# ============================================================
# us_cpi_bootstrap_r2q.R -- bootstrap CIs for the HEADLINE inflation-broadening
#   results (six BLS CPI subcomponents, NO oil -- matches us_cpi_6cpi_driver.R).
#
# ADAPTED from the original 7-series (6 groups + WTI oil) version: the oil merge
# is removed and the system is the same six no-oil columns as us_cpi_6cpi_driver.R
# (Food, Gasoline, HHEnergy, CoreGoods, Shelter, CoreServ_xS), so the resulting CIs
# bracket the numbers actually reported in tables/h1.tex. The hardcoded setwd() to
# the old project path was also removed -- run this from the v2 project root.
#
# Method: FIXED-DESIGN quantile-wild bootstrap (in the Caporin-Bonaccolto-Shahzad
#   wild family, consistent with CtR2QBoot's residual resampling, but NON-recursive).
#   Why fixed-design: the recursive QVAR simulation explodes at tau=0.9 here -- the
#   tau=0.9 quantile-VAR is non-stationary, so a recursive path diverges (~1e24) and
#   collapses every draw to a constant. Fixed-design forms bootstrap values as
#   (fitted + wild residual) using the ACTUAL lags -- no feedback, no explosion.
#
#   For each equation j: fit a quantile regression of y_j on all lags at tau (giving
#   fitted_j and residuals); wild-resample the residuals (Caporin et al. construction:
#   leverage H, density-at-zero via quantreg::akj, delta-split weights {-2tau,2(1-tau)});
#   form y*_j = fitted_j + wild residual; pair y* with the ACTUAL lags and recompute the
#   connectedness via conn_from_Z() (the core of nearPD R2QConnectedness on an embedded
#   matrix). STATISTIC = the NEW headline estimator (nearPD R2QConnectedness), NOT clipping.
#
# Produces (full sample, nlag=2): bias-corrected percentile CIs for TCI at each tau
#   (Overall/Contemp/Lagged), a CI + one-sided p for the broadening gap TCI(0.9)-TCI(0.5),
#   and CIs for the tau=0.9 directional NET of each series.
# Output: results_boot_r2q/, plots_boot_r2q/
#
# NOTE: B is set to 10 below for a quick smoke test (~seconds, not a real CI).
#   Set B <- 1000 (or BOOT_B=1000 env var) for a publication-quality run.
# ============================================================
suppressMessages({
  library(quantmod); library(zoo); library(xts); library(quantreg)
  library(MASS); library(corpcor); library(Matrix); library(ConnectednessApproach)
  library(parallel)
})
source("use_R2Q.R")   # QuantileCorrelation, .ShrinkCorrelation (PSD-targeted), .IsPSD, R2QConnectedness
W <- 0.60

## ---- config (set B for the full run) ----
TAUS  <- c(0.1, 0.3, 0.5, 0.7, 0.9)
NLAG  <- 2
B     <- 10                                     # <-- quick smoke test (was 1000 for the full run)
if (nzchar(Sys.getenv("BOOT_B"))) B <- as.integer(Sys.getenv("BOOT_B"))
DELTA <- 0.3
ALPHA <- 0.05
NWORK <- max(1L, parallel::detectCores() - 2L)
if (nzchar(Sys.getenv("BOOT_NWORK"))) NWORK <- as.integer(Sys.getenv("BOOT_NWORK"))
SEED  <- 20260625
set.seed(SEED)
dir.create("results_boot_r2q", showWarnings=FALSE); dir.create("plots_boot_r2q", showWarnings=FALSE)

## ---- data: six-CPI headline system (no oil; matches us_cpi_6cpi_driver.R) ----
codes <- c(Food="CPIUFDSL", Gasoline="CUSR0000SETB01", HHEnergy="CUSR0000SEHF",
           CoreGoods="CUSR0000SACL1E", Shelter="CUSR0000SAH1", CoreServices="CUSR0000SASLE")
cat("Downloading CPI series from FRED...\n")
sl <- list(); for (nm in names(codes)) { s <- getSymbols(codes[[nm]], src="FRED", auto.assign=FALSE); colnames(s)<-nm; sl[[nm]]<-s }
idx <- do.call(merge, sl)[, names(codes)]
pcz <- zoo(coredata(100*(idx/lag.xts(idx,k=1)-1)), as.yearmon(index(idx))); colnames(pcz)<-names(codes)
exS <- (pcz[,"CoreServices"] - W*pcz[,"Shelter"])/(1-W)
cpi6 <- cbind(Food=pcz[,"Food"], Gasoline=pcz[,"Gasoline"], HHEnergy=pcz[,"HHEnergy"],
              CoreGoods=pcz[,"CoreGoods"], Shelter=pcz[,"Shelter"], CoreServ_xS=exS)
dat <- zoo(coredata(cpi6), index(cpi6))
dat <- dat[, c("Food","Gasoline","HHEnergy","CoreGoods","Shelter","CoreServ_xS")]
dat <- na.omit(dat); storage.mode(dat) <- "double"; X0 <- as.matrix(dat); NMS <- colnames(X0); K <- ncol(X0)
cat(sprintf("Sample: %d months, %s to %s; %d series. B=%d, NWORK=%d, nlag=%d\n\n",
            nrow(dat), format(start(dat)), format(end(dat)), K, B, NWORK, NLAG))

## ---- connectedness from an embedded matrix Z (core of nearPD R2QConnectedness) ----
# Z columns: [current k | lag1 k | ... | lag_nlag k]. Returns flat named stat vector.
conn_from_Z <- function(Z, k, nlag, tau, nms, method="fn"){
  R <- QuantileCorrelation(Z, tau=tau, method=method)
  R <- .ShrinkCorrelation(R, Z, verbose=FALSE); R <- 0.5*(R+t(R))
  if (!.IsPSD(R, tol=1e-8)) { R <- as.matrix(Matrix::nearPD(R, corr=TRUE)$mat); R <- R + 1e-10*diag(nrow(R)) }
  C <- matrix(0,k,k); L <- matrix(0,k,k)
  for (i in 1:k){
    ryx <- R[-i, i, drop=FALSE]; rxx <- R[-i,-i,drop=FALSE]
    e <- eigen(rxx, symmetric=TRUE); e$values <- pmax(round(e$values,3),0)
    root <- e$vectors %*% diag(sqrt(e$values), nrow=length(e$values)) %*% t(e$vectors)
    cd <- root^2 %*% (MASS::ginv(root) %*% ryx)^2
    C[i,-i] <- cd[1:(k-1)]
    if (nlag>0) L[i,] <- apply(array(cd[-(1:(k-1))], c(1,k,nlag)), 1:2, sum)
  }
  C <- C*100; L <- L*100; O <- C+L
  offd <- function(M){diag(M)<-0;M}; Co<-offd(C); Lo<-offd(L); Oo<-offd(O)
  net <- colSums(Oo) - rowSums(Oo); names(net) <- nms
  c(at=mean(colSums(Oo)), ct=mean(colSums(Co)), lt=mean(colSums(Lo)),
    setNames(net, paste0("NET_", nms)))
}

## ---- fixed-design fit: per-equation QR at tau -> fitted, residuals, leverage ----
fd_fit <- function(X, p, tau, method="fn"){
  X <- as.matrix(X); N <- ncol(X); z <- embed(X, p+1)
  Xreg <- z[,-(1:N),drop=FALSE]; Xhat <- cbind(1, Xreg)
  Y <- z[,1:N,drop=FALSE]; fitted <- matrix(0,nrow(Y),N); Res <- matrix(0,nrow(Y),N)
  for (j in 1:N){ f <- quantreg::rq(Y[,j]~Xreg, tau=tau, method=method)
    fitted[,j] <- Xhat %*% coef(f); Res[,j] <- Y[,j]-fitted[,j] }
  H <- pmin(pmax(rowSums((Xhat %*% qr.solve(crossprod(Xhat), diag(ncol(Xhat)))) * Xhat),0),1)
  list(Z=z, lags=Xreg, fitted=fitted, Res=Res, H=H, N=N)
}
## ---- one fixed-design wild resample -> embedded matrix Z_star ----
fd_sample <- function(fit, tau, delta, f0_floor=1e-8){
  Res <- fit$Res; H <- fit$H; Teff <- nrow(Res); N <- fit$N
  n_common <- floor(Teff*delta); pos_all <- seq_len(Teff)
  pos_com <- if(n_common>0) sort(sample(pos_all,n_common)) else integer(0); pos_ind <- setdiff(pos_all,pos_com)
  draw_w <- function(m) if(m<=0) numeric(0) else sample(c(-2*tau,2*(1-tau)), prob=c(tau,1-tau), size=m, replace=TRUE)
  w_com <- draw_w(length(pos_com)); Res_star <- matrix(0,Teff,N)
  for (j in 1:N){ w_ind <- draw_w(length(pos_ind)); w <- numeric(Teff)
    if(length(pos_com)) w[pos_com] <- w_com; if(length(pos_ind)) w[pos_ind] <- w_ind
    f0 <- quantreg::akj(Res[,j], z=0)$dens; if(!is.finite(f0)||f0<=0) f0 <- f0_floor
    RR <- Res[,j] + H*(tau-as.integer(Res[,j]<0))/f0; Res_star[,j] <- w*abs(RR) }
  cbind(fit$fitted + Res_star, fit$lags)        # bootstrap current + ACTUAL lags (no recursion)
}

## ---- bootstrap one tau ----
boot_one <- function(cl, fit, k, p, tau, delta, nms, B){
  base <- conn_from_Z(fit$Z, k, p, tau, nms)
  draws <- parLapply(cl, seq_len(B), function(b, fit,k,p,tau,delta,nms){
    Zs <- tryCatch(fd_sample(fit, tau, delta), error=function(e) NULL)
    if (is.null(Zs)) return(NULL)
    tryCatch(conn_from_Z(Zs, k, p, tau, nms), error=function(e) NULL)
  }, fit=fit, k=k, p=p, tau=tau, delta=delta, nms=nms)
  draws <- draws[!vapply(draws, is.null, logical(1))]
  D <- do.call(rbind, draws); D <- D[stats::complete.cases(D), , drop=FALSE]
  bias <- base - colMeans(D); Dbc <- sweep(D, 2, bias, FUN="+")
  list(base=base, Dbc=Dbc, B_eff=nrow(D))
}
pctl_ci <- function(v, alpha) stats::quantile(v, c(alpha/2, 1-alpha/2), na.rm=TRUE)

## ---- run ----
cl <- makeCluster(NWORK); on.exit(stopCluster(cl), add=TRUE)
clusterCall(cl, setwd, getwd())
clusterEvalQ(cl, { suppressMessages({library(zoo);library(quantreg);library(MASS);
  library(corpcor);library(Matrix);library(ConnectednessApproach)}); source("use_R2Q.R"); TRUE })
clusterExport(cl, c("conn_from_Z","fd_sample"))
clusterSetRNGStream(cl, SEED)

store <- list(); ci_rows <- list()
for (tau in TAUS){
  fit <- fd_fit(X0, NLAG, tau)
  t0 <- Sys.time(); bo <- boot_one(cl, fit, K, NLAG, tau, DELTA, NMS, B); el <- as.numeric(Sys.time()-t0, units="secs")
  store[[as.character(tau)]] <- bo
  for (ty in c("at","ct","lt")){ ci <- pctl_ci(bo$Dbc[,ty], ALPHA)
    ci_rows[[length(ci_rows)+1]] <- data.frame(nlag=NLAG, tau=tau, type=ty,
      base=round(bo$base[ty],3), lower=round(ci[1],3), upper=round(ci[2],3), B_eff=bo$B_eff) }
  cat(sprintf("  tau=%.1f  TCI=%.2f [%.2f, %.2f]  (%.1fs, B_eff=%d)\n", tau, bo$base["at"],
      pctl_ci(bo$Dbc[,"at"],ALPHA)[1], pctl_ci(bo$Dbc[,"at"],ALPHA)[2], el, bo$B_eff))
}
stopCluster(cl)
ci_df <- do.call(rbind, ci_rows); write.csv(ci_df, "results_boot_r2q/TCI_CI_by_quantile.csv", row.names=FALSE)

## ---- broadening gap ----
hi <- store[["0.9"]]; lo <- store[["0.5"]]; nb <- min(nrow(hi$Dbc), nrow(lo$Dbc)); gap_rows <- list()
for (ty in c("at","ct","lt")){
  diffb <- hi$Dbc[1:nb,ty] - lo$Dbc[1:nb,ty]; gap_base <- hi$base[ty] - lo$base[ty]
  ci <- pctl_ci(diffb, ALPHA); pval <- mean(diffb <= 0)
  gap_rows[[length(gap_rows)+1]] <- data.frame(nlag=NLAG, type=ty, gap_base=round(gap_base,3),
    lower=round(ci[1],3), upper=round(ci[2],3), p_one_sided=signif(pval,3), B_eff=nb) }
gap_df <- do.call(rbind, gap_rows); write.csv(gap_df, "results_boot_r2q/broadening_gap_CI.csv", row.names=FALSE)

## ---- tau=0.9 NET CIs ----
b9 <- store[["0.9"]]; netcols <- grep("^NET_", colnames(b9$Dbc), value=TRUE)
net_df <- do.call(rbind, lapply(netcols, function(nc){ ci <- pctl_ci(b9$Dbc[,nc], ALPHA)
  data.frame(series=sub("^NET_","",nc), base=round(b9$base[nc],2), lower=round(ci[1],2), upper=round(ci[2],2)) }))
net_df <- net_df[order(-net_df$base),]; write.csv(net_df, "results_boot_r2q/NET_CI_tau0.9.csv", row.names=FALSE)

## ---- plot ----
at <- ci_df[ci_df$type=="at",]; at <- at[order(at$tau),]
pdf("plots_boot_r2q/TCI_bands.pdf", width=7, height=5)
plot(at$tau, at$base, type="n", ylim=c(0, max(at$upper)*1.05), xlab=expression(tau),
     ylab="Total Connectedness Index (%)", main="Inflation broadening with 95% bootstrap CIs")
polygon(c(at$tau, rev(at$tau)), c(at$lower, rev(at$upper)), col=adjustcolor("#2166ac",0.18), border=NA)
lines(at$tau, at$base, type="b", pch=19, lwd=2); dev.off()

cat("\n=== TCI CIs (Overall) ===\n"); print(ci_df[ci_df$type=="at",], row.names=FALSE)
cat("\n=== Broadening gap ===\n"); print(gap_df, row.names=FALSE)
cat("\n=== tau=0.9 NET CIs ===\n"); print(net_df, row.names=FALSE)
cat("\nSaved results_boot_r2q/ and plots_boot_r2q/\nDone.\n")
print(Sys.time())
# or
cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")