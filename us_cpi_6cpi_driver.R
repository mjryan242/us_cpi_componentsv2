# ============================================================
# us_cpi_6cpi_driver.R -- HEADLINE six-CPI system (NO oil).
#   Six BLS expenditure-weighted CPI subcomponents (a complete,
#   non-overlapping decomposition of headline CPI):
#     Food          CPIUFDSL
#     Gasoline      CUSR0000SETB01
#     HHEnergy      CUSR0000SEHF
#     CoreGoods     CUSR0000SACL1E
#     Shelter       CUSR0000SAH1
#     CoreServ_xS   = core services (CUSR0000SASLE) less shelter, w=0.60
#
# Estimator: R2QConnectedness (nearPD) via use_R2Q.R.
# Main: nlag=2.  Robustness: nlag=1.  Output: results_6cpi/, nlag1_6cpi/, plots_6cpi/
# ============================================================
pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R"); W <- 0.60

codes <- c(Food="CPIUFDSL", Gasoline="CUSR0000SETB01", HHEnergy="CUSR0000SEHF",
           CoreGoods="CUSR0000SACL1E", Shelter="CUSR0000SAH1", CoreServices="CUSR0000SASLE")
cat("Downloading", length(codes), "CPI series from FRED...\n")
sl <- list(); for (nm in names(codes)) { s <- getSymbols(codes[[nm]], src="FRED", auto.assign=FALSE); colnames(s)<-nm; sl[[nm]]<-s }
idx <- do.call(merge, sl)[, names(codes)]
pcz <- zoo(coredata(100*(idx/lag.xts(idx,k=1)-1)), as.yearmon(index(idx))); colnames(pcz)<-names(codes)
exS <- (pcz[,"CoreServices"]-W*pcz[,"Shelter"])/(1-W)
dat <- cbind(Food=pcz[,"Food"], Gasoline=pcz[,"Gasoline"], HHEnergy=pcz[,"HHEnergy"],
             CoreGoods=pcz[,"CoreGoods"], Shelter=pcz[,"Shelter"], CoreServ_xS=exS)
dat <- na.omit(dat); storage.mode(dat) <- "double"
cat(sprintf("\nSix-CPI system (NO oil): %d months, %s to %s, %d series (w=%.2f)\n",
            nrow(dat), format(start(dat)), format(end(dat)), ncol(dat), W))

make_conn_table <- function(res){ C<-res$CT[,,1,1]*100; L<-res$CT[,,1,2]*100; O<-C+L
  od<-function(M){diag(M)<-0;M}; Co<-od(C);Lo<-od(L);Oo<-od(O)
  list(Cp=C,Lp=L,Op=O,FROM_C=rowSums(Co),FROM_L=rowSums(Lo),FROM_T=rowSums(Oo),
       TO_C=colSums(Co),TO_L=colSums(Lo),TO_T=colSums(Oo),NET_T=colSums(Oo)-rowSums(Oo),
       TCI_C=mean(colSums(Co)),TCI_L=mean(colSums(Lo)),TCI=mean(colSums(Oo))) }
tci_by_quantile <- function(rl) do.call(rbind, lapply(names(rl), function(lbl){ tb<-make_conn_table(rl[[lbl]])
  data.frame(Tau=as.numeric(sub("tau_","",lbl)), Overall=round(tb$TCI,3),
             Contemporaneous=round(tb$TCI_C,3), Lagged=round(tb$TCI_L,3)) }))
write_dir <- function(rl, outdir){ dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  write.csv(tci_by_quantile(rl), file.path(outdir,"TCI_by_quantile.csv"), row.names=FALSE)
  dr <- list(); for (lbl in names(rl)){ tau<-as.numeric(sub("tau_","",lbl)); tb<-make_conn_table(rl[[lbl]]); nm<-rownames(tb$Op)
    for (i in seq_along(nm)) dr[[length(dr)+1]] <- data.frame(Tau=tau, Series=nm[i],
      FROM_Total=round(tb$FROM_T[i],3), TO_C=round(tb$TO_C[i],3), TO_L=round(tb$TO_L[i],3),
      TO_Total=round(tb$TO_T[i],3), NET_Total=round(tb$NET_T[i],3)) }
  write.csv(do.call(rbind,dr), file.path(outdir,"conn_directional.csv"), row.names=FALSE) }

taus <- c(0.1,0.3,0.5,0.7,0.9)
est <- function(nlag){ o<-list(); for (tau in taus){ lbl<-sprintf("tau_%.1f",tau)
  o[[lbl]] <- R2ConnectednessQ2(dat, window.size=NULL, nlag=nlag, tau=tau, shrink=TRUE, drop_own_lags=FALSE, progbar=FALSE) }; o }
res2 <- est(2); res1 <- est(1)
write_dir(res2,"results_6cpi"); write_dir(res1,"nlag1_6cpi"); dir.create("plots_6cpi",showWarnings=FALSE)
t2 <- tci_by_quantile(res2); t1 <- tci_by_quantile(res1)
cat("\n=== TCI by quantile (nlag=2, main) ===\n"); print(t2, row.names=FALSE)
cat("\n=== TCI by quantile (nlag=1) ===\n"); print(t1, row.names=FALSE)
cat(sprintf("\n>> Broadening (nlag=2): TCI(0.9)-TCI(0.5) = %.2f  (%.2fx)\n",
            t2$Overall[t2$Tau==0.9]-t2$Overall[t2$Tau==0.5], t2$Overall[t2$Tau==0.9]/t2$Overall[t2$Tau==0.5]))

# decomposition plot (Overall black, Contemp blue, Lagged red)
pdf("plots_6cpi/TCI_decomp.pdf", width=8, height=5)
yr <- range(c(t2$Overall,t2$Contemporaneous,t2$Lagged))
plot(t2$Tau,t2$Overall,type="b",pch=19,lwd=2,ylim=yr,xlab=expression(tau),ylab="TCI (%)",
     main="Six-CPI system: TCI decomposition (nlag=2)")
lines(t2$Tau,t2$Contemporaneous,type="b",pch=17,lwd=2,col="#2166ac")
lines(t2$Tau,t2$Lagged,type="b",pch=15,lwd=2,col="#d6604d")
legend("topleft",legend=c("Overall","Contemporaneous","Lagged"),
       col=c("black","#2166ac","#d6604d"),pch=c(19,17,15),lwd=2,bty="n")
dev.off()
cat("\nSaved results_6cpi/, nlag1_6cpi/, plots_6cpi/TCI_decomp.pdf\nDone.\n")
