# ============================================================
# us_cpi_6cpi_era_q.R -- QUARTERLY episode mechanism, six-CPI (NO oil).
#   Quarterly analogue of us_cpi_6cpi_era.R.
# Output: results_6cpi_era_q/episode_marginal_decomposition.csv
# ============================================================
pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R"); source("quarterly_utils.R"); W <- 0.60
codes <- c(Food="CPIUFDSL", Gasoline="CUSR0000SETB01", HHEnergy="CUSR0000SEHF",
           CoreGoods="CUSR0000SACL1E", Shelter="CUSR0000SAH1", CoreServices="CUSR0000SASLE")
sl <- list(); for (nm in names(codes)) { s <- getSymbols(codes[[nm]], src="FRED", auto.assign=FALSE); colnames(s)<-nm; sl[[nm]]<-s }
idx <- do.call(merge, sl)[, names(codes)]
pcz <- to_qpc(idx); colnames(pcz) <- names(codes)
exS <- (pcz[,"CoreServices"]-W*pcz[,"Shelter"])/(1-W)
dat <- cbind(Food=pcz[,"Food"], Gasoline=pcz[,"Gasoline"], HHEnergy=pcz[,"HHEnergy"],
             CoreGoods=pcz[,"CoreGoods"], Shelter=pcz[,"Shelter"], CoreServ_xS=exS)
dat <- na.omit(dat); storage.mode(dat) <- "double"
mct <- function(res){ C<-res$CT[,,1,1]*100; L<-res$CT[,,1,2]*100; O<-C+L; od<-function(M){diag(M)<-0;M}
  c(Overall=mean(colSums(od(O))), Contemp=mean(colSums(od(C))), Lagged=mean(colSums(od(L)))) }
est <- function(d, tau) mct(R2ConnectednessQ2(d,window.size=NULL,nlag=2,tau=tau,shrink=TRUE,drop_own_lags=FALSE,progbar=FALSE))
# episode quarter counts (for the small-sample point in the paper's H3 discussion)
nq <- c(GI=nrow(window(dat, end=as.yearqtr("1982 Q4"))),
        Core=nrow(window(dat, start=as.yearqtr("1983 Q1"), end=as.yearqtr("2019 Q4"))),
        COVID=nrow(window(dat, start=as.yearqtr("2020 Q1"))))
cat(sprintf("episode quarter counts: GI %d, Core %d, COVID %d\n", nq["GI"], nq["Core"], nq["COVID"]))
out <- NULL; base <- NULL
for (tau in c(0.7, 0.9)) {
  c0 <- est(window(dat, start=as.yearqtr("1983 Q1"), end=as.yearqtr("2019 Q4")), tau)
  gi <- est(window(dat, end=as.yearqtr("2019 Q4")), tau)
  cv <- est(window(dat, start=as.yearqtr("1983 Q1")), tau)
  mk <- function(x){ d<-x-c0; data.frame(d_Overall=round(d["Overall"],2), d_Contemp=round(d["Contemp"],2),
    d_Lagged=round(d["Lagged"],2), Lagged_share_pct=round(100*d["Lagged"]/d["Overall"],1)) }
  out <- rbind(out,
    cbind(Tau=tau, Episode="GreatInflation (1967-82)", mk(gi)),
    cbind(Tau=tau, Episode="COVID (2020-26)",          mk(cv)))
  base <- rbind(base, data.frame(Tau=tau, Overall=round(c0["Overall"],2),
    Contemp=round(c0["Contemp"],2), Lagged=round(c0["Lagged"],2)))
  cat(sprintf("core (1983-2019) tau%.1f: Overall %.1f, Contemp %.1f, Lagged %.1f\n",
      tau, c0["Overall"],c0["Contemp"],c0["Lagged"]))
}
rownames(out) <- NULL; rownames(base) <- NULL
dir.create("results_6cpi_era_q", showWarnings=FALSE)
write.csv(out,  "results_6cpi_era_q/episode_marginal_decomposition.csv", row.names=FALSE)
write.csv(base, "results_6cpi_era_q/core_base_levels.csv", row.names=FALSE)
print(out, row.names=FALSE)
cat("\nSaved results_6cpi_era_q/episode_marginal_decomposition.csv + core_base_levels.csv\nDone.\n")
