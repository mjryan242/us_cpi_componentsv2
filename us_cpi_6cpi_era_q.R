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
est <- function(d) mct(R2ConnectednessQ2(d,window.size=NULL,nlag=2,tau=0.9,shrink=TRUE,drop_own_lags=FALSE,progbar=FALSE))
c0 <- est(window(dat, start=as.yearqtr("1983 Q1"), end=as.yearqtr("2019 Q4")))
gi <- est(window(dat, end=as.yearqtr("2019 Q4")))
cv <- est(window(dat, start=as.yearqtr("1983 Q1")))
mk <- function(x){ d<-x-c0; data.frame(d_Overall=round(d["Overall"],2), d_Contemp=round(d["Contemp"],2),
  d_Lagged=round(d["Lagged"],2), Lagged_share_pct=round(100*d["Lagged"]/d["Overall"],1)) }
out <- rbind(cbind(Episode="GreatInflation (1967-82)", mk(gi)),
             cbind(Episode="COVID (2020-26)",          mk(cv)))
rownames(out) <- NULL
dir.create("results_6cpi_era_q", showWarnings=FALSE)
write.csv(out, "results_6cpi_era_q/episode_marginal_decomposition.csv", row.names=FALSE)
cat(sprintf("core (1983-2019) tau0.9: Overall %.1f, Contemp %.1f, Lagged %.1f\n", c0["Overall"],c0["Contemp"],c0["Lagged"]))
print(out, row.names=FALSE)
cat("\nSaved results_6cpi_era_q/episode_marginal_decomposition.csv\nDone.\n")
