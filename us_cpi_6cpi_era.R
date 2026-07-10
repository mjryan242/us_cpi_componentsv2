# ============================================================
# us_cpi_6cpi_era.R -- episode mechanism for the six-CPI system (NO oil).
#   Marginal dTCI(0.9) vs the 1983-2019 core (nested subtraction):
#     GreatInflation marginal = dropCOVID(1967-2019) - dropBoth(1983-2019)
#     COVID marginal          = dropGI(1983-2026)    - dropBoth(1983-2019)
# Output: results_6cpi_era/episode_marginal_decomposition.csv
# ============================================================
pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R"); W <- 0.60
codes <- c(Food="CPIUFDSL", Gasoline="CUSR0000SETB01", HHEnergy="CUSR0000SEHF",
           CoreGoods="CUSR0000SACL1E", Shelter="CUSR0000SAH1", CoreServices="CUSR0000SASLE")
sl <- list(); for (nm in names(codes)) { s <- getSymbols(codes[[nm]], src="FRED", auto.assign=FALSE); colnames(s)<-nm; sl[[nm]]<-s }
idx <- do.call(merge, sl)[, names(codes)]
pcz <- zoo(coredata(100*(idx/lag.xts(idx,k=1)-1)), as.yearmon(index(idx))); colnames(pcz)<-names(codes)
exS <- (pcz[,"CoreServices"]-W*pcz[,"Shelter"])/(1-W)
dat <- cbind(Food=pcz[,"Food"], Gasoline=pcz[,"Gasoline"], HHEnergy=pcz[,"HHEnergy"],
             CoreGoods=pcz[,"CoreGoods"], Shelter=pcz[,"Shelter"], CoreServ_xS=exS)
dat <- na.omit(dat); storage.mode(dat) <- "double"
mct <- function(res){ C<-res$CT[,,1,1]*100; L<-res$CT[,,1,2]*100; O<-C+L; od<-function(M){diag(M)<-0;M}
  c(Overall=mean(colSums(od(O))), Contemp=mean(colSums(od(C))), Lagged=mean(colSums(od(L)))) }
est <- function(d) mct(R2ConnectednessQ2(d,window.size=NULL,nlag=2,tau=0.9,shrink=TRUE,drop_own_lags=FALSE,progbar=FALSE))
c0 <- est(window(dat, start=as.yearmon("Jan 1983"), end=as.yearmon("Dec 2019")))  # dropBoth core
gi <- est(window(dat, end=as.yearmon("Dec 2019")))                                # 1967-2019 (dropCOVID)
cv <- est(window(dat, start=as.yearmon("Jan 1983")))                              # 1983-2026 (dropGI)
mk <- function(x){ d<-x-c0; data.frame(d_Overall=round(d["Overall"],2), d_Contemp=round(d["Contemp"],2),
  d_Lagged=round(d["Lagged"],2), Lagged_share_pct=round(100*d["Lagged"]/d["Overall"],1)) }
out <- rbind(cbind(Episode="GreatInflation (1967-82)", mk(gi)),
             cbind(Episode="COVID (2020-26)",          mk(cv)))
rownames(out) <- NULL
dir.create("results_6cpi_era", showWarnings=FALSE)
write.csv(out, "results_6cpi_era/episode_marginal_decomposition.csv", row.names=FALSE)
cat(sprintf("core (1983-2019) tau0.9: Overall %.1f, Contemp %.1f, Lagged %.1f\n", c0["Overall"],c0["Contemp"],c0["Lagged"]))
print(out, row.names=FALSE)
cat("\nSaved results_6cpi_era/episode_marginal_decomposition.csv\nDone.\n")
