# ============================================================
# us_cpi_6cpi_driver_q.R -- QUARTERLY six-CPI system (NO oil).
#   Quarterly analogue of us_cpi_6cpi_driver.R; each series aggregated
#   monthly index -> quarterly period-average -> QoQ % (quarterly_utils.R).
# Output: results_6cpi_q/, nlag1_6cpi_q/, plots_6cpi_q/
# ============================================================
pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R"); source("quarterly_utils.R"); W <- 0.60

codes <- c(Food="CPIUFDSL", Gasoline="CUSR0000SETB01", HHEnergy="CUSR0000SEHF",
           CoreGoods="CUSR0000SACL1E", Shelter="CUSR0000SAH1", CoreServices="CUSR0000SASLE")
cat("Downloading", length(codes), "CPI series from FRED...\n")
sl <- list(); for (nm in names(codes)) { s <- getSymbols(codes[[nm]], src="FRED", auto.assign=FALSE); colnames(s)<-nm; sl[[nm]]<-s }
idx <- do.call(merge, sl)[, names(codes)]
pcz <- to_qpc(idx); colnames(pcz) <- names(codes)
exS <- (pcz[,"CoreServices"]-W*pcz[,"Shelter"])/(1-W)
dat <- cbind(Food=pcz[,"Food"], Gasoline=pcz[,"Gasoline"], HHEnergy=pcz[,"HHEnergy"],
             CoreGoods=pcz[,"CoreGoods"], Shelter=pcz[,"Shelter"], CoreServ_xS=exS)
dat <- na.omit(dat); storage.mode(dat) <- "double"
cat(sprintf("\nSix-CPI QUARTERLY: %d quarters, %s to %s, %d series\n",
            nrow(dat), as.character(start(dat)), as.character(end(dat)), ncol(dat)))

make_conn_table <- function(res){ C<-res$CT[,,1,1]*100; L<-res$CT[,,1,2]*100; O<-C+L
  od<-function(M){diag(M)<-0;M}; Co<-od(C);Lo<-od(L);Oo<-od(O)
  list(TCI_C=mean(colSums(Co)),TCI_L=mean(colSums(Lo)),TCI=mean(colSums(Oo))) }
tci_by_quantile <- function(rl) do.call(rbind, lapply(names(rl), function(lbl){ tb<-make_conn_table(rl[[lbl]])
  data.frame(Tau=as.numeric(sub("tau_","",lbl)), Overall=round(tb$TCI,3),
             Contemporaneous=round(tb$TCI_C,3), Lagged=round(tb$TCI_L,3)) }))
taus <- c(0.1,0.3,0.5,0.7,0.9)
est <- function(nlag){ o<-list(); for (tau in taus){ lbl<-sprintf("tau_%.1f",tau)
  o[[lbl]] <- R2ConnectednessQ2(dat, window.size=NULL, nlag=nlag, tau=tau, shrink=TRUE, drop_own_lags=FALSE, progbar=FALSE) }; o }
res2 <- est(2); res1 <- est(1)
dir.create("results_6cpi_q",showWarnings=FALSE); dir.create("nlag1_6cpi_q",showWarnings=FALSE)
write.csv(tci_by_quantile(res2), "results_6cpi_q/TCI_by_quantile.csv", row.names=FALSE)
write.csv(tci_by_quantile(res1), "nlag1_6cpi_q/TCI_by_quantile.csv", row.names=FALSE)
t2 <- tci_by_quantile(res2)
cat("\n=== TCI by quantile (quarterly, nlag=2) ===\n"); print(t2, row.names=FALSE)
cat(sprintf(">> Broadening (quarterly nlag=2): TCI(0.9)-TCI(0.5) = %.2f  (%.2fx)\n",
            t2$Overall[t2$Tau==0.9]-t2$Overall[t2$Tau==0.5], t2$Overall[t2$Tau==0.9]/t2$Overall[t2$Tau==0.5]))
cat("\nSaved results_6cpi_q/, nlag1_6cpi_q/\nDone.\n")
