# ============================================================
# us_cpi_wage_taulevel.R -- level-matched wage-CPI connectedness (H2 analogue
#   of the issue-#4 exercise). Table 2 compares the GI-inclusive (1964-2019)
#   and post-GI (1983-2026) subperiods at the SAME tau, which is a different
#   inflation LEVEL. Here we instead ask: at the tau in the 1964-2019 sample
#   whose CPI-inflation level matches tau=0.9 of the 1983-2026 sample, what is
#   the wage-CPI connectedness?
#
# Data build + windows replicate us_cpi_wage_spiral.R exactly.
# ============================================================
pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R")

cpi <- getSymbols("CPIAUCSL", src="FRED", auto.assign=FALSE)
cpi_pc <- zoo(coredata(100*(cpi/lag.xts(cpi,k=1)-1)), as.yearmon(index(cpi))); colnames(cpi_pc) <- "CPI"
w_raw <- read.csv("AHETPI.csv", stringsAsFactors=FALSE); w_raw$observation_date <- as.Date(w_raw$observation_date)
wz <- zoo(w_raw$AHETPI, as.yearmon(w_raw$observation_date)); w_pc <- 100*(wz/stats::lag(wz,-1)-1)
dat <- merge(Wages=w_pc, CPI=cpi_pc, all=FALSE); dat <- dat[, c("Wages","CPI")]; dat <- na.omit(dat)

pre  <- window(dat, end=as.yearmon("Dec 2019"))     # 1964-2019 (with Great Inflation)
post <- window(dat, start=as.yearmon("Jan 1983"))   # 1983-2026 (post-Great-Inflation)

## ---- connectedness (same metrics as us_cpi_wage_spiral.R) ----
cpi_pre  <- as.numeric(pre[,"CPI"]); cpi_post <- as.numeric(post[,"CPI"])
ann <- function(mm) 100*((1+mm/100)^12 - 1)                            # m/m % -> annualised %
metrics <- function(C, L){ O <- C+L; od <- function(M){diag(M)<-0;M}
  list(TCI=mean(colSums(od(O))), TCI_C=mean(colSums(od(C))), TCI_L=mean(colSums(od(L))),
       w2p_L=L["CPI","Wages"], p2w_L=L["Wages","CPI"]) }   # lagged Wages->CPI ; lagged CPI->Wages
genizi <- function(d, tau){ r <- R2ConnectednessQ2(d, window.size=NULL, nlag=2, tau=tau,
    shrink=TRUE, drop_own_lags=FALSE, progbar=FALSE); metrics(r$CT[,,1,1]*100, r$CT[,,1,2]*100) }
report <- function(lab, d, tau){ m <- genizi(d, tau)
  cat(sprintf("%-46s tau=%.3f  TCI=%5.2f  Wages->CPI(L)=%5.2f  CPI->Wages(L)=%5.2f\n",
      lab, tau, m$TCI, m$w2p_L, m$p2w_L)); invisible(m) }

cat("\n--- Table 2 cells (same-tau), for reference ---\n")
report("1983-2026 @ tau=0.7 (Table 2 as printed)", post, 0.7)
report("1983-2026 @ tau=0.9 (Table 2 as printed)", post, 0.9)
report("1964-2019 @ tau=0.7 (Table 2 as printed)", pre,  0.7)
report("1964-2019 @ tau=0.9 (Table 2 as printed)", pre,  0.9)

## ---- level-match the 1964-2019 sample to each post-1983 tau in turn ----
rows <- list()
for (post_tau in c(0.7, 0.9)) {
  lev <- as.numeric(quantile(cpi_post, post_tau, na.rm=TRUE, type=8))  # target level (m/m %)
  tau_star <- mean(cpi_pre <= lev, na.rm=TRUE)                         # empirical CDF in pre sample
  cat(sprintf("\n=== match to post-1983 tau=%.1f  (CPI %.3f%% m/m = %.2f%% ann.) -> pre tau*=%.3f ===\n",
              post_tau, lev, ann(lev), tau_star))
  m <- report("1964-2019 @ level-matched tau", pre, tau_star)
  cat(sprintf("  RESULT: TCI=%.2f, lagged Wages->CPI=%.2f, lagged CPI->Wages=%.2f\n",
      m$TCI, m$w2p_L, m$p2w_L))
  rows[[length(rows)+1]] <- data.frame(post_tau=post_tau, level_mm=round(lev,3),
    level_ann=round(ann(lev),2), tau_matched=round(tau_star,3),
    TCI=round(m$TCI,2), Wages_to_CPI_L=round(m$w2p_L,2), CPI_to_Wages_L=round(m$p2w_L,2))
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, "results_wage_spiral/tau_level_matched.csv", row.names=FALSE)
cat("\nSaved results_wage_spiral/tau_level_matched.csv\n")
