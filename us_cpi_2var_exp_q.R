# ============================================================
# us_cpi_2var_exp_q.R -- QUARTERLY 2-variable system:
#   MICH_QTR (Michigan 1yr inflation expectations, quarterly mean, %)
#     -- used AS-IS, already a % expectation, NO transform --
#   vs Overall CPI (CPIAUCSL) aggregated to quarterly QoQ % change.
#
# This is the quarterly analogue of us_cpi_2var_exp.R.  Its PAYOFF is
# span: MICH_QTR begins 1960Q1, so unlike monthly MICH (1978+) it
# covers the full Great Inflation INCLUDING the 1973-77 un-anchoring.
#
# Estimator: R2QConnectedness (nearPD) via use_R2Q.R.
# nlag=2 main (6-month horizon), nlag=1 robustness. tau in {0.1..0.9}.
# In a 2-var system NET_MICH = -NET_CPI: NET_MICH>0 => expectations
# LEAD CPI (transmitter); NET_MICH<0 => expectations FOLLOW.
# Output: results_2var_exp_q/
# ============================================================

pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R")
source("quarterly_utils.R")

# CPI -> quarterly QoQ %
cpi   <- getSymbols("CPIAUCSL", src="FRED", auto.assign=FALSE)
cpi_q <- to_qpc(cpi); colnames(cpi_q) <- "CPI"
# MICH_QTR: quarterly Michigan inflation expectations, already %, as-is
mich_q <- read_mich_qtr("MICH_QTR.csv")

dat <- merge(MICH=mich_q, CPI=cpi_q, all=FALSE); dat <- dat[, c("MICH","CPI")]
dat <- na.omit(dat); storage.mode(dat) <- "double"
cat(sprintf("Full quarterly sample: %d quarters, %s to %s\n",
            nrow(dat), as.character(start(dat)), as.character(end(dat))))

dir.create("results_2var_exp_q", showWarnings=FALSE)
samples <- list(
  Full           = dat,
  GreatInflation = window(dat, start=as.yearqtr("1967 Q2"), end=as.yearqtr("1982 Q4")),
  Core           = window(dat, start=as.yearqtr("1983 Q1"), end=as.yearqtr("2019 Q4")),
  COVID          = window(dat, start=as.yearqtr("2020 Q1")))
taus <- c(0.1,0.3,0.5,0.7,0.9)

dir_tab <- function(res){ O <- (res$CT[,,1,1]+res$CT[,,1,2])*100; diag(O)<-0
  list(TO=colSums(O), FROM=rowSums(O), NET=colSums(O)-rowSums(O), TCI=mean(colSums(O))) }

rows <- list()
for (nl in c(2,1)) for (nm in names(samples)){
  sub <- samples[[nm]]
  cat(sprintf("\n== %s, nlag=%d : %d quarters, %s to %s ==\n", nm, nl, nrow(sub),
              as.character(start(sub)), as.character(end(sub))))
  cat("tau   MICH_TO  MICH_FROM  MICH_NET   TCI   (NET>0 = expectations lead)\n")
  for (tau in taus){
    r <- R2ConnectednessQ2(sub, window.size=NULL, nlag=nl, tau=tau, shrink=TRUE, drop_own_lags=FALSE, progbar=FALSE)
    d <- dir_tab(r)
    rows[[length(rows)+1]] <- data.frame(nlag=nl, Sample=nm, Tau=tau,
      MICH_TO=round(d$TO["MICH"],2), MICH_FROM=round(d$FROM["MICH"],2),
      MICH_NET=round(d$NET["MICH"],2), TCI=round(d$TCI,2))
    cat(sprintf("%.1f  %8.2f  %9.2f  %8.2f  %6.2f\n", tau, d$TO["MICH"], d$FROM["MICH"], d$NET["MICH"], d$TCI))
  }
}
out <- do.call(rbind, rows)
write.csv(out, "results_2var_exp_q/MICH_vs_CPI_directional.csv", row.names=FALSE)

cat("\n>> MICH NET (lead +/follow -) by sample, nlag=2:\n")
w <- out[out$nlag==2, c("Sample","Tau","MICH_NET","TCI")]
print(w, row.names=FALSE)
cat("\n>> H6 episode contrast (nlag=2, tau=0.9): expectations should LEAD more in Great Inflation than COVID\n")
gi <- out[out$nlag==2 & out$Tau==0.9 & out$Sample=="GreatInflation","MICH_NET"]
cv <- out[out$nlag==2 & out$Tau==0.9 & out$Sample=="COVID","MICH_NET"]
cat(sprintf("   GreatInflation NET_MICH = %+.2f ;  COVID NET_MICH = %+.2f\n", gi, cv))
cat("\nSaved results_2var_exp_q/MICH_vs_CPI_directional.csv\nDone.\n")
