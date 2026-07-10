# ============================================================
# us_cpi_wage_spiral_q.R -- QUARTERLY wage-price spiral test.
#   2-variable connectedness between WAGES (AHETPI, avg hourly earnings)
#   and overall CPI (CPIAUCSL), both aggregated to quarterly QoQ %
#   (level -> quarterly period-average -> QoQ %; see quarterly_utils.R),
#   estimated separately on Great Inflation vs COVID windows.
#
#   A spiral is LAGGED feedback, so the quarterly horizon is the point:
#   nlag=2 = 6 months, nlag=4 = 1 year.  Spiral metric = lagged
#   Wages->CPI + lagged CPI->Wages.
#
#   Estimators: Genizi (R2ConnectednessQ2) and LASSO predictive
#   (R2Q_lasso_CT).  AHETPI starts 1964 -> full Great-Inflation cover.
# Output: results_wage_spiral_q/
# ============================================================

pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R"); source("R2Q_lasso_dir.R")
source("quarterly_utils.R")

cpi   <- getSymbols("CPIAUCSL", src="FRED", auto.assign=FALSE)
cpi_q <- to_qpc(cpi); colnames(cpi_q) <- "CPI"
w_raw <- read.csv("AHETPI.csv", stringsAsFactors=FALSE); w_raw$observation_date <- as.Date(w_raw$observation_date)
wz    <- zoo(w_raw$AHETPI, w_raw$observation_date)      # Date-indexed monthly wage level
w_q   <- to_qpc(wz); colnames(w_q) <- "Wages"
dat <- merge(Wages=w_q, CPI=cpi_q, all=FALSE); dat <- dat[, c("Wages","CPI")]; dat <- na.omit(dat)
cat(sprintf("Full quarterly sample: %d quarters, %s to %s\n\n", nrow(dat),
            as.character(start(dat)), as.character(end(dat))))
dir.create("results_wage_spiral_q", showWarnings=FALSE)

windows <- list(
  Full           = dat,
  GreatInflation = window(dat, start=as.yearqtr("1967 Q2"), end=as.yearqtr("1982 Q4")),
  Core           = window(dat, start=as.yearqtr("1983 Q1"), end=as.yearqtr("2019 Q4")),
  COVID          = window(dat, start=as.yearqtr("2020 Q1")))

metrics <- function(C, L){ O <- C+L; od <- function(M){diag(M)<-0;M}
  list(TCI=mean(colSums(od(O))), TCI_C=mean(colSums(od(C))), TCI_L=mean(colSums(od(L))),
       w2p_L=L["CPI","Wages"], p2w_L=L["Wages","CPI"],          # lagged wages->CPI ; CPI->wages
       spiral_L=L["CPI","Wages"]+L["Wages","CPI"]) }            # total lagged bidirectional feedback

genizi <- function(d, tau, nl){ r <- R2ConnectednessQ2(d, window.size=NULL, nlag=nl, tau=tau,
    shrink=TRUE, drop_own_lags=FALSE, progbar=FALSE); metrics(r$CT[,,1,1]*100, r$CT[,,1,2]*100) }
lasso <- function(d, tau, nl){ ct <- R2Q_lasso_CT(as.matrix(d), nlag=nl, tau=tau); metrics(ct$C, ct$L) }

rows <- list()
for (tau in c(0.5, 0.7, 0.9)) for (nm in names(windows)) {
  d <- windows[[nm]]
  for (est in list(c("Genizi","2"), c("Genizi","4"), c("LASSO","2"), c("LASSO","4"))) {
    lab <- est[1]; nl <- as.integer(est[2])
    m <- tryCatch(if (lab=="Genizi") genizi(d, tau, nl) else lasso(d, tau, nl), error=function(e) NULL)
    if (is.null(m)) next
    rows[[length(rows)+1]] <- data.frame(tau=tau, window=nm, estimator=lab, nlag=nl, n=nrow(d),
      TCI=round(m$TCI,2), TCI_C=round(m$TCI_C,2), TCI_L=round(m$TCI_L,2),
      Wages_to_CPI_L=round(m$w2p_L,2), CPI_to_Wages_L=round(m$p2w_L,2), spiral_L=round(m$spiral_L,2))
  }
}
out <- do.call(rbind, rows)
write.csv(out, "results_wage_spiral_q/wage_cpi_connectedness.csv", row.names=FALSE)

cat("=== Spiral metric (total LAGGED wage<->CPI feedback) by window ===\n")
for (tau in c(0.5,0.7,0.9)) { cat(sprintf("\n-- tau=%.1f --\n", tau))
  s <- out[out$tau==tau, ]
  cat("window          est     nlag   n   TCI  TCI_L  W->CPI(L) CPI->W(L) spiral_L\n")
  for (i in seq_len(nrow(s))) cat(sprintf("%-14s %-7s %d   %3d  %5.1f %5.1f   %6.2f   %6.2f   %6.2f\n",
    s$window[i], s$estimator[i], s$nlag[i], s$n[i], s$TCI[i], s$TCI_L[i],
    s$Wages_to_CPI_L[i], s$CPI_to_Wages_L[i], s$spiral_L[i])) }

cat("\n>> Hypothesis (spiral stronger in Great Inflation than COVID), tau=0.5 lagged feedback:\n")
for (est in list(c("Genizi",2),c("Genizi",4),c("LASSO",2),c("LASSO",4))) {
  gi <- out$spiral_L[out$tau==0.5 & out$window=="GreatInflation" & out$estimator==est[1] & out$nlag==est[2]]
  cv <- out$spiral_L[out$tau==0.5 & out$window=="COVID"          & out$estimator==est[1] & out$nlag==est[2]]
  if (length(gi)&&length(cv)) cat(sprintf("  %-7s nlag=%s : GreatInflation=%.2f  COVID=%.2f  -> %s\n", est[1], est[2], gi, cv,
              ifelse(gi>cv,"GI > COVID (supports)","GI <= COVID (does not)")))
}
cat("\nSaved results_wage_spiral_q/wage_cpi_connectedness.csv\nDone.\n")
