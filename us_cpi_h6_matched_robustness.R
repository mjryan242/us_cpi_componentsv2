# ============================================================
# us_cpi_h6_matched_robustness.R -- matched-horizon robustness for table h6 (H4).
#
# Table h6 pairs MICH_QTR (Michigan 1-YEAR-ahead expected inflation) with
# quarter-on-quarter CPI inflation -- a horizon mismatch. This script re-runs
# the identical directional decomposition with YoY (4-quarter) CPI inflation,
# matched to the survey's one-year horizon, and prints both side by side.
#
# QoQ columns replicate make_tables_hyp.R's h6 block (sanity check against the
# published table: Full tau=0.5 should give ~24.2/23.7/27.4/25.x; Great
# Inflation tau=0.9 lagged Infl->Exp ~50.7 vs Exp->Infl ~27.5, up to FRED
# vintage changes). NET_MICH replicates the h6net logic (positive = expectations
# lead / transmit; negative = expectations follow).
#
# Output: results_2var_exp_q/h6_matched_robustness.csv  (new file only)
# ============================================================

pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R")
source("quarterly_utils.R")

## ---- data: quarterly-average CPI level -> QoQ and YoY inflation ----
cpi    <- getSymbols("CPIAUCSL", src="FRED", auto.assign=FALSE)
qcpi_x <- apply.quarterly(cpi, function(v) mean(v, na.rm = TRUE))   # period-average level (to_qpc convention)
cpi_qoq <- zoo(as.numeric(100 * (qcpi_x / lag.xts(qcpi_x, k = 1) - 1)), as.yearqtr(index(qcpi_x)))
cpi_yoy <- zoo(as.numeric(100 * (qcpi_x / lag.xts(qcpi_x, k = 4) - 1)), as.yearqtr(index(qcpi_x)))

mich_q <- read_mich_qtr("MICH_QTR.csv")   # 1-year-ahead expectation, %, as-is

build <- function(cpiz) {
  d <- na.omit(merge(MICH = mich_q, CPI = cpiz, all = FALSE))
  d <- d[, c("MICH", "CPI")]; storage.mode(d) <- "double"; d
}
de <- list(QoQ = build(cpi_qoq), YoY = build(cpi_yoy))
for (m in names(de)) cat(sprintf("%s system: %d quarters, %s to %s\n",
  m, nrow(de[[m]]), as.character(start(de[[m]])), as.character(end(de[[m]]))))

## ---- samples exactly as in make_tables_hyp.R ----
make_samples <- function(d) list(
  Full           = d,
  GreatInflation = window(d, start = as.yearqtr("1967 Q2"), end = as.yearqtr("1982 Q4")),
  Core           = window(d, start = as.yearqtr("1983 Q1"), end = as.yearqtr("2019 Q4")),
  COVID          = window(d, start = as.yearqtr("2020 Q1")))

## ---- directional decomposition (same as h6) + NET (same as h6net) ----
dec <- function(d, tau) {
  r <- R2ConnectednessQ2(d, window.size = NULL, nlag = 2, tau = tau, shrink = TRUE, progbar = FALSE)
  C <- r$CT[,,1,1] * 100; L <- r$CT[,,1,2] * 100      # [receiver, source]
  O <- C + L; diag(O) <- 0
  c(cCE = C["MICH","CPI"], lCE = L["MICH","CPI"],     # Inflation -> Expectations
    cEC = C["CPI","MICH"], lEC = L["CPI","MICH"],     # Expectations -> Inflation
    NET_MICH = unname(colSums(O)["MICH"] - rowSums(O)["MICH"]))
}

taus <- c(0.5, 0.7, 0.9)
rows <- list()
for (meas in names(de)) {
  smp <- make_samples(de[[meas]])
  for (s in names(smp)) for (tau in taus) {
    m <- dec(smp[[s]], tau)
    rows[[length(rows) + 1]] <- data.frame(measure = meas, sample = s, tau = tau,
      InflToExp_C = round(m["cCE"], 1), InflToExp_L = round(m["lCE"], 1),
      ExpToInfl_C = round(m["cEC"], 1), ExpToInfl_L = round(m["lEC"], 1),
      NET_MICH = round(m["NET_MICH"], 1))
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
dir.create("results_2var_exp_q", showWarnings = FALSE)
write.csv(out, "results_2var_exp_q/h6_matched_robustness.csv", row.names = FALSE)

cat("\n================ QoQ (current table h6) ================\n")
print(out[out$measure == "QoQ", -1], row.names = FALSE)
cat("\n================ YoY (matched horizon) =================\n")
print(out[out$measure == "YoY", -1], row.names = FALSE)
cat("\nKey H4 checks:\n")
gq <- out[out$measure=="QoQ" & out$sample=="GreatInflation" & out$tau==0.9, ]
gy <- out[out$measure=="YoY" & out$sample=="GreatInflation" & out$tau==0.9, ]
cat(sprintf("  GI tau=0.9 lagged Infl->Exp vs Exp->Infl:  QoQ %.1f vs %.1f | YoY %.1f vs %.1f\n",
            gq$InflToExp_L, gq$ExpToInfl_L, gy$InflToExp_L, gy$ExpToInfl_L))
for (s in c("GreatInflation","Core","COVID")) {
  nq <- out[out$measure=="QoQ" & out$sample==s & out$tau==0.9, "NET_MICH"]
  ny <- out[out$measure=="YoY" & out$sample==s & out$tau==0.9, "NET_MICH"]
  cat(sprintf("  NET_MICH tau=0.9 %-15s QoQ %6.1f | YoY %6.1f\n", s, nq, ny))
}
cat("\nSaved results_2var_exp_q/h6_matched_robustness.csv\nDone.\n")
