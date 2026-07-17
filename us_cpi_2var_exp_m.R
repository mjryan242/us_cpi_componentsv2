# ============================================================
# us_cpi_2var_exp_m.R -- MONTHLY inflation-expectations directional
#   connectedness (GitHub issue #1: the paper's tables h6/h6net are
#   quarterly; referees will want the monthly analogue in the appendix).
#
# System: Michigan 1-year-ahead expected inflation (FRED MICH, monthly,
#   from Jan 1978) paired with month-on-month CPI inflation (the paper's
#   working short-horizon / "mismatched" measure -- see DECISIONS.md).
# Samples: Full (1978-), late Great Inflation (1978-82; only ~60 months,
#   caveated in the table notes), Core (1983-2019), COVID (2020-).
# tau in {0.5, 0.7, 0.9}; nlag=2 (a two-month lag horizon).
#
# Outputs:
#   results_2var_exp_m/MICH_vs_CPI_directional_monthly.csv
#   paper/tables/h6_monthly.tex   (appendix fragment, h6 house style)
# ============================================================
pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R")

## ---- data ----
mich <- getSymbols("MICH", src="FRED", auto.assign=FALSE)
cat(sprintf("MICH monthly starts %s\n", as.character(start(mich))))   # verify the paper's "1978" claim
cpi  <- getSymbols("CPIAUCSL", src="FRED", auto.assign=FALSE)
cpi_mm <- 100 * (cpi / lag.xts(cpi, k = 1) - 1)

michz <- zoo(as.numeric(mich),   as.yearmon(index(mich)))
cpiz  <- zoo(as.numeric(cpi_mm), as.yearmon(index(cpi_mm)))
de <- na.omit(merge(MICH = michz, CPI = cpiz, all = FALSE))
de <- de[, c("MICH","CPI")]; storage.mode(de) <- "double"
cat(sprintf("system: %d months, %s to %s\n", nrow(de),
    as.character(start(de)), as.character(end(de))))

samp <- list(
  Full   = de,
  LateGI = window(de, end = as.yearmon("Dec 1982")),
  Core   = window(de, start = as.yearmon("Jan 1983"), end = as.yearmon("Dec 2019")),
  COVID  = window(de, start = as.yearmon("Jan 2020")))
lab <- c(Full="Full (1978--2026)", LateGI="Late Great Inflation (1978--82)",
         Core="Core (1983--2019)", COVID="COVID (2020--26)")

## ---- directional decomposition (same convention as make_tables_hyp.R h6) ----
dec <- function(d, tau) {
  r <- R2ConnectednessQ2(d, window.size=NULL, nlag=2, tau=tau, shrink=TRUE,
                         drop_own_lags=FALSE, progbar=FALSE)
  C <- r$CT[,,1,1]*100; L <- r$CT[,,1,2]*100     # [receiver, source]
  O <- C + L; diag(O) <- 0
  c(cCE=C["MICH","CPI"], lCE=L["MICH","CPI"],    # Inflation -> Expectations
    cEC=C["CPI","MICH"], lEC=L["CPI","MICH"],    # Expectations -> Inflation
    NET_MICH=unname(colSums(O)["MICH"] - rowSums(O)["MICH"]))
}

taus <- c(0.5, 0.7, 0.9)
rows <- list()
for (s in names(samp)) for (tau in taus) {
  m <- dec(samp[[s]], tau)
  rows[[length(rows)+1]] <- data.frame(Sample=s, Tau=tau, n=nrow(samp[[s]]),
    InflToExp_C=round(m["cCE"],1), InflToExp_L=round(m["lCE"],1),
    ExpToInfl_C=round(m["cEC"],1), ExpToInfl_L=round(m["lEC"],1),
    NET_MICH=round(m["NET_MICH"],1))
}
out <- do.call(rbind, rows); rownames(out) <- NULL
dir.create("results_2var_exp_m", showWarnings=FALSE)
write.csv(out, "results_2var_exp_m/MICH_vs_CPI_directional_monthly.csv", row.names=FALSE)
print(out, row.names=FALSE)

## ---- appendix table fragment (h6 house style + NET column) ----
f1 <- function(x) sprintf("%.1f", x)
blocks <- sapply(names(samp), function(s)
  paste(sapply(taus, function(tau){
    m <- out[out$Sample==s & out$Tau==tau, ]
    lead <- if (tau==0.5) lab[s] else ""
    sprintf("%s & %.1f & %s & %s & %s & %s & %s \\\\", lead, tau,
      f1(m$InflToExp_C), f1(m$InflToExp_L), f1(m$ExpToInfl_C), f1(m$ExpToInfl_L),
      sprintf("%+.1f", m$NET_MICH)) }), collapse="\n"))
body <- paste(blocks, collapse="\n\\addlinespace\n")
dir.create("paper/tables", showWarnings=FALSE, recursive=TRUE)
writeLines(paste0(
"\\begin{table}[!ht]\n\\centering\n",
"\\caption{Monthly robustness --- connectedness between inflation and expectations}\n",
"\\label{tab:h6monthly}\n",
"\\setlength{\\tabcolsep}{6pt}\\footnotesize\n\\begin{threeparttable}\n",
"\\begin{tabular}{l c c c c c c}\n\\toprule\n",
" & & \\multicolumn{2}{c}{Inflation $\\to$ Expectations} & \\multicolumn{2}{c}{Expectations $\\to$ Inflation} & \\\\\n",
"\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n",
"Sample & $\\tau$ & Contemp. & Lagged & Contemp. & Lagged & NET \\\\\n\\midrule\n",
body, "\n\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\\footnotesize\n",
"\\item \\textit{Notes:} Monthly analogue of Table~\\ref{tab:h6}: pairwise pseudo-quantile $R^2$ ",
"connectedness (\\%) in a two-variable system of the Michigan one-year inflation expectation ",
"(FRED \\texttt{MICH}, available monthly from January 1978) and month-on-month CPI inflation ",
"(CPIAUCSL), $n_{\\text{lag}}=2$ (a two-month lag horizon). NET is the net directional ",
"connectedness of expectations (TO $-$ FROM); $>0$ means expectations lead, $<0$ that they ",
"follow. Because the monthly survey begins in 1978, the Great-Inflation sample covers only ",
"its final years (1978--82, 60 months) and should be read with caution. Dependence, not ",
"identified causation.\n",
"\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n"),
"paper/tables/h6_monthly.tex")
cat("\nwrote paper/tables/h6_monthly.tex\nDone.\n")
