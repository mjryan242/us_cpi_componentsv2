# ============================================================
# us_cpi_2var_exp_q_marginal.R -- one-in/one-out episode marginals for the
#   quarterly inflation-expectations system (the exact analogue of the
#   six-CPI episode mechanism in us_cpi_6cpi_era.R / Table h3).
#
#   GI marginal    = (1960-2019, COVID out) - (1983-2019 core)
#   COVID marginal = (1983-2026, GI out)    - (1983-2019 core)
#
# System: MICH_QTR (Michigan 1-yr expectation, quarterly from 1960) and
# quarter-on-quarter CPI inflation, as in table h6. tau in {0.7, 0.9}.
#
# Outputs:
#   results_2var_exp_q/episode_marginal_decomposition.csv
#   paper/tables/h6_marginal.tex   (appendix fragment)
# ============================================================
pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R"); source("quarterly_utils.R")

cpi_q  <- to_qpc(getSymbols("CPIAUCSL", src="FRED", auto.assign=FALSE)); colnames(cpi_q) <- "CPI"
mich_q <- read_mich_qtr("MICH_QTR.csv")
de <- merge(MICH=mich_q, CPI=cpi_q, all=FALSE); de <- na.omit(de[, c("MICH","CPI")]); storage.mode(de) <- "double"
cat(sprintf("system: %d quarters, %s to %s\n", nrow(de), as.character(start(de)), as.character(end(de))))

mct <- function(res){ C<-res$CT[,,1,1]*100; L<-res$CT[,,1,2]*100; O<-C+L; od<-function(M){diag(M)<-0;M}
  c(Overall=mean(colSums(od(O))), Contemp=mean(colSums(od(C))), Lagged=mean(colSums(od(L)))) }
est <- function(d, tau) mct(R2ConnectednessQ2(d,window.size=NULL,nlag=2,tau=tau,shrink=TRUE,drop_own_lags=FALSE,progbar=FALSE))

out <- NULL
for (tau in c(0.7, 0.9)) {
  c0 <- est(window(de, start=as.yearqtr("1983 Q1"), end=as.yearqtr("2019 Q4")), tau)
  gi <- est(window(de, end=as.yearqtr("2019 Q4")), tau)              # GI in, COVID out (1960-2019)
  cv <- est(window(de, start=as.yearqtr("1983 Q1")), tau)            # COVID in, GI out (1983-2026)
  mk <- function(x){ d<-x-c0; data.frame(d_Overall=round(d["Overall"],2), d_Contemp=round(d["Contemp"],2),
    d_Lagged=round(d["Lagged"],2), Lagged_share_pct=round(100*d["Lagged"]/d["Overall"],1)) }
  out <- rbind(out,
    cbind(Tau=tau, Episode="GreatInflation (1960-82 in)", mk(gi)),
    cbind(Tau=tau, Episode="COVID (2020-26 in)",          mk(cv)))
  cat(sprintf("core (1983-2019) tau%.1f: Overall %.1f, Contemp %.1f, Lagged %.1f\n",
      tau, c0["Overall"],c0["Contemp"],c0["Lagged"]))
}
rownames(out) <- NULL
write.csv(out, "results_2var_exp_q/episode_marginal_decomposition.csv", row.names=FALSE)
print(out, row.names=FALSE)

## ---- appendix fragment (h3 house style) ----
erow <- function(i){ r<-out[i,]; sprintf("%s & %+.1f & %+.1f & %+.1f & %.1f \\\\",
  ifelse(grepl("GreatInflation", r$Episode), "Great Inflation (pre-1983 in)", "COVID (2020-- in)"),
  r$d_Overall, r$d_Contemp, r$d_Lagged, r$Lagged_share_pct) }
blk <- function(idx, tau) paste0(
  "\\multicolumn{5}{l}{\\textit{$\\tau=", tau, "$}}\\\\\n",
  paste(sapply(idx, erow), collapse="\n"), "\n\\addlinespace\n")
tex <- paste0(
"\\begin{table}[!ht]\n\\centering\n",
"\\caption{Episode marginals for the inflation--expectations system (one-in/one-out)}\n",
"\\label{tab:h6marginal}\n",
"\\setlength{\\tabcolsep}{6pt}\\footnotesize\n\\begin{threeparttable}\n",
"\\begin{tabular}{l c c c c}\n\\toprule\n",
"Episode & $\\Delta$Overall & $\\Delta$Contemp. & $\\Delta$Lagged & Lagged \\% \\\\\n\\midrule\n",
blk(1:2, "0.7"), blk(3:4, "0.9"),
"\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\\footnotesize\n",
"\\item \\textit{Notes:} One-in/one-out episode marginals for the quarterly two-variable ",
"system of Michigan inflation expectations (\\texttt{MICH\\_QTR}) and quarter-on-quarter CPI ",
"inflation, constructed exactly as in Table~\\ref{tab:h3}: each row is total connectedness ",
"estimated with that episode included (and the other excluded) minus the calm 1983--2019 core, ",
"split into contemporaneous and lagged; $n_{\\text{lag}}=2$ (a six-month lag horizon). The ",
"pre-1983 sample runs from 1960, so the ``Great Inflation'' marginal includes the 1960s ",
"run-up. Dependence, not identified causation.\n",
"\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n")
writeLines(tex, "paper/tables/h6_marginal.tex")
cat("\nSaved results_2var_exp_q/episode_marginal_decomposition.csv and paper/tables/h6_marginal.tex\nDone.\n")
