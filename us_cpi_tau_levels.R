# ============================================================
# us_cpi_tau_levels.R -- GitHub issue #4: "taus and levels".
#
# Problem: comparing the Great Inflation and COVID at the SAME tau compares
# DIFFERENT inflation levels, because the two episodes' inflation
# distributions differ. This script equalises the level instead:
#
#   Design 1 (common level -> episode-specific taus): take the headline
#     inflation rate at tau=0.9 of the COVID episode's distribution; find
#     the tau in the Great-Inflation distribution with the same headline
#     inflation rate; re-run the six-CPI episode comparison at those
#     episode-specific taus. Also the reverse (GI tau=0.9 level mapped
#     into the COVID distribution; if the level lies beyond the COVID
#     support this is reported, not silently capped).
#
#   Design 2 (common absolute threshold): fix 5% annualised headline
#     inflation; each episode's tau at that level; same comparison there.
#
# Interpretive note: tau in the connectedness estimator is a quantile of
# each series' own distribution; we map tau <-> inflation level via the
# within-episode empirical CDF of headline m/m CPI inflation (annualised
# for readability), as an interpretable common yardstick.
#
# The connectedness object compared is the one-in/one-out marginal vs the
# 1983-2019 core (as in us_cpi_6cpi_era.R), evaluated at each episode's
# level-matched tau.
#
# Outputs:
#   results_6cpi_era/tau_levels.csv
#   paper/tables/tau_levels.tex   (appendix fragment)
# ============================================================
pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R"); W <- 0.60

## ---- six-CPI panel (as in us_cpi_6cpi_era.R) ----
codes <- c(Food="CPIUFDSL", Gasoline="CUSR0000SETB01", HHEnergy="CUSR0000SEHF",
           CoreGoods="CUSR0000SACL1E", Shelter="CUSR0000SAH1", CoreServices="CUSR0000SASLE")
sl <- list(); for (nm in names(codes)) { s <- getSymbols(codes[[nm]], src="FRED", auto.assign=FALSE); colnames(s)<-nm; sl[[nm]]<-s }
idx <- do.call(merge, sl)[, names(codes)]
pcz <- zoo(coredata(100*(idx/lag.xts(idx,k=1)-1)), as.yearmon(index(idx))); colnames(pcz)<-names(codes)
exS <- (pcz[,"CoreServices"]-W*pcz[,"Shelter"])/(1-W)
dat <- cbind(Food=pcz[,"Food"], Gasoline=pcz[,"Gasoline"], HHEnergy=pcz[,"HHEnergy"],
             CoreGoods=pcz[,"CoreGoods"], Shelter=pcz[,"Shelter"], CoreServ_xS=exS)
dat <- na.omit(dat); storage.mode(dat) <- "double"

## ---- headline inflation within each episode (annualised m/m, %) ----
cpi <- getSymbols("CPIAUCSL", src="FRED", auto.assign=FALSE)
head_mm <- zoo(as.numeric(100*((cpi/lag.xts(cpi,k=1))^12 - 1)), as.yearmon(index(cpi)))
gi_infl <- window(head_mm, start=as.yearmon("Jan 1967"), end=as.yearmon("Dec 1982"))
cv_infl <- window(head_mm, start=as.yearmon("Jan 2020"))
lev <- function(x, tau) as.numeric(quantile(x, tau, na.rm=TRUE, type=8))
tau_of <- function(x, level) mean(x <= level, na.rm=TRUE)   # empirical CDF

## ---- Design 1: common level -> episode taus ----
lev_cv09 <- lev(cv_infl, 0.9)          # COVID tau=0.9 headline level
tau_gi_star <- tau_of(gi_infl, lev_cv09)
lev_gi09 <- lev(gi_infl, 0.9)          # GI tau=0.9 headline level
tau_cv_star <- tau_of(cv_infl, lev_gi09)
cat(sprintf("COVID tau=0.9 level: %.2f%% ann. -> GI tau* = %.3f\n", lev_cv09, tau_gi_star))
cat(sprintf("GI    tau=0.9 level: %.2f%% ann. -> COVID tau* = %.3f%s\n", lev_gi09, tau_cv_star,
            if (tau_cv_star >= 0.995) "  [beyond COVID support]" else ""))

## ---- Design 2: common absolute threshold (5% annualised) ----
THR <- 5
tau_gi_thr <- tau_of(gi_infl, THR); tau_cv_thr <- tau_of(cv_infl, THR)
cat(sprintf("5%% ann. threshold -> GI tau = %.3f, COVID tau = %.3f\n", tau_gi_thr, tau_cv_thr))

## ---- one-in/one-out marginal at a given tau (as in us_cpi_6cpi_era.R) ----
mct <- function(res){ C<-res$CT[,,1,1]*100; L<-res$CT[,,1,2]*100; O<-C+L; od<-function(M){diag(M)<-0;M}
  c(Overall=mean(colSums(od(O))), Contemp=mean(colSums(od(C))), Lagged=mean(colSums(od(L)))) }
est <- function(d, tau) mct(R2ConnectednessQ2(d,window.size=NULL,nlag=2,tau=tau,shrink=TRUE,drop_own_lags=FALSE,progbar=FALSE))
core  <- window(dat, start=as.yearmon("Jan 1983"), end=as.yearmon("Dec 2019"))
s_gi  <- window(dat, end=as.yearmon("Dec 2019"))       # GI in, COVID out
s_cv  <- window(dat, start=as.yearmon("Jan 1983"))     # COVID in, GI out
marg <- function(smp, tau) { d <- est(smp, tau) - est(core, tau)
  data.frame(d_Overall=round(d["Overall"],2), d_Contemp=round(d["Contemp"],2), d_Lagged=round(d["Lagged"],2)) }

## clamp taus to an estimable range; record the raw value alongside
cl <- function(tau) min(max(tau, 0.05), 0.95)
runs <- list(
  list(Design="1: common level (COVID 0.9 ref)", Episode="GreatInflation", TauRaw=tau_gi_star, Level=lev_cv09, smp=s_gi),
  list(Design="1: common level (COVID 0.9 ref)", Episode="COVID",          TauRaw=0.9,         Level=lev_cv09, smp=s_cv),
  list(Design="1: common level (GI 0.9 ref)",    Episode="GreatInflation", TauRaw=0.9,         Level=lev_gi09, smp=s_gi),
  list(Design="1: common level (GI 0.9 ref)",    Episode="COVID",          TauRaw=tau_cv_star, Level=lev_gi09, smp=s_cv),
  list(Design="2: 5% annualised threshold",      Episode="GreatInflation", TauRaw=tau_gi_thr,  Level=THR,      smp=s_gi),
  list(Design="2: 5% annualised threshold",      Episode="COVID",          TauRaw=tau_cv_thr,  Level=THR,      smp=s_cv))
out <- NULL
for (r in runs) {
  tau_used <- cl(r$TauRaw)
  m <- marg(r$smp, tau_used)
  out <- rbind(out, cbind(data.frame(Design=r$Design, Episode=r$Episode,
    Level_ann_pct=round(r$Level,2), Tau_raw=round(r$TauRaw,3), Tau_used=round(tau_used,3)), m))
  cat(sprintf("%-32s %-15s tau=%.3f  dOverall %+6.2f  dContemp %+6.2f  dLagged %+6.2f\n",
      r$Design, r$Episode, tau_used, m$d_Overall, m$d_Contemp, m$d_Lagged))
}
rownames(out) <- NULL
dir.create("results_6cpi_era", showWarnings=FALSE)
write.csv(out, "results_6cpi_era/tau_levels.csv", row.names=FALSE)

## ---- appendix table fragment ----
fr <- function(i) { r <- out[i,]
  sprintf("%s & %.2f & %.3f%s & %+.1f & %+.1f & %+.1f \\\\",
    ifelse(r$Episode=="GreatInflation","Great Inflation","COVID"),
    r$Level_ann_pct, r$Tau_used,
    ifelse(abs(r$Tau_raw - r$Tau_used) > 1e-9, "$^{\\dagger}$", ""),
    r$d_Overall, r$d_Contemp, r$d_Lagged) }
pan <- function(rows_idx, title) paste0(
  "\\multicolumn{6}{l}{\\textit{", title, "}}\\\\\n",
  paste(sapply(rows_idx, fr), collapse="\n"), "\n\\addlinespace\n")
tex <- paste0(
"\\begin{table}[!ht]\n\\centering\n",
"\\caption{Level-matched episode comparison: same inflation rate, episode-specific $\\tau$}\n",
"\\label{tab:taulevels}\n",
"\\setlength{\\tabcolsep}{6pt}\\footnotesize\n\\begin{threeparttable}\n",
"\\begin{tabular}{l c c c c c}\n\\toprule\n",
"Episode & Level (\\% ann.) & $\\tau$ used & $\\Delta$Overall & $\\Delta$Contemp. & $\\Delta$Lagged \\\\\n\\midrule\n",
pan(1:2, "Panel A. Common level: COVID $\\tau=0.9$ headline rate as reference"),
pan(3:4, "Panel B. Common level: Great-Inflation $\\tau=0.9$ headline rate as reference"),
pan(5:6, "Panel C. Common absolute threshold: 5\\% annualised headline inflation"),
"\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\\footnotesize\n",
"\\item \\textit{Notes:} Comparing the two episodes at the same $\\tau$ compares different inflation ",
"levels, because their inflation distributions differ. Here the comparison is level-matched: within ",
"each episode, $\\tau$ is chosen so that the corresponding quantile of annualised month-on-month ",
"headline CPI inflation equals the common reference level shown. $\\Delta$ columns are the episode's ",
"one-in/one-out marginal connectedness relative to the 1983--2019 core (as in Table~\\ref{tab:h3}), ",
"evaluated at that episode-specific $\\tau$; six-CPI system, monthly, $n_{\\text{lag}}=2$. ",
"$^{\\dagger}$~The implied $\\tau$ lies outside $[0.05,0.95]$ (the reference level sits at or beyond ",
"the edge of that episode's inflation distribution) and is clamped to the boundary; the comparison ",
"at that reference should be read as one-sided.\n",
"\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n")
writeLines(tex, "paper/tables/tau_levels.tex")
cat("\nSaved results_6cpi_era/tau_levels.csv and paper/tables/tau_levels.tex\nDone.\n")
