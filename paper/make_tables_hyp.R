# ============================================================
# make_tables_hyp.R -- ONE publication-ready table per hypothesis
#   (new H1--H6 numbering), JMCB house style (booktabs/threeparttable).
#   Reads the result CSVs; recomputes only the H5 sequencing piece.
#   Writes paper/tables/h1.tex ... h6.tex.
#
#   RUN FROM the us_cpi_components/ directory:
#     Rscript paper/make_tables_hyp.R
# ============================================================
pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R"); source("R2Q_lasso_dir.R"); source("quarterly_utils.R")
W <- 0.60
TD <- "paper/tables"; dir.create(TD, showWarnings=FALSE, recursive=TRUE)
wr <- function(name, txt){ writeLines(txt, file.path(TD, paste0(name,".tex"))); cat("wrote", name, "\n") }
f1 <- function(x) sprintf("%.1f", x); f2 <- function(x) sprintf("%.2f", x)
ci <- function(b,l,h) sprintf("%.1f [%.1f, %.1f]", b, l, h)
ratio <- function(d) d$Overall[d$Tau==0.9]/d$Overall[d$Tau==0.5]

## ============================================================
## H1 -- CPI broadening
## ============================================================
cat("building h1\n")
taus <- c(0.1,0.3,0.5,0.7,0.9)
tm <- read.csv("results_6cpi/TCI_by_quantile.csv")     # monthly six-CPI
tq <- read.csv("results_6cpi_q/TCI_by_quantile.csv")   # quarterly six-CPI
prow <- function(df) paste(sapply(taus, function(t){ r <- df[df$Tau==t, ]
  sprintf("%s & %s & %s & %s \\\\", f1(t), f1(r$Overall), f1(r$Contemporaneous), f1(r$Lagged)) }), collapse="\n")
bodyA <- prow(tm); bodyB <- prow(tq); rb_m <- ratio(tm); rb_q <- ratio(tq)
wr("h1", paste0(
"\\begin{table}[!ht]\n\\centering\n\\caption{H1 --- Inflation broadens into the high tail}\n\\label{tab:h1}\n",
"\\setlength{\\tabcolsep}{6pt}\\footnotesize\n\\begin{threeparttable}\n",
"\\begin{tabular}{c c c c}\n",
"\\multicolumn{4}{l}{\\textit{Panel A. Monthly total connectedness by quantile}}\\\\\n\\toprule\n",
"$\\tau$ & Overall & Contemporaneous & Lagged \\\\\n\\midrule\n", bodyA,
"\n\\bottomrule\n\\end{tabular}\n\n\\vspace{4pt}\n",
"\\begin{tabular}{c c c c}\n",
"\\multicolumn{4}{l}{\\textit{Panel B. Quarterly total connectedness by quantile}}\\\\\n\\toprule\n",
"$\\tau$ & Overall & Contemporaneous & Lagged \\\\\n\\midrule\n", bodyB,
"\n\\bottomrule\n\\end{tabular}\n",
"\\begin{tablenotes}[flushleft]\\footnotesize\n",
"\\item \\textit{Notes:} Pseudo-quantile $R^2$ (Genizi/nearPD) connectedness (\\%) among the six CPI subcomponents; ",
"Overall $=$ Contemporaneous $+$ Lagged, $n_{\\text{lag}}=2$. Panel A is monthly 1967--2026 (broadening ratio ",
"$\\text{TCI}(0.9)/\\text{TCI}(0.5)=", f2(rb_m), "\\times$); Panel B quarterly 1967Q2--2026Q2 (ratio ", f2(rb_q),
"$\\times$). Both frequencies show the tail rise, increasingly driven by the lagged component. Point estimates; ",
"bootstrap confidence intervals for the six-CPI system are to be added.\n",
"\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n"))

## ============================================================
## H2 -- Wage-price interaction
## ============================================================
cat("building h2\n")
sp  <- read.csv("results_wage_spiral/wage_cpi_connectedness.csv")
spq <- read.csv("results_wage_spiral_q/wage_cpi_connectedness.csv")
tausH2 <- c(0.5,0.7,0.9)
vrow <- function(df, win, label, metric){ s <- df[df$window==win & df$estimator=="Genizi" & df$nlag==2, ]
  v <- sapply(tausH2, function(t) s[s$tau==t, metric])
  sprintf("%s & %s & %s & %s \\\\", label, f1(v[1]), f1(v[2]), f1(v[3])) }
panel3 <- function(df, win) paste(
  vrow(df,win,"Wage--CPI total connectedness (TCI)","TCI"),
  vrow(df,win,"\\quad lagged Wages$\\to$CPI","Wages_to_CPI_L"),
  vrow(df,win,"\\quad lagged CPI$\\to$Wages","CPI_to_Wages_L"), sep="\n")
hdr <- "$\\tau=0.5$ & $\\tau=0.7$ & $\\tau=0.9$"
# Panel D: 1964-2019, level-matched taus (from us_cpi_wage_taulevel.R)
tl <- read.csv("results_wage_spiral/tau_level_matched.csv")
tl7 <- tl[tl$post_tau==0.7, ]; tl9 <- tl[tl$post_tau==0.9, ]
panelD <- paste(
  sprintf("Wage--CPI total connectedness (TCI) & %s & %s \\\\", f1(tl7$TCI), f1(tl9$TCI)),
  sprintf("\\quad lagged Wages$\\to$CPI & %s & %s \\\\", f1(tl7$Wages_to_CPI_L), f1(tl9$Wages_to_CPI_L)),
  sprintf("\\quad lagged CPI$\\to$Wages & %s & %s \\\\", f1(tl7$CPI_to_Wages_L), f1(tl9$CPI_to_Wages_L)), sep="\n")
hdrD <- sprintf("Matched to $\\tau{=}0.7$ & Matched to $\\tau{=}0.9$")
wr("h2", paste0(
"\\begin{table}[!ht]\n\\centering\n\\caption{H2 --- Wage--price connectedness is stronger when inflation is high}\n\\label{tab:h2}\n",
"\\setlength{\\tabcolsep}{6pt}\\footnotesize\n\\begin{threeparttable}\n",
"\\begin{tabular}{l c c c}\n",
"\\multicolumn{4}{l}{\\textit{Panel A. Two-variable wage--CPI system (monthly, full sample), by quantile}}\\\\\n\\toprule\n",
" & ", hdr, " \\\\\n\\midrule\n", panel3(sp,"Full"),
"\n\\bottomrule\n\\end{tabular}\n\n\\vspace{4pt}\n",
"\\begin{tabular}{l c c c}\n",
"\\multicolumn{4}{l}{\\textit{Panel B. Quarterly (full sample), by quantile}}\\\\\n\\toprule\n",
" & ", hdr, " \\\\\n\\midrule\n", panel3(spq,"Full"),
"\n\\bottomrule\n\\end{tabular}\n\n\\vspace{4pt}\n",
"\\begin{tabular}{l c c c}\n",
"\\multicolumn{4}{l}{\\textit{Panel C. By subperiod (monthly), by quantile}}\\\\\n\\toprule\n",
" & ", hdr, " \\\\\n\\midrule\n",
"\\multicolumn{4}{l}{\\textit{\\;1964--2019 (with Great Inflation)}}\\\\\n", panel3(sp,"Pre_1964_2019"), "\n\\addlinespace\n",
"\\multicolumn{4}{l}{\\textit{\\;1983--2026 (post-Great-Inflation)}}\\\\\n", panel3(sp,"Post_1983_2026"),
"\n\\bottomrule\n\\end{tabular}\n\n\\vspace{4pt}\n",
"\\begin{tabular}{l c c}\n",
"\\multicolumn{3}{l}{\\textit{Panel D. 1964--2019 (with Great Inflation), level-matched $\\tau$ (monthly)}}\\\\\n\\toprule\n",
" & ", hdrD, " \\\\\n\\midrule\n", panelD,
"\n\\bottomrule\n\\end{tabular}\n",
"\\begin{tablenotes}[flushleft]\\footnotesize\n",
"\\item \\textit{Notes:} Two-variable system of wage growth (AHETPI) and overall CPI (CPIAUCSL); pseudo-quantile ",
"$R^2$ connectedness (\\%), Genizi estimator, $n_{\\text{lag}}=2$. Wage--price connectedness rises sharply into ",
"the high-inflation tail at all frequencies (Panels A--B), and is far stronger in the Great-Inflation-inclusive ",
"subperiod (1964--2019) than post-1983 (Panel C). In the high tail the lagged components dominate. Quarterly ",
"($n_{\\text{lag}}=2$) spans a six-month horizon. Panel~D addresses the fact that comparing the two subperiods ",
"at the same $\\tau$ compares different inflation \\emph{levels}: it reports the 1964--2019 connectedness at the ",
"quantiles ($\\tau=", f2(tl7$tau_matched), "$ and $", f2(tl9$tau_matched), "$) whose month-on-month CPI-inflation ",
"level matches that of the post-1983 sample at $\\tau=0.7$ and $\\tau=0.9$ (about ", sprintf("%.1f",tl7$level_ann),

" and ", sprintf("%.1f",tl9$level_ann), " per cent annualised). Even level-matched, 1964--2019 connectedness ",
"(", f1(tl7$TCI), " and ", f1(tl9$TCI), ") remains far above the post-1983 values in Panel~C (", f1(sp$TCI[sp$window=="Post_1983_2026" & sp$estimator=="Genizi" & sp$nlag==2 & sp$tau==0.7]),
" and ", f1(sp$TCI[sp$window=="Post_1983_2026" & sp$estimator=="Genizi" & sp$nlag==2 & sp$tau==0.9]),
"), so the gap is not a level artefact; the level-matched spillovers also stay roughly symmetric, the balanced ",
"feedback characteristic of a spiral.\n",
"\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n"))

## ============================================================
## H3 -- Episode mechanism
## ============================================================
cat("building h3\n")
mm <- read.csv("results_6cpi_era/episode_marginal_decomposition.csv")   # same-quantile monthly marginals
tl <- read.csv("results_6cpi_era/tau_levels.csv")                        # level-matched marginals
# --- Panel A: same-quantile monthly marginals (no Lagged% col, no quarterly panel) ---
erow <- function(df, tau, ep) { r <- df[df$Tau==tau & grepl(ep, df$Episode), ]
  sprintf("%s & %s & %s & %s \\\\", gsub("_","\\\\_",r$Episode),
    sprintf("%+.1f",r$d_Overall), sprintf("%+.1f",r$d_Contemp), sprintf("%+.1f",r$d_Lagged)) }
tblock <- function(df, tau) paste0(
  "\\multicolumn{4}{l}{\\textit{\\;$\\tau=", tau, "$}}\\\\\n",
  erow(df, tau, "GreatInflation"), "\n", erow(df, tau, "COVID"))
epanel <- function(df) paste0(tblock(df, 0.7), "\n\\addlinespace\n", tblock(df, 0.9))
# --- level-matched panels (fold in tau_levels; drop its GI-reference design) ---
lrow <- function(design, ep) { r <- tl[tl$Design==design & tl$Episode==ep, ]
  lab <- if (ep=="GreatInflation") "Great Inflation" else "COVID"
  clamp <- if (abs(r$Tau_raw - r$Tau_used) > 1e-9) "$^{\\dagger}$" else ""
  sprintf("%s ($\\tau=%.2f$)%s & %s & %s & %s \\\\", lab, r$Tau_used, clamp,
    sprintf("%+.1f",r$d_Overall), sprintf("%+.1f",r$d_Contemp), sprintf("%+.1f",r$d_Lagged)) }
lpanel <- function(design) paste0(lrow(design,"GreatInflation"), "\n", lrow(design,"COVID"))
covid_ref <- "1: common level (COVID 0.9 ref)"; thr5 <- "2: 5% annualised threshold"
ehdr  <- "Episode & $\\Delta$Overall & $\\Delta$Contemp. & $\\Delta$Lagged"
wr("h3", paste0(
"\\begin{table}[!ht]\n\\centering\n\\caption{H3 --- COVID broadened contemporaneously; the Great Inflation through lagged propagation}\n\\label{tab:h3}\n",
"\\setlength{\\tabcolsep}{6pt}\\footnotesize\n\\begin{threeparttable}\n",
"\\begin{tabular}{l c c c}\n",
"\\multicolumn{4}{l}{\\textit{Panel A. Marginal $\\Delta$TCI$(\\tau)$ vs.\\ the 1983--2019 core, same $\\tau$}}\\\\\n\\toprule\n",
ehdr, " \\\\\n\\midrule\n", epanel(mm),
"\n\\bottomrule\n\\end{tabular}\n\n\\vspace{4pt}\n",
"\\begin{tabular}{l c c c}\n",
"\\multicolumn{4}{l}{\\textit{Panel B. Level-matched to COVID $\\tau=0.9$ (8.7\\% annualised)}}\\\\\n\\toprule\n",
ehdr, " \\\\\n\\midrule\n", lpanel(covid_ref),
"\n\\bottomrule\n\\end{tabular}\n\n\\vspace{4pt}\n",
"\\begin{tabular}{l c c c}\n",
"\\multicolumn{4}{l}{\\textit{Panel C. Level-matched to a common 5\\% annualised threshold}}\\\\\n\\toprule\n",
ehdr, " \\\\\n\\midrule\n", lpanel(thr5),
"\n\\bottomrule\n\\end{tabular}\n",
"\\begin{tablenotes}[flushleft]\\footnotesize\n",
"\\item \\textit{Notes:} Each episode's marginal contribution to monthly connectedness relative to the calm ",
"1983--2019 core (nested subtraction), split into contemporaneous and lagged; six-CPI system, $n_{\\text{lag}}=2$. ",
"Panel~A compares the episodes at the same quantile $\\tau\\in\\{0.7,0.9\\}$. Because a given $\\tau$ corresponds to a ",
"higher inflation level in the Great Inflation, Panels~B and~C instead match on the inflation \\emph{level}: each ",
"episode is evaluated at the $\\tau$ whose annualised month-on-month CPI inflation equals a common reference --- ",
"the COVID $\\tau=0.9$ level of $8.7$ per cent in Panel~B, and $5$ per cent in Panel~C. Throughout, the ",
"Great Inflation broadened mainly through the \\emph{lagged} (propagation) channel, whereas COVID broadened ",
"contemporaneously (its lagged marginal is negative). The level-matched panels show the Great Inflation's larger ",
"broadening is not merely a level artefact. $^{\\dagger}$~The implied $\\tau$ falls outside $[0.05,0.95]$ (the ",
"reference level sits beyond that episode's inflation distribution) and is clamped to the boundary.\n",
"\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n"))

## ============================================================
## H4 -- Energy second-round effects
## ============================================================
## (H4 energy / oil->core removed: oil dropped from the paper.)

## ============================================================
## H5 -- Goods-to-services propagation (recompute on ex-GI sample)
## ============================================================
cat("building h5 (recompute)\n")
codes <- c(Food="CPIUFDSL", Gasoline="CUSR0000SETB01", HHEnergy="CUSR0000SEHF",
           CoreGoods="CUSR0000SACL1E", Shelter="CUSR0000SAH1", CoreServices="CUSR0000SASLE")
sl <- list(); for (nm in names(codes)) { s <- getSymbols(codes[[nm]], src="FRED", auto.assign=FALSE); colnames(s)<-nm; sl[[nm]]<-s }
idx <- do.call(merge, sl)[, names(codes)]
pcz <- zoo(coredata(100*(idx/lag.xts(idx,k=1)-1)), as.yearmon(index(idx))); colnames(pcz)<-names(codes)
oil_raw <- read.csv("WTISPLC.csv"); oil_raw$observation_date <- as.Date(oil_raw$observation_date)
oilz <- zoo(oil_raw$WTISPLC, as.yearmon(oil_raw$observation_date)); oil_pc <- 100*(oilz/stats::lag(oilz,-1)-1)
exS <- (pcz[,"CoreServices"]-W*pcz[,"Shelter"])/(1-W)
cpi6 <- cbind(Food=pcz[,"Food"],Gasoline=pcz[,"Gasoline"],HHEnergy=pcz[,"HHEnergy"],
              CoreGoods=pcz[,"CoreGoods"],Shelter=pcz[,"Shelter"],CoreServ_xS=exS)
dat <- merge(zoo(coredata(cpi6),index(cpi6)), Oil=oil_pc, all=FALSE)
dat <- dat[, c("Oil","Food","Gasoline","HHEnergy","CoreGoods","Shelter","CoreServ_xS")]; dat <- na.omit(dat)
datGI <- window(dat, start=as.yearmon("Jan 1983")); YGI <- as.matrix(datGI)
goods <- "CoreGoods"; svcs <- c("Shelter","CoreServ_xS")
seq_g <- function(nl){ r<-R2ConnectednessQ2(datGI,window.size=NULL,nlag=nl,tau=0.9,shrink=TRUE,drop_own_lags=FALSE,progbar=FALSE)
  L<-r$CT[,,1,2]*100; diag(L)<-0; net<-sum(L[svcs,goods])-sum(L[goods,svcs])
  O<-(r$CT[,,1,1]+r$CT[,,1,2])*100; diag(O)<-0; c(net=net, cg=(colSums(O)-rowSums(O))["CoreGoods"]) }
seq_l <- function(nl){ ct<-R2Q_lasso_CT(YGI,nlag=nl,tau=0.9); L<-ct$L; diag(L)<-0
  net<-sum(L[svcs,goods])-sum(L[goods,svcs]); O<-ct$C+ct$L; diag(O)<-0
  c(net=net, cg=(colSums(O)-rowSums(O))["CoreGoods"]) }
sg<-seq_g(2); l2<-seq_l(2); l6<-seq_l(6)
bodyS <- paste(c(
  sprintf("Genizi dependence ($n_{\\text{lag}}{=}2$) & %s & %s \\\\", sprintf("%+.1f",sg["net"]), sprintf("%+.1f",sg["cg.CoreGoods"])),
  sprintf("LASSO predictive ($n_{\\text{lag}}{=}2$) & %s & %s \\\\", sprintf("%+.1f",l2["net"]), sprintf("%+.1f",l2["cg.CoreGoods"])),
  sprintf("LASSO predictive ($n_{\\text{lag}}{=}6$) & %s & %s \\\\", sprintf("%+.1f",l6["net"]), sprintf("%+.1f",l6["cg.CoreGoods"]))), collapse="\n")
wr("h5", paste0(
"\\begin{table}[!ht]\n\\centering\n\\caption{Goods-to-services sequencing under alternative metrics and lag lengths}\n\\label{tab:h5}\n",
"\\setlength{\\tabcolsep}{6pt}\\footnotesize\n\\begin{threeparttable}\n",
"\\begin{tabular}{l c c}\n\\toprule\n",
"Metric & Goods$\\to$services (lagged, net) & Core goods NET \\\\\n\\midrule\n",
bodyS, "\n\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\\footnotesize\n",
"\\item \\textit{Notes:} ``Goods$\\to$services (lagged, net)'' is lagged connectedness from core goods to core ",
"services (shelter $+$ core services ex-shelter) minus the reverse, at $\\tau=0.9$ on the ex-Great-Inflation sample ",
"(1983--2026); $>0$ means goods lead. The two-lag Genizi dependence metric cannot see a multi-month sequence and ",
"tilts negative; the LASSO predictive metric confirms the sequencing, increasingly at longer lags. ",
"Predictive ordering, not identified causation.\n",
"\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n"))

## ============================================================
## H6 -- Expectations: contemporaneous & lagged, both directions
## ============================================================
cat("building h6 (directional decomposition)\n")
cpi_q  <- to_qpc(getSymbols("CPIAUCSL", src="FRED", auto.assign=FALSE)); colnames(cpi_q) <- "CPI"
mich_q <- read_mich_qtr("MICH_QTR.csv")
de <- merge(MICH=mich_q, CPI=cpi_q, all=FALSE); de <- na.omit(de[, c("MICH","CPI")]); storage.mode(de) <- "double"
sampE <- list(Full=de,
  Core=window(de,start=as.yearqtr("1983 Q1"),end=as.yearqtr("2019 Q4")),
  Y1960_2019=window(de,end=as.yearqtr("2019 Q4")),
  Y1970_2019=window(de,start=as.yearqtr("1970 Q1"),end=as.yearqtr("2019 Q4")),
  Y1983_2026=window(de,start=as.yearqtr("1983 Q1")),
  Y1993_2026=window(de,start=as.yearqtr("1993 Q1")))
labE <- c(Full="Full (1960--2026)", Core="Core (1983--2019)",
          Y1960_2019="1960--2019", Y1970_2019="1970--2019", Y1983_2026="1983--2026",
          Y1993_2026="1993--2026")
dec <- function(d,tau){ r<-R2ConnectednessQ2(d,window.size=NULL,nlag=2,tau=tau,shrink=TRUE,drop_own_lags=FALSE,progbar=FALSE)
  C<-r$CT[,,1,1]*100; L<-r$CT[,,1,2]*100  # [receiver, source]
  O<-C+L; diag(O)<-0
  c(cCE=C["MICH","CPI"], lCE=L["MICH","CPI"], cEC=C["CPI","MICH"], lEC=L["CPI","MICH"],  # CPI->Exp ; Exp->CPI
    NET=unname(colSums(O)["MICH"]-rowSums(O)["MICH"])) }                                  # NET of expectations
# tau within each sample whose quarter-on-quarter CPI inflation equals 5% annualised
lev5_qoq <- 100*((1.05)^(1/4) - 1)                       # QoQ % consistent with 5% ann.
tau5 <- function(d) mean(as.numeric(d[,"CPI"]) <= lev5_qoq, na.rm=TRUE)   # empirical CDF
rowE <- function(s, tau, lead, taulab){ m <- dec(sampE[[s]], tau)
  sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", lead, taulab, f1(m["cCE"]), f1(m["lCE"]),
    f1(m["cEC"]), f1(m["lEC"]), sprintf("%+.1f", m["NET"])) }
blocks <- sapply(names(sampE), function(s){
  ts_raw <- tau5(sampE[[s]]); ts <- min(max(ts_raw, 0.05), 0.95)
  dag <- if (abs(ts_raw - ts) > 1e-9) "$^{\\dagger}$" else ""
  paste(rowE(s, 0.5, labE[s], "0.5"),
        rowE(s, 0.7, "", "0.7"),
        rowE(s, 0.9, "", "0.9"),
        rowE(s, ts, "", paste0(sprintf("%.2f", ts), "$^{\\ddagger}$", dag)), sep="\n") })
bodyH6 <- paste(blocks, collapse="\n\\addlinespace\n")
wr("h6", paste0(
"\\begin{table}[!ht]\n\\centering\n\\caption{H4 --- Contemporaneous and lagged connectedness between inflation and expectations}\n\\label{tab:h6}\n",
"\\setlength{\\tabcolsep}{6pt}\\footnotesize\n\\begin{threeparttable}\n",
"\\begin{tabular}{l c c c c c c}\n\\toprule\n",
" & & \\multicolumn{2}{c}{Inflation $\\to$ Expectations} & \\multicolumn{2}{c}{Expectations $\\to$ Inflation} & \\\\\n",
"\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n",
"Sample & $\\tau$ & Contemp. & Lagged & Contemp. & Lagged & NET \\\\\n\\midrule\n",
bodyH6, "\n\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\\footnotesize\n",
"\\item \\textit{Notes:} Pairwise pseudo-quantile $R^2$ connectedness (\\%) in a two-variable quarterly system of ",
"Michigan inflation expectations (\\texttt{MICH\\_QTR}, from 1960) and overall CPI (CPIAUCSL), $n_{\\text{lag}}=2$ ",
"(a six-month lag horizon). ``Inflation $\\to$ Expectations'' is the share of \\emph{expectations'} variation ",
"explained by CPI, split into contemporaneous and lagged; ``Expectations $\\to$ Inflation'' is the reverse. NET ",
"is the net directional connectedness of expectations (TO $-$ FROM); $>0$ means expectations lead (transmit), ",
"$<0$ that they follow. Dependence, ",
"not identified causation. $^{\\ddagger}$~The fourth row of each block reports connectedness at the quantile ",
"whose quarter-on-quarter CPI inflation equals $5$ per cent annualised within that sample (a common-level ",
"comparison across samples). $^{\\dagger}$~The $5$ per cent level lies outside $[0.05,0.95]$ for that sample and ",
"the quantile is clamped to the boundary.\n",
"\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n"))

## ============================================================
## H6net -- NET directional connectedness of expectations (lead/follow)
##   RETIRED: the NET column is now folded into table h6 above (user
##   request, 2026-07-17). Kept dormant for reference; not written.
## ============================================================
if (FALSE) {
cat("building h6net\n")
qd <- read.csv("results_2var_exp_q/MICH_vs_CPI_directional.csv")
sB <- c("Full","GreatInflation","Core","COVID")
lB <- c(Full="Full (1960--2026)", GreatInflation="Great Inflation (1967--82)",
        Core="Core (1983--2019)", COVID="COVID (2020--26)")
qn <- function(s,t) qd$MICH_NET[qd$nlag==2 & qd$Sample==s & qd$Tau==t]
bodyN <- paste(sapply(sB, function(s)
  sprintf("%s & %s & %s & %s \\\\", lB[s], sprintf("%+.1f",qn(s,0.5)), sprintf("%+.1f",qn(s,0.7)), sprintf("%+.1f",qn(s,0.9)))), collapse="\n")
wr("h6net", paste0(
"\\begin{table}[!ht]\n\\centering\n\\caption{H4 --- Do inflation expectations lead or follow? Net directional connectedness}\n\\label{tab:h6net}\n",
"\\setlength{\\tabcolsep}{6pt}\\footnotesize\n\\begin{threeparttable}\n",
"\\begin{tabular}{l c c c}\n\\toprule\n",
"Sample & $\\tau=0.5$ & $\\tau=0.7$ & $\\tau=0.9$ \\\\\n\\midrule\n",
bodyN, "\n\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}[flushleft]\\footnotesize\n",
"\\item \\textit{Notes:} Net directional connectedness of quarterly Michigan inflation expectations ",
"(\\texttt{MICH\\_QTR}, from 1960) versus overall CPI in a two-variable system, $n_{\\text{lag}}=2$. ",
"$\\text{NET}=\\text{TO}-\\text{FROM}$; $>0$ means expectations lead (transmit), $<0$ means they follow. ",
"Pooled, expectations lead in the upper tail --- but this reflects the anchored 1983--2019 core (NET $+29.1$ ",
"at $\\tau=0.9$); conditional on a high-inflation episode the sign reverses (Great Inflation $-20.2$, COVID ",
"$-6.7$), consistent with adaptive expectations. The COVID window (25 quarters) is short.\n",
"\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n"))
}

cat("\nAll per-hypothesis tables written to", TD, "\n")
