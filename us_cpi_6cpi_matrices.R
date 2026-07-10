# ============================================================
# us_cpi_6cpi_matrices.R -- 6x6 contemporaneous & lagged connectedness
#   matrices (TO/FROM/NET) for the six-CPI system (NO oil) at tau=0.9,
#   nlag=2, for three samples: full (1967-2026), 1967-1982, 1983-2026.
#   Writes six LaTeX fragments to paper/tables/conn_{contemp,lagged}_{full,gi,core}.tex.
# ============================================================
pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R"); W <- 0.60; TAU <- 0.9; TD <- "paper/tables"
codes <- c(Food="CPIUFDSL", Gasoline="CUSR0000SETB01", HHEnergy="CUSR0000SEHF",
           CoreGoods="CUSR0000SACL1E", Shelter="CUSR0000SAH1", CoreServices="CUSR0000SASLE")
sl <- list(); for (nm in names(codes)) { s <- getSymbols(codes[[nm]], src="FRED", auto.assign=FALSE); colnames(s)<-nm; sl[[nm]]<-s }
idx <- do.call(merge, sl)[, names(codes)]
pcz <- zoo(coredata(100*(idx/lag.xts(idx,k=1)-1)), as.yearmon(index(idx))); colnames(pcz)<-names(codes)
exS <- (pcz[,"CoreServices"]-W*pcz[,"Shelter"])/(1-W)
DAT <- cbind(Food=pcz[,"Food"], Gasoline=pcz[,"Gasoline"], HHEnergy=pcz[,"HHEnergy"],
             CoreGoods=pcz[,"CoreGoods"], Shelter=pcz[,"Shelter"], CoreServ_xS=exS)
DAT <- na.omit(DAT); storage.mode(DAT) <- "double"
lab <- c("Food","Gasoline","HH energy","Core goods","Shelter","Core svc.\\ ex-sh.")
f1 <- function(x) sprintf("%.1f", x)

emit <- function(M, name, caption, label, diag_dash){
  Mo<-M; diag(Mo)<-0; FROM<-rowSums(Mo); TO<-colSums(Mo); NET<-TO-FROM
  body <- paste(sapply(1:6, function(i){
    cells <- sapply(1:6, function(j) if (i==j && diag_dash) "---" else f1(M[i,j]))
    sprintf("%s & %s & %s \\\\", lab[i], paste(cells, collapse=" & "), f1(FROM[i])) }), collapse="\n")
  tor <- paste("TO", paste(f1(TO), collapse=" & "), "", sep=" & ")
  netr<- paste("NET", paste(sprintf("%+.1f",NET), collapse=" & "), "", sep=" & ")
  txt <- paste0(
"\\begin{table}[!ht]\n\\centering\n\\caption{", caption, "}\n\\label{", label, "}\n",
"\\setlength{\\tabcolsep}{4pt}\\footnotesize\n\\begin{threeparttable}\n",
"\\begin{tabular}{l c c c c c c c}\n\\toprule\n",
" & Food & Gas. & HH en. & Core gds & Shelter & Core svc. & FROM \\\\\n\\midrule\n",
body, "\n\\midrule\n", tor, " \\\\\n", netr, " \\\\\n\\bottomrule\n\\end{tabular}\n",
"\\begin{tablenotes}[flushleft]\\footnotesize\n",
"\\item \\textit{Notes:} Pseudo-quantile $R^2$ connectedness (\\%) among the six CPI subcomponents at ",
"$\\tau=0.9$, $n_{\\text{lag}}=2$. Cell $(i,j)$ is the share of row $i$'s variation attributed to ",
"column $j$. FROM $=$ row sum (received); TO $=$ column sum (transmitted); NET $=$ TO$-$FROM ",
"($>0$ net transmitter). Diagonal (own) terms are excluded from FROM/TO", if(!diag_dash) " (and shown for the lagged block as the own-lag/persistence term)" else "", ".\n",
"\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n")
  writeLines(txt, file.path(TD, paste0(name,".tex"))); cat("wrote", name, "\n")
}

samples <- list(
  full = list(d=DAT, per="full sample, 1967--2026", tag="full"),
  gi   = list(d=window(DAT, end=as.yearmon("Dec 1982")),   per="1967--1982", tag="gi"),
  core = list(d=window(DAT, start=as.yearmon("Jan 1983")), per="1983--2026", tag="core"))
for (s in samples){
  r <- R2ConnectednessQ2(s$d, window.size=NULL, nlag=2, tau=TAU, shrink=TRUE, drop_own_lags=FALSE, progbar=FALSE)
  C <- r$CT[,,1,1]*100; L <- r$CT[,,1,2]*100
  emit(C, paste0("conn_contemp_",s$tag), paste0("Contemporaneous connectedness among CPI subcomponents (", s$per, ", $\\tau=0.9$)"), paste0("tab:connc",s$tag), TRUE)
  emit(L, paste0("conn_lagged_", s$tag), paste0("Lagged connectedness among CPI subcomponents (", s$per, ", $\\tau=0.9$)"),         paste0("tab:connl",s$tag), FALSE)
}
cat("\nDone.\n")
