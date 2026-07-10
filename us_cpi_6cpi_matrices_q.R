# ============================================================
# us_cpi_6cpi_matrices_q.R -- QUARTERLY analogue of the full-sample
#   contemporaneous & lagged connectedness matrices (TO/FROM/NET) for the
#   six-CPI system (NO oil) at tau=0.9, nlag=2.
#   Quarterly build via to_qpc() (monthly index -> quarterly period-average
#   -> QoQ %), matching us_cpi_6cpi_era_q.R. Full sample only.
#   Writes: paper/tables/conn_contemp_full_q.tex, conn_lagged_full_q.tex
#   Labels: tab:conncfullq, tab:connlfullq.
# ============================================================
pkgs <- c("quantmod","zoo","xts","quantreg","MASS","corpcor","Matrix","ConnectednessApproach")
for (p in pkgs) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
source("use_R2Q.R"); source("quarterly_utils.R"); W <- 0.60; TAU <- 0.9; TD <- "paper/tables"
codes <- c(Food="CPIUFDSL", Gasoline="CUSR0000SETB01", HHEnergy="CUSR0000SEHF",
           CoreGoods="CUSR0000SACL1E", Shelter="CUSR0000SAH1", CoreServices="CUSR0000SASLE")
sl <- list(); for (nm in names(codes)) { s <- getSymbols(codes[[nm]], src="FRED", auto.assign=FALSE); colnames(s)<-nm; sl[[nm]]<-s }
idx <- do.call(merge, sl)[, names(codes)]
pcz <- to_qpc(idx); colnames(pcz) <- names(codes)                 # quarterly QoQ %
exS <- (pcz[,"CoreServices"]-W*pcz[,"Shelter"])/(1-W)
DAT <- cbind(Food=pcz[,"Food"], Gasoline=pcz[,"Gasoline"], HHEnergy=pcz[,"HHEnergy"],
             CoreGoods=pcz[,"CoreGoods"], Shelter=pcz[,"Shelter"], CoreServ_xS=exS)
DAT <- na.omit(DAT); storage.mode(DAT) <- "double"
per <- sprintf("quarterly, %s--%s", format(as.yearqtr(start(DAT))), format(as.yearqtr(end(DAT))))
cat(sprintf("Quarterly six-CPI: %d quarters (%s)\n", nrow(DAT), per))
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
"$\\tau=0.9$, $n_{\\text{lag}}=2$, quarterly data. Cell $(i,j)$ is the share of row $i$'s variation ",
"attributed to column $j$. FROM $=$ row sum (received); TO $=$ column sum (transmitted); ",
"NET $=$ TO$-$FROM ($>0$ net transmitter). Diagonal (own) terms are excluded from FROM/TO",
if(!diag_dash) " (and shown for the lagged block as the own-lag/persistence term)" else "", ".\n",
"\\end{tablenotes}\n\\end{threeparttable}\n\\end{table}\n")
  writeLines(txt, file.path(TD, paste0(name,".tex"))); cat("wrote", name, "\n")
}

r <- R2ConnectednessQ2(DAT, window.size=NULL, nlag=2, tau=TAU, shrink=TRUE, drop_own_lags=FALSE, progbar=FALSE)
C <- r$CT[,,1,1]*100; L <- r$CT[,,1,2]*100
emit(C, "conn_contemp_full_q", paste0("Contemporaneous connectedness among CPI subcomponents (", per, ", $\\tau=0.9$)"), "tab:conncfullq", TRUE)
emit(L, "conn_lagged_full_q",  paste0("Lagged connectedness among CPI subcomponents (", per, ", $\\tau=0.9$)"),          "tab:connlfullq", FALSE)

od<-function(X){diag(X)<-0;X}
cat("\n--- Contemporaneous: TO / FROM / NET ---\n")
print(round(rbind(FROM=rowSums(od(C)), TO=colSums(od(C)), NET=colSums(od(C))-rowSums(od(C))),1))
cat("\n--- Lagged: TO / FROM / NET ---\n")
print(round(rbind(FROM=rowSums(od(L)), TO=colSums(od(L)), NET=colSums(od(L))-rowSums(od(L))),1))
cat(sprintf("\nTCI: contemp %.1f, lagged %.1f, overall %.1f\n",
            mean(colSums(od(C))), mean(colSums(od(L))), mean(colSums(od(C+L)))))
cat("\nDone.\n")
