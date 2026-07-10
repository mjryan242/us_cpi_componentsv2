# ================================================================================
# us_cpi_msq_exp12_addcol.R
#
# Rebuilds regime_tex/msq_exp12_120.tex (matched YoY expectations, quarterly
# connectedness, 120-month window) with an EXTRA left-hand column for the EXTENDED
# sample. The three original columns use the monthly Michigan survey (from 1978)
# aggregated to quarters (earliest ~1987); the extra column uses the native
# QUARTERLY Michigan survey back to 1960 (earliest ~1969).
#
# Columns:  Extended (1969) | Earliest (1987) | 1993 | 2003
#   - Extended : long series (MICH_QTR), full sample
#   - other 3  : short series (aggregated monthly MICH), as before
#
# Only regime_tex/msq_exp12_120.tex is written.
# ================================================================================
suppressMessages({ library(quantmod); library(zoo); library(xts); library(sandwich)
                   library(lmtest); library(MSwM); library(dplyr); library(lubridate); library(tidyr) })
source("use_R2Q.R"); source("quarterly_utils.R")
NLAG <- 2; PISTAR <- 2; K <- 2; WIN <- 40; MIN_Q <- 20; TAUS <- c("0.5", "0.7", "0.9")
set.seed(20260707)

fetch <- function(code){for(k in 1:6){x<-tryCatch(getSymbols(code,src="FRED",auto.assign=FALSE),error=function(e)NULL);if(!is.null(x))return(x);Sys.sleep(4)};stop("fail ",code)}
fred_zoo <- function(code){x<-fetch(code);zoo(as.numeric(coredata(x)),as.yearmon(index(x)))}
lag_one  <- function(z) zoo(c(NA, coredata(z)[-length(z)]), index(z))
quarter_months <- function(dq) round(dq + (0:2)/12, 3)

roll_tci <- function(Sm, win, tau){
  n <- nrow(Sm); out <- rep(NA_real_, n)
  for (t in win:n){ sub <- zoo(Sm[(t-win+1):t,,drop=FALSE], order.by=seq_len(win))
    r <- R2ConnectednessQ2(sub, window.size=NULL, nlag=NLAG, tau=tau, shrink=TRUE, progbar=FALSE)
    tot <- r$CT[,,1,1]*100 + r$CT[,,1,2]*100; diag(tot) <- 0; out[t] <- mean(colSums(tot)) }
  out
}
build_conn <- function(Sm, yq, win){
  data.frame(date_q=round(as.numeric(yq),3),
             `Conn_0.5`=roll_tci(Sm,win,0.5), `Conn_0.7`=roll_tci(Sm,win,0.7), `Conn_0.9`=roll_tci(Sm,win,0.9),
             check.names=FALSE)
}

## ---- data ----
cat("Building short and long matched (YoY) panels...\n")
cpi <- fetch("CPIAUCSL")
qcpi_x <- apply.quarterly(cpi, function(v) mean(v, na.rm=TRUE))
cpi_yoy <- zoo(as.numeric(100*(qcpi_x/lag.xts(qcpi_x,k=4)-1)), as.yearqtr(index(qcpi_x)))

mich <- read.csv("MICH.csv", stringsAsFactors=FALSE); mich$observation_date <- as.Date(mich$observation_date)
mich_q_short <- to_qmean(zoo(mich$MICH, as.yearmon(mich$observation_date)))   # monthly survey -> quarterly (from 1978)
mich_q_long  <- read_mich_qtr("MICH_QTR.csv")                                 # native quarterly survey (from 1960)

sys_short <- na.omit(merge(Exp=mich_q_short, CPI=cpi_yoy)); Sm_short <- coredata(sys_short); colnames(Sm_short)<-c("Exp","CPI"); yq_short <- as.yearqtr(index(sys_short))
sys_long  <- na.omit(merge(Exp=mich_q_long,  CPI=cpi_yoy)); Sm_long  <- coredata(sys_long);  colnames(Sm_long) <-c("Exp","CPI"); yq_long  <- as.yearqtr(index(sys_long))

## ---- monthly macro frame ----
ffz <- fred_zoo("FEDFUNDS")
raw <- read.csv("WuXiaShadowRate.csv", skip=1, header=FALSE, stringsAsFactors=FALSE)[,1:3]; names(raw)<-c("date","effr","shadow"); raw<-raw[!is.na(raw$date)&raw$date!="",]
pp<-strsplit(raw$date,"-"); m2<-substr(sapply(pp,`[`,1),1,3); y2<-as.integer(sapply(pp,`[`,2)); yr2<-ifelse(y2<=30,2000+y2,1900+y2)
ish<-ifelse(!is.na(as.numeric(raw$shadow)),as.numeric(raw$shadow),as.numeric(raw$effr)); isf<-zoo(ish,as.yearmon(yr2+(match(m2,month.abb)-1)/12)); spx<-merge(isf,ffz)
irate<-spx[,"isf"]; irate[is.na(irate)]<-spx[is.na(irate),"ffz"]; irate<-irate[!is.na(irate)]
piz<-zoo(as.numeric(coredata(100*(cpi/lag.xts(cpi,k=12)-1))),as.yearmon(index(cpi)))
un<-merge(u=fred_zoo("UNRATE"),nrou=fred_zoo("NROU")); un[,"nrou"]<-na.approx(un[,"nrou"],na.rm=FALSE,rule=2); ugapz<-un[,"u"]-un[,"nrou"]
macro<-merge(i=irate,pi=piz,ugap=ugapz); macro<-cbind(macro,pigap=macro[,"pi"]-PISTAR,i_L1=lag_one(macro[,"i"]))
macro<-macro[complete.cases(macro[,c("i","i_L1","pigap","ugap")]),]
hold<-data.frame(date_n=round(as.numeric(index(macro)),3), i=as.numeric(macro[,"i"]), i_L1=as.numeric(macro[,"i_L1"]), pigap=as.numeric(macro[,"pigap"]), ugap=as.numeric(macro[,"ugap"]))

## ---- MS on quarterly connectedness -> broadcast -> monthly Taylor rule ----
est_ms_q <- function(conn, taucol, start){
  d <- data.frame(date_q=conn$date_q, Conn=conn[[taucol]]); d<-d[!is.na(d$Conn)&d$date_q>=start,]; d<-d[order(d$date_q),]
  d$Conn_L1<-c(NA,head(d$Conn,-1)); d<-d[!is.na(d$Conn_L1),]; if(nrow(d)<MIN_Q) return(list(converged=FALSE))
  ms<-tryCatch(msmFit(lm(Conn~Conn_L1,d),k=K,p=0,sw=rep(TRUE,3),control=list(parallel=FALSE)),error=function(e)NULL)
  if(is.null(ms)) return(list(converged=FALSE))
  smo<-ms@Fit@smoProb; if(nrow(smo)==nrow(d)+1) smo<-smo[-1,,drop=FALSE]; if(nrow(smo)!=nrow(d)) return(list(converged=FALSE))
  state<-max.col(smo,ties.method="first"); hs<-as.integer(names(which.max(tapply(d$Conn,state,mean)))); d$High<-as.integer(state==hs)
  mh<-do.call(rbind,lapply(seq_len(nrow(d)),function(k) data.frame(date_n=quarter_months(d$date_q[k]),High=d$High[k])))
  reg<-merge(hold,mh,by="date_n"); reg$R1_inf<-reg$High*reg$pigap; reg<-reg[complete.cases(reg[,c("i","i_L1","pigap","ugap","R1_inf")]),]
  if(length(unique(reg$High))<2) return(list(converged=FALSE)); if(nrow(reg)<30) return(list(converged=FALSE))
  m1<-lm(i~i_L1+pigap+R1_inf+ugap,reg); ct<-tryCatch(coeftest(m1,vcov.=NeweyWest(m1)),error=function(e)NULL); if(is.null(ct)) return(list(converged=FALSE))
  b<-coef(m1); rho<-unname(b["i_L1"])
  list(converged=TRUE, ct=ct, n=nrow(reg), high_share=mean(reg$High), lr_low=unname(b["pigap"])/(1-rho), lr_high=unname(b["pigap"]+b["R1_inf"])/(1-rho))
}
stars <- function(p){ if(is.na(p)) return(""); if(p<0.01) return("***"); if(p<0.05) return("**"); if(p<0.10) return("*"); "" }
cell <- function(res,cn){ if(!isTRUE(res$converged)) return("n.c."); ct<-res$ct; if(!(cn%in%rownames(ct))) return("--")
  s<-stars(ct[cn,4]); if(nzchar(s)) sprintf("$%.2f^{%s}$\\,{\\scriptsize(%.1f)}",ct[cn,1],s,ct[cn,3]) else sprintf("$%.2f$\\,{\\scriptsize(%.1f)}",ct[cn,1],ct[cn,3]) }

## ---- estimate: extra column (long, full) + original three (short) ----
cat("Estimating...\n")
conn_long  <- build_conn(Sm_long,  yq_long,  WIN)
conn_short <- build_conn(Sm_short, yq_short, WIN)
yr_long  <- floor(min(conn_long$date_q [!is.na(conn_long$`Conn_0.5`)]))
yr_short <- floor(min(conn_short$date_q[!is.na(conn_short$`Conn_0.5`)]))

R <- list()
for (tau in TAUS) R[[tau]] <- list(
  ext = est_ms_q(conn_long,  paste0("Conn_",tau), -Inf),   # EXTENDED (long, full)
  e   = est_ms_q(conn_short, paste0("Conn_",tau), -Inf),   # short earliest
  a   = est_ms_q(conn_short, paste0("Conn_",tau), 1993),
  b   = est_ms_q(conn_short, paste0("Conn_",tau), 2003)
)

## ---- emit 4-column table ----
coef_rows <- list(c("Intercept","(Intercept)"), c("$\\rho$ ($i_{t-1}$)","i_L1"),
                  c("$\\beta_\\pi$ ($\\pi-\\pi^*$)","pigap"), c("$\\Delta\\beta_\\pi^{\\text{High}}$","R1_inf"),
                  c("$\\beta_x$ (ugap)","ugap"))
block <- function(tau){
  r <- R[[tau]]; ord <- list(r$ext, r$e, r$a, r$b)
  L <- sprintf("\\multicolumn{5}{l}{\\textit{$\\tau=%s$}} \\\\", tau)
  for (cr in coef_rows) L <- c(L, sprintf("%s & %s \\\\", cr[1], paste(sapply(ord,function(x) cell(x,cr[2])), collapse=" & ")))
  nc <- sapply(ord, function(x) if(isTRUE(x$converged)) sprintf("%d",x$n) else "n.c.")
  sc <- sapply(ord, function(x) if(isTRUE(x$converged)) sprintf("%.2f",x$high_share) else "n.c.")
  lc <- sapply(ord, function(x) if(isTRUE(x$converged)) sprintf("$%.2f/%.2f$",x$lr_low,x$lr_high) else "n.c.")
  L <- c(L, "\\addlinespace",
         sprintf("$n$ (months) & %s \\\\", paste(nc,collapse=" & ")),
         sprintf("high-conn.\\ share & %s \\\\", paste(sc,collapse=" & ")),
         sprintf("long-run $\\beta_\\pi$ (low/high) & %s \\\\", paste(lc,collapse=" & ")))
  paste(L, collapse="\n")
}
header <- sprintf(" & Extended$^{\\dagger}$ (%d) & Earliest (%d) & 1993 & 2003 \\\\", yr_long, yr_short)
body <- paste("\\setlength{\\tabcolsep}{4pt}\\scriptsize", "\\begin{tabular}{lcccc}", "\\toprule", header, "\\midrule",
              block("0.5"), "\\midrule", block("0.7"), "\\midrule", block("0.9"), "\\bottomrule", "\\end{tabular}", sep="\n")
writeLines(body, "regime_tex/msq_exp12_120.tex")
cat("wrote regime_tex/msq_exp12_120.tex  (extended col =", yr_long, "; short cols from", yr_short, ")\n")
