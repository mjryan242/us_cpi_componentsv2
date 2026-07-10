# ================================================================================
# us_cpi_regime_ms_quarterly.R
#
# QUARTERLY-CONNECTEDNESS variant of the Markov-switching regime Taylor rule.
#
# Motivation: with monthly data and nlag = 2, connectedness only lets spillovers
# propagate two months. Here we aggregate the series to QUARTERLY frequency and
# compute connectedness on the quarterly data with nlag = 2, so spillovers can
# propagate two quarters (~6 months). We then:
#   1. fit the two-state Markov-switching model to the QUARTERLY connectedness
#      series and classify each quarter as high- or low-connectedness;
#   2. broadcast that label to months -- a month is "High" if it falls in a
#      high-connectedness quarter;
#   3. re-estimate the MONTHLY Taylor rule (unchanged) with the monthly High
#      indicator.
#
# Quarterly aggregation follows quarterly_utils.R: aggregate the price INDEX to a
# quarterly period-average then take the % change; survey RATES (MICH) take the
# quarterly mean. Two inflation measures for the expectations system:
#   * exp12 (matched)  : MICH vs 4-quarter (year-over-year) CPI inflation
#   * expmm (mismatch) : MICH vs 1-quarter (quarter-over-quarter) CPI inflation
#
# Windows: 40 quarters (= 120 months) and 20 quarters (= 60 months).
#
# *** WRITES ONLY NEW FILES *** (no existing file is overwritten):
#   tables -> regime_tex/msq_{six,exp12,expmm}_{120,60}.tex
#   plots  -> plots_regime/ms_q/msq_regimes_{six,exp12,expmm}.pdf
# ================================================================================


# --------------------------------------------------------------------------------
# 0. PACKAGES, HELPERS, KNOBS
# --------------------------------------------------------------------------------
suppressMessages({
  library(quantmod); library(zoo); library(xts); library(sandwich); library(lmtest)
  library(MSwM); library(dplyr); library(lubridate); library(tidyr)
})
source("use_R2Q.R")
source("quarterly_utils.R")     # to_qpc(), to_qmean(), .as_xts_daily() -- READ ONLY

WIN_Q   <- c(40, 20)                 # rolling windows in QUARTERS (= 120 and 60 months)
NLAG    <- 2                         # VAR lags (now quarters -> ~6 months of propagation)
PISTAR  <- 2                         # inflation target pi*
K       <- 2                         # Markov states
W_SHELTER <- 0.60                    # shelter weight in core-services-ex-shelter split
TAUS    <- c("0.5", "0.7", "0.9")
MIN_Q   <- 20                        # minimum quarterly obs to attempt an MS fit
SHADE   <- "grey85"

set.seed(20260707)
dir.create("regime_tex", showWarnings = FALSE)
dir.create("plots_regime/ms_q", showWarnings = FALSE, recursive = TRUE)

fetch <- function(code) {
  for (k in 1:6) {
    x <- tryCatch(getSymbols(code, src = "FRED", auto.assign = FALSE), error = function(e) NULL)
    if (!is.null(x)) return(x)
    Sys.sleep(4)
  }
  stop("FRED download failed for ", code)
}
fred_zoo <- function(code) { x <- fetch(code); zoo(as.numeric(coredata(x)), as.yearmon(index(x))) }
lag_one  <- function(z) zoo(c(NA, coredata(z)[-length(z)]), index(z))

# months (as date_n = year+(month-1)/12) belonging to a quarter whose date_q is given
quarter_months <- function(date_q) round(date_q + (0:2) / 12, 3)


# --------------------------------------------------------------------------------
# 1. CONNECTEDNESS ROLLER (total TCI, quarterly panel)
# --------------------------------------------------------------------------------
roll_tci <- function(Sm, win, tau) {
  n <- nrow(Sm)
  out <- rep(NA_real_, n)
  for (t in win:n) {
    sub <- zoo(Sm[(t - win + 1):t, , drop = FALSE], order.by = seq_len(win))
    r <- R2ConnectednessQ2(sub, window.size = NULL, nlag = NLAG, tau = tau, shrink = TRUE, progbar = FALSE)
    total <- r$CT[, , 1, 1] * 100 + r$CT[, , 1, 2] * 100
    diag(total) <- 0
    out[t] <- mean(colSums(total))
  }
  out
}

# roll all three taus for one quarterly panel/window -> tidy data frame keyed by date_q
build_conn_q_df <- function(Sm, yq, win, label) {
  cat(sprintf("  [%s, win=%dq] rolling quarterly connectedness at tau 0.5/0.7/0.9 ...\n", label, win))
  data.frame(
    date_q     = round(as.numeric(yq), 3),           # year + (quarter-1)/4
    `Conn_0.5` = roll_tci(Sm, win, 0.5),
    `Conn_0.7` = roll_tci(Sm, win, 0.7),
    `Conn_0.9` = roll_tci(Sm, win, 0.9),
    check.names = FALSE
  )
}


# --------------------------------------------------------------------------------
# 2. BUILD THE THREE QUARTERLY PANELS
# --------------------------------------------------------------------------------
cat("Downloading data and building quarterly connectedness panels...\n")

cpi <- fetch("CPIAUCSL")                                     # monthly index level (xts)

# quarterly CPI inflation, two horizons, from the quarterly period-average level
qcpi_x <- apply.quarterly(cpi, function(v) mean(v, na.rm = TRUE))
cpi_qoq <- zoo(as.numeric(100 * (qcpi_x / lag.xts(qcpi_x, k = 1) - 1)), as.yearqtr(index(qcpi_x)))  # 1-quarter %
cpi_yoy <- zoo(as.numeric(100 * (qcpi_x / lag.xts(qcpi_x, k = 4) - 1)), as.yearqtr(index(qcpi_x)))  # 4-quarter %

# Michigan expectations -> quarterly mean (already a rate)
mich <- read.csv("MICH.csv", stringsAsFactors = FALSE)
mich$observation_date <- as.Date(mich$observation_date)
michz  <- zoo(mich$MICH, as.yearmon(mich$observation_date))
mich_q <- to_qmean(michz)

# ---- (A) six-CPI subcomponents: index -> quarterly QoQ % of each ----
six_codes <- c(Food = "CPIUFDSL", Gasoline = "CUSR0000SETB01", HHEnergy = "CUSR0000SEHF",
               CoreGoods = "CUSR0000SACL1E", Shelter = "CUSR0000SAH1", CoreServices = "CUSR0000SASLE")
six_list <- list()
for (nm in names(six_codes)) { s <- fetch(six_codes[[nm]]); colnames(s) <- nm; six_list[[nm]] <- s }
six_idx <- do.call(merge, six_list)[, names(six_codes)]
six_qpc <- to_qpc(six_idx)                                   # zoo yearqtr, 6 cols, QoQ %
cs_ex_shelter_q <- (six_qpc[, "CoreServices"] - W_SHELTER * six_qpc[, "Shelter"]) / (1 - W_SHELTER)
six_panel_q <- na.omit(merge(Food = six_qpc[, "Food"], Gasoline = six_qpc[, "Gasoline"],
                             HHEnergy = six_qpc[, "HHEnergy"], CoreGoods = six_qpc[, "CoreGoods"],
                             Shelter = six_qpc[, "Shelter"], CoreServ_xS = cs_ex_shelter_q))
yq_six <- as.yearqtr(index(six_panel_q)); Sm_six <- coredata(six_panel_q); storage.mode(Sm_six) <- "double"

# ---- (B) exp12: MICH vs 4-quarter (YoY) CPI inflation  [matched horizon] ----
sys_exp12 <- na.omit(merge(Exp = mich_q, CPI = cpi_yoy)); yq_exp12 <- as.yearqtr(index(sys_exp12))
Sm_exp12  <- coredata(sys_exp12); colnames(Sm_exp12) <- c("Exp", "CPI")

# ---- (C) expmm: MICH vs 1-quarter (QoQ) CPI inflation  [mismatched horizon] ----
sys_expmm <- na.omit(merge(Exp = mich_q, CPI = cpi_qoq)); yq_expmm <- as.yearqtr(index(sys_expmm))
Sm_expmm  <- coredata(sys_expmm); colnames(Sm_expmm) <- c("Exp", "CPI")


# --------------------------------------------------------------------------------
# 3. MACRO TAYLOR-RULE FRAME (MONTHLY, unchanged; 12-month CPI inflation)
# --------------------------------------------------------------------------------
ffz <- fred_zoo("FEDFUNDS")
raw <- read.csv("WuXiaShadowRate.csv", skip = 1, header = FALSE, stringsAsFactors = FALSE)[, 1:3]
names(raw) <- c("date", "effr", "shadow"); raw <- raw[!is.na(raw$date) & raw$date != "", ]
parts   <- strsplit(raw$date, "-")
mon_abb <- substr(sapply(parts, `[`, 1), 1, 3)
yr_two  <- as.integer(sapply(parts, `[`, 2))
yr_full <- ifelse(yr_two <= 30, 2000 + yr_two, 1900 + yr_two)
shadow_or_effr <- ifelse(!is.na(as.numeric(raw$shadow)), as.numeric(raw$shadow), as.numeric(raw$effr))
shadow_z <- zoo(shadow_or_effr, as.yearmon(yr_full + (match(mon_abb, month.abb) - 1) / 12))
spliced <- merge(shadow_z, ffz)
irate <- spliced[, "shadow_z"]; irate[is.na(irate)] <- spliced[is.na(irate), "ffz"]; irate <- irate[!is.na(irate)]

piz <- zoo(as.numeric(coredata(100 * (cpi / lag.xts(cpi, k = 12) - 1))), as.yearmon(index(cpi)))  # 12-month inflation
un <- merge(u = fred_zoo("UNRATE"), nrou = fred_zoo("NROU"))
un[, "nrou"] <- na.approx(un[, "nrou"], na.rm = FALSE, rule = 2)
ugapz <- un[, "u"] - un[, "nrou"]

macro <- merge(i = irate, pi = piz, ugap = ugapz)
macro <- cbind(macro, pigap = macro[, "pi"] - PISTAR, i_L1 = lag_one(macro[, "i"]))
macro <- macro[complete.cases(macro[, c("i", "i_L1", "pigap", "ugap")]), ]
hold <- data.frame(
  date_n = round(as.numeric(index(macro)), 3),
  i = as.numeric(macro[, "i"]), i_L1 = as.numeric(macro[, "i_L1"]),
  pigap = as.numeric(macro[, "pigap"]), ugap = as.numeric(macro[, "ugap"])
)


# --------------------------------------------------------------------------------
# 4. MS ON QUARTERLY CONNECTEDNESS -> BROADCAST TO MONTHS -> MONTHLY TAYLOR RULE
# --------------------------------------------------------------------------------
est_ms_q <- function(conn_q_df, taucol, start_year) {
  # (a) quarterly connectedness sample from start onward
  d <- data.frame(date_q = conn_q_df$date_q, Conn = conn_q_df[[taucol]])
  d <- d[!is.na(d$Conn) & d$date_q >= start_year, ]
  d <- d[order(d$date_q), ]
  d$Conn_L1 <- c(NA, head(d$Conn, -1))
  d <- d[!is.na(d$Conn_L1), ]
  if (nrow(d) < MIN_Q) return(list(converged = FALSE))

  # (b) 2-state Markov switch on the quarterly connectedness
  ms <- tryCatch(msmFit(lm(Conn ~ Conn_L1, d), k = K, p = 0, sw = rep(TRUE, 3),
                        control = list(parallel = FALSE)), error = function(e) NULL)
  if (is.null(ms)) return(list(converged = FALSE))
  smo <- ms@Fit@smoProb
  if (nrow(smo) == nrow(d) + 1) smo <- smo[-1, , drop = FALSE]
  if (nrow(smo) != nrow(d))     return(list(converged = FALSE))
  state <- max.col(smo, ties.method = "first")

  # (c) high = higher-mean-connectedness quarter
  high_state <- as.integer(names(which.max(tapply(d$Conn, state, mean))))
  d$High <- as.integer(state == high_state)

  # (d) broadcast the quarterly regime to its three months
  monthly_high <- do.call(rbind, lapply(seq_len(nrow(d)), function(k)
    data.frame(date_n = quarter_months(d$date_q[k]), High = d$High[k])))

  # (e) merge onto the MONTHLY macro frame; interacted monthly Taylor rule
  reg <- merge(hold, monthly_high, by = "date_n")
  reg$R1_inf <- reg$High * reg$pigap
  reg <- reg[complete.cases(reg[, c("i", "i_L1", "pigap", "ugap", "R1_inf")]), ]
  if (length(unique(reg$High)) < 2) return(list(converged = FALSE))
  if (nrow(reg) < 30)               return(list(converged = FALSE))

  m1 <- lm(i ~ i_L1 + pigap + R1_inf + ugap, data = reg)
  ct <- tryCatch(coeftest(m1, vcov. = NeweyWest(m1)), error = function(e) NULL)
  if (is.null(ct)) return(list(converged = FALSE))
  b <- coef(m1); rho <- unname(b["i_L1"])
  list(converged = TRUE, ct = ct, n = nrow(reg), high_share = mean(reg$High),
       lr_low = unname(b["pigap"]) / (1 - rho),
       lr_high = unname(b["pigap"] + b["R1_inf"]) / (1 - rho))
}


# --------------------------------------------------------------------------------
# 5. LATEX TABLE (identical format to the monthly script)
# --------------------------------------------------------------------------------
stars <- function(p) { if (is.na(p)) return(""); if (p < 0.01) return("***"); if (p < 0.05) return("**"); if (p < 0.10) return("*"); "" }
coef_cell <- function(res, cn) {
  if (!isTRUE(res$converged)) return("n.c.")
  ct <- res$ct; if (!(cn %in% rownames(ct))) return("--")
  est <- ct[cn, 1]; tval <- ct[cn, 3]; sig <- stars(ct[cn, 4])
  if (nzchar(sig)) sprintf("$%.2f^{%s}$\\,{\\scriptsize(%.1f)}", est, sig, tval)
  else             sprintf("$%.2f$\\,{\\scriptsize(%.1f)}", est, tval)
}
emit_ms <- function(conn_q_df, sys_label, file, fixed_starts) {
  available   <- conn_q_df$date_q[!is.na(conn_q_df$`Conn_0.5`)]
  earliest_yr <- floor(min(available))
  starts <- c(list(list(label = sprintf("Earliest (%d)", earliest_yr), year = -Inf)),
              lapply(fixed_starts, function(y) list(label = as.character(y), year = y)))
  cat(sprintf("  estimating quarterly-connectedness MS models for %s...\n", sys_label))
  results <- list()
  for (tau in TAUS) results[[tau]] <- lapply(starts, function(s) est_ms_q(conn_q_df, paste0("Conn_", tau), s$year))

  ncol_start <- length(starts)
  header <- paste0(" & ", paste(sapply(starts, function(s) s$label), collapse = " & "), " \\\\")
  coef_rows <- list(c("Intercept", "(Intercept)"), c("$\\rho$ ($i_{t-1}$)", "i_L1"),
                    c("$\\beta_\\pi$ ($\\pi-\\pi^*$)", "pigap"),
                    c("$\\Delta\\beta_\\pi^{\\text{High}}$", "R1_inf"), c("$\\beta_x$ (ugap)", "ugap"))
  block_for_tau <- function(tau) {
    rl <- results[[tau]]
    lines <- sprintf("\\multicolumn{%d}{l}{\\textit{$\\tau=%s$}} \\\\", ncol_start + 1, tau)
    for (cr in coef_rows) lines <- c(lines, sprintf("%s & %s \\\\", cr[1],
                                     paste(sapply(rl, function(r) coef_cell(r, cr[2])), collapse = " & ")))
    n_c <- sapply(rl, function(r) if (isTRUE(r$converged)) sprintf("%d", r$n) else "n.c.")
    s_c <- sapply(rl, function(r) if (isTRUE(r$converged)) sprintf("%.2f", r$high_share) else "n.c.")
    l_c <- sapply(rl, function(r) if (isTRUE(r$converged)) sprintf("$%.2f/%.2f$", r$lr_low, r$lr_high) else "n.c.")
    lines <- c(lines, "\\addlinespace",
               sprintf("$n$ (months) & %s \\\\", paste(n_c, collapse = " & ")),
               sprintf("high-conn.\\ share & %s \\\\", paste(s_c, collapse = " & ")),
               sprintf("long-run $\\beta_\\pi$ (low/high) & %s \\\\", paste(l_c, collapse = " & ")))
    paste(lines, collapse = "\n")
  }
  colspec <- paste0("l", paste(rep("c", ncol_start), collapse = ""))
  body <- paste("\\setlength{\\tabcolsep}{4pt}\\scriptsize",
                sprintf("\\begin{tabular}{%s}", colspec), "\\toprule", header, "\\midrule",
                block_for_tau("0.5"), "\\midrule", block_for_tau("0.7"), "\\midrule", block_for_tau("0.9"),
                "\\bottomrule", "\\end{tabular}", sep = "\n")
  writeLines(body, file.path("regime_tex", file))
  cat(sprintf("  wrote regime_tex/%s\n", file))
}


# --------------------------------------------------------------------------------
# 6. REGIME PLOTS (quarterly connectedness line; high quarters shaded)
# --------------------------------------------------------------------------------
regime_series_q <- function(conn_vec, yq) {
  d <- data.frame(date_q = as.numeric(yq), Conn = conn_vec)
  d <- d[!is.na(d$Conn), ]; d <- d[order(d$date_q), ]
  d$Conn_L1 <- c(NA, head(d$Conn, -1)); d$High <- NA_integer_
  fit_rows <- !is.na(d$Conn_L1); df <- d[fit_rows, ]
  if (nrow(df) >= MIN_Q) {
    ms <- tryCatch(msmFit(lm(Conn ~ Conn_L1, df), k = K, p = 0, sw = rep(TRUE, 3),
                          control = list(parallel = FALSE)), error = function(e) NULL)
    if (!is.null(ms)) {
      smo <- ms@Fit@smoProb
      if (nrow(smo) == nrow(df) + 1) smo <- smo[-1, , drop = FALSE]
      if (nrow(smo) == nrow(df)) {
        state <- max.col(smo, ties.method = "first")
        hs <- as.integer(names(which.max(tapply(df$Conn, state, mean))))
        d$High[fit_rows] <- as.integer(state == hs)
      }
    }
  }
  d[, c("date_q", "Conn", "High")]
}
draw_panel <- function(pd, win_months, tau) {
  pd <- pd[order(pd$date_q), ]; tt <- pd$date_q; y <- pd$Conn
  ylim <- range(y, na.rm = TRUE); ylim[2] <- ylim[2] + 0.06 * diff(ylim); xlim <- range(tt)
  plot(tt, y, type = "n", xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE)
  dt <- 0.125                                            # half a quarter, in years
  hi <- ifelse(is.na(pd$High), 0L, pd$High)
  r <- rle(hi == 1L); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  for (k in which(r$values)) rect(tt[starts[k]] - dt, ylim[1], tt[ends[k]] + dt, ylim[2], col = SHADE, border = NA)
  lines(tt, y, lwd = 1.3, col = "black")
  yr0 <- ceiling(xlim[1] / 10) * 10
  axis(1, at = seq(yr0, floor(xlim[2]), by = 10), cex.axis = 0.8, mgp = c(1.6, 0.4, 0), tcl = -0.3)
  axis(2, cex.axis = 0.8, las = 1, mgp = c(1.9, 0.5, 0), tcl = -0.3); box(lwd = 0.8)
  share <- mean(pd$High, na.rm = TRUE)
  title(main = bquote(.(win_months) * "-month,  " * tau == .(tau)), cex.main = 0.95, font.main = 1, line = 0.5)
  if (is.finite(share)) mtext(sprintf("high %.0f%%", 100 * share), side = 3, line = -1.1, adj = 0.97, cex = 0.62, col = "grey30")
}
make_system_pdf <- function(Sm, yq, sys_label, file) {
  cat(sprintf("  [%s] computing quarterly regimes and drawing...\n", sys_label))
  panels <- list()
  for (win in WIN_Q) for (tau in as.numeric(TAUS)) panels[[paste(win, tau)]] <- regime_series_q(roll_tci(Sm, win, tau), yq)
  pdf(file.path("plots_regime/ms_q", file), width = 9.5, height = 5.6)
  par(mfrow = c(length(WIN_Q), length(TAUS)), mar = c(2.4, 3.0, 1.9, 0.7), oma = c(1.4, 1.6, 2.4, 0.4))
  for (win in WIN_Q) for (tau in as.numeric(TAUS)) draw_panel(panels[[paste(win, tau)]], win * 3, tau)
  mtext(sys_label, side = 3, outer = TRUE, cex = 1.05, font = 2, line = 0.7)
  mtext("Quarterly connectedness (TCI, %)", side = 2, outer = TRUE, cex = 0.8, line = 0.2)
  mtext("Grey = high-connectedness Markov regime (quarters); unshaded = low", side = 1, outer = TRUE, cex = 0.72, line = 0.2)
  dev.off()
  cat(sprintf("  wrote plots_regime/ms_q/%s\n", file))
}


# --------------------------------------------------------------------------------
# 7. RUN: tables (both windows) + one plot, for each system
# --------------------------------------------------------------------------------
systems <- list(
  list(key = "six",   plab = "Six-CPI subcomponents (quarterly connectedness)",
       tlab = "Six-CPI subcomponents", Sm = Sm_six,   yq = yq_six,   starts = c(1983, 1993, 2003)),
  list(key = "exp12", plab = "Inflation-expectations, matched (MICH & YoY quarterly CPI)",
       tlab = "matched (YoY)",  Sm = Sm_exp12, yq = yq_exp12, starts = c(1993, 2003)),
  list(key = "expmm", plab = "Inflation-expectations, mismatched (MICH & QoQ quarterly CPI)",
       tlab = "mismatched (QoQ)", Sm = Sm_expmm, yq = yq_expmm, starts = c(1993, 2003))
)

for (win in WIN_Q) {
  cat(sprintf("\n=== window = %d quarters (= %d months) ===\n", win, win * 3))
  for (s in systems) {
    conn_q_df <- build_conn_q_df(s$Sm, s$yq, win, s$key)
    emit_ms(conn_q_df, sprintf("%s, %dq window", s$tlab, win), sprintf("msq_%s_%d.tex", s$key, win * 3), s$starts)
  }
}

cat("\n=== quarterly regime figures ===\n")
for (s in systems) make_system_pdf(s$Sm, s$yq, s$plab, sprintf("msq_regimes_%s.pdf", s$key))

cat("\nDone. Tables regime_tex/msq_*, figures plots_regime/ms_q/. No existing files overwritten.\n")
