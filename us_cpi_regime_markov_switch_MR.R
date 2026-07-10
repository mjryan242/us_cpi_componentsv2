# ================================================================================
# us_cpi_regime_markov_switch_MR.R
#
# ***  SWITCH TO A MARKOV-SWITCHING MODEL  ***
#
#let the connectedness DATA decide the regime: we fit a 2-state Markov-switching AR(1)
# to the connectedness series itself, and each month is assigned to the latent
# state (low- or high-connectedness) it most probably belongs to.
#
# 
#   does the Taylor-rule inflation response differ between low- and high-
#   connectedness months?
# Use regression that has a regime x inflation-gap
# interaction term:
#
#     i_t = a + rho*i_{t-1} + beta_pi*(pi_t - pi*) + delta*(High_t * (pi_t-pi*))
#             + beta_x*x_t + e_t
#
#   i        = Wu-Xia shadow policy rate
#   pi       = 12-month CPI inflation, pi* = 2
#   x        = unemployment gap
#   High_t   = 1 if month t is in the high-connectedness Markov state, else 0
#   delta    = EXTRA inflation response in the high-connectedness regime   
#
# The regression coefficients are reported with Newey-West (HAC) standard errors,
# exactly via  coeftest(m1, vcov. = NeweyWest(m1)).
#
# The script is FLEXIBLE over:
#   * window length       (WIN_SET at the top; connectedness recomputed each window)
#   * quantile tau         (0.5, 0.7, 0.9)
#   * connectedness system (six-CPI subcomponents; inflation<->expectations)
#   * sample start year    (earliest available, 1983, 1993, 2003)
#
# For the inflation<->expectations system we build BOTH pairings:
#   * matched   : MICH 1-yr expected inflation  <->  12-month CPI inflation
#   * mismatched: MICH 1-yr expected inflation  <->  month-on-month CPI inflation
#
# OUTPUT: one LaTeX table fragment per (system x window) in regime_tex/:
#   ms_six_120.tex   ms_expM_120.tex   ms_expX_120.tex
#   ms_six_60.tex    ms_expM_60.tex    ms_expX_60.tex
# Each table has three tau subsections; the columns are the four start dates.
# ================================================================================


# --------------------------------------------------------------------------------
# 0. PACKAGES, HELPERS, AND THE KNOBS YOU CAN TURN
# --------------------------------------------------------------------------------

# Load packages quietly (the libraries still load; only the chatter is hidden).
suppressMessages({
  library(quantmod)   # getSymbols() to download FRED data
  library(zoo)        # time-indexed series (yearmon index)
  library(xts)        # lag.xts() etc.
  library(sandwich)   # NeweyWest() HAC variance
  library(lmtest)     # coeftest()
  library(MSwM)       # msmFit() Markov-switching estimation
  library(dplyr)      # readable data wrangling
  library(lubridate)  # year(), month()
  library(tidyr)      # drop_na() etc.
})

# The pseudo-quantile R^2 connectedness estimator lives here (R2ConnectednessQ2()).
source("use_R2Q.R")

# ---- The knobs -----------------------------------------------------------------
WIN_SET     <- c(120, 60)          # rolling-window lengths (months) to make tables for
NLAG        <- 2                   # VAR lags inside the connectedness estimator
PISTAR      <- 2                   # assumed inflation target pi* (only shifts the intercept)
K           <- 2                   # number of Markov states (low / high connectedness)
W_SHELTER   <- 0.60                # shelter weight in the core-services-ex-shelter split
TAUS        <- c("0.5", "0.7", "0.9")   # the three quantiles
# (the fixed sample-start columns are set per system near the bottom of the file)

set.seed(20260707)   # msmFit uses random starting values -> fix them for reproducibility

# Where the LaTeX fragments go.
dir.create("regime_tex", showWarnings = FALSE)

# Robust FRED download: FRED sometimes returns a transient error, so retry up to 6x.
fetch <- function(code) {
  for (k in 1:6) {
    x <- tryCatch(getSymbols(code, src = "FRED", auto.assign = FALSE),
                  error = function(e) NULL)
    if (!is.null(x)) return(x)
    Sys.sleep(4)
  }
  stop("FRED download failed for ", code)
}

# Download a FRED series and return it as a plain monthly zoo.
fred_zoo <- function(code) {
  x <- fetch(code)
  zoo(as.numeric(coredata(x)), as.yearmon(index(x)))
}

# Shift a zoo back one month (previous month's value, same dates).
lag_one <- function(z) {
  zoo(c(NA, coredata(z)[-length(z)]), index(z))
}


# --------------------------------------------------------------------------------
# 1. THE CONNECTEDNESS ROLLER (used for every system and every window)
# --------------------------------------------------------------------------------

# roll_tci(): total connectedness index over a trailing window of length `win`,
# at quantile `tau`, for the numeric panel `Sm` (rows = months, cols = variables).
# Returns a numeric vector the same length as the number of rows in Sm; the first
# (win - 1) entries are NA because a full window is not yet available.
roll_tci <- function(Sm, win, tau) {
  n_months <- nrow(Sm)
  out <- rep(NA_real_, n_months)

  for (t in win:n_months) {
    # the trailing window ending at month t
    window_rows <- (t - win + 1):t
    sub <- Sm[window_rows, , drop = FALSE]
    # the estimator wants a zoo; the dummy 1..win index is fine (values matter, not dates)
    sub <- zoo(sub, order.by = seq_len(win))

    r <- R2ConnectednessQ2(sub, window.size = NULL, nlag = NLAG,
                           tau = tau, shrink = TRUE, progbar = FALSE)

    # r$CT[,,1,1] = contemporaneous block, r$CT[,,1,2] = lagged block (R^2 shares)
    contemp <- r$CT[, , 1, 1] * 100
    lagged  <- r$CT[, , 1, 2] * 100
    total   <- contemp + lagged
    diag(total) <- 0                      # drop "own" connectedness
    out[t] <- mean(colSums(total))        # Total Connectedness Index for this window
  }

  out
}

# build_conn_df(): run the roller at all three taus for one panel/window and
# return a tidy data frame with a numeric date and the three connectedness columns.
build_conn_df <- function(Sm, ym, win, label) {
  cat(sprintf("  [%s, win=%d] rolling connectedness at tau 0.5 / 0.7 / 0.9 ...\n",
              label, win))

  c50 <- roll_tci(Sm, win, 0.5)
  c70 <- roll_tci(Sm, win, 0.7)
  c90 <- roll_tci(Sm, win, 0.9)

  # date_n = year + (month-1)/12, rounded so it matches the macro frame exactly.
  data.frame(
    date_n     = round(as.numeric(ym), 3),
    `Conn_0.5` = c50,
    `Conn_0.7` = c70,
    `Conn_0.9` = c90,
    check.names = FALSE
  )
}


# --------------------------------------------------------------------------------
# 2. BUILD THE THREE CONNECTEDNESS PANELS (six-CPI, exp-matched, exp-mismatched)
# --------------------------------------------------------------------------------

cat("Downloading data and building connectedness panels...\n")

# ---- shared CPI pieces ----
cpi  <- fetch("CPIAUCSL")                           # headline CPI index level
cpc  <- 100 * (cpi / lag.xts(cpi, k = 1)  - 1)     # month-on-month %
cp12 <- 100 * (cpi / lag.xts(cpi, k = 12) - 1)     # 12-month %

cpc_z  <- zoo(as.numeric(coredata(cpc)),  as.yearmon(index(cpi)))
cp12_z <- zoo(as.numeric(coredata(cp12)), as.yearmon(index(cpi)))

# ---- Michigan 1-year inflation expectations (local CSV) ----
mich <- read.csv("MICH.csv", stringsAsFactors = FALSE)
mich$observation_date <- as.Date(mich$observation_date)
michz <- zoo(mich$MICH, as.yearmon(mich$observation_date))

# ---- (A) six-CPI subcomponent panel ----
# Six seasonally-adjusted component indexes, turned into m/m % changes, with core
# services measured EX-shelter (weight W_SHELTER on shelter). Same panel the main
# six-CPI paper uses.
six_codes <- c(Food = "CPIUFDSL", Gasoline = "CUSR0000SETB01", HHEnergy = "CUSR0000SEHF",
               CoreGoods = "CUSR0000SACL1E", Shelter = "CUSR0000SAH1",
               CoreServices = "CUSR0000SASLE")

six_list <- list()
for (nm in names(six_codes)) {
  series <- fetch(six_codes[[nm]])
  colnames(series) <- nm
  six_list[[nm]] <- series
}
six_idx <- do.call(merge, six_list)[, names(six_codes)]

six_pc <- 100 * (six_idx / lag.xts(six_idx, k = 1) - 1)     # m/m % of each index
six_pc <- zoo(coredata(six_pc), as.yearmon(index(six_idx)))
colnames(six_pc) <- names(six_codes)

core_serv_ex_shelter <- (six_pc[, "CoreServices"] - W_SHELTER * six_pc[, "Shelter"]) /
                        (1 - W_SHELTER)

six_panel <- cbind(
  Food        = six_pc[, "Food"],
  Gasoline    = six_pc[, "Gasoline"],
  HHEnergy    = six_pc[, "HHEnergy"],
  CoreGoods   = six_pc[, "CoreGoods"],
  Shelter     = six_pc[, "Shelter"],
  CoreServ_xS = core_serv_ex_shelter
)
six_panel <- na.omit(six_panel)

ym_six <- as.yearmon(index(six_panel))
Sm_six <- coredata(six_panel)
storage.mode(Sm_six) <- "double"

# ---- (B) inflation<->expectations, MATCHED (MICH <-> 12-month CPI inflation) ----
sys_expM <- na.omit(merge(Exp = michz, CPI = cp12_z))
ym_expM  <- as.yearmon(index(sys_expM))
Sm_expM  <- coredata(sys_expM)

# ---- (C) inflation<->expectations, MISMATCHED (MICH <-> month-on-month CPI) ----
sys_expX <- na.omit(merge(Exp = michz, CPI = cpc_z))
ym_expX  <- as.yearmon(index(sys_expX))
Sm_expX  <- coredata(sys_expX)


# --------------------------------------------------------------------------------
# 3. MACRO DATA FOR THE TAYLOR RULE (policy rate, inflation gap, slack)
# --------------------------------------------------------------------------------

# Effective federal funds rate (used to fill the shadow rate off the lower bound).
ffz <- fred_zoo("FEDFUNDS")

# Wu-Xia shadow short rate from the local CSV.
raw <- read.csv("WuXiaShadowRate.csv", skip = 1, header = FALSE, stringsAsFactors = FALSE)[, 1:3]
names(raw) <- c("date", "effr", "shadow")
raw <- raw[!is.na(raw$date) & raw$date != "", ]

# Parse "Jan-60" style dates, locale-free (two-digit year <=30 -> 2000s).
parts   <- strsplit(raw$date, "-")
mon_abb <- substr(sapply(parts, `[`, 1), 1, 3)
yr_two  <- as.integer(sapply(parts, `[`, 2))
yr_full <- ifelse(yr_two <= 30, 2000 + yr_two, 1900 + yr_two)

# shadow rate where available, else the effective funds rate; then fill remaining gaps.
shadow_or_effr <- ifelse(!is.na(as.numeric(raw$shadow)),
                         as.numeric(raw$shadow), as.numeric(raw$effr))
shadow_z <- zoo(shadow_or_effr, as.yearmon(yr_full + (match(mon_abb, month.abb) - 1) / 12))

spliced <- merge(shadow_z, ffz)
irate <- spliced[, "shadow_z"]
irate[is.na(irate)] <- spliced[is.na(irate), "ffz"]
irate <- irate[!is.na(irate)]

# 12-month CPI inflation is the "inflation" in the rule.
piz <- cp12_z

# Unemployment gap = UNRATE - CBO natural rate (NROU, quarterly -> interpolated).
un <- merge(u = fred_zoo("UNRATE"), nrou = fred_zoo("NROU"))
un[, "nrou"] <- na.approx(un[, "nrou"], na.rm = FALSE, rule = 2)
ugapz <- un[, "u"] - un[, "nrou"]

# Assemble the macro Taylor-rule frame `hold` (one row per month).
macro <- merge(i = irate, pi = piz, ugap = ugapz)
macro <- cbind(macro, pigap = macro[, "pi"] - PISTAR, i_L1 = lag_one(macro[, "i"]))
macro <- macro[complete.cases(macro[, c("i", "i_L1", "pigap", "ugap")]), ]

hold <- data.frame(
  date_n = round(as.numeric(index(macro)), 3),
  i      = as.numeric(macro[, "i"]),
  i_L1   = as.numeric(macro[, "i_L1"]),
  pigap  = as.numeric(macro[, "pigap"]),
  ugap   = as.numeric(macro[, "ugap"])
)


# --------------------------------------------------------------------------------
# 4. ONE MARKOV-SWITCHING REGIME + INTERACTED TAYLOR RULE
# --------------------------------------------------------------------------------

# est_ms(): for one connectedness series (taucol) from one start year, fit the
# 2-state Markov-switching AR(1), label the high-connectedness state, and run the
# interacted Taylor rule. Returns a small list of everything the table needs.
# If the Markov fit fails or the regime is not identified, returns converged=FALSE.
est_ms <- function(conn_df, taucol, start_year) {

  # ---- (a) the connectedness sample from `start_year` onward ----
  d <- data.frame(date_n = conn_df$date_n, Conn = conn_df[[taucol]])
  d <- d[!is.na(d$Conn) & d$date_n >= start_year, ]
  d <- d[order(d$date_n), ]
  d$Conn_L1 <- c(NA, head(d$Conn, -1))       # AR(1) regressor = previous month's connectedness
  d <- d[!is.na(d$Conn_L1), ]

  if (nrow(d) < 40) return(list(converged = FALSE))   # too short to fit a switching model

  # ---- (b) AR(1), then a 2-state Markov switch on it ----
  ar1 <- lm(Conn ~ Conn_L1, data = d)
  ms <- tryCatch(
    msmFit(ar1, k = K, p = 0, sw = rep(TRUE, 3), control = list(parallel = FALSE)),
    error = function(e) NULL
  )
  if (is.null(ms)) return(list(converged = FALSE))

  # ---- (c) smoothed probabilities -> a hard regime for each month ----
  smo <- ms@Fit@smoProb
  if (nrow(smo) == nrow(d) + 1) smo <- smo[-1, , drop = FALSE]   # some versions add a leading row
  if (nrow(smo) != nrow(d))     return(list(converged = FALSE))
  state <- max.col(smo, ties.method = "first")

  # ---- (d) label the state with the higher MEAN connectedness as "high" ----
  mean_conn_by_state <- tapply(d$Conn, state, mean)
  high_state <- as.integer(names(which.max(mean_conn_by_state)))
  d$High <- as.integer(state == high_state)

  # ---- (e) merge the regime into the macro frame; build the interaction ----
  reg <- merge(hold, d[, c("date_n", "High")], by = "date_n")
  reg$R1_inf <- reg$High * reg$pigap                    # regime x inflation gap
  reg <- reg[complete.cases(reg[, c("i", "i_L1", "pigap", "ugap", "R1_inf")]), ]

  # need both regimes present to identify the interaction
  if (length(unique(reg$High)) < 2) return(list(converged = FALSE))
  if (nrow(reg) < 30)               return(list(converged = FALSE))

  # ---- (f) the interacted Taylor rule with Newey-West HAC t-stats ----
  m1 <- lm(i ~ i_L1 + pigap + R1_inf + ugap, data = reg)
  ct <- tryCatch(coeftest(m1, vcov. = NeweyWest(m1)), error = function(e) NULL)
  if (is.null(ct)) return(list(converged = FALSE))

  # long-run inflation response = beta_pi / (1 - rho), low vs high regime
  b   <- coef(m1)
  rho <- unname(b["i_L1"])
  lr_low  <- unname(b["pigap"]) / (1 - rho)
  lr_high <- unname(b["pigap"] + b["R1_inf"]) / (1 - rho)

  list(
    converged  = TRUE,
    ct         = ct,
    n          = nrow(reg),
    high_share = mean(reg$High),
    lr_low     = lr_low,
    lr_high    = lr_high
  )
}


# --------------------------------------------------------------------------------
# 5. FORMAT ONE LATEX TABLE FOR ONE SYSTEM / WINDOW
# --------------------------------------------------------------------------------

# significance stars from a (HAC) p-value: *** 1%, ** 5%, * 10%.
stars <- function(p) {
  if (is.na(p))    return("")
  if (p < 0.01)    return("***")
  if (p < 0.05)    return("**")
  if (p < 0.10)    return("*")
  ""
}

# a coefficient cell "estimate*** (t)"; "n.c." if the model did not converge.
# Stars use the Newey-West HAC p-value (column 4 of the coeftest matrix).
coef_cell <- function(res, coef_name) {
  if (!isTRUE(res$converged)) return("n.c.")
  ct <- res$ct
  if (!(coef_name %in% rownames(ct))) return("--")
  est  <- ct[coef_name, 1]
  tval <- ct[coef_name, 3]
  sig  <- stars(ct[coef_name, 4])
  if (nzchar(sig)) {
    sprintf("$%.2f^{%s}$\\,{\\scriptsize(%.1f)}", est, sig, tval)
  } else {
    sprintf("$%.2f$\\,{\\scriptsize(%.1f)}", est, tval)
  }
}

# emit_ms(): estimate all (tau x start) models for one system/window and write
# the LaTeX fragment. Columns are the start dates; each tau is a subsection.
# `fixed_starts` = the fixed start years to add as columns beyond "earliest"
# (the exp systems drop 1983, which is redundant since their data start later).
emit_ms <- function(conn_df, win, sys_label, file, fixed_starts) {

  # first month actually available for this system/window (the "earliest" column)
  available   <- conn_df$date_n[!is.na(conn_df$`Conn_0.5`)]
  earliest    <- min(available)
  earliest_yr <- floor(earliest)

  # the samples: earliest (no lower bound), then the fixed years for this system
  starts <- c(list(list(label = sprintf("Earliest (%d)", earliest_yr), year = -Inf)),
              lapply(fixed_starts, function(y) list(label = as.character(y), year = y)))

  # estimate every (tau, start) once and store in results[[tau]][[start_index]]
  cat(sprintf("  estimating MS models for %s (win=%d)...\n", sys_label, win))
  results <- list()
  for (tau in TAUS) {
    results[[tau]] <- lapply(starts, function(s)
      est_ms(conn_df, paste0("Conn_", tau), s$year))
  }

  # ---- build the table body ----
  ncol_start <- length(starts)
  header <- paste0(" & ",
                   paste(sapply(starts, function(s) s$label), collapse = " & "),
                   " \\\\")

  # the coefficient rows for one tau block: (row label, coefficient name)
  coef_rows <- list(
    c("Intercept",                          "(Intercept)"),
    c("$\\rho$ ($i_{t-1}$)",                 "i_L1"),
    c("$\\beta_\\pi$ ($\\pi-\\pi^*$)",       "pigap"),
    c("$\\Delta\\beta_\\pi^{\\text{High}}$", "R1_inf"),
    c("$\\beta_x$ (ugap)",                   "ugap")
  )

  block_for_tau <- function(tau) {
    res_list <- results[[tau]]
    lines <- character(0)

    # subsection heading = this tau
    lines <- c(lines, sprintf("\\multicolumn{%d}{l}{\\textit{$\\tau=%s$}} \\\\",
                              ncol_start + 1, tau))

    # coefficient rows
    for (cr in coef_rows) {
      cells <- sapply(res_list, function(r) coef_cell(r, cr[2]))
      lines <- c(lines, sprintf("%s & %s \\\\", cr[1], paste(cells, collapse = " & ")))
    }

    # memo rows: n, high-regime share, long-run beta_pi (low/high)
    n_cells     <- sapply(res_list, function(r) if (isTRUE(r$converged)) sprintf("%d", r$n) else "n.c.")
    share_cells <- sapply(res_list, function(r) if (isTRUE(r$converged)) sprintf("%.2f", r$high_share) else "n.c.")
    lr_cells    <- sapply(res_list, function(r) if (isTRUE(r$converged)) sprintf("$%.2f/%.2f$", r$lr_low, r$lr_high) else "n.c.")

    lines <- c(lines, "\\addlinespace")
    lines <- c(lines, sprintf("$n$ & %s \\\\",                             paste(n_cells,     collapse = " & ")))
    lines <- c(lines, sprintf("high-conn.\\ share & %s \\\\",              paste(share_cells, collapse = " & ")))
    lines <- c(lines, sprintf("long-run $\\beta_\\pi$ (low/high) & %s \\\\", paste(lr_cells,    collapse = " & ")))

    paste(lines, collapse = "\n")
  }

  # column spec: one label column + one per start date
  colspec <- paste0("l", paste(rep("c", ncol_start), collapse = ""))

  body <- paste(
    "\\setlength{\\tabcolsep}{4pt}\\scriptsize",
    sprintf("\\begin{tabular}{%s}", colspec),
    "\\toprule",
    header,
    "\\midrule",
    block_for_tau("0.5"),
    "\\midrule",
    block_for_tau("0.7"),
    "\\midrule",
    block_for_tau("0.9"),
    "\\bottomrule",
    "\\end{tabular}",
    sep = "\n"
  )

  writeLines(body, file.path("regime_tex", file))
  cat(sprintf("  wrote regime_tex/%s\n", file))
}


# --------------------------------------------------------------------------------
# 6. PRODUCE ALL SIX FRAGMENTS (three systems x two windows)
# --------------------------------------------------------------------------------

# `starts` = the fixed start years used as columns (besides "earliest").
# Six-CPI data reach back to the 1970s, so 1983 is a genuine subsample. The two
# expectations systems only start in the 1980s, so a 1983 column would just
# duplicate the "earliest" one -> we drop it there.
systems <- list(
  list(key = "six",  label = "Six-CPI subcomponents",
       Sm = Sm_six,  ym = ym_six,  starts = c(1983, 1993, 2003)),
  list(key = "expM", label = "Inflation$\\leftrightarrow$expectations (matched)",
       Sm = Sm_expM, ym = ym_expM, starts = c(1993, 2003)),
  list(key = "expX", label = "Inflation$\\leftrightarrow$expectations (m/m)",
       Sm = Sm_expX, ym = ym_expX, starts = c(1993, 2003))
)

for (win in WIN_SET) {
  cat(sprintf("\n=== window = %d months ===\n", win))
  for (s in systems) {
    conn_df <- build_conn_df(s$Sm, s$ym, win, s$key)
    file    <- sprintf("ms_%s_%d.tex", s$key, win)
    emit_ms(conn_df, win, s$label, file, s$starts)
  }
}

cat("\nDone. Six LaTeX fragments written to regime_tex/.\n")
