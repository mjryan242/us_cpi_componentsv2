# ============================================================
# quarterly_utils.R -- shared helpers for the quarterly-frequency
# robustness analysis.
#
# PRINCIPLE: aggregate the monthly INDEX / price LEVEL to a quarterly
# period-average, THEN take the quarter-on-quarter % change.  Never
# aggregate monthly % changes.  Survey RATE series (e.g. expectations
# already expressed in %) take the quarterly mean of the level, with
# no differencing.
#
# Aggregation uses xts::apply.quarterly + lag.xts (k=1 = previous
# quarter) to avoid the zoo lag() sign ambiguity.  Output is a zoo
# object indexed by yearqtr, matching what the estimators consume.
# ============================================================
suppressMessages({ library(xts); library(zoo) })

# coerce a monthly zoo/xts (yearmon or Date index) to a Date-indexed xts
.as_xts_daily <- function(z) {
  if (xts::is.xts(z)) return(z)
  idx <- zoo::index(z)
  if (inherits(idx, "yearmon")) idx <- as.Date(idx)
  xts::xts(zoo::coredata(z), order.by = as.Date(idx))
}

# monthly index/level (1+ cols) -> quarterly QoQ % change (zoo, yearqtr)
to_qpc <- function(z) {
  x  <- .as_xts_daily(z)
  qx <- xts::apply.quarterly(x, function(v) colMeans(v, na.rm = TRUE))  # period-average level
  pc <- 100 * (qx / xts::lag.xts(qx, k = 1) - 1)                        # QoQ %; k=1 = previous qtr
  out <- zoo::zoo(zoo::coredata(pc), zoo::as.yearqtr(zoo::index(pc)))
  if (!is.null(colnames(pc))) colnames(out) <- colnames(pc)
  out
}

# monthly RATE series (already a %) -> quarterly mean of the level (zoo, yearqtr)
to_qmean <- function(z) {
  x  <- .as_xts_daily(z)
  qx <- xts::apply.quarterly(x, function(v) colMeans(v, na.rm = TRUE))
  out <- zoo::zoo(zoo::coredata(qx), zoo::as.yearqtr(zoo::index(qx)))
  if (!is.null(colnames(qx))) colnames(out) <- colnames(qx)
  out
}

# read MICH_QTR.csv (cols: Quarter, Year, Mean) -> zoo of expectations (%), yearqtr
read_mich_qtr <- function(path = "MICH_QTR.csv") {
  m <- read.csv(path, stringsAsFactors = FALSE)
  yq <- zoo::as.yearqtr(paste0(m$Year, " Q", m$Quarter))
  o  <- order(yq)
  out <- zoo::zoo(as.numeric(m$Mean)[o], yq[o])
  out  # already inflation expectations in % -- enter as-is, NO transform
}
