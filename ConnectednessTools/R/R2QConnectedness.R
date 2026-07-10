#' R² Decomposed Connectedness via Quantile Correlation
#'
#' Extends the R² decomposed connectedness framework by replacing linear
#' correlations with the quantile correlation coefficient of Choi and Shin (2022),
#' enabling dependence measurement at a specified quantile level \code{tau}.
#'
#' Internally, for each (rolling) window, the function:
#' \enumerate{
#'   \item Computes a quantile correlation matrix via \code{\link{QuantileCorrelation}}.
#'   \item Applies shrinkage to improve conditioning (internal helper \code{.ShrinkCorrelation}).
#'   \item Ensures positive semi-definiteness (PSD) via an eigenvalue tolerance check
#'         (internal helper \code{.IsPSD}), with a fallback to \code{Matrix::nearPD()}.
#'   \item Computes R² connectedness components and summarizes them using
#'         \code{ConnectednessApproach::ConnectednessTable()}.
#' }
#'
#' @param x A \code{zoo} object containing multivariate time series data
#'   (columns represent variables).
#' @param window.size Integer or \code{NULL}. Rolling window length for dynamic
#'   connectedness estimation. If \code{NULL}, the full sample is used.
#' @param nlag Integer (default \code{1}). Number of lags for the decomposition.
#' @param tau Numeric scalar in \code{(0,1)} (default \code{0.5}). Quantile level
#'   for the quantile correlation coefficient.
#' @param method Character string for quantile regression solver, passed to
#'   \code{\link{QuantileCorrelation}} (e.g., \code{"fn"}, \code{"sfn"}, \code{"br"}).
#' @param progbar Logical (default \code{TRUE}). Display a text progress bar during
#'   rolling estimation.
#' @param return_CT_only Logical (default \code{FALSE}). If \code{TRUE}, returns only
#'   the \code{CT} array and \code{config}, which is useful for bootstrap procedures.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{CT}: 4-D array of component contributions
#'     (\eqn{k \times k \times T \times (1+nlag)}).
#'   \item \code{TO}, \code{FROM}, \code{NET}: directional connectedness measures.
#'   \item \code{TCI}: total connectedness index over time.
#'   \item \code{NPDC}: normalized pairwise directional connectedness.
#'   \item \code{TABLE}: summary tables (overall, contemporaneous, and/or lagged).
#'   \item \code{config}: list of configuration arguments used.
#' }
#'
#' @references
#' Balli, F., Balli, H. O., Dang, T. H. N., & Gabauer, D. (2023).
#' Contemporaneous and lagged R2 decomposed connectedness approach:
#' New evidence from the energy futures market.
#' Finance Research Letters, 57, 104168.
#'
#' Choi, J. E., & Shin, D. W. (2022).
#' Quantile correlation coefficient: a new tail dependence measure.
#' Statistical Papers, 63(4), 1075-1104.
#'
#' @seealso \code{\link{QuantileCorrelation}},
#'   \code{\link[ConnectednessApproach]{ConnectednessTable}},
#'   \code{\link[Matrix]{nearPD}}
#'
#' @examples
#' \donttest{
#' data("dy2012", package = "ConnectednessApproach")
#' dca <- R2QConnectedness(dy2012, window.size = NULL, nlag = 0, 
#'                         tau = 0.5, method = "fn", progbar = FALSE)
#' dca$TABLE
#' }
#'
#' @export
R2QConnectedness <- function(x,
                             window.size = NULL,
                             nlag = 1,
                             tau = 0.5,
                             method = "fn",
                             progbar = TRUE,
                             return_CT_only = FALSE) {

  if (!inherits(x, "zoo")) {
    stop("`x` must be a 'zoo' object.", call. = FALSE)
  }
  if (!is.numeric(nlag) || length(nlag) != 1L || !is.finite(nlag) || nlag < 0) {
    stop("`nlag` must be a single finite number >= 0.", call. = FALSE)
  }
  nlag <- as.integer(nlag)

  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("`tau` must be a finite numeric scalar in (0, 1).", call. = FALSE)
  }

  DATE <- as.character(zoo::index(x))
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  k <- ncol(x)
  if (k < 2L) stop("`x` must have at least 2 columns (variables).", call. = FALSE)

  NAMES <- colnames(x)
  if (is.null(NAMES)) NAMES <- paste0("V", seq_len(k))
  colnames(x) <- NAMES

  Z <- stats::embed(x, nlag + 1L)

  if (is.null(window.size)) {
    window.size_eff <- nrow(Z)
    t0 <- 1L
  } else {
    if (!is.numeric(window.size) || length(window.size) != 1L || !is.finite(window.size)) {
      stop("`window.size` must be NULL or a single finite number.", call. = FALSE)
    }
    window.size <- as.integer(window.size)
    window.size_eff <- window.size - nlag
    if (window.size_eff <= 0L) stop("`window.size` must be > nlag.", call. = FALSE)
    t0 <- nrow(Z) - window.size_eff + 1L
    if (t0 <= 0L) stop("`window.size` is too large for the available sample.", call. = FALSE)
  }

  date <- utils::tail(DATE, t0)

  CT <- array(
    0,
    dim = c(k, k, t0, nlag + 1L),
    dimnames = list(NAMES, NAMES, date, 0:nlag)
  )

  pb <- NULL
  if (isTRUE(progbar)) {
    pb <- utils::txtProgressBar(max = t0, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  for (j in seq_len(t0)) {

    Xwin <- Z[j:(j + window.size_eff - 1L), , drop = FALSE]

    R <- QuantileCorrelation(Xwin, tau = tau, method = method)
    R <- .ShrinkCorrelation(R, Xwin, verbose = FALSE)

    # numeric symmetry before PSD test
    R <- 0.5 * (R + t(R))

    if (!.IsPSD(R, tol = 1e-8)) {
      R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
      R <- R + 1e-10 * diag(nrow(R))
    }

    for (i in seq_len(k)) {
      ryx <- R[-i, i, drop = FALSE]
      rxx <- R[-i, -i, drop = FALSE]

      eigcovx <- eigen(rxx, symmetric = TRUE)
      eigcovx$values <- round(eigcovx$values, 3)
      eigcovx$values <- pmax(eigcovx$values, 0)

      rootcovx <- eigcovx$vectors %*%
        diag(sqrt(eigcovx$values), nrow = length(eigcovx$values)) %*%
        t(eigcovx$vectors)

      cd <- rootcovx^2 %*% (MASS::ginv(rootcovx) %*% ryx)^2

      CT[i, -i, j, 1] <- cd[1:(k - 1)]

      if (nlag > 0) {
        CT[i, , j, 2] <- apply(
          array(cd[-(1:(k - 1))], c(1, k, nlag)),
          1:2, sum
        )
      }
    }

    if (isTRUE(progbar)) utils::setTxtProgressBar(pb, j)
  }

  config <- list(nlag = nlag, approach = "R2", window.size = window.size_eff)

  if (isTRUE(return_CT_only)) {
    return(list(CT = CT, config = config))
  }

  kl <- 1L
  dimensions <- "TCI"
  if (nlag > 0) {
    kl <- 3L
    dimensions <- c("Overall", "Contemporaneous", "Lagged")
  }

  TCI <- array(0, c(t0, kl), dimnames = list(date, dimensions))
  TO <- FROM <- NET <- array(0, c(t0, k, kl), dimnames = list(date, NAMES, dimensions))
  NPDC <- array(0, c(k, k, t0, kl), dimnames = list(NAMES, NAMES, date, dimensions))

  for (tt in seq_len(t0)) {

    if (nlag > 0) {
      ct <- ConnectednessApproach::ConnectednessTable(CT[, , tt, 1])
      lt <- ConnectednessApproach::ConnectednessTable(CT[, , tt, 2])
      at <- ConnectednessApproach::ConnectednessTable(CT[, , tt, 2] + CT[, , tt, 1])

      TO[tt, , 1] <- at$TO;   FROM[tt, , 1] <- at$FROM; NET[tt, , 1] <- at$NET
      NPDC[, , tt, 1] <- at$NPDC; TCI[tt, 1] <- at$TCI

      TO[tt, , 2] <- ct$TO;   FROM[tt, , 2] <- ct$FROM; NET[tt, , 2] <- ct$NET
      NPDC[, , tt, 2] <- ct$NPDC; TCI[tt, 2] <- ct$TCI

      TO[tt, , 3] <- lt$TO;   FROM[tt, , 3] <- lt$FROM; NET[tt, , 3] <- lt$NET
      NPDC[, , tt, 3] <- lt$NPDC; TCI[tt, 3] <- lt$TCI

    } else {
      ct <- ConnectednessApproach::ConnectednessTable(CT[, , tt, 1])

      TO[tt, , 1] <- ct$TO
      FROM[tt, , 1] <- ct$FROM
      NET[tt, , 1] <- ct$NET
      NPDC[, , tt, 1] <- ct$NPDC
      TCI[tt, 1] <- ct$TCI
    }
  }

  TABLE <- ConnectednessApproach::ConnectednessTable(CT[, , , 1])$TABLE
  if (nlag > 0) {
    lt_tab <- ConnectednessApproach::ConnectednessTable(CT[, , , 2])$TABLE
    at_tab <- ConnectednessApproach::ConnectednessTable(CT[, , , 1] + CT[, , , 2])$TABLE
    TABLE <- list(Overall = at_tab, Contemporaneous = TABLE, Lagged = lt_tab)
  }

  if (nlag == 0) {
  TO   <- TO[, , 1, drop = FALSE]
  FROM <- FROM[, , 1, drop = FALSE]
  NET  <- NET[, , 1, drop = FALSE]
  NPDC <- NPDC[, , , 1, drop = FALSE]

  # If t0 == 1, collapse time dimension to plain matrices for TO/FROM/NET
  if (dim(TO)[1] == 1L) {
    TO   <- TO[1, , , drop = TRUE]
    FROM <- FROM[1, , , drop = TRUE]
    NET  <- NET[1, , , drop = TRUE]
    TO   <- as.matrix(TO)
    FROM <- as.matrix(FROM)
    NET  <- as.matrix(NET)

    # NPDC: keep as k x k matrix if single time point
    if (length(dim(NPDC)) == 4L && dim(NPDC)[3] == 1L) {
      NPDC <- NPDC[, , 1, 1, drop = TRUE]
    }
  }
}


  list(
    CT = CT,
    TO = TO,
    FROM = FROM,
    NET = NET,
    TCI = TCI,
    NPDC = NPDC,
    TABLE = TABLE,
    config = config
  )
}
