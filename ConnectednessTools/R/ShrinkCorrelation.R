#' Shrink a correlation matrix using data-driven lambda
#'
#' Internal helper that shrinks a correlation matrix using a data-driven
#' shrinkage intensity estimated from the data via corpcor::estimate.lambda().
#'
#' @param R A p x p numeric matrix (intended to be a correlation matrix).
#' @param x Numeric data matrix (n x p) aligned with R.
#' @param verbose Logical; passed to corpcor::estimate.lambda().
#'
#' @return A p x p numeric matrix with diagonal fixed at 1. Attributes:
#'   "lambda" and "lambda.estimated".
#'
#' @keywords internal
.ShrinkCorrelation <- function(R, x, verbose = FALSE) {
  R <- as.matrix(R)
  x <- as.matrix(x)

  if (!is.numeric(R) || !is.numeric(x)) {
    stop("`R` and `x` must be numeric.", call. = FALSE)
  }
  if (nrow(R) != ncol(R)) {
    stop("`R` must be a square matrix.", call. = FALSE)
  }
  if (ncol(x) != ncol(R)) {
    stop("Number of columns in `x` must match dimensions of `R`.", call. = FALSE)
  }
  if (nrow(x) < 2L) {
    stop("`x` must have at least two observations (rows).", call. = FALSE)
  }
  if (!all(is.finite(R))) {
    stop("`R` must contain only finite values.", call. = FALSE)
  }

  # Defensive correlation structure enforcement (avoid brittle exact checks)
  R <- 0.5 * (R + t(R))
  diag(R) <- 1

  n <- nrow(x)
  w <- rep.int(1 / n, n)

  lambda <- corpcor::estimate.lambda(x, w = w, verbose = verbose)

  # Keep your original shrink form (equivalent to (1-lambda)R + lambda I after diag fix)
  R_shrunk <- (1 - lambda) * R
  diag(R_shrunk) <- 1

  attr(R_shrunk, "lambda") <- lambda
  attr(R_shrunk, "lambda.estimated") <- TRUE

  R_shrunk
}
