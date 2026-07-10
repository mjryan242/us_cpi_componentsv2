#' Check if a matrix is positive semi-definite (PSD)
#'
#' Internal helper that tests whether a square numeric matrix is PSD within a
#' numerical tolerance. The matrix is symmetrized and then eigenvalues are checked
#' against \code{-tol}.
#'
#' @param M A square numeric matrix.
#' @param tol Numeric tolerance allowing small negative eigenvalues from numerical
#'   rounding (default \code{1e-8}).
#'
#' @return Logical scalar. \code{TRUE} if PSD within tolerance, otherwise \code{FALSE}.
#'
#' @keywords internal
.IsPSD <- function(M, tol = 1e-8) {
  if (!is.matrix(M)) stop("`M` must be a matrix.", call. = FALSE)
  if (!is.numeric(M)) stop("`M` must be numeric.", call. = FALSE)
  if (length(tol) != 1L || !is.numeric(tol) || !is.finite(tol) || tol < 0) {
    stop("`tol` must be a single non-negative finite number.", call. = FALSE)
  }
  if (nrow(M) != ncol(M)) stop("`M` must be square.", call. = FALSE)

  # symmetrize to ensure a real spectrum (up to numerical error)
  M <- (M + t(M)) / 2

  eigvals <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  all(eigvals >= -tol)
}
