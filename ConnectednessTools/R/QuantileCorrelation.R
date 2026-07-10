#' Quantile Correlation Matrix
#'
#' Computes the symmetric matrix of pairwise quantile correlations at a given
#' quantile level `tau`, following Choi and Shin (2022).
#'
#' @param x A numeric matrix or data frame (columns represent variables). `zoo`/`xts`
#'   objects are accepted and will be coerced via `zoo::coredata()`.
#' @param tau Numeric scalar in (0, 1); the quantile level at which the quantile
#'   correlation is evaluated.
#' @param method Quantile regression fitting method: one of `"fn"`, `"sfn"`, `"br"`.
#' @return A symmetric numeric matrix of dimension `ncol(x)` by `ncol(x)`:
#' \itemize{
#'   \item Diagonal elements are equal to 1.
#'   \item Off-diagonal elements are pairwise quantile correlations (possibly `NA`).
#' }
#'
#' @details
#' For two variables \eqn{X} and \eqn{Y}, the quantile correlation at level \eqn{\tau}
#' is based on the square root of the product of the quantile regression slopes from
#' \eqn{Y \sim X} and \eqn{X \sim Y} at quantile \eqn{\tau}. The sign is taken from
#' the \eqn{Y \sim X} slope. The result is truncated to lie in \eqn{[-1, 1]}.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200), ncol = 4)
#' colnames(x) <- paste0("Var", 1:4)
#' QuantileCorrelation(x, tau = 0.5)
#'
#' @references
#' Choi, J. E., & Shin, D. W. (2022).
#' Quantile correlation coefficient: a new tail dependence measure.
#' Statistical Papers, 63(4), 1075-1104.
#'
#' @export
QuantileCorrelation <- function(x, tau = 0.5, method = "fn") {
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("'tau' must be a finite numeric scalar in (0, 1).", call. = FALSE)
  }

  # accept zoo/xts without requiring xts as a hard dependency
  if (inherits(x, "zoo") || inherits(x, "xts")) x <- zoo::coredata(x)

  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (ncol(x) < 2L) stop("'x' must have at least 2 columns.", call. = FALSE)

  n  <- ncol(x)
  nm <- colnames(x)
  if (is.null(nm)) nm <- paste0("V", seq_len(n))

  corQ <- matrix(NA_real_, n, n, dimnames = list(nm, nm))
  diag(corQ) <- 1

  rq_fit <- switch(
    method,
    br  = function(y, X) quantreg::rq.fit.br (X, y, tau = tau)$coefficients,
    sfn = function(y, X) quantreg::rq.fit.sfn(X, y, tau = tau)$coefficients,
    fn  = function(y, X) quantreg::rq.fit.fnb(X, y, tau = tau)$coefficients,
    stop("Unknown 'method': use one of 'br', 'sfn', 'fn'.", call. = FALSE)
  )

  inds <- which(upper.tri(corQ), arr.ind = TRUE)

  vals <- vapply(seq_len(nrow(inds)), function(ii) {
    i <- inds[ii, 1]
    j <- inds[ii, 2]
    xi <- x[, i]
    xj <- x[, j]

    ok <- is.finite(xi) & is.finite(xj)
    xi <- xi[ok]
    xj <- xj[ok]

    if (length(xi) < 3L) return(NA_real_)
    if (stats::sd(xi) == 0 || stats::sd(xj) == 0) return(NA_real_)

    X1 <- cbind(1, xi)
    b2.1 <- rq_fit(y = xj, X = X1)[2]

    X2 <- cbind(1, xj)
    b1.2 <- rq_fit(y = xi, X = X2)[2]

    qcor <- if ((b2.1 * b1.2) > 0) sign(b2.1) * sqrt(b2.1 * b1.2) else 0
    qcor <- max(-1, min(1, qcor))
    qcor
  }, numeric(1))

  corQ[upper.tri(corQ)] <- vals
  corQ[lower.tri(corQ)] <- t(corQ)[lower.tri(corQ)]
  corQ
}
