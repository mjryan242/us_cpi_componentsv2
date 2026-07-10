# ============================================================
# use_R2Q.R  --  shared estimator switch.
#
# Sources Jawad's package estimator (ConnectednessTools::R2QConnectedness)
# and exposes it under the name the drivers already call
# (R2ConnectednessQ2), so every driver uses the SAME estimator by
# sourcing this one file.
#
# Difference from our old clipping estimator: R2QConnectedness repairs a
# non-PSD quantile-correlation matrix with Matrix::nearPD (corr=TRUE)
# before the Genizi decomposition, instead of per-block eigenvalue
# clipping. The two agree wherever the matrix is already PSD (tau<=0.7
# here) and diverge only in the extreme tail (tau=0.9).
#
# Shrinkage: Schafer-Strimmer (.ShrinkCorrelation), but here OVERRIDDEN with
# a PSD-targeted variant -- lambda is raised just enough to make the shrunk
# quantile-correlation matrix positive semi-definite, so nearPD becomes a
# rarely-firing safety net rather than a large tail repair. This de-saturates
# the extreme-tail directional cells. Where the data-driven lambda already
# yields a PSD matrix (tau <= 0.7 here), lambda is unchanged.
# ============================================================

.CT_PKG <- "ConnectednessTools/R"
source(file.path(.CT_PKG, "QuantileCorrelation.R"))
source(file.path(.CT_PKG, "ShrinkCorrelation.R"))
source(file.path(.CT_PKG, "IsPSD.R"))
source(file.path(.CT_PKG, "R2QConnectedness.R"))

# ---- PSD-targeted shrinkage override --------------------------------------
# Redefines .ShrinkCorrelation in the global env; R2QConnectedness looks the
# name up at call time, so it picks up this version.
#
# For a correlation matrix R (unit diagonal), the shrink S = (1-lambda)R + lambda*I
# has eigenvalues (1-lambda)*mu_i + lambda. The smallest is >= psd_tol iff
#   lambda >= (psd_tol - mu_min) / (1 - mu_min),  mu_min = min eigenvalue of R.
# We use lambda = max(lambda_hat, lambda_psd), clamped to [0, 0.999].
.ShrinkCorrelation <- function(R, x, verbose = FALSE, psd_tol = 1e-6) {
  R <- as.matrix(R); x <- as.matrix(x)
  R <- 0.5 * (R + t(R)); diag(R) <- 1
  n <- nrow(x); w <- rep.int(1 / n, n)
  lambda_hat <- corpcor::estimate.lambda(x, w = w, verbose = verbose)

  mu_min <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  lambda_psd <- if (mu_min >= psd_tol) 0 else (psd_tol - mu_min) / (1 - mu_min)

  lambda <- min(max(lambda_hat, lambda_psd), 0.999)

  R_shrunk <- (1 - lambda) * R
  diag(R_shrunk) <- 1
  attr(R_shrunk, "lambda") <- lambda
  attr(R_shrunk, "lambda.hat") <- lambda_hat
  attr(R_shrunk, "lambda.psd") <- lambda_psd
  attr(R_shrunk, "lambda.estimated") <- TRUE
  R_shrunk
}

# Drop-in replacement: same call signature the drivers use.
# `shrink` is accepted for compatibility but ignored (R2QConnectedness
# always shrinks). `drop_own_lags` is not supported by this estimator.
R2ConnectednessQ2 <- function(x, window.size = NULL, nlag = 1, tau = 0.5,
                              shrink = TRUE, method = "fn",
                              drop_own_lags = FALSE, progbar = TRUE) {
  if (isTRUE(drop_own_lags)) {
    stop("drop_own_lags is not supported by R2QConnectedness (it uses ",
         "nearPD PSD repair, not own-lag dropping). Use the clipping ",
         "estimator R_files for the drop-own-lags robustness instead.")
  }
  R2QConnectedness(x, window.size = window.size, nlag = nlag, tau = tau,
                   method = method, progbar = progbar, return_CT_only = FALSE)
}
