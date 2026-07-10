# ============================================================
# R2Q_lasso_dir.R -- DIRECTIONAL LASSO-QR pseudo-R^2 connectedness (BH22 style)
#   with an exact Shapley decomposition over the variable-blocks, so we can
#   read TO/FROM/NET and the lagged source->receiver pass-through.
#
# Same metric as R2Q_lasso.R (Koenker-Machado pseudo-R^2 via LASSO-penalised
# quantile regression), now with directional attribution. Exact Shapley is
# over the k variable-blocks (NOT columns), so it stays feasible at any nlag.
# lambda is selected once per equation (SIC) and held fixed across coalitions.
# ============================================================

suppressMessages(library(quantreg))
.rho <- function(u, tau) u * (tau - (u < 0))

# select lambda* (SIC) on the full standardised design; return lam, L0, full-fit coef
.select_lambda <- function(y, Xstd, tau, ngrid = 15L) {
  n <- length(y); p <- ncol(Xstd); pc <- 2:p
  b0 <- as.numeric(stats::quantile(y, tau, type = 1)); L0 <- sum(.rho(y - b0, tau))
  if (!is.finite(L0) || L0 <= 0) L0 <- 1
  psi <- tau - ((y - b0) < 0)
  lam_max <- max(abs(crossprod(Xstd[, pc, drop = FALSE], psi)))
  if (!is.finite(lam_max) || lam_max <= 0) lam_max <- 1
  grid <- exp(seq(log(lam_max), log(lam_max * 1e-3), length.out = ngrid))
  best <- list(sic = Inf, lam = grid[1], coef = NULL)
  for (lam in grid) {
    fit <- tryCatch(quantreg::rq.fit.lasso(Xstd, y, tau = tau, lambda = c(0, rep(lam, p - 1))),
                    error = function(e) NULL)
    if (is.null(fit)) next
    u <- y - as.vector(Xstd %*% fit$coefficients); L1 <- sum(.rho(u, tau))
    df <- sum(abs(fit$coefficients[pc]) > 1e-7)
    sic <- log(L1 / n) + (df + 1) * log(n) / (2 * n)
    if (is.finite(sic) && sic < best$sic) best <- list(sic = sic, lam = lam, coef = fit$coefficients)
  }
  if (is.null(best$coef)) best$coef <- quantreg::rq.fit(Xstd, y, tau = tau)$coefficients
  list(lam = best$lam, L0 = L0, coef = best$coef)
}

# pseudo-R^2 of a column subset (cols includes intercept index 1) at fixed lambda
.pr2_subset <- function(y, Xstd, cols, tau, lamstar, L0) {
  if (length(cols) <= 1L) return(0)
  Xs <- Xstd[, cols, drop = FALSE]
  fit <- tryCatch(quantreg::rq.fit.lasso(Xs, y, tau = tau, lambda = c(0, rep(lamstar, length(cols) - 1L))),
                  error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  u <- y - as.vector(Xs %*% fit$coefficients)
  1 - sum(.rho(u, tau)) / L0
}

# exact Shapley over `blocks` (each a vector of column indices into Xstd)
.shapley_exact <- function(y, Xstd, blocks, tau, lamstar, L0) {
  K <- length(blocks)
  pr2 <- numeric(2^K)
  for (s in 0:(2^K - 1L)) {
    sel <- as.logical(intToBits(s))[seq_len(K)]
    cols <- c(1L, unlist(blocks[sel], use.names = FALSE))
    pr2[s + 1L] <- .pr2_subset(y, Xstd, cols, tau, lamstar, L0)
  }
  phi <- numeric(K)
  for (k in seq_len(K)) for (s in 0:(2^K - 1L)) {
    sel <- as.logical(intToBits(s))[seq_len(K)]; if (sel[k]) next
    sel2 <- sel; sel2[k] <- TRUE; s2 <- sum(2^(which(sel2) - 1L))
    m <- sum(sel); w <- factorial(m) * factorial(K - m - 1L) / factorial(K)
    phi[k] <- phi[k] + w * (pr2[s2 + 1L] - pr2[s + 1L])
  }
  phi
}

# Returns contemporaneous (C) and lagged (L) matrices, [receiver, source] x100.
R2Q_lasso_CT <- function(Y, nlag, tau) {
  Y <- as.matrix(Y); k <- ncol(Y); NM <- colnames(Y)
  Z <- stats::embed(Y, nlag + 1L)
  yend <- Z[, 1:k, drop = FALSE]
  lagcols <- Z[, (k + 1L):(k * (nlag + 1L)), drop = FALSE]
  sds <- apply(lagcols, 2, stats::sd); sds[sds == 0 | !is.finite(sds)] <- 1
  Xfull <- cbind(1, sweep(lagcols, 2, sds, "/"))
  blocks_lag <- lapply(seq_len(k), function(j) 1L + (j + k * (0:(nlag - 1L))))  # var j's lags -> Xfull cols

  L <- matrix(0, k, k, dimnames = list(NM, NM)); resid <- matrix(NA_real_, nrow(yend), k)
  for (i in seq_len(k)) {
    y <- yend[, i]; sl <- .select_lambda(y, Xfull, tau)
    L[i, ] <- 100 * .shapley_exact(y, Xfull, blocks_lag, tau, sl$lam, sl$L0)
    resid[, i] <- y - as.vector(Xfull %*% sl$coef)
  }
  rsd <- apply(resid, 2, stats::sd); rsd[rsd == 0 | !is.finite(rsd)] <- 1
  C <- matrix(0, k, k, dimnames = list(NM, NM))
  for (i in seq_len(k)) {
    y <- resid[, i]; others <- setdiff(seq_len(k), i)
    Xc <- cbind(1, sweep(resid[, others, drop = FALSE], 2, rsd[others], "/"))
    blocks_c <- as.list(seq_along(others) + 1L)
    sl <- .select_lambda(y, Xc, tau)
    C[i, others] <- 100 * .shapley_exact(y, Xc, blocks_c, tau, sl$lam, sl$L0)
  }
  list(C = C, L = L)
}
