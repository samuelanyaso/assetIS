## ============================================================
## Internal utility functions for assetIS
## None of these are exported; they are called by the four main
## estimators and the Newton solver.
## ============================================================

# ----------------------------------------------------------------
# Subset enumeration
# ----------------------------------------------------------------

#' @keywords internal
.generate_subsets <- function(C) {
  lapply(seq_len(2^C - 1L),
         function(mask) which(as.logical(intToBits(mask))[seq_len(C)]))
}

# ----------------------------------------------------------------
# Log-sum-exp (numerically stable)
# ----------------------------------------------------------------

#' @keywords internal
.logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# ----------------------------------------------------------------
# Weight matrix for independent Gaussian studies
# W[m, j] = sqrt(n_m) / sqrt(sum_{l in A_j} n_l)  for m in A_j
# ----------------------------------------------------------------

#' @keywords internal
.build_weight_matrix_ind <- function(n, subsets) {
  M      <- length(n)
  J      <- length(subsets)
  sqrt_n <- sqrt(n)
  W      <- matrix(0.0, nrow = M, ncol = J)
  denom  <- vapply(subsets, function(A) sqrt(sum(n[A])), numeric(1L))
  for (j in seq_len(J)) {
    A       <- subsets[[j]]
    W[A, j] <- sqrt_n[A] / denom[j]
  }
  W
}

# ----------------------------------------------------------------
# GLS weight matrix for correlated Gaussian / non-Gaussian studies
# W[c, j] = [Sigma_A^{-1} 1_A]_c / sqrt(1_A' Sigma_A^{-1} 1_A)
# ----------------------------------------------------------------

#' @keywords internal
.build_weight_matrix_gls <- function(Sigma, subsets, C, ridge = 1e-10) {
  J <- length(subsets)
  W <- matrix(0.0, nrow = C, ncol = J)
  for (j in seq_len(J)) {
    A   <- subsets[[j]]
    S_A <- Sigma[A, A, drop = FALSE]
    S_A <- (S_A + t(S_A)) / 2 + diag(ridge, length(A))
    one <- rep(1.0, length(A))
    a   <- tryCatch(solve(S_A, one),
                    error = function(e)
                      stop(sprintf("Sigma_A is singular for subset %d: %s",
                                   j, conditionMessage(e))))
    denom  <- sqrt(as.numeric(crossprod(one, a)))
    W[A, j] <- a / denom
  }
  W
}

# ----------------------------------------------------------------
# Robust Cholesky (for Gaussian correlated estimator)
# ----------------------------------------------------------------

#' @keywords internal
.safe_chol <- function(S) {
  U <- tryCatch(chol(S), error = function(e) NULL)
  if (!is.null(U)) return(U)
  if (requireNamespace("Matrix", quietly = TRUE)) {
    S_pd <- as.matrix(Matrix::nearPD(S, corr = FALSE)$mat)
    S_pd <- (S_pd + t(S_pd)) / 2
    U <- tryCatch(chol(S_pd), error = function(e) NULL)
    if (!is.null(U)) {
      warning("'Sigma' was not positive definite; nearPD regularisation applied.")
      return(U)
    }
    S <- S_pd
  }
  for (k in seq_len(8L)) {
    jitter <- 1e-10 * 10^(k - 1L)
    U <- tryCatch(chol(S + diag(jitter, nrow(S))), error = function(e) NULL)
    if (!is.null(U)) {
      warning(sprintf("'Sigma' regularised with diagonal jitter %.2e.", jitter))
      return(U)
    }
  }
  stop("Failed to compute Cholesky of 'Sigma' after regularisation.")
}

# ----------------------------------------------------------------
# Newton--Raphson solver: phi'_A(xi) = b0
#
# Works for any per-subject allele frequency vector f_subj, so it
# handles both common-frequency and multi-ancestry settings.
#
# phi'(xi) = sum_i 2 f_i omega_i exp(xi omega_i) /
#                ((1 - f_i) + f_i exp(xi omega_i))
#           - 2 sum_i f_i omega_i
# Setting phi'(xi) = b is equivalent to solving f(xi) = 0 where
# f(xi) = phi'(xi) - b.
# ----------------------------------------------------------------

#' @keywords internal
.newton_xi <- function(b0, omega, f_subj,
                       tol      = 1e-10,
                       max_iter = 100L) {

  B_const <- b0 + 2 * sum(f_subj * omega)

  fval <- function(xi) {
    x     <- xi * omega
    denom <- (1 - f_subj) + f_subj * exp(x)
    sum(2 * f_subj * omega * exp(x) / denom) - B_const
  }

  fpval <- function(xi) {
    x     <- xi * omega
    denom <- (1 - f_subj) + f_subj * exp(x)
    sum(2 * f_subj * (1 - f_subj) * omega^2 * exp(x) / denom^2)
  }

  xi <- 0; fx <- fval(xi)

  for (iter in seq_len(max_iter)) {
    if (abs(fx) < tol)
      return(list(xi = xi, converged = TRUE, iter = iter - 1L))
    fp <- fpval(xi)
    if (!is.finite(fp) || fp <= 0)
      stop(sprintf("Newton step invalid at iter %d (fp = %.3e).", iter, fp))
    step   <- -fx / fp
    xi_new <- xi + step
    fx_new <- fval(xi_new)
    n_halv <- 0L
    while ((!is.finite(fx_new) || abs(fx_new) >= abs(fx)) && n_halv < 50L) {
      step <- step / 2; xi_new <- xi + step
      fx_new <- fval(xi_new); n_halv <- n_halv + 1L
    }
    xi <- xi_new; fx <- fx_new
    if (abs(step) < tol * (1 + abs(xi)))
      return(list(xi = xi, converged = TRUE, iter = iter))
  }
  list(xi = xi, converged = FALSE, iter = max_iter)
}

# ----------------------------------------------------------------
# CGF phi_A(xi) for the trinomial genotype model
# phi_A(xi) = sum_i log[(1-f_i)^2 exp(-2f_i xi w_i)
#                      + 2f_i(1-f_i) exp((1-2f_i) xi w_i)
#                      + f_i^2 exp((2-2f_i) xi w_i)]
# ----------------------------------------------------------------

#' @keywords internal
.phi_nonnorm <- function(xi, omega, f_subj) {
  xw <- xi * omega
  t0 <- (1 - f_subj)^2 * exp(-2 * f_subj * xw)
  t1 <- 2 * f_subj * (1 - f_subj) * exp((1 - 2 * f_subj) * xw)
  t2 <- f_subj^2 * exp((2 - 2 * f_subj) * xw)
  sum(log(t0 + t1 + t2))
}

# ----------------------------------------------------------------
# Precompute xi_vec and phi_vec for all subsets (non-Gaussian IS)
# ----------------------------------------------------------------

#' @keywords internal
.precompute_xi_phi <- function(omega_list, f_subj, b0,
                               tol, max_iter) {
  J       <- length(omega_list)
  xi_vec  <- numeric(J)
  phi_vec <- numeric(J)
  for (j in seq_len(J)) {
    res <- .newton_xi(b0, omega_list[[j]], f_subj, tol, max_iter)
    if (!res$converged)
      warning(sprintf("Newton's method did not converge for subset %d.", j))
    xi_vec[j]  <- res$xi
    phi_vec[j] <- .phi_nonnorm(res$xi, omega_list[[j]], f_subj)
  }
  list(xi_vec = xi_vec, phi_vec = phi_vec)
}

# ----------------------------------------------------------------
# Split K simulations into per-core chunks
# ----------------------------------------------------------------

#' @keywords internal
.make_chunks <- function(K, n_cores) {
  base_c <- K %/% n_cores
  rem    <- K %%  n_cores
  cs     <- rep(base_c, n_cores)
  if (rem > 0L) cs[seq_len(rem)] <- cs[seq_len(rem)] + 1L
  cs
}

# ----------------------------------------------------------------
# Derive per-chunk seeds from a master seed
# ----------------------------------------------------------------

#' @keywords internal
.make_seeds <- function(seed, n_chunks) {
  if (is.null(seed)) return(rep(list(NULL), n_chunks))
  set.seed(seed)
  as.list(sample.int(.Machine$integer.max, n_chunks))
}

# ----------------------------------------------------------------
# Combine accumulator lists by elementwise addition
# ----------------------------------------------------------------

#' @keywords internal
.comb_acc <- function(a, b)
  list(sum_w = a$sum_w + b$sum_w, sum_w2 = a$sum_w2 + b$sum_w2)

# ----------------------------------------------------------------
# Convert raw accumulator totals to p-value estimates + CI
# ----------------------------------------------------------------

#' @keywords internal
## NOTE: two.sided correction is now handled inside each IS simulation loop
## (proper mixture tilt with 2J denominator). This function simply maps
## raw accumulator sums back to the original b scale and computes CIs.
.summarise_acc <- function(acc, K, b_orig, b_sorted, b0) {
  p_sorted   <- acc$sum_w  / K
  E2_sorted  <- acc$sum_w2 / K
  var_sorted <- pmax(0, E2_sorted - p_sorted^2)
  se_sorted  <- sqrt(var_sorted / K)

  inv   <- match(b_orig, b_sorted)
  p_est <- p_sorted[inv]
  p_se  <- se_sorted[inv]

  data.frame(
    p.est    = p_est,
    p.se     = p_se,
    p.lower  = pmax(0, p_est - 1.96 * p_se),
    p.upper  = pmin(1, p_est + 1.96 * p_se),
    b_anchor = b0,
    stringsAsFactors = FALSE
  )
}
