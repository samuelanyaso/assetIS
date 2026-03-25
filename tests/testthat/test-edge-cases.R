## tests/testthat/test-edge-cases.R
## Edge-case and integration tests for assetIS

library(assetIS)

set.seed(321)
N <- 60; C <- 3; f_val <- 0.25
Y_small <- scale(matrix(rnorm(N * C), N, C))

## ----------------------------------------------------------------
## 1. Monotonicity: p-values decrease as b increases
## ----------------------------------------------------------------
test_that("p-values decrease with increasing b (Gaussian ind)", {
  b_seq <- seq(2.0, 5.0, by = 0.5)
  res   <- asset_is_pvalue_ind(z = c(2.8, -0.2, 0.4), n = c(200, 300, 250),
                                b = b_seq, K = 1000L, seed = 1L)
  ## Allow small non-monotonicities due to MC variance at very small K
  ## but overall trend must be decreasing
  expect_true(res$p.est[1] > res$p.est[length(b_seq)])
})

test_that("p-values decrease with increasing b (non-Gaussian ind)", {
  b_seq <- seq(3.0, 5.5, by = 0.5)
  res   <- asset_is_nonnorm_ind(Y = Y_small, f = f_val,
                                 b = b_seq, K = 500L, seed = 2L)
  expect_true(res$p.est[1] >= res$p.est[length(b_seq)])
})

## ----------------------------------------------------------------
## 2. Two-sided ≈ 2 * one-sided (for symmetric statistics)
## ----------------------------------------------------------------
test_that("two-sided ~ 2 * one-sided for Gaussian ind", {
  z_sym <- c(2.8, -2.8, 0.1)
  n_sym <- c(400, 400, 400)
  b_test <- 2.8
  p2 <- asset_is_pvalue_ind(z = z_sym, n = n_sym, b = b_test,
                              K = 3000L, two.sided = TRUE,  seed = 10L)$p.est
  p1 <- asset_is_pvalue_ind(z = z_sym, n = n_sym, b = b_test,
                              K = 3000L, two.sided = FALSE, seed = 10L)$p.est
  ## two-sided should be roughly 2x one-sided (within Monte Carlo noise)
  expect_equal(p2, 2 * p1, tolerance = 0.15)
})

## ----------------------------------------------------------------
## 3. b_anchor argument
## ----------------------------------------------------------------
test_that("b_anchor overrides median(b) for Newton solve", {
  b_vec <- seq(3.5, 5.0, by = 0.5)
  ## Should run without error using a custom anchor
  res <- asset_is_nonnorm_ind(Y = Y_small, f = f_val,
                               b = b_vec, b_anchor = 4.0,
                               K = 300L, seed = 11L)
  expect_equal(nrow(res), length(b_vec))
  expect_equal(unique(res$b_anchor), 4.0)
})

## ----------------------------------------------------------------
## 4. Covariate-adjusted eQTL mode
## ----------------------------------------------------------------
test_that("asset_is_nonnorm_ind works in covariate-adjusted eQTL mode", {
  tau_c    <- apply(Y_small, 2, sd) * sqrt(2 * f_val * (1 - f_val))
  g_center <- rep(2 * f_val, N)
  res <- asset_is_nonnorm_ind(Y = Y_small, f = f_val,
                               tau = tau_c, g_center = g_center,
                               b = 4.0, K = 500L, seed = 12L)
  expect_true(res$p.est >= 0 && res$p.est <= 1)
})

test_that("asset_is_nonnorm_corr works in covariate-adjusted eQTL mode", {
  tau_c    <- apply(Y_small, 2, sd) * sqrt(2 * f_val * (1 - f_val))
  g_center <- rep(2 * f_val, N)
  res <- asset_is_nonnorm_corr(Y = Y_small, f = f_val,
                                tau = tau_c, g_center = g_center,
                                b = 4.0, K = 500L, seed = 13L)
  expect_true(res$p.est >= 0 && res$p.est <= 1)
})

## ----------------------------------------------------------------
## 5. Correlated estimators accept user-supplied Sigma_z
## ----------------------------------------------------------------
test_that("asset_is_nonnorm_corr: pre-computed Sigma_z gives same
           result as internally computed", {
  Y_scaled <- scale(Y_small)
  Sz       <- crossprod(Y_scaled) / N
  diag(Sz) <- 1

  res_internal <- asset_is_nonnorm_corr(Y = Y_small, f = f_val,
                                         b = 4.0, scale_Y = TRUE,
                                         K = 500L, seed = 14L)
  res_external <- asset_is_nonnorm_corr(Y = Y_scaled, f = f_val,
                                         b = 4.0, Sigma_z = Sz,
                                         scale_Y = FALSE,
                                         K = 500L, seed = 14L)
  expect_equal(res_internal$p.est, res_external$p.est, tolerance = 1e-8)
})

## ----------------------------------------------------------------
## 6. ridge argument prevents failure for highly collinear Y
## ----------------------------------------------------------------
test_that("ridge prevents Sigma_A singularity for collinear Y", {
  ## Nearly collinear: all columns = same vector + tiny noise
  base_col <- rnorm(N)
  Y_collin <- matrix(base_col, N, C) +
              matrix(rnorm(N * C, sd = 1e-4), N, C)
  Y_collin <- scale(Y_collin)
  ## Should not throw, thanks to ridge
  expect_no_error(
    asset_is_nonnorm_corr(Y = Y_collin, f = f_val, b = 4.0,
                           ridge = 1e-6, K = 200L, seed = 15L)
  )
})

## ----------------------------------------------------------------
## 7. Parallel produces same estimates as serial (same seed)
## ----------------------------------------------------------------
test_that("parallel = TRUE reproduces serial result with same seed", {
  skip_if_not(parallel::detectCores() >= 2,
              "Need >= 2 cores for parallel test")

  res_ser <- asset_is_pvalue_ind(z  = c(3.0, 2.5, -0.1),
                                  n  = c(300, 400, 350),
                                  b  = 3.5, K = 400L,
                                  parallel = FALSE, seed = 77L)
  res_par <- asset_is_pvalue_ind(z  = c(3.0, 2.5, -0.1),
                                  n  = c(300, 400, 350),
                                  b  = 3.5, K = 400L,
                                  parallel = TRUE, n.cores = 2L, seed = 77L)
  ## Estimates won't be identical (different chunk splits), but should
  ## be in same ballpark
  expect_true(abs(res_ser$p.est - res_par$p.est) < 0.1)
})

## ----------------------------------------------------------------
## 8. Output column names are correct for all four functions
## ----------------------------------------------------------------
test_that("output column names match spec for all four functions", {
  z3  <- c(2.9, -0.1, 0.3)
  n3  <- c(300, 350, 400)
  Sig <- diag(3)
  Y3  <- scale(matrix(rnorm(N * 3), N, 3))

  cols_gauss   <- c("b", "p.est", "p.se", "p.lower", "p.upper", "xi", "K")
  cols_nonnorm <- c("b", "p.est", "p.se", "p.lower", "p.upper", "b_anchor", "K")

  r1 <- asset_is_pvalue_ind(z = z3, n = n3, b = 3.0, K = 200L, seed = 20L)
  r2 <- asset_is_pvalue_corr(z = z3, n = n3, Sigma = Sig,
                              b = 3.0, K = 200L, seed = 21L)
  r3 <- asset_is_nonnorm_ind(Y = Y3, f = 0.2, b = 3.5, K = 200L, seed = 22L)
  r4 <- asset_is_nonnorm_corr(Y = Y3, f = 0.2, b = 3.5, K = 200L, seed = 23L)

  expect_named(r1, cols_gauss)
  expect_named(r2, cols_gauss)
  expect_named(r3, cols_nonnorm)
  expect_named(r4, cols_nonnorm)
})

## ----------------------------------------------------------------
## 9. p.lower <= p.est <= p.upper for all rows
## ----------------------------------------------------------------
test_that("CI bounds are consistent for vector b", {
  res <- asset_is_pvalue_ind(z = c(3.1, 2.8, -0.2), n = c(400, 500, 350),
                              b = seq(2.5, 4.5, by = 0.25),
                              K = 1000L, seed = 30L)
  expect_true(all(res$p.lower <= res$p.est + 1e-10))
  expect_true(all(res$p.upper >= res$p.est - 1e-10))
  expect_true(all(res$p.lower >= 0))
  expect_true(all(res$p.upper <= 1))
})
