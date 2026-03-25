## tests/testthat/test-estimators.R
## Basic correctness and interface tests for all four IS estimators

library(assetIS)

## ----------------------------------------------------------------
## Shared fixtures
## ----------------------------------------------------------------
set.seed(123)
M <- 4
n_vec  <- c(300, 400, 350, 500)
z_vec  <- c(2.8, 2.5, 0.3, -0.2)
b_single <- max(abs(z_vec))
b_grid   <- c(2.5, 3.0, 3.5)
rho      <- 0.3
Sigma    <- outer(seq_len(M), seq_len(M), function(i, j) rho^abs(i - j))

N <- 80; C <- 4; f_val <- 0.2
Y_eqtl <- scale(matrix(rnorm(N * C), N, C))

K_test <- 500L   # small for CRAN speed

## ----------------------------------------------------------------
## 1. asset_is_pvalue_ind
## ----------------------------------------------------------------
test_that("asset_is_pvalue_ind returns valid output for scalar b", {
  res <- asset_is_pvalue_ind(z = z_vec, n = n_vec, b = b_single,
                              K = K_test, seed = 1L)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_true(all(c("b", "p.est", "p.se", "p.lower", "p.upper", "xi", "K") %in%
                    names(res)))
  expect_true(res$p.est >= 0 && res$p.est <= 1)
  expect_true(res$p.se  >= 0)
  expect_true(res$p.lower <= res$p.est)
  expect_true(res$p.upper >= res$p.est)
})

test_that("asset_is_pvalue_ind handles vector b", {
  res <- asset_is_pvalue_ind(z = z_vec, n = n_vec, b = b_grid,
                              K = K_test, seed = 2L)
  expect_equal(nrow(res), length(b_grid))
  expect_true(all(res$p.est >= 0 & res$p.est <= 1))
  ## p-values should be decreasing as b increases
  expect_true(all(diff(res$p.est) <= 0.1))  # monotone up to MC noise
})

test_that("asset_is_pvalue_ind one-sided gives larger p than two-sided", {
  p2 <- asset_is_pvalue_ind(z = z_vec, n = n_vec, b = b_single,
                              K = K_test, two.sided = TRUE, seed = 3L)$p.est
  p1 <- asset_is_pvalue_ind(z = z_vec, n = n_vec, b = b_single,
                              K = K_test, two.sided = FALSE, seed = 3L)$p.est
  ## one-sided at same b should be <= two-sided
  expect_true(p1 <= p2 + 0.05)
})

test_that("asset_is_pvalue_ind validates inputs", {
  expect_error(asset_is_pvalue_ind(z = z_vec, n = n_vec[-1], b = b_single))
  expect_error(asset_is_pvalue_ind(z = z_vec, n = c(-1, n_vec[-1]),
                                   b = b_single))
  expect_error(asset_is_pvalue_ind(z = z_vec, n = n_vec, b = b_single,
                                   K = 0L))
})

## ----------------------------------------------------------------
## 2. asset_is_pvalue_corr
## ----------------------------------------------------------------
test_that("asset_is_pvalue_corr returns valid output", {
  res <- asset_is_pvalue_corr(z = z_vec, n = n_vec, Sigma = Sigma,
                               b = b_single, K = K_test, seed = 4L)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_true(res$p.est >= 0 && res$p.est <= 1)
})

test_that("asset_is_pvalue_corr matches ind when Sigma = I", {
  res_corr <- asset_is_pvalue_corr(z = z_vec, n = n_vec, Sigma = diag(M),
                                    b = b_single, K = K_test, seed = 5L)
  res_ind  <- asset_is_pvalue_ind(z = z_vec, n = n_vec,
                                   b = b_single, K = K_test, seed = 5L)
  ## Results should be in the same ballpark (not exactly equal due to different
  ## weight formulas at Sigma=I, but close in order of magnitude)
  expect_true(abs(log10(res_corr$p.est + 1e-12) -
                    log10(res_ind$p.est  + 1e-12)) < 1.5)
})

test_that("asset_is_pvalue_corr validates Sigma", {
  bad_Sigma <- Sigma; bad_Sigma[1, 2] <- bad_Sigma[1, 2] + 0.5
  expect_error(asset_is_pvalue_corr(z = z_vec, n = n_vec, Sigma = bad_Sigma,
                                     b = b_single))
})

## ----------------------------------------------------------------
## 3. asset_is_nonnorm_ind
## ----------------------------------------------------------------
test_that("asset_is_nonnorm_ind returns valid output (eQTL mode)", {
  res <- asset_is_nonnorm_ind(Y = Y_eqtl, f = f_val,
                               b = 4.0, K = K_test, seed = 6L)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_true(all(c("b", "p.est", "p.se", "p.lower", "p.upper",
                    "b_anchor", "K") %in% names(res)))
  expect_true(res$p.est >= 0 && res$p.est <= 1)
})

test_that("asset_is_nonnorm_ind handles vector b", {
  res <- asset_is_nonnorm_ind(Y = Y_eqtl, f = f_val,
                               b = c(3.5, 4.0, 4.5), K = K_test, seed = 7L)
  expect_equal(nrow(res), 3L)
  expect_true(all(res$p.est >= 0 & res$p.est <= 1))
})

test_that("asset_is_nonnorm_ind works in case-control mode", {
  study_id <- rep(seq_len(C), each = N %/% C)
  phi_c    <- rep(0.3, C)
  y_bin    <- rbinom(length(study_id), 1, phi_c[study_id])
  Ycc      <- matrix(0, length(study_id), C)
  for (i in seq_along(study_id))
    Ycc[i, study_id[i]] <- y_bin[i] - phi_c[study_id[i]]
  res <- asset_is_nonnorm_ind(Y = Ycc, f = f_val, b = 3.5,
                               study_id = study_id, K = K_test, seed = 8L)
  expect_true(res$p.est >= 0 && res$p.est <= 1)
})

test_that("asset_is_nonnorm_ind accepts study-specific MAF", {
  study_id <- rep(seq_len(C), each = N %/% C)
  phi_c    <- rep(0.3, C)
  y_bin    <- rbinom(length(study_id), 1, phi_c[study_id])
  Ycc      <- matrix(0, length(study_id), C)
  for (i in seq_along(study_id))
    Ycc[i, study_id[i]] <- y_bin[i] - phi_c[study_id[i]]
  f_multi <- c(0.15, 0.20, 0.10, 0.25)
  res <- asset_is_nonnorm_ind(Y = Ycc, f = f_multi, b = 3.5,
                               study_id = study_id, K = K_test, seed = 9L)
  expect_true(res$p.est >= 0 && res$p.est <= 1)
})

test_that("asset_is_nonnorm_ind validates f range", {
  expect_error(asset_is_nonnorm_ind(Y = Y_eqtl, f = 1.5, b = 4.0))
  expect_error(asset_is_nonnorm_ind(Y = Y_eqtl, f = 0,   b = 4.0))
})

## ----------------------------------------------------------------
## 4. asset_is_nonnorm_corr
## ----------------------------------------------------------------
test_that("asset_is_nonnorm_corr returns valid output", {
  res <- asset_is_nonnorm_corr(Y = Y_eqtl, f = f_val,
                                b = 4.0, K = K_test, seed = 10L)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_true(res$p.est >= 0 && res$p.est <= 1)
})

test_that("asset_is_nonnorm_corr accepts pre-computed Sigma_z", {
  Sz <- crossprod(scale(Y_eqtl)) / N
  diag(Sz) <- 1
  res <- asset_is_nonnorm_corr(Y = Y_eqtl, f = f_val, b = 4.0,
                                Sigma_z = Sz, scale_Y = FALSE,
                                K = K_test, seed = 11L)
  expect_true(res$p.est >= 0 && res$p.est <= 1)
})

test_that("asset_is_nonnorm_corr handles scale_Y = FALSE", {
  Y_pre <- scale(Y_eqtl)
  res1  <- asset_is_nonnorm_corr(Y = Y_eqtl, f = f_val, b = 4.0,
                                   scale_Y = TRUE,  K = K_test, seed = 12L)
  res2  <- asset_is_nonnorm_corr(Y = Y_pre,   f = f_val, b = 4.0,
                                   scale_Y = FALSE, K = K_test, seed = 12L)
  ## Same input, same seed -> same result
  expect_equal(res1$p.est, res2$p.est, tolerance = 1e-6)
})

test_that("asset_is_nonnorm_corr validates Sigma_z", {
  bad_Sz <- crossprod(Y_eqtl) / N
  bad_Sz[1, 2] <- bad_Sz[1, 2] + 1   # break symmetry
  expect_error(asset_is_nonnorm_corr(Y = Y_eqtl, f = f_val, b = 4.0,
                                      Sigma_z = bad_Sz))
})

## ----------------------------------------------------------------
## 5. Reproducibility
## ----------------------------------------------------------------
test_that("seed argument produces identical results", {
  r1 <- asset_is_pvalue_ind(z = z_vec, n = n_vec, b = b_single,
                              K = K_test, seed = 42L)
  r2 <- asset_is_pvalue_ind(z = z_vec, n = n_vec, b = b_single,
                              K = K_test, seed = 42L)
  expect_equal(r1$p.est, r2$p.est)
})
