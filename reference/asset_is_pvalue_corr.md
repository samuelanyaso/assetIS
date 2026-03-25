# Gaussian IS Estimator of ASSET Tail Probabilities for Correlated or Overlapping Studies

Estimates \\p(b) = P_0(\max_A Z_A \> b)\\ (one-sided) or \\p(b) =
P_0(\max_A \|Z_A\| \> b)\\ (two-sided) for the ASSET statistic when the
study-level score statistics \\(z_1,\ldots,z_M)\\ are correlated, e.g.,
due to overlapping subjects across studies. The optimal meta-analysis
statistic for subset \\A\\ is \$\$Z_A = \frac{N_A^\top \Sigma_A^{-1}
z_A} {\sqrt{N_A^\top \Sigma_A^{-1} N_A}},\$\$ where \\N_A =
(\sqrt{n_m})\_{m \in A}\\ and \\\Sigma_A\\ is the principal submatrix of
the \\M \times M\\ inter-study correlation matrix \\\Sigma\\. The IS
algorithm tilts the joint MVN distribution via \\dP\_{\xi,A} = \exp(\xi
Z_A - \xi^2/2)\\dP_0\\.

## Usage

``` r
asset_is_pvalue_corr(
  z,
  n,
  Sigma,
  b,
  K = 5000L,
  two.sided = TRUE,
  xi = NULL,
  parallel = FALSE,
  n.cores = parallel::detectCores() - 1L,
  seed = NULL
)
```

## Arguments

- z:

  Numeric vector of length \\M\\. Observed study-level score statistics.

- n:

  Numeric vector of length \\M\\. Study sample sizes. For case-control
  studies supply effective sample sizes \\n_m^{\mathrm{eff}} = n_m
  \phi_m (1-\phi_m)\\.

- Sigma:

  Numeric \\M \times M\\ inter-study correlation matrix under the null.
  Must be symmetric positive definite with unit diagonal.

- b:

  Numeric scalar or vector of threshold values.

- K:

  Positive integer. Number of IS simulations. Default `5000L`.

- two.sided:

  Logical. If `TRUE` (default) compute the two-sided p-value.

- xi:

  Numeric scalar. Tilting parameter. Defaults to `median(b)` (after
  two-sided halving).

- parallel:

  Logical. Default `FALSE`.

- n.cores:

  Positive integer. Default `parallel::detectCores()-1`.

- seed:

  Integer or `NULL`. Default `NULL`.

## Value

A data frame with one row per value in `b` and columns:

- `b`:

  Threshold value (as supplied by the user).

- `p.est`:

  IS point estimate of \\p(b)\\.

- `p.se`:

  Estimated standard error of `p.est`.

- `p.lower`:

  Lower bound of 95\\ \\\hat{p} - 1.96\\\widehat{\mathrm{SE}}\\,
  truncated to \\\[0,1\]\\.

- `p.upper`:

  Upper bound of 95\\ \\\[0,1\]\\.

- `xi`:

  Tilting parameter used.

- `K`:

  Number of IS simulations used.

## Details

Importance Sampling for ASSET: Correlated Studies, Gaussian
Approximation

**Case-control studies.** Replace \\n_m\\ with the effective sample size
\\n_m^{\mathrm{eff}} = n_m \phi_m (1-\phi_m)\\.

**Numerical stability.** If \\\Sigma\\ is not numerically positive
definite, the function attempts regularisation via
[`Matrix::nearPD`](https://rdrr.io/pkg/Matrix/man/nearPD.html) followed
by diagonal jitter, issuing a warning if applied.

For independent studies use
[`asset_is_pvalue_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_ind.md).
For non-Gaussian settings use
[`asset_is_nonnorm_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_ind.md)
or
[`asset_is_nonnorm_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_corr.md).

## References

Bhattacharjee, S. et al. (2012). A subset-based approach improves power
and interpretation for the combined analysis of genetic association
studies of heterogeneous traits. *American Journal of Human Genetics*,
**90**(5), 821–835.

Shi, J. et al. (2007). Importance sampling for estimating p-values in
linkage analysis. *Journal of the American Statistical Association*,
**102**(479), 929–937.

Siegmund, D. (1976). Importance sampling in the Monte Carlo study of
sequential tests. *The Annals of Statistics*, **4**(4), 673–684.

## See also

[`asset_is_pvalue_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_ind.md),
[`asset_is_nonnorm_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_ind.md),
[`asset_is_nonnorm_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_corr.md)

## Examples

``` r
set.seed(42)
M   <- 5
rho <- 0.4
Sigma <- outer(1:M, 1:M, function(i, j) rho^abs(i - j))
z   <- rnorm(M, mean = c(3, 3, 0, 0, 0))
n   <- c(500, 600, 450, 700, 550)

## Two-sided p-value at a single threshold
b_obs <- max(abs(z))
asset_is_pvalue_corr(z = z, n = n, Sigma = Sigma, b = b_obs,
                     K = 1e4, two.sided = TRUE)
#>          b        p.est         p.se     p.lower      p.upper       xi     K
#> 1 4.370958 0.0001826234 4.315504e-06 0.000174165 0.0001910818 4.370958 10000

## One-sided p-values over a grid (range approximation)
b_grid <- seq(3, 6, by = 0.1)
asset_is_pvalue_corr(z = z, n = n, Sigma = Sigma, b = b_grid,
                     K = 5e4, two.sided = FALSE)
#>      b        p.est         p.se      p.lower      p.upper  xi     K
#> 1  3.0 1.277276e-02 2.249459e-04 1.233186e-02 1.321365e-02 4.5 50000
#> 2  3.1 9.510179e-03 1.559336e-04 9.204549e-03 9.815809e-03 4.5 50000
#> 3  3.2 6.918790e-03 1.057077e-04 6.711603e-03 7.125977e-03 4.5 50000
#> 4  3.3 5.066806e-03 7.369355e-05 4.922367e-03 5.211246e-03 4.5 50000
#> 5  3.4 3.682224e-03 5.050924e-05 3.583226e-03 3.781222e-03 4.5 50000
#> 6  3.5 2.632632e-03 3.460781e-05 2.564801e-03 2.700464e-03 4.5 50000
#> 7  3.6 1.877647e-03 2.365651e-05 1.831281e-03 1.924014e-03 4.5 50000
#> 8  3.7 1.322997e-03 1.604062e-05 1.291558e-03 1.354437e-03 4.5 50000
#> 9  3.8 8.955185e-04 1.050624e-05 8.749263e-04 9.161107e-04 4.5 50000
#> 10 3.9 6.060622e-04 6.946086e-06 5.924478e-04 6.196765e-04 4.5 50000
#> 11 4.0 4.087319e-04 4.574842e-06 3.997652e-04 4.176986e-04 4.5 50000
#> 12 4.1 2.781990e-04 3.093825e-06 2.721351e-04 2.842629e-04 4.5 50000
#> 13 4.2 1.827885e-04 1.995726e-06 1.788769e-04 1.867001e-04 4.5 50000
#> 14 4.3 1.197530e-04 1.298324e-06 1.172083e-04 1.222977e-04 4.5 50000
#> 15 4.4 7.916720e-05 8.465111e-07 7.750804e-05 8.082636e-05 4.5 50000
#> 16 4.5 5.151148e-05 5.555746e-07 5.042255e-05 5.260040e-05 4.5 50000
#> 17 4.6 3.300152e-05 3.614495e-07 3.229307e-05 3.370996e-05 4.5 50000
#> 18 4.7 2.090717e-05 2.324954e-07 2.045148e-05 2.136286e-05 4.5 50000
#> 19 4.8 1.295570e-05 1.477401e-07 1.266613e-05 1.324527e-05 4.5 50000
#> 20 4.9 8.109483e-06 9.478165e-08 7.923711e-06 8.295255e-06 4.5 50000
#> 21 5.0 4.905760e-06 5.954331e-08 4.789055e-06 5.022464e-06 4.5 50000
#> 22 5.1 2.931661e-06 3.693001e-08 2.859278e-06 3.004044e-06 4.5 50000
#> 23 5.2 1.790689e-06 2.355934e-08 1.744512e-06 1.836865e-06 4.5 50000
#> 24 5.3 1.056758e-06 1.469950e-08 1.027947e-06 1.085569e-06 4.5 50000
#> 25 5.4 6.226540e-07 9.083527e-09 6.048502e-07 6.404577e-07 4.5 50000
#> 26 5.5 3.579976e-07 5.572035e-09 3.470764e-07 3.689188e-07 4.5 50000
#> 27 5.6 2.070728e-07 3.458422e-09 2.002943e-07 2.138514e-07 4.5 50000
#> 28 5.7 1.152491e-07 2.076838e-09 1.111785e-07 1.193197e-07 4.5 50000
#> 29 5.8 6.443466e-08 1.255790e-09 6.197331e-08 6.689601e-08 4.5 50000
#> 30 5.9 3.595407e-08 7.585844e-10 3.446725e-08 3.744090e-08 4.5 50000
#> 31 6.0 1.891556e-08 4.397079e-10 1.805373e-08 1.977739e-08 4.5 50000

## Case-control: supply effective sample sizes
phi   <- c(0.3, 0.4, 0.5, 0.35, 0.45)
n_eff <- n * phi * (1 - phi)
asset_is_pvalue_corr(z = z, n = n_eff, Sigma = Sigma, b = b_obs,
                     K = 1e4, two.sided = TRUE)
#>          b        p.est         p.se     p.lower      p.upper       xi     K
#> 1 4.370958 0.0001789857 4.321798e-06 0.000170515 0.0001874564 4.370958 10000

## Parallel execution
if (FALSE) { # \dontrun{
asset_is_pvalue_corr(z = z, n = n, Sigma = Sigma, b = b_grid,
                     K = 1e5, parallel = TRUE, n.cores = 4)
} # }
```
