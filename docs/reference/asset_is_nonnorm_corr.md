# Conditional IS Estimator of ASSET Tail Probabilities for Non-Gaussian Settings with Correlated Score Statistics

Estimates \\p(b \mid Y) = P_0(\max_A Z_A \> b \mid Y)\\ (one-sided) or
\\p(b \mid Y) = P_0(\max_A \|Z_A\| \> b \mid Y)\\ (two-sided) when both
the normality assumption may fail and score statistics are correlated,
e.g., due to correlated gene expression across cell types or overlapping
subjects in case-control meta-analysis.

The conditional covariance of the score vector \\\mathbf{z}\\ is
estimated as \\\widehat\Sigma_Y = Y^\top Y / N\\ (after scaling). The
optimal GLS statistic for subset \\A\\ is \$\$Z_A =
\frac{\mathbf{1}\_A^\top \Sigma_A^{-1} \mathbf{z}\_A}
{\sqrt{\mathbf{1}\_A^\top \Sigma_A^{-1} \mathbf{1}\_A}},\$\$ which
admits a linear genotype representation \\Z_A = \sum_i
\omega\_{i,A}(g_i - \tilde g_i)\\ with GLS-adjusted weights. The IS
algorithm then follows
[`asset_is_nonnorm_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_ind.md):
Newton solve for \\\xi_A\\, trinomial genotype simulation under the
tilted HWE distribution, and the mixture IS weight \\(2^C-1) / \sum_A
\exp\\\xi_A Z_A - \phi_A(\xi_A)\\\\.

## Usage

``` r
asset_is_nonnorm_corr(
  Y,
  f,
  b,
  study_id = NULL,
  g_center = NULL,
  tau = NULL,
  Sigma_z = NULL,
  scale_Y = TRUE,
  ridge = 1e-10,
  two.sided = TRUE,
  K = 5000L,
  b_anchor = NULL,
  newton_tol = 1e-10,
  newton_max_iter = 100L,
  parallel = FALSE,
  n.cores = parallel::detectCores() - 1L,
  seed = NULL
)
```

## Arguments

- Y:

  Numeric \\N \times C\\ matrix. See
  [`asset_is_nonnorm_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_ind.md)
  for mode-specific conventions. Scaled internally if `scale_Y = TRUE`.

- f:

  Numeric scalar or length-\\C\\ vector. Minor allele frequency (MAF). A
  scalar applies a common MAF to all cell types / studies; a vector
  applies study-specific MAFs (multi-ancestry case-control). All values
  must be in \\(0, 0.5\]\\.

- b:

  Numeric scalar or vector of threshold values.

- study_id:

  Integer vector of length \\N\\. Activates case-control mode. Default
  `NULL`.

- g_center:

  Numeric vector of length \\N\\. Per-subject genotype centering.
  Default `NULL` uses \\2f\_{c_i}\\.

- tau:

  Numeric positive vector of length \\C\\. Scale factors for
  covariate-adjusted eQTL. Default `NULL`.

- Sigma_z:

  Numeric \\C \times C\\ matrix. Conditional covariance (correlation)
  matrix of \\(z_1,\ldots,z_C)\\ given \\Y\\. Defaults to `NULL`, in
  which case \\\widehat{\Sigma}\_Y = Y^\top Y / N\\ is computed
  internally after scaling. Supply a pre-computed or externally
  estimated correlation matrix to override the default (e.g., from a
  previous run, or estimated from a larger reference panel). Must be
  symmetric positive definite with unit diagonal.

- scale_Y:

  Logical. If `TRUE` (default), `Y` is centered and scaled to mean zero
  and unit variance per column. Set to `FALSE` only if `Y` has already
  been pre-processed.

- ridge:

  Positive numeric. Ridge added to \\\Sigma_A\\ before inversion.
  Default `1e-10`.

- two.sided:

  Logical. Default `TRUE`.

- K:

  Positive integer. Number of IS simulations. Default `5000L`.

- b_anchor:

  Numeric scalar. The anchor threshold \\b_0\\ at which \\\xi\_{A,\pm}\\
  are solved. Defaults to `NULL`, using \\b_0 = \mathrm{median}(b)\\.
  Ignored for scalar `b`.

- newton_tol:

  Numeric. Newton tolerance. Default `1e-10`.

- newton_max_iter:

  Integer. Maximum Newton iterations. Default `100L`.

- parallel:

  Logical. Default `FALSE`.

- n.cores:

  Positive integer. Default `parallel::detectCores()-1`.

- seed:

  Integer or `NULL`. Default `NULL`.

## Value

A data frame with one row per threshold in `b` and columns:

- `b`:

  Threshold value as supplied.

- `p.est`:

  IS point estimate of \\p(b \mid Y)\\.

- `p.se`:

  Estimated standard error of `p.est`.

- `p.lower`:

  Lower bound of 95\\ \\\[0,1\]\\.

- `p.upper`:

  Upper bound of 95\\ \\\[0,1\]\\.

- `b_anchor`:

  Anchor threshold used for Newton solve.

- `K`:

  Number of IS simulations used.

## Details

Conditional IS for ASSET When Normality Fails: Correlated Cell Types or
Overlapping Studies

**Covariate adjustment.** When `tau` and `g_center` are supplied, the
adjusted omega formula \\\omega\_{i,A} = N^{-1/2}\sum\_{c \in
A}a\_{c,A}\\y'\_{ic}/\tau_c / \sqrt{\mathbf{1}\_A^\top a_A}\\ is used,
where \\y'\_{ic}\\ are expression residuals after regressing out
covariates.

**Scaling.** The formula \\\widehat\Sigma_Y = Y^\top Y/N\\ requires
zero-mean, unit-variance columns. Set `scale_Y = TRUE` (default) to
apply [`base::scale()`](https://rdrr.io/r/base/scale.html) internally.

**Supplying `Sigma_z`.** A pre-computed or reference-panel correlation
matrix can be passed directly via `Sigma_z`, bypassing the empirical
estimate.

**Modes.** The same three analysis modes as
[`asset_is_nonnorm_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_ind.md)
are supported (eQTL, covariate-adjusted eQTL, case-control), now with
GLS weights derived from \\\widehat\Sigma_Y\\.

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

[`asset_is_nonnorm_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_ind.md),
[`asset_is_pvalue_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_corr.md)

## Examples

``` r
set.seed(42)
N <- 200; C <- 5; f <- 0.2; rho <- 0.4

## Correlated expression data (compound symmetry)
Sigma_expr <- (1 - rho) * diag(C) + rho
Y <- MASS::mvrnorm(N, mu = rep(0, C), Sigma = Sigma_expr)

## --- eQTL mode (correlated cell types) ---
b_obs <- 4.5
asset_is_nonnorm_corr(Y = Y, f = f, b = b_obs, K = 2000L, two.sided = TRUE)
#>     b        p.est         p.se      p.lower      p.upper b_anchor    K
#> 1 4.5 7.987572e-05 4.866892e-06 7.033661e-05 8.941483e-05      4.5 2000

## Range of thresholds
b_grid <- seq(3.5, 6.0, by = 0.1)
asset_is_nonnorm_corr(Y = Y, f = f, b = b_grid, K = 1e4, two.sided = TRUE)
#>      b        p.est         p.se      p.lower      p.upper b_anchor     K
#> 1  3.5 4.054243e-03 1.522211e-04 3.755890e-03 4.352597e-03     4.75 10000
#> 2  3.6 2.945766e-03 1.049116e-04 2.740140e-03 3.151393e-03     4.75 10000
#> 3  3.7 1.965829e-03 6.366422e-05 1.841047e-03 2.090611e-03     4.75 10000
#> 4  3.8 1.410119e-03 4.392195e-05 1.324032e-03 1.496206e-03     4.75 10000
#> 5  3.9 9.632054e-04 2.861356e-05 9.071228e-04 1.019288e-03     4.75 10000
#> 6  4.0 6.800352e-04 1.929261e-05 6.422217e-04 7.178487e-04     4.75 10000
#> 7  4.1 4.486800e-04 1.272169e-05 4.237455e-04 4.736145e-04     4.75 10000
#> 8  4.2 2.926937e-04 8.085187e-06 2.768468e-04 3.085407e-04     4.75 10000
#> 9  4.3 1.878238e-04 5.076626e-06 1.778736e-04 1.977740e-04     4.75 10000
#> 10 4.4 1.286189e-04 3.373030e-06 1.220078e-04 1.352301e-04     4.75 10000
#> 11 4.5 8.155001e-05 2.086892e-06 7.745970e-05 8.564032e-05     4.75 10000
#> 12 4.6 5.349119e-05 1.359545e-06 5.082648e-05 5.615590e-05     4.75 10000
#> 13 4.7 3.448797e-05 8.640196e-07 3.279449e-05 3.618144e-05     4.75 10000
#> 14 4.8 2.145417e-05 5.499406e-07 2.037629e-05 2.253205e-05     4.75 10000
#> 15 4.9 1.326585e-05 3.459355e-07 1.258781e-05 1.394388e-05     4.75 10000
#> 16 5.0 8.288391e-06 2.162345e-07 7.864572e-06 8.712211e-06     4.75 10000
#> 17 5.1 5.114105e-06 1.351706e-07 4.849171e-06 5.379040e-06     4.75 10000
#> 18 5.2 2.986760e-06 8.327817e-08 2.823535e-06 3.149985e-06     4.75 10000
#> 19 5.3 1.768281e-06 5.031768e-08 1.669658e-06 1.866903e-06     4.75 10000
#> 20 5.4 1.040238e-06 3.081215e-08 9.798458e-07 1.100629e-06     4.75 10000
#> 21 5.5 6.070955e-07 1.900854e-08 5.698387e-07 6.443522e-07     4.75 10000
#> 22 5.6 3.599241e-07 1.196466e-08 3.364733e-07 3.833748e-07     4.75 10000
#> 23 5.7 2.000200e-07 7.034092e-09 1.862332e-07 2.138068e-07     4.75 10000
#> 24 5.8 1.069890e-07 3.995038e-09 9.915870e-08 1.148193e-07     4.75 10000
#> 25 5.9 6.202141e-08 2.499719e-09 5.712197e-08 6.692086e-08     4.75 10000
#> 26 6.0 3.314947e-08 1.485981e-09 3.023695e-08 3.606199e-08     4.75 10000

## Supply pre-computed Sigma_z (e.g., from a reference panel)
Sigma_z_ref <- cor(Y)
asset_is_nonnorm_corr(Y = Y, f = f, b = b_obs, Sigma_z = Sigma_z_ref, K = 2000L)
#>     b       p.est         p.se      p.lower      p.upper b_anchor    K
#> 1 4.5 8.11448e-05 4.540486e-06 7.224545e-05 9.004416e-05      4.5 2000

## --- eQTL with covariate adjustment ---
tau_c     <- apply(scale(Y), 2, sd) * sqrt(2 * f * (1 - f))
g_center  <- rep(2 * f, N)
asset_is_nonnorm_corr(Y = Y, f = f, b = b_obs,
                      tau = tau_c, g_center = g_center, K = 2000L)
#>     b        p.est         p.se      p.lower      p.upper b_anchor    K
#> 1 4.5 8.425665e-05 4.826406e-06 7.479689e-05 9.371641e-05      4.5 2000

## --- Case-control with sample overlap (same ancestry) ---
study_id <- sample(1:C, N, replace = TRUE)
phi_c    <- runif(C, 0.2, 0.4)
y_raw    <- rbinom(N, 1, prob = phi_c[study_id])
Ycc      <- matrix(0, N, C)
for (i in seq_len(N)) Ycc[i, study_id[i]] <- y_raw[i] - phi_c[study_id[i]]
asset_is_nonnorm_corr(Y = Ycc, f = f, b = b_obs,
                      study_id = study_id, K = 2000L)
#>     b        p.est         p.se      p.lower      p.upper b_anchor    K
#> 1 4.5 1.095769e-05 4.931541e-06 1.291874e-06 2.062352e-05      4.5 2000

## --- Multi-ancestry case-control ---
f_c <- runif(C, 0.1, 0.4)
asset_is_nonnorm_corr(Y = Ycc, f = f_c, b = b_obs,
                      study_id = study_id, K = 2000L)
#>     b        p.est         p.se      p.lower      p.upper b_anchor    K
#> 1 4.5 6.936307e-06 2.938049e-06 1.177731e-06 1.269488e-05      4.5 2000

## --- Parallel execution ---
if (FALSE) { # \dontrun{
asset_is_nonnorm_corr(Y = Y, f = f, b = b_grid,
                      K = 1e5, parallel = TRUE, n.cores = 4)
} # }
```
