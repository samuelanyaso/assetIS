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
#> 1 4.5 8.610092e-05 4.997332e-06 7.630615e-05 9.589569e-05      4.5 2000

## Range of thresholds
b_grid <- seq(3.5, 6.0, by = 0.1)
asset_is_nonnorm_corr(Y = Y, f = f, b = b_grid, K = 1e4, two.sided = TRUE)
#>      b        p.est         p.se      p.lower      p.upper b_anchor     K
#> 1  3.5 4.335015e-03 1.550822e-04 4.031054e-03 4.638976e-03     4.75 10000
#> 2  3.6 3.085065e-03 1.073403e-04 2.874678e-03 3.295452e-03     4.75 10000
#> 3  3.7 2.067795e-03 6.795803e-05 1.934598e-03 2.200993e-03     4.75 10000
#> 4  3.8 1.435726e-03 4.436274e-05 1.348775e-03 1.522677e-03     4.75 10000
#> 5  3.9 1.013705e-03 3.045905e-05 9.540057e-04 1.073405e-03     4.75 10000
#> 6  4.0 6.799915e-04 2.019476e-05 6.404098e-04 7.195732e-04     4.75 10000
#> 7  4.1 4.356949e-04 1.237558e-05 4.114388e-04 4.599510e-04     4.75 10000
#> 8  4.2 2.895161e-04 8.071528e-06 2.736959e-04 3.053363e-04     4.75 10000
#> 9  4.3 1.948310e-04 5.205050e-06 1.846291e-04 2.050329e-04     4.75 10000
#> 10 4.4 1.255303e-04 3.305402e-06 1.190517e-04 1.320089e-04     4.75 10000
#> 11 4.5 8.177838e-05 2.139314e-06 7.758532e-05 8.597143e-05     4.75 10000
#> 12 4.6 5.322616e-05 1.351750e-06 5.057673e-05 5.587559e-05     4.75 10000
#> 13 4.7 3.448518e-05 8.928358e-07 3.273522e-05 3.623514e-05     4.75 10000
#> 14 4.8 2.082457e-05 5.385499e-07 1.976901e-05 2.188012e-05     4.75 10000
#> 15 4.9 1.287301e-05 3.381560e-07 1.221023e-05 1.353580e-05     4.75 10000
#> 16 5.0 8.119459e-06 2.221589e-07 7.684027e-06 8.554890e-06     4.75 10000
#> 17 5.1 4.783573e-06 1.344944e-07 4.519964e-06 5.047183e-06     4.75 10000
#> 18 5.2 2.882920e-06 8.209352e-08 2.722017e-06 3.043823e-06     4.75 10000
#> 19 5.3 1.741594e-06 5.037503e-08 1.642859e-06 1.840329e-06     4.75 10000
#> 20 5.4 1.069962e-06 3.288892e-08 1.005500e-06 1.134424e-06     4.75 10000
#> 21 5.5 6.285259e-07 1.981317e-08 5.896921e-07 6.673598e-07     4.75 10000
#> 22 5.6 3.685298e-07 1.236049e-08 3.443032e-07 3.927563e-07     4.75 10000
#> 23 5.7 1.954922e-07 7.042299e-09 1.816893e-07 2.092951e-07     4.75 10000
#> 24 5.8 1.133610e-07 4.309343e-09 1.049147e-07 1.218073e-07     4.75 10000
#> 25 5.9 6.160586e-08 2.509946e-09 5.668637e-08 6.652535e-08     4.75 10000
#> 26 6.0 3.566073e-08 1.579081e-09 3.256573e-08 3.875573e-08     4.75 10000

## Supply pre-computed Sigma_z (e.g., from a reference panel)
Sigma_z_ref <- cor(Y)
asset_is_nonnorm_corr(Y = Y, f = f, b = b_obs, Sigma_z = Sigma_z_ref, K = 2000L)
#>     b        p.est         p.se      p.lower    p.upper b_anchor    K
#> 1 4.5 8.673704e-05 4.882125e-06 7.716807e-05 9.6306e-05      4.5 2000

## --- eQTL with covariate adjustment ---
tau_c     <- apply(scale(Y), 2, sd) * sqrt(2 * f * (1 - f))
g_center  <- rep(2 * f, N)
asset_is_nonnorm_corr(Y = Y, f = f, b = b_obs,
                      tau = tau_c, g_center = g_center, K = 2000L)
#>     b       p.est         p.se      p.lower     p.upper b_anchor    K
#> 1 4.5 8.80634e-05 5.446577e-06 7.738811e-05 9.87387e-05      4.5 2000

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
