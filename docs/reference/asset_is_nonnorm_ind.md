# Conditional IS Estimator of ASSET Tail Probabilities for Non-Gaussian, Independent Settings

Estimates the conditional tail probability \\p(b \mid Y) = P_0(\max_A
Z_A \> b \mid Y)\\ (one-sided) or \\p(b \mid Y) = P_0(\max_A \|Z_A\| \>
b \mid Y)\\ (two-sided) for the ASSET statistic when the normality
assumption on score statistics may be violated and cell types / studies
are independent. This arises in single-cell eQTL mapping with small
sample sizes or meta-analysis of rare diseases.

Inference is conditional on the observed phenotype matrix \\Y\\; the
only source of randomness is the genotype vector \\(g_1,\ldots,g_N)\\.
Genotypes are simulated from a discrete (trinomial) distribution under
Hardy–Weinberg equilibrium (HWE), with per-subject exponential tilting.
The tilting parameter \\\xi_A\\ for each subset is chosen to satisfy
\\\phi_A'(\xi_A) = b_0\\ via Newton's method.

The function supports three analysis modes selected automatically:

- eQTL:

  All \\N\\ subjects in all \\C\\ cell types (`study_id = NULL`,
  `tau = NULL`).

- eQTL with covariate adjustment:

  `tau` and `g_center` supplied; `Y` contains expression residuals.

- Case-control meta-analysis:

  `study_id` supplied; subjects nested within \\C\\ studies;
  multi-ancestry via a length-\\C\\ vector `f`.

## Usage

``` r
asset_is_nonnorm_ind(
  Y,
  f,
  b,
  study_id = NULL,
  g_center = NULL,
  tau = NULL,
  scale_Y = TRUE,
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

  Numeric matrix of dimension \\N \times C\\, centered and scaled (mean
  zero, unit variance per column) or raw (if `scale_Y = TRUE`). For
  eQTL: the expression matrix (or residual matrix after covariate
  adjustment, with `tau` supplied). For case-control: a matrix whose
  entry `Y[i, study_id[i]]` contains the centered binary outcome
  \\y\_{ic} - \bar{y}\_c\\ (or the covariate residual \\y\_{ic} -
  \hat{\mu}\_{ic}\\) for subject \\i\\; entries in other columns are
  ignored.

- f:

  Numeric scalar or length-\\C\\ vector. Minor allele frequency (MAF).
  Scalar for common MAF; length-\\C\\ vector for study-specific MAFs
  (multi-ancestry case-control). Must be in \\(0, 0.5\]\\.

- b:

  Numeric scalar or vector of threshold values. When a vector, \\\xi_A\\
  is solved once at the anchor \\b_0 = \mathrm{median}(\mathtt{b})\\ (or
  `b_anchor`).

- study_id:

  Integer vector of length \\N\\ with values in \\\\1,\ldots,C\\\\.
  Activates case-control mode. Default `NULL` (eQTL mode).

- g_center:

  Numeric vector of length \\N\\. Per-subject genotype centering values
  \\\tilde g_i\\. Default `NULL` uses \\2f\_{c_i}\\.

- tau:

  Numeric positive vector of length \\C\\. Per-cell-type scale factors
  for covariate-adjusted eQTL weights. Default `NULL`.

- scale_Y:

  Logical. If `TRUE` (default), `Y` is centered and scaled to mean zero
  and unit variance per column. Set to `FALSE` only if `Y` has already
  been pre-processed.

- two.sided:

  Logical. Default `TRUE`.

- K:

  Positive integer. Number of IS simulations. Default `5000L`.

- b_anchor:

  Numeric scalar. The anchor threshold \\b_0\\ at which \\\xi\_{A,\pm}\\
  are solved. Defaults to `NULL`, using \\b_0 = \mathrm{median}(b)\\.
  Ignored for scalar `b`.

- newton_tol:

  Positive numeric. Newton convergence tolerance. Default `1e-10`.

- newton_max_iter:

  Positive integer. Maximum Newton iterations. Default `100L`.

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

  Lower bound of 95\\ \\\hat{p} - 1.96\\\widehat{\mathrm{SE}}\\,
  truncated to \\\[0,1\]\\.

- `p.upper`:

  Upper bound of 95\\ \\\[0,1\]\\.

- `b_anchor`:

  Anchor threshold used for Newton solve.

- `K`:

  Number of IS simulations used.

## Details

Conditional IS for ASSET When Normality Fails: Independent Cell Types or
Studies

**Two-sided test.** A proper two-sided mixture tilt is used rather than
the heuristic of halving the threshold and doubling the estimate.
Because the CGF \\\phi_A(\xi)\\ is not symmetric in \\\xi\\ for general
(data-dependent) weights \\\omega\_{i,A}\\, two per-subset tilting
parameters are solved: \\\xi\_{A,+}\\ satisfying \\\phi_A'(\xi\_{A,+}) =
+b_0\\ (tilting toward the upper tail) and \\\xi\_{A,-}\\ satisfying
\\\phi_A'(\xi\_{A,-}) = -b_0\\ (tilting toward the lower tail). For each
simulation, a driving subset \\A_k\\ and a sign \\s_k \in \\+1,-1\\\\
are chosen uniformly at random. Genotypes are drawn under the tilted
distribution corresponding to \\\xi\_{A_k,s_k}\\, and the IS weight
denominator involves \\2(2^C-1)\\ terms: \$\$\frac{2(2^C-1)}{\sum_A
\left\[ \exp\\\xi\_{A,+} Z_A - \phi_A(\xi\_{A,+})\\ + \exp\\\xi\_{A,-}
Z_A - \phi_A(\xi\_{A,-})\\ \right\]} \mathcal{I}\\\left\\\max_A \|Z_A\|
\> b\right\\.\$\$

**Range of thresholds.** When `b` is a vector, \\\xi\_{A,\pm}\\ are
solved once at the anchor \\b_0 = \mathrm{median}(b)\\ (or `b_anchor`)
at the original scale of \\b\\.

**Scaling of \\Y\\.** Each column of \\Y\\ should have mean zero and
unit variance. If `scale_Y = TRUE` (default), the function applies
[`base::scale`](https://rdrr.io/r/base/scale.html) to \\Y\\ internally.

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

[`asset_is_nonnorm_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_corr.md),
[`asset_is_pvalue_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_ind.md),
[`asset_is_pvalue_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_corr.md)

## Examples

``` r
set.seed(42)
N <- 200; C <- 5
f <- 0.2

## --- eQTL mode ---
Y <- matrix(rnorm(N * C), N, C)
b_obs <- 4.5
asset_is_nonnorm_ind(Y = Y, f = f, b = b_obs, K = 2000L, two.sided = TRUE)
#>     b        p.est         p.se      p.lower      p.upper b_anchor    K
#> 1 4.5 0.0001317644 7.299017e-06 0.0001174583 0.0001460705      4.5 2000

## Range of thresholds (range approximation)
b_grid <- seq(3.5, 6.0, by = 0.1)
asset_is_nonnorm_ind(Y = Y, f = f, b = b_grid, K = 1e4, two.sided = TRUE)
#>      b        p.est         p.se      p.lower      p.upper b_anchor     K
#> 1  3.5 6.848820e-03 2.241423e-04 6.409501e-03 7.288139e-03     4.75 10000
#> 2  3.6 4.809340e-03 1.482995e-04 4.518673e-03 5.100007e-03     4.75 10000
#> 3  3.7 3.448998e-03 1.004270e-04 3.252161e-03 3.645835e-03     4.75 10000
#> 4  3.8 2.417772e-03 6.898457e-05 2.282562e-03 2.552982e-03     4.75 10000
#> 5  3.9 1.594580e-03 4.316903e-05 1.509969e-03 1.679192e-03     4.75 10000
#> 6  4.0 1.115857e-03 2.977806e-05 1.057492e-03 1.174222e-03     4.75 10000
#> 7  4.1 7.463551e-04 1.970268e-05 7.077378e-04 7.849723e-04     4.75 10000
#> 8  4.2 4.823754e-04 1.222078e-05 4.584227e-04 5.063282e-04     4.75 10000
#> 9  4.3 3.201975e-04 8.098006e-06 3.043254e-04 3.360696e-04     4.75 10000
#> 10 4.4 2.023438e-04 5.148871e-06 1.922521e-04 2.124356e-04     4.75 10000
#> 11 4.5 1.273127e-04 3.262404e-06 1.209183e-04 1.337070e-04     4.75 10000
#> 12 4.6 8.329940e-05 2.160534e-06 7.906475e-05 8.753404e-05     4.75 10000
#> 13 4.7 5.029585e-05 1.301642e-06 4.774463e-05 5.284707e-05     4.75 10000
#> 14 4.8 3.150545e-05 8.341452e-07 2.987053e-05 3.314038e-05     4.75 10000
#> 15 4.9 1.946267e-05 5.395319e-07 1.840519e-05 2.052016e-05     4.75 10000
#> 16 5.0 1.175900e-05 3.363267e-07 1.109980e-05 1.241820e-05     4.75 10000
#> 17 5.1 6.934379e-06 2.036932e-07 6.535140e-06 7.333617e-06     4.75 10000
#> 18 5.2 4.077761e-06 1.247328e-07 3.833284e-06 4.322237e-06     4.75 10000
#> 19 5.3 2.339265e-06 7.639016e-08 2.189540e-06 2.488990e-06     4.75 10000
#> 20 5.4 1.365974e-06 4.753852e-08 1.272799e-06 1.459150e-06     4.75 10000
#> 21 5.5 8.127734e-07 2.986985e-08 7.542284e-07 8.713183e-07     4.75 10000
#> 22 5.6 4.757698e-07 1.856078e-08 4.393906e-07 5.121489e-07     4.75 10000
#> 23 5.7 2.703817e-07 1.124850e-08 2.483347e-07 2.924288e-07     4.75 10000
#> 24 5.8 1.525255e-07 6.798568e-09 1.392003e-07 1.658507e-07     4.75 10000
#> 25 5.9 8.460472e-08 4.165577e-09 7.644019e-08 9.276925e-08     4.75 10000
#> 26 6.0 4.743095e-08 2.680086e-09 4.217798e-08 5.268392e-08     4.75 10000

## --- eQTL with covariate adjustment ---
## Y contains expression residuals y'_ic after regressing out covariates
tau_c <- apply(scale(Y), 2, sd) * sqrt(2 * f * (1 - f))
g_center <- rep(2 * f, N)   # replace with predicted g_i0 in real applications
asset_is_nonnorm_ind(Y = Y, f = f, b = b_obs, tau = tau_c,
                 g_center = g_center, K = 2000L)
#>     b      p.est        p.se    p.lower    p.upper b_anchor    K
#> 1 4.5 0.08311106 0.003436058 0.07637638 0.08984573      4.5 2000

## --- Case-control meta-analysis (same ancestry) ---
study_id <- sample(1:C, N, replace = TRUE)
phi_c    <- runif(C, 0.2, 0.4)   # case fractions
y_raw    <- rbinom(N, 1, prob = phi_c[study_id])
Ycc      <- matrix(0, N, C)
for (i in seq_len(N)) {
  c_i <- study_id[i]
  Ycc[i, c_i] <- y_raw[i] - phi_c[c_i]
}
asset_is_nonnorm_ind(Y = Ycc, f = f, b = b_obs, study_id = study_id, K = 2000L)
#>     b      p.est        p.se    p.lower    p.upper b_anchor    K
#> 1 4.5 0.04954485 0.002294865 0.04504692 0.05404279      4.5 2000

## --- Case-control, different ancestries (study-specific MAF) ---
f_c <- runif(C, 0.1, 0.4)
asset_is_nonnorm_ind(Y = Ycc, f = f_c, b = b_obs, study_id = study_id, K = 2000L)
#>     b      p.est     p.se    p.lower    p.upper b_anchor    K
#> 1 4.5 0.07372001 0.003624 0.06661697 0.08082305      4.5 2000
```
