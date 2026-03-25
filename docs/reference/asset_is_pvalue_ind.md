# Gaussian IS Estimator of ASSET Tail Probabilities for Independent Studies

Estimates \\p(b) = P_0(\max_A Z_A \> b)\\ (one-sided) or \\p(b) =
P_0(\max_A \|Z_A\| \> b)\\ (two-sided) for the ASSET statistic when the
\\M\\ study-level score statistics \\(z_1,\ldots,z_M)\\ are independent
and approximately normally distributed under the global null. The IS
algorithm uses exponential tilting of the joint Gaussian distribution
toward the target region, substantially reducing variance compared with
naive Monte Carlo for extreme thresholds.

For each non-empty subset \\A \subseteq \\1,\ldots,M\\\\, the
fixed-effect meta-analysis statistic is \$\$Z_A = \frac{\sum\_{m \in A}
\sqrt{n_m}\\ z_m}{\sqrt{\sum\_{l \in A} n_l}}.\$\$ The IS estimator is
based on the likelihood-ratio identity \$\$P_0(\max_A Z_A \> b) =
\frac{1}{2^M-1} \sum\_{A_0} E\_{\xi,A_0}\\\left( \frac{2^M-1}{\sum_A
\exp\\\xi Z_A - \xi^2/2\\} \\\mathcal{I}\\\max_A Z_A \> b\\ \right),\$\$
with the near-optimal choice \\\xi = \mathrm{median}(\mathtt{b})\\.

## Usage

``` r
asset_is_pvalue_ind(
  z,
  n,
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

  Numeric vector of length \\M\\. Observed study-level score statistics,
  approximately \\N(0,1)\\ under the global null.

- n:

  Numeric vector of length \\M\\. Study sample sizes. For case-control
  studies supply effective sample sizes \\n_m^{\mathrm{eff}} = n_m
  \phi_m (1-\phi_m)\\.

- b:

  Numeric scalar or vector of threshold values. When a vector, a single
  \\\xi = \mathrm{median}(\mathtt{b})\\ is reused across all thresholds.

- K:

  Positive integer. Number of IS simulations. Default `5000L`.

- two.sided:

  Logical. If `TRUE` (default) compute \\P_0(\max_A \|Z_A\| \> b)\\; if
  `FALSE` compute \\P_0(\max_A Z_A \> b)\\.

- xi:

  Numeric scalar. Tilting parameter. Defaults to `median(b)` (after
  two-sided halving). Override to supply a custom value.

- parallel:

  Logical. If `TRUE` distribute simulations across `n.cores` workers.
  Default `FALSE`.

- n.cores:

  Positive integer. Workers used when `parallel = TRUE`. Default
  `parallel::detectCores() - 1`.

- seed:

  Integer or `NULL`. Master random seed. Default `NULL`.

## Value

A data frame with one row per threshold value in `b` and the following
columns:

- `b`:

  The threshold value.

- `p.est`:

  IS point estimate of \\p(b)\\.

- `p.se`:

  Estimated standard error of `p.est`.

- `p.lower`:

  Lower bound of a 95\\ \\\hat{p} - 1.96\\\widehat{\mathrm{SE}}\\,
  truncated to \\\[0,1\]\\.

- `p.upper`:

  Upper bound of a 95\\ truncated to \\\[0,1\]\\.

- `xi`:

  Tilting parameter used for the corresponding threshold.

- `K`:

  Number of IS simulations used.

## Details

Importance Sampling for ASSET: Independent Studies, Gaussian
Approximation

**Case-control studies.** Replace each \\n_m\\ with the effective sample
size \\n_m^{\mathrm{eff}} = n_m \phi_m (1-\phi_m)\\, where \\\phi_m\\ is
the fraction of cases in study \\m\\.

**Range of thresholds.** When `b` is a vector, a single tilting
parameter \\\xi = \mathrm{median}(\mathtt{b})\\ is shared, and all
thresholds are evaluated in one pass via `findInterval`.

**Two-sided test.** Uses the symmetry \\P_0(\max_A \|Z_A\| \> b) = 2
P_0(\max_A Z_A \> b/2)\\.

**Parallel computation.** When `parallel = TRUE` the simulations are
distributed across `n.cores` workers using `foreach` with a `doParallel`
backend.

For correlated studies use
[`asset_is_pvalue_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_corr.md).
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

[`asset_is_pvalue_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_corr.md),
[`asset_is_nonnorm_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_ind.md),
[`asset_is_nonnorm_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_corr.md)

## Examples

``` r
set.seed(42)
M <- 5
z <- rnorm(M, mean = c(3, 3, 0, 0, 0))   # two studies have signal
n <- c(500, 600, 450, 700, 550)

## Two-sided p-value at a single threshold
b_obs <- max(abs(z))
asset_is_pvalue_ind(z = z, n = n, b = b_obs, K = 1e4, two.sided = TRUE)
#>          b        p.est         p.se      p.lower      p.upper       xi     K
#> 1 4.370958 0.0002649982 5.833978e-06 0.0002535636 0.0002764328 4.370958 10000

## One-sided p-values over a grid (range approximation)
b_grid <- seq(3, 6, by = 0.1)
asset_is_pvalue_ind(z = z, n = n, b = b_grid, K = 5e4, two.sided = FALSE)
#>      b        p.est         p.se      p.lower      p.upper  xi     K
#> 1  3.0 1.812728e-02 2.906282e-04 1.755765e-02 1.869691e-02 4.5 50000
#> 2  3.1 1.331330e-02 1.986254e-04 1.292399e-02 1.370260e-02 4.5 50000
#> 3  3.2 9.955434e-03 1.397337e-04 9.681556e-03 1.022931e-02 4.5 50000
#> 4  3.3 7.354234e-03 9.682617e-05 7.164455e-03 7.544013e-03 4.5 50000
#> 5  3.4 5.309203e-03 6.611327e-05 5.179621e-03 5.438785e-03 4.5 50000
#> 6  3.5 3.829000e-03 4.597845e-05 3.738882e-03 3.919117e-03 4.5 50000
#> 7  3.6 2.662120e-03 3.061344e-05 2.602117e-03 2.722122e-03 4.5 50000
#> 8  3.7 1.863929e-03 2.068715e-05 1.823382e-03 1.904475e-03 4.5 50000
#> 9  3.8 1.290747e-03 1.388426e-05 1.263534e-03 1.317960e-03 4.5 50000
#> 10 3.9 8.792855e-04 9.208833e-06 8.612362e-04 8.973348e-04 4.5 50000
#> 11 4.0 5.988095e-04 6.167188e-06 5.867218e-04 6.108971e-04 4.5 50000
#> 12 4.1 4.029991e-04 4.078569e-06 3.950051e-04 4.109931e-04 4.5 50000
#> 13 4.2 2.661005e-04 2.665589e-06 2.608760e-04 2.713251e-04 4.5 50000
#> 14 4.3 1.750540e-04 1.756997e-06 1.716103e-04 1.784977e-04 4.5 50000
#> 15 4.4 1.126198e-04 1.129806e-06 1.104054e-04 1.148342e-04 4.5 50000
#> 16 4.5 7.250396e-05 7.388301e-07 7.105585e-05 7.395206e-05 4.5 50000
#> 17 4.6 4.610281e-05 4.733966e-07 4.517495e-05 4.703066e-05 4.5 50000
#> 18 4.7 2.865625e-05 3.013733e-07 2.806556e-05 2.924694e-05 4.5 50000
#> 19 4.8 1.799323e-05 1.937781e-07 1.761343e-05 1.837304e-05 4.5 50000
#> 20 4.9 1.118614e-05 1.247378e-07 1.094165e-05 1.143062e-05 4.5 50000
#> 21 5.0 6.714722e-06 7.791024e-08 6.562018e-06 6.867426e-06 4.5 50000
#> 22 5.1 4.033499e-06 4.884330e-08 3.937766e-06 4.129232e-06 4.5 50000
#> 23 5.2 2.406219e-06 3.056879e-08 2.346304e-06 2.466134e-06 4.5 50000
#> 24 5.3 1.435335e-06 1.928637e-08 1.397534e-06 1.473136e-06 4.5 50000
#> 25 5.4 8.250405e-07 1.176461e-08 8.019819e-07 8.480991e-07 4.5 50000
#> 26 5.5 4.774874e-07 7.243902e-09 4.632894e-07 4.916855e-07 4.5 50000
#> 27 5.6 2.730047e-07 4.417578e-09 2.643463e-07 2.816632e-07 4.5 50000
#> 28 5.7 1.588988e-07 2.750385e-09 1.535080e-07 1.642895e-07 4.5 50000
#> 29 5.8 8.871214e-08 1.669888e-09 8.543916e-08 9.198512e-08 4.5 50000
#> 30 5.9 4.855330e-08 9.979491e-10 4.659732e-08 5.050928e-08 4.5 50000
#> 31 6.0 2.566336e-08 5.856667e-10 2.451545e-08 2.681126e-08 4.5 50000

## Case-control: supply effective sample sizes
phi <- c(0.3, 0.4, 0.5, 0.35, 0.45)     # case fractions
n_eff <- n * phi * (1 - phi)
asset_is_pvalue_ind(z = z, n = n_eff, b = b_obs, K = 1e4, two.sided = TRUE)
#>          b        p.est        p.se      p.lower      p.upper       xi     K
#> 1 4.370958 0.0002595515 5.79281e-06 0.0002481976 0.0002709054 4.370958 10000

## Parallel execution
if (FALSE) { # \dontrun{
asset_is_pvalue_ind(z = z, n = n, b = b_grid, K = 1e5,
                parallel = TRUE, n.cores = 4)
} # }
```
