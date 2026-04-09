# assetIS ![](reference/figures/logo.png)

## Overview

**assetIS** provides efficient importance sampling (IS) estimators of
extreme tail probabilities for the
[ASSET](https://doi.org/10.1016/j.ajhg.2012.03.015) subset-based
meta-analysis statistic (Bhattacharjee et al., 2012).

The ASSET statistic is defined as

$$Z_{\text{ASSET}} = \max\limits_{A \neq \varnothing}\left| Z_{A} \right|,$$

where the maximum is over all $2^{M} - 1$ non-empty subsets $A$ of $M$
studies (or cell types). Evaluating its significance requires an extreme
tail probability, which is computationally intractable by naive Monte
Carlo. The IS estimators in **assetIS** use exponential tilting to
concentrate simulation effort in the target region, achieving
orders-of-magnitude variance reductions.

## Four estimators — one unified interface

| Setting                         | Independent                                                                                            | Correlated                                                                                               |
|---------------------------------|--------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
| **Large-sample** (Gaussian)     | [`asset_is_pvalue_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_ind.md)   | [`asset_is_pvalue_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_corr.md)   |
| **Small-sample** (Non-Gaussian) | [`asset_is_nonnorm_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_ind.md) | [`asset_is_nonnorm_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_corr.md) |

All four functions share the same design:

- **One-sided or two-sided** p-values
- **Scalar or vector thresholds** (efficient range approximation)
- **Case-control studies** (effective sample sizes or binary outcome
  matrices)
- **Multi-ancestry** analyses with study-specific allele frequencies
- **Covariate adjustment** (non-Gaussian estimators)
- **Parallel computation** via `doParallel`
- **Reproducibility** via a master `seed` argument

## When to use which estimator

Use the **Gaussian** estimators (`asset_is_pvalue_*`) when sample sizes
are large enough to justify asymptotic normality of the study-level
score statistics.

Use the **non-Gaussian** estimators (`asset_is_nonnorm_*`) when: -
Sample sizes are small (e.g., single-cell eQTL across cell types) - The
trait is rare (rare-disease meta-analysis) - You wish to condition on
the observed phenotype data rather than rely on a Gaussian approximation

Use the **correlated** variants (`*_corr`) when: - Expression
measurements are correlated across cell types (eQTL) - Studies share
subjects (overlapping case-control cohorts)

## Installation

Install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("samuelanyaso/assetIS")
```

## Quick start

``` r
library(assetIS)

## --- Gaussian IS, independent studies ---
set.seed(1)
M <- 5
n <- c(500, 600, 450, 700, 550)
z <- c(3.1, 2.9, 0.2, -0.1, 0.4)

asset_is_pvalue_ind(z = z, n = n, b = max(abs(z)), K = 5000L)

## --- Gaussian IS, correlated studies ---
rho   <- 0.35
Sigma <- outer(1:M, 1:M, function(i, j) rho^abs(i - j))
asset_is_pvalue_corr(z = z, n = n, Sigma = Sigma, b = max(abs(z)), K = 5000L)

## --- Non-Gaussian IS, single-cell eQTL ---
N <- 100; C <- 5; f <- 0.15
Y <- scale(matrix(rnorm(N * C), N, C))
asset_is_nonnorm_ind(Y = Y, f = f, b = 4.5, K = 2000L)

## --- Non-Gaussian IS, correlated cell types ---
asset_is_nonnorm_corr(Y = Y, f = f, b = 4.5, K = 2000L)
```

## Case-control and multi-ancestry examples

``` r
## Gaussian IS with effective sample sizes (case-control)
phi   <- c(0.30, 0.40, 0.50, 0.35, 0.45)
n_eff <- n * phi * (1 - phi)
asset_is_pvalue_ind(z = z, n = n_eff, b = max(abs(z)), K = 5000L)

## Non-Gaussian IS, multi-ancestry (study-specific MAF)
f_multi  <- c(0.15, 0.22, 0.10, 0.30, 0.18)
study_id <- sample(1:C, N, replace = TRUE)
Ycc      <- matrix(0, N, C)
# ... fill Ycc with centered binary outcomes ...
asset_is_nonnorm_ind(Y = Ycc, f = f_multi, b = 4.5,
                     study_id = study_id, K = 2000L)
```

## Parallel computation

``` r
asset_is_pvalue_ind(z = z, n = n, b = seq(3, 6, by = 0.1),
                    K = 1e5L, parallel = TRUE, n.cores = 4L, seed = 42L)
```

## Citation

If you use **assetIS** in published work, please cite:

> Anyaso-Samuel S (2025). *assetIS: Importance Sampling for ASSET Tail
> Probability Estimation*. R package version 0.1.0.
> <https://github.com/samuelanyaso/assetIS>

## References

- Bhattacharjee S et al. (2012). *Am J Hum Genet*, **90**(5), 821–835.
- Shi J et al. (2007). *Journal of the American Statistical
  Association*, **102**(479), 929–937.
- Siegmund D (1976). *Ann Statist*, **4**(4), 673–684.

## License

GPL (\>= 3)
