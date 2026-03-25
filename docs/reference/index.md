# Package index

## Gaussian settings (large samples)

Use these when study-level score statistics are approximately normal.
Appropriate for large-sample continuous-trait or case-control
meta-analysis.

- [`asset_is_pvalue_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_ind.md)
  : Gaussian IS Estimator of ASSET Tail Probabilities for Independent
  Studies
- [`asset_is_pvalue_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_corr.md)
  : Gaussian IS Estimator of ASSET Tail Probabilities for Correlated or
  Overlapping Studies

## Non-Gaussian settings (small samples)

Use these when normality cannot be assumed — single-cell eQTL,
rare-disease meta-analysis, or any setting with small to moderate sample
sizes. Genotypes are simulated from a discrete Hardy–Weinberg
distribution conditional on the observed phenotype matrix.

- [`asset_is_nonnorm_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_ind.md)
  : Conditional IS Estimator of ASSET Tail Probabilities for
  Non-Gaussian, Independent Settings
- [`asset_is_nonnorm_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_corr.md)
  : Conditional IS Estimator of ASSET Tail Probabilities for
  Non-Gaussian Settings with Correlated Score Statistics

## Package

- [`assetIS`](https://anyasosamuelcs.github.io/assetIS/reference/assetIS-package.md)
  [`assetIS-package`](https://anyasosamuelcs.github.io/assetIS/reference/assetIS-package.md)
  : assetIS: Importance Sampling for ASSET Tail Probability Estimation
