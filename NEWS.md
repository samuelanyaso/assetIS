# assetIS 0.1.0

## Initial release

### New functions

* `asset_is_pvalue_ind()`: Gaussian IS for the ASSET statistic assuming
  independent studies and asymptotic normality of score statistics.

* `asset_is_pvalue_corr()`: Gaussian IS for the ASSET statistic when
  study-level score statistics are correlated (e.g., due to overlapping
  subjects).

* `asset_is_nonnorm_ind()`: Conditional IS for the ASSET statistic in
  non-Gaussian settings (single-cell eQTL, rare-disease meta-analysis)
  assuming independent cell types or studies. Genotypes are simulated
  from a discrete Hardy--Weinberg distribution conditional on the
  observed phenotype matrix. Supports eQTL, covariate-adjusted eQTL,
  and case-control modes, with common or study-specific allele
  frequencies.

* `asset_is_nonnorm_corr()`: Conditional IS for the ASSET statistic in
  non-Gaussian settings with correlated score statistics. GLS weights
  derived from the empirical conditional covariance
  $\widehat\Sigma_Y = Y^\top Y / N$ replace the equal-weight formula.

### Features shared across all four functions

* One-sided and two-sided p-values.
* Scalar or vector thresholds with the range approximation (single
  tilting parameter reused over a neighbourhood of thresholds).
* Case-control support via effective sample sizes (Gaussian estimators)
  or binary outcome matrices (non-Gaussian estimators).
* Multi-ancestry analyses with study-specific allele frequencies
  (non-Gaussian estimators).
* Covariate adjustment (non-Gaussian estimators).
* Optional parallel computation via **doParallel**.
* Reproducibility via a master `seed` argument with per-worker derived
  seeds.
* 95% Wald confidence intervals for all p-value estimates.
