# assetIS: Importance Sampling for ASSET Tail Probability Estimation

The assetIS package provides efficient importance sampling (IS)
estimators of extreme tail probabilities for the ASSET subset-based
meta-analysis statistic (Bhattacharjee et al., 2012).

## Details

The ASSET statistic is defined as \$\$Z\_{\mathrm{ASSET}} = \max\_{A
\neq \varnothing} \|Z_A\|,\$\$ where the maximum is taken over all
\\2^M - 1\\ non-empty subsets \\A\\ of \\M\\ studies (or cell types),
and \\Z_A\\ is a fixed-effect meta-analysis statistic for the studies in
\\A\\. Evaluating the significance of \\Z\_{\mathrm{ASSET}}\\ requires
estimating an extreme tail probability, which is computationally
prohibitive by naive Monte Carlo. The package implements
exponential-tilting IS estimators that achieve large variance reductions
by concentrating simulation effort in the target region.

Four estimators are provided, covering two design dimensions:

- **Gaussian vs. non-Gaussian:**:

  Gaussian estimators
  ([`asset_is_pvalue_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_ind.md),
  [`asset_is_pvalue_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_corr.md))
  simulate study-level \\z\\-statistics from a (possibly correlated)
  normal distribution; they are appropriate when sample sizes are large
  enough to justify asymptotic normality. Non-Gaussian estimators
  ([`asset_is_nonnorm_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_ind.md),
  [`asset_is_nonnorm_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_corr.md))
  simulate genotypes directly from a trinomial (Hardy–Weinberg)
  distribution and condition on the observed phenotype data \\Y\\; they
  are designed for small-sample settings such as single-cell eQTL or
  rare-disease meta-analysis.

- **Independent vs. correlated:**:

  Independent estimators
  ([`asset_is_pvalue_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_ind.md),
  [`asset_is_nonnorm_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_ind.md))
  assume score statistics are uncorrelated across studies / cell types.
  Correlated estimators
  ([`asset_is_pvalue_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_corr.md),
  [`asset_is_nonnorm_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_corr.md))
  account for inter-study correlation via a known or estimated
  correlation matrix.

All four functions share a common interface and support:

- One-sided and two-sided p-values.

- Scalar or vector thresholds \\b\\ (range approximation for vectors,
  reusing a single tilting parameter over a neighbourhood).

- Case-control studies (effective sample sizes or binary outcome
  matrices).

- Multi-ancestry analyses with study-specific allele frequencies.

- Covariate adjustment (non-Gaussian estimators).

- Optional parallel computation via doParallel.

- Reproducibility via a master `seed` argument.

## Choosing an estimator

|  |  |  |
|----|----|----|
| **Setting** | **Independent** | **Correlated** |
| Large-sample (Gaussian) | [`asset_is_pvalue_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_ind.md) | [`asset_is_pvalue_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_pvalue_corr.md) |
| Small-sample (non-Gaussian) | [`asset_is_nonnorm_ind()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_ind.md) | [`asset_is_nonnorm_corr()`](https://anyasosamuelcs.github.io/assetIS/reference/asset_is_nonnorm_corr.md) |

## References

Bhattacharjee S, Rajaraman P, Bhattacharjee S, et al. (2012). A
subset-based approach improves power and interpretation for the combined
analysis of genetic association studies of heterogeneous traits.
*American Journal of Human Genetics*, **90**(5), 821–835.

Shi J and Siegmund D (2007). Importance sampling for estimating p-values
in linkage analysis. *Biostatistics*, **8**(3), 587–600.

Siegmund D (1976). Importance sampling in the Monte Carlo study of
sequential tests. *The Annals of Statistics*, **4**(4), 673–684.

## See also

Useful links:

- <https://github.com/samuelanyaso/assetIS>

- Report bugs at <https://github.com/samuelanyaso/assetIS/issues>

## Author

**Maintainer**: Samuel Anyaso-Samuel <samuel.anyaso-samuel@nih.gov>
