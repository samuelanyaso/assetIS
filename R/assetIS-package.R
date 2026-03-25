#' assetIS: Importance Sampling for ASSET Tail Probability Estimation
#'
#' @description
#' The \pkg{assetIS} package provides efficient importance sampling (IS)
#' estimators of extreme tail probabilities for the ASSET subset-based
#' meta-analysis statistic (Bhattacharjee et al., 2012).
#'
#' @details
#' The ASSET statistic is defined as
#' \deqn{Z_{\mathrm{ASSET}} = \max_{A \neq \varnothing} |Z_A|,}
#' where the maximum is taken over all \eqn{2^M - 1} non-empty subsets \eqn{A}
#' of \eqn{M} studies (or cell types), and \eqn{Z_A} is a fixed-effect
#' meta-analysis statistic for the studies in \eqn{A}. Evaluating the
#' significance of \eqn{Z_{\mathrm{ASSET}}} requires estimating an extreme
#' tail probability, which is computationally prohibitive by naive Monte
#' Carlo. The package implements exponential-tilting IS estimators that
#' achieve large variance reductions by concentrating simulation effort in
#' the target region.
#'
#' Four estimators are provided, covering two design dimensions:
#'
#' \describe{
#'   \item{\strong{Gaussian vs. non-Gaussian:}}{
#'     Gaussian estimators ([asset_is_pvalue_ind()], [asset_is_pvalue_corr()])
#'     simulate study-level \eqn{z}-statistics from a (possibly correlated)
#'     normal distribution; they are appropriate when sample sizes are large
#'     enough to justify asymptotic normality.  Non-Gaussian estimators
#'     ([asset_is_nonnorm_ind()], [asset_is_nonnorm_corr()]) simulate
#'     genotypes directly from a trinomial (Hardy--Weinberg) distribution and
#'     condition on the observed phenotype data \eqn{Y}; they are designed
#'     for small-sample settings such as single-cell eQTL or rare-disease
#'     meta-analysis.}
#'   \item{\strong{Independent vs. correlated:}}{
#'     Independent estimators ([asset_is_pvalue_ind()],
#'     [asset_is_nonnorm_ind()]) assume score statistics are uncorrelated
#'     across studies / cell types.  Correlated estimators
#'     ([asset_is_pvalue_corr()], [asset_is_nonnorm_corr()]) account for
#'     inter-study correlation via a known or estimated correlation matrix.}
#' }
#'
#' All four functions share a common interface and support:
#' \itemize{
#'   \item One-sided and two-sided p-values.
#'   \item Scalar or vector thresholds \eqn{b} (range approximation for
#'     vectors, reusing a single tilting parameter over a neighbourhood).
#'   \item Case-control studies (effective sample sizes or binary outcome
#'     matrices).
#'   \item Multi-ancestry analyses with study-specific allele frequencies.
#'   \item Covariate adjustment (non-Gaussian estimators).
#'   \item Optional parallel computation via \pkg{doParallel}.
#'   \item Reproducibility via a master \code{seed} argument.
#' }
#'
#' @section Choosing an estimator:
#' \tabular{lll}{
#'   \strong{Setting}          \tab \strong{Independent}          \tab \strong{Correlated} \cr
#'   Large-sample (Gaussian)   \tab [asset_is_pvalue_ind()]        \tab [asset_is_pvalue_corr()] \cr
#'   Small-sample (non-Gaussian) \tab [asset_is_nonnorm_ind()]     \tab [asset_is_nonnorm_corr()] \cr
#' }
#'
#' @references
#' Bhattacharjee S, Rajaraman P, Bhattacharjee S, et al. (2012). A
#' subset-based approach improves power and interpretation for the combined
#' analysis of genetic association studies of heterogeneous traits.
#' \emph{American Journal of Human Genetics}, \bold{90}(5), 821--835.
#'
#' Shi J and Siegmund D (2007). Importance sampling for estimating p-values
#' in linkage analysis. \emph{Biostatistics}, \bold{8}(3), 587--600.
#'
#' Siegmund D (1976). Importance sampling in the Monte Carlo study of
#' sequential tests. \emph{The Annals of Statistics}, \bold{4}(4), 673--684.
#'
#' @keywords internal
"_PACKAGE"
