#' Importance Sampling for ASSET: Independent Studies, Gaussian Approximation
#'
#' @title Gaussian IS Estimator of ASSET Tail Probabilities for
#'   Independent Studies
#'
#' @description
#' Estimates \eqn{p(b) = P_0(\max_A Z_A > b)} (one-sided) or
#' \eqn{p(b) = P_0(\max_A |Z_A| > b)} (two-sided) for the ASSET statistic
#' when the \eqn{M} study-level score statistics \eqn{(z_1,\ldots,z_M)} are
#' independent and approximately normally distributed under the global null.
#' The IS algorithm uses exponential tilting of the joint Gaussian distribution
#' toward the target region, substantially reducing variance compared with
#' naive Monte Carlo for extreme thresholds.
#'
#' For each non-empty subset \eqn{A \subseteq \{1,\ldots,M\}}, the
#' fixed-effect meta-analysis statistic is
#' \deqn{Z_A = \frac{\sum_{m \in A} \sqrt{n_m}\, z_m}{\sqrt{\sum_{l \in A} n_l}}.}
#' The IS estimator is based on the likelihood-ratio identity
#' \deqn{P_0(\max_A Z_A > b) = \frac{1}{2^M-1} \sum_{A_0}
#'   E_{\xi,A_0}\!\left(
#'     \frac{2^M-1}{\sum_A \exp\{\xi Z_A - \xi^2/2\}}
#'     \,\mathcal{I}\{\max_A Z_A > b\}
#'   \right),}
#' with the near-optimal choice \eqn{\xi = \mathrm{median}(\mathtt{b})}.
#'
#' @details
#' \strong{Case-control studies.} Replace each \eqn{n_m} with the effective
#' sample size \eqn{n_m^{\mathrm{eff}} = n_m \phi_m (1-\phi_m)}, where
#' \eqn{\phi_m} is the fraction of cases in study \eqn{m}.
#'
#' \strong{Range of thresholds.} When \code{b} is a vector, a single
#' tilting parameter \eqn{\xi = \mathrm{median}(\mathtt{b})} is shared, and
#' all thresholds are evaluated in one pass via \code{findInterval}.
#'
#' \strong{Two-sided test.} Uses the symmetry
#' \eqn{P_0(\max_A |Z_A| > b) = 2 P_0(\max_A Z_A > b/2)}.
#'
#'
#' \strong{Parallel computation.}  When \code{parallel = TRUE} the
#' simulations are distributed across \code{n.cores} workers using
#' \code{foreach} with a \code{doParallel} backend. 
#'
#' For correlated studies use [asset_is_pvalue_corr()]. For non-Gaussian
#' settings use [asset_is_nonnorm_ind()] or [asset_is_nonnorm_corr()].
#'
#' @param z Numeric vector of length \eqn{M}. Observed study-level score
#'   statistics, approximately \eqn{N(0,1)} under the global null.
#' @param n Numeric vector of length \eqn{M}. Study sample sizes. For
#'   case-control studies supply effective sample sizes
#'   \eqn{n_m^{\mathrm{eff}} = n_m \phi_m (1-\phi_m)}.
#' @param b Numeric scalar or vector of threshold values. When a vector, a
#'   single \eqn{\xi = \mathrm{median}(\mathtt{b})} is reused across all
#'   thresholds.
#' @param K Positive integer. Number of IS simulations. Default \code{5000L}.
#' @param two.sided Logical. If \code{TRUE} (default) compute
#'   \eqn{P_0(\max_A |Z_A| > b)}; if \code{FALSE} compute
#'   \eqn{P_0(\max_A Z_A > b)}.
#' @param xi Numeric scalar. Tilting parameter. Defaults to
#'   \code{median(b)} (after two-sided halving). Override to supply a
#'   custom value.
#' @param parallel Logical. If \code{TRUE} distribute simulations across
#'   \code{n.cores} workers. Default \code{FALSE}.
#' @param n.cores Positive integer. Workers used when
#'   \code{parallel = TRUE}. Default \code{parallel::detectCores() - 1}.
#' @param seed Integer or \code{NULL}. Master random seed. Default
#'   \code{NULL}.
#'
#' @return A data frame with one row per threshold value in \code{b} and the
#'   following columns:
#'   \describe{
#'     \item{\code{b}}{The threshold value.}
#'     \item{\code{p.est}}{IS point estimate of \eqn{p(b)}.}
#'     \item{\code{p.se}}{Estimated standard error of \code{p.est}.}
#'     \item{\code{p.lower}}{Lower bound of a 95\% Wald confidence interval,
#'       \eqn{\hat{p} - 1.96\,\widehat{\mathrm{SE}}}, truncated to \eqn{[0,1]}.}
#'     \item{\code{p.upper}}{Upper bound of a 95\% Wald confidence interval,
#'       truncated to \eqn{[0,1]}.}
#'     \item{\code{xi}}{Tilting parameter used for the corresponding threshold.}
#'     \item{\code{K}}{Number of IS simulations used.}
#'   }
#'
#' @seealso [asset_is_pvalue_corr()], [asset_is_nonnorm_ind()],
#'   [asset_is_nonnorm_corr()]
#'
#' @references
#' Bhattacharjee, S. et al. (2012). A subset-based approach improves power
#' and interpretation for the combined analysis of genetic association studies
#' of heterogeneous traits. \emph{American Journal of Human Genetics},
#' \bold{90}(5), 821--835.
#' 
#' Shi, J. et al. (2007). Importance sampling for estimating
#' p-values in linkage analysis. \emph{Journal of the American Statistical Association}, 
#' \bold{102}(479), 929--937.
#'
#' Siegmund, D. (1976). Importance sampling in the Monte Carlo study of
#' sequential tests. \emph{The Annals of Statistics}, \bold{4}(4), 673--684.
#'
#' @examples
#' set.seed(42)
#' M <- 5
#' z <- rnorm(M, mean = c(3, 3, 0, 0, 0))   # two studies have signal
#' n <- c(500, 600, 450, 700, 550)
#'
#' ## Two-sided p-value at a single threshold
#' b_obs <- max(abs(z))
#' asset_is_pvalue_ind(z = z, n = n, b = b_obs, K = 1e4, two.sided = TRUE)
#'
#' ## One-sided p-values over a grid (range approximation)
#' b_grid <- seq(3, 6, by = 0.1)
#' asset_is_pvalue_ind(z = z, n = n, b = b_grid, K = 5e4, two.sided = FALSE)
#'
#' ## Case-control: supply effective sample sizes
#' phi <- c(0.3, 0.4, 0.5, 0.35, 0.45)     # case fractions
#' n_eff <- n * phi * (1 - phi)
#' asset_is_pvalue_ind(z = z, n = n_eff, b = b_obs, K = 1e4, two.sided = TRUE)
#'
#' ## Parallel execution
#' \dontrun{
#' asset_is_pvalue_ind(z = z, n = n, b = b_grid, K = 1e5,
#'                 parallel = TRUE, n.cores = 4)
#' }
#'
#' @importFrom stats median rnorm runif
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @export
asset_is_pvalue_ind <- function(z,
                                n,
                                b,
                                K         = 5000L,
                                two.sided = TRUE,
                                xi        = NULL,
                                parallel  = FALSE,
                                n.cores   = parallel::detectCores() - 1L,
                                seed      = NULL) {

  ## ---- validation ----
  M <- length(z)
  if (length(n) != M) stop("'z' and 'n' must have the same length.")
  if (any(n <= 0))    stop("All entries of 'n' must be positive.")
  if (!is.numeric(b) || length(b) == 0) stop("'b' must be a non-empty numeric vector.")
  K <- as.integer(K); if (K < 1L) stop("'K' must be a positive integer.")

  ## ---- tilting parameter at the original scale of b (no halving) ----
  if (is.null(xi)) xi <- median(b)

  ## ---- pre-compute subsets and weight matrix ----
  subsets   <- .generate_subsets(M)
  J         <- length(subsets)
  W         <- .build_weight_matrix_ind(n, subsets)
  b_sorted  <- sort(b)
  B         <- length(b_sorted)
  sqrt_n    <- sqrt(n)
  xi2_half  <- xi^2 / 2
  two_sided <- two.sided          # local copy for closure

  ## ---- parallel setup ----
  if (parallel) {
    n.cores <- min(as.integer(n.cores), K)
    cl <- parallel::makeCluster(n.cores); doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    chunk_sizes <- .make_chunks(K, n.cores)
  } else {
    chunk_sizes <- K
  }
  seeds <- .make_seeds(seed, length(chunk_sizes))

  ## ---- core IS chunk ----
  ## One-sided: pick A_k; stat = max(ZA); denom J terms.
  ## Two-sided (proper mixture tilt): pick A_k and s_k in {+1,-1};
  ##   stat = max(|ZA|); denom 2J terms combining both signs.
  .run <- function(n_this, seed_chunk) {
    if (!is.null(seed_chunk)) set.seed(seed_chunk)
    sum_w <- numeric(B); sum_w2 <- numeric(B)
    for (ii in seq_len(n_this)) {
      Aind <- sample.int(J, 1L)
      A    <- subsets[[Aind]]
      s    <- if (two_sided) (if (runif(1L) < 0.5) 1L else -1L) else 1L
      mu   <- numeric(M)
      mu[A] <- s * xi * sqrt_n[A] / sqrt(sum(n[A]))
      z_sim <- rnorm(M) + mu
      ZA    <- as.numeric(crossprod(z_sim, W))
      if (two_sided) {
        log_d <- .logsumexp(c(xi * ZA - xi2_half, -xi * ZA - xi2_half))
        w     <- 2L * J * exp(-log_d)
        stat  <- max(abs(ZA))
      } else {
        log_d <- .logsumexp(xi * ZA - xi2_half)
        w     <- J * exp(-log_d)
        stat  <- max(ZA)
      }
      k_hit <- findInterval(stat, b_sorted)
      if (k_hit > 0L) {
        idx <- seq_len(k_hit)
        sum_w[idx]  <- sum_w[idx]  + w
        sum_w2[idx] <- sum_w2[idx] + w * w
      }
    }
    list(sum_w = sum_w, sum_w2 = sum_w2)
  }

  ## ---- run ----
  if (parallel) {
    acc <- foreach::foreach(cc = seq_along(chunk_sizes),
                            .combine = .comb_acc,
                            .export  = c(".logsumexp"),
                            .inorder = FALSE) %dopar%
      .run(chunk_sizes[cc], seeds[[cc]])
  } else {
    acc <- .run(chunk_sizes, seeds[[1L]])
  }

  ## ---- summarise ----
  res <- .summarise_acc(acc, K, b, b_sorted, xi)
  df  <- data.frame(b = b, res, xi = xi, K = K,
                    stringsAsFactors = FALSE)
  df$b_anchor <- NULL
  df
}
