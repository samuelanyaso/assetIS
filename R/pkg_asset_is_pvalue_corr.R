#' Importance Sampling for ASSET: Correlated Studies, Gaussian Approximation
#'
#' @title Gaussian IS Estimator of ASSET Tail Probabilities for
#'   Correlated or Overlapping Studies
#'
#' @description
#' Estimates \eqn{p(b) = P_0(\max_A Z_A > b)} (one-sided) or
#' \eqn{p(b) = P_0(\max_A |Z_A| > b)} (two-sided) for the ASSET statistic
#' when the study-level score statistics \eqn{(z_1,\ldots,z_M)} are
#' correlated, e.g., due to overlapping subjects across studies. The optimal
#' meta-analysis statistic for subset \eqn{A} is
#' \deqn{Z_A = \frac{N_A^\top \Sigma_A^{-1} z_A}
#'             {\sqrt{N_A^\top \Sigma_A^{-1} N_A}},}
#' where \eqn{N_A = (\sqrt{n_m})_{m \in A}} and \eqn{\Sigma_A} is the
#' principal submatrix of the \eqn{M \times M} inter-study correlation matrix
#' \eqn{\Sigma}. The IS algorithm tilts the joint MVN distribution via
#' \eqn{dP_{\xi,A} = \exp(\xi Z_A - \xi^2/2)\,dP_0}.
#'
#' @details
#' \strong{Case-control studies.} Replace \eqn{n_m} with the effective
#' sample size \eqn{n_m^{\mathrm{eff}} = n_m \phi_m (1-\phi_m)}.
#'
#' \strong{Numerical stability.} If \eqn{\Sigma} is not numerically positive
#' definite, the function attempts regularisation via \code{Matrix::nearPD}
#' followed by diagonal jitter, issuing a warning if applied.
#'
#' For independent studies use [asset_is_pvalue_ind()]. For non-Gaussian
#' settings use [asset_is_nonnorm_ind()] or [asset_is_nonnorm_corr()].
#'
#' @param z Numeric vector of length \eqn{M}. Observed study-level score
#'   statistics.
#' @param n Numeric vector of length \eqn{M}. Study sample sizes. For
#'   case-control studies supply effective sample sizes
#'   \eqn{n_m^{\mathrm{eff}} = n_m \phi_m (1-\phi_m)}.
#' @param Sigma Numeric \eqn{M \times M} inter-study correlation matrix
#'   under the null. Must be symmetric positive definite with unit diagonal.
#' @param b Numeric scalar or vector of threshold values.
#' @param K Positive integer. Number of IS simulations. Default \code{5000L}.
#' @param two.sided Logical. If \code{TRUE} (default) compute the two-sided
#'   p-value.
#' @param xi Numeric scalar. Tilting parameter. Defaults to
#'   \code{median(b)} (after two-sided halving).
#' @param parallel Logical. Default \code{FALSE}.
#' @param n.cores Positive integer. Default \code{parallel::detectCores()-1}.
#' @param seed Integer or \code{NULL}. Default \code{NULL}.
#'
#' @return A data frame with one row per value in \code{b} and columns:
#'   \describe{
#'     \item{\code{b}}{Threshold value (as supplied by the user).}
#'     \item{\code{p.est}}{IS point estimate of \eqn{p(b)}.}
#'     \item{\code{p.se}}{Estimated standard error of \code{p.est}.}
#'     \item{\code{p.lower}}{Lower bound of 95\% Wald CI,
#'       \eqn{\hat{p} - 1.96\,\widehat{\mathrm{SE}}}, truncated to \eqn{[0,1]}.}
#'     \item{\code{p.upper}}{Upper bound of 95\% Wald CI, truncated to
#'       \eqn{[0,1]}.}
#'     \item{\code{xi}}{Tilting parameter used.}
#'     \item{\code{K}}{Number of IS simulations used.}
#'   }
#' 
#'
#' @seealso [asset_is_pvalue_ind()], [asset_is_nonnorm_ind()],
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
#' M   <- 5
#' rho <- 0.4
#' Sigma <- outer(1:M, 1:M, function(i, j) rho^abs(i - j))
#' z   <- rnorm(M, mean = c(3, 3, 0, 0, 0))
#' n   <- c(500, 600, 450, 700, 550)
#'
#' ## Two-sided p-value at a single threshold
#' b_obs <- max(abs(z))
#' asset_is_pvalue_corr(z = z, n = n, Sigma = Sigma, b = b_obs,
#'                      K = 1e4, two.sided = TRUE)
#'
#' ## One-sided p-values over a grid (range approximation)
#' b_grid <- seq(3, 6, by = 0.1)
#' asset_is_pvalue_corr(z = z, n = n, Sigma = Sigma, b = b_grid,
#'                      K = 5e4, two.sided = FALSE)
#'
#' ## Case-control: supply effective sample sizes
#' phi   <- c(0.3, 0.4, 0.5, 0.35, 0.45)
#' n_eff <- n * phi * (1 - phi)
#' asset_is_pvalue_corr(z = z, n = n_eff, Sigma = Sigma, b = b_obs,
#'                      K = 1e4, two.sided = TRUE)
#'
#' ## Parallel execution
#' \dontrun{
#' asset_is_pvalue_corr(z = z, n = n, Sigma = Sigma, b = b_grid,
#'                      K = 1e5, parallel = TRUE, n.cores = 4)
#' }
#'
#' @importFrom stats median rnorm runif
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @export
asset_is_pvalue_corr <- function(z,
                                 n,
                                 Sigma,
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
  if (!is.matrix(Sigma) || nrow(Sigma) != M || ncol(Sigma) != M)
    stop("'Sigma' must be an M x M matrix.")
  if (!isTRUE(all.equal(Sigma, t(Sigma), tolerance = 1e-8)))
    stop("'Sigma' must be symmetric.")
  if (!is.numeric(b) || length(b) == 0) stop("'b' must be non-empty numeric.")
  K <- as.integer(K); if (K < 1L) stop("'K' must be a positive integer.")

  ## ---- tilting parameter at the original scale of b (no halving) ----
  if (is.null(xi)) xi <- median(b)

  ## ---- pre-compute ----
  subsets   <- .generate_subsets(M)
  J         <- length(subsets)
  W         <- .build_weight_matrix_gls(Sigma, subsets, M, ridge = 0)
  SigW      <- Sigma %*% W
  U         <- .safe_chol(Sigma)
  b_sorted  <- sort(b)
  B         <- length(b_sorted)
  xi2_half  <- xi^2 / 2
  two_sided <- two.sided

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
  ## One-sided: pick A_k; mu = xi * SigW[, A_k]; stat = max(ZA); denom J terms.
  ## Two-sided (proper mixture tilt): pick A_k and s_k in {+1,-1};
  ##   mu = s_k * xi * SigW[, A_k]; stat = max(|ZA|); denom 2J terms.
  .run <- function(n_this, seed_chunk) {
    if (!is.null(seed_chunk)) set.seed(seed_chunk)
    sum_w <- numeric(B); sum_w2 <- numeric(B)
    for (ii in seq_len(n_this)) {
      Aind  <- sample.int(J, 1L)
      s     <- if (two_sided) (if (runif(1L) < 0.5) 1L else -1L) else 1L
      mu    <- s * xi * SigW[, Aind]
      z_sim <- mu + as.numeric(crossprod(U, rnorm(M)))
      ZA    <- as.numeric(crossprod(W, z_sim))
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
  df  <- data.frame(b = b, res, xi = xi, K = K, stringsAsFactors = FALSE)
  df$b_anchor <- NULL
  df
}
