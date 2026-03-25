#' Conditional IS for ASSET When Normality Fails: Independent Cell Types
#'   or Studies
#'
#' @title Conditional IS Estimator of ASSET Tail Probabilities for
#'   Non-Gaussian, Independent Settings
#'
#' @description
#' Estimates the conditional tail probability
#' \eqn{p(b \mid Y) = P_0(\max_A Z_A > b \mid Y)} (one-sided) or
#' \eqn{p(b \mid Y) = P_0(\max_A |Z_A| > b \mid Y)} (two-sided) for the
#' ASSET statistic when the normality assumption on score statistics may be
#' violated and cell types / studies are independent. This arises in
#' single-cell eQTL mapping with small sample sizes or meta-analysis of
#' rare diseases.
#'
#' Inference is conditional on the observed phenotype matrix \eqn{Y};
#' the only source of randomness is the genotype vector
#' \eqn{(g_1,\ldots,g_N)}. Genotypes are simulated from a discrete
#' (trinomial) distribution under Hardy--Weinberg equilibrium (HWE), with
#' per-subject exponential tilting. The tilting parameter \eqn{\xi_A} for
#' each subset is chosen to satisfy \eqn{\phi_A'(\xi_A) = b_0} via
#' Newton's method.
#'
#' The function supports three analysis modes selected automatically:
#' \describe{
#'   \item{eQTL}{All \eqn{N} subjects in all \eqn{C} cell types
#'     (\code{study_id = NULL}, \code{tau = NULL}).}
#'   \item{eQTL with covariate adjustment}{\code{tau} and \code{g_center}
#'     supplied; \code{Y} contains expression residuals.}
#'   \item{Case-control meta-analysis}{\code{study_id} supplied; subjects
#'     nested within \eqn{C} studies; multi-ancestry via a length-\eqn{C}
#'     vector \code{f}.}
#' }
#'
#' @details
#' \strong{Two-sided test.}  A proper two-sided mixture tilt is used rather
#' than the heuristic of halving the threshold and doubling the estimate.
#' Because the CGF \eqn{\phi_A(\xi)} is not symmetric in \eqn{\xi} for
#' general (data-dependent) weights \eqn{\omega_{i,A}}, two per-subset
#' tilting parameters are solved:
#' \eqn{\xi_{A,+}} satisfying \eqn{\phi_A'(\xi_{A,+}) = +b_0} (tilting
#' toward the upper tail) and \eqn{\xi_{A,-}} satisfying
#' \eqn{\phi_A'(\xi_{A,-}) = -b_0} (tilting toward the lower tail).
#' For each simulation, a driving subset \eqn{A_k} and a sign
#' \eqn{s_k \in \{+1,-1\}} are chosen uniformly at random.  Genotypes are
#' drawn under the tilted distribution corresponding to \eqn{\xi_{A_k,s_k}},
#' and the IS weight denominator involves \eqn{2(2^C-1)} terms:
#' \deqn{\frac{2(2^C-1)}{\sum_A \left[
#'   \exp\{\xi_{A,+} Z_A - \phi_A(\xi_{A,+})\} +
#'   \exp\{\xi_{A,-} Z_A - \phi_A(\xi_{A,-})\}
#' \right]} \mathcal{I}\!\left\{\max_A |Z_A| > b\right\}.}
#'
#' \strong{Range of thresholds.}  When \code{b} is a vector, \eqn{\xi_{A,\pm}}
#' are solved once at the anchor \eqn{b_0 = \mathrm{median}(b)} (or
#' \code{b_anchor}) at the original scale of \eqn{b}.
#'
#' \strong{Scaling of \eqn{Y}.}  Each column of \eqn{Y} should have mean
#' zero and unit variance.  If \code{scale_Y = TRUE} (default), the function
#' applies \code{base::scale} to \eqn{Y} internally.
#' 
#' @param Y Numeric matrix of dimension \eqn{N \times C}, centered and
#'   scaled (mean zero, unit variance per column) or raw (if
#'   \code{scale_Y = TRUE}).  For eQTL: the expression matrix
#'   (or residual matrix after covariate adjustment, with \code{tau}
#'   supplied). For case-control: a matrix whose entry \code{Y[i,
#'   study_id[i]]} contains the centered binary outcome
#'   \eqn{y_{ic} - \bar{y}_c} (or the covariate residual
#'   \eqn{y_{ic} - \hat{\mu}_{ic}}) for subject \eqn{i}; entries in other
#'   columns are ignored.
#' @param f Numeric scalar or length-\eqn{C} vector. Minor allele
#'   frequency (MAF). Scalar for common MAF; length-\eqn{C} vector for
#'   study-specific MAFs (multi-ancestry case-control). Must be in
#'   \eqn{(0, 0.5]}.
#' @param b Numeric scalar or vector of threshold values. When a vector,
#'   \eqn{\xi_A} is solved once at the anchor
#'   \eqn{b_0 = \mathrm{median}(\mathtt{b})} (or \code{b_anchor}).
#' @param study_id Integer vector of length \eqn{N} with values in
#'   \eqn{\{1,\ldots,C\}}. Activates case-control mode. Default \code{NULL}
#'   (eQTL mode).
#' @param g_center Numeric vector of length \eqn{N}. Per-subject genotype
#'   centering values \eqn{\tilde g_i}. Default \code{NULL} uses
#'   \eqn{2f_{c_i}}.
#' @param tau Numeric positive vector of length \eqn{C}. Per-cell-type
#'   scale factors for covariate-adjusted eQTL weights. Default \code{NULL}.
#' @param scale_Y Logical.  If \code{TRUE} (default), \code{Y} is centered
#'   and scaled to mean zero and unit variance per column.  Set to \code{FALSE}
#'   only if \code{Y} has already been pre-processed.
#' @param two.sided Logical. Default \code{TRUE}.
#' @param K Positive integer. Number of IS simulations. Default \code{5000L}.
#' @param b_anchor Numeric scalar.  The anchor threshold \eqn{b_0} at which
#'   \eqn{\xi_{A,\pm}} are solved.  Defaults to \code{NULL}, using
#'   \eqn{b_0 = \mathrm{median}(b)}.  Ignored for scalar \code{b}.
#' @param newton_tol Positive numeric. Newton convergence tolerance.
#'   Default \code{1e-10}.
#' @param newton_max_iter Positive integer. Maximum Newton iterations.
#'   Default \code{100L}.
#' @param parallel Logical. Default \code{FALSE}.
#' @param n.cores Positive integer. Default \code{parallel::detectCores()-1}.
#' @param seed Integer or \code{NULL}. Default \code{NULL}.
#'
#' @return A data frame with one row per threshold in \code{b} and columns:
#'   \describe{
#'     \item{\code{b}}{Threshold value as supplied.}
#'     \item{\code{p.est}}{IS point estimate of \eqn{p(b \mid Y)}.}
#'     \item{\code{p.se}}{Estimated standard error of \code{p.est}.}
#'     \item{\code{p.lower}}{Lower bound of 95\% Wald CI,
#'       \eqn{\hat{p} - 1.96\,\widehat{\mathrm{SE}}}, truncated to
#'       \eqn{[0,1]}.}
#'     \item{\code{p.upper}}{Upper bound of 95\% Wald CI, truncated to
#'       \eqn{[0,1]}.}
#'     \item{\code{b_anchor}}{Anchor threshold used for Newton solve.}
#'     \item{\code{K}}{Number of IS simulations used.}
#'   }
#'
#' @seealso [asset_is_nonnorm_corr()], [asset_is_pvalue_ind()],
#'   [asset_is_pvalue_corr()]
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
#' N <- 200; C <- 5
#' f <- 0.2
#'
#' ## --- eQTL mode ---
#' Y <- matrix(rnorm(N * C), N, C)
#' b_obs <- 4.5
#' asset_is_nonnorm_ind(Y = Y, f = f, b = b_obs, K = 2000L, two.sided = TRUE)
#'
#' ## Range of thresholds (range approximation)
#' b_grid <- seq(3.5, 6.0, by = 0.1)
#' asset_is_nonnorm_ind(Y = Y, f = f, b = b_grid, K = 1e4, two.sided = TRUE)
#'
#' ## --- eQTL with covariate adjustment ---
#' ## Y contains expression residuals y'_ic after regressing out covariates
#' tau_c <- apply(scale(Y), 2, sd) * sqrt(2 * f * (1 - f))
#' g_center <- rep(2 * f, N)   # replace with predicted g_i0 in real applications
#' asset_is_nonnorm_ind(Y = Y, f = f, b = b_obs, tau = tau_c,
#'                  g_center = g_center, K = 2000L)
#'
#' ## --- Case-control meta-analysis (same ancestry) ---
#' study_id <- sample(1:C, N, replace = TRUE)
#' phi_c    <- runif(C, 0.2, 0.4)   # case fractions
#' y_raw    <- rbinom(N, 1, prob = phi_c[study_id])
#' Ycc      <- matrix(0, N, C)
#' for (i in seq_len(N)) {
#'   c_i <- study_id[i]
#'   Ycc[i, c_i] <- y_raw[i] - phi_c[c_i]
#' }
#' asset_is_nonnorm_ind(Y = Ycc, f = f, b = b_obs, study_id = study_id, K = 2000L)
#'
#' ## --- Case-control, different ancestries (study-specific MAF) ---
#' f_c <- runif(C, 0.1, 0.4)
#' asset_is_nonnorm_ind(Y = Ycc, f = f_c, b = b_obs, study_id = study_id, K = 2000L)
#'
#' @importFrom stats median rnorm runif
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @export
asset_is_nonnorm_ind <- function(Y,
                                 f,
                                 b,
                                 study_id        = NULL,
                                 g_center        = NULL,
                                 tau             = NULL,
                                 scale_Y         = TRUE,
                                 two.sided       = TRUE,
                                 K               = 5000L,
                                 b_anchor        = NULL,
                                 newton_tol      = 1e-10,
                                 newton_max_iter = 100L,
                                 parallel        = FALSE,
                                 n.cores         = parallel::detectCores() - 1L,
                                 seed            = NULL) {

  ## ---- validation ----
  if (!is.matrix(Y) || !is.numeric(Y)) stop("'Y' must be a numeric matrix.")
  N <- nrow(Y); C <- ncol(Y)
  f <- as.numeric(f); if (length(f) == 1L) f <- rep(f, C)
  if (length(f) != C) stop("'f' must be a scalar or length ncol(Y) vector.")
  if (any(f <= 0) || any(f >= 1)) stop("'f' must be in (0, 1).")
  if (!is.null(study_id)) {
    study_id <- as.integer(study_id)
    if (length(study_id) != N) stop("'study_id' must have length nrow(Y).")
    if (any(study_id < 1L) || any(study_id > C))
      stop("'study_id' values must be in {1,...,ncol(Y)}.")
  }
  if (!is.null(g_center)) {
    g_center <- as.numeric(g_center)
    if (length(g_center) != N) stop("'g_center' must have length nrow(Y).")
  }
  if (!is.null(tau)) {
    tau <- as.numeric(tau)
    if (length(tau) != C || any(tau <= 0))
      stop("'tau' must be a positive vector of length ncol(Y).")
  }
  if (!is.numeric(b) || length(b) == 0) stop("'b' must be non-empty numeric.")
  K <- as.integer(K); if (K < 1L) stop("'K' must be a positive integer.")

  ## ---- Scale Y ----
  if (scale_Y) Y <- scale(Y)
  
  ## ---- anchor threshold at the original scale of b (no halving) ----
  b0 <- if (length(b) == 1L) b
        else if (!is.null(b_anchor)) b_anchor
        else median(b)

  ## ---- per-subject frequency and centering ----
  f_subj   <- if (is.null(study_id)) rep(f[1L], N) else f[study_id]
  if (is.null(g_center)) g_center <- 2 * f_subj

  ## ---- subsets and omega vectors ----
  subsets <- .generate_subsets(C)
  J       <- length(subsets)

  mode <- if (!is.null(study_id)) "cc"
          else if (!is.null(tau))  "eqtl_cov"
          else                     "eqtl"

  if (mode == "eqtl") {
    sigma_g    <- sqrt(2 * f[1L] * (1 - f[1L]))
    omega_list <- lapply(subsets, function(A) {
      s_i <- if (length(A) == 1L) Y[, A] else rowSums(Y[, A, drop = FALSE])
      s_i / (sqrt(N * length(A)) * sigma_g)
    })
  } else if (mode == "eqtl_cov") {
    Y_sc       <- sweep(Y, 2, tau, FUN = "/")
    omega_list <- lapply(subsets, function(A) {
      s_i <- if (length(A) == 1L) Y_sc[, A] else rowSums(Y_sc[, A, drop = FALSE])
      s_i / sqrt(N)
    })
  } else {
    N_c       <- tabulate(study_id, nbins = C)
    sigma2_g  <- 2 * f * (1 - f)
    resid_var <- vapply(seq_len(C), function(c) {
      idx <- which(study_id == c)
      if (length(idx) == 0L) return(0)
      mean(Y[idx, c]^2)
    }, numeric(1L))
    V_c     <- N_c * resid_var * sigma2_g
    numer_i <- vapply(seq_len(N), function(i) Y[i, study_id[i]], numeric(1L))
    omega_list <- lapply(subsets, function(A) {
      denom <- sqrt(sum(V_c[A]))
      if (denom <= 0) return(rep(0, N))
      numer_i / denom
    })
  }

  ## ---- Newton solve for +b0 (all modes) and -b0 (two-sided only) ----
  xphi_pos <- .precompute_xi_phi(omega_list, f_subj, +b0,
                                  newton_tol, newton_max_iter)
  xi_vec_pos  <- xphi_pos$xi_vec
  phi_vec_pos <- xphi_pos$phi_vec

  if (two.sided) {
    xphi_neg <- .precompute_xi_phi(omega_list, f_subj, -b0,
                                    newton_tol, newton_max_iter)
    xi_vec_neg  <- xphi_neg$xi_vec
    phi_vec_neg <- xphi_neg$phi_vec
  } else {
    xi_vec_neg  <- numeric(J)
    phi_vec_neg <- numeric(J)
  }

  ## ---- sorted thresholds at the original scale of b ----
  b_sorted  <- sort(b); B_len <- length(b_sorted)
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
  ## One-sided: pick A_k; use xi_{A_k,+}; stat = max(ZA); denom J terms.
  ## Two-sided (proper mixture tilt): pick A_k and s_k in {+1,-1};
  ##   use xi_{A_k,s_k}; stat = max(|ZA|); denom 2J terms (pos + neg).
  .run <- function(n_this, seed_chunk) {
    if (!is.null(seed_chunk)) set.seed(seed_chunk)
    sum_w <- numeric(B_len); sum_w2 <- numeric(B_len)
    for (ii in seq_len(n_this)) {
      Aind <- sample.int(J, 1L)
      if (two_sided) {
        xi_A <- if (runif(1L) < 0.5) xi_vec_pos[Aind] else xi_vec_neg[Aind]
      } else {
        xi_A <- xi_vec_pos[Aind]
      }
      omega_A <- omega_list[[Aind]]
      xw  <- xi_A * omega_A
      p0  <- (1 - f_subj)^2
      p1  <- 2 * f_subj * (1 - f_subj) * exp(xw)
      p2  <- f_subj^2 * exp(2 * xw)
      s   <- p0 + p1 + p2
      u   <- runif(N) * s
      g   <- integer(N)
      g[u > p0]        <- 1L
      g[u > (p0 + p1)] <- 2L
      gc_vec <- g - g_center
      ZA  <- vapply(seq_len(J), function(j)
        sum(omega_list[[j]] * gc_vec), numeric(1L))
      if (two_sided) {
        log_terms <- c(xi_vec_pos * ZA - phi_vec_pos,
                       xi_vec_neg * ZA - phi_vec_neg)
        log_d     <- .logsumexp(log_terms)
        w         <- 2L * J * exp(-log_d)
        stat      <- max(abs(ZA))
      } else {
        log_terms <- xi_vec_pos * ZA - phi_vec_pos
        log_d     <- .logsumexp(log_terms)
        w         <- J * exp(-log_d)
        stat      <- max(ZA)
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
  res <- .summarise_acc(acc, K, b, b_sorted, b0)
  data.frame(b = b, res, K = K, stringsAsFactors = FALSE)
}
