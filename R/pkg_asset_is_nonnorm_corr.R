#' Conditional IS for ASSET When Normality Fails: Correlated Cell Types
#'   or Overlapping Studies
#'
#' @title Conditional IS Estimator of ASSET Tail Probabilities for
#'   Non-Gaussian Settings with Correlated Score Statistics
#'
#' @description
#' Estimates \eqn{p(b \mid Y) = P_0(\max_A Z_A > b \mid Y)} (one-sided) or
#' \eqn{p(b \mid Y) = P_0(\max_A |Z_A| > b \mid Y)} (two-sided) when both
#' the normality assumption may fail and score statistics are correlated,
#' e.g., due to correlated gene expression across cell types or overlapping
#' subjects in case-control meta-analysis.
#'
#' The conditional covariance of the score vector \eqn{\mathbf{z}} is
#' estimated as \eqn{\widehat\Sigma_Y = Y^\top Y / N} (after scaling). The
#' optimal GLS statistic for subset \eqn{A} is
#' \deqn{Z_A = \frac{\mathbf{1}_A^\top \Sigma_A^{-1} \mathbf{z}_A}
#'             {\sqrt{\mathbf{1}_A^\top \Sigma_A^{-1} \mathbf{1}_A}},}
#' which admits a linear genotype representation
#' \eqn{Z_A = \sum_i \omega_{i,A}(g_i - \tilde g_i)} with GLS-adjusted
#' weights. The IS algorithm then follows [asset_is_nonnorm_ind()]:
#' Newton solve for \eqn{\xi_A}, trinomial genotype simulation under the
#' tilted HWE distribution, and the mixture IS weight
#' \eqn{(2^C-1) / \sum_A \exp\{\xi_A Z_A - \phi_A(\xi_A)\}}.
#'
#' @details
#' \strong{Covariate adjustment.}  When \code{tau} and \code{g_center} are
#' supplied, the adjusted omega formula
#' \eqn{\omega_{i,A} = N^{-1/2}\sum_{c \in A}a_{c,A}\,y'_{ic}/\tau_c /
#' \sqrt{\mathbf{1}_A^\top a_A}} is used, where \eqn{y'_{ic}} are
#' expression residuals after regressing out covariates.
#'
#' \strong{Scaling.} The formula \eqn{\widehat\Sigma_Y = Y^\top Y/N}
#' requires zero-mean, unit-variance columns. Set \code{scale_Y = TRUE}
#' (default) to apply \code{base::scale()} internally.
#'
#' \strong{Supplying \code{Sigma_z}.} A pre-computed or reference-panel
#' correlation matrix can be passed directly via \code{Sigma_z}, bypassing
#' the empirical estimate.
#'
#' \strong{Modes.} The same three analysis modes as [asset_is_nonnorm_ind()]
#' are supported (eQTL, covariate-adjusted eQTL, case-control), now with
#' GLS weights derived from \eqn{\widehat\Sigma_Y}.
#'
#' @param Y Numeric \eqn{N \times C} matrix. See [asset_is_nonnorm_ind()]
#'   for mode-specific conventions. Scaled internally if
#'   \code{scale_Y = TRUE}.
#' @param f Numeric scalar or length-\eqn{C} vector.  Minor allele
#'   frequency (MAF).  A scalar applies a common MAF to all cell types /
#'   studies; a vector applies study-specific MAFs (multi-ancestry
#'   case-control).  All values must be in \eqn{(0, 0.5]}.
#' @param b Numeric scalar or vector of threshold values.
#' @param study_id Integer vector of length \eqn{N}. Activates
#'   case-control mode. Default \code{NULL}.
#' @param g_center Numeric vector of length \eqn{N}. Per-subject genotype
#'   centering. Default \code{NULL} uses \eqn{2f_{c_i}}.
#' @param tau Numeric positive vector of length \eqn{C}. Scale factors for
#'   covariate-adjusted eQTL. Default \code{NULL}.
#' @param Sigma_z Numeric \eqn{C \times C} matrix.  Conditional covariance
#'   (correlation) matrix of \eqn{(z_1,\ldots,z_C)} given \eqn{Y}.
#'   Defaults to \code{NULL}, in which case \eqn{\widehat{\Sigma}_Y =
#'   Y^\top Y / N} is computed internally after scaling.  Supply a
#'   pre-computed or externally estimated correlation matrix to override the
#'   default (e.g., from a previous run, or estimated from a larger
#'   reference panel).  Must be symmetric positive definite with unit
#'   diagonal.
#' @param scale_Y Logical.  If \code{TRUE} (default), \code{Y} is centered
#'   and scaled to mean zero and unit variance per column.  Set to \code{FALSE}
#'   only if \code{Y} has already been pre-processed.
#' @param ridge Positive numeric. Ridge added to \eqn{\Sigma_A} before
#'   inversion. Default \code{1e-10}.
#' @param two.sided Logical. Default \code{TRUE}.
#' @param K Positive integer. Number of IS simulations. Default \code{5000L}.
#' @param b_anchor Numeric scalar.  The anchor threshold \eqn{b_0} at which
#'   \eqn{\xi_{A,\pm}} are solved.  Defaults to \code{NULL}, using
#'   \eqn{b_0 = \mathrm{median}(b)}.  Ignored for scalar \code{b}.
#' @param newton_tol Numeric. Newton tolerance. Default \code{1e-10}.
#' @param newton_max_iter Integer. Maximum Newton iterations. Default
#'   \code{100L}.
#' @param parallel Logical. Default \code{FALSE}.
#' @param n.cores Positive integer. Default \code{parallel::detectCores()-1}.
#' @param seed Integer or \code{NULL}. Default \code{NULL}.
#'
#' @return A data frame with one row per threshold in \code{b} and columns:
#'   \describe{
#'     \item{\code{b}}{Threshold value as supplied.}
#'     \item{\code{p.est}}{IS point estimate of \eqn{p(b \mid Y)}.}
#'     \item{\code{p.se}}{Estimated standard error of \code{p.est}.}
#'     \item{\code{p.lower}}{Lower bound of 95\% Wald CI, truncated to
#'       \eqn{[0,1]}.}
#'     \item{\code{p.upper}}{Upper bound of 95\% Wald CI, truncated to
#'       \eqn{[0,1]}.}
#'     \item{\code{b_anchor}}{Anchor threshold used for Newton solve.}
#'     \item{\code{K}}{Number of IS simulations used.}
#'   }
#'
#' @seealso [asset_is_nonnorm_ind()], [asset_is_pvalue_corr()]
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
#' N <- 200; C <- 5; f <- 0.2; rho <- 0.4
#'
#' ## Correlated expression data (compound symmetry)
#' Sigma_expr <- (1 - rho) * diag(C) + rho
#' Y <- MASS::mvrnorm(N, mu = rep(0, C), Sigma = Sigma_expr)
#'
#' ## --- eQTL mode (correlated cell types) ---
#' b_obs <- 4.5
#' asset_is_nonnorm_corr(Y = Y, f = f, b = b_obs, K = 2000L, two.sided = TRUE)
#'
#' ## Range of thresholds
#' b_grid <- seq(3.5, 6.0, by = 0.1)
#' asset_is_nonnorm_corr(Y = Y, f = f, b = b_grid, K = 1e4, two.sided = TRUE)
#'
#' ## Supply pre-computed Sigma_z (e.g., from a reference panel)
#' Sigma_z_ref <- cor(Y)
#' asset_is_nonnorm_corr(Y = Y, f = f, b = b_obs, Sigma_z = Sigma_z_ref, K = 2000L)
#'
#' ## --- eQTL with covariate adjustment ---
#' tau_c     <- apply(scale(Y), 2, sd) * sqrt(2 * f * (1 - f))
#' g_center  <- rep(2 * f, N)
#' asset_is_nonnorm_corr(Y = Y, f = f, b = b_obs,
#'                       tau = tau_c, g_center = g_center, K = 2000L)
#'
#' ## --- Case-control with sample overlap (same ancestry) ---
#' study_id <- sample(1:C, N, replace = TRUE)
#' phi_c    <- runif(C, 0.2, 0.4)
#' y_raw    <- rbinom(N, 1, prob = phi_c[study_id])
#' Ycc      <- matrix(0, N, C)
#' for (i in seq_len(N)) Ycc[i, study_id[i]] <- y_raw[i] - phi_c[study_id[i]]
#' asset_is_nonnorm_corr(Y = Ycc, f = f, b = b_obs,
#'                       study_id = study_id, K = 2000L)
#'
#' ## --- Multi-ancestry case-control ---
#' f_c <- runif(C, 0.1, 0.4)
#' asset_is_nonnorm_corr(Y = Ycc, f = f_c, b = b_obs,
#'                       study_id = study_id, K = 2000L)
#'
#' ## --- Parallel execution ---
#' \dontrun{
#' asset_is_nonnorm_corr(Y = Y, f = f, b = b_grid,
#'                       K = 1e5, parallel = TRUE, n.cores = 4)
#' }
#'
#' @importFrom stats median rnorm runif
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @export
asset_is_nonnorm_corr <- function(Y,
                                  f,
                                  b,
                                  study_id        = NULL,
                                  g_center        = NULL,
                                  tau             = NULL,
                                  Sigma_z         = NULL,
                                  scale_Y         = TRUE,
                                  ridge           = 1e-10,
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
  if (length(f) != C) stop("'f' must be scalar or length ncol(Y).")
  if (any(f <= 0) || any(f >= 1)) stop("'f' must be in (0,1).")
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
      stop("'tau' must be positive with length ncol(Y).")
  }
  if (!is.null(Sigma_z)) {
    if (!is.matrix(Sigma_z) || nrow(Sigma_z) != C || ncol(Sigma_z) != C)
      stop("'Sigma_z' must be C x C.")
    if (!isTRUE(all.equal(Sigma_z, t(Sigma_z), tolerance = 1e-8)))
      stop("'Sigma_z' must be symmetric.")
  }
  if (!is.numeric(b) || length(b) == 0) stop("'b' must be non-empty numeric.")
  K <- as.integer(K); if (K < 1L) stop("'K' must be a positive integer.")

  ## ---- scale Y, compute Sigma_z ----
  if (scale_Y) Y <- scale(Y)
  if (is.null(Sigma_z)) {
    Sigma_z <- crossprod(Y) / N
    Sigma_z <- (Sigma_z + t(Sigma_z)) / 2
    diag(Sigma_z) <- 1
  }

  ## ---- anchor threshold at the original scale of b (no halving) ----
  b0 <- if (length(b) == 1L) b
        else if (!is.null(b_anchor)) b_anchor
        else median(b)

  ## ---- per-subject frequency and centering ----
  f_subj <- if (is.null(study_id)) rep(f[1L], N) else f[study_id]
  if (is.null(g_center)) g_center <- 2 * f_subj

  ## ---- subsets, GLS weights, projected matrix ----
  subsets <- .generate_subsets(C)
  J       <- length(subsets)
  Wmat    <- .build_weight_matrix_gls(Sigma_z, subsets, C, ridge)
  Sproj   <- Y %*% Wmat   # N x J

  ## ---- omega vectors ----
  mode <- if (!is.null(study_id)) "cc"
          else if (!is.null(tau))  "eqtl_cov"
          else                     "eqtl"

  if (mode == "eqtl") {
    sigma_g    <- sqrt(2 * f[1L] * (1 - f[1L]))
    denom_z    <- sqrt(N) * sigma_g
    omega_list <- lapply(seq_len(J), function(j) Sproj[, j] / denom_z)

  } else if (mode == "eqtl_cov") {
    Y_sc       <- sweep(Y, 2, tau, FUN = "/")
    Sproj_cov  <- Y_sc %*% Wmat
    omega_list <- lapply(seq_len(J), function(j) Sproj_cov[, j] / sqrt(N))

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
    omega_list <- lapply(seq_len(J), function(j) {
      A     <- subsets[[j]]
      denom <- sqrt(sum(V_c[A]))
      if (denom <= 0) return(rep(0, N))
      w_ci  <- Wmat[cbind(study_id, rep(j, N))]
      w_ci * numer_i / denom
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
      z_vec  <- as.numeric(crossprod(Y, gc_vec)) /
                  (sqrt(N) * sqrt(2 * f * (1 - f)))
      ZA     <- as.numeric(crossprod(Wmat, z_vec))
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
