# Sensitivity analysis for AIPE on a mixed-effects fixed effect.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Mixed-Effects Fixed Effect
#'
#' @description
#' Quantifies how much misspecification of the variance components
#' (\eqn{\sigma^2_Y}, \eqn{\sigma^2_X}, and the intraclass correlation
#' \eqn{\mathrm{icc}}) distorts an AIPE-based sample size plan for a
#' cluster-level fixed effect under a two-level random-intercept model.
#' On each replication the function simulates \emph{K} clusters of
#' \code{cluster_size} observations each from
#' \deqn{Y_{ki} = \beta\,X_k + u_k + \epsilon_{ki},}
#' with \eqn{X_k \sim N(0, \sigma^2_X)}, \eqn{u_k \sim N(0,
#' \mathrm{icc}\cdot\sigma^2_Y)}, and \eqn{\epsilon_{ki} \sim N(0,
#' (1 - \mathrm{icc})\sigma^2_Y)}. The model is then refit by either
#' \code{lme4::lmer} (if available) or by GLS-by-cluster aggregation,
#' and a Wald CI on the fixed effect is recorded.
#'
#' @param true_sigma2_y Population total variance of \emph{Y}.
#' @param true_sigma2_x Population variance of the cluster-level
#'   predictor.
#' @param true_icc Population intraclass correlation (between-cluster
#'   share of total variance).
#' @param true_beta Population fixed-effect slope (default \code{0}).
#' @param estimated_sigma2_y,estimated_sigma2_x,estimated_icc Planning
#'   values passed to \code{\link{ss_aipe_mixed_effects}}; supply
#'   all three or \code{specified_K} but not both.
#' @param width Desired full width of the CI on the fixed effect.
#' @param cluster_size Number of observations per cluster (assumed
#'   balanced).
#' @param specified_K Number of clusters to evaluate.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param save Logical. Save per-replication CSV.
#' @param filename Path used when \code{save = TRUE}.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   the realized fixed-effect estimate and CI width, the proportion of
#'   intervals at or below \code{width}, tail-specific and overall
#'   non-coverage of \code{true_beta}, and the input echoes.
#'
#' @references
#' McNeish, D., & Kelley, K. (2019). Fixed effects versus mixed effects
#'   models for clustered data: Reviewing the approaches, disentangling
#'   the differences, and making recommendations. \emph{Psychological
#'   Methods, 24}, 20--35. \doi{10.1037/met0000182}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapters 15 and 16 on mixed-effects
#'   models.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_mixed_effects}}, \code{\link{icc_lmer}}
#'
#' @examples
#' # Monte Carlo sensitivity check, reduced sizes for a fast example.
#' set.seed(113)
#' ss_aipe_mixed_effects_sensitivity(
#'   true_sigma2_y = 1, true_sigma2_x = 1, true_icc = 0.10,
#'   true_beta = 0.30,
#'   specified_K = 25, cluster_size = 10,
#'   width = 0.40,
#'   G = 25, print_iter = FALSE
#' )
#'
#' @keywords design multivariate
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family AIPE sample size planning
#' @family mixed models
#'
#' @export
ss_aipe_mixed_effects_sensitivity <- function(true_sigma2_y = NULL,
                                                   true_sigma2_x = NULL,
                                                   true_icc = NULL,
                                                   true_beta = 0,
                                                   estimated_sigma2_y = NULL,
                                                   estimated_sigma2_x = NULL,
                                                   estimated_icc = NULL,
                                                   width,
                                                   cluster_size = 20L,
                                                   specified_K = NULL,
                                                   conf_level = 0.95,
                                                   G = 1000, print_iter = FALSE,
                                                   save = FALSE,
                                                   filename = "ss_aipe_mixed_effects_sensitivity_result.csv") {
  est_supplied <- !is.null(estimated_sigma2_y) && !is.null(estimated_sigma2_x) &&
                  !is.null(estimated_icc)
  if (!est_supplied && is.null(specified_K))
    stop("Supply all three estimated variance components, or 'specified_K'.", call. = FALSE)
  if (est_supplied && !is.null(specified_K))
    stop("Supply estimated variance components or 'specified_K', not both.", call. = FALSE)
  if (is.null(true_sigma2_y) || true_sigma2_y <= 0)
    stop("'true_sigma2_y' must be a positive number.", call. = FALSE)
  if (is.null(true_sigma2_x) || true_sigma2_x <= 0)
    stop("'true_sigma2_x' must be a positive number.", call. = FALSE)
  if (is.null(true_icc) || true_icc < 0 || true_icc >= 1)
    stop("'true_icc' must be in [0, 1).", call. = FALSE)

  if (est_supplied) {
    plan <- ss_aipe_mixed_effects(sigma2_y = estimated_sigma2_y,
                                  sigma2_x = estimated_sigma2_x,
                                  icc = estimated_icc,
                                  width = width, cluster_size = cluster_size,
                                  conf_level = conf_level)
    K <- plan$value[plan$term == "necessary_n_clusters"]
  } else {
    K <- specified_K
  }
  K <- as.integer(K)
  if (K < 5) stop("Resolved number of clusters is < 5.", call. = FALSE)
  m <- as.integer(cluster_size)

  has_lme4 <- requireNamespace("lme4", quietly = TRUE)

  sigma2_u <- true_icc * true_sigma2_y
  sigma2_e <- (1 - true_icc) * true_sigma2_y
  sd_x     <- sqrt(true_sigma2_x)
  sd_u     <- sqrt(sigma2_u)
  sd_e     <- sqrt(sigma2_e)

  beta_hat <- numeric(G)
  ci_lo    <- numeric(G)
  ci_hi    <- numeric(G)
  ci_w     <- numeric(G)
  tI_lower <- logical(G)
  tI_upper <- logical(G)

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    X_k <- stats::rnorm(K, mean = 0, sd = sd_x)
    u_k <- stats::rnorm(K, mean = 0, sd = sd_u)
    X_long  <- rep(X_k, each = m)
    u_long  <- rep(u_k, each = m)
    eps     <- stats::rnorm(K * m, mean = 0, sd = sd_e)
    Y       <- true_beta * X_long + u_long + eps
    cluster <- factor(rep(seq_len(K), each = m))

    if (has_lme4) {
      fit <- tryCatch(
        lme4::lmer(Y ~ X_long + (1 | cluster), REML = TRUE),
        error = function(e) NULL,
        warning = function(w) suppressWarnings(
          lme4::lmer(Y ~ X_long + (1 | cluster), REML = TRUE))
      )
      if (is.null(fit) || !inherits(fit, "merMod")) {
        bhat <- NA_real_; lo <- NA_real_; hi <- NA_real_
      } else {
        coefs <- lme4::fixef(fit)
        bhat  <- unname(coefs["X_long"])
        vc    <- as.matrix(stats::vcov(fit))
        se    <- sqrt(vc["X_long", "X_long"])
        # Conservative cluster-level df = K - 2 (number of clusters minus
        # the two fixed effects). A more refined Satterthwaite df would
        # require the lmerTest package; the conservative df errs on the
        # side of wider intervals, which is the correct direction for a
        # sensitivity analysis.
        crit <- stats::qt(1 - (1 - conf_level) / 2, df = K - 2)
        lo <- bhat - crit * se
        hi <- bhat + crit * se
      }
    } else {
      # GLS-by-aggregation fallback: collapse to cluster means and regress.
      Y_bar  <- tapply(Y, cluster, mean)
      X_clus <- as.numeric(tapply(X_long, cluster, mean))
      fit <- stats::lm(Y_bar ~ X_clus)
      bhat <- unname(stats::coef(fit)[2])
      se   <- sqrt(stats::vcov(fit)[2, 2])
      crit <- stats::qt(1 - (1 - conf_level) / 2, df = K - 2)
      lo <- bhat - crit * se
      hi <- bhat + crit * se
    }
    beta_hat[g] <- bhat
    ci_lo[g]    <- lo
    ci_hi[g]    <- hi
    ci_w[g]     <- hi - lo
    tI_lower[g] <- true_beta < lo
    tI_upper[g] <- true_beta > hi
  }

  if (isTRUE(save)) {
    per_rep <- data.frame(beta = beta_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_w,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_beta", "median_beta", "sd_beta",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "n_clusters", "cluster_size", "total_N",
              "true_sigma2_y", "true_sigma2_x", "true_icc", "true_beta",
              "estimated_sigma2_y", "estimated_sigma2_x", "estimated_icc",
              "width", "conf_level"),
    value = c(mean(beta_hat, na.rm = TRUE),
              stats::median(beta_hat, na.rm = TRUE),
              stats::sd(beta_hat, na.rm = TRUE),
              mean(ci_w, na.rm = TRUE),
              stats::median(ci_w, na.rm = TRUE),
              stats::sd(ci_w, na.rm = TRUE),
              mean(ci_w <= width, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE),
              mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              K, m, K * m,
              true_sigma2_y, true_sigma2_x, true_icc, true_beta,
              if (is.null(estimated_sigma2_y)) NA_real_ else estimated_sigma2_y,
              if (is.null(estimated_sigma2_x)) NA_real_ else estimated_sigma2_x,
              if (is.null(estimated_icc))      NA_real_ else estimated_icc,
              width, conf_level),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
