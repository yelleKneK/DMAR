# Sensitivity analysis for AIPE on a mixed-effects fixed effect.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Mixed-Effects Fixed Effect
#'
#' @description
#' Simulates the design that \code{\link{ss_aipe_mixed_effects}} plans and
#' reports what its confidence interval for the fixed-effect slope
#' delivers: the mean width, the share of intervals no wider than desired,
#' and the coverage. The two functions take the same description of the
#' predictor (\code{predictor_level}, \code{centering}, the intraclass
#' correlation of the predictor), use the same expressions, and form the
#' same interval, so a plan made by the first can be checked by the
#' second, at the planning values or at population values that differ
#' from them.
#'
#' On each replication the function draws \code{n_clusters} clusters
#' (i.e., Level 2 units) of \code{cluster_size} individuals (i.e., Level 1
#' units) from
#' \deqn{y_{ki} = \beta x_{ki} + u_k + e_{ki},}
#' with \eqn{u_k \sim N(0, \rho V)} and \eqn{e_{ki} \sim N(0, (1 - \rho) V)},
#' where \eqn{\beta} is \code{true_reg_coef}, \eqn{V} is
#' \code{true_total_var}, and \eqn{\rho} is \code{true_icc}. An
#' individual-level (i.e., Level 1) predictor is
#' \eqn{x_{ki} = b_k + w_{ki}}, with variance \eqn{\rho_X \sigma^2_X}
#' between clusters and \eqn{(1 - \rho_X)\sigma^2_X} within them
#' (\code{true_var_X}, \code{true_icc_X}); a cluster-level (i.e., Level 2)
#' predictor is \eqn{x_k \sim N(0, \sigma^2_X)}, the same for every
#' individual of a cluster. The predictor is then centered as
#' \code{centering} says, the model is fit, and the interval
#' \eqn{\hat\beta \pm t_{1 - \alpha/2,\, \nu}\,\widehat{SE}(\hat\beta)} is
#' recorded, with the degrees of freedom \eqn{\nu} of
#' \code{\link{ss_aipe_mixed_effects}}.
#'
#' @param true_total_var Population variance of the outcome given the
#'   predictor, the variance of the random intercept plus the residual
#'   variance: the value the data are generated from.
#' @param true_var_X Population variance of the predictor the data are
#'   generated from: over individuals for an individual-level predictor
#'   (between and within clusters together), over clusters for a
#'   cluster-level predictor.
#' @param true_icc Population intraclass correlation of the outcome given
#'   the predictor (the share of \code{true_total_var} between clusters).
#' @param true_reg_coef Population value of the fixed-effect slope, the
#'   regression coefficient of the outcome on the predictor (\eqn{\beta}
#'   in the model above; default \code{0}). The width of the interval does
#'   not depend on it; the coverage is evaluated against it.
#' @param supposed_total_var,supposed_var_X,supposed_icc Supposed
#'   values, the planning values passed to
#'   \code{\link{ss_aipe_mixed_effects}} (not estimates computed from
#'   data), which may differ from the population values above; supply all
#'   three, or
#'   \code{specified_n_clusters}, but not both.
#' @param width Desired full width of the confidence interval for the
#'   slope.
#' @param cluster_size Number of individuals in each cluster (the same in
#'   every cluster).
#' @param specified_n_clusters Number of clusters to evaluate, in place of
#'   the number the planning values give.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param predictor_level,centering The level of the predictor (\code{1},
#'   the default, or \code{2}) and how it is centered
#'   (\code{"cluster_mean"}, \code{"grand_mean"}, or \code{"none"}), as in
#'   \code{\link{ss_aipe_mixed_effects}}. The same values are used to plan
#'   and to analyze each simulated data set.
#' @param true_icc_X Population intraclass correlation of an
#'   individual-level predictor: the share of \code{true_var_X} between
#'   clusters in the data generated. Default \code{0}. Not used with
#'   \code{predictor_level = 2}.
#' @param supposed_icc_X The supposed value of \code{true_icc_X}, the
#'   planning value passed to \code{\link{ss_aipe_mixed_effects}}.
#'   Default: \code{true_icc_X}.
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param filename Optional path for a comma separated file recording
#'   every replication (the fixed effect estimate, the two confidence
#'   limits, the interval width, and two indicators of whether the
#'   interval missed \code{true_reg_coef} below or above): nothing is written
#'   when \code{filename} is \code{NULL} (the default), a new file with a
#'   header row is created otherwise, an existing file at that path is
#'   appended to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @details
#' \strong{How each case is fit.} Two of the three analyses have closed
#' forms in a balanced design and are computed directly. The slope of an
#' individual-level predictor centered at its cluster mean is the
#' within-cluster regression slope, with the residual mean square on
#' \eqn{K(m - 1) - 1} degrees of freedom. The slope of a cluster-level
#' predictor is the slope of the regression of the \eqn{K} cluster means
#' of the outcome on the predictor, on \eqn{K - 2} degrees of freedom. In
#' both, restricted maximum likelihood for the random-intercept model
#' gives the same estimate and standard error when its estimate of the
#' intercept variance is positive. An individual-level predictor that is
#' centered at its grand mean or not centered has no closed form, because
#' its weights depend on the estimated variance components; it is fit by
#' restricted maximum likelihood with \pkg{lme4}, which must then be
#' installed, and its interval uses \eqn{K(m - 1) - 1} degrees of freedom.
#'
#' \strong{What to look at.} \code{mean_ci_width} against \code{width}
#' says whether the plan delivers the width on average, and
#' \code{pct_ci_less_width} how often. \code{total_type_I_error} against
#' \eqn{1 - } \code{conf_level} says whether the interval covers as it
#' should; with few clusters and a cluster-level predictor this is where
#' a problem would show. Replications in which no interval could be formed
#' are counted in the record of the simulation.
#'
#' @return A \code{data.frame} with rows for the mean, median, and
#'   standard deviation of the estimated slope (\code{mean_reg_coef},
#'   \code{median_reg_coef}, \code{sd_reg_coef}) and of the interval width,
#'   the share of intervals no wider than \code{width}
#'   (\code{pct_ci_less_width}), the shares that missed \code{true_reg_coef}
#'   below and above and their sum, the design simulated
#'   (\code{n_clusters}, \code{cluster_size}, \code{total_N},
#'   \code{predictor_level}), the \code{expected_width} that
#'   \code{\link{ss_aipe_mixed_effects}} gives for that design at the
#'   population values, and the inputs. The centering is recorded in
#'   \code{attr(x, "centering")}.
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
#' # Plan, then check the plan (few replications, for a fast example).
#' ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1, true_icc = 0.10,
#'                       width = 0.40, cluster_size = 10)
#' ss_aipe_mixed_effects_sensitivity(
#'   true_total_var = 1, true_var_X = 1, true_icc = 0.10,
#'   true_reg_coef = 0.30,
#'   supposed_total_var = 1, supposed_var_X = 1, supposed_icc = 0.10,
#'   width = 0.40, cluster_size = 10,
#'   G = 50, print_iter = FALSE, seed = 113
#' )
#'
#' # A cluster-level predictor, at a number of clusters chosen by hand.
#' ss_aipe_mixed_effects_sensitivity(
#'   true_total_var = 1, true_var_X = 1, true_icc = 0.10,
#'   true_reg_coef = 0.30,
#'   specified_n_clusters = 25, cluster_size = 10, width = 0.40,
#'   predictor_level = 2, centering = "grand_mean",
#'   G = 50, print_iter = FALSE, seed = 113
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
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export
ss_aipe_mixed_effects_sensitivity <- function(true_total_var = NULL,
                                                   true_var_X = NULL,
                                                   true_icc = NULL,
                                                   true_reg_coef = 0,
                                                   supposed_total_var = NULL,
                                                   supposed_var_X = NULL,
                                                   supposed_icc = NULL,
                                                   width,
                                                   cluster_size = 20L,
                                                   specified_n_clusters = NULL,
                                                   conf_level = 0.95,
                                                   predictor_level = 1,
                                                   centering = c("cluster_mean",
                                                                 "grand_mean",
                                                                 "none"),
                                                   true_icc_X = 0,
                                                   supposed_icc_X = true_icc_X,
                                                   G = 1000, print_iter = FALSE,
                                                   filename = NULL,
                                                   seed = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!missing(width)) .dmar_check_width(width)
  .dmar_check_numbers("ss_aipe_mixed_effects_sensitivity", c("true_total_var", "true_var_X", "true_icc"))
  est_supplied <- !is.null(supposed_total_var) && !is.null(supposed_var_X) &&
                  !is.null(supposed_icc)
  if (!est_supplied && is.null(specified_n_clusters))
    .dmar_refuse_argument(
      "Supply all three supposed values ('supposed_total_var', 'supposed_var_X', and 'supposed_icc'), or 'specified_n_clusters'.",
      arg = "specified_n_clusters")
  if (est_supplied && !is.null(specified_n_clusters))
    .dmar_refuse_argument(
      "Supply the supposed values or 'specified_n_clusters', not both.",
      arg = "specified_n_clusters")
  if (is.null(true_total_var) || true_total_var <= 0)
    .dmar_refuse_argument(
      "'true_total_var' must be a positive number.", arg = "true_total_var")
  if (is.null(true_var_X) || true_var_X <= 0)
    .dmar_refuse_argument(
      "'true_var_X' must be a positive number.", arg = "true_var_X")
  if (is.null(true_icc) || true_icc < 0 || true_icc >= 1)
    .dmar_refuse_argument("'true_icc' must be in [0, 1).", arg = "true_icc")
  if (!is.numeric(true_reg_coef) || length(true_reg_coef) != 1L ||
      !is.finite(true_reg_coef))
    .dmar_refuse_argument(
      "'true_reg_coef' must be a single finite number.", arg = "true_reg_coef")
  centering <- if (missing(centering) && is.numeric(predictor_level) &&
                   length(predictor_level) == 1L &&
                   isTRUE(predictor_level == 2)) {
    "grand_mean"
  } else {
    match.arg(centering)
  }
  design <- .me_design(predictor_level, centering, true_icc_X)
  .me_design(predictor_level, centering, supposed_icc_X,
             icc_arg = "supposed_icc_X")
  needs_lme4 <- design$level == 1L && design$centering != "cluster_mean"
  if (needs_lme4 && !requireNamespace("lme4", quietly = TRUE)) {
    .dmar_refuse(c(
      "An individual-level predictor that is not centered at its cluster ",
      "mean (centering = \"", design$centering, "\") is fit with the lme4 ",
      "package, which is not installed. Install lme4, or use ",
      "centering = \"cluster_mean\", which is computed without it."),
      "package_required", package = "lme4")
  }
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)

  if (est_supplied) {
    plan <- ss_aipe_mixed_effects(true_total_var = supposed_total_var,
                                  true_var_X = supposed_var_X,
                                  true_icc = supposed_icc,
                                  width = width, cluster_size = cluster_size,
                                  conf_level = conf_level,
                                  predictor_level = design$level,
                                  centering = design$centering,
                                  true_icc_X = supposed_icc_X)
    K <- plan$value[plan$term == "necessary_n_clusters"]
  } else {
    K <- .dmar_check_count(specified_n_clusters, "specified_n_clusters")
  }
  K <- as.integer(K)
  if (K < .me_min_clusters(design$level)) .dmar_refuse_argument(sprintf(paste0(
    "The number of clusters is %d, and the interval for %s needs at least ",
    "%d (its t quantile has %s degrees of freedom)."), K,
    if (design$level == 2L) "a cluster-level predictor" else
      "an individual-level predictor",
    .me_min_clusters(design$level),
    if (design$level == 2L) "n_clusters - 2" else
      "n_clusters * (cluster_size - 1) - 1"),
    arg = "specified_n_clusters")
  m <- .dmar_check_count(cluster_size, "cluster_size", min = 2L)

  tau2   <- true_icc * true_total_var
  sigma2 <- (1 - true_icc) * true_total_var

  beta_hat <- numeric(G)
  ci_lo    <- numeric(G)
  ci_hi    <- numeric(G)
  ci_w     <- numeric(G)
  tI_lower <- logical(G)
  tI_upper <- logical(G)

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    one <- .me_one_interval(K, m, true_reg_coef, tau2, sigma2, true_var_X,
                            true_icc_X, design$level, design$centering,
                            conf_level)
    beta_hat[g] <- one[1L]
    ci_lo[g]    <- one[2L]
    ci_hi[g]    <- one[3L]
    ci_w[g]     <- one[3L] - one[2L]
    tI_lower[g] <- true_reg_coef < one[2L]
    tI_upper[g] <- true_reg_coef > one[3L]
  }

  if (!is.null(filename)) {
    per_rep <- data.frame(reg_coef = beta_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_w,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  expected <- .me_expected_width(K, m, true_total_var, true_icc, true_var_X,
                                 true_icc_X, design$level, design$centering,
                                 conf_level)
  na_if_null <- function(x) if (is.null(x)) NA_real_ else x
  out <- data.frame(
    term  = c("mean_reg_coef", "median_reg_coef", "sd_reg_coef",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "n_clusters", "cluster_size", "total_N", "predictor_level",
              "expected_width",
              "true_total_var", "true_var_X", "true_icc", "true_icc_X",
              "true_reg_coef",
              "supposed_total_var", "supposed_var_X", "supposed_icc",
              "supposed_icc_X",
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
              K, m, K * m, design$level,
              expected,
              true_total_var, true_var_X, true_icc, true_icc_X, true_reg_coef,
              na_if_null(supposed_total_var), na_if_null(supposed_var_X),
              na_if_null(supposed_icc),
              if (est_supplied) supposed_icc_X else NA_real_,
              width, conf_level),
    stringsAsFactors = FALSE, row.names = NULL
  )
  out <- .ss_sens_record(.as_dmar_tbl(out, conf_level = conf_level), G,
                         !is.na(ci_w), conf_level)
  attr(out, "centering") <- design$centering
  out
}
