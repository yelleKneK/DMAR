#' @title Sensitivity Analysis for Sample Size Planning for the Standardized ANOVA Contrast From
#' the Accuracy in Parameter Estimation (AIPE) Perspective
#'
#' @description Performs a sensitivity analysis when planning sample size from the Accuracy in Parameter Estimation (AIPE)
#' Perspective for the standardized ANOVA contrast.
#'
#' @param true_std_contrast population standardized contrast
#' @param supposed_std_contrast supposed standardized contrast, the planning value given to \code{\link{ss_aipe_sc}} to choose the sample size
#' @param c_weights the contrast weights, fractions whose positive values sum to 1 and whose sum is zero, the rule of \code{\link{ci_sc}} applied to the values as \code{\link{ss_aipe_sc}} applies it: weights \code{\link{ss_aipe_sc}} refuses are refused here, with a condition of class \code{dmar_argument_invalid}, and integer weights such as \code{-1:1} are accepted
#' @param width the desired full width of the obtained confidence interval
#' @param specified_n_per_group selected sample size to use in order to determine distributional properties of at a given value of sample size
#' @param assurance parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity)
#' @param conf_level the desired confidence interval coverage, (i.e., 1 - Type I error rate)
#' @param G number of generations (i.e., replications) of the simulation
#' @param print_iter to print the current value of the iterations
#' @param filename an optional path for a comma separated file recording
#'   every replication (the realized standardized contrast,
#'   \code{std_contrast}, the full and
#'   the two one-sided interval widths, the three non-coverage indicators,
#'   and the two confidence limits): nothing is written when
#'   \code{filename} is \code{NULL} (the default), a new file with a header
#'   row is created otherwise, an existing file at that path is appended
#'   to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}
#' @param assurance_method passed to \code{\link{ss_aipe_sc}} when the
#'   sample size is planned from \code{supposed_std_contrast} with
#'   \code{assurance}. Only \code{"analytic"} (the default) is offered: the
#'   plan is the smallest sample size per group, at or above the plan
#'   without \code{assurance}, at which the exact probability of an
#'   interval no wider than \code{width} is at least \code{assurance}.
#'   \code{"simulation"} is refused with a condition of class
#'   \code{dmar_argument_invalid}, as it is by \code{\link{ss_aipe_sc}},
#'   since that probability is exact, whether or not \code{assurance} or
#'   \code{specified_n_per_group} is given.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @details
#' Each replication draws \emph{J} groups of \emph{n} normal observations
#' with standard deviation 1 and forms the \code{\link{ci_sc}} interval
#' from the group means and the pooled standard deviation, the square root
#' of the mean square error that \code{\link{ci_sc}} takes as
#' \code{s_anova}.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}
#' summarizing the Monte Carlo sensitivity analysis across the \code{G}
#' replications. The \code{term} entries are: \code{mean_std_contrast},
#' \code{median_std_contrast}, \code{sd_std_contrast} (summaries of the realized
#' standardized contrast); \code{mean_ci_width}, \code{median_ci_width},
#' \code{sd_ci_width} (summaries of the full interval widths);
#' \code{mean_ci_width_lower} and \code{mean_ci_width_upper} (mean
#' one-sided widths, measured from the observed contrast to each limit);
#' \code{pct_ci_less_width} (proportion of intervals at or below the target
#' width); \code{pct_ci_miss_low} and \code{pct_ci_miss_high}
#' (tail-specific empirical non-coverage of \code{true_std_contrast});
#' \code{total_type_I_error} (overall empirical non-coverage, the sum of
#' the two tails); and the input echoes \code{n_per_group} (the size simulated at, planned or specified),
#' \code{total_N}, \code{true_std_contrast}, \code{supposed_std_contrast} (NA when
#' \code{specified_n_per_group} was supplied instead), \code{width},
#' \code{conf_level}, and \code{assurance} (present only when an
#' assurance was supplied). The proportion and Type I error rows are proportions
#' on the 0 to 1 scale, not percentages.
#'
#' @references
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#'   calculation of confidence intervals that are based on central and
#'   noncentral distributions. \emph{Educational and Psychological
#'   Measurement, 61}(4), 532--574. \doi{10.1177/0013164401614002}
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators.
#' \emph{Journal of Educational Statistics, 6}(2), 107--128.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application,
#' and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the standardized mean difference: Accuracy in Parameter
#' Estimation via narrow confidence intervals. \emph{Psychological Methods, 11}(4), 363--385.
#'   \doi{10.1037/1082-989X.11.4.363}
#'
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
#' ANCOVA and ANOVA contrasts: Sample size planning via narrow
#' confidence intervals.
#' \emph{British Journal of Mathematical and Statistical Psychology, 65},
#' 350--370. \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical models. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_sc}}, \code{\link{ss_aipe_c}}, \code{\link{ci_nc_t}}
#'
#' @examples
#' # Sensitivity analysis for a standardized three-group ANOVA contrast
#' # (-1, 0, 1) with a standardized contrast of 0.5 and target full width 0.40. G is kept small
#' # here so the example runs quickly; raise it for a stable sweep.
#' set.seed(113)
#' ss_aipe_sc_sensitivity(
#'   true_std_contrast = 0.5, supposed_std_contrast = 0.5,
#'   c_weights = c(-1, 0, 1),
#'   width = 0.40,
#'   conf_level = 0.95, G = 50, print_iter = FALSE
#' )
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export


ss_aipe_sc_sensitivity <- function(true_std_contrast = NULL, supposed_std_contrast = NULL, c_weights, width = NULL,
                                   specified_n_per_group = NULL, assurance = NULL, conf_level = .95, G = 10000,
                                   print_iter = TRUE, filename = NULL,
                                   assurance_method = "analytic", seed = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!is.null(width)) .dmar_check_width(width)
  .dmar_check_assurance_single(assurance)
  # The planner's method, refused here in this function's name when it is
  # not one ss_aipe_sc() offers (R/ss_aipe_assurance_internals.R).
  assurance_method <- .ss_aipe_assurance_method(
    assurance_method, "ss_aipe_sc_sensitivity", simulation = FALSE,
    both = FALSE)
  .dmar_check_numbers("ss_aipe_sc_sensitivity", c("true_std_contrast", "supposed_std_contrast"))
  if (is.null(supposed_std_contrast) && is.null(specified_n_per_group)) .dmar_refuse_argument(
    "You must specify either \'supposed_std_contrast\' or \'specified_n_per_group\' (i.e., the sample size per group ).",
    arg = c("supposed_std_contrast", "specified_n_per_group"))
  if (!is.null(supposed_std_contrast) && !is.null(specified_n_per_group)) .dmar_refuse_argument(
    "You must specify either \'supposed_std_contrast\' or \'specified_n_per_group\' (i.e., the per group sample size), but not both.",
    arg = c("supposed_std_contrast", "specified_n_per_group"))
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)

  # The weights rule of ci_sc() and ss_aipe_sc() (R/ss_aipe_sc_internals.R),
  # whether or not the sample size is planned here.
  .sc_check_c_weights(c_weights)

  if (!is.null(supposed_std_contrast)) {
    n <- ss_aipe_sc(
      true_std_contrast = supposed_std_contrast, conf_level = conf_level,
      c_weights = c_weights, width = width, assurance = assurance,
      assurance_method = assurance_method)[1, 2]
  }
  else {
    n <- specified_n_per_group
  }

  J <- length(c_weights)

  G <- G
  Full_Width <- rep(NA, G)
  Width_from_psi_obs_Lower <- rep(NA, G)
  Width_from_psi_obs_Upper <- rep(NA, G)
  Type_I_Error_Upper <- rep(NA, G)
  Type_I_Error_Lower <- rep(NA, G)
  Type_I_Error <- rep(NA, G)
  Lower_Limit <- rep(NA, G)
  Upper_Limit <- rep(NA, G)
  psi_obs <- rep(NA, G)
  psi_pop <- true_std_contrast

  # Create population group means for group 1 to group J
  signal <- TRUE
  mu_c <- rep(NA, J)
  for (cw in 1:J)
  {
    if (c_weights[cw] > 0 && signal) {
      mu_c[cw] <- psi_pop / c_weights[cw]
      signal <- FALSE
    }

    else {
      mu_c[cw] <- 0
    }
  }

  for (i in 1:G)
  {
    if (print_iter) cat(c(i), "\n")

    # generate random samples
    group_data <- array(NA, dim = c(n, J))
    x_bar <- rep(NA, J)
    sd_group <- rep(NA, J)

    for (p in 1:J)
    {
      group_data[, p] <- rnorm(n, mean = mu_c[p], sd = 1)
      x_bar[p] <- mean(group_data[, p])
      sd_group[p] <- sd(group_data[, p])
    }

    # The standardized contrast and its ci_sc() interval: the divisor is
    # the pooled standard deviation, the root mean square error ci_sc()
    # takes as s_anova (DMAR 1.0.x used mean(sd_group); see Details).
    s_pooled <- sqrt(sum((n - 1) * sd_group^2) / (n * J - J))
    psi_obs[i] <- sum(c_weights * x_bar) / s_pooled
    lambda <- psi_obs[i] / sqrt(sum(c_weights^2) / n)
    lambda_limits <- ci_nc_t(ncp = lambda, df = n * J - J, conf_level = conf_level)
    psi_limit_upper <- lambda_limits[2, 2] * sqrt(sum(c_weights^2) / n)
    psi_limit_lower <- lambda_limits[1, 2] * sqrt(sum(c_weights^2) / n)

    Full_Width[i] <- psi_limit_upper - psi_limit_lower
    Width_from_psi_obs_Lower[i] <- psi_obs[i] - psi_limit_lower
    Width_from_psi_obs_Upper[i] <- psi_limit_upper - psi_obs[i]

    Type_I_Error_Upper[i] <- psi_pop > psi_limit_upper
    Type_I_Error_Lower[i] <- psi_pop < psi_limit_lower
    Type_I_Error[i] <- Type_I_Error_Upper[i] | Type_I_Error_Lower[i]

    Lower_Limit[i] <- psi_limit_lower
    Upper_Limit[i] <- psi_limit_upper
  }

  Results <- data.frame(
    std_contrast = psi_obs, full_width = Full_Width, width_lower = Width_from_psi_obs_Lower,
    width_upper = Width_from_psi_obs_Upper, type_I_error_upper = Type_I_Error_Upper,
    type_I_error_lower = Type_I_Error_Lower, type_I_error = Type_I_Error, ci_lower = Lower_Limit,
    ci_upper = Upper_Limit
  )

  if (!is.null(filename)) {
    if (file.exists(filename) && file.size(filename) > 0) {
      message("The file '", filename, "' already exists; the simulation results are appended to it.")
      utils::write.table(Results, filename, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    } else {
      utils::write.table(Results, filename, sep = ",", row.names = FALSE, append = FALSE)
    }
  }

  Summary <- data.frame(
    term = c(
      "mean_std_contrast", "median_std_contrast", "sd_std_contrast",
      "mean_ci_width", "median_ci_width", "sd_ci_width",
      "mean_ci_width_lower", "mean_ci_width_upper",
      "pct_ci_less_width",
      "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
      "n_per_group", "total_N",
      "true_std_contrast", "supposed_std_contrast",
      "width", "conf_level",
      if (!is.null(assurance)) "assurance"
    ),
    value = c(
      mean(psi_obs), median(psi_obs), sd(psi_obs),
      mean(Full_Width), median(Full_Width), sd(Full_Width),
      mean(Width_from_psi_obs_Lower), mean(Width_from_psi_obs_Upper),
      mean(Full_Width <= width),
      mean(Type_I_Error_Lower), mean(Type_I_Error_Upper), mean(Type_I_Error),
      n, n * J,
      true_std_contrast,
      if (is.null(supposed_std_contrast)) NA_real_ else supposed_std_contrast,
      if (is.null(width)) NA_real_ else width,
      conf_level,
      if (!is.null(assurance)) assurance
    )
  )

  return(.ss_sens_record(.as_dmar_tbl(Summary, conf_level = conf_level), G,
                         !is.na(Full_Width), conf_level))
}
