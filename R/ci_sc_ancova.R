#' Confidence Interval for a Standardized Contrast in ANCOVA With One Covariate
#'
#' @description
#' Calculate the confidence interval for a standardized contrast in ANCOVA with one covariate. The standardizer
#' (i.e., the divisor) can be either the error standard deviation of the ANOVA model (i.e., the model excluding the covariate)
#' or of the ANCOVA model.
#'
#' @param contrast Unstandardized contrast of adjusted means
#' @param adj_means The vector that contains the adjusted mean of each group on the dependent variable
#' @param s_anova The standard deviation of the errors from the ANOVA model (i.e., the square root of the mean square error from ANOVA); required, as a single positive number, when \code{standardizer = "s_anova"}
#' @param s_ancova The standard deviation of the errors from the ANCOVA model (i.e., the square root of the mean square error from ANCOVA); required, as a single positive number, with either standardizer
#' @param standardizer Which error standard deviation the user wants to use, the value of which can be either \code{"s_ancova"} or \code{"s_anova"}
#' @param c_weights The contrast weights (chose weights so that the positive \emph{c}-weights sum to 1 and the negative \emph{c}-weights sum to -1; i.e., use fractional values not integers).
#' @param n_per_group Either a single number that indicates the sample size per group, or a vector that contains the sample size of each group
#' @param cov_means A vector that contains the group means of the covariate
#' @param SSwithin_x The sum of squares within groups obtained from the summary table for ANOVA on the covariate
#' @param conf_level The desired confidence interval coverage, (i.e., 1 - Type I error rate). A level above
#'   .999999998, whose tail areas are below 1e-9, is refused (see "Small tail areas" in \code{\link{ci_nc_t}})
#'
#' @details
#' The argument \code{SSwithin_x} is the sum of squares within groups for the
#' covariate, taken from the ANOVA source table in which the covariate (not the
#' outcome) is the dependent variable. Published reports do not always print
#' this quantity directly. When a report gives the covariate group means, the
#' group sample sizes, and the \emph{F} statistic from the one-way ANOVA on the
#' covariate, \code{SSwithin_x} can be recovered algebraically. The worked
#' example below follows Lai and Kelley (2012): three groups of sizes 19, 18,
#' and 19 (so \eqn{N = 56}) have covariate means 60.08, 57.08, and 57.97, and
#' the covariate ANOVA reports \eqn{F = 0.756} with 2 and 53 degrees of freedom. The sum of squares between
#' groups for the covariate, computed from the group means and sample sizes, is
#' approximately 88.5, so the mean square between groups is approximately
#' \eqn{88.5 / 2 = 44.3}. Because \emph{F} is the ratio of the mean square
#' between groups to the mean square within groups, the mean square within
#' groups is approximately \eqn{44.3 / 0.756 = 58.6}, and the sum of squares
#' within groups is that mean square times its degrees of freedom, approximately
#' \eqn{58.6 \times 53 = 3103}. That recovered value is what you would pass to
#' \code{SSwithin_x}. The \dQuote{Examples} section reproduces this computation
#' in code.
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term} and \code{value}
#' (numeric). The \code{term} values are \code{"ci_lower"} (the lower
#' confidence limit on the standardized ANCOVA contrast), \code{"std_contrast"}
#' (the standardized contrast), and \code{"ci_upper"} (the upper
#' limit). The divisor used in standardization (either \code{"s_anova"}
#' or \code{"s_ancova"}) is attached as the \code{"standardizer"}
#' attribute of the returned data.frame.
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized
#'   effect sizes: Theory, application, and implementation.
#'   \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals.
#'   \emph{Psychological Methods, 11}, 363--385.
#'   \doi{10.1037/1082-989X.11.4.363}
#'
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
#'   ANCOVA and ANOVA contrasts: Sample size planning via narrow
#'   confidence intervals.
#'   \emph{British Journal of Mathematical and Statistical Psychology, 65},
#'   350--370. \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9.)
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical models. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' Be sure to use the standard deviations and not the error variances for \code{s_anova} and \code{s_ancova},
#' not the squares of these values which would come from the source tables
#' (i.e., do not use the variance of the errors but rather use its square root, the standard deviation).
#'
#' If \code{n_per_group} receives a single number, that number is considered as the sample size per group.
#' If \code{n_per_group} is assigned to a vector, the vector is considered as the sample size of each group.
#'
#' Be sure to use fractional \emph{c}-weights when doing complex contrasts (not integers) to specify
#' \code{c_weights}. For example, in an ANCOVA of four groups, if the user wants to compare the mean of
#' group 1 and 2 with the mean of group 3 and 4, \code{c_weights} should be specified as c(0.5, 0.5, -0.5, -0.5)
#' rather than c(1, 1, -1, -1). Make sure the sum of the contrast weights are zero.
#'
#' The argument to be assigned to \code{standardizer} must be either \code{"s_ancova"} or \code{"s_anova"}.
#'
#' @seealso \code{\link{ci_c_ancova}}, \code{\link{ci_sc}}
#'
#' @examples
#' # Maxwell, Delaney, & Kelley (2027) offer an example that 30 depressive
#' # individuals are randomly assigned to three groups, 10 in each, and ANCOVA
#' # is performed on the posttest scores using the participants' pretest
#' # scores as the covariate. The means of pretest scores of group 1, 2, and 3 are
#' # 17, 17.7, and 17.4, respectively, whereas the adjusted means of groups 1, 2, and 3
#' # are 7.5, 12, and 14, respectively. The error variance in ANCOVA is 29 and thus
#' # 5.385165 is the error standard deviation, with the sum of squares within groups
#' # from an ANOVA on the covariate is 752.5.
#'
#' # To obtained the confidence interval for the standardized adjusted mean difference
#' # between group 1 and 2, using the ANCOVA error standard deviation:
#' ci_sc_ancova(adj_means = c(7.5, 12, 14), s_ancova = 5.385165, c_weights = c(1, -1, 0),
#'              n_per_group = 10, cov_means = c(17, 17.7, 17.4), SSwithin_x = 752.5)
#'
#' # Or, with less error in rounding:
#' ci_sc_ancova(adj_means = c(7.54, 11.98, 13.98), s_ancova = 5.393, c_weights = c(-1, 0, 1),
#'              n_per_group = 10, cov_means = c(17, 17.7, 17.4), SSwithin_x = 752.5)
#'
#' # Now, using the standard deviation from ANOVA (and not ANCOVA as above), we have:
#' ci_sc_ancova(adj_means = c(7.54, 11.98, 13.98), s_anova = 6.294, s_ancova = 5.393,
#'              c_weights = c(-1, 0, 1),n_per_group = 10, cov_means = c(17, 17.7, 17.4),
#'              SSwithin_x = 752.5, standardizer = "s_anova", conf_level = .95)
#'
#' # Recovering SSwithin_x from a covariate ANOVA F when a report does not print
#' # it directly (see the Details section). This example follows Lai and Kelley
#' # (2012): three groups of sizes 19, 18, and 19 have covariate means 60.08,
#' # 57.08, and 57.97, and the one-way ANOVA on the covariate reports F = 0.756.
#' cov_means_ex  <- c(60.08, 57.08, 57.97)
#' n_ex          <- c(19, 18, 19)
#' grand_x       <- sum(n_ex * cov_means_ex) / sum(n_ex)
#' ss_between_x  <- sum(n_ex * (cov_means_ex - grand_x)^2)
#' ms_between_x  <- ss_between_x / (length(n_ex) - 1)
#' ms_within_x   <- ms_between_x / 0.756
#' SSwithin_x_ex <- ms_within_x * (sum(n_ex) - length(n_ex))
#'
#' ci_sc_ancova(adj_means = c(63.88, 62.39, 56.48), s_ancova = 20.267,
#'              c_weights = c(0.5, 0.5, -1), n_per_group = n_ex, cov_means = cov_means_ex,
#'              SSwithin_x = SSwithin_x_ex)
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_sc_ancova <- function(contrast = NULL, adj_means = NULL, s_anova = NULL, s_ancova = NULL, standardizer = "s_ancova", c_weights, n_per_group, cov_means, SSwithin_x, conf_level = .95) {
  .dmar_check_conf_level(conf_level)
  .dmar_check_numbers("ci_sc_ancova", "SSwithin_x", required = "SSwithin_x")
  if (!is.character(standardizer) || length(standardizer) != 1L ||
      !standardizer %in% c("s_ancova", "s_anova"))
    .dmar_refuse_argument(
      "The standardizer must be either 's_anova' or 's_ancova'.",
      arg = c("s_anova", "s_ancova"))
  # The noncentrality parameter is the contrast over s_ancova times
  # sqrt(sum(c_weights^2 / n_per_group) + (sum(c_weights * cov_means))^2 /
  # SSwithin_x) whichever standardizer is chosen, so s_ancova is required
  # with either; s_anova is required when it is the standardizer (the
  # interval is then rescaled by s_ancova / s_anova). Each must be a single
  # positive, finite number (.check_error_sd(), shared with ci_c()).
  # Through 1.0.1 s_ancova was checked only when s_anova was the
  # standardizer, so the default call without it failed with "argument is
  # of length zero", and neither was refused at 0 or below.
  .check_error_sd(s_ancova, "s_ancova", "ANCOVA")
  if (standardizer == "s_anova") .check_error_sd(s_anova, "s_anova", "ANOVA")

  if (missing(cov_means)) .dmar_refuse_argument(
    "The mean of the covariate within each group (i.e., the vector 'cov_means') is missing.",
    arg = "cov_means")
  if (is.null(contrast) && is.null(adj_means)) .dmar_refuse_argument(
    "Supply one of 'contrast' and 'adj_means'; neither was given.", arg = c("contrast", "adj_means"))
  if (!is.null(contrast) && !is.null(adj_means)) .dmar_refuse_argument(
    "Supply one of 'contrast' and 'adj_means', not both.", arg = c("contrast", "adj_means"))

  if (missing(SSwithin_x)) .dmar_refuse_argument(
    "The sum of squares within from the ANOVA on the covariate is missing (i.e., 'SSwithin_x').",
    arg = "SSwithin_x")

  if (is.null(contrast)) contrast <- sum(adj_means * c_weights)

  if (abs(sum(c_weights)) > 1e-8) .dmar_refuse_argument(
    "'c_weights' must sum to zero.", arg = "c_weights")
  # All-zero weights pass both rules and define no contrast; they stopped
  # with "missing value where TRUE/FALSE needed" (ss_aipe_sc_ancova()) or a
  # message about 'ncp' (ci_sc_ancova()).
  .dmar_check_contrast_nonzero(c_weights)
  if (sum(c_weights[c_weights > 0]) > 1) .dmar_refuse_argument(
    "The positive values of 'c_weights' must not sum to more than 1; give the weights as fractions, such as c(1/2, 1/2, -1).", arg = "c_weights")

  J <- length(c_weights)
  if (length(n_per_group) > 1) {
    if (length(n_per_group) != J) .dmar_refuse_argument(
      "'c_weights' and 'n_per_group' imply different numbers of groups.",
      arg = c("c_weights", "n_per_group"))
  }
  if (length(n_per_group) == 1) n_per_group <- rep(n_per_group, J)

  if (length(cov_means) != J) .dmar_refuse_argument(
    "'c_weights' and 'cov_means' imply different numbers of groups.",
    arg = c("c_weights", "cov_means"))

  f_x_numerator <- (sum(c_weights * cov_means))^2
  f_x_denominator <- SSwithin_x
  sample_size_weighted <- sum(c_weights^2 / n_per_group)
  ratio <- s_ancova / s_anova
  alpha <- 1 - conf_level
  nu <- sum(n_per_group) - J - 1
  # The floor on a tail area (R/dmar_conditions.R): the interval has equal
  # tails, each half of 1 - conf_level.
  .dmar_check_tail_floor(alpha / 2, alpha / 2, conf_level = conf_level)

  if (standardizer == "s_ancova") {
    contrast <- contrast / s_ancova
    lambda_obs <- contrast / sqrt(sample_size_weighted + (f_x_numerator / f_x_denominator))
    lambda_limits <- ci_nc_t(ncp = lambda_obs, df = nu, conf_level = 1 - alpha)

    psi_limit_upper <- lambda_limits[which(lambda_limits$term == "ci_upper"), 2] * sqrt(sample_size_weighted + (f_x_numerator / f_x_denominator))
    psi_limit_lower <- lambda_limits[which(lambda_limits$term == "ci_lower"), 2] * sqrt(sample_size_weighted + (f_x_numerator / f_x_denominator))
  }

  if (standardizer == "s_anova") {
    contrast <- contrast / s_anova
    lambda_obs <- contrast / (ratio * sqrt(sample_size_weighted + (f_x_numerator / f_x_denominator)))
    lambda_limits <- ci_nc_t(ncp = lambda_obs, df = nu, conf_level = 1 - alpha)

    psi_limit_upper <- lambda_limits[which(lambda_limits$term == "ci_upper"), 2] * ratio * sqrt(sample_size_weighted + (f_x_numerator / f_x_denominator))
    psi_limit_lower <- lambda_limits[which(lambda_limits$term == "ci_lower"), 2] * ratio * sqrt(sample_size_weighted + (f_x_numerator / f_x_denominator))
  }

  out <- data.frame(
    term  = c("ci_lower", "std_contrast", "ci_upper"),
    value = c(psi_limit_lower, contrast, psi_limit_upper)
  )
  attr(out, "standardizer") <- standardizer
  .as_dmar_tbl(out, conf_level = conf_level)
}
