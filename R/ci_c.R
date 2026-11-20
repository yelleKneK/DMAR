#' Confidence Interval for a Contrast in a Fixed Effects ANOVA
#'
#' Computes the confidence interval for an unstandardized contrast of means in
#' a fixed effects analysis of variance, so a focused comparison among groups
#' (a pairwise difference or any weighted combination of the means) is
#' reported with its precision and in the units of the response. Homogeneity
#' of variance is assumed, as in the ANOVA on which \code{s_anova} is based.
#'
#' @param means A vector of the group means or the means of the particular level of the effect (for fixed effect designs)
#' @param s_anova The standard deviation of the errors from the ANOVA model (i.e., the square root of the mean square error); required, as a single positive number, since the standard error of the contrast is formed from it and the means and sample sizes cannot supply it
#' @param c_weights The contrast weights (choose weights so that the positive \emph{c}-weights sum to 1 and the negative \emph{c}-weights sum to -1; i.e., use fractional values not integers)
#' @param n_per_group Sample sizes \emph{per group} or level of the particular factor (if length 1 it is assumed that the per group/level sample sizes are equal)
#' @param N Total sample size
#' @param contrast The contrast effect, obtained by multiplying the \emph{j}th mean by the \emph{j}th contrast weight.
#' @param conf_level Confidence interval coverage (i.e., 1- Type I error rate); default is .95. A level above .999999998,
#'   whose tail areas are below 1e-9, is refused
#' @param alpha_lower Type I error for the lower confidence limit, given
#'   together with \code{alpha_upper}; one without the other is refused.
#' @param alpha_upper Type I error for the upper confidence limit.
#'   The table records the coverage, \code{1 - alpha_lower - alpha_upper},
#'   as its level and, when the two tail areas differ, the areas
#'   themselves, which the line beneath the printed table,
#'   \code{\link{results_sentence}}, and the \code{alpha_lower} and
#'   \code{alpha_upper} columns of \code{tidy()} state; with a tail area
#'   of 0 the interval is a one-sided confidence bound. A nonzero tail
#'   area below 1e-9 is refused (see "Small tail areas" in
#'   \code{\link{ci_nc_t}}).
#' @param df_error The degrees of freedom for the error. In one-way designs it is \emph{N} minus the number of groups, \code{N - length(means)}, and need not be specified; it must be specified if the design has multiple factors.
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term} and \code{value}. The
#' \code{term} values are \code{"ci_lower"} (the lower confidence limit
#' on the population contrast), \code{"contrast"} (the estimated
#' unstandardized contrast), and \code{"ci_upper"} (the upper limit).
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
#'   ANCOVA and ANOVA contrasts: Sample size planning via narrow
#'   confidence intervals.
#'   \emph{British Journal of Mathematical and Statistical Psychology, 65},
#'   350--370. \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons of means.)
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the analysis of variance and contrast analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' Be sure to use the standard deviation and not the error variance for \code{s_anova}, not the square of this value (the error variance) which would come from the source table
#' (i.e., use the root mean square error, not the mean square error).
#'
#' Be sure to use fractional \emph{c}-weights when doing complex contrasts (not integers) to specify \code{c_weights}.
#' For example, in an ANCOVA of four groups, if the user wants to compare the mean of group 1 and 2 with the mean of group 3 and 4,
#' \code{c_weights} should be specified as c(0.5, 0.5, -0.5, -0.5) rather than c(1, 1, -1, -1). Make sure the sum of the contrast weights is zero.
#'
#' @seealso
#' \code{\link{ci_sc}}
#'
#' @examples
#' # Here is a four group example. Suppose that the means of groups 1--4 are 2, 4, 9,
#' # and 13, respectively. Further, let the error variance be .64 and thus the standard
#' # deviation would be .80 (note we use the standard deviation in the function, not the
#' # variance). The contrast of interest here is the average of groups 1 and 4 versus the
#' # average of groups 2 and 3.
#' ci_c(means = c(2, 4, 9, 13), s_anova = .80, c_weights = c(.5, -.5, -.5, .5),
#' n_per_group = c(3, 3, 3, 3), N = 12, conf_level = .95)
#'
#' # Here is an example with two groups.
#' ci_c(means = c(1.6, 0), s_anova = .80, c_weights = c(1, -1),
#' n_per_group = c(10, 10), N = 20, conf_level = .95)
#'
#' # An example given by Maxwell, Delaney, & Kelley (2027) :
#' # 20 subjects of mild hypertensives are assigned to one of four treatments: drug
#' # therapy, biofeedback, dietary modification, and a treatment combining all the
#' # three previous treatments. Subjects' blood pressure is measured two weeks
#' # after the termination of treatment. Now we want to form a 95% level
#' # confidence interval for the difference in blood pressure between subjects
#' # who received drug treatment and those who received biofeedback treatment
#'
#' ## Drug group's mean = 94; group size=4
#' ## Biofeedback group's mean = 91; group size=6
#' ## Diet group's mean = 92; group size=5
#' ## Combination group's mean = 83; group size=5
#' ## Mean Square Within (i.e., 'error_variance') = 67.375
#'
#' ci_c(means = c(94, 91, 92, 83), s_anova = sqrt(67.375), c_weights = c(1, -1, 0, 0),
#' n_per_group = c(4, 6, 5, 5), N = 20, conf_level = .95)
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_c <- function(means = NULL, s_anova = NULL, c_weights = NULL, n_per_group = NULL, N = NULL,
                 contrast = NULL, conf_level = .95, alpha_lower = NULL, alpha_upper = NULL, df_error = NULL) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  # The inputs are checked before either path is taken. The standard error
  # of the contrast is s_anova * sqrt(sum(c_weights^2 / n_per_group)), and s_anova,
  # the root mean square error (the square root of MS within), cannot be
  # recovered from the means and the sample sizes (Maxwell, Delaney, &
  # Kelley, 2027, ch. 4), so it is required on both paths and must be a
  # positive standard deviation. Through 1.0.1 it was checked only when contrast
  # was given: means without s_anova gave a zero-width interval at the
  # estimate (each limit added contrast to a product with NULL, which is
  # numeric(0), and data.frame() recycled contrast into all three rows), an
  # s_anova of 0 gave the same, and a negative one swapped the limits.
  if (is.null(c_weights)) .dmar_refuse_argument(
    "You must specify the vector of contrast weights ('c_weights').",
    arg = "c_weights")
  if (is.null(means) && is.null(contrast)) .dmar_refuse_argument(
    "Specify either the vector of means ('means') or the contrast effect ('contrast').",
    arg = c("means", "contrast"))
  if (!is.null(contrast) && !is.null(means)) .dmar_refuse_argument(
    "Since the contrast effect ('contrast') was specified, you should not specify the vector of means ('means').",
    arg = c("contrast", "means"))
  .check_error_sd(s_anova)

  if (!is.null(means) && length(means) != length(c_weights)) .dmar_refuse_argument(
    "The lengths of 'means' and 'c_weights' differ, which should not be the case.",
    arg = c("means", "c_weights"))

  # Recycle a scalar n_per_group by the number of contrast weights, not the number of
  # means: on the contrast-only path means is NULL, so length(means) is 0 and a
  # scalar n_per_group would recycle to length zero and fail the length check below.
  if (length(n_per_group) == 1) {
    n_per_group <- rep(n_per_group, length(c_weights))
  }

  if (length(n_per_group) != length(c_weights)) .dmar_refuse_argument(
    "The lengths of 'n_per_group' and 'c_weights' differ, which should not be the case.",
    arg = c("n_per_group", "c_weights"))

  # if(!identical(sum(c_weights[c_weights>0]), 1)) stop("Please use fractions to specify the contrast weights")

  if (!identical(round(sum(c_weights), 5), 0)) .dmar_refuse_argument(
    "The sum of the contrast weights ('c_weights') should equal zero.",
    arg = "c_weights")
  # All-zero weights define no contrast; they gave an interval of zero
  # width around the contrast ([0, 0] from the means, [.5, .5] with
  # contrast = .5).
  .dmar_check_contrast_nonzero(c_weights)

  part_of_se <- sqrt(sum((c_weights^2) / n_per_group))


  if (!is.null(contrast)) {
    if (is.null(n_per_group)) .dmar_refuse_argument(
      "You must specify the vector per group/level sample size ('n_per_group').",
      arg = "n_per_group")
  }


  if (!is.null(means)) {
    contrast <- sum(c_weights * means)
  }

  # A supplied pair of alphas sets the coverage to 1 - alpha_lower - alpha_upper,
  # so the default conf_level must not travel to the footer (it would mislabel a,
  # say, 90% interval as 95%). Drop conf_level when the alphas are given
  # directly; the table records the coverage they define instead.
  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (alphas_supplied) conf_level <- NULL
  # One tail area without the other stopped inside qt() or data.frame()
  # with R's own message ("Non-numeric argument to mathematical function"
  # or "arguments imply differing number of rows: 3, 2").
  if (alphas_supplied && (is.null(alpha_lower) || is.null(alpha_upper))) {
    given <- if (is.null(alpha_lower)) "alpha_upper" else "alpha_lower"
    .dmar_refuse_argument(c(
      "Only '", given, "' was given. 'alpha_lower' and 'alpha_upper' go ",
      "together: each is the tail area beyond one limit, so the interval ",
      "needs both. Give both (for example alpha_lower = .025 and ",
      "alpha_upper = .025, or 0 for a one-sided interval), or neither and ",
      "a 'conf_level'."),
      arg = c("alpha_lower", "alpha_upper"))
  }

  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }
  # The floor on a tail area (R/dmar_conditions.R), which ci_sc(), the
  # standardized interval, also keeps.
  .dmar_check_tail_floor(alpha_lower, alpha_upper, conf_level = conf_level)

  if (is.null(N) && is.null(n_per_group)) .dmar_refuse_argument(
    "You must specify the either total sample size ('N'), or sample sizes per group('n_per_group').",
    arg = c("N", "n_per_group"))
  if (is.null(N) && !is.null(n_per_group)) N <- sum(n_per_group)
  # One-way error df is N - (number of groups). Use length(c_weights), which is
  # always the group count, rather than length(means): on the contrast-only path
  # means is NULL, so length(means) is 0 and df would collapse to N, giving a
  # too-narrow (anti-conservative) interval.
  if (is.null(df_error)) df_2 <- N - length(c_weights)
  if (!is.null(df_error)) df_2 <- df_error

  # The central t quantiles. (Through 3338679 these passed ncp = 0, which
  # sends qt() to its noncentral search; the result is the same for
  # df_error >= 1, and below 1 the central quantile is exact.) The upper
  # one is taken in the upper tail (lower.tail = FALSE): 1 - alpha_upper
  # rounds a small tail area away, and DMAR 1.0.0 gave an upper limit of
  # Inf at alpha_upper = 1e-20.
  cv_upper <- qt(alpha_upper, df = df_2, lower.tail = FALSE)
  cv_lower <- qt(alpha_lower, df = df_2)

  term <- c("ci_lower", "contrast", "ci_upper")
  value <- c(contrast + cv_lower * part_of_se * s_anova, contrast, contrast + cv_upper * part_of_se * s_anova)

  # print(paste("The", 1 - (alpha_lower + alpha_upper), "confidence limits for the contrast are given as:"))
  return(.as_dmar_tbl(data.frame(term, value), conf_level = conf_level %||% (1 - alpha_lower - alpha_upper),
                      tail_areas = .tail_areas(alpha_lower, alpha_upper)))
}

# The standard deviation of the errors that a contrast interval scales the
# contrast by (s_anova of ci_c() and ci_sc(), s_ancova of ci_c_ancova(),
# ci_sc_ancova() and ci_c_ancova_bp(), and s_anova of ci_sc_ancova() when
# it is the standardizer): required, and a single positive, finite number,
# since it cannot be recovered from the means and sample sizes, a standard
# deviation of 0 gives a zero-width interval, and a negative one swaps the
# limits (or the sign of a standardized contrast). `arg` names the
# argument and `model` the model whose errors it describes. Not exported.
.check_error_sd <- function(value, arg = "s_anova", model = "ANOVA") {
  if (is.null(value))
    .dmar_refuse_argument(
      sprintf("You must specify the standard deviation of the errors of the %s model ('%s'), the square root of its mean square error.", model, arg),
      arg = arg)
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0)
    .dmar_refuse_argument(
      sprintf("'%s' must be a single positive, finite number: the standard deviation of the errors of the %s model, the square root of its mean square error (not the error variance).", arg, model),
      arg = arg)
  invisible(TRUE)
}
