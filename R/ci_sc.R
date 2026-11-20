#' Confidence Interval for a Standardized Contrast in a Fixed Effects ANOVA
#'
#' @description
#' Computes the exact noncentral \emph{t}-based confidence interval for
#' a standardized contrast of means in a fixed effects analysis of
#' variance: the contrast of interest divided by the error standard
#' deviation, so groups measured in raw units are compared on a common
#' standardized scale.
#'
#' @param means A vector of the group means or the means of the particular level of the effect (for fixed effect designs)
#' @param s_anova The standard deviation of the errors from the ANOVA model (i.e., the square root of the mean square error); required, as a single positive number, on every path
#' @param c_weights The contrast weights (choose weights so that the positive \emph{c}-weights sum to 1 and the negative \emph{c}-weights sum to -1; i.e., use fractional values such as \code{c(.5, .5, -1)}, not \code{c(1, 1, -2)}). The rule is applied to the values, so weights stored as integers, such as \code{-1:1}, are accepted.
#' @param n_per_group Sample sizes \emph{per group} or sample sizes for the level of the particular factor (if length 1 it is assumed that the sample size \emph{per group} or for the level of the particular factor are are equal)
#' @param N Total sample size
#' @param contrast The (unstandardized) contrast effect, obtained by multiplying the \emph{j}th mean by the \emph{j}th contrast weight (this is the unstandardized effect)
#' @param ncp The noncentrality parameter from the \emph{t}-distribution
#' @param conf_level Desired level of confidence for the computed interval (i.e., 1 - the Type I error rate). A level
#'   above .999999998, whose tail areas are below 1e-9, is refused
#' @param alpha_lower The Type I error rate for the lower confidence interval limit
#' @param alpha_upper The Type I error rate for the upper confidence interval limit.
#'   The table records the coverage, \code{1 - alpha_lower - alpha_upper},
#'   as its level and, when the two tail areas differ, the areas
#'   themselves, which the line beneath the printed table,
#'   \code{\link{results_sentence}}, and the \code{alpha_lower} and
#'   \code{alpha_upper} columns of \code{tidy()} state; with a tail area
#'   of 0 the interval is a one-sided confidence bound. A nonzero tail
#'   area below 1e-9 is refused (see "Small tail areas" in
#'   \code{\link{ci_nc_t}}).
#' @param df_error The degrees of freedom for the error. In one-way designs it is \emph{N} minus the number of groups, \code{N - length(means)}, and need not be specified; it must be specified if the design has multiple factors.
#' @param \dots Optional additional specifications for nested functions
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term} and \code{value}. The
#' \code{term} values are \code{"ci_lower"} (the lower confidence limit
#' on the standardized contrast), \code{"std_contrast"} (the standardized
#' contrast), and \code{"ci_upper"} (the upper limit).
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
#' ANCOVA and ANOVA contrasts: Sample size planning via narrow
#' confidence intervals.
#' \emph{British Journal of Mathematical and Statistical Psychology, 65},
#' 350--370. \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the
#' Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' Be sure to use the standard deviation and not the error variance for \code{s_anova},
#' not the square of this value (the error variance) which would come from the source table
#' (i.e., do not use the variance of the error but rather use its square root, the standard deviation).
#'
#' Be sure to use fractional \emph{c}-weights when doing complex contrasts (not integer values such as \code{c(1, 1, -2)}) to specify \code{c_weights}.
#' For example, in an ANOVA of four groups, if the user wants to compare the mean of group 1 and 2 with the mean of
#' group 3 and 4, \code{c_weights} should be specified as c(0.5, 0.5, -0.5, -0.5) rather than c(1, 1, -1, -1).
#' Make sure the sum of the contrast weights are zero.
#'
#' @seealso
#' \code{\link{ci_nc_t}}, \code{\link{ci_src}}, \code{\link{ci_smd}}, \code{\link{ci_smd_c}}, \code{\link{ci_sm}}, \code{\link{ci_c}} \code{\link{ci_c_ancova}}
#'
#' @examples
#' # Here is a four group example. Suppose that the means of groups 1--4 are 2, 4, 9,
#' # and 13, respectively. Further, let the error variance be .64 and thus the standard
#' # deviation would be .80 (note we use the standard deviation in the function, not the
#' # variance). The standardized contrast of interest here is the average of groups 1 and 4
#' # versus the average of groups 2 and 3.
#'
#' ci_sc(means = c(2, 4, 9, 13), s_anova = .80, c_weights = c(.5, -.5, -.5, .5),
#'       n_per_group = c(3, 3, 3, 3), N = 12, conf_level = .95)
#'
#' # Here is an example with two groups.
#' ci_sc(means = c(1.6, 0), s_anova = .80, c_weights = c(1, -1),
#'       n_per_group = c(10, 10), N = 20, conf_level = .95)
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_sc <- function(means = NULL, s_anova = NULL, c_weights = NULL, n_per_group = NULL, N = NULL,
                  contrast = NULL, ncp = NULL, conf_level = .95, alpha_lower = NULL, alpha_upper = NULL, df_error = NULL, ...) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  # The unstandardized contrast was 'psi' through DMAR 1.0.x. The function
  # passes `...` to ci_nc_t(), so that name would not stop with "unused
  # argument" here; it is refused instead.
  # The name is read from the call, so the arguments in `...` are not
  # evaluated here.
  if ("psi" %in% names(match.call(expand.dots = FALSE)$...))
    .dmar_refuse_argument(c(
      "'psi' is not an argument of ci_sc(): give the unstandardized contrast ",
      "as 'contrast'."),
      arg = "psi")
  .dmar_check_numbers("ci_sc", "N")
  # The inputs are checked before any path is taken: exactly one of means,
  # contrast, and ncp, and on every path a positive standard deviation of the
  # errors (.check_error_sd(), shared with ci_c()). Through 1.0.1 s_anova
  # was checked only on the contrast and ncp paths, so means without it failed
  # with "argument is of length zero", and an s_anova of 0 or below was
  # not refused.
  if (is.null(c_weights)) .dmar_refuse_argument(
    "You must specify the vector of contrast weights ('c_weights').",
    arg = "c_weights")
  if (is.null(means) && is.null(contrast) && is.null(ncp)) .dmar_refuse_argument(
    "Specify one of the vector of means ('means'), the contrast effect ('contrast'), or the noncentral parameter ('ncp').",
    arg = c("means", "contrast", "ncp"))
  .check_error_sd(s_anova)

  # The weights rule, judged by value (R/ss_aipe_sc_internals.R), which
  # ss_aipe_sc() and ss_aipe_sc_sensitivity() share. Through 1.0.x the
  # positive weights were compared with 1 by identical(), so integer
  # weights such as -1:1 were refused.
  .sc_check_c_weights(c_weights)

  if (!is.null(means) && length(means) != length(c_weights)) .dmar_refuse_argument(
    "The lengths of 'means' and 'c_weights' differ, which should not be the case.",
    arg = c("means", "c_weights"))

  # Recycle a scalar n_per_group by the number of contrast weights, not the number of
  # means: on the contrast-only and ncp-only paths means is NULL, so length(means)
  # is 0 and a scalar n_per_group would recycle to length zero and fail the length
  # check below.
  if (length(n_per_group) == 1) {
    n_per_group <- rep(n_per_group, length(c_weights))
  }

  if (length(n_per_group) != length(c_weights)) .dmar_refuse_argument(
    "The lengths of 'n_per_group' and 'c_weights' differ, which should not be the case.",
    arg = c("n_per_group", "c_weights"))


  part_of_se <- sqrt(sum((c_weights^2) / n_per_group))

  if (!is.null(contrast)) {
    if (!is.null(means)) .dmar_refuse_argument(
      "Since the contrast effect ('contrast') was specified, you should not specify the vector of means ('means').",
      arg = c("contrast", "means"))
    if (!is.null(ncp)) .dmar_refuse_argument(
      "Since the contrast effect ('contrast') was specified, you should not specify the noncentral parameter ('ncp').",
      arg = c("contrast", "ncp"))
    if (is.null(n_per_group)) .dmar_refuse_argument(
      "You must specify the vector per group/level sample size ('n_per_group').",
      arg = "n_per_group")
    if (is.null(c_weights)) .dmar_refuse_argument(
      "You must specify the vector of contrast weights ('c_weights').",
      arg = "c_weights")

    contrast <- contrast / s_anova
    lambda <- contrast / part_of_se
  }

  if (!is.null(ncp)) {
    if (!is.null(means)) .dmar_refuse_argument(
      "Since the noncentral parameter was specified directly, you should not specify the vector of means ('means').",
      arg = "means")
    if (!is.null(contrast)) .dmar_refuse_argument(
      "Since the noncentral parameter was specified directly, you should not specify the the contrast effect ('contrast').",
      arg = "contrast")
    if (is.null(n_per_group)) .dmar_refuse_argument(
      "You must specify the vector per group/level sample size ('n_per_group').",
      arg = "n_per_group")
    if (is.null(c_weights)) .dmar_refuse_argument(
      "You must specify the vector of contrast weights ('c_weights'.",
      arg = "c_weights")

    lambda <- ncp
    # The means and contrast branches define the noncentrality as the standardized
    # contrast divided by part_of_se, so a directly supplied noncentrality
    # parameter implies a standardized contrast of lambda * part_of_se.
    # Without this inversion the std_contrast row was left undefined and the
    # closing data.frame() errored (three terms against two values).
    contrast <- lambda * part_of_se
  }

  if (!is.null(means)) {
    contrast <- sum(c_weights * means)
    contrast <- contrast / s_anova

    lambda <- contrast / part_of_se
  }

  # A supplied pair of alphas sets the coverage to 1 - alpha_lower - alpha_upper,
  # so the default conf_level must not also be attached (the footer would then
  # mislabel a, say, 90% interval as 95%); the table records the coverage the
  # alphas define. Reject mixing an explicit conf_level with the alphas.
  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (alphas_supplied) {
    if (!missing(conf_level) && !is.null(conf_level)) {
      .dmar_refuse_argument(
        "Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'; you cannot mix them.",
        arg = c("conf_level", "alpha_lower", "alpha_upper"))
    }
    conf_level <- NULL
  }

  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }
  # The floor on a tail area (R/dmar_conditions.R), here so that a level
  # that leaves too small a tail is refused in its own name; ci_nc_t()
  # checks the tails it is given against the same floor. With one tail
  # given and the other not, ci_nc_t() refuses first, asking for both.
  if (!is.null(alpha_lower) && !is.null(alpha_upper))
    .dmar_check_tail_floor(alpha_lower, alpha_upper,
                           conf_level = if (!alphas_supplied) conf_level)

  # Mirror ci_c: when only the per-group sample sizes are given, the total
  # sample size is their sum. This keeps the two sibling functions, which share
  # the (means, s_anova, c_weights, n_per_group, N) interface, resolving N the same way.
  if (is.null(N) && !is.null(n_per_group)) N <- sum(n_per_group)
  if (is.null(N)) .dmar_refuse_argument(
    "You must specify the total sample size ('N').", arg = "N")
  # One-way error df is N - (number of groups). Use length(c_weights), which is
  # always supplied and equals the group count, rather than length(means): on
  # the contrast-only and ncp-only paths means is NULL, so length(means) is 0 and df
  # would collapse to N, giving too-narrow noncentral t limits.
  if (is.null(df_error)) df_2 <- N - length(c_weights)
  if (!is.null(df_error)) df_2 <- df_error

  Lims <- ci_nc_t(
    ncp = lambda, df = df_2, conf_level = NULL, alpha_lower = alpha_lower,
    alpha_upper = alpha_upper, ...
  )

  term <- c("ci_lower", "std_contrast", "ci_upper")
  value <- c(Lims[which(Lims$term == "ci_lower"), 2] * part_of_se, contrast, Lims[which(Lims$term == "ci_upper"), 2] * part_of_se)

  # print(paste("The", 1 - (alpha_lower + alpha_upper), "confidence limits for the standardized contrast are given as:"))

  return(.as_dmar_tbl(data.frame(term, value), conf_level = conf_level %||% (1 - alpha_lower - alpha_upper),
                      tail_areas = .tail_areas(alpha_lower, alpha_upper)))
}
