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
#' @param s_anova The standard deviation of the errors from the ANOVA model (i.e., the square root of the mean square error)
#' @param c_weights The contrast weights (chose weights so that the positive \emph{c}-weights sum to 1 and the negative \emph{c}-weights sum to -1; i.e., use fractional values not integers).
#' @param n Sample sizes \emph{per group} or sample sizes for the level of the particular factor (if length 1 it is assumed that the sample size \emph{per group} or for the level of the particular factor are are equal)
#' @param N Total sample size
#' @param psi The (unstandardized) contrast effect, obtained by multiplying the \emph{j}th mean by the \emph{j}th contrast weight (this is the unstandardized effect)
#' @param ncp The noncentrality parameter from the \emph{t}-distribution
#' @param conf_level Desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)
#' @param alpha_lower The Type I error rate for the lower confidence interval limit
#' @param alpha_upper The Type I error rate for the upper confidence interval limit
#' @param df_error The degrees of freedom for the error. In one-way designs, this is simply \emph{N}-length (means) and need not be specified; it must be specified if the design has multiple factors.
#' @param \dots Optional additional specifications for nested functions
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term} and \code{value}. The
#' \code{term} values are \code{"lower_limit"} (the lower confidence limit
#' on the standardized contrast), \code{"std_contrast"} (the standardized
#' contrast), and \code{"upper_limit"} (the upper limit).
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
#' Be sure to use fractional \emph{c}-weights when doing complex contrasts (not integers) to specify \code{c_weights}.
#' For example, in an ANCOVA of four groups, if the user wants to compare the mean of group 1 and 2 with the mean of
#' group 3 and 4, \code{c_weights} should be specified as c(0.5, 0.5, -0.5, -0.5) rather than c(1, 1, -1, -1).
#' Make sure the sum of the contrast weights are zero.
#'
#' @seealso
#' \code{\link{conf_limits_nct}}, \code{\link{ci_src}}, \code{\link{ci_smd}}, \code{\link{ci_smd_c}}, \code{\link{ci_sm}}, \code{\link{ci_c}} \code{\link{ci_c_ancova}}
#'
#' @examples
#' # Here is a four group example. Suppose that the means of groups 1--4 are 2, 4, 9,
#' # and 13, respectively. Further, let the error variance be .64 and thus the standard
#' # deviation would be .80 (note we use the standard deviation in the function, not the
#' # variance). The standardized contrast of interest here is the average of groups 1 and 4
#' # versus the average of groups 2 and 3.
#'
#' ci_sc(means = c(2, 4, 9, 13), s_anova = .80, c_weights = c(.5, -.5, -.5, .5),
#'       n = c(3, 3, 3, 3), N = 12, conf_level = .95)
#'
#' # Here is an example with two groups.
#' ci_sc(means = c(1.6, 0), s_anova = .80, c_weights = c(1, -1),
#'       n = c(10, 10), N = 20, conf_level = .95)
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_sc <- function(means = NULL, s_anova = NULL, c_weights = NULL, n = NULL, N = NULL,
                  psi = NULL, ncp = NULL, conf_level = .95, alpha_lower = NULL, alpha_upper = NULL, df_error = NULL, ...) {
  if (!identical(sum(c_weights[c_weights > 0]), 1)) stop("Please use fractions to specify the contrast weights")
  if (!identical(round(sum(c_weights), 5), 0)) stop("The sum of the contrast weights ('c_weights') should equal zero.")

  if (!is.null(means) && length(means) != length(c_weights)) stop("The lengths of 'means' and 'c_weights' differ, which should not be the case.")

  # Recycle a scalar n by the number of contrast weights, not the number of
  # means: on the psi-only and ncp-only paths means is NULL, so length(means)
  # is 0 and a scalar n would recycle to length zero and fail the length
  # check below.
  if (length(n) == 1) {
    n <- rep(n, length(c_weights))
  }

  if (length(n) != length(c_weights)) stop("The lengths of 'n' and 'c_weights' differ, which should not be the case.")


  part_of_se <- sqrt(sum((c_weights^2) / n))

  if (!is.null(psi)) {
    if (!is.null(means)) stop("Since the contrast effect ('psi') was specified, you should not specify the vector of means ('means').")
    if (!is.null(ncp)) stop("Since the contrast effect ('psi') was specified, you should not specify the noncentral parameter ('ncp').")
    if (is.null(s_anova)) stop("You must specify the standard deviation of the errors (i.e., the square root of the error variance).")
    if (is.null(n)) stop("You must specify the vector per group/level sample size ('n').")
    if (is.null(c_weights)) stop("You must specify the vector of contrast weights ('c_weights').")

    psi <- psi / s_anova
    lambda <- psi / part_of_se
  }

  if (!is.null(ncp)) {
    if (!is.null(means)) stop("Since the noncentral parameter was specified directly, you should not specify the vector of means ('means').")
    if (!is.null(psi)) stop("Since the noncentral parameter was specified directly, you should not specify the the contrast effect ('psi').")
    if (is.null(s_anova)) stop("You must specify the standard deviation of the errors (i.e., the square root of the error variance).")
    if (is.null(n)) stop("You must specify the vector per group/level sample size ('n').")
    if (is.null(c_weights)) stop("You must specify the vector of contrast weights ('c_weights'.")

    lambda <- ncp
    # The means and psi branches define the noncentrality as the standardized
    # contrast divided by part_of_se, so a directly supplied noncentrality
    # parameter implies a standardized contrast of lambda * part_of_se.
    # Without this inversion the std_contrast row was left undefined and the
    # closing data.frame() errored (three terms against two values).
    psi <- lambda * part_of_se
  }

  if (!is.null(means)) {
    psi <- sum(c_weights * means)
    psi <- psi / s_anova

    lambda <- psi / part_of_se
  }

  # A supplied pair of alphas sets the coverage to 1 - alpha_lower - alpha_upper,
  # so the default conf_level must not also be attached (the footer would then
  # mislabel a, say, 90% interval as 95%). Reject mixing an explicit conf_level
  # with the alphas.
  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (alphas_supplied) {
    if (!missing(conf_level) && !is.null(conf_level)) {
      stop("Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'; you cannot mix them.", call. = FALSE)
    }
    conf_level <- NULL
  }

  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }

  # Mirror ci_c: when only the per-group sample sizes are given, the total
  # sample size is their sum. This keeps the two sibling functions, which share
  # the (means, s_anova, c_weights, n, N) interface, resolving N the same way.
  if (is.null(N) && !is.null(n)) N <- sum(n)
  if (is.null(N)) stop("You must specify the total sample size ('N').")
  # One-way error df is N - (number of groups). Use length(c_weights), which is
  # always supplied and equals the group count, rather than length(means): on
  # the psi-only and ncp-only paths means is NULL, so length(means) is 0 and df
  # would collapse to N, giving too-narrow noncentral t limits.
  if (is.null(df_error)) df_2 <- N - length(c_weights)
  if (!is.null(df_error)) df_2 <- df_error

  Lims <- conf_limits_nct(
    ncp = lambda, df = df_2, conf_level = NULL, alpha_lower = alpha_lower,
    alpha_upper = alpha_upper, ...
  )

  term <- c("lower_limit", "std_contrast", "upper_limit")
  value <- c(Lims[which(Lims$term == "lower_limit"), 2] * part_of_se, psi, Lims[which(Lims$term == "upper_limit"), 2] * part_of_se)

  # print(paste("The", 1 - (alpha_lower + alpha_upper), "confidence limits for the standardized contrast are given as:"))

  return(.as_dmar_tbl(data.frame(term, value), conf_level = conf_level))
}
