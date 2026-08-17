#' Confidence Interval for a Contrast in a Fixed Effects ANOVA
#'
#' Computes the confidence interval for an unstandardized contrast of means in
#' a fixed effects analysis of variance, so a focused comparison among groups
#' (a pairwise difference or any weighted combination of the means) is
#' reported with its precision and in the units of the response. Homogeneity
#' of variance is assumed, as in the ANOVA on which \code{s_anova} is based.
#'
#' @param means A vector of the group means or the means of the particular level of the effect (for fixed effect designs)
#' @param s_anova The standard deviation of the errors from the ANOVA model (i.e., the square root of the mean square error)
#' @param c_weights The contrast weights (choose weights so that the positive \emph{c}-weights sum to 1 and the negative \emph{c}-weights sum to -1; i.e., use fractional values not integers)
#' @param n Sample sizes \emph{per group} or level of the particular factor (if length 1 it is assumed that the per group/level sample sizes are equal)
#' @param N Total sample size
#' @param psi The contrast effect, obtained by multiplying the \emph{j}th mean by the \emph{j}th contrast weight.
#' @param conf_level Confidence interval coverage (i.e., 1- Type I error rate); default is .95
#' @param alpha_lower Type I error for the lower confidence limit
#' @param alpha_upper Type I error for the upper confidence limit
#' @param df_error The degrees of freedom for the error. In one-way designs, this is simply \emph{N}-length (means) and need not be specified; it must be specified if the design has multiple factors.
#' @param \dots Allows one to potentially include parameter values for inner functions
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term} and \code{value}. The
#' \code{term} values are \code{"lower_limit"} (the lower confidence limit
#' on the population contrast), \code{"contrast"} (the estimated
#' unstandardized contrast), and \code{"upper_limit"} (the upper limit).
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
#' n = c(3, 3, 3, 3), N = 12, conf_level = .95)
#'
#' # Here is an example with two groups.
#' ci_c(means = c(1.6, 0), s_anova = .80, c_weights = c(1, -1),
#' n = c(10, 10), N = 20, conf_level = .95)
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
#' n = c(4, 6, 5, 5), N = 20, conf_level = .95)
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_c <- function(means = NULL, s_anova = NULL, c_weights = NULL, n = NULL, N = NULL,
                 psi = NULL, conf_level = .95, alpha_lower = NULL, alpha_upper = NULL, df_error = NULL, ...) {
  if (!is.null(means) && length(means) != length(c_weights)) stop("The lengths of 'means' and 'c_weights' differ, which should not be the case.")

  # Recycle a scalar n by the number of contrast weights, not the number of
  # means: on the psi-only path means is NULL, so length(means) is 0 and a
  # scalar n would recycle to length zero and fail the length check below.
  if (length(n) == 1) {
    n <- rep(n, length(c_weights))
  }

  if (length(n) != length(c_weights)) stop("The lengths of 'n' and 'c_weights' differ, which should not be the case.")

  # if(!identical(sum(c_weights[c_weights>0]), 1)) stop("Please use fractions to specify the contrast weights")

  if (!identical(round(sum(c_weights), 5), 0)) stop("The sum of the contrast weights ('c_weights') should equal zero.")

  part_of_se <- sqrt(sum((c_weights^2) / n))


  if (!is.null(psi)) {
    if (!is.null(means)) stop("Since the contrast effect ('psi') was specified, you should not specify the vector of means ('means').")
    if (is.null(s_anova)) stop("You must specify the standard deviation of the errors (i.e., the square root of the error variance).")
    if (is.null(n)) stop("You must specify the vector per group/level sample size ('n').")
    if (is.null(c_weights)) stop("You must specify the vector of contrast weights ('c_weights').")
  }


  if (!is.null(means)) {
    psi <- sum(c_weights * means)
  }

  # A supplied pair of alphas sets the coverage to 1 - alpha_lower - alpha_upper,
  # so the default conf_level must not travel to the footer (it would mislabel a,
  # say, 90% interval as 95%). Drop conf_level when the alphas are given directly.
  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (alphas_supplied) conf_level <- NULL

  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }

  if (is.null(N) && is.null(n)) stop("You must specify the either total sample size ('N'), or sample sizes per group('n').")
  if (is.null(N) && !is.null(n)) N <- sum(n)
  # One-way error df is N - (number of groups). Use length(c_weights), which is
  # always the group count, rather than length(means): on the psi-only path
  # means is NULL, so length(means) is 0 and df would collapse to N, giving a
  # too-narrow (anti-conservative) interval.
  if (is.null(df_error)) df_2 <- N - length(c_weights)
  if (!is.null(df_error)) df_2 <- df_error

  cv_upper <- qt(1 - alpha_upper, df = df_2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  cv_lower <- qt(alpha_lower, df = df_2, ncp = 0, lower.tail = TRUE, log.p = FALSE)

  term <- c("lower_limit", "contrast", "upper_limit")
  value <- c(psi + cv_lower * part_of_se * s_anova, psi, psi + cv_upper * part_of_se * s_anova)

  # print(paste("The", 1 - (alpha_lower + alpha_upper), "confidence limits for the contrast are given as:"))
  return(.as_dmar_tbl(data.frame(term, value), conf_level = conf_level))
}
