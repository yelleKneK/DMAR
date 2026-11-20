#' Plan Sample Size to Make the Test of a Change in the Squared Multiple Correlation Coefficient Sufficiently Powerful
#'
#' @description
#' Determine the necessary sample size, or the statistical power at a
#' specified sample size, for the model comparison \emph{F}-test that a set of
#' predictors adds nothing to a regression model that already contains other
#' predictors. The restricted model has \code{p_restricted} predictors and the
#' full model adds one or more predictors to it, for \code{p_full} in all. The
#' null hypothesis is that the two models have the same population squared
#' multiple correlation coefficient, so that the change in \eqn{R^2} is zero.
#' This is the test of an increment in \eqn{R^2} in hierarchical regression
#' (Cohen, Cohen, West, & Aiken, 2003), with the effect size \eqn{f^2} of
#' Cohen (1988). \code{\link{ss_power_R2}} plans the omnibus test of all
#' predictors at once, and \code{\link{ss_power_reg_coef}} the test of one
#' coefficient with fixed predictors.
#'
#' @param true_R2_full Population squared multiple correlation coefficient of
#'   the outcome with all \code{p_full} predictors of the full model, at least
#'   0 and below 1.
#' @param true_R2_restricted Population squared multiple correlation
#'   coefficient of the outcome with the \code{p_restricted} predictors of the
#'   restricted model. The restricted model's predictors are a subset of the
#'   full model's, so \code{true_R2_restricted} cannot exceed
#'   \code{true_R2_full}.
#' @param true_R2_change Population change in \eqn{R^2} when the added
#'   predictors join the restricted model, \code{true_R2_full} minus
#'   \code{true_R2_restricted}. It is the squared semipartial correlation of
#'   the added predictor when one predictor is added, and the squared multiple
#'   semipartial correlation of the added set when several are. Give it with
#'   \code{true_R2_full} or with \code{true_R2_restricted}.
#' @param true_f2 Cohen's (1988) effect size for the increment,
#'   \eqn{f^2 = (\rho^2_{full} - \rho^2_{restricted}) / (1 - \rho^2_{full})},
#'   a nonnegative number. Give it alone, in place of the \eqn{R^2} values.
#' @param p_full Number of predictors in the full model, not counting the
#'   intercept. A whole number of at least 1.
#' @param p_restricted Number of predictors in the restricted model, not
#'   counting the intercept. A whole number from 0 to \code{p_full - 1}. Zero
#'   is the model with an intercept only, for which the test is the omnibus
#'   test of \code{\link{ss_power_R2}}.
#' @param alpha_level Type I error rate.
#' @param desired_power Desired statistical power, a single number strictly
#'   between 0 and 1 (0 and 1 themselves are refused, with a condition of
#'   class \code{dmar_desired_power_invalid}); used when a sample size is
#'   sought. With \code{specified_N} the power is computed and
#'   \code{desired_power} is not used, but a value outside (0, 1) is
#'   still refused.
#' @param specified_N Total sample size at which to report power, in place of
#'   a search for the necessary sample size. A whole number of at least
#'   \code{p_full + 2}, the smallest size that leaves an error degree of
#'   freedom.
#' @param random_predictors Whether the predictors are treated as random
#'   (\code{TRUE}, the default) or fixed (\code{FALSE}). See Details.
#'
#' @details
#' \strong{The effect.} Write \eqn{\rho^2_{full}} and
#' \eqn{\rho^2_{restricted}} (\code{true_R2_full} and
#' \code{true_R2_restricted}) for the two population squared multiple
#' correlation coefficients, \eqn{q} for the number of predictors removed
#' from the full model to give the restricted model
#' (\eqn{q = p_{full} - p_{restricted}}, the predictors being tested; row
#' \code{q} of the result, and \eqn{u} in Cohen, 1988), and \eqn{w} for the
#' number in the restricted model, \code{p_restricted}. In the partial and
#' semipartial correlation functions (\code{\link{var_partial_r}},
#' \code{\link{var_semipartial_r}}, and their planners) \eqn{p} is
#' \code{p_full}, the number of variables partialed out is \eqn{p - q}, and
#' one predictor is of interest, so \eqn{q = 1}. The change
#' \eqn{\rho^2_{full} - \rho^2_{restricted}} is the squared multiple
#' correlation of the outcome with the part of the added predictors that the
#' restricted model's predictors do not explain: the squared multiple
#' semipartial correlation of the added set (Cohen et al., 2003, Section
#' 5.4.3), and with one added predictor its squared semipartial correlation
#' (Kelley & Maxwell, 2008, Equation 20, p. 176; Aloe & Becker, 2012,
#' Equation 1).
#' So the answer to whether the change in \eqn{R^2} is the squared
#' semipartial correlation is yes, for one added predictor, and for a set it
#' is the multiple version of the same quantity. Dividing the change by
#' \eqn{1 - \rho^2_{restricted}} gives the squared multiple partial
#' correlation (the squared partial correlation for one added predictor).
#' Cohen's effect size divides the change by the full model's unexplained
#' variance,
#' \deqn{f^2 = \frac{\rho^2_{full} - \rho^2_{restricted}}{1 - \rho^2_{full}},}
#' which is also the squared multiple partial correlation divided by one minus
#' itself. Power depends on the two \eqn{R^2} values only through \eqn{f^2},
#' so the ways of giving the effect are equivalent and exactly one is given:
#' any two of the three \eqn{R^2} quantities (\code{true_R2_full},
#' \code{true_R2_restricted}, \code{true_R2_change}), or \code{true_f2}
#' alone. The same change is a larger effect when the full model explains
#' more: a change of .05 is \eqn{f^2 = .071} when the full model explains
#' .30 and \eqn{f^2 = .25} when it explains .80.
#'
#' \strong{The test.} The statistic is
#' \deqn{F = \frac{(R^2_{full} - R^2_{restricted}) / q}{(1 - R^2_{full}) / (N - p_{full} - 1)},}
#' with \eqn{q} and \eqn{N - p_{full} - 1} degrees of freedom (the form of
#' Cohen et al., 2003, Equation 5.5.1), the general model comparison
#' \emph{F} of Maxwell, Delaney, and Kelley (2027, Chapter 3) applied to the
#' two regression models. With one added predictor it is the square of the
#' \emph{t}-statistic of that predictor's coefficient. In the notation of
#' Maxwell, Delaney, and Kelley (2027), the restricted model has
#' \eqn{df_R = N - p_{restricted} - 1} error degrees of freedom and the full
#' model \eqn{df_F = N - p_{full} - 1}; the numerator degrees of freedom are
#' \eqn{df_R - df_F = q} (rows \code{q} and \code{df_1}), and
#' the denominator degrees of freedom are \eqn{df_F} (row \code{df_2}). The
#' function takes the numbers of predictors, since the error degrees of
#' freedom depend on the sample size being planned. A categorical predictor
#' with \emph{a} levels counts as \emph{a} - 1 predictors.
#'
#' \strong{Fixed predictors.} When the predictor values are held constant in
#' hypothetical replications of the study, \emph{F} follows a noncentral
#' \emph{F} distribution with \eqn{q} and \eqn{N - p_{full} - 1} degrees of
#' freedom and noncentrality \eqn{\lambda = N f^2}, where the \eqn{R^2}
#' values are those of the design itself (computed from the design's sums of
#' squares, with divisor \eqn{N}). The calculation is exact under normal
#' errors. With one added predictor this is the square of the noncentrality
#' of Kelley and Maxwell (2008, Equations 16 and 19, p. 175, and Equation 31,
#' p. 181), who give \eqn{\lambda_j = f_j \sqrt{N}} for the \emph{t}-test of
#' one coefficient with fixed predictors. That is the noncentrality of
#' \code{\link{ss_power_reg_coef}}, whose power this function reproduces when
#' one predictor is added: the noncentral \emph{t} of that function has
#' noncentrality \eqn{\sqrt{N} f}, whose square is \eqn{N f^2}. For any
#' number of added predictors \eqn{N f^2} is also the noncentrality of
#' \code{\link{ss_power_R2}} with fixed predictors.
#'
#' \strong{Random predictors.} When the predictors are sampled along with the
#' outcome from a multivariate normal population, the sample squared multiple
#' partial correlation of the outcome with the \eqn{q} added predictors, given
#' the \eqn{w} others, has the distribution that a sample squared multiple
#' correlation has with \eqn{q} predictors and a sample size of \eqn{N - w}
#' (Gatsonis & Sampson, 1989, p. 518). The critical value is the same as with
#' fixed predictors, since under the null hypothesis the statistic has the
#' central \emph{F} distribution in both models; only the power differs. The
#' power is computed from that distribution with the noncentral \emph{F}
#' approximation of Lee (1971, Section 5.1) that \code{\link{ss_power_R2}} and
#' \code{\link{ci_R2}} use. Against the exact distribution, the power by
#' Lee's approximation differed by at most .0043 over the combinations
#' examined (squared partial correlations up to .95, \code{alpha_level} from
#' .001 to .10, and one or more error degrees of freedom) and by at most .0004 with 25
#' or more error degrees of freedom; the necessary sample size differed from
#' the one the exact distribution gives, by one case, in 20 of 1,500 searches,
#' and in those the power by the exact distribution at the size returned was
#' within .0011 of \code{desired_power}. At the levels of power used in planning, power is
#' lower and the necessary sample size larger than with fixed predictors, and
#' the difference grows with the number of predictors in the restricted model
#' relative to the sample size. When power is very low the order can reverse
#' by a small amount (at most .0005 in 1,920 combinations examined, all with
#' no predictors in the restricted model, \code{alpha_level} of .001, and
#' power below .09). Predictors are random in most studies in the behavioral,
#' educational, and social sciences, so \code{random_predictors = TRUE} is
#' the default.
#'
#' \strong{Relation to the test of one regression coefficient.} With one
#' added predictor and \code{random_predictors = FALSE} this function and
#' \code{\link{ss_power_reg_coef}} plan the same test and return the same
#' sample size and power; \code{true_R2_Y_X} there is \code{true_R2_full} here
#' and \code{true_R2_Y_X_without_j} is \code{true_R2_restricted}, with
#' \code{p} equal to \code{p_full} and \code{p_restricted = p_full - 1}.
#' \code{ss_power_reg_coef()} has fixed predictors only, so with the default
#' of this function (random predictors) the sample size here is larger: 56
#' against 50 for a change from .30 to .40 with five predictors and power
#' .80 (50 is the example of Kelley & Maxwell, 2008, pp. 181--182). The two functions report the effect size differently:
#' \code{true_f2} here is \eqn{f^2}, and \code{true_f} there is
#' \eqn{f}.
#'
#' \strong{Cohen's (1988) noncentrality.} Cohen computes power from
#' \eqn{\lambda = f^2 (u + v + 1)}, with \eqn{u} the number of predictors
#' tested (\eqn{q} here) and \eqn{v} the error degrees of freedom
#' (as given by Maxwell, 2000, Equation 1). For the omnibus test
#' \eqn{q + v + 1 = N}. For an increment \eqn{v = N - q - w - 1}, so Cohen's
#' noncentrality is \eqn{f^2 (N - w)}, smaller than \eqn{N f^2} by
#' \eqn{w f^2}. Cohen's value equals the fixed-predictor calculation for
#' \eqn{q} predictors, no restricted model, and \code{specified_N} set to
#' \eqn{N - w}, and the last example reproduces it that way. That
#' \eqn{N - w} is only the number that makes the noncentrality agree; it is
#' not the study's sample size, which remains \eqn{N}. The necessary sample size by Cohen's
#' noncentrality exceeds the fixed-predictor answer of this function by at
#' most \eqn{w} cases. Maxwell (2000, p. 436) notes that \eqn{N f^2} is the
#' noncentrality for fixed predictors and reads Cohen's form as an adjustment
#' toward random predictors; the random-predictor calculation here makes that
#' adjustment from the sampling distribution itself. MBESS's
#' \code{ss.power.R2()} with a nonzero \code{Null.R2} computes this
#' quantity when \code{p} is given as the number of predictors added: the
#' sample size it reports is Cohen's total sample size less \eqn{w}.
#'
#' \strong{Error from a larger model.} The error variance of the test can be
#' taken from a model with further predictors that are not part of the
#' comparison (the Model 2 error of Cohen et al., 2003, Section 5.5.4). With
#' fixed predictors that calculation is obtained by giving \code{true_f2}
#' with the larger model's unexplained variance in its denominator,
#' \code{p_full} as the number of predictors in that larger model, and
#' \code{p_restricted} as \code{p_full} minus the number of predictors
#' tested. Set \code{random_predictors = FALSE} for this calculation. The
#' random-predictor calculation is for the comparison of two nested models
#' only.
#'
#' \strong{Edge cases.} When the change is zero the null hypothesis is true:
#' at a \code{specified_N} the power returned is \code{alpha_level}, and a
#' search stops with an error of class \code{dmar_effect_zero}, because no
#' sample size attains \code{desired_power}. A \code{specified_N} below \code{p_full + 2} leaves
#' no error degrees of freedom and is an error, as is one that is not a whole
#' number. With \code{p_restricted = 0} the restricted model explains
#' nothing, so \code{true_R2_restricted}, if given, must be 0 and
#' \code{true_R2_change}, if given, must equal \code{true_R2_full}; the
#' result is then that of \code{\link{ss_power_R2}}. A search starts at
#' \code{p_full + 2} and returns that size when it already attains
#' \code{desired_power}. The search does not go beyond a sample size of ten
#' million and stops with an error when the power is not attained by then. R
#' computes the noncentral \emph{F} distribution accurately for a
#' noncentrality parameter up to about one million and warns beyond it; no
#' power is returned from an evaluation at which R gave that warning, and
#' the error says so.
#'
#' \strong{Estimation.} To plan for a narrow confidence interval rather than
#' for power when one predictor is added, see
#' \code{\link{ss_aipe_semipartial_r}}; its \code{true_r_sp} squared is
#' \code{true_R2_change} and its \code{p} is \code{p_full} (that function
#' counts the predictor of interest and the variables partialed from it, and
#' one predictor is added, so \code{p_restricted + 1 = p_full}). DMAR has no
#' accuracy in parameter estimation planner for the change in \eqn{R^2} from
#' a set of predictors.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value} and nine
#' rows, in this order:
#' \describe{
#' \item{\code{necessary_N}}{The smallest total sample size at which power is
#'   no less than \code{desired_power}. When \code{specified_N} was given the
#'   row is named \code{specified_N} and holds that value.}
#' \item{\code{actual_power}}{Power at that sample size.}
#' \item{\code{noncentral_f_parm}}{The noncentrality parameter \eqn{N f^2}
#'   under fixed predictors; \code{NA} under random predictors, where the
#'   power does not come from a single noncentral \emph{F} with that
#'   parameter.}
#' \item{\code{true_f2}}{The planning value of Cohen's \eqn{f^2} for the
#'   increment: \code{true_f2} as given, or the value the \eqn{R^2}
#'   planning values imply.}
#' \item{\code{true_partial_R2}}{The population squared multiple partial
#'   correlation the planning values imply, \eqn{f^2 / (1 + f^2)}.}
#' \item{\code{true_R2_change}}{The population change in \eqn{R^2};
#'   \code{NA} when \code{true_f2} was given, which does not determine it.}
#' \item{\code{q}}{The number of predictors removed from the full model to
#'   give the restricted model, \code{p_full} minus \code{p_restricted}: the
#'   predictors being tested.}
#' \item{\code{df_1}}{Numerator degrees of freedom, equal to \code{q}.}
#' \item{\code{df_2}}{Denominator degrees of freedom at the sample size
#'   reported, \eqn{N - p_{full} - 1}.}
#' }
#' The result carries the \code{dmar_ss_power} class, so
#' \code{\link[generics]{tidy}} and \code{\link[generics]{glance}} apply, and
#' the attribute \code{random_predictors} records the predictor model.
#'
#' @note
#' The sample size returned is the smallest whole number at which power is no
#' less than \code{desired_power}, so \code{actual_power} is slightly above
#' the target. With random predictors that power is the one Lee's
#' approximation gives; the paragraph on random predictors in Details says
#' how closely that approximation follows the exact distribution.
#'
#' @references
#' Aloe, A. M., & Becker, B. J. (2012). An effect size for regression
#'   predictors in meta-analysis. \emph{Journal of Educational and Behavioral
#'   Statistics, 37}(2), 278--297. \doi{10.3102/1076998610396901}
#'
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum. (See Chapter 9.)
#'
#' Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). \emph{Applied
#'   multiple regression/correlation analysis for the behavioral sciences}
#'   (3rd ed.). Lawrence Erlbaum. (See Sections 5.4 and 5.5.)
#'
#' Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact power
#'   and sample size calculations. \emph{Psychological Bulletin, 106}(3),
#'   516--524.
#'
#' Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
#'   applications to multiple regression: Power and accuracy for omnibus and
#'   targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
#'   \emph{The Sage handbook of social research methods} (pp. 166--192). Sage.
#'
#' Lee, Y. S. (1971). Some results on the sampling distribution of the
#'   multiple correlation coefficient. \emph{Journal of the Royal Statistical
#'   Society, Series B, 33}(1), 117--130.
#'
#' Maxwell, S. E. (2000). Sample size and multiple regression analysis.
#'   \emph{Psychological Methods, 5}(4), 434--458.
#'   \doi{10.1037/1082-989X.5.4.434}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on the general form of the
#'   model comparison \emph{F}-test.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_R2}}, \code{\link{ss_power_reg_coef}},
#'   \code{\link{ss_aipe_semipartial_r}}, \code{\link{var_semipartial_r}},
#'   \code{\link{ci_R2}}
#'
#' @examples
#' # Five predictors are expected to explain .25 of the variance and four
#' # added predictors to raise that to .30. Necessary sample size for power
#' # .80 at alpha .01, with random predictors (the default).
#' ss_power_R2_change(true_R2_restricted = .25, true_R2_full = .30,
#'                    p_restricted = 5, p_full = 9, alpha_level = .01,
#'                    desired_power = .80)
#'
#' # The same plan with predictors fixed by design.
#' ss_power_R2_change(true_R2_restricted = .25, true_R2_full = .30,
#'                    p_restricted = 5, p_full = 9, alpha_level = .01,
#'                    desired_power = .80, random_predictors = FALSE)
#'
#' # The effect given as the change with the full model's value, and as
#' # Cohen's f squared. Power at a sample size of 90, fixed predictors.
#' ss_power_R2_change(true_R2_change = .05, true_R2_full = .30,
#'                    p_restricted = 5, p_full = 9, alpha_level = .01,
#'                    specified_N = 90, random_predictors = FALSE)
#' ss_power_R2_change(true_f2 = .05 / .70, p_restricted = 5, p_full = 9,
#'                    alpha_level = .01, specified_N = 90,
#'                    random_predictors = FALSE)
#'
#' # The effect given as the restricted model's value and the change, the
#' # way a hierarchical regression is often described. Random predictors.
#' ss_power_R2_change(true_R2_restricted = .25, true_R2_change = .05,
#'                    p_restricted = 5, p_full = 9, alpha_level = .01,
#'                    desired_power = .80)
#'
#' # One added predictor with fixed predictors is the test of its regression
#' # coefficient: the two functions return the same sample size.
#' ss_power_R2_change(true_R2_restricted = .30, true_R2_full = .40,
#'                    p_restricted = 4, p_full = 5, desired_power = .80,
#'                    random_predictors = FALSE)
#' ss_power_reg_coef(true_R2_Y_X = .40, true_R2_Y_X_without_j = .30, p = 5,
#'                   desired_power = .80)
#'
#' # Cohen's (1988) noncentrality for the design of the third example,
#' # f squared times (N - w) with w = 5: four predictors and no restricted
#' # model, with specified_N set to 90 - 5 so the noncentrality matches
#' # Cohen's. The study's sample size is still 90.
#' ss_power_R2_change(true_f2 = .05 / .70, p_restricted = 0, p_full = 4,
#'                    alpha_level = .01, specified_N = 85,
#'                    random_predictors = FALSE)
#'
#' @keywords design
#'
#' @family sample size for power
#'
#' @export

ss_power_R2_change <- function(true_R2_full = NULL, true_R2_restricted = NULL,
                               true_R2_change = NULL, true_f2 = NULL,
                               p_full, p_restricted,
                               alpha_level = .05, desired_power = .85,
                               specified_N = NULL, random_predictors = TRUE) {
  # A single finite number, the form every scalar argument below must have
  # before its range is checked.
  .is_number <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  .is_whole <- function(x) .is_number(x) && x == round(x)
  num <- function(x) format(x, digits = 15, scientific = FALSE)

  if (!.is_number(alpha_level) || alpha_level <= 0 || alpha_level >= 1) {
    .dmar_refuse_argument(
      "'alpha_level' must be a single numeric value in (0, 1).",
      arg = "alpha_level")
  }
  if (missing(p_full) || !.is_whole(p_full) || p_full < 1) {
    .dmar_refuse_argument(
      "'p_full' must be a single whole number of at least 1.", arg = "p_full")
  }
  if (missing(p_restricted) || !.is_whole(p_restricted) ||
      p_restricted < 0 || p_restricted >= p_full) {
    .dmar_refuse_argument(
      c("'p_restricted' must be a single whole number from 0 to p_full - 1. The restricted model has fewer predictors than the full model; 0 is the model with an intercept only."),
      arg = c("p_restricted", "p_full"))
  }
  if (!isTRUE(random_predictors) && !isFALSE(random_predictors)) {
    .dmar_refuse_argument(
      "'random_predictors' must be TRUE or FALSE.", arg = "random_predictors")
  }

  # The effect: two of the three R-squared quantities, or true_f2 alone.
  has_full       <- !is.null(true_R2_full)
  has_restricted <- !is.null(true_R2_restricted)
  has_change     <- !is.null(true_R2_change)
  n_R2 <- has_full + has_restricted + has_change
  if (!((is.null(true_f2) && n_R2 == 2L) || (!is.null(true_f2) && n_R2 == 0L))) {
    .dmar_refuse_argument(
      c("Give the effect in exactly one way: two of 'true_R2_full', 'true_R2_restricted', and 'true_R2_change', or 'true_f2' alone."),
      arg = c("true_R2_full", "true_R2_restricted", "true_R2_change", "true_f2"))
  }
  .check_R2 <- function(x, arg) {
    if (!.is_number(x) || x < 0 || x >= 1) {
      .dmar_refuse_argument(c(
        "'", arg, "' must be a single numeric value in [0, 1)."), arg = NULL)
    }
  }
  if (has_full)       .check_R2(true_R2_full, "true_R2_full")
  if (has_restricted) .check_R2(true_R2_restricted, "true_R2_restricted")
  if (has_change)     .check_R2(true_R2_change, "true_R2_change")

  # R2_restricted_implied is the restricted model's value the input implies;
  # it is not defined by true_f2 alone. tol_equal is the difference below
  # which two planning values are the same number entered two ways.
  tol_equal <- 1e-12
  if (has_full && has_restricted) {
    if (true_R2_restricted - true_R2_full > tol_equal) {
      .dmar_refuse_argument(
        c("'true_R2_restricted' (", num(true_R2_restricted), ") exceeds 'true_R2_full' (", num(true_R2_full), "). The restricted model's predictors are a subset of the full model's, so its squared multiple correlation cannot be the larger one."),
        arg = c("true_R2_restricted", "true_R2_full"))
    }
    # Values that differ only in the last bits (.1 + .2 against .3) are equal.
    R2_change <- max(true_R2_full - true_R2_restricted, 0)
    if (R2_change <= tol_equal) R2_change <- 0
    R2_restricted_implied <- true_R2_restricted
    f2 <- R2_change / (1 - true_R2_full)
  } else if (has_full && has_change) {
    if (true_R2_change - true_R2_full > tol_equal) {
      .dmar_refuse_argument(
        c("'true_R2_change' (", num(true_R2_change), ") exceeds 'true_R2_full' (", num(true_R2_full), "). The change is the part of 'true_R2_full' that the added predictors contribute, so it cannot be the larger one."),
        arg = c("true_R2_change", "true_R2_full"))
    }
    R2_change <- true_R2_change
    R2_restricted_implied <- true_R2_full - true_R2_change
    if (abs(R2_restricted_implied) <= tol_equal) R2_restricted_implied <- 0
    f2 <- R2_change / (1 - true_R2_full)
  } else if (has_restricted && has_change) {
    if (true_R2_restricted + true_R2_change >= 1) {
      .dmar_refuse_argument(
        c("'true_R2_restricted' plus 'true_R2_change' is 'true_R2_full' and must be below 1; the values given sum to ", num(true_R2_restricted + true_R2_change), "."),
        arg = c("true_R2_restricted", "true_R2_change"))
    }
    R2_change <- true_R2_change
    R2_restricted_implied <- true_R2_restricted
    f2 <- R2_change / (1 - (true_R2_restricted + true_R2_change))
  } else {
    if (!.is_number(true_f2) || true_f2 < 0) {
      .dmar_refuse_argument(
        "'true_f2' must be a single nonnegative numeric value.",
        arg = "true_f2")
    }
    R2_change <- NA_real_
    R2_restricted_implied <- NA_real_
    f2 <- true_f2
  }
  if (p_restricted == 0 && !is.na(R2_restricted_implied) &&
      abs(R2_restricted_implied) > 1e-12) {
    .dmar_refuse_argument(
      c("With 'p_restricted' = 0 the restricted model has an intercept only and explains nothing, so 'true_R2_restricted' must be 0 and 'true_R2_change' must equal 'true_R2_full'."),
      arg = c("p_restricted", "true_R2_restricted", "true_R2_change"))
  }
  # desired_power is used only by the sample size search, but an invalid
  # value is refused on either branch, before any work.
  if (is.null(specified_N) || !is.null(desired_power)) {
    .dmar_check_desired_power(desired_power)
  }
  # A total N below p_full + 2 leaves no error degrees of freedom; a
  # fractional N is not a sample size.
  if (!is.null(specified_N)) {
    specified_N <- .check_whole_n(specified_N, "specified_N", p_full + 2L)
  }

  # Power at a given N, NA where R's noncentral F routine warned of lost
  # accuracy (R/R2_internals.R); a missing power is refused below.
  .power_at_N <- function(N_now) {
    .R2_change_power_guarded(N_now, f2, p_full = p_full,
                             p_restricted = p_restricted,
                             alpha_level = alpha_level,
                             random_predictors = random_predictors)
  }
  .imprecise_message <- function(N_now) .R2_imprecise_message(N_now, f2)

  if (is.null(specified_N)) {
    if (f2 == 0) {
      .dmar_refuse(c(
        "The change in R-squared is 0, so the null hypothesis is true: ",
        "power equals 'alpha_level' at every sample size and no sample ",
        "size attains 'desired_power'. Give the nonzero change the study ",
        "is planned to detect, or give 'specified_N' to have the power at ",
        "that sample size (which is 'alpha_level') reported."),
        "effect_zero")
    }
    found <- .R2_power_search_N(.power_at_N, N_min = p_full + 2,
                                desired_power = desired_power,
                                imprecise_message = .imprecise_message)
    N_out <- found$N
    power_out <- found$power
    size_term <- "necessary_N"
  } else {
    N_out <- specified_N
    power_out <- .power_at_N(N_out)
    if (is.na(power_out)) {
      .dmar_refuse(.imprecise_message(N_out), "computation_failed")
    }
    size_term <- "specified_N"
  }

  out <- data.frame(
    term  = c(size_term, "actual_power", "noncentral_f_parm", "true_f2",
              "true_partial_R2", "true_R2_change", "q", "df_1", "df_2"),
    value = c(N_out, power_out,
              if (random_predictors) NA_real_ else N_out * f2,
              f2, f2 / (1 + f2), R2_change,
              p_full - p_restricted,
              p_full - p_restricted, N_out - p_full - 1),
    stringsAsFactors = FALSE
  )
  attr(out, "random_predictors") <- random_predictors
  .as_dmar_tbl(out, subclass = "dmar_ss_power",
               measured_terms = c("true_f2", "true_partial_R2",
                                  "true_R2_change"))
}
