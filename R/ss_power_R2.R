#' Plan Sample Size to Make the Test of the Squared Multiple Correlation Coefficient Sufficiently Powerful
#'
#' @description
#' Determine the necessary sample size for the omnibus test of the squared
#' multiple correlation coefficient (\eqn{R^2}), or the realized statistical
#' power given a specified sample size, under either fixed or random
#' predictors. The fixed-predictors path uses Cohen's (1988) noncentral
#' \emph{F} formulation; the random-predictors path uses the three-moment
#' noncentral \emph{F} approximation of Lee (1971, Section 5.1) to the
#' sampling distribution of the sample \eqn{R^2} under joint multivariate
#' normality.
#'
#' @param true_R2 Population squared multiple correlation coefficient the
#'   sample size or power is planned for. The planner treats this value as
#'   the truth; \code{\link{ss_power_R2_sensitivity}} shows what the plan
#'   delivers when the truth (\code{true_R2} there) differs from the value
#'   planned with (\code{supposed_R2} there).
#' @param alpha_level Type I error rate
#' @param desired_power Desired degree of statistical power, a single
#'   number strictly between 0 and 1. The values 0 and 1 are refused, with
#'   a condition of class \code{dmar_desired_power_invalid}.
#' @param p The number of predictor variables
#' @param specified_N The sample size used to calculate power (rather than determine necessary sample size). This is the \emph{total} sample size across all groups or observations.
#' @param true_f2 Cohen's (1988) effect size for multiple regression: \code{true_R2}/(1 - \code{true_R2}).
#' @param null_R2 The null value of the omnibus test, which must be 0. Any other value stops with an error of class \code{dmar_null_R2_unsupported} that explains why and gives the call to \code{\link{ss_power_R2_change}}; see Details.
#' @param random_predictors Whether the predictor variables are treated as random (\code{TRUE}, the default) or fixed (\code{FALSE}). See Details.
#' @param print_iter If the progress of the iterative procedure is printed to the screen as the iterations are occurring
#'
#' @details
#' Determine the necessary sample size given a particular
#' \code{true_R2}, \code{alpha_level}, \code{p}, and
#' \code{desired_power}. Alternatively, given \code{true_R2},
#' \code{alpha_level}, \code{p}, and \code{specified_N}, the function can
#' be used to determine the statistical power.
#'
#' \strong{Fixed vs.\ random predictors.} The two regression models give
#' \emph{different} sampling distributions for the omnibus \eqn{F}-statistic,
#' and so different power. Under fixed predictors the design matrix is
#' treated as constant in hypothetical replications of the study, and
#' \eqn{F} follows a noncentral \emph{F} with \eqn{p} and \eqn{N-p-1}
#' degrees of freedom and noncentrality \eqn{\lambda = N \cdot f^2},
#' where \eqn{f^2 = \rho^2 / (1 - \rho^2)} and \eqn{\rho^2} is \code{true_R2}
#' (Cohen, 1988). Cohen's
#' noncentrality is \eqn{\lambda = f^2 (u + v + 1)}, with \eqn{u} and
#' \eqn{v} the numerator and denominator degrees of freedom (as given by
#' Maxwell, 2000, Equation 1); for the omnibus test \eqn{u = p} and
#' \eqn{v = N - p - 1}, so \eqn{u + v + 1 = N} and the two forms are the
#' same. Under random
#' predictors the design matrix is itself a draw from a joint multivariate
#' normal distribution, and the unconditional distribution of the sample
#' \eqn{R^2} is due to Fisher (1928); \code{ss_power_R2()} uses the
#' three-moment noncentral \emph{F} approximation of Lee (1971, Section 5.1)
#' to that distribution, the same approximation \code{ci_R2()} uses for
#' random-predictor confidence intervals. Gatsonis and Sampson (1989) give
#' the exact calculation for random predictors and compare it with Cohen's
#' tables (the 1977 edition of the book cited below), whose sample sizes were
#' slightly smaller than the exact ones for most designs with up to 10
#' predictors. The sample size this function returns under random predictors
#' is the smallest at which the power by Lee's approximation attains
#' \code{desired_power}. Against the exact distribution, the power by
#' Lee's approximation differed by at most .0043 over the combinations
#' examined (squared correlations up to .95, \code{alpha_level} from .001 to
#' .10, and one or more error degrees of freedom) and by at most .0004 with 25
#' or more error degrees of freedom; the necessary sample size differed from
#' the one the exact distribution gives, by one case, in 20 of 1,500 searches,
#' and in those the power by the exact distribution at the size returned was
#' within .0011 of \code{desired_power}. In the behavioral,
#' educational, and social sciences predictor variables are almost always
#' random, so the default is \code{random_predictors = TRUE}; pass
#' \code{random_predictors = FALSE} for designs in which the predictor
#' variables are fixed by design (for example, planned dosing levels).
#'
#' \strong{Tests other than the omnibus test.} \code{ss_power_R2()} plans
#' the test that the population squared multiple correlation coefficient of
#' all \code{p} predictors is zero. The test that predictors added to a
#' model contribute nothing, a test of the change in \eqn{R^2}, is planned
#' with \code{\link{ss_power_R2_change}}, which takes the number of
#' predictors in both models.
#'
#' \code{null_R2} exists so that a call written for MBESS with a nonzero
#' null value is answered with an explanation; it must be 0. A nonzero
#' value stops with an error of class \code{dmar_null_R2_unsupported} (see
#' \code{\link{dmar_metadata}}), for these reasons. In MBESS a nonzero
#' null value \eqn{\rho_0^2} forms the effect size
#' \eqn{(\rho^2 - \rho_0^2) / (1 - \rho^2)}, which is Cohen's (1988)
#' \eqn{f^2} for the increment in \eqn{R^2} when predictors are added to a
#' model that already explains \eqn{\rho_0^2}, and refers it to an
#' \emph{F} distribution with \eqn{p} and \eqn{N - p - 1} degrees of
#' freedom and noncentrality \eqn{N f^2}. What that calculation is
#' depends on how \code{p} is read. With \code{p} read as the number of
#' predictors added, it is Cohen's calculation,
#' \eqn{\lambda = f^2 (u + v + 1)}, for a study whose total sample size is
#' \eqn{N} plus the number of predictors in the restricted model: the
#' sample size it reports is Cohen's total sample size less that number.
#' With \code{p} read as all of the predictors, as this function defines
#' it, the calculation is not the test of the increment, whose numerator
#' degrees of freedom are the number of predictors added and not \code{p}.
#' The argument alone cannot say which calculation
#' is meant, because the test of an increment needs the number of
#' predictors in the restricted model, so the value is refused and the
#' message writes the call to \code{ss_power_R2_change()}.
#'
#' The number that MBESS returns for a population value \eqn{a} and a
#' null value \eqn{b} is \code{ss_power_R2(true_f2 = (a - b) / (1 - a),
#' p = p, random_predictors = FALSE)}, the omnibus test of \code{p} predictors with
#' an effect of that size; the last example shows it. Whether the population
#' squared multiple correlation coefficient of all \code{p} predictors
#' exceeds a stated value above zero is a third question. DMAR does not
#' plan sample size for that test; the test itself rejects
#' when the one-sided lower confidence limit from \code{\link{ci_R2}}
#' (with \code{alpha_lower = alpha_level} and \code{alpha_upper = 0})
#' exceeds the stated value.
#'
#' \strong{Limit of the calculation.} R computes the noncentral \emph{F}
#' distribution accurately for a noncentrality parameter up to about one
#' million, and warns beyond it. \code{ss_power_R2()} returns no power from
#' an evaluation at which R gave that warning: a power requested at such a
#' \code{specified_N} is refused, and a search whose answer lies there
#' stops and says so. A search that only passes through such sizes on its
#' way to a smaller answer returns that answer, which is evaluated without
#' the warning. The search does not go beyond a sample size of ten million,
#' and stops with an error when \code{desired_power} is not attained by
#' then. With \code{true_R2} of 0 the null hypothesis is true: the
#' power at a \code{specified_N} is \code{alpha_level}, and a search stops
#' with an error of class \code{dmar_effect_zero} because no sample size
#' attains \code{desired_power}.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}.
#' For an \eqn{N} search the rows are \code{necessary_N},
#' \code{actual_power}, \code{noncentral_f_parm} (only meaningful for
#' \code{random_predictors = FALSE}; \code{NA} otherwise), and
#' \code{true_f2} (Cohen's \eqn{f^2 = \rho^2 / (1 - \rho^2)} at
#' \code{true_R2}). For a
#' power-at-specified-\code{N} computation the first row is
#' \code{specified_N} instead of \code{necessary_N}.
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Fisher, R. A. (1928). The general sampling distribution of the multiple
#'   correlation coefficient. \emph{Proceedings of the Royal Society of
#'   London, Series A, 121}(788), 654--673. \doi{10.1098/rspa.1928.0224}
#'
#' Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact
#'   power and sample size calculations. \emph{Psychological Bulletin,
#'   106}(3), 516--524.
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple
#'   correlation coefficient: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{Multivariate Behavioral Research, 43}(4),
#'   524--555. \doi{10.1080/00273170802490632}
#'
#' Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
#'   applications to multiple regression: Power and accuracy for omnibus
#'   and targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen
#'   (Eds.), \emph{The Sage handbook of social research methods}
#'   (pp. 166--192). Sage.
#'
#' Lee, Y. S. (1971). Some results on the sampling distribution of the
#'   multiple correlation coefficient. \emph{Journal of the Royal
#'   Statistical Society, Series B, 33}(1), 117--130.
#'
#' Maxwell, S. E. (2000). Sample size and multiple regression analysis.
#'   \emph{Psychological Methods, 5}(4), 434--458.
#'   \doi{10.1037/1082-989X.5.4.434}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{R^2} as a model comparison
#'   effect size.)
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#'   for statistical power and accuracy in parameter estimation.
#'   \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
#'   planning for more accurate statistical power: A method adjusting
#'   sample effect sizes for publication bias and uncertainty.
#'   \emph{Psychological Science, 28}(11), 1547--1562.
#'   \doi{10.1177/0956797617723724}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' When determining sample size for a desired degree of power, there will
#' always be a slightly larger degree of actual power. This is the case
#' because the algorithm employed determines sample size until the actual
#' power is no less than the desired power (given sample size is a whole
#' number power will almost certainly not be exactly the specified value).
#' This is the same as other statistical power procedures that return
#' whole numbers for necessary sample size.
#'
#' @seealso
#' \code{\link{ss_power_R2_change}} for the test of predictors added to a
#' model, \code{\link{ss_aipe_R2}}, \code{\link{ss_power_R2_sensitivity}},
#' \code{\link{ss_power_reg_coef}}, \code{\link{ci_nc_F}},
#' \code{\link{ci_R2}}
#'
#' @examples
#' # Random predictors (default; appropriate for most behavioral / social
#' # science applications).
#' ss_power_R2(true_R2 = .5, alpha_level = .05, desired_power = .85, p = 5)
#'
#' # Fixed predictors (Cohen 1988): predictor variables fixed by design.
#' ss_power_R2(true_R2 = .5, alpha_level = .05, desired_power = .85,
#'             p = 5, random_predictors = FALSE)
#'
#' # Effect size input (Cohen's f^2).
#' ss_power_R2(true_f2 = 1, alpha_level = .05, desired_power = .85, p = 5)
#'
#' # Realized power at a specified N.
#' ss_power_R2(true_R2 = .5, specified_N = 15, alpha_level = .05,
#'             desired_power = .85, p = 5)
#'
#' # The number MBESS reported for a population value of .30 against a null
#' # value of .25 with p = 4: the omnibus test of four fixed predictors with
#' # Cohen's f squared for that increment. The test of the increment itself
#' # is planned with ss_power_R2_change().
#' ss_power_R2(true_f2 = (.30 - .25) / (1 - .30), p = 4, specified_N = 85,
#'             alpha_level = .01, random_predictors = FALSE)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export

ss_power_R2 <- function(true_R2 = NULL, alpha_level = .05,
                        desired_power = .85, p, specified_N = NULL,
                        true_f2 = NULL, null_R2 = 0,
                        random_predictors = TRUE,
                        print_iter = FALSE) {
  # A single finite number, the form every scalar argument below must have
  # before its range is checked.
  .is_number <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)

  if (!.is_number(alpha_level) || alpha_level <= 0 || alpha_level >= 1) {
    .dmar_refuse_argument(
      "'alpha_level' must be a single numeric value in (0, 1).",
      arg = "alpha_level")
  }
  if (!.is_number(p) || p != round(p) || p < 1) {
    .dmar_refuse_argument(
      "'p' must be a single whole number of at least 1.", arg = "p")
  }
  if (!isTRUE(random_predictors) && !isFALSE(random_predictors)) {
    .dmar_refuse_argument(
      "'random_predictors' must be TRUE or FALSE.", arg = "random_predictors")
  }
  if (!.is_number(null_R2) || null_R2 < 0 || null_R2 >= 1) {
    .dmar_refuse_argument(
      "'null_R2' must be a single numeric value in [0, 1).", arg = "null_R2")
  }
  if (is.null(true_R2) && is.null(true_f2)) {
    .dmar_refuse_argument(
      "Specify either 'true_R2' or 'true_f2'.", arg = c("true_R2", "true_f2"))
  }
  if (!is.null(true_R2) && !is.null(true_f2)) {
    .dmar_refuse_argument(
      "Since you specified 'true_R2', do not specify 'true_f2'.",
      arg = c("true_R2", "true_f2"))
  }
  if (!is.null(true_R2) &&
      (!.is_number(true_R2) || true_R2 < 0 || true_R2 >= 1)) {
    .dmar_refuse_argument(
      "'true_R2' must be a single numeric value in [0, 1).", arg = "true_R2")
  }
  if (!is.null(true_f2) && (!.is_number(true_f2) || true_f2 < 0)) {
    .dmar_refuse_argument(
      "'true_f2' must be a single nonnegative numeric value.", arg = "true_f2")
  }
  # The function plans the omnibus test, whose null value is 0. A nonzero
  # null_R2 had three possible readings (the MBESS arithmetic, the test of
  # an increment, a minimum effect test) that give different numbers, and
  # the argument alone cannot say which is meant. It is refused with the
  # reason and the calls that answer each reading.
  if (null_R2 != 0) {
    .dmar_refuse(.ss_power_R2_null_message(null_R2, true_R2, true_f2, p),
                 "null_R2_unsupported", null_R2 = null_R2,
                 true_R2 = true_R2, true_f2 = true_f2, p = p)
  }
  if (!is.null(true_f2)) {
    true_R2 <- true_f2 / (1 + true_f2)
  }
  # desired_power is used only by the sample size search, but an invalid
  # value is refused on either branch, before any work.
  if (is.null(specified_N) || !is.null(desired_power)) {
    .dmar_check_desired_power(desired_power)
  }
  # A total N below p + 2 leaves no residual degrees of freedom; a fractional
  # N is not a sample size. Reject both rather than returning 0, NaN, or a
  # power table for a size the design cannot have.
  if (!is.null(specified_N)) {
    specified_N <- .check_whole_n(specified_N, "specified_N", p + 2L)
  }

  # Cohen's f^2, rho^2 / (1 - rho^2).
  f_Squared <- true_R2 / (1 - true_R2)

  # Power at a given N: the omnibus test is the model comparison test of all
  # p predictors against a model with none, so the calculation is the shared
  # one of R/R2_internals.R with no predictors in the restricted model. It
  # returns NA where R's noncentral F routine warned of lost accuracy, and
  # the callers below refuse a missing power.
  .power_at_N <- function(N_now) {
    .R2_change_power_guarded(N_now, f_Squared, p_full = p, p_restricted = 0,
                             alpha_level = alpha_level,
                             random_predictors = random_predictors)
  }
  .imprecise_message <- function(N_now) .R2_imprecise_message(N_now, f_Squared)

  if (is.null(specified_N)) {
    # With no effect the power is alpha_level at every sample size. Stop
    # before the search.
    if (true_R2 == 0) {
      .dmar_refuse(c(
        "With 'true_R2' of 0 the null hypothesis is true, so power equals ",
        "'alpha_level' at every sample size and no sample size attains ",
        "'desired_power'. Give the nonzero squared multiple correlation ",
        "the study is planned to detect, or give 'specified_N' to have ",
        "the power at that sample size (which is 'alpha_level') reported."),
        "effect_zero")
    }

    .power_at_N_reported <- function(N_now) {
      power_now <- .power_at_N(N_now)
      if (print_iter && !is.na(power_now)) {
        cat(c(Current.Power = power_now,
              Current.NC.F.Parm = ifelse(random_predictors, NA, N_now * f_Squared),
              Current.N = N_now), "\n")
      }
      power_now
    }

    # The smallest admissible size leaves one error degree of freedom.
    found <- .R2_power_search_N(.power_at_N_reported, N_min = p + 2,
                                desired_power = desired_power,
                                imprecise_message = .imprecise_message)
    N_i <- found$N
    Actual_Power <- found$power
    VALUE <- data.frame(
      term  = c("necessary_N", "actual_power", "noncentral_f_parm", "true_f2"),
      value = c(N_i, Actual_Power,
                ifelse(random_predictors, NA_real_, N_i * f_Squared),
                f_Squared)
    )
  } else {
    actual_power_specified_N <- .power_at_N(specified_N)
    if (is.na(actual_power_specified_N)) {
      .dmar_refuse(.imprecise_message(specified_N), "computation_failed")
    }
    VALUE <- data.frame(
      term  = c("specified_N", "actual_power", "noncentral_f_parm", "true_f2"),
      value = c(specified_N, actual_power_specified_N,
                ifelse(random_predictors, NA_real_, specified_N * f_Squared),
                f_Squared)
    )
  }
  class(VALUE) <- c("dmar_ss_power", "dmar_tbl", "data.frame")
  return(VALUE)
}

# The message of the dmar_null_R2_unsupported refusal of ss_power_R2(). The
# values of the call are written into the calls the message suggests, so that
# each can be completed and run. With 'true_f2' given in place of
# 'true_R2' the suggested calls carry it in place of the R-squared
# arguments: DMAR 1.0.x with fixed predictors took a 'true_f2' given with a
# nonzero null as the effect size of the increment.
.ss_power_R2_null_message <- function(null_R2, true_R2, true_f2, p) {
  num <- function(x) format(x, digits = 15, scientific = FALSE)
  b <- num(null_R2)
  p_text <- num(p)
  no_increment <- FALSE
  if (is.null(true_f2)) {
    a <- num(true_R2)
    effect_args <- paste0("true_R2_full = ", a, ", true_R2_restricted = ", b)
    f2_text <- paste0("(", a, " - ", b, ") / (1 - ", a, ")")
    no_increment <- true_R2 <= null_R2
  } else {
    f2_num <- num(true_f2)
    effect_args <- paste0("true_f2 = ", f2_num)
    f2_text <- f2_num
  }
  # The calls to make instead. With 'true_R2' at or below 'null_R2' there is
  # no increment, and the two calls would stop on a negative effect size.
  calls_text <- if (no_increment) {
    paste0(
      "'true_R2' (", a, ") is not above 'null_R2' (", b, "), so there is no ",
      "increment to test. The test of an increment in R-squared, for a full ",
      "model that explains more than the restricted one, has its own ",
      "function, ss_power_R2_change(), which takes the number of predictors ",
      "in both models."
    )
  } else {
    paste0(
      "The test of an increment needs the number of predictors in ",
      "both models and has its own function: ss_power_R2_change(", effect_args,
      ", p_full = (predictors in the full model), p_restricted = (predictors ",
      "in the restricted model)), with 'alpha_level', 'desired_power' or ",
      "'specified_N', and 'random_predictors' as given here. If ", p_text,
      " counted all of the predictors, 'p_full' is ", p_text, "; if it counted ",
      "the predictors added, 'p_full' is ", p_text, " plus 'p_restricted'. ",
      "The number MBESS returned is ss_power_R2(true_f2 = ", f2_text,
      ", p = ", p_text, ", random_predictors = FALSE)."
    )
  }
  paste0(
    "'null_R2' must be 0. ss_power_R2() plans the omnibus test that the ",
    "population squared multiple correlation coefficient of all 'p' ",
    "predictors is zero, and a 'null_R2' of ", b, " was given. In MBESS a ",
    "nonzero 'null_R2' formed the effect size (true_R2 - null_R2) / ",
    "(1 - true_R2), which is Cohen's (1988) f-squared for the increment ",
    "in R-squared when predictors are added to a model that already ",
    "explains 'null_R2', and referred it to an F distribution with 'p' and ",
    "N - p - 1 degrees of freedom. With 'p' read as the number of predictors added, that is ",
    "Cohen's calculation for a total sample size of N plus the number of ",
    "predictors in the restricted model. With 'p' read as all of the ",
    "predictors, as this function defines it, it is not the test of the ",
    "increment, whose numerator degrees of freedom are the number of ",
    "predictors added. ", calls_text,
    " If the question is ",
    "instead whether the population squared multiple correlation ",
    "coefficient exceeds ", b, ", DMAR does not plan sample size for that ",
    "test; the test rejects when the one-sided lower ",
    "confidence limit from ci_R2() exceeds ", b, "."
  )
}
