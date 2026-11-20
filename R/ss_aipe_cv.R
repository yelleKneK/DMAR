#' Sample Size Planning for the Coefficient of Variation Given the Goal of Accuracy in Parameter Estimation Approach to Sample Size Planning
#'
#' @description
#' Takes a planning value of the population coefficient of variation and a
#' desired confidence interval width, and returns the smallest sample size
#' at which the interval of \code{\link{ci_cv}} is expected to be that
#' narrow (Kelley, 2007). With \code{assurance} it returns the sample size
#' of the published procedure of that article, at which the realized
#' two-sided interval is that narrow with about the stated probability.
#' The interval and the plan assume a normal population (see Details).
#'
#' The population coefficient of variation may be given directly as
#' \code{true_cv} or through \code{true_mean} and \code{true_sd}, in which case
#' \code{true_cv} is taken as \code{true_sd / true_mean}. It must be positive.
#'
#' @param true_cv Population coefficient of variation the sample size is
#'   planned for (kappa in Kelley, 2007). The planner treats this value as
#'   the truth; \code{\link{ss_aipe_cv_sensitivity}} shows what the plan
#'   delivers when the truth (\code{true_cv} there) differs from the value
#'   planned with (\code{supposed_cv} there).
#' @param width Desired (full) width of the confidence interval
#' @param conf_level Confidence interval coverage; 1-Type I error rate. A
#'   level above .999999998, whose tail areas are below 1e-9, is refused
#'   (see "Small tail areas" in \code{\link{ci_nc_t}})
#' @param assurance Optional. The probability that the realized interval
#'   is no wider than \code{width}: a single number strictly between 0 and
#'   1 (for example .80, .90, or .99). A value such as .50 may be used; it
#'   asks that the realized interval be no wider than \code{width} in about
#'   half of the studies that could be run, where the plan without
#'   \code{assurance} asks that the width at the expected value of the
#'   sample coefficient be no larger than \code{width}. Near .50 the
#'   realized probability can be several points below the request and the
#'   sample size can be one smaller than without \code{assurance}; see
#'   Details for the procedure and the numbers. Any other value is
#'   refused with a condition of class \code{dmar_assurance_invalid}.
#'   \code{assurance} cannot be given with \code{alpha_lower = 0}, a
#'   one-sided upper confidence bound (see \strong{Edge cases}).
#' @param true_mean Population mean, a positive number (specified with \code{true_sd} when \code{true_cv} is not specified)
#' @param true_sd Population standard deviation, a positive number (specified with \code{true_mean} when \code{true_cv} is not specified)
#' @param alpha_lower Type I error for the lower confidence limit, given
#'   with \code{alpha_upper} when \code{conf_level} is \code{NULL}: a
#'   single number that is at least 0 and below 1
#' @param alpha_upper Type I error for the upper confidence limit. A
#'   nonzero tail area below 1e-9 is refused, as in \code{\link{ci_cv}}
#'
#' @details
#' \strong{Estimand and interval.} The quantity planned for is the
#' population coefficient of variation \eqn{\kappa = \sigma / \mu} of a
#' normal population with a positive mean, estimated from \eqn{N}
#' independent observations by \eqn{k = s / M}. The interval is the one
#' \code{\link{ci_cv}} reports, from the noncentral \emph{t} distribution
#' of \eqn{\sqrt{N} / k}, two-sided at \code{conf_level} (or with the tail
#' areas \code{alpha_lower} and \code{alpha_upper} when \code{conf_level}
#' is \code{NULL}).
#'
#' \strong{The expected-width plan.} At each candidate sample size \eqn{N}
#' the width is evaluated at the value Kelley (2007, Equation 3) gives
#' for the expected value of the sample coefficient of variation given
#' \code{true_cv} and \eqn{N}, \deqn{E[k \mid \kappa, N] \approx \kappa
#' \left[1 - \frac{1}{4(N - 1)} + \frac{\kappa^2}{N} +
#' \frac{1}{2(N - 1)^2}\right],}
#' and the smallest \eqn{N} whose width at that value is no larger than
#' \code{width} is returned. The expression is the approximation of Sokal
#' and Braumann (1980, p. 51, Expression 3) for normal data. It stops at
#' the third-order term and is above the expected value at small
#' \eqn{N}: for a coefficient of .25, by 3.1 percent at \eqn{N = 5}, 0.6
#' percent at \eqn{N = 10}, and 0.1 percent at \eqn{N = 20}. The function
#' keeps it because the published tables were computed with it. The
#' search starts at \eqn{N = 4}. The interval is computed exactly; for
#' how this relates to the sample sizes tabled in Kelley (2007), see
#' \strong{The published tables} below.
#'
#' The width of the interval is not a linear function of the sample
#' coefficient of variation, so the width at the value of Equation 3 is
#' not exactly the expected width. Over 62 conditions (coefficients of
#' variation of .05 to .50, widths of .025 to .20, and confidence levels
#' of .90, .95, and .99, with plans from 5 to 1,500) the returned \eqn{N}
#' was the smallest \eqn{N} with an expected width no larger than
#' \code{width} in 47 and one smaller in 15; where it was one smaller the
#' expected width exceeded \code{width} by at most 0.7 percent.
#'
#' \strong{The assurance step.} With \code{assurance} the function runs
#' the computation that produced the tables of Kelley (2007). At the
#' expected-width sample size \eqn{N_0} it finds \eqn{\kappa_\gamma}, the
#' value that the sample coefficient exceeds with probability
#' \code{1 - assurance}, from the noncentral \emph{t} distribution of
#' Equation 5 of the article. It then takes the value of Equation 3 at
#' \eqn{\kappa_\gamma} and \eqn{N_0} as the planning value and runs the
#' expected-width search for it, which applies Equation 3 to that value
#' again at each \eqn{N}. The quantile and the intervals are computed
#' exactly (see \code{\link{ci_nc_t}}). The article's text (p. 759)
#' describes the search as run for \eqn{\kappa_\gamma} itself.
#' Recomputed with R's own noncentral \emph{t} functions, as in 2007, the
#' computation described here gives the printed value in 524 of the 528
#' assurance cells of Tables 1 to 3 (the other four are in the two rows
#' that show a misprint), where the search for \eqn{\kappa_\gamma} itself
#' gives the printed value in 365. The two computations give different
#' plans in 160 cells.
#'
#' \strong{The published tables.} Kelley (2007) computed Tables 1 to 3
#' (pp. 760--762) with R's noncentral \emph{t} distribution, which is a
#' normal approximation above a noncentrality of 37.62. Here the
#' noncentrality is \eqn{\sqrt{N}/k}, which narrow intervals pass, and the
#' function computes there exactly, so it does not return every printed
#' value. Two rows of the tables show a misprint. Table 1 (90 percent
#' confidence), assurance .80, width .125, coefficient of variation .40
#' prints 127, between 69 and 122 in its row; the procedure gives 93.
#' Table 3 (99 percent confidence), assurance .80, width .010 prints
#' 14,279; 28,321; 50,209; 38,120; and 82,817 for coefficients of variation
#' of .30 to .50, which do not increase; the procedure gives 14,277;
#' 20,473; 28,320; 38,118; and 50,207, so the printed 28,321 and 50,209 are
#' its values for .40 and .50 (to within the approximation), one and two
#' columns to the left of where they belong. Apart from the four
#' misprinted cells (127; 28,321; 50,209; and 82,817), the function
#' returns the printed value in 468 of the 788 cells, one to three fewer
#' cases in 232, and one to three more in 88. It returns every printed
#' value with
#' a width of .125 or more at 90 and 95 percent confidence, the misprint
#' aside, and 9 of the 120 with a width of .075 or less at 99 percent.
#' The worked example of the article (pp. 759--760), a coefficient of
#' variation of .25, a width of .10, and 99 percent confidence, gives 101,
#' and 142 with an assurance of .99, where the article gives 99 and 141:
#' at 99 the interval at the expected value is .1006 wide, and .0999 by the
#' approximation.
#'
#' The procedure is a close approximation and not an exact search, so the
#' probability that the realized interval is no wider than \code{width} at
#' the returned sample size is near \code{assurance} and not exactly it. In
#' the 488 printed cells with an assurance and at most 5,000 cases that
#' probability was .744 to .902 where .80 was requested (below .80 in 37 of
#' 244 cells) and .9916 to .9997 where .99 was requested. Away from the
#' published tables, and more so as the request falls toward .50, it can
#' fall below the request: for a coefficient of variation of .22, a width
#' of .09, a 95 percent interval, and an assurance of .60 the function
#' returns 57, at which the probability is .587, and one more case gives
#' .621. In 105 plans computed at 95 percent confidence (coefficients of
#' variation from .05 to .50, widths of a quarter, a half, and three
#' quarters of the coefficient) the probability was below the request by
#' at most .038 for a request of .50, .036 for .55, .058 for .60, .032 for
#' .70, and .001 for .80, and it was not below the request in any plan for
#' .90 and .99. A request of .50 does not reproduce the plan without
#' \code{assurance}: in 12 of the 15 plans it returned a sample size one
#' smaller than the plan without \code{assurance} does, and for a
#' coefficient of variation of .05 and a width of .0375 it returns 17, at
#' which the probability is .462, where the plan without \code{assurance}
#' returns 18. In the other direction, when \code{width} is close to
#' \code{true_cv} and \code{assurance} is high the returned sample size
#' can be well above the smallest sufficient one: for a coefficient of
#' variation of .5, a width of .5, a 95 percent interval, and an assurance
#' of .99 the function returns 57, and the probability reaches .99 at 37.
#' \code{\link{ss_aipe_cv_sensitivity}} estimates the probability for any
#' plan by simulation.
#'
#' \strong{Populations that are not normal.} The interval and the plan
#' assume a normal population. With a skewed population the coverage of
#' the interval departs from the stated level, and a larger sample does
#' not repair it (see \code{\link{ci_cv}}). In 40,000 samples at each of
#' the three plans for a coefficient of variation of .50, a width of .20,
#' and 95 percent confidence (79 without \code{assurance}, 97 with an
#' assurance of .80, and 142 with .99), coverage was .97 for a gamma
#' population with that coefficient of variation, and for a lognormal one
#' it was .915, .913, and .905, lowest at the largest plan. The proportion
#' of intervals no wider than \code{width} was .55, .88, and 1.00 for the
#' gamma population and .59, .84, and .99 for the lognormal one. For a
#' normal population that probability, computed from the noncentral
#' \emph{t} distribution, is .53, .84, and 1.00.
#'
#' \strong{Edge cases.} When \code{width} is already met at \eqn{N = 4}
#' the function returns 4. A \code{width} that no sample size reaches
#' within 100,000 steps of the search stops with an error. With
#' \code{alpha_upper = 0} the interval is a one-sided lower confidence bound
#' whose upper limit is infinite at every sample size, so no plan exists,
#' and the call is refused at once with a condition of class
#' \code{dmar_argument_invalid}. With \code{alpha_lower = 0}, an upper
#' bound, the lower limit is 0 and the width is the bound itself, which
#' comes down toward the coefficient of variation as \eqn{N} grows, so a
#' \code{width} larger than \code{true_cv} can be planned without
#' \code{assurance}. With \code{assurance} the call is refused, with
#' class \code{dmar_argument_invalid}: the assurance step was written
#' for a two-sided interval and does not give the requested probability
#' for a bound. For a coefficient of variation of .2, a width of .25,
#' \code{conf_level = NULL}, \code{alpha_lower = 0}, and
#' \code{alpha_upper = .05}, where the plan without \code{assurance} is
#' 37, the step gives 33, 91, and 219 for an assurance of .50, .80, and
#' .90. The probability that the bound is no larger than the width is
#' .470, .878, and .998 at those sizes, and the smallest sample sizes at
#' which it reaches the request are 36, 74, and 98. With
#' \code{alpha_upper} below .5 the bound at the value of Equation 3
#' stays above \code{true_cv} at every sample size the search can reach:
#' over \code{true_cv} from .001 to 5 and \eqn{N} from 4 to 100,004 the
#' smallest bound was 1.000006 times \code{true_cv}. (The value of
#' Equation 3 is itself below \code{true_cv} when \code{true_cv} is
#' below \eqn{\sqrt{N(N - 3)} / (2(N - 1))}, which is .33 at
#' \eqn{N = 4}, .46 at \eqn{N = 10}, and .49 at \eqn{N = 20}, and
#' approaches .5.) A \code{width} at or below \code{true_cv} is therefore
#' never reached, and the call is refused at once, also with class
#' \code{dmar_argument_invalid}. With a large
#' coefficient of variation the first sample sizes of the search can be so
#' small that the upper limit of the interval is infinite (at 99 percent
#' confidence and a coefficient of .45, \eqn{N = 4} is one). An infinite
#' interval is never narrow enough, so the search goes on, and the
#' interval at the returned sample size is finite.
#'
#' \strong{Comparison with MBESS.} MBESS 3.4.0 and later evaluate the
#' width at \code{true_cv} multiplied by \eqn{(1 + 1/(4N))}, the nearly
#' unbiased estimate of Kelley (2007, Equation 4) applied to the
#' population value, which is larger than \code{true_cv} where the
#' expected value is smaller, and use R's approximation to the noncentral
#' \emph{t} beyond a noncentrality of 37.62. Their sample sizes are one
#' larger than the tables in Kelley (2007) in many conditions, and up to
#' three larger with \code{assurance}; for the worked example they are 100
#' and 142, where this function gives 101 and 142.
#'
#' @return
#' A \code{data.frame} with the columns \code{term} and \code{value} and
#' one row, \code{necessary_N}, the planned sample size, with or without
#' \code{assurance}.
#'
#' @references
#' Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
#'   variation with minimum risk: A sequential method for minimizing
#'   sampling error and study cost.
#'   \emph{Multivariate Behavioral Research, 51}(5), 627--648.
#'   \doi{10.1080/00273171.2016.1203279}
#'
#' Kelley, K. (2007). Sample size planning for the coefficient of variation
#'   from the accuracy in parameter estimation approach.
#'   \emph{Behavior Research Methods, 39}(4), 755--766.
#'   \doi{10.3758/BF03192966}
#'
#' Sokal, R. R., & Braumann, C. A. (1980). Significance tests for
#'   coefficients of variation and variability profiles.
#'   \emph{Systematic Zoology, 29}(1), 50--66. \doi{10.1093/sysbio/29.1.50}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on confidence intervals for
#'   effect sizes.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{ss_aipe_cv_sensitivity}}, \code{\link{cv}}
#'
#' @examples
#' # Suppose one wishes to have a confidence interval with an expected width of .10
#' # for a 99% confidence interval when the population coefficient of variation is .10.
#' ss_aipe_cv(true_cv = .1, width = .1, conf_level = .99)
#'
#' # The same planning problem parameterized by the population mean and standard
#' # deviation: true_mean = 10 and true_sd = 1 imply the same coefficient of variation, .10.
#' ss_aipe_cv(true_mean = 10, true_sd = 1, width = .1, conf_level = .99)
#'
#' # Ensuring that the confidence interval will be sufficiently narrow with a 99%
#' # certainty for the situation above.
#' ss_aipe_cv(true_cv = .1, width = .1, conf_level = .99, assurance = .99)
#'
#' @keywords design htest
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_cv <- function(true_cv = NULL, width = NULL, conf_level = .95, assurance = NULL, true_mean = NULL, true_sd = NULL, alpha_lower = NULL, alpha_upper = NULL) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  .dmar_check_numbers("ss_aipe_cv", "true_cv")
  # The coefficient of variation can be supplied directly as 'true_cv' or
  # through the population mean and standard deviation, in which case
  # true_cv = true_sd / true_mean. Exactly one of the two parameterizations is allowed;
  # accepting both would let them disagree silently.
  if (!is.null(true_cv) && (!is.null(true_mean) || !is.null(true_sd)))
    .dmar_refuse_argument(
      "Specify either 'true_cv' or both 'true_mean' and 'true_sd'; you cannot mix the two parameterizations.",
      arg = c("true_cv", "true_mean", "true_sd"))
  if (is.null(true_cv)) {
    if (is.null(true_mean) || is.null(true_sd))
      .dmar_refuse_argument(
        "Specify either 'true_cv' or both 'true_mean' and 'true_sd'.",
        arg = c("true_cv", "true_mean", "true_sd"))
    # Each must be positive on its own: a negative mean with a negative
    # standard deviation has a positive ratio and no population behind it.
    for (nm in c("true_mean", "true_sd")) {
      x <- get(nm)
      if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
        .dmar_refuse_argument(c(
          "'", nm, "' must be a single finite positive number; ",
          .dmar_given(x), ". The plan is for a population with a positive ",
          "mean and a positive standard deviation."),
          arg = nm, value = x)
    }
    true_cv <- true_sd / true_mean
  }

  if (!is.numeric(true_cv) || length(true_cv) != 1L || !is.finite(true_cv) ||
      true_cv <= 0)
    .dmar_refuse_argument(c(
      "'true_cv' must be a single finite positive number; ",
      .dmar_given(true_cv), ". The coefficient of variation the sample ",
      "size is planned for is that of a population with a positive mean ",
      "(see Chattopadhyay & Kelley, 2016)."),
      arg = "true_cv", value = true_cv)

  # A width that is not positive or not finite is never reached by the
  # interval, and the search would run to its cap.
  if (is.null(width) || !is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0)
    .dmar_refuse_argument(
      "'width' must be a single finite positive number.", arg = "width")

  if (!is.null(conf_level) && (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1))
    .dmar_refuse_argument(
      "'conf_level' must be NULL or a single number strictly between 0 and 1.",
      arg = "conf_level")

  # Without a level each tail area is given, and each is one number that
  # is at least 0 and below 1. A tail left out (NULL), NA, or a vector is
  # refused here, in the argument's name.
  if (is.null(conf_level)) {
    for (nm in c("alpha_lower", "alpha_upper")) {
      a <- get(nm)
      if (!is.numeric(a) || length(a) != 1L || !is.finite(a) || a < 0 || a >= 1)
        .dmar_refuse_argument(c(
          "'", nm, "' must be a single number that is at least 0 and below ",
          "1 when 'conf_level' is NULL; ", .dmar_given(a), ". Give both ",
          "tail areas, 'alpha_lower' and 'alpha_upper', or a 'conf_level'."),
          arg = nm, value = a)
    }
  }

  if (!is.null(conf_level)) {
    if (!is.null(alpha_lower) || !is.null(alpha_upper)) .dmar_refuse_argument(c(
      "'alpha_lower' and 'alpha_upper' are used only when 'conf_level' is ",
      "NULL, and 'conf_level' is ", format(conf_level), " (its default is ",
      ".95). To plan with the two tail areas, give them with ",
      "conf_level = NULL; otherwise leave both out."),
      arg = c("conf_level", "alpha_lower", "alpha_upper"))
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }
  # The floor on a tail area (R/dmar_conditions.R), which ci_nc_t(), whose
  # limits ci_cv() takes, also checks; here a level that leaves too small a
  # tail is refused in its own name. The refusal does not offer a tail of
  # 0: alpha_upper = 0 is refused below, and alpha_lower = 0 only suits
  # some widths.
  .dmar_check_tail_floor(alpha_lower, alpha_upper, conf_level = conf_level,
                         zero_ok = FALSE)

  # With alpha_upper = 0 the interval of ci_cv() is a one-sided lower
  # confidence bound: its upper limit is infinite at every N, so no N makes
  # the width finite. alpha_lower = 0 is an upper bound: its lower limit is
  # 0, and its width, the bound itself, comes down toward the coefficient
  # of variation as N grows, so a width above true_cv can be planned.
  if (alpha_upper == 0)
    .dmar_refuse_argument(c(
      "'alpha_upper' is 0, so the interval is a one-sided lower confidence ",
      "bound for the coefficient of variation: its upper limit is infinite ",
      "at every sample size, and no sample size makes it as narrow as ",
      "'width'. Give 'alpha_upper' a value above 0. (With 'alpha_lower' 0 ",
      "the interval is an upper confidence bound: its lower limit is 0, and ",
      "its width is the bound itself, which comes down toward the ",
      "coefficient of variation as the sample size grows.)"),
      arg = "alpha_upper", value = alpha_upper)

  # With alpha_lower = 0 and alpha_upper below .5 the width the search
  # compares with `width` is the upper bound at the planning value of
  # Equation 3, and it stays above true_cv at every N the search can
  # visit, coming down toward it as N grows: over true_cv from .001 to 5,
  # alpha_upper from .001 to .499, and N from 4 to 100,004 (every N to 400),
  # the smallest bound was 1.000006 times true_cv (work/round6/evidence/
  # g1_final_fixes/fix1/26). This is a finding and not a consequence of
  # the definition, since the value of Equation 3 is below true_cv when
  # true_cv is below sqrt(N (N - 3)) / (2 (N - 1)). A width at or below
  # true_cv is therefore never reached, and the search would run its
  # 100,000 steps before stopping. With alpha_upper above .5 the bound can
  # fall below true_cv, and the search runs.
  if (alpha_lower == 0 && alpha_upper < .5 && width <= true_cv)
    .dmar_refuse_argument(c(
      "'alpha_lower' is 0, so the interval is a one-sided upper confidence ",
      "bound for the coefficient of variation, and its width is the bound ",
      "itself. With 'alpha_upper' below .5 that bound, at the expected ",
      "sample coefficient, is above the coefficient of variation at every ",
      "sample size (it comes down toward it as the sample size grows), so ",
      "a 'width' of ", format(width),
      ", at or below the coefficient of variation of ", format(true_cv),
      ", is never reached. Ask for a 'width' above the coefficient of ",
      "variation, or give 'alpha_lower' a value above 0."),
      arg = c("width", "alpha_lower"), value = width)

  # One check and one classed refusal (dmar_assurance_invalid). Any
  # probability strictly inside (0, 1) is a request the published procedure
  # can be run for.
  .dmar_check_assurance(assurance)

  # The published assurance step was written for a two-sided interval. For
  # an upper confidence bound (alpha_lower = 0) it does not give the
  # probability asked for: with true_cv = .2, width = .25, and
  # alpha_upper = .05 the step gives 33, 91, and 219 for assurances of
  # .50, .80, and .90, where the probability that the bound is no larger
  # than the width is .470, .878, and .998 and the smallest sufficient
  # sample sizes are 36, 74, and 98. The combination is refused, after the
  # check of `assurance` itself and after the refusal of a width that no
  # bound reaches.
  if (!is.null(assurance) && alpha_lower == 0)
    .dmar_refuse_argument(c(
      "'assurance' cannot be given with 'alpha_lower' = 0. With ",
      "'alpha_lower' = 0 the interval is a one-sided upper confidence ",
      "bound, and the assurance step of Kelley (2007), which was written ",
      "for a two-sided interval, does not give the requested probability ",
      "for it: the sample size it returns can be too small or far too ",
      "large. Plan without 'assurance'."),
      arg = c("assurance", "alpha_lower"), value = assurance)

  # The planning value at each N is the value of Kelley's (2007) Equation
  # 3 for the sample coefficient of variation, not the nearly unbiased
  # estimate of Equation 4 applied to the population value.
  minimal_N <- 4

  # The interval at the planning value for a candidate N. At a small N
  # with a large coefficient its upper limit is infinite (at 99 percent
  # with a coefficient of .45, N = 4 is such a size); an infinite width is
  # never narrow enough, so the search goes on.
  interval_width <- function(N) {
    ci <- ci_cv(cv = .cv_expected_value(true_cv, N), N = N, alpha_lower = alpha_lower, alpha_upper = alpha_upper, conf_level = NULL)
    ci$value[ci$term == "ci_upper"] - ci$value[ci$term == "ci_lower"]
  }

  dif <- interval_width(minimal_N) - width
  N_0 <- minimal_N
  # N starts at the minimal sample size, so that it is defined when the
  # requested width is already attained there (dif <= 0) and the loop
  # never runs.
  N <- minimal_N

  # An explicit cap on the increment-by-one search so a pathological input
  # cannot loop without terminating.
  .iter <- 0L
  max_iter <- 100000L
  while (dif > 0) {
    .iter <- .iter + 1L
    if (.iter > max_iter)
      .dmar_refuse(c(
        "The sample size search did not converge within ", max_iter, " iterations; the specified target may be infeasible."),
        "computation_failed")
    N <- N_0 + 1
    dif <- interval_width(N) - width
    N_0 <- N
  }

  # The assurance step, the computation that produced the tables of Kelley
  # (2007). At the expected-width sample size N it finds kappa_gamma, the
  # value the sample coefficient exceeds with probability 1 - assurance,
  # takes the value of Equation 3 at kappa_gamma and N as the planning
  # value, and runs the expected-width search for that value, which
  # applies Equation 3 to it again at each candidate size. The article's
  # text (p. 759) describes the search as run for kappa_gamma itself.
  if (!is.null(assurance)) {
    # The noncentral t value that sqrt(N) / k falls below with probability
    # 1 - assurance. The quantile is exact (.nc_t_quantile(),
    # R/nc_t_internals.R); R's qt() with ncp inverts pt(), a normal
    # approximation above a noncentrality of 37.62, which
    # sqrt(N) / true_cv commonly passes.
    beyond_CV_NCP <- .nc_t_quantile(1 - assurance, df = N - 1, ncp = sqrt(N) / true_cv)

    # That value as a coefficient of variation, kappa_gamma.
    Lim_for_Certainty <- sqrt(N) / beyond_CV_NCP

    N_gamma <- ss_aipe_cv(true_cv = .cv_expected_value(Lim_for_Certainty, N), width = width, alpha_lower = alpha_lower, alpha_upper = alpha_upper, conf_level = NULL, assurance = NULL)[1, 2]
    return(.as_dmar_tbl(data.frame(term = "necessary_N", value = N_gamma), conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }

  .as_dmar_tbl(data.frame(term = "necessary_N", value = N), conf_level = conf_level, subclass = "dmar_ss_aipe")
}


# The planning value for the sample coefficient of variation of normal
# data, given the population value 'kappa' and the sample size 'N': Kelley
# (2007), Equation 3, the approximation of Sokal and Braumann (1980, p. 51,
# Expression 3) to its expected value. The approximation stops at the
# third-order term and is above the expected value at small N (by 3.1
# percent at N = 5 for kappa = .25); the planner keeps it because the
# published tables were computed with it. MBESS through 3.3.3 returned
# this from cv(C.of.V, N, unbiased = TRUE), which ss.aipe.cv() called;
# cv() later changed to the (1 + 1/(4N)) multiplier of Equation 4 and the
# planner followed it without being edited.
.cv_expected_value <- function(kappa, N) {
  kappa * (1 - 1 / (4 * (N - 1)) + kappa^2 / N + 1 / (2 * (N - 1)^2))
}
