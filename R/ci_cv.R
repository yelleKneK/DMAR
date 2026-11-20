#' Confidence Interval for the Coefficient of Variation
#'
#' Computes the confidence interval for the population coefficient of
#' variation, the standard deviation relative to the mean, from the
#' noncentral \emph{t} distribution, so that variability can be reported
#' on a scale that is free of the measurement units. The interval is exact
#' under normality and is not robust to departures from it (see Details).
#'
#' @param cv Sample coefficient of variation, the sample standard
#'   deviation divided by the sample mean, given with \code{N}. Must be
#'   positive.
#' @param mean Sample mean; must be positive. It is given with \code{sd}
#'   and \code{N}.
#' @param sd Sample standard deviation (square root of the unbiased estimate of the variance). Must be positive.
#' @param N Total sample size, given with \code{cv} or with \code{mean}
#'   and \code{sd}. Must exceed 1; it need not be a whole number.
#' @param data Numeric vector of at least two finite scores, given in
#'   place of the summary statistics: the coefficient of variation and the
#'   sample size are taken from it.
#' @param conf_level Desired confidence level (1-Type I error rate). A level above .999999998, whose tail areas
#'   are below 1e-9, is refused
#' @param alpha_lower The lower tail area. An omitted confidence level yields to supplied tail areas; an explicit level must agree with them.
#' @param alpha_upper The upper tail area, supplied together with \code{alpha_lower}.
#'   The table records the coverage, \code{1 - alpha_lower - alpha_upper},
#'   as its level and, when the two tail areas differ, the areas
#'   themselves, which the line beneath the printed table,
#'   \code{\link{results_sentence}}, and the \code{alpha_lower} and
#'   \code{alpha_upper} columns of \code{tidy()} state; with a tail area
#'   of 0 the interval is a one-sided confidence bound. A nonzero tail
#'   area below 1e-9 is refused (see "Small tail areas" in
#'   \code{\link{ci_nc_t}}).
#'
#' @details
#' \strong{The interval.} For a random sample of \eqn{N} scores from a
#' normal population with a positive mean, \eqn{\sqrt{N} / k}, where
#' \eqn{k = s / M} is the sample coefficient of variation, has a
#' noncentral \emph{t} distribution on \eqn{N - 1} degrees of freedom
#' whose noncentrality parameter is \eqn{\sqrt{N} / \kappa}, with
#' \eqn{\kappa = \sigma / \mu} the population coefficient of variation
#' (Johnson & Welch, 1940; Kelley, 2007, Equations 5 to 8; McKay, 1932).
#' The function finds the confidence limits for the noncentrality
#' parameter (\code{\link{ci_nc_t}}) and divides \eqn{\sqrt{N}} by each:
#' the lower limit for \eqn{\kappa} comes from the upper limit for the
#' noncentrality parameter, and the upper limit from the lower.
#'
#' \strong{Normality.} The interval is exact when the scores are a random
#' sample from a normal distribution with a positive mean: its coverage is
#' then the stated level at every sample size. It is not robust to
#' departures from normality, and a larger sample does not repair it,
#' because the distribution that is inverted is that of a normal sample
#' whatever its size. In simulations of a nominal 95 percent interval,
#' with 10,000 samples in each condition (population coefficients of
#' variation of .3 and .5, sample sizes of 30 and 200), coverage was .81
#' to .85 for Laplace data, .90 to .94 for lognormal data, and .96 to .97
#' for gamma data, where it was .945 to .951 for normal data, the stated
#' level within the error of the simulation. For the Laplace and the
#' lognormal data it was lower at 200 than at 30. A
#' large-sample interval that does not assume normality follows from the
#' asymptotic normal distribution of \eqn{k} (Kelley, Darku, &
#' Chattopadhyay, 2018, Equations 29 to 33). It needs a large sample:
#' with 400,000 samples in each of the same conditions it covered .806 to
#' .905 at \eqn{N = 30} and .894 to .942 at \eqn{N = 200}, below the
#' stated level in every condition, normal data included. The function
#' does not compute it.
#'
#' \strong{An infinite upper limit.} The upper limit is infinite when
#' \eqn{\sqrt{N} / k} is at most the \code{1 - alpha_upper} quantile of
#' the central \emph{t} distribution on \eqn{N - 1} degrees of freedom,
#' which happens when the sample coefficient is large for the sample
#' size. The lower confidence limit of the noncentrality parameter is
#' then not positive: at that tail area the data do not rule out a
#' population mean of zero, and no finite coefficient of variation bounds
#' the interval above. The infinite limit is part of the exact interval.
#' Over repeated samples from a normal population with a positive mean,
#' the interval covers the population coefficient with the stated
#' probability, counting the intervals that are bounded above and those
#' that are not. In 20,000 simulated samples with a population
#' coefficient of .8 and \eqn{N = 6}, the upper limit of the 95 percent
#' interval was infinite in 31 percent and coverage was .952; with a
#' coefficient of .5 and \eqn{N = 5} it was infinite in 9 percent and
#' coverage was .948. The stated level describes all of the intervals
#' together and neither kind alone (computed from the noncentral \emph{t}
#' distribution, coverage is .92 and .73 among the intervals with an
#' infinite limit and .96 and .97 among the others).
#' Nothing is signaled for an infinite limit, and the \code{prob_less}
#' and \code{prob_greater} columns hold the requested tail areas there as
#' at a finite limit. A sample whose mean is not positive has no interval
#' and is refused; with a population coefficient of .5 or less and five
#' or more scores the probability of such a sample is below .00001.
#' Where such samples are not rare, coverage among the samples that have
#' an interval is somewhat above the stated level: it is .962 for a 95
#' percent interval with a coefficient of 1 and \eqn{N = 5}, where 1.3
#' percent of samples have a mean that is not positive. Of the 20,000
#' samples with a coefficient of .8 above, 25 had such a mean, and the
#' figures given describe the other 19,975.
#'
#' @return
#' A 4-row \code{data.frame} with columns \code{term}, \code{value},
#' \code{prob_less}, and \code{prob_greater}. The rows are ordered so the two
#' point estimates sit between the confidence limits:
#' \code{"ci_lower"} (lower confidence limit on the coefficient of
#' variation), \code{"cv"} (the sample coefficient of variation),
#' \code{"cv_unbiased"} (the nearly unbiased estimate,
#' \eqn{(1 + 1/(4 N))} times the sample value; see \code{\link{cv}}), and
#' \code{"ci_upper"} (upper confidence limit). The \code{prob_less}
#' and \code{prob_greater} columns hold the requested tail areas:
#' \code{alpha_lower} and \code{1 - alpha_lower} at the lower limit, and
#' \code{1 - alpha_upper} and \code{alpha_upper} at the upper limit, also
#' when the upper limit is infinite. They are \code{NA} for the
#' point-estimate rows.
#'
#' @references
#' Johnson, N. L., & Welch, B. L. (1940). Applications of the non-central
#'   \emph{t}-distribution. \emph{Biometrika, 31}(3--4), 362--389.
#'   \doi{10.1093/biomet/31.3-4.362}
#'
#' Kelley, K. (2007). Sample size planning for the coefficient of variation from the accuracy in parameter estimation approach. \emph{Behavior Research Methods, 39}(4), 755--766.
#'   \doi{10.3758/BF03192966}
#'
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#'   parameter estimation for a general class of effect sizes: A
#'   sequential approach. \emph{Psychological Methods, 23}(2), 226--243.
#'   \doi{10.1037/met0000127}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on confidence intervals for
#'   effect sizes.)
#'
#' McKay, A. T. (1932). Distribution of the coefficient of variation and
#'   the extended \emph{t} distribution. \emph{Journal of the Royal
#'   Statistical Society, 95}(4), 695--698.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cv}}
#'
#' @examples
#' # Fifteen scores from a normal population with a mean of 5 and a
#' # standard deviation of 1, whose coefficient of variation is .20.
#' set.seed(113)
#' x <- rnorm(15, mean = 5, sd = 1)
#'
#' # The interval from the scores, from their mean and standard deviation,
#' # and from the sample coefficient itself: three routes to one interval.
#' ci_cv(data = x)
#' ci_cv(mean = mean(x), sd = sd(x), N = 15)
#' ci_cv(cv = sd(x) / mean(x), N = 15)
#'
#' # Unequal tail areas, with .01 below the interval and .04 above it.
#' ci_cv(data = x, alpha_lower = .01, alpha_upper = .04)
#'
#' @keywords models htest
#'
#' @import stats
#' @family confidence intervals for effect sizes
#'
#' @export

ci_cv <- function(cv = NULL, mean = NULL, sd = NULL, N = NULL, data = NULL,
                  conf_level = 0.95, alpha_lower = NULL, alpha_upper = NULL) {
  tails_given <- !is.null(alpha_lower) || !is.null(alpha_upper)
  tails <- .ci_resolve_tails(conf_level, alpha_lower, alpha_upper,
    level_missing = missing(conf_level))
  conf_level <- tails$conf_level
  alpha_lower <- tails$alpha_lower
  alpha_upper <- tails$alpha_upper
  # The floor on a tail area (R/dmar_conditions.R). ci_nc_t() checks it
  # too, but the interval for the coefficient passes it the two tails
  # exchanged (a larger noncentrality is a smaller coefficient), so its
  # refusal would name the other tail, and a level as two tails.
  .dmar_check_tail_floor(alpha_lower, alpha_upper,
                         conf_level = if (!tails_given) conf_level)

  # Three input routes: the scores, the coefficient with the sample size,
  # or the mean and the standard deviation with the sample size. Each
  # resolves the sample coefficient k and the sample size N, and the
  # interval is formed in one place, .ci_cv_interval(). The interval
  # inverts the noncentral t distribution of sqrt(N) / k on N - 1 degrees
  # of freedom, which needs more than one score, a positive mean, and a
  # positive standard deviation.
  check_N <- function() {
    if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 1)
      .dmar_refuse_argument(
        "'N' must be a single finite sample size greater than 1.", arg = "N")
  }
  if (!is.null(data)) {
    if (!is.null(cv))
      .dmar_refuse_argument(
        "Give either 'data', or 'cv' with 'N', not both.",
        arg = c("data", "cv"))
    if (!is.null(mean))
      .dmar_refuse_argument(
        "Since 'data' is specified, do not specify the 'mean'.",
        arg = c("data", "mean"))
    if (!is.null(sd))
      .dmar_refuse_argument(
        "Since 'data' is specified, do not specify the 'sd'.",
        arg = c("data", "sd"))
    if (!is.null(N))
      .dmar_refuse_argument(
        "Since 'data' is specified, do not specify the 'N'.",
        arg = c("data", "N"))
    if (!is.numeric(data) || length(data) < 2L || any(!is.finite(data)))
      .dmar_refuse_argument(
        "'data' must be a numeric vector of at least 2 finite observations.",
        arg = "data")
    N <- length(data)
    mean_data <- base::mean(data)
    sd_data <- stats::sd(data)
    if (mean_data <= 0) .ci_cv_refuse_not_positive("mean", mean_data, "data")
    if (sd_data <= 0)
      .ci_cv_refuse_not_positive("standard deviation", sd_data, "data")
    k <- sd_data / mean_data
  } else if (!is.null(cv)) {
    if (!is.null(mean))
      .dmar_refuse_argument(
        "Since you specified the coefficient of variation ('cv') directly, do not specify the mean ('mean').",
        arg = c("cv", "mean"))
    if (!is.null(sd))
      .dmar_refuse_argument(
        "Since you specified the coefficient of variation ('cv') directly, do not specify the standard deviation ('sd').",
        arg = c("cv", "sd"))
    if (is.null(N))
      .dmar_refuse_argument(
        "'N' is required with 'cv': the interval depends on the sample size.",
        arg = "N")
    check_N()
    if (!is.numeric(cv) || length(cv) != 1L || !is.finite(cv) || cv <= 0)
      .dmar_refuse_argument(c(
        "'cv' must be a single positive number; ", .dmar_given(cv), ". ",
        "The confidence interval for the coefficient of variation requires ",
        "a positive mean and a positive standard deviation."),
        arg = "cv", value = cv)
    k <- cv
  } else {
    if (is.null(mean) && is.null(sd))
      .dmar_refuse_argument(c(
        "Give the scores as 'data', or the coefficient of variation as ",
        "'cv' with 'N', or 'mean' and 'sd' with 'N'."),
        arg = c("data", "cv", "mean", "sd"))
    if (is.null(mean))
      .dmar_refuse_argument(
        "'mean' is required with 'sd' and 'N' (or give 'cv' with 'N', or the scores as 'data').",
        arg = "mean")
    if (is.null(sd))
      .dmar_refuse_argument(
        "'sd' is required with 'mean' and 'N' (or give 'cv' with 'N', or the scores as 'data').",
        arg = "sd")
    if (is.null(N))
      .dmar_refuse_argument(
        "'N' is required with 'mean' and 'sd': the interval depends on the sample size.",
        arg = "N")
    check_N()
    if (!is.numeric(mean) || length(mean) != 1L || !is.finite(mean))
      .dmar_refuse_argument(c(
        "'mean' must be a single finite number, the sample mean; ",
        .dmar_given(mean), "."), arg = "mean", value = mean)
    if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd))
      .dmar_refuse_argument(c(
        "'sd' must be a single finite number, the sample standard ",
        "deviation; ", .dmar_given(sd), "."), arg = "sd", value = sd)
    if (mean <= 0) .ci_cv_refuse_not_positive("mean", mean, "mean")
    if (sd <= 0) .ci_cv_refuse_not_positive("standard deviation", sd, "sd")
    k <- sd / mean
  }
  .ci_cv_interval(k, N, alpha_lower, alpha_upper, conf_level)
}

# .ci_cv_refuse_not_positive(what, value, arg)
#
# The one refusal of a sample mean, or of a sample standard deviation, that
# is not positive, whichever argument it came from: `what` is "mean" or
# "standard deviation", `value` the statistic, and `arg` the argument that
# held it ('mean', 'sd', or 'data').
.ci_cv_refuse_not_positive <- function(what, value, arg) {
  .dmar_refuse_argument(c(
    "The sample ", what, " is ", format(value), ", and it must be ",
    "positive: the confidence interval for the coefficient of variation ",
    "requires a positive mean and a positive standard deviation."),
    arg = arg, value = value)
}

# .ci_cv_interval(k, N, alpha_lower, alpha_upper, conf_level)
#
# The interval of ci_cv() for a sample coefficient of variation k from N
# scores, with the tail areas alpha_lower and alpha_upper already
# resolved. sqrt(N) / k has a noncentral t distribution on N - 1 degrees
# of freedom with noncentrality parameter sqrt(N) / kappa, so the limits
# for kappa are sqrt(N) over the limits for the noncentrality parameter,
# exchanged: the tail area below the interval for kappa is the one above
# the interval for the noncentrality parameter. The limits of ci_nc_t()
# are read by term name.
#
# A lower noncentrality limit that is not positive gives an infinite
# upper limit for kappa (with alpha_upper = 0 it is -Inf). That limit is
# part of the exact interval, so the tail columns hold the requested
# areas there as at a finite limit and nothing is signaled. With
# alpha_lower = 0 the upper noncentrality limit is Inf and the lower
# limit for kappa is 0.
.ci_cv_interval <- function(k, N, alpha_lower, alpha_upper, conf_level) {
  ncp <- ci_nc_t(ncp = sqrt(N) / k, df = N - 1, alpha_lower = alpha_upper,
                 alpha_upper = alpha_lower, conf_level = NULL)
  ncp_lower <- ncp$value[ncp$term == "ci_lower"]
  ncp_upper <- ncp$value[ncp$term == "ci_upper"]
  out <- data.frame(
    term = c("ci_lower", "cv", "cv_unbiased", "ci_upper"),
    value = c(sqrt(N) / ncp_upper, k, k * (1 + 1 / (4 * N)),
              if (ncp_lower <= 0) Inf else sqrt(N) / ncp_lower),
    prob_less = c(alpha_lower, NA, NA, 1 - alpha_upper),
    prob_greater = c(1 - alpha_lower, NA, NA, alpha_upper)
  )
  .as_dmar_tbl(out, conf_level = conf_level,
               tail_areas = .tail_areas(alpha_lower, alpha_upper))
}
