# Large-sample variance of the coefficient of variation under normality.
#' Asymptotic Variance of the Coefficient of Variation
#'
#' Computes the large-sample variance of the sample coefficient of
#' variation \eqn{k = s / M}, the sample standard deviation divided by the
#' sample mean, for a random sample from a normal population:
#' \eqn{\kappa^2 (1/2 + \kappa^2) / (N - 1)}. It takes the population
#' coefficient of variation and the sample size and returns the variance.
#' Companion to \code{\link{ci_cv}} and \code{\link{ss_aipe_cv}}.
#'
#' @param true_cv Population coefficient of variation
#'   \eqn{\kappa = \sigma / \mu}. A single finite number greater than 0.
#'   For \eqn{\kappa} above about 0.4 the value returned is too small, and
#'   more so the larger \eqn{\kappa} is (see Details); for an interval use
#'   \code{\link{ci_cv}}, which does not use this variance.
#' @param N Total sample size. A single finite number of at least 3.
#'
#' @return A \code{data.frame} with the columns \code{term} and
#'   \code{value} and one row, \code{var_cv_normal}: the first-order
#'   variance of the sample coefficient of variation under normality (see
#'   Details).
#'
#' @details
#' \strong{Estimand.} The variance, over repeated samples of \eqn{N}
#' independent observations from a normal population with mean
#' \eqn{\mu > 0} and standard deviation \eqn{\sigma}, of the sample
#' coefficient of variation \eqn{k = s / M}, evaluated at the population
#' value \eqn{\kappa = \sigma / \mu} given as \code{true_cv}. Strictly, the
#' sample coefficient has no finite variance, because the sample mean can
#' be arbitrarily close to zero; the expression below describes its spread
#' when that event has negligible probability, that is, when
#' \eqn{\sqrt{N} / \kappa} is large. Substituting a sample coefficient for
#' \code{true_cv} gives a plug-in estimate of the variance, which carries
#' the sampling error of the coefficient it was computed from.
#'
#' \strong{The variance.} The first-order (delta method) variance is
#' \deqn{\mathrm{Var}(k) \;\approx\;
#'   \frac{\kappa^2}{N - 1}
#'   \left(\frac{1}{2} + \kappa^2\right).}
#' It follows from \eqn{\mathrm{Var}(s) \approx \sigma^2 / \{2 (N - 1)\}},
#' \eqn{\mathrm{Var}(M) = \sigma^2 / N}, and the independence of \eqn{s}
#' and \eqn{M} under normality, which give
#' \eqn{\mathrm{Var}(k) \approx \kappa^2 / \{2 (N - 1)\} + \kappa^4 / N}.
#' That expression is Expression 5 of Sokal and Braumann (1980, p. 52)
#' without its third-order term. The function uses \eqn{N - 1} under both
#' terms, which agrees with it to first order in \eqn{1/N} and is larger by
#' \eqn{\kappa^4 / \{N (N - 1)\}}. The form usually given, their
#' Expression 7, has \eqn{N} under both terms:
#' \eqn{\kappa^2 (1/2 + \kappa^2) / N}.
#'
#' \strong{Accuracy under normality.} The comparisons in this paragraph
#' are with the variance of \eqn{k} when the sample mean is restricted to
#' lie within 8 standard errors of \eqn{\mu} and to be at least
#' \eqn{\mu / 50}. The restriction excludes a probability below
#' \eqn{10^{-9}} in every case reported, and the figures given are the
#' same for any lower bound from \eqn{\mu / 1000} to \eqn{\mu / 5}. For
#' \eqn{\kappa \le 0.3} and \eqn{N \ge 10} the value is too large. The
#' excess is largest for a small \eqn{\kappa}, where it is that of
#' \eqn{\sigma^2 / \{2 (N - 1)\}} over the exact variance of \eqn{s}: to
#' one decimal place, 3.0 percent at \eqn{N = 10}, 1.4 percent at
#' \eqn{N = 20}, 0.9 percent at \eqn{N = 30}, and 0.5 percent at
#' \eqn{N = 50}. For \eqn{\kappa} above about 0.4 the value is too small,
#' because of the terms of higher order that a first-order variance omits:
#' by 5.7 percent at \eqn{\kappa = 0.5} and \eqn{N = 10}, and by 1.7
#' percent at \eqn{\kappa = 0.5} and \eqn{N = 30}.
#'
#' \strong{Normality.} The value is for a normal population. For any
#' population with a finite fourth moment, the large-sample variance of
#' \eqn{k} is \eqn{\xi^2 / N}, with
#' \deqn{\xi^2 \;=\; \frac{\kappa^2}{2} + \kappa^4 - \kappa^3 \gamma_1
#'   + \frac{\kappa^2 \gamma_2}{4},}
#' where \eqn{\gamma_1} is the skewness and \eqn{\gamma_2} the excess
#' kurtosis of the population (Kelley, Darku, & Chattopadhyay, 2018,
#' Equation 30, which is written in central moments). Under normality both
#' are 0 and \eqn{\xi^2 = \kappa^2 (1/2 + \kappa^2)}. For other populations
#' the normal-theory value can be far off in either direction:
#' \eqn{\kappa^2 (1/2 + \kappa^2)} is 39 percent larger than \eqn{\xi^2}
#' for a gamma population with \eqn{\kappa = 0.8}, and 56 percent smaller
#' for a Laplace (double exponential) population with \eqn{\kappa = 0.3}.
#' These compare the limiting values, so a larger sample does not remove
#' the difference. \code{\link{mr_cv}} estimates \eqn{\xi^2} from the
#' scores, without assuming a normal population.
#'
#' \strong{Use.} Neither \code{\link{ci_cv}} nor \code{\link{ss_aipe_cv}}
#' uses this variance: the interval inverts the noncentral \emph{t}
#' distribution of \eqn{\sqrt{N} / k}, and the planner (Kelley, 2007)
#' chooses the sample size for that interval. The value serves as a quick
#' description of the sampling variability of \eqn{k} and for comparison
#' with a simulated variance. McKay (1932) and Vangel (1996) approximate
#' the distribution of a function of \eqn{k} by a chi square distribution,
#' from which confidence limits follow; the variance here is not taken
#' from those approximations.
#'
#' @references
#' Kelley, K. (2007). Sample size planning for the coefficient of variation
#'   from the accuracy in parameter estimation approach.
#'   \emph{Behavior Research Methods, 39}(4), 755--766.
#'   \doi{10.3758/BF03192966}
#'
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#'   parameter estimation for a general class of effect sizes: A
#'   sequential approach. \emph{Psychological Methods, 23}(2), 226--243.
#'   \doi{10.1037/met0000127}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' McKay, A. T. (1932). Distribution of the coefficient of variation and
#'   the extended \emph{t} distribution. \emph{Journal of the Royal
#'   Statistical Society, 95}(4), 695--698.
#'
#' Sokal, R. R., & Braumann, C. A. (1980). Significance tests for
#'   coefficients of variation and variability profiles.
#'   \emph{Systematic Zoology, 29}(1), 50--66.
#'   \doi{10.1093/sysbio/29.1.50}
#'
#' Vangel, M. G. (1996). Confidence intervals for a normal coefficient
#'   of variation. \emph{The American Statistician, 50}(1), 21--26.
#'   \doi{10.1080/00031305.1996.10473537}
#'
#' @seealso \code{\link{cv}}, \code{\link{ci_cv}}, \code{\link{ss_aipe_cv}},
#'   \code{\link{ss_aipe_cv_sensitivity}}, \code{\link{mr_cv}}
#'
#' @examples
#' # The variance of the sample coefficient of variation in samples of 30
#' # from a normal population whose coefficient is .20. The variance is
#' # small, so it is printed to six decimal places.
#' v <- var_cv(true_cv = 0.20, N = 30)
#' print(v, digits = 6)
#'
#' # Its square root is the standard error of the sample coefficient, shown
#' # here beside the standard deviation of 2000 simulated coefficients.
#' set.seed(113)
#' k <- replicate(2000, {
#'   y <- rnorm(30, mean = 10, sd = 2)
#'   sd(y) / mean(y)
#' })
#' c(normal_theory = sqrt(v$value), simulated = sd(k))
#'
#' # The same comparison for a skewed population, a gamma distribution
#' # with a coefficient of .80 and samples of 200: the normal-theory
#' # standard error is too large.
#' k_gamma <- replicate(2000, {
#'   y <- rgamma(200, shape = 1 / 0.80^2)
#'   sd(y) / mean(y)
#' })
#' c(normal_theory = sqrt(var_cv(true_cv = 0.80, N = 200)$value),
#'   simulated = sd(k_gamma))
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family variance utilities
#'
#' @export

var_cv <- function(true_cv, N) {
  if (!is.numeric(true_cv) || length(true_cv) != 1L || !is.finite(true_cv) ||
      true_cv <= 0)
    .dmar_refuse_argument(c(
      "'true_cv' must be a single finite number greater than 0, the ",
      "population coefficient of variation; ", .dmar_given(true_cv), "."),
      arg = "true_cv", value = true_cv)
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 3)
    .dmar_refuse_argument(c(
      "'N' must be a single finite sample size of at least 3; ",
      .dmar_given(N), "."),
      arg = "N", value = N)

  out <- data.frame(
    term  = "var_cv_normal",
    value = (true_cv^2 / (N - 1)) * (0.5 + true_cv^2),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out)
}
