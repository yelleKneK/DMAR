#' Coefficient of Variation (Sample Value or Nearly Unbiased Estimate)
#'
#' Computes the sample coefficient of variation \eqn{k = s / M}, the
#' sample standard deviation divided by the sample mean, or, optionally,
#' the nearly unbiased estimate of Kelley (2007, Equation 4) under
#' normality. Either supply a precomputed \code{cv} or supply the raw
#' \code{mean} and \code{sd}; with \code{unbiased = TRUE} the value is
#' multiplied by \eqn{(1 + 1/(4 N))}. The (biased) sample
#' coefficient of variation (the default, \code{unbiased = FALSE}) is the
#' form usually reported. To accompany it with a confidence interval, use
#' \code{\link{ci_cv}}.
#'
#' @param cv The sample coefficient of variation, \eqn{k = s / M}, a single
#'   finite number. Optional; if supplied, \code{mean} and \code{sd} must
#'   not be.
#' @param mean Sample mean. A single finite number other than 0. A negative
#'   mean gives a negative coefficient; \code{\link{ci_cv}} requires a
#'   positive one.
#' @param sd Sample standard deviation, using \eqn{N - 1} in the variance
#'   denominator. A single finite number that is not negative.
#' @param N Sample size. Required when \code{unbiased = TRUE}, where it
#'   must exceed 1.
#' @param unbiased Logical. If \code{TRUE}, multiplies the sample value by
#'   \eqn{(1 + 1/(4 N))}. The result is nearly, not exactly, unbiased (see
#'   Details). Default \code{FALSE}.
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} and \code{value}.
#' The \code{term} value is \code{"cv"} and \code{value} is either the
#' sample coefficient (default) or the nearly unbiased estimate
#' (when \code{unbiased = TRUE}).
#'
#' @details
#' The sample coefficient of variation \eqn{k = s / M} is the form
#' usually reported in applied work. It is what this function returns by
#' default (\code{unbiased = FALSE}). A point estimate is most informative
#' when paired with a confidence interval; \code{\link{ci_cv}} computes one
#' for the population value \eqn{\kappa}.
#'
#' Under normality the sample coefficient is biased. Its expected value is
#' approximately \eqn{\kappa [1 - 1/(4 (N - 1)) + \kappa^2/N]}, the first
#' and second order terms of Expression 3 of Sokal and Braumann (1980,
#' p. 51), which Kelley (2007) reproduces as his Equation 3. The term
#' \eqn{-1/(4 (N - 1))} comes from the sample standard deviation and the
#' term \eqn{\kappa^2/N} from the sample mean in the denominator, so
#' \eqn{k} underestimates \eqn{\kappa} on average when \eqn{\kappa} is
#' below about one half and overestimates it when \eqn{\kappa} is larger.
#' Multiplying \eqn{k} by \eqn{(1 + 1/(4 N))}, a correction that Sokal and
#' Braumann (1980) and Kelley (2007) credit to Haldane (1955), offsets the
#' first term and not the second. The relative bias that remains is about
#' \eqn{\kappa^2/N - 9/(32 N^2)}; the second of these terms is derived
#' from the expected value of \eqn{s} under normality and is not taken
#' from those sources. Kelley (2007, Equation 4) calls the result a nearly
#' unbiased estimate and recommends it when a point estimate is reported;
#' the sample value \eqn{k} is the one used to form the confidence
#' interval.
#'
#' At \eqn{\kappa = 0.1} and \eqn{N = 10} the relative bias is about
#' \eqn{-2.6} percent for \eqn{k} and \eqn{-0.2} percent for the corrected
#' value; at \eqn{\kappa = 0.3} and \eqn{N = 10} it is about \eqn{-1.8}
#' percent and \eqn{+0.6} percent. With 10 cases the two are about
#' equally biased, in opposite directions, at a coefficient near 0.38
#' (near 0.35 in large samples, and above 0.4 with 5 cases or fewer),
#' and for a larger \eqn{\kappa} the corrected value is the more biased of
#' the two. These statements describe a mean that stays away from 0: under
#' an unrestricted normal distribution the expectation of
#' \eqn{1/M} does not exist as an ordinary integral.
#'
#' For confidence intervals on \eqn{\kappa}, \code{\link{ci_cv}} inverts
#' the noncentral \emph{t} distribution of \eqn{\sqrt{N}/k} (Johnson &
#' Welch, 1940; McKay, 1932; Kelley, 2007). McKay (1932) and Vangel (1996)
#' also give approximations based on the chi square distribution, which
#' \code{ci_cv} does not use. The first-order variance of \eqn{k} under
#' normality, \eqn{\kappa^2 (1/2 + \kappa^2) / (N - 1)}, is in
#' \code{\link{var_cv}}, whose page relates it to Expression 5 of Sokal and
#' Braumann (1980, p. 52).
#'
#' @examples
#' # 1. Point estimate from raw mean and SD.
#' cv(mean = 100, sd = 15)
#'
#' # 2. Nearly unbiased estimate at N = 50; the multiplier is close to 1
#' #    at this sample size.
#' cv(mean = 100, sd = 15, N = 50, unbiased = TRUE)
#'
#' # 3. The same at N = 10, where the multiplier is larger.
#' cv(cv = .15, N = 10, unbiased = TRUE)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_cv}}, \code{\link{var_cv}},
#'   \code{\link{ss_aipe_cv}}
#'
#' @references
#' Haldane, J. B. S. (1955). The measurement of variation.
#'   \emph{Evolution, 9}(4), 484.
#'   \doi{10.1111/j.1558-5646.1955.tb01559.x}
#'
#' Johnson, N. L., & Welch, B. L. (1940). Applications of the non-central
#'   \emph{t}-distribution. \emph{Biometrika, 31}(3--4), 362--389.
#'   \doi{10.1093/biomet/31.3-4.362}
#'
#' Kelley, K. (2007). Sample size planning for the coefficient of variation
#'   from the accuracy in parameter estimation approach.
#'   \emph{Behavior Research Methods, 39}(4), 755--766.
#'   \doi{10.3758/BF03192966}
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
#' @keywords design htest
#'
#' @export

cv <- function(cv = NULL, mean = NULL, sd = NULL, N = NULL, unbiased = FALSE) {
  if (!is.logical(unbiased) || length(unbiased) != 1L || is.na(unbiased))
    .dmar_refuse_argument("'unbiased' must be TRUE or FALSE.", arg = "unbiased")
  single_finite <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)

  if (!is.null(cv)) {
    if (!is.null(mean) && !is.null(sd))
      .dmar_refuse_argument(
        "Because 'mean' and 'sd' were specified, do not specify 'cv'.",
        arg = c("mean", "sd", "cv"))
    if (!is.null(mean) || !is.null(sd))
      .dmar_refuse_argument(
        "Because 'cv' was specified, do not specify 'mean' or 'sd'.",
        arg = c("cv", "mean", "sd"))
    if (!single_finite(cv))
      .dmar_refuse_argument("'cv' must be a single finite number.", arg = "cv")
    k <- cv
  } else if (!is.null(mean) && !is.null(sd)) {
    if (!single_finite(mean) || mean == 0)
      .dmar_refuse_argument(c(
        "'mean' must be a single finite number other than 0; the ",
        "coefficient of variation is undefined at a mean of 0."), arg = "mean")
    if (!single_finite(sd) || sd < 0)
      .dmar_refuse_argument(
        "'sd' must be a single finite number that is not negative.", arg = "sd")
    k <- sd / mean
  } else {
    .dmar_refuse_argument(
      "You must specify either 'cv' or both 'mean' and 'sd'.",
      arg = c("cv", "mean", "sd"))
  }

  if (unbiased) {
    if (is.null(N)) .dmar_refuse_argument(
      "You must specify 'N' when 'unbiased = TRUE'.", arg = "N")
    if (!single_finite(N) || N <= 1)
      .dmar_refuse_argument(
        "'N' must be a single finite sample size greater than 1.", arg = "N")
    k <- k * (1 + 1 / (4 * N))
  }

  return(.as_dmar_tbl(data.frame(term = 'cv', value = k)))
}
