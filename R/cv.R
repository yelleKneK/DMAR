#' Coefficient of Variation (Biased or Unbiased Estimator)
#'
#' Computes the sample coefficient of variation \eqn{\hat\kappa = s / \bar Y}
#' or, optionally, its first-order bias-corrected counterpart under normality.
#' Either supply a precomputed \code{cv} or supply the raw \code{mean} and
#' \code{sd}; with \code{unbiased = TRUE} the value is multiplied by the
#' small-sample correction \eqn{(1 + 1/(4 N))}. The (biased) sample
#' coefficient of variation (the default, \code{unbiased = FALSE}) is the
#' form usually reported. To accompany it with a confidence interval, use
#' \code{\link{ci_cv}}.
#'
#' @param cv The sample coefficient of variation, \eqn{s/\bar Y}.
#'   Optional; if supplied, \code{mean} and \code{sd} must not be.
#' @param mean Sample mean. Numeric scalar.
#' @param sd Sample standard deviation, using \eqn{N - 1} in the variance
#'   denominator. Numeric scalar.
#' @param N Sample size. Required when \code{unbiased = TRUE}.
#' @param unbiased Logical. If \code{TRUE}, applies the first-order
#'   small-sample bias correction \eqn{(1 + 1/(4 N))} to the plug-in
#'   estimator. Default \code{FALSE}.
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} and \code{value}.
#' The \code{term} value is \code{"cv"} and \code{value} is either the
#' plug-in estimator (default) or the first-order bias-corrected estimator
#' (when \code{unbiased = TRUE}).
#'
#' @details
#' The plug-in estimator \eqn{\hat\kappa = s / \bar Y} is the workhorse
#' coefficient of variation in applied work and is the form usually
#' reported. It is what this function returns by default
#' (\code{unbiased = FALSE}). A point estimate is most informative when
#' paired with a confidence interval; \code{\link{ci_cv}} computes one for
#' \eqn{\kappa}. Under normality, however, the plug-in estimator
#' is biased downward, and the leading-order expansion is
#' \eqn{E[\hat\kappa] = \kappa (1 - 1/(4 N)) + O(N^{-2})}
#' (Sokal & Rohlf, 1995). The bias is negligible at \eqn{N} above about
#' 100 but is non-trivial in small samples, where multiplying the plug-in
#' value by \eqn{(1 + 1/(4 N))} removes the leading-order term. The
#' \code{unbiased = TRUE} option applies that correction.
#'
#' For confidence intervals on \eqn{\kappa}, the McKay (1932) noncentral
#' \emph{t} based interval is implemented in \code{\link{ci_cv}}; the
#' corresponding asymptotic variances are in \code{\link{var_cv}}. The
#' Vangel (1996) small-sample refinement of McKay's interval is a further
#' option described in the literature; it matters more than the bias
#' correction in this function when \eqn{\kappa} is larger than about
#' 0.3 (Kelley, 2007).
#'
#' @examples
#' # 1. Point estimate from raw mean and SD.
#' cv(mean = 100, sd = 15)
#'
#' # 2. Bias-corrected estimate at N = 50; the correction is small but
#' #    non-negligible at this sample size.
#' cv(mean = 100, sd = 15, N = 50, unbiased = TRUE)
#'
#' # 3. Bias correction at N = 10; the correction is larger here.
#' cv(cv = .15, N = 10, unbiased = TRUE)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_cv}}, \code{\link{var_cv}},
#'   \code{\link{ss_aipe_cv}}
#'
#' @references
#' Kelley, K. (2007). Sample size planning for the coefficient of variation
#'   from the accuracy in parameter estimation approach.
#'   \emph{Behavior Research Methods, 39}(4), 755--766.
#'   \doi{10.3758/BF03192966}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3.)
#'
#' McKay, A. T. (1932). Distribution of the coefficient of variation and
#'   the extended \emph{t} distribution. \emph{Journal of the Royal
#'   Statistical Society, 95}(4), 695--698.
#'
#' Sokal, R. R., & Rohlf, F. J. (1995). \emph{Biometry: The principles and
#'   practice of statistics in biological research} (3rd ed.). W. H. Freeman.
#'
#' Vangel, M. G. (1996). Confidence intervals for a normal coefficient
#'   of variation. \emph{The American Statistician, 50}(1), 21--26.
#'   \doi{10.1080/00031305.1996.10473537}
#'
#' @keywords design htest
#'
#' @export

cv <- function(cv = NULL, mean = NULL, sd = NULL, N = NULL, unbiased = FALSE) {
  if (!is.null(mean) && !is.null(sd)) {
    k <- sd / mean
    if (!is.null(cv)) stop("Because \'mean\' and \'sd\' were specified, do not specify \'cv\'.")
  } else if (!is.null(cv)) {
    k <- cv
  } else {
    stop("You must specify either 'cv' or both 'mean' and 'sd'.", call. = FALSE)
  }

  if (unbiased == TRUE) {
    if (is.null(N)) stop("You must specify 'N' when 'unbiased = TRUE'.", call. = FALSE)
    k <- k * (1 + 1 / (4 * N))
  }

  return(.as_dmar_tbl(data.frame(term = 'cv', value = k)))
}
