# Critical value(s) for a chi square distribution.
#' Provides the Critical Value(s) for a Chi Square Distribution
#'
#' @param alpha_level Type I error rate (i.e., the false positive rate).
#' @param df The number of degrees of freedom (a positive number).
#' @param alternative The type of alternative hypothesis of interest. The
#'   default, \code{"greater"}, puts the whole of \code{alpha_level} in the upper
#'   tail, which is how the chi square distribution is used to test a model
#'   or an association (see Details).
#' @param alpha_lower The error rate in the lower tail of the distribution.
#' @param alpha_upper The error rate in the upper tail of the distribution.
#' @param ncp The noncentral parameter (if zero, the default, it is the
#'   central chi square distribution).
#' @param verbose Provides extra information about areas under the curve.
#'
#' @return Returns the critical value(s), based on the input specifications,
#'   in a output style (a \code{data.frame} with a row for the lower and the
#'   upper critical value, following the format used by \code{\link{cv_t}}).
#'
#' @details Like the \emph{F} distribution and unlike \emph{t} and \emph{z},
#'   the chi square distribution is not symmetric and takes only non-negative
#'   values. Its common uses are one-sided in the upper tail: a test of
#'   association in a contingency table, a likelihood ratio test, and a test
#'   of model fit all reject for large values, because a poorly fitting model
#'   produces a large discrepancy, never a small one. That is why
#'   \code{alternative} defaults to \code{"greater"} here whereas it defaults
#'   to \code{"not_equal"} in \code{\link{cv_t}}. Maxwell, Delaney, and Kelley
#'   (2027) tabulate these upper-tail values in their Appendix Table A.9.
#'
#'   Both tails remain available for the situations that need them, most
#'   commonly an interval for a variance, which uses an upper and a lower chi
#'   square quantile. Set \code{alternative = "not_equal"}, or give
#'   \code{alpha_lower} and \code{alpha_upper} directly. When a tail is given
#'   zero area its critical value is the boundary of the support, so
#'   \code{lower_cv} is 0 under the default.
#'
#'   A noncentral parameter can be supplied, which is what a power analysis
#'   for a test of model fit needs, though it would not be used for a standard
#'   null hypothesis significance test. See \code{\link{conf_limits_nc_chisq}}
#'   for confidence limits on the noncentral parameter itself.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (Appendix Table A.9 reports these critical
#'   values.)
#'
#' @examples
#' # The critical value for a test on 3 degrees of freedom at the .05 level.
#' cv_chisq(alpha_level = .05, df = 3)
#'
#' # Simple output.
#' cv_chisq(alpha_level = .05, df = 3, verbose = FALSE)
#'
#' # Both tails, as an interval for a variance would need.
#' cv_chisq(alpha_level = .05, df = 10, alternative = "not_equal")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cv_t}}, \code{\link{cv_f}},
#'   \code{\link{conf_limits_nc_chisq}}
#'
#' @keywords distribution htest
#'
#' @family critical values
#'
#' @export
#' @import stats

cv_chisq <- function(alpha_level, df, alternative = "greater",
                     alpha_lower, alpha_upper, ncp = 0, verbose = TRUE) {
  if (missing(df)) stop("You must specify the degrees of freedom (i.e., 'df'), which must be a positive number.")
  if (!df > 0) stop("You must specify a positive value for the degrees of freedom.")

  # A negative 'ncp' would route qchisq()/pchisq() into the noncentral algorithm
  # with an out-of-domain argument and return a polished NaN rather than a
  # critical value; require a finite, non-negative noncentral parameter (zero
  # being the central chi square).
  if (!is.numeric(ncp) || length(ncp) != 1L || !is.finite(ncp) || ncp < 0)
    stop("'ncp' must be a single finite, non-negative number (zero gives the central chi square distribution).")

  if (missing(alpha_level) && missing(alpha_lower) && missing(alpha_upper)) stop("You must specify either 'alpha_level' or 'alpha_lower' and 'alpha_upper'.")

  if (missing(alpha_level)) alpha_level <- NULL
  if (missing(alpha_lower)) alpha_lower <- NULL
  if (missing(alpha_upper)) alpha_upper <- NULL

  if (!is.null(alpha_level)) {
    if (!is.null(alpha_lower)) stop("You have specified both 'alpha_level' and 'alpha_lower'; only use one approach.")
    if (!is.null(alpha_upper)) stop("You have specified both 'alpha_level' and 'alpha_upper'; only use one approach.")
    if (alpha_level <= 0 || alpha_level >= 1) stop("Specify 'alpha_level' to be greater than zero and less than 1.")

    if (tolower(alternative) %in% c("greater than", "greater-than", "greater", "greater.than", "gt", "g", ">", ">=")) {
      alpha_lower <- 0
      alpha_upper <- alpha_level
    }

    if (tolower(alternative) %in% c("less than", "less-than", "less", "lesser", "less.than", "lt", "l", "<", "<=")) {
      alpha_lower <- alpha_level
      alpha_upper <- 0
    }

    if (tolower(alternative) %in% c("ne", "2s", "not equal", "two.sided", "two sided", "two-sided", "!=", "not_equal")) {
      alpha_lower <- alpha_upper <- alpha_level / 2
    }
  }

  if (is.null(alpha_level)) {
    if (is.null(alpha_lower)) stop("With 'alpha_level=NULL' you need to specify 'alpha_lower' (and alpha_upper)")
    if (is.null(alpha_upper)) stop("With 'alpha_level=NULL' you need to specify 'alpha_upper' (and alpha_lower)")
  }

  if (is.null(alpha_lower) || is.null(alpha_upper)) stop("'alternative' must be one of \"not_equal\", \"greater\", or \"less\".")
  if (alpha_lower < 0 || alpha_lower >= .5) stop("Specify 'alpha_lower' to be greater than or equal to zero but less than .50")
  if (alpha_upper < 0 || alpha_upper >= .5) stop("Specify 'alpha_upper' to be greater than or equal to zero but less than .50.")

  # Outputs
  term <- c("lower_cv", "upper_cv")
  # As in cv_f(), pass 'ncp' only when it is nonzero, so the central algorithm
  # serves the central case rather than the noncentral one standing in for it.
  if (ncp == 0) {
    value <- qchisq(p = c(alpha_lower, 1 - alpha_upper), df = df)
    area_less <- pchisq(value, df = df, lower.tail = TRUE)
    area_greater <- pchisq(value, df = df, lower.tail = FALSE)
  } else {
    value <- qchisq(p = c(alpha_lower, 1 - alpha_upper), df = df, ncp = ncp)
    area_less <- pchisq(value, df = df, ncp = ncp, lower.tail = TRUE)
    area_greater <- pchisq(value, df = df, ncp = ncp, lower.tail = FALSE)
  }

  if (verbose == TRUE) {
    return(.as_dmar_tbl(data.frame(term, value, area_less, area_greater)))
  }
  .as_dmar_tbl(data.frame(term, value))
}
