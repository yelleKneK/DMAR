# Critical value(s) for an F distribution.
#' Provides the Critical Value(s) for an \emph{F} Distribution
#'
#' @param alpha_level Type I error rate (i.e., the false positive rate).
#' @param df_numerator The numerator degrees of freedom (a positive number).
#'   In a model comparison this is the difference in the number of parameters
#'   between the two models.
#' @param df_denominator The denominator (error) degrees of freedom (a
#'   positive number).
#' @param alternative The type of alternative hypothesis of interest. The
#'   default, \code{"greater"}, puts the whole of \code{alpha_level} in the upper
#'   tail, which is how the \emph{F} distribution is used to test a model
#'   comparison (see Details).
#' @param alpha_lower The error rate in the lower tail of the distribution.
#' @param alpha_upper The error rate in the upper tail of the distribution.
#' @param ncp The noncentral parameter (if zero, the default, it is the
#'   central \emph{F} distribution).
#' @param verbose Provides extra information about areas under the curve.
#'
#' @return Returns the critical value(s), based on the input specifications,
#'   in a output style (a \code{data.frame} with a row for the lower and the
#'   upper critical value, following the format used by \code{\link{cv_t}}).
#'
#' @details Unlike the \emph{t} and \emph{z} distributions, the \emph{F}
#'   distribution is not symmetric and takes only non-negative values, and the
#'   usual test of a model comparison is one-sided: a restricted model fits
#'   worse than a full model, so evidence against the restriction shows up as
#'   a large \emph{F}, never a small one. That is why \code{alternative}
#'   defaults to \code{"greater"} here whereas it defaults to
#'   \code{"not_equal"} in \code{\link{cv_t}}. Maxwell, Delaney, and Kelley
#'   (2027) tabulate these upper-tail values in their Appendix Table A.2.
#'
#'   Both tails remain available for the situations that need them, such as an
#'   interval for a ratio of variances, either by setting
#'   \code{alternative = "not_equal"} or by giving \code{alpha_lower} and
#'   \code{alpha_upper} directly. When a tail is given zero area its critical
#'   value is the boundary of the support, so \code{lower_cv} is 0 under the
#'   default.
#'
#'   A noncentral parameter can be supplied, which is what a power analysis
#'   needs, though it would not be used for a standard null hypothesis
#'   significance test. See \code{\link{conf_limits_ncf}} for confidence
#'   limits on the noncentral parameter itself.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3, where the \emph{F} test of a model
#'   comparison is developed; Appendix Table A.2 reports these critical
#'   values.)
#'
#' @examples
#' # The critical value for a model comparison with 3 numerator and 20
#' # denominator degrees of freedom, at the .05 level.
#' cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20)
#'
#' # An omnibus test of four groups with 24 participants: a - 1 = 3 and
#' # N - a = 20 degrees of freedom; simple output.
#' cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20, verbose = FALSE)
#'
#' # Both tails, as an interval for a ratio of variances would need.
#' cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20,
#'      alternative = "not_equal")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{cv_t}}, \code{\link{cv_chisq}},
#'   \code{\link{cv_bonferroni_f}}, \code{\link{cv_scheffe}},
#'   \code{\link{conf_limits_ncf}}
#'
#' @keywords distribution htest
#'
#' @family critical values
#'
#' @export
#' @import stats

cv_f <- function(alpha_level, df_numerator, df_denominator, alternative = "greater",
                 alpha_lower, alpha_upper, ncp = 0, verbose = TRUE) {
  if (missing(df_numerator)) stop("You must specify the numerator degrees of freedom (i.e., 'df_numerator'), which must be a positive number.")
  if (missing(df_denominator)) stop("You must specify the denominator degrees of freedom (i.e., 'df_denominator'), which must be a positive number.")
  if (!df_numerator > 0) stop("You must specify a positive value for 'df_numerator'.")
  if (!df_denominator > 0) stop("You must specify a positive value for 'df_denominator'.")

  # A negative 'ncp' would route qf()/pf() into the noncentral algorithm with an
  # out-of-domain argument and return a polished NaN rather than a critical
  # value; require a finite, non-negative noncentral parameter (zero being the
  # central F).
  if (!is.numeric(ncp) || length(ncp) != 1L || !is.finite(ncp) || ncp < 0)
    stop("'ncp' must be a single finite, non-negative number (zero gives the central F distribution).")

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
  # Pass 'ncp' only when it is nonzero. Supplying ncp = 0 explicitly sends qf()
  # and pf() down their noncentral algorithm, which does not admit an infinite
  # numerator df and returns NaN there, whereas the central algorithm handles it
  # (the infinite-numerator column of Appendix Table A.2 is the case in point).
  if (ncp == 0) {
    value <- qf(p = c(alpha_lower, 1 - alpha_upper), df1 = df_numerator, df2 = df_denominator)
    area_less <- pf(value, df1 = df_numerator, df2 = df_denominator, lower.tail = TRUE)
    area_greater <- pf(value, df1 = df_numerator, df2 = df_denominator, lower.tail = FALSE)
  } else {
    value <- qf(p = c(alpha_lower, 1 - alpha_upper), df1 = df_numerator, df2 = df_denominator, ncp = ncp)
    area_less <- pf(value, df1 = df_numerator, df2 = df_denominator, ncp = ncp, lower.tail = TRUE)
    area_greater <- pf(value, df1 = df_numerator, df2 = df_denominator, ncp = ncp, lower.tail = FALSE)
  }

  if (verbose == TRUE) {
    return(.as_dmar_tbl(data.frame(term, value, area_less, area_greater)))
  }
  .as_dmar_tbl(data.frame(term, value))
}
