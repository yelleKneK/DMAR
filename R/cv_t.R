# Critical value(s) for a \emph{t}-distribution
#' Provides the Critical Value(s) for a \emph{t}-distribution
#'
#' @param alpha_level Type I error rate (i.e., the false positive rate).
#' @param df The number of degrees of freedom (a positive number)
#' @param alternative The type of alternative hypothesis of interest.
#' @param alpha_lower The error rate on the lower (negative) side of the distribution.
#' @param alpha_upper The error rate on the upper (positive) side of the distribution.
#' @param ncp The noncentral parameter (if zero, the default, it is the central \emph{t}-distribution).
#' @param verbose Provides extra information about areas under the curve.
#'
#' @return Returns the critical value(s), based on the input specifications, in a output style.
#'
#' @details Though a noncentral parameter can be included, that would not be done for a standard null hypothesis significance test.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @examples
#' # A basic call for finding critical values with equal area in the two tails.
#' cv_t(alpha_level = .05, df = 13)
#'
#' # A basic call for a single-sided confidence interval (for "a greater than" alternative hypothesis)
#' cv_t(alpha_level = .05, df = 13, alternative = "greater")
#'
#' # A single-sided confidence interval (for "a greater than" alternative hypothesis); simple output.
#' cv_t(alpha_lower = 0, alpha_upper = .05, df = 13, verbose = FALSE)
#'
#' # For a nonsymmetric 95% confidence interval.
#' cv_t(alpha_lower = .01, alpha_upper = .04, df = 13)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @family critical values
#'
#' @export
#' @import stats

cv_t <- function(alpha_level, df, alternative = "not_equal", alpha_lower, alpha_upper, ncp = 0, verbose = TRUE) {
  if (missing(df)) stop("You must specify the degrees of freedom (i.e., 'df'), which must be a positive number.")
  if (!df > 0) stop("You must specify a positive value for the degrees of freedom.")

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

  if (alpha_lower < 0 || alpha_lower >= .5) stop("Specify 'alpha_lower' to be greater than or equal to zero but less than .50")
  if (alpha_upper < 0 || alpha_upper >= .5) stop("Specify 'alpha_upper' to be greater than or equal to zero but less than .50.")

  # Outputs
  term <- c("lower_cv", "upper_cv")
  value <- qt(p = c(alpha_lower, 1 - alpha_upper), df = df, ncp = ncp)
  area_less <- c(pt(value[1], df = df, ncp = ncp, lower.tail = TRUE), pt(value[2], df = df, ncp = ncp, lower.tail = TRUE))
  area_greater <- c(pt(value[1], df = df, ncp = ncp, lower.tail = FALSE), pt(value[2], df = df, ncp = ncp, lower.tail = FALSE))

  if (verbose == TRUE) {
    return(.as_dmar_tbl(data.frame(term, value, area_less, area_greater)))
  }
  if (verbose == FALSE) {
    return(.as_dmar_tbl(data.frame(term, value)))
  }
}

