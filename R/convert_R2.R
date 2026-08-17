#' @rdname convert_R2
#' @name convert_R2
#' @aliases convert_R2_f
#' @aliases convert_f_R2
#' @aliases convert_lambda_R2
#' @aliases convert_R2_lambda
#'
#' @title Convert Between \emph{F}, \eqn{R^2}, and Their Noncentral
#'   Parameters
#'
#' @description
#' Given values of test statistics (and the appropriate additional information) the value of the noncentral
#' values can be obtained. Likewise, given noncentral values (and the appropriate additional information)
#' the value of the test statistic can be obtained.
#'
#' @param R2 Squared multiple correlation coefficient (population or observed)
#' @param df_1 Degrees of freedom for the numerator of the \emph{F}-distribution
#' @param df_2 Degrees of freedom for the denominator of the \emph{F}-distribution
#' @param p Number of predictor variables for \code{R2}
#' @param N Sample size
#' @param F_value The obtained \emph{F} value from a test of significance for the squared multiple correlation coefficient
#' @param lambda The noncentral parameter from an \emph{F}-distribution
#'
#' @details These functions are especially helpful in the search for confidence intervals for noncentral parameters, as they convert to and from related quantities.
#'
#' @return Each of the four functions returns a 1-row \code{data.frame}
#'   with columns \code{term} and \code{value}.  The \code{term} entry
#'   identifies the conversion performed
#'   (\code{"r2_f"}, \code{"f_r2"}, \code{"lambda_r2"}, or
#'   \code{"r2_lambda"}) and \code{value} is the converted scalar. The
#'   conversions are exact inverses of one another (with the appropriate
#'   degrees-of-freedom / sample size inputs supplied), which is what
#'   makes them useful inside the noncentrality-parameter confidence
#'   interval machinery of \code{\link{ci_R2}} and
#'   \code{\link{conf_limits_ncf}}.
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. \emph{Journal of Statistical
#' Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_R2}}, \code{\link{ci_R2}}, \code{\link{conf_limits_nct}}, \code{\link{conf_limits_ncf}}
#'
#' @examples
#' convert_R2_lambda(R2 = .5, N = 100)
#'
#' @keywords multivariate design
#'
#' @family parameterization conversions
#'
#' @export

convert_R2_f <- function(R2 = NULL, df_1 = NULL, df_2 = NULL, p = NULL, N = NULL) {
  if (is.null(df_1) && is.null(df_2) && !is.null(N) && !is.null(p)) {
    df_1 <- p
    df_2 <- N - p - 1
  }
  if (is.null(df_1) || is.null(df_2)) stop("You have not specified \'df_1\', \'df_2\', \'N\', and/or \'p\' correctly.")
  .convert_scalar_check(R2, "R2", lower = 0, upper_open = 1)
  .convert_scalar_check(df_1, "df_1", lower_open = 0)
  .convert_scalar_check(df_2, "df_2", lower_open = 0)
  return(.as_dmar_tbl(data.frame(term = 'r2_f', value = (R2 / df_1) / ((1 - R2) / df_2))))
}

#' @rdname convert_R2
#' @export

convert_f_R2 <- function(F_value = NULL, df_1 = NULL, df_2 = NULL) {
  if (is.null(df_1) || is.null(df_2)) stop("You have not specified \'df_1\' and/or \'df_2\'.")
  .convert_scalar_check(F_value, "F_value", lower = 0)
  .convert_scalar_check(df_1, "df_1", lower_open = 0)
  .convert_scalar_check(df_2, "df_2", lower_open = 0)
  return(.as_dmar_tbl(data.frame(term = 'f_r2', value = F_value * df_1 / (F_value * df_1 + df_2))))
}

#' @rdname convert_R2
#' @export

convert_lambda_R2 <- function(lambda = NULL, N = NULL) {
  if (is.null(lambda) || is.null(N)) stop("You must specify \'lambda\' (i.e., a noncentral parameter) and \'N\' (i.e., sample size) in order to calculate the noncentrality parameter.")
  .convert_scalar_check(lambda, "lambda", lower = 0)
  .convert_scalar_check(N, "N", lower_open = 0)
  return(.as_dmar_tbl(data.frame(term = 'lambda_r2', value = lambda / (lambda + N))))
}

#' @rdname convert_R2
#' @export

convert_R2_lambda <- function(R2 = NULL, N = NULL) {
  if (is.null(R2) || is.null(N)) stop("You must specify \'R2\' (i.e., R Square) and \'N\' (i.e., sample size) in order to calculate the noncentrality parameter.")
  .convert_scalar_check(R2, "R2", lower = 0, upper_open = 1)
  .convert_scalar_check(N, "N", lower_open = 0)
  return(.as_dmar_tbl(data.frame(term = 'r2_lambda', value = R2 / (1 - R2) * N)))
}



# -----------------------------------------------------------------
# Internal fast-path companions to the convert_R2_* family.
# Return bare numerics, skipping the per-call data.frame allocation
# that dominates the runtime of iterative sample size search loops
# (ss_aipe_R2 calls these helpers thousands of times per planning
# call). Used by ci_R2's hot path; public callers continue to use
# the data.frame-returning forms.
# -----------------------------------------------------------------

.convert_R2_f_fast <- function(R2, df_1, df_2) {
  (R2 / df_1) / ((1 - R2) / df_2)
}

.convert_f_R2_fast <- function(F_value, df_1, df_2) {
  F_value * df_1 / (F_value * df_1 + df_2)
}

.convert_lambda_R2_fast <- function(lambda, N) {
  lambda / (lambda + N)
}

.convert_R2_lambda_fast <- function(R2, N) {
  R2 / (1 - R2) * N
}

# Scalar guard shared by the convert_R2_* quartet: every argument is a
# single number, so vector input stops here instead of recycling the
# term column into duplicated rows (for a vector of values, apply the
# arithmetic directly). Bounds are checked where the map's domain ends.
.convert_scalar_check <- function(x, name, lower = -Inf, upper = Inf,
                                  lower_open = NULL, upper_open = NULL) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop("'", name, "' must be a single number. For a vector of values, ",
         "apply the conversion arithmetic directly.", call. = FALSE)
  }
  if (!is.null(lower_open) && x <= lower_open) {
    stop("'", name, "' must be greater than ", lower_open, "; got ", x, ".",
         call. = FALSE)
  }
  if (x < lower) {
    stop("'", name, "' must be at least ", lower, "; got ", x, ".",
         call. = FALSE)
  }
  if (!is.null(upper_open) && x >= upper_open) {
    stop("'", name, "' must be less than ", upper_open, "; got ", x, ".",
         call. = FALSE)
  }
  if (x > upper) {
    stop("'", name, "' must be at most ", upper, "; got ", x, ".",
         call. = FALSE)
  }
  invisible(TRUE)
}
