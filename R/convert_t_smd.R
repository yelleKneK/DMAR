#' @rdname convert_t_smd
#' @name convert_t_smd
#' @aliases convert_lambda_delta
#' @aliases convert_delta_lambda
#'
#' @title Conversion Functions for Noncentral \emph{t}-distribution
#'
#' @description Functions useful for converting a standardized mean difference to a noncentrality parameter, and vice versa.
#'
#' @param lambda noncentral value from a \emph{t}-distribution
#' @param delta Population value of the standardized mean difference
#' @param n_1 Sample size in group 1
#' @param n_2 Sample size in group 2
#'
#' @details
#' Although \code{lambda} is the population noncentral value, an estimate of it is the observed value of a
#' \emph{t}-statistic. Likewise, delta can be estimated as the observed standardized mean difference. Thus, the observed
#' standardized mean difference can be converted to the observed \emph{t}-value. These functions are especially helpful in the
#' context of forming confidence intervals for the population standardized mean difference.
#'
#' @return Each function returns a 1-row \code{data.frame} with columns
#'   \code{term} and \code{value}. The \code{term} entry identifies the
#'   conversion (\code{"delta_lambda"} or \code{"lambda_delta"}) and
#'   \code{value} is the converted scalar. The two functions are exact
#'   inverses given the \emph{per-group} sample sizes
#'   \code{n_1} and \code{n_2}.
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. \emph{Journal of Statistical
#' Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{smd}}, \code{\link{ci_smd}}, \code{\link{ss_aipe_smd}}
#'
#' @examples
#' convert_lambda_delta(lambda = 2, n_1 = 113, n_2 = 113)
#' convert_delta_lambda(delta = .266076, n_1 = 113, n_2 = 113)
#'
#' @keywords design
#'
#' @family parameterization conversions
#'
#' @export

convert_delta_lambda <- function(delta, n_1, n_2) {
  if (!is.numeric(delta) || length(delta) != 1L || is.na(delta)) {
    stop("'delta' must be a single number.", call. = FALSE)
  }
  .convert_scalar_check(n_1, "n_1", lower_open = 0)
  .convert_scalar_check(n_2, "n_2", lower_open = 0)
  term <- "delta_lambda"
  value <- delta * sqrt((n_2 * n_1) / (n_1 + n_2))
  return(.as_dmar_tbl(data.frame(term, value)))
}

#' @rdname convert_t_smd
#' @export

convert_lambda_delta <- function(lambda, n_1, n_2) {
  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda)) {
    stop("'lambda' must be a single number.", call. = FALSE)
  }
  .convert_scalar_check(n_1, "n_1", lower_open = 0)
  .convert_scalar_check(n_2, "n_2", lower_open = 0)
  term <- "lambda_delta"
  value <- lambda * sqrt((n_2 + n_1) / (n_1 * n_2))
  return(.as_dmar_tbl(data.frame(term, value)))
}
