#' Two Factor Confirmatory Factor Analysis Model
#'
#' Fits a two factor confirmatory factor analysis model to raw item data
#' or a sample covariance matrix. This is the two factor special case of
#' \code{\link{cfa_k}}: the function is a convenience wrapper that only
#' requires the items of each factor, builds the two factor
#' specification, and forwards everything else to \code{cfa_k()}. The
#' factors are named \code{f1} and \code{f2}, so the rows of the
#' returned table are \code{lambda_f1_1}, \code{lambda_f2_1},
#' \code{phi_f1_f2} (the factor correlation), \code{omega_f1},
#' \code{omega_f2}, and so on, exactly as a two factor \code{cfa_k()}
#' call would report them (the \code{syntax} column names the item
#' behind each number). To
#' name the factors substantively, or for three or more factors, call
#' \code{\link{cfa_k}} directly.
#'
#' Each factor is identified by fixing its variance to 1 and estimating
#' every loading; each item loads on exactly one factor (simple
#' structure). The factor correlation is estimated by default and
#' \code{correlated_factors = FALSE} fixes it to zero. Per-factor
#' constraint vectors use the factor names, for example
#' \code{equal_loading = c(f1 = TRUE, f2 = FALSE)}.
#'
#' @param data A raw data matrix or data frame, rows are respondents and
#'   columns include the items named in \code{factor_1} and
#'   \code{factor_2}. Supply exactly one of \code{data} or \code{S}.
#' @param factor_1 Character vector naming the items of the first
#'   factor (two or more).
#' @param factor_2 Character vector naming the items of the second
#'   factor (two or more). No item may appear in both factors.
#' @param S A symmetric covariance matrix of the items, with dimnames
#'   naming the items; \code{N} is then required. Supply exactly one of
#'   \code{data} or \code{S}.
#' @param N Total sample size. Required with \code{S}; ignored (inferred
#'   from the rows) with \code{data}.
#' @inheritParams cfa_k
#' @param \dots Additional arguments forwarded to \code{\link{cfa_k}}
#'   and, through it, to \code{\link[lavaan]{lavaan}} (for example
#'   \code{ordered}, \code{equal_intercept}, or \code{M}).
#'
#' @return The value of the corresponding \code{\link{cfa_k}} call: a
#'   \code{data.frame} (classes \code{dmar_cfa_k}, \code{dmar_tbl}) with
#'   one row per parameter (\code{estimate}, \code{se}, \code{z_value},
#'   \code{p_value}, \code{ci_lower}, \code{ci_upper}) followed by the
#'   fit rows, or the alternative shapes selected by \code{output}
#'   (see \code{?cfa_k}).
#'
#' @seealso \code{\link{cfa_k}} for the general function this wraps;
#'   \code{\link{cfa_1}} for the one factor wrapper;
#'   \code{\link{htmt}} and \code{\link{average_variance_extracted}}
#'   for the discriminant and convergent validity summaries the
#'   \code{output = "measurement"} table reports alongside omega.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' data(holzinger_swineford)
#'
#' # Two factors, each named by its items.
#' cfa_2(holzinger_swineford,
#'       factor_1 = c("t6_paragraph_comprehension", "t7_sentence",
#'                    "t9_word_meaning"),
#'       factor_2 = c("t20_deduction", "t22_problem_reasoning",
#'                    "t23_series_completion"))
#'
#' # The measurement properties: omega, ave, and H per factor, the
#' # factor correlation, and the htmt ratio.
#' cfa_2(holzinger_swineford,
#'       factor_1 = c("t6_paragraph_comprehension", "t7_sentence",
#'                    "t9_word_meaning"),
#'       factor_2 = c("t20_deduction", "t22_problem_reasoning",
#'                    "t23_series_completion"),
#'       output = "measurement")
#'
#' @export
cfa_2 <- function(data = NULL, factor_1, factor_2, S = NULL, N = NULL,
                  equal_loading = FALSE, equal_error = FALSE,
                  correlated_factors = TRUE,
                  estimator = "ML", missing = "listwise",
                  se = "standard", conf_level = 0.95,
                  output = c("verbose", "measurement", "summary",
                             "standardized", "fit"),
                  ...) {
  if (!is.character(factor_1) || length(factor_1) < 2L ||
      !is.character(factor_2) || length(factor_2) < 2L) {
    stop("'factor_1' and 'factor_2' must each name two or more items.",
         call. = FALSE)
  }
  cfa_k(data = data, factors = list(f1 = factor_1, f2 = factor_2),
        S = S, N = N,
        equal_loading = equal_loading, equal_error = equal_error,
        correlated_factors = correlated_factors,
        estimator = estimator, missing = missing, se = se,
        conf_level = conf_level, output = output, ...)
}
