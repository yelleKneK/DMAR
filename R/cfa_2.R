#' Two Factor Confirmatory Factor Analysis Model
#'
#' Fits a two factor confirmatory factor analysis model to raw item data
#' or a sample covariance matrix. This is the two factor special case of
#' \code{\link{fa_confirmatory}}: the function is a convenience wrapper that only
#' requires the items of each factor, builds the two factor
#' specification, and forwards everything else to \code{fa_confirmatory()}. The
#' factors are named \code{f1} and \code{f2}, so the rows of the
#' returned table are \code{lambda_f1_1}, \code{lambda_f2_1},
#' \code{phi_f1_f2} (the factor correlation), \code{omega_f1},
#' \code{omega_f2}, and so on, exactly as a two factor \code{fa_confirmatory()}
#' call would report them (the \code{syntax} column names the item
#' behind each number). To
#' name the factors substantively, or for three or more factors, call
#' \code{\link{fa_confirmatory}} directly.
#'
#' Each factor is identified by fixing its variance to 1 and estimating
#' every loading; each item loads on exactly one factor (simple
#' structure). The factor correlation is estimated by default and
#' \code{correlated_factors = FALSE} fixes it to zero. Per-factor
#' constraint vectors use the factor names, for example
#' \code{equal_loading = c(f1 = TRUE, f2 = FALSE)}. A fit whose likelihood
#' has no maximum (an error variance that runs off below zero while the
#' likelihood is flat or still rising) stops with an error of class
#' \code{dmar_no_maximum} that prints the \pkg{lavaan} code of the fit,
#' as \code{\link{fa_confirmatory}} describes; there,
#' \code{bounds = "pos.var"} (which \code{\dots} also passes on from
#' here) or a different model are the choices.
#'
#' @param data A raw data matrix or data frame, rows are respondents and
#'   columns include the items named in \code{factor_1} and
#'   \code{factor_2}. Supply exactly one of \code{data} or \code{S}.
#' @param factor_1 Character vector naming the items of the first
#'   factor (two or more).
#' @param factor_2 Character vector naming the items of the second
#'   factor (two or more). No item may appear in both factors. An item of
#'   either factor that is not a column of \code{data} or an item of
#'   \code{S} is refused with an error of class \code{dmar_items_invalid}
#'   (see \code{\link{dmar_metadata}}).
#' @param S A symmetric covariance matrix of the items, with dimnames
#'   naming the items; \code{N} is then required. Supply exactly one of
#'   \code{data} or \code{S}.
#' @param N Total sample size. Required with \code{S}; ignored (inferred
#'   from the rows) with \code{data}.
#' @inheritParams fa_confirmatory
#' @param \dots Additional arguments forwarded to \code{\link{fa_confirmatory}}
#'   and, through it, to \code{\link[lavaan]{lavaan}} (for example
#'   \code{ordered}, \code{equal_intercept}, or \code{M}).
#'
#' @return The value of the corresponding \code{\link{fa_confirmatory}} call: a
#'   \code{data.frame} (classes \code{dmar_fa_confirmatory}, \code{dmar_tbl}) with
#'   one row per parameter (\code{estimate}, \code{se}, \code{test_statistic},
#'   \code{p_value}, \code{ci_lower}, \code{ci_upper}) followed by the
#'   fit rows, or the alternative shapes selected by \code{output}
#'   (see \code{?fa_confirmatory}).
#'
#' @references
#' Browne, M. W. (1974). Generalized least squares estimators in the
#'   analysis of covariance structures. \emph{South African Statistical
#'   Journal, 8}, 1--24.
#'
#' Browne, M. W. (1984). Asymptotically distribution-free methods for the
#'   analysis of covariance structures. \emph{British Journal of
#'   Mathematical and Statistical Psychology, 37}(1), 62--83.
#'   \doi{10.1111/j.2044-8317.1984.tb00789.x}
#'
#' Muthén, B. (1984). A general structural equation model with
#'   dichotomous, ordered categorical, and continuous latent variable
#'   indicators. \emph{Psychometrika, 49}(1), 115--132.
#'
#' Muthén, B., du Toit, S. H. C., & Spisic, D. (1997). \emph{Robust
#'   inference using weighted least squares and quadratic estimating
#'   equations in latent variable modeling with categorical and continuous
#'   outcomes}. Unpublished technical report.
#'
#' Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics
#'   and standard errors in covariance structure analysis. In A. von Eye &
#'   C. C. Clogg (Eds.), \emph{Latent variables analysis: Applications for
#'   developmental research} (pp. 399--419). Thousand Oaks, CA: Sage.
#'
#' Savalei, V. (2010). Expected versus observed information in SEM with
#'   incomplete normal and nonnormal data. \emph{Psychological Methods,
#'   15}(4), 352--367. \doi{10.1037/a0020143}
#'
#' von Hippel, P. T. (2016). New confidence intervals and bias comparisons
#'   show that maximum likelihood can beat multiple imputation in small
#'   samples. \emph{Structural Equation Modeling, 23}(3), 422--437.
#'   \doi{10.1080/10705511.2015.1047931}
#'
#' Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for
#'   mean and covariance structure analysis with nonnormal missing data.
#'   \emph{Sociological Methodology, 30}, 165--200.
#'
#' @seealso \code{\link{fa_confirmatory}} for the general function this wraps;
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
                  estimator = "ML", missing = "fiml",
                  se = NULL, conf_level = 0.95,
                  output = c("verbose", "measurement", "summary",
                             "standardized", "fit"),
                  ...,
                  wald_df = Inf, rmsea_conf_level = 0.90) {
  if (!is.character(factor_1) || length(factor_1) < 2L ||
      !is.character(factor_2) || length(factor_2) < 2L) {
    .dmar_refuse_argument(
      "'factor_1' and 'factor_2' must each name two or more items.",
      arg = c("factor_1", "factor_2"))
  }
  # An unsupplied 'missing' reaches fa_confirmatory() unsupplied, so that
  # its default resolves there as it does for a direct call. The refusal of
  # a fit whose likelihood has no maximum names this function, the one in
  # the user's call (.no_max_as(), R/no_maximum_internals.R).
  missing_default <- base::missing(missing)
  .no_max_as(if (missing_default) {
    fa_confirmatory(data = data, factors = list(f1 = factor_1, f2 = factor_2),
          S = S, N = N,
          equal_loading = equal_loading, equal_error = equal_error,
          correlated_factors = correlated_factors,
          estimator = estimator, se = se,
          conf_level = conf_level, rmsea_conf_level = rmsea_conf_level,
          output = output, ..., wald_df = wald_df)
  } else {
    fa_confirmatory(data = data, factors = list(f1 = factor_1, f2 = factor_2),
          S = S, N = N,
          equal_loading = equal_loading, equal_error = equal_error,
          correlated_factors = correlated_factors,
          estimator = estimator, missing = missing, se = se,
          conf_level = conf_level, rmsea_conf_level = rmsea_conf_level,
          output = output, ..., wald_df = wald_df)
  }, "cfa_2")
}
