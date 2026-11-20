#' One Factor Confirmatory Factor Analysis Model
#'
#' Fits a single factor (congeneric by default) confirmatory factor
#' analysis model to raw item data or a sample covariance matrix. This
#' is the one factor special case of \code{\link{fa_confirmatory}}: the function
#' is a convenience wrapper that only requires the data (and, when the
#' data hold more than the items, a vector of item names), builds the
#' one factor specification, and forwards everything else to
#' \code{fa_confirmatory()}. The factor is named \code{f1}, so the rows of the
#' returned table are \code{lambda_f1_1}, \code{lambda_f1_2}, ...,
#' \code{psi_f1_1}, ..., and \code{omega_f1}, exactly as a one factor
#' \code{fa_confirmatory()} call would report them (the \code{syntax} column
#' names the item behind each number).
#'
#' The model is identified by fixing the factor variance to 1 and
#' estimating every loading. \code{equal_loading} and \code{equal_error}
#' impose the classical measurement structures on the single factor, and
#' the header of the printed table names the structure implied by the
#' constraints. For the composite reliability coefficient with the
#' observed total variance in the denominator, use
#' \code{\link{reliability_omega}} with \code{denominator = "observed"};
#' for the model implied omega, the \code{omega_f1} row of this
#' function's output and \code{\link{reliability_omega}} agree.
#'
#' Ordered categorical items are not supported here; use
#' \code{\link{fa_confirmatory}}, whose \code{ordered} argument fits WLSMV with
#' the theta parameterization and reports categorical omega on the sum
#' score metric, or
#' \code{\link{reliability_omega_categorical}}.
#'
#' A fit whose likelihood has no maximum (an error variance that runs off
#' below zero while the likelihood is flat or still rising, as it can on
#' three items whose correlations have a negative product) stops with an
#' error of class \code{dmar_no_maximum} that prints the \pkg{lavaan}
#' code of the fit, as \code{\link{fa_confirmatory}} describes; there,
#' \code{bounds = "pos.var"} (which \code{\dots} also passes on from
#' here) or a different model are the choices.
#'
#' @param data A raw data matrix or data frame, rows are respondents and
#'   columns include the items. A matrix without column names is given
#'   the names \code{y1}, \code{y2}, ... Supply exactly one of
#'   \code{data} or \code{S}.
#' @param items Character vector naming the items of the factor (three
#'   or more; two are accepted with \code{equal_loading = TRUE}, the
#'   just identified tau-equivalent case). The default \code{NULL} uses every column of \code{data}
#'   (or of \code{S}), so a data set that holds only the items needs no
#'   \code{items} at all. An item that is not a column of \code{data} or
#'   an item of \code{S} is refused with an error of class
#'   \code{dmar_items_invalid} (see \code{\link{dmar_metadata}}).
#' @param S A symmetric covariance matrix of the items; \code{N} is then
#'   required. Dimnames are optional here: a matrix without them is
#'   given the item names \code{y1}, \code{y2}, ... Supply exactly one
#'   of \code{data} or \code{S}.
#' @param N Total sample size. Required with \code{S}; ignored (inferred
#'   from the rows) with \code{data}.
#' @inheritParams fa_confirmatory
#' @param \dots Additional arguments forwarded to \code{\link{fa_confirmatory}}
#'   and, through it, to \code{\link[lavaan]{lavaan}}.
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
#' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores
#'   using structural equation modeling: An alternative to coefficient
#'   alpha. \emph{Psychometrika, 74}(1), 155--167.
#'   \doi{10.1007/s11336-008-9099-3}
#'
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
#'   population reliability coefficients: Evaluation of methods,
#'   recommendations, and software for composite measures.
#'   \emph{Psychological Methods, 21}(1), 69--92. \doi{10.1037/a0040086}
#'
#' McDonald, R. P. (1999). \emph{Test theory: A unified treatment}.
#'   Lawrence Erlbaum Associates.
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
#' @seealso \code{\link{fa_confirmatory}} for the general function this wraps
#'   (factor analysis with any number of factors, intercept constraints,
#'   ordered categorical items, and the measurement output);
#'   \code{\link{cfa_2}} for the two factor wrapper;
#'   \code{\link{reliability_omega}} for coefficient omega with
#'   confidence intervals.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' set.seed(113)
#' f <- rnorm(200)
#' loadings <- c(0.5, 0.6, 0.65, 0.7, 0.8)
#' X <- sapply(loadings, function(l) l * f + rnorm(200, sd = sqrt(1 - l^2)))
#' colnames(X) <- paste0("y", 1:5)
#'
#' # All columns are items, so the data are all that is needed.
#' cfa_1(X)
#'
#' # Equal loadings (essentially tau-equivalent), named in the header.
#' cfa_1(X, equal_loading = TRUE)
#'
#' # From a covariance matrix and sample size, as a paper reports them,
#' # here with a subset of its variables named in \'items\' (which selects
#' # columns of raw data in the same way).
#' cfa_1(S = cov(X), N = 200, items = c("y1", "y2", "y3", "y4"))
#'
#' @export
cfa_1 <- function(data = NULL, items = NULL, S = NULL, N = NULL,
                  equal_loading = FALSE, equal_error = FALSE,
                  estimator = "ML", missing = "fiml",
                  se = NULL, conf_level = 0.95,
                  output = c("verbose", "measurement", "summary",
                             "standardized", "fit"),
                  ...,
                  wald_df = Inf, rmsea_conf_level = 0.90) {
  if (is.null(data) == is.null(S)) {
    .dmar_refuse_argument(c(
      "Supply exactly one of 'data' (raw data) or 'S' (a symmetric ",
      "covariance matrix, with 'N')."), arg = c("data", "S", "N"))
  }
  if ("ordered" %in% names(list(...))) {
    .dmar_refuse_argument(c(
      "cfa_1() does not support ordered categorical items. Use ",
      "fa_confirmatory(), whose 'ordered' argument fits WLSMV with the theta ",
      "parameterization and reports categorical omega on the sum ",
      "score metric, or ",
      "reliability_omega_categorical()."), arg = NULL)
  }

  # The wrapper's convenience over fa_confirmatory(): unnamed input is auto-named,
  # and 'items' defaults to every column, so cfa_1(X) is a complete call.
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data)) {
      .dmar_refuse_argument(c(
        "'data' must be a data frame or matrix of raw data (rows are ",
        "respondents, columns include the items)."), arg = "data")
    }
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("y", seq_len(ncol(data)))
    }
    all_names <- colnames(data)
  } else {
    if (!is.matrix(S) || !isSymmetric(unname(S), tol = 1e-05)) {
      .dmar_refuse_argument(c(
        "'S' must be a symmetric covariance matrix; raw data are ",
        "passed through 'data'."), arg = c("S", "data"))
    }
    if (is.null(colnames(S))) {
      colnames(S) <- rownames(S) <- paste0("y", seq_len(ncol(S)))
    }
    all_names <- colnames(S)
  }
  if (is.null(items)) items <- all_names
  if (!is.character(items) || length(items) < 2L) {
    .dmar_refuse_argument(
      "'items' must name the factor's columns.", arg = "items")
  }

  # An unsupplied 'missing' reaches fa_confirmatory() unsupplied, so that
  # its default resolves there as it does for a direct call. The refusal of
  # a fit whose likelihood has no maximum names this function, the one in
  # the user's call (.no_max_as(), R/no_maximum_internals.R).
  missing_default <- base::missing(missing)
  .no_max_as(if (missing_default) {
    fa_confirmatory(data = data, factors = list(f1 = items), S = S, N = N,
          equal_loading = equal_loading, equal_error = equal_error,
          estimator = estimator, se = se,
          conf_level = conf_level, rmsea_conf_level = rmsea_conf_level,
          output = output, ..., wald_df = wald_df)
  } else {
    fa_confirmatory(data = data, factors = list(f1 = items), S = S, N = N,
          equal_loading = equal_loading, equal_error = equal_error,
          estimator = estimator, missing = missing, se = se,
          conf_level = conf_level, rmsea_conf_level = rmsea_conf_level,
          output = output, ..., wald_df = wald_df)
  }, "cfa_1")
}
