#' One Factor Confirmatory Factor Analysis Model
#'
#' Fits a single factor (congeneric by default) confirmatory factor
#' analysis model to raw item data or a sample covariance matrix. This
#' is the one factor special case of \code{\link{cfa_k}}: the function
#' is a convenience wrapper that only requires the data (and, when the
#' data hold more than the items, a vector of item names), builds the
#' one factor specification, and forwards everything else to
#' \code{cfa_k()}. The factor is named \code{f1}, so the rows of the
#' returned table are \code{lambda_f1_1}, \code{lambda_f1_2}, ...,
#' \code{psi_f1_1}, ..., and \code{omega_f1}, exactly as a one factor
#' \code{cfa_k()} call would report them (the \code{syntax} column
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
#' \code{\link{cfa_k}}, whose \code{ordered} argument fits WLSMV with
#' the theta parameterization and reports the Green and Yang (2009)
#' categorical sum score omega, or
#' \code{\link{reliability_omega_categorical}}.
#'
#' @param data A raw data matrix or data frame, rows are respondents and
#'   columns include the items. A matrix without column names is given
#'   the names \code{y1}, \code{y2}, ... Supply exactly one of
#'   \code{data} or \code{S}.
#' @param items Character vector naming the items of the factor (three
#'   or more; two are accepted with \code{equal_loading = TRUE}, the
#'   just identified tau-equivalent case). The default \code{NULL} uses every column of \code{data}
#'   (or of \code{S}), so a data set that holds only the items needs no
#'   \code{items} at all.
#' @param S A symmetric covariance matrix of the items; \code{N} is then
#'   required. Dimnames are optional here: a matrix without them is
#'   given the item names \code{y1}, \code{y2}, ... Supply exactly one
#'   of \code{data} or \code{S}.
#' @param N Total sample size. Required with \code{S}; ignored (inferred
#'   from the rows) with \code{data}.
#' @inheritParams cfa_k
#' @param \dots Additional arguments forwarded to \code{\link{cfa_k}}
#'   and, through it, to \code{\link[lavaan]{lavaan}}.
#'
#' @return The value of the corresponding \code{\link{cfa_k}} call: a
#'   \code{data.frame} (classes \code{dmar_cfa_k}, \code{dmar_tbl}) with
#'   one row per parameter (\code{estimate}, \code{se}, \code{z_value},
#'   \code{p_value}, \code{ci_lower}, \code{ci_upper}) followed by the
#'   fit rows, or the alternative shapes selected by \code{output}
#'   (see \code{?cfa_k}).
#'
#' @references
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
#' @seealso \code{\link{cfa_k}} for the general function this wraps
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
#' # From a covariance matrix and sample size, as a paper reports them.
#' cfa_1(S = cov(X), N = 200)
#'
#' # A subset of columns via \'items\'.
#' cfa_1(X, items = c("y1", "y2", "y3", "y4"))
#'
#' @export
cfa_1 <- function(data = NULL, items = NULL, S = NULL, N = NULL,
                  equal_loading = FALSE, equal_error = FALSE,
                  estimator = "ML", missing = "listwise",
                  se = "standard", conf_level = 0.95,
                  output = c("verbose", "measurement", "summary",
                             "standardized", "fit"),
                  ...) {
  if (is.null(data) == is.null(S)) {
    stop("Supply exactly one of 'data' (raw data) or 'S' (a symmetric ",
         "covariance matrix, with 'N').", call. = FALSE)
  }
  if ("ordered" %in% names(list(...))) {
    stop("cfa_1() does not support ordered categorical items. Use ",
         "cfa_k(), whose 'ordered' argument fits WLSMV with the theta ",
         "parameterization and reports the Green and Yang (2009) ",
         "categorical sum score omega, or ",
         "reliability_omega_categorical().", call. = FALSE)
  }

  # The wrapper's convenience over cfa_k(): unnamed input is auto-named,
  # and 'items' defaults to every column, so cfa_1(X) is a complete call.
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data)) {
      stop("'data' must be a data frame or matrix of raw data (rows are ",
           "respondents, columns include the items).", call. = FALSE)
    }
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("y", seq_len(ncol(data)))
    }
    all_names <- colnames(data)
  } else {
    if (!is.matrix(S) || !isSymmetric(unname(S), tol = 1e-05)) {
      stop("'S' must be a symmetric covariance matrix; raw data are ",
           "passed through 'data'.", call. = FALSE)
    }
    if (is.null(colnames(S))) {
      colnames(S) <- rownames(S) <- paste0("y", seq_len(ncol(S)))
    }
    all_names <- colnames(S)
  }
  if (is.null(items)) items <- all_names
  if (!is.character(items) || length(items) < 2L) {
    stop("'items' must name the factor's columns.", call. = FALSE)
  }

  cfa_k(data = data, factors = list(f1 = items), S = S, N = N,
        equal_loading = equal_loading, equal_error = equal_error,
        estimator = estimator, missing = missing, se = se,
        conf_level = conf_level, output = output, ...)
}
