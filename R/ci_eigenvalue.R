# Confidence interval on the largest eigenvalue of a sample covariance.
#' Confidence Interval on the Largest Eigenvalue of a Sample Covariance Matrix
#'
#' Computes an asymptotic confidence interval on the largest population
#' eigenvalue \eqn{\lambda_1} of a population covariance matrix, given a
#' sample covariance matrix from \eqn{n} observations on \eqn{p}
#' variables under multivariate normality. Useful in principal-components
#' analysis and dimension-reduction settings to gauge whether the largest
#' eigenvalue is well-separated from the second.
#'
#' @param cov_matrix Sample covariance matrix (a symmetric, positive-
#'   semidefinite numeric matrix), or a data frame whose columns are
#'   the variables (in which case \code{cov_matrix} is computed
#'   internally).
#' @param n Sample size (number of rows of the original data). Required
#'   when \code{cov_matrix} is supplied as a matrix; ignored when
#'   \code{cov_matrix} is a data frame (\code{nrow(cov_matrix)} is
#'   used instead).
#' @param conf_level Confidence level. Default \code{0.95}.
#' @param k Which eigenvalue (1 = largest, 2 = next, ...) to bracket.
#'   Default \code{1}.
#'
#' @return A 3-row \code{data.frame} with rows ordered
#'   \code{"lower_limit"}, \code{"eigenvalue"} (the sample eigenvalue point
#'   estimate), and \code{"upper_limit"}, so the point estimate sits between
#'   its confidence limits.
#'
#' @details
#' \strong{Asymptotic distribution.} Under multivariate normality with
#' eigenvalues \eqn{\lambda_1 > \lambda_2 \ge \cdots \ge \lambda_p}, the
#' sample eigenvalues \eqn{\hat\lambda_j} are asymptotically independent
#' and approximately normal with mean \eqn{\lambda_j} and variance
#' \eqn{2 \lambda_j^2 / (n - 1)} when the eigenvalues are simple
#' (well-separated) (Anderson, 2003, Theorem 13.3.1; Muirhead, 1982,
#' Section 9.7). The asymptotic CI is therefore
#' \deqn{\hat\lambda_j
#'   \cdot \exp\!\left(\pm z_{1 - \alpha/2} \sqrt{\frac{2}{n - 1}}\right),}
#' on the multiplicative scale (equivalently, a Wald CI on
#' \eqn{\log \lambda_j} with variance \eqn{2/(n - 1)}). The log scale is
#' the natural variance-stabilizing transformation for an eigenvalue.
#'
#' \strong{Caveats.} The asymptotic CI assumes well-separated population
#' eigenvalues. When the largest two eigenvalues are close, the sample
#' eigenvalue exhibits a "repulsion" phenomenon and the CI is biased
#' (typically too narrow). Diagnostic: if \eqn{\hat\lambda_1 /
#' \hat\lambda_2} is close to 1, the asymptotic CI should not be relied
#' upon; a bootstrap is preferable.
#'
#' @references
#' Anderson, T. W. (2003). \emph{An introduction to multivariate
#'   statistical analysis} (4th ed.). Wiley. (See Chapter 13.)
#'
#' Muirhead, R. J. (1982). \emph{Aspects of multivariate statistical
#'   theory}. Wiley. (See Section 9.7.)
#'
#' @seealso \code{\link[stats]{prcomp}}, \code{\link[base]{eigen}}
#'
#' @examples
#' # 1. From a data frame:
#' set.seed(113)
#' X <- data.frame(matrix(rnorm(200), nrow = 50))
#' ci_eigenvalue(X, k = 1)
#'
#' # 2. From an explicit covariance matrix:
#' S <- cov(X)
#' ci_eigenvalue(S, n = nrow(X), k = 1)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family multivariate and latent variable methods
#'
#' @export

ci_eigenvalue <- function(cov_matrix, n = NULL,
                          conf_level = 0.95, k = 1) {
  if (is.data.frame(cov_matrix)) {
    if (is.null(n)) n <- nrow(cov_matrix)
    cov_matrix <- stats::cov(cov_matrix)
  }
  if (!is.matrix(cov_matrix) || !is.numeric(cov_matrix))
    stop("'cov_matrix' must be a numeric matrix or data.frame.")
  if (is.null(n) || !is.numeric(n) || n < 3)
    stop("'n' must be a numeric scalar >= 3.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")
  if (!is.numeric(k) || k < 1 || k > nrow(cov_matrix))
    stop("'k' must be between 1 and ncol(cov_matrix).")

  ev <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
  lambda_k <- ev[k]
  if (lambda_k <= 0) stop("Selected eigenvalue is non-positive; CI undefined.")

  z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)
  log_se  <- sqrt(2 / (n - 1))
  lo <- lambda_k * exp(-z_alpha * log_se)
  hi <- lambda_k * exp( z_alpha * log_se)

  out <- data.frame(
    term  = c("lower_limit", "eigenvalue", "upper_limit"),
    value = c(lo, lambda_k, hi),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
