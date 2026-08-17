#' Single-Common-Factor Screen for Common Method Variance
#'
#' This function implements Harman's single-factor test, the most widely used
#' (and weakest) screen for common method variance: fit a one-factor model to
#' all of the items by maximum likelihood and inspect how much of their
#' variance the common factor accounts for. The rationale is that if a single
#' method factor dominated the responses, one common factor would capture a
#' large share of the variance. A factor accounting for more than half of the
#' variance is the customary red flag (Podsakoff, MacKenzie, Lee, & Podsakoff,
#' 2003). The screen is coarse and cannot by itself rule
#' method variance in or out; the marker-variable and latent method factor
#' approaches are stronger (see \code{\link{common_method_marker}}).
#'
#' Harman's single-factor test (the proportion of variance explained by one
#' common factor) is related to but distinct from the marker-variable
#' technique. A \emph{marker variable} (or common-method marker) is a variable
#' chosen to be theoretically unrelated to the substantive constructs under
#' study, so that any observed correlation between it and the substantive items
#' is attributable to shared method rather than to a true relationship; it is
#' used to estimate or partial out common method variance (Lindell & Whitney,
#' 2001). The single-factor test uses no such marker, it asks only whether a
#' single dimension dominates the item set, so it can flag a strong common
#' factor but cannot identify whether that factor is method or substance.
#'
#' @param data A \code{data.frame} or numeric matrix of item responses.
#'   Supply this, a covariance matrix \code{S}, or a correlation matrix
#'   \code{R} (exactly one).
#' @param S A symmetric covariance matrix among the items, when raw data are
#'   not available but the summary statistics a paper reports are. It is
#'   converted to a correlation matrix internally, so the test acts on the
#'   same scale-free quantity regardless of which input is supplied.
#' @param R A correlation matrix among the items, when raw data are not
#'   available.
#'
#' @details
#' The one-factor model is fit to the item correlation matrix by maximum
#' likelihood with \code{\link[stats]{factanal}}, and the statistic is the
#' proportion of total variance the common factor accounts for: the sum of the
#' squared standardized loadings divided by the number of items (equivalently,
#' the mean communality). Much of the applied literature computes the screen
#' from the largest eigenvalue of the correlation matrix, which describes the
#' first principal component, not a factor; the test is implemented factor
#' analytically here, in the psychometric tradition, because a principal
#' component absorbs unique as well as common variance and so overstates the
#' share a common factor accounts for. Correlations from raw data use
#' pairwise-complete observations. A supplied covariance matrix is first
#' standardized to a correlation matrix with \code{\link[stats]{cov2cor}}.
#' The one-factor model requires at least three items.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows
#'   \code{variance_explained} (the proportion of total variance the single
#'   common factor accounts for) and \code{n_items} in the \code{value}
#'   column.
#'
#' @references
#' Harman, H. H. (1976). \emph{Modern factor analysis} (3rd ed.). University
#'   of Chicago Press.
#'
#' Lindell, M. K., & Whitney, D. J. (2001). Accounting for common method
#'   variance in cross-sectional research designs. \emph{Journal of Applied
#'   Psychology, 86}(1), 114--121. \doi{10.1037/0021-9010.86.1.114}
#'
#' Podsakoff, P. M., MacKenzie, S. B., Lee, J.-Y., & Podsakoff, N. P.
#'   (2003). Common method biases in behavioral research: A critical review
#'   of the literature and recommended remedies. \emph{Journal of Applied
#'   Psychology, 88}(5), 879--903. \doi{10.1037/0021-9010.88.5.879}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{common_method_marker}} for the marker-variable
#'   adjustment.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' set.seed(113)
#' f <- rnorm(200)
#' d <- data.frame(
#'   x1 = f + rnorm(200), x2 = f + rnorm(200), x3 = f + rnorm(200),
#'   x4 = rnorm(200),     x5 = rnorm(200),     x6 = rnorm(200))
#' common_method_single_factor(d)
#'
#' # The same screen from the summary statistics a paper reports.
#' common_method_single_factor(S = cov(d))
#'
#' @export
#' @importFrom stats cor cov2cor
common_method_single_factor <- function(data = NULL, S = NULL, R = NULL) {
  n_supplied <- (!is.null(data)) + (!is.null(S)) + (!is.null(R))
  if (n_supplied != 1L) {
    stop("Supply exactly one of 'data', 'S', or 'R'.", call. = FALSE)
  }
  if (!is.null(data)) {
    if (is.data.frame(data)) data <- as.matrix(data)
    if (!is.numeric(data) || ncol(data) < 3L) {
      stop("'data' must be numeric with three or more item columns.",
           call. = FALSE)
    }
    R <- stats::cor(data, use = "pairwise.complete.obs")
  }
  if (!is.null(S)) {
    if (!is.matrix(S) || !isSymmetric(unname(S)) || nrow(S) < 3L ||
        any(diag(S) <= 0)) {
      stop("'S' must be a symmetric covariance matrix among three or more ",
           "items with positive variances.", call. = FALSE)
    }
    R <- stats::cov2cor(S)
  }
  if (!is.matrix(R) || !isSymmetric(unname(R)) || nrow(R) < 3L ||
      any(abs(R) > 1 + 1e-8)) {
    stop("'R' must be a symmetric correlation matrix among three or more ",
         "items.", call. = FALSE)
  }
  p <- nrow(R)
  fit <- stats::factanal(covmat = R, factors = 1L)
  Lambda <- fit$loadings[, 1L]
  out <- data.frame(
    term  = c("variance_explained", "n_items"),
    value = c(sum(Lambda^2) / p, p),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out)
}
