#' Marker-Variable Adjustment for Common Method Variance
#'
#' The marker-variable technique of Lindell and Whitney (2001) estimates
#' common method variance from the correlation of a \emph{marker} variable
#' that is theoretically unrelated to at least one of the substantive
#' variables: any non-zero correlation it shows with that variable is
#' attributed to shared method, and that amount is partialled out of the
#' substantive correlations. When no a priori marker is available, the
#' smallest positive correlation among the substantive items is used as a
#' proxy, the common marker-free variant of the method. A correlation that
#' remains statistically significant after the adjustment, by the paper's
#' \emph{t} test of the adjusted correlation with \eqn{N - 3} degrees of
#' freedom (their Equation 5), is evidence that the relationship is not an
#' artifact of method variance; the test is applied by the user, since this
#' function works from the correlation matrix alone and does not take
#' \eqn{N}.
#'
#' @param R A correlation matrix among the substantive items.
#' @param marker_r The marker variable's (CMV) correlation. When
#'   \code{NULL} (default) the smallest positive off-diagonal correlation in
#'   \code{R} is used as the proxy marker.
#'
#' @details
#' Writing \eqn{r_M} for the marker (or proxy) correlation, each substantive
#' correlation is adjusted as
#' \eqn{r^{A}_{ij} = (r_{ij} - r_M) / (1 - r_M)} (Lindell & Whitney, 2001,
#' Equation 4). The CMV-adjusted correlation matrix is returned as the
#' \code{"adjusted"} attribute; the reported table summarizes the marker
#' correlation and the average absolute correlation before and after
#' adjustment.
#'
#' The method presumes the variables are reflected so that their
#' intercorrelations are positive; a negative substantive correlation is
#' pushed further from zero by the adjustment rather than attenuated, so
#' reverse-code as needed before adjusting.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with rows
#'   \code{marker_correlation}, \code{mean_abs_r_unadjusted}, and
#'   \code{mean_abs_r_adjusted} in the \code{value} column. The full adjusted
#'   correlation matrix is the \code{"adjusted"} attribute.
#'
#' @references
#' Lindell, M. K., & Whitney, D. J. (2001). Accounting for common method
#'   variance in cross-sectional research designs. \emph{Journal of Applied
#'   Psychology, 86}(1), 114--121. \doi{10.1037/0021-9010.86.1.114}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{common_method_single_factor}} for the single-factor
#'   screen.
#'
#' @family multivariate and latent variable methods
#'
#' @keywords multivariate
#'
#' @examples
#' R <- matrix(c(1, .5, .4, .5, 1, .45, .4, .45, 1), 3, 3,
#'             dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
#' res <- common_method_marker(R, marker_r = 0.10)
#' res
#' attr(res, "adjusted")
#'
#' @export
common_method_marker <- function(R, marker_r = NULL) {
  if (!is.matrix(R) || !isSymmetric(unname(R)) || nrow(R) < 2L ||
      any(abs(R) > 1 + 1e-8)) {
    stop("'R' must be a symmetric correlation matrix.", call. = FALSE)
  }
  off <- R[upper.tri(R)]
  if (is.null(marker_r)) {
    pos <- off[off > 0]
    if (!length(pos)) {
      stop("No positive correlation is available to serve as a proxy ",
           "marker; supply 'marker_r'.", call. = FALSE)
    }
    marker_r <- min(pos)
  }
  if (!is.numeric(marker_r) || length(marker_r) != 1L || is.na(marker_r) ||
      marker_r <= -1 || marker_r >= 1) {
    stop("'marker_r' must be a single number in (-1, 1).", call. = FALSE)
  }
  adj <- (R - marker_r) / (1 - marker_r)
  diag(adj) <- 1
  dimnames(adj) <- dimnames(R)
  out <- data.frame(
    term  = c("marker_correlation", "mean_abs_r_unadjusted", "mean_abs_r_adjusted"),
    value = c(marker_r, mean(abs(off)), mean(abs(adj[upper.tri(adj)]))),
    stringsAsFactors = FALSE, row.names = NULL
  )
  res <- .as_dmar_tbl(out)
  attr(res, "adjusted") <- adj
  res
}
