#' Correlation Matrix to Covariance Matrix Conversion
#'
#' Rescales a correlation matrix into the covariance matrix implied by
#' a set of standard deviations, the inverse of the standardization
#' that produces a correlation matrix from a covariance matrix. Useful
#' when a published article reports correlations and standard
#' deviations but an analysis needs the covariances.
#'
#' @param cor_mat The correlation matrix to be converted
#' @param sd A vector that contains the standard deviations of the variables in the correlation matrix
#' @param discrepancy A small nonnegative tolerance (near 0; default \code{1e-5}). A value on the main diagonal of the correlation matrix is treated as equal to 1 when it is within \code{discrepancy} of 1, that is, when \eqn{|d - 1| \le} \code{discrepancy}
#'
#' @details The correlation matrix to convert can be either symmetric or triangular. The covariance matrix returned is always a symmetric matrix.
#'
#' @return A square numeric matrix giving the covariance matrix
#'   implied by the supplied correlation matrix and standard
#'   deviations, with the same row / column names as \code{cor_mat}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' Cor.Mat <- rbind(c(1.0000, 0.8254, 0.4261, 0.6237, 0.5901, 0.1564, 0.1551),
#'               c(0.8254, 1.0000, 0.5583, 0.5967, 0.6692, 0.1877, 0.2246),
#'               c(0.4261, 0.5583, 1.0000, 0.4933, 0.4455, 0.1472, 0.3433),
#'               c(0.6237, 0.5967, 0.4933, 1.0000, 0.6403, 0.1160, 0.5316),
#'               c(0.5901, 0.6692, 0.4455, 0.6403, 1.0000, 0.3769, 0.5742),
#'               c(0.1564, 0.1877, 0.1472, 0.1160, 0.3769, 1.0000, 0.2833),
#'               c(0.1551, 0.2246, 0.3433, 0.5316, 0.5742, 0.2833, 1.0000))
#' colnames(Cor.Mat) <- rownames(Cor.Mat) <- c("rating", "complaints", "privileges",
#' "learning", "raises", "critical", "advance")
#'
#' SDs <- c(12.172562, 13.314757, 12.235430, 11.737013, 10.397226, 9.894908, 10.288706)
#' convert_cor_cov(cor_mat=Cor.Mat, sd=SDs)
#'
#' @note
#' The correlation matrix input should be a square matrix, and the length of \code{sd} should be equal to the number of variables in the correlation matrix (i.e., the number of rows/columns).
#' Sometimes the correlation matrix input may not have exactly 1's on the main diagonal, due to, e.g., rounding; \code{discrepancy} specifies the allowable discrepancy so that the function still considers the input as a correlation matrix and can proceed
#' (but the function does not change the numbers on the main diagonal).
#'
#' @keywords design
#'
#' @family parameterization conversions
#'
#' @export

convert_cor_cov <- function(cor_mat, sd, discrepancy = 1e-5) {
  if (dim(cor_mat)[1] != dim(cor_mat)[2]) stop("'cor_mat' should be a square matrix")

  n <- sqrt(length(cor_mat))
  if (n != length(sd)) stop("The length of 'sd' should be the same as the number of rows of 'cor_mat'")

  if (length(sd[sd > 0]) != n) stop("The elements in 'sd' should all be non-negative")

  if (isSymmetric(cor_mat)) {
    IS_symmetric <- TRUE
  } else {
    IS_symmetric <- FALSE
  }
  p <- dim(cor_mat)[1]
  q <- p * (p - 1) / 2
  if (isTRUE(all.equal(cor_mat[lower.tri(cor_mat)], rep(0, q))) || isTRUE(all.equal(cor_mat[upper.tri(cor_mat)], rep(0, q)))) {
    IS_triangular <- TRUE
  } else {
    IS_triangular <- FALSE
  }
  if (!IS_symmetric && !IS_triangular) stop("The object 'cor_mat' should be either a symmetric or a triangular matrix")

  cov_mat <- diag(sd) %*% cor_mat %*% diag(sd)
  colnames(cov_mat) <- rownames(cov_mat) <- colnames(cor_mat)
  return(cov_mat)
}
