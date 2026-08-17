# Check whether a contrast-coefficient matrix is mutually orthogonal.
#' Check Whether a Set of Contrasts Is Mutually Orthogonal
#'
#' Tests whether every pair of columns in a contrast-coefficient matrix
#' is orthogonal under either the equal-\eqn{n} convention
#' \eqn{\sum_i c_{ik} c_{ij} = 0} or the unequal-\eqn{n} convention
#' \eqn{\sum_i c_{ik} c_{ij} / n_i = 0} (Maxwell, Delaney, & Kelley,
#' 2027, Sec. 4.10; Kirk, 2013). Also checks that each column sums to
#' zero (the contrast property).
#'
#' @param contrasts A numeric \eqn{a \times m} matrix or
#'   \code{data.frame}, where \eqn{a} is the number of groups and
#'   \eqn{m} is the number of contrasts.
#' @param n Optional integer vector of length \eqn{a} giving the per-
#'   group sample sizes. If supplied, the unequal-\eqn{n} convention is
#'   used; otherwise the equal-\eqn{n} convention is assumed.
#' @param tol Numerical tolerance for declaring orthogonality. Default
#'   \code{1e-8}.
#'
#' @return A \code{data.frame} with rows for the overall
#'   orthogonality flag (\code{1} = all pairs orthogonal, \code{0} =
#'   not), the contrast-sum-to-zero flag, the number of contrasts
#'   tested, and one row per pairwise dot-product, named by contrast
#'   pair.
#'
#' @details
#' \strong{Equal-\eqn{n}.} Two contrasts \eqn{\mathbf c, \mathbf d} on
#' \eqn{a} groups of equal size are orthogonal iff
#' \eqn{\sum_{i=1}^{a} c_i d_i = 0}.
#'
#' \strong{Unequal-\eqn{n}.} With sample sizes \eqn{n_1, \ldots, n_a},
#' the orthogonality condition that yields uncorrelated sample
#' contrasts is \eqn{\sum_{i=1}^{a} c_i d_i / n_i = 0}.
#'
#' \strong{Useful for design checks.} Before performing planned
#' comparisons or partitioning the omnibus sums of squares, the user
#' typically wants confirmation that the chosen contrast set is
#' orthogonal so that its component SS sum to the omnibus SS.
#'
#' @references
#' Kirk, R. E. (2013). \emph{Experimental design: Procedures for the
#'   behavioral sciences} (4th ed.). Sage.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Sec. 4.10.)
#'
#' @seealso \code{\link{effects_coding}}, \code{\link{helmert_coding}},
#'   \code{\link{ci_scheffe}}
#'
#' @examples
#' # 1. Two orthogonal contrasts on a 4-group design (equal n):
#' cmat <- cbind(
#'   c_linear  = c(-3, -1,  1,  3),
#'   c_quad    = c( 1, -1, -1,  1)
#' )
#' is_orthogonal_set(cmat)
#'
#' # 2. Same contrasts under unequal sample sizes:
#' is_orthogonal_set(cmat, n = c(10, 8, 12, 9))
#'
#' # 3. Non-orthogonal pair:
#' cmat_bad <- cbind(
#'   c_diff_1  = c( 1, -1,  0,  0),
#'   c_diff_2  = c( 1,  0, -1,  0)
#' )
#' is_orthogonal_set(cmat_bad)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @family design utilities
#'
#' @export

is_orthogonal_set <- function(contrasts, n = NULL, tol = 1e-8) {
  if (is.data.frame(contrasts)) contrasts <- as.matrix(contrasts)
  if (!is.numeric(contrasts) || !is.matrix(contrasts))
    stop("'contrasts' must be a numeric matrix or data.frame.")
  a <- nrow(contrasts); m <- ncol(contrasts)
  if (a < 2L)  stop("Need at least 2 groups (rows).")
  if (m < 2L)  stop("Need at least 2 contrasts (columns).")
  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != a)
      stop(sprintf("'n' must be a numeric vector of length %d.", a))
    if (any(n <= 0)) stop("'n' must contain positive integers.")
  }
  if (is.null(colnames(contrasts)))
    colnames(contrasts) <- paste0("C", seq_len(m))

  sums_to_zero <- abs(colSums(contrasts)) < tol
  pairs <- utils::combn(m, 2L)

  dots <- numeric(ncol(pairs))
  labels <- character(ncol(pairs))
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1L, k]; j <- pairs[2L, k]
    c_i <- contrasts[, i]; c_j <- contrasts[, j]
    dot <- if (is.null(n)) sum(c_i * c_j) else sum(c_i * c_j / n)
    dots[k] <- dot
    labels[k] <- paste(colnames(contrasts)[i], colnames(contrasts)[j],
                       sep = " . ")
  }

  all_orthogonal <- all(abs(dots) < tol)
  all_sum_zero   <- all(sums_to_zero)

  pair_rows <- data.frame(
    term  = paste0("dot[", labels, "]"),
    value = dots, stringsAsFactors = FALSE)
  out <- rbind(
    data.frame(term  = c("all_orthogonal", "all_contrasts_sum_to_zero",
                          "n_contrasts"),
               value = c(as.integer(all_orthogonal),
                         as.integer(all_sum_zero), m),
               stringsAsFactors = FALSE),
    pair_rows
  )
  .as_dmar_tbl(out)
}
