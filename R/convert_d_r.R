#' Convert Between the Standardized Mean Difference and the Correlation
#'
#' Invertible conversions between a two-group standardized mean difference
#' (Cohen's \emph{d}) and the (point-biserial) correlation between the
#' outcome and group membership. \code{convert_d_r()} maps \emph{d} to
#' \emph{r}; \code{convert_r_d()} maps \emph{r} back to \emph{d}. These are
#' the standard conversions used to bring effect sizes reported in different
#' metrics onto a common scale, for example when synthesizing a literature in
#' which some studies report mean differences and others report correlations
#' (Borenstein, Hedges, Higgins, & Rothstein, 2009, Chapter 7).
#'
#' @param d The standardized mean difference.
#' @param r The point-biserial correlation, in \eqn{(-1, 1)}.
#' @param n_1,n_2 Optional per-group sample sizes. When supplied, the
#'   conversion uses the unequal-group factor
#'   \eqn{a = (n_1 + n_2)^2 / (n_1 n_2)}; when omitted, equal group sizes are
#'   assumed, for which \eqn{a = 4}.
#'
#' @details
#' With \eqn{a = (n_1 + n_2)^2/(n_1 n_2)} (equal to 4 for equal groups), the
#' two directions are
#' \deqn{r = \frac{d}{\sqrt{d^2 + a}}, \qquad
#'       d = \frac{\sqrt{a}\, r}{\sqrt{1 - r^2}},}
#' exact inverses of one another for a given \eqn{a}. The same \code{n_1} and
#' \code{n_2} must be supplied to both directions for the round trip to be
#' exact.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with a single row:
#'   term \code{r} (for \code{convert_d_r}) or \code{smd} (for
#'   \code{convert_r_d}) and its \code{value}.
#'
#' @references
#' Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R.
#'   (2009). \emph{Introduction to meta-analysis}. Wiley.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{convert_d_or}} / \code{\link{convert_or_d}} for the
#'   odds ratio leg of the same triangle; \code{\link{smd}} and
#'   \code{\link{ci_r}} for estimating the quantities being converted.
#'
#' @family parameterization conversions
#'
#' @keywords design
#'
#' @examples
#' # Equal groups: d = 0.5 corresponds to r about .243.
#' convert_d_r(d = 0.5)
#'
#' # And back, exactly.
#' convert_r_d(r = convert_d_r(d = 0.5)$value)
#'
#' # Unequal groups change the conversion factor.
#' convert_d_r(d = 0.5, n_1 = 20, n_2 = 80)
#'
#' @export
convert_d_r <- function(d, n_1 = NULL, n_2 = NULL) {
  if (!is.numeric(d) || length(d) != 1L || is.na(d)) {
    stop("'d' must be a single number.", call. = FALSE)
  }
  a <- .convert_d_r_factor(n_1, n_2)
  term  <- "r"
  value <- d / sqrt(d^2 + a)
  return(.as_dmar_tbl(data.frame(term, value)))
}

#' @rdname convert_d_r
#' @export
convert_r_d <- function(r, n_1 = NULL, n_2 = NULL) {
  if (!is.numeric(r) || length(r) != 1L || is.na(r) || abs(r) >= 1) {
    stop("'r' must be a single correlation in (-1, 1).", call. = FALSE)
  }
  a <- .convert_d_r_factor(n_1, n_2)
  term  <- "smd"
  value <- sqrt(a) * r / sqrt(1 - r^2)
  return(.as_dmar_tbl(data.frame(term, value)))
}

# Shared unequal-group factor a = (n_1 + n_2)^2 / (n_1 n_2); 4 when group
# sizes are not supplied (equal groups). Not exported.
.convert_d_r_factor <- function(n_1, n_2) {
  if (is.null(n_1) != is.null(n_2)) {
    stop("Supply both 'n_1' and 'n_2', or neither.", call. = FALSE)
  }
  if (is.null(n_1)) return(4)
  for (nm in c("n_1", "n_2")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val < 1 || val != round(val)) {
      stop(sprintf("'%s' must be a single positive integer.", nm),
           call. = FALSE)
    }
  }
  (n_1 + n_2)^2 / (n_1 * n_2)
}
