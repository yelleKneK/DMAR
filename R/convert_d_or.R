#' Convert Between the Standardized Mean Difference and the Odds Ratio
#'
#' Invertible conversions between a two-group standardized mean difference
#' (Cohen's \emph{d}) and an odds ratio, by the logistic-distribution method
#' of Hasselblad and Hedges (1995): a continuous outcome split at a threshold
#' under logistic errors implies
#' \deqn{d = \log(\mathrm{OR}) \cdot \frac{\sqrt{3}}{\pi}, \qquad
#'       \mathrm{OR} = \exp\!\bigl(d \cdot \pi / \sqrt{3}\bigr).}
#' These conversions let binary-outcome studies enter a synthesis on the
#' standardized mean difference scale, or mean-difference studies enter one
#' on the odds ratio scale (Borenstein, Hedges, Higgins, & Rothstein, 2009,
#' Chapter 7).
#'
#' @param d The standardized mean difference.
#' @param or The odds ratio, a single positive number.
#'
#' @return A \code{data.frame} (class \code{dmar_tbl}) with a single
#'   row: term \code{odds_ratio} (for \code{convert_d_or}) or \code{smd}
#'   (for \code{convert_or_d}) and its \code{value}.
#'
#' @references
#' Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R.
#'   (2009). \emph{Introduction to meta-analysis}. Wiley.
#'
#' Hasselblad, V., & Hedges, L. V. (1995). Meta-analysis of screening and
#'   diagnostic tests. \emph{Psychological Bulletin, 117}(1), 167--178.
#'   \doi{10.1037/0033-2909.117.1.167}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{convert_d_r}} / \code{\link{convert_r_d}} for the
#'   correlation leg of the same triangle.
#'
#' @family parameterization conversions
#'
#' @keywords design
#'
#' @examples
#' # d = 0.5 corresponds to an odds ratio of about 2.48.
#' convert_d_or(d = 0.5)
#'
#' # And back, exactly.
#' convert_or_d(or = convert_d_or(d = 0.5)$value)
#'
#' # The null maps to the null: d = 0 is an odds ratio of 1.
#' convert_d_or(d = 0)
#'
#' @export
convert_d_or <- function(d) {
  if (!is.numeric(d) || length(d) != 1L || is.na(d)) {
    stop("'d' must be a single number.", call. = FALSE)
  }
  term  <- "odds_ratio"
  value <- exp(d * pi / sqrt(3))
  return(.as_dmar_tbl(data.frame(term, value)))
}

#' @rdname convert_d_or
#' @export
convert_or_d <- function(or) {
  if (!is.numeric(or) || length(or) != 1L || is.na(or) || or <= 0) {
    stop("'or' must be a single positive number.", call. = FALSE)
  }
  term  <- "smd"
  value <- log(or) * sqrt(3) / pi
  return(.as_dmar_tbl(data.frame(term, value)))
}
