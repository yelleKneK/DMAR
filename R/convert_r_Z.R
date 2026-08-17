#' Convert a Correlation Coefficient (\emph{r}) Into the Scale of Fisher's \emph{Z}
#'
#' This function converts a correlation coefficient into the scale of
#' Fisher's \emph{Z}, the variance-stabilizing transformation of a
#' correlation. Many authors call this map the \emph{z}-prime transform
#' and write the transformed value as \emph{z}'. The capital \emph{Z} is
#' meaningful: Fisher's \emph{Z} is not a \emph{z}-score (it is not a
#' standardized variate, that is, an observation centered and divided by a
#' standard deviation). It is the transform \eqn{Z = \mathrm{atanh}(r)} of
#' a correlation coefficient, applied because the sampling distribution of
#' \emph{Z} is approximately normal with a variance that does not depend on
#' the population correlation, which makes \emph{Z} convenient for forming
#' confidence intervals.
#'
#' @param r Correlation coefficient (between two variables)
#'
#' @return A 1-row \code{data.frame} with columns \code{term} and
#'   \code{value}. The \code{term} is \code{"Z_from_r"} and \code{value}
#'   is Fisher's \emph{Z} corresponding to the supplied correlation
#'   coefficient. The inverse direction is \code{\link{convert_Z_r}}.
#'
#' @details This function is typically used in the context of forming a confidence interval for a population correlation coefficient. Note that, in that situation, the two variables are assumed to follow a bivariate normal distribution (e.g., Hays, 1994).
#'
#' @examples
#' # From Hays (1994, pp. 649--650)
#' convert_r_Z(.35)
#'
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Hays, W. L. (1994). \emph{Statistics} (5th ed.). Fort Worth, TX:
#'   Harcourt Brace College Publishers.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{convert_Z_r}}, \code{\link{ci_r}}
#'
#' @keywords design
#'
#' @family parameterization conversions
#'
#' @rdname convert_r_Z
#'
#' @export

convert_r_Z <- function(r) {
  if (!is.numeric(r) || length(r) != 1L || is.na(r)) {
    stop("'r' must be a single numeric value. For a vector of ",
         "correlations, call atanh() directly.", call. = FALSE)
  }
  if (any(!is.na(r) & abs(r) >= 1)) {
    stop("'r' must satisfy |r| < 1; Fisher's Z is undefined at the ",
         "boundary and beyond. Got r = ",
         paste(r[abs(r) >= 1], collapse = ", "), ".", call. = FALSE)
  }
  term <- "Z_from_r"
  # atanh(r) is mathematically equivalent to (1/2) * log((1+r)/(1-r))
  # and is the numerically stable form.
  value <- atanh(r)
  .as_dmar_tbl(data.frame(term, value))
}

