#' Convert Fisher's \emph{Z} Into the Scale of a Correlation Coefficient (\emph{r})
#'
#' Converts Fisher's \emph{Z} back into the scale of a correlation
#' coefficient (\emph{r}). Fisher's \emph{Z} is the variance-stabilizing
#' transformation of a correlation; many authors call it the \emph{z}-prime
#' transform and write the transformed value as \emph{z}'. The capital
#' \emph{Z} is meaningful: Fisher's \emph{Z} is not a \emph{z}-score (it is
#' not a standardized variate, that is, an observation centered and divided
#' by a standard deviation). This function applies the inverse transform
#' \eqn{r = \mathrm{tanh}(Z)} to return to the scale of a correlation
#' coefficient.
#'
#' @param Z Fisher's \emph{Z} (the variance-stabilizing transform of a
#'   correlation, which many authors call \emph{z}')
#'
#' @return A 1-row \code{data.frame} with columns \code{term} and
#'   \code{value}. The \code{term} is \code{"r_from_Z"} and \code{value}
#'   is the correlation coefficient corresponding to the supplied
#'   Fisher's \emph{Z}. The inverse direction is
#'   \code{\link{convert_r_Z}}.
#'
#' @details This function is typically used in the context of forming a confidence interval for a population correlation coefficient. Note that, in that situation, the two variables are assumed to follow a bivariate normal distribution (e.g., Hays, 1994).
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Hays, W. L. (1994). \emph{Statistics} (5th ed.). Fort Worth, TX:
#'   Harcourt Brace College Publishers.
#'
#' @examples
#' # From Hays (1994, pp. 649--650)
#' convert_Z_r(0.3654438)
#'
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{convert_r_Z}}, \code{\link{ci_r}}
#'
#' @keywords design
#'
#' @family parameterization conversions
#'
#' @rdname convert_Z_r
#'
#' @export

convert_Z_r <- function(Z) {
  if (!is.numeric(Z) || length(Z) != 1L || is.na(Z)) {
    stop("'Z' must be a single numeric value. For a vector of ",
         "Fisher's Z values, call tanh() directly.", call. = FALSE)
  }
  term <- "r_from_Z"
  # tanh(Z) is mathematically equivalent to (exp(2Z) - 1) / (exp(2Z) + 1)
  # and is numerically stable for large |Z|, where the explicit form
  # overflows to NaN around |Z| >= 350.
  value <- tanh(Z)
  .as_dmar_tbl(data.frame(term, value))
}

