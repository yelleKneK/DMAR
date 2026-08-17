# Asymptotic variance of the coefficient of variation.
#' Asymptotic Variance of the Coefficient of Variation
#'
#' Computes the asymptotic variance of the sample coefficient of
#' variation \eqn{\hat\kappa = s / \bar Y} under normality, using
#' McKay's (1932) original noncentral \emph{t}-based approximation and
#' Vangel's (1996) refinement. Companion to \code{\link{ci_cv}} and
#' \code{\link{ss_aipe_cv}}.
#'
#' @param cv Population coefficient of variation
#'   \eqn{\kappa = \sigma / \mu}. Numeric scalar in \eqn{(0, \infty)}.
#'   When \eqn{\kappa} is very large (> 0.5), both McKay's and Vangel's
#'   approximations degrade and exact methods (\code{\link{ci_cv}})
#'   should be preferred.
#' @param n Sample size.
#'
#' @return A \code{data.frame} with rows for the McKay (1932) and
#'   Vangel (1996) approximations; columns are \code{term} and
#'   \code{value}.
#'
#' @details
#' \strong{McKay (1932).} The classical large-sample variance of the
#' sample CV under normality is
#' \deqn{\mathrm{Var}(\hat\kappa) \;\approx\;
#'   \frac{\kappa^2}{n - 1}
#'   \cdot \left(\frac{1}{2} + \kappa^2\right).}
#' This is exact up to \eqn{O(1/n)} and is what most planning tables
#' use. It begins to drift when \eqn{\kappa > 0.3} or so.
#'
#' \strong{Vangel (1996).} Vangel showed that a small-sample correction
#' that adjusts the McKay form for the noncentral \emph{t} mean factor
#' gives substantially better coverage of CIs derived from the variance:
#' \deqn{\mathrm{Var}_{\mathrm{Vangel}}(\hat\kappa) \;\approx\;
#'   \frac{\kappa^2}{n - 1}
#'   \cdot \left(\frac{1}{2} + \kappa^2 \cdot
#'   \frac{n + 1}{n - 1}\right).}
#' The two forms coincide in the large-\eqn{n} limit. We report both so
#' the user can see the magnitude of the small-sample correction.
#'
#' @references
#' Kelley, K. (2007). Sample size planning for the coefficient of variation
#'   from the accuracy in parameter estimation approach.
#'   \emph{Behavior Research Methods, 39}(4), 755--766.
#'   \doi{10.3758/BF03192966}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3.)
#'
#' McKay, A. T. (1932). Distribution of the coefficient of variation and
#'   the extended \emph{t} distribution. \emph{Journal of the Royal
#'   Statistical Society, 95}(4), 695--698.
#'
#' Vangel, M. G. (1996). Confidence intervals for a normal coefficient
#'   of variation. \emph{The American Statistician, 50}(1), 21--26.
#'   \doi{10.1080/00031305.1996.10473537}
#'
#' @seealso \code{\link{ci_cv}}, \code{\link{ss_aipe_cv}},
#'   \code{\link{ss_aipe_cv_sensitivity}}
#'
#' @examples
#' # 1. CV = 0.20 in a sample of 30:
#' var_cv(cv = 0.20, n = 30)
#'
#' # 2. The Vangel correction grows with cv (becomes non-trivial
#' #        for kappa > 0.3):
#' var_cv(cv = 0.50, n = 30)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family variance utilities
#'
#' @export

var_cv <- function(cv, n) {
  if (!is.numeric(cv) || length(cv) != 1L || cv <= 0)
    stop("'cv' must be a single positive number.")
  if (!is.numeric(n) || length(n) != 1L || n < 3)
    stop("'n' must be a single integer >= 3.")

  mckay  <- (cv^2 / (n - 1)) * (0.5 + cv^2)
  vangel <- (cv^2 / (n - 1)) * (0.5 + cv^2 * (n + 1) / (n - 1))

  out <- data.frame(
    term  = c("var_cv_mckay", "var_cv_vangel"),
    value = c(mckay, vangel),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out)
}
