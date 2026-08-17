# Asymptotic variance of the Pearson correlation.
#' Asymptotic Variance of the Pearson Correlation Coefficient
#'
#' Computes the asymptotic (large-sample) variance of the sample Pearson
#' product-moment correlation \eqn{r} under bivariate normality (Fisher,
#' 1915) and, optionally, the Bonett-Wright (2000) kurtosis-corrected
#' variance for non-normal margins. Stand-alone variance utility:
#' surprisingly absent from CRAN despite being a building block for
#' AIPE planning, meta-analytic weighting, and Wald-style inference on
#' \eqn{r}.
#'
#' @param rho Population correlation coefficient. Numeric scalar or
#'   vector in \eqn{(-1, 1)}.
#' @param n Total sample size on which the Pearson \eqn{r} would be
#'   computed. Scalar or vector.
#' @param kurtosis_x Optional excess kurtosis of the marginal
#'   distribution of \eqn{X}. When \code{kurtosis_x} and
#'   \code{kurtosis_y} are both supplied (along with \code{rho_xy_2x2y}
#'   if available), the Bonett-Wright (2000) corrected variance is
#'   returned in addition to the normal-theory variance.
#' @param kurtosis_y Optional excess kurtosis of \eqn{Y}; see
#'   \code{kurtosis_x}.
#'
#' @return A \code{data.frame} with the rows
#'   \itemize{
#'     \item \code{var_r_normal}, the normal-theory
#'       asymptotic variance \eqn{(1 - \rho^2)^2 / (n - 1)}
#'       (Hotelling, 1953, Section 7),
#'     \item \code{var_r_bonett_wright} (when kurtoses are supplied)
#'       the Bonett & Wright (2000) kurtosis-corrected variance,
#'     \item \code{var_fisher_z}, the variance of the Fisher
#'       \eqn{Z} transform, \eqn{1/(n - 3)}.
#'   }
#'
#' @details
#' \strong{Normal-theory variance (default).} Under bivariate normality
#' the large-sample variance is
#' \deqn{\mathrm{Var}(\hat r) \;\approx\; (1 - \rho^2)^2 / (n - 1),}
#' the leading term of the exact moment expansion (Hotelling, 1953,
#' Section 7; the exact density of \eqn{r} is Fisher's, 1915).
#' This is the workhorse variance and is exact in the limit; it is also
#' what Fisher's \eqn{Z} CI \code{\link{ci_r}} uses on the
#' transformed scale.
#'
#' \strong{Bonett-Wright kurtosis correction.} When the marginals are
#' \emph{not} normal, the asymptotic variance picks up a kurtosis-
#' dependent correction (Bonett & Wright, 2000):
#' \deqn{\mathrm{Var}(\hat r) \;\approx\; (1 - \rho^2)^2 / (n - 1)
#'   \cdot \bigl(1 + \rho^2 (\gamma_2^{(X)} + \gamma_2^{(Y)}) / 4\bigr),}
#' where \eqn{\gamma_2^{(X)}}, \eqn{\gamma_2^{(Y)}} are the excess kurtoses
#' of the two marginals. The correction is exact when the joint
#' distribution is elliptical; for non-elliptical joints it is a
#' first-order approximation. When the kurtosis arguments are
#' \code{NULL}, only the normal-theory variance is returned.
#'
#' \strong{Connection to Fisher's Z transform.} On the variance-
#' stabilized scale \eqn{Z = \tanh^{-1}(r)} the asymptotic variance is
#' \eqn{1/(n-3)} regardless of \eqn{\rho} (Fisher, 1921). This is
#' reported alongside the raw-scale variance because it is the natural
#' working scale for CI construction (\code{\link{ci_r}}).
#'
#' @references
#' Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for
#'   estimating Pearson, Kendall and Spearman correlations.
#'   \emph{Psychometrika, 65}(1), 23--28. \doi{10.1007/BF02294183}
#'
#' Fisher, R. A. (1915). Frequency distribution of the values of the
#'   correlation coefficient in samples from an indefinitely large
#'   population. \emph{Biometrika, 10}(4), 507--521.
#'
#' Fisher, R. A. (1921). On the "probable error" of a coefficient of
#'   correlation deduced from a small sample. \emph{Metron, 1}, 3--32.
#'
#' Hotelling, H. (1953). New light on the correlation coefficient and its
#'   transforms. \emph{Journal of the Royal Statistical Society, Series B,
#'   15}(2), 193--232. (Section 7 gives the exact moments of \eqn{r}.)
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge.
#'
#' Olkin, I., & Finn, J. D. (1995). Correlations redux.
#'   \emph{Psychological Bulletin, 118}(1), 155--164.
#'   \doi{10.1037/0033-2909.118.1.155}
#'
#' @seealso \code{\link{expected_r}}, \code{\link{ci_r}},
#'   \code{\link{var_partial_r}}, \code{\link{var_semipartial_r}}
#'
#' @examples
#' # 1. Normal-theory variance:
#' var_r(rho = 0.30, n = 50)
#'
#' # 2. With Bonett-Wright correction for leptokurtic margins
#' #        (excess kurtosis = 3 for each variable):
#' var_r(rho = 0.30, n = 50, kurtosis_x = 3, kurtosis_y = 3)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family variance utilities
#'
#' @export

var_r <- function(rho, n, kurtosis_x = NULL, kurtosis_y = NULL) {
  if (!is.numeric(rho)) stop("'rho' must be numeric.")
  if (any(abs(rho) >= 1, na.rm = TRUE))
    stop("'rho' must lie strictly inside (-1, 1).")
  if (!is.numeric(n) || any(n < 4, na.rm = TRUE))
    stop("'n' must be numeric and >= 4.")

  base <- (1 - rho^2)^2 / (n - 1)
  var_z <- 1 / (n - 3)

  rows <- data.frame(
    term  = c("var_r_normal", "var_fisher_z"),
    value = c(base, var_z),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (!is.null(kurtosis_x) && !is.null(kurtosis_y)) {
    if (!is.numeric(kurtosis_x) || !is.numeric(kurtosis_y))
      stop("'kurtosis_x' and 'kurtosis_y' must be numeric (excess kurtosis).")
    bw <- base * (1 + rho^2 * (kurtosis_x + kurtosis_y) / 4)
    rows <- rbind(rows,
                  data.frame(term  = "var_r_bonett_wright",
                              value = bw,
                              stringsAsFactors = FALSE,
                              row.names = NULL))
  }

  .as_dmar_tbl(rows)
}
