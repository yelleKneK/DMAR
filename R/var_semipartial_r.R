#' Asymptotic Variance of the Semipartial (Part) Correlation Coefficient
#'
#' Computes the large-sample variance of the sample semipartial (also known
#' as the part) correlation coefficient \eqn{r_{Y(X \cdot Z_1 \cdots Z_J)}}
#' under multivariate normality. The function provides two variances: the
#' asymptotic variance under the alternative (the partial-correlation
#' form applied by analogy) and, when the full-model coefficient of
#' multiple determination \eqn{R^2_{Y \cdot X Z_1 \cdots Z_J}} is supplied,
#' the exact null-hypothesis variance derived from the multiple-regression
#' \emph{F}-test for the unique contribution of \eqn{X} (Cohen, Cohen, West,
#' & Aiken, 2003).
#'
#' @param r_sp Sample semipartial correlation coefficient
#'   \eqn{r_{Y(X \cdot Z_1 \cdots Z_J)}}, with the controlled variables
#'   partialled out of \eqn{X} only (not \eqn{Y}). Must be in \eqn{[-1, 1]}.
#' @param n Total sample size.
#' @param J Number of variables partialled out of \eqn{X} (i.e., the count
#'   of \eqn{Z_1, \ldots, Z_J}); must be at least 1. Defaults to 1.
#' @param R2_full Optional coefficient of multiple determination
#'   \eqn{R^2_{Y \cdot X Z_1 \cdots Z_J}} from the full model that includes
#'   \eqn{X} and all controls. When supplied, the null-hypothesis variance
#'   \eqn{(1 - R^2)/(n - J - 2)} is returned instead of the alternative-side
#'   asymptotic variance. Must be in \eqn{[0, 1]}.
#'
#' @return A one-row \code{data.frame} with columns \code{term}
#'   (either \code{"var_semipartial_r"} or
#'   \code{"var_semipartial_r_under_null"}) and \code{value} (the requested
#'   variance).
#'
#' @details
#' \strong{Background.} The squared semipartial \eqn{r^2_{Y(X \cdot Z)}}
#' equals the increase in \eqn{R^2} when \eqn{X} is added to a model already
#' containing the controls \eqn{Z_1, \ldots, Z_J}, i.e., the unique variance
#' in \eqn{Y} attributable to \eqn{X}. Unlike the partial, the semipartial
#' is on the original scale of \eqn{Y} rather than on the partialled scale,
#' which makes it the natural effect size companion to standardized
#' regression coefficients in multiple-regression reports (Cohen et al.,
#' 2003).
#'
#' \strong{Asymptotic variance (default).} Under multivariate normality the
#' semipartial admits the same large-sample form as the partial
#' (Fisher, 1924, applied by analogy):
#' \deqn{\mathrm{Var}(\hat r_{Y(X \cdot Z)}) \;\approx\;
#'         \frac{(1 - \rho^2_{Y(X \cdot Z)})^2}{n - J - 1}.}
#' The function evaluates this with \eqn{\hat r_{sp}} substituted for
#' \eqn{\rho_{sp}}. This is the appropriate quantity for Wald-style
#' inference and for AIPE-style precision planning analogous to that of the
#' partial correlation. Aloe and Becker (2012) develop the asymptotic
#' variance of the semipartial as a function of the full population
#' correlation structure, and Yuan and Chan (2011) give exact higher-order results
#' for the closely related standardized regression coefficients; the
#' present approximation matches the leading \eqn{1/n} behavior.
#'
#' \strong{Null-hypothesis variance (when \code{R2_full} is supplied).} In
#' multiple regression the unique contribution of \eqn{X} is tested with
#' \deqn{F \;=\; \frac{r^2_{Y(X \cdot Z)}\,(n - J - 2)}{1 - R^2_{Y \cdot X Z}}
#'         \;\sim\; F(1,\, n - J - 2)}
#' under \eqn{H_0\!: \rho_{Y(X \cdot Z)} = 0} (Cohen et al., 2003,
#' equation 3.7.3). Equivalently
#' \eqn{t = \hat r_{sp}\,\sqrt{(n - J - 2)/(1 - R^2_{Y \cdot X Z})}} is a
#' \emph{t}-statistic on \eqn{n - J - 2} degrees of freedom, so the
#' \emph{under-the-null} variance of \eqn{\hat r_{sp}} is
#' \deqn{\mathrm{Var}_0(\hat r_{sp}) \;=\;
#'         \frac{1 - R^2_{Y \cdot X Z}}{n - J - 2}.}
#' Supplying \code{R2_full} returns this null variance, which is the
#' standard ingredient for testing the significance of \eqn{X}'s unique
#' contribution.
#'
#' @references
#' Aloe, A. M., & Becker, B. J. (2012). An effect size for regression
#'   predictors in meta-analysis. \emph{Journal of Educational and Behavioral
#'   Statistics, 37}(2), 278--297. \doi{10.3102/1076998610396901}
#'
#' Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003).
#'   \emph{Applied multiple regression/correlation analysis for the
#'   behavioral sciences} (3rd ed.). Lawrence Erlbaum.
#'
#' Fisher, R. A. (1924). The distribution of the partial correlation
#'   coefficient. \emph{Metron, 3}, 329--332.
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 4 on contrasts,
#'   Chapter 5 on multiple comparisons, and Chapter 9 on ANCOVA.)
#'
#' Yuan, K.-H., & Chan, W. (2011). Biases and standard errors of standardized
#'   regression coefficients. \emph{Psychometrika, 76}(4), 670--690.
#'   \doi{10.1007/s11336-011-9224-6}
#'
#' @examples
#' # Olkin-Finn-style asymptotic variance of a semipartial r = .25 with
#' #     n = 100 and J = 3 controls in X.
#' var_semipartial_r(r_sp = 0.25, n = 100, J = 3)
#'
#' # With R^2 of the full model also supplied, the function returns the
#' #     null-hypothesis variance used in the F-test for X's unique
#' #     contribution.
#' var_semipartial_r(r_sp = 0.25, n = 100, J = 3, R2_full = 0.42)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{var_partial_r}}, \code{\link{var_R2}},
#'   \code{\link{ci_reg_coef}}
#'
#' @keywords design regression
#'
#' @export
var_semipartial_r <- function(r_sp, n, J = 1, R2_full = NULL) {
  if (!is.numeric(r_sp) || length(r_sp) != 1L || is.na(r_sp) ||
      abs(r_sp) > 1) {
    stop("'r_sp' must be a single number in [-1, 1].", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || n <= 0) {
    stop("'n' must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(J) || length(J) != 1L || J < 1L) {
    stop("'J' (number of partialled variables) must be a single integer >= 1.",
         call. = FALSE)
  }
  if (!is.null(R2_full)) {
    if (!is.numeric(R2_full) || length(R2_full) != 1L ||
        is.na(R2_full) || R2_full < 0 || R2_full > 1) {
      stop("'R2_full' must be a single number in [0, 1].", call. = FALSE)
    }
    if (n - J - 2 <= 0) {
      stop("Need n > J + 2 for the null-hypothesis variance.", call. = FALSE)
    }
    value <- (1 - R2_full) / (n - J - 2)
    term  <- "var_semipartial_r_under_null"
  } else {
    if (n - J - 1 <= 0) {
      stop("Need n > J + 1 for the asymptotic variance.", call. = FALSE)
    }
    value <- (1 - r_sp^2)^2 / (n - J - 1)
    term  <- "var_semipartial_r"
  }
  out <- data.frame(term = term, value = value,
                    stringsAsFactors = FALSE, row.names = NULL)
  .as_dmar_tbl(out)
}
