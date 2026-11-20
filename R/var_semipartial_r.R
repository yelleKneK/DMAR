#' Asymptotic Variance of the Semipartial (Part) Correlation Coefficient
#'
#' Computes a variance for the sample semipartial (also known as the part)
#' correlation coefficient \eqn{r_{Y(X \cdot Z_1 \cdots Z_{p-1})}} under
#' multivariate normality. The function provides two quantities, both of
#' which need the full-model coefficient of multiple determination
#' \eqn{R^2_{Y \cdot X Z_1 \cdots Z_{p-1}}}: the squared standard error in the
#' \emph{t} test of \eqn{X}'s coefficient (the default when \code{R2_full}
#' is supplied; Cohen, Cohen, West, & Aiken, 2003), and the first-order
#' normal-theory variance of the sample semipartial correlation
#' (\code{under_null = FALSE}), which is the quantity for a Wald confidence
#' interval and for sample size planning. Without \code{R2_full} the
#' function returns the first-order variance for the special case in which
#' the variables partialed out explain none of the outcome, and warns.
#'
#' @param r_sp Sample semipartial correlation coefficient
#'   \eqn{r_{Y(X \cdot Z_1 \cdots Z_{p-1})}}, with the other variables
#'   partialed out of \eqn{X} only (not \eqn{Y}). Must be in \eqn{[-1, 1]}.
#' @param N Total sample size.
#' @param p Number of predictors in the regression of \emph{Y} on \emph{X}
#'   and the variables partialed out, so \code{p - 1} variables are
#'   partialed out and \code{p} is at least 2. A categorical variable with
#'   \emph{a} levels counts as \emph{a} - 1 predictors. It is the \code{p}
#'   of \code{\link{ss_power_R2}}, \code{\link{ci_R2}}, and
#'   \code{\link{ci_reg_coef}}: the correlation concerns one predictor of a
#'   regression with \code{p} predictors. In the notation of a model
#'   comparison, \eqn{q} is the number of predictors removed from the full
#'   model to give the restricted model (\eqn{q = p_{full} - p_{restricted}},
#'   the predictors being tested), and the number of variables partialed out
#'   is \eqn{p - q}; this correlation concerns one predictor, so \eqn{q = 1}
#'   and \eqn{p - 1} variables are partialed out. Defaults to 2 (one variable
#'   partialed out).
#' @param R2_full Coefficient of multiple determination
#'   \eqn{R^2_{Y \cdot X Z_1 \cdots Z_{p-1}}} from the full model that includes
#'   \eqn{X} and every variable partialed out (\code{R2_Y_X} in
#'   \code{\link{ci_reg_coef}}; the squared semipartial correlation is
#'   \code{R2_full} minus the squared multiple correlation of the model
#'   without \emph{X}). Must be in \eqn{[0, 1]} and
#'   cannot be smaller than \code{r_sp} squared, because the squared
#'   semipartial correlation is the part of the full-model \eqn{R^2} that
#'   \eqn{X} adds; a smaller value is refused with an error of class
#'   \code{dmar_semipartial_r_incoherent}, unless the two values could
#'   have come from rounding to two decimals (see Details). When
#'   \code{NULL} (the default),
#'   the value returned is the first-order variance that holds when the
#'   variables partialed out explain none of the outcome, with a warning of
#'   class \code{dmar_semipartial_r_no_R2_full}.
#' @param under_null Logical, used only when \code{R2_full} is supplied.
#'   If \code{TRUE} (the default), the function returns
#'   \eqn{(1 - R^2)/(N - p - 1)}, the squared standard error in the
#'   \emph{t} test of \eqn{X}'s coefficient. If \code{FALSE}, it returns
#'   the first-order normal-theory variance of the sample semipartial
#'   correlation, evaluated at \code{r_sp} and \code{R2_full}.
#'
#' @return A one-row \code{data.frame} with columns \code{term}
#'   (either \code{"var_semipartial_r"} or
#'   \code{"var_semipartial_r_under_null"}) and \code{value} (the requested
#'   variance).
#'
#' @details
#' \strong{Background.} The squared semipartial \eqn{r^2_{Y(X \cdot Z)}}
#' equals the increase in \eqn{R^2} when \eqn{X} is added to a model already
#' containing the variables partialed out, \eqn{Z_1, \ldots, Z_{p-1}}, that is, the unique
#' variance in \eqn{Y} attributable to \eqn{X}. Unlike the partial, the
#' semipartial is on the original scale of \eqn{Y} rather than on the
#' partialed scale, which makes it the natural effect size companion to
#' standardized regression coefficients in multiple regression reports
#' (Cohen et al., 2003). Aloe and Becker (2012) propose it as an effect
#' size for regression predictors in meta-analysis.
#' \code{\link{ss_power_R2_change}} plans the power of the test of that
#' change in \eqn{R^2}.
#'
#' \strong{First-order variance (\code{under_null = FALSE}).} Write
#' \eqn{\rho_{sp}} for the population semipartial correlation,
#' \eqn{P^2} for the population full-model squared multiple correlation,
#' and \eqn{P_r^2 = P^2 - \rho_{sp}^2}, the squared multiple correlation of
#' \eqn{Y} with \eqn{Z_1, \ldots, Z_{p-1}} alone. Under multivariate normality,
#' to first order in \eqn{1/N},
#' \deqn{N\,\mathrm{Var}(\hat r_{sp}) \;\to\;
#'         (1 - \rho_{sp}^2)^2 \;-\; P_r^2\,(1 - 2 \rho_{sp}^2).}
#' The variance depends on the population through \eqn{\rho_{sp}} and
#' \eqn{P_r^2} only. Three cases show how it behaves. When the variables
#' partialed out explain none of the outcome (\eqn{P_r^2 = 0}) it is
#' \eqn{(1 - \rho_{sp}^2)^2}, the form of a simple or partial correlation.
#' When \eqn{\rho_{sp} = 0} it is \eqn{1 - P^2}. For \eqn{P_r^2 > 0}, it is
#' smaller than \eqn{(1 - \rho_{sp}^2)^2} when \eqn{\rho_{sp}^2 < 1/2} and
#' larger when \eqn{\rho_{sp}^2 > 1/2}. Because the variance depends on
#' \eqn{P_r^2}, no function of \eqn{r_{sp}} alone stabilizes it, and Fisher's
#' \emph{Z} transformation does not apply to a semipartial correlation as
#' it does to a partial correlation. The expression is the first-order
#' variance of Aloe and Becker (2012, Equation 4), which they obtain from
#' Alf and Graf (1999), with \eqn{N - p} in place of their \eqn{n}: their
#' numerator, \eqn{P^4 - 2 P^2 + P_{(f)}^2 + 1 - P_{(f)}^4} with
#' \eqn{P_{(f)}^2} the squared multiple correlation of the model without
#' \eqn{X}, equals the expression above with \eqn{P_r^2 = P_{(f)}^2}. The
#' result follows from the delta
#' method variance of the difference between two nested squared multiple
#' correlations (Model A of Olkin & Finn, 1995, for one variable partialed
#' out; Alf & Graf, 1999, for any number), because \eqn{r_{sp}^2} is that difference and
#' so \eqn{\mathrm{Var}(\hat r_{sp}) \approx
#'   \mathrm{Var}(\widehat{\Delta R^2}) / (4 \rho_{sp}^2)}.
#'
#' The function evaluates the expression at \code{r_sp} and
#' \code{R2_full} and divides by \eqn{N - p}. The number of variables
#' partialed out does not enter the first-order variance. The divisor
#' \eqn{N - p} is a convention, kept for agreement with
#' \code{\link{var_partial_r}} and because, with estimated values, a Wald
#' interval formed with it has coverage closer to the nominal level than
#' one formed with \eqn{N}; it is not the most accurate estimate of the
#' variance itself (in simulation neither \eqn{N - 1} nor \eqn{N - p}
#' was closer to the observed variance in every condition).
#'
#' The approximation has limits that simulation shows and that no theory
#' given here corrects: it understates the variance when
#' \eqn{|\rho_{sp}|} is large and \eqn{N - p} is small (in one simulated
#' condition, \eqn{\rho_{sp} = .95}, \eqn{P_r^2 = 0}, \eqn{p = 2}, and
#' \eqn{N = 40}, it was about 39 percent of the variance observed), and the
#' sample semipartial correlation is biased toward zero by roughly
#' \eqn{\rho_{sp} (p - 1) / (2N)}. A Wald interval formed from it also misses
#' the population value too often when the full-model \eqn{R^2} is near 1
#' and \eqn{N} is small, for any \eqn{\rho_{sp}} (in one simulated
#' condition, \eqn{\rho_{sp} = .1}, a full-model \eqn{R^2} of .99,
#' \eqn{p = 4}, and \eqn{N = 12}, it missed in .13 of samples). When
#' \eqn{|\rho_{sp}|} is above about .9 and the variables partialed out
#' explain some of the outcome, the interval misses more often on one side
#' than the other even with several hundred observations (in one simulated
#' condition, \eqn{\rho_{sp} = .95}, a full-model \eqn{R^2} of .97,
#' \eqn{p = 2}, and \eqn{N = 395}, it missed in .064 of samples, .051 of
#' them on one side).
#'
#' \strong{Values rounded to two decimals.} Reported values are often
#' rounded, and a rounded pair can have \code{R2_full} slightly below
#' \code{r_sp} squared: a semipartial correlation reported as .46 has a
#' square of .2116, and the full-model \eqn{R^2} beside it may be
#' reported as .21. When \code{R2_full + .005} is at least
#' \eqn{(|r_{sp}| - .005)^2} (with the difference in parentheses not below
#' 0), the pair is consistent with rounding. The function then gives a
#' warning of class \code{dmar_semipartial_r_rounded}, which carries the
#' two values, and proceeds: the squared semipartial correlation cannot
#' exceed the full-model \eqn{R^2}, so the variables partialed out are
#' taken to explain none of the outcome (\eqn{P_r^2 = 0}) in the first-order
#' variance, and the squared standard error of the \emph{t} test uses
#' \code{R2_full} as given. A pair further apart than rounding allows is
#' refused. \code{\link{ss_aipe_semipartial_r}} and
#' \code{\link{ss_aipe_semipartial_r_sensitivity}} apply the same rule to
#' their values.
#'
#' \strong{Without \code{R2_full}.} The function returns
#' \eqn{(1 - r_{sp}^2)^2 / (N - p)}, the first-order variance at
#' \eqn{P_r^2 = 0}. When the variables partialed out explain some of the
#' outcome, this value is too large if \eqn{\rho_{sp}^2 < 1/2} and too
#' small if \eqn{\rho_{sp}^2 > 1/2}, which is why the call warns. At
#' \eqn{r_{sp} = 0} it is \eqn{1/(N - p)} whatever the full-model
#' \eqn{R^2}, where the first-order variance is \eqn{(1 - P^2)/N}.
#'
#' \strong{Squared standard error of the \emph{t} test (\code{R2_full}
#' supplied, \code{under_null = TRUE}).} In multiple regression the unique
#' contribution of \eqn{X} is tested with
#' \deqn{F \;=\; \frac{r^2_{Y(X \cdot Z)}\,(N - p - 1)}{1 - R^2_{Y \cdot X Z}}
#'         \;\sim\; F(1,\, N - p - 1)}
#' under \eqn{H_0\!: \rho_{Y(X \cdot Z)} = 0} (Cohen et al., 2003).
#' Equivalently
#' \eqn{t = \hat r_{sp}\,\sqrt{(N - p - 1)/(1 - R^2_{Y \cdot X Z})}} is a
#' \emph{t}-statistic on \eqn{N - p - 1} degrees of freedom, and the
#' function returns the square of its denominator,
#' \deqn{\frac{1 - R^2_{Y \cdot X Z}}{N - p - 1},}
#' the squared standard error in the \emph{t} test of \eqn{X}'s
#' coefficient. It uses the sample \eqn{R^2} and equals the first-order
#' variance at \eqn{\rho_{sp} = 0}; it is the ingredient for testing the
#' significance of \eqn{X}'s unique contribution, not a variance for an
#' interval about a nonzero \eqn{r_{sp}}.
#'
#' With \code{under_null = TRUE} the value is
#' \eqn{(1 - R^2_{full})/(N - p - 1)}, which equals \eqn{r_{sp}^2/t^2} for
#' the \emph{t} statistic of the predictor's coefficient, by the identity
#' \eqn{r_{sp} = t\sqrt{(1 - R^2_{full})/(N - p - 1)}} (Aloe & Becker,
#' 2012, Equation 2).
#'
#' @references
#' Alf, E. F., Jr., & Graf, R. G. (1999). Asymptotic confidence limits for
#'   the difference between two squared multiple correlations: A simplified
#'   approach. \emph{Psychological Methods, 4}(1), 70--75.
#'   \doi{10.1037/1082-989X.4.1.70}
#'
#' Aloe, A. M., & Becker, B. J. (2012). An effect size for regression
#'   predictors in meta-analysis. \emph{Journal of Educational and Behavioral
#'   Statistics, 37}(2), 278--297. \doi{10.3102/1076998610396901}
#'
#' Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003).
#'   \emph{Applied multiple regression/correlation analysis for the
#'   behavioral sciences} (3rd ed.). Lawrence Erlbaum.
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
#' Olkin, I., & Finn, J. D. (1995). Correlations redux.
#'   \emph{Psychological Bulletin, 118}(1), 155--164.
#'   \doi{10.1037/0033-2909.118.1.155}
#'
#' @examples
#' # First-order variance of a semipartial r = .25 with N = 100, three
#' #     variables partialed out of X (p = 4 predictors in the full model),
#' #     and a full-model R^2 of .42.
#' var_semipartial_r(r_sp = 0.25, N = 100, p = 4, R2_full = 0.42,
#'                   under_null = FALSE)
#'
#' # By default, with R^2 of the full model supplied, the function returns
#' #     the squared standard error in the t test of X's unique
#' #     contribution.
#' var_semipartial_r(r_sp = 0.25, N = 100, p = 4, R2_full = 0.42)
#'
#' # A reported pair, r_sp = .46 with a full-model R^2 of .21, in which
#' #     rounding has left R^2 below .46^2 = .2116: the value is returned
#' #     with a warning that says how the pair was read. The handler keeps
#' #     the message and lets the call finish, so both are shown.
#' noted <- new.env()
#' rounded <- withCallingHandlers(
#'   var_semipartial_r(r_sp = 0.46, N = 100, p = 2, R2_full = 0.21),
#'   dmar_semipartial_r_rounded = function(w) {
#'     noted$message <- conditionMessage(w)
#'     invokeRestart("muffleWarning")
#'   })
#' rounded
#' noted$message
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{var_partial_r}}, \code{\link{var_R2}},
#'   \code{\link{ci_reg_coef}}, \code{\link{ss_aipe_semipartial_r}},
#'   \code{\link{ss_power_R2_change}}
#'
#' @keywords design regression
#'
#' @export
var_semipartial_r <- function(r_sp, N, p = 2, R2_full = NULL,
                              under_null = TRUE) {
  p <- .check_partial_p(p)
  # The arithmetic below is written in J, the number of variables partialed
  # out: N - J - 1 = N - p and N - J - 2 = N - p - 1.
  J <- p - 1
  if (!is.numeric(r_sp) || length(r_sp) != 1L || is.na(r_sp) ||
      abs(r_sp) > 1) {
    .dmar_refuse_argument(
      "'r_sp' must be a single number in [-1, 1].", arg = "r_sp")
  }
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0) {
    .dmar_refuse_argument(
      "'N' must be a single finite positive number.", arg = "N")
  }
  if (!is.logical(under_null) || length(under_null) != 1L ||
      is.na(under_null)) {
    .dmar_refuse_argument(
      "'under_null' must be TRUE or FALSE.", arg = "under_null")
  }
  if (is.null(R2_full)) {
    if (N - J - 1 <= 0) {
      .dmar_refuse_argument(
        "Need N > p for the asymptotic variance.", arg = NULL)
    }
    value <- (1 - r_sp^2)^2 / (N - J - 1)
    term  <- "var_semipartial_r"
    .dmar_warn(c(
      "Without 'R2_full' the value is the variance that holds when the ",
      "variables partialed out explain none of the outcome. When they ",
      "explain some of it, this value is too large if r_sp squared is ",
      "below 1/2 and too small if it is above. Supply 'R2_full' with ",
      "under_null = FALSE."),
      "semipartial_r_no_R2_full")
  } else {
    .check_semipartial_R2_full(r_sp, R2_full,
                               uses = if (under_null) "as_given" else
                                 "q_zero")
    if (under_null) {
      if (N - J - 2 <= 0) {
        .dmar_refuse_argument(
          "Need N > p + 1 for the squared standard error of the t test.",
          arg = NULL)
      }
      value <- (1 - R2_full) / (N - J - 2)
      term  <- "var_semipartial_r_under_null"
    } else {
      if (N - J - 1 <= 0) {
        .dmar_refuse_argument(
          "Need N > p for the asymptotic variance.", arg = NULL)
      }
      value <- .semipartial_r_var_numerator(r_sp, R2_full) / (N - J - 1)
      term  <- "var_semipartial_r"
    }
  }
  out <- data.frame(term = term, value = value,
                    stringsAsFactors = FALSE, row.names = NULL)
  .as_dmar_tbl(out)
}
