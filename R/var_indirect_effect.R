# Variance of the mediated (indirect) effect ab.
#' Variance of the Mediated (Indirect) Effect \eqn{ab}
#'
#' Computes the asymptotic variance of the product of two regression
#' coefficients \eqn{\hat a \hat b} (the mediated/indirect effect in a
#' simple three-variable mediator model: \eqn{X \to M \to Y}) under four
#' competing formulas: Sobel (1982) first-order, Aroian (1947) /
#' Goodman (1960) second-order, and the full second-order delta method
#' with optional cross-product covariance. All four are reported in a
#' single output so the user can see the relative contributions of
#' the higher-order terms.
#'
#' @param reg_coef_a,reg_coef_b The regression coefficients \eqn{a} of
#'   \eqn{X \to M} and \eqn{b} of \eqn{M \to Y}: anticipated population
#'   values or estimates, on whatever scale the analysis uses (raw or
#'   standardized); the variances below must be on the same scale.
#' @param var_reg_coef_a,var_reg_coef_b Variances (squared standard errors)
#'   of \eqn{\hat a} and \eqn{\hat b}. For standardized regression with
#'   no covariates,
#'   \eqn{\mathrm{Var}(\hat a) \approx (1 - a^2)/(n - 2)} and
#'   \eqn{\mathrm{Var}(\hat b) \approx (1 - b^2)/\{(n - 3)(1 - a^2)\}},
#'   the standardized-model variance accounting for the correlation the
#'   \eqn{a} path induces between the predictors of \eqn{Y}.
#' @param cov_reg_coef_a_b Optional covariance between \eqn{\hat a} and
#'   \eqn{\hat b}. Defaults to 0 (the assumption underlying the standard
#'   Sobel formula). In practice the two estimators are nearly
#'   uncorrelated when the controls in the \eqn{Y}-on-\eqn{M} regression
#'   are uncorrelated with the \eqn{M}-on-\eqn{X} regression's
#'   predictors.
#'
#' @return A \code{data.frame} with rows for the four variance
#'   formulas; columns are \code{term} and \code{value}.
#'
#' @details
#' In the formulas below \eqn{a} and \eqn{b} are \code{reg_coef_a} and
#' \code{reg_coef_b}, \eqn{\mathrm{Var}(\hat a)} and
#' \eqn{\mathrm{Var}(\hat b)} are \code{var_reg_coef_a} and
#' \code{var_reg_coef_b}, and \eqn{\mathrm{Cov}(\hat a, \hat b)} is
#' \code{cov_reg_coef_a_b}.
#'
#' \strong{Sobel (1982) first-order.} The delta method variance of
#' \eqn{\hat a \hat b} under independent \eqn{\hat a}, \eqn{\hat b} is
#' \deqn{\mathrm{Var}_{\mathrm{Sobel}}(\hat a \hat b) \;=\;
#'   a^2 \mathrm{Var}(\hat b) + b^2 \mathrm{Var}(\hat a).}
#' This is the most cited form and is the variance used by the standard
#' Sobel \emph{z}-test (Sobel, 1982).
#'
#' \strong{Aroian (1947).} Aroian retains the second-order term:
#' \deqn{\mathrm{Var}_{\mathrm{Aroian}}(\hat a \hat b) \;=\;
#'   a^2 \mathrm{Var}(\hat b) + b^2 \mathrm{Var}(\hat a)
#'   + \mathrm{Var}(\hat a)\,\mathrm{Var}(\hat b).}
#' Aroian shows this is exact under joint normality of the two
#' independent estimators.
#'
#' \strong{Goodman (1960).} Goodman's "unbiased" variance subtracts the
#' second-order term instead of adding it:
#' \deqn{\mathrm{Var}_{\mathrm{Goodman}}(\hat a \hat b) \;=\;
#'   a^2 \mathrm{Var}(\hat b) + b^2 \mathrm{Var}(\hat a)
#'   - \mathrm{Var}(\hat a)\,\mathrm{Var}(\hat b).}
#' For small variances the three forms agree to leading order; they
#' diverge for noisy \eqn{\hat a}, \eqn{\hat b}.
#'
#' \strong{Second-order delta method (with covariance).} When
#' \eqn{\hat a} and \eqn{\hat b} share variability (e.g., they are both
#' estimated from the same regression of \eqn{Y} on \eqn{X} and \eqn{M}),
#' the cross-product covariance term enters:
#' \deqn{\mathrm{Var}(\hat a \hat b) \;\approx\;
#'   a^2 \mathrm{Var}(\hat b) + b^2 \mathrm{Var}(\hat a)
#'   + 2 a b \cdot \mathrm{Cov}(\hat a, \hat b).}
#' MacKinnon et al.\ (2002) show this matters in models with covariates
#' that simultaneously load on \eqn{M} and \eqn{Y}.
#'
#' \strong{Connection to \code{\link{ss_aipe_indirect_effect}}.} The
#' Sobel (delta method) variance is what
#' \code{\link{ss_aipe_indirect_effect}} builds on under
#' \code{method = "closed_form"} for AIPE planning; this function makes
#' the alternative formulas available for explicit comparison.
#'
#' @references
#' Aroian, L. A. (1947). The probability function of the product of two
#'   normally distributed variables. \emph{The Annals of Mathematical
#'   Statistics, 18}(2), 265--271.
#'
#' Goodman, L. A. (1960). On the exact variance of products.
#'   \emph{Journal of the American Statistical Association, 55}(292),
#'   708--713.
#'
#' Lachowicz, M. J., Preacher, K. J., & Kelley, K. (2018). A novel measure
#'   of effect size for mediation analysis.
#'   \emph{Psychological Methods, 23}(2), 244--261. \doi{10.1037/met0000165}
#'
#' MacKinnon, D. P., Lockwood, C. M., Hoffman, J. M., West, S. G., &
#'   Sheets, V. (2002). A comparison of methods to test mediation and
#'   other intervening variable effects. \emph{Psychological Methods,
#'   7}(1), 83--104. \doi{10.1037/1082-989X.7.1.83}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
#'   models: Quantitative strategies for communicating indirect effects.
#'   \emph{Psychological Methods, 16}(2), 93--115. \doi{10.1037/a0022658}
#'
#' Sobel, M. E. (1982). Asymptotic confidence intervals for indirect
#'   effects in structural equation models. \emph{Sociological
#'   Methodology, 13}, 290--312.
#'
#' Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
#'   analysis: Introducing the model-based constrained optimization
#'   procedure. \emph{Psychological Methods, 25}(4), 496--515.
#'   \doi{10.1037/met0000259}
#'
#' @seealso \code{\link{ss_aipe_indirect_effect}}
#'
#' @examples
#' # 1. Coefficients of 0.40 for both paths, each with a sampling variance
#' #    of 0.02, and no covariance:
#' var_indirect_effect(reg_coef_a = 0.40, reg_coef_b = 0.40,
#'                     var_reg_coef_a = 0.02, var_reg_coef_b = 0.02)
#'
#' # 2. With a positive covariance between a-hat and b-hat:
#' var_indirect_effect(reg_coef_a = 0.40, reg_coef_b = 0.40,
#'                     var_reg_coef_a = 0.02, var_reg_coef_b = 0.02,
#'                     cov_reg_coef_a_b = 0.005)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family variance utilities
#'
#' @export

var_indirect_effect <- function(reg_coef_a, reg_coef_b, var_reg_coef_a,
                                var_reg_coef_b, cov_reg_coef_a_b = 0) {
  if (!is.numeric(reg_coef_a) || length(reg_coef_a) != 1L)
    .dmar_refuse_argument(
      "'reg_coef_a' must be a single number.", arg = "reg_coef_a")
  if (!is.numeric(reg_coef_b) || length(reg_coef_b) != 1L)
    .dmar_refuse_argument(
      "'reg_coef_b' must be a single number.", arg = "reg_coef_b")
  if (!is.numeric(var_reg_coef_a) || length(var_reg_coef_a) != 1L ||
      var_reg_coef_a < 0)
    .dmar_refuse_argument(
      "'var_reg_coef_a' must be a single non-negative number.",
      arg = "var_reg_coef_a")
  if (!is.numeric(var_reg_coef_b) || length(var_reg_coef_b) != 1L ||
      var_reg_coef_b < 0)
    .dmar_refuse_argument(
      "'var_reg_coef_b' must be a single non-negative number.",
      arg = "var_reg_coef_b")
  if (!is.numeric(cov_reg_coef_a_b) || length(cov_reg_coef_a_b) != 1L)
    .dmar_refuse_argument("'cov_reg_coef_a_b' must be a single number.",
                          arg = "cov_reg_coef_a_b")

  sobel   <- reg_coef_a^2 * var_reg_coef_b + reg_coef_b^2 * var_reg_coef_a
  aroian  <- sobel + var_reg_coef_a * var_reg_coef_b
  goodman <- sobel - var_reg_coef_a * var_reg_coef_b
  delta_2 <- sobel + 2 * reg_coef_a * reg_coef_b * cov_reg_coef_a_b

  out <- data.frame(
    term  = c("var_sobel", "var_aroian", "var_goodman",
              "var_delta_second_order"),
    value = c(sobel, aroian, goodman, delta_2),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out)
}
