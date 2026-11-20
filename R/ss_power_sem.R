#' Sample Size Planning for Structural Equation Modeling From the Power Analysis Perspective
#'
#' @description
#' Calculate the necessary sample size for an SEM study, so as to have enough power to reject the
#' null hypothesis that (a) the model has perfect fit, or (b) the difference in fit between two
#' nested models equal some specified amount.
#'
#' @param true_F_ML The true maximum likelihood fit function value in the population for the model of interest. Leave this argument NULL if you are doing nested model significance tests
#' @param df The degrees of freedom of the model of interest. Leave this argument NULL if you are doing nested model significance tests
#' @param null_rmsea The model's population RMSEA under the null hypothesis. Leave this argument NULL if you are doing nested model significance tests
#' @param true_rmsea The model's population RMSEA under the alternative hypothesis. This should be the model's true population RMSEA value. Leave this argument NULL if you are doing nested model significance tests
#' @param true_F_full The population maximum likelihood fit function value for the full model
#' @param true_F_restricted The population maximum likelihood fit function value for the restricted model
#' @param true_rmsea_full The population RMSEA value for the full model
#' @param true_rmsea_restricted The population RMSEA value for the restricted model
#' @param df_full The degrees of freedom for the full model
#' @param df_restricted The degrees of freedom for the restricted model
#' @param alpha_level The Type I error rate. Defaults to 0.05.
#' @param desired_power The desired power. Defaults to 0.85,
#'   matching the rest of the \code{ss_power_*} family.
#'
#' @return A \code{data.frame} with a \code{necessary_N} row, the smallest
#'   integer \emph{N} whose power reaches \code{desired_power} under the
#'   supplied fit-function or RMSEA alternative, and an \code{actual_power} row
#'   giving the realized power at that \emph{N}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # One-model test: necessary N to reject H0: RMSEA = 0 in favor of
#' # a model whose population RMSEA is 0.05 at 80% power, alpha = .05,
#' # with df = 20.
#' ss_power_sem(null_rmsea = 0, true_rmsea = 0.05, df = 20,
#'              alpha_level = 0.05, desired_power = 0.80)
#'
#' # Equivalent input via the population fit function, which is the
#' # degrees of freedom times the squared population RMSEA.
#' ss_power_sem(true_F_ML = 20 * 0.05^2, df = 20, alpha_level = 0.05, desired_power = 0.80)
#'
#' # Two-model nested test: necessary N to detect the difference
#' # between a full model (RMSEA = 0.04, df = 18) and a restricted
#' # model (RMSEA = 0.06, df = 22) at 80% power.
#' ss_power_sem(true_rmsea_full = 0.04, df_full = 18,
#'              true_rmsea_restricted = 0.06, df_restricted = 22,
#'              alpha_level = 0.05, desired_power = 0.80)
#'
#' @references
#' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
#'   analysis and determination of sample size for covariance structure
#'   modeling. \emph{Psychological Methods, 1}(2), 130--149.
#'   \doi{10.1037/1082-989X.1.2.130}
#'
#' Lai, K., & Kelley, K. (2011). Accuracy in parameter estimation for
#'   targeted effects in structural equation modeling: Sample size
#'   planning for narrow confidence intervals.
#'   \emph{Psychological Methods, 16}(2), 127--148. \doi{10.1037/a0021764}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @keywords multivariate design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export


ss_power_sem <- function(true_F_ML = NULL, df = NULL, null_rmsea = NULL, true_rmsea = NULL, true_F_full = NULL, true_F_restricted = NULL, true_rmsea_full = NULL, true_rmsea_restricted = NULL,
                         df_full = NULL, df_restricted = NULL, alpha_level = 0.05, desired_power = 0.85) {
  .dmar_check_numbers("ss_power_sem", c("df", "null_rmsea", "true_rmsea"))
  .dmar_check_desired_power(desired_power)
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L ||
      !is.finite(alpha_level) || alpha_level <= 0 || alpha_level >= 1)
    .dmar_refuse_argument(
      "'alpha_level' must be a single number strictly between 0 and 1.",
      arg = "alpha_level")

  input_F_ml <- input_RMSEA_true <- input_F_full <- input_RMSEA_full <- 0
  One_model <- Two_model <- FALSE

  if (!is.null(true_F_ML)) {
    if (is.null(df)) .dmar_refuse_argument(
      "Because 'true_F_ML' is specified, the function requires the model's degrees of freedom 'df' as input as well.",
      arg = c("true_F_ML", "df"))
    rmsea0 <- 0
    rmseaa <- sqrt(true_F_ML / df)
    d <- df
    input_F_ml <- 1
    One_model <- TRUE
  } else {
    input_F_ml <- 0
  }

  if (!is.null(true_rmsea)) {
    if (is.null(null_rmsea) || is.null(df)) .dmar_refuse_argument(
      "Because 'true_rmsea' is specified, the function requires 'null_rmsea' and 'df' as input as well.",
      arg = c("true_rmsea", "null_rmsea", "df"))
    rmsea0 <- null_rmsea
    rmseaa <- true_rmsea
    d <- df
    input_RMSEA_true <- 1
    One_model <- TRUE
  } else {
    input_RMSEA_true <- 0
  }

  if (!is.null(true_F_full)) {
    if (is.null(df_full) || is.null(df_restricted) || is.null(true_F_restricted)) .dmar_refuse_argument(
      "Because 'true_F_full' is specified, the function requires 'true_F_restricted', 'df_full', and 'df_restricted' as input as well",
      arg = c("true_F_full", "true_F_restricted", "df_full"))
    if (!is.null(true_F_ML) || !is.null(true_rmsea_full) || !is.null(true_rmsea)) .dmar_refuse_argument(
      "Because 'true_F_full' is specified, do not put in 'true_F_ML' or any RMSEA values",
      arg = c("true_F_full", "true_F_ML"))
    if (true_F_full > true_F_restricted) .dmar_refuse_argument(
      "The full model's fit-function value 'true_F_full' must be smaller than or equal to the restricted model's 'true_F_restricted'.",
      arg = c("true_F_full", "true_F_restricted"))
    rmseaa <- sqrt(true_F_restricted / df_restricted)
    rmseab <- sqrt(true_F_full / df_full)
    da <- df_restricted
    db <- df_full
    input_F_full <- 1
    Two_model <- TRUE
  } else {
    input_F_full <- 0
  }

  if (!is.null(true_rmsea_full)) {
    if (is.null(true_rmsea_restricted) || is.null(df_full) || is.null(df_restricted)) .dmar_refuse_argument(
      "Because 'true_rmsea_full' is specified, the function also requires 'true_rmsea_restricted', 'df_full', and 'df_restricted' as input",
      arg = c("true_rmsea_full", "true_rmsea_restricted", "df_full"))
    if (true_rmsea_full > true_rmsea_restricted) .dmar_refuse_argument(
      "The full model's RMSEA must be smaller than or equal to the restricted model's RMSEA",
      arg = NULL)
    rmseaa <- true_rmsea_restricted
    rmseab <- true_rmsea_full
    da <- df_restricted
    db <- df_full
    input_RMSEA_full <- 1
    Two_model <- TRUE
  } else {
    input_RMSEA_full <- 0
  }

  if (input_F_ml + input_RMSEA_true + input_F_full + input_RMSEA_full == 0) .dmar_refuse_argument(
    "Please specify at least one valid input set: either ('true_F_ML', 'df') or ('null_rmsea', 'true_rmsea', 'df') for a one-model test, or ('true_F_full', 'true_F_restricted', 'df_full', 'df_restricted') or ('true_rmsea_full', 'true_rmsea_restricted', 'df_full', 'df_restricted') for a nested two-model test.",
    arg = c("true_F_ML", "df", "null_rmsea"))

  alpha   <- alpha_level
  desired <- desired_power

  # Power is monotone increasing in N under both frameworks, so the necessary
  # sample size is the smallest integer N whose power reaches 'desired'. Given a
  # power function of N, bracket the crossing with a coarse jump, then step back
  # to the exact minimal integer. This returns an N that actually delivers the
  # requested power (the earlier bisect-to-tolerance-then-ceiling scheme could
  # stop up to .001 short of 'desired' and then ceiling to a non-minimal or
  # under-powered N).
  min_n_for_power <- function(power_at) {
    n <- 100
    while (power_at(n) < desired) {
      n <- n + 100
      if (n > 1e7) .dmar_refuse(
        "Failed to reach 'desired_power' within a reasonable sample size.",
        "computation_failed")
    }
    # 'n' now meets the target; walk down to the smallest integer that still does.
    while (n > 1 && power_at(n - 1) >= desired) n <- n - 1
    n
  }

  if (One_model) {
    power_at <- function(n) {
      ncp0 <- (n - 1) * d * rmsea0^2
      ncpa <- (n - 1) * d * rmseaa^2
      if (rmsea0 < rmseaa) {
        cval <- qchisq(alpha, d, ncp = ncp0, lower.tail = FALSE)
        pchisq(cval, d, ncp = ncpa, lower.tail = FALSE)
      } else {
        cval <- qchisq(1 - alpha, d, ncp = ncp0, lower.tail = FALSE)
        1 - pchisq(cval, d, ncp = ncpa, lower.tail = FALSE)
      }
    }
    minn <- min_n_for_power(power_at)
  }

  ###########################################
  if (Two_model) {
    ddiff <- da - db
    fa <- da * rmseaa^2
    fb <- db * rmseab^2

    power_at <- function(n) {
      ncp <- (n - 1) * (fa - fb)
      cval <- qchisq(alpha, ddiff, ncp = 0, lower.tail = FALSE)
      pchisq(cval, ddiff, ncp = ncp, lower.tail = FALSE)
    }
    minn <- min_n_for_power(power_at)
  }

  # Report the realized power at the resolved N alongside the sample size, so
  # the achieved power travels with the plan and broom's tidy()/glance() find a
  # power row (the "actual_power" term the ss_power_* schema recognizes).
  out <- data.frame(term  = c("necessary_N", "actual_power"),
                    value = c(minn, power_at(minn)))
  class(out) <- c("dmar_ss_power", "dmar_tbl", "data.frame")
  return(out)
}
