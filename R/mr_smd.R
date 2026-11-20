#' Minimum Risk Point Estimation of the Population Standardized Mean Difference
#'
#' @description
#' A function for the sequential estimation of the standardized mean
#' difference with minimum risk. The function implements the procedure of
#' Chattopadhyay and Kelley (2017), which considers study cost and the
#' accuracy of the estimated standardized mean difference together.
#' \code{mr_smd} assumes normally distributed data, equal sample sizes in
#' the two groups, and an equal cost of sampling per observation in each
#' group.
#'
#' @param A The structural cost divided by the square of \code{epsilon}: the price one is willing to pay for a maximum allowable difference of \code{epsilon} between the estimate of the standardized mean difference and its parameter, a single positive number. Give either \code{A} or both \code{structural_cost} and \code{epsilon}.
#' @param structural_cost The structural cost of what one is willing to pay in a study, a single positive number.
#' @param epsilon The maximum desired difference between the estimated standardized mean difference and the population value, a single positive number.
#' @param d The current estimate of the standardized mean difference, a single finite number. It is needed when \code{pilot = FALSE}.
#' @param n_per_group Current sample size \emph{per group} (thus total sample size is \eqn{2n}), a whole number of at least 2; requires equal sample size \emph{per group}. It is needed when \code{pilot = FALSE}.
#' @param sampling_cost The cost of collecting one more observation, a single positive number. For example, if each survey costs 10 dollars to distribute and score, \code{sampling_cost} would be 10 dollars per additional observation.
#' @param pilot \code{TRUE} to plan the pilot sample size per group, or \code{FALSE} (the default) to evaluate the risk function and the stopping rule at the current sample size.
#' @param m0 The smallest pilot sample size per group allowed: a whole number of at least 4 (the default).
#' @param gamma The exponent \eqn{\gamma} of the correction term \eqn{n^{-\gamma}}{n^(-gamma)} in the stopping rule, which keeps sampling from stopping too early: a single number strictly between 0 and 1/2. Chattopadhyay and Kelley (2017) take \eqn{0 < \gamma < 1/2} and suggest .49 for practical purposes, the default.
#'
#' @details
#' Chattopadhyay and Kelley (2017) estimate the population standardized
#' mean difference sequentially, minimizing a risk function that adds the
#' mean square error of the estimator, weighted by \code{A}, to the total
#' sampling cost. With \eqn{n} observations per group and
#' \eqn{V_n^2 = 2 + d^2/4}, sampling stops at the first \eqn{n} for which
#' \eqn{n \ge \sqrt{A / (2c)}\,(V_n + n^{-\gamma})}{n >= sqrt(A / (2c)) (V_n + n^(-gamma))}, where \eqn{c} is
#' \code{sampling_cost}. With \code{pilot = FALSE} the function evaluates
#' that risk and the stopping rule at the current sample size; with
#' \code{pilot = TRUE} it returns the pilot sample size per group,
#' \eqn{\max(m_0, \lceil (A / (2c))^{1/(2 + 2\gamma)} \rceil)}{max(m0, ceiling((A / (2c))^(1 / (2 + 2 gamma))))}.
#' A formula value that is itself a whole number is that number, although
#' R can compute it a little above (with \eqn{A = 1{,}000{,}000}{A = 1,000,000},
#' \eqn{c = 5}, and \eqn{\gamma = .25} it is 100, not 101; see
#' \code{\link{mr_cv}}).
#' See Kelley and
#' Rausch (2006) for additional information on the standardized mean
#' difference.
#'
#' @return
#' A \code{dmar_tbl} with columns \code{term} and \code{value}. With
#' \code{pilot = FALSE} the rows are
#' \item{risk}{The estimated risk at the current sample size,
#'   \eqn{A V_n^2 / n + 2cn}{A V_n^2 / n + 2 c n}.}
#' \item{n_1}{Sample size for group 1 (the input \code{n_per_group}).}
#' \item{n_2}{Sample size for group 2 (the input \code{n_per_group}).}
#' \item{d}{The input estimate of the standardized mean difference.}
#' The attribute \code{is_satisfied} is \code{TRUE} when the stopping rule
#' is satisfied, so that sampling can stop, and the attribute
#' \code{criterion} holds
#' \eqn{\sqrt{A / (2c)}\,(V_n + n^{-\gamma})}{sqrt(A / (2c)) (V_n + n^(-gamma))},
#' the value that the sample size per group is compared with. The decision
#' is printed under the table: whether the stopping rule is met, the two
#' numbers compared, and whether to stop or to collect another observation
#' in each group. With \code{pilot = TRUE} the single row is
#' \code{pilot_n_per_group} (see Note).
#'
#' @references
#' Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
#'   variation with minimum risk: A sequential method for minimizing
#'   sampling error and study cost.
#'   \emph{Multivariate Behavioral Research, 51}(5), 627--648.
#'   \doi{10.1080/00273171.2016.1203279}
#'
#' Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
#'   difference with minimum risk: Maximizing accuracy and minimizing cost
#'   with sequential estimation.
#'   \emph{Psychological Methods, 22}(1), 94--113. \doi{10.1037/met0000089}
#'
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#'   parameter estimation for a general class of effect sizes: A
#'   sequential approach. \emph{Psychological Methods, 23}(2), 226--243.
#'   \doi{10.1037/met0000127}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via narrow
#'   confidence intervals.
#'   \emph{Psychological Methods, 11}(4), 363--385.
#'   \doi{10.1037/1082-989X.11.4.363}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' When \code{pilot=TRUE} the function returns one row, \code{pilot_n_per_group}: the size of the pilot sample \emph{per group}
#' that should be used (thus, the total sample size is twice the pilot sample size).
#'
#' @seealso
#' \code{\link{ci_smd}}, \code{\link{mr_cv}}
#'
#' @examples
#' # To obtain pilot sample size in a situation in which A=10000. Note that 'A' is
#' # 'structural_cost' divided by the square of 'epsilon'.
#'
#' # From Chattopadhyay and Kelley (2017)
#' mr_smd(pilot = TRUE, A = 10000, sampling_cost = 2.4, gamma = .49)
#'
#' High.SLS <- c(11, 7, 22, 13, 6, 9, 11, 16, 12, 17, 14, 8, 16)
#' Low.SLS <- c(3, 6, 10, 8, 14, 5, 12, 10, 6, 8, 13, 5, 9)
#'
#' mr_smd(d = 1.021484, n_per_group = 13, A = 10000, sampling_cost = 2.40, gamma = .49)
#'
#' # Or, using the smd() function:
#' mr_smd(d = smd(group_1 = High.SLS, group_2 = Low.SLS)[1,2], n_per_group = 13, A = 10000,
#'        sampling_cost = 2.40, gamma = .49)
#'
#' # Here, for this situation, the stopping rule is satisfied:
#' mr_smd(d = 1.00, n_per_group = 75, A = 10000, sampling_cost = 2.40, gamma = .49)
#'
#' @keywords design misc htest
#'
#' @export


mr_smd <- function(A, structural_cost, epsilon, d, n_per_group, sampling_cost, pilot = FALSE, m0 = 4, gamma = .49) {
  .dmar_check_numbers(
    "mr_smd", c("A", "structural_cost", "epsilon", "sampling_cost", "gamma"),
    required = c("A", "structural_cost", "epsilon", "sampling_cost", "gamma"))
  .mr_check_flag(pilot, "pilot")
  A <- .mr_resolve_A(A, structural_cost, epsilon,
                     given = c(A = !missing(A),
                               structural_cost = !missing(structural_cost),
                               epsilon = !missing(epsilon)))
  if (missing(sampling_cost)) .dmar_refuse_argument(c(
    "'sampling_cost' must be given: the cost of collecting one more ",
    "observation, a single positive finite number. The stopping rule and ",
    "the pilot sample size both divide 'A' by twice its value."),
    arg = "sampling_cost")
  .mr_check_positive(sampling_cost, "sampling_cost",
                     "the cost of collecting one more observation")
  .mr_check_gamma(gamma, paste0(
    "The stopping rule is defined for an exponent in that open interval ",
    "(Chattopadhyay & Kelley, 2017, Equation 15)."))
  m0 <- .dmar_check_count(m0, "m0", min = 4L, reason = paste0(
    "The procedure takes a pilot sample of at least 4 per group ",
    "(Chattopadhyay & Kelley, 2017, Equation 14)."))
  .mr_check_ratio(A, sampling_cost)

  if (pilot) {
    out <- data.frame(
      term = "pilot_n_per_group",
      value = .mr_pilot_size((A / (2 * sampling_cost))^(1 / (2 + 2 * gamma)), m0))
    return(.as_dmar_tbl(out, subclass = "dmar_mr"))
  }

  if (missing(d)) .dmar_refuse_argument(c(
    "'d' must be given when 'pilot' is FALSE: the current estimate of the ",
    "standardized mean difference, from which the stopping rule is ",
    "evaluated. For the pilot sample size, which needs no data, set ",
    "pilot = TRUE."), arg = "d")
  if (!(is.numeric(d) && length(d) == 1L && is.finite(d)))
    .dmar_refuse_argument(c(
      "'d' must be a single finite number; ", .dmar_given(d), ". It is the ",
      "current estimate of the standardized mean difference."),
      arg = "d", value = d)
  if (missing(n_per_group)) .dmar_refuse_argument(c(
    "'n_per_group' must be given when 'pilot' is FALSE: the number of ",
    "observations collected so far in each of the two groups."),
    arg = "n_per_group")
  n <- .dmar_check_count(n_per_group, "n_per_group", min = 2L, reason = paste0(
    "It is the number of observations collected so far in each of the two ",
    "groups, which are of equal size."))

  # V_n^2 = 2 + d^2 / 4 is at least 2, so the rule never meets a variance
  # estimate it cannot use.
  criterion <- sqrt(A / (2 * sampling_cost)) * (sqrt(2 + d^2 / 4) + n^(-gamma))
  risk <- A * (2 + d^2 / 4) / n + 2 * sampling_cost * n
  out <- data.frame(term  = c("risk", "n_1", "n_2", "d"),
                    value = c(risk, n, n, d))
  # The decision and the criterion are attributes, so that the value column
  # stays numeric.
  attr(out, "is_satisfied") <- n >= criterion
  attr(out, "criterion") <- criterion
  .as_dmar_tbl(out, subclass = "dmar_mr")
}
