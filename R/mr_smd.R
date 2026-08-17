#' Minimum Risk Point Estimation of the Population Standardized Mean Difference
#'
#' @description
#' A function for the sequential estimation of the standardized mean difference with minimum risk. The
#' function implements the ideas of Chattopadhyay and Kelley (2017), which
#' considers study cost and accuracy of the estimated standardized mean difference simultaneously. This
#' is important to specify that \code{mr_smd} was developed under the assumption of normally distributed
#' data with equal sample size and equal cost of sampling per observation for each group.
#'
#' @param A The price one is willing to pay in order to have a maximum allowable difference of \eqn{epsilon^2} between the estimate of the standardized mean difference and its corresponding parameter
#' @param structural_cost The structural cost of what one is willing to pay in a study
#' @param epsilon The maximum desired difference between the estimated standardized mean difference and the population value
#' @param d The current estimate of the standardized mean difference
#' @param n Current sample size \emph{per group} (thus total sample size is \eqn{2n}); requires equal sample size \emph{per group}
#' @param sampling_cost The sampling cost to collect an additional observation. For example, if each survey costs 10 dollars to distribute and score, \code{sampling_cost} would be 10 dollars per additional observation
#' @param pilot \code{TRUE} or \code{FALSE} based on whether the users is using the function to plan a pilot sample size (TRUE) or if it is being used to assess if the optimization criterion has been satisfied (FALSE)
#' @param m0 The minimum bound on the initial pilot sample size
#' @param gamma A correction factor in which we suggest .49; see the two Chattopadhyay & Kelley articles for more details (ignorable for most users)
#'
#' @details
#' The standardized mean difference is a widely used measure effect size. In this article, we developed a
#' general theory for estimating the population standardized mean difference by minimizing both the mean
#' square error of the estimator and the total sampling cost. This function implements our ideas discussed
#' in Chattopadhyay and Kelley (2017). See also Kelley and Rausch (2006) for additional information on
#' the standardized mean difference.
#'
#' @return
#' \item{risk}{The value of the risk function.}
#' \item{n1}{Sample size for group 1 (echos the input value)}
#' \item{n2}{Sample size for group 2 (echos the input value)}
#' \item{d}{Observed value of the standardized mean difference (i.e., \emph{d}; echos the input value)}
#' \item{is_satisfied}{A \code{TRUE} or \code{FALSE} statement of that evaluates a stopping rule using the risk function to determine if the optimization criterion has been satisfied (based on the goals of the researcher and current information available)}
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
#'   sequential approach. \emph{Psychological Methods, 23}, 226--243.
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
#' When \code{pilot=TRUE} the function returns the size of the pilot sample size, \emph{per group},
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
#' mr_smd(d = 1.021484, n = 13, A = 10000, sampling_cost = 2.40, gamma = .49)
#'
#' # Or, using the smd() function:
#' mr_smd(d = smd(group_1 = High.SLS, group_2 = Low.SLS)[1,2], n = 13, A = 10000,
#'        sampling_cost = 2.40, gamma = .49)
#'
#' # Here, for this situation, the stopping rule is satisfied:
#' mr_smd(d = 1.00, n = 75, A = 10000, sampling_cost = 2.40, gamma = .49)
#'
#' @keywords design misc htest
#'
#' @export


mr_smd <- function(A, structural_cost, epsilon, d, n, sampling_cost, pilot = FALSE, m0 = 4, gamma = .49) {
  # d current estimate of the standardized mean difference.
  # n current sample size per group.
  # A the cost that the researcher would be willing to pay to have a maximum allowable difference of the population standardized mean difference and its estimate of $\epsilon$.
  # sampling_cost the cost of collecting an observation; thus the cost is 2n for increasing sample size by 1 in each group.
  # gamma
  # m pilot sample size.
  # Here we check to see if the current sample size is sufficiently large

  if (!missing(A)) {
    if (!missing(structural_cost)) stop("Because you specified \'A\' directly, you should not also specify \'structural_cost\'.")
    if (!missing(epsilon)) stop("Because you specified \'A\' directly, you should not also specify \'epsilon\'.")
    if (A <= 0) stop("A should be a non-zero positive value")
  }

  if (missing(A)) {
    if (missing(structural_cost)) stop("Because you did not specificy \'A\' directly, you must specify \'structural_cost\'.")
    if (missing(epsilon)) stop("Because you did not specificy \'A\' directly, you must specify \'epsilon\'.")
    A <- structural_cost / (epsilon^2)
    if (A <= 0) stop("A should be a non-zero positive value")
  }

  if (pilot == FALSE) {
    Stop <- FALSE

    Criterion <- sqrt(A / (2 * sampling_cost)) * (sqrt(2 + d^2 / 4) + n^(-gamma))

    Rk <- A * ((1 / n + 1 / n) + (d^2) / (2 * (n + n))) + sampling_cost * (n + n)
    Outcome <- data.frame(
      term  = c("risk", "n1", "n2", "d"),
      value = c(Rk, n, n, d)
    )
    # is_satisfied (logical) lives on an attribute so that the value
    # column stays numeric.
    attr(Outcome, "is_satisfied") <- (n >= Criterion)
    attr(Outcome, "criterion") <- Criterion
  }


  if (pilot == TRUE) {
    if (m0 < 4) stop("The value of 'm0' must be 4 or greater.")
    Outcome <- data.frame(term = "pilot_ss", value = max(m0, ceiling((A / (2 * sampling_cost))^(1 / (2 + 2 * gamma)))))
  }

  return(.as_dmar_tbl(Outcome))
}
