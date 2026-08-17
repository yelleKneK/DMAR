#' Minimum Risk Point Estimation of the Population Coefficient of Variation
#'
#' @description
#' A function for the sequential estimation of the coefficient of variations with minimum risk. The function
#' implements the ideas of Chattopadhyay and Kelley (2016), which considers study cost and accuracy of
#' the estimated coefficient of variation simultaneously.
#'
#' @param data the data for which to evaluate the function
#' @param A \eqn{structural_cost/(epsilon^2)}; this is the structural cost that one is willing to pay in a study to estimate the coefficient of variation divided by the square of the desired difference (between the estimate and the parameter)
#' @param structural_cost The structural cost of what one is willing to pay in a study (see note below)
#' @param epsilon The maximum desired difference between the estimated coefficient of variation and the population value)
#' @param sampling_cost The sampling cost to collect an additional observation. For example, if each survey costs 10 dollars to distribute and score, \code{sampling_cost} would be 10 dollars per additional observation
#' @param pilot \code{TRUE} or \code{FALSE} based on whether the users is using the function to plan a pilot sample size (\code{TRUE}) or if it is being used to assess if the optimization criterion has been satisfied (\code{FALSE})
#' @param m0 The minimum bound on the initial pilot sample size
#' @param gamma A correction factor in which we suggest .49; see the two Chattopadhyay & Kelley articles for more details (ignorable for most users)
#' @param verbose If \code{TRUE}, extra information is printed; defaults to \code{FALSE}
#'
#' @details
#' The value of \code{epsilon} is context specific; the smaller the value the closer the estimated value will tend to be to the population value.
#'
#' @return
#' \item{risk}{The value of the risk function }
#' \item{n}{The current sample size}
#' \item{cv}{The current coefficient of variation}
#' \item{is_satisfied}{A TRUE/FALSE statement of whether or not the risk function has been satisfied. If TRUE then sampling can stop as the stopping rule has been satisfied}
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
#' Kelley, K. (2007). Sample size planning for the coefficient of variation
#'   from the accuracy in parameter estimation approach.
#'   \emph{Behavior Research Methods, 39}(4), 755--766.
#'   \doi{10.3758/BF03192966}
#'
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#'   parameter estimation for a general class of effect sizes: A
#'   sequential approach. \emph{Psychological Methods, 23}, 226--243.
#'   \doi{10.1037/met0000127}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' When a study's aim is to estimate a parameter accurately, such as the coefficient of variation, the
#' structural costs and the maximum probable error of the estimate (i.e., \eqn{\epsilon}) are combined
#' to form \eqn{A}. When we say "what the researcher is willing to pay", we literally mean the structural
#' cost (\eqn{c}) the researcher is willing to invest in a study in order to estimate the parameter of
#' interest with the desired degree of accuracy. This value is implicitly included (along with anticipated
#' sampling cost) in grant applications for empirical studies when a certain amount of money is requested
#' to conduct a study. If a researcher is willing to pay more and/or desire a smaller value of
#' \eqn{\epsilon}, \eqn{A} is larger than it would have been. A larger \eqn{A} value will translate into
#' a more expensive study, holding everything else constant. Notice that \eqn{A} is a fixed value in any
#' investigation, as the researcher specifies \eqn{A} directly or by specifying its two components
#' (structural cost and \eqn{\epsilon}) individually. However, what is not fixed but rather evaluated in
#' multiple steps throughout the process is the sampling cost, as it is unknown the necessary sample size
#' in order to accomplish the study's goal of achieving a sufficiently accurate estimate of the coefficient
#' of variation. This is the core of our contributions: minimizing sampling cost, and thereby study cost,
#' by using a sequential procedure that evaluates a stopping rule using the risk function to determine if
#' the optimization criterion has been satisfied (based on the goals of the researcher and current
#' information available). This function implements the ideas of  sampling error and the study costs are
#' considered simultaneously, so that the cost is not higher than necessary for the tolerable sampling error.
#'
#' @seealso
#' \code{\link{ci_cv}}, \code{\link{cv}}, \code{\link{mr_smd}}
#'
#' @examples
#' # Determine pilot sample size:
#' mr_cv(pilot = TRUE, A = 400000, sampling_cost = 75, gamma = .49)
#'
#' # Collect data (the size of which is the pilot sample size)
#' Data <- c(36, 53, 19, 11, 10, 24, 14, 65, 18, 48, 25, 35, 13, 18, 3, 41, 5, 3)
#'
#' # Use mr_cv() to assess if the criterion for stopping the sequential study has been satisfied:
#' mr_cv(data = Data, A = 400000, sampling_cost = 75, gamma = .49)
#'
#' # Collect another data (m=1 here) and perform another check:
#' Data <- c(Data, 44)
#' mr_cv(data = Data, A = 400000, sampling_cost = 75, gamma = .49)
#'
#' # Continue adding obervations, checking each time if m=1, until the minimum risk criteria
#' # are satisfied:
#' Data <- c(Data, 26, 13, 39, 2, 3, 26, 22, 8, 15, 12, 22, 5, 21, 23, 40, 18)
#' mr_cv(data = Data, A = 400000, sampling_cost = 75, gamma = .49)
#'
#' @keywords misc design htest
#'
#' @export


mr_cv <- function(data, A, structural_cost, epsilon, sampling_cost, pilot = FALSE, m0 = 4, gamma = .49, verbose = FALSE) {

  # Internal functions (could be stand alone if one desired by simply
  # putting them in the workspace directly as their own functions).
  ########################################################################################################################
  V2_n <- function(data) {
    if (is.data.frame(data)) {
      data <- as.vector(data)
    }
    n <- length(data)
    t1 <- (sqrt(var(data)))^4 / (4 * mean(data)^4)
    t2 <- mu_4_n(data) / (4 * mean(data)^4)
    t3 <- var(data) / (2 * mean(data)^2)
    t4 <- mu_3_n(data) / (mean(data)^3)
    return(t1 + t2 + t3 - t4)
  }


  mu_4_n <- function(data) {
    if (is.data.frame(data)) {
      data <- as.vector(data)
    }
    n <- length(data)
    y2 <- data^2
    y3 <- data^3
    y4 <- data^4
    w1 <- n^2 * sum((data - mean(data))^4)
    w2 <- (-2 * n + 3) * sum(y4)
    w3 <- (8 * n - 12) * mean(data) * sum(y3)
    w4 <- (-6 + 9 / n) * (sum(y2))^2
    return((w1 + w2 + w3 + w4) / ((n - 1) * (n - 2) * (n - 3)))
  }

  mu_3_n <- function(data) {
    if (is.data.frame(data)) {
      data <- as.vector(data)
    }
    n <- length(data)
    w <- sum((data - mean(data))^3)
    return(n * w / ((n - 1) * (n - 2)))
  }
  ########################################################################################################################

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
    if (is.data.frame(data)) {
      data <- as.vector(data)
    }
    n <- length(data)

    Criterion <- (A / sampling_cost) * (V2_n(data) + n^(-2 * gamma))

    CV <- cv(mean = mean(data), sd = sqrt(var(data)))[1, 2]

    Rk <- (1 / n) * A * V2_n(data) + (sampling_cost * n) # This is the risk function.

    is_sat <- (n^2 >= Criterion)
    if (verbose == FALSE) {
      Outcome <- data.frame(
        term  = c("risk", "n", "cv"),
        value = c(Rk, n, CV)
      )
    }
    if (verbose == TRUE) {
      Outcome <- data.frame(
        term  = c("risk", "n", "cv", "v2_n", "mu_3_n", "mu_4_n",
                  "criterion"),
        value = c(Rk, n, CV, V2_n(data), mu_3_n(data), mu_4_n(data),
                  Criterion)
      )
    }
    # is_satisfied (logical) and the criterion live on attributes so
    # that the value column stays numeric.
    attr(Outcome, "is_satisfied") <- is_sat
    attr(Outcome, "criterion") <- Criterion
  }

  if (pilot == TRUE) {
    Outcome <- data.frame(term = "pilot_ss", value = max(m0, ceiling((A / sampling_cost)^(1 / ((2 + 2 * gamma))))))
  }

  return(.as_dmar_tbl(Outcome))
}
