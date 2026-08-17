#' Variance of the Squared Multiple Correlation Coefficient
#'
#' @description
#' Computes the sampling variance of the squared multiple correlation
#' coefficient from the population value, the sample size, and the
#' number of predictors, the quantity that governs how precisely
#' \eqn{R^2} is estimated at a given design size.
#'
#' @param population_R2 Population squared multiple correlation coefficient
#' @param N Sample size
#' @param p The number of predictor variables
#'
#' @details
#' Uses the hypergeometric function as discussed in and section 28 of Stuart, Ord, and Arnold (1999) in
#' order to obtain the \emph{correct} value for the variance of the squared multiple correlation coefficient.
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} and \code{value}.
#' The \code{term} value is \code{"var_R2"} and \code{value} is the
#' asymptotic variance of \eqn{R^2}.
#'
#' @references
#' Kelley, K. (2008). Sample size planning for the squared multiple
#'   correlation coefficient: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{Multivariate Behavioral Research, 43},
#'   524--555. \doi{10.1080/00273170802490632}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{R^2} as a model comparison
#'   effect size.)
#'
#' Stuart, A., Ord, J. K., & Arnold, S. (1999). \emph{Kendall's advanced
#'   theory of statistics, volume 2A: Classical inference and the linear
#'   model} (6th ed.). Arnold.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' The Gauss hypergeometric function \eqn{{}_2F_1} is computed in base R (see
#' the internal \code{.hyperg_2F1}); no GSL system library is required.
#'
#' @seealso
#' \code{\link{expected_R2}}, \code{\link{ci_R2}}, \code{\link{ss_aipe_R2}}
#'
#' @examples
#' var_R2(.5, 10, 5)
#' var_R2(.5, 25, 5)
#' var_R2(.5, 50, 5)
#' var_R2(.5, 100, 5)
#'
#' @keywords design
#'
#' @export

var_R2 <- function(population_R2, N, p) {
  # gsl::hyperg_2F1(2, 2, .5 * (N + 3), population_R2) is another way to obtain
  # this value; DMAR now computes the 2F1 in base R via .hyperg_2F1() (no GSL
  # system dependency). See R/R2_internals.R.
  result <- (((N - p - 1) * (N - p + 1)) / (N^2 - 1)) * ((1 - population_R2)^2) * (.hyperg_2F1(2, 2, .5 * (N + 3), population_R2)) - ((expected_R2(population_R2 = population_R2, N = N, p = p)[1,2] - 1)^2)
  return(.as_dmar_tbl(data.frame(term = 'var_R2', value = result)))
}
