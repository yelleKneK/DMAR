#' Signal to Noise Estimators for the Squared Multiple Correlation Coefficient
#'
#' Computes five estimators of the population signal to noise ratio
#' \eqn{\phi^2 = \rho^2 / (1 - \rho^2)} associated with the squared multiple
#' correlation coefficient \eqn{\rho^2}. Two are functions of the unadjusted
#' and the Wherry-adjusted sample \eqn{R^2}; the other three are the
#' Muirhead (1985) unique minimum variance unbiased estimators that improve
#' substantially on the plug-in estimator at small \emph{N} and modest
#' numbers of predictors.
#'
#' @param R2 The usual sample estimate of the squared multiple correlation
#'   coefficient (no degrees of freedom adjustment). Numeric scalar in
#'   \eqn{(0, 1)}.
#' @param N Sample size.
#' @param p Number of predictor variables.
#'
#' @details
#' The signal to noise ratio \eqn{\phi^2 = \rho^2 / (1 - \rho^2)} is a
#' natural reparameterization of \eqn{\rho^2} that is bounded only below
#' (at zero) and so behaves more like a variance ratio than a proportion.
#' It is also the noncentrality parameter (up to a factor of \emph{N}) for
#' the omnibus \emph{F}-test of \eqn{\rho^2 = 0} under fixed predictors;
#' see \code{\link{convert_R2_f}}.
#'
#' The five estimators returned, in increasing order of bias-correction
#' machinery, are:
#' \itemize{
#'   \item \code{phi2_hat}: the plug-in estimator
#'     \eqn{\hat\phi^2 = R^2 / (1 - R^2)}. Biased upward in small samples
#'     because the sample \eqn{R^2} is itself biased upward.
#'   \item \code{phi2_adj_hat}: the plug-in estimator applied to the
#'     Wherry-adjusted \eqn{R^2}. Removes the leading-order bias in
#'     \eqn{R^2} but is not itself unbiased for \eqn{\phi^2}.
#'   \item \code{phi2_umvue}: Muirhead's (1985) unique minimum variance
#'     unbiased estimator (his \eqn{\theta_U}, their Eq. 4); equivalent to
#'     Stuart, Ord, and Arnold's (1999) equation 28.97. Requires
#'     \eqn{N \ge p + 6} (the gate on all three Muirhead
#'     estimators); for smaller \eqn{N} the value is \code{NA}.
#'   \item \code{phi2_umvue_l}: Muirhead's (1985) linearly-improved unique
#'     minimum variance unbiased estimator (his \eqn{\theta_L}); equivalent
#'     to Stuart et al.\ (1999) equation 28.98. Dominates \code{phi2_umvue}
#'     in mean squared error.
#'   \item \code{phi2_umvue_nl}: Muirhead's (1985) nonlinearly-improved
#'     estimator (his \eqn{\theta_{NL}}). Dominates the linear improvement
#'     in MSE but requires \eqn{p \ge 5}; for smaller \eqn{p} the value is
#'     \code{NA}.
#' }
#' The nonlinear estimator dominates the linear one in risk when
#' \eqn{p \ge 5}, though Muirhead notes the two perform very similarly in
#' practice; for smaller \eqn{p} the linear estimator is preferred over
#' the plug-in and adjusted-\eqn{R^2} forms. All three Muirhead estimators
#' are reported truncated at zero, so the returned value is on the same
#' scale as \eqn{\phi^2} (the truncation introduces a negligible bias
#' only when the population \eqn{\phi^2} is near zero).
#'
#' As \eqn{N} grows with \eqn{p} fixed, the five estimators converge to a
#' common value (the population \eqn{\phi^2}); the difference between them
#' is the small-sample bias machinery in operation. The \code{@examples}
#' block illustrates that convergence.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value} and
#' one row per estimator: \code{phi2_hat}, \code{phi2_adj_hat},
#' \code{phi2_umvue}, \code{phi2_umvue_l}, and \code{phi2_umvue_nl}.
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
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
#' Muirhead, R. J. (1985). Estimating a particular function of the multiple correlation coefficient.
#' \emph{Journal of the American Statistical Association,  80}, 923--925.
#'
#' Stuart, A., Ord, J. K., & Arnold, S. (1999). \emph{Kendall's advanced
#'   theory of statistics, volume 2A: Classical inference and the linear
#'   model} (6th ed.). Arnold.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_R2}}, \code{\link{ss_aipe_R2}},
#'   \code{\link{convert_R2_f}}
#'
#' @examples
#' # 1. Fixed R^2 = 0.5 and p = 2, growing N: the five estimators agree
#' #    to within a small fraction once N is moderate.
#' signal_to_noise_R2(R2 = .5, N = 50,   p = 2)
#' signal_to_noise_R2(R2 = .5, N = 100,  p = 2)
#' signal_to_noise_R2(R2 = .5, N = 500,  p = 2)
#'
#' # 2. With p = 5 the nonlinear estimator is available; it differs
#' #    most from the plug-in at small N.
#' signal_to_noise_R2(R2 = .5, N = 50,  p = 5)
#' signal_to_noise_R2(R2 = .5, N = 500, p = 5)
#'
#' @keywords models htest
#'
#' @export


signal_to_noise_R2 <- function(R2, N, p) {
  Y <- Y_Adj <- Theta_U <- Theta_L <- Theta_NL <- NA

  # From Muirhead, 1985
  # Uses his notation but define his parameters in terms of N and p.
  n <- N - 1
  m <- p + 1

  Y <- R2 / (1 - R2)

  # Often used method of adjustment.
  R_Square_Adj <- ((N - 1) / (N - p - 1)) * R2 - (p) / (N - p - 1)

  Y_Adj <- R_Square_Adj / (1 - R_Square_Adj)

  if (n - m - 3 >= 1) {
    # Theta_U is Muirhead's (1985) Eq. 4, equal to Stuart, Ord & Arnold's
    # (1999) Equation 28.97. The untruncated value must be kept: Theta_NL
    # below is defined from the untruncated linear estimate (Muirhead's
    # Eq. 10 adds c/Y to theta_L before any truncation), so truncating
    # early would change Theta_NL whenever the sample R2 is small enough
    # to make Theta_U negative.
    Theta_U_raw <- ((n - m - 1) / n) * Y - (m - 1) / n
    Theta_U <- max(Theta_U_raw, 0)

    # Theta_L is Muirhead's Eq. 8, equal to Stuart, Ord & Arnold's (1999)
    # Equation 28.98. The shrinkage factor is positive, so truncating
    # before or after the multiplication reports the same value; the raw
    # version feeds Theta_NL.
    shrink <- (n * (n - m - 3)) / ((n + 2) * (n - m - 1))
    Theta_L_raw <- shrink * Theta_U_raw
    Theta_L <- max(Theta_L_raw, 0)

    # Muirhead's Eq. 10, only for m >= 6 (i.e., p >= 5); otherwise
    # Theta_NL stays NA from its initialization above and
    # 'phi2_umvue_nl' is reported as NA.
    if (m >= 6) {
      Theta_NL <- Theta_L_raw + ((2 * (n - 2) * (n - m - 3) * (m - 5)) / ((n + 2) * (n - m - 1) * (n - m + 1) * (n - m + 3))) * (1 / Y)
      Theta_NL <- max(Theta_NL, 0)
    }
  }
  term <- c('phi2_hat', 'phi2_adj_hat', 'phi2_umvue', 'phi2_umvue_l', 'phi2_umvue_nl')
  value <- c(Y, Y_Adj, Theta_U, Theta_L, Theta_NL)

  return(.as_dmar_tbl(data.frame(term, value)))
}
