# Asymptotic variance of Cronbach's coefficient alpha.
#' Asymptotic Variance of Coefficient Alpha (Cronbach, Guttman)
#'
#' Computes two large-sample approximations to the variance of the sample
#' coefficient alpha (Guttman, 1945; Cronbach, 1951) under multivariate
#' normality of the \emph{J} item scores and compound symmetry of their
#' covariance matrix (equal item variances and equal covariances, the
#' parallel-items model): the asymptotic variance of van Zyl, Neudecker,
#' and Nel (2000), \eqn{2J(1 - \alpha)^2 / ((J - 1)(N - 1))}, and the
#' variance that follows from Bonett (2002),
#' \eqn{2J(1 - \alpha)^2 / ((J - 1)(N - 2))}.
#'
#' For an interval estimate of alpha from data see
#' \code{\link{reliability_alpha}} and \code{\link{reliability}}; for
#' sample size planning see \code{\link{ss_aipe_reliability}}.
#'
#' @param alpha Population coefficient alpha. Numeric scalar in
#'   \eqn{[0, 1)}.
#' @param N Total sample size (number of \emph{respondents}).
#' @param n_items Number of items contributing to the composite alpha
#'   coefficient (\eqn{J} in the formulas); at least 2.
#'
#' @return A \code{data.frame} with columns \code{term} and \code{value}
#'   and two rows: \code{var_alpha_van_zyl}, the asymptotic variance of
#'   van Zyl, Neudecker, and Nel (2000), and \code{var_alpha_bonett}, the
#'   variance that follows from Bonett (2002).
#'
#' @details
#' \strong{van Zyl, Neudecker, and Nel (2000).} Under multivariate
#' normality and compound symmetry, van Zyl et al. show that
#' \eqn{n\,\mathrm{Var}(\hat\alpha)} tends to
#' \eqn{2J(1 - \alpha)^2 / (J - 1)} as \eqn{n} grows, where
#' \eqn{n = N - 1} is the degrees of freedom of the item covariance
#' matrix (their Equation 13). The asymptotic variance is therefore
#' \deqn{\mathrm{Var}(\hat\alpha) \;\approx\;
#'   \frac{2 J (1 - \alpha)^2}{(J - 1)(N - 1)}.}
#' Without compound symmetry the asymptotic variance is a function of the
#' whole item covariance matrix (their Equation 21);
#' \code{\link{reliability_alpha}} uses that general form for its
#' \code{"ml"} interval, and at a compound symmetric covariance matrix it
#' equals the form above.
#'
#' \strong{Bonett (2002).} Bonett (2003, Equations 2 and 4), citing
#' Bonett (2002), uses \eqn{2J / ((J - 1)(N - 2))} as the large-sample
#' variance of \eqn{\ln(1 - \hat\alpha)} under the same assumptions.
#' Carried to the alpha scale by the delta method it gives
#' \deqn{\mathrm{Var}(\hat\alpha) \;\approx\;
#'   \frac{2 J (1 - \alpha)^2}{(J - 1)(N - 2)}.}
#'
#' \strong{Comparing the two.} The two forms differ only in the
#' denominator: the form with \eqn{N - 2} is the form with \eqn{N - 1}
#' times \eqn{(N - 1)/(N - 2)}, and a sample size chosen to reach a target
#' variance is one larger with the \eqn{N - 2} denominator than with
#' \eqn{N - 1}. Both are large-sample approximations. Under multivariate
#' normality and compound symmetry the exact variance of \eqn{\hat\alpha}
#' for \eqn{N > 5} is
#' \deqn{\frac{2 (1 - \alpha)^2 n (nJ - 2)}{(J - 1)(n - 2)^2 (n - 4)},
#'   \quad n = N - 1}
#' (van Zyl, Neudecker, & Nel, 2000, Equation 12). It is larger than both
#' forms, and the form with \eqn{N - 2} is the closer of the two: with
#' \eqn{J = 5} the exact variance is 2.53 times that form at \eqn{N = 10},
#' 1.27 times at \eqn{N = 30}, and 1.07 times at \eqn{N = 100}.
#'
#' @references
#' Bonett, D. G. (2002). Sample size requirements for testing and
#'   estimating coefficient alpha. \emph{Journal of Educational and
#'   Behavioral Statistics, 27}(4), 335--340. \doi{10.3102/10769986027004335}
#'
#' Bonett, D. G. (2003). Sample size requirements for comparing two alpha
#'   coefficients. \emph{Applied Psychological Measurement, 27}(1),
#'   72--74. \doi{10.1177/0146621602239477}
#'
#' Cronbach, L. J. (1951). Coefficient alpha and the internal structure
#'   of tests. \emph{Psychometrika, 16}(3), 297--334.
#'
#' Guttman, L. (1945). A basis for analyzing test-retest reliability.
#'   \emph{Psychometrika, 10}(4), 255--282.
#'
#' Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
#'   formation for reliability coefficients of homogeneous measurement
#'   instruments. \emph{Methodology, 8}, 39--50.
#'   \doi{10.1027/1614-2241/a000036}
#'
#' Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
#'   population reliability coefficients: Evaluation of methods,
#'   recommendations, and software for composite measures.
#'   \emph{Psychological Methods, 21}, 69--92. \doi{10.1037/a0040086}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' McDonald, R. P. (1999). \emph{Test theory: A unified treatment}.
#'   Lawrence Erlbaum.
#'
#' Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
#'   reliability coefficients: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{British Journal of Mathematical and
#'   Statistical Psychology, 65}, 371--401.
#'   \doi{10.1111/j.2044-8317.2011.02030.x}
#'
#' van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000). On the
#'   distribution of the maximum likelihood estimator of Cronbach's
#'   alpha. \emph{Psychometrika, 65}(3), 271--280. \doi{10.1007/BF02296146}
#'
#' @seealso \code{\link{reliability_alpha}}, \code{\link{reliability}},
#'   \code{\link{ss_aipe_reliability}}
#'
#' @examples
#' # 1. Variance of alpha = 0.80 from a 10-item test with N = 100.
#' var_alpha(alpha = 0.80, N = 100, n_items = 10)
#'
#' # 2. Variance shrinks with N and grows as alpha moves away from 1:
#' var_alpha(alpha = 0.90, N = 50,  n_items = 5)
#' var_alpha(alpha = 0.90, N = 500, n_items = 5)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family variance utilities
#'
#' @export

var_alpha <- function(alpha, N, n_items) {
  # R would match a leftover `n = ` to `n_items` and send the count given
  # by position to `N`, exchanging the two. The name as written is read
  # from the call and refused.
  if ("n" %in% names(sys.call()))
    .dmar_refuse_argument(c(
      "'n' is not an argument of var_alpha(). The sample size is 'N' and ",
      "the number of items is 'n_items'. R would match 'n' to 'n_items' ",
      "by its leading characters and take the other count as the sample ",
      "size, which exchanges the two and gives a wrong variance, so the ",
      "call is refused. Write the sample size as N = , for example ",
      "var_alpha(alpha = .80, N = 100, n_items = 10)."),
      arg = "N")
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      alpha < 0 || alpha >= 1)
    .dmar_refuse_argument(
      "'alpha' must be a single value in [0, 1).", arg = "alpha")
  if (!is.numeric(N) || length(N) != 1L || N < 3)
    .dmar_refuse_argument("'N' must be a single integer >= 3.", arg = "N")
  if (!is.numeric(n_items) || length(n_items) != 1L || n_items < 2)
    .dmar_refuse_argument(
      "'n_items' must be a single integer >= 2.", arg = "n_items")

  # van Zyl, Neudecker, and Nel (2000, Eq. 13): n Var(alpha_hat) tends to
  # 2 J (1 - alpha)^2 / (J - 1) with n = N - 1, the Wishart degrees of
  # freedom. Bonett (2003, Eqs. 2 and 4, citing Bonett, 2002): the variance
  # of ln(1 - alpha_hat) is 2 J / ((J - 1)(N - 2)), carried to the alpha
  # scale by the delta method. DMAR 1.0.x had the two denominators in each
  # other's rows.
  van_zyl <- 2 * n_items * (1 - alpha)^2 / ((n_items - 1) * (N - 1))
  bonett  <- 2 * n_items * (1 - alpha)^2 / ((n_items - 1) * (N - 2))

  out <- data.frame(
    term  = c("var_alpha_van_zyl", "var_alpha_bonett"),
    value = c(van_zyl, bonett),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out)
}
