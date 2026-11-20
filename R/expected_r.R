# Exact expected value of the sample Pearson r given rho and N.
#' Exact Expected Value of the Sample Pearson Correlation Given \eqn{\rho} and \eqn{N}
#'
#' Computes \eqn{\mathrm{E}[r \mid \rho, N]}, the exact expected value of the
#' sample Pearson product-moment correlation coefficient under bivariate
#' normality, using the Olkin-Pratt (1958) closed form built on Hotelling's
#' (1953) exact density. The sample \eqn{r} is downward-biased as an
#' estimator of the population \eqn{\rho}, a fact known since Soper (1913)
#' and Fisher (1915). The exact expected value requires summing a Gauss
#' hypergeometric series (see Details), which this function does, so that
#' users can see the exact bias and decide whether to apply the
#' Olkin-Pratt (1958) unbiased estimator of \eqn{\rho} (see Details).
#'
#' \code{expected_r()} is especially useful at the \emph{design stage} of a
#' study, where sample size planning typically proceeds from an assumed
#' value of the population correlation \eqn{\rho}. Because the value of
#' \eqn{r} a researcher should expect to observe on average is
#' \emph{smaller in absolute value} than the assumed \eqn{\rho}, plugging
#' \eqn{\rho} directly into sampling-distribution machinery (precision of
#' \eqn{r}, width of a confidence interval, power of a test of
#' \eqn{H_0\!: \rho = 0}) systematically over-promises on the realized
#' precision or power. Substituting \code{expected_r(true_cor, N)} for the bare
#' \eqn{\rho} in such planning calculations corrects the leading-order
#' over-promise. See \strong{Examples} for a design-stage walk-through.
#'
#' @param true_cor Population correlation coefficient \eqn{\rho}. A numeric scalar or
#'   vector in the interval \eqn{[-1, 1]}.
#' @param N Sample size on which the Pearson \eqn{r} would be computed.
#'   A scalar or vector of integers with \eqn{N \ge 4}. (At \eqn{N = 3}
#'   the bias is defined but the exact formula is numerically delicate;
#'   we require \eqn{N \ge 4} for well-conditioned computation.) When
#'   \code{true_cor} and \code{N} are both vectors they must have the same length or recycle cleanly.
#'
#' @return A \code{data.frame} with one row per (\code{true_cor}, \code{N})
#'   input and the columns
#'   \itemize{
#'     \item \code{true_cor}, the input population correlation,
#'     \item \code{N}, the input sample size,
#'     \item \code{expected_r}, \eqn{\mathrm{E}[r \mid \rho, N]},
#'     \item \code{bias}, \eqn{\rho - \mathrm{E}[r \mid \rho, N]} (the
#'       amount by which the sample \eqn{r} underestimates \eqn{\rho} on
#'       average),
#'     \item \code{relative_bias}, \code{bias / true_cor} when \eqn{\rho \ne 0},
#'       and \code{NA} when \eqn{\rho = 0}.
#'   }
#'
#' @details
#' \strong{The exact formula.} Under bivariate normality, if
#' \eqn{r} is the sample Pearson correlation in a sample of size \eqn{N},
#' \deqn{\mathrm{E}[r \mid \rho, N] \;=\;
#'   \rho \,\cdot\,
#'   {}_2F_1\!\left(\tfrac{1}{2},\, \tfrac{1}{2};\, \tfrac{N+1}{2};\, \rho^2 \right)
#'   \,\cdot\,
#'   \frac{\Gamma\!\left(\tfrac{N}{2}\right)^2}
#'        {\Gamma\!\left(\tfrac{N-1}{2}\right)\,\Gamma\!\left(\tfrac{N+1}{2}\right)}.}
#' Here \eqn{{}_2F_1(a, b; c; z) = \sum_{k=0}^{\infty} \frac{(a)_k (b)_k}{(c)_k\, k!}\, z^k}
#' is the Gauss hypergeometric function with Pochhammer symbols
#' \eqn{(x)_k = x(x+1)\cdots(x+k-1)} (Hotelling, 1953, Section 7 for the
#' moments of \eqn{r}, Section 3, equation 25, for the exact density;
#' Olkin & Pratt, 1958, equation 3.2, for this closed form).
#' For \eqn{\rho^2 < 1} the series converges absolutely; this implementation
#' sums by a numerically stable forward recurrence and stops when the
#' relative contribution of the next term falls below a tolerance.
#' The \code{bias} column is computed directly, as \eqn{\rho} times one
#' minus the product of the series and the gamma ratio, formed from their
#' logarithms, and not as the difference \eqn{\rho - \mathrm{E}[r]}. At a
#' large \eqn{N} the two numbers share their leading digits (about eight of
#' them at \eqn{N = 10^8}), and the subtraction would leave the bias only
#' the digits after those.
#'
#' \strong{Tuning the series convergence (rarely needed).} Two internal
#' tuning constants control the \eqn{{}_2F_1} series summation:
#' \itemize{
#'   \item \code{dmar.expected_r.tol}, relative tolerance for stopping
#'     the series (default \code{1e-15}). The summation stops once a term
#'     is below this fraction of the sum of the terms after the leading 1,
#'     the only part of the series on which the bias depends.
#'   \item \code{dmar.expected_r.max_iter}, maximum number of series
#'     terms (default \code{5000L}). If the last term summed is still not
#'     below \code{dmar.expected_r.tol} times the whole series, a
#'     non-convergence warning is issued.
#' }
#' Users who need different values
#' (e.g., for very near-boundary \eqn{|\rho|} where the series converges
#' slowly) can set them via \code{\link[base]{options}()}, for example
#' \code{options(dmar.expected_r.tol = 1e-12)}. They are deliberately
#' hidden from the function signature so as not to clutter the everyday
#' user's view of the call.
#'
#' \strong{Why the bias is present even though \eqn{s^2} is unbiased for
#' \eqn{\sigma^2}.} The downward bias arises because \eqn{r} is a nonlinear
#' function of unbiased sample moments. By Jensen's inequality and the
#' concavity of the square root in the denominator of \eqn{r}, the expectation
#' of the ratio is not the ratio of the expectations. The bias is largest
#' when \eqn{N} is small or \eqn{|\rho|} is moderate; as \eqn{N \to \infty},
#' \eqn{\mathrm{E}[r] \to \rho}.
#'
#' \strong{Sign and magnitude.} The bias \eqn{\rho - \mathrm{E}[r]} has the
#' same sign as \eqn{\rho} and is approximately
#' \eqn{\rho(1 - \rho^2)/[2(N - 1)]} to leading order (Fisher, 1915;
#' Hotelling, 1953; Ghosh, 1966); the exact formula above incorporates all
#' higher-order corrections. For \eqn{\rho = 0.5}, \eqn{N = 10}, the bias
#' is about \eqn{+0.021}; for \eqn{\rho = 0.5}, \eqn{N = 30}, about
#' \eqn{+0.0065}.
#'
#' \strong{Olkin-Pratt (1958) unbiased estimator of \eqn{\rho}.} The
#' companion to this expected-value calculation is the Olkin-Pratt
#' unbiased estimator of \eqn{\rho} given an observed \eqn{r}:
#' \deqn{\tilde{\rho}_{\text{OP}}(r, N) \;=\;
#'   r \,\cdot\, {}_2F_1\!\left(\tfrac{1}{2},\, \tfrac{1}{2};\, \tfrac{N-2}{2};\, 1 - r^2 \right).}
#' Olkin & Pratt (1958) prove
#' \eqn{\mathrm{E}[\tilde{\rho}_{\text{OP}}(r, N) \mid \rho, N] = \rho}
#' exactly, for every \eqn{\rho \in (-1, 1)} and \eqn{N \ge 4}. The
#' commonly quoted first-order approximation
#' \eqn{\tilde{\rho} \approx r\,[1 + (1 - r^2)/(2(N - 3))]} (Olkin, 1967)
#' is the truncation of the OP series at \eqn{k = 1}; the full series is
#' what makes the estimator exactly unbiased.
#'
#' \strong{Size of the bias relative to sampling variability.} The
#' standard deviation of \eqn{r} is approximately
#' \eqn{(1 - \rho^2)/\sqrt{N - 1}}, so the bias of a single \eqn{r} is
#' about \eqn{\rho/(2\sqrt{N - 1})} of its standard deviation, whether
#' the coefficient is tested, given a confidence interval, or used as a
#' reliability, validity, or planning value; at \eqn{\rho = 0.8} and
#' \eqn{N = 20} the exact ratio is 0.09. When \eqn{K} correlations from
#' samples of the same \eqn{N} and the same \eqn{\rho} are averaged, as in
#' a meta-analysis, the bias of the average equals the bias of a single
#' \eqn{r}, while the standard error of the average falls as
#' \eqn{1/\sqrt{K}}; the bias is then about
#' \eqn{\rho\sqrt{K}/(2\sqrt{N - 1})} standard errors of the average, and
#' at \eqn{\rho = 0.8}, \eqn{N = 20}, and \eqn{K = 50} the exact ratio is
#' 0.63.
#'
#' \strong{Connection to Fisher's Z transform.} The variance-stabilizing
#' transform \eqn{z = \tanh^{-1}(r)}, proposed in passing in Fisher (1915,
#' p. 521) and developed in Fisher (1921), has approximate variance
#' \eqn{1/(N-3)} regardless of \eqn{\rho}, but \eqn{\mathrm{E}[z]} also
#' carries a small-sample bias of order \eqn{1/N} (Hotelling, 1953,
#' Section 8). Hotelling (1953, Sections 9--10) gives bias-adjusted and
#' variance-stabilized refinements of \eqn{Z}, e.g.
#' \eqn{z - (3z + r)/(4n)}.
#'
#' @references
#' Anderson, T. W. (2003). \emph{An introduction to multivariate
#'   statistical analysis} (3rd ed.), Section 4.2. Wiley.
#'
#' Fisher, R. A. (1915). Frequency distribution of the values of the
#'   correlation coefficient in samples from an indefinitely large
#'   population. \emph{Biometrika, 10}(4), 507--521.
#'
#' Fisher, R. A. (1921). On the "probable error" of a coefficient of
#'   correlation deduced from a small sample. \emph{Metron, 1}, 3--32.
#'
#' Ghosh, B. K. (1966). Asymptotic expansions for the moments of the
#'   distribution of correlation coefficient. \emph{Biometrika, 53}(1/2),
#'   258--262.
#'
#' Hotelling, H. (1953). New light on the correlation coefficient and its
#'   transforms. \emph{Journal of the Royal Statistical Society, Series B,
#'   15}(2), 193--232. (Discussion, pp.\ 225--232.)
#'
#' Olkin, I. (1967). Correlations revisited. In J. C. Stanley (Ed.),
#'   \emph{Improving experimental design and statistical analysis}
#'   (pp.\ 102--128). Rand McNally.
#'
#' Olkin, I., & Pratt, J. W. (1958). Unbiased estimation of certain
#'   correlation coefficients. \emph{The Annals of Mathematical Statistics,
#'   29}(1), 201--211.
#'
#' Soper, H. E. (1913). On the probable error of the correlation
#'   coefficient to a second approximation. \emph{Biometrika, 9}(1/2),
#'   91--115.
#'
#' Soper, H. E., Young, A. W., Cave, B. M., Lee, A., & Pearson, K. (1917).
#'   On the distribution of the correlation coefficient in small samples.
#'   Appendix II to the papers of "Student" and R. A. Fisher: A cooperative
#'   study. \emph{Biometrika, 11}(4), 328--413.
#'
#' Stuart, A., & Ord, J. K. (1994). \emph{Kendall's advanced theory of
#'   statistics, Vol.\ 1: Distribution theory} (6th ed.), Section 16.32.
#'   Edward Arnold.
#'
#' @seealso \code{\link{ci_r}}, \code{\link[stats]{cor}},
#'   \code{\link[stats]{cor.test}}
#'
#' @examples
#' # 1. A single value: true_cor = 0.5, N = 10. The sample r is downwardly
#' #        biased by about 0.021, roughly 4% of true_cor.
#' expected_r(true_cor = 0.5, N = 10)
#'
#' # 2. Bias as a function of N for fixed true_cor. The bias is roughly
#' #        true_cor * (1 - true_cor^2) / (2(N - 1)) to leading order; as N
#' #        grows it shrinks toward zero.
#' expected_r(true_cor = 0.5, N = c(5, 10, 20, 50, 100, 500))
#'
#' # 3. Bias as a function of true_cor for fixed N. The bias is zero at
#' #        true_cor = 0 and true_cor = +/- 1, and largest near
#' #        true_cor = +/- 0.6.
#' expected_r(true_cor = seq(0, 0.95, by = 0.05), N = 10)
#'
#' # 4. The Olkin-Pratt unbiased estimator: invert the bias for an
#' #        observed sample r.
#' set.seed(113)
#' x <- rnorm(20); y <- 0.4 * x + sqrt(1 - 0.4^2) * rnorm(20)
#' r_obs <- cor(x, y)
#' r_obs
#'
#' # Olkin-Pratt unbiased estimator of rho:
#' op_unbiased <- function(r, N, tol = 1e-15, max_iter = 5000) {
#'   z <- 1 - r^2; c_par <- (N - 2) / 2
#'   s <- 1; term <- 1
#'   for (k in seq_len(max_iter)) {
#'     term <- term * ((k - 0.5)^2) / ((c_par + k - 1) * k) * z
#'     s <- s + term
#'     if (abs(term) < tol * abs(s)) break
#'   }
#'   r * s
#' }
#' op_unbiased(r_obs, N = 20)
#'
#' # 5. Design-stage use: a study planned around true_cor = 0.4 with N = 30.
#' #        Naive plug-in says we expect to observe r = 0.40 on average,
#' #        but the realized expected r is smaller, and that gap matters
#' #        for any precision- or power-based sample size calculation that
#' #        plugs in true_cor as if it were the expected sample r.
#' cor_planned <- 0.4
#' n_planned   <- 30
#' expected_r(true_cor = cor_planned, N = n_planned)
#'
#' # 6. Tuning constants are hidden from the signature but tunable
#' #        through options() for the rare cases that need them (e.g.,
#' #        very near-boundary |true_cor| where the 2F1 series converges
#' #        slowly). The defaults rarely need to be changed.
#' options(dmar.expected_r.tol = 1e-12,
#'         dmar.expected_r.max_iter = 20000L)
#' expected_r(true_cor = 0.999, N = 5)
#' options(dmar.expected_r.tol = NULL,
#'         dmar.expected_r.max_iter = NULL)   # restore defaults
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest correlation
#'
#' @family effect size estimates
#'
#' @export

expected_r <- function(true_cor, N) {
  if (!is.numeric(true_cor)) .dmar_refuse_argument(
    "'true_cor' must be numeric.", arg = "true_cor")
  if (!is.numeric(N)) .dmar_refuse_argument("'N' must be numeric.", arg = "N")
  if (any(abs(true_cor) > 1, na.rm = TRUE))
    .dmar_refuse_argument("'true_cor' must lie in [-1, 1].", arg = "true_cor")
  if (any(N < 4, na.rm = TRUE))
    .dmar_refuse_argument(
      "'N' must be >= 4 for well-conditioned computation of E[r].", arg = "N")

  # Hidden tuning constants for the 2F1 series; documented in @details.
  # Pulled from options() so they are controllable but not visible in the
  # signature; sensible defaults that almost never need to be changed.
  tol      <- getOption("dmar.expected_r.tol",      1e-15)
  max_iter <- as.integer(getOption("dmar.expected_r.max_iter", 5000L))

  # .expected_r_one() returns the relative bias w = 1 - E[r] / rho, from
  # which E[r] = rho (1 - w) and the bias rho - E[r] = rho w, so the bias
  # keeps the relative precision of w at any N.
  args <- data.frame(true_cor = true_cor, N = N)
  w <- vapply(seq_len(nrow(args)), function(i) {
    .expected_r_one(args$true_cor[i], args$N[i], tol = tol, max_iter = max_iter)
  }, numeric(1))

  .as_dmar_tbl(data.frame(
    true_cor      = args$true_cor,
    N             = args$N,
    expected_r    = args$true_cor * (1 - w),
    bias          = args$true_cor * w,
    relative_bias = ifelse(args$true_cor == 0, NA_real_, w),
    stringsAsFactors = FALSE,
    row.names     = NULL
  ))
}


# Single-(rho, N) worker for the exact Olkin-Pratt (Hotelling, 1953)
# expectation E[r | rho, N] = rho * 2F1(1/2, 1/2; (N+1)/2; rho^2) *
# Gamma(N/2)^2 / (Gamma((N-1)/2) * Gamma((N+1)/2)) under bivariate normality.
# It returns the relative bias w = 1 - E[r] / rho = 1 - 2F1 * gamma ratio,
# which does not depend on the sign of rho. At rho = 0 it returns the limit
# of w as rho goes to 0, which expected_r() multiplies by rho = 0.
.expected_r_one <- function(true_cor, N, tol, max_iter) {
  if (is.na(true_cor) || is.na(N)) return(NA_real_)
  if (abs(true_cor) == 1) return(0)  # boundary: r is degenerate at +/-1

  z <- true_cor^2

  # 2F1(1/2, 1/2; (N+1)/2; rho^2) = 1 + t by forward recurrence:
  # the standard 2F1 ratio of consecutive terms is
  #   term_{k+1} / term_k = (a + k)(b + k) / ((c + k)(k + 1)) * z,
  # which for a = b = 1/2, c = (N+1)/2 simplifies to
  #   (k - 1/2)^2 / (((N+1)/2 + k - 1) * k) * z
  # when iterating from term_{k-1} to term_k. The terms after the leading 1
  # are summed apart from it, as t: at a large N, t is about rho^2 / (2N),
  # and adding it to 1 would keep only its leading digits. The terms are
  # positive and decreasing (zero when z is 0). The sum stops once a term
  # is below tol times t, the precision the bias needs, or after max_iter
  # terms. The warning is given, as in 1.0.x, only when the last term
  # is still not below tol times the whole sum 1 + t.
  c_par <- (N + 1) / 2
  t <- 0
  term <- 1
  for (k in seq_len(max_iter)) {
    term <- term * ((k - 0.5)^2) / ((c_par + k - 1) * k) * z
    t <- t + term
    if (term <= tol * t) break
  }
  if (term >= tol * (1 + t)) {
    .dmar_warn(
      sprintf(
"expected_r(): 2F1 series did not converge to tol = %.1e in %d iterations (true_cor = %.3f, N = %d). Returning current partial sum.",
tol, max_iter, true_cor, N
), "result_qualified")
  }

  # Gamma ratio: Gamma(N/2)^2 / (Gamma((N-1)/2) * Gamma((N+1)/2)). Since
  # Gamma((N+1)/2) = ((N-1)/2) Gamma((N-1)/2), it equals
  # 2/(N-1) {Gamma(N/2) / Gamma((N-1)/2)}^2 = c4(N)^2, the square of
  # E(S) / sigma for a normal sample of size N. Its logarithm comes from
  # .chi_log_mean_factor() (R/c4_internals.R) on N - 1 df, and
  #   w = 1 - (1 + t) c4(N)^2 = -expm1(log1p(t) + 2 log c4(N)).
  # The only cancellation left is between the two logarithms, about
  # rho^2 / (2N) and -1 / (2N) at a large N, and it is severe only as rho^2
  # nears 1, where the bias itself goes to 0. Through 1.0.x the ratio was
  # the sum of three lgamma() values, each about (N/2) log(N/2), which lost
  # the digits they share (it was off by 1e-7 at N = 1e8, so that E[r]
  # could come out above rho), and the bias was rho - E[r], which keeps
  # only the digits the two do not share: with the ratio right to a few
  # units in the last place, that difference is off by 3e-7 at N = 1e8
  # and 7e-6 at N = 1e9, relative, at rho = .5.
  -expm1(log1p(t) + 2 * .chi_log_mean_factor(N - 1))
}
