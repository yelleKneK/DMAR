# Exact expected value of the sample Pearson r given rho and n.
#' Exact Expected Value of the Sample Pearson Correlation Given \eqn{\rho} and \eqn{n}
#'
#' Computes \eqn{\mathrm{E}[r \mid \rho, n]}, the exact expected value of the
#' sample Pearson product-moment correlation coefficient under bivariate
#' normality, using the Olkin-Pratt (1958) closed form built on Hotelling's
#' (1953) exact density. The sample \eqn{r} is downward-biased as an
#' estimator of the population \eqn{\rho}, a fact known since Soper (1913)
#' and Fisher (1915); although this exact-bias equation has been available
#' for more than half a century, it has historically not been easy to
#' implement in general-purpose statistical software and is correspondingly
#' rarely used in applied work. This function makes the exact bias visible
#' so that users can decide whether to apply the Olkin-Pratt (1958)
#' unbiased estimator of \eqn{\rho} (see Details).
#'
#' \code{expected_r()} is especially useful at the \emph{design stage} of a
#' study, where sample size planning typically proceeds from an assumed
#' value of the population correlation \eqn{\rho}. Because the value of
#' \eqn{r} a researcher should expect to observe on average is
#' \emph{smaller in absolute value} than the assumed \eqn{\rho}, plugging
#' \eqn{\rho} directly into sampling-distribution machinery (precision of
#' \eqn{r}, width of a confidence interval, power of a test of
#' \eqn{H_0\!: \rho = 0}) systematically over-promises on the realized
#' precision or power. Substituting \code{expected_r(rho, n)} for the bare
#' \eqn{\rho} in such planning calculations corrects the leading-order
#' over-promise. See \strong{Examples} for a design-stage walk-through.
#'
#' @param rho Population correlation coefficient. A numeric scalar or
#'   vector in the interval \eqn{[-1, 1]}.
#' @param n Sample size on which the Pearson \eqn{r} would be computed.
#'   A scalar or vector of integers with \eqn{n \ge 4}. (At \eqn{n = 3}
#'   the bias is defined but the exact formula is numerically delicate;
#'   we require \eqn{n \ge 4} for well-conditioned computation.) When \code{rho} and \code{n} are both
#'   vectors they must have the same length or recycle cleanly.
#'
#' @return A \code{data.frame} with one row per (\code{rho}, \code{n})
#'   input and the columns
#'   \itemize{
#'     \item \code{rho}, the input population correlation,
#'     \item \code{n}, the input sample size,
#'     \item \code{expected_r}, \eqn{\mathrm{E}[r \mid \rho, n]},
#'     \item \code{bias}, \eqn{\rho - \mathrm{E}[r \mid \rho, n]} (the
#'       amount by which the sample \eqn{r} underestimates \eqn{\rho} on
#'       average),
#'     \item \code{relative_bias}, \code{bias / rho} when \eqn{\rho \ne 0},
#'       and \code{NA} when \eqn{\rho = 0}.
#'   }
#'
#' @details
#' \strong{The exact formula.} Under bivariate normality, if
#' \eqn{r} is the sample Pearson correlation in a sample of size \eqn{n},
#' \deqn{\mathrm{E}[r \mid \rho, n] \;=\;
#'   \rho \,\cdot\,
#'   {}_2F_1\!\left(\tfrac{1}{2},\, \tfrac{1}{2};\, \tfrac{n+1}{2};\, \rho^2 \right)
#'   \,\cdot\,
#'   \frac{\Gamma\!\left(\tfrac{n}{2}\right)^2}
#'        {\Gamma\!\left(\tfrac{n-1}{2}\right)\,\Gamma\!\left(\tfrac{n+1}{2}\right)}.}
#' Here \eqn{{}_2F_1(a, b; c; z) = \sum_{k=0}^{\infty} \frac{(a)_k (b)_k}{(c)_k\, k!}\, z^k}
#' is the Gauss hypergeometric function with Pochhammer symbols
#' \eqn{(x)_k = x(x+1)\cdots(x+k-1)} (Hotelling, 1953, Section 7 for the
#' moments of \eqn{r}, Section 3, equation 25, for the exact density;
#' Olkin & Pratt, 1958, equation 3.2, for this closed form).
#' For \eqn{\rho^2 < 1} the series converges absolutely; this implementation
#' sums by a numerically stable forward recurrence and stops when the
#' relative contribution of the next term falls below a tolerance.
#'
#' \strong{Tuning the series convergence (rarely needed).} Two internal
#' tuning constants control the \eqn{{}_2F_1} series summation:
#' \itemize{
#'   \item \code{DMAR.expected_r.tol}, relative tolerance for stopping
#'     the series (default \code{1e-15}).
#'   \item \code{DMAR.expected_r.max_iter}, maximum number of series
#'     terms before issuing a non-convergence warning (default
#'     \code{5000L}).
#' }
#' Both have sensible defaults; advanced users who need different values
#' (e.g., for very near-boundary \eqn{|\rho|} where the series converges
#' slowly) can set them via \code{\link[base]{options}()}, for example
#' \code{options(DMAR.expected_r.tol = 1e-12)}. They are deliberately
#' hidden from the function signature so as not to clutter the everyday
#' user's view of the call.
#'
#' \strong{Why the bias is present even though \eqn{s^2} is unbiased for
#' \eqn{\sigma^2}.} The downward bias arises because \eqn{r} is a nonlinear
#' function of unbiased sample moments. By Jensen's inequality and the
#' concavity of the square root in the denominator of \eqn{r}, the expectation
#' of the ratio is not the ratio of the expectations. The bias is largest
#' when \eqn{n} is small or \eqn{|\rho|} is moderate; as \eqn{n \to \infty},
#' \eqn{\mathrm{E}[r] \to \rho}.
#'
#' \strong{Sign and magnitude.} The bias \eqn{\rho - \mathrm{E}[r]} has the
#' same sign as \eqn{\rho} and is approximately
#' \eqn{\rho(1 - \rho^2)/[2(n - 1)]} to leading order (Fisher, 1915;
#' Hotelling, 1953; Ghosh, 1966); the exact formula above incorporates all
#' higher-order corrections. For \eqn{\rho = 0.5}, \eqn{n = 10}, the bias
#' is about \eqn{+0.021}; for \eqn{\rho = 0.5}, \eqn{n = 30}, about
#' \eqn{+0.0065}.
#'
#' \strong{Olkin-Pratt (1958) unbiased estimator of \eqn{\rho}.} The
#' companion to this expected-value calculation is the Olkin-Pratt
#' unbiased estimator of \eqn{\rho} given an observed \eqn{r}:
#' \deqn{\tilde{\rho}_{\text{OP}}(r, n) \;=\;
#'   r \,\cdot\, {}_2F_1\!\left(\tfrac{1}{2},\, \tfrac{1}{2};\, \tfrac{n-2}{2};\, 1 - r^2 \right).}
#' Olkin & Pratt (1958) prove
#' \eqn{\mathrm{E}[\tilde{\rho}_{\text{OP}}(r, n) \mid \rho, n] = \rho}
#' exactly, for every \eqn{\rho \in (-1, 1)} and \eqn{n \ge 4}. The
#' commonly quoted first-order approximation
#' \eqn{\tilde{\rho} \approx r\,[1 + (1 - r^2)/(2(n - 3))]} (Olkin, 1967)
#' is the truncation of the OP series at \eqn{k = 1}; the full series is
#' what makes the estimator exactly unbiased.
#'
#' Although the unbiased estimator is straightforward to apply, it is
#' rarely used because in most downstream uses (significance testing,
#' Fisher's \eqn{Z} confidence intervals, structural-equation models) the
#' bias is small relative to other sources of uncertainty. Where it does
#' matter, meta-analyses with many small samples, reliability /
#' validity coefficients estimated from short calibration samples,
#' design-stage estimates feeding into AIPE sample size machinery
#' the correction is well worth applying.
#'
#' \strong{Connection to Fisher's Z transform.} The variance-stabilizing
#' transform \eqn{z = \tanh^{-1}(r)}, proposed in passing in Fisher (1915,
#' p. 521) and developed in Fisher (1921), has approximate variance
#' \eqn{1/(n-3)} regardless of \eqn{\rho}, but \eqn{\mathrm{E}[z]} also
#' carries a small-sample bias of order \eqn{1/n} (Hotelling, 1953,
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
#' # 1. A single value: rho = 0.5, n = 10. The sample r is downwardly
#' #        biased by about 0.021, roughly 4% of rho.
#' expected_r(rho = 0.5, n = 10)
#'
#' # 2. Bias as a function of n for fixed rho. The bias is roughly
#' #        rho * (1 - rho^2) / (2(n - 1)) to leading order; as n grows it
#' #        shrinks toward zero.
#' expected_r(rho = 0.5, n = c(5, 10, 20, 50, 100, 500))
#'
#' # 3. Bias as a function of rho for fixed n. The bias is zero at
#' #        rho = 0 and rho = +/- 1, and largest near rho = +/- 0.6.
#' expected_r(rho = seq(0, 0.95, by = 0.05), n = 10)
#'
#' # 4. The Olkin-Pratt unbiased estimator: invert the bias for an
#' #        observed sample r.
#' set.seed(113)
#' x <- rnorm(20); y <- 0.4 * x + sqrt(1 - 0.4^2) * rnorm(20)
#' r_obs <- cor(x, y)
#' r_obs
#'
#' # Olkin-Pratt unbiased estimator of rho:
#' op_unbiased <- function(r, n, tol = 1e-15, max_iter = 5000) {
#'   z <- 1 - r^2; c_par <- (n - 2) / 2
#'   s <- 1; term <- 1
#'   for (k in seq_len(max_iter)) {
#'     term <- term * ((k - 0.5)^2) / ((c_par + k - 1) * k) * z
#'     s <- s + term
#'     if (abs(term) < tol * abs(s)) break
#'   }
#'   r * s
#' }
#' op_unbiased(r_obs, n = 20)
#'
#' # 5. Design-stage use: a study planned around rho = 0.4 with n = 30.
#' #        Naive plug-in says we expect to observe r = 0.40 on average,
#' #        but the realized expected r is smaller, and that gap matters
#' #        for any precision- or power-based sample size calculation that
#' #        plugs in rho as if it were the expected sample r.
#' rho_planned <- 0.4
#' n_planned   <- 30
#' expected_r(rho = rho_planned, n = n_planned)
#'
#' # 6. Tuning constants are hidden from the signature but tunable
#' #        through options() for the rare cases that need them (e.g.,
#' #        very near-boundary |rho| where the 2F1 series converges
#' #        slowly). The defaults rarely need to be changed.
#' options(DMAR.expected_r.tol = 1e-12,
#'         DMAR.expected_r.max_iter = 20000L)
#' expected_r(rho = 0.999, n = 5)
#' options(DMAR.expected_r.tol = NULL,
#'         DMAR.expected_r.max_iter = NULL)   # restore defaults
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest correlation
#'
#' @family effect size estimates
#'
#' @export

expected_r <- function(rho, n) {
  if (!is.numeric(rho)) stop("'rho' must be numeric.")
  if (!is.numeric(n)) stop("'n' must be numeric.")
  if (any(abs(rho) > 1, na.rm = TRUE))
    stop("'rho' must lie in [-1, 1].")
  if (any(n < 4, na.rm = TRUE))
    stop("'n' must be >= 4 for well-conditioned computation of E[r].")

  # Hidden tuning constants for the 2F1 series; documented in @details.
  # Pulled from options() so they are controllable but not visible in the
  # signature; sensible defaults that almost never need to be changed.
  tol      <- getOption("DMAR.expected_r.tol",      1e-15)
  max_iter <- as.integer(getOption("DMAR.expected_r.max_iter", 5000L))

  args <- data.frame(rho = rho, n = n)
  out <- vapply(seq_len(nrow(args)), function(i) {
    .expected_r_one(args$rho[i], args$n[i], tol = tol, max_iter = max_iter)
  }, numeric(1))

  .as_dmar_tbl(data.frame(
    rho           = args$rho,
    n             = args$n,
    expected_r    = out,
    bias          = args$rho - out,
    relative_bias = ifelse(args$rho == 0, NA_real_,
                           (args$rho - out) / args$rho),
    stringsAsFactors = FALSE,
    row.names     = NULL
  ))
}


# Single-(rho, n) worker. Computes the exact Olkin-Pratt (Hotelling, 1953)
# expectation E[r | rho, n] = rho * 2F1(1/2, 1/2; (n+1)/2; rho^2) *
# Gamma(n/2)^2 / (Gamma((n-1)/2) * Gamma((n+1)/2)) under bivariate normality.
.expected_r_one <- function(rho, n, tol, max_iter) {
  if (is.na(rho) || is.na(n)) return(NA_real_)
  if (abs(rho) == 1) return(rho)  # boundary: r is degenerate at +/-1

  z <- rho^2
  if (z == 0) return(0)  # rho = 0 => E[r] = 0 exactly

  # 2F1(1/2, 1/2; (n+1)/2; rho^2) by forward recurrence:
  # the standard 2F1 ratio of consecutive terms is
  #   term_{k+1} / term_k = (a + k)(b + k) / ((c + k)(k + 1)) * z,
  # which for a = b = 1/2, c = (n+1)/2 simplifies to
  #   (k - 1/2)^2 / (((n+1)/2 + k - 1) * k) * z
  # when iterating from term_{k-1} to term_k.
  c_par <- (n + 1) / 2
  s <- 1
  term <- 1
  converged <- FALSE
  for (k in seq_len(max_iter)) {
    term <- term * ((k - 0.5)^2) / ((c_par + k - 1) * k) * z
    s <- s + term
    if (abs(term) < tol * abs(s)) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    warning(sprintf(
      "expected_r(): 2F1 series did not converge to tol = %.1e in %d iterations (rho = %.3f, n = %d). Returning current partial sum.",
      tol, max_iter, rho, n
    ))
  }

  # Gamma ratio: Gamma(n/2)^2 / (Gamma((n-1)/2) * Gamma((n+1)/2)).
  # Compute on log scale for numerical stability at large n.
  gamma_factor <- exp(2 * lgamma(n / 2) -
                      lgamma((n - 1) / 2) -
                      lgamma((n + 1) / 2))

  rho * s * gamma_factor
}
