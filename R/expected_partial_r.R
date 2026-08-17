# Exact expected value of the sample partial correlation given rho, n, J.
#' Exact Expected Value of the Sample Partial Correlation
#'
#' Computes \eqn{\mathrm{E}[r_{XY \cdot Z_1 \cdots Z_J} \mid \rho, n, J]},
#' the exact expected value of the sample partial Pearson correlation
#' coefficient under multivariate normality, applying the Olkin-Pratt
#' (1958) / Hotelling (1953) result for the simple Pearson correlation to
#' the partial-correlation setting, an extension Olkin and Pratt (1958,
#' Section 2.3) give themselves. Like the simple
#' \emph{r}, the partial \emph{r} is downward-biased as an estimator of
#' the population partial correlation \eqn{\rho_{XY \cdot Z}}, with the
#' magnitude of the bias growing in the number of controls \eqn{J} and
#' shrinking in the sample size \eqn{n}.
#'
#' The formula is the same as for the simple \emph{r} but with the
#' degrees of freedom reduced from \eqn{n} to \eqn{n - J}:
#' \deqn{\mathrm{E}[r_{XY \cdot Z} \mid \rho_{XY \cdot Z}, n, J] \;=\;
#'   \rho_{XY \cdot Z} \,\cdot\,
#'   {}_2F_1\!\left(\tfrac{1}{2},\, \tfrac{1}{2};\,
#'                  \tfrac{n - J + 1}{2};\,
#'                  \rho_{XY \cdot Z}^{\,2}\right)
#'   \,\cdot\,
#'   \frac{\Gamma\!\left(\tfrac{n - J}{2}\right)^2}
#'        {\Gamma\!\left(\tfrac{n - J - 1}{2}\right)\,
#'         \Gamma\!\left(\tfrac{n - J + 1}{2}\right)}.}
#'
#' \code{expected_partial_r()} is especially useful at the \emph{design
#' stage}: sample size plans that assume \eqn{\rho_{XY \cdot Z}} as the
#' value to be observed will under-deliver on the expected width or power
#' of a CI on the partial correlation by an amount that grows in \eqn{J}.
#'
#' @param rho Population partial correlation
#'   \eqn{\rho_{XY \cdot Z_1 \cdots Z_J}}. Scalar or vector in \eqn{[-1, 1]}.
#' @param n Total sample size; \eqn{n - J - 1 \ge 3} is required for
#'   well-conditioned computation.
#' @param J Number of variables partialled out (the count of
#'   \eqn{Z_1, \ldots, Z_J}); at least 1.
#'
#' @return A \code{data.frame} with one row per (\code{rho},
#'   \code{n}, \code{J}) input and the columns \code{rho}, \code{n},
#'   \code{J}, \code{expected_partial_r}, \code{bias}, and
#'   \code{relative_bias}.
#'
#' @details
#' \strong{Generalization of Olkin-Pratt.} Under multivariate normality,
#' the partial correlation \eqn{r_{XY \cdot Z}} computed from a sample of
#' size \eqn{n} has the same sampling distribution as a simple Pearson
#' correlation from a sample of size \eqn{n - J} (Anderson, 2003, Theorem
#' 4.3.5, Section 4.3.2, p. 143, who credits the derivation to Fisher,
#' 1924). The bias formula for the simple \emph{r} (Hotelling, 1953;
#' Olkin & Pratt, 1958) therefore applies directly to the partial \emph{r}
#' with the substitution \eqn{n \to n - J}. The Olkin-Pratt (1958) unbiased
#' estimator extends in the same way: given an observed \eqn{r_{XY \cdot Z}},
#' the unbiased estimator of \eqn{\rho_{XY \cdot Z}} is
#' \eqn{r \cdot {}_2F_1(1/2, 1/2; (n - J - 2)/2; 1 - r^2)} (Olkin &
#' Pratt, 1958, Section 2.3, where the partial correlation is shown to
#' have the simple correlation's density at the reduced sample size, so
#' their estimator applies with that substitution).
#'
#' \strong{Magnitude of the bias.} To leading order, from Hotelling's
#' (1953) expansion of the simple-correlation expectation applied at the
#' reduced sample size:
#' \deqn{\rho_{XY \cdot Z} - \mathrm{E}[r_{XY \cdot Z}] \;\approx\;
#'   \rho_{XY \cdot Z} (1 - \rho_{XY \cdot Z}^2) / [2 (n - J - 1)],}
#' matching the sign convention of the returned \code{bias} column,
#' \eqn{\rho_{XY \cdot Z} - \mathrm{E}[r_{XY \cdot Z}]}, which is
#' positive for positive \eqn{\rho_{XY \cdot Z}}. For
#' \eqn{\rho_{XY \cdot Z} = 0.4}, \eqn{n = 30}, \eqn{J = 2}, the exact
#' bias is about \eqn{+0.00624}; for \eqn{n = 30}, \eqn{J = 10}, about
#' \eqn{+0.00888}.
#'
#' \strong{Tuning the series convergence.} The underlying \eqn{{}_2F_1}
#' series is summed by forward recurrence and stops when the relative
#' contribution falls below \code{getOption("DMAR.expected_r.tol",
#' 1e-15)}, with a maximum of
#' \code{getOption("DMAR.expected_r.max_iter", 5000L)} terms. These are
#' the same options as for \code{\link{expected_r}}; users almost never
#' need to change them.
#'
#' @references
#' Anderson, T. W. (2003). \emph{An introduction to multivariate
#'   statistical analysis} (3rd ed.), Sections 4.2 and 4.3. Wiley.
#'   (Theorem 4.3.5, Section 4.3.2, p. 143: the sample partial correlation
#'   based on \eqn{N} observations, with \eqn{J} variables partialled out,
#'   has cdf \eqn{F[r \mid N - J, \rho]}; Anderson writes the count of
#'   partialled variables as \eqn{p - q}.)
#'
#' Hotelling, H. (1953). New light on the correlation coefficient and its
#'   transforms. \emph{Journal of the Royal Statistical Society, Series B,
#'   15}(2), 193--232.
#'
#' Olkin, I., & Pratt, J. W. (1958). Unbiased estimation of certain
#'   correlation coefficients. \emph{The Annals of Mathematical Statistics,
#'   29}(1), 201--211.
#'
#' @seealso \code{\link{expected_r}}, \code{\link{var_partial_r}},
#'   \code{\link{ss_aipe_partial_r}}
#'
#' @examples
#' # 1. A single value: rho = 0.4, n = 30, J = 3 controls.
#' expected_partial_r(rho = 0.4, n = 30, J = 3)
#'
#' # 2. Bias grows with J at fixed (rho, n):
#' expected_partial_r(rho = 0.4, n = 30, J = c(1, 2, 5, 10, 20))
#'
#' # 3. Partialing J = 1 variable at n = 11 is distributed exactly as a
#' #    simple correlation at n = 10 (Anderson, 2003, Theorem 4.3.5):
#' expected_partial_r(rho = 0.5, n = 11, J = 1)$expected_partial_r
#' expected_r(rho = 0.5, n = 10)$expected_r
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest correlation
#'
#' @family effect size estimates
#'
#' @export

expected_partial_r <- function(rho, n, J) {
  if (!is.numeric(rho)) stop("'rho' must be numeric.")
  if (!is.numeric(n))   stop("'n' must be numeric.")
  if (!is.numeric(J))   stop("'J' must be numeric.")
  if (any(abs(rho) > 1, na.rm = TRUE))
    stop("'rho' must lie in [-1, 1].")
  if (any(J < 1, na.rm = TRUE))
    stop("'J' (number of controls) must be >= 1.")
  if (any(n - J < 4, na.rm = TRUE))
    stop("'n - J' must be >= 4 for well-conditioned computation.")

  # The partial-r sampling distribution at sample size n with J controls
  # is identical to the simple-r sampling distribution at sample size n - J
  # (Anderson, 2003). So delegate to expected_r() with n' = n - J.
  args <- data.frame(rho = rho, n = n, J = J)
  res  <- expected_r(rho = args$rho, n = args$n - args$J)
  .as_dmar_tbl(data.frame(
    rho                 = args$rho,
    n                   = args$n,
    J                   = args$J,
    expected_partial_r  = res$expected_r,
    bias                = res$bias,
    relative_bias       = ifelse(args$rho == 0, NA_real_,
                                 res$bias / args$rho),
    stringsAsFactors = FALSE,
    row.names    = NULL
  ))
}
