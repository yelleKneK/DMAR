# Exact expected value of the sample partial correlation given rho, N, p.
#' Exact Expected Value of the Sample Partial Correlation
#'
#' Computes \eqn{\mathrm{E}[r_{XY \cdot Z_1 \cdots Z_{p-1}} \mid \rho, N, p]},
#' the exact expected value of the sample partial Pearson correlation
#' coefficient under multivariate normality, applying the Olkin-Pratt
#' (1958) / Hotelling (1953) result for the simple Pearson correlation to
#' the partial-correlation setting, an extension Olkin and Pratt (1958,
#' Section 2.3) give themselves. Like the simple
#' \emph{r}, the partial \emph{r} is downward-biased as an estimator of
#' the population partial correlation \eqn{\rho_{XY \cdot Z}}, with the
#' magnitude of the bias growing in the number of predictors \eqn{p} and
#' shrinking in the sample size \eqn{N}.
#'
#' The formula is the same as for the simple \emph{r} but with the
#' degrees of freedom reduced from \eqn{N} to \eqn{N - (p - 1)}:
#' \deqn{\mathrm{E}[r_{XY \cdot Z} \mid \rho_{XY \cdot Z}, N, p] \;=\;
#'   \rho_{XY \cdot Z} \,\cdot\,
#'   {}_2F_1\!\left(\tfrac{1}{2},\, \tfrac{1}{2};\,
#'                  \tfrac{N - p + 2}{2};\,
#'                  \rho_{XY \cdot Z}^{\,2}\right)
#'   \,\cdot\,
#'   \frac{\Gamma\!\left(\tfrac{N - p + 1}{2}\right)^2}
#'        {\Gamma\!\left(\tfrac{N - p}{2}\right)\,
#'         \Gamma\!\left(\tfrac{N - p + 2}{2}\right)}.}
#'
#' \code{expected_partial_r()} is especially useful at the \emph{design
#' stage}: sample size plans that assume \eqn{\rho_{XY \cdot Z}} as the
#' value to be observed will under-deliver on the expected width or power
#' of a CI on the partial correlation by an amount that grows in \eqn{p}.
#'
#' @param true_cor Population partial correlation
#'   \eqn{\rho_{XY \cdot Z_1 \cdots Z_{p-1}}}. Scalar or vector in \eqn{[-1, 1]}.
#' @param N Total sample size; \eqn{N - p \ge 3} is required for
#'   well-conditioned computation. Scalar or vector.
#' @param p Number of predictors in the regression of \emph{Y} on \emph{X}
#'   and the variables partialed out, so \code{p - 1} variables are
#'   partialed out and \code{p} is at least 2. A categorical variable with
#'   \emph{a} levels counts as \emph{a} - 1 predictors. It is the \code{p}
#'   of \code{\link{ss_power_R2}}, \code{\link{ci_R2}}, and
#'   \code{\link{ci_reg_coef}}: the correlation concerns one predictor of a
#'   regression with \code{p} predictors. In the notation of a model
#'   comparison, \eqn{q} is the number of predictors removed from the full
#'   model to give the restricted model (\eqn{q = p_{full} - p_{restricted}},
#'   the predictors being tested), and the number of variables partialed out
#'   is \eqn{p - q}; this correlation concerns one predictor, so \eqn{q = 1}
#'   and \eqn{p - 1} variables are partialed out. Scalar or vector.
#'
#' @return A \code{data.frame} with one row per (\code{true_cor},
#'   \code{N}, \code{p}) input and the columns \code{true_cor}, \code{N},
#'   \code{p}, \code{expected_partial_r}, \code{bias}, and
#'   \code{relative_bias}.
#'
#' @details
#' \strong{Generalization of Olkin-Pratt.} Under multivariate normality,
#' the partial correlation \eqn{r_{XY \cdot Z}} computed from a sample of
#' size \eqn{N} has the same sampling distribution as a simple Pearson
#' correlation from a sample of size \eqn{N - (p - 1)} (Anderson, 2003, Theorem
#' 4.3.5, Section 4.3.2, p. 143, who credits the derivation to Fisher).
#' The bias formula for the simple \emph{r} (Hotelling, 1953;
#' Olkin & Pratt, 1958) therefore applies directly to the partial \emph{r}
#' with the substitution \eqn{N \to N - (p - 1)}. The Olkin-Pratt (1958) unbiased
#' estimator extends in the same way: given an observed \eqn{r_{XY \cdot Z}},
#' the unbiased estimator of \eqn{\rho_{XY \cdot Z}} is
#' \eqn{r \cdot {}_2F_1(1/2, 1/2; (N - p - 1)/2; 1 - r^2)} (Olkin &
#' Pratt, 1958, Section 2.3, where the partial correlation is shown to
#' have the simple correlation's density at the reduced sample size, so
#' their estimator applies with that substitution).
#'
#' \strong{Magnitude of the bias.} To leading order, from Hotelling's
#' (1953) expansion of the simple-correlation expectation applied at the
#' reduced sample size:
#' \deqn{\rho_{XY \cdot Z} - \mathrm{E}[r_{XY \cdot Z}] \;\approx\;
#'   \rho_{XY \cdot Z} (1 - \rho_{XY \cdot Z}^2) / [2 (N - p)],}
#' matching the sign convention of the returned \code{bias} column,
#' \eqn{\rho_{XY \cdot Z} - \mathrm{E}[r_{XY \cdot Z}]}, which is
#' positive for positive \eqn{\rho_{XY \cdot Z}}. For
#' \eqn{\rho_{XY \cdot Z} = 0.4}, \eqn{N = 30}, \eqn{p = 3}, the exact
#' bias is about \eqn{+0.00624}; for \eqn{N = 30}, \eqn{p = 11}, about
#' \eqn{+0.00888}.
#'
#' \strong{Tuning the series convergence.} The underlying \eqn{{}_2F_1}
#' series is summed by forward recurrence and stops when the relative
#' contribution falls below \code{getOption("dmar.expected_r.tol",
#' 1e-15)}, with a maximum of
#' \code{getOption("dmar.expected_r.max_iter", 5000L)} terms. These are
#' the same options as for \code{\link{expected_r}}; users almost never
#' need to change them.
#'
#' @references
#' Anderson, T. W. (2003). \emph{An introduction to multivariate
#'   statistical analysis} (3rd ed.), Sections 4.2 and 4.3. Wiley.
#'   (Theorem 4.3.5, Section 4.3.2, p. 143: the sample partial correlation
#'   based on \eqn{N} observations, with \eqn{k} variables partialed out,
#'   has cdf \eqn{F[r \mid N - k, \rho]}; Anderson writes the count of
#'   partialed variables as \eqn{p - q} in his notation, and here it is
#'   \eqn{p - 1}.)
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
#' # 1. A single value: true_cor = 0.4, N = 30, three variables partialed out
#' #    (p = 4 predictors: X and the three variables).
#' expected_partial_r(true_cor = 0.4, N = 30, p = 4)
#'
#' # 2. Bias grows with p at fixed (true_cor, N):
#' expected_partial_r(true_cor = 0.4, N = 30, p = c(2, 3, 6, 11, 21))
#'
#' # 3. Partialing one variable (p = 2) at N = 11 is distributed exactly as
#' #    a simple correlation at N = 10 (Anderson, 2003, Theorem 4.3.5):
#' expected_partial_r(true_cor = 0.5, N = 11, p = 2)$expected_partial_r
#' expected_r(true_cor = 0.5, N = 10)$expected_r
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest correlation
#'
#' @family effect size estimates
#'
#' @export

expected_partial_r <- function(true_cor, N, p) {
  if (missing(p))
    .dmar_refuse_argument(
      "'p' (the number of predictors) is required.", arg = "p")
  if (!is.numeric(p)) .dmar_refuse_argument("'p' must be numeric.", arg = "p")
  p <- .check_partial_p(p, scalar = FALSE)
  if (!is.numeric(true_cor)) .dmar_refuse_argument(
    "'true_cor' must be numeric.", arg = "true_cor")
  if (!is.numeric(N))   .dmar_refuse_argument("'N' must be numeric.", arg = "N")
  if (any(abs(true_cor) > 1, na.rm = TRUE))
    .dmar_refuse_argument("'true_cor' must lie in [-1, 1].", arg = "true_cor")
  if (any(N - p < 3, na.rm = TRUE))
    .dmar_refuse_argument(
      "'N - p' must be at least 3 for well-conditioned computation.",
      arg = NULL)

  # The partial-r sampling distribution at sample size N with p - 1
  # variables partialed out is identical to the simple-r sampling
  # distribution at sample size N - (p - 1) (Anderson, 2003). So delegate
  # to expected_r() at that sample size.
  args <- data.frame(true_cor = true_cor, N = N, p = p)
  res  <- expected_r(true_cor = args$true_cor, N = args$N - (args$p - 1))
  .as_dmar_tbl(data.frame(
    true_cor            = args$true_cor,
    N                   = args$N,
    p                   = args$p,
    expected_partial_r  = res$expected_r,
    bias                = res$bias,
    relative_bias       = res$relative_bias,
    stringsAsFactors = FALSE,
    row.names    = NULL
  ))
}
