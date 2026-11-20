# TOST on the Pearson correlation.
#' Equivalence Test for the Pearson Correlation via Two One-Sided Tests (TOST)
#'
#' Performs a two one-sided tests procedure for equivalence of a
#' Pearson correlation \eqn{\rho} to zero against user-specified
#' equivalence bounds \eqn{[-\rho_L, \rho_U]} (Counsell & Cribbie,
#' 2015; Goertzen & Cribbie, 2010). Uses the Fisher's \eqn{Z}
#' transformation throughout, a large-sample approximation whose
#' accuracy under bivariate normality improves quickly with \emph{N}.
#' Equivalence is declared when the 100(1 - 2\eqn{\alpha})\%
#' Fisher's \eqn{Z} CI on \eqn{\rho} lies entirely inside the
#' equivalence region.
#'
#' @param r,N Observed sample correlation \emph{r} and sample size.
#'   Alternatively supply \code{x} and \code{y} to compute \emph{r}
#'   from raw data.
#' @param x,y Numeric vectors of paired observations. If supplied,
#'   \code{r} and \code{N} are computed from the data and the
#'   \code{r}/\code{N} arguments are ignored.
#' @param cor_lower,cor_upper Equivalence bounds on the correlation
#'   scale, both positive. The equivalence region is
#'   \eqn{[-\rho_L, +\rho_U]}, where \eqn{\rho_L} is \code{cor_lower} and
#'   \eqn{\rho_U} is \code{cor_upper}. If only \code{cor_upper} is
#'   supplied, the bounds are symmetric.
#' @param alpha_level One-sided significance level. Default \code{0.05}.
#'
#' @return A \code{data.frame} with rows for the observed \emph{r},
#'   the two one-sided test statistics on the Fisher's \eqn{Z} scale,
#'   their \emph{p}-values, the joint TOST \emph{p}-value, the
#'   100(1 - 2\eqn{\alpha})\% CI on \eqn{\rho}, the equivalence bounds,
#'   a binary equivalence flag, and the sample size (\code{N}). The
#'   table names the TOST \emph{p}-value, \code{p_tost}, in its
#'   \code{p_sentence} attribute, so that \code{\link{results_sentence}}
#'   pairs the interval with the test of equivalence rather than with
#'   \code{p_lower}, the first \emph{p}-value row.
#'
#' @details
#' \strong{Fisher's \eqn{Z} transformation.}
#' \deqn{Z = \tfrac{1}{2} \log\left(\frac{1 + r}{1 - r}\right), \quad
#'   \mathrm{Var}(Z) = \frac{1}{N - 3}.}
#' The TOST is run on the Fisher's \eqn{Z} scale:
#' \itemize{
#'   \item Lower test: \eqn{(Z - Z_{-\rho_L}) \sqrt{N - 3}} compared
#'         against the upper \eqn{\alpha} of \eqn{N(0, 1)}.
#'   \item Upper test: \eqn{(Z - Z_{\rho_U}) \sqrt{N - 3}} compared
#'         against the lower \eqn{\alpha} of \eqn{N(0, 1)}.
#' }
#' The CI bounds are back-transformed from the \eqn{Z} scale via
#' \eqn{r = \tanh(Z)} so they remain in \eqn{[-1, 1]}.
#'
#' \strong{Choosing \eqn{\rho_L} and \eqn{\rho_U}.} Common choices in
#' psychology are \eqn{0.1}, or domain-specific meaningfulness
#' thresholds (e.g., \eqn{0.2} for cognitive task correlations). The
#' bounds must be set \emph{before} data collection.
#'
#' @references
#' Counsell, A., & Cribbie, R. A. (2015). Equivalence tests for
#'   comparing correlation and regression coefficients. \emph{British
#'   Journal of Mathematical and Statistical Psychology, 68}(2),
#'   292--309. \doi{10.1111/bmsp.12045}
#'
#' Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of
#'   association: An equivalence testing approach. \emph{British
#'   Journal of Mathematical and Statistical Psychology, 63}(3),
#'   527--537. \doi{10.1348/000711009X475853}
#'
#' Lakens, D. (2017). Equivalence tests: A practical primer for
#'   \emph{t} tests, correlations, and meta-analyses. \emph{Social
#'   Psychological and Personality Science, 8}(4), 355--362.
#'   \doi{10.1177/1948550617697177}
#'
#' @seealso \code{\link{equivalence_smd}}, \code{\link{ci_r}}
#'
#' @examples
#' # 1. Equivalence test that |rho| < 0.10 with N = 200 and r = 0.05:
#' equivalence_r(r = 0.05, N = 200, cor_upper = 0.10)
#'
#' # 2. From raw data:
#' set.seed(113)
#' x <- rnorm(150); y <- 0.04 * x + rnorm(150)
#' equivalence_r(x = x, y = y, cor_upper = 0.15)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords htest
#'
#' @family hypothesis tests
#'
#' @family equivalence testing
#'
#' @export

equivalence_r <- function(r = NULL, N = NULL,
                   x = NULL, y = NULL,
                   cor_lower = NULL, cor_upper = NULL,
                   alpha_level = 0.05) {
  .dmar_check_alpha_level(alpha_level)
  if (is.null(cor_upper))
    .dmar_refuse_argument("'cor_upper' must be specified.", arg = "cor_upper")
  if (!is.numeric(cor_upper) || cor_upper <= 0 || cor_upper >= 1)
    .dmar_refuse_argument("'cor_upper' must be in (0, 1).", arg = "cor_upper")
  if (is.null(cor_lower)) cor_lower <- cor_upper
  if (!is.numeric(cor_lower) || cor_lower <= 0 || cor_lower >= 1)
    .dmar_refuse_argument("'cor_lower' must be in (0, 1).", arg = "cor_lower")
  if (alpha_level <= 0 || alpha_level >= 0.5)
    .dmar_refuse_argument(
      "'alpha_level' must be in (0, 0.5).", arg = "alpha_level")

  if (!is.null(x) && !is.null(y)) {
    if (!is.numeric(x) || !is.numeric(y))
      .dmar_refuse_argument(
        "'x' and 'y' must be numeric vectors.", arg = c("x", "y"))
    if (length(x) != length(y))
      .dmar_refuse_argument(
        "'x' and 'y' must have the same length.", arg = c("x", "y"))
    ok <- !is.na(x) & !is.na(y)
    x <- x[ok]; y <- y[ok]
    N <- length(x)
    if (N < 4)
      .dmar_refuse_argument(
        "'x' and 'y' must supply at least 4 complete pairs.", arg = c("x", "y"))
    r <- stats::cor(x, y)
  } else if (is.null(r) || is.null(N)) {
    .dmar_refuse_argument(
      "Supply either ('x', 'y') or ('r', 'N').", arg = c("x", "y", "r"))
  } else {
    if (!is.numeric(r) || abs(r) >= 1)
      .dmar_refuse_argument("'r' must be in (-1, 1).", arg = "r")
    if (!is.numeric(N) || N < 4)
      .dmar_refuse_argument("'N' must be >= 4.", arg = "N")
  }

  z   <- atanh(r)
  se  <- 1 / sqrt(N - 3)
  z_l <- atanh(-cor_lower); z_u <- atanh(cor_upper)

  z_test_lower <- (z - z_l) / se   # against H0: rho <= -cor_lower
  z_test_upper <- (z - z_u) / se   # against H0: rho >= cor_upper
  p_lower <- stats::pnorm(z_test_lower, lower.tail = FALSE)
  p_upper <- stats::pnorm(z_test_upper, lower.tail = TRUE)
  p_tost  <- max(p_lower, p_upper)

  z_crit  <- stats::qnorm(1 - alpha_level)
  ci_z_lo <- z - z_crit * se
  ci_z_hi <- z + z_crit * se
  ci_lo   <- tanh(ci_z_lo)
  ci_hi   <- tanh(ci_z_hi)

  equivalent <- as.integer((ci_lo > -cor_lower) && (ci_hi < cor_upper))

  out <- data.frame(
    term  = c("r", "z_lower_test", "z_upper_test",
              "p_lower", "p_upper", "p_tost",
              "ci_lower", "ci_upper",
              "cor_lower", "cor_upper",
              "equivalent", "N"),
    value = c(r, z_test_lower, z_test_upper,
              p_lower, p_upper, p_tost,
              ci_lo, ci_hi,
              -cor_lower, cor_upper,
              equivalent, N),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- .as_dmar_tbl(out, p_terms = c("p_lower", "p_upper", "p_tost"),
                      conf_level = 1 - 2 * alpha_level)
  # The p-value a results sentence reports beside the 1 - 2 alpha
  # interval is that of the equivalence test, not the first p_terms row.
  attr(out, "p_sentence") <- "p_tost"
  out
}
