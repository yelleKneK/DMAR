# TOST on the Pearson correlation.
#' Equivalence Test for the Pearson Correlation via Two One-Sided Tests (TOST)
#'
#' Performs a two one-sided tests procedure for equivalence of a
#' Pearson correlation \eqn{\rho} to zero against user-specified
#' equivalence bounds \eqn{[-\rho_L, \rho_U]} (Counsell & Cribbie,
#' 2015; Goertzen & Cribbie, 2010). Uses the Fisher's \eqn{Z}
#' transformation throughout, a large-sample approximation whose
#' accuracy under bivariate normality improves quickly with \emph{n}.
#' Equivalence is declared when the 100(1 - 2\eqn{\alpha})\%
#' Fisher's \eqn{Z} CI on \eqn{\rho} lies entirely inside the
#' equivalence region.
#'
#' @param r,n Observed sample correlation \emph{r} and sample size.
#'   Alternatively supply \code{x} and \code{y} to compute \emph{r}
#'   from raw data.
#' @param x,y Numeric vectors of paired observations. If supplied,
#'   \code{r} and \code{n} are computed from the data and the
#'   \code{r}/\code{n} arguments are ignored.
#' @param rho_lower,rho_upper Equivalence bounds on the correlation
#'   scale, both positive. The equivalence region is
#'   \eqn{[-\rho_L, +\rho_U]}. If only \code{rho_upper} is supplied,
#'   the bounds are symmetric.
#' @param alpha_level One-sided significance level. Default \code{0.05}.
#'
#' @return A \code{data.frame} with rows for the observed \emph{r},
#'   the two one-sided test statistics on the Fisher's \eqn{Z} scale,
#'   their \emph{p}-values, the joint TOST \emph{p}-value, the
#'   100(1 - 2\eqn{\alpha})\% CI on \eqn{\rho}, the equivalence bounds,
#'   a binary equivalence flag, and the sample size (\code{n}).
#'
#' @details
#' \strong{Fisher's \eqn{Z} transformation.}
#' \deqn{Z = \tfrac{1}{2} \log\left(\frac{1 + r}{1 - r}\right), \quad
#'   \mathrm{Var}(Z) = \frac{1}{n - 3}.}
#' The TOST is run on the Fisher's \eqn{Z} scale:
#' \itemize{
#'   \item Lower test: \eqn{(Z - Z_{-\rho_L}) \sqrt{n - 3}} compared
#'         against the upper \eqn{\alpha} of \eqn{N(0, 1)}.
#'   \item Upper test: \eqn{(Z - Z_{\rho_U}) \sqrt{n - 3}} compared
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
#' # 1. Equivalence test that |rho| < 0.10 with n = 200 and r = 0.05:
#' equivalence_r(r = 0.05, n = 200, rho_upper = 0.10)
#'
#' # 2. From raw data:
#' set.seed(113)
#' x <- rnorm(150); y <- 0.04 * x + rnorm(150)
#' equivalence_r(x = x, y = y, rho_upper = 0.15)
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

equivalence_r <- function(r = NULL, n = NULL,
                   x = NULL, y = NULL,
                   rho_lower = NULL, rho_upper = NULL,
                   alpha_level = 0.05) {
  if (is.null(rho_upper))
    stop("'rho_upper' must be specified.")
  if (!is.numeric(rho_upper) || rho_upper <= 0 || rho_upper >= 1)
    stop("'rho_upper' must be in (0, 1).")
  if (is.null(rho_lower)) rho_lower <- rho_upper
  if (!is.numeric(rho_lower) || rho_lower <= 0 || rho_lower >= 1)
    stop("'rho_lower' must be in (0, 1).")
  if (alpha_level <= 0 || alpha_level >= 0.5)
    stop("'alpha_level' must be in (0, 0.5).")

  if (!is.null(x) && !is.null(y)) {
    if (!is.numeric(x) || !is.numeric(y))
      stop("'x' and 'y' must be numeric vectors.")
    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length.")
    ok <- !is.na(x) & !is.na(y)
    x <- x[ok]; y <- y[ok]
    n <- length(x)
    if (n < 4)
      stop("'x' and 'y' must supply at least 4 complete pairs.")
    r <- stats::cor(x, y)
  } else if (is.null(r) || is.null(n)) {
    stop("Supply either ('x', 'y') or ('r', 'n').")
  } else {
    if (!is.numeric(r) || abs(r) >= 1)
      stop("'r' must be in (-1, 1).")
    if (!is.numeric(n) || n < 4)
      stop("'n' must be >= 4.")
  }

  z   <- atanh(r)
  se  <- 1 / sqrt(n - 3)
  z_l <- atanh(-rho_lower); z_u <- atanh(rho_upper)

  z_test_lower <- (z - z_l) / se   # against H0: rho <= -rho_lower
  z_test_upper <- (z - z_u) / se   # against H0: rho >= rho_upper
  p_lower <- stats::pnorm(z_test_lower, lower.tail = FALSE)
  p_upper <- stats::pnorm(z_test_upper, lower.tail = TRUE)
  p_tost  <- max(p_lower, p_upper)

  z_crit  <- stats::qnorm(1 - alpha_level)
  ci_z_lo <- z - z_crit * se
  ci_z_hi <- z + z_crit * se
  ci_lo   <- tanh(ci_z_lo)
  ci_hi   <- tanh(ci_z_hi)

  equivalent <- as.integer((ci_lo > -rho_lower) && (ci_hi < rho_upper))

  out <- data.frame(
    term  = c("r", "z_lower_test", "z_upper_test",
              "p_lower", "p_upper", "p_tost",
              "lower_limit", "upper_limit",
              "rho_lower", "rho_upper",
              "equivalent", "n"),
    value = c(r, z_test_lower, z_test_upper,
              p_lower, p_upper, p_tost,
              ci_lo, ci_hi,
              -rho_lower, rho_upper,
              equivalent, n),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, p_terms = c("p_lower", "p_upper", "p_tost"),
               conf_level = 1 - 2 * alpha_level)
}
