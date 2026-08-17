#' Sample Size or Power for a Pearson Correlation Coefficient (Fisher Z Transformation)
#'
#' Determine the necessary sample size to achieve a desired level of statistical power for the
#' test of a Pearson correlation against a null value (typically zero), or, given a sample size,
#' return the realized statistical power. The computation uses the Fisher's \emph{Z} transformation,
#' which has a near-normal sampling distribution.
#'
#' @param rho The population correlation coefficient under the alternative hypothesis
#' @param rho_0 The null hypothesis value of the correlation (default 0)
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param N Sample size (\emph{number of pairs}); if specified, returns the realized power (the ss_power_* family is not uniform here: \code{\link{ss_power_contrast}} takes a \emph{per-group} size)
#' @param directional Logical: \code{TRUE} for a one-sided test (in the same direction as the difference \code{rho - rho_0}), \code{FALSE} (default) for a two-sided test
#'
#' @details
#' Under the alternative the Fisher-transformed correlation \eqn{Z_r = \tanh^{-1}(r)} is
#' approximately normal with mean \eqn{Z_\rho = \tanh^{-1}(\rho)} and variance \eqn{1 / (N - 3)}.
#' Power is computed from this normal approximation.
#'
#' For sample size, a closed-form expression is used as the starting point,
#' \deqn{N = ((z_{\alpha} + z_{\beta}) / (Z_\rho - Z_{\rho_0}))^2 + 3,}
#' which is then verified iteratively to ensure power exactly meets or exceeds \code{desired_power}.
#' The search is bounded at \eqn{N = 10^7}. When \code{rho} and \code{rho_0} are so close that
#' \code{desired_power} is unreachable within that bound, the function stops with an error rather
#' than searching indefinitely.
#'
#' @return A \code{data.frame} with rows for \code{necessary_N} (or \code{specified_N}),
#'   \code{actual_power}, \code{rho}, and \code{rho_0}.
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Fisher, R. A. (1921). On the "probable error" of a coefficient of correlation deduced from a small sample. \emph{Metron, 1}, 3--32.
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 3 on the one-way
#'   ANOVA and Chapter 4 on contrasts.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_r}}, \code{\link{convert_r_Z}}, \code{\link{convert_Z_r}}
#'
#' @examples
#' # Population r = 0.30, null r = 0, desired power = .80, two-sided
#' ss_power_r(rho = 0.30, desired_power = 0.80)
#'
#' # Same with a directional alternative
#' ss_power_r(rho = 0.30, desired_power = 0.80, directional = TRUE)
#'
#' # Realized power for N = 100 pairs
#' ss_power_r(rho = 0.30, N = 100)
#'
#' # Test against a non-zero null (rho_0 = 0.20) -- looking for evidence rho > 0.20
#' ss_power_r(rho = 0.40, rho_0 = 0.20, desired_power = 0.80, directional = TRUE)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export
ss_power_r <- function(rho, rho_0 = 0, desired_power = 0.85, alpha_level = 0.05,
                       N = NULL, directional = FALSE) {
  if (!is.numeric(rho) || length(rho) != 1L || rho <= -1 || rho >= 1) {
    stop("'rho' must be a single numeric value in (-1, 1).", call. = FALSE)
  }
  if (!is.numeric(rho_0) || length(rho_0) != 1L || rho_0 <= -1 || rho_0 >= 1) {
    stop("'rho_0' must be a single numeric value in (-1, 1).", call. = FALSE)
  }
  if (rho == rho_0) {
    stop("'rho' equals 'rho_0'; there is no effect to detect.", call. = FALSE)
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be in (0, 1).", call. = FALSE)
  }

  z_rho   <- atanh(rho)
  z_rho_0 <- atanh(rho_0)
  z_diff  <- z_rho - z_rho_0

  power_at <- function(N) {
    if (N <= 3) return(NA_real_)
    se <- 1 / sqrt(N - 3)
    if (directional) {
      crit <- qnorm(1 - alpha_level)
      pnorm((abs(z_diff) - crit * se) / se)
    } else {
      crit <- qnorm(1 - alpha_level / 2)
      # Two-sided: chance of exceeding +crit OR falling below -crit, under the alternative
      pnorm((abs(z_diff) - crit * se) / se) + pnorm((-abs(z_diff) - crit * se) / se)
    }
  }

  if (!is.null(N)) {
    if (N < 4) stop("'N' must be at least 4.", call. = FALSE)
    out <- data.frame(
      term  = c("specified_N", "actual_power", "rho", "rho_0"),
      value = c(N, power_at(N), rho, rho_0)
    )
    class(out) <- c("dmar_ss_power", "dmar_tbl", "data.frame")
    return(out)
  }

  if (desired_power <= 0 || desired_power >= 1) {
    stop("'desired_power' must be in (0, 1).", call. = FALSE)
  }

  z_alpha <- if (directional) qnorm(1 - alpha_level) else qnorm(1 - alpha_level / 2)
  z_beta  <- qnorm(desired_power)
  N_start <- ceiling(((z_alpha + z_beta) / abs(z_diff))^2 + 3)
  N_start <- max(N_start, 4)

  # Fail fast when the target cannot be reached at any realistic sample size,
  # for example when 'rho' sits within rounding error of 'rho_0'. Power is
  # monotone in N, so a bound that misses the target rules out every smaller
  # size. The bound also keeps the search in the range where it is meaningful:
  # a degenerate effect sends the closed-form start into the 1e24 range, where
  # N - 1 is not a distinct double and the walk-down step can never advance.
  N_max <- 1e7
  if (power_at(N_max) < desired_power) {
    stop("Could not reach 'desired_power' at a sample size up to ",
         format(N_max, scientific = FALSE), "; the difference between 'rho' ",
         "and 'rho_0' may be too small to detect at this power.", call. = FALSE)
  }

  N_i <- min(N_start, N_max)
  pwr <- power_at(N_i)
  # Walk up if needed (the closed form can be slightly low for two-sided)
  while (pwr < desired_power) {
    N_i <- N_i + 1
    pwr <- power_at(N_i)
    if (N_i > N_max) {
      stop("The search for 'desired_power' did not converge below ",
           format(N_max, scientific = FALSE), ".", call. = FALSE)
    }
  }
  # Walk down to find the smallest N satisfying power, under the same bound
  steps <- 0
  while (N_i > 4) {
    pwr_below <- power_at(N_i - 1)
    if (pwr_below < desired_power) break
    N_i <- N_i - 1
    pwr <- pwr_below
    steps <- steps + 1
    if (steps > N_max) {
      stop("The search for 'desired_power' did not converge within ",
           format(N_max, scientific = FALSE), " steps.", call. = FALSE)
    }
  }

  out <- data.frame(
    term  = c("necessary_N", "actual_power", "rho", "rho_0"),
    value = c(N_i, pwr, rho, rho_0)
  )
  class(out) <- c("dmar_ss_power", "dmar_tbl", "data.frame")
  out
}
