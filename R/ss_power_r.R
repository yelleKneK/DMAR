#' Sample Size or Power for a Pearson Correlation Coefficient (Fisher Z Transformation)
#'
#' Determine the necessary sample size to achieve a desired level of statistical power for the
#' test of a Pearson correlation against a null value (typically zero), or, given a sample size,
#' return the realized statistical power. The computation uses the Fisher's \emph{Z} transformation,
#' which has a near-normal sampling distribution.
#'
#' @param true_cor Population correlation coefficient the sample size or
#'   power is planned for (the value under the alternative hypothesis), in
#'   \eqn{(-1, 1)}. The planner treats this value as the truth.
#' @param null_cor The null hypothesis value of the correlation (default 0)
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param specified_N Total sample size (the \emph{number of pairs}); if
#'   specified, the function returns the power at that sample size.
#' @param directional Logical: \code{TRUE} for a one-sided test (in the same direction as the difference \code{true_cor - null_cor}), \code{FALSE} (default) for a two-sided test
#'
#' @details
#' Under the alternative the Fisher-transformed correlation \eqn{Z_r = \tanh^{-1}(r)} is
#' approximately normal with mean \eqn{Z_\rho = \tanh^{-1}(\rho)} and variance \eqn{1 / (N - 3)}.
#' Power is computed from this normal approximation.
#'
#' For sample size, a closed-form expression is used as the starting point,
#' \deqn{N = ((z_{\alpha} + z_{\beta}) / (Z_\rho - Z_{\rho_0}))^2 + 3,}
#' where \eqn{\rho} is \code{true_cor}, \eqn{\rho_0} is \code{null_cor},
#' \eqn{z_{\alpha}} is the standard normal critical value at \code{alpha_level}
#' (one- or two-sided, as \code{directional} says), and \eqn{z_{\beta}} is the
#' standard normal quantile at \code{desired_power},
#' which is then verified iteratively to ensure power exactly meets or exceeds \code{desired_power}.
#' The search is bounded at \eqn{N = 10^7}. When \code{true_cor} and \code{null_cor} are so close that
#' \code{desired_power} is unreachable within that bound, the function stops with an error rather
#' than searching indefinitely.
#'
#' @return A \code{data.frame} with rows for \code{necessary_N} (or \code{specified_N}),
#'   \code{actual_power}, \code{true_cor}, and \code{null_cor}.
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
#' # true_cor = 0.30, null_cor = 0, desired power = .80, two-sided
#' ss_power_r(true_cor = 0.30, desired_power = 0.80)
#'
#' # Same with a directional alternative
#' ss_power_r(true_cor = 0.30, desired_power = 0.80, directional = TRUE)
#'
#' # Realized power for N = 100 pairs
#' ss_power_r(true_cor = 0.30, specified_N = 100)
#'
#' # Test against a non-zero null (null_cor = 0.20), looking for evidence that rho > 0.20
#' ss_power_r(true_cor = 0.40, null_cor = 0.20, desired_power = 0.80, directional = TRUE)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export
ss_power_r <- function(true_cor, null_cor = 0, desired_power = 0.85,
                       alpha_level = 0.05, specified_N = NULL,
                       directional = FALSE) {
  .dmar_check_alpha_level(alpha_level)
  # The arithmetic below is written in rho, the planning value, and N.
  rho <- true_cor
  N <- specified_N
  if (!is.numeric(rho) || length(rho) != 1L || rho <= -1 || rho >= 1) {
    .dmar_refuse_argument(
      "'true_cor' must be a single numeric value in (-1, 1).", arg = "true_cor")
  }
  if (!is.numeric(null_cor) || length(null_cor) != 1L || null_cor <= -1 || null_cor >= 1) {
    .dmar_refuse_argument(
      "'null_cor' must be a single numeric value in (-1, 1).", arg = "null_cor")
  }
  if (rho == null_cor) {
    .dmar_refuse_argument(
      "'true_cor' equals 'null_cor'; there is no effect to detect.",
      arg = c("true_cor", "null_cor"))
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    .dmar_refuse_argument(
      "'alpha_level' must be in (0, 1).", arg = "alpha_level")
  }

  z_rho   <- atanh(rho)
  z_rho_0 <- atanh(null_cor)
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

  # An invalid desired_power is refused whether or not it is used: with a
  # specified sample size the power is computed, not planned for.
  if (!is.null(desired_power)) .dmar_check_desired_power(desired_power)
  if (!is.null(N)) {
    if (N < 4) .dmar_refuse_argument(
      "'specified_N' must be at least 4.", arg = "specified_N")
    out <- data.frame(
      term  = c("specified_N", "actual_power", "true_cor", "null_cor"),
      value = c(N, power_at(N), rho, null_cor)
    )
    class(out) <- c("dmar_ss_power", "dmar_tbl", "data.frame")
    return(out)
  }

  .dmar_check_desired_power(desired_power)

  z_alpha <- if (directional) qnorm(1 - alpha_level) else qnorm(1 - alpha_level / 2)
  z_beta  <- qnorm(desired_power)
  N_start <- ceiling(((z_alpha + z_beta) / abs(z_diff))^2 + 3)
  N_start <- max(N_start, 4)

  # Fail fast when the target cannot be reached at any realistic sample size,
  # for example when 'true_cor' sits within rounding error of 'null_cor'. Power is
  # monotone in N, so a bound that misses the target rules out every smaller
  # size. The bound also keeps the search in the range where it is meaningful:
  # a degenerate effect sends the closed-form start into the 1e24 range, where
  # N - 1 is not a distinct double and the walk-down step can never advance.
  N_max <- 1e7
  if (power_at(N_max) < desired_power) {
    .dmar_refuse(c(
      "Could not reach 'desired_power' at a sample size up to ",
      format(N_max, scientific = FALSE), "; the difference between 'true_cor' ",
      "and 'null_cor' may be too small to detect at this power."),
      "computation_failed")
  }

  N_i <- min(N_start, N_max)
  pwr <- power_at(N_i)
  # Walk up if needed (the closed form can be slightly low for two-sided)
  while (pwr < desired_power) {
    N_i <- N_i + 1
    pwr <- power_at(N_i)
    if (N_i > N_max) {
      .dmar_refuse(c(
        "The search for 'desired_power' did not converge below ",
        format(N_max, scientific = FALSE), "."), "computation_failed")
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
      .dmar_refuse(c(
        "The search for 'desired_power' did not converge within ",
        format(N_max, scientific = FALSE), " steps."), "computation_failed")
    }
  }

  out <- data.frame(
    term  = c("necessary_N", "actual_power", "true_cor", "null_cor"),
    value = c(N_i, pwr, rho, null_cor)
  )
  class(out) <- c("dmar_ss_power", "dmar_tbl", "data.frame")
  out
}
