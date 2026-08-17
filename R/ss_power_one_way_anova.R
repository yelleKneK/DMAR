#' Sample Size or Power for a One-Way Between-Subjects ANOVA Omnibus \emph{F} Test
#'
#' Determine the necessary total sample size to achieve a desired level of statistical power
#' for the omnibus \emph{F} test in a one-way between-subjects analysis of variance, or, given a total
#' sample size, return the realized statistical power.
#'
#' @param a Number of groups (levels of the between-subjects factor)
#' @param f Cohen's \emph{f} effect size (the population value); supply this or \code{eta_squared}, but not both
#' @param eta_squared Population eta squared (proportion of total variance accounted for by group membership); supply this or \code{f}
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param N Total sample size; if specified, the function returns the realized power (the ss_power_* family is not uniform here: \code{\link{ss_power_contrast}} takes a \emph{per-group} size)
#'
#' @details
#' Under the alternative hypothesis, the omnibus \emph{F} statistic follows a noncentral \emph{F} distribution with
#' numerator df \eqn{a - 1}, denominator df \eqn{N - a}, and noncentrality parameter
#' \eqn{\lambda = N f^2}. Cohen's \emph{f} relates to eta squared via
#' \eqn{f = \sqrt{\eta^2 / (1 - \eta^2)}}.
#'
#' The function searches over total sample sizes \eqn{N} (treating per-group \eqn{N/a} as balanced)
#' until power first reaches \code{desired_power}. When \code{N} is supplied it instead reports the
#' realized power at that \code{N}.
#'
#' @return A \code{data.frame}. When a sample size is being planned (\code{N} not
#'   supplied) the rows are \code{necessary_N}, \code{n_per_group}, \code{a},
#'   \code{noncentrality}, and \code{actual_power}; the search constructs the total as
#'   a balanced design, so \code{n_per_group} is a whole-number per-group count. When
#'   \code{N} is supplied, power is evaluated at that total \code{N} directly and the
#'   rows are \code{specified_N}, \code{a}, \code{noncentrality}, and
#'   \code{actual_power} (no \code{n_per_group} row, since balance is not assumed).
#'   The result carries the \code{dmar_ss_power} class, so
#'   \code{\link[generics]{tidy}} and \code{\link[generics]{glance}} summarize it
#'   in broom convention; the summarized sample size is the total \code{N}.
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing experiments and analyzing data: A model comparison perspective} (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_factorial_anova}}, \code{\link{ss_power_c}}, \code{\link{ss_power_sc}}, \code{\link{conf_limits_ncf}}
#'
#' @examples
#' # Three groups, f = 0.25, power = .80
#' ss_power_one_way_anova(a = 3, f = 0.25, desired_power = 0.80)
#'
#' # Same effect specified via eta squared
#' ss_power_one_way_anova(a = 3, eta_squared = 0.0588, desired_power = 0.80)
#'
#' # Realized power at N = 60 across 3 groups
#' ss_power_one_way_anova(a = 3, f = 0.25, N = 60)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export
ss_power_one_way_anova <- function(a, f = NULL, eta_squared = NULL,
                                   desired_power = 0.85, alpha_level = 0.05, N = NULL) {
  if (!is.numeric(a) || length(a) != 1L || a < 2) {
    stop("'a' must be a single integer >= 2 (number of groups).", call. = FALSE)
  }
  if (!is.null(f) && !is.null(eta_squared)) {
    stop("Specify either 'f' or 'eta_squared', not both.", call. = FALSE)
  }
  if (is.null(f) && is.null(eta_squared)) {
    stop("Specify exactly one of 'f' or 'eta_squared'.", call. = FALSE)
  }
  if (!is.null(eta_squared)) {
    if (eta_squared <= 0 || eta_squared >= 1) {
      stop("'eta_squared' must be in (0, 1).", call. = FALSE)
    }
    f <- sqrt(eta_squared / (1 - eta_squared))
  }
  if (f <= 0) stop("'f' must be positive.", call. = FALSE)
  if (alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be in (0, 1).", call. = FALSE)
  }

  power_at <- function(N) {
    df_1 <- a - 1
    df_2 <- N - a
    if (df_2 <= 0) return(NA_real_)
    ncp  <- N * f^2
    crit <- qf(1 - alpha_level, df_1, df_2)
    1 - pf(crit, df_1, df_2, ncp = ncp)
  }

  if (!is.null(N)) {
    if (N <= a) stop("'N' must exceed a (need at least a + 1 observations).", call. = FALSE)
    pwr <- power_at(N)
    out <- data.frame(
      term  = c("specified_N", "a", "noncentrality", "actual_power"),
      value = c(N, a, N * f^2, pwr)
    )
    class(out) <- c("dmar_ss_power", class(out))
    return(.as_dmar_tbl(out))
  }

  if (desired_power <= 0 || desired_power >= 1) {
    stop("'desired_power' must be in (0, 1).", call. = FALSE)
  }

  n_i <- 1
  pwr <- 0
  while (pwr < desired_power) {
    n_i <- n_i + 1
    pwr <- power_at(n_i * a)
    if (n_i > 1e6) stop("Failed to converge within reasonable sample size.", call. = FALSE)
  }
  N_i <- n_i * a
  out <- data.frame(
    term  = c("necessary_N", "n_per_group", "a", "noncentrality", "actual_power"),
    value = c(N_i, n_i, a, N_i * f^2, pwr)
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
