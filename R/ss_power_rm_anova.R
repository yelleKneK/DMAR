#' Sample Size or Power for a One-Way Repeated Measures ANOVA Omnibus \emph{F} Test
#'
#' Determine the necessary number of subjects to achieve a desired level of statistical power for
#' the omnibus \emph{F} test of the within-subjects factor in a one-way repeated measures ANOVA, or, given
#' a number of subjects, return the realized statistical power.
#'
#' @param a Number of measurement occasions (levels of the within-subjects factor)
#' @param f Cohen's \emph{f} effect size for the within-subjects factor (the population value); supply this or \code{eta_squared}, but not both
#' @param eta_squared Population eta squared (proportion of variance, on the relevant scale, accounted for by the within-subjects factor); supply this or \code{f}
#' @param rho Average correlation among the repeated measures (default 0). With \code{rho > 0}, the within-subjects test gains efficiency relative to a between-subjects analogue
#' @param epsilon Greenhouse-Geisser / Huynh-Feldt sphericity adjustment in (0, 1] (default 1, sphericity assumed). When \code{epsilon < 1}, both numerator and denominator degrees of freedom are multiplied by \code{epsilon}
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param n Number of subjects (each measured at all \code{a} occasions); if specified, returns the realized power
#'
#' @details
#' Under the alternative hypothesis with sphericity (\code{epsilon = 1}), the within-subjects \emph{F}
#' statistic follows a noncentral \emph{F} distribution with numerator df \eqn{a - 1}, denominator df
#' \eqn{(n - 1)(a - 1)}, and noncentrality parameter
#' \eqn{\lambda = n a f^2 / (1 - \rho)}, where \eqn{f} is Cohen's \emph{f} for the within-subjects
#' effect and \eqn{\rho} is the average correlation across the repeated measures (Maxwell, Delaney,
#' & Kelley, 2027). Setting \code{rho = 0} reduces to the between-subjects expression.
#'
#' When sphericity is violated, supplying \code{epsilon} (e.g., a Greenhouse-Geisser estimate)
#' rescales the test using the Muller-Barton convention: both numerator and denominator
#' degrees of freedom are multiplied by \code{epsilon}, and the noncentrality parameter is
#' likewise multiplied by \code{epsilon}. Smaller \code{epsilon} therefore reduces power and
#' increases the necessary sample size.
#'
#' @return A \code{data.frame} with rows for \code{necessary_n_subjects} (or \code{specified_n_subjects}),
#'   \code{a}, \code{effect_df}, \code{error_df}, \code{noncentrality}, and \code{actual_power}.
#'   The result carries the \code{dmar_ss_power} class, so
#'   \code{\link[generics]{tidy}} and \code{\link[generics]{glance}} summarize
#'   it in broom convention (the reported size is the number of subjects).
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing experiments and analyzing data: A model comparison perspective} (4th ed.). Routledge.
#'
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_one_way_anova}}, \code{\link{ss_power_pcm}}
#'
#' @examples
#' # 4 measurement occasions, f = 0.25, average within-subject correlation 0.5, power = .80
#' ss_power_rm_anova(a = 4, f = 0.25, rho = 0.5, desired_power = 0.80)
#'
#' # Same but with Greenhouse-Geisser epsilon = 0.75
#' ss_power_rm_anova(a = 4, f = 0.25, rho = 0.5, epsilon = 0.75, desired_power = 0.80)
#'
#' # Realized power at n = 20 subjects, a = 4
#' ss_power_rm_anova(a = 4, f = 0.25, rho = 0.5, n = 20)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export
ss_power_rm_anova <- function(a, f = NULL, eta_squared = NULL, rho = 0, epsilon = 1,
                              desired_power = 0.85, alpha_level = 0.05, n = NULL) {
  if (!is.numeric(a) || length(a) != 1L || a < 2 || a != as.integer(a)) {
    stop("'a' must be a single integer >= 2 (number of measurement occasions).", call. = FALSE)
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
  if (!is.numeric(rho) || length(rho) != 1L || rho <= -1 || rho >= 1) {
    stop("'rho' must be a single numeric value in (-1, 1).", call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0 || epsilon > 1) {
    stop("'epsilon' must be a single numeric value in (0, 1].", call. = FALSE)
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be in (0, 1).", call. = FALSE)
  }

  power_at <- function(n) {
    df_1 <- (a - 1) * epsilon
    df_2 <- (n - 1) * (a - 1) * epsilon
    if (df_2 <= 0) return(NA_real_)
    ncp <- n * a * f^2 / (1 - rho) * epsilon
    crit <- qf(1 - alpha_level, df_1, df_2)
    1 - pf(crit, df_1, df_2, ncp = ncp)
  }

  if (!is.null(n)) {
    if (n < 2) stop("'n' (number of subjects) must be at least 2.", call. = FALSE)
    pwr  <- power_at(n)
    df_1 <- (a - 1) * epsilon
    df_2 <- (n - 1) * (a - 1) * epsilon
    out <- data.frame(
      term  = c("specified_n_subjects", "a", "effect_df", "error_df", "noncentrality", "actual_power"),
      value = c(n, a, df_1, df_2, n * a * f^2 / (1 - rho) * epsilon, pwr)
    )
    class(out) <- c("dmar_ss_power", class(out))
    return(.as_dmar_tbl(out))
  }

  if (desired_power <= 0 || desired_power >= 1) {
    stop("'desired_power' must be in (0, 1).", call. = FALSE)
  }

  # Evaluate the smallest admissible size before incrementing, so a target that
  # the minimum already attains returns that minimum rather than one above it.
  n_i <- 2
  pwr <- power_at(n_i)
  while (pwr < desired_power) {
    n_i <- n_i + 1
    if (n_i > 1e6) stop("Failed to converge within reasonable sample size.", call. = FALSE)
    pwr <- power_at(n_i)
  }
  out <- data.frame(
    term  = c("necessary_n_subjects", "a", "effect_df", "error_df", "noncentrality", "actual_power"),
    value = c(n_i, a, (a - 1) * epsilon, (n_i - 1) * (a - 1) * epsilon, n_i * a * f^2 / (1 - rho) * epsilon, pwr)
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
