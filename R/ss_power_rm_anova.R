#' Sample Size or Power for a One-Way Repeated Measures ANOVA Omnibus \emph{F} Test
#'
#' Determine the necessary number of subjects (the total sample size \eqn{N}) to achieve a desired level of statistical power for
#' the omnibus \emph{F} test of the within-subjects factor in a one-way repeated measures ANOVA, or, given
#' a number of subjects, return the realized statistical power.
#'
#' @param a Number of measurement occasions (levels of the within-subjects factor)
#' @param true_f Cohen's \emph{f} effect size for the within-subjects factor (the population value); supply this or \code{true_eta_squared}, but not both
#' @param true_eta_squared Population eta squared (proportion of variance, on the relevant scale, accounted for by the within-subjects factor); supply this or \code{true_f}
#' @param true_cor Average correlation among the repeated measures (default 0). With \code{true_cor > 0}, the within-subjects test gains efficiency relative to a between-subjects analog
#' @param true_epsilon Greenhouse-Geisser / Huynh-Feldt sphericity adjustment in (0, 1] (default 1, sphericity assumed). When \code{true_epsilon < 1}, both numerator and denominator degrees of freedom are multiplied by \code{true_epsilon}
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param specified_N Total sample size: the number of subjects, each measured at all \code{a} occasions; if specified, the function returns the realized power
#'
#' @details
#' Under the alternative hypothesis with sphericity (\code{true_epsilon = 1}), the within-subjects \emph{F}
#' statistic follows a noncentral \emph{F} distribution with numerator df \eqn{a - 1}, denominator df
#' \eqn{(N - 1)(a - 1)}, and noncentrality parameter
#' \eqn{\lambda = N a f^2 / (1 - \rho)}, where \eqn{N} is the number of subjects, \eqn{f} is Cohen's \emph{f} (\code{true_f}) for the within-subjects
#' effect and \eqn{\rho} is the average correlation across the repeated measures
#' (\code{true_cor}; Maxwell, Delaney,
#' & Kelley, 2027). Setting \code{true_cor = 0} reduces to the between-subjects expression.
#'
#' When sphericity is violated, supplying \code{true_epsilon} (e.g., a Greenhouse-Geisser estimate)
#' rescales the test using the Muller-Barton convention: both numerator and denominator
#' degrees of freedom are multiplied by \code{true_epsilon}, and the noncentrality parameter is
#' likewise multiplied by \code{true_epsilon}. Smaller \code{true_epsilon} therefore reduces power and
#' increases the necessary sample size.
#'
#' @return A \code{data.frame} with rows for \code{necessary_N} (or \code{specified_N}), the number of subjects,
#'   \code{a}, \code{df_effect}, \code{df_error}, \code{noncentral_f_parm}, and \code{actual_power}.
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
#' ss_power_rm_anova(a = 4, true_f = 0.25, true_cor = 0.5, desired_power = 0.80)
#'
#' # Same but with Greenhouse-Geisser epsilon = 0.75
#' ss_power_rm_anova(a = 4, true_f = 0.25, true_cor = 0.5, true_epsilon = 0.75, desired_power = 0.80)
#'
#' # Realized power at N = 20 subjects, a = 4
#' ss_power_rm_anova(a = 4, true_f = 0.25, true_cor = 0.5, specified_N = 20)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export
ss_power_rm_anova <- function(a, true_f = NULL, true_eta_squared = NULL, true_cor = 0, true_epsilon = 1,
                              desired_power = 0.85, alpha_level = 0.05, specified_N = NULL) {
  .dmar_check_alpha_level(alpha_level)
  .dmar_check_numbers("ss_power_rm_anova", "true_f")
  if (!is.numeric(a) || length(a) != 1L || a < 2 || a != as.integer(a)) {
    .dmar_refuse_argument(
      "'a' must be a single integer >= 2 (number of measurement occasions).",
      arg = "a")
  }
  if (!is.null(true_f) && !is.null(true_eta_squared)) {
    .dmar_refuse_argument(
      "Specify either 'true_f' or 'true_eta_squared', not both.",
      arg = c("true_f", "true_eta_squared"))
  }
  if (is.null(true_f) && is.null(true_eta_squared)) {
    .dmar_refuse_argument(
      "Specify exactly one of 'true_f' or 'true_eta_squared'.",
      arg = c("true_f", "true_eta_squared"))
  }
  if (!is.null(true_eta_squared)) {
    if (true_eta_squared <= 0 || true_eta_squared >= 1) {
      .dmar_refuse_argument(
        "'true_eta_squared' must be in (0, 1).", arg = "true_eta_squared")
    }
    true_f <- sqrt(true_eta_squared / (1 - true_eta_squared))
  }
  if (true_f <= 0) .dmar_refuse_argument(
    "'true_f' must be positive.", arg = "true_f")
  if (!is.numeric(true_cor) || length(true_cor) != 1L || true_cor <= -1 || true_cor >= 1) {
    .dmar_refuse_argument(
      "'true_cor' must be a single numeric value in (-1, 1).", arg = "true_cor")
  }
  if (!is.numeric(true_epsilon) || length(true_epsilon) != 1L || true_epsilon <= 0 || true_epsilon > 1) {
    .dmar_refuse_argument(
      "'true_epsilon' must be a single numeric value in (0, 1].",
      arg = "true_epsilon")
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    .dmar_refuse_argument(
      "'alpha_level' must be in (0, 1).", arg = "alpha_level")
  }

  power_at <- function(N) {
    df_1 <- (a - 1) * true_epsilon
    df_2 <- (N - 1) * (a - 1) * true_epsilon
    if (df_2 <= 0) return(NA_real_)
    ncp <- N * a * true_f^2 / (1 - true_cor) * true_epsilon
    crit <- qf(1 - alpha_level, df_1, df_2)
    1 - pf(crit, df_1, df_2, ncp = ncp)
  }

  # An invalid desired_power is refused whether or not it is used: with a
  # specified sample size the power is computed, not planned for.
  if (!is.null(desired_power)) .dmar_check_desired_power(desired_power)
  if (!is.null(specified_N)) {
    if (specified_N < 2) .dmar_refuse_argument(
      "'specified_N' (number of subjects) must be at least 2.",
      arg = "specified_N")
    pwr  <- power_at(specified_N)
    df_1 <- (a - 1) * true_epsilon
    df_2 <- (specified_N - 1) * (a - 1) * true_epsilon
    out <- data.frame(
      term  = c("specified_N", "a", "df_effect", "df_error", "noncentral_f_parm", "actual_power"),
      value = c(specified_N, a, df_1, df_2, specified_N * a * true_f^2 / (1 - true_cor) * true_epsilon, pwr)
    )
    class(out) <- c("dmar_ss_power", class(out))
    return(.as_dmar_tbl(out))
  }

  .dmar_check_desired_power(desired_power)

  # Evaluate the smallest admissible size before incrementing, so a target that
  # the minimum already attains returns that minimum rather than one above it.
  N_i <- 2
  pwr <- power_at(N_i)
  while (pwr < desired_power) {
    N_i <- N_i + 1
    if (N_i > 1e6) .dmar_refuse_no_sample_size(
      "ss_power_rm_anova", desired_power, size = "sample size")
    pwr <- power_at(N_i)
  }
  out <- data.frame(
    term  = c("necessary_N", "a", "df_effect", "df_error", "noncentral_f_parm", "actual_power"),
    value = c(N_i, a, (a - 1) * true_epsilon, (N_i - 1) * (a - 1) * true_epsilon, N_i * a * true_f^2 / (1 - true_cor) * true_epsilon, pwr)
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
