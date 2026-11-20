#' Sample Size or Power for a Mixed-Effects ANOVA (Between X Within Design)
#'
#' Determine the necessary per-group sample size to achieve a desired level of statistical power
#' for one of the three \emph{F} tests in a mixed-effects ANOVA with one between-subjects factor and one
#' within-subjects factor (the between-subjects main effect, the within-subjects main effect, or the
#' between x within interaction), or, given a per-group sample size, return the realized
#' statistical power. (This design is also commonly called a split-plot factorial.)
#'
#' @param a Number of levels of the between-subjects factor (i.e., number of groups)
#' @param b Number of levels of the within-subjects factor (i.e., number of measurement occasions)
#' @param effect Which \emph{F} test to compute power for: \code{"between"}, \code{"within"}, or \code{"interaction"}
#' @param true_f Cohen's \emph{f} effect size for the chosen effect (the population value); supply this or \code{true_eta_squared_partial}, but not both
#' @param true_eta_squared_partial Population partial eta squared for the chosen effect; supply this or \code{true_f}
#' @param true_cor Average correlation among the repeated measures within a subject (must lie in (-1, 1)). Higher \code{true_cor} reduces power for the between-subjects test (because subject means contain more redundant information) and increases power for the within-subjects and interaction tests
#' @param true_epsilon Greenhouse-Geisser / Huynh-Feldt sphericity adjustment in (0, 1] (default 1, sphericity assumed). Applied to the within-subjects and interaction tests but not the between-subjects test. Both numerator and denominator df, and the noncentrality, are multiplied by \code{true_epsilon} (Muller-Barton convention)
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param specified_n_per_group Per-group (between-subjects) sample size; if specified, returns the realized power
#'
#' @details
#' This is a two-factor mixed-effects design: one between-subjects factor with \eqn{a} levels and
#' one within-subjects factor with \eqn{b} levels; \eqn{n} subjects are randomly assigned to each
#' between-subjects level and each subject is measured at all \eqn{b} within-subjects levels, for
#' \eqn{N = na} subjects total and \eqn{Nb} observations. The covariance among the \eqn{b}
#' within-subject observations is summarized by \code{true_cor}, the average pairwise correlation.
#'
#' The three \emph{F} tests have noncentrality parameters
#' \deqn{\lambda_{B} = N b f^2 / (1 + (b - 1) \rho)}
#' for the between-subjects test (numerator df \eqn{a - 1}, denominator df \eqn{N - a}),
#' \deqn{\lambda_{W} = N b f^2 \, \epsilon / (1 - \rho)}
#' for the within-subjects test (numerator df \eqn{(b - 1)\epsilon}, denominator df \eqn{(N - a)(b - 1)\epsilon}),
#' and the same form as \eqn{\lambda_W} for the interaction (numerator df \eqn{(a - 1)(b - 1)\epsilon},
#' same denominator df). Cohen's \emph{f} relates to partial eta squared via
#' \eqn{f = \sqrt{\eta_p^2 / (1 - \eta_p^2)}}. In these formulas \eqn{f} is
#' \code{true_f} (or the value implied by \code{true_eta_squared_partial},
#' which is \eqn{\eta_p^2}), \eqn{\rho} is \code{true_cor}, \eqn{\epsilon}
#' is \code{true_epsilon}, and \eqn{n} is the per-group size
#' (\code{specified_n_per_group}, or the size the function plans).
#'
#' This design is the compound-symmetry (random intercept) special case of the
#' two-level linear mixed-effects model: \code{true_cor} is the intraclass
#' correlation and the \eqn{b} occasions are the observations within a
#' subject (i.e., the Level 1 units).
#' For two between-subjects groups the between-subjects \emph{F}(1, .) test is
#' therefore the two-level treatment \emph{t} test of
#' \code{\link{ss_power_mixed_effects}} squared, so the two planners agree on
#' that shared case.
#'
#' @return A \code{data.frame} with rows for \code{necessary_n_per_group} (or
#'   \code{specified_n_per_group}), \code{total_N}, \code{df_effect}, \code{df_error},
#'   \code{noncentral_f_parm}, and \code{actual_power}. The result carries the
#'   \code{dmar_ss_power} class, so \code{\link[generics]{tidy}} and
#'   \code{\link[generics]{glance}} summarize it in broom convention (the
#'   reported size is the per-group count).
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing experiments and analyzing data: A model comparison perspective} (4th ed.). Routledge.
#'
#' Muller, K. E., & Barton, C. N. (1989). Approximate power for repeated measures ANOVA lacking sphericity. \emph{Journal of the American Statistical Association, 84}, 549--555.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_one_way_anova}}, \code{\link{ss_power_factorial_anova}}, \code{\link{ss_power_rm_anova}}, \code{\link{ss_power_mixed_effects}}
#'
#' @examples
#' # 2 groups, 4 occasions, between-subjects effect, f = 0.25,
#' # average within-subject correlation 0.5, power = .80
#' ss_power_split_plot_anova(a = 2, b = 4, effect = "between", true_f = 0.25,
#'                     true_cor = 0.5, desired_power = 0.80)
#'
#' # Same design, within-subjects (occasion) main effect, f = 0.25
#' ss_power_split_plot_anova(a = 2, b = 4, effect = "within", true_f = 0.25,
#'                     true_cor = 0.5, desired_power = 0.80)
#'
#' # Same design, between x within interaction, f = 0.25
#' ss_power_split_plot_anova(a = 2, b = 4, effect = "interaction", true_f = 0.25,
#'                     true_cor = 0.5, desired_power = 0.80)
#'
#' # Realized power for n = 25 per group on the interaction test, partial eta^2 = 0.06
#' ss_power_split_plot_anova(a = 2, b = 4, effect = "interaction",
#'                     true_eta_squared_partial = 0.06, true_cor = 0.5, specified_n_per_group = 25)
#'
#' # Greenhouse-Geisser correction with epsilon = 0.7 on the within-subjects test
#' ss_power_split_plot_anova(a = 2, b = 4, effect = "within", true_f = 0.25,
#'                     true_cor = 0.5, true_epsilon = 0.7, desired_power = 0.80)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#' @family mixed models
#'
#' @export
ss_power_split_plot_anova <- function(a, b, effect, true_f = NULL, true_eta_squared_partial = NULL,
                                true_cor, true_epsilon = 1, desired_power = 0.85,
                                alpha_level = 0.05, specified_n_per_group = NULL) {
  .dmar_check_alpha_level(alpha_level)
  .dmar_check_numbers("ss_power_split_plot_anova", "true_f")
  if (!is.numeric(a) || length(a) != 1L || a < 2 || a != as.integer(a)) {
    .dmar_refuse_argument(
      "'a' must be a single integer >= 2 (between-subjects factor levels).",
      arg = "a")
  }
  if (!is.numeric(b) || length(b) != 1L || b < 2 || b != as.integer(b)) {
    .dmar_refuse_argument(
      "'b' must be a single integer >= 2 (within-subjects factor levels).",
      arg = "b")
  }
  if (missing(effect) || !(effect %in% c("between", "within", "interaction"))) {
    .dmar_refuse_argument(
      "'effect' must be one of \"between\", \"within\", or \"interaction\".",
      arg = "effect")
  }
  if (!is.null(true_f) && !is.null(true_eta_squared_partial)) {
    .dmar_refuse_argument(
      "Specify either 'true_f' or 'true_eta_squared_partial', not both.",
      arg = c("true_f", "true_eta_squared_partial"))
  }
  if (is.null(true_f) && is.null(true_eta_squared_partial)) {
    .dmar_refuse_argument(
      "Specify exactly one of 'true_f' or 'true_eta_squared_partial'.",
      arg = c("true_f", "true_eta_squared_partial"))
  }
  if (!is.null(true_eta_squared_partial)) {
    if (true_eta_squared_partial <= 0 || true_eta_squared_partial >= 1) {
      .dmar_refuse_argument(
        "'true_eta_squared_partial' must be in (0, 1).",
        arg = "true_eta_squared_partial")
    }
    true_f <- sqrt(true_eta_squared_partial / (1 - true_eta_squared_partial))
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

  components <- function(specified_n_per_group) {
    N <- specified_n_per_group * a
    if (effect == "between") {
      df_1 <- a - 1
      df_2 <- N - a
      ncp  <- N * b * true_f^2 / (1 + (b - 1) * true_cor)
    } else if (effect == "within") {
      df_1 <- (b - 1) * true_epsilon
      df_2 <- (N - a) * (b - 1) * true_epsilon
      ncp  <- N * b * true_f^2 / (1 - true_cor) * true_epsilon
    } else {
      df_1 <- (a - 1) * (b - 1) * true_epsilon
      df_2 <- (N - a) * (b - 1) * true_epsilon
      ncp  <- N * b * true_f^2 / (1 - true_cor) * true_epsilon
    }
    list(df_1 = df_1, df_2 = df_2, ncp = ncp, N = N)
  }

  power_at <- function(specified_n_per_group) {
    cc <- components(specified_n_per_group)
    if (cc$df_2 <= 0) return(NA_real_)
    crit <- qf(1 - alpha_level, cc$df_1, cc$df_2)
    1 - pf(crit, cc$df_1, cc$df_2, ncp = cc$ncp)
  }

  # An invalid desired_power is refused whether or not it is used: with a
  # specified sample size the power is computed, not planned for.
  if (!is.null(desired_power)) .dmar_check_desired_power(desired_power)
  if (!is.null(specified_n_per_group)) {
    if (specified_n_per_group < 2) .dmar_refuse_argument(
      "'specified_n_per_group' (per-group sample size) must be at least 2.",
      arg = "specified_n_per_group")
    cc <- components(specified_n_per_group)
    pwr <- power_at(specified_n_per_group)
    out <- data.frame(
      term  = c("specified_n_per_group", "total_N", "df_effect", "df_error", "noncentral_f_parm", "actual_power"),
      value = c(specified_n_per_group, cc$N, cc$df_1, cc$df_2, cc$ncp, pwr)
    )
    class(out) <- c("dmar_ss_power", class(out))
    return(.as_dmar_tbl(out))
  }

  .dmar_check_desired_power(desired_power)

  n_i <- 1
  pwr <- 0
  while (pwr < desired_power) {
    n_i <- n_i + 1
    pwr <- power_at(n_i)
    if (n_i > 1e6) .dmar_refuse_no_sample_size(
      "ss_power_split_plot_anova", desired_power, size = "per-group sample size")
  }
  cc <- components(n_i)
  out <- data.frame(
    term  = c("necessary_n_per_group", "total_N", "df_effect", "df_error", "noncentral_f_parm", "actual_power"),
    value = c(n_i, cc$N, cc$df_1, cc$df_2, cc$ncp, pwr)
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
