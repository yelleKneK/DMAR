#' Sample Size or Power for a Mixed-Effects ANOVA (Between X Within Design)
#'
#' Determine the necessary per-group sample size to achieve a desired level of statistical power
#' for one of the three \emph{F} tests in a mixed-effects ANOVA with one between-subjects factor and one
#' within-subjects factor -- between-subjects main effect, within-subjects main effect, or the
#' between x within interaction -- or, given a per-group sample size, return the realized
#' statistical power. (This design is also commonly called a split-plot factorial.)
#'
#' @param a Number of levels of the between-subjects factor (i.e., number of groups)
#' @param b Number of levels of the within-subjects factor (i.e., number of measurement occasions)
#' @param effect Which \emph{F} test to compute power for: \code{"between"}, \code{"within"}, or \code{"interaction"}
#' @param f Cohen's \emph{f} effect size for the chosen effect (the population value); supply this or \code{partial_eta_squared}, but not both
#' @param partial_eta_squared Partial eta squared for the chosen effect; supply this or \code{f}
#' @param rho Average correlation among the repeated measures within a subject (must lie in (-1, 1)). Higher \code{rho} reduces power for the between-subjects test (because subject means contain more redundant information) and increases power for the within-subjects and interaction tests
#' @param epsilon Greenhouse-Geisser / Huynh-Feldt sphericity adjustment in (0, 1] (default 1, sphericity assumed). Applied to the within-subjects and interaction tests but not the between-subjects test. Both numerator and denominator df, and the noncentrality, are multiplied by \code{epsilon} (Muller-Barton convention)
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param n Per-group (between-subjects) sample size; if specified, returns the realized power
#'
#' @details
#' This is a two-factor mixed-effects design: one between-subjects factor with \eqn{a} levels and
#' one within-subjects factor with \eqn{b} levels; \eqn{n} subjects are randomly assigned to each
#' between-subjects level and each subject is measured at all \eqn{b} within-subjects levels, for
#' \eqn{N = na} subjects total and \eqn{Nb} observations. The covariance among the \eqn{b}
#' within-subject observations is summarized by \code{rho}, the average pairwise correlation.
#'
#' The three \emph{F} tests have noncentrality parameters
#' \deqn{\lambda_{B} = N b f^2 / (1 + (b - 1) \rho)}
#' for the between-subjects test (numerator df \eqn{a - 1}, denominator df \eqn{N - a}),
#' \deqn{\lambda_{W} = N b f^2 \, \epsilon / (1 - \rho)}
#' for the within-subjects test (numerator df \eqn{(b - 1)\epsilon}, denominator df \eqn{(N - a)(b - 1)\epsilon}),
#' and the same form as \eqn{\lambda_W} for the interaction (numerator df \eqn{(a - 1)(b - 1)\epsilon},
#' same denominator df). Cohen's \emph{f} relates to partial eta squared via
#' \eqn{f = \sqrt{\eta_p^2 / (1 - \eta_p^2)}}.
#'
#' This design is the compound-symmetry (random intercept) special case of the
#' two-level linear mixed-effects model: \code{rho} is the intraclass
#' correlation and the \eqn{b} occasions are the level-1 units of a subject.
#' For two between-subjects groups the between-subjects \emph{F}(1, .) test is
#' therefore the two-level treatment \emph{t} test of
#' \code{\link{ss_power_mixed_effects}} squared, so the two planners agree on
#' that shared case.
#'
#' @return A \code{data.frame} with rows for \code{necessary_n_per_group} (or
#'   \code{specified_n_per_group}), \code{total_N}, \code{effect_df}, \code{error_df},
#'   \code{noncentrality}, and \code{actual_power}. The result carries the
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
#' ss_power_split_plot_anova(a = 2, b = 4, effect = "between", f = 0.25,
#'                     rho = 0.5, desired_power = 0.80)
#'
#' # Same design, within-subjects (occasion) main effect, f = 0.25
#' ss_power_split_plot_anova(a = 2, b = 4, effect = "within", f = 0.25,
#'                     rho = 0.5, desired_power = 0.80)
#'
#' # Same design, between x within interaction, f = 0.25
#' ss_power_split_plot_anova(a = 2, b = 4, effect = "interaction", f = 0.25,
#'                     rho = 0.5, desired_power = 0.80)
#'
#' # Realized power for n = 25 per group on the interaction test, partial eta^2 = 0.06
#' ss_power_split_plot_anova(a = 2, b = 4, effect = "interaction",
#'                     partial_eta_squared = 0.06, rho = 0.5, n = 25)
#'
#' # Greenhouse-Geisser correction with epsilon = 0.7 on the within-subjects test
#' ss_power_split_plot_anova(a = 2, b = 4, effect = "within", f = 0.25,
#'                     rho = 0.5, epsilon = 0.7, desired_power = 0.80)
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
ss_power_split_plot_anova <- function(a, b, effect, f = NULL, partial_eta_squared = NULL,
                                rho, epsilon = 1, desired_power = 0.85,
                                alpha_level = 0.05, n = NULL) {
  if (!is.numeric(a) || length(a) != 1L || a < 2 || a != as.integer(a)) {
    stop("'a' must be a single integer >= 2 (between-subjects factor levels).", call. = FALSE)
  }
  if (!is.numeric(b) || length(b) != 1L || b < 2 || b != as.integer(b)) {
    stop("'b' must be a single integer >= 2 (within-subjects factor levels).", call. = FALSE)
  }
  if (missing(effect) || !(effect %in% c("between", "within", "interaction"))) {
    stop("'effect' must be one of \"between\", \"within\", or \"interaction\".", call. = FALSE)
  }
  if (!is.null(f) && !is.null(partial_eta_squared)) {
    stop("Specify either 'f' or 'partial_eta_squared', not both.", call. = FALSE)
  }
  if (is.null(f) && is.null(partial_eta_squared)) {
    stop("Specify exactly one of 'f' or 'partial_eta_squared'.", call. = FALSE)
  }
  if (!is.null(partial_eta_squared)) {
    if (partial_eta_squared <= 0 || partial_eta_squared >= 1) {
      stop("'partial_eta_squared' must be in (0, 1).", call. = FALSE)
    }
    f <- sqrt(partial_eta_squared / (1 - partial_eta_squared))
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

  components <- function(n) {
    N <- n * a
    if (effect == "between") {
      df_1 <- a - 1
      df_2 <- N - a
      ncp  <- N * b * f^2 / (1 + (b - 1) * rho)
    } else if (effect == "within") {
      df_1 <- (b - 1) * epsilon
      df_2 <- (N - a) * (b - 1) * epsilon
      ncp  <- N * b * f^2 / (1 - rho) * epsilon
    } else {
      df_1 <- (a - 1) * (b - 1) * epsilon
      df_2 <- (N - a) * (b - 1) * epsilon
      ncp  <- N * b * f^2 / (1 - rho) * epsilon
    }
    list(df_1 = df_1, df_2 = df_2, ncp = ncp, N = N)
  }

  power_at <- function(n) {
    cc <- components(n)
    if (cc$df_2 <= 0) return(NA_real_)
    crit <- qf(1 - alpha_level, cc$df_1, cc$df_2)
    1 - pf(crit, cc$df_1, cc$df_2, ncp = cc$ncp)
  }

  if (!is.null(n)) {
    if (n < 2) stop("'n' (per-group sample size) must be at least 2.", call. = FALSE)
    cc <- components(n)
    pwr <- power_at(n)
    out <- data.frame(
      term  = c("specified_n_per_group", "total_N", "effect_df", "error_df", "noncentrality", "actual_power"),
      value = c(n, cc$N, cc$df_1, cc$df_2, cc$ncp, pwr)
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
    pwr <- power_at(n_i)
    if (n_i > 1e6) stop("Failed to converge within reasonable sample size.", call. = FALSE)
  }
  cc <- components(n_i)
  out <- data.frame(
    term  = c("necessary_n_per_group", "total_N", "effect_df", "error_df", "noncentrality", "actual_power"),
    value = c(n_i, cc$N, cc$df_1, cc$df_2, cc$ncp, pwr)
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
