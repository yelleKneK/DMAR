#' Sample Size or Power for an Unstandardized Contrast in a One-Way Between-Subjects ANOVA
#'
#' Determine the necessary per-group sample size to achieve a desired level of statistical power
#' for the test of a single planned (unstandardized) contrast in a one-way between-subjects
#' analysis of variance, or, given a per-group sample size, return the realized statistical power.
#'
#' @param true_contrast The population unstandardized contrast effect, \eqn{\psi = \sum c_j \mu_j}
#' @param c_weights Vector of contrast weights (must sum to zero); use fractional weights so that the positive weights sum to 1 (e.g., \code{c(0.5, 0.5, -0.5, -0.5)})
#' @param true_sd Within-group population standard deviation, \eqn{\sigma} in Details
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param specified_n_per_group Per-group sample size (assumed balanced); if specified, returns the realized power
#' @param directional Logical: \code{TRUE} for a one-sided test (in the same sign as \code{true_contrast}), \code{FALSE} (default) for a two-sided test
#'
#' @details
#' Under the alternative hypothesis the contrast \emph{t}-statistic follows a noncentral \emph{t}-distribution with
#' degrees of freedom \eqn{N - J} (where \eqn{N = n J} is the total sample size and \eqn{J} the number
#' of groups, taken as \code{length(c_weights)}) and noncentrality parameter
#' \eqn{\lambda = \psi / (\sigma \sqrt{\sum c_j^2 / n})}.
#'
#' The function searches over per-group sample sizes \eqn{n} until power first reaches
#' \code{desired_power}; when \code{specified_n_per_group} is supplied it instead returns the realized power.
#'
#' @return A \code{data.frame} with rows for \code{necessary_n_per_group} (or \code{specified_n_per_group}),
#'   \code{actual_power}, and \code{noncentral_t_parm}. The result carries the
#'   \code{dmar_ss_power} class, so \code{\link[generics]{tidy}} and
#'   \code{\link[generics]{glance}} summarize it in broom convention.
#'
#' @references
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
#' ANCOVA and ANOVA contrasts: Sample size planning via narrow
#' confidence intervals.
#' \emph{British Journal of Mathematical and Statistical Psychology, 65},
#' 350--370. \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing experiments and analyzing data: A model comparison perspective} (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_sc}}, \code{\link{ss_power_one_way_anova}}, \code{\link{ss_power_c_ancova}}, \code{\link{ci_c}}, \code{\link{ss_aipe_c}}
#'
#' @examples
#' # Power for the contrast (Group 1 + Group 2) / 2 vs (Group 3 + Group 4) / 2
#' # with population contrast = 0.5, within-group sigma = 1, desired power = .80
#' ss_power_c(true_contrast = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5), true_sd = 1,
#'            desired_power = 0.80)
#'
#' # Realized power for n = 30 per group
#' ss_power_c(true_contrast = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5), true_sd = 1,
#'            specified_n_per_group = 30)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export
ss_power_c <- function(true_contrast, c_weights, true_sd, desired_power = 0.85,
                       alpha_level = 0.05, specified_n_per_group = NULL, directional = FALSE) {
  .dmar_check_alpha_level(alpha_level)
  if (!is.numeric(true_contrast) || length(true_contrast) != 1L || !is.finite(true_contrast)) {
    .dmar_refuse_argument(
      "'true_contrast' must be a single finite numeric value.", arg = "true_contrast")
  }
  if (!is.numeric(c_weights) || length(c_weights) < 2) {
    .dmar_refuse_argument(
      "'c_weights' must be a numeric vector of length >= 2.", arg = "c_weights")
  }
  if (round(sum(c_weights), 8) != 0) {
    .dmar_refuse_argument("'c_weights' must sum to zero.", arg = "c_weights")
  }
  # All-zero weights define no contrast; they gave a power of 1 with an
  # infinite noncentrality parameter.
  .dmar_check_contrast_nonzero(c_weights)
  if (!is.numeric(true_sd) || length(true_sd) != 1L || true_sd <= 0) {
    .dmar_refuse_argument(
      "'true_sd' must be a single positive numeric value.",
      arg = "true_sd")
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    .dmar_refuse_argument(
      "'alpha_level' must be in (0, 1).", arg = "alpha_level")
  }

  J <- length(c_weights)
  ssq_c <- sum(c_weights^2)

  power_at <- function(specified_n_per_group) {
    df  <- specified_n_per_group * J - J
    if (df <= 0) return(NA_real_)
    ncp <- true_contrast / (true_sd * sqrt(ssq_c / specified_n_per_group))
    # A directional test is taken in the direction of true_contrast; the
    # power is that of |ncp| (R/ss_power_internals.R).
    .power_noncentral_t(ncp, df, alpha_level, directional)
  }

  # An invalid desired_power is refused whether or not it is used: with a
  # specified sample size the power is computed, not planned for.
  if (!is.null(desired_power)) .dmar_check_desired_power(desired_power)
  if (!is.null(specified_n_per_group)) {
    if (specified_n_per_group < 2) .dmar_refuse_argument(
      "'specified_n_per_group' (per-group sample size) must be at least 2.",
      arg = "specified_n_per_group")
    pwr <- power_at(specified_n_per_group)
    ncp <- true_contrast / (true_sd * sqrt(ssq_c / specified_n_per_group))
    out <- data.frame(
      term  = c("specified_n_per_group", "actual_power", "noncentral_t_parm"),
      value = c(specified_n_per_group, pwr, ncp)
    )
    class(out) <- c("dmar_ss_power", class(out))
    return(.as_dmar_tbl(out))
  }

  .dmar_check_desired_power(desired_power)
  if (true_contrast == 0) .dmar_refuse_zero_effect(
    "ss_power_c", "true_contrast", "a population contrast")

  # Evaluate the smallest admissible size before incrementing, so a target that
  # the minimum already attains returns that minimum rather than one above it.
  n_i <- 2
  pwr <- power_at(n_i)
  while (pwr < desired_power) {
    n_i <- n_i + 1
    if (n_i > 1e6) .dmar_refuse_no_sample_size(
      "ss_power_c", desired_power, size = "per-group sample size")
    pwr <- power_at(n_i)
  }
  out <- data.frame(
    term  = c("necessary_n_per_group", "actual_power", "noncentral_t_parm"),
    value = c(n_i, pwr, true_contrast / (true_sd * sqrt(ssq_c / n_i)))
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
