#' Sample Size or Power for an Unstandardized Contrast in a One-Way Between-Subjects ANOVA
#'
#' Determine the necessary per-group sample size to achieve a desired level of statistical power
#' for the test of a single planned (unstandardized) contrast in a one-way between-subjects
#' analysis of variance, or, given a per-group sample size, return the realized statistical power.
#'
#' @param psi The population unstandardized contrast effect, \eqn{\psi = \sum c_j \mu_j}
#' @param c_weights Vector of contrast weights (must sum to zero); use fractional weights so that the positive weights sum to 1 (e.g., \code{c(0.5, 0.5, -0.5, -0.5)})
#' @param sigma Within-group population standard deviation
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param n Per-group sample size (assumed balanced); if specified, returns the realized power
#' @param directional Logical: \code{TRUE} for a one-sided test (in the same sign as \code{psi}), \code{FALSE} (default) for a two-sided test
#'
#' @details
#' Under the alternative hypothesis the contrast \emph{t}-statistic follows a noncentral \emph{t}-distribution with
#' degrees of freedom \eqn{N - J} (where \eqn{N = n J} is the total sample size and \eqn{J} the number
#' of groups, taken as \code{length(c_weights)}) and noncentrality parameter
#' \eqn{\lambda = \psi / (\sigma \sqrt{\sum c_j^2 / n})}.
#'
#' The function searches over per-group sample sizes \eqn{n} until power first reaches
#' \code{desired_power}; when \code{n} is supplied it instead returns the realized power.
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
#' ss_power_c(psi = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5), sigma = 1,
#'            desired_power = 0.80)
#'
#' # Realized power for n = 30 per group
#' ss_power_c(psi = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5), sigma = 1, n = 30)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export
ss_power_c <- function(psi, c_weights, sigma, desired_power = 0.85,
                       alpha_level = 0.05, n = NULL, directional = FALSE) {
  if (!is.numeric(psi) || length(psi) != 1L || !is.finite(psi)) {
    stop("'psi' must be a single finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(c_weights) || length(c_weights) < 2) {
    stop("'c_weights' must be a numeric vector of length >= 2.", call. = FALSE)
  }
  if (round(sum(c_weights), 8) != 0) {
    stop("'c_weights' must sum to zero.", call. = FALSE)
  }
  if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0) {
    stop("'sigma' must be a single positive numeric value.", call. = FALSE)
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be in (0, 1).", call. = FALSE)
  }

  J <- length(c_weights)
  ssq_c <- sum(c_weights^2)

  power_at <- function(n) {
    df  <- n * J - J
    if (df <= 0) return(NA_real_)
    ncp <- psi / (sigma * sqrt(ssq_c / n))
    if (directional) {
      crit <- qt(1 - alpha_level, df = df)
      if (psi >= 0) 1 - pt(crit, df = df, ncp = ncp) else pt(-crit, df = df, ncp = ncp)
    } else {
      cu <- qt(1 - alpha_level / 2, df = df)
      cl <- qt(alpha_level / 2,     df = df)
      (1 - pt(cu, df = df, ncp = ncp)) + pt(cl, df = df, ncp = ncp)
    }
  }

  if (!is.null(n)) {
    if (n < 2) stop("'n' (per-group sample size) must be at least 2.", call. = FALSE)
    pwr <- power_at(n)
    ncp <- psi / (sigma * sqrt(ssq_c / n))
    out <- data.frame(
      term  = c("specified_n_per_group", "actual_power", "noncentral_t_parm"),
      value = c(n, pwr, ncp)
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
    term  = c("necessary_n_per_group", "actual_power", "noncentral_t_parm"),
    value = c(n_i, pwr, psi / (sigma * sqrt(ssq_c / n_i)))
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
