#' Sample Size or Power for a Standardized Contrast in a One-Way Between-Subjects ANOVA
#'
#' Determine the necessary per-group sample size to achieve a desired level of statistical power
#' for the test of a single planned standardized contrast in a one-way between-subjects analysis
#' of variance, or, given a per-group sample size, return the realized statistical power.
#'
#' @param psi_standardized The population standardized contrast effect, \eqn{\psi / \sigma}, where \eqn{\sigma} is the within-group standard deviation
#' @param c_weights Vector of contrast weights (must sum to zero); use fractional weights so that the positive weights sum to 1 (e.g., \code{c(0.5, 0.5, -0.5, -0.5)})
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param n Per-group sample size (assumed balanced); if specified, returns the realized power
#' @param directional Logical: \code{TRUE} for a one-sided test (in the same sign as \code{psi_standardized}), \code{FALSE} (default) for a two-sided test
#'
#' @details
#' Under the alternative hypothesis the contrast \emph{t}-statistic follows a noncentral \emph{t}-distribution with
#' degrees of freedom \eqn{N - J} (\eqn{J = }\code{length(c_weights)}) and noncentrality parameter
#' \eqn{\lambda = \psi^* / \sqrt{\sum c_j^2 / n}}, where \eqn{\psi^*} is the standardized contrast.
#'
#' The function searches over per-group sample sizes \eqn{n} until power first reaches
#' \code{desired_power}; when \code{n} is supplied it returns the realized power.
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
#' @seealso \code{\link{ss_power_c}}, \code{\link{ss_power_one_way_anova}}, \code{\link{ci_sc}}, \code{\link{ss_aipe_sc}}
#'
#' @examples
#' # Power for a standardized contrast of 0.5 across 4 groups,
#' # contrast (G1 + G2)/2 vs (G3 + G4)/2, desired power = .80
#' ss_power_sc(psi_standardized = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
#'             desired_power = 0.80)
#'
#' # Realized power at n = 30 per group
#' ss_power_sc(psi_standardized = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5), n = 30)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export
ss_power_sc <- function(psi_standardized, c_weights, desired_power = 0.85,
                        alpha_level = 0.05, n = NULL, directional = FALSE) {
  if (!is.numeric(psi_standardized) || length(psi_standardized) != 1L || !is.finite(psi_standardized)) {
    stop("'psi_standardized' must be a single finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(c_weights) || length(c_weights) < 2) {
    stop("'c_weights' must be a numeric vector of length >= 2.", call. = FALSE)
  }
  if (round(sum(c_weights), 8) != 0) {
    stop("'c_weights' must sum to zero.", call. = FALSE)
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be in (0, 1).", call. = FALSE)
  }

  J <- length(c_weights)
  ssq_c <- sum(c_weights^2)

  power_at <- function(n) {
    df  <- n * J - J
    if (df <= 0) return(NA_real_)
    ncp <- psi_standardized / sqrt(ssq_c / n)
    if (directional) {
      crit <- qt(1 - alpha_level, df = df)
      if (psi_standardized >= 0) 1 - pt(crit, df = df, ncp = ncp) else pt(-crit, df = df, ncp = ncp)
    } else {
      cu <- qt(1 - alpha_level / 2, df = df)
      cl <- qt(alpha_level / 2,     df = df)
      (1 - pt(cu, df = df, ncp = ncp)) + pt(cl, df = df, ncp = ncp)
    }
  }

  if (!is.null(n)) {
    if (n < 2) stop("'n' (per-group sample size) must be at least 2.", call. = FALSE)
    pwr <- power_at(n)
    ncp <- psi_standardized / sqrt(ssq_c / n)
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
    value = c(n_i, pwr, psi_standardized / sqrt(ssq_c / n_i))
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
