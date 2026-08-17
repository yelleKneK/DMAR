#' Sample Size or Power for an Unstandardized Contrast in a One-Way ANCOVA
#'
#' Determine the necessary per-group sample size to achieve a desired level of statistical power
#' for the test of a single planned (unstandardized) contrast on the adjusted means in a one-way
#' analysis of covariance, or, given a per-group sample size, return the realized statistical power.
#'
#' @param psi The population unstandardized contrast effect on the adjusted means, \eqn{\psi = \sum c_j \mu^{(adj)}_j}
#' @param c_weights Vector of contrast weights (must sum to zero); use fractional weights so the positive weights sum to 1
#' @param sigma Within-group population standard deviation of the response (the same \eqn{\sigma} as in a one-way ANOVA on the response)
#' @param rho Within-group population correlation between the response and the covariate; must lie in (-1, 1)
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param n Per-group sample size (assumed balanced); if specified, returns the realized power
#' @param directional Logical: \code{TRUE} for a one-sided test (in the same sign as \code{psi}), \code{FALSE} (default) for a two-sided test
#'
#' @details
#' This function uses the standard large-sample formulation in which the ANCOVA error variance is
#' \eqn{\sigma^2_{adj} = \sigma^2 (1 - \rho^2)}, the contrast \emph{t}-statistic has degrees of freedom
#' \eqn{N - J - 1} (one less than the corresponding ANOVA contrast because of the covariate), and the
#' noncentrality parameter is
#' \eqn{\lambda = \psi / (\sigma \sqrt{1 - \rho^2} \sqrt{\sum c_j^2 / n})}. This assumes the
#' covariate means are equal across groups (the typical assumption under random assignment); for
#' designs with substantial group differences in the covariate, the small-sample correction
#' \eqn{1 + (\bar X_{j} - \bar X_{\cdot})^2 / SS^{(within)}_X} would slightly inflate the standard
#' error and reduce power, an effect that is negligible for moderate or large \eqn{n}.
#'
#' @return A \code{data.frame} with rows for \code{necessary_n_per_group} (or
#'   \code{specified_n_per_group}), \code{actual_power}, and \code{noncentral_t_parm}.
#'   The result carries the \code{dmar_ss_power} class, so
#'   \code{\link[generics]{tidy}} and \code{\link[generics]{glance}} summarize it
#'   in broom convention.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing experiments and analyzing data: A model comparison perspective} (4th ed.). Routledge.
#'
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence intervals. \emph{British Journal of Mathematical and Statistical Psychology, 65}, 350--370.
#'   \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_c}}, \code{\link{ci_c_ancova}}, \code{\link{ss_aipe_c_ancova}}
#'
#' @examples
#' # Same population contrast as in the ANOVA example, with rho = 0.5 between
#' # outcome and covariate; ANCOVA is more efficient than ANOVA here.
#' ss_power_c_ancova(psi = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
#'                   sigma = 1, rho = 0.5, desired_power = 0.80)
#'
#' # Realized power for n = 30 per group
#' ss_power_c_ancova(psi = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
#'                   sigma = 1, rho = 0.5, n = 30)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export
ss_power_c_ancova <- function(psi, c_weights, sigma, rho,
                              desired_power = 0.85, alpha_level = 0.05,
                              n = NULL, directional = FALSE) {
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
  if (!is.numeric(rho) || length(rho) != 1L || rho <= -1 || rho >= 1) {
    stop("'rho' must be a single numeric value in (-1, 1).", call. = FALSE)
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be in (0, 1).", call. = FALSE)
  }

  J <- length(c_weights)
  ssq_c <- sum(c_weights^2)
  sigma_adj <- sigma * sqrt(1 - rho^2)

  power_at <- function(n) {
    df  <- n * J - J - 1
    if (df <= 0) return(NA_real_)
    ncp <- psi / (sigma_adj * sqrt(ssq_c / n))
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
    ncp <- psi / (sigma_adj * sqrt(ssq_c / n))
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
    value = c(n_i, pwr, psi / (sigma_adj * sqrt(ssq_c / n_i)))
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
