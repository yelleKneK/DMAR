#' Sample Size or Power for an Unstandardized Contrast in a One-Way ANCOVA
#'
#' Determine the necessary per-group sample size to achieve a desired level of statistical power
#' for the test of a single planned (unstandardized) contrast on the adjusted means in a one-way
#' analysis of covariance, or, given a per-group sample size, return the realized statistical power.
#'
#' @param true_contrast The population unstandardized contrast effect on the adjusted means, \eqn{\psi = \sum c_j \mu^{(adj)}_j}
#' @param c_weights Vector of contrast weights (must sum to zero); use fractional weights so the positive weights sum to 1
#' @param true_sd Within-group population standard deviation of the response (the same \eqn{\sigma} as in a one-way ANOVA on the response)
#' @param true_cor Within-group population correlation between the response and the covariate, \eqn{\rho} in Details; must lie in (-1, 1)
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param specified_n_per_group Per-group sample size (assumed balanced); if specified, returns the realized power
#' @param directional Logical: \code{TRUE} for a one-sided test (in the same sign as \code{true_contrast}), \code{FALSE} (default) for a two-sided test
#'
#' @details
#' This function uses the standard large-sample formulation in which the ANCOVA error variance is
#' \eqn{\sigma^2_{adj} = \sigma^2 (1 - \rho^2)} (\eqn{\sigma} is \code{true_sd} and \eqn{\rho} is \code{true_cor}), the contrast \emph{t}-statistic has degrees of freedom
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
#' # Same population contrast as in the ANOVA example, with true_cor = 0.5 between
#' # outcome and covariate; ANCOVA is more efficient than ANOVA here.
#' ss_power_c_ancova(true_contrast = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
#'                   true_sd = 1, true_cor = 0.5, desired_power = 0.80)
#'
#' # Realized power for n = 30 per group
#' ss_power_c_ancova(true_contrast = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
#'                   true_sd = 1, true_cor = 0.5, specified_n_per_group = 30)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export
ss_power_c_ancova <- function(true_contrast, c_weights, true_sd, true_cor,
                              desired_power = 0.85, alpha_level = 0.05,
                              specified_n_per_group = NULL, directional = FALSE) {
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
  if (!is.numeric(true_cor) || length(true_cor) != 1L || true_cor <= -1 || true_cor >= 1) {
    .dmar_refuse_argument(
      "'true_cor' must be a single numeric value in (-1, 1).", arg = "true_cor")
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    .dmar_refuse_argument(
      "'alpha_level' must be in (0, 1).", arg = "alpha_level")
  }

  J <- length(c_weights)
  ssq_c <- sum(c_weights^2)
  sigma_adj <- true_sd * sqrt(1 - true_cor^2)

  power_at <- function(n) {
    df  <- n * J - J - 1
    if (df <= 0) return(NA_real_)
    ncp <- true_contrast / (sigma_adj * sqrt(ssq_c / n))
    # A directional test is taken in the direction of true_contrast; the
    # power is that of |ncp| (R/ss_power_internals.R).
    .power_noncentral_t(ncp, df, alpha_level, directional)
  }

  # An invalid desired_power is refused whether or not it is used: with a
  # specified sample size the power is computed, not planned for.
  if (!is.null(desired_power)) .dmar_check_desired_power(desired_power)
  if (!is.null(specified_n_per_group)) {
    if (!is.numeric(specified_n_per_group) || length(specified_n_per_group) != 1L ||
        !is.finite(specified_n_per_group) || specified_n_per_group < 2) {
      .dmar_refuse_argument(paste0(
        "'specified_n_per_group' must be a single number of at least 2: it is the ",
        "number of individuals in each group, and the error variance of ",
        "the analysis of covariance cannot be estimated from fewer. Give ",
        "the per-group sample size, or leave 'specified_n_per_group' out to have it ",
        "planned for 'desired_power'."), arg = "specified_n_per_group")
    }
    pwr <- power_at(specified_n_per_group)
    ncp <- true_contrast / (sigma_adj * sqrt(ssq_c / specified_n_per_group))
    out <- data.frame(
      term  = c("specified_n_per_group", "actual_power", "noncentral_t_parm"),
      value = c(specified_n_per_group, pwr, ncp)
    )
    class(out) <- c("dmar_ss_power", class(out))
    return(.as_dmar_tbl(out))
  }

  .dmar_check_desired_power(desired_power)

  # Evaluate the smallest admissible size before incrementing, so a target that
  # the minimum already attains returns that minimum rather than one above it.
  n_i <- 2
  pwr <- power_at(n_i)
  while (pwr < desired_power) {
    n_i <- n_i + 1
    if (n_i > 1e6) .dmar_refuse_no_sample_size(
      "ss_power_c_ancova", desired_power, size = "per-group sample size")
    pwr <- power_at(n_i)
  }
  out <- data.frame(
    term  = c("necessary_n_per_group", "actual_power", "noncentral_t_parm"),
    value = c(n_i, pwr, true_contrast / (sigma_adj * sqrt(ssq_c / n_i)))
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
