#' Sample Size or Power for a Treatment Effect in a Two-Level Mixed-Effects Model
#'
#' Determine the necessary number of level-2 units per arm to achieve a desired level of
#' statistical power for a treatment-versus-control comparison in a two-level mixed-effects model
#' with a random intercept (e.g., individuals nested within clusters in a cluster-randomized trial,
#' or repeated measurements nested within subjects in a person-randomized longitudinal study).
#' Alternatively, given a number of level-2 units per arm, return the realized statistical power.
#'
#' @param d Standardized treatment effect, defined as the population mean difference divided by the population standard deviation of the level-1 outcome
#' @param n Number of level-1 units per level-2 unit (e.g., individuals per cluster, or measurements per subject); assumed equal across level-2 units
#' @param rho Intra-class correlation (the proportion of total outcome variance attributable to differences between level-2 units); must be in [0, 1)
#' @param J Number of level-2 units per arm (i.e., \code{J} treatment clusters and \code{J} control clusters); if specified, returns the realized power
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param directional Logical: \code{TRUE} for a one-sided test (in the same sign as \code{d}), \code{FALSE} (default) for a two-sided test
#'
#' @details
#' This function computes power for the fixed treatment effect at the higher level of a two-level
#' mixed-effects model with random intercept,
#' \deqn{y_{ij} = \beta_0 + \beta_1 T_j + u_j + \epsilon_{ij},}
#' where \eqn{u_j \sim N(0, \sigma_u^2)} is the level-2 random intercept and
#' \eqn{\epsilon_{ij} \sim N(0, \sigma_e^2)} is the level-1 residual. The treatment indicator
#' \eqn{T_j} varies between level-2 units (i.e., entire clusters or entire subjects are assigned
#' to treatment or control). The intra-class correlation is \eqn{\rho = \sigma_u^2 / (\sigma_u^2
#' + \sigma_e^2)} and the total outcome variance is \eqn{\sigma_y^2 = \sigma_u^2 + \sigma_e^2}.
#'
#' The standard error of the estimated treatment effect is
#' \eqn{SE(\hat\beta_1) = \sigma_y \sqrt{2 (1 + (n - 1)\rho) / (J n)}}, giving a noncentrality
#' parameter of
#' \deqn{\lambda = d \sqrt{J n / (2 (1 + (n - 1)\rho))}}
#' under a two-sample \emph{t}-test with \eqn{2J - 2} degrees of freedom.
#'
#' The factor \eqn{1 + (n - 1)\rho} is the design effect: as the within-cluster correlation grows,
#' the effective information per level-2 unit shrinks, so more level-2 units are needed for a given
#' level of power.
#'
#' @return A \code{data.frame} with rows for \code{necessary_J_per_arm} (or \code{specified_J_per_arm}),
#'   \code{total_N}, \code{noncentrality}, and \code{actual_power}. The result
#'   carries the \code{dmar_ss_power} class, so \code{\link[generics]{tidy}} and
#'   \code{\link[generics]{glance}} summarize it in broom convention (the
#'   reported size is the number of clusters per arm).
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing experiments and analyzing data: A model comparison perspective} (4th ed.). Routledge.
#'
#' Raudenbush, S. W. (1997). Statistical analysis and optimal design for cluster randomized trials. \emph{Psychological Methods, 2}, 173--185.
#'   \doi{10.1037/1082-989X.2.2.173}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_split_plot_anova}}, \code{\link{ss_power_smd}}, \code{\link{ss_power_pcm}}
#'
#' @examples
#' # 30 individuals per cluster, ICC = 0.05, standardized effect d = 0.30, power = .80
#' ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.05, desired_power = 0.80)
#'
#' # Same effect but with much higher ICC (e.g., schools or therapists)
#' ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.20, desired_power = 0.80)
#'
#' # Realized power with 25 level-2 units per arm
#' ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.05, J = 25)
#'
#' # directional test
#' ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.05, desired_power = 0.80,
#'                          directional = TRUE)
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
ss_power_mixed_effects <- function(d, n, rho, J = NULL,
                                     desired_power = 0.85, alpha_level = 0.05,
                                     directional = FALSE) {
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d == 0) {
    stop("'d' must be a single finite nonzero numeric value.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || n < 2 || n != as.integer(n)) {
    stop("'n' must be a single integer >= 2 (level-1 units per level-2 unit).", call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1L || rho < 0 || rho >= 1) {
    stop("'rho' must be a single numeric value in [0, 1).", call. = FALSE)
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be in (0, 1).", call. = FALSE)
  }

  power_at <- function(J) {
    df  <- 2 * J - 2
    if (df <= 0) return(NA_real_)
    ncp <- d * sqrt(J * n / (2 * (1 + (n - 1) * rho)))
    if (directional) {
      crit <- qt(1 - alpha_level, df = df)
      if (d >= 0) 1 - pt(crit, df = df, ncp = ncp) else pt(-crit, df = df, ncp = ncp)
    } else {
      cu <- qt(1 - alpha_level / 2, df = df)
      cl <- qt(alpha_level / 2,     df = df)
      (1 - pt(cu, df = df, ncp = ncp)) + pt(cl, df = df, ncp = ncp)
    }
  }

  if (!is.null(J)) {
    if (!is.numeric(J) || length(J) != 1L || J < 2 || J != as.integer(J)) {
      stop("'J' must be a single integer >= 2 (level-2 units per arm).", call. = FALSE)
    }
    pwr <- power_at(J)
    ncp <- d * sqrt(J * n / (2 * (1 + (n - 1) * rho)))
    out <- data.frame(
      term  = c("specified_J_per_arm", "total_N", "noncentrality", "actual_power"),
      value = c(J, 2 * J * n, ncp, pwr)
    )
    class(out) <- c("dmar_ss_power", class(out))
    return(.as_dmar_tbl(out))
  }

  if (desired_power <= 0 || desired_power >= 1) {
    stop("'desired_power' must be in (0, 1).", call. = FALSE)
  }

  J_i <- 1
  pwr <- 0
  while (pwr < desired_power) {
    J_i <- J_i + 1
    pwr <- power_at(J_i)
    if (J_i > 1e6) stop("Failed to converge within reasonable sample size.", call. = FALSE)
  }
  ncp <- d * sqrt(J_i * n / (2 * (1 + (n - 1) * rho)))
  out <- data.frame(
    term  = c("necessary_J_per_arm", "total_N", "noncentrality", "actual_power"),
    value = c(J_i, 2 * J_i * n, ncp, pwr)
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
