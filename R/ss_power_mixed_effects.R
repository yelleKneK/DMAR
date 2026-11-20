#' Sample Size or Power for a Treatment Effect in a Two-Level Mixed-Effects Model
#'
#' Determine the necessary number of clusters (i.e., Level 2 units) per arm to achieve a desired level of
#' statistical power for a treatment-versus-control comparison in a two-level mixed-effects model
#' with a random intercept (e.g., individuals nested within clusters in a cluster-randomized trial,
#' or repeated measurements nested within subjects in a person-randomized longitudinal study).
#' Alternatively, given a number of clusters per arm, return the realized statistical power.
#'
#' @param true_smd Population standardized treatment effect the test is powered for (\eqn{d} in the formulas below), defined as the population mean difference divided by the population standard deviation of the outcome among individual observations (i.e., Level 1 units)
#' @param cluster_size Number of observations (i.e., Level 1 units) in each cluster (i.e., Level 2 unit); for example, individuals per cluster or measurements per subject, \eqn{n} in the formulas below; a single whole number of at least 2, assumed equal across clusters. It has no default
#' @param true_icc Population intraclass correlation (the proportion of total outcome variance attributable to differences between clusters (i.e., Level 2 units)), \eqn{\rho} in the formulas below; a planning value that must be in [0, 1)
#' @param specified_n_clusters_per_arm Number of clusters (i.e., Level 2 units) in each arm, \eqn{J} in the formulas below (i.e., \code{specified_n_clusters_per_arm} treatment clusters and as many control clusters, so the total number of clusters is twice this value); if specified, returns the realized power
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param directional Logical: \code{TRUE} for a one-sided test (in the same sign as \code{true_smd}), \code{FALSE} (default) for a two-sided test
#'
#' @details
#' This function computes power for the fixed treatment effect at the higher level of a two-level
#' mixed-effects model with random intercept,
#' \deqn{y_{ij} = \beta_0 + \beta_1 T_j + u_j + \epsilon_{ij},}
#' where \eqn{u_j \sim N(0, \sigma_u^2)} is the cluster (i.e., Level 2) random intercept and
#' \eqn{\epsilon_{ij} \sim N(0, \sigma_e^2)} is the within-cluster (i.e., Level 1) residual. The treatment indicator
#' \eqn{T_j} varies between clusters (i.e., entire clusters or entire subjects are assigned
#' to treatment or control). The intra-class correlation is \eqn{\rho = \sigma_u^2 / (\sigma_u^2
#' + \sigma_e^2)} and the total outcome variance is \eqn{\sigma_y^2 = \sigma_u^2 + \sigma_e^2}.
#'
#' With \eqn{n} observations in each cluster (\code{cluster_size}) and \eqn{J} clusters in each arm
#' (\code{specified_n_clusters_per_arm}), the standard error of the estimated treatment effect is
#' \eqn{SE(\hat\beta_1) = \sigma_y \sqrt{2 (1 + (n - 1)\rho) / (J n)}}, giving a noncentrality
#' parameter of
#' \deqn{\lambda = d \sqrt{J n / (2 (1 + (n - 1)\rho))}}
#' under a two-sample \emph{t}-test with \eqn{2J - 2} degrees of freedom.
#'
#' The factor \eqn{1 + (n - 1)\rho} is the design effect: as the within-cluster correlation grows,
#' the effective information per cluster shrinks, so more clusters are needed for a given
#' level of power.
#'
#' @return A \code{data.frame} with rows for \code{necessary_n_clusters_per_arm} (or
#'   \code{specified_n_clusters_per_arm} when \code{specified_n_clusters_per_arm} is given), \code{total_N} (the total
#'   number of observations (i.e., Level 1 units), \eqn{2 J n}), \code{noncentral_t_parm} (the noncentrality parameter
#'   \eqn{\lambda} of the \emph{t} distribution), and \code{actual_power}. The result
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
#' # 30 individuals per cluster, ICC = 0.05, standardized effect 0.30, power = .80
#' ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.05, desired_power = 0.80)
#'
#' # Same effect but with much higher ICC (e.g., schools or therapists)
#' ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.20, desired_power = 0.80)
#'
#' # Realized power with 25 clusters per arm
#' ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.05,
#'                        specified_n_clusters_per_arm = 25)
#'
#' # directional test
#' ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.05,
#'                        desired_power = 0.80, directional = TRUE)
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
ss_power_mixed_effects <- function(true_smd, cluster_size, true_icc,
                                   specified_n_clusters_per_arm = NULL,
                                   desired_power = 0.85, alpha_level = 0.05,
                                   directional = FALSE) {
  .dmar_check_alpha_level(alpha_level)
  # cluster_size is checked first and has no default: a call written for
  # 1.0.x that gave the intraclass correlation by position leaves it
  # missing, and must stop here and say so. (A 1.0.x call that named `n`
  # stops earlier, with R's "unused argument".)
  if (missing(cluster_size)) {
    .dmar_refuse_argument(c(
      "'cluster_size' is missing. It is the number of observations in ",
      "each cluster (individuals per cluster, or measurements per ",
      "subject; i.e., Level 1 units in each Level 2 unit) and has no default, because the design effect ",
      "1 + (cluster_size - 1) * true_icc depends on it. Give a whole ",
      "number of at least 2, such as cluster_size = 30. The intraclass ",
      "correlation is 'true_icc', and the number of clusters per arm is ",
      "'specified_n_clusters_per_arm'."),
      arg = "cluster_size")
  }
  if (!is.numeric(true_smd) || length(true_smd) != 1L || !is.finite(true_smd) || true_smd == 0) {
    .dmar_refuse_argument(
      "'true_smd' must be a single finite nonzero numeric value.",
      arg = "true_smd")
  }
  if (!is.numeric(cluster_size) || length(cluster_size) != 1L ||
      !is.finite(cluster_size) || cluster_size < 2 ||
      cluster_size != round(cluster_size)) {
    .dmar_refuse_argument(c(
      "'cluster_size' must be a single whole number of at least 2: the ",
      "number of observations (i.e., Level 1 units) in each cluster. With one per ",
      "cluster there is no clustering to plan for; use ss_power_smd() ",
      "for two independent groups."),
      arg = "cluster_size", value = cluster_size)
  }
  if (missing(true_icc) || !is.numeric(true_icc) || length(true_icc) != 1L ||
      !is.finite(true_icc) || true_icc < 0 || true_icc >= 1) {
    .dmar_refuse_argument(c(
      "'true_icc' must be a single number in [0, 1): the population ",
      "intraclass correlation, the share of the outcome variance that ",
      "lies between clusters. A value of 1 would leave no variance ",
      "within clusters. Give a planning value such as true_icc = 0.05."),
      arg = "true_icc", value = if (missing(true_icc)) NULL else true_icc)
  }
  if (alpha_level <= 0 || alpha_level >= 1) {
    .dmar_refuse_argument(
      "'alpha_level' must be in (0, 1).", arg = "alpha_level")
  }

  ncp_at <- function(J) {
    true_smd * sqrt(J * cluster_size / (2 * (1 + (cluster_size - 1) * true_icc)))
  }
  power_at <- function(J) {
    df  <- 2 * J - 2
    if (df <= 0) return(NA_real_)
    ncp <- ncp_at(J)
    # A directional test is taken in the direction of true_smd; the power
    # is that of |ncp| (R/ss_power_internals.R).
    .power_noncentral_t(ncp, df, alpha_level, directional)
  }

  # An invalid desired_power is refused whether or not it is used: with a
  # specified sample size the power is computed, not planned for.
  if (!is.null(desired_power)) .dmar_check_desired_power(desired_power)
  if (!is.null(specified_n_clusters_per_arm)) {
    if (!is.numeric(specified_n_clusters_per_arm) || length(specified_n_clusters_per_arm) != 1L ||
        !is.finite(specified_n_clusters_per_arm) || specified_n_clusters_per_arm < 2 ||
        specified_n_clusters_per_arm != round(specified_n_clusters_per_arm)) {
      .dmar_refuse_argument(c(
        "'specified_n_clusters_per_arm' must be a single whole number of at least 2: ",
        "the number of clusters in each of the two arms. The test has ",
        "2 * specified_n_clusters_per_arm - 2 degrees of freedom, which is 0 with one ",
        "cluster per arm. Leave it NULL to have the number found for ",
        "'desired_power'."),
        arg = "specified_n_clusters_per_arm", value = specified_n_clusters_per_arm)
    }
    out <- data.frame(
      term  = c("specified_n_clusters_per_arm", "total_N", "noncentral_t_parm", "actual_power"),
      value = c(specified_n_clusters_per_arm, 2 * specified_n_clusters_per_arm * cluster_size,
                ncp_at(specified_n_clusters_per_arm), power_at(specified_n_clusters_per_arm))
    )
    class(out) <- c("dmar_ss_power", class(out))
    return(.as_dmar_tbl(out))
  }

  .dmar_check_desired_power(desired_power)

  J_i <- 1
  pwr <- 0
  while (pwr < desired_power) {
    J_i <- J_i + 1
    pwr <- power_at(J_i)
    if (J_i > 1e6) .dmar_refuse_no_sample_size(
      "ss_power_mixed_effects", desired_power, size = "number of clusters per arm")
  }
  out <- data.frame(
    term  = c("necessary_n_clusters_per_arm", "total_N", "noncentral_t_parm", "actual_power"),
    value = c(J_i, 2 * J_i * cluster_size, ncp_at(J_i), pwr)
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
