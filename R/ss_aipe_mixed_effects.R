# AIPE sample size planning for a fixed effect in a two-level mixed-effects model.
#' AIPE Sample Size Planning for a Fixed Effect in a Two-Level Mixed-Effects Model
#'
#' Computes the minimum number of clusters (level-2 units) needed so
#' that the confidence interval on a level-1 fixed-effect slope has
#' expected full width no larger than \eqn{\omega} (Kelley, 2007;
#' Raudenbush & Liu, 2001; Snijders & Bosker, 2012). The function
#' inverts the closed-form approximation for the variance of a
#' fixed-effect slope in a balanced two-level random-intercept model.
#'
#' @param sigma2_y Total variance of the outcome variable.
#' @param sigma2_x Variance of the level-1 predictor (covariate).
#' @param icc Intraclass correlation of the outcome.
#' @param width Target full CI width on the slope.
#' @param cluster_size Per-cluster sample size (number of level-1
#'   units per level-2 unit). Default \code{20L}.
#' @param conf_level Confidence level. Default \code{0.95}.
#'
#' @return A \code{data.frame} with rows for the recommended
#'   number of clusters \code{necessary_n_clusters}, the implied total
#'   sample size \code{total_N} (\code{necessary_n_clusters *
#'   cluster_size}), the target \code{width}, the intraclass
#'   correlation \code{icc}, the \code{cluster_size}, and the resulting
#'   \code{ci_width_expected} (the expected full CI width at the
#'   recommended size).
#'
#' @details
#' \strong{Variance of the slope.} For a level-1 predictor centered
#' within cluster, the asymptotic variance of \eqn{\hat\beta} is
#' approximately
#' \deqn{\mathrm{Var}(\hat\beta) \;\approx\;
#'   \frac{\sigma^2_y (1 - \rho_I)}{N \sigma^2_x},}
#' where \eqn{N = n_{\mathrm{clusters}} \cdot m} is the total number of
#' level-1 units, \eqn{m} is the cluster size, and \eqn{\rho_I} the
#' intraclass correlation. Because within-cluster centering removes the
#' cluster-level variation from the predictor, the design effect
#' \eqn{1 + (m - 1) \rho_I} that inflates the variance of a
#' cluster-level estimand does not appear here; clustering enters only
#' through the residual variance \eqn{\sigma^2_y (1 - \rho_I)}. The
#' function inverts this expression for \eqn{N}. No anticipated slope
#' value is needed: \eqn{\beta} does not appear in the variance, so the
#' recommended number of clusters is the same whatever the slope.
#'
#' \strong{Scope.} Planning is for the most common single-level
#' covariate case (random intercept, fixed slope, level-1 predictor
#' centered within cluster). For cross-level interactions or random
#' slopes, the variance formula changes and a Monte Carlo planner
#' should be used instead (Schoemann, Boulton, & Short, 2017).
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for
#'   standardized effect sizes: Theory, application, and
#'   implementation. \emph{Journal of Statistical Software, 20}(8),
#'   1--24. \doi{10.18637/jss.v020.i08}
#'
#' Raudenbush, S. W., & Liu, X.-F. (2001). Effects of study duration,
#'   frequency of observation, and sample size on power in studies of
#'   group differences in polynomial change. \emph{Psychological
#'   Methods, 6}(4), 387--401. \doi{10.1037/1082-989X.6.4.387}
#'
#' Schoemann, A. M., Boulton, A. J., & Short, S. D. (2017). Determining
#'   power and sample size for simple and complex mediation models.
#'   \emph{Social Psychological and Personality Science, 8}(4), 379--386.
#'   \doi{10.1177/1948550617715068}
#'
#' Snijders, T. A. B., & Bosker, R. J. (2012). \emph{Multilevel
#'   analysis: An introduction to basic and advanced multilevel
#'   modeling} (2nd ed.). Sage.
#'
#' @seealso \code{\link{ss_power_mixed_effects}}, \code{\link{var_icc}},
#'   \code{\link{ss_aipe_icc}}
#'
#' @examples
#' # 1. Plan a two-level study with cluster size 20, ICC = 0.10,
#' #        sigma_y = 1, sigma_x = 1, target CI width = 0.20:
#' ss_aipe_mixed_effects(sigma2_y = 1, sigma2_x = 1, icc = 0.10,
#'                            width = 0.20, cluster_size = 20)
#'
#' # 2. The same study with stronger clustering (ICC = 0.20):
#' ss_aipe_mixed_effects(sigma2_y = 1, sigma2_x = 1, icc = 0.20,
#'                            width = 0.20, cluster_size = 20)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#' @family mixed models
#'
#' @export

ss_aipe_mixed_effects <- function(sigma2_y, sigma2_x, icc, width,
                                  cluster_size = 20L,
                                  conf_level = 0.95) {
  for (nm in c("sigma2_y", "sigma2_x", "icc", "width", "cluster_size",
               "conf_level")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L)
      stop(sprintf("'%s' must be a single numeric value.", nm))
  }
  if (sigma2_y <= 0) stop("'sigma2_y' must be positive.")
  if (sigma2_x <= 0) stop("'sigma2_x' must be positive.")
  if (icc < 0 || icc >= 1) stop("'icc' must be in [0, 1).")
  if (width <= 0) stop("'width' must be positive.")
  if (cluster_size < 2) stop("'cluster_size' must be at least 2.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be in (0, 1).")

  m <- cluster_size
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  half_target <- width / 2

  # Variance of beta for level-1 predictor (cluster-mean-centered):
  half_width <- function(n_cl) {
    N <- n_cl * m
    var_beta <- sigma2_y * (1 - icc) / (N * sigma2_x)
    z * sqrt(var_beta)
  }

  n_cl <- 2L
  while (half_width(n_cl) > half_target) {
    n_cl <- n_cl + 1L
    if (n_cl > 1e6L)
      stop("Required number of clusters exceeded 1e6; check inputs.")
  }

  out <- data.frame(
    term  = c("necessary_n_clusters", "total_N", "width", "icc", "cluster_size",
              "ci_width_expected"),
    value = c(n_cl, n_cl * m, width, icc, m, 2 * half_width(n_cl)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
