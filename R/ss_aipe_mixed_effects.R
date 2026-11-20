# AIPE sample size planning for a fixed effect in a two-level mixed-effects model.
#' AIPE Sample Size Planning for a Fixed Effect in a Two-Level Mixed-Effects Model
#'
#' Computes the smallest number of clusters (i.e., Level 2 units) for which
#' the confidence interval for the fixed-effect slope of one predictor has
#' an expected full width no larger than \code{width}, in a balanced
#' two-level model with a random intercept. The predictor may be measured
#' on the individuals (i.e., Level 1) or on the clusters (i.e., Level 2),
#' and may be centered at its cluster mean, centered at its grand mean, or
#' left as it is. \code{\link{ss_aipe_mixed_effects_sensitivity}} simulates
#' the same design, analysis, and interval, so the two can be used
#' together: plan here, then check the plan there, including under values
#' other than the ones planned with.
#'
#' @param true_total_var Population variance of the outcome given the
#'   predictor (a planning value): the variance of the random intercept
#'   (between clusters) plus the residual variance (within clusters). It
#'   is the name the total variance has in the cluster randomized planners
#'   (\code{\link{ss_aipe_crd}}).
#' @param true_var_X Population variance of the predictor (a planning
#'   value). For an individual-level (i.e., Level 1) predictor it is the
#'   total variance over individuals, between and within clusters together
#'   (see \code{true_icc_X}); for a cluster-level (i.e., Level 2)
#'   predictor it is the variance over clusters.
#' @param true_icc Population intraclass correlation of the outcome given
#'   the predictor (a planning value): the share of \code{true_total_var}
#'   that lies between clusters. In [0, 1).
#' @param width Desired full width of the confidence interval for the
#'   slope.
#' @param cluster_size Number of individuals (i.e., Level 1 units) in each
#'   cluster (i.e., Level 2 unit), the same in every cluster. Default
#'   \code{20L}.
#' @param conf_level Confidence level. Default \code{0.95}.
#' @param predictor_level \code{1} (the default) for an individual-level
#'   predictor, which varies within a cluster, or \code{2} for a
#'   cluster-level predictor, which is the same for every individual of a
#'   cluster.
#' @param centering How the predictor enters the model.
#'   \code{"cluster_mean"} (the default) subtracts each cluster's mean
#'   from its individuals' values, so the slope is the within-cluster
#'   slope; it applies to an individual-level predictor only.
#'   \code{"grand_mean"} subtracts the overall mean, and \code{"none"}
#'   leaves the predictor as it is. These two give the same slope and the
#'   same interval (only the intercept differs), and both are offered so
#'   that the call can say what the analysis will do.
#' @param true_icc_X Population intraclass correlation of an
#'   individual-level predictor (a planning value): the share of
#'   \code{true_var_X} that lies between clusters. Default \code{0}, a
#'   predictor whose cluster means differ by sampling alone. It is not
#'   used with \code{predictor_level = 2}.
#'
#' @return A \code{data.frame} with rows for the necessary number of
#'   clusters \code{necessary_n_clusters}, the implied total sample size
#'   \code{total_N} (\code{necessary_n_clusters * cluster_size}), the
#'   desired \code{width}, the \code{expected_width} at the necessary
#'   number of clusters, the degrees of freedom \code{df} of the interval's
#'   \emph{t} quantile, and the inputs \code{true_icc}, \code{true_icc_X},
#'   \code{cluster_size}, and \code{predictor_level}. The centering is
#'   recorded in \code{attr(x, "centering")}.
#'
#' @details
#' \strong{The model.} For individual \eqn{i} of cluster \eqn{k},
#' \deqn{y_{ki} = \gamma_0 + \beta x_{ki} + u_k + e_{ki},}
#' with \eqn{u_k \sim N(0, \tau^2)} the random intercept of the cluster
#' and \eqn{e_{ki} \sim N(0, \sigma^2)} the residual of the individual,
#' in \eqn{K} clusters of \eqn{m} individuals each. With \eqn{V} for
#' \code{true_total_var} and \eqn{\rho} for \code{true_icc},
#' \eqn{\tau^2 = \rho V} and \eqn{\sigma^2 = (1 - \rho) V}. The slope has
#' the same value within and between clusters. Writing \eqn{\sigma^2_X}
#' for \code{true_var_X} and \eqn{\rho_X} for \code{true_icc_X}, an
#' individual-level predictor has variance \eqn{\rho_X \sigma^2_X}
#' between clusters and \eqn{\sigma^2_W = (1 - \rho_X)\sigma^2_X} within
#' them.
#'
#' \strong{The interval planned for.} The interval is
#' \eqn{\hat\beta \pm t_{1 - \alpha/2,\, \nu}\, \widehat{SE}(\hat\beta)},
#' with \eqn{\nu = K(m - 1) - 1} for an individual-level predictor and
#' \eqn{\nu = K - 2} for a cluster-level predictor (the regression of
#' \eqn{K} cluster means). Its expected width is
#' \eqn{2\, t_{1 - \alpha/2,\, \nu}\, E[\widehat{SE}(\hat\beta)]}, and
#' the function returns the smallest \eqn{K} at which that is no larger
#' than \code{width}. No value of the slope is needed: \eqn{\beta} does
#' not appear in any of the expressions below.
#'
#' \strong{The expected standard error.} Each expression follows from the
#' sums of squares of a balanced design.
#' \itemize{
#'   \item \emph{Individual-level predictor centered at its cluster mean.}
#'     The estimate is the within-cluster slope, and for a normally
#'     distributed predictor the expected standard error is exactly
#'     \deqn{E[\widehat{SE}(\hat\beta)] =
#'       \frac{\sigma}{\sigma_W \sqrt{K (m - 1) - 1}}.}
#'     Clustering of the outcome enters only through \eqn{\sigma}: the
#'     design effect \eqn{1 + (m - 1)\rho} that inflates the variance of a
#'     cluster-level quantity does not appear, because centering removes
#'     the differences among clusters from both variables.
#'   \item \emph{Cluster-level predictor.} The estimate is the slope of
#'     the regression of the cluster means of the outcome on the
#'     predictor, and for a normally distributed predictor the expected
#'     standard error is exactly
#'     \deqn{E[\widehat{SE}(\hat\beta)] =
#'       \frac{\sqrt{\tau^2 + \sigma^2 / m}}{\sigma_X \sqrt{K - 2}}.}
#'     Here the number of clusters is what buys precision: as \eqn{m}
#'     grows the numerator approaches \eqn{\tau} and not 0.
#'   \item \emph{Individual-level predictor centered at its grand mean, or
#'     not centered.} The estimate pools the within-cluster and the
#'     between-cluster information about the one slope, each weighted by
#'     its precision, and to first order
#'     \deqn{E[\widehat{SE}(\hat\beta)]^2 \approx \left[
#'       \frac{\{K (m - 1) - 1\}\, \sigma^2_W}{\sigma^2} +
#'       \frac{(K - 2)\,(\rho_X \sigma^2_X + \sigma^2_W / m)}
#'            {\tau^2 + \sigma^2 / m} \right]^{-1},}
#'     which reduces to the two expressions above when either part is
#'     absent. It is never larger than the value for the within-cluster
#'     slope, and the gain is what the assumption of one slope within and
#'     between clusters buys; when that assumption is in doubt, plan for
#'     the cluster-mean-centered analysis.
#' }
#' In 4,000 simulated data sets for each of six designs (clusters of 5 to
#' 20, 9 to 44 clusters, intraclass correlations of .05 to .30), the mean
#' width of the interval was within 0.3 percent of the expected width in
#' the four designs with an exact expression, and 0.5 and 1.3 percent
#' narrower than it in the two with a pooled slope, so the first-order
#' expression errs toward a slightly larger number of clusters.
#' \code{\link{ss_aipe_mixed_effects_sensitivity}} reports the mean width
#' and the coverage for any design.
#'
#' \strong{Scope.} One predictor, a random intercept, a fixed slope, and
#' clusters of equal size. For random slopes, cross-level interactions, or
#' unequal cluster sizes the variance changes and a Monte Carlo planner
#' should be used (Schoemann, Boulton, & Short, 2017). The two-level
#' random-intercept model is treated at length by Snijders and Bosker
#' (2012).
#'
#' @references
#' Schoemann, A. M., Boulton, A. J., & Short, S. D. (2017). Determining
#'   power and sample size for simple and complex mediation models.
#'   \emph{Social Psychological and Personality Science, 8}(4), 379--386.
#'   \doi{10.1177/1948550617715068}
#'
#' Snijders, T. A. B., & Bosker, R. J. (2012). \emph{Multilevel
#'   analysis: An introduction to basic and advanced multilevel
#'   modeling} (2nd ed.). Sage.
#'
#' @seealso \code{\link{ss_aipe_mixed_effects_sensitivity}} to check a
#'   plan by simulation; \code{\link{ss_power_mixed_effects}},
#'   \code{\link{var_icc}}, \code{\link{ss_aipe_icc}}
#'
#' @examples
#' # 1. An individual-level predictor centered at its cluster mean (the
#' #    default): clusters of 20, intraclass correlation .10, unit
#' #    variances, and a desired width of .20.
#' ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1, true_icc = 0.10,
#'                       width = 0.20, cluster_size = 20)
#'
#' # 2. The same predictor not centered, when 30 percent of its variance
#' #    lies between clusters.
#' ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1, true_icc = 0.10,
#'                       width = 0.20, cluster_size = 20,
#'                       centering = "none", true_icc_X = 0.30)
#'
#' # 3. A cluster-level predictor: far more clusters are needed for the
#' #    same width, because each cluster supplies one value of it.
#' ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1, true_icc = 0.10,
#'                       width = 0.20, cluster_size = 20,
#'                       predictor_level = 2, centering = "grand_mean")
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

ss_aipe_mixed_effects <- function(true_total_var, true_var_X, true_icc, width,
                                  cluster_size = 20L,
                                  conf_level = 0.95,
                                  predictor_level = 1,
                                  centering = c("cluster_mean", "grand_mean",
                                                "none"),
                                  true_icc_X = 0) {
  .dmar_check_conf_level(conf_level)
  if (!missing(width)) .dmar_check_width(width)
  for (nm in c("true_total_var", "true_var_X", "true_icc", "width", "cluster_size",
               "conf_level")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L)
      .dmar_refuse_argument(
        sprintf("'%s' must be a single numeric value.", nm), arg = nm)
  }
  if (true_total_var <= 0) .dmar_refuse_argument(
    "'true_total_var' must be positive.", arg = "true_total_var")
  if (true_var_X <= 0) .dmar_refuse_argument(
    "'true_var_X' must be positive.", arg = "true_var_X")
  if (true_icc < 0 || true_icc >= 1) .dmar_refuse_argument(
    "'true_icc' must be in [0, 1).", arg = "true_icc")
  if (width <= 0) .dmar_refuse_argument(
    "'width' must be positive.", arg = "width")
  if (cluster_size < 2) .dmar_refuse_argument(
    "'cluster_size' must be at least 2.", arg = "cluster_size")
  if (conf_level <= 0 || conf_level >= 1)
    .dmar_refuse_argument("'conf_level' must be in (0, 1).", arg = "conf_level")
  # With predictor_level = 2 the default centering (the cluster mean) has no
  # meaning, so an unset centering becomes the grand mean there; a centering
  # the user names is checked as given.
  centering <- if (missing(centering) && is.numeric(predictor_level) &&
                   length(predictor_level) == 1L &&
                   isTRUE(predictor_level == 2)) {
    "grand_mean"
  } else {
    match.arg(centering)
  }
  design <- .me_design(predictor_level, centering, true_icc_X)

  m <- cluster_size
  ew <- function(K) {
    .me_expected_width(K, m, true_total_var, true_icc, true_var_X,
                       true_icc_X, design$level, design$centering, conf_level)
  }

  n_cl <- .me_min_clusters(design$level)
  while (ew(n_cl) > width) {
    n_cl <- n_cl + 1L
    if (n_cl > 1e6L)
      .dmar_refuse_argument(
        "Required number of clusters exceeded 1e6; check inputs.", arg = NULL)
  }

  out <- data.frame(
    term  = c("necessary_n_clusters", "total_N", "width", "expected_width",
              "df", "true_icc", "true_icc_X", "cluster_size",
              "predictor_level"),
    value = c(n_cl, n_cl * m, width, ew(n_cl),
              .me_df(n_cl, m, design$level), true_icc, true_icc_X, m,
              design$level),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
  attr(out, "centering") <- design$centering
  out
}
