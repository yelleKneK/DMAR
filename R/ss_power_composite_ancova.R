#' Sample Size or Composite Power for a One-Way or Factorial ANCOVA
#'
#' Determine the necessary per-cell sample size to achieve a desired level of
#' composite statistical power in a balanced analysis of covariance with
#' \eqn{a} groups (a one-way design) or a factorial arrangement of factors, or,
#' given a per-cell sample size, return the realized composite power. Composite
#' power is the probability that every effect named in \code{effects} is
#' statistically significant in the same study, the quantity a design must be
#' planned against when its conclusion requires more than one result to hold at
#' once.
#'
#' This is the general entry point for the ANCOVA composite. The two-group
#' design has its own simpler interface in
#' \code{\link{ss_power_composite_ancova_2group}}, stated through a single
#' \code{smd} and one or two correlations; use that when there are exactly two
#' groups. Use this function for more than two groups or for a factorial design.
#' With no covariate, name the effects through \code{\link{ss_power_composite_anova}}
#' instead.
#'
#' \strong{Homogeneous or heterogeneous slopes.} The \code{slopes} argument
#' chooses the model. \code{"homogeneous"} (the default) assumes one common
#' covariate slope across the cells; the covariate is a variance reducer and the
#' composite is over the factorial mean effects. \code{"heterogeneous"} lets the
#' slope differ across cells, which makes the average covariate slope and the
#' factor-by-covariate slope heterogeneity into testable effects that can join
#' the mean effects in the composite. The heterogeneous one-way case is the
#' \eqn{a}-group generalization of the two-group ANCOVA composite: a group mean
#' effect, the covariate effect, and the group-by-covariate slope heterogeneity.
#'
#' This function is a thin dispatcher. \code{slopes = "homogeneous"} forwards to
#' \code{\link{ss_power_composite_factorial_ancova}} and
#' \code{slopes = "heterogeneous"} to
#' \code{\link{ss_power_composite_factorial_ancova_het}}; see those for the full
#' method, the exactness discussion, and the shape of the returned object. The
#' population effects can be stated as effect sizes or as population values
#' (cell means, and for heterogeneous slopes a covariate-outcome correlation per
#' cell) with a common within-cell standard deviation, from which the
#' \code{plot()} method draws the population pattern.
#'
#' @param factor_levels Integer vector of the number of levels of each factor,
#'   one entry per factor (each at least 2). A single value \code{a} is a
#'   one-way design with \eqn{a} groups; \code{c(2, 3)} is a 2 by 3 factorial.
#' @param effects A non-empty list naming the effects in the composite. For
#'   \code{slopes = "homogeneous"} each element has \code{factors} (the factor
#'   indices the effect spans) and, unless \code{means} is supplied, its effect
#'   size \code{f} or \code{partial_eta_squared}. For
#'   \code{slopes = "heterogeneous"} each element also has a \code{type}, one of
#'   \code{"mean"} (a factorial mean effect, the default), \code{"covariate"}
#'   (the average slope, spanning no factors), or \code{"slope"} (a
#'   factor-by-covariate slope heterogeneity). An optional \code{label} names the
#'   effect. See the forwarded functions for the exact grammar.
#' @param slopes The covariate-slope model, \code{"homogeneous"} (one common
#'   slope, the default) or \code{"heterogeneous"} (the slope may differ across
#'   cells, making the covariate and slope-heterogeneity effects testable).
#' @param means Optional array of population cell means (dimensions
#'   \code{factor_levels}, or a vector in array order) from which the mean
#'   effects' sizes are read; enables the mean-pattern figure.
#' @param sigma The common within-cell population standard deviation of the
#'   outcome, required with the population-values interface.
#' @param covariate_R2 For \code{slopes = "homogeneous"} only: proportion of the
#'   outcome's within-cell variance the covariate or covariates explain, in
#'   \eqn{[0, 1)}, which raises every effect's noncentrality through
#'   \eqn{f / \sqrt{1 - R^2}}. Defaults to 0.
#' @param n_covariates For \code{slopes = "homogeneous"} only: number of
#'   covariates, each spending one residual degree of freedom. Must be positive
#'   when \code{covariate_R2} is. Defaults to 0.
#' @param correlations For \code{slopes = "heterogeneous"} only: optional array
#'   of the population covariate-outcome correlation within each cell
#'   (dimensions \code{factor_levels}), which supplies the covariate and slope
#'   effects' sizes with the population-values interface.
#' @param sd_cov For \code{slopes = "heterogeneous"} only: population standard
#'   deviation of the covariate, the units the slopes and the figure are drawn
#'   in. Scale free, so it changes no power. Default 1.
#' @param desired_power Desired composite statistical power (default 0.85). Used
#'   only when \code{n_per_cell} is \code{NULL}.
#' @param alpha_level Type I error rate for each individual \emph{F} test
#'   (default 0.05), the per-test rate, not a rate for the composite event.
#' @param n_per_cell Per-cell sample size (balanced); if supplied, the realized
#'   composite power is returned rather than a sample size planned.
#'
#' @return The \code{data.frame} the forwarded planner returns, with \code{term}
#'   and \code{value} columns and the \code{dmar_ss_power} class for
#'   \code{\link[generics]{tidy}} and \code{\link[generics]{glance}}, plus the
#'   \code{dmar_composite_power_factorial} (homogeneous) or
#'   \code{dmar_composite_power_factorial_het} (heterogeneous) class for
#'   \code{plot()}. See \code{\link{ss_power_composite_factorial_ancova}} for the
#'   row-by-row description.
#'
#' @references
#' Maxwell, S. E. (2004). The persistence of underpowered studies in
#'   psychological research: Causes, consequences, and remedies.
#'   \emph{Psychological Methods, 9}, 147--163.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective} (4th ed.).
#'   Routledge. (See Chapter 9 on the analysis of covariance, Chapter 7 on
#'   factorial designs, and Chapter 3 on statistical power.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_composite_ancova_2group}} for the two-group
#'   special case with the simple \code{smd}/\code{rho} interface;
#'   \code{\link{ss_power_composite_anova}} for the no-covariate design;
#'   \code{\link{ss_power_composite_factorial_ancova}} and
#'   \code{\link{ss_power_composite_factorial_ancova_het}}, the planners this
#'   forwards to; \code{\link{ss_power_factorial_ancova}} for a single effect
#'
#' @examples
#' # A four-group (one-way) ANCOVA with heterogeneous slopes: the a-group
#' # generalization of the two-group composite. The conclusion needs the group
#' # mean effect, a covariate effect, and evidence that the covariate slope
#' # differs across the groups, so the design is planned against all three.
#' ss_power_composite_ancova(
#'   factor_levels = 4, slopes = "heterogeneous",
#'   effects = list(list(type = "mean",      factors = 1, f = 0.30),
#'                  list(type = "covariate",              f = 0.40),
#'                  list(type = "slope",      factors = 1, f = 0.20)),
#'   desired_power = 0.80)
#'
#' # A 2 by 3 factorial ANCOVA with one common slope. Both main effects must
#' # hold; the covariate explains 25 percent of the within-cell variance.
#' ss_power_composite_ancova(
#'   factor_levels = c(2, 3),
#'   effects = list(list(factors = 1, f = 0.25),
#'                  list(factors = 2, f = 0.20)),
#'   covariate_R2 = 0.25, n_covariates = 1,
#'   desired_power = 0.80)
#'
#' # The population effects can instead be a full pattern of cell means with a
#' # common within-cell SD; plot() then draws the mean pattern itself.
#' cell_means <- matrix(c(10, 12, 11,
#'                        13, 12, 16), nrow = 2, byrow = TRUE)
#' fit <- ss_power_composite_ancova(
#'   factor_levels = c(2, 3), means = cell_means, sigma = 4,
#'   effects = list(list(factors = 1, label = "A"),
#'                  list(factors = 2, label = "B")),
#'   n_per_cell = 30)
#' fit
#' plot(fit)
#'
#' # A one-row broom summary of a plan.
#' generics::tidy(ss_power_composite_ancova(
#'   factor_levels = c(2, 3),
#'   effects = list(list(factors = 1, f = 0.25),
#'                  list(factors = 2, f = 0.20)),
#'   desired_power = 0.80))
#'
#' # The two-group case matches the dedicated two-group planner. Here factor 1
#' # has two levels, so the heterogeneous one-way composite reproduces
#' # ss_power_composite_ancova_2group with one correlation per group.
#' ss_power_composite_ancova(
#'   factor_levels = 2, slopes = "heterogeneous",
#'   means = c(0, 0.5), correlations = c(0.1, 0.5), sigma = 1,
#'   effects = list(list(type = "mean", factors = 1),
#'                  list(type = "covariate"),
#'                  list(type = "slope", factors = 1)),
#'   n_per_cell = 95)
#' ss_power_composite_ancova_2group(smd = 0.5, rho = c(0.1, 0.5), n = 95)
#'
#' @keywords design htest
#'
#' @family sample size for power
#'
#' @family composite power
#'
#' @export
ss_power_composite_ancova <- function(factor_levels, effects,
                                      slopes = c("homogeneous",
                                                 "heterogeneous"),
                                      means = NULL, sigma = NULL,
                                      covariate_R2 = 0, n_covariates = 0,
                                      correlations = NULL, sd_cov = 1,
                                      desired_power = 0.85, alpha_level = 0.05,
                                      n_per_cell = NULL) {
  slopes <- match.arg(slopes)

  if (slopes == "heterogeneous") {
    if (!identical(covariate_R2, 0) && any(covariate_R2 != 0)) {
      stop("'covariate_R2' applies only to slopes = \"homogeneous\"; a ",
           "heterogeneous-slope model estimates the covariate effect rather ",
           "than treating it as a fixed variance reduction.", call. = FALSE)
    }
    if (!identical(n_covariates, 0) && any(n_covariates != 0)) {
      stop("'n_covariates' applies only to slopes = \"homogeneous\".",
           call. = FALSE)
    }
    return(ss_power_composite_factorial_ancova_het(
      factor_levels = factor_levels, effects = effects,
      means = means, correlations = correlations, sigma = sigma,
      sd_cov = sd_cov, desired_power = desired_power,
      alpha_level = alpha_level, n_per_cell = n_per_cell))
  }

  if (!is.null(correlations)) {
    stop("'correlations' applies only to slopes = \"heterogeneous\"; a ",
         "homogeneous-slope model summarizes the covariate through ",
         "'covariate_R2'.", call. = FALSE)
  }
  if (!identical(sd_cov, 1) && any(sd_cov != 1)) {
    stop("'sd_cov' applies only to slopes = \"heterogeneous\".", call. = FALSE)
  }
  ss_power_composite_factorial_ancova(
    factor_levels = factor_levels, effects = effects,
    means = means, sigma = sigma, covariate_R2 = covariate_R2,
    n_covariates = n_covariates, desired_power = desired_power,
    alpha_level = alpha_level, n_per_cell = n_per_cell)
}
