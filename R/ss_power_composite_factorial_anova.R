#' Sample Size or Composite Power for a Factorial ANOVA
#'
#' Determine the necessary per-cell sample size to achieve a desired level of
#' composite statistical power in a balanced factorial analysis of variance, or,
#' given a per-cell sample size, return the realized composite power. Composite
#' power is the probability that every effect named in \code{effects} is
#' statistically significant in the same study, the quantity a design must be
#' planned against when its conclusion requires more than one result to hold at
#' once.
#'
#' This is the no-covariate case of \code{\link{ss_power_composite_factorial_ancova}},
#' named directly. It does not take a covariate; when a covariate belongs in the
#' model, use \code{\link{ss_power_composite_factorial_ancova}}, which raises
#' every effect's power for the variance the covariate explains. With no
#' covariate the composite is exact to quadrature precision, since the balanced
#' factorial \emph{F} tests are exactly noncentral \emph{F} and exactly
#' independent given the shared error estimate.
#'
#' @param factor_levels Integer vector of the number of levels of each factor,
#'   one entry per factor (each at least 2). A \code{c(2, 3, 2)} argument is a
#'   2 by 3 by 2 design.
#' @param effects A non-empty list naming the effects in the composite; see
#'   \code{\link{ss_power_composite_factorial_ancova}} for the format. Each
#'   element gives the \code{factors} an effect spans and its \code{f} or
#'   \code{partial_eta_squared}, unless \code{means} is supplied.
#' @param means Optional array of population cell means (dimensions
#'   \code{factor_levels}), or a numeric vector of length
#'   \code{prod(factor_levels)} in array order, from which each effect's Cohen's
#'   \emph{f} is read given \code{sigma}. When supplied, \code{plot()} draws the
#'   mean pattern. See \code{\link{ss_power_composite_factorial_ancova}}.
#' @param sigma The common within-cell standard deviation of the outcome,
#'   required with \code{means} and used only there.
#' @param desired_power Desired composite statistical power (default 0.85). Used
#'   only when \code{n_per_cell} is \code{NULL}.
#' @param alpha_level Type I error rate for each individual \emph{F} test
#'   (default 0.05), the per-test rate rather than a rate for the composite
#'   event.
#' @param n_per_cell Per-cell sample size, assumed balanced; if supplied, the
#'   realized composite power is returned rather than a sample size planned.
#'
#' @details
#' See \code{\link{ss_power_composite_factorial_ancova}} for the model, the
#' one-dimensional integral that evaluates the composite over the shared error
#' estimate, and the figure the \code{plot()} method draws. Naming one effect
#' reproduces \code{\link{ss_power_factorial_anova}} exactly.
#'
#' @return A \code{data.frame} with \code{term} and \code{value} columns, as
#'   \code{\link{ss_power_composite_factorial_ancova}} returns but with
#'   \code{covariate_R2} 0 and \code{n_covariates} 0. It carries the
#'   \code{dmar_ss_power} class for \code{\link[generics]{tidy}} /
#'   \code{\link[generics]{glance}} and the
#'   \code{dmar_composite_power_factorial} class for \code{plot()}.
#'
#' @references
#' Maxwell, S. E. (2004). The persistence of underpowered studies in
#'   psychological research: Causes, consequences, and remedies.
#'   \emph{Psychological Methods, 9}, 147--163.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective} (4th ed.).
#'   Routledge. (See Chapter 7 on factorial designs and Chapter 3 on statistical
#'   power.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_composite_factorial_ancova}} for the version
#'   that admits a covariate; \code{\link{ss_power_factorial_anova}} for a single
#'   effect
#'
#' @examples
#' # A 2 by 3 factorial ANOVA whose conclusion needs both main effects. Plan
#' # against their composite, not against either one alone.
#' ss_power_composite_factorial_anova(
#'   factor_levels = c(2, 3),
#'   effects = list(list(factors = 1, f = 0.25),
#'                  list(factors = 2, f = 0.20)),
#'   desired_power = 0.80)
#'
#' # Realized composite power at 25 per cell for a main effect and the
#' # three-way interaction of a 2 by 2 by 2 design.
#' ss_power_composite_factorial_anova(
#'   factor_levels = c(2, 2, 2),
#'   effects = list(list(factors = 1,          f = 0.30, label = "A"),
#'                  list(factors = c(1, 2, 3), f = 0.25, label = "AxBxC")),
#'   n_per_cell = 25)
#'
#' # Naming one effect reproduces the single-effect planner exactly.
#' ss_power_composite_factorial_anova(
#'   factor_levels = c(2, 3), effects = list(list(factors = 2, f = 0.25)),
#'   n_per_cell = 20)
#' ss_power_factorial_anova(factor_levels = c(2, 3), effect_indices = 2,
#'                          f = 0.25, n_per_cell = 20)
#'
#' # The figure of the purported population effect sizes.
#' plot(ss_power_composite_factorial_anova(
#'   factor_levels = c(2, 3),
#'   effects = list(list(factors = 1, f = 0.25),
#'                  list(factors = 2, f = 0.20)),
#'   n_per_cell = 30))
#'
#' # Stating the effects as population cell means with a common within-cell SD.
#' # plot() then draws the mean pattern, with error bars of one SD.
#' cell_means <- matrix(c(10, 12, 11,
#'                        13, 12, 16), nrow = 2, byrow = TRUE)
#' plot(ss_power_composite_factorial_anova(
#'   factor_levels = c(2, 3), means = cell_means, sigma = 4,
#'   effects = list(list(factors = 1, label = "A"),
#'                  list(factors = 2, label = "B")),
#'   n_per_cell = 30))
#'
#' @keywords design htest
#'
#' @family sample size for power
#'
#' @family composite power
#'
#' @export
ss_power_composite_factorial_anova <- function(factor_levels, effects,
                                               means = NULL, sigma = NULL,
                                               desired_power = 0.85,
                                               alpha_level = 0.05,
                                               n_per_cell = NULL) {
  ss_power_composite_factorial_ancova(
    factor_levels = factor_levels, effects = effects,
    means = means, sigma = sigma, covariate_R2 = 0, n_covariates = 0,
    desired_power = desired_power, alpha_level = alpha_level,
    n_per_cell = n_per_cell)
}
