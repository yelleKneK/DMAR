#' Sample Size or Composite Power for a One-Way or Factorial ANOVA
#'
#' Determine the necessary per-cell sample size to achieve a desired level of
#' composite statistical power in a balanced analysis of variance with \eqn{a}
#' groups (a one-way design) or a factorial arrangement of factors, or, given a
#' per-cell sample size, return the realized composite power. Composite power is
#' the probability that every effect named in \code{effects} is statistically
#' significant in the same study, the quantity a design must be planned against
#' when its conclusion requires more than one result to hold at once. Each
#' effect is a main effect or an interaction, tested by its own \emph{F} test,
#' and any subset of them can make up the composite.
#'
#' The population effects can be stated two ways: as effect sizes, a Cohen's
#' \emph{f} or partial eta squared per effect, or as a full array of population
#' cell means together with a common within-cell standard deviation, from which
#' each named effect's \emph{f} is read off the analysis of variance
#' decomposition of the means. Supplying means lets the \code{plot()} method draw
#' the mean pattern itself, so the size of the population effects is on the page.
#'
#' If the design includes one or more covariates, use
#' \code{\link{ss_power_composite_ancova}}.
#'
#' @param factor_levels Integer vector of the number of levels of each factor,
#'   one entry per factor (each at least 2). A single value \code{a} is a
#'   one-way design with \eqn{a} groups; \code{c(2, 3)} is a 2 by 3 factorial.
#' @param effects A non-empty list naming the effects in the composite. Each
#'   element is a list with \code{factors} (a vector of factor indices into
#'   \code{factor_levels}: one index for a main effect, several for an
#'   interaction) and, unless \code{true_means} is supplied, exactly one of \code{true_f}
#'   (Cohen's \emph{f}) or \code{true_eta_squared_partial}. An optional \code{label}
#'   names the effect in the output and the figure; the default label is the
#'   factor indices joined by \code{x}. The purported effect sizes are
#'   population values the researcher posits, never sample estimates. A
#'   field under any other name is refused.
#' @param true_means Optional array of population cell means whose dimensions are
#'   \code{factor_levels} (a matrix for two factors), or a numeric vector of
#'   length \code{prod(factor_levels)} in array order. When supplied, each named
#'   effect's Cohen's \emph{f} is computed from the means and \code{true_sd}, and
#'   \code{plot()} draws the mean pattern.
#' @param true_sd The common within-cell population standard deviation of the
#'   outcome, required with \code{true_means} and used only there.
#' @param desired_power Desired composite statistical power (default 0.85). Used
#'   only when \code{specified_n_per_cell} is \code{NULL}.
#' @param alpha_level Type I error rate for each individual \emph{F} test
#'   (default 0.05), the per-test rate, not a rate for the composite event.
#' @param specified_n_per_cell Per-cell sample size (balanced); if supplied, the realized
#'   composite power is returned rather than a sample size planned.
#'
#' @return A \code{data.frame} with \code{term} and \code{value} columns: the
#'   \code{necessary_n_per_cell} (or the supplied \code{specified_n_per_cell})
#'   and the total sample size it implies, \code{total_N}, the
#'   \code{composite_power}, the \code{df_error}, then for each effect its
#'   marginal \code{power_<label>}, purported \code{true_f_<label>}, numerator
#'   \code{df_<label>}, and \code{noncentral_f_parm_<label>}, followed by
#'   \code{cells} and \code{alpha_level}. The result carries the
#'   \code{dmar_ss_power} class for \code{\link[generics]{tidy}} and
#'   \code{\link[generics]{glance}}, and a \code{dmar_composite_power_factorial}
#'   class so \code{plot()} draws the figure.
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
#' @seealso \code{\link{ss_power_composite_ancova}} for the design with one or
#'   more covariates; \code{\link{ss_power_factorial_anova}} and
#'   \code{\link{ss_power_one_way_anova}} for a single effect
#'
#' @examples
#' # A 2 by 3 factorial ANOVA whose conclusion needs both main effects, so the
#' # design is planned against their composite.
#' ss_power_composite_anova(
#'   factor_levels = c(2, 3),
#'   effects = list(list(factors = 1, true_f = 0.25),
#'                  list(factors = 2, true_f = 0.20)),
#'   desired_power = 0.80)
#'
#' # Realized composite power at 40 per cell for a main effect and the
#' # interaction of a 2 by 2 design, sizes given as partial eta squared.
#' ss_power_composite_anova(
#'   factor_levels = c(2, 2),
#'   effects = list(list(factors = 1,       true_eta_squared_partial = 0.06),
#'                  list(factors = c(1, 2), true_eta_squared_partial = 0.04)),
#'   specified_n_per_cell = 40)
#'
#' # The effects can instead be a full pattern of population cell means with a
#' # common within-cell SD; plot() then draws the mean pattern itself.
#' cell_means <- matrix(c(10, 12, 11,
#'                        13, 12, 16), nrow = 2, byrow = TRUE)
#' fit <- ss_power_composite_anova(
#'   factor_levels = c(2, 3), true_means = cell_means, true_sd = 4,
#'   effects = list(list(factors = 1, label = "A"),
#'                  list(factors = 2, label = "B")),
#'   specified_n_per_cell = 30)
#' fit
#' plot(fit)
#'
#' @keywords design htest
#'
#' @family sample size for power
#'
#' @family composite power
#'
#' @export
ss_power_composite_anova <- function(factor_levels, effects,
                                     true_means = NULL, true_sd = NULL,
                                     desired_power = 0.85, alpha_level = 0.05,
                                     specified_n_per_cell = NULL) {
  out <- ss_power_composite_factorial_anova(
    factor_levels = factor_levels, effects = effects,
    true_means = true_means, true_sd = true_sd, desired_power = desired_power,
    alpha_level = alpha_level, specified_n_per_cell = specified_n_per_cell)
  # The ANCOVA engine echoes true_covariate_R2 and n_covariates rows (both 0 here).
  # This is the ANOVA entry point, which surfaces no covariate, so drop those
  # two rows while preserving the class, the attributes plot() reads, and the
  # numeric value column. Subsetting resets the attributes, so re-apply them.
  keep <- !out$term %in% c("true_covariate_R2", "n_covariates")
  if (all(keep)) return(out)
  ats <- attributes(out)
  out <- out[keep, , drop = FALSE]
  ats$row.names <- seq_len(nrow(out))
  attributes(out) <- ats
  out
}
