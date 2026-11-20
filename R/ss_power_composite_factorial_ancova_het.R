#' Composite Power for a Factorial ANCOVA With Heterogeneous Slopes
#'
#' Determine the necessary per-cell sample size, or the realized composite power
#' at a supplied per-cell size, for a balanced factorial analysis of covariance
#' in which the covariate's slope may differ across the cells. When the slopes
#' differ, the covariate main effect (the average slope) and the
#' factor-by-covariate slope heterogeneity are themselves testable effects, and
#' any of them, together with the factorial mean effects, can make up the
#' composite. This is the factorial generalization of
#' \code{\link{ss_power_composite_ancova_2group}}, whose two-group model with a
#' correlation per group is the one-factor, two-level case here.
#'
#' Kept separate from \code{\link{ss_power_composite_factorial_ancova}}, which
#' assumes one common slope and treats the covariate only as a variance reducer.
#' Use this function when the covariate slope is expected to differ across
#' conditions, or when the test of that difference is part of the design.
#'
#' @param factor_levels Integer vector of the number of levels of each factor
#'   (each at least 2).
#' @param effects A non-empty list naming the effects in the composite. Each
#'   element has a \code{type}, one of \code{"mean"} (a factorial main effect or
#'   interaction on the cell means; the default), \code{"covariate"} (the
#'   average covariate slope), or \code{"slope"} (a factor-by-covariate slope
#'   heterogeneity). A \code{"mean"} or \code{"slope"} effect also gives
#'   \code{factors}, the factor indices it spans; a \code{"covariate"} effect
#'   spans none. Each element carries its effect size, \code{true_f} or
#'   \code{true_eta_squared_partial}, unless the population values (\code{true_means} and
#'   \code{true_cor}) are supplied, in which case the sizes come from them.
#'   An optional \code{label} names the effect; defaults are the factor indices
#'   joined by \code{x} for a mean effect, \code{"covariate"}, and
#'   \code{"cov_x_<factors>"} for a slope effect. A field under any other
#'   name is refused with an error of class \code{dmar_argument_invalid}.
#' @param true_means Optional array of population cell means (dimensions
#'   \code{factor_levels}), needed when a \code{"mean"} effect is in the
#'   composite. Supplies the mean effects' sizes.
#' @param true_cor Optional array of the population covariate-outcome
#'   correlation within each cell (dimensions \code{factor_levels}, each in
#'   \eqn{(-1, 1)}); the correlation in cell \emph{c} is \eqn{\rho_c} in
#'   the formulas below, and \eqn{\bar{\rho^2}} is the mean of their
#'   squares. The cell slopes it implies supply the covariate and slope
#'   effects' sizes, and its spread across cells sets the slope heterogeneity.
#'   Required with the population-values interface, since the pooled error
#'   depends on every cell's correlation.
#' @param true_sd The common within-cell population standard deviation of the
#'   outcome (before adjustment), \eqn{\sigma} in the formulas below;
#'   required with the population-values interface.
#' @param true_sd_X Population standard deviation of the covariate. Default 1. A
#'   correlation is scale free, so \code{true_sd_X} does not change any power; it
#'   sets the units the slopes and the figure are drawn in.
#' @param desired_power Desired composite statistical power (default 0.85). Used
#'   only when \code{specified_n_per_cell} is \code{NULL}.
#' @param alpha_level Type I error rate for each individual test (default 0.05).
#' @param specified_n_per_cell Per-cell sample size (balanced); if supplied, the realized
#'   composite power is returned rather than a sample size planned.
#'
#' @details
#' The model fits the factorial mean structure, the covariate, and every
#' factor-by-covariate slope term, so it has \eqn{2 \times \mathrm{cells}}
#' parameters and residual degrees of freedom \eqn{N - 2\,\mathrm{cells}}. Under
#' balance and a covariate with a common distribution across the cells, those
#' terms are mutually orthogonal, so the tests share only the pooled residual
#' and the composite is the shared-error integral of
#' \code{\link{ss_power_composite_factorial_ancova}}. A "mean" effect on factor
#' set \eqn{S} has numerator df \eqn{\prod(a - 1)} and is read from the cell
#' means against the pooled adjusted error \eqn{\sigma^2(1 - \bar{\rho^2})}; a
#' "covariate" effect is the grand slope on 1 df; a "slope" effect on \eqn{S} is
#' the slope heterogeneity across those factors, with the same df as the matching
#' mean effect, read from the cell slopes.
#'
#' The two approximations are those of
#' \code{\link{ss_power_composite_ancova_2group}}. The covariate-related tests
#' condition on the covariate cross-products, overstating power by an amount of
#' order \eqn{1 / N}, which is negligible past a few dozen per cell. Separately,
#' correlations that differ in absolute value across cells make the pooled error
#' a mixture of scaled chi squares rather than the single one the integral
#' assumes; that one does not shrink with \emph{N} and is the subject of the
#' section below.
#'
#' @section Unequal Residual Variances When the Cell Correlations Differ:
#' In cell \emph{c} the residual variance is
#' \eqn{\tau_c^2 = \sigma^2 (1 - \rho_c^2)}, so cell correlations that differ in
#' absolute value leave the cells with different residual variances. Averaging
#' the squared correlations, which is what \eqn{\sigma^2(1 - \bar{\rho^2})}
#' does, gets the expected pooled error right and its distribution wrong: the
#' residual sum of squares is \eqn{\sum_c \tau_c^2 Q_c} with the \eqn{Q_c}
#' independent chi square variables, a mixture of scaled chi squares with the
#' same mean and a larger spread than the single scaled chi square the
#' shared-error integral integrates over. The numerators are affected too. The
#' covariate effect is the average cell slope and a slope effect is a contrast
#' among the cell slopes, and cell slopes estimated with different residual
#' variances give estimators of those two that are correlated, so the tests are
#' not independent even given the error estimate.
#'
#' Every power the function reports is then an approximation, and the output
#' says so. \code{composite_power} is reported as
#' \code{approximate_composite_power}, each \code{power_<label>} as
#' \code{approximate_power_<label>}, and a planned sample size as
#' \code{approximate_n_per_cell} and \code{approximate_total_N}, because that size is
#' the smallest one at which the approximation reaches \code{desired_power}
#' rather than a size known to attain it. The numbers do not change; the names
#' do, and only in this case. Equal absolute cell correlations, including cells
#' whose correlations share a magnitude and differ in sign, are the exact case
#' and keep the ordinary names.
#'
#' The error has no guaranteed sign, and it grows with the spread of the cell
#' correlations rather than shrinking with \emph{N}. Treat the number as an
#' approximation and confirm a design you intend to run by simulating it: draw
#' each cell's errors with its own residual standard deviation
#' \eqn{\sigma \sqrt{1 - \rho_c^2}}, fit the same model, and count the
#' replications in which every effect in the composite is significant. The
#' two-group help page reports the size of the departure over a range of
#' correlation gaps.
#'
#' The effect size interface is a special case worth naming. Supplying \code{true_f}
#' or \code{true_eta_squared_partial} for each effect states the effects directly and
#' says nothing about the cell correlations, so the function has nothing to
#' detect unequal residual variances from and labels the table exact. If the
#' design those effect sizes came from has cell correlations that differ in
#' absolute value, the same approximation applies and the labels will not tell
#' you; supply \code{true_cor} instead when you want the function to keep
#' track of it.
#'
#' @section The figure:
#' When the population values are supplied, \code{plot()} draws the population
#' regression line in each cell over the covariate, so heterogeneous slopes show
#' as lines of different angle and mean effects as vertical separation, colored
#' by the first factor and faceted by any others. When effect sizes are supplied
#' instead, it draws the effect size lollipop of
#' \code{\link{ss_power_composite_factorial_ancova}}. Either way the composite
#' power is in the subtitle. Requires \pkg{ggplot2}.
#'
#' @return A \code{data.frame} with \code{term} and \code{value} columns: the
#'   \code{necessary_n_per_cell} (or the supplied \code{specified_n_per_cell})
#'   and the total sample size it implies, \code{total_N}, the
#'   \code{composite_power}, the \code{df_error}, then for each effect its
#'   marginal \code{power_<label>}, purported \code{true_f_<label>}, numerator
#'   \code{df_<label>}, and \code{noncentral_f_parm_<label>}, followed by
#'   \code{cells} and \code{alpha_level}. Carries the \code{dmar_ss_power} class
#'   for \code{\link[generics]{tidy}} / \code{\link[generics]{glance}} and a
#'   \code{dmar_composite_power_factorial_het} class for \code{plot()}.
#'
#'   When the supplied cell correlations differ in absolute value the powers are
#'   approximations and the row names say so, with an \code{approximate_} prefix
#'   on every power and on a planned sample size; see the section on unequal
#'   residual variances. An \code{approximate} attribute carries the same flag
#'   for a program to test without parsing row names, and
#'   \code{\link[generics]{tidy}} and \code{\link[generics]{glance}} read both
#'   sets of names, so a relabeled table still summarizes to its per-cell size
#'   and its composite power.
#'
#' @references
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective} (4th ed.).
#'   Routledge. (See Chapter 9 on the analysis of covariance and heterogeneity of
#'   regression, and Chapter 7 on factorial designs.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_composite_factorial_ancova}} for the common-slope
#'   version; \code{\link{ss_power_composite_ancova_2group}} for the two-group case
#'
#' @examples
#' # A 2 by 2 design whose conclusion needs a main effect and evidence that the
#' # covariate slope differs across the first factor. Effect sizes stated
#' # directly: the main effect (mean), the average covariate effect, and the
#' # factor 1 by covariate slope heterogeneity.
#' ss_power_composite_factorial_ancova_het(
#'   factor_levels = c(2, 2),
#'   effects = list(list(type = "mean",      factors = 1,  true_f = 0.25),
#'                  list(type = "covariate",               true_f = 0.40),
#'                  list(type = "slope",      factors = 1,  true_f = 0.20)),
#'   desired_power = 0.80)
#'
#' # The same design from population values: cell means, a covariate-outcome
#' # correlation per cell (they differ across factor 1, so the slopes do), and a
#' # common within-cell SD. plot() then draws the per-cell regression lines.
#' cell_means <- matrix(c(10, 12,
#'                        11, 13), nrow = 2, byrow = TRUE)
#' cell_rho   <- matrix(c(0.55, 0.55,
#'                        0.15, 0.20), nrow = 2, byrow = TRUE)
#' fit <- ss_power_composite_factorial_ancova_het(
#'   factor_levels = c(2, 2), true_means = cell_means, true_cor = cell_rho,
#'   true_sd = 4, true_sd_X = 2,
#'   effects = list(list(type = "mean",      factors = 1),
#'                  list(type = "covariate"),
#'                  list(type = "slope",      factors = 1)),
#'   specified_n_per_cell = 40)
#' fit
#' plot(fit)
#'
#' # The one-factor, two-level case is the two-group composite ANCOVA.
#' ss_power_composite_factorial_ancova_het(
#'   factor_levels = 2,
#'   true_means = c(-0.25, 0.25), true_cor = c(0.1, 0.5), true_sd = 1,
#'   effects = list(list(type = "mean", factors = 1),
#'                  list(type = "slope", factors = 1)),
#'   specified_n_per_cell = 100)
#'
#' @keywords design htest
#'
#' @family sample size for power
#'
#' @family composite power
#'
#' @export
ss_power_composite_factorial_ancova_het <- function(factor_levels, effects,
                                                    true_means = NULL,
                                                    true_cor = NULL,
                                                    true_sd = NULL, true_sd_X = 1,
                                                    desired_power = 0.85,
                                                    alpha_level = 0.05,
                                                    specified_n_per_cell = NULL) {
  if (!is.numeric(factor_levels) || !length(factor_levels) ||
      any(is.na(factor_levels)) || any(factor_levels < 2) ||
      any(factor_levels != round(factor_levels))) {
    .dmar_refuse_argument(
      "'factor_levels' must be a vector of integers, each at least 2.",
      arg = "factor_levels")
  }
  factor_levels <- as.integer(factor_levels)
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L ||
      is.na(alpha_level) || alpha_level <= 0 || alpha_level >= 1) {
    .dmar_refuse_argument(
      "'alpha_level' must be a single number in (0, 1).", arg = "alpha_level")
  }
  if (!is.numeric(true_sd_X) || length(true_sd_X) != 1L || is.na(true_sd_X) ||
      true_sd_X <= 0) {
    .dmar_refuse_argument(
      "'true_sd_X' must be a single positive number.", arg = "true_sd_X")
  }

  from_pop <- !is.null(true_means) || !is.null(true_cor)
  if (from_pop) {
    if (is.null(true_cor)) {
      .dmar_refuse_argument(c(
        "The population-values interface needs 'true_cor', a covariate ",
        "outcome correlation for every cell, since the pooled error depends ",
        "on all of them."), arg = "true_cor")
    }
    if (is.null(true_sd) || !is.numeric(true_sd) || length(true_sd) != 1L ||
        is.na(true_sd) || true_sd <= 0) {
      .dmar_refuse_argument(c(
        "With population values, supply a single positive 'true_sd', the ",
        "common within-cell standard deviation."), arg = "true_sd")
    }
    true_cor <- .as_factorial_means(true_cor, factor_levels)
    if (any(abs(true_cor) >= 1)) {
      .dmar_refuse_argument(
        "Every entry of 'true_cor' must lie in (-1, 1).",
        arg = "true_cor")
    }
    # Correlations that differ in absolute value across cells make the cells'
    # residual variances differ, so the pooled error the tests share is a
    # mixture of scaled chi squares rather than the single one the composite
    # integral assumes, and the average slope and slope heterogeneity
    # estimators become correlated. Every power is then an approximation, and a
    # resolved sample size is the smallest one at which the approximation
    # reaches the target rather than a size known to attain it. Warn, and carry
    # the same statement into the row names.
    approximate <- length(unique(round(abs(as.vector(true_cor)), 12))) > 1L
    if (approximate) {
      .dmar_warn(c(
        "The cell correlations differ in absolute value, so the cells' ",
        "residual variances differ and the shared-error composite power ",
        "is an approximation whose error grows with the spread of the ",
        "correlations. The powers are reported as approximate_* rows, ",
        "and a planned sample size as approximate_n_per_cell, because ",
        "neither is exact; confirm a final design by simulation."),
        "result_qualified")
    }
    effects <- .parse_factorial_het_effects(effects, length(factor_levels),
                                            require_f = FALSE)
    if (any(vapply(effects, function(e) e$type == "mean", logical(1)))) {
      if (is.null(true_means)) {
        .dmar_refuse_argument(
          "A \"mean\" effect is in the composite, so supply 'true_means'.",
          arg = "true_means")
      }
      true_means <- .as_factorial_means(true_means, factor_levels)
    } else {
      true_means <- array(0, dim = factor_levels)   # unused when no mean effect
    }
    effects <- lapply(effects, function(e) {
      e$f <- .factorial_het_effect_f(e, true_means, true_cor,
                                     true_sd, true_sd_X)
      e
    })
    zero <- vapply(effects, function(e) e$f < 1e-9, logical(1))
    if (any(zero)) {
      .dmar_refuse_argument(c(
        "These effects do not vary in the supplied population values ",
        "(Cohen's f is 0), so no sample size attains a power above the Type ",
        "I error rate: ",
        paste(vapply(effects[zero], function(e) e$label, character(1)),
              collapse = ", "),
        ". Drop them, or supply values in which they vary."), arg = NULL)
    }
  } else {
    if (!is.null(true_sd)) {
      .dmar_refuse_argument(c(
        "'true_sd' is only used with the population-values interface. Supply ",
        "effect sizes in 'effects', or supply 'true_means' and 'true_cor'."),
        arg = c("true_sd", "effects", "true_means"))
    }
    effects <- .parse_factorial_het_effects(effects, length(factor_levels),
                                            require_f = TRUE)
    # Effect sizes state the effects directly and say nothing about the cells'
    # correlations, so there is nothing here to detect unequal residual
    # variances from. The table is labeled exact because the function cannot
    # tell otherwise; the documentation says what to do about that.
    approximate <- FALSE
  }

  cells <- prod(factor_levels)

  evaluate_at <- function(n) {
    d   <- .factorial_het_design(factor_levels, effects, n)
    res <- .composite_power_shared_sigma_F(d$df_num, d$lambda, d$df_error,
                                           alpha_level)
    list(composite = res$composite, marginal = res$marginal, design = d)
  }

  # An invalid desired_power is refused whether or not it is used: with a
  # specified sample size the power is computed, not planned for.
  if (!is.null(desired_power)) .dmar_check_desired_power(desired_power)
  if (!is.null(specified_n_per_cell)) {
    if (!is.numeric(specified_n_per_cell) || length(specified_n_per_cell) != 1L ||
        is.na(specified_n_per_cell) || specified_n_per_cell < 2 || specified_n_per_cell != round(specified_n_per_cell)) {
      .dmar_refuse_argument(
        "'specified_n_per_cell' must be a single whole number of at least 2.",
        arg = "specified_n_per_cell")
    }
    specified_n_per_cell <- as.integer(specified_n_per_cell)
    res <- evaluate_at(specified_n_per_cell)
    if (res$design$df_error < 1) {
      .dmar_refuse_argument(c(
        "'specified_n_per_cell' leaves no residual degrees of freedom for this model ",
        "(it has 2 * cells parameters); increase it."), arg = "specified_n_per_cell")
    }
    return(.finish_factorial_het(
      size_term = "specified_n_per_cell", size = specified_n_per_cell,
      res = res, cells = cells,
      alpha_level = alpha_level, desired_power = NULL,
      pop = if (from_pop) list(true_means = true_means,
                               true_cor = true_cor,
                               true_sd = true_sd,
                               true_sd_X = true_sd_X) else NULL,
      approximate = approximate))
  }

  .dmar_check_desired_power(desired_power)
  n_i <- max(2L, as.integer(ceiling((2 * cells + 1) / cells)) + 1L)
  res <- evaluate_at(n_i)
  while (is.na(res$composite) || res$composite < desired_power) {
    n_i <- n_i + 1L
    res <- evaluate_at(n_i)
    if (n_i > 1e6) {
      .dmar_refuse(
        "Failed to converge within a reasonable sample size.",
        "computation_failed")
    }
  }
  .finish_factorial_het(
    size_term = "necessary_n_per_cell", size = n_i,
    res = res, cells = cells,
    alpha_level = alpha_level, desired_power = desired_power,
    pop = if (from_pop) list(true_means = true_means,
                             true_cor = true_cor,
                             true_sd = true_sd,
                             true_sd_X = true_sd_X) else NULL,
    approximate = approximate)
}


# Assemble the tidy table for the heterogeneous-slope planner.
#
# When the cell correlations leave the cells with unequal residual variances,
# the rows carrying a power, and a sample size resolved against one, are renamed
# with an approximate_ prefix. The numbers are the same either way; the names
# are what stop the table from claiming more than the method delivers.
.finish_factorial_het <- function(size_term, size, res, cells,
                                  alpha_level, desired_power, pop,
                                  approximate = FALSE) {
  d      <- res$design
  labels <- d$label
  out <- data.frame(
    term = c(size_term, "total_N", "composite_power", "df_error",
             paste0("power_", labels), paste0("true_f_", labels),
             paste0("df_", labels), paste0("noncentral_f_parm_", labels),
             "cells", "alpha_level"),
    value = c(size, d$N, res$composite, d$df_error,
              unname(res$marginal), d$f, d$df_num, unname(d$lambda),
              cells, alpha_level),
    stringsAsFactors = FALSE)
  if (!is.null(desired_power)) {
    out <- rbind(out, data.frame(term = "desired_power", value = desired_power))
  }
  if (approximate) out$term <- .composite_approximate_terms(out$term)
  class(out) <- c("dmar_composite_power_factorial_het", "dmar_ss_power",
                  "data.frame")
  attr(out, "composite_terms") <- labels
  attr(out, "effect_types")    <- d$type
  attr(out, "approximate")     <- approximate
  if (!is.null(pop)) {
    attr(out, "true_means")        <- pop$true_means
    attr(out, "true_cor") <- pop$true_cor
    attr(out, "true_sd")        <- pop$true_sd
    attr(out, "true_sd_X")      <- pop$true_sd_X
  }
  .as_dmar_tbl(out)
}


#' @describeIn ss_power_composite_factorial_ancova_het Draw the purported
#'   population values: the per-cell regression lines when population values were
#'   supplied, or the effect size lollipop otherwise.
#'
#' @param x An object returned by \code{ss_power_composite_factorial_ancova_het}.
#' @param \dots Further arguments to the figure: \code{palette}, \code{title},
#'   and, for the regression-line figure, \code{cov_range} (covariate range in
#'   SDs either side of the mean, default 2).
#'
#' @export
plot.dmar_composite_power_factorial_het <- function(x, ...) {
  dots <- list(...)
  arg  <- function(nm, default) if (is.null(dots[[nm]])) default else dots[[nm]]
  # A table planned under unequal residual variances carries approximate_ row
  # names. Restoring the exact names here lets the lookups below name one thing.
  x    <- .composite_as_exact(x)
  val  <- function(term) { v <- x$value[x$term == term]
    if (!length(v)) NA_real_ else v[1L] }
  labels <- attr(x, "composite_terms")
  n <- val("specified_n_per_cell"); if (is.na(n)) n <- val("necessary_n_per_cell")
  N <- val("total_N")

  true_means <- attr(x, "true_means")
  if (!is.null(true_means)) {
    return(.plot_factorial_het(
      true_means, attr(x, "true_cor"), attr(x, "true_sd"),
      attr(x, "true_sd_X"),
      factor_levels = dim(true_means), composite = val("composite_power"),
      n_per_cell = n, N = N, alpha_level = val("alpha_level"),
      palette = arg("palette", "okabe_ito"), cov_range = arg("cov_range", 2),
      title = arg("title", "Purported Population Regression Lines by Cell")))
  }

  effects_df <- data.frame(
    label          = labels,
    f              = vapply(labels, function(l) val(paste0("true_f_", l)), numeric(1)),
    df_num         = vapply(labels, function(l) val(paste0("df_", l)), numeric(1)),
    marginal_power = vapply(labels, function(l) val(paste0("power_", l)),
                            numeric(1)),
    stringsAsFactors = FALSE)
  .plot_factorial_composite(
    effects_df, composite = val("composite_power"), n_per_cell = n, N = N,
    alpha_level = val("alpha_level"), palette = arg("palette", "okabe_ito"),
    title = arg("title", "Purported Population Effects the Design Is Planned On"))
}
