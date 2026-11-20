#' Sample Size or Composite Power for a Factorial ANCOVA
#'
#' Determine the necessary per-cell sample size to achieve a desired level of
#' composite statistical power in a balanced factorial analysis of covariance,
#' or, given a per-cell sample size, return the realized composite power.
#' Composite power is the probability that every effect named in \code{effects}
#' is statistically significant in the same study, the quantity a design must be
#' planned against when its conclusion requires more than one result to hold at
#' once. Each effect is a main effect or an interaction of the factorial design,
#' tested by its own \emph{F} test, and any subset of them can make up the
#' composite.
#'
#' The population effects can be stated two ways: as effect sizes, a Cohen's
#' \emph{f} or partial eta squared per effect, or as a full array of population
#' cell means together with a common within-cell standard deviation, from which
#' each named effect's \emph{f} is read off the analysis of variance
#' decomposition of the means. Supplying means lets the \code{plot()} method draw
#' the mean pattern itself, so the size of the population effects is on the page.
#'
#' With no covariate (the defaults \code{true_covariate_R2 = 0} and
#' \code{n_covariates = 0}) this is a factorial ANOVA; the wrapper
#' \code{\link{ss_power_composite_factorial_anova}} is that case named directly.
#'
#' @param factor_levels Integer vector of the number of levels of each factor,
#'   one entry per factor (each at least 2). A \code{c(2, 3, 2)} argument is a
#'   2 by 3 by 2 design.
#' @param effects A non-empty list naming the effects in the composite. Each
#'   element is itself a list with \code{factors} (a vector of factor indices
#'   into \code{factor_levels}: one index for a main effect, several for an
#'   interaction) and exactly one of \code{true_f} (Cohen's \emph{f} for that effect)
#'   or \code{true_eta_squared_partial}. An optional \code{label} names the effect in
#'   the output and the figure; the default label is the factor indices joined
#'   by \code{x} (for example \code{"1x2"} for the interaction of factors 1 and
#'   2). Each element also carries the effect size, one of \code{true_f} or
#'   \code{true_eta_squared_partial}, unless \code{true_means} is supplied, in which case
#'   the sizes come from the means and no effect size is given here. The
#'   purported effect sizes are population values the researcher posits, never
#'   sample estimates, which is why the two fields carry \code{true_}. Every
#'   field of an element is named, and a field under any other name is
#'   refused with an error of class \code{dmar_argument_invalid}, so that an
#'   effect size is never dropped without a word.
#' @param true_means Optional array of population cell means whose dimensions are
#'   \code{factor_levels} (a matrix for two factors), or a numeric vector of
#'   length \code{prod(factor_levels)} in array order (the first factor varying
#'   fastest). When supplied, each named effect's Cohen's \emph{f} is computed
#'   from the means and \code{true_sd}, and \code{plot()} draws the mean pattern.
#'   The means are the population values the researcher posits, on the raw scale
#'   of the outcome.
#' @param true_sd The common within-cell population standard deviation of the
#'   outcome (the square root of the error variance), required with \code{true_means}
#'   and used only there. Cohen's \emph{f} for an effect is the spread of its
#'   cell-mean component relative to \code{true_sd}.
#' @param true_covariate_R2 Proportion of the outcome's within-cell variance the
#'   covariate or covariates explain, in \eqn{[0, 1)}. The covariate removes
#'   that fraction of the error variance, which raises every effect's
#'   noncentrality through \eqn{f / \sqrt{1 - R^2}}. Defaults to 0.
#' @param n_covariates Number of covariates, a non-negative integer. Each spends
#'   one residual degree of freedom. Must be positive when \code{true_covariate_R2}
#'   is. Defaults to 0.
#' @param desired_power Desired composite statistical power (default 0.85). Used
#'   only when \code{specified_n_per_cell} is \code{NULL}.
#' @param alpha_level Type I error rate for each individual \emph{F} test
#'   (default 0.05). This is the per-test rate, not a rate for the composite
#'   event.
#' @param specified_n_per_cell Per-cell sample size, assumed balanced across cells; if
#'   supplied, the realized composite power is returned rather than a sample
#'   size planned.
#'
#' @details
#' In a balanced factorial design the effect sums of squares are mutually
#' orthogonal, so the effects are uncorrelated. Orthogonal effects do not give
#' independent tests: every \emph{F} test divides by the same error mean square,
#' so an error estimate that lands low inflates all of the test statistics
#' together. The tests are positively dependent, and composite power is strictly
#' larger than the product of the marginal powers, bounded above by the least
#' powerful test in the set, so the weakest effect governs the design.
#'
#' Conditional on the error estimate the tests are independent, which reduces the
#' composite to a one dimensional integral over the chi square distribution of
#' that estimate. Effect \eqn{j} has numerator degrees of freedom
#' \eqn{\prod (a - 1)} over the factors it spans and noncentrality
#' \eqn{N f_{\mathrm{adj}}^2} with \eqn{f_{\mathrm{adj}} = f / \sqrt{1 - R^2}} and
#' \eqn{N} the total sample size, the same convention and covariate adjustment
#' \code{\link{ss_power_factorial_ancova}} uses. The residual degrees of freedom
#' are \eqn{N - \mathrm{cells} - \mathrm{covariates}}. Adaptive quadrature
#' evaluates the integral, so nothing is simulated and the result is
#' deterministic to quadrature precision. A single-effect composite reproduces
#' the ordinary noncentral \emph{F} power, and naming one effect reproduces
#' \code{\link{ss_power_factorial_anova}} (or \code{\link{ss_power_factorial_ancova}}
#' with a covariate) exactly.
#'
#' \strong{Exactness.} With no covariate the composite is exact to quadrature
#' precision: the balanced factorial \emph{F} tests are exactly noncentral
#' \emph{F} and exactly independent given the error estimate. A covariate
#' introduces the one approximation \code{\link{ss_power_factorial_ancova}}
#' already carries, treating \code{true_covariate_R2} as a fixed reduction of the
#' error variance rather than an estimated one; the departure is of order
#' \eqn{1 / N} and is negligible at the sample sizes a multi-effect composite
#' usually needs.
#'
#' Because every test divides by the same error estimate, composite power is not
#' strictly monotone in \code{specified_n_per_cell} at the smallest residual degrees of
#' freedom. \code{necessary_n_per_cell} is the smallest per-cell size attaining
#' \code{desired_power}.
#'
#' @section The figure:
#' The \code{plot()} method draws the purported population values. When
#' \code{true_means} were supplied it draws the mean pattern itself: a profile of the
#' population cell means over the first factor, one line per level of the second,
#' faceted by any further factors, with an error bar of plus or minus one
#' within-cell standard deviation at each mean so the effect sizes read against
#' the noise. When effect sizes were supplied instead, the cell means are not
#' pinned (many mean patterns share one Cohen's \emph{f}), so it draws the effect
#' sizes: one lollipop per named effect at its partial eta squared, colored and
#' labeled by that effect's marginal power. Either way the composite power is in
#' the subtitle and nothing is simulated. Requires \pkg{ggplot2}.
#'
#' @return A \code{data.frame} with \code{term} and \code{value} columns: the
#'   \code{necessary_n_per_cell} (or the supplied \code{specified_n_per_cell})
#'   and the total sample size it implies, \code{total_N}, the
#'   \code{composite_power}, the \code{df_error}, then for each effect its
#'   marginal \code{power_<label>}, purported \code{true_f_<label>}, numerator
#'   \code{df_<label>}, and \code{noncentral_f_parm_<label>}, followed by rows
#'   echoing \code{true_covariate_R2}, \code{n_covariates}, \code{cells}, and
#'   \code{alpha_level}. The result carries the \code{dmar_ss_power} class, so
#'   \code{\link[generics]{tidy}} and \code{\link[generics]{glance}} summarize
#'   the per-cell size and the composite power in broom convention, and a
#'   \code{dmar_composite_power_factorial} class so \code{plot()} draws the
#'   figure.
#'
#' @references
#' Maxwell, S. E. (2004). The persistence of underpowered studies in
#'   psychological research: Causes, consequences, and remedies.
#'   \emph{Psychological Methods, 9}, 147--163.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective} (4th ed.).
#'   Routledge. (See Chapter 7 on factorial designs, Chapter 9 on the analysis
#'   of covariance, and Chapter 3 on statistical power.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_composite_factorial_anova}} for the no-covariate
#'   case; \code{\link{ss_power_composite_ancova_2group}} for the two-group ANCOVA
#'   composite of the group effect, the covariate effect, and their interaction;
#'   \code{\link{ss_power_factorial_ancova}} and
#'   \code{\link{ss_power_factorial_anova}} for a single effect
#'
#' @examples
#' # A 2 by 3 factorial ANCOVA. The conclusion needs both main effects to hold,
#' # so the design is planned against their composite. The covariate explains
#' # 25 percent of the within-cell variance (one covariate).
#' ss_power_composite_factorial_ancova(
#'   factor_levels = c(2, 3),
#'   effects = list(list(factors = 1, true_f = 0.25),
#'                  list(factors = 2, true_f = 0.20)),
#'   true_covariate_R2 = 0.25, n_covariates = 1,
#'   desired_power = 0.80)
#'
#' # Realized composite power at 40 per cell for a main effect and the
#' # interaction of a 2 by 2 design, effect sizes given as partial eta squared.
#' ss_power_composite_factorial_ancova(
#'   factor_levels = c(2, 2),
#'   effects = list(list(factors = 1,        true_eta_squared_partial = 0.06),
#'                  list(factors = c(1, 2),  true_eta_squared_partial = 0.04)),
#'   specified_n_per_cell = 40)
#'
#' # Naming one effect reproduces the single-effect planner exactly.
#' ss_power_composite_factorial_ancova(
#'   factor_levels = c(2, 3), effects = list(list(factors = 2, true_f = 0.25)),
#'   specified_n_per_cell = 20)
#' ss_power_factorial_anova(factor_levels = c(2, 3), effect_indices = 2,
#'                          true_f = 0.25, specified_n_per_cell = 20)
#'
#' # The composite is not the product of the marginal powers. The tests share
#' # one error estimate, so they are positively dependent and the composite is
#' # the larger of the two.
#' plan <- ss_power_composite_factorial_ancova(
#'   factor_levels = c(2, 2, 3),
#'   effects = list(list(factors = 1, true_f = 0.30, label = "A"),
#'                  list(factors = c(1, 3), true_f = 0.25, label = "AxC")),
#'   specified_n_per_cell = 15)
#' plan$value[plan$term == "composite_power"]
#' prod(plan$value[plan$term %in% c("power_A", "power_AxC")])
#'
#' # A one-row broom summary of the plan.
#' generics::tidy(ss_power_composite_factorial_ancova(
#'   factor_levels = c(2, 3),
#'   effects = list(list(factors = 1, true_f = 0.25),
#'                  list(factors = 2, true_f = 0.20)),
#'   desired_power = 0.80))
#'
#' # The figure of the purported population values, annotated with the power a
#' # given sample size delivers.
#' plot(ss_power_composite_factorial_ancova(
#'   factor_levels = c(2, 3),
#'   effects = list(list(factors = 1, true_f = 0.25),
#'                  list(factors = 2, true_f = 0.20)),
#'   specified_n_per_cell = 30))
#'
#' # The effects can instead be stated as a full pattern of population cell means
#' # with a common within-cell SD. Rows are the 2-level factor, columns the
#' # 3-level factor. plot() then draws the mean pattern itself.
#' cell_means <- matrix(c(10, 12, 11,
#'                        13, 12, 16), nrow = 2, byrow = TRUE)
#' fit <- ss_power_composite_factorial_ancova(
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
ss_power_composite_factorial_ancova <- function(factor_levels, effects,
                                                true_means = NULL, true_sd = NULL,
                                                true_covariate_R2 = 0,
                                                n_covariates = 0,
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
  if (!is.numeric(true_covariate_R2) || length(true_covariate_R2) != 1L ||
      is.na(true_covariate_R2) || true_covariate_R2 < 0 || true_covariate_R2 >= 1) {
    .dmar_refuse_argument(
      "'true_covariate_R2' must be a single number in [0, 1).",
      arg = "true_covariate_R2")
  }
  if (!is.numeric(n_covariates) || length(n_covariates) != 1L ||
      is.na(n_covariates) || n_covariates < 0 ||
      n_covariates != round(n_covariates)) {
    .dmar_refuse_argument(
      "'n_covariates' must be a single non-negative integer.",
      arg = "n_covariates")
  }
  if (true_covariate_R2 > 0 && n_covariates == 0) {
    .dmar_refuse_argument(c(
      "'true_covariate_R2' is positive but 'n_covariates' is 0; say how many ",
      "covariates carry that variance."),
      arg = c("true_covariate_R2", "n_covariates"))
  }
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L ||
      is.na(alpha_level) || alpha_level <= 0 || alpha_level >= 1) {
    .dmar_refuse_argument(
      "'alpha_level' must be a single number in (0, 1).", arg = "alpha_level")
  }
  n_covariates <- as.integer(n_covariates)

  # Two ways to state the effects: their sizes directly in 'effects', or a full
  # array of population cell means with a common within-cell 'true_sd', from which
  # each named effect's Cohen's f is read off the orthogonal ANOVA decomposition
  # of the means.
  from_means <- !is.null(true_means)
  if (from_means) {
    if (is.null(true_sd) || !is.numeric(true_sd) || length(true_sd) != 1L ||
        is.na(true_sd) || true_sd <= 0) {
      .dmar_refuse_argument(c(
        "With 'true_means', supply a single positive 'true_sd', the common ",
        "within-cell standard deviation."), arg = c("true_means", "true_sd"))
    }
    true_means <- .as_factorial_means(true_means, factor_levels)
    effects <- .parse_factorial_effects(effects, length(factor_levels),
                                        require_f = FALSE)
    effects <- lapply(effects, function(e) {
      e$f <- .factorial_effect_f(true_means, e$factors, true_sd); e
    })
    zero <- vapply(effects, function(e) e$f < 1e-9, logical(1))
    if (any(zero)) {
      .dmar_refuse_argument(c(
        "These effects do not vary in the supplied means (Cohen's f is 0), ",
        "so no sample size attains a power above the Type I error rate: ",
        paste(vapply(effects[zero], function(e) e$label, character(1)),
              collapse = ", "),
        ". Drop them from 'effects', or supply means in which they vary."),
        arg = "effects")
    }
  } else {
    if (!is.null(true_sd)) {
      .dmar_refuse_argument(c(
        "'true_sd' is only used with 'true_means'. Supply the effect sizes in ",
        "'effects', or supply both 'true_means' and 'true_sd'."),
        arg = c("true_sd", "true_means", "effects"))
    }
    effects <- .parse_factorial_effects(effects, length(factor_levels),
                                        require_f = TRUE)
  }

  cells <- prod(factor_levels)

  evaluate_at <- function(n) {
    d   <- .factorial_composite_design(factor_levels, effects, true_covariate_R2,
                                       n_covariates, n)
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
        "'specified_n_per_cell' leaves no residual degrees of freedom for this model; ",
        "increase it."), arg = "specified_n_per_cell")
    }
    return(.finish_factorial_composite(
      size_term = "specified_n_per_cell", size = specified_n_per_cell,
      res = res, effects = effects,
      true_covariate_R2 = true_covariate_R2, n_covariates = n_covariates,
      cells = cells, alpha_level = alpha_level, desired_power = NULL,
      true_means = if (from_means) true_means else NULL,
      true_sd = if (from_means) true_sd else NULL))
  }

  .dmar_check_desired_power(desired_power)
  # Start where the model has at least one residual degree of freedom, then walk
  # up to the smallest per-cell size attaining the target. Every effect has a
  # positive f, so each marginal power tends to 1 and the composite does too,
  # and the search terminates.
  n_i <- max(2L, as.integer(ceiling((cells + n_covariates + 1) / cells)) + 1L)
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
  .finish_factorial_composite(
    size_term = "necessary_n_per_cell", size = n_i,
    res = res, effects = effects,
    true_covariate_R2 = true_covariate_R2, n_covariates = n_covariates,
    cells = cells, alpha_level = alpha_level, desired_power = desired_power,
    true_means = if (from_means) true_means else NULL,
    true_sd = if (from_means) true_sd else NULL)
}


# Assemble the tidy table: the design result, then each effect's marginal power,
# purported f, numerator df, and noncentrality, then the planning echoes. The
# per-effect labels come from the parsed effects. dmar_composite_power_factorial
# leads so plot() dispatches to the factorial figure; dmar_ss_power brings the
# shared tidy()/glance(), which read composite_power and the per-cell size.
.finish_factorial_composite <- function(size_term, size, res, effects,
                                        true_covariate_R2, n_covariates, cells,
                                        alpha_level, desired_power,
                                        true_means = NULL, true_sd = NULL) {
  d      <- res$design
  labels <- d$label
  out <- data.frame(
    term = c(size_term, "total_N", "composite_power", "df_error",
             paste0("power_", labels),
             paste0("true_f_", labels),
             paste0("df_", labels),
             paste0("noncentral_f_parm_", labels),
             "true_covariate_R2", "n_covariates", "cells", "alpha_level"),
    value = c(size, d$N, res$composite, d$df_error,
              unname(res$marginal),
              d$f,
              d$df_num,
              unname(d$lambda),
              true_covariate_R2, n_covariates, cells, alpha_level),
    stringsAsFactors = FALSE)
  if (!is.null(desired_power)) {
    out <- rbind(out, data.frame(term = "desired_power", value = desired_power))
  }
  class(out) <- c("dmar_composite_power_factorial", "dmar_ss_power",
                  "data.frame")
  attr(out, "composite_terms") <- labels
  # When the effects came from cell means, carry the means and the within-cell
  # SD so plot() can draw the mean pattern rather than the effect sizes.
  if (!is.null(true_means)) {
    attr(out, "true_means") <- true_means
    attr(out, "true_sd") <- true_sd
  }
  .as_dmar_tbl(out)
}


#' @describeIn ss_power_composite_factorial_ancova Draw the purported population
#'   effect sizes a result was planned on, each annotated with its marginal
#'   power, and the composite power in the subtitle.
#'
#' @param x An object returned by \code{ss_power_composite_factorial_ancova} or
#'   \code{\link{ss_power_composite_factorial_anova}}.
#' @param \dots Further arguments to the figure: \code{palette} (a palette
#'   name, default \code{"okabe_ito"}) and \code{title}.
#'
#' @export
plot.dmar_composite_power_factorial <- function(x, ...) {
  dots <- list(...)
  arg  <- function(nm, default) if (is.null(dots[[nm]])) default else dots[[nm]]

  val <- function(term) {
    v <- x$value[x$term == term]
    if (!length(v)) NA_real_ else v[1L]
  }
  labels <- attr(x, "composite_terms")
  n <- val("specified_n_per_cell")
  if (is.na(n)) n <- val("necessary_n_per_cell")
  N <- val("total_N")

  true_means <- attr(x, "true_means")
  if (!is.null(true_means)) {
    # Cell means were supplied: show the population mean pattern itself.
    return(.plot_factorial_means(
      true_means, factor_levels = dim(true_means), composite = val("composite_power"),
      n_per_cell = n, N = N, true_sd = attr(x, "true_sd"),
      alpha_level = val("alpha_level"), palette = arg("palette", "okabe_ito"),
      title = arg("title", "Purported Population Cell Means")))
  }

  effects_df <- data.frame(
    label          = labels,
    f              = vapply(labels, function(l) val(paste0("true_f_", l)), numeric(1)),
    df_num         = vapply(labels, function(l) val(paste0("df_", l)), numeric(1)),
    marginal_power = vapply(labels, function(l) val(paste0("power_", l)),
                            numeric(1)),
    stringsAsFactors = FALSE)

  .plot_factorial_composite(
    effects_df, composite = val("composite_power"), n_per_cell = n,
    N = N, alpha_level = val("alpha_level"),
    palette = arg("palette", "okabe_ito"),
    title = arg("title", "Purported Population Effects the Design Is Planned On"))
}
