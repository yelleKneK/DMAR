#' Sample Size Planning for Power in Factorial ANCOVA
#'
#' Power and sample size for any effect (a main effect or any interaction) in
#' a between-subjects factorial design with covariates: the analysis of
#' covariance generalization of \code{\link{ss_power_factorial_anova}}.
#' Covariates earn their keep by absorbing error variance: with a joint
#' squared multiple correlation \eqn{R^2} between the covariates and the
#' outcome within cells, the error variance falls by the factor
#' \eqn{1 - R^2}, so an effect size \eqn{f} defined on the original (ANOVA)
#' metric grows to \eqn{f / \sqrt{1 - R^2}} in the covariate-adjusted
#' analysis, at the price of one error degree of freedom per covariate.
#'
#' @param factor_levels Integer vector giving the number of levels of each
#'   factor, for example \code{c(2, 4, 3)} for a 2 x 4 x 3 design.
#' @param effect_indices Integer vector identifying the factors that define
#'   the effect of interest: \code{1} for the first factor's main effect,
#'   \code{c(2, 3)} for the B x C interaction, and so on.
#' @param f Cohen's \eqn{f} for the chosen effect \emph{on the unadjusted
#'   (ANOVA) metric}, that is, with the within-cell standard deviation of
#'   the outcome in its denominator. Supply this or
#'   \code{partial_eta_squared}, not both.
#' @param partial_eta_squared Partial eta squared for the chosen effect on
#'   the unadjusted metric.
#' @param covariate_R2 Joint squared multiple correlation between the
#'   covariates and the outcome within cells, in \eqn{[0, 1)}. \code{0}
#'   reproduces the ANOVA analysis (with the covariate degrees of freedom
#'   still spent if \code{n_covariates > 0}).
#' @param n_covariates Number of covariates, a non-negative integer.
#' @param desired_power Desired power; the per-cell sample size is solved
#'   when \code{n_per_cell} is \code{NULL}. Defaults to 0.85 to match
#'   \code{\link{ss_power_factorial_anova}}.
#' @param alpha_level Type I error rate.
#' @param n_per_cell Per-cell sample size; when supplied, the realized power
#'   at that size is returned instead of solving for size.
#'
#' @details
#' The test of an effect with numerator degrees of freedom
#' \eqn{\mathit{df}_h} (the product of the involved factors' levels each
#' minus one) is a noncentral \emph{F} with noncentrality
#' \eqn{\lambda = N f_{\mathrm{adj}}^2}, where \eqn{N} is the total sample
#' size, \eqn{f_{\mathrm{adj}} = f / \sqrt{1 - R^2}}, and error degrees of
#' freedom \eqn{N - (\prod \mathrm{levels}) - q} for \eqn{q} covariates
#' (the standard one-line ANCOVA adjustment; Maxwell, Delaney, & Kelley,
#' 2027, Chapter 9). The covariate slopes are assumed homogeneous across
#' cells and the covariates measured at baseline, so that adjusting does not
#' bias the treatment effects in a randomized design.
#'
#' The complete worked example for this function, a 2 x 4 x 3 ANCOVA with
#' two baseline covariates, planned effect by effect and then analyzed with
#' Type III sums of squares, interaction plots, and focused follow-up
#' contrasts, is the \dQuote{Power for factorial ANCOVA} vignette:
#' \code{vignette("ancova_2x4x3_power", package = "DMAR")}.
#'
#' @return A \code{data.frame} with the
#'   per-cell and total sample sizes (or the supplied ones), the realized
#'   \code{actual_power}, the numerator and error degrees of freedom, the
#'   unadjusted and covariate-adjusted effect sizes (\code{f},
#'   \code{f_adjusted}), \code{covariate_R2}, \code{n_covariates}, the
#'   noncentrality parameter, and \code{alpha_level}. The result carries the
#'   \code{dmar_ss_power} class, so \code{\link[generics]{tidy}} and
#'   \code{\link[generics]{glance}} summarize it in broom convention (the
#'   reported size is the per-cell count).
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9 on designs with covariates.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_factorial_anova}} for the no-covariate
#'   case this wraps; \code{\link{ancova}} and \code{\link{ci_sc_ancova}}
#'   for the analysis side;
#'   \code{vignette("ancova_2x4x3_power", package = "DMAR")} for the full
#'   worked design.
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @keywords design
#'
#' @examples
#' # A 2 x 4 x 3 design, planning the f = .10 main effect of the
#' # two-level factor A. Two baseline covariates with a modest joint
#' # R^2 = .25 cut the required total N by roughly a quarter:
#' ss_power_factorial_ancova(factor_levels = c(2, 4, 3), effect_indices = 1,
#'                           f = 0.10, covariate_R2 = 0,    n_covariates = 0,
#'                           desired_power = 0.80)
#' ss_power_factorial_ancova(factor_levels = c(2, 4, 3), effect_indices = 1,
#'                           f = 0.10, covariate_R2 = 0.25, n_covariates = 2,
#'                           desired_power = 0.80)
#'
#' # Realized power for the three-way interaction at 6 per cell.
#' ss_power_factorial_ancova(factor_levels = c(2, 4, 3),
#'                           effect_indices = c(1, 2, 3), f = 0.15,
#'                           covariate_R2 = 0.25, n_covariates = 2,
#'                           n_per_cell = 6)
#'
#' @export
#' @importFrom stats pf qf
ss_power_factorial_ancova <- function(factor_levels, effect_indices,
                                      f = NULL, partial_eta_squared = NULL,
                                      covariate_R2 = 0, n_covariates = 0,
                                      desired_power = 0.85,
                                      alpha_level = 0.05,
                                      n_per_cell = NULL) {
  if (!is.numeric(factor_levels) || length(factor_levels) < 1 ||
      any(factor_levels < 2) || any(factor_levels != round(factor_levels))) {
    stop("'factor_levels' must be a vector of integers each >= 2.",
         call. = FALSE)
  }
  if (!is.numeric(effect_indices) || length(effect_indices) < 1 ||
      any(effect_indices < 1) ||
      any(effect_indices > length(factor_levels)) ||
      any(duplicated(effect_indices))) {
    stop("'effect_indices' must be unique integer indices into ",
         "'factor_levels'.", call. = FALSE)
  }
  if (is.null(f) == is.null(partial_eta_squared)) {
    stop("Specify exactly one of 'f' or 'partial_eta_squared'.",
         call. = FALSE)
  }
  if (!is.null(partial_eta_squared)) {
    if (!is.numeric(partial_eta_squared) ||
        length(partial_eta_squared) != 1L || is.na(partial_eta_squared) ||
        partial_eta_squared <= 0 || partial_eta_squared >= 1) {
      stop("'partial_eta_squared' must be a single number in (0, 1).",
           call. = FALSE)
    }
    f <- sqrt(partial_eta_squared / (1 - partial_eta_squared))
  }
  if (!is.numeric(f) || length(f) != 1L || is.na(f) || f <= 0) {
    stop("'f' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(covariate_R2) || length(covariate_R2) != 1L ||
      is.na(covariate_R2) || covariate_R2 < 0 || covariate_R2 >= 1) {
    stop("'covariate_R2' must be a single number in [0, 1).", call. = FALSE)
  }
  if (!is.numeric(n_covariates) || length(n_covariates) != 1L ||
      is.na(n_covariates) || n_covariates < 0 ||
      n_covariates != round(n_covariates)) {
    stop("'n_covariates' must be a single non-negative integer.",
         call. = FALSE)
  }
  if (covariate_R2 > 0 && n_covariates == 0) {
    stop("'covariate_R2' is positive but 'n_covariates' is 0; say how many ",
         "covariates carry that R2.", call. = FALSE)
  }
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L ||
      is.na(alpha_level) || alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be a single number in (0, 1).", call. = FALSE)
  }

  cells <- prod(factor_levels)
  df_h  <- prod(factor_levels[effect_indices] - 1)
  f_adj <- f / sqrt(1 - covariate_R2)

  power_at <- function(n_cell) {
    N    <- n_cell * cells
    df_e <- N - cells - n_covariates
    if (df_e < 1) return(0)
    lambda <- N * f_adj^2
    crit   <- qf(1 - alpha_level, df_h, df_e)
    1 - pf(crit, df_h, df_e, ncp = lambda)
  }

  if (is.null(n_per_cell)) {
    if (!is.numeric(desired_power) || length(desired_power) != 1L ||
        is.na(desired_power) || desired_power <= 0 || desired_power >= 1) {
      stop("'desired_power' must be a single number in (0, 1).",
           call. = FALSE)
    }
    n_cell <- max(2L, ceiling((n_covariates + 1) / cells) + 1L)
    # Fail fast when the target cannot be reached at any realistic per-cell
    # size (for example a near-zero effect) rather than incrementing without
    # bound. Power is monotone in the per-cell size.
    n_max <- 1e7
    if (power_at(n_max) < desired_power)
      stop("Could not reach 'desired_power' at a per-cell size up to ",
           format(n_max, scientific = FALSE), "; the effect may be too small ",
           "to detect at this power.", call. = FALSE)
    while (power_at(n_cell) < desired_power) {
      n_cell <- n_cell + 1L
      if (n_cell > n_max)
        stop("The search for 'desired_power' did not converge below ",
             format(n_max, scientific = FALSE), " per cell.", call. = FALSE)
    }
  } else {
    if (!is.numeric(n_per_cell) || length(n_per_cell) != 1L ||
        is.na(n_per_cell) || n_per_cell < 2 ||
        n_per_cell != round(n_per_cell)) {
      stop("'n_per_cell' must be a single integer of at least 2.",
           call. = FALSE)
    }
    n_cell <- as.integer(n_per_cell)
  }

  N <- n_cell * cells
  # The size row is named for its role, matching ss_power_factorial_anova:
  # a solved size is the planner's answer (necessary_), a user-supplied one
  # is an echo (specified_).
  size_term <- if (is.null(n_per_cell)) "necessary_n_per_cell" else "specified_n_per_cell"
  out <- data.frame(
    term  = c(size_term, "total_N", "actual_power", "df_effect",
              "df_error", "f", "f_adjusted", "covariate_R2", "n_covariates",
              "noncentrality", "alpha_level"),
    value = c(n_cell, N, power_at(n_cell), df_h,
              N - cells - n_covariates, f, f_adj, covariate_R2,
              n_covariates, N * f_adj^2, alpha_level),
    stringsAsFactors = FALSE
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
