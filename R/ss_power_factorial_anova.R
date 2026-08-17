#' Sample Size or Power for a Factorial Between-Subjects ANOVA Effect
#'
#' Determine the necessary per-cell sample size to achieve a desired level of statistical power
#' for a single \emph{F} test (main effect or interaction) in a between-subjects factorial ANOVA, or,
#' given a per-cell sample size, return the realized statistical power. The function handles
#' two-way and higher-order factorial designs.
#'
#' @param factor_levels Integer vector giving the number of levels of each factor (e.g., \code{c(2, 3)} for a 2 x 3 design, \code{c(2, 2, 4)} for a 2 x 2 x 4 design)
#' @param effect_indices Integer vector identifying the factors that define the effect of interest. For example, \code{1} requests the main effect of the first factor; \code{c(1, 2)} requests the AxB two-way interaction; \code{c(1, 2, 3)} requests the three-way interaction
#' @param f Cohen's \emph{f} effect size for the chosen effect; supply this or \code{partial_eta_squared}, but not both
#' @param partial_eta_squared Partial eta squared for the chosen effect; supply this or \code{f}
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param n_per_cell Per-cell sample size; if specified, returns the realized power
#'
#' @details
#' For a between-subjects factorial design, the \emph{F} statistic for the chosen effect follows a
#' noncentral \emph{F} distribution under the alternative with numerator degrees of freedom
#' \eqn{\prod_{i \in S} (k_i - 1)} (where \eqn{S} is \code{effect_indices} and \eqn{k_i} is
#' \code{factor_levels[i]}), denominator degrees of freedom \eqn{N - K} (where
#' \eqn{N = n_{cell} \prod k_i} and \eqn{K = \prod k_i} is the number of cells), and noncentrality
#' parameter \eqn{\lambda = N f^2}. Cohen's \emph{f} relates to partial eta squared via
#' \eqn{f = \sqrt{\eta_p^2 / (1 - \eta_p^2)}}.
#'
#' The function searches over per-cell sample sizes until power reaches \code{desired_power}; when
#' \code{n_per_cell} is supplied it returns the realized power.
#'
#' @return A \code{data.frame} with rows for \code{necessary_n_per_cell} (or \code{specified_n_per_cell}),
#'   \code{total_N}, \code{effect_df}, \code{error_df}, \code{noncentrality}, and \code{actual_power}.
#'   The result carries the \code{dmar_ss_power} class, so
#'   \code{\link[generics]{tidy}} and \code{\link[generics]{glance}} summarize
#'   it in broom convention (the reported size is the per-cell count).
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing experiments and analyzing data: A model comparison perspective} (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_one_way_anova}}, \code{\link{ss_power_c}}, \code{\link{conf_limits_ncf}}
#'
#' @examples
#' # 2 x 3 design, main effect of factor B (the 3-level factor), f = 0.25, power = .80
#' ss_power_factorial_anova(factor_levels = c(2, 3), effect_indices = 2,
#'                          f = 0.25, desired_power = 0.80)
#'
#' # 2 x 2 design, AxB interaction, partial eta squared = 0.06, power = .80
#' ss_power_factorial_anova(factor_levels = c(2, 2), effect_indices = c(1, 2),
#'                          partial_eta_squared = 0.06, desired_power = 0.80)
#'
#' # 2 x 2 x 3 design, three-way interaction, f = 0.20
#' ss_power_factorial_anova(factor_levels = c(2, 2, 3), effect_indices = c(1, 2, 3),
#'                          f = 0.20, desired_power = 0.80)
#'
#' # Realized power for n_per_cell = 20 in a 2x3 design, AxB interaction, f = 0.25
#' ss_power_factorial_anova(factor_levels = c(2, 3), effect_indices = c(1, 2),
#'                          f = 0.25, n_per_cell = 20)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @details
#' For covariates, see \code{\link{ss_power_factorial_ancova}}. A complete
#' worked three-factor example (a 2 x 4 x 3 design planned effect by
#' effect, simulated, analyzed with Type III sums of squares, plotted,
#' and followed up with focused contrasts) is the vignette
#' \code{vignette("ancova_2x4x3_power", package = "DMAR")}.
#'
#' @export
ss_power_factorial_anova <- function(factor_levels, effect_indices, f = NULL,
                                     partial_eta_squared = NULL,
                                     desired_power = 0.85, alpha_level = 0.05,
                                     n_per_cell = NULL) {
  if (!is.numeric(factor_levels) || length(factor_levels) < 1 || any(factor_levels < 2) || any(factor_levels != as.integer(factor_levels))) {
    stop("'factor_levels' must be a vector of integers each >= 2.", call. = FALSE)
  }
  if (!is.numeric(effect_indices) || length(effect_indices) < 1 || any(effect_indices < 1) || any(effect_indices > length(factor_levels)) || any(duplicated(effect_indices))) {
    stop("'effect_indices' must be unique integer indices into 'factor_levels'.", call. = FALSE)
  }
  if (!is.null(f) && !is.null(partial_eta_squared)) {
    stop("Specify either 'f' or 'partial_eta_squared', not both.", call. = FALSE)
  }
  if (is.null(f) && is.null(partial_eta_squared)) {
    stop("Specify exactly one of 'f' or 'partial_eta_squared'.", call. = FALSE)
  }
  if (!is.null(partial_eta_squared)) {
    if (partial_eta_squared <= 0 || partial_eta_squared >= 1) {
      stop("'partial_eta_squared' must be in (0, 1).", call. = FALSE)
    }
    f <- sqrt(partial_eta_squared / (1 - partial_eta_squared))
  }
  if (f <= 0) stop("'f' must be positive.", call. = FALSE)
  if (alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be in (0, 1).", call. = FALSE)
  }

  K <- prod(factor_levels)
  effect_df <- prod(factor_levels[effect_indices] - 1)

  power_at <- function(n_cell) {
    N <- n_cell * K
    df_2 <- N - K
    if (df_2 <= 0) return(NA_real_)
    ncp <- N * f^2
    crit <- qf(1 - alpha_level, effect_df, df_2)
    1 - pf(crit, effect_df, df_2, ncp = ncp)
  }

  if (!is.null(n_per_cell)) {
    if (n_per_cell < 2) stop("'n_per_cell' must be at least 2.", call. = FALSE)
    pwr  <- power_at(n_per_cell)
    N    <- n_per_cell * K
    df_2 <- N - K
    out <- data.frame(
      term  = c("specified_n_per_cell", "total_N", "effect_df", "error_df", "noncentrality", "actual_power"),
      value = c(n_per_cell, N, effect_df, df_2, N * f^2, pwr)
    )
    class(out) <- c("dmar_ss_power", class(out))
    return(.as_dmar_tbl(out))
  }

  if (desired_power <= 0 || desired_power >= 1) {
    stop("'desired_power' must be in (0, 1).", call. = FALSE)
  }

  n_i <- 1
  pwr <- 0
  while (pwr < desired_power) {
    n_i <- n_i + 1
    pwr <- power_at(n_i)
    if (n_i > 1e6) stop("Failed to converge within reasonable sample size.", call. = FALSE)
  }
  N <- n_i * K
  out <- data.frame(
    term  = c("necessary_n_per_cell", "total_N", "effect_df", "error_df", "noncentrality", "actual_power"),
    value = c(n_i, N, effect_df, N - K, N * f^2, pwr)
  )
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
