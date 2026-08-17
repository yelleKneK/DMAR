#' Sample Size or Power for a Standardized Mean Difference (Two Independent Groups)
#'
#' Determine the necessary per-group sample size to achieve a desired level of statistical power
#' for the two-sample (independent groups) \emph{t}-test on a standardized mean difference (Cohen's
#' \emph{d}; equivalently Hedges' \emph{g} and Glass's \emph{g} for sample size purposes).
#' Alternatively, given a per-group sample size, return the realized statistical power.
#'
#' @param smd Supposed standardized mean difference (Cohen's \emph{d}) the design is planned
#'   against: a value the researcher posits for the population, either a minimally important
#'   effect or a value believed to be true in the population, never a sample estimate. Echoed
#'   in the returned table as the \code{supposed_smd} row.
#' @param desired_power Desired statistical power (default 0.85)
#' @param alpha_level Type I error rate (default 0.05)
#' @param n_1 Sample size for group 1 (if specified, the function returns the realized power; assumes \code{n_2 = n_1} unless \code{n_2} is also given)
#' @param n_2 Sample size for group 2 (defaults to \code{n_1} when \code{n_1} is supplied)
#' @param directional Logical: \code{TRUE} for a one-sided test (in the same sign as \code{smd}), \code{FALSE} (default) for a two-sided test
#'
#' @details
#' The two-sample \emph{t}-statistic with pooled standard deviation follows a noncentral \emph{t}-distribution with
#' \eqn{n_1 + n_2 - 2} degrees of freedom and noncentrality parameter
#' \eqn{\lambda = \delta \sqrt{n_1 n_2 / (n_1 + n_2)}}, where \eqn{\delta} is the population
#' standardized mean difference. For balanced designs (\eqn{n_1 = n_2 = n}) this simplifies to
#' \eqn{\lambda = \delta \sqrt{n / 2}}.
#'
#' Power is computed as the probability that the absolute value of the test statistic exceeds the
#' critical value(s) under the alternative; the function returns the per-group sample size for which
#' power first reaches \code{desired_power}.
#'
#' Kelley and Rausch (2006) develop the accuracy in parameter estimation
#' approach to planning the sample size for the standardized mean
#' difference, implemented in \code{\link{ss_aipe_smd}}.
#'
#' @return
#' A \code{data.frame} with \code{term} and \code{value} columns. The design
#' result comes first, followed by rows that echo the user-supplied planning
#' inputs, so the assumptions the power was evaluated under travel with the
#' result. The \code{supposed_smd} row is the supposed effect the plan is built
#' on: a value the researcher posits, either a minimally important effect or a
#' value believed to be true in the population, never a sample estimate. The
#' \code{tails} row is 2 for a nondirectional test and 1 for a directional test.
#' \describe{
#'   \item{When \code{n_1} is \code{NULL}}{Result rows \code{necessary_n_per_group},
#'     \code{actual_power}, and \code{noncentral_t_parm}, then the planning inputs
#'     \code{supposed_smd}, \code{desired_power}, \code{alpha_level}, and \code{tails}.}
#'   \item{When \code{n_1} is specified}{Result rows \code{specified_n_1},
#'     \code{specified_n_2}, \code{actual_power}, and \code{noncentral_t_parm}, then
#'     the planning inputs \code{supposed_smd}, \code{alpha_level}, and \code{tails}
#'     (the supplied group sizes are the \code{specified_n_1} / \code{specified_n_2}
#'     rows).}
#' }
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
#'   obtaining precision: Delineating methods of sample size planning.
#'   \emph{Evaluation and the Health Professions, 26}(3), 258--287.
#'   \doi{10.1177/0163278703255242}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing experiments and analyzing data: A model comparison perspective} (4th ed.). Routledge.
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#'   for statistical power and accuracy in parameter estimation.
#'   \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_smd}}, \code{\link{ci_smd}}, \code{\link{smd}}, \code{\link{conf_limits_nct}}
#'
#' @examples
#' # Per-group sample size for d = 0.5, alpha = .05, power = .80, two-sided
#' ss_power_smd(smd = 0.5, desired_power = 0.80)
#'
#' # Same with a directional (one-sided) test
#' ss_power_smd(smd = 0.5, desired_power = 0.80, directional = TRUE)
#'
#' # Realized power given balanced n = 30 per group
#' ss_power_smd(smd = 0.5, n_1 = 30)
#'
#' # Realized power for unbalanced (n_1 = 30, n_2 = 50)
#' ss_power_smd(smd = 0.5, n_1 = 30, n_2 = 50)
#'
#' @keywords design htest
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family sample size for power
#'
#' @export
ss_power_smd <- function(smd, desired_power = 0.85, alpha_level = 0.05,
                         n_1 = NULL, n_2 = NULL, directional = FALSE) {
  if (!is.numeric(smd) || length(smd) != 1L || !is.finite(smd)) {
    stop("'smd' must be a single finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L || alpha_level <= 0 || alpha_level >= 1) {
    stop("'alpha_level' must be a single numeric value in (0, 1).", call. = FALSE)
  }
  if (is.null(n_1)) {
    if (!is.numeric(desired_power) || length(desired_power) != 1L || desired_power <= 0 || desired_power >= 1) {
      stop("'desired_power' must be a single numeric value in (0, 1).", call. = FALSE)
    }
    # A null effect has power equal to alpha_level at every sample size, so the
    # search for an N attaining desired_power would never terminate. Stop cleanly
    # rather than iterate to the sample size cap and report a convergence
    # failure. (Power at a *specified* N with smd = 0 is a valid query, equal to
    # alpha_level, so the specified-N branch below is left unguarded.)
    if (smd == 0) {
      stop("'smd' must be nonzero to plan a sample size; a null effect has power equal to 'alpha_level' at every N.", call. = FALSE)
    }
  }

  power_at <- function(n_a, n_b) {
    df  <- n_a + n_b - 2
    ncp <- smd * sqrt(n_a * n_b / (n_a + n_b))
    if (directional) {
      crit <- qt(1 - alpha_level, df = df)
      if (smd >= 0) 1 - pt(crit, df = df, ncp = ncp) else pt(-crit, df = df, ncp = ncp)
    } else {
      cu <- qt(1 - alpha_level / 2, df = df)
      cl <- qt(alpha_level / 2,     df = df)
      (1 - pt(cu, df = df, ncp = ncp)) + pt(cl, df = df, ncp = ncp)
    }
  }

  if (!is.null(n_1)) {
    if (is.null(n_2)) n_2 <- n_1
    # Reject fractional or too-small group sizes rather than returning a power
    # table for a design that cannot have them.
    n_1 <- .check_whole_n(n_1, "n_1", 2L)
    n_2 <- .check_whole_n(n_2, "n_2", 2L)
    pwr <- power_at(n_1, n_2)
    ncp <- smd * sqrt(n_1 * n_2 / (n_1 + n_2))
    # The design result is followed by the user-supplied planning inputs, each
    # echoed as its own row so the assumptions the power was evaluated under
    # travel with the result. supposed_smd is the supposed effect the plan is
    # built on, a value the researcher posits (a minimally important effect or
    # a value believed to be true in the population), never a sample estimate.
    # The supplied group sizes are already reported as specified_n_1 /
    # specified_n_2; tails is 2 for a nondirectional test and 1 for a
    # directional test, keeping the value column numeric.
    out <- data.frame(
      term  = c("specified_n_1", "specified_n_2", "actual_power",
                "noncentral_t_parm",
                "supposed_smd", "alpha_level", "tails"),
      value = c(n_1, n_2, pwr, ncp,
                smd, alpha_level, if (directional) 1 else 2)
    )
    class(out) <- c("dmar_ss_power", "dmar_tbl", "data.frame")
    return(out)
  }

  # Start at the smallest admissible per-group size (n = 2) and evaluate before
  # incrementing, so a design whose power already meets the target at n = 2 is
  # returned rather than skipped.
  n <- 2
  pwr <- power_at(n, n)
  while (pwr < desired_power) {
    n <- n + 1
    pwr <- power_at(n, n)
    if (n > 1e6) stop("Failed to converge within reasonable sample size.", call. = FALSE)
  }
  # The design result is followed by the user-supplied planning inputs, each
  # echoed as its own row so the assumptions the design was planned against
  # travel with the result. supposed_smd is the supposed effect the plan is
  # built on, a value the researcher posits (a minimally important effect or a
  # value believed to be true in the population), never a sample estimate;
  # tails is 2 for a nondirectional test and 1 for a directional test, keeping
  # the value column numeric.
  out <- data.frame(
    term  = c("necessary_n_per_group", "actual_power", "noncentral_t_parm",
              "supposed_smd", "desired_power", "alpha_level", "tails"),
    value = c(n, pwr, smd * sqrt(n / 2),
              smd, desired_power, alpha_level, if (directional) 1 else 2)
  )
  class(out) <- c("dmar_ss_power", "dmar_tbl", "data.frame")
  out
}
