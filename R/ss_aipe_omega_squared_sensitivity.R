# Sensitivity analysis for AIPE on omega squared.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for Omega Squared
#'
#' @description
#' Quantifies how much misspecification of the population
#' \eqn{\omega^2} distorts an AIPE-based sample size plan. The planner
#' \code{\link{ss_aipe_omega_squared}} solves for the smallest \emph{N}
#' that yields an expected CI width below the target at the planning
#' value. Here we generate \code{G} datasets from a balanced one-way
#' ANOVA with population \eqn{\omega^2 = }\code{true_omega_squared} and
#' \code{df_effect + 1} groups at the planner-recommended \emph{N},
#' compute the noncentral \emph{F} confidence interval on each
#' replication via \code{\link{ci_omega_squared}}, and summarize the
#' realized widths and coverage of \code{true_omega_squared}.
#'
#' @param true_omega_squared Population \eqn{\omega^2} (the
#'   data generating value); in \eqn{[0, 1)}.
#' @param estimated_omega_squared \eqn{\omega^2} used to plan the study;
#'   supply this or \code{specified_N} but not both.
#' @param df_effect Numerator degrees of freedom for the omnibus
#'   \emph{F}, equal to the number of groups minus 1.
#' @param width Desired full width of the confidence interval on
#'   \eqn{\omega^2}.
#' @param specified_N Total sample size to evaluate (incompatible with
#'   \code{estimated_omega_squared}).
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance Optional assurance probability passed to
#'   \code{\link{ss_aipe_omega_squared}} when resolving the planned
#'   sample size.
#' @param G Number of Monte Carlo replications (default 1000).
#' @param print_iter Logical. Print iteration index per replication.
#' @param save Logical. If \code{TRUE} write per-replication results to
#'   \code{filename}.
#' @param filename Path used when \code{save = TRUE}.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   the realized \eqn{\hat\omega^2} and interval width, the proportion
#'   of intervals at or below \code{width}, tail-specific and overall
#'   empirical non-coverage of \code{true_omega_squared}, and the input
#'   echoes, including \code{assurance} (present only when an
#' assurance was supplied).
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect
#'   sizes: Theory, application, and implementation. \emph{Journal of
#'   Statistical Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on effect size measures.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_omega_squared}}, \code{\link{ci_omega_squared}}
#'
#' @examples
#' # Well-specified: planner used omega^2 = 0.10, truth is 0.10.
#' # G is kept small here so the example runs quickly; raise it for a
#' # stable sensitivity estimate.
#' set.seed(113)
#' ss_aipe_omega_squared_sensitivity(
#'   true_omega_squared      = 0.10,
#'   estimated_omega_squared = 0.10,
#'   df_effect = 2, width = 0.10,
#'   G = 25, print_iter = FALSE
#' )
#'
#' @keywords design htest
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family AIPE sample size planning
#'
#' @export
ss_aipe_omega_squared_sensitivity <- function(true_omega_squared = NULL,
                                              estimated_omega_squared = NULL,
                                              df_effect, width,
                                              specified_N = NULL,
                                              conf_level = 0.95,
                                              assurance = NULL,
                                              G = 1000, print_iter = FALSE,
                                              save = FALSE,
                                              filename = "ss_aipe_omega_squared_sensitivity_result.csv") {
  if (is.null(estimated_omega_squared) && is.null(specified_N))
    stop("You must specify either 'estimated_omega_squared' or 'specified_N'.", call. = FALSE)
  if (!is.null(estimated_omega_squared) && !is.null(specified_N))
    stop("You must specify 'estimated_omega_squared' or 'specified_N', but not both.", call. = FALSE)
  if (is.null(true_omega_squared) || !is.numeric(true_omega_squared) ||
      true_omega_squared < 0 || true_omega_squared >= 1)
    stop("'true_omega_squared' must be a single value in [0, 1).", call. = FALSE)
  if (!is.numeric(df_effect) || df_effect < 1)
    stop("'df_effect' must be a positive integer (groups - 1).", call. = FALSE)

  if (!is.null(estimated_omega_squared)) {
    plan <- suppressWarnings(
      ss_aipe_omega_squared(population_omega_squared = estimated_omega_squared,
                            df_effect = df_effect, width = width,
                            conf_level = conf_level, assurance = assurance)
    )
    N <- plan$value[plan$term == "necessary_N"]
  } else {
    N <- specified_N
  }
  N <- as.integer(N)
  J <- as.integer(df_effect + 1)
  n_per_group <- N %/% J
  if (n_per_group < 2)
    stop("Resolved per-group sample size is < 2.", call. = FALSE)
  N <- n_per_group * J

  # Group means chosen to deliver the requested population omega^2.
  # For a balanced design with total variance 1, omega^2 = sigma2_alpha /
  # (sigma2_alpha + sigma2_e), so set sigma2_alpha = true_omega_squared and
  # sigma2_e = 1 - true_omega_squared; spread means symmetrically with
  # variance equal to sigma2_alpha.
  sigma2_e   <- 1 - true_omega_squared
  sigma2_a   <- true_omega_squared
  mu_centers <- seq(-1, 1, length.out = J)
  mu_centers <- mu_centers - mean(mu_centers)
  if (sigma2_a > 0) {
    mu <- mu_centers * sqrt(sigma2_a / stats::var(mu_centers))
  } else {
    mu <- rep(0, J)
  }

  omega_hat <- numeric(G)
  ci_lo     <- numeric(G)
  ci_hi     <- numeric(G)
  ci_width  <- numeric(G)
  tI_lower  <- logical(G)
  tI_upper  <- logical(G)

  sd_e <- sqrt(sigma2_e)
  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    y <- unlist(lapply(seq_len(J), function(j)
      stats::rnorm(n_per_group, mean = mu[j], sd = sd_e)))
    grp <- factor(rep(seq_len(J), each = n_per_group))
    fit <- stats::aov(y ~ grp)
    # ci_omega_squared() returns a wide data.frame with columns named
    # "omega_squared", "lower_limit", "upper_limit", "F_value", ...
    ci <- suppressWarnings(ci_omega_squared(object = fit, conf_level = conf_level))
    om  <- as.numeric(ci$omega_squared[1])
    lo  <- as.numeric(ci$lower_limit[1])
    hi  <- as.numeric(ci$upper_limit[1])
    omega_hat[g] <- om
    ci_lo[g]     <- lo
    ci_hi[g]     <- hi
    ci_width[g]  <- hi - lo
    tI_lower[g]  <- true_omega_squared < lo
    tI_upper[g]  <- true_omega_squared > hi
  }

  if (isTRUE(save)) {
    per_rep <- data.frame(omega_hat = omega_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_omega_squared", "median_omega_squared", "sd_omega_squared",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_w",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "total_N", "n_per_group",
              "true_omega_squared", "estimated_omega_squared",
              "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(omega_hat, na.rm = TRUE),
              stats::median(omega_hat, na.rm = TRUE),
              stats::sd(omega_hat, na.rm = TRUE),
              mean(ci_width, na.rm = TRUE),
              stats::median(ci_width, na.rm = TRUE),
              stats::sd(ci_width, na.rm = TRUE),
              mean(ci_width <= width, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE), mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              N, n_per_group,
              true_omega_squared,
              if (is.null(estimated_omega_squared)) NA_real_ else estimated_omega_squared,
              width, conf_level,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level)
}
