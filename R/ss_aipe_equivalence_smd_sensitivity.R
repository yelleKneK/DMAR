# Sensitivity analysis for AIPE on the equivalence-test (TOST) SMD.
#' Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Equivalence-Test SMD
#'
#' @description
#' Quantifies how much misspecification of the population standardized
#' mean difference distorts an AIPE-based sample size plan for the
#' two-one-sided-tests (TOST) confidence interval on the standardized mean
#' difference. On each replication the function simulates two normal groups
#' of size \emph{n} per group with population standardized mean difference
#' \code{true_smd}, computes the sample standardized mean difference and the
#' interval \code{\link{equivalence_smd}} reports for it, and summarizes the
#' realized widths, the proportion of replications in which the interval
#' falls entirely inside the equivalence bounds (\code{equivalent}), and the
#' proportions that miss \code{true_smd}.
#'
#' @details
#' \strong{The interval simulated.} The interval is the one
#' \code{\link{equivalence_smd}} reports and
#' \code{\link{ss_aipe_equivalence_smd}} plans for,
#' \deqn{\hat d \pm t_{1 - \alpha,\, 2n - 2}
#'   \sqrt{(8 + \hat d^{\,2}) / (4n)},}
#' a large-sample standard error of \eqn{\hat d} with a Student \emph{t}
#' critical value, at confidence level \code{conf_level}, so that
#' \eqn{\alpha} is (1 - \code{conf_level}) / 2: the default 0.95 is the
#' interval of the two one-sided tests at \eqn{\alpha = .025}, and
#' \code{conf_level = 0.90} is the interval of the tests at .05. A sample
#' size from \code{supposed_smd} is planned by
#' \code{ss_aipe_equivalence_smd(alpha_level = (1 - conf_level) / 2)}, so
#' the plan, the analysis, and the simulation describe one interval at one
#' level: a simulated sample with standardized mean difference \eqn{\hat d}
#' gives the limits \code{equivalence_smd(smd = d, n_1 = n, n_2 = n,
#' alpha_level = (1 - conf_level) / 2)} reports. When \code{true_smd}
#' equals \code{supposed_smd}, \code{pct_ci_less_width} estimates the
#' probability that the interval at the planned size is no wider than
#' \code{width}, which the planner computes and reports as
#' \code{assurance_attained} when \code{assurance} is given; the
#' \code{simulation} attribute (below) gives the Monte Carlo limits of the
#' estimate.
#'
#' The interval is a large-sample approximation, so its coverage of
#' \code{true_smd} (the \code{pct_ci_miss_low}, \code{pct_ci_miss_high},
#' and \code{total_type_I_error} rows) can differ from
#' 1 - \code{conf_level}, most in small samples; the noncentral
#' \emph{t} interval of \code{\link{ci_smd}}, which
#' \code{\link{ss_aipe_smd_sensitivity}} simulates, has exact coverage for
#' the normal model.
#'
#' @param true_smd Population standardized mean difference (the
#'   data generating value), a single finite number. Defaults to
#'   \code{0} (perfect equivalence).
#' @param supposed_smd The supposed population SMD, the planning value
#'   passed to \code{\link{ss_aipe_equivalence_smd}} (not an estimate
#'   computed from data), a single finite number;
#'   supply this or \code{specified_n_per_group} but not both.
#' @param width Desired full width of the two-sided CI on the SMD.
#' @param smd_lower,smd_upper Equivalence bounds on the SMD, as
#'   positive magnitudes with the same meaning as in
#'   \code{\link{equivalence_smd}}: the region is
#'   \eqn{(-\delta_L, +\delta_U)}, with \code{smd_lower} as
#'   \eqn{\delta_L} and \code{smd_upper} as \eqn{\delta_U}.
#'   \code{smd_upper} is required; \code{smd_lower} defaults to
#'   \code{smd_upper} (a symmetric region). The simulator records
#'   whether the realized CI falls entirely inside the region.
#' @param specified_n_per_group Per-group sample size to evaluate, a whole
#'   number of at least 2.
#' @param conf_level Confidence level of the interval (default
#'   \code{0.95}); the two one-sided tests it decides are at level
#'   (1 - \code{conf_level}) / 2.
#' @param assurance Optional assurance probability, passed to
#'   \code{\link{ss_aipe_equivalence_smd}} when the sample size is planned
#'   from \code{supposed_smd}.
#' @param G Number of Monte Carlo replications.
#' @param print_iter Logical.
#' @param filename Optional path for a comma separated file recording
#'   every replication (the sample standardized mean difference, the two
#'   confidence limits, the interval width, whether the interval fell
#'   inside the equivalence region, and two indicators of whether the
#'   interval missed \code{true_smd} below or above): nothing is written when
#'   \code{filename} is \code{NULL} (the default), a new file with a
#'   header row is created otherwise, an existing file at that path is
#'   appended to, and a throwaway run should point it at
#'   \code{tempfile(fileext = ".csv")}.
#' @param seed \code{NULL} (the default) to draw from the current random
#'   number stream, or a single number that makes the call reproducible; a
#'   supplied seed is set for this call only, and the caller's random
#'   number state is restored when the function returns.
#'
#' @return A \code{data.frame} with rows for mean / median / SD of
#'   the realized SMD and CI width, the proportion of intervals at or
#'   below \code{width}, tail-specific and overall non-coverage of
#'   \code{true_smd}, the proportion of intervals classified as
#'   \code{equivalent} (CI fully inside the bounds), and the input
#'   echoes, including \code{assurance} (present only when an
#' assurance was supplied).
#'
#' @references
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11},
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_equivalence_smd}}, \code{\link{equivalence_smd}}, \code{\link{ss_aipe_smd_sensitivity}}
#'
#' @examples
#' # Reduced Monte Carlo sweep (small G) for a fast, illustrative run.
#' set.seed(113)
#' ss_aipe_equivalence_smd_sensitivity(
#'   true_smd      = 0.0,
#'   supposed_smd  = 0.0,
#'   width         = 0.30,
#'   smd_upper     = 0.20,
#'   G = 50, print_iter = FALSE
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
#' @inheritSection ss_aipe_cv_sensitivity The record of the simulation
#' @export
ss_aipe_equivalence_smd_sensitivity <- function(true_smd = 0,
                                         supposed_smd = NULL,
                                         width,
                                         smd_lower = NULL,
                                         smd_upper = NULL,
                                         specified_n_per_group = NULL,
                                         conf_level = 0.95,
                                         assurance = NULL,
                                         G = 1000, print_iter = FALSE,
                                         filename = NULL,
                                         seed = NULL) {
  .dmar_check_conf_level(conf_level)
  # A width left out is refused here, with the helper's classed message
  # (NULL stands for it there), rather than with R's "argument is missing"
  # after the simulation has run.
  .dmar_check_width(if (missing(width)) NULL else width,
                    "the standardized mean difference")
  .dmar_check_assurance(assurance)
  .dmar_check_numbers("ss_aipe_equivalence_smd_sensitivity",
                      c("true_smd", "supposed_smd", "specified_n_per_group",
                        "smd_lower", "smd_upper"))
  if (is.null(supposed_smd) && is.null(specified_n_per_group))
    .dmar_refuse_argument(
      "You must specify either 'supposed_smd' or 'specified_n_per_group'.",
      arg = c("supposed_smd", "specified_n_per_group"))
  if (!is.null(supposed_smd) && !is.null(specified_n_per_group))
    .dmar_refuse_argument(
      "You must specify 'supposed_smd' or 'specified_n_per_group', but not both.",
      arg = c("supposed_smd", "specified_n_per_group"))
  # true_smd generates the data and supposed_smd plans the size, so each
  # must be a single finite number: unchecked, an infinite or missing
  # true_smd makes every simulated interval NaN, and an infinite
  # supposed_smd reaches the planner's search.
  .equivalence_smd_check_smd(true_smd, "true_smd")
  if (!is.null(supposed_smd))
    .equivalence_smd_check_smd(supposed_smd, "supposed_smd")
  # A given size is used as it is, so it must be a whole number of at least
  # 2 (the smallest equivalence_smd() accepts); a fractional size had been
  # truncated without a word.
  if (!is.null(specified_n_per_group) &&
      !(is.numeric(specified_n_per_group) &&
        length(specified_n_per_group) == 1L &&
        is.finite(specified_n_per_group) &&
        specified_n_per_group >= 2 &&
        specified_n_per_group == round(specified_n_per_group)))
    .dmar_refuse_argument(c(
      "'specified_n_per_group' must be a single whole number of at least 2, ",
      "the size of each group; ", .dmar_given(specified_n_per_group), "."),
      arg = "specified_n_per_group", value = specified_n_per_group)

  if (is.null(smd_upper))
    .dmar_refuse_argument(
      "'smd_upper' must be specified (the upper equivalence bound).",
      arg = "smd_upper")
  if (!is.numeric(smd_upper) || length(smd_upper) != 1L ||
      smd_upper <= 0)
    .dmar_refuse_argument(
      "'smd_upper' must be a positive number.", arg = "smd_upper")
  if (is.null(smd_lower)) smd_lower <- smd_upper
  if (!is.numeric(smd_lower) || length(smd_lower) != 1L ||
      smd_lower <= 0)
    .dmar_refuse_argument(
      "'smd_lower' must be a positive number.", arg = "smd_lower")
  .check_filename(filename)
  G <- .dmar_check_count(G, "G")
  .dmar_local_seed(seed)

  # ss_aipe_equivalence_smd() plans the width of a (1 - 2 * alpha) confidence
  # interval, and the simulation below forms the interval at conf_level, so the
  # two match only when alpha = (1 - conf_level) / 2. Using 1 - conf_level
  # planned a 90% width while the 95% interval was scored, underplanning the
  # sample size.
  alpha <- (1 - conf_level) / 2
  if (!is.null(supposed_smd)) {
    plan <- ss_aipe_equivalence_smd(true_smd = supposed_smd, width = width,
                             alpha_level = alpha, assurance = assurance,
                             balanced = TRUE)
    n_per <- plan$value[plan$term == "necessary_n_per_group"]
  } else {
    n_per <- specified_n_per_group
  }
  # Both sources are whole numbers of at least 2: the planner's search
  # starts at 2, and a given size is checked above. The size stays a double,
  # so products of sizes cannot overflow R's integer range.
  n_per <- as.numeric(n_per)

  smd_hat    <- numeric(G)
  ci_lo      <- numeric(G)
  ci_hi      <- numeric(G)
  ci_width   <- numeric(G)
  equivalent <- logical(G)
  tI_lower   <- logical(G)
  tI_upper   <- logical(G)

  for (g in seq_len(G)) {
    if (isTRUE(print_iter)) cat(g, "\n")
    x1 <- stats::rnorm(n_per, mean = true_smd, sd = 1)
    x2 <- stats::rnorm(n_per, mean = 0,        sd = 1)
    m_diff <- mean(x1) - mean(x2)
    s_pool <- sqrt(((n_per - 1) * stats::var(x1) +
                    (n_per - 1) * stats::var(x2)) / (2 * n_per - 2))
    d <- m_diff / s_pool

    # The interval equivalence_smd() reports for this sample, at
    # conf_level = 1 - 2 * alpha, from the helper the analysis and the
    # planner also use (R/equivalence_internals.R).
    ci <- .equivalence_smd_interval(d, n_per, n_per, alpha)
    lo <- ci$lower
    hi <- ci$upper
    smd_hat[g]    <- d
    ci_lo[g]      <- lo
    ci_hi[g]      <- hi
    ci_width[g]   <- ci$width
    equivalent[g] <- (lo > -smd_lower) && (hi < smd_upper)
    tI_lower[g]   <- true_smd < lo
    tI_upper[g]   <- true_smd > hi
  }

  if (!is.null(filename)) {
    per_rep <- data.frame(smd = smd_hat, ci_lower = ci_lo,
                          ci_upper = ci_hi, ci_width = ci_width,
                          equivalent = equivalent,
                          type_I_lower = tI_lower, type_I_upper = tI_upper)
    .write_sensitivity_csv(per_rep, filename)
  }

  out <- data.frame(
    term  = c("mean_smd", "median_smd", "sd_smd",
              "mean_ci_width", "median_ci_width", "sd_ci_width",
              "pct_ci_less_width", "pct_equivalent",
              "pct_ci_miss_low", "pct_ci_miss_high", "total_type_I_error",
              "n_per_group", "total_N",
              "true_smd", "supposed_smd",
              "width", "conf_level",
              "smd_lower", "smd_upper",
              if (!is.null(assurance)) "assurance"),
    value = c(mean(smd_hat, na.rm = TRUE),
              stats::median(smd_hat, na.rm = TRUE),
              stats::sd(smd_hat, na.rm = TRUE),
              mean(ci_width, na.rm = TRUE),
              stats::median(ci_width, na.rm = TRUE),
              stats::sd(ci_width, na.rm = TRUE),
              mean(ci_width <= width, na.rm = TRUE),
              mean(equivalent, na.rm = TRUE),
              mean(tI_lower, na.rm = TRUE),
              mean(tI_upper, na.rm = TRUE),
              mean(tI_lower | tI_upper, na.rm = TRUE),
              n_per, 2 * n_per,
              true_smd,
              if (is.null(supposed_smd)) NA_real_ else supposed_smd,
              width, conf_level,
              -smd_lower, smd_upper,
              if (!is.null(assurance)) assurance),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .ss_sens_record(.as_dmar_tbl(out, conf_level = conf_level), G,
                 !is.na(ci_width), conf_level)
}
