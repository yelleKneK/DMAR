# Sample size for the TOST or noninferiority test on a linear contrast.
#' Sample Size for Equivalence or Noninferiority of a Linear Contrast
#'
#' Computes the smallest sample size (per group, or the number of units
#' in a paired design) at which the two one-sided tests procedure
#' (Schuirmann, 1987) for a linear contrast
#' \eqn{\psi = \sum_j c_j \mu_j}, or the companion one-sided
#' noninferiority test, attains a desired power. Power is computed
#' exactly through \code{\link{power_equivalence_c}}. This is the
#' declaration-probability route to planning; the accuracy in
#' parameter estimation (AIPE) route, which targets the confidence
#' interval width directly, is \code{\link{ss_aipe_c}}, and the two
#' answer the same question whenever \code{true_contrast = 0} and the width
#' target is calibrated to the bounds.
#'
#' @param c_weights The contrast weights, required for
#'   \code{design = "independent"}. The weights must sum to zero with
#'   the positive weights summing to 1 and the negative weights to -1,
#'   so that the bounds are on the raw scale of the response. Leave
#'   \code{NULL} for \code{design = "paired"}, where the weights are
#'   already in the per-unit contrast scores.
#' @param true_sd For \code{design = "independent"}, the anticipated error
#'   standard deviation (the square root of the mean square error). For
#'   \code{design = "paired"}, the anticipated standard deviation of the
#'   per-unit contrast scores (for two conditions, of the difference
#'   scores). It is \eqn{\sigma} in the formulas below.
#' @param contrast_lower,contrast_upper Bounds on the raw scale of the
#'   response, each positive; the equivalence region is
#'   \eqn{(-\delta_L, +\delta_U)}, where \eqn{\delta_L} is
#'   \code{contrast_lower} and \eqn{\delta_U} is \code{contrast_upper}.
#'   For \code{side = "equivalence"}, \code{contrast_upper} is required
#'   and \code{contrast_lower} defaults to it (symmetric bounds). For
#'   \code{side = "noninferiority"}, the test is against \eqn{-\delta_L} alone, so \code{contrast_lower} is the
#'   argument to give; a call that gives only the upper bound uses it as
#'   the margin.
#' @param true_contrast The population contrast (\eqn{\psi}) the design
#'   should be able to detect as equivalent (or noninferior). Default
#'   \code{0}. For equivalence it must lie strictly inside the bounds; for
#'   noninferiority, strictly above \eqn{-\delta_L}. The farther
#'   \code{true_contrast} sits from the center of the region, the larger
#'   the required sample size.
#' @param desired_power The target probability of declaring
#'   equivalence (or noninferiority) at \code{true_contrast}. Default
#'   \code{0.85}, matching \code{\link{ss_power_contrast}}.
#' @param alpha_level One-sided significance level for each test. Default
#'   \code{0.05}.
#' @param side \code{"equivalence"} (default) or
#'   \code{"noninferiority"}.
#' @param design \code{"independent"} (default) for groups of different
#'   units sharing one pooled error term, or \code{"paired"} for a
#'   contrast measured within each unit, analyzed as the mean of the
#'   per-unit contrast scores (see \code{\link{power_equivalence_c}}).
#'
#' @return A \code{data.frame} with columns \code{term} and
#'   \code{value}. For \code{design = "independent"} the rows are
#'   \code{necessary_n_per_group} (the recommended sample size of each
#'   group the contrast involves), \code{total_N} (the number of units
#'   the design needs, \eqn{J \times n} for the \eqn{J} groups named by
#'   \code{c_weights}), and \code{actual_power} (the exact power achieved
#'   at the recommendation). For \code{design = "paired"} the rows are
#'   \code{necessary_N} and \code{actual_power}: each condition is
#'   observed on every unit, so the answer is the number of units, which
#'   is the total sample size, and there is no per-group size to report
#'   beside it. The row is named for the quantity it holds; it is the
#'   value \code{\link{power_equivalence_c}} takes as \code{N}, as
#'   \code{necessary_n_per_group} is the value it takes as
#'   \code{n_per_group}. The result carries the \code{dmar_ss_power}
#'   class, so \code{\link[generics]{tidy}} and
#'   \code{\link[generics]{glance}} summarize it in broom convention;
#'   \code{tidy()} reports the planned size (per group, or the number of
#'   units) and the power in the same columns for both designs.
#'
#' @details
#' \strong{Design.} Planning for independent groups assumes equal
#' allocation across the \eqn{J} groups named by \code{c_weights} and a
#' pooled error term on \eqn{N - J} degrees of freedom. Groups with zero
#' weight still contribute error degrees of freedom, which is why they
#' belong in \code{c_weights} when the fitted model will include them.
#' Planning for a paired design counts units. Each condition is observed
#' on every unit, so the answer is the total sample size \eqn{N}, and
#' the mean of the \eqn{N} per-unit contrast scores has standard error
#' \eqn{\sigma/\sqrt{N}} on \eqn{N - 1} degrees of freedom, where
#' \eqn{\sigma} is the standard deviation of the scores.
#' That standard deviation carries the correlation between the measures
#' (see the paired designs paragraph of
#' \code{\link{power_equivalence_c}}), and a positive correlation is
#' what lets a paired design reach a given power with fewer observations
#' than independent groups need.
#'
#' \strong{The search.} The function starts from the normal-theory
#' approximation and moves to the smallest integer size (per group, or
#' the number of units of a paired design) whose exact power reaches
#' \code{desired_power}. Power is monotone in the size once the design
#' is feasible, so the search is a short walk.
#'
#' \strong{Relation to the half-width rule.} With symmetric bounds
#' \eqn{\pm\delta}, \code{true_contrast = 0}, and \eqn{\alpha = .05},
#' targeting a 90\% CI half-width of \eqn{\delta/2} yields a
#' declaration probability of about .90; this function makes the
#' probability the target directly rather than through the width.
#'
#' @references
#' Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal,
#'   J. J. (2025). A sequential approach for noninferiority or
#'   equivalence of a linear contrast under cost constraints.
#'   \emph{Psychological Methods, 30}(2), 425--439. \doi{10.1037/met0000570}
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size
#'   planning for statistical power and accuracy in parameter
#'   estimation. \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' @seealso \code{\link{power_equivalence_c}}, \code{\link{ss_aipe_c}},
#'   \code{\link{ss_power_contrast}}, \code{\link{equivalence_c}}
#'
#' @examples
#' # 1. Two-group equivalence with bounds of 5 raw-scale points and an
#' #    anticipated error SD of 15.67: n per group for 90% power at
#' #    true equivalence.
#' ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
#'                        contrast_upper = 5, desired_power = 0.90)
#'
#' # 2. Noninferiority is cheaper than equivalence at the same bound.
#' ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
#'                        contrast_upper = 5, desired_power = 0.90,
#'                        side = "noninferiority")
#'
#' # 3. A true contrast off center raises the requirement.
#' ss_power_equivalence_c(c_weights = c(1, -1), true_sd = 15.67,
#'                        contrast_upper = 5, true_contrast = 2,
#'                        desired_power = 0.90)
#'
#' # 4. A paired noninferiority design: difference scores with an
#' #    anticipated SD of 1 and a margin of 0.5. Each unit is measured
#' #    under both conditions, so the answer is the number of units,
#' #    necessary_N.
#' ss_power_equivalence_c(true_sd = 1, contrast_lower = 0.5,
#'                        desired_power = 0.90, side = "noninferiority",
#'                        design = "paired")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @family equivalence testing
#'
#' @family sample size for power
#'
#' @export

ss_power_equivalence_c <- function(c_weights = NULL, true_sd,
                                   contrast_lower = NULL, contrast_upper = NULL,
                                   true_contrast = 0, desired_power = 0.85,
                                   alpha_level = 0.05,
                                   side = c("equivalence", "noninferiority"),
                                   design = c("independent", "paired")) {
  .dmar_check_alpha_level(alpha_level)
  side   <- match.arg(side)
  design <- match.arg(design)

  # The bounds resolve as in power_equivalence_c(): contrast_upper is required
  # for equivalence only, and a noninferiority call that gives contrast_upper
  # alone, as DMAR 1.0.x required, uses it as the margin.
  if (side == "equivalence" && is.null(contrast_upper))
    .dmar_refuse_argument(
      "'contrast_upper' must be specified (the upper equivalence bound).",
      arg = "contrast_upper")
  if (side == "noninferiority" && is.null(contrast_lower) && is.null(contrast_upper))
    .dmar_refuse_argument(c(
      "'contrast_lower' must be specified (the noninferiority margin; the ",
      "test is against -contrast_lower)."), arg = "contrast_lower")
  if (!is.null(contrast_upper) &&
      (!is.numeric(contrast_upper) || length(contrast_upper) != 1L ||
       is.na(contrast_upper) || contrast_upper <= 0))
    .dmar_refuse_argument(
      "'contrast_upper' must be a single positive number.", arg = "contrast_upper")
  if (is.null(contrast_lower)) contrast_lower <- contrast_upper
  if (!is.numeric(contrast_lower) || length(contrast_lower) != 1L ||
      is.na(contrast_lower) || contrast_lower <= 0)
    .dmar_refuse_argument(
      "'contrast_lower' must be a single positive number.", arg = "contrast_lower")
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L || alpha_level <= 0 || alpha_level >= 0.5)
    .dmar_refuse_argument(
      "'alpha_level' must be a single number in (0, 0.5).", arg = "alpha_level")
  if (!is.numeric(true_sd) || length(true_sd) != 1L || true_sd <= 0)
    .dmar_refuse_argument(
      "'true_sd' must be a single positive number.", arg = "true_sd")
  .dmar_check_desired_power(desired_power)
  if (!is.numeric(true_contrast) || length(true_contrast) != 1L)
    .dmar_refuse_argument(
      "'true_contrast' must be a single numeric value.", arg = "true_contrast")
  if (design == "independent") {
    if (is.null(c_weights))
      .dmar_refuse_argument(
        "'c_weights' must be supplied for design = \"independent\".",
        arg = "c_weights")
    if (!identical(round(sum(c_weights), 5), 0))
      .dmar_refuse_argument(
        "The sum of the contrast weights ('c_weights') should equal zero.",
        arg = "c_weights")
    .dmar_check_contrast_nonzero(c_weights)
    pos <- sum(c_weights[c_weights > 0])
    neg <- sum(c_weights[c_weights < 0])
    if (!isTRUE(all.equal(pos, 1)) || !isTRUE(all.equal(neg, -1)))
      .dmar_refuse_argument(c(
        "The positive weights must sum to 1 and the negative weights to ",
        "-1, so that the bounds are on the raw scale of the response."),
        arg = "c_weights")
  } else if (!is.null(c_weights)) {
    .dmar_refuse_argument(c(
      "Leave 'c_weights' NULL for design = \"paired\": 'true_sd' is the ",
      "standard deviation of the per-unit contrast scores, which ",
      "already carry the weights."), arg = c("c_weights", "true_sd"))
  }

  if (side == "equivalence" &&
      (true_contrast <= -contrast_lower || true_contrast >= contrast_upper))
    .dmar_refuse_argument(c(
      "For equivalence, 'true_contrast' must lie strictly inside ",
      "(-contrast_lower, contrast_upper); no sample size can reach the ",
      "desired power at a population contrast outside the bounds."),
      arg = "true_contrast")
  if (side == "noninferiority" && true_contrast <= -contrast_lower)
    .dmar_refuse_argument(c(
      "For noninferiority, 'true_contrast' must exceed -contrast_lower; no ",
      "sample size can reach the desired power at a population ",
      "contrast at or below the bound."), arg = "true_contrast")

  # A paired design is the one-group case of the search below: the mean of
  # n per-unit contrast scores has standard error true_sd * sqrt(1 / n) on
  # n - 1 degrees of freedom, that is, J = 1 and a weight sum of squares of
  # 1. The search and the rows reported are shared.
  if (design == "independent") {
    J   <- length(c_weights)
    ssq <- sum(c_weights^2)
  } else {
    J   <- 1L
    ssq <- 1
  }

  # Normal-theory start: half-width h = z * true_sd * sqrt(ssq / n) must fit
  # between true_contrast and the nearer bound with z_power to spare. The exact
  # search below corrects the approximation.
  z_a <- stats::qnorm(1 - alpha_level)
  z_p <- stats::qnorm(desired_power)
  gap <- if (side == "noninferiority") true_contrast + contrast_lower else
    min(true_contrast + contrast_lower, contrast_upper - true_contrast)
  n <- max(2, ceiling(((z_a + z_p) * true_sd)^2 * ssq / gap^2))

  power_at <- function(n) {
    .power_equivalence_c_fast(k = sqrt(ssq / n), df_error = n * J - J,
                              true_sd = true_sd,
                              contrast_lower = contrast_lower,
                              contrast_upper = contrast_upper,
                              true_contrast = true_contrast, alpha_level = alpha_level,
                              side = side)
  }

  # Walk down while the target still holds, then up until it is reached.
  # Power is monotone in n once feasible, so both walks terminate.
  while (n > 2 && power_at(n - 1) >= desired_power) n <- n - 1
  n_max <- n * 1000
  while (power_at(n) < desired_power) {
    n <- n + 1
    if (n > n_max)
      .dmar_refuse_argument(c(
        "The search did not reach the desired power; check that the ",
        "bounds, 'true_sd', and 'true_contrast' describe a feasible design."),
        arg = c("true_sd", "true_contrast"))
  }

  # Each size row is named for the quantity it holds. Independent groups
  # need n units in each of J groups: necessary_n_per_group, with the total
  # beside it. A paired design observes each condition on every unit, so its
  # answer is the number of units, the total sample size: necessary_N, with
  # no per-group size to report. Both names are in the shared ss_power_*
  # vocabulary (.SS_POWER_SIZE_TERMS in R/dmar_tidiers.R), so tidy() and
  # glance() report the planned size and the power for either design.
  out <- if (design == "independent") {
    data.frame(
      term  = c("necessary_n_per_group", "total_N", "actual_power"),
      value = c(n, n * J, power_at(n)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  } else {
    data.frame(
      term  = c("necessary_N", "actual_power"),
      value = c(n, power_at(n)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  class(out) <- c("dmar_ss_power", class(out))
  .as_dmar_tbl(out)
}
