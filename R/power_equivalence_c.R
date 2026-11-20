#' Power of the TOST or Noninferiority Test for a Linear Contrast
#'
#' Computes the exact power of the Schuirmann (1987) two one-sided
#' tests procedure, or of the one-sided noninferiority test, for a
#' linear contrast of group means \eqn{\psi = \sum_j c_j \mu_j} with
#' one pooled error term, or for a contrast measured within each unit
#' of a paired design. For equivalence, the power is the
#' probability that the \eqn{(1 - 2\alpha)} confidence interval for
#' \eqn{\psi} lies entirely inside \eqn{(-\delta_L, \delta_U)},
#' computed by numerical integration over the chi distribution of the
#' estimated error standard deviation; for noninferiority, the power
#' is a noncentral \emph{t} probability in closed form. This is the
#' contrast generalization of \code{\link{power_equivalence_md}}.
#'
#' @param c_weights The contrast weights, required for
#'   \code{design = "independent"}. The weights must sum to zero with
#'   the positive weights summing to 1 and the negative weights to -1,
#'   so that the bounds are on the raw scale of the response. Leave
#'   \code{NULL} for \code{design = "paired"}, where the weights are
#'   already in the per-unit contrast scores.
#' @param n_per_group For \code{design = "independent"}, the sample
#'   sizes per group (if length 1, equal group sizes are assumed);
#'   together with \code{c_weights}, \code{n_per_group} determines the
#'   standard error factor \eqn{\sqrt{\sum_j c_j^2 / n_j}}. Leave
#'   \code{NULL} for \code{design = "paired"}, which has no groups of
#'   different units and takes \code{N}; a paired call that gives
#'   \code{n_per_group} is refused.
#' @param N For \code{design = "paired"}, the number of units, each
#'   measured under every condition the contrast involves. Every unit
#'   contributes one contrast score, so this is the total sample size.
#'   Leave \code{NULL} for \code{design = "independent"}, where the
#'   sizes are given per group in \code{n_per_group}; an independent
#'   call that gives \code{N} is refused.
#' @param true_sd For \code{design = "independent"}, the error standard
#'   deviation (the square root of the mean square error). For
#'   \code{design = "paired"}, the standard deviation of the per-unit
#'   contrast scores (for two conditions, of the difference scores). It is
#'   \eqn{\sigma} in the formulas below.
#' @param contrast_lower,contrast_upper Bounds on the raw scale of the
#'   response, each positive; the equivalence region is
#'   \eqn{(-\delta_L, +\delta_U)}, where \eqn{\delta_L} is
#'   \code{contrast_lower} and \eqn{\delta_U} is \code{contrast_upper}. For \code{side = "equivalence"},
#'   \code{contrast_upper} is required and \code{contrast_lower} defaults to
#'   it (symmetric bounds). For \code{side = "noninferiority"}, the test
#'   is against \eqn{-\delta_L} alone, so \code{contrast_lower} is the
#'   argument to give; a call that gives only the upper bound uses it as
#'   the margin.
#' @param true_contrast The population value of the contrast at which the
#'   power is evaluated (\eqn{\psi} in the formulas). Default
#'   \code{0}, the most favorable point for an equivalence declaration.
#' @param alpha_level One-sided significance level for each test. Default
#'   \code{0.05}.
#' @param side \code{"equivalence"} (default) for the TOST power, or
#'   \code{"noninferiority"} for the one-sided test against
#'   \eqn{-\delta_L}.
#' @param df_error The error degrees of freedom, \eqn{\nu} in the
#'   formulas below. Defaults to
#'   \eqn{\sum_j n_j - J} for \code{design = "independent"}, the sum of
#'   \code{n_per_group} less the number of groups, and to \eqn{N - 1}
#'   for \code{design = "paired"}; supply it directly when the error
#'   term comes from a model with more groups or additional predictors
#'   than the contrast involves (for example, a five-group model
#'   supplying the pooled error for a two-group contrast).
#' @param design \code{"independent"} (default) for groups of different
#'   units sharing one pooled error term, or \code{"paired"} for a
#'   contrast measured within each unit (two conditions, or a pre and a
#'   post measurement, on the same units), analyzed as the mean of the
#'   \eqn{N} per-unit contrast scores. The design decides which size
#'   argument is read: \code{n_per_group} for independent groups,
#'   \code{N} for a paired design. See Details.
#'
#' @return A one-row \code{data.frame} with columns \code{term}
#'   (\code{"power"}) and \code{value} (the computed power, in
#'   \eqn{[0, 1]}).
#'
#' @details
#' \strong{Equivalence power.} Conditional on the estimated error
#' standard deviation \eqn{S}, the \eqn{(1 - 2\alpha)} CI fits inside
#' the bounds on a computable event, and the unconditional power
#' integrates that event over the scaled chi distribution of \eqn{S}
#' on \code{df_error} degrees of freedom. With \code{c_weights =
#' c(1, -1)} and equal \code{n_per_group}, the result reproduces
#' \code{\link{power_equivalence_md}} exactly.
#'
#' \strong{Noninferiority power.} The one-sided test rejects when
#' \eqn{t = (\hat\psi + \delta_L)/\mathrm{SE}(\hat\psi)} exceeds
#' \eqn{t_{1-\alpha,\nu}}, so the power is
#' \eqn{\Pr(T'_{\nu}(\lambda) > t_{1-\alpha,\nu})} with noncentrality
#' \eqn{\lambda = (\psi + \delta_L)/(\sigma \sqrt{\sum_j c_j^2/n_j})}.
#'
#' \strong{Paired designs.} With \code{design = "paired"}, each unit
#' contributes one contrast score, \eqn{d_i = \sum_j c_j y_{ij}}, and
#' the analysis is of their mean: the standard error is
#' \eqn{\sigma / \sqrt{N}} on \eqn{N - 1} degrees of freedom, where
#' \eqn{N} is the number of units (the argument \code{N}), so
#' \eqn{\sqrt{\sum_j c_j^2 / n_j}} above becomes \eqn{1 / \sqrt{N}}
#' and both powers are computed as for the independent design. The
#' weights that form the scores should follow the same rule as
#' \code{c_weights} (positive weights summing to 1, negative to -1), so
#' that the bounds are on the raw scale; for two conditions the score
#' is the difference. Its standard deviation \code{true_sd} carries the
#' correlation between the measures: for two measures with a common
#' standard deviation \eqn{\sigma_y} correlated \eqn{\rho}, it is
#' \eqn{\sigma_y \sqrt{2(1 - \rho)}}, smaller than
#' \eqn{\sigma_y \sqrt{2}}, the standard deviation of a difference
#' between independent observations, whenever \eqn{\rho > 0}; that is
#' the precision a paired design buys.
#'
#' \strong{The feasibility condition.} If the expected half-width of
#' the CI is not smaller than the bounds allow, the equivalence power
#' is zero or near zero regardless of \code{true_contrast}: an imprecise
#' design cannot declare equivalence even when the arms are truly
#' identical. Planning should target a half-width of about half the
#' bound; see \code{\link{ss_power_equivalence_c}} and
#' \code{\link{ss_aipe_c}}.
#'
#' @references
#' Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal,
#'   J. J. (2025). A sequential approach for noninferiority or
#'   equivalence of a linear contrast under cost constraints.
#'   \emph{Psychological Methods, 30}(2), 425--439. \doi{10.1037/met0000570}
#'
#' Phillips, K. F. (1990). Power of the two one-sided tests procedure
#'   in bioequivalence. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 18}(2), 137--144. \doi{10.1007/BF01063556}
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' @seealso \code{\link{power_equivalence_md}},
#'   \code{\link{ss_power_equivalence_c}}, \code{\link{equivalence_c}},
#'   \code{\link{ss_aipe_c}}
#'
#' @examples
#' # 1. Two groups of 61 and 113 sharing a five-group pooled error term
#' #    (so df_error = 404 - 5 = 399), bounds of 5 raw-scale points:
#' #    the design's probability of declaring equivalence when the
#' #    groups are truly identical.
#' power_equivalence_c(c_weights = c(1, -1), n_per_group = c(61, 113),
#'                     true_sd = 15.67, contrast_upper = 5,
#'                     true_contrast = 0, df_error = 399)
#'
#' # 2. The same design's noninferiority power at the same point.
#' power_equivalence_c(c_weights = c(1, -1), n_per_group = c(61, 113),
#'                     true_sd = 15.67, contrast_upper = 5,
#'                     true_contrast = 0, df_error = 399,
#'                     side = "noninferiority")
#'
#' # 3. Agreement with power_equivalence_md() in the two-group case
#' #    (Phillips, 1990, Table 1: expected 0.8029678).
#' power_equivalence_c(c_weights = c(1, -1), n_per_group = 24, true_sd = 0.20,
#'                     contrast_lower = 0.2, contrast_upper = 0.2,
#'                     true_contrast = 0.05, df_error = 22)
#'
#' # 4. A paired design: 50 units measured under both conditions, with
#' #    difference scores whose SD is 1, and a noninferiority margin of
#' #    0.5. No upper bound is needed for a noninferiority question.
#' power_equivalence_c(N = 50, true_sd = 1, contrast_lower = 0.5,
#'                     side = "noninferiority", design = "paired")
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @family equivalence testing
#'
#' @export

power_equivalence_c <- function(c_weights = NULL, n_per_group = NULL, N = NULL,
                                true_sd,
                                contrast_lower = NULL, contrast_upper = NULL,
                                true_contrast = 0, alpha_level = 0.05,
                                side = c("equivalence", "noninferiority"),
                                df_error = NULL,
                                design = c("independent", "paired")) {
  side   <- match.arg(side)
  design <- match.arg(design)

  # The equivalence test uses both bounds and the noninferiority test
  # -contrast_lower alone. DMAR 1.0.x required the upper bound on either side
  # and took it for both bounds, so a noninferiority call that gives only
  # contrast_upper still uses it as the margin. A bound that is given is
  # checked even when the test does not use it.
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
       !is.finite(contrast_upper) || contrast_upper <= 0))
    .dmar_refuse_argument(
      "'contrast_upper' must be a single positive number.", arg = "contrast_upper")
  if (is.null(contrast_lower)) contrast_lower <- contrast_upper
  if (!is.numeric(contrast_lower) || length(contrast_lower) != 1L ||
      !is.finite(contrast_lower) || contrast_lower <= 0)
    .dmar_refuse_argument(
      "'contrast_lower' must be a single positive number.", arg = "contrast_lower")
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L || is.na(alpha_level) ||
      alpha_level <= 0 || alpha_level >= 0.5)
    .dmar_refuse_argument(
      "'alpha_level' must be a single number in (0, 0.5).", arg = "alpha_level")
  if (!is.numeric(true_sd) || length(true_sd) != 1L || !is.finite(true_sd) || true_sd <= 0)
    .dmar_refuse_argument(
      "'true_sd' must be a single positive number.", arg = "true_sd")
  if (!is.numeric(true_contrast) || length(true_contrast) != 1L || !is.finite(true_contrast))
    .dmar_refuse_argument(
      "'true_contrast' must be a single numeric value.", arg = "true_contrast")

  # Each design has its own size argument, because the two sizes are
  # different quantities: sizes per group for independent groups, and the
  # number of units (the total sample size) for a paired design. The
  # argument of the other design is refused, so a size is never read as
  # the wrong quantity.
  if (design == "independent") {
    if (!is.null(N))
      .dmar_refuse_argument(c(
        "'N' was given with design = \"independent\". 'N' is the number of ",
        "units of a paired design, in which every unit is measured under ",
        "every condition the contrast involves. Independent groups are ",
        "sized per group: give 'n_per_group' (one number for equal groups, ",
        "or one size for each contrast weight) and leave 'N' NULL. If each ",
        "unit is measured under every condition, set design = \"paired\"."),
        arg = c("N", "design"))
    if (is.null(n_per_group))
      .dmar_refuse_argument(c(
        "'n_per_group' must be given for design = \"independent\": the ",
        "sample size of each group, as one number for equal groups or one ",
        "size for each contrast weight. The standard error of the contrast ",
        "is the error standard deviation times sqrt(sum(c_weights^2 / ",
        "n_per_group)), so the power cannot be computed without it."),
        arg = "n_per_group")
    if (!is.numeric(n_per_group) || anyNA(n_per_group) || any(!is.finite(n_per_group)))
      .dmar_refuse_argument(
        "'n_per_group' must be finite and not missing.", arg = "n_per_group")
    if (is.null(c_weights))
      .dmar_refuse_argument(
        "'c_weights' must be supplied for design = \"independent\".",
        arg = "c_weights")
    if (length(n_per_group) == 1L) n_per_group <- rep(n_per_group, length(c_weights))
    if (length(n_per_group) != length(c_weights))
      .dmar_refuse_argument(
        "The lengths of 'n_per_group' and 'c_weights' differ, which should not be the case.",
        arg = c("n_per_group", "c_weights"))
    if (any(n_per_group < 2)) .dmar_refuse_argument(
      "Each group needs at least 2 observations.", arg = NULL)
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

    k  <- sqrt(sum(c_weights^2 / n_per_group))
    nu <- if (is.null(df_error)) sum(n_per_group) - length(c_weights) else df_error
  } else {
    if (!is.null(n_per_group))
      .dmar_refuse_argument(c(
        "'n_per_group' was given with design = \"paired\". A paired design ",
        "has one set of units, each measured under every condition the ",
        "contrast involves, so there are no groups of different units to ",
        "size. Give the number of units as 'N' and leave 'n_per_group' ",
        "NULL. If the conditions are observed on different units, set ",
        "design = \"independent\"."),
        arg = c("n_per_group", "design"))
    if (!is.null(c_weights))
      .dmar_refuse_argument(c(
        "Leave 'c_weights' NULL for design = \"paired\": 'true_sd' is the ",
        "standard deviation of the per-unit contrast scores, which ",
        "already carry the weights."), arg = c("c_weights", "true_sd"))
    if (is.null(N) || !is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 2)
      .dmar_refuse_argument(c(
        "For design = \"paired\", 'N' must be the number of units, a ",
        "single finite number of at least 2. Each unit contributes one ",
        "contrast score, and the standard deviation of the scores is ",
        "estimated on N - 1 degrees of freedom, which needs at least 2 ",
        "units."),
        arg = "N")

    # The mean of N per-unit contrast scores: the standard error factor is
    # 1 / sqrt(N), on N - 1 degrees of freedom.
    k  <- 1 / sqrt(N)
    nu <- if (is.null(df_error)) N - 1 else df_error
  }
  if (!is.numeric(nu) || length(nu) != 1L || is.na(nu) || !is.finite(nu) || nu <= 0)
    .dmar_refuse_argument(
      "The error degrees of freedom must be a single positive number.",
      arg = NULL)

  power <- .power_equivalence_c_fast(k = k, df_error = nu, true_sd = true_sd,
                                     contrast_lower = contrast_lower,
                                     contrast_upper = contrast_upper,
                                     true_contrast = true_contrast, alpha_level = alpha_level,
                                     side = side)

  .as_dmar_tbl(data.frame(term = "power", value = power,
                          stringsAsFactors = FALSE))
}


# --- Fast path: bare-numeric power for the search loops -------------------
# Skips the argument validation and data.frame construction of the public
# power_equivalence_c(); consumed by ss_power_equivalence_c()'s sample size
# search, which evaluates the power at many candidate sizes per invocation.
# 'k' is the standard error factor: sqrt(sum(c_weights^2 / n_per_group)) for
# independent groups, 1 / sqrt(N) for a paired design. The noninferiority
# branch never reads contrast_upper, which is NULL when the caller planned
# noninferiority without one.
.power_equivalence_c_fast <- function(k, df_error, true_sd, contrast_lower,
                                      contrast_upper, true_contrast,
                                      alpha_level, side) {
  if (side == "noninferiority") {
    # The noncentral t upper tail, exact beyond a noncentrality of 37.62
    # (.nc_t_prob(), R/nc_t_internals.R).
    ncp <- (true_contrast + contrast_lower) / (true_sd * k)
    return(.nc_t_prob(stats::qt(1 - alpha_level, df_error), df_error, ncp,
                      lower_tail = FALSE))
  }

  # Equivalence: integrate the conditional probability that the interval fits
  # inside the bounds over the sampling distribution of the estimated error
  # standard deviation. The bounds enter only through their distance from
  # 'true_contrast' in standard errors, which is the unit-free scale the shared
  # quadrature works on. That quadrature (.power_equivalence_tost_quad(), in
  # power_equivalence_md.R) documents the integrand and explains why its
  # limits are quantiles of the sampling distribution of the error standard
  # deviation rather than the raw interval that deviation must fall in.
  one_over_se <- 1 / (true_sd * k)

  .power_equivalence_tost_quad(
    power_t  = stats::qt(1 - alpha_level, df = df_error),
    z_lower  = (-contrast_lower - true_contrast) * one_over_se,
    z_upper  = ( contrast_upper - true_contrast) * one_over_se,
    df_error = df_error
  )
}
