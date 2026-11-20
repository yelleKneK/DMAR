# AIPE sample size planning for TOST on the SMD.
#' AIPE Sample Size Planning for an Equivalence Test on the Standardized Mean Difference
#'
#' Computes the smallest per-group sample size at which the confidence
#' interval that \code{\link{equivalence_smd}} reports, the
#' 100(1 - 2\eqn{\alpha})\% interval on the standardized mean difference
#' whose position against the equivalence bounds decides the two one-sided
#' tests (Lakens, 2017), is no wider than \code{width} when the sample
#' standardized mean difference equals the planning value \code{true_smd}
#' (accuracy in parameter estimation, AIPE; Kelley & Rausch, 2006). With
#' the conventional \eqn{\alpha = .05} the interval is the 90\% interval.
#' With \code{assurance} the function returns the smallest per-group sample
#' size at which the realized interval is no wider than \code{width} with
#' at least that probability.
#'
#' @param true_smd Planning value of the population standardized mean
#'   difference \eqn{\delta}, a single finite number. Default \code{0}. The interval is narrowest
#'   when the sample standardized mean difference is 0, so for a given
#'   width the plan at 0 is the smallest over planning values.
#' @param width Target full width of the interval on the standardized mean
#'   difference scale (e.g., \code{0.20} for a 90\% interval of width
#'   0.20).
#' @param alpha_level One-sided significance level of each of the two
#'   tests. The interval planned for is at confidence level
#'   \eqn{1 - 2\alpha}, as in \code{\link{equivalence_smd}}. Default
#'   \code{0.05} (90\% interval).
#' @param assurance Optional probability in \eqn{(0, 1)} that the realized
#'   interval is no wider than \code{width}. The function returns the
#'   smallest per-group sample size, at or above the one planned without
#'   it, at which that probability is at least \code{assurance} (Kelley,
#'   Maxwell, & Rausch, 2003; Kelley & Rausch, 2006). Default \code{NULL}
#'   (the plan is for the width at \code{true_smd} only).
#' @param balanced Logical; \code{TRUE} (default) plans equal-\eqn{n}
#'   groups. Unequal-\eqn{n} planning is not yet supported.
#'
#' @return A \code{data.frame} with columns \code{term} and \code{value}
#'   and one row per term:
#'   \describe{
#'     \item{\code{necessary_n_per_group}}{the planned sample size of each
#'       group.}
#'     \item{\code{total_N}}{twice \code{necessary_n_per_group}.}
#'     \item{\code{width}}{the target width, as given.}
#'     \item{\code{true_smd}}{the planning value.}
#'     \item{\code{expected_width}}{the full width of the interval
#'       \code{\link{equivalence_smd}} reports with
#'       \code{necessary_n_per_group} in each group when the sample
#'       standardized mean difference equals \code{true_smd}, that is,
#'       \code{ci_upper} minus \code{ci_lower} of
#'       \code{equivalence_smd(smd = true_smd, n_1 = n, n_2 = n,
#'       alpha_level = alpha_level)} with \code{n} the planned size. It is
#'       no larger than \code{width}. Without \code{assurance}, the same
#'       interval with one fewer case in each group is wider than
#'       \code{width}, unless the plan is the smallest size,
#'       2 per group.}
#'     \item{\code{assurance}}{the requested probability; present only when
#'       \code{assurance} was given.}
#'     \item{\code{assurance_attained}}{the probability, with
#'       \code{necessary_n_per_group} in each group and population
#'       standardized mean difference \code{true_smd}, that the realized
#'       interval is no wider than \code{width}; present only when
#'       \code{assurance} was given. It is computed, not simulated (see
#'       Details), and it is at least \code{assurance}.}
#'   }
#'   The \code{conf_level} attribute records the confidence level of the
#'   interval, \eqn{1 - 2\alpha}.
#'
#' @details
#' \strong{The interval planned for.} \code{\link{equivalence_smd}} reports
#' the interval
#' \deqn{\hat d \pm t_{1 - \alpha,\, n_1 + n_2 - 2}
#'   \sqrt{\frac{n_1 + n_2}{n_1 n_2} + \frac{\hat d^{\,2}}{2 (n_1 + n_2)}},}
#' a large-sample standard error of \eqn{\hat d} with a critical value from
#' Student's \emph{t} distribution on \eqn{n_1 + n_2 - 2} degrees of
#' freedom, and it declares equivalence when that interval lies inside the
#' bounds. This function plans for the same interval, because it is the
#' interval the analysis reports; a plan for any other interval would
#' describe a width the analysis does not show. With \eqn{n} in each group
#' the width is
#' \deqn{w(\hat d, n) \;=\; t_{1 - \alpha,\, 2n - 2}
#'   \sqrt{(8 + \hat d^{\,2}) / n}.}
#' The plan is the smallest \eqn{n} of at least 2 with
#' \eqn{w(\delta, n) \le \omega}, where \eqn{\delta} is \code{true_smd} and
#' \eqn{\omega} is \code{width}. Because \eqn{w} decreases in \eqn{n}, the
#' interval with \eqn{n - 1} in each group is wider than \eqn{\omega}. The
#' search starts at \eqn{z_{1 - \alpha}^2 (8 + \delta^2) / \omega^2}, the
#' answer with the standard normal quantile in place of the \emph{t}
#' quantile, which is not above the answer because the \emph{t} quantile is
#' the larger.
#'
#' \strong{Assurance.} Without \code{assurance} the plan makes the interval
#' no wider than the target when the sample standardized mean difference
#' equals \code{true_smd}; a sample with a larger \eqn{|\hat d|} gives a
#' wider interval. With \code{assurance} the function finds the sample size
#' at which the realized interval is no wider than the target with at least
#' that probability, from the sampling distribution of \eqn{\hat d}, as
#' Kelley and Rausch (2006, pp. 372--373) do for the noncentral \emph{t}
#' interval:
#' \enumerate{
#'   \item At a given \eqn{n} the width depends on the sample only through
#'     \eqn{\hat d} and increases with \eqn{|\hat d|}, so the interval is no
#'     wider than \eqn{\omega} exactly when \eqn{|\hat d| \le c(n)}, where
#'     \eqn{c(n)^2 = n \omega^2 / t_{1 - \alpha,\, 2n - 2}^2 - 8}. When
#'     \eqn{c(n)^2 \le 0}, even the interval at \eqn{\hat d = 0} is wider
#'     than \eqn{\omega} and the probability is 0.
#'   \item For two independent samples from normal populations with a
#'     common variance, \eqn{T = \hat d \sqrt{n / 2}} is the two-group
#'     \emph{t} statistic, distributed as noncentral \emph{t} on
#'     \eqn{\nu = 2n - 2} degrees of freedom with noncentrality
#'     \eqn{\lambda = \delta \sqrt{n / 2}}. The probability is
#'     \eqn{P_n = P(|T| \le c(n) \sqrt{n / 2})}, the form of Kelley and
#'     Rausch (2006, Equation 20). With \eqn{T = (Z + \lambda) / \sqrt{V / \nu}},
#'     where \eqn{Z} is standard normal and \eqn{V} is chi square on
#'     \eqn{\nu} degrees of freedom, independent of \eqn{Z},
#'     \deqn{P(|T| \le k) \;=\; E_V\Big[\Phi\big(k \sqrt{V / \nu} - \lambda\big)
#'       - \Phi\big(-k \sqrt{V / \nu} - \lambda\big)\Big],}
#'     which is computed by numerical integration over \eqn{V}. Only the
#'     normal distribution function and the chi square density enter,
#'     because R's noncentral \emph{t} distribution function is an
#'     approximation at noncentrality values above 37.62, which plans for
#'     narrow intervals reach.
#'   \item The function returns the smallest \eqn{n}, at or above the sample
#'     size planned without \code{assurance}, with \eqn{P_n \ge}
#'     \code{assurance}, found by doubling a step upward and then bisecting.
#'     \code{assurance_attained} is \eqn{P_n} at that \eqn{n}.
#' }
#' The probability is exact for the normal model, up to numerical
#' integration. Nothing is simulated, so the function draws no random
#' numbers and returns the same plan on every call. The probability agreed
#' with R's noncentral \emph{t} distribution function to within
#' \eqn{10^{-10}} wherever that function is accurate, and with simulation
#' from normal scores: in 24 planning conditions of 20,000 simulated
#' studies each, the proportion of intervals no wider than the target at
#' the planned size differed from \eqn{P_n} by at most 2.0 Monte Carlo
#' standard errors.
#'
#' Three properties of the step are worth knowing.
#' \itemize{
#'   \item The answer is never below the sample size planned without
#'     \code{assurance}, and an \code{assurance} at or below the probability
#'     at that size returns it unchanged. At or above that size \eqn{P_n}
#'     did not decrease in any of 48 planning conditions scanned over
#'     every whole \eqn{n} (6,639 sizes in all).
#'   \item At \code{true_smd = 0} the plan without \code{assurance} often
#'     meets a moderate request already. The probability at that size
#'     depends on how far rounding up to a whole number narrows the
#'     interval: over widths from 0.05 to 3 at \eqn{\alpha = .05} it ranged
#'     from .11 to .97, with a median of .85, and for a width of 0.50 it is
#'     .839 with 88 per group.
#'   \item As \eqn{|\delta|} grows, the probability at the plan without
#'     \code{assurance} falls toward one half, because \eqn{|\hat d|}
#'     exceeds \eqn{|\delta|} about half the time, and the step adds more
#'     cases, as it does when the target narrows: for
#'     \code{true_smd = 0.5}, a width of 0.30, and \code{assurance = 0.99}
#'     the plan moves from 249 to 257 per group, and for
#'     \code{true_smd = 1}, a width of 0.20, and \code{assurance = 0.90} from
#'     610 to 621.
#' }
#'
#' \strong{The noncentral t interval.} Kelley and Rausch (2006) plan for the
#' interval from the noncentral \emph{t} distribution, which
#' \code{\link{ci_smd}} computes (Kelley, 2007) and
#' \code{\link{ss_aipe_smd}} plans for; they use the interval with the
#' standard normal quantile only as the starting value of their search
#' (pp. 371--372). The noncentral \emph{t} interval has exact coverage for the
#' normal model, where the interval of \code{\link{equivalence_smd}} is a
#' large-sample approximation. At the same \eqn{\hat d}, \eqn{n}, and level
#' the interval of \code{equivalence_smd()} is usually the wider of the two:
#' for 90\% intervals with \eqn{|\hat d|} up to 3, by 18 to 30 percent with
#' 3 per group and by less than 1 percent with 100 per group. At the 80
#' percent level it can be the narrower, by less than 1 percent, when
#' \eqn{|\hat d|} is 5 or more (\code{\link{equivalence_smd}} says where the
#' comparison was computed). A study
#' that will report \code{ci_smd()} at level \eqn{1 - 2\alpha} is planned
#' with \code{ss_aipe_smd(conf_level = 1 - 2 * alpha_level)}.
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for
#'   standardized effect sizes: Theory, application, and
#'   implementation. \emph{Journal of Statistical Software, 20}(8),
#'   1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
#'   obtaining precision: Delineating methods of sample-size planning.
#'   \emph{Evaluation & the Health Professions, 26}(3), 258--287.
#'   \doi{10.1177/0163278703255242}
#'
#' Lakens, D. (2017). Equivalence tests: A practical primer for
#'   \emph{t} tests, correlations, and meta-analyses. \emph{Social
#'   Psychological and Personality Science, 8}(4), 355--362.
#'   \doi{10.1177/1948550617697177}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons and
#'   Chapter 3 on one-way ANOVA.)
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests
#'   procedure and the power approach for assessing the equivalence of
#'   average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics, 15}(6), 657--680.
#'
#' @seealso \code{\link{equivalence_smd}}, the analysis planned for;
#'   \code{\link{ss_aipe_equivalence_smd_sensitivity}}, which simulates the
#'   plan; \code{\link{ss_aipe_smd}} and \code{\link{ci_smd}}, for the
#'   noncentral \emph{t} interval.
#'
#' @examples
#' # 1. Plan for a 90% interval on the standardized mean difference of
#' #    width 0.20, with a planning value of 0.
#' ss_aipe_equivalence_smd(true_smd = 0, width = 0.20)
#'
#' # 2. The analysis at the planned size, when the sample standardized mean
#' #    difference equals the planning value, reports an interval no wider
#' #    than the target (expected_width is its width).
#' plan <- ss_aipe_equivalence_smd(true_smd = 0.05, width = 0.20)
#' n <- plan$value[plan$term == "necessary_n_per_group"]
#' fit <- equivalence_smd(smd = 0.05, n_1 = n, n_2 = n, smd_upper = 0.3)
#' fit$value[fit$term == "ci_upper"] - fit$value[fit$term == "ci_lower"]
#'
#' # 3. With 80% assurance that the realized interval is no wider than 0.20.
#' #    The last row is the probability at the returned size.
#' ss_aipe_equivalence_smd(true_smd = 0.05, width = 0.20, assurance = 0.80)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @family AIPE sample size planning
#'
#' @export

ss_aipe_equivalence_smd <- function(true_smd = 0,
                             width,
                             alpha_level = 0.05,
                             assurance = NULL,
                             balanced = TRUE) {
  .dmar_check_alpha_level(alpha_level)
  # A width left out is refused with the helper's classed message (NULL
  # stands for it there) rather than R's "argument is missing".
  .dmar_check_width(if (missing(width)) NULL else width,
                    "the standardized mean difference")
  .dmar_check_assurance(assurance)
  .dmar_check_numbers("ss_aipe_equivalence_smd", "true_smd")
  .equivalence_smd_check_smd(true_smd, "true_smd")
  if (alpha_level <= 0 || alpha_level >= 0.5)
    .dmar_refuse_argument(
      "'alpha_level' must be in (0, 0.5).", arg = "alpha_level")
  if (!isTRUE(balanced))
    .dmar_refuse_argument(
      "Only balanced two-group designs are supported.",
      arg = NULL)

  # The plan is for the interval equivalence_smd() reports, formed by the
  # same helper (.equivalence_smd_interval() in R/equivalence_internals.R):
  # the smallest n per group whose interval at a sample SMD equal to
  # true_smd is no wider than the target.
  n <- .equivalence_smd_plan_n(true_smd, width, alpha_level)

  # The assurance step: the smallest n, at or above that plan, at which the
  # probability that the realized interval is no wider than the target is
  # at least `assurance`, computed from the noncentral t distribution of
  # the two-group t statistic (.equivalence_smd_assurance_n()).
  attained <- NULL
  if (!is.null(assurance)) {
    step <- .equivalence_smd_assurance_n(n, true_smd, width, alpha_level,
                                         assurance)
    n <- step$n
    attained <- step$attained
  }

  expected_w <- .equivalence_smd_interval(true_smd, n, n, alpha_level)$width

  out <- data.frame(
    term  = c("necessary_n_per_group", "total_N", "width",
              "true_smd", "expected_width",
              if (!is.null(assurance)) c("assurance", "assurance_attained")),
    value = c(n, 2 * n, width, true_smd, expected_w,
              if (!is.null(assurance)) c(assurance, attained)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = 1 - 2 * alpha_level,
              subclass = "dmar_ss_aipe")
}
