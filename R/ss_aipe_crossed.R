# Sample size planning for the mean of a fully crossed units by raters design.
#' Sample Size Planning From the AIPE Perspective for the Mean of a Crossed Design of Units and Raters
#'
#' @description
#' Plans the numbers of units and raters for a fully crossed design, in
#' which every rater scores every unit once, so that the confidence
#' interval for the mean score is sufficiently narrow at the least cost.
#' The units are the objects of measurement (persons, essays, recordings,
#' or the topics of an information retrieval evaluation) and the raters
#' are whoever or whatever produces the scores (judges, coders, assessors,
#' or automated scoring models). When the raters are a sample from a
#' population of raters, as the units are a sample of units, the interval
#' for the population mean carries the variability of both samples, and
#' no number of units can make up for too few raters. The interval planned
#' for is the one \code{\link{ci_crossed}} computes from the data.
#'
#' @param true_var_unit,true_var_rater,true_var_residual Planning values of the three
#'   variance components: of the unit effects (\eqn{\sigma^2_u}, how much
#'   units differ in their mean score over the population of raters), of
#'   the rater effects (\eqn{\sigma^2_r}, how much raters differ in
#'   leniency), and of the residuals (\eqn{\sigma^2_e}, the unit by rater
#'   interaction together with error, which one score per cell cannot
#'   separate). The first two must be nonnegative and the third positive.
#'   The \code{var_unit}, \code{var_rater}, and \code{var_residual} rows of
#'   \code{\link{ci_crossed}} applied to a pilot array estimate them.
#' @param width Desired full width of the confidence interval for the
#'   mean.
#' @param conf_level Confidence level of the interval (default
#'   \code{0.95}).
#' @param assurance Optional probability, strictly between 0 and 1, that
#'   the realized interval is no wider than \code{width}. The default
#'   \code{NULL} plans for the width at the planning values. A plan with
#'   an assurance is never smaller than the plan without one (see
#'   Details), so an assurance at or below one half returns the plan
#'   without one. An assurance requires \code{quantile = "t"}, since it
#'   concerns the sampling variability of the estimated variance that the
#'   \emph{t} quantile accounts for.
#' @param n_raters Optional fixed number of raters (a whole number of at
#'   least 2), for a study whose panel of raters is already set; only the
#'   number of units is then planned. The default \code{NULL} plans both.
#' @param unit_cost,rater_cost,judgment_cost Costs of adding a unit, of
#'   recruiting a rater, and of one judgment (one rater's score of one
#'   unit), in any common currency. A design with \eqn{n_u} units and
#'   \eqn{n_r} raters costs
#'   \eqn{c_u n_u + c_r n_r + c_j n_u n_r}. The defaults \code{0},
#'   \code{0}, and \code{1} make the cost the number of judgments, so the
#'   plan is the design that collects the fewest scores. The costs must be
#'   nonnegative and not all zero.
#' @param quantile \code{"t"} (the default) or \code{"normal"}: the
#'   quantile the planned interval is formed with (see Details).
#' @param max_units,max_raters The largest numbers of units and of raters
#'   the search considers. Defaults \code{100000} and \code{500}.
#'
#' @return A \code{data.frame} (class \code{dmar_ss_aipe}) with rows
#'   \code{necessary_n_units} and \code{necessary_n_raters} (the design;
#'   the second is \code{specified_n_raters} when \code{n_raters} was
#'   supplied), \code{total_judgments} (\eqn{n_u n_r}, the number of
#'   scores to collect), \code{cost}, \code{expected_width} (the full
#'   width of the interval at the planning values, at most \code{width};
#'   a plan with an assurance guards against an unlucky sample, so its
#'   expected width is no wider, and narrower once the assurance asks for
#'   more than the plan without one delivers), \code{df} (the
#'   Satterthwaite degrees of freedom the interval will have at the
#'   planned design, computed from the planning values),
#'   \code{rater_floor} (see Details), the echoes
#'   \code{true_var_unit}, \code{true_var_rater}, \code{true_var_residual}, \code{width},
#'   and \code{conf_level}, and \code{assurance} when one was supplied.
#'   The quantile the plan used is the \code{"quantile"} attribute.
#'   \code{tidy()} reports the number of units as the sample size, and
#'   \code{glance()} keeps every row.
#'
#' @details
#' \strong{The model and the variance of the mean.} The score of unit
#' \eqn{u} by rater \eqn{r} is
#' \eqn{y_{ur} = \mu + a_u + b_r + e_{ur}}, with independent unit
#' effects, rater effects, and residuals of mean zero and variances
#' \eqn{\sigma^2_u}, \eqn{\sigma^2_r}, and \eqn{\sigma^2_e}. Both units
#' and raters are random, so the target is \eqn{\mu}, the mean over the
#' population of units and the population of raters, and the grand mean
#' of a complete array of \eqn{n_u} units by \eqn{n_r} raters has
#' variance
#' \deqn{V = \frac{\sigma^2_u}{n_u} + \frac{\sigma^2_r}{n_r} +
#'   \frac{\sigma^2_e}{n_u n_r}.}
#' In the vocabulary of generalizability theory (Brennan, 2001) the units
#' are the objects of measurement and the raters a random facet crossed
#' with them.
#'
#' \strong{The rater floor.} Whatever the number of units, \eqn{V} is at
#' least \eqn{\sigma^2_r / n_r}, so with \eqn{n_r} raters the half-width
#' cannot fall below \eqn{z\,\sigma_r / \sqrt{n_r}}, with \eqn{z} the
#' normal quantile. A design meeting a target half-width \eqn{h} (half of
#' \code{width}) under the normal quantile therefore needs more than
#' \eqn{z^2 \sigma^2_r / h^2} raters; that number is the
#' \code{rater_floor} row. The \emph{t} quantile raises the requirement
#' further. When raters vary, adding units cannot make up for a panel
#' that is too small, and a plan with \code{n_raters} at or below the
#' floor is refused with that explanation.
#'
#' \strong{The interval planned for.} \code{\link{ci_crossed}} estimates
#' \eqn{V} by the ANOVA estimates of the three components, a negative
#' estimate truncated at zero, and refers the grand mean to a \emph{t}
#' distribution whose degrees of freedom are Satterthwaite's (1946) for
#' that combination of mean squares. With \code{quantile = "t"}, the
#' default, the plan forms the interval the same way at the planning
#' values: the planned half-width is \eqn{t_{\nu} \sqrt{V}}, where the
#' degrees of freedom \eqn{\nu} are those of the plug-in estimate
#' \eqn{(MS_u + MS_r - MS_e) / (n_u n_r)} computed from the expected
#' mean squares. With \code{quantile = "normal"} the planned half-width
#' is \eqn{z \sqrt{V}}, the half-width of an interval whose variance is
#' known; \code{quantile = "normal"} reproduces planning that uses the
#' normal quantile. Its designs are smaller, and because the analysis
#' interval uses the \emph{t} quantile, its realized intervals exceed the
#' target width more often than those of a plan made with the \emph{t}
#' quantile. The gap is largest when the rater component is at or near
#' zero. The normal plan then uses two raters, since raters add nothing
#' to \eqn{V}, but the analysis still estimates the rater component from
#' those two. With no rater variance, \eqn{MS_r} exceeds \eqn{MS_e} by
#' chance in about a third of samples and the estimate is kept; the
#' Satterthwaite degrees of freedom then fall, to one or two when the
#' unit component is small beside the residual one, and the interval
#' can be several times wider than the target. The \emph{t} plan sees the
#' few degrees of freedom two raters give at the planning values and adds
#' raters or units to make up for them.
#' \code{\link{ss_aipe_crossed_sensitivity}} shows how a given plan
#' fares.
#'
#' \strong{Assurance.} The plan without an assurance gives an interval
#' no wider than the target at the planning values, and the realized
#' interval is narrower or wider as the estimated variance falls below
#' or above \eqn{V}. An assurance \eqn{\gamma} (Kelley & Maxwell, 2003)
#' asks instead for a realized width at most \code{width} with
#' probability at least \eqn{\gamma}. The plan approximates the
#' distribution of the plug-in variance by a scaled chi square on the
#' Satterthwaite degrees of freedom, as Satterthwaite (1946) does, and
#' requires
#' \deqn{t_{\nu} \sqrt{V \max(1, q_{\gamma} / \nu)} \le h,}
#' with \eqn{q_{\gamma}} the \eqn{\gamma} quantile of the chi square
#' distribution on \eqn{\nu} degrees of freedom. Under the same
#' approximation the plan without an assurance already meets the target
#' with probability at least \eqn{P(\chi^2_{\nu} \le \nu)}, which
#' exceeds one half because the median of a chi square lies below its
#' mean: about .58 at \eqn{\nu = 5}, .53 at \eqn{\nu = 30}, and nearer one
#' half as \eqn{\nu} grows. An assurance at or below that probability
#' asks for nothing more, and \eqn{q_{\gamma} / \nu} is then at most 1;
#' taking the larger of 1 and \eqn{q_{\gamma} / \nu} keeps such an
#' assurance from shrinking the design. The rule is that a plan with an
#' assurance is never smaller than the plan without one: it costs at
#' least as much, and with a fixed panel it has at least as many units.
#' An assurance at or below one half always returns the plan without
#' one.
#'
#' \strong{The search.} For every number of raters from 2 to
#' \code{max_raters} (or for \code{n_raters} alone) the planner finds the
#' smallest number of units that meets the target, prices each such
#' design, and returns the cheapest; designs of equal cost are resolved
#' toward fewer raters. With the \emph{t} quantile the planned width is
#' not monotone in the number of units: as units are added, \eqn{V}
#' approaches its rater floor while the degrees of freedom eventually
#' fall toward \eqn{n_r - 1}, and past some point the larger quantile
#' outweighs the smaller variance, so more units lengthen the planned
#' interval. The search therefore looks for the smallest feasible number
#' of units on the stretch where the planned width is still falling, and
#' it counts a number of raters as infeasible when the shortest interval
#' that number allows is too wide. A warning says so when the cheapest
#' design found uses \code{max_raters} raters, since a larger limit might
#' find a cheaper one. When no design within the limits meets the target,
#' the error names the limit that binds: if more units would do, it
#' gives the design they would allow and points to \code{max_units};
#' otherwise more raters are needed, and it says whether the reason is
#' the rater floor or the loss of degrees of freedom under the \emph{t}
#' quantile.
#'
#' \strong{A closed form for the allocation.} The search needs no closed
#' form, but one special case has one, which serves as a check. When only
#' judgments cost (\code{unit_cost = 0}, \code{rater_cost = 0}), the
#' problem of minimizing \eqn{V} for a fixed number of judgments has the
#' form Marcoulides and Goldstein (1990) solved for a two-facet
#' generalizability study, whose error variance is
#' \eqn{\sigma^2_1 / n_1 + \sigma^2_2 / n_2 + \sigma^2_{12} / (n_1 n_2)}
#' under a budget on \eqn{n_1 n_2}: the optimum has
#' \eqn{n_1 / n_2 = \sigma^2_1 / \sigma^2_2} (their Equations 5 and 6).
#' Here that is \eqn{n_u / n_r = \sigma^2_u / \sigma^2_r}. The cost is
#' flat near the optimum, so the least costly design in whole numbers
#' can stray from that ratio while collecting barely more judgments than
#' the continuous optimum.
#'
#' \strong{An information retrieval evaluation.} In a comparison of two
#' retrieval systems over a set of topics, each judged by a set of
#' relevance assessors (people, or automated judges such as language
#' models), the units are the topics, the raters are the assessors, and
#' each cell holds the difference between the two systems' effectiveness
#' scores on the topic as computed from that assessor's judgments. The
#' grand mean is then the mean difference between the systems over the
#' populations of topics and assessors, and \code{true_var_unit} and
#' \code{true_var_rater} describe how much the difference between the systems
#' varies over topics and over assessors (the topic by system and
#' assessor by system interactions).
#'
#' \strong{When the planning values are uncertain.} The plan is only as
#' good as the planning values, and a rater component that is
#' underestimated is costly, since it governs the floor.
#' \code{\link{ss_aipe_crossed_sensitivity}} shows by Monte Carlo
#' simulation what a plan delivers when the true components differ from
#' the planning values. When no trustworthy planning values exist, a
#' two-stage sequential procedure, in which a pilot array determines the
#' final numbers of units and raters, can take their place.
#'
#' @references
#' Brennan, R. L. (2001). \emph{Generalizability theory}. Springer.
#'   \doi{10.1007/978-1-4757-3456-0}
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Marcoulides, G. A., & Goldstein, Z. (1990). The optimization of
#'   generalizability studies with resource constraints.
#'   \emph{Educational and Psychological Measurement, 50}(4), 761--768.
#'   \doi{10.1177/0013164490504004}
#'
#' Satterthwaite, F. E. (1946). An approximate distribution of
#'   estimates of variance components. \emph{Biometrics Bulletin, 2}(6),
#'   110--114. \doi{10.2307/3002019}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_crossed}} for the interval the plan is for,
#'   \code{\link{ss_aipe_crossed_sensitivity}} for what a plan delivers
#'   when the planning values are wrong, and \code{\link{icc}} for the
#'   reliability of a single rater's scores in the same design.
#'
#' @examples
#' # Essays scored by raters. Planning values: essays differ with
#' # variance 0.50, raters differ in leniency with variance 0.10, and the
#' # residual variance is 0.40. The goal is a 95% interval no wider than
#' # 0.40 for the mean score; by default the plan collects the fewest
#' # scores.
#' ss_aipe_crossed(true_var_unit = 0.50, true_var_rater = 0.10, true_var_residual = 0.40,
#'                 width = 0.40)
#'
#' # Recruiting and training a rater costs as much as 40 scores, so the
#' # cheapest design uses fewer raters and more essays.
#' ss_aipe_crossed(true_var_unit = 0.50, true_var_rater = 0.10, true_var_residual = 0.40,
#'                 width = 0.40, rater_cost = 40, judgment_cost = 1)
#'
#' # A panel of 30 raters is already set: plan the essays only, with 80%
#' # assurance that the interval is no wider than 0.40.
#' ss_aipe_crossed(true_var_unit = 0.50, true_var_rater = 0.10, true_var_residual = 0.40,
#'                 width = 0.40, n_raters = 30, assurance = 0.80)
#'
#' @keywords design
#'
#' @family AIPE sample size planning
#'
#' @export
ss_aipe_crossed <- function(true_var_unit, true_var_rater, true_var_residual, width,
                            conf_level = 0.95, assurance = NULL,
                            n_raters = NULL, unit_cost = 0, rater_cost = 0,
                            judgment_cost = 1, quantile = c("t", "normal"),
                            max_units = 100000, max_raters = 500) {
  quantile <- match.arg(quantile)
  .crossed_check_number(true_var_unit, "true_var_unit", lower = 0)
  .crossed_check_number(true_var_rater, "true_var_rater", lower = 0)
  .crossed_check_number(true_var_residual, "true_var_residual", lower = 0,
                        open = TRUE)
  .crossed_check_number(width, "width", lower = 0, open = TRUE)
  .crossed_check_number(conf_level, "conf_level", lower = 0, upper = 1,
                        open = TRUE)
  if (!is.null(assurance)) {
    .crossed_check_number(assurance, "assurance", lower = 0, upper = 1,
                          open = TRUE)
    if (quantile == "normal")
      .dmar_refuse_argument(c(
        "An 'assurance' concerns the sampling variability of the ",
        "estimated variance, which a plan with the normal quantile ",
        "treats as known; use quantile = \"t\" with an assurance."),
        arg = "assurance")
  }
  if (!is.null(n_raters)) .crossed_check_count(n_raters, "n_raters")
  .crossed_check_number(unit_cost, "unit_cost", lower = 0)
  .crossed_check_number(rater_cost, "rater_cost", lower = 0)
  .crossed_check_number(judgment_cost, "judgment_cost", lower = 0)
  if (unit_cost + rater_cost + judgment_cost == 0)
    .dmar_refuse_argument(c(
      "At least one of 'unit_cost', 'rater_cost', and 'judgment_cost' ",
      "must be positive, or every design costs nothing."),
      arg = c("unit_cost", "rater_cost", "judgment_cost"))
  .crossed_check_count(max_units, "max_units")
  .crossed_check_count(max_raters, "max_raters")

  h <- width / 2
  rater_floor <- stats::qnorm(1 - (1 - conf_level) / 2)^2 * true_var_rater / h^2

  plan <- .crossed_plan(true_var_unit, true_var_rater, true_var_residual, width,
                        conf_level, assurance, n_raters, unit_cost,
                        rater_cost, judgment_cost, quantile, max_units,
                        max_raters)
  if (is.na(plan$n_units))
    .dmar_refuse_argument(
      .crossed_infeasible_message(true_var_unit, true_var_rater, true_var_residual,
                                  width, conf_level, assurance,
                                  n_raters, unit_cost, rater_cost,
                                  judgment_cost, quantile, max_units,
                                  max_raters, rater_floor), arg = NULL)
  if (plan$at_limit)
    .dmar_warn(
      sprintf("The least costly design found uses %d raters, the search limit 'max_raters'; a larger 'max_raters' may find a cheaper design.",
              as.integer(max_raters)), "result_qualified")

  # The expected width leaves out the assurance factor: it is the width at
  # the planning values, which an assurance plan makes narrower than
  # 'width' (or leaves at the plan without an assurance; see Details).
  q <- if (quantile == "normal") stats::qnorm(1 - (1 - conf_level) / 2) else
    stats::qt(1 - (1 - conf_level) / 2, plan$df)
  V <- .crossed_design(plan$n_units, plan$n_raters, true_var_unit, true_var_rater,
                       true_var_residual)$variance

  out <- data.frame(
    term  = c("necessary_n_units",
              if (is.null(n_raters)) "necessary_n_raters" else "specified_n_raters",
              "total_judgments", "cost", "expected_width", "df",
              "rater_floor", "true_var_unit", "true_var_rater", "true_var_residual",
              "width", "conf_level",
              if (!is.null(assurance)) "assurance"),
    value = c(plan$n_units, plan$n_raters, plan$n_units * plan$n_raters,
              plan$cost, 2 * q * sqrt(V), plan$df, rater_floor,
              true_var_unit, true_var_rater, true_var_residual, width, conf_level,
              assurance),
    stringsAsFactors = FALSE
  )
  attr(out, "quantile") <- quantile
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
