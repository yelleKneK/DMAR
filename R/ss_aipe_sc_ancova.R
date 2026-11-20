#' Sample Size Planning From the AIPE Perspective for Standardized ANCOVA Contrasts
#'
#' Sample size planning from the accuracy in parameter estimation (AIPE) perspective for standardized ANCOVA contrasts.
#'
#' @param true_contrast The population unstandardized ANCOVA (adjusted) contrast
#' @param true_sd_anova The population error standard deviation of the ANOVA model
#' @param true_sd_ancova The population error standard deviation of the ANCOVA model
#' @param true_std_contrast The population standardized ANCOVA (adjusted) contrast
#' @param true_sd_ratio The ratio of \code{true_sd_ancova} over \code{true_sd_anova}, greater than 0 and at most 1
#' @param true_cor The population correlation coefficient between the response and the covariate
#' @param divisor Which error standard deviation to be used in standardizing the contrast; the value can be either \code{"s_ancova"} or \code{"s_anova"}
#' @param c_weights Contrast weights
#' @param width The desired full width of the obtained confidence interval
#' @param conf_level The desired confidence interval coverage (i.e., 1 - Type I error rate). Default is \code{.95}, which gives a symmetric two-sided interval, the interval \code{\link{ci_sc_ancova}} reports. A level above .999999998, whose tail areas are below 1e-9, is refused.
#' @param alpha_lower,alpha_upper The lower-tail and upper-tail Type I error
#'   rates, an alternative to \code{conf_level} (give \code{conf_level =
#'   NULL} with them). \code{\link{ci_sc_ancova}} reports intervals with
#'   equal tail areas only, so the two must be equal, and
#'   \code{alpha_lower = alpha_upper = a} plans the interval of
#'   \code{conf_level = 1 - 2a}. Unequal values are refused with a
#'   condition of class \code{dmar_argument_invalid}: they would plan for
#'   an interval no DMAR function reports. Each must be at least 1e-9:
#'   below it the noncentral \emph{t} limits lose accuracy (see "Small tail
#'   areas" in \code{\link{ci_nc_t}}).
#' @param assurance Parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be \code{NULL} or greater than .5 and less than 1). An assurance above .999999998 is refused, with either \code{assurance_method}: the published assurance step finds noncentral \emph{t} limits with tail areas of \code{(1 - assurance) / 2}, which must be at least 1e-9.
#' @param assurance_method How the plan with \code{assurance} is made, used
#'   only with \code{assurance}. \code{"analytic"} (the default) is the
#'   published plan of Lai and Kelley (2012), computed without simulation,
#'   which does not compute the probability it attains and can attain less
#'   than \code{assurance}; \code{assurance_attained} is then \code{NA}.
#'   \code{"simulation"} returns the smallest sample size per group, at or
#'   above the plan without \code{assurance}, at which the probability that
#'   the \code{\link{ci_sc_ancova}} interval is no wider than \code{width}
#'   is at least \code{assurance}, with that probability estimated from
#'   \code{G} simulated samples and reported as \code{assurance_attained}
#'   (see Details).
#' @param G Number of simulated samples at each sample size the assurance
#'   step tries, and in its estimate of the attained probability; a whole
#'   number of at least 100 (default 10000, which gives a Monte Carlo
#'   standard error of .004 for a probability of .80). Used and checked only
#'   with \code{assurance} and \code{assurance_method = "simulation"}.
#' @param seed Optional single number. When given, the assurance step is
#'   reproducible, and the random number stream of the calling session is
#'   left as it was. With \code{NULL} (the default) the draws come from the
#'   session's current stream. Used only where \code{G} is.
#'
#' @details
#' The sample size planning method this function is based on is developed in the context of simple (i.e., one-response-one-covariate)
#' ANCOVA model and randomized design (i.e., same population covariate mean across groups).
#'
#' An ANCOVA contrast can be standardized in at least two ways: (a) divided by the error standard deviation of the
#' ANOVA model, (b) divided by the error standard deviation of the ANCOVA model. This function can be used to analyze
#' both types of standardized ANCOVA contrasts.
#'
#' Not all of the effect size arguments need to be specified. When
#' \code{divisor="s_ancova"} the input is either (a) \code{true_std_contrast},
#' or (b) \code{true_contrast} (the unstandardized ANCOVA contrast) and
#' \code{true_sd_ancova}. When \code{divisor="s_anova"}, the valid input
#' combinations are (a) \code{true_std_contrast} and \code{true_sd_ratio};
#' (b) \code{true_std_contrast} and \code{true_cor}; or
#' (c) \code{true_contrast}, \code{true_sd_anova}, and \code{true_sd_ancova}.
#' With \code{divisor = "s_anova"} the ratio of the two standard deviations,
#' \eqn{\sqrt{1 - \rho^2}}, must be greater than 0 and at most 1, as it is in
#' every population (a covariate cannot raise the error variance), so
#' \code{true_cor} of 1 or -1 is refused there.
#'
#' \strong{The interval planned for.} The plan is for the interval of
#' \code{\link{ci_sc_ancova}}: the limits for the noncentrality parameter of
#' a noncentral \emph{t} on \eqn{N - J - 1} degrees of freedom at the
#' observed \emph{t} of the adjusted contrast, carried back to the scale of
#' the standardized contrast by the standard error
#' \eqn{\sqrt{C/n + (\sum_j c_j \bar{X}_j)^2 / SS_{W,X}}}, with
#' \eqn{C = \sum_j c_j^2}, \eqn{\bar{X}_j} the covariate means, and
#' \eqn{SS_{W,X}} the within-group sum of squares of the covariate (and, for
#' the divisor \code{"s_anova"}, by the ratio \eqn{s_{ANCOVA}/s_{ANOVA}}).
#' The plan without \code{assurance} is for the criterion that this
#' interval, evaluated at the planning value, is no wider than \code{width}
#' (Lai & Kelley, 2012, Equations 21 and 28); it uses
#' \eqn{\sigma_{ANCOVA}\sqrt{C/n}} as the standard error of the contrast,
#' leaving out the covariate means term, which in a randomized design is
#' small (Lai & Kelley, 2012, pp. 354--355). The plan is the smallest
#' \emph{n} per group that meets the criterion, searched from 2 per group,
#' the smallest sample size with positive error degrees of freedom. The
#' search doubles a step and then bisects, which finds the smallest
#' \emph{n} because the interval at the planning value narrows as \emph{n}
#' grows from 2 (checked at every \emph{n} up to 80, and in steps of 10 up
#' to 600, for 420 combinations of the contrast, the ratio
#' \eqn{\sigma_{ANCOVA}/\sigma_{ANOVA}}, weights for two to five groups,
#' and the confidence level). At \code{true_std_contrast = 0},
#' \code{c_weights = c(1, -1)}, and \code{width = .5} the plan is 123 per
#' group, where the interval is .49985 wide (.50190 at 122), and at
#' \code{true_std_contrast = .5}, \code{true_cor = .5},
#' \code{divisor = "s_anova"}, \code{c_weights = c(1, -1)}, and
#' \code{width = 3.6} it is 2 per group, where the interval is 3.596 wide.
#' A \code{width} that would need more than \eqn{2^{53}} (about \eqn{9.0 \times 10^{15}}, the largest number up to which R holds every whole number exactly) per group is refused
#' with a condition of class \code{dmar_argument_invalid}. The intervals
#' have equal tail areas, as those of \code{\link{ci_sc_ancova}} do and as
#' Lai and Kelley (2012, p. 354) assume.
#'
#' Past about \eqn{10^{11}} per group the plan and
#' \code{assurance_attained} are limited by the precision of doubles, as
#' for \code{\link{ss_aipe_sc}}: the width changes by about \eqn{1/(2n)}
#' of itself from one \emph{n} to the next, and for a nonzero contrast it
#' is the difference of two limits near the noncentrality. Against a
#' 512-bit computation of the interval, with \code{c_weights = c(1, -1)},
#' \code{conf_level = .95}, and the divisor \code{"s_ancova"}, and with
#' \code{width} and the tail areas as R holds them (each tail area,
#' (1 - .95) / 2, is 0.025000000000000022 in doubles): at
#' \code{true_std_contrast = .5} and \code{width} of 1e-5, 1e-6, and 1e-7
#' the plan was 4, 1,035, and 525,379 above the smallest \emph{n} that
#' meets the criterion (about 3.2e11, 3.2e13, and 3.2e15 per group). At a
#' contrast of 0, where the width is computed from the normal quantiles,
#' the plan was the smallest \emph{n} at \code{width} of 1e-3, 1e-5, and
#' 1e-6 (about 3.1e7, 3.1e11, and 3.1e13 per group) and 2 below it at 1e-7
#' (about 3.1e15). Past about \eqn{10^{15}} per group the plan at a
#' contrast of 0 can be below the smallest \emph{n}, by up to about
#' \eqn{10^{-15} n}: at \code{width = 6e-8} it is 8,536,575,157,098,048
#' per group, 8 below the smallest \emph{n}, 8,536,575,157,098,056. With
#' \code{assurance = .9}, \code{assurance_method = "simulation"},
#' \code{G = 200}, and \code{seed = 1}, \code{assurance_attained} was the
#' share of the same simulated
#' statistics whose interval, computed in 512 bits, is no wider than
#' \code{width}, up to 3.1e11 per group; at about 3.1e13 per group it was
#' .950 where that share is .940 (contrast 0) and .900 where it is .895
#' (contrast .5), and at about 3.1e15 it was 1.000 where it is 0 (the plan
#' is below the smallest \emph{n}) and .905 where it is .930.
#'
#' \strong{The assurance step.} With \code{assurance} the plan is by
#' default (\code{assurance_method = "analytic"}) the published one (Lai &
#' Kelley, 2012, pp. 358--359, a generalization of Kelley & Rausch, 2006):
#' the value \eqn{t^*} that \eqn{|T|} does not exceed with probability
#' \code{assurance} under the planner's model, carried back to the scale of
#' the standardized contrast, is planned for without assurance. It is
#' computed without simulation and does not compute the probability it
#' attains, so \code{assurance_attained} is \code{NA}, and the printed
#' footer says so. That model leaves out the covariate means term and
#' treats \eqn{s_{ANCOVA}/s_{ANOVA}} as known, and the widths vary little,
#' so the plan attains less than the request: at
#' \code{true_std_contrast = .5}, \code{true_cor = .5},
#' \code{c_weights = c(1, -1)}, \code{width = .5}, and
#' \code{assurance = .80}, the published plans (129 per group with
#' \code{divisor = "s_ancova"}, 98 with \code{"s_anova"}) give an interval
#' no wider than .5 with probability .77 and .55.
#'
#' With \code{assurance_method = "simulation"} the plan is the smallest
#' sample size per group, not below the plan without \code{assurance}, at
#' which the
#' \code{\link{ci_sc_ancova}} interval is no wider than \code{width} with at
#' least that probability, for normal within-group errors and a normal
#' covariate with the same distribution in every group (the randomized
#' design). The width depends on the sample through five independent
#' statistics: the contrast of the covariate means, the within-group sum of
#' squares of the covariate, the contrast of the group means of the errors,
#' the error of the pooled within-group slope, and the ANCOVA error sum of
#' squares. The probability is estimated from \code{G} draws of these
#' statistics from their exact sampling distributions, the same draws at
#' every sample size the search tries, so that the estimate changes
#' smoothly as the sample size grows, and both the covariate means term
#' and, with \code{divisor = "s_anova"}, the sampling of
#' \eqn{s_{ANCOVA}/s_{ANOVA}} enter it. \code{assurance_attained} is the
#' estimate at the returned sample size, at least \code{assurance}, with its
#' Monte Carlo standard error in the attribute \code{assurance_mc_se}. Give
#' \code{seed} for a plan that can be reproduced. For the divisor
#' \code{"s_anova"}, Lai and Kelley (2012, pp. 360--361) wrote that their
#' assurance method is not appropriate, because the width varies with
#' \eqn{s_{ANCOVA}/s_{ANOVA}} as well as with \emph{t}, and gave as a
#' practical solution an a priori Monte Carlo simulation that raises the
#' sample size until the proportion of intervals no wider than the target
#' reaches the assurance; this step is that simulation, made from the
#' sufficient statistics rather than from raw data, and applied to both
#' divisors. In the example above the simulated plans (130 and 104 per
#' group, \code{seed = 1}) give an interval no wider than .5 with
#' probability .87 and .84. Because the probability is estimated, a plan
#' can differ by one from the smallest sample size at which the probability
#' itself reaches \code{assurance}: with \code{"s_anova"} above, 103 gives
#' .803. \code{\link{ss_aipe_sc_ancova_sensitivity}} simulates raw data at
#' a plan, and the attained probability of any plan can be checked there.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}:
#' \code{necessary_n_per_group}, the per-group sample size needed for the
#' planned ANCOVA contrast; with \code{assurance}, the rows
#' \code{assurance} (the request) and \code{assurance_attained}: \code{NA}
#' for the published plan of \code{assurance_method = "analytic"}, which
#' does not compute it, and for \code{assurance_method = "simulation"} the
#' simulated probability at \code{necessary_n_per_group} that the
#' \code{\link{ci_sc_ancova}} interval is no wider than \code{width}, with
#' its Monte Carlo standard error in the attribute \code{assurance_mc_se}
#' and the row \code{G}. With \code{assurance} the printed footer names the
#' assurance method and says how \code{assurance_attained} was found, or
#' that the published plan does not compute it; the attributes
#' \code{assurance_method} and \code{assurance_probability}
#' (\code{"simulation"} or \code{"not computed"}) record the same.
#'
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized
#'   effect sizes: Theory, application, and implementation.
#'   \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals.
#'   \emph{Psychological Methods, 11}(4), 363--385.
#'   \doi{10.1037/1082-989X.11.4.363}
#'
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
#'   ANCOVA and ANOVA contrasts: Sample size planning via narrow
#'   confidence intervals.
#'   \emph{British Journal of Mathematical and Statistical Psychology, 65},
#'   350--370. \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9.)
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical models. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_aipe_sc}}, \code{\link{ss_aipe_sc_ancova_sensitivity}}
#'
#' @examples
#' ss_aipe_sc_ancova(true_std_contrast = .8, width = .5, c_weights = c(.5, .5, 0, -1))
#'
#' ss_aipe_sc_ancova(true_std_contrast = .8, true_sd_ratio = .6, width = .5,
#'                   c_weights = c(.5, .5, 0, -1), divisor = "s_anova")
#'
#' ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .4, width = .3,
#'                c_weights = c(.5, .5, 0, -1), divisor = "s_anova")
#'
#' # With assurance, the published plan of Lai and Kelley (2012), which
#' # does not compute the probability it attains.
#' ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5, width = .5,
#'                   c_weights = c(1, -1), assurance = .90)
#'
#' # By simulation: the smallest sample size per group at which the
#' # interval is no wider than .5 with probability .90. G = 2000 keeps the
#' # example quick; a reported plan deserves the default G = 10000.
#' plan <- ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5, width = .5,
#'                           c_weights = c(1, -1), assurance = .90,
#'                           assurance_method = "simulation",
#'                           G = 2000, seed = 113)
#' plan
#' attr(plan, "assurance_mc_se")
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_sc_ancova <- function(true_contrast = NULL, true_sd_anova = NULL, true_sd_ancova = NULL, true_std_contrast = NULL, true_sd_ratio = NULL, true_cor = NULL, divisor = "s_ancova", c_weights, width,
                              conf_level = .95,
                              alpha_lower = NULL, alpha_upper = NULL,
                              assurance = NULL,
                              assurance_method = c("analytic", "simulation"),
                              G = 10000, seed = NULL) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  if (!missing(width)) .dmar_check_width(width)
  .dmar_check_assurance_single(assurance)
  # An assurance in (.5, 1), as in ss_aipe_sc(): no sample size attains 1.
  if (!is.null(assurance)) {
    if ((assurance <= 0) || (assurance >= 1)) .dmar_refuse_argument(
      "The 'assurance' must either be NULL or some value greater than zero and less than unity.",
      arg = "assurance")
    if (assurance <= .50) .dmar_refuse_argument(
      "The 'assurance' should be larger than 0.5 (but less than 1).",
      arg = "assurance")
    # The published step (assurance_method = "analytic") finds noncentral t
    # limits with tail areas of (1 - assurance) / 2 and 1 - assurance, held
    # to the floor on a tail area (R/dmar_conditions.R). The check is made
    # for either method, so that the two refuse the same assurances.
    .dmar_check_assurance_tail_floor(assurance)
  }
  .dmar_check_numbers("ss_aipe_sc_ancova", "true_std_contrast")
  # "analytic" (the default), the published plan, or "simulation"
  # (R/ss_aipe_assurance_internals.R); used only with assurance.
  assurance_method <- .ss_aipe_assurance_method(
    assurance_method, "ss_aipe_sc_ancova", simulation = TRUE)
  if (divisor != "s_ancova" && divisor != "s_anova") .dmar_refuse_argument(
    "The argument 'divisor' must be either 's_ancova' or 's_anova'",
    arg = "divisor")
  if (!is.null(true_sd_ratio)) {
    if (true_sd_ratio > 1 || true_sd_ratio < 0) .dmar_refuse_argument(
      "'true_sd_ratio' must be larger than 0 and smaller than or equal to 1.",
      arg = "true_sd_ratio")
    if (!is.null(true_cor)) .dmar_refuse_argument(
      "You just need to specify either 'true_sd_ratio' or 'true_cor', not both.",
      arg = c("true_sd_ratio", "true_cor"))
  }
  if (!is.null(true_cor)) {
    if (true_cor > 1 || true_cor < -1) .dmar_refuse_argument(
      "'true_cor' must be no smaller than -1 and no larger than 1.",
      arg = "true_cor")
    true_sd_ratio <- sqrt(1 - true_cor^2)
  }

  # Resolve interval bounds (mirrors ss_aipe_sc / ci_nc_t contract):
  # supply either a symmetric conf_level (default .95) or a pair
  # alpha_lower/alpha_upper, but not both.
  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (alphas_supplied) {
    if (!missing(conf_level) && !is.null(conf_level)) {
      .dmar_refuse_argument(
        "Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'; you cannot mix them.",
        arg = c("conf_level", "alpha_lower", "alpha_upper"))
    }
    if (is.null(alpha_lower) || is.null(alpha_upper)) {
      .dmar_refuse_argument(
        "Supply both 'alpha_lower' and 'alpha_upper' together.",
        arg = c("alpha_lower", "alpha_upper"))
    }
    if (alpha_lower <= 0 || alpha_upper <= 0 || alpha_lower + alpha_upper >= 1) {
      .dmar_refuse_argument(
        "'alpha_lower' and 'alpha_upper' must be strictly positive and sum to less than 1.",
        arg = c("alpha_lower", "alpha_upper"))
    }
    # ci_sc_ancova() reports equal tails only (Lai & Kelley, 2012, p. 354).
    if (!isTRUE(all.equal(alpha_lower, alpha_upper, tolerance = 1e-10))) {
      .dmar_refuse_argument(c(
        "'alpha_lower' and 'alpha_upper' must be equal: ss_aipe_sc_ancova() ",
        "plans for the interval of ci_sc_ancova(), which has equal tail ",
        "areas (it takes conf_level only), so a plan for unequal tails would ",
        "be for an interval no DMAR function reports. Give conf_level, or ",
        "alpha_lower = alpha_upper."),
        arg = c("alpha_lower", "alpha_upper"),
        alpha_lower = alpha_lower, alpha_upper = alpha_upper)
    }
  } else {
    if (is.null(conf_level) || !is.numeric(conf_level) || length(conf_level) != 1 ||
        conf_level <= 0 || conf_level >= 1) {
      .dmar_refuse_argument(
        "'conf_level' must be a single number strictly between 0 and 1.",
        arg = "conf_level")
    }
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }
  # The floor on a tail area (R/dmar_conditions.R).
  .dmar_check_tail_floor(alpha_lower, alpha_upper,
                         conf_level = if (!alphas_supplied) conf_level,
                         zero_ok = FALSE)

  if (abs(sum(c_weights)) > 1e-8) .dmar_refuse_argument(
    "'c_weights' must sum to zero.", arg = "c_weights")
  # All-zero weights pass both rules and define no contrast; they stopped
  # with "missing value where TRUE/FALSE needed" (ss_aipe_sc_ancova()) or a
  # message about 'ncp' (ci_sc_ancova()).
  .dmar_check_contrast_nonzero(c_weights)
  if (sum(c_weights[c_weights > 0]) > 1) .dmar_refuse_argument(
    "The positive values of 'c_weights' must not sum to more than 1; give the weights as fractions, such as c(1/2, 1/2, -1).", arg = "c_weights")

  # The planning value on the scale of the divisor and the ratio
  # sigma_ancova / sigma_anova (1 for the divisor s_ancova).
  if (divisor == "s_ancova") {
    if (is.null(true_std_contrast)) {
      if (is.null(true_contrast) || is.null(true_sd_ancova)) .dmar_refuse_argument(
        "You must specify either 'true_std_contrast', or both 'true_contrast' and 'true_sd_ancova'.",
        arg = c("true_std_contrast", "true_contrast", "true_sd_ancova"))
      true_std_contrast <- true_contrast / true_sd_ancova
    }
    ratio <- 1
  } else {
    if (!is.null(true_contrast) || !is.null(true_sd_anova) || !is.null(true_sd_ancova)) {
      if (is.null(true_contrast) || is.null(true_sd_anova) || is.null(true_sd_ancova)) .dmar_refuse_argument(
        "'true_contrast', 'true_sd_anova', and 'true_sd_ancova' must be all specified.",
        arg = c("true_contrast", "true_sd_anova", "true_sd_ancova"))
    }
    if (!is.null(true_std_contrast) || !is.null(true_sd_ratio)) {
      if (is.null(true_std_contrast) || is.null(true_sd_ratio)) .dmar_refuse_argument(
        "Both 'true_std_contrast' and 'true_sd_ratio' (or 'true_cor') must be specified",
        arg = c("true_std_contrast", "true_sd_ratio", "true_cor"))
    }
    if (is.null(true_sd_ratio)) {
      true_sd_ratio <- true_sd_ancova / true_sd_anova
      true_std_contrast <- true_contrast / true_sd_anova
    }
    # sqrt(1 - rho^2) is in (0, 1] in every population.
    if (!(true_sd_ratio > 0 && true_sd_ratio <= 1)) .dmar_refuse_argument(c(
      "With divisor = \"s_anova\" the ratio of the ANCOVA to the ANOVA error ",
      "standard deviation, sqrt(1 - rho^2), must be greater than 0 and at ",
      "most 1; ", format(true_sd_ratio), " was given or implied by ",
      "'true_sd_ratio', 'true_cor', or 'true_sd_ancova' / 'true_sd_anova'."),
      arg = c("true_sd_ratio", "true_cor", "true_sd_ancova", "true_sd_anova"),
      value = true_sd_ratio)
    ratio <- true_sd_ratio
  }

  n_0 <- .sc_ancova_plan_width(true_std_contrast, ratio, c_weights, width,
                               alpha_lower, alpha_upper)
  if (is.null(assurance)) {
    return(.as_dmar_tbl(data.frame(term = "necessary_n_per_group", value = n_0),
                        conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }

  # The published plan of Lai and Kelley (2012), assurance_method =
  # "analytic", which does not compute the probability it attains; with
  # equal tails a contrast and its negative get one plan.
  if (assurance_method == "analytic") {
    n <- .sc_ancova_published_assurance(true_std_contrast, ratio, c_weights,
                                        width, alpha_lower, alpha_upper,
                                        assurance)
    out <- .as_dmar_tbl(
      data.frame(term = c("necessary_n_per_group", "assurance",
                          "assurance_attained"),
                 value = c(n, assurance, NA_real_)),
      conf_level = conf_level, subclass = "dmar_ss_aipe")
    return(.ss_aipe_assurance_record(out, "analytic", "not computed"))
  }

  G <- .dmar_check_count(G, "G", min = 100L, reason = paste0(
    "The assurance step estimates a probability from 'G' simulated ",
    "samples at each sample size it tries; with 100 the standard error ",
    "of that estimate is at most .05."))
  .dmar_local_seed(seed)
  # The search starts at the plan without assurance for |contrast|.
  if (true_std_contrast < 0)
    n_0 <- .sc_ancova_plan_width(-true_std_contrast, ratio, c_weights, width,
                                 alpha_lower, alpha_upper)
  # The simulated step.
  step <- .sc_ancova_assurance_n(
    n_0, c_weights, delta_e = true_std_contrast / ratio, divisor = divisor,
    ratio = ratio, width = width, alpha_lower = alpha_lower,
    alpha_upper = alpha_upper, assurance = assurance, G = G)
  out <- .as_dmar_tbl(
    data.frame(term = c("necessary_n_per_group", "assurance",
                        "assurance_attained", "G"),
               value = c(step$n, assurance, step$attained, G)),
    conf_level = conf_level, subclass = "dmar_ss_aipe")
  attr(out, "assurance_mc_se") <- step$mc_se
  .ss_aipe_assurance_record(out, "simulation", "simulation")
}
