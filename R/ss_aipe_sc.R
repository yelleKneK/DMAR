#' Sample Size Planning for Accuracy in Parameter Estimation (AIPE) of the Standardized Contrast in ANOVA
#'
#' @description
#' Plans the sample size per group so that the confidence interval for a
#' standardized contrast of means in a fixed effects analysis of variance,
#' the interval computed by \code{\link{ci_sc}}, is sufficiently narrow, an
#' application of the accuracy in parameter estimation (AIPE) approach to the
#' standardized contrast.
#'
#' @param true_std_contrast Population standardized contrast
#' @param c_weights The contrast weights, as \code{\link{ci_sc}} takes them:
#'   fractions whose positive values sum to 1 and whose sum is zero (for
#'   example \code{c(1/2, 1/2, -1)}, not \code{c(1, 1, -2)}), so that the
#'   standardized contrast is a difference of averages of means in units of
#'   the error standard deviation. Weights that break the rule are refused
#'   with the message of \code{\link{ci_sc}} and a condition of class
#'   \code{dmar_argument_invalid}. The rule is applied to the values, so
#'   integer weights such as \code{-1:1} are accepted, as they are by
#'   \code{\link{ci_sc}}.
#' @param width The desired full width of the obtained confidence interval
#' @param conf_level The desired confidence interval coverage (i.e., 1 - Type I error rate). Default is \code{.95}, which gives a symmetric two-sided interval. Specify either \code{conf_level} or both of \code{alpha_lower} and \code{alpha_upper}, not both. A level above .999999998, whose tail areas are below 1e-9, is refused.
#' @param alpha_lower Lower-tail Type I error rate, used to plan an asymmetric confidence interval. When supplied together with \code{alpha_upper}, the planned interval has lower-tail probability \code{alpha_lower} and upper-tail probability \code{alpha_upper}. Set \code{conf_level = NULL} when supplying these.
#' @param alpha_upper Upper-tail Type I error rate, used together with \code{alpha_lower} to plan an asymmetric confidence interval. Each tail area must be at least 1e-9: below it the noncentral \emph{t} limits lose accuracy (see "Small tail areas" in \code{\link{ci_nc_t}}).
#' @param assurance Parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or greater than .5 and less than 1). An assurance above .999999998 is refused: the floor of 1e-9 on a tail area (see "Small tail areas" in \code{\link{ci_nc_t}}) applies to \code{(1 - assurance) / 2}.
#' @param assurance_method How the probability that the
#'   \code{\link{ci_sc}} interval is no wider than \code{width} is found,
#'   used only with \code{assurance}. Only \code{"analytic"} (the default)
#'   is offered here: the probability is computed exactly, and the plan is
#'   the smallest sample size per group, not below the plan without
#'   \code{assurance}, at which it reaches \code{assurance}.
#'   \code{"simulation"}, which \code{\link{ss_aipe_sc_ancova}},
#'   \code{\link{ss_aipe_reg_coef}} with \code{noncentral = TRUE}, and
#'   \code{\link{ss_aipe_src}} offer, where the probability has no exact
#'   form, is refused with a condition of class
#'   \code{dmar_argument_invalid}, whether or not \code{assurance} is given,
#'   since a simulated estimate of an exact probability would only add Monte
#'   Carlo error.
#'
#' @details
#' The plan is for the interval of \code{\link{ci_sc}}, the limits for the
#' noncentrality parameter of a noncentral \emph{t} on \eqn{J(n - 1)}
#' degrees of freedom at the observed \emph{t} times
#' \eqn{\sqrt{\sum_j c_j^2 / n}}. Without \code{assurance} it is the
#' smallest \emph{n} per group, and at least 4, at which that interval,
#' evaluated at the planning value, is no wider than \code{width}; a
#' \code{width} that would need more than \eqn{2^{53}} (about \eqn{9.0 \times 10^{15}}, the largest number up to which R holds every whole number exactly) per group is refused
#' with a condition of class \code{dmar_argument_invalid}.
#'
#' With \code{assurance} the plan is the smallest sample size per group,
#' not below the plan without assurance, at which the probability that the
#' interval is no wider than \code{width} is at least \code{assurance}
#' (\code{assurance_method = "analytic"}), and that probability, computed
#' exactly, is \code{assurance_attained}. The observed \emph{t} is
#' noncentral \emph{t} on \eqn{J(n - 1)} degrees of freedom with
#' noncentrality \eqn{\delta / \sqrt{\sum_j c_j^2 / n}}, the width is
#' smallest at \eqn{t = 0} and rises on either side, so the event is
#' \eqn{t_L \le T \le t_U} for the two values at which the width equals the
#' target, and its probability is that of the noncentral \emph{t}. The plan
#' is found by doubling a step upward from the plan without assurance and
#' then bisecting, which assumes that the probability does not decrease as
#' \emph{n} grows from there. In 96 planning sets (contrasts of 0, .5, -1,
#' and 2; weights for two, three, and five groups; widths of .5 and 1;
#' assurances of .80 and .99; and tail areas of .025 and .025 or .01 and
#' .04) it did not decrease at any \emph{n} from the plan without
#' assurance to ten past the plan, and the plan was the smallest \emph{n}
#' found by raising \emph{n} one at a time.
#'
#' The step of Kelley and Rausch (2006), as Lai and Kelley (2012,
#' pp. 358--359 and 361) apply it to a standardized contrast, plans
#' without assurance for the value \eqn{t^*} that \eqn{|T|} does not
#' exceed with probability \code{assurance} at the plan without
#' assurance, carried back to the scale of the standardized contrast; it
#' does not compute the probability. That step can fall short of the
#' request: at \code{true_std_contrast = 0}, \code{c_weights = c(1, -1)},
#' \code{width = .8}, and \code{assurance = .95} it gives 49 per group,
#' which attains .9499, where the plan here is 50, which attains .9946. It
#' can also be larger than the smallest plan that attains the request: in
#' the 96 sets above it attained the request in every set and was one case
#' or more larger in 17.
#'
#' Past about \eqn{10^{11}} per group the plan and
#' \code{assurance_attained} are limited by the precision of doubles. The
#' width changes by about \eqn{1/(2n)} of itself from one \emph{n} to the
#' next, and for a nonzero contrast it is the difference of two limits near
#' the noncentrality \eqn{\delta / \sqrt{\sum_j c_j^2 / n}}, which a
#' double holds to about \eqn{1.1 \times 10^{-16}} of its size. Against a
#' 512-bit computation of the interval, with \code{c_weights = c(1, -1)}
#' and \code{conf_level = .95}, and with \code{width} and the tail areas
#' as R holds them (each tail area, (1 - .95) / 2, is
#' 0.025000000000000022 in doubles): the plans checked up to 4.6e9 per
#' group were the smallest \emph{n}; at \code{true_std_contrast = 2} and
#' \code{width} of 1e-5, 1e-6, 3e-7, and 1e-7 the plan was 2, 1,850,
#' 487,182, and 6,485,787 below the smallest \emph{n} (about 4.6e11,
#' 4.6e13, 5.1e14, and 4.6e15), and at .5 it was 4 above, 1,956 below,
#' 129 above, and 4,854,160 above it (about 3.2e11, 3.2e13, 3.5e14, and
#' 3.2e15). At a contrast of 0 the width is computed from the normal
#' quantiles, and the plan is below the smallest \emph{n} only when the
#' \emph{n} at which the width equals \code{width} is within about
#' \eqn{10^{-15} n} above a whole number, and then by at most the larger
#' of 1 and about \eqn{10^{-15} n}. At \code{width = 1e-6} that \emph{n} is
#' 30,731,670,565,552.9985, and the plan, 30,731,670,565,553 per group, is
#' the smallest \emph{n}; at 3e-7 that \emph{n} is
#' 341,463,006,283,922.21, and the plan, 341,463,006,283,922, is 1 below
#' the smallest \emph{n}; and the plan was 2 below at 1e-7 (about 3.1e15
#' per group), 4 below at 7e-8 (6.3e15), and 8 below at 6e-8 (8.5e15).
#' With \code{assurance = .9}, \code{assurance_attained} was within 1e-6
#' of the exact probability up to 4.6e9 per group.
#'
#' With unequal tails the width at an observed \emph{t} is not the width at
#' \eqn{-t}: the interval at \eqn{-t} with tail areas
#' (\code{alpha_lower}, \code{alpha_upper}) is the reflection of the
#' interval at \emph{t} with the two exchanged. A negative
#' \code{true_std_contrast} therefore gets the plan of its magnitude with
#' \code{alpha_lower} and \code{alpha_upper} exchanged, so that reversing
#' the contrast and the tails together leaves the plan unchanged. Without
#' \code{assurance} the signed planning value already gives that plan.
#' With equal tails the plan does not depend on the sign.
#'
#' @return
#' \item{necessary_n_per_group}{Necessary sample size \emph{per group}}
#' \item{assurance, assurance_attained}{With \code{assurance} only: the
#'   request, and the exact probability at \code{necessary_n_per_group}
#'   that the \code{\link{ci_sc}} interval is no wider than \code{width}
#'   (at least \code{assurance}). The printed footer names the assurance
#'   method, which the attribute \code{assurance_method} records, with
#'   \code{assurance_probability} \code{"exact"}.}
#'
#' @references
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#'   calculation of confidence intervals that are based on central and
#'   noncentral distributions. \emph{Educational and Psychological
#'   Measurement, 61}(4), 532--574. \doi{10.1177/0013164401614002}
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators.
#' \emph{Journal of Educational Statistics, 6}(2), 107--128.
#'
#' Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the
#' standardized mean difference: Bootstrap and parametric confidence intervals,
#' \emph{Educational and Psychological Measurement, 65}, 51--69.
#'   \doi{10.1177/0013164404264850}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
#' ANCOVA and ANOVA contrasts: Sample size planning via narrow
#' confidence intervals.
#' \emph{British Journal of Mathematical and Statistical Psychology, 65},
#' 350--370. \doi{10.1111/j.2044-8317.2011.02029.x}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical models. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{ci_sc}}, \code{\link{ci_nc_t}}, \code{\link{ss_aipe_c}}
#'
#' @examples
#' # Suppose the population standardized contrast is believed to be .6
#' # in some 5-group ANOVA model. The researcher is interested in comparing
#' # the average of means of group 1 and 2 with the average of group 3 and 4.
#'
#' # To calculate the necessary sample size per group such that the width
#' # of 95 percent confidence interval of the standardized
#' # contrast is, with 90 percent assurance, no wider than .4 (the last row
#' # is the exact probability at that sample size):
#'
#' ss_aipe_sc(true_std_contrast=.6, c_weights=c(.5, .5, -.5, -.5, 0), width=.4, assurance=.90)
#'
#' # Asymmetric confidence interval: most of the alpha goes in the upper tail
#' # (e.g., when a one-sided concern dominates). Pass alpha_lower and
#' # alpha_upper instead of conf_level.
#' ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
#'            conf_level = NULL, alpha_lower = .01, alpha_upper = .04)
#'

#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_sc <- function(true_std_contrast, c_weights, width,
                       conf_level = .95,
                       alpha_lower = NULL, alpha_upper = NULL,
                       assurance = NULL, assurance_method = "analytic") {
    .dmar_check_conf_level(conf_level, null_ok = TRUE)
    if (!missing(width)) .dmar_check_width(width)
    .dmar_check_assurance_single(assurance)
    # The probability is exact here, so "simulation" is not offered
    # (R/ss_aipe_assurance_internals.R).
    assurance_method <- .ss_aipe_assurance_method(
      assurance_method, "ss_aipe_sc", simulation = FALSE, both = FALSE)
    .dmar_check_numbers("ss_aipe_sc", "true_std_contrast", required = "true_std_contrast")
    # The weights rule of ci_sc(), whose interval is planned for
    # (R/ss_aipe_sc_internals.R).
    .sc_check_c_weights(c_weights)

    # Resolve interval bounds. Two equivalent specifications: a symmetric
    # conf_level (default .95) OR an explicit pair alpha_lower/alpha_upper
    # (possibly asymmetric). The two cannot be mixed, matching ci_nc_t.
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

    J <- length(c_weights)

    if (is.null(assurance)) {
        # Initial sample size guess from the z distribution. For an asymmetric
        # split the width is bracketed by the upper alpha_lower quantile of
        # the normal on the lower tail and the upper alpha_upper quantile on
        # the upper tail; for a symmetric 95% interval this reduces to the
        # standard 2 * qnorm(.975) / width. The interval at t = 0 has exactly
        # this width, and the width is smallest at t = 0, so no n below n_0
        # meets the target. Each quantile is taken in the upper tail
        # (lower.tail = FALSE), which keeps a small tail area exact where
        # 1 - alpha would round it away.
        z_sum <- qnorm(alpha_lower, lower.tail = FALSE) +
            qnorm(alpha_upper, lower.tail = FALSE)
        n_0 <- sum(c_weights^2) * z_sum^2 / width^2

        n <- ceiling(n_0)

        # To ensure that the initial n is not too big.
        n <- max(4, n - 5)

        # The width of the ci_sc() interval at the planning value with n per
        # group. At a contrast of 0 the limits are the normal quantiles
        # -qnorm(1 - alpha_lower) and qnorm(1 - alpha_upper) at any df
        # (P(T <= 0) is pnorm(-ncp)), and the width is computed from them:
        # ci_nc_t() finds them only to its tolerance, and from about 3e13
        # per group on that error, not the width, could decide the plan (at
        # width 6e-7 the plan was 85,365,751,570,980, one below the
        # smallest, and at width 1e-7 it was 3,073,167,056,555,294, the
        # start, where the interval is 9.5e-16 of its width too wide).
        z_tails <- qnorm(alpha_lower, lower.tail = FALSE) +
            qnorm(alpha_upper, lower.tail = FALSE)
        width_at <- function(n) {
            if (true_std_contrast == 0)
                return(sqrt(sum(c_weights^2) / n) * z_tails)
            lambda <- true_std_contrast / sqrt(sum(c_weights^2) / n)
            lambda_limits <- ci_nc_t(ncp = lambda, df = n * J - J,
                                     conf_level = NULL,
                                     alpha_lower = alpha_lower,
                                     alpha_upper = alpha_upper)
            psi_limit_upper <- lambda_limits[which(lambda_limits$term == 'ci_upper'),2] * sqrt(sum(c_weights^2) / n)
            psi_limit_lower <- lambda_limits[which(lambda_limits$term == 'ci_lower'),2] * sqrt(sum(c_weights^2) / n)
            abs(psi_limit_upper - psi_limit_lower)
        }

        # The smallest n at or above the start whose interval is no wider
        # than `width`, found by doubling a step and bisecting, which gives
        # what raising n by one gives (R/ss_aipe_sc_internals.R). A width
        # that needs more than 2^53 per group (the largest number up to
        # which a double holds every whole number exactly) is refused, as in
        # ss_aipe_sc_ancova().
        n_max <- 2^53
        if (n > n_max) .sc_width_unreachable(width, n_max)
        n <- .sc_smallest_n(n, function(m) width_at(m) - width <= 0, n_max,
                            function() .sc_width_unreachable(width, n_max))
        return(.as_dmar_tbl(data.frame(term = 'necessary_n_per_group', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }

    if (!is.null(assurance)) {
        if ((assurance <= 0) || (assurance >= 1)) .dmar_refuse_argument(
          "The 'assurance' must either be NULL or some value greater than zero and less than unity.",
          arg = "assurance")
        if (assurance <= .50) .dmar_refuse_argument(
          "The 'assurance' should be larger than 0.5 (but less than 1).",
          arg = "assurance")
        # The floor on a tail area (R/dmar_conditions.R) applies to
        # (1 - assurance) / 2, as in ss_aipe_sc_ancova(); an assurance that
        # would need a smaller tail is refused in its own name.
        .dmar_check_assurance_tail_floor(assurance)

        # The ci_sc() interval at -t with tails (a_l, a_u) is the reflection
        # of the interval at t with (a_u, a_l), so a negative contrast is
        # planned as its magnitude with the tails exchanged (see Details;
        # DMAR 1.0.x kept the tails as given).
        tails_given <- c(alpha_lower, alpha_upper)
        reflected <- true_std_contrast < 0
        if (reflected) {
            true_std_contrast <- -true_std_contrast
            alpha_lower <- tails_given[2]
            alpha_upper <- tails_given[1]
        }

        # The plan without assurance, where the search starts.
        n0 <- ss_aipe_sc(true_std_contrast = true_std_contrast,
                         conf_level = NULL,
                         alpha_lower = alpha_lower, alpha_upper = alpha_upper,
                         width = width, c_weights = c_weights,
                         assurance = NULL)[1,2]

        # The smallest n at or above n0 whose exact probability that the
        # ci_sc() interval is no wider than `width` is at least `assurance`,
        # as ss_aipe_smd() plans. The reflected problem has the same
        # probability.
        step <- .sc_assurance_search(n0, function(m) .sc_width_probability(
            m, true_std_contrast, c_weights, width, alpha_lower,
            alpha_upper), assurance)
        out <- .as_dmar_tbl(
            data.frame(term = c("necessary_n_per_group", "assurance",
                                "assurance_attained"),
                       value = c(step$n, assurance, step$attained)),
            conf_level = conf_level, subclass = "dmar_ss_aipe")
        return(.ss_aipe_assurance_record(out, assurance_method, "exact"))
    }
}
