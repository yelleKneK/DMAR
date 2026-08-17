#' Sample Size Planning for Polynomial Change Models in Longitudinal Study
#'
#' @description
#' This function plans sample size with respect to the group-by-time interaction in the context of a
#' longitudinal design with two groups. It plans sample size from the accuracy in parameter estimation
#' (AIPE) perspective, where the goal is to obtain a sufficiently narrow confidence interval for the fixed
#' effect polynomial change coefficient parameter (e.g., linear, quadratic, etc.). The sample size returned
#' can be one such that (a) the expected confidence interval width is sufficiently narrow, or (b) the
#' observed confidence interval will be sufficiently narrow with a specified high degree of assurance
#' (e.g., .99, .95, .90, etc.). This function accompanies Kelley and Rausch (2011).
#'
#' @param variance_trend The variance of the individuals' true change coefficients (i.e., \eqn{\sigma^2_{\upsilon_m}} in Kelley & Rausch, 2011) for the polynomial trend (e.g., linear, quadratic, etc.) of interest
#' @param error_variance The true level one error variance (i.e., \eqn{\sigma^2_{\epsilon}} in Kelley & Rausch, 2011). Either \code{error_variance} or \code{variance_true_minus_estimated_trend} must be supplied; if \code{variance_true_minus_estimated_trend} is given directly, \code{error_variance} may be omitted.
#' @param variance_true_minus_estimated_trend The variance of the difference between the \eqn{m}th true change coefficient minus the \eqn{m}th estimated change coefficient (i.e., \eqn{\sigma^2_{\hat{\pi}_{m} - \pi_{m}}} from Equation 19 in Kelley & Rausch, 2011). When derived from \code{error_variance} this equals \eqn{\sigma^2_{\epsilon} f^{2p} / \sum_t c_{mt}^2}, where \eqn{f} is the frequency, \eqn{p} the polynomial order, and \eqn{\sum_t c_{mt}^2} the sum of squared orthogonal polynomial contrast weights over the measurement occasions. A user who already has this variance may supply it directly and omit \code{error_variance}.
#' @param duration The duration of the study
#' @param frequency The number of times measurement occurs within each unit of time
#' @param width Width of the confidence interval
#' @param conf_level The desired level of confidence for the confidence interval that will be computed at the completion of the study
#' @param trend The polynomial trend (1st-3rd) of interest specified as "linear", "quadratic", or "cubic"
#' @param assurance Value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g, .80, .90, .95)
#'
#' @return
#' A \code{data.frame} (class \code{dmar_tbl}) with a single row,
#' \code{necessary_n_per_group}, giving the necessary number of subjects \emph{per group}
#' (the total study size is twice this value) for the combination of the desired
#' confidence interval width, confidence level, optional assurance, and the
#' population parameters at the specified design.
#'
#' @references
#' Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
#'   applications to multiple regression: Power and accuracy for omnibus
#'   and targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen
#'   (Eds.), \emph{The Sage handbook of social research methods}
#'   (pp. 166--192). Sage.
#'
#' Kelley, K., & Rausch, J. R. (2011). Sample size planning for longitudinal
#'   models: Accuracy in parameter estimation for polynomial change parameters.
#'   \emph{Psychological Methods, 16}(4), 391--405. \doi{10.1037/a0023352}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapters 11, 15.)
#'
#' Raudenbush, S. W., & Liu, X.-F. (2001). Effects of study duration,
#'   frequency of observation, and sample size on power in studies of
#'   group differences in polynomial change. \emph{Psychological
#'   Methods, 6}(4), 387--401. \doi{10.1037/1082-989X.6.4.387}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ss_power_pcm}} for the power analytic analog (planning
#'   to detect the group-by-time change difference rather than to estimate it
#'   precisely) on the same model, and \code{\link{ss_aipe_pcm_sensitivity}} for
#'   a Monte Carlo check of how parameter misspecification affects the plan.
#'
#' @note
#' Like in all formal sample size planning methods that require the value of one or more population
#' parameter(s), if the population parameters are incorrectly specified, there is no guarantee that the
#' sample size this function returns will be accurate. Of course, the further away from the true values,
#' the further away the true sample size will tend to be.
#'
#' The number of timepoints in a study (say \eqn{M}) is defined by \eqn{f \times D + 1}, where \eqn{f} is
#' the frequency and \eqn{D} is the duration.
#'
#' @examples
#' # The examples reproduce the tolerance-of-antisocial-thinking illustration
#' # of Kelley and Rausch (2011, Tables 1 and 2), which draws on the National
#' # Youth Survey data also used by Raudenbush and Liu (2001). The level-one
#' # error variance is 0.0262 and the between-subject slope variance is 0.003.
#' # The planner finds the sample size needed for a confidence interval on the
#' # group-by-time slope difference that is no wider than `width`. The returned
#' # necessary_n_per_group is per group, so the total study size is twice that. Unlike
#' # power analysis, the value of the slope is not needed here: the confidence
#' # interval width does not depend on it.
#'
#' # (1) Expected-width planning. With five measurement occasions
#' #     (M = frequency * duration + 1 = 1 * 4 + 1) and a target width of
#' #     0.025, the expected 95% confidence interval is sufficiently narrow at
#' #     278 subjects per group (Kelley & Rausch, 2011, Table 1, T = 5).
#' ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
#'             duration = 4, frequency = 1, width = 0.025, conf_level = .95)
#'
#' # (2) More measurement occasions sharpen the estimate. Extending the study
#' #     so that M = 10 (duration = 9, frequency = 1) cuts the expected-width
#' #     requirement from 278 to 165 per group (Kelley & Rausch, 2011, Table 1,
#' #     T = 10).
#' ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
#'             duration = 9, frequency = 1, width = 0.025, conf_level = .95)
#'
#' # (3) A wider tolerated interval costs less. Relaxing the target width from
#' #     0.025 to 0.05 at M = 5 drops the requirement from 278 to 71 per group
#' #     (Kelley & Rausch, 2011, Table 1, T = 5).
#' ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
#'             duration = 4, frequency = 1, width = 0.05, conf_level = .95)
#'
#' # (4) Adding an assurance parameter. Requiring 85% assurance that the
#' #     realized confidence interval will be no wider than 0.025 raises the
#' #     M = 5 requirement from 278 to 295 per group (Kelley & Rausch, 2011,
#' #     Table 2, T = 5). Assurance guards against the expected-width plan being
#' #     too small for the particular sample obtained.
#' ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
#'             duration = 4, frequency = 1, width = 0.025, conf_level = .95,
#'             assurance = .85)
#'
#' # (5) A higher assurance costs more. Demanding 99% assurance rather than 85%
#' #     raises the per-group requirement further, from 295 to 316.
#' ss_aipe_pcm(variance_trend = 0.003, error_variance = 0.0262,
#'             duration = 4, frequency = 1, width = 0.025, conf_level = .95,
#'             assurance = .99)
#'
#' @keywords multivariate design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_pcm <- function(variance_trend, error_variance = NULL, variance_true_minus_estimated_trend = NULL, duration, frequency, width, conf_level = .95, trend = "linear", assurance = NULL) {

    # Validate the planning targets at entry so a boundary input errors clearly
    # rather than driving the per-group n down to 1 (nu = 2n - 2 = 0), which
    # sends qt() out of domain and crashes the fixed-point search with a missing
    # value.
    if (missing(width) || !is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0)
        stop("'width' must be a single finite positive number.", call. = FALSE)
    if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
        stop("'conf_level' must be a single number strictly between 0 and 1.", call. = FALSE)
    if (!is.null(assurance) && (!is.numeric(assurance) || length(assurance) != 1L || !is.finite(assurance) || assurance <= 0 || assurance >= 1))
        stop("'assurance' must be NULL or a single number strictly between 0 and 1.", call. = FALSE)
    if (missing(variance_trend) || !is.numeric(variance_trend) || length(variance_trend) != 1L || !is.finite(variance_trend) || variance_trend < 0)
        stop("'variance_trend' must be a single finite non-negative number.", call. = FALSE)

    # The admissible minimum per-group sample size: the confidence interval uses
    # 2n - 2 error degrees of freedom, so n must be at least 2 for a positive
    # df. A target so wide that the closed-form n falls below 2 is met at n = 2;
    # return that rather than a degenerate n = 1 (or a crash).
    min_n <- 2L

    # variance_trend: the true variance of the individuals' change coefficients.
    # error_variance: the level one variance of the errors.
    # variance_true_minus_estimated_trend: the variance of the (true minus
    #   estimated) change coefficient. The planner ordinarily derives it from
    #   error_variance as
    #     variance_true_minus_estimated_trend = error_variance * frequency^(2p)
    #                                            / sum_c2_pm,
    #   the same per-subject sampling variance ss_power_pcm() calls V. A user
    #   who has that converted variance in hand may supply it directly and omit
    #   error_variance.


    M <- frequency * duration + 1



    if (trend == "intercept" || trend == "Intercept" || trend == 0 || trend == "INTERCEPT") {
        Trend <- 0
        K_p <- 1
    }

    if (trend == "linear" || trend == "Linear" || trend == 1 || trend == "LINEAR") {
        Trend <- 1
        K_p <- 1 / 12
    }
    if (trend == "quadratic" || trend == "Quadratic" || trend == 2 || trend == "QUADRATIC") {
        Trend <- 2
        K_p <- 1 / 720
    }
    if (trend == "cubic" || trend == "Cubic" || trend == 3 || trend == "CUBIC") {
        Trend <- 3
        K_p <- 1 / 100800
    }


    # Ratio of factorials on the log scale so it does not overflow: factorial()
    # returns Inf for arguments beyond 170, whereas lfactorial is finite for all
    # M and the exponentiated difference recovers the exact ratio. K_p and this
    # expression are from Raudenbush and Liu (2001).
    sum_c2_pm <- K_p * exp(lfactorial(M + Trend) - lfactorial(M - Trend - 1))

    # Convert the level one error variance to the variance of the (true minus
    # estimated) change coefficient. The minimum-variance estimate of the
    # degree-p change coefficient for one subject has sampling variance
    #   error_variance * frequency^(2p) / sum_c2_pm
    # (Raudenbush & Liu, 2001, Equations 4-6 and 17; Kelley & Rausch, 2011,
    # Equation 19). The frequency^(2p) factor arises because the change
    # coefficient is expressed per unit of time, whereas sum_c2_pm is the sum
    # of squared orthogonal polynomial contrast weights over the M measurement
    # positions; it is the same factor ss_power_pcm() carries in V, and reduces
    # to 1 when frequency == 1 (the case of every Kelley & Rausch table). A user
    # may instead supply the converted variance directly through
    # `variance_true_minus_estimated_trend` and omit `error_variance`.
    if (is.null(error_variance) && is.null(variance_true_minus_estimated_trend))
        stop("Specify either 'error_variance' or ",
             "'variance_true_minus_estimated_trend'.", call. = FALSE)

    error_to_trend_factor <- frequency^(2 * Trend) / sum_c2_pm

    if (!is.null(variance_true_minus_estimated_trend)) {
        # Supplying both is allowed only when they agree; an exact-equality
        # check would reject pairs differing only by floating point noise, so
        # the cross-check uses all.equal().
        if (!is.null(error_variance) &&
            !isTRUE(all.equal(error_variance * error_to_trend_factor,
                              variance_true_minus_estimated_trend)))
            stop("'error_variance' and 'variance_true_minus_estimated_trend' ",
                 "are inconsistent: 'variance_true_minus_estimated_trend' should ",
                 "equal error_variance * frequency^(2p) / (sum of the squared ",
                 "polynomial weights). Supply one or the other, or make them ",
                 "agree.", call. = FALSE)
    } else {
        variance_true_minus_estimated_trend <- error_variance * error_to_trend_factor
    }

    n_i_plus_1 <- M + 2
    dif <- 1
    counter <- 0
    past_ss <- 0
    while (dif != 0) {
        counter <- counter + 1
        past_ss <- c(past_ss, n_i_plus_1)

        n_i <- n_i_plus_1
        nu_i <- 2 * n_i - 2
        critival_t_i <- qt(df = nu_i, (1 - (1 - conf_level) / 2))
        # Floor the candidate per-group n at the admissible minimum so nu never
        # drops to zero on the next pass (which would send qt() out of domain);
        # a very wide target settles at min_n rather than a degenerate n = 1.
        n_i_plus_1 <- max(ceiling((8 * (variance_trend + variance_true_minus_estimated_trend) * critival_t_i^2) / (width^2)), min_n)
        dif <- n_i_plus_1 - n_i

        if (counter == 100) {
            n_i_plus_1 <- max(past_ss[95:100])
            dif <- 0
        }
    }

    if (is.null(assurance)) {
        return(.as_dmar_tbl(data.frame(term = 'necessary_n_per_group', value = n_i_plus_1), conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }

    if (!is.null(assurance)) {
        n_i_plus_1 <- n_i_plus_1
        dif <- 1
        counter <- 0
        past_ss <- 0
        while (dif != 0) {
            counter <- counter + 1
            past_ss <- c(past_ss, n_i_plus_1)
            n_i <- n_i_plus_1
            nu_i <- 2 * n_i - 2
            # The critical t tracks the candidate n's degrees of freedom, as in
            # the expected-width loop above; reusing the frozen value from that
            # loop is very slightly conservative but not what the width
            # equation states.
            critival_t_i <- qt(df = nu_i, (1 - (1 - conf_level) / 2))
            critival_chi_square_i <- qchisq(df = nu_i, assurance)
            variance_pi_hat_gamma <- (critival_chi_square_i * (variance_trend + variance_true_minus_estimated_trend)) / nu_i
            # As above, floor at the admissible minimum so nu stays positive.
            n_i_plus_1 <- max(ceiling((8 * (variance_pi_hat_gamma) * critival_t_i^2) / (width^2)), min_n)
            dif <- n_i_plus_1 - n_i


            if (counter == 100) {
                n_i_plus_1 <- max(past_ss[95:100])
                dif <- 0
            }
        }

        return(.as_dmar_tbl(data.frame(term = 'necessary_n_per_group', value = n_i_plus_1), conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }
}
