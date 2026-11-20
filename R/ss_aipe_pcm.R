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
#' @param true_var_trend The variance of the individuals' true change coefficients (i.e., \eqn{\sigma^2_{\upsilon_m}} in Kelley & Rausch, 2011) for the polynomial trend (e.g., linear, quadratic, etc.) of interest
#' @param true_var_error The true within-subject (i.e., Level 1) error variance (\eqn{\sigma^2_{\epsilon}} in Kelley & Rausch, 2011). Either \code{true_var_error} or \code{true_var_trend_error} must be supplied; if \code{true_var_trend_error} is given directly, \code{true_var_error} may be omitted.
#' @param true_var_trend_error The sampling-error variance of one individual's estimated change coefficient: the variance, across individuals, of the difference between the \eqn{m}th estimated and the \eqn{m}th true change coefficient of an individual, \eqn{\sigma^2_{\hat{\pi}_{m} - \pi_{m}}} (Kelley & Rausch, 2011, Equation 19). It is the part of the variance of the estimated change coefficients that measurement error contributes, beside \code{true_var_trend}, the variance of the true ones. Kelley and Rausch (2011, Equation 19) give it as \eqn{\sigma^2_{\epsilon} / \sum_t c_{mt}^2}, where \eqn{\sum_t c_{mt}^2} is the sum of the squared orthogonal polynomial contrast weights over the measurement occasions; with the change coefficient expressed per unit of time it is \eqn{\sigma^2_{\epsilon} f^{2p} / \sum_t c_{mt}^2}, where \eqn{f} is the frequency and \eqn{p} the polynomial order, which is how the function computes it from \code{true_var_error} (the two agree at \eqn{f = 1}). A user who already has this variance may supply it directly and omit \code{true_var_error}.
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
#' # Youth Survey data also used by Raudenbush and Liu (2001). The
#' # within-subject (i.e., Level 1) error variance is 0.0262 and the
#' # between-subject slope variance is 0.003.
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
#' ss_aipe_pcm(true_var_trend = 0.003, true_var_error = 0.0262,
#'             duration = 4, frequency = 1, width = 0.025, conf_level = .95)
#'
#' # (2) More measurement occasions sharpen the estimate. Extending the study
#' #     so that M = 10 (duration = 9, frequency = 1) cuts the expected-width
#' #     requirement from 278 to 165 per group (Kelley & Rausch, 2011, Table 1,
#' #     T = 10).
#' ss_aipe_pcm(true_var_trend = 0.003, true_var_error = 0.0262,
#'             duration = 9, frequency = 1, width = 0.025, conf_level = .95)
#'
#' # (3) A wider tolerated interval costs less. Relaxing the target width from
#' #     0.025 to 0.05 at M = 5 drops the requirement from 278 to 71 per group
#' #     (Kelley & Rausch, 2011, Table 1, T = 5).
#' ss_aipe_pcm(true_var_trend = 0.003, true_var_error = 0.0262,
#'             duration = 4, frequency = 1, width = 0.05, conf_level = .95)
#'
#' # (4) Adding an assurance parameter. Requiring 85% assurance that the
#' #     realized confidence interval will be no wider than 0.025 raises the
#' #     M = 5 requirement from 278 to 295 per group (Kelley & Rausch, 2011,
#' #     Table 2, T = 5). Assurance guards against the expected-width plan being
#' #     too small for the particular sample obtained.
#' ss_aipe_pcm(true_var_trend = 0.003, true_var_error = 0.0262,
#'             duration = 4, frequency = 1, width = 0.025, conf_level = .95,
#'             assurance = .85)
#'
#' # (5) A higher assurance costs more. Demanding 99% assurance rather than 85%
#' #     raises the per-group requirement further, from 295 to 316.
#' ss_aipe_pcm(true_var_trend = 0.003, true_var_error = 0.0262,
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


ss_aipe_pcm <- function(true_var_trend, true_var_error = NULL, true_var_trend_error = NULL, duration, frequency, width, conf_level = .95, trend = "linear", assurance = NULL) {
    .dmar_check_numbers("ss_aipe_pcm", c("duration", "frequency"), required = c("duration", "frequency"))

    # Validate the planning targets at entry so a boundary input errors clearly
    # rather than driving the per-group n down to 1 (nu = 2n - 2 = 0), which
    # sends qt() out of domain and crashes the fixed-point search with a missing
    # value.
    if (missing(width) || !is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0)
        .dmar_refuse_argument(
          "'width' must be a single finite positive number.", arg = "width")
    if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
        .dmar_refuse_argument(
          "'conf_level' must be a single number strictly between 0 and 1.",
          arg = "conf_level")
    if (!is.null(assurance) && (!is.numeric(assurance) || length(assurance) != 1L || !is.finite(assurance) || assurance <= 0 || assurance >= 1))
        .dmar_refuse_argument(
          "'assurance' must be NULL or a single number strictly between 0 and 1.",
          arg = "assurance")
    if (missing(true_var_trend) || !is.numeric(true_var_trend) || length(true_var_trend) != 1L || !is.finite(true_var_trend) || true_var_trend < 0)
        .dmar_refuse_argument(
          "'true_var_trend' must be a single finite non-negative number.",
          arg = "true_var_trend")

    # The admissible minimum per-group sample size: the confidence interval uses
    # 2n - 2 error degrees of freedom, so n must be at least 2 for a positive
    # df. A target so wide that the closed-form n falls below 2 is met at n = 2;
    # return that rather than a degenerate n = 1 (or a crash).
    min_n <- 2L

    # true_var_trend: the true variance of the individuals' change coefficients.
    # true_var_error: the within-subject (Level 1) variance of the errors.
    # true_var_trend_error: the sampling-error variance of one individual's
    #   estimated change coefficient, the variance of (estimated minus true)
    #   change coefficient (Kelley & Rausch, 2011, Equation 19). The planner
    #   ordinarily derives it from true_var_error as
    #     true_var_trend_error = true_var_error * frequency^(2p) / sum_c2_pm,
    #   the same per-subject sampling variance ss_power_pcm() calls V and
    #   reports in its row true_var_trend_error. A user who has that converted
    #   variance in hand may supply it directly and omit true_var_error.


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

    # Convert the within-subject error variance to the variance of the (true minus
    # estimated) change coefficient. The minimum-variance estimate of the
    # degree-p change coefficient for one subject has sampling variance
    #   true_var_error * frequency^(2p) / sum_c2_pm
    # (Raudenbush & Liu, 2001, Equations 4-6 and 17; Kelley & Rausch, 2011,
    # Equation 19). The frequency^(2p) factor arises because the change
    # coefficient is expressed per unit of time, whereas sum_c2_pm is the sum
    # of squared orthogonal polynomial contrast weights over the M measurement
    # positions; it is the same factor ss_power_pcm() carries in V, and reduces
    # to 1 when frequency == 1 (the case of every Kelley & Rausch table). A user
    # may instead supply the converted variance directly through
    # `true_var_trend_error` and omit `true_var_error`.
    if (is.null(true_var_error) && is.null(true_var_trend_error))
        .dmar_refuse_argument(c(
          "Specify either 'true_var_error' or ",
          "'true_var_trend_error'."),
          arg = c("true_var_error", "true_var_trend_error"))

    error_to_trend_factor <- frequency^(2 * Trend) / sum_c2_pm

    if (!is.null(true_var_trend_error)) {
        # Supplying both is allowed only when they agree; an exact-equality
        # check would reject pairs differing only by floating point noise, so
        # the cross-check uses all.equal().
        if (!is.null(true_var_error) &&
            !isTRUE(all.equal(true_var_error * error_to_trend_factor,
                              true_var_trend_error)))
            .dmar_refuse_argument(c(
              "'true_var_error' and 'true_var_trend_error' ",
              "are inconsistent: 'true_var_trend_error' should ",
              "equal true_var_error * frequency^(2p) / (sum of the squared ",
              "polynomial weights). Supply one or the other, or make them ",
              "agree."),
              arg = c("true_var_error", "true_var_trend_error"))
    } else {
        true_var_trend_error <- true_var_error * error_to_trend_factor
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
        n_i_plus_1 <- max(ceiling((8 * (true_var_trend + true_var_trend_error) * critival_t_i^2) / (width^2)), min_n)
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
            variance_pi_hat_gamma <- (critival_chi_square_i * (true_var_trend + true_var_trend_error)) / nu_i
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
