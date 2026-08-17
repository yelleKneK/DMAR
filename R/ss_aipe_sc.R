#' Sample Size Planning for Accuracy in Parameter Estimation (AIPE) of the Standardized Contrast in ANOVA
#'
#' @description
#' Plans the sample size per group so that the confidence interval for a
#' standardized contrast of means in a fixed effects analysis of variance,
#' the interval computed by \code{\link{ci_sc}}, is sufficiently narrow, an
#' application of the accuracy in parameter estimation (AIPE) approach to the
#' standardized contrast.
#'
#' @param psi_standardized Population standardized contrast
#' @param c_weights The contrast weights
#' @param width The desired full width of the obtained confidence interval
#' @param conf_level The desired confidence interval coverage (i.e., 1 - Type I error rate). Default is \code{.95}, which gives a symmetric two-sided interval. Specify either \code{conf_level} or both of \code{alpha_lower} and \code{alpha_upper}, not both.
#' @param alpha_lower Lower-tail Type I error rate, used to plan an asymmetric confidence interval. When supplied together with \code{alpha_upper}, the planned interval has lower-tail probability \code{alpha_lower} and upper-tail probability \code{alpha_upper}. Set \code{conf_level = NULL} when supplying these.
#' @param alpha_upper Upper-tail Type I error rate, used together with \code{alpha_lower} to plan an asymmetric confidence interval.
#' @param assurance Parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity)
#' @param \dots Allows one to potentially include parameter values for inner functions
#'
#' @return
#' \item{necessary_n_per_group}{Necessary sample size \emph{per group}}
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
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{ci_sc}}, \code{\link{conf_limits_nct}}, \code{\link{ss_aipe_c}}
#'
#' @examples
#' # Suppose the population standardized contrast is believed to be .6
#' # in some 5-group ANOVA model. The researcher is interested in comparing
#' # the average of means of group 1 and 2 with the average of group 3 and 4.
#'
#' # To calculate the necessary sample size per group such that the width
#' # of 95 percent confidence interval of the standardized
#' # contrast is, with 90 percent assurance, no wider than .4:
#'
#' ss_aipe_sc(psi_standardized=.6, c_weights=c(.5, .5, -.5, -.5, 0), width=.4, assurance=.90)
#'
#' # Asymmetric confidence interval: most of the alpha goes in the upper tail
#' # (e.g., when a one-sided concern dominates). Pass alpha_lower and
#' # alpha_upper instead of conf_level.
#' ss_aipe_sc(psi_standardized = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
#'            conf_level = NULL, alpha_lower = .01, alpha_upper = .04)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_sc <- function(psi_standardized, c_weights, width,
                       conf_level = .95,
                       alpha_lower = NULL, alpha_upper = NULL,
                       assurance = NULL, ...) {
    if (abs(sum(c_weights)) > 1e-8) stop("The sum of the coefficients must be zero")
    if (sum(c_weights[c_weights > 0]) > 1) stop("Please use fractions to specify the contrast weights")

    # Resolve interval bounds. Two equivalent specifications: a symmetric
    # conf_level (default .95) OR an explicit pair alpha_lower/alpha_upper
    # (possibly asymmetric). The two cannot be mixed, matching conf_limits_nct.
    alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
    if (alphas_supplied) {
        if (!missing(conf_level) && !is.null(conf_level)) {
            stop("Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'; you cannot mix them.", call. = FALSE)
        }
        if (is.null(alpha_lower) || is.null(alpha_upper)) {
            stop("Supply both 'alpha_lower' and 'alpha_upper' together.", call. = FALSE)
        }
        if (alpha_lower <= 0 || alpha_upper <= 0 || alpha_lower + alpha_upper >= 1) {
            stop("'alpha_lower' and 'alpha_upper' must be strictly positive and sum to less than 1.", call. = FALSE)
        }
    } else {
        if (is.null(conf_level) || !is.numeric(conf_level) || length(conf_level) != 1 ||
            conf_level <= 0 || conf_level >= 1) {
            stop("'conf_level' must be a single number strictly between 0 and 1.", call. = FALSE)
        }
        alpha_lower <- (1 - conf_level) / 2
        alpha_upper <- (1 - conf_level) / 2
    }

    J <- length(c_weights)

    if (is.null(assurance)) {
        # Initial sample size guess from the z distribution. For an asymmetric
        # split the width is bracketed by qnorm(1 - alpha_lower) on the lower
        # tail and qnorm(1 - alpha_upper) on the upper tail; for a symmetric
        # 95% interval this reduces to the standard 2 * qnorm(.975) / width.
        z_sum <- qnorm(1 - alpha_lower) + qnorm(1 - alpha_upper)
        n_0 <- sum(c_weights^2) * z_sum^2 / width^2

        n <- ceiling(n_0)

        # To ensure that the initial n is not too big.
        n <- max(4, n - 5)

        # Initial estimate of noncentral value.
        # This is literally the theoretical t-value given psi_standardized and the initial estimate of sample size.
        lambda_0 <- psi_standardized / sqrt(sum(c_weights^2) / n)

        # Initial confidence limits.
        lambda_limits_0 <- conf_limits_nct(ncp = lambda_0, df = n * J - J,
                                            conf_level = NULL,
                                            alpha_lower = alpha_lower,
                                            alpha_upper = alpha_upper)
        psi_limit_upper_0 <- lambda_limits_0[which(lambda_limits_0$term == 'upper_limit'),2] * sqrt(sum(c_weights^2) / n)
        psi_limit_lower_0 <- lambda_limits_0[which(lambda_limits_0$term == 'lower_limit'),2] * sqrt(sum(c_weights^2) / n)

        # Initial full-width for confidence interval.
        Diff_width_Full <- abs(psi_limit_upper_0 - psi_limit_lower_0) - width

        while (Diff_width_Full > 0) {
            n <- n + 1
            lambda <- psi_standardized / sqrt(sum(c_weights^2) / n)
            lambda_limits <- conf_limits_nct(ncp = lambda, df = n * J - J,
                                              conf_level = NULL,
                                              alpha_lower = alpha_lower,
                                              alpha_upper = alpha_upper)

            psi_limit_upper <- lambda_limits[which(lambda_limits$term == 'upper_limit'),2] * sqrt(sum(c_weights^2) / n)
            psi_limit_lower <- lambda_limits[which(lambda_limits$term == 'lower_limit'),2] * sqrt(sum(c_weights^2) / n)

            Current_width <- abs(psi_limit_upper - psi_limit_lower)
            Diff_width_Full <- Current_width - width
        }
        return(.as_dmar_tbl(data.frame(term = 'necessary_n_per_group', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }

    if (!is.null(assurance)) {
        psi_sign <- psi_standardized
        psi_standardized <- abs(psi_standardized)

        if ((assurance <= 0) || (assurance >= 1)) stop("The 'assurance' must either be NULL or some value greater than zero and less than unity.", call. = FALSE)
        if (assurance <= .50) stop("The 'assurance' should be larger than 0.5 (but less than 1).", call. = FALSE)

        n0 <- ss_aipe_sc(psi_standardized = psi_standardized,
                         conf_level = NULL,
                         alpha_lower = alpha_lower, alpha_upper = alpha_upper,
                         width = width, c_weights = c_weights,
                         assurance = NULL, ...)[1,2]

        Lim_2 <- conf_limits_nct(
            ncp = psi_standardized / sqrt(sum(c_weights^2) / n0), df = n0 * J - J,
            conf_level = NULL,
            alpha_upper = (1 - assurance) / 2, alpha_lower = (1 - assurance) / 2
        )
        limit_2_sided <- sqrt(sum(c_weights^2) / n0) * Lim_2[which(Lim_2$term == 'upper_limit'),2]

        Lim_1 <- conf_limits_nct(
            ncp = psi_standardized / sqrt(sum(c_weights^2) / n0), df = n0 * J - J,
            conf_level = NULL,
            alpha_upper = 1 - assurance, alpha_lower = 0
        )
        limit_1_sided <- sqrt(sum(c_weights^2) / n0) * Lim_1[which(Lim_1$term == 'upper_limit'),2]

        determine_limit <- function(current_psi_limit = current_psi_limit, sample_size = n0, psi_standardized = psi_standardized,
                                    assurance = assurance) {
            Less <- pt(
                q = -current_psi_limit / sqrt(sum(c_weights^2) / sample_size),
                df = J * sample_size - J, ncp = psi_standardized / sqrt(sum(c_weights^2) / sample_size)
            )

            Greater <- 1 - pt(
                q = current_psi_limit / sqrt(sum(c_weights^2) / sample_size),
                df = J * sample_size - J, ncp = psi_standardized / sqrt(sum(c_weights^2) / sample_size)
            )

            expected_widths_too_large <- Less + Greater
            return((expected_widths_too_large - (1 - assurance))^2)
        }
        optimize_result <- optimize(
            f = determine_limit, interval = c(limit_1_sided, limit_2_sided),
            psi_standardized = psi_standardized, assurance = assurance
        )

        n <- ss_aipe_sc(
            psi_standardized = optimize_result$minimum,
            conf_level = NULL,
            alpha_lower = alpha_lower, alpha_upper = alpha_upper,
            width = width, c_weights = c_weights,
            assurance = NULL, ...
        )[1,2]
        return(.as_dmar_tbl(data.frame(term = 'necessary_n_per_group', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }
}
