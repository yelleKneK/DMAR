#' Sample Size Planning From the AIPE Perspective for Standardized ANCOVA Contrasts
#'
#' Sample size planning from the accuracy in parameter estimation (AIPE) perspective for standardized ANCOVA contrasts.
#'
#' @param psi The population unstandardized ANCOVA (adjusted) contrast
#' @param sigma_anova The population error standard deviation of the ANOVA model
#' @param sigma_ancova The population error standard deviation of the ANCOVA model
#' @param psi_standardized The population standardized ANCOVA (adjusted) contrast
#' @param ratio The ratio of \code{sigma_ancova} over \code{sigma_anova}
#' @param rho The population correlation coefficient between the response and the covariate
#' @param divisor Which error standard deviation to be used in standardizing the contrast; the value can be either \code{"s_ancova"} or \code{"s_anova"}
#' @param c_weights Contrast weights
#' @param width The desired full width of the obtained confidence interval
#' @param conf_level The desired confidence interval coverage (i.e., 1 - Type I error rate). Default is \code{.95}, which gives a symmetric two-sided interval. Specify either \code{conf_level} or both of \code{alpha_lower} and \code{alpha_upper}, not both.
#' @param alpha_lower Lower-tail Type I error rate, used to plan an asymmetric confidence interval. When supplied together with \code{alpha_upper}, the planned interval has lower-tail probability \code{alpha_lower} and upper-tail probability \code{alpha_upper}. Set \code{conf_level = NULL} when supplying these.
#' @param alpha_upper Upper-tail Type I error rate, used together with \code{alpha_lower} to plan an asymmetric confidence interval.
#' @param assurance Parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be \code{NULL} or between zero and unity)
#' @param \dots Allows one to potentially include parameter values for inner functions
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
#' \code{divisor="s_ancova"} the input is either (a) \code{psi_standardized},
#' or (b) \code{psi} (the unstandardized ANCOVA contrast) and
#' \code{sigma_ancova}. When \code{divisor="s_anova"}, the valid input
#' combinations are (a) \code{psi_standardized} and \code{ratio};
#' (b) \code{psi_standardized} and \code{rho}; or
#' (c) \code{psi}, \code{sigma_anova}, and \code{sigma_ancova}.
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} and \code{value}.
#' The \code{term} is \code{"necessary_n_per_group"} and \code{value} is the
#' per-group sample size needed for the planned ANCOVA contrast.
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
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @note
#' When \code{divisor="s_anova"} and the argument \code{assurance} is specified, the necessary
#' sample size \emph{per group} returned by the function with \code{assurance} specified is slightly underestimated.
#' The method to obtain exact sample size in the above situation has not been developed yet. A practical solution is
#' to use the sample size returned as the starting value to conduct a priori Monte Carlo simulations with
#' function \code{\link{ss_aipe_sc_ancova_sensitivity}}, as discussed in Lai & Kelley (2012).
#'
#' @seealso \code{\link{ss_aipe_sc}}, \code{\link{ss_aipe_sc_ancova_sensitivity}}
#'
#' @examples
#' ss_aipe_sc_ancova(psi_standardized = .8, width = .5, c_weights = c(.5, .5, 0, -1))
#'
#' ss_aipe_sc_ancova(psi_standardized = .8, ratio = .6, width = .5,
#'                   c_weights = c(.5, .5, 0, -1), divisor = "s_anova")
#'
#' ss_aipe_sc_ancova(psi_standardized = .5, rho = .4, width = .3,
#'                c_weights = c(.5, .5, 0, -1), divisor = "s_anova")
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_sc_ancova <- function(psi = NULL, sigma_anova = NULL, sigma_ancova = NULL, psi_standardized = NULL, ratio = NULL, rho = NULL, divisor = "s_ancova", c_weights, width,
                              conf_level = .95,
                              alpha_lower = NULL, alpha_upper = NULL,
                              assurance = NULL, ...) {
  if (divisor != "s_ancova" && divisor != "s_anova") stop("The argument 'divisor' must be either 's_ancova' or 's_anova'")
  if (!is.null(ratio)) {
    if (ratio > 1 || ratio < 0) stop("'ratio' must be larger than 0 and smaller than or equal to 1.")
    if (!is.null(rho)) stop("You just need to specify either 'ratio' or 'rho', not both.")
  }
  if (!is.null(rho)) {
    if (rho > 1 || rho < -1) stop("'rho' must be no smaller than -1 and no larger than 1.")
    ratio <- sqrt(1 - rho^2)
  }

  # Resolve interval bounds (mirrors ss_aipe_sc / conf_limits_nct contract):
  # supply either a symmetric conf_level (default .95) or a possibly asymmetric
  # pair alpha_lower/alpha_upper, but not both.
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

  if (divisor == "s_ancova") {
    if (abs(sum(c_weights)) > 1e-8) stop("The sum of the contrast weights must be zero")
    if (sum(c_weights[c_weights > 0]) > 1) stop("Please use fractions to specify the contrast weights")
    if (is.null(psi_standardized)) {
      if (is.null(psi) || is.null(sigma_ancova)) stop("You must specify either 'psi_standardized', or both 'psi' and 'sigma_ancova'.")
      psi_standardized <- psi / sigma_ancova
    }

    J <- length(c_weights)

    if (is.null(assurance)) { # Initial starting value for n using the z distribution.
      # Asymmetric-aware width formula: each tail contributes its own quantile.
      z_sum <- qnorm(1 - alpha_lower) + qnorm(1 - alpha_upper)
      n_0 <- sum(c_weights^2) * z_sum^2 / width^2

      # Second starting value for n using the central t distribution.
      t_sum_0 <- qt(1 - alpha_lower, n_0 * J - J - 1) + qt(1 - alpha_upper, n_0 * J - J - 1)
      n <- sum(c_weights^2) * t_sum_0^2 / width^2

      # measures the discrepancy between the initial and second starting values.
      Difference <- abs(n - n_0)

      while (Difference > .000001) {
        n_p <- n
        t_sum <- qt(1 - alpha_lower, n * J - J - 1) + qt(1 - alpha_upper, n * J - J - 1)
        n <- sum(c_weights^2) * t_sum^2 / width^2
        Difference <- abs(n - n_p)
      }
      n <- ceiling(n)

      # Initial estimate of noncentral value.
      lambda_0 <- psi_standardized / sqrt(sum(c_weights^2) / n)

      lambda_limits_0 <- conf_limits_nct(ncp = lambda_0, df = n * J - J - 1,
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
        lambda_limits <- conf_limits_nct(ncp = lambda, df = n * J - J - 1,
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
      psi_standardized <- abs(psi_standardized)

      if ((assurance <= 0) || (assurance >= 1)) stop("The 'assurance' must either be NULL or some value greater than zero and less than unity.", call. = FALSE)
      if (assurance <= .50) stop("The 'assurance' should be larger than 0.5 (but less than 1).", call. = FALSE)

      n0 <- ss_aipe_sc_ancova(psi_standardized = psi_standardized,
                              conf_level = NULL,
                              alpha_lower = alpha_lower, alpha_upper = alpha_upper,
                              width = width, c_weights = c_weights,
                              assurance = NULL, ...)[1,2]

      lambda_2_sided <- conf_limits_nct(ncp = psi_standardized / sqrt(sum(c_weights^2) / n0), df = n0 * J - J - 1, conf_level = NULL, alpha_upper = (1 - assurance) / 2, alpha_lower = (1 - assurance) / 2)

      limit_2_sided <- sqrt(sum(c_weights^2) / n0) * lambda_2_sided[which(lambda_2_sided$term == 'upper_limit'),2]

      lambda_1_sided <- conf_limits_nct(ncp = psi_standardized / sqrt(sum(c_weights^2) / n0), df = n0 * J - J - 1, conf_level = NULL, alpha_upper = 1 - assurance, alpha_lower = 0)

      limit_1_sided <- sqrt(sum(c_weights^2) / n0) * lambda_1_sided[which(lambda_1_sided$term == 'upper_limit'),2]

      determine_limit <- function(current_psi_limit = current_psi_limit, sample_size = n0, psi_standardized = psi_standardized, assurance = assurance) {
        Less <- pt(
          q = -current_psi_limit / sqrt(sum(c_weights^2) / sample_size),
          df = J * sample_size - J - 1, ncp = psi_standardized / sqrt(sum(c_weights^2) / sample_size)
        )

        Greater <- 1 - pt(
          q = current_psi_limit / sqrt(sum(c_weights^2) / sample_size),
          df = J * sample_size - J - 1, ncp = psi_standardized / sqrt(sum(c_weights^2) / sample_size)
        )

        expected_widths_too_large <- Less + Greater
        return((expected_widths_too_large - (1 - assurance))^2)
      }
      optimize_result <- optimize(f = determine_limit, interval = c(limit_1_sided, limit_2_sided), psi_standardized = psi_standardized, assurance = assurance)

      n <- ss_aipe_sc_ancova(psi_standardized = optimize_result$minimum,
                             conf_level = NULL,
                             alpha_lower = alpha_lower, alpha_upper = alpha_upper,
                             width = width, c_weights = c_weights,
                             assurance = NULL, ...)[1,2]

      return(.as_dmar_tbl(data.frame(term = 'necessary_n_per_group', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }
  }
  ####################################################################################
  ####################################################################################

  if (divisor == "s_anova") {
    if (abs(sum(c_weights)) > 1e-8) stop("The sum of the contrast weights must be zero")
    if (sum(c_weights[c_weights > 0]) > 1) stop("Please use fractions to specify the contrast weights")

    if (!is.null(psi) || !is.null(sigma_anova) || !is.null(sigma_ancova)) {
      if (is.null(psi) || is.null(sigma_anova) || is.null(sigma_ancova)) stop("'psi', 'sigma_anova', and 'sigma_ancova' must be all specified.")
    }

    if (!is.null(psi_standardized) || !is.null(ratio)) {
      if (is.null(psi_standardized) || is.null(ratio)) stop("Both 'psi_standardized' and 'ratio' (or 'rho') must be specified")
    }

    J <- length(c_weights)

    if (is.null(ratio)) {
      ratio <- sigma_ancova / sigma_anova
      psi_standardized <- psi / sigma_anova
    }

    if (is.null(assurance)) { # Initial starting value for n using the z distribution.
      z_sum <- qnorm(1 - alpha_lower) + qnorm(1 - alpha_upper)
      n_0 <- sum(c_weights^2) * ratio^2 * z_sum^2 / width^2

      # Second starting value for n using the central t distribution.
      t_sum_0 <- qt(1 - alpha_lower, n_0 * J - J - 1) + qt(1 - alpha_upper, n_0 * J - J - 1)
      n <- sum(c_weights^2) * ratio^2 * t_sum_0^2 / width^2

      # measures the discrepancy between the initial and second starting values.
      Difference <- abs(n - n_0)

      while (Difference > .000001) {
        n_p <- n
        t_sum <- qt(1 - alpha_lower, n * J - J - 1) + qt(1 - alpha_upper, n * J - J - 1)
        n <- sum(c_weights^2) * ratio^2 * t_sum^2 / width^2
        Difference <- abs(n - n_p)
      }
      n <- ceiling(n)

      # Initial estimate of noncentral value.
      lambda_0 <- psi_standardized / (ratio * sqrt(sum(c_weights^2) / n))

      lambda_limits_0 <- conf_limits_nct(ncp = lambda_0, df = n * J - J - 1,
                                          conf_level = NULL,
                                          alpha_lower = alpha_lower,
                                          alpha_upper = alpha_upper)
      psi_limit_upper_0 <- lambda_limits_0[which(lambda_limits_0$term == 'upper_limit'),2] * ratio * sqrt(sum(c_weights^2) / n)
      psi_limit_lower_0 <- lambda_limits_0[which(lambda_limits_0$term == 'lower_limit'),2] * ratio * sqrt(sum(c_weights^2) / n)
      # Initial full-width for confidence interval.
      Diff_width_Full <- abs(psi_limit_upper_0 - psi_limit_lower_0) - width

      while (Diff_width_Full > 0) {
        n <- n + 1
        lambda <- psi_standardized / (ratio * sqrt(sum(c_weights^2) / n))
        lambda_limits <- conf_limits_nct(ncp = lambda, df = n * J - J - 1,
                                          conf_level = NULL,
                                          alpha_lower = alpha_lower,
                                          alpha_upper = alpha_upper)
        psi_limit_upper <- lambda_limits[which(lambda_limits$term == 'upper_limit'),2] * ratio * sqrt(sum(c_weights^2) / n)
        psi_limit_lower <- lambda_limits[which(lambda_limits$term == 'lower_limit'),2] * ratio * sqrt(sum(c_weights^2) / n)
        Current_width <- abs(psi_limit_upper - psi_limit_lower)
        Diff_width_Full <- Current_width - width
      }
      return(.as_dmar_tbl(data.frame(term = 'necessary_n_per_group', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }

    if (!is.null(assurance)) {
      psi_standardized <- abs(psi_standardized)

      if ((assurance <= 0) || (assurance >= 1)) stop("The 'assurance' must either be NULL or some value greater than zero and less than unity.", call. = FALSE)
      if (assurance <= .50) stop("The 'assurance' should be larger than 0.5 (but less than 1).", call. = FALSE)

      n0 <- ss_aipe_sc_ancova(psi_standardized = psi_standardized, ratio = ratio, divisor = "s_anova",
                              conf_level = NULL,
                              alpha_lower = alpha_lower, alpha_upper = alpha_upper,
                              width = width, c_weights = c_weights,
                              assurance = NULL, ...)[1,2]

      lambda_2_sided <- conf_limits_nct(ncp = psi_standardized / (ratio * sqrt(sum(c_weights^2) / n0)), df = n0 * J - J - 1, conf_level = NULL, alpha_upper = (1 - assurance) / 2, alpha_lower = (1 - assurance) / 2)

      limit_2_sided <- sqrt(sum(c_weights^2) / n0) * ratio * lambda_2_sided[which(lambda_2_sided$term == 'upper_limit'),2]

      lambda_1_sided <- conf_limits_nct(ncp = psi_standardized / (ratio * sqrt(sum(c_weights^2) / n0)), df = n0 * J - J - 1, conf_level = NULL, alpha_upper = 1 - assurance, alpha_lower = 0)

      limit_1_sided <- sqrt(sum(c_weights^2) / n0) * ratio * lambda_1_sided[which(lambda_1_sided$term == 'upper_limit'),2]

      determine_limit <- function(current_psi_limit = current_psi_limit, sample_size = n0, psi_standardized = psi_standardized, assurance = assurance) {
        Less <- pt(q = -current_psi_limit / (ratio * sqrt(sum(c_weights^2) / sample_size)), df = J * sample_size - J - 1, ncp = psi_standardized / (ratio * sqrt(sum(c_weights^2) / sample_size)))

        Greater <- 1 - pt(q = current_psi_limit / (ratio * sqrt(sum(c_weights^2) / sample_size)), df = J * sample_size - J - 1, ncp = psi_standardized / (ratio * sqrt(sum(c_weights^2) / sample_size)))

        expected_widths_too_large <- Less + Greater
        return((expected_widths_too_large - (1 - assurance))^2)
      }
      optimize_result <- optimize(f = determine_limit, interval = c(limit_1_sided, limit_2_sided), psi_standardized = psi_standardized, assurance = assurance)

      n <- ss_aipe_sc_ancova(psi_standardized = optimize_result$minimum, ratio = ratio, divisor = "s_anova",
                             conf_level = NULL,
                             alpha_lower = alpha_lower, alpha_upper = alpha_upper,
                             width = width, c_weights = c_weights,
                             assurance = NULL, ...)[1,2]

      return(.as_dmar_tbl(data.frame(term = 'necessary_n_per_group', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
    }
  }
}
