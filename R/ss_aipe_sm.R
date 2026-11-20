#' Sample Size Planning for Accuracy in Parameter Estimation (AIPE) of the Standardized Mean
#'
#' @description
#' Plans the sample size needed for a sufficiently narrow confidence interval
#' for the population standardized mean, the mean divided by the standard
#' deviation, from the accuracy in parameter estimation (AIPE) perspective.
#'
#' @param true_sm The planning value of the population standardized mean (the population mean divided by the population standard deviation)
#' @param width The desired full width of the obtained confidence interval
#' @param conf_level The desired confidence interval coverage, (i.e., 1 - Type I error rate). A level above
#'   .999999998, whose tail areas are below 1e-9, is refused (see "Small tail areas" in \code{\link{ci_nc_t}})
#' @param assurance Optional. The probability that the confidence interval
#'   of the planned study is no wider than \code{width}: \code{NULL} (the
#'   default), which plans for the expected width, or a single number
#'   strictly between 0 and 1. Any such value may be used. Values such as
#'   \code{.80}, \code{.90}, and \code{.99} are the usual ones, and a value
#'   at or below \code{.50} is accepted as well (see \emph{Details}).
#'   Anything else is refused with a condition of class
#'   \code{dmar_assurance_invalid}. An assurance above .999999998 is
#'   refused with a condition of class \code{dmar_argument_invalid}: the
#'   assurance step finds noncentral \emph{t} limits with tail areas of
#'   \code{(1 - assurance) / 2}, which must be at least 1e-9.
#'
#' @details
#' \strong{The plan for the expected width.} Without \code{assurance} the
#' function finds the smallest sample size at which the noncentral
#' \emph{t} confidence interval for the standardized mean
#' (\code{\link{ci_sm}}), computed at \code{true_sm}, is no wider than
#' \code{width}. The interval of the study itself is computed at the
#' sample standardized mean, which varies from sample to sample, and the
#' interval is wider the larger that value is in magnitude. A plan for the
#' expected width therefore leaves a probability near one half that the
#' interval obtained is wider than \code{width}.
#'
#' \strong{The assurance step.} With \code{assurance} the function follows
#' the one-sample analog of the procedure Kelley and Rausch (2006,
#' Equations 18 to 20 and 25 to 31) give for the standardized mean
#' difference. At the sample size of the plan for the expected width it
#' finds the value that the sample standardized mean exceeds in magnitude
#' with probability \code{1 - assurance}, by a search between the one-sided
#' and the two-sided confidence limits their Equations 25 and 28 name as
#' its bounds, and then plans for the expected width with that value in
#' place of \code{true_sm}. That value is defined for every probability
#' strictly between 0 and 1, so every such \code{assurance} is accepted. A
#' value above about one half gives a larger sample size than the plan for
#' the expected width, and a value below it gives a smaller one, which is
#' the sample size at which the interval is no wider than \code{width}
#' with only that probability: for \code{true_sm = 1}, \code{width = 0.5},
#' and 95 percent confidence the function returns 93 without
#' \code{assurance}, 93 at \code{.50}, 87 at \code{.20}, 100 at
#' \code{.80}, and 105 at \code{.90}. Because the value planned for is found
#' at the first sample size and the returned sample size differs from it,
#' the probability at the returned sample size is near \code{assurance}
#' and not exactly it; in that example it is .519, .228, .812, and .925 at
#' the four values of \code{assurance}, computed from the noncentral
#' \emph{t} distribution.
#'
#' @return
#' A 1-row \code{data.frame} with columns \code{term} and \code{value}:
#' \item{necessary_N}{The necessary total sample size in order to achieve the desired degree of accuracy (i.e., the sufficiently narrow confidence interval)}
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
#' Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the standardized
#' mean difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological
#' Measurement, 65}, 51--69. \doi{10.1177/0013164404264850}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application,
#' and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical models. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_nc_t}}, \code{\link{ci_sm}}
#'
#' @examples
#' # Suppose the population mean is believed to be 20, and the population
#' # standard deviation is believed to be 2; thus the population standardized
#' # mean is believed to be 10. To determine the necessary sample size for a
#' # study so that the full width of the 95 percent confidence interval
#' # obtained in the study will be, with 90% assurance, no wider than 2.5,
#' # the function should be specified as follows.
#'
#' ss_aipe_sm(true_sm = 10, width = 2.5, conf_level = .95, assurance = .90)
#'
#' # The plan for the expected width, and the sample size at which the
#' # interval is no wider than 0.5 with probability .50, for a standardized
#' # mean of 1.
#' ss_aipe_sm(true_sm = 1, width = 0.5)
#' ss_aipe_sm(true_sm = 1, width = 0.5, assurance = .50)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_sm <- function(true_sm, width, conf_level = .95, assurance = NULL) {
  .dmar_check_conf_level(conf_level)
  if (!missing(width)) .dmar_check_width(width)
  .dmar_check_numbers("ss_aipe_sm", "true_sm", required = "true_sm")
  # One check and one classed refusal (dmar_assurance_invalid). Any
  # probability strictly inside (0, 1) is accepted: the value the assurance
  # step plans for (the standardized mean that the sample value exceeds in
  # magnitude with probability 1 - assurance) is defined for each of them.
  .dmar_check_assurance(assurance)
  # The assurance step finds noncentral t limits with tail areas of
  # (1 - assurance) / 2 and 1 - assurance, held to the floor on a tail area
  # (R/dmar_conditions.R).
  .dmar_check_assurance_tail_floor(assurance)
  alpha <- 1 - conf_level

  # Because the properties of the distribution does not change for negative or positive values, other than sign changes,
  # this is added to simplify later code.
  true_sm <- abs(true_sm)
  # Thanks to Guy Prochilo <University of Melbourne> for asking a question about negative standardized means that led to a fix in this function.
  # Namely, the addition of the abs() here, where the lack of this caused problems later for assurance parameters.

  if (is.null(assurance)) { # Initial starting value for n using the z distribution.
    n_0 <- (qnorm(1 - alpha / 2) / (width / 2))^2

    n <- ceiling(n_0)

    # To ensure that the initial n is not too big.
    n <- max(4, n - 5)

    # Initial estimate of noncentral value.
    # This is literally the theoretical t-value given delta and the initial estimate of sample size.
    lambda_0 <- true_sm * sqrt(n)

    # Initial confidence limits.
    lambda_limits_0 <- ci_nc_t(ncp = lambda_0, df = n - 1, conf_level = 1 - alpha)
    sm_limit_upper_0 <- lambda_limits_0[which(lambda_limits_0$term == 'ci_upper'),2] / sqrt(n)
    sm_limit_lower_0 <- lambda_limits_0[which(lambda_limits_0$term == 'ci_lower'),2] / sqrt(n)

    # Initial full-width for confidence interval.
    Diff_width_Full <- abs(sm_limit_upper_0 - sm_limit_lower_0) - width

    while (Diff_width_Full > 0) {
      n <- n + 1
      lambda <- true_sm * sqrt(n)
      lambda_limits <- ci_nc_t(ncp = lambda, df = n - 1, conf_level = 1 - alpha)
      sm_limit_upper <- lambda_limits[which(lambda_limits$term == 'ci_upper'),2] / sqrt(n)
      sm_limit_lower <- lambda_limits[which(lambda_limits$term == 'ci_lower'),2] / sqrt(n)
      Current_width <- abs(sm_limit_upper - sm_limit_lower)
      Diff_width_Full <- Current_width - width
    }

    return(.as_dmar_tbl(data.frame(term = 'necessary_N', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }

  if (!is.null(assurance)) {
    n0 <- ss_aipe_sm(true_sm = true_sm, conf_level = conf_level, width = width, assurance = NULL)[1,2]

    Lim_2 <- ci_nc_t(
      ncp = true_sm * sqrt(n0), df = n0 - 1, conf_level = NULL,
      alpha_upper = (1 - assurance) / 2, alpha_lower = (1 - assurance) / 2
    )
    limit_2_sided <- (1 / sqrt(n0)) * Lim_2[which(Lim_2$term == 'ci_upper'),2]

    Lim_1 <- ci_nc_t(
      ncp = true_sm * sqrt(n0), df = n0 - 1, conf_level = NULL,
      alpha_upper = 1 - assurance, alpha_lower = 0
    )
    limit_1_sided <- (1 / sqrt(n0)) * Lim_1[which(Lim_1$term == 'ci_upper'),2]

    determine_limit <- function(current_sm_limit = current_sm_limit, samp_size = n0, true_sm = true_sm,
                                    assurance = assurance) {
      # The noncentral t probabilities are exact (.nc_t_prob() uses pt()
      # only where it is accurate; see R/nc_t_internals.R).
      Less <- .nc_t_prob(-current_sm_limit * sqrt(samp_size), samp_size - 1, true_sm * sqrt(samp_size))
      Greater <- .nc_t_prob(current_sm_limit * sqrt(samp_size), samp_size - 1, true_sm * sqrt(samp_size), lower_tail = FALSE)
      expected_widths_too_large <- Less + Greater
      return((expected_widths_too_large - (1 - assurance))^2)
    }
    optimize_result <- optimize(
      f = determine_limit, interval = c(limit_1_sided, limit_2_sided),
      true_sm = true_sm, assurance = assurance
    )

    n <- ss_aipe_sm(true_sm = optimize_result$minimum, conf_level = 1 - alpha, width = width, assurance = NULL)[1,2]
    return(.as_dmar_tbl(data.frame(term = 'necessary_N', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }
}
