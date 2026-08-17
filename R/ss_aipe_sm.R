#' Sample Size Planning for Accuracy in Parameter Estimation (AIPE) of the Standardized Mean
#'
#' @description
#' Plans the sample size needed for a sufficiently narrow confidence interval
#' for the population standardized mean, the mean divided by the standard
#' deviation, from the accuracy in parameter estimation (AIPE) perspective.
#'
#' @param sm The population standardized mean
#' @param width The desired full width of the obtained confidence interval
#' @param conf_level The desired confidence interval coverage, (i.e., 1 - Type I error rate)
#' @param assurance Parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be \code{NULL} or between zero and unity)
#' @param \dots Allows one to potentially include parameter values for inner functions
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
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{conf_limits_nct}}, \code{\link{ci_sm}}
#'
#' @examples
#' # Suppose the population mean is believed to be 20, and the population
#' # standard deviation is believed to be 2; thus the population standardized
#' # mean is believed to be 10. To determine the necessary sample size for a
#' # study so that the full width of the 95 percent confidence interval
#' # obtained in the study will be, with 90% assurance, no wider than 2.5,
#' # the function should be specified as follows.
#'
#' ss_aipe_sm(sm = 10, width = 2.5, conf_level = .95, assurance = .90)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export


ss_aipe_sm <- function(sm, width, conf_level = .95, assurance = NULL, ...) {
  alpha <- 1 - conf_level

  # Because the properties of the distribution does not change for negative or positive values, other than sign changes,
  # this is added to simplify later code.
  sm <- abs(sm)
  # Thanks to Guy Prochilo <University of Melbourne> for asking a question about negative standardized means that led to a fix in this function.
  # Namely, the addition of the abs() here, where the lack of this caused problems later for assurance parameters.

  # The iterative search calls conf_limits_nct() repeatedly with intermediate
  # trial values of N. When the population standardized mean is large, the
  # implied noncentrality parameter routinely exceeds the 37.62 magnitude at
  # which pt()/qt() lose accuracy, and conf_limits_nct() signals this with a
  # warning. Surfacing that warning once per iteration produces dozens of
  # identical messages (see ?conf_limits_nct). We muffle the per-iteration
  # warnings with withCallingHandlers, count them, and emit a single summary
  # warning at the end (via on.exit so it fires on any return path). The
  # assurance branch calls ss_aipe_sm() recursively; a nested call detects the
  # active handler through .ncp_env and defers to the outer call so exactly one
  # summary warning is emitted per top-level invocation.
  .ncp_env <- getOption("DMAR.ss_aipe_sm_ncp_env", NULL)
  .top_level <- is.null(.ncp_env)
  if (.top_level) {
    .ncp_env <- new.env(parent = emptyenv())
    .ncp_env$count <- 0L
    old_opt <- options(DMAR.ss_aipe_sm_ncp_env = .ncp_env)
    on.exit(options(old_opt), add = TRUE)
    on.exit({
      if (.ncp_env$count > 0L) {
        warning(sprintf(
          "During the iterative sample size search, conf_limits_nct() reported a noncentrality parameter exceeding 37.62 in magnitude in %d intermediate evaluations, the limit at which R's pt()/qt() can return accurate noncentral t probabilities. The returned sample size may be affected; see ?conf_limits_nct.",
          .ncp_env$count
        ), call. = FALSE)
      }
    }, add = TRUE)
  }

  withCallingHandlers({

  if (is.null(assurance)) { # Initial starting value for n using the z distribution.
    n_0 <- (qnorm(1 - alpha / 2) / (width / 2))^2

    n <- ceiling(n_0)

    # To ensure that the initial n is not too big.
    n <- max(4, n - 5)

    # Initial estimate of noncentral value.
    # This is literally the theoretical t-value given delta and the initial estimate of sample size.
    lambda_0 <- sm * sqrt(n)

    # Initial confidence limits.
    lambda_limits_0 <- conf_limits_nct(ncp = lambda_0, df = n - 1, conf_level = 1 - alpha)
    sm_limit_upper_0 <- lambda_limits_0[which(lambda_limits_0$term == 'upper_limit'),2] / sqrt(n)
    sm_limit_lower_0 <- lambda_limits_0[which(lambda_limits_0$term == 'lower_limit'),2] / sqrt(n)

    # Initial full-width for confidence interval.
    Diff_width_Full <- abs(sm_limit_upper_0 - sm_limit_lower_0) - width

    while (Diff_width_Full > 0) {
      n <- n + 1
      lambda <- sm * sqrt(n)
      lambda_limits <- conf_limits_nct(ncp = lambda, df = n - 1, conf_level = 1 - alpha)
      sm_limit_upper <- lambda_limits[which(lambda_limits$term == 'upper_limit'),2] / sqrt(n)
      sm_limit_lower <- lambda_limits[which(lambda_limits$term == 'lower_limit'),2] / sqrt(n)
      Current_width <- abs(sm_limit_upper - sm_limit_lower)
      Diff_width_Full <- Current_width - width
    }

    return(.as_dmar_tbl(data.frame(term = 'necessary_N', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }

  if (!is.null(assurance)) {
    if ((assurance <= 0) || (assurance >= 1)) stop("The 'assurance' must either be NULL or some value greater than zero and less than unity.", call. = FALSE)
    if (assurance <= .50) stop("The 'assurance' should be larger than 0.5 (but less than 1).", call. = FALSE)

    n0 <- ss_aipe_sm(sm = sm, conf_level = conf_level, width = width, assurance = NULL, ...)[1,2]

    Lim_2 <- conf_limits_nct(
      ncp = sm * sqrt(n0), df = n0 - 1, conf_level = NULL,
      alpha_upper = (1 - assurance) / 2, alpha_lower = (1 - assurance) / 2
    )
    limit_2_sided <- (1 / sqrt(n0)) * Lim_2[which(Lim_2$term == 'upper_limit'),2]

    Lim_1 <- conf_limits_nct(
      ncp = sm * sqrt(n0), df = n0 - 1, conf_level = NULL,
      alpha_upper = 1 - assurance, alpha_lower = 0
    )
    limit_1_sided <- (1 / sqrt(n0)) * Lim_1[which(Lim_1$term == 'upper_limit'),2]

    determine_limit <- function(current_sm_limit = current_sm_limit, samp_size = n0, sm = sm,
                                    assurance = assurance) {
      Less <- pt(q = -current_sm_limit * sqrt(samp_size), df = samp_size - 1, ncp = sm * sqrt(samp_size))
      Greater <- 1 - pt(q = current_sm_limit * sqrt(samp_size), df = samp_size - 1, ncp = sm * sqrt(samp_size))
      expected_widths_too_large <- Less + Greater
      return((expected_widths_too_large - (1 - assurance))^2)
    }
    optimize_result <- optimize(
      f = determine_limit, interval = c(limit_1_sided, limit_2_sided),
      sm = sm, assurance = assurance
    )

    n <- ss_aipe_sm(sm = optimize_result$minimum, conf_level = 1 - alpha, width = width, assurance = NULL, ...)[1,2]
    return(.as_dmar_tbl(data.frame(term = 'necessary_N', value = n), conf_level = conf_level, subclass = "dmar_ss_aipe"))
  }

  }, warning = function(w) {
    if (grepl("noncentrality parameter exceeds 37.62", conditionMessage(w))) {
      .ncp_env$count <- .ncp_env$count + 1L
      invokeRestart("muffleWarning")
    }
  })
}
