#' @rdname ss_aipe_smd
#' @name ss_aipe_smd
#'
#' @title Sample Size Planning for the Standardized Mean Difference (AIPE)
#'
#' @description
#' Determines the per-group sample size needed for a two-independent-groups
#' design so that the (expected) confidence interval for Cohen's
#' \emph{d}, the population standardized mean difference, denoted
#' \eqn{\delta}, is no wider than a user-specified value.  This is the
#' Accuracy in Parameter Estimation (AIPE) framework of Kelley and Rausch
#' (2006), the standardized-mean-difference companion to power-based
#' planning via \code{\link{ss_power_smd}}.  Optionally, supplying
#' \code{assurance} returns the larger sample size needed so that the
#' realized interval will be at or below the target width with that
#' probability rather than just on average.
#'
#' @param delta The supposed value of the population standardized mean
#'   difference \eqn{\delta} the sample size is planned against: a value the
#'   researcher posits, either a minimally important effect or a value believed
#'   to be true in the population, never a sample estimate. Echoed in the
#'   returned table as the \code{supposed_smd} row.
#' @param conf_level Desired confidence level (i.e., \eqn{1-\alpha}, where
#'   \eqn{\alpha} is the Type I error rate). Default \code{0.95}.
#' @param width Desired (full) width of the two-sided confidence interval
#'   on \eqn{\delta}.
#' @param assurance Optional probability with which the realized
#'   confidence interval is to be no wider than \code{width}.  When
#'   \code{NULL} (the default), the planning targets the
#'   \emph{expected} width; when supplied (e.g., 0.80, 0.90, 0.99), the
#'   procedure returns the larger \emph{N} that guarantees the desired
#'   width with that assurance.  Must be \code{NULL} or strictly between
#'   0.50 and 1.
#'
#' @return
#' A \code{data.frame} with columns \code{term} and \code{value}. The first
#' row, \code{necessary_n_per_group}, is the necessary per-group sample size
#' \emph{N}; the remaining rows echo the user-supplied planning inputs
#' \code{supposed_smd} and \code{width} (and \code{assurance} when supplied), so
#' the assumptions the sample size was planned under travel with the result. The
#' \code{supposed_smd} row is the supposed effect the plan is built on: a value
#' the researcher posits, either a minimally important effect or a value
#' believed to be true in the population, never a sample estimate. The
#' confidence level is reported in the printed footer.
#'
#' @references
#' Anderson, S. F., & Kelley, K. (2024). Sample size planning for
#'   replication studies: The devil is in the design. \emph{Psychological
#'   Methods, 29}(5), 844--867. \doi{10.1037/met0000520}
#'
#' Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
#'   planning for more accurate statistical power: A method adjusting
#'   sample effect sizes for publication bias and uncertainty.
#'   \emph{Psychological Science, 28}(11), 1547--1562.
#'   \doi{10.1177/0956797617723724}
#'
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#'   calculation of confidence intervals that are based on central and
#'   noncentral distributions. \emph{Educational and Psychological
#'   Measurement, 61}(4), 532--574. \doi{10.1177/0013164401614002}
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 6}(2), 107--128.
#'
#' Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the standardized mean
#' difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological Measurement, 65}, 51--69.
#'   \doi{10.1177/0013164404264850}
#'
#' Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
#'   obtaining precision: Delineating methods of sample size planning.
#'   \emph{Evaluation and the Health Professions, 26}(3), 258--287.
#'   \doi{10.1177/0163278703255242}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 4 on individual comparisons and
#'   Chapter 3 on one-way ANOVA.)
#'
#' Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
#'   for statistical power and accuracy in parameter estimation.
#'   \emph{Annual Review of Psychology, 59}, 537--563.
#'   \doi{10.1146/annurev.psych.59.103006.093735}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @section Warning: The returned value is the sample size \emph{per group}.
#'
#' @seealso
#' \code{\link{smd}}, \code{\link{smd_c}}, \code{\link{ci_smd}}, \code{\link{ci_smd_c}}, 
#' \code{\link{conf_limits_nct}}, \code{\link[stats:power.t.test]{stats::power.t.test()}}
#'
#' @examples
#' ss_aipe_smd(delta = .5, conf_level = .95, width = .30)
#' ss_aipe_smd(delta = .5, conf_level = .95, width = .30, assurance = .8)
#' ss_aipe_smd(delta = .5, conf_level = .95, width = .30, assurance = .95)
#'
#' @keywords design
#'
#' @seealso \code{\link{design_consequences}} for what a chosen design delivers:
#'   power, the Type S (sign) and Type M (exaggeration) errors of the
#'   significance filter, and the expected confidence interval width.
#'
#' @export

ss_aipe_smd <- function(delta, conf_level = 0.95, width, assurance = NULL) {
  # Validate the planning targets at entry so a boundary input errors clearly
  # rather than sending qnorm()/qt() out of domain and crashing the fixed-point
  # search with a missing value.
  if (missing(delta) || !is.numeric(delta) || length(delta) != 1L || !is.finite(delta))
    stop("'delta' must be a single finite number.", call. = FALSE)
  if (missing(width) || !is.numeric(width) || length(width) != 1L || !is.finite(width) || width <= 0)
    stop("'width' must be a single finite positive number.", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be a single number strictly between 0 and 1.", call. = FALSE)

  alpha <- 1 - conf_level

  # An explicit cap on each fixed-point search so a pathological input cannot
  # loop without terminating.
  max_iter <- 100000L

  ss_aipe_smd_expected_width <- function(delta, conf_level = 1 - alpha, width) {
    alpha <- 1 - conf_level

    # Initial starting value for n using the z distribution.
    n_0 <- 2 * (qnorm(1 - alpha / 2) / (width / 2))^2

    # Second starting value for n using the central t distribution. Floor the
    # error df at 1 so that a very wide target (whose z-based start can be below
    # 1 per group) never passes a nonpositive df to qt() and yields NaN; the
    # final n is floored at 4 below regardless.
    n <- 2 * ((qt(1 - alpha / 2, max(2 * n_0 - 2, 1))) / (width / 2))^2

    # Measures the discrepancy between the initial and second starting values.
    Difference <- abs(n - n_0)

    iter <- 0L
    while (Difference > .000001) {
      iter <- iter + 1L
      if (iter > max_iter)
        stop("The sample size search did not converge within ", max_iter, " iterations; the specified target may be infeasible.", call. = FALSE)
      n_p <- n
      n <- 2 * ((qt(1 - alpha / 2, max(2 * n - 2, 1))) / (width / 2))^2
      Difference <- abs(n - n_p)
    }
    n <- ceiling(n)

    # To ensure that the initial n is not too small. Four per group is the
    # admissible minimum this planner returns, so an impossibly wide target
    # settles here rather than at a degenerate value.
    n <- max(4, n - 5)

    # Initial estimate of noncentral value.
    # This is literally the theoretical t-value given delta and the initial estimate of sample size.
    lambda_0 <- delta * sqrt(n / 2)

    # Initial confidence limits.
    Limits_0 <- ci_smd(ncp = lambda_0, n_1 = n, n_2 = n, conf_level = 1 - alpha)

    # Initial (full-width) for confidence interval.
    diff_width <- abs(Limits_0[Limits_0$term == "upper_limit", "value"] - Limits_0[Limits_0$term == "lower_limit", "value"]) - width

    iter <- 0L
    while (diff_width > 0) {
      iter <- iter + 1L
      if (iter > max_iter)
        stop("The sample size search did not converge within ", max_iter, " iterations; the specified target may be infeasible.", call. = FALSE)
      n <- n + 1
      lambda <- delta * sqrt(n / 2)
      Limits <- ci_smd(ncp = lambda, n_1 = n, n_2 = n, conf_level = 1 - alpha)
      Current_width <- abs(Limits[Limits$term == "upper_limit", "value"] - Limits[Limits$term == "lower_limit", "value"])
      diff_width <- Current_width - width
    }
    return(n)
  }

  n <- ss_aipe_smd_expected_width(delta = delta, conf_level = conf_level, width = width)

  # This part is for the optional assurance (that the CI will be narrow enough).
  if (!is.null(assurance)) {
    if ((assurance <= 0) || (assurance >= 1)) stop("The 'assurance' must either be NULL or some value greater than zero and less than unity.", call. = FALSE)
    if (assurance <= .50) stop("The 'assurance' should be > 0.50 (but less than 1).", call. = FALSE)

    # Start with the original n, if only the expected width is specified.
    n0 <- n

    limit_2_sided <- ci_smd(smd = delta, n_1 = n0, n_2 = n0, conf_level = NULL, alpha_lower = (1 - assurance) / 2, alpha_upper = (1 - assurance) / 2)
    limit_2_sided <- limit_2_sided[limit_2_sided$term == "upper_limit", "value"]

    limit_1_sided <- ci_smd(smd = delta, n_1 = n0, n_2 = n0, conf_level = NULL, alpha_lower = 0, alpha_upper = 1 - assurance)
    limit_1_sided <- limit_1_sided[limit_1_sided$term == "upper_limit", "value"]

    determine_limit <- function(current_delta_limit = current_delta_limit, samp_size = n0, delta = delta, assurance = assurance) {
      Less <- pt(q = convert_delta_lambda(delta = -current_delta_limit, n_1 = samp_size, n_2 = samp_size)[1, 2], df = 2 * samp_size - 2, ncp = convert_delta_lambda(delta = delta, n_1 = samp_size, n_2 = samp_size)[1, 2])
      Greater <- 1 - pt(q = convert_delta_lambda(delta = current_delta_limit, n_1 = samp_size, n_2 = samp_size)[1, 2], df = 2 * samp_size - 2, ncp = convert_delta_lambda(delta = delta, n_1 = samp_size, n_2 = samp_size)[1, 2])
      expected_widths_too_large <- Less + Greater
      return((expected_widths_too_large - (1 - assurance))^2)
    }

    optimize_result <- optimize(f = determine_limit, interval = c(limit_1_sided, limit_2_sided), delta = delta, assurance = assurance)
    # This uses the basic approach, with no assurance (i.e., the expected width) but with a new value of delta (based on the above). This overwrites the first n.
    n <- ss_aipe_smd_expected_width(delta = optimize_result$minimum, conf_level = conf_level, width = width)
  }

  # The necessary per-group sample size comes first, followed by rows that echo
  # the user-supplied planning inputs, so the assumptions the sample size was
  # planned under travel with the result. supposed_smd is the supposed effect
  # the plan is built on, a value the researcher posits (a minimally important
  # effect or a value believed to be true in the population), never a sample
  # estimate. The confidence level is reported in the printed footer; the value
  # column stays numeric.
  terms  <- c("necessary_n_per_group", "supposed_smd", "width")
  values <- c(n, delta, width)
  if (!is.null(assurance)) {
    terms  <- c(terms, "assurance")
    values <- c(values, assurance)
  }
  return(.as_dmar_tbl(data.frame(term = terms, value = values),
                      conf_level = conf_level, subclass = "dmar_ss_aipe"))
}

#' @rdname ss_aipe_smd
#' @export
