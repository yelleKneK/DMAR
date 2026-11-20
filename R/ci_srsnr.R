#' Confidence Interval for the Square Root of the Signal-to-Noise Ratio
#'
#' Computes the exact confidence interval for the square root of the
#' signal-to-noise ratio, the standard deviation of the group means
#' relative to the error standard deviation. On this root scale the
#' quantity is an effect size in standard deviation units, the
#' multi-group analog of the standardized mean difference.
#'
#' @param F_value Observed \emph{F}-value from the analysis of variance.
#'   Use this argument when re-analyzing existing data.
#' @param df_1 Numerator degrees of freedom
#' @param df_2 Denominator degrees of freedom
#' @param N Sample size
#' @param means Numeric vector of group means, one per group and at least
#'   two, given with \code{var_within} and \code{n_per_group} in place of
#'   \code{F_value}: the observed sample means of a one-way design, from
#'   which the function computes the observed \emph{F}, or population means
#'   for a hypothetical interval (see Details).
#' @param var_within The within-group variance that goes with \code{means},
#'   a single positive number: the pooled within-group mean square (the
#'   mean square error of the one-way analysis of variance), or a
#'   population within-group variance for a hypothetical interval.
#' @param n_per_group The group sizes that go with \code{means}: a single
#'   size for equal groups, or one size per group, each a whole number of
#'   at least 2.
#' @param conf_level Confidence interval coverage (i.e., 1 - Type I error rate); default is .95
#' @param alpha_lower Type I error for the lower confidence limit
#' @param alpha_upper Type I error for the upper confidence limit.
#'   The table records the coverage, \code{1 - alpha_lower - alpha_upper},
#'   as its level and, when the two tail areas differ, the areas
#'   themselves, which the line beneath the printed table,
#'   \code{\link{results_sentence}}, and the \code{alpha_lower} and
#'   \code{alpha_upper} columns of \code{tidy()} state; with a tail area
#'   of 0 the interval is a one-sided confidence bound.
#' @param \dots Allows one to potentially include parameter values for inner functions
#'
#' @details
#' The confidence level is given in one of two ways: as the coverage,
#' \code{conf_level}, or as the two tail areas, \code{alpha_lower} and
#' \code{alpha_upper}.
#'
#' The square root of the signal-to-noise ratio is the standard deviation
#' of the group means, each weighted by its share of the total sample size
#' \eqn{n_j / N}, divided by the within-group standard deviation; it is
#' \eqn{\sqrt{\lambda / N}} for the noncentrality parameter \eqn{\lambda}
#' of the \emph{F} test, and it is Cohen's \emph{f} with those weights (see
#' \code{\link{cohen_f}}). This function uses the confidence interval
#' transformation principle (Steiger, 2004) to transform the confidence
#' limits for the noncentrality parameter to the confidence limits for the
#' square root of the signal-to-noise ratio. The confidence interval for the
#' noncentrality parameter itself is given by \code{\link{ci_nc_F}}.
#'
#' The interval has two routes. The first takes the observed \emph{F}, its
#' degrees of freedom, and \emph{N}. The second takes the group means
#' \eqn{\bar{y}_j}{ybar_j}, the pooled within-group mean square
#' \eqn{MS_{within}}{MS_within}, and the group sizes \eqn{n_j} of a one-way design
#' with \emph{k} groups, and computes
#' \deqn{F = \frac{\sum_j n_j (\bar{y}_j - \bar{y})^2 / (k - 1)}{MS_{within}},
#'   \qquad \bar{y} = \frac{1}{N} \sum_j n_j \bar{y}_j,}{F = [sum_j n_j
#'   (ybar_j - ybar)^2 / (k - 1)] / MS_within, with ybar = sum_j n_j ybar_j / N,}
#' with \eqn{df_1 = k - 1}, \eqn{df_2 = N - k}, and
#' \eqn{N = \sum_j n_j}. With the observed sample means and mean square,
#' this is the observed \emph{F}, and the interval is the one the first
#' route gives for the same data. Values of \code{df_1}, \code{df_2}, or
#' \code{N} given with \code{means} must agree with these. Group means that
#' are all equal give \eqn{F = 0}, which is refused, as it is on the first
#' route: at \eqn{F = 0} the upper confidence limit is undefined.
#'
#' Population means and a population within-group variance can be given in
#' the same places. The result is then a hypothetical interval: the interval
#' a study with these group sizes would report if its sample means and its
#' pooled within-group mean square equaled the values given. It is not an
#' interval from observed data, and it is not the expected interval of such
#' a study: the limits of the interval vary from sample to sample, and the
#' interval computed at the population values need not equal their average.
#'
#' @return
#' A 2-row \code{data.frame} with columns \code{term} and \code{value}. The
#' \code{term} values are \code{"ci_lower"} and \code{"ci_upper"},
#' giving the square roots of the corresponding signal-to-noise-ratio
#' confidence limits.
#'
#' @references
#' Bargmann, R. E. (1970). Interpretation and use of a generalized
#'   discriminant function. In R. C. Bose, I. M. Chakravarti, P. C.
#'   Mahalanobis, C. R. Rao, & K. J. C. Smith (Eds.), \emph{Essays in
#'   probability and statistics} (pp. 35--60). University of North
#'   Carolina Press.
#'
#' Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
#'   \emph{Educational and Psychological Measurement, 40}(3), 659--670.
#'   \doi{10.1177/001316448004000309}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' Venables, W. (1975). Calculation of confidence intervals for
#'   noncentrality parameters. \emph{Journal of the Royal Statistical
#'   Society, Series B, 37}(3), 406--412.
#'   \doi{10.1111/j.2517-6161.1975.tb01554.x}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_snr}}, \code{\link{ci_nc_F}}
#'
#' @examples
#' ## Bargmann (1970) gave an example of a one-way ANOVA with an observed
#' ## F value of 11.2213 on 4 and 50 degrees of freedom, which Venables (1975),
#' ## Fleishman (1980), and Steiger (2004) used in turn. Steiger describes it
#' ## as five groups of 11 and reports F as 11.221, the value used here.
#' ## The exact confidence interval for the square root of the signal-to-noise ratio of that example:
#'
#' ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
#'
#' ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, conf_level = .90)
#'
#' ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, alpha_lower = .02, alpha_upper = .03)
#'
#' # The interval from the summaries of an observed one-way design: the
#' # group means, the pooled within-group mean square, and the group sizes
#' # of the three arms of drinks_trial. The F route gives the same interval
#' # for the same data.
#' fit <- stats::lm(log_drinks ~ treatment, data = drinks_trial)
#' ci_srsnr(means = tapply(drinks_trial$log_drinks, drinks_trial$treatment, mean),
#'          var_within = sum(stats::residuals(fit)^2) / stats::df.residual(fit),
#'          n_per_group = as.vector(table(drinks_trial$treatment)))
#' ci_srsnr(F_value = stats::anova(fit)[1, "F value"], df_1 = 2, df_2 = 85, N = 88)
#'
#' # A hypothetical interval: with population means of 94, 91, 92, and 83, a
#' # within-group variance of 67.375, and 20 per group, the interval a study
#' # would report if its summaries equaled these values. It is not the
#' # expected interval of that study.
#' ci_srsnr(means = c(94, 91, 92, 83), var_within = 67.375, n_per_group = 20)
#'
#' @keywords design
#'
#' @family confidence intervals for effect sizes
#'
#' @export

ci_srsnr <- function(F_value = NULL, df_1 = NULL, df_2 = NULL, N = NULL,
                     means = NULL, var_within = NULL, n_per_group = NULL,
                     conf_level = .95,
                     alpha_lower = NULL, alpha_upper = NULL, ...) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  .dmar_check_numbers("ci_srsnr", c("df_1", "df_2", "N", "var_within",
                                    "n_per_group"))

  if (is.null(means) && (!is.null(var_within) || !is.null(n_per_group))) {
    .dmar_refuse_argument(c(
      "'var_within' and 'n_per_group' are used only with 'means', which was ",
      "not given. Give the group means as well, to compute F from the ",
      "summaries, or give 'F_value', 'df_1', 'df_2', and 'N' without them."),
      arg = c("var_within", "n_per_group", "means"))
  }

  # Summary route: F, df_1, df_2, and N from the group means, the pooled
  # within-group mean square, and the group sizes.
  if (!is.null(means)) {
    if (is.null(var_within) || is.null(n_per_group)) {
      .dmar_refuse_argument(c(
        "When 'means' is given, 'var_within' and 'n_per_group' must be ",
        "given with it: F is computed from the group means, the pooled ",
        "within-group mean square, and the group sizes."),
        arg = c("means", "var_within", "n_per_group"))
    }
    if (!is.null(F_value)) {
      .dmar_refuse_argument(c(
        "Give either 'F_value' or 'means' (with 'var_within' and ",
        "'n_per_group'), not both: each determines F, and the two could ",
        "disagree."),
        arg = c("F_value", "means"))
    }
    if (!is.numeric(means) || length(means) < 2L || !all(is.finite(means))) {
      .dmar_refuse_argument(c(
        "'means' must be a numeric vector of at least two finite group ",
        "means; ", .dmar_given(means), "."), arg = "means", value = means)
    }
    means <- as.vector(means)
    if (!is.numeric(var_within) || length(var_within) != 1L ||
        !is.finite(var_within) || var_within <= 0) {
      .dmar_refuse_argument(c(
        "'var_within' must be a single positive number, the pooled ",
        "within-group mean square; ", .dmar_given(var_within), "."),
        arg = "var_within", value = var_within)
    }
    if (length(n_per_group) == 1L) n_per_group <- rep(n_per_group, length(means))
    if (length(n_per_group) != length(means)) {
      .dmar_refuse_argument(
        "'n_per_group' must be a scalar or a vector the same length as 'means'.",
        arg = c("n_per_group", "means"))
    }
    n_per_group <- as.vector(n_per_group)
    # Group sizes set df_2 and N, so they are counts: a fractional size
    # would give a fractional df_2 that no observed design has.
    if (!all(is.finite(n_per_group)) || any(n_per_group <= 1) ||
        any(n_per_group != round(n_per_group))) .dmar_refuse_argument(c(
      "Each value of 'n_per_group' must be a whole number of at least 2, ",
      "the size of a group; the group sizes were ",
      paste(as.character(n_per_group), collapse = ", "), "."),
      arg = "n_per_group", value = n_per_group)
    J          <- length(means)
    N_summary  <- sum(n_per_group)
    implied    <- c(df_1 = J - 1L, df_2 = N_summary - J, N = N_summary)
    supplied   <- list(df_1 = df_1, df_2 = df_2, N = N)
    disagree   <- vapply(names(supplied), function(nm) {
      x <- supplied[[nm]]
      !is.null(x) && !(length(x) == 1L && abs(x - implied[[nm]]) <=
                         sqrt(.Machine$double.eps) * max(1, abs(implied[[nm]])))
    }, logical(1))
    if (any(disagree)) {
      bad <- names(supplied)[disagree]
      .dmar_refuse_argument(c(
        paste0("'", bad, "'", collapse = ", "), " disagree", if (length(bad) == 1L) "s",
        " with the group sizes: with ", J, " groups and n_per_group summing ",
        "to ", format(N_summary), ", df_1 is ", format(implied[["df_1"]]),
        ", df_2 is ", format(implied[["df_2"]]), ", and N is ",
        format(implied[["N"]]), ". Leave them out when 'means' is given, ",
        "or give values that agree."),
        arg = bad)
    }
    N          <- N_summary
    df_1       <- implied[["df_1"]]
    df_2       <- implied[["df_2"]]
    grand_mean <- sum(n_per_group * means) / N
    SS_between <- sum(n_per_group * (means - grand_mean)^2)
    F_value    <- (SS_between / df_1) / var_within
    if (F_value <= 0) {
      .dmar_refuse_argument(c(
        "The group means in 'means' are all equal, so the F computed from ",
        "the summaries is 0. ci_srsnr() needs an F larger than 0, on this ",
        "route as on the 'F_value' route: at F = 0 the upper confidence ",
        "limit is undefined."),
        arg = "means", value = means)
    }
  }

  if (is.null(alpha_lower) && is.null(alpha_upper)) {
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }

  if (!is.null(alpha_lower) && !is.null(alpha_upper)) {
    conf_level <- 1 - alpha_upper - alpha_lower
  }

  if (!is.null(alpha_lower) && is.null(alpha_upper)) .dmar_refuse_argument(
    "This is a problem with the desired confidence level ('alpha_lower' specified but not 'alpha_upper').",
    arg = c("alpha_lower", "alpha_upper"))
  if (is.null(alpha_lower) && !is.null(alpha_upper)) .dmar_refuse_argument(
    "This is a problem with the desired confidence level ('alpha_upper' specified but not 'alpha_lower').",
    arg = c("alpha_upper", "alpha_lower"))

  if (is.null(F_value)) .dmar_refuse_argument(
    "Specify either 'F_value' (with 'df_1', 'df_2', and 'N') or 'means' with 'var_within' and 'n_per_group'.",
    arg = c("F_value", "means", "var_within"))
  if (is.null(df_1) || is.null(df_2) || is.null(N)) .dmar_refuse_argument(
    "You need to specify 'df_1', 'df_2', and 'N'.",
    arg = c("df_1", "df_2", "N"))

  if (alpha_lower > .5 || alpha_lower < 0) .dmar_refuse_argument(
    "'alpha_lower' must be between 0 and .5.",
    arg = "alpha_lower")
  if (conf_level > 1 || conf_level < 0) .dmar_refuse_argument(
    "'conf_level' must be between 0 and 1.",
    arg = "conf_level")
  if (F_value <= 0) .dmar_refuse_argument(
    "'F_value' must be larger than 0.", arg = "F_value")
  if (N <= 0 || N <= df_1 + df_2) .dmar_refuse_argument(
    "'N' must be larger than 'df_1' + 'df_2'.", arg = c("N", "df_1", "df_2"))

  Lims <- .ci_nc_F_for(
    caller = "ci_srsnr",
    quantity = "the square root of the signal-to-noise ratio",
    F_value = F_value, conf_level = NULL, df_1 = df_1,
    df_2 = df_2, alpha_lower = alpha_lower, alpha_upper = alpha_upper, ...
  )

  term <- c("ci_lower", "ci_upper")
  value <- c(sqrt(max(0, Lims[which(Lims$term == "ci_lower"), 2] / N)), sqrt(Lims[which(Lims$term == "ci_upper"), 2] / N))
  length(value) <- length(term)
  return(.as_dmar_tbl(data.frame(term, value), conf_level = conf_level,
                      tail_areas = .tail_areas(alpha_lower, alpha_upper)))
}
