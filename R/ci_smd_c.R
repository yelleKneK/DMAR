#' Confidence Limits for the Standardized Mean Difference Using the Control Group Standard Deviation as the Divisor
#'
#' Computes the exact noncentral \emph{t}-based confidence limits for the
#' standardized mean difference that uses the control group standard deviation
#' as the divisor (Glass's \emph{g}). Standardizing by the control group alone
#' keeps the scale of the effect anchored in the untreated population, which
#' matters when the treatment may alter variability as well as the mean.
#'
#' @param ncp The estimated noncentrality parameter, this is generally the observed \emph{t}-statistic from comparing the control and experimental group (assuming homogeneity of variance)
#' @param smd_c The standardized mean difference (using the control group standard deviation in the denominator)
#' @param n_C The sample size for the control group
#' @param n_E The sample size for experimental group
#' @param conf_level The confidence level (1-Type I error rate). A level above .999999998, whose tail areas are
#'   below 1e-9, is refused
#' @param alpha_lower The Type I error rate for the lower tail
#' @param alpha_upper The Type I error rate for the upper tail.
#'   The table records the coverage, \code{1 - alpha_lower - alpha_upper},
#'   as its level and, when the two tail areas differ, the areas
#'   themselves, which the line beneath the printed table,
#'   \code{\link{results_sentence}}, and the \code{alpha_lower} and
#'   \code{alpha_upper} columns of \code{tidy()} state; with a tail area
#'   of 0 the interval is a one-sided confidence bound. A nonzero tail
#'   area below 1e-9 is refused (see "Small tail areas" in
#'   \code{\link{ci_nc_t}}).
#' @param tol The tolerance of the iterative method for determining the critical values
#' @param \dots Potentially include parameter for inner functions
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term} and \code{value}. The
#' \code{term} values are \code{"ci_lower"} (the lower bound of the
#' confidence interval), \code{"smd_c"} (the standardized mean difference
#' standardized by the control group standard deviation), and
#' \code{"ci_upper"} (the upper bound).
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#'   calculation of confidence intervals that are based on central and
#'   noncentral distributions. \emph{Educational and Psychological
#'   Measurement, 61}(4), 532--574. \doi{10.1177/0013164401614002}
#'
#' Glass, G. V. (1976). Primary, secondary, and meta-analysis of research. \emph{Educational Researcher, 5}, 3--8.
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 6}(2), 107--128.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
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
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical models. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @section Warning:
#' This function uses \code{ci_nc_t}, which has as one of its arguments \code{tol} (and can be modified with \code{tol} of the present function).
#' If the present function fails to converge (i.e., if it runs but does not report a solution), it is likely that the \code{tol} value is too restrictive and should be increased by a factor of 10, but probably by no more than 100.
#' Running the function \code{ci_nc_t} directly will report the actual probability values of the limits found. This should be done if any modification to \code{tol} is necessary in order to ensure acceptable confidence limits for the noncentral \emph{t} parameter have been achieved.
#'
#' @seealso
#' \code{\link{smd_c}}, \code{\link{smd}}, \code{\link{ci_smd}}, \code{\link{ci_nc_t}}
#'
#' @examples
#' ci_smd_c(smd_c = .5, n_C = 100, n_E = 100, conf_level = .95)
#'
#' @concept Glass's delta
#'
#' @keywords models htest
#'
#' @family confidence intervals for effect sizes
#'
#' @export

ci_smd_c <- function(ncp = NULL, smd_c = NULL, n_C = NULL, n_E = NULL, conf_level = .95, alpha_lower = NULL, alpha_upper = NULL, tol = 1e-9, ...) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  .dmar_check_numbers("ci_smd_c", c("n_C", "n_E", "tol"))
  if (is.null(ncp) && is.null(smd_c)) .dmar_refuse_argument(
    "You must specify either the estimated noncentral parameter 'ncp' (generally the observed t-statistic) or the standardized mean difference 'smd_c' (as might be obtained from the 'smd_c()' function).",
    arg = c("ncp", "smd_c"))
  if (length(ncp) == 1 && length(smd_c) == 1) .dmar_refuse_argument(
    "You only need to specify either 'ncp' or 'smd_c', not both.",
    arg = c("ncp", "smd_c"))
  if (is.null(n_C) || is.null(n_E)) .dmar_refuse_argument(
    "You must specify sample size per group in order to determine confidence limits.",
    arg = NULL)
  # A size given as an integer (50000L, or a count from length()) is used
  # as a double: R forms the product of two integers as an integer, which
  # overflows to NA once n_C * n_E passes 2^31 - 1. Only integers are
  # converted. as.numeric() of a factor gives its codes, not its labels,
  # so a factor is left as it is: the arithmetic below warns that it is
  # not meaningful for factors, and the call is refused.
  if (is.integer(n_C)) n_C <- as.numeric(n_C)
  if (is.integer(n_E)) n_E <- as.numeric(n_E)
  if (!is.null(conf_level) && conf_level >= 1) .dmar_refuse_argument(
    "'conf_level' must be below 1.", arg = "conf_level")
  if (is.null(conf_level) && sum(alpha_lower, alpha_upper) >= 1) .dmar_refuse_argument(
    "'alpha_lower' and 'alpha_upper' must sum to less than 1.",
    arg = c("alpha_lower", "alpha_upper"))

  # Resolve conf_level vs explicit alpha bounds. Mirrors the ci_nc_t
  # contract: supply one or the other, never both.
  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (alphas_supplied) {
    if (!missing(conf_level) && !is.null(conf_level)) {
      .dmar_refuse_argument(
        "Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'; you cannot mix them.",
        arg = c("conf_level", "alpha_lower", "alpha_upper"))
    }
    conf_level <- NULL
  }

  if (is.null(smd_c)) smd_c <- ncp * sqrt((n_C + n_E) / (n_C * n_E))

  df <- n_C - 1

  if (length(ncp) == 1) {
    Limits <- ci_nc_t(ncp, df, conf_level = conf_level, alpha_lower = alpha_lower, alpha_upper = alpha_upper, tol = tol, ...)
    Limits_L <- Limits[which(Limits$term == "ci_lower"), 2]
    Limits_U <- Limits[which(Limits$term == "ci_upper"), 2]
    Lower_Conf_Limit <- Limits_L * sqrt((n_C + n_E) / (n_C * n_E))
    Upper_Conf_Limit <- Limits_U * sqrt((n_C + n_E) / (n_C * n_E))
    term <- c("ci_lower", "smd_c", "ci_upper")
    value <- c(Lower_Conf_Limit, smd_c, Upper_Conf_Limit)
    out <- data.frame(term, value); attr(out, "conf_level") <- conf_level %||% (1 - alpha_lower - alpha_upper); attr(out, "tail_areas") <- .tail_areas(alpha_lower, alpha_upper); class(out) <- c("dmar_ci_long", "dmar_tbl", "data.frame"); return(out)
  }


  if (length(smd_c) == 1) {
    ncp <- smd_c * sqrt((n_C * n_E) / (n_C + n_E))
    Limits <- ci_nc_t(ncp, df, conf_level = conf_level, alpha_lower = alpha_lower, alpha_upper = alpha_upper, tol = tol, ...)
    Limits_L <- Limits[which(Limits$term == "ci_lower"), 2]
    Limits_U <- Limits[which(Limits$term == "ci_upper"), 2]
    Lower_Conf_Limit <- Limits_L * sqrt((n_C + n_E) / (n_C * n_E))
    Upper_Conf_Limit <- Limits_U * sqrt((n_C + n_E) / (n_C * n_E))
    term <- c("ci_lower", "smd_c", "ci_upper")
    value <- c(Lower_Conf_Limit, smd_c, Upper_Conf_Limit)
    out <- data.frame(term, value); attr(out, "conf_level") <- conf_level %||% (1 - alpha_lower - alpha_upper); attr(out, "tail_areas") <- .tail_areas(alpha_lower, alpha_upper); class(out) <- c("dmar_ci_long", "dmar_tbl", "data.frame"); return(out)
  }
}
