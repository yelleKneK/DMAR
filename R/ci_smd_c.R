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
#' @param conf_level The confidence level (1-Type I error rate)
#' @param alpha_lower The Type I error rate for the lower tail
#' @param alpha_upper The Type I error rate for the upper tail
#' @param tol The tolerance of the iterative method for determining the critical values
#' @param \dots Potentially include parameter for inner functions
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term} and \code{value}. The
#' \code{term} values are \code{"lower_limit"} (the lower bound of the
#' confidence interval), \code{"smd_c"} (the standardized mean difference
#' standardized by the control group standard deviation), and
#' \code{"upper_limit"} (the upper bound).
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
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @section Warning:
#' This function uses \code{conf_limits_nct}, which has as one of its arguments \code{tol} (and can be modified with \code{tol} of the present function).
#' If the present function fails to converge (i.e., if it runs but does not report a solution), it is likely that the \code{tol} value is too restrictive and should be increased by a factor of 10, but probably by no more than 100.
#' Running the function \code{conf_limits_nct} directly will report the actual probability values of the limits found. This should be done if any modification to \code{tol} is necessary in order to ensure acceptable confidence limits for the noncentral \emph{t} parameter have been achieved.
#'
#' @seealso
#' \code{\link{smd_c}}, \code{\link{smd}}, \code{\link{ci_smd}}, \code{\link{conf_limits_nct}}
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
  if (is.null(ncp) && is.null(smd_c)) stop("You must specify either the estimated noncentral parameter 'ncp' (generally the observed t-statistic) or the standardized mean difference 'smd_c' (as might be obtained from the 'smd_c()' function).", call. = FALSE)
  if (length(ncp) == 1 && length(smd_c) == 1) stop("You only need to specify either 'ncp' or 'smd_c', not both.", call. = FALSE)
  if (is.null(n_C) || is.null(n_E)) stop("You must specify sample size per group in order to determine confidence limits.", call. = FALSE)
  if (!is.null(conf_level) && conf_level >= 1) stop("There is a problem with your confidence level.", call. = FALSE)
  if (is.null(conf_level) && sum(alpha_lower, alpha_upper) >= 1) stop("There is a problem with your upper and or lower confidence limits.", call. = FALSE)

  # Resolve conf_level vs explicit alpha bounds. Mirrors the conf_limits_nct
  # contract: supply one or the other, never both.
  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (alphas_supplied) {
    if (!missing(conf_level) && !is.null(conf_level)) {
      stop("Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'; you cannot mix them.", call. = FALSE)
    }
    conf_level <- NULL
  }

  if (is.null(smd_c)) smd_c <- ncp * sqrt((n_C + n_E) / (n_C * n_E))

  df <- n_C - 1

  if (length(ncp) == 1) {
    Limits <- conf_limits_nct(ncp, df, conf_level = conf_level, alpha_lower = alpha_lower, alpha_upper = alpha_upper, tol = tol, ...)
    Limits_L <- Limits[which(Limits$term == "lower_limit"), 2]
    Limits_U <- Limits[which(Limits$term == "upper_limit"), 2]
    Lower_Conf_Limit <- Limits_L * sqrt((n_C + n_E) / (n_C * n_E))
    Upper_Conf_Limit <- Limits_U * sqrt((n_C + n_E) / (n_C * n_E))
    term <- c("lower_limit", "smd_c", "upper_limit")
    value <- c(Lower_Conf_Limit, smd_c, Upper_Conf_Limit)
    out <- data.frame(term, value); attr(out, "conf_level") <- conf_level; class(out) <- c("dmar_ci_long", "dmar_tbl", "data.frame"); return(out)
  }


  if (length(smd_c) == 1) {
    ncp <- smd_c * sqrt((n_C * n_E) / (n_C + n_E))
    Limits <- conf_limits_nct(ncp, df, conf_level = conf_level, alpha_lower = alpha_lower, alpha_upper = alpha_upper, tol = tol, ...)
    Limits_L <- Limits[which(Limits$term == "lower_limit"), 2]
    Limits_U <- Limits[which(Limits$term == "upper_limit"), 2]
    Lower_Conf_Limit <- Limits_L * sqrt((n_C + n_E) / (n_C * n_E))
    Upper_Conf_Limit <- Limits_U * sqrt((n_C + n_E) / (n_C * n_E))
    term <- c("lower_limit", "smd_c", "upper_limit")
    value <- c(Lower_Conf_Limit, smd_c, Upper_Conf_Limit)
    out <- data.frame(term, value); attr(out, "conf_level") <- conf_level; class(out) <- c("dmar_ci_long", "dmar_tbl", "data.frame"); return(out)
  }
}
