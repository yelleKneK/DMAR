#' Confidence Interval for the Standardized Mean Difference (Two Independent Groups)
#'
#' Constructs an exact-coverage confidence interval for the population
#' standardized mean difference \eqn{\delta = (\mu_1 - \mu_2)/\sigma}
#' (Cohen's \emph{d} when expressed as a sample quantity) for two
#' independent groups under bivariate normality with equal variances. The
#' interval is obtained by inverting the noncentral \emph{t} sampling
#' distribution of the rescaled statistic
#' \eqn{t = \hat d \sqrt{n_1 n_2 / (n_1 + n_2)}}, which under the
#' independent groups equal variances model is exactly noncentral
#' \emph{t} with \eqn{n_1 + n_2 - 2} degrees of freedom and noncentrality
#' parameter \eqn{\lambda = \delta \sqrt{n_1 n_2 / (n_1 + n_2)}} (Hedges,
#' 1981). The confidence limits are then rescaled back to the \eqn{\delta}
#' metric. This is the same construction Steiger and Fouladi (1997) and
#' Kelley (2007) describe for noncentral effect size CIs.
#'
#' @param ncp The estimated noncentrality parameter, this is generally the observed \emph{t}-statistic from comparing the two groups and assumes homogeneity of variance
#' @param smd The standardized mean difference (using the pooled standard deviation in the denominator)
#' @param n_1 The sample size for Group 1
#' @param n_2 The sample size for Group 2
#' @param conf_level The confidence level (1-Type I error rate)
#' @param alpha_lower The Type I error rate for the lower tail
#' @param alpha_upper The Type I error rate for the upper tail
#' @param tol The tolerance of the iterative method for determining the critical values
#' @param \dots Allows one to potentially include parameter values for inner functions
#'
#' @details
#' \strong{ncp-input vs. smd-input paths.} The function accepts the effect
#' size in either of two equivalent metrics: the observed
#' \emph{t}-statistic (via \code{ncp}) or the sample standardized mean
#' difference (via \code{smd}). The two paths are mathematically
#' equivalent under the equal variances assumption (since
#' \eqn{t = \hat d \sqrt{n_1 n_2 / (n_1 + n_2)}}); pick whichever is
#' easier to obtain. Supply exactly one. Both paths internally call
#' \code{\link{conf_limits_nct}} to invert the noncentral \emph{t}
#' distribution at the specified two-tailed (or asymmetric, via
#' \code{alpha_lower} / \code{alpha_upper}) confidence level.
#'
#' \strong{Independent vs.\ paired comparison.} \code{ci_smd} assumes two
#' \emph{independent} groups with a common variance. DMAR does not
#' currently provide a confidence interval for the standardized mean
#' difference in a paired or within-subject design, whose sampling
#' distribution depends on the correlation between the paired
#' measurements; applying the independent groups interval to paired data
#' gives the wrong coverage. (\code{\link{ci_smd_c}} is not a paired
#' interval either; it is the interval for Glass's estimator, which
#' standardizes the difference between two independent groups by the
#' control group standard deviation.)
#'
#' \strong{Bias correction (Hedges' g).} \code{ci_smd} reports the CI on
#' \emph{d}; if the bias-corrected \emph{g} is desired, multiply the
#' bounds by the Hedges and Olkin (1985) correction factor
#' \eqn{J(\nu) = 1 - 3/(4 \nu - 1)} (with \eqn{\nu = n_1 + n_2 - 2}).
#' Because \eqn{J(\nu)} is a constant, the rescaling preserves coverage.
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term} and \code{value}. The
#' \code{term} values are \code{"lower_limit"} (the lower bound of the
#' confidence interval on the standardized mean difference), \code{"smd"} (the
#' point estimate), and \code{"upper_limit"} (the upper bound).
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
#' Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 6}(2), 107--128.
#'
#' Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for
#'   meta-analysis}. Academic Press.
#'
#' Kelley, K. (2005). The effects of nonnormal distributions on confidence
#'   intervals around the standardized mean difference: Bootstrap and
#'   parametric confidence intervals.
#'   \emph{Educational and Psychological Measurement, 65}(1), 51--69.
#'   \doi{10.1177/0013164404264850}
#'
#' Kelley, K. (2007). Confidence intervals for standardized
#'   effect sizes: Theory, application, and implementation.
#'   \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
#'   obtaining precision: Delineating methods of sample size planning.
#'   \emph{Evaluation and the Health Professions, 26}(3), 258--287.
#'   \doi{10.1177/0163278703255242}
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
#' @section Warning:
#' This function uses \code{conf_limits_nct}, which has as one of its arguments \code{tol} (and can be modified with \code{tol} of the present function).
#' If the present function fails to converge (i.e., if it runs but does not report a solution), it is likely that the \code{tol} value is too restrictive and should be increased by a factor of 10, but probably by no more than 100.
#' Running the function \code{conf_limits_nct} directly will report the actual probability values of the limits found. This should be done if any modification to \code{tol} is necessary in order to ensure acceptable confidence limits for the noncentral \emph{t} parameter have been achieved.
#'
#' @seealso
#' \code{\link{smd}}, \code{\link{smd_c}}, \code{\link{ci_smd_c}},
#' \code{\link{ss_aipe_smd}}, \code{\link{ss_power_smd}},
#' \code{\link{plot_smd}}, \code{\link{conf_limits_nct}}
#'
#' @examples
#' # Steiger and Fouladi (1997) example values.
#' ci_smd(ncp = 2.6, n_1 = 10, n_2 = 10, conf_level = 1 - .05)
#' ci_smd(ncp = 2.4, n_1 = 300, n_2 = 300, conf_level = 1 - .05)
#'
#' @concept Cohen's d
#' @concept Hedges' g
#'
#' @keywords univar htest
#'
#' @family confidence intervals for effect sizes
#'
#' @export


ci_smd <- function(ncp = NULL, smd = NULL, n_1 = NULL, n_2 = NULL, conf_level = .95, alpha_lower = NULL, alpha_upper = NULL, tol = 1e-9, ...) {
  if (is.null(ncp) && is.null(smd)) stop("You must specify either the estimated noncentral parameter 'ncp' (generally the observed t-statistic) or the standardized mean difference 'smd' (as might be obtained from the 'smd' function.).", call. = FALSE)
  if (length(ncp) == 1 && length(smd) == 1) stop("You only need to specify either 'ncp' or 'smd', not both.", call. = FALSE)
  if (is.null(n_1) || is.null(n_2)) stop("You must specify sample size per group in order to determine confidence limits.", call. = FALSE)
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

  df <- n_1 + n_2 - 2

  if (length(ncp) == 1) {
    # if(ncp==0) stop("You need not use a noncentral method since the noncentrality parameter is zero; use the critical value from the central t-distribution.", call.=FALSE)
    smd <- ncp * sqrt((n_1 + n_2) / (n_1 * n_2))
    Limits <- conf_limits_nct(ncp, df, conf_level = conf_level, alpha_lower = alpha_lower, alpha_upper = alpha_upper, tol = tol, ...)

    Limits_L <- Limits[which(Limits$term == "lower_limit"), 2]
    Limits_U <- Limits[which(Limits$term == "upper_limit"), 2]
    Lower_Conf_Limit <- Limits_L * sqrt((n_1 + n_2) / (n_1 * n_2))
    Upper_Conf_Limit <- Limits_U * sqrt((n_1 + n_2) / (n_1 * n_2))
    term <- c("lower_limit", "smd", "upper_limit")
    value <- c(Lower_Conf_Limit, smd, Upper_Conf_Limit)
    out <- data.frame(term, value)
    attr(out, "conf_level") <- conf_level
    class(out) <- c("dmar_ci_smd", "dmar_tbl", "data.frame")
    return(out)
  }


  if (length(smd) == 1) {
    # if(smd==0) stop("You need not use a noncentral method since the effect size is zero; use the critical value from the central t-distribution.", call.=FALSE)
    ncp <- smd * sqrt((n_1 * n_2) / (n_1 + n_2))
    Limits <- conf_limits_nct(ncp, df, conf_level = conf_level, alpha_lower = alpha_lower, alpha_upper = alpha_upper, tol = tol, ...)
    Limits_L <- Limits[which(Limits$term == "lower_limit"), 2]
    Limits_U <- Limits[which(Limits$term == "upper_limit"), 2]
    Lower_Conf_Limit <- Limits_L * sqrt((n_1 + n_2) / (n_1 * n_2))
    Upper_Conf_Limit <- Limits_U * sqrt((n_1 + n_2) / (n_1 * n_2))
    term <- c("lower_limit", "smd", "upper_limit")
    value <- c(Lower_Conf_Limit, smd, Upper_Conf_Limit)
    out <- data.frame(term, value)
    attr(out, "conf_level") <- conf_level
    class(out) <- c("dmar_ci_smd", "dmar_tbl", "data.frame")
    return(out)
  }
}


#' Tidy / Glance Methods for ci_smd Output
#'
#' Returns the standard broom-style one-row summary of the SMD
#' confidence interval: \code{term} (always \code{"smd"}),
#' \code{estimate}, \code{ci_lower}, \code{ci_upper},
#' \code{conf_level}.
#'
#' @param x A \code{dmar_ci_smd} object returned by \code{ci_smd}.
#' @param \dots Unused.
#' @return A one-row \code{data.frame} with columns \code{term} (always
#'   \code{"smd"}), \code{estimate} (the point estimate of \emph{d}),
#'   \code{ci_lower} (the lower confidence limit on \eqn{\delta}),
#'   \code{ci_upper} (the upper confidence limit on \eqn{\delta}), and
#'   \code{conf_level} (the confidence level used in construction; \code{NA}
#'   when asymmetric \code{alpha_lower} / \code{alpha_upper} were supplied
#'   instead). Column names follow the broom convention so the result
#'   composes with the broom ecosystem.
#' @author Ken Kelley \email{kkelley@@nd.edu}
#' @examples
#' res <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
#' generics::tidy(res)
#' generics::glance(res)
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_ci_smd <- function(x, ...) {
  data.frame(
    term       = "smd",
    estimate   = x$value[x$term == "smd"],
    ci_lower   = x$value[x$term == "lower_limit"],
    ci_upper  = x$value[x$term == "upper_limit"],
    conf_level = attr(x, "conf_level") %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

#' @rdname tidy.dmar_ci_smd
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_ci_smd <- function(x, ...) {
  tidy.dmar_ci_smd(x, ...)
}

# Local %||% helper (also defined elsewhere for the reliability methods;
# redefining here keeps the file self-contained).
`%||%` <- function(a, b) if (is.null(a)) b else a
