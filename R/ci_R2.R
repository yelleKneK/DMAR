#' Confidence Interval for the Population Squared Multiple Correlation Coefficient
#'
#' Constructs a confidence interval for the population squared multiple
#' correlation coefficient \eqn{\rho^2} by inverting the sampling
#' distribution of the sample \eqn{R^2}. The confidence interval is for the
#' population value \eqn{\rho^2}; the required input is the corresponding
#' sample value, the observed sample squared multiple correlation coefficient
#' \eqn{R^2} (or, equivalently, the observed \emph{F}-statistic and degrees of
#' freedom). The right choice of sampling
#' distribution, and so the right interval, depends on whether the
#' predictors are treated as random draws from a joint multivariate normal
#' distribution (the default and the typical case in the behavioral,
#' educational, and social sciences) or as fixed by design (planned dosing
#' levels, factorial covariates, etc.). The function selects the sampling
#' distribution via \code{random_predictors} and inverts the corresponding
#' noncentral distribution; the construction is the regression analogue of
#' a noncentral distribution based CI on a standardized effect size
#' (Steiger & Fouladi, 1992; Kelley, 2007).
#'
#' @param R2 Observed value of the sample squared multiple correlation coefficient
#' @param df_1 Numerator degrees of freedom
#' @param df_2 Denominator degrees of freedom
#' @param conf_level Confidence interval coverage; 1-Type I error rate
#' @param random_predictors Whether or not the predictor variables are random or fixed (random is default)
#' @param F_value Obtained \emph{F}-value
#' @param N Sample size
#' @param p Number of predictors
#' @param alpha_lower Type I error for the lower confidence limit
#' @param alpha_upper Type I error for the upper confidence limit
#' @param tol The convergence tolerance passed to \code{\link[stats]{uniroot}}
#'   when \code{random_predictors = FALSE} and the confidence limits are found
#'   by inverting the noncentral \emph{F} distribution (see
#'   \code{\link{conf_limits_ncf}}); ignored when
#'   \code{random_predictors = TRUE}, where the Lee (1971) bisection uses its
#'   own fixed tolerance
#'
#' @details
#' \strong{Fixed vs.\ random predictors.} The two regression models give
#' \emph{different} sampling distributions for the sample \eqn{R^2}, and
#' so different confidence intervals. Under fixed predictors the design
#' matrix is treated as constant in hypothetical replications of the
#' study, and the omnibus \eqn{F}-statistic
#' \eqn{F = (R^2 / p) / ((1 - R^2) / (N - p - 1))} follows a noncentral
#' \emph{F} with \eqn{p} and \eqn{N - p - 1} degrees of freedom and
#' noncentrality \eqn{\lambda = N \rho^2 / (1 - \rho^2)} (Cohen, 1988); the
#' CI is obtained by inverting that distribution at the supplied
#' confidence level (see \code{\link{conf_limits_ncf}}). Under random
#' predictors the design matrix is itself a draw from a joint multivariate
#' normal distribution and the unconditional sampling distribution of the
#' sample \eqn{R^2} is given by Lee (1971); \code{ci_R2} uses the Lee
#' (1971) bisection (the same construction Algina and Olejnik 2000
#' implemented in SAS) to invert that distribution. Gatsonis and Sampson
#' (1989) document the comparison and show that treating random predictors
#' as fixed tends to over-state precision (and so under-state the CI
#' width); the discrepancy is modest at moderate to large \eqn{N} but
#' non-trivial at small \eqn{N} with moderate-to-large effects. In the
#' behavioral, educational, and social sciences predictor variables are
#' almost always random, so the default is
#' \code{random_predictors = TRUE}; pass \code{random_predictors = FALSE}
#' for designs in which the predictor variables are fixed by design.
#'
#' @return
#' A 3-row \code{data.frame} with columns \code{term}, \code{value},
#' \code{prob_less}, and \code{prob_greater}. The rows are ordered
#' \code{"lower_limit"} (lower confidence limit on the population
#' \eqn{\rho^2}), \code{"R2"} (the sample squared multiple correlation
#' coefficient supplied by the user, the point estimate), and
#' \code{"upper_limit"} (upper confidence limit on the population
#' \eqn{\rho^2}), so the point estimate sits between its confidence limits.
#' The \code{prob_less} and \code{prob_greater} columns report the achieved
#' lower-tail and upper-tail error probabilities at each limit (they are
#' \code{NA} for the \code{"R2"} estimate row). For random-predictor mode
#' (\code{random_predictors = TRUE}) the limits are computed via the Lee
#' (1971) bisection over the multiple-correlation sampling distribution; for
#' fixed-predictor mode they are computed by inversion of the noncentral
#' \emph{F} distribution.
#'
#' @references
#' Algina, J. & Olejnik, S. (2000). Determining sample size for accurate estimation of the squared multiple correlation coefficient. \emph{Multivariate Behavioral Research, 35}, 119--137.
#'   \doi{10.1207/s15327906mbr3501_5}
#'
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#'   sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact
#'   power and sample size calculations. \emph{Psychological Bulletin,
#'   106}(3), 516--524.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple
#'   correlation coefficient: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{Multivariate Behavioral Research, 43},
#'   524--555. \doi{10.1080/00273170802490632}
#'
#' Lee, Y. S. (1971). Some results on the sampling distribution of the
#'   multiple correlation coefficient. \emph{Journal of the Royal
#'   Statistical Society, Series B, 33}(1), 117--130.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 3 on \eqn{R^2} as a model comparison
#'   effect size.)
#'
#' Smithson, M. (2003). \emph{Confidence intervals}. Thousand Oaks, CA: Sage Publications.
#'
#' Steiger, J. H., & Fouladi, R. T. (1992). R2: A computer program for
#'   interval estimation, power calculations, sample size estimation, and
#'   hypothesis testing in multiple regression. \emph{Behavior Research
#'   Methods, Instruments, & Computers, 24}(4), 581--582.
#'   \doi{10.3758/BF03203611}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{ss_aipe_R2}}, \code{\link{conf_limits_ncf}}
#'
#' @examples
#' # For random predictor variables.
#' ci_R2(R2 = .25, N = 100, p = 5, conf_level = .95, random_predictors = TRUE)
#'
#' ci_R2(F_value = 6.266667, N = 100, p = 5, conf_level = .95, random_predictors = TRUE)
#'
#' # For fixed predictor variables.
#' ci_R2(R2 = .25, N = 100, p = 5, conf_level = .95, random_predictors = FALSE)
#'
#' ci_R2(F_value = 6.266667, N = 100, p = 5, conf_level = .95, random_predictors = FALSE)
#'
#' # One sided confidence intervals when predictors are random.
#' ci_R2(R2 = .25, N = 100, p = 5, alpha_lower = .05, alpha_upper = 0,
#'       conf_level = NULL, random_predictors = TRUE)
#'
#' ci_R2(R2 = .25, N = 100, p = 5, alpha_lower = 0, alpha_upper = .05,
#'       conf_level = NULL, random_predictors = TRUE)
#'
#' # One sided confidence intervals when predictors are fixed.
#' ci_R2(R2 = .25, N = 100, p = 5, alpha_lower = .05, alpha_upper = 0,
#'       conf_level = NULL, random_predictors = FALSE)
#'
#' ci_R2(R2 = .25, N = 100, p = 5, alpha_lower = 0, alpha_upper = .05,
#'       conf_level = NULL, random_predictors = FALSE)
#'
#' @keywords multivariate htest regression
#'
#' @family confidence intervals for effect sizes
#'
#' @export

ci_R2 <- function(R2 = NULL, df_1 = NULL, df_2 = NULL, conf_level = .95, random_predictors = TRUE,
                  F_value = NULL, N = NULL, p = NULL, alpha_lower = NULL, alpha_upper = NULL, tol = 1e-9) {

  # Local helper: tag every return value with the dmar_ci_R2 class
  # so generics::tidy / generics::glance dispatch. Carries R2 and
  # conf_level on attributes so the methods can reach them.
  # conf_level is nulled out below once converted to symmetric
  # alpha_lower / alpha_upper, so recover it from those when needed.
  add_class <- function(df) {
    attr(df, "R2") <- R2
    cl <- conf_level
    if (is.null(cl) && !is.null(alpha_lower) && !is.null(alpha_upper)) {
      cl <- 1 - alpha_lower - alpha_upper
    }
    attr(df, "conf_level") <- cl
    class(df) <- c("dmar_ci_R2", "dmar_tbl", "data.frame")
    df
  }


  if ((!is.null(N) || !is.null(p)) && (!is.null(df_1) || !is.null(df_2))) stop("Either specify \'df_1\' and \'df_2\' or \'N\' and \'p,\' but not both combinations.")

  if (!is.null(N) && !is.null(p) && is.null(df_1) && is.null(df_2)) {
    df_1 <- p
    df_2 <- N - p - 1
  }

  if (!is.null(df_1) && !is.null(df_2) && is.null(N) && is.null(p)) {
    N <- df_1 + df_2 + 1
    p <- df_1
  }

  # Supplying alpha_lower / alpha_upper defines the coverage as
  # 1 - alpha_lower - alpha_upper, so the default conf_level must not also be
  # attached (its footer would mislabel, for example, a 90% interval as 95%).
  # Reject an explicitly supplied conf_level mixed with the alphas; otherwise
  # drop the default conf_level and let add_class() recover coverage from them.
  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (alphas_supplied) {
    if (!missing(conf_level) && !is.null(conf_level)) {
      stop("Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'; you cannot mix them.", call. = FALSE)
    }
    conf_level <- NULL
  }

  if (!is.null(conf_level) && is.null(alpha_lower) && is.null(alpha_upper)) {
    if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level) ||
        conf_level <= 0 || conf_level >= 1) {
      stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
    }
    alpha_lower <- alpha_upper <- (1 - conf_level) / 2
    conf_level <- NULL
  }

  # Use the fast-path internals on the hot path to skip the per-call
  # data.frame allocation in convert_R2_f / convert_f_R2 /
  # convert_lambda_R2; behavior is identical, the public functions
  # are unchanged. ci_R2 is called thousands of times per
  # ss_aipe_R2() invocation, so the difference matters.
  if (is.null(F_value)) {
    F_value <- .convert_R2_f_fast(R2, df_1, df_2)
  }

  if (is.null(R2)) {
    R2 <- .convert_f_R2_fast(F_value, df_1, df_2)
  }

  if (random_predictors == FALSE) {
    Limits <- .conf_limits_ncf_for(
      caller = "ci_R2",
      quantity = "the population squared multiple correlation coefficient",
      F_value = F_value, df_1 = df_1, df_2 = df_2, conf_level = NULL,
      tol = tol, alpha_lower = alpha_lower, alpha_upper = alpha_upper
    )

    lower_lambda <- Limits$value[Limits$term == "lower_limit"]
    upper_lambda <- Limits$value[Limits$term == "upper_limit"]
    achieved_alpha_lower <- Limits$prob_greater[Limits$term == "lower_limit"]
    achieved_alpha_upper <- Limits$prob_less[Limits$term == "upper_limit"]

    LL <- .convert_lambda_R2_fast(lower_lambda, N)
    UL <- if (is.infinite(upper_lambda)) 1 else .convert_lambda_R2_fast(upper_lambda, N)

    return(add_class(data.frame(
      term = c("lower_limit", "R2", "upper_limit"),
      value = c(LL, R2, UL),
      prob_less = c(achieved_alpha_lower, NA_real_, 1 - achieved_alpha_upper),
      prob_greater = c(1 - achieved_alpha_lower, NA_real_, achieved_alpha_upper)
    )))
  }


  if (random_predictors == TRUE) {
    pul <- alpha_upper
    pll <- 1 - alpha_lower

    # P(R^2_sample <= R2_obs | rho^2) is monotone decreasing in rho^2, so its
    # largest value is at rho^2 = 0. When a limit's target probability exceeds
    # that maximum, the inversion has no solution in [0, 1) and the rho^2 limit
    # is 0 (the observed R^2 lies at or below that quantile of the null). This
    # single test guards both limits symmetrically; the previous code guarded
    # only the lower limit and let the upper limit walk to the boundary
    # artifact rho^2 ~ 1 whenever R^2 fell below its alpha_upper null quantile.
    C0 <- .lee_random_R2_cdf(R2_obs = R2, rho2 = 0, N = N, p = p)

    ulrhosq <- if (pul == 0) 1 else if (C0 <= pul) 0 else
      .lee_random_R2_bisect(R2_obs = R2, p_target = pul, N = N, p = p)

    llrhosq <- if (pll == 1) 0 else if (C0 <= pll) 0 else
      .lee_random_R2_bisect(R2_obs = R2, p_target = pll, N = N, p = p)

    if (llrhosq > ulrhosq) warning("There is a problem; the lower limit is greater than the upper limit (Are you at one of the boundaries of Rho^2 or is alpha very large?).")

    if (pll == 1) {
      return(add_class(data.frame(
        term = c("lower_limit", "R2", "upper_limit"), value = c(0, R2, ulrhosq),
        prob_less = c(0, NA_real_, 1 - pul), prob_greater = c(1, NA_real_, pul)
      )))
    }
    if (pul == 0) {
      return(add_class(data.frame(
        term = c("lower_limit", "R2", "upper_limit"), value = c(llrhosq, R2, 1),
        prob_less = c(1 - pll, NA_real_, 1), prob_greater = c(pll, NA_real_, 0)
      )))
    }
    return(add_class(data.frame(
      term = c("lower_limit", "R2", "upper_limit"), value = c(llrhosq, R2, ulrhosq),
      prob_less = c(1 - pll, NA_real_, 1 - pul), prob_greater = c(pll, NA_real_, pul)
    )))
  }
}


#' Tidy / Glance Methods for ci_R2 Output
#'
#' Returns the standard broom-style one-row summary of the
#' \eqn{R^2} confidence interval: \code{term} (always
#' \code{"R2"}), \code{estimate}, \code{ci_lower},
#' \code{ci_upper}, \code{conf_level}.
#'
#' @param x A \code{dmar_ci_R2} object returned by \code{ci_R2}.
#' @param \dots Unused.
#' @return A one-row \code{data.frame} in broom convention.
#' @author Ken Kelley \email{kkelley@@nd.edu}
#' @examples
#' res <- ci_R2(R2 = 0.25, N = 100, p = 5)
#' generics::tidy(res)
#' generics::glance(res)
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_ci_R2 <- function(x, ...) {
  data.frame(
    term       = "R2",
    estimate   = attr(x, "R2"),
    ci_lower   = x$value[x$term == "lower_limit"],
    ci_upper  = x$value[x$term == "upper_limit"],
    conf_level = attr(x, "conf_level"),
    stringsAsFactors = FALSE
  )
}

#' @rdname tidy.dmar_ci_R2
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_ci_R2 <- function(x, ...) {
  tidy.dmar_ci_R2(x, ...)
}
