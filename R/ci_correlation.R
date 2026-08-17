#' Confidence Intervals for the Population Correlation and Multiple Correlation
#'
#' @description
#' Two confidence intervals on the correlation scale share this page, named
#' by the convention that lowercase \emph{r} is the Pearson product-moment
#' correlation between two variables and capital \emph{R} is the multiple
#' correlation between an outcome and a set of predictors.
#'
#' \code{ci_r()} forms a confidence interval for the population correlation
#' coefficient \eqn{\rho}. The confidence interval is for the population
#' value \eqn{\rho}; the required input is the corresponding sample value,
#' the observed sample correlation coefficient \emph{r}. This approach
#' assumes that the two variables on which the correlation is based are
#' bivariate normally distributed (e.g., Hays, 1994, Chapter 14).
#'
#' \code{ci_R()} constructs a confidence interval for the population
#' multiple correlation coefficient \eqn{\rho = \sqrt{\rho^2}} from the
#' sample multiple correlation coefficient (or, equivalently, from the
#' observed \emph{F}-statistic and degrees of freedom). The interval is
#' obtained by inverting the sampling distribution of the sample \eqn{R^2}
#' and propagating the limits through the monotone (square root) transform.
#'
#' The two estimands meet at a single predictor: the multiple correlation
#' from a regression on one predictor is the absolute value of the Pearson
#' correlation between the outcome and that predictor.
#'
#' @param r Observed value of the sample correlation coefficient
#'   (specifically the zero-order Pearson product-moment correlation
#'   coefficient), for \code{ci_r()}
#' @param n Sample size for \code{ci_r()}, which must be at least 4 (see
#'   Details)
#' @param conf_level Confidence interval coverage (i.e., 1 - Type I error
#'   rate); default is .95
#' @param alpha_lower The Type I error rate for the lower confidence
#'   interval limit
#' @param alpha_upper The Type I error rate for the upper confidence
#'   interval limit
#' @param R Observed value of the sample multiple correlation coefficient,
#'   for \code{ci_R()}
#' @param df_1 Numerator degrees of freedom
#' @param df_2 Denominator degrees of freedom
#' @param random_predictors Whether or not the predictor variables are
#'   random or fixed (random is default)
#' @param F_value Obtained \emph{F}-value
#' @param N Sample size
#' @param p Number of predictors
#' @param \dots Allows one to potentially include parameter values for inner
#'   functions
#'
#' @details
#' \strong{The Pearson correlation interval (\code{ci_r}).} This approach
#' will not generally lead to a symmetric confidence interval. The function
#' first transforms \eqn{r} into \emph{Z}', forms a confidence interval for
#' the population value (i.e., \eqn{\zeta}), and then transforms the
#' confidence limits for \eqn{\zeta} into the scale of the correlation
#' coefficient. The interval requires a sample size of at least 4. The
#' variance of \emph{Z}' is \eqn{1/(n - 3)}, which is infinite at
#' \eqn{n = 3}; there the interval would be vacuous, covering \eqn{[-1, 1]}
#' regardless of \emph{r}, and for smaller \emph{n} the variance is
#' undefined. The function therefore stops with an error when \eqn{n < 4}.
#'
#' \strong{Fixed vs. random predictors (\code{ci_R}).} The two regression
#' models give \emph{different} sampling distributions for the sample
#' \eqn{R^2}, and so different confidence intervals on \eqn{\rho}. Under
#' fixed predictors the design matrix is treated as constant in
#' hypothetical replications of the study, and the omnibus
#' \eqn{F}-statistic follows a noncentral \emph{F} with \eqn{p} and
#' \eqn{N - p - 1} degrees of freedom and noncentrality
#' \eqn{\lambda = N \rho^2 / (1 - \rho^2)} (Cohen, 1988); the CI on
#' \eqn{\rho^2} is obtained by inverting that distribution and then taking
#' the square root (see \code{\link{conf_limits_ncf}}). Under random
#' predictors the design matrix is itself a draw from a joint multivariate
#' normal distribution and the unconditional sampling distribution of the
#' sample \eqn{R^2} is given by Lee (1971); the same Lee bisection that
#' \code{\link{ci_R2}} uses for the random-predictor CI on \eqn{\rho^2} is
#' applied here and the limits are mapped to \eqn{\rho}. Gatsonis and
#' Sampson (1989) document the comparison; in the behavioral, educational,
#' and social sciences predictor variables are almost always random, so the
#' default is \code{random_predictors = TRUE}. Pass
#' \code{random_predictors = FALSE} for designs in which the predictor
#' variables are fixed by design.
#'
#' @return
#' \code{ci_r()} returns a 3-row \code{data.frame} with columns \code{term}
#' and \code{value}. The \code{term} values are \code{"lower_limit"} (the
#' lower confidence limit on the population correlation \eqn{\rho}),
#' \code{"r"} (the observed sample correlation coefficient), and
#' \code{"upper_limit"} (the upper limit on \eqn{\rho}).
#'
#' \code{ci_R()} returns a 3-row \code{data.frame} with columns
#' \code{term}, \code{value}, \code{prob_less}, and \code{prob_greater}.
#' The rows are ordered \code{"lower_limit"}, \code{"R"} (the sample
#' multiple correlation coefficient supplied by the user, the point
#' estimate), and \code{"upper_limit"}, so the point estimate sits between
#' its confidence limits. The lower and upper limits are the confidence
#' limits on the population multiple correlation coefficient \eqn{\rho}
#' (square roots of the corresponding limits on \eqn{\rho^2}). The
#' \code{prob_less} and \code{prob_greater} columns report the achieved
#' lower-tail and upper-tail error probabilities at each limit (they are
#' \code{NA} for the \code{"R"} estimate row).
#'
#' @note The \code{ci_r()} confidence interval assumes that the two
#'   variables the correlation is based on are bivariate normal. See Hays
#'   (1994, Chapter 14) for details.
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
#' Hays, W. L. (1994). \emph{Statistics} (5th ed.). Fort Worth, TX:
#'   Harcourt Brace College Publishers.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#'   Theory, application, and implementation. \emph{Journal of Statistical
#'   Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' Kelley, K. (2008). Sample size planning for the squared multiple
#'   correlation coefficient: Accuracy in parameter estimation via narrow
#'   confidence intervals. \emph{Multivariate Behavioral Research, 43},
#'   524--555. \doi{10.1080/00273170802490632}
#'
#' Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
#'   regression: Obtaining regression coefficients that are accurate,
#'   not simply significant. \emph{Psychological Methods, 8}(3),
#'   305--321. \doi{10.1037/1082-989X.8.3.305}
#'
#' Lee, Y. S. (1971). Some results on the sampling distribution of the
#'   multiple correlation coefficient. \emph{Journal of the Royal
#'   Statistical Society, Series B, 33}(1), 117--130.
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027).
#'   \emph{Designing experiments and analyzing data: A model comparison
#'   perspective} (4th ed.). Routledge. (See Chapter 3 on \eqn{R^2} as a
#'   model comparison effect size.)
#'
#' Smithson, M. (2003). \emph{Confidence intervals}. Thousand Oaks, CA: Sage Publications.
#'
#' Steiger, J. H. (2004). Beyond the \emph{F} test: Effect size confidence intervals and tests of close fit in the analysis of variance and contrast analysis. \emph{Psychological Methods, 9}(2), 164--182.
#'   \doi{10.1037/1082-989X.9.2.164}
#'
#' Steiger, J. H., & Fouladi, R. T. (1992). R2: A computer program for
#'   interval estimation, power calculations, sample size estimation, and
#'   hypothesis testing in multiple regression. \emph{Behavior Research
#'   Methods, Instruments, & Computers, 24}(4), 581--582.
#'   \doi{10.3758/BF03203611}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # Pearson correlation, from Hays. Suppose n = 100 and r = .35.
#' ci_r(r = .35, n = 100, conf_level = .95)
#'
#' # Here is another way to enter the above example.
#' ci_r(r = .35, n = 100, conf_level = NULL,
#'      alpha_lower = .025, alpha_upper = .025)
#'
#' # Here are examples of one-sided confidence intervals.
#' ci_r(r = .35, n = 100, conf_level = NULL, alpha_lower = 0, alpha_upper = .05)
#' ci_r(r = .35, n = 100, conf_level = NULL, alpha_lower = .05, alpha_upper = 0)
#'
#' # Multiple correlation from a five-predictor regression.
#' ci_R(R = .7071, df_1 = 5, df_2 = 50, conf_level = .95,
#'      random_predictors = TRUE)
#'
#' @seealso \code{\link{ci_R2}}, \code{\link{ss_aipe_r}},
#'   \code{\link{ss_power_r}}, \code{\link{var_r}}, \code{\link{ss_aipe_R2}},
#'   \code{\link{convert_r_Z}}, \code{\link{convert_Z_r}},
#'   \code{\link{conf_limits_nct}}
#'
#' @keywords design regression
#'
#' @family confidence intervals for effect sizes
#'
#' @name ci_correlation
NULL

#' @rdname ci_correlation
#' @export
ci_r <- function(r, n, conf_level = .95, alpha_lower = NULL, alpha_upper = NULL) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 4) {
    if (is.numeric(n) && length(n) == 1L && !is.na(n) && n == 3) {
      stop("'n' must be at least 4. At n = 3 the Fisher's Z interval is ",
           "vacuous: the variance of Z is 1/(n - 3), which is infinite at ",
           "n = 3, so the interval is [-1, 1] regardless of 'r' and carries ",
           "no information about the population correlation.", call. = FALSE)
    }
    stop("'n' must be a single sample size of at least 4; the Fisher's Z ",
         "interval is based on the variance 1/(n - 3).", call. = FALSE)
  }

  if (!is.null(conf_level)) {
    if (conf_level >= 1 || conf_level <= 0) stop("You have not properly specified 'conf_level'", call. = FALSE)
    if (!is.null(alpha_lower)) stop("You specified both 'conf_level' and 'alpha_lower', specify confidence level using only one approach.", call. = FALSE)
    if (!is.null(alpha_upper)) stop("You specified both 'conf_level' and 'alpha_upper', specify confidence level using only one approach.", call. = FALSE)
    alpha_lower <- alpha_upper <- (1 - conf_level) / 2
  }

  if (is.null(conf_level)) {
    if (is.null(alpha_lower) || is.null(alpha_upper)) stop("You need to specify either 'conf_level', or both 'alpha_lower' and 'alpha_upper'.", call. = FALSE)
    if (alpha_lower > 0.5 || alpha_lower < 0) stop("You have not properly specified 'alpha_lower' correctly.", call. = FALSE)
    if (alpha_upper > 0.5 || alpha_upper < 0) stop("You have not properly specified 'alpha_upper' correctly.", call. = FALSE)
  }

  CV_Lower <- qnorm(1 - alpha_lower)
  CV_Upper <- qnorm(1 - alpha_upper)

  Z <- convert_r_Z(r)[1, 2]
  SE_Z <- sqrt(1 / (n - 3))

  CI_Lower_Zeta <- Z - CV_Lower * SE_Z
  CI_Upper_Zeta <- Z + CV_Upper * SE_Z

  if (alpha_lower > 0) CI_Lower_rho <- convert_Z_r(CI_Lower_Zeta)[1, 2]
  if (alpha_upper > 0) CI_Upper_rho <- convert_Z_r(CI_Upper_Zeta)[1, 2]

  if (alpha_lower == 0) CI_Lower_rho <- -1
  if (alpha_upper == 0) CI_Upper_rho <- 1

  term <- c("lower_limit", "r", "upper_limit")
  value <- c(CI_Lower_rho, r, CI_Upper_rho)

  out <- data.frame(term, value); attr(out, "conf_level") <- conf_level; class(out) <- c("dmar_ci_long", "dmar_tbl", "data.frame"); return(out)
}

#' @rdname ci_correlation
#' @export
ci_R <- function(R = NULL, df_1 = NULL, df_2 = NULL, conf_level = .95, random_predictors = TRUE,
                 F_value = NULL, N = NULL, p = NULL, alpha_lower = NULL, alpha_upper = NULL, ...) {
  if (!is.null(R)) {
    if (R < 0) stop("Your multiple correlation coefficient ('R') cannot be less than zero.")
    if (R > 1) stop("Your multiple correlation coefficient ('R') cannot be greater than one.")
  }

  if (is.null(R)) {
    if (is.null(df_1)) {
      if (is.null(p)) stop("You need to specify 'p' or 'df_1'.")
      df_1 <- p
    }

    if (is.null(df_2)) {
      if (is.null(N)) stop("You need to specify 'N' or 'df_2'.")
      if (is.null(p)) stop("You need to specify 'p' or 'df_2'.")
      df_2 <- N - p - 1
    }

    R <- sqrt(convert_f_R2(F_value = F_value, df_1 = df_1, df_2 = df_2)[1, 2])

    N <- NULL
    p <- NULL
  }


  Limits <- ci_R2(
    R2 = R^2, df_1 = df_1, df_2 = df_2, conf_level = conf_level, random_predictors = random_predictors,
    F_value = NULL, N = N, p = p, alpha_lower = alpha_lower, alpha_upper = alpha_upper, ...
  )

  ll_row <- which(Limits$term == "lower_limit")
  ul_row <- which(Limits$term == "upper_limit")
  return(.as_dmar_tbl(data.frame(
    term = c("lower_limit", "R", "upper_limit"),
    value = c(sqrt(Limits$value[ll_row]), R, sqrt(Limits$value[ul_row])),
    prob_less = c(Limits[ll_row, 3], NA_real_, 1 - Limits[ul_row, 4]),
    prob_greater = c(1 - Limits[ll_row, 3], NA_real_, Limits[ul_row, 4])
  ), conf_level = conf_level))
}
