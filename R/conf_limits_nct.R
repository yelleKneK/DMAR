#' Confidence Limits for a Noncentrality Parameter From a \emph{t}-distribution
#'
#' Finds the noncentrality parameters of a noncentral \emph{t}-distribution
#' that bracket an observed \emph{t}-value with the requested tail
#' probabilities, giving a confidence interval on the population noncentrality
#' parameter. Together with \code{\link{conf_limits_ncf}} and
#' \code{\link{conf_limits_nc_chisq}}, this is one of the low-level
#' noncentral distribution workhorses on which the \code{ci_*} confidence
#' interval functions (e.g., \code{\link{ci_smd}}, \code{\link{ci_smd_c}},
#' \code{\link{ci_cv}}) are built; most analyses reach it through those
#' functions rather than calling it directly.
#'
#' @param ncp The noncentrality parameter (e.g., observed \emph{t}-value) of interest
#' @param df The degrees of freedom
#' @param conf_level The level of confidence for a symmetric confidence interval
#' @param alpha_lower The proportion of values beyond the lower limit of the confidence interval (cannot be used with \code{conf_level})
#' @param alpha_upper The proportion of values beyond the upper limit of the confidence interval (cannot be used with \code{conf_level})
#' @param t_value Alias for \code{ncp}
#' @param tol The convergence tolerance passed to \code{\link[stats]{uniroot}} when locating each limit
#' @param verbose If \code{TRUE} (the default), the returned data frame additionally reports the achieved tail probabilities at each limit; if \code{FALSE}, only \code{term} and \code{value} are returned
#' @param \dots Additional arguments forwarded to \code{\link[stats]{uniroot}}
#'
#' @details
#' Each confidence limit is the noncentrality parameter of a noncentral
#' \emph{t}-distribution with \code{df} degrees of freedom whose appropriate tail
#' at the observed \code{ncp} contains the requested probability:
#' \itemize{
#'   \item the lower limit satisfies \eqn{P(T \ge \mathrm{ncp}) = \mathtt{alpha\_lower}};
#'   \item the upper limit satisfies \eqn{P(T \le \mathrm{ncp}) = \mathtt{alpha\_upper}}.
#' }
#' Each tail probability is continuous and strictly monotone in the
#' noncentrality parameter, so each limit is the unique root of a
#' one-dimensional equation. The roots are located with
#' \code{\link[stats]{uniroot}} starting from a bracket centered on \code{ncp}
#' with half-width scaled by the asymptotic standard error of the noncentrality
#' estimator; \code{extendInt} is used to widen the bracket if needed.
#'
#' This function is especially useful for forming confidence intervals around
#' standardized mean differences (Cohen's \emph{d}, Glass's \emph{g}, Hedges' \emph{g}),
#' standardized regression coefficients, and coefficients of variation.
#'
#' @return
#' A \code{data.frame} with one row per confidence limit and the columns:
#' \item{term}{Either \code{"lower_limit"} or \code{"upper_limit"}.}
#' \item{value}{The noncentrality parameter at that limit. \code{-Inf} when \code{alpha_lower = 0}; \code{Inf} when \code{alpha_upper = 0}.}
#' \item{prob_less}{(\code{verbose = TRUE}) The probability \eqn{P(T \le \mathrm{ncp})} that a \emph{t}-statistic from the noncentral \emph{t}-distribution centered at the row's limit falls at or below the observed \code{ncp}. By construction this equals \code{alpha_upper} on the \code{upper_limit} row and \eqn{1 - \mathtt{alpha\_lower}} on the \code{lower_limit} row.}
#' \item{prob_greater}{(\code{verbose = TRUE}) The complementary probability \eqn{P(T \ge \mathrm{ncp})}. By construction this equals \code{alpha_lower} on the \code{lower_limit} row and \eqn{1 - \mathtt{alpha\_upper}} on the \code{upper_limit} row.}
#'
#' @references
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#'   calculation of confidence intervals that are based on central and
#'   noncentral distributions. \emph{Educational and Psychological
#'   Measurement, 61}(4), 532--574. \doi{10.1177/0013164401614002}
#'
#' Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the standardized mean difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological Measurement, 65}, 51--69.
#'   \doi{10.1177/0013164404264850}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20}(8), 1--24.
#'   \doi{10.18637/jss.v020.i08}
#'
#' Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
#'   estimation and the evaluation of statistical methods. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @section Warning:
#' As of R 4.0.0, the largest \code{ncp} that R can accurately handle is 37.62.
#'
#' @seealso
#' \code{\link[stats:TDist]{stats::pt()}}, \code{\link[stats:TDist]{stats::qt()}}, \code{\link[stats]{uniroot}}, \code{\link{ci_smd}}, \code{\link{ci_smd_c}}, \code{\link{conf_limits_ncf}}, \code{\link{conf_limits_nc_chisq}}
#'
#' @examples
#'# Suppose observed t-value based on 'df'=126 is 2.83. Finding the lower
#'# and upper critical values for the population noncentrality parameter
#'# with a symmetric confidence interval with 95\% confidence is given as:
#'conf_limits_nct(ncp = 2.83, df = 126, conf_level = .95)
#'
#'# Modifying the above example so that a nonsymmetric 95% confidence interval
#'# can be formed:
#'conf_limits_nct(ncp = 2.83, df = 126, alpha_lower = .01, alpha_upper = .04, conf_level = NULL)
#'
#'# Modifying the above example so that a single-sided 95% confidence interval
#'# can be formed:
#'conf_limits_nct(ncp = 2.83, df = 126, alpha_lower = 0, alpha_upper = .05, conf_level = NULL)
#'
#' @keywords models htest
#'
#' @family noncentral distribution confidence limits
#'
#' @export
#' @importFrom stats pt uniroot

conf_limits_nct <- function(ncp, df, conf_level = .95, alpha_lower = NULL, alpha_upper = NULL,
                            t_value, tol = 1e-9, verbose = TRUE, ...) {
  if (missing(ncp)) {
    if (missing(t_value)) {
      stop("You need to specify either 'ncp' or its alias, 't_value'; you have not specified either.", call. = FALSE)
    }
    ncp <- t_value
  }

  if (df <= 0) stop("The degrees of freedom must be some positive value.", call. = FALSE)

  if (abs(ncp) > 37.62) {
    warning(
      "The observed noncentrality parameter exceeds 37.62 in magnitude, which is the limit at which R's pt()/qt() can return accurate noncentral t probabilities. Results may be inaccurate; use caution.",
      call. = FALSE
    )
  }

  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (!is.null(conf_level) && alphas_supplied) {
    stop("Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'; you cannot mix them.", call. = FALSE)
  }
  if (is.null(conf_level)) {
    if (is.null(alpha_lower) || is.null(alpha_upper)) {
      stop("When 'conf_level' is NULL, both 'alpha_lower' and 'alpha_upper' must be specified.", call. = FALSE)
    }
  } else {
    if (conf_level <= 0 || conf_level >= 1) stop("'conf_level' must be between 0 and 1.", call. = FALSE)
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }

  if (alpha_lower < 0 || alpha_upper < 0 || alpha_lower + alpha_upper >= 1) {
    stop("'alpha_lower' and 'alpha_upper' must be non-negative and sum to less than 1.", call. = FALSE)
  }

  # Initial bracket half-width: a generous multiple of the asymptotic SE of the
  # noncentrality estimator, sqrt(1 + ncp^2 / (2 * df)). uniroot's extendInt
  # widens the bracket if a sign change is not initially present.
  bracket_half_width <- max(10, 8 * sqrt(1 + ncp^2 / (2 * df)))

  solve_for_lower_limit <- function(target_alpha) {
    objective <- function(candidate_ncp) suppressWarnings(pt(ncp, df = df, ncp = candidate_ncp, lower.tail = FALSE)) - target_alpha
    suppressWarnings(uniroot(
      objective,
      interval = c(ncp - bracket_half_width, ncp),
      tol = tol,
      extendInt = "upX",
      ...
    ))$root
  }

  solve_for_upper_limit <- function(target_alpha) {
    objective <- function(candidate_ncp) suppressWarnings(pt(ncp, df = df, ncp = candidate_ncp, lower.tail = TRUE)) - target_alpha
    suppressWarnings(uniroot(
      objective,
      interval = c(ncp, ncp + bracket_half_width),
      tol = tol,
      extendInt = "downX",
      ...
    ))$root
  }

  if (alpha_lower == 0) {
    lower_ncp <- -Inf
    achieved_alpha_lower <- 0
  } else {
    lower_ncp <- solve_for_lower_limit(alpha_lower)
    achieved_alpha_lower <- suppressWarnings(pt(ncp, df = df, ncp = lower_ncp, lower.tail = FALSE))
  }

  if (alpha_upper == 0) {
    upper_ncp <- Inf
    achieved_alpha_upper <- 0
  } else {
    upper_ncp <- solve_for_upper_limit(alpha_upper)
    achieved_alpha_upper <- suppressWarnings(pt(ncp, df = df, ncp = upper_ncp, lower.tail = TRUE))
  }

  # prob_less and prob_greater report the noncentral t probability mass on each
  # side of the observed `ncp` for the distribution centered at that limit.
  # By construction, prob_greater at the lower_limit equals alpha_lower and
  # prob_less at the upper_limit equals alpha_upper.
  term <- c("lower_limit", "upper_limit")
  value <- c(lower_ncp, upper_ncp)
  prob_less <- c(1 - achieved_alpha_lower, achieved_alpha_upper)
  prob_greater <- c(achieved_alpha_lower, 1 - achieved_alpha_upper)

  if (verbose) {
    data.frame(term, value, prob_less, prob_greater)
  } else {
    data.frame(term, value)
  }
}
