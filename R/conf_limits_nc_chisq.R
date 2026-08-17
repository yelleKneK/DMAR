#' Confidence Limits for the Noncentrality Parameter of a Noncentral Chi Square Distribution
#'
#' Finds the noncentrality parameters of a noncentral chi square distribution
#' that bracket an observed chi square value with the requested tail
#' probabilities, giving a confidence interval on the population noncentrality
#' parameter. Together with \code{\link{conf_limits_nct}} and
#' \code{\link{conf_limits_ncf}}, this is one of the low-level noncentral
#' distribution workhorses on which the \code{ci_*} confidence interval
#' functions are built; most analyses reach it through those functions rather
#' than calling it directly.
#'
#' @param chi_square The observed chi square value
#' @param conf_level The desired degree of confidence for a symmetric interval
#' @param df The degrees of freedom
#' @param alpha_lower The proportion of values beyond the lower limit (cannot be used with \code{conf_level})
#' @param alpha_upper The proportion of values beyond the upper limit (cannot be used with \code{conf_level})
#' @param tol The convergence tolerance passed to \code{\link[stats]{uniroot}}
#' @param verbose If \code{TRUE} (the default), the returned data frame additionally reports the achieved tail probabilities at each limit; if \code{FALSE}, only \code{term} and \code{value} are returned
#' @param \dots Additional arguments forwarded to \code{\link[stats]{uniroot}}
#'
#' @details
#' Each confidence limit is the noncentrality parameter \eqn{\lambda \ge 0} of a
#' noncentral chi square distribution with \code{df} degrees of freedom whose
#' appropriate tail at the observed \code{chi_square} contains the requested
#' probability:
#' \itemize{
#'   \item the lower limit satisfies \eqn{P(X \ge \mathtt{chi\_square}) = \mathtt{alpha\_lower}};
#'   \item the upper limit satisfies \eqn{P(X \le \mathtt{chi\_square}) = \mathtt{alpha\_upper}}.
#' }
#' The two conditions run in opposite directions in \eqn{\lambda}: the
#' lower-tail probability \eqn{P(X \le \mathtt{chi\_square})} is continuous
#' and strictly decreasing in the noncentrality parameter, so the upper-tail
#' probability \eqn{P(X \ge \mathtt{chi\_square})} is continuous and strictly
#' increasing in it. The lower limit is the \eqn{\lambda} at which the upper
#' tail has grown to \code{alpha_lower}, and the upper limit is the
#' \eqn{\lambda} at which the lower tail has shrunk to \code{alpha_upper}.
#' Each is therefore the unique non-negative root of a one-dimensional
#' equation, and both are located with \code{\link[stats]{uniroot}} on the
#' decreasing lower-tail scale; \code{extendInt} is used to widen the search
#' bracket if needed.
#'
#' Because the noncentrality parameter is bounded below by zero, the lower limit
#' is set to zero whenever the observed \code{chi_square} is smaller than the
#' \code{alpha_lower} critical value of the central chi square distribution
#' (i.e., the data is consistent with \eqn{\lambda = 0} at the requested
#' confidence level). A warning is issued in that case, and the achieved
#' probabilities reported in the output reflect the actual values at
#' \eqn{\lambda = 0}.
#'
#' Symmetrically, when the observed \code{chi_square} is so small that even
#' at \eqn{\lambda = 0} the lower-tail probability is already at or below
#' \code{alpha_upper}, no \eqn{\lambda \ge 0} places as much as
#' \code{alpha_upper} mass at or below \code{chi_square}; the upper limit is
#' undefined and is returned as \code{NA}, with a warning.
#'
#' @return
#' A \code{data.frame} with one row per confidence limit and the columns:
#' \item{term}{Either \code{"lower_limit"} or \code{"upper_limit"}.}
#' \item{value}{The noncentrality parameter at that limit. \code{0} when \code{alpha_lower = 0} or the lower limit is unattainable; \code{Inf} when \code{alpha_upper = 0}; \code{NA} when the observed \code{chi_square} is so small that even at \eqn{\lambda = 0} the lower-tail probability is already at or below \code{alpha_upper}, so the upper limit is undefined (a warning is issued).}
#' \item{prob_less}{(\code{verbose = TRUE}) The probability \eqn{P(X \le \mathtt{chi\_square})} that an observation from the noncentral chi square distribution centered at the row's limit falls at or below the observed \code{chi_square}.}
#' \item{prob_greater}{(\code{verbose = TRUE}) The complementary probability \eqn{P(X \ge \mathtt{chi\_square})}. By construction this equals \code{alpha_lower} on the \code{lower_limit} row and \eqn{1 - \mathtt{alpha\_upper}} on the \code{upper_limit} row.}
#'
#' @references
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#'   calculation of confidence intervals that are based on central and
#'   noncentral distributions. \emph{Educational and Psychological
#'   Measurement, 61}(4), 532--574. \doi{10.1177/0013164401614002}
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. \emph{Journal of Statistical
#' Software, 20}(8), 1--24. \doi{10.18637/jss.v020.i08}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{conf_limits_nct}}, \code{\link{conf_limits_ncf}}, \code{\link[stats:Chisquare]{stats::pchisq()}}, \code{\link[stats:Chisquare]{stats::qchisq()}}, \code{\link[stats]{uniroot}}
#'
#' @examples
#' # A typical call to the function.
#' conf_limits_nc_chisq(chi_square = 30, conf_level = .95, df = 15)
#'
#' # A one-sided (upper) confidence interval.
#' conf_limits_nc_chisq(chi_square = 30, alpha_lower = 0, alpha_upper = .05,
#'                      conf_level = NULL, df = 15)
#'
#' @keywords design multivariate regression
#'
#' @family noncentral distribution confidence limits
#'
#' @export
#' @importFrom stats pchisq uniroot

conf_limits_nc_chisq <- function(chi_square = NULL, conf_level = .95, df = NULL,
                                 alpha_lower = NULL, alpha_upper = NULL, tol = 1e-9,
                                 verbose = TRUE, ...) {
  if (is.null(chi_square)) stop("You must specify 'chi_square'.", call. = FALSE)
  if (chi_square < 0) stop("'chi_square' must be non-negative.", call. = FALSE)
  if (is.null(df)) stop("You must specify the degrees of freedom ('df').", call. = FALSE)

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
    alpha_lower <- alpha_upper <- (1 - conf_level) / 2
  }
  if (alpha_lower < 0 || alpha_upper < 0 || alpha_lower + alpha_upper >= 1) {
    stop("'alpha_lower' and 'alpha_upper' must be non-negative and sum to less than 1.", call. = FALSE)
  }

  # The mean of a noncentral chi square is df + lambda, so a bracket of order
  # max(chi_square, df) * 4 puts the root well inside; uniroot's extendInt
  # widens it if not.
  bracket_right <- max(chi_square * 4, df * 4, 100)

  central_lower_tail <- suppressWarnings(pchisq(chi_square, df = df, ncp = 0))

  if (alpha_lower == 0) {
    lower_ncp <- 0
    achieved_alpha_lower <- 0
  } else if (1 - central_lower_tail >= alpha_lower) {
    # chi_square is not significant at alpha_lower under H0; no lambda > 0 puts
    # as little as alpha_lower mass beyond chi_square, since the noncentral
    # chi square's upper-tail mass at chi_square increases with lambda.
    warning(
      "The observed chi_square is below the alpha_lower critical value of the central chi square distribution; the lower noncentrality limit has been clamped to 0 and the reported 'prob_greater' on the lower_limit row reflects the actual upper-tail probability at lambda = 0.",
      call. = FALSE
    )
    lower_ncp <- 0
    achieved_alpha_lower <- 1 - central_lower_tail
  } else {
    lower_ncp <- suppressWarnings(uniroot(
      function(candidate_ncp) suppressWarnings(pchisq(chi_square, df = df, ncp = candidate_ncp)) - (1 - alpha_lower),
      interval = c(0, bracket_right),
      tol = tol,
      extendInt = "downX",
      ...
    ))$root
    achieved_alpha_lower <- 1 - suppressWarnings(pchisq(chi_square, df = df, ncp = lower_ncp))
  }

  if (alpha_upper == 0) {
    upper_ncp <- Inf
    achieved_alpha_upper <- 0
  } else if (central_lower_tail <= alpha_upper) {
    warning(
      "The observed chi_square is so small that even at lambda = 0 the lower-tail probability is already at or below alpha_upper; the upper noncentrality limit is undefined and reported as NA.",
      call. = FALSE
    )
    upper_ncp <- NA_real_
    achieved_alpha_upper <- NA_real_
  } else {
    upper_ncp <- suppressWarnings(uniroot(
      function(candidate_ncp) suppressWarnings(pchisq(chi_square, df = df, ncp = candidate_ncp)) - alpha_upper,
      interval = c(0, bracket_right),
      tol = tol,
      extendInt = "downX",
      ...
    ))$root
    achieved_alpha_upper <- suppressWarnings(pchisq(chi_square, df = df, ncp = upper_ncp))
  }

  # prob_less and prob_greater report the noncentral chi square probability
  # mass on each side of the observed chi_square for the distribution centered
  # at that limit. By construction, prob_greater at the lower_limit equals
  # alpha_lower and prob_less at the upper_limit equals alpha_upper.
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
