#' Confidence Limits for the Noncentrality Parameter of a Noncentral \emph{F}-distribution
#'
#' Finds the noncentrality parameters of a noncentral \emph{F}-distribution
#' that bracket an observed \emph{F}-value with the requested tail
#' probabilities, giving a confidence interval on the population noncentrality
#' parameter. Together with \code{\link{conf_limits_nct}} and
#' \code{\link{conf_limits_nc_chisq}}, this is one of the low-level
#' noncentral distribution workhorses on which the \code{ci_*} confidence
#' interval functions (e.g., \code{\link{ci_pvaf}}, \code{\link{ci_snr}},
#' \code{\link{ci_R2}}) are built; most analyses reach it through those
#' functions rather than calling it directly.
#'
#' @param F_value The observed \emph{F}-value
#' @param conf_level The desired degree of confidence for a symmetric interval
#' @param df_1 The numerator degrees of freedom
#' @param df_2 The denominator degrees of freedom
#' @param alpha_lower The proportion of values beyond the lower limit (cannot be used with \code{conf_level})
#' @param alpha_upper The proportion of values beyond the upper limit (cannot be used with \code{conf_level})
#' @param tol The convergence tolerance passed to \code{\link[stats]{uniroot}}
#' @param verbose If \code{TRUE} (the default), the returned data frame additionally reports the achieved tail probabilities at each limit; if \code{FALSE}, only \code{term} and \code{value} are returned
#' @param \dots Additional arguments forwarded to \code{\link[stats]{uniroot}}
#'
#' @details
#' Each confidence limit is the noncentrality parameter \eqn{\lambda \ge 0} of a
#' noncentral \emph{F}-distribution with \code{df_1} and \code{df_2} degrees of
#' freedom whose appropriate tail at the observed \code{F_value} contains the
#' requested probability:
#' \itemize{
#'   \item the lower limit satisfies \eqn{P(F \ge \mathtt{F\_value}) = \mathtt{alpha\_lower}};
#'   \item the upper limit satisfies \eqn{P(F \le \mathtt{F\_value}) = \mathtt{alpha\_upper}}.
#' }
#' The two conditions run in opposite directions in \eqn{\lambda}: the
#' lower-tail probability \eqn{P(F \le \mathtt{F\_value})} is continuous and
#' strictly decreasing in the noncentrality parameter, so the upper-tail
#' probability \eqn{P(F \ge \mathtt{F\_value})} is continuous and strictly
#' increasing in it. The lower limit is the \eqn{\lambda} at which the upper
#' tail has grown to \code{alpha_lower}, and the upper limit is the
#' \eqn{\lambda} at which the lower tail has shrunk to \code{alpha_upper}.
#' Each is therefore the unique non-negative root of a one-dimensional
#' equation, and both are located with \code{\link[stats]{uniroot}} on the
#' decreasing lower-tail scale; \code{extendInt} is used to widen the search
#' bracket if needed.
#'
#' Because the noncentrality parameter is bounded below by zero, the lower limit
#' is set to zero whenever the observed \code{F_value} is smaller than the
#' \code{alpha_lower} critical value of the central \emph{F}-distribution (i.e.,
#' the data is consistent with \eqn{\lambda = 0} at the requested confidence
#' level). A warning is issued in that case, and the achieved probabilities
#' reported in the output reflect the actual values at \eqn{\lambda = 0} rather
#' than the requested \code{alpha_lower}. The warning carries the condition
#' class \code{dmar_ncf_clamp}, so a caller that inverts the noncentral
#' \emph{F} repeatedly can muffle or deduplicate it by class.
#'
#' @return
#' A \code{data.frame} with one row per confidence limit and the columns:
#' \item{term}{Either \code{"lower_limit"} or \code{"upper_limit"}.}
#' \item{value}{The noncentrality parameter at that limit. \code{0} when \code{alpha_lower = 0} or the lower limit is unattainable; \code{Inf} when \code{alpha_upper = 0}; \code{NA} on the \code{upper_limit} row when the observed \code{F_value} is so small that even at \eqn{\lambda = 0} the lower-tail probability is already at or below \code{alpha_upper}, leaving the upper noncentrality limit undefined (a warning is issued).}
#' \item{prob_less}{(\code{verbose = TRUE}) The probability \eqn{P(F \le \mathtt{F\_value})} that an \emph{F}-statistic from the noncentral \emph{F}-distribution centered at the row's limit falls at or below the observed \code{F_value}.}
#' \item{prob_greater}{(\code{verbose = TRUE}) The complementary probability \eqn{P(F \ge \mathtt{F\_value})}. By construction this equals \code{alpha_lower} on the \code{lower_limit} row and \eqn{1 - \mathtt{alpha\_upper}} on the \code{upper_limit} row.}
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
#' \code{\link{ss_aipe_R2}}, \code{\link{ci_R2}}, \code{\link{conf_limits_nct}}, \code{\link{conf_limits_nc_chisq}}, \code{\link[stats:FDist]{stats::pf()}}, \code{\link[stats:FDist]{stats::qf()}}, \code{\link[stats]{uniroot}}
#'
#' @examples
#' conf_limits_ncf(F_value = 5, conf_level = .95, df_1 = 5, df_2 = 100)
#'
#' # A one-sided (upper) confidence interval.
#' conf_limits_ncf(F_value = 5, conf_level = NULL, df_1 = 5, df_2 = 100,
#'                 alpha_lower = 0, alpha_upper = .05)
#'
#' @keywords design multivariate regression
#'
#' @family noncentral distribution confidence limits
#'
#' @export
#' @importFrom stats pf uniroot

conf_limits_ncf <- function(F_value = NULL, conf_level = .95, df_1 = NULL, df_2 = NULL,
                            alpha_lower = NULL, alpha_upper = NULL, tol = 1e-9,
                            verbose = TRUE, ...) {
  if (is.null(F_value)) stop("You must specify 'F_value'.", call. = FALSE)
  if (F_value < 0) stop("'F_value' must be non-negative.", call. = FALSE)
  if (is.null(df_1) || is.null(df_2)) {
    stop("You must specify both numerator and denominator degrees of freedom ('df_1' and 'df_2').", call. = FALSE)
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
    alpha_lower <- alpha_upper <- (1 - conf_level) / 2
  }
  if (alpha_lower < 0 || alpha_upper < 0 || alpha_lower + alpha_upper >= 1) {
    stop("'alpha_lower' and 'alpha_upper' must be non-negative and sum to less than 1.", call. = FALSE)
  }

  # The mean of a noncentral F is roughly (df_1 + lambda) / df_1 * df_2 / (df_2 - 2).
  # An initial right bracket scaled by F_value * df_1 puts the root well inside,
  # and uniroot's extendInt widens it if not.
  bracket_right <- max(F_value * df_1 * 4, 100)

  central_lower_tail <- suppressWarnings(pf(F_value, df1 = df_1, df2 = df_2, ncp = 0))

  if (alpha_lower == 0) {
    lower_ncp <- 0
    achieved_alpha_lower <- 0
  } else if (1 - central_lower_tail >= alpha_lower) {
    # F_value is not significant at alpha_lower under H0; no lambda > 0 puts as
    # little as alpha_lower mass beyond F_value, since the noncentral F's
    # upper-tail mass at F_value increases with lambda from its already-too-
    # large value at lambda = 0. This is a normal consequence of a small
    # observed F, not a failure, so it is a warning, not an error. The
    # condition class "dmar_ncf_clamp" lets the ci_* callers restate the
    # consequence for their own effect size and lets the iterative planners
    # deduplicate the warning by class.
    warning(warningCondition(
      "The observed F_value is below the alpha_lower critical value of the central F-distribution, so the lower confidence limit on the noncentrality parameter is 0.",
      class = "dmar_ncf_clamp"
    ))
    lower_ncp <- 0
    achieved_alpha_lower <- 1 - central_lower_tail
  } else {
    lower_ncp <- tryCatch(
      suppressWarnings(uniroot(
        function(candidate_ncp) suppressWarnings(pf(F_value, df1 = df_1, df2 = df_2, ncp = candidate_ncp)) - (1 - alpha_lower),
        interval = c(0, bracket_right),
        tol = tol,
        extendInt = "downX",
        ...
      ))$root,
      error = function(e) {
        stop(sprintf(
          "In conf_limits_ncf(), the root search for the lower noncentrality limit failed for F_value = %s with df_1 = %s and df_2 = %s (alpha_lower = %s). The inner solver reported: %s. Verify that 'F_value' and the degrees of freedom describe the intended analysis; if they do, adjust 'tol' or request one-sided limits through 'alpha_lower' and 'alpha_upper'.",
          format(F_value), format(df_1), format(df_2), format(alpha_lower), conditionMessage(e)
        ), call. = FALSE)
      }
    )
    achieved_alpha_lower <- 1 - suppressWarnings(pf(F_value, df1 = df_1, df2 = df_2, ncp = lower_ncp))
  }

  if (alpha_upper == 0) {
    upper_ncp <- Inf
    achieved_alpha_upper <- 0
  } else if (central_lower_tail <= alpha_upper) {
    warning(
      "The observed F_value is so small that even at lambda = 0 the lower-tail probability is already at or below alpha_upper; the upper noncentrality limit is undefined and reported as NA.",
      call. = FALSE
    )
    upper_ncp <- NA_real_
    achieved_alpha_upper <- NA_real_
  } else {
    upper_ncp <- tryCatch(
      suppressWarnings(uniroot(
        function(candidate_ncp) suppressWarnings(pf(F_value, df1 = df_1, df2 = df_2, ncp = candidate_ncp)) - alpha_upper,
        interval = c(0, bracket_right),
        tol = tol,
        extendInt = "downX",
        ...
      ))$root,
      error = function(e) {
        stop(sprintf(
          "In conf_limits_ncf(), the root search for the upper noncentrality limit failed for F_value = %s with df_1 = %s and df_2 = %s (alpha_upper = %s). The inner solver reported: %s. Verify that 'F_value' and the degrees of freedom describe the intended analysis; if they do, adjust 'tol' or request one-sided limits through 'alpha_lower' and 'alpha_upper'.",
          format(F_value), format(df_1), format(df_2), format(alpha_upper), conditionMessage(e)
        ), call. = FALSE)
      }
    )
    achieved_alpha_upper <- suppressWarnings(pf(F_value, df1 = df_1, df2 = df_2, ncp = upper_ncp))
  }

  # prob_less and prob_greater report the noncentral F probability mass on each
  # side of the observed F_value for the distribution centered at that limit.
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
