#' Confidence Interval for the Noncentrality Parameter of a Noncentral \emph{t} Distribution
#'
#' Finds the noncentrality parameters of a noncentral \emph{t}-distribution
#' that bracket an observed \emph{t}-value with the requested tail
#' probabilities, giving a confidence interval on the population noncentrality
#' parameter. Together with \code{\link{ci_nc_F}} and
#' \code{\link{ci_nc_chisq}}, this is one of the low-level
#' noncentral distribution workhorses on which the \code{ci_*} confidence
#' interval functions (e.g., \code{\link{ci_smd}}, \code{\link{ci_smd_c}},
#' \code{\link{ci_cv}}) are built; most analyses reach it through those
#' functions rather than calling it directly.
#' The function corresponds to \code{conf.limits.nct()} in MBESS; its name
#' places it in the \code{ci_*} family because a confidence interval is what
#' it computes.
#'
#' @param ncp The observed \emph{t}-value, which is the estimate of the noncentrality parameter the interval is formed for
#' @param df The degrees of freedom
#' @param conf_level The level of confidence for a symmetric confidence interval. A level above .999999998, whose tail
#'   areas are below 1e-9, is refused (see "Small tail areas" in Details)
#' @param alpha_lower The proportion of values beyond the lower limit of the confidence interval (cannot be used with
#'   \code{conf_level}): 0, for an interval with no lower limit, or at least 1e-9 (see "Small tail areas" in Details)
#' @param alpha_upper The proportion of values beyond the upper limit of the confidence interval (cannot be used with
#'   \code{conf_level}): 0, for an interval with no upper limit, or at least 1e-9
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
#' \strong{Computation.} R's \code{\link[stats:TDist]{pt}()} sums a series for
#' the noncentral \emph{t} distribution up to a noncentrality of 37.62 and
#' with at most 400,000 degrees of freedom, and returns a normal
#' approximation otherwise. Compared with the exact probability over a grid
#' of degrees of freedom, noncentralities, and quantiles (from \eqn{10^{-6}}
#' to \eqn{1 - 10^{-6}}), it was within \eqn{10^{-9}} up to a noncentrality
#' of 37.62 with at most 2,000 degrees of freedom and up to 33 with more (the
#' series stops short in the upper tail near the top of its range: the error
#' was \eqn{10^{-6}} at 34 and .001 at 36 with 100,000 degrees of freedom),
#' while the normal approximation was off by as much as .14 with one degree
#' of freedom, .008 with 38, and .001 with 1,000; below one degree of
#' freedom it is not accurate at any noncentrality other than 0. Each limit
#' is therefore found with \code{pt()} where the limit lies in that range
#' (with at least one degree of freedom), and otherwise with the exact
#' probability
#' \deqn{P(T \le t) = \int_0^\infty \Phi(t s - \lambda)\, f(s)\, ds,}
#' where \eqn{\Phi} is the standard normal distribution function, \eqn{f} the
#' density of \eqn{S = \sqrt{V/\nu}} for \eqn{V} chi square on \eqn{\nu}
#' degrees of freedom, and \eqn{T = (Z + \lambda)/S} with \eqn{Z} standard
#' normal, independent of \eqn{V}; the integral is computed numerically. The
#' limits agreed with bisection on a second exact form, which conditions on
#' the normal variate instead, to within \eqn{3.5 \times 10^{-9}} of the
#' width of the interval in 400 random cases with 1 to 1,000,000 degrees of
#' freedom. The tail probabilities in \code{prob_less} and
#' \code{prob_greater} are computed the same way as the limit in their row.
#' A limit outside the range takes about a millisecond; one inside it, a
#' tenth of that.
#'
#' \strong{Small tail areas.} A nonzero tail area below 1e-9, or a
#' \code{conf_level} above .999999998, is refused with an error of class
#' \code{dmar_argument_invalid} whose \code{reason} is
#' \code{"tail_below_floor"}. A limit is the root of a tail probability,
#' and \code{pt()} sums its series to an absolute tolerance of
#' \eqn{10^{-12}} and forms one tail as 1 minus the other, so the error of
#' the probability becomes a larger part of the tail as the tail shrinks.
#' Against the exact probability, the limits with equal tails at
#' \eqn{t} = 5 on 47 degrees of freedom, .5 on 10, 2 on 200, 10 on 20, -3 on
#' 30, and 1 on 1,000, at 66 tail areas from \eqn{10^{-6}} down to
#' \eqn{10^{-9}}, were off by at most \eqn{1.7 \times 10^{-5}} at tail
#' areas of \eqn{10^{-8}} and above, \eqn{7.3 \times 10^{-5}} at
#' \eqn{2 \times 10^{-9}} and above, and \eqn{1.4 \times 10^{-4}} at
#' \eqn{10^{-9}} (the upper limit at \eqn{t} = 10 on 20 degrees of freedom,
#' \eqn{6.0 \times 10^{-6}} of the width of the interval); at
#' \eqn{10^{-10}} the largest error was \eqn{1.1 \times 10^{-3}}, and at
#' \eqn{10^{-20}} the upper limit at \eqn{t} = 5 on 47 degrees of freedom
#' was 7.6 percent too small. Every interval and plan whose limits come
#' from this function keeps the floor, among them \code{\link{ci_smd}},
#' \code{\link{ci_smd_c}}, \code{\link{ci_sm}}, \code{\link{ci_cv}},
#' \code{\link{ci_sc}}, \code{\link{ci_sc_ancova}}, \code{\link{ci_src}},
#' and \code{\link{ci_reg_coef}} with \code{noncentral = TRUE}, and the
#' planners of those intervals. The central intervals of a contrast, a
#' correlation, and a regression coefficient (\code{\link{ci_c}},
#' \code{\link{ci_r}}, \code{\link{ci_rc}}, and \code{\link{ci_reg_coef}}
#' with \code{noncentral = FALSE}) and their planners
#' (\code{\link{ss_aipe_c}}, \code{\link{ss_aipe_r}}, and
#' \code{\link{ss_aipe_rc}}) keep it too, so that an interval and the plan
#' for it accept the same tail areas.
#'
#' \strong{A limit past 37.62.} Each limit is found by its own rule, so an
#' interval can have its lower limit from \code{pt()} and its upper limit
#' from the exact probability. The interval for the standardized mean
#' difference of two groups of 20 at \eqn{d = 10} with 99 percent
#' confidence (\code{\link{ci_smd}}) is such an interval: it is
#' [7.003, 13.099], 6.096 wide. In 100,000 samples of two groups of 20 with
#' \eqn{\delta = 13}, this 99 percent interval fell below \eqn{\delta} in
#' .52 percent and above it in .48 percent (.50 percent each by design).
#'
#' This function is especially useful for forming confidence intervals around
#' standardized mean differences (Cohen's \emph{d}, Glass's \emph{g}, Hedges' \emph{g}),
#' standardized regression coefficients, and coefficients of variation.
#'
#' @return
#' A \code{data.frame} with one row per confidence limit and the columns:
#' \item{term}{Either \code{"ci_lower"} or \code{"ci_upper"}.}
#' \item{value}{The noncentrality parameter at that limit. \code{-Inf} when \code{alpha_lower = 0}; \code{Inf} when \code{alpha_upper = 0}.}
#' \item{prob_less}{(\code{verbose = TRUE}) The probability \eqn{P(T \le \mathrm{ncp})} that a \emph{t}-statistic from the noncentral \emph{t}-distribution centered at the row's limit falls at or below the observed \code{ncp}. By construction this equals \code{alpha_upper} on the \code{ci_upper} row and \eqn{1 - \mathtt{alpha\_lower}} on the \code{ci_lower} row.}
#' \item{prob_greater}{(\code{verbose = TRUE}) The complementary probability \eqn{P(T \ge \mathrm{ncp})}. By construction this equals \code{alpha_lower} on the \code{ci_lower} row and \eqn{1 - \mathtt{alpha\_upper}} on the \code{ci_upper} row.}
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
#'   estimation and the evaluation of statistical models. In L. L. Harlow,
#'   S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were no
#'   significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link[stats:TDist]{stats::pt()}}, \code{\link[stats:TDist]{stats::qt()}}, \code{\link[stats]{uniroot}}, \code{\link{ci_smd}}, \code{\link{ci_smd_c}}, \code{\link{ci_nc_F}}, \code{\link{ci_nc_chisq}}
#'
#' @examples
#'# Suppose observed t-value based on 'df'=126 is 2.83. Finding the lower
#'# and upper critical values for the population noncentrality parameter
#'# with a symmetric confidence interval with 95\% confidence is given as:
#'ci_nc_t(ncp = 2.83, df = 126, conf_level = .95)
#'
#'# Modifying the above example so that a nonsymmetric 95% confidence interval
#'# can be formed:
#'ci_nc_t(ncp = 2.83, df = 126, alpha_lower = .01, alpha_upper = .04, conf_level = NULL)
#'
#'# Modifying the above example so that a single-sided 95% confidence interval
#'# can be formed:
#'ci_nc_t(ncp = 2.83, df = 126, alpha_lower = 0, alpha_upper = .05, conf_level = NULL)
#'
#' @keywords models htest
#'
#' @family noncentral distribution confidence intervals
#'
#' @export
#' @importFrom stats pt uniroot

ci_nc_t <- function(ncp, df, conf_level = .95, alpha_lower = NULL, alpha_upper = NULL,
                            tol = 1e-9, verbose = TRUE, ...) {
  .dmar_check_conf_level(conf_level, null_ok = TRUE)
  .dmar_check_numbers("ci_nc_t", c("ncp", "df", "tol"), required = c("ncp", "df"))
  if (missing(ncp) || "t_value" %in% names(list(...))) {
    .dmar_refuse_argument(c(
      "ci_nc_t() takes the observed t value as 'ncp', and it cannot be ",
      "left out. 'ncp' is the one name for the value the interval is ",
      "formed for; 't_value' is not an argument. Write, for example, ",
      "ci_nc_t(ncp = 2.83, df = 126)."),
      arg = "ncp")
  }

  if (df <= 0) .dmar_refuse_argument(
    "'df', the degrees of freedom, must be positive.", arg = "df")

  alphas_supplied <- !is.null(alpha_lower) || !is.null(alpha_upper)
  if (!is.null(conf_level) && alphas_supplied) {
    .dmar_refuse_argument(
      "Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'; you cannot mix them.",
      arg = c("conf_level", "alpha_lower", "alpha_upper"))
  }
  if (is.null(conf_level)) {
    if (is.null(alpha_lower) || is.null(alpha_upper)) {
      .dmar_refuse_argument(
        "When 'conf_level' is NULL, both 'alpha_lower' and 'alpha_upper' must be specified.",
        arg = c("conf_level", "alpha_lower", "alpha_upper"))
    }
  } else {
    if (conf_level <= 0 || conf_level >= 1) .dmar_refuse_argument(
      "'conf_level' must be between 0 and 1.", arg = "conf_level")
    alpha_lower <- (1 - conf_level) / 2
    alpha_upper <- (1 - conf_level) / 2
  }

  if (alpha_lower < 0 || alpha_upper < 0 || alpha_lower + alpha_upper >= 1) {
    .dmar_refuse_argument(
      "'alpha_lower' and 'alpha_upper' must be non-negative and sum to less than 1.",
      arg = c("alpha_lower", "alpha_upper"))
  }
  # A nonzero tail below 1e-9 is refused (R/dmar_conditions.R): the tail
  # probability a limit is found from carries an absolute error of up to
  # about 1e-12 (pt() sums its series to that tolerance, and .nc_t_pt()
  # forms one tail as 1 minus the other), and below the floor the limit
  # loses accuracy. Every interval and plan whose limits come from this
  # function is held to it here.
  .dmar_check_tail_floor(alpha_lower, alpha_upper, conf_level = conf_level)

  # Initial bracket half-width: a generous multiple of the asymptotic SE of the
  # noncentrality estimator, sqrt(1 + ncp^2 / (2 * df)). uniroot's extendInt
  # widens the bracket if a sign change is not initially present.
  bracket_half_width <- max(10, 8 * sqrt(1 + ncp^2 / (2 * df)))

  # Each limit is found first with R's pt() (.nc_t_pt()), which is fast. A
  # limit inside the range where pt() is accurate (.nc_t_pt_limit()) by more
  # than twice the tolerance, so that a root uniroot() puts at the jump
  # where pt() changes method is not taken for one, is then exact. A limit
  # outside it is found again with the probability of .nc_t_prob(): pt()
  # inside the range and the exact integral of .nc_t_cdf() outside it. The
  # search starts at the first value with a Newton step that takes its slope
  # from pt() there and continues by secant steps, which settle in three or
  # four evaluations of the integral; should they not settle within 30
  # steps, uniroot() takes over from a bracket about the first value. The
  # tail probability at a limit (prob_less and prob_greater) comes from the
  # same computation as the limit. `lower_tail` names the tail of T at the
  # observed ncp: P(T > ncp) for the lower limit, which increases with the
  # noncentrality ("upX"), and P(T <= ncp) for the upper limit, which
  # decreases ("downX").
  inside_end <- .nc_t_pt_limit(df) - 2 * tol - 1e-6
  find_limit <- function(target, lower_tail, interval, extend) {
    first <- function(candidate_ncp) .nc_t_pt(ncp, df, candidate_ncp, lower_tail)
    root <- uniroot(function(candidate_ncp) first(candidate_ncp) - target,
                    interval = interval, tol = tol, extendInt = extend, ...)$root
    if (abs(root) < inside_end) return(c(root, first(root)))
    exact <- function(candidate_ncp)
      .nc_t_prob(ncp, df, candidate_ncp, lower_tail) - target
    h <- 1e-4 * max(1, abs(root))
    slope <- (first(root + h) - first(root - h)) / (2 * h)
    x_0 <- root
    f_0 <- exact(x_0)
    x_1 <- x_0 - f_0 / slope
    if (is.finite(x_1) && (slope > 0) == (extend == "upX")) {
      for (i in seq_len(30L)) {
        f_1 <- exact(x_1)
        if (f_1 == 0) return(c(x_1, target))
        if (f_1 == f_0) break
        step <- f_1 * (x_1 - x_0) / (f_1 - f_0)
        if (!is.finite(step)) break
        if (abs(step) <= tol) return(c(x_1, f_1 + target))
        x_0 <- x_1
        f_0 <- f_1
        x_1 <- x_1 - step
      }
    }
    step <- max(1, bracket_half_width / 20)
    root <- uniroot(exact, interval = root + c(-step, step), tol = tol,
                    extendInt = extend, ...)$root
    c(root, exact(root) + target)
  }

  if (alpha_lower == 0) {
    lower_ncp <- -Inf
    achieved_alpha_lower <- 0
  } else {
    found <- find_limit(alpha_lower, lower_tail = FALSE,
                        interval = c(ncp - bracket_half_width, ncp),
                        extend = "upX")
    lower_ncp <- found[1]
    achieved_alpha_lower <- found[2]
  }

  if (alpha_upper == 0) {
    upper_ncp <- Inf
    achieved_alpha_upper <- 0
  } else {
    found <- find_limit(alpha_upper, lower_tail = TRUE,
                        interval = c(ncp, ncp + bracket_half_width),
                        extend = "downX")
    upper_ncp <- found[1]
    achieved_alpha_upper <- found[2]
  }

  # prob_less and prob_greater report the noncentral t probability mass on each
  # side of the observed `ncp` for the distribution centered at that limit.
  # By construction, prob_greater at the ci_lower equals alpha_lower and
  # prob_less at the ci_upper equals alpha_upper.
  term <- c("ci_lower", "ci_upper")
  value <- c(lower_ncp, upper_ncp)
  prob_less <- c(1 - achieved_alpha_lower, achieved_alpha_upper)
  prob_greater <- c(achieved_alpha_lower, 1 - achieved_alpha_upper)

  if (verbose) {
    data.frame(term, value, prob_less, prob_greater)
  } else {
    data.frame(term, value)
  }
}
