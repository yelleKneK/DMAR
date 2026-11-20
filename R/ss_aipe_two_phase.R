#' Sample Size Planning From the AIPE Perspective for a Mean Under Two-Phase Sampling
#'
#' @description
#' Plans the number of labeled units, \eqn{n}, on which a costly target
#' measurement must be taken so that the confidence interval for the target
#' mean from \code{\link{ci_two_phase}} is sufficiently narrow, when an
#' inexpensive proxy measurement is available on every unit of a pool of
#' \eqn{M} units. The plan follows the accuracy in parameter estimation
#' (AIPE) approach: the goal is an interval no wider than \code{width},
#' either at the planning values (the expected width) or with a stated
#' assurance that the realized interval will be that narrow (Kelley &
#' Rausch, 2006). Beside the answer it reports the number of labeled units the
#' same goal needs without the proxy, so the saving is in view.
#'
#' @param true_sd_Y Planning value of the standard deviation of the target
#'   measurement, \eqn{\sigma_Y}.
#' @param true_cor Planning value of the correlation between the target and the
#'   proxy, \eqn{\rho}.
#' @param width The desired full width of the confidence interval (the
#'   upper limit minus the lower limit), in the units of the target.
#' @param n_pool The number of units in the pool, \eqn{M}, every one of
#'   which carries the proxy. The default \code{Inf} plans as though the pool
#'   were unlimited.
#' @param conf_level Confidence level (default \code{0.95}).
#' @param assurance \code{NULL} (default) to plan for the expected width, or
#'   a probability in (0, 1) that the realized width will be no wider than
#'   \code{width}, as in \code{\link{ss_aipe_c}}.
#' @param true_kurtosis Excess kurtosis (the scale of \code{\link{kurtosis}}) of
#'   the residual of the target on the proxy; default \code{0}, its normal
#'   value. Used only with \code{assurance}, where it adjusts the degrees of
#'   freedom of the variance estimate (see Details).
#'
#' @details
#' \strong{The interval planned for.} The regression estimator of double
#' sampling (Tamhane, 1978, Section 2.1, p. 420; the design is that of
#' Neyman, 1938, who formed strata from the inexpensive variable) that \code{\link{ci_two_phase}}
#' reports has, to first order, variance
#' \eqn{\sigma_Y^2(1 - \rho^2)/n + \rho^2\sigma_Y^2/M}, and its interval
#' uses the \emph{t} distribution on \eqn{n - 2} degrees of freedom. With
#' \eqn{\omega} the desired width (\code{width}), \eqn{1 - \alpha} the
#' \code{conf_level}, and \eqn{t_n = t_{1-\alpha/2,\, n-2}},
#' the plan for the expected width is the smallest \eqn{n} with
#' \deqn{2\, t_n \left\{\frac{\sigma_Y^2(1 - \rho^2)}{n} +
#'   \frac{\rho^2\sigma_Y^2}{M}\right\}^{1/2} \le \omega.}
#' In an evaluation of a system whose outputs an automated judge labels,
#' for example, \eqn{n} is the number of outputs a human assessor must
#' label when the judge has labeled \eqn{M} of them.
#'
#' \strong{Assurance.} The realized width replaces \eqn{\sigma_Y^2(1 -
#' \rho^2)} by the residual mean square \eqn{S_e^2} and \eqn{\rho^2
#' \sigma_Y^2} by \eqn{\hat b^2 S_X^2}. Treating the second, which is
#' divided by the large \eqn{M}, as fixed at its planning value, and
#' \eqn{\nu S_e^2 / \sigma_e^2} as chi square on \eqn{\nu} degrees of
#' freedom, the probability that the realized width is at most \eqn{\omega}
#' is
#' \deqn{F_{\chi^2_\nu}\left(\frac{\nu\, n\, \{\omega^2 / (4 t_n^2) -
#'   \rho^2 \sigma_Y^2 / M\}}{\sigma_Y^2 (1 - \rho^2)}\right),}
#' and the plan is the smallest \eqn{n} from which it is at least
#' \code{assurance} at every larger \eqn{n}. Under normal residuals
#' \eqn{\nu = n - 2}. For residuals with excess kurtosis \eqn{\kappa}
#' (\code{true_kurtosis}) the
#' variance of \eqn{S_e^2} is approximately
#' \eqn{\sigma_e^4\{2/(n - 2) + \kappa/n\}}, and matching it to a scaled chi
#' square gives \eqn{\nu = 2 / \{2/(n - 2) + \kappa/n\}}: for an assurance
#' of .80 or more, a lighter-tailed residual (\eqn{\kappa < 0}) needs fewer
#' units and a heavier-tailed one more. A binary target is not light-tailed
#' as a rule: a 0/1 variable with proportion \eqn{p} has excess kurtosis
#' \eqn{\{1 - 6p(1 - p)\}/\{p(1 - p)\}}, which is negative only for \eqn{p}
#' between about .21 and .79 and is about 5.1 at \eqn{p = .90}. The
#' planning value is best taken from the residuals of a pilot regression of
#' the target on the proxy, through \code{\link{kurtosis}}. With a large
#' \eqn{\kappa} and a small \eqn{n}, \eqn{\nu} is small, and the two-moment
#' approximation can reach the assurance at a few small \eqn{n}, fall
#' below it over a stretch of larger ones, and reach it again; the plan is
#' the start of the stretch over which it holds for good, which the search
#' finds exactly for any \code{assurance} in (0, 1). Under normal
#' residuals and with an assurance of one half or more the two readings
#' agree.
#'
#' \strong{The pool sets a floor.} The labeled units are drawn from the
#' pool, and \code{\link{ci_two_phase}} needs at least one unit that
#' carries the proxy alone, so \eqn{n \le M - 1}. No number of labeled
#' units can therefore bring the width below
#' \eqn{2 z_{1-\alpha/2} \{\sigma_Y^2(1 - \rho^2)/(M - 1) +
#' \rho^2\sigma_Y^2/M\}^{1/2}} (the \code{minimum_width} row), which is
#' close to \eqn{2 z_{1-\alpha/2}\, \sigma_Y / \sqrt{M}}, the width of the
#' target mean were every unit of the pool labeled; a target at or below
#' it is an error. A target just above it can still be out of reach once
#' the \emph{t} quantile and the assurance are accounted for; when the plan
#' needs more than \eqn{M - 1} labeled units, that too is an error, naming
#' the pool as the constraint, rather than a missing value.
#'
#' \strong{The target-only comparison.} \code{necessary_n_target_only} is
#' the planned size of the usual \emph{t} interval for the target mean from
#' labeled units alone, the plan \code{ss_aipe_mean()} makes for a single
#' mean, so the two functions agree. Without an assurance it is the
#' smallest \eqn{n} with
#' \eqn{2\, t_{1-\alpha/2,\, n-1}\, c_4(n)\, \sigma_Y / \sqrt{n} \le
#' \omega}, where \eqn{c_4(n)\sigma_Y} is the expected sample standard
#' deviation under normality, so that the left side is the exact expected
#' width; with \code{assurance} it is the same chi square argument on
#' \eqn{n - 1} degrees of freedom. It is not bounded by \code{n_pool},
#' because a design without the proxy has no use for a pool. The
#' \code{true_kurtosis} supplied is the residual's, and the comparison uses it
#' for the target as well, an approximation: the target's own kurtosis can
#' differ, and the comparison for it is \code{ss_aipe_mean()} with that
#' value. The labeled-unit plan, by contrast, sets the two estimates in its
#' standard error at their planning values, since the expected value of
#' the two-phase standard error has no closed form once the pool term
#' enters; this overstates its expected width slightly, by a factor of at
#' most about \eqn{1 + 1/\{4(n - 2)\}}, and can add a unit to
#' \code{necessary_n_labeled}.
#'
#' @return A \code{data.frame} (classes \code{dmar_ss_aipe} and
#'   \code{dmar_tbl}) with rows
#'   \item{necessary_n_labeled}{the planned number of labeled units.}
#'   \item{necessary_n_target_only}{the number the same goal needs without
#'     the proxy.}
#'   \item{ratio}{the first divided by the second.}
#'   \item{expected_width}{the full width at \code{necessary_n_labeled} and
#'     the planning values.}
#'   \item{minimum_width}{the narrowest width the pool allows, the
#'     normal-quantile width with all but one unit of the pool labeled
#'     (0 for an unlimited pool; see Details).}
#'   \item{width, true_sd_Y, true_cor, n_pool, conf_level}{the inputs.}
#'   \item{assurance, true_kurtosis}{the inputs, present only when an assurance
#'     was supplied.}
#'
#' @references
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
#'   standardized mean difference: Accuracy in parameter estimation via
#'   narrow confidence intervals. \emph{Psychological Methods, 11}(4),
#'   363--385. \doi{10.1037/1082-989X.11.4.363}
#'
#' Neyman, J. (1938). Contribution to the theory of sampling human
#'   populations. \emph{Journal of the American Statistical Association,
#'   33}(201), 101--116. \doi{10.1080/01621459.1938.10503378}
#'
#' Tamhane, A. C. (1978). Inference based on regression estimator in
#'   double sampling. \emph{Biometrika, 65}(2), 419--427.
#'   \doi{10.1093/biomet/65.2.419}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{ci_two_phase}} for the interval,
#'   \code{\link{ss_aipe_two_phase_cost}} for the least costly numbers of
#'   labeled and pool units, \code{\link{ss_aipe_two_phase_sensitivity}} to
#'   evaluate a plan by simulation, and \code{\link{ss_seq_two_phase}} for
#'   the sequential rule.
#'
#' @family AIPE sample size planning
#'
#' @keywords design
#'
#' @examples
#' # Planning values from a pilot: a target standard deviation of 0.43 and
#' # a correlation of 0.55 with the proxy, which is available on a pool of
#' # 20,000 units. A 90% interval no wider than 0.10:
#' ss_aipe_two_phase(true_sd_Y = 0.43, true_cor = 0.55, width = 0.10,
#'                   n_pool = 20000,
#'                   conf_level = 0.90)
#'
#' # The same goal with 80% assurance that the realized interval is no
#' # wider than 0.10:
#' ss_aipe_two_phase(true_sd_Y = 0.43, true_cor = 0.55, width = 0.10,
#'                   n_pool = 20000,
#'                   conf_level = 0.90, assurance = 0.80)
#'
#' # A smaller pool raises the floor and the number of labeled units.
#' ss_aipe_two_phase(true_sd_Y = 0.43, true_cor = 0.55, width = 0.10,
#'                   n_pool = 1000,
#'                   conf_level = 0.90, assurance = 0.80)
#'
#' @export
ss_aipe_two_phase <- function(true_sd_Y, true_cor, width, n_pool = Inf,
                              conf_level = 0.95, assurance = NULL,
                              true_kurtosis = 0) {
  .dmar_check_conf_level(conf_level)
  if (!missing(width)) .dmar_check_width(width)
  .two_phase_check_planning(sd_y = true_sd_Y, rho = true_cor, width = width,
                            conf_level = conf_level, assurance = assurance,
                            true_kurtosis = true_kurtosis)
  if (!is.numeric(n_pool) || length(n_pool) != 1L || is.na(n_pool) ||
      n_pool < 4 || (is.finite(n_pool) && n_pool != round(n_pool)))
    .dmar_refuse_argument(
      "'n_pool' must be Inf or a single whole number of at least 4.",
      arg = "n_pool")
  if (is.null(assurance) && true_kurtosis != 0)
    .dmar_warn(c(
      "'true_kurtosis' affects only a plan with an 'assurance'; it is ",
      "ignored when 'assurance' is NULL."), "argument_ignored")

  half <- width / 2
  q <- 1 - (1 - conf_level) / 2
  z <- stats::qnorm(q)
  var_error <- true_sd_Y^2 * (1 - true_cor^2)
  var_fixed <- if (is.finite(n_pool)) true_cor^2 * true_sd_Y^2 / n_pool else 0
  # ci_two_phase() needs at least one unit with the proxy alone, so at most
  # n_pool - 1 units can be labeled, and the floor is the normal-quantile
  # width there.
  n_max <- n_pool - 1
  minimum_width <- if (is.finite(n_pool))
    2 * z * sqrt(var_error / n_max + var_fixed) else 0
  if (width <= minimum_width)
    .dmar_refuse_argument(c(
      "The pool is too small for this width: with n_pool = ", n_pool,
      " no number of labeled units gives a width below ",
      signif(minimum_width, 4), " (the width when every unit of the ",
      "pool but one is labeled)."), arg = NULL)

  # The smallest n from which the goal holds at every larger n, for an
  # interval whose variance estimate spends df_lost degrees of freedom
  # (2 for the regression, 1 for the mean alone), whose variance is
  # var_err / n + var_fix, and whose expected width at n is
  # expected_width(n), the criterion when no assurance is given. With an
  # assurance the realized width is at most `width` exactly when
  # S^2 / sigma^2 <= bound(n), and the criterion is the scaled chi square
  # probability of that event; see .two_phase_search() for the
  # certificates. bound(n) and nu(n) both rise with n (nu for any kurtosis
  # above -2, the derivative of its denominator being
  # -2 / (n - df_lost)^2 - true_kurtosis / n^2 < 0), and qchisq(assurance, nu)
  # rises with nu, which gives holds_on(); Cantelli's inequality gives
  # qchisq(assurance, nu) / nu <= 1 + sqrt(2 assurance / ((1 - assurance) nu)),
  # a bound that falls as nu rises, which gives holds_from(). The factor tol
  # keeps rounding in qt() and qchisq() from certifying an n at the edge of
  # the criterion.
  plan_n <- function(var_err, var_fix, df_lost, expected_width) {
    n_min <- df_lost + 1
    start <- var_err / ((half / z)^2 - var_fix)
    if (is.null(assurance)) {
      # The expected width falls as n rises, so meeting the goal at n
      # certifies every larger n.
      meets <- function(n) expected_width(n) <= width
      return(.two_phase_search(meets, meets, function(a, b) meets(a),
                               start, n_min))
    }
    bound <- function(n)
      n * ((half / stats::qt(q, n - df_lost))^2 - var_fix) / var_err
    nu <- function(n) 2 / (2 / (n - df_lost) + true_kurtosis / n)
    meets <- function(n) {
      b <- bound(n)
      b > 0 && stats::pchisq(nu(n) * b, nu(n)) >= assurance
    }
    tol <- 1 + 1e-12
    holds_from <- function(n) {
      bound(n) >= tol * (1 + sqrt(2 * assurance / ((1 - assurance) * nu(n))))
    }
    holds_on <- function(a, b) {
      nu(a) * bound(a) >= tol * stats::qchisq(assurance, nu(b))
    }
    .two_phase_search(meets, holds_from, holds_on, start, n_min)
  }

  # The labeled units: the estimates in the standard error at their
  # planning values, on the n - 2 degrees of freedom of the residual mean
  # square.
  n_labeled <- plan_n(var_error, var_fixed, 2, function(n)
    2 * stats::qt(q, n - 2) * sqrt(var_error / n + var_fixed))
  if (n_labeled > n_max)
    .dmar_refuse_argument(c(
      "No number of labeled units up to n_pool - 1 = ", n_max,
      " (which leaves one unit of the pool with the proxy alone) attains ",
      "this width",
      if (!is.null(assurance)) " with this assurance",
      ": the requirement passes the pool once the t quantile",
      if (!is.null(assurance)) " and the assurance are" else " is",
      " taken into account. Enlarge the pool or widen the target."), arg = NULL)
  # The target alone: the t interval for one mean, whose expected width
  # under normality is exact through c4(n), the plan ss_aipe_mean() makes;
  # both use the one .c4() so the two give the same answer.
  n_target_only <- plan_n(true_sd_Y^2, 0, 1, function(n)
    2 * stats::qt(q, n - 1) * .c4(n) * true_sd_Y / sqrt(n))
  expected_width <- 2 * stats::qt(q, n_labeled - 2) *
    sqrt(var_error / n_labeled + var_fixed)

  out <- data.frame(
    term  = c("necessary_n_labeled", "necessary_n_target_only", "ratio",
              "expected_width", "minimum_width", "width", "true_sd_Y",
              "true_cor", "n_pool", "conf_level",
              if (!is.null(assurance)) c("assurance", "true_kurtosis")),
    value = c(n_labeled, n_target_only, n_labeled / n_target_only,
              expected_width, minimum_width, width, true_sd_Y,
              true_cor, n_pool, conf_level,
              if (!is.null(assurance)) c(assurance, true_kurtosis)),
    stringsAsFactors = FALSE, row.names = NULL
  )
  .as_dmar_tbl(out, conf_level = conf_level, subclass = "dmar_ss_aipe")
}
