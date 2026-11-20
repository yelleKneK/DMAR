#' The Bryant-Paulson Generalized Studentized Range Distribution
#'
#' @description
#' Distribution function (\code{pbryant_paulson}), quantile/critical-value
#' function (\code{qbryant_paulson}), and density (\code{dbryant_paulson})
#' for the Bryant-Paulson generalized studentized range, the sampling
#' distribution of the studentized range of covariate-\emph{adjusted} means
#' in the analysis of covariance (ANCOVA) when the covariate(s) are
#' \emph{random}. These are the analysis-of-covariance analogs of
#' \code{\link[stats]{ptukey}} / \code{\link[stats]{qtukey}} and supply the
#' critical values needed for Tukey-Kramer-type simultaneous confidence
#' intervals on (and tests of) contrasts of adjusted means.
#'
#' @param q Vector of quantiles (values of the generalized studentized
#'   range statistic).
#' @param prob Vector of probabilities. For \code{qbryant_paulson} this is
#'   the cumulative probability (e.g., \code{0.95} returns the upper 5\%
#'   critical value).
#' @param n_covariates The number of random covariates, \eqn{p}
#'   (\eqn{p \ge 0}). With \code{n_covariates = 0} the distribution is
#'   exactly the ordinary studentized range, computed as in
#'   \code{\link{cv_tukey_hsd}} (\code{\link[stats]{ptukey}} and
#'   \code{\link[stats]{qtukey}} where those are accurate, the exact
#'   distribution where they are not).
#' @param n_groups The number of groups (treatments) being compared,
#'   \eqn{k} (\eqn{k \ge 2}); this is the \dQuote{sample size} of the range.
#' @param df The error degrees of freedom of the ANCOVA model, \eqn{\nu}.
#'   In a one-way ANCOVA with \eqn{N} total observations, \eqn{k} groups,
#'   and \eqn{p} covariates, \eqn{\nu = N - k - p}. \code{Inf} gives the
#'   limiting distribution (see Details).
#' @param lower_tail Logical; if \code{TRUE} (default) probabilities are
#'   \eqn{P(Q \le q)}, otherwise \eqn{P(Q > q)}.
#'
#' @details
#' \strong{The statistic.} In a balanced ANCOVA with \eqn{k} groups and
#' \eqn{p} random covariates, let \eqn{\hat\theta_i} be the adjusted group
#' means and \eqn{\hat\sigma_{y \mid x}} the square root of the ANCOVA error
#' mean square (on \eqn{\nu} degrees of freedom). The Bryant-Paulson
#' statistic is the studentized range of the adjusted means,
#' \deqn{Q \;=\; \frac{\max_i \hat\theta_i - \min_i \hat\theta_i}{\hat\sigma_{y\mid x}\sqrt{K_1 - K_2}},}
#' where \eqn{K_1 - K_2} is the design constant that scales the variance of a
#' single adjusted mean (for a one-way design with \eqn{n} per group,
#' \eqn{K_1 - K_2 = 1/n}). The studentizer uses only this
#' \dQuote{between-only} standard error: the extra sampling variability
#' induced by having to \emph{estimate} the covariate adjustment from random
#' covariates is carried by the distribution of \eqn{Q} itself, not by a
#' per-comparison standard-error correction. This is what distinguishes the
#' procedure from naively applying Tukey's method to adjusted means.
#'
#' \strong{The distribution.} Bryant and Paulson (1976) give the exact CDF of
#' \eqn{Q_p} in their Equation (17), a single integral over a variable that
#' combines the \eqn{\chi^2_\nu} error estimate with a random
#' covariate-shrinkage factor \eqn{\delta} that, by their Equations (11)--(12),
#' has a \eqn{\mathrm{Beta}((\nu+1)/2,\, p/2)} distribution. Carrying out the
#' error integral with the studentized-range routine \code{ptukey} reduces
#' Equation (17) to the equivalent one-dimensional form
#' \deqn{P(Q_p \le q) \;=\; \int_0^1 \mathrm{ptukey}\!\left(q\sqrt{\delta};\, k,\, \nu\right)\, f_{\mathrm{Beta}}\!\left(\delta;\, \tfrac{\nu+1}{2},\, \tfrac{p}{2}\right) d\delta,}
#' which this package evaluates (the reduction is exact; see the source-code
#' comments in \file{R/bryant_paulson.R} for the one-line derivation from
#' Bryant and Paulson's p. 634 conditioning argument). When \eqn{p = 0} the
#' factor \eqn{\delta} degenerates at 1 and \eqn{Q_p} is exactly the ordinary
#' studentized range (Bryant and Paulson, 1976, Sec. 1), so the code
#' short-circuits to the studentized range distribution function. The integral
#' is evaluated with
#' \code{\link[stats]{integrate}} after the change of variables
#' \eqn{\delta = 1 - u^2}, which removes the endpoint singularity of the Beta
#' weight at \eqn{\delta = 1} when \eqn{p = 1} and makes the quadrature
#' converge in a few subdivisions; \code{qbryant_paulson} inverts it with
#' \code{\link[stats]{uniroot}}. Bryant and Bruvold (1980) later showed the
#' same distribution and critical values remain valid when the covariates are
#' \emph{not} identically distributed across groups (their grouped-covariate
#' model, Eq. 1.3), and added the Duncan multiple-range extension.
#'
#' \strong{Accuracy at small df.} The studentized range distribution function
#' inside the integral is \code{\link[stats]{ptukey}} where that is accurate and
#' is evaluated directly elsewhere, without \code{ptukey}, by integrating the
#' probability integral of the range against the \eqn{\chi^2_\nu} error
#' density. \code{ptukey} loses accuracy at small \eqn{\nu}, increasingly so as
#' \eqn{k} grows (at \eqn{\nu = 3}, \eqn{k = 20} its probability error reaches
#' \eqn{\approx 3\times10^{-4}}, enough to move the critical value by about
#' 0.2, and it is larger at \eqn{\nu = 2}), and again above 10000 df; the
#' direct evaluation is used for fewer than 9 error df with up to 6 groups,
#' fewer than 11 with up to 10, fewer than 15 with up to 20, fewer than 21 with
#' up to 50, fewer than 31 with up to 100, fewer than 41 with up to 250, and,
#' with more groups, below a bound that rises by steps to 111 with more than
#' 3000 groups; above 10000 and up to \eqn{10^8} df; and, with more than 50 groups,
#' wherever the probability is below .75 (.80 with more than 1000 groups;
#' \code{\link{cv_tukey_hsd}} gives the
#' whole rule and how it was measured). At probabilities above .999 or below
#' .001 the quantile also comes from the exact distribution, with or without
#' covariates (with covariates every node of the integral over the
#' covariate-shrinkage weight then uses the exact studentized range: through
#' \code{ptukey} the quantile at .9999 with one covariate, 10 groups, and 11
#' df was 12.977072 instead of 12.977142), and with no covariates also
#' wherever \code{qtukey} stops short of its root. The direct path costs a
#' fraction of a second per probability.
#'
#' \strong{Large df.} As \eqn{\nu} grows, \eqn{\delta} tends to 1 and
#' \eqn{Q_p} to the studentized range with infinite df, whatever the number of
#' covariates; with \code{df = Inf} the functions return that limit, computed
#' as for \code{n_covariates = 0}. At large \eqn{\nu} the weight of the
#' integral lies within a few multiples of \eqn{\sqrt{p/\nu}} of
#' \eqn{u = 0}, where \eqn{\delta = 1 - u^2}, so the integral over \eqn{u}
#' stops at the square root of the upper \eqn{10^{-15}} quantile of the
#' \eqn{\mathrm{Beta}(p/2,\, (\nu+1)/2)} distribution of \eqn{u^2} rather
#' than at 1, which leaves out weight \eqn{10^{-15}}. Over the whole unit
#' interval the quadrature missed the weight from about \eqn{4 \times 10^7}
#' df with one covariate (\eqn{5.6 \times 10^7} with six), and the
#' distribution function came out near 0. The weight is computed on the log
#' scale, since with many covariates its factors overflow and underflow
#' separately. At \eqn{10^8} df the distribution function is still below its
#' infinite-df limit by about \eqn{10^{-8}} (\eqn{7 \times 10^{-9}} at 3.7
#' with one covariate and four groups).
#'
#' \strong{Validation.} The implementation reproduces Bryant and Paulson's
#' (1976) Table 1 exactly, to the two decimal places tabled, over the whole of
#' its range: both tail areas (\eqn{\alpha = .05} and \eqn{\alpha = .01}), all
#' three covariate counts (\eqn{p = 1, 2, 3}), every tabled number of groups
#' (\eqn{k = 2, \ldots, 8, 10, 12, 16, 20}), and every tabled error degrees of
#' freedom (\eqn{\nu = 2, \ldots, 8, 10, 12, 14, 16, 18, 20, 24, 30, 40, 60,
#' 120}), which is 1188 critical values in all. The corners of the table are
#' included: \eqn{q_{.01;\,1,2,2} = 19.09}, \eqn{q_{.01;\,3,20,2} = 73.01},
#' \eqn{q_{.01;\,3,20,3} = 33.13}, and \eqn{q_{.05;\,1,6,14} = 4.83} (the value
#' used in the Bryant and Bruvold, 1980 worked example). Two of the 1188
#' entries, \eqn{q_{.01;\,2,8,3}} and \eqn{q_{.01;\,2,20,4}}, have exact values
#' of 23.165013 and 19.745008, each roughly \eqn{10^{-5}} above the 23.165 and
#' 19.745 half-way points, so they round to 23.17 and 19.75; the 1976 table
#' rounds them down, to 23.16 and 19.74. Both values were confirmed to fourteen
#' significant figures by two independent high-order quadrature engines that
#' share no code with the package implementation. A large-scale simulation of
#' the Bryant and Paulson statistic confirms the computed values independently,
#' and the Bryant and Bruvold (1980) Table 2 Duncan ranges are reproduced as
#' well. See the package tests.
#'
#' @return
#' Numeric vectors. \code{pbryant_paulson} returns cumulative (or upper-tail)
#' probabilities, \code{dbryant_paulson} returns density values, and
#' \code{qbryant_paulson} returns critical values (quantiles) of the
#' generalized studentized range. Results are recycled to the length of the
#' longest of \code{q}/\code{prob} and the parameter arguments.
#'
#' @references
#' Bryant, J. L., & Paulson, A. S. (1976). An extension of Tukey's method of
#'   multiple comparisons to experimental designs with random concomitant
#'   variables. \emph{Biometrika, 63}, 631--638.
#'   \doi{10.1093/biomet/63.3.631}
#'
#' Bryant, J. L., & Bruvold, N. T. (1980). Multiple comparison procedures in
#'   the analysis of covariance. \emph{Journal of the American Statistical
#'   Association, 75}(372), 874--880. \doi{10.2307/2287175}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#'   experiments and analyzing data: A model comparison perspective}
#'   (4th ed.). Routledge. (See Chapter 9.)
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso
#' \code{\link{ci_c_ancova_bp}} for the simultaneous confidence intervals
#'   these critical values produce; \code{\link[stats]{ptukey}} and
#'   \code{\link[stats]{qtukey}} for the ordinary (fixed-covariate or
#'   no-covariate) studentized range.
#'
#' @examples
#' # Critical value from the worked example of Bryant and Bruvold (1980):
#' # k = 6 panels, p = 1 covariate, nu = 14 error df, alpha = .05. The quantile
#' # is found by inverting the distribution function with uniroot; it is 4.83,
#' # the entry in Table 1 of Bryant and Paulson (1976).
#' qbryant_paulson(0.95, n_covariates = 1, n_groups = 6, df = 14)
#'
#' # The ordinary Tukey value (ignoring that the covariate is random and
#' # estimated) is smaller, 4.64, so it yields intervals that are too narrow:
#' qtukey(0.95, nmeans = 6, df = 14)
#'
#' # How much too narrow: the Bryant-Paulson area beyond the Tukey value is the
#' # familywise error rate Tukey's method actually delivers in this design.
#' # With no covariate the two distributions coincide, so the area is .0499,
#' # the nominal .05 up to the rounding of 4.64 itself. Each additional random
#' # covariate carries more estimation uncertainty, stretches the distribution
#' # to the right, and pushes the rate up, to .063, .076, and .091.
#' pbryant_paulson(4.64, n_covariates = 0:3, n_groups = 6, df = 14,
#'                 lower_tail = FALSE)
#'
#' # The p = 0 entry above is exactly the ordinary studentized range:
#' ptukey(4.64, nmeans = 6, df = 14, lower.tail = FALSE)
#'
#' @keywords distribution design
#'
#' @name bryant_paulson
NULL


# Internal engine: the cumulative distribution function of the Bryant-Paulson
# Q_p statistic at a single (q, p, k, nu). The formula is derived in Bryant and
# Paulson (1976); equation numbers below refer to that paper.
#
# Bryant and Paulson (1976, Equation 17) give the exact distribution function
# of Q_p as a single integral over a variable A that bundles together two
# independent sources of variability: the error estimate g, distributed as
# chi squared on nu degrees of freedom (Bryant & Paulson, 1976, Equation 9,
# with v/v estimating sigma^2); and a random covariate-shrinkage factor
# delta = {1 + t' Sigma_xx^{-1} t / (nu + 1)}^{-1}, which has a Beta
# distribution with parameters (nu + 1)/2 and p/2 (Bryant & Paulson, 1976,
# Equations 11 and 12). The two are combined as A = (g * delta / 2)^{1/2}
# (Bryant & Paulson, 1976, Equation 13).
#
# This implementation evaluates the factored form of Equation (17). In the
# conditioning argument of Bryant and Paulson (1976, p. 634), the event
# {Q_p <= q} is the event that the range of k independent standard normal
# variates is at most (2/nu)^{1/2} * q * A. Substituting A = (g * delta / 2)^{1/2}
# with g = nu * S^2 (so that S^2 is the error variance estimate and nu * S^2 is
# distributed as chi squared on nu degrees of freedom) turns (2/nu)^{1/2} * q * A
# into q * S * sqrt(delta). Writing R_k for the range of k independent standard
# normals, Q_p is distributed as R_k / (sqrt(delta) * S), with delta and S
# independent. Conditioning on delta, and writing G(x; k, nu) = P(R_k / S <= x)
# for the studentized-range CDF, gives
#
#     P(Q_p <= q) = E_delta[ G(q * sqrt(delta); k, nu) ]
#                 = integral over delta in (0, 1) of
#                       G(q * sqrt(delta); k, nu) *
#                       dbeta(delta, (nu + 1)/2, p/2)  d(delta).
#
# The studentized-range CDF G is .sr_cdf() (R/studentized_range_internals.R):
# stats::ptukey() where it is accurate, and an exact evaluation that does not
# call ptukey where it is not. ptukey's algorithm (Copenhaver & Holland, 1988)
# loses accuracy at small nu: at nu = 3, k = 20 its error reaches ~3e-4 in
# probability, enough to shift the tabled critical value by ~0.2, and it is
# larger still at nu = 2. Bryant and Paulson's (1976) own hand-tabulated values
# are accurate there; a direct simulation of Q_p confirms them and disagrees
# with ptukey. The exact path makes qbryant_paulson()/pbryant_paulson()
# reproduce the original Bryant and Paulson (1976) Table 1 across the whole df
# range, including nu = 2 and 3; see tests/testthat/test-cv_bryant_paulson.R.
# Until 1.1.0 the exact path was used for nu < 7 and lived in this file. It now
# lives in R/studentized_range_internals.R and serves every studentized range
# caller, over the region measured there (ptukey's errors reach past nu = 7
# when k is large).
#
# When p = 0 the factor delta is degenerate at 1 and Q_p is exactly the ordinary
# studentized range (Bryant & Paulson, 1976, Section 1); p = 0 returns .sr_cdf()
# directly so the documented reduction to sqrt(2) * cv_tukey_hsd() holds exactly.

# The integral over the covariate-shrinkage factor is carried out after the
# change of variables delta = 1 - u^2. The Beta((nu + 1)/2, p/2) density of delta
# (Bryant & Paulson, 1976, Equation 12) carries the factor (1 - delta)^(p/2 - 1),
# which for p = 1 is an integrable singularity at delta = 1, right where the
# distribution puts most of its mass; integrate() copes with it only by
# subdividing heavily near that endpoint, so the CDF cost several hundred
# studentized-range evaluations. Under delta = 1 - u^2 the weight becomes
#     2 u^(p - 1) (1 - u^2)^((nu - 1)/2) / B((nu + 1)/2, p/2),
# which is smooth on [0, 1] for every p >= 1, and the same integral converges in
# a handful of subdivisions: a critical value costs about a twentieth of what it
# did, and the two forms agree to about 1e-10. The change of variables is exact,
# so the Table 1 reproduction in the tests is unaffected. The weight is taken
# on the log scale: the direct power loses six digits at 10^10 df, and with
# many covariates u^(p - 1) underflows where 1 / B() overflows (NaN with 200
# covariates at 10^5 df).
.bp_shrink_weight <- function(u, p, nu) {
  a <- (nu + 1) / 2
  log_u <- if (p == 1) 0 else (p - 1) * log(u)
  log_power <- if (a == 1) 0 else (a - 1) * log1p(-u^2)
  2 * exp(log_u + log_power - lbeta(a, p / 2))
}

# Upper limit of the integral over u: u^2 = 1 - delta is Beta(p/2, (nu + 1)/2),
# so the weight above this limit is .BP_U_TAIL. At large nu the weight lies
# within a few multiples of sqrt(p / nu) of 0, between the nodes of a
# quadrature over [0, 1]: from about 4e7 df with one covariate the CDF came
# out near 0 (5.6e-25 at 1e8 df, where it is .956). The limit is above .998
# below 10 df and about 8 / sqrt(nu) for one covariate at large nu.
.BP_U_TAIL <- 1e-15
.bp_u_upper <- function(p, nu)
  min(1, sqrt(stats::qbeta(.BP_U_TAIL, p / 2, (nu + 1) / 2,
                           lower.tail = FALSE)))

.pbp_one <- function(q, p, k, nu, exact = FALSE) {
  if (is.na(q) || is.na(p) || is.na(k) || is.na(nu)) return(NA_real_)
  if (q <= 0) return(0)
  # p = 0: ordinary studentized range (Bryant & Paulson, 1976, Section 1).
  # nu = Inf: delta = {1 + t' Sigma_xx^{-1} t / (nu + 1)}^{-1} is 1 with
  # probability 1 and the error variance is known, so Q_p is the studentized
  # range at infinite df, whatever the number of covariates.
  if (p == 0 || is.infinite(nu)) return(.sr_cdf(q, k, nu))
  u_upper <- .bp_u_upper(p, nu)
  # exact = TRUE integrates the exact studentized-range CDF at every node.
  # The far tails need it: ptukey()'s small absolute error there (about
  # 9e-9 at q = 12 with 10 groups and 11 df) moves a quantile beyond .999
  # by as much as 7e-5 (qbryant_paulson(.9999, 1, 10, 11)).
  if (exact) {
    rs <- .sr_range_spline(k)
    integrand_exact <- function(u)
      vapply(q * sqrt(1 - u^2), .sr_cdf_exact_one, numeric(1),
             k = k, nu = nu, rs = rs) * .bp_shrink_weight(u, p, nu)
    val <- stats::integrate(integrand_exact, 0, u_upper, rel.tol = 1e-10,
                            stop.on.error = FALSE)$value
    return(min(max(val, 0), 1))
  }
  # Studentized-range CDF from .sr_cdf() (exact where ptukey is not accurate,
  # ptukey otherwise; see the note above), integrated against the
  # covariate-shrinkage weight (Eq. 12) in the substituted variable u, with
  # delta = 1 - u^2. The quadrature tolerances are those the two paths had
  # before 1.1.0: 1e-10 where the CDF is exact, 1e-9 where it is ptukey.
  integrand <- function(u)
    .sr_cdf(q * sqrt(1 - u^2), k, nu) * .bp_shrink_weight(u, p, nu)
  val <- stats::integrate(integrand, 0, u_upper,
                          rel.tol = if (.sr_needs_exact(k, nu)) 1e-10 else 1e-9,
                          stop.on.error = FALSE)$value
  min(max(val, 0), 1)
}


#' @rdname bryant_paulson
#' @export
pbryant_paulson <- function(q, n_covariates, n_groups, df,
                            lower_tail = TRUE) {
  if (any(n_covariates < 0, na.rm = TRUE))
    .dmar_refuse_argument(
      "'n_covariates' (p) must be non-negative.", arg = "n_covariates")
  if (any(n_groups < 2, na.rm = TRUE))
    .dmar_refuse_argument(
      "'n_groups' (k) must be at least 2.", arg = "n_groups")
  if (any(df <= 0, na.rm = TRUE))
    .dmar_refuse_argument("'df' (nu) must be positive.", arg = "df")
  n <- max(length(q), length(n_covariates), length(n_groups), length(df))
  q   <- rep_len(q, n)
  p   <- rep_len(n_covariates, n)
  k   <- rep_len(n_groups, n)
  nu  <- rep_len(df, n)
  out <- vapply(seq_len(n), function(i) {
    val <- .pbp_one(q[i], p[i], k[i], nu[i])
    # In the far tails, recompute with the exact studentized range.
    if (!is.na(val) && p[i] > 0 && (val > 0.999 || val < 0.001))
      val <- .pbp_one(q[i], p[i], k[i], nu[i], exact = TRUE)
    val
  }, numeric(1))
  if (!lower_tail) out <- 1 - out
  out
}


#' @rdname bryant_paulson
#' @export
qbryant_paulson <- function(prob, n_covariates, n_groups, df,
                            lower_tail = TRUE) {
  if (any(n_covariates < 0, na.rm = TRUE))
    .dmar_refuse_argument(
      "'n_covariates' (p) must be non-negative.", arg = "n_covariates")
  if (any(n_groups < 2, na.rm = TRUE))
    .dmar_refuse_argument(
      "'n_groups' (k) must be at least 2.", arg = "n_groups")
  if (any(df <= 0, na.rm = TRUE))
    .dmar_refuse_argument("'df' (nu) must be positive.", arg = "df")
  if (!lower_tail) prob <- 1 - prob
  if (any(prob <= 0 | prob >= 1, na.rm = TRUE))
    .dmar_refuse_argument(
      "'prob' must be strictly between 0 and 1.", arg = "prob")
  n <- max(length(prob), length(n_covariates), length(n_groups),
           length(df))
  prob <- rep_len(prob, n)
  p    <- rep_len(n_covariates, n)
  k    <- rep_len(n_groups, n)
  nu   <- rep_len(df, n)
  vapply(seq_len(n), function(i) {
    if (is.na(prob[i]) || is.na(p[i]) || is.na(k[i]) || is.na(nu[i]))
      return(NA_real_)
    # No covariates, or infinite df: the ordinary studentized range (see
    # .pbp_one()).
    if (p[i] == 0 || is.infinite(nu[i]))
      return(.sr_quantile(prob[i], k[i], nu[i]))
    # Beyond .999 and below .001 the CDF integrates the exact studentized
    # range, as .sr_quantile() does for the ordinary studentized range.
    far <- prob[i] > 0.999 || prob[i] < 0.001
    f <- function(qq) .pbp_one(qq, p[i], k[i], nu[i], exact = far) - prob[i]
    # Bracket: the BP quantile exceeds the Tukey quantile (p > 0 widens the
    # distribution); use the Tukey value as a lower bound and grow an upper
    # bound until the CDF clears prob, tested with f itself (at 10^10 df a
    # bound found with the other CDF missed the far-tail root).
    lo <- .sr_quantile(prob[i], k[i], nu[i])
    hi <- lo
    f_hi <- f(hi)
    while (f_hi < 0 && hi < 1e4) {
      hi <- hi * 2
      f_hi <- f(hi)
    }
    stats::uniroot(f, lower = max(lo * 0.5, 1e-6), upper = hi,
                   f.upper = f_hi, tol = 1e-10)$root
  }, numeric(1))
}


#' @rdname bryant_paulson
#' @export
dbryant_paulson <- function(q, n_covariates, n_groups, df) {
  # Numerical derivative of the CDF; adequate for plotting/illustration.
  h <- 1e-4
  (pbryant_paulson(q + h, n_covariates, n_groups, df) -
     pbryant_paulson(q - h, n_covariates, n_groups, df)) / (2 * h)
}
