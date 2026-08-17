#' The Bryant--Paulson Generalized Studentized Range Distribution
#'
#' @description
#' Distribution function (\code{pbryant_paulson}), quantile/critical-value
#' function (\code{qbryant_paulson}), and density (\code{dbryant_paulson})
#' for the Bryant--Paulson generalized studentized range, the sampling
#' distribution of the studentized range of covariate-\emph{adjusted} means
#' in the analysis of covariance (ANCOVA) when the covariate(s) are
#' \emph{random}. These are the analysis-of-covariance analogues of
#' \code{\link[stats]{ptukey}} / \code{\link[stats]{qtukey}} and supply the
#' critical values needed for Tukey--Kramer-type simultaneous confidence
#' intervals on (and tests of) contrasts of adjusted means.
#'
#' @param q Vector of quantiles (values of the generalized studentized
#'   range statistic).
#' @param prob Vector of probabilities. For \code{qbryant_paulson} this is
#'   the cumulative probability (e.g., \code{0.95} returns the upper 5\%
#'   critical value).
#' @param num_covariates The number of random covariates, \eqn{p}
#'   (\eqn{p \ge 0}). With \code{num_covariates = 0} the distribution is
#'   exactly the ordinary studentized range and the functions reduce to
#'   \code{\link[stats]{ptukey}} / \code{\link[stats]{qtukey}}.
#' @param num_groups The number of groups (treatments) being compared,
#'   \eqn{k} (\eqn{k \ge 2}); this is the \dQuote{sample size} of the range.
#' @param df The error degrees of freedom of the ANCOVA model, \eqn{\nu}.
#'   In a one-way ANCOVA with \eqn{N} total observations, \eqn{k} groups,
#'   and \eqn{p} covariates, \eqn{\nu = N - k - p}.
#' @param lower_tail Logical; if \code{TRUE} (default) probabilities are
#'   \eqn{P(Q \le q)}, otherwise \eqn{P(Q > q)}.
#' @param \dots Additional arguments (currently unused; for extensibility).
#'
#' @details
#' \strong{The statistic.} In a balanced ANCOVA with \eqn{k} groups and
#' \eqn{p} random covariates, let \eqn{\hat\theta_i} be the adjusted group
#' means and \eqn{\hat\sigma_{y \mid x}} the square root of the ANCOVA error
#' mean square (on \eqn{\nu} degrees of freedom). The Bryant--Paulson
#' statistic is the studentized range of the adjusted means,
#' \deqn{Q \;=\; \frac{\max_i \hat\theta_i - \min_i \hat\theta_i}{\hat\sigma_{y\mid x}\sqrt{K_1 - K_2}},}
#' where \eqn{K_1 - K_2} is the design constant that scales the variance of a
#' single adjusted mean (for a one-way design with \eqn{n} per group,
#' \eqn{K_1 - K_2 = 1/n}). Crucially, the studentizer uses only this
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
#' short-circuits to \code{ptukey}. The integral is evaluated with
#' \code{\link[stats]{integrate}}; \code{qbryant_paulson} inverts it with
#' \code{\link[stats]{uniroot}}. Bryant and Bruvold (1980) later showed the
#' same distribution and critical values remain valid when the covariates are
#' \emph{not} identically distributed across groups (their grouped-covariate
#' model, Eq. 1.3), and added the Duncan multiple-range extension.
#'
#' \strong{Accuracy at small df.} The error integral is carried out with
#' \code{\link[stats]{ptukey}} for \eqn{\nu \ge 7}, where it is accurate to
#' about \eqn{10^{-9}}. For \eqn{\nu < 7} \code{ptukey}'s algorithm loses
#' accuracy (at \eqn{\nu = 3}, \eqn{k = 20} its probability error reaches
#' \eqn{\approx 3\times10^{-4}}, enough to move the critical value by about
#' 0.2, and it is larger at \eqn{\nu = 2}), so the studentized-range
#' distribution is instead evaluated directly, without \code{ptukey}, by
#' integrating the probability integral of the range against the
#' \eqn{\chi^2_\nu} error density. The small-\eqn{\nu} path costs a fraction of
#' a second.
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
#' # k = 6 panels, p = 1 covariate, nu = 14 error df, alpha = .05. Getting a
#' # quantile means inverting the distribution function with uniroot, and every
#' # step of that root search evaluates the integral over the covariate-shrinkage
#' # factor, so the call takes about half a second and is shown here rather than
#' # run.
#' # qbryant_paulson(0.95, num_covariates = 1, num_groups = 6, df = 14)
#' # It returns 4.83, the entry in Table 1 of Bryant and Paulson (1976). The
#' # distribution function itself is a single integral and is quick, so the
#' # pbryant_paulson calls below do run.
#'
#' # The ordinary Tukey value (ignoring that the covariate is random and
#' # estimated) is smaller, so it yields intervals that are too narrow:
#' qtukey(0.95, nmeans = 6, df = 14)                                   # 4.64
#'
#' # How much too narrow: the Bryant-Paulson area beyond the Tukey value is the
#' # familywise error rate Tukey's method actually delivers in this design.
#' # With no covariate the two distributions coincide, so the area is .0499,
#' # the nominal .05 up to the rounding of 4.64 itself. Each additional random
#' # covariate carries more estimation uncertainty, stretches the distribution
#' # to the right, and pushes the rate up, to .063, .076, and .091.
#' pbryant_paulson(4.64, num_covariates = 0:3, num_groups = 6, df = 14,
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
# The studentized-range CDF G is stats::ptukey(). For nu >= 7 ptukey is accurate
# to ~1e-9 and is used directly. For small nu, however, ptukey's algorithm
# (AS 190; Copenhaver & Holland, 1988) loses accuracy: at nu = 3, k = 20 its
# error reaches ~3e-4 in probability, enough to shift the tabled critical value
# by ~0.2, and it is larger still at nu = 2. Bryant and Paulson's (1976) own
# hand-tabulated values are accurate there; a direct simulation of Q_p confirms
# them and disagrees with ptukey. So for nu < 7 the studentized-range CDF is
# evaluated exactly by .bp_qsr_cdf() below, which does not call ptukey. This
# makes qbryant_paulson()/pbryant_paulson() reproduce the original Bryant and
# Paulson (1976) Table 1 across the whole df range, including nu = 2 and 3; see
# tests/testthat/test-cv_bryant_paulson.R.
#
# When p = 0 the factor delta is degenerate at 1 and Q_p is exactly the ordinary
# studentized range (Bryant & Paulson, 1976, Section 1); p = 0 returns ptukey()
# directly so the documented reduction to sqrt(2) * cv_tukey_hsd() holds exactly.

# Degrees of freedom below which ptukey is bypassed for the exact CDF.
.BP_NU_EXACT <- 7L

# CDF of the range of k iid N(0,1) at w >= 0:
#   P(W <= w) = k * integral phi(u) [Phi(u) - Phi(u - w)]^{k-1} du,
# the "probability integral of the range" (David & Nagaraja, Order Statistics).
# phi(u) confines the integrand to about [-8.5, 8.5] for any w.
.bp_range_cdf <- function(w, k) {
  if (w <= 0) return(0)
  f <- function(u) {
    d <- stats::pnorm(u) - stats::pnorm(u - w); d[d < 0] <- 0
    k * stats::dnorm(u) * d^(k - 1)
  }
  # rel.tol 1e-12 (subdivisions 1000) makes the knot values machine-accurate, so
  # the interpolating spline built from them, not the quadrature, is the only
  # remaining approximation. The dnorm(u) weight confines the integrand to about
  # [-8.5, 8.5] for any w; the mass beyond those bounds is below 1e-15.
  min(1, stats::integrate(f, -8.5, 8.5, rel.tol = 1e-12,
                          subdivisions = 1000L)$value)
}

# Monotone spline of the range CDF for a given k, up to where it saturates at 1.
# The CDF depends only on k, so splines are memoized by k in .bp_spline_cache:
# building one (a few thousand range integrals) is paid once per k per session,
# and every subsequent CDF evaluation is a spline lookup. The spline stands in
# for the (machine-accurate) direct integral inside the nested quadrature of
# .bp_qsr_cdf(), which would otherwise be evaluated tens of thousands of times
# per critical value; so the spline's interpolation error is what limits the
# accuracy of the whole small-nu path. Knots are placed at a fixed spacing of
# about 0.003 in w (a few thousand of them), which holds the monoH.FC
# interpolation error below ~1e-9 across k = 2, ..., 20; that maps to a critical
# value accurate to about 1e-7, matching the direct integral to seven figures.
# A coarser 400-knot grid (the earlier default) left a ~2e-7 range-CDF error
# that inflated to ~4e-5 in q_{.01; 2,8,3}. The monoH.FC method is Fritsch and
# Carlson's monotone cubic, so the interpolant is a genuine (non-decreasing) CDF,
# which keeps uniroot's bracketing in qbryant_paulson robust.
.bp_spline_cache <- new.env(parent = emptyenv())
.bp_get_range_spline <- function(k) {
  key <- as.character(k)
  hit <- get0(key, envir = .bp_spline_cache, inherits = FALSE)
  if (!is.null(hit)) return(hit)
  wsat <- 2
  while (.bp_range_cdf(wsat, k) < 1 - 1e-13 && wsat < 80) wsat <- wsat + 0.5
  n_knots <- max(1500L, as.integer(ceiling(wsat / 0.003)) + 1L)
  wg <- seq(0, wsat, length.out = n_knots)
  rc <- vapply(wg, function(w) .bp_range_cdf(w, k), numeric(1))
  sp <- stats::splinefun(wg, rc, method = "monoH.FC")
  rs <- list(f = function(w) pmin(1, pmax(0, sp(w))), wsat = wsat)
  assign(key, rs, envir = .bp_spline_cache)
  rs
}

# Exact studentized-range CDF G(x; k, nu) = P(R_k / S <= x), S = sqrt(chi^2_nu/nu),
# without ptukey: integrate the range CDF against the chi squared error density.
# The range CDF saturates at 1 once its argument exceeds wsat, i.e. once
# cc > (wsat/x)^2 * nu; beyond that point the integrand is exactly the chi squared
# density, so that tail is added analytically (1 - pchisq). Splitting there keeps
# both pieces smooth (no derivative kink) and the quadrature well-behaved.
.bp_qsr_cdf <- function(x, k, nu, rs) {
  if (x <= 0) return(0)
  if (is.infinite(nu)) return(rs$f(x))
  chi_hi <- stats::qchisq(1 - 1e-12, nu)
  ccstar <- (rs$wsat / x)^2 * nu
  upper  <- min(ccstar, chi_hi)
  below  <- stats::integrate(function(cc)
      rs$f(x * sqrt(cc / nu)) * stats::dchisq(cc, nu),
      0, upper, rel.tol = 1e-11, stop.on.error = FALSE)$value
  tail <- if (ccstar < chi_hi) 1 - stats::pchisq(ccstar, nu) else 0
  min(1, max(0, below + tail))
}

.pbp_one <- function(q, p, k, nu) {
  if (is.na(q) || is.na(p) || is.na(k) || is.na(nu)) return(NA_real_)
  if (q <= 0) return(0)
  # p = 0: ordinary studentized range (Bryant & Paulson, 1976, Section 1).
  if (p == 0) return(stats::ptukey(q, nmeans = k, df = nu))
  # Studentized-range CDF: exact for small nu, ptukey otherwise (see the note
  # above). dbeta((nu+1)/2, p/2) is the covariate-shrinkage factor (Eq. 12).
  if (nu < .BP_NU_EXACT) {
    rs <- .bp_get_range_spline(k)
    integrand <- function(delta)
      vapply(delta, function(dd) .bp_qsr_cdf(q * sqrt(dd), k, nu, rs), numeric(1)) *
        stats::dbeta(delta, (nu + 1) / 2, p / 2)
    val <- stats::integrate(integrand, 0, 1, rel.tol = 1e-10,
                            stop.on.error = FALSE)$value
    return(min(max(val, 0), 1))
  }
  integrand <- function(delta)
    stats::ptukey(q * sqrt(delta), nmeans = k, df = nu) *
      stats::dbeta(delta, (nu + 1) / 2, p / 2)
  val <- stats::integrate(integrand, 0, 1, rel.tol = 1e-9,
                          stop.on.error = FALSE)$value
  min(max(val, 0), 1)
}


#' @rdname bryant_paulson
#' @export
pbryant_paulson <- function(q, num_covariates, num_groups, df,
                            lower_tail = TRUE, ...) {
  if (any(num_covariates < 0, na.rm = TRUE))
    stop("'num_covariates' (p) must be non-negative.")
  if (any(num_groups < 2, na.rm = TRUE))
    stop("'num_groups' (k) must be at least 2.")
  if (any(df <= 0, na.rm = TRUE))
    stop("'df' (nu) must be positive.")
  n <- max(length(q), length(num_covariates), length(num_groups), length(df))
  q   <- rep_len(q, n)
  p   <- rep_len(num_covariates, n)
  k   <- rep_len(num_groups, n)
  nu  <- rep_len(df, n)
  out <- vapply(seq_len(n),
                function(i) .pbp_one(q[i], p[i], k[i], nu[i]),
                numeric(1))
  if (!lower_tail) out <- 1 - out
  out
}


#' @rdname bryant_paulson
#' @export
qbryant_paulson <- function(prob, num_covariates, num_groups, df,
                            lower_tail = TRUE, ...) {
  if (any(num_covariates < 0, na.rm = TRUE))
    stop("'num_covariates' (p) must be non-negative.")
  if (any(num_groups < 2, na.rm = TRUE))
    stop("'num_groups' (k) must be at least 2.")
  if (any(df <= 0, na.rm = TRUE))
    stop("'df' (nu) must be positive.")
  if (!lower_tail) prob <- 1 - prob
  if (any(prob <= 0 | prob >= 1, na.rm = TRUE))
    stop("'prob' must be strictly between 0 and 1.")
  n <- max(length(prob), length(num_covariates), length(num_groups),
           length(df))
  prob <- rep_len(prob, n)
  p    <- rep_len(num_covariates, n)
  k    <- rep_len(num_groups, n)
  nu   <- rep_len(df, n)
  vapply(seq_len(n), function(i) {
    if (is.na(prob[i]) || is.na(p[i]) || is.na(k[i]) || is.na(nu[i]))
      return(NA_real_)
    if (p[i] == 0)
      return(stats::qtukey(prob[i], nmeans = k[i], df = nu[i]))
    f <- function(qq) .pbp_one(qq, p[i], k[i], nu[i]) - prob[i]
    # Bracket: the BP quantile exceeds the Tukey quantile (p > 0 widens the
    # distribution); use the Tukey value as a lower bound and grow an upper
    # bound until the CDF clears prob.
    lo <- stats::qtukey(prob[i], nmeans = k[i], df = nu[i])
    hi <- lo
    while (.pbp_one(hi, p[i], k[i], nu[i]) < prob[i] && hi < 1e4) hi <- hi * 2
    stats::uniroot(f, lower = max(lo * 0.5, 1e-6), upper = hi,
                   tol = 1e-10)$root
  }, numeric(1))
}


#' @rdname bryant_paulson
#' @export
dbryant_paulson <- function(q, num_covariates, num_groups, df, ...) {
  # Numerical derivative of the CDF; adequate for plotting/illustration.
  h <- 1e-4
  (pbryant_paulson(q + h, num_covariates, num_groups, df) -
     pbryant_paulson(q - h, num_covariates, num_groups, df)) / (2 * h)
}
