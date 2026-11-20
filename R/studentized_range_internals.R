# The studentized range distribution: one engine for every caller.
#
# Q = R_k / S, where R_k is the range of k independent standard normal
# variates and nu * S^2 is an independent chi square variate on nu degrees of
# freedom. Every DMAR function that needs a quantile or a probability of Q
# calls .sr_quantile() or .sr_cdf() below, never stats::qtukey() or
# stats::ptukey() directly: cv_tukey_hsd(), ci_tukey_kramer(),
# ci_games_howell(), contrast_test(adjust = "tukey"), and the Bryant-Paulson
# functions (whose n_covariates = 0 case is Q itself, and whose general case
# integrates the distribution function of Q).
#
# Why not stats::ptukey() everywhere. R's documentation states that ptukey()
# uses a 16-point Legendre formula for its integral and that qtukey(), a
# secant search on ptukey(), is accurate to the fourth decimal place. Two
# independent quadratures (the 1.1.0 review's F3 and G2 evidence) show where
# that holds and where it does not. The G2 quadrature was compared with R's
# values at 52,892 cells: k = 3 to 10000 means at error df nu from just below
# each small-nu bound to 20000; nu = Inf, the distribution R uses for every
# nu above 25000, for k = 2 to 20000; nu = 1e8 + 1, 1e9, and 1e10 for k = 3
# to 10000; lower-tail probabilities p = .00001 to .99999. The F3 quadrature
# covered k = 2 to 50 at nu = 1 to 1000 and Inf, and k up to 500 at nu up to
# 10^6 and Inf.
#   - small nu, increasingly so as k and p grow: at k = 20, nu = 2,
#     qtukey(.99) is 44.526 where the exact quantile is 37.943 (17 percent
#     high), at nu = 3 it is 0.9 percent high, at k = 50, nu = 12 and
#     p = .999 it is still 4.8e-4 high; below nu = 2 both functions return
#     NaN. The df at which the errors fall below the tolerances grows with k
#     (.SR_EXACT_DF_BOUNDS below), to 100 df at 3000 to 10000 means;
#   - nu above 10000: ptukey()'s probability error grows from about 1e-8 at
#     15000 to 3e-8 at 18000, which moves qtukey(.999) by 1.2e-5, and from
#     25001 on ptukey() returns the infinite-df distribution, so that at
#     nu = 25001 qtukey(.99, 50) is 1.1e-3 low and ptukey() is off by up to
#     9.5e-5 in probability. The gap closes as 1/nu, and above 10^8 the
#     infinite-df values are accurate again (.SR_EXACT_UP_TO_DF below);
#   - far tails: ptukey()'s upper tail carries an absolute error of about
#     3e-8 (more with several hundred means), which qtukey()'s secant search
#     turns into a quantile error of 1e-4 at p = .9999 (qtukey(.9999, 10, 11)
#     is 12.0084624 where the exact quantile is 12.0085798) and of up to 4 at
#     p = .99999 (qtukey(.99999, 2000, 1000) is 12.93 where the exact
#     quantile is 9.87; at nu = Inf it is 13.75 against 9.75). In the lower
#     tail qtukey() returns values far from the root (qtukey(.00001, 3, 150)
#     is 9.2e-5 where the exact quantile is 6.0e-3) or NaN;
#   - many means below the upper quarter: with more than 50 means ptukey()
#     is off by more than 1e-6 at lower-tail probabilities below .75 (up to
#     about 8.6e-6 with 2000 means, and with 300 or more means it returns 0 where the
#     probability is 1e-5, and with 1500 or more where it is 1e-3); with
#     thousands of means it is off just above .75 as well: at nu = Inf from
#     about 4500 means on (1.04e-6 at .75 with 5000 means, 1.18e-6 with
#     10000, within 1e-6 again from about .765), and at nu = 1e8 + 1, where
#     the infinite-df value is returned, from 2500 means on (1.07e-6 at .75
#     with 3000 means);
#   - anywhere: qtukey()'s secant search can stop short of ptukey()'s root:
#     qtukey(.975, 10, 12) is 4.1e-5 below the exact 6.0038073,
#     qtukey(.9975, 5, 23) 1.2e-4 above, qtukey(.975, 150, 300) 8.3e-5 above.
# The exact computation is used in each of these places and nowhere else:
# .sr_needs_exact() marks the (k, nu) regions, .sr_needs_exact_p() the
# probabilities, and .sr_qtukey_offset() checks each value qtukey() returns
# against ptukey() and sends one that stopped short to the exact root.
# Everywhere else the results are exactly what stats::qtukey() and
# stats::ptukey() give. There, over the G2 cells, qtukey() was within 9.6e-6
# of the exact quantile and ptukey() within 9.1e-7 of the exact probability
# (k = 50, nu = Inf, p = .1). The largest quantile error came from a dense
# probe of 16,422 cells between the grid's cells near the df bounds, where the
# kept errors come closest to 1e-5 (k = 700, nu = 57, p = .9965); at the
# grid's own cells it was about 8.6e-6. Larger errors between the measured cells
# cannot be ruled out, and beyond 10000 means (20000 at nu = Inf) the df bound
# and the probability cuts are extrapolated, not measured.

# The small-nu region, as the measurements placed it. Below the bound for
# k's row, qtukey() was more than 1e-5 from the exact quantile (after the
# check in .sr_qtukey_offset()), or ptukey() more than 1e-6 from the exact
# probability, at some measured p from .001 to .999 (.75 to .999 with more
# than 50 means). Each bound sits at least one measured df above the largest
# failing df: k = 2 to 6 fail up to nu = 6 or 7; k = 8 and 10 up to 8 and 9;
# k = 15 and 20 up to 10 and 13; k = 30 and 50 up to 15 and 19; k = 100 up to
# 20; k = 200 up to 30 (k = 250: none from 39, the first df measured there);
# k = 300 and 400 up to 42 and 44; 500 and 700 up to 52 and 54; 1000 up to 65;
# 1500 up to 75; 2000 up to 85; 2500 up to 90; 3000, 4000, 5000, 7000, and
# 10000 up to 100 (the G2 measurements, which also confirmed the rows up to
# 250 means with a finer grid of probabilities). With 2500 to 10000 means
# R's errors fall by one to two orders of magnitude between 100 and 100.25
# df, and 100.25, 100.5, 101, and 102 df all passed. The last measured row is
# k = 10000, in the row for more than 5000 means; beyond 10000 means its bound
# is extrapolated.
.SR_EXACT_DF_BOUNDS <- data.frame(
  max_k    = c(6, 10, 20, 50, 100, 250, 400, 700, 1000, 1500, 2000, 3000,
               5000, Inf),
  below_df = c(9, 11, 15, 21, 31, 41, 47, 57, 71, 81, 91, 101, 111, 111))

# Above this many degrees of freedom (finite) ptukey() is not accurate either:
# measured failing at 15000 (qtukey's secant, caught by the check as well) and
# 18000 (ptukey's own error, which no check on ptukey can see), accurate at
# 5000 and 10000.
.SR_EXACT_ABOVE_DF <- 10000

# The large-nu region ends here. From 25001 df on, ptukey() and qtukey()
# return the infinite-df distribution, which differs from the distribution at
# nu by about c / nu, with c up to 81 in the quantile and 9 in the probability
# over k = 2 to 500 and p = .05 to .9999. The infinite-df values meet the
# tolerances that set the small-nu bounds from about 5.7e6 df on (k up to 200,
# p = .90 to .999; 8.1e6 for the quantile over every measured cell), and at
# 1e8 they are within 8.1e-7 of the exact quantile and 9.2e-8 of the exact
# probability (k up to 500; with up to 10000 means, G2 measured the gap at
# 1e8 + 1 as at most 1.4e-6 in the quantile and 2.3e-7 in the probability).
# Above 1e8, as at nu = Inf, the results are therefore R's infinite-df values,
# as before 1.1.0, wherever R's values are accurate there (the probability
# rules and the check below still apply, with the exact infinite-df
# distribution). R's own infinite-df error adds to that gap, which is why
# the probability rule for more than 1000 means reaches .80 (see
# .sr_needs_exact_p()). The cap also keeps the exact computation away
# from nu at which double precision cannot place the chi square bounds of its
# integral:
# qchisq(1e-13, nu) / nu rounds to 1 from about 1e30 on, the integral then
# vanishes, and the distribution function would collapse to a step at the
# spline's saturation point (at k = 4 and p = .95 the relative error of the
# uncapped exact quantile was 1e-11 up to 1e12, 1.1e-8 at 1e16, and 2e-3 at
# 1e30, and from 1e40 on the quantile came out 11 rather than 3.633).
.SR_EXACT_UP_TO_DF <- 1e8

# TRUE where the exact computation replaces stats::ptukey() and
# stats::qtukey() for every probability: nu below the bound for k, or nu above
# 10000 and at most 1e8.
.sr_needs_exact <- function(k, nu) {
  row <- findInterval(k, .SR_EXACT_DF_BOUNDS$max_k, left.open = TRUE) + 1L
  nu < .SR_EXACT_DF_BOUNDS$below_df[row] |
    (nu > .SR_EXACT_ABOVE_DF & nu <= .SR_EXACT_UP_TO_DF)
}

# Probabilities (lower tail) at which R is not accurate, at any nu:
# quantiles below .001 or above .999, and, with more than 50 means, every
# probability below .75 (below .80 with more than 1000 means). Of the 1020 G2
# grid cells per probability outside the F3 region (24 values of k from 3 to
# 2000, df from just below each F3 bound to 20000), qtukey() missed the exact
# quantile by more than 1e-5 in 83 at p = .9995, 177 at .9999, 190 at .99995,
# and 291 at .99999, with failures at every k and at df from 9 to 20000 (and,
# with 2000 means, at Inf); in 236 at .00001 (NaN in 636 more) and 6 at .0001
# (NaN in 576), none at .001 (NaN in 466, which the fallback handles); with
# 2500 to 10000 means the same held. ptukey() failed the 1e-6 rule below its
# .75 point from k = 75 on; with at most 50 means it stayed within 9.1e-7 at
# every measured nu, Inf and 1e8 + 1 to 1e10 included. Above .75 it failed
# only with thousands of means and only at large nu: at Inf from about 4500
# means on, and at 1e8 + 1 (R's infinite-df value plus the gap to the
# distribution at that nu) from 2500 on, by up to 1.07e-6 at 3000 means
# (G2 round 2: 21_inf_df.log, 24b_above_cap.log). Hence the second cut:
# with more than 1000 means, .80. Where R's distribution function is kept,
# it was within 7.0e-7 at .80 and above with 1001 to 10000 means (4.7e-7 at
# Inf with 20000) and at .75 and above with 51 to 1000 means, at every
# measured nu (29b_claims.log).
.SR_TAIL_P <- 0.001
.SR_MANY_MEANS <- 50
.SR_MANY_MEANS_BELOW_P <- 0.75
.SR_VERY_MANY_MEANS <- 1000
.SR_VERY_MANY_MEANS_BELOW_P <- 0.80
.sr_needs_exact_p <- function(p, k, tails = TRUE) {
  cut <- ifelse(k > .SR_VERY_MANY_MEANS, .SR_VERY_MANY_MEANS_BELOW_P,
                .SR_MANY_MEANS_BELOW_P)
  many <- k > .SR_MANY_MEANS & p < cut
  if (tails) many | p < .SR_TAIL_P | p > 1 - .SR_TAIL_P else many
}

# The check on qtukey()'s value q: how far q is from ptukey()'s own root,
# estimated as the residual ptukey(q) - p over ptukey()'s density at q (a
# central difference with step 1e-3 q). On the measured grid the estimate was
# within 2.4e-7 of the true distance. A value more than 2e-6 from ptukey()'s
# root (the secant search stopped short) goes to the exact root; NaN (ptukey()
# warned) does too.
.SR_QTUKEY_CHECK_TOL <- 2e-6
.sr_qtukey_offset <- function(q, p, k, nu) {
  h <- 1e-3 * q
  f <- function(x) .sr_r_values(stats::ptukey, x, k, nu)
  (f(q) - p) / ((f(q + h) - f(q - h)) / (2 * h))
}

# CDF of the range of k iid N(0, 1) at w >= 0:
#   P(W <= w) = k * integral phi(u) [Phi(u) - Phi(u - w)]^{k-1} du,
# the "probability integral of the range" (David & Nagaraja, Order Statistics).
# phi(u) confines the integrand to about [-8.5, 8.5] for any w.
.sr_range_cdf <- function(w, k) {
  if (w <= 0) return(0)
  f <- function(u) {
    d <- stats::pnorm(u) - stats::pnorm(u - w); d[d < 0] <- 0
    k * stats::dnorm(u) * d^(k - 1)
  }
  # The integrand has a single peak, between u = 0 (the mode of phi) and
  # u = w / 2 (the mode of the window probability Phi(u) - Phi(u - w), which is
  # symmetric about w / 2); both factors are log-concave. With many means and
  # a small w the power k - 1 makes the peak narrow and puts it just below
  # w / 2 (at k = 7000, w = 6.2 it sits 0.02 below w / 2 with a standard
  # deviation of 0.08), and integrate() over the whole of [-8.5, 8.5] can miss
  # it: at k = 7000, w = 6.218556 it returned 2.3e-17 where the probability is
  # 9.5e-6. Splitting the range at w / 2 puts such a peak at the end of each
  # piece, where the Gauss-Kronrod nodes are closest together. rel.tol 1e-12
  # (subdivisions 1000) then makes the knot values accurate to about 1e-12, so
  # the interpolating spline built from them, not the quadrature, is the
  # remaining approximation (the 1.1.0 review's G2 evidence,
  # 26_range_cdf_audit.log). The dnorm(u) weight confines the integrand to
  # about [-8.5, 8.5] for any w; the mass beyond those bounds is below 1e-15.
  br <- c(-8.5, if (w / 2 < 8.5) w / 2, 8.5)
  val <- 0
  for (j in seq_len(length(br) - 1L))
    val <- val + stats::integrate(f, br[j], br[j + 1L], rel.tol = 1e-12,
                                  subdivisions = 1000L)$value
  min(1, val)
}

# Monotone spline of the range CDF for a given k, up to where it saturates at 1.
# The CDF depends only on k, so splines are memoized by k in .sr_spline_cache:
# building one (a few thousand range integrals) is paid once per k per session,
# and every subsequent CDF evaluation is a spline lookup. The spline stands in
# for the (machine-accurate) direct integral inside the outer quadrature of
# .sr_cdf_exact_one(), which would otherwise evaluate the range integral
# hundreds of times per probability; so the spline's interpolation error is
# what limits the accuracy of the exact path. Knots are placed at a fixed
# spacing of about 0.003 in w (a few thousand of them), which holds the
# monoH.FC interpolation error below ~1e-9 across k = 2, ..., 200; against the
# independent quadrature the exact quantiles agree to 5e-10 (relative) and the
# probabilities to 1.5e-10 (F3). With the range integral split at w / 2 the
# knot values agree with the G2 quadrature to 1.5e-12 at every knot checked
# (22 values of k from 3 to 20000, about 9500 points each), the range CDF
# reaches 1 - 1e-13 by w = 15.5 at all of them (in one piece it never did
# from about 8000 means on, which pushed wsat to 80 and the spline to 26,667
# knots), and the exact quantiles agree to 2.0e-8 at p = .01 to .99999,
# 4.6e-8 at .001, 1.1e-7 at .0001, and 4.8e-7 at .00001
# (26_range_cdf_audit.log, 29_summary.log). A coarser 400-knot grid (an earlier
# default) left a ~2e-7 range-CDF error that inflated to ~4e-5 in the
# Bryant-Paulson q_{.01; 2,8,3}. The monoH.FC method is Fritsch and Carlson's
# monotone cubic, so the interpolant is a genuine (non-decreasing) CDF, which
# keeps the root bracketing of the quantile functions robust.
.sr_spline_cache <- new.env(parent = emptyenv())
.sr_range_spline <- function(k) {
  key <- as.character(k)
  hit <- get0(key, envir = .sr_spline_cache, inherits = FALSE)
  if (!is.null(hit)) return(hit)
  wsat <- 2
  while (.sr_range_cdf(wsat, k) < 1 - 1e-13 && wsat < 80) wsat <- wsat + 0.5
  n_knots <- max(1500L, as.integer(ceiling(wsat / 0.003)) + 1L)
  wg <- seq(0, wsat, length.out = n_knots)
  rc <- vapply(wg, function(w) .sr_range_cdf(w, k), numeric(1))
  sp <- stats::splinefun(wg, rc, method = "monoH.FC")
  rs <- list(f = function(w) pmin(1, pmax(0, sp(w))), wsat = wsat)
  assign(key, rs, envir = .sr_spline_cache)
  rs
}

# Exact studentized-range CDF G(x; k, nu) = P(R_k / S <= x), S = sqrt(C / nu),
# without ptukey: integrate the range CDF against the chi square density of C.
# The integral runs over the central 1 - 2e-13 of the chi square distribution,
# so that a large nu, whose density is a narrow peak far from zero, is
# integrated where its mass is. The range CDF saturates at 1 once its argument
# exceeds wsat, i.e. once c > (wsat / x)^2 * nu; beyond that point the
# integrand is exactly the chi square density, so that tail is added
# analytically (1 - pchisq). Splitting there keeps both pieces smooth (no
# derivative kink) and the quadrature well-behaved. Against the independent
# quadrature this one interval was as accurate as splitting the range at a
# dozen chi square quantiles (722 cells, nu = 1 to 10^6), at a third of the
# cost (the 1.1.0 review's F3 evidence, cdf_variants.rds).
.sr_cdf_exact_one <- function(x, k, nu, rs = .sr_range_spline(k)) {
  if (x <= 0) return(0)
  # Beyond .SR_EXACT_UP_TO_DF the distribution is the infinite-df one to within
  # 1e-7 in probability, and at far larger nu the chi square bounds below
  # collapse in double precision. The dispatchers send such nu to R; this
  # guard serves the fallback for a cell where R's function failed there.
  if (nu > .SR_EXACT_UP_TO_DF) nu <- Inf
  if (k == 2) {
    # The range of two standard normals is sqrt(2) |Z|, so Q = sqrt(2) |T|.
    return(if (is.infinite(nu)) 2 * stats::pnorm(x / sqrt(2)) - 1
           else 2 * stats::pt(x / sqrt(2), nu) - 1)
  }
  if (is.infinite(nu)) return(rs$f(x))
  chi_lo <- stats::qchisq(1e-13, nu)
  chi_hi <- stats::qchisq(1e-13, nu, lower.tail = FALSE)
  ccstar <- (rs$wsat / x)^2 * nu
  upper  <- min(ccstar, chi_hi)
  below  <- if (upper > chi_lo)
    stats::integrate(function(cc) rs$f(x * sqrt(cc / nu)) * stats::dchisq(cc, nu),
                     chi_lo, upper, rel.tol = 1e-11, abs.tol = 1e-15,
                     stop.on.error = FALSE)$value else 0
  tail <- if (ccstar < chi_hi)
    stats::pchisq(ccstar, nu, lower.tail = FALSE) else 0
  min(1, max(0, below + tail))
}

# Exact quantile: the root of the exact CDF, on the log scale, bracketed by the
# k = 2 quantile (below) and the Bonferroni bound (above).
.sr_quantile_exact_one <- function(p, k, nu) {
  if (nu > .SR_EXACT_UP_TO_DF) nu <- Inf  # as in .sr_cdf_exact_one()
  a <- 1 - p
  t_quantile <- function(prob)
    if (is.infinite(nu)) stats::qnorm(prob) else stats::qt(prob, nu)
  lo <- sqrt(2) * t_quantile(1 - a / 2)
  if (k == 2) return(lo)
  hi <- sqrt(2) * t_quantile(1 - a / (k * (k - 1)))
  rs <- .sr_range_spline(k)
  f <- function(lq) .sr_cdf_exact_one(exp(lq), k, nu, rs) - p
  # The bracket holds in exact arithmetic; extendInt only guards against the
  # spline's last digits at probabilities within about 1e-9 of 0 or 1.
  exp(stats::uniroot(f, c(log(lo), log(hi)), tol = 1e-11,
                     extendInt = "upX")$root)
}

# Values of stats::ptukey() or stats::qtukey() (fun) at cells inside the
# domain, with NaN for every cell at which the function warned. The two warn
# only where they fail: qtukey() with "convergence failed in 'qtukey'" and
# NaN, ptukey() with NaN or with "full precision may not have been achieved"
# beside a value it does not vouch for. A NaN cell is recomputed exactly by the
# caller, so the warning would describe a value that is never returned, and it
# is muffled here. The vectorized call comes first; only if it warned are the
# cells evaluated one at a time to find which failed, so that the others keep
# R's value. The handler records that it fired in an environment made for it.
.sr_r_values <- function(fun, x, k, nu, ...) {
  quiet <- function(x, k, nu) {
    state <- new.env(parent = emptyenv())
    state$warned <- FALSE
    val <- withCallingHandlers(fun(x, k, nu, ...), warning = function(w) {
      state$warned <- TRUE
      invokeRestart("muffleWarning")
    })
    if (state$warned) NA_real_ * val else val
  }
  val <- quiet(x, k, nu)
  if (!anyNA(val) || length(x) == 1L) return(val)
  mapply(quiet, x, k, nu, USE.NAMES = FALSE)
}

# Distribution function of Q, vectorized with recycling like stats::ptukey():
# the exact computation where .sr_needs_exact(), stats::ptukey() elsewhere,
# and the exact computation again wherever stats::ptukey() fails inside the
# domain (k >= 2, nu > 0; see .sr_r_values()) or, with more than 50 means,
# where its lower-tail probability is below .75 (.80 with more than 1000
# means; .sr_needs_exact_p()).
# Outside the domain the call goes to stats::ptukey(), which returns NaN with
# its own warning.
.sr_cdf <- function(q, k, nu, lower_tail = TRUE) {
  n <- max(length(q), length(k), length(nu))
  q <- rep_len(q, n); k <- rep_len(k, n); nu <- rep_len(nu, n)
  out <- rep(NA_real_, n)
  ok  <- !(is.na(q) | is.na(k) | is.na(nu))
  dom <- ok & k >= 2 & nu > 0
  ex  <- dom & .sr_needs_exact(k, nu)
  rd  <- dom & !ex
  rn  <- ok & !dom
  if (any(rd)) {
    out[rd] <- .sr_r_values(stats::ptukey, q[rd], k[rd], nu[rd],
                            lower.tail = lower_tail)
    cum <- if (lower_tail) out else 1 - out
    ex <- ex | (rd & (!is.finite(out) |
                        .sr_needs_exact_p(cum, k, tails = FALSE)))
  }
  if (any(rn))
    out[rn] <- stats::ptukey(q[rn], nmeans = k[rn], df = nu[rn],
                             lower.tail = lower_tail)
  for (i in which(ex)) {
    p <- .sr_cdf_exact_one(q[i], k[i], nu[i])
    out[i] <- if (lower_tail) p else 1 - p
  }
  out
}

# Quantile function of Q (lower-tail probability p), vectorized with recycling:
# the root of the exact distribution function where .sr_needs_exact() or
# .sr_needs_exact_p(), stats::qtukey() elsewhere, unless the check
# (.sr_qtukey_offset()) finds qtukey()'s value more than 2e-6 from ptukey()'s
# root. stats::qtukey() also fails to converge, and returns NaN with a
# warning, at lower probabilities when k is large (for example p = .25 with
# k = 51 and nu = 32, or p = .05 with k = 20 and nu = 20000); inside the
# domain (0 < p < 1, k >= 2, nu > 0) the exact root is used there as well.
# Outside the domain the call goes to stats::qtukey(), which returns NaN with
# its own warning.
.sr_quantile <- function(p, k, nu) {
  n <- max(length(p), length(k), length(nu))
  p <- rep_len(p, n); k <- rep_len(k, n); nu <- rep_len(nu, n)
  out <- rep(NA_real_, n)
  ok  <- !(is.na(p) | is.na(k) | is.na(nu))
  dom <- ok & p > 0 & p < 1 & k >= 2 & nu > 0
  ex  <- dom & (.sr_needs_exact(k, nu) | .sr_needs_exact_p(p, k))
  rd  <- dom & !ex
  rn  <- ok & !dom
  if (any(rd)) {
    out[rd] <- .sr_r_values(stats::qtukey, p[rd], k[rd], nu[rd])
    ex <- ex | (rd & !is.finite(out))
    chk <- which(rd & is.finite(out))
    if (length(chk)) {
      off <- .sr_qtukey_offset(out[chk], p[chk], k[chk], nu[chk])
      ex[chk] <- !is.finite(off) | abs(off) > .SR_QTUKEY_CHECK_TOL
    }
  }
  if (any(rn))
    out[rn] <- stats::qtukey(p[rn], nmeans = k[rn], df = nu[rn])
  # A root find costs tens of milliseconds, and callers often ask for the same
  # quantile many times (contrast_test() asks once per contrast), so each
  # distinct (p, k, nu) is solved once.
  if (any(ex)) {
    key <- sprintf("%.17g|%.17g|%.17g", p[ex], k[ex], nu[ex])
    first <- !duplicated(key)
    idx <- which(ex)[first]
    val <- vapply(idx, function(i) .sr_quantile_exact_one(p[i], k[i], nu[i]),
                  numeric(1))
    out[ex] <- val[match(key, key[first])]
  }
  out
}
