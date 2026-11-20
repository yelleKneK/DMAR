# Shared numerics of ci_ratio_of_means() and ci_ratio_of_means_paired().
#
# Both functions estimate a ratio num / den of two means (or of two
# advantages over a baseline) and report one of four sets for it: the
# Fieller confidence set, in closed form when the degrees of freedom are
# fixed (.fieller_set()) and by root finding when they change with the
# ratio, as under unequal variances (.fieller_set_df()); the log-scale
# interval of the response ratio (.log_ratio_limits()); and the
# percentile and BCa bootstrap intervals (.ratio_boot_means(),
# .ratio_boot_totals(), .ratio_boot_keep(), the jackknife accelerations,
# and .ratio_boot_limits()). The exported functions validate their
# inputs, form the means and their covariance matrix, call these, and
# assemble the table (.ratio_table()). The footer line, tidy(), and
# glance() of the shared class dmar_ci_ratio_of_means live here too.

# ---- The Fieller set ------------------------------------------------------
#
# Let num and den estimate mu_num and mu_den with covariance matrix V
# (V[1, 1] = var(num), V[2, 2] = var(den), V[1, 2] = cov) and let
#   T(rho) = (num - rho den) / sqrt(V11 - 2 rho V12 + rho^2 V22),
# the t statistic of the linear hypothesis mu_num - rho mu_den = 0
# (Fieller, 1954, p. 176). The confidence set is every rho the test does
# not reject: {rho : -q_upper <= T(rho) <= q_lower}, where q_lower and
# q_upper are the critical values of the two tail areas (Inf for a tail
# of zero, a one-sided bound).
#
# Negating num and den together leaves the ratio and V unchanged, so the
# denominator is oriented positive first. Then T(est) = 0 at the estimate
# est = num / den, T > 0 below it and T < 0 above it, and the set splits
# at est into a part below, governed by q_lower alone, and a part above,
# governed by q_upper alone. On each side the constraint is
# T(rho)^2 <= q^2, that is, the quadratic
#   Q_q(rho) = A rho^2 + B rho + C <= 0, with
#   A = den^2 - q^2 V22, B = -2 (num den - q^2 V12), C = num^2 - q^2 V11.
# Its shape is decided by A, which is positive exactly when the
# denominator's t, den / sqrt(V22), exceeds q (after the orientation den is
# positive, so this is the absolute value of the t of the denominator as
# the data give it):
#   A > 0          {Q <= 0} is the bounded interval between the roots;
#   A < 0, D > 0   it is the line outside the roots;
#   A < 0, D <= 0  it is the whole line;
#   A = 0          it is a half-line (an edge case of probability zero),
# with D = B^2 - 4 A C. When A > 0, D >= 0 in exact arithmetic, because
# Q_q(est) = -q^2 var(num - est den) <= 0 puts est between the roots, so
# a negative D there is rounding and is clamped at zero; when A <= 0, D
# is used as computed. With equal tail areas the two sides give Fieller's
# three cases: a bounded interval, the complement of an interval, or the
# whole line. With unequal tail areas a side can be bounded while the
# other is not, which yields an interval unbounded on one side or, when
# one side is a complement, the union of two intervals.
.fieller_set <- function(num, den, V, q_lower, q_upper) {
  if (den < 0) {
    num <- -num
    den <- -den
  }
  est <- num / den
  pieces <- rbind(.fieller_side(num, den, V, q_lower, est, below = TRUE),
                  .fieller_side(num, den, V, q_upper, est, below = FALSE))
  .ratio_set_summary(pieces)
}

# The part of {Q_q <= 0} on one side of the estimate, as a two-column
# matrix of closed pieces (from, to). The estimate always belongs to the
# set; the clamps keep a rounding error in a root from detaching it.
.fieller_side <- function(num, den, V, q, est, below) {
  half <- if (below) c(-Inf, est) else c(est, Inf)
  if (is.infinite(q)) return(matrix(half, ncol = 2L))
  A <- den^2 - q^2 * V[2L, 2L]
  B <- -2 * (num * den - q^2 * V[1L, 2L])
  C <- num^2 - q^2 * V[1L, 1L]
  pieces <- if (A > 0) {
    r <- .quadratic_roots(A, B, C, max(B^2 - 4 * A * C, 0))
    matrix(c(min(r[1L], est), max(r[2L], est)), ncol = 2L)
  } else if (A < 0) {
    D <- B^2 - 4 * A * C
    if (D <= 0) {
      matrix(c(-Inf, Inf), ncol = 2L)
    } else {
      r <- .quadratic_roots(A, B, C, D)
      if (est <= mean(r)) {
        rbind(c(-Inf, max(r[1L], est)), c(r[2L], Inf))
      } else {
        rbind(c(-Inf, r[1L]), c(min(r[2L], est), Inf))
      }
    }
  } else if (B > 0) {
    matrix(c(-Inf, max(-C / B, est)), ncol = 2L)
  } else if (B < 0) {
    matrix(c(min(-C / B, est), Inf), ncol = 2L)
  } else {
    matrix(c(-Inf, Inf), ncol = 2L)
  }
  from <- pmax(pieces[, 1L], half[1L])
  to <- pmin(pieces[, 2L], half[2L])
  keep <- from <= to
  cbind(from[keep], to[keep])
}

# The two real roots of A r^2 + B r + C, given the discriminant D >= 0,
# in increasing order. The form that avoids subtracting nearly equal
# numbers (B and sqrt(D)) keeps both roots accurate when B^2 dominates
# 4 A C.
.quadratic_roots <- function(A, B, C, D) {
  h <- -0.5 * (B + if (B >= 0) sqrt(D) else -sqrt(D))
  if (h == 0) return(c(0, 0))
  sort(c(h / A, C / h))
}

# Reads a set given as closed pieces: merges pieces that touch, then
# describes the result by its two outer limits and, when it has a gap,
# the gap's ends. The set is every value from `lower` to `upper` except
# those strictly between `excluded_lower` and `excluded_upper`.
.ratio_set_summary <- function(pieces) {
  pieces <- pieces[order(pieces[, 1L]), , drop = FALSE]
  merged <- pieces[1L, , drop = FALSE]
  for (k in seq_len(nrow(pieces))[-1L]) {
    last <- nrow(merged)
    if (pieces[k, 1L] <= merged[last, 2L]) {
      merged[last, 2L] <- max(merged[last, 2L], pieces[k, 2L])
    } else {
      merged <- rbind(merged, pieces[k, ])
    }
  }
  if (nrow(merged) > 2L) {
    .dmar_refuse(c(
      "The confidence set came out in more than two pieces, which a ",
      "ratio of means cannot produce; please report this with the data."),
      "internal_error")
  }
  lower <- merged[1L, 1L]
  upper <- merged[nrow(merged), 2L]
  unbounded <- is.infinite(lower) && is.infinite(upper)
  if (nrow(merged) == 1L) {
    excluded <- c(NA_real_, NA_real_)
    shape <- if (unbounded) "whole line" else "interval"
  } else {
    excluded <- c(merged[1L, 2L], merged[2L, 1L])
    shape <- if (unbounded) "complement of an interval" else
      "union of two intervals"
  }
  list(lower = lower, upper = upper, excluded_lower = excluded[1L],
       excluded_upper = excluded[2L], shape = shape)
}

# The Fieller set when the degrees of freedom of the test change with the
# hypothesized ratio, as they do for independent groups with unequal
# variances: the test of mu_num - rho mu_den = 0 is a Welch-type test of
# the contrast m_1 - rho m_2 - (1 - rho) m_b, whose Satterthwaite degrees
# of freedom depend on rho. Inverting that test at every rho gives the
# set that agrees with the test, so that each finite limit L is a value
# at which the test's p-value equals its tail area (for two groups,
# stats::t.test(x, L * y) has p-value alpha_lower + alpha_upper when the
# tails are equal).
#
# `m` holds the group means (numerator group, denominator group, and the
# baseline group if there is one), `v` the variances of those means, and
# `nu` the degrees of freedom of each variance. No closed form exists, so
# the roots are found on the angle theta = atan(rho), which maps the whole
# line onto [-pi/2, pi/2]: with c = cos(theta) and s = sin(theta), the
# contrast c (m_1 - m_b) - s (m_2 - m_b) has the sign and the t statistic
# of num - rho den for every finite rho, and at theta = -pi/2 and pi/2 its
# t statistic is the denominator's t and its negative, so the ends of the
# line are evaluated rather than approached. The sides below and above the
# estimate are scanned on a grid that includes the point at which the t
# statistic is extreme, and then the point at which each side's boundary
# function g peaks, since a narrow excluded gap sits around that peak; it
# is near the peak of t but not at it, because the critical value moves
# with the degrees of freedom. Every sign change is refined by
# stats::uniroot(). Returns the set as .ratio_set_summary() reads it, with
# the degrees of freedom of the test at each finite end of the set
# (lower, upper, excluded_lower, excluded_upper; NA at an infinite end
# and where there is no gap).
.fieller_set_df <- function(m, v, nu, alpha_lower, alpha_upper,
                            grid = 2000L) {
  k <- length(m)
  base <- if (k == 3L) m[3L] else 0
  if (m[2L] - base < 0) {
    m <- -m
    base <- -base
  }
  num <- m[1L] - base
  den <- m[2L] - base
  est <- num / den
  th_hat <- atan(est)
  weights <- function(th) {
    w <- cbind(cos(th), -sin(th))
    if (k == 3L) w <- cbind(w, sin(th) - cos(th))
    w
  }
  t_df <- function(th) {
    w <- weights(th)
    vc <- drop(w^2 %*% v)
    list(t = drop(w %*% m) / sqrt(vc),
         df = vc^2 / drop(w^4 %*% (v^2 / nu)))
  }
  # Each side's set is {g <= 0}.
  g_below <- function(th) {
    s <- t_df(th)
    s$t - .ratio_crit(alpha_lower, s$df)
  }
  g_above <- function(th) {
    s <- t_df(th)
    -s$t - .ratio_crit(alpha_upper, s$df)
  }
  # The stationary point of the t statistic in rho, where |t| is largest
  # (Fieller, 1954, p. 177); its angle joins the grid of its side.
  V12 <- if (k == 3L) v[3L] else 0
  V11 <- v[1L] + V12
  V22 <- v[2L] + V12
  th_star <- atan((den * V11 - num * V12) / (den * V12 - num * V22))
  side <- function(g, alpha, from, to) {
    if (alpha == 0) return(matrix(c(from, to), ncol = 2L))
    th <- seq(from, to, length.out = grid)
    if (is.finite(th_star) && th_star > from && th_star < to) {
      th <- sort(c(th, th_star))
    }
    gv <- g(th)
    # g can peak above zero between two grid points at which it is
    # negative, and the gap around that peak would then be missed. Since g
    # rises to its peak and falls from it, the peak lies within one
    # spacing of the largest grid value; it is also near th_star, the peak
    # of t. Each of those two grid points is bracketed by its neighbors,
    # the maximum of g in the bracket is found by stats::optimize(), and
    # the maximizer joins the grid.
    n_th <- length(th)
    for (j in unique(stats::na.omit(c(which.max(gv), match(th_star, th))))) {
      bracket <- th[c(max(j - 1L, 1L), min(j + 1L, n_th))]
      top <- stats::optimize(g, bracket, maximum = TRUE, tol = 1e-12)
      if (!top$maximum %in% th) {
        th <- c(th, top$maximum)
        gv <- c(gv, top$objective)
      }
    }
    o <- order(th)
    th <- th[o]
    gv <- gv[o]
    roots <- th[gv == 0]
    change <- which(gv[-length(gv)] * gv[-1L] < 0)
    for (j in change) {
      roots <- c(roots, stats::uniroot(g, c(th[j], th[j + 1L]),
                                       tol = 1e-14, maxiter = 1000L)$root)
    }
    b <- sort(unique(c(from, roots, to)))
    mid <- (b[-1L] + b[-length(b)]) / 2
    keep <- g(mid) <= 0
    cbind(b[-length(b)][keep], b[-1L][keep])
  }
  pieces <- rbind(side(g_below, alpha_lower, -pi / 2, th_hat),
                  side(g_above, alpha_upper, th_hat, pi / 2))
  to_ratio <- function(th) {
    ifelse(th <= -pi / 2, -Inf, ifelse(th >= pi / 2, Inf,
                                       ifelse(th == th_hat, est, tan(th))))
  }
  out <- .ratio_set_summary(matrix(to_ratio(pieces), ncol = 2L))
  df_at <- function(r) if (is.finite(r)) t_df(atan(r))$df else NA_real_
  out$df_limits <- c(lower = df_at(out$lower), upper = df_at(out$upper),
                     excluded_lower = df_at(out$excluded_lower),
                     excluded_upper = df_at(out$excluded_upper))
  out
}

# The Satterthwaite degrees of freedom of a linear combination with
# weights w of independent means whose variances are v, each on nu
# degrees of freedom.
.satterthwaite_df <- function(w, v, nu) {
  sum(w^2 * v)^2 / sum(w^4 * v^2 / nu)
}

# The critical value of a tail area on df degrees of freedom; a tail of
# zero (a one-sided bound) has an infinite critical value.
.ratio_crit <- function(alpha, df) {
  if (alpha == 0) Inf else stats::qt(1 - alpha, df)
}

# ---- The log ratio ------------------------------------------------------
#
# The interval of the response ratio (Hedges, Gurevitch, & Curtis, 1999):
# log(num / den) plus or minus a critical value times its delta method
# standard error,
#   se^2 = V11 / num^2 + V22 / den^2 - 2 V12 / (num den),
# mapped back by exp(). Defined only when both num and den are positive,
# which the caller checks.
.log_ratio_limits <- function(num, den, V, q_lower, q_upper) {
  se <- sqrt(V[1L, 1L] / num^2 + V[2L, 2L] / den^2 -
               2 * V[1L, 2L] / (num * den))
  log_est <- log(num / den)
  c(lower = exp(log_est - q_lower * se), upper = exp(log_est + q_upper * se))
}

# The classed refusal of the log interval for a ratio whose terms are not
# both positive. `what` names the two terms ("means" or "advantages over
# the baseline").
.log_ratio_refuse <- function(num, den, what) {
  .dmar_refuse(c(
    "ci_method = \"log\" forms the interval for log(ratio), which is defined ",
    "only when both ", what, " are positive (here ",
    format(num, digits = 4), " and ", format(den, digits = 4), "). Use ",
    "ci_method = \"fieller\", which needs no sign condition."),
    "log_ratio_undefined")
}

# ---- The bootstrap --------------------------------------------------------
#
# Written in base R, as cohen_kappa() is, so that the bootstrap runs
# without the boot package. Replicates are drawn in blocks: each block
# draws n * k indices at once and averages the columns of an n by k
# matrix, which is far faster than one call per replicate and holds at
# most about a million indices in memory.

# B bootstrap means of v, resampled from v alone (one group of an
# independent groups design).
.ratio_boot_means <- function(v, B) {
  n <- length(v)
  per <- max(1L, floor(1e6 / n))
  out <- numeric(B)
  done <- 0L
  while (done < B) {
    k <- min(per, B - done)
    idx <- sample.int(n, n * k, replace = TRUE)
    out[done + seq_len(k)] <- colMeans(matrix(v[idx], nrow = n, ncol = k))
    done <- done + k
  }
  out
}

# B bootstrap replicates of the ratio sum(a) / sum(b) when the pairs
# (a_i, b_i) are resampled together: the per-unit values of a paired
# design, or the cluster totals of a clustered one. The ratio of the
# means of the resampled units equals the ratio of their sums, because
# both means divide by the same number of units. Returns the replicate
# ratios and their denominators.
.ratio_boot_totals <- function(a, b, B) {
  n <- length(a)
  per <- max(1L, floor(1e6 / n))
  num <- den <- numeric(B)
  done <- 0L
  while (done < B) {
    k <- min(per, B - done)
    idx <- sample.int(n, n * k, replace = TRUE)
    at <- done + seq_len(k)
    num[at] <- colSums(matrix(a[idx], nrow = n, ncol = k))
    den[at] <- colSums(matrix(b[idx], nrow = n, ncol = k))
    done <- done + k
  }
  list(ratio = num / den, den = den)
}

# What the rule of R/dmar_bootstrap.R needs to know about a bootstrap
# interval with the tail areas alpha_lower and alpha_upper: the kind
# ("quantile" for the percentile interval, "bca"), the smallest tail area
# a reported limit cuts off (for a limit at quantile position p,
# min(p, 1 - p), over the limits asked for), and how the result is named
# in a message: an interval with equal tails, two limits with tail areas
# of their own, or a one-sided bound (a tail area of zero on the other
# side) with the quantile it sits at.
.ratio_boot_rule <- function(alpha_lower, alpha_upper, conf_level, method) {
  areas <- c(alpha_lower, alpha_upper)
  asked <- areas > 0
  tail <- min(pmin(areas[asked], 1 - areas[asked]))
  what <- if (!all(asked)) "bound" else
    if (isTRUE(all.equal(alpha_lower, alpha_upper))) "interval" else "limits"
  position <- if (what != "bound") NULL else
    if (alpha_lower > 0) alpha_lower else 1 - alpha_upper
  list(kind = if (method == "bca") "bca" else "quantile", tail = tail,
       what = what, position = position, conf_level = conf_level)
}

# The check of B for a bootstrap interval of a ratio (.dmar_check_B()):
# 10 replications beyond each limit for the percentile interval, 50
# beyond each nominal limit for BCa, at the smallest tail area asked for.
.ratio_check_B <- function(B, rule) {
  .dmar_check_B(B, kind = rule$kind, tail = rule$tail,
                conf_level = rule$conf_level, what = rule$what,
                position = rule$position)
}

# Drops the replicates whose ratio is undefined (a denominator of exactly
# zero, which discrete data can produce) and returns the kept draws with
# the verdict of the rule on them (.dmar_boot_settle(), not yet signaled:
# the BCa position rule may still withhold the interval). Through DMAR
# 1.0.x the call stopped when fewer than 100 replicates remained; the
# limits are now NA under a warning of class dmar_bootstrap_insufficient.
.ratio_boot_keep <- function(ratio, den, B, rule, fun) {
  keep <- is.finite(ratio) & den != 0
  verdict <- .dmar_boot_settle(
    B, sum(keep), kind = rule$kind, tail = rule$tail,
    conf_level = rule$conf_level, fun = fun, quantity = "the ratio",
    dropped_because = paste0("had a denominator of zero, so their ratio ",
                             "was undefined,"),
    what = rule$what, position = rule$position, signal = FALSE)
  list(draws = ratio[keep], verdict = verdict)
}

# The BCa acceleration from the jackknife,
# a = sum(l^3) / (6 (sum(l^2))^(3/2)), where l_i = (n - 1)(est - est_(-i))
# are the jackknife influence values centered at the full-sample estimate,
# as boot's usual jackknife (boot:::usual.jack, behind
# boot::empinf(type = "jack")) centers them.
# Efron and Tibshirani (1993, Equation 14.15, p. 186) center the
# leave-one-out values at their own mean instead; the two differ at a higher
# order (for the ratio of the means of their Table 14.1 data, .0100678
# against .0100692). For a paired or clustered design the units deleted are
# the pairs or the clusters, whose totals are a and b.
.ratio_acceleration_totals <- function(a, b) {
  n <- length(a)
  est <- sum(a) / sum(b)
  loo <- (sum(a) - a) / (sum(b) - b)
  if (!all(is.finite(loo))) .ratio_jackknife_refuse()
  l <- (n - 1) * (est - loo)
  .ratio_acceleration(list(l), n)
}

# The same for independent groups, whose bootstrap resamples each group
# from itself: the jackknife deletes one observation at a time within its
# group, as boot::empinf(type = "jack") does for stratified resampling,
# and the influence values of group g enter the acceleration scaled by
# its size n_g,
#   a = sum_g n_g^-3 sum_i l_gi^3 / (6 (sum_g n_g^-2 sum_i l_gi^2)^(3/2)),
# the skewness of the linear approximation to the bootstrap distribution
# of the ratio. With equal group sizes the scaling cancels and a is the
# unstratified formula applied to all the influence values.
# boot::boot.ci() applies that unscaled formula to stratified influence
# values whatever the group sizes, so the two agree only when the groups
# are the same size.
.ratio_acceleration_groups <- function(x, y, z = NULL) {
  loo <- function(v) (sum(v) - v) / (length(v) - 1L)
  mx <- mean(x)
  my <- mean(y)
  mz <- if (is.null(z)) 0 else mean(z)
  est <- (mx - mz) / (my - mz)
  jack <- list((loo(x) - mz) / (my - mz), (mx - mz) / (loo(y) - mz))
  n <- c(length(x), length(y))
  if (!is.null(z)) {
    jack <- c(jack, list((mx - loo(z)) / (my - loo(z))))
    n <- c(n, length(z))
  }
  if (!all(is.finite(unlist(jack)))) .ratio_jackknife_refuse()
  l <- Map(function(t, ng) (ng - 1) * (est - t), jack, n)
  .ratio_acceleration(l, n)
}

.ratio_acceleration <- function(l, n) {
  num <- sum(mapply(function(li, ng) sum(li^3) / ng^3, l, n))
  den <- sum(mapply(function(li, ng) sum(li^2) / ng^2, l, n))
  if (den == 0) return(0)
  num / (6 * den^1.5)
}

.ratio_jackknife_refuse <- function() {
  .dmar_refuse_argument(c(
    "Deleting one observation leaves a denominator of zero, so the ",
    "jackknife that estimates the BCa acceleration is undefined for ",
    "these data. Use ci_method = \"percentile\"."), arg = NULL)
}

# Percentile and BCa limits from the kept replicates of the estimate
# `est`. `kept` is the list .ratio_boot_keep() returns. A tail area of
# zero gives the one-sided bound and an infinite limit on that side. The
# BCa levels are
#   Phi(z0 + (z0 + z_a) / (1 - acc (z0 + z_a))), z0 = Phi^-1(P(boot < est)),
# the bias correction taken from the bootstrap distribution, as
# boot::boot.ci() takes it. `acc` is a function that returns the
# acceleration, called only when a BCa interval is to be formed. Returns
# list(limits, verdict) after signaling the conditions of the verdict: a
# limit that was asked for is NA when the rule withholds the interval
# (more than half of the replicates dropped, or for BCa a B that puts
# fewer than 10 replicates beyond a limit where the adjustment places it).
.ratio_boot_limits <- function(kept, est, alpha_lower, alpha_upper,
                               method, acc = function() 0) {
  boot <- kept$draws
  verdict <- kept$verdict
  probs <- c(alpha_lower, 1 - alpha_upper)
  asked <- c(alpha_lower > 0, alpha_upper > 0)
  lim <- c(-Inf, Inf)
  lim[asked] <- NA_real_
  if (verdict$formed && method == "bca") {
    p0 <- mean(boot < est)
    if (p0 == 0 || p0 == 1) {
      .dmar_refuse_argument(c(
        "Every bootstrap replicate fell on one side of the estimate, so ",
        "the BCa bias correction is infinite. Use ci_method = ",
        "\"percentile\"."), arg = NULL)
    }
    z0 <- stats::qnorm(p0)
    a <- acc()
    open <- probs > 0 & probs < 1
    z <- stats::qnorm(probs[open])
    probs[open] <- stats::pnorm(z0 + (z0 + z) / (1 - a * (z0 + z)))
    # The count of 10 beyond each limit, where the adjustment places it;
    # a side with no limit asked for sets no requirement.
    verdict <- .dmar_boot_bca_positions(
      verdict, c(if (asked[1L]) probs[1L] else 0.5,
                 if (asked[2L]) probs[2L] else 0.5))
  }
  if (verdict$formed) {
    if (asked[1L]) lim[1L] <- stats::quantile(boot, probs[1L], names = FALSE)
    if (asked[2L]) lim[2L] <- stats::quantile(boot, probs[2L], names = FALSE)
  }
  .dmar_boot_signal(verdict)
  list(limits = lim, verdict = verdict)
}

# ---- Validation -------------------------------------------------------

# A numeric vector of observations with at least two values; NAs are
# allowed here and handled by the caller.
.ratio_check_data <- function(v, name) {
  if (!is.numeric(v) || !is.null(dim(v))) {
    .dmar_refuse_argument(c(
      "'", name, "' must be a numeric vector."), arg = NULL)
  }
  if (any(is.infinite(v))) {
    .dmar_refuse_argument(c(
      "'", name, "' contains an infinite value."), arg = NULL)
  }
  invisible(v)
}

.ratio_check_scalar <- function(v, name, positive = FALSE) {
  if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
      (positive && v <= 0)) {
    .dmar_refuse_argument(c(
      "'", name, "' must be a single finite number",
      if (positive) " greater than zero", "."), arg = NULL)
  }
  invisible(v)
}

.ratio_check_n <- function(v, name) {
  .ratio_check_scalar(v, name)
  if (v < 2 || v != round(v)) {
    .dmar_refuse_argument(c(
      "'", name, "' must be a whole number of at least 2."), arg = NULL)
  }
  invisible(v)
}

.ratio_bootstrap_refuse <- function(method) {
  .dmar_refuse(c(
    "ci_method = \"", method, "\" resamples the observations, so it needs ",
    "the raw data in 'x' and 'y'; from summary statistics use ",
    "ci_method = \"fieller\" or \"log\"."), "bootstrap_needs_data")
}

# ---- The table, its footer, and its tidiers -------------------------------

# Assembles the long table both functions return: the limits and the
# estimate, the gap of a set that has one, then the rows the caller
# supplies in `rows` (a named numeric vector), as a dmar_ci_long table of
# class dmar_ci_ratio_of_means. `set` is a .ratio_set_summary() list.
.ratio_table <- function(set, ratio, rows, conf_level, alpha_lower,
                         alpha_upper, ci_method, attrs) {
  head <- c(ci_lower = set$lower, ratio = ratio, ci_upper = set$upper,
            excluded_lower = set$excluded_lower,
            excluded_upper = set$excluded_upper)
  value <- c(head, rows)
  out <- data.frame(term = names(value), value = unname(value),
                    stringsAsFactors = FALSE)
  out <- .as_dmar_tbl(out, conf_level = conf_level,
                      subclass = "dmar_ci_long",
                      tail_areas = .tail_areas(alpha_lower, alpha_upper))
  class(out) <- c("dmar_ci_ratio_of_means", class(out))
  # A bootstrap interval withheld by the rule of R/dmar_bootstrap.R has no
  # limits and so no form: the attribute is NA, and the footer says that
  # no interval is reported.
  withheld <- is.na(set$lower) && is.na(set$upper) &&
    isTRUE(attrs$bootstrap_status %in%
             c("insufficient_count", "insufficient_share"))
  attr(out, "confidence_set") <- if (withheld) NA_character_ else set$shape
  attr(out, "ci_method") <- ci_method
  for (nm in names(attrs)) attr(out, nm) <- attrs[[nm]]
  out
}

# The sentence that names the form of the confidence set. It reads the
# rows of the table, so it states the gap of a complement in the same
# numbers the table prints.
.ratio_set_sentence <- function(x) {
  shape <- attr(x, "confidence_set", exact = TRUE)
  method <- attr(x, "ci_method", exact = TRUE)
  if (is.null(shape) || is.null(method) || is.null(x$term)) return(NULL)
  v <- function(term) x$value[x$term == term][1L]
  lo <- v("ci_lower")
  hi <- v("ci_upper")
  label <- switch(method,
    fieller = "Fieller confidence set",
    log = "Log ratio interval",
    percentile = "Percentile bootstrap interval",
    bca = "BCa bootstrap interval")
  # A withheld bootstrap interval: the reason is in the bootstrap line.
  if (is.na(shape)) return(paste0(label, ": no interval is reported."))
  form <- switch(shape,
    "interval" = if (is.finite(lo) && is.finite(hi)) {
      "a bounded interval"
    } else if (is.finite(lo)) {
      "an interval unbounded above"
    } else if (is.finite(hi)) {
      "an interval unbounded below"
    } else {
      "the whole line"
    },
    "complement of an interval" = paste0(
      "the complement of an interval, every value except those strictly ",
      "between excluded_lower and excluded_upper"),
    "whole line" = "the whole line; no value of the ratio is excluded at this level",
    "union of two intervals" = paste0(
      "the union of two intervals, from ci_lower to excluded_lower and ",
      "from excluded_upper to ci_upper"),
    shape)
  line <- paste0(label, ": ", form, ".")
  # A side is open when its limit is infinite although its tail area is
  # not zero; a zero tail asks for a one-sided bound, unbounded on that side
  # by design. Every open side, whatever the form of the set, is explained:
  # an interval unbounded above can come from a small upper tail area as
  # well as from a one-sided bound, and the reader must be able to tell
  # the two apart.
  tails <- .recorded_tail_areas(x)
  open <- .ratio_open_sides(lo, hi, if (is.null(tails)) 1 else tails[[1L]],
                            if (is.null(tails)) 1 else tails[[2L]])
  if (identical(method, "fieller") && any(open)) {
    line <- paste0(line, " The set is unbounded", if (!all(open)) {
      if (open[["lower"]]) " below" else " above"
    }, " because |t_denominator| does not exceed the critical value",
    if (!is.null(tails) && !all(open)) {
      paste0(" of the ", if (open[["lower"]]) "lower" else "upper",
             " tail area")
    }, ".")
  }
  fieller <- attr(x, "fieller_set", exact = TRUE)
  if (!identical(method, "fieller") && is.character(fieller)) {
    line <- paste0(line, " The Fieller set at this level is ", fieller,
                   "; this interval is bounded regardless and can ",
                   "undercover when the denominator is near zero.")
  }
  line
}

# The counts behind a bootstrap interval (.dmar_boot_footer(): requested,
# used, dropped, and required, with the reason when no interval is
# reported and the mark of an exploratory run).
.ratio_boot_sentence <- function(x) {
  B <- attr(x, "B", exact = TRUE)
  used <- attr(x, "B_used", exact = TRUE)
  if (is.null(B) || is.null(used)) return(NULL)
  line <- .dmar_boot_footer(x)
  if (!is.null(line)) {
    if (used < B) line <- paste0(
      line, if (!grepl("\\.$", line)) ".",
      " A dropped replication had a denominator of zero")
    return(line)
  }
  line <- sprintf("Bootstrap: %d replications, %d used", B, used)
  if (used < B) {
    line <- paste0(line, sprintf(
      " (%d had a denominator of zero and were dropped)", B - used))
  }
  line
}

# The form of the Fieller set at the same tail areas, for the note a log
# or bootstrap interval carries when that set is unbounded on a side whose
# tail area is not zero (a zero tail asks for a one-sided bound, which is
# unbounded on that side by design). NULL when the Fieller set is bounded
# wherever a limit was asked for.
.ratio_fieller_note <- function(set, alpha_lower, alpha_upper) {
  open <- any(.ratio_open_sides(set$lower, set$upper, alpha_lower,
                                alpha_upper)) ||
    !is.na(set$excluded_lower)
  if (!open) return(NULL)
  switch(set$shape,
    "interval" = "an interval unbounded on one side",
    "complement of an interval" = "the complement of an interval",
    "whole line" = "the whole line",
    "the union of two intervals")
}

# The sides of a set that are unbounded although their tail area is not
# zero, c(lower = , upper = ). Such a side is unbounded exactly when the
# absolute value of the denominator's t does not exceed that side's
# critical value. A gap always comes with one, since a side whose critical
# value |t| exceeds is a bounded piece with no gap.
.ratio_open_sides <- function(lower, upper, alpha_lower, alpha_upper) {
  c(lower = alpha_lower > 0 && is.infinite(lower),
    upper = alpha_upper > 0 && is.infinite(upper))
}

# The lines under the table: those every table carries (the confidence
# level and any unequal tail areas), the form of the confidence set, and
# for a bootstrap interval the replication counts.
#' @exportS3Method .dmar_tbl_footer dmar_ci_ratio_of_means
#' @noRd
.dmar_tbl_footer.dmar_ci_ratio_of_means <- function(x, ...) {
  lines <- .dmar_footer_add(NextMethod(), .ratio_set_sentence(x),
                            "confidence_set")
  .dmar_footer_add(lines, .ratio_boot_sentence(x), "bootstrap")
}

# The broom verbs for the two ratio of means functions. The dmar_ci_long
# row (the ratio, its limits, and the level) is kept, and the ends of an
# excluded gap and the form of the set follow the limits, so that a
# complement of an interval is not read from tidy() as the whole line.
#' @rdname dmar_tidiers
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.dmar_ci_ratio_of_means <- function(x, ...) {
  out <- tidy.dmar_ci_long(x, ...)
  v <- function(term) {
    r <- x$value[x$term == term]
    if (length(r)) r[1L] else NA_real_
  }
  out$excluded_lower <- v("excluded_lower")
  out$excluded_upper <- v("excluded_upper")
  out$confidence_set <- attr(x, "confidence_set", exact = TRUE) %||%
    NA_character_
  first <- c("term", "estimate", "ci_lower", "ci_upper", "excluded_lower",
             "excluded_upper", "confidence_set")
  out[c(first, setdiff(names(out), first))]
}

#' @rdname dmar_tidiers
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.dmar_ci_ratio_of_means <- function(x, ...) {
  tidy.dmar_ci_ratio_of_means(x, ...)
}
