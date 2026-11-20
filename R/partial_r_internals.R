# Internal helpers shared by the partial and semipartial correlation family:
# var_partial_r(), var_semipartial_r(), expected_partial_r(),
# ss_aipe_partial_r(), ss_aipe_semipartial_r(), and the two sensitivity
# functions. None is exported.

# .check_partial_p(p, scalar = TRUE)
#
# p is the number of predictors in the regression of Y on X and the
# variables partialed out, the same quantity in every function of the
# family, so p - 1 variables are partialed out and p must be a whole number
# of at least 2. (DMAR 1.0.x took the count of variables partialed out, J,
# with p = J + 1; that name is not accepted.) In the notation of a model
# comparison q predictors are removed from the full model to give the
# restricted model, p - q variables are partialed out, and q = 1 for the
# correlation of one predictor. The simple correlation (p = 1, nothing partialed out) has
# functions of its own (var_r(), expected_r(), ss_aipe_r(),
# ss_aipe_r_sensitivity()). With scalar = FALSE a vector is accepted
# (expected_partial_r() is vectorized over p) and every element is checked;
# a missing element passes and gives a missing row. Returns p unchanged.
#
# Inside the family's functions and in the helpers below the arithmetic is
# written in J = p - 1, the number of variables partialed out, as it was in
# 1.0.x: each exported function sets J <- p - 1 once, after
# .check_partial_p(), and the helpers below take that J.
.check_partial_p <- function(p, scalar = TRUE) {
  msg <- paste0(
    "'p' (the number of predictors in the regression of Y on X and the ",
    "variables partialed out) must be a ",
    if (scalar) "single " else "", "whole number of at least 2. With ",
    "p = 1 nothing is partialed out and the correlation is the simple ",
    "correlation, which var_r(), expected_r(), ss_aipe_r(), and ",
    "ss_aipe_r_sensitivity() cover.")
  # The refusal has class dmar_argument_invalid and carries arg = "p" and
  # the value given.
  if (scalar) {
    if (!is.numeric(p) || length(p) != 1L || !is.finite(p) ||
        p != round(p) || p < 2) .dmar_refuse_argument(msg, arg = "p", value = p)
  } else {
    if (!is.numeric(p) || any(p != round(p) | p < 2, na.rm = TRUE))
      .dmar_refuse_argument(msg, arg = "p", value = p)
  }
  p
}

# .semipartial_r_var_numerator(r_sp, R2_full)
#
# First-order normal-theory variance of the sample semipartial correlation,
# times n: with q = R2_full - r_sp^2, the squared multiple correlation of
# the outcome with the variables partialed out,
#   n Var(r_sp) -> (1 - r_sp^2)^2 - q (1 - 2 r_sp^2).
# The number of variables partialed out does not enter at this order. Both
# pmax() calls guard a sample pair: a rounded pair can have q slightly
# below 0, and the expression is 0 at its lower boundary (r_sp = 0,
# R2_full = 1).
.semipartial_r_var_numerator <- function(r_sp, R2_full) {
  q <- pmax(R2_full - r_sp^2, 0)
  pmax((1 - r_sp^2)^2 - q * (1 - 2 * r_sp^2), 0)
}

# .semipartial_r_var_numerator_max(r_sp)
#
# The largest value of that numerator over every full-model R2 that is
# possible with r_sp (r_sp^2 <= R2_full <= 1), for a plan made without a
# planning value of R2_full. The numerator is linear in q: it decreases in
# q when r_sp^2 < 1/2, so its largest value is (1 - r_sp^2)^2 at q = 0, and
# increases when r_sp^2 > 1/2, so its supremum is r_sp^2 (1 - r_sp^2) as q
# approaches 1 - r_sp^2. J is at least 1 throughout the family, so q > 0 is
# always possible.
.semipartial_r_var_numerator_max <- function(r_sp) {
  max((1 - r_sp^2)^2, r_sp^2 * (1 - r_sp^2))
}

# .check_semipartial_R2_full(r_sp, R2_full, arg = "R2_full", r_arg = "r_sp")
#
# R2_full must be a single number in [0, 1] and cannot be below r_sp^2,
# because the squared semipartial correlation is the part of the full-model
# R2 that X adds. A pair that differs by less than 1e-8 passes, since
# R2_full = r_sp^2 computed in floating point may fall a few ulps short.
#
# A pair that two-decimal reports can produce also passes, with a warning.
# A reported R2_full stands for a value up to .005 above it and a reported
# r_sp for a value up to .005 nearer 0, so the pair is consistent with
# rounding when R2_full + .005 >= max(0, |r_sp| - .005)^2 (r_sp = .46 with
# R2_full = .21: .46^2 = .2116, and .215 >= .455^2 = .207). The callers
# then proceed as they did in DMAR 1.0.x. Where the result is computed from
# q = max(R2_full - r_sp^2, 0) (the first-order variance, the planner, the
# sensitivity function; uses = "q_zero") that is q = 0, variables
# partialed out that explain none of the outcome, and the warning says so.
# The squared standard error of the t test (var_semipartial_r() with
# under_null = TRUE; uses = "as_given") is (1 - R2_full) / (N - p - 1),
# which does not involve r_sp, so there R2_full is used as given and the
# warning says that instead. The warning is classed
# (dmar_semipartial_r_rounded) and carries the two values.
#
# The refusal of any other incoherent pair is classed
# (dmar_semipartial_r_incoherent, a dmar_refusal) and carries the two
# values.
.check_semipartial_R2_full <- function(r_sp, R2_full, arg = "R2_full",
                                       r_arg = "r_sp", uses = "q_zero") {
  if (!is.numeric(R2_full) || length(R2_full) != 1L || is.na(R2_full) ||
      R2_full < 0 || R2_full > 1)
    .dmar_refuse_argument(c(
      "'", arg, "' must be a single number in [0, 1]."), arg = NULL)
  if (R2_full < r_sp^2 - 1e-8) {
    if (R2_full + 0.005 >= max(0, abs(r_sp) - 0.005)^2) {
      .dmar_warn(c(
        "'", arg, "' (", format(R2_full), ") is below '", r_arg,
        "' squared (", format(r_sp^2), ") by no more than rounding to two ",
        "decimals allows. The squared semipartial correlation cannot ",
        "exceed the full-model R2, so ",
        if (identical(uses, "as_given")) {
          paste0("one of the two was rounded. The value returned does not ",
                 "depend on ", r_arg, " and uses ", arg, " as given.")
        } else {
          paste0("the variables partialed out are taken to explain none ",
                 "of the outcome (", arg, " = ", r_arg, "^2).")
        }),
        "semipartial_r_rounded", R2_full = R2_full, r_sp = r_sp)
    } else {
      .dmar_refuse(c(
        "'", arg, "' (", format(R2_full), ") is smaller than '", r_arg,
        "' squared (", format(r_sp^2), "). The squared semipartial ",
        "correlation is the part of the full-model R2 that X adds, so '", arg,
        "' cannot be below it. If the two values were rounded, use ", arg,
        " = ", r_arg, "^2."),
        "semipartial_r_incoherent", R2_full = R2_full, r_sp = r_sp)
    }
  }
  invisible(R2_full)
}

# ---- The assurance step of ss_aipe_partial_r() --------------------------

# .partial_r_cdf(c, m, rho)
#
# P(r <= c) for the sample correlation r of m bivariate normal observations
# with population correlation rho, and so, with m = n - J, for a sample
# partial correlation given J variables from n observations.
#
# With U chi square on m - 1 degrees of freedom, V chi square on m - 2, and
# Z standard normal, all independent, and theta = rho / sqrt(1 - rho^2),
#   r / sqrt(1 - r^2) = (Z + theta sqrt(U)) / sqrt(V)
# in distribution (regress Y on X: given X the slope is normal about its
# population value, U is the sum of squares of X over its variance, and V
# the residual sum of squares over the residual variance). Put S = U + V,
# chi square on 2m - 3 degrees of freedom, and B = U / S, which is
# Beta((m - 1)/2, (m - 2)/2) and independent of S. With
# a = c / sqrt(1 - c^2), the event r <= c is
#   Z <= sqrt(S) {a sqrt(1 - B) - theta sqrt(B)},
# and for a fixed B the probability of that is a central t probability,
# because Z / sqrt(S / (2m - 3)) is Student's t on 2m - 3 degrees of
# freedom. So
#   P(r <= c) = integral over b in (0, 1) of
#     pt(sqrt(2m - 3) {a sqrt(1 - b) - theta sqrt(b)}, 2m - 3) *
#     dbeta(b, (m - 1)/2, (m - 2)/2) db.
# Only the central t and a beta density enter. The better known form, a
# noncentral t probability integrated over a chi square, is not used: R's
# pt() with ncp is a normal approximation above a noncentrality of 37.62,
# which plans for narrow intervals reach (at m = 400 and rho = .95 it was
# off by .00037, 4.6 standard errors of a simulation of 40,000,000
# correlations, where this form was within 0.8).
#
# The integral runs between the 1e-14 and 1 - 1e-14 quantiles of the beta
# density, which leaves out at most 2e-14 of probability. m must be at
# least 3 (2.5 or more works numerically; the callers pass whole numbers).
# At rho = 0 the value equals pt(sqrt(m - 2) a, m - 2).
.partial_r_cdf <- function(c, m, rho) {
  if (c <= -1) return(0)
  if (c >= 1) return(1)
  a <- c / sqrt(1 - c^2)
  theta <- rho / sqrt(1 - rho^2)
  k <- 2 * m - 3
  shape_1 <- (m - 1) / 2
  shape_2 <- (m - 2) / 2
  integrand <- function(b)
    stats::pt(sqrt(k) * (a * sqrt(1 - b) - theta * sqrt(b)), k) *
      stats::dbeta(b, shape_1, shape_2)
  value <- stats::integrate(
    integrand,
    lower = stats::qbeta(1e-14, shape_1, shape_2),
    upper = stats::qbeta(1 - 1e-14, shape_1, shape_2),
    rel.tol = 1e-11, subdivisions = 2000L)$value
  min(max(value, 0), 1)
}

# .partial_r_half_width(r, n, J, z, fisher_z)
#
# Half the width of the interval ss_aipe_partial_r() plans, at a sample
# partial correlation r: z (1 - r^2) / sqrt(n - J - 1) on the raw scale, and
# half the distance between the back-transformed limits of Fisher's Z
# interval, whose standard error is 1 / sqrt(n - J - 3), otherwise. Both
# decrease as |r| grows.
.partial_r_half_width <- function(r, n, J, z, fisher_z) {
  if (!fisher_z) return(z * (1 - r^2) / sqrt(n - J - 1))
  half_z <- z / sqrt(n - J - 3)
  (tanh(atanh(r) + half_z) - tanh(atanh(r) - half_z)) / 2
}

# .partial_r_width_probability(n, rho, J, half_width, z, fisher_z)
#
# The probability, at sample size n and population partial correlation
# rho, that the realized half-width is no larger than half_width. Because
# the half-width decreases in |r|, that event is |r| >= r*, with r* the
# |r| at which the half-width equals the target (closed form on the raw
# scale, a root on Fisher's), and the probability is
# 1 - F(r*) + F(-r*) with F = .partial_r_cdf(., n - J, rho). When the
# half-width at r = 0, the widest the interval can be, already meets the
# target, the probability is 1.
.partial_r_width_probability <- function(n, rho, J, half_width, z,
                                         fisher_z) {
  if (.partial_r_half_width(0, n, J, z, fisher_z) <= half_width) return(1)
  r_star <- if (!fisher_z) {
    sqrt(1 - half_width * sqrt(n - J - 1) / z)
  } else {
    stats::uniroot(
      function(r) .partial_r_half_width(r, n, J, z, TRUE) - half_width,
      lower = 0, upper = 1 - 1e-12, tol = 1e-13)$root
  }
  m <- n - J
  1 - .partial_r_cdf(r_star, m, rho) + .partial_r_cdf(-r_star, m, rho)
}

# .partial_r_assurance_n(n_0, rho, J, half_width, z, fisher_z, assurance)
#
# The smallest sample size, at or above n_0 (the plan without assurance)
# and the smallest n the interval exists for (J + 3 on the raw scale,
# J + 4 on Fisher's), at which .partial_r_width_probability() is at least
# `assurance`. The search does not go below n_0: below it the probability
# is not monotone in n (adjacent sample sizes differed by as much as .165
# in a scan), so a smaller n that happens to satisfy the request is not a
# plan; at or above n_0 the probability did not decrease from one n to the
# next in any of the 72 planning sets scanned over every integer n. The
# step upward doubles until the request is met, and a bisection then finds
# the first n that meets it. Returns list(n, attained).
.partial_r_assurance_n <- function(n_0, rho, J, half_width, z, fisher_z,
                                   assurance) {
  prob <- function(n)
    .partial_r_width_probability(n, rho, J, half_width, z, fisher_z)
  lo <- max(n_0, J + if (fisher_z) 4 else 3)
  p_lo <- prob(lo)
  if (p_lo >= assurance) return(list(n = lo, attained = p_lo))
  step <- 1
  repeat {
    hi <- lo + step
    p_hi <- prob(hi)
    if (p_hi >= assurance) break
    lo <- hi
    step <- 2 * step
    if (hi > 1e9)
      .dmar_refuse_argument(c(
        "The sample size search for 'assurance' passed 1e9 without ",
        "reaching the requested probability."), arg = "assurance")
  }
  while (hi - lo > 1) {
    mid <- lo + (hi - lo) %/% 2
    p_mid <- prob(mid)
    if (p_mid >= assurance) {
      hi <- mid
      p_hi <- p_mid
    } else {
      lo <- mid
    }
  }
  list(n = hi, attained = p_hi)
}

# ---- The assurance step of ss_aipe_semipartial_r() ----------------------

# .semipartial_r_half_width_draws(G, n, r_sp, q, J, z)
#
# G independent draws of the half-width of the interval
# ss_aipe_semipartial_r() plans, z sqrt(V-hat / (n - J - 1)), for n
# multivariate normal observations from a population with semipartial
# correlation r_sp and full-model R2 r_sp^2 + q. V-hat is
# .semipartial_r_var_numerator() at the sample semipartial correlation and
# the sample full-model R2, the same function and clamping
# ss_aipe_semipartial_r_sensitivity() applies to each simulated data set.
#
# The population is the canonical one (the comment in
# ss_aipe_semipartial_r_sensitivity() gives the argument that nothing is
# lost): Z_1, ..., Z_J, X, and an error E are independent standard normal
# and
#   Y = r_sp X + sqrt(q) Z_1 + e E,   e = sqrt(1 - r_sp^2 - q),
# which also covers a full-model R2 of 1 (e = 0). The matrix A of sums of
# squares and products of (Z_1, ..., Z_J, X, E) about their means is
# Wishart with n - 1 degrees of freedom and identity scale, and the sums
# of squares and products that involve Y follow from the three
# coefficients. A is not formed. Its Cholesky factor L (A = L L', lower
# triangular) has independent entries, chi (not chi square) on
# n - 1, n - 2, ... degrees of freedom down the diagonal and standard
# normal below it (the Bartlett decomposition), and with the Z's first the
# quantities needed are functions of a few of those entries. Write
# l_x = L[X, X], chi on n - 1 - J degrees of freedom; l_ex = L[E, X],
# standard normal; l_e = L[E, E], chi on n - 2 - J; and l_1 = L[Z_1, Z_1],
# chi on n - 1. Then, with residuals taken on the Z's,
#   sum of squares of the residual of X          l_x^2
#   cross product of that residual with Y        l_x (r_sp l_x + e l_ex)
#   residual sum of squares of Y on the Z's      (r_sp l_x + e l_ex)^2 +
#                                                e^2 l_e^2
#   total sum of squares of Y                    the line above, plus
#     (sqrt(q) l_1 + s N)^2 + s^2 W,
# where s = sqrt(r_sp^2 + e^2), N is standard normal, and W is chi square
# on J - 1 degrees of freedom (the entries of the rows of X and E in the
# columns of the Z's enter only through r_sp L[X, Z_j] + e L[E, Z_j],
# which is s times a standard normal for each j). The sample semipartial
# correlation is (r_sp l_x + e l_ex) / sqrt(total), and the sample
# full-model R2 is its square plus 1 - residual / total. This is the
# distribution of the statistics computed from stats::rWishart(G, n - 1,
# diag(J + 2)); a test forms A = L L' from the same entries and confirms
# the algebra against the matrix expressions. n must be at least J + 3.
.semipartial_r_half_width_draws <- function(G, n, r_sp, q, J, z) {
  e <- sqrt(max(0, 1 - r_sp^2 - q))
  s <- sqrt(r_sp^2 + e^2)
  l_1  <- sqrt(stats::rchisq(G, n - 1))
  l_x  <- sqrt(stats::rchisq(G, n - 1 - J))
  l_e  <- sqrt(stats::rchisq(G, n - 2 - J))
  l_ex <- stats::rnorm(G)
  n_1  <- stats::rnorm(G)
  w    <- if (J > 1) stats::rchisq(G, J - 1) else 0
  .semipartial_r_half_width_from_bartlett(l_1, l_x, l_e, l_ex, n_1, w,
                                          r_sp, q, n, J, z)
}

# The arithmetic of .semipartial_r_half_width_draws(), apart from the
# draws, so that a test can pass entries of a Cholesky factor it built.
.semipartial_r_half_width_from_bartlett <- function(l_1, l_x, l_e, l_ex,
                                                    n_1, w, r_sp, q, n, J,
                                                    z) {
  e <- sqrt(max(0, 1 - r_sp^2 - q))
  s <- sqrt(r_sp^2 + e^2)
  y_x <- r_sp * l_x + e * l_ex
  ss_resid <- y_x^2 + e^2 * l_e^2
  ss_total <- ss_resid + (sqrt(q) * l_1 + s * n_1)^2 + s^2 * w
  r_hat <- y_x / sqrt(ss_total)
  R2_hat <- r_hat^2 + 1 - ss_resid / ss_total
  z * sqrt(.semipartial_r_var_numerator(r_hat, R2_hat) / (n - J - 1))
}

# .semipartial_r_assurance_n(n_0, r_sp, q, J, half_width, z, assurance, G)
#
# The Monte Carlo assurance step of ss_aipe_semipartial_r(): the smallest
# sample size found, at or above n_0 (the plan without assurance) and
# J + 3, at which the share of G simulated half-widths no larger than
# half_width is at least `assurance`. The step upward doubles until the
# share reaches the request and a bisection follows, with fresh draws at
# every sample size tried. The share at the returned n is then estimated
# again from G draws that had no part in choosing n, because the share
# that chose n is biased upward by the choice. Returns
# list(n, attained, mc_se), the last two from the fresh draws. The caller
# sets the seed.
.semipartial_r_assurance_n <- function(n_0, r_sp, q, J, half_width, z,
                                       assurance, G) {
  share <- function(n)
    mean(.semipartial_r_half_width_draws(G, n, r_sp, q, J, z) <= half_width)
  lo <- max(n_0, J + 3)
  if (share(lo) >= assurance) {
    hi <- lo
  } else {
    step <- 1
    repeat {
      hi <- lo + step
      if (share(hi) >= assurance) break
      lo <- hi
      step <- 2 * step
      if (hi > 1e9)
        .dmar_refuse_argument(c(
          "The sample size search for 'assurance' passed 1e9 without ",
          "reaching the requested probability."), arg = "assurance")
    }
    while (hi - lo > 1) {
      mid <- lo + (hi - lo) %/% 2
      if (share(mid) >= assurance) hi <- mid else lo <- mid
    }
  }
  attained <- share(hi)
  list(n = hi, attained = attained,
       mc_se = sqrt(attained * (1 - attained) / G))
}
