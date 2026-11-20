# Helpers of ci_sc(), ss_aipe_sc(), ss_aipe_sc_sensitivity(),
# ss_aipe_sc_ancova(), and ss_aipe_sc_ancova_sensitivity(). None is
# exported.

# .sc_check_c_weights(c_weights): the weights rule of ci_sc(), with its
# messages, shared by ci_sc(), the planner of its interval, and that
# planner's sensitivity function. The weights are judged by their values:
# integer weights such as -1:1 pass as c(-1, 0, 1) does, and the positive
# values must sum to 1 within 1e-8, so the order of the weights does not
# matter (1.0.x compared the sum with 1 by identical(), which refused
# c(.7, .2, .1, -1), whose sum is 1 - 1.1e-16, and accepted c(.1, .2, .7, -1)).
# Weights with no nonzero value are refused first, with the message of
# .dmar_check_contrast_nonzero() (R/dmar_conditions.R), rather than with
# the message about the positive values, which they also break.
.sc_check_c_weights <- function(c_weights) {
  if (!is.numeric(c_weights)) .dmar_refuse_argument(
    "'c_weights' must be a numeric vector of contrast weights, such as c(1/2, 1/2, -1).",
    arg = "c_weights")
  .dmar_check_contrast_nonzero(c_weights)
  w <- as.numeric(c_weights)
  if (!isTRUE(abs(sum(w[w > 0]) - 1) <= 1e-8)) .dmar_refuse_argument(
    "The positive values of 'c_weights' must sum to 1; give the weights as fractions, such as c(1/2, 1/2, -1).",
    arg = "c_weights")
  if (!identical(round(sum(w), 5), 0)) .dmar_refuse_argument(
    "The sum of the contrast weights ('c_weights') should equal zero.",
    arg = "c_weights")
  invisible(c_weights)
}

# .sc_width_probability(n, delta, c_weights, width, alpha_lower, alpha_upper)
#
# The exact probability, with n per group and the standardized contrast
# delta, that the ci_sc() interval with tail areas alpha_lower and
# alpha_upper is no wider than `width`. The observed t is noncentral t on
# J(n - 1) df with noncentrality delta / se, se = sqrt(sum(c^2) / n), and
# the width is se times the distance between the ci_nc_t() limits at t.
# At t = 0 the limits are -qnorm(1 - alpha_lower) and qnorm(1 - alpha_upper)
# (P(T <= 0) is pnorm(-ncp)), and both have slope E(sqrt(V / df)) there, so
# the distance is smallest at t = 0 and rises on either side. The event is
# t_lo <= t <= t_hi, and its probability is one less the two noncentral t
# tails, each exact (.nc_t_prob()); 0 when the interval at t = 0 is already
# wider than `width`. The width at t = 0 is computed from the normal
# quantiles, not from ci_nc_t(), which finds them only to its tolerance: at
# delta = 0, c(1, -1), width 1e-7, and n = 3,073,167,056,555,294, where the
# interval at t = 0 is 9.5e-16 of its width too wide, the limits of
# ci_nc_t() made it 3.3e-22 narrower than `width`, and the probability came
# out 1 where it is 0.
#
# Precision: with many degrees of freedom the distance changes slowly with
# t near t_lo and t_hi (its slope there is about (z_l + z_u) t / (2 df)),
# so a small error in the ci_nc_t() limits moves t_lo and t_hi much more.
# At delta = 0, c(1, -1), width .03, and n = 34148 (df 68294), the width
# from the ci_nc_t() limits at t = 2.6072304 is 3.3e-12 below the width
# from an oracle that conditions on the normal variate, so t_hi comes out
# 6.0e-6 high and the probability .99087039 against .99087022 (1.6e-7). A
# smaller `tol` for ci_nc_t() does not change this (the error is that of
# pt() there), and it is far below the error of any planning value.
.sc_width_probability <- function(n, delta, c_weights, width, alpha_lower,
                                  alpha_upper) {
  se <- sqrt(sum(c_weights^2) / n)
  df <- length(c_weights) * (n - 1)
  z_tails <- stats::qnorm(alpha_lower, lower.tail = FALSE) +
    stats::qnorm(alpha_upper, lower.tail = FALSE)
  over <- function(t) {
    if (t == 0) return(se * z_tails - width)
    lim <- ci_nc_t(ncp = t, df = df, conf_level = NULL,
                   alpha_lower = alpha_lower, alpha_upper = alpha_upper,
                   verbose = FALSE)
    se * (lim$value[lim$term == "ci_upper"] -
            lim$value[lim$term == "ci_lower"]) - width
  }
  if (over(0) > 0) return(0)
  ncp <- delta / se
  edge <- function(direction) {
    b <- direction * max(1, abs(ncp))
    doublings <- 0L
    while (over(b) <= 0) {
      b <- 2 * b
      doublings <- doublings + 1L
      if (doublings > 60L)
        .dmar_refuse(c(
          "The search for the observed t at which the interval reaches ",
          "'width' did not end; the target may be infeasible."),
          "computation_failed")
    }
    stats::uniroot(over, sort(c(0, b)), tol = 1e-10)$root
  }
  t_lo <- edge(-1)
  t_hi <- edge(1)
  p <- 1 - .nc_t_prob(t_lo, df, ncp) -
    .nc_t_prob(t_hi, df, ncp, lower_tail = FALSE)
  min(max(p, 0), 1)
}

# ---- ss_aipe_sc_ancova(): the plans of Lai and Kelley (2012) -----------
# delta is the standardized contrast on the scale of the divisor and
# `ratio` = sigma_ancova / sigma_anova = sqrt(1 - rho^2); the planner's
# noncentrality is delta / (ratio sqrt(sum(c^2) / n)), and ratio = 1 gives
# the divisor s_ancova (their Equations 17 and 24). The helpers serve both
# divisors.

# The plan without assurance, for the criterion that the interval at
# t = lambda is no wider than `width` (Equations 21 and 28): the smallest n
# per group at which it holds, searched from n_min = 2, the smallest n with
# a positive error df J(n - 1) - 1. The width at the planning value falls
# as n grows from n_min (checked at every n from n_min to 80 and in steps
# of 10 to 600, for 420 cases of contrast, ratio, weights for two to five
# groups, and level), so the search of .sc_smallest_n() returns the
# smallest n that meets it. The interval at t = 0 is the narrowest, z_l +
# z_u wide on the scale of t, so no n below n_0 = C ratio^2 (z_l + z_u)^2 /
# width^2 meets the criterion, and a width whose n_0 is past n_max is
# refused before the search, as is one the search does not reach by n_max.
.sc_ancova_plan_width <- function(delta, ratio, c_weights, width,
                                  alpha_lower, alpha_upper) {
  J <- length(c_weights)
  C <- sum(c_weights^2)
  n_min <- ceiling(1 + 2 / J)
  n_max <- 2^53
  # At delta = 0 the width is computed from the normal quantiles, as in
  # ss_aipe_sc() (see .sc_width_probability()).
  z_tails <- stats::qnorm(alpha_lower, lower.tail = FALSE) +
    stats::qnorm(alpha_upper, lower.tail = FALSE)
  n_0 <- C * ratio^2 * z_tails^2 / width^2
  if (!is.finite(n_0) || n_0 > n_max) .sc_width_unreachable(width, n_max)
  width_at <- function(n) {
    if (delta == 0) return(ratio * sqrt(C / n) * z_tails)
    lim <- ci_nc_t(ncp = delta / (ratio * sqrt(C / n)), df = n * J - J - 1,
                   conf_level = NULL, alpha_lower = alpha_lower,
                   alpha_upper = alpha_upper)
    upper <- lim[which(lim$term == "ci_upper"), 2] * ratio * sqrt(C / n)
    lower <- lim[which(lim$term == "ci_lower"), 2] * ratio * sqrt(C / n)
    abs(upper - lower)
  }
  .sc_smallest_n(n_min, function(m) width_at(m) - width <= 0, n_max,
                 function() .sc_width_unreachable(width, n_max))
}

# .sc_smallest_n(n, meets, n_max, unreachable)
#
# The smallest whole number at or above n at which meets() is TRUE, for a
# meets() that stays TRUE as the sample size grows (the interval at the
# planning value narrows as n grows). A step upward from n is doubled until
# meets() holds and the last step is then bisected, so the result is the
# one that raising n by one until meets() holds would give, from about
# 2 log2(m - n) evaluations in place of m - n (one ci_nc_t() call each, so
# raising n by one took minutes for a plan of millions per group). When
# no n up to n_max meets it, unreachable() is called, which refuses.
.sc_smallest_n <- function(n, meets, n_max, unreachable) {
  if (meets(n)) return(n)
  lo <- n
  step <- 1
  repeat {
    hi <- lo + step
    if (hi >= n_max) {
      hi <- n_max
      if (!meets(hi)) unreachable()
      break
    }
    if (meets(hi)) break
    lo <- hi
    step <- 2 * step
  }
  while (hi - lo > 1) {
    mid <- lo + (hi - lo) %/% 2
    if (meets(mid)) hi <- mid else lo <- mid
  }
  hi
}

# The refusal of a width that the plan of ss_aipe_sc() or
# ss_aipe_sc_ancova() cannot reach within n_max per group.
.sc_width_unreachable <- function(width, n_max) {
  .dmar_refuse_argument(c(
    "'width' of ", format(width), " needs more than 2^53 per group (about ",
    "9.0e15, the largest number up to which R holds every whole number ",
    "exactly) at these planning values, so no plan is returned; give a ",
    "wider 'width'."),
    arg = "width", value = width, n_max = n_max)
}

# The published assurance step (pp. 358--359; Kelley & Rausch, 2006): t*
# with P(|T| <= t*) = assurance under the planner's model at the plan
# without assurance, carried to the scale of the contrast and planned for.
# The model leaves out the covariate means and the sampling of
# s_ancova / s_anova, so the plan can attain less than the request. The
# plan of assurance_method = "analytic", the default.
.sc_ancova_published_assurance <- function(delta, ratio, c_weights, width,
                                           alpha_lower, alpha_upper,
                                           assurance) {
  delta <- abs(delta)
  J <- length(c_weights)
  C <- sum(c_weights^2)
  n0 <- .sc_ancova_plan_width(delta, ratio, c_weights, width, alpha_lower,
                              alpha_upper)
  ncp0 <- delta / (ratio * sqrt(C / n0))
  lim_2 <- ci_nc_t(ncp = ncp0, df = n0 * J - J - 1, conf_level = NULL,
                   alpha_upper = (1 - assurance) / 2,
                   alpha_lower = (1 - assurance) / 2)
  limit_2_sided <- sqrt(C / n0) * ratio *
    lim_2[which(lim_2$term == "ci_upper"), 2]
  lim_1 <- ci_nc_t(ncp = ncp0, df = n0 * J - J - 1, conf_level = NULL,
                   alpha_upper = 1 - assurance, alpha_lower = 0)
  limit_1_sided <- sqrt(C / n0) * ratio *
    lim_1[which(lim_1$term == "ci_upper"), 2]
  determine_limit <- function(current_psi_limit, sample_size = n0) {
    se <- ratio * sqrt(C / sample_size)
    less <- .nc_t_prob(-current_psi_limit / se, J * sample_size - J - 1,
                       delta / se)
    greater <- .nc_t_prob(current_psi_limit / se, J * sample_size - J - 1,
                          delta / se, lower_tail = FALSE)
    (less + greater - (1 - assurance))^2
  }
  limit <- stats::optimize(f = determine_limit,
                           interval = c(limit_1_sided, limit_2_sided))$minimum
  .sc_ancova_plan_width(limit, ratio, c_weights, width, alpha_lower,
                        alpha_upper)
}

# ---- ss_aipe_sc_ancova(): the simulated assurance step -----------------
# The event is that the ci_sc_ancova() interval is no wider than `width`.
# Model: J groups of n, X ~ N(mu_X, sigma_X^2) in every group, Y = mu_j +
# beta X + e, e ~ N(0, sigma_e^2) (sigma_e = sigma_ancova). The interval
# is Delta(t) sqrt(V) (divisor s_ancova) or Delta(t) (s_ancova / s_anova)
# sqrt(V) (s_anova), with Delta(t) the distance between the ci_nc_t()
# limits at t on N - J - 1 df, V = C/n + (c'xbar)^2 / SSW_x, C = sum(c^2),
# and t = c'm_adj / (s_ancova sqrt(V)), m_adj the adjusted means. Five
# independent statistics determine the width:
#   z_x = c'xbar / (sigma_X sqrt(C/n)) ~ N(0, 1);
#   W = SSW_x / sigma_X^2 ~ chi square on N - J;
#   z_e = c'ebar / (sigma_e sqrt(C/n)) ~ N(0, 1);
#   z_b = (b_w - beta) sqrt(SSW_x) / sigma_e ~ N(0, 1), b_w the pooled
#         within-group slope;
#   Q = SSE / sigma_e^2 ~ chi square on N - J - 1.
# (Given X, z_e, z_b, and Q are functions of orthogonal projections of the
# errors and do not depend on X.) Then V = (C/n)(1 + z_x^2 / W),
# c'm_adj - c'mu = sigma_e sqrt(C/n) (z_e - z_b z_x / sqrt(W)), so
# t = (delta_e / sqrt(V) + z) / sqrt(Q / (N - J - 1)) with
# z = sqrt(C/n) (z_e - z_b z_x / sqrt(W)) / sqrt(V) and delta_e = c'mu /
# sigma_e; and SSW_y / sigma_e^2 = Q + (kappa sqrt(W) + z_b)^2 with
# kappa = sqrt(1 - ratio^2) / ratio. The chi squares are drawn as uniforms
# mapped by qchisq() at each n, so one set of draws serves every n tried
# (common random numbers). The tests rebuild ci_sc_ancova() widths of raw
# data from these statistics.
.sc_ancova_base_draws <- function(G) {
  list(z_x = stats::rnorm(G), u_w = stats::runif(G), z_e = stats::rnorm(G),
       z_b = stats::rnorm(G), u_q = stats::runif(G))
}

# The widths of the simulated ci_sc_ancova() intervals at n per group;
# delta_e is the standardized contrast on the scale of sigma_ancova.
.sc_ancova_widths <- function(base, n, c_weights, delta_e, divisor, ratio,
                              alpha_lower, alpha_upper) {
  J <- length(c_weights)
  C <- sum(c_weights^2)
  df_e <- J * (n - 1) - 1
  W <- stats::qchisq(base$u_w, J * (n - 1))
  Q <- stats::qchisq(base$u_q, df_e)
  V <- C / n * (1 + base$z_x^2 / W)
  z <- sqrt(C / n) * (base$z_e - base$z_b * base$z_x / sqrt(W)) / sqrt(V)
  t <- (delta_e / sqrt(V) + z) / sqrt(Q / df_e)
  width <- .sc_nc_distance(t, df_e, alpha_lower, alpha_upper)(t) * sqrt(V)
  if (divisor == "s_anova") {
    kappa <- sqrt(1 - ratio^2) / ratio
    ssw_y <- Q + (kappa * sqrt(W) + base$z_b)^2
    width <- width * sqrt((Q / df_e) / (ssw_y / (J * (n - 1))))
  }
  width
}

# The distance between the ci_nc_t() limits as a function of t over the
# range of `t`: ci_nc_t() at Chebyshev nodes and a cubic spline through the
# distance over sqrt(1 + t^2 / (2 df)), as .reg_coef_nc_delta() does (a
# copy, so the contrast planners do not depend on that file).
.sc_nc_distance <- function(t, df, alpha_lower, alpha_upper, nodes = 41L) {
  lo <- min(t)
  hi <- max(t)
  if (!(hi > lo)) {
    lo <- lo - 1
    hi <- hi + 1
  }
  x <- lo + (hi - lo) * (1 - cos(seq(0, pi, length.out = nodes))) / 2
  dist <- vapply(x, function(v) {
    lim <- ci_nc_t(ncp = v, df = df, conf_level = NULL,
                   alpha_lower = alpha_lower, alpha_upper = alpha_upper,
                   verbose = FALSE)
    lim$value[lim$term == "ci_upper"] - lim$value[lim$term == "ci_lower"]
  }, numeric(1))
  scale <- function(v) sqrt(1 + v^2 / (2 * df))
  f <- stats::splinefun(x, dist / scale(x), method = "fmm")
  function(v) f(v) * scale(v)
}

# The smallest n at or above n_0 with prob(n) >= assurance: doubling, then
# bisection, so the n below the one returned falls short (when above n_0).
.sc_assurance_search <- function(n_0, prob, assurance) {
  lo <- n_0
  p_lo <- prob(lo)
  if (p_lo >= assurance) return(list(n = lo, attained = p_lo))
  step <- 1
  repeat {
    hi <- lo + step
    p_hi <- prob(hi)
    if (p_hi >= assurance) break
    lo <- hi
    step <- 2 * step
    if (hi > 2^53)
      .dmar_refuse_argument(c(
        "The sample size search for 'assurance' passed 2^53 per group ",
        "without reaching the requested probability."), arg = "assurance")
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

# The simulated step of assurance_method = "simulation": the smallest n at
# or above n_0 at which the share of G simulated intervals no wider than
# `width` is at least `assurance`, from one set of draws. With equal tails
# the width depends on the sign only through t's reflection, so |delta_e|
# is used. Returns list(n, attained, mc_se). The caller sets the seed.
.sc_ancova_assurance_n <- function(n_0, c_weights, delta_e, divisor, ratio,
                                   width, alpha_lower, alpha_upper,
                                   assurance, G) {
  base <- .sc_ancova_base_draws(G)
  share <- function(n)
    mean(.sc_ancova_widths(base, n, c_weights, abs(delta_e), divisor, ratio,
                           alpha_lower, alpha_upper) <= width)
  step <- .sc_assurance_search(n_0, share, assurance)
  list(n = step$n, attained = step$attained,
       mc_se = sqrt(step$attained * (1 - step$attained) / G))
}
