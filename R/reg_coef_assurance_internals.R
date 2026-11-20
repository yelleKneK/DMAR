# The assurance steps of ss_aipe_reg_coef() and its wrappers ss_aipe_rc()
# and ss_aipe_src() that compute the probability that the realized
# confidence interval is no wider than the target: the smallest sample
# size, at or above the plan without assurance, at which it is at least
# `assurance`, exactly (assurance_method = "analytic" for the central
# interval in the units of the data) or by simulation (assurance_method =
# "simulation" for the intervals for the standardized coefficient). None
# is exported.
#
# Three intervals are planned for, and each planner names the one it
# targets:
#   * the central t interval for the coefficient in the units of the data
#     (ss_aipe_reg_coef(noncentral = FALSE), ss_aipe_rc()), whose width has
#     an exact F law (.reg_coef_raw_width_probability());
#   * the noncentral t interval for the standardized coefficient
#     (ss_aipe_reg_coef(noncentral = TRUE), ss_aipe_src()), simulated
#     (.reg_coef_std_assurance_n() with interval = "noncentral");
#   * the central t interval for the standardized coefficient, the interval
#     of Kelley and Maxwell (2003, Equation 1), whose standard error is
#     formed from the two sample squared multiple correlations
#     (ss_aipe_src(noncentral = FALSE)), simulated from the same draws
#     (interval = "central"). Its width varies less than that of the
#     interval in the units of the data, so the F law of the first
#     overstates the sample size it needs.
#
# Notation. N observations of Y and p predictors from a multivariate normal
# population; X_j is the targeted predictor and X_{-j} the other p - 1.
# R2 is the squared multiple correlation of Y with all p predictors, R2_j
# that of X_j with X_{-j}, and the analysis fits Y on all p predictors with
# an intercept. RSS is the residual sum of squares of Y and SSR_j the
# residual sum of squares of X_j on X_{-j} (with an intercept). The
# interval has N - p - 1 degrees of freedom.

# ---- The interval for the coefficient in the units of the data ---------

# .reg_coef_raw_width_probability(n, p, R2, R2_j, sd_ratio, width,
#                                 alpha_lower, alpha_upper)
#
# The probability, at sample size n, that the central t interval of
# ci_reg_coef() for the unstandardized coefficient (the interval
# ss_aipe_reg_coef() plans for when noncentral = FALSE) is no wider than
# `width`, with random predictors. sd_ratio is sd_Y / sd_X_j.
#
# The interval is b_j plus and minus a t quantile times the standard error
#   SE^2 = RSS / {(n - p - 1) SSR_j},
# so its width is (q_upper + q_lower) SE with q_upper = qt(1 - alpha_upper)
# and q_lower = qt(1 - alpha_lower) on n - p - 1 degrees of freedom (each
# computed in the upper tail, qt(alpha, lower.tail = FALSE), which keeps a
# small tail area that 1 - alpha would round away). Given
# the predictors, RSS / {sd_Y^2 (1 - R2)} is chi square on n - p - 1
# degrees of freedom, whatever the predictors are, so RSS is independent
# of the predictors and of SSR_j. SSR_j / {sd_X^2 (1 - R2_j)} is chi square
# on n - p degrees of freedom (the residual of a regression of X_j on the
# p - 1 others with an intercept; with p = 1 it is the sum of squares of
# X_j about its mean). So
#   SE^2 = K F / (n - p),  K = sd_ratio^2 (1 - R2) / (1 - R2_j),
# with F the ratio of the two chi squares, each over its degrees of
# freedom, an F variate on n - p - 1 and n - p degrees of freedom, and
#   P(width <= w) = P{F <= w^2 (n - p) / [(q_upper + q_lower)^2 K]}.
# With fixed predictors SSR_j is a constant and the same event is a chi
# square probability (RSS alone varies); the planner is for random
# predictors, the model ss_aipe_reg_coef_sensitivity() simulates.
.reg_coef_raw_width_probability <- function(n, p, R2, R2_j, sd_ratio, width,
                                            alpha_lower, alpha_upper) {
  nu <- n - p - 1
  q <- stats::qt(alpha_upper, nu, lower.tail = FALSE) +
    stats::qt(alpha_lower, nu, lower.tail = FALSE)
  K <- sd_ratio^2 * (1 - R2) / (1 - R2_j)
  stats::pf(width^2 * (n - p) / (q^2 * K), nu, n - p)
}

# .reg_coef_first_n(meets, lower, guess = lower, limit = Inf, beyond = NULL)
#
# The smallest whole n at or above `lower` at which meets(n) is TRUE, for a
# meets() that is FALSE below some n and TRUE from it on. Every search of
# the planners goes through it: the plans without assurance (lower = p + 2,
# the smallest sample with a residual degree of freedom, and `guess` near
# the answer), the analytic approximation for the noncentral interval,
# and the assurance steps that compute the probability
# (.reg_coef_assurance_search() below). From `guess` (raised to `lower`)
# the steps double, downward while meets() holds and upward while it does
# not, until meets() holds at the upper end of a bracket and not at the
# lower end (or the lower end is `lower`), and a bisection then finds the
# first n that meets. A guess a few cases from the answer costs a few
# evaluations, and the cost grows with the logarithm of the distance.
# beyond() is called, and must stop, when an upward step passes `limit`
# without meeting. DMAR 1.0.x and MBESS scanned upward one case at a time
# from a start below the answer (for the central plan, the larger of ten
# below the normal approximation and 2p + 1), and so returned the start
# when a smaller sample size already met the width.
.reg_coef_first_n <- function(meets, lower, guess = lower, limit = Inf,
                              beyond = NULL) {
  n <- max(lower, ceiling(guess))
  if (meets(n)) {
    hi <- n
    step <- 1
    repeat {
      if (hi == lower) return(hi)
      lo <- max(lower, hi - step)
      if (!meets(lo)) break
      hi <- lo
      step <- 2 * step
    }
  } else {
    lo <- n
    step <- 1
    repeat {
      hi <- lo + step
      if (meets(hi)) break
      if (hi > limit) beyond()
      lo <- hi
      step <- 2 * step
    }
  }
  while (hi - lo > 1) {
    mid <- lo + (hi - lo) %/% 2
    if (meets(mid)) hi <- mid else lo <- mid
  }
  hi
}

# .reg_coef_assurance_search(n_0, prob, assurance)
#
# The smallest whole n at or above n_0 at which prob(n) is at least
# `assurance` (.reg_coef_first_n() from n_0: the step upward doubles until
# the request is met, and a bisection then finds the first n that meets
# it, as .partial_r_assurance_n() does). The search does not go below n_0,
# the plan without assurance. Returns list(n, attained), with attained the
# value of prob() at that n. prob() is the exact probability for the
# central interval, and for the noncentral one the share of a fixed set
# of simulated samples (.reg_coef_std_assurance_n()); both are
# deterministic functions of n, and each value is computed once.
.reg_coef_assurance_search <- function(n_0, prob, assurance) {
  seen <- new.env(parent = emptyenv())
  key <- function(n) format(n, scientific = FALSE)
  meets <- function(n) {
    p_n <- prob(n)
    assign(key(n), p_n, envir = seen)
    p_n >= assurance
  }
  n <- .reg_coef_first_n(meets, n_0, limit = 1e9, beyond = function()
    .dmar_refuse_argument(c(
      "The sample size search for 'assurance' passed 1e9 without ",
      "reaching the requested probability."), arg = "assurance"))
  list(n = n, attained = get(key(n), envir = seen))
}

# ---- The two intervals for the standardized coefficient ----------------

# The interval ss_aipe_reg_coef() plans for with noncentral = TRUE (and so
# ss_aipe_src()) is the one ci_reg_coef(noncentral = TRUE) gives for
# standardized data: with t = b_j / SE the usual t statistic (it does not
# depend on the units) and
#   SE_std = sqrt{(1 - R2_hat) / [(1 - R2_j_hat)(n - p - 1)]},
# the limits are the limits of the noncentrality of a noncentral t on
# n - p - 1 degrees of freedom at the observed t (ci_nc_t()), times SE_std.
# The width, Delta(t) SE_std with Delta(t) the distance between the two
# noncentrality limits, depends on t, R2_hat, and R2_j_hat, and no closed
# form for its distribution is used: the probability is simulated from the
# exact sampling distribution of the sums of squares and products.
#
# The central interval ss_aipe_src(noncentral = FALSE) plans for is the one
# ci_reg_coef(noncentral = FALSE) gives for standardized data, Equation 1
# of Kelley and Maxwell (2003): the estimate plus and minus t quantiles on
# n - p - 1 degrees of freedom times the same SE_std. Its width is
# q SE_std with q = qt(1 - alpha_lower) + qt(1 - alpha_upper), the
# expression above with Delta(t) replaced by the constant q, so the same
# draws give its distribution. It depends on the sample through R2_hat and
# R2_j_hat alone (Kelley and Maxwell, 2003, p. 309).
#
# The population. In correlation units (unit variances of Y and X_j) write
#   X_j = a'Z + s_u U,   Y = beta s_u U + d'Z + s_e E,
# where Z (p - 1 standard normals) stands for X_{-j} whitened, U is the
# standardized residual of X_j on X_{-j}, E the standardized error, all
# independent; beta is the standardized coefficient, s_u^2 = 1 - R2_j,
# s_e^2 = 1 - R2, |a|^2 = R2_j, |d|^2 = R2 - beta^2 s_u^2 (the squared
# multiple correlation of Y with X_{-j}), and a'd = rho_YXj - beta s_u^2,
# with rho_YXj the correlation of Y with X_j. The distribution of the sums
# of squares and products of Z is unchanged by a rotation of Z, so a can
# be put on the first axis and d in the plane of the first two; the
# distribution of the width therefore depends on the population only
# through p, beta, R2, R2_j, and rho_YXj.
#
# The draws. The matrix of sums of squares and products of
# (Z_1, ..., Z_{p-1}, U, E) about their means is Wishart with n - 1 degrees
# of freedom and identity scale; its Cholesky factor L has independent
# entries, chi on n - 1, n - 2, ... degrees of freedom down the diagonal
# and standard normal below it (the Bartlett decomposition). Rows p (U) and
# p + 1 (E) and the first two rows of Z are all that the statistics need:
#   SSR_j = s_u^2 L[p, p]^2                      (residual of X_j)
#   RSS   = s_e^2 L[p+1, p+1]^2                  (residual of Y)
#   b_j   = beta + s_e L[p+1, p] / (s_u L[p, p]) (its coefficient)
#   SS_X  = (|a| L[1, 1] + s_u L[p, 1])^2 + s_u^2 sum_{m >= 2} L[p, m]^2
#   SS_Y  = sum over the columns m of L of
#           (d_1 L[1, m] + d_2 L[2, m] + beta s_u L[p, m] + s_e L[p+1, m])^2,
# where SS_X and SS_Y are the sums of squares of X_j and Y about their
# means and (d_1, d_2) are the coordinates of d. Then
# 1 - R2_hat = RSS / SS_Y and 1 - R2_j_hat = SSR_j / SS_X, and
#   SE_std = sqrt{RSS / [(n - p - 1) SSR_j]} sqrt(SS_X / SS_Y).
# A test forms the statistics from data by least squares and compares the
# distributions.
#
# The draws are made in two parts so that one set of samples serves every
# sample size the search tries (common random numbers).
# .reg_coef_std_base_draws(G, p) draws what does not depend on n: the
# standard normal entries, and, for each chi entry, a uniform number u
# whose chi square quantile at the degrees of freedom of a given n,
# qchisq(u, df), is that entry squared. .reg_coef_std_draws_at(base, n,
# beta, R2, R2_j) maps them to the statistics at sample size n. For a
# fixed u, qchisq(u, df) increases with df, so each simulated sample
# changes smoothly with n. In seven sets checked (p from 1 to 10), 3208 of
# 95 million comparisons of a simulated interval at n and n + 1, all with n
# below 60, were wider at n + 1, and the estimated probability never
# decreased from one n to the next. The search assumes that it does not,
# and where that holds it finds the smallest n at which the estimate
# reaches the request. With fresh draws at every n
# the estimates are independent from one n to the next, and a search can
# stop at an n where a lucky estimate reached the request though the next
# smaller or larger one did not. The population enters through
# .reg_coef_std_width() so that one set of draws serves every value of
# rho_YXj that is tried. .reg_coef_std_draws(G, n, p, beta, R2, R2_j) is
# the two in one, for one sample size.
.reg_coef_std_base_draws <- function(G, p) {
  base <- list(
    p = p,
    u_u = stats::runif(G),
    u_e = stats::runif(G),
    l_eu = stats::rnorm(G))
  if (p >= 2) {
    base$u_11 <- stats::runif(G)
    base$l_u1 <- stats::rnorm(G)
    base$l_e1 <- stats::rnorm(G)
  }
  if (p >= 3) {
    base$l_21 <- stats::rnorm(G)
    base$u_22 <- stats::runif(G)
    base$l_u2 <- stats::rnorm(G)
    base$l_e2 <- stats::rnorm(G)
  }
  if (p >= 4) {
    base$l_um <- matrix(stats::rnorm(G * (p - 3)), G)
    base$l_em <- matrix(stats::rnorm(G * (p - 3)), G)
  }
  base
}

.reg_coef_std_draws_at <- function(base, n, beta, R2, R2_j) {
  p <- base$p
  s_u <- sqrt(1 - R2_j)
  s_e <- sqrt(1 - R2)
  draws <- base[setdiff(names(base), c("u_u", "u_e", "u_11", "u_22"))]
  draws$l_u <- sqrt(stats::qchisq(base$u_u, n - p))
  draws$l_e <- sqrt(stats::qchisq(base$u_e, n - p - 1))
  if (p >= 2) draws$l_11 <- sqrt(stats::qchisq(base$u_11, n - 1))
  if (p >= 3) draws$l_22 <- sqrt(stats::qchisq(base$u_22, n - 2))
  b <- beta + s_e * draws$l_eu / (s_u * draws$l_u)
  draws$se_raw <- s_e * draws$l_e / (s_u * draws$l_u * sqrt(n - p - 1))
  draws$t <- b / draws$se_raw
  ss_x <- s_u^2 * draws$l_u^2
  if (p >= 2) ss_x <- ss_x + (sqrt(R2_j) * draws$l_11 + s_u * draws$l_u1)^2
  if (p >= 3) ss_x <- ss_x + s_u^2 * draws$l_u2^2
  if (p >= 4) ss_x <- ss_x + s_u^2 * rowSums(draws$l_um^2)
  draws$ss_x <- ss_x
  # The parts of SS_Y that do not involve d.
  ss_y <- (beta * s_u * draws$l_u + s_e * draws$l_eu)^2 + s_e^2 * draws$l_e^2
  if (p >= 4) ss_y <- ss_y + rowSums((beta * s_u * draws$l_um +
                                        s_e * draws$l_em)^2)
  draws$ss_y_fixed <- ss_y
  draws$n <- n
  draws$beta <- beta
  draws$R2 <- R2
  draws$R2_j <- R2_j
  draws
}

.reg_coef_std_draws <- function(G, n, p, beta, R2, R2_j) {
  .reg_coef_std_draws_at(.reg_coef_std_base_draws(G, p), n, beta, R2, R2_j)
}

# .reg_coef_std_d(R2, R2_j, beta, ad, p)
#
# The coordinates (d_1, d_2) of d for a'd = ad: d_1 = ad / |a| and d_2 the
# rest of |d|; with |a| = 0 all of d is put on the first axis, which then
# does not enter SS_X. With p = 2 there is one axis and d_2 is 0.
.reg_coef_std_d <- function(R2, R2_j, beta, ad, p) {
  a <- sqrt(R2_j)
  d_sq <- max(0, R2 - beta^2 * (1 - R2_j))
  if (a > 0) {
    d_1 <- ad / a
    d_2 <- if (p >= 3) sqrt(max(0, d_sq - d_1^2)) else 0
  } else {
    d_1 <- sqrt(d_sq)
    d_2 <- 0
  }
  c(d_1, d_2)
}

# .reg_coef_std_width(draws, ad, delta)
#
# The widths of the noncentral standardized interval for the draws, at
# a'd = ad, with delta() the distance between the noncentrality limits as a
# function of t (.reg_coef_nc_delta()).
.reg_coef_std_width <- function(draws, ad, delta) {
  p <- draws$p
  s_u <- sqrt(1 - draws$R2_j)
  s_e <- sqrt(1 - draws$R2)
  ss_y <- draws$ss_y_fixed
  if (p >= 2) {
    d <- .reg_coef_std_d(draws$R2, draws$R2_j, draws$beta, ad, p)
    col_1 <- d[1] * draws$l_11 + draws$beta * s_u * draws$l_u1 +
      s_e * draws$l_e1
    if (p >= 3) {
      col_1 <- col_1 + d[2] * draws$l_21
      ss_y <- ss_y + (d[2] * draws$l_22 + draws$beta * s_u * draws$l_u2 +
                        s_e * draws$l_e2)^2
    }
    ss_y <- ss_y + col_1^2
  }
  delta(draws$t) * draws$se_raw * sqrt(draws$ss_x / ss_y)
}

# .reg_coef_nc_delta(t, df, alpha_lower, alpha_upper, nodes = 41)
#
# The distance between the noncentrality limits that ci_nc_t() gives at an
# observed t on df degrees of freedom, as a function of t over the range of
# the values in `t`. ci_nc_t() is evaluated at `nodes` Chebyshev points of
# that range and the ratio of the distance to sqrt(1 + t^2 / (2 df)), which
# is nearly constant (the distance grows like that square root), is
# interpolated by a cubic spline. Against ci_nc_t() at the draws themselves
# the relative error was below 1e-9 in the checks recorded in the tests.
# Every limit therefore comes from ci_nc_t(), the function the interval of
# the analysis uses.
.reg_coef_nc_delta <- function(t, df, alpha_lower, alpha_upper, nodes = 41L) {
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
                   verbose = FALSE)$value
    lim[2] - lim[1]
  }, numeric(1))
  scale <- function(v) sqrt(1 + v^2 / (2 * df))
  f <- stats::splinefun(x, dist / scale(x), method = "fmm")
  function(v) f(v) * scale(v)
}

# .reg_coef_central_delta(df, alpha_lower, alpha_upper)
#
# The counterpart of .reg_coef_nc_delta() for the central interval: the
# width of the central t interval in units of its standard error,
# qt(1 - alpha_lower) + qt(1 - alpha_upper) on df degrees of freedom
# (each computed in the upper tail), whatever t is.
.reg_coef_central_delta <- function(df, alpha_lower, alpha_upper) {
  q <- stats::qt(alpha_lower, df, lower.tail = FALSE) +
    stats::qt(alpha_upper, df, lower.tail = FALSE)
  function(v) rep(q, length(v))
}

# .reg_coef_std_share(draws, ad_values, width, alpha_lower, alpha_upper,
#                     interval = "noncentral")
#
# The proportion of the draws whose interval for the standardized
# coefficient is no wider than `width`, for each value of a'd in ad_values
# (one set of draws serves them all: t and the residual sums of squares do
# not depend on a'd). `interval` is "noncentral" (ci_reg_coef(noncentral =
# TRUE) on standardized data) or "central" (ci_reg_coef(noncentral =
# FALSE) on standardized data, Kelley and Maxwell, 2003, Equation 1).
.reg_coef_std_share <- function(draws, ad_values, width, alpha_lower,
                                alpha_upper,
                                interval = c("noncentral", "central")) {
  interval <- match.arg(interval)
  df <- draws$n - draws$p - 1
  delta <- if (interval == "noncentral")
    .reg_coef_nc_delta(draws$t, df, alpha_lower, alpha_upper) else
      .reg_coef_central_delta(df, alpha_lower, alpha_upper)
  vapply(ad_values, function(ad)
    mean(.reg_coef_std_width(draws, ad, delta) <= width), numeric(1))
}

# .reg_coef_std_ad_values(R2, R2_j, beta, p, cor_YXj)
#
# The values of a'd = rho_YXj - beta (1 - R2_j) the step evaluates. With
# the correlation of Y with X_j known (the covariance structure path) it is
# that one value. Otherwise (the squared multiple correlations path) it is
# `grid` (21) equally spaced values across the range that is possible with
# R2, R2_j, and beta: a'd lies in [-|a||d|, |a||d|]. The smallest share
# over them is an approximate minimum over the range, since the share
# between two grid values can be a little lower than at either. With
# p = 2, where a and d are numbers, a'd is one of the two ends, and both
# are evaluated. With p = 1 there is no Z and a'd does not enter.
.reg_coef_std_ad_values <- function(R2, R2_j, beta, p, cor_YXj = NULL,
                                    grid = 21L) {
  if (p == 1) return(0)
  if (!is.null(cor_YXj)) return(cor_YXj - beta * (1 - R2_j))
  bound <- sqrt(R2_j) * sqrt(max(0, R2 - beta^2 * (1 - R2_j)))
  if (bound == 0) return(0)
  if (p == 2) return(c(-bound, bound))
  seq(-bound, bound, length.out = grid)
}

# .reg_coef_std_assurance_n(n_0, p, beta, R2, R2_j, ad_values, width,
#                           alpha_lower, alpha_upper, assurance, G,
#                           interval = "noncentral")
#
# The simulated assurance step for an interval for the standardized
# coefficient, the noncentral one or the central one (`interval`, as in
# .reg_coef_std_share()):
# the smallest sample size, at or above n_0 (the plan without assurance),
# at which the share of G simulated intervals no wider than `width` is at
# least `assurance`, the share being the smallest over ad_values (one
# value on the covariance structure path; on the other path the smallest
# over the grid of .reg_coef_std_ad_values(), an approximate minimum over
# the populations the inputs allow). The same G samples serve every n
# tried (.reg_coef_std_base_draws()), so the share changes smoothly with
# n. The search of .reg_coef_assurance_search() (doubling, then
# bisection) assumes that the share does not decrease with n, which held
# in every set checked though a single simulated interval can be wider at
# n + 1 (see .reg_coef_std_base_draws() above), and where it holds the
# search returns the smallest n at which the share reaches the request.
# Returns list(n, attained, mc_se): attained is the share at n,
# the estimate the plan rests on, so it is at least `assurance`; mc_se is
# its Monte Carlo standard error, sqrt{attained (1 - attained) / G}. The
# exact probability at n differs from the share by Monte Carlo error, so
# n can differ by a case or two from the n the exact probability gives.
# The caller sets the seed.
.reg_coef_std_assurance_n <- function(n_0, p, beta, R2, R2_j, ad_values,
                                      width, alpha_lower, alpha_upper,
                                      assurance, G,
                                      interval = c("noncentral", "central")) {
  interval <- match.arg(interval)
  base <- .reg_coef_std_base_draws(G, p)
  share <- function(n)
    min(.reg_coef_std_share(.reg_coef_std_draws_at(base, n, beta, R2, R2_j),
                            ad_values, width, alpha_lower, alpha_upper,
                            interval))
  step <- .reg_coef_assurance_search(n_0, share, assurance)
  list(n = step$n, attained = step$attained,
       mc_se = sqrt(step$attained * (1 - step$attained) / G))
}

# .reg_coef_std_population(true_R2_Y_X, true_R2_j_X_without_j, p,
#                          true_std_reg_coef_j, true_cor_XX, true_cor_YX,
#                          which_predictor)
#
# What the simulated step needs to know of the population, from the
# arguments of ss_aipe_src() after ss_aipe_reg_coef() has checked them:
# list(p, R2, R2_j, beta, cor_YXj). On the covariance structure path
# (true_cor_XX and true_cor_YX) they are derived as ss_aipe_reg_coef()
# derives them, and cor_YXj is the correlation of Y with the targeted
# predictor; on the squared multiple correlations path they are the values
# given, and cor_YXj is NULL (.reg_coef_std_ad_values() then takes every
# value the three allow).
.reg_coef_std_population <- function(true_R2_Y_X, true_R2_j_X_without_j, p,
                                     true_std_reg_coef_j, true_cor_XX,
                                     true_cor_YX, which_predictor) {
  if (!is.null(true_cor_XX) && !is.null(true_cor_YX)) {
    R <- as.matrix(true_cor_XX)
    r <- as.numeric(true_cor_YX)
    p <- ncol(R)
    j <- if (is.null(which_predictor)) 1 else which_predictor
    coef <- solve(R, r)
    list(p = p, R2 = sum(r * coef),
         R2_j = if (p > 1) 1 - 1 / solve(R)[j, j] else 0,
         beta = coef[j], cor_YXj = r[j])
  } else {
    list(p = p, R2 = as.numeric(true_R2_Y_X),
         R2_j = as.numeric(true_R2_j_X_without_j),
         beta = as.numeric(true_std_reg_coef_j), cor_YXj = NULL)
  }
}
