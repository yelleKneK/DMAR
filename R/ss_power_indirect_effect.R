#' Sample Size Planning for Power for the Indirect (Mediation) Effect
#'
#' Power, or the sample size required for a desired power, for the test of
#' the indirect effect \eqn{a b} in the simple mediation model, with paths
#' specified in standardized metric (unit-variance \eqn{X}, \eqn{M}, and
#' \eqn{Y}). The default test is joint significance (the indirect effect is
#' declared when both \eqn{\hat a} and \eqn{\hat b} are individually
#' significant), which tracks the resampling tests' power closely and far
#' exceeds the Sobel test in small samples (Fritz & MacKinnon, 2007); the
#' Sobel normal-theory test is available for comparison. Power comes from
#' the exact sampling distribution of the two path \emph{t} statistics,
#' with \eqn{X} random by default. This is the power
#' counterpart of the accuracy in parameter estimation (AIPE) planner
#' \code{\link{ss_aipe_indirect_effect}}, and the planning complement of
#' the analysis function \code{\link{mediate}}.
#'
#' @param true_std_reg_coef_a Population standardized regression
#'   coefficient of the \eqn{X \to M} path \eqn{a}.
#' @param true_std_reg_coef_b Population standardized regression
#'   coefficient of the \eqn{M \to Y} path \eqn{b}, holding \eqn{X}.
#' @param true_std_reg_coef_direct Population standardized regression
#'   coefficient of the direct effect \eqn{c'} of \eqn{X} on \eqn{Y}
#'   holding \eqn{M}. Defaults to 0; it enters only through the residual
#'   variance of \eqn{Y}.
#' @param desired_power Desired power; supply this to solve for \eqn{N}.
#' @param specified_N Total sample size \eqn{N}; supply this to evaluate the
#'   realized power. Specify exactly one of \code{desired_power} and
#'   \code{specified_N}.
#' @param alpha_level Two-sided Type I error rate for each component test.
#'   Defaults to 0.05.
#' @param method \code{"joint_significance"} (default) or \code{"sobel"}.
#' @param random_predictors \code{TRUE} (default) when \eqn{X} is random,
#'   as when it is measured, and \code{FALSE} when it is fixed by design, as
#'   when groups are assigned. The mediator is random under either setting.
#'
#' @details
#' With unit-variance variables, \eqn{M = a X + e_M} and
#' \eqn{Y = c' X + b M + e_Y}, with normal errors of variances \eqn{1 - a^2}
#' and \eqn{\sigma^2_{e_Y} = 1 - (b^2 + c'^2 + 2 a b c')}, which must be
#' positive, or the function stops. As in \code{\link{mediate}}, \eqn{t_a}
#' tests the slope of \eqn{M} on \eqn{X} (\eqn{N - 2} degrees of freedom)
#' and \eqn{t_b} the coefficient of \eqn{M} in the regression of \eqn{Y} on
#' \eqn{X} and \eqn{M} (\eqn{N - 3}). Given \eqn{X}, \eqn{t_a} is noncentral
#' \emph{t} with noncentrality \eqn{a \sqrt{SS_X / (1 - a^2)}}; given
#' \eqn{X} and \eqn{M}, \eqn{t_b} is noncentral \emph{t} with noncentrality
#' \eqn{b \sqrt{SS_{M \cdot X}} / \sigma_{e_Y}}, where \eqn{SS_{M \cdot X}}
#' is the residual sum of squares of \eqn{M} on \eqn{X}.
#' \eqn{SS_{M \cdot X} / (1 - a^2)} is chi square on \eqn{N - 2} degrees of
#' freedom in every design, because \eqn{M} is random given \eqn{X};
#' \eqn{SS_X} is chi square on \eqn{N - 1} with \eqn{X} random, and \eqn{N}
#' with \eqn{X} fixed at variance 1 (divisor \eqn{N}; two groups of equal
#' size, for example). Each power averages the noncentral \emph{t} power
#' over these chi squares, by Gauss-Hermite sums on their normal scores.
#'
#' \eqn{SS_{M \cdot X}} is the error sum of squares of one test and the
#' predictor sum of squares of the other, so the tests are dependent:
#' \code{power_ab} is the probability that both reject, which is below
#' \code{power_a} times \code{power_b}. The first-order Sobel statistic of
#' \code{\link{mediate}} equals \eqn{t_a t_b / \sqrt{t_a^2 + t_b^2}}, and
#' its power is the probability, integrated numerically, that it exceeds
#' the normal critical value; \code{power_a} and \code{power_b} are then
#' \code{NA}. Because it refers \emph{t} statistics to the normal, its
#' power need not rise with \eqn{N}: with paths of .9 and .05 it falls from
#' .0729 at \eqn{N = 10} to .0562 at 59, and with paths of .8 and .02 and
#' \code{alpha_level = .5} it is .5072 at 12 and .5035 at 59, and does not
#' reach .507 again until 216. The joint significance power rose with
#' \eqn{N} in each of 1,000 settings checked at every \eqn{N} from 10 to
#' 1,500 (paths from .001 to .99 in absolute value, \code{alpha_level} from
#' .0001 to .997, \eqn{X} random or fixed), which is not a proof. So the
#' search assumes neither: \code{necessary_N} is the smallest \eqn{N} of at
#' least 10 whose power reaches \code{desired_power}. It scans upward,
#' passing over a range of \eqn{N} only where a bound shows that none
#' reaches it: written as its quantile at a fixed normal score, each chi
#' square rises with \eqn{N}, so the power at every \eqn{N} of the range is
#' at most the power with the noncentralities at the top of the range and
#' each test's threshold (critical value times scale) at its least over the
#' range. Where a path is near zero, the power is also at most the rate at
#' which the test of that path rejects when the path is zero plus
#' \eqn{\phi(1)} (the standard normal density at 1) times the expected
#' square of its noncentrality (for the Sobel method, the test of
#' \eqn{|t|} against the normal critical value, since the Sobel statistic
#' is at most the smaller of \eqn{|t_a|} and \eqn{|t_b|}), which lets the
#' search pass at once over a long range in which the power stays just
#' below \code{desired_power}. It stops at \eqn{N} = 10,000,000 with an
#' error of class \code{dmar_computation_failed}. The tests check
#' \code{power_ab} against raw-data simulation within Monte Carlo error at
#' the plan for power .80 and at one case fewer.
#'
#' Fixing each noncentrality at its population standard error and
#' multiplying the two powers overstates the joint significance power: with
#' paths of .39 that approximation gives .8017 at \eqn{N = 65}, where the
#' power is .7675, and the plan for .80 is 69.
#'
#' @return A tidy \code{data.frame} with
#'   \code{necessary_N} (or \code{specified_N}), \code{power_ab} (the
#'   power to detect the indirect effect \eqn{a b}, the quantity the sample
#'   size is planned against), the powers of the tests of the two paths
#'   (\code{power_a}, \code{power_b}; \code{NA} for the Sobel method), the
#'   paths (\code{true_std_reg_coef_a}, \code{true_std_reg_coef_b},
#'   \code{true_std_reg_coef_direct}), their implied standardized indirect
#'   effect \code{true_std_indirect_effect}, and \code{alpha_level}. The
#'   method is recorded in the \code{"method"} attribute and the predictor
#'   model in the \code{"random_predictors"} attribute. The result carries
#'   the \code{dmar_ss_power} class, so
#'   \code{\link[generics]{tidy}} reports the sample size and the power to
#'   detect the indirect effect (the row \code{power_ab}), and
#'   \code{\link[generics]{glance}} adds the component powers and the
#'   planning inputs.
#'
#' @references
#' Fritz, M. S., & MacKinnon, D. P. (2007). Required sample size to detect
#'   the mediated effect. \emph{Psychological Science, 18}(3), 233--239.
#'   \doi{10.1111/j.1467-9280.2007.01882.x}
#'
#' MacKinnon, D. P., Lockwood, C. M., Hoffman, J. M., West, S. G., &
#'   Sheets, V. (2002). A comparison of methods to test mediation and
#'   other intervening variable effects. \emph{Psychological Methods,
#'   7}(1), 83--104. \doi{10.1037/1082-989X.7.1.83}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{mediate}} to analyze the study this plans;
#'   \code{\link{ss_aipe_indirect_effect}} to plan for confidence interval
#'   width instead of detection; \code{\link{design_consequences}} for what
#'   the chosen design delivers.
#'
#' @family mediation
#' @family sample size for power
#'
#' @keywords design
#'
#' @examples
#' # Fritz and MacKinnon's (2007) running scenario (a = b = .39), X random.
#' ss_power_indirect_effect(true_std_reg_coef_a = .39,
#'                          true_std_reg_coef_b = .39, desired_power = .80)
#'
#' # X fixed by design (two groups of equal size); the mediator is random.
#' ss_power_indirect_effect(true_std_reg_coef_a = .39,
#'                          true_std_reg_coef_b = .39, desired_power = .80,
#'                          random_predictors = FALSE)
#'
#' # A near-zero a path against a larger b: the weak link drives the requirement.
#' ss_power_indirect_effect(true_std_reg_coef_a = .14,
#'                          true_std_reg_coef_b = .39, desired_power = .80)
#'
#' # Realized power at a given N, and the Sobel test at the same N.
#' ss_power_indirect_effect(true_std_reg_coef_a = .39,
#'                          true_std_reg_coef_b = .39, specified_N = 75)
#' ss_power_indirect_effect(true_std_reg_coef_a = .39,
#'                          true_std_reg_coef_b = .39, specified_N = 75,
#'                          method = "sobel")
#'
#' @export
#' @importFrom stats pnorm qnorm
ss_power_indirect_effect <- function(true_std_reg_coef_a, true_std_reg_coef_b,
                                     true_std_reg_coef_direct = 0,
                                     desired_power = NULL, specified_N = NULL,
                                     alpha_level = 0.05,
                                     method = c("joint_significance",
                                                "sobel"),
                                     random_predictors = TRUE) {
  method <- match.arg(method)
  for (nm in c("true_std_reg_coef_a", "true_std_reg_coef_b",
               "true_std_reg_coef_direct")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        abs(val) >= 1) {
      .dmar_refuse_argument(
        sprintf("'%s' must be a single standardized path in (-1, 1).",
                nm), arg = nm)
    }
  }
  if (true_std_reg_coef_a == 0 || true_std_reg_coef_b == 0) {
    .dmar_refuse_argument(c(
      "'true_std_reg_coef_a' and 'true_std_reg_coef_b' must be nonzero: a ",
      "zero path has no indirect effect to detect."),
      arg = c("true_std_reg_coef_a", "true_std_reg_coef_b"))
  }
  if (is.null(desired_power) == is.null(specified_N)) {
    .dmar_refuse_argument(
      "Specify exactly one of 'desired_power' or 'specified_N'.",
      arg = c("desired_power", "specified_N"))
  }
  if (!is.null(desired_power)) .dmar_check_desired_power(desired_power)
  if (!is.null(specified_N) && (!is.numeric(specified_N) || length(specified_N) != 1L || is.na(specified_N) ||
                      specified_N < 10 || specified_N != round(specified_N))) {
    .dmar_refuse_argument(
      "'specified_N' must be a single integer of at least 10.",
      arg = "specified_N")
  }
  if (!is.numeric(alpha_level) || length(alpha_level) != 1L ||
      is.na(alpha_level) || alpha_level <= 0 || alpha_level >= 1) {
    .dmar_refuse_argument(
      "'alpha_level' must be a single number in (0, 1).", arg = "alpha_level")
  }
  if (!isTRUE(random_predictors) && !isFALSE(random_predictors)) {
    .dmar_refuse_argument(
      "'random_predictors' must be TRUE or FALSE.", arg = "random_predictors")
  }

  resid_M <- 1 - true_std_reg_coef_a^2
  resid_Y <- 1 - (true_std_reg_coef_b^2 + true_std_reg_coef_direct^2 +
                    2 * true_std_reg_coef_a * true_std_reg_coef_b *
                    true_std_reg_coef_direct)
  if (resid_Y <= 0) {
    .dmar_refuse_argument(c(
      "The supplied standardized paths imply a non-positive residual ",
      "variance for the outcome; the combination is not admissible."),
      arg = NULL)
  }
  delta_a <- abs(true_std_reg_coef_a) / sqrt(resid_M)
  delta_b <- abs(true_std_reg_coef_b) * sqrt(resid_M) / sqrt(resid_Y)
  power_at <- function(n, n_lo = n) {
    .indirect_power(n, delta_a, delta_b, alpha_level, method,
                    random_predictors, n_lo)
  }

  if (is.null(specified_N)) {
    # The smallest N from 10 whose power reaches desired_power, which need
    # not rise with N (Details). Every N up to n_short falls short. First,
    # .indirect_null_reach() passes over every N up to the last at which its
    # bound for a near-zero path stays below the target. Then the block
    # n_short + 1 to n_short + h is passed over when the bound of
    # .indirect_power() falls short by more than its numerical error, and h
    # doubles; otherwise h is halved, down to one N, whose power is exact.
    n_max <- 10000000L
    n_short <- 9L
    h <- 1L
    repeat {
      reach <- .indirect_null_reach(n_short + 1L, desired_power - 1e-8,
                                    delta_a, delta_b, alpha_level, method)
      if (reach > n_short) n_short <- as.integer(min(reach, n_max))
      if (n_short >= n_max)
        .dmar_refuse(c(
          "Could not reach 'desired_power' at a sample size up to ",
          format(n_max, scientific = FALSE), "; the indirect effect may be ",
          "too small to detect at this power."), "computation_failed")
      n_i <- n_short + 1L
      n_hi <- min(n_short + h, n_max)
      pw <- power_at(n_hi, n_i)
      if (n_hi == n_i && pw[1] >= desired_power) break
      if (n_hi > n_i && pw[1] >= desired_power - 1e-8) {
        h <- (n_hi - n_short) %/% 2L
        next
      }
      n_short <- n_hi
      h <- 2L * h
    }
    label <- "necessary_N"
  } else {
    n_i <- as.integer(specified_N)
    label <- "specified_N"
    pw <- power_at(n_i)
  }

  out <- data.frame(
    term  = c(label, "power_ab", "power_a", "power_b",
              "true_std_reg_coef_a", "true_std_reg_coef_b",
              "true_std_reg_coef_direct", "true_std_indirect_effect",
              "alpha_level"),
    value = c(n_i, pw[1], pw[2], pw[3], true_std_reg_coef_a,
              true_std_reg_coef_b, true_std_reg_coef_direct,
              true_std_reg_coef_a * true_std_reg_coef_b, alpha_level),
    stringsAsFactors = FALSE
  )
  # tidy() reports power_ab with the sample size through the power registry
  # in R/dmar_tidiers.R; glance() adds power_a and power_b.
  class(out) <- c("dmar_ss_power", class(out))
  out <- .as_dmar_tbl(out)
  attr(out, "method") <- method
  attr(out, "random_predictors") <- random_predictors
  out
}

# The last N (Inf for every N; -Inf for none) up to which no N from n_lo has
# power_ab above `target`, by the bound for a near-zero path (Details). Let
# c be the critical value of a path test: qt(1 - alpha_level / 2, df) for
# joint significance and z for Sobel, whose statistic
# |t_a t_b| / sqrt(t_a^2 + t_b^2) is at most the smaller of |t_a| and |t_b|,
# so power_ab is at most P(|t_a| > c) and at most P(|t_b| > c). With
# g(m) = P(|Z + m| > T), g'(0) = 0 and |g''| <= 2 dnorm(1) (the largest
# |x| dnorm(x)), so g(m) <= g(0) + dnorm(1) m^2; averaged over the chi
# squares, P(|t_a| > c) is at most P(|t| > c), t central on N - 2 df, plus
# dnorm(1) delta_a^2 E(SS_X) (N - 1, or N with X fixed), and P(|t_b| > c)
# at most the same on N - 3 df plus dnorm(1) delta_b^2 (N - 2). The first
# term is alpha_level for joint significance and 2 pt(-z, df) for Sobel,
# which falls as df rises (it is E 2 pnorm(-z sqrt(Q)), Q chi square over
# its df, a convex function of Q, and Q falls in the convex order as df
# rises, a mean of more squares), so from n_lo to N it is at most its value
# at n_lo.
.indirect_null_reach <- function(n_lo, target, delta_a, delta_b, alpha_level,
                                 method) {
  rate <- function(df) if (method == "sobel")
    2 * stats::pt(stats::qnorm(alpha_level / 2), df) else alpha_level
  last <- function(room, slope, shift) {
    if (room <= 0) return(-Inf)
    floor(room / slope) + shift
  }
  slope <- stats::dnorm(1)
  max(last(target - rate(n_lo - 2), slope * delta_a^2, 0),
      last(target - rate(n_lo - 3), slope * delta_b^2, 2))
}

# c(power_ab, power_a, power_b) at n (Details). With V = SS_X (n with X
# fixed), W = SS_M.X / (1 - a^2), U = SSE_Y / var(e_Y), and Z_a, Z_b normal,
# all independent, t_a = (Z_a + delta_a sqrt(V)) / sqrt(W / (n - 2)) and
# t_b = (Z_b + delta_b sqrt(W)) / sqrt(U / (n - 3)), so the tests are
# independent given W. With X fixed, power_a is .nc_t_power()'s.
# With n_lo < n, c(bound, NA, NA), a bound of power_ab at every N from n_lo
# to n. Each chi square is its quantile at the normal score of a node, so at
# a node it rises with N, as the means delta sqrt(V) and delta sqrt(W) do,
# and the critical values fall. Over a run of N, the threshold of |Z + mean|
# (critical value times sqrt(chi square / df); for Sobel, the scale of t) is
# at least its value with the chi square at the run's first N and the rest
# at its last. The bound takes the means at n and each threshold at its
# least over the runs, so at each node the probability that either test, or
# the Sobel test (an upper set in |t_a| and |t_b|), rejects is at least that
# at any N from n_lo to n. More runs make the bound closer to the power, at
# 48 chi square quantiles a run: 256 for Sobel, whose integral costs far
# more, and 16 for joint significance.
.indirect_power <- function(n, delta_a, delta_b, alpha_level, method,
                            random_predictors, n_lo = n) {
  W <- .indirect_chi_nodes(n - 2)
  V <- if (random_predictors) .indirect_chi_nodes(n - 1) else
    list(s = sqrt(n), w = 1)
  sobel <- method == "sobel"
  crit <- function(df) if (sobel) 1 else stats::qt(1 - alpha_level / 2, df)
  scale <- function(d) {
    one <- function(lo, hi)
      crit(hi - d) * .indirect_chi_nodes(lo - d)$s / sqrt(hi - d)
    if (n_lo == n) return(one(n, n))
    e <- unique(round(seq(n_lo - 1, n,
                          length.out = if (sobel) 257L else 17L)))
    Reduce(pmin, Map(one, e[-length(e)] + 1, e[-1]))
  }
  crit_a <- scale(2)
  crit_b <- scale(3)
  if (sobel) {
    z <- stats::qnorm(1 - alpha_level / 2)
    return(c(.indirect_sobel_power(n, delta_a, delta_b, z, V, W, crit_a,
                                   crit_b, n_lo), NA_real_, NA_real_))
  }
  mean_a <- delta_a * V$s
  reject_a <- drop((stats::pnorm(outer(-crit_a, mean_a, "+")) +
                      stats::pnorm(outer(-crit_a, -mean_a, "+"))) %*% V$w)
  mean_b <- delta_b * W$s
  reject_b <- drop((stats::pnorm(outer(mean_b, -crit_b, "+")) +
                      stats::pnorm(outer(-mean_b, -crit_b, "+"))) %*% W$w)
  if (n_lo < n) return(c(sum(W$w * reject_a * reject_b), NA_real_, NA_real_))
  power_a <- if (random_predictors) sum(W$w * reject_a) else
    .nc_t_power(crit(n - 2), n - 2, delta_a * sqrt(n))
  c(sum(W$w * reject_a * reject_b), power_a, sum(W$w * reject_b))
}

# The n-point Gauss-Hermite rule for the standard normal, built as
# .gauss_legendre() (R/reliability_internals.R) builds its rule.
.gauss_hermite <- function(n) {
  i <- seq_len(n - 1L)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(i, i + 1L)] <- sqrt(i)
  jacobi[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(jacobi, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1L, ]^2)
}

.INDIRECT_GAUSS_HERMITE <- .gauss_hermite(24L)

# sqrt(V) at the nodes, V chi square on df (upper tail for x > 0).
.indirect_chi_nodes <- function(df, rule = .INDIRECT_GAUSS_HERMITE) {
  x <- rule$x
  v <- numeric(length(x))
  up <- x > 0
  v[up] <- stats::qchisq(stats::pnorm(-x[up]), df, lower.tail = FALSE)
  v[!up] <- stats::qchisq(stats::pnorm(x[!up]), df)
  list(s = sqrt(v), w = rule$w)
}

# The Sobel test rejects when |t_a| > z and |t_b| > k(|t_a|), k(y) =
# y z / sqrt(y^2 - z^2), and never when |t_a| <= z, where k is infinite.
# y_lo = k(k_top) is above z, but with k_top more than about 6e7 times z
# it is within rounding of z, and y^2 - z^2 can come out negative near
# y_lo; k_of() holds y^2 - z^2 at 0 or above, so k is infinite there and
# the integrand is 0 rather than NaN. Given W, t_a and t_b are independent
# normal mixtures over V and U, and the power sums over W the integral over
# y = |t_a| of f(y | W) P(|t_b| > k(y) | W): negligible below y_lo and past
# max(mu1) + 10 s1, P(|t_a| > y_one | W) past y_one, where the second factor
# is 1 within 1e-18, and between by integrate() on pieces split toward z
# and, where f(y | W) has separate peaks, around each of them.
# sigma_a[i] and tau are the scales of t_a at W node i and of t_b at the U
# nodes, sqrt(W / (n - 2)) and sqrt(U / (n - 3)), or the bound's for N
# from n_lo to n.
.indirect_sobel_power <- function(n, delta_a, delta_b, z, V, W, sigma_a, tau,
                                  n_lo = n) {
  what <- if (n_lo < n) {
    paste0("A bound of the power of the Sobel test for N from ",
           format(n_lo, scientific = FALSE), " to ",
           format(n, scientific = FALSE))
  } else {
    paste0("The power of the Sobel test at N = ", format(n, scientific = FALSE))
  }
  k_of <- function(y) y * z / sqrt(pmax(y^2 - z^2, 0))
  s2 <- 1 / tau
  total <- 0
  for (i in seq_along(W$s)) {
    sigma <- sigma_a[i]
    s1 <- 1 / sigma
    mu1 <- delta_a * V$s / sigma
    mu2 <- delta_b * W$s[i] / tau
    k_top <- max(mu2 + 9 * s2)
    if (k_top <= z) next
    k_one <- min(mu2 - 9 * s2)
    y_lo <- k_of(k_top)
    y_one <- if (k_one > z) k_of(k_one) else Inf
    y_top <- max(mu1) + 10 * s1
    y_end <- min(y_one, y_top)
    part <- 0
    if (y_one < y_top) {
      part <- sum(V$w * (stats::pnorm((mu1 - y_one) / s1) +
                           stats::pnorm((-mu1 - y_one) / s1)))
    }
    if (y_lo < y_end) {
      integrand <- function(y) {
        k <- k_of(y)
        m <- length(y)
        dens <- drop((stats::dnorm(outer(y, mu1, "-") / s1) +
                        stats::dnorm(outer(y, -mu1, "-") / s1)) %*% V$w) / s1
        mean2 <- rep(mu2, each = m)
        sd2 <- rep(s2, each = m)
        reject_b <- stats::pnorm((mean2 - k) / sd2) +
          stats::pnorm((-mean2 - k) / sd2)
        dim(reject_b) <- c(m, length(mu2))
        dens * drop(reject_b %*% W$w)
      }
      breaks <- c(y_lo, y_end, z + (y_lo - z) * 4^(1:40),
                  k_of(seq(max(k_one, z * (1 + 1e-6)), k_top,
                           length.out = 6)),
                  min(mu1) + s1 * c(-10, -5, -2, 0),
                  max(mu1) + s1 * c(0, 2, 5), 10 * s1)
      # With X random, f(y | W) is a mixture over the V nodes of normal
      # densities with sd s1 at mu1. Where two neighboring means stand
      # more than 2 s1 apart, which takes |a| above about .9, the mixture
      # has separate peaks, up to about 200 s1 apart at |a| = .99999, and
      # integrate() on a piece much longer than s1 can step over one: with
      # a = .999999, b = .5, N = 10, and alpha_level = .5, the power was
      # .352 where it is .522. Each mean then gets a break at it and 10 s1
      # either side, where its density is 2e-22 of its height; at 6 s1 a
      # piece that ended on the flank of a peak came back "the integral
      # is probably divergent".
      if (length(mu1) > 1L && max(diff(sort(mu1))) > 2 * s1) {
        breaks <- c(breaks, outer(mu1, s1 * c(-10, 0, 10), "+"))
      }
      breaks <- sort(unique(breaks[breaks >= y_lo & breaks <= y_end]))
      # The budget of absolute error for a piece that integrate() returns
      # with a roundoff status (.indirect_sobel_piece()): 1e-9 shared
      # equally by the pieces of this node. A piece returned "OK" has an
      # estimated error of at most the larger of 1e-15 and 1e-10 of its
      # value, and the pieces of a node, at most 127, sum to at most 1, so
      # the node's integral carries an estimated error below 1.2e-9; the W
      # weights sum to 1, so the power does too. That is well under the
      # 1e-8 by which the search requires a bound to fall short before it
      # passes over N.
      budget <- 1e-9 / (length(breaks) - 1L)
      for (j in seq_len(length(breaks) - 1L)) {
        part <- part + .indirect_sobel_piece(integrand, breaks[j],
                                             breaks[j + 1L], budget, what)
      }
    }
    total <- total + W$w[i] * part
  }
  min(max(total, 0), 1)
}

# One piece of the integral of the Sobel power, from lower to upper.
# integrate() runs with stop.on.error = FALSE, so that its status can be
# read (the strings are those of R's integrate()). A piece returned "OK"
# is accepted whatever its error estimate. A piece returned with either
# roundoff status, "roundoff error was detected" or "roundoff error is
# detected in the extrapolation table", is accepted when its estimated
# absolute error is finite and at most `budget` and its value, the
# integral of a nonnegative function, is at least -budget. QUADPACK gives
# these statuses when subdividing stops reducing its error estimate, which
# happens on short pieces next to z when alpha_level is near 1, pieces
# worth little: at alpha_level = .999 and N = 10, less than 1e-45, with
# estimated errors near 1e-47. A roundoff piece over its budget is refused
# with the status, the error, and the budget, and the message names the
# joint significance test, which does not need the integral. Every other
# status, and an error of integrate() itself, is refused with the status
# in the message.
.indirect_sobel_piece <- function(integrand, lower, upper, budget, what) {
  refuse <- function(detail, remedy = NULL) .dmar_refuse(c(
    what, " could not ",
    "be computed reliably: the numerical integral over the ",
    "distribution of the path t statistics failed (", detail, ").",
    remedy), "computation_failed")
  q <- tryCatch(stats::integrate(
    integrand, lower, upper, rel.tol = 1e-10, abs.tol = 1e-15,
    subdivisions = 1000L, stop.on.error = FALSE),
    error = function(e) refuse(conditionMessage(e)))
  if (identical(q$message, "OK")) return(q$value)
  roundoff <- c("roundoff error was detected",
                "roundoff error is detected in the extrapolation table")
  if (!q$message %in% roundoff) refuse(q$message)
  below <- is.finite(q$value) && q$value < -budget
  if (is.finite(q$value) && is.finite(q$abs.error) &&
      q$abs.error <= budget && !below) {
    return(q$value)
  }
  refuse(c(q$message, ", with an estimated absolute error of ",
           format(q$abs.error, digits = 3), " against a budget of ",
           format(budget, digits = 3), " for this piece",
           if (below) c(" and a value of ", format(q$value, digits = 3),
                        " for the integral of a nonnegative function")),
         c(" The joint significance test (method = ",
           "\"joint_significance\") does not need this integral."))
}
