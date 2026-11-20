# The limits of .ratio_boot_limits() from a vector of replicates, through
# the two steps the ratio functions take (.ratio_boot_keep(), which settles
# the rule of R/dmar_bootstrap.R, then the limits). The mechanics are
# tested here, so the rule's counts are lowered for this file.
local_small_B()
boot_limits <- function(reps, est, alpha_lower, alpha_upper, method,
                        acc = 0) {
  rule <- DMAR:::.ratio_boot_rule(alpha_lower, alpha_upper,
                                  1 - alpha_lower - alpha_upper, method)
  kept <- DMAR:::.ratio_boot_keep(reps, rep(1, length(reps)), length(reps),
                                  rule, "ci_ratio_of_means")
  DMAR:::.ratio_boot_limits(kept, est, alpha_lower, alpha_upper, method,
                            acc = function() acc)$limits
}

# The shared numerics of ci_ratio_of_means() and ci_ratio_of_means_paired():
# the closed-form Fieller set, the root-finding set for degrees of freedom
# that change with the ratio, and the bootstrap pieces.

# The t statistic of the hypothesis mu_num - rho mu_den = 0.
fieller_t <- function(rho, num, den, V) {
  (num - rho * den) / sqrt(V[1, 1] - 2 * rho * V[1, 2] + rho^2 * V[2, 2])
}

# Fieller's (1954, p. 180) Cushny-Peebles summary: numerator 0.75,
# denominator 2.33, variances of the means 0.3201 and 0.4009, covariance
# 0.2848, so the denominator's t is 3.68 and the largest |t|, at the
# stationary point of t (Fieller, 1954, p. 177), is 4.53.
V_cp <- matrix(c(0.3201, 0.2848, 0.2848, 0.4009), 2)

test_that("a bounded Fieller set has limits where |t| equals the critical value", {
  q <- 2.5
  s <- DMAR:::.fieller_set(0.75, 2.33, V_cp, q, q)
  expect_identical(s$shape, "interval")
  expect_true(is.na(s$excluded_lower) && is.na(s$excluded_upper))
  expect_lt(s$lower, 0.75 / 2.33)
  expect_gt(s$upper, 0.75 / 2.33)
  expect_equal(fieller_t(s$lower, 0.75, 2.33, V_cp), q, tolerance = 1e-10)
  expect_equal(fieller_t(s$upper, 0.75, 2.33, V_cp), -q, tolerance = 1e-10)
})

test_that("the complement of an interval keeps its gap and the whole line has none", {
  # sqrt(13.54) < q < sqrt(20.51): the denominator is not significant but
  # some ratios are still rejected.
  q <- sqrt(15)
  s <- DMAR:::.fieller_set(0.75, 2.33, V_cp, q, q)
  expect_identical(s$shape, "complement of an interval")
  expect_identical(c(s$lower, s$upper), c(-Inf, Inf))
  expect_lt(s$excluded_lower, s$excluded_upper)
  expect_gt(s$excluded_lower, 0.75 / 2.33)
  expect_equal(abs(fieller_t(c(s$excluded_lower, s$excluded_upper),
                             0.75, 2.33, V_cp)), c(q, q), tolerance = 1e-10)
  # Inside the gap the test rejects; just outside it does not.
  mid <- mean(c(s$excluded_lower, s$excluded_upper))
  expect_gt(abs(fieller_t(mid, 0.75, 2.33, V_cp)), q)

  w <- DMAR:::.fieller_set(0.75, 2.33, V_cp, 5, 5)
  expect_identical(w$shape, "whole line")
  expect_identical(c(w$lower, w$upper), c(-Inf, Inf))
  expect_true(is.na(w$excluded_lower))
})

test_that("negating the numerator and denominator together leaves the set unchanged", {
  for (q in c(2, sqrt(15), 5)) {
    expect_equal(DMAR:::.fieller_set(-0.75, -2.33, V_cp, q, q),
                 DMAR:::.fieller_set(0.75, 2.33, V_cp, q, q))
  }
  # Unequal tails keep their sides under the negation as well.
  expect_equal(DMAR:::.fieller_set(-0.75, -2.33, V_cp, 1.7, 3),
               DMAR:::.fieller_set(0.75, 2.33, V_cp, 1.7, 3))
})

test_that("each tail governs its own side of the estimate", {
  # A zero tail (an infinite critical value) leaves its side unbounded.
  one <- DMAR:::.fieller_set(0.75, 2.33, V_cp, 1.833, Inf)
  two <- DMAR:::.fieller_set(0.75, 2.33, V_cp, 1.833, 1.833)
  expect_identical(one$shape, "interval")
  expect_identical(one$upper, Inf)
  expect_equal(one$lower, two$lower)
  # A small upper tail whose critical value exceeds the denominator's t,
  # beside a bounded lower side, gives the union of two intervals.
  u <- DMAR:::.fieller_set(0.75, 2.33, V_cp, 2.4, 4.3)
  expect_identical(u$shape, "union of two intervals")
  expect_true(is.finite(u$lower))
  expect_identical(u$upper, Inf)
  expect_equal(fieller_t(u$lower, 0.75, 2.33, V_cp), 2.4, tolerance = 1e-10)
  expect_equal(fieller_t(c(u$excluded_lower, u$excluded_upper), 0.75, 2.33,
                         V_cp), c(-4.3, -4.3), tolerance = 1e-10)
})

test_that("the quadratic roots stay accurate when one root is tiny", {
  r <- DMAR:::.quadratic_roots(1, -1e8, 1, 1e16 - 4)
  expect_equal(r[1L], 1e-8, tolerance = 1e-12)
  expect_equal(r[2L], 1e8, tolerance = 1e-12)
})

test_that("the root-finding engine reproduces the closed form at fixed degrees of freedom", {
  # With every group on infinite degrees of freedom the Satterthwaite
  # degrees of freedom are infinite at every ratio, the critical values are
  # normal ones, and the two engines must agree in every shape.
  v <- c(0.04, 0.09, 0.05)
  m <- c(1.3, 1.0, 0.6)
  for (alphas in list(c(0.025, 0.025), c(0.0005, 0.0005), c(0.05, 0),
                      c(0.01, 0.04))) {
    for (den in c(0.4, 1.0)) {
      mm <- m
      mm[2L] <- mm[3L] + den
      V <- matrix(c(v[1] + v[3], v[3], v[3], v[2] + v[3]), 2)
      closed <- DMAR:::.fieller_set(mm[1] - mm[3], mm[2] - mm[3], V,
                                    DMAR:::.ratio_crit(alphas[1], Inf),
                                    DMAR:::.ratio_crit(alphas[2], Inf))
      roots <- DMAR:::.fieller_set_df(mm, v, rep(Inf, 3), alphas[1],
                                      alphas[2])
      expect_identical(roots$shape, closed$shape)
      expect_equal(unlist(roots[1:4]), unlist(closed[1:4]),
                   tolerance = 1e-9)
    }
  }
})

test_that("the unequal-variance set inverts the Welch test at every finite end", {
  set.seed(3)
  x <- rnorm(15, 5, 2)
  y <- rnorm(12, 0.8, 2)
  m <- c(mean(x), mean(y))
  v <- c(var(x) / 15, var(y) / 12)
  for (level in c(0.80, 0.95, 0.99)) {
    a <- (1 - level) / 2
    s <- DMAR:::.fieller_set_df(m, v, c(14, 11), a, a)
    ends <- c(s$lower, s$upper, s$excluded_lower, s$excluded_upper)
    ends <- ends[is.finite(ends)]
    expect_gt(length(ends), 0L)
    p <- vapply(ends, function(L) stats::t.test(x, L * y)$p.value, 0)
    expect_equal(p, rep(1 - level, length(ends)), tolerance = 1e-8)
  }
})

test_that("bootstrap limits from a zero tail are infinite on that side", {
  set.seed(1)
  b <- rnorm(500, 1, 0.1)
  lim <- boot_limits(b, 1, 0.05, 0, "percentile")
  expect_identical(lim[2L], Inf)
  expect_equal(lim[1L], unname(stats::quantile(b, 0.05)))
  lim2 <- boot_limits(b, 1, 0, 0.05, "bca", acc = 0.01)
  expect_identical(lim2[1L], -Inf)
  expect_true(is.finite(lim2[2L]))
})

test_that("the jackknife accelerations match boot's influence values", {
  skip_if_not_installed("boot")
  stat_groups <- function(d, i) {
    d <- d[i, ]
    (mean(d$v[d$g == 1]) - mean(d$v[d$g == 3])) /
      (mean(d$v[d$g == 2]) - mean(d$v[d$g == 3]))
  }
  set.seed(9)
  x <- rexp(20, 1 / 10)
  y <- rexp(45, 1 / 6)
  z <- rexp(31, 1 / 2)
  d <- data.frame(v = c(x, y, z), g = rep(1:3, c(20, 45, 31)))
  bo <- boot::boot(d, stat_groups, R = 10, strata = d$g)
  L <- boot::empinf(bo, type = "jack")
  ng <- c(20, 45, 31)[d$g]
  expect_equal(DMAR:::.ratio_acceleration_groups(x, y, z),
               sum(L^3 / ng^3) / (6 * sum(L^2 / ng^2)^1.5),
               tolerance = 1e-10)

  a <- rnorm(40, 3)
  b <- rnorm(40, 5)
  bp <- boot::boot(data.frame(a, b), function(d, i) sum(d$a[i]) / sum(d$b[i]),
                   R = 10)
  Lp <- boot::empinf(bp, type = "jack")
  expect_equal(DMAR:::.ratio_acceleration_totals(a, b),
               sum(Lp^3) / (6 * sum(Lp^2)^1.5), tolerance = 1e-10)
})

test_that("the BCa levels agree with boot.ci() on the same replicates", {
  skip_if_not_installed("boot")
  # Equal group sizes, where boot.ci()'s unscaled acceleration and the
  # size-scaled one coincide.
  set.seed(5)
  x <- rexp(30, 1 / 10)
  y <- rexp(30, 1 / 6)
  B <- 999L
  reps <- DMAR:::.ratio_boot_means(x, B) / DMAR:::.ratio_boot_means(y, B)
  est <- mean(x) / mean(y)
  acc <- DMAR:::.ratio_acceleration_groups(x, y)
  d <- data.frame(v = c(x, y), g = rep(1:2, c(30, 30)))
  bo <- boot::boot(d, function(d, i) {
    d <- d[i, ]
    mean(d$v[d$g == 1]) / mean(d$v[d$g == 2])
  }, R = 10, strata = d$g)
  bo$t <- matrix(reps, ncol = 1L)
  bo$R <- B
  L <- boot::empinf(bo, type = "jack")
  ci <- boot::boot.ci(bo, conf = 0.95, type = c("bca", "perc"), L = L)
  # boot.ci() reports the order statistic positions (B + 1) p it uses,
  # rounded to two places.
  z0 <- stats::qnorm(mean(reps < est))
  zq <- stats::qnorm(c(0.025, 0.975))
  adj <- stats::pnorm(z0 + (z0 + zq) / (1 - acc * (z0 + zq)))
  expect_equal(round(adj * (B + 1), 2), unname(ci$bca[2:3]))
  # The limits differ only by the quantile rule (stats::quantile() against
  # boot's interpolation between order statistics): within one spacing.
  sorted <- sort(reps)
  spacing <- function(pos) {
    max(vapply(floor(pos), function(k) max(diff(sorted[(k - 1):(k + 2)])),
               numeric(1L)))
  }
  ours <- boot_limits(reps, est, 0.025, 0.025, "bca", acc)
  expect_lt(max(abs(ours - unname(ci$bca[4:5]))), 2 * spacing(adj * (B + 1)))
  pct <- boot_limits(reps, est, 0.025, 0.025, "percentile")
  expect_lt(max(abs(pct - unname(ci$percent[4:5]))),
            2 * spacing(c(0.025, 0.975) * (B + 1)))
})
