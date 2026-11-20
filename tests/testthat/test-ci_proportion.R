# Tests for ci_proportion(). Oracle: the closed form of the Wilson (1927)
# score interval, written out here independently of the package code,
#   center = (x + z^2 / 2) / (N + z^2),
#   half   = z / (N + z^2) * sqrt(x (N - x) / N + z^2 / 4),
# and the Wald interval p_hat +/- z sqrt(p_hat (1 - p_hat) / N).

wilson_oracle <- function(x, N, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level) / 2)
  center <- (x + z^2 / 2) / (N + z^2)
  half <- z / (N + z^2) * sqrt(x * (N - x) / N + z^2 / 4)
  c(center - half, center + half)
}

limits <- function(out) {
  c(out$value[out$term == "ci_lower"], out$value[out$term == "ci_upper"])
}

test_that("ci_proportion() reproduces the closed-form Wilson interval", {
  for (x in c(0, 3, 7, 10)) {
    expect_equal(limits(ci_proportion(x, 10)), wilson_oracle(x, 10),
                 tolerance = 1e-12)
  }
  expect_equal(limits(ci_proportion(17, 40, conf_level = 0.90)),
               wilson_oracle(17, 40, conf_level = 0.90), tolerance = 1e-12)
})

test_that("the Wilson interval has positive width at 0 and N successes and the Wald interval does not (help claim)", {
  # The help page says the Wilson interval stays inside [0, 1] and has
  # positive width when the number of successes is 0 or N, unlike the Wald
  # interval. At x = 0 the Wilson upper limit is z^2 / (N + z^2).
  z <- qnorm(0.975)
  w0 <- limits(ci_proportion(0, 10))
  wN <- limits(ci_proportion(10, 10))
  expect_equal(w0, c(0, z^2 / (10 + z^2)), tolerance = 1e-12)
  expect_equal(wN, c(10 / (10 + z^2), 1), tolerance = 1e-12)
  expect_gt(diff(w0), 0)
  expect_gt(diff(wN), 0)
  expect_equal(diff(limits(ci_proportion(0, 10, method = "wald"))), 0)
  expect_equal(diff(limits(ci_proportion(10, 10, method = "wald"))), 0)
  # Inside [0, 1] exactly (see the next test).
  for (x in 0:10) {
    w <- limits(ci_proportion(x, 10))
    expect_gte(w[1], 0)
    expect_lte(w[2], 1)
  }
  # The Wald interval leaves [0, 1] at a small count; the Wilson does not.
  wald1 <- limits(ci_proportion(1, 10, method = "wald"))
  expect_lt(wald1[1], 0)
  expect_equal(wald1, 0.1 + c(-1, 1) * z * sqrt(0.1 * 0.9 / 10),
               tolerance = 1e-12)
})

test_that("the Wilson limits are 0 exactly at 0 successes and 1 exactly at N", {
  # Through the first builds of 1.1.0 the lower limit at 0 successes was a
  # difference of two equal terms and came back a few ulps either side of
  # 0 (-5.6e-17 to 2.8e-17), and the upper limit at N successes a few ulps
  # either side of 1, over N = 1 to 500 and these five levels. The exact
  # values are those of the closed form: x = 0 gives a center equal to the
  # half width, so the lower limit is 0, and x = N gives center + half = 1.
  # Every other limit agrees with the oracle and with stats::prop.test()
  # without continuity correction, the same score interval.
  for (conf_level in c(.80, .90, .95, .99, .999)) {
    for (N in c(1:60, 97, 250, 500)) {
      lo <- limits(ci_proportion(0, N, conf_level))
      hi <- limits(ci_proportion(N, N, conf_level))
      expect_identical(lo[1], 0)
      expect_identical(hi[2], 1)
      expect_equal(lo[2], wilson_oracle(0, N, conf_level)[2], tolerance = 1e-12)
      expect_equal(hi[1], wilson_oracle(N, N, conf_level)[1], tolerance = 1e-12)
    }
  }
  for (x in c(1, 5, 13, 29)) {
    w <- limits(ci_proportion(x, 30, conf_level = .95))
    ref <- stats::prop.test(x, 30, conf.level = .95, correct = FALSE)$conf.int
    expect_equal(w, as.numeric(ref), tolerance = 1e-10)
    expect_true(w[1] > 0 && w[2] < 1)
  }
})

test_that("the coverage figures stated on ?ci_proportion", {
  # The help page gives, for a 95% interval, the coverage averaged over a
  # uniformly distributed true proportion (Wilson .954 at N = 10 and .953 at
  # N = 20; Wald .769 and .846) and says the Wilson coverage falls to about
  # .84 for some true proportions near 0 or 1. Oracles, with no package code:
  # the average coverage in closed form,
  #   int_0^1 sum_x dbinom(x, N, pi) 1{L_x <= pi <= U_x} dpi
  #     = sum_x [pbeta(U_x; x + 1, N - x + 1) - pbeta(L_x; x + 1, N - x + 1)] / (N + 1),
  # and the exact coverage sum_x dbinom(x, N, pi) 1{L_x <= pi <= U_x} at
  # points just outside each limit, where the infimum is approached.
  wald_oracle <- function(x, N) {
    z <- qnorm(0.975); p <- x / N
    c(p - z * sqrt(p * (1 - p) / N), p + z * sqrt(p * (1 - p) / N))
  }
  lims <- function(f, N) t(sapply(0:N, f, N = N))
  avg_cov <- function(L, U, N) {
    x <- 0:N
    sum(pbeta(pmin(pmax(U, 0), 1), x + 1, N - x + 1) -
          pbeta(pmin(pmax(L, 0), 1), x + 1, N - x + 1)) / (N + 1)
  }
  for (N in c(10, 20)) {
    w <- lims(wilson_oracle, N)
    pk <- t(sapply(0:N, function(x) limits(ci_proportion(x, N))))
    expect_equal(pk, w, tolerance = 1e-12)
    d <- lims(wald_oracle, N)
    stated <- if (N == 10) c(.954, .769) else c(.953, .846)
    expect_equal(round(avg_cov(w[, 1], w[, 2], N), 3), stated[1])
    expect_equal(round(avg_cov(d[, 1], d[, 2], N), 3), stated[2])
    # Infimum of the Wilson coverage, approached just outside a limit.
    e <- sort(unique(c(w)))
    pts <- c(e - 1e-9, e + 1e-9)
    pts <- pts[pts > 0 & pts < 1]
    cover <- sapply(pts, function(pp)
      sum(dbinom(0:N, N, pp) * (w[, 1] <= pp & pp <= w[, 2])))
    # .8350 at N = 10 and .8366 at N = 20: "about .84" on the page.
    expect_gt(min(cover), 0.83)
    expect_lt(min(cover), 0.845)
    expect_true(min(abs(pts[which.min(cover)] - c(0, 1))) < 0.05)
  }
})
