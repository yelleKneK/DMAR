# ci_proportion_paired() and the paired proportion internals it shares
# with equivalence_proportion_paired().

lims_of <- function(res) {
  v <- stats::setNames(res$value, res$term)
  c(v[["ci_lower"]], v[["ci_upper"]])
}

test_that("ci_proportion_paired() returns the documented rows and attributes", {
  res <- ci_proportion_paired(n_10 = 32, n_01 = 27, N = 344)
  expect_s3_class(res, "dmar_tbl")
  expect_identical(res$term, c("ci_lower", "difference", "ci_upper",
                               "n_10", "n_01", "N"))
  expect_type(res$value, "double")
  expect_equal(res$value[res$term == "difference"], 5 / 344)
  expect_equal(res$value[4:6], c(32, 27, 344))
  expect_identical(attr(res, "conf_level"), 0.95)
  expect_identical(attr(res, "method"), "score")
  expect_identical(attr(ci_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
                                             method = "wald"), "method"),
                   "wald")
})

test_that("the restricted estimate of p_01 maximizes the restricted likelihood", {
  # Under p_10 - p_01 = psi0 the log likelihood in q = p_01 is
  # n_10 log(q + psi0) + n_01 log(q) + (n - n_10 - n_01) log(1 - 2q - psi0).
  loglik <- function(q, n_10, n_01, n, psi0)
    n_10 * log(q + psi0) + n_01 * log(q) +
      (n - n_10 - n_01) * log(1 - 2 * q - psi0)
  cases <- list(c(38, 52, 400, -0.05), c(38, 52, 400, 0.03),
                c(12, 5, 100, -0.10), c(3, 1, 20, 0.20), c(7, 9, 30, -0.30),
                c(32, 27, 344, -0.05))
  for (cs in cases) {
    n_10 <- cs[1]; n_01 <- cs[2]; n <- cs[3]; psi0 <- cs[4]
    lo <- max(0, -psi0) + 1e-10
    hi <- (1 - psi0) / 2 - 1e-10
    num <- stats::optimize(loglik, c(lo, hi), n_10 = n_10, n_01 = n_01,
                           n = n, psi0 = psi0, maximum = TRUE,
                           tol = 1e-12)$maximum
    expect_equal(.paired_rd_p01(n_10, n_01, n, psi0), num, tolerance = 1e-6)
  }
})

test_that("the restricted estimate stays in the parameter space", {
  # p_01 >= max(0, -psi0) keeps p_10 = p_01 + psi0 nonnegative, and
  # p_01 <= (1 - psi0) / 2 keeps the concordant probability nonnegative.
  # The cancellation-free root form must hold this at every psi0,
  # including where B > 0 and the textbook form loses digits.
  grid <- seq(-0.99, 0.99, by = 0.01)
  for (cs in list(c(38, 52, 400), c(0, 0, 50), c(0, 7, 30), c(40, 0, 40),
                  c(1, 0, 1000), c(0, 1, 1000))) {
    p01 <- .paired_rd_p01(cs[1], cs[2], cs[3], grid)
    expect_true(all(p01 >= pmax(0, -grid) - 1e-12))
    expect_true(all(p01 <= (1 - grid) / 2 + 1e-12))
  }
})

test_that("Nam's (1997) restricted estimates are reproduced", {
  # Nam (1997, p. 1427): the liver scan data, n = 344 with 32 and 27
  # discordant pairs, under psi0 = -.05 give restricted estimates
  # p_01 = .1179 and p_10 = .0679.
  p01 <- .paired_rd_p01(32, 27, 344, -0.05)
  expect_equal(round(p01, 4), 0.1179)
  expect_equal(round(p01 - 0.05, 4), 0.0679)
})

test_that("the score statistic decreases in psi0 and the limits solve Z = +/- q", {
  for (cs in list(c(38, 52, 400), c(12, 5, 100), c(3, 0, 20), c(0, 7, 30),
                  c(0, 0, 50))) {
    grid <- seq(-0.99, 0.99, by = 0.01)
    z <- .paired_rd_score(cs[1], cs[2], cs[3], grid)
    expect_true(all(diff(z) < 0))
    res <- ci_proportion_paired(n_10 = cs[1], n_01 = cs[2], N = cs[3],
                                conf_level = 0.90)
    l <- lims_of(res)
    q <- stats::qnorm(0.95)
    expect_equal(.paired_rd_score(cs[1], cs[2], cs[3], l[1]), q,
                 tolerance = 1e-8)
    expect_equal(.paired_rd_score(cs[1], cs[2], cs[3], l[2]), -q,
                 tolerance = 1e-8)
  }
})

test_that("the score interval reproduces PropCIs::scoreci.mp", {
  # Pinned from PropCIs::scoreci.mp(b = n_01, c = n_10, n, conf.level)
  # (PropCIs 0.3.0, 2026-11-05), Tango's (1998) interval; live comparison
  # over a wider sweep in tools/oracle_checks.R. PropCIs stops its
  # bisection at a step of 1e-7, and in these tables its limits sit within
  # 7e-8 of DMAR's, whose limits solve Z = +/- q (checked above); the
  # tolerance is set above that and far below any difference that would
  # matter. The last case is the example on PropCIs's own help page.
  ref <- rbind(
    c(38, 52, 400, 0.90, -0.074377383,  0.004056882),
    c(38, 52, 400, 0.95, -0.082102288,  0.011661878),
    c(12,  5, 100, 0.90,  0.002337099,  0.142046661),
    c( 0,  0,  50, 0.95, -0.071347594,  0.071347654),
    c( 3,  0,  20, 0.90,  0.012968540,  0.321828046),
    c( 0,  7,  30, 0.90, -0.379026000, -0.131306730),
    c( 1,  9,  15, 0.90, -0.743214194, -0.210461104),
    c(40,  0,  40, 0.95,  0.824756861,  1.000000000),
    c( 0, 40,  40, 0.95, -1.000000000, -0.824756801),
    c(32, 27, 344, 0.95, -0.029898173,  0.059375877),
    c(20, 40, 160, 0.95, -0.218064234, -0.030864589))
  for (i in seq_len(nrow(ref))) {
    r <- ref[i, ]
    got <- lims_of(ci_proportion_paired(n_10 = r[1], n_01 = r[2], N = r[3],
                                        conf_level = r[4]))
    expect_lt(max(abs(got - r[5:6])), 2e-7)
  }
})

test_that("the limits reach the edges of [-1, 1] only when the data do", {
  expect_identical(lims_of(ci_proportion_paired(n_10 = 40, n_01 = 0,
                                                N = 40))[2], 1)
  expect_identical(lims_of(ci_proportion_paired(n_10 = 0, n_01 = 40,
                                                N = 40))[1], -1)
  l <- lims_of(ci_proportion_paired(n_10 = 39, n_01 = 0, N = 40))
  expect_true(l[1] > -1 && l[2] < 1)
})

test_that("with no discordant pairs the limits are -/+ z^2 / (n + z^2)", {
  z <- stats::qnorm(0.975)
  expect_equal(lims_of(ci_proportion_paired(n_10 = 0, n_01 = 0, N = 50)),
               c(-1, 1) * z^2 / (50 + z^2), tolerance = 1e-10)
  # The Wald interval collapses to a point there.
  expect_identical(lims_of(ci_proportion_paired(n_10 = 0, n_01 = 0, N = 50,
                                                method = "wald")),
                   c(0, 0))
})

test_that("the score interval excludes 0 exactly when McNemar's test rejects", {
  # At psi0 = 0 the score statistic is McNemar's (n_10 - n_01) /
  # sqrt(n_10 + n_01), whose square is stats::mcnemar.test()'s statistic
  # without the continuity correction.
  for (cs in list(c(32, 27, 344), c(12, 5, 100), c(20, 40, 160),
                  c(9, 2, 60), c(15, 5, 80))) {
    tab <- matrix(c(cs[3] - cs[1] - cs[2], cs[2], cs[1], 0), nrow = 2)
    mc <- stats::mcnemar.test(tab, correct = FALSE)
    expect_equal(.paired_rd_score(cs[1], cs[2], cs[3], 0)^2,
                 unname(mc$statistic), tolerance = 1e-12)
    l <- lims_of(ci_proportion_paired(n_10 = cs[1], n_01 = cs[2], N = cs[3]))
    expect_identical(l[1] > 0 || l[2] < 0, mc$p.value < 0.05)
  }
})

test_that("the Wald option follows its formula", {
  w <- lims_of(ci_proportion_paired(n_10 = 38, n_01 = 52, N = 400,
                                    method = "wald"))
  hw <- stats::qnorm(0.975) * sqrt((90 / 400 - (-14 / 400)^2) / 400)
  expect_equal(w, -14 / 400 + c(-1, 1) * hw)
})

test_that("exact coverage favors the score interval where the Wald interval fails", {
  skip_on_cran()
  # Coverage summed over every possible table (no simulation), computed
  # without the package's own statistic or limits. The score interval is
  # the set of differences the score test does not reject, so it covers
  # the true psi exactly when |Z(psi)| <= q, and Z(psi) is computed here
  # from Nam's (1997) equation 3 with the restricted estimate of p_01
  # found by optimize() on the restricted likelihood. The help page
  # quotes .986 and .779 for p_10 = .03, p_01 = .01, n = 50; the limits
  # of PropCIs::scoreci.mp() (PropCIs 0.3.0) give the same .98612 in
  # tools/oracle_snippets/propcis_paired.R. The second configuration has
  # more discordant pairs.
  restricted_p01 <- function(n_10, n_01, n, psi0) {
    loglik <- function(q) n_10 * log(q + psi0) + n_01 * log(q) +
      (n - n_10 - n_01) * log(1 - 2 * q - psi0)
    stats::optimize(loglik, c(max(0, -psi0) + 1e-12, (1 - psi0) / 2 - 1e-12),
                    maximum = TRUE, tol = 1e-12)$maximum
  }
  exact_coverage <- function(p10, p01, n, level = 0.95) {
    psi <- p10 - p01
    q <- stats::qnorm(1 - (1 - level) / 2)
    g <- expand.grid(n_10 = 0:n, n_01 = 0:n)
    g <- g[g$n_10 + g$n_01 <= n, ]
    pr <- mapply(function(a, b) {
      stats::dmultinom(c(a, b, n - a - b),
                       prob = c(p10, p01, 1 - p10 - p01))
    }, g$n_10, g$n_01)
    z <- mapply(function(a, b) {
      (a - b - n * psi) /
        sqrt(n * (2 * restricted_p01(a, b, n, psi) + psi * (1 - psi)))
    }, g$n_10, g$n_01)
    score <- abs(z) <= q
    # Table by table, the package's limits cover psi exactly when the
    # independently computed test does not reject it.
    ours <- mapply(function(a, b) {
      l <- lims_of(ci_proportion_paired(n_10 = a, n_01 = b, N = n,
                                        conf_level = level))
      l[1] <= psi && psi <= l[2]
    }, g$n_10, g$n_01)
    expect_identical(ours, score)
    est <- (g$n_10 - g$n_01) / n
    hw <- q * sqrt(pmax((g$n_10 + g$n_01) / n - est^2, 0) / n)
    c(score = sum(pr * score),
      wald = sum(pr * (est - hw <= psi & psi <= est + hw)))
  }
  a <- exact_coverage(0.03, 0.01, 50)
  expect_equal(round(a[["score"]], 3), 0.986)
  expect_equal(round(a[["wald"]], 3), 0.779)
  b <- exact_coverage(0.10, 0.02, 40)
  expect_gt(b[["score"]], 0.95)
  expect_lt(b[["wald"]], 0.92)
})

test_that("raw vectors and counts give the same result", {
  set.seed(113)
  y <- stats::rbinom(200, 1, 0.70)
  x <- ifelse(stats::runif(200) < 0.15, 1 - y, y)
  a <- ci_proportion_paired(x = x, y = y)
  b <- ci_proportion_paired(n_10 = sum(x == 1 & y == 0),
                            n_01 = sum(x == 0 & y == 1), N = 200)
  expect_equal(a$value, b$value)
  # Logical vectors are read as 1 and 0.
  expect_equal(ci_proportion_paired(x = x == 1, y = y == 1)$value, a$value)
  # A pair with a missing member is dropped as a pair.
  x_na <- c(x, NA, 1)
  y_na <- c(y, 1, NA)
  expect_equal(ci_proportion_paired(x = x_na, y = y_na)$value, a$value)
})

test_that("ci_proportion_paired() refuses malformed input", {
  expect_error(ci_proportion_paired(n_10 = 5, n_01 = 7, N = 10), "at least")
  expect_error(ci_proportion_paired(n_10 = 5, n_01 = 7), "Supply both")
  expect_error(ci_proportion_paired(n_10 = Inf, n_01 = 7, N = 10),
               "non-negative integer")
  expect_error(ci_proportion_paired(n_10 = 2.5, n_01 = 7, N = 10),
               "non-negative integer")
  expect_error(ci_proportion_paired(x = c(1, 0, 2), y = c(1, 1, 0)),
               "only 0 and 1")
  expect_error(ci_proportion_paired(x = c("1", "0"), y = c("1", "1")),
               "numeric or logical")
  expect_error(ci_proportion_paired(x = c(1, 0, 1), y = c(1, 1)),
               "same length")
  expect_error(ci_proportion_paired(x = c(1, 0)), "Supply both")
  expect_error(ci_proportion_paired(x = c(1, 0), y = c(0, 1), N = 2),
               "not both")
  expect_error(ci_proportion_paired(x = c(NA, 1), y = c(1, NA)),
               "At least one unit")
  expect_error(ci_proportion_paired(n_10 = 1, n_01 = 2, N = 10,
                                    conf_level = 95), "conf_level")
  expect_error(ci_proportion_paired(n_10 = 1, n_01 = 2, N = 10,
                                    method = "exact"))
})
