test_that("ss_aipe_semipartial_r() plans from the first-order variance at R2_full", {
  N <- function(...) {
    res <- ss_aipe_semipartial_r(...)
    res$value[res$term == "necessary_N"]
  }
  # r_sp = .8, p = 2, width = .10: V = .2136 at R2_full = .94 and .1296 at
  # R2_full = r_sp^2; without R2_full the largest V, r^2 (1 - r^2) = .2304.
  expect_equal(N(true_r_sp = 0.8, p = 2, width = 0.10), 357)
  expect_equal(N(true_r_sp = 0.8, p = 2, width = 0.10, true_R2_full = 0.94), 331)
  expect_equal(N(true_r_sp = 0.8, p = 2, width = 0.10, true_R2_full = 0.64), 202)
  expect_equal(N(true_r_sp = 0.8, p = 2, width = 0.10, true_R2_full = 0.8^2), 202)
  # r_sp = .2, p = 4: the largest V is (1 - r^2)^2, at R2_full = r_sp^2.
  expect_equal(N(true_r_sp = 0.2, p = 4, width = 0.10), 1421)
  expect_equal(N(true_r_sp = 0.2, p = 4, width = 0.10, true_R2_full = 0.54), 714)
  # The closed form, computed here.
  z <- stats::qnorm(0.975)
  V <- (1 - 0.2^2)^2 - (0.54 - 0.2^2) * (1 - 2 * 0.2^2)
  expect_equal(N(true_r_sp = 0.2, p = 4, width = 0.10, true_R2_full = 0.54),
               3 + 1 + ceiling(z^2 * V / 0.05^2))

  # Expected width uses the same variance at the returned N.
  res <- ss_aipe_semipartial_r(true_r_sp = 0.8, p = 2, width = 0.10,
                               true_R2_full = 0.94)
  expect_equal(res$value[res$term == "expected_width"],
               2 * z * sqrt(0.2136 / (331 - 1 - 1)), tolerance = 1e-10)
  expect_lte(res$value[res$term == "expected_width"], 0.10)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_aipe_semipartial_r() returns documented rows and width <= target", {
  res <- ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15)
  expect_true(all(c("necessary_N", "expected_width", "true_r_sp", "p",
                    "width", "conf_level") %in% res$term))
  expect_lte(res$value[res$term == "expected_width"], 0.15 + 1e-6)
})

# ---- The assurance step (a simulation) ------------------------------------

test_that("ss_aipe_semipartial_r() with assurance never plans fewer cases and appends its three rows", {
  res_50 <- ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15)
  res_80 <- ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15,
                                   assurance = 0.80, G = 2000, seed = 113)
  expect_gt(res_80$value[1], res_50$value[1])
  expect_identical(res_80$term, c(res_50$term, "assurance",
                                  "assurance_attained", "G"))
  expect_equal(res_80$value[res_80$term == "assurance"], 0.80)
  expect_equal(res_80$value[res_80$term == "G"], 2000)
  expect_type(res_80$value, "double")
  p <- res_80$value[res_80$term == "assurance_attained"]
  expect_equal(attr(res_80, "assurance_mc_se"), sqrt(p * (1 - p) / 2000))
  # Without assurance there is no row, no attribute, and no draw.
  expect_null(attr(res_50, "assurance_mc_se"))
  set.seed(4)
  before <- .Random.seed
  ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15, G = 10)
  expect_identical(.Random.seed, before)
})

test_that("a seeded ss_aipe_semipartial_r() assurance plan is reproducible and leaves the caller's stream alone", {
  skip_on_cran()
  plan <- function(seed) ss_aipe_semipartial_r(true_r_sp = 0.5, p = 2, width = 0.2,
                                               assurance = 0.90,
                                               true_R2_full = 0.25, G = 2000,
                                               seed = seed)
  set.seed(99)
  before <- .Random.seed
  a <- plan(113)
  expect_identical(.Random.seed, before)
  b <- plan(113)
  expect_identical(a, b)
  # Unseeded, the draws come from the session's stream and advance it.
  set.seed(99)
  c1 <- ss_aipe_semipartial_r(true_r_sp = 0.5, p = 2, width = 0.2,
                              assurance = 0.90, true_R2_full = 0.25, G = 2000)
  expect_false(identical(.Random.seed, before))
  set.seed(99)
  c2 <- ss_aipe_semipartial_r(true_r_sp = 0.5, p = 2, width = 0.2,
                              assurance = 0.90, true_R2_full = 0.25, G = 2000)
  expect_identical(c1, c2)
})

test_that("ss_aipe_semipartial_r() attains about the requested assurance in fresh draws", {
  skip_on_cran()
  # The plan for this condition is 253 or 254 at G = 20000; the plan
  # without assurance is 219.
  res <- ss_aipe_semipartial_r(true_r_sp = 0.5, p = 2, width = 0.2,
                               assurance = 0.90, true_R2_full = 0.25, G = 20000,
                               seed = 5150)
  n <- res$value[1]
  expect_gte(n, 251)
  expect_lte(n, 256)
  p <- res$value[res$term == "assurance_attained"]
  se <- attr(res, "assurance_mc_se")
  expect_gt(p, 0.90 - 4 * se)
  # Several conditions, a smaller G: within 4 standard errors or above.
  # The third entry of a cell is the number of variables partialed out.
  cells <- list(c(-.6, .70, 4, .25, .80), c(.2, .04, 2, .30, .95),
                c(.8, .90, 2, .15, .80), c(.3, NA, 2, .20, .80),
                c(.85, NA, 1, .15, .80))
  for (cc in cells) {
    res <- suppressWarnings(ss_aipe_semipartial_r(
      true_r_sp = cc[1], p = cc[3] + 1, width = cc[4], assurance = cc[5],
      true_R2_full = if (is.na(cc[2])) NULL else cc[2], G = 4000, seed = 31))
    p <- res$value[res$term == "assurance_attained"]
    expect_gt(p, cc[5] - 4 * attr(res, "assurance_mc_se"),
              label = paste(cc, collapse = " "))
    n_0 <- suppressWarnings(ss_aipe_semipartial_r(
      true_r_sp = cc[1], p = cc[3] + 1, width = cc[4],
      true_R2_full = if (is.na(cc[2])) NULL else cc[2]))$value[1]
    expect_gte(res$value[1], n_0)
  }
})

test_that("ss_aipe_semipartial_r() treats width as a half-width under which_width = 'Lower' in the assurance step too", {
  full  <- ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15,
                                 assurance = 0.80, G = 2000, seed = 113)
  lower <- ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.075,
                                 which_width = "Lower",
                                 assurance = 0.80, G = 2000, seed = 113)
  expect_identical(lower$value[1], full$value[1])
  expect_identical(lower$value[lower$term == "assurance_attained"],
                   full$value[full$term == "assurance_attained"])
})

test_that("the attained assurance is the share in fresh draws, above or below the request", {
  skip_on_cran()
  # At a small G the search is noisy and the size it returns can be one at
  # which the probability is below the request; the value reported is then
  # below it too, and agrees with an independent large draw at that size.
  res <- ss_aipe_semipartial_r(true_r_sp = 0.5, p = 2, width = 0.2,
                               assurance = 0.90, true_R2_full = 0.25,
                               G = 400, seed = 11)
  n <- res$value[1]
  attained <- res$value[res$term == "assurance_attained"]
  expect_equal(n, 250)
  expect_equal(attained, 0.855)
  expect_lt(attained, 0.90)
  set.seed(77)
  own <- DMAR:::.semipartial_r_half_width_draws(
    200000, n, 0.5, 0, 1, stats::qnorm(0.975))
  truth <- mean(own <= 0.1)
  expect_lt(truth, 0.90)
  expect_lt(abs(attained - truth),
            4 * sqrt(attr(res, "assurance_mc_se")^2 +
                       truth * (1 - truth) / 200000))
})

test_that("without true_R2_full and r_sp squared above 1/2 the simulated population has a full-model R2 of 1", {
  skip_on_cran()
  # The plan without assurance is 382. With a full-model R2 of 1 the
  # probability of an interval no wider than .10 is about .73 at 385 and
  # .82 at 386; with variables partialed out that explain nothing (q = 0)
  # it is already about .95 at 382, and the step would return 382.
  expect_equal(suppressWarnings(
    ss_aipe_semipartial_r(true_r_sp = .75, p = 3, width = .10))$value[1], 382)
  expect_equal(suppressWarnings(
    ss_aipe_semipartial_r(true_r_sp = .75, p = 3, width = .10,
                          assurance = .80, seed = 905))$value[1], 386)
})

test_that("ss_aipe_semipartial_r() refuses fewer than 100 simulated samples for the assurance step", {
  for (bad in c(1, 99)) {
    err <- expect_error(
      ss_aipe_semipartial_r(true_r_sp = .3, p = 2, width = .2,
                            assurance = .9, seed = 1, G = bad),
      class = "dmar_count_invalid")
    expect_identical(err$min, 100L)
    expect_match(conditionMessage(err),
                 "estimates a probability from 'G' simulated samples",
                 fixed = TRUE)
  }
  expect_no_error(ss_aipe_semipartial_r(true_r_sp = .3, p = 2, width = .2,
                                        assurance = .9, seed = 1, G = 100))
})

test_that("ss_aipe_semipartial_r() checks G and seed only when assurance is given", {
  expect_no_error(ss_aipe_semipartial_r(true_r_sp = 0.3, p = 3, width = 0.2,
                                        G = 4.5, seed = "a"))
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.3, p = 3, width = 0.2,
                                     assurance = 0.8, G = 4.5),
               class = "dmar_count_invalid")
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.3, p = 3, width = 0.2,
                                     assurance = 0.8, G = 0),
               class = "dmar_count_invalid")
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.3, p = 3, width = 0.2,
                                     assurance = 0.8, G = 100, seed = "a"),
               "'seed' must be a single number or NULL")
})

test_that("the draws of the semipartial assurance step are the statistics of a Wishart matrix", {
  skip_on_cran()
  # The step draws entries of the Cholesky factor of the Wishart matrix of
  # (Z_1, ..., Z_J, X, E). Here the matrix is formed from such a factor and
  # the sample semipartial correlation, the sample R2, and the half-width
  # are computed from it by matrix algebra.
  set.seed(1)
  J <- 3; n <- 25; r <- .4; q <- .3; z <- stats::qnorm(.975); d <- J + 2
  e <- sqrt(1 - r^2 - q); s <- sqrt(r^2 + e^2)
  for (i in 1:25) {
    L <- matrix(0, d, d)
    diag(L) <- sqrt(stats::rchisq(d, n - 1 - (0:(d - 1))))
    L[lower.tri(L)] <- stats::rnorm(d * (d - 1) / 2)
    A <- L %*% t(L)
    b <- c(sqrt(q), rep(0, J - 1), r, e)
    zi <- 1:J; xi <- J + 1
    s_yy <- sum(b * (A %*% b)); s_y <- A %*% b
    S_zz <- A[zi, zi]; s_yz <- s_y[zi]; s_xz <- A[xi, zi]
    s_xx <- A[xi, xi] - sum(s_xz * solve(S_zz, s_xz))
    s_xy <- s_y[xi] - sum(s_xz * solve(S_zz, s_yz))
    r_hat <- s_xy / sqrt(s_xx * s_yy)
    R2_hat <- r_hat^2 + sum(s_yz * solve(S_zz, s_yz)) / s_yy
    by_matrix <- z * sqrt(DMAR:::.semipartial_r_var_numerator(r_hat, R2_hat) /
                            (n - J - 1))
    by_entries <- DMAR:::.semipartial_r_half_width_from_bartlett(
      l_1 = L[1, 1], l_x = L[xi, xi], l_e = L[d, d], l_ex = L[d, xi],
      n_1 = (r * L[xi, 1] + e * L[d, 1]) / s,
      w = sum((r * L[xi, 2:J] + e * L[d, 2:J])^2) / s^2,
      r_sp = r, q = q, n = n, J = J, z = z)
    expect_equal(by_entries, by_matrix, tolerance = 1e-12)
  }
})

test_that("the draws of the semipartial assurance step agree with stats::rWishart() draws", {
  skip_on_cran()
  # The share of narrow intervals agrees with the same share computed from
  # stats::rWishart() draws (two independent sets of 40,000; the standard
  # error of the difference is about .0035).
  set.seed(2)
  z <- stats::qnorm(.975)
  n <- 90; J <- 2; r <- .6; q <- .2; hw <- .11
  e <- sqrt(1 - r^2 - q)
  W <- stats::rWishart(40000, n - 1, diag(J + 2))   # order X, Z_1, Z_2, E
  b <- c(r, sqrt(q), 0, e); zi <- 2:3
  half <- apply(W, 3, function(A) {
    s_yy <- sum(b * (A %*% b)); s_y <- A %*% b
    S_zz <- A[zi, zi]; s_yz <- s_y[zi]; s_xz <- A[1, zi]
    s_xx <- A[1, 1] - sum(s_xz * solve(S_zz, s_xz))
    s_xy <- s_y[1] - sum(s_xz * solve(S_zz, s_yz))
    r_hat <- s_xy / sqrt(s_xx * s_yy)
    R2_hat <- r_hat^2 + sum(s_yz * solve(S_zz, s_yz)) / s_yy
    z * sqrt(DMAR:::.semipartial_r_var_numerator(r_hat, R2_hat) / (n - J - 1))
  })
  own <- DMAR:::.semipartial_r_half_width_draws(40000, n, r, q, J, z)
  expect_lt(abs(mean(half <= hw) - mean(own <= hw)), 4 * .0035)
  expect_equal(mean(half), mean(own), tolerance = .002)
  expect_equal(stats::sd(half), stats::sd(own), tolerance = .02)
})

test_that("ss_aipe_semipartial_r_sensitivity() confirms the simulated assurance plan from data", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  # Two conditions with R2_full given and below 1 (the sensitivity function
  # refuses a population R2 of 1). G = 4000: the standard error of a
  # proportion near .8 is .0063, near .9 it is .0047.
  cells <- list(c(.5, .25, 1, .20, .90), c(-.6, .70, 4, .25, .80))
  for (cc in cells) {
    plan <- ss_aipe_semipartial_r(true_r_sp = cc[1], p = cc[3] + 1, width = cc[4],
                                  assurance = cc[5], true_R2_full = cc[2],
                                  seed = 2718)
    sens <- ss_aipe_semipartial_r_sensitivity(
      true_r_sp = cc[1], specified_N = plan$value[1], p = cc[3] + 1,
      width = cc[4], G = 4000, true_R2_full = cc[2], seed = 314)
    share <- sens$value[sens$term == "pct_ci_less_width"]
    expect_lt(abs(share - cc[5]), 3 * sqrt(cc[5] * (1 - cc[5]) / 4000),
              label = paste(cc, collapse = " "))
  }
})

test_that("ss_aipe_semipartial_r() rejects bad inputs", {
  expect_error(ss_aipe_semipartial_r(true_r_sp = 1.5, p = 4, width = 0.2), "in \\(-1, 1\\)")
  expect_error(ss_aipe_semipartial_r(true_r_sp = 1.5, p = 4, width = 0.2),
               "'true_r_sp' must be a single value in (-1, 1).", fixed = TRUE)
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.3, p = 1, width = 0.2),
               "whole number of at least 2")
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.3, p = 2.5, width = 0.2),
               "whole number of at least 2")
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.3, p = NA_real_, width = 0.2),
               "whole number of at least 2")
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.3, p = 0, width = 0.2),
               "whole number of at least 2")
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.3, width = 0.2),
               "'p' (the number of predictors) is required.", fixed = TRUE)
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.3, p = 3, width = 0.2,
                                     conf_level = c(0.9, 0.95)),
               "conf_level")
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.3, p = 3, width = 0.2,
                                     assurance = 0.4), "assurance")
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.3, p = 3, width = 0.2,
                                     true_R2_full = 1.2),
               "'true_R2_full' must be a single number in [0, 1].",
               fixed = TRUE)
})

test_that("ss_aipe_semipartial_r() without R2_full equals the 1.0.x plan when r_sp squared is at most 1/2", {
  z <- stats::qnorm(0.975)
  for (r in c(0.1, 0.3, 0.5, 0.7)) {
    for (J in c(1, 4)) {
      res <- ss_aipe_semipartial_r(true_r_sp = r, p = J + 1, width = 0.10)
      expect_equal(res$value[res$term == "necessary_N"],
                   J + 1 + ceiling(z^2 * (1 - r^2)^2 / 0.05^2),
                   info = paste("r_sp =", r, "J =", J))
      expect_equal(
        res$value[res$term == "necessary_N"],
        ss_aipe_semipartial_r(true_r_sp = r, p = J + 1, width = 0.10,
                              true_R2_full = r^2)$value[1])
    }
  }
})

test_that("ss_aipe_semipartial_r() sample size moves with R2_full as the variance does", {
  N_at <- function(r, R2) vapply(R2, function(v)
    ss_aipe_semipartial_r(true_r_sp = r, p = 3, width = 0.10, true_R2_full = v)$value[1],
    numeric(1))
  # Below r_sp^2 = 1/2 a larger R2_full never needs more; above, never fewer.
  n_low  <- N_at(0.3, seq(0.09, 0.99, by = 0.10))
  n_high <- N_at(0.8, seq(0.64, 0.99, by = 0.05))
  expect_true(all(diff(n_low) <= 0))
  expect_true(all(diff(n_high) >= 0))
  expect_lt(n_low[length(n_low)], n_low[1])
  expect_gt(n_high[length(n_high)], n_high[1])
  # The plan without R2_full is sufficient for every R2_full.
  expect_true(all(n_low  <= ss_aipe_semipartial_r(0.3, p = 3, 0.10)$value[1]))
  expect_true(all(n_high <= ss_aipe_semipartial_r(0.8, p = 3, 0.10)$value[1]))
})

test_that("ss_aipe_semipartial_r() reports true_R2_full as its last row", {
  res <- ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15)
  expect_equal(res$term, c("necessary_N", "expected_width", "true_r_sp", "p",
                           "width", "conf_level", "true_R2_full"))
  expect_identical(res$value[res$term == "true_r_sp"], 0.25)
  expect_identical(res$value[res$term == "p"], 4)
  expect_false(any(c("r_sp", "J", "R2_full", "n") %in% res$term))
  expect_true(is.na(res$value[res$term == "true_R2_full"]))
  res <- ss_aipe_semipartial_r(true_r_sp = 0.25, p = 4, width = 0.15,
                               true_R2_full = 0.42)
  expect_equal(res$term[nrow(res)], "true_R2_full")
  expect_equal(res$value[nrow(res)], 0.42)
  expect_type(res$value, "double")
})

test_that("ss_aipe_semipartial_r() refuses an R2_full below r_sp squared, by class", {
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.9, p = 2, width = 0.2,
                                     true_R2_full = 0.1),
               class = "dmar_semipartial_r_incoherent")
  expect_error(ss_aipe_semipartial_r(true_r_sp = 0.9, p = 2, width = 0.2,
                                     true_R2_full = 0.1),
               class = "dmar_refusal")
})

test_that("ss_aipe_semipartial_r() plans without R2_full above r_sp squared of 1/2 without a warning", {
  expect_no_warning(ss_aipe_semipartial_r(true_r_sp = 0.7, p = 2, width = 0.10))
  expect_no_warning(ss_aipe_semipartial_r(true_r_sp = 0.8, p = 2, width = 0.10,
                                          true_R2_full = 0.94))
  expect_no_warning(
    res <- ss_aipe_semipartial_r(true_r_sp = 0.8, p = 2, width = 0.10))
  expect_equal(res$value[res$term == "necessary_N"], 357)
  expect_no_warning(ss_aipe_semipartial_r(true_r_sp = -0.8, p = 2, width = 0.10))
})

test_that("ss_aipe_semipartial_r() refuses the plan with no variance", {
  # At r_sp = 0 with R2_full = 1 the first-order variance is 0; the planner
  # returned necessary_N = p with an expected width of NaN.
  expect_error(ss_aipe_semipartial_r(0, p = 2, .2, true_R2_full = 1),
               "leaves no variance to plan for")
  near <- ss_aipe_semipartial_r(0, p = 2, .2, true_R2_full = .99)
  expect_true(is.finite(near$value[near$term == "expected_width"]))
})

test_that("ss_aipe_semipartial_r() plans at r_sp squared of 1/2 without a warning", {
  # The two worst-case numerators are equal there, and the plan is 100;
  # sqrt(.5)^2 > 1/2 in floating point.
  res <- expect_no_warning(ss_aipe_semipartial_r(sqrt(.5), p = 3, .2))
  expect_equal(res$value[res$term == "necessary_N"], 100)
  expect_no_warning(ss_aipe_semipartial_r(.75, p = 3, .2))
})

# ---- The count in its place, and the 1.0.x names --------------------------

test_that("ss_aipe_semipartial_r() takes p in its place, by position or by name", {
  named <- ss_aipe_semipartial_r(true_r_sp = .25, p = 4, width = .15)
  expect_identical(ss_aipe_semipartial_r(.25, 4, .15), named)
  expect_identical(do.call(ss_aipe_semipartial_r, list(.25, 4, .15)), named)
  z <- stats::qnorm(.975)
  for (r in c(.05, .25, .5)) for (k in c(1, 3, 6, 10)) {
    # N = p + ceiling(z^2 V / h^2), with the worst-case V when
    # true_R2_full is not given.
    expect_identical(
      ss_aipe_semipartial_r(true_r_sp = r, p = k + 1, width = .15)$value[1],
      k + 1 + ceiling(z^2 * (1 - r^2)^2 / .075^2))
  }
})

test_that("ss_aipe_semipartial_r() does not accept the DMAR 1.0.x names", {
  expect_error(ss_aipe_semipartial_r(r_sp = .25, p = 4, width = .15),
               "unused argument")
  expect_error(ss_aipe_semipartial_r(true_r_sp = .25, J = 3, width = .15),
               "unused argument")
  expect_error(ss_aipe_semipartial_r(true_r_sp = .25, p = 4, width = .15,
                                     R2_full = .4), "unused argument")
  expect_identical(names(formals(ss_aipe_semipartial_r)),
                   c("true_r_sp", "p", "width", "which_width", "conf_level",
                     "assurance", "true_R2_full", "G", "seed"))
})

test_that("ss_aipe_semipartial_r() gives no warning with assurance", {
  expect_no_warning(ss_aipe_semipartial_r(true_r_sp = .3, p = 2, width = .2,
                                          assurance = .8, G = 200, seed = 113))
})
