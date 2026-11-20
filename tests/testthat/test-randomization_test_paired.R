test_that("randomization_test_paired() returns the documented rows", {
  res <- randomization_test_paired(c(95, 102, 98, 107, 105),
                                   c(102, 108, 100, 112, 109))
  expect_setequal(res$term,
                  c("statistic", "p_value", "n_pairs",
                    "n_evaluated", "exact"))
})

test_that("randomization_test_paired() exact enumeration matches 2^n", {
  res <- randomization_test_paired(c(1, 2, 3, 4, 5),
                                   c(2, 3, 4, 5, 6))
  expect_equal(res$value[res$term == "n_evaluated"], 2^5)
  expect_equal(res$value[res$term == "exact"], 1)
})

test_that("randomization_test_paired() finds significance when all signs are positive", {
  res <- randomization_test_paired(1:6, 11:16)
  expect_lt(res$value[res$term == "p_value"], 0.05)
})

test_that("randomization_test_paired() Monte-Carlo branch fires for large n", {
  set.seed(113)
  x <- rnorm(25); y <- rnorm(25)
  res <- randomization_test_paired(x, y, exact = FALSE, n_resamples = 500L)
  expect_equal(res$value[res$term == "exact"], 0)
  expect_equal(res$value[res$term == "n_evaluated"], 500)
})

test_that("randomization_test_paired() flips p-value direction for one-sided alternatives", {
  set.seed(113)
  d <- rnorm(15, 1, 0.5)
  x <- runif(15); y <- x + d
  p_less    <- randomization_test_paired(x, y, alternative = "less")$value[2]
  p_greater <- randomization_test_paired(x, y, alternative = "greater")$value[2]
  expect_gt(p_less,    0.5)  # y > x, so "less" is anti-direction
  expect_lt(p_greater, 0.5)
})

test_that("randomization_test_paired() rejects mismatched lengths", {
  expect_error(randomization_test_paired(1:5, 1:6), "same length")
})


# ---- Ties, boundary values, and the count of resamples ------------------

rtp_value <- function(x, term) x$value[x$term == term]

# Exact p-value for whole-number scores: every signed sum of the
# differences is a whole number, so ties are exact.
whole_number_sign_flip_p <- function(k_x, k_y, alternative) {
  d <- k_y - k_x
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(d))))
  s <- as.vector(signs %*% d)
  obs <- sum(d)
  switch(alternative,
         two_sided = mean(abs(s) >= abs(obs)),
         greater   = mean(s >= obs),
         less      = mean(s <= obs))
}


test_that("randomization_test_paired() gives constant differences the smallest p-value under the mean and refuses the paired t", {
  # Differences that are all the same nonzero value, as stored (whole
  # numbers) or only in exact arithmetic (0.3 - 0.1 and 0.9 - 0.7 differ in
  # their last binary place). The paired t has a zero standard error there
  # and no value; through round 1 of 1.1.0 development it was reported as
  # Inf, and in 1.0.0 as 0 with p = 1.
  tenths <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  designs <- list(list(1:5, 1:5 + 2), list(tenths, tenths + 0.2))
  p_exact <- c(two_sided = 2 / 32, greater = 1 / 32, less = 1)
  for (d in designs) for (alternative in names(p_exact)) {
    res <- randomization_test_paired(d[[1]], d[[2]], alternative = alternative)
    expect_identical(rtp_value(res, "p_value"), p_exact[[alternative]])
    e <- expect_error(
      randomization_test_paired(d[[1]], d[[2]], statistic = "t",
                                alternative = alternative),
      class = "dmar_studentized_statistic_undefined")
    expect_s3_class(e, "dmar_refusal")
    expect_identical(e$reason, "constant_differences")
    expect_identical(e$alternative, alternative)
    expect_identical(e$n, 5L)
    expect_identical(e$p_value_mean, p_exact[[alternative]])
    expect_match(conditionMessage(e), "zero denominator", fixed = TRUE)
    expect_match(conditionMessage(e), "Use statistic = \"mean\".",
                 fixed = TRUE)
    # Negative differences: the two one-sided values trade places.
    mirror <- c(two_sided = "two_sided", greater = "less",
                less = "greater")[[alternative]]
    e <- expect_error(
      randomization_test_paired(d[[2]], d[[1]], statistic = "t",
                                alternative = mirror),
      class = "dmar_studentized_statistic_undefined")
    expect_identical(e$p_value_mean, p_exact[[alternative]])
  }
  # The message states the value the mean statistic returns, and says why
  # a p-value of 1 is 1.
  e <- expect_error(randomization_test_paired(1:5, 1:5 + 2, statistic = "t"),
                    class = "dmar_studentized_statistic_undefined")
  expect_match(conditionMessage(e), "All 5 paired differences equal 2,",
               fixed = TRUE)
  expect_match(conditionMessage(e), "of the 32 sign patterns", fixed = TRUE)
  expect_match(conditionMessage(e),
               "the exact p-value 2/32 = 0.0625 (two_sided)", fixed = TRUE)
  e <- expect_error(randomization_test_paired(1:5, 1:5 + 2, statistic = "t",
                                              alternative = "less"),
                    class = "dmar_studentized_statistic_undefined")
  expect_match(conditionMessage(e),
               "the exact p-value 1 (less), because every sign pattern",
               fixed = TRUE)
  # Six equal negative differences.
  e <- expect_error(randomization_test_paired(rep(3, 6), rep(1, 6),
                                              statistic = "t",
                                              alternative = "less"),
                    class = "dmar_studentized_statistic_undefined")
  expect_identical(e$p_value_mean, 1 / 64)
})

test_that("differences equal only to within rounding are described as such", {
  # The stored differences are 1, 1.25, and 1.5, and the scores are 1e15
  # times their spread, so the comparison cannot tell their standard
  # deviation from 0. statistic = "mean" on the stored values gives .75,
  # not the 2/8 of differences that are all equal.
  x <- 1e15 + c(0, 0, 0)
  y <- 1e15 + c(1, 1.25, 1.5)
  e <- expect_error(randomization_test_paired(x, y, statistic = "t"),
                    class = "dmar_studentized_statistic_undefined")
  msg <- conditionMessage(e)
  expect_match(msg, "to within the rounding of scores of this size (they span 1 to 1.5)",
               fixed = TRUE)
  expect_match(msg, "cannot be told from 0", fixed = TRUE)
  expect_match(msg, "The p-value of differences that are all equal is 2/8 = 0.25",
               fixed = TRUE)
  expect_no_match(msg, "the exact p-value", fixed = TRUE)
  expect_no_match(msg, "so their standard deviation is 0", fixed = TRUE)
  res <- randomization_test_paired(x, y, statistic = "mean")
  expect_equal(res$value[res$term == "p_value"], 0.75)
  # Differences that are exactly equal keep the earlier wording.
  e <- expect_error(randomization_test_paired(1:5, 1:5 + 2, statistic = "t"),
                    class = "dmar_studentized_statistic_undefined")
  expect_match(conditionMessage(e), "so their standard deviation is 0",
               fixed = TRUE)
  expect_match(conditionMessage(e), "with the exact p-value", fixed = TRUE)
})


test_that("randomization_test_paired() refuses the paired t in the Monte Carlo branch without drawing", {
  set.seed(8128)
  before <- .Random.seed
  e <- expect_error(
    randomization_test_paired(rep(1, 30), rep(3, 30), statistic = "t"),
    class = "dmar_studentized_statistic_undefined")
  expect_identical(.Random.seed, before)
  expect_identical(e$p_value_mean, 2 / 2^30)
  expect_match(conditionMessage(e),
               paste0("With 30 pairs the test is run by Monte Carlo sampling, ",
                      "which estimates this value and cannot report less ",
                      "than 1/(n_resamples + 1) = 0.00009999."),
               fixed = TRUE)
  # Sampling asked for where the default would enumerate.
  e <- expect_error(
    randomization_test_paired(1:5, 1:5 + 2, statistic = "t", exact = FALSE,
                              n_resamples = 999, seed = 4),
    class = "dmar_studentized_statistic_undefined")
  expect_identical(.Random.seed, before)
  expect_match(conditionMessage(e),
               "With exact = FALSE the test is run by Monte Carlo sampling",
               fixed = TRUE)
  expect_match(conditionMessage(e), "1/(n_resamples + 1) = 0.001.",
               fixed = TRUE)
  # The mean statistic runs, and reports the smallest value the run can.
  res <- randomization_test_paired(rep(1, 30), rep(3, 30), seed = 4)
  expect_identical(rtp_value(res, "statistic"), 2)
  expect_equal(rtp_value(res, "p_value"), 1 / 10001)
})


test_that("randomization_test_paired() keeps a finite paired t when the differences are equal in size and mixed in sign", {
  # The pattern with every sign positive (and its mirror image) has a zero
  # standard error, but it is not the observed one: the observed paired t
  # is defined and the count is the mean difference's.
  y <- c(2, 2, 2, -2, 2)
  res_t <- randomization_test_paired(rep(0, 5), y, statistic = "t")
  expect_equal(rtp_value(res_t, "statistic"), mean(y) * sqrt(5) / sd(y))
  expect_equal(rtp_value(res_t, "statistic"), 1.5)
  expect_identical(rtp_value(res_t, "p_value"), 12 / 32)
  res_m <- randomization_test_paired(rep(0, 5), y)
  expect_identical(rtp_value(res_m, "p_value"), 12 / 32)
})


test_that("randomization_test_paired() returns 0 and p = 1 when every difference is zero", {
  for (statistic in c("mean", "t")) {
    res <- randomization_test_paired(c(4, 7, 9, 2), c(4, 7, 9, 2),
                                     statistic = statistic)
    expect_identical(rtp_value(res, "statistic"), 0)
    expect_identical(rtp_value(res, "p_value"), 1)
  }
})


test_that("randomization_test_paired() counts ties, and only ties, under either statistic", {
  x <- rep(0, 6)
  for (statistic in c("mean", "t")) {
    p_of <- function(y, alternative) {
      rtp_value(randomization_test_paired(x, y, statistic = statistic,
                                          alternative = alternative),
                "p_value")
    }
    # Five differences of 1 and one of -1. The six patterns with one
    # negative sign tie with the observed one: 7 of 64 sums are at least
    # the observed 4, and 14 are at least as large in absolute value.
    y <- c(1, 1, 1, 1, 1, -1)
    expect_equal(p_of(y, "two_sided"), 14 / 64, tolerance = 1e-12)
    expect_equal(p_of(y, "greater"), 7 / 64, tolerance = 1e-12)

    # One difference larger by 1e-12, which is a difference in the data
    # and far above the rounding of sums of this size. The pattern that
    # puts the negative sign on that difference now falls 2e-12 below the
    # observed sum, so 6 of 64 sums are at least the observed one.
    y <- c(1, 1, 1, 1, 1 + 1e-12, -1)
    expect_equal(p_of(y, "two_sided"), 12 / 64, tolerance = 1e-12)
    expect_equal(p_of(y, "greater"), 6 / 64, tolerance = 1e-12)
  }
})


test_that("randomization_test_paired() counts ties among thousandths on an offset of a million", {
  x <- c(1000000.002, 1000000.005, 1000000.000, 1000000.002, 1000000.001,
         1000000.001)
  y <- c(1000000.003, 1000000.007, 1000000.002, 1000000.001,  999999.999,
         1000000.001)
  p_offset <- rtp_value(randomization_test_paired(x, y), "p_value")
  expect_equal(p_offset, 0.8125, tolerance = 1e-12)
  # The same scores without the offset, entered as thousandths.
  x_0 <- c(0.002, 0.005, 0.000, 0.002, 0.001, 0.001)
  y_0 <- c(0.003, 0.007, 0.002, 0.001, -0.001, 0.001)
  expect_equal(rtp_value(randomization_test_paired(x_0, y_0), "p_value"),
               p_offset, tolerance = 1e-12)
  expect_equal(p_offset,
               whole_number_sign_flip_p(c(2, 5, 0, 2, 1, 1),
                                        c(3, 7, 2, 1, -1, 1), "two_sided"),
               tolerance = 1e-12)
})


test_that("randomization_test_paired() gives one p-value for both statistics", {
  set.seed(20261109)
  for (i in 1:50) {
    n <- sample(4:10, 1)
    x <- rnorm(n, 50, 10)
    y <- x + rt(n, 3)
    alternative <- c("two_sided", "less", "greater")[i %% 3 + 1]
    exact_mean <- randomization_test_paired(x, y, alternative = alternative)
    exact_t <- randomization_test_paired(x, y, statistic = "t",
                                         alternative = alternative)
    expect_identical(rtp_value(exact_t, "p_value"),
                     rtp_value(exact_mean, "p_value"))
    mc_mean <- randomization_test_paired(x, y, alternative = alternative,
                                         exact = FALSE, n_resamples = 300,
                                         seed = i)
    mc_t <- randomization_test_paired(x, y, statistic = "t",
                                      alternative = alternative,
                                      exact = FALSE, n_resamples = 300,
                                      seed = i)
    expect_identical(rtp_value(mc_t, "p_value"), rtp_value(mc_mean, "p_value"))
  }
})


test_that("randomization_test_paired() Monte Carlo draws and counts as DMAR 1.1.0 did for the mean statistic", {
  # The p-value recorded from the reviewed 1.1.0 build (df120a0d) for this
  # seeded call. It pins the order of the draws.
  set.seed(113)
  x <- rnorm(30, 100, 15)
  y <- x + rnorm(30, 3, 12)
  res <- randomization_test_paired(x, y, exact = FALSE, n_resamples = 2000,
                                   seed = 113)
  expect_equal(rtp_value(res, "p_value"), 0.025487256371814093,
               tolerance = 1e-14)
  expect_equal(rtp_value(res, "statistic"), mean(y - x), tolerance = 1e-14)
})


test_that("randomization_test_paired() equals whole-number arithmetic in four unit systems", {
  set.seed(20261109)
  for (i in 1:20) {
    n <- sample(4:8, 1)
    k_x <- sample(0:9, n, replace = TRUE)
    k_y <- sample(0:9, n, replace = TRUE)
    k_x_outlier <- k_x
    k_x_outlier[1] <- 3e9
    for (alternative in c("two_sided", "less", "greater")) {
      exact_p <- whole_number_sign_flip_p(k_x, k_y, alternative)
      versions <- list(
        tenths            = list(k_x / 10, k_y / 10),
        tenths_offset     = list((10003 + k_x) / 10, (10003 + k_y) / 10),
        thousandths_at_1e6 = list((1e9 + k_x) / 1000, (1e9 + k_y) / 1000)
      )
      for (v in versions) {
        expect_equal(rtp_value(randomization_test_paired(
          v[[1]], v[[2]], alternative = alternative), "p_value"),
          exact_p, tolerance = 1e-12)
      }
      # One score about 1e9 times the rest.
      expect_equal(rtp_value(randomization_test_paired(
        k_x_outlier, k_y, alternative = alternative), "p_value"),
        whole_number_sign_flip_p(k_x_outlier, k_y, alternative),
        tolerance = 1e-12)
    }
  }
})


test_that("randomization_test_paired() requires n_resamples to be a whole number of at least 1", {
  set.seed(113)
  x <- rnorm(12)
  y <- x + rnorm(12, 0.5)
  for (bad in list(2.5, 0, -1, NA_real_, Inf, "10", c(10, 20))) {
    expect_error(randomization_test_paired(x, y, exact = FALSE,
                                           n_resamples = bad),
                 "'n_resamples'", class = "dmar_count_invalid")
  }
  a <- randomization_test_paired(x, y, exact = FALSE, n_resamples = 200,
                                 seed = 113)
  b <- randomization_test_paired(x, y, exact = FALSE, n_resamples = 200L,
                                 seed = 113)
  expect_identical(a$value, b$value)
  # The denominator of the Monte Carlo p-value is the number of draws made.
  expect_identical(rtp_value(a, "n_evaluated"), 200)
  expect_equal(201 * rtp_value(a, "p_value"),
               round(201 * rtp_value(a, "p_value")), tolerance = 1e-10)
})


test_that("randomization_test_paired() refuses scores that are not finite", {
  expect_error(randomization_test_paired(c(1, 2, 3, Inf), c(2, 3, 4, 5)),
               "must be finite")
  expect_error(randomization_test_paired(c(1, 2, 3, 4), c(2, -Inf, 4, 5)),
               "must be finite")
})

test_that("randomization_test_paired() reads scores stored as integers as numbers", {
  # The rounding bounds sum absolute scores, which overflows in integer
  # arithmetic: these calls warned ("NAs produced by integer overflow") and
  # the first stopped.
  pop_before <- c(1500000000L, 1200000000L, 1100000000L, 1300000000L,
                  900000000L)
  pop_after  <- c(1400000000L, 1250000000L, 1000000000L, 1350000000L,
                  950000000L)
  for (statistic in c("t", "mean")) {
    as_int <- expect_no_warning(
      randomization_test_paired(pop_before, pop_after, statistic = statistic))
    as_dbl <- randomization_test_paired(as.numeric(pop_before),
                                        as.numeric(pop_after),
                                        statistic = statistic)
    expect_equal(as_int$value, as_dbl$value)
  }
  expect_equal(as_int$value[as_int$term == "p_value"], 1)

  set.seed(41)
  price_before <- as.integer(round(rnorm(500, 300000, 60000), -2))
  price_after  <- price_before + as.integer(round(rnorm(500, 4000, 20000), -2))
  as_int <- expect_no_warning(
    randomization_test_paired(price_before, price_after, seed = 2))
  as_dbl <- randomization_test_paired(as.numeric(price_before),
                                      as.numeric(price_after), seed = 2)
  expect_equal(as_int$value, as_dbl$value)
})

test_that("randomization_test_paired() checks n_resamples under exact enumeration too", {
  x <- c(5, 7, 6, 9, 8, 7, 10, 6)
  y <- c(6, 9, 6, 11, 9, 9, 12, 7)
  expect_error(randomization_test_paired(x, y, n_resamples = NULL),
               class = "dmar_count_invalid")
})

test_that("randomization_test_paired() has the attainable size for continuous differences and no more with ties", {
  # Under the null each of the 2^n sign patterns of the absolute
  # differences is equally likely to be the one observed, so the exact
  # size of the test at .05, given the absolute differences, is the share
  # of patterns whose p-value is .05 or less. When no two patterns give
  # the same |S| (continuous differences, almost surely) the two-sided
  # p-values are k / 2^(n - 1), k = 1, ..., 2^(n - 1), each taken by two
  # patterns, so the size is floor(.05 * 2^(n - 1)) / 2^(n - 1) = 6/128 for
  # n = 8, whatever the absolute differences (light-tailed, heavy-tailed,
  # or one outlier); one-sided it is floor(.05 * 2^n) / 2^n = 12/256.
  # With differences of equal size many patterns give the same |S| and the
  # size is lower (the page: "does not exceed the nominal level").
  n <- 8L
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  size_given <- function(abs_d, alternative) {
    p <- apply(signs, 1L, function(s) {
      tab <- randomization_test_paired(rep(0, n), s * abs_d,
                                       alternative = alternative)
      tab$value[tab$term == "p_value"]
    })
    mean(p <= .05)
  }
  set.seed(113)
  for (abs_d in list(abs(rnorm(n)), abs(rt(n, df = 1)),
                     c(abs(rnorm(n - 1)), 40))) {
    expect_identical(size_given(abs_d, "two_sided"), 6 / 128)
    expect_identical(size_given(abs_d, "greater"), 12 / 256)
  }
  ties <- size_given(rep(1, n), "two_sided")
  expect_lt(ties, 6 / 128)
  # Equal differences: |S| = |2k - n| with k positive signs; p <= .05 only
  # for all signs alike, 2 patterns of 256.
  expect_identical(ties, 2 / 256)
})
