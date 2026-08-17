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
