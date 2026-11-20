test_that("expected_R2() returns its value in the row expected_R2, as var_R2() does in var_R2", {
  result <- expected_R2(true_R2 = 0.3, N = 100, p = 5)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_identical(result$term, "expected_R2")
  expect_identical(var_R2(true_R2 = 0.3, N = 100, p = 5)$term, "var_R2")
  # The DMAR 1.0.x row, expected_value_population_R2, matches nothing.
  expect_length(result$value[result$term == "expected_value_population_R2"], 0L)
})

test_that("expected_R2() expected value is at least the population R2 (positive bias)", {
  pop <- 0.3
  result <- expected_R2(true_R2 = pop, N = 100, p = 5)
  expect_gte(result$value, pop)
})

test_that("expected_R2() bias shrinks toward zero as N grows", {
  pop <- 0.3
  bias_small <- expected_R2(true_R2 = pop, N = 50,    p = 5)$value - pop
  bias_large <- expected_R2(true_R2 = pop, N = 5000,  p = 5)$value - pop
  expect_gt(bias_small, bias_large)
  expect_lt(bias_large, 0.01)
})

test_that("expected_R2() takes true_R2; the DMAR 1.0.x name population_R2 stops", {
  expect_identical(names(formals(expected_R2)), c("true_R2", "N", "p"))
  expect_error(expected_R2(population_R2 = 0.3, N = 100, p = 5),
               "unused argument")
})

test_that("expected_R2() is the mean of the exact distribution of R^2", {
  # Oracle, with no package code: the negative binomial mixture of Beta
  # distributions described in test-var_R2.R.
  mixture_mean <- function(rho2, N, p) {
    m <- (N - 1) / 2
    b <- (N - p - 1) / 2
    j <- 0:20000
    w <- exp(lgamma(m + j) - lgamma(m) - lgamma(j + 1) + j * log(rho2) +
               m * log1p(-rho2))
    a <- p / 2 + j
    sum(w * a / (a + b)) / sum(w)
  }
  for (cell in list(c(0.3, 100, 5), c(0.05, 20, 2), c(0.6, 40, 6),
                    c(0.9, 500, 10), c(0.99, 50, 3))) {
    expect_equal(expected_R2(true_R2 = cell[1], N = cell[2], p = cell[3])$value,
                 mixture_mean(cell[1], cell[2], cell[3]), tolerance = 1e-12,
                 info = paste(cell, collapse = " "))
  }
})
