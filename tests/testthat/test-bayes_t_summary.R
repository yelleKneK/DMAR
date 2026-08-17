# The Bayes factor depends on the data only through the t statistic and
# the sample sizes, so the summary-statistics form must reproduce the raw
# form exactly, not approximately.

test_that("summary statistics reproduce the raw-data result exactly", {
  set.seed(113)
  x <- rnorm(24, 0.6, 1.1)
  y <- rnorm(31, 0.0, 0.9)

  raw <- bayes_independent_t(x, y)
  summ <- bayes_independent_t(mean_1 = mean(x), sd_1 = sd(x), n_1 = length(x),
                              mean_2 = mean(y), sd_2 = sd(y), n_2 = length(y))
  expect_equal(raw, summ, tolerance = 1e-12)

  raw <- bayes_one_sample_t(x, mu_0 = 0.3)
  summ <- bayes_one_sample_t(mean = mean(x), sd = sd(x), n = length(x),
                             mu_0 = 0.3)
  expect_equal(raw, summ, tolerance = 1e-12)

  z <- rnorm(24, 0.2, 1)
  raw <- bayes_paired_t(x, z)
  d <- x - z
  summ <- bayes_paired_t(mean_diff = mean(d), sd_diff = sd(d), n = length(d))
  expect_equal(raw, summ, tolerance = 1e-12)
})

test_that("a standardized effect size enters through the summary form", {
  # d = 0.5 with n = 30 per group: mean difference of 0.5 in pooled SD
  # units, so mean_1 = 0.5, mean_2 = 0, both SDs 1.
  res <- bayes_independent_t(mean_1 = 0.5, sd_1 = 1, n_1 = 30,
                             mean_2 = 0.0, sd_2 = 1, n_2 = 30)
  expect_equal(res$value[res$term == "t"],
               0.5 * sqrt(30 * 30 / 60), tolerance = 1e-12)
})

test_that("the two entry paths are mutually exclusive and complete", {
  expect_error(bayes_one_sample_t(), "Supply either")
  expect_error(bayes_one_sample_t(x = rnorm(10), mean = 1), "not both")
  expect_error(bayes_one_sample_t(mean = 1, sd = 1), "'n' must be")
  expect_error(bayes_one_sample_t(mean = 1, sd = -1, n = 10),
               "'sd' must be positive")
  expect_error(bayes_independent_t(mean_1 = 1, sd_1 = 1, n_1 = 10),
               "'mean_2' must be")
  expect_error(bayes_paired_t(x = rnorm(5), y = rnorm(5), mean_diff = 1),
               "not both")
})
