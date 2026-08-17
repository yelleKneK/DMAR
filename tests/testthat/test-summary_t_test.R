test_that("summary_t_test() returns the documented rows", {
  res <- summary_t_test(mean_1 = 100, sd_1 = 15, n_1 = 30,
                        mean_2 = 108, sd_2 = 18, n_2 = 25)
  expect_setequal(res$term,
                  c("mean_difference", "t_statistic", "df",
                    "p_value", "lower_limit", "upper_limit"))
})

test_that("summary_t_test() matches stats::t.test() pooled-variance", {
  set.seed(113)
  x <- rnorm(20, 100, 15); y <- rnorm(20, 110, 15)
  ref <- stats::t.test(x, y, var.equal = TRUE)
  res <- summary_t_test(mean_1 = mean(x), sd_1 = sd(x), n_1 = length(x),
                        mean_2 = mean(y), sd_2 = sd(y), n_2 = length(y))
  expect_equal(res$value[res$term == "t_statistic"],
               as.numeric(ref$statistic), tolerance = 1e-10)
  expect_equal(res$value[res$term == "df"],
               as.numeric(ref$parameter), tolerance = 1e-10)
  expect_equal(res$value[res$term == "p_value"],
               ref$p.value, tolerance = 1e-10)
})

test_that("summary_t_test() Welch matches stats::t.test() Welch", {
  set.seed(113)
  x <- rnorm(20, 0, 1); y <- rnorm(20, 0.5, 2)
  ref <- stats::t.test(x, y, var.equal = FALSE)
  res <- summary_t_test(mean_1 = mean(x), sd_1 = sd(x), n_1 = length(x),
                        mean_2 = mean(y), sd_2 = sd(y), n_2 = length(y),
                        var_equal = FALSE)
  expect_equal(res$value[res$term == "t_statistic"],
               as.numeric(ref$statistic), tolerance = 1e-10)
  expect_equal(res$value[res$term == "df"],
               as.numeric(ref$parameter), tolerance = 1e-10)
})

test_that("summary_t_test() rejects invalid input", {
  expect_error(summary_t_test(0, -1, 30, 0, 1, 30), "non-negative")
  expect_error(summary_t_test(0, 1,  1, 0, 1, 30), "at least 2")
})

test_that("tidy()/glance() return broom columns matching the source table", {
  res <- summary_t_test(mean_1 = 100, sd_1 = 15, n_1 = 30,
                        mean_2 = 108, sd_2 = 18, n_2 = 25,
                        conf_level = 0.95)

  td <- generics::tidy(res)
  expect_equal(nrow(td), 1L)
  expect_named(td, c("term", "estimate", "ci_lower", "ci_upper",
                     "statistic", "df", "p_value", "conf_level"))
  expect_equal(td$term, "mean_difference")
  expect_equal(td$estimate,   res$value[res$term == "mean_difference"])
  expect_equal(td$ci_lower,   res$value[res$term == "lower_limit"])
  expect_equal(td$ci_upper,  res$value[res$term == "upper_limit"])
  expect_equal(td$statistic,  res$value[res$term == "t_statistic"])
  expect_equal(td$df,         res$value[res$term == "df"])
  expect_equal(td$p_value,    res$value[res$term == "p_value"])
  expect_equal(td$conf_level, 0.95)

  expect_equal(generics::glance(res), td)
})
