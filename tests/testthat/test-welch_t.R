test_that("welch_t() returns the documented rows", {
  set.seed(113)
  x <- rnorm(20, 0, 1); y <- rnorm(20, 0.5, 2)
  res <- welch_t(x, y)
  expect_setequal(res$term,
                  c("mean_difference", "t_statistic", "df", "p_value",
                    "lower_limit", "upper_limit",
                    "mean_x", "mean_y", "sd_x", "sd_y", "n_x", "n_y"))
})

test_that("welch_t() matches stats::t.test(var.equal = FALSE)", {
  set.seed(113)
  x <- rnorm(15, 0, 1); y <- rnorm(20, 0.5, 2)
  ref <- stats::t.test(x, y, var.equal = FALSE)
  res <- welch_t(x, y)
  expect_equal(res$value[res$term == "t_statistic"],
               as.numeric(ref$statistic), tolerance = 1e-10)
  expect_equal(res$value[res$term == "df"],
               as.numeric(ref$parameter), tolerance = 1e-10)
  expect_equal(res$value[res$term == "p_value"],
               ref$p.value, tolerance = 1e-10)
  expect_equal(res$value[res$term == "lower_limit"],
               ref$conf.int[1], tolerance = 1e-10)
})

test_that("welch_t() one-sided alternatives flip the p-value", {
  set.seed(113)
  x <- rnorm(20, 0); y <- rnorm(20, 1)
  p_less    <- welch_t(x, y, alternative = "less")$value[
                  welch_t(x, y, alternative = "less")$term == "p_value"]
  p_greater <- welch_t(x, y, alternative = "greater")$value[
                  welch_t(x, y, alternative = "greater")$term == "p_value"]
  expect_lt(p_less, 0.5)
  expect_gt(p_greater, 0.5)
  expect_equal(p_less + p_greater, 1, tolerance = 1e-10)
})

test_that("welch_t() rejects degenerate input", {
  expect_error(welch_t(1, 1:5),       "at least 2")
  expect_error(welch_t(rep(1, 5), rep(2, 5)), "zero within-group variance")
})

test_that("tidy()/glance() return broom columns matching the source table", {
  set.seed(113)
  x <- rnorm(20, mean = 100, sd = 15)
  y <- rnorm(20, mean = 110, sd = 25)
  res <- welch_t(x, y, conf_level = 0.95)

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

  # glance() coincides with tidy() for a single-estimand test.
  expect_equal(generics::glance(res), td)
})
