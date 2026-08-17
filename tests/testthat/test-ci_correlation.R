## ci_r() and ci_R() -- confidence intervals for the population Pearson
## correlation and the population multiple correlation. The two functions
## share R/ci_correlation.R and one help page because their natural file
## names differ only by case, which case-insensitive filesystems and the
## R CMD check portable-file-names check both refuse.

test_that("ci_r() returns a tidy data.frame with lower_limit, r, upper_limit", {
  res <- ci_r(r = 0.8, n = 50, conf_level = 0.95)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_setequal(res$term, c("lower_limit", "r", "upper_limit"))
})

test_that("ci_r() interval brackets the point estimate", {
  res <- ci_r(r = 0.8, n = 50, conf_level = 0.95)
  lo <- res$value[res$term == "lower_limit"]
  hi <- res$value[res$term == "upper_limit"]
  expect_lt(lo, 0.8)
  expect_gt(hi, 0.8)
})

test_that("ci_r() narrows as n grows", {
  small <- ci_r(0.8, n = 30,   0.95)
  large <- ci_r(0.8, n = 5000, 0.95)
  w_small <- small$value[small$term == "upper_limit"] - small$value[small$term == "lower_limit"]
  w_large <- large$value[large$term == "upper_limit"] - large$value[large$term == "lower_limit"]
  expect_gt(w_small, w_large)
})

test_that("ci_r() widens as confidence rises", {
  ci90 <- ci_r(0.8, 100, 0.90)
  ci99 <- ci_r(0.8, 100, 0.99)
  w90 <- ci90$value[ci90$term == "upper_limit"] - ci90$value[ci90$term == "lower_limit"]
  w99 <- ci99$value[ci99$term == "upper_limit"] - ci99$value[ci99$term == "lower_limit"]
  expect_gt(w99, w90)
})

test_that("ci_r() stops for n < 4", {
  expect_error(ci_r(r = 0.35, n = 2), "at least 4")
  expect_error(ci_r(r = 0.35, n = 0), "at least 4")
  expect_error(ci_r(r = 0.35, n = -5), "at least 4")
})

test_that("ci_r() at n = 3 explains that the Fisher's Z interval is vacuous", {
  expect_error(ci_r(r = 0.35, n = 3), "vacuous")
  expect_error(ci_r(r = 0.35, n = 3), "1/\\(n - 3\\)")
  expect_error(ci_r(r = 0.35, n = 3), "\\[-1, 1\\]")
})

test_that("ci_r() rejects a non-scalar or non-numeric 'n'", {
  expect_error(ci_r(r = 0.35, n = c(10, 20)), "single")
  expect_error(ci_r(r = 0.35, n = "100"), "single")
  expect_error(ci_r(r = 0.35, n = NA), "single")
})

test_that("ci_r() is defined at the n = 4 floor", {
  res <- ci_r(r = 0.35, n = 4, conf_level = 0.95)
  lo <- res$value[res$term == "lower_limit"]
  hi <- res$value[res$term == "upper_limit"]
  expect_true(is.finite(lo))
  expect_true(is.finite(hi))
  expect_gt(lo, -1)
  expect_lt(hi, 1)
})

test_that("ci_R() returns the documented 3-row table", {
  result <- ci_R(R = 0.4, df_1 = 5, df_2 = 100, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("term", "value", "prob_less", "prob_greater") %in%
                    names(result)))
  expect_equal(result$term, c("lower_limit", "R", "upper_limit"))
})

test_that("ci_R() limits bracket the sample R on the correlation scale", {
  R <- 0.4
  result <- ci_R(R = R, df_1 = 5, df_2 = 100, conf_level = 0.95)
  ll <- result$value[result$term == "lower_limit"]
  ul <- result$value[result$term == "upper_limit"]
  expect_lte(ll, R)
  expect_gte(ul, R)
  expect_gte(ll, 0)
  expect_lte(ul, 1)
})

test_that("ci_R() limits are the square roots of the ci_R2() limits", {
  R <- 0.4
  r_ci  <- ci_R(R = R, df_1 = 5, df_2 = 100, conf_level = 0.95)
  r2_ci <- ci_R2(R2 = R^2, df_1 = 5, df_2 = 100, conf_level = 0.95)
  expect_equal(r_ci$value[r_ci$term == "lower_limit"],
               sqrt(r2_ci$value[r2_ci$term == "lower_limit"]))
  expect_equal(r_ci$value[r_ci$term == "upper_limit"],
               sqrt(r2_ci$value[r2_ci$term == "upper_limit"]))
})

test_that("ci_R() at one predictor agrees with |r| from ci_r() at the estimate", {
  # The multiple correlation from a one-predictor regression is the
  # absolute value of the Pearson correlation, so the two point-estimate
  # rows must coincide; the intervals differ by construction (Lee's
  # R^2 distribution versus Fisher's Z), so only the estimand is checked.
  r_obs <- 0.35
  res_R <- ci_R(R = r_obs, df_1 = 1, df_2 = 98, conf_level = 0.95)
  res_r <- ci_r(r = r_obs, n = 100, conf_level = 0.95)
  expect_equal(res_R$value[res_R$term == "R"],
               abs(res_r$value[res_r$term == "r"]))
})
