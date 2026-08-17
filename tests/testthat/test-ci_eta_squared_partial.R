test_that("partial CI point estimate matches the F-and-df formula", {
  res <- ci_eta_squared_partial(F_value = 11.221, df_effect = 4,
                                   df_error = 50, N = 55)
  expected <- 4 * 11.221 / (4 * 11.221 + 50)
  expect_equal(res$eta_squared_partial, expected, tolerance = 1e-12)
})

test_that("returns the documented columns (with partial value-name)", {
  res <- ci_eta_squared_partial(F_value = 11.221, df_effect = 4,
                                   df_error = 50, N = 55)
  expect_named(res, c("effect", "eta_squared_partial", "lower_limit",
                      "upper_limit", "F_value", "df_effect", "df_error", "N"))
})

test_that("partial and total CI bounds are numerically identical", {
  res_total   <- ci_eta_squared(F_value = 11.221, df_effect = 4,
                                   df_error = 50, N = 55)
  res_partial <- ci_eta_squared_partial(F_value = 11.221, df_effect = 4,
                                           df_error = 50, N = 55)
  expect_equal(res_partial$eta_squared_partial, res_total$eta_squared,
               tolerance = 1e-12)
  expect_equal(res_partial$lower_limit, res_total$lower_limit, tolerance = 1e-12)
  expect_equal(res_partial$upper_limit, res_total$upper_limit, tolerance = 1e-12)
})

test_that("model-based factorial returns one row per effect", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  # The wool effect's low F legitimately clamps its lower limit to 0.
  expect_warning(res <- ci_eta_squared_partial(fit),
                 "below the alpha_lower critical value")
  expect_equal(nrow(res), 3L)
  expect_setequal(res$effect, c("wool", "tension", "wool:tension"))
})

test_that("errors on missing or invalid arguments", {
  expect_error(ci_eta_squared_partial(), "Either provide 'object'")
  expect_error(ci_eta_squared_partial(F_value = -1, df_effect = 2,
                                         df_error = 10, N = 20),
               "'F_value' must be")
})

test_that("ci_eta_squared_partial() reports the lower-limit clamp once, in terms of partial eta squared, and still returns the interval", {
  msgs <- capture_warnings(
    result <- ci_eta_squared_partial(F_value = 1.2, df_effect = 4, df_error = 50, N = 55)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on partial eta squared is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_equal(result$lower_limit, 0)
  expect_gt(result$upper_limit, 0)
})
