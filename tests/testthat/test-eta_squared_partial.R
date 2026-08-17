test_that("partial point estimate equals the F-and-df formula", {
  res <- eta_squared_partial(F_value = 11.221, df_effect = 4, df_error = 50)
  expected <- 4 * 11.221 / (4 * 11.221 + 50)
  expect_equal(res$eta_squared_partial, expected, tolerance = 1e-12)
})

test_that("returns the documented columns (with the partial value-name)", {
  res <- eta_squared_partial(F_value = 11.221, df_effect = 4, df_error = 50)
  expect_named(res, c("effect", "eta_squared_partial", "F_value",
                      "df_effect", "df_error"))
})

test_that("partial and total agree in one-way (PlantGrowth)", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  res_total   <- eta_squared(fit)
  res_partial <- eta_squared_partial(fit)
  expect_equal(res_total$eta_squared, res_partial$eta_squared_partial,
               tolerance = 1e-12)
})

test_that("partial and total return identical numeric values in factorial", {
  # By construction the two functions share the F-and-df formula.
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  res_total   <- eta_squared(fit)
  res_partial <- eta_squared_partial(fit)
  expect_equal(res_total$eta_squared, res_partial$eta_squared_partial,
               tolerance = 1e-12)
  expect_equal(res_total$effect, res_partial$effect)
})

test_that("errors on missing or invalid arguments", {
  expect_error(eta_squared_partial(), "Either provide 'object'")
  expect_error(eta_squared_partial(F_value = -1, df_effect = 2,
                                      df_error = 10),
               "'F_value' must be")
})
