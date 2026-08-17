test_that("omega_squared_partial() raw-argument interface returns the documented columns", {
  res <- omega_squared_partial(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("effect", "omega_squared_partial",
                      "F_value", "df_effect", "df_error", "N"))
  expect_equal(res$effect, "overall")
})

test_that("omega_squared_partial() reproduces df_eff(F-1)/(df_eff(F-1)+N)", {
  F_value <- 11.221; df_effect <- 4; df_error <- 50; N <- 55
  num <- df_effect * (F_value - 1)
  expected <- num / (num + N)
  res <- omega_squared_partial(F_value = F_value, df_effect = df_effect,
                               df_error = df_error, N = N)
  expect_equal(res$omega_squared_partial, expected, tolerance = 1e-12)
})

test_that("omega_squared_partial() truncates at zero when F < 1", {
  res <- omega_squared_partial(F_value = 0.5, df_effect = 2, df_error = 50, N = 55)
  expect_equal(res$omega_squared_partial, 0)
})

test_that("omega_squared_partial() Bargman (1970) one-way matches Steiger (2004) Table 3", {
  res <- omega_squared_partial(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  # Steiger (2004) Table 3: omega-hat-squared = 0.426. In a one-way design,
  # partial omega squared equals total omega squared.
  expect_equal(round(res$omega_squared_partial, 3), 0.426)
})

test_that("omega_squared_partial() one-way model interface matches the raw interface", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  tbl <- summary(fit)[[1]]
  F_value   <- tbl[["F value"]][1]
  df_effect <- tbl[["Df"]][1]
  df_error  <- tbl[["Df"]][2]
  N         <- nobs(fit)
  via_model <- omega_squared_partial(fit)
  via_raw   <- omega_squared_partial(F_value = F_value, df_effect = df_effect,
                                     df_error = df_error, N = N)
  expect_equal(via_model$omega_squared_partial, via_raw$omega_squared_partial,
               tolerance = 1e-12)
  expect_equal(via_model$F_value, via_raw$F_value, tolerance = 1e-12)
})

test_that("omega_squared_partial() agrees row-by-row with omega_squared() (named pair)", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  pt_partial <- omega_squared_partial(fit)
  pt_omega   <- omega_squared(fit)
  expect_equal(pt_partial$effect, pt_omega$effect)
  expect_equal(pt_partial$omega_squared_partial, pt_omega$omega_squared,
               tolerance = 1e-12)
  expect_equal(pt_partial$F_value, pt_omega$F_value, tolerance = 1e-12)
})

test_that("omega_squared_partial() agrees row-by-row with ci_omega_squared() point estimate", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  pt  <- omega_squared_partial(fit)
  # The wool effect's low F legitimately clamps its lower limit to 0.
  expect_warning(ci <- ci_omega_squared(fit),
                 "below the alpha_lower critical value")
  expect_equal(pt$effect, ci$effect)
  expect_equal(pt$omega_squared_partial, ci$omega_squared, tolerance = 1e-12)
  expect_equal(pt$F_value, ci$F_value, tolerance = 1e-12)
  expect_equal(pt$df_effect, ci$df_effect)
  expect_equal(pt$df_error,  ci$df_error)
  expect_equal(pt$N, ci$N)
})

test_that("omega_squared_partial() factorial model returns one row per non-Residuals effect", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  res <- omega_squared_partial(fit)
  expect_equal(nrow(res), 3L)
  expect_setequal(res$effect, c("wool", "tension", "wool:tension"))
  expect_true(all(res$omega_squared_partial >= 0))
  expect_true(all(res$N == 54))
})

test_that("omega_squared_partial() rejects bad inputs", {
  expect_error(omega_squared_partial(F_value = 11.221, df_effect = 4, df_error = 50),
               "all of 'F_value'")
  expect_error(omega_squared_partial(F_value = -1, df_effect = 4, df_error = 50, N = 55),
               "positive number")
  expect_error(omega_squared_partial(F_value = 11.221, df_effect = 4, df_error = 50, N = 10),
               "greater than df_effect")
  expect_error(omega_squared_partial(object = lm(mpg ~ wt, data = mtcars)) |> length(),
               NA)  # lm dispatch should work
  expect_error(omega_squared_partial(object = list(not = "a model")),
               "aov, lm, or aovlist fit")
})
