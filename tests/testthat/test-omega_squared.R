test_that("omega_squared() raw-argument interface returns the documented columns", {
  res <- omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("effect", "omega_squared",
                      "F_value", "df_effect", "df_error", "N"))
  expect_equal(res$effect, "overall")
})

test_that("omega_squared() reproduces the Hays (1994) closed-form: df_eff(F-1)/(df_eff(F-1)+N)", {
  F_value <- 11.221; df_effect <- 4; df_error <- 50; N <- 55
  num <- df_effect * (F_value - 1)
  expected <- num / (num + N)
  res <- omega_squared(F_value = F_value, df_effect = df_effect,
                       df_error = df_error, N = N)
  expect_equal(res$omega_squared, expected, tolerance = 1e-12)
})

test_that("omega_squared() truncates at zero when F < 1", {
  res <- omega_squared(F_value = 0.5, df_effect = 2, df_error = 50, N = 55)
  expect_equal(res$omega_squared, 0)
})

test_that("omega_squared() Bargman (1970) example matches Steiger (2004) Table 3 value", {
  res <- omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  # Steiger (2004) Table 3 reports omega-hat-squared = 0.426 for this F.
  expect_equal(round(res$omega_squared, 3), 0.426)
})

test_that("omega_squared() one-way model interface matches the raw interface", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  tbl <- summary(fit)[[1]]
  F_value   <- tbl[["F value"]][1]
  df_effect <- tbl[["Df"]][1]
  df_error  <- tbl[["Df"]][2]
  N         <- nobs(fit)
  via_model <- omega_squared(fit)
  via_raw   <- omega_squared(F_value = F_value, df_effect = df_effect,
                             df_error = df_error, N = N)
  expect_equal(via_model$omega_squared, via_raw$omega_squared,
               tolerance = 1e-12)
  expect_equal(via_model$F_value, via_raw$F_value, tolerance = 1e-12)
})

test_that("omega_squared() and ci_omega_squared() agree on the point estimate row-by-row", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  pt  <- omega_squared(fit)
  # The wool effect's low F legitimately clamps its lower limit to 0.
  expect_warning(ci <- ci_omega_squared(fit),
                 "below the alpha_lower critical value")
  expect_equal(pt$effect, ci$effect)
  expect_equal(pt$omega_squared, ci$omega_squared, tolerance = 1e-12)
  expect_equal(pt$F_value, ci$F_value, tolerance = 1e-12)
  expect_equal(pt$df_effect, ci$df_effect)
  expect_equal(pt$df_error,  ci$df_error)
  expect_equal(pt$N, ci$N)
})

test_that("omega_squared() factorial model returns one row per non-Residuals effect", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  res <- omega_squared(fit)
  expect_equal(nrow(res), 3L)
  expect_setequal(res$effect, c("wool", "tension", "wool:tension"))
  expect_true(all(res$omega_squared >= 0))
  expect_true(all(res$N == 54))
})

test_that("omega_squared() rejects bad inputs", {
  expect_error(omega_squared(F_value = 11.221, df_effect = 4, df_error = 50),
               "all of 'F_value'")
  expect_error(omega_squared(F_value = -1, df_effect = 4, df_error = 50, N = 55),
               "positive number")
  expect_error(omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 10),
               "greater than df_effect")
  expect_error(omega_squared(object = lm(mpg ~ wt, data = mtcars)) |> length(),
               NA)  # lm dispatch should work
  expect_error(omega_squared(object = list(not = "a model")),
               "aov or lm fit")
})
