test_that("point estimate matches the closed-form (F-and-df) formula", {
  res <- eta_squared(F_value = 11.221, df_effect = 4, df_error = 50)
  expected <- 4 * 11.221 / (4 * 11.221 + 50)
  expect_equal(res$eta_squared, expected, tolerance = 1e-12)
  expect_equal(res$effect, "overall")
})

test_that("raw interface returns the documented columns", {
  res <- eta_squared(F_value = 11.221, df_effect = 4, df_error = 50)
  expect_named(res, c("effect", "eta_squared", "F_value", "df_effect", "df_error"))
  expect_equal(nrow(res), 1L)
})

test_that("point estimate equals SS_effect / SS_total in one-way (PlantGrowth)", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  tbl <- anova(fit)
  ss_eff <- tbl["group",     "Sum Sq"]
  ss_err <- tbl["Residuals", "Sum Sq"]
  expected_eta_sq <- ss_eff / (ss_eff + ss_err)   # one-way: SS_total = SS_eff + SS_err

  res <- eta_squared(fit)
  expect_equal(res$eta_squared, expected_eta_sq, tolerance = 1e-12)
  expect_equal(res$effect, "group")
})

test_that("model-based one-way matches raw interface with the same F/df", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  tbl <- anova(fit)
  F_val <- tbl["group", "F value"]
  df_e  <- tbl["group", "Df"]
  df_r  <- tbl["Residuals", "Df"]

  res_model <- eta_squared(fit)
  res_raw   <- eta_squared(F_value = F_val, df_effect = df_e, df_error = df_r)

  expect_equal(res_model$eta_squared, res_raw$eta_squared, tolerance = 1e-12)
})

test_that("factorial returns one row per non-Residuals effect", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  res <- eta_squared(fit)
  expect_equal(nrow(res), 3L)
  expect_setequal(res$effect, c("wool", "tension", "wool:tension"))
  expect_false("Residuals" %in% res$effect)
})

test_that("factorial: each row matches the partial eta squared formula", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  tbl <- anova(fit)
  df_err <- tbl["Residuals", "Df"]
  res <- eta_squared(fit)

  for (eff in res$effect) {
    F_val <- tbl[eff, "F value"]
    df_e  <- tbl[eff, "Df"]
    expected <- (df_e * F_val) / (df_e * F_val + df_err)   # partial
    expect_equal(res$eta_squared[res$effect == eff], expected, tolerance = 1e-12)
  }
})

test_that("supports an lm() fit as well as aov()", {
  fit_lm  <- lm(weight ~ group,  data = PlantGrowth)
  fit_aov <- aov(weight ~ group, data = PlantGrowth)
  expect_equal(eta_squared(fit_lm)$eta_squared,
               eta_squared(fit_aov)$eta_squared,
               tolerance = 1e-12)
})

test_that("F = 0 yields eta_squared = 0", {
  res <- eta_squared(F_value = 0, df_effect = 2, df_error = 30)
  expect_equal(res$eta_squared, 0)
})

test_that("eta squared is bounded in [0, 1]", {
  set.seed(113)
  Fs <- runif(50, 0.1, 100)
  df1 <- sample.int(10, 50, replace = TRUE)
  df2 <- sample(20:200, 50, replace = TRUE)
  for (i in seq_along(Fs)) {
    res <- eta_squared(F_value = Fs[i], df_effect = df1[i], df_error = df2[i])
    expect_gte(res$eta_squared, 0)
    expect_lte(res$eta_squared, 1)
  }
})

test_that("errors on missing or invalid arguments", {
  expect_error(eta_squared(), "Either provide 'object'")
  expect_error(eta_squared(F_value = 1.5, df_effect = 2), "Either provide")
  expect_error(eta_squared(F_value = -1, df_effect = 2, df_error = 10),
               "'F_value' must be")
  expect_error(eta_squared(F_value = 1, df_effect = 0, df_error = 10),
               "'df_effect' must be")
  expect_error(eta_squared(F_value = 1, df_effect = 2, df_error = -1),
               "'df_error' must be")
  expect_error(eta_squared(object = list(foo = 1)), "aov, lm, or aovlist")
})
