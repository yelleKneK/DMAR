test_that("point estimate matches the closed-form formula for Bargman (1970)", {
  res <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  expected <- 4 * (11.221 - 1) / (4 * (11.221 - 1) + 55)
  expect_equal(res$omega_squared, expected, tolerance = 1e-6)
  expect_equal(res$effect, "overall")
})

test_that("raw-interface CI agrees with the direct Steiger transformation of ci_pvaf", {
  res <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  pvaf <- ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
  # Both apply ncp / (ncp + N); the bounds should match to high precision.
  expect_equal(res$lower_limit,
               pvaf$value[pvaf$term == "lower_limit"], tolerance = 1e-6)
  expect_equal(res$upper_limit,
               pvaf$value[pvaf$term == "upper_limit"], tolerance = 1e-6)
})

test_that("CI brackets the point estimate", {
  res <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  expect_gte(res$omega_squared, res$lower_limit)
  expect_lte(res$omega_squared, res$upper_limit)
})

test_that("conf_level narrows with smaller coverage and widens with larger", {
  r90 <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
                          conf_level = 0.90)
  r99 <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
                          conf_level = 0.99)
  expect_gt(r99$upper_limit - r99$lower_limit,
            r90$upper_limit - r90$lower_limit)
})

test_that("low F (CI includes 0) yields lower_limit == 0", {
  # Deliberately non-significant: F close to 1
  expect_warning(
    res <- ci_omega_squared(F_value = 0.5, df_effect = 2, df_error = 30, N = 33),
    "below the alpha_lower critical value")
  expect_equal(res$lower_limit, 0)
  expect_equal(res$omega_squared, 0)  # truncated at 0 as F < 1
})

test_that("model-based one-way matches raw interface with the same F/df/N", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  tbl <- anova(fit)
  F_val <- tbl["group", "F value"]
  df_e  <- tbl["group", "Df"]
  df_r  <- tbl["Residuals", "Df"]
  N_val <- nobs(fit)

  res_model <- ci_omega_squared(fit)
  res_raw   <- ci_omega_squared(F_value = F_val, df_effect = df_e,
                                df_error = df_r, N = N_val)

  expect_equal(res_model$omega_squared, res_raw$omega_squared)
  expect_equal(res_model$lower_limit,   res_raw$lower_limit)
  expect_equal(res_model$upper_limit,   res_raw$upper_limit)
  expect_equal(res_model$effect, "group")
})

test_that("model-based factorial returns one row per effect (not Residuals)", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  # The wool effect's low F legitimately clamps its lower limit to 0.
  expect_warning(res <- ci_omega_squared(fit),
                 "below the alpha_lower critical value")
  expect_equal(nrow(res), 3L)
  expect_setequal(res$effect, c("wool", "tension", "wool:tension"))
  expect_false("Residuals" %in% res$effect)
})

test_that("model-based factorial matches the closed-form point estimates", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  tbl <- anova(fit)
  N_val <- nobs(fit)
  # The wool effect's low F legitimately clamps its lower limit to 0.
  expect_warning(res <- ci_omega_squared(fit),
                 "below the alpha_lower critical value")

  for (eff in res$effect) {
    F_val <- tbl[eff, "F value"]
    df_e  <- tbl[eff, "Df"]
    expected <- df_e * (F_val - 1) / (df_e * (F_val - 1) + N_val)
    expected <- max(0, expected)
    expect_equal(res$omega_squared[res$effect == eff], expected, tolerance = 1e-9)
  }
})

test_that("supports an lm() fit as well as aov()", {
  fit_lm  <- lm(weight ~ group, data = PlantGrowth)
  fit_aov <- aov(weight ~ group, data = PlantGrowth)
  res_lm  <- ci_omega_squared(fit_lm)
  res_aov <- ci_omega_squared(fit_aov)
  expect_equal(res_lm$omega_squared, res_aov$omega_squared)
  expect_equal(res_lm$lower_limit,   res_aov$lower_limit)
  expect_equal(res_lm$upper_limit,   res_aov$upper_limit)
})

test_that("asymmetric alpha arguments are respected", {
  res_sym  <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                               N = 55, conf_level = 0.90)
  res_asym <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                               N = 55, alpha_lower = 0, alpha_upper = 0.10)
  # With alpha_lower = 0 the lower bound on the NCP, and thus omega^2, is 0
  expect_equal(res_asym$lower_limit, 0)
  # Upper limit of asym (alpha_upper = .10) should match upper of sym (alpha_upper = .05)... no.
  # sym has alpha_upper = .05, asym has alpha_upper = .10, so asym_upper < sym_upper.
  expect_lt(res_asym$upper_limit, res_sym$upper_limit)
})

test_that("errors on missing or invalid arguments", {
  expect_error(ci_omega_squared(), "Either provide 'object'")
  expect_error(ci_omega_squared(F_value = 1.5, df_effect = 2), "Either provide")
  expect_error(ci_omega_squared(F_value = -1, df_effect = 2,
                                df_error = 10, N = 20),
               "'F_value' must be")
  expect_error(ci_omega_squared(F_value = 3, df_effect = 2,
                                df_error = 10, N = 5),
               "greater than")
  expect_error(ci_omega_squared(object = list(foo = 1)), "aov or lm")
})

test_that("returned data.frame has the documented columns", {
  res <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  expect_named(res, c("effect", "omega_squared", "lower_limit", "upper_limit",
                      "F_value", "df_effect", "df_error", "N"))
})

test_that("ci_omega_squared() reports the lower-limit clamp once, in terms of omega squared, and still returns the interval", {
  msgs <- capture_warnings(
    result <- ci_omega_squared(F_value = 1.2, df_effect = 4, df_error = 50, N = 55)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on omega squared is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_equal(result$lower_limit, 0)
  expect_gt(result$upper_limit, 0)
})

test_that("ci_omega_squared() on a model with several clamped effects warns once, not once per effect", {
  set.seed(113)
  d <- data.frame(
    y = rnorm(60),
    a = factor(rep(1:3, each = 20)),
    b = factor(rep(1:2, 30))
  )
  fit <- aov(y ~ a + b, data = d)
  msgs <- capture_warnings(result <- ci_omega_squared(fit))
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on omega squared is 0")
  expect_equal(result$lower_limit, c(0, 0))
})

test_that("ci_omega_squared() names itself and the inputs when the noncentral root search fails", {
  expect_error(
    ci_omega_squared(F_value = Inf, df_effect = 4, df_error = 50, N = 55),
    "In ci_omega_squared\\(\\).*F-statistic Inf"
  )
})
