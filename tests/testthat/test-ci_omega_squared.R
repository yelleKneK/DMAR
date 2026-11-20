test_that("point estimate matches the closed-form formula for Bargmann (1970)", {
  res <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  expected <- 4 * (11.221 - 1) / (4 * (11.221 - 1) + 55)
  expect_equal(res$omega_squared_partial, expected, tolerance = 1e-6)
  expect_equal(res$effect, "overall")
})

test_that("raw-interface CI agrees with the direct Steiger transformation of ci_pvaf", {
  res <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  pvaf <- ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
  # Both apply ncp / (ncp + N); the bounds should match to high precision.
  expect_equal(res$omega_squared_partial_ci_lower,
               pvaf$value[pvaf$term == "ci_lower"], tolerance = 1e-6)
  expect_equal(res$omega_squared_partial_ci_upper,
               pvaf$value[pvaf$term == "ci_upper"], tolerance = 1e-6)
})

test_that("CI brackets the point estimate", {
  res <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  expect_gte(res$omega_squared_partial, res$omega_squared_partial_ci_lower)
  expect_lte(res$omega_squared_partial, res$omega_squared_partial_ci_upper)
})

test_that("conf_level narrows with smaller coverage and widens with larger", {
  r90 <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
                          conf_level = 0.90)
  r99 <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
                          conf_level = 0.99)
  expect_gt(r99$omega_squared_partial_ci_upper - r99$omega_squared_partial_ci_lower,
            r90$omega_squared_partial_ci_upper - r90$omega_squared_partial_ci_lower)
})

test_that("low F (CI includes 0) yields ci_lower == 0", {
  # Deliberately non-significant: F close to 1
  expect_warning(
    res <- ci_omega_squared(F_value = 0.5, df_effect = 2, df_error = 30, N = 33),
    "below the alpha_lower critical value")
  expect_equal(res$omega_squared_partial_ci_lower, 0)
  expect_equal(res$omega_squared_partial, 0)  # truncated at 0 as F < 1
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

  expect_equal(res_model$omega_squared_partial, res_raw$omega_squared_partial)
  expect_equal(res_model$omega_squared_partial_ci_lower,   res_raw$omega_squared_partial_ci_lower)
  expect_equal(res_model$omega_squared_partial_ci_upper,   res_raw$omega_squared_partial_ci_upper)
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
    expect_equal(res$omega_squared_partial[res$effect == eff], expected, tolerance = 1e-9)
  }
})

test_that("supports an lm() fit as well as aov()", {
  fit_lm  <- lm(weight ~ group, data = PlantGrowth)
  fit_aov <- aov(weight ~ group, data = PlantGrowth)
  res_lm  <- ci_omega_squared(fit_lm)
  res_aov <- ci_omega_squared(fit_aov)
  expect_equal(res_lm$omega_squared_partial, res_aov$omega_squared_partial)
  expect_equal(res_lm$omega_squared_partial_ci_lower,   res_aov$omega_squared_partial_ci_lower)
  expect_equal(res_lm$omega_squared_partial_ci_upper,   res_aov$omega_squared_partial_ci_upper)
})

test_that("asymmetric alpha arguments are respected", {
  res_sym  <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                               N = 55, conf_level = 0.90)
  res_asym <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                               N = 55, alpha_lower = 0, alpha_upper = 0.10)
  # With alpha_lower = 0 the lower bound on the NCP, and thus omega^2, is 0
  expect_equal(res_asym$omega_squared_partial_ci_lower, 0)
  # Upper limit of asym (alpha_upper = .10) should match upper of sym (alpha_upper = .05)... no.
  # sym has alpha_upper = .05, asym has alpha_upper = .10, so asym_upper < sym_upper.
  expect_lt(res_asym$omega_squared_partial_ci_upper, res_sym$omega_squared_partial_ci_upper)
})

test_that("errors on missing or invalid arguments", {
  expect_error(ci_omega_squared(), "Give a fitted aov, lm, or aovlist model")
  expect_error(ci_omega_squared(F_value = 1.5, df_effect = 2),
               "'F_value' was given without 'df_error' and 'N'")
  expect_error(ci_omega_squared(F_value = -1, df_effect = 2,
                                df_error = 10, N = 20),
               "'F_value' must be")
  expect_error(ci_omega_squared(F_value = 3, df_effect = 2,
                                df_error = 10, N = 11),
               "at least df_effect \\+ df_error")
  expect_error(ci_omega_squared(object = list(foo = 1)), "aov, lm, or aovlist")
})

test_that("returned data.frame has the documented columns", {
  res <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  expect_named(res, c("effect", "omega_squared_total",
                      "omega_squared_total_ci_lower",
                      "omega_squared_total_ci_upper", "omega_squared_partial",
                      "omega_squared_partial_ci_lower",
                      "omega_squared_partial_ci_upper",
                      "F_value", "df_effect", "df_error", "N"))
  expect_true(all(is.na(unlist(res[, 2:4]))))
  part <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                           N = 55, estimand = "partial")
  expect_named(part, c("effect", "omega_squared_partial",
                       "omega_squared_partial_ci_lower",
                       "omega_squared_partial_ci_upper",
                       "F_value", "df_effect", "df_error", "N"))
})

test_that("ci_omega_squared() reports the lower-limit clamp once, in terms of omega squared, and still returns the interval", {
  msgs <- capture_warnings(
    result <- ci_omega_squared(F_value = 1.2, df_effect = 4, df_error = 50, N = 55)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on partial omega squared is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_equal(result$omega_squared_partial_ci_lower, 0)
  expect_gt(result$omega_squared_partial_ci_upper, 0)
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
  expect_match(msgs, "lower confidence limit on partial omega squared is 0")
  expect_equal(result$omega_squared_partial_ci_lower, c(0, 0))
})

test_that("ci_omega_squared() refuses an F statistic that is not finite", {
  expect_error(
    ci_omega_squared(F_value = Inf, df_effect = 4, df_error = 50, N = 55),
    "'F_value' must be one or more finite non-negative numbers",
    class = "dmar_argument_invalid"
  )
})

test_that("ci_omega_squared() records the confidence level of its intervals", {
  # The 1.1.0 review: the table recorded no level, so results_sentence()
  # wrote a 90% interval as "95% CI".
  fit <- aov(len ~ supp * factor(dose), data = ToothGrowth)
  o90 <- ci_omega_squared(fit, conf_level = 0.90)
  expect_identical(attr(o90, "conf_level"), 0.90)
  expect_match(results_sentence(o90[1, ], digits = 3),
               "^omega_squared_partial = 0\\.[0-9]{3}, 90% CI \\[0\\.073, 0\\.350\\]$")
  expect_true("Confidence level: 90%" %in% utils::capture.output(print(o90)))
  expect_identical(attr(ci_omega_squared(F_value = 11.221, df_effect = 4,
                                         df_error = 50, N = 55),
                        "conf_level"), 0.95)
  tails <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                            N = 55, alpha_lower = 0.05, alpha_upper = 0.05)
  expect_equal(attr(tails, "conf_level"), 0.90, tolerance = 1e-12)
})


# ---------- The total and the partial estimand (1.1.0) ----------

test_that("the one-way mood-induction example of Maxwell, Delaney, and Kelley (2027) is reproduced", {
  # Chapter 3 (pp. 120-121 of the 7 July 2026 draft): E_R = 72.67,
  # E_F = 26 on 27 df, a = 3, N = 30; omega-hat squared .608 (Eq. 3.96),
  # F(2, 27) = 24.23, and the 95% interval for omega squared [.362, .753].
  # A rounded comparison at the printed precision.
  F_value <- ((72.67 - 26) / 2) / (26 / 27)
  expect_equal(round(F_value, 2), 24.23)
  res <- ci_omega_squared(F_value = F_value, df_effect = 2, df_error = 27,
                          N = 30, estimand = "partial")
  expect_equal(round(res$omega_squared_partial, 3), 0.608)
  expect_equal(round(c(res$omega_squared_partial_ci_lower,
                       res$omega_squared_partial_ci_upper), 3),
               c(0.362, 0.753))
  # Eq. 3.96 from the sums of squares: the same number.
  expect_equal(res$omega_squared_partial,
               ((72.67 - 26) - 2 * (26 / 27)) / (72.67 + 26 / 27),
               tolerance = 1e-12)
})

test_that("a one-way fit reports the total with the partial value's interval; a factorial fit gives the total no limits", {
  one <- ci_omega_squared(aov(iq_gain ~ factor(grade), data = pygmalion))
  ss <- anova(aov(iq_gain ~ factor(grade), data = pygmalion))
  ms_e <- ss[["Mean Sq"]][2]
  expect_equal(one$omega_squared_total,
               (ss[["Sum Sq"]][1] - 5 * ms_e) / (sum(ss[["Sum Sq"]]) + ms_e),
               tolerance = 1e-12)
  expect_equal(one$omega_squared_total, one$omega_squared_partial,
               tolerance = 1e-12)
  expect_identical(one$omega_squared_total_ci_upper,
                   one$omega_squared_partial_ci_upper)
  fit <- aov(y ~ A * B, data = .mdk_blood_pressure())
  res <- suppressWarnings(ci_omega_squared(fit))
  # Maxwell, Delaney, and Kelley (2027), Table 7.8, to two decimals.
  expect_equal(round(res$omega_squared_total, 2), c(0.10, 0.34, 0.05))
  expect_equal(round(res$omega_squared_partial, 2), c(0.17, 0.40, 0.09))
  expect_true(all(is.na(res$omega_squared_total_ci_lower)))
  expect_error(ci_omega_squared(fit, estimand = "total"),
               class = "dmar_estimand_undetermined")
  # The partial interval is the one of ci_eta_squared().
  eta <- suppressWarnings(ci_eta_squared(fit, estimand = "partial"))
  expect_identical(res$omega_squared_partial_ci_upper,
                   eta$eta_squared_partial_ci_upper)
})

test_that("ci_omega_squared() takes an aovlist fit, computing each effect with its stratum's size", {
  n <- nrow(pygmalion)
  occasions <- c("pretest", "intermediate", "end_of_study")
  iq_long <- data.frame(
    id       = factor(rep(seq_len(n), times = 3)),
    occasion = factor(rep(occasions, each = n), levels = occasions),
    iq       = c(pygmalion$iq_pre, pygmalion$iq_4, pygmalion$iq_8))
  fit <- aov(iq ~ occasion + Error(id / occasion), data = iq_long)
  res <- ci_omega_squared(fit)
  expect_identical(res$stratum, "id:occasion")
  expect_identical(res$N, 930)
  expect_true(is.na(res$omega_squared_total))
  num <- 2 * (res$F_value - 1)
  expect_equal(res$omega_squared_partial, num / (num + 620),
               tolerance = 1e-12)
  expect_equal(round(res$omega_squared_partial, 3), 0.209)
  eta <- ci_eta_squared(fit, estimand = "partial")
  expect_identical(res$omega_squared_partial_ci_lower,
                   eta$eta_squared_partial_ci_lower)
  raw <- ci_omega_squared(F_value = res$F_value, df_effect = 2,
                          df_error = 618, N = 620, estimand = "partial")
  expect_identical(raw$omega_squared_partial, res$omega_squared_partial)
  expect_identical(raw$omega_squared_partial_ci_upper,
                   res$omega_squared_partial_ci_upper)
  expect_error(ci_omega_squared(fit, estimand = "total"),
               class = "dmar_estimand_undetermined")
})
