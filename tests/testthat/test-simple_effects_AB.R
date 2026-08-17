test_that("simple_effects_AB() returns a tidy data.frame with the documented columns", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  res <- suppressWarnings(simple_effects_AB(fit))
  expect_s3_class(res, "data.frame")
  expect_named(
    res,
    c("effect", "focal_factor", "conditioning_factor", "conditioning_level",
      "F_value", "df_effect", "df_error", "p_value", "p_adjusted",
      "partial_eta_squared", "lower_limit", "upper_limit", "n_at_level")
  )
  expect_identical(attr(res, "factor_A"), "wool")
  expect_identical(attr(res, "factor_B"), "tension")
  expect_identical(attr(res, "error_term"), "pooled")
  expect_identical(attr(res, "adjust"), "none")
  expect_equal(attr(res, "conf_level"), 0.95)
})

test_that("which = 'both' returns a + b rows in A_at_B then B_at_A order", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  res <- suppressWarnings(simple_effects_AB(fit))
  expect_equal(nrow(res), 2 + 3)  # 3 wool|tension + 2 tension|wool
  expect_equal(res$focal_factor,
               c("wool", "wool", "wool", "tension", "tension"))
  expect_equal(res$conditioning_factor,
               c("tension", "tension", "tension", "wool", "wool"))
})

test_that("which = 'A_at_B' / 'B_at_A' filter to the right family", {
  fit  <- aov(breaks ~ wool * tension, data = warpbreaks)
  AatB <- suppressWarnings(simple_effects_AB(fit, which = "A_at_B"))
  BatA <- suppressWarnings(simple_effects_AB(fit, which = "B_at_A"))
  expect_equal(nrow(AatB), 3L)
  expect_equal(nrow(BatA), 2L)
  expect_true(all(AatB$focal_factor == "wool"))
  expect_true(all(BatA$focal_factor == "tension"))
})

test_that("pooled F matches MS_simple / MS_W from the full model by hand", {
  fit       <- aov(breaks ~ wool * tension, data = warpbreaks)
  ms_w_full <- anova(fit)["Residuals", "Mean Sq"]
  df_w_full <- anova(fit)["Residuals", "Df"]

  # Hand-compute F for tension | wool = A.
  sub      <- subset(warpbreaks, wool == "A")
  ss_T_A   <- anova(aov(breaks ~ tension, data = sub))["tension", "Sum Sq"]
  F_T_A    <- (ss_T_A / 2) / ms_w_full

  # And for wool | tension = L.
  sub      <- subset(warpbreaks, tension == "L")
  ss_W_L   <- anova(aov(breaks ~ wool,    data = sub))["wool",    "Sum Sq"]
  F_W_L    <- (ss_W_L / 1) / ms_w_full

  res <- suppressWarnings(simple_effects_AB(fit))
  expect_equal(res$F_value[res$effect == "tension | wool = A"], F_T_A,
               tolerance = 1e-10)
  expect_equal(res$F_value[res$effect == "wool | tension = L"], F_W_L,
               tolerance = 1e-10)
  expect_true(all(res$df_error == df_w_full))
})

test_that("welch F matches stats::oneway.test on the subset", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  res <- suppressWarnings(simple_effects_AB(fit, error_term = "welch"))

  sub <- subset(warpbreaks, tension == "L")
  ow  <- oneway.test(breaks ~ wool, data = sub, var.equal = FALSE)
  row <- res[res$effect == "wool | tension = L", ]
  expect_equal(row$F_value,   as.numeric(ow$statistic),    tolerance = 1e-10)
  expect_equal(row$df_effect, as.numeric(ow$parameter[1]), tolerance = 1e-10)
  expect_equal(row$df_error,  as.numeric(ow$parameter[2]), tolerance = 1e-10)
  expect_equal(row$p_value,   ow$p.value,                  tolerance = 1e-10)
})

test_that("partial_eta_squared equals (df_eff * F) / (df_eff * F + df_err)", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  res <- suppressWarnings(simple_effects_AB(fit))
  expect_equal(res$partial_eta_squared,
               (res$df_effect * res$F_value) /
                 (res$df_effect * res$F_value + res$df_error),
               tolerance = 1e-12)
})

test_that("bonferroni p_adjusted = min(1, m * p_value) across the full family", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  res <- suppressWarnings(simple_effects_AB(fit, adjust = "bonferroni"))
  m   <- nrow(res)
  expect_equal(res$p_adjusted, pmin(1, res$p_value * m))
})

test_that("sequential adjustments match stats::p.adjust over the family", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  for (m in c("holm", "hochberg", "BH", "BY")) {
    res <- suppressWarnings(simple_effects_AB(fit, adjust = m))
    expect_equal(res$p_adjusted, p.adjust(res$p_value, method = m),
                 tolerance = 1e-12, info = paste("adjust =", m))
  }
})

test_that("which family is the adjustment scope (not the full a + b family)", {
  fit  <- aov(breaks ~ wool * tension, data = warpbreaks)
  full <- suppressWarnings(simple_effects_AB(fit, adjust = "bonferroni"))
  AatB <- suppressWarnings(simple_effects_AB(fit, adjust = "bonferroni",
                                             which = "A_at_B"))
  # In the full call there are 5 tests; in the A_at_B-only call there are 3.
  # Multiplier differs, so p_adjusted differs.
  p_full_row <- full$p_adjusted[full$effect == "wool | tension = L"]
  p_AB_row   <- AatB$p_adjusted[AatB$effect == "wool | tension = L"]
  expect_equal(p_full_row, pmin(1, full$p_value[full$effect == "wool | tension = L"] * 5))
  expect_equal(p_AB_row,   pmin(1, AatB$p_value[AatB$effect == "wool | tension = L"] * 3))
})

test_that("CI lower and upper limits are in [0, 1] and bracket the point estimate", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  res <- suppressWarnings(simple_effects_AB(fit))
  ok  <- !is.na(res$F_value)
  expect_true(all(res$lower_limit[ok] >= 0))
  expect_true(all(res$upper_limit[ok] <= 1))
  expect_true(all(res$lower_limit[ok] <= res$partial_eta_squared[ok]))
  expect_true(all(res$partial_eta_squared[ok] <= res$upper_limit[ok]))
})

test_that("conf_level controls the interval width", {
  fit  <- aov(breaks ~ wool * tension, data = warpbreaks)
  r90  <- suppressWarnings(simple_effects_AB(fit, conf_level = 0.90))
  r99  <- suppressWarnings(simple_effects_AB(fit, conf_level = 0.99))
  w90  <- r90$upper_limit - r90$lower_limit
  w99  <- r99$upper_limit - r99$lower_limit
  ok   <- !is.na(w90) & !is.na(w99)
  expect_true(all(w99[ok] >= w90[ok] - 1e-9))
})

test_that("no-interaction model issues a warning naming the two factors", {
  fit <- aov(breaks ~ wool + tension, data = warpbreaks)
  expect_warning(
    suppressWarnings(simple_effects_AB(fit), classes = "warning") |>
      invisible(),
    NA  # placeholder: real check below
  )
  # Capture the first warning's message.
  msgs <- character()
  withCallingHandlers(
    simple_effects_AB(fit),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("no interaction term", msgs)))
  expect_true(any(grepl("wool \\* tension", msgs)))
})

test_that("non-aov/lm input is rejected with a clear error", {
  expect_error(simple_effects_AB(lm(breaks ~ 1, data = warpbreaks)),
               "two-factor design")
  expect_error(simple_effects_AB("not a model"), "fitted aov or lm")
})

test_that("conf_level outside (0, 1) is rejected", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  expect_error(simple_effects_AB(fit, conf_level = 0),  "conf_level")
  expect_error(simple_effects_AB(fit, conf_level = 1),  "conf_level")
  expect_error(simple_effects_AB(fit, conf_level = -1), "conf_level")
})

test_that("unbalanced design uses correct n_at_level per row", {
  d   <- warpbreaks[-c(1:5, 28:31), ]   # drop rows to create imbalance
  fit <- aov(breaks ~ wool * tension, data = d)
  res <- suppressWarnings(simple_effects_AB(fit))
  # n_at_level == nrow of the conditioning subset
  for (i in seq_len(nrow(res))) {
    cond_f <- res$conditioning_factor[i]
    cond_l <- res$conditioning_level[i]
    expect_equal(res$n_at_level[i], sum(d[[cond_f]] == cond_l))
  }
})

test_that("simple effect tests are equivalent to one-way ANOVAs (pooled error)", {
  # For a balanced design, F_pooled(A | b_j) computed by simple_effects_AB
  # should equal SS_A(subset) / (a-1) / MS_W_full.
  fit       <- aov(breaks ~ wool * tension, data = warpbreaks)
  res       <- suppressWarnings(simple_effects_AB(fit))
  ms_w_full <- anova(fit)["Residuals", "Mean Sq"]
  for (lvl in levels(warpbreaks$tension)) {
    sub      <- subset(warpbreaks, tension == lvl)
    F_hand   <- (anova(aov(breaks ~ wool, data = sub))["wool", "Sum Sq"]) / ms_w_full
    F_func   <- res$F_value[res$effect == paste0("wool | tension = ", lvl)]
    expect_equal(F_func, F_hand, tolerance = 1e-10)
  }
})

test_that("warning deduplication emits one summary warning, not one per row", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  msgs <- character()
  withCallingHandlers(
    simple_effects_AB(fit),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  clamp_msgs <- msgs[grepl("clamp", msgs)]
  expect_length(clamp_msgs, 1L)
  expect_match(clamp_msgs, "[0-9]+ of the simple effect rows")
})
