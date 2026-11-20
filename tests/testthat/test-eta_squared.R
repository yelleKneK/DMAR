# Tests of eta_squared(): the total and the partial estimand from a fitted
# model, from sums of squares, and from F. The tests of
# eta_squared_partial(), removed in 1.1.0, are folded in under
# estimand = "partial". Every expected value is computed here from the
# formulas, independently of the package, or is a published value.

# ---------- The F interface (partial only) ----------

test_that("the F interface gives the partial value of the F-and-df formula", {
  res <- eta_squared(F_value = 11.221, df_effect = 4, df_error = 50)
  expected <- 4 * 11.221 / (4 * 11.221 + 50)
  expect_equal(res$eta_squared_partial, expected, tolerance = 1e-12)
  expect_equal(res$effect, "overall")
  expect_true(is.na(res$eta_squared_total))
  expect_named(res, c("effect", "eta_squared_total", "eta_squared_partial",
                      "F_value", "df_effect", "df_error"))
  expect_s3_class(res, c("dmar_anova_es", "dmar_tbl", "data.frame"),
                  exact = TRUE)
  expect_identical(attr(res, "estimand"), "all")
  expect_named(attr(res, "estimand_notes"), "eta_squared_total")
  expect_true(any(grepl("^NA eta_squared_total: F and its degrees of freedom",
                        utils::capture.output(print(res)))))
  part <- eta_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                      estimand = "partial")
  expect_named(part, c("effect", "eta_squared_partial", "F_value",
                       "df_effect", "df_error"))
  expect_identical(part$eta_squared_partial, res$eta_squared_partial)
  expect_null(attr(part, "estimand_notes"))
})

test_that("the total is refused when asked for by name from F", {
  e <- tryCatch(eta_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                            estimand = "total"), error = function(e) e)
  expect_s3_class(e, c("dmar_estimand_undetermined", "dmar_refusal"),
                  exact = FALSE)
  expect_identical(e$code, "estimand_undetermined")
  expect_identical(e$estimand, "total")
  expect_match(conditionMessage(e), "do not determine the total")
})

test_that("several F statistics against one error term give one row each", {
  Fs <- c(A = 6.934, B = 10.979, "A:B" = 2.504)
  res <- eta_squared(F_value = Fs, df_effect = c(1, 2, 2), df_error = 24)
  expect_identical(res$effect, names(Fs))
  num <- c(1, 2, 2) * Fs
  expect_equal(res$eta_squared_partial, unname(num / (num + 24)),
               tolerance = 1e-12)
  expect_identical(eta_squared(F_value = c(2, 3), df_effect = 1,
                               df_error = 20)$effect,
                   c("effect_1", "effect_2"))
})

test_that("F = 0 yields 0, and the partial value lies in [0, 1]", {
  expect_equal(eta_squared(F_value = 0, df_effect = 2,
                           df_error = 30)$eta_squared_partial, 0)
  set.seed(113)
  Fs <- runif(50, 0.1, 100)
  df1 <- sample.int(10, 50, replace = TRUE)
  res <- eta_squared(F_value = Fs, df_effect = df1, df_error = 100)
  expect_true(all(res$eta_squared_partial >= 0 &
                    res$eta_squared_partial <= 1))
})


# ---------- Sums of squares: the Kim (2016) chain ----------

test_that("Kim (2016) Table 2: the one-way eta squared from the sums of squares", {
  # One-way ANOVA of MATERIAL: SS 1996.998 (3 df), error 3866.717 (76 df),
  # total 5863.715 = their sum. Total and partial coincide.
  res <- eta_squared(SS_effect = 1996.998, SS_error = 3866.717,
                     SS_total = 5863.715, df_effect = 3, df_error = 76)
  expect_equal(res$eta_squared_total, 1996.998 / 5863.715, tolerance = 1e-12)
  expect_equal(res$eta_squared_partial, 1996.998 / (1996.998 + 3866.717),
               tolerance = 1e-12)
  expect_equal(res$eta_squared_total, 0.340568735008, tolerance = 1e-11)
  expect_equal(res$eta_squared_total, res$eta_squared_partial,
               tolerance = 1e-12)
  # The F statistic, from the sums of squares; Kim prints 13.084.
  expect_equal(res$F_value, (1996.998 / 3) / (3866.717 / 76),
               tolerance = 1e-12)
  expect_equal(round(res$F_value, 3), 13.084)
  # The F interface with the printed F gives the rounded value.
  expect_equal(round(eta_squared(F_value = 13.084, df_effect = 3,
                                 df_error = 76)$eta_squared_partial, 3),
               0.341)
})

test_that("Kim (2016) Table 3: total and partial eta squared of a 2 by 4 design", {
  k <- .kim_table3
  res <- eta_squared(SS_effect = k$SS_effect, SS_error = k$SS_error,
                     SS_total = k$SS_total, df_effect = k$df_effect,
                     df_error = k$df_error)
  expect_identical(res$effect, names(k$SS_effect))
  ss <- unname(k$SS_effect)
  expect_equal(res$eta_squared_total, ss / 5863.715, tolerance = 1e-12)
  expect_equal(res$eta_squared_partial, ss / (ss + 2261.612),
               tolerance = 1e-12)
  # The values of the planning record (numerical_verification.txt), to
  # twelve digits.
  expect_equal(res$eta_squared_total,
               c(0.00592047874087, 0.34056873500844, 0.26781468744644),
               tolerance = 1e-11)
  expect_equal(res$eta_squared_partial,
               c(0.0151180493379, 0.468931881529, 0.409809131052),
               tolerance = 1e-11)
  # Kim's printed partial eta squared, .015, .469, .410 (a rounded
  # comparison at the printed precision).
  expect_equal(round(res$eta_squared_partial, 3), c(0.015, 0.469, 0.410))
  # The F interface, with F from the sums of squares, gives the partial
  # values and only them.
  Fs <- (ss / k$df_effect) / (2261.612 / 72)
  viaF <- eta_squared(F_value = Fs, df_effect = k$df_effect, df_error = 72)
  expect_equal(viaF$eta_squared_partial, res$eta_squared_partial,
               tolerance = 1e-12)
  expect_true(all(is.na(viaF$eta_squared_total)))
  # In this balanced design the totals of the effects and the error's
  # share add up to 1.
  expect_equal(sum(res$eta_squared_total) + 2261.612 / 5863.715, 1,
               tolerance = 1e-12)
})

test_that("Olejnik and Algina (2000) Table 10: the treatment factor of a fixed-effects factorial", {
  # SS_treatment = 2 x 54.061 (the table prints the numerator as 88.122,
  # a misprint: 88.122 / 748.621 is .118, not the .144 it reports),
  # SS_error = 431.273, SS_total = 748.621. A rounded comparison.
  res <- eta_squared(SS_effect = 2 * 54.061, SS_error = 431.273,
                     SS_total = 748.621)
  expect_equal(round(res$eta_squared_total, 3), 0.144)
  expect_equal(round(res$eta_squared_partial, 3), 0.200)
  # With only one of SS_error and SS_total, the other estimand is NA.
  tot <- eta_squared(SS_effect = 108.122, SS_total = 748.621)
  expect_true(is.na(tot$eta_squared_partial))
  expect_match(attr(tot, "estimand_notes")[["eta_squared_partial"]],
               "SS_error was not given")
  expect_error(eta_squared(SS_effect = 108.122, SS_total = 748.621,
                           estimand = "partial"),
               class = "dmar_estimand_undetermined")
  part <- eta_squared(SS_effect = 108.122, SS_error = 431.273)
  expect_true(is.na(part$eta_squared_total))
  expect_true(is.na(part$F_value))
  expect_equal(eta_squared(SS_effect = 108.122, SS_total = 748.621,
                           estimand = "total")$eta_squared_total,
               108.122 / 748.621, tolerance = 1e-12)
})


# ---------- Fitted models ----------

test_that("Maxwell, Delaney, and Kelley (2027) Tables 7.5 to 7.8: the blood pressure data", {
  # The raw data of Table 7.5 reproduce the sums of squares of Table 7.7
  # exactly, and the R squared columns of Table 7.8 to two decimals.
  fit <- aov(y ~ A * B, data = .mdk_blood_pressure())
  tab <- anova(fit)
  expect_equal(unname(tab[["Sum Sq"]]), c(1080, 3420, 780, 3738),
               tolerance = 1e-10)
  res <- eta_squared(fit)
  expect_identical(res$effect, c("A", "B", "A:B"))
  ss <- c(1080, 3420, 780)
  expect_equal(res$eta_squared_total, ss / 9018, tolerance = 1e-12)
  expect_equal(res$eta_squared_partial, ss / (ss + 3738), tolerance = 1e-12)
  expect_equal(round(res$eta_squared_total, 2), c(0.12, 0.38, 0.09))
  expect_equal(round(res$eta_squared_partial, 2), c(0.22, 0.48, 0.17))
  # Eq. 7.43: the total from the F statistics of every effect.
  dfF <- c(1, 2, 2) * res$F_value
  expect_equal(res$eta_squared_total, dfF / (sum(dfF) + 24),
               tolerance = 1e-12)
  # The sums of squares interface with the same values agrees.
  ss_res <- eta_squared(SS_effect = c(A = 1080, B = 3420, "A:B" = 780),
                        SS_error = 3738, SS_total = 9018,
                        df_effect = c(1, 2, 2), df_error = 24)
  expect_equal(ss_res$eta_squared_total, res$eta_squared_total,
               tolerance = 1e-12)
  expect_equal(ss_res$eta_squared_partial, res$eta_squared_partial,
               tolerance = 1e-12)
  expect_equal(ss_res$F_value, res$F_value, tolerance = 1e-10)
})

test_that("a one-way fit gives the same total and partial value, SS_effect / SS_total", {
  fit <- aov(bdi_post ~ condition, data = depression_bdi)
  tbl <- anova(fit)
  res <- eta_squared(fit)
  expect_equal(res$effect, "condition")
  expect_equal(res$eta_squared_total,
               tbl[["Sum Sq"]][1] / sum(tbl[["Sum Sq"]]), tolerance = 1e-12)
  expect_equal(res$eta_squared_partial, res$eta_squared_total,
               tolerance = 1e-12)
  raw <- eta_squared(F_value = tbl[["F value"]][1], df_effect = 2,
                     df_error = 27)
  expect_equal(raw$eta_squared_partial, res$eta_squared_partial,
               tolerance = 1e-12)
  expect_null(attr(res, "estimand_notes"))
})

test_that("an unbalanced factorial fit: sequential sums of squares, whose totals sum with the error's share to 1", {
  fit <- aov(iq_8 ~ treatment * factor(grade), data = pygmalion)
  tbl <- anova(fit)
  res <- eta_squared(fit)
  expect_equal(nrow(res), 3L)
  ss <- tbl[["Sum Sq"]]
  expect_equal(res$eta_squared_total, ss[1:3] / sum(ss), tolerance = 1e-12)
  expect_equal(res$eta_squared_partial, ss[1:3] / (ss[1:3] + ss[4]),
               tolerance = 1e-12)
  expect_equal(sum(res$eta_squared_total) + ss[4] / sum(ss), 1,
               tolerance = 1e-12)
  # The order of the terms changes the Type I values in an unbalanced
  # design.
  rev_fit <- aov(iq_8 ~ factor(grade) * treatment, data = pygmalion)
  expect_false(isTRUE(all.equal(
    eta_squared(rev_fit)$eta_squared_total[2], res$eta_squared_total[1])))
})

test_that("an lm fit, with a covariate, and an lm fit without an intercept", {
  fit_lm <- lm(bdi_post ~ condition, data = depression_bdi)
  fit_aov <- aov(bdi_post ~ condition, data = depression_bdi)
  expect_equal(eta_squared(fit_lm)$eta_squared_total,
               eta_squared(fit_aov)$eta_squared_total, tolerance = 1e-12)
  cov_fit <- lm(iq_8 ~ iq_pre + treatment, data = pygmalion)
  res <- eta_squared(cov_fit)
  tbl <- anova(cov_fit)
  expect_identical(res$effect, c("iq_pre", "treatment"))
  expect_equal(res$eta_squared_total, tbl[["Sum Sq"]][1:2] /
                 sum(tbl[["Sum Sq"]]), tolerance = 1e-12)
  # Without an intercept the sums of squares do not partition the
  # corrected total.
  noint <- lm(bdi_post ~ 0 + condition, data = depression_bdi)
  r0 <- eta_squared(noint)
  expect_true(is.na(r0$eta_squared_total))
  expect_match(attr(r0, "estimand_notes")[["eta_squared_total"]],
               "no intercept")
  expect_error(eta_squared(noint, estimand = "total"),
               class = "dmar_estimand_undetermined")
})

test_that("an aovlist fit: the partial value per stratum, the total over every stratum", {
  set.seed(113)
  d <- data.frame(id = factor(rep(1:16, each = 3)),
                  time = factor(rep(1:3, 16)),
                  group = factor(rep(c("a", "b"), each = 24)))
  d$y <- rnorm(16, sd = 1.5)[as.integer(d$id)] + 0.5 * as.integer(d$time) +
    rnorm(48)
  fit <- aov(y ~ group * time + Error(id / time), data = d)
  res <- eta_squared(fit)
  expect_identical(res$effect, c("group", "time", "group:time"))
  expect_identical(res$stratum, c("id", "id:time", "id:time"))
  s <- summary(fit)
  sub_t <- s[["Error: id"]][[1]]
  wit <- s[["Error: id:time"]][[1]]
  rownames(sub_t) <- trimws(rownames(sub_t))
  rownames(wit) <- trimws(rownames(wit))
  ss_y <- sum((d$y - mean(d$y))^2)
  expect_equal(res$eta_squared_total,
               c(sub_t["group", "Sum Sq"], wit["time", "Sum Sq"],
                 wit["group:time", "Sum Sq"]) / ss_y, tolerance = 1e-12)
  expect_equal(res$eta_squared_partial[2],
               wit["time", "Sum Sq"] /
                 (wit["time", "Sum Sq"] + wit["Residuals", "Sum Sq"]),
               tolerance = 1e-12)
  expect_true(all(res$eta_squared_total < res$eta_squared_partial))
  # The corrected total over the strata is sum((y - mean(y))^2), balanced
  # or not (Bakeman, 2005, p. 380: the denominator of eta squared holds
  # every component).
  expect_equal(DMAR:::.aovlist_total_ss(fit), ss_y, tolerance = 1e-10)
  d2 <- d[-c(2, 9, 17, 30, 44), ]
  fit2 <- suppressWarnings(aov(y ~ group * time + Error(id / time),
                               data = d2))
  expect_equal(DMAR:::.aovlist_total_ss(fit2),
               sum((d2$y - mean(d2$y))^2), tolerance = 1e-10)
})


# ---------- Inputs that conflict or do not suffice ----------

test_that("incomplete and conflicting inputs are refused", {
  k <- .kim_table3
  expect_error(eta_squared(), "Give a fitted aov, lm, or aovlist model",
               class = "dmar_argument_invalid")
  expect_error(eta_squared(F_value = 1.5, df_effect = 2),
               "'F_value' was given without 'df_error'")
  expect_error(eta_squared(F_value = -1, df_effect = 2, df_error = 10),
               "'F_value' must be")
  expect_error(eta_squared(F_value = 1, df_effect = 0, df_error = 10),
               "'df_effect' must be")
  expect_error(eta_squared(F_value = 1, df_effect = 2, df_error = -1),
               "'df_error' must be")
  expect_error(eta_squared(object = list(foo = 1)), "aov, lm, or aovlist")
  expect_error(eta_squared(SS_effect = 10), "Give 'SS_total'")
  expect_error(eta_squared(SS_error = 10, SS_total = 20),
               "needs 'SS_effect'")
  expect_error(eta_squared(SS_effect = 10, SS_error = 20, F_value = 2),
               "were both given")
  expect_error(eta_squared(aov(bdi_post ~ condition, data = depression_bdi),
                           SS_total = 3),
               "were both given")
  expect_error(eta_squared(SS_effect = 10, SS_error = 20, SS_total = 25),
               "smaller than SS_effect \\+ SS_error",
               class = "dmar_argument_invalid")
  expect_error(eta_squared(SS_effect = 10, SS_total = 20, N = 30),
               "unused argument")
  expect_error(eta_squared(SS_effect = k$SS_effect, SS_error = k$SS_error,
                           df_effect = c(1, 3)),
               "one value for all effects or one for each")
  expect_error(eta_squared(SS_effect = 1, SS_total = 2, estimand = "both"),
               "'estimand' must be one of")
  expect_error(eta_squared(SS_effect = 1, SS_total = 2,
                           estimand = c("total", "partial")),
               "'estimand' must be one of")
  # A published table whose rounded sums of squares add up exactly is
  # accepted (Kim's one-way table: 1996.998 + 3866.717 = 5863.715).
  expect_error(eta_squared(SS_effect = 1996.998, SS_error = 3866.717,
                           SS_total = 5863.715), NA)
})


# ---------- tidy() and glance() ----------

test_that("tidy() gives one row per effect and estimand; glance() one row", {
  fit <- aov(y ~ A * B, data = .mdk_blood_pressure())
  res <- eta_squared(fit)
  td <- generics::tidy(res)
  expect_named(td, c("term", "effect", "estimate"))
  expect_identical(td$term, rep(c("eta_squared_total", "eta_squared_partial"),
                                each = 3))
  expect_identical(td$estimate, c(res$eta_squared_total,
                                  res$eta_squared_partial))
  gl <- generics::glance(res)
  expect_identical(gl$term, "eta_squared_total, eta_squared_partial")
  expect_identical(gl$n_effects, 3L)
  one <- eta_squared(F_value = 3, df_effect = 2, df_error = 20,
                     estimand = "partial")
  expect_identical(generics::glance(one), generics::tidy(one))
  # A fit with an Error() term names the stratum.
  set.seed(113)
  d <- data.frame(s = factor(rep(1:10, each = 3)),
                  w = factor(rep(1:3, 10)), y = rnorm(30))
  tw <- generics::tidy(eta_squared(aov(y ~ w + Error(s/w), data = d)))
  expect_identical(names(tw), c("term", "effect", "stratum", "estimate"))
})
