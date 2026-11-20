# Tests of omega_squared(): the total and the partial estimand from a
# fitted model, from sums of squares, and from F. The tests of
# omega_squared_partial(), removed in 1.1.0, are folded in under
# estimand = "partial". Every expected value is computed here from the
# formulas (Olejnik & Algina, 2000, Table 10; Maxwell, Delaney, & Kelley,
# 2027, Eqs. 7.44 to 7.47), independently of the package, or is a
# published value.

# The two estimators, coded from the published formulas.
.omega_total <- function(ss, df, ms_e, ss_t) {
  pmax(0, (ss - df * ms_e) / (ss_t + ms_e))
}
.omega_partial <- function(ss, df, ms_e, N) {
  pmax(0, (ss - df * ms_e) / (ss + (N - df) * ms_e))
}

# ---------- The F interface (partial only) ----------

test_that("the F interface gives df (F - 1) / (df (F - 1) + N), the total NA", {
  res <- omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  num <- 4 * (11.221 - 1)
  expect_equal(res$omega_squared_partial, num / (num + 55), tolerance = 1e-12)
  expect_true(is.na(res$omega_squared_total))
  expect_named(res, c("effect", "omega_squared_total",
                      "omega_squared_partial", "F_value", "df_effect",
                      "df_error", "N"))
  expect_s3_class(res, c("dmar_anova_es", "dmar_tbl", "data.frame"),
                  exact = TRUE)
  expect_equal(res$effect, "overall")
  # Bargmann's example: .426.
  expect_equal(round(res$omega_squared_partial, 3), 0.426)
  expect_error(omega_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                             N = 55, estimand = "total"),
               class = "dmar_estimand_undetermined")
  part <- omega_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                        N = 55, estimand = "partial")
  expect_named(part, c("effect", "omega_squared_partial", "F_value",
                       "df_effect", "df_error", "N"))
})

test_that("a negative estimate is reported as 0, and F = 0 is accepted", {
  res <- omega_squared(F_value = 0.5, df_effect = 2, df_error = 50, N = 55)
  expect_equal(res$omega_squared_partial, 0)
  expect_equal(omega_squared(F_value = 0, df_effect = 2, df_error = 50,
                             N = 55)$omega_squared_partial, 0)
  # The total too: SS_effect below df_effect MS_error.
  tot <- omega_squared(SS_effect = 1, df_effect = 2, MS_error = 1,
                       SS_total = 50, estimand = "total")
  expect_equal(tot$omega_squared_total, 0)
})

test_that("the size of a one-effect stratum is accepted as N, and nothing smaller", {
  # A within-subjects stratum that tests one effect holds df_effect +
  # df_error degrees of freedom. DMAR 1.0.x refused N = df_effect + df_error
  # in omega_squared().
  at_size <- omega_squared(F_value = 83.11, df_effect = 2, df_error = 618,
                           N = 620)
  expect_equal(at_size$omega_squared_partial, 2 * 82.11 / (2 * 82.11 + 620),
               tolerance = 1e-12)
  expect_error(omega_squared(F_value = 83.11, df_effect = 2, df_error = 618,
                             N = 619), "at least df_effect \\+ df_error")
  expect_error(omega_squared(F_value = 83.11, df_effect = 2, df_error = 618,
                             N = NA_real_), "at least df_effect \\+ df_error")
  expect_error(omega_squared(F_value = 3, df_effect = 2, df_error = 10,
                             N = 12), NA)
  expect_error(omega_squared(F_value = 3, df_effect = 2, df_error = 10,
                             N = 11), "at least df_effect \\+ df_error")
})


# ---------- Sums of squares: the Kim (2016) chain ----------

test_that("Kim (2016) Table 3: total omega squared with the plus denominator, and partial omega squared", {
  k <- .kim_table3
  res <- omega_squared(SS_effect = k$SS_effect, df_effect = k$df_effect,
                       SS_error = k$SS_error, df_error = k$df_error,
                       SS_total = k$SS_total, N = k$N)
  ss <- unname(k$SS_effect)
  ms_e <- 2261.612 / 72
  expect_equal(res$omega_squared_total,
               .omega_total(ss, c(1, 3, 3), ms_e, 5863.715),
               tolerance = 1e-12)
  expect_equal(res$omega_squared_partial,
               .omega_partial(ss, c(1, 3, 3), ms_e, 80), tolerance = 1e-12)
  # The values of the planning record (numerical_verification.txt).
  expect_equal(res$omega_squared_total,
               c(0.000560585484772, 0.322769025973084, 0.250402637214264),
               tolerance = 1e-11)
  expect_equal(res$omega_squared_partial,
               c(0.00131337457339, 0.430912120412, 0.370051123049),
               tolerance = 1e-11)
  # A source discrepancy, not an oracle: Kim's Table 4 divides by
  # SS_total - MS_error (5832.31) and prints .001, .326, .253; the
  # estimator's denominator is SS_total + MS_error (Olejnik & Algina,
  # 2000, p. 266 and Table 10; Maxwell, Delaney, & Kelley, 2027, Eq. 7.45),
  # which gives .001, .323, .250. Olejnik and Algina (2000, p. 262) print
  # the minus sign in the formula for a contrast, while their Table 6 adds.
  minus <- (ss - c(1, 3, 3) * ms_e) / (5863.715 - ms_e)
  expect_equal(round(minus, 3), c(0.001, 0.326, 0.253))
  expect_equal(round(res$omega_squared_total, 3), c(0.001, 0.323, 0.250))
  expect_false(isTRUE(all.equal(res$omega_squared_total, minus)))
  # MS_error in place of SS_error and df_error gives the same values, and
  # F, from the sums of squares, gives the partial values.
  ms <- omega_squared(SS_effect = k$SS_effect, df_effect = k$df_effect,
                      MS_error = ms_e, SS_total = k$SS_total, N = 80)
  expect_equal(ms$omega_squared_total, res$omega_squared_total,
               tolerance = 1e-12)
  expect_true(all(is.na(ms$df_error)))
  Fs <- (ss / c(1, 3, 3)) / ms_e
  viaF <- omega_squared(F_value = Fs, df_effect = c(1, 3, 3),
                        df_error = 72, N = 80)
  expect_equal(viaF$omega_squared_partial, res$omega_squared_partial,
               tolerance = 1e-12)
})

test_that("Olejnik and Algina (2000) Table 10 and Maxwell, Delaney, and Kelley (2027) Table 7.8", {
  # Olejnik and Algina, treatment factor: total .124, partial .165
  # (N = 66). A rounded comparison.
  oa <- omega_squared(SS_effect = 2 * 54.061, df_effect = 2,
                      MS_error = 7.188, SS_total = 748.621, N = 66)
  expect_equal(round(oa$omega_squared_total, 3), 0.124)
  expect_equal(round(oa$omega_squared_partial, 3), 0.165)
  # Maxwell, Delaney, and Kelley, the blood pressure data (exact raw data,
  # Table 7.5): total .10, .34, .05 and partial .17, .40, .09.
  fit <- aov(y ~ A * B, data = .mdk_blood_pressure())
  res <- omega_squared(fit)
  ss <- c(1080, 3420, 780)
  expect_equal(res$omega_squared_total,
               .omega_total(ss, c(1, 2, 2), 155.75, 9018), tolerance = 1e-12)
  expect_equal(res$omega_squared_partial,
               .omega_partial(ss, c(1, 2, 2), 155.75, 30), tolerance = 1e-12)
  expect_equal(round(res$omega_squared_total, 2), c(0.10, 0.34, 0.05))
  expect_equal(round(res$omega_squared_partial, 2), c(0.17, 0.40, 0.09))
  # Eq. 7.47: the total from the F statistics of every effect.
  dfF <- c(1, 2, 2) * res$F_value
  expect_equal(res$omega_squared_total,
               c(1, 2, 2) * (res$F_value - 1) / (sum(dfF) + 24 + 1),
               tolerance = 1e-12)
})

test_that("the sums-of-squares interface takes exactly one error mean square and the inputs each estimand needs", {
  expect_error(omega_squared(SS_effect = 10, df_effect = 2, MS_error = 1,
                             SS_error = 30, df_error = 30, SS_total = 50),
               "two ways of\\s+giving the error mean square")
  expect_error(omega_squared(SS_effect = 10, df_effect = 2, SS_error = 30,
                             SS_total = 50), "give 'MS_error', or")
  expect_error(omega_squared(SS_effect = 10, MS_error = 1, SS_total = 50),
               "needs 'df_effect'")
  expect_error(omega_squared(SS_effect = 10, df_effect = 2, MS_error = 1),
               "Give 'SS_total'")
  expect_error(omega_squared(SS_effect = 10, df_effect = 2, MS_error = 1,
                             N = 2), "'N' must be greater than df_effect")
  # SS_total only: the partial value is NA, and refused by name.
  tot <- omega_squared(SS_effect = 10, df_effect = 2, MS_error = 1,
                       SS_total = 50)
  expect_true(is.na(tot$omega_squared_partial))
  expect_match(attr(tot, "estimand_notes")[["omega_squared_partial"]],
               "N was not given")
  expect_error(omega_squared(SS_effect = 10, df_effect = 2, MS_error = 1,
                             SS_total = 50, estimand = "partial"),
               class = "dmar_estimand_undetermined")
  expect_error(omega_squared(SS_effect = 10, df_effect = 2, MS_error = 1,
                             SS_error = 30, SS_total = 50), "both given")
  expect_error(omega_squared(SS_effect = 40, df_effect = 2, SS_error = 30,
                             df_error = 30, SS_total = 50),
               "smaller than SS_effect \\+ SS_error")
  expect_error(omega_squared(SS_effect = 10, df_effect = 2, MS_error = 1,
                             SS_total = 50, F_value = 3), "were both given")
  expect_error(omega_squared(F_value = 11.221, df_effect = 4,
                             df_error = 50), "'F_value' was given without 'N'")
  expect_error(omega_squared(F_value = -1, df_effect = 4, df_error = 50,
                             N = 55), "'F_value' must be")
  expect_error(omega_squared(object = list(not = "a model")),
               "aov, lm, or aovlist fit")
})


# ---------- Fitted models ----------

test_that("a one-way fit gives the same total and partial value, and matches the F interface", {
  fit <- aov(bdi_post ~ condition, data = depression_bdi)
  tbl <- anova(fit)
  res <- omega_squared(fit)
  ms_e <- tbl[["Mean Sq"]][2]
  expect_equal(res$omega_squared_total,
               .omega_total(tbl[["Sum Sq"]][1], 2, ms_e, sum(tbl[["Sum Sq"]])),
               tolerance = 1e-12)
  expect_equal(res$omega_squared_total, res$omega_squared_partial,
               tolerance = 1e-12)
  raw <- omega_squared(F_value = tbl[["F value"]][1], df_effect = 2,
                       df_error = 27, N = 30)
  expect_identical(raw$omega_squared_partial, res$omega_squared_partial)
  expect_identical(res$N, 30L)
})

test_that("the partial values of fits without an Error() term are those of DMAR 1.0.0", {
  # DMAR 1.0.0 (bcc0624), df (F - 1) / (df (F - 1) + N) with N = 310.
  pyg <- pygmalion
  pyg$grade <- factor(pyg$grade)
  one <- omega_squared(aov(iq_gain ~ grade, data = pyg), estimand = "partial")
  expect_named(one, c("effect", "omega_squared_partial", "F_value",
                      "df_effect", "df_error", "N"))
  expect_equal(one$omega_squared_partial, 0.10351136851412224,
               tolerance = 1e-10)
  fac <- omega_squared(aov(iq_8 ~ treatment * grade, data = pyg))
  expect_identical(fac$effect, c("treatment", "grade", "treatment:grade"))
  expect_equal(fac$omega_squared_partial,
               c(0.0175551142359686094, 0.0281424532197332052,
                 0.0030662003403285538), tolerance = 1e-10)
  expect_identical(fac$N, rep(310L, 3L))
  tbl <- anova(aov(iq_8 ~ treatment * grade, data = pyg))
  expect_equal(fac$omega_squared_total,
               .omega_total(tbl[["Sum Sq"]][1:3], tbl[["Df"]][1:3],
                            tbl[["Mean Sq"]][4], sum(tbl[["Sum Sq"]])),
               tolerance = 1e-12)
  expect_true(all(fac$omega_squared_total <= fac$omega_squared_partial))
})

test_that("omega_squared() and ci_omega_squared() agree on the point estimates row by row", {
  fit <- aov(y ~ A * B, data = .mdk_blood_pressure())
  pt <- omega_squared(fit)
  ci <- suppressWarnings(ci_omega_squared(fit))
  for (nm in names(pt)) expect_identical(ci[[nm]], pt[[nm]], info = nm)
})


# ---------- Fits with an Error() term ----------
#
# The partial value of an effect of an aovlist fit is computed with the
# size of its stratum, the size with which ci_eta_squared() converts the
# effect's limits; the N column reports the number of observations. DMAR
# 1.0.x computed every effect with the number of rows. The total is NA:
# the total omega squared is defined for between-subjects designs.

# The three IQ occasions of pygmalion in long form, as vignette("DMAR")
# analyzes them, with the expectancy treatment beside each child's rows.
.pygmalion_long <- function() {
  n <- nrow(pygmalion)
  occasions <- c("pretest", "intermediate", "end_of_study")
  data.frame(
    id        = factor(rep(seq_len(n), times = 3)),
    occasion  = factor(rep(occasions, each = n), levels = occasions),
    treatment = rep(pygmalion$treatment, times = 3),
    iq        = c(pygmalion$iq_pre, pygmalion$iq_4, pygmalion$iq_8)
  )
}

test_that("the pygmalion occasions are computed with the stratum's 620, inside the interval of their partial eta squared", {
  skip_on_cran()
  fit <- aov(iq ~ occasion + Error(id / occasion), data = .pygmalion_long())
  res <- omega_squared(fit)
  expect_named(res, c("effect", "omega_squared_total",
                      "omega_squared_partial", "stratum", "F_value",
                      "df_effect", "df_error", "N"))
  expect_identical(res$stratum, "id:occasion")
  expect_identical(res$N, 930)
  expect_true(is.na(res$omega_squared_total))
  expect_match(attr(res, "estimand_notes")[["omega_squared_total"]],
               "Error\\(\\) term")
  expect_error(omega_squared(fit, estimand = "total"),
               class = "dmar_estimand_undetermined")
  num <- 2 * (res$F_value - 1)
  expect_equal(res$omega_squared_partial, num / (num + 620),
               tolerance = 1e-12)
  expect_equal(round(res$omega_squared_partial, 3), 0.209)
  eta <- ci_eta_squared(fit, estimand = "partial")
  expect_gt(res$omega_squared_partial, eta$eta_squared_partial_ci_lower)
  expect_lt(res$omega_squared_partial, eta$eta_squared_partial_ci_upper)
  # The row's F, degrees of freedom, and the stratum's size reproduce the
  # estimate through the F interface.
  raw <- omega_squared(F_value = res$F_value, df_effect = 2, df_error = 618,
                       N = 620)
  expect_identical(raw$omega_squared_partial, res$omega_squared_partial)
  # With the 930 rows, the value of DMAR 1.0.x and of the partial omega
  # squared of Olejnik and Algina (2003) for a repeated measures factor,
  # the estimate lies below the interval.
  rows <- omega_squared(F_value = res$F_value, df_effect = 2, df_error = 618,
                        N = 930)
  expect_equal(round(rows$omega_squared_partial, 3), 0.150)
  expect_lt(rows$omega_squared_partial, eta$eta_squared_partial_ci_lower)
})

test_that("an effect in the subjects stratum of a mixed design is computed with the number of subjects", {
  skip_on_cran()
  fit <- aov(iq ~ treatment * occasion + Error(id / occasion),
             data = .pygmalion_long())
  res <- omega_squared(fit, estimand = "partial")
  expect_identical(res$effect, c("treatment", "occasion",
                                 "treatment:occasion"))
  expect_identical(res$stratum, c("id", "id:occasion", "id:occasion"))
  expect_identical(res$N, rep(930, 3))
  size <- c(310, 620, 620)
  num <- res$df_effect * (res$F_value - 1)
  expect_equal(res$omega_squared_partial, num / (num + size),
               tolerance = 1e-12)
  # Treatment, F(1, 308): .0039 with the 930 rows in 1.0.x, .0117 with the
  # 310 children.
  expect_equal(round(res$omega_squared_partial[1], 4), 0.0117)
  eta <- suppressWarnings(ci_eta_squared(fit, estimand = "partial"))
  for (k in seq_len(nrow(res))) {
    expect_gte(res$omega_squared_partial[k],
               eta$eta_squared_partial_ci_lower[k])
    expect_lte(res$omega_squared_partial[k],
               eta$eta_squared_partial_ci_upper[k])
  }
})


# ---------- tidy() and glance() ----------

test_that("tidy() and glance() of an omega squared table", {
  res <- omega_squared(aov(y ~ A * B, data = .mdk_blood_pressure()))
  td <- generics::tidy(res)
  expect_named(td, c("term", "effect", "estimate"))
  expect_identical(td$estimate, c(res$omega_squared_total,
                                  res$omega_squared_partial))
  gl <- generics::glance(res)
  expect_identical(gl$n_effects, 3L)
  expect_identical(gl$N, 30L)
  # Through 1.0.x the table was a plain data.frame, and tidy() had no
  # method for it.
  expect_s3_class(res, "dmar_tbl")
})
