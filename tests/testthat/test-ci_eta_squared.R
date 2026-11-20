test_that("point estimate matches the F-and-df formula (Bargmann 1970)", {
  res <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  expected <- 4 * 11.221 / (4 * 11.221 + 50)
  expect_equal(res$eta_squared_partial, expected, tolerance = 1e-12)
  expect_equal(res$effect, "overall")
})

test_that("returns the documented columns, the total NA from F", {
  res <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  expect_named(res, c("effect", "eta_squared_total",
                      "eta_squared_total_ci_lower",
                      "eta_squared_total_ci_upper", "eta_squared_partial",
                      "eta_squared_partial_ci_lower",
                      "eta_squared_partial_ci_upper",
                      "F_value", "df_effect", "df_error", "N"))
  expect_true(all(is.na(unlist(res[, 2:4]))))
  expect_identical(attr(res, "estimand"), "all")
  expect_named(attr(res, "estimand_notes"), "eta_squared_total")
  expect_s3_class(res, c("dmar_ci_anova", "dmar_tbl", "data.frame"),
                  exact = TRUE)
  part <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                         N = 55, estimand = "partial")
  expect_named(part, c("effect", "eta_squared_partial",
                       "eta_squared_partial_ci_lower",
                       "eta_squared_partial_ci_upper",
                       "F_value", "df_effect", "df_error", "N"))
  expect_identical(unclass(part)[[2]], unclass(res)$eta_squared_partial)
  expect_null(attr(part, "estimand_notes"))
})

test_that("raw CI matches ci_pvaf (shared Steiger transformation)", {
  res  <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  pvaf <- ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
  expect_equal(res$eta_squared_partial_ci_lower,
               pvaf$value[pvaf$term == "ci_lower"], tolerance = 1e-6)
  expect_equal(res$eta_squared_partial_ci_upper,
               pvaf$value[pvaf$term == "ci_upper"], tolerance = 1e-6)
})

test_that("CI brackets the point estimate", {
  res <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
  expect_gte(res$eta_squared_partial, res$eta_squared_partial_ci_lower)
  expect_lte(res$eta_squared_partial, res$eta_squared_partial_ci_upper)
})

test_that("CI widens with higher conf_level", {
  r90 <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
                           conf_level = 0.90)
  r99 <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
                           conf_level = 0.99)
  expect_gt(r99$eta_squared_partial_ci_upper - r99$eta_squared_partial_ci_lower,
            r90$eta_squared_partial_ci_upper - r90$eta_squared_partial_ci_lower)
})

test_that("low F yields ci_lower = 0", {
  expect_warning(
    res <- ci_eta_squared(F_value = 0.5, df_effect = 2, df_error = 30, N = 33),
    "below the alpha_lower critical value")
  expect_equal(res$eta_squared_partial_ci_lower, 0)
})

test_that("model-based one-way matches raw interface", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  tbl <- anova(fit)
  F_val <- tbl["group", "F value"]
  df_e  <- tbl["group", "Df"]
  df_r  <- tbl["Residuals", "Df"]
  N_val <- nobs(fit)

  res_model <- ci_eta_squared(fit)
  res_raw   <- ci_eta_squared(F_value = F_val, df_effect = df_e,
                                 df_error = df_r, N = N_val)

  expect_equal(res_model$eta_squared_partial, res_raw$eta_squared_partial, tolerance = 1e-12)
  expect_equal(res_model$eta_squared_partial_ci_lower, res_raw$eta_squared_partial_ci_lower, tolerance = 1e-6)
  expect_equal(res_model$eta_squared_partial_ci_upper, res_raw$eta_squared_partial_ci_upper, tolerance = 1e-6)
  expect_equal(res_model$effect, "group")
})

test_that("model-based factorial returns one row per effect", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  # The wool effect's low F legitimately clamps its lower limit to 0.
  expect_warning(res <- ci_eta_squared(fit),
                 "below the alpha_lower critical value")
  expect_equal(nrow(res), 3L)
  expect_setequal(res$effect, c("wool", "tension", "wool:tension"))
})

test_that("CI machinery agrees with ci_omega_squared (shared transformation)", {
  # Both functions use the same Steiger ncp -> [bound]/(bound + N) mapping;
  # only the point estimators differ.
  res_eta   <- ci_eta_squared(F_value = 11.221, df_effect = 4,
                                 df_error = 50, N = 55)
  res_omega <- ci_omega_squared(F_value = 11.221, df_effect = 4,
                                df_error = 50, N = 55)
  expect_equal(res_eta$eta_squared_partial_ci_lower,
               res_omega$omega_squared_partial_ci_lower, tolerance = 1e-9)
  expect_equal(res_eta$eta_squared_partial_ci_upper,
               res_omega$omega_squared_partial_ci_upper, tolerance = 1e-9)
})

test_that("supports an lm() fit as well as aov()", {
  fit_lm  <- lm(weight ~ group,  data = PlantGrowth)
  fit_aov <- aov(weight ~ group, data = PlantGrowth)
  res_lm  <- ci_eta_squared(fit_lm)
  res_aov <- ci_eta_squared(fit_aov)
  expect_equal(res_lm$eta_squared_partial,  res_aov$eta_squared_partial,  tolerance = 1e-12)
  expect_equal(res_lm$eta_squared_partial_ci_lower,  res_aov$eta_squared_partial_ci_lower,  tolerance = 1e-9)
  expect_equal(res_lm$eta_squared_partial_ci_upper,  res_aov$eta_squared_partial_ci_upper,  tolerance = 1e-9)
})

test_that("asymmetric alpha arguments are respected", {
  res_sym  <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                                N = 55, conf_level = 0.90)
  res_asym <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                                N = 55, alpha_lower = 0, alpha_upper = 0.10)
  expect_equal(res_asym$eta_squared_partial_ci_lower, 0)
  expect_lt(res_asym$eta_squared_partial_ci_upper, res_sym$eta_squared_partial_ci_upper)
})

test_that("errors on missing or invalid arguments", {
  expect_error(ci_eta_squared(), "Give a fitted aov, lm, or aovlist model",
               class = "dmar_argument_invalid")
  expect_error(ci_eta_squared(F_value = 1.5, df_effect = 2),
               "'F_value' was given without 'df_error' and 'N'")
  expect_error(ci_eta_squared(estimand = "both"), "'estimand' must be one of",
               class = "dmar_argument_invalid")
  expect_error(ci_eta_squared(aov(weight ~ group, data = PlantGrowth),
                              F_value = 2),
               "were both given", class = "dmar_argument_invalid")
  expect_error(ci_eta_squared(F_value = 3, df_effect = 2, df_error = 10,
                              N = 20, SS_effect = 4), "unused argument")
  expect_error(ci_eta_squared(F_value = -1, df_effect = 2,
                                 df_error = 10, N = 20),
               "'F_value' must be")
  expect_error(ci_eta_squared(F_value = 3, df_effect = 2,
                                 df_error = 10, N = 5),
               "at least df_effect \\+ df_error")
  expect_error(ci_eta_squared(object = list(foo = 1)), "aov, lm, or aovlist")
})

test_that("ci_eta_squared() reports the lower-limit clamp once, in terms of eta squared, and still returns the interval", {
  msgs <- capture_warnings(
    result <- ci_eta_squared(F_value = 1.2, df_effect = 4, df_error = 50, N = 55)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on partial eta squared is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_equal(result$eta_squared_partial_ci_lower, 0)
  expect_gt(result$eta_squared_partial_ci_upper, 0)
})

test_that("ci_eta_squared() on a model with several clamped effects warns once, not once per effect", {
  set.seed(113)
  d <- data.frame(
    y = rnorm(60),
    a = factor(rep(1:3, each = 20)),
    b = factor(rep(1:2, 30))
  )
  fit <- aov(y ~ a + b, data = d)
  msgs <- capture_warnings(result <- ci_eta_squared(fit))
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on partial eta squared is 0")
  expect_equal(result$eta_squared_partial_ci_lower, c(0, 0))
})

test_that("ci_eta_squared() refuses an F statistic that is not finite", {
  # Through the first builds of 1.1.0 an infinite F reached the noncentral
  # root search, which failed; it is refused on entry since the summary
  # interfaces were added.
  expect_error(
    ci_eta_squared(F_value = Inf, df_effect = 4, df_error = 50, N = 55),
    "'F_value' must be one or more finite non-negative numbers",
    class = "dmar_argument_invalid"
  )
})

test_that("ci_eta_squared() records the confidence level of its intervals", {
  # The 1.1.0 review: the table recorded no level, so the console printed
  # none and results_sentence() wrote a 90% interval as "95% CI" (the true
  # 95% interval for the first effect here is [0.053, 0.377]).
  fit <- aov(len ~ supp * factor(dose), data = ToothGrowth)
  e90 <- ci_eta_squared(fit, conf_level = 0.90)
  expect_identical(attr(e90, "conf_level"), 0.90)
  expect_s3_class(e90, c("dmar_ci_anova", "dmar_tbl", "data.frame"),
                  exact = TRUE)
  # The row holds no interval for the total, so the sentence is the
  # partial value's.
  expect_identical(results_sentence(e90[1, ], digits = 3),
                   "eta_squared_partial = 0.224, 90% CI [0.073, 0.350]")
  expect_identical(results_sentence(ci_eta_squared(fit)[1, ], digits = 3),
                   "eta_squared_partial = 0.224, 95% CI [0.053, 0.377]")
  expect_true("Confidence level: 90%" %in% utils::capture.output(print(e90)))
  expect_identical(generics::tidy(e90)$conf_level, rep(0.90, 6L))
  # The raw interface of the help page, and the level the tail areas
  # define.
  raw <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                        N = 55, conf_level = 0.90)
  expect_identical(attr(raw, "conf_level"), 0.90)
  tails <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                          N = 55, alpha_lower = 0.05, alpha_upper = 0.05)
  expect_equal(attr(tails, "conf_level"), 0.90, tolerance = 1e-12)
  expect_equal(tails$eta_squared_partial_ci_lower, raw$eta_squared_partial_ci_lower, tolerance = 1e-12)
  # A within-subjects fit records it too.
  set.seed(113)
  rm <- data.frame(subject = factor(rep(1:12, each = 3)),
                   time = factor(rep(1:3, 12)), y = rnorm(36))
  rm$y <- rm$y + as.numeric(rm$time)
  w <- ci_eta_squared(aov(y ~ time + Error(subject / time), data = rm),
                      conf_level = 0.90)
  expect_identical(attr(w, "conf_level"), 0.90)
})


# ---------- The size with which the limits are converted (1.1.0) ----------
#
# A fit without an Error() term and the raw interface convert the
# noncentrality limits with the total sample size, lambda / (lambda + N), as
# DMAR 1.0.x did; an aovlist fit converts each effect with its stratum's size
# (test-eta_squared_aovlist.R and test-ci_eta_squared_partial.R). The checks
# below invert the noncentral F by uniroot() on pf(), independently of
# ci_nc_F().

.ncp_limits_by_pf <- function(F_value, df_1, df_2, conf_level = 0.95) {
  a <- (1 - conf_level) / 2
  gap <- function(ncp, p) stats::pf(F_value, df_1, df_2, ncp = ncp) - p
  lower <- if (stats::pf(F_value, df_1, df_2) < 1 - a) 0 else
    stats::uniroot(gap, c(0, 1e4), p = 1 - a, tol = 1e-12)$root
  upper <- stats::uniroot(gap, c(0, 1e4), p = a, tol = 1e-12)$root
  c(lower, upper)
}

test_that("a fit without an Error() term converts with N, so one-way and factorial limits are those of 1.0.x", {
  one <- ci_eta_squared(aov(iq_gain ~ factor(grade), data = pygmalion))
  # DMAR 1.0.x (e093b1af).
  expect_equal(c(one$eta_squared_partial_ci_lower, one$eta_squared_partial_ci_upper),
               c(0.0480719904184997, 0.175788966752462), tolerance = 1e-10)
  L <- .ncp_limits_by_pf(one$F_value, one$df_effect, one$df_error)
  expect_equal(c(one$eta_squared_partial_ci_lower, one$eta_squared_partial_ci_upper), L / (L + 310),
               tolerance = 1e-6)
  fit <- aov(iq_gain ~ treatment * factor(grade), data = pygmalion)
  res <- suppressWarnings(ci_eta_squared(fit))
  expect_identical(res$N, rep(310L, 3L))
  expect_equal(res$eta_squared_partial_ci_upper,
               c(0.0498375103356816, 0.180634872665218, 0.0608023036598636),
               tolerance = 1e-10)
  for (k in seq_len(nrow(res))) {
    L <- .ncp_limits_by_pf(res$F_value[k], res$df_effect[k], res$df_error[k])
    expect_equal(c(res$eta_squared_partial_ci_lower[k], res$eta_squared_partial_ci_upper[k]), L / (L + 310),
                 tolerance = 1e-6, info = res$effect[k])
    expect_gte(res$eta_squared_partial[k], res$eta_squared_partial_ci_lower[k])
    expect_lte(res$eta_squared_partial[k], res$eta_squared_partial_ci_upper[k])
    raw <- suppressWarnings(ci_eta_squared(F_value = res$F_value[k],
                                           df_effect = res$df_effect[k],
                                           df_error = res$df_error[k],
                                           N = 310))
    expect_identical(raw$eta_squared_partial_ci_upper, res$eta_squared_partial_ci_upper[k])
  }
  # estimand = "partial", factorial_anova(), and ci_omega_squared()
  # convert with the same N, so they share the interval.
  part <- suppressWarnings(ci_eta_squared(fit, estimand = "partial"))
  expect_identical(part$eta_squared_partial_ci_upper, res$eta_squared_partial_ci_upper)
  om <- suppressWarnings(ci_omega_squared(fit))
  expect_equal(om$omega_squared_partial_ci_upper,
               res$eta_squared_partial_ci_upper, tolerance = 1e-12)
})

test_that("the total sample size and Smithson's degrees of freedom conversion reproduce Steiger (2004), Smithson (2001), and Fidler and Thompson (2001)", {
  # The 2 by 7 design with four observations per cell (N = 56), 90%.
  # Steiger (2004, p. 171): interaction, F = 2.50 on 6 and 42 df,
  # [.0082, .3160] with N = 56. Smithson (2001, p. 619) converts with
  # u + v + 1: interaction [.0093, .3455] from a lower noncentrality limit
  # of 0.4621 (ci_nc_F() and Steiger give 0.4628, hence .0094 here), row
  # effect (F = 6.00 on 1 and 42 df) [.0117, .2801], column effect
  # (F = 2.65 on 6 and 42 df) [.0174, .3577].
  ab_56 <- ci_eta_squared(F_value = 2.50, df_effect = 6, df_error = 42,
                          N = 56, conf_level = 0.90)
  expect_equal(round(c(ab_56$eta_squared_partial_ci_lower, ab_56$eta_squared_partial_ci_upper), 4),
               c(0.0082, 0.3160))
  ab_49 <- ci_eta_squared(F_value = 2.50, df_effect = 6, df_error = 42,
                          N = 49, conf_level = 0.90)
  expect_equal(round(ab_49$eta_squared_partial_ci_upper, 4), 0.3455)
  expect_equal(round(ab_49$eta_squared_partial_ci_lower, 3), 0.009)
  a_44 <- ci_eta_squared(F_value = 6.00, df_effect = 1, df_error = 42,
                         N = 44, conf_level = 0.90)
  expect_equal(round(c(a_44$eta_squared_partial_ci_lower, a_44$eta_squared_partial_ci_upper), 4),
               c(0.0117, 0.2801))
  b_49 <- ci_eta_squared(F_value = 2.65, df_effect = 6, df_error = 42,
                         N = 49, conf_level = 0.90)
  expect_equal(round(c(b_49$eta_squared_partial_ci_lower, b_49$eta_squared_partial_ci_upper), 4),
               c(0.0174, 0.3577))
  # The values the help page quotes for the two conversions.
  a_56 <- ci_eta_squared(F_value = 6.00, df_effect = 1, df_error = 42,
                         N = 56, conf_level = 0.90)
  expect_equal(round(c(a_56$eta_squared_partial_ci_lower, a_56$eta_squared_partial_ci_upper,
                       a_44$eta_squared_partial_ci_lower, a_44$eta_squared_partial_ci_upper), 3),
               c(0.009, 0.234, 0.012, 0.280))
  expect_equal(round(c(ab_56$eta_squared_partial_ci_lower, ab_56$eta_squared_partial_ci_upper,
                       ab_49$eta_squared_partial_ci_lower, ab_49$eta_squared_partial_ci_upper), 3),
               c(0.008, 0.316, 0.009, 0.346))
  # Fidler and Thompson (2001, Table 8, p. 596): the noncentral intervals
  # of their 2 by 4 fixed-effects ANOVA (N = 24, 16 error df), computed
  # with Smithson's syntax from F and the df alone, are the 90% limits
  # converted with u + v + 1.
  ft <- function(F_value, df_effect)
    suppressWarnings(ci_eta_squared(F_value = F_value, df_effect = df_effect,
                                    df_error = 16, N = df_effect + 17,
                                    conf_level = 0.90))
  a <- ft(1.5, 1); b <- ft(4.0, 3); ab <- ft(1.5, 3)
  got <- c(a$eta_squared_partial_ci_lower, a$eta_squared_partial_ci_upper, b$eta_squared_partial_ci_lower, b$eta_squared_partial_ci_upper,
           ab$eta_squared_partial_ci_lower, ab$eta_squared_partial_ci_upper)
  expect_true(all(abs(got - c(0, .31669, .03572, .56709, 0, .37779)) < 5e-5))
})

test_that("the coverage the help page quotes for Smithson's conversion of a factorial effect", {
  # Converted with u + v + 1 = 44 instead of N = 56, the 95% interval for
  # the row effect of the 2 by 7 design covers the population partial eta
  # squared lambda / (lambda + 56) less often than 95%. Its coverage is the
  # probability that F falls between the .025 and .975 quantiles of the
  # noncentral F at the noncentrality (44 / 56) lambda.
  coverage <- function(eta2, M, D = 56, d1 = 1, d2 = 42) {
    lambda <- D * eta2 / (1 - eta2)
    lp <- (M / D) * lambda
    stats::pf(stats::qf(0.975, d1, d2, ncp = lp), d1, d2, ncp = lambda) -
      stats::pf(stats::qf(0.025, d1, d2, ncp = lp), d1, d2, ncp = lambda)
  }
  expect_equal(coverage(0.20, M = 56), 0.95, tolerance = 1e-8)
  expect_equal(round(c(coverage(0.02, M = 44), coverage(0.40, M = 44)), 3),
               c(0.945, 0.919))
})

test_that("the raw interface converts with N as given and accepts the size of a one-effect within-subjects stratum", {
  # Occasion in the three-occasion pygmalion analysis: F(2, 618), 310
  # children, 930 rows, a stratum of 2 + 618 = 620 degrees of freedom.
  L <- .ncp_limits_by_pf(83.11, 2, 618)
  stratum <- ci_eta_squared(estimand = "partial", F_value = 83.11, df_effect = 2,
                                    df_error = 618, N = 620)
  expect_equal(c(stratum$eta_squared_partial_ci_lower, stratum$eta_squared_partial_ci_upper), L / (L + 620),
               tolerance = 1e-6)
  expect_identical(stratum$N, 620)
  rows <- ci_eta_squared(estimand = "partial", F_value = 83.11, df_effect = 2,
                                 df_error = 618, N = 930)
  expect_equal(c(rows$eta_squared_partial_ci_lower, rows$eta_squared_partial_ci_upper), L / (L + 930),
               tolerance = 1e-6)
  expect_error(ci_eta_squared(F_value = 3, df_effect = 2, df_error = 10,
                              N = 11), "at least df_effect \\+ df_error")
  expect_error(ci_eta_squared(F_value = 3, df_effect = 2, df_error = 10,
                              N = NA), "at least df_effect \\+ df_error")
  expect_error(ci_eta_squared(F_value = 3, df_effect = 2, df_error = 10),
               "'F_value' was given without 'N'")
  expect_error(ci_eta_squared(estimand = "partial", F_value = 3,
                              df_effect = 2, df_error = 10),
               "'F_value' was given without 'N'")
})

# ---------- estimand = "partial" (ci_eta_squared_partial() through 1.0.x) ----------
#
# The tests of ci_eta_squared_partial(), removed in 1.1.0, folded in: its
# interval is ci_eta_squared(estimand = "partial").

test_that("estimand = \"partial\" refuses missing or invalid arguments", {
  expect_error(ci_eta_squared(estimand = "partial"), "Give a fitted aov")
  expect_error(ci_eta_squared(estimand = "partial", F_value = -1, df_effect = 2,
                                         df_error = 10, N = 20),
               "'F_value' must be")
})

test_that("estimand = \"partial\" reports the lower-limit clamp once, in terms of partial eta squared, and still returns the interval", {
  msgs <- capture_warnings(
    result <- ci_eta_squared(estimand = "partial", F_value = 1.2, df_effect = 4, df_error = 50, N = 55)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on partial eta squared is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_equal(result$eta_squared_partial_ci_lower, 0)
  expect_gt(result$eta_squared_partial_ci_upper, 0)
})

test_that("estimand = \"partial\" records the confidence level of its intervals", {
  # The 1.1.0 review: the table recorded no level, so results_sentence()
  # wrote a 90% interval as "95% CI".
  fit <- aov(len ~ supp * factor(dose), data = ToothGrowth)
  p90 <- ci_eta_squared(estimand = "partial", fit, conf_level = 0.90)
  expect_identical(attr(p90, "conf_level"), 0.90)
  expect_match(results_sentence(p90[1, ], digits = 3),
               "^eta_squared_partial = 0\\.224, 90% CI \\[0\\.073, 0\\.350\\]$")
  expect_true("Confidence level: 90%" %in% utils::capture.output(print(p90)))
  raw <- ci_eta_squared(estimand = "partial", F_value = 11.221, df_effect = 4,
                                df_error = 50, N = 55, conf_level = 0.90)
  expect_identical(attr(raw, "conf_level"), 0.90)
})

test_that("the three-occasion pygmalion interval brackets its estimate and is the pf() inversion converted with the stratum's 620 degrees of freedom", {
  # The within-subjects analysis of vignette("DMAR"): 310 children, three
  # IQ occasions, F(2, 618). DMAR 1.0.x converted the noncentrality limits
  # with the 930 rows and put the interval, [.111, .194], below the
  # estimate of .212.
  n <- nrow(pygmalion)
  occasions <- c("pretest", "intermediate", "end_of_study")
  iq_long <- data.frame(
    id       = factor(rep(seq_len(n), times = 3)),
    occasion = factor(rep(occasions, each = n), levels = occasions),
    iq       = c(pygmalion$iq_pre, pygmalion$iq_4, pygmalion$iq_8)
  )
  fit <- aov(iq ~ occasion + Error(id / occasion), data = iq_long)
  res <- ci_eta_squared(estimand = "partial", fit)
  expect_identical(res$effect, "occasion")
  expect_identical(res$stratum, "id:occasion")
  expect_equal(c(res$df_effect, res$df_error), c(2, 618))
  expect_identical(res$N, 930)
  expect_gt(res$eta_squared_partial, res$eta_squared_partial_ci_lower)
  expect_lt(res$eta_squared_partial, res$eta_squared_partial_ci_upper)

  # Independently: the noncentral F inverted by uniroot() on pf(), and the
  # limits converted with the stratum's 2 + 618 = 620 degrees of freedom.
  a <- 0.025
  gap <- function(ncp, p) stats::pf(res$F_value, 2, 618, ncp = ncp) - p
  L <- c(stats::uniroot(gap, c(0, 1e4), p = 1 - a, tol = 1e-12)$root,
         stats::uniroot(gap, c(0, 1e4), p = a, tol = 1e-12)$root)
  expect_equal(c(res$eta_squared_partial_ci_lower, res$eta_squared_partial_ci_upper), L / (L + 620),
               tolerance = 1e-6)
  # The values the help page, NEWS, and vignette("DMAR") quote.
  expect_equal(round(c(res$eta_squared_partial, res$eta_squared_partial_ci_lower,
                       res$eta_squared_partial_ci_upper), 3), c(0.212, 0.157, 0.265))
  # The raw interface with the stratum's size reproduces the row.
  raw <- ci_eta_squared(estimand = "partial", F_value = res$F_value, df_effect = 2,
                                df_error = 618, N = 620)
  expect_identical(c(raw$eta_squared_partial_ci_lower, raw$eta_squared_partial_ci_upper),
                   c(res$eta_squared_partial_ci_lower, res$eta_squared_partial_ci_upper))

  # DMAR 1.0.x (e093b1af): the same estimate, the interval below it; the
  # raw interface with the 930 rows reproduces that interval.
  expect_equal(res$eta_squared_partial, 0.211961634219325, tolerance = 1e-10)
  old <- ci_eta_squared(estimand = "partial", F_value = res$F_value, df_effect = 2,
                                df_error = 618, N = 930)
  expect_equal(c(old$eta_squared_partial_ci_lower, old$eta_squared_partial_ci_upper),
               c(0.110669058375401, 0.193552508221394), tolerance = 1e-10)
  expect_gt(res$eta_squared_partial, old$eta_squared_partial_ci_upper)

  # ci_eta_squared() gives the same row.
  tot <- ci_eta_squared(fit)
  expect_identical(tot$eta_squared_partial_ci_lower, res$eta_squared_partial_ci_lower)
  expect_identical(tot$eta_squared_partial_ci_upper, res$eta_squared_partial_ci_upper)

  # With the expectancy treatment as a between-subjects factor, its effect
  # is tested among the 310 children and converted with them; NEWS quotes
  # its upper limit, 0.018 when converted with the 930 rows in 1.0.x.
  iq_long$treatment <- rep(pygmalion$treatment, times = 3)
  mixed <- suppressWarnings(ci_eta_squared(estimand = "partial", 
    aov(iq ~ treatment * occasion + Error(id / occasion), data = iq_long)))
  tr <- mixed[mixed$effect == "treatment", ]
  expect_identical(tr$stratum, "id")
  expect_equal(round(tr$eta_squared_partial_ci_upper, 3), 0.052)
  tr_rows <- suppressWarnings(ci_eta_squared(estimand = "partial", 
    F_value = tr$F_value, df_effect = 1, df_error = tr$df_error, N = 930))
  expect_equal(round(tr_rows$eta_squared_partial_ci_upper, 3), 0.018)
  expect_equal(mixed$eta_squared_partial_ci_upper[mixed$effect == "occasion"],
               suppressWarnings(ci_eta_squared(estimand = "partial", 
                 F_value = mixed$F_value[mixed$effect == "occasion"],
                 df_effect = 2, df_error = 616, N = 620))$eta_squared_partial_ci_upper,
               tolerance = 1e-12)
})

test_that("the within-subjects example of the help page gives the interval NEWS quotes", {
  set.seed(113)
  n <- 20
  rm_data <- data.frame(
    subject = factor(rep(seq_len(n), each = 3)),
    time    = factor(rep(c("Pre", "Mid", "Post"), n),
                     levels = c("Pre", "Mid", "Post")),
    y       = rnorm(n, sd = 1.5)[rep(seq_len(n), each = 3)] +
              0.7 * rep(1:3, n) + rnorm(n * 3, sd = 1.2)
  )
  res <- ci_eta_squared(estimand = "partial", aov(y ~ time + Error(subject/time),
                                    data = rm_data))
  # Converted with the 2 + 38 = 40 degrees of freedom of the subject:time
  # stratum; DMAR 1.0.x converted with the 60 rows.
  expect_equal(round(c(res$eta_squared_partial_ci_lower, res$eta_squared_partial_ci_upper), 3),
               c(0.023, 0.417))
  old <- ci_eta_squared(estimand = "partial", F_value = res$F_value, df_effect = 2,
                                df_error = 38, N = 60)
  expect_equal(round(c(old$eta_squared_partial_ci_lower, old$eta_squared_partial_ci_upper), 3),
               c(0.015, 0.322))
})


# ---------- The total and the partial estimand (1.1.0) ----------
#
# The noncentral F interval from an effect's F is an interval for its
# partial value. The total has an interval only in a one-way design, where
# the two coincide; elsewhere its limits are NA under estimand = "all" and
# estimand = "total" is refused.

test_that("a one-way fit reports the total and the partial value with the same interval", {
  fit <- aov(iq_gain ~ factor(grade), data = pygmalion)
  res <- ci_eta_squared(fit)
  ss <- anova(fit)[["Sum Sq"]]
  expect_equal(res$eta_squared_total, ss[1] / sum(ss), tolerance = 1e-12)
  expect_equal(res$eta_squared_total, res$eta_squared_partial,
               tolerance = 1e-12)
  expect_identical(res$eta_squared_total_ci_lower,
                   res$eta_squared_partial_ci_lower)
  expect_identical(res$eta_squared_total_ci_upper,
                   res$eta_squared_partial_ci_upper)
  L <- .ncp_by_pf(res$F_value, res$df_effect, res$df_error)
  expect_equal(c(res$eta_squared_total_ci_lower,
                 res$eta_squared_total_ci_upper), L / (L + 310),
               tolerance = 1e-6)
  expect_null(attr(res, "estimand_notes"))
  tot <- ci_eta_squared(fit, estimand = "total")
  expect_named(tot, c("effect", "eta_squared_total",
                      "eta_squared_total_ci_lower",
                      "eta_squared_total_ci_upper", "F_value", "df_effect",
                      "df_error", "N"))
  expect_identical(attr(tot, "estimand"), "total")
  expect_identical(tot$eta_squared_total_ci_upper,
                   res$eta_squared_total_ci_upper)
  # Two intervals in one row: the sentence names the one it reports.
  expect_identical(results_sentence(res, term = "eta_squared_total"),
                   results_sentence(tot))
  expect_error(results_sentence(res), "intervals of several estimates",
               class = "dmar_argument_invalid")
})

test_that("a factorial fit reports the total estimate without limits and refuses estimand = \"total\"", {
  fit <- aov(y ~ A * B, data = .mdk_blood_pressure())
  res <- suppressWarnings(ci_eta_squared(fit))
  # Maxwell, Delaney, and Kelley (2027), Table 7.7: SS 1,080, 3,420, 780,
  # and 3,738 within cells, a corrected total of 9,018.
  expect_equal(res$eta_squared_total, c(1080, 3420, 780) / 9018,
               tolerance = 1e-12)
  expect_equal(res$eta_squared_partial,
               c(1080, 3420, 780) / (c(1080, 3420, 780) + 3738),
               tolerance = 1e-12)
  expect_true(all(is.na(res$eta_squared_total_ci_lower)))
  expect_true(all(is.na(res$eta_squared_total_ci_upper)))
  for (k in 1:3) {
    L <- .ncp_by_pf(res$F_value[k], res$df_effect[k], 24)
    expect_equal(c(res$eta_squared_partial_ci_lower[k],
                   res$eta_squared_partial_ci_upper[k]), L / (L + 30),
                 tolerance = 1e-6, info = res$effect[k])
  }
  expect_named(attr(res, "estimand_notes"), "eta_squared_total_ci")
  expect_true(any(grepl("^NA eta_squared_total limits: No confidence interval",
                        utils::capture.output(print(res)))))
  # Refused before any interval is computed, so no clamp warning (the AB
  # interaction's lower limit is clamped) comes with the refusal.
  seen <- new.env(parent = emptyenv())
  seen$warnings <- character(0)
  e <- withCallingHandlers(
    tryCatch(ci_eta_squared(fit, estimand = "total"), error = function(e) e),
    warning = function(w_) {
      seen$warnings <- c(seen$warnings, conditionMessage(w_))
      invokeRestart("muffleWarning")
    })
  expect_s3_class(e, c("dmar_estimand_undetermined", "dmar_refusal"),
                  exact = FALSE)
  expect_identical(e$estimand, "total")
  expect_match(e$reason, "interval for the partial value")
  expect_length(seen$warnings, 0L)
  part <- suppressWarnings(ci_eta_squared(fit, estimand = "partial"))
  expect_identical(names(part),
                   setdiff(names(res), grep("_total", names(res),
                                            value = TRUE)))
  expect_identical(part$eta_squared_partial_ci_upper,
                   res$eta_squared_partial_ci_upper)
})

test_that("the F interface determines neither the total nor its interval", {
  res <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                        N = 55)
  expect_true(all(is.na(c(res$eta_squared_total,
                          res$eta_squared_total_ci_lower,
                          res$eta_squared_total_ci_upper))))
  expect_match(attr(res, "estimand_notes")[["eta_squared_total"]],
               "use estimand = \"partial\"")
  expect_error(ci_eta_squared(F_value = 11.221, df_effect = 4,
                              df_error = 50, N = 55, estimand = "total"),
               class = "dmar_estimand_undetermined")
  # The partial interval is the one of 1.0.x, Steiger's (2004) 90% limits
  # for Bargmann's example being .261 and .565.
  r90 <- ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                        N = 55, estimand = "partial", conf_level = 0.90)
  expect_equal(round(c(r90$eta_squared_partial_ci_lower,
                       r90$eta_squared_partial_ci_upper), 3),
               c(0.261, 0.565))
})

test_that("several F statistics give one row each, as separate calls do", {
  Fs <- c(A = 6.934, B = 10.979, "A:B" = 2.504)
  res <- suppressWarnings(ci_eta_squared(F_value = Fs, df_effect = c(1, 2, 2),
                                         df_error = 24, N = 30,
                                         estimand = "partial"))
  expect_identical(res$effect, names(Fs))
  for (k in 1:3) {
    one <- suppressWarnings(ci_eta_squared(F_value = Fs[[k]],
                                           df_effect = c(1, 2, 2)[k],
                                           df_error = 24, N = 30,
                                           estimand = "partial"))
    expect_identical(unlist(one[1, -1]), unlist(res[k, -1]),
                     info = names(Fs)[k])
  }
  expect_identical(suppressWarnings(ci_eta_squared(
    F_value = c(3, 4), df_effect = 2, df_error = 20, N = 23))$effect,
    c("effect_1", "effect_2"))
  expect_error(ci_eta_squared(F_value = c(3, 4), df_effect = c(1, 2, 3),
                              df_error = 20, N = 30),
               "one value for all effects or one for each")
})

test_that("an aovlist fit reports the total over every stratum, without limits", {
  set.seed(113)
  d <- data.frame(id = factor(rep(1:15, each = 3)),
                  time = factor(rep(1:3, 15)))
  d$y <- rnorm(15, sd = 1.5)[as.integer(d$id)] + 0.5 * as.integer(d$time) +
    rnorm(45)
  fit <- aov(y ~ time + Error(id / time), data = d)
  res <- ci_eta_squared(fit)
  ss_time <- summary(fit)[["Error: id:time"]][[1]]["time", "Sum Sq"]
  expect_equal(res$eta_squared_total, ss_time / sum((d$y - mean(d$y))^2),
               tolerance = 1e-12)
  expect_lt(res$eta_squared_total, res$eta_squared_partial)
  expect_true(is.na(res$eta_squared_total_ci_lower))
  expect_identical(res$stratum, "id:time")
  expect_error(ci_eta_squared(fit, estimand = "total"),
               class = "dmar_estimand_undetermined")
})

test_that("tidy() and glance() give one row per effect and estimand", {
  fit <- aov(y ~ A * B, data = .mdk_blood_pressure())
  res <- suppressWarnings(ci_eta_squared(fit))
  td <- generics::tidy(res)
  expect_named(td, c("term", "effect", "estimate", "ci_lower", "ci_upper",
                     "conf_level"))
  expect_identical(td$term, rep(c("eta_squared_total",
                                  "eta_squared_partial"), each = 3))
  expect_identical(td$estimate,
                   c(res$eta_squared_total, res$eta_squared_partial))
  expect_identical(td$ci_upper[4:6], res$eta_squared_partial_ci_upper)
  expect_true(all(is.na(td$ci_upper[1:3])))
  expect_identical(results_sentence(td[5, ]), results_sentence(res[2, ]))
  gl <- generics::glance(res)
  expect_identical(gl$term, "eta_squared_total, eta_squared_partial")
  expect_identical(gl$n_effects, 3L)
  expect_identical(gl$N, 30L)
  p1 <- suppressWarnings(ci_eta_squared(fit, estimand = "partial"))[2, ]
  expect_identical(generics::glance(p1), generics::tidy(p1))
})
