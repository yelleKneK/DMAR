test_that(".format_dmar_value keeps whole numbers integer and others at 3 sig figs", {
  fmt <- DMAR:::.format_dmar_value
  # Whole numbers (degrees of freedom, sample sizes) print without decimals.
  expect_identical(fmt(c(1, 27, 0, -3, 30, 250)),
                   c("1", "27", "0", "-3", "30", "250"))
  # Non-integers print to 3 significant figures by default.
  expect_identical(fmt(29.28272), "29.3")
  expect_identical(fmt(0.5202791), "0.52")
  expect_identical(fmt(0.2246837), "0.225")
  # Very small magnitudes use scientific notation; whole numbers never do.
  expect_identical(fmt(1.010542e-05), "1.01e-05")
  expect_false(grepl("e", fmt(27), fixed = TRUE))
  # Missing values are preserved as NA, not "NA" coerced through format.
  expect_true(is.na(fmt(NA_real_)))
})

test_that(".format_dmar_value honors the digits argument", {
  fmt <- DMAR:::.format_dmar_value
  expect_identical(fmt(29.28272, digits = 5L), "29.283")
  expect_identical(fmt(29.28272, digits = 2L), "29")
})

test_that(".format_dmar_pvalue uses fixed decimals with a < 0.0001 floor", {
  fmt <- DMAR:::.format_dmar_pvalue
  # Four decimal places by default, padding trailing zeros.
  expect_identical(fmt(c(0.0234, 0.3361016, 0.5, 1)),
                   c("0.0234", "0.3361", "0.5000", "1.0000"))
  # Below the 0.0001 floor (including exact zero) reads as "< 0.0001",
  # never scientific notation and never a rounded "0.0000".
  expect_identical(fmt(1.010542e-05), "< 0.0001")
  expect_identical(fmt(0), "< 0.0001")
  # The digits_p argument moves both the decimals and the floor together.
  expect_identical(fmt(0.3361016, digits_p = 2L), "0.34")
  expect_identical(fmt(1e-04, digits_p = 2L), "< 0.01")
  # Missing values are preserved.
  expect_true(is.na(fmt(NA_real_)))
})

test_that(".format_dmar_fixed uses fixed decimals with no floor or scientific fallback", {
  fmt <- DMAR:::.format_dmar_fixed
  # Three decimal places by default, padding trailing zeros, so an AIC of
  # 2284.830 keeps the decimals that three significant figures would discard.
  expect_identical(fmt(c(2284.83, -1134.415, 0, 5)),
                   c("2284.830", "-1134.415", "0.000", "5.000"))
  # Unlike the p-value formatter there is no "< ..." floor: a tiny value
  # rounds to the decimals rather than reading as below a threshold.
  expect_identical(fmt(1e-05), "0.000")
  # The digits_fixed argument moves the number of decimals.
  expect_identical(fmt(2284.830284, digits_fixed = 2L), "2284.83")
  # Missing values are preserved.
  expect_true(is.na(fmt(NA_real_)))
})

test_that("ancova returns a dmar_tbl with full-precision numeric values", {
  set.seed(113)
  n   <- 30
  grp <- rep(c("A", "B"), each = n / 2)
  x   <- rnorm(n, 50, 10)
  y   <- 0.6 * x + 12 * (grp == "B") + rnorm(n, 0, 4)
  d   <- data.frame(y = y, group = grp, x = x)
  res <- ancova(d, outcome = "y", treatment = "group", covariates = "x")

  # Class is additive: still a data.frame, now also a dmar_tbl.
  expect_s3_class(res, "dmar_tbl")
  expect_s3_class(res, "data.frame")
  # Stored value column is numeric and keeps full precision (display only is
  # rounded). Degrees of freedom are stored as exact whole numbers.
  expect_type(res$value, "double")
  expect_identical(res$value[res$term == "df_1"], 1)
  expect_identical(res$value[res$term == "df_2"], 27)
  expect_gt(res$value[res$term == "F_value"], 50)  # not rounded in storage
  expect_false(res$value[res$term == "F_value"] == 57.6)
  # Confidence level recorded for the print footer.
  expect_identical(attr(res, "conf_level"), 0.95)
})

test_that("print.dmar_tbl shows integers for df and avoids scientific notation there", {
  set.seed(113)
  n   <- 30
  grp <- rep(c("A", "B"), each = n / 2)
  x   <- rnorm(n, 50, 10)
  y   <- 0.6 * x + 12 * (grp == "B") + rnorm(n, 0, 4)
  d   <- data.frame(y = y, group = grp, x = x)
  res <- ancova(d, outcome = "y", treatment = "group", covariates = "x")

  out <- capture.output(print(res))
  txt <- paste(out, collapse = "\n")

  # Degrees of freedom print as bare integers, never in scientific notation.
  expect_true(any(grepl("^\\s*df_1\\s+1\\s*$", out)))
  expect_true(any(grepl("^\\s*df_2\\s+27\\s*$", out)))
  expect_false(grepl("e\\+00", txt))
  # F to 3 significant figures, not full precision.
  expect_true(any(grepl("F_value\\s+57\\.6\\b", out)))
  expect_false(grepl("57.566071", txt, fixed = TRUE))
  # p-values print to 4 decimal places; one below the 0.0001 floor reads as
  # "< 0.0001" rather than rounding to 0.0000, and the homogeneity-of-
  # regression p-value reads to four decimals.
  expect_true(any(grepl("p_value\\s+< 0\\.0001", out)))
  expect_true(any(grepl("p_homogeneity_of_regression\\s+0\\.3361", out)))
  # Confidence-level footer present.
  expect_true(any(grepl("Confidence level: 95%", out)))
})

test_that("print.dmar_tbl digits argument and dmar.digits option control precision", {
  set.seed(113)
  n   <- 30
  grp <- rep(c("A", "B"), each = n / 2)
  x   <- rnorm(n, 50, 10)
  y   <- 0.6 * x + 12 * (grp == "B") + rnorm(n, 0, 4)
  d   <- data.frame(y = y, group = grp, x = x)
  res <- ancova(d, outcome = "y", treatment = "group", covariates = "x")

  out5 <- capture.output(print(res, digits = 5))
  expect_true(any(grepl("F_value\\s+57\\.566\\b", out5)))

  old <- getOption("dmar.digits")
  on.exit(options(dmar.digits = old), add = TRUE)
  options(dmar.digits = 2)
  out2 <- capture.output(print(res))
  expect_true(any(grepl("F_value\\s+58\\b", out2)))
})

test_that("format.dmar_tbl formats a dedicated p_value column to 4 decimals", {
  # A wide-format dmar_tbl: p-values live in their own column (the ancova
  # path above is long-format and exercises the p_terms attribute instead).
  x <- data.frame(
    term    = c("a", "b", "c"),
    estimate = c(1.23456, 2, 0.25),
    p_value  = c(0.0234, 1.010542e-05, 0.5),
    stringsAsFactors = FALSE
  )
  class(x) <- c("dmar_tbl", "data.frame")
  fx <- format(x)
  expect_identical(fx$p_value, c("0.0234", "< 0.0001", "0.5000"))
  # The non-p_value numeric column still uses significant figures, and a
  # whole number there still prints without a decimal part.
  expect_identical(fx$estimate, c("1.23", "2", "0.25"))
})

test_that("format.dmar_tbl returns formatted character columns without mutating x", {
  set.seed(113)
  n   <- 30
  grp <- rep(c("A", "B"), each = n / 2)
  x   <- rnorm(n, 50, 10)
  y   <- 0.6 * x + 12 * (grp == "B") + rnorm(n, 0, 4)
  d   <- data.frame(y = y, group = grp, x = x)
  res <- ancova(d, outcome = "y", treatment = "group", covariates = "x")

  fres <- format(res)
  expect_type(fres$value, "character")
  expect_identical(fres$value[fres$term == "df_2"], "27")
  # The original object is untouched (still numeric, full precision).
  expect_type(res$value, "double")
})

# ---- Package-wide rollout of the dmar_tbl class --------------------------
# The pretty printer is opted into by ~115 tidy-returning functions through
# .as_dmar_tbl(). These tests guard the rollout invariants on a representative
# sample: the class is applied, the stored value column stays numeric and
# full-precision (only the display rounds), and any leading broom-dispatch
# subclass is preserved ahead of dmar_tbl so tidy() / glance() still work.

test_that(".as_dmar_tbl tags a bare table and is idempotent", {
  out <- data.frame(term = c("a", "df"), value = c(1.23456, 7),
                    stringsAsFactors = FALSE)
  tagged <- DMAR:::.as_dmar_tbl(out, conf_level = 0.95)
  expect_identical(class(tagged), c("dmar_tbl", "data.frame"))
  expect_identical(attr(tagged, "conf_level"), 0.95)
  # Calling again does not duplicate the class or reorder it.
  again <- DMAR:::.as_dmar_tbl(tagged)
  expect_identical(class(again), c("dmar_tbl", "data.frame"))
  # A leading subclass is kept ahead of dmar_tbl, before data.frame.
  sub <- out
  class(sub) <- c("dmar_ci_long", "data.frame")
  sub <- DMAR:::.as_dmar_tbl(sub)
  expect_identical(class(sub), c("dmar_ci_long", "dmar_tbl", "data.frame"))
})

test_that("rolled-out functions return dmar_tbl with full-precision numeric values", {
  set.seed(113)
  res <- welch_t(rnorm(20, 100, 15), rnorm(20, 110, 25))
  expect_s3_class(res, "dmar_tbl")
  expect_type(res$value, "double")
  # The display rounds to 3 significant figures, but the stored value does not.
  disp <- format(res)
  df_stored <- res$value[res$term == "df"]
  expect_identical(disp$value[disp$term == "df"], "28.3")
  expect_false(df_stored == as.numeric(disp$value[disp$term == "df"]))
  # Whole-number rows (sample sizes) print without a decimal part.
  expect_identical(disp$value[disp$term == "n_x"], "20")
})

test_that("classed families keep their leading subclass ahead of dmar_tbl", {
  # ci_smd carries dmar_ci_smd for its tidy()/glance() methods; dmar_tbl is
  # inserted just before data.frame so print falls through to print.dmar_tbl
  # while tidy()/glance() still dispatch on the leading class.
  cs <- ci_smd(smd = 0.5, n_1 = 30, n_2 = 30)
  expect_identical(class(cs), c("dmar_ci_smd", "dmar_tbl", "data.frame"))
  expect_s3_class(generics::tidy(cs), "data.frame")

  sp <- ss_power_smd(smd = 0.5, desired_power = 0.8)
  expect_identical(class(sp), c("dmar_ss_power", "dmar_tbl", "data.frame"))
  expect_s3_class(generics::glance(sp), "data.frame")
})

test_that("p_terms rows of a rolled-out long table format to 4 decimals", {
  # equivalence_smd reports three p-values in a shared value column and names them
  # via the p_terms attribute, so they print to fixed decimals, not sig figs.
  ts <- equivalence_smd(smd = 0.1, n_1 = 50, n_2 = 50,
                 delta_lower = 0.5, delta_upper = 0.5)
  expect_s3_class(ts, "dmar_tbl")
  expect_identical(attr(ts, "p_terms"), c("p_lower", "p_upper", "p_tost"))
  expect_type(ts$value, "double")
  disp <- format(ts)
  expect_identical(disp$value[ts$term == "p_lower"], "0.0017")
  expect_identical(disp$value[ts$term == "p_upper"], "0.0242")
  expect_identical(disp$value[ts$term == "p_tost"],  "0.0242")
})

test_that("fixed_terms rows of cfa_1 format to fixed decimals, not significant figures", {
  # cfa_1(output = "verbose") reports AIC, BIC, and the H0 / H1 log-likelihoods
  # in its `estimate` column and names them via the fixed_terms attribute, so
  # they print to a fixed number of decimal places. Under three significant
  # figures a model-comparison difference of a few points would round away.
  skip_if_not_installed("lavaan")

  cov_mat <- matrix(
    c(1.384, 1.484, 1.988, 2.429, 3.031,
      1.484, 2.756, 2.874, 3.588, 4.390,
      1.988, 2.874, 4.845, 4.894, 6.080,
      2.429, 3.588, 4.894, 6.951, 7.476,
      3.031, 4.390, 6.080, 7.476, 10.313),
    nrow = 5)
  fit <- cfa_1(N = 300, S = cov_mat)

  expect_s3_class(fit, "dmar_tbl")
  expect_identical(attr(fit, "fixed_terms"), c("AIC", "BIC", "H0", "H1"))

  disp <- format(fit)
  ic <- disp$estimate[disp$term %in% c("AIC", "BIC", "H0", "H1")]
  # Every information-criterion row reads as a number with exactly three
  # decimal places, never collapsed into three significant figures (which
  # would turn 4793.531 into "4790").
  expect_true(all(grepl("^-?[0-9]+\\.[0-9]{3}$", ic)))
  expect_false(any(grepl("e", ic, fixed = TRUE)))
  # The primary (estimate) column is still numeric and the stored AIC is not
  # the three-significant-figure display value.
  expect_type(fit$estimate, "double")
  aic_stored <- fit$estimate[fit$term == "AIC"]
  expect_false(identical(aic_stored, signif(aic_stored, 3)))
})

test_that("wide-format dmar_tbl formats each numeric column on its own terms (ci_mahalanobis)", {
  # ci_mahalanobis returns a wide table: a leading sample_type label beside
  # several typed numeric columns. Degrees of freedom, sample sizes, and the
  # dimensionality count print as whole numbers; D2 and the interval limits
  # print to significant figures. The `p` column is the dimensionality (a count
  # of variables), NOT a p-value, so it keeps the whole-number rule.
  m <- ci_mahalanobis(D2 = 103.2, n_1 = 50, n_2 = 50, p = 4)
  expect_s3_class(m, "dmar_tbl")
  expect_identical(attr(m, "conf_level"), 0.95)

  fm <- format(m)
  # Whole-number columns print without a decimal part.
  expect_identical(fm$df_1, "4")
  expect_identical(fm$df_2, "95")
  expect_identical(fm$n_1,  "50")
  expect_identical(fm$n_2,  "50")
  # The `p` dimensionality count prints as a bare integer, never as a p-value
  # (which would read "4.0000").
  expect_identical(fm$p, "4")
  # D2 stays numeric (full precision) and rounds only for display.
  expect_type(m$D2, "double")
  expect_identical(fm$D2, DMAR:::.format_dmar_value(m$D2))
  # The leading character label column is left untouched.
  expect_identical(fm$sample_type, "two-sample")
})

test_that("wide-format dmar_tbl tags ci_eta_squared_generalized; footer only when a CI exists", {
  fit <- aov(len ~ supp * dose, data = ToothGrowth)

  # method = "none": the class is applied but no conf_level is attached, since
  # no interval is computed and a confidence-level footer would be misleading.
  e0 <- suppressMessages(ci_eta_squared_generalized(fit, observed = "supp"))
  expect_s3_class(e0, "dmar_tbl")
  expect_null(attr(e0, "conf_level"))
  # The leading effect label stays character; the effect-size column stays
  # numeric and is formatted column-wise like any other dmar_tbl numeric.
  expect_type(e0$effect, "character")
  expect_type(e0$eta_squared_generalized, "double")
  fe <- format(e0)
  expect_identical(fe$eta_squared_generalized,
                   DMAR:::.format_dmar_value(e0$eta_squared_generalized))

  # method = "parametric": a CI is computed, so the conf_level footer applies.
  ep <- suppressWarnings(
    ci_eta_squared_generalized(fit, observed = "supp", method = "parametric"))
  expect_s3_class(ep, "dmar_tbl")
  expect_identical(attr(ep, "conf_level"), 0.95)
})

test_that("print.dmar_tbl honors caller-supplied row.names / right without colliding", {
  set.seed(113)
  res <- welch_t(rnorm(20, 100, 15), rnorm(20, 110, 25))
  # print() forwards row.names / right to print.data.frame (the documented
  # `...` contract). These used to collide with hard-coded values and error.
  expect_error(capture.output(print(res, row.names = TRUE)), NA)
  expect_error(capture.output(print(res, right = TRUE)), NA)
  # The default still suppresses row names for the tidy look.
  out <- capture.output(print(res))
  expect_false(any(grepl("^[0-9]+ ", out)))
})
