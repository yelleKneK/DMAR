test_that(".format_dmar_value keeps whole numbers integer and others at 3 decimal places", {
  fmt <- DMAR:::.format_dmar_value
  # Whole numbers (degrees of freedom, sample sizes) print without decimals.
  expect_identical(fmt(c(1, 27, 0, -3, 30, 250)),
                   c("1", "27", "0", "-3", "30", "250"))
  # A count computed in floating point is whole up to its rounding error
  # and prints as the whole number; a value 1e-9 away is not whole.
  expect_identical(fmt(c(0.1 * 3 * 10 - 1, 22 + 4e-15, 1 - 1e-16, 1 + 1e-9)),
                   c("2", "22", "1", "1.000"))
  # Other values print to 3 decimal places by default, trailing zeros kept,
  # the examples of the 1.1.0 NEWS entry; through 1.0.x they printed to 3
  # significant figures (13, 0.56, 0.56, 233, 0.0412, 1.5).
  expect_identical(fmt(c(13.018, 0.5596, 0.56, 233.4107, 0.041234, 1.5)),
                   c("13.018", "0.560", "0.560", "233.411", "0.041", "1.500"))
  expect_identical(fmt(c(29.28272, 0.5202791, 0.2246837)),
                   c("29.283", "0.520", "0.225"))
  # Never scientific notation: a value too small for the decimals reads
  # 0.000 (1.01e-05 through 1.0.x), and a large one is written out.
  expect_identical(fmt(c(0.000123, 1.010542e-05)), c("0.000", "0.000"))
  expect_false(any(grepl("e", fmt(c(1e-20, 27, 1e15, 123456789.123)),
                         fixed = TRUE)))
  # A value that rounds to zero prints without a minus sign; one that
  # rounds away from zero keeps it. Only an exact zero is a whole zero: a
  # tiny nonzero value, a difference left by rounding error (-6.66e-16) or
  # a relative weight (3.35e-52, which printed as itself through 1.0.x),
  # is a measured value that reads 0.000.
  expect_identical(fmt(c(-0.0001, -0.0004, -6.66e-16, 3.35e-52, -0, 0)),
                   c("0.000", "0.000", "0.000", "0.000", "0", "0"))
  expect_identical(fmt(c(-6.66e-16, -0), whole = FALSE), c("0.000", "0.000"))
  expect_identical(fmt(-0.0006), "-0.001")
  # Missing values (NaN among them) are preserved as NA, not "NA" coerced
  # through format; infinities print as such.
  expect_true(all(is.na(fmt(c(NA_real_, NaN)))))
  expect_identical(fmt(c(Inf, -Inf)), c("Inf", "-Inf"))
})

test_that(".format_dmar_value prints a measured whole number to the decimals", {
  fmt <- DMAR:::.format_dmar_value
  # whole = FALSE marks a measured value, which keeps its decimals when it
  # lands on a whole number (a limit of 0, a CFI of 1).
  expect_identical(fmt(c(0, 1, -2, 0.5, 1 - 1e-16), whole = FALSE),
                   c("0.000", "1.000", "-2.000", "0.500", "1.000"))
  # `whole` is recycled along the values, so a caller types each element.
  expect_identical(fmt(c(98, 0, 0.083), whole = c(TRUE, FALSE, FALSE)),
                   c("98", "0.000", "0.083"))
})

test_that(".format_dmar_value reads digits as decimal places and checks it", {
  fmt <- DMAR:::.format_dmar_value
  expect_identical(fmt(29.28272, digits = 5L), "29.28272")
  expect_identical(fmt(29.28272, digits = 2L), "29.28")
  expect_identical(fmt(c(29.28272, -0.4), digits = 0L), c("29", "0"))
  # A number of decimal places is a whole number of 0 or more; signif()
  # took a fractional or negative value through 1.0.x.
  for (bad in list(-1, 2.5, NA_real_, "3", c(2, 3), Inf)) {
    expect_error(fmt(1.5, digits = bad), "whole number of 0 or more")
  }
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
  # rounds to the decimals rather than reading as below a threshold, and a
  # tiny negative one prints without a minus sign.
  expect_identical(fmt(c(1e-05, -1e-05)), c("0.000", "0.000"))
  expect_identical(fmt(c(Inf, -Inf)), c("Inf", "-Inf"))
  # The digits_fixed argument moves the number of decimals.
  expect_identical(fmt(2284.830284, digits_fixed = 2L), "2284.83")
  # Missing values are preserved.
  expect_true(is.na(fmt(NA_real_)))
  # The text helper behind it drops the sign of a rounded zero only, so an
  # infinite limit keeps its sign, under either decimal mark.
  expect_identical(DMAR:::.dmar_fixed_text(c(-Inf, -0.0004, -0.0006, Inf), 3L),
                   c("-Inf", "0.000", "-0.001", "Inf"))
  withr::local_options(OutDec = ",")
  expect_identical(DMAR:::.dmar_fixed_text(c(-0.0004, -1.5), 2L),
                   c("0,00", "-1,50"))
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
  # F to 3 decimal places, not full precision.
  expect_true(any(grepl("F_value\\s+57\\.566\\s*$", out)))
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

  # digits counts decimal places, and whole numbers keep none.
  out5 <- capture.output(print(res, digits = 5))
  expect_true(any(grepl("F_value\\s+57\\.56607\\s*$", out5)))
  expect_true(any(grepl("^\\s*df_2\\s+27\\s*$", out5)))

  old <- getOption("dmar.digits")
  on.exit(options(dmar.digits = old), add = TRUE)
  options(dmar.digits = 2)
  out2 <- capture.output(print(res))
  expect_true(any(grepl("F_value\\s+57\\.57\\s*$", out2)))
  # An option that is not a number of decimal places stops the print.
  options(dmar.digits = 2.5)
  expect_error(print(res), "whole number of 0 or more")
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
  # The non-p_value numeric column prints to 3 decimal places, and a whole
  # number there still prints without a decimal part.
  expect_identical(fx$estimate, c("1.235", "2", "0.250"))
})

test_that("the limits of an interval and its estimate keep their decimals when whole", {
  # A long table: the ci_lower and ci_upper rows with the estimate
  # they belong to (its first row), a limit pair named for its estimate
  # (tau2_ci_lower, tau2_ci_upper) with that estimate, the heterogeneity index tau
  # (a measured kind), and beside them counts and a whole planning value
  # (width), which keep the whole-number rule. A p-value row named like a
  # limit (p_lower) prints as a p-value.
  long <- DMAR:::.as_dmar_tbl(data.frame(
    term = c("R2", "ci_lower", "ci_upper", "df_1", "N", "tau2",
             "tau2_ci_lower", "tau2_ci_upper", "tau", "width", "p_lower"),
    value = c(0.05, 0, 0.112, 5, 50, 0, 0, 0.048, 0, 1, 1),
    stringsAsFactors = FALSE), p_terms = "p_lower")
  expect_identical(format(long)$value,
                   c("0.050", "0.000", "0.112", "5", "50", "0.000", "0.000",
                     "0.048", "0.000", "1", "1.0000"))

  # A wide table: the limit columns, the column tidy() takes as their
  # estimate (eta_squared), and the column a limit is named for
  # (eta_squared_partial). A whole F beside one that is not reads 4.000,
  # and the counts print as whole numbers.
  wide <- DMAR:::.as_dmar_tbl(data.frame(
    effect = c("a", "b"), eta_squared = c(0, 0.0141),
    ci_lower = c(0, 0), ci_upper = c(0.0498, 0.0608),
    F_value = c(4, 4.27), df_1 = c(1, 5), N = c(310, 310),
    eta_squared_partial = c(0.2, 0), eta_squared_partial_ci_lower = c(0, 0.01),
    stringsAsFactors = FALSE))
  fw <- format(wide)
  expect_identical(fw$eta_squared, c("0.000", "0.014"))
  expect_identical(fw$ci_lower, c("0.000", "0.000"))
  expect_identical(fw$ci_upper, c("0.050", "0.061"))
  expect_identical(fw$F_value, c("4.000", "4.270"))
  expect_identical(fw$df_1, c("1", "5"))
  expect_identical(fw$N, c("310", "310"))
  expect_identical(fw$eta_squared_partial, c("0.200", "0.000"))
  expect_identical(fw$eta_squared_partial_ci_lower, c("0.000", "0.010"))

  # A coefficient table whose estimate column also holds fit rows: a row
  # with limits beside it (a loading fixed at 1, with its interval [1, 1])
  # is an estimate, and the fit rows are counts. The standard error column
  # holds values that are not whole, so a fixed parameter's 0 reads 0.000.
  coef <- DMAR:::.as_dmar_tbl(data.frame(
    term = c("f =~ y1", "f =~ y2", "df", "N"),
    estimate = c(1, 0.8, 24, 301), se = c(0, 0.05, NA, NA),
    ci_lower = c(1, 0.7, NA, NA), ci_upper = c(1, 0.9, NA, NA),
    stringsAsFactors = FALSE))
  fc <- format(coef)
  expect_identical(fc$estimate, c("1.000", "0.800", "24", "301"))
  expect_identical(fc$se, c("0.000", "0.050", NA, NA))
  expect_identical(fc$ci_lower, c("1.000", "0.700", NA, NA))
  expect_identical(fc$ci_upper, c("1.000", "0.900", NA, NA))

  # A real interval truncated at zero: the lower limit of R2 reads 0.000,
  # where it read 0 through 1.0.x.
  r2 <- format(ci_R2(R2 = 0.05, N = 50, p = 5))
  expect_identical(r2$value[r2$term == "ci_lower"], "0.000")
  expect_identical(r2$value[r2$term == "R2"], "0.050")
})

test_that("a measured kind, a column of measured values, and named terms keep their decimals", {
  # A fit index is measured wherever it appears: a CFI of 1 reads 1.000
  # beside a TLI of 1.012 in the long table of sem_fit_measures(), while
  # the degrees of freedom, N, and groups stay counts.
  fm <- format(sem_fit_measures(chisq = 2.109, df = 4, N = 301, p = 6,
                                baseline_chisq = 600, baseline_df = 15))
  shown <- stats::setNames(fm$value, fm$term)
  expect_identical(unname(shown[c("cfi", "tli", "df", "N", "n_groups")]),
                   c("1.000", "1.012", "4", "301", "1"))
  expect_identical(unname(shown[c("rmsea", "rmsea_ci_lower")]),
                   c("0.000", "0.000"))

  # Powers, the proportions of a sensitivity analysis, a fraction of
  # missing information, a model comparison weight, a standard deviation, a
  # moment, a test statistic, a standard error, and a noncentrality
  # parameter, as rows of a long table or as columns of a wide one whose
  # values are all whole. A name that is also a count (delta_df, median_N)
  # stays a count.
  long <- DMAR:::.as_dmar_tbl(data.frame(
    term = c("power_a", "actual_power", "pct_ci_less_width",
             "pct_ci_less_width_a", "total_type_I_error", "fmi",
             "relative_weight", "sd_ci_width", "skewness", "F_value", "se",
             "ncp", "delta_df", "median_N", "necessary_N"),
    value = c(1, 1, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 3, 430, 120),
    stringsAsFactors = FALSE))
  expect_identical(format(long)$value,
                   c("1.000", "1.000", "1.000", "1.000", "0.000", "0.000",
                     "1.000", "0.000", "0.000", "2.000", "0.000", "0.000",
                     "3", "430", "120"))
  # A one-row wide table has no column of values that are not whole, and
  # its test statistic is measured by name.
  one <- DMAR:::.as_dmar_tbl(data.frame(
    effect = "overall", omega_squared = 0.1, ci_lower = 0,
    ci_upper = 0.456, F_value = 2, df_effect = 2, df_error = 15, N = 18,
    stringsAsFactors = FALSE))
  expect_identical(unlist(format(one)[-1L], use.names = FALSE),
                   c("0.100", "0.000", "0.456", "2.000", "2", "15", "18"))
  wide <- DMAR:::.as_dmar_tbl(data.frame(
    level = c("configural", "metric"), cfi = c(1, 1), rmsea = c(0, 0),
    delta_df = c(NA, 2), delta_cfi = c(NA, 0), stringsAsFactors = FALSE))
  fw <- format(wide)
  expect_identical(fw$cfi, c("1.000", "1.000"))
  expect_identical(fw$rmsea, c("0.000", "0.000"))
  expect_identical(fw$delta_df, c(NA, "2"))
  expect_identical(fw$delta_cfi, c(NA, "0.000"))

  # A column that holds values that are not whole is measured throughout,
  # unless it is named for a count (degrees of freedom adjusted for
  # sphericity) or the row records the sum-of-squares type.
  typed <- DMAR:::.as_dmar_tbl(data.frame(
    effect = c("A", "A", "B", "sum_of_squares_type"),
    statistic_value = c(0.25, 7.5, 2, 3), df_1 = c(2, 1.83, 1, NA),
    relative_weight = c(1, 1.23e-06, 3.35e-52, NA), stringsAsFactors = FALSE))
  ft <- format(typed)
  expect_identical(ft$statistic_value, c("0.250", "7.500", "2.000", "3"))
  expect_identical(ft$df_1, c("2", "1.830", "1", NA))
  expect_identical(ft$relative_weight, c("1.000", "0.000", "0.000", NA))

  # A function names columns of its own kind through measured_terms.
  own <- DMAR:::.as_dmar_tbl(data.frame(
    factor_x = c("A", "B"), congruence = c(1, 1 - 1e-16),
    stringsAsFactors = FALSE), measured_terms = "congruence")
  expect_identical(format(own)$congruence, c("1.000", "1.000"))
  attr(own, "measured_terms") <- NULL
  expect_identical(format(own)$congruence, c("1", "1"))
  # dmacs() names its loadings and intercepts: an intercept of 2 reads
  # 2.000 beside one of 2.300, and a dMACS of exactly 0 reads 0.000.
  dm <- format(dmacs(factor_loading_reference = c(0.8, 0.75),
                     factor_loading_focal = c(0.8, 0.75), intercept_reference = c(2, 2.3),
                     intercept_focal = c(2, 2), sd_pooled = c(1.2, 1.1),
                     item_names = c("optimism", "worry")))
  expect_identical(dm$intercept_focal, c("2.000", "2.000"))
  expect_identical(dm$dmacs, c("0.000", "0.273"))
})

test_that("the limits of a p-value's interval print to its decimals without the floor", {
  x <- DMAR:::.as_dmar_tbl(data.frame(
    term = c("p_bootstrap", "p_bootstrap_ci_lower", "p_bootstrap_ci_upper",
             "B"),
    value = c(0.02, 0.0041232, 0.105, 200), stringsAsFactors = FALSE),
    p_terms = "p_bootstrap")
  expect_identical(format(x)$value, c("0.0200", "0.0041", "0.1050", "200"))
  x$value[1:3] <- c(0, 0, 0.0189)
  expect_identical(format(x)$value, c("< 0.0001", "0.0000", "0.0189", "200"))
  expect_identical(format(x, digits_p = 3)$value[2:3], c("0.000", "0.019"))
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
  # The display rounds to 3 decimal places, but the stored value does not.
  disp <- format(res)
  df_stored <- res$value[res$term == "df"]
  expect_identical(disp$value[disp$term == "df"], sprintf("%.3f", df_stored))
  expect_match(disp$value[disp$term == "df"], "^28\\.[0-9]{3}$")
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

  sp <- ss_power_smd(true_smd = 0.5, desired_power = 0.8)
  expect_identical(class(sp), c("dmar_ss_power", "dmar_tbl", "data.frame"))
  expect_s3_class(generics::glance(sp), "data.frame")
})

test_that("p_terms rows of a rolled-out long table format to 4 decimals", {
  # equivalence_smd reports three p-values in a shared value column and names them
  # via the p_terms attribute, so they print to digits_p decimal places.
  ts <- equivalence_smd(smd = 0.1, n_1 = 50, n_2 = 50,
                 smd_lower = 0.5, smd_upper = 0.5)
  expect_s3_class(ts, "dmar_tbl")
  expect_identical(attr(ts, "p_terms"), c("p_lower", "p_upper", "p_tost"))
  expect_type(ts$value, "double")
  disp <- format(ts)
  expect_identical(disp$value[ts$term == "p_lower"], "0.0017")
  expect_identical(disp$value[ts$term == "p_upper"], "0.0242")
  expect_identical(disp$value[ts$term == "p_tost"],  "0.0242")
})

test_that("fixed_terms rows of cfa_1 keep digits_fixed decimals whatever digits is", {
  skip_on_cran()
  # cfa_1(output = "verbose") reports AIC, BIC, and the H0 / H1 log-likelihoods
  # in its `estimate` column and names them via the fixed_terms attribute, so
  # they print to digits_fixed decimal places. Under three significant
  # figures (the rule through 1.0.x) a model comparison difference of a few
  # points would round away, and it must not when digits is lowered.
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
  # would turn 4793.531 into "4790"), and keeps them when the rest of the
  # table is shown to one decimal place.
  expect_true(all(grepl("^-?[0-9]+\\.[0-9]{3}$", ic)))
  expect_false(any(grepl("e", ic, fixed = TRUE)))
  disp1 <- format(fit, digits = 1)
  expect_identical(disp1$estimate[disp1$term %in% c("AIC", "BIC", "H0", "H1")],
                   ic)
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
  # print to 3 decimal places. The `p` column is the dimensionality (a count
  # of variables), not a p-value, so it keeps the whole-number rule.
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
  skip_on_cran()  # The extended validation grid runs in the full local suite.
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

# ---- rbind() of DMAR tables ----------------------------------------------

test_that("rbind() of DMAR tables with the same attributes keeps the class and the attributes", {
  x <- ci_smd(smd = 0.5, n_1 = 30, n_2 = 30)
  s <- rbind(x, x)
  expect_identical(class(s), class(x))
  expect_identical(attr(s, "conf_level"), 0.95)
  expect_identical(nrow(s), 2L * nrow(x))
  expect_identical(s$value, c(x$value, x$value))
  # Every attribute but the row names is that of the tables stacked
  # (rbind.data.frame() returns them in another order).
  facts <- function(d) {
    a <- attributes(d)
    a <- a[setdiff(names(a), "row.names")]
    a[order(names(a))]
  }
  expect_identical(facts(s), facts(x))
  # A stack can be stacked again.
  expect_identical(class(rbind(s, x)), class(x))
  # One table, and NULL arguments, which rbind() sets aside.
  expect_identical(class(rbind(x)), class(x))
  expect_identical(class(rbind(NULL, x, NULL, x)), class(x))
  expect_identical(attr(rbind(NULL, x), "conf_level"), 0.95)
  # The named arguments of rbind.data.frame() are not tables.
  r <- rbind(x, x, make.row.names = FALSE)
  expect_s3_class(r, "dmar_tbl")
  expect_identical(nrow(r), 2L * nrow(x))
})

test_that("rbind() of DMAR tables at different confidence levels gives a plain data frame with no level", {
  x90 <- ci_smd(smd = 0.5, n_1 = 30, n_2 = 30, conf_level = 0.90)
  x95 <- ci_smd(smd = 0.5, n_1 = 30, n_2 = 30, conf_level = 0.95)
  m <- rbind(x90, x95)
  expect_identical(class(m), "data.frame")
  expect_null(attr(m, "conf_level"))
  expect_identical(sort(names(attributes(m))), c("class", "names", "row.names"))
  # Nothing is lost from the rows: both sets of limits are there.
  expect_identical(nrow(m), nrow(x90) + nrow(x95))
  expect_identical(m$term, c(x90$term, x95$term))
  expect_identical(m$value, c(x90$value, x95$value))
  # No footer can state one level for all rows, and no condition is given.
  expect_no_condition(rbind(x90, x95))
  expect_false(any(grepl("Confidence level", capture.output(print(m)))))
  # The level travels with tidy() rows instead.
  td <- rbind(generics::tidy(x90), generics::tidy(x95))
  expect_identical(td$conf_level, c(0.90, 0.95))
})

test_that("rbind() of a DMAR table with a plain data frame gives a plain data frame, in either order", {
  x <- ci_smd(smd = 0.5, n_1 = 30, n_2 = 30)
  plain <- data.frame(term = "note", value = 1, stringsAsFactors = FALSE)
  a <- rbind(x, plain)
  b <- rbind(plain, x)
  expect_identical(class(a), "data.frame")
  expect_identical(class(b), "data.frame")
  expect_null(attr(a, "conf_level"))
  expect_null(attr(b, "conf_level"))
  expect_identical(a$term, c(x$term, "note"))
  expect_identical(b$term, c("note", x$term))
})

test_that("rbind() of tables from different functions gives a plain data frame", {
  a <- ci_smd(smd = 0.5, n_1 = 30, n_2 = 30)
  b <- moments_nc_F(3, 40, 8)
  s <- rbind(a, b)
  expect_identical(class(s), "data.frame")
  expect_identical(nrow(s), nrow(a) + nrow(b))
})

test_that("do.call(rbind, ) on a list of three follows the same rule", {
  x90 <- ci_smd(smd = 0.5, n_1 = 30, n_2 = 30, conf_level = 0.90)
  x95 <- ci_smd(smd = 0.5, n_1 = 30, n_2 = 30, conf_level = 0.95)
  same <- do.call(rbind, list(x95, x95, x95))
  expect_identical(class(same), class(x95))
  expect_identical(attr(same, "conf_level"), 0.95)
  expect_identical(nrow(same), 3L * nrow(x95))
  mixed <- do.call(rbind, list(x90, x95, x90))
  expect_identical(class(mixed), "data.frame")
  expect_null(attr(mixed, "conf_level"))
  expect_identical(nrow(mixed), 3L * nrow(x95))
})
