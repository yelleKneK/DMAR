# Reads one row of the term / value table.
nnt_row <- function(res, term) res$value[res$term == term]
# The standardized mean difference whose success rate difference is `s`.
smd_of_srd <- function(s) sqrt(2) * qnorm((s + 1) / 2)

test_that("nnt_from_smd() point estimate matches Kraemer-Kupfer formula", {
  res <- nnt_from_smd(smd = 0.5)
  srd_expected <- 2 * pnorm(0.5 / sqrt(2)) - 1
  expect_equal(res$value[res$term == "srd"], srd_expected, tolerance = 1e-12)
  expect_equal(res$value[res$term == "nnt"], 1 / srd_expected, tolerance = 1e-12)
})

test_that("nnt_from_smd() returns Inf NNT when d = 0", {
  res <- nnt_from_smd(smd = 0)
  expect_equal(res$value[res$term == "srd"], 0)
  expect_true(is.infinite(res$value[res$term == "nnt"]))
})

test_that("nnt_from_smd() reproduces Kraemer and Kupfer (2006) Table 1 NNTs", {
  # Kraemer & Kupfer (2006), Table 1 (p. 992) and text (p. 994): d = 0.2,
  # 0.5, 0.8 give NNT = 8.892, 3.619, 2.334, quoted in the text as 8.9,
  # 3.6, 2.3.
  nnt <- vapply(c(0.2, 0.5, 0.8),
                function(d) nnt_from_smd(d)$value[3], numeric(1L))
  expect_equal(round(nnt, 1), c(8.9, 3.6, 2.3))
  expect_equal(round(nnt, 3), c(8.892, 3.619, 2.334))
})

test_that("nnt_from_smd() with CI on d propagates monotonically (smaller d -> larger NNT bound)", {
  res <- nnt_from_smd(smd = 0.5, smd_ci_lower = 0.20, smd_ci_upper = 0.80)
  nnt_lo <- res$value[res$term == "nnt_ci_lower"]
  nnt_hi <- res$value[res$term == "nnt_ci_upper"]
  expect_true(nnt_hi > nnt_lo)
  # Bound at d_upper = 0.8 should match d=0.8 point estimate:
  expect_equal(nnt_lo, nnt_from_smd(smd = 0.80)$value[3], tolerance = 1e-12)
})

test_that("nnt_from_smd() with sample sizes uses ci_smd() internally", {
  res <- nnt_from_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.95)
  expect_true("nnt_ci_lower" %in% res$term)
  expect_true("nnt_ci_upper" %in% res$term)
})

test_that("nnt_from_smd() rows keep their 1.0.0 positions, with the limits named for their estimate and the SRD limits appended", {
  point <- nnt_from_smd(smd = 0.5)
  expect_identical(point$term, c("smd", "srd", "nnt"))
  expect_true(is.numeric(point$value))
  expect_null(attr(point, "nnt_set"))
  expect_null(attr(point, "nnt_confidence_set"))

  nine <- c("smd", "srd", "nnt", "smd_ci_lower", "smd_ci_upper",
            "nnt_ci_lower", "nnt_ci_upper", "srd_ci_lower", "srd_ci_upper")
  res <- nnt_from_smd(smd = 0.5, n_1 = 50, n_2 = 50)
  expect_identical(res$term, nine)
  expect_true(is.numeric(res$value))
  # The values DMAR 1.0.0 returned in rows 1 to 7 for this call.
  expect_equal(res$value[1:7],
               c(0.5, 0.2763263902, 3.6189087817, 0.1005857103,
                 0.8969414275, 2.1093866118, 17.6361867918),
               tolerance = 1e-9)
  expect_equal(res$value[6:7], c(2.109, 17.636), tolerance = 1e-3)
  expect_identical(attr(res, "nnt_set"), "interval")

  expect_identical(nnt_from_smd(0.3, smd_ci_lower = -0.1, smd_ci_upper = 0.7)$term,
                   nine)
})

test_that("nnt_from_smd() reports the SRD at the two limits as an ordered interval", {
  for (lim in list(c(0.2, 0.8), c(-0.1, 0.7), c(-0.8, -0.2))) {
    res <- nnt_from_smd(smd = mean(lim), smd_ci_lower = lim[1], smd_ci_upper = lim[2])
    expect_equal(nnt_row(res, "srd_ci_lower"), 2 * pnorm(lim[1] / sqrt(2)) - 1,
                 tolerance = 1e-12)
    expect_equal(nnt_row(res, "srd_ci_upper"), 2 * pnorm(lim[2] / sqrt(2)) - 1,
                 tolerance = 1e-12)
    expect_lt(nnt_row(res, "srd_ci_lower"), nnt_row(res, "srd_ci_upper"))
  }
  res <- nnt_from_smd(smd = 0.5, n_1 = 50, n_2 = 50)
  expect_equal(nnt_row(res, "srd_ci_lower"),
               2 * pnorm(nnt_row(res, "smd_ci_lower") / sqrt(2)) - 1,
               tolerance = 1e-12)
  expect_equal(nnt_row(res, "srd_ci_upper"),
               2 * pnorm(nnt_row(res, "smd_ci_upper") / sqrt(2)) - 1,
               tolerance = 1e-12)
})

test_that("nnt_from_smd() returns two intervals when the interval for d includes zero", {
  res <- nnt_from_smd(smd = 0.3, smd_ci_lower = -0.1, smd_ci_upper = 0.7)
  expect_true(is.na(nnt_row(res, "nnt_ci_lower")))
  expect_true(is.na(nnt_row(res, "nnt_ci_upper")))
  expect_identical(attr(res, "nnt_set"), "two_intervals")
  set <- attr(res, "nnt_confidence_set")
  expect_equal(set,
               data.frame(lower = c(-Inf, 2.635865),
                          upper = c(-17.73931, Inf)),
               tolerance = 1e-6)
  # The two finite ends are the reciprocals of the SRD limits.
  expect_equal(set$upper[1], 1 / nnt_row(res, "srd_ci_lower"), tolerance = 1e-12)
  expect_equal(set$lower[2], 1 / nnt_row(res, "srd_ci_upper"), tolerance = 1e-12)
  # The estimate lies in the benefit interval.
  nnt <- nnt_row(res, "nnt")
  expect_equal(nnt, 5.9525, tolerance = 1e-4)
  expect_true(nnt >= set$lower[2] && nnt <= set$upper[2])

  # The noncentral t route reaches the same case.
  res_n <- nnt_from_smd(smd = 0.1, n_1 = 30, n_2 = 30)
  expect_identical(attr(res_n, "nnt_set"), "two_intervals")
  expect_true(all(is.na(res_n$value[6:7])))
  set_n <- attr(res_n, "nnt_confidence_set")
  expect_equal(c(set_n$upper[1], set_n$lower[2]), c(-4.417, 3.015),
               tolerance = 1e-3)
})

test_that("nnt_from_smd() places the infinite end when a limit for d is exactly zero", {
  res <- nnt_from_smd(smd = 0.3, smd_ci_lower = 0, smd_ci_upper = 0.7)
  expect_equal(nnt_row(res, "nnt_ci_lower"), 2.635865, tolerance = 1e-6)
  expect_identical(nnt_row(res, "nnt_ci_upper"), Inf)
  expect_identical(attr(res, "nnt_set"), "interval")

  res <- nnt_from_smd(smd = -0.3, smd_ci_lower = -0.7, smd_ci_upper = 0)
  expect_identical(nnt_row(res, "nnt_ci_lower"), -Inf)
  expect_equal(nnt_row(res, "nnt_ci_upper"), -2.635865, tolerance = 1e-6)
  expect_identical(attr(res, "nnt_set"), "interval")
  expect_equal(attr(res, "nnt_confidence_set"),
               data.frame(lower = -Inf, upper = -2.635865), tolerance = 1e-6)

  res <- nnt_from_smd(smd = 0, smd_ci_lower = 0, smd_ci_upper = 0)
  expect_identical(nnt_row(res, "nnt_ci_lower"), Inf)
  expect_identical(nnt_row(res, "nnt_ci_upper"), Inf)
  expect_identical(attr(res, "nnt_set"), "interval")
})

test_that("nnt_from_smd() gives an ordered negative interval when both limits for d are negative", {
  res <- nnt_from_smd(smd = -0.5, smd_ci_lower = -0.8, smd_ci_upper = -0.2)
  lo <- nnt_row(res, "nnt_ci_lower")
  hi <- nnt_row(res, "nnt_ci_upper")
  expect_equal(c(lo, hi), c(-8.892, -2.334), tolerance = 1e-3)
  expect_lt(lo, hi)
  expect_true(nnt_row(res, "nnt") > lo && nnt_row(res, "nnt") < hi)
  expect_identical(attr(res, "nnt_set"), "interval")
})

test_that("nnt_from_smd() decides the case on the limits for d, not on an SRD that underflows", {
  res <- nnt_from_smd(smd = 0.3, smd_ci_lower = -1e-320, smd_ci_upper = 0.7)
  expect_identical(attr(res, "nnt_set"), "two_intervals")
  expect_true(is.na(nnt_row(res, "nnt_ci_lower")))
  expect_true(is.na(nnt_row(res, "nnt_ci_upper")))
})

test_that("nnt_from_smd() reproduces Altman's (1998) reciprocal limits", {
  # Altman (1998): a difference of 10% with limits 5% and 15% gives an NNT
  # interval of 6.7 to 20.
  res <- nnt_from_smd(smd = smd_of_srd(0.10),
                      smd_ci_lower = smd_of_srd(0.05),
                      smd_ci_upper = smd_of_srd(0.15))
  expect_equal(c(nnt_row(res, "nnt_ci_lower"), nnt_row(res, "nnt_ci_upper")),
               c(1 / 0.15, 20), tolerance = 1e-8)
  expect_equal(round(nnt_row(res, "nnt_ci_lower"), 1), 6.7)

  # Limits -5% and 25%: NNTH 20 to infinity to NNTB 4.
  res <- nnt_from_smd(smd = smd_of_srd(0.10),
                      smd_ci_lower = smd_of_srd(-0.05),
                      smd_ci_upper = smd_of_srd(0.25))
  set <- attr(res, "nnt_confidence_set")
  expect_identical(attr(res, "nnt_set"), "two_intervals")
  expect_equal(set, data.frame(lower = c(-Inf, 4), upper = c(-20, Inf)),
               tolerance = 1e-8)

  # Limits -0.27 and 0.11: NNTH 3.7 to infinity to NNTB 9.1.
  res <- nnt_from_smd(smd = smd_of_srd(-0.08),
                      smd_ci_lower = smd_of_srd(-0.27),
                      smd_ci_upper = smd_of_srd(0.11))
  set <- attr(res, "nnt_confidence_set")
  expect_equal(c(set$upper[1], set$lower[2]), c(-1 / 0.27, 1 / 0.11),
               tolerance = 1e-8)
  expect_equal(round(c(set$upper[1], set$lower[2]), 3), c(-3.704, 9.091))
})

test_that("nnt_from_smd() confidence set holds the true NNT exactly when the interval for d holds delta", {
  skip_on_cran()
  in_set <- function(nnt, set) any(nnt >= set$lower & nnt <= set$upper)
  true_nnt <- function(delta) {
    if (delta == 0) Inf else 1 / (2 * pnorm(delta / sqrt(2)) - 1)
  }
  sizes <- list(c(10, 10), c(20, 35), c(60, 60))
  checked <- 0L
  for (n in sizes) for (d in seq(-1, 1, by = 0.25)) {
    res <- nnt_from_smd(smd = d, n_1 = n[1], n_2 = n[2])
    set <- attr(res, "nnt_confidence_set")
    lower <- nnt_row(res, "smd_ci_lower")
    upper <- nnt_row(res, "smd_ci_upper")
    for (delta in c(-0.8, -0.3, -0.05, 0, 0.05, 0.4, 0.9)) {
      expect_identical(in_set(true_nnt(delta), set),
                       delta >= lower && delta <= upper)
      checked <- checked + 1L
    }
  }
  expect_identical(checked, 189L)
})

test_that("nnt_from_smd() prints the two intervals under the table, and only then", {
  two <- nnt_from_smd(smd = 0.3, smd_ci_lower = -0.1, smd_ci_upper = 0.7)
  expect_output(print(two), "two intervals")
  expect_output(print(two),
                "NNT <= -17.739 (harm) or NNT >= 2.636 (benefit), joined at infinity (Altman, 1998)",
                fixed = TRUE)
  expect_output(print(two, digits = 1), "NNT <= -17.7 (harm) or NNT >= 2.6 (benefit)",
                fixed = TRUE)
  one <- nnt_from_smd(smd = 0.5, smd_ci_lower = 0.2, smd_ci_upper = 0.8)
  expect_false(any(grepl("two intervals", capture.output(print(one)))))
  expect_false(any(grepl("two intervals",
                         capture.output(print(nnt_from_smd(smd = 0.5))))))
})

test_that("nnt_from_smd() states a confidence level only for limits whose level it knows", {
  supplied <- nnt_from_smd(smd = 0.5, smd_ci_lower = 0.2, smd_ci_upper = 0.8)
  expect_null(attr(supplied, "conf_level"))
  expect_false(any(grepl("Confidence level", capture.output(print(supplied)))))

  stated <- nnt_from_smd(smd = 0.5, smd_ci_lower = 0.2, smd_ci_upper = 0.8,
                         conf_level = 0.90)
  expect_identical(attr(stated, "conf_level"), 0.90)
  expect_output(print(stated), "Confidence level: 90%", fixed = TRUE)

  computed <- nnt_from_smd(smd = 0.5, n_1 = 50, n_2 = 50)
  expect_identical(attr(computed, "conf_level"), 0.95)
  expect_output(print(computed), "Confidence level: 95%", fixed = TRUE)

  expect_null(attr(nnt_from_smd(smd = 0.5), "conf_level"))
})

test_that("nnt_from_smd() refuses arguments it cannot use", {
  skip_on_cran()
  smd_msg <- "'smd' must be a single finite number."
  expect_error(nnt_from_smd(NA_real_), smd_msg, fixed = TRUE)
  expect_error(nnt_from_smd(Inf), smd_msg, fixed = TRUE)
  expect_error(nnt_from_smd(c(0.2, 0.5)), smd_msg, fixed = TRUE)
  expect_error(nnt_from_smd("0.5"), smd_msg, fixed = TRUE)

  level_msg <- "'conf_level' must be a single number in (0, 1)."
  expect_error(nnt_from_smd(0.3, smd_ci_lower = 0.1, smd_ci_upper = 0.5,
                            conf_level = 7), level_msg, fixed = TRUE)
  expect_error(nnt_from_smd(0.3, n_1 = 20, n_2 = 20, conf_level = 1),
               level_msg, fixed = TRUE)
  expect_error(nnt_from_smd(0.3, conf_level = NA_real_), level_msg,
               fixed = TRUE)

  both_limits <- "Supply both 'smd_ci_lower' and 'smd_ci_upper', or neither."
  expect_error(nnt_from_smd(0.3, smd_ci_lower = 0.1), both_limits, fixed = TRUE)
  expect_error(nnt_from_smd(0.3, smd_ci_upper = 0.5), both_limits, fixed = TRUE)

  limit_msg <- "'smd_ci_lower' and 'smd_ci_upper' must each be a single finite number."
  expect_error(nnt_from_smd(0.3, smd_ci_lower = NA_real_, smd_ci_upper = 0.5),
               limit_msg, fixed = TRUE)
  expect_error(nnt_from_smd(0.3, smd_ci_lower = 0.1, smd_ci_upper = Inf),
               limit_msg, fixed = TRUE)
  expect_error(nnt_from_smd(0.3, smd_ci_lower = c(0.1, 0.2), smd_ci_upper = 0.5),
               limit_msg, fixed = TRUE)
  expect_error(nnt_from_smd(0.3, smd_ci_lower = 0.7, smd_ci_upper = 0.1),
               "'smd_ci_lower' must not exceed 'smd_ci_upper'.", fixed = TRUE)

  both_n <- "Supply both 'n_1' and 'n_2', or neither."
  expect_error(nnt_from_smd(0.3, n_1 = 10), both_n, fixed = TRUE)
  expect_error(nnt_from_smd(0.3, n_2 = 10), both_n, fixed = TRUE)

  n_msg <- paste0("'n_1' and 'n_2' must each be a single number of at least 1, ",
                  "with n_1 + n_2 greater than 2.")
  expect_error(nnt_from_smd(0.3, n_1 = NA, n_2 = 10), n_msg, fixed = TRUE)
  expect_error(nnt_from_smd(0.3, n_1 = "a", n_2 = 10), n_msg, fixed = TRUE)
  expect_error(nnt_from_smd(0.3, n_1 = 0, n_2 = 10), n_msg, fixed = TRUE)
  expect_error(nnt_from_smd(0.3, n_1 = 1, n_2 = 1), n_msg, fixed = TRUE)
  expect_error(nnt_from_smd(0.3, n_1 = c(10, 12), n_2 = 10), n_msg,
               fixed = TRUE)
})

test_that("nnt_from_smd() warns about an estimate outside its limits and still returns", {
  expect_warning(
    res <- nnt_from_smd(2, smd_ci_lower = 0.1, smd_ci_upper = 0.5),
    "'smd' lies outside the interval from 'smd_ci_lower' to 'smd_ci_upper'.",
    fixed = TRUE)
  expect_equal(nnt_row(res, "nnt_ci_lower"), nnt_from_smd(0.5)$value[3],
               tolerance = 1e-12)
  expect_no_warning(nnt_from_smd(0.1, smd_ci_lower = 0.1, smd_ci_upper = 0.5))
})

test_that("nnt_from_smd() uses supplied limits and leaves n_1, n_2 unused", {
  limits_only <- nnt_from_smd(0.5, smd_ci_lower = 0.2, smd_ci_upper = 0.8)
  expect_no_error(
    with_n <- nnt_from_smd(0.5, smd_ci_lower = 0.2, smd_ci_upper = 0.8,
                           n_1 = 10, n_2 = 10))
  expect_identical(with_n, limits_only)
})

test_that("nnt_from_smd() results work with tidy() and glance()", {
  res <- nnt_from_smd(smd = 0.3, smd_ci_lower = -0.1, smd_ci_upper = 0.7)
  td <- generics::tidy(res)
  expect_identical(nrow(td), 9L)
  expect_identical(td$term, res$term)
  expect_identical(nrow(generics::glance(res)), 1L)
})

test_that("nnt_from_smd() keeps the sign of an infinite end, the SRD near zero, and the warning boundary", {
  # A lower limit so small that its SRD rounds to zero from below: the
  # harm interval ends at -Inf, not Inf.
  set <- attr(nnt_from_smd(0, smd_ci_lower = -1e-320, smd_ci_upper = 0.7),
              "nnt_confidence_set")
  expect_identical(set$upper[1], -Inf)
  # The SRD near zero is d / sqrt(pi) to first order; 2 pnorm(d / sqrt(2))
  # - 1 loses digits there (5.642153e-13 at d = 1e-12).
  expect_equal(nnt_from_smd(1e-12)$value[2], 1e-12 / sqrt(pi),
               tolerance = 1e-12)
  # An estimate equal to a limit is not outside its limits.
  expect_no_warning(nnt_from_smd(0.5, smd_ci_lower = 0.1, smd_ci_upper = 0.5))
  expect_no_warning(nnt_from_smd(0.1, smd_ci_lower = 0.1, smd_ci_upper = 0.5))
  expect_warning(nnt_from_smd(0.6, smd_ci_lower = 0.1, smd_ci_upper = 0.5))
})

test_that("nnt_from_smd() results have the class dmar_nnt and tidy() keeps the SRD limits", {
  res <- nnt_from_smd(0.3, smd_ci_lower = -0.1, smd_ci_upper = 0.7)
  expect_identical(class(res), c("dmar_nnt", "dmar_tbl", "data.frame"))
  expect_identical(class(nnt_from_smd(0.5)),
                   c("dmar_nnt", "dmar_tbl", "data.frame"))
  td <- generics::tidy(res)
  set <- attr(res, "nnt_confidence_set")
  srd <- td$estimate[td$term %in% c("srd_ci_lower", "srd_ci_upper")]
  expect_equal(1 / srd, c(set$upper[1], set$lower[2]))
})

test_that("nnt_from_smd() reproduces the numbers its help page states", {
  # The table of the conversion (SRD to three decimals, NNT to three).
  d <- seq(0.1, 1, by = 0.1)
  srd <- vapply(d, function(x) nnt_row(nnt_from_smd(x), "srd"), numeric(1))
  nnt <- vapply(d, function(x) nnt_row(nnt_from_smd(x), "nnt"), numeric(1))
  expect_equal(round(srd, 3), c(0.056, 0.112, 0.168, 0.223, 0.276, 0.329,
                                0.379, 0.428, 0.475, 0.520))
  expect_equal(round(nnt, 3), c(17.739, 8.892, 5.953, 4.490, 3.619, 3.043,
                                2.636, 2.334, 2.103, 1.921))
  # The worked example with an interval that excludes zero.
  res <- nnt_from_smd(smd = 0.5, n_1 = 50, n_2 = 50)
  expect_equal(round(res$value, 3),
               c(0.5, 0.276, 3.619, 0.101, 0.897, 2.109, 17.636, 0.057,
                 0.474))
  expect_identical(attr(res, "nnt_set"), "interval")
  # The worked example with an interval that includes zero.
  res <- nnt_from_smd(smd = 0.2, n_1 = 50, n_2 = 50)
  expect_equal(round(res$value[c(1:5, 8:9)], 3),
               c(0.2, 0.112, 8.892, -0.193, 0.592, -0.109, 0.325))
  expect_true(all(is.na(res$value[6:7])))
  expect_identical(attr(res, "nnt_set"), "two_intervals")
  expect_equal(round(attr(res, "nnt_confidence_set")$upper[1], 3), -9.189)
  expect_equal(round(attr(res, "nnt_confidence_set")$lower[2], 3), 3.079)
  expect_output(print(res),
                "NNT <= -9.189 (harm) or NNT >= 3.079 (benefit)", fixed = TRUE)
  # The threshold NNT under the same normal model, at d = 0.5: larger than
  # the Kraemer and Kupfer NNT at every control response proportion, with
  # its smallest value at p_0 = pnorm(-d / 2).
  threshold_nnt <- function(d, p_0) 1 / (pnorm(qnorm(p_0) + d) - p_0)
  kk <- nnt_row(nnt_from_smd(0.5), "nnt")
  expect_equal(round(kk, 2), 3.62)
  expect_equal(round(1 / (2 * pnorm(0.25) - 1), 2), 5.07)
  expect_equal(round(pnorm(-0.25), 2), 0.40)
  expect_equal(threshold_nnt(0.5, pnorm(-0.25)), 1 / (2 * pnorm(0.25) - 1))
  expect_equal(round(threshold_nnt(0.5, c(0.5, 0.2)), 2), c(5.22, 6.01))
  p_0 <- seq(0.01, 0.99, by = 0.01)
  expect_true(all(threshold_nnt(0.5, p_0) >= 1 / (2 * pnorm(0.25) - 1) - 1e-12))
  expect_true(all(threshold_nnt(0.5, p_0) > kk))
  # The edge cases the page lists.
  neg <- nnt_from_smd(smd = -0.5, n_1 = 50, n_2 = 50)
  expect_equal(round(neg$value[6:7], 3), c(-17.636, -2.109))
  expect_equal(nnt_row(nnt_from_smd(1e-8), "nnt"), 1 / (1e-8 / sqrt(pi)),
               tolerance = 1e-12)
  tiny <- nnt_from_smd(smd = 0.5, n_1 = 1, n_2 = 2)
  expect_equal(round(tiny$value[4:5], 2), c(-2.07, 2.88))
  expect_identical(attr(tiny, "nnt_set"), "two_intervals")
  expect_error(nnt_from_smd(0.3, conf_level = NULL), "conf_level")
  expect_warning(nnt_from_smd(0.9, smd_ci_lower = 0.1, smd_ci_upper = 0.5),
                 "lies outside")
})

test_that("nnt_from_smd() takes and reports the limits on d as smd_ci_lower and smd_ci_upper", {
  # DMAR 1.0.x called the limits smd_lower and smd_upper; those names are
  # the equivalence bounds of equivalence_smd() and stop here with
  # "unused argument".
  res <- nnt_from_smd(smd = 0.5, smd_ci_lower = 0.20, smd_ci_upper = 0.80)
  expect_identical(res$value[res$term == "smd_ci_lower"], 0.20)
  expect_identical(res$value[res$term == "smd_ci_upper"], 0.80)
  expect_equal(res$value[res$term == "srd_ci_lower"],
               2 * pnorm(0.20 / sqrt(2)) - 1, tolerance = 1e-15)
  expect_equal(res$value[res$term == "nnt_ci_upper"],
               1 / (2 * pnorm(0.20 / sqrt(2)) - 1), tolerance = 1e-12)
  expect_error(nnt_from_smd(smd = 0.5, smd_lower = 0.20, smd_upper = 0.80),
               "unused argument")
  err <- expect_error(nnt_from_smd(smd = 0.5, smd_ci_lower = 0.20),
                      class = "dmar_argument_invalid")
  expect_identical(err$arg, c("smd_ci_lower", "smd_ci_upper"))
})

test_that("results_sentence() reports the SRD on the table, its tidy() rows, stacks and subsets of them", {
  skip_on_cran()
  # The limits on d are now named for d (smd_ci_lower, smd_ci_upper), so
  # the verbs see the interval for d beside those of the SRD and the NNT.
  # The table names d as the source of the other two (ci_source), and the
  # sentence is the SRD's wherever the rows go, as it was in DMAR 1.0.x,
  # when the limits on d were smd_lower and smd_upper. Oracle: the limits
  # on d by inverting the noncentral t with uniroot() and pt(), and the
  # SRD as 2 pnorm(d / sqrt(2)) - 1 at the estimate and at each limit.
  srd_of <- function(d) 2 * pnorm(d / sqrt(2)) - 1
  n_1 <- 50; n_2 <- 40; d <- 0.5
  scale <- sqrt(n_1 * n_2 / (n_1 + n_2))
  ncp_at <- function(p) uniroot(function(ncp) pt(d * scale, n_1 + n_2 - 2, ncp) - p,
                                c(-3, 8), tol = 1e-12)$root
  s <- srd_of(c(d, c(ncp_at(.975), ncp_at(.025)) / scale))
  want <- sprintf("srd = %.2f, 95%% CI [%.2f, %.2f]", s[1], s[2], s[3])
  x <- nnt_from_smd(smd = d, n_1 = n_1, n_2 = n_2)
  td <- generics::tidy(x)
  expect_identical(attr(x, "ci_source"), "smd")
  expect_identical(results_sentence(x), want)
  expect_identical(results_sentence(td), want)
  expect_identical(results_sentence(rbind(td, td)), want)
  expect_identical(results_sentence(rbind(td, td)[seq_len(nrow(td)), ]), want)
  expect_identical(results_sentence(td[-1, ]), want)
  bare <- td
  attr(bare, "ci_source") <- NULL
  expect_identical(results_sentence(bare), want)
  # tidy() gives the level on the rows of the SRD and the NNT and their
  # limits, and NA on the rows of d, as in DMAR 1.0.x.
  expect_identical(td$conf_level,
                   c(NA, 0.95, 0.95, NA, NA, 0.95, 0.95, 0.95, 0.95))
  # The glance() row holds the intervals of the SRD and the NNT, and asks
  # which one (as in DMAR 1.0.x); the interval for d is not counted.
  gl <- generics::glance(x)
  expect_error(results_sentence(gl, conf_level = 0.95),
               "several estimates (\"nnt\", \"srd\")", fixed = TRUE,
               class = "dmar_argument_invalid")
  expect_identical(results_sentence(gl, term = "srd", conf_level = 0.95),
                   want)
  # Supplied limits at a stated level, and with no level: the table then
  # records none, its tidy() column is NA throughout, and the sentence at
  # the level the caller gives is still the SRD's.
  s <- srd_of(c(0.3, -0.1, 0.7))
  want90 <- sprintf("srd = %.2f, 90%% CI [%.2f, %.2f]", s[1], s[2], s[3])
  y <- nnt_from_smd(smd = 0.3, smd_ci_lower = -0.1, smd_ci_upper = 0.7,
                    conf_level = 0.90)
  expect_identical(results_sentence(y), want90)
  expect_identical(results_sentence(generics::tidy(y)), want90)
  z <- nnt_from_smd(smd = 0.3, smd_ci_lower = -0.1, smd_ci_upper = 0.7)
  expect_true(all(is.na(generics::tidy(z)$conf_level)))
  expect_identical(results_sentence(generics::tidy(z), conf_level = 0.90),
                   want90)
  # The glance() row of an interval for d that contains 0 has no NNT
  # interval, so the SRD's is the one it reports.
  w <- nnt_from_smd(smd = 0.3, smd_ci_lower = -0.1, smd_ci_upper = 0.7,
                    conf_level = 0.90)
  expect_identical(results_sentence(generics::glance(w), conf_level = 0.90),
                   want90)
})
