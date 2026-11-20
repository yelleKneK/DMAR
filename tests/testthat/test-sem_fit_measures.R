sem_value <- function(x, term) x$value[match(term, x$term)]

fit_terms <- c("chi_square", "df", "p_chi_square",
               "rmsea", "rmsea_ci_lower", "rmsea_ci_upper",
               "close_fit_bound", "p_close_fit",
               "not_close_fit_bound", "p_not_close_fit",
               "cfi", "tli", "nfi", "ifi", "rni", "srmr",
               "gamma_hat", "gamma_hat_ci_lower", "gamma_hat_ci_upper",
               "adjusted_gamma_hat", "adjusted_gamma_hat_ci_lower",
               "adjusted_gamma_hat_ci_upper", "N", "n_groups")
gamma_terms <- c("gamma_hat", "gamma_hat_ci_lower", "gamma_hat_ci_upper",
                 "adjusted_gamma_hat", "adjusted_gamma_hat_ci_lower",
                 "adjusted_gamma_hat_ci_upper")

# The statistics lavaan reports for the three-factor model of the nine
# Holzinger and Swineford tests (N = 301), formed with N.
hs_T <- 85.30552; hs_df <- 24; hs_N <- 301; hs_Tb <- 918.85159; hs_dfb <- 36

# The noncentrality that puts T at the given lower-tail probability, by
# uniroot on pchisq directly, independent of ci_nc_chisq().
ncp_at <- function(T, df, prob) {
  if (stats::pchisq(T, df) <= prob) return(0)
  stats::uniroot(function(l) stats::pchisq(T, df, ncp = l) - prob,
                 c(0, 10 * T + 100), tol = 1e-12)$root
}

test_that("the table has the same rows in the same order whatever the input", {
  a <- sem_fit_measures(chisq = 30, df = 20, N = 200)
  b <- sem_fit_measures(chisq = 30, df = 20, N = 200, baseline_chisq = 500,
                        baseline_df = 28, p = 8)
  c2 <- sem_fit_measures(chisq = 60, df = 40, N = 300, n_groups = 2)
  for (x in list(a, b, c2)) {
    expect_s3_class(x, c("dmar_sem_fit", "dmar_tbl", "data.frame"))
    expect_identical(x$term, fit_terms)
    expect_type(x$value, "double")
  }
  expect_true(all(is.na(sem_value(a, c("cfi", "tli", "nfi", "ifi", "rni",
                                       "srmr", gamma_terms)))))
  expect_false(anyNA(sem_value(b, c("cfi", "tli", "nfi", "ifi", "rni",
                                    gamma_terms))))
  expect_equal(attr(a, "conf_level"), 0.90)
  expect_equal(sem_value(c2, c("N", "n_groups")), c(300, 2))
})

test_that("the default level is 0.90, the RMSEA's conventional level, and the table states it", {
  # A 90 percent interval decides the one-sided tests of close and not
  # close fit at .05 (MacCallum, Browne, & Sugawara, 1996, p. 136), the
  # tests the table reports beside it; the gamma hat intervals come from
  # the same noncentrality limits and share the level.
  expect_equal(eval(formals(sem_fit_measures)$conf_level), 0.90)
  x <- sem_fit_measures(chisq = hs_T, df = hs_df, N = hs_N,
                        baseline_chisq = hs_Tb, baseline_df = hs_dfb, p = 9,
                        sample_size_convention = "N")
  y <- sem_fit_measures(chisq = hs_T, df = hs_df, N = hs_N,
                        baseline_chisq = hs_Tb, baseline_df = hs_dfb, p = 9,
                        conf_level = 0.90, sample_size_convention = "N")
  z <- sem_fit_measures(chisq = hs_T, df = hs_df, N = hs_N,
                        baseline_chisq = hs_Tb, baseline_df = hs_dfb, p = 9,
                        conf_level = 0.95, sample_size_convention = "N")
  expect_identical(x$value, y$value)
  expect_true("Confidence level: 90%" %in% utils::capture.output(print(x)))
  expect_true("Confidence level: 95%" %in% utils::capture.output(print(z)))
  # The lower limit exceeds the close fit bound of .05 exactly when the test
  # of close fit rejects at .05, at the default level.
  expect_identical(sem_value(x, "rmsea_ci_lower") > 0.05,
                   sem_value(x, "p_close_fit") < 0.05)
  for (g in c("gamma_hat_ci_lower", "adjusted_gamma_hat_ci_upper")) {
    expect_false(isTRUE(all.equal(sem_value(x, g), sem_value(z, g))))
  }
  # Where the two levels disagree: on 24 df with N = 301, close fit is
  # rejected at .05 (p = .040) while the 95 percent interval still holds
  # .05; the 90 percent interval excludes it, as the test does.
  w90 <- sem_fit_measures(chisq = 62.85, df = 24, N = 301)
  w95 <- sem_fit_measures(chisq = 62.85, df = 24, N = 301, conf_level = 0.95)
  expect_lt(sem_value(w90, "p_close_fit"), 0.05)
  expect_gt(sem_value(w90, "rmsea_ci_lower"), 0.05)
  expect_lt(sem_value(w95, "rmsea_ci_lower"), 0.05)
})

test_that("summary input reproduces the closed forms under both conventions", {
  for (conv in c("N-1", "N")) {
    a <- if (conv == "N") hs_N else hs_N - 1
    x <- sem_fit_measures(chisq = hs_T, df = hs_df, N = hs_N,
                          baseline_chisq = hs_Tb, baseline_df = hs_dfb, p = 9,
                          sample_size_convention = conv)
    expect_equal(sem_value(x, "p_chi_square"),
                 stats::pchisq(hs_T, hs_df, lower.tail = FALSE))
    expect_equal(sem_value(x, "rmsea"), sqrt((hs_T - hs_df) / (a * hs_df)))
    expect_equal(sem_value(x, "cfi"), 1 - (hs_T - hs_df) / (hs_Tb - hs_dfb))
    expect_equal(sem_value(x, "rni"), 1 - (hs_T - hs_df) / (hs_Tb - hs_dfb))
    expect_equal(sem_value(x, "tli"),
                 (hs_Tb / hs_dfb - hs_T / hs_df) / (hs_Tb / hs_dfb - 1))
    expect_equal(sem_value(x, "nfi"), (hs_Tb - hs_T) / hs_Tb)
    expect_equal(sem_value(x, "ifi"), (hs_Tb - hs_T) / (hs_Tb - hs_df))
    g <- 9 / (9 + 2 * (hs_T - hs_df) / a)
    expect_equal(sem_value(x, "gamma_hat"), g)
    expect_equal(sem_value(x, "adjusted_gamma_hat"),
                 1 - 9 * 10 / (2 * hs_df) * (1 - g))
    expect_equal(attr(x, "sem_fit")$multiplier, a)
    expect_equal(attr(x, "sem_fit")$sample_size_convention, conv)
  }
  # The default is N - 1, the convention of ci_rmsea() and Kelley and Lai (2011).
  y <- sem_fit_measures(chisq = 30, df = 20, N = 200)
  expect_equal(sem_value(y, "rmsea"), sqrt(10 / (199 * 20)))
})

test_that("the RMSEA limits invert the noncentral chi square at the requested level", {
  for (level in c(0.95, 0.90)) {
    x <- sem_fit_measures(chisq = hs_T, df = hs_df, N = hs_N, conf_level = level,
                          sample_size_convention = "N")
    lo <- ncp_at(hs_T, hs_df, 1 - (1 - level) / 2)
    hi <- ncp_at(hs_T, hs_df, (1 - level) / 2)
    expect_equal(sem_value(x, "rmsea_ci_lower"), sqrt(lo / (hs_N * hs_df)),
                 tolerance = 1e-7)
    expect_equal(sem_value(x, "rmsea_ci_upper"), sqrt(hi / (hs_N * hs_df)),
                 tolerance = 1e-7)
    expect_equal(attr(x, "conf_level"), level)
  }
  # A statistic below df: the point estimate and the lower limit are 0,
  # while the upper limit of the 95 percent interval remains positive (at
  # 90 percent it is 0; see the next test).
  y <- sem_fit_measures(chisq = 10, df = 20, N = 200, conf_level = 0.95,
                        sample_size_convention = "N")
  expect_equal(sem_value(y, c("rmsea", "rmsea_ci_lower")), c(0, 0))
  expect_equal(sem_value(y, "rmsea_ci_upper"),
               sqrt(ncp_at(10, 20, 0.025) / (200 * 20)), tolerance = 1e-7)
  expect_length(attr(y, "sem_fit")$boundary, 0)
})

test_that("a statistic in the lower tail has an upper limit of 0, without a warning", {
  # pchisq(10, 20) is .032: above .025, below .05.
  expect_no_warning(x <- sem_fit_measures(chisq = 10, df = 20, N = 200,
                                          conf_level = 0.90))
  expect_equal(sem_value(x, c("rmsea_ci_lower", "rmsea_ci_upper")), c(0, 0))
  expect_match(attr(x, "sem_fit")$boundary, "upper RMSEA limit is 0")
  # The upper limit and the test of not close fit agree: the test rejects
  # for every positive bound at .05.
  expect_lt(sem_value(x, "p_not_close_fit"), 0.05)
  e <- equivalence_rmsea(10, 20, 200, rmsea_bound = 0.05)
  expect_equal(sem_value(x, "p_not_close_fit"), sem_value(e, "p_value"))
  expect_equal(sem_value(e, "ci_upper"), 0)
})

test_that("the tests of close and not close fit follow MacCallum et al. (1996)", {
  x <- sem_fit_measures(chisq = hs_T, df = hs_df, N = hs_N, close_fit = 0.05,
                        not_close_fit = 0.08, sample_size_convention = "N")
  expect_equal(sem_value(x, "p_close_fit"),
               stats::pchisq(hs_T, hs_df, ncp = hs_N * hs_df * 0.05^2,
                             lower.tail = FALSE))
  expect_equal(sem_value(x, "p_not_close_fit"),
               stats::pchisq(hs_T, hs_df, ncp = hs_N * hs_df * 0.08^2))
  # lavaan's values for this model at its bounds of .05 and .08.
  expect_equal(sem_value(x, "p_close_fit"), 0.0006612368, tolerance = 1e-6)
  expect_equal(sem_value(x, "p_not_close_fit"), 0.8395529, tolerance = 1e-6)
  expect_equal(sem_value(x, c("close_fit_bound", "not_close_fit_bound")),
               c(0.05, 0.08))
  # The test of not close fit is the equivalence test of equivalence_rmsea().
  e <- equivalence_rmsea(hs_T, hs_df, hs_N, rmsea_bound = 0.08,
                         sample_size_convention = "N")
  expect_equal(sem_value(x, "p_not_close_fit"), sem_value(e, "p_value"))
  # The defaults are .05 for both bounds (MacCallum et al., 1996).
  y <- sem_fit_measures(chisq = hs_T, df = hs_df, N = hs_N)
  expect_equal(sem_value(y, c("close_fit_bound", "not_close_fit_bound")),
               c(0.05, 0.05))
})

test_that("a multiple group statistic gets Steiger's multiple group RMSEA", {
  # lavaan, the two schools: chisq 115.8513 on 48 df, RMSEA 0.09691486, 95
  # percent limits 0.06989664 and 0.12372288.
  x <- sem_fit_measures(chisq = 115.8513, df = 48, N = 301, n_groups = 2,
                        conf_level = 0.95, sample_size_convention = "N")
  expect_equal(sem_value(x, "rmsea"), sqrt(2) * sqrt((115.8513 - 48) / (301 * 48)))
  expect_equal(sem_value(x, "rmsea"), 0.09691486, tolerance = 1e-6)
  expect_equal(sem_value(x, "rmsea_ci_lower"), 0.06989664, tolerance = 1e-6)
  expect_equal(sem_value(x, "rmsea_ci_upper"), 0.12372288, tolerance = 1e-6)
  expect_equal(sem_value(x, "p_not_close_fit"),
               stats::pchisq(115.8513, 48, ncp = 301 * 48 * 0.05^2 / 2))
  # Under N - 1 each group loses one: the multiplier is N - G.
  y <- sem_fit_measures(chisq = 115.8513, df = 48, N = 301, n_groups = 2)
  expect_equal(attr(y, "sem_fit")$multiplier, 299)
  expect_equal(sem_value(y, "rmsea"), sqrt(2 * (115.8513 - 48) / (299 * 48)))
  # The incremental indices do not depend on the number of groups.
  z1 <- sem_fit_measures(chisq = 115.8513, df = 48, N = 301, n_groups = 2,
                         baseline_chisq = 957.769, baseline_df = 72)
  z2 <- sem_fit_measures(chisq = 115.8513, df = 48, N = 301,
                         baseline_chisq = 957.769, baseline_df = 72)
  expect_equal(sem_value(z1, c("cfi", "tli", "nfi", "ifi", "rni")),
               sem_value(z2, c("cfi", "tli", "nfi", "ifi", "rni")))
})

test_that("zero degrees of freedom uses saturated conventions with explicit unavailable tests", {
  expect_warning(x <- sem_fit_measures(chisq = 0, df = 0, N = 100,
                                       baseline_chisq = 50, baseline_df = 3),
                 "zero degrees of freedom")
  expect_true(all(is.na(sem_value(x, c("p_chi_square", "p_close_fit", "p_not_close_fit")))))
  expect_equal(sem_value(x, c("cfi", "nfi", "rni")), c(1, 1, 1))
  expect_match(attr(x, "sem_fit")$boundary, "zero degrees of freedom")
})

test_that("a degenerate baseline leaves the affected indices NA, with a warning", {
  expect_warning(z <- sem_fit_measures(chisq = 20, df = 10, N = 100,
                                       baseline_chisq = 20, baseline_df = 20),
                 "tli, rni")
  expect_true(all(is.na(sem_value(z, c("tli", "rni")))))
  expect_equal(sem_value(z, "nfi"), 0)
})

test_that("incoherent summary input is refused with a message that says why", {
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 200, baseline_chisq = 500,
                                baseline_df = 0), "baseline_df")
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 200, baseline_chisq = 500,
                                baseline_df = 10), "nested")
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 200, p = 3), "moments")
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 1.5), "integer")
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 1.5,
                                sample_size_convention = "N"), "integer")
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 200.5), "integer")
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 2, n_groups = 2), "exceed")
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 200, n_groups = 2, p = 5),
               "single")
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 200, baseline_chisq = 500),
               "both")
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 200, version = "robust"),
               "standard measures only")
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 200, version = "scaled"),
               "standard measures only")
  expect_error(sem_fit_measures(chisq = 30, df = 20), "chisq")
  expect_error(sem_fit_measures(chisq = 30, df = 20.5, N = 200), "integer")
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 200, conf_level = 1),
               "below 1")
  expect_error(sem_fit_measures(chisq = 30, df = 20, N = 200, close_fit = 0),
               "greater")
})

test_that("named inputs, such as values from fitMeasures(), keep the row names", {
  x <- sem_fit_measures(chisq = c(chisq = 30), df = c(df = 20), N = c(ntotal = 200),
                        baseline_chisq = c(baseline.chisq = 500),
                        baseline_df = c(baseline.df = 28))
  expect_identical(x$term, fit_terms)
  expect_null(names(x$value))
})

test_that("the printed table names the version and the sample size convention", {
  x <- sem_fit_measures(chisq = hs_T, df = hs_df, N = hs_N)
  out <- utils::capture.output(print(x))
  expect_true(any(grepl("Version: standard", out)))
  expect_true(any(grepl("multiplier a = 300", out)))
})

hs_model <- 'visual =~ x1 + x2 + x3
             textual =~ x4 + x5 + x6
             speed =~ x7 + x8 + x9'

test_that("a fitted ML model reports lavaan's own values at the requested level", {
  skip_if_not_installed("lavaan")
  f <- lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939, std.lv = TRUE)
  for (level in c(0.95, 0.90)) {
    fm <- lavaan::fitMeasures(f, fm.args = list(rmsea.ci.level = level))
    x <- sem_fit_measures(f, conf_level = level)
    expect_identical(x$term, fit_terms)
    expect_equal(sem_value(x, c("chi_square", "df", "p_chi_square", "rmsea",
                                "rmsea_ci_lower", "rmsea_ci_upper", "p_close_fit",
                                "cfi", "tli", "nfi", "ifi", "rni")),
                 unname(fm[c("chisq", "df", "pvalue", "rmsea", "rmsea.ci.lower",
                             "rmsea.ci.upper", "rmsea.pvalue", "cfi", "tli",
                             "nfi", "ifi", "rni")]))
    expect_equal(sem_value(x, "srmr"), unname(fm[["srmr_bentler_nomean"]]))
    expect_equal(sem_value(x, c("N", "n_groups")), c(301, 1))
  }
  # Gamma hat against its closed form with lavaan's N.
  fm <- lavaan::fitMeasures(f)
  g <- 9 / (9 + 2 * (fm[["chisq"]] - 24) / 301)
  expect_equal(sem_value(x, "gamma_hat"), g)
  expect_equal(sem_value(x, "adjusted_gamma_hat"), 1 - 90 / 48 * (1 - g))
  details <- attr(x, "sem_fit")
  expect_equal(details$version, "standard")
  expect_equal(details$sample_size_convention, "N")
  expect_equal(details$multiplier, 301)
  # Summary input with lavaan's statistics and its convention gives the
  # same table, but for the SRMR, which needs the residuals.
  y <- sem_fit_measures(chisq = fm[["chisq"]], df = fm[["df"]], N = 301,
                        baseline_chisq = fm[["baseline.chisq"]],
                        baseline_df = fm[["baseline.df"]], p = 9,
                        conf_level = 0.90, sample_size_convention = "N")
  keep <- setdiff(fit_terms, "srmr")
  expect_equal(sem_value(y, keep), sem_value(x, keep), tolerance = 1e-7)
  # The test of not close fit at lavaan's bound of .08.
  z <- sem_fit_measures(f, not_close_fit = 0.08)
  expect_equal(sem_value(z, "p_not_close_fit"),
               unname(lavaan::fitMeasures(f, "rmsea.notclose.pvalue")))
})

test_that("semTools agrees on gamma hat and adjusted gamma hat", {
  skip_if_not_installed("lavaan")
  f <- lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939, std.lv = TRUE)
  x <- sem_fit_measures(f)
  # Pinned from semTools::moreFitIndices(f, c("gammaHat", "adjGammaHat"))
  # (semTools 0.5.9, 2026-11-01); live comparison in
  # tools/oracle_snippets/semtools.R.
  expect_equal(sem_value(x, c("gamma_hat", "adjusted_gamma_hat")),
               c(0.95669919291192251, 0.91881098670985473), tolerance = 1e-10)
})

test_that("a robust ML fit reports the robust measures unless told otherwise", {
  skip_if_not_installed("lavaan")
  f <- lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939, std.lv = TRUE,
                   estimator = "MLR")
  fm <- lavaan::fitMeasures(f, fm.args = list(rmsea.ci.level = 0.90))
  x <- sem_fit_measures(f, conf_level = 0.90)
  expect_equal(attr(x, "sem_fit")$version, "robust")
  expect_match(attr(x, "sem_fit")$version_reason, "yuan.bentler.mplus")
  expect_equal(sem_value(x, c("chi_square", "df", "p_chi_square", "rmsea",
                              "rmsea_ci_lower", "rmsea_ci_upper", "p_close_fit",
                              "cfi", "tli", "rni")),
               unname(fm[c("chisq.scaled", "df.scaled", "pvalue.scaled",
                           "rmsea.robust", "rmsea.ci.lower.robust",
                           "rmsea.ci.upper.robust", "rmsea.pvalue.robust",
                           "cfi.robust", "tli.robust", "rni.robust")]))
  expect_true(all(is.na(sem_value(x, c("nfi", "ifi")))))
  # Gamma hat has a robust version: lavaan's gfi_lrt.robust and its limits.
  # (lavaan's gamma hat interval is at its default level, 90 percent.)
  fg <- lavaan::fitMeasures(f, c("gfi_lrt.robust", "gfi_lrt.ci.lower.robust",
                                 "gfi_lrt.ci.upper.robust"))
  expect_equal(sem_value(x, c("gamma_hat", "gamma_hat_ci_lower",
                              "gamma_hat_ci_upper")), unname(fg),
               tolerance = 1e-7)
  expect_match(attr(x, "sem_fit")$correction, "Brosseau-Liard")
  expect_true(any(grepl("Version: robust", utils::capture.output(print(x)))))
  s <- sem_fit_measures(f, conf_level = 0.90, version = "scaled")
  expect_equal(sem_value(s, c("chi_square", "rmsea", "cfi", "nfi")),
               unname(fm[c("chisq.scaled", "rmsea.scaled", "cfi.scaled",
                           "nfi.scaled")]))
  expect_true(all(is.na(sem_value(s, gamma_terms))))
  st <- sem_fit_measures(f, conf_level = 0.90, version = "standard")
  expect_equal(sem_value(st, c("chi_square", "rmsea", "cfi")),
               unname(fm[c("chisq", "rmsea", "cfi")]))
  expect_false(is.na(sem_value(st, "gamma_hat")))
  # The SRMR has no test statistic in it and is the same in every version.
  expect_equal(sem_value(x, "srmr"), sem_value(st, "srmr"))
  # A plain ML fit has no robust or scaled measures.
  g <- lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939, std.lv = TRUE)
  expect_error(sem_fit_measures(g, version = "robust"), "corrected test statistic")
  expect_error(sem_fit_measures(g, version = "scaled"), "no scaled test statistic")
})

test_that("a FIML fit reports the FIML corrected measures of Zhang and Savalei (2023)", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  set.seed(113)
  d$x2[sample(nrow(d), 40)] <- NA
  d$x9[sample(nrow(d), 30)] <- NA
  f <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, missing = "ml")
  fm <- lavaan::fitMeasures(f)
  expect_no_warning(x <- sem_fit_measures(f, conf_level = 0.90))
  expect_equal(attr(x, "sem_fit")$version, "robust")
  expect_match(attr(x, "sem_fit")$correction, "Zhang and Savalei.*version 6")
  # No scaled statistic: the chi square row is the standard FIML statistic.
  expect_equal(sem_value(x, c("chi_square", "rmsea", "cfi", "tli")),
               unname(fm[c("chisq", "rmsea.robust", "cfi.robust", "tli.robust")]))
  expect_error(sem_fit_measures(f, version = "scaled"), "no scaled")
  expect_no_warning(sem_fit_measures(f, version = "standard"))
  fr <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, missing = "ml",
                    estimator = "MLR")
  fmr <- lavaan::fitMeasures(fr)
  y <- sem_fit_measures(fr)
  expect_match(attr(y, "sem_fit")$correction, "version 3")
  expect_equal(sem_value(y, c("chi_square", "rmsea", "cfi")),
               unname(fmr[c("chisq.scaled", "rmsea.robust", "cfi.robust")]))
  # The SRMR keeps the covariance definition, which lavaan's srmr does not
  # under a mean structure.
  expect_equal(sem_value(y, "srmr"), unname(fmr[["srmr_bentler_nomean"]]))
})

test_that("ordered data under WLSMV report Savalei's (2021) corrected measures", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  d <- lavaan::HolzingerSwineford1939
  for (v in paste0("x", 1:9)) {
    d[[v]] <- cut(d[[v]], stats::quantile(d[[v]], c(0, .25, .5, .75, 1)),
                  include.lowest = TRUE, labels = FALSE)
  }
  f <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, ordered = TRUE)
  fm <- lavaan::fitMeasures(f)
  x <- sem_fit_measures(f)
  expect_equal(attr(x, "sem_fit")$version, "robust")
  expect_match(attr(x, "sem_fit")$correction, "Savalei \\(2021\\)")
  expect_false(is.na(sem_value(x, "p_chi_square")))
  expect_equal(sem_value(x, c("chi_square", "rmsea", "cfi")),
               unname(fm[c("chisq.scaled", "rmsea.robust", "cfi.robust")]))
  # The least squares statistic is formed with N - 1.
  expect_equal(attr(x, "sem_fit")$sample_size_convention, "N-1")
  expect_warning(sem_fit_measures(f, version = "standard"), "not chi square")
  expect_warning(sem_fit_measures(f, version = "scaled"), "Savalei")
  # A fit by DWLS has only the standard statistic, lavaan's default test
  # for that estimator, which lavaan gives no reference distribution:
  # "auto" reports the standard measures, says so from lavaan's options
  # (not as a test the user set), and names the remedy.
  g <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, ordered = TRUE,
                   estimator = "DWLS")
  expect_identical(lavaan::lavInspect(g, "options")$test, "standard")
  expect_warning(y <- sem_fit_measures(g), "refit with a corrected test")
  expect_identical(attr(y, "sem_fit")$version_reason, paste0(
    "a fit to ordered data by DWLS with the standard test alone, so it has ",
    "no corrected statistic (estimator = \"WLSMV\" is DWLS with one)"))
  expect_true(is.na(sem_value(y, "p_chi_square")))
  expect_error(sem_fit_measures(g, version = "robust"),
               "not chi square distributed. Refit it with a corrected test")
  # Full WLS is the exception (the 1.1.0 review): its statistic is
  # asymptotically chi square (Browne, 1984; Muthen, 1984) and lavaan
  # reports its p-value, which the warning had called NA beside a table
  # that printed it, and the refusal of the robust version had called not
  # chi square distributed.
  w <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, ordered = TRUE,
                   estimator = "WLS")
  expect_no_warning(z <- sem_fit_measures(w))
  expect_equal(attr(z, "sem_fit")$version, "standard")
  expect_identical(attr(z, "sem_fit")$version_reason, paste0(
    "a fit to ordered data by full WLS, whose statistic is asymptotically ",
    "chi square without a correction"))
  expect_equal(sem_value(z, c("chi_square", "p_chi_square")),
               unname(lavaan::fitMeasures(w, c("chisq", "pvalue"))))
  expect_false(is.na(sem_value(z, "p_chi_square")))
  expect_error(sem_fit_measures(w, version = "robust"),
               "that of full WLS, is asymptotically chi square")
  wr <- tryCatch(sem_fit_measures(w, version = "robust"),
                 error = conditionMessage)
  expect_false(grepl("not chi square distributed", wr, fixed = TRUE))
})

test_that("the reason for the standard version names the estimator of a continuous fit", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939
  reason <- function(est) {
    attr(sem_fit_measures(lavaan::cfa(hs_model, data = hs, std.lv = TRUE,
                                      estimator = est)), "sem_fit")$version_reason
  }
  expect_identical(reason("ML"),
                   "a normal theory ML fit of complete data with no corrected statistic")
  # Full WLS of continuous data is the asymptotically distribution free
  # estimator, which the reason had called a normal theory fit.
  expect_identical(reason("WLS"), paste0(
    "an asymptotically distribution free (WLS) fit of complete data, whose ",
    "statistic needs no correction"))
  expect_identical(reason("ULS"), paste0(
    "a ULS fit of complete data with no corrected statistic, tested by ",
    "Browne's residual based statistic"))
})

test_that("a least squares fit of continuous data with the standard test alone is not treated as chi square", {
  skip_if_not_installed("lavaan")
  # lavaan's lav_model_test() gives the standard statistic of DWLS and ULS
  # no reference distribution whatever the data; with continuous data the
  # default test adds Browne's residual based statistic, which is chi
  # square under normality and is then the one reported, so only a test
  # set to "standard" leaves the fit without one. The 1.1.0 review found
  # the continuous DWLS fit reported with an RMSEA interval and close fit
  # tests and no warning, and the continuous ULS fit's warning offering a
  # scaled version it does not have.
  hs <- lavaan::HolzingerSwineford1939
  w <- lavaan::cfa(hs_model, data = hs, std.lv = TRUE, estimator = "DWLS",
                   ordered = FALSE, test = "standard")
  expect_true(is.na(lavaan::fitMeasures(w, "pvalue")))
  expect_warning(x <- sem_fit_measures(w),
                 "standard statistic of a DWLS fit is not chi square")
  expect_identical(attr(x, "sem_fit")$version_reason, paste0(
    "a DWLS fit of complete data with the standard test alone, whose ",
    "statistic lavaan gives no reference distribution"))
  expect_true(is.na(sem_value(x, "p_chi_square")))
  expect_error(sem_fit_measures(w, version = "robust"),
               "that of DWLS, is not chi square distributed")
  wr <- tryCatch(sem_fit_measures(w, version = "robust"),
                 error = conditionMessage)
  expect_false(grepl("the ones to report", wr, fixed = TRUE))
  u <- lavaan::cfa(hs_model, data = hs, std.lv = TRUE, estimator = "ULS",
                   test = "standard")
  expect_error(sem_fit_measures(u, version = "robust"),
               "ML-derived robust fit-index correction is not supported")
  # The default test of a DWLS fit of continuous data is Browne's, and
  # nothing is said.
  wb <- lavaan::cfa(hs_model, data = hs, std.lv = TRUE, estimator = "DWLS",
                    ordered = FALSE)
  expect_no_warning(xb <- sem_fit_measures(wb))
  expect_match(attr(xb, "sem_fit")$version_reason,
               "tested by Browne's residual based statistic", fixed = TRUE)
  expect_false(is.na(sem_value(xb, "p_chi_square")))
})

test_that("a multiple group fit reports lavaan's multiple group RMSEA", {
  skip_if_not_installed("lavaan")
  f <- lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939,
                   std.lv = TRUE, group = "school")
  # lavaan's own default level, 0.90, is DMAR's.
  fm <- lavaan::fitMeasures(f)
  x <- sem_fit_measures(f)
  expect_equal(sem_value(x, c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper")),
               unname(fm[c("rmsea", "rmsea.ci.lower", "rmsea.ci.upper")]))
  expect_equal(sem_value(x, c("N", "n_groups")), c(301, 2))
  expect_true(all(is.na(sem_value(x, c("gamma_hat", "adjusted_gamma_hat")))))
  y <- sem_fit_measures(chisq = fm[["chisq"]], df = fm[["df"]], N = 301, n_groups = 2,
                        sample_size_convention = "N")
  expect_equal(sem_value(y, c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper")),
               sem_value(x, c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper")),
               tolerance = 1e-7)
})

test_that("Wishart and least squares fits keep lavaan's multiplier N - G", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  f <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, likelihood = "wishart")
  fm <- lavaan::fitMeasures(f)
  x <- sem_fit_measures(f)
  expect_equal(sem_value(x, "rmsea"), unname(fm[["rmsea"]]))
  expect_equal(sem_value(x, "rmsea"), sqrt((fm[["chisq"]] - 24) / (300 * 24)))
  expect_equal(attr(x, "sem_fit")$sample_size_convention, "N-1")
  g <- 9 / (9 + 2 * (fm[["chisq"]] - 24) / 300)
  expect_equal(sem_value(x, "gamma_hat"), g)
  # The RMSEA and gamma hat imply the same population discrepancy.
  expect_equal(sem_value(x, "rmsea")^2 * 24, (9 / g - 9) / 2, tolerance = 1e-10)
  # The same statistic as summary input under the default convention.
  y <- sem_fit_measures(chisq = fm[["chisq"]], df = 24, N = 301)
  expect_equal(sem_value(y, c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper")),
               sem_value(x, c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper")),
               tolerance = 1e-7)
  h <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, estimator = "GLS")
  z <- sem_fit_measures(h)
  expect_equal(attr(z, "sem_fit")$multiplier, 300)
  expect_equal(sem_value(z, "rmsea"), unname(lavaan::fitMeasures(h, "rmsea")))
  expect_true(is.na(sem_value(z, "gamma_hat")))
})

test_that("saturated fitted and summary input use the same conventions", {
  skip_if_not_installed("lavaan")
  f <- lavaan::cfa("f =~ x1 + x2 + x3", data = lavaan::HolzingerSwineford1939,
                   std.lv = TRUE)
  fm <- lavaan::fitMeasures(f)
  expect_warning(x <- sem_fit_measures(f), "zero degrees of freedom")
  expect_warning(y <- sem_fit_measures(chisq = fm[["chisq"]], df = 0, N = 301,
                                       baseline_chisq = fm[["baseline.chisq"]],
                                       baseline_df = fm[["baseline.df"]],
                                       sample_size_convention = "N"),
                 "zero degrees of freedom")
  na_rows <- c("p_chi_square", "p_close_fit", "p_not_close_fit")
  expect_true(all(is.na(sem_value(x, na_rows))))
  expect_true(all(is.na(sem_value(y, na_rows))))
  expect_equal(sem_value(x, "cfi"), sem_value(y, "cfi"))
  expect_match(attr(x, "sem_fit")$boundary, "zero degrees of freedom")
})

test_that("a fitted statistic in the lower tail gets an upper limit of 0 in both paths", {
  skip_if_not_installed("lavaan")
  # A correctly specified one-factor model whose statistic lands in the lower
  # tail (T = 1.74 on 9 df, pchisq = .005).
  L <- rep(0.7, 6)
  R <- chol(tcrossprod(L) + diag(1 - L^2))
  set.seed(65)
  X <- matrix(stats::rnorm(100 * 6), 100) %*% R
  colnames(X) <- paste0("y", 1:6)
  f <- lavaan::cfa("f =~ y1 + y2 + y3 + y4 + y5 + y6", data = as.data.frame(X),
                   std.lv = TRUE)
  T <- unname(lavaan::fitMeasures(f, "chisq"))
  expect_lt(stats::pchisq(T, 9), 0.025)
  for (level in c(0.90, 0.95)) {
    x <- sem_fit_measures(f, conf_level = level)
    y <- sem_fit_measures(chisq = T, df = 9, N = 100, conf_level = level,
                          sample_size_convention = "N")
    expect_equal(sem_value(x, c("rmsea_ci_lower", "rmsea_ci_upper")), c(0, 0))
    expect_equal(sem_value(y, c("rmsea_ci_lower", "rmsea_ci_upper")), c(0, 0))
    expect_match(attr(x, "sem_fit")$boundary, "upper RMSEA limit is 0")
    expect_match(attr(y, "sem_fit")$boundary, "upper RMSEA limit is 0")
  }
})

test_that("unconverged, improper, multilevel, and doubly specified fits are refused", {
  skip_if_not_installed("lavaan")
  P <- diag(5)
  P[1, 2] <- P[2, 1] <- P[1, 3] <- P[3, 1] <- 0.8
  P[2, 3] <- P[3, 2] <- 0.3
  P[4, 5] <- P[5, 4] <- 0.5
  dimnames(P) <- list(paste0("y", 1:5), paste0("y", 1:5))
  heywood <- suppressWarnings(lavaan::cfa("f1 =~ y1 + y2 + y3\n f2 =~ y4 + y5",
                                          sample.cov = P, sample.nobs = 300))
  # lavaan's own check warns as it finds the negative variance.
  expect_error(suppressWarnings(sem_fit_measures(heywood)), "inadmissible")
  # The refusal is classed and names what failed.
  e <- tryCatch(suppressWarnings(sem_fit_measures(heywood)), dmar_refusal = function(e) e)
  expect_s3_class(e, c("dmar_inadmissible", "dmar_refusal", "error"))
  expect_identical(e$code, "inadmissible")
  expect_identical(e$reason, "negative_variance")
  expect_identical(e$variables, "y1")
  expect_identical(e$auxiliary, character(0))
  expect_match(conditionMessage(e), paste0(
    "has an inadmissible solution (lavaan's post-estimation check finds a ",
    "negative variance estimate for y1); resolve it before evaluating fit."),
    fixed = TRUE)
  stalled <- suppressWarnings(lavaan::cfa(hs_model,
                                          data = lavaan::HolzingerSwineford1939,
                                          control = list(iter.max = 1)))
  expect_error(sem_fit_measures(stalled), "did not converge")
  expect_error(sem_fit_measures(stalled), class = "dmar_not_converged")
  expect_error(sem_fit_measures(stats::lm(x1 ~ x2,
                                          data = lavaan::HolzingerSwineford1939)),
               "lavaan model")
  f <- lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939)
  expect_error(sem_fit_measures(f, N = 301), "not both")
  expect_error(sem_fit_measures(f, n_groups = 1), "'n_groups'")
  expect_error(sem_fit_measures(f, sample_size_convention = "N"),
               "sample_size_convention")
  skip_on_cran()
  two_level <- lavaan::sem("level: 1
                              fw =~ y1 + y2 + y3
                            level: 2
                              fb =~ y1 + y2 + y3",
                           data = lavaan::Demo.twolevel, cluster = "cluster")
  expect_error(sem_fit_measures(two_level), "multilevel")
})

test_that("the check of admissibility is lavaan's, and its refusal names what failed", {
  skip_if_not_installed("lavaan")
  # Without auxiliary variables the verdict is lavaan's post-estimation
  # check, and the conditions found here are what that check flags, in its
  # order: the variances, the latent covariance matrix, theta, and psi.
  P <- diag(5)
  P[1, 2] <- P[2, 1] <- P[1, 3] <- P[3, 1] <- 0.8
  P[2, 3] <- P[3, 2] <- 0.3
  P[4, 5] <- P[5, 4] <- 0.5
  dimnames(P) <- list(paste0("y", 1:5), paste0("y", 1:5))
  P2 <- diag(5)
  P2[1:3, 1:3] <- 0.5
  P2[4, 5] <- P2[5, 4] <- 0.5
  diag(P2) <- 1
  dimnames(P2) <- dimnames(P)
  R <- diag(6)
  R[1:3, 1:3] <- R[4:6, 4:6] <- 0.5
  R[1:3, 4:6] <- R[4:6, 1:3] <- 0.55
  diag(R) <- 1
  dimnames(R) <- list(paste0("z", 1:6), paste0("z", 1:6))
  # Three items that correlate .6 with an observed variable and .25 with
  # each other: the factor's correlation with it is 1.2.
  Q <- matrix(0.25, 4, 4)
  Q[4, ] <- Q[, 4] <- 0.6
  diag(Q) <- 1
  dimnames(Q) <- list(c("w1", "w2", "w3", "x"), c("w1", "w2", "w3", "x"))
  quiet <- function(expr) suppressWarnings(expr)
  fits <- list(
    heywood = quiet(lavaan::cfa("f1 =~ y1 + y2 + y3\n f2 =~ y4 + y5", sample.cov = P,
                                sample.nobs = 300)),
    factor_correlation = quiet(lavaan::cfa("a =~ z1 + z2 + z3\n b =~ z4 + z5 + z6",
                                           sample.cov = R, sample.nobs = 200,
                                           std.lv = TRUE)),
    latent_variance = quiet(lavaan::cfa("f =~ 0.5*z1 + 0.5*z2 + 0.5*z3\n f ~~ -0.2*f",
                                        sample.cov = R[1:3, 1:3], sample.nobs = 200)),
    structural = quiet(lavaan::cfa("f =~ w1 + w2 + w3\n f ~~ x", sample.cov = Q,
                                   sample.nobs = 300, std.lv = TRUE)),
    group = quiet(lavaan::cfa("f1 =~ y1 + y2 + y3\n f2 =~ y4 + y5",
                              sample.cov = list(P2, P), sample.nobs = c(300, 300))),
    proper = lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939),
    proper_groups = lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939,
                                group = "school")
  )
  for (f in fits) {
    expect_true(lavaan::lavInspect(f, "converged"))
    lavaan_ok <- isTRUE(suppressWarnings(lavaan::lavInspect(f, "post.check")))
    expect_identical(nrow(DMAR:::.sem_inadmissible_conditions(f)) == 0L, lavaan_ok)
    expect_identical(as.logical(suppressWarnings(DMAR:::.sem_admissible(f))), lavaan_ok)
  }
  reasons <- vapply(fits, function(f) {
    paste(DMAR:::.sem_inadmissible_conditions(f)$reason, collapse = " ")
  }, character(1))
  expect_identical(unname(reasons), c("negative_variance", "cov.lv", "negative_variance",
                                      "psi", "negative_variance", "", ""))
  refusal <- function(f) {
    tryCatch(suppressWarnings(sem_fit_measures(f)), dmar_refusal = function(e) e)
  }
  e <- refusal(fits$factor_correlation)
  expect_identical(e$reason, "cov.lv")
  expect_identical(e$variables, character(0))
  expect_match(conditionMessage(e), paste0("finds a covariance matrix of the latent ",
                                           "variables (cov.lv) that is not positive ",
                                           "definite)"), fixed = TRUE)
  expect_lt(e$problems$value, 0)
  e <- refusal(fits$latent_variance)
  expect_identical(e$variables, "f")
  expect_equal(e$problems$value, -0.2)
  e <- refusal(fits$structural)
  expect_identical(e$reason, "psi")
  expect_match(conditionMessage(e), "structural part of the model (psi)", fixed = TRUE)
  e <- refusal(fits$group)
  expect_identical(e$problems$group, "Group 2")
  expect_match(conditionMessage(e), "negative variance estimate for y1 in group Group 2",
               fixed = TRUE)
  # compare_sem() names the model at fault.
  e <- tryCatch(suppressWarnings(compare_sem(list(good = fits$proper,
                                                  improper = fits$heywood))),
                dmar_refusal = function(e) e)
  expect_s3_class(e, "dmar_inadmissible")
  # The test of a matrix is lavaan's, tolerance included.
  psd <- tryCatch(utils::getFromNamespace("lav_object_post_check_psd", "lavaan"),
                  error = function(e) NULL)
  skip_if(is.null(psd), "lavaan's own test of a matrix is not available")
  for (shift in c(-1e-3, -1e-9, -1e-13, 0, 1e-13)) {
    m <- diag(c(2, 1, shift))
    expect_identical(DMAR:::.sem_psd(m), psd(m))
  }
})

test_that("the refusal of an inadmissible solution places each failure in its group and at its level", {
  # The message names the group and, in a multilevel fit, the level.
  p <- data.frame(reason = c("negative_variance", "theta", "theta"),
                  group = c("A", "A", "B"), level = c(2L, 1L, 2L),
                  variable = c("y1", NA, NA), value = c(-0.05, -0.1, -0.2),
                  stringsAsFactors = FALSE)
  e <- tryCatch(DMAR:::.sem_refuse_inadmissible(p, character(0), "'fit'"),
                dmar_refusal = function(e) e)
  expect_s3_class(e, "dmar_inadmissible")
  expect_identical(e$reason, c("negative_variance", "theta"))
  expect_match(conditionMessage(e), paste0(
    "finds a negative variance estimate for y1 in group A at level 2; a ",
    "covariance matrix of the residuals of the observed variables (theta) ",
    "that is not positive definite in group A at level 1 and in group B at ",
    "level 2)"), fixed = TRUE)
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # A two-level fit reaches the check through sem_information_criteria(),
  # which takes it for AIC, and through sem_fit_measures(), whose check
  # comes before its refusal of a multilevel fit. lavaan names the second
  # level's block after the cluster variable; the refusal gives its number.
  m <- "level: 1
          fw =~ y1 + y2 + y3
        level: 2
          fb =~ y1 + y2 + y3
          y1 ~~ -0.05*y1"
  two_level <- suppressWarnings(lavaan::sem(m, data = lavaan::Demo.twolevel,
                                            cluster = "cluster"))
  e <- tryCatch(suppressWarnings(sem_information_criteria(two_level, measures = "AIC")),
                dmar_refusal = function(e) e)
  expect_s3_class(e, "dmar_inadmissible")
  expect_identical(e$problems$group, NA_character_)
  expect_identical(e$problems$level, 2L)
  expect_identical(e$variables, "y1")
  expect_match(conditionMessage(e), "negative variance estimate for y1 at level 2)",
               fixed = TRUE)
  expect_error(suppressWarnings(sem_fit_measures(two_level)), class = "dmar_inadmissible")
  # With two groups lavaan numbers the blocks level within group.
  d <- lavaan::Demo.twolevel
  d$grp <- ifelse(d$cluster %% 2 == 0, "even", "odd")
  m2 <- "group: odd
         level: 1
           fw =~ y1 + y2 + y3
         level: 2
           fb =~ y1 + y2 + y3
         group: even
         level: 1
           fw =~ y1 + y2 + y3
         level: 2
           fb =~ y1 + y2 + y3
           y1 ~~ -0.05*y1"
  two_groups <- suppressWarnings(lavaan::sem(m2, data = d, cluster = "cluster",
                                             group = "grp"))
  expect_identical(lavaan::lavInspect(two_groups, "group.label"), c("odd", "even"))
  e <- tryCatch(suppressWarnings(sem_information_criteria(two_groups, measures = "AIC")),
                dmar_refusal = function(e) e)
  expect_identical(e$problems$group, "even")
  expect_identical(e$problems$level, 2L)
  expect_match(conditionMessage(e), "for y1 in group even at level 2)", fixed = TRUE)
})

test_that("listwise deletion is named in a warning, and a weighted fit is robust", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  d$x1[1:20] <- NA
  expect_warning(x <- sem_fit_measures(lavaan::cfa(hs_model, data = d)),
                 "removed 20 of 301 cases by listwise deletion")
  expect_equal(sem_value(x, "N"), 281)
  w <- lavaan::HolzingerSwineford1939
  w$wt <- rep(c(0.5, 1.5), length.out = nrow(w))
  fw <- lavaan::cfa(hs_model, data = w, sampling.weights = "wt")
  y <- sem_fit_measures(fw)
  expect_equal(attr(y, "sem_fit")$version, "robust")
  expect_warning(sem_fit_measures(fw, version = "standard"), "sampling weighted")
})

test_that("gamma hat has the limits of the RMSEA interval, as in Steiger and Fouladi (1997)", {
  # Steiger and Fouladi (1997, p. 252): a circumplex model of six variables,
  # T = 27.05 on 12 df, n = 710. The 90 percent interval for the adjusted
  # population gamma is .972 to .997, and the RMSEA interval begins at .021.
  x <- sem_fit_measures(chisq = 27.05, df = 12, N = 710, p = 6, conf_level = 0.90)
  expect_equal(round(sem_value(x, c("adjusted_gamma_hat_ci_lower",
                                    "adjusted_gamma_hat_ci_upper")), 3),
               c(0.972, 0.997))
  expect_equal(round(sem_value(x, "rmsea_ci_lower"), 3), 0.021)
  # Each limit is p/(p + 2 df RMSEA^2) at the opposite RMSEA limit, and the
  # adjusted limits are the adjustment of those.
  lo <- sem_value(x, "rmsea_ci_lower"); hi <- sem_value(x, "rmsea_ci_upper")
  g <- 6 / (6 + 2 * 12 * c(hi, lo)^2)
  expect_equal(sem_value(x, c("gamma_hat_ci_lower", "gamma_hat_ci_upper")), g)
  expect_equal(sem_value(x, c("adjusted_gamma_hat_ci_lower",
                              "adjusted_gamma_hat_ci_upper")),
               1 - 42 / 24 * (1 - g))
  # A statistic below df: the estimate exceeds 1 and lies above its
  # interval, whose upper limit is 1.
  y <- sem_fit_measures(chisq = 10, df = 20, N = 200, p = 8)
  expect_gt(sem_value(y, "gamma_hat"), 1)
  expect_equal(sem_value(y, c("gamma_hat_ci_upper", "adjusted_gamma_hat_ci_upper")),
               c(1, 1))
})

test_that("a fitted model's gamma hat limits are lavaan's gfi_lrt limits", {
  skip_if_not_installed("lavaan")
  f <- lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939, std.lv = TRUE)
  x <- sem_fit_measures(f, conf_level = 0.90)
  fg <- lavaan::fitMeasures(f, c("gfi_lrt", "gfi_lrt.ci.lower", "gfi_lrt.ci.upper"))
  expect_equal(sem_value(x, c("gamma_hat", "gamma_hat_ci_lower",
                              "gamma_hat_ci_upper")), unname(fg),
               tolerance = 1e-7)
  # lavaan's default gfi is the same function of Browne's residual based
  # statistic, which differs from gamma hat.
  expect_false(isTRUE(all.equal(sem_value(x, "gamma_hat"),
                                unname(lavaan::fitMeasures(f, "gfi")))))
})

# Rows of the 90 percent interval and the tests of close and not close fit.
rmsea_terms <- c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper", "p_close_fit",
                 "p_not_close_fit")

test_that("the scaled version puts the scaled statistic in the standard formulas", {
  skip_if_not_installed("lavaan")
  # Heavy tailed data (multivariate t on 4 df), for which the Yuan-Bentler
  # scaling constant is far from 1: the defect showed only then. lavaan's
  # default scaled interval refers the unscaled statistic to c df degrees of
  # freedom, and the table once paired it with the scaled statistic on df.
  d <- lavaan::HolzingerSwineford1939[paste0("x", 1:9)]
  set.seed(1)
  X <- scale(as.matrix(d), scale = FALSE)
  d <- as.data.frame(sweep(X, 1, sqrt(4 / stats::rchisq(301, 4)), "*"))
  f <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, estimator = "MLR")
  s <- sem_fit_measures(f, version = "scaled", conf_level = 0.90)
  expect_gt(unname(lavaan::fitMeasures(f, "chisq.scaling.factor")), 1.5)
  # The same rows as summary input with the chi_square and df rows.
  plug <- sem_fit_measures(chisq = sem_value(s, "chi_square"),
                           df = sem_value(s, "df"), N = 301, conf_level = 0.90,
                           sample_size_convention = "N")
  expect_equal(sem_value(s, rmsea_terms), sem_value(plug, rmsea_terms))
  # lavaan's values under rmsea.scaled.mplus = TRUE, as Mplus reports them.
  g <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, estimator = "MLR",
                   rmsea.scaled.mplus = TRUE)
  fm <- lavaan::fitMeasures(g, fm.args = list(rmsea.ci.level = 0.90,
                                              rmsea.h0.notclosefit = 0.05))
  lav <- unname(fm[c("rmsea.scaled", "rmsea.ci.lower.scaled",
                     "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled",
                     "rmsea.notclose.pvalue.scaled")])
  expect_equal(sem_value(s, rmsea_terms), lav, tolerance = 1e-6)
  # The values differ from lavaan's default, whose interval begins at .056
  # where the one the rows imply begins at .051.
  expect_equal(sem_value(s, c("rmsea_ci_lower", "p_close_fit")), c(0.0506, 0.0459),
               tolerance = 1e-3)
  fd <- lavaan::fitMeasures(f, fm.args = list(rmsea.ci.level = 0.90))
  expect_gt(abs(sem_value(s, "rmsea_ci_lower") - fd[["rmsea.ci.lower.scaled"]]), 0.005)
  # The scaled CFI and the point estimate were already the substitution.
  expect_equal(sem_value(s, c("rmsea", "cfi")),
               unname(fd[c("rmsea.scaled", "cfi.scaled")]))
})

test_that("a scaled and shifted statistic gives lavaan's own scaled values", {
  skip_if_not_installed("lavaan")
  f <- lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939, std.lv = TRUE,
                   estimator = "MLMV")
  s <- sem_fit_measures(f, version = "scaled", conf_level = 0.90)
  fm <- lavaan::fitMeasures(f, fm.args = list(rmsea.ci.level = 0.90,
                                              rmsea.h0.notclosefit = 0.05))
  expect_equal(sem_value(s, rmsea_terms),
               unname(fm[c("rmsea.scaled", "rmsea.ci.lower.scaled",
                           "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled",
                           "rmsea.notclose.pvalue.scaled")]),
               tolerance = 1e-6)
  # Its robust version reports gamma hat (Savalei, 2018), lavaan's
  # gfi_lrt.robust.
  r <- sem_fit_measures(f, conf_level = 0.90)
  expect_equal(sem_value(r, c("gamma_hat", "gamma_hat_ci_lower", "gamma_hat_ci_upper")),
               unname(lavaan::fitMeasures(f, c("gfi_lrt.robust", "gfi_lrt.ci.lower.robust",
                                               "gfi_lrt.ci.upper.robust"))),
               tolerance = 1e-7)
})

test_that("a mean and variance adjusted statistic reports the robust measures, which lavaan does not compute", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939
  f <- lavaan::cfa(hs_model, data = hs, std.lv = TRUE, estimator = "MLMVS")
  expect_false("rmsea.robust" %in% names(lavaan::fitMeasures(f)))
  # The robust measures subtract c df from the ML statistic, c = tr(U
  # Gamma)/df, whichever corrected statistic accompanies them (Savalei,
  # 2018); the Satterthwaite statistic's s df_s is tr(U Gamma), so they are
  # lavaan's robust measures for the same model fitted by MLM, whose
  # information is the same. Through the 1.1.0 review "auto" reported the
  # scaled measures here, with a warning.
  expect_no_warning(r <- sem_fit_measures(f, conf_level = 0.90))
  expect_equal(attr(r, "sem_fit")$version, "robust")
  expect_match(attr(r, "sem_fit")$correction, "computed by DMAR")
  g <- lavaan::cfa(hs_model, data = hs, std.lv = TRUE, estimator = "MLM")
  y <- sem_fit_measures(g, conf_level = 0.90)
  rows <- c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper", "p_close_fit",
            "p_not_close_fit", "cfi", "tli", "rni", gamma_terms)
  expect_equal(sem_value(r, rows), sem_value(y, rows), tolerance = 1e-7)
  expect_equal(sem_value(r, c("rmsea", "cfi")),
               unname(lavaan::fitMeasures(g, c("rmsea.robust", "cfi.robust"))),
               tolerance = 1e-10)
  # The test of exact fit is the fit's own mean and variance adjusted one.
  expect_equal(sem_value(r, c("chi_square", "df")),
               unname(lavaan::fitMeasures(f, c("chisq.scaled", "df.scaled"))))
  # The scaled measures remain on request: the scaled statistic, with its
  # fractional degrees of freedom, inverted as given.
  x <- sem_fit_measures(f, conf_level = 0.90, version = "scaled")
  expect_equal(attr(x, "sem_fit")$version, "scaled")
  df_s <- sem_value(x, "df")
  expect_false(df_s == round(df_s))
  T_s <- sem_value(x, "chi_square")
  lo <- ncp_at(T_s, df_s, 0.95); hi <- ncp_at(T_s, df_s, 0.05)
  expect_equal(sem_value(x, c("rmsea_ci_lower", "rmsea_ci_upper")),
               sqrt(c(lo, hi) / (301 * df_s)), tolerance = 1e-7)
  expect_equal(sem_value(x, "p_close_fit"),
               stats::pchisq(T_s, df_s, ncp = 301 * df_s * 0.05^2, lower.tail = FALSE))
  # Were the robust measures of a continuous ML fit not to be had, "auto"
  # would report its standard measures, which estimate the population
  # values of the ML measures, and not its scaled ones, which do not
  # (Savalei, 2018), as fa_confirmatory() and measurement_invariance() do.
  expect_warning(s <- with_mocked_bindings(
    sem_fit_measures(f), .sem_robust_mva = function(...) numeric(0)),
    "so the standard versions are reported")
  expect_equal(attr(s, "sem_fit")$version, "standard")
  expect_equal(sem_value(s, c("chi_square", "cfi")),
               unname(lavaan::fitMeasures(f, c("chisq", "cfi"))))
})

test_that("gamma hat is NA for a model that restricts the means", {
  skip_if_not_installed("lavaan")
  set.seed(1)
  n <- 400
  i <- stats::rnorm(n); sl <- stats::rnorm(n, 0.5, 0.4)
  Y <- sapply(0:3, function(t) 2 + i + sl * t + 0.06 * t^2 + stats::rnorm(n, 0, 0.7))
  colnames(Y) <- paste0("y", 1:4)
  Y <- as.data.frame(Y)
  growth <- lavaan::growth("i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
                            s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4", data = Y)
  x <- sem_fit_measures(growth)
  expect_true(all(is.na(sem_value(x, gamma_terms))))
  expect_match(attr(x, "sem_fit")$gamma, "restricts the means")
  expect_false(is.na(sem_value(x, "rmsea")))
  # Intercepts constrained equal restrict the means too.
  eq <- lavaan::cfa("f =~ y1 + y2 + y3 + y4\n y1 ~ a*1\n y2 ~ a*1", data = Y,
                    meanstructure = TRUE)
  expect_true(is.na(sem_value(sem_fit_measures(eq), "gamma_hat")))
  # So do the restrictions a count of free mean parameters missed (the
  # 1.1.0 review): one shared label under ceq.simple = TRUE, which lavaan
  # keeps as one parameter with no == row; an intercept constrained to a
  # constant, or entered in a linear constraint; and an intercept tied to
  # a loading. The means of each fit miss the sample means, and each had a
  # gamma hat recorded as computed (0.873 under ceq.simple = TRUE, 0.759
  # with the loading).
  restricted <- list(
    ceq_simple = lavaan::cfa("f =~ y1 + y2 + y3 + y4\n y1 ~ a*1\n y2 ~ a*1",
                             data = Y, meanstructure = TRUE, ceq.simple = TRUE),
    constant = lavaan::cfa("f =~ y1 + y2 + y3 + y4\n y1 ~ a*1\n a == 4.9",
                           data = Y, meanstructure = TRUE),
    linear = lavaan::cfa(paste("f =~ y1 + y2 + y3 + y4\n y1 ~ a*1\n y2 ~ b*1",
                               "y3 ~ c*1\n a == b + c - 2", sep = "\n"),
                         data = Y, meanstructure = TRUE),
    loading = lavaan::cfa("f =~ l1*y1 + y2 + y3 + y4\n y1 ~ a*1\n a == l1",
                          data = Y, meanstructure = TRUE))
  for (nm in names(restricted)) {
    fr <- restricted[[nm]]
    expect_gt(max(abs(lavaan::lavInspect(fr, "residuals")$mean)), 0.1)
    xr <- sem_fit_measures(fr)
    expect_true(all(is.na(sem_value(xr, gamma_terms))), info = nm)
    expect_match(attr(xr, "sem_fit")$gamma, "restricts the means", info = nm)
  }
  # Free (saturated) means change neither index.
  m1 <- sem_fit_measures(lavaan::cfa("f =~ y1 + y2 + y3 + y4", data = Y,
                                     meanstructure = TRUE))
  m0 <- sem_fit_measures(lavaan::cfa("f =~ y1 + y2 + y3 + y4", data = Y))
  expect_equal(sem_value(m1, gamma_terms), sem_value(m0, gamma_terms))
  # Nor do free means parameterized otherwise (one intercept fixed at zero
  # with the latent mean free), or an equality between two loadings, which
  # restricts the covariances only.
  ml <- sem_fit_measures(lavaan::cfa("f =~ y1 + y2 + y3 + y4\n y1 ~ 0*1\n f ~ 1",
                                     data = Y, meanstructure = TRUE))
  expect_identical(attr(ml, "sem_fit")$gamma, "computed")
  expect_equal(sem_value(ml, gamma_terms), sem_value(m0, gamma_terms),
               tolerance = 1e-6)
  lc <- sem_fit_measures(lavaan::cfa("f =~ a*y1 + a*y2 + y3 + y4", data = Y,
                                     meanstructure = TRUE))
  expect_identical(attr(lc, "sem_fit")$gamma, "computed")
  expect_equal(sem_value(lc, gamma_terms),
               sem_value(sem_fit_measures(lavaan::cfa("f =~ a*y1 + a*y2 + y3 + y4",
                                                      data = Y)), gamma_terms),
               tolerance = 1e-8)
  # Exogenous covariates under fixed.x keep their means out of both counts.
  xs <- sem_fit_measures(lavaan::sem("y4 ~ y1 + y2\n y3 ~ y4", data = Y,
                                     meanstructure = TRUE))
  expect_false(is.na(sem_value(xs, "gamma_hat")))
  # An inequality constraint restricts the means when it binds at the
  # estimates, though lavaan's degrees of freedom do not count it: held
  # above 5, y1's intercept sits at the bound and its fitted mean misses the
  # sample mean, and the first version of the 1.1.0 fix, which set every
  # inequality aside, recorded gamma hat as computed. An inequality that
  # does not bind restricts nothing.
  bind <- lavaan::cfa("f =~ y1 + y2 + y3 + y4\n y1 ~ a*1\n a > 5", data = Y,
                      meanstructure = TRUE)
  expect_equal(as.numeric(lavaan::fitMeasures(bind, "df")), 2)
  expect_gt(max(abs(lavaan::lavInspect(bind, "residuals")$mean)), 1)
  xb <- sem_fit_measures(bind)
  expect_true(all(is.na(sem_value(xb, gamma_terms))))
  expect_match(attr(xb, "sem_fit")$gamma, "restricts the means")
  slack <- lavaan::cfa("f =~ y1 + y2 + y3 + y4\n y1 ~ a*1\n a > 0", data = Y,
                       meanstructure = TRUE)
  expect_identical(attr(sem_fit_measures(slack), "sem_fit")$gamma, "computed")
  expect_equal(sem_value(sem_fit_measures(slack), gamma_terms),
               sem_value(m0, gamma_terms), tolerance = 1e-6)
  # Summary input whose df exceeds the p(p + 1)/2 covariances.
  expect_warning(z <- sem_fit_measures(chisq = 30, df = 12, N = 400, p = 4),
                 "restricts the means")
  expect_true(all(is.na(sem_value(z, gamma_terms))))
  expect_false(is.na(sem_value(z, "rmsea")))
})

test_that("a weighted fit without a corrected test says so and names the remedy", {
  skip_if_not_installed("lavaan")
  w <- lavaan::HolzingerSwineford1939
  set.seed(1)
  w$wt <- stats::runif(nrow(w), 0.5, 2)
  fw <- lavaan::cfa(hs_model, data = w, std.lv = TRUE, sampling.weights = "wt",
                    estimator = "ML", test = "standard", se = "standard")
  expect_warning(x <- sem_fit_measures(fw), "refit with a corrected test")
  expect_match(attr(x, "sem_fit")$version_reason, "sampling weighted fit")
  expect_true(all(is.na(sem_value(x, gamma_terms))))
  expect_error(sem_fit_measures(fw, version = "robust"),
               "not chi square distributed. Refit it with a corrected test")
  # With lavaan's default test the robust measures are there.
  fd <- lavaan::cfa(hs_model, data = w, std.lv = TRUE, sampling.weights = "wt")
  expect_warning(sem_fit_measures(fd, version = "standard"),
                 "the default version = \"auto\" reports the robust measures")
})

test_that("a statistic in the millions gets limits that bracket the estimate", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # lavaan's limits come from R's pchisq(), which does not converge at a
  # noncentrality of 2.8 million: they collapsed onto one value.
  S <- stats::cov(lavaan::HolzingerSwineford1939[paste0("x", 1:9)])
  f <- lavaan::cfa(hs_model, sample.cov = S, sample.nobs = 1e7, std.lv = TRUE)
  expect_no_warning(x <- sem_fit_measures(f, conf_level = 0.90))
  y <- sem_fit_measures(chisq = sem_value(x, "chi_square"), df = 24, N = 1e7,
                        p = 9, conf_level = 0.90, sample_size_convention = "N")
  keep <- c(rmsea_terms, gamma_terms)
  expect_equal(sem_value(x, keep), sem_value(y, keep))
  expect_lt(sem_value(x, "rmsea_ci_lower"), sem_value(x, "rmsea"))
  expect_gt(sem_value(x, "rmsea_ci_upper"), sem_value(x, "rmsea"))
  # The normal approximation to the noncentral chi square, accurate at this
  # size: the half width is z sqrt(2 (df + 2 lambda)) on the noncentrality.
  T <- sem_value(x, "chi_square")
  h <- stats::qnorm(0.95) * sqrt(2 * (24 + 2 * (T - 24)))
  expect_equal(sem_value(x, c("rmsea_ci_lower", "rmsea_ci_upper")),
               sqrt((T - 24 + c(-h, h)) / (1e7 * 24)), tolerance = 1e-5)
})

test_that("the robust limits are recomputed from the scaled statistic when lavaan's cannot be used", {
  skip_if_not_installed("lavaan")
  # A not close fit bound of 4 puts its noncentrality above 1e5, where the
  # robust rows of a complete data robust ML fit are recomputed by
  # Brosseau-Liard et al.'s (2012) Eq. 10. They agree with lavaan's.
  f <- lavaan::cfa(hs_model, data = lavaan::HolzingerSwineford1939, std.lv = TRUE,
                   estimator = "MLR")
  x <- sem_fit_measures(f, conf_level = 0.90, not_close_fit = 4)
  fm <- lavaan::fitMeasures(f, fm.args = list(rmsea.ci.level = 0.90,
                                              rmsea.h0.notclosefit = 4))
  expect_equal(sem_value(x, rmsea_terms),
               unname(fm[c("rmsea.robust", "rmsea.ci.lower.robust",
                           "rmsea.ci.upper.robust", "rmsea.pvalue.robust",
                           "rmsea.notclose.pvalue.robust")]),
               tolerance = 1e-7)
  # A FIML fit's robust statistic cannot be rebuilt: the test whose
  # noncentrality exceeds 1e6 is NA, with a warning, and the rest stand.
  d <- lavaan::HolzingerSwineford1939
  set.seed(113)
  d$x2[sample(nrow(d), 40)] <- NA
  ff <- lavaan::cfa(hs_model, data = d, std.lv = TRUE, missing = "ml")
  expect_warning(z <- sem_fit_measures(ff, not_close_fit = 12),
                 "robust test of not close fit")
  expect_true(is.na(sem_value(z, "p_not_close_fit")))
  expect_false(anyNA(sem_value(z, c("rmsea_ci_lower", "rmsea_ci_upper",
                                    "p_close_fit"))))
})

test_that("a saturated fit has a CFI of 1 in every version", {
  skip_if_not_installed("lavaan")
  # The statistic of a saturated model is zero in every version, so the
  # model term of its incremental indices is zero whatever the scaling
  # constant, which is 0/0 there; lavaan's robust CFI and RNI of this one
  # are NA by rounding.
  f <- lavaan::cfa("f =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
                   data = holzinger_swineford, estimator = "MLR")
  expect_true(is.na(lavaan::fitMeasures(f, "cfi.robust")))
  expect_warning(r <- sem_fit_measures(f), "zero degrees of freedom")
  expect_equal(attr(r, "sem_fit")$version, "robust")
  expect_identical(sem_value(r, "cfi"), 1)
  expect_equal(sem_value(r, "rni"), 1, tolerance = 1e-10)
  expect_true(all(is.na(sem_value(r, c("nfi", "ifi")))))
  expect_warning(s <- sem_fit_measures(f, version = "scaled"),
                 "zero degrees of freedom")
  expect_equal(sem_value(s, c("cfi", "rni", "nfi", "ifi")), rep(1, 4),
               tolerance = 1e-10)
  # With the mean and variance adjusted statistic, whose robust measures
  # DMAR computes, a saturated fit is reported in the robust version too,
  # with no warning that robust measures are missing (the first version of
  # the 1.1.0 fix reported the standard version and asked for a refit).
  fv <- lavaan::cfa("f =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
                    data = holzinger_swineford, estimator = "MLMVS")
  seen <- new.env(parent = emptyenv())
  seen$w <- character(0)
  v <- withCallingHandlers(sem_fit_measures(fv), warning = function(cond) {
    seen$w <- c(seen$w, conditionMessage(cond))
    invokeRestart("muffleWarning")
  })
  expect_length(seen$w, 1L)
  expect_match(seen$w, "zero degrees of freedom")
  expect_equal(attr(v, "sem_fit")$version, "robust")
  expect_identical(sem_value(v, "cfi"), 1)
  expect_equal(sem_value(v, c("tli", "rmsea")), c(1, 0))
  # Nothing was computed from its scaling constant, which is 0/0, so the
  # footer names the corrections as for the same model fitted by MLM; the
  # first version of the rule said DMAR had computed them from it.
  fm <- lavaan::cfa("f =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
                    data = holzinger_swineford, estimator = "MLM")
  m <- suppressWarnings(sem_fit_measures(fm))
  expect_false(grepl("computed by DMAR", attr(v, "sem_fit")$correction,
                     fixed = TRUE))
  expect_false(any(grepl("computed by DMAR", attr(v, "sem_fit")$footer,
                         fixed = TRUE)))
  expect_identical(attr(v, "sem_fit")$correction,
                   attr(m, "sem_fit")$correction)
})

test_that("a model with no degrees of freedom that does not fit exactly keeps its CFI", {
  skip_if_not_installed("lavaan")
  # An inequality bound that binds fixes its parameter as an equality
  # would, but lavaan does not count it among the degrees of freedom: one
  # factor of three variables correlated .8, .8, and .5 has a negative
  # error variance, which bounds = "pos.var" holds at zero, so the
  # statistic is 32.8 on 0 df. Its CFI is .925 (lavaan's, and fa_confirmatory()'s
  # through 1.0.1), not the 1 of a model that fits exactly, which the first
  # version of the rule for saturated fits reported; lavaan's NFI of 1 at
  # zero degrees of freedom is replaced by 1 - T/T_b.
  S <- matrix(c(1, .8, .8, .8, 1, .5, .8, .5, 1), 3, 3,
              dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  f <- lavaan::cfa("f =~ a + b + c", sample.cov = S, sample.nobs = 200,
                   std.lv = TRUE, bounds = "pos.var")
  fm <- lavaan::fitMeasures(f, c("chisq", "df", "baseline.chisq", "cfi",
                                 "rni", "ifi", "nfi"))
  expect_identical(unname(fm[["df"]]), 0)
  expect_gt(fm[["chisq"]], 30)
  expect_identical(unname(fm[["nfi"]]), 1)
  expect_false(DMAR:::.sem_exact_fit(fm))
  expect_true(DMAR:::.sem_exact_fit(c(chisq = 1e-12, df = 0)))
  expect_true(DMAR:::.sem_exact_fit(c(df = 0)))
  expect_false(DMAR:::.sem_exact_fit(c(chisq = 0, df = 1)))
  expect_warning(r <- sem_fit_measures(f), "does not fit exactly",
                 fixed = TRUE)
  expect_equal(sem_value(r, "cfi"), 0.925, tolerance = 1e-3)
  expect_equal(sem_value(r, c("cfi", "rni", "ifi")),
               unname(fm[c("cfi", "rni", "ifi")]), tolerance = 1e-12)
  expect_equal(sem_value(r, "nfi"),
               unname(1 - fm[["chisq"]] / fm[["baseline.chisq"]]),
               tolerance = 1e-12)
  expect_true(all(is.na(sem_value(r, c("tli", "rmsea", "p_chi_square")))))
  expect_match(attr(r, "sem_fit")$boundary, "does not fit exactly",
               fixed = TRUE)
  # The summary statistics give the same CFI and say the same.
  expect_warning(
    s <- sem_fit_measures(chisq = unname(fm[["chisq"]]), df = 0, N = 200,
                          baseline_chisq = unname(fm[["baseline.chisq"]]),
                          baseline_df = 3),
    "does not fit exactly", fixed = TRUE)
  expect_equal(sem_value(s, "cfi"), sem_value(r, "cfi"), tolerance = 1e-12)
})

test_that("a ULS fit of continuous data reports no measure that depends on the units", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # eight ULS fits
  # The ULS statistic, (N - G) times the ULS fit function, is in the units of the
  # squared covariances. lavaan's robust RMSEA of a ULSM fit puts it into
  # the correction of Brosseau-Liard et al. (2012) for the ML fit function:
  # dividing every score by 10 divides it by 100 (16.1 against 0.161 here).
  # The standard RMSEA, CFI, TLI, RNI, and IFI set the degrees of freedom
  # against it and move too; the robust CFI, the scaled measures, and the
  # NFI do not, while the means are free.
  m <- "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
        textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
        speed   =~ t10_addition + t12_counting_groups_of_dots +
                   t13_straight_and_curved_capitals"
  hs <- holzinger_swineford
  h10 <- hs
  items <- lavaan::lavNames(lavaan::lavaanify(m), "ov")
  h10[items] <- h10[items] / 10
  f <- lavaan::cfa(m, data = hs, std.lv = TRUE, estimator = "ULSM")
  f10 <- lavaan::cfa(m, data = h10, std.lv = TRUE, estimator = "ULSM")
  lr <- c(lavaan::fitMeasures(f, "rmsea.robust"),
          lavaan::fitMeasures(f10, "rmsea.robust"))
  expect_equal(unname(lr[1] / lr[2]), 100, tolerance = 1e-6)
  expect_gt(lr[[1]], 16)
  x <- sem_fit_measures(f)
  x10 <- sem_fit_measures(f10)
  expect_equal(attr(x, "sem_fit")$version, "scaled")
  rmsea_rows <- c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper", "p_close_fit",
                  "p_not_close_fit")
  expect_error(sem_fit_measures(f, version = "robust"),
               "ML-derived robust.*not supported.*continuous ULS")
  expect_match(attr(x, "sem_fit")$uls, "conventional ML cutoffs")
  expect_output(print(x), "Scaled descriptive indices", fixed = TRUE)
  s <- sem_fit_measures(f, version = "scaled")
  s10 <- sem_fit_measures(f10, version = "scaled")
  expect_identical(x$value, s$value)
  expect_identical(x10$value, s10$value)
  expect_equal(sem_value(s, "rmsea"), 0.1595, tolerance = 1e-3)
  keep <- c("rmsea", "rmsea_ci_lower", "cfi", "tli", "rni", "nfi", "ifi")
  expect_equal(sem_value(s, keep), sem_value(s10, keep), tolerance = 1e-6)
  # The standard version: only the NFI, a ratio of two ULS statistics,
  # keeps its value (the CFI was .9767 and .9793, the first version of the
  # 1.1.0 fix reported it).
  expect_warning(st <- sem_fit_measures(f, version = "standard"),
                 "(N - G) times the ULS fit function", fixed = TRUE)
  expect_warning(st10 <- sem_fit_measures(f10, version = "standard"),
                 "eligible scaled")
  std_rows <- c(rmsea_rows, "cfi", "tli", "rni", "ifi")
  expect_true(all(is.na(sem_value(st, std_rows))))
  expect_equal(sem_value(st, "nfi"), unname(lavaan::fitMeasures(f, "nfi")))
  expect_equal(sem_value(st, "nfi"), sem_value(st10, "nfi"), tolerance = 1e-6)
  expect_false(isTRUE(all.equal(unname(lavaan::fitMeasures(f, "cfi")),
                                unname(lavaan::fitMeasures(f10, "cfi")),
                                tolerance = 1e-4)))
  # estimator = "ULS" with test = "standard" reports the same version under
  # "auto", and its remedy is not a scaled version it does not have.
  u0 <- lavaan::cfa(m, data = hs, std.lv = TRUE, estimator = "ULS",
                    test = "standard")
  expect_warning(xs <- sem_fit_measures(u0), "Leave out test = \"standard\"")
  wu <- tryCatch(sem_fit_measures(u0), warning = conditionMessage)
  expect_false(grepl("version = \"scaled\"", wu, fixed = TRUE))
  expect_true(all(is.na(sem_value(xs, std_rows))))
  # ULS alone is tested by Browne's residual based statistic, which has no
  # units, and keeps every measure.
  u <- lavaan::cfa(m, data = hs, std.lv = TRUE, estimator = "ULS")
  expect_no_warning(xu <- sem_fit_measures(u))
  expect_equal(sem_value(xu, c("rmsea", "cfi")),
               unname(lavaan::fitMeasures(u, c("rmsea", "cfi"))))
  expect_false(anyNA(sem_value(xu, c("rmsea", "cfi", "nfi"))))
  # A model that restricts the means adds the squared mean residuals to the
  # squared covariance residuals, in other units: the robust CFI was .979
  # for the scores and .972 for the scores divided by 10, and the scaled
  # RMSEA .160 and .181. No version reports an index; the chi_square row
  # and the SRMR are kept.
  mr <- paste(m, "t6_paragraph_comprehension ~ a*1",
              "t9_word_meaning ~ a*1", sep = "\n")
  r <- lavaan::cfa(mr, data = hs, std.lv = TRUE, estimator = "ULSM",
                   meanstructure = TRUE)
  r10 <- lavaan::cfa(mr, data = h10, std.lv = TRUE, estimator = "ULSM",
                     meanstructure = TRUE)
  expect_equal(unname(lavaan::fitMeasures(r, "cfi.robust")), 0.9787, tolerance = 1e-3)
  expect_equal(unname(lavaan::fitMeasures(r10, "cfi.robust")), 0.9725, tolerance = 1e-3)
  all_rows <- c(rmsea_rows, "cfi", "tli", "rni", "nfi", "ifi")
  expect_error(sem_fit_measures(r, version = "robust"), "not supported")
  for (v in c("scaled", "standard")) {
    expect_warning(xr <- sem_fit_measures(r, version = v), "restricts the means",
                   info = v)
    expect_true(all(is.na(sem_value(xr, all_rows))), info = v)
    expect_false(anyNA(sem_value(xr, c("chi_square", "srmr"))), info = v)
    expect_match(attr(xr, "sem_fit")$uls, "restricts the means", fixed = TRUE,
                 info = v)
  }
  # Browne's statistic of the same model has no units, and the equal
  # intercepts restrict the means through a parameter the covariances do
  # not share (the latent means are fixed), so the ULS estimates do not
  # change with the units either and the measures are kept.
  rb <- lavaan::cfa(mr, data = hs, std.lv = TRUE, estimator = "ULS",
                    meanstructure = TRUE)
  rb10 <- lavaan::cfa(mr, data = h10, std.lv = TRUE, estimator = "ULS",
                      meanstructure = TRUE)
  expect_true(DMAR:::.sem_means_restricted(rb))
  expect_false(DMAR:::.sem_means_linked(rb))
  expect_no_warning(xb <- sem_fit_measures(rb))
  expect_equal(sem_value(xb, c("rmsea", "cfi", "tli", "nfi")),
               sem_value(sem_fit_measures(rb10), c("rmsea", "cfi", "tli", "nfi")),
               tolerance = 1e-6)
  expect_true(is.na(attr(xb, "sem_fit")$uls))
})

test_that("Browne's measures of a ULS fit whose estimates change with the units are NA", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # six ULS fits
  # Browne's residual based statistic has no units, but it is evaluated at
  # the ULS estimates, which change with the units when the model restricts
  # the means through parameters the covariances share: the ULS fit
  # function is k^4 times its covariance part plus k^2 times its mean part.
  # Through the first version of the ULS rule the measures were reported:
  # on the scalar rung of the first example of measurement_invariance(), a
  # CFI of .8874 for the scores and .8861 for the scores divided by 10, and
  # an RMSEA of .0680 and .0684 (marker identification, as there).
  items <- c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning",
             "t20_deduction", "t22_problem_reasoning", "t23_series_completion")
  m <- paste("verbal =~", paste(items[1:3], collapse = " + "),
             "\ndeduction =~", paste(items[4:6], collapse = " + "))
  hs <- holzinger_swineford
  h10 <- hs
  h10[items] <- h10[items] / 10
  scalar <- function(d) {
    lavaan::cfa(m, data = d, group = "school", estimator = "ULS",
                group.equal = c("loadings", "intercepts"))
  }
  s <- scalar(hs)
  s10 <- scalar(h10)
  expect_true(any(startsWith(lavaan::lavInspect(s, "options")$test,
                             "browne.residual")))
  expect_equal(unname(lavaan::fitMeasures(s, "cfi")), 0.8874, tolerance = 1e-3)
  expect_equal(unname(lavaan::fitMeasures(s10, "cfi")), 0.8861, tolerance = 1e-3)
  expect_true(DMAR:::.sem_means_linked(s))
  expect_true(DMAR:::.sem_uls_estimates_units(s))
  rmsea_rows <- c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper", "p_close_fit",
                  "p_not_close_fit")
  all_rows <- c(rmsea_rows, "cfi", "tli", "rni", "nfi", "ifi")
  expect_warning(x <- sem_fit_measures(s), "and this model does",
                 fixed = TRUE)
  expect_equal(attr(x, "sem_fit")$version, "standard")
  expect_true(all(is.na(sem_value(x, all_rows))))
  expect_false(anyNA(sem_value(x, c("chi_square", "p_chi_square", "srmr"))))
  expect_equal(sem_value(x, "chi_square"),
               unname(lavaan::fitMeasures(s, "chisq")))
  expect_match(attr(x, "sem_fit")$uls, "Browne's residual based statistic",
               fixed = TRUE)
  expect_output(print(x), "ULS of continuous data: Browne's", fixed = TRUE)
  # The same model with the means free (loadings equal only) keeps them,
  # and they do not change with the units.
  metric <- function(d) {
    lavaan::cfa(m, data = d, group = "school", estimator = "ULS",
                group.equal = "loadings")
  }
  expect_no_warning(xm <- sem_fit_measures(metric(hs)))
  xm10 <- sem_fit_measures(metric(h10))
  expect_false(anyNA(sem_value(xm, all_rows)))
  expect_equal(sem_value(xm, all_rows), sem_value(xm10, all_rows),
               tolerance = 1e-5)
  # In one group: a latent mean freed beside two intercepts fixed at zero
  # puts the loadings into the means.
  m1 <- paste("verbal =~", paste(items[1:3], collapse = " + "),
              "\nt6_paragraph_comprehension ~ 0*1\nt7_sentence ~ 0*1",
              "\nverbal ~ NA*1")
  one <- function(d) {
    lavaan::cfa(m1, data = d, std.lv = TRUE, estimator = "ULS",
                meanstructure = TRUE)
  }
  o <- one(hs)
  o10 <- one(h10)
  expect_false(isTRUE(all.equal(unname(lavaan::fitMeasures(o, "cfi")),
                                unname(lavaan::fitMeasures(o10, "cfi")),
                                tolerance = 1e-5)))
  expect_warning(xo <- sem_fit_measures(o), "Browne's residual based",
                 fixed = TRUE)
  expect_true(all(is.na(sem_value(xo, all_rows))))
})


test_that("H1 screening respects fixed-predictor likelihood scope", {
  skip_if_not_installed("lavaan")
  set.seed(31)
  d <- data.frame(x1 = rnorm(100), x2 = rnorm(100),
                  y1 = rnorm(100), y2 = rnorm(100))
  d$x2 <- d$x1
  model <- paste("y1 ~ x1 + x2", "y2 ~ x1 + x2", "y1 ~~ 1*y1",
                  "y2 ~~ 1*y2", "y1 ~~ 0*y2", sep = "\n")
  fit <- suppressWarnings(lavaan::sem(model, data = d, missing = "ml",
                                        fixed.x = TRUE))
  # Redundant X can give nonunique slopes, but finite conditional
  # likelihood. The joint-normal collapse path would change fixed X.
  expect_true(lavaan::lavInspect(fit, "converged"))
  expect_true(is.finite(as.numeric(lavaan::logLik(fit))))
  expect_null(.fiml_fit_h1_diagnostic(fit))
  expect_false(is.null(.fiml_likelihood_obstacle(d)))
  d$y2 <- 2 * d$y1 + 1
  response_fit <- suppressWarnings(lavaan::sem(model, data = d,
    missing = "ml", fixed.x = TRUE))
  why <- .fiml_fit_h1_diagnostic(response_fit)
  expect_match(why, "y1 and y2", fixed = TRUE)
  expect_match(why, "response variables only", fixed = TRUE)
  expect_match(why, "not assessed", fixed = TRUE)
  # An exact Y-on-X relation illustrates an intentionally missed
  # obstacle, not a claim that a conditional maximum exists.
  d$y1 <- d$x1
  d$y2 <- rnorm(100)
  exact_fit <- suppressWarnings(lavaan::sem(model, data = d,
    missing = "ml", fixed.x = TRUE))
  expect_null(.fiml_fit_h1_diagnostic(exact_fit))
})


test_that("H1-dependent bootstrap points screen the selected resample", {
  set.seed(73)
  d <- as.data.frame(matrix(rnorm(120), 30, 4))
  names(d) <- paste0("y", 1:4)
  # Four distinct complete rows are on one hyperplane; rows missing y4
  # add distinct patterns but cannot prevent that full-support collapse.
  d$y4[1:4] <- rowSums(d[1:4, 1:3])
  d$y4[5:15] <- NA_real_
  expect_null(.fiml_likelihood_obstacle(d))
  selected <- d[1:15, ]
  expect_gt(nrow(unique(selected)), ncol(selected))
  expect_false(is.null(.fiml_likelihood_obstacle(selected)))
  for (missing in c("fiml", "two_stage")) {
    called <- FALSE
    checked <- .relia_checked_point(function(x) {
      called <<- TRUE
      1
    }, missing = missing, items = names(d), fitted = TRUE,
    needs_saturated = TRUE)
    expect_identical(checked(d), 1)
    called <- FALSE
    expect_error(checked(selected), class = "dmar_fiml_unbounded")
    expect_false(called)
    # H0-only points retain their separate scope.
    h0 <- .relia_checked_point(function(x) 1, missing = missing,
                               items = names(d), fitted = TRUE)
    expect_identical(h0(selected), 1)
  }
})

test_that("with auxiliary variables the SRMR is over the model's own variables", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939
  hs$agemo2 <- hs$agemo + 12 * hs$ageyr
  miss <- hs
  set.seed(1)
  for (j in paste0("x", 1:9)) miss[sample(nrow(miss), 40), j] <- NA
  # The SRMR of the nine items, computed apart from DMAR from lavaan's h1
  # and implied moments with the Bentler standardization.
  own <- function(fit) {
    S <- lavaan::lavInspect(fit, "sampstat.h1")$cov
    E <- lavaan::lavInspect(fit, "implied")$cov
    k <- paste0("x", 1:9)
    r <- (S[k, k] - E[k, k]) / tcrossprod(sqrt(diag(S)[k]))
    sqrt(mean(r[lower.tri(r, diag = TRUE)]^2))
  }
  srmr <- function(fit) {
    d <- sem_fit_measures(fit)
    d$value[d$term == "srmr"]
  }
  fiml <- function(dat, aux = NULL) {
    suppressMessages(lavaan::cfa(hs_model, data = dat, std.lv = TRUE,
                                 missing = "fiml", aux = aux))
  }
  for (dat in list(hs, miss)) {
    f0 <- fiml(dat)
    f2 <- fiml(dat, c("ageyr", "agemo2"))
    # Without auxiliary variables the row is lavaan's srmr_bentler_nomean.
    expect_equal(srmr(f0), lavaan::fitMeasures(f0, "srmr_bentler_nomean")[[1L]],
                 tolerance = 1e-12)
    expect_equal(srmr(f2), own(f2), tolerance = 1e-10)
    # lavaan's own measure averages the auxiliary variables in and is lower.
    expect_lt(lavaan::fitMeasures(f2, "srmr_bentler_nomean")[[1L]], own(f2) - 0.005)
    expect_match(attr(sem_fit_measures(f2), "sem_fit")$srmr,
                 "auxiliary variables set aside: ageyr, agemo2", fixed = TRUE)
  }
  # With complete data an auxiliary variable changes nothing in the model.
  expect_equal(srmr(fiml(hs, "ageyr")), srmr(fiml(hs)), tolerance = 1e-4)
  # Two groups are combined as lavaan combines them, by their numbers of cases.
  g0 <- suppressMessages(lavaan::cfa(hs_model, data = hs, group = "school",
                                     missing = "fiml"))
  g1 <- suppressMessages(lavaan::cfa(hs_model, data = hs, group = "school",
                                     missing = "fiml", aux = "ageyr"))
  expect_equal(.sem_srmr_own(g0, character(0)),
               lavaan::fitMeasures(g0, "srmr_bentler_nomean")[[1L]], tolerance = 1e-12)
  expect_equal(srmr(g1), srmr(g0), tolerance = 1e-4)
})
