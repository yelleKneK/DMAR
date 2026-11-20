test_that("cv_t() returns the expected tidy data frame for two-sided alpha", {
  result <- cv_t(alpha_level = 0.05, df = 30)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "area_less", "area_greater"))
  expect_equal(result$term, c("lower_cv", "upper_cv"))
})

test_that("cv_t() returns a dmar_tbl", {
  expect_s3_class(cv_t(alpha_level = .05, df = 13), "dmar_tbl")
  expect_s3_class(cv_t(alpha_lower = 0, alpha_upper = .05, df = 13, verbose = FALSE), "dmar_tbl")
})

test_that("cv_t() at alpha = 0.05, df = 30 matches qt(0.025, 30) and qt(0.975, 30)", {
  result <- cv_t(alpha_level = 0.05, df = 30)
  expect_equal(result$value[result$term == "lower_cv"], qt(0.025, df = 30), tolerance = 1e-12)
  expect_equal(result$value[result$term == "upper_cv"], qt(0.975, df = 30), tolerance = 1e-12)
})

test_that("cv_t() with ncp != 0 returns the noncentral t critical values", {
  result <- cv_t(alpha_level = 0.05, df = 30, ncp = 2)
  expect_equal(result$value[result$term == "lower_cv"], qt(0.025, df = 30, ncp = 2), tolerance = 1e-10)
  expect_equal(result$value[result$term == "upper_cv"], qt(0.975, df = 30, ncp = 2), tolerance = 1e-10)
})

test_that("the lower-tail alternative accepts 'less', the counterpart of 'greater'", {
  # The upper-tail accept list carried "greater" while the lower-tail list
  # carried "lesser" but not "less", so the natural counterpart of the
  # spelling that worked did not. It did not error cleanly either: the call
  # fell through both branches and died on an unset condition with
  # "missing value where TRUE/FALSE needed".
  expect_equal(cv_t(alpha_level = .05, df = 13, alternative = "less"),
               cv_t(alpha_level = .05, df = 13, alternative = "lesser"))
  expect_equal(cv_f(alpha_level = .05, df_numerator = 2, df_denominator = 20,
                    alternative = "less"),
               cv_f(alpha_level = .05, df_numerator = 2, df_denominator = 20,
                    alternative = "lesser"))
  expect_equal(cv_chisq(alpha_level = .05, df = 3, alternative = "less"),
               cv_chisq(alpha_level = .05, df = 3, alternative = "lesser"))
  # And "less" really is routed to the lower tail: it spends the whole of
  # alpha there, so its lower critical value sits closer to zero than the
  # two-sided one, and nothing is cut off above.
  one <- cv_t(alpha_level = .05, df = 13, alternative = "less")
  two <- cv_t(alpha_level = .05, df = 13, alternative = "not_equal")
  expect_gt(one$value[one$term == "lower_cv"],
            two$value[two$term == "lower_cv"])
  expect_identical(one$value[one$term == "upper_cv"], Inf)
})

test_that("cv_t() gives exact noncentral t critical values past a noncentrality of 37.62", {
  # R's qt() with ncp is a normal approximation there; DMAR through
  # 3338679 returned it: 429.17833 for the upper value at ncp = 50 and
  # df = 3 (exactly 186.49387) and 83.93688 at ncp = 45 and df = 10
  # (exactly 79.08391).
  r <- expect_no_warning(cv_t(alpha_level = .05, df = 3, ncp = 50))
  expect_equal(r$value, c(28.2828377277, 186.4938651548), tolerance = 1e-10)
  expect_equal(r$area_less, c(.025, .975), tolerance = 1e-10)
  expect_equal(r$area_greater, c(.975, .025), tolerance = 1e-10)
  expect_equal(.oracle_nct(r$value[2], 3, 50, lower_tail = FALSE), .025,
               tolerance = 1e-10)
  r <- cv_t(alpha_level = .05, df = 10, ncp = 45)
  expect_equal(r$value, c(31.3532061443, 79.0839051910), tolerance = 1e-10)
  # One-sided, a negative noncentrality: the lower value is -Inf.
  r <- cv_t(alpha_level = .05, df = 5, ncp = -40, alternative = "less")
  expect_identical(r$value[2], Inf)
  expect_identical(c(r$area_less[2], r$area_greater[2]), c(1, 0))
  expect_equal(.oracle_nct(r$value[1], 5, -40), .05, tolerance = 1e-10)
  # Inside the range of pt(): unchanged, and without the warnings qt() gave.
  r <- expect_no_warning(cv_t(alpha_level = .05, df = 2000, ncp = 37))
  expect_equal(r$value, suppressWarnings(qt(c(.025, .975), 2000, 37)),
               tolerance = 1e-11)
})

test_that("cv_t() with fewer than one degree of freedom is exact, and does not stop with integrate()'s error", {
  # Through the first 1.1.0 builds these stopped with an unclassed error
  # from stats::integrate(), "the integral is probably divergent".
  for (x in list(list(.001, .5, 38), list(.05, .2, 45))) {
    r <- expect_no_error(cv_t(x[[1]], x[[2]], ncp = x[[3]]))
    a <- x[[1]] / 2
    expect_equal(.oracle_nct(r$value[1], x[[2]], x[[3]]), a, tolerance = 1e-9)
    expect_equal(.oracle_nct(r$value[2], x[[2]], x[[3]], lower_tail = FALSE), a,
                 tolerance = 1e-9)
    expect_equal(.oracle_nct_series(r$value[2], x[[2]], x[[3]], lower_tail = FALSE), a,
                 tolerance = 1e-9)
    expect_equal(r$area_less, c(a, 1 - a), tolerance = 1e-9)
  }
  # Inside 37.62 the values were pt()'s, off in the far tail: the lower
  # critical value at df = 0.5, ncp = -37.62, and alpha .001 was -85697561.7
  # where the .0005 quantile is -91565073.6 (by both oracles).
  r <- cv_t(.001, .5, ncp = -37.62)
  expect_equal(r$value[1], -91565073.64, tolerance = 1e-9)
  expect_equal(.oracle_nct_series(r$value[1], .5, -37.62), 5e-4, tolerance = 1e-9)
  # The central t (ncp = 0) is qt() at any df.
  expect_equal(cv_t(.05, .5)$value, stats::qt(c(.025, .975), .5))
})

test_that("cv_t() returns a critical value past the largest double as -Inf or Inf, and never stops with uniroot()'s error", {
  # Through the first 1.1.0 builds these stopped with the unclassed
  # "no sign change found in 1000 iterations" of stats::uniroot(): with so
  # few degrees of freedom the quantile can be past the largest double.
  big <- .Machine$double.xmax
  # df = 0.001, ncp = 5: the .025 quantile is finite, the .975 one infinite.
  r <- cv_t(.05, df = .001, ncp = 5)
  expect_equal(r$value[1], 1.444368535e10, tolerance = 1e-9)
  expect_equal(.oracle_nct_far(r$value[1], .001, 5), 1 - .025, tolerance = 1e-9)
  expect_identical(r$value[2], Inf)
  expect_gte(.oracle_nct_far(big, .001, 5), .025)
  # df = 0.001, ncp = -1: both, as qt() gives them with ncp = 0.
  r <- cv_t(.05, .001, ncp = -1)
  expect_identical(r$value, c(-Inf, Inf))
  expect_gte(.oracle_nct_far(-big, .001, -1), .025)
  expect_gte(.oracle_nct_far(big, .001, -1), .025)
  expect_identical(cv_t(.05, .001)$value, c(-Inf, Inf))
  # alpha 1e-6 at df = 0.01 and 0.001.
  r <- cv_t(1e-6, .01, ncp = -100)
  expect_identical(r$value[1], -Inf)
  expect_gte(.oracle_nct_far(-big, .01, -100), 5e-7)
  expect_equal(.oracle_nct_far(r$value[2], .01, -100), 1 - 5e-7, tolerance = 1e-9)
  # A tail near 1e-6 carries the absolute tolerance of the integral, 1e-13,
  # so it is checked to 1e-7 of itself (it is within 2e-8 here; the series
  # oracle agrees with .oracle_nct() to 1e-15).
  for (df in c(.001, .01)) {
    r <- cv_t(1e-6, df, ncp = 5)
    expect_equal(.oracle_nct(r$value[1], df, 5), 5e-7, tolerance = 1e-7)
    expect_equal(.oracle_nct_series(r$value[1], df, 5), 5e-7, tolerance = 1e-7)
    expect_identical(r$value[2], Inf)
    expect_gte(.oracle_nct_far(big, df, 5), 5e-7)
  }
  for (x in list(list(2e-6, .01, 5), list(.002, .005, 5))) {
    r <- expect_no_error(cv_t(x[[1]], x[[2]], ncp = x[[3]]))
    expect_equal(.oracle_nct(r$value[1], x[[2]], x[[3]]), x[[1]] / 2, tolerance = 1e-7)
    expect_identical(r$value[2], Inf)
    expect_gte(.oracle_nct_far(big, x[[2]], x[[3]]), x[[1]] / 2)
  }
  # Finite quantiles beyond the reach of uniroot()'s extension of its
  # bracket (about 1e301 times its first step), found between it and the
  # largest double: by the far-tail oracle each tail is at its target.
  r <- cv_t(.05, .005, ncp = -1)
  expect_true(is.finite(r$value[1]) && r$value[1] < -1e303)
  expect_equal(.oracle_nct_far(r$value[1], .005, -1), .025, tolerance = 1e-9)
  r <- cv_t(.001, .01, ncp = .3)
  expect_true(is.finite(r$value[2]) && r$value[2] > 1e307)
  expect_equal(.oracle_nct_far(r$value[2], .01, .3), 5e-4, tolerance = 1e-9)
})
