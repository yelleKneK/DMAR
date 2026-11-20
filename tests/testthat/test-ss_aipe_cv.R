test_that("ss_aipe_cv() with assurance returns one row, the published sample size", {
  plain <- ss_aipe_cv(true_cv = .25, width = .10, conf_level = .99)
  res <- ss_aipe_cv(true_cv = .25, width = .10, conf_level = .99, assurance = .99)
  expect_identical(plain$term, "necessary_N")
  expect_identical(res$term, "necessary_N")
  expect_named(res, c("term", "value"))
  expect_s3_class(res, "dmar_tbl")
  expect_s3_class(res, "dmar_ss_aipe")
  # Kelley (2007, pp. 759-760) gives 99 and, with assurance .99, 141, from
  # R's approximation to the noncentral t; the exact plans are 101 and 142.
  expect_equal(plain$value, 101)
  expect_equal(res$value, 142)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_aipe_cv() returns a 1-row tidy data frame with necessary_N", {
  result <- ss_aipe_cv(true_cv = 0.2, width = 0.1, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, "necessary_N")
})

test_that("ss_aipe_cv() requires a larger n for a narrower width", {
  skip_on_cran()
  # width = 0.05 forces n large enough that the noncentrality passes 37.62,
  # where R's pt() is an approximation; the interval is exact there and
  # nothing warns.
  n_wide   <- ss_aipe_cv(true_cv = 0.2, width = 0.20, conf_level = 0.95)$value
  n_narrow <- expect_no_warning(ss_aipe_cv(true_cv = 0.2, width = 0.05, conf_level = 0.95))$value
  expect_lt(n_wide, n_narrow)
})

test_that("ss_aipe_cv() errors on non-positive coefficient of variation", {
  expect_error(ss_aipe_cv(true_cv = 0,  width = 0.1, conf_level = 0.95), "positive")
  expect_error(ss_aipe_cv(true_cv = -1, width = 0.1, conf_level = 0.95), "positive")
})

test_that("ss_aipe_cv() errors on a non-positive width instead of never terminating", {
  # A non-positive width makes the search target unreachable, and the
  # increment-by-one search would not terminate. It must error at entry
  # rather than hang.
  expect_error(ss_aipe_cv(true_cv = 0.1, width = -1, conf_level = 0.99), "width")
  expect_error(ss_aipe_cv(true_cv = 0.1, width = 0, conf_level = 0.99), "width")
})

test_that("ss_aipe_cv() validates conf_level", {
  expect_error(ss_aipe_cv(true_cv = 0.1, width = 0.1, conf_level = 1.2),
               "conf_level")
})

test_that("ss_aipe_cv() accepts the mu/sigma parameterization", {
  skip_on_cran()
  # Regression: neither formal was read, so the documented mu/sigma call
  # died on "argument is of length zero". mu = 10 and sigma = 1 imply
  # true_cv = .1; both parameterizations of the same population must plan
  # the same sample size, anchored at N = 20 for a 99 percent interval of
  # width .10. Table 3 of Kelley (2007) prints 19, computed with R's
  # approximation to the noncentral t past a noncentrality of 37.62
  # (sqrt(19) / .1 = 43.6); with the exact interval 19 is .1003 wide at the
  # expected value.
  n_musigma <- expect_no_warning(
    ss_aipe_cv(true_mean = 10, true_sd = 1, width = .1, conf_level = .99)
  )
  n_cv <- ss_aipe_cv(true_cv = .1, width = .1, conf_level = .99)
  expect_named(n_musigma, c("term", "value"))
  expect_equal(n_musigma$term, "necessary_N")
  expect_equal(n_musigma$value, n_cv$value)
  expect_equal(n_musigma$value, 20)
})

test_that("ss_aipe_cv() rejects conflicting or incomplete parameterizations", {
  expect_error(
    ss_aipe_cv(true_cv = .1, true_mean = 10, true_sd = 1, width = .1),
    "either 'true_cv' or both 'true_mean' and 'true_sd'"
  )
  expect_error(ss_aipe_cv(true_mean = 10, width = .1),
               "both 'true_mean' and 'true_sd'")
  expect_error(ss_aipe_cv(true_sd = 1, width = .1),
               "both 'true_mean' and 'true_sd'")
})

test_that("ss_aipe_cv() plans on the expected value of k, as Kelley (2007) does", {
  skip_on_cran()
  # Kelley (2007, p. 758) substitutes E[k | kappa, N] of Equation 3 for kappa
  # at each N. Worked example, pp. 759-760 and the Appendix: kappa = .25,
  # width = .10, 99 percent: the article gives N = 99, and 141 with
  # assurance .99, computed with R's approximation to the noncentral t
  # (sqrt(99) / .25 = 39.8, past 37.62). With the exact interval the plans
  # are 101 and 142 (the interval at the expected value is .1006 wide at
  # 99 and .0995 at 101). MBESS plans on kappa (1 + 1/(4N)) with the
  # approximation and returns 100 and 142.
  expect_equal(ss_aipe_cv(true_cv = .25, width = .10, conf_level = .99)$value, 101)
  expect_equal(ss_aipe_cv(true_cv = .25, width = .10, conf_level = .99, assurance = .99)$value[1], 142)
  # Table 1 (90 percent), expected-width subtable, width = .10: every
  # printed value.
  kappa <- seq(.05, .50, .05)
  n <- vapply(kappa, function(k) ss_aipe_cv(true_cv = k, width = .10, conf_level = .90)$value, 0)
  expect_equal(n, c(5, 9, 17, 27, 42, 62, 87, 119, 159, 209))
  # The width at the expected value meets the target at N and not at N - 1.
  w <- function(k, N) { ci <- ci_cv(cv = DMAR:::.cv_expected_value(k, N), N = N, conf_level = .90); diff(ci$value[ci$term %in% c("ci_lower", "ci_upper")]) }
  expect_lte(w(.30, 62), .10)
  expect_gt(w(.30, 61), .10)
})

test_that("ss_aipe_cv() against the rest of Table 1 and Table 3 of Kelley (2007)", {
  skip_on_cran()
  # The tables were computed with R's approximation to the noncentral t
  # past a noncentrality of 37.62; the plans here use the exact interval
  # and differ from print by one where a narrow width puts the
  # noncentrality past it. Each printed row is in the comment above its
  # check. Of the 792 printed cells the function returns 468 (the
  # evidence of the change that made the intervals exact,
  # work/round6/evidence/e1_nct_exact/12_kelley2007_summary.log).
  plan <- function(k, w, ...) ss_aipe_cv(true_cv = k, width = w, conf_level = .90, ...)$value[1]
  # Printed: 26, 93, 208, 379, 614.
  expect_equal(vapply(c(.05, .10, .15, .20, .25), plan, 0, w = .025), c(26, 92, 208, 378, 613))
  # Printed: 9, 27, 56, 98, 157, 235, 335, 463, 622, 818.
  expect_equal(vapply(seq(.05, .50, .05), plan, 0, w = .05), c(9, 26, 55, 98, 156, 234, 335, 462, 621, 817))
  # Printed and returned: 5, 7, 10, 14, 19, 25, 33, 44, 56.
  expect_equal(vapply(seq(.10, .50, .05), plan, 0, w = .20), c(5, 7, 10, 14, 19, 25, 33, 44, 56))
  # Table 3 (99 percent), expected-width subtable, width = .10. Printed:
  # 9, 19, 37, 63, 99, 147, 209, 287, 389, 510. At a coefficient of .45
  # or .50 the search starts at N = 4, where the upper limit of ci_cv()'s
  # interval is infinite; an infinite width is not narrow enough, the
  # search goes on, and nothing is signaled.
  plan_99 <- function(k) ss_aipe_cv(true_cv = k, width = .10, conf_level = .99)$value
  expect_equal(expect_no_warning(vapply(seq(.05, .50, .05), plan_99, 0)),
               c(9, 20, 38, 64, 101, 149, 211, 289, 387, 508))
  # Assurance .99 subtable. Printed: 16, 43, 82, 137, 209, 303.
  expect_equal(vapply(seq(.05, .30, .05), plan, 0, w = .05, assurance = .99), c(17, 43, 82, 137, 209, 303))
  # Printed and returned: 8, 18, 31, 49, 71, 100, 136, 180, 235, 301.
  expect_equal(vapply(seq(.05, .50, .05), plan, 0, w = .10, assurance = .99), c(8, 18, 31, 49, 71, 100, 136, 180, 235, 301))
})

test_that("ss_aipe_cv() does not accept the DMAR 1.0.x name C_of_V", {
  # C_of_V is not a formal and the function takes no `...`, so the name is
  # an unused argument.
  expect_false("C_of_V" %in% names(formals(ss_aipe_cv)))
  expect_error(ss_aipe_cv(C_of_V = .25, width = .10, conf_level = .99),
               "unused argument")
  expect_identical(names(formals(ss_aipe_cv)),
                   c("true_cv", "width", "conf_level", "assurance",
                     "true_mean", "true_sd", "alpha_lower", "alpha_upper"))
})

test_that(".cv_expected_value() is Equation 3 of Kelley (2007)", {
  expect_equal(DMAR:::.cv_expected_value(.2, 98), .2 * (1 - 1 / (4 * 97) + .04 / 98 + 1 / (2 * 97^2)))
  # Agreement with the mean of the sampling distribution of k (quantile average).
  p <- (seq_len(20000) - .5) / 20000
  k <- sqrt(200) / suppressWarnings(qt(1 - p, df = 199, ncp = sqrt(200) / .2))
  expect_equal(DMAR:::.cv_expected_value(.2, 200), mean(k), tolerance = 1e-4)
})

# ---- The assurance step: the published procedure of Kelley (2007) ---------

test_that("ss_aipe_cv() runs the published procedure for any assurance strictly between 0 and 1", {
  skip_on_cran()
  plan <- function(a) ss_aipe_cv(true_cv = .22, width = .09, conf_level = .95, assurance = a)$value
  # The cell named on the help page: 57 for a requested .60 (55 without
  # assurance). With R's approximation to the noncentral t the plan for
  # .60 was 56.
  expect_equal(ss_aipe_cv(true_cv = .22, width = .09, conf_level = .95)$value, 55)
  expect_equal(plan(.60), 57)
  # The sample size does not decrease as the request rises, through .50.
  n <- vapply(c(.20, .40, .50, .60, .80, .95, .99), plan, 0)
  expect_false(is.unsorted(n))
  expect_true(all(n == round(n)))
  expect_gt(n[7], n[3])
})

test_that("ss_aipe_cv() refuses a malformed assurance with one message", {
  skip_on_cran()
  for (bad in list(c(.8, .9), "0.8", NA, NA_real_, 1, 0, -1, 1.2, Inf)) {
    err <- tryCatch(
      suppressWarnings(ss_aipe_cv(true_cv = .2, width = .1, assurance = bad)),
      error = function(e) e)
    expect_s3_class(err, "dmar_assurance_invalid")
    expect_s3_class(err, "dmar_refusal")
    expect_match(conditionMessage(err),
                 "'assurance' must be NULL or a single number in (0, 1)",
                 fixed = TRUE)
    expect_identical(err$arg, "assurance")
    expect_identical(err$value, bad)
  }
})

test_that("ss_aipe_cv() has no dots: an unknown or 1.0.x name is an unused argument", {
  expect_error(ss_aipe_cv(true_cv = .2, width = .1, assurance = .9,
                          assurance_method = "exact"),
               "unused argument")
  expect_error(ss_aipe_cv(true_cv = .22, C_of_V = .5, width = .09),
               "unused argument")
  expect_false("..." %in% names(formals(ss_aipe_cv)))
})

test_that("ss_aipe_cv() with assurance .50 is the published procedure, not the plan without assurance", {
  skip_on_cran()
  # The help states this pair: a request of .50 can return one fewer than
  # the plan without assurance. (Another pair, true_cv = .10 and width =
  # .075, gives 18 both ways.)
  with_half <- ss_aipe_cv(true_cv = .05, width = .0375, assurance = .50)
  without <- ss_aipe_cv(true_cv = .05, width = .0375)
  expect_equal(with_half$value, 17)
  expect_equal(without$value, 18)
})

test_that("ss_aipe_cv() plans past a noncentrality of 37.62 without a warning", {
  skip_on_cran()
  # Past a noncentrality of 37.62 R's noncentral t is an approximation;
  # the intervals and the quantile of the assurance step are exact there,
  # and nothing is signaled.
  res <- expect_no_warning(
    ss_aipe_cv(true_cv = .1, width = .1, conf_level = .99, assurance = .99))
  expect_equal(res$value, 32)
  expect_no_warning(ss_aipe_cv(true_cv = .3, width = .3))
})

test_that("ss_aipe_cv() searches past an infinite upper limit without a warning", {
  skip_on_cran()  # numerical oracle check, run in the full local suite
  # At 99 percent with a coefficient of .45 the search starts at N = 4,
  # where the upper limit of ci_cv()'s interval is infinite. That interval
  # is valid and ci_cv() signals nothing for it; its width is not narrow
  # enough, so the search goes on to a plan whose interval is finite.
  first <- expect_no_condition(
    ci_cv(cv = DMAR:::.cv_expected_value(.45, 4), N = 4, conf_level = .99))
  expect_identical(first$value[first$term == "ci_upper"], Inf)
  plan <- expect_no_condition(ss_aipe_cv(true_cv = .45, width = .10, conf_level = .99))
  expect_identical(plan$value, 387)
  at_plan <- ci_cv(cv = DMAR:::.cv_expected_value(.45, 387), N = 387, conf_level = .99)
  expect_true(is.finite(at_plan$value[at_plan$term == "ci_upper"]))
  expect_no_condition(ss_aipe_cv(true_cv = .50, width = .10, conf_level = .99, assurance = .90))
  # The search muffles nothing: a warning from the interval reaches the
  # user.
  local_mocked_bindings(ci_cv = function(cv, N, ...) {
    warning("another warning")
    data.frame(term = c("ci_lower", "ci_upper"), value = c(0, 1 / N))
  })
  expect_warning(ss_aipe_cv(true_cv = .2, width = .5), "another warning")
})

test_that("ss_aipe_cv() refuses alpha_upper = 0 at once, and plans alpha_lower = 0", {
  # With alpha_upper = 0 the interval is a lower confidence bound whose
  # upper limit is infinite at every N: ci_cv() gives Inf at N = 1000.
  expect_identical(ci_cv(cv = .2, N = 1000, conf_level = NULL, alpha_lower = .05,
                         alpha_upper = 0)$value[4], Inf)
  # The search would run its 100,000 steps (about 100 s) before stopping,
  # so the request is refused before the search, with or without
  # assurance.
  for (a in list(NULL, .9)) {
    t0 <- proc.time()[["elapsed"]]
    e <- tryCatch(ss_aipe_cv(true_cv = .2, width = .1, conf_level = NULL,
                             alpha_lower = .05, alpha_upper = 0, assurance = a),
                  error = function(e) e)
    expect_lt(proc.time()[["elapsed"]] - t0, 2)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, "alpha_upper")
    expect_match(conditionMessage(e), "infinite at every sample size", fixed = TRUE)
  }
  # With alpha_lower = 0 the lower limit is 0 and the width is the upper
  # bound itself, which comes down toward the coefficient of variation: a
  # width above it is planned, the smallest N whose bound at the expected
  # coefficient is no more than the width.
  plan <- ss_aipe_cv(true_cv = .2, width = .3, conf_level = NULL, alpha_lower = 0,
                     alpha_upper = .05)
  N <- plan$value
  up <- function(n) ci_cv(cv = DMAR:::.cv_expected_value(.2, n), N = n, conf_level = NULL,
                          alpha_lower = 0, alpha_upper = .05)$value
  expect_equal(up(N)[1], 0)
  expect_lte(up(N)[4], .3)
  expect_gt(up(N - 1)[4], .3)
})

test_that("ss_aipe_cv() refuses alpha_lower = 0 with a width at or below true_cv at once", {
  # The upper bound at the expected coefficient stays above true_cv at every
  # N the search visits (smallest 1.000006 times true_cv over the grid of
  # work/round6/evidence/g1_final_fixes/fix1/26), here at N = 4, 100, and
  # 100,004, although E[k_N] is below true_cv = .2.
  up <- function(n, au = .05) ci_cv(cv = DMAR:::.cv_expected_value(.2, n), N = n,
                                    conf_level = NULL, alpha_lower = 0,
                                    alpha_upper = au)$value[4]
  expect_lt(DMAR:::.cv_expected_value(.2, 100), .2)
  for (n in c(4, 100, 100004)) expect_gt(up(n), .2)
  # The search would run its 100,000 steps (87 s) and stop with
  # dmar_computation_failed, so the request is refused before the search.
  for (w in c(.2, .15)) for (a in list(NULL, .9)) {
    t0 <- proc.time()[["elapsed"]]
    e <- tryCatch(ss_aipe_cv(true_cv = .2, width = w, conf_level = NULL, alpha_lower = 0,
                             alpha_upper = .05, assurance = a), error = function(e) e)
    expect_lt(proc.time()[["elapsed"]] - t0, 2)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, c("width", "alpha_lower"))
    expect_match(conditionMessage(e), "is never reached", fixed = TRUE)
  }
  # The same through true_mean and true_sd.
  expect_error(ss_aipe_cv(true_mean = 10, true_sd = 2, width = .2, conf_level = NULL,
                          alpha_lower = 0, alpha_upper = .05),
               class = "dmar_argument_invalid")
  # With alpha_upper above .5 the bound can fall below true_cv (.996 of it
  # at N = 7 with alpha_upper = .55, and above it at N = 6), so that request
  # is still planned: 7.
  expect_lt(up(7, .55), .2)
  expect_gt(up(6, .55), .2)
  expect_identical(ss_aipe_cv(true_cv = .2, width = .2, conf_level = NULL, alpha_lower = 0,
                              alpha_upper = .55)$value, 7)
})

# ---- Assurance with a one-sided upper confidence bound --------------------

test_that("ss_aipe_cv() refuses assurance with alpha_lower = 0, and advises planning without it", {
  one_sided <- function(...) ss_aipe_cv(true_cv = .2, width = .25, conf_level = NULL,
                                       alpha_lower = 0, alpha_upper = .05, ...)
  # The assurance step of Kelley (2007) does not give the requested
  # probability for an upper bound, so the combination is refused for
  # every assurance.
  for (a in c(.50, .80, .90, .99)) {
    e <- tryCatch(one_sided(assurance = a), error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_s3_class(e, "dmar_refusal")
    expect_identical(e$arg, c("assurance", "alpha_lower"))
    expect_identical(e$value, a)
    expect_match(conditionMessage(e), "Plan without 'assurance'.", fixed = TRUE)
    # The only advice is to plan without assurance: a small positive lower
    # tail is a different interval, not a way around the refusal.
    expect_no_match(conditionMessage(e), "above 0", fixed = TRUE)
    expect_no_match(conditionMessage(e), "give 'alpha_lower'", fixed = TRUE)
  }
  # Also where the upper tail is above one half.
  expect_error(ss_aipe_cv(true_cv = .2, width = .3, conf_level = NULL, alpha_lower = 0,
                          alpha_upper = .55, assurance = .9),
               class = "dmar_argument_invalid")
  # Without assurance the bound is planned.
  expect_identical(one_sided()$value, 37)
  # The check of the assurance itself comes first, and so does the
  # refusal of a width that no bound reaches (tested above).
  e <- tryCatch(one_sided(assurance = 1.2), error = function(e) e)
  expect_s3_class(e, "dmar_assurance_invalid")
  e <- tryCatch(ss_aipe_cv(true_cv = .2, width = .2, conf_level = NULL, alpha_lower = 0,
                           alpha_upper = .05, assurance = .9), error = function(e) e)
  expect_identical(e$arg, c("width", "alpha_lower"))
  # The refusal is of a lower tail of exactly 0. Positive tails, however
  # unequal, are planned with assurance.
  expect_identical(ss_aipe_cv(true_cv = .25, width = .10, conf_level = NULL, alpha_lower = .01,
                              alpha_upper = .04, assurance = .80)$value, 68)
  expect_identical(ss_aipe_cv(true_cv = .2, width = .25, conf_level = NULL, alpha_lower = 1e-9,
                              alpha_upper = .05, assurance = .80)$value, 14)
})

# ---- The probability the plans attain, computed apart from the package ----

# The probability that the interval of ci_cv() is no wider than `width`
# when the population coefficient is `kappa` and the sample size `N`. The
# width increases with the sample coefficient k, so the interval is narrow
# enough exactly when k is at most k*, the coefficient whose interval has
# that width, which is when sqrt(N) / k is at least sqrt(N) / k*. That
# statistic is noncentral t on N - 1 degrees of freedom with noncentrality
# sqrt(N) / kappa. The probability is taken from the oracle of
# helper-nc_t_oracle.R, which integrates over the normal variate in base R
# and shares no code with the package; the limits at k* are checked by the
# same oracle, so that the interval the plan is for is right as well.
attained_probability <- function(kappa, width, N, alpha_lower, alpha_upper = alpha_lower) {
  width_at <- function(k) {
    ci <- ci_cv(cv = k, N = N, conf_level = NULL, alpha_lower = alpha_lower,
                alpha_upper = alpha_upper)
    ci$value[ci$term == "ci_upper"] - ci$value[ci$term == "ci_lower"]
  }
  k_star <- stats::uniroot(function(k) width_at(k) - width,
                           c(kappa / 50, kappa * 3), tol = 1e-13)$root
  ci <- ci_cv(cv = k_star, N = N, conf_level = NULL, alpha_lower = alpha_lower,
              alpha_upper = alpha_upper)
  t_star <- sqrt(N) / k_star
  expect_equal(.oracle_nct(t_star, N - 1, sqrt(N) / ci$value[1], lower_tail = TRUE),
               alpha_lower, tolerance = 1e-8)
  expect_equal(.oracle_nct(t_star, N - 1, sqrt(N) / ci$value[4], lower_tail = FALSE),
               alpha_upper, tolerance = 1e-8)
  .oracle_nct(t_star, N - 1, sqrt(N) / kappa, lower_tail = FALSE)
}

test_that("the plans of the worked example attain .5307 and .9966, computed apart from the package", {
  # Kelley (2007, pp. 759-760): a coefficient of .25, a width of .10, 99
  # percent confidence. The plan without assurance, 101, makes the interval
  # narrow enough a little more than half the time, and the plan for an
  # assurance of .99, 142, with a probability above .99.
  plain <- ss_aipe_cv(true_cv = .25, width = .10, conf_level = .99)$value
  assured <- ss_aipe_cv(true_cv = .25, width = .10, conf_level = .99, assurance = .99)$value
  expect_identical(c(plain, assured), c(101, 142))
  p_plain <- attained_probability(.25, .10, plain, .005)
  p_assured <- attained_probability(.25, .10, assured, .005)
  expect_equal(p_plain, .5307, tolerance = 5e-5)
  expect_equal(p_assured, .9966, tolerance = 5e-5)
  expect_gte(p_assured, .99)
  # The probability rises with the sample size, so one case fewer attains
  # less.
  expect_lt(attained_probability(.25, .10, 141, .005), p_assured)
  # The same event by simulation of raw scores through ci_cv(): 2,000
  # samples of 101, all with an interval, and a standard error of .011.
  set.seed(113)
  narrow <- vapply(seq_len(2000), function(i) {
    ci <- ci_cv(data = rnorm(101, 10, 2.5), conf_level = .99)
    ci$value[4] - ci$value[1] <= .10
  }, logical(1))
  expect_lt(abs(mean(narrow) - p_plain), 4 * sqrt(p_plain * (1 - p_plain) / 2000))
})

test_that("a plan of the assurance step can attain less than the request, as the help says", {
  # Table 1 of Kelley (2007, p. 760), 90 percent confidence, assurance
  # .80, width .10, coefficient .05: the function returns 5 where the
  # table prints 6. At 5 the probability is .7436, below the .80
  # requested (the procedure is an approximation and not an exact
  # search), and at 6 it is above .80.
  N <- ss_aipe_cv(true_cv = .05, width = .10, conf_level = .90, assurance = .80)$value
  expect_identical(N, 5)
  expect_equal(attained_probability(.05, .10, 5, .05), .7436, tolerance = 5e-5)
  expect_gt(attained_probability(.05, .10, 6, .05), .80)
  # With unequal positive tails the step attains about the request.
  expect_equal(attained_probability(.25, .10, 68, .01, .04), .8073, tolerance = 5e-5)
})

test_that("ss_aipe_cv() gives 93 where Table 1 of Kelley (2007) prints 127", {
  # Table 1 (90 percent confidence), assurance .80, width .125,
  # coefficient .40 prints 127, between 69 and 122 in its row: a misprint.
  row <- vapply(c(.35, .40, .45), function(k)
    ss_aipe_cv(true_cv = k, width = .125, conf_level = .90, assurance = .80)$value, 0)
  expect_identical(row, c(69, 93, 122))
  expect_gte(attained_probability(.40, .125, 93, .05), .80)
})

# ---- Input checks ----------------------------------------------------------

test_that("ss_aipe_cv() refuses tail areas that are missing, NA, or not one number", {
  tails <- function(...) ss_aipe_cv(true_cv = .2, width = .25, conf_level = NULL, ...)
  refused <- function(expr, arg) {
    e <- tryCatch(expr, error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, arg)
    expect_match(conditionMessage(e), "at least 0 and below 1 when 'conf_level' is NULL",
                 fixed = TRUE)
  }
  refused(tails(alpha_upper = .05), "alpha_lower")
  refused(tails(alpha_lower = .05), "alpha_upper")
  refused(tails(), "alpha_lower")
  refused(tails(alpha_lower = NA, alpha_upper = .05), "alpha_lower")
  refused(tails(alpha_lower = .05, alpha_upper = NA_real_), "alpha_upper")
  refused(tails(alpha_lower = c(.01, .02), alpha_upper = .05), "alpha_lower")
  refused(tails(alpha_lower = .05, alpha_upper = "0.05"), "alpha_upper")
  refused(tails(alpha_lower = -.01, alpha_upper = .05), "alpha_lower")
  refused(tails(alpha_lower = .05, alpha_upper = 1), "alpha_upper")
  # The function keeps what it accepted: both tails with a NULL level, a
  # lower tail of 0, and a level alone.
  expect_identical(tails(alpha_lower = .025, alpha_upper = .025)$value,
                   ss_aipe_cv(true_cv = .2, width = .25)$value)
  expect_identical(tails(alpha_lower = 0, alpha_upper = .05)$value, 37)
  expect_null(attr(tails(alpha_lower = .01, alpha_upper = .04), "conf_level"))
  expect_error(ss_aipe_cv(true_cv = .2, width = .25, alpha_lower = .01, alpha_upper = .04),
               "conf_level = NULL", fixed = TRUE, class = "dmar_argument_invalid")
})

test_that("ss_aipe_cv() refuses a planning value that is not one positive number", {
  refused <- function(expr, arg) {
    e <- tryCatch(expr, error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, arg)
  }
  refused(ss_aipe_cv(true_cv = c(.2, .3), width = .1), "true_cv")
  refused(ss_aipe_cv(true_cv = Inf, width = .1), "true_cv")
  refused(ss_aipe_cv(true_cv = 0, width = .1), "true_cv")
  # A negative mean with a negative standard deviation has a positive
  # ratio; each of the two must be positive on its own.
  refused(ss_aipe_cv(true_mean = -10, true_sd = -1, width = .1), "true_mean")
  refused(ss_aipe_cv(true_mean = 10, true_sd = -1, width = .1), "true_sd")
  refused(ss_aipe_cv(true_mean = 0, true_sd = 1, width = .1), "true_mean")
  refused(ss_aipe_cv(true_mean = 10, true_sd = 0, width = .1), "true_sd")
  refused(ss_aipe_cv(true_mean = c(10, 20), true_sd = 1, width = .1), "true_mean")
  expect_identical(ss_aipe_cv(true_mean = 10, true_sd = 1, width = .1)$value,
                   ss_aipe_cv(true_cv = .1, width = .1)$value)
})

# ---- The sensitivity function at a level other than the default -----------

test_that("ss_aipe_cv_sensitivity() plans and forms its intervals at the level it is given", {
  # The worked example at 99 percent: the plans are 101 and, with an
  # assurance of .99, 142 (at 95 percent they are 59 and 92).
  out <- ss_aipe_cv_sensitivity(true_cv = .25, supposed_cv = .25, width = .10,
                                conf_level = .99, G = 2, seed = 1)
  value <- function(x, t) x$value[x$term == t]
  expect_identical(value(out, "total_N"), 101)
  expect_identical(value(out, "conf_level"), .99)
  assured <- ss_aipe_cv_sensitivity(true_cv = .25, supposed_cv = .25, width = .10,
                                    conf_level = .99, assurance = .99, G = 2, seed = 1)
  expect_identical(value(assured, "total_N"), 142)
  expect_identical(value(assured, "assurance"), .99)
  # The intervals inside the replications are at that level too: the mean
  # width equals a loop over the same draws at 99 percent, and not the one
  # at 95 percent.
  hand <- function(level) {
    set.seed(1)
    mean(vapply(1:2, function(i) {
      ci <- ci_cv(data = rnorm(101, mean = 100, sd = 25), conf_level = level)
      ci$value[4] - ci$value[1]
    }, numeric(1)))
  }
  expect_equal(value(out, "mean_ci_width"), hand(.99), tolerance = 1e-12)
  expect_equal(value(out, "mean_ci_width"), 0.09116115, tolerance = 1e-7)
  expect_gt(abs(value(out, "mean_ci_width") - hand(.95)), .02)
})
