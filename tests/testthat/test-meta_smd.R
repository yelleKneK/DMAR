v <- function(tab, t) tab$value[tab$term == t]

test_that("meta_smd() matches metafor's escalc(measure='SMD') pipeline", {
  d  <- teacher_expectancy$d
  ne <- teacher_expectancy$n_experimental
  nc <- teacher_expectancy$n_control

  ours <- meta_smd(smd = d, n_1 = ne, n_2 = nc, hartung_knapp = FALSE)

  # metafor's SMD route: g = J * d with J = 1 - 3/(4m - 1) (their
  # approximation); DMAR uses the exact gamma-function J, so allow the
  # tiny approximation gap at these sample sizes.
  # Pinned from metafor::rma() on the exact-J g values and their
  # large-sample variances (metafor 5.0.1, 2026-08-09); live comparison
  # in tools/oracle_checks.R.
  expect_equal(v(ours, "estimate"), 0.05438657047124985, tolerance = 1e-3)
  expect_lt(abs(v(ours, "tau2") - 3.883947625026544e-07), 1e-4)
})

test_that("meta_smd(unbiased = FALSE) pools the raw d values", {
  d  <- teacher_expectancy$d
  ne <- teacher_expectancy$n_experimental
  nc <- teacher_expectancy$n_control
  raw <- meta_smd(d, ne, nc, unbiased = FALSE, hartung_knapp = FALSE)
  via_es <- meta_es(d, (ne + nc) / (ne * nc) + d^2 / (2 * (ne + nc)),
                    hartung_knapp = FALSE)
  expect_equal(v(raw, "estimate"), v(via_es, "estimate"))
  expect_equal(v(raw, "tau2"), v(via_es, "tau2"))
  # The correction shrinks every study, so the pooled estimate shrinks too.
  cor <- meta_smd(d, ne, nc, unbiased = TRUE, hartung_knapp = FALSE)
  expect_lt(abs(v(cor, "estimate")), abs(v(raw, "estimate")) + 1e-12)
})

test_that("meta_smd() validates its arguments", {
  expect_error(meta_smd(c(.2, .3), n_1 = c(10, 20), n_2 = c(10, 1.5)),
               "integer sample size")
  expect_error(meta_smd(c(.2, .3), n_1 = 10, n_2 = c(10, 20)),
               "integer sample size")
  expect_error(meta_smd(.2, 10, 10), "two or more")
})

test_that("meta_smd() warns when the studies are identical under Hartung-Knapp", {
  # The same standardized mean difference from the same group sizes gives
  # the same bias-corrected effect size in every study.
  expect_warning(res <- meta_smd(smd = rep(.4, 3), n_1 = rep(30, 3),
                                 n_2 = rep(30, 3)),
                 class = "dmar_meta_hk_zero_se")
  # Zero or numerically zero: whether the weighted mean of three equal
  # numbers reproduces that number to the last bit depends on the number.
  expect_lt(abs(v(res, "se")), 1e-14)
  expect_equal(v(res, "ci_lower"), v(res, "estimate"))
  # The warning's conventional interval is the one hartung_knapp = FALSE
  # reports, and the table carries q and the conventional standard error.
  w <- expect_warning(meta_smd(smd = rep(.4, 3), n_1 = rep(30, 3),
                               n_2 = rep(30, 3)),
                      class = "dmar_meta_hk_zero_se")
  plain <- meta_smd(smd = rep(.4, 3), n_1 = rep(30, 3), n_2 = rep(30, 3),
                    hartung_knapp = FALSE)
  expect_identical(w$n_studies, 3L)
  expect_equal(w$estimate, v(res, "estimate"))
  expect_equal(w$se_conventional, v(plain, "se"))
  expect_equal(c(w$lower_conventional, w$upper_conventional),
               c(v(plain, "ci_lower"), v(plain, "ci_upper")))
  # q is zero or numerically zero, as the standard error is: with the
  # correction factor computed without cancellation the common effect
  # size is not reproduced to the last bit by its weighted mean.
  expect_gte(attr(res, "hartung_knapp_q"), 0)
  expect_lt(attr(res, "hartung_knapp_q"), 1e-20)
  expect_equal(attr(res, "se_conventional"), v(plain, "se"))
  expect_no_warning(meta_smd(smd = rep(.4, 3), n_1 = rep(30, 3),
                             n_2 = rep(30, 3), hartung_knapp = FALSE))
})

test_that("meta_smd() refuses an infinite sample size in its own words", {
  expect_error(meta_smd(c(0.1, 0.2, 0.3), c(20, Inf, 20), rep(20, 3)),
               "'n_1' must give an integer sample size")
  expect_error(meta_smd(c(0.1, 0.2, 0.3), rep(20, 3), c(20, Inf, 20)),
               "'n_2' must give an integer sample size")
})
