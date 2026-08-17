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
