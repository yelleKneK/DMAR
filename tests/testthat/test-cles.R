test_that("cles() point estimate matches McGraw-Wong formula", {
  expect_equal(cles(smd = 0)$value[cles(smd = 0)$term == "cles"], 0.5,
               tolerance = 1e-12)
  expect_equal(cles(smd = 0.5)$value[2], pnorm(0.5 / sqrt(2)),
               tolerance = 1e-12)
})

test_that("cles() values match published reference (Brooks et al., 2014)", {
  # Standardized mean differences of 0.2 / 0.5 / 0.8 -> CL of ≈ 0.56 / 0.64 / 0.71
  expect_equal(round(cles(smd = 0.2)$value[2], 2), 0.56)
  expect_equal(round(cles(smd = 0.5)$value[2], 2), 0.64)
  expect_equal(round(cles(smd = 0.8)$value[2], 2), 0.71)
})

test_that("cles() is symmetric: CL(-d) = 1 - CL(d)", {
  expect_equal(cles(smd = -0.5)$value[2],
               1 - cles(smd = 0.5)$value[2],
               tolerance = 1e-12)
})

test_that("cles() with CI on d propagates monotonically", {
  res <- cles(smd = 0.5, smd_ci_lower = 0.20, smd_ci_upper = 0.80)
  cl_lo <- res$value[res$term == "cles_ci_lower"]
  cl_hi <- res$value[res$term == "cles_ci_upper"]
  expect_true(cl_hi > cl_lo)
  expect_equal(cl_lo, pnorm(0.20 / sqrt(2)), tolerance = 1e-12)
})

test_that("cles() with sample sizes uses ci_smd() internally", {
  res <- cles(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.95)
  expect_true("cles_ci_lower" %in% res$term)
  expect_true("cles_ci_upper" %in% res$term)
})

test_that("cles() takes and reports the limits on d as smd_ci_lower and smd_ci_upper", {
  # One name, one meaning: smd_lower and smd_upper are the equivalence
  # bounds of equivalence_smd(), so the confidence limits on d that cles()
  # takes and echoes are named for their estimate (DMAR 1.0.x called them
  # smd_lower and smd_upper). The old names stop with "unused argument".
  res <- cles(smd = 0.5, smd_ci_lower = 0.20, smd_ci_upper = 0.80)
  expect_identical(res$term, c("smd", "cles", "smd_ci_lower", "smd_ci_upper",
                               "cles_ci_lower", "cles_ci_upper"))
  expect_identical(res$value[3:4], c(0.20, 0.80))
  expect_equal(res$value[5:6], pnorm(c(0.20, 0.80) / sqrt(2)),
               tolerance = 1e-15)
  expect_error(cles(smd = 0.5, smd_lower = 0.20, smd_upper = 0.80),
               "unused argument")
  # The limits computed from the sample sizes are those of the noncentral
  # t interval: the t statistic d sqrt(n_1 n_2 / (n_1 + n_2)) and the ncp
  # at which it is the 97.5th and 2.5th percentile, found here by uniroot().
  n_1 <- 12; n_2 <- 15; d <- -0.84
  t_obs <- d * sqrt(n_1 * n_2 / (n_1 + n_2))
  df <- n_1 + n_2 - 2
  ncp_at <- function(p) uniroot(function(ncp) pt(t_obs, df, ncp) - p,
                                c(-20, 20), tol = 1e-12)$root
  d_lim <- c(ncp_at(.975), ncp_at(.025)) / sqrt(n_1 * n_2 / (n_1 + n_2))
  res <- cles(smd = d, n_1 = n_1, n_2 = n_2)
  expect_equal(res$value[res$term %in% c("smd_ci_lower", "smd_ci_upper")],
               d_lim, tolerance = 1e-7)
  expect_equal(res$value[res$term %in% c("cles_ci_lower", "cles_ci_upper")],
               pnorm(d_lim / sqrt(2)), tolerance = 1e-7)
  expect_false(any(c("smd_lower", "smd_upper") %in% res$term))
})

test_that("results_sentence() reports CL on the table, its tidy() rows, stacks and subsets of them, and its glance() row", {
  # The limits on d are now named for d (smd_ci_lower, smd_ci_upper), as
  # the limits on CL are named for CL, so the verbs see two intervals. The
  # table names d as the source of the interval for CL (ci_source), and the
  # sentence is CL's wherever the rows go, as it was in DMAR 1.0.x, when
  # the limits on d were smd_lower and smd_upper. Oracle: the limits on d
  # by inverting the noncentral t with uniroot() and pt(), and CL as
  # pnorm(d / sqrt(2)) at the estimate and at each limit.
  n_1 <- 50; n_2 <- 40; d <- 0.5
  scale <- sqrt(n_1 * n_2 / (n_1 + n_2))
  ncp_at <- function(p) uniroot(function(ncp) pt(d * scale, n_1 + n_2 - 2, ncp) - p,
                                c(-3, 8), tol = 1e-12)$root
  d_lim <- c(ncp_at(.975), ncp_at(.025)) / scale
  cl <- pnorm(c(d, d_lim) / sqrt(2))
  want <- sprintf("cles = %.2f, 95%% CI [%.2f, %.2f]", cl[1], cl[2], cl[3])
  x <- cles(smd = d, n_1 = n_1, n_2 = n_2)
  td <- generics::tidy(x)
  expect_identical(results_sentence(x), want)
  expect_identical(results_sentence(td), want)
  expect_identical(results_sentence(rbind(td, td)), want)
  expect_identical(results_sentence(rbind(td, td)[seq_len(nrow(td)), ]), want)
  expect_identical(results_sentence(td[-1, ]), want)
  expect_identical(results_sentence(generics::glance(x), conf_level = 0.95),
                   want)
  # Rows that have lost the ci_source attribute (a stack whose first frame
  # is another table's) are read from their conf_level column.
  bare <- td
  attr(bare, "ci_source") <- NULL
  expect_identical(results_sentence(bare), want)
  # tidy() gives the level on the rows of CL and its limits alone, as in
  # DMAR 1.0.x; the rows of d keep NA.
  expect_identical(td$conf_level, c(NA, 0.95, NA, NA, 0.95, 0.95))
  # The interval for d is reported when term names it.
  expect_identical(results_sentence(x, term = "smd"),
                   sprintf("smd = 0.50, 95%% CI [%.2f, %.2f]", d_lim[1],
                           d_lim[2]))
  # The same with limits on d supplied.
  y <- cles(smd = 0.3, smd_ci_lower = 0.1, smd_ci_upper = 0.7)
  cl <- pnorm(c(0.3, 0.1, 0.7) / sqrt(2))
  want <- sprintf("cles = %.2f, 95%% CI [%.2f, %.2f]", cl[1], cl[2], cl[3])
  expect_identical(results_sentence(y), want)
  expect_identical(results_sentence(generics::tidy(y)), want)
  expect_identical(attr(y, "ci_source"), "smd")
  expect_null(attr(cles(smd = 0.3), "ci_source"))
})

test_that("the common language effect size is the row cles, in cles() and randomization_test()", {
  # Ruling 21 of 2026-11-11: one name for the quantity in both functions.
  # DMAR 1.0.x named it cl (limits cl_lower, cl_upper) here and cles
  # (limits cles_lower_limit, cles_upper_limit) in randomization_test().
  # A lookup by an old name finds nothing.
  res <- cles(smd = 0.5, n_1 = 50, n_2 = 50)
  rows <- c("cles", "cles_ci_lower", "cles_ci_upper")
  expect_true(all(rows %in% res$term))
  expect_false(any(c("cl", "cl_lower", "cl_upper", "cl_ci_lower",
                     "cl_ci_upper") %in% res$term))
  expect_length(res$value[res$term == "cl"], 0L)
  # The broom verbs and the sentence carry the name.
  expect_true(all(rows %in% generics::tidy(res)$term))
  expect_true(all(rows %in% names(generics::glance(res))))
  expect_match(results_sentence(res), "^cles = 0\\.64, 95% CI")
  # randomization_test() reports the same rows, and its values are CL at
  # its own d and the limits of d's interval, Phi(d / sqrt(2)) (the
  # oracle), and those of cles() given that d and those limits.
  set.seed(113)
  y_1 <- rnorm(15, .5); y_2 <- rnorm(15)
  rt <- randomization_test(group_1 = y_1, group_2 = y_2, n_resamples = 999,
                           seed = 113)
  expect_true(all(rows %in% rt$term))
  expect_false(any(c("cl", "cles_lower_limit", "cles_upper_limit") %in% rt$term))
  v <- function(tab, t) tab$value[match(t, tab$term)]
  d3 <- v(rt, c("smd", "smd_ci_lower", "smd_ci_upper"))
  expect_equal(v(rt, rows), pnorm(d3 / sqrt(2)), tolerance = 1e-14)
  same <- cles(smd = d3[1], smd_ci_lower = d3[2], smd_ci_upper = d3[3])
  expect_identical(v(same, rows), v(rt, rows))
})
