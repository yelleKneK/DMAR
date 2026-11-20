## Tests for proportion_of_superiority() (sometimes called Cohen's U3).

test_that("proportion_of_superiority() point estimate equals pnorm(d)", {
  for (d in c(0.2, 0.5, 0.8, -0.5)) {
    res <- proportion_of_superiority(d)
    expect_equal(res$value[res$term == "proportion_of_superiority"],
                 pnorm(d), tolerance = 1e-12)
  }
})

test_that("proportion_of_superiority() values match Cohen (1988) Table 2.2.1", {
  # Cohen (1988) Table 2.2.1: d = 0.2, 0.5, 0.8 -> 58%, 69%, 79%.
  expect_equal(round(proportion_of_superiority(0.2)$value[2], 2), 0.58)
  expect_equal(round(proportion_of_superiority(0.5)$value[2], 2), 0.69)
  expect_equal(round(proportion_of_superiority(0.8)$value[2], 2), 0.79)
})

test_that("proportion_of_superiority() CI from sample sizes is monotone", {
  res <- proportion_of_superiority(smd = 0.5, n_1 = 50, n_2 = 50)
  expect_gt(res$value[res$term == "ci_upper"],
            res$value[res$term == "ci_lower"])
})

test_that("proportion_of_superiority() agrees with cles relationship for positive d", {
  # The proportion of superiority is pnorm(d); CLES is pnorm(d / sqrt(2)).
  # For positive d the proportion of superiority strictly exceeds CLES.
  ps <- proportion_of_superiority(0.5)$value[2]
  cl <- cles(0.5)$value[cles(0.5)$term == "cles"]
  expect_gt(ps, cl)
})

test_that("proportion_of_superiority() returns a tidy numeric data.frame", {
  res <- proportion_of_superiority(smd = 0.5, n_1 = 50, n_2 = 50)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_type(res$value, "double")
  expect_true(all(c("smd", "proportion_of_superiority",
                    "smd_ci_lower", "smd_ci_upper",
                    "ci_lower", "ci_upper") %in% res$term))
})

test_that("proportion_of_superiority() names its own row as its interval's estimate", {
  # The 1.1.0 review: results_sentence() took the first row, the
  # standardized mean difference, as the estimate of the interval for the
  # proportion of superiority, and wrote "smd = 0.50, 95% CI [0.54, 0.82]",
  # an interval that does not contain the estimate it names.
  res <- proportion_of_superiority(smd = .5, n_1 = 50, n_2 = 50)
  expect_identical(attr(res, "ci_estimate"), "proportion_of_superiority")
  expect_identical(results_sentence(res),
                   "proportion_of_superiority = 0.69, 95% CI [0.54, 0.82]")
})

test_that("proportion_of_superiority() takes and reports the limits on d as smd_ci_lower and smd_ci_upper", {
  # DMAR 1.0.x called the limits smd_lower and smd_upper; those names are
  # the equivalence bounds of equivalence_smd() and stop here with
  # "unused argument".
  res <- proportion_of_superiority(smd = 0.5, smd_ci_lower = 0.2,
                                   smd_ci_upper = 0.8)
  expect_identical(res$term, c("smd", "proportion_of_superiority",
                               "smd_ci_lower", "smd_ci_upper",
                               "ci_lower", "ci_upper"))
  expect_identical(res$value[3:6], c(0.2, 0.8, pnorm(0.2), pnorm(0.8)))
  expect_error(proportion_of_superiority(smd = 0.5, smd_lower = 0.2,
                                         smd_upper = 0.8),
               "unused argument")
})

test_that("tidy() of proportion_of_superiority() gives the level on its own interval, and the sentence follows its rows", {
  # The limits on d are named for d (smd_ci_lower, smd_ci_upper); the
  # table names d as the source of its interval (ci_source), so tidy()
  # gives the level on the proportion of superiority and its limits alone,
  # as in DMAR 1.0.x. Oracle: pnorm() of d and of its limits, computed
  # from the noncentral t inversion with uniroot() and pt().
  n_1 <- 50; n_2 <- 40; d <- 0.5
  scale <- sqrt(n_1 * n_2 / (n_1 + n_2))
  ncp_at <- function(p) uniroot(function(ncp) pt(d * scale, n_1 + n_2 - 2, ncp) - p,
                                c(-3, 8), tol = 1e-12)$root
  ps <- pnorm(c(d, c(ncp_at(.975), ncp_at(.025)) / scale))
  want <- sprintf("proportion_of_superiority = %.2f, 95%% CI [%.2f, %.2f]",
                  ps[1], ps[2], ps[3])
  x <- proportion_of_superiority(smd = d, n_1 = n_1, n_2 = n_2)
  td <- generics::tidy(x)
  expect_identical(td$conf_level, c(NA, 0.95, NA, NA, 0.95, 0.95))
  expect_identical(attr(x, "ci_source"), "smd")
  expect_identical(results_sentence(x), want)
  expect_identical(results_sentence(td), want)
  expect_identical(results_sentence(rbind(td, td)[seq_len(nrow(td)), ]), want)
  expect_identical(results_sentence(generics::glance(x)), want)
})
