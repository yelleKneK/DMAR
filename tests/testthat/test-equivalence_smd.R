test_that("equivalence_smd() returns documented rows", {
  res <- equivalence_smd(smd = 0.05, n_1 = 40, n_2 = 40, smd_upper = 0.4)
  expect_setequal(res$term,
                  c("smd", "t_lower", "t_upper", "df",
                    "p_lower", "p_upper", "p_tost",
                    "ci_lower", "ci_upper",
                    "smd_lower", "smd_upper", "equivalent"))
})

test_that("equivalence_smd() declares equivalence for a tiny d and large n", {
  res <- equivalence_smd(smd = 0.02, n_1 = 400, n_2 = 400, smd_upper = 0.3)
  expect_equal(res$value[res$term == "equivalent"], 1)
})

test_that("equivalence_smd() does not declare equivalence for a large d", {
  res <- equivalence_smd(smd = 0.7, n_1 = 50, n_2 = 50, smd_upper = 0.3)
  expect_equal(res$value[res$term == "equivalent"], 0)
})

test_that("equivalence_smd() raw-data interface agrees with summary interface", {
  set.seed(113)
  x <- rnorm(40, 100, 15); y <- rnorm(40, 101, 15)
  s_p <- sqrt(((40 - 1) * var(x) + (40 - 1) * var(y)) / (40 + 40 - 2))
  d_obs <- (mean(x) - mean(y)) / s_p
  r1 <- equivalence_smd(x = x, y = y, smd_upper = 0.3)
  r2 <- equivalence_smd(smd = d_obs, n_1 = 40, n_2 = 40, smd_upper = 0.3)
  expect_equal(r1$value[r1$term == "smd"],
               r2$value[r2$term == "smd"], tolerance = 1e-10)
  expect_equal(r1$value[r1$term == "p_tost"],
               r2$value[r2$term == "p_tost"], tolerance = 1e-10)
})

test_that("equivalence_smd() errors on bad bounds", {
  expect_error(equivalence_smd(smd = 0, n_1 = 30, n_2 = 30, smd_upper = -1),
               "positive number")
})

test_that("results_sentence() reports the TOST p-value, on the table and on its tidy()", {
  # The 1.1.0 review: the sentence reported p_lower (here 3e-08) beside the
  # 90% interval of a result that is not equivalent (TOST p = .081).
  s <- equivalence_smd(smd = .3, n_1 = 100, n_2 = 100, smd_lower = .5,
                       smd_upper = .5)
  p <- setNames(s$value, s$term)
  expect_lt(p[["p_lower"]], 1e-4)
  want <- "smd = 0.30, 90% CI [0.06, 0.54], TOST p = 0.0806"
  expect_identical(want, paste0("smd = 0.30, 90% CI [0.06, 0.54], TOST p = ",
                                format_p(p[["p_tost"]])))
  expect_identical(results_sentence(s), want)
  expect_identical(results_sentence(generics::tidy(s)), want)
})

test_that("the interval is d -/+ qt(1 - alpha, n_1 + n_2 - 2) SE(d), the one the planner plans for", {
  # Independent algebra: the large-sample standard error of d with a Student
  # t critical value, at confidence level 1 - 2 alpha. ss_aipe_equivalence_smd()
  # and ss_aipe_equivalence_smd_sensitivity() plan for and simulate this
  # interval (their tests check that against this function).
  grid <- expand.grid(smd = c(-1.2, 0, 0.35, 2), n_1 = c(2, 9, 120),
                      n_2 = c(3, 9, 400), alpha = c(0.05, 0.01))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    se <- sqrt((g$n_1 + g$n_2) / (g$n_1 * g$n_2) +
                 g$smd^2 / (2 * (g$n_1 + g$n_2)))
    crit <- stats::qt(1 - g$alpha, g$n_1 + g$n_2 - 2)
    r <- equivalence_smd(smd = g$smd, n_1 = g$n_1, n_2 = g$n_2,
                         smd_lower = 0.4, smd_upper = 0.6,
                         alpha_level = g$alpha)
    lo <- r$value[r$term == "ci_lower"]; hi <- r$value[r$term == "ci_upper"]
    expect_equal(lo, g$smd - crit * se, tolerance = 1e-14)
    expect_equal(hi, g$smd + crit * se, tolerance = 1e-14)
    expect_equal(attr(r, "conf_level"), 1 - 2 * g$alpha)
    # Both one-sided tests reject at alpha exactly when the interval lies
    # inside the bounds.
    both <- max(r$value[r$term == "p_lower"], r$value[r$term == "p_upper"]) <
      g$alpha
    expect_identical(both, lo > -0.4 && hi < 0.6)
    expect_identical(r$value[r$term == "equivalent"], as.numeric(both))
  }
})
