test_that("var_smd_trimmed() returns a one-row tidy data.frame", {
  res <- var_smd_trimmed(true_smd_trimmed = 0.5, n_1 = 30, n_2 = 30)
  expect_equal(nrow(res), 1)
  expect_equal(res$term, "var_smd_trimmed")
  expect_gt(res$value, 0)
})

test_that("var_smd_trimmed() trim = 0 reduces to var_smd-like value", {
  v_t <- var_smd_trimmed(0.5, 30, 30, trim = 0)$value
  v_s <- (30 + 30) / (30 * 30) + 0.5^2 / (2 * (30 + 30))
  expect_equal(v_t, v_s, tolerance = 1e-12)
})

test_that("var_smd_trimmed() trim > 0 increases variance (smaller h)", {
  v0 <- var_smd_trimmed(0.5, 30, 30, trim = 0)$value
  v2 <- var_smd_trimmed(0.5, 30, 30, trim = 0.20)$value
  expect_gt(v2, v0)
})

test_that("var_smd_trimmed() rejects invalid n / trim", {
  expect_error(var_smd_trimmed(0.5, 3, 30), ">= 4")
  expect_error(var_smd_trimmed(0.5, 30, 30, trim = 0.6), "in \\[0, 0.5\\)")
})

test_that("var_smd_trimmed() takes true_smd_trimmed; the DMAR 1.0.x name stops", {
  expect_identical(names(formals(var_smd_trimmed)),
                   c("true_smd_trimmed", "n_1", "n_2", "trim"))
  expect_error(var_smd_trimmed(population_smd_trimmed = 0.5, n_1 = 30, n_2 = 30),
               "unused argument")
})

test_that("var_smd_trimmed() evaluates the variance on the help page", {
  # By hand from the formula on ?var_smd_trimmed: with 20 percent trimming
  # of 30 per group, h = 30 - 2 * 6 = 18 in each group, so the variance is
  # 36 / 324 + 0.25 / 72; with 10 percent trimming of 25 and 41, h = 21 and
  # 33, and at 0.8 the variance is 54 / 693 + 0.64 / 108.
  expect_equal(var_smd_trimmed(true_smd_trimmed = 0.5, n_1 = 30, n_2 = 30)$value,
               36 / 324 + 0.25 / 72, tolerance = 1e-14)
  expect_equal(var_smd_trimmed(true_smd_trimmed = 0.8, n_1 = 25, n_2 = 41,
                               trim = 0.1)$value,
               54 / 693 + 0.64 / 108, tolerance = 1e-14)
})
