test_that("var_smd_trimmed() returns a one-row tidy data.frame", {
  res <- var_smd_trimmed(population_smd_trimmed = 0.5, n_1 = 30, n_2 = 30)
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
