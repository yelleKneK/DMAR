test_that("cles() point estimate matches McGraw-Wong formula", {
  expect_equal(cles(smd = 0)$value[cles(smd = 0)$term == "cl"], 0.5,
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
  res <- cles(smd = 0.5, smd_lower = 0.20, smd_upper = 0.80)
  cl_lo <- res$value[res$term == "cl_lower"]
  cl_hi <- res$value[res$term == "cl_upper"]
  expect_true(cl_hi > cl_lo)
  expect_equal(cl_lo, pnorm(0.20 / sqrt(2)), tolerance = 1e-12)
})

test_that("cles() with sample sizes uses ci_smd() internally", {
  res <- cles(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.95)
  expect_true("cl_lower" %in% res$term)
  expect_true("cl_upper" %in% res$term)
})
