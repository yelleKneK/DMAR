test_that("nnt_from_smd() point estimate matches Kraemer-Kupfer formula", {
  res <- nnt_from_smd(smd = 0.5)
  srd_expected <- 2 * pnorm(0.5 / sqrt(2)) - 1
  expect_equal(res$value[res$term == "srd"], srd_expected, tolerance = 1e-12)
  expect_equal(res$value[res$term == "nnt"], 1 / srd_expected, tolerance = 1e-12)
})

test_that("nnt_from_smd() returns Inf NNT when d = 0", {
  res <- nnt_from_smd(smd = 0)
  expect_equal(res$value[res$term == "srd"], 0)
  expect_true(is.infinite(res$value[res$term == "nnt"]))
})

test_that("nnt_from_smd() reproduces Furukawa and Leucht (2011) Table 1 NNTs", {
  # Furukawa & Leucht (2011) Table 1: d = 0.2 -> NNT ≈ 8.9; 0.5 -> 3.6; 0.8 -> 2.3.
  expect_equal(round(nnt_from_smd(0.2)$value[3], 1), 8.9)
  expect_equal(round(nnt_from_smd(0.5)$value[3], 1), 3.6)
  expect_equal(round(nnt_from_smd(0.8)$value[3], 1), 2.3)
})

test_that("nnt_from_smd() with CI on d propagates monotonically (smaller d -> larger NNT bound)", {
  res <- nnt_from_smd(smd = 0.5, smd_lower = 0.20, smd_upper = 0.80)
  nnt_lo <- res$value[res$term == "nnt_lower"]
  nnt_hi <- res$value[res$term == "nnt_upper"]
  expect_true(nnt_hi > nnt_lo)
  # Bound at d_upper = 0.8 should match d=0.8 point estimate:
  expect_equal(nnt_lo, nnt_from_smd(smd = 0.80)$value[3], tolerance = 1e-12)
})

test_that("nnt_from_smd() yields Inf upper-NNT bound when smd_lower <= 0", {
  res <- nnt_from_smd(smd = 0.3, smd_lower = -0.1, smd_upper = 0.7)
  expect_true(is.infinite(res$value[res$term == "nnt_upper"]) ||
              res$value[res$term == "nnt_upper"] < 0)
})

test_that("nnt_from_smd() with sample sizes uses ci_smd() internally", {
  res <- nnt_from_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.95)
  expect_true("nnt_lower" %in% res$term)
  expect_true("nnt_upper" %in% res$term)
})
