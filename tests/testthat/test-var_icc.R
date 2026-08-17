# Tests for var_icc() against the closed-form Smith (1956) / Spearman-Brown
# delta method expressions.

test_that("var_icc() ICC(1,1) reduces to 0 at rho = 1 and is positive otherwise", {
  expect_equal(var_icc(rho = 1, n = 30, k = 4, type = "ICC(1,1)")$value, 0)
  expect_gt(var_icc(rho = 0.6, n = 30, k = 4, type = "ICC(1,1)")$value, 0)
})

test_that("var_icc() ICC(1,1) matches Smith (1956) hand calculation", {
  # 2 (1 - .6)^2 (1 + 3 * .6)^2 / (30 * 4 * 3)
  expected <- 2 * (1 - 0.6)^2 * (1 + 3 * 0.6)^2 / (30 * 4 * 3)
  expect_equal(var_icc(rho = 0.6, n = 30, k = 4, type = "ICC(1,1)")$value,
               expected, tolerance = 1e-12)
})

test_that("var_icc() ICC(1,k) closed-form matches the delta method derivation", {
  # For rho_k = .857 at the average-of-4 level, the closed form is
  # 2 k (1 - rho_k)^2 / (n (k - 1)).
  expected <- 2 * 4 * (1 - 0.857)^2 / (30 * 3)
  expect_equal(var_icc(rho = 0.857, n = 30, k = 4, type = "ICC(1,k)")$value,
               expected, tolerance = 1e-12)
})

test_that("var_icc() ICC(1,k) reduces correctly via Spearman-Brown to ICC(1,1)", {
  # If rho_k is obtained from rho_single = .6 by Spearman-Brown then the
  # closed-form Var(rho_k) should equal the delta method transform of
  # Var(rho_single).
  rho_single <- 0.6; k <- 4; n <- 30
  rho_k <- k * rho_single / (1 + (k - 1) * rho_single)
  var_k <- var_icc(rho = rho_k, n = n, k = k, type = "ICC(1,k)")$value

  var_single <- var_icc(rho = rho_single, n = n, k = k, type = "ICC(1,1)")$value
  jac        <- k / (1 + (k - 1) * rho_single)^2
  expect_equal(var_k, jac^2 * var_single, tolerance = 1e-10)
})

test_that("var_icc() single-rater forms share the closed form across ICC types", {
  v_one   <- var_icc(rho = 0.5, n = 40, k = 3, type = "ICC(1,1)")$value
  v_two   <- var_icc(rho = 0.5, n = 40, k = 3, type = "ICC(2,1)")$value
  v_three <- var_icc(rho = 0.5, n = 40, k = 3, type = "ICC(3,1)")$value
  expect_equal(v_two,   v_one, tolerance = 1e-12)
  expect_equal(v_three, v_one, tolerance = 1e-12)
})

test_that("var_icc() accepts the icc() shorthand aliases", {
  v_full  <- var_icc(rho = 0.5, n = 40, k = 3, type = "ICC(1,1)")$value
  v_alias <- var_icc(rho = 0.5, n = 40, k = 3, type = "1")$value
  expect_equal(v_alias, v_full, tolerance = 1e-12)
})

test_that("var_icc() rejects invalid input", {
  expect_error(var_icc(rho = -0.1, n = 30, k = 4), "rho")
  expect_error(var_icc(rho = 1.1,  n = 30, k = 4), "rho")
  expect_error(var_icc(rho = 0.5,  n = 1,  k = 4), "n")
  expect_error(var_icc(rho = 0.5,  n = 30, k = 1), "k")
  expect_error(var_icc(rho = 0.5,  n = 30, k = 4, type = "ICC(4,1)"),
               "Unrecognized")
})

test_that("var_icc() variance scales like 1/n", {
  v_30  <- var_icc(rho = 0.4, n = 30,  k = 4)$value
  v_300 <- var_icc(rho = 0.4, n = 300, k = 4)$value
  expect_equal(v_30 / v_300, 10, tolerance = 1e-10)
})
