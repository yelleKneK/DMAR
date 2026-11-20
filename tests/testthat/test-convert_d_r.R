test_that("convert_d_r() and convert_r_d() match the textbook formulas", {
  # Equal groups (a = 4): d = 0.5 -> r = 0.5 / sqrt(0.25 + 4).
  expect_equal(convert_d_r(0.5)$value, 0.5 / sqrt(0.5^2 + 4))
  expect_equal(convert_r_d(0.243)$value, 2 * 0.243 / sqrt(1 - 0.243^2))
  expect_identical(convert_d_r(0.5)$term, "r")
  expect_identical(convert_r_d(0.2)$term, "smd")
})

test_that("convert_d_r() / convert_r_d() are exact inverses (property test)", {
  set.seed(113)
  for (d in c(-2, -0.5, 0, 0.1, 0.8, 3, runif(20, -2, 2))) {
    expect_equal(convert_r_d(convert_d_r(d)$value)$value, d, tolerance = 1e-12)
  }
  for (r in c(-0.9, -0.3, 0, 0.25, 0.7, runif(20, -0.95, 0.95))) {
    expect_equal(convert_d_r(convert_r_d(r)$value)$value, r, tolerance = 1e-12)
  }
})

test_that("the unequal-group factor is honored and round-trips", {
  a <- (20 + 80)^2 / (20 * 80)   # 6.25
  expect_equal(convert_d_r(0.5, n_1 = 20, n_2 = 80)$value,
               0.5 / sqrt(0.5^2 + a))
  d_back <- convert_r_d(convert_d_r(0.5, n_1 = 20, n_2 = 80)$value,
                        n_1 = 20, n_2 = 80)$value
  expect_equal(d_back, 0.5, tolerance = 1e-12)
  # Equal explicit n reproduces the default a = 4.
  expect_equal(convert_d_r(0.5, n_1 = 30, n_2 = 30)$value,
               convert_d_r(0.5)$value)
})

test_that("convert_d_r family validates its arguments", {
  expect_error(convert_d_r("x"), "single number")
  expect_error(convert_r_d(1), "\\(-1, 1\\)")
  expect_error(convert_d_r(0.5, n_1 = 20), "both")
  expect_error(convert_d_r(0.5, n_1 = 20, n_2 = 2.5), "positive integer")
})
