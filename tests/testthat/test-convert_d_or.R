test_that("convert_d_or() and convert_or_d() match the Hasselblad-Hedges forms", {
  expect_equal(convert_d_or(0.5)$value, exp(0.5 * pi / sqrt(3)))
  expect_equal(convert_or_d(2)$value, log(2) * sqrt(3) / pi)
  expect_identical(convert_d_or(0.5)$term, "odds_ratio")
  expect_identical(convert_or_d(2)$term, "smd")
  # The null maps to the null.
  expect_equal(convert_d_or(0)$value, 1)
  expect_equal(convert_or_d(1)$value, 0)
})

test_that("convert_d_or() / convert_or_d() are exact inverses (property test)", {
  set.seed(113)
  for (d in c(-1.5, -0.2, 0, 0.3, 1, runif(20, -2, 2))) {
    expect_equal(convert_or_d(convert_d_or(d)$value)$value, d, tolerance = 1e-12)
  }
  for (or in c(0.1, 0.5, 1, 2.5, 10, exp(runif(20, -2, 2)))) {
    expect_equal(convert_d_or(convert_or_d(or)$value)$value, or, tolerance = 1e-12)
  }
})

test_that("convert_d_or family validates its arguments", {
  expect_error(convert_d_or("x"), "single number")
  expect_error(convert_or_d(0), "positive")
  expect_error(convert_or_d(-2), "positive")
})
