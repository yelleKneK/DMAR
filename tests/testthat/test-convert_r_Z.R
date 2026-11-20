test_that("convert_r_Z() returns a tidy 1-row data.frame", {
  res <- convert_r_Z(0.5)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("term", "value"))
  expect_equal(res$term, "Z_from_r")
})

test_that("convert_r_Z() matches the closed-form definition Z = atanh(r)", {
  for (r in c(-0.9, -0.5, -0.1, 0, 0.1, 0.5, 0.9)) {
    expect_equal(convert_r_Z(r)$value, atanh(r), tolerance = 1e-12)
  }
})

test_that("convert_r_Z() Hays (1994) example: r = .35 yields Z = .3654", {
  expect_equal(round(convert_r_Z(.35)$value, 4), 0.3654)
})

test_that("convert_Z_r() returns a tidy 1-row data.frame", {
  res <- convert_Z_r(0.5)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("term", "value"))
  expect_equal(res$term, "r_from_Z")
})

test_that("convert_Z_r() inverts convert_r_Z() to floating-point precision", {
  for (r in c(-0.9, -0.5, -0.1, 0.1, 0.5, 0.9)) {
    Z <- convert_r_Z(r)$value
    expect_equal(convert_Z_r(Z)$value, r, tolerance = 1e-12)
  }
})

test_that("convert_Z_r() matches the closed-form definition r = tanh(Z)", {
  for (Z in c(-2, -1, -0.5, 0, 0.5, 1, 2)) {
    expect_equal(convert_Z_r(Z)$value, tanh(Z), tolerance = 1e-12)
  }
})

test_that("convert_r_Z() alias is identical to convert_r_Z()", {
  expect_identical(convert_r_Z, convert_r_Z)
  for (r in c(-0.9, -0.5, 0, 0.5, 0.9)) {
    expect_equal(convert_r_Z(r)$value, convert_r_Z(r)$value, tolerance = 1e-12)
  }
})

test_that("convert_Z_r() alias is identical to convert_Z_r()", {
  expect_identical(convert_Z_r, convert_Z_r)
  for (Z in c(-2, -0.5, 0, 0.5, 2)) {
    expect_equal(convert_Z_r(Z)$value, convert_Z_r(Z)$value, tolerance = 1e-12)
  }
})

test_that("both converts insist on a single value and point at the vector idiom", {
  expect_error(convert_r_Z(c(0.1, 0.2)), "single numeric value")
  expect_error(convert_r_Z(NA_real_), "single numeric value")
  expect_error(convert_Z_r(c(0.1, 0.2)), "single numeric value")
  expect_error(convert_Z_r(numeric(0)), "single numeric value")
})
