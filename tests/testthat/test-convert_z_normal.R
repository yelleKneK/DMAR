test_that("convert_z_normal() returns a tidy 1-row data.frame", {
  res <- convert_z_normal(z = 1.96, mean = 100, sd = 15)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("term", "value"))
  expect_equal(res$term, "value_from_z")
})

test_that("convert_z_normal() matches value = mean + z * sd", {
  for (z in c(-2, -1, -0.5, 0, 0.5, 1, 2)) {
    expect_equal(
      convert_z_normal(z, mean = 100, sd = 15)$value,
      100 + z * 15,
      tolerance = 1e-12
    )
  }
})

test_that("convert_z_normal() matches qnorm(pnorm(z), mean, sd)", {
  for (z in c(-2, -1, -0.5, 0, 0.5, 1, 2)) {
    expect_equal(
      convert_z_normal(z, mean = 100, sd = 15)$value,
      qnorm(pnorm(z), mean = 100, sd = 15),
      tolerance = 1e-9
    )
  }
})

test_that("convert_z_normal() returns z unchanged under the standard normal defaults", {
  for (z in c(-2, -1, 0, 1, 1.96, 2)) {
    expect_equal(convert_z_normal(z)$value, z, tolerance = 1e-12)
  }
})

test_that("convert_z_normal() preserves the percentile of z", {
  for (z in c(-1.5, 0, 0.8, 1.96)) {
    value <- convert_z_normal(z, mean = 50, sd = 10)$value
    expect_equal(pnorm(value, mean = 50, sd = 10), pnorm(z), tolerance = 1e-12)
  }
})

test_that("convert_z_normal() rejects nonnumeric input and negative sd", {
  expect_error(convert_z_normal("a"), "'z' must be a single number")
  expect_error(convert_z_normal(1, mean = "a"), "'mean' must be a single number")
  expect_error(convert_z_normal(1, sd = "a"), "'sd' must be a single number")
  expect_error(convert_z_normal(1, sd = -1), "'sd' must be nonnegative")
})
