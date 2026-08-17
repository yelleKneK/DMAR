test_that("format_p() rounds to 4 dp by default", {
  expect_equal(format_p(0.5), "0.5000")
  expect_equal(format_p(0.0234), "0.0234")
  expect_equal(format_p(0.0001234), "0.0001")
  expect_equal(format_p(0.05), "0.0500")
})

test_that("format_p() uses '< 0.0001' for values below the default floor", {
  expect_equal(format_p(1e-10), "< 0.0001")
  expect_equal(format_p(0.00001), "< 0.0001")
  expect_equal(format_p(2.2e-16), "< 0.0001")
})

test_that("format_p() returns NA_character_ for NA input, preserving length", {
  expect_identical(format_p(NA_real_), NA_character_)
  expect_identical(format_p(c(0.05, NA, 0.5)),
                   c("0.0500", NA_character_, "0.5000"))
})

test_that("format_p() respects digits_p and tracks the floor with it", {
  expect_equal(format_p(0.0001234, digits_p = 6), "0.000123")
  expect_equal(format_p(1e-10,     digits_p = 6), "< 0.000001")
  expect_equal(format_p(0.5,       digits_p = 2), "0.50")
  expect_equal(format_p(0.03,      digits_p = 2), "0.03")
  expect_equal(format_p(0.005,     digits_p = 2), "< 0.01")
})

test_that("format_p() vectorizes", {
  out <- format_p(c(0.5, 0.0234, 1e-10, NA))
  expect_length(out, 4L)
  expect_type(out, "character")
  expect_equal(out, c("0.5000", "0.0234", "< 0.0001", NA_character_))
})

test_that("format_p() errors on bad input", {
  expect_error(format_p("0.05"), "must be numeric")
  expect_error(format_p(0.05, digits_p = 0), "positive integer")
  expect_error(format_p(0.05, digits_p = c(2, 3)), "single positive integer")
  expect_error(format_p(0.05, digits_p = NA_integer_), "single positive integer")
})

test_that("format_p() never returns scientific-notation strings", {
  out <- format_p(c(2.2e-16, 1e-300, 0.001, 0.5))
  expect_false(any(grepl("e[+-]?[0-9]", out)))
})
