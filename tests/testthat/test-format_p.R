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

test_that("format_p() refuses a fractional digits_p rather than cutting it (DR46)", {
  # Through the first builds of 1.1.0 digits_p = 3.9 was cut to 3 by
  # as.integer() and showed "0.123".
  for (bad in list(3.9, 2.5, 1.0000001)) {
    err <- tryCatch(format_p(.123456, digits_p = bad), error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "digits_p")
    expect_match(conditionMessage(err), "single positive integer", fixed = TRUE)
  }
  # A whole number given as a double or an integer is accepted, and the
  # formatting is that of sprintf() at that many decimal places.
  expect_identical(format_p(.123456, digits_p = 3), "0.123")
  expect_identical(format_p(.123456, digits_p = 3L), "0.123")
  expect_identical(format_p(.123456, digits_p = 3),
                   sprintf("%.3f", .123456))
  expect_identical(format_p(4e-4, digits_p = 3), "< 0.001")
})

test_that("format_p() is display only: values outside [0, 1] are shown as given (DR51)", {
  # The function does not check that p is a probability; the help page
  # says what each kind of value shows as.
  expect_identical(format_p(c(-.1, 1.1, NA, NaN, -Inf)),
                   c("< 0.0001", "1.1000", NA_character_, NA_character_,
                     "< 0.0001"))
  expect_identical(trimws(format_p(Inf)), "Inf")
})

test_that("a refused digits_p is named as given, and the largest is 324 (review)", {
  # A value that seven significant digits round to a whole number was
  # named as that number ("1 was given" for 1 + 1e-12).
  err <- tryCatch(format_p(.5, digits_p = 1 + 1e-12), error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_match(conditionMessage(err), "1.000000000001 was given", fixed = TRUE)
  # 324 is the most decimal places formatC() writes; 1e10 became NA in
  # as.integer() with a warning.
  expect_identical(format_p(.5, digits_p = 324),
                   formatC(.5, format = "f", digits = 324))
  for (bad in list(325, 1e10)) {
    err <- tryCatch(format_p(.5, digits_p = bad),
                    error = function(e) e, warning = function(w) w)
    expect_s3_class(err, "dmar_argument_invalid")
  }
})
