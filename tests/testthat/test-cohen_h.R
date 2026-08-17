test_that("cohen_h is the arcsine difference and matches known values", {
  phi <- function(p) 2 * asin(sqrt(p))
  expect_equal(cohen_h(0.55, 0.40)$value, phi(0.55) - phi(0.40), tolerance = 1e-12)
  # Cohen (1988) example: p1 = .5, p2 = .5 gives h = 0.
  expect_equal(cohen_h(0.5, 0.5)$value, 0)
  # signed: p1 > p2 is positive; swapping flips the sign.
  expect_gt(cohen_h(0.6, 0.4)$value, 0)
  expect_equal(cohen_h(0.6, 0.4)$value, -cohen_h(0.4, 0.6)$value, tolerance = 1e-12)
})

test_that("the same raw difference is a larger h nearer the floor", {
  near_floor  <- abs(cohen_h(0.20, 0.05)$value)
  near_middle <- abs(cohen_h(0.55, 0.40)$value)   # both raw differences are .15
  expect_gt(near_floor, near_middle)
})

test_that("cohen_h returns the tidy term/value shape", {
  out <- cohen_h(0.3, 0.5)
  expect_s3_class(out, "data.frame")
  expect_equal(names(out), c("term", "value"))
  expect_equal(out$term, "cohen_h")
})

test_that("cohen_h validates its inputs", {
  expect_error(cohen_h(1.2, 0.5), "in \\[0, 1\\]")
  expect_error(cohen_h(0.5, -0.1), "in \\[0, 1\\]")
  expect_error(cohen_h(c(0.3, 0.4), 0.5), "single proportion")
})
