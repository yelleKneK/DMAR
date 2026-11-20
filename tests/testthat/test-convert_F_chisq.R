# Tests for convert_F_chisq() and convert_chisq_F(). Both conversions are
# exact once the degrees of freedom are known, so the assertions are
# equalities against independently computed quantities.

test_that("the default (Inf denominator) is the scaling conversion", {
  # chi square = df_numerator * F, and back F = chi square / df.
  expect_equal(convert_F_chisq(2.75, df_numerator = 3)$value, 3 * 2.75)
  expect_equal(convert_chisq_F(8.25, df = 3)$value, 8.25 / 3)
  # Inf may also be passed explicitly; it is the documented default.
  expect_equal(convert_F_chisq(2.75, df_numerator = 3,
                               df_denominator = Inf)$value, 3 * 2.75)
  # Scaling round-trips exactly.
  expect_equal(convert_chisq_F(convert_F_chisq(4.1, df_numerator = 2)$value,
                               df = 2)$value, 4.1)
})

test_that("a finite denominator matches the upper-tail probability", {
  for (nu1 in c(1, 3, 8)) {
    for (nu2 in c(5, 20, 200)) {
      f <- 2.75
      x <- convert_F_chisq(f, df_numerator = nu1, df_denominator = nu2)$value
      expect_equal(pchisq(x, nu1, lower.tail = FALSE),
                   pf(f, nu1, nu2, lower.tail = FALSE))
    }
  }
})

test_that("the two directions are inverses at a finite denominator", {
  set.seed(113)
  for (i in seq_len(50)) {
    nu1 <- sample(1:12, 1)
    nu2 <- sample(3:400, 1)
    f   <- runif(1, 0.01, 40)
    x    <- convert_F_chisq(f, df_numerator = nu1, df_denominator = nu2)$value
    back <- convert_chisq_F(x, df = nu1, df_denominator = nu2)$value
    expect_equal(back, f)
  }
})

test_that("probability matching approaches the scaling value as nu_2 grows", {
  # The two conversions are the same family: the finite-denominator value
  # rises to nu_1 * F from below as the denominator degrees of freedom
  # grow, reaching it in the limit.
  x <- vapply(c(10, 50, 1000, 100000),
              function(nu2) convert_F_chisq(2.75, 3, nu2)$value, numeric(1))
  expect_true(all(diff(x) > 0))
  expect_true(all(x < 3 * 2.75))
  expect_equal(x[4], 3 * 2.75, tolerance = 1e-3)
})

test_that("with one numerator degree of freedom this is t and z, squared", {
  # F(1, nu) = t(nu)^2 and chi square(1) = z^2, so probability matching a
  # squared t gives the square of the matched normal deviate.
  tv <- 2.2; nu <- 30
  expect_equal(convert_F_chisq(tv^2, df_numerator = 1,
                               df_denominator = nu)$value,
               qnorm(pt(tv, nu))^2)
})

test_that("probability matching stays finite in the far upper tail", {
  # Referring the LOWER-tail probability to qchisq() loses the answer:
  # pf() reaches 1 in double precision by about F = 500 at these degrees
  # of freedom, so the lower-tail composition returns Inf. The upper-tail
  # p-value used here (no logarithms) stays accurate far past that.
  expect_true(is.infinite(qchisq(pf(500, 3, 20), 3)))

  for (f in c(500, 1e5, 1e10, 1e20)) {
    x <- convert_F_chisq(f, df_numerator = 3, df_denominator = 20)$value
    expect_true(is.finite(x))
    expect_gt(x, 0)
  }
  x <- vapply(c(500, 1e5, 1e10),
              function(f) convert_F_chisq(f, 3, 20)$value, numeric(1))
  expect_true(all(diff(x) > 0))
})

test_that("the return honors the package output contract", {
  r <- convert_F_chisq(2.75, df_numerator = 3, df_denominator = 50)
  expect_s3_class(r, "data.frame")
  expect_s3_class(r, "dmar_tbl")
  expect_named(r, c("term", "value"))
  expect_identical(r$term, "chi_square_from_F")
  expect_type(r$value, "double")

  r2 <- convert_chisq_F(7.710814, df = 3, df_denominator = 50)
  expect_identical(r2$term, "F_from_chi_square")
  expect_type(r2$value, "double")
})

test_that("inputs are validated", {
  expect_error(convert_F_chisq("a", 3), "single non-missing number")
  expect_error(convert_F_chisq(c(1, 2), 3), "single non-missing number")
  expect_error(convert_F_chisq(NA_real_, 3), "single non-missing number")
  expect_error(convert_F_chisq(-1, 3), "nonnegative")
  expect_error(convert_F_chisq(2.75, 0), "greater than zero")
  expect_error(convert_F_chisq(2.75, 3, 0), "greater than zero")
  # The numerator degrees of freedom cannot be infinite; the denominator
  # can, and that is the documented default.
  expect_error(convert_F_chisq(2.75, Inf), "must be finite")
  expect_error(convert_F_chisq(2.75, 3, Inf), NA)

  expect_error(convert_chisq_F(-1, 3), "nonnegative")
  expect_error(convert_chisq_F(7.7, 0), "greater than zero")
  expect_error(convert_chisq_F(7.7, Inf), "must be finite")
})

test_that("a zero statistic converts to zero", {
  expect_equal(convert_F_chisq(0, df_numerator = 3)$value, 0)
  expect_equal(convert_chisq_F(0, df = 3)$value, 0)
  expect_equal(convert_F_chisq(0, df_numerator = 3, df_denominator = 50)$value, 0)
})
