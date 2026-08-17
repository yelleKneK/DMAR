test_that("average_variance_extracted() computes the Fornell-Larcker mean", {
  res <- average_variance_extracted(loadings = c(.8, .7, .6))
  expect_s3_class(res, "dmar_tbl")
  expect_equal(res$ave, mean(c(.8, .7, .6)^2))
})

test_that("average_variance_extracted() on a fit is the mean squared standardized loading", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  f <- rnorm(250)
  d <- data.frame(y1 = .8 * f + rnorm(250, 0, .6),
                  y2 = .7 * f + rnorm(250, 0, .7),
                  y3 = .6 * f + rnorm(250, 0, .8))
  fit <- lavaan::cfa("f =~ y1 + y2 + y3", data = d, std.lv = TRUE)
  res <- average_variance_extracted(fit)
  expect_equal(res$factor, "f")

  # The oracle is the Fornell and Larcker (1981) definition computed from
  # scratch off the fitted solution: the mean of the squared completely
  # standardized loadings. Deriving it here rather than calling another
  # package's AVE() keeps the assertion a test of the definition, and it
  # holds to machine precision rather than to a third decimal.
  std <- lavaan::standardizedSolution(fit)
  lam <- std$est.std[std$op == "=~" & std$lhs == "f"]
  expect_equal(res$ave, mean(lam^2))
})

test_that("average_variance_extracted() validates", {
  expect_error(average_variance_extracted(), "exactly one")
  expect_error(average_variance_extracted(loadings = .8), "two or more")
  expect_error(average_variance_extracted(fit = lm(mpg ~ wt, mtcars)),
               "lavaan fit")
})

test_that("the percentile bootstrap brackets the estimate and reproduces", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # 2 x 300 CFA refits; the cheap-path tests cover CRAN
  set.seed(113)
  f <- rnorm(250)
  d <- data.frame(y1 = .8 * f + rnorm(250, 0, .6),
                  y2 = .7 * f + rnorm(250, 0, .7),
                  y3 = .6 * f + rnorm(250, 0, .8))
  fit <- lavaan::cfa("f =~ y1 + y2 + y3", data = d, std.lv = TRUE)

  plain <- average_variance_extracted(fit)
  expect_true(all(is.na(plain$ci_lower)) && all(is.na(plain$ci_upper)))

  b1 <- average_variance_extracted(fit, ci_method = "percentile",
                                   B = 300, seed = 113)
  b2 <- average_variance_extracted(fit, ci_method = "percentile",
                                   B = 300, seed = 113)
  expect_identical(as.data.frame(b1), as.data.frame(b2))
  expect_equal(b1$ave, plain$ave)
  expect_lt(b1$ci_lower, b1$ave)
  expect_gt(b1$ci_upper, b1$ave)
  expect_identical(attr(b1, "B_used"), 300L)
  expect_identical(attr(b1, "conf_level"), 0.95)
})

test_that("loadings alone cannot carry an interval", {
  expect_error(
    average_variance_extracted(loadings = c(.8, .7, .6),
                               ci_method = "percentile"),
    "loadings' alone")
  # The loadings path keeps the canonical four columns.
  r <- average_variance_extracted(loadings = c(.8, .7, .6))
  expect_named(r, c("factor", "ave", "ci_lower", "ci_upper"))
})
