test_that("ss_aipe_cv() returns a 1-row tidy data frame with necessary_N", {
  result <- ss_aipe_cv(C_of_V = 0.2, width = 0.1, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, "necessary_N")
})

test_that("ss_aipe_cv() requires a larger n for a narrower width", {
  # width = 0.05 forces n large enough that the noncentral t's ncp exceeds R's
  # 37.62 stable limit and the function emits the documented precision warning
  # at that boundary; suppress only the expected warning.
  n_wide   <- ss_aipe_cv(C_of_V = 0.2, width = 0.20, conf_level = 0.95)$value
  n_narrow <- suppressWarnings(ss_aipe_cv(C_of_V = 0.2, width = 0.05, conf_level = 0.95))$value
  expect_lt(n_wide, n_narrow)
})

test_that("ss_aipe_cv() errors on non-positive coefficient of variation", {
  expect_error(ss_aipe_cv(C_of_V = 0,  width = 0.1, conf_level = 0.95), "positive")
  expect_error(ss_aipe_cv(C_of_V = -1, width = 0.1, conf_level = 0.95), "positive")
})

test_that("ss_aipe_cv() errors on a non-positive width instead of never terminating", {
  # A non-positive width makes the search target unreachable, so the
  # increment-by-one search previously ran without terminating. It must error
  # at entry rather than hang.
  expect_error(ss_aipe_cv(C_of_V = 0.1, width = -1, conf_level = 0.99), "width")
  expect_error(ss_aipe_cv(C_of_V = 0.1, width = 0, conf_level = 0.99), "width")
})

test_that("ss_aipe_cv() validates conf_level", {
  expect_error(ss_aipe_cv(C_of_V = 0.1, width = 0.1, conf_level = 1.2),
               "conf_level")
})

test_that("ss_aipe_cv() accepts the mu/sigma parameterization", {
  # Regression: neither formal was read, so the documented mu/sigma call
  # died on "argument is of length zero". mu = 10 and sigma = 1 imply
  # C_of_V = .1; both parameterizations of the same population must plan
  # the same sample size, anchored at the N = 20 the C_of_V path has
  # always returned for this input. suppressWarnings() covers only the
  # documented noncentrality-precision warning the .99 search emits.
  n_musigma <- suppressWarnings(
    ss_aipe_cv(mu = 10, sigma = 1, width = .1, conf_level = .99)
  )
  n_cv <- suppressWarnings(
    ss_aipe_cv(C_of_V = .1, width = .1, conf_level = .99)
  )
  expect_named(n_musigma, c("term", "value"))
  expect_equal(n_musigma$term, "necessary_N")
  expect_equal(n_musigma$value, n_cv$value)
  expect_equal(n_musigma$value, 20)
})

test_that("ss_aipe_cv() rejects conflicting or incomplete parameterizations", {
  expect_error(
    ss_aipe_cv(C_of_V = .1, mu = 10, sigma = 1, width = .1),
    "either 'C_of_V' or both 'mu' and 'sigma'"
  )
  expect_error(ss_aipe_cv(mu = 10, width = .1),
               "both 'mu' and 'sigma'")
  expect_error(ss_aipe_cv(sigma = 1, width = .1),
               "both 'mu' and 'sigma'")
})
