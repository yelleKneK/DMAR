test_that("ss_aipe_c() returns a 1-row tidy data frame with necessary_n_per_group", {
  result <- ss_aipe_c(error_variance = 4, c_weights = c(1, -.5, -.5), width = 0.5)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, "necessary_n_per_group")
})

test_that("ss_aipe_c() requires a larger n for larger error variance", {
  n_small <- ss_aipe_c(error_variance = 1, c_weights = c(1, -.5, -.5), width = 0.5)$value
  n_large <- ss_aipe_c(error_variance = 9, c_weights = c(1, -.5, -.5), width = 0.5)$value
  expect_lt(n_small, n_large)
})

test_that("ss_aipe_c() floors an impossibly wide target at the admissible minimum", {
  # An enormous target width was previously met at a degenerate n = 1, which
  # cannot estimate the within-group error variance (df = J * (n - 1) = 0). The
  # smallest admissible design is n = 2 per group.
  r <- ss_aipe_c(error_variance = 40, c_weights = c(1, -.5, -.5), width = 100)
  expect_equal(r$value[r$term == "necessary_n_per_group"], 2)
})

test_that("ss_aipe_c() validates its boundary inputs", {
  expect_error(ss_aipe_c(error_variance = 40, c_weights = c(1, -.5, -.5),
                         width = 3, conf_level = 1.2), "conf_level")
  expect_error(ss_aipe_c(error_variance = 40, c_weights = c(1, -.5, -.5),
                         width = 0), "width")
  expect_error(ss_aipe_c(error_variance = 40, c_weights = c(1, -.5, -.5),
                         width = Inf), "width")
  expect_error(ss_aipe_c(error_variance = -1, c_weights = c(1, -.5, -.5),
                         width = 3), "error variance")
  expect_error(ss_aipe_c(error_variance = 40, c_weights = c(1, -.5, -.5),
                         width = 3, assurance = 1.5), "assurance")
})
