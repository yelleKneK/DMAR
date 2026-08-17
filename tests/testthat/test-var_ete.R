# var_ete(): variance of the estimated treatment effect at selected
# covariate values (Li, McLouth, & Delaney, 2020). The numerical anchor
# is MBESS::var.ete(), the reference implementation contributed by
# Li Li, checked across every branch (three covariate values by two
# variance types).

test_that("var_ete() matches MBESS::var.ete() on every branch", {
  grid <- expand.grid(
    covariate_value = c("sample_mean", "sd", "fixed"),
    type = c("sample", "population"),
    stringsAsFactors = FALSE
  )
  # Pinned from MBESS::var.ete (MBESS 4.9.3, 2026-08-09), one value per grid
  # row in order; live comparison in tools/oracle_checks.R.
  ref <- c(3.009106514479761, 6.083023321637939, 6.269164859483777,
           3.018998698280908, 6.097873582984068, 6.300455054358817)

  for (i in seq_len(nrow(grid))) {
    ours <- var_ete(
      sigma2 = 150.5, sigma2_Z = 210.7, n_1 = 64, n_2 = 246,
      beta_1 = 0.968, beta_2 = 0.778, mu_Z = 100, fixed_value = 115,
      type = grid$type[i], covariate_value = grid$covariate_value[i])
    expect_equal(ours$value, ref[i], tolerance = 1e-12,
                 info = paste(grid$covariate_value[i], grid$type[i]))
  }
})

test_that("var_ete() returns the house shape and validates input", {
  res <- var_ete(sigma2 = 100, sigma2_Z = 200, n_1 = 50, n_2 = 60,
                 beta_1 = 1, beta_2 = 0.5)
  expect_s3_class(res, "dmar_tbl")
  expect_identical(res$term, "var_ete")
  expect_true(is.numeric(res$value) && res$value > 0)
  expect_identical(attr(res, "type"), "sample")
  expect_identical(attr(res, "covariate_value"), "sample_mean")

  expect_error(var_ete(sigma2 = -1, sigma2_Z = 1, n_1 = 10, n_2 = 10,
                       beta_1 = 1, beta_2 = 1), "positive")
  expect_error(var_ete(sigma2 = 1, sigma2_Z = 1, n_1 = 3, n_2 = 10,
                       beta_1 = 1, beta_2 = 1), "greater than 3")
})
