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
      error_variance = 150.5, var_Z = 210.7, n_1 = 64, n_2 = 246,
      reg_coef_1 = 0.968, reg_coef_2 = 0.778, mean_Z = 100, fixed_value = 115,
      type = grid$type[i], covariate_value = grid$covariate_value[i])
    expect_equal(ours$value, ref[i], tolerance = 1e-12,
                 info = paste(grid$covariate_value[i], grid$type[i]))
  }
})

test_that("var_ete() with covariate_value = 'sd' is right at a very large N", {
  # c_1 = 1/N + 1 - c4(N)^2 there, and 1 - c4(N)^2 = Var(S) / sigma^2 comes
  # from .chi_var_factor() (R/c4_internals.R). Through 1.0.x it was formed
  # from lgamma() values and the subtraction lost the leading digits: 3e-6
  # of relative error in the variance at n_1 = n_2 = 5e5, and c_1 itself
  # negative at N = 1e8. The pins evaluate the same expression in Rmpfr at
  # 320 bits, with every input taken as the double it is.
  sizes <- c(5e5, 5e7)
  pins <- list(c(sample = 0.001215411514865371,
                 population = 0.0012154124178692088),
               c(sample = 1.2154094260985249e-05,
                 population = 1.2154094351285252e-05))
  for (i in seq_along(sizes)) {
    for (ty in c("sample", "population")) {
      ours <- var_ete(
        error_variance = 150.5, var_Z = 210.7, n_1 = sizes[i], n_2 = sizes[i],
        reg_coef_1 = 0.968, reg_coef_2 = 0.778, mean_Z = 100, fixed_value = 115,
        type = ty, covariate_value = "sd")
      expect_equal(ours$value, pins[[i]][[ty]], tolerance = 1e-12,
                   info = paste(sizes[i], "per group,", ty))
    }
  }
})

test_that("var_ete() returns the house shape and validates input", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  res <- var_ete(error_variance = 100, var_Z = 200, n_1 = 50, n_2 = 60,
                 reg_coef_1 = 1, reg_coef_2 = 0.5)
  expect_s3_class(res, "dmar_tbl")
  expect_identical(res$term, "var_ete")
  expect_true(is.numeric(res$value) && res$value > 0)
  expect_identical(attr(res, "type"), "sample")
  expect_identical(attr(res, "covariate_value"), "sample_mean")

  expect_error(var_ete(error_variance = -1, var_Z = 1, n_1 = 10, n_2 = 10,
                       reg_coef_1 = 1, reg_coef_2 = 1), "positive")
  expect_error(var_ete(error_variance = 1, var_Z = 1, n_1 = 3, n_2 = 10,
                       reg_coef_1 = 1, reg_coef_2 = 1), "greater than 3")
})
