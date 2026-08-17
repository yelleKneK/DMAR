test_that("var_r() returns the Fisher (1915) and Fisher's Z variances", {
  res <- var_r(rho = 0.30, n = 50)
  expect_setequal(res$term, c("var_r_normal", "var_fisher_z"))
  expect_equal(res$value[res$term == "var_r_normal"],
               (1 - 0.30^2)^2 / 49, tolerance = 1e-12)
  expect_equal(res$value[res$term == "var_fisher_z"],
               1 / 47, tolerance = 1e-12)
})

test_that("var_r() returns the Bonett-Wright row when kurtoses are supplied", {
  res <- var_r(rho = 0.30, n = 50, kurtosis_x = 3, kurtosis_y = 3)
  expect_true("var_r_bonett_wright" %in% res$term)
  # Bonett-Wright variance >= normal-theory variance for positive kurtosis:
  expect_gt(res$value[res$term == "var_r_bonett_wright"],
            res$value[res$term == "var_r_normal"])
})

test_that("var_r() rejects bad inputs", {
  expect_error(var_r(1.5, 50),  "\\(-1, 1\\)")
  expect_error(var_r(0.3, 3),   ">= 4")
})
