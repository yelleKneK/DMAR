test_that("var_r() returns the Fisher (1915) and Fisher's Z variances", {
  res <- var_r(true_cor = 0.30, N = 50)
  expect_setequal(res$term, c("var_r_normal", "var_fisher_z"))
  expect_equal(res$value[res$term == "var_r_normal"],
               (1 - 0.30^2)^2 / 49, tolerance = 1e-12)
  expect_equal(res$value[res$term == "var_fisher_z"],
               1 / 47, tolerance = 1e-12)
})

test_that("var_r() no longer takes the marginal kurtoses", {
  expect_error(var_r(true_cor = 0.30, N = 50, kurtosis_x = 3, kurtosis_y = 3),
               "unused argument")
})

test_that("var_r() rejects bad inputs", {
  expect_error(var_r(1.5, 50),  "\\(-1, 1\\)")
  expect_error(var_r(0.3, 3),   ">= 4")
})

test_that("var_r() help statements: delta method transfer, Z variance, sign", {
  # The help page says the large-sample variance of r, carried to the
  # Fisher Z scale by the delta method (dZ/dr = 1 / (1 - rho^2)), is
  # 1 / (N - 1), that ci_r() uses 1 / (N - 3) on that scale, and that the
  # variance of Z depends on N alone. Oracles: the algebra of the delta
  # method and Fisher's (1921) variance, written out here.
  for (rho in c(-0.8, -0.3, 0, 0.45, 0.9)) {
    for (N in c(10, 50, 400)) {
      res <- var_r(true_cor = rho, N = N)
      v_r <- res$value[res$term == "var_r_normal"]
      v_z <- res$value[res$term == "var_fisher_z"]
      expect_equal(v_r / (1 - rho^2)^2, 1 / (N - 1), tolerance = 1e-12)
      expect_equal(v_z, 1 / (N - 3), tolerance = 1e-12)
      # Reversing the sign of the correlation leaves both variances alone.
      expect_equal(var_r(true_cor = -rho, N = N)$value, res$value,
                   tolerance = 1e-14)
    }
  }
})
