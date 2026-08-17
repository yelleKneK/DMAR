test_that("equivalence_r() returns documented rows", {
  res <- equivalence_r(r = 0.05, n = 200, rho_upper = 0.10)
  expect_setequal(res$term,
                  c("r", "z_lower_test", "z_upper_test",
                    "p_lower", "p_upper", "p_tost",
                    "lower_limit", "upper_limit",
                    "rho_lower", "rho_upper",
                    "equivalent", "n"))
})

test_that("equivalence_r() declares equivalence for r near 0 with large n", {
  res <- equivalence_r(r = 0.01, n = 1000, rho_upper = 0.10)
  expect_equal(res$value[res$term == "equivalent"], 1)
})

test_that("equivalence_r() does not declare equivalence for r near the bound", {
  res <- equivalence_r(r = 0.18, n = 100, rho_upper = 0.20)
  expect_equal(res$value[res$term == "equivalent"], 0)
})

test_that("equivalence_r() raw-data interface matches summary interface", {
  set.seed(113)
  x <- rnorm(150); y <- 0.04 * x + rnorm(150)
  r_obs <- cor(x, y)
  r1 <- equivalence_r(x = x, y = y, rho_upper = 0.15)
  r2 <- equivalence_r(r = r_obs, n = 150, rho_upper = 0.15)
  expect_equal(r1$value[r1$term == "p_tost"],
               r2$value[r2$term == "p_tost"], tolerance = 1e-10)
})
