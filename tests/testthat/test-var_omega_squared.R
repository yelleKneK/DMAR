test_that("var_omega_squared() returns a non-negative variance", {
  res <- var_omega_squared(0.10, 2, 57, 60)
  expect_named(res, c("term", "value"))
  expect_equal(res$term, "var_omega_squared")
  expect_gt(res$value, 0)
})

test_that("var_omega_squared() variance grows with smaller df_error", {
  high_df <- var_omega_squared(0.10, 2, 200, 203)$value
  low_df  <- var_omega_squared(0.10, 2, 20,   23)$value
  expect_gt(low_df, high_df)
})

test_that("var_omega_squared() rejects bad inputs", {
  expect_error(var_omega_squared(1.2, 2, 57, 60), "in \\[0, 1\\)")
  expect_error(var_omega_squared(0.1, 0, 57, 60), ">=")
  expect_error(var_omega_squared(0.1, 2, 3,  60), ">=")
})

test_that("the variance tracks a Monte Carlo of the sampling distribution", {
  # The sampling distribution of omega^2-hat is a function of the noncentral
  # F, so it can be simulated directly from chi-square draws without
  # generating data. Fleishman (1980, Eq. 22) transferred to the omega^2
  # scale should match that variance.
  mc_var <- function(omega2, J, n, B = 200000L, seed = 113L) {
    set.seed(seed)
    N <- J * n; df_e <- J - 1; df_r <- N - J
    lambda <- omega2 * N / (1 - omega2)
    F_v <- (rchisq(B, df_e, ncp = lambda) / df_e) / (rchisq(B, df_r) / df_r)
    o2 <- df_e * (F_v - 1) / (df_e * (F_v - 1) + N)
    var(o2)
  }
  for (om2 in c(0.06, 0.25)) {
    for (J in c(3, 5)) {
      for (n in c(30, 100)) {
        N <- J * n
        got <- var_omega_squared(om2, J - 1, N - J, N)$value
        expect_equal(got / mc_var(om2, J, n), 1, tolerance = 0.10)
      }
    }
  }
})

test_that("the variance is consistent: N times the variance stays bounded", {
  # The regression guard for the defect this replaced, whose N * var grew
  # without bound (0.8 at N = 30 to 51 at N = 5100) because it converged to
  # a positive constant instead of shrinking like 1/N.
  Nv <- vapply(c(30, 120, 300, 1200, 5100), function(N) {
    N * var_omega_squared(0.10, 2, N - 3, N)$value
  }, numeric(1))
  expect_true(all(diff(Nv) < 0))          # decreasing, not growing
  expect_true(all(Nv < 1))                # and bounded well below 1
  expect_equal(Nv[5] / Nv[3], 1, tolerance = 0.15)   # settling to a limit
})
