test_that("signal_to_noise_R2() returns a tidy data frame with the expected estimators", {
  result <- signal_to_noise_R2(R2 = 0.3, N = 100, p = 5)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_setequal(
    result$term,
    c("phi2_hat", "phi2_adj_hat", "phi2_umvue", "phi2_umvue_l", "phi2_umvue_nl")
  )
})

test_that("signal_to_noise_R2() returns a dmar_tbl", {
  expect_s3_class(signal_to_noise_R2(R2 = .5, N = 50, p = 2), "dmar_tbl")
})

test_that("signal_to_noise_R2() phi2_hat equals R2 / (1 - R2)", {
  R2 <- 0.3
  result <- signal_to_noise_R2(R2 = R2, N = 100, p = 5)
  expect_equal(
    result$value[result$term == "phi2_hat"],
    R2 / (1 - R2),
    tolerance = 1e-12
  )
})

test_that("signal_to_noise_R2() adjusted estimator is smaller than the naive plug-in", {
  result <- signal_to_noise_R2(R2 = 0.3, N = 100, p = 5)
  expect_lt(
    result$value[result$term == "phi2_adj_hat"],
    result$value[result$term == "phi2_hat"]
  )
})

test_that("phi2_umvue_nl adds Muirhead's c/Y term to the untruncated linear estimate", {
  # Muirhead (1985, Eq. 10) defines theta_NL from the untruncated theta_L.
  # With sample R2 below p/(N-1) the raw theta_U is negative, and building
  # from a truncated theta_L would overstate theta_NL.
  out <- signal_to_noise_R2(R2 = 0.05, N = 40, p = 6)
  v <- setNames(out$value, out$term)
  n <- 39; m <- 7; Y <- 0.05 / 0.95
  theta_U_raw <- ((n - m - 1) / n) * Y - (m - 1) / n
  expect_lt(theta_U_raw, 0)
  shrink <- (n * (n - m - 3)) / ((n + 2) * (n - m - 1))
  c_term <- (2 * (n - 2) * (n - m - 3) * (m - 5)) /
    ((n + 2) * (n - m - 1) * (n - m + 1) * (n - m + 3))
  expect_equal(unname(v[["phi2_umvue_nl"]]),
               max(shrink * theta_U_raw + c_term / Y, 0), tolerance = 1e-12)
})
