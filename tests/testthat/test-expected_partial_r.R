test_that("expected_partial_r() returns documented columns and is downward-biased", {
  res <- expected_partial_r(rho = 0.4, n = 30, J = 3)
  expect_named(res, c("rho", "n", "J", "expected_partial_r",
                      "bias", "relative_bias"))
  expect_lt(res$expected_partial_r, res$rho)  # downward bias
  expect_gt(res$bias, 0)
})

test_that("expected_partial_r() returns a dmar_tbl", {
  expect_s3_class(expected_partial_r(rho = 0.4, n = 30, J = 3), "dmar_tbl")
})

test_that("expected_partial_r() reduces to expected_r() at J = 1 with adjusted n", {
  # partial-r at (n, J = 1) is the same sampling distribution as
  # simple-r at n - 1 (Anderson, 2003).
  pr <- expected_partial_r(rho = 0.5, n = 11, J = 1)$expected_partial_r
  sr <- expected_r(rho = 0.5, n = 10)$expected_r
  expect_equal(pr, sr, tolerance = 1e-12)
})

test_that("expected_partial_r() bias grows with k at fixed (rho, n)", {
  res <- expected_partial_r(rho = 0.4, n = 30, J = c(1, 2, 5, 10, 15))
  expect_true(all(diff(res$bias) > 0))  # bias monotone increasing in k
})

test_that("expected_partial_r() rejects bad inputs", {
  expect_error(expected_partial_r(rho = 1.5, n = 30, J = 2),  "lie in")
  expect_error(expected_partial_r(rho = 0.4, n = 5,  J = 5),  ">=")
  expect_error(expected_partial_r(rho = 0.4, n = 30, J = 0),  ">=")
})

# Doc-code agreement: the help page's displayed expectation is
# E[r | rho, n, J] = rho * 2F1(1/2, 1/2; (n - J + 1)/2; rho^2)
#   * Gamma((n - J)/2)^2 / (Gamma((n - J - 1)/2) Gamma((n - J + 1)/2)).
# The formula previously printed carried every index one lower and, as
# printed, exceeds 1 (1.0498 at rho = 0.99, n - J = 5), impossible for
# the expectation of a correlation. This test computes the corrected
# expression independently of the package's series code.
test_that("expected_partial_r() matches the documented formula under the n to n - J substitution", {
  hyp2f1_half <- function(c_par, z) {
    s <- 1; term <- 1
    for (k in 1:20000) {
      term <- term * ((k - 0.5)^2) / ((c_par + k - 1) * k) * z
      s <- s + term
      if (abs(term) < 1e-15 * abs(s)) break
    }
    s
  }
  doc_formula <- function(rho, n, J) {
    m <- n - J
    rho * hyp2f1_half((m + 1) / 2, rho^2) *
      exp(2 * lgamma(m / 2) - lgamma((m - 1) / 2) - lgamma((m + 1) / 2))
  }
  grid <- expand.grid(rho = c(-0.7, 0.2, 0.4, 0.9),
                      n   = c(15, 30, 100),
                      J   = c(1, 3, 8))
  fn  <- expected_partial_r(grid$rho, grid$n, grid$J)$expected_partial_r
  doc <- mapply(doc_formula, grid$rho, grid$n, grid$J)
  expect_equal(fn, doc, tolerance = 1e-12)
  expect_true(all(abs(fn) < 1))
})

test_that("expected_partial_r() bias magnitudes match the values quoted in the help page", {
  # The page quotes +0.00624 (rho = 0.4, n = 30, J = 2) and +0.00888
  # (J = 10); the exact values are 0.00624153 and 0.00887618.
  expect_equal(round(expected_partial_r(0.4, 30, 2)$bias, 5), 0.00624)
  expect_equal(round(expected_partial_r(0.4, 30, 10)$bias, 5), 0.00888)
})
