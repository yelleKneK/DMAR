test_that("expected_partial_r() returns documented columns and is downward-biased", {
  res <- expected_partial_r(true_cor = 0.4, N = 30, p = 4)
  expect_named(res, c("true_cor", "N", "p", "expected_partial_r",
                      "bias", "relative_bias"))
  expect_identical(res$N, 30)
  expect_identical(res$p, 4)
  expect_lt(res$expected_partial_r, res$true_cor)  # downward bias
  expect_gt(res$bias, 0)
})

test_that("expected_partial_r() returns a dmar_tbl", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  expect_s3_class(expected_partial_r(true_cor = 0.4, N = 30, p = 4), "dmar_tbl")
})

test_that("expected_partial_r() reduces to expected_r() at p = 2 with adjusted N", {
  # partial-r at (N, p = 2), one variable partialed out, has the same
  # sampling distribution as simple-r at N - 1 (Anderson, 2003).
  pr <- expected_partial_r(true_cor = 0.5, N = 11, p = 2)$expected_partial_r
  sr <- expected_r(true_cor = 0.5, N = 10)$expected_r
  expect_equal(pr, sr, tolerance = 1e-12)
})

test_that("expected_partial_r() bias grows with p at fixed (rho, N)", {
  res <- expected_partial_r(true_cor = 0.4, N = 30, p = c(2, 3, 6, 11, 16))
  expect_true(all(diff(res$bias) > 0))  # bias monotone increasing in p
})

test_that("expected_partial_r() rejects bad inputs", {
  expect_error(expected_partial_r(true_cor = 1.5, N = 30, p = 3),  "lie in")
  expect_error(expected_partial_r(true_cor = 0.4, N = 5,  p = 6),
               "'N - p' must be at least 3", fixed = TRUE)
  expect_error(expected_partial_r(true_cor = 0.4, N = 30, p = 1),
               "whole number of at least 2")
  expect_error(expected_partial_r(true_cor = 0.4, N = 30, p = 2.5),
               "whole number of at least 2")
  expect_error(expected_partial_r(true_cor = 0.4, N = 30, p = c(2, 3.5, 4)),
               "whole number of at least 2")
  expect_error(expected_partial_r(true_cor = 0.4, N = 30, p = 0),
               "whole number of at least 2")
  expect_error(expected_partial_r(true_cor = 0.4, N = 30),
               "'p' (the number of predictors) is required.", fixed = TRUE)
  expect_error(expected_partial_r(true_cor = 0.4, N = "30", p = 3),
               "'N' must be numeric.", fixed = TRUE)
  expect_error(expected_partial_r(true_cor = 0.4, N = 30, p = "3"),
               "'p' must be numeric.", fixed = TRUE)
  # The boundary the message states: N - p = 2 stops and N - p = 3 runs.
  expect_error(expected_partial_r(true_cor = 0.4, N = 6, p = 4), "at least 3")
  expect_no_error(expected_partial_r(true_cor = 0.4, N = 7, p = 4))
})

test_that("expected_partial_r() takes p in its place and does not accept the DMAR 1.0.x names", {
  new <- expected_partial_r(true_cor = 0.4, N = 30, p = c(2, 4))
  expect_identical(expected_partial_r(0.4, 30, c(2, 4)), new)
  expect_identical(names(formals(expected_partial_r)), c("true_cor", "N", "p"))
  expect_error(expected_partial_r(true_cor = 0.4, N = 30, J = 1), "unused argument")
  expect_error(expected_partial_r(true_cor = 0.4, n = 30, p = 2), "unused argument")
  expect_error(expected_partial_r(true_cor = 0.4, N = 30),
               "'p' (the number of predictors) is required.", fixed = TRUE)
})

# Doc-code agreement: the help page's displayed expectation is, with
# m = N - (p - 1) and J = p - 1 the number of variables partialed out,
# E[r | rho, N, p] = rho * 2F1(1/2, 1/2; (m + 1)/2; rho^2)
#   * Gamma(m/2)^2 / (Gamma((m - 1)/2) Gamma((m + 1)/2)).
# The formula previously printed carried every index one lower and, as
# printed, exceeds 1 (1.0498 at rho = 0.99, m = 5), impossible for
# the expectation of a correlation. This test computes the corrected
# expression independently of the package's series code.
test_that("expected_partial_r() matches the documented formula under the N to N - (p - 1) substitution", {
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
  fn  <- expected_partial_r(grid$rho, grid$n, p = grid$J + 1)$expected_partial_r
  doc <- mapply(doc_formula, grid$rho, grid$n, grid$J)
  expect_equal(fn, doc, tolerance = 1e-12)
  expect_true(all(abs(fn) < 1))
})

test_that("expected_partial_r() bias magnitudes match the values quoted in the help page", {
  # The page quotes +0.00624 (rho = 0.4, N = 30, p = 3) and +0.00888
  # (p = 11); the exact values are 0.00624153 and 0.00887618.
  expect_equal(round(expected_partial_r(0.4, 30, p = 3)$bias, 5), 0.00624)
  expect_equal(round(expected_partial_r(0.4, 30, p = 11)$bias, 5), 0.00888)
})

test_that("expected_partial_r() keeps a positive bias of the right size at a very large n", {
  # It calls expected_r() at N - (p - 1), whose gamma ratio lost its digits
  # at a large N through 1.0.x (the bias was 31 times its value at 1e8,
  # rho = .5, and negative at 3e8), and takes its bias and relative bias,
  # now formed directly rather than as rho - E[r]; see test-expected_r.R.
  rho <- 0.5
  for (m in c(1e8, 3e8)) {
    res <- expected_partial_r(true_cor = rho, N = m + 2, p = 3)
    simple <- expected_r(true_cor = rho, N = m)
    expect_gt(res$bias, 0)
    expect_identical(res$expected_partial_r, simple$expected_r)
    expect_identical(res$bias, simple$bias)
    expect_identical(res$relative_bias, simple$relative_bias)
    expect_equal(res$bias, rho * (1 - rho^2) / (2 * (m - 1)),
                 tolerance = 2e-6, info = paste("N - (p - 1) =", m))
  }
})

test_that("expected_partial_r() bias agrees with 320-bit values at N - (p - 1) = 1e8", {
  # The bias of the simple r at N = 1e8 and rho = .5, pinned in
  # test-expected_r.R from Rmpfr at 320 bits.
  res <- expected_partial_r(true_cor = 0.5, N = 1e8 + 3, p = 4)
  expect_equal(res$bias, 1.8750000246093751e-9, tolerance = 1e-14)
  expect_equal(res$relative_bias, 3.7500000492187503e-9, tolerance = 1e-14)
  expect_true(is.na(expected_partial_r(true_cor = 0, N = 20, p = 3)$relative_bias))
})
