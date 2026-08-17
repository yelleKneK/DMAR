test_that("expected_smd() returns the documented columns and is upward-biased", {
  res <- expected_smd(delta = 0.5, n_1 = 10)
  expect_named(res, c("delta", "n_1", "n_2", "expected_smd",
                      "bias", "j_correction"))
  expect_gt(res$expected_smd, res$delta)  # upward bias
  expect_lt(res$j_correction, 1)
})

test_that("expected_smd() returns a dmar_tbl", {
  expect_s3_class(expected_smd(delta = 0.5, n_1 = 10), "dmar_tbl")
})

test_that("expected_smd() bias relationship: E[d] = delta / J(df)", {
  res <- expected_smd(delta = 0.5, n_1 = 10, n_2 = 10)
  df <- res$n_1 + res$n_2 - 2
  J  <- gamma(df / 2) / (sqrt(df / 2) * gamma((df - 1) / 2))
  expect_equal(res$expected_smd, 0.5 / J, tolerance = 1e-12)
  expect_equal(res$j_correction, J, tolerance = 1e-12)
})

test_that("expected_smd() handles n_1 = n_2 default (balanced)", {
  res_a <- expected_smd(delta = 0.5, n_1 = 20)
  res_b <- expected_smd(delta = 0.5, n_1 = 20, n_2 = 20)
  expect_equal(res_a$expected_smd, res_b$expected_smd, tolerance = 1e-12)
})

test_that("expected_smd() bias shrinks to zero as n grows", {
  res <- expected_smd(delta = 0.5, n_1 = c(5, 10, 20, 50, 100, 500))
  expect_true(all(diff(res$bias) < 0))  # monotone toward zero
  expect_lt(res$bias[length(res$bias)], 1e-3)
})

test_that("expected_smd() handles negative delta symmetrically", {
  pos <- expected_smd(delta =  0.5, n_1 = 10)$expected_smd
  neg <- expected_smd(delta = -0.5, n_1 = 10)$expected_smd
  expect_equal(pos, -neg, tolerance = 1e-12)
})

test_that("expected_smd() rejects bad inputs", {
  expect_error(expected_smd(delta = "a", n_1 = 10), "numeric")
  expect_error(expected_smd(delta = 0.5, n_1 = 1),  ">=")
})
