test_that("expected_smd() returns the documented columns and is upward-biased", {
  res <- expected_smd(true_smd = 0.5, n_1 = 10)
  expect_named(res, c("true_smd", "n_1", "n_2", "expected_smd",
                      "bias", "j_correction"))
  expect_gt(res$expected_smd, res$true_smd)  # upward bias
  expect_lt(res$j_correction, 1)
})

test_that("expected_smd() returns a dmar_tbl", {
  expect_s3_class(expected_smd(true_smd = 0.5, n_1 = 10), "dmar_tbl")
})

test_that("expected_smd() bias relationship: E[d] = true_smd / J(df)", {
  res <- expected_smd(true_smd = 0.5, n_1 = 10, n_2 = 10)
  df <- res$n_1 + res$n_2 - 2
  J  <- gamma(df / 2) / (sqrt(df / 2) * gamma((df - 1) / 2))
  expect_equal(res$expected_smd, 0.5 / J, tolerance = 1e-12)
  expect_equal(res$j_correction, J, tolerance = 1e-12)
})

test_that("expected_smd() handles n_1 = n_2 default (balanced)", {
  res_a <- expected_smd(true_smd = 0.5, n_1 = 20)
  res_b <- expected_smd(true_smd = 0.5, n_1 = 20, n_2 = 20)
  expect_equal(res_a$expected_smd, res_b$expected_smd, tolerance = 1e-12)
})

test_that("expected_smd() bias shrinks to zero as n grows", {
  res <- expected_smd(true_smd = 0.5, n_1 = c(5, 10, 20, 50, 100, 500))
  expect_true(all(diff(res$bias) < 0))  # monotone toward zero
  expect_lt(res$bias[length(res$bias)], 1e-3)
})

test_that("expected_smd() handles a negative true_smd symmetrically", {
  pos <- expected_smd(true_smd =  0.5, n_1 = 10)$expected_smd
  neg <- expected_smd(true_smd = -0.5, n_1 = 10)$expected_smd
  expect_equal(pos, -neg, tolerance = 1e-12)
})

test_that("expected_smd() rejects bad inputs", {
  expect_error(expected_smd(true_smd = "a", n_1 = 10), "numeric")
  expect_error(expected_smd(true_smd = 0.5, n_1 = 1),  ">=")
})

test_that("expected_smd() keeps the bias at very large samples", {
  # Reference values from a 60-digit evaluation of true_smd (1 / J - 1). The
  # form used through 1.0.x returned -1.19e-07 (the wrong sign, with a J
  # above 1) in place of 1.875e-09 at n = 1e8 per group.
  r <- expected_smd(true_smd = 0.5, n_1 = 1e8, n_2 = 1e8)
  expect_equal(r$bias, 1.8750000285156254e-9, tolerance = 1e-13)
  expect_equal(r$expected_smd, 0.50000000187500003, tolerance = 1e-15)
  expect_equal(r$j_correction, 0.99999999624999996, tolerance = 1e-15)
  expect_equal(r$bias / 1.875e-09, 1, tolerance = 1e-6)
  r5 <- expected_smd(true_smd = 0.5, n_1 = 1e5, n_2 = 2e5)
  expect_equal(r5$bias, 1.2500126737397293e-6, tolerance = 1e-13)
  # Small samples are as before.
  r20 <- expected_smd(true_smd = 0.5, n_1 = 20)
  expect_equal(r20$expected_smd, 0.51014661094597607, tolerance = 1e-14)
  expect_equal(r20$bias, 0.010146610945976066, tolerance = 1e-13)
  expect_equal(r20$bias, r20$expected_smd - 0.5, tolerance = 1e-12)
  expect_identical(expected_smd(true_smd = 0, n_1 = 1e8)$bias, 0)
})
