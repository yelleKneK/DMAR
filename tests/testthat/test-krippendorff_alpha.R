test_that("krippendorff_alpha() returns documented rows", {
  ratings <- matrix(c(
    1, 2, 3, 3, 2, 1, 4, 1, 2, NA, NA, NA,
    1, 2, 3, 3, 2, 2, 4, 1, 2, 5,  NA, 3,
    NA, 3, 3, 3, 2, 3, 4, 2, 2, 5,  1,  NA,
    1, 2, 3, 3, 2, 4, 4, 1, 2, 5,  1,  NA
  ), nrow = 12, ncol = 4)
  res <- krippendorff_alpha(ratings, level = "nominal", boot = FALSE)
  expect_setequal(res$term,
                  c("krippendorff_alpha", "D_observed", "D_expected",
                    "n_pairable"))
})

test_that("krippendorff_alpha() is 1.0 on perfect agreement", {
  # All four raters give identical scores on every unit.
  R <- matrix(rep(1:5, times = 4), 5, 4)
  res <- krippendorff_alpha(R, boot = FALSE)
  expect_equal(res$value[res$term == "krippendorff_alpha"], 1,
               tolerance = 1e-10)
})

test_that("krippendorff_alpha() interval level returns finite alpha", {
  set.seed(113)
  r1 <- rnorm(30); r2 <- r1 + rnorm(30, 0, 0.2)
  res <- krippendorff_alpha(cbind(r1, r2), level = "interval",
                            boot = FALSE)
  expect_true(is.finite(res$value[1]))
})

test_that("krippendorff_alpha() bootstrap returns CI bounds", {
  set.seed(113)
  R <- matrix(sample(1:4, 40, TRUE), 20, 2)
  res <- krippendorff_alpha(R, boot = TRUE, B = 200L)
  expect_true("lower_limit" %in% res$term)
  expect_true("upper_limit" %in% res$term)
})

test_that("krippendorff_alpha() ratio level handles an identical zero pair (HIGH-05)", {
  # An identical zero pair (0, 0) makes the ratio distance ((a - b) / (a + b))^2
  # evaluate 0 / 0 = NaN, which used to crash the comparison. The oracle is
  # irr's independent implementation.
  # Pinned from irr::kripp.alpha (irr 0.85, 2026-08-09); live comparison in
  # tools/oracle_checks.R.
  ratings <- rbind(c(0, 0), c(1, 1), c(0, 1), c(1, 0))
  oracle  <- 0.125
  res     <- krippendorff_alpha(ratings, level = "ratio", boot = FALSE)
  alpha   <- res$value[res$term == "krippendorff_alpha"]
  expect_false(is.nan(alpha))
  expect_equal(alpha, oracle, tolerance = 1e-8)
})

test_that("the bootstrap is opt-in: no analysis runs one unless requested", {
  set.seed(113)
  r1 <- rnorm(20)
  r2 <- r1 + rnorm(20, 0, 0.4)
  res <- krippendorff_alpha(cbind(r1, r2), level = "interval")
  expect_false(any(c("lower_limit", "upper_limit", "B_used") %in% res$term))
})
