test_that("ci_rmsea() returns the expected 3-row tidy data frame", {
  rmsea_pt <- sqrt(30 / (15 * 199))
  result <- ci_rmsea(rmsea = rmsea_pt, df = 15, N = 200, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, c("lower_limit", "rmsea", "upper_limit"))
})

test_that("ci_rmsea() lower_limit <= point estimate <= upper_limit", {
  rmsea_pt <- sqrt(30 / (15 * 199))
  result <- ci_rmsea(rmsea = rmsea_pt, df = 15, N = 200, conf_level = 0.95)
  ll <- result$value[result$term == "lower_limit"]
  est <- result$value[result$term == "rmsea"]
  ul <- result$value[result$term == "upper_limit"]
  expect_gte(ll, 0)
  expect_lte(ll, est)
  expect_lte(est, ul)
})

test_that("ci_rmsea() clamps lower_limit to 0 for a well-fitting (small chi^2) model", {
  rmsea_pt <- sqrt(5 / (15 * 199))
  expect_message(
    result <- ci_rmsea(rmsea = rmsea_pt, df = 15, N = 200, conf_level = 0.95),
    "lower RMSEA limit is set to 0"
  )
  expect_equal(result$value[result$term == "lower_limit"], 0)
})

test_that("the help page example's 90 percent upper limit sits just above 0.05", {
  # Anchors the ?ci_rmsea example prose: for rmsea = .035, df = 40, N = 425
  # the 90 percent upper limit is 0.052, just above the Browne and Cudeck
  # close fit threshold of 0.05 (an earlier version of the page described it
  # as landing below 0.05). Recomputed independently by inverting the
  # noncentral chi square lower tail with uniroot on pchisq directly.
  suppressMessages(
    res <- ci_rmsea(rmsea = .035, df = 40, N = 425, conf_level = .90)
  )
  ul <- res$value[res$term == "upper_limit"]
  chi_sq <- .035^2 * 40 * 424 + 40
  lam <- uniroot(function(l) pchisq(chi_sq, df = 40, ncp = l) - 0.05,
                 c(0, 500), tol = 1e-10)$root
  expect_equal(ul, sqrt(lam / (40 * 424)), tolerance = 1e-6)
  expect_gt(ul, 0.05)
  expect_lt(ul, 0.053)
})
