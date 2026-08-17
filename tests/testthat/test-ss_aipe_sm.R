test_that("ss_aipe_sm() returns a 1-row tidy data frame with necessary_N", {
  result <- ss_aipe_sm(sm = 0.5, width = 0.4)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, "necessary_N")
})

test_that("ss_aipe_sm() requires a larger n for a narrower width", {
  n_wide   <- ss_aipe_sm(sm = 0.5, width = 0.5)$value
  n_narrow <- ss_aipe_sm(sm = 0.5, width = 0.3)$value
  expect_lt(n_wide, n_narrow)
})

test_that("ss_aipe_sm() is symmetric about zero in sm (sign doesn't matter)", {
  pos <- ss_aipe_sm(sm =  0.5, width = 0.4)$value
  neg <- ss_aipe_sm(sm = -0.5, width = 0.4)$value
  expect_equal(pos, neg)
})

test_that("ss_aipe_sm() assurance branch returns a sensible, finite sample size", {
  # The assurance branch's noncentrality was previously the reciprocal form
  # (sm / sqrt(n)), which made the assurance target unreachable and pinned the
  # search to a bracket endpoint. With the correct sm * sqrt(n) it solves and
  # an assurance >= 0.5 requires at least as large an n as the no-assurance plan.
  n_plain <- ss_aipe_sm(sm = 0.5, width = 0.4)$value
  n_85    <- ss_aipe_sm(sm = 0.5, width = 0.4, assurance = 0.85)$value
  expect_true(is.finite(n_85) && n_85 > 1)
  expect_gte(n_85, n_plain)
})
