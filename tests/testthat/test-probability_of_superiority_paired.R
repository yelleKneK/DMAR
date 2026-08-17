test_that("probability_of_superiority_paired() returns documented rows", {
  set.seed(113)
  pre <- rnorm(30, 100, 15); post <- pre + rnorm(30, 5, 10)
  res <- probability_of_superiority_paired(pre, post)
  expect_setequal(res$term,
                   c("probability_of_superiority", "lower_limit",
                     "upper_limit", "var_ps",
                     "wins_y_over_x", "losses_y_under_x", "ties"))
})

test_that("probability_of_superiority_paired() PS bounded in [0, 1]", {
  set.seed(113)
  res <- probability_of_superiority_paired(rnorm(30), rnorm(30, 0.5))
  expect_gte(res$value[res$term == "probability_of_superiority"], 0)
  expect_lte(res$value[res$term == "probability_of_superiority"], 1)
  expect_gte(res$value[res$term == "lower_limit"], 0)
  expect_lte(res$value[res$term == "upper_limit"], 1)
})

test_that("probability_of_superiority_paired() PS = 0.5 when paired diffs are symmetric around 0", {
  set.seed(113)
  d <- rnorm(1000, 0, 1)
  pre <- runif(1000); post <- pre + d
  res <- probability_of_superiority_paired(pre, post)
  expect_equal(res$value[res$term == "probability_of_superiority"],
               0.5, tolerance = 0.05)
})

test_that("probability_of_superiority_paired() handles all wins (PS = 1)", {
  expect_equal(
    probability_of_superiority_paired(1:10, 11:20)$value[1],
    1, tolerance = 1e-12
  )
})

test_that("probability_of_superiority_paired() rejects mismatched lengths", {
  expect_error(probability_of_superiority_paired(1:5, 1:6), "same length")
})
