## Tests for proportion_of_superiority() (sometimes called Cohen's U3).

test_that("proportion_of_superiority() point estimate equals pnorm(d)", {
  for (d in c(0.2, 0.5, 0.8, -0.5)) {
    res <- proportion_of_superiority(d)
    expect_equal(res$value[res$term == "proportion_of_superiority"],
                 pnorm(d), tolerance = 1e-12)
  }
})

test_that("proportion_of_superiority() values match Cohen (1988) Table 2.2.1", {
  # Cohen (1988) Table 2.2.1: d = 0.2, 0.5, 0.8 -> 58%, 69%, 79%.
  expect_equal(round(proportion_of_superiority(0.2)$value[2], 2), 0.58)
  expect_equal(round(proportion_of_superiority(0.5)$value[2], 2), 0.69)
  expect_equal(round(proportion_of_superiority(0.8)$value[2], 2), 0.79)
})

test_that("proportion_of_superiority() CI from sample sizes is monotone", {
  res <- proportion_of_superiority(smd = 0.5, n_1 = 50, n_2 = 50)
  expect_gt(res$value[res$term == "upper_limit"],
            res$value[res$term == "lower_limit"])
})

test_that("proportion_of_superiority() agrees with cles relationship for positive d", {
  # The proportion of superiority is pnorm(d); CLES is pnorm(d / sqrt(2)).
  # For positive d the proportion of superiority strictly exceeds CLES.
  ps <- proportion_of_superiority(0.5)$value[2]
  cl <- cles(0.5)$value[cles(0.5)$term == "cl"]
  expect_gt(ps, cl)
})

test_that("proportion_of_superiority() returns a tidy numeric data.frame", {
  res <- proportion_of_superiority(smd = 0.5, n_1 = 50, n_2 = 50)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_type(res$value, "double")
  expect_true(all(c("smd", "proportion_of_superiority",
                    "smd_lower", "smd_upper",
                    "lower_limit", "upper_limit") %in% res$term))
})
