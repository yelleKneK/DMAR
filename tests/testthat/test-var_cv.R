test_that("var_cv() returns McKay and Vangel rows", {
  res <- var_cv(0.20, 30)
  expect_setequal(res$term, c("var_cv_mckay", "var_cv_vangel"))
  expect_true(all(res$value > 0))
})

test_that("var_cv() Vangel exceeds McKay for non-trivial cv (small-n correction)", {
  res <- var_cv(0.50, 30)
  expect_gt(res$value[res$term == "var_cv_vangel"],
            res$value[res$term == "var_cv_mckay"])
})

test_that("var_cv() McKay row equals its closed form at cv = 0.20, n = 30", {
  cv <- 0.20; n <- 30
  mckay <- (cv^2 / (n - 1)) * (0.5 + cv^2)  # 0.000744827586206897
  res <- var_cv(cv, n)
  expect_equal(res$value[res$term == "var_cv_mckay"], mckay,
               tolerance = 1e-12)
})

test_that("var_cv() rejects bad inputs", {
  expect_error(var_cv(-0.1, 30), "positive")
  expect_error(var_cv(0.2,  2),  ">= 3")
})
