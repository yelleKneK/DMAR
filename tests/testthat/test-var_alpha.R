test_that("var_alpha() returns van Zyl and Bonett rows", {
  res <- var_alpha(alpha = 0.80, n = 100, p_items = 10)
  expect_setequal(res$term, c("var_alpha_van_zyl", "var_alpha_bonett"))
  expect_true(all(res$value > 0))
})

test_that("var_alpha() variance shrinks with larger n", {
  small <- var_alpha(0.90, n = 30,  p_items = 5)$value[1]
  large <- var_alpha(0.90, n = 300, p_items = 5)$value[1]
  expect_gt(small, large)
})

test_that("var_alpha() approximations converge as n grows", {
  out <- var_alpha(0.90, n = 1000, p_items = 5)
  expect_equal(out$value[1], out$value[2], tolerance = 1e-3)
})

test_that("var_alpha() rejects bad inputs", {
  expect_error(var_alpha(1.2, 100, 10), "in \\[0, 1\\)")
  expect_error(var_alpha(0.8, 2,   10), ">= 3")
  expect_error(var_alpha(0.8, 100, 1),  ">= 2")
})
