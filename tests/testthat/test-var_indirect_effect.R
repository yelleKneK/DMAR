test_that("var_indirect_effect() returns all four formulas", {
  res <- var_indirect_effect(a = 0.40, b = 0.40, var_a = 0.02, var_b = 0.02)
  expect_setequal(res$term,
                   c("var_sobel", "var_aroian", "var_goodman",
                     "var_delta_second_order"))
  expect_true(all(res$value > 0))
})

test_that("var_indirect_effect() Sobel = Aroian - var_a*var_b = Goodman + var_a*var_b", {
  va <- 0.02; vb <- 0.02
  res <- var_indirect_effect(0.4, 0.4, va, vb)
  s   <- res$value[res$term == "var_sobel"]
  ar  <- res$value[res$term == "var_aroian"]
  g   <- res$value[res$term == "var_goodman"]
  expect_equal(ar - s, va * vb, tolerance = 1e-12)
  expect_equal(s - g, va * vb, tolerance = 1e-12)
})

test_that("var_indirect_effect() cov_ab affects only the second-order delta row", {
  base <- var_indirect_effect(0.4, 0.4, 0.02, 0.02, cov_ab = 0)
  with_cov <- var_indirect_effect(0.4, 0.4, 0.02, 0.02, cov_ab = 0.01)
  expect_equal(base$value[base$term == "var_sobel"],
               with_cov$value[with_cov$term == "var_sobel"])
  expect_gt(with_cov$value[with_cov$term == "var_delta_second_order"],
            base$value[base$term == "var_delta_second_order"])
})
