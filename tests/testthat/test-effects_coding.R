test_that("effects_coding() default reference is the last level", {
  M <- effects_coding(c("a", "b", "c"))
  expect_equal(rownames(M), c("a", "b", "c"))
  expect_equal(colnames(M), c("a", "b"))
  expect_equal(M["c", ], c(a = -1, b = -1))
})

test_that("effects_coding() with explicit reference", {
  M <- effects_coding(c("low", "med", "high"), reference = "low")
  expect_equal(M["low", ], c(med = -1, high = -1))
})

test_that("effects_coding() column sums are zero", {
  M <- effects_coding(4)
  expect_true(all(abs(colSums(M)) < 1e-12))
})

test_that("helmert_coding() yields an orthogonal contrast matrix", {
  M <- helmert_coding(5)
  res <- is_orthogonal_set(M)
  expect_equal(res$value[res$term == "all_orthogonal"], 1)
  expect_equal(res$value[res$term == "all_contrasts_sum_to_zero"], 1)
})

test_that("helmert_coding() level labels preserved", {
  M <- helmert_coding(c("base", "wk1", "wk2", "wk3"))
  expect_equal(rownames(M), c("base", "wk1", "wk2", "wk3"))
  expect_match(colnames(M), "_vs_prior")
})
