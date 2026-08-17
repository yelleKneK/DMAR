test_that("ss_aipe_equivalence_smd() returns documented rows", {
  res <- ss_aipe_equivalence_smd(population_smd = 0, width = 0.20)
  expect_setequal(res$term,
                  c("necessary_n_per_group", "total_N", "width",
                    "population_smd", "ci_width_expected"))
})

test_that("ss_aipe_equivalence_smd() recommends a finite n", {
  res <- ss_aipe_equivalence_smd(population_smd = 0, width = 0.20)
  n   <- res$value[res$term == "necessary_n_per_group"]
  expect_gt(n, 2)
  expect_lt(n, 1e5)
})

test_that("ss_aipe_equivalence_smd() expected width is at most target", {
  res <- ss_aipe_equivalence_smd(population_smd = 0, width = 0.20)
  w   <- res$value[res$term == "ci_width_expected"]
  expect_lte(w, 0.20 + 1e-8)
})

test_that("ss_aipe_equivalence_smd() smaller width => larger n", {
  n1 <- ss_aipe_equivalence_smd(population_smd = 0, width = 0.30)$value[1]
  n2 <- ss_aipe_equivalence_smd(population_smd = 0, width = 0.10)$value[1]
  expect_lt(n1, n2)
})
