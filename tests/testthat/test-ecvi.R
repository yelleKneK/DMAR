test_that("ecvi() matches the Browne-Cudeck / lavaan point estimate", {
  # (chisq + 2*npar)/n.
  res <- ecvi(chisq = 24.361, df = 8, npar = 13, n = 301)
  expect_s3_class(res, "dmar_tbl")
  expect_equal(res$value[res$term == "ecvi"], (24.361 + 2 * 13) / 301)
  # The interval brackets the point estimate.
  lo <- res$value[res$term == "lower_limit"]
  hi <- res$value[res$term == "upper_limit"]
  est <- res$value[res$term == "ecvi"]
  expect_lte(lo, est)
  expect_gte(hi, est)
  expect_identical(attr(res, "conf_level"), 0.95)
})

test_that("ecvi(fit) agrees with lavaan::fitMeasures", {
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa(
    "visual =~ t1_visual_perception + t2_cubes + t4_lozenges
     verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    data = holzinger_swineford, std.lv = TRUE)
  res <- ecvi(fit)
  oracle <- unname(lavaan::fitMeasures(fit, "ecvi"))
  expect_equal(res$value[res$term == "ecvi"], oracle, tolerance = 1e-6)
})

test_that("ecvi() validates input", {
  expect_error(ecvi(chisq = 10, df = 0, npar = 5, n = 100), "at least 1")
  expect_error(ecvi(chisq = 10, df = 5, npar = 5, n = -1), "non-negative")
  expect_error(ecvi(chisq = 10, df = 5, npar = 5, n = 100, conf_level = 1),
               "\\(0, 1\\)")
})
