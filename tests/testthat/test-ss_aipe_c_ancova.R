## ss_aipe_c_ancova() -- sample size for an unstandardized ANCOVA contrast.

test_that("ss_aipe_c_ancova() returns a tidy data.frame with necessary_n_per_group", {
  res <- ss_aipe_c_ancova(error_var_anova = 1, rho = 0.4,
                          c_weights = c(-1, 0, 1), width = 0.20)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_n_per_group" %in% res$term)
  expect_gt(res$value[res$term == "necessary_n_per_group"], 0)
})

test_that("ss_aipe_c_ancova() requires more N for tighter widths", {
  wide  <- ss_aipe_c_ancova(error_var_anova = 1, rho = 0.4,
                            c_weights = c(-1, 0, 1), width = 0.40)$value[1]
  tight <- ss_aipe_c_ancova(error_var_anova = 1, rho = 0.4,
                            c_weights = c(-1, 0, 1), width = 0.10)$value[1]
  expect_gt(tight, wide)
})

test_that("ss_aipe_c_ancova() requires fewer N when the covariate explains more variance", {
  weak_rho   <- ss_aipe_c_ancova(error_var_anova = 1, rho = 0.10,
                                 c_weights = c(-1, 0, 1), width = 0.20)$value[1]
  strong_rho <- ss_aipe_c_ancova(error_var_anova = 1, rho = 0.80,
                                 c_weights = c(-1, 0, 1), width = 0.20)$value[1]
  expect_gt(weak_rho, strong_rho)
})

test_that("ss_aipe_c_ancova() rejects both ANCOVA-and-ANOVA error specification simultaneously", {
  expect_error(
    ss_aipe_c_ancova(error_var_ancova = 0.8, error_var_anova = 1, rho = 0.4,
                     c_weights = c(-1, 0, 1), width = 0.20),
    "do not input"
  )
})

test_that("ss_aipe_c_ancova() floors an impossibly wide target at the admissible minimum", {
  # A huge target width previously crashed the fixed-point search ("missing
  # value where TRUE/FALSE needed") because the error df J * (n - 1) - 1 went
  # negative. The smallest admissible design is n = 2 per group.
  r <- ss_aipe_c_ancova(error_var_anova = 40, rho = 0.22,
                        c_weights = c(1, -.5, -.5), width = 100)
  expect_equal(r$value[r$term == "necessary_n_per_group"], 2)
})

test_that("ss_aipe_c_ancova() validates its boundary inputs", {
  expect_error(ss_aipe_c_ancova(error_var_anova = 40, rho = 0.22,
                                c_weights = c(1, -.5, -.5), width = 3,
                                conf_level = 1.2), "conf_level")
  expect_error(ss_aipe_c_ancova(error_var_anova = 40, rho = 0.22,
                                c_weights = c(1, -.5, -.5), width = 0), "width")
  expect_error(ss_aipe_c_ancova(error_var_anova = 40, rho = 1,
                                c_weights = c(1, -.5, -.5), width = 3), "rho")
  expect_error(ss_aipe_c_ancova(error_var_anova = -1, rho = 0.22,
                                c_weights = c(1, -.5, -.5), width = 3),
               "error_var_anova")
})
