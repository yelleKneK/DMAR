test_that("ss_aipe_sc_ancova() returns a tidy single-row data frame", {
  res <- ss_aipe_sc_ancova(psi_standardized = .8, width = .5, c_weights = c(.5, .5, 0, -1))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("term", "value"))
  expect_equal(res$term, "necessary_n_per_group")
  expect_true(is.numeric(res$value) && res$value > 0)
})

test_that("ss_aipe_sc_ancova() ratio path returns a smaller N than the unconditioned plan", {
  # ratio < 1 represents an effective reduction in residual variance from a
  # covariate, so fewer subjects are needed.
  base   <- ss_aipe_sc_ancova(psi_standardized = .8, width = .5, c_weights = c(.5, .5, 0, -1))
  ratio  <- ss_aipe_sc_ancova(psi_standardized = .8, ratio = .6, width = .5,
                              c_weights = c(.5, .5, 0, -1), divisor = "s_anova")
  expect_lt(ratio$value, base$value)
})

test_that("ss_aipe_sc_ancova() accepts explicit alpha_lower/alpha_upper (s_ancova)", {
  ref <- ss_aipe_sc_ancova(psi_standardized = .8, width = .5, c_weights = c(.5, .5, 0, -1))
  alt <- ss_aipe_sc_ancova(psi_standardized = .8, width = .5, c_weights = c(.5, .5, 0, -1),
                           conf_level = NULL,
                           alpha_lower = .025, alpha_upper = .025)
  expect_equal(alt$value, ref$value)
})

test_that("ss_aipe_sc_ancova() accepts asymmetric alpha bounds (s_anova path)", {
  res <- ss_aipe_sc_ancova(psi_standardized = .8, ratio = .6, width = .5,
                           c_weights = c(.5, .5, 0, -1), divisor = "s_anova",
                           conf_level = NULL,
                           alpha_lower = .01, alpha_upper = .04)
  expect_s3_class(res, "data.frame")
  expect_equal(res$term, "necessary_n_per_group")
})

test_that("ss_aipe_sc_ancova() rejects mixing conf_level and alphas", {
  expect_error(
    ss_aipe_sc_ancova(psi_standardized = .8, width = .5, c_weights = c(.5, .5, 0, -1),
                      conf_level = .95, alpha_lower = .025, alpha_upper = .025),
    "cannot mix them"
  )
})

test_that("ss_aipe_sc_ancova() s_anova assurance path runs end-to-end", {
  # The s_anova + assurance branch previously also forwarded conf_level alongside
  # alphas to conf_limits_nct, which is the same bug pattern as s_ancova.
  res <- ss_aipe_sc_ancova(psi_standardized = .8, ratio = .6, width = .5,
                           c_weights = c(.5, .5, 0, -1), divisor = "s_anova",
                           assurance = .90)
  expect_s3_class(res, "data.frame")
  expect_equal(res$term, "necessary_n_per_group")
})

test_that("ss_aipe_sc_ancova() assurance path no longer errors on inner conf_limits_nct calls", {
  # Regression test: same bug as ss_aipe_sc() -- the inner conf_limits_nct
  # calls forwarded the default conf_level = .95 alongside alpha_lower /
  # alpha_upper, which conf_limits_nct refuses.
  res <- ss_aipe_sc_ancova(psi_standardized = .8, width = .5,
                           c_weights = c(.5, .5, 0, -1), assurance = .90)
  expect_s3_class(res, "data.frame")
  expect_equal(res$term, "necessary_n_per_group")
})
