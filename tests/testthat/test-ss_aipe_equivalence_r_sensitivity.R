test_that("ss_aipe_equivalence_r_sensitivity() reports pct_equivalent inside the bounds", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_equivalence_r_sensitivity(
    true_r = 0, estimated_r = 0, width = 0.30,
    rho_upper = 0.20,
    G = 60, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("mean_r", "mean_ci_width", "pct_equivalent",
                    "total_N") %in% res$term))
  p_eq <- res$value[res$term == "pct_equivalent"]
  expect_gte(p_eq, 0)
  expect_lte(p_eq, 1)
})

test_that("ss_aipe_equivalence_r_sensitivity maps alpha to (1 - conf_level)/2", {
  skip_on_cran()
  # The planner sizes the width of a (1 - 2*alpha) interval, and the
  # simulation scores ci_r() at conf_level, so the resolved sample size
  # must equal the planner called with alpha = (1 - conf_level)/2.
  for (cl in c(0.90, 0.95, 0.99)) {
    oracle <- ss_aipe_equivalence_r(population_r = 0, width = 0.20,
                             alpha_level = (1 - cl) / 2)
    oracle_N <- oracle$value[oracle$term == "necessary_N"]
    s <- ss_aipe_equivalence_r_sensitivity(true_r = 0, estimated_r = 0,
                                    width = 0.20, rho_upper = 0.20,
                                    conf_level = cl, G = 1)
    expect_equal(s$value[s$term == "total_N"], oracle_N)
  }
})

test_that("ss_aipe_equivalence_r_sensitivity() validates its inputs", {
  expect_error(ss_aipe_equivalence_r_sensitivity(width = 0.2),
               "either 'estimated_r' or 'specified_N'")
  expect_error(ss_aipe_equivalence_r_sensitivity(estimated_r = 0, specified_N = 50,
                                          width = 0.2),
               "but not both")
  expect_error(ss_aipe_equivalence_r_sensitivity(estimated_r = 0, width = 0.2),
               "'rho_upper' must be specified")
  expect_error(ss_aipe_equivalence_r_sensitivity(estimated_r = 0, width = 0.2,
                                          rho_upper = 1.2),
               "must be in \\(0, 1\\)")
})
