
test_that("ss_aipe_equivalence_smd_sensitivity maps alpha to (1 - conf_level)/2 (HIGH-03)", {
  skip_on_cran()
  # The planner ss_aipe_equivalence_smd() sizes the width of a (1 - 2*alpha) interval,
  # and the simulation scores ci_smd() at conf_level, so the resolved sample
  # size must equal the planner called with alpha = (1 - conf_level)/2. The bug
  # used alpha = 1 - conf_level, planning a 90% width while scoring a 95% one.
  for (cl in c(0.90, 0.95, 0.99)) {
    oracle <- ss_aipe_equivalence_smd(population_smd = 0, width = 0.20,
                               alpha_level = (1 - cl) / 2, balanced = TRUE)
    oracle_n <- oracle$value[oracle$term == "necessary_n_per_group"]
    s <- ss_aipe_equivalence_smd_sensitivity(true_smd = 0, estimated_smd = 0,
                                      width = 0.20, delta_upper = 0.20,
                                      conf_level = cl, G = 1)
    expect_equal(s$value[s$term == "n_per_group"], oracle_n)
  }
})
