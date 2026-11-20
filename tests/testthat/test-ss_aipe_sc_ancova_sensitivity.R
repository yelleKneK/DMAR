## ss_aipe_sc_ancova_sensitivity() plans the sample size from the supposed
## value through ss_aipe_sc_ancova(), and passes it assurance_method (new in
## 1.1.0, "analytic" by default, as in the planner; Ken's decision of
## 2026-11-13). The planner's simulated step draws from the stream that the
## sensitivity function's seed sets, with the planner's default G of 10000;
## the G of the sensitivity function counts the replications it simulates.

test_that("ss_aipe_sc_ancova_sensitivity() plans through the planner's assurance method", {
  expect_identical(eval(formals(ss_aipe_sc_ancova_sensitivity)$assurance_method),
                   c("analytic", "simulation"))
  for (div in c("s_ancova", "s_anova")) {
    for (m in list(NULL, "analytic", "simulation")) {
      res <- do.call(ss_aipe_sc_ancova_sensitivity, c(list(
        true_std_contrast = .5, supposed_std_contrast = .5, c_weights = c(1, -1),
        width = .5, true_cor = .5, divisor = div, assurance = .8, G = 3,
        print_iter = FALSE, seed = 1), if (!is.null(m)) list(assurance_method = m)))
      planned <- do.call(ss_aipe_sc_ancova, c(list(
        true_std_contrast = .5, true_cor = .5, divisor = div,
        c_weights = c(1, -1), width = .5, assurance = .8, seed = 1),
        if (!is.null(m)) list(assurance_method = m)))
      v <- stats::setNames(res$value, res$term)
      expect_identical(unname(v["n_per_group"]), planned$value[1])
      expect_identical(unname(v["assurance"]), .8)
    }
  }
  # The default plans the published plan.
  res <- ss_aipe_sc_ancova_sensitivity(
    true_std_contrast = .5, supposed_std_contrast = .5, c_weights = c(1, -1),
    width = .5, true_cor = .5, divisor = "s_anova", assurance = .8, G = 3,
    print_iter = FALSE, seed = 1)
  expect_identical(res$value[res$term == "n_per_group"], 98)
  # It is not used with a specified sample size.
  res <- ss_aipe_sc_ancova_sensitivity(
    true_std_contrast = .5, specified_n_per_group = 20, c_weights = c(1, -1),
    width = .5, true_cor = .5, G = 3, print_iter = FALSE, seed = 1,
    assurance_method = "simulation")
  expect_identical(res$value[res$term == "n_per_group"], 20)
  err <- tryCatch(ss_aipe_sc_ancova_sensitivity(
    true_std_contrast = .5, specified_n_per_group = 20, c_weights = c(1, -1),
    width = .5, true_cor = .5, G = 3, print_iter = FALSE, assurance_method = NA),
    error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "assurance_method")
})

test_that("ss_aipe_sc_ancova_sensitivity() plans with assurance without a warning", {
  for (m in c("analytic", "simulation"))
    expect_no_warning(ss_aipe_sc_ancova_sensitivity(
      true_std_contrast = .5, supposed_std_contrast = .5, c_weights = c(1, -1),
      width = .5, true_cor = .5, divisor = "s_anova", assurance = .8, G = 3,
      print_iter = FALSE, seed = 1, assurance_method = m))
})
