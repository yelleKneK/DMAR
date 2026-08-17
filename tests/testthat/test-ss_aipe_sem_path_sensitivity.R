test_that("ss_aipe_sem_path_sensitivity runs a lavaan-backed Monte Carlo study", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")

  pop_model <- "
    f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
    f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
    f2 ~ 0.5*f1
    f1 ~~ 1*f1
    f2 ~~ 0.75*f2
    y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
    y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
  "
  Sigma <- cov_sem(pop_model)$sigma_theta

  analysis_model <- "
    f1 =~ y1 + y2 + y3
    f2 =~ y4 + y5 + y6
    f2 ~ b*f1
  "

  set.seed(113)
  out <- ss_aipe_sem_path_sensitivity(
    model = analysis_model, est_Sigma = Sigma, true_Sigma = Sigma,
    which_path = "b", desired_width = 0.30, N = 150, G = 15
  )

  expect_s3_class(out, "dmar_tbl")

  suc_rep <- out$value[out$term == "suc_rep"]
  expect_gt(suc_rep, 0)

  mean_width <- out$value[out$term == "mean_ci_width"]
  expect_true(is.finite(mean_width))
  expect_gt(mean_width, 0)
})
