test_that("ss_aipe_sem_path plans the expected sample size for a labeled path", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

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

  out <- ss_aipe_sem_path(model = analysis_model, Sigma = Sigma,
                          desired_width = 0.30, which_path = "b")

  expect_s3_class(out, "dmar_tbl")

  necessary_N <- out$value[out$term == "necessary_N"]
  expect_equal(necessary_N, 264)

  path_index <- out$value[out$term == "path_index"]
  expect_true(path_index > 0)
  expect_equal(path_index, round(path_index))

  var_theta_j <- out$value[out$term == "var_theta_j"]
  expect_true(is.finite(var_theta_j))
  expect_true(var_theta_j > 0)
})

test_that("ss_aipe_sem_path errors when which_path is not a labeled parameter", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

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

  expect_error(
    ss_aipe_sem_path(model = analysis_model, Sigma = Sigma,
                     desired_width = 0.30, which_path = "not_a_label"),
    "not a labeled parameter"
  )
})
