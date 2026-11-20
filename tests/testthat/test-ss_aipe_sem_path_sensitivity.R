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
    model = analysis_model, supposed_cov_matrix = Sigma, true_cov_matrix = Sigma,
    which_path = "b", width = 0.30, specified_N = 150, G = 15
  )

  expect_s3_class(out, "dmar_tbl")

  expect_gt(attr(out, "simulation")$returned, 0)
  expect_false("suc_rep" %in% out$term)

  mean_width <- out$value[out$term == "mean_ci_width"]
  expect_true(is.finite(mean_width))
  expect_gt(mean_width, 0)
})

test_that("ss_aipe_sem_path_sensitivity() writes per-replication results only to a named file", {
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

  # Nothing is written unless the caller names a file; the old `save`
  # switch and its default path in the working directory are gone, and a
  # malformed filename is refused before any replication runs.
  expect_null(eval(formals(ss_aipe_sem_path_sensitivity)$filename))
  expect_false("save" %in% names(formals(ss_aipe_sem_path_sensitivity)))
  expect_error(
    ss_aipe_sem_path_sensitivity(model = analysis_model, supposed_cov_matrix = Sigma,
                                 which_path = "b", width = 0.30,
                                 specified_N = 150, G = 3,
                                 filename = c("a.csv", "b.csv")),
    "'filename' must be NULL or a single character string", fixed = TRUE)

  out_file <- tempfile(fileext = ".csv")
  on.exit(unlink(out_file), add = TRUE)
  set.seed(113)
  ss_aipe_sem_path_sensitivity(model = analysis_model, supposed_cov_matrix = Sigma,
                               which_path = "b", width = 0.30,
                               specified_N = 150, G = 3, filename = out_file)
  written <- utils::read.csv(out_file)
  expect_named(written, c("theta_hat", "se_theta_hat", "ci_lower", "ci_upper", "width"))
  expect_equal(nrow(written), 3L)
  expect_equal(written$width, written$ci_upper - written$ci_lower)

  # A second run appends to the existing file and says so.
  set.seed(113)
  expect_message(
    ss_aipe_sem_path_sensitivity(model = analysis_model, supposed_cov_matrix = Sigma,
                                 which_path = "b", width = 0.30,
                                 specified_N = 150, G = 3, filename = out_file),
    "already exists", fixed = TRUE)
  expect_equal(nrow(utils::read.csv(out_file)), 6L)
})

test_that("ss_aipe_sem_path_sensitivity() passes the arguments in '...' to every fit", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")

  # With std.lv = TRUE the path is 0.5 / sqrt(0.75) = 0.577, against the
  # 0.5 of the marker parameterization: the fit to true_cov_matrix, which gives
  # the population value covered, and the fit to each sample estimate it
  # in the parameterization of the plan.
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
    model = analysis_model, supposed_cov_matrix = Sigma, which_path = "b",
    width = 0.30, specified_N = 1000, G = 10, std.lv = TRUE
  )
  v <- function(term) out$value[out$term == term]
  expect_equal(v("true_path"), 0.5 / sqrt(0.75), tolerance = 1e-4)
  expect_lt(abs(v("mean_path") - 0.5 / sqrt(0.75)), abs(v("mean_path") - 0.5))
})
