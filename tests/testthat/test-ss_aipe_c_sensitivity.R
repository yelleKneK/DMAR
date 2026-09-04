## ss_aipe_c_sensitivity() and the per-replication CSV writer it shares with
## the rest of the ss_aipe_*_sensitivity() family. The Monte Carlo schema
## and the monotonicity in the error variance are covered in
## test-ss_aipe_new_sensitivities.R and test-ss_aipe_sensitivity_family.R;
## this file covers the filename switch, the only route to a file since the
## 2026-09-04 CRAN review retired the save argument and its default path.

test_that(".write_sensitivity_csv() writes a header on a new file and appends rows without one", {
  per_rep <- data.frame(a = c(1.5, 2.5), b = c(TRUE, FALSE))
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)

  DMAR:::.write_sensitivity_csv(per_rep, path)
  back <- utils::read.csv(path)
  expect_named(back, c("a", "b"))
  expect_equal(nrow(back), 2L)
  expect_equal(back$a, per_rep$a)
  expect_equal(back$b, per_rep$b)

  # A second call appends under the existing header rather than repeating it.
  DMAR:::.write_sensitivity_csv(per_rep, path)
  back <- utils::read.csv(path)
  expect_named(back, c("a", "b"))
  expect_equal(nrow(back), 4L)
  expect_equal(back$a, rep(per_rep$a, 2))
  expect_equal(length(readLines(path)), 5L)
})

test_that("ss_aipe_c_sensitivity() writes one row per replication to filename", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  set.seed(113)
  res <- ss_aipe_c_sensitivity(true_error_variance = 4, n_per_group = 10,
                               c_weights = c(-1, 0, 1), width = 1,
                               G = 7, print_iter = FALSE, filename = path)
  expect_true(file.exists(path))
  per_rep <- utils::read.csv(path)
  expect_named(per_rep, c("psi_hat", "ci_lower", "ci_upper", "ci_width",
                          "type_I_lower", "type_I_upper"))
  expect_equal(nrow(per_rep), 7L)
  # The file holds the replications the summary rows were computed from.
  v <- stats::setNames(res$value, res$term)
  expect_equal(mean(per_rep$ci_width), unname(v["mean_ci_width"]))
  expect_equal(mean(per_rep$psi_hat), unname(v["mean_psi"]))
  expect_equal(mean(per_rep$type_I_lower | per_rep$type_I_upper),
               unname(v["total_type_I_error"]))

  # Running again against the same path appends the new replications.
  set.seed(113)
  ss_aipe_c_sensitivity(true_error_variance = 4, n_per_group = 10,
                        c_weights = c(-1, 0, 1), width = 1,
                        G = 7, print_iter = FALSE, filename = path)
  expect_equal(nrow(utils::read.csv(path)), 14L)
})

test_that("ss_aipe_c_sensitivity() writes nothing by default and rejects a bad filename", {
  before <- list.files(getwd(), all.files = TRUE)
  set.seed(113)
  res <- ss_aipe_c_sensitivity(true_error_variance = 4, n_per_group = 10,
                               c_weights = c(-1, 0, 1), width = 1,
                               G = 3, print_iter = FALSE)
  expect_s3_class(res, "data.frame")
  expect_identical(list.files(getwd(), all.files = TRUE), before)

  for (bad in list(TRUE, 1, c("a.csv", "b.csv"), NA_character_, "")) {
    expect_error(
      ss_aipe_c_sensitivity(true_error_variance = 4, n_per_group = 10,
                            c_weights = c(-1, 0, 1), width = 1,
                            G = 3, print_iter = FALSE, filename = bad),
      "'filename' must be NULL or a single character string"
    )
  }
})
