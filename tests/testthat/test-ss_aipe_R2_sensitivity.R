## The per-replication file of ss_aipe_R2_sensitivity(): its columns are
## named as the rows that summarize them (ci_width_lower, ci_width_upper,
## ci_width; DMAR 1.0.x wrote lower_width_ci, upper_width_ci, width_ci),
## and the rows are the summaries of those columns.

test_that("the per-replication widths carry the names of the rows that summarize them", {
  skip_on_cran()
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  res <- suppressMessages(ss_aipe_R2_sensitivity(
    true_R2 = 0.3, specified_N = 60, width = 0.3, p = 2, G = 40,
    print_iter = FALSE, filename = path, seed = 113))
  f <- utils::read.csv(path)
  expect_named(f, c("ci_lower", "R2", "ci_upper", "ci_width_lower",
                    "ci_width_upper", "ci_width"))
  v <- stats::setNames(res$value, res$term)
  for (col in c("ci_width_lower", "ci_width_upper", "ci_width")) {
    expect_equal(v[[paste0("mean_", col)]], mean(f[[col]], na.rm = TRUE))
    expect_equal(v[[paste0("median_", col)]], median(f[[col]], na.rm = TRUE))
    expect_equal(v[[paste0("sd_", col)]], sd(f[[col]], na.rm = TRUE))
  }
  # The widths are the distances from the observed R2 to each limit.
  expect_equal(f$ci_width_lower, f$R2 - f$ci_lower)
  expect_equal(f$ci_width_upper, f$ci_upper - f$R2)
  expect_equal(f$ci_width, f$ci_upper - f$ci_lower)
  expect_false(any(c("lower_width_ci", "upper_width_ci", "width_ci") %in% names(f)))
})

test_that("true_R2, width, or p left out is refused by name", {
  skip_on_cran()
  # Each is needed whether the sample size is planned from supposed_R2 or
  # given as specified_N, so a call without it is refused in its own name
  # before anything is planned or simulated.
  refused <- function(expr, arg) {
    e <- tryCatch(expr, error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_s3_class(e, "dmar_refusal")
    expect_identical(e$arg, arg)
    expect_match(conditionMessage(e), sprintf("You need to specify '%s'", arg),
                 fixed = TRUE)
  }
  refused(ss_aipe_R2_sensitivity(true_R2 = .5, supposed_R2 = .5, p = 5,
                                 G = 5, print_iter = FALSE), "width")
  refused(ss_aipe_R2_sensitivity(true_R2 = .5, specified_N = 100, p = 5,
                                 G = 5, print_iter = FALSE), "width")
  refused(ss_aipe_R2_sensitivity(supposed_R2 = .5, width = .1, p = 5,
                                 G = 5, print_iter = FALSE), "true_R2")
  refused(ss_aipe_R2_sensitivity(specified_N = 100, width = .1, p = 5,
                                 G = 5, print_iter = FALSE), "true_R2")
  refused(ss_aipe_R2_sensitivity(true_R2 = .5, specified_N = 100, width = .1,
                                 G = 5, print_iter = FALSE), "p")
  refused(ss_aipe_R2_sensitivity(true_R2 = .5, supposed_R2 = .5, width = .1,
                                 G = 5, print_iter = FALSE), "p")
})
