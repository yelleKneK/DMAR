test_that("ss_aipe_sm() returns a 1-row tidy data frame with necessary_N", {
  result <- ss_aipe_sm(true_sm = 0.5, width = 0.4)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, "necessary_N")
})

test_that("ss_aipe_sm() requires a larger n for a narrower width", {
  n_wide   <- ss_aipe_sm(true_sm = 0.5, width = 0.5)$value
  n_narrow <- ss_aipe_sm(true_sm = 0.5, width = 0.3)$value
  expect_lt(n_wide, n_narrow)
})

test_that("ss_aipe_sm() is symmetric about zero in sm (sign doesn't matter)", {
  pos <- ss_aipe_sm(true_sm =  0.5, width = 0.4)$value
  neg <- ss_aipe_sm(true_sm = -0.5, width = 0.4)$value
  expect_equal(pos, neg)
})

test_that("ss_aipe_sm() assurance branch returns a sensible, finite sample size", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # The assurance branch's noncentrality was previously the reciprocal form
  # (sm / sqrt(n)), which made the assurance target unreachable and pinned the
  # search to a bracket endpoint. With the correct sm * sqrt(n) it solves and
  # an assurance of .85 requires at least as large an n as the no-assurance plan.
  n_plain <- ss_aipe_sm(true_sm = 0.5, width = 0.4)$value
  n_85    <- ss_aipe_sm(true_sm = 0.5, width = 0.4, assurance = 0.85)$value
  expect_true(is.finite(n_85) && n_85 > 1)
  expect_gte(n_85, n_plain)
})

test_that("ss_aipe_sm() takes any assurance strictly between 0 and 1", {
  skip_on_cran()
  # The value the assurance step plans for, the standardized mean that the
  # sample value exceeds in magnitude with probability 1 - assurance
  # (Kelley & Rausch, 2006, Equation 31, in its one-sample form), is
  # defined for every probability in (0, 1). The sample sizes are those the
  # function returned when the range was opened (2026-11-11); the realized
  # probabilities beside them, from the noncentral t distribution, are in
  # the validation record: .228, .519, .812, .925.
  n <- function(a) {
    ss_aipe_sm(true_sm = 1, width = 0.5, assurance = a)$value
  }
  expect_identical(ss_aipe_sm(true_sm = 1, width = 0.5)$value, 93)
  expect_identical(n(.20), 87)
  expect_identical(n(.50), 93)
  expect_identical(n(.80), 100)
  expect_identical(n(.90), 105)
  # The sample size does not fall as the assurance rises.
  grid <- vapply(c(.05, .20, .40, .50, .60, .80, .99), n, numeric(1))
  expect_false(is.unsorted(grid))
  # The probability that the interval is no wider than the width, at the
  # returned sample size: the width depends on the sample only through
  # |t| and rises with it, so the probability is that of |t| below the
  # value at which the width equals the target.
  realized <- function(N, true_sm = 1, width = 0.5) {
    w <- function(tv) {
      l <- ci_nc_t(ncp = tv, df = N - 1, conf_level = .95)
      diff(l$value[l$term %in% c("ci_lower", "ci_upper")]) / sqrt(N)
    }
    cut <- stats::uniroot(function(tv) w(tv) - width,
                          c(0, 3 * true_sm * sqrt(N)), tol = 1e-8)$root
    stats::pt(cut, N - 1, true_sm * sqrt(N)) -
      stats::pt(-cut, N - 1, true_sm * sqrt(N))
  }
  expect_equal(realized(87), .228, tolerance = .005)
  expect_equal(realized(93), .519, tolerance = .005)
  expect_equal(realized(100), .812, tolerance = .005)
  expect_equal(realized(105), .925, tolerance = .005)
})

test_that("ss_aipe_sm() refuses an assurance outside (0, 1) with a classed condition", {
  for (bad in list(0, 1, -.1, 1.2, NA_real_, c(.8, .9), "0.9")) {
    expect_error(ss_aipe_sm(true_sm = 1, width = 0.5, assurance = bad),
                 class = "dmar_assurance_invalid")
  }
})

test_that("ss_aipe_sm_sensitivity() writes per-replication results only when 'filename' is supplied", {
  skip_on_cran()
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  args <- list(true_sm = 0.5, specified_N = 30, width = 0.5,
               G = 3, print_iter = FALSE)
  set.seed(113)
  res_none <- do.call(ss_aipe_sm_sensitivity, args)
  expect_false(file.exists(csv))
  set.seed(113)
  res_file <- do.call(ss_aipe_sm_sensitivity, c(args, list(filename = csv)))
  expect_equal(res_file, res_none)
  written <- utils::read.csv(csv)
  expect_equal(nrow(written), 3L)
  expect_named(written, c("sm_obs", "full_width", "width_lower", "width_upper",
                          "type_I_error_upper", "type_I_error_lower",
                          "type_I_error", "ci_lower", "ci_upper"))
  # A second run appends to the existing file.
  set.seed(113)
  do.call(ss_aipe_sm_sensitivity, c(args, list(filename = csv)))
  expect_equal(nrow(utils::read.csv(csv)), 6L)
  expect_error(
    do.call(ss_aipe_sm_sensitivity, c(args, list(filename = 1))),
    "'filename' must be NULL or a single character string"
  )
})
