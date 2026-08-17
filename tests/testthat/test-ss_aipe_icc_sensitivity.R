## ss_aipe_icc_sensitivity() -- Monte Carlo sensitivity check for the AIPE
## sample size planner on the intraclass correlation coefficient.

test_that("ss_aipe_icc_sensitivity() returns a tidy data.frame with the expected schema", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_icc_sensitivity(
    true_rho      = 0.70,
    estimated_rho = 0.70,
    k             = 3,
    width         = 0.20,
    conf_level    = 0.95,
    G             = 100,
    print_iter    = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_type(res$value, "double")
  expect_true(all(c("mean_icc", "median_icc", "sd_icc",
                    "mean_ci_width", "median_ci_width", "sd_ci_width",
                    "pct_ci_less_w", "pct_ci_miss_low", "pct_ci_miss_high",
                    "total_type_I_error",
                    "total_N", "k", "true_rho",
                    "estimated_rho", "width", "conf_level") %in% res$term))
  expect_identical(attr(res, "icc_type"), "ICC(1,1)")
})

test_that("ss_aipe_icc_sensitivity() mean realized ICC tracks true_rho when well-specified", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_icc_sensitivity(
    true_rho      = 0.50,
    estimated_rho = 0.50,
    k             = 3,
    width         = 0.30,
    G             = 200,
    print_iter    = FALSE
  )
  mean_icc <- res$value[res$term == "mean_icc"]
  expect_gt(mean_icc, 0.30)
  expect_lt(mean_icc, 0.70)
})

test_that("ss_aipe_icc_sensitivity() validates 'type'", {
  # The message pattern pins the match.arg validation: the pre-fix code
  # ignored 'type' until icc() failed inside the simulation loop with a
  # different message ("Unrecognized ICC type(s): ...").
  expect_error(
    ss_aipe_icc_sensitivity(true_rho = 0.5, specified_N = 40,
                            k = 3, width = 0.30, G = 5, type = "garbage"),
    "should be one of"
  )
})

test_that("ss_aipe_icc_sensitivity() generates data at the average-of-k level", {
  # skip_on_cran: the fast anchor that stays on CRAN is the deterministic
  # single-rater/average-of-k planning split in test-ss_aipe_icc.R
  # ("planning is shared within the single-rater and average-of-k
  # families"); this test confirms the generator side by Monte Carlo.
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_icc_sensitivity(
    true_rho    = 0.70,
    specified_N = 60,
    k           = 3,
    width       = 0.30,
    type        = "ICC(1,k)",
    G           = 400,
    print_iter  = FALSE
  )
  mean_icc <- res$value[res$term == "mean_icc"]
  # true_rho is on the ICC(1,k) scale. With the inverse Spearman-Brown
  # mapping in place, the mean realized ICC(1,k) sits near .70; generating
  # at the single-rater level instead would put it near
  # 3 * .70 / (1 + 2 * .70) = .875.
  expect_lt(abs(mean_icc - 0.70), 0.05)
  expect_identical(attr(res, "icc_type"), "ICC(1,k)")
})

test_that("ss_aipe_icc_sensitivity() rejects bad inputs", {
  expect_error(
    ss_aipe_icc_sensitivity(true_rho = 0.5, k = 3, width = 0.20, G = 10),
    "either 'estimated_rho' or 'specified_N'"
  )
  expect_error(
    ss_aipe_icc_sensitivity(true_rho = 0.5, estimated_rho = 0.5, specified_N = 50,
                            k = 3, width = 0.20, G = 10),
    "but not both"
  )
  expect_error(
    ss_aipe_icc_sensitivity(estimated_rho = 0.5, k = 3, width = 0.20, G = 10),
    "'true_rho' must be"
  )
  expect_error(
    ss_aipe_icc_sensitivity(true_rho = 1.2, estimated_rho = 0.5,
                            k = 3, width = 0.20, G = 10),
    "'true_rho' must be"
  )
})

test_that("ss_aipe_icc_sensitivity() supports specified_N mode", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_icc_sensitivity(
    true_rho    = 0.50,
    specified_N = 60,
    k           = 3,
    width       = 0.30,
    G           = 100,
    print_iter  = FALSE
  )
  expect_equal(res$value[res$term == "total_N"], 60)
  expect_true(is.na(res$value[res$term == "estimated_rho"]))
})

test_that("ss_aipe_icc_sensitivity() supports the two-way random and mixed types", {
  skip_on_cran()
  set.seed(113)
  res2 <- ss_aipe_icc_sensitivity(true_rho = 0.50, specified_N = 50,
                                  k = 3, width = 0.30, G = 50,
                                  type = "ICC(2,1)", print_iter = FALSE)
  expect_identical(attr(res2, "icc_type"), "ICC(2,1)")
  expect_gt(res2$value[res2$term == "mean_icc"], 0.10)

  set.seed(113)
  res3 <- ss_aipe_icc_sensitivity(true_rho = 0.50, specified_N = 50,
                                  k = 3, width = 0.30, G = 50,
                                  type = "ICC(3,1)", print_iter = FALSE)
  expect_identical(attr(res3, "icc_type"), "ICC(3,1)")
  expect_gt(res3$value[res3$term == "mean_icc"], 0.10)
})
