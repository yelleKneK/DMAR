test_that("ci_smd() returns a 3-row tidy data frame with the expected terms", {
  result <- ci_smd(ncp = 2.83, n_1 = 60, n_2 = 70, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, c("lower_limit", "smd", "upper_limit"))
})

test_that("ci_smd() lower_limit < smd < upper_limit for a positive ncp", {
  result <- ci_smd(ncp = 2.83, n_1 = 60, n_2 = 70, conf_level = 0.95)
  ll  <- result$value[result$term == "lower_limit"]
  est <- result$value[result$term == "smd"]
  ul  <- result$value[result$term == "upper_limit"]
  expect_lt(ll, est)
  expect_lt(est, ul)
})

test_that("ci_smd() accepts smd in place of ncp and gives consistent results", {
  smd_value <- 2.83 * sqrt((60 + 70) / (60 * 70))
  r_from_smd <- ci_smd(smd = smd_value, n_1 = 60, n_2 = 70, conf_level = 0.95)
  r_from_ncp <- ci_smd(ncp = 2.83,      n_1 = 60, n_2 = 70, conf_level = 0.95)
  expect_equal(r_from_smd$value, r_from_ncp$value, tolerance = 1e-6)
})

test_that("ci_smd() asymmetric alphas produce a wider lower / narrower upper CI", {
  sym  <- ci_smd(ncp = 2.83, n_1 = 60, n_2 = 70, conf_level = NULL, alpha_lower = 0.025, alpha_upper = 0.025)
  asym <- ci_smd(ncp = 2.83, n_1 = 60, n_2 = 70, conf_level = NULL, alpha_lower = 0.04,  alpha_upper = 0.01)
  expect_gt(asym$value[asym$term == "lower_limit"], sym$value[sym$term == "lower_limit"])
  expect_gt(asym$value[asym$term == "upper_limit"], sym$value[sym$term == "upper_limit"])
})

test_that("ci_smd() rejects mixing conf_level and alpha bounds explicitly", {
  expect_error(
    ci_smd(ncp = 2.83, n_1 = 60, n_2 = 70,
           conf_level = 0.95, alpha_lower = 0.025, alpha_upper = 0.025),
    "cannot mix them"
  )
})

test_that("ci_smd() accepts alpha bounds without requiring conf_level = NULL explicitly", {
  # With the default conf_level = .95 in place but alphas explicitly supplied,
  # the function should silently treat the alphas as the source of truth
  # rather than erroring on "mix" -- only an *explicit* conf_level mixed with
  # alphas is rejected.
  expect_silent(
    res <- ci_smd(ncp = 2.83, n_1 = 60, n_2 = 70,
                  alpha_lower = 0.025, alpha_upper = 0.025)
  )
  expect_equal(res$term, c("lower_limit", "smd", "upper_limit"))
})

test_that("ci_smd_c() is not a paired-design interval (no correlation argument)", {
  # The ci_smd help page states that a paired-design SMD interval is
  # not currently provided; ci_smd_c() is the interval for Glass's
  # estimator (control group SD, two independent groups) and takes no
  # correlation between paired measurements. This pins the fact the
  # page now states.
  expect_false(any(c("rho", "correlation") %in% names(formals(ci_smd_c))))
  expect_true(all(c("n_C", "n_E") %in% names(formals(ci_smd_c))))
})
