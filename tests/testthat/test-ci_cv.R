test_that("ci_cv() returns the expected 4-row tidy data frame", {
  result <- ci_cv(cv = 0.2, N = 30, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
  expect_setequal(result$term, c("ci_lower", "ci_upper", "cv", "cv_unbiased"))
  # The rows in order, under the names of the argument and of cv()'s row.
  expect_identical(result$term,
                   c("ci_lower", "cv", "cv_unbiased", "ci_upper"))
  expect_false(any(c("c_of_v", "c_of_v_unbiased") %in% result$term))
})

test_that("ci_cv() ci_lower < cv < ci_upper for a typical positive cv", {
  result <- ci_cv(cv = 0.2, N = 30, conf_level = 0.95)
  ll  <- result$value[result$term == "ci_lower"]
  est <- result$value[result$term == "cv"]
  ul  <- result$value[result$term == "ci_upper"]
  expect_lt(ll, est)
  expect_lt(est, ul)
})

test_that("ci_cv() accepts (mean, sd) in place of cv and gives consistent limits", {
  m <- 10; s <- 2; n <- 30
  r1 <- ci_cv(cv = s / m, N = n, conf_level = 0.95)
  r2 <- ci_cv(mean = m, sd = s, N = n, conf_level = 0.95)
  expect_equal(
    r1$value[r1$term %in% c("ci_lower", "ci_upper")],
    r2$value[r2$term %in% c("ci_lower", "ci_upper")],
    tolerance = 1e-8
  )
})

limits_cv <- function(x) x$value[x$term %in% c("ci_lower", "ci_upper")]

test_that("ci_cv() gives the interval of Kelley (2007, p. 758)", {
  # The coefficient of variation of glucose metabolism in 35 patients was
  # .271, and the article prints the 95 percent interval as [.216, .363].
  lims <- limits_cv(ci_cv(cv = .271, N = 35))
  expect_identical(round(lims, 3), c(.216, .363))
  # The package's own values, to more figures than the article prints.
  expect_equal(lims, c(0.2164671, 0.3633569), tolerance = 1e-6)
  # The two point estimates between the limits: the coefficient as given,
  # and (1 + 1/(4 N)) times it, here .271 times 141/140, which is
  # 0.2729357.
  out <- ci_cv(cv = .271, N = 35)
  expect_equal(out$value[out$term == "cv"], .271, tolerance = 1e-12)
  expect_equal(out$value[out$term == "cv_unbiased"], .271 * 141 / 140,
               tolerance = 1e-12)
})

test_that("ci_cv() agrees with limits computed to 30 digits", {
  # Limits by quadrature of the noncentral t distribution in 30-digit
  # arithmetic (mpmath), computed apart from the package: equal tails,
  # small tails with a noncentrality the range of pt() does not reach
  # exactly, and unequal tails.
  expect_equal(limits_cv(ci_cv(cv = .2, N = 30, conf_level = .95)),
               c(0.15811520714959514, 0.27270129792523299), tolerance = 1e-9)
  expect_equal(limits_cv(ci_cv(cv = .15, N = 15, conf_level = NULL,
                               alpha_lower = .005, alpha_upper = .005)),
               c(0.099680411185073852, 0.2844686771994545), tolerance = 1e-9)
  expect_equal(limits_cv(ci_cv(cv = .3, N = 20, conf_level = NULL,
                               alpha_lower = .01, alpha_upper = .04)),
               c(0.21293711196397681, 0.43508015068140117), tolerance = 1e-9)
})

test_that("ci_cv() gives each one-sided bound from the right tail", {
  # A lower bound takes the upper limit of the noncentrality parameter and
  # an upper bound the lower limit. Each bound is the matching limit of
  # the 90 percent two-sided interval, and the root of pt() in the
  # noncentrality parameter (sqrt(20) / .3 is 14.9, where pt() is exact).
  lower <- ci_cv(cv = .3, N = 20, conf_level = NULL, alpha_lower = .05, alpha_upper = 0)
  upper <- ci_cv(cv = .3, N = 20, conf_level = NULL, alpha_lower = 0, alpha_upper = .05)
  expect_equal(lower$value[lower$term == "ci_lower"], 0.234279967, tolerance = 1e-8)
  expect_equal(upper$value[upper$term == "ci_upper"], 0.424624431, tolerance = 1e-8)
  expect_identical(lower$value[lower$term == "ci_upper"], Inf)
  expect_identical(upper$value[upper$term == "ci_lower"], 0)
  two_sided <- limits_cv(ci_cv(cv = .3, N = 20, conf_level = .90))
  expect_equal(c(lower$value[1], upper$value[4]), two_sided, tolerance = 1e-10)
  t_obs <- sqrt(20) / .3
  by_pt <- sqrt(20) / c(
    stats::uniroot(function(d) stats::pt(t_obs, 19, d) - .05,
                   c(t_obs, t_obs + 20), tol = 1e-13)$root,
    stats::uniroot(function(d) stats::pt(t_obs, 19, d) - .95,
                   c(0, t_obs), tol = 1e-13)$root)
  expect_equal(c(lower$value[1], upper$value[4]), by_pt, tolerance = 1e-8)
  # The tail columns of a bound are the requested areas, 0 on its open side.
  expect_identical(lower$prob_less, c(.05, NA, NA, 1))
  expect_identical(lower$prob_greater, c(.95, NA, NA, 0))
  expect_identical(upper$prob_less, c(0, NA, NA, .95))
  expect_identical(upper$prob_greater, c(1, NA, NA, .05))
})

test_that("ci_cv() reports the requested tail areas, and its three routes agree with unequal tails", {
  set.seed(113)
  x <- rnorm(15, 5, 1)
  by_data <- ci_cv(data = x, alpha_lower = .01, alpha_upper = .04)
  by_stats <- ci_cv(mean = mean(x), sd = sd(x), N = 15, alpha_lower = .01,
                    alpha_upper = .04)
  by_cv <- ci_cv(cv = sd(x) / mean(x), N = 15, alpha_lower = .01,
                 alpha_upper = .04)
  expect_identical(by_data$prob_less, c(.01, NA, NA, .96))
  expect_identical(by_data$prob_greater, c(.99, NA, NA, .04))
  expect_identical(as.data.frame(by_stats), as.data.frame(by_data))
  expect_identical(as.data.frame(by_cv), as.data.frame(by_data))
  expect_equal(attr(by_data, "conf_level"), .95)
  expect_identical(attr(by_data, "tail_areas"), c(lower = .01, upper = .04))
  # More of the area above the interval moves both limits down.
  equal_tails <- limits_cv(ci_cv(data = x))
  expect_true(all(limits_cv(by_data) < equal_tails))
})

test_that("ci_cv() gives an infinite upper limit without a warning and with the requested tail areas", {
  # sqrt(6) / 1.2 is 2.04, below the .975 quantile of the central t
  # distribution on 5 degrees of freedom, 2.57, so the lower limit of the
  # noncentrality parameter is not positive and the upper limit of the
  # interval is infinite. That interval has the stated coverage.
  expect_lt(sqrt(6) / 1.2, stats::qt(.975, 5))
  out <- expect_no_condition(ci_cv(cv = 1.2, N = 6))
  expect_identical(out$value[out$term == "ci_upper"], Inf)
  expect_equal(out$value[out$term == "ci_lower"], 0.5714759, tolerance = 1e-6)
  expect_equal(out$prob_greater[out$term == "ci_upper"], .025, tolerance = 1e-15)
  expect_equal(out$prob_less[out$term == "ci_upper"], .975, tolerance = 1e-15)
  expect_equal(out$prob_less[out$term == "ci_lower"], .025, tolerance = 1e-15)
  expect_equal(attr(out, "conf_level"), .95)
  # Just above the quantile the limit is finite.
  finite <- expect_no_condition(ci_cv(cv = sqrt(6) / 2.6, N = 6))
  expect_true(is.finite(finite$value[finite$term == "ci_upper"]))
  # With unequal tails the rule uses the upper tail area.
  cut <- stats::qt(1 - .04, 5)
  expect_identical(ci_cv(cv = sqrt(6) / (cut - .01), N = 6, alpha_lower = .01,
                         alpha_upper = .04)$value[4], Inf)
  expect_true(is.finite(ci_cv(cv = sqrt(6) / (cut + .01), N = 6, alpha_lower = .01,
                              alpha_upper = .04)$value[4]))
})

test_that("ci_cv() covers the population value at the stated rate when many upper limits are infinite", {
  skip_on_cran()
  # Normal samples of 6 with a coefficient of variation of .8: about 31
  # percent of the 95 percent intervals have an infinite upper limit, and
  # the coverage over all of them is the stated level. With 2,000 samples
  # the standard error of the coverage is .0049.
  set.seed(113)
  res <- t(vapply(seq_len(2000), function(i) {
    x <- rnorm(6, 1, .8)
    if (mean(x) <= 0) return(c(NA_real_, NA_real_))
    limits_cv(ci_cv(data = x))
  }, numeric(2)))
  res <- res[!is.na(res[, 1]), ]
  expect_gt(nrow(res), 1990)
  share_infinite <- mean(is.infinite(res[, 2]))
  expect_gt(share_infinite, .27)
  expect_lt(share_infinite, .36)
  coverage <- mean(res[, 1] <= .8 & .8 <= res[, 2])
  expect_gt(coverage, .95 - 4 * .0049)
  expect_lt(coverage, .95 + 4 * .0049)
})

test_that("ci_cv() names the argument that is missing, malformed, or given twice", {
  refusal <- function(expr) {
    e <- tryCatch(expr, error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    e
  }
  # The sample size is what is missing when the coefficient is given alone.
  e <- refusal(ci_cv(cv = .2))
  expect_identical(e$arg, "N")
  expect_match(conditionMessage(e), "'N' is required with 'cv'", fixed = TRUE)
  expect_identical(refusal(ci_cv(mean = 5, sd = 1))$arg, "N")
  expect_identical(refusal(ci_cv(mean = 5, N = 20))$arg, "sd")
  expect_identical(refusal(ci_cv(sd = 1, N = 20))$arg, "mean")
  expect_identical(refusal(ci_cv())$arg, c("data", "cv", "mean", "sd"))
  # The scores and the coefficient are two routes to one interval.
  e <- refusal(ci_cv(data = 1:10, cv = .2))
  expect_identical(e$arg, c("data", "cv"))
  expect_identical(conditionMessage(e),
                   "Give either 'data', or 'cv' with 'N', not both.")
  expect_identical(refusal(ci_cv(data = 1:10, cv = .2, N = 10))$arg, c("data", "cv"))
  expect_identical(refusal(ci_cv(data = 1:10, N = 10))$arg, c("data", "N"))
  expect_identical(refusal(ci_cv(cv = .2, N = 10, sd = 1))$arg, c("cv", "sd"))
  # A standard deviation or a mean that is not a number is refused by name.
  e <- refusal(ci_cv(mean = 5, sd = "a", N = 30))
  expect_identical(e$arg, "sd")
  expect_match(conditionMessage(e), "'sd' must be a single finite number", fixed = TRUE)
  expect_identical(refusal(ci_cv(mean = "a", sd = 1, N = 30))$arg, "mean")
  expect_identical(refusal(ci_cv(mean = c(5, 6), sd = 1, N = 30))$arg, "mean")
  # One message for a mean that is not positive, from 'mean' or from the
  # scores, and the same for a standard deviation.
  by_mean <- refusal(ci_cv(mean = -2.5, sd = 1.2, N = 4))
  by_data <- refusal(ci_cv(data = c(-1, -2, -3, -4)))
  expect_identical(conditionMessage(by_mean), conditionMessage(by_data))
  expect_identical(conditionMessage(by_data), paste0(
    "The sample mean is -2.5, and it must be positive: the confidence ",
    "interval for the coefficient of variation requires a positive mean ",
    "and a positive standard deviation."))
  expect_identical(by_mean$arg, "mean")
  expect_identical(by_data$arg, "data")
  expect_no_match(conditionMessage(by_mean), "absolute value")
  expect_identical(refusal(ci_cv(mean = 0, sd = 1, N = 10))$arg, "mean")
  expect_identical(refusal(ci_cv(data = c(-2, -1, 1, 2)))$arg, "data")
  expect_identical(refusal(ci_cv(mean = 10, sd = 0, N = 20))$arg, "sd")
  expect_identical(refusal(ci_cv(data = c(5, 5, 5, 5)))$arg, "data")
  expect_identical(refusal(ci_cv(cv = -.2, N = 20))$arg, "cv")
})
