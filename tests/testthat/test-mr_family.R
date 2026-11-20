## mr_smd() and mr_cv(): minimum risk point estimation by a purely
## sequential procedure, in which one observation (one per group for
## mr_smd()) is added after the pilot sample until the stopping rule is met
## (Chattopadhyay & Kelley, 2016, p. 633).

# The 35 observations of the worked example of Chattopadhyay and Kelley
# (2016, pp. 637-638), with A = 400000 and a sampling cost of 75; the first
# 18 are the pilot sample.
mr_books <- c(36, 53, 19, 11, 10, 24, 14, 65, 18, 48, 25, 35, 13, 18, 3, 41,
              5, 3, 44, 26, 13, 39, 2, 3, 26, 22, 8, 15, 12, 22, 5, 21, 23,
              40, 18)
mr_term <- function(x, t) x$value[x$term == t]
# The refusal a call gives, or NULL when it gives none.
mr_refusal <- function(expr) tryCatch({ expr; NULL }, error = function(e) e)
expect_refused <- function(expr, arg, class = "dmar_argument_invalid") {
  e <- mr_refusal(expr)
  expect_s3_class(e, class)
  expect_s3_class(e, "dmar_refusal")
  expect_identical(e$arg, arg)
}

test_that("mr_smd() returns a tidy data.frame with numeric value column", {
  res <- mr_smd(structural_cost = 0.0001, epsilon = 0.05, d = 0.5,
                n_per_group = 30, sampling_cost = 0.001)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("risk", "n_1", "n_2", "d") %in% res$term))
  expect_true(is.numeric(res$value))
  # is_satisfied is non-numeric metadata and lives on an attribute,
  # not as a row with a stringly-typed value.
  expect_false("is_satisfied" %in% res$term)
  expect_true(is.logical(attr(res, "is_satisfied")))
  expect_true(is.numeric(attr(res, "criterion")))
})

test_that("mr_smd() risk equals its closed form at a pinned design", {
  # A = structural_cost / epsilon^2 = 0.0001 / 0.05^2 = 0.04.
  # Risk = A * (2/n + d^2 / (4n)) + sampling_cost * 2n
  #      = 0.04 * (2/30 + 0.25/120) + 0.001 * 60 = 0.06275.
  res <- mr_smd(structural_cost = 0.0001, epsilon = 0.05, d = 0.5,
                n_per_group = 30, sampling_cost = 0.001)
  A <- 0.0001 / 0.05^2
  risk_hand <- A * (2 / 30 + 0.5^2 / (4 * 30)) + 0.001 * 60
  expect_equal(res$value[res$term == "risk"], risk_hand,
               tolerance = 1e-12)
  expect_equal(res$value[res$term == "risk"], 0.06275, tolerance = 1e-10)
})

test_that("mr_smd() accepts A directly and skips epsilon/structural_cost", {
  res <- mr_smd(A = 1, d = 0.5, n_per_group = 30, sampling_cost = 0.001)
  expect_s3_class(res, "data.frame")
  expect_true("risk" %in% res$term)
})

test_that("mr_smd() raises an error when A and epsilon are both supplied", {
  expect_error(
    mr_smd(A = 1, epsilon = 0.05, d = 0.5, n_per_group = 30, sampling_cost = 0.001),
    "specified 'A'"
  )
})

test_that("mr_cv() pilot mode returns a small pilot sample size", {
  set.seed(113)
  x <- rnorm(40, mean = 10, sd = 2)
  res <- mr_cv(data = x, A = 1, sampling_cost = 0.001, pilot = TRUE, m0 = 4)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("pilot_N" %in% res$term)
  expect_gt(res$value[res$term == "pilot_N"], 0)
})


test_that("mr_cv() non-pilot mode returns numeric value column with attributes", {
  set.seed(113)
  x <- rlnorm(60, meanlog = 1, sdlog = 0.3)
  res <- mr_cv(data = x, A = 1, sampling_cost = 0.001, gamma = 0.49,
               pilot = FALSE, m0 = 60, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_true(is.numeric(res$value))
  expect_false("is_satisfied" %in% res$term)
  expect_true(is.logical(attr(res, "is_satisfied")))
  expect_true(is.numeric(attr(res, "criterion")))
  expect_false(attr(res, "variance_unusable"))
})


test_that("mr_cv() verbose mode keeps numeric value column with extra term rows", {
  set.seed(113)
  x <- rlnorm(60)
  res <- mr_cv(data = x, A = 1, sampling_cost = 0.001, gamma = 0.49,
               pilot = FALSE, m0 = 60, verbose = TRUE)
  expect_true(is.numeric(res$value))
  expect_true(all(c("risk", "N", "cv", "v2_n", "mu_3_n", "mu_4_n",
                    "criterion") %in% res$term))
})


test_that("mr_cv() differs from the worked example of Chattopadhyay and Kelley (2016) by the correction of Equation 22 and by nothing else", {
  # The article estimates xi^2 by its Equation 22; mr_cv() uses Equation 32
  # of Kelley, Darku, and Chattopadhyay (2018). On any sample the second
  # minus the first is D = (k^2 / 4)(1 - k^2)(mu4_hat / s^4 - 3), with k
  # the sample coefficient of variation. Taking A D / n out of the
  # function's risk and (A / c) D out of its criterion must therefore give
  # the numbers the article prints (pp. 637-638): risks of 5964.345,
  # 6224.861, and 4891.284 at 18, 19, and 35 observations, coefficients of
  # variation of 0.7391157, 0.7113385, and 0.7013904, and a stopping rule
  # first met at 35. The function's own risk and criterion enter each
  # expectation, so a variance that is not exactly Equation 32 fails here.
  A <- 400000; cc <- 75
  risk_2016 <- crit_2016 <- cv_n <- numeric(0)
  for (n in 18:35) {
    x <- mr_books[seq_len(n)]
    res <- mr_cv(data = x, A = A, sampling_cost = cc, verbose = TRUE)
    k <- mr_term(res, "cv")
    D <- (k^2 / 4) * (1 - k^2) * (mr_term(res, "mu_4_n") / var(x)^2 - 3)
    risk_2016 <- c(risk_2016, mr_term(res, "risk") - (A / n) * D)
    crit_2016 <- c(crit_2016, attr(res, "criterion") - (A / cc) * D)
    cv_n <- c(cv_n, k)
    expect_equal(mr_term(res, "criterion"), attr(res, "criterion"))
  }
  at <- match(c(18, 19, 35), 18:35)
  expect_identical(round(risk_2016[at], 3), c(5964.345, 6224.861, 4891.284))
  expect_identical(round(cv_n[at], 7), c(0.7391157, 0.7113385, 0.7013904))
  expect_identical((18:35)[which((18:35)^2 >= crit_2016)[1L]], 35L)
  # At 35 the rule of Equation 22 is met by 3.79, and at 34 it is short by
  # 123.02.
  expect_equal(35^2 - crit_2016[18], 3.79, tolerance = 2e-3)
  expect_equal(34^2 - crit_2016[17], -123.02, tolerance = 1e-4)
})

test_that("mr_cv() does not stop on the article's 35 observations, and its values at 35 are pinned", {
  A <- 400000; cc <- 75
  met <- vapply(18:35, function(n)
    attr(mr_cv(data = mr_books[seq_len(n)], A = A, sampling_cost = cc),
         "is_satisfied"), logical(1))
  expect_false(any(met))
  # Regression values of this function (not published ones): the estimate
  # of Equation 32, the risk, and the criterion at n = 35.
  res <- mr_cv(data = mr_books, A = A, sampling_cost = cc, verbose = TRUE)
  expect_equal(mr_term(res, "v2_n"), 0.20582264894, tolerance = 1e-9)
  expect_equal(mr_term(res, "risk"), 4977.258845025, tolerance = 1e-9)
  expect_equal(mr_term(res, "criterion"), 1261.331623081, tolerance = 1e-9)
  expect_equal(attr(res, "criterion"), 1261.331623081, tolerance = 1e-9)
  expect_equal(mr_term(res, "N"), 35)
  expect_equal(mr_term(res, "cv"), sd(mr_books) / mean(mr_books),
               tolerance = 1e-12)
  expect_false(attr(res, "is_satisfied"))
  expect_false(attr(res, "variance_unusable"))
  # The risk and the criterion are those of the article's Equations 14 and
  # 23 at the function's own estimate.
  v2 <- mr_term(res, "v2_n")
  expect_equal(mr_term(res, "risk"), A * v2 / 35 + cc * 35, tolerance = 1e-12)
  expect_equal(attr(res, "criterion"), (A / cc) * (v2 + 35^(-2 * .49)),
               tolerance = 1e-12)
})

test_that("mr_cv() estimates the variance of the sample coefficient of variation in a large sample", {
  # For gamma scores with a coefficient of variation k the skewness is 2k
  # and the excess kurtosis 6k^2, so xi^2 = k^2/2 + k^4 - k^3 (2k) +
  # k^2 (6k^2) / 4, which is 0.5248 at k = .8. (Equation 22 of the 2016
  # article estimates 0.3036 there.)
  k <- .8
  xi2 <- k^2 / 2 + k^4 - k^3 * (2 * k) + k^2 * (6 * k^2) / 4
  expect_equal(xi2, 0.5248, tolerance = 1e-12)
  set.seed(113)
  z <- rgamma(2e6, shape = 1 / k^2)
  res <- mr_cv(data = z, A = 1, sampling_cost = 1, verbose = TRUE)
  expect_lt(abs(mr_term(res, "v2_n") - xi2), 0.02)
  # An oracle that does not retype Equation 32: the mean squared influence
  # value of s / M, k ((x - M)^2 - s^2) / (2 s^2) - k (x - M) / M, at the
  # sample moments.
  m <- mean(z); s2 <- var(z); kn <- sqrt(s2) / m
  infl <- kn * (((z - m)^2 - s2) / (2 * s2) - (z - m) / m)
  expect_equal(mr_term(res, "v2_n"), mean(infl^2), tolerance = 1e-4)
})

test_that("mr_cv() gives the article's pilot sample size by both routes", {
  by_A <- mr_cv(pilot = TRUE, A = 400000, sampling_cost = 75, gamma = .49)
  by_parts <- mr_cv(pilot = TRUE, structural_cost = 1000, epsilon = .05,
                    sampling_cost = 75, gamma = .49)
  expect_identical(by_A$term, "pilot_N")
  expect_equal(by_A$value, 18)
  expect_equal(by_parts$value, 18)
  expect_equal(by_A$value, max(4, ceiling((400000 / 75)^(1 / (2 + 2 * .49)))))
  # m0 is the floor of the pilot sample size.
  expect_equal(mr_cv(pilot = TRUE, A = 400000, sampling_cost = 75,
                     m0 = 30)$value, 30)
})

test_that("mr_cv() reports the rule as not met when the variance estimate is zero, negative, or not finite", {
  unusable <- function(res, n) {
    expect_false(attr(res, "is_satisfied"))
    expect_true(attr(res, "variance_unusable"))
    expect_identical(attr(res, "criterion"), NA_real_)
    expect_identical(mr_term(res, "risk"), NA_real_)
    expect_identical(mr_term(res, "criterion"), NA_real_)
    expect_equal(mr_term(res, "N"), n)
    expect_true(is.numeric(res$value))
  }
  # Eight scores at the pilot sample size of A = 485 and c = 1: a negative
  # estimate, which is the value the v2_n row holds.
  x8 <- c(2.3, 0.61, 1.6, 0.54, 0.51, 0.41, 0.67, 0.37)
  expect_equal(mr_cv(pilot = TRUE, A = 485, sampling_cost = 1)$value, 8)
  res <- mr_cv(data = x8, A = 485, sampling_cost = 1, verbose = TRUE)
  unusable(res, 8)
  expect_equal(mr_term(res, "v2_n"), -0.0249923, tolerance = 1e-5)
  expect_identical(mr_term(mr_cv(data = x8, A = 485, sampling_cost = 1),
                           "risk"), NA_real_)
  # Here n^-3 in place of the estimate would meet the rule: the criterion
  # would be 400 (8^-3 + 8^-.98) = 52.9, below 64.
  expect_lt(400 * (8^-3 + 8^-.98), 64)
  res <- mr_cv(data = c(1, 1, 1, 1, 1, 1, 3, 3), A = 400, sampling_cost = 1,
               verbose = TRUE)
  unusable(res, 8)
  expect_lt(mr_term(res, "v2_n"), 0)
  # Scores that do not vary: the estimate is 0, and again n^-3 would meet
  # the rule (314.8 against 324).
  expect_lt((400000 / 75) * (18^-3 + 18^-.98), 18^2)
  res <- mr_cv(data = rep(25, 18), A = 400000, sampling_cost = 75,
               verbose = TRUE)
  unusable(res, 18)
  expect_identical(mr_term(res, "v2_n"), 0)
  expect_identical(mr_term(res, "cv"), 0)
  # An estimate that is not finite is unusable as well. Scores a user would
  # give do not produce one, so the rule is given the value itself.
  for (v2 in c(NaN, Inf)) {
    expect_identical(DMAR:::.mr_cv_rule(v2, 10, 400, 1, .49),
                     list(unusable = TRUE, risk = NA_real_,
                          criterion = NA_real_, is_satisfied = FALSE))
  }
})

test_that("mr_cv() uses n^-3 in place of a positive estimate below it, and reports the estimate itself", {
  # Ten scores with a coefficient of variation of .018: the estimate is
  # 9.30e-5, below 10^-3, so the risk is A 10^-3 / 10 + 75 (10) = 790 and
  # the criterion (A / c)(10^-3 + 10^-.98) = 563.80. Without the bound the
  # risk would be 753.72.
  x <- c(98, 99, 100, 101, 102, 100, 99, 101, 103, 97)
  res <- mr_cv(data = x, A = 400000, sampling_cost = 75, verbose = TRUE)
  expect_equal(mr_term(res, "v2_n"), 9.30278e-5, tolerance = 1e-5)
  expect_lt(mr_term(res, "v2_n"), 10^-3)
  expect_equal(mr_term(res, "risk"), 790, tolerance = 1e-12)
  expect_equal(mr_term(res, "criterion"),
               (400000 / 75) * (10^-3 + 10^(-.98)), tolerance = 1e-12)
  expect_equal(attr(res, "criterion"), 563.80, tolerance = 1e-5)
  expect_false(attr(res, "is_satisfied"))
  expect_false(attr(res, "variance_unusable"))
})

test_that("mr_cv() computes the central moments from centered scores", {
  set.seed(113)
  z <- rnorm(30)
  mu4 <- function(x) mr_term(mr_cv(data = x, A = 400000, sampling_cost = 75,
                                   verbose = TRUE), "mu_4_n")
  mu3 <- function(x) mr_term(mr_cv(data = x, A = 400000, sampling_cost = 75,
                                   verbose = TRUE), "mu_3_n")
  # The unbiased estimates from the centered scores, written out.
  d <- z - mean(z); n <- 30
  mu4_hand <- ((n^2 - 2 * n + 3) * sum(d^4) - 3 * (2 * n - 3) / n * sum(d^2)^2) /
    ((n - 1) * (n - 2) * (n - 3))
  mu3_hand <- n * sum(d^3) / ((n - 1) * (n - 2))
  expect_equal(mu4(z), mu4_hand, tolerance = 1e-12)
  expect_equal(mu3(z), mu3_hand, tolerance = 1e-12)
  # The same in the raw power sums of the published form of the estimator
  # (Chattopadhyay & Kelley, 2016, Equation 21), which is accurate when the
  # mean is not large beside the standard deviation.
  mu4_raw <- (n^2 * sum(d^4) - (2 * n - 3) * sum(z^4) +
                (8 * n - 12) * mean(z) * sum(z^3) -
                (6 * n - 9) / n * sum(z^2)^2) / ((n - 1) * (n - 2) * (n - 3))
  expect_equal(mu4(z), mu4_raw, tolerance = 1e-10)
  # A central moment does not depend on where the scores are centered.
  expect_equal(mu4(z + 1e5), mu4(z), tolerance = 1e-6)
  expect_equal(mu4(z + 1e7), mu4(z), tolerance = 1e-6)
  expect_equal(mu3(z + 1e7), mu3(z), tolerance = 1e-6)
})

test_that("mr_cv() refuses each argument it cannot use, by name", {
  b <- mr_books[1:18]
  ok <- function(...) mr_cv(data = b, A = 400000, sampling_cost = 75, ...)
  # data: a numeric vector of at least 4 finite observations with a
  # positive mean. Centered scores are refused whether their mean comes
  # out as exactly 0 (b - mean(b)) or as a positive rounding error: the
  # mean of c(-1, -1, 1, 1 + 2^-50) is exactly 2^-52, and that of
  # b / 10 - mean(b / 10) is near 2e-16.
  tiny_mean <- c(-1, -1, 1, 1 + 2^-50)
  expect_identical(mean(tiny_mean), 2^-52)
  for (bad in list(data.frame(b), matrix(b, ncol = 1), as.character(b),
                   c(b, NA), c(b, Inf), b[1:3], numeric(0), -b,
                   b - mean(b), b / 10 - mean(b / 10), tiny_mean)) {
    expect_refused(mr_cv(data = bad, A = 400000, sampling_cost = 75), "data")
  }
  expect_refused(mr_cv(A = 400000, sampling_cost = 75), "data")
  expect_match(conditionMessage(mr_refusal(
    mr_cv(data = -b, A = 400000, sampling_cost = 75))), "positive mean")
  expect_match(conditionMessage(mr_refusal(
    mr_cv(data = tiny_mean, A = 400000, sampling_cost = 75))),
    "positive mean; the mean of the 4 observations is 2.22e-16, which is zero to rounding",
    fixed = TRUE)
  # A positive mean that is small beside the largest score is accepted:
  # the mean of these 500 lognormal scores is 0.045 of their maximum.
  set.seed(113)
  skewed <- rlnorm(500, meanlog = 0, sdlog = 1.5)
  expect_lt(mean(skewed) / max(skewed), 0.05)
  expect_identical(mr_term(mr_cv(data = skewed, A = 400000,
                                 sampling_cost = 75), "N"), 500)
  expect_match(conditionMessage(mr_refusal(
    mr_cv(data = b[1:3], A = 400000, sampling_cost = 75))),
    "at least 4 finite observations; 3 were given", fixed = TRUE)
  expect_match(conditionMessage(mr_refusal(
    mr_cv(data = c(b, NA), A = 400000, sampling_cost = 75))),
    "1 of the 19 values is missing or infinite", fixed = TRUE)
  expect_match(conditionMessage(mr_refusal(
    mr_cv(data = c(b, NA, Inf), A = 400000, sampling_cost = 75))),
    "2 of the 20 values are missing or infinite", fixed = TRUE)
  # gamma: strictly between 0 and 1/2, so .5 itself is refused.
  for (bad in list(.5, 0, -1, .7, NA, "0.49", c(.3, .4), NULL)) {
    expect_refused(ok(gamma = bad), "gamma")
    expect_refused(mr_cv(pilot = TRUE, A = 400000, sampling_cost = 75,
                         gamma = bad), "gamma")
  }
  expect_no_error(ok(gamma = .499))
  # sampling_cost, A, structural_cost, epsilon: single positive finite
  # numbers, with the data and for the pilot sample size alike.
  for (bad in list(0, -75, NA, Inf, "75", c(75, 80), NULL)) {
    expect_refused(mr_cv(data = b, A = 400000, sampling_cost = bad),
                   "sampling_cost")
    expect_refused(mr_cv(pilot = TRUE, A = 400000, sampling_cost = bad),
                   "sampling_cost")
  }
  expect_refused(mr_cv(data = b, A = 400000), "sampling_cost")
  expect_refused(mr_cv(pilot = TRUE, A = 400000), "sampling_cost")
  for (bad in list(0, -1, NA, Inf, "a", c(1, 2), NULL)) {
    expect_refused(mr_cv(data = b, A = bad, sampling_cost = 75), "A")
    expect_refused(mr_cv(pilot = TRUE, A = bad, sampling_cost = 75), "A")
    expect_refused(mr_cv(pilot = TRUE, structural_cost = bad, epsilon = .05,
                         sampling_cost = 75), "structural_cost")
    expect_refused(mr_cv(pilot = TRUE, structural_cost = 1000, epsilon = bad,
                         sampling_cost = 75), "epsilon")
  }
  expect_refused(mr_cv(pilot = TRUE, A = 400000, structural_cost = 1000,
                       sampling_cost = 75), c("A", "structural_cost"))
  expect_refused(mr_cv(pilot = TRUE, A = 400000, epsilon = .05,
                       sampling_cost = 75), c("A", "epsilon"))
  expect_refused(mr_cv(pilot = TRUE, sampling_cost = 75),
                 c("A", "structural_cost"))
  expect_refused(mr_cv(pilot = TRUE, structural_cost = 1000,
                       sampling_cost = 75), c("A", "epsilon"))
  # m0: a whole number of at least 4.
  for (bad in list(4.5, 3, 2, NA, "4", c(4, 5))) {
    expect_refused(mr_cv(pilot = TRUE, A = 400000, sampling_cost = 75,
                         m0 = bad), "m0", class = "dmar_count_invalid")
  }
  # pilot and verbose: TRUE or FALSE.
  for (bad in list("yes", 1, 0, NA, c(TRUE, FALSE), NULL)) {
    expect_refused(mr_cv(pilot = bad, A = 400000, sampling_cost = 75),
                   "pilot")
    expect_refused(ok(verbose = bad), "verbose")
  }
})

test_that("mr_smd() refuses each argument it cannot use, by name", {
  ok <- function(...) mr_smd(A = 10000, sampling_cost = 2.4, ...)
  # d: a single finite number.
  for (bad in list(c(.2, .5), NA, Inf, "a", NULL)) {
    expect_refused(ok(d = bad, n_per_group = 13), "d")
  }
  expect_refused(ok(n_per_group = 13), "d")
  expect_no_error(ok(d = -.5, n_per_group = 13))
  expect_no_error(ok(d = 0, n_per_group = 13))
  # n_per_group: a whole number of at least 2.
  for (bad in list(-5, 0, 1, 74.5, c(13, 75), NA, "13", NULL)) {
    expect_refused(ok(d = .5, n_per_group = bad), "n_per_group",
                   class = "dmar_count_invalid")
  }
  expect_refused(ok(d = .5), "n_per_group")
  # d and n_per_group are not needed, and not checked, for the pilot size.
  expect_equal(ok(pilot = TRUE)$value, 13)
  expect_equal(ok(pilot = TRUE, d = "a", n_per_group = -5)$value, 13)
  # gamma, the costs, m0, and pilot, as in mr_cv().
  for (bad in list(.5, 0, -1, NA, "0.49", c(.3, .4))) {
    expect_refused(ok(d = .5, n_per_group = 13, gamma = bad), "gamma")
    expect_refused(ok(pilot = TRUE, gamma = bad), "gamma")
  }
  for (bad in list(0, -2.4, NA, Inf, "a", c(1, 2))) {
    expect_refused(mr_smd(A = 10000, sampling_cost = bad, d = .5,
                          n_per_group = 13), "sampling_cost")
    expect_refused(mr_smd(A = bad, sampling_cost = 2.4, d = .5,
                          n_per_group = 13), "A")
    expect_refused(mr_smd(structural_cost = bad, epsilon = .1,
                          sampling_cost = 2.4, pilot = TRUE),
                   "structural_cost")
    expect_refused(mr_smd(structural_cost = 100, epsilon = bad,
                          sampling_cost = 2.4, pilot = TRUE), "epsilon")
  }
  expect_refused(mr_smd(A = 10000, d = .5, n_per_group = 13), "sampling_cost")
  for (bad in list(3, 4.5, NA)) {
    expect_refused(ok(pilot = TRUE, m0 = bad), "m0",
                   class = "dmar_count_invalid")
  }
  for (bad in list("yes", 1, NA)) {
    expect_refused(ok(pilot = bad, d = .5, n_per_group = 13), "pilot")
  }
})

test_that("mr_cv() prints its decision under the table in every state", {
  line <- function(x, ...) unname(DMAR:::.dmar_tbl_footer(x, ...)["stopping_rule"])
  not_met <- mr_cv(data = mr_books, A = 400000, sampling_cost = 75)
  expect_s3_class(not_met, c("dmar_mr", "dmar_tbl", "data.frame"),
                  exact = TRUE)
  expect_identical(line(not_met), paste0(
    "Stopping rule: not met (N^2 = 1225 is below the criterion, 1261.332); ",
    "collect another observation."))
  # With A = 380000 the rule is met on the same 35 observations (a
  # criterion of 1198.265 against 1225).
  met <- mr_cv(data = mr_books, A = 380000, sampling_cost = 75)
  expect_true(attr(met, "is_satisfied"))
  expect_identical(line(met), paste0(
    "Stopping rule: met (N^2 = 1225 is at least the criterion, 1198.265); ",
    "sampling can stop."))
  unusable <- mr_cv(data = c(2.3, 0.61, 1.6, 0.54, 0.51, 0.41, 0.67, 0.37),
                    A = 485, sampling_cost = 1)
  expect_identical(line(unusable), paste0(
    "The variance estimate is zero or negative at this sample size, so the ",
    "stopping rule is reported as not met; collect another observation."))
  # The console shows the line, and the criterion follows print()'s digits.
  expect_output(print(not_met), "Stopping rule: not met (N^2 = 1225 is below",
                fixed = TRUE)
  expect_output(print(not_met, digits = 1), "the criterion, 1261.3);",
                fixed = TRUE)
  expect_identical(line(not_met, digits = 1), paste0(
    "Stopping rule: not met (N^2 = 1225 is below the criterion, 1261.3); ",
    "collect another observation."))
  expect_output(print(unusable), "reported as not met", fixed = TRUE)
  # The verbose table carries the same line.
  expect_identical(line(mr_cv(data = mr_books, A = 400000, sampling_cost = 75,
                              verbose = TRUE)), line(not_met))
  # Nothing is added for a pilot sample size, or for a table that has lost
  # the decision or the row it compares.
  pilot <- mr_cv(pilot = TRUE, A = 400000, sampling_cost = 75)
  expect_false("stopping_rule" %in% names(DMAR:::.dmar_tbl_footer(pilot)))
  bare <- not_met
  attr(bare, "is_satisfied") <- NULL
  expect_false("stopping_rule" %in% names(DMAR:::.dmar_tbl_footer(bare)))
  expect_output(print(bare), "risk", fixed = TRUE)
  cut <- not_met[not_met$term == "cv", ]
  expect_false("stopping_rule" %in% names(DMAR:::.dmar_tbl_footer(cut)))
})

test_that("mr_smd() prints its decision under the table in both states", {
  line <- function(x) unname(DMAR:::.dmar_tbl_footer(x)["stopping_rule"])
  # A = 10000, c = 2.4, d = .5: the criterion is 78.539 at 13 per group and
  # 71.054 at 75.
  not_met <- mr_smd(d = .5, n_per_group = 13, A = 10000, sampling_cost = 2.4)
  expect_s3_class(not_met, c("dmar_mr", "dmar_tbl", "data.frame"),
                  exact = TRUE)
  expect_false(attr(not_met, "is_satisfied"))
  expect_identical(line(not_met), paste0(
    "Stopping rule: not met (the sample size per group, 13, is below the ",
    "criterion, 78.539); collect another observation in each group."))
  met <- mr_smd(d = .5, n_per_group = 75, A = 10000, sampling_cost = 2.4)
  expect_true(attr(met, "is_satisfied"))
  expect_identical(line(met), paste0(
    "Stopping rule: met (the sample size per group, 75, is at least the ",
    "criterion, 71.054); sampling can stop."))
  expect_output(print(met), "sampling can stop.", fixed = TRUE)
  pilot <- mr_smd(pilot = TRUE, A = 10000, sampling_cost = 2.4)
  expect_false("stopping_rule" %in% names(DMAR:::.dmar_tbl_footer(pilot)))
})

test_that("mr_smd() stopping rule and pilot size are those of Chattopadhyay and Kelley (2017)", {
  # Oracle: V_n^2 = 2 + d^2 / 4 (Eq. 13), the stopping rule
  # n >= sqrt(A / (2c)) (V_n + n^(-gamma)) (Eq. 15), and the pilot size
  # per group max(m0, ceiling((A / (2c))^(1 / (2 + 2 gamma)))) (Eq. 14).
  A <- 10000; cc <- 2.4; g <- .49
  for (d in c(0, .4, 1.2)) for (n in c(13, 40, 120)) {
    res <- mr_smd(A = A, d = d, n_per_group = n, sampling_cost = cc,
                  gamma = g)
    crit <- sqrt(A / (2 * cc)) * (sqrt(2 + d^2 / 4) + n^(-g))
    expect_equal(attr(res, "criterion"), crit, tolerance = 1e-12)
    expect_identical(attr(res, "is_satisfied"), n >= crit)
    expect_equal(res$value[res$term == "risk"],
                 A * (2 + d^2 / 4) / n + 2 * cc * n, tolerance = 1e-12)
  }
  pilot <- mr_smd(pilot = TRUE, A = A, sampling_cost = cc, gamma = g)
  expect_equal(pilot$value, max(4, ceiling((A / (2 * cc))^(1 / (2 + 2 * g)))))
  # The article's worked example: a pilot sample of 13 per group.
  expect_equal(pilot$value, 13)
})

test_that("a pilot sample size whose formula is a whole number is that number", {
  # (1e5)^0.4 and (1e6 / 10)^0.4 are exactly 100; R computes
  # 100.00000000000003, whose ceiling is 101.
  expect_gt((1e5)^0.4, 100)
  expect_identical(mr_cv(pilot = TRUE, A = 1e5, sampling_cost = 1,
                         gamma = .25)$value, 100)
  expect_identical(mr_smd(pilot = TRUE, A = 1e6, sampling_cost = 5,
                          gamma = .25)$value, 100)
  # A formula value that is not a whole number is rounded up as before.
  expect_identical(mr_cv(pilot = TRUE, A = 4e5, sampling_cost = 75)$value, 18)
  expect_identical(mr_cv(pilot = TRUE, A = 101^2.5, sampling_cost = 1,
                         gamma = .25)$value, 101)
  expect_identical(mr_cv(pilot = TRUE, A = 100.5^2.5, sampling_cost = 1,
                         gamma = .25)$value, 101)
})

test_that("A divided by sampling_cost must be finite", {
  for (f in list(function() mr_cv(pilot = TRUE, A = 1e300, sampling_cost = 1e-10),
                 function() mr_smd(pilot = TRUE, A = 1e300, sampling_cost = 1e-10))) {
    e <- tryCatch(f(), error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, c("A", "sampling_cost"))
  }
})

test_that("scores that differ only by rounding error give an unusable estimate", {
  # Equal gains computed as posttest minus pretest differ in the last bits.
  pre <- c(1.1, 2.3, 3.7, 4.2, 5.9, 6.1, 7.3, 8.8, 9.4, 10.6, 11.2, 12.7)
  gains <- (pre + 0.3) - pre
  expect_gt(var(gains), 0)
  g <- mr_cv(gains, A = 4e5, sampling_cost = 75, verbose = TRUE)
  expect_true(attr(g, "variance_unusable"))
  expect_false(attr(g, "is_satisfied"))
  expect_identical(mr_term(g, "v2_n"), 0)
  # An estimate that is zero to rounding at the zero crossing of Equation 32.
  t <- 0.504532046616077
  z <- (1 - t) * c(1, 1, 1, 1, 1, 10) + t * (1:6)
  zz <- mr_cv(z, A = 100, sampling_cost = 1, verbose = TRUE)
  expect_true(attr(zz, "variance_unusable"))
  expect_identical(mr_term(zz, "v2_n"), 0)
  # Scores with real variation are unaffected.
  expect_false(attr(mr_cv(mr_books, A = 4e5, sampling_cost = 75),
                    "variance_unusable"))
})
