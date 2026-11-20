test_that("ss_power_R2() solves for the necessary sample size", {
  result <- ss_power_R2(true_R2 = 0.3, p = 5, alpha_level = 0.05, desired_power = 0.8)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_true("necessary_N" %in% result$term)
  expect_true("actual_power" %in% result$term)
  expect_gt(result$value[result$term == "necessary_N"], 0)
  expect_gte(result$value[result$term == "actual_power"], 0.8)
})

test_that("ss_power_R2() with specified_N returns realized power", {
  result <- ss_power_R2(true_R2 = 0.3, p = 5, alpha_level = 0.05, specified_N = 100)
  expect_true("specified_N" %in% result$term)
  expect_equal(result$value[result$term == "specified_N"], 100)
  expect_gt(result$value[result$term == "actual_power"], 0.99)
})

test_that("ss_power_R2() requires fewer subjects for a larger effect", {
  small <- ss_power_R2(true_R2 = 0.10, p = 5, alpha_level = 0.05, desired_power = 0.8)$value[1]
  large <- ss_power_R2(true_R2 = 0.30, p = 5, alpha_level = 0.05, desired_power = 0.8)$value[1]
  expect_gt(small, large)
})

test_that("ss_power_R2() random predictors requires at least as many subjects as fixed predictors", {
  # Cohen (1988) noncentral-F overstates power under random predictors and so
  # understates required N. The Lee (1971) random-predictor path corrects this,
  # so N_random >= N_fixed across plausible scenarios.
  scenarios <- expand.grid(rho2 = c(0.10, 0.25, 0.50),
                           p    = c(3, 5),
                           pow  = c(0.80, 0.90))
  for (i in seq_len(nrow(scenarios))) {
    s <- scenarios[i, ]
    n_fixed  <- ss_power_R2(true_R2 = s$rho2, p = s$p, desired_power = s$pow,
                            random_predictors = FALSE)$value[1]
    n_random <- ss_power_R2(true_R2 = s$rho2, p = s$p, desired_power = s$pow,
                            random_predictors = TRUE)$value[1]
    expect_gte(n_random, n_fixed)
  }
})

test_that("ss_power_R2() random predictors at specified N has noncentral_f_parm NA", {
  res <- ss_power_R2(true_R2 = 0.3, p = 5, specified_N = 50,
                     random_predictors = TRUE)
  expect_true(is.na(res$value[res$term == "noncentral_f_parm"]))
})

test_that("ss_power_R2() respects random_predictors", {
  # The fixed-predictors and random-predictors paths give different
  # noncentralities at the same population R^2 and N; verify the
  # function actually branches.
  r_fixed  <- ss_power_R2(true_R2 = 0.3, p = 5, specified_N = 50,
                          random_predictors = FALSE)$value[2]
  r_random <- ss_power_R2(true_R2 = 0.3, p = 5, specified_N = 50,
                          random_predictors = TRUE)$value[2]
  expect_false(isTRUE(all.equal(r_fixed, r_random, tolerance = 1e-3)))
})

power_R2_value <- function(...) {
  res <- ss_power_R2(...)
  res$value[res$term == "actual_power"]
}

test_that("ss_power_R2() has size alpha_level when true_R2 is 0", {
  for (rp in c(FALSE, TRUE)) {
    expect_equal(power_R2_value(true_R2 = 0, null_R2 = 0, p = 5,
                                specified_N = 200, random_predictors = rp),
                 0.05, tolerance = 1e-10)
  }
})

# --- A nonzero null value is refused ------------------------------------------
#
# ss_power_R2() plans the omnibus test. DMAR 1.0.x accepted a nonzero null_R2;
# it now stops with a classed refusal that says what the earlier calculation
# was and writes the call to ss_power_R2_change().

test_that("ss_power_R2() refuses a nonzero null_R2 with the replacement call", {
  calls <- list(
    function() ss_power_R2(true_R2 = 0.3, null_R2 = 0.25, p = 4,
                           specified_N = 85, alpha_level = 0.01),
    function() ss_power_R2(true_R2 = 0.3, null_R2 = 0.25, p = 4,
                           desired_power = 0.80),
    function() ss_power_R2(true_R2 = 0.3, null_R2 = 0.25, p = 4,
                           specified_N = 85, random_predictors = FALSE),
    function() ss_power_R2(true_R2 = 0.3, null_R2 = 0.25, p = 4,
                           desired_power = 0.80, random_predictors = FALSE)
  )
  for (call_i in calls) {
    cnd <- expect_error(call_i(), class = "dmar_null_R2_unsupported")
    expect_s3_class(cnd, "dmar_refusal")
    expect_identical(class(cnd), c("dmar_null_R2_unsupported", "dmar_refusal",
                                   "error", "condition"))
    msg <- conditionMessage(cnd)
    expect_match(msg, "'null_R2' must be 0", fixed = TRUE)
    expect_match(msg, "a 'null_R2' of 0.25 was given", fixed = TRUE)
    expect_match(msg, paste0("ss_power_R2_change(true_R2_full = 0.3, ",
                             "true_R2_restricted = 0.25, p_full = "),
                 fixed = TRUE)
    expect_match(msg, "If 4 counted all of the predictors, 'p_full' is 4",
                 fixed = TRUE)
    expect_match(msg, paste0("ss_power_R2(true_f2 = (0.3 - 0.25) / (1 - 0.3), ",
                             "p = 4, random_predictors = FALSE)"), fixed = TRUE)
    expect_match(msg, "ci_R2()", fixed = TRUE)
    # The fields a program reads.
    expect_identical(cnd$code, "null_R2_unsupported")
    expect_identical(cnd$null_R2, 0.25)
    expect_identical(cnd$true_R2, 0.3)
    expect_null(cnd$true_f2)
    expect_identical(cnd$p, 4)
    expect_null(conditionCall(cnd))
  }
  # R cuts an error message at 8170 bytes; this one is well inside.
  expect_lt(nchar(conditionMessage(cnd), "bytes"), 4000)
})

test_that("ss_power_R2() refuses a nonzero null_R2 on the true_f2 entry point", {
  for (N in list(NULL, 50)) {
    for (rp in c(TRUE, FALSE)) {
      cnd <- expect_error(
        ss_power_R2(true_f2 = 0.1, null_R2 = 0.2, p = 5, specified_N = N,
                    desired_power = 0.80, random_predictors = rp),
        class = "dmar_null_R2_unsupported")
      expect_s3_class(cnd, "dmar_refusal")
      msg <- conditionMessage(cnd)
      expect_match(msg, "a 'null_R2' of 0.2 was given", fixed = TRUE)
      expect_match(msg, "ss_power_R2_change(true_f2 = 0.1, p_full = ",
                   fixed = TRUE)
      expect_match(msg, paste0("ss_power_R2(true_f2 = 0.1, p = 5, ",
                               "random_predictors = FALSE)"), fixed = TRUE)
      expect_identical(cnd$true_f2, 0.1)
      expect_null(cnd$true_R2)
      expect_identical(cnd$null_R2, 0.2)
      expect_identical(cnd$p, 5)
    }
  }
})

test_that("the null_R2 refusal says there is no increment when true_R2 is not above null_R2", {
  # The two calls the message otherwise writes would stop on a negative
  # effect size. DMAR 1.0.0 returned alpha_level (.05) for this call.
  for (rp in c(TRUE, FALSE)) {
    for (a in c(0.10, 0.25)) {
      cnd <- expect_error(
        ss_power_R2(true_R2 = a, null_R2 = 0.25, p = 4, specified_N = 85,
                    random_predictors = rp),
        class = "dmar_null_R2_unsupported")
      msg <- conditionMessage(cnd)
      expect_match(msg, paste0("'true_R2' (", a, ") is not above 'null_R2' ",
                               "(0.25), so there is no increment to test."),
                   fixed = TRUE)
      expect_no_match(msg, "ss_power_R2(true_f2 = ", fixed = TRUE)
      expect_no_match(msg, "true_R2_restricted = ", fixed = TRUE)
      expect_match(msg, "ss_power_R2_change()", fixed = TRUE)
      expect_match(msg, "ci_R2()", fixed = TRUE)
    }
  }
})

test_that("the true_f2 refusal on random predictors names the call MBESS's number comes from", {
  cnd <- expect_error(
    ss_power_R2(true_f2 = 0.1, null_R2 = 0.2, p = 5, specified_N = 100),
    class = "dmar_null_R2_unsupported")
  expect_match(conditionMessage(cnd),
               "ss_power_R2(true_f2 = 0.1, p = 5, random_predictors = FALSE)",
               fixed = TRUE)
  expect_no_match(conditionMessage(cnd), "1.0.x", fixed = TRUE)
  # DMAR 1.0.0 returned 0.9975493 for the refused call (random predictors)
  # and 0.6377504 for the same call without null_R2.
  expect_equal(power_R2_value(true_R2 = (0.1 + 0.2) / (1 + 0.1), p = 5,
                              specified_N = 100),
               0.9975493, tolerance = 1e-6)
  expect_equal(ss_power_R2(true_f2 = 0.1, p = 5,
                           specified_N = 100)$value[2],
               0.6377504, tolerance = 1e-6)
})

test_that("the calls the null_R2 refusal writes run and give the stated numbers", {
  # The number MBESS (5.0.1, ss.power.R2(Population.R2 = .30, Null.R2 = .25,
  # p = 4, Specified.N = 85, alpha.level = .01)) and the fixed-predictor
  # calculation of DMAR 1.0.0 returned: power .2225083 with
  # noncentrality 6.071429, which is Cohen's f^2 (u + v + 1) at N = 90 with
  # five predictors in the restricted model.
  mbess <- ss_power_R2(true_f2 = (0.3 - 0.25) / (1 - 0.3), p = 4,
                       specified_N = 85, alpha_level = 0.01,
                       random_predictors = FALSE)
  expect_equal(mbess$value[mbess$term == "actual_power"], 0.222508251,
               tolerance = 1e-8)
  expect_equal(mbess$value[mbess$term == "noncentral_f_parm"], 6.07142857,
               tolerance = 1e-8)
  # .25 against .20 with p = 5 at N = 200, the value the candidate returned
  # under fixed predictors.
  expect_equal(power_R2_value(true_f2 = 0.05 / 0.75, p = 5, specified_N = 200,
                              random_predictors = FALSE),
               0.8046769, tolerance = 1e-6)
  # The suggested call, completed: four predictors added to five.
  inc <- ss_power_R2_change(true_R2_full = 0.3, true_R2_restricted = 0.25,
                            p_full = 9, p_restricted = 5, specified_N = 90,
                            alpha_level = 0.01, random_predictors = FALSE)
  expect_equal(inc$value[inc$term == "actual_power"], 0.241296492,
               tolerance = 1e-8)
})

test_that("ss_power_R2() results with null_R2 = 0 are those of DMAR 1.0.0", {
  # Values returned by the released DMAR 1.0.0 (and by the 1.1.0 candidate
  # df120a0d) with the argument named population_R2.
  expect_equal(power_R2_value(true_R2 = 0.3, p = 5, specified_N = 50,
                              random_predictors = TRUE),
               0.91700234132275138, tolerance = 1e-10)
  expect_equal(power_R2_value(true_R2 = 0.3, p = 5, specified_N = 50,
                              random_predictors = FALSE),
               0.93965427440688798, tolerance = 1e-10)

  cells <- expand.grid(rho2 = c(0.02, 0.1, 0.3, 0.6), p = c(1, 5, 12),
                       pow = c(0.5, 0.8, 0.95))
  n_random <- c(192, 39, 13, 6, 350, 70, 24, 12, 498, 102, 36, 21,
                390, 76, 23, 10, 637, 124, 39, 17, 863, 170, 55, 26,
                644, 124, 37, 15, 980, 189, 57, 23, 1283, 249, 76, 33)
  n_fixed  <- c(191, 37, 12, 6, 349, 69, 23, 12, 496, 101, 35, 20,
                387, 73, 21, 8, 635, 122, 36, 16, 861, 167, 52, 25,
                639, 119, 33, 12, 975, 184, 53, 20, 1278, 244, 72, 30)
  n_of <- function(i, rp) {
    ss_power_R2(true_R2 = cells$rho2[i], p = cells$p[i],
                desired_power = cells$pow[i], random_predictors = rp)$value[1]
  }
  expect_equal(vapply(seq_len(nrow(cells)), n_of, numeric(1), rp = TRUE), n_random)
  expect_equal(vapply(seq_len(nrow(cells)), n_of, numeric(1), rp = FALSE), n_fixed)
})

test_that("ss_power_R2() is finite at a very small nonzero true_R2", {
  # The discriminant of Lee's approximation, computed by subtraction, came out
  # negative near rho^2 = 1e-8, and the power was NaN with a warning.
  expect_no_warning(
    tiny <- power_R2_value(true_R2 = 1e-8, p = 10, specified_N = 12)
  )
  expect_true(is.finite(tiny))
  expect_equal(tiny, 0.05, tolerance = 1e-6)
  expect_true(is.finite(DMAR:::.lee_random_R2_cdf(0.5, 1e-8, 12, 10)))

  # At rho^2 = 0 the quantile is the central F critical value of the omnibus
  # test, and it is continuous there.
  f_crit  <- qf(0.95, 5, 194)
  R2_crit <- 5 * f_crit / (5 * f_crit + 194)
  expect_equal(DMAR:::.lee_random_R2_quantile(0.95, 0, 200, 5), R2_crit,
               tolerance = 1e-14)
  expect_equal(DMAR:::.lee_random_R2_quantile(0.95, 1e-10, 200, 5), R2_crit,
               tolerance = 1e-8)
})

test_that("the quantile of Lee's approximation inverts its distribution function", {
  cells <- expand.grid(q = c(0.9, 0.95, 0.99), rho2 = c(0.1, 0.5, 0.9),
                       design = 1:3)
  N <- c(15, 200, 60)
  p <- c(1, 5, 20)
  for (i in seq_len(nrow(cells))) {
    s <- cells[i, ]
    x <- DMAR:::.lee_random_R2_quantile(s$q, s$rho2, N[s$design], p[s$design])
    expect_equal(DMAR:::.lee_random_R2_cdf(x, s$rho2, N[s$design], p[s$design]),
                 s$q, tolerance = 1e-9)
  }
})

test_that("the upper point of the distribution of R^2 at a nonzero rho^2 is the dual of the one-sided lower limit of ci_R2()", {
  # No exported function plans a test against a nonzero null value of rho^2
  # in this release. The critical values such a test would use are checked
  # here against ci_R2(): the sample R^2 just above the upper alpha point of
  # its distribution at rho^2 = null has a one-sided lower limit above the
  # null, and the one just below has a limit below it. Under random
  # predictors the two rules agree to the tolerance of the bisection ci_R2()
  # makes for its limit (1e-5 in probability), well inside the offset of
  # .001 used here.
  designs <- data.frame(N = c(200, 40, 25), p = c(5, 8, 2),
                        null = c(0.20, 0.40, 0.60))
  ci_lower <- function(R2, N, p, rp) {
    ci <- ci_R2(R2 = R2, N = N, p = p, alpha_lower = 0.05, alpha_upper = 0,
                random_predictors = rp)
    ci$value[ci$term == "ci_lower"]
  }
  for (i in seq_len(nrow(designs))) {
    s <- designs[i, ]
    df_2 <- s$N - s$p - 1
    f_crit <- qf(0.95, s$p, df_2, ncp = s$N * s$null / (1 - s$null))
    crit_fixed  <- s$p * f_crit / (s$p * f_crit + df_2)
    crit_random <- DMAR:::.lee_random_R2_quantile(0.95, rho2 = s$null,
                                                  N = s$N, p = s$p)
    expect_lt(ci_lower(crit_fixed - 0.001, s$N, s$p, FALSE), s$null)
    expect_gt(ci_lower(crit_fixed + 0.001, s$N, s$p, FALSE), s$null)
    expect_lt(ci_lower(crit_random - 0.001, s$N, s$p, TRUE), s$null)
    expect_gt(ci_lower(crit_random + 0.001, s$N, s$p, TRUE), s$null)
  }
})

test_that("ss_power_R2() refuses arguments outside their ranges", {
  expect_error(ss_power_R2(true_R2 = 1, p = 5, specified_N = 50),
               "true_R2")
  expect_error(ss_power_R2(p = 5, specified_N = 50),
               "'true_R2' or 'true_f2'")
  expect_error(ss_power_R2(true_R2 = 0.3, null_R2 = NA, p = 5,
                           specified_N = 50), "null_R2")
  expect_error(ss_power_R2(true_R2 = 0.3, null_R2 = c(0, 0.2), p = 5,
                           specified_N = 50), "null_R2")
  expect_error(ss_power_R2(true_R2 = 0.3, null_R2 = 1, p = 5,
                           specified_N = 50), "null_R2")
  expect_error(ss_power_R2(true_R2 = 0.3, null_R2 = -0.1, p = 5,
                           specified_N = 50), "null_R2")
  expect_error(ss_power_R2(true_R2 = 0.3, p = 5.5, specified_N = 50),
               "'p' must be")
  expect_error(ss_power_R2(true_R2 = 0.3, p = 0, specified_N = 50),
               "'p' must be")
  expect_error(ss_power_R2(true_R2 = 0.3, p = 5, specified_N = 50,
                           alpha_level = 0), "alpha_level")
  expect_error(ss_power_R2(true_R2 = 0.3, p = 5, specified_N = 50,
                           alpha_level = 1), "alpha_level")
  expect_error(ss_power_R2(true_f2 = -1, p = 5, specified_N = 50),
               "true_f2")
  expect_error(ss_power_R2(true_R2 = 0.3, true_f2 = 0.1, p = 5,
                           specified_N = 50), "do not specify 'true_f2'")
  expect_error(ss_power_R2(true_R2 = 0.3, p = 5, specified_N = 50,
                           random_predictors = NA), "random_predictors")
})

test_that("ss_power_R2() keeps the validation message for an invalid null_R2", {
  # A value that is not a single number in [0, 1) is a validation error, not
  # the refusal of a nonzero null.
  for (bad in list(NA, c(0, 0.2), 1, -0.1, "0.2")) {
    cnd <- expect_error(ss_power_R2(true_R2 = 0.3, null_R2 = bad, p = 5,
                                    specified_N = 50),
                        "'null_R2' must be a single numeric value in [0, 1).",
                        fixed = TRUE)
    expect_false(inherits(cnd, "dmar_null_R2_unsupported"))
  }
  # The other arguments are checked first, so their own message is the one
  # given when they are invalid too.
  expect_error(ss_power_R2(true_R2 = 1.2, null_R2 = 0.2, p = 5,
                           specified_N = 50), "'true_R2' must be")
  expect_error(ss_power_R2(null_R2 = 0.2, p = 5, specified_N = 50),
               "'true_R2' or 'true_f2'")
})

test_that("ss_power_R2() true_f2 and true_R2 entry points agree at null_R2 = 0", {
  for (rp in c(FALSE, TRUE)) {
    expect_equal(ss_power_R2(true_f2 = 0.3 / 0.7, p = 5, desired_power = 0.80,
                             random_predictors = rp),
                 ss_power_R2(true_R2 = 0.3, p = 5, desired_power = 0.80,
                             random_predictors = rp))
    expect_equal(ss_power_R2(true_f2 = 0.3 / 0.7, p = 5, specified_N = 50,
                             random_predictors = rp),
                 ss_power_R2(true_R2 = 0.3, p = 5, specified_N = 50,
                             random_predictors = rp))
  }
})

test_that("ss_power_R2() refuses a desired_power of 0 or 1, as ss_power_R2_change() does", {
  # Through 1.0.x both ends were accepted (desired_power = 1 returned 346
  # for true_R2 = .3 and p = 5). One rule now holds in every power planner.
  for (rp in c(TRUE, FALSE)) {
    for (bad in c(0, 1)) {
      err <- expect_error(
        ss_power_R2(true_R2 = 0.3, p = 5, desired_power = bad,
                    random_predictors = rp),
        class = "dmar_desired_power_invalid")
      expect_identical(err$value, bad)
      expect_match(conditionMessage(err), "strictly between 0 and 1",
                   fixed = TRUE)
    }
  }
  # Not used when a sample size is given, and refused all the same.
  expect_error(ss_power_R2(true_R2 = 0.3, p = 5, specified_N = 60,
                           desired_power = 1),
               class = "dmar_desired_power_invalid")
  expect_s3_class(ss_power_R2(true_R2 = 0.3, p = 5, specified_N = 60,
                              desired_power = .5), "dmar_tbl")
  # A value just inside is planned.
  expect_gt(ss_power_R2(true_R2 = 0.3, p = 5, desired_power = 0.999)$value[1],
            ss_power_R2(true_R2 = 0.3, p = 5, desired_power = 0.99)$value[1])
})

test_that("ss_power_R2() validates desired_power on the true_f2 entry point", {
  # Previously, desired_power validation was inside the true_R2 branch,
  # so an invalid desired_power supplied alongside true_f2 alone would drive
  # the iterative search into an infinite loop.
  expect_error(
    ss_power_R2(true_f2 = 1, desired_power = 1.5, p = 5,
                random_predictors = FALSE),
    "desired_power"
  )
  expect_error(
    ss_power_R2(true_f2 = 1, desired_power = -0.5, p = 5,
                random_predictors = TRUE),
    "desired_power"
  )
})

test_that("ss_power_R2() rejects a specified_N with no residual degrees of freedom (MEDIUM-03)", {
  # df_2 = N - p - 1 = 0 at specified_N = p + 1: the F test is undefined, so the
  # function errors rather than returning a numeric-looking power of 0.
  expect_error(ss_power_R2(true_R2 = 0.30, p = 5, specified_N = 6,
                           random_predictors = TRUE), "whole number")
  expect_error(ss_power_R2(true_R2 = 0.30, p = 5, specified_N = 6,
                           random_predictors = FALSE), "whole number")
})

test_that("ss_power_R2() random-predictor power matches a Monte Carlo simulation", {
  skip_on_cran()
  # The random-predictor path uses the Lee (1971) approximation to the sampling
  # distribution of R^2 when the predictors are themselves a draw from a joint
  # multivariate normal. This pins the analytic power to ground truth: generate
  # data with the predictors RE-DRAWN each replication (the random-predictor
  # model), reject by the usual overall F-test, and compare the empirical
  # rejection rate to ss_power_R2()'s analytic power.
  set.seed(113)
  p <- 3; rho2 <- 0.20; N <- 60; alpha <- 0.05; G <- 4000

  analytic <- ss_power_R2(true_R2 = rho2, p = p, specified_N = N,
                          random_predictors = TRUE,
                          alpha_level = alpha)$value[2]

  # Equal coefficients on independent standard-normal predictors give a
  # population R^2 of exactly rho2 when the error variance is 1.
  b <- sqrt(rho2 / (p * (1 - rho2)))
  reject <- logical(G)
  for (g in seq_len(G)) {
    X    <- matrix(rnorm(N * p), N, p)              # predictors random each rep
    y    <- as.numeric(X %*% rep(b, p)) + rnorm(N)
    fobs <- summary(lm(y ~ X))$fstatistic
    reject[g] <- pf(fobs[1], fobs[2], fobs[3], lower.tail = FALSE) < alpha
  }
  empirical <- mean(reject)

  # The Lee approximation is excellent here; the Monte Carlo standard error is
  # about 0.005, so a 0.02 band is comfortable and still a real test.
  expect_lt(abs(analytic - empirical), 0.02)

  # And the fixed-predictor (Cohen, 1988) formula over-states the random-design
  # power, the discrepancy that motivates the random-predictor default.
  fixed <- ss_power_R2(true_R2 = rho2, p = p, specified_N = N,
                       random_predictors = FALSE, alpha_level = alpha)$value[2]
  expect_gt(fixed, analytic)
})

test_that("ss_power_R2 returns the true minimum sample size (MEDIUM-02)", {
  pw <- function(n) { r <- ss_power_R2(true_R2 = 0.3, p = 5, specified_N = n)
    r$value[r$term == "actual_power"] }
  target <- pw(7) - 1e-4                        # minimum admissible is p + 2 = 7 (df_2 >= 1)
  r <- ss_power_R2(true_R2 = 0.3, p = 5, desired_power = target)
  expect_equal(r$value[r$term == "necessary_N"], 7)
})

test_that("ss_power_R2 rejects invalid fixed sample sizes (MEDIUM-03)", {
  expect_error(ss_power_R2(true_R2 = 0.3, p = 5, specified_N = 5), "whole number")   # no residual df
  expect_error(ss_power_R2(true_R2 = 0.3, p = 2, specified_N = 20.5), "whole number") # fractional
  expect_silent(ss_power_R2(true_R2 = 0.3, p = 2, specified_N = 20))
})

test_that("ss_power_R2() stops before the search when true_R2 is 0", {
  for (rp in c(TRUE, FALSE)) {
    expect_error(ss_power_R2(true_R2 = 0, p = 4, random_predictors = rp),
                 paste0("With 'true_R2' of 0 the null hypothesis is true, ",
                        "so power equals 'alpha_level' at every sample size and ",
                        "no sample size attains 'desired_power'."),
                 fixed = TRUE)
    expect_error(ss_power_R2(true_f2 = 0, p = 4, random_predictors = rp),
                 "the null hypothesis is true")
    expect_error(ss_power_R2(true_f2 = 0, p = 4, random_predictors = rp),
                 class = "dmar_effect_zero")
  }
})

test_that("ss_power_R2() searches up to ten million and stops there", {
  # A necessary size above 100,000, found without a warning and without a
  # search in steps of one.
  big <- expect_no_warning(
    ss_power_R2(true_R2 = 5e-5, p = 5, desired_power = 0.80,
                random_predictors = FALSE)
  )
  n_big <- big$value[big$term == "necessary_N"]
  expect_equal(n_big, 256545)
  expect_gte(big$value[big$term == "actual_power"], 0.80)
  expect_lt(power_R2_value(true_R2 = 5e-5, p = 5,
                           specified_N = n_big - 1,
                           random_predictors = FALSE), 0.80)
  big_random <- ss_power_R2(true_R2 = 5e-5, p = 5, desired_power = 0.80)
  expect_equal(big_random$value[big_random$term == "necessary_N"], 256548)
  # An effect so small that ten million cases do not attain the power.
  expect_error(
    ss_power_R2(true_R2 = 1e-7, p = 5, desired_power = 0.80,
                random_predictors = FALSE),
    "Failed to reach 'desired_power'"
  )
})

test_that("the shared search returns no sample size where the power is missing", {
  # R's noncentral F routine warns of lost accuracy beyond a noncentrality of
  # about one million, and the power functions then return NA in place of a
  # number. With a null value of zero no call was found that reaches the
  # warning on this build of R, so the handling of a missing power is tested
  # on the search itself with a power function that is missing above a size.
  msg <- function(N) DMAR:::.R2_imprecise_message(N, 0.5)
  power_fn <- function(limit) function(N) {
    if (N > limit) NA_real_ else stats::pnorm(N / 100 - 3)
  }
  # The answer (385) lies below the sizes without a power: it is returned.
  found <- DMAR:::.R2_power_search_N(power_fn(400), N_min = 7,
                                     desired_power = 0.80,
                                     imprecise_message = msg)
  expect_equal(found$N, 385)
  expect_gte(found$power, 0.80)
  expect_lt(stats::pnorm(384 / 100 - 3), 0.80)
  # The answer lies among the sizes without a power: the search stops and
  # says above which size it lies and why.
  expect_error(
    DMAR:::.R2_power_search_N(power_fn(300), N_min = 7, desired_power = 0.80,
                              imprecise_message = msg),
    "is above 300\\. At a sample size of 301 the noncentrality parameter"
  )
  # No power at the smallest size.
  expect_error(
    DMAR:::.R2_power_search_N(power_fn(5), N_min = 7, desired_power = 0.80,
                              imprecise_message = msg),
    "This happens only at very large sample sizes."
  )
  # A target met at the smallest size returns that size.
  expect_equal(DMAR:::.R2_power_search_N(function(N) 0.9, N_min = 7,
                                         desired_power = 0.80,
                                         imprecise_message = msg)$N, 7)
  # The guard passes a value through and turns the warning into NA.
  expect_equal(DMAR:::.R2_change_power_guarded(50, 0.3 / 0.7, 5, 0, 0.05, FALSE),
               0.93965427440688798, tolerance = 1e-10)
  # At an N of ten million the power is computed without a warning.
  expect_no_warning(expect_equal(
    power_R2_value(true_f2 = 0.5, p = 5, specified_N = 1e7,
                   random_predictors = FALSE), 1))
})
