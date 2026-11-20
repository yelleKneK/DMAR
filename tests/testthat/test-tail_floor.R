# The floor on a nonzero tail area (.dmar_min_tail_area, 1e-9, and
# .dmar_check_tail_floor() and .dmar_check_assurance_tail_floor() in
# R/dmar_conditions.R), checked by ci_nc_t() for every interval and plan
# whose limits come from it and by the intervals and plans of the
# regression coefficient (ci_reg_coef(), ci_rc(), ci_src(),
# ss_aipe_reg_coef(), ss_aipe_rc(), ss_aipe_src()), the contrast (ci_c(),
# ci_sc(), ci_sc_ancova(), ss_aipe_c(), ss_aipe_sc(), ss_aipe_sc_ancova()),
# and the correlation (ci_r(), ss_aipe_r()), and by ci_sm() and ci_cv(),
# which give ci_nc_t() tails of their own; and the upper quantiles these
# take with lower.tail = FALSE, which keep a small tail area that
# 1 - alpha rounded away (DMAR 1.0.0 gave ci_c() an upper limit of Inf and
# ci_r() a limit of 1 at a tail of 1e-20).

# The refusal: class dmar_argument_invalid, reason "tail_below_floor", the
# argument or arguments named, and the floor carried.
expect_floor_refusal <- function(object, arg) {
  e <- expect_error(object, class = "dmar_argument_invalid")
  expect_identical(e$reason, "tail_below_floor")
  expect_identical(e$arg, arg)
  expect_identical(e$floor, 1e-9)
  invisible(e)
}

test_that("the floor is 1e-9, and the check refuses a nonzero tail below it in the name of that tail", {
  check <- DMAR:::.dmar_check_tail_floor
  expect_identical(DMAR:::.dmar_min_tail_area, 1e-9)
  # Zero (a one-sided interval), the floor itself, and anything the
  # function's own checks deal with pass.
  expect_null(check(0, .05))
  expect_null(check(1e-9, 1e-9))
  expect_null(check(.025, .025))
  expect_null(check(NULL, NA))
  expect_null(check(-1, "a"))
  e <- expect_floor_refusal(check(.025, 1e-20), "alpha_upper")
  expect_identical(e$value, c(alpha_upper = 1e-20))
  expect_match(conditionMessage(e), "^'alpha_upper' is 1e-20, a tail area above 0 and below 1e-9")
  expect_match(conditionMessage(e), "Give a tail area of at least 1e-9, or 0 for a one-sided interval.", fixed = TRUE)
  e <- expect_floor_refusal(check(5e-10, .025, zero_ok = FALSE), "alpha_lower")
  expect_match(conditionMessage(e), "Give a tail area of at least 1e-9.", fixed = TRUE)
  e <- expect_floor_refusal(check(1e-12, 1e-20), c("alpha_lower", "alpha_upper"))
  expect_match(conditionMessage(e), "^'alpha_lower' and 'alpha_upper' are 1e-12 and 1e-20")
  e <- expect_floor_refusal(check(1e-12, 1e-12), c("alpha_lower", "alpha_upper"))
  expect_match(conditionMessage(e), "^'alpha_lower' and 'alpha_upper' are both 1e-12")
  # With the other tail not given (NULL), the tail below the floor is the
  # one named and carried (c() dropped the NULL, and the message read
  # "'alpha_upper' is NA").
  e <- expect_floor_refusal(check(NULL, 1e-20), "alpha_upper")
  expect_identical(e$value, c(alpha_upper = 1e-20))
  expect_match(conditionMessage(e), "^'alpha_upper' is 1e-20, a tail area")
  e <- expect_floor_refusal(check(1e-20, NULL), "alpha_lower")
  expect_identical(e$value, c(alpha_lower = 1e-20))
  # Tails formed from a level: the level is named, with the largest one
  # accepted, which leaves tails at the floor in double precision.
  cl <- 1 - 1e-12
  e <- expect_floor_refusal(check((1 - cl) / 2, (1 - cl) / 2, conf_level = cl),
                            "conf_level")
  expect_match(conditionMessage(e), "Give a conf_level of at most 0.999999998.", fixed = TRUE)
  expect_null(check((1 - .999999998) / 2, (1 - .999999998) / 2,
                    conf_level = .999999998))
  # No em dash in the messages.
  expect_false(grepl("\u2014", conditionMessage(e)))
})

test_that("an assurance whose tail areas fall below the floor is refused in its own name", {
  check <- DMAR:::.dmar_check_assurance_tail_floor
  expect_null(check(NULL))
  expect_null(check(.99))
  expect_null(check(.999999998))
  e <- expect_error(check(1 - 1e-10), class = "dmar_argument_invalid")
  expect_identical(e$arg, "assurance")
  expect_identical(e$reason, "tail_below_floor")
  expect_match(conditionMessage(e), "Give an assurance of at most 0.999999998.", fixed = TRUE)
})

test_that("ci_nc_t() refuses a nonzero tail below the floor, and every interval built on it does", {
  expect_floor_refusal(ci_nc_t(ncp = 5, df = 47, conf_level = NULL,
                               alpha_lower = .025, alpha_upper = 1e-20),
                       "alpha_upper")
  expect_floor_refusal(ci_nc_t(ncp = 5, df = 47, conf_level = 1 - 1e-12),
                       "conf_level")
  # A tail of 0 is still a one-sided interval, beside a tail at the floor.
  res <- ci_nc_t(ncp = 5, df = 47, conf_level = NULL, alpha_lower = 0,
                 alpha_upper = 1e-9)
  expect_identical(res$value[1], -Inf)
  expect_true(is.finite(res$value[2]))
  expect_floor_refusal(ci_smd(smd = .5, n_1 = 20, n_2 = 20, conf_level = 1 - 1e-12),
                       "conf_level")
  expect_floor_refusal(ci_sm(sm = .5, N = 20, conf_level = NULL,
                             alpha_lower = 1e-12, alpha_upper = .025),
                       "alpha_lower")
  expect_floor_refusal(ci_smd_c(smd_c = .5, n_C = 20, n_E = 20, conf_level = NULL,
                                alpha_lower = .025, alpha_upper = 1e-12),
                       "alpha_upper")
})

test_that("at the floor the ci_nc_t() limits put the requested tails at the observed t", {
  # The tails at the limits by the oracle that conditions on the normal
  # variate (helper-nc_t_oracle.R). At 1e-9 a limit is found from a tail
  # whose absolute error is up to about 1e-12, so the tails are compared
  # to 1e-3 of themselves (the lower one at t = 5 is 3.6e-4 off).
  for (x in list(c(5, 47), c(.5, 10))) {
    res <- ci_nc_t(ncp = x[1], df = x[2], conf_level = NULL,
                   alpha_lower = 1e-9, alpha_upper = 1e-9)
    expect_equal(unname(.oracle_tails(res, x[1], x[2])), c(1e-9, 1e-9),
                 tolerance = 1e-3)
  }
})

test_that("ci_reg_coef(), ci_rc(), and ci_src() refuse a tail below the floor, and at it the central limits are the upper-tail quantiles", {
  for (nc in c(FALSE, TRUE)) {
    expect_floor_refusal(ci_reg_coef(reg_coef_j = .5, se_reg_coef_j = .1, N = 50, p = 2,
                                     conf_level = NULL, alpha_lower = .025,
                                     alpha_upper = 1e-20, noncentral = nc),
                         "alpha_upper")
    expect_floor_refusal(ci_reg_coef(reg_coef_j = .5, se_reg_coef_j = .1, N = 50, p = 2,
                                     conf_level = NULL, alpha_lower = 1e-20,
                                     alpha_upper = .025, noncentral = nc),
                         "alpha_lower")
  }
  expect_floor_refusal(ci_reg_coef(reg_coef_j = .5, se_reg_coef_j = .1, N = 50, p = 2,
                                   conf_level = 1 - 1e-12), "conf_level")
  expect_floor_refusal(ci_rc(reg_coef_j = .3, se_reg_coef_j = .1, N = 50, p = 2,
                             conf_level = NULL, alpha_lower = 1e-12, alpha_upper = .025),
                       "alpha_lower")
  expect_floor_refusal(ci_src(std_reg_coef_j = .3, se_std_reg_coef_j = .1, N = 50, p = 2,
                              conf_level = NULL, alpha_lower = 1e-12, alpha_upper = .025),
                       "alpha_lower")
  # At the floor, in either tail.
  ci <- ci_reg_coef(reg_coef_j = .5, se_reg_coef_j = .1, N = 50, p = 2,
                    conf_level = NULL, alpha_lower = 1e-9, alpha_upper = 1e-9)
  expect_equal(ci$value, c(.5 + stats::qt(1e-9, 47) * .1, .5,
                           .5 + stats::qt(1e-9, 47, lower.tail = FALSE) * .1),
               tolerance = 1e-14)
  # A zero tail is still a one-sided bound.
  ci <- ci_reg_coef(reg_coef_j = .5, se_reg_coef_j = .1, N = 50, p = 2,
                    conf_level = NULL, alpha_lower = 1e-9, alpha_upper = 0)
  expect_identical(ci$value[3], Inf)
})

test_that("ci_c(), ci_sc(), and ci_sc_ancova() refuse a tail below the floor, and at it ci_c() takes the upper-tail quantile", {
  expect_floor_refusal(ci_c(means = c(1, 0), s_anova = 1, c_weights = c(1, -1),
                            n_per_group = c(20, 20), conf_level = NULL,
                            alpha_lower = .025, alpha_upper = 1e-20),
                       "alpha_upper")
  expect_floor_refusal(ci_c(means = c(1, 0), s_anova = 1, c_weights = c(1, -1),
                            n_per_group = 20, conf_level = 1 - 1e-12),
                       "conf_level")
  expect_floor_refusal(ci_sc(means = c(1, 0), s_anova = 1, c_weights = c(1, -1),
                             n_per_group = c(20, 20), conf_level = NULL,
                             alpha_lower = .025, alpha_upper = 1e-20),
                       "alpha_upper")
  expect_floor_refusal(ci_sc(means = c(1, 0), s_anova = 1, c_weights = c(1, -1),
                             n_per_group = 20, conf_level = 1 - 1e-12),
                       "conf_level")
  expect_floor_refusal(ci_sc_ancova(contrast = .5, s_ancova = 1, c_weights = c(1, -1),
                                    n_per_group = c(20, 20), cov_means = c(0, 0),
                                    SSwithin_x = 10, conf_level = 1 - 1e-12),
                       "conf_level")
  se <- sqrt(2 / 20)
  ci <- ci_c(means = c(1, 0), s_anova = 1, c_weights = c(1, -1),
             n_per_group = c(20, 20), conf_level = NULL, alpha_lower = 0,
             alpha_upper = 1e-9)
  expect_identical(ci$value[1], -Inf)
  expect_equal(ci$value[3], 1 + stats::qt(1e-9, 38, lower.tail = FALSE) * se,
               tolerance = 1e-14)
  expect_true(all(is.finite(ci_sc(means = c(1, 0), s_anova = 1, c_weights = c(1, -1),
                                  n_per_group = c(20, 20), conf_level = NULL,
                                  alpha_lower = 1e-9, alpha_upper = 1e-9)$value)))
})

test_that("ci_r() refuses a tail below the floor, and at it takes the upper-tail normal quantile", {
  expect_floor_refusal(ci_r(r = .3, N = 50, conf_level = NULL, alpha_lower = .025,
                            alpha_upper = 1e-20), "alpha_upper")
  expect_floor_refusal(ci_r(r = .3, N = 50, conf_level = NULL, alpha_lower = 1e-20,
                            alpha_upper = .025), "alpha_lower")
  expect_floor_refusal(ci_r(r = .3, N = 50, conf_level = 1 - 1e-12), "conf_level")
  ci <- ci_r(r = .3, N = 50, conf_level = NULL, alpha_lower = 1e-9,
             alpha_upper = .025)
  z <- atanh(.3)
  se <- 1 / sqrt(47)
  expect_equal(ci$value[c(1, 3)],
               tanh(z + c(-1, 1) * c(stats::qnorm(1e-9, lower.tail = FALSE),
                                     stats::qnorm(.025, lower.tail = FALSE)) * se),
               tolerance = 1e-14)
  # qnorm(1 - 1e-9) is 4.6e-9 too large; the limit is the exact one.
  expect_gt(abs(stats::qnorm(1 - 1e-9) - stats::qnorm(1e-9, lower.tail = FALSE)), 1e-9)
})

test_that("ss_aipe_sc() and ss_aipe_sc_ancova() refuse tails and assurances below the floor", {
  expect_floor_refusal(ss_aipe_sc(true_std_contrast = .5, c_weights = c(1, -1), width = .5,
                                  conf_level = NULL, alpha_lower = 1e-20,
                                  alpha_upper = .025), "alpha_lower")
  expect_floor_refusal(ss_aipe_sc(true_std_contrast = .5, c_weights = c(1, -1), width = .5,
                                  conf_level = 1 - 1e-12), "conf_level")
  expect_floor_refusal(ss_aipe_sc_ancova(true_std_contrast = .5, c_weights = c(1, -1),
                                         width = .5, conf_level = NULL,
                                         alpha_lower = 1e-20, alpha_upper = 1e-20),
                       c("alpha_lower", "alpha_upper"))
  expect_floor_refusal(ss_aipe_sc_ancova(true_std_contrast = .5, c_weights = c(1, -1),
                                         width = .5, conf_level = 1 - 1e-12),
                       "conf_level")
  # The published assurance step finds limits with tails of
  # (1 - assurance) / 2: at an assurance of 1 - 1e-10, 5e-11.
  a <- 1 - 1e-10
  for (call in list(
    function() ss_aipe_sc(true_std_contrast = .5, c_weights = c(1, -1), width = .5,
                          assurance = a),
    function() ss_aipe_sc_ancova(true_std_contrast = .5, c_weights = c(1, -1),
                                 width = .5, assurance = a),
    function() ss_aipe_sc_ancova(true_std_contrast = .5, c_weights = c(1, -1),
                                 width = .5, assurance = a,
                                 assurance_method = "simulation", G = 100,
                                 seed = 113),
    function() ss_aipe_sm(true_sm = .5, width = .5, assurance = a))) {
    e <- expect_error(call(), class = "dmar_argument_invalid")
    expect_identical(e$arg, "assurance")
    expect_identical(e$reason, "tail_below_floor")
  }
})

test_that("at the floor ss_aipe_sc() plans a contrast of 0 from the upper-tail normal quantiles", {
  # At a contrast of 0 the limits are the normal quantiles at any df, so the
  # plan is the smallest n with sqrt(2 / n) (z_lower + z_upper) <= width.
  z <- stats::qnorm(1e-9, lower.tail = FALSE) + stats::qnorm(.025, lower.tail = FALSE)
  plan <- ss_aipe_sc(true_std_contrast = 0, c_weights = c(1, -1), width = .5,
                     conf_level = NULL, alpha_lower = 1e-9, alpha_upper = .025)$value
  expect_lte(sqrt(2 / plan) * z, .5)
  expect_gt(sqrt(2 / (plan - 1)) * z, .5)
  expect_true(is.finite(ss_aipe_sc_ancova(true_std_contrast = .5, c_weights = c(1, -1),
                                          width = .5, conf_level = NULL,
                                          alpha_lower = 1e-9,
                                          alpha_upper = 1e-9)$value))
})

test_that("ss_aipe_reg_coef(), ss_aipe_rc(), and ss_aipe_src() refuse a tail below the floor, and at it the central plan takes the upper-tail quantiles", {
  args <- list(true_R2_Y_X = .1, true_R2_j_X_without_j = .2, p = 2,
               true_reg_coef_j = .1, width = .4)
  for (nc in c(FALSE, TRUE)) {
    expect_floor_refusal(do.call(ss_aipe_reg_coef, c(args, list(
      noncentral = nc, conf_level = NULL, alpha_lower = 1e-20, alpha_upper = .05))),
      "alpha_lower")
  }
  expect_floor_refusal(do.call(ss_aipe_reg_coef, c(args, list(conf_level = 1 - 1e-12))),
                       "conf_level")
  expect_floor_refusal(ss_aipe_rc(true_R2_Y_X = .1, true_R2_j_X_without_j = .2, p = 2,
                                  true_reg_coef_j = .1, width = .4, conf_level = NULL,
                                  alpha_lower = 1e-20, alpha_upper = .05),
                       "alpha_lower")
  expect_floor_refusal(ss_aipe_src(true_R2_Y_X = .1, true_R2_j_X_without_j = .2, p = 2,
                                   true_std_reg_coef_j = .1, width = .4, conf_level = NULL,
                                   alpha_lower = .05, alpha_upper = 1e-20),
                       "alpha_upper")
  # The plan at the floor with one predictor: the smallest N whose interval
  # at the expected squared multiple correlations, with the t quantiles of
  # the two tails taken in the upper tail, is no wider than the target.
  # With both population values 0 the expected squared multiple
  # correlation of Y is 1 / (N - 1) and that of X_j, on no other
  # predictor, is 0, so the standard error is 1 / sqrt(N - 1).
  width_at <- function(N) {
    (stats::qt(1e-9, N - 2, lower.tail = FALSE) +
       stats::qt(.025, N - 2, lower.tail = FALSE)) / sqrt(N - 1)
  }
  plan <- ss_aipe_reg_coef(true_R2_Y_X = 0, true_R2_j_X_without_j = 0,
                           p = 1, true_reg_coef_j = 0, width = .5,
                           conf_level = NULL, alpha_lower = 1e-9,
                           alpha_upper = .025)$value
  expect_lte(width_at(plan), .5)
  expect_gt(width_at(plan - 1), .5)
  # The largest level accepted is planned: 1,032 for the arguments above
  # (1,033 in DMAR 1.0.0, before the search from p + 2 and the expectation
  # on p - 1 predictors).
  expect_identical(do.call(ss_aipe_reg_coef, c(args, list(conf_level = .999999998)))$value[1],
                   1032)
})

test_that("ss_aipe_cv() names conf_level when the level leaves tails below the floor", {
  expect_floor_refusal(ss_aipe_cv(true_cv = .2, width = .1, conf_level = 1 - 1e-12),
                       "conf_level")
})

test_that("ss_aipe_cv() does not offer a tail of 0, which it refuses as alpha_upper", {
  e <- expect_floor_refusal(ss_aipe_cv(true_cv = .2, width = .1, conf_level = NULL,
                                       alpha_lower = .025, alpha_upper = 5e-10),
                            "alpha_upper")
  expect_match(conditionMessage(e), "Give a tail area of at least 1e-9.", fixed = TRUE)
  expect_no_match(conditionMessage(e), "one-sided")
})

test_that("ci_cv() names the tail the caller set below the floor, and a level as the level", {
  # The interval passes ci_nc_t() the two tails exchanged, so the refusal
  # there named the other tail.
  e <- expect_floor_refusal(ci_cv(cv = .2, N = 30, conf_level = NULL,
                                  alpha_lower = .025, alpha_upper = 5e-10),
                            "alpha_upper")
  expect_identical(e$value, c(alpha_upper = 5e-10))
  expect_floor_refusal(ci_cv(cv = .2, N = 30, conf_level = NULL,
                             alpha_lower = 5e-10, alpha_upper = .025),
                       "alpha_lower")
  expect_floor_refusal(ci_cv(mean = 5, sd = 1, N = 30, conf_level = NULL,
                             alpha_lower = .025, alpha_upper = 5e-10),
                       "alpha_upper")
  expect_floor_refusal(ci_cv(cv = .2, N = 50, conf_level = 1 - 1e-12), "conf_level")
  # At the floor, and with a tail of 0, the interval is formed.
  expect_true(all(is.finite(ci_cv(cv = .2, N = 30, conf_level = NULL,
                                  alpha_lower = 1e-9, alpha_upper = 1e-9)$value)))
  expect_identical(ci_cv(cv = .2, N = 30, conf_level = NULL, alpha_lower = 0,
                         alpha_upper = .05)$value[1], 0)
})

test_that("ci_sm() names conf_level when the level leaves tails below the floor", {
  expect_floor_refusal(ci_sm(sm = .5, N = 50, conf_level = 1 - 1e-12), "conf_level")
  expect_floor_refusal(ci_sm(sm = -.5, N = 30, conf_level = NULL, alpha_lower = .025,
                             alpha_upper = 5e-10), "alpha_upper")
})

test_that("ci_sc() with one tail given asks for both before it checks the floor", {
  e <- expect_error(ci_sc(means = c(1, 0), s_anova = 1, c_weights = c(1, -1),
                          n_per_group = c(20, 20), conf_level = NULL,
                          alpha_upper = 1e-20),
                    class = "dmar_argument_invalid")
  expect_null(e$reason)
  expect_match(conditionMessage(e), "both 'alpha_lower' and 'alpha_upper' must be specified",
               fixed = TRUE)
})

test_that("ss_aipe_smd() and ss_aipe_smd_sensitivity() refuse an assurance whose tail areas fall below the floor in its own name", {
  # The floor applies to (1 - assurance) / 2, and the refusal names
  # 'assurance', the argument given.
  a <- 1 - 1e-10
  for (m in list(NULL, "analytic")) {
    e <- expect_error(do.call(ss_aipe_smd, c(list(true_smd = .5, width = .5,
                                                 assurance = a),
                                            assurance_method = m)),
                      class = "dmar_argument_invalid")
    expect_identical(e$arg, "assurance")
    expect_identical(e$reason, "tail_below_floor")
  }
  e <- expect_error(ss_aipe_smd_sensitivity(true_smd = .5, supposed_smd = .5, width = .5,
                                            assurance = a, G = 5, seed = 113),
                    class = "dmar_argument_invalid")
  expect_identical(e$arg, "assurance")
  # The largest assurance accepted is planned.
  expect_true(is.finite(ss_aipe_smd(true_smd = .5, width = .5,
                                    assurance = .999999998)$value[1]))
  expect_floor_refusal(ss_aipe_smd(true_smd = .5, width = .5, conf_level = 1 - 1e-12),
                       "conf_level")
})

test_that("ss_aipe_c(), ss_aipe_r(), and their sensitivity functions refuse a level that leaves tails below the floor, as ci_c() and ci_r() do", {
  expect_floor_refusal(ss_aipe_c(true_var_error = 1, c_weights = c(1, -1), width = .5,
                                 conf_level = 1 - 1e-12), "conf_level")
  expect_floor_refusal(ss_aipe_r(true_cor = .3, width = .2, conf_level = 1 - 1e-12),
                       "conf_level")
  expect_floor_refusal(ss_aipe_r_sensitivity(true_cor = .3, specified_N = 50, width = .2,
                                             conf_level = 1 - 1e-12, G = 5, seed = 113),
                       "conf_level")
  expect_floor_refusal(ss_aipe_r_sensitivity(true_cor = .3, supposed_cor = .3, width = .2,
                                             conf_level = 1 - 1e-12, G = 5, seed = 113),
                       "conf_level")
  expect_floor_refusal(ss_aipe_c_sensitivity(true_var_error = 1, supposed_var_error = 1,
                                             c_weights = c(1, -1), width = .5,
                                             conf_level = 1 - 1e-12, G = 5, seed = 113),
                       "conf_level")
  # The largest level accepted is planned.
  expect_true(is.finite(ss_aipe_c(true_var_error = 1, c_weights = c(1, -1), width = .5,
                                  conf_level = .999999998)$value[1]))
  expect_true(is.finite(ss_aipe_r(true_cor = .3, width = .2,
                                  conf_level = .999999998)$value[1]))
})

test_that("the planners take their upper quantiles in the upper tail", {
  # qt(1 - 1e-9, df) and qnorm(1 - 1e-9) are about 7.7e-10 of themselves
  # off; each check below fails with them.
  expect_identical(DMAR:::.reg_coef_central_delta(47, 1e-9, .025)(0),
                   stats::qt(1e-9, 47, lower.tail = FALSE) +
                     stats::qt(.025, 47, lower.tail = FALSE))
  q <- stats::qt(.025, 46, lower.tail = FALSE) + stats::qt(1e-9, 46, lower.tail = FALSE)
  K <- 1.5^2 * (1 - .3) / (1 - .2)
  expect_identical(DMAR:::.reg_coef_raw_width_probability(50, 3, .3, .2, 1.5, 2, 1e-9, .025),
                   stats::pf(2^2 * 47 / (q^2 * K), 46, 47))
  # Near 1e12 per group the error of 1 - 1e-9 moves a plan by about a
  # thousand: at a contrast of 0, ss_aipe_sc() and ss_aipe_sc_ancova()
  # returned plans above the smallest n, and the central plan of
  # ss_aipe_reg_coef() was 508 above the smallest N. (ss_aipe_sc_ancova()
  # takes equal tails only; its search is called here with the two tails.)
  w <- 1.2e-5
  z <- stats::qnorm(1e-9, lower.tail = FALSE) + stats::qnorm(.025, lower.tail = FALSE)
  plan <- ss_aipe_sc(true_std_contrast = 0, c_weights = c(1, -1), width = w,
                     conf_level = NULL, alpha_lower = 1e-9, alpha_upper = .025)$value
  expect_lte(sqrt(2 / plan) * z, w)
  expect_gt(sqrt(2 / (plan - 1)) * z, w)
  plan <- DMAR:::.sc_ancova_plan_width(0, 1, c(1, -1), w, 1e-9, .025)
  expect_lte(sqrt(2 / plan) * z, w)
  expect_gt(sqrt(2 / (plan - 1)) * z, w)
  width_at <- function(N) {
    (stats::qt(1e-9, N - 2, lower.tail = FALSE) +
       stats::qt(.025, N - 2, lower.tail = FALSE)) / sqrt(N - 1)
  }
  plan <- ss_aipe_reg_coef(true_R2_Y_X = 0, true_R2_j_X_without_j = 0,
                           p = 1, true_reg_coef_j = 0, width = w,
                           conf_level = NULL, alpha_lower = 1e-9,
                           alpha_upper = .025)$value
  expect_lte(width_at(plan), w)
  expect_gt(width_at(plan - 1), w)
})
