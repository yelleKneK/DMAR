# ss_power_R2_change(): the model comparison F test of predictors added to a
# regression model. Fixed predictors: noncentral F with noncentrality N f^2.
# Random predictors: the distribution of the sample squared multiple partial
# correlation (Gatsonis & Sampson, 1989, p. 518) through Lee's (1971)
# approximation.

row_value <- function(x, term) x$value[x$term == term]
power_change <- function(...) row_value(ss_power_R2_change(...), "actual_power")
n_change <- function(...) row_value(ss_power_R2_change(...), "necessary_N")

test_that("fixed-predictor power reproduces the G*Power 3.1 manual, Section 14.3.1", {
  # Five predictors explain .25, four added predictors raise that to .30,
  # alpha .01, N = 90. The manual prints power 0.241297, noncentrality
  # 6.428574, critical F 3.563110 on 4 and 80 degrees of freedom, and a
  # sample size of 242 for power .80.
  res <- ss_power_R2_change(true_R2_restricted = 0.25, true_R2_full = 0.30,
                            p_restricted = 5, p_full = 9, alpha_level = 0.01,
                            specified_N = 90, random_predictors = FALSE)
  expect_equal(row_value(res, "actual_power"), 0.241296492, tolerance = 1e-8)
  expect_equal(row_value(res, "noncentral_f_parm"), 6.42857143, tolerance = 1e-8)
  expect_equal(row_value(res, "df_1"), 4)
  expect_equal(row_value(res, "df_2"), 80)
  expect_equal(qf(0.99, 4, 80), 3.56310963, tolerance = 1e-8)
  expect_equal(row_value(res, "true_f2"), 0.05 / 0.70)
  expect_equal(row_value(res, "true_partial_R2"), 0.05 / 0.75)
  expect_equal(row_value(res, "true_R2_change"), 0.05)

  found <- ss_power_R2_change(true_R2_restricted = 0.25, true_R2_full = 0.30,
                              p_restricted = 5, p_full = 9, alpha_level = 0.01,
                              desired_power = 0.80, random_predictors = FALSE)
  expect_equal(row_value(found, "necessary_N"), 242)
  expect_gte(row_value(found, "actual_power"), 0.80)
  expect_equal(row_value(found, "df_2"), 232)
  expect_equal(power_change(true_R2_restricted = 0.25, true_R2_full = 0.30,
                            p_restricted = 5, p_full = 9, alpha_level = 0.01,
                            specified_N = 241, random_predictors = FALSE),
               0.799252437, tolerance = 1e-8)
})

test_that("error from a larger model reproduces the G*Power 3.1 manual, Section 14.3.2", {
  # f^2 = .075 with three predictors tested and twelve in the model that
  # supplies the error, alpha .01, N = 200: the manual prints power 0.766990,
  # noncentrality 15, and 187 denominator degrees of freedom.
  res <- ss_power_R2_change(true_f2 = 0.075, p_full = 12, p_restricted = 9,
                            alpha_level = 0.01, specified_N = 200,
                            random_predictors = FALSE)
  expect_equal(row_value(res, "actual_power"), 0.766989564, tolerance = 1e-8)
  expect_equal(row_value(res, "noncentral_f_parm"), 15)
  expect_equal(row_value(res, "df_1"), 3)
  expect_equal(row_value(res, "df_2"), 187)
})

# On CRAN the file stops here, after the anchors above; the blocks below run
# in the full local suite (see cran_files in tests/testthat.R).
skip_on_cran()

test_that("Cohen's noncentrality f^2 (N - w) is the call with no restricted model at N - w", {
  # The design of the first test under Cohen's noncentrality: u = 4, w = 5,
  # N = 90, so f^2 (N - w) = (.05 / .70) * 85.
  res <- ss_power_R2_change(true_f2 = 0.05 / 0.70, p_restricted = 0, p_full = 4,
                            alpha_level = 0.01, specified_N = 85,
                            random_predictors = FALSE)
  expect_equal(row_value(res, "noncentral_f_parm"), 6.07142857, tolerance = 1e-8)
  expect_equal(row_value(res, "actual_power"), 0.222508251, tolerance = 1e-8)
  expect_equal(1 - pf(qf(0.99, 4, 80), 4, 80, ncp = (0.05 / 0.70) * (90 - 5)),
               0.222508251, tolerance = 1e-8)
  # The design of the second test: u = 3, w = 9, N = 200.
  res2 <- ss_power_R2_change(true_f2 = 0.075, p_restricted = 0, p_full = 3,
                             alpha_level = 0.01, specified_N = 191,
                             random_predictors = FALSE)
  expect_equal(row_value(res2, "noncentral_f_parm"), 14.325)
  expect_equal(row_value(res2, "actual_power"), 0.740721888, tolerance = 1e-8)
  # The sample size by Cohen's noncentrality is at most w above the
  # fixed-predictor answer: 247 against 242 here.
  cohen_N <- n_change(true_f2 = 0.05 / 0.70, p_restricted = 0, p_full = 4,
                      alpha_level = 0.01, desired_power = 0.80,
                      random_predictors = FALSE) + 5
  expect_equal(cohen_N, 247)
})

test_that("one added predictor with fixed predictors is ss_power_reg_coef()", {
  grid <- expand.grid(R2_full = c(0.2, 0.4, 0.6), change = c(0.02, 0.05, 0.10),
                      p = c(2, 5, 10), N = c(30, 100, 400))
  expect_equal(nrow(grid), 81)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    ours <- power_change(true_R2_full = g$R2_full, true_R2_change = g$change,
                         p_full = g$p, p_restricted = g$p - 1,
                         specified_N = g$N, random_predictors = FALSE)
    theirs <- ss_power_reg_coef(true_R2_Y_X = g$R2_full,
                                true_R2_Y_X_without_j = g$R2_full - g$change,
                                p = g$p, specified_N = g$N)
    expect_equal(ours, row_value(theirs, "actual_power"), tolerance = 1e-8)
  }
  search <- expand.grid(R2_full = c(0.2, 0.4, 0.6), change = c(0.02, 0.05, 0.10))
  for (i in seq_len(nrow(search))) {
    g <- search[i, ]
    ours <- n_change(true_R2_full = g$R2_full, true_R2_change = g$change,
                     p_full = 5, p_restricted = 4, desired_power = 0.80,
                     random_predictors = FALSE)
    theirs <- ss_power_reg_coef(true_R2_Y_X = g$R2_full,
                                true_R2_Y_X_without_j = g$R2_full - g$change,
                                p = 5, desired_power = 0.80)
    expect_equal(ours, row_value(theirs, "necessary_N"))
  }
  expect_equal(n_change(true_R2_restricted = 0.30, true_R2_full = 0.40,
                        p_restricted = 4, p_full = 5, desired_power = 0.80,
                        random_predictors = FALSE), 50)
})

test_that("no predictors in the restricted model is ss_power_R2()", {
  grid <- expand.grid(R2 = c(0.1, 0.3, 0.6), p = c(1, 5, 12), N = c(20, 60, 200),
                      rp = c(TRUE, FALSE))
  expect_equal(nrow(grid), 54)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    ours <- ss_power_R2_change(true_R2_full = g$R2, true_R2_restricted = 0,
                               p_full = g$p, p_restricted = 0,
                               specified_N = g$N, random_predictors = g$rp)
    omni <- ss_power_R2(true_R2 = g$R2, p = g$p, specified_N = g$N,
                        random_predictors = g$rp)
    expect_equal(row_value(ours, "actual_power"),
                 row_value(omni, "actual_power"), tolerance = 1e-12)
    expect_equal(row_value(ours, "true_f2"),
                 row_value(omni, "true_f2"), tolerance = 1e-12)
    expect_equal(row_value(ours, "noncentral_f_parm"),
                 row_value(omni, "noncentral_f_parm"), tolerance = 1e-12)
  }
  for (rp in c(TRUE, FALSE)) {
    expect_equal(n_change(true_R2_full = 0.3, true_R2_change = 0.3, p_full = 5,
                          p_restricted = 0, desired_power = 0.80,
                          random_predictors = rp),
                 row_value(ss_power_R2(true_R2 = 0.3, p = 5,
                                       desired_power = 0.80,
                                       random_predictors = rp), "necessary_N"))
  }
})

test_that("random-predictor power matches the exact distribution of the squared multiple partial correlation", {
  # Under multivariate normality the sample squared multiple partial
  # correlation of the outcome with the u added predictors given the w others
  # is a negative binomial mixture of beta variates: given K it is
  # Beta(u/2 + K, (N - p_full - 1)/2), with K negative binomial (size
  # (N - 1 - w)/2, probability 1 - the population squared partial). The
  # critical value is that of the central F test.
  exact_power <- function(N, p_full, p_restricted, f2, alpha) {
    u <- p_full - p_restricted
    df_2 <- N - p_full - 1
    f_crit <- qf(1 - alpha, u, df_2)
    R2_crit <- u * f_crit / (u * f_crit + df_2)
    pR2 <- f2 / (1 + f2)
    size <- (N - 1 - p_restricted) / 2
    k <- 0:qnbinom(1e-17, size, 1 - pR2, lower.tail = FALSE)
    sum(dnbinom(k, size, 1 - pR2) *
          pbeta(R2_crit, u / 2 + k, df_2 / 2, lower.tail = FALSE))
  }
  cells <- data.frame(
    N            = c(60, 40, 31, 200, 15, 90),
    p_full       = c(5, 13, 26, 8, 3, 9),
    p_restricted = c(3, 10, 25, 2, 2, 5),
    f2           = c(0.10 / 0.70, 0.15 / 0.55, 1.5, 0.125, 0.25, 0.05 / 0.70),
    alpha        = c(0.05, 0.05, 0.05, 0.01, 0.05, 0.01),
    exact        = c(0.6825505, 0.5659049, 0.5014549, 0.8997967, 0.3477531,
                     0.2201927))
  for (i in seq_len(nrow(cells))) {
    s <- cells[i, ]
    exact <- exact_power(s$N, s$p_full, s$p_restricted, s$f2, s$alpha)
    expect_equal(exact, s$exact, tolerance = 1e-6)
    expect_lt(abs(power_change(true_f2 = s$f2, p_full = s$p_full,
                               p_restricted = s$p_restricted,
                               alpha_level = s$alpha, specified_N = s$N) -
                    exact), 1.5e-3)
  }
  # The largest difference found in 1,500 further combinations: three error
  # degrees of freedom, a squared partial correlation of .95, alpha .10.
  widest <- exact_power(55, 51, 50, 19, 0.10)
  expect_equal(widest, 0.97796, tolerance = 1e-4)
  expect_lt(abs(power_change(true_f2 = 19, p_full = 51, p_restricted = 50,
                             alpha_level = 0.10, specified_N = 55) - widest),
            4e-3)
  # The sample size for the last of the six designs at power .80 is the one
  # the exact distribution gives.
  expect_equal(n_change(true_f2 = 0.05 / 0.70, p_full = 9, p_restricted = 5,
                        alpha_level = 0.01, desired_power = 0.80), 250)
  expect_lt(exact_power(249, 9, 5, 0.05 / 0.70, 0.01), 0.80)
  expect_gte(exact_power(250, 9, 5, 0.05 / 0.70, 0.01), 0.80)
  expect_lt(power_change(true_f2 = 0.05 / 0.70, p_full = 9, p_restricted = 5,
                         alpha_level = 0.01, specified_N = 249), 0.80)
})

test_that("random predictors give lower power than fixed at the power used in planning", {
  # Not a universal ordering: at very low power the order can reverse by a
  # small amount. The six designs here have fixed-predictor power of at
  # least .50.
  designs <- data.frame(N            = c(60, 50, 200, 40, 25, 120),
                        p_full       = c(5, 5, 8, 13, 4, 6),
                        p_restricted = c(3, 4, 2, 10, 1, 2),
                        f2           = c(0.10 / 0.70, 0.10 / 0.60, 0.125,
                                         0.15 / 0.55, 0.35 / 0.60, 0.10))
  for (i in seq_len(nrow(designs))) {
    s <- designs[i, ]
    args <- list(true_f2 = s$f2, p_full = s$p_full,
                 p_restricted = s$p_restricted)
    fixed  <- do.call(power_change, c(args, specified_N = s$N,
                                      random_predictors = FALSE))
    random <- do.call(power_change, c(args, specified_N = s$N,
                                      random_predictors = TRUE))
    expect_gte(fixed, 0.50)
    expect_lt(random, fixed)
    expect_gte(do.call(n_change, c(args, desired_power = 0.80,
                                   random_predictors = TRUE)),
               do.call(n_change, c(args, desired_power = 0.80,
                                   random_predictors = FALSE)))
  }
  expect_equal(n_change(true_R2_restricted = 0.30, true_R2_full = 0.40,
                        p_restricted = 4, p_full = 5, desired_power = 0.80), 56)
  expect_equal(n_change(true_R2_restricted = 0.30, true_R2_full = 0.40,
                        p_restricted = 4, p_full = 5, desired_power = 0.80,
                        random_predictors = FALSE), 50)
})

test_that("the four ways of giving the effect return the same result", {
  for (rp in c(TRUE, FALSE)) {
    common <- list(p_full = 9, p_restricted = 5, alpha_level = 0.01,
                   specified_N = 90, random_predictors = rp)
    a <- do.call(ss_power_R2_change,
                 c(list(true_R2_full = 0.30, true_R2_restricted = 0.25), common))
    b <- do.call(ss_power_R2_change,
                 c(list(true_R2_full = 0.30, true_R2_change = 0.05), common))
    c2 <- do.call(ss_power_R2_change,
                  c(list(true_R2_restricted = 0.25, true_R2_change = 0.05), common))
    d <- do.call(ss_power_R2_change, c(list(true_f2 = 0.05 / 0.70), common))
    expect_equal(b$value, a$value, tolerance = 1e-12)
    expect_equal(c2$value, a$value, tolerance = 1e-12)
    keep <- a$term != "true_R2_change"
    expect_equal(d$value[keep], a$value[keep], tolerance = 1e-12)
    expect_true(is.na(row_value(d, "true_R2_change")))
    expect_identical(d$term, a$term)
  }
  # The same holds for a search.
  expect_equal(n_change(true_R2_restricted = 0.25, true_R2_change = 0.05,
                        p_restricted = 5, p_full = 9, alpha_level = 0.01,
                        desired_power = 0.80), 250)
})

test_that("power equals alpha_level when the change is zero, and the search stops", {
  for (rp in c(TRUE, FALSE)) {
    res <- ss_power_R2_change(true_R2_full = 0.2, true_R2_restricted = 0.2,
                              p_full = 6, p_restricted = 3, specified_N = 100,
                              random_predictors = rp)
    expect_equal(row_value(res, "actual_power"), 0.05, tolerance = 1e-10)
    expect_equal(row_value(res, "true_f2"), 0)
    expect_equal(power_change(true_f2 = 0, p_full = 6, p_restricted = 3,
                              alpha_level = 0.01, specified_N = 100,
                              random_predictors = rp), 0.01, tolerance = 1e-10)
    expect_error(
      ss_power_R2_change(true_R2_full = 0.2, true_R2_restricted = 0.2,
                         p_full = 6, p_restricted = 3, random_predictors = rp),
      paste0("The change in R-squared is 0, so the null hypothesis is true: ",
             "power equals 'alpha_level' at every sample size and no sample ",
             "size attains 'desired_power'."), fixed = TRUE)
    # The refusal is classed, so a caller can read it without its wording.
    expect_error(
      ss_power_R2_change(true_f2 = 0, p_full = 6, p_restricted = 3,
                         random_predictors = rp),
      class = "dmar_effect_zero")
  }
})

test_that("the search returns the smallest sample size that attains the power", {
  grid <- expand.grid(f2 = c(0.03, 0.15, 0.60), u = c(1, 4), w = c(2, 9),
                      rp = c(TRUE, FALSE))
  expect_equal(nrow(grid), 24)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    args <- list(true_f2 = g$f2, p_full = g$u + g$w, p_restricted = g$w,
                 random_predictors = g$rp)
    res <- do.call(ss_power_R2_change, c(args, desired_power = 0.90))
    N <- row_value(res, "necessary_N")
    expect_gte(row_value(res, "actual_power"), 0.90)
    expect_equal(row_value(res, "df_2"), N - (g$u + g$w) - 1)
    if (N > g$u + g$w + 2) {
      expect_lt(do.call(power_change, c(args, specified_N = N - 1)), 0.90)
    }
    expect_equal(row_value(res, "actual_power"),
                 do.call(power_change, c(args, specified_N = N)))
  }
  # A target the smallest admissible size already attains returns that size.
  for (rp in c(TRUE, FALSE)) {
    at_min <- power_change(true_f2 = 0.3, p_full = 5, p_restricted = 3,
                           specified_N = 7, random_predictors = rp)
    res <- ss_power_R2_change(true_f2 = 0.3, p_full = 5, p_restricted = 3,
                              desired_power = at_min - 1e-4,
                              random_predictors = rp)
    expect_equal(row_value(res, "necessary_N"), 7)
    expect_equal(row_value(res, "df_2"), 1)
  }
})

test_that("planning values equal up to the last bits are taken as equal", {
  # .1 + .2 exceeds .3 by 5.6e-17; neither is a restricted value above the
  # full one or a change above the full value.
  whole <- ss_power_R2_change(true_R2_full = 0.3, true_R2_change = 0.1 + 0.2,
                              p_full = 3, p_restricted = 0, specified_N = 50)
  swapped <- ss_power_R2_change(true_R2_full = 0.1 + 0.2, true_R2_change = 0.3,
                                p_full = 3, p_restricted = 0, specified_N = 50)
  omnibus <- ss_power_R2(true_R2 = 0.3, p = 3, specified_N = 50)
  expect_equal(whole$value[2], swapped$value[2], tolerance = 1e-12)
  expect_equal(whole$value[2], omnibus$value[2], tolerance = 1e-12)
  expect_equal(whole$value[2], 0.95687, tolerance = 1e-5)

  none <- ss_power_R2_change(true_R2_full = 0.3, true_R2_restricted = 0.1 + 0.2,
                             p_full = 5, p_restricted = 3, specified_N = 50)
  expect_identical(none$value[none$term == "true_R2_change"], 0)
  expect_equal(none$value[none$term == "actual_power"], 0.05, tolerance = 1e-12)

  # A real excess is still refused.
  expect_error(
    ss_power_R2_change(true_R2_full = 0.3, true_R2_change = 0.3 + 1e-9,
                       p_full = 5, p_restricted = 3, specified_N = 50),
    "'true_R2_change' .* exceeds 'true_R2_full'")
  expect_error(
    ss_power_R2_change(true_R2_full = 0.3, true_R2_restricted = 0.3 + 1e-9,
                       p_full = 5, p_restricted = 3, specified_N = 50),
    "'true_R2_restricted' .* exceeds 'true_R2_full'")
})

test_that("a restricted value and a change that sum to exactly 1 are refused", {
  expect_error(
    ss_power_R2_change(true_R2_restricted = 0.5, true_R2_change = 0.5,
                       p_full = 5, p_restricted = 3, specified_N = 50),
    "must be below 1")
})

test_that("ss_power_R2_change() returns the nine rows, the class, and the attribute", {
  terms <- c("actual_power", "noncentral_f_parm", "true_f2", "true_partial_R2",
             "true_R2_change", "q", "df_1", "df_2")
  found <- ss_power_R2_change(true_R2_restricted = 0.25, true_R2_full = 0.30,
                              p_restricted = 5, p_full = 9, desired_power = 0.80)
  expect_identical(found$term, c("necessary_N", terms))
  expect_identical(class(found), c("dmar_ss_power", "dmar_tbl", "data.frame"))
  expect_named(found, c("term", "value"))
  expect_type(found$value, "double")
  expect_true(is.na(row_value(found, "noncentral_f_parm")))
  expect_true(attr(found, "random_predictors"))

  at_N <- ss_power_R2_change(true_R2_restricted = 0.25, true_R2_full = 0.30,
                             p_restricted = 5, p_full = 9, specified_N = 150,
                             random_predictors = FALSE)
  expect_identical(at_N$term, c("specified_N", terms))
  expect_false(attr(at_N, "random_predictors"))
  expect_equal(row_value(at_N, "noncentral_f_parm"), 150 * 0.05 / 0.70)
  expect_equal(row_value(at_N, "df_2"), 140)
  # q is the number of predictors removed from the full model to give the
  # restricted model, p_full - p_restricted, and equals df_1.
  expect_identical(row_value(found, "q"), 4)
  expect_identical(row_value(at_N, "q"), row_value(at_N, "df_1"))
  expect_identical(
    row_value(ss_power_R2_change(true_f2 = .1, p_full = 7, p_restricted = 0,
                                 specified_N = 60), "q"), 7)

  td <- generics::tidy(found)
  expect_equal(td$estimate, row_value(found, "necessary_N"))
  expect_equal(td$power, row_value(found, "actual_power"))
  gl <- generics::glance(at_N)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$estimate, 150)
  expect_true(all(c("true_f2", "true_partial_R2", "true_R2_change", "df_1",
                    "df_2") %in% names(gl)))
  expect_output(print(found), "necessary_N")
})

test_that("ss_power_R2_change() refuses invalid arguments with the reason", {
  ok <- function(...) {
    args <- utils::modifyList(
      list(true_R2_full = 0.30, true_R2_restricted = 0.25, p_full = 9,
           p_restricted = 5, specified_N = 90), list(...), keep.null = TRUE)
    do.call(ss_power_R2_change, args)
  }
  expect_s3_class(ok(), "dmar_ss_power")

  expect_error(ok(alpha_level = 0), "'alpha_level' must be a single numeric value in (0, 1).", fixed = TRUE)
  expect_error(ok(alpha_level = 1), "'alpha_level' must be")
  expect_error(ok(alpha_level = c(0.05, 0.01)), "'alpha_level' must be")

  expect_error(ok(p_full = 0), "'p_full' must be a single whole number of at least 1.", fixed = TRUE)
  expect_error(ok(p_full = 9.5), "'p_full' must be")
  expect_error(ss_power_R2_change(true_R2_full = 0.3, true_R2_restricted = 0.25,
                                  p_restricted = 5, specified_N = 90),
               "'p_full' must be")

  p_restricted_msg <- "'p_restricted' must be a single whole number from 0 to p_full - 1. The restricted model has fewer predictors than the full model; 0 is the model with an intercept only."
  expect_error(ok(p_restricted = 9), p_restricted_msg, fixed = TRUE)
  expect_error(ok(p_restricted = 5.5), p_restricted_msg, fixed = TRUE)
  expect_error(ok(p_restricted = -1), p_restricted_msg, fixed = TRUE)
  expect_error(ok(p_restricted = 12), p_restricted_msg, fixed = TRUE)
  expect_error(ss_power_R2_change(true_R2_full = 0.3, true_R2_restricted = 0.25,
                                  p_full = 9, specified_N = 90),
               p_restricted_msg, fixed = TRUE)

  expect_error(ok(random_predictors = NA), "'random_predictors' must be TRUE or FALSE.", fixed = TRUE)
  expect_error(ok(random_predictors = "yes"), "'random_predictors' must be")

  one_way <- "Give the effect in exactly one way: two of 'true_R2_full', 'true_R2_restricted', and 'true_R2_change', or 'true_f2' alone."
  # All three R-squared values, true_f2 with one or two of them, one value
  # alone, and nothing.
  expect_error(ok(true_R2_change = 0.05), one_way, fixed = TRUE)
  expect_error(ok(true_f2 = 0.1), one_way, fixed = TRUE)
  expect_error(ok(true_R2_restricted = NULL, true_f2 = 0.1), one_way, fixed = TRUE)
  expect_error(ok(true_R2_restricted = NULL), one_way, fixed = TRUE)
  expect_error(ok(true_R2_full = NULL, true_R2_restricted = NULL), one_way, fixed = TRUE)
  expect_error(ok(true_R2_full = NULL, true_R2_restricted = NULL,
                  true_R2_change = 0.05), one_way, fixed = TRUE)

  expect_error(ok(true_R2_full = 1), "'true_R2_full' must be a single numeric value in [0, 1).", fixed = TRUE)
  expect_error(ok(true_R2_full = -0.1), "'true_R2_full' must be")
  expect_error(ok(true_R2_restricted = NA_real_), "'true_R2_restricted' must be a single numeric value in [0, 1).", fixed = TRUE)
  expect_error(ok(true_R2_restricted = NULL, true_R2_change = c(0.05, 0.1)),
               "'true_R2_change' must be a single numeric value in [0, 1).", fixed = TRUE)

  expect_error(ok(true_R2_restricted = 0.35),
               "'true_R2_restricted' (0.35) exceeds 'true_R2_full' (0.3). The restricted model's predictors are a subset of the full model's, so its squared multiple correlation cannot be the larger one.",
               fixed = TRUE)
  expect_error(ok(true_R2_restricted = NULL, true_R2_change = 0.4),
               "'true_R2_change' (0.4) exceeds 'true_R2_full' (0.3). The change is the part of 'true_R2_full' that the added predictors contribute, so it cannot be the larger one.",
               fixed = TRUE)
  expect_error(ok(true_R2_full = NULL, true_R2_restricted = 0.6, true_R2_change = 0.5),
               "'true_R2_restricted' plus 'true_R2_change' is 'true_R2_full' and must be below 1; the values given sum to 1.1.",
               fixed = TRUE)

  expect_error(ok(true_R2_full = NULL, true_R2_restricted = NULL, true_f2 = -0.1),
               "'true_f2' must be a single nonnegative numeric value.", fixed = TRUE)
  expect_error(ok(true_R2_full = NULL, true_R2_restricted = NULL, true_f2 = Inf),
               "'true_f2' must be")

  intercept_only <- "With 'p_restricted' = 0 the restricted model has an intercept only and explains nothing, so 'true_R2_restricted' must be 0 and 'true_R2_change' must equal 'true_R2_full'."
  expect_error(ok(p_restricted = 0, true_R2_restricted = 0.1), intercept_only, fixed = TRUE)
  expect_error(ok(p_restricted = 0, true_R2_restricted = NULL, true_R2_change = 0.2),
               intercept_only, fixed = TRUE)
  expect_s3_class(ok(p_restricted = 0, true_R2_restricted = 0), "dmar_ss_power")
  expect_s3_class(ok(p_restricted = 0, true_R2_restricted = NULL,
                     true_R2_change = 0.30), "dmar_ss_power")

  power_msg <- "'desired_power' must be a single number in (0, 1)"
  expect_error(ok(specified_N = NULL, desired_power = 1), power_msg, fixed = TRUE)
  expect_error(ok(specified_N = NULL, desired_power = 0), power_msg, fixed = TRUE)
  expect_error(ok(specified_N = NULL, desired_power = 1.5), power_msg, fixed = TRUE)
  expect_error(ok(specified_N = NULL, desired_power = 1),
               class = "dmar_desired_power_invalid")
  # desired_power is not used when a sample size is given, and a value
  # outside (0, 1) is refused all the same.
  expect_error(ok(desired_power = 1.5), class = "dmar_desired_power_invalid")
  expect_s3_class(ok(desired_power = .5), "dmar_ss_power")

  # p_full + 1 leaves no error degrees of freedom; p_full + 2 leaves one.
  expect_error(ok(specified_N = 10), "'specified_N' must be a single whole number of at least 11")
  expect_error(ok(specified_N = 90.5), "'specified_N' must be a single whole number")
  expect_equal(row_value(ok(specified_N = 11), "df_2"), 1)
})

test_that("ss_power_R2_change() size and power match a raw-data Monte Carlo", {
  skip_on_cran()
  # Three predictors explain .20 and two added predictors raise that to .30;
  # N = 60, alpha .05. The predictors are equicorrelated (.3) and the added
  # predictors' coefficients differ in size and sign. 4,000 replications:
  # the Monte Carlo standard error is .0034 for the size and at most .0074
  # for the power.
  set.seed(113)
  N <- 60; w <- 3; u <- 2; p <- w + u; G <- 4000; alpha <- 0.05
  Sxx <- matrix(0.3, p, p); diag(Sxx) <- 1
  A <- seq_len(w); B <- w + seq_len(u)
  population <- function(R2_restricted, R2_full) {
    S_BA <- Sxx[B, B] - Sxx[B, A] %*% solve(Sxx[A, A], Sxx[A, B])
    gA <- rev(A)
    gA <- gA * sqrt(R2_restricted / drop(t(gA) %*% Sxx[A, A] %*% gA))
    gB <- c(1, -2)
    change <- R2_full - R2_restricted
    gB <- gB * sqrt(change / drop(t(gB) %*% S_BA %*% gB))
    beta <- c(gA - solve(Sxx[A, A], Sxx[A, B] %*% gB), gB)
    # The construction gives the two population R^2 values exactly.
    stopifnot(abs(drop(t(beta) %*% Sxx %*% beta) - R2_full) < 1e-10)
    list(beta = beta, sigma = sqrt(1 - R2_full))
  }
  R <- chol(Sxx)
  df_2 <- N - p - 1
  f_crit <- qf(1 - alpha, u, df_2)
  one <- rep(1, N)
  reject_rate <- function(pop, random) {
    if (random) {
      e_full <- e_restricted <- numeric(G)
      for (g in seq_len(G)) {
        X <- matrix(rnorm(N * p), N) %*% R
        y <- drop(X %*% pop$beta) + pop$sigma * rnorm(N)
        e_full[g] <- sum(qr.resid(qr(cbind(one, X)), y)^2)
        e_restricted[g] <- sum(qr.resid(qr(cbind(one, X[, A])), y)^2)
      }
    } else {
      # One design held fixed, standardized so that its covariance matrix
      # with divisor N is Sxx: the R^2 values of the design itself are then
      # the population values, and the noncentrality is N f^2.
      Z <- scale(matrix(rnorm(N * p), N), scale = FALSE)
      X <- Z %*% solve(chol(crossprod(Z) / N)) %*% R
      Y <- drop(X %*% pop$beta) + pop$sigma * matrix(rnorm(N * G), N)
      e_full <- colSums(qr.resid(qr(cbind(one, X)), Y)^2)
      e_restricted <- colSums(qr.resid(qr(cbind(one, X[, A])), Y)^2)
    }
    mean(((e_restricted - e_full) / u) / (e_full / df_2) > f_crit)
  }
  alternative <- population(0.20, 0.30)
  null <- population(0.20, 0.20)

  random_power <- power_change(true_R2_restricted = 0.20, true_R2_full = 0.30,
                               p_restricted = w, p_full = p, specified_N = N)
  fixed_power <- power_change(true_R2_restricted = 0.20, true_R2_full = 0.30,
                              p_restricted = w, p_full = p, specified_N = N,
                              random_predictors = FALSE)
  expect_equal(random_power, 0.6825, tolerance = 1e-3)
  expect_equal(fixed_power, 0.7236, tolerance = 1e-3)

  expect_lt(abs(reject_rate(null, random = TRUE) - alpha), 0.012)
  expect_lt(abs(reject_rate(alternative, random = TRUE) - random_power), 0.025)
  expect_lt(abs(reject_rate(null, random = FALSE) - alpha), 0.012)
  expect_lt(abs(reject_rate(alternative, random = FALSE) - fixed_power), 0.025)
})
