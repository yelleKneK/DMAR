# Tests for vargha_delaney_A() against the Mann-Whitney U formulation and
# the DeLong et al. (1988) variance.

test_that("vargha_delaney_A() recovers A from the rank/U formulation", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(3, 4, 5, 6, 7)
  res <- vargha_delaney_A(x, y)
  # Hand calculation: sum of indicators (ties at half) is 4.5 over 25 pairs.
  expect_equal(res$A, 4.5 / 25, tolerance = 1e-10)
})

test_that("vargha_delaney_A() equals 0.5 for identical samples", {
  set.seed(113)
  x <- stats::rnorm(40)
  res <- vargha_delaney_A(x, x)
  expect_equal(res$A, 0.5, tolerance = 1e-10)
})

test_that("vargha_delaney_A() and 1 - A are symmetric across group order", {
  set.seed(113)
  x <- stats::rnorm(40, mean = 0.6)
  y <- stats::rnorm(40, mean = 0)
  res_xy <- vargha_delaney_A(x, y)
  res_yx <- vargha_delaney_A(y, x)
  expect_equal(res_xy$A, 1 - res_yx$A, tolerance = 1e-10)
  expect_equal(res_xy$se, res_yx$se,   tolerance = 1e-10)
})

test_that("vargha_delaney_A() matches the Mann-Whitney U / (n_1 n_2) computation", {
  set.seed(113)
  x <- stats::rnorm(50, mean = 0.5)
  y <- stats::rnorm(50)
  res <- vargha_delaney_A(x, y)
  # Independent computation: U_x via wilcox.test (R reports U for x as the
  # MWU "statistic" when 'x' is the first argument).
  U_x <- as.numeric(stats::wilcox.test(x, y, exact = FALSE)$statistic)
  expect_equal(res$A, U_x / (length(x) * length(y)), tolerance = 1e-8)
})

test_that("vargha_delaney_A() DeLong variance reproduces a hand calculation on a small sample", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(3, 4, 5, 6, 7)
  res <- vargha_delaney_A(x, y, ci_method = "wald")
  # V_10 placements per X_i:
  #   X=1: 0; X=2: 0; X=3: 0.1; X=4: 0.3; X=5: 0.5
  V_10 <- c(0, 0, 0.1, 0.3, 0.5)
  V_01 <- c(0.5, 0.3, 0.1, 0, 0)
  var_check <- stats::var(V_10) / 5 + stats::var(V_01) / 5
  expect_equal(res$se^2, var_check, tolerance = 1e-10)
})

test_that("vargha_delaney_A() formula interface matches two-vector interface", {
  set.seed(113)
  res_form <- vargha_delaney_A(len ~ supp, data = ToothGrowth)
  x_oj <- ToothGrowth$len[ToothGrowth$supp == "OJ"]
  x_vc <- ToothGrowth$len[ToothGrowth$supp == "VC"]
  res_vec <- vargha_delaney_A(x_oj, x_vc)
  expect_equal(res_form$A,  res_vec$A,  tolerance = 1e-10)
  expect_equal(res_form$se, res_vec$se, tolerance = 1e-10)
})

test_that("vargha_delaney_A() CI lies in [0, 1] and brackets the point estimate", {
  set.seed(113)
  x <- stats::rnorm(40, mean = 0.6)
  y <- stats::rnorm(40, mean = 0)
  res <- vargha_delaney_A(x, y)
  expect_gte(res$ci_lower, 0)
  expect_lte(res$ci_upper, 1)
  expect_lt(res$ci_lower, res$A)
  expect_gt(res$ci_upper, res$A)
})

test_that("vargha_delaney_A() wald CI agrees with normal approximation", {
  set.seed(113)
  x <- stats::rnorm(40, mean = 0.6)
  y <- stats::rnorm(40, mean = 0)
  res <- vargha_delaney_A(x, y, ci_method = "wald")
  z   <- stats::qnorm(0.975)
  expect_equal(res$ci_lower, max(0, res$A - z * res$se),
               tolerance = 1e-10)
  expect_equal(res$ci_upper, min(1, res$A + z * res$se),
               tolerance = 1e-10)
})

test_that("vargha_delaney_A() is a DMAR table that records its confidence level", {
  # The 1.1.0 review: a plain data.frame that recorded no level, so
  # results_sentence() wrote a 90% interval as "95% CI" (the 95% interval
  # here is [0.543, 0.782]).
  set.seed(1)
  u <- rnorm(40, .5)
  v <- rnorm(40)
  a90 <- vargha_delaney_A(u, v, conf_level = 0.90)
  expect_s3_class(a90, c("dmar_vargha_delaney_A", "dmar_tbl", "data.frame"),
                  exact = TRUE)
  expect_identical(attr(a90, "conf_level"), 0.90)
  # The sentence reports the p-value of the test that was run: Student's t
  # by default from 1.1.0, and the 1.0.x normal reference (p = 0.0051,
  # which these two sentences pinned before the default changed) under
  # test = "z".
  expect_identical(results_sentence(a90, digits = 3),
                   "A = 0.674, 90% CI [0.565, 0.767], unadjusted p = 0.0064")
  expect_identical(results_sentence(vargha_delaney_A(u, v), digits = 3),
                   "A = 0.674, 95% CI [0.543, 0.782], unadjusted p = 0.0064")
  expect_identical(
    results_sentence(vargha_delaney_A(u, v, conf_level = 0.90, test = "z"),
                     digits = 3),
    "A = 0.674, 90% CI [0.565, 0.767], unadjusted p = 0.0051")
  expect_identical(
    results_sentence(vargha_delaney_A(u, v, test = "z"), digits = 3),
    "A = 0.674, 95% CI [0.543, 0.782], unadjusted p = 0.0051")
  expect_true("Confidence level: 90%" %in% utils::capture.output(print(a90)))
  # The first nine columns are those of 1.0.x, in their order; df and test
  # follow them.
  expect_identical(names(a90), c("A", "se", "ci_lower", "ci_upper",
                                 "test_statistic", "p_value", "n_1", "n_2",
                                 "ci_method", "df", "test"))
  # tidy() names the row for the estimate, not for the method label, and
  # gives the same sentence.
  td <- generics::tidy(a90)
  expect_identical(td$term, "A")
  expect_identical(td$conf_level, 0.90)
  expect_equal(td$estimate, a90$A, tolerance = 0)
  expect_equal(c(td$ci_lower, td$ci_upper),
               c(a90$ci_lower, a90$ci_upper), tolerance = 0)
  expect_identical(results_sentence(td, digits = 3),
                   results_sentence(a90, digits = 3))
  expect_identical(generics::glance(a90), td)
})

# ---- The test of A = 1/2 (1.1.0: test = "t", "permutation", "z") ---------

test_that("test = \"t\" reproduces two published Brunner-Munzel results", {
  # The pain score data in the help of lawstat::brunner.munzel.test():
  # statistic 3.1375, df 17.683, p 0.005786, estimate 0.788961.
  pain_n <- c(3, 3, 4, 3, 1, 2, 3, 1, 1, 5, 4)
  pain_y <- c(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1)
  res <- vargha_delaney_A(pain_n, pain_y, test = "t")
  expect_equal(res$test_statistic, 3.137467, tolerance = 1e-6)
  expect_equal(res$df, 17.682842, tolerance = 1e-7)
  expect_equal(res$p_value, 0.005786, tolerance = 1e-4)
  expect_equal(res$A, 0.788961, tolerance = 1e-6)
  expect_identical(res$test, "t")
  # Table 2 of Schuurhuis, Konietschke, and Brunner (2025): estimate 0.8,
  # p = 0.0239.
  g1 <- c(1956, 3828, 2051, 3721, 3233, 2000, 4000, 4428, 2603, 2370)
  g2 <- c(820, 3364, 1957, 1851, 2984, 744, 2044)
  res <- vargha_delaney_A(g1, g2)
  expect_equal(res$A, 0.8, tolerance = 1e-12)
  expect_equal(res$df, 11.2192, tolerance = 1e-5)
  expect_equal(round(res$p_value, 4), 0.0239)
  # The degrees of freedom are Satterthwaite's, from the two placement
  # variances.
  psi <- outer(g1, g2, ">") + 0.5 * outer(g1, g2, "==")
  v_1 <- stats::var(rowMeans(psi)) / 10
  v_2 <- stats::var(colMeans(psi)) / 7
  expect_equal(res$df, (v_1 + v_2)^2 / (v_1^2 / 9 + v_2^2 / 6),
               tolerance = 1e-12)
  expect_equal(res$p_value, 2 * stats::pt(-abs(res$test_statistic), res$df),
               tolerance = 1e-12)
})

test_that("the three tests share the estimate, standard error, statistic, and interval", {
  x <- c(3.5, 5, 6, 7, 8)
  y <- c(1, 2, 3, 4, 4.2)
  t_ref <- vargha_delaney_A(x, y, test = "t")
  z_ref <- vargha_delaney_A(x, y, test = "z")
  perm  <- vargha_delaney_A(x, y, test = "permutation")
  same <- c("A", "se", "ci_lower", "ci_upper", "test_statistic", "n_1", "n_2",
            "ci_method")
  expect_identical(as.list(z_ref[same]), as.list(t_ref[same]))
  expect_identical(as.list(perm[same]), as.list(t_ref[same]))
  # test = "z" is the 1.0.x p-value.
  expect_equal(z_ref$p_value, 2 * stats::pnorm(-abs(z_ref$test_statistic)),
               tolerance = 1e-14)
  expect_true(is.na(z_ref$df) && is.na(perm$df))
  expect_identical(c(t_ref$test, perm$test, z_ref$test),
                   c("t", "permutation", "z"))
  # The worked example of the help: A = .92, logit interval [.486, .993],
  # t test p = .0033, exact permutation p = 8/252.
  expect_equal(t_ref$A, 0.92)
  expect_equal(round(c(t_ref$ci_lower, t_ref$ci_upper), 3),
               c(0.486, 0.993))
  expect_equal(round(t_ref$p_value, 4), 0.0033)
  expect_equal(perm$p_value, 8 / 252, tolerance = 1e-14)
  expect_gt(t_ref$p_value, z_ref$p_value)
  # Only the permutation test records its reassignments, and prints them.
  expect_null(attr(t_ref, "resampling"))
  expect_identical(attr(perm, "resampling"),
                   list(exact = TRUE, n_arrangements = 252,
                        n_resamples_used = 252L))
  expect_true("Permutation test: all 252 reassignments" %in%
                utils::capture.output(print(perm)))
  # tidy() and glance() give df and test after p_value.
  td <- generics::tidy(t_ref)
  expect_identical(names(td)[5:7], c("p_value", "df", "test"))
  expect_equal(td$df, t_ref$df)
  expect_identical(td$interval_degenerate, FALSE)
  expect_identical(generics::glance(t_ref), td)
})

test_that("the permutation p-value equals a brute-force enumeration, ties included", {
  W_reference <- function(x, y) {
    psi <- outer(x, y, function(a, b) (a > b) + 0.5 * (a == b))
    A <- mean(psi)
    v <- stats::var(rowMeans(psi)) / length(x) +
      stats::var(colMeans(psi)) / length(y)
    if (v < 1e-14) {
      if (abs(A - 0.5) < 1e-12) 0 else sign(A - 0.5) * Inf
    } else {
      (A - 0.5) / sqrt(v)
    }
  }
  cases <- list(list(c(1, 2, 2, 3), c(2, 3, 3, 4, 5)),
                list(c(0.3, 1.9, 2.2, 4.1, 5), c(-1, 0.2, 0.4)),
                list(c(1, 1, 0, 1), c(0, 0, 1, 0, 0)),
                list(c(7, 8, 9), c(1, 2, 3, 4)),
                list(c(2, 2, 2), c(2, 2, 2)))
  for (case in cases) {
    x <- case[[1L]]
    y <- case[[2L]]
    pooled <- c(x, y)
    W <- W_reference(x, y)
    alloc <- utils::combn(length(pooled), length(x))
    W_b <- apply(alloc, 2L, function(i) W_reference(pooled[i], pooled[-i]))
    got <- suppressWarnings(vargha_delaney_A(x, y, test = "permutation"))
    expect_equal(got$p_value, mean(abs(W_b) >= abs(W) - 1e-9),
                 tolerance = 1e-14)
    if (is.finite(W)) {
      expect_equal(got$test_statistic, W, tolerance = 1e-12)
    } else {
      expect_identical(got$test_statistic, W)
    }
  }
})

test_that("with an estimated variance of 0 only the permutation test gives a p-value", {
  for (test in c("t", "z")) {
    res <- suppressWarnings(vargha_delaney_A(6:10, 1:5, test = test))
    expect_true(is.na(res$test_statistic) && is.na(res$p_value) && is.na(res$df))
  }
  expect_warning(up <- vargha_delaney_A(6:10, 1:5, test = "permutation"),
                 class = "dmar_interval_degenerate")
  expect_identical(up$test_statistic, Inf)
  expect_equal(up$p_value, 2 / 252, tolerance = 1e-14)
  down <- suppressWarnings(vargha_delaney_A(1:5, 6:10, test = "permutation"))
  expect_identical(down$test_statistic, -Inf)
  expect_equal(down$p_value, 2 / 252, tolerance = 1e-14)
  flat <- suppressWarnings(vargha_delaney_A(c(2, 2, 2), c(2, 2, 2),
                                            test = "permutation"))
  expect_identical(c(flat$test_statistic, flat$p_value), c(0, 1))
  # The warning says which tests are missing.
  w <- tryCatch(vargha_delaney_A(6:10, 1:5),
                dmar_interval_degenerate = function(w) w)
  expect_match(conditionMessage(w),
               "No test is reported (test_statistic and p_value are NA).",
               fixed = TRUE)
  w <- tryCatch(vargha_delaney_A(6:10, 1:5, test = "permutation"),
                dmar_interval_degenerate = function(w) w)
  expect_match(conditionMessage(w), "The permutation test is reported",
               fixed = TRUE)
})

test_that("exact, n_resamples, and seed are read only under test = \"permutation\"", {
  skip_on_cran()
  x <- c(1, 4, 2, 8, 5, 7, 3)
  y <- c(2, 9, 6, 10, 12, 11)
  base <- vargha_delaney_A(x, y, test = "t")
  expect_identical(vargha_delaney_A(x, y, test = "t", n_resamples = 4.5,
                                    exact = "yes", seed = "a"), base)
  expect_identical(vargha_delaney_A(x, y, test = "z", n_resamples = 4.5)$A,
                   base$A)
  expect_error(vargha_delaney_A(x, y, test = "permutation",
                                n_resamples = 4.5),
               class = "dmar_count_invalid")
  expect_error(vargha_delaney_A(x, y, test = "permutation", exact = "yes"),
               "'exact' must be NULL, TRUE, or FALSE.", fixed = TRUE)
  expect_error(vargha_delaney_A(x, y, test = "exact"), "should be one of")

  # Random reassignments: reproducible with a seed, which leaves the
  # caller's random number state as it was, and of the documented form.
  set.seed(42)
  before <- .Random.seed
  first  <- vargha_delaney_A(x, y, test = "permutation", exact = FALSE,
                             n_resamples = 500, seed = 113)
  expect_identical(.Random.seed, before)
  second <- vargha_delaney_A(x, y, test = "permutation", exact = FALSE,
                             n_resamples = 500, seed = 113)
  expect_identical(first, second)
  expect_identical(attr(first, "resampling"),
                   list(exact = FALSE, n_arrangements = choose(13, 7),
                        n_resamples_used = 500L))
  expect_true("Permutation test: 500 random reassignments" %in%
                utils::capture.output(print(first)))
  expect_equal((first$p_value * 501 - 1) %% 1, 0, tolerance = 1e-9)
  # Without a seed the draws come from the caller's stream and advance it.
  set.seed(7)
  unseeded <- vargha_delaney_A(x, y, test = "permutation", exact = FALSE,
                               n_resamples = 500)
  expect_false(identical(.Random.seed, {set.seed(7); .Random.seed}))
  # The exact p-value is within Monte Carlo error of the sampled one.
  exact <- vargha_delaney_A(x, y, test = "permutation")
  expect_identical(attr(exact, "resampling")$n_resamples_used, 1716L)
  expect_lt(abs(first$p_value - exact$p_value),
            4 * sqrt(exact$p_value * (1 - exact$p_value) / 500) + 2 / 501)
  # Enumeration is refused above the cap of randomization_test().
  expect_error(vargha_delaney_A(seq(1, 29, 2), seq(2, 30, 2),
                                test = "permutation", exact = TRUE),
               "Exact enumeration would require")
})

test_that("the t reference holds the level at 10 per group where the normal reference does not", {
  skip_on_cran()
  set.seed(20261110)
  reps <- 4000
  reject_t <- reject_z <- logical(reps)
  for (r in seq_len(reps)) {
    res <- vargha_delaney_A(stats::rnorm(10), stats::rnorm(10), test = "t")
    reject_t[r] <- isTRUE(res$p_value < 0.05)
    reject_z[r] <- isTRUE(2 * stats::pnorm(-abs(res$test_statistic)) < 0.05)
  }
  # Monte Carlo standard error of a rate of .05 over 4,000 samples: .0034.
  expect_gt(mean(reject_t), 0.035)
  expect_lt(mean(reject_t), 0.07)
  expect_gt(mean(reject_z), mean(reject_t))
})

test_that("vargha_delaney_A() gives no warning with the default test", {
  expect_no_warning(vargha_delaney_A(c(3.5, 5, 6, 7, 8), c(1, 2, 3, 4, 4.2)))
})
