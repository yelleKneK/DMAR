# The oracle for the limits is a noncentral F interval computed by hand in
# this file: each limit of the noncentrality parameter solves
# pf(F, df_1, df_2, ncp) = 1 - alpha_lower or alpha_upper with
# stats::uniroot(), is 0 when the central F already puts F beyond that
# quantile, is infinite for a tail area of 0, and becomes a limit of the
# square root of the signal-to-noise ratio as sqrt(ncp / N). It does not
# call DMAR's noncentral F code.
.srsnr_by_hand <- function(F_value, df_1, df_2, N, alpha_lower, alpha_upper) {
  solve_ncp <- function(target) {
    if (target == 0) return(Inf)   # a one-sided bound
    if (stats::pf(F_value, df_1, df_2, ncp = 0) <= target) return(0)
    upper <- max(10, 2 * F_value * df_1)
    while (stats::pf(F_value, df_1, df_2, ncp = upper) > target) upper <- 2 * upper
    stats::uniroot(function(ncp) stats::pf(F_value, df_1, df_2, ncp = ncp) - target,
                   c(0, upper), tol = 1e-14)$root
  }
  sqrt(c(solve_ncp(1 - alpha_lower), solve_ncp(alpha_upper)) / N)
}

test_that("ci_srsnr() returns the expected 2-row tidy data frame", {
  result <- ci_srsnr(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, c("ci_lower", "ci_upper"))
})

test_that("ci_srsnr() takes var_within with means, and no earlier name", {
  expect_identical(names(formals(ci_srsnr)),
                   c("F_value", "df_1", "df_2", "N", "means", "var_within",
                     "n_per_group", "conf_level", "alpha_lower",
                     "alpha_upper", "..."))
  # ci_srsnr() has `...`, which passes on to the noncentral F search, so an
  # old name does not stop at the call itself. With 'means' the call stops
  # because 'var_within' is missing; with 'F_value' the old name reaches the
  # search, which stops with R's "unused argument".
  for (old in c("true_var_within", "sigma_squared")) {
    args <- list(means = c(1, 2, 3), n_per_group = 5, 1)
    names(args)[3] <- old
    expect_error(do.call(ci_srsnr, args),
                 "'var_within' and 'n_per_group' must be given",
                 class = "dmar_argument_invalid")
    args <- list(F_value = 5, df_1 = 2, df_2 = 27, N = 30, 1)
    names(args)[5] <- old
    expect_error(do.call(ci_srsnr, args), paste0("unused argument \\(", old))
  }
})

test_that("ci_srsnr() limits from F match the noncentral F interval computed by hand", {
  # The Bargmann (1970) example as Steiger (2004) reports it: F = 11.221 on
  # 4 and 50 degrees of freedom, N = 55.
  for (tails in list(c(.025, .025), c(.05, .05), c(.02, .03), c(.10, 0))) {
    result <- ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55,
                       alpha_lower = tails[1], alpha_upper = tails[2])
    expect_equal(result$value,
                 .srsnr_by_hand(11.221, 4, 50, 55, tails[1], tails[2]),
                 tolerance = 1e-10)
  }
  expect_equal(ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55,
                        conf_level = .90)$value,
               .srsnr_by_hand(11.221, 4, 50, 55, .05, .05), tolerance = 1e-10)
})

test_that("ci_srsnr() limits are the square roots of the corresponding ci_snr() limits", {
  snr_result    <- ci_snr(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  srsnr_result  <- ci_srsnr(F_value = 5, df_1 = 5, df_2 = 100, N = 110, conf_level = 0.95)
  expect_equal(srsnr_result$value, sqrt(pmax(0, snr_result$value)), tolerance = 1e-8)
})

test_that("ci_srsnr() from the summaries of observed data gives the interval of the observed F", {
  set.seed(2207)
  n <- c(9, 14, 11, 16)
  group <- factor(rep(seq_along(n), n))
  y <- stats::rnorm(sum(n), mean = c(0, .5, .3, 1)[group], sd = 1.1)
  tab <- stats::anova(stats::lm(y ~ group))
  by_summary <- ci_srsnr(means = tapply(y, group, mean),
                         var_within = tab[2, "Mean Sq"],
                         n_per_group = as.vector(table(group)))
  by_F <- ci_srsnr(F_value = tab[1, "F value"], df_1 = tab[1, "Df"],
                   df_2 = tab[2, "Df"], N = sum(n))
  expect_equal(by_summary$value, by_F$value, tolerance = 1e-12)
  expect_equal(by_summary$value,
               .srsnr_by_hand(tab[1, "F value"], 3, sum(n) - 4, sum(n),
                              .025, .025),
               tolerance = 1e-10)
  # The plug-in Cohen's f of the same summaries lies inside the interval.
  f_hat <- cohen_f(means = tapply(y, group, mean),
                   var_within = tab[2, "Mean Sq"], n_per_group = n)$value
  expect_gt(f_hat, by_summary$value[1])
  expect_lt(f_hat, by_summary$value[2])
  # Values of df_1, df_2, and N that agree with the summaries are accepted.
  expect_equal(ci_srsnr(means = tapply(y, group, mean),
                        var_within = tab[2, "Mean Sq"], n_per_group = n,
                        df_1 = 3, df_2 = sum(n) - 4, N = sum(n))$value,
               by_F$value, tolerance = 1e-12)
})

test_that("ci_srsnr() with population values gives the hypothetical interval at the implied F", {
  means <- c(94, 91, 92, 83)
  var_within <- 67.375
  n <- 20
  # The F a study of 20 per group would report if its summaries equaled
  # these values: SS_between = 20 * 70, on 3 and 76 degrees of freedom.
  F_implied <- (n * sum((means - mean(means))^2) / 3) / var_within
  result <- expect_silent(ci_srsnr(means = means, var_within = var_within,
                                   n_per_group = n))
  expect_equal(result$value, .srsnr_by_hand(F_implied, 3, 76, 80, .025, .025),
               tolerance = 1e-10)
})

test_that("ci_srsnr() summary route matches the implied F route when the lower limit is clamped", {
  means <- c(94, 91, 92, 83)
  sig2  <- 67.375
  n     <- 6
  J     <- length(means)
  N     <- J * n
  df1   <- J - 1
  df2   <- N - J
  mu_bar      <- mean(means)
  SS_between  <- n * sum((means - mu_bar)^2)
  F_implied   <- (SS_between / df1) / sig2

  # The small implied F legitimately clamps the lower limit to 0 on both paths.
  expect_warning(
    r_design <- ci_srsnr(means = means, var_within = sig2, n_per_group = n,
                         conf_level = 0.95),
    "below the alpha_lower critical value")
  expect_warning(
    r_F <- ci_srsnr(F_value = F_implied, df_1 = df1, df_2 = df2, N = N,
                    conf_level = 0.95),
    "below the alpha_lower critical value")
  expect_equal(r_design$value, r_F$value, tolerance = 1e-12)
  expect_equal(r_design$value, .srsnr_by_hand(F_implied, df1, df2, N, .025, .025),
               tolerance = 1e-10)
})

test_that("ci_srsnr() refuses an incomplete or conflicting specification", {
  invalid <- "dmar_argument_invalid"
  expect_error(
    ci_srsnr(F_value = 5, means = c(1, 2, 3),
             var_within = 1, n_per_group = 5),
    "either 'F_value'", class = invalid
  )
  expect_error(ci_srsnr(means = c(1, 2, 3), var_within = 1),
               "'n_per_group' must be", class = invalid)
  expect_error(ci_srsnr(means = c(1, 2, 3), n_per_group = 5),
               "'var_within' and 'n_per_group' must be", class = invalid)
  err <- expect_error(ci_srsnr(F_value = 5, df_1 = 2, df_2 = 27, N = 30,
                               var_within = 1), class = invalid)
  expect_match(conditionMessage(err), "used only with 'means'")
  expect_error(ci_srsnr(F_value = 5, df_1 = 2, df_2 = 27, N = 30,
                        n_per_group = 10),
               "used only with 'means'", class = invalid)
  # df_1, df_2, or N that disagree with the group sizes
  err <- expect_error(ci_srsnr(means = c(1, 2, 3), var_within = 1,
                               n_per_group = 10, N = 40), class = invalid)
  expect_identical(err$arg, "N")
  expect_match(conditionMessage(err), "N is 30")
  err <- expect_error(ci_srsnr(means = c(1, 2, 3), var_within = 1,
                               n_per_group = 10, df_1 = 3, df_2 = 26),
                      class = invalid)
  expect_identical(err$arg, c("df_1", "df_2"))
})

test_that("ci_srsnr() refuses invalid summary values with classed conditions", {
  invalid <- "dmar_argument_invalid"
  check <- function(args, arg, pattern) {
    err <- expect_error(do.call(ci_srsnr, args), class = invalid)
    expect_true(arg %in% err$arg)
    expect_match(conditionMessage(err), pattern)
  }
  check(list(means = c(1, 2, 3), var_within = 0, n_per_group = 5),
        "var_within", "positive")
  check(list(means = c(1, 2, 3), var_within = c(1, 2), n_per_group = 5),
        "var_within", "single")
  check(list(means = c(1, 2, 3), var_within = NA, n_per_group = 5),
        "var_within", "not missing")
  check(list(means = 2, var_within = 1, n_per_group = 5),
        "means", "at least two")
  check(list(means = c(1, NA, 3), var_within = 1, n_per_group = 5),
        "means", "finite")
  check(list(means = c(1, 2, 3), var_within = 1, n_per_group = c(5, 6)),
        "n_per_group", "same length")
  check(list(means = c(1, 2, 3), var_within = 1, n_per_group = c(5, 1, 5)),
        "n_per_group", "at least 2")
  # Group sizes are counts: they set df_2 and N.
  check(list(means = c(1, 2, 3), var_within = 1, n_per_group = c(7.5, 8, 9)),
        "n_per_group", "whole number of at least 2")
  check(list(means = c(1, 2, 3), var_within = 1, n_per_group = 10.5),
        "n_per_group", "group sizes were 10.5, 10.5, 10.5")
  # A whole number stored as a double, and an integer, are accepted.
  expect_equal(ci_srsnr(means = c(1, 2, 3), var_within = 1,
                        n_per_group = c(7, 8, 9))$value,
               ci_srsnr(means = c(1, 2, 3), var_within = 1,
                        n_per_group = c(7L, 8L, 9L))$value)
})

test_that("ci_srsnr() refuses equal group means as F = 0 and names 'means'", {
  invalid <- "dmar_argument_invalid"
  err <- expect_error(ci_srsnr(means = c(5, 5, 5), var_within = 2,
                               n_per_group = 10), class = invalid)
  expect_identical(err$arg, "means")
  expect_match(conditionMessage(err), "all equal, so the F computed from the summaries is 0")
  # The F route refuses the same F, naming 'F_value'.
  err <- expect_error(ci_srsnr(F_value = 0, df_1 = 2, df_2 = 27, N = 30),
                      class = invalid)
  expect_identical(err$arg, "F_value")
})

test_that("ci_srsnr() reports the lower-limit clamp once, in terms of the root signal-to-noise ratio, and still returns the interval", {
  msgs <- capture_warnings(
    result <- ci_srsnr(F_value = 1.2, df_1 = 4, df_2 = 50, N = 55)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "lower confidence limit on the square root of the signal-to-noise ratio is 0")
  expect_false(grepl("prob_greater", msgs))
  expect_equal(result$value[result$term == "ci_lower"], 0)
  expect_gt(result$value[result$term == "ci_upper"], 0)
})

test_that("ci_srsnr() names itself and the inputs when the noncentral root search fails", {
  expect_error(
    ci_srsnr(F_value = Inf, df_1 = 4, df_2 = 50, N = 55),
    "In ci_srsnr\\(\\).*F-statistic Inf"
  )
})
