# Oracles in this file are computed separately from R/cohen_f.R: the
# variance of the means through stats::var() rescaled, through matrix
# products, or through the F of stats::anova() on raw data, and the
# published values of Kim (2016) and Cohen (1988) at the precision they
# print. Every value is computed at full precision; a value recorded to 12
# decimals in the evidence is compared after rounding to 12 decimals.

test_that("cohen_f() returns a one-row dmar_tbl with the term cohen_f", {
  result <- cohen_f(eta_squared = 0.10)
  expect_s3_class(result, "dmar_tbl")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term, "cohen_f")
  expect_type(result$value, "double")
})

test_that("cohen_f() takes the neutral quantity names and no earlier ones", {
  expect_identical(names(formals(cohen_f)),
                   c("means", "var_within", "n_per_group", "eta_squared",
                     "sd_means", "sd_within"))
  unused <- "unused argument"
  # DMAR 1.0.x names
  expect_error(cohen_f(mu = c(1, 2), sigma_squared = 1), unused)
  expect_error(cohen_f(sigma_m = 1, sigma = 2), unused)
  # Names of the 1.1.0 development builds
  expect_error(cohen_f(mu = c(1, 2), true_var_within = 1), unused)
  expect_error(cohen_f(true_sd_means = 1, true_sd = 2), unused)
})

test_that("cohen_f() from means with equal weights uses the variance of the means with divisor k", {
  means <- c(94, 91, 92, 83)
  var_within <- 67.375
  k <- length(means)
  # stats::var() divides by k - 1; rescale to the divisor k.
  oracle <- sqrt(stats::var(means) * (k - 1) / k / var_within)
  expect_equal(cohen_f(means = means, var_within = var_within)$value,
               oracle, tolerance = 1e-14)
  # A single n_per_group states equal groups and changes nothing.
  expect_equal(cohen_f(means = means, var_within = var_within,
                       n_per_group = 6)$value, oracle, tolerance = 1e-14)
})

test_that("cohen_f() with unequal n_per_group weights the means by n_j / N", {
  means <- c(94, 91, 92, 83)
  var_within <- 67.375
  n <- c(4, 6, 5, 5)
  # Weighted variance of the means by matrix algebra: W = diag(n / N), with
  # the means centered at the weighted grand mean.
  w <- n / sum(n)
  centered <- means - drop(crossprod(w, means))
  oracle <- sqrt(drop(t(centered) %*% diag(w) %*% centered) / var_within)
  expect_equal(cohen_f(means = means, var_within = var_within,
                       n_per_group = n)$value, oracle, tolerance = 1e-14)
  # Equal weighting gives a different value for these unequal groups.
  expect_false(isTRUE(all.equal(
    cohen_f(means = means, var_within = var_within)$value, oracle)))
  # Only the proportions n_j / N matter.
  expect_equal(cohen_f(means = means, var_within = var_within,
                       n_per_group = 10 * n)$value, oracle, tolerance = 1e-14)
})

test_that("cohen_f() from sample summaries is the plug-in sqrt(df_1 F / N) of the analysis of variance", {
  set.seed(4127)
  n <- c(7, 12, 9, 15)
  group <- factor(rep(seq_along(n), n))
  y <- stats::rnorm(sum(n), mean = c(0, .4, .2, .9)[group], sd = 1.3)
  tab <- stats::anova(stats::lm(y ~ group))
  F_value <- tab[1, "F value"]
  df_1 <- tab[1, "Df"]
  N <- sum(n)
  k <- length(n)
  oracle <- sqrt(df_1 * F_value / N)
  result <- cohen_f(means = tapply(y, group, mean),
                    var_within = tab[2, "Mean Sq"],
                    n_per_group = as.vector(table(group)))
  expect_equal(result$value, oracle, tolerance = 1e-12)

  # A sample eta squared, SS_between / SS_total, pairs with the within-group
  # variance SS_within / N rather than the mean square SS_within / (N - k),
  # so its plug-in value is larger by sqrt(N / (N - k)).
  ss_b <- tab[1, "Sum Sq"]
  ss_w <- tab[2, "Sum Sq"]
  f_eta <- cohen_f(eta_squared = ss_b / (ss_b + ss_w))$value
  expect_equal(f_eta, sqrt(ss_b / ss_w), tolerance = 1e-12)
  expect_equal(f_eta / result$value, sqrt(N / (N - k)), tolerance = 1e-12)
  expect_equal(cohen_f(means = tapply(y, group, mean), var_within = ss_w / N,
                       n_per_group = n)$value, f_eta, tolerance = 1e-12)
})

test_that("cohen_f() on drinks_trial: sample eta squared exceeds the means route by sqrt(N / (N - k))", {
  # The help example. Sums of squares by hand from the raw scores, without
  # anova() or lm().
  y <- drinks_trial$log_drinks
  g <- drinks_trial$treatment
  groups <- split(y, g)
  n <- vapply(groups, length, numeric(1))
  ybar_j <- vapply(groups, mean, numeric(1))
  N <- sum(n)
  k <- length(n)
  ybar <- sum(y) / N
  ss_b <- sum(n * (ybar_j - ybar)^2)
  ss_w <- sum(vapply(groups, function(x) sum((x - mean(x))^2), numeric(1)))
  expect_equal(c(N, k), c(88, 3))

  f_means <- cohen_f(means = ybar_j, var_within = ss_w / (N - k),
                     n_per_group = n)$value
  f_eta <- cohen_f(eta_squared = ss_b / (ss_b + ss_w))$value
  expect_equal(f_means, sqrt(ss_b / N / (ss_w / (N - k))), tolerance = 1e-13)
  expect_equal(f_eta, sqrt(ss_b / ss_w), tolerance = 1e-13)
  expect_gt(f_eta, f_means)
  expect_equal(f_eta / f_means, sqrt(88 / 85), tolerance = 1e-13)
  # The two standard deviations give either value, by which within-group
  # variance they use.
  expect_equal(cohen_f(sd_means = sqrt(ss_b / N),
                       sd_within = sqrt(ss_w / (N - k)))$value,
               f_means, tolerance = 1e-13)
  expect_equal(cohen_f(sd_means = sqrt(ss_b / N),
                       sd_within = sqrt(ss_w / N))$value,
               f_eta, tolerance = 1e-13)
  # Printed values in the help example, at the printed precision
  expect_equal(round(f_means, 3), .311)
  expect_equal(round(f_eta, 3), .317)
})

test_that("cohen_f() gives one value by its three specifications on consistent inputs", {
  means <- c(12, 15, 11, 20, 17)
  var_within <- 30
  n <- c(8, 12, 10, 6, 14)
  for (weights in list(NULL, n)) {
    w <- if (is.null(weights)) rep(1 / length(means), length(means)) else
      weights / sum(weights)
    var_means <- sum(w * (means - sum(w * means))^2)
    eta2 <- var_means / (var_means + var_within)
    f_means <- cohen_f(means = means, var_within = var_within,
                       n_per_group = weights)$value
    f_eta <- cohen_f(eta_squared = eta2)$value
    f_sd <- cohen_f(sd_means = sqrt(var_means),
                    sd_within = sqrt(var_within))$value
    expect_equal(f_eta, f_means, tolerance = 1e-13)
    expect_equal(f_sd, f_means, tolerance = 1e-13)
    expect_equal(f_means, sqrt(var_means / var_within), tolerance = 1e-14)
  }
  # The boundary of each specification
  expect_equal(cohen_f(eta_squared = 0)$value, 0)
  expect_equal(cohen_f(sd_means = 0, sd_within = 2)$value, 0)
  expect_equal(cohen_f(means = c(5, 5, 5), var_within = 2)$value, 0)
})

test_that("cohen_f() reproduces Kim (2016, Table 1) and its power plan", {
  # Kim (2016), Restorative Dentistry & Endodontics 41(3), p. 231, Table 1:
  # means 18.4, 22.2, 25.1, 32.1 with a common SD of 8.6; the article prints
  # f as 0.58 (rounded comparison) and, on p. 232, a total N of 40 with an
  # actual power of 0.847 from G*Power for power .80 at alpha .05 (rounded
  # comparison). The full-precision values are DMAR's, checked against
  # direct algebra and base R's noncentral F.
  means <- c(18.4, 22.2, 25.1, 32.1)
  f <- cohen_f(means = means, var_within = 8.6^2)$value
  f_algebra <- sqrt(sum((means - 24.45)^2) / 4) / 8.6   # grand mean 24.45
  expect_equal(f, f_algebra, tolerance = 1e-14)
  expect_equal(round(f, 12), 0.583165908653)
  expect_equal(round(f, 2), 0.58)

  plan <- ss_power_one_way_anova(a = 4, true_f = f, desired_power = .80,
                                 alpha_level = .05)
  power_at <- function(N, f) {
    1 - stats::pf(stats::qf(.95, 3, N - 4), 3, N - 4, ncp = N * f^2)
  }
  expect_equal(plan$value[plan$term == "necessary_N"], 40)
  expect_equal(plan$value[plan$term == "n_per_group"], 10)
  expect_equal(plan$value[plan$term == "actual_power"], power_at(40, f),
               tolerance = 1e-12)
  expect_equal(round(power_at(40, f), 12), 0.846617726244)
  expect_equal(round(power_at(40, f), 3), 0.847)
  # 40 is the smallest balanced total that reaches .80.
  expect_lt(power_at(36, f), .80)
  expect_equal(round(power_at(36, f), 12), 0.797040832003)
  # Rounding f to the printed .58 before planning gives a lower power at
  # N = 40, an illustration of why the full precision is kept.
  expect_equal(round(power_at(40, .58), 12), 0.842326131984)
})

test_that("cohen_f() reproduces Cohen (1988, p. 277) for seven means in Pattern 1", {
  # Cohen (1988), Section 8.2.1, p. 277, Equation 8.2.8: f = d sqrt(1 / (2k))
  # for k means, one at each end of the range d and the rest at the
  # midpoint. Seven means spanning d = .5 within-population standard
  # deviations give f = .267(.5) = .13 as printed (rounded comparison);
  # Table 8.2.1 (p. 278) prints the constant .267 for k = 7.
  f <- cohen_f(means = c(-.25, 0, 0, 0, 0, 0, .25), var_within = 1)$value
  expect_equal(f, .5 * sqrt(1 / (2 * 7)), tolerance = 1e-14)
  expect_equal(round(f, 12), 0.133630620956)
  expect_equal(round(f, 2), 0.13)
  expect_equal(round(sqrt(1 / 14), 3), 0.267)
  # The same means on another scale give the same f.
  expect_equal(cohen_f(means = 100 + 15 * c(-.25, 0, 0, 0, 0, 0, .25),
                       var_within = 15^2)$value, f, tolerance = 1e-13)
})

test_that("cohen_f() refuses an empty, incomplete, or conflicting specification", {
  invalid <- "dmar_argument_invalid"
  expect_error(cohen_f(), class = invalid)
  expect_error(cohen_f(), "needs one complete specification")
  # Conflicting: arguments of two specifications
  conflicting <- list(
    list(eta_squared = .1, means = c(1, 2)),
    list(eta_squared = .1, means = c(1, 2), var_within = 1),
    list(eta_squared = .1, n_per_group = 5),
    list(eta_squared = .1, sd_within = 2),
    list(eta_squared = .1, sd_means = 1, sd_within = 2),
    list(means = c(1, 2), var_within = 1, sd_means = 1, sd_within = 2),
    list(var_within = 1, sd_within = 2))
  for (args in conflicting) {
    err <- expect_error(do.call(cohen_f, args), class = invalid)
    expect_match(conditionMessage(err), "different specifications of f")
    expect_setequal(err$arg, names(args))
  }
  # Incomplete: one half of a pair, or n_per_group without the means
  err <- expect_error(cohen_f(means = c(1, 2)), class = invalid)
  expect_match(conditionMessage(err), "'means' was given without 'var_within'")
  err <- expect_error(cohen_f(var_within = 1), class = invalid)
  expect_match(conditionMessage(err), "'var_within' was given without 'means'")
  err <- expect_error(cohen_f(means = c(1, 2), n_per_group = 5), class = invalid)
  expect_match(conditionMessage(err), "without 'var_within'")
  err <- expect_error(cohen_f(n_per_group = c(5, 6)), class = invalid)
  expect_match(conditionMessage(err), "'n_per_group' was given without")
  err <- expect_error(cohen_f(sd_means = 1), class = invalid)
  expect_match(conditionMessage(err), "'sd_means' was given without 'sd_within'")
  err <- expect_error(cohen_f(sd_within = 1), class = invalid)
  expect_match(conditionMessage(err), "'sd_within' was given without 'sd_means'")
})

test_that("cohen_f() refuses invalid values with classed conditions that name the argument", {
  invalid <- "dmar_argument_invalid"
  check <- function(args, arg, pattern) {
    err <- expect_error(do.call(cohen_f, args), class = invalid)
    expect_true(arg %in% err$arg)
    expect_match(conditionMessage(err), pattern)
  }
  check(list(eta_squared = -0.1), "eta_squared", "in \\[0, 1\\)")
  check(list(eta_squared = 1), "eta_squared", "in \\[0, 1\\)")
  check(list(eta_squared = c(.1, .2)), "eta_squared", "single number")
  check(list(eta_squared = NA), "eta_squared", "not missing")
  check(list(eta_squared = "0.1"), "eta_squared", "not missing")
  check(list(means = 3, var_within = 1), "means", "at least two")
  check(list(means = c(1, NA), var_within = 1), "means", "finite")
  check(list(means = c(1, 2), var_within = 0), "var_within", "positive")
  check(list(means = c(1, 2), var_within = -1), "var_within", "positive")
  check(list(means = c(1, 2), var_within = c(1, 2)), "var_within", "single")
  check(list(means = c(1, 2), var_within = NA), "var_within", "not missing")
  check(list(means = c(1, 2, 3), var_within = 1, n_per_group = c(5, 6)),
        "n_per_group", "same length")
  check(list(means = c(1, 2), var_within = 1, n_per_group = c(5, 0)),
        "n_per_group", "positive")
  check(list(means = c(1, 2), var_within = 1, n_per_group = c(5, Inf)),
        "n_per_group", "positive")
  check(list(sd_means = 1, sd_within = 0), "sd_within", "positive")
  check(list(sd_means = -1, sd_within = 1), "sd_means", "at least 0")
  check(list(sd_means = c(1, 2), sd_within = 1), "sd_means", "single")
})
