# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# ci_ratio_of_means(): the ratio of two independent means, or of two
# groups' advantages over a third, baseline group.

# Two independent groups (rounded normal draws), used for the oracle pins.
x <- c(37.4, 56.9, 53.1, 44.7, 56.1, 52.1, 52.5, 62.1, 58.9, 40.2, 64.6,
       48.2, 58, 54.7, 49.3, 52.7, 72.4, 56.9, 51.9, 65.9, 55.2, 36.1, 48.2,
       43.6)
y <- c(31.4, 43.2, 47.7, 44.7, 56, 65.8, 66, 46.3, 41, 34.5, 43.2, 32.6,
       44.6, 46.1, 40.2, 54.5, 39, 27.5, 20, 50.5, -0.5, 40.1, 37, 41.9, 61.3,
       21, 35.5, 40.6, 52.8, 50.8)

value_of <- function(out, term) out$value[match(term, out$term)]
limits_of <- function(out) value_of(out, c("ci_lower", "ci_upper"))

ratio_terms <- c("ci_lower", "ratio", "ci_upper", "excluded_lower",
                 "excluded_upper", "mean_1", "mean_2", "mean_baseline",
                 "se_1", "se_2", "se_baseline", "t_denominator", "df", "n_1",
                 "n_2", "n_baseline")

test_that("pooled variances give Fieller's closed-form interval", {
  out <- ci_ratio_of_means(x, y, var_equal = TRUE)
  n1 <- length(x)
  n2 <- length(y)
  sp2 <- ((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2)
  q <- qt(0.975, n1 + n2 - 2)
  est <- mean(x) / mean(y)
  v11 <- sp2 / n1
  v22 <- sp2 / n2
  g <- q^2 * v22 / mean(y)^2
  half <- q / mean(y) * sqrt(v11 + est^2 * v22 - g * v11)
  expect_equal(limits_of(out), (est + c(-1, 1) * half) / (1 - g),
               tolerance = 1e-12)
  expect_equal(value_of(out, "ratio"), est)
  expect_equal(value_of(out, "df"), n1 + n2 - 2)
  expect_equal(value_of(out, "se_1"), sqrt(v11))
  expect_identical(attr(out, "confidence_set"), "interval")
  # Pinned from mratios::ttestratio(x, y, var.equal = TRUE) (mratios
  # 1.4.4, 2026-11-05); live comparison in tools/oracle_snippets/mratios.R.
  expect_equal(limits_of(out), c(1.10356403989741, 1.45658125339170),
               tolerance = 1e-12)
  expect_equal(limits_of(ci_ratio_of_means(x, y, var_equal = TRUE,
                                           conf_level = 0.90)),
               c(1.12894005237528, 1.42280131534973), tolerance = 1e-12)
})

test_that("unequal variances invert the Welch test, so the interval agrees with t.test()", {
  out <- ci_ratio_of_means(x, y)
  L <- limits_of(out)
  p <- vapply(L, function(l) t.test(x, l * y)$p.value, numeric(1))
  expect_equal(p, c(0.05, 0.05), tolerance = 1e-9)
  df_at <- vapply(L, function(l) unname(t.test(x, l * y)$parameter), 0)
  df_limits <- attr(out, "df_limits")
  expect_identical(names(df_limits), c("lower", "upper", "excluded_lower",
                                       "excluded_upper"))
  expect_equal(unname(df_limits[1:2]), df_at, tolerance = 1e-8)
  expect_true(all(is.na(df_limits[3:4])))
  # The df row is that of the denominator's t, the y group alone.
  expect_equal(value_of(out, "df"), length(y) - 1)
  expect_equal(value_of(out, "t_denominator"),
               mean(y) / (sd(y) / sqrt(length(y))))
  # Pinned from mratios::ttestratio(x, y, var.equal = FALSE,
  # iterativeCI = TRUE) (mratios 1.4.4); its uniroot() runs at the default
  # tolerance, which leaves its limits about 1e-5 from the exact inversion.
  expect_equal(L, c(1.10710046951399, 1.46395875091004), tolerance = 1e-4)
  one <- ci_ratio_of_means(x, y, alpha_lower = 0.05, alpha_upper = 0)
  expect_equal(value_of(one, "ci_lower"), 1.13142695018883,
               tolerance = 1e-4)
  expect_identical(value_of(one, "ci_upper"), Inf)
  expect_equal(t.test(x, value_of(one, "ci_lower") * y,
                      alternative = "greater")$p.value, 0.05,
               tolerance = 1e-9)
})

test_that("with a baseline group the limits invert the test of the three-arm contrast", {
  set.seed(2003)
  placebo <- rnorm(25, 50, 8)
  reference <- rnorm(30, 62, 10)
  candidate <- rnorm(28, 60, 13)
  contrast_p <- function(L, pooled) {
    groups <- list(candidate, reference, placebo)
    w <- c(1, -L, L - 1)
    n <- lengths(groups)
    m <- vapply(groups, mean, 0)
    s2 <- vapply(groups, var, 0)
    if (pooled) {
      v <- sum((n - 1) * s2) / sum(n - 1) / n
      df <- sum(n) - 3
    } else {
      v <- s2 / n
      df <- sum(w^2 * v)^2 / sum(w^4 * v^2 / (n - 1))
    }
    2 * pt(-abs(sum(w * m) / sqrt(sum(w^2 * v))), df)
  }
  for (pooled in c(TRUE, FALSE)) {
    out <- ci_ratio_of_means(candidate, reference, baseline = placebo,
                             var_equal = pooled)
    expect_equal(value_of(out, "ratio"), (mean(candidate) - mean(placebo)) /
                   (mean(reference) - mean(placebo)))
    expect_equal(vapply(limits_of(out), contrast_p, 0, pooled = pooled),
                 c(0.05, 0.05), tolerance = 1e-9)
    expect_equal(value_of(out, "n_baseline"), 25)
  }
  expect_equal(value_of(ci_ratio_of_means(candidate, reference,
                                          baseline = placebo,
                                          var_equal = TRUE), "df"),
               25 + 30 + 28 - 3)
})

test_that("noninferiority at a fraction holds exactly below the one-sided lower bound", {
  set.seed(7)
  placebo <- rnorm(80, 0, 1)
  reference <- rnorm(80, 1, 1)
  candidate <- rnorm(80, 0.85, 1)
  out <- ci_ratio_of_means(candidate, reference, baseline = placebo,
                           var_equal = TRUE, alpha_lower = 0.05,
                           alpha_upper = 0)
  lower <- value_of(out, "ci_lower")
  expect_identical(value_of(out, "ci_upper"), Inf)
  # The one-sided bound is the lower limit of the two-sided 90% interval.
  expect_equal(lower, value_of(ci_ratio_of_means(candidate, reference,
                                                 baseline = placebo,
                                                 var_equal = TRUE,
                                                 conf_level = 0.90),
                               "ci_lower"))
  # The one-sided test of (mu_c - mu_p) - f (mu_r - mu_p) <= 0 rejects at .05
  # for every f below the bound and for none above it.
  sp2 <- (var(candidate) + var(reference) + var(placebo)) / 3
  rejects <- function(f) {
    est <- mean(candidate) - f * mean(reference) - (1 - f) * mean(placebo)
    se <- sqrt(sp2 * (1 + f^2 + (1 - f)^2) / 80)
    est / se > qt(0.95, 237)
  }
  f <- seq(0.3, 1.2, by = 0.01)
  f <- f[abs(f - lower) > 1e-6]
  expect_identical(vapply(f, rejects, logical(1)), f < lower)
})

test_that("summary statistics give the same interval as the raw data", {
  for (pooled in c(TRUE, FALSE)) {
    raw <- ci_ratio_of_means(x, y, var_equal = pooled)
    smry <- ci_ratio_of_means(mean_1 = mean(x), sd_1 = sd(x),
                              n_1 = length(x), mean_2 = mean(y),
                              sd_2 = sd(y), n_2 = length(y),
                              var_equal = pooled)
    expect_equal(smry$value, raw$value, tolerance = 1e-10)
    expect_identical(smry$term, ratio_terms)
  }
})

test_that("a baseline group from summary statistics gives the raw data's interval", {
  set.seed(2003)
  placebo <- rnorm(25, 50, 8)
  reference <- rnorm(30, 62, 10)
  candidate <- rnorm(28, 60, 13)
  for (pooled in c(TRUE, FALSE)) {
    for (method in c("fieller", "log")) {
      raw <- ci_ratio_of_means(candidate, reference, baseline = placebo,
                               var_equal = pooled, ci_method = method)
      smry <- ci_ratio_of_means(mean_1 = mean(candidate),
                                sd_1 = sd(candidate), n_1 = 28,
                                mean_2 = mean(reference),
                                sd_2 = sd(reference), n_2 = 30,
                                mean_baseline = mean(placebo),
                                sd_baseline = sd(placebo), n_baseline = 25,
                                var_equal = pooled, ci_method = method)
      expect_equal(smry$value, raw$value, tolerance = 1e-10)
      expect_identical(smry$term, ratio_terms)
      expect_equal(attr(smry, "df_limits"), attr(raw, "df_limits"),
                   tolerance = 1e-10)
    }
  }
  expect_error(ci_ratio_of_means(mean_1 = 1, sd_1 = 1, n_1 = 5, mean_2 = 2,
                                 sd_2 = 1, n_2 = 6, mean_baseline = 0.5),
               "all three")
  expect_error(ci_ratio_of_means(x, y, mean_baseline = 0.5, sd_baseline = 1,
                                 n_baseline = 5), "not both")
  e <- tryCatch(ci_ratio_of_means(mean_1 = 1, sd_1 = 1, n_1 = 5, mean_2 = 2,
                                  sd_2 = 1, n_2 = 6, mean_baseline = 0.5,
                                  sd_baseline = 1, n_baseline = 5,
                                  ci_method = "percentile"),
                dmar_refusal = function(e) e)
  expect_s3_class(e, c("dmar_bootstrap_needs_data", "dmar_refusal"))
})

test_that("missing values are removed from each group separately", {
  expect_equal(ci_ratio_of_means(c(x, NA), c(NA, y))$value,
               ci_ratio_of_means(x, y)$value)
})

test_that("an unbounded set is reported as it is, with its gap", {
  set.seed(3)
  x2 <- rnorm(15, 5, 2)
  y2 <- rnorm(12, 0.8, 2)
  out <- ci_ratio_of_means(x2, y2)
  expect_identical(attr(out, "confidence_set"), "complement of an interval")
  expect_identical(limits_of(out), c(-Inf, Inf))
  gap <- value_of(out, c("excluded_lower", "excluded_upper"))
  expect_true(all(is.finite(gap)) && gap[1L] < gap[2L])
  expect_equal(vapply(gap, function(l) t.test(x2, l * y2)$p.value, 0),
               c(0.05, 0.05), tolerance = 1e-9)
  # The degrees of freedom of the test at the ends of the gap, the set's
  # only finite ends.
  df_limits <- attr(out, "df_limits")
  expect_true(all(is.na(df_limits[c("lower", "upper")])))
  expect_equal(unname(df_limits[c("excluded_lower", "excluded_upper")]),
               vapply(gap, function(l) unname(t.test(x2, l * y2)$parameter),
                      0), tolerance = 1e-8)
  expect_lt(abs(value_of(out, "t_denominator")),
            qt(0.975, value_of(out, "df")))
  expect_output(print(out), "complement of an interval")
  td <- generics::tidy(out)
  expect_equal(c(td$excluded_lower, td$excluded_upper), gap)
  expect_identical(td$confidence_set, "complement of an interval")

  x3 <- c(-1.2, 0.4, 1.1, 2.0, -0.6, 0.9)
  y3 <- c(0.3, -1.1, 1.2, -0.4, 0.6, 0.2)
  for (pooled in c(TRUE, FALSE)) {
    whole <- ci_ratio_of_means(x3, y3, var_equal = pooled)
    expect_identical(attr(whole, "confidence_set"), "whole line")
    expect_true(all(is.na(value_of(whole, c("excluded_lower",
                                            "excluded_upper")))))
  }
})

test_that("a narrow excluded gap under unequal variances is found", {
  # The Welch test's one-sided p-value is smallest, 0.000694789, near
  # rho = -20.52, so at a tail area just above it the set excludes a narrow
  # gap there. Its boundary peaks near, not at, the peak of t, because the
  # degrees of freedom change with the ratio, and the gap is narrower in
  # the angle atan(rho) than the search grid's spacing.
  set.seed(8)
  x8 <- rnorm(5, 3, 6)
  y8 <- rnorm(40, 0.5, 1)
  a <- 0.000694789 * (1 + 1e-4)
  out <- ci_ratio_of_means(x8, y8, alpha_lower = a, alpha_upper = a)
  expect_identical(attr(out, "confidence_set"), "complement of an interval")
  gap <- value_of(out, c("excluded_lower", "excluded_upper"))
  expect_true(gap[1L] > -20.7 && gap[2L] < -20.3)
  expect_equal(vapply(gap, function(l) t.test(x8, l * y8)$p.value, 0),
               c(2 * a, 2 * a), tolerance = 1e-8)
  expect_lt(t.test(x8, mean(gap) * y8)$p.value, 2 * a)
  # Below the smallest p-value no ratio is rejected.
  a_below <- 0.000694789 * (1 - 1e-3)
  expect_identical(attr(ci_ratio_of_means(x8, y8, alpha_lower = a_below,
                                          alpha_upper = a_below),
                        "confidence_set"), "whole line")
})

test_that("the log interval is the response ratio interval and needs positive terms", {
  out <- ci_ratio_of_means(x, y, ci_method = "log", var_equal = TRUE)
  sp2 <- ((23) * var(x) + (29) * var(y)) / 52
  se <- sqrt(sp2 / 24 / mean(x)^2 + sp2 / 30 / mean(y)^2)
  expect_equal(limits_of(out),
               exp(log(mean(x) / mean(y)) + c(-1, 1) * qt(0.975, 52) * se),
               tolerance = 1e-12)
  # Pinned from metafor::escalc("ROM", correct = FALSE) (metafor 5.0.1,
  # 2026-11-05): yi is the log ratio, and vi is the squared standard error,
  # vtype "LS" for unequal and "HO" for pooled variances. Live comparison
  # in tools/oracle_snippets/metafor.R.
  expect_equal(log(value_of(out, "ratio")), 0.23612355151674830,
               tolerance = 1e-12)
  expect_equal(se^2, 0.00474425205250613, tolerance = 1e-10)
  welch <- ci_ratio_of_means(x, y, ci_method = "log")
  se_w <- sqrt(var(x) / 24 / mean(x)^2 + var(y) / 30 / mean(y)^2)
  expect_equal(se_w^2, 0.00476047839886595, tolerance = 1e-10)
  df_w <- se_w^4 / ((var(x) / 24 / mean(x)^2)^2 / 23 +
                      (var(y) / 30 / mean(y)^2)^2 / 29)
  expect_equal(unname(attr(welch, "df_limits")), c(df_w, df_w, NA, NA))
  expect_equal(limits_of(welch),
               exp(log(mean(x) / mean(y)) + c(-1, 1) * qt(0.975, df_w) * se_w),
               tolerance = 1e-12)

  e <- tryCatch(ci_ratio_of_means(x, -y, ci_method = "log"),
                dmar_refusal = function(e) e)
  expect_s3_class(e, c("dmar_log_ratio_undefined", "dmar_refusal", "error"))
  expect_identical(e$code, "log_ratio_undefined")
})

test_that("the log interval with a baseline uses the linearized log ratio", {
  set.seed(2003)
  placebo <- rnorm(25, 50, 8)
  reference <- rnorm(30, 62, 10)
  candidate <- rnorm(28, 60, 13)
  out <- ci_ratio_of_means(candidate, reference, baseline = placebo,
                           ci_method = "log")
  num <- mean(candidate) - mean(placebo)
  den <- mean(reference) - mean(placebo)
  v <- c(var(candidate) / 28, var(reference) / 30, var(placebo) / 25)
  # log(num) - log(den) moves by dm_1 / num - dm_2 / den + dm_b (1 / den -
  # 1 / num) when the three means move by dm.
  w <- c(1 / num, -1 / den, 1 / den - 1 / num)
  se <- sqrt(sum(w^2 * v))
  expect_equal(se^2, (v[1] + v[3]) / num^2 + (v[2] + v[3]) / den^2 -
                 2 * v[3] / (num * den))
  df <- sum(w^2 * v)^2 / sum(w^4 * v^2 / c(27, 29, 24))
  expect_equal(unname(attr(out, "df_limits")), c(df, df, NA, NA))
  expect_equal(limits_of(out),
               exp(log(num / den) + c(-1, 1) * qt(0.975, df) * se),
               tolerance = 1e-12)
  pooled <- ci_ratio_of_means(candidate, reference, baseline = placebo,
                              ci_method = "log", var_equal = TRUE)
  sp2 <- sum(c(27, 29, 24) * c(var(candidate), var(reference),
                               var(placebo))) / 80
  vp <- sp2 / c(28, 30, 25)
  sep <- sqrt(sum(w^2 * vp))
  expect_equal(limits_of(pooled),
               exp(log(num / den) + c(-1, 1) * qt(0.975, 80) * sep),
               tolerance = 1e-12)
  expect_null(attr(pooled, "df_limits"))
})

test_that("the bootstrap is reproducible from a seed and leaves the caller's stream alone", {
  set.seed(1)
  before <- .Random.seed
  a <- ci_ratio_of_means(x, y, ci_method = "percentile", B = 200, seed = 11)
  expect_identical(.Random.seed, before)
  b <- ci_ratio_of_means(x, y, ci_method = "percentile", B = 200, seed = 11)
  expect_identical(a$value, b$value)
  expect_identical(attr(a, "B"), 200L)
  expect_identical(attr(a, "B_used"), 200L)
  expect_identical(a$term, ratio_terms)
  # Each group is resampled from itself, the numerator group first, so the
  # same seed redraws the same replicates.
  set.seed(11)
  bx <- DMAR:::.ratio_boot_means(x, 200L)
  by <- DMAR:::.ratio_boot_means(y, 200L)
  expect_equal(limits_of(a), unname(quantile(bx / by, c(0.025, 0.975))))
  bca <- ci_ratio_of_means(x, y, ci_method = "bca", B = 200, seed = 11)
  expect_true(all(is.finite(limits_of(bca))))
  expect_identical(attr(bca, "ci_method"), "bca")
  md <- dmar_metadata(bca)
  expect_identical(md$bootstrap$B, 200L)
  expect_identical(md$ci_method, "bca")
})

test_that("replications with a zero denominator are dropped with one warning", {
  # Discrete scores: the resampled reference and baseline means often tie.
  xx <- c(2, 3, 3, 4)
  yy <- c(1, 1, 2, 2)
  zz <- c(1, 1, 1, 2)
  expect_warning(
    out <- ci_ratio_of_means(xx, yy, baseline = zz,
                             ci_method = "percentile", B = 200, seed = 4),
    "of 200 bootstrap replications \\([0-9.]+%\\) had a denominator of zero",
    class = "dmar_bootstrap_replications_dropped")
  expect_lt(attr(out, "B_used"), 200L)
  expect_identical(attr(out, "n_dropped"), 200L - attr(out, "B_used"))
  expect_gte(attr(out, "B_used"), 100L)
  expect_output(print(out), "replications, [0-9]+ used")
  # The Fieller set for these data is unbounded, and the line under the
  # bootstrap interval says so.
  expect_identical(attr(out, "fieller_set"), "the complement of an interval")
  expect_output(print(out), "The Fieller set at this level is")
})

test_that("a larger bootstrap is close to Fieller's interval for well-determined means", {
  skip_on_cran()
  out <- ci_ratio_of_means(x, y, ci_method = "bca", B = 10000, seed = 2)
  expect_equal(limits_of(out), limits_of(ci_ratio_of_means(x, y)),
               tolerance = 0.05)
})

test_that("the table has one shape and tidies to one row", {
  out <- ci_ratio_of_means(x, y)
  expect_identical(class(out), c("dmar_ci_ratio_of_means", "dmar_ci_long",
                                 "dmar_tbl", "data.frame"))
  expect_identical(out$term, ratio_terms)
  expect_true(is.numeric(out$value))
  expect_true(is.na(value_of(out, "mean_baseline")))
  expect_identical(attr(out, "ci_method"), "fieller")
  expect_false(attr(out, "var_equal"))
  expect_output(print(out), "Fieller confidence set: a bounded interval")
  td <- generics::tidy(out)
  expect_identical(names(td), c("term", "estimate", "ci_lower", "ci_upper",
                                "excluded_lower", "excluded_upper",
                                "confidence_set", "conf_level"))
  expect_identical(td$term, "ratio")
  expect_equal(c(td$ci_lower, td$ci_upper), limits_of(out))
  expect_identical(generics::glance(out), td)
  one <- generics::tidy(ci_ratio_of_means(x, y, alpha_lower = 0.05,
                                          alpha_upper = 0))
  expect_identical(c(one$alpha_lower, one$alpha_upper), c(0.05, 0))
})

test_that("inputs are checked", {
  expect_error(ci_ratio_of_means(x, y, mean_1 = 1), "not both")
  expect_error(ci_ratio_of_means(x), "Supply both 'x'")
  expect_error(ci_ratio_of_means(mean_1 = 1, sd_1 = 1, n_1 = 5),
               "all six")
  expect_error(ci_ratio_of_means(x, rep(3, 10)), "no variability")
  expect_error(ci_ratio_of_means(x, c(1, NA)), "at least two")
  expect_error(ci_ratio_of_means(x, y, var_equal = NA), "TRUE or FALSE")
  expect_error(ci_ratio_of_means(x, y, ci_method = "percentile", B = 1),
               "single whole number", class = "dmar_count_invalid")
  expect_error(ci_ratio_of_means(mean_1 = 1, sd_1 = 1, n_1 = 5, mean_2 = 2,
                                 sd_2 = 1, n_2 = 1.5), "whole number")
  expect_error(ci_ratio_of_means(x, y, conf_level = 0.9, alpha_lower = 0.05,
                                 alpha_upper = 0.01), "inconsistent")
  expect_error(ci_ratio_of_means(x, c(-1, 1, -2, 2)), "exactly zero")
  e <- tryCatch(ci_ratio_of_means(mean_1 = 1, sd_1 = 1, n_1 = 5, mean_2 = 2,
                                  sd_2 = 1, n_2 = 6, ci_method = "bca"),
                dmar_refusal = function(e) e)
  expect_s3_class(e, c("dmar_bootstrap_needs_data", "dmar_refusal"))
})
