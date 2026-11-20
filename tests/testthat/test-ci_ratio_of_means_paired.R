# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# ci_ratio_of_means_paired(): the ratio of two means measured on the same
# units, or of their advantages over a third, aligned variable.

value_of <- function(out, term) out$value[match(term, out$term)]
limits_of <- function(out) value_of(out, c("ci_lower", "ci_upper"))

paired_terms <- c("ci_lower", "ratio", "ci_upper", "excluded_lower",
                  "excluded_upper", "mean_1", "mean_2", "mean_baseline",
                  "se_numerator", "se_denominator", "r", "t_denominator",
                  "df", "N", "n_clusters")

# The Cushny-Peebles hours of extra sleep under two hypnotics, the same ten
# patients under each (the values of datasets::sleep).
drug_a <- c(0.7, -1.6, -0.2, -1.2, -0.1, 3.4, 3.7, 0.8, 0.0, 2.0)
drug_b <- c(1.9, 0.8, 1.1, 0.1, -0.1, 4.4, 5.5, 1.6, 4.6, 3.4)

test_that("Fieller's (1954) three forms are reproduced from his example", {
  # Fieller (1954, p. 180): at P = .95 the limits are -0.485 and 0.657,
  # inclusive; at .99, -2.920 and 0.808, inclusive; at .995, 0.896 and
  # 146.3, exclusive; at .999 none exist. He worked from summary values
  # rounded to four places and the rounded t^2 = 13.61 next to the
  # asymptote t^2 = 13.54, where the outer root moves fast, so 146.3 is
  # matched to within one percent.
  ci <- function(level) ci_ratio_of_means_paired(drug_a, drug_b,
                                                 conf_level = level)
  p95 <- ci(0.95)
  expect_identical(attr(p95, "confidence_set"), "interval")
  expect_equal(limits_of(p95), c(-0.485, 0.657), tolerance = 1e-3)
  p99 <- ci(0.99)
  expect_equal(limits_of(p99), c(-2.920, 0.808), tolerance = 1e-3)
  p995 <- ci(0.995)
  expect_identical(attr(p995, "confidence_set"), "complement of an interval")
  expect_identical(limits_of(p995), c(-Inf, Inf))
  expect_equal(value_of(p995, "excluded_lower"), 0.896, tolerance = 1e-3)
  expect_equal(value_of(p995, "excluded_upper"), 146.3, tolerance = 0.01)
  p999 <- ci(0.999)
  expect_identical(attr(p999, "confidence_set"), "whole line")
  expect_true(all(is.na(value_of(p999, c("excluded_lower",
                                         "excluded_upper")))))
  expect_equal(value_of(p95, "ratio"), 0.75 / 2.33)
  expect_equal(value_of(p95, "t_denominator"), 2.33 / sqrt(var(drug_b) / 10))
  # The same from the summary values Fieller printed, as the help shows.
  smry <- ci_ratio_of_means_paired(mean_1 = 0.75, sd_1 = sqrt(10 * 0.3201),
                                   mean_2 = 2.33, sd_2 = sqrt(10 * 0.4009),
                                   r = 0.2848 / sqrt(0.3201 * 0.4009),
                                   N = 10)
  expect_equal(limits_of(smry), c(-0.485, 0.657), tolerance = 1e-3)
})

test_that("each finite limit is where the one-sample t test of x - L y has p = alpha", {
  for (level in c(0.90, 0.95, 0.995)) {
    out <- ci_ratio_of_means_paired(drug_a, drug_b, conf_level = level)
    ends <- value_of(out, c("ci_lower", "ci_upper", "excluded_lower",
                            "excluded_upper"))
    ends <- ends[is.finite(ends)]
    p <- vapply(ends, function(L) t.test(drug_a - L * drug_b)$p.value, 0)
    expect_equal(p, rep(1 - level, length(ends)), tolerance = 1e-9)
  }
})

test_that("the retained fraction of the TOIS application is reproduced", {
  # Summary statistics of the per-user advantages over the MostPop
  # recommender on MovieLens 100K (nDCG@10, 943 users): PureSVD minus
  # MostPop and EASE minus MostPop. The per-user data do not ship (their
  # redistribution terms are unverified); these six numbers determine the
  # paired Fieller interval. The manuscript reports rho-hat .877 with 90
  # percent limits .849 and .905; the 95 percent limits are .844 and .911.
  ml <- list(mean_1 = 0.161900190775728, sd_1 = 0.173855817899603,
             mean_2 = 0.184617167557725, sd_2 = 0.179862516950324,
             r = 0.834233262646623, N = 943)
  p95 <- do.call(ci_ratio_of_means_paired, ml)
  expect_equal(value_of(p95, c("ci_lower", "ratio", "ci_upper")),
               c(0.8437, 0.8770, 0.9108), tolerance = 1e-4)
  expect_equal(value_of(p95, c("ci_lower", "ratio", "ci_upper")),
               c(0.8436845, 0.8769509, 0.9107666), tolerance = 1e-6)
  p90 <- do.call(ci_ratio_of_means_paired, c(ml, conf_level = 0.90))
  expect_equal(limits_of(p90), c(0.8490, 0.9053), tolerance = 1e-4)
  expect_equal(limits_of(p90), c(0.8490205, 0.9052675), tolerance = 1e-6)
  expect_equal(value_of(p95, "t_denominator"), 31.52008, tolerance = 1e-6)
  # The one-sided 95 percent lower bound is the 90 percent lower limit.
  one <- do.call(ci_ratio_of_means_paired,
                 c(ml, alpha_lower = 0.05, alpha_upper = 0))
  expect_equal(value_of(one, "ci_lower"), limits_of(p90)[1L])
  expect_identical(value_of(one, "ci_upper"), Inf)
  # The response ratio interval at 90 percent, for comparison.
  lg <- do.call(ci_ratio_of_means_paired, c(ml, conf_level = 0.90,
                                            ci_method = "log"))
  expect_equal(limits_of(lg), c(0.8493114, 0.9054899), tolerance = 1e-6)
})

test_that("a baseline gives the ratio of the per-unit advantages", {
  set.seed(42)
  z <- runif(150, 0, 0.3)
  y <- z + rnorm(150, 0.2, 0.1)
  x <- z + rnorm(150, 0.17, 0.1)
  raw <- ci_ratio_of_means_paired(x, y, baseline = z)
  a <- x - z
  b <- y - z
  smry <- ci_ratio_of_means_paired(mean_1 = mean(a), sd_1 = sd(a),
                                   mean_2 = mean(b), sd_2 = sd(b),
                                   r = cor(a, b), N = 150)
  expect_equal(value_of(raw, c("ci_lower", "ratio", "ci_upper")),
               value_of(smry, c("ci_lower", "ratio", "ci_upper")),
               tolerance = 1e-12)
  expect_equal(value_of(raw, "mean_baseline"), mean(z))
  expect_equal(value_of(raw, "r"), cor(a, b))
  expect_equal(vapply(limits_of(raw), function(L) t.test(a - L * b)$p.value,
                      0), c(0.05, 0.05), tolerance = 1e-9)
  # Noninferiority at f holds exactly for the f below the lower bound.
  one <- ci_ratio_of_means_paired(x, y, baseline = z, alpha_lower = 0.05,
                                  alpha_upper = 0)
  lower <- value_of(one, "ci_lower")
  f <- seq(0.5, 0.99, by = 0.01)
  f <- f[abs(f - lower) > 1e-6]
  rejects <- vapply(f, function(ff) {
    t.test(a - ff * b, alternative = "greater")$p.value < 0.05
  }, logical(1))
  expect_identical(rejects, f < lower)
})

# Units in twelve clusters of three to six, with a cluster effect shared by
# all three variables.
cluster <- rep(1:12, times = c(4, 5, 3, 6, 4, 5, 4, 3, 6, 5, 4, 4))
zc <- c(12.32, 12.62, 11.95, 10.27, 9.6, 7.93, 9.61, 10.68, 11.26, 10.78,
        9.42, 11.82, 8.96, 9.64, 7.92, 8.95, 9.65, 9.05, 10.09, 8.32, 9.33,
        7.84, 12.22, 14.54, 13.3, 12.66, 12.7, 11.65, 11.18, 11.36, 10.68,
        8.53, 9.88, 10.14, 12.1, 11.21, 13.67, 11.08, 11.13, 11.24, 12.7,
        11.38, 10.99, 11.45, 12.1, 11.51, 12.43, 11.87, 10.62, 7.88, 8.98,
        8.41, 9.59)
yc <- c(14.87, 15.58, 16, 13.93, 12.06, 8.89, 11.73, 15.13, 14.19, 15.57,
        13.22, 15.76, 12.95, 11.29, 10.55, 11.9, 12.53, 12.15, 11.96, 10.98,
        12.23, 11.98, 18.06, 17.94, 19.98, 16.53, 15.53, 15.66, 15.58, 13.44,
        12.7, 12.75, 12.92, 11.83, 17.04, 15.71, 18.49, 15.88, 15.63, 15.47,
        14.16, 14.72, 14.35, 15.45, 15.83, 15.08, 16.4, 16.58, 13.16, 8.93,
        9.98, 10.41, 12.07)
xc <- c(13.87, 15.51, 15.6, 12.1, 12.61, 11.62, 11.95, 13.65, 13.16, 15.69,
        12.72, 15.17, 9.83, 12.07, 8.68, 10.25, 12.77, 9.7, 12.37, 9.6,
        10.52, 8.65, 15.43, 16.29, 15.61, 14.51, 17.15, 14.98, 12.16, 14.75,
        12.58, 10.78, 13.3, 12.16, 15.35, 15.27, 17.77, 14.53, 15.95, 16.37,
        17.5, 14.24, 14.35, 13.04, 14.56, 15.27, 14.25, 15.09, 13.05, 11.06,
        10.09, 10.4, 11.79)

test_that("clusters give the cluster-robust covariance on Q - 1 degrees of freedom", {
  out <- ci_ratio_of_means_paired(xc, yc, baseline = zc, cluster = cluster)
  expect_equal(value_of(out, "df"), 11)
  expect_equal(value_of(out, "n_clusters"), 12)
  expect_equal(value_of(out, "N"), 53)
  se_n <- value_of(out, "se_numerator")
  se_d <- value_of(out, "se_denominator")
  # Pinned from sandwich::vcovCL(lm(v ~ 1), cluster = cluster, type = "HC1")
  # for v = x - z, y - z, and their difference (sandwich 3.1.3,
  # 2026-11-05): the variances of the two intercepts and, by
  # var(a - b) = var(a) + var(b) - 2 cov(a, b), their covariance. Live
  # comparison in tools/oracle_snippets/sandwich.R.
  expect_equal(se_n^2, 0.0675418874878823, tolerance = 1e-12)
  expect_equal(se_d^2, 0.0692726405168686, tolerance = 1e-12)
  expect_equal(value_of(out, "r") * se_n * se_d, 0.0469585608693910,
               tolerance = 1e-12)
  expect_equal(value_of(out, c("ci_lower", "ratio", "ci_upper")),
               c(0.6712287747, 0.7961137548, 0.9283755926), tolerance = 1e-9)
  # Each limit inverts the test of the contrast with those variances: the
  # t statistic of mean(a) - L mean(b), over its cluster-robust standard
  # error, is the critical value of t on Q - 1 = 11 degrees of freedom.
  a <- xc - zc
  b <- yc - zc
  va <- 0.0675418874878823
  vb <- 0.0692726405168686
  cab <- 0.0469585608693910
  t_at <- function(L) (mean(a) - L * mean(b)) / sqrt(va - 2 * L * cab +
                                                        L^2 * vb)
  expect_equal(t_at(limits_of(out)), c(1, -1) * qt(0.975, 11),
               tolerance = 1e-10)
  # Without clusters the same data use the per-unit covariance on n - 1.
  plain <- ci_ratio_of_means_paired(xc, yc, baseline = zc)
  expect_equal(value_of(plain, "df"), 52)
  expect_true(is.na(value_of(plain, "n_clusters")))
  expect_lt(diff(limits_of(plain)), diff(limits_of(out)))
  # Cluster labels of any type, in any order, give the same answer.
  shuffled <- ci_ratio_of_means_paired(xc, yc, baseline = zc,
                                       cluster = letters[13 - cluster])
  expect_equal(shuffled$value, out$value)
})

test_that("the log interval adds the covariance term of correlated means", {
  out <- ci_ratio_of_means_paired(xc, yc, ci_method = "log")
  V <- cov(cbind(xc, yc)) / 53
  m <- c(mean(xc), mean(yc))
  se2 <- V[1, 1] / m[1]^2 + V[2, 2] / m[2]^2 - 2 * V[1, 2] / prod(m)
  # Pinned from metafor::escalc("ROMC", correct = FALSE) (metafor 5.0.1,
  # 2026-11-05); live comparison in tools/oracle_snippets/metafor.R.
  expect_equal(log(value_of(out, "ratio")), -0.049602638656155065,
               tolerance = 1e-12)
  expect_equal(se2, 0.000231707617997902, tolerance = 1e-10)
  expect_equal(limits_of(out),
               exp(log(m[1] / m[2]) + c(-1, 1) * qt(0.975, 52) * sqrt(se2)),
               tolerance = 1e-12)
  e <- tryCatch(ci_ratio_of_means_paired(drug_a, drug_b - 5,
                                         ci_method = "log"),
                dmar_refusal = function(e) e)
  expect_s3_class(e, c("dmar_log_ratio_undefined", "dmar_refusal", "error"))
})

test_that("the bootstrap resamples units or clusters, reproducibly", {
  set.seed(1)
  before <- .Random.seed
  a <- ci_ratio_of_means_paired(xc, yc, baseline = zc, cluster = cluster,
                                ci_method = "bca", B = 200, seed = 3)
  expect_identical(.Random.seed, before)
  b <- ci_ratio_of_means_paired(xc, yc, baseline = zc, cluster = cluster,
                                ci_method = "bca", B = 200, seed = 3)
  expect_identical(a$value, b$value)
  expect_identical(a$term, paired_terms)
  expect_identical(attr(a, "B_used"), 200L)
  expect_true(all(is.finite(limits_of(a))))
  # Resampling clusters is resampling their totals: a seed that draws the
  # same cluster indices gives the same replicates.
  tot_a <- rowsum(xc - zc, cluster, reorder = FALSE)[, 1]
  tot_b <- rowsum(yc - zc, cluster, reorder = FALSE)[, 1]
  set.seed(3)
  reps <- DMAR:::.ratio_boot_totals(tot_a, tot_b, 200L)$ratio
  pct <- ci_ratio_of_means_paired(xc, yc, baseline = zc, cluster = cluster,
                                  ci_method = "percentile", B = 200,
                                  seed = 3)
  expect_equal(limits_of(pct), unname(quantile(reps, c(0.025, 0.975))))
})

test_that("replications with a zero denominator are dropped with one warning", {
  xp <- c(3, 2, 4, 3, 5, 2)
  yp <- c(2, 2, 1, 1, 2, 3)
  zp <- c(1, 1, 2, 1, 2, 2)
  expect_warning(
    out <- ci_ratio_of_means_paired(xp, yp, baseline = zp,
                                    ci_method = "percentile", B = 200,
                                    seed = 4),
    "of 200 bootstrap replications \\([0-9.]+%\\) had a denominator of zero",
    class = "dmar_bootstrap_replications_dropped")
  expect_lt(attr(out, "B_used"), 200L)
  expect_identical(attr(out, "n_dropped"), 200L - attr(out, "B_used"))
  expect_output(print(out), "dropped")
  expect_output(print(out), "A dropped replication had a denominator of zero")
})

test_that("a larger paired bootstrap is close to Fieller's interval", {
  skip_on_cran()
  out <- ci_ratio_of_means_paired(xc, yc, baseline = zc, ci_method = "bca",
                                  B = 10000, seed = 5)
  expect_equal(limits_of(out),
               limits_of(ci_ratio_of_means_paired(xc, yc, baseline = zc)),
               tolerance = 0.05)
})

test_that("the table has one shape and tidies to one row", {
  out <- ci_ratio_of_means_paired(drug_a, drug_b, conf_level = 0.995)
  expect_identical(class(out), c("dmar_ci_ratio_of_means", "dmar_ci_long",
                                 "dmar_tbl", "data.frame"))
  expect_identical(out$term, paired_terms)
  td <- generics::tidy(out)
  expect_identical(td$confidence_set, "complement of an interval")
  expect_equal(c(td$excluded_lower, td$excluded_upper),
               value_of(out, c("excluded_lower", "excluded_upper")))
  expect_identical(generics::glance(out), td)
  expect_output(print(out),
                "|t_denominator| does not exceed the critical value.",
                fixed = TRUE)
  smry <- ci_ratio_of_means_paired(mean_1 = 1, sd_1 = 1, mean_2 = 2,
                                   sd_2 = 1, r = 0.5, N = 20)
  expect_identical(smry$term, paired_terms)
  expect_true(all(is.na(value_of(smry, c("mean_baseline", "n_clusters")))))
})

test_that("results_sentence() writes the set as what it is", {
  ci <- function(...) ci_ratio_of_means_paired(drug_a, drug_b, ...)
  # Fieller's (1954, p. 180) three forms.
  expect_identical(results_sentence(ci()),
                   "ratio = 0.32, 95% CI [-0.48, 0.66]")
  complement <- ci(conf_level = 0.995)
  s <- paste0("ratio = 0.32, 99.5% confidence set: every value except ",
              "those between 0.90 and 147.05")
  expect_identical(results_sentence(complement), s)
  expect_identical(results_sentence(generics::tidy(complement)), s)
  expect_identical(results_sentence(complement[1:5, ]), s)
  whole <- ci(conf_level = 0.999)
  s <- paste0("ratio = 0.32, 99.9% confidence set: the whole line, which ",
              "excludes no value")
  expect_identical(results_sentence(whole), s)
  expect_identical(results_sentence(generics::tidy(whole)), s)
  # With unequal tail areas: the union of two intervals, and an interval
  # unbounded above although its upper tail area is not zero, which is not
  # a one-sided bound.
  union <- ci(alpha_lower = 0.02, alpha_upper = 0.002)
  s <- paste0("ratio = 0.32, 97.8% confidence set: the union of [-0.61, ",
              "0.93] and [9.55, Inf) (unequal tail areas: 0.02 below, ",
              "0.002 above)")
  expect_identical(results_sentence(union), s)
  expect_identical(results_sentence(generics::tidy(union)), s)
  open_above <- ci(alpha_lower = 0.025, alpha_upper = 0.0005)
  expect_identical(results_sentence(open_above),
                   paste0("ratio = 0.32, 97.45% CI [-0.48, Inf) (unequal ",
                          "tail areas: 0.025 below, 0.0005 above)"))
  # A one-sided bound is written as one, and a set that is not an interval
  # says that it is one-sided.
  expect_identical(results_sentence(ci(alpha_lower = 0.005, alpha_upper = 0)),
                   "ratio = 0.32, one-sided 99.5% lower confidence bound -2.92")
  expect_identical(results_sentence(ci(alpha_lower = 0, alpha_upper = 0.0005)),
                   paste0("ratio = 0.32, 99.95% confidence set: the whole ",
                          "line, which excludes no value (one-sided: tail ",
                          "area 0.0005 above, none below)"))
  # The sentence of every other table is unchanged.
  expect_identical(results_sentence(ci_smd(smd = 0.5, n_1 = 50, n_2 = 50),
                                    label = "Cohen's d"),
                   "Cohen's d = 0.50, 95% CI [0.10, 0.90]")
})

test_that("the line under the table explains every side left open by the data", {
  ci <- function(...) ci_ratio_of_means_paired(drug_a, drug_b, ...)
  expect_output(print(ci(alpha_lower = 0.025, alpha_upper = 0.0005)),
                paste0("an interval unbounded above. The set is unbounded ",
                       "above because |t_denominator| does not exceed the ",
                       "critical value of the upper tail area."),
                fixed = TRUE)
  expect_output(print(ci(alpha_lower = 0.02, alpha_upper = 0.002)),
                "unbounded above because |t_denominator|", fixed = TRUE)
  # A one-sided bound is unbounded above by design, and nothing is
  # explained.
  one <- utils::capture.output(print(ci(alpha_lower = 0.005,
                                        alpha_upper = 0)))
  expect_true(any(grepl("an interval unbounded above.", one, fixed = TRUE)))
  expect_false(any(grepl("because", one, fixed = TRUE)))
})

test_that("inputs are checked", {
  expect_error(ci_ratio_of_means_paired(drug_a, drug_b[-1]), "same length")
  expect_error(ci_ratio_of_means_paired(drug_a, drug_b, cluster = rep(1, 10)),
               "at least two clusters")
  expect_error(ci_ratio_of_means_paired(mean_1 = 1, sd_1 = 1, mean_2 = 2,
                                        sd_2 = 1, r = 1.2, N = 20),
               "between -1 and 1")
  expect_error(ci_ratio_of_means_paired(mean_1 = 1, sd_1 = 1, mean_2 = 2),
               "all six")
  expect_error(ci_ratio_of_means_paired(drug_a, drug_b, r = 0.5), "not both")
  expect_error(ci_ratio_of_means_paired(drug_a, rep(2, 10)),
               "denominator has no variability")
  e <- tryCatch(ci_ratio_of_means_paired(mean_1 = 1, sd_1 = 1, mean_2 = 2,
                                         sd_2 = 1, r = 0.5, N = 20,
                                         ci_method = "percentile"),
                dmar_refusal = function(e) e)
  expect_s3_class(e, c("dmar_bootstrap_needs_data", "dmar_refusal"))
  # A unit missing any of its values is dropped whole.
  with_na <- ci_ratio_of_means_paired(c(drug_a, NA, 1), c(drug_b, 2, NA))
  expect_equal(with_na$value, ci_ratio_of_means_paired(drug_a, drug_b)$value)
})
