## ss_aipe_mean_sensitivity(): the Monte Carlo sibling of ss_aipe_mean().
## The family-wide return schema is asserted in
## test-ss_aipe_sensitivity_family.R; this file covers the member's own
## arguments (the resampled population, the seed, the file), its echoes,
## and, behind skip_on_cran(), what the simulation should find: a correct
## normal-theory plan delivers its expected width, its assurance, and its
## coverage, and the kurtosis adjustment delivers the assurance for a
## population whose tails are not normal.

sens_value <- function(res, term) res$value[res$term == term]

test_that("ss_aipe_mean_sensitivity() refuses inconsistent or malformed inputs before simulating", {
  expect_error(ss_aipe_mean_sensitivity(true_sd = 1, width = 0.5, G = 2),
               "exactly one of 'supposed_sd'")
  expect_error(ss_aipe_mean_sensitivity(true_sd = 1, supposed_sd = 1,
                                        specified_N = 30, width = 0.5, G = 2),
               "exactly one of 'supposed_sd'")
  expect_error(ss_aipe_mean_sensitivity(supposed_sd = 1, width = 0.5, G = 2),
               "Supply 'true_sd'")
  expect_error(ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 30,
                                        width = 0.5, population = 1:10, G = 2),
               "not both")
  expect_error(ss_aipe_mean_sensitivity(specified_N = 30, width = 0.5,
                                        population = 1:10, true_mean = 5,
                                        G = 2),
               "'true_mean' describes a normal population")
  expect_error(ss_aipe_mean_sensitivity(specified_N = 30, width = 0.5,
                                        population = c(1, NA, 3), G = 2),
               "missing or infinite")
  expect_error(ss_aipe_mean_sensitivity(specified_N = 30, width = 0.5,
                                        population = rep(4, 10), G = 2),
               "constant")
  expect_error(ss_aipe_mean_sensitivity(true_sd = -1, specified_N = 30,
                                        width = 0.5, G = 2), "'true_sd'")
  expect_error(ss_aipe_mean_sensitivity(true_sd = 1, supposed_sd = 0,
                                        width = 0.5, G = 2), "'supposed_sd'")
  expect_error(ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 30.5,
                                        width = 0.5, G = 2), "'specified_N'")
  expect_error(ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 1,
                                        width = 0.5, G = 2), "'specified_N'")
  expect_error(ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 30,
                                        width = 0.5, G = 0), "'G'")
  expect_error(ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 30,
                                        width = 0.5, G = 2, filename = 42),
               "'filename' must be NULL or a single character string")
  expect_error(ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 30,
                                        width = 0.5, G = 2, seed = "a"),
               "'seed'")
  # The kurtosis is checked as the planner checks it, whether or not the
  # planner is called.
  expect_error(ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 30,
                                        width = 0.5, true_kurtosis = -2, G = 2),
               "-2 or less")
  expect_error(ss_aipe_mean_sensitivity(true_sd = 1, supposed_sd = 1,
                                        width = 0.5, assurance = 0.9,
                                        true_kurtosis = -2.1, G = 2),
               "-2 or less")
  # A kurtosis that cannot reach the plan is reported, not silently dropped.
  expect_warning(ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 30,
                                          width = 0.5, true_kurtosis = 2, G = 2),
                 "unused when 'specified_N' is supplied")
})

test_that("the echoes report the plan's size and the population that generated the data", {
  planned <- ss_aipe_mean(true_sd = 2, width = 1, assurance = 0.80)
  res <- ss_aipe_mean_sensitivity(true_sd = 2.5, supposed_sd = 2, width = 1,
                                  true_mean = 10, assurance = 0.80, G = 3,
                                  seed = 113)
  expect_equal(sens_value(res, "total_N"),
               planned$value[planned$term == "necessary_N"])
  expect_equal(sens_value(res, "true_sd"), 2.5)
  expect_equal(sens_value(res, "true_mean"), 10)
  expect_equal(sens_value(res, "supposed_sd"), 2)
  expect_equal(sens_value(res, "assurance"), 0.80)
  expect_equal(attr(res, "population"), "normal")
  expect_equal(attr(res, "conf_level"), 0.95)

  # A specified size echoes NA for the planning value, and no assurance row
  # appears when none was supplied.
  res_n <- ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 40, width = 1,
                                    G = 3, seed = 113)
  expect_equal(sens_value(res_n, "total_N"), 40)
  expect_true(is.na(sens_value(res_n, "supposed_sd")))
  expect_false("assurance" %in% res_n$term)

  # With a resampled population the true values are those of the
  # distribution resampled: its mean, and its standard deviation with
  # divisor equal to the number of scores.
  pop <- c(0, 0, 1, 3, 7)
  res_p <- ss_aipe_mean_sensitivity(specified_N = 10, width = 4,
                                    population = pop, G = 3, seed = 113)
  expect_equal(sens_value(res_p, "true_mean"), mean(pop))
  expect_equal(sens_value(res_p, "true_sd"), sqrt(mean((pop - mean(pop))^2)))
  expect_equal(attr(res_p, "population"), "resampled")
})

test_that("a seed reproduces the run and leaves the caller's random stream as it was", {
  set.seed(42)
  before <- .Random.seed
  a <- ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 25, width = 0.8,
                                G = 20, seed = 113)
  expect_identical(.Random.seed, before)
  b <- ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 25, width = 0.8,
                                G = 20, seed = 113)
  expect_identical(a, b)
  # Without a seed the draws come from, and advance, the caller's stream.
  set.seed(7)
  before_unseeded <- .Random.seed
  c1 <- ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 25, width = 0.8,
                                 G = 20)
  expect_false(identical(.Random.seed, before_unseeded))
  set.seed(7)
  c2 <- ss_aipe_mean_sensitivity(true_sd = 1, specified_N = 25, width = 0.8,
                                 G = 20)
  expect_identical(c1, c2)
})

test_that("the per-replication file holds the replications the summary was computed from", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  res <- ss_aipe_mean_sensitivity(specified_N = 12, width = 6,
                                  population = pygmalion$iq_gain, G = 7,
                                  filename = path, seed = 113)
  per_rep <- utils::read.csv(path)
  expect_named(per_rep, c("mean_hat", "sd_hat", "ci_lower", "ci_upper",
                          "ci_width", "type_I_lower", "type_I_upper"))
  expect_equal(nrow(per_rep), 7L)
  expect_equal(mean(per_rep$mean_hat), sens_value(res, "mean_mean"))
  expect_equal(mean(per_rep$ci_width), sens_value(res, "mean_ci_width"))
  expect_equal(mean(per_rep$ci_width <= 6), sens_value(res, "pct_ci_less_width"))
  expect_equal(mean(per_rep$type_I_lower), sens_value(res, "pct_ci_miss_low"))
  expect_equal(mean(per_rep$type_I_upper), sens_value(res, "pct_ci_miss_high"))
  # The intervals are the t intervals at the size evaluated.
  half <- stats::qt(0.975, 11) * per_rep$sd_hat / sqrt(12)
  expect_equal(per_rep$ci_upper - per_rep$ci_lower, 2 * half)

  # A second run appends to the same file.
  ss_aipe_mean_sensitivity(specified_N = 12, width = 6,
                           population = pygmalion$iq_gain, G = 7,
                           filename = path, seed = 113)
  expect_equal(nrow(utils::read.csv(path)), 14L)
})

test_that("a correct normal-theory plan delivers its expected width, its assurance, and its coverage", {
  skip_on_cran()
  plan <- ss_aipe_mean(true_sd = 1, width = 0.4)
  res <- ss_aipe_mean_sensitivity(true_sd = 1, supposed_sd = 1, width = 0.4,
                                  G = 4000, seed = 20261105)
  # The mean realized width estimates the exact expected width the planner
  # reports; its Monte Carlo standard error is sd_ci_width / sqrt(G).
  se_width <- sens_value(res, "sd_ci_width") / sqrt(4000)
  expect_lt(abs(sens_value(res, "mean_ci_width") -
                  plan$value[plan$term == "expected_width"]), 4 * se_width)
  # Each tail misses at alpha / 2 under normality (Monte Carlo SE about
  # .0025 at G = 4000).
  expect_lt(abs(sens_value(res, "pct_ci_miss_low") - 0.025), 0.012)
  expect_lt(abs(sens_value(res, "pct_ci_miss_high") - 0.025), 0.012)

  res_a <- ss_aipe_mean_sensitivity(true_sd = 1, supposed_sd = 1,
                                    width = 0.4, assurance = 0.80, G = 4000,
                                    seed = 20261105)
  # Monte Carlo SE of the attained assurance about .0063 at G = 4000.
  expect_lt(abs(sens_value(res_a, "pct_ci_less_width") - 0.80), 0.03)

  # A planning value below the truth widens every interval in proportion:
  # with the same seed the normal draws are the same deviates scaled by
  # true_sd, so the ratio of the mean widths is exactly the ratio of the
  # standard deviations.
  res_m <- ss_aipe_mean_sensitivity(true_sd = 1.2, supposed_sd = 1,
                                    width = 0.4, G = 4000, seed = 20261105)
  expect_equal(sens_value(res_m, "mean_ci_width") / sens_value(res, "mean_ci_width"),
               1.2, tolerance = 1e-12)
  expect_lt(sens_value(res_m, "pct_ci_less_width"), sens_value(res, "pct_ci_less_width"))
})

test_that("the kurtosis adjustment delivers the assurance for a population with heavy tails", {
  skip_on_cran()
  # An exponential population, represented by its quantiles: excess
  # kurtosis about 6. The normal-theory plan undershoots its assurance;
  # the plan made with the population's kurtosis attains it within Monte
  # Carlo error of the two-moment approximation.
  pop   <- stats::qexp(stats::ppoints(20001))
  sigma <- sqrt(mean((pop - mean(pop))^2))
  kappa <- mean(((pop - mean(pop)) / sigma)^4) - 3
  normal_plan <- ss_aipe_mean_sensitivity(supposed_sd = sigma, width = 0.4,
                                          population = pop, assurance = 0.80,
                                          G = 4000, seed = 20261105)
  kurtosis_plan <- ss_aipe_mean_sensitivity(supposed_sd = sigma, width = 0.4,
                                            population = pop,
                                            assurance = 0.80, true_kurtosis = kappa,
                                            G = 4000, seed = 20261105)
  expect_gt(sens_value(kurtosis_plan, "total_N"),
            sens_value(normal_plan, "total_N"))
  expect_lt(sens_value(normal_plan, "pct_ci_less_width"), 0.76)
  expect_lt(abs(sens_value(kurtosis_plan, "pct_ci_less_width") - 0.80), 0.035)
})
