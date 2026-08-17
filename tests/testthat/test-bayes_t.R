v <- function(tab, t) tab$value[tab$term == t]

test_that("bayes_one_sample_t() returns the documented schema with coherent values", {
  set.seed(113)
  x <- rnorm(40, 0.4, 1)
  res <- bayes_one_sample_t(x)
  expect_s3_class(res, "dmar_tbl")
  expect_true(all(c("delta_posterior_median", "delta_posterior_mean",
                    "delta_lower", "delta_upper", "p_delta_positive",
                    "raw_posterior_median", "bf_10", "bf_01", "t", "df",
                    "prior_scale", "n") %in% res$term))
  # Coherence: median inside the credible interval, bf_01 = 1/bf_10,
  # the posterior shrinks the sample effect toward zero, and the raw rows
  # are the delta rows times the sample SD.
  expect_gt(v(res, "delta_posterior_median"), v(res, "delta_lower"))
  expect_lt(v(res, "delta_posterior_median"), v(res, "delta_upper"))
  expect_equal(v(res, "bf_01"), 1 / v(res, "bf_10"))
  d_hat <- v(res, "t") / sqrt(v(res, "n"))
  expect_lt(abs(v(res, "delta_posterior_median")), abs(d_hat))
  expect_equal(v(res, "raw_posterior_median"),
               v(res, "delta_posterior_median") * sd(x))
  expect_gt(v(res, "p_delta_positive"), 0.5)
})

test_that("the JZS Bayes factor matches pinned ttestBF values in all three designs", {
  set.seed(113)
  x <- rnorm(35, 0.3, 1); y <- rnorm(40, 0, 1)

  # Pinned from BayesFactor::ttestBF (BayesFactor 0.9.12.4.8, 2026-08-09);
  # live comparison in tools/oracle_checks.R.
  bf_one <- 2.997079992617161
  expect_equal(v(bayes_one_sample_t(x), "bf_10"), bf_one, tolerance = 1e-4)

  bf_two <- 0.6494000559300029
  expect_equal(v(bayes_independent_t(x, y), "bf_10"), bf_two,
               tolerance = 1e-4)

  y2 <- x + rnorm(35, 0.2, 0.8)
  bf_pair <- 0.3103826542565354
  expect_equal(v(bayes_paired_t(y2, x), "bf_10"), bf_pair, tolerance = 1e-4)
})

test_that("the posterior summaries match pinned MCMC values within Monte Carlo error", {
  set.seed(113)
  x <- rnorm(50, 0.5, 1)
  res <- bayes_one_sample_t(x)
  # Pinned from 40000 draws of BayesFactor::posterior under the seed above
  # (BayesFactor 0.9.12.4.8, 2026-08-09); live comparison in
  # tools/oracle_checks.R.
  expect_equal(v(res, "delta_posterior_median"), 0.6501664738341795,
               tolerance = 0.02)
  expect_equal(v(res, "delta_posterior_mean"), 0.6518518146352267,
               tolerance = 0.02)
  expect_equal(v(res, "delta_lower"), 0.3500982113270965,
               tolerance = 0.04)
  expect_equal(v(res, "delta_upper"), 0.962897091299841,
               tolerance = 0.04)
  expect_equal(v(res, "p_delta_positive"), 1,
               tolerance = 0.01)
})

test_that("the prior scale behaves as it should", {
  set.seed(113)
  x <- rnorm(25, 0.2, 1)
  med  <- bayes_one_sample_t(x, prior_scale = sqrt(2) / 2)
  wide <- bayes_one_sample_t(x, prior_scale = 2)
  # A wider prior spreads mass over big effects the data rule out, so it
  # penalizes H1: the Bayes factor for the effect drops.
  expect_lt(v(wide, "bf_10"), v(med, "bf_10"))
})

test_that("paired analysis equals the one-sample analysis of differences", {
  set.seed(113)
  a <- rnorm(20, 10, 2); b <- rnorm(20, 9, 2)
  expect_equal(bayes_paired_t(a, b)$value,
               bayes_one_sample_t(a - b)$value)
})

test_that("the bayes_* t functions validate their arguments", {
  expect_error(bayes_one_sample_t(c(1, NA, 2)), "non-missing")
  expect_error(bayes_one_sample_t(1), "at least 2")
  expect_error(bayes_one_sample_t(rnorm(10), prior_scale = -1), "positive")
  expect_error(bayes_one_sample_t(rnorm(10), conf_level = 1.1),
               "\\(0, 1\\)")
  expect_error(bayes_paired_t(rnorm(10), rnorm(11)), "same length")
})
