## .ss_sens_simulation() and .ss_sens_record(): the record of a
## sensitivity function's Monte Carlo run that every *_sensitivity()
## function attaches as attr(out, "simulation"). One denominator: the
## width proportion is over the replications with an interval.

sim <- function(...) DMAR:::.ss_sens_simulation(...)
sim_names <- c("attempted", "returned", "failed", "denominator",
               "width_target_probability", "width_target_mc_ci")

test_that("every replication with an interval: the denominator is G", {
  met <- rep(c(TRUE, FALSE), c(37L, 13L))
  s <- sim(50, rep(TRUE, 50), met, 0.95)
  expect_identical(names(s), sim_names)
  expect_identical(s$attempted, 50L)
  expect_identical(s$returned, 50L)
  expect_identical(s$failed, 0L)
  expect_identical(s$denominator, 50L)
  expect_identical(s$width_target_probability, 37 / 50)
  expect_identical(s$width_target_mc_ci,
                   unname(stats::binom.test(37, 50, conf.level = 0.95)$conf.int))
})

test_that("the Monte Carlo limits are the exact binomial limits at conf_level", {
  s <- sim(50, rep(TRUE, 50), rep(c(TRUE, FALSE), c(37L, 13L)), 0.90)
  # Clopper-Pearson limits from the beta quantiles, written out.
  expect_equal(as.numeric(s$width_target_mc_ci),
               c(stats::qbeta(.05, 37, 14), stats::qbeta(.95, 38, 13)),
               tolerance = 1e-12)
  none <- sim(20, rep(TRUE, 20), rep(FALSE, 20), 0.95)
  expect_equal(as.numeric(none$width_target_mc_ci),
               c(0, 1 - .025^(1 / 20)), tolerance = 1e-12)
  all_met <- sim(20, rep(TRUE, 20), rep(TRUE, 20), 0.95)
  expect_equal(as.numeric(all_met$width_target_mc_ci),
               c(.025^(1 / 20), 1), tolerance = 1e-12)
})

test_that("a replication without an interval is left out of the width proportion and counted", {
  has <- c(TRUE, TRUE, FALSE, FALSE)
  # NA where there is no interval (NA <= width), and TRUE by mistake: the
  # entry of a replication without an interval is not read.
  for (met in list(c(TRUE, FALSE, NA, NA), c(TRUE, FALSE, TRUE, TRUE))) {
    s <- sim(4L, has, met, 0.95)
    expect_identical(s$attempted, 4L)
    expect_identical(s$returned, 2L)
    expect_identical(s$failed, 2L)
    expect_identical(s$denominator, 2L)
    expect_identical(s$width_target_probability, 0.5)
    # The limits are over the two replications with an interval, not G.
    expect_identical(s$width_target_mc_ci,
                     unname(stats::binom.test(1, 2)$conf.int))
  }
  # NA in met_width for a replication with an interval counts as not met.
  expect_identical(
    sim(3, rep(TRUE, 3), c(TRUE, NA, TRUE), 0.95)$width_target_probability,
    2 / 3)
})

test_that("no interval in any replication: the width entries are NA", {
  s <- sim(3, rep(FALSE, 3), rep(NA, 3), 0.95)
  expect_identical(s$attempted, 3L)
  expect_identical(s$returned, 0L)
  expect_identical(s$failed, 3L)
  expect_identical(s$denominator, 0L)
  expect_identical(s$width_target_probability, NA_real_)
  expect_identical(s$width_target_mc_ci, c(NA_real_, NA_real_))
})

test_that("a function with no width target has NA width entries", {
  s <- sim(200, rep(c(TRUE, FALSE), c(198L, 2L)), NULL, 0.95)
  expect_identical(names(s), sim_names)
  expect_identical(s$attempted, 200L)
  expect_identical(s$returned, 198L)
  expect_identical(s$failed, 2L)
  expect_identical(s$denominator, 198L)
  expect_identical(s$width_target_probability, NA_real_)
  expect_identical(s$width_target_mc_ci, c(NA_real_, NA_real_))
})

test_that("arguments of the wrong shape stop", {
  expect_error(sim(4, rep(TRUE, 3), rep(TRUE, 4), 0.95), "has_interval")
  expect_error(sim(4, c(TRUE, NA, TRUE, TRUE), rep(TRUE, 4), 0.95),
               "has_interval")
  expect_error(sim(4, c(1, 1, 1, 1), rep(TRUE, 4), 0.95), "has_interval")
  expect_error(sim(4, rep(TRUE, 4), rep(TRUE, 3), 0.95), "met_width")
  expect_error(sim(4, rep(TRUE, 4), c(1, 0, 1, 1), 0.95), "met_width")
  expect_error(sim(2.5, rep(TRUE, 2), rep(TRUE, 2), 0.95),
               class = "dmar_count_invalid")
})

test_that(".ss_sens_record() reads the width count from the table's own row", {
  record <- function(pct, G, has, conf_level = 0.95) {
    tab <- data.frame(term = c("mean_ci_width", "pct_ci_less_width"),
                      value = c(1, pct))
    attr(DMAR:::.ss_sens_record(tab, G, has, conf_level), "simulation")
  }
  has <- rep(c(TRUE, FALSE), c(8L, 2L))
  # Six of the eight intervals met the target: the same list as from the
  # vectors, whatever the two replications without an interval hold.
  expect_identical(record(6 / 8, 10, has),
                   sim(10, has, seq_len(10) <= 6, 0.95))
  expect_identical(record(0, 10, has)$width_target_probability, 0)
  expect_identical(record(1, 10, has)$width_target_probability, 1)
  # A row over all ten attempted replications is not a whole number of the
  # eight with an interval, which breaks the rule and stops.
  expect_error(record(6 / 10, 10, has), "replications with an interval")
  # No width target (NA or NaN in the row) and no row: the counts alone.
  for (pct in c(NA_real_, NaN)) {
    s <- record(pct, 10, has)
    expect_identical(s, sim(10, has, NULL, 0.95))
  }
  none <- data.frame(term = "coverage", value = .9)
  s <- attr(DMAR:::.ss_sens_record(none, 10, rep(TRUE, 10)), "simulation")
  expect_identical(s, sim(10, rep(TRUE, 10), NULL, 0.95))
  # No replication with an interval.
  expect_identical(record(NA_real_, 3, rep(FALSE, 3)),
                   sim(3, rep(FALSE, 3), rep(NA, 3), 0.95))
  # The table is returned with its rows and other attributes as they were.
  tab <- data.frame(term = "pct_ci_less_width", value = .5)
  attr(tab, "conf_level") <- 0.9
  out <- DMAR:::.ss_sens_record(tab, 4, rep(TRUE, 4), 0.9)
  expect_identical(out$value, .5)
  expect_identical(attr(out, "conf_level"), 0.9)
  expect_error(DMAR:::.ss_sens_record(tab, 4, rep(TRUE, 3), 0.9),
               "has_interval")
  expect_error(DMAR:::.ss_sens_record(tab, 0, logical(0), 0.9),
               class = "dmar_count_invalid")
})

test_that("ss_aipe_reliability_sensitivity() attaches the helper's list", {
  skip_on_cran()
  set.seed(113)
  out <- ss_aipe_reliability_sensitivity(true_reliability = .8,
    specified_N = 60, n_items = 4, width = .15, G = 40)
  s <- attr(out, "simulation")
  expect_identical(names(s), append(sim_names, "first_failure", after = 5L))
  expect_identical(s$attempted, 40L)
  expect_identical(s$failed, 0L)
  p <- out$value[out$term == "pct_ci_less_width"]
  expect_identical(s[sim_names],
                   sim(40, rep(TRUE, 40), seq_len(40) <= round(40 * p), 0.95))
})
