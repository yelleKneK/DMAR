# ss_aipe_cliff_delta_sensitivity(): replications without an interval, the
# record of the simulation, and the seed. (The family contract is in
# test-ss_aipe_sensitivity_family.R and the recovery of the estimate in
# test-ss_aipe_new_sensitivities.R.)

# Runs `expr`, collecting the warnings of class
# dmar_replications_without_interval and muffling them.
without_interval <- function(expr) {
  seen <- new.env(parent = emptyenv())
  seen$w <- list()
  value <- withCallingHandlers(expr,
    dmar_replications_without_interval = function(w) {
      seen$w[[length(seen$w) + 1L]] <- w
      invokeRestart("muffleWarning")
    })
  list(value = value, warnings = seen$w)
}
row <- function(res, term) res$value[res$term == term]

test_that("a completely separated replication has no interval and is in no width or coverage row", {
  skip_on_cran()
  # In DMAR 1.0.x this call returned pct_ci_less_width = .405, all of it from
  # replications whose two limits were at the estimate (width 0).
  set.seed(113)
  run <- expect_no_warning(without_interval(
    ss_aipe_cliff_delta_sensitivity(true_delta = 0.9, specified_N = 12,
                                    width = 0.4, G = 400)),
    class = "dmar_interval_degenerate")
  res <- run$value
  expect_length(run$warnings, 1L)
  w <- run$warnings[[1L]]
  expect_identical(w$n_without_interval, 162L)
  expect_identical(w$G, 400L)
  expect_match(conditionMessage(w),
               "162 of 400 replications had completely separated groups",
               fixed = TRUE)
  expect_match(conditionMessage(w),
               "describe the other 238 replications", fixed = TRUE)
  expect_identical(row(res, "pct_ci_less_width"), 0)
  sim <- attr(res, "simulation")
  expect_identical(sim$attempted, 400L)
  expect_identical(sim$failed, 162L)
  expect_identical(sim$returned, 238L)
  expect_identical(sim$denominator, 238L)
  expect_identical(sim$width_target_probability, 0)
  expect_equal(sim$width_target_mc_ci,
               as.numeric(stats::binom.test(0, 238)$conf.int),
               ignore_attr = TRUE, tolerance = 1e-12)
  # The estimates of all 400 replications enter the first three rows.
  expect_gt(row(res, "mean_cliff_delta"), 0.85)
  # Every interval that exists is wider than the target here.
  expect_gt(row(res, "mean_ci_width"), 0.4)
  expect_true(all(is.finite(res$value[res$term != "supposed_delta"])))
})

test_that("the rows and the per-replication file agree, with NA where there is no interval", {
  skip_on_cran()
  file <- tempfile(fileext = ".csv")
  on.exit(unlink(file), add = TRUE)
  run <- without_interval(
    ss_aipe_cliff_delta_sensitivity(true_delta = 0.85, specified_N = 14,
                                    width = 1.2, G = 300, filename = file,
                                    seed = 2718))
  res <- run$value
  rep <- utils::read.csv(file)
  expect_identical(nrow(rep), 300L)
  none <- is.na(rep$ci_width)
  expect_identical(sum(none), attr(res, "simulation")$failed)
  expect_gt(sum(none), 0L)
  # Without an interval: the estimate is at the end of its range and the
  # limits, width, and indicators are missing.
  expect_true(all(abs(rep$cliff_delta[none]) == 1))
  expect_true(all(is.na(rep$ci_lower[none]) & is.na(rep$ci_upper[none]) &
                    is.na(rep$type_I_lower[none]) &
                    is.na(rep$type_I_upper[none])))
  expect_true(all(abs(rep$cliff_delta[!none]) < 1))
  # The rows from the file.
  expect_equal(row(res, "mean_cliff_delta"), mean(rep$cliff_delta),
               tolerance = 1e-12)
  expect_equal(row(res, "sd_cliff_delta"), stats::sd(rep$cliff_delta),
               tolerance = 1e-12)
  expect_equal(row(res, "mean_ci_width"), mean(rep$ci_width[!none]),
               tolerance = 1e-12)
  expect_equal(row(res, "median_ci_width"),
               stats::median(rep$ci_width[!none]), tolerance = 1e-12)
  expect_equal(row(res, "pct_ci_less_width"),
               mean(rep$ci_width[!none] <= 1.2), tolerance = 1e-12)
  expect_gt(row(res, "pct_ci_less_width"), 0)
  expect_equal(row(res, "pct_ci_miss_low"), mean(rep$type_I_lower[!none]),
               tolerance = 1e-12)
  expect_equal(row(res, "total_type_I_error"),
               mean(rep$type_I_lower[!none] | rep$type_I_upper[!none]),
               tolerance = 1e-12)
  sim <- attr(res, "simulation")
  expect_identical(sim$denominator, sum(!none))
  expect_equal(sim$width_target_probability,
               row(res, "pct_ci_less_width"), tolerance = 1e-12)
  expect_equal(sim$width_target_mc_ci,
               as.numeric(stats::binom.test(sum(rep$ci_width[!none] <= 1.2),
                                            sum(!none))$conf.int),
               ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("a run in which every replication has an interval is as it was, with the record attached", {
  set.seed(113)
  res <- expect_no_warning(ss_aipe_cliff_delta_sensitivity(
    true_delta = 0.30, supposed_delta = 0.30, width = 0.50, G = 30))
  sim <- attr(res, "simulation")
  expect_identical(c(sim$attempted, sim$returned, sim$failed),
                   c(30L, 30L, 0L))
  expect_equal(sim$width_target_probability,
               row(res, "pct_ci_less_width"), tolerance = 1e-12)
  expect_identical(res$term,
                   c("mean_cliff_delta", "median_cliff_delta",
                     "sd_cliff_delta", "mean_ci_width", "median_ci_width",
                     "sd_ci_width", "pct_ci_less_width", "pct_ci_miss_low",
                     "pct_ci_miss_high", "total_type_I_error", "n_1", "n_2",
                     "total_N", "true_delta", "supposed_delta", "ratio",
                     "width", "conf_level"))
})

test_that("seed is the last formal, makes the call reproducible, and leaves the caller's stream alone", {
  skip_on_cran()
  f <- formals(ss_aipe_cliff_delta_sensitivity)
  expect_identical(names(f)[length(f)], "seed")
  expect_null(f$seed)
  call <- function(seed = NULL) {
    ss_aipe_cliff_delta_sensitivity(true_delta = 0.3, specified_N = 40,
                                    width = 0.5, G = 20, seed = seed)
  }
  set.seed(99)
  before <- .Random.seed
  first <- call(seed = 113)
  expect_identical(.Random.seed, before)
  expect_identical(call(seed = 113), first)
  expect_false(identical(call(seed = 114)$value, first$value))
  # seed = NULL draws from the caller's stream, as the function always did.
  set.seed(113)
  from_stream <- call()
  expect_false(identical(.Random.seed, {set.seed(113); .Random.seed}))
  expect_identical(from_stream$value, first$value)
  expect_error(call(seed = "a"), "'seed' must be a single number or NULL.",
               fixed = TRUE)
})

test_that("the table prints the number of replications, and the three counts when some had no interval", {
  out <- ss_aipe_cliff_delta_sensitivity(true_delta = 0.3, specified_N = 60,
                                         width = 0.5, G = 10, seed = 113)
  lines <- DMAR:::.dmar_tbl_footer(out)
  expect_identical(unname(lines["simulation"]), "Replications: 10")
  expect_output(print(out), "Replications: 10", fixed = TRUE)
  sep <- suppressWarnings(ss_aipe_cliff_delta_sensitivity(
    true_delta = 0.9, specified_N = 12, width = 0.4, G = 40, seed = 113))
  n_failed <- attr(sep, "simulation")$failed
  expect_gt(n_failed, 0L)
  expect_identical(unname(DMAR:::.dmar_tbl_footer(sep)["simulation"]),
                   sprintf(paste0("Replications: 40 attempted, %d with an ",
                                  "interval, %d without"),
                           40L - n_failed, n_failed))
})

test_that("a run in which every replication has an interval gives no warning", {
  expect_no_warning(ss_aipe_cliff_delta_sensitivity(
    true_delta = 0.3, specified_N = 60, width = 0.5, G = 10, seed = 113))
})
