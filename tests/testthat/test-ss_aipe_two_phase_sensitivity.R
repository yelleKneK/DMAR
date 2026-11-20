## ss_aipe_two_phase_sensitivity(): Monte Carlo evaluation of a two-phase
## design. Every block simulates, so every block is skipped on CRAN; the
## closed-form planner tests (test-ss_aipe_two_phase.R) run there. The
## family-wide schema is asserted in test-ss_aipe_sensitivity_family.R.

.v <- function(res) stats::setNames(res$value, res$term)

test_that("coverage is near nominal under normality", {
  skip_on_cran()
  s <- .v(ss_aipe_two_phase_sensitivity(true_sd_Y = 1, true_cor = 0.6,
                                        specified_n_labeled = 150,
                                        n_pool = 5000, width = 0.30,
                                        conf_level = 0.90, G = 4000,
                                        seed = 5))
  # Monte Carlo SE of a .90 coverage at G = 4000 is about .0047.
  expect_lt(abs(s[["total_type_I_error"]] - 0.10), 0.02)
  expect_lt(abs(s[["pct_ci_miss_low"]] - s[["pct_ci_miss_high"]]), 0.02)
  expect_lt(abs(s[["mean_mean_Y"]]), 0.01)
  expect_equal(s[["total_type_I_error"]],
               s[["pct_ci_miss_low"]] + s[["pct_ci_miss_high"]])
})

test_that("drawing the unlabeled proxies through their mean and sum of squares matches drawing every unit", {
  skip_on_cran()
  # The normal generator draws the unlabeled part of the pool from the exact
  # joint distribution of its mean and sum of squares; a full simulation
  # through ci_two_phase() must give the same width distribution.
  n <- 60; M <- 1500; rho <- 0.7; G <- 1500
  set.seed(29)
  w_full <- vapply(seq_len(G), function(g) {
    x <- rnorm(M)
    y <- 2 * (rho * x[1:n] + sqrt(1 - rho^2) * rnorm(n))
    r <- ci_two_phase(y, x[1:n], x[-(1:n)])
    diff(r$value[r$term %in% c("ci_lower", "ci_upper")])
  }, numeric(1))
  s <- .v(ss_aipe_two_phase_sensitivity(true_sd_Y = 2, true_cor = rho,
                                        specified_n_labeled = n, n_pool = M,
                                        width = 1, G = G, seed = 31))
  # Width SD is about 0.06 here, so the Monte Carlo SE of each mean is about
  # 0.0016; a difference of 0.008 is five of them.
  expect_lt(abs(s[["mean_ci_width"]] - mean(w_full)), 0.008)
  expect_lt(abs(s[["sd_ci_width"]] / sd(w_full) - 1), 0.1)
})

test_that("a planned design evaluates the planner's answer and echoes its inputs", {
  skip_on_cran()
  plan <- ss_aipe_two_phase(true_sd_Y = 0.43, true_cor = 0.55, width = 0.10,
                            n_pool = 20000, conf_level = 0.90,
                            assurance = 0.80)
  s <- .v(ss_aipe_two_phase_sensitivity(true_sd_Y = 0.43, true_cor = 0.55,
                                        supposed_sd_Y = 0.43,
                                        supposed_cor = 0.55, n_pool = 20000,
                                        width = 0.10, conf_level = 0.90,
                                        assurance = 0.80, G = 400, seed = 1))
  expect_equal(s[["n_labeled"]],
               plan$value[plan$term == "necessary_n_labeled"])
  expect_equal(s[["total_N"]], 20000)
  expect_equal(c(s[["supposed_sd_Y"]], s[["supposed_cor"]]), c(0.43, 0.55))
  expect_equal(c(s[["assurance"]], s[["true_kurtosis"]]), c(0.80, 0))

  fixed <- .v(ss_aipe_two_phase_sensitivity(true_sd_Y = 0.43, true_cor = 0.55,
                                            specified_n_labeled = 120,
                                            n_pool = 20000, width = 0.10,
                                            G = 20, seed = 1))
  expect_true(is.na(fixed[["supposed_sd_Y"]]))
  expect_true(is.na(fixed[["supposed_cor"]]))
  expect_false(any(c("assurance", "true_kurtosis") %in% names(fixed)))
})

test_that("the realized assurance matches the probability the plan computed", {
  skip_on_cran()
  # In the population it was planned for, the proportion of intervals no
  # wider than the target estimates the probability that the planner's
  # criterion computes at the planned size (at least the assurance asked
  # for, and close to it). The Monte Carlo SE at G = 10000 is about .004, so
  # a miscalibration of a few points fails here.
  planned_prob <- function(n, sd_y, rho, width, M, conf_level) {
    t <- qt(1 - (1 - conf_level) / 2, n - 2)
    room <- (width / (2 * t))^2 - rho^2 * sd_y^2 / M
    pchisq((n - 2) * n * room / (sd_y^2 * (1 - rho^2)), n - 2)
  }
  for (cfg in list(c(0.43, 0.55, 0.10, 20000, 0.90, 0.80),
                   c(1, 0.6, 0.30, 5000, 0.95, 0.90))) {
    s <- .v(ss_aipe_two_phase_sensitivity(
      true_sd_Y = cfg[1], true_cor = cfg[2], supposed_sd_Y = cfg[1],
      supposed_cor = cfg[2], n_pool = cfg[4], width = cfg[3],
      conf_level = cfg[5], assurance = cfg[6], G = 10000, seed = 3))
    p <- planned_prob(s[["n_labeled"]], cfg[1], cfg[2], cfg[3], cfg[4], cfg[5])
    expect_gte(p, cfg[6])
    expect_lt(abs(s[["pct_ci_less_width"]] - p), 0.02, label = paste(cfg, collapse = " "))
  }
})

test_that("a resampled population supplies the true values", {
  skip_on_cran()
  set.seed(7)
  x <- rbinom(400, 1, 0.6)
  y <- ifelse(runif(400) < 0.8, x, rbinom(400, 1, 0.6))
  pop <- data.frame(y = y, x = x)
  s <- .v(ss_aipe_two_phase_sensitivity(specified_n_labeled = 150,
                                        n_pool = 3000, width = 0.10,
                                        population = pop, G = 300, seed = 2))
  expect_equal(s[["true_mean_Y"]], mean(y))
  expect_equal(s[["true_sd_Y"]], sqrt(mean((y - mean(y))^2)))
  expect_equal(s[["true_cor"]], cor(y, x))
  expect_lt(abs(s[["mean_mean_Y"]] - mean(y)), 0.01)
})

test_that("a supplied seed reproduces the run and leaves the caller's stream alone", {
  skip_on_cran()
  call <- function() ss_aipe_two_phase_sensitivity(true_sd_Y = 1, true_cor = 0.5,
                                                   specified_n_labeled = 50,
                                                   n_pool = 1000, width = 0.6,
                                                   G = 25, seed = 113)
  set.seed(1); before <- runif(1)
  set.seed(1); a <- call(); after <- runif(1)
  expect_equal(before, after)
  expect_equal(a$value, call()$value)
})

test_that("replications without a regression interval are dropped with one warning", {
  skip_on_cran()
  # A proxy that is 1 for one unit in ten: with 4 labeled units about two
  # replications in three see a constant proxy.
  set.seed(17)
  pop <- data.frame(y = rnorm(50), x = rep(c(1, 0), c(5, 45)))
  expect_warning(
    res <- ss_aipe_two_phase_sensitivity(specified_n_labeled = 4, n_pool = 200,
                                         width = 2, population = pop, G = 40,
                                         seed = 3),
    "drew labeled units that all share one proxy value")
  expect_null(attr(res, "G_used"))
  expect_identical(attr(res, "simulation")$attempted, 40L)
  expect_lt(attr(res, "simulation")$returned, 40)
  expect_gt(attr(res, "simulation")$returned, 0)
})

test_that("replications whose labeled units share one target value are dropped, not scored at zero width", {
  skip_on_cran()
  # A binary target that is 1 for 96 units in 100: with 10 labeled units
  # about two replications in three see a constant target.
  set.seed(19)
  pop <- data.frame(y = rep(c(0, 1), c(4, 96)), x = rnorm(100))
  expect_warning(
    res <- ss_aipe_two_phase_sensitivity(specified_n_labeled = 10,
                                         n_pool = 500, width = 0.2,
                                         population = pop, G = 60, seed = 8),
    "drew labeled units that all share one target value")
  v <- .v(res)
  expect_lt(attr(res, "simulation")$returned, 60)
  expect_gt(attr(res, "simulation")$returned, 0)
  # No interval of zero width enters the summaries.
  expect_gt(v[["median_ci_width"]], 0)
  # A population in which every draw is constant has nothing to summarize.
  expect_error(ss_aipe_two_phase_sensitivity(specified_n_labeled = 5,
                                             n_pool = 50, width = 0.2,
                                             population = data.frame(
                                               y = rep(c(0, 1), c(1, 99999)),
                                               x = rnorm(1e5)),
                                             G = 3, seed = 1),
               "No replication produced an interval",
               class = "dmar_argument_invalid")
})

test_that("per-replication results are written only to a named file", {
  skip_on_cran()
  expect_null(eval(formals(ss_aipe_two_phase_sensitivity)$filename))
  run <- function(f) ss_aipe_two_phase_sensitivity(true_sd_Y = 1, true_cor = 0.5,
                                                   specified_n_labeled = 40,
                                                   n_pool = 500, width = 0.8,
                                                   G = 3, filename = f,
                                                   seed = 4)
  expect_error(run(42), "'filename' must be NULL or a single character string",
               fixed = TRUE)
  out_file <- tempfile(fileext = ".csv")
  run(out_file)
  written <- utils::read.csv(out_file)
  expect_named(written, c("mean_Y", "ci_lower", "ci_upper", "ci_width",
                          "type_I_lower", "type_I_upper"))
  expect_equal(nrow(written), 3L)
  run(out_file)
  expect_equal(nrow(utils::read.csv(out_file)), 6L)
  unlink(out_file)
})

test_that("contradictory or incomplete specifications are refused", {
  skip_on_cran()
  base <- list(n_pool = 1000, width = 0.5, G = 2)
  call <- function(...) do.call(ss_aipe_two_phase_sensitivity,
                                c(list(...), base))
  expect_error(call(true_sd_Y = 1, true_cor = 0.5), "Supply either")
  expect_error(call(true_sd_Y = 1, true_cor = 0.5, supposed_sd_Y = 1,
                    supposed_cor = 0.5, specified_n_labeled = 50),
               "not both")
  expect_error(call(true_sd_Y = 1, true_cor = 0.5, supposed_sd_Y = 1),
               "Supply both")
  expect_error(call(specified_n_labeled = 50), "true_sd_Y")
  expect_error(call(true_sd_Y = 1, true_cor = 0.5, specified_n_labeled = 50,
                    population = data.frame(y = 1:5, x = c(1, 3, 2, 5, 4))),
               "leave them NULL")
  expect_error(call(specified_n_labeled = 50,
                    population = data.frame(target = 1:5, x = 1:5)),
               "columns 'y'")
  expect_error(call(true_sd_Y = 1, true_cor = 0.5,
                    specified_n_labeled = 1000), "must exceed")
  expect_error(call(true_sd_Y = 1, true_cor = 1.5, specified_n_labeled = 50),
               "true_cor")
})

test_that("a replication without an interval for another reason is counted as such, not put down to too few labeled units", {
  # Integer sizes once left every interval with a missing limit (an
  # overflow, now fixed), and the refusal read "of 5 replications, . The
  # design has too few labeled units for this population." A replication
  # whose interval has a missing limit for a reason other than a constant
  # draw is counted as its own kind. The arguments were accepted and the
  # computation gave no interval, so when no replication is left the
  # refusal has class dmar_computation_failed, not dmar_argument_invalid,
  # which stays the class when every dropped replication was a constant
  # draw.
  local_mocked_bindings(.ci_two_phase_fast = function(...)
    c(mean_Y = NA_real_, se = NA_real_, df = NA_real_, ci_lower = NA_real_,
      ci_upper = NA_real_, slope = NA_real_))
  err <- expect_error(
    ss_aipe_two_phase_sensitivity(true_sd_Y = 1, true_cor = 0.7,
                                  specified_n_labeled = 50, n_pool = 500,
                                  width = 0.5, G = 5, seed = 1),
    class = "dmar_computation_failed")
  expect_false(inherits(err, "dmar_argument_invalid"))
  expect_identical(conditionMessage(err), paste0(
    "No replication produced an interval: of 5 replications, 5 gave an ",
    "interval with a missing limit."))
  # Constant draws and missing limits together: both are counted, and the
  # refusal is not put down to the number of labeled units. A binary target
  # that is 1 for 87 units in 100 gives five labeled units that all share
  # one value in about half the replications; the mock gives the others a
  # missing limit.
  set.seed(2)
  pop <- data.frame(y = rep(c(0, 1), c(13, 87)), x = rnorm(100))
  err <- expect_error(
    ss_aipe_two_phase_sensitivity(specified_n_labeled = 5, n_pool = 50,
                                  width = 0.5, population = pop, G = 10,
                                  seed = 1),
    class = "dmar_computation_failed")
  msg <- conditionMessage(err)
  expect_match(msg, "drew labeled units that all share one target value",
               fixed = TRUE)
  expect_match(msg, "gave an interval with a missing limit.", fixed = TRUE)
  expect_false(grepl("too few labeled units", msg, fixed = TRUE))
})

test_that("a replication with a missing limit is counted in the warning when others are kept", {
  calls <- 0L
  local_mocked_bindings(.ci_two_phase_fast = function(...) {
    calls <<- calls + 1L
    v <- c(mean_Y = 0, se = 0.1, df = 48, ci_lower = -0.2, ci_upper = 0.2,
           slope = 0.5)
    if (calls %% 2L == 1L) v[c("ci_lower", "ci_upper")] <- NA_real_
    v
  })
  expect_warning(
    res <- ss_aipe_two_phase_sensitivity(true_sd_Y = 1, true_cor = 0.7,
                                         specified_n_labeled = 50,
                                         n_pool = 500, width = 0.5, G = 4,
                                         seed = 1),
    paste0("Of 4 replications, 2 gave an interval with a missing limit; ",
           "they were dropped, and the summaries use the remaining 2."),
    fixed = TRUE, class = "dmar_result_qualified")
  expect_equal(attr(res, "simulation")$returned, 2)
})
