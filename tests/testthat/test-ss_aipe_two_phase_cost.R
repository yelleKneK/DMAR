## ss_aipe_two_phase_cost(): the least costly whole-number design for the
## expected width of the ci_two_phase() interval, checked against an
## exhaustive search over pairs (n, M), against fixed anchors, against the
## companion planner ss_aipe_two_phase(), and at the edges.

.v <- function(res) stats::setNames(res$value, res$term)

# The planned width of ci_two_phase() at the planning values, and the
# expected width of the t interval for one mean.
.tp_width <- function(n, M, sd_y, rho, q)
  2 * qt(q, n - 2) * sd_y * sqrt((1 - rho^2) / n + rho^2 / M)
.tp_width_target <- function(n, sd_y, q)
  2 * qt(q, n - 1) * sqrt(2 / (n - 1)) *
    exp(lgamma(n / 2) - lgamma((n - 1) / 2)) * sd_y / sqrt(n)

# An exhaustive search that uses neither the normal-quantile bound nor its
# convexity: every n from 3 for as long as n labeled units with the
# smallest admissible pool, n + 1, would cost less than the best design so
# far, and for each n the smallest pool that meets the width (the width
# falls as the pool grows, so the first M that meets it is found by
# bisection between n + 1 and the largest pool that could still cost
# less). Returns c(n, M, cost); M = 0 is the target alone.
.tp_brute <- function(sd_y, rho, width, r, conf_level) {
  q <- 1 - (1 - conf_level) / 2
  n_t <- 2
  while (.tp_width_target(n_t, sd_y, q) > width) n_t <- n_t + 1
  best <- c(n_t, 0, n_t)
  n <- 3
  while (n + r * (n + 1) < best[3]) {
    M_top <- floor((best[3] - n) / r) + 1
    if (M_top >= n + 1 && .tp_width(n, M_top, sd_y, rho, q) <= width) {
      lo <- n + 1
      hi <- M_top
      if (.tp_width(n, lo, sd_y, rho, q) <= width) {
        hi <- lo
      } else {
        while (hi - lo > 1) {
          mid <- floor((lo + hi) / 2)
          if (.tp_width(n, mid, sd_y, rho, q) <= width) hi <- mid else lo <- mid
        }
      }
      if (n + r * hi < best[3] * (1 - 1e-12)) best <- c(n, hi, n + r * hi)
    }
    n <- n + 1
  }
  best
}

# Sixty fixed planning sets: correlations of either sign, among them ones
# with rho^2 below r / (1 + r), where the lower bound of the search has its
# minimum at the point at which the pool is n + 1; cost ratios from .002 to
# .9; widths from a quarter of a standard deviation to three; three
# confidence levels.
.cost_sets <- function() {
  g <- expand.grid(rho = c(0.15, -0.45, 0.7, -0.9, 0.98),
                   r = c(0.002, 0.05, 0.3, 0.9),
                   width = c(0.25, 0.8, 3))
  g$conf_level <- rep(c(0.90, 0.95, 0.99), length.out = nrow(g))
  g$sd_y <- rep(c(1, 0.02, 40, 7), length.out = nrow(g))
  g$width <- g$width * g$sd_y
  g
}

test_that("the designs reproduce the anchors", {
  a <- .v(ss_aipe_two_phase_cost(true_sd_Y = 1, true_cor = 0.9, width = 0.2,
                                 cost_ratio = 0.05))
  expect_equal(c(a[["necessary_n_labeled"]], a[["necessary_n_pool"]]),
               c(109, 1011))
  expect_equal(a[["necessary_n_labeled"]] + 0.05 * a[["necessary_n_pool"]],
               159.55)
  expect_equal(a[["necessary_n_target_only"]], 387)
  expect_equal(a[["relative_cost"]], 159.55 / 387)
  expect_equal(a[["relative_cost_continuous"]],
               (sqrt(1 - 0.81) + 0.9 * sqrt(0.05))^2)
  expect_equal(a[["threshold_cor"]], 2 * sqrt(0.05) / 1.05)
  # The examples of the help page.
  b <- .v(ss_aipe_two_phase_cost(0.43, 0.55, 0.10, 0.01, conf_level = 0.90))
  expect_equal(c(b[["necessary_n_labeled"]], b[["necessary_n_pool"]],
                 b[["necessary_n_target_only"]]), c(151, 957, 202))
  d <- .v(ss_aipe_two_phase_cost(0.43, 0.55, 0.10, 0.10, conf_level = 0.90))
  expect_equal(d[["threshold_cor"]], 0.57496, tolerance = 1e-4)
  expect_equal(c(d[["necessary_n_labeled"]], d[["necessary_n_pool"]]),
               c(202, 0))
  expect_equal(d[["relative_cost"]], 1)
  expect_equal(d[["relative_cost_continuous"]], 1)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("the design equals an exhaustive search over whole-number designs on 60 planning sets", {
  skip_on_cran()
  g <- .cost_sets()
  expect_identical(nrow(g), 60L)
  expect_gt(sum(g$rho^2 < g$r / (1 + g$r)), 10)
  n_two_phase <- 0
  for (i in seq_len(nrow(g))) {
    info <- sprintf("set %d: rho %.2f, r %.3f, width %.3g, level %.2f", i,
                    g$rho[i], g$r[i], g$width[i], g$conf_level[i])
    plan <- .v(ss_aipe_two_phase_cost(g$sd_y[i], g$rho[i], g$width[i], g$r[i],
                                      conf_level = g$conf_level[i]))
    b <- .tp_brute(g$sd_y[i], g$rho[i], g$width[i], g$r[i], g$conf_level[i])
    expect_equal(c(plan[["necessary_n_labeled"]], plan[["necessary_n_pool"]]),
                 b[1:2], info = info)
    n_two_phase <- n_two_phase + (b[2] > 0)
  }
  expect_gt(n_two_phase, 15)
  expect_lt(n_two_phase, 55)
})

test_that("every returned design meets the width, and the next smaller pool and the next smaller n do not", {
  skip_on_cran()
  g <- .cost_sets()
  for (i in seq_len(nrow(g))) {
    info <- sprintf("set %d", i)
    q <- 1 - (1 - g$conf_level[i]) / 2
    plan <- .v(ss_aipe_two_phase_cost(g$sd_y[i], g$rho[i], g$width[i], g$r[i],
                                      conf_level = g$conf_level[i]))
    n <- plan[["necessary_n_labeled"]]; M <- plan[["necessary_n_pool"]]
    n_t <- plan[["necessary_n_target_only"]]
    expect_true(n == round(n) && M == round(M), info = info)
    expect_lte(.tp_width_target(n_t, g$sd_y[i], q), g$width[i], label = info)
    if (n_t > 2)
      expect_gt(.tp_width_target(n_t - 1, g$sd_y[i], q), g$width[i],
                label = info)
    if (M > 0) {
      expect_true(n >= 3 && M >= n + 1, info = info)
      expect_lte(.tp_width(n, M, g$sd_y[i], g$rho[i], q), g$width[i],
                 label = info)
      if (M > n + 1)
        expect_gt(.tp_width(n, M - 1, g$sd_y[i], g$rho[i], q), g$width[i],
                  label = info)
      if (n > 3)
        expect_gt(.tp_width(n - 1, M, g$sd_y[i], g$rho[i], q), g$width[i],
                  label = info)
      expect_lt(n + g$r[i] * M, n_t, label = info)
      expect_equal(plan[["relative_cost"]], (n + g$r[i] * M) / n_t,
                   info = info)
      expect_equal(plan[["expected_width"]],
                   .tp_width(n, M, g$sd_y[i], g$rho[i], q), info = info)
    } else {
      expect_identical(n, n_t, info = info)
      expect_identical(plan[["relative_cost"]], 1, info = info)
      expect_equal(plan[["expected_width"]],
                   .tp_width_target(n, g$sd_y[i], q), info = info)
    }
    expect_lte(plan[["expected_width"]], g$width[i], label = info)
    # The fractional, normal-quantile figure, and the bound it gives.
    rc <- if (abs(g$rho[i]) > 2 * sqrt(g$r[i]) / (1 + g$r[i]))
      (sqrt(1 - g$rho[i]^2) + abs(g$rho[i]) * sqrt(g$r[i]))^2 else 1
    expect_equal(plan[["relative_cost_continuous"]], rc, info = info)
    s <- (2 * qnorm(q) * g$sd_y[i] / g$width[i])^2
    expect_gte(n + g$r[i] * M, rc * s * (1 - 1e-12), label = info)
  }
})

test_that("ss_aipe_two_phase() returns the labeled units for the returned pool, and the same target-alone size", {
  skip_on_cran()
  g <- .cost_sets()
  checked <- 0
  for (i in seq_len(nrow(g))) {
    plan <- .v(ss_aipe_two_phase_cost(g$sd_y[i], g$rho[i], g$width[i], g$r[i],
                                      conf_level = g$conf_level[i]))
    M <- plan[["necessary_n_pool"]]
    if (M < 4) next
    comp <- .v(ss_aipe_two_phase(true_sd_Y = g$sd_y[i], true_cor = g$rho[i],
                                 width = g$width[i], n_pool = M,
                                 conf_level = g$conf_level[i]))
    info <- sprintf("set %d", i)
    expect_equal(comp[["necessary_n_labeled"]], plan[["necessary_n_labeled"]],
                 info = info)
    expect_equal(comp[["necessary_n_target_only"]],
                 plan[["necessary_n_target_only"]], info = info)
    expect_equal(comp[["expected_width"]], plan[["expected_width"]],
                 info = info)
    checked <- checked + 1
  }
  expect_gt(checked, 15)
})

test_that("a weak proxy, a proxy with no correlation, and a wide width return the target alone", {
  a <- .v(ss_aipe_two_phase_cost(1, 0.1, 0.2, 0.05))
  expect_equal(a[["necessary_n_pool"]], 0)
  expect_equal(a[["necessary_n_labeled"]], a[["necessary_n_target_only"]])
  z <- .v(ss_aipe_two_phase_cost(1, 0, 0.2, 0.05))
  expect_equal(c(z[["necessary_n_labeled"]], z[["necessary_n_pool"]]),
               c(387, 0))
  # The target alone needs at least 2 units, and at 2 units the t quantile
  # on 1 degree of freedom is part of the requirement.
  w <- .v(ss_aipe_two_phase_cost(1, 0.5, 100, 0.5))
  expect_equal(c(w[["necessary_n_labeled"]], w[["necessary_n_pool"]]),
               c(2, 0))
  w3 <- .v(ss_aipe_two_phase_cost(1, 0.99, 10, 0.01))
  expect_equal(c(w3[["necessary_n_labeled"]], w3[["necessary_n_pool"]]),
               c(3, 0))
  expect_equal(w3[["relative_cost_continuous"]], 0.0576, tolerance = 1e-3)
  # The sign of the correlation does not matter.
  expect_equal(.v(ss_aipe_two_phase_cost(1, -0.9, 0.2, 0.05))[1:7],
               .v(ss_aipe_two_phase_cost(1, 0.9, 0.2, 0.05))[1:7])
})

test_that("the target alone can be the design above the threshold correlation", {
  # The threshold is for fractional sizes under the normal quantile. Here
  # the correlation is above it, and the least costly whole-number design
  # is still the target alone.
  r <- 0.3
  thr <- 2 * sqrt(r) / (1 + r)
  rho <- thr + 0.02
  plan <- .v(ss_aipe_two_phase_cost(1, rho, 0.8, r))
  expect_lt(plan[["relative_cost_continuous"]], 1)
  expect_equal(plan[["necessary_n_pool"]], 0)
  b <- .tp_brute(1, rho, 0.8, r, 0.95)
  expect_equal(b[2], 0)
})

test_that("the cost planner joins the AIPE tidiers and refuses a costly proxy", {
  skip_on_cran()
  res <- ss_aipe_two_phase_cost(0.43, 0.55, 0.10, 0.01, conf_level = 0.90)
  expect_s3_class(res, "dmar_ss_aipe")
  expect_s3_class(res, "dmar_tbl")
  expect_equal(attr(res, "conf_level"), 0.90)
  expect_equal(generics::tidy(res)$estimate,
               res$value[res$term == "necessary_n_labeled"])
  expect_equal(generics::glance(res)$necessary_n_pool,
               res$value[res$term == "necessary_n_pool"])
  expect_error(ss_aipe_two_phase_cost(1, 0.5, 0.2, cost_ratio = 1),
               "below 1")
  expect_error(ss_aipe_two_phase_cost(1, 0.5, 0.2, cost_ratio = 0),
               "positive")
  # Both refusals are classed and carry the value given.
  for (bad in c(1, 2.5, 0, -0.1)) {
    err <- expect_error(ss_aipe_two_phase_cost(1, 0.5, 0.2, cost_ratio = bad),
                        class = "dmar_cost_ratio_invalid")
    expect_s3_class(err, "dmar_refusal")
    expect_identical(err$value, bad)
  }
  expect_error(ss_aipe_two_phase_cost(1, 0.5, 0.2, cost_ratio = c(0.1, 0.2)),
               class = "dmar_cost_ratio_invalid")
  expect_error(ss_aipe_two_phase_cost(1, 1.2, 0.2, cost_ratio = 0.1), "'true_cor'")
})

test_that("the table has its twelve rows in the documented order", {
  res <- ss_aipe_two_phase_cost(1, 0.8, 0.2, 0.01)
  expect_identical(res$term,
                   c("necessary_n_labeled", "necessary_n_pool",
                     "relative_cost", "relative_cost_continuous",
                     "expected_width", "necessary_n_target_only",
                     "threshold_cor", "width", "true_sd_Y", "true_cor", "cost_ratio",
                     "conf_level"))
  expect_true(all(c("relative_cost", "relative_cost_continuous",
                    "expected_width", "necessary_n_target_only") %in%
                    names(generics::glance(res))))
})

test_that("sizes that cannot be represented are refused, and the scale of sd_y does not matter", {
  expect_error(ss_aipe_two_phase_cost(1e200, 0.5, 1, 0.1),
               class = "dmar_two_phase_no_design")
  expect_error(ss_aipe_two_phase_cost(1, 0.9, 1, 1e-320),
               class = "dmar_two_phase_no_design")
  expect_error(ss_aipe_two_phase_cost(1, 0.5, 1e-8, 0.1),
               class = "dmar_two_phase_no_design")
  tiny <- .v(ss_aipe_two_phase_cost(1e-200, 0.9, 1, 0.1))
  expect_equal(c(tiny[["necessary_n_labeled"]], tiny[["necessary_n_pool"]]),
               c(2, 0))
  # The design depends on sd_y and width through their ratio alone.
  p1 <- .v(ss_aipe_two_phase_cost(1, 0.8, 0.25, 0.01))
  p2 <- .v(ss_aipe_two_phase_cost(1e-150, 0.8, 0.25e-150, 0.01))
  expect_equal(p2[c(1:4, 6)], p1[c(1:4, 6)])
  expect_equal(p2[["expected_width"]] * 1e150, p1[["expected_width"]])
  # A very large design is still found: the target alone at a size where
  # the t quantile and the normal quantile agree to the last digit.
  wide <- .v(ss_aipe_two_phase_cost(1, 0.5, 1e-7, 0.1))
  expect_lte(abs(wide[["necessary_n_labeled"]] - 1536583528277650), 3)
  expect_equal(wide[["necessary_n_pool"]], 0)
})

test_that("a large design is found quickly", {
  # About a million labeled units; the search examines the n whose lower
  # bound is below the best cost, a few thousand of them.
  t0 <- proc.time()[["elapsed"]]
  big <- .v(ss_aipe_two_phase_cost(1, 0.9, 0.002, 0.05))
  expect_lt(proc.time()[["elapsed"]] - t0, 5)
  n <- big[["necessary_n_labeled"]]; M <- big[["necessary_n_pool"]]
  expect_lte(.tp_width(n, M, 1, 0.9, 0.975), 0.002)
  expect_gt(.tp_width(n, M - 1, 1, 0.9, 0.975), 0.002)
  # Whole-number cost within a few units of the fractional bound.
  s <- (2 * qnorm(0.975) / 0.002)^2
  expect_lt(n + 0.05 * M - big[["relative_cost_continuous"]] * s, 10)
})

test_that("the design equals an exhaustive search on 600 random planning sets", {
  skip_on_cran()
  set.seed(31207)
  K <- 600
  same <- 0
  two_phase <- 0
  for (k in seq_len(K)) {
    rho <- runif(1, 0.02, 0.999) * sample(c(-1, 1), 1)
    r <- 10^runif(1, -3, log10(0.98))
    sd_y <- 10^runif(1, -2, 2)
    width <- 10^runif(1, -0.9, 0.6) * sd_y
    cl <- sample(c(0.80, 0.90, 0.95, 0.99), 1)
    plan <- .v(ss_aipe_two_phase_cost(sd_y, rho, width, r, conf_level = cl))
    b <- .tp_brute(sd_y, rho, width, r, cl)
    same <- same + (plan[["necessary_n_labeled"]] == b[1] &&
                      plan[["necessary_n_pool"]] == b[2])
    two_phase <- two_phase + (b[2] > 0)
  }
  expect_equal(same, K)
  expect_gt(two_phase, 150)
})

test_that("the sensitivity function accepts a returned design, and its median width is near the plan", {
  skip_on_cran()
  plan <- .v(ss_aipe_two_phase_cost(1, 0.9, 0.5, 0.05))
  res <- .v(ss_aipe_two_phase_sensitivity(
    true_sd_Y = 1, true_cor = 0.9,
    specified_n_labeled = plan[["necessary_n_labeled"]],
    n_pool = plan[["necessary_n_pool"]], width = 0.5, G = 400, seed = 7))
  expect_true(is.finite(res[["mean_ci_width"]]))
  expect_lt(abs(res[["median_ci_width"]] - plan[["expected_width"]]),
            0.1 * plan[["expected_width"]])
})

test_that("the planning checks of the two-phase family name the argument they refuse", {
  # .two_phase_check_planning() and .two_phase_prepare() built the message
  # from the argument's name but stored arg = NULL; the condition now
  # carries the name.
  arg_of <- function(expr) tryCatch(expr, error = function(e) e)$arg
  expect_identical(arg_of(ss_aipe_two_phase_cost(true_sd_Y = -1, true_cor = .5,
                                                 width = .2, cost_ratio = .1)),
                   "true_sd_Y")
  expect_identical(arg_of(ss_aipe_two_phase_cost(true_sd_Y = 1, true_cor = 1.2,
                                                 width = .2, cost_ratio = .1)),
                   "true_cor")
  expect_identical(arg_of(ss_aipe_two_phase(true_sd_Y = 1, true_cor = c(.3, .4),
                                            width = .2)),
                   "true_cor")
  expect_identical(arg_of(ss_aipe_two_phase_sensitivity(
    true_sd_Y = 1, true_cor = .5, supposed_sd_Y = 0, supposed_cor = .5,
    width = .3, n_pool = 500, G = 2)), "supposed_sd_Y")
  expect_identical(arg_of(ss_aipe_two_phase_sensitivity(
    true_sd_Y = 1, true_cor = .5, supposed_sd_Y = 1, supposed_cor = -1,
    width = .3, n_pool = 500, G = 2)), "supposed_cor")
  expect_identical(arg_of(ci_two_phase(y = letters[1:5], x = 1:5,
                                       x_unlabeled = 1:10)), "y")
  expect_identical(arg_of(ci_two_phase(y = 1:5, x = c(1:4, Inf),
                                       x_unlabeled = 1:10)), "x")
})
