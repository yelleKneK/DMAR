
test_that("ss_aipe_equivalence_smd_sensitivity maps alpha to (1 - conf_level)/2 (HIGH-03)", {
  skip_on_cran()
  # The planner ss_aipe_equivalence_smd() sizes the width of a (1 - 2*alpha)
  # interval, and the simulation scores the interval of equivalence_smd() at
  # conf_level (ci_smd() through DMAR 1.0.x), so the resolved sample size must
  # equal the planner called with alpha = (1 - conf_level)/2. The bug used
  # alpha = 1 - conf_level, planning a 90% width while scoring a 95% one.
  for (cl in c(0.90, 0.95, 0.99)) {
    oracle <- ss_aipe_equivalence_smd(true_smd = 0, width = 0.20,
                               alpha_level = (1 - cl) / 2, balanced = TRUE)
    oracle_n <- oracle$value[oracle$term == "necessary_n_per_group"]
    s <- ss_aipe_equivalence_smd_sensitivity(true_smd = 0, supposed_smd = 0,
                                      width = 0.20, smd_upper = 0.20,
                                      conf_level = cl, G = 1)
    expect_equal(s$value[s$term == "n_per_group"], oracle_n)
  }
})

# ---- The simulation scores the interval equivalence_smd() reports (DR17) ----

test_that("each replication's interval is the one equivalence_smd() reports for the same scores", {
  skip_on_cran()
  # The planner draws nothing, so the call's stream is the seed followed by
  # rnorm(n, true_smd) and rnorm(n) for each replication; replaying it here
  # and analyzing the scores with equivalence_smd(x, y) must give the
  # recorded estimate and limits (to the 15 digits the file keeps).
  for (cs in list(list(true = 0.3, est = 0.3, cl = 0.90, seed = 12),
                  list(true = -0.5, est = 0.2, cl = 0.95, seed = 13))) {
    f <- tempfile(fileext = ".csv")
    s <- ss_aipe_equivalence_smd_sensitivity(
      true_smd = cs$true, supposed_smd = cs$est, width = 0.8,
      smd_upper = 0.8, conf_level = cs$cl, G = 15, seed = cs$seed,
      filename = f)
    reps <- utils::read.csv(f)
    unlink(f)
    n <- s$value[s$term == "n_per_group"]
    plan <- ss_aipe_equivalence_smd(true_smd = cs$est, width = 0.8,
                                    alpha_level = (1 - cs$cl) / 2)
    expect_identical(n, plan$value[plan$term == "necessary_n_per_group"])
    set.seed(cs$seed)
    for (g in seq_len(15)) {
      x <- stats::rnorm(n, mean = cs$true)
      y <- stats::rnorm(n)
      e <- equivalence_smd(x = x, y = y, smd_upper = 50,
                           alpha_level = (1 - cs$cl) / 2)
      expect_equal(reps$smd[g], e$value[e$term == "smd"], tolerance = 1e-13)
      expect_equal(reps$ci_lower[g], e$value[e$term == "ci_lower"],
                   tolerance = 1e-13)
      expect_equal(reps$ci_upper[g], e$value[e$term == "ci_upper"],
                   tolerance = 1e-13)
    }
  }
})

test_that("pct_ci_less_width estimates the exact probability for the reported interval", {
  skip_on_cran()
  # Exact probability from R's pt() (every noncentrality here is below
  # 37.62): the interval is no wider than the target exactly when
  # d^2 <= n w^2 / t^2 - 8, and d sqrt(n/2) is noncentral t.
  p_exact <- function(n, d, w, a) {
    k <- sqrt((n * w^2 / stats::qt(1 - a, 2 * n - 2)^2 - 8) * n / 2)
    stats::pt(k, 2 * n - 2, abs(d) * sqrt(n / 2)) -
      stats::pt(-k, 2 * n - 2, abs(d) * sqrt(n / 2))
  }
  G <- 400L
  s <- ss_aipe_equivalence_smd_sensitivity(
    true_smd = 0.5, supposed_smd = 0.5, width = 0.6, smd_upper = 1,
    conf_level = 0.90, assurance = 0.80, G = G, seed = 7)
  n <- s$value[s$term == "n_per_group"]
  ex <- p_exact(n, 0.5, 0.6, 0.05)
  plan <- ss_aipe_equivalence_smd(true_smd = 0.5, width = 0.6,
                                  alpha_level = 0.05, assurance = 0.80)
  expect_equal(plan$value[plan$term == "assurance_attained"], ex,
               tolerance = 1e-9)
  sim <- attr(s, "simulation")
  expect_identical(c(sim$attempted, sim$returned, sim$failed), c(G, G, 0L))
  expect_lt(abs(s$value[s$term == "pct_ci_less_width"] - ex),
            4 * sqrt(ex * (1 - ex) / G))
})

# ---- Arguments that are not finite (review of the DR17 change) ----
# With the interval formed by .equivalence_smd_interval() in place of
# ci_smd(), a missing or infinite true_smd no longer stopped inside the
# loop: the call returned a table of NaN (with only R's "NAs produced"
# warnings for NA, and nothing for Inf). The arguments are now checked on
# entry, as ss_aipe_smd_sensitivity() checks true_smd.
test_that("a true_smd that is not a single finite number is refused by name", {
  for (v in list(NA_real_, NaN, Inf, -Inf, NA, c(0, 1))) {
    e <- expect_error(ss_aipe_equivalence_smd_sensitivity(
      true_smd = v, specified_n_per_group = 10, width = 1,
      smd_upper = 0.5, G = 5, seed = 1),
      class = "dmar_argument_invalid")
    expect_identical(e$arg, "true_smd")
  }
})

test_that("supposed_smd, specified_n_per_group, the bounds, and width are checked on entry", {
  refused_arg <- function(...) {
    e <- expect_error(ss_aipe_equivalence_smd_sensitivity(..., G = 5, seed = 1),
                      class = "dmar_argument_invalid")
    e$arg
  }
  for (v in list(NA_real_, Inf, -Inf))
    expect_identical(refused_arg(supposed_smd = v, width = 1,
                                 smd_upper = 0.5), "supposed_smd")
  for (v in list(NA_real_, Inf, 1, 10.7, c(10, 12)))
    expect_identical(refused_arg(specified_n_per_group = v, width = 1,
                                 smd_upper = 0.5), "specified_n_per_group")
  expect_identical(refused_arg(specified_n_per_group = 10, width = 1,
                               smd_upper = NA_real_), "smd_upper")
  expect_identical(refused_arg(specified_n_per_group = 10, width = 1,
                               smd_upper = c(0.5, 0.6)), "smd_upper")
  expect_identical(refused_arg(specified_n_per_group = 10, width = 1,
                               smd_lower = NA_real_, smd_upper = 0.5),
                   "smd_lower")
  expect_identical(refused_arg(specified_n_per_group = 10,
                               smd_upper = 0.5), "width")
  expect_identical(refused_arg(supposed_smd = 0, smd_upper = 0.5),
                   "width")
  # The smallest size equivalence_smd() accepts is still accepted.
  s <- ss_aipe_equivalence_smd_sensitivity(specified_n_per_group = 2,
                                           width = 10, smd_upper = 0.5,
                                           G = 5, seed = 1)
  expect_identical(s$value[s$term == "n_per_group"], 2)
})

test_that("a per-group size above 46,340 gives finite intervals, as equivalence_smd() does", {
  # n_1 * n_2 overflowed R's integer range above 46,340 per group when the
  # sizes arrived as integers, and every simulated interval was NA.
  out <- ss_aipe_equivalence_smd_sensitivity(
    true_smd = 0, specified_n_per_group = 46341L, width = 0.02,
    conf_level = 0.90, smd_upper = 1, G = 3, seed = 1)
  v <- function(term) out$value[out$term == term]
  expect_true(is.finite(v("mean_ci_width")))
  expect_identical(v("pct_equivalent"), 1)
  iv <- DMAR:::.equivalence_smd_interval(0.1, 46341L, 46341L, 0.05)
  n <- 46341
  se <- sqrt(2 / n + 0.1^2 / (4 * n))
  expect_equal(iv$se, se, tolerance = 1e-14)
  expect_equal(iv$width, 2 * stats::qt(0.95, 2 * n - 2) * se,
               tolerance = 1e-12)
})

test_that("ss_aipe_equivalence_smd_sensitivity() gives no warning", {
  expect_no_warning(ss_aipe_equivalence_smd_sensitivity(
    true_smd = 0, specified_n_per_group = 20, width = 0.8, smd_upper = 0.5,
    G = 2, seed = 1))
})
