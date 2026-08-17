test_that("ss_aipe_mixed_effects() returns documented rows", {
  res <- ss_aipe_mixed_effects(sigma2_y = 1, sigma2_x = 1, icc = 0.10,
                                    width = 0.20, cluster_size = 20)
  expect_setequal(res$term,
                  c("necessary_n_clusters", "total_N", "width", "icc",
                    "cluster_size", "ci_width_expected"))
})

test_that("ss_aipe_mixed_effects() expected width matches target", {
  res <- ss_aipe_mixed_effects(sigma2_y = 1, sigma2_x = 1, icc = 0.10,
                                    width = 0.20)
  expect_lte(res$value[res$term == "ci_width_expected"], 0.20 + 1e-8)
})

test_that("ss_aipe_mixed_effects() smaller width => more clusters", {
  r1 <- ss_aipe_mixed_effects(1, 1, icc = 0.10, width = 0.30)
  r2 <- ss_aipe_mixed_effects(1, 1, icc = 0.10, width = 0.10)
  expect_lt(r1$value[1], r2$value[1])
})

test_that("ss_aipe_mixed_effects() rejects bad inputs", {
  expect_error(ss_aipe_mixed_effects(1, 1, icc = 1.2, width = 0.2),
               "in \\[0, 1\\)")
  expect_error(ss_aipe_mixed_effects(-1, 1, icc = 0.1, width = 0.2),
               "positive")
})

test_that("ss_aipe_mixed_effects() fires a type guard for every formal", {
  # Non-numeric or non-scalar values are rejected by name; no argument
  # is silently accepted.
  expect_error(ss_aipe_mixed_effects("1", 1, icc = 0.1, width = 0.2),
               "'sigma2_y' must be a single numeric value")
  expect_error(ss_aipe_mixed_effects(1, c(1, 2), icc = 0.1, width = 0.2),
               "'sigma2_x' must be a single numeric value")
  expect_error(ss_aipe_mixed_effects(1, 1, icc = "0.1", width = 0.2),
               "'icc' must be a single numeric value")
  expect_error(ss_aipe_mixed_effects(1, 1, icc = 0.1, width = TRUE),
               "'width' must be a single numeric value")
  expect_error(ss_aipe_mixed_effects(1, 1, icc = 0.1, width = 0.2,
                                     cluster_size = "20"),
               "'cluster_size' must be a single numeric value")
  expect_error(ss_aipe_mixed_effects(1, 1, icc = 0.1, width = 0.2,
                                     conf_level = "0.95"),
               "'conf_level' must be a single numeric value")
})

test_that("ss_aipe_mixed_effects() fires a range guard for every formal", {
  expect_error(ss_aipe_mixed_effects(-1, 1, icc = 0.1, width = 0.2),
               "'sigma2_y' must be positive")
  expect_error(ss_aipe_mixed_effects(1, -1, icc = 0.1, width = 0.2),
               "'sigma2_x' must be positive")
  expect_error(ss_aipe_mixed_effects(1, 1, icc = 1.2, width = 0.2),
               "'icc' must be in \\[0, 1\\)")
  expect_error(ss_aipe_mixed_effects(1, 1, icc = 0.1, width = -0.2),
               "'width' must be positive")
  expect_error(ss_aipe_mixed_effects(1, 1, icc = 0.1, width = 0.2,
                                     cluster_size = 1),
               "'cluster_size' must be at least 2")
  expect_error(ss_aipe_mixed_effects(1, 1, icc = 0.1, width = 0.2,
                                     conf_level = 1.5),
               "'conf_level' must be in \\(0, 1\\)")
})

test_that("ss_aipe_mixed_effects() no longer accepts a 'beta' argument", {
  # The closed-form variance of the slope does not involve the slope
  # value, so the never-read formal was removed; a supplied beta fails
  # loudly instead of being silently ignored.
  expect_error(ss_aipe_mixed_effects(sigma2_y = 1, sigma2_x = 1, icc = 0.1,
                                     width = 0.2, beta = 0.3),
               "unused argument")
})

test_that("expected half width follows the documented Var(beta) with no design-effect factor", {
  # Var(beta-hat) = sigma2_y (1 - rho_I) / (N sigma2_x) for the
  # cluster-mean-centered level-1 slope. The trailing design-effect
  # factor the page previously printed would give a half width of
  # 0.0373 on these inputs, not the reported 0.0980, so this
  # recomputation discriminates the two formulas.
  res <- ss_aipe_mixed_effects(sigma2_y = 1, sigma2_x = 1, icc = 0.10,
                               width = 0.20, cluster_size = 20)
  v <- function(t) res$value[res$term == t]
  N <- v("necessary_n_clusters") * 20
  hw_doc <- stats::qnorm(0.975) * sqrt(1 * (1 - 0.10) / (N * 1))
  expect_equal(v("ci_width_expected"), 2 * hw_doc, tolerance = 1e-12)
  # And the planner is minimal: one fewer cluster misses the target.
  N_less <- (v("necessary_n_clusters") - 1) * 20
  expect_gt(stats::qnorm(0.975) * sqrt(0.9 / N_less), 0.10)
})
