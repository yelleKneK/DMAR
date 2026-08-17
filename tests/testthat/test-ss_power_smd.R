test_that("ss_power_smd() with smd alone solves for the per-group n", {
  result <- ss_power_smd(smd = 0.5)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_true("necessary_n_per_group" %in% result$term)
  expect_true("actual_power"          %in% result$term)
  expect_gt(result$value[result$term == "necessary_n_per_group"], 0)
  expect_gte(result$value[result$term == "actual_power"], 0.85)
})

test_that("ss_power_smd() with both smd and n_1 returns realized power", {
  result <- ss_power_smd(smd = 0.5, n_1 = 30)
  expect_true("specified_n_1" %in% result$term)
  expect_true("actual_power"  %in% result$term)
  expect_equal(result$value[result$term == "specified_n_1"], 30)
})

test_that("ss_power_smd() requires fewer subjects for a larger effect", {
  small <- ss_power_smd(smd = 0.3)$value[1]
  large <- ss_power_smd(smd = 0.8)$value[1]
  expect_gt(small, large)
})

test_that("ss_power_smd() echoes the planning inputs as numeric rows", {
  r <- ss_power_smd(smd = 0.5, desired_power = 0.80, alpha_level = 0.05)
  expect_equal(r$term[1], "necessary_n_per_group")      # result stays first
  # the supposed effect is denoted as such, not as a bare "smd"
  expect_equal(r$value[r$term == "supposed_smd"], 0.5)
  expect_false("smd" %in% r$term)
  expect_equal(r$value[r$term == "desired_power"], 0.80)
  expect_equal(r$value[r$term == "alpha_level"], 0.05)
  expect_equal(r$value[r$term == "tails"], 2)           # nondirectional
  expect_type(r$value, "double")                        # value column numeric
  d <- ss_power_smd(smd = 0.5, desired_power = 0.80, directional = TRUE)
  expect_equal(d$value[d$term == "tails"], 1)           # directional
  # the specified-n branch echoes supposed_smd, alpha_level, and tails
  # alongside the supplied group sizes
  s <- ss_power_smd(smd = 0.5, n_1 = 30)
  expect_equal(s$value[s$term == "supposed_smd"], 0.5)
  expect_equal(s$value[s$term == "tails"], 2)
})

test_that("ss_power_smd() rejects a null effect in the sample-size search", {
  # smd = 0 has power equal to alpha at every N, so the search cannot converge;
  # it must stop cleanly and quickly, not iterate to the sample-size cap.
  expect_error(ss_power_smd(smd = 0, desired_power = 0.85), "nonzero")
  # power at a *specified* N with a null effect is a valid query (= alpha).
  sp <- ss_power_smd(smd = 0, n_1 = 30)
  expect_equal(sp$value[sp$term == "actual_power"], 0.05, tolerance = 1e-9)
})

test_that("ss_power_smd rejects fractional group sizes (MEDIUM-03)", {
  expect_error(ss_power_smd(smd = 0.5, n_1 = 20.5, n_2 = 21.25), "whole number")
  expect_error(ss_power_smd(smd = 0.5, n_1 = 1), "whole number")
})
