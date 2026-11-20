## Tests for ss_power_mixed_effects(), the two-level (random-intercept) power
## planner for a treatment effect at the higher level.

test_that("ss_power_mixed_effects() the necessary number of clusters per arm grows with ICC", {
  low_icc  <- ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.01, desired_power = 0.80)
  high_icc <- ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.20, desired_power = 0.80)
  expect_gt(high_icc[high_icc$term == "necessary_n_clusters_per_arm", "value"],
            low_icc[low_icc$term   == "necessary_n_clusters_per_arm", "value"])
})

test_that("ss_power_mixed_effects() the necessary number of clusters per arm shrinks with larger d", {
  small <- ss_power_mixed_effects(true_smd = 0.20, cluster_size = 30, true_icc = 0.05, desired_power = 0.80)
  large <- ss_power_mixed_effects(true_smd = 0.50, cluster_size = 30, true_icc = 0.05, desired_power = 0.80)
  expect_gt(small[small$term == "necessary_n_clusters_per_arm", "value"],
            large[large$term == "necessary_n_clusters_per_arm", "value"])
})

test_that("ss_power_mixed_effects() with true_icc = 0 reduces toward an unclustered two-sample t-test", {
  # With true_icc = 0, the design effect 1 + (n-1)*rho = 1, so each cluster contributes n
  # independent observations. The total sample per arm is J*n, and the test reduces to
  # a standard two-sample t-test on J*n per arm.
  cluster_res  <- ss_power_mixed_effects(true_smd = 0.30, cluster_size = 5, true_icc = 0, desired_power = 0.80)
  J_required <- cluster_res[cluster_res$term == "necessary_n_clusters_per_arm", "value"]
  total_n_per_arm <- J_required * 5
  smd_res <- ss_power_smd(true_smd = 0.30, n_1 = total_n_per_arm)
  realized_power <- smd_res[smd_res$term == "actual_power", "value"]
  # The cluster-level test uses df = 2J - 2 instead of 2*Jn - 2, so it's a bit more conservative.
  # Realized power from the smd unclustered test should at least equal the desired_power level.
  expect_gte(realized_power, 0.80)
})

test_that("ss_power_mixed_effects() power at the returned number of clusters equals or exceeds desired_power", {
  res <- ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.05, desired_power = 0.80)
  J   <- res[res$term == "necessary_n_clusters_per_arm", "value"]
  pwr <- ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.05, specified_n_clusters_per_arm = J)
  expect_gte(pwr[pwr$term == "actual_power", "value"], 0.80)
})

test_that("ss_power_mixed_effects() validates inputs", {
  expect_error(ss_power_mixed_effects(true_smd = 0,    cluster_size = 30, true_icc = 0.05, desired_power = 0.80),
               "'true_smd' must be a single finite nonzero numeric")
  expect_error(ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 1.0,  desired_power = 0.80),
               "'true_icc' must be a single number in", class = "dmar_argument_invalid")
  expect_error(ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = -0.1, desired_power = 0.80),
               "'true_icc' must be a single number in", class = "dmar_argument_invalid")
})

test_that("ss_power_mixed_effects() names its sizes, its planning value, and its rows as the package does", {
  expect_identical(names(formals(ss_power_mixed_effects)),
                   c("true_smd", "cluster_size", "true_icc", "specified_n_clusters_per_arm",
                     "desired_power", "alpha_level", "directional"))
  planned <- ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.05,
                                    desired_power = 0.80)
  expect_identical(planned$term, c("necessary_n_clusters_per_arm", "total_N",
                                   "noncentral_t_parm", "actual_power"))
  J <- planned$value[1]
  expect_identical(planned$value[2], 2 * J * 30)
  expect_equal(planned$value[3], 0.30 * sqrt(J * 30 / (2 * (1 + 29 * 0.05))))
  given <- ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.05,
                                  specified_n_clusters_per_arm = 25)
  expect_identical(given$term, c("specified_n_clusters_per_arm", "total_N",
                                 "noncentral_t_parm", "actual_power"))
  expect_identical(generics::tidy(planned)$estimate, J)
  expect_identical(generics::tidy(given)$estimate, 25)
})

test_that("ss_power_mixed_effects() refuses the 1.0.x names and a missing or inadmissible size by class", {
  # A complete 1.0.x call names rho (and perhaps J), which R no longer finds.
  # (suppressWarnings: under warnPartialMatchArgs R warns that `n` matches
  # n_clusters_per_arm before it stops for `rho`.)
  expect_error(suppressWarnings(
    do.call(ss_power_mixed_effects,
            list(true_smd = 0.30, n = 30, rho = 0.05, desired_power = 0.80))),
    "unused argument")
  expect_error(do.call(ss_power_mixed_effects,
                       list(true_smd = 0.30, cluster_size = 30, true_icc = 0.05, J = 25)),
               "unused argument")
  # `n` matches no argument now that the size is
  # specified_n_clusters_per_arm: R stops the call.
  expect_error(do.call(ss_power_mixed_effects, list(true_smd = 0.30, n = 30)),
               "unused argument")
  # A call that leaves cluster_size out stops and names the new arguments.
  err <- tryCatch(ss_power_mixed_effects(true_smd = 0.30, true_icc = 0.05),
                  error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "cluster_size")
  expect_match(conditionMessage(err), "'cluster_size' is missing", fixed = TRUE)
  expect_match(conditionMessage(err), "n_clusters_per_arm", fixed = TRUE)
  for (bad in list(1, 2.5, c(10, 20), "30", NA_real_, Inf)) {
    err <- tryCatch(ss_power_mixed_effects(true_smd = 0.30, cluster_size = bad, true_icc = 0.05),
                    error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "cluster_size")
  }
  for (bad in list(1, 2.5, c(10, 20), "30")) {
    err <- tryCatch(ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30, true_icc = 0.05,
                                           specified_n_clusters_per_arm = bad),
                    error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "specified_n_clusters_per_arm")
  }
  err <- tryCatch(ss_power_mixed_effects(true_smd = 0.30, cluster_size = 30),
                  error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "true_icc")
})
