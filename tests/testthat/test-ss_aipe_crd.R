## Sample-size planning for cluster-randomized designs (Lai & Kelley, 2011/2018)
## --------------------------------------------------------------------------
## The non-effect size CRD family targets the population mean DIFFERENCE in
## treatment minus control. Six exported functions handle the three "what
## you control" axes (n_clusters, n_individuals, both) crossed with two "what you
## constrain" axes (fixed_width, fixed_budget).

test_that("ss_aipe_crd_n_clusters_fixed_width() returns a tidy data.frame with n_clusters + expected width", {
  res <- ss_aipe_crd_n_clusters_fixed_width(width = 0.3, n_individuals = 20, pr_treat = 0.5,
                                      icc_Y = 0.1, total_var = 1)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_n_clusters" %in% res$term)
  expect_true(any(grepl("exp_width", res$term)))
  expect_true(res$value[res$term == "necessary_n_clusters"] > 0)
})

test_that("ss_aipe_crd_n_clusters_fixed_width() achieves expected width <= target width", {
  width_target <- 0.30
  res <- ss_aipe_crd_n_clusters_fixed_width(width = width_target, n_individuals = 20,
                                      pr_treat = 0.5, icc_Y = 0.1, total_var = 1)
  expect_lte(res$value[grepl("exp_width", res$term)], width_target + 1e-6)
})

test_that("ss_aipe_crd_n_clusters_fixed_width() necessary n_clusters grows with higher ICC", {
  low_icc  <- ss_aipe_crd_n_clusters_fixed_width(width = 0.3, n_individuals = 20,
                                           pr_treat = 0.5, icc_Y = 0.05,
                                           total_var = 1)$value[1]
  high_icc <- ss_aipe_crd_n_clusters_fixed_width(width = 0.3, n_individuals = 20,
                                           pr_treat = 0.5, icc_Y = 0.30,
                                           total_var = 1)$value[1]
  expect_gt(high_icc, low_icc)
})

test_that("ss_aipe_crd_n_clusters_fixed_width() smaller width target requires more clusters", {
  wide <- ss_aipe_crd_n_clusters_fixed_width(width = 0.50, n_individuals = 20,
                                       pr_treat = 0.5, icc_Y = 0.1,
                                       total_var = 1)$value[1]
  tight <- ss_aipe_crd_n_clusters_fixed_width(width = 0.20, n_individuals = 20,
                                        pr_treat = 0.5, icc_Y = 0.1,
                                        total_var = 1)$value[1]
  expect_gt(tight, wide)
})

test_that("ss_aipe_crd_n_individuals_fixed_width() errors when n_clusters is too small to reach the target width", {
  expect_error(
    ss_aipe_crd_n_individuals_fixed_width(width = 0.3, n_clusters = 30, pr_treat = 0.5,
                                  icc_Y = 0.1, total_var = 1),
    "impossible to achieve the target width"
  )
})

test_that("ss_aipe_crd_n_individuals_fixed_width() returns a tidy data.frame when n_clusters is large enough", {
  res <- ss_aipe_crd_n_individuals_fixed_width(width = 0.30, n_clusters = 200,
                                       pr_treat = 0.5, icc_Y = 0.1,
                                       total_var = 1)
  expect_s3_class(res, "data.frame")
  expect_true("cluster_size" %in% res$term)
  expect_lte(res$value[grepl("exp_width", res$term)], 0.30 + 1e-6)
})

test_that("ss_aipe_crd_both_fixed_width() jointly returns n_clusters, cluster_size, and a budget total", {
  res <- ss_aipe_crd_both_fixed_width(width = 0.30, clus_cost = 10,
                                     indiv_cost = 1, pr_treat = 0.5,
                                     icc_Y = 0.1, total_var = 1)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("necessary_n_clusters", "cluster_size", "budget") %in% res$term))
  n_clusters  <- res$value[res$term == "necessary_n_clusters"]
  csize  <- res$value[res$term == "cluster_size"]
  budget <- res$value[res$term == "budget"]
  expect_equal(budget, n_clusters * 10 + n_clusters * csize * 1, tolerance = 1)
})

test_that("ss_aipe_crd_n_clusters_fixed_budget() respects the budget constraint", {
  res <- ss_aipe_crd_n_clusters_fixed_budget(budget = 1000, n_individuals = 20,
                                       clus_cost = 10, indiv_cost = 1,
                                       pr_treat = 0.5, icc_Y = 0.1, total_var = 1)
  expect_s3_class(res, "data.frame")
  expect_lte(res$value[res$term == "budget"], 1000)
  expect_true(res$value[res$term == "necessary_n_clusters"] > 0)
})

test_that("ss_aipe_crd_n_individuals_fixed_budget() respects the budget constraint", {
  res <- ss_aipe_crd_n_individuals_fixed_budget(budget = 1000, n_clusters = 30,
                                        clus_cost = 10, indiv_cost = 1,
                                        pr_treat = 0.5, icc_Y = 0.1, total_var = 1)
  expect_s3_class(res, "data.frame")
  expect_lte(res$value[res$term == "budget"], 1000)
  expect_true(res$value[res$term == "cluster_size"] > 0)
})

test_that("ss_aipe_crd_both_fixed_budget() returns a feasible joint (n_clusters, cluster_size) under budget", {
  res <- ss_aipe_crd_both_fixed_budget(budget = 1000, clus_cost = 10,
                                      indiv_cost = 1, pr_treat = 0.5,
                                      icc_Y = 0.1, total_var = 1)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("necessary_n_clusters", "cluster_size", "budget") %in% res$term))
  expect_lte(res$value[res$term == "budget"], 1000)
})

test_that("ss_aipe_crd_*_fixed_budget() expected widths shrink as the budget grows", {
  small <- ss_aipe_crd_both_fixed_budget(budget =  500, clus_cost = 10,
                                        indiv_cost = 1, pr_treat = 0.5,
                                        icc_Y = 0.1,
                                        total_var = 1)$value[grepl("exp_width", c("necessary_n_clusters","cluster_size","exp_width_of_unstd_conditions_diff","budget"))]
  big   <- ss_aipe_crd_both_fixed_budget(budget = 5000, clus_cost = 10,
                                        indiv_cost = 1, pr_treat = 0.5,
                                        icc_Y = 0.1,
                                        total_var = 1)$value[grepl("exp_width", c("necessary_n_clusters","cluster_size","exp_width_of_unstd_conditions_diff","budget"))]
  expect_gt(small, big)
})
