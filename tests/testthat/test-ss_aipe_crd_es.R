## Sample-size planning for cluster-randomized designs with a standardized
## effect size target. These functions use a priori Monte Carlo simulation
## (nrep replications inside an iterative search) and are slow, so we cap
## nrep at a small value and skip the heaviest runs on CRAN.

test_that("ss_aipe_crd_es_n_clusters_fixed_width() returns a tidy table that reaches the target width", {
  skip_on_cran()
  # One planner call carries both the shape assertions and the width
  # assertion: an identical second call used to double the runtime of the
  # slowest file in the suite for no additional coverage.
  target <- 0.40
  res <- suppressMessages(
    ss_aipe_crd_es_n_clusters_fixed_width(width = target, n_individuals = 20, es = 0.5,
                                    icc_Y = 0.10, pr_treat = 0.5,
                                    nrep = 12, seed = 113)
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_n_clusters" %in% res$term)
  expect_true(any(grepl("exp_width", res$term)))
  expect_lte(res$value[grepl("exp_width", res$term)], target + 0.05)
})

test_that("ss_aipe_crd_es_n_clusters_fixed_budget() returns a tidy data.frame within budget", {
  skip_on_cran()
  res <- suppressMessages(
    ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, n_individuals = 20,
                                     clus_cost = 10, indiv_cost = 1,
                                     pr_treat = 0.5, icc_Y = 0.10, es = 0.5,
                                     nrep = 12, seed = 113)
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_n_clusters" %in% res$term)
  expect_lte(res$value[res$term == "budget"], 1000)
})

test_that("ss_aipe_crd_es_*() produces deterministic output for a fixed seed", {
  skip_on_cran()
  res1 <- suppressMessages(
    ss_aipe_crd_es_n_clusters_fixed_width(width = 0.40, n_individuals = 20, es = 0.5,
                                    icc_Y = 0.10, pr_treat = 0.5,
                                    nrep = 6, seed = 113)
  )
  res2 <- suppressMessages(
    ss_aipe_crd_es_n_clusters_fixed_width(width = 0.40, n_individuals = 20, es = 0.5,
                                    icc_Y = 0.10, pr_treat = 0.5,
                                    nrep = 6, seed = 113)
  )
  expect_equal(res1, res2)
})

test_that("ss_aipe_crd_es_n_individuals_fixed_budget() returns a tidy data.frame within budget", {
  skip_on_cran()
  res <- suppressMessages(
    ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 30,
                                      clus_cost = 10, indiv_cost = 1,
                                      pr_treat = 0.5, icc_Y = 0.10, es = 0.5,
                                      nrep = 12, seed = 113)
  )
  expect_s3_class(res, "data.frame")
  expect_true("cluster_size" %in% res$term)
  expect_lte(res$value[res$term == "budget"], 1000)
})

test_that("the CRD effect-size Monte Carlo internals default seed to NULL, not a baked-in constant", {
  # The package's reproducibility-seed discipline: a function that exposes a
  # seed argument defaults to NULL ("use the caller's current RNG state"), never
  # to a constant such as 113. A baked-in default would silently make a default
  # planning call look reproducible and hide the sampling variability of the a
  # priori Monte Carlo. The public ss_aipe_crd_es_*() planners already default
  # seed = NULL and forward it down, so these non-exported helpers must too;
  # this is the fast guard for the seed-discipline grep over R/.
  internals <- c(".find_n_clus_crd_es", ".find_width_crd_es",
                 ".find_n_indiv_crd_es", ".find_min_width_crd_es",
                 ".find_min_cost_crd_es")
  for (fn in internals) {
    expect_null(formals(getFromNamespace(fn, "DMAR"))$seed, info = fn)
  }
})

test_that("the public ss_aipe_crd_es_*() planners default seed to NULL", {
  planners <- c("ss_aipe_crd_es_n_clusters_fixed_width",
                "ss_aipe_crd_es_n_individuals_fixed_width",
                "ss_aipe_crd_es_n_clusters_fixed_budget",
                "ss_aipe_crd_es_n_individuals_fixed_budget",
                "ss_aipe_crd_es_both_fixed_budget",
                "ss_aipe_crd_es_both_fixed_width")
  for (fn in planners) {
    expect_null(formals(get(fn))$seed, info = fn)
  }
})

test_that("ss_aipe_crd_es_both_fixed_budget() reports a feasible n_clusters/cluster_size combination", {
  skip_on_cran()
  res <- suppressMessages(
    ss_aipe_crd_es_both_fixed_budget(budget = 1000, clus_cost = 10,
                                    indiv_cost = 1, es = 0.5, icc_Y = 0.10,
                                    pr_treat = 0.5, nrep = 12, seed = 113)
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("necessary_n_clusters", "cluster_size", "budget") %in% res$term))
  expect_lte(res$value[res$term == "budget"], 1000)
})

test_that("ss_aipe_crd_es_n_individuals_fixed_width() returns a tidy cluster size for a fixed number of clusters", {
  skip_on_cran()
  res <- suppressMessages(
    ss_aipe_crd_es_n_individuals_fixed_width(width = 0.60, n_clusters = 60, es = 0.5,
                                      icc_Y = 0.10, pr_treat = 0.5,
                                      nrep = 12, seed = 113)
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("cluster_size" %in% res$term)
  expect_gte(res$value[res$term == "cluster_size"], 1)
  expect_true(any(grepl("exp_width", res$term)))
})

test_that("ss_aipe_crd_es_both_fixed_width() reports a feasible n_clusters/cluster_size combination", {
  skip_on_cran()
  res <- suppressMessages(
    ss_aipe_crd_es_both_fixed_width(width = 0.60, clus_cost = 10, indiv_cost = 1,
                                    es = 0.5, icc_Y = 0.10, pr_treat = 0.5,
                                    nrep = 12, seed = 113)
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("necessary_n_clusters", "cluster_size") %in% res$term))
  expect_gte(res$value[res$term == "necessary_n_clusters"], 1)
  expect_gte(res$value[res$term == "cluster_size"], 1)
})
