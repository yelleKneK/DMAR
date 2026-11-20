## Sample-size planning for cluster-randomized designs (Lai & Kelley, 2011/2018)
## --------------------------------------------------------------------------
## The non-effect size CRD family targets the population mean DIFFERENCE in
## treatment minus control. Six exported functions handle the three "what
## you control" axes (n_clusters, cluster_size, both) crossed with two "what you
## constrain" axes (fixed_width, fixed_budget).

test_that("ss_aipe_crd_n_clusters_fixed_width() returns a tidy data.frame with n_clusters + expected width", {
  res <- ss_aipe_crd_n_clusters_fixed_width(width = 0.3, cluster_size = 20, pr_treat = 0.5,
                                      true_icc = 0.1, true_total_var = 1)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_n_clusters" %in% res$term)
  expect_true("expected_width" %in% res$term)
  expect_true(res$value[res$term == "necessary_n_clusters"] > 0)
})

test_that("ss_aipe_crd_n_clusters_fixed_width() achieves expected width <= target width", {
  width_target <- 0.30
  res <- ss_aipe_crd_n_clusters_fixed_width(width = width_target, cluster_size = 20,
                                      pr_treat = 0.5, true_icc = 0.1, true_total_var = 1)
  expect_lte(res$value[res$term == "expected_width"], width_target + 1e-6)
})

test_that("ss_aipe_crd_n_clusters_fixed_width() necessary n_clusters grows with higher ICC", {
  low_icc  <- ss_aipe_crd_n_clusters_fixed_width(width = 0.3, cluster_size = 20,
                                           pr_treat = 0.5, true_icc = 0.05,
                                           true_total_var = 1)$value[1]
  high_icc <- ss_aipe_crd_n_clusters_fixed_width(width = 0.3, cluster_size = 20,
                                           pr_treat = 0.5, true_icc = 0.30,
                                           true_total_var = 1)$value[1]
  expect_gt(high_icc, low_icc)
})

test_that("ss_aipe_crd_n_clusters_fixed_width() smaller width target requires more clusters", {
  wide <- ss_aipe_crd_n_clusters_fixed_width(width = 0.50, cluster_size = 20,
                                       pr_treat = 0.5, true_icc = 0.1,
                                       true_total_var = 1)$value[1]
  tight <- ss_aipe_crd_n_clusters_fixed_width(width = 0.20, cluster_size = 20,
                                        pr_treat = 0.5, true_icc = 0.1,
                                        true_total_var = 1)$value[1]
  expect_gt(tight, wide)
})

test_that("ss_aipe_crd_n_individuals_fixed_width() errors when n_clusters is too small to reach the target width", {
  # With 30 clusters the width falls toward 2 qt(.975, 28) sqrt(.1 / 7.5),
  # its value with every cluster mean known exactly, which is above 0.3.
  err <- tryCatch(
    ss_aipe_crd_n_individuals_fixed_width(width = 0.3, n_clusters = 30, pr_treat = 0.5,
                                          true_icc = 0.1, true_total_var = 1),
    error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$reason, "width_not_reachable")
  expect_identical(err$arg, "width")
  expect_equal(err$width_limit, 2 * stats::qt(0.975, 28) * sqrt(0.1 / 7.5))
  expect_gt(err$width_at_max, err$width_limit)
  expect_match(conditionMessage(err), "no cluster size reaches the target", fixed = TRUE)
})

test_that("ss_aipe_crd_n_individuals_fixed_width() returns a tidy data.frame when n_clusters is large enough", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  res <- ss_aipe_crd_n_individuals_fixed_width(width = 0.30, n_clusters = 200,
                                       pr_treat = 0.5, true_icc = 0.1,
                                       true_total_var = 1)
  expect_s3_class(res, "data.frame")
  expect_true("necessary_cluster_size" %in% res$term)
  expect_lte(res$value[res$term == "expected_width"], 0.30 + 1e-6)
})

test_that("ss_aipe_crd_both_fixed_width() jointly returns n_clusters, cluster_size, and a budget total", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  res <- ss_aipe_crd_both_fixed_width(width = 0.30, cluster_cost = 10,
                                     individual_cost = 1, pr_treat = 0.5,
                                     true_icc = 0.1, true_total_var = 1)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("necessary_n_clusters", "necessary_cluster_size", "budget") %in% res$term))
  n_clusters  <- res$value[res$term == "necessary_n_clusters"]
  csize  <- res$value[res$term == "necessary_cluster_size"]
  budget <- res$value[res$term == "budget"]
  expect_equal(budget, n_clusters * 10 + n_clusters * csize * 1, tolerance = 1)
})

test_that("ss_aipe_crd_n_clusters_fixed_budget() respects the budget constraint", {
  res <- ss_aipe_crd_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
                                       cluster_cost = 10, individual_cost = 1,
                                       pr_treat = 0.5, true_icc = 0.1, true_total_var = 1)
  expect_s3_class(res, "data.frame")
  expect_lte(res$value[res$term == "budget"], 1000)
  expect_true(res$value[res$term == "necessary_n_clusters"] > 0)
})

test_that("ss_aipe_crd_n_individuals_fixed_budget() respects the budget constraint", {
  res <- ss_aipe_crd_n_individuals_fixed_budget(budget = 1000, n_clusters = 30,
                                        cluster_cost = 10, individual_cost = 1,
                                        pr_treat = 0.5, true_icc = 0.1, true_total_var = 1)
  expect_s3_class(res, "data.frame")
  expect_lte(res$value[res$term == "budget"], 1000)
  expect_true(res$value[res$term == "necessary_cluster_size"] > 0)
})

test_that("ss_aipe_crd_both_fixed_budget() returns a feasible joint (n_clusters, cluster_size) under budget", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  res <- ss_aipe_crd_both_fixed_budget(budget = 1000, cluster_cost = 10,
                                      individual_cost = 1, pr_treat = 0.5,
                                      true_icc = 0.1, true_total_var = 1)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("necessary_n_clusters", "necessary_cluster_size", "budget") %in% res$term))
  expect_lte(res$value[res$term == "budget"], 1000)
})

.crd_row <- function(x, term) x$value[x$term == term]

test_that("ss_aipe_crd_*_fixed_budget() expected widths shrink as the budget grows", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  small <- .crd_row(ss_aipe_crd_both_fixed_budget(budget =  500, cluster_cost = 10,
                                        individual_cost = 1, pr_treat = 0.5,
                                        true_icc = 0.1,
                                        true_total_var = 1), "expected_width")
  big   <- .crd_row(ss_aipe_crd_both_fixed_budget(budget = 5000, cluster_cost = 10,
                                        individual_cost = 1, pr_treat = 0.5,
                                        true_icc = 0.1,
                                        true_total_var = 1), "expected_width")
  expect_gt(small, big)
})

test_that("the width row is expected_width, or assured_width with the assurance in its own row", {
  skip_on_cran()
  plain <- ss_aipe_crd_n_clusters_fixed_width(width = 0.3, cluster_size = 20, pr_treat = 0.5,
                                              true_icc = 0.1, true_total_var = 1)
  expect_identical(plain$term, c("necessary_n_clusters", "width", "expected_width"))
  expect_identical(.crd_row(plain, "width"), 0.3)
  assured <- ss_aipe_crd_n_clusters_fixed_width(width = 0.3, cluster_size = 20, pr_treat = 0.5,
                                                true_icc = 0.1, true_total_var = 1,
                                                assurance = 0.8)
  expect_identical(assured$term, c("necessary_n_clusters", "width", "assured_width", "assurance"))
  expect_identical(.crd_row(assured, "assurance"), 0.8)
  expect_lte(.crd_row(assured, "assured_width"), 0.3 + 1e-6)
  expect_gte(.crd_row(assured, "necessary_n_clusters"), .crd_row(plain, "necessary_n_clusters"))
  # Every member reports the same row names; no row name carries a value.
  both <- ss_aipe_crd_both_fixed_budget(budget = 1000, cluster_cost = 10, individual_cost = 1,
                                        pr_treat = 0.5, true_icc = 0.1, true_total_var = 1,
                                        assurance = 0.9)
  expect_identical(both$term, c("necessary_n_clusters", "necessary_cluster_size", "assured_width",
                                "assurance", "budget"))
  expect_identical(.crd_row(both, "assurance"), 0.9)
  gl <- generics::glance(both)
  expect_true(all(c("assured_width", "assurance", "budget") %in% names(gl)))
})

test_that("the planning values of the cluster randomized planners carry true_ and the 1.0.x names are refused", {
  skip_on_cran()
  crd <- c("ss_aipe_crd_n_clusters_fixed_width", "ss_aipe_crd_n_individuals_fixed_width",
           "ss_aipe_crd_n_clusters_fixed_budget", "ss_aipe_crd_n_individuals_fixed_budget",
           "ss_aipe_crd_both_fixed_budget", "ss_aipe_crd_both_fixed_width")
  old <- c("num_predictors", "n_individuals", "icc_Y", "tau_Y", "sigma2_Y", "total_var",
           "R2_between", "R2_within", "clus_cost", "indiv_cost")
  for (fn in crd) {
    fm <- names(formals(get(fn, asNamespace("DMAR"))))
    expect_false(any(old %in% fm), info = fn)
    expect_true(all(c("n_covariates", "true_icc", "true_var_between", "true_var_within", "true_total_var",
                      "true_R2_between", "true_R2_within", "cluster_cost", "individual_cost") %in% fm),
                info = fn)
  }
  expect_true("cluster_size" %in% names(formals(ss_aipe_crd_n_clusters_fixed_width)))
  expect_true("cluster_size" %in% names(formals(ss_aipe_crd_n_clusters_fixed_budget)))
  # do.call() with the 1.0.x name in the list: R must find no such argument.
  expect_error(do.call(ss_aipe_crd_n_clusters_fixed_width,
                       list(width = 0.3, n_individuals = 20, pr_treat = 0.5,
                            true_icc = 0.1, true_total_var = 1)),
               "unused argument")
  expect_error(do.call(ss_aipe_crd_n_clusters_fixed_width,
                       list(width = 0.3, cluster_size = 20, pr_treat = 0.5,
                            icc_Y = 0.1, true_total_var = 1)),
               "unused argument")
  expect_error(do.call(ss_aipe_crd_n_clusters_fixed_budget,
                       list(budget = 1000, cluster_size = 20, clus_cost = 10,
                            individual_cost = 1)),
               "unused argument")
})

test_that("n_covariates counts the covariates beside the treatment indicator", {
  # The degrees of freedom are n_clusters - 2 - n_covariates, so with an assurance
  # the width that is not exceeded grows as covariates use up degrees of freedom
  # when they explain nothing.
  w0 <- .crd_row(ss_aipe_crd_n_individuals_fixed_budget(budget = 1000, n_clusters = 20,
          cluster_cost = 10, individual_cost = 1, pr_treat = 0.5, true_icc = 0.1,
          true_total_var = 1, assurance = 0.9, n_covariates = 0), "assured_width")
  w3 <- .crd_row(ss_aipe_crd_n_individuals_fixed_budget(budget = 1000, n_clusters = 20,
          cluster_cost = 10, individual_cost = 1, pr_treat = 0.5, true_icc = 0.1,
          true_total_var = 1, assurance = 0.9, n_covariates = 3), "assured_width")
  expect_gt(w3, w0)
})

test_that("tidy() and glance() find the planned size and the desired width of every cluster randomized plan", {
  skip_on_cran()
  # The size solved for is a necessary_ row; a size that was given is not
  # echoed. width is the width asked for, NA for a plan made from a budget.
  pop <- list(pr_treat = 0.5, true_var_between = 0.25, true_var_within = 0.75)
  plans <- list(
    n_clusters_width = do.call(ss_aipe_crd_n_clusters_fixed_width,
                               c(list(width = 0.3, cluster_size = 20), pop)),
    n_clusters_width_assured = do.call(ss_aipe_crd_n_clusters_fixed_width,
                                       c(list(width = 0.3, cluster_size = 20, assurance = 0.9), pop)),
    n_individuals_width = do.call(ss_aipe_crd_n_individuals_fixed_width,
                                  c(list(width = 0.3, n_clusters = 250), pop)),
    n_individuals_width_assured = do.call(ss_aipe_crd_n_individuals_fixed_width,
                                          c(list(width = 0.3, n_clusters = 250, assurance = 0.9), pop)),
    n_clusters_budget = do.call(ss_aipe_crd_n_clusters_fixed_budget,
                                c(list(budget = 10000, cluster_size = 20, cluster_cost = 10,
                                       individual_cost = 1), pop)),
    n_individuals_budget = do.call(ss_aipe_crd_n_individuals_fixed_budget,
                                   c(list(budget = 10000, n_clusters = 100, cluster_cost = 10,
                                          individual_cost = 1), pop)),
    both_width = do.call(ss_aipe_crd_both_fixed_width,
                         c(list(width = 0.3, cluster_cost = 10, individual_cost = 1), pop))
  )
  size_row <- c(n_clusters_width = "necessary_n_clusters",
                n_clusters_width_assured = "necessary_n_clusters",
                n_individuals_width = "necessary_cluster_size",
                n_individuals_width_assured = "necessary_cluster_size",
                n_clusters_budget = "necessary_n_clusters",
                n_individuals_budget = "necessary_cluster_size",
                both_width = "necessary_n_clusters")
  for (nm in names(plans)) {
    x <- plans[[nm]]
    td <- generics::tidy(x)
    gl <- generics::glance(x)
    expect_false("cluster_size" %in% x$term, info = nm)
    expect_identical(td$estimate, .crd_row(x, size_row[[nm]]), info = nm)
    expect_false(is.na(td$estimate), info = nm)
    expect_identical(gl$estimate, td$estimate, info = nm)
    if (grepl("width", nm)) {
      expect_identical(td$width, 0.3, info = nm)
      expect_identical(gl$width, 0.3, info = nm)
    } else {
      expect_true(is.na(td$width), info = nm)
    }
    # The width the plan gives stays under its own name.
    wrow <- if (grepl("assured", nm)) "assured_width" else "expected_width"
    expect_true(wrow %in% names(gl), info = nm)
    expect_identical(gl[[wrow]], .crd_row(x, wrow), info = nm)
  }
  expect_true("necessary_cluster_size" %in% names(generics::glance(plans$both_width)))
})

## Boundaries, validation, and the variance model (the round 6 review of
## November 2026: SQ01 to SQ07 and SQ10, and the problems its evaluation
## found). Every call here is closed form and runs at once.

.crd_pop <- list(true_var_between = 0.2, true_var_within = 0.8)
.crd_call <- function(fn, ...) do.call(fn, c(list(...), .crd_pop))
.crd_err <- function(expr) tryCatch(expr, error = function(e) e)

test_that("a budget equal to the cost of the design at the end of a bracket buys that design", {
  # 1000 clusters of 20 at a cost of 1 each cost 20000, and 10 clusters of
  # 1000 cost 10000. Through the first builds of 1.1.0 both stopped with
  # "'from' must be a finite number".
  a <- ss_aipe_crd_n_clusters_fixed_budget(budget = 20000, cluster_size = 20,
                                           cluster_cost = 0, individual_cost = 1)
  expect_identical(.crd_row(a, "necessary_n_clusters"), 1000)
  expect_identical(.crd_row(a, "budget"), 20000)
  b <- ss_aipe_crd_n_individuals_fixed_budget(budget = 10000, n_clusters = 10,
                                              cluster_cost = 0, individual_cost = 1)
  expect_identical(.crd_row(b, "necessary_cluster_size"), 1000)
  expect_identical(.crd_row(b, "budget"), 10000)
  # With a cost per cluster as well, and one more cluster or individual is
  # beyond the budget.
  c2 <- ss_aipe_crd_n_clusters_fixed_budget(budget = 2000 * 30, cluster_size = 20,
                                            cluster_cost = 10, individual_cost = 1)
  expect_identical(.crd_row(c2, "necessary_n_clusters"), 2000)
})

test_that("a target width equal to the width at the end of a bracket is met there", {
  ns <- asNamespace("DMAR")
  w1 <- ns$.find_width_crd_diff(1000, 20, 0.5, true_var_between = 0.2, true_var_within = 0.8)
  a <- .crd_call(ss_aipe_crd_n_clusters_fixed_width, width = w1, cluster_size = 20, pr_treat = 0.5)
  expect_identical(.crd_row(a, "necessary_n_clusters"), 1000)
  expect_identical(.crd_row(a, "expected_width"), w1)
  w2 <- ns$.find_width_crd_diff(100, 1000, 0.5, true_var_between = 0.2, true_var_within = 0.8)
  b <- .crd_call(ss_aipe_crd_n_individuals_fixed_width, width = w2, n_clusters = 100, pr_treat = 0.5)
  expect_identical(.crd_row(b, "necessary_cluster_size"), 1000)
  expect_identical(.crd_row(b, "expected_width"), w2)
})

test_that("the search for the number of clusters keeps two clusters in each condition and returns a size it tested", {
  skip_on_cran()
  # A wide target is met by every number of clusters from 100 to 1000. The
  # search returned 3 + n_covariates, which it had not tested (at pr_treat
  # = .1, 3 clusters of which none is treated, with an infinite width), and
  # stopped with "wrong sign in 'by' argument" at pr_treat = .01, where the
  # floor ceiling(2 / .01) = 200 is above 100.
  for (p in c(0.1, 0.01, 0.005, 0.5, 0.9)) for (k in 0:1) {
    lo <- ceiling(2 * (1 / min(p, 1 - p))) + k
    res <- .crd_call(ss_aipe_crd_n_clusters_fixed_width, width = 100, cluster_size = 20,
                     pr_treat = p, n_covariates = k)
    info <- paste("pr_treat", p, "n_covariates", k)
    expect_identical(.crd_row(res, "necessary_n_clusters"), lo, info = info)
    expect_lte(.crd_row(res, "expected_width"), 100)
    expect_true(is.finite(.crd_row(res, "expected_width")), info = info)
  }
  # Unequal cluster sizes too: the treated clusters' terms were summed over
  # 1:0 with no treated cluster, which gave NaN.
  # (The cluster sizes are printed ahead of the table.)
  utils::capture.output(res <- .crd_call(ss_aipe_crd_n_clusters_fixed_width, width = 100,
                                         cluster_size = 20, pr_treat = 0.1,
                                         diff_size = c(-2, 0, 2)))
  expect_identical(.crd_row(res, "necessary_n_clusters"), 20)
  expect_true(is.finite(.crd_row(res, "expected_width")))
  ns <- asNamespace("DMAR")
  expect_identical(ns$.var_crd_diff(5, ns$.find_n_indiv_vec(20, c(-2, 0, 2), 5), 0, 0.2, 0.8), Inf)
})

test_that("the plan for the number of clusters is the first number from the floor that meets the target", {
  skip_on_cran()
  # Enumeration of every number of clusters from the floor, with equal
  # cluster sizes, where the width falls as clusters are added.
  ns <- asNamespace("DMAR")
  for (p in c(0.02, 0.1, 0.3, 0.9)) for (k in c(0, 2)) for (a in list(NULL, 0.8)) {
    lo <- ceiling(2 * (1 / min(p, 1 - p))) + k
    n <- as.numeric(lo:20000)
    w_n <- ns$.find_width_crd_diff(n, 20, p, true_var_between = 0.2, true_var_within = 0.8,
                                   n_covariates = k, assurance = a)
    for (target in c(0.15, 0.4, 1.5, 6)) {
      res <- .crd_call(ss_aipe_crd_n_clusters_fixed_width, width = target, cluster_size = 20,
                       pr_treat = p, n_covariates = k, assurance = a)
      expect_identical(.crd_row(res, "necessary_n_clusters"), n[which(w_n <= target)[1]],
                       info = paste(p, k, format(a), target))
    }
  }
})

test_that("the planners that search over both sizes keep the floor and stop with a reason where nothing is affordable", {
  skip_on_cran()
  # pr_treat = .1 puts the floor at 20 clusters; 20 clusters of 2 cost
  # 20,040 at 1000 per cluster, more than a budget of 20,000. Through the
  # first builds of 1.1.0 this stopped with "invalid 'type' (list) of
  # argument", and pr_treat = .005 (a floor of 400) with "wrong sign in
  # 'by' argument".
  err <- .crd_err(.crd_call(ss_aipe_crd_both_fixed_budget, budget = 20000, cluster_cost = 1000,
                            individual_cost = 1, pr_treat = 0.1))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$reason, "budget_infeasible")
  expect_identical(err$arg, "budget")
  expect_identical(c(err$n_clusters, err$cluster_size, err$min_cost), c(20, 2, 20040))
  expect_match(conditionMessage(err), "at least two clusters in each condition at pr_treat = 0.1",
               fixed = TRUE)
  err <- .crd_err(.crd_call(ss_aipe_crd_both_fixed_budget, budget = 20000, cluster_cost = 1000,
                            individual_cost = 1, pr_treat = 0.005))
  expect_identical(err$reason, "budget_infeasible")
  expect_identical(err$min_cost, 400800)
  err <- .crd_err(.crd_call(ss_aipe_crd_both_fixed_budget, budget = 1, cluster_cost = 10,
                            individual_cost = 1, pr_treat = 0.5))
  expect_identical(err$reason, "budget_infeasible")
  expect_identical(err$min_cost, 48)
  # A budget that buys the floor's design plans at or above it.
  res <- .crd_call(ss_aipe_crd_both_fixed_budget, budget = 30000, cluster_cost = 1000,
                   individual_cost = 1, pr_treat = 0.1)
  expect_gte(.crd_row(res, "necessary_n_clusters"), 20)
  expect_lte(.crd_row(res, "budget"), 30000)
  res <- .crd_call(ss_aipe_crd_both_fixed_width, width = 3, cluster_cost = 1000,
                   individual_cost = 1, pr_treat = 0.005)
  expect_identical(.crd_row(res, "necessary_n_clusters"), 400)
  expect_lte(.crd_row(res, "expected_width"), 3)
})

test_that("every cluster randomized planner returns a numeric value column", {
  skip_on_cran()
  # .find_n_indiv_crd_diff() returned the string "> 100000" for a number
  # of clusters that could not reach the width, which made the value
  # column of ss_aipe_crd_both_fixed_width() character ("80", "7", ...).
  x <- .crd_call(ss_aipe_crd_both_fixed_width, width = 0.5, cluster_cost = 10,
                 individual_cost = 1, pr_treat = 0.5)
  expect_type(x$value, "double")
  expect_identical(x$value[c(1, 2, 5)], c(80, 7, 1360))
  expect_type(generics::tidy(x)$estimate, "double")
  expect_type(generics::glance(x)$budget, "double")
  expect_identical(asNamespace("DMAR")$.find_n_indiv_crd_diff(
    width = 0.5, n_clusters = 10, pr_treat = 0.5, true_var_between = 0.2, true_var_within = 0.8), Inf)
  plans <- list(
    .crd_call(ss_aipe_crd_n_clusters_fixed_width, width = 0.5, cluster_size = 20, pr_treat = 0.5,
              cluster_cost = 10, individual_cost = 1),
    .crd_call(ss_aipe_crd_n_individuals_fixed_width, width = 0.5, n_clusters = 80, pr_treat = 0.5,
              cluster_cost = 10, individual_cost = 1),
    .crd_call(ss_aipe_crd_n_clusters_fixed_budget, budget = 2000, cluster_size = 20, cluster_cost = 10,
              individual_cost = 1, pr_treat = 0.5),
    .crd_call(ss_aipe_crd_n_individuals_fixed_budget, budget = 2000, n_clusters = 40, cluster_cost = 10,
              individual_cost = 1, pr_treat = 0.5),
    .crd_call(ss_aipe_crd_both_fixed_budget, budget = 2000, cluster_cost = 10, individual_cost = 1,
              pr_treat = 0.5),
    x,
    ss_aipe_crd_es_n_clusters_fixed_budget(budget = 2000, cluster_size = 20, cluster_cost = 10,
                                           individual_cost = 1),
    ss_aipe_crd_es_n_individuals_fixed_budget(budget = 2000, n_clusters = 40, cluster_cost = 10,
                                              individual_cost = 1))
  for (i in seq_along(plans)) expect_type(plans[[i]]$value, "double")
})

test_that("a target no cluster size reaches with the number of clusters is refused with its reason", {
  err <- .crd_err(.crd_call(ss_aipe_crd_n_individuals_fixed_width, width = 0.5, n_clusters = 10,
                            pr_treat = 0.5))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$reason, "width_not_reachable")
  expect_equal(err$width_limit, 2 * stats::qt(0.975, 8) * sqrt(0.2 / 2.5))
  expect_identical(err$max_cluster_size, 1e5)
  # A target above the limit and below the width at 100,000 per cluster:
  # a larger cluster size would reach it.
  lim <- 2 * stats::qt(0.975, 8) * sqrt(0.002 / 2.5)
  err <- .crd_err(ss_aipe_crd_n_individuals_fixed_width(width = lim * 1.000001, n_clusters = 10,
                    pr_treat = 0.5, true_var_between = 0.002, true_var_within = 0.998))
  expect_identical(err$reason, "width_not_reachable")
  expect_match(conditionMessage(err), "a cluster size above 100,000 would reach the target", fixed = TRUE)
})

test_that("a budget that buys a design with a condition without clusters, or too few clusters, is refused", {
  # A budget of 100 buys 5 clusters of 20, and round(5 * .1) = 0 are
  # treated; the table had an expected width of Inf.
  err <- .crd_err(.crd_call(ss_aipe_crd_n_clusters_fixed_budget, budget = 100, cluster_size = 20,
                            cluster_cost = 0, individual_cost = 1, pr_treat = 0.1))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$reason, "arm_without_clusters")
  expect_identical(c(err$n_clusters, err$n_treated, err$fewest, err$fewest_cost), c(5, 0, 6, 120))
  expect_match(conditionMessage(err), "a budget of 100 buys 5 clusters of 20 individuals (cost 100)",
               fixed = TRUE)
  # Without pr_treat no width is computed and the design is returned.
  ok <- ss_aipe_crd_n_clusters_fixed_budget(budget = 100, cluster_size = 20, cluster_cost = 0,
                                            individual_cost = 1)
  expect_identical(.crd_row(ok, "necessary_n_clusters"), 5)
  # Four clusters with two covariates leave no degrees of freedom.
  err <- .crd_err(.crd_call(ss_aipe_crd_n_clusters_fixed_budget, budget = 80, cluster_size = 20,
                            cluster_cost = 0, individual_cost = 1, pr_treat = 0.5, n_covariates = 2))
  expect_identical(err$reason, "too_few_clusters")
  expect_identical(c(err$fewest, err$fewest_cost), c(5, 100))
  # A fixed number of clusters with an empty condition is refused on entry.
  for (fn in list(ss_aipe_crd_n_individuals_fixed_budget, ss_aipe_crd_n_individuals_fixed_width)) {
    args <- c(list(n_clusters = 5, pr_treat = 0.1), .crd_pop)
    args <- c(args, if ("budget" %in% names(formals(fn))) list(budget = 1000) else list(width = 1))
    err <- .crd_err(do.call(fn, args))
    expect_identical(err$reason, "arm_without_clusters")
    expect_identical(err$arg, c("n_clusters", "pr_treat"))
    expect_identical(err$fewest, 6)
  }
})

test_that("costs and budgets are checked, and costs of 0 are refused where they decide the design", {
  skip_on_cran()
  # Costs of 0 searched without end in the budget planners, and
  # ss_aipe_crd_both_fixed_width() returned a table of NA.
  zero <- list(
    function(cc, ic) ss_aipe_crd_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
                                                         cluster_cost = cc, individual_cost = ic),
    function(cc, ic) ss_aipe_crd_n_individuals_fixed_budget(budget = 1000, n_clusters = 20,
                                                            cluster_cost = cc, individual_cost = ic),
    function(cc, ic) .crd_call(ss_aipe_crd_both_fixed_budget, budget = 10000, cluster_cost = cc,
                               individual_cost = ic, pr_treat = 0.5),
    function(cc, ic) .crd_call(ss_aipe_crd_both_fixed_width, width = 0.5, cluster_cost = cc,
                               individual_cost = ic, pr_treat = 0.5))
  for (i in seq_along(zero)) {
    err <- .crd_err(zero[[i]](0, 0))
    expect_identical(err$reason, "costs_all_zero", info = i)
    expect_identical(err$arg, c("cluster_cost", "individual_cost"), info = i)
    if (i > 1) {
      err <- .crd_err(zero[[i]](10, 0))
      expect_identical(err$reason, "individual_cost_zero", info = i)
      expect_identical(err$arg, "individual_cost", info = i)
    }
  }
  # Clusters of a fixed size cost only their clusters: a finite plan.
  expect_identical(.crd_row(zero[[1]](10, 0), "necessary_n_clusters"), 100)
  # Where the costs only price the plan, 0 is a price, and a budget of 0 is
  # reported.
  res <- .crd_call(ss_aipe_crd_n_clusters_fixed_width, width = 0.5, cluster_size = 20,
                   pr_treat = 0.5, cluster_cost = 0, individual_cost = 0)
  expect_identical(.crd_row(res, "budget"), 0)
  for (bad in list(-1, NA_real_, c(1, 2), "1", Inf)) {
    expect_error(.crd_call(ss_aipe_crd_n_clusters_fixed_width, width = 0.5, cluster_size = 20,
                           pr_treat = 0.5, cluster_cost = bad, individual_cost = 1),
                 class = "dmar_argument_invalid")
    expect_error(ss_aipe_crd_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
                                                     individual_cost = bad),
                 class = "dmar_argument_invalid")
    expect_error(ss_aipe_crd_n_clusters_fixed_budget(budget = bad, cluster_size = 20),
                 class = "dmar_argument_invalid")
  }
})

test_that("a budget below the cost of the smallest design is refused with that cost", {
  # The table had NA for the size and the budget.
  err <- .crd_err(ss_aipe_crd_n_individuals_fixed_budget(budget = 1, n_clusters = 10,
                                                         cluster_cost = 10, individual_cost = 1))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$reason, "budget_infeasible")
  expect_identical(c(err$n_clusters, err$cluster_size, err$min_cost), c(10, 2, 120))
  expect_match(conditionMessage(err), "Give a budget of at least 120", fixed = TRUE)
  err <- .crd_err(ss_aipe_crd_n_clusters_fixed_budget(budget = 1, cluster_size = 10,
                                                      cluster_cost = 10, individual_cost = 1))
  expect_identical(err$reason, "budget_infeasible")
  expect_identical(c(err$n_clusters, err$cluster_size, err$min_cost), c(4, 10, 80))
  # The smallest design itself is bought.
  ok <- ss_aipe_crd_n_clusters_fixed_budget(budget = 80, cluster_size = 10, cluster_cost = 10,
                                            individual_cost = 1)
  expect_identical(.crd_row(ok, "necessary_n_clusters"), 4)
})

test_that("counts, proportions, and the allocation are checked on entry", {
  skip_on_cran()
  base <- list(width = 1, cluster_size = 20, pr_treat = 0.5)
  refused <- function(changes, arg) {
    err <- .crd_err(do.call(ss_aipe_crd_n_clusters_fixed_width,
                            utils::modifyList(c(base, .crd_pop), changes)))
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, arg, info = paste(names(changes), format(changes[[1]])))
  }
  # n_covariates = .5 planned 17.5 clusters and true_R2_between = 1.1
  # planned 5; cluster_size = 20.5 was planned for.
  for (v in list(0.5, -1, NA_real_, c(0, 1))) refused(list(n_covariates = v), "n_covariates")
  for (nm in c("true_R2_between", "true_R2_within")) {
    for (v in list(1.1, -0.5, 1, NA_real_)) refused(stats::setNames(list(v), nm), nm)
  }
  for (v in list(20.5, 1, 0, c(20, 30))) refused(list(cluster_size = v), "cluster_size")
  # pr_treat = 1 (and 0) searched without end.
  for (v in list(0, 1, -0.1, 1.5)) refused(list(pr_treat = v), "pr_treat")
  # n_clusters: whole, and at least 3 + n_covariates (2 stopped with
  # "missing value where TRUE/FALSE needed").
  for (v in list(60.5, 2)) {
    err <- .crd_err(.crd_call(ss_aipe_crd_n_individuals_fixed_width, width = 1, n_clusters = v,
                              pr_treat = 0.5))
    expect_identical(err$arg, "n_clusters")
  }
  err <- .crd_err(.crd_call(ss_aipe_crd_n_individuals_fixed_budget, budget = 1000, n_clusters = 4,
                            pr_treat = 0.5, n_covariates = 2))
  expect_identical(err$arg, "n_clusters")
  expect_match(conditionMessage(err), "at least 5 (3 + n_covariates)", fixed = TRUE)
  expect_s3_class(.crd_call(ss_aipe_crd_n_individuals_fixed_budget, budget = 1000, n_clusters = 5,
                            pr_treat = 0.5, n_covariates = 2), "dmar_tbl")
})

test_that("a between-level R2 without a cluster-level covariate is refused by the six planners", {
  skip_on_cran()
  # With n_covariates = 0, true_R2_between = .9 shrank the plan for a width
  # of 1 at clusters of 20 from 18 clusters to 7.
  calls <- list(
    ss_aipe_crd_n_clusters_fixed_width = list(width = 1, cluster_size = 20, pr_treat = 0.5),
    ss_aipe_crd_n_individuals_fixed_width = list(width = 1, n_clusters = 20, pr_treat = 0.5),
    ss_aipe_crd_n_clusters_fixed_budget = list(budget = 1000, cluster_size = 20, pr_treat = 0.5),
    ss_aipe_crd_n_individuals_fixed_budget = list(budget = 1000, n_clusters = 20, pr_treat = 0.5),
    ss_aipe_crd_both_fixed_budget = list(budget = 1000, cluster_cost = 10, pr_treat = 0.5),
    ss_aipe_crd_both_fixed_width = list(width = 1, cluster_cost = 10, pr_treat = 0.5))
  for (fn in names(calls)) {
    args <- c(calls[[fn]], .crd_pop)
    err <- .crd_err(do.call(fn, c(args, list(true_R2_between = 0.9))))
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "true_R2_between", info = fn)
    expect_identical(err$reason, "R2_between_without_covariates", info = fn)
    expect_identical(err$value, 0.9, info = fn)
    # With a covariate the share is planned for, and true_R2_within needs none.
    expect_s3_class(do.call(fn, c(args, list(true_R2_between = 0.9, n_covariates = 1))), "dmar_tbl")
    expect_s3_class(do.call(fn, c(args, list(true_R2_within = 0.5))), "dmar_tbl")
  }
  msg <- conditionMessage(err)
  expect_match(msg, "share of the between-cluster variance", fixed = TRUE)
  expect_match(msg, "n_covariates = 0 says there are none", fixed = TRUE)
  expect_match(msg, "Give n_covariates", fixed = TRUE)
  expect_match(msg, "set true_R2_between = 0", fixed = TRUE)
  # A budget planner that computes no width checks it as well.
  err <- .crd_err(ss_aipe_crd_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
                                                      true_R2_between = 0.9))
  expect_identical(err$reason, "R2_between_without_covariates")
})

test_that("any two of the four variance inputs give the same plan, and inputs that disagree are refused", {
  skip_on_cran()
  # The documented pairs (between, ICC) and (within, ICC) stopped with
  # "missing value where TRUE/FALSE needed".
  pairs <- list(list(true_var_between = 0.2, true_var_within = 0.8),
                list(true_total_var = 1, true_icc = 0.2),
                list(true_total_var = 1, true_var_between = 0.2),
                list(true_total_var = 1, true_var_within = 0.8),
                list(true_var_between = 0.2, true_icc = 0.2),
                list(true_var_within = 0.8, true_icc = 0.2),
                list(true_var_between = 0.2, true_var_within = 0.8, true_total_var = 1, true_icc = 0.2))
  for (v in pairs) {
    res <- do.call(ss_aipe_crd_n_clusters_fixed_width,
                   c(list(width = 1, cluster_size = 20, pr_treat = 0.5), v))
    expect_identical(.crd_row(res, "necessary_n_clusters"), 18, info = paste(names(v), collapse = " "))
  }
  plan <- function(...) .crd_err(ss_aipe_crd_n_clusters_fixed_width(width = 1, cluster_size = 20,
                                                                    pr_treat = 0.5, ...))
  # A negative variance is refused before the search (ten warnings and an
  # error from R).
  err <- expect_no_warning(plan(true_var_between = -0.2, true_var_within = 0.8))
  expect_identical(err$reason, "variance_out_of_range")
  expect_identical(err$arg, "true_var_between")
  # Total 1, ICC .2, and a between-cluster variance of .5 disagree; the
  # call planned 36 clusters with between .5 and within .8.
  err <- plan(true_total_var = 1, true_icc = 0.2, true_var_between = 0.5)
  expect_identical(err$reason, "inconsistent_variances")
  expect_identical(err$arg, c("true_total_var", "true_icc", "true_var_between"))
  expect_identical(err$supplied$true_var_between, 0.5)
  expect_match(conditionMessage(err), "imply true_var_between = 0.2, but 0.5 was given", fixed = TRUE)
  # A rounded ICC beside two variances shows the value they imply.
  err <- plan(true_var_between = 1, true_var_within = 2, true_icc = 0.333)
  expect_match(conditionMessage(err), "imply true_icc = 0.333333333333333, but 0.333 was given",
               fixed = TRUE)
  expect_identical(plan(true_var_between = 0.2)$reason, "underidentified_variances")
  expect_identical(plan(true_var_between = 0, true_icc = 0)$reason, "underidentified_variances")
  expect_identical(plan(true_var_within = 0.5, true_icc = 1)$reason, "inconsistent_variances")
  expect_identical(plan(true_total_var = 1, true_var_between = 2)$reason, "inconsistent_variances")
  expect_identical(plan(true_total_var = 0, true_icc = 0.2)$reason, "variance_out_of_range")
  expect_identical(plan(true_icc = 1.2, true_total_var = 1)$reason, "variance_out_of_range")
  expect_identical(plan(true_var_between = NA_real_, true_var_within = 1)$arg, "true_var_between")
  # With unequal cluster sizes a within-cluster variance of 0 gave NaN.
  err <- plan(true_var_between = 0.2, true_var_within = 0, diff_size = c(-2, 0, 2))
  expect_identical(err$reason, "within_variance_zero")
  expect_s3_class(plan(true_var_between = 0.2, true_var_within = 0), "dmar_tbl")
  # The variances are resolved only where a width is computed.
  expect_s3_class(ss_aipe_crd_n_clusters_fixed_budget(budget = 1000, cluster_size = 20), "dmar_tbl")
})

test_that("the variance resolver checks redundant inputs on the scale of the total variance", {
  skip_on_cran()
  resolve <- getFromNamespace(".crd_resolve_variances", "DMAR")
  all4 <- list(true_var_between = 0.2, true_var_within = 0.8, true_total_var = 1, true_icc = 0.2)
  for (pair in utils::combn(names(all4), 2, simplify = FALSE)) {
    got <- do.call(resolve, all4[pair])
    expect_equal(unclass(got)[names(all4)], all4, ignore_attr = TRUE)
    expect_setequal(attr(got, "identified_by"), pair)
  }
  for (s in c(1e-200, 1e-6, 1, 1e6, 1e200)) {
    got <- resolve(0.2 * s, 0.8 * s, s, 0.2)
    expect_equal(got$true_total_var / s, 1)
    expect_error(resolve(0.2 * s, 0.8 * s, 1.001 * s, 0.2), class = "dmar_argument_invalid")
  }
  # Within 64 machine epsilons a redundant total agrees; the variances
  # given are used as given.
  got <- resolve(0.2, 0.8, 1 + 4 * .Machine$double.eps, 0.2)
  expect_identical(c(got$true_var_between, got$true_var_within), c(0.2, 0.8))
  expect_error(resolve(0.2, 0.8, 1 + 4 * .Machine$double.eps, 0.2, rel_tol = 0),
               class = "dmar_argument_invalid")
  for (rho in c(0, 1)) {
    got <- resolve(true_total_var = 1, true_icc = rho)
    expect_identical(c(got$true_var_between, got$true_var_within), c(rho, 1 - rho))
  }
  # A total that is not finite, or 0 (two variances of 0, which planned 3
  # clusters with an expected width of 0), is out of range, as a total of 0
  # given as true_total_var is.
  expect_identical(tryCatch(resolve(1e308, 1e308), error = function(e) e$reason),
                   "variance_out_of_range")
  expect_identical(tryCatch(resolve(0, 0), error = function(e) e$reason),
                   "variance_out_of_range")
  expect_identical(tryCatch(resolve(true_var_within = 0, true_icc = 0.2), error = function(e) e$reason),
                   "variance_out_of_range")
  for (bad in list(-1, NA, NaN, Inf, TRUE, "1", numeric(0), c(0.2, 0.8), list(0.2))) {
    expect_error(resolve(true_total_var = 1, true_icc = bad), class = "dmar_argument_invalid")
  }
})

test_that("a target equal to the width of a plan is met by that plan", {
  skip_on_cran()
  # A width meets the target when it is no larger than it, at the end of a
  # bracket and elsewhere. The candidate asked for a width below the target
  # inside a bracket, so the width a plan reported, given back as the
  # target, gave the next size: 19 clusters, a cluster size of 5, 23
  # clusters, and 19 clusters of 13 at a budget of 817.
  again <- function(fn, args, term, width_row = "expected_width") {
    first <- do.call(fn, c(args, .crd_pop))
    args$width <- .crd_row(first, width_row)
    second <- do.call(fn, c(args, .crd_pop))
    c(.crd_row(first, term), .crd_row(second, term))
  }
  expect_identical(again(ss_aipe_crd_n_clusters_fixed_width,
                         list(width = 1, cluster_size = 20, pr_treat = 0.5), "necessary_n_clusters"),
                   c(18, 18))
  expect_identical(again(ss_aipe_crd_n_individuals_fixed_width,
                         list(width = 1, n_clusters = 30, pr_treat = 0.5), "necessary_cluster_size"),
                   c(4, 4))
  expect_identical(again(ss_aipe_crd_n_clusters_fixed_width,
                         list(width = 1, cluster_size = 20, pr_treat = 0.5, assurance = 0.8),
                         "necessary_n_clusters", "assured_width"), c(22, 22))
  first <- .crd_call(ss_aipe_crd_both_fixed_width, width = 1, cluster_cost = 30, individual_cost = 1,
                     pr_treat = 0.5)
  second <- .crd_call(ss_aipe_crd_both_fixed_width, width = .crd_row(first, "expected_width"),
                      cluster_cost = 30, individual_cost = 1, pr_treat = 0.5)
  expect_identical(second$value, first$value[c(1, 2, 4, 4, 5)])
  expect_identical(.crd_row(second, "budget"), 800)
})

test_that("a width that is not a number stops a search with its reason", {
  ns <- asNamespace("DMAR")
  # The public planners check the variances first; a search given a
  # negative variance had gone on without end in the first build of this
  # repair, reading the widths of NaN as widths that do not meet the target.
  # The square root of a negative variance gives R's warning "NaNs
  # produced", which is not at issue here.
  nan_args <- list(pr_treat = 0.5, true_var_between = -0.5, true_var_within = 0.8)
  checks <- list(
    function() do.call(ns$.find_n_clus_crd_diff, c(list(width = 1, cluster_size = 20, fn = "planner_x"), nan_args)),
    function() do.call(ns$.find_n_indiv_crd_diff, c(list(width = 1, n_clusters = 30), nan_args)),
    function() do.call(ns$.find_min_width_crd_diff, c(list(budget = 3000, cluster_cost = 10, individual_cost = 1), nan_args)),
    function() do.call(ns$.find_min_cost_crd_diff, c(list(width = 1, cluster_cost = 10, individual_cost = 1), nan_args)))
  for (i in seq_along(checks)) {
    err <- suppressWarnings(.crd_err(checks[[i]]()))
    expect_s3_class(err, "dmar_computation_failed")
    expect_identical(err$reason, "width_not_a_number", info = i)
  }
  expect_match(conditionMessage(suppressWarnings(.crd_err(checks[[1]]()))), "^planner_x\\(\\): the width")
})

test_that("a target no design of up to a million clusters meets is refused", {
  skip_on_cran()
  # A width of 1e-300 was searched for without end (DMAR 1.0.0 and the
  # candidate), and one of .003 planned 1,639,025 clusters of 20.
  for (w in c(1e-300, 0.003)) {
    err <- .crd_err(.crd_call(ss_aipe_crd_n_clusters_fixed_width, width = w, cluster_size = 20,
                              pr_treat = 0.5))
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$reason, "n_clusters_beyond_reach")
    expect_identical(err$max_n_clusters, 1e6)
    expect_gt(err$width_at_max, w)
  }
  err <- .crd_err(.crd_call(ss_aipe_crd_both_fixed_width, width = 0.003, cluster_cost = 10,
                            individual_cost = 1, pr_treat = 0.5))
  expect_identical(err$reason, "n_clusters_beyond_reach")
  expect_match(conditionMessage(err), "clusters of 100,000 individuals", fixed = TRUE)
  # A target met below the bound is planned as before.
  res <- .crd_call(ss_aipe_crd_n_clusters_fixed_width, width = 0.004, cluster_size = 20, pr_treat = 0.5)
  expect_identical(.crd_row(res, "necessary_n_clusters"), 921953)
})

test_that("diff_size is NULL or a vector of finite numbers", {
  skip_on_cran()
  # An empty vector or one with Inf searched without end in the first
  # build of this repair; DMAR 1.0.0 and the candidate stopped with
  # "missing value where TRUE/FALSE needed".
  for (bad in list(numeric(0), c(Inf, 0), NA, NA_real_, "a", list(1))) {
    for (fn in list(
      function(d) .crd_call(ss_aipe_crd_n_clusters_fixed_width, width = 1, cluster_size = 20, pr_treat = 0.5, diff_size = d),
      function(d) .crd_call(ss_aipe_crd_both_fixed_budget, budget = 3000, cluster_cost = 10, individual_cost = 1, pr_treat = 0.5, diff_size = d),
      function(d) ss_aipe_crd_n_clusters_fixed_budget(budget = 1000, cluster_size = 10, diff_size = d))) {
      err <- .crd_err(fn(bad))
      expect_s3_class(err, "dmar_argument_invalid")
      expect_identical(err$arg, "diff_size")
    }
  }
})

test_that("the searches over both sizes name the planner when no candidate is found", {
  skip_on_cran()
  ns <- asNamespace("DMAR")
  # The public planners refuse such a budget first (budget_infeasible); the
  # internal search refuses it too.
  err <- .crd_err(ns$.find_min_width_crd_diff(budget = 10, cluster_cost = 10, individual_cost = 1,
                                             pr_treat = 0.5, true_var_between = 0.2,
                                             true_var_within = 0.8, fn = "ss_aipe_crd_both_fixed_budget"))
  expect_s3_class(err, "dmar_computation_failed")
  expect_identical(err$reason, "no_candidate")
  expect_identical(err$floor, 4)
  expect_match(conditionMessage(err), "^ss_aipe_crd_both_fixed_budget\\(\\): the search")
  expect_match(conditionMessage(err), "Give a larger budget, lower costs, or a pr_treat nearer 0.5",
               fixed = TRUE)
  # The floor as the refusals state it, with covariates.
  for (k in 1:2) {
    err <- .crd_err(.crd_call(ss_aipe_crd_both_fixed_budget, budget = 10, cluster_cost = 10,
                              individual_cost = 1, pr_treat = 0.5, n_covariates = k))
    expect_identical(err$n_clusters, 4 + k)
    expect_match(conditionMessage(err),
                 if (k == 1) "plus one for the covariate" else "plus one for each of the 2 covariates",
                 fixed = TRUE)
  }
})

test_that("the searches over the number of clusters reach past the floor and past 1100 clusters", {
  ns <- asNamespace("DMAR")
  # At pr_treat = .005 the floor is 400. A target first met at 200 (one
  # treated cluster) was planned there; the search goes on from the floor.
  w <- ns$.find_width_crd_diff(200, 20, 0.005, true_var_between = 0.2, true_var_within = 0.8) * (1 + 1e-9)
  expect_equal(ns$.find_n_clus_crd_diff(width = w, cluster_size = 20, pr_treat = 0.005,
                                        true_var_between = 0.2, true_var_within = 0.8), 400)
  # Every number from 100 to 1000 meets the target and the floor, above
  # 100 and off the grid, does not: the first from the floor that does.
  lo <- ns$.crd_min_n_clusters(0.008)
  expect_gt(lo %% 100, 0)
  local_mocked_bindings(.find_width_crd_diff = function(n_clusters, ...) if (n_clusters == lo) 2 else 0.5)
  expect_equal(ns$.find_n_clus_crd_diff(width = 1, cluster_size = 20, pr_treat = 0.008,
                                        true_var_between = 0.2, true_var_within = 0.8), lo + 1)
})

test_that(".crd_first_meeting() looks past its first block", {
  meets <- function(n) n >= 250
  expect_identical(asNamespace("DMAR")$.crd_first_meeting(5, meets), 250L)
})

test_that("the searches over both sizes plan past 1100 clusters as before", {
  skip_on_cran()
  # Plans of 2,000 to 17,864 clusters, from the blocks above 1100; the
  # candidate gave the same.
  res <- .crd_call(ss_aipe_crd_both_fixed_budget, budget = 30000, cluster_cost = 10, individual_cost = 1,
                   pr_treat = 0.5)
  expect_identical(res$value[1:2], c(2000, 5))
  expect_identical(.crd_row(res, "budget"), 30000)
  # Past a second block of 1000.
  res <- .crd_call(ss_aipe_crd_both_fixed_budget, budget = 60000, cluster_cost = 10, individual_cost = 1,
                   pr_treat = 0.5)
  expect_identical(res$value[1:2], c(3529, 7))
  expect_identical(.crd_row(res, "budget"), 59993)
  res <- .crd_call(ss_aipe_crd_both_fixed_width, width = 0.05, cluster_cost = 10, individual_cost = 1,
                   pr_treat = 0.5)
  expect_identical(res$value[1:2], c(6887, 10))
  expect_identical(.crd_row(res, "budget"), 137740)
})

test_that("the variance resolver refuses the remaining contradictions and keeps a difference within the tolerance", {
  resolve <- getFromNamespace(".crd_resolve_variances", "DMAR")
  rs <- function(...) tryCatch(resolve(...), error = function(e) e)
  # A between-cluster variance with an intraclass correlation of 0.
  err <- rs(true_var_between = 0.2, true_icc = 0)
  expect_identical(err$reason, "inconsistent_variances")
  expect_identical(err$arg, c("true_var_between", "true_icc"))
  # No within-cluster variance with an intraclass correlation of 1 leaves
  # the scale open.
  err <- rs(true_var_within = 0, true_icc = 1)
  expect_identical(err$reason, "underidentified_variances")
  # A total a hair below the between-cluster variance gives a within
  # variance of 0.
  got <- resolve(true_total_var = 1, true_var_between = 1 + 2 * .Machine$double.eps)
  expect_identical(got$true_var_within, 0)
  # Every refusal says what to give.
  for (e in list(rs(true_total_var = 1, true_icc = NA), rs(true_var_between = -1, true_var_within = 1),
                 rs(true_total_var = 0, true_icc = 0.2), rs(true_total_var = 1, true_icc = 2))) {
    expect_match(conditionMessage(e), "Give ", fixed = TRUE)
  }
})

test_that("a budget below the smallest design with a width names the budget that buys one", {
  # The budget of 60 was refused for not buying 4 clusters (80), and 80 then
  # bought 4 clusters with none treated.
  err <- .crd_err(.crd_call(ss_aipe_crd_n_clusters_fixed_budget, budget = 60, cluster_size = 20,
                            cluster_cost = 0, individual_cost = 1, pr_treat = 0.1))
  expect_identical(err$reason, "budget_infeasible")
  expect_identical(c(err$n_clusters, err$min_cost), c(6, 120))
  expect_match(conditionMessage(err), "the fewest that put a cluster in each condition", fixed = TRUE)
  expect_s3_class(.crd_call(ss_aipe_crd_n_clusters_fixed_budget, budget = 120, cluster_size = 20,
                            cluster_cost = 0, individual_cost = 1, pr_treat = 0.1), "dmar_tbl")
  err <- .crd_err(.crd_call(ss_aipe_crd_n_clusters_fixed_budget, budget = 60, cluster_size = 20,
                            cluster_cost = 0, individual_cost = 1, pr_treat = 0.5, n_covariates = 2))
  expect_identical(c(err$n_clusters, err$min_cost), c(5, 100))
  # Without pr_treat no width is computed, and 4 clusters are the smallest.
  err <- .crd_err(ss_aipe_crd_n_clusters_fixed_budget(budget = 60, cluster_size = 20, cluster_cost = 0,
                                                      individual_cost = 1))
  expect_identical(c(err$n_clusters, err$min_cost), c(4, 80))
  expect_match(conditionMessage(err), "4 clusters of 20 individuals, which costs 80", fixed = TRUE)
})

test_that("the cost and budget refusals say what to give", {
  for (e in list(.crd_err(ss_aipe_crd_n_clusters_fixed_budget(budget = -5, cluster_size = 20)),
                 .crd_err(ss_aipe_crd_n_clusters_fixed_budget(budget = 100, cluster_size = 20, cluster_cost = -1)),
                 .crd_err(ss_aipe_crd_n_clusters_fixed_budget(budget = 100, cluster_size = 20,
                                                             individual_cost = c(1, 2))))) {
    expect_s3_class(e, "dmar_argument_invalid")
    expect_match(conditionMessage(e), "Give one ", fixed = TRUE)
  }
})

test_that("a search over both sizes without a candidate states the target width", {
  # The public planner refuses such a target first (n_clusters_beyond_reach);
  # the refusal of the search itself names the planner and the target.
  err <- .crd_err(asNamespace("DMAR")$.crd_refuse_no_candidate("ss_aipe_crd_both_fixed_width", 4, 0.5, 0,
                                                               c(4, 200), width = 1))
  expect_s3_class(err, "dmar_computation_failed")
  expect_identical(err$reason, "no_candidate")
  expect_match(conditionMessage(err), "reaches the target width of 1 with clusters of at most 100,000",
               fixed = TRUE)
  expect_match(conditionMessage(err), "Give a wider target or a pr_treat nearer 0.5", fixed = TRUE)
})
