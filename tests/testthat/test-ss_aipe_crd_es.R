## Sample-size planning for cluster-randomized designs with a standardized
## effect size target. These functions use a priori Monte Carlo simulation
## (G replications inside an iterative search) and are slow, so we cap
## G at a small value and skip the heaviest runs on CRAN.

test_that("ss_aipe_crd_es_n_clusters_fixed_width() returns a tidy table that reaches the target width", {
  skip_on_cran()
  # One planner call carries both the shape assertions and the width
  # assertion: an identical second call used to double the runtime of the
  # slowest file in the suite for no additional coverage.
  target <- 0.40
  res <- suppressMessages(
    ss_aipe_crd_es_n_clusters_fixed_width(width = target, cluster_size = 20, true_es = 0.5,
                                    true_icc = 0.10, pr_treat = 0.5,
                                    G = 12, seed = 113)
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_n_clusters" %in% res$term)
  expect_true("expected_width" %in% res$term)
  expect_lte(res$value[res$term == "expected_width"], target + 0.05)
})

test_that("ss_aipe_crd_es_n_clusters_fixed_budget() returns a tidy data.frame within budget", {
  skip_on_cran()
  res <- suppressMessages(
    ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
                                     cluster_cost = 10, individual_cost = 1,
                                     pr_treat = 0.5, true_icc = 0.10, true_es = 0.5,
                                     G = 12, seed = 113)
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_n_clusters" %in% res$term)
  expect_lte(res$value[res$term == "budget"], 1000)
})

test_that("ss_aipe_crd_es_*() produces deterministic output for a fixed seed", {
  skip_on_cran()
  res1 <- suppressMessages(
    ss_aipe_crd_es_n_clusters_fixed_width(width = 0.40, cluster_size = 20, true_es = 0.5,
                                    true_icc = 0.10, pr_treat = 0.5,
                                    G = 6, seed = 113)
  )
  res2 <- suppressMessages(
    ss_aipe_crd_es_n_clusters_fixed_width(width = 0.40, cluster_size = 20, true_es = 0.5,
                                    true_icc = 0.10, pr_treat = 0.5,
                                    G = 6, seed = 113)
  )
  expect_equal(res1, res2)
})

test_that("ss_aipe_crd_es_n_individuals_fixed_budget() returns a tidy data.frame within budget", {
  skip_on_cran()
  res <- suppressMessages(
    ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 30,
                                      cluster_cost = 10, individual_cost = 1,
                                      pr_treat = 0.5, true_icc = 0.10, true_es = 0.5,
                                      G = 12, seed = 113)
  )
  expect_s3_class(res, "data.frame")
  expect_true("necessary_cluster_size" %in% res$term)
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
    ss_aipe_crd_es_both_fixed_budget(budget = 1000, cluster_cost = 10,
                                    individual_cost = 1, true_es = 0.5, true_icc = 0.10,
                                    pr_treat = 0.5, G = 12, seed = 113)
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("necessary_n_clusters", "necessary_cluster_size", "budget") %in% res$term))
  expect_lte(res$value[res$term == "budget"], 1000)
})

test_that("ss_aipe_crd_es_n_individuals_fixed_width() returns a tidy cluster size for a fixed number of clusters", {
  skip_on_cran()
  res <- suppressMessages(
    ss_aipe_crd_es_n_individuals_fixed_width(width = 0.60, n_clusters = 60, true_es = 0.5,
                                      true_icc = 0.10, pr_treat = 0.5,
                                      G = 12, seed = 113)
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_cluster_size" %in% res$term)
  expect_gte(res$value[res$term == "necessary_cluster_size"], 1)
  expect_true("expected_width" %in% res$term)
})

test_that("ss_aipe_crd_es_both_fixed_width() reports a feasible n_clusters/cluster_size combination", {
  skip_on_cran()
  res <- suppressMessages(
    ss_aipe_crd_es_both_fixed_width(width = 0.60, cluster_cost = 10, individual_cost = 1,
                                    true_es = 0.5, true_icc = 0.10, pr_treat = 0.5,
                                    G = 12, seed = 113)
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("necessary_n_clusters", "necessary_cluster_size") %in% res$term))
  expect_gte(res$value[res$term == "necessary_n_clusters"], 1)
  expect_gte(res$value[res$term == "necessary_cluster_size"], 1)
})

test_that("the six ss_aipe_crd_es_*() planners refuse a G that is not a count", {
  skip_on_cran()
  # The check is made at entry, before any simulation, so these run at once.
  bad <- list(fractional = 12.5, zero = 0, negative = -3, missing = NA,
              character = "12", two = c(10, 20), infinite = Inf)
  planners <- list(
    ss_aipe_crd_es_n_clusters_fixed_width = function(G)
      ss_aipe_crd_es_n_clusters_fixed_width(width = 0.4, cluster_size = 20,
        true_es = 0.5, true_icc = 0.1, pr_treat = 0.5, G = G),
    ss_aipe_crd_es_n_individuals_fixed_width = function(G)
      ss_aipe_crd_es_n_individuals_fixed_width(width = 0.4, n_clusters = 40,
        true_es = 0.5, true_icc = 0.1, pr_treat = 0.5, G = G),
    ss_aipe_crd_es_n_clusters_fixed_budget = function(G)
      ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
        cluster_cost = 10, individual_cost = 1, pr_treat = 0.5, true_icc = 0.1,
        true_es = 0.5, G = G),
    ss_aipe_crd_es_n_individuals_fixed_budget = function(G)
      ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 30,
        cluster_cost = 10, individual_cost = 1, pr_treat = 0.5, true_icc = 0.1,
        true_es = 0.5, G = G),
    ss_aipe_crd_es_both_fixed_budget = function(G)
      ss_aipe_crd_es_both_fixed_budget(budget = 1000, cluster_cost = 10,
        individual_cost = 1, true_es = 0.5, true_icc = 0.1, pr_treat = 0.5, G = G),
    ss_aipe_crd_es_both_fixed_width = function(G)
      ss_aipe_crd_es_both_fixed_width(width = 0.5, cluster_cost = 5,
        individual_cost = 1, true_es = 0.5, true_icc = 0.25, pr_treat = 0.5, G = G))
  for (p in names(planners)) {
    for (b in names(bad)) {
      err <- tryCatch(planners[[p]](bad[[b]]), error = function(e) e)
      expect_s3_class(err, "dmar_count_invalid")
      expect_s3_class(err, "dmar_refusal")
      expect_match(conditionMessage(err), "'G' must be a single whole number",
                   fixed = TRUE, info = paste(p, b))
      expect_identical(err$arg, "G", info = paste(p, b))
    }
  }
})

test_that("G = NULL still skips the simulation in the two budget planners", {
  a <- ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
                                              cluster_cost = 0, individual_cost = 1)
  expect_equal(a$value[a$term == "necessary_n_clusters"], 50)
  b <- ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 200,
                                                 cluster_cost = 0, individual_cost = 1)
  expect_s3_class(b, "dmar_tbl")
  expect_false(any(c("expected_width", "assured_width") %in% b$term))
})

test_that("the standardized planners name the width row and the planning values as their siblings do", {
  skip_on_cran()
  es <- c("ss_aipe_crd_es_n_clusters_fixed_width", "ss_aipe_crd_es_n_individuals_fixed_width",
          "ss_aipe_crd_es_n_clusters_fixed_budget", "ss_aipe_crd_es_n_individuals_fixed_budget",
          "ss_aipe_crd_es_both_fixed_budget", "ss_aipe_crd_es_both_fixed_width")
  old <- c("num_predictors", "n_individuals", "icc_Y", "icc_Z", "R2_between", "R2_within",
           "clus_cost", "indiv_cost", "nrep")
  for (fn in es) {
    fm <- names(formals(get(fn, asNamespace("DMAR"))))
    expect_false(any(old %in% fm), info = fn)
    expect_true(all(c("n_covariates", "true_icc", "true_icc_covariate", "true_R2_between",
                      "true_R2_within", "cluster_cost", "individual_cost", "G") %in% fm), info = fn)
  }
  expect_error(do.call(ss_aipe_crd_es_n_clusters_fixed_budget,
                       list(budget = 1000, cluster_size = 20, cluster_cost = 0,
                            individual_cost = 1, nrep = 2)),
               "unused argument")
  skip_if_not_installed("OpenMx")
  res <- ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
                                                cluster_cost = 0, individual_cost = 1,
                                                true_es = 0.5, es_type = 1, true_icc = 0.25,
                                                pr_treat = 0.5, G = 4, seed = 113,
                                                assurance = 0.75)
  expect_identical(res$term, c("necessary_n_clusters", "assured_width", "assurance", "budget"))
  expect_identical(res$value[res$term == "assurance"], 0.75)
})

test_that("a width asked for without the effect size is refused by name in all six planners", {
  skip_on_cran()
  # The four planners that take or return a width need true_es always.
  # No simulation starts, so OpenMx is not needed for the refusal.
  always <- list(
    ss_aipe_crd_es_n_clusters_fixed_width = function(...)
      ss_aipe_crd_es_n_clusters_fixed_width(width = 0.4, cluster_size = 20,
        true_icc = 0.1, pr_treat = 0.5, G = 2, ...),
    ss_aipe_crd_es_n_individuals_fixed_width = function(...)
      ss_aipe_crd_es_n_individuals_fixed_width(width = 0.4, n_clusters = 40,
        true_icc = 0.1, pr_treat = 0.5, G = 2, ...),
    ss_aipe_crd_es_both_fixed_budget = function(...)
      ss_aipe_crd_es_both_fixed_budget(budget = 1000, cluster_cost = 10,
        individual_cost = 1, true_icc = 0.1, pr_treat = 0.5, G = 2, ...),
    ss_aipe_crd_es_both_fixed_width = function(...)
      ss_aipe_crd_es_both_fixed_width(width = 0.5, cluster_cost = 5,
        individual_cost = 1, true_icc = 0.25, pr_treat = 0.5, G = 2, ...))
  for (fn in names(always)) {
    err <- tryCatch(always[[fn]](), error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_s3_class(err, "dmar_refusal")
    expect_identical(err$arg, "true_es", info = fn)
    expect_match(conditionMessage(err),
                 paste0(fn, "() needs 'true_es'"), fixed = TRUE, info = fn)
    expect_match(conditionMessage(err), "none was given", fixed = TRUE)
    expect_match(conditionMessage(err), "such as true_es = 0.5", fixed = TRUE)
    for (bad in list(NULL, NA_real_, c(0.2, 0.5), "0.5", Inf)) {
      expect_error(always[[fn]](true_es = bad),
                   class = "dmar_argument_invalid", info = fn)
    }
  }
  # The two planners that fix a size and a budget need it once G asks
  # for the expected width, and say which population values are missing.
  budget <- list(
    ss_aipe_crd_es_n_clusters_fixed_budget = function(...)
      ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 20,
        cluster_cost = 10, individual_cost = 1, ...),
    ss_aipe_crd_es_n_individuals_fixed_budget = function(...)
      ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 30,
        cluster_cost = 10, individual_cost = 1, ...))
  for (fn in names(budget)) {
    err <- tryCatch(budget[[fn]](G = 2, pr_treat = 0.5, true_icc = 0.1),
                    error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "true_es", info = fn)
    expect_match(conditionMessage(err), paste0(
      fn, "() was asked for the expected width of the interval (G = 2 ",
      "replications), but 'true_es' (the population standardized ",
      "difference between the condition means) was not given."),
      fixed = TRUE, info = fn)
    expect_match(conditionMessage(err), "or leave G = NULL", fixed = TRUE)
    # Every missing population value is named, where the table used to
    # come back without the width.
    all3 <- tryCatch(budget[[fn]](G = 2), error = function(e) e)
    expect_s3_class(all3, "dmar_argument_invalid")
    expect_identical(all3$arg, c("pr_treat", "true_icc", "true_es"))
    expect_match(conditionMessage(all3), "were not given", fixed = TRUE)
    expect_error(budget[[fn]](G = 2, pr_treat = 0.5, true_icc = 0.1,
                              true_es = "0.5"),
                 class = "dmar_argument_invalid", info = fn)
    # Without G no width is asked for, and nothing more is needed.
    expect_s3_class(budget[[fn]](), "dmar_tbl")
  }
})

test_that("the two budget planners describe a covariate as the width planners do, and refuse a bad value by class", {
  skip_on_cran()
  args <- list(budget = 1500, cluster_size = 20, cluster_cost = 2, individual_cost = 1,
               G = 3, pr_treat = 0.5, true_icc = 0.2, true_es = 0.4, es_type = 1,
               n_covariates = 1, seed = 7)
  # The defaults: R2 of 0 at both levels, and the outcome's intraclass
  # correlation for the covariate.
  expect_identical(formals(ss_aipe_crd_es_n_clusters_fixed_budget)$true_R2_within, 0)
  expect_identical(formals(ss_aipe_crd_es_n_individuals_fixed_budget)$true_R2_between, 0)
  plain <- do.call(ss_aipe_crd_es_n_clusters_fixed_budget, args)
  stated <- do.call(ss_aipe_crd_es_n_clusters_fixed_budget,
                    c(args, list(true_icc_covariate = 0.2, true_R2_within = 0,
                                 true_R2_between = 0)))
  expect_identical(plain$value, stated$value)
  expect_true(is.finite(plain$value[plain$term == "expected_width"]))
  for (nm in c("true_R2_within", "true_R2_between", "true_icc_covariate")) {
    bad <- args
    bad[[nm]] <- if (nm == "true_icc_covariate") 1.5 else NA_real_
    err <- tryCatch(do.call(ss_aipe_crd_es_n_clusters_fixed_budget, bad), error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, nm)
  }
  ind <- args
  ind$cluster_size <- NULL
  ind$n_clusters <- 100
  err <- tryCatch(do.call(ss_aipe_crd_es_n_individuals_fixed_budget,
                          c(ind, list(true_R2_within = NULL))),
                  error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "true_R2_within")
  two <- args
  two$n_covariates <- 2
  err <- tryCatch(do.call(ss_aipe_crd_es_n_clusters_fixed_budget, two), error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "n_covariates")
})

## The simulated data sets that give no usable interval, the warnings the
## fits give, the socket cluster, the bounds of the searches, and the seed
## (DR30 to DR35, DR47, and DR50 of the November 2026 review). Where the
## behavior under test is the summary of the replications, the replication
## worker is replaced by a fake that hands out known widths, so that the
## expected values follow from arithmetic and no OpenMx fit is run.

# A fake replication worker: the g-th call returns widths[g] (NA with
# reasons[g] for a failed draw) and the warnings in warnings[[g]].
crd_fake_worker <- function(widths, reasons = rep(NA_character_, length(widths)),
                            warnings = NULL, messages = NULL) {
  st <- new.env(parent = emptyenv())
  st$g <- 0L
  function(seed, ...) {
    st$g <- st$g + 1L
    g <- (st$g - 1L) %% length(widths) + 1L
    list(width = widths[g], reason = reasons[g],
         message = if (is.null(messages)) NULL else messages[g],
         warnings = if (is.null(warnings)) character(0) else warnings[[g]])
  }
}

# Collects the warnings an expression gives, muffled, in an environment
# made for the purpose rather than by assignment to an enclosing scope.
crd_collect_warnings <- function(expr) {
  st <- new.env(parent = emptyenv())
  st$seen <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    st$seen[[length(st$seen) + 1L]] <- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = st$seen)
}

crd_budget_call <- function(G, assurance = NULL, ncpus = 1, seed = 11) {
  ss_aipe_crd_es_n_clusters_fixed_budget(
    budget = 1000, cluster_size = 20, cluster_cost = 0, individual_cost = 1,
    true_es = 0.5, es_type = 1, true_icc = 0.25, pr_treat = 0.5, G = G,
    assurance = assurance, seed = seed, ncpus = ncpus)
}

test_that("draws without a usable interval are counted, and the expected width is the mean of the others", {
  skip_on_cran()
  local_mocked_bindings(.runrep_width_es_crd = crd_fake_worker(
    widths = c(0.2, NA, 0.4, NA, 0.6, NA),
    reasons = c(NA, "limit_not_found", NA, "optimizer_status_6", NA, "limit_not_found")))
  w <- expect_warning(res <- crd_budget_call(G = 6),
                      class = "dmar_replications_without_interval")
  # The mean of the three usable widths, by arithmetic.
  expect_equal(res$value[res$term == "expected_width"], (0.2 + 0.4 + 0.6) / 3)
  sim <- attr(res, "simulation")
  expect_identical(c(sim$attempted, sim$returned, sim$failed, sim$denominator),
                   c(6L, 3L, 3L, 3L))
  expect_identical(sim$returned + sim$failed, sim$attempted)
  expect_identical(sim$failure_reasons,
                   c(limit_not_found = 2L, optimizer_status_6 = 1L))
  expect_identical(c(w$attempted, w$returned, w$failed, w$denominator),
                   c(6L, 3L, 3L, 3L))
  expect_match(conditionMessage(w), "3 of the 6 data sets simulated", fixed = TRUE)
  expect_match(conditionMessage(w), "mean over the 3 that gave one", fixed = TRUE)
  expect_true(any(grepl("Replications: 6 attempted, 3 with an interval, 3 without",
                        capture.output(print(res)), fixed = TRUE)))
})

test_that("with an assurance a draw without an interval counts as wider than any width", {
  widths <- c(0.31, 0.12, NA, 0.55, 0.27, 0.48, NA, 0.39, 0.22, 0.60)
  reasons <- ifelse(is.na(widths), "limit_not_found", NA_character_)
  local_mocked_bindings(.runrep_width_es_crd = crd_fake_worker(widths, reasons))
  # The type 7 quantile of the ten widths with the two failed draws at Inf,
  # by hand: sorted 0.12 0.22 0.27 0.31 0.39 0.48 0.55 0.60 Inf Inf; the
  # 0.7 quantile sits at position 1 + 9 * 0.7 = 7.3, so it is
  # 0.7 * 0.55 + 0.3 * 0.60.
  w <- expect_warning(res <- crd_budget_call(G = 10, assurance = 0.7),
                      class = "dmar_replications_without_interval")
  expect_equal(res$value[res$term == "assured_width"], 0.7 * 0.55 + 0.3 * 0.60)
  expect_identical(attr(res, "simulation")$denominator, 10L)
  expect_match(conditionMessage(w), "wider than any width", fixed = TRUE)
  # At 0.8 the position is 8.2, between 0.60 and a failed draw: no width is
  # reached with that probability, and the call is refused with the counts.
  err <- tryCatch(suppressWarnings(crd_budget_call(G = 10, assurance = 0.8)),
                  error = function(e) e)
  expect_s3_class(err, "dmar_computation_failed")
  expect_identical(c(err$attempted, err$returned, err$failed), c(10L, 8L, 2L))
  expect_match(conditionMessage(err), "no width is reached with probability 0.8",
               fixed = TRUE)
})

test_that("a design at which no draw gave an interval is refused with the counts, not NaN", {
  local_mocked_bindings(.runrep_width_es_crd = crd_fake_worker(
    widths = rep(NA_real_, 4), reasons = rep("fit_error", 4)))
  err <- tryCatch(crd_budget_call(G = 4), error = function(e) e)
  expect_s3_class(err, "dmar_computation_failed")
  expect_identical(c(err$attempted, err$returned, err$failed), c(4L, 0L, 4L))
  expect_identical(err$failure_reasons, c(fit_error = 4L))
  expect_match(conditionMessage(err), "none of the 4 data sets", fixed = TRUE)
  # The internal summary gives Inf, which a search reads as a width that
  # does not meet the target.
  f <- getFromNamespace(".find_width_crd_es", "DMAR")
  v <- suppressWarnings(f(G = 4, n_clusters = 20, n_treat_clus = 10, cluster_size = 10,
                          true_icc = 0.2, es = 0.5, seed = 1))
  expect_identical(as.numeric(v), Inf)
  expect_identical(attr(v, "simulation")$failed, 4L)
})

test_that("the search for the number of clusters reaches a design far above its start, as the one-step search would", {
  ns <- asNamespace("DMAR")
  # Width 0.2 + 100.1 / K against 0.5: the smallest K with
  # 0.2 + 100.1 / K < 0.5 is 334 (K > 100.1 / 0.3 = 333.67), by algebra,
  # 322 steps above a start of 12. The first build of this repair refused
  # it after 100 steps (the reviewer's r05c gives a natural case).
  st <- new.env(parent = emptyenv())
  st$n <- 0L
  local_mocked_bindings(
    .find_n_clus_crd_diff = function(...) 12,
    .find_width_crd_es = function(..., n_clusters) {
      st$n <- st$n + 1L
      0.2 + 100.1 / n_clusters
    })
  res <- ns$.find_n_clus_crd_es(width = 0.5, cluster_size = 10, es = 0.5, true_icc = 0.2,
                                pr_treat = 0.5, G = 4)
  expect_identical(res[1], 334)
  expect_equal(res[2], 0.2 + 100.1 / 334)
  # The start, a bracket of doubling steps (13, 15, 19, ..., 523), and a
  # bisection: at most 1 + 2 * ceiling(log2(334 - 12 + 1)) + 1 designs, where
  # the one-step search would take 323.
  expect_lte(st$n, 1 + 2 * ceiling(log2(334 - 12 + 1)) + 1)
})

test_that("an upward search stops after three designs in a row with no finite width, with the counts and the error", {
  st <- new.env(parent = emptyenv())
  st$n <- 0L
  fake <- crd_fake_worker(widths = rep(NA_real_, 3), reasons = rep("fit_error", 3),
                          messages = rep("object 'zz' not found", 3))
  local_mocked_bindings(.runrep_width_es_crd = function(...) {
    st$n <- st$n + 1L
    fake(...)
  })
  err <- tryCatch(ss_aipe_crd_es_n_clusters_fixed_width(
    width = 0.8, cluster_size = 10, true_es = 0.5, es_type = 1, true_icc = 0.2,
    pr_treat = 0.5, G = 3, seed = 1), error = function(e) e)
  expect_s3_class(err, "dmar_computation_failed")
  expect_identical(err$reason, "no_finite_width")
  # Three designs of three data sets each: the start and two designs each
  # twice the size of the one before (the first build of this repair fitted
  # 101 designs before it stopped, and said a larger G could help).
  expect_identical(st$n, 9L)
  expect_identical(err$sizes[2:3], err$sizes[1] * c(2, 4))
  expect_identical(c(err$attempted, err$returned, err$failed), c(9L, 0L, 9L))
  expect_identical(err$failure_reasons, c(fit_error = 9L))
  expect_identical(err$failure_messages, c(fit_error = "object 'zz' not found"))
  expect_identical(err$search$attempted, 9L)
  expect_match(conditionMessage(err), "no simulated data set gave a usable confidence interval",
               fixed = TRUE)
  expect_match(conditionMessage(err), "object 'zz' not found", fixed = TRUE)
  expect_no_match(conditionMessage(err), "larger G", fixed = TRUE)

  # With an assurance that the usable draws cannot reach (two of four
  # failed, more than 1 - 0.8 allows) the message says so and points to a
  # smaller assurance.
  local_mocked_bindings(.runrep_width_es_crd = crd_fake_worker(
    widths = c(0.2, NA, 0.2, NA), reasons = c(NA, "limit_not_found", NA, "limit_not_found")))
  err <- tryCatch(ss_aipe_crd_es_n_clusters_fixed_width(
    width = 0.8, cluster_size = 10, true_es = 0.5, es_type = 1, true_icc = 0.2,
    pr_treat = 0.5, G = 4, seed = 1, assurance = 0.8), error = function(e) e)
  expect_s3_class(err, "dmar_computation_failed")
  expect_identical(err$reason, "no_finite_width")
  expect_identical(c(err$attempted, err$returned, err$failed), c(12L, 6L, 6L))
  expect_match(conditionMessage(err), "A smaller assurance", fixed = TRUE)
})

test_that("the cluster size search stops at its reach and says what the reach, not the target, decided", {
  ns <- asNamespace("DMAR")
  # es_type 0 has no known limit. A width that levels off at 0.6 against a
  # target of 0.5, its Monte Carlo interval 0.01 on each side.
  st <- new.env(parent = emptyenv())
  st$n <- 0L
  local_mocked_bindings(
    .find_n_indiv_crd_diff = function(...) 5,
    .find_width_crd_es = function(..., cluster_size) {
      st$n <- st$n + 1L
      w <- 0.6 + 1 / cluster_size
      structure(w, simulation = list(width_mc_ci = w + c(-0.01, 0.01)))
    })
  err <- tryCatch(ns$.find_n_indiv_crd_es(width = 0.5, n_clusters = 40, es = 0.5, es_type = 0,
                                          true_icc = 0.2, pr_treat = 0.5, G = 4),
                  error = function(e) e)
  expect_s3_class(err, "dmar_computation_failed")
  expect_identical(err$reason, "width_not_reached")
  expect_identical(c(err$from, err$to, err$narrowest_at), c(5, 105, 105))
  expect_equal(err$narrowest, 0.6 + 1 / 105)
  # The start and the sizes 6, 8, 12, 20, 36, 68, and 105 (the reach): 8
  # designs, where the first build of this repair simulated 101.
  expect_identical(st$n, 8L)
  expect_match(conditionMessage(err), "was not determined", fixed = TRUE)
  # The target is outside the Monte Carlo interval of the narrowest width,
  # so a larger G is not suggested.
  expect_no_match(conditionMessage(err), "larger G", fixed = TRUE)
  # When the target lies within that interval, it is.
  local_mocked_bindings(
    .find_n_indiv_crd_diff = function(...) 5,
    .find_width_crd_es = function(..., cluster_size) {
      w <- 0.6 + 1 / cluster_size
      structure(w, simulation = list(width_mc_ci = w + c(-0.2, 0.2)))
    })
  err <- tryCatch(ns$.find_n_indiv_crd_es(width = 0.5, n_clusters = 40, es = 0.5, es_type = 0,
                                          true_icc = 0.2, pr_treat = 0.5, G = 4),
                  error = function(e) e)
  expect_match(conditionMessage(err), "a larger G can also tell", fixed = TRUE)
})

test_that("the width the cluster size search approaches with es_type = 1 agrees with an independent computation", {
  ns <- asNamespace("DMAR")
  # The limit as every cluster grows: K cluster means with a common unknown
  # variance tau and the within-cluster standard deviation known. Computed
  # here by another route: the profile likelihood multiplier by root
  # finding (K log(1 + c^2) = qchisq(0.95, 1)), and E sqrt(chi^2_{K - 2}) by
  # numerical integration.
  for (cs in list(c(K = 10, K1 = 5, icc = 0.2), c(K = 7, K1 = 3, icc = 0.05),
                  c(K = 40, K1 = 20, icc = 0.3))) {
    K <- cs[["K"]]; K1 <- cs[["K1"]]; icc <- cs[["icc"]]
    q <- stats::qchisq(0.95, 1)
    c_lr <- stats::uniroot(function(c) K * log(1 + c^2) - q, c(0, 10), tol = 1e-12)$root
    e_root <- stats::integrate(function(x) sqrt(x) * stats::dchisq(x, K - 2), 0, Inf,
                               rel.tol = 1e-10)$value
    s_xx <- K1 * (K - K1) / K
    oracle <- 2 * c_lr * sqrt(icc / (1 - icc)) * e_root / sqrt(s_xx)
    expect_equal(ns$.crd_es_width_limit(K, K1, icc, 0.95), oracle, tolerance = 1e-7)
    oracle_q <- 2 * c_lr * sqrt(icc / (1 - icc)) * sqrt(stats::qchisq(0.8, K - 2)) / sqrt(s_xx)
    expect_equal(ns$.crd_es_width_limit(K, K1, icc, 0.95, assurance = 0.8), oracle_q,
                 tolerance = 1e-7)
  }
  expect_true(is.na(ns$.crd_es_width_limit(2, 1, 0.2)))
  # And against simulation: the mean width of the likelihood-based interval
  # found by root finding on the profile deviance of simulated cluster means
  # (no closed form), 2000 data sets, within four Monte Carlo standard errors.
  # evidence/crd/verify_after_round2.log has 20000 per design.
  skip_on_cran()
  withr::local_seed(113)
  K <- 10; x <- rep(c(1, 0), each = 5); tau <- 0.2 / 0.8
  w <- vapply(seq_len(2000), function(i) {
    y <- 0.5 * x + stats::rnorm(K, 0, sqrt(tau))
    rss <- function(d) { r <- y - d * x; sum((r - mean(r))^2) }
    dhat <- mean(y[x == 1]) - mean(y[x == 0])
    f <- function(d) K * log(rss(d) / rss(dhat)) - stats::qchisq(0.95, 1)
    stats::uniroot(f, c(dhat, dhat + 20), tol = 1e-10)$root -
      stats::uniroot(f, c(dhat - 20, dhat), tol = 1e-10)$root
  }, numeric(1))
  expect_lt(abs(mean(w) - ns$.crd_es_width_limit(10, 5, 0.2)), 4 * stats::sd(w) / sqrt(2000))
})

test_that("a target at or below that width is refused before anything is simulated, and one above it is searched", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  lim <- getFromNamespace(".crd_es_width_limit", "DMAR")(10, 5, 0.2, 0.95)
  err <- tryCatch(ss_aipe_crd_es_n_individuals_fixed_width(
    width = lim, n_clusters = 10, true_es = 0.5, es_type = 1, true_icc = 0.2, pr_treat = 0.5,
    G = 4, seed = 1), error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$reason, "width_not_reachable")
  expect_equal(err$width_limit, lim)
  expect_match(conditionMessage(err), "no cluster size reaches it", fixed = TRUE)
  # With the assurance the limit is that of the quantile, which is larger.
  lim_a <- getFromNamespace(".crd_es_width_limit", "DMAR")(10, 5, 0.2, 0.95, assurance = 0.8)
  err <- tryCatch(ss_aipe_crd_es_n_individuals_fixed_width(
    width = (lim + lim_a) / 2, n_clusters = 10, true_es = 0.5, es_type = 1, true_icc = 0.2,
    pr_treat = 0.5, G = 4, seed = 1, assurance = 0.8), error = function(e) e)
  expect_identical(err$reason, "width_not_reachable")

  # A target of 1.30 lies above the limit (1.187) but below the
  # unstandardized width at 100,000 individuals per cluster (about 1.46,
  # from the t quantile with 8 degrees of freedom), which refused it
  # through the first builds of 1.1.0. The search now starts from the
  # smallest m with lim * sqrt(1 + 4 / m) < 1.30 (4 = (1 - 0.2) / 0.2),
  # which is m = 20 by algebra (4 / ((1.30 / lim)^2 - 1) = 19.4). The
  # unstandardized search gives Inf for a width no cluster size up to
  # 100,000 reaches (the string "> 100000" through the first builds of
  # 1.1.0).
  ns <- asNamespace("DMAR")
  expect_identical(ns$.find_n_indiv_crd_diff(width = 1.30, n_clusters = 10, pr_treat = 0.5,
                                             true_total_var = 1 / 0.8, true_icc = 0.2),
                   Inf)
  m_start <- floor(4 / ((1.30 / lim)^2 - 1)) + 1
  expect_identical(ns$.crd_es_cluster_size_start(1.30, lim, 0.2), m_start)
  st <- new.env(parent = emptyenv())
  st$sizes <- numeric(0)
  local_mocked_bindings(.find_width_crd_es = function(..., cluster_size) {
    st$sizes <- c(st$sizes, cluster_size)
    lim * sqrt(1 + 4 / cluster_size) + 0.05
  })
  res <- ns$.find_n_indiv_crd_es(width = 1.30, n_clusters = 10, es = 0.5, es_type = 1,
                                 true_icc = 0.2, pr_treat = 0.5, G = 4)
  expect_identical(st$sizes[1], m_start)
  # The smallest m with lim * sqrt(1 + 4 / m) + 0.05 < 1.30, by algebra.
  expect_identical(res[1], floor(4 / (((1.30 - 0.05) / lim)^2 - 1)) + 1)
})

test_that("a call stops after its budget of designs, with the totals", {
  ns <- asNamespace("DMAR")
  tracker <- ns$.crd_tracker("ss_aipe_crd_es_both_fixed_width")
  record <- list(width = 0.7, attempted = 4L, returned = 3L, failed = 1L,
                 failure_reasons = c(limit_not_found = 1L), failure_messages = character(0),
                 n_clusters = 30, cluster_size = 5)
  tracker$records <- rep(list(record), ns$.CRD_ES_MAX_DESIGNS)
  err <- tryCatch(ns$.find_width_crd_es(G = 4, n_clusters = 31, n_treat_clus = 15, cluster_size = 5,
                                        true_icc = 0.2, es = 0.5, seed = 1, tracker = tracker),
                  error = function(e) e)
  expect_s3_class(err, "dmar_computation_failed")
  expect_identical(err$reason, "search_budget")
  expect_identical(err$n_designs, 500L)
  expect_identical(c(err$attempted, err$returned, err$failed), c(2000L, 1500L, 500L))
  expect_identical(err$failure_reasons, c(limit_not_found = 500L))
  expect_match(conditionMessage(err), "not a finding about the target", fixed = TRUE)
})

test_that("the both-size budget planner refuses three designs with no finite width with their counts", {
  skip_on_cran()
  local_mocked_bindings(.runrep_width_es_crd = crd_fake_worker(
    widths = rep(NA_real_, 3), reasons = rep("limit_not_found", 3)))
  err <- tryCatch(ss_aipe_crd_es_both_fixed_budget(
    budget = 1000, cluster_cost = 5, individual_cost = 1, true_es = 0.5, es_type = 1,
    true_icc = 0.25, pr_treat = 0.5, G = 3, seed = 1), error = function(e) e)
  expect_s3_class(err, "dmar_computation_failed")
  expect_identical(err$reason, "no_finite_width")
  expect_identical(c(err$attempted, err$returned, err$failed), c(9L, 0L, 9L))
  expect_identical(err$failure_reasons, c(limit_not_found = 9L))
})

test_that("failed draws at other designs of a search are reported when the design reported has none", {
  # The starting design (37 clusters) has no usable interval in any of its
  # four data sets; every later data set gives a width of 0.3, below the
  # target of 0.8, so the search brackets at 74 and bisects down to 38.
  st <- new.env(parent = emptyenv())
  st$g <- 0L
  local_mocked_bindings(
    .find_n_clus_crd_diff = function(...) 37,
    .runrep_width_es_crd = function(...) {
      st$g <- st$g + 1L
      if (st$g <= 4L) list(width = NA_real_, reason = "optimizer_status_5", warnings = character(0))
      else list(width = 0.3, reason = NA_character_, warnings = character(0))
    })
  got <- crd_collect_warnings(ss_aipe_crd_es_n_clusters_fixed_width(
    width = 0.8, cluster_size = 10, true_es = 0.5, es_type = 1, true_icc = 0.2,
    pr_treat = 0.5, G = 4, seed = 1))
  res <- got$value
  expect_identical(res$value[res$term == "necessary_n_clusters"], 38)
  w <- Filter(function(w) inherits(w, "dmar_replications_without_interval"), got$warnings)
  expect_length(w, 1L)
  expect_identical(w[[1]]$scope, "search")
  expect_identical(w[[1]]$failed, 4L)
  expect_identical(w[[1]]$failure_reasons, c(optimizer_status_5 = 4L))
  sim <- attr(res, "simulation")
  expect_identical(sim$failed, 0L)
  expect_identical(sim$search$failed, 4L)
  expect_identical(sim$search$failure_reasons, c(optimizer_status_5 = 4L))
})

test_that("warnings whose messages differ only in their numbers are given once, with every message kept", {
  local_mocked_bindings(.runrep_width_es_crd = crd_fake_worker(
    widths = c(0.2, 0.3, 0.4, 0.5),
    warnings = list(c(simpleWarning = "step size 1.25e-3 at iteration 12"),
                    c(simpleWarning = "step size 0.5 at iteration 3"),
                    character(0),
                    c(simpleWarning = "step size 2 at iteration 40", my_warning = "other"))))
  got <- crd_collect_warnings(crd_budget_call(G = 4))
  w <- Filter(function(w) inherits(w, "dmar_result_qualified"), got$warnings)
  expect_length(w, 2L)
  step <- Filter(function(w) identical(w$warning_class, "simpleWarning"), w)[[1]]
  expect_identical(step$n_replications, 3L)
  expect_identical(step$messages, c("step size 1.25e-3 at iteration 12",
                                    "step size 0.5 at iteration 3",
                                    "step size 2 at iteration 40"))
  expect_match(conditionMessage(step), "3 messages of this form", fixed = TRUE)
  other <- Filter(function(w) identical(w$warning_class, "my_warning"), w)[[1]]
  expect_identical(other$n_replications, 1L)
  sw <- attr(got$value, "simulation")$search$warnings
  expect_length(sw, 2L)
})

test_that("a fit that stops with an error keeps its message in the record", {
  ns <- asNamespace("DMAR")
  got <- ns$.crd_fit_interval(simpleError("model building failed"))
  expect_identical(got$reason, "fit_error")
  expect_identical(got$message, "model building failed")
  rec <- ns$.crd_mc_record(list(
    list(width = NA_real_, reason = "fit_error", message = "model building failed"),
    list(width = 0.4, reason = NA_character_),
    list(width = NA_real_, reason = "limit_not_found")))
  expect_identical(rec$failure_messages, c(fit_error = "model building failed"))
  expect_identical(rec$failure_reasons, c(fit_error = 1L, limit_not_found = 1L))
})

test_that("a reachable width gives the smallest design that meets it, from below and from above", {
  ns <- asNamespace("DMAR")
  # Width 0.5 + 4.2 / m against a target of 1: the smallest whole m with
  # 0.5 + 4.2 / m < 1 is 9, whether the search starts below or above it.
  # es_type 0 has no closed-form start, so the search starts where the mock
  # of the unstandardized start puts it.
  for (start in c(5, 20)) {
    local_mocked_bindings(
      .find_n_indiv_crd_diff = function(...) start,
      .find_width_crd_es = function(..., cluster_size) 0.5 + 4.2 / cluster_size)
    res <- ns$.find_n_indiv_crd_es(width = 1, n_clusters = 40, es = 0.5, es_type = 0,
                                   true_icc = 0.2, pr_treat = 0.5, G = 4)
    expect_identical(res[1], 9, info = paste("start", start))
    expect_equal(res[2], 0.5 + 4.2 / 9)
  }
  # Width 0.2 + 6.1 / K against 0.5: the smallest K with 0.2 + 6.1 / K < 0.5
  # is 21.
  for (start in c(12, 40)) {
    local_mocked_bindings(
      .find_n_clus_crd_diff = function(...) start,
      .find_width_crd_es = function(..., n_clusters) 0.2 + 6.1 / n_clusters)
    res <- ns$.find_n_clus_crd_es(width = 0.5, cluster_size = 10, es = 0.5, true_icc = 0.2,
                                  pr_treat = 0.5, G = 4)
    expect_identical(res[1], 21, info = paste("start", start))
  }
})

test_that("a width equal to the target meets it in both directions of both searches", {
  ns <- asNamespace("DMAR")
  # Widths that equal the target exactly in floating point at one size:
  # 0.5 + 4 / 8 = 1 and 0.25 + 5 / 20 = 0.5. A width meets the target when
  # it is no larger than the target, so the smallest size that meets it is
  # 8 (and 20) from a start below it, at it, or above it. Through the first
  # builds of 1.1.0 the upward search asked for a width below the target and
  # returned 9 (and 21) from below, while the downward search returned 8
  # (and 20).
  expect_identical(0.5 + 4 / 8, 1)
  expect_identical(0.25 + 5 / 20, 0.5)
  for (start in c(5, 8, 20)) {
    local_mocked_bindings(
      .find_n_indiv_crd_diff = function(...) start,
      .find_width_crd_es = function(..., cluster_size) 0.5 + 4 / cluster_size)
    res <- ns$.find_n_indiv_crd_es(width = 1, n_clusters = 40, es = 0.5, es_type = 0,
                                   true_icc = 0.2, pr_treat = 0.5, G = 4)
    expect_identical(res, c(8, 1), info = paste("start", start))
  }
  for (start in c(12, 20, 40)) {
    local_mocked_bindings(
      .find_n_clus_crd_diff = function(...) start,
      .find_width_crd_es = function(..., n_clusters) 0.25 + 5 / n_clusters)
    res <- ns$.find_n_clus_crd_es(width = 0.5, cluster_size = 10, es = 0.5, true_icc = 0.2,
                                  pr_treat = 0.5, G = 4)
    expect_identical(res, c(20, 0.5), info = paste("start", start))
  }
})

test_that("the interval of a fit is read by the column names, and its status is checked", {
  skip_on_cran()
  ns <- asNamespace("DMAR")
  tab <- matrix(c(-0.1, 0.3, 0.9), 1, dimnames = list("es", c("lbound", "estimate", "ubound")))
  row <- ns$.crd_interval_row(tab)
  expect_identical(row, c(-0.1, 0.3, 0.9))
  ok <- ns$.crd_interval_status(status = 0, ci = row, codes = c(0, 0),
                                diagnostic = c("success", "success"))
  # The width is the upper limit minus the lower (1.0), not the estimate
  # minus the lower limit (0.4), which is what the first builds of 1.1.0
  # and DMAR 1.0.x used.
  expect_equal(ok$width, 0.9 - (-0.1))
  expect_true(is.na(ok$reason))
  expect_true(is.na(ns$.crd_interval_status(1, row, c(1, 0))$reason))
  cases <- list(
    list(args = list(5, row, c(0, 0)), reason = "optimizer_status_5"),
    list(args = list(6, row, c(0, 0)), reason = "optimizer_status_6"),
    list(args = list(4, row, c(0, 0)), reason = "optimizer_status_4"),
    list(args = list(NA_real_, row, c(0, 0)), reason = "optimizer_status_NA"),
    list(args = list(0, row, c(0, 3)), reason = "limit_status_3"),
    list(args = list(0, row, c(NA, 0)), reason = "limit_not_found"),
    list(args = list(0, row, c(0, 0), c("success", "alpha level not reached")),
         reason = "limit_not_found"),
    list(args = list(0, NULL, c(0, 0)), reason = "limit_not_found"),
    list(args = list(0, c(-0.1, 0.3, NA), c(0, 0)), reason = "limit_not_finite"),
    list(args = list(0, c(0.5, 0.3, 0.9), c(0, 0)), reason = "limits_not_ordered"),
    list(args = list(0, c(0.9, 0.9, 0.9), c(0, 0)), reason = "limits_not_ordered"))
  for (cs in cases) {
    got <- do.call(ns$.crd_interval_status, cs$args)
    expect_identical(got$reason, cs$reason)
    expect_true(is.na(got$width))
  }
  expect_null(ns$.crd_interval_row(matrix(1:3, 1)))
})

test_that("a natural fit with an unacceptable optimizer status is a failed draw, and a usable fit agrees with the profile likelihood", {
  skip_on_cran()
  skip_if_not_installed("OpenMx")
  ns <- asNamespace("DMAR")
  # Six clusters of two to four individuals, intraclass correlation 0.05:
  # with OpenMx 2.22.11 this data set's fit ends with optimizer status 5
  # (a Hessian that is not convex), with the between-cluster variance
  # estimated below 0, and the earlier code took a width from it.
  bad <- ns$.runrep_width_es_crd(seed = 13, n_clusters = 6, n_treat_clus = 3, cluster_size = 3,
                                 true_icc = 0.05, es = 0.5, es_type = 1,
                                 true_total_var = 1 / 0.95, diff_size = c(-1, 0, 1))
  expect_match(bad$reason, "^optimizer_status_")
  expect_false(bad$reason %in% c("optimizer_status_0", "optimizer_status_1"))
  expect_true(is.na(bad$width))
  # OpenMx's warning of the status is muffled, and no other warning is given
  # (the deprecated objective functions are no longer used).
  expect_length(bad$warnings, 0L)

  # An ordinary data set: the limits against an independent profile
  # likelihood interval for the effect standardized by the within-cluster
  # standard deviation, computed from the closed form of the
  # random-intercept likelihood (no OpenMx).
  set.seed(2)
  K <- 20; m <- 5
  dw <- ns$.create_data_crd_wide(K, K / 2, m, 0.2, 0.5, 1, 1 / 0.8)
  fit <- ns$.lik_ci_es_crd(dw, paste0("y", 1:m), "x", es_type = 1, true_icc = 0.2,
                           es = 0.5, true_total_var = 1 / 0.8)
  expect_true(is.na(fit$reason))
  y <- as.matrix(dw[, paste0("y", 1:m)]); x <- dw$x
  ybar <- rowMeans(y); ssw <- sum((y - ybar)^2)
  dev <- function(mu, gamma, s2, lambda) K * (m - 1) * log(s2) + K * log(lambda) +
    ssw / s2 + m * sum((ybar - mu - gamma * x)^2) / lambda
  profile <- function(delta) stats::optimize(function(ls) {
    s2 <- exp(2 * ls); gamma <- delta * sqrt(s2); mu <- mean(ybar - gamma * x)
    dev(mu, gamma, s2, m * sum((ybar - mu - gamma * x)^2) / K)
  }, c(-5, 5), tol = 1e-12)$objective
  b <- stats::lm.fit(cbind(1, x), ybar)
  s2_hat <- ssw / (K * (m - 1))
  dmin <- dev(b$coefficients[1], b$coefficients[2], s2_hat, m * sum(b$residuals^2) / K)
  delta_hat <- b$coefficients[[2]] / sqrt(s2_hat)
  f <- function(d) profile(d) - dmin - stats::qchisq(0.95, 1)
  lower <- stats::uniroot(f, c(delta_hat - 3, delta_hat), tol = 1e-10)$root
  upper <- stats::uniroot(f, c(delta_hat, delta_hat + 3), tol = 1e-10)$root
  expect_equal(fit$estimate, delta_hat, tolerance = 1e-4)
  # OpenMx's limit search stops a little past the exact limit (the deviance
  # there has risen by up to about 3.89 rather than 3.84), so the widths
  # agree to about one percent.
  expect_equal(fit$width, upper - lower, tolerance = 0.01)
  expect_equal(c(fit$lower, fit$upper), c(lower, upper), tolerance = 0.02)
})

test_that("a warning from a fit reaches the user once, with the number of data sets that gave it", {
  skip_if_not_installed("OpenMx")
  local_mocked_bindings(.lik_ci_es_crd = function(...) {
    warning("injected fit warning")
    list(width = 0.3, reason = NA_character_, lower = 0.1, estimate = 0.25, upper = 0.4)
  })
  got <- crd_collect_warnings(crd_budget_call(G = 5))
  res <- got$value
  seen <- got$warnings
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "dmar_result_qualified")
  expect_identical(seen[[1]]$warning_class, "simpleWarning")
  expect_identical(seen[[1]]$n_replications, 5L)
  expect_match(conditionMessage(seen[[1]]), "5 of the 5 simulated data sets", fixed = TRUE)
  expect_match(conditionMessage(seen[[1]]), "injected fit warning", fixed = TRUE)
  expect_equal(res$value[res$term == "expected_width"], 0.3)
  # Forked workers do not pass a warning back; the worker records it, so the
  # same single warning is given with two processes.
  skip_on_cran()
  skip_on_os("windows")
  seen <- crd_collect_warnings(crd_budget_call(G = 5, ncpus = 2))$warnings
  expect_length(seen, 1L)
  expect_identical(seen[[1]]$n_replications, 5L)
})

test_that("a socket cluster is stopped when a worker stops with an error, and a clean run still works", {
  skip_on_cran()
  ns <- asNamespace("DMAR")
  keep <- new.env(parent = emptyenv())
  real <- parallel::makeCluster
  local_mocked_bindings(makeCluster = function(...) {
    keep$cl <- real(...)
    keep$cl
  }, .package = "parallel")
  # The worker functions have the global environment, so that a socket
  # worker does not need DMAR to unserialize them.
  bad <- function(s) if (s == 3) stop("injected worker error") else s
  environment(bad) <- globalenv()
  expect_error(ns$.crd_apply_replications(as.list(1:4), bad, ncpus = 2L, cluster_type = "PSOCK"),
               "injected worker error")
  still_open <- vapply(keep$cl, function(node)
    tryCatch(isOpen(node$con), error = function(e) FALSE), logical(1))
  expect_false(any(still_open))
  good <- function(s, add) s + add
  environment(good) <- globalenv()
  out <- ns$.crd_apply_replications(as.list(1:4), good, ncpus = 2L, add = 10,
                                    cluster_type = "PSOCK")
  expect_identical(unlist(out), c(11, 12, 13, 14))
  still_open <- vapply(keep$cl, function(node)
    tryCatch(isOpen(node$con), error = function(e) FALSE), logical(1))
  expect_false(any(still_open))
})

test_that("seed = NULL advances the random number stream, and a given seed restores it", {
  local_mocked_bindings(.runrep_width_es_crd = crd_fake_worker(widths = c(0.2, 0.3)))
  set.seed(239)
  before <- .Random.seed
  crd_budget_call(G = 4, seed = NULL)
  expect_false(identical(before, .Random.seed))
  set.seed(239)
  before <- .Random.seed
  crd_budget_call(G = 4, seed = 113)
  expect_identical(before, .Random.seed)
})

test_that("a call that simulates, and one that simulates nothing, give no warning", {
  local_mocked_bindings(.runrep_width_es_crd = crd_fake_worker(widths = c(0.2, 0.3)))
  expect_no_warning(crd_budget_call(G = 4))
  expect_no_warning(ss_aipe_crd_es_n_clusters_fixed_budget(
    budget = 1000, cluster_size = 20, cluster_cost = 0, individual_cost = 1))
})

test_that("a real call gives no warning when every data set has an interval", {
  skip_on_cran()
  skip_if_not_installed("OpenMx")
  # Through the first builds of 1.1.0 every fit gave a deprecation warning
  # from OpenMx, which the planners discarded with everything else.
  expect_no_warning(res <- crd_budget_call(G = 3, seed = 113))
  expect_identical(attr(res, "simulation")$failed, 0L)
  expect_true(is.finite(res$value[res$term == "expected_width"]))
})

test_that("a cluster size search that would start above the largest simulated size is refused before anything is simulated", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  lim <- getFromNamespace(".crd_es_width_limit", "DMAR")(10, 5, 0.2, 0.95)
  start_fn <- getFromNamespace(".crd_es_cluster_size_start", "DMAR")
  cap <- getFromNamespace(".CRD_ES_MAX_CLUSTER_SIZE", "DMAR")
  # Just above the limit the start is far above the cap: with 10 clusters
  # and an intraclass correlation of .2, (1 - .2) / .2 = 4, and the start is
  # floor(4 / ((1.19 / lim)^2 - 1)) + 1.
  target <- 1.19
  expect_gt(target, lim)
  start <- start_fn(target, lim, 0.2)
  expect_identical(start, floor(4 / ((target / lim)^2 - 1)) + 1)
  expect_gt(start, cap)
  err <- tryCatch(ss_aipe_crd_es_n_individuals_fixed_width(
    width = target, n_clusters = 10, true_es = 0.5, es_type = 1, true_icc = 0.2,
    pr_treat = 0.5, G = 4, seed = 1), error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$reason, "cluster_size_beyond_reach")
  expect_identical(err$start, start)
  expect_identical(err$max_cluster_size, cap)
  expect_equal(err$width_limit, lim)
  expect_match(conditionMessage(err), "ss_aipe_crd_es_n_clusters_fixed_width()", fixed = TRUE)
  # With no covariate, equal cluster sizes, and the expected width, the
  # start is a lower bound on the size the target needs, and the refusal
  # says so.
  expect_true(err$lower_bound)
  expect_match(conditionMessage(err), "would start at 711 individuals per cluster, a lower bound",
               fixed = TRUE)
  expect_match(conditionMessage(err), "This search simulates cluster sizes up to 500", fixed = TRUE)
})

test_that("a start beyond the largest simulated size is called a lower bound only where it is one", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  ns <- asNamespace("DMAR")
  refuse <- function(...) tryCatch(ss_aipe_crd_es_n_individuals_fixed_width(
    n_clusters = 10, true_es = 0.5, true_icc = 0.2, pr_treat = 0.5, G = 4, seed = 1, ...),
    error = function(e) e)
  check <- function(err, start, info) {
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$reason, "cluster_size_beyond_reach", info = info)
    expect_identical(err$start, start, info = info)
    expect_false(err$lower_bound, info = info)
    expect_no_match(conditionMessage(err), "lower bound", fixed = TRUE)
  }
  # With an assurance: just above the limit of the 0.8 quantile, the
  # closed-form start, from the quantile's limit, is far above 500.
  lim_a <- ns$.crd_es_width_limit(10, 5, 0.2, 0.95, assurance = 0.8)
  target <- lim_a * 1.0005
  err <- refuse(width = target, es_type = 1, assurance = 0.8)
  check(err, floor(4 / ((target / lim_a)^2 - 1)) + 1, "assurance")
  expect_match(conditionMessage(err), "(the assured width falls toward 1.438", fixed = TRUE)
  # A covariate that varies within clusters only (the limit is known, and
  # the closed-form start uses the variance it leaves, 0.7 of 4).
  lim <- ns$.crd_es_width_limit(10, 5, 0.2, 0.95)
  target <- lim * 1.0005
  check(refuse(width = target, es_type = 1, n_covariates = 1, true_icc_covariate = 0,
               true_R2_within = 0.3),
        floor(0.7 * 4 / ((target / lim)^2 - 1)) + 1, "within covariate")
  # Unequal cluster sizes.
  check(refuse(width = target, es_type = 1, diff_size = c(-1, 0, 1)),
        floor(4 / ((target / lim)^2 - 1)) + 1, "unequal sizes")
  # es_type 0 has no closed-form start: the unstandardized start, 1000 by
  # construction (its width with the t quantile on 8 degrees of freedom and
  # a total variance of 1 falls just below the target there), is refused
  # without a limit.
  t_width <- function(m) 2 * stats::qt(0.975, 8) * sqrt((0.2 + 0.8 / m) / 2.5)
  target <- (t_width(999) + t_width(1000)) / 2
  err <- refuse(width = target, es_type = 0)
  check(err, 1000, "es_type 0")
  expect_true(is.na(err$width_limit))
  expect_no_match(conditionMessage(err), "falls toward", fixed = TRUE)
})

# The width of every simulated data set at cluster size m, for a mock of
# .runrep_width_es_crd(): w(m) = L sqrt(1 + r / m) (1 + c / (K m)), the
# exact expected width with the within-cluster standard deviation known
# (L the limit, r = (1 - icc) / icc), widened by a term that falls with m,
# as estimating the standard deviation widens it (c of about 2.5 at a
# standardized effect of 0.5 with 10 clusters, and some hundreds at an
# effect of 3; evidence/j4_crd_start/02_oracle_band.log and
# 03_reach_risk.log of the November 2026 review). The cluster size of
# every simulated data set is recorded in `st`.
crd_mock_width <- function(L, r, K, c, st) {
  function(seed, n_clusters, n_treat_clus, cluster_size, ...) {
    st$sizes <- c(st$sizes, cluster_size)
    list(width = L * sqrt(1 + r / cluster_size) * (1 + c / (K * cluster_size)),
         reason = NA_character_, warnings = character(0))
  }
}

test_that("es_type 1 starts from the smaller start and plans the smallest size that meets the target", {
  skip_on_cran()  # numerical oracle check, run in the full local suite
  # Pinned to algebra: the limit from its formula, the closed-form start,
  # the unstandardized start by enumeration of its width with the t
  # quantile, and the size the target needs by enumeration of the mocked
  # width. The cases put the unstandardized start below 500, between 501
  # and 100,000 (where the planner refused through the first builds of
  # 1.1.0 although a small cluster size reaches the target), and beyond
  # 100,000. A search from the unstandardized start, which stepped down to
  # the same size, gives the same plan wherever it was at most 500
  # (evidence/j4_crd_start/05_grid_base_vs_new.log has the comparison).
  q <- stats::qchisq(0.95, 1)
  for (K in c(6, 10, 40)) for (icc in c(0.05, 0.2)) for (cc in c(0, 25)) {
    K1 <- K / 2; s_xx <- K1 * (K - K1) / K; r <- (1 - icc) / icc; tau <- icc / (1 - icc)
    e_root <- sqrt(2) * exp(lgamma((K - 1) / 2) - lgamma((K - 2) / 2))
    L <- 2 * sqrt((exp(q / K) - 1) / s_xx) * sqrt(tau) * e_root
    t_width <- function(m) 2 * stats::qt(0.975, K - 2) * sqrt((tau + 1 / m) / s_xx)
    ratio <- t_width(1) / (L * sqrt(1 + r))
    m <- as.numeric(2:100000)
    for (target in L * c(1.05, 1.3, ratio * sqrt(1 + r / c(300, 2000)) * 1.000001)) {
      info <- sprintf("K %d, icc %.2f, c %d, target %.4f", K, icc, cc, target)
      t_ok <- which(t_width(m) < target)
      t_start <- if (length(t_ok)) m[t_ok[1]] else Inf
      cf <- max(2, floor(r / ((target / L)^2 - 1)) + 1)
      start <- min(cf, t_start)
      expect_identical(start, cf, info = info)
      w <- L * sqrt(1 + r / m) * (1 + cc / (K * m))
      m_star <- m[which(w < target)[1]]
      st <- new.env(parent = emptyenv())
      st$sizes <- numeric(0)
      local_mocked_bindings(.runrep_width_es_crd = crd_mock_width(L, r, K, cc, st))
      res <- tryCatch(ss_aipe_crd_es_n_individuals_fixed_width(
        width = target, n_clusters = K, true_es = 0.5, es_type = 1, true_icc = icc,
        pr_treat = 0.5, G = 2, seed = 1), error = function(e) e)
      if (start > 500) {
        expect_identical(res$reason, "cluster_size_beyond_reach", info = info)
        expect_true(res$lower_bound, info = info)
        expect_length(st$sizes, 0L)
        next
      }
      expect_identical(st$sizes[1], start, info = info)
      expect_identical(res$value[res$term == "necessary_cluster_size"], m_star, info = info)
      expect_lte(attr(res, "simulation")$search$n_designs, 16L)
    }
  }
  # With 10 clusters and an intraclass correlation of 0.2, the target 1.46
  # starts at 8 (its unstandardized start is 1,875), and 1.465 at 8 (444,
  # from which the search stepped down one size at a time).
  for (target in c(1.46, 1.465)) {
    st <- new.env(parent = emptyenv())
    st$sizes <- numeric(0)
    local_mocked_bindings(.runrep_width_es_crd = crd_mock_width(
      getFromNamespace(".crd_es_width_limit", "DMAR")(10, 5, 0.2, 0.95), 4, 10, 2.5, st))
    res <- ss_aipe_crd_es_n_individuals_fixed_width(
      width = target, n_clusters = 10, true_es = 0.5, es_type = 1, true_icc = 0.2,
      pr_treat = 0.5, G = 2, seed = 1)
    expect_identical(st$sizes[1], 8)
    expect_lte(attr(res, "simulation")$search$n_designs, 4L)
  }
})

test_that("the cluster size search reaches 100 above the unstandardized start", {
  # A large effect: the mocked width is widened by 1 + 40 / m (c = 400 with
  # 10 clusters), and the size the target needs, 192 by enumeration, lies
  # more than 100 above the closed-form start (32) but below the
  # unstandardized start (300), from which the search reached it by
  # stepping down. The search from the closed-form start reaches 400.
  K <- 10; icc <- 0.05; r <- 19; q <- stats::qchisq(0.95, 1)
  e_root <- sqrt(2) * exp(lgamma(4.5) - lgamma(4))
  L <- 2 * sqrt((exp(q / K) - 1) / 2.5) * sqrt(icc / (1 - icc)) * e_root
  ratio <- stats::qt(0.975, 8) / (sqrt(exp(q / K) - 1) * e_root)
  target <- ratio * L * sqrt(1 + r / 300) * 1.000001
  ns <- asNamespace("DMAR")
  expect_identical(ns$.find_n_indiv_crd_diff(width = target, n_clusters = K, pr_treat = 0.5,
                                             true_total_var = 1 / (1 - icc), true_icc = icc), 300)
  m <- 2:5000
  m_star <- m[which(L * sqrt(1 + r / m) * (1 + 400 / (K * m)) < target)[1]]
  expect_identical(m_star, 192L)
  st <- new.env(parent = emptyenv())
  st$sizes <- numeric(0)
  local_mocked_bindings(.runrep_width_es_crd = crd_mock_width(L, r, K, 400, st))
  res <- ss_aipe_crd_es_n_individuals_fixed_width(
    width = target, n_clusters = K, true_es = 3, es_type = 1, true_icc = icc,
    pr_treat = 0.5, G = 2, seed = 1)
  expect_identical(st$sizes[1], floor(r / ((target / L)^2 - 1)) + 1)
  expect_identical(res$value[res$term == "necessary_cluster_size"], 192)
  expect_gt(max(st$sizes), st$sizes[1] + 100)
  # A width that never meets the target is refused at the reach, 500 here
  # (100 above 444, the unstandardized start of 1.465 with 10 clusters and
  # an intraclass correlation of 0.2), and the message says from where.
  local_mocked_bindings(.runrep_width_es_crd = function(...)
    list(width = 2, reason = NA_character_, warnings = character(0)))
  err <- tryCatch(ss_aipe_crd_es_n_individuals_fixed_width(
    width = 1.465, n_clusters = 10, true_es = 0.5, es_type = 1, true_icc = 0.2,
    pr_treat = 0.5, G = 2, seed = 1), error = function(e) e)
  expect_s3_class(err, "dmar_computation_failed")
  expect_identical(err$reason, "width_not_reached")
  expect_identical(c(err$from, err$to), c(8, 500))
  expect_match(conditionMessage(err), paste0("(100 above 444, the cluster size the ",
               "unstandardized difference between the condition means needs, and at most 500"),
               fixed = TRUE)
})

test_that("without a covariate, true_R2_within does not move the closed-form start", {
  # The simulated data use true_R2_within only with a covariate, so the
  # closed-form start leaves it out without one, and the start stays a
  # lower bound. Through the first build of this repair it entered the
  # start: with 10 clusters, an intraclass correlation of 0.2, a target of
  # 1.46, and true_R2_within = -0.5, the search started at 12, above the
  # size the target needs (about 10), and a start beyond 500 was not called
  # a lower bound. A proportion of variance explained below 0 is now
  # refused on entry, as in the two budget planners, so the check is made
  # at values in [0, 1): had 0.9 entered the start, the start would be
  # floor(0.1 * 4 / ...) + 1 rather than 8.
  ns <- asNamespace("DMAR")
  lim <- ns$.crd_es_width_limit(10, 5, 0.2, 0.95)
  call <- function(width, r2w) ss_aipe_crd_es_n_individuals_fixed_width(
    width = width, n_clusters = 10, true_es = 0.5, es_type = 1, true_icc = 0.2,
    pr_treat = 0.5, true_R2_within = r2w, G = 2, seed = 1)
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  err <- tryCatch(call(1.46, -0.5), error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "true_R2_within")
  expect_match(conditionMessage(err), "must be a single number in [0, 1)", fixed = TRUE)
  for (r2w in c(0.3, 0.9)) {
    st <- new.env(parent = emptyenv())
    st$sizes <- numeric(0)
    local_mocked_bindings(.runrep_width_es_crd = crd_mock_width(lim, 4, 10, 2.5, st))
    call(1.46, r2w)
    expect_identical(st$sizes[1], 8, info = paste("true_R2_within", r2w))
    local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
    target <- lim * 1.0005
    err <- tryCatch(call(target, r2w), error = function(e) e)
    expect_identical(err$reason, "cluster_size_beyond_reach", info = paste("true_R2_within", r2w))
    expect_identical(err$start, floor(4 / ((target / lim)^2 - 1)) + 1, info = paste("true_R2_within", r2w))
    expect_true(err$lower_bound, info = paste("true_R2_within", r2w))
    expect_match(conditionMessage(err), "a lower bound on the size the target needs", fixed = TRUE)
  }
})

test_that("a target the t-based start put out of reach is planned with a small cluster size", {
  skip_on_cran()
  skip_if_not_installed("OpenMx")
  # 10 clusters, an intraclass correlation of 0.2, and a target of 1.46:
  # the unstandardized start is 1,875, which was refused as beyond reach
  # and called a lower bound. Nothing is mocked. The closed-form start is
  # 8 by algebra; with seed 11 its simulated width is above the target and
  # that of 9 below it, so the plan is 9 after two designs.
  q <- stats::qchisq(0.95, 1)
  L <- 2 * sqrt((exp(q / 10) - 1) / 2.5) * 0.5 * sqrt(2) * exp(lgamma(4.5) - lgamma(4))
  expect_identical(floor(4 / ((1.46 / L)^2 - 1)) + 1, 8)
  res <- suppressWarnings(ss_aipe_crd_es_n_individuals_fixed_width(
    width = 1.46, n_clusters = 10, true_es = 0.5, es_type = 1, true_icc = 0.2,
    pr_treat = 0.5, G = 20, seed = 11))
  expect_identical(res$value[res$term == "necessary_cluster_size"], 9)
  sim <- attr(res, "simulation")
  expect_identical(sim$search$n_designs, 2L)
  expect_identical(c(sim$attempted, sim$returned, sim$failed), c(20L, 20L, 0L))
  expect_lt(res$value[res$term == "expected_width"], 1.46)
  # A sanity check only, since 20 data sets give a Monte Carlo standard
  # error near 0.08: the width reported at 9 lies within four standard
  # errors of the oracle's mean width there (the next test pins the width
  # at 9 to the oracle with 200 data sets).
  se <- diff(sim$width_mc_ci) / (2 * stats::qnorm(0.975))
  expect_lt(abs(res$value[res$term == "expected_width"] - 1.4688), 4 * sqrt(se^2 + 0.0038^2))
})

test_that("the simulated width at a design of the band agrees with an independent oracle", {
  skip_on_cran()
  skip_if_not_installed("OpenMx")
  # 10 clusters (5 treated) of 9, an intraclass correlation of 0.2, and a
  # standardized effect of 0.5. The oracle computes the profile likelihood
  # interval for the effect from raw individual-level data with no DMAR
  # code (the intercept and the variance of the cluster means profiled in
  # closed form, the within-cluster standard deviation by a one-dimensional
  # optimization, the limits by root finding): mean width 1.4688 with Monte
  # Carlo standard error 0.0038 over 10,000 data sets, all usable (2,000 at
  # 1.4699 and 8,000 at 1.4685; evidence/j4_crd_start/02_oracle_band.log and
  # round2/r2_02_oracle.log of the November 2026 review). 200 data sets give a
  # standard error near 0.027, so the check is within three standard errors
  # of the difference, about 0.08.
  w <- suppressWarnings(getFromNamespace(".find_width_crd_es", "DMAR")(
    G = 200, n_clusters = 10, n_treat_clus = 5, cluster_size = 9, true_icc = 0.2, es = 0.5,
    es_type = 1, true_total_var = 1.25, seed = 11))
  sim <- attr(w, "simulation")
  expect_identical(c(sim$attempted, sim$returned, sim$failed), c(200L, 200L, 0L))
  se <- diff(sim$width_mc_ci) / (2 * stats::qnorm(0.975))
  expect_lt(se, 0.035)
  expect_lt(abs(as.numeric(w) - 1.4688), 3 * sqrt(se^2 + 0.0038^2))
})

## Entry checks, budgets, and the start of the searches (the round 6
## review of November 2026: SQ02 to SQ06 and the problems its evaluation
## found). Nothing here is simulated with OpenMx: the refusals come before
## any simulation, and the searches run on mocks.

crd_es_err <- function(expr) tryCatch(expr, error = function(e) e)

test_that("G above the pool of seeds is refused on entry, and the seeds drawn are unchanged", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  calls <- list(
    function(G) ss_aipe_crd_es_n_clusters_fixed_width(width = 0.4, cluster_size = 20,
      true_es = 0.5, true_icc = 0.1, pr_treat = 0.5, G = G),
    function(G) ss_aipe_crd_es_n_individuals_fixed_width(width = 0.4, n_clusters = 40,
      true_es = 0.5, true_icc = 0.1, pr_treat = 0.5, G = G),
    function(G) ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 10,
      cluster_cost = 0, individual_cost = 1, G = G, pr_treat = 0.5, true_icc = 0.2, true_es = 0.5),
    function(G) ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 30,
      cluster_cost = 0, individual_cost = 1, G = G, pr_treat = 0.5, true_icc = 0.2, true_es = 0.5),
    function(G) ss_aipe_crd_es_both_fixed_budget(budget = 1000, cluster_cost = 5,
      individual_cost = 1, true_es = 0.5, true_icc = 0.2, pr_treat = 0.5, G = G),
    function(G) ss_aipe_crd_es_both_fixed_width(width = 0.5, cluster_cost = 5,
      individual_cost = 1, true_es = 0.5, true_icc = 0.25, pr_treat = 0.5, G = G))
  for (i in seq_along(calls)) {
    # A million replications stopped with "cannot take a sample larger than
    # the population when 'replace = FALSE'".
    err <- crd_es_err(calls[[i]](1000000))
    expect_s3_class(err, "dmar_count_invalid")
    expect_identical(err$arg, "G")
    expect_identical(err$max, 999999L)
    expect_match(conditionMessage(err), "at most 999999", fixed = TRUE)
  }
  expect_identical(DMAR:::.dmar_check_count(999999, "G", max = DMAR:::.CRD_ES_MAX_G), 999999L)
  # Each simulated data set takes the seed sample(1:999999, G) gives.
  st <- new.env(parent = emptyenv())
  st$seeds <- integer(0)
  local_mocked_bindings(.runrep_width_es_crd = function(seed, ...) {
    st$seeds <- c(st$seeds, seed)
    list(width = 0.5, reason = NA_character_, warnings = character(0))
  })
  set.seed(113)
  expected <- sample(1:999999, 5)
  getFromNamespace(".find_width_crd_es", "DMAR")(G = 5, n_clusters = 20, n_treat_clus = 10,
                                                cluster_size = 10, true_icc = 0.2, es = 0.5, seed = 113)
  expect_identical(st$seeds, expected)
})

test_that("costs of 0 that leave the design open are refused before anything is simulated", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  # These searched without end (the budget planners) or stopped with
  # "missing value where TRUE/FALSE needed" (the width planner).
  calls <- list(
    n_clusters_budget = function(cc, ic) ss_aipe_crd_es_n_clusters_fixed_budget(
      budget = 1000, cluster_size = 20, cluster_cost = cc, individual_cost = ic),
    n_individuals_budget = function(cc, ic) ss_aipe_crd_es_n_individuals_fixed_budget(
      budget = 1000, n_clusters = 20, cluster_cost = cc, individual_cost = ic),
    both_budget = function(cc, ic) ss_aipe_crd_es_both_fixed_budget(
      budget = 10000, cluster_cost = cc, individual_cost = ic, true_es = 0.5, true_icc = 0.2,
      pr_treat = 0.5, G = 2, seed = 113),
    both_width = function(cc, ic) ss_aipe_crd_es_both_fixed_width(
      width = 0.5, cluster_cost = cc, individual_cost = ic, true_es = 0.5, true_icc = 0.2,
      pr_treat = 0.5, G = 2, seed = 113))
  for (nm in names(calls)) {
    expect_identical(crd_es_err(calls[[nm]](0, 0))$reason, "costs_all_zero", info = nm)
    if (nm != "n_clusters_budget") {
      expect_identical(crd_es_err(calls[[nm]](10, 0))$reason, "individual_cost_zero", info = nm)
    }
  }
  expect_identical(calls$n_clusters_budget(10, 0)$value[1], 100)
  # The planners that take a width and a fixed size only price the plan.
  expect_error(ss_aipe_crd_es_n_individuals_fixed_width(width = 0.4, n_clusters = 40,
    true_es = 0.5, true_icc = 0.1, pr_treat = 0.5, G = 2, cluster_cost = -1, individual_cost = 1),
    class = "dmar_argument_invalid")
})

test_that("a budget that buys no design the planner considers is refused with the cost of the smallest", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  err <- crd_es_err(ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1, cluster_size = 10,
                                                           cluster_cost = 10, individual_cost = 1))
  expect_identical(err$reason, "budget_infeasible")
  expect_identical(err$min_cost, 80)
  err <- crd_es_err(ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1, n_clusters = 10,
                                                              cluster_cost = 10, individual_cost = 1))
  expect_identical(err$reason, "budget_infeasible")
  expect_identical(err$min_cost, 120)
  # It stopped with "invalid 'type' (list) of argument".
  err <- crd_es_err(ss_aipe_crd_es_both_fixed_budget(budget = 1, cluster_cost = 10,
    individual_cost = 1, true_es = 0.5, true_icc = 0.2, pr_treat = 0.5, G = 2, seed = 113))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$reason, "budget_infeasible")
  expect_identical(c(err$n_clusters, err$min_cost), c(4, 48))
  # A budget that buys a design with no treated cluster is refused before
  # the width is simulated.
  err <- crd_es_err(ss_aipe_crd_es_n_clusters_fixed_budget(budget = 100, cluster_size = 20,
    cluster_cost = 0, individual_cost = 1, G = 2, pr_treat = 0.1, true_icc = 0.2, true_es = 0.5))
  expect_identical(err$reason, "arm_without_clusters")
})

test_that("the standardized planners check their design arguments and width on entry", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  w <- function(...) crd_es_err(do.call(ss_aipe_crd_es_n_clusters_fixed_width, utils::modifyList(
    list(width = 0.4, cluster_size = 20, true_es = 0.5, true_icc = 0.1, pr_treat = 0.5, G = 2),
    list(...))))
  expect_identical(w(cluster_size = 20.5)$arg, "cluster_size")
  expect_identical(w(pr_treat = 1)$arg, "pr_treat")
  expect_identical(w(true_R2_between = 1.2)$arg, "true_R2_between")
  expect_identical(w(n_covariates = 2)$arg, "n_covariates")
  # A width of 0 or below searched without end.
  err <- crd_es_err(ss_aipe_crd_es_n_clusters_fixed_width(width = -1, cluster_size = 20,
    true_es = 0.5, true_icc = 0.1, pr_treat = 0.5, G = 2))
  expect_identical(err$arg, "width")
  expect_identical(crd_es_err(ss_aipe_crd_es_both_fixed_width(width = 0, cluster_cost = 5,
    individual_cost = 1, true_es = 0.5, true_icc = 0.25, pr_treat = 0.5, G = 2))$arg, "width")
  expect_identical(crd_es_err(ss_aipe_crd_es_n_individuals_fixed_width(width = 0.4,
    n_clusters = 2, true_es = 0.5, true_icc = 0.1, pr_treat = 0.5, G = 2))$arg, "n_clusters")
  expect_identical(crd_es_err(ss_aipe_crd_es_n_individuals_fixed_width(width = 0.4,
    n_clusters = 5, true_es = 0.5, true_icc = 0.1, pr_treat = 0.1, G = 2))$reason,
    "arm_without_clusters")
  expect_identical(crd_es_err(ss_aipe_crd_es_n_individuals_fixed_budget(
    budget = 1000, n_clusters = 30.5, cluster_cost = 0, individual_cost = 1))$arg, "n_clusters")
  # With es_type = 0 and 10 clusters no cluster size up to 100,000 brings
  # the unstandardized width to 1 (it falls toward 1.30), so the search has
  # no start.
  err <- crd_es_err(ss_aipe_crd_es_n_individuals_fixed_width(width = 1, n_clusters = 10,
    true_es = 0.5, es_type = 0, true_icc = 0.2, pr_treat = 0.5, G = 2))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$reason, "cluster_size_beyond_reach")
  expect_identical(err$start, Inf)
})

test_that("the search for the number of clusters starts at the floor of two clusters in each condition", {
  # Every number of clusters meets a wide target, and the unstandardized
  # start was 3 (at pr_treat = .1, no treated cluster); it is the floor,
  # 20, and the downward search stops above one treated cluster.
  st <- new.env(parent = emptyenv())
  st$n <- numeric(0)
  local_mocked_bindings(.find_width_crd_es = function(..., n_clusters) {
    st$n <- c(st$n, n_clusters)
    0.1
  })
  res <- asNamespace("DMAR")$.find_n_clus_crd_es(width = 50, cluster_size = 20, es = 0.5,
                                                 true_icc = 0.2, pr_treat = 0.1, G = 2)
  expect_identical(st$n[1], 20)
  expect_identical(res[1], 15)
  expect_identical(round(res[1] * 0.1), 2)
})

test_that("the both-size budget search does not simulate a number of clusters the budget does not buy", {
  # At an intraclass correlation of .9 the narrowest unstandardized design
  # is the most clusters of 2 the budget buys, 142, and one more cluster is
  # beyond it. That neighbor was simulated with a cluster size of NA,
  # which stopped with "invalid 'times' argument".
  local_mocked_bindings(.find_width_crd_es = function(..., n_clusters, cluster_size) {
    if (is.na(cluster_size)) stop("a design the budget does not buy was simulated")
    1 / n_clusters
  })
  res <- asNamespace("DMAR")$.find_min_width_crd_es(budget = 1000, cluster_cost = 5,
    individual_cost = 1, es = 0.5, true_icc = 0.9, pr_treat = 0.5, G = 2)
  expect_identical(res, c(142, 2, 1 / 142))
})

test_that("es_type and true_icc are checked wherever a width is simulated", {
  skip_on_cran()
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  calls <- list(
    n_clusters_width = function(icc, type) ss_aipe_crd_es_n_clusters_fixed_width(width = 1,
      cluster_size = 20, true_es = 0.5, true_icc = icc, es_type = type, pr_treat = 0.5, G = 2, seed = 113),
    n_individuals_width = function(icc, type) ss_aipe_crd_es_n_individuals_fixed_width(width = 1,
      n_clusters = 30, true_es = 0.5, true_icc = icc, es_type = type, pr_treat = 0.5, G = 2, seed = 113),
    n_clusters_budget = function(icc, type) ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000,
      cluster_size = 10, cluster_cost = 0, individual_cost = 1, G = 2, pr_treat = 0.5, true_icc = icc,
      true_es = 0.5, es_type = type, seed = 113),
    n_individuals_budget = function(icc, type) ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000,
      n_clusters = 20, cluster_cost = 0, individual_cost = 1, G = 2, pr_treat = 0.5, true_icc = icc,
      true_es = 0.5, es_type = type, seed = 113),
    both_budget = function(icc, type) ss_aipe_crd_es_both_fixed_budget(budget = 3000, cluster_cost = 10,
      individual_cost = 1, true_es = 0.5, true_icc = icc, es_type = type, pr_treat = 0.5, G = 2, seed = 113),
    both_width = function(icc, type) ss_aipe_crd_es_both_fixed_width(width = 1, cluster_cost = 10,
      individual_cost = 1, true_es = 0.5, true_icc = icc, es_type = type, pr_treat = 0.5, G = 2, seed = 113))
  # true_icc = 1 or 1.5, and 0 with es_type = 2, stopped with "missing
  # value where TRUE/FALSE needed" (the first build of this repair searched
  # without end); true_icc = 1 with es_type = 0 simulated three designs at
  # which every fit failed.
  cases <- list(list(1, 1, "true_icc"), list(1.5, 1, "true_icc"), list(-0.2, 1, "true_icc"),
                list(1, 0, "true_icc"), list(NA_real_, 1, "true_icc"), list(0, 2, c("true_icc", "es_type")),
                list(0.2, 7, "es_type"))
  for (nm in names(calls)) {
    for (cs in cases) {
      err <- crd_es_err(calls[[nm]](cs[[1]], cs[[2]]))
      expect_s3_class(err, "dmar_argument_invalid")
      expect_identical(err$arg, cs[[3]], info = paste(nm, cs[[1]], cs[[2]]))
    }
  }
  # Without G nothing is simulated and true_icc is not used.
  expect_identical(ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 10,
    cluster_cost = 0, individual_cost = 1, true_icc = 1.5)$value[1], 100)
})

test_that("diff_size and the covariate are checked on entry by the standardized planners", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  w <- function(...) crd_es_err(do.call(ss_aipe_crd_es_n_clusters_fixed_width, utils::modifyList(
    list(width = 1, cluster_size = 20, true_es = 0.5, true_icc = 0.2, pr_treat = 0.5, G = 2), list(...))))
  for (bad in list(numeric(0), c(Inf, 0), NA)) expect_identical(w(diff_size = bad)$arg, "diff_size")
  expect_identical(crd_es_err(ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 10,
    cluster_cost = 0, individual_cost = 1, diff_size = numeric(0)))$arg, "diff_size")
  # true_icc_covariate = 1.5 simulated and stopped with no_finite_width.
  err <- w(n_covariates = 1, true_icc_covariate = 1.5, true_R2_within = 0.1, true_R2_between = 0.1)
  expect_identical(err$arg, "true_icc_covariate")
})

test_that("a between-level R2 without a covariate is refused by the six standardized planners", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  # With n_covariates = 0, true_R2_between = .9 moved the start of the
  # search for a width of 1 at clusters of 20 from 22 clusters to 8, while
  # the simulated data left it out: both values planned 18 clusters.
  pop <- list(true_es = 0.5, true_icc = 0.2, pr_treat = 0.5, G = 2, true_R2_between = 0.9)
  calls <- list(
    ss_aipe_crd_es_n_clusters_fixed_width = list(width = 1, cluster_size = 20),
    ss_aipe_crd_es_n_individuals_fixed_width = list(width = 1, n_clusters = 20),
    ss_aipe_crd_es_n_clusters_fixed_budget = list(budget = 1000, cluster_size = 20,
                                                  cluster_cost = 0, individual_cost = 1),
    ss_aipe_crd_es_n_individuals_fixed_budget = list(budget = 1000, n_clusters = 20,
                                                     cluster_cost = 0, individual_cost = 1),
    ss_aipe_crd_es_both_fixed_budget = list(budget = 1000, cluster_cost = 10, individual_cost = 1),
    ss_aipe_crd_es_both_fixed_width = list(width = 1, cluster_cost = 10, individual_cost = 1))
  for (fn in names(calls)) {
    err <- crd_es_err(do.call(fn, c(calls[[fn]], pop)))
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "true_R2_between", info = fn)
    expect_identical(err$reason, "R2_between_without_covariates", info = fn)
  }
  # Without G the budget planners simulate nothing and still refuse it.
  err <- crd_es_err(ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 20,
    cluster_cost = 0, individual_cost = 1, true_R2_between = 0.9))
  expect_identical(err$reason, "R2_between_without_covariates")
  # One covariate passes the check; so does true_R2_within without one.
  ns <- asNamespace("DMAR")
  expect_null(ns$.crd_check_design(1, 0.9, 0, pr_treat = 0.5))
  expect_null(ns$.crd_check_design(0, 0, 0.5, pr_treat = 0.5))
})

test_that("the standardized budget planners with a fixed size check assurance on entry", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  # An assurance of 2 stopped with "'probs' outside [0,1]" in the budget
  # planners with G, and was not checked without it.
  budget_calls <- list(
    function(...) ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, cluster_size = 10,
      cluster_cost = 10, individual_cost = 1, ...),
    function(...) ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 20,
      cluster_cost = 10, individual_cost = 1, ...))
  for (f in budget_calls) {
    expect_s3_class(crd_es_err(f(assurance = 2)), "dmar_assurance_invalid")
    expect_s3_class(crd_es_err(f(assurance = 2, G = 2, pr_treat = 0.5, true_icc = 0.2, true_es = 0.5,
                                 seed = 113)), "dmar_assurance_invalid")
  }
})

test_that("a fixed number of clusters with an empty condition is refused only where a width is simulated", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  # Without G pr_treat is not used, and the size the budget buys is
  # returned, as by ss_aipe_crd_es_n_clusters_fixed_budget().
  res <- ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 5, cluster_cost = 0,
                                                   individual_cost = 1, pr_treat = 0.1)
  expect_identical(res$value[1], 200)
  err <- crd_es_err(ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 5,
    cluster_cost = 0, individual_cost = 1, pr_treat = 0.1, G = 2, true_icc = 0.2, true_es = 0.5))
  expect_identical(err$reason, "arm_without_clusters")
  expect_identical(err$fewest, 6)
})

test_that("a budget below the smallest design with a width names the budget that buys one", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  err <- crd_es_err(ss_aipe_crd_es_n_clusters_fixed_budget(budget = 60, cluster_size = 20,
    cluster_cost = 0, individual_cost = 1, G = 2, pr_treat = 0.1, true_icc = 0.2, true_es = 0.5))
  expect_identical(err$reason, "budget_infeasible")
  expect_identical(c(err$n_clusters, err$min_cost), c(6, 120))
  err <- crd_es_err(ss_aipe_crd_es_n_clusters_fixed_budget(budget = 60, cluster_size = 20,
    cluster_cost = 0, individual_cost = 1))
  expect_identical(c(err$n_clusters, err$min_cost), c(4, 80))
})

test_that("a start beyond a million clusters is refused before anything is simulated", {
  local_mocked_bindings(.runrep_width_es_crd = function(...) stop("nothing should be simulated"))
  err <- crd_es_err(ss_aipe_crd_es_n_clusters_fixed_width(width = 1e-300, cluster_size = 20,
    true_es = 0.5, true_icc = 0.2, pr_treat = 0.5, G = 2))
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$reason, "n_clusters_beyond_reach")
  expect_match(conditionMessage(err), "^ss_aipe_crd_es_n_clusters_fixed_width\\(\\): no design")
  err <- crd_es_err(ss_aipe_crd_es_both_fixed_width(width = 1e-300, cluster_cost = 10,
    individual_cost = 1, true_es = 0.5, true_icc = 0.2, pr_treat = 0.5, G = 2))
  expect_identical(err$reason, "n_clusters_beyond_reach")
})
