# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# Tests for ci_bifactor_indices(): bootstrap intervals for the indices
# bifactor_indices() reports. The fits are of simulated nine-item
# bifactor data (one general and three group factors), chosen so that a
# replication takes a few hundredths of a second: a clean model whose
# resamples are proper (general loadings .6, group loadings .55; with
# .5 and .5 at n = 300 a resample now and then gives a solution that
# runs away, which takes seconds to fit), a model with one item near the
# boundary whose resamples are often improper without running away
# (n = 300), and the same model at n = 150, whose point estimate is
# itself improper. Each replication is checked against a hand refit of
# the same resample, and each interval against the engine run on that
# index alone.

skip_if_not_installed("lavaan")
skip_if_not_installed("boot")

bf_model <- "g  =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
             f1 =~ x1 + x2 + x3
             f2 =~ x4 + x5 + x6
             f3 =~ x7 + x8 + x9"

bf_data <- function(n, seed, lg = rep(0.6, 9), ls = rep(0.55, 9)) {
  set.seed(seed)
  g <- rnorm(n)
  grp <- list(rnorm(n), rnorm(n), rnorm(n))
  X <- vapply(1:9, function(i) {
    lg[i] * g + ls[i] * grp[[ceiling(i / 3)]] +
      sqrt(1 - lg[i]^2 - ls[i]^2) * rnorm(n)
  }, numeric(n))
  colnames(X) <- paste0("x", 1:9)
  as.data.frame(X)
}

bf_fit <- function(data, ...) {
  lavaan::cfa(bf_model, data = data, orthogonal = TRUE, std.lv = TRUE, ...)
}

# x1 near the boundary: a general loading of .9 and a group loading of
# .42 leave it a residual variance of about .01.
near_lg <- c(0.9, rep(0.6, 8))
near_ls <- c(0.42, rep(0.5, 8))

# The resamples .bootstrap_ci() draws for a seed: boot draws its whole
# index array before any replication, whatever the statistic is.
bf_indices <- function(data, B, seed) {
  bo <- withr::with_seed(seed, boot::boot(data, function(d, i) 0, R = B,
                                          stype = "i"))
  boot::boot.array(bo, indices = TRUE)
}

bf_pool <- function(fit) {
  X <- as.data.frame(lavaan::lavInspect(fit, "data"))
  X[rowSums(!is.na(X)) > 0L, , drop = FALSE]
}

bf_quietly <- function(expr) suppressWarnings(expr)


test_that("one row per index, with the point estimates of bifactor_indices()", {
  fit <- bf_fit(bf_data(300, 113))
  ci <- bf_quietly(ci_bifactor_indices(fit, B = 20, seed = 113))
  bi <- bifactor_indices(fit)
  expect_s3_class(ci, c("dmar_ci_bifactor", "dmar_tbl", "data.frame"))
  expect_named(ci, c("factor", "index", "estimate", "se", "ci_lower",
                     "ci_upper"))
  expect_identical(ci$factor, c(rep("g", 6), rep(c("f1", "f2", "f3"),
                                                 each = 5)))
  expect_identical(ci$index, c("ECV", "omega", "omega_H", "PUC", "H",
                               "maximal_reliability",
                               rep(c("ECV", "omega", "omega_HS", "H",
                                     "maximal_reliability"), 3)))
  maximal <- ci$index == "maximal_reliability"
  for (k in which(!maximal)) {
    expect_identical(ci$estimate[k],
                     bi[[ci$index[k]]][bi$factor == ci$factor[k]])
  }
  expect_identical(ci$estimate[maximal],
                   unname(attr(bi, "maximal_reliability")[ci$factor[maximal]]))
  puc <- ci$index == "PUC"
  expect_true(all(is.na(c(ci$se[puc], ci$ci_lower[puc],
                          ci$ci_upper[puc]))))
  expect_true(all(is.finite(ci$se[!puc])))
  expect_true(all(ci$ci_lower[!puc] < ci$ci_upper[!puc]))
  expect_identical(attr(ci, "B"), 20L)
  expect_identical(attr(ci, "B_used"), attr(ci, "bootstrap")$returned)
  expect_identical(attr(ci, "conf_level"), 0.95)
  expect_identical(attr(ci, "ci_method_requested"), "percentile")
  expect_identical(unname(attr(ci, "ci_method")), rep("percentile", 20L))
  expect_identical(names(attr(ci, "ci_method")),
                   paste(ci$factor, ci$index, sep = ":")[!puc])
  expect_identical(attr(ci, "improper_rule"), "drop")
  expect_identical(attr(ci, "denominator"), "model_implied")
  expect_identical(attr(ci, "metric"), "observed")
  expect_identical(attr(ci, "general"), "g")
  expect_false(attr(ci, "improper"))
  expect_identical(attr(ci, "N"), 300L)
  expect_identical(attr(ci, "N_incomplete"), 0L)
  # bifactor_indices() itself is unchanged by any of this.
  expect_identical(names(bi), c("factor", "ECV", "omega", "omega_H",
                                "omega_HS", "PUC", "H"))
})

test_that("every replication is a refit of the model to its resample", {
  skip_on_cran()
  data <- bf_data(300, 113)
  fit <- bf_fit(data)
  ci <- bf_quietly(ci_bifactor_indices(fit, B = 20, seed = 113))
  bs <- attr(ci, "bootstrap")
  expect_identical(dim(bs$t), c(20L, 20L))
  expect_identical(colnames(bs$t), names(attr(ci, "ci_method")))
  expect_identical(levels(bs$status),
                   c("proper", "improper", "error", "non_finite",
                     "undetermined", "not_converged"))
  expect_identical(bs$attempted, 20L)
  expect_identical(bs$returned + bs$failed, 20L)
  expect_identical(sum(bs$status == "proper"), bs$returned)
  expect_identical(bs$returned, 20L)
  idx <- bf_indices(bf_pool(fit), B = 20, seed = 113)
  for (r in c(1L, 7L, 20L)) {
    refit <- bf_fit(bf_pool(fit)[idx[r, ], ])
    hand <- bf_quietly(bifactor_indices(refit, general = "g"))
    cells <- DMAR:::.bifactor_cells(hand)
    cells <- cells[cells$index != "PUC", ]
    expect_equal(unname(bs$t[r, ]), DMAR:::.bifactor_cell_values(hand, cells),
                 tolerance = 1e-10)
  }
})

test_that("each interval is the engine's interval for that index alone, on the same replications", {
  skip_on_cran()
  fit <- bf_fit(bf_data(300, 113, near_lg, near_ls))
  ci <- bf_quietly(ci_bifactor_indices(fit, B = 30, seed = 113))
  pool <- bf_pool(fit)
  pt <- lavaan::parTable(fit)
  fn <- DMAR:::.bifactor_point_fn(
    model = pt[setdiff(names(pt), c("est", "se", "start"))],
    options = DMAR:::.bifactor_refit_options(fit, "model_implied"),
    ordered = character(0), categories = list(), general = "g",
    denominator = "model_implied",
    cells = DMAR:::.bifactor_cells(bifactor_indices(fit))[-4, ])
  stat_names <- names(attr(ci, "ci_method"))
  for (j in match(c("g:omega_H", "f2:omega_HS"), stat_names)) {
    one <- function(d) {
      v <- fn(d)
      out <- v[j]
      attr(out, "dmar_reason") <- attr(v, "dmar_reason")
      out
    }
    r <- bf_quietly(DMAR:::.bootstrap_ci(pool, one, B = 30, conf_level = 0.95,
                                         kind = "percentile", seed = 113,
                                         improper = "drop"))
    row <- which(paste(ci$factor, ci$index, sep = ":") == stat_names[j])
    expect_identical(c(ci$ci_lower[row], ci$ci_upper[row], ci$se[row]),
                     c(r$lower, r$upper, r$se))
    expect_identical(attr(ci, "bootstrap")$returned, r$bootstrap$returned)
    expect_identical(attr(ci, "bootstrap")$status, r$bootstrap$status)
  }
})

test_that("the maximal reliability of each factor has an interval from the same replications", {
  skip_on_cran()
  # bifactor_indices() gives the maximal reliability of each factor in an
  # attribute, and ci_bifactor_indices() computed it on every replication
  # and dropped it: the table had no row for it and the footer pointed to
  # the bare attribute.
  fit <- bf_fit(bf_data(300, 113))
  ci <- bf_quietly(ci_bifactor_indices(fit, B = 20, seed = 113))
  maximal <- ci$index == "maximal_reliability"
  expect_identical(ci$factor[maximal], c("g", "f1", "f2", "f3"))
  expect_identical(ci$estimate[maximal],
                   unname(attr(bifactor_indices(fit), "maximal_reliability")))
  expect_true(all(is.finite(ci$se[maximal])))
  expect_true(all(ci$ci_lower[maximal] < ci$ci_upper[maximal]))
  expect_true(all(ci$ci_lower[maximal] > 0 & ci$ci_upper[maximal] < 1))
  # Each replication's value is that of a hand refit of its resample.
  bs <- attr(ci, "bootstrap")
  cols <- paste0(c("g", "f1", "f2", "f3"), ":maximal_reliability")
  expect_true(all(cols %in% colnames(bs$t)))
  idx <- bf_indices(bf_pool(fit), B = 20, seed = 113)
  for (r in c(1L, 20L)) {
    refit <- bf_fit(bf_pool(fit)[idx[r, ], ])
    hand <- attr(bf_quietly(bifactor_indices(refit, general = "g")),
                 "maximal_reliability")
    expect_equal(unname(bs$t[r, cols]), unname(hand), tolerance = 1e-10)
  }
  # Each interval is the engine's interval for that value alone.
  pool <- bf_pool(fit)
  pt <- lavaan::parTable(fit)
  cells <- DMAR:::.bifactor_cells(bifactor_indices(fit))
  fn <- DMAR:::.bifactor_point_fn(
    model = pt[setdiff(names(pt), c("est", "se", "start"))],
    options = DMAR:::.bifactor_refit_options(fit, "model_implied"),
    ordered = character(0), categories = list(), general = "g",
    denominator = "model_implied", cells = cells[cells$index != "PUC", ])
  stat_names <- names(attr(ci, "ci_method"))
  for (j in match(c("g:maximal_reliability", "f2:maximal_reliability"),
                  stat_names)) {
    one <- function(d) {
      v <- fn(d)
      out <- v[j]
      attr(out, "dmar_reason") <- attr(v, "dmar_reason")
      out
    }
    r <- bf_quietly(DMAR:::.bootstrap_ci(pool, one, B = 20, conf_level = 0.95,
                                         kind = "percentile", seed = 113,
                                         improper = "drop"))
    row <- which(paste(ci$factor, ci$index, sep = ":") == stat_names[j])
    expect_identical(c(ci$ci_lower[row], ci$ci_upper[row], ci$se[row]),
                     c(r$lower, r$upper, r$se))
  }
  # The bootstrap estimate of each index's bias is the mean of the
  # replicates its interval uses minus the estimate, and the printed table
  # gives it for H and the maximal reliability, which are biased upward.
  key <- paste(ci$factor, ci$index, sep = ":")
  expect_identical(names(bs$bias), colnames(bs$t))
  expect_equal(unname(bs$bias),
               unname(colMeans(bs$t[bs$status == "proper", , drop = FALSE]) -
                        ci$estimate[match(colnames(bs$t), key)]),
               tolerance = 1e-12)
  foot <- DMAR:::.dmar_tbl_footer(ci)
  expect_match(foot[["bifactor_bias"]], paste0(
    "^Bootstrap estimates of bias \\(mean of the replicates used minus ",
    "the estimate\\): H g [-0-9.e]+, f1 .*; maximal_reliability g ",
    "[-0-9.e]+, f1 [-0-9.e]+, f2 [-0-9.e]+, f3 [-0-9.e]+$"))
  expect_no_match(foot[["bifactor_bias"]], "omega|ECV")
  # The printed table says what the rows are, and the table of
  # bifactor_indices() names this function for their intervals.
  expect_no_match(paste(foot, collapse = "\n"),
                  "attribute of bifactor_indices")
  bi_foot <- DMAR:::.dmar_tbl_footer(bifactor_indices(fit))
  expect_match(bi_foot[names(bi_foot) == "bifactor_definitions"][2],
               "ci_bifactor_indices\\(\\) gives each an interval")
})

test_that("a statistic the point function does not define for a replication leaves the replication to the others", {
  # The engine sets a replication aside whole when any statistic is not
  # finite, unless the point function marks that statistic not defined
  # for the replication (dmar_undefined), as the bifactor point function
  # marks the maximal reliability of an improper solution kept under
  # improper = "keep".
  v <- structure(c(0.8, NA), dmar_reason = "improper",
                 dmar_undefined = c(FALSE, TRUE))
  expect_identical(DMAR:::.boot_classify(v, "keep", 2L), c(0.8, NA, 2, 0.8, NA))
  expect_identical(DMAR:::.boot_classify(v, "drop", 2L),
                   c(NA, NA, 2, 0.8, NA))
  unmarked <- structure(c(0.8, NA), dmar_reason = "improper")
  expect_identical(DMAR:::.boot_classify(unmarked, "keep", 2L),
                   c(NA, NA, 2, NA, NA))
  # The other statistics still decide, and the first is never excused.
  bad <- structure(c(NA, 0.5), dmar_undefined = c(FALSE, TRUE))
  expect_identical(DMAR:::.boot_classify(bad, "keep", 2L),
                   c(NA, NA, 4, NA, NA))
  first <- structure(c(NA, 0.5), dmar_undefined = c(TRUE, FALSE))
  expect_identical(DMAR:::.boot_classify(first, "keep", 2L),
                   c(NA, NA, 4, NA, NA))
  # A marked statistic is NA even if a value came with it.
  marked <- structure(c(0.8, 0.3), dmar_undefined = c(FALSE, TRUE))
  expect_identical(DMAR:::.boot_classify(marked, "keep", 2L),
                   c(0.8, NA, 1, 0.8, NA))
})

test_that("improper replications are set aside by default, kept on request, and counted either way", {
  skip_on_cran()
  fit <- bf_fit(bf_data(300, 113, near_lg, near_ls))
  w_drop <- testthat::capture_warnings(
    d <- ci_bifactor_indices(fit, B = 40, seed = 113))
  w_keep <- testthat::capture_warnings(
    k <- ci_bifactor_indices(fit, B = 40, seed = 113, improper = "keep"))
  bd <- attr(d, "bootstrap")
  bk <- attr(k, "bootstrap")
  n_imp <- sum(bd$status == "improper")
  expect_gt(n_imp, 0L)
  # The same resamples and the same refits under either rule.
  expect_identical(bd$status, bk$status)
  expect_identical(bd$t, bk$t)
  # "drop" sets every improper replication aside.
  expect_identical(bd$improper, n_imp)
  expect_identical(bd$improper_kept, 0L)
  expect_identical(bd$returned, sum(bd$status == "proper"))
  expect_identical(attr(d, "B_used"), bd$returned)
  expect_identical(bd$improper_rule, "drop")
  # "keep" uses those whose indices are proportions, which here is all.
  kept <- bk$status == "improper" & is.finite(bk$t[, 1])
  expect_identical(bk$improper_kept, sum(kept))
  expect_identical(bk$improper, n_imp - sum(kept))
  expect_identical(bk$returned, sum(bk$status == "proper") + sum(kept))
  expect_false(identical(d$ci_lower, k$ci_lower))
  expect_identical(d$estimate, k$estimate)
  # The maximal reliability is not defined for an improper solution, so
  # the improper replications kept have none, and its intervals rest on
  # the proper replications under either rule: the same intervals.
  maximal <- d$index == "maximal_reliability"
  expect_true(all(is.na(bk$t[kept, endsWith(colnames(bk$t),
                                            ":maximal_reliability")])))
  expect_identical(d[maximal, c("se", "ci_lower", "ci_upper")],
                   k[maximal, c("se", "ci_lower", "ci_upper")])
  # So is the bootstrap estimate of its bias, while the other indices'
  # biases are taken over the replications each interval uses.
  mr <- endsWith(names(bd$bias), ":maximal_reliability")
  expect_identical(bd$bias[mr], bk$bias[mr])
  used <- bk$status == "proper" | kept
  expect_equal(unname(bk$bias[!mr]),
               unname(colMeans(bk$t[used, !mr, drop = FALSE]) -
                        k$estimate[match(names(bk$bias)[!mr],
                                         paste(k$factor, k$index,
                                               sep = ":"))]),
               tolerance = 1e-12)
  # The warnings and the printed lines say what was done.
  expect_match(w_drop, sprintf(paste0(
    "%d of 40 bootstrap replications were set aside ",
    "\\(%d gave an improper \\(Heywood\\) solution\\)"), n_imp, n_imp),
    all = FALSE)
  expect_match(w_drop, "conditions the intervals on admissible solutions",
               all = FALSE)
  foot_d <- DMAR:::.dmar_tbl_footer(d)
  foot_k <- DMAR:::.dmar_tbl_footer(k)
  expect_match(foot_d[["bifactor_interval"]], "improper replications set aside")
  expect_match(foot_k[["bifactor_interval"]], paste0(
    "improper replications kept where every index but the maximal ",
    "reliability is a proportion"))
  expect_match(DMAR:::.bifactor_improper_note("keep"), paste0(
    "set aside only when an index other than the maximal reliability is ",
    "not finite"))
  expect_match(foot_k[["bifactor_interval"]], sprintf(paste0(
    "maximal reliability is not defined for an improper solution, so its ",
    "intervals rest on the %d proper replications of the %d used"),
    sum(bk$status == "proper"), bk$returned))
  expect_no_match(foot_d[["bifactor_interval"]], "proper replications of")
  expect_match(foot_k[["bootstrap"]],
               sprintf(paste0("%d of the %d came from an improper ",
                              "\\(Heywood\\) solution and were kept"),
                       sum(kept), bk$returned))
  expect_match(w_keep, "Caution:", all = FALSE)
  expect_match(foot_d[["bootstrap"]], sprintf("%d gave an improper", n_imp))
  # A share of improper replications this large reaches the tail area.
  expect_true(bd$caution)
  expect_match(foot_d[["bootstrap"]], "caution")
  expect_match(w_drop, "Caution:", all = FALSE)
})

test_that("an improper point estimate stops under drop, and under keep is reported with a warning", {
  skip_on_cran()
  fit <- bf_quietly(bf_fit(bf_data(150, 113, near_lg, near_ls)))
  expect_true(attr(bf_quietly(bifactor_indices(fit)), "improper"))
  expect_error(ci_bifactor_indices(fit, B = 10, seed = 113),
               "improper = \"keep\" to bootstrap it as it is")
  w <- testthat::capture_warnings(
    k <- ci_bifactor_indices(fit, B = 30, seed = 113, improper = "keep"))
  expect_match(w, "^Improper \\(Heywood\\) solution", all = FALSE)
  expect_true(attr(k, "improper"))
  bk <- attr(k, "bootstrap")
  # A replication whose indices are not all proportions yields no
  # coefficient and is set aside even under "keep"; its row of t is NA.
  aside <- bk$status == "improper" & !is.finite(bk$t[, 1])
  expect_identical(bk$improper, sum(aside))
  expect_true(all(is.na(bk$t[aside, ])))
  kept <- bk$status == "improper" & is.finite(bk$t[, 1])
  expect_true(all(bk$t[kept, ] >= 0 & bk$t[kept, ] <= 1))
  # The improper point estimate has no maximal reliability, so its rows
  # carry no estimate and no interval, and nothing of it is bootstrapped.
  maximal <- k$index == "maximal_reliability"
  expect_true(all(is.na(as.matrix(k[maximal, c("estimate", "se",
                                              "ci_lower",
                                              "ci_upper")]))))
  expect_false(any(endsWith(colnames(bk$t), ":maximal_reliability")))
  expect_match(DMAR:::.dmar_tbl_footer(k)[["bifactor_interval"]],
               "that of the fit included, so its rows have no estimate")
})

test_that("a replication whose H is not defined is set aside under either rule", {
  skip_on_cran()
  # On these 150 cases the general loading of x1 is 1.06, so the general
  # factor's H is not defined. bifactor_indices() left that loading out
  # and gave the H of the other eight, a proportion, so under "keep" such
  # a replication was used: 31 of the 72 kept on the Holzinger and
  # Swineford model of the help page were of this kind.
  fit <- bf_quietly(bf_fit(bf_data(150, 7, near_lg, near_ls)))
  tab <- bf_quietly(bifactor_indices(fit))
  expect_true(is.na(tab$H[1]))
  cells <- DMAR:::.bifactor_cells(tab)
  cells <- cells[cells$index != "PUC", , drop = FALSE]
  pt <- lavaan::parTable(fit)
  pt <- pt[setdiff(names(pt), c("est", "se", "start"))]
  point_fn <- DMAR:::.bifactor_point_fn(
    model = pt, options = DMAR:::.bifactor_refit_options(fit, "model_implied"),
    ordered = character(0), categories = list(), general = "g",
    denominator = "model_implied", cells = cells)
  v <- point_fn(bf_pool(fit))
  expect_identical(attr(v, "dmar_reason"), "improper")
  expect_true(all(is.na(v)))
})

test_that("full information maximum likelihood resamples every row the fit used, with its pattern", {
  skip_on_cran()
  data <- bf_data(300, 113)
  data[1:40, "x1"] <- NA
  data[41:60, c("x4", "x8")] <- NA
  data[61, ] <- NA
  fit <- bf_quietly(bf_fit(data, missing = "ml", estimator = "MLR"))
  expect_identical(lavaan::lavInspect(fit, "nobs"), 299L)
  for (den in c("model_implied", "observed")) {
    ci <- bf_quietly(ci_bifactor_indices(fit, denominator = den, B = 12,
                                         seed = 113))
    expect_identical(attr(ci, "N"), 299L)
    expect_identical(attr(ci, "N_incomplete"), 60L)
    bi <- bifactor_indices(fit, denominator = den)
    expect_identical(ci$estimate[1:3], c(bi$ECV[1], bi$omega[1], bi$omega_H[1]))
    # A replication is a full information refit of its resample, with the
    # observed denominator taken from that refit's own saturated model.
    pool <- bf_pool(fit)
    expect_identical(nrow(pool), 299L)
    idx <- bf_indices(pool, B = 12, seed = 113)
    bs <- attr(ci, "bootstrap")
    for (r in c(2L, 11L)) {
      refit <- bf_quietly(bf_fit(pool[idx[r, ], ], missing = "ml",
                                 estimator = "MLR"))
      hand <- bf_quietly(bifactor_indices(refit, general = "g",
                                          denominator = den))
      cells <- DMAR:::.bifactor_cells(hand)
      cells <- cells[cells$index != "PUC", ]
      if (bs$status[r] == "proper") {
        expect_equal(unname(bs$t[r, ]),
                     DMAR:::.bifactor_cell_values(hand, cells),
                     tolerance = 1e-10)
      }
    }
    expect_match(DMAR:::.dmar_tbl_footer(ci)[["bifactor_interval"]],
                 "299 rows the fit used \\(60 of them incomplete")
  }
})

test_that("a fit with an auxiliary variable, and each replication, is judged on the model's own parameters", {
  skip_on_cran()
  # The help page's simulated items with an auxiliary variable that
  # correlates .8 with the general factor, given which x2 and x5 are
  # missing at random in 187 of the 600 rows. lavaan models it as a
  # saturated correlate whose residual covariances with the items make
  # theta indefinite, so its own check fails for the fit and for every
  # replication although the bifactor solution is proper.
  set.seed(113)
  n <- 600
  g <- rnorm(n)
  grp <- list(rnorm(n), rnorm(n), rnorm(n))
  X <- vapply(1:9, function(i) {
    0.5 * g + 0.5 * grp[[ceiling(i / 3)]] + sqrt(0.5) * rnorm(n)
  }, numeric(n))
  colnames(X) <- paste0("x", 1:9)
  data <- as.data.frame(X)
  data$a <- 0.8 * g + 0.6 * rnorm(n)
  drop <- runif(n) < stats::plogis(2 * data$a - 1.2)
  expect_identical(sum(drop), 187L)
  data[drop, c("x2", "x5")] <- NA
  aux_fit <- function(d) {
    suppressMessages(bf_quietly(bf_fit(d, missing = "ml", aux = "a")))
  }
  fit <- aux_fit(data)
  items <- paste0("x", 1:9)
  smallest <- function(m) min(eigen(m, symmetric = TRUE, only.values = TRUE)$values)
  expect_false(bf_quietly(lavaan::lavInspect(fit, "post.check")))
  expect_lt(smallest(lavaan::lavInspect(fit, "theta")), 0)
  expect_gt(smallest(lavaan::lavInspect(fit, "theta")[items, items]), 0)
  ci <- bf_quietly(ci_bifactor_indices(fit, B = 30, seed = 113))
  expect_false(attr(ci, "improper"))
  expect_identical(attr(ci, "N_incomplete"), 187L)
  # Each resample refit by hand: lavaan's own check fails every time, and a
  # replication is improper exactly when the model's own solution is, a
  # variance of the model below zero or a residual covariance matrix of
  # the items that is not positive definite.
  bs <- attr(ci, "bootstrap")
  pool <- bf_pool(fit)
  idx <- bf_indices(pool, B = 30, seed = 113)
  own_improper <- vapply(seq_len(30), function(r) {
    refit <- aux_fit(pool[idx[r, ], ])
    expect_false(bf_quietly(lavaan::lavInspect(refit, "post.check")))
    pe <- lavaan::parameterEstimates(refit)
    variances <- pe$est[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs != "a"]
    any(variances < 0) ||
      smallest(lavaan::lavInspect(refit, "theta")[items, items]) <= 0
  }, logical(1))
  expect_identical(bs$status == "improper", own_improper)
  expect_identical(bs$improper, sum(own_improper))
  expect_lt(bs$improper, 3L)
})

test_that("ordered items are refit with their thresholds, and a resample that loses a category is counted with its reason", {
  skip_on_cran()
  cont <- bf_data(300, 113)
  ord <- as.data.frame(lapply(cont, function(v) {
    as.integer(cut(v, c(-Inf, -0.8, 0, 0.8, Inf)))
  }))
  # One respondent alone in a fifth category of x2: about e^-1 of the
  # resamples do not draw that respondent.
  ord$x2[which.max(cont$x2)] <- 5L
  fit <- bf_fit(ord, ordered = TRUE)
  w <- testthat::capture_warnings(
    ci <- ci_bifactor_indices(fit, B = 15, seed = 113))
  bs <- attr(ci, "bootstrap")
  expect_identical(attr(ci, "metric"), "latent_response")
  expect_gt(bs$reasons[["error"]], 0L)
  expect_identical(sum(bs$status == "error"), bs$reasons[["error"]])
  expect_match(names(bs$error_messages),
               "drew no response in category 5 of the ordered item x2",
               all = FALSE)
  expect_match(w, paste0("raised an error \\(most often: the resample ",
                         "drew no response"), all = FALSE)
  # A replication that kept every category is a refit with the fit's
  # estimator (DWLS on the polychoric correlations).
  pool <- bf_pool(fit)
  idx <- bf_indices(pool, B = 15, seed = 113)
  r <- which(bs$status == "proper")[1]
  refit <- bf_fit(pool[idx[r, ], ], ordered = TRUE)
  hand <- bf_quietly(bifactor_indices(refit, general = "g"))
  cells <- DMAR:::.bifactor_cells(hand)
  cells <- cells[cells$index != "PUC", ]
  expect_equal(unname(bs$t[r, ]), DMAR:::.bifactor_cell_values(hand, cells),
               tolerance = 1e-8)
})

test_that("a refit that does not converge is named as such", {
  skip_on_cran()
  fit <- bf_fit(bf_data(300, 113))
  pt <- lavaan::parTable(fit)
  opts <- DMAR:::.bifactor_refit_options(fit, "model_implied")
  opts$control <- list(iter.max = 2L, eval.max = 4L)
  cells <- DMAR:::.bifactor_cells(bifactor_indices(fit))[-4, ]
  fn <- DMAR:::.bifactor_point_fn(
    model = pt[setdiff(names(pt), c("est", "se", "start"))], options = opts,
    ordered = character(0), categories = list(), general = "g",
    denominator = "model_implied", cells = cells)
  v <- fn(bf_pool(fit))
  expect_identical(attr(v, "dmar_reason"), "not_converged")
  expect_true(all(is.na(v)))
  expect_length(v, nrow(cells))
  code <- DMAR:::.boot_classify(v, "drop", nrow(cells))[nrow(cells) + 1L]
  expect_identical(code, 6)
})

test_that("a seed reproduces the intervals and leaves the caller's random stream as it was", {
  skip_on_cran()
  fit <- bf_fit(bf_data(300, 113))
  set.seed(42)
  before <- .Random.seed
  a <- bf_quietly(ci_bifactor_indices(fit, B = 10, seed = 113))
  expect_identical(.Random.seed, before)
  b <- bf_quietly(ci_bifactor_indices(fit, B = 10, seed = 113))
  expect_identical(a, b)
  unseeded <- bf_quietly(ci_bifactor_indices(fit, B = 10))
  expect_false(identical(.Random.seed, before))
})

test_that("ncpus spreads the replications without changing them", {
  skip_on_cran()
  skip_on_os("windows")
  fit <- bf_fit(bf_data(300, 113, near_lg, near_ls))
  a <- bf_quietly(ci_bifactor_indices(fit, B = 12, seed = 113))
  b <- bf_quietly(ci_bifactor_indices(fit, B = 12, seed = 113, ncpus = 2))
  ba <- attr(a, "bootstrap")
  bb <- attr(b, "bootstrap")
  expect_identical(bb$ncpus, 2L)
  ba$ncpus <- bb$ncpus <- NULL
  expect_identical(ba, bb)
  for (column in names(a)) expect_identical(a[[column]], b[[column]])
})

test_that("the point function a socket worker receives carries its values", {
  skip_on_cran()
  skip_on_os("windows")
  # A fork cluster sends functions serialized, as the socket cluster the
  # engine uses on Windows does, without needing an installed DMAR.
  fit <- bf_fit(bf_data(300, 113))
  pt <- lavaan::parTable(fit)
  cells <- DMAR:::.bifactor_cells(bifactor_indices(fit))[-4, ]
  fn <- DMAR:::.bifactor_point_fn(
    model = pt[setdiff(names(pt), c("est", "se", "start"))],
    options = DMAR:::.bifactor_refit_options(fit, "model_implied"),
    ordered = character(0), categories = list(), general = "g",
    denominator = "model_implied", cells = cells)
  expect_setequal(ls(environment(fn)),
                  c("model", "options", "ordered", "categories", "general",
                    "denominator", "cells"))
  st <- DMAR:::.boot_statistic(fn, "drop", nrow(cells))
  cl <- parallel::makeForkCluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  pool <- bf_pool(fit)
  b1 <- withr::with_seed(113, boot::boot(pool, st, R = 6, stype = "i"))
  b2 <- withr::with_seed(113, boot::boot(pool, st, R = 6, stype = "i",
                                         parallel = "snow", ncpus = 2L,
                                         cl = cl))
  expect_identical(b1$t, b2$t)
})

test_that("the BCa interval reads its acceleration from the jackknife when replications are fewer than cases", {
  skip_on_cran()
  fit <- bf_fit(bf_data(200, 113, rep(0.6, 9), rep(0.6, 9)))
  expect_false(isTRUE(attr(bifactor_indices(fit), "improper")))
  # Fewer replications used than cases: the jackknife, once for every index.
  j <- bf_quietly(ci_bifactor_indices(fit, ci_method = "bca", B = 20,
                                      seed = 113))
  bs <- attr(j, "bootstrap")
  expect_identical(bs$jackknife$attempted, 200L)
  expect_length(bs$acceleration, 20L)
  ok <- !is.na(bs$acceleration_method)
  expect_true(all(bs$acceleration_method[ok] == "jackknife"))
  expect_true(all(is.finite(bs$acceleration[ok])))
  expect_identical(attr(j, "ci_method_requested"), "bca")
  expect_true(all(attr(j, "ci_method")[ok] == "bca"))
  expect_match(DMAR:::.dmar_tbl_footer(j)[["bootstrap"]],
               "jackknife for the BCa acceleration: 200 attempted")
  p <- bf_quietly(ci_bifactor_indices(fit, B = 20, seed = 113))
  expect_identical(j$se, p$se)
  expect_false(identical(j$ci_lower, p$ci_lower))
})

test_that("what a bootstrap of cases cannot describe is refused", {
  data <- bf_data(300, 113)
  fit <- bf_fit(data)
  expect_error(ci_bifactor_indices(lm(x1 ~ x2, data = data)), "lavaan fit")
  expect_error(ci_bifactor_indices(fit, B = 1), "'B' must be")
  expect_error(ci_bifactor_indices(fit, B = 10.5), "'B' must be")
  expect_error(ci_bifactor_indices(fit, conf_level = 95), "'conf_level'")
  expect_error(ci_bifactor_indices(fit, ncpus = 0), "'ncpus' must be")
  expect_error(ci_bifactor_indices(fit, ci_method = "bootstrap_se"),
               "should be one of")
  # Correlated factors, with bifactor_indices()'s own message (the
  # oblique bifactor model is not identified, and lavaan says so on the
  # way).
  oblique <- bf_quietly(lavaan::cfa(bf_model, data = data, std.lv = TRUE))
  expect_error(bf_quietly(ci_bifactor_indices(oblique, B = 10)),
               "needs orthogonal general and group factors")
  # A multiple-group fit.
  grouped <- data
  grouped$grp <- rep(c("a", "b"), each = 150)
  mg <- bf_quietly(bf_fit(grouped, group = "grp"))
  expect_error(ci_bifactor_indices(mg, B = 10), "single-group fit")
  # A cluster variable.
  clustered <- data
  clustered$cl <- rep(1:30, each = 10)
  cf <- bf_quietly(bf_fit(clustered, cluster = "cl"))
  expect_error(ci_bifactor_indices(cf, B = 10), "cluster variable \\(cl\\)")
  # Sampling weights.
  weighted <- data
  weighted$w <- rep(c(0.5, 1.5), 150)
  wf <- bf_quietly(bf_fit(weighted, sampling.weights = "w"))
  expect_error(ci_bifactor_indices(wf, B = 10), "sampling weights")
  # A fit to a covariance matrix has no cases.
  sf <- bf_quietly(lavaan::cfa(bf_model, sample.cov = stats::cov(data),
                               sample.nobs = 300, orthogonal = TRUE,
                               std.lv = TRUE))
  expect_error(ci_bifactor_indices(sf, B = 10), "carries no raw data")
  # A fit that did not converge.
  nc <- bf_quietly(bf_fit(data, control = list(iter.max = 2L, eval.max = 4L)))
  expect_false(lavaan::lavInspect(nc, "converged"))
  expect_error(ci_bifactor_indices(nc, B = 10), "did not converge")
})

test_that("the printed table names its intervals, and the footer is the same for a report", {
  skip_on_cran()
  fit <- bf_fit(bf_data(300, 113))
  ci <- bf_quietly(ci_bifactor_indices(fit, B = 10, seed = 113))
  foot <- DMAR:::.dmar_tbl_footer(ci)
  expect_identical(names(foot), c("conf_level", "denominator",
                                  "bifactor_interval", "bifactor_bias",
                                  "bifactor_definitions",
                                  "bifactor_definitions", "bootstrap"))
  expect_match(foot[["bifactor_interval"]],
               "^Intervals: percentile bootstrap over the 300 rows the fit used")
  expect_match(foot[["bifactor_interval"]],
               "PUC depends only on which items share a group factor")
  # What ECV and omega mean for a group factor, and which H this is.
  expect_match(foot[5], "ECV is the general factor's share of the common")
  expect_match(foot[6], paste0("index of Rodriguez, Reise, and Haviland ",
                               "\\(2016\\).*maximal_reliability is the ",
                               "reliability the optimally weighted items"))
  printed <- utils::capture.output(print(ci))
  for (line in foot) expect_true(any(printed == line))
  expect_identical(DMAR:::.ci_bifactor_footer(bifactor_indices(fit)),
                   character(0))
  skip_if_not_installed("knitr")
  expect_match(paste(as_kable(ci, format = "pipe"), collapse = "\n"),
               "Intervals: percentile bootstrap over the 300 rows")
})

test_that("tidy() gives one row per index and glance() one row of counts and settings", {
  skip_on_cran()
  fit <- bf_fit(bf_data(300, 113, near_lg, near_ls))
  ci <- bf_quietly(ci_bifactor_indices(fit, B = 10, seed = 113))
  td <- generics::tidy(ci)
  expect_named(td, c("term", "estimate", "se", "ci_lower", "ci_upper",
                     "conf_level"))
  expect_identical(td$term, paste(ci$factor, ci$index, sep = ":"))
  expect_identical(td$ci_lower, ci$ci_lower)
  gl <- generics::glance(ci)
  expect_identical(nrow(gl), 1L)
  bs <- attr(ci, "bootstrap")
  expect_identical(gl$B, 10L)
  expect_identical(gl$B_used, bs$returned)
  expect_identical(gl$set_aside, bs$failed)
  expect_identical(gl$improper_set_aside, bs$improper)
  expect_identical(gl$not_converged, 0L)
  expect_identical(gl$n_factors, 4L)
  expect_identical(gl$improper_rule, "drop")
})

test_that("the bootstrap attribute records require_half = TRUE, as the help says", {
  skip_on_cran()
  fit <- bf_fit(bf_data(300, 113))
  ci <- bf_quietly(ci_bifactor_indices(fit, B = 10, seed = 113))
  bs <- attr(ci, "bootstrap")
  expect_true(bs$require_half)
  expect_identical(names(bs$target_returned), colnames(bs$t))
  expect_identical(names(bs$target_withheld), colnames(bs$t))
})


test_that("the help describes released history rather than an interim implementation", {
  path <- file.path("..", "..", "man", "ci_bifactor_indices.Rd")
  text <- if (file.exists(path)) {
    paste(readLines(path, warn = FALSE), collapse = " ")
  } else {
    rd <- tools::Rd_db("DMAR")[["ci_bifactor_indices.Rd"]]
    paste(unlist(rd), collapse = " ")
  }
  text <- gsub("[[:space:]]+", " ", text)
  expect_false(grepl("rose from 0.852 to 0.876", text, fixed = TRUE))
})
