# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

bs_test_items <- function() c("t1_visual_perception", "t2_cubes", "t4_lozenges",
  "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
bs_test_model <- function() paste(
  "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges",
  "textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
  sep = "\n")
bs_incomplete <- function() {
  d <- holzinger_swineford[bs_test_items()]
  set.seed(113)
  drop <- runif(nrow(d)) < 0.3
  d[drop, c("t2_cubes", "t7_sentence")] <- NA
  d[sample(nrow(d), 20), "t9_word_meaning"] <- NA
  d
}

test_that("with complete data the bootstrap is lavaan's Bollen-Stine bootstrap, resample for resample", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[bs_test_items()]
  fit <- lavaan::cfa(bs_test_model(), data = d, std.lv = TRUE)
  x <- sem_exact_fit_bootstrap(fit, B = 25, seed = 7)
  oracle <- suppressWarnings(lavaan::bootstrapLavaan(
    fit, R = 25, type = "bollen.stine", iseed = 7,
    FUN = function(f) lavaan::lavTech(f, "test")[[1L]]$stat))
  boot <- attr(x, "bootstrap")
  ok <- !is.na(oracle)
  expect_true(sum(ok) >= 20L)
  expect_equal(boot$chi_square[ok], as.numeric(oracle)[ok], tolerance = 1e-8)
  expect_identical(levels(boot$status), c("proper", "improper", "not_converged", "error"))
  v <- function(term) x$value[x$term == term]
  test <- lavaan::lavTech(fit, "test")[["standard"]]
  expect_equal(v("chi_square"), test$stat)
  expect_equal(v("df"), test$df)
  expect_equal(v("p_chi_square"), test$pvalue)
  used <- boot$chi_square[boot$status %in% c("proper", "improper")]
  expect_equal(v("p_bootstrap"), mean(used >= test$stat))
  expect_equal(v("bootstrap_mean"), mean(used))
  mc <- ci_proportion(sum(used >= test$stat), length(used), conf_level = 0.95)
  expect_equal(v("p_bootstrap_ci_lower"), mc$value[mc$term == "ci_lower"])
  expect_equal(v("p_bootstrap_ci_upper"), mc$value[mc$term == "ci_upper"])
  expect_equal(v("B"), 25)
  expect_equal(v("B_used") + v("B_failed"), 25)
  expect_equal(v("N"), nrow(d))
  expect_equal(v("n_groups"), 1)
  record <- attr(x, "transformation")[[1L]]
  expect_identical(record$method, "bollen_stine")
  expect_true(is.na(record$mean_residual))   # no mean structure to check
  expect_output(print(x), "Bollen and Stine")
  expect_output(print(x), "Wilson interval")
  expect_s3_class(x, "dmar_sem_fit")
})

test_that("the complete data transformation gives the model implied moments exactly", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[bs_test_items()]
  fit <- lavaan::cfa(bs_test_model(), data = d, std.lv = TRUE, meanstructure = TRUE)
  Sigma <- unclass(lavaan::lavInspect(fit, "cov.ov"))
  Mu <- lavaan::lavInspect(fit, "mean.ov")
  S <- lavaan::lavInspect(fit, "sampstat")$cov
  tr <- DMAR:::.bs_transform_group(as.matrix(d), Sigma, Mu, S, "common")
  Y <- tr$Y
  expect_equal(unname(crossprod(scale(Y, scale = FALSE)) / nrow(Y)), unname(Sigma),
               tolerance = 1e-10)
  expect_equal(unname(colMeans(Y)), unname(as.numeric(Mu)), tolerance = 1e-10)
  refit <- lavaan::cfa(bs_test_model(), data = as.data.frame(Y), std.lv = TRUE,
                       meanstructure = TRUE)
  expect_lt(lavaan::lavTech(refit, "test")[[1L]]$stat, 1e-6)
  # With a mean structure the exported function still matches lavaan's
  # bootstrap resample for resample, and an MLR fit contributes its ML
  # estimates and standard statistic.
  x <- sem_exact_fit_bootstrap(fit, B = 10, seed = 3)
  oracle <- suppressWarnings(lavaan::bootstrapLavaan(
    fit, R = 10, type = "bollen.stine", iseed = 3,
    FUN = function(f) lavaan::lavTech(f, "test")[[1L]]$stat))
  ok <- !is.na(oracle)
  expect_equal(attr(x, "bootstrap")$chi_square[ok], as.numeric(oracle)[ok], tolerance = 1e-8)
  expect_lt(attr(x, "transformation")[[1L]]$mean_residual, 1e-10)
  mlr <- lavaan::cfa(bs_test_model(), data = d, std.lv = TRUE, estimator = "MLR")
  y <- sem_exact_fit_bootstrap(mlr, B = 6, seed = 3)
  expect_equal(y$value[y$term == "chi_square"],
               lavaan::lavTech(mlr, "test")[["standard"]]$stat)
  expect_identical(attr(y, "sem_fit")$estimator, "MLR")
  expect_equal(attr(y, "sem_fit")$chi_square, y$value[y$term == "chi_square"])
  expect_identical(attr(y, "sem_fit")$transformation, "bollen_stine")
})

test_that("both incomplete data transformations make the saturated FIML estimates the model implied moments", {
  skip_if_not_installed("lavaan")
  d <- bs_incomplete()
  fit <- lavaan::cfa(bs_test_model(), data = d, missing = "ml", std.lv = TRUE)
  X <- as.matrix(lavaan::lavInspect(fit, "data"))
  Sigma <- unclass(lavaan::lavInspect(fit, "cov.ov"))
  Mu <- lavaan::lavInspect(fit, "mean.ov")
  S <- lavaan::lavInspect(fit, "sampstat")$cov
  for (method in c("common", "pattern")) {
    tr <- DMAR:::.bs_transform_group(X, Sigma, Mu, S, method)
    expect_identical(tr$method, method)
    expect_identical(is.na(tr$Y), is.na(X))
    expect_lt(tr$score, 1e-9 * nrow(X))
    expect_lt(tr$mean_residual, 1e-10)
    # cov_fiml() is an independent implementation of the saturated FIML
    # estimates (EM); on the transformed data they are the fitted moments.
    em <- cov_fiml(as.data.frame(tr$Y))
    expect_equal(unname(em$S[colnames(Sigma), colnames(Sigma)]), unname(Sigma),
                 tolerance = 1e-6)
    expect_equal(unname(em$mu[names(Mu)]), unname(as.numeric(Mu)), tolerance = 1e-8)
    refit <- lavaan::cfa(bs_test_model(), data = as.data.frame(tr$Y),
                         missing = "ml", std.lv = TRUE)
    expect_lt(lavaan::lavTech(refit, "test")[[1L]]$stat, 1e-5)
  }
  common <- DMAR:::.bs_transform_group(X, Sigma, Mu, S, "common")
  expect_gt(common$iterations, 0)
  expect_gt(common$patterns, 1)
})

test_that("a FIML fit is bootstrapped with its missing values, and the two transformations are selectable", {
  skip_if_not_installed("lavaan")
  d <- bs_incomplete()
  fit <- lavaan::cfa(bs_test_model(), data = d, missing = "ml", std.lv = TRUE)
  x <- sem_exact_fit_bootstrap(fit, B = 12, seed = 113)
  expect_identical(attr(x, "transformation")[[1L]]$method, "common")
  expect_equal(x$value[x$term == "N"], sum(lavaan::lavInspect(fit, "nobs")))
  expect_true(x$value[x$term == "p_bootstrap"] >= 0 && x$value[x$term == "p_bootstrap"] <= 1)
  expect_output(print(x), "under FIML")
  y <- sem_exact_fit_bootstrap(fit, B = 12, seed = 113, transformation = "pattern")
  expect_identical(attr(y, "transformation")[[1L]]$method, "pattern")
  expect_output(print(y), "pattern by pattern")
  # The same seed gives the same resamples whatever ncpus is.
  skip_on_os("windows")
  z <- sem_exact_fit_bootstrap(fit, B = 12, seed = 113, ncpus = 2L)
  expect_equal(attr(z, "bootstrap")$chi_square, attr(x, "bootstrap")$chi_square)
})

test_that("the pattern transformation refuses a pattern with too few rows", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[bs_test_items()]
  d[1, c("t1_visual_perception", "t2_cubes")] <- NA
  fit <- lavaan::cfa(bs_test_model(), data = d, missing = "ml", std.lv = TRUE)
  expect_error(sem_exact_fit_bootstrap(fit, B = 3, transformation = "pattern"),
               "more rows than observed variables")
  # The refusal is classed and carries the pattern (the first, by its
  # first row), its rows, its observed variables, and its group.
  e <- tryCatch(sem_exact_fit_bootstrap(fit, B = 3, transformation = "pattern"),
                dmar_refusal = function(e) e)
  expect_s3_class(e, c("dmar_pattern_too_small", "dmar_refusal", "error"))
  expect_identical(e$code, "pattern_too_small")
  expect_identical(c(e$pattern, e$n, e$p, e$group), c(1L, 1L, 4L, 1L))
  expect_identical(e$observed, bs_test_items()[3:6])
  expect_no_error(sem_exact_fit_bootstrap(fit, B = 3, seed = 1))
  # With several groups the message names the group and the condition
  # carries its index.
  hs <- holzinger_swineford[c(bs_test_items(), "school")]
  hs[c(1, 200), c("t1_visual_perception", "t2_cubes")] <- NA
  mg <- lavaan::cfa(bs_test_model(), data = hs, missing = "ml", std.lv = TRUE,
                    group = "school")
  e <- tryCatch(sem_exact_fit_bootstrap(mg, B = 3, transformation = "pattern"),
                dmar_refusal = function(e) e)
  expect_s3_class(e, "dmar_pattern_too_small")
  expect_match(conditionMessage(e), paste0("pattern 1 of group Pasteur has 1 rows on ",
                                           "4 variables"), fixed = TRUE)
  expect_identical(e$group, 1L)
})

test_that("a transformation undefined for the data is refused with its reason and group", {
  X <- as.matrix(holzinger_swineford[bs_test_items()])
  S <- stats::cov(X) * (nrow(X) - 1) / nrow(X)
  undefined <- function(expr, reason) {
    e <- tryCatch(expr, dmar_refusal = function(e) e)
    expect_s3_class(e, c("dmar_transformation_undefined", "dmar_refusal", "error"))
    expect_identical(e$reason, reason)
    invisible(e)
  }
  # A model implied covariance matrix that is not positive definite.
  Sigma <- S
  Sigma[1, 2] <- Sigma[2, 1] <- 1.5 * sqrt(S[1, 1] * S[2, 2])
  e <- undefined(DMAR:::.bs_transform_group(X, Sigma, NULL, S, group = 2L),
                 "not_positive_definite")
  expect_identical(e$matrix, "model implied covariance matrix")
  expect_identical(e$group, 2L)
  expect_match(conditionMessage(e), "is not positive definite, so the model based")
  # Incomplete data without a mean structure, and a saturated variance of 0.
  Xm <- X
  Xm[1:20, 2] <- NA
  undefined(DMAR:::.bs_transform_group(Xm, S, NULL, S), "no_mean_structure")
  S0 <- S
  S0[3, 3] <- 0
  e <- undefined(DMAR:::.bs_transform_group(Xm, S, colMeans(X), S0),
                 "nonpositive_variance")
  expect_identical(e$group, 1L)
  # A model implied covariance matrix that is singular on the variables a
  # pattern observes: u and v are the same variable in Sigma. The first
  # pattern (first row) observes u and w, whose submatrix is proper; the
  # second observes all three.
  set.seed(1)
  Xs <- matrix(rnorm(30), 10, 3, dimnames = list(NULL, c("u", "v", "w")))
  Xs[1, "v"] <- NA
  Sig <- matrix(c(1, 1, 0.5, 1, 1, 0.5, 0.5, 0.5, 1), 3,
                dimnames = list(colnames(Xs), colnames(Xs)))
  e <- undefined(DMAR:::.bs_transform_group(Xs, Sig, c(0, 0, 0), diag(3),
                                            group = 3L),
                 "singular_pattern")
  expect_identical(e$pattern, 2L)
  expect_identical(e$group, 3L)
  expect_match(conditionMessage(e), "observed in one missing data pattern is singular",
               fixed = TRUE)
})

test_that("the incomplete data transformations do not depend on the units, unlike their raw unit versions", {
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa(bs_test_model(), data = bs_incomplete(), missing = "ml",
                     std.lv = TRUE)
  X <- as.matrix(lavaan::lavInspect(fit, "data"))
  Sigma <- unclass(lavaan::lavInspect(fit, "cov.ov"))
  Mu <- as.numeric(lavaan::lavInspect(fit, "mean.ov"))
  S <- unclass(lavaan::lavInspect(fit, "sampstat")$cov)
  # The first item in tenths of its unit: its transformed values are ten
  # times the others', and nothing else changes.
  L <- c(10, rep(1, ncol(X) - 1L))
  for (method in c("common", "pattern")) {
    a <- DMAR:::.bs_transform_group(X, Sigma, Mu, S, method)$Y
    b <- DMAR:::.bs_transform_group(sweep(X, 2L, L, `*`), Sigma * tcrossprod(L),
                                    Mu * L, S * tcrossprod(L), method)$Y
    expect_equal(sweep(b, 2L, L, `/`), a, tolerance = 1e-8)
  }
  # The pattern transformation with the symmetric square roots taken in raw
  # units also makes the saturated FIML estimates the model implied moments,
  # but gives other values than the one solved on the standardized scale.
  sym_sqrt <- function(M, inverse = FALSE) {
    e <- eigen(M, symmetric = TRUE)
    e$vectors %*% ((if (inverse) 1 / sqrt(e$values) else sqrt(e$values)) *
                     t(e$vectors))
  }
  key <- apply(is.na(X), 1L, paste, collapse = "")
  raw <- X
  for (k in unique(key)) {
    r <- which(key == k)
    o <- !is.na(X[r[1L], ])
    xc <- scale(X[r, o, drop = FALSE], scale = FALSE)
    raw[r, o] <- sweep(xc %*% sym_sqrt(crossprod(xc) / length(r), inverse = TRUE) %*%
                         sym_sqrt(Sigma[o, o]), 2L, Mu[o], `+`)
  }
  standardized <- DMAR:::.bs_transform_group(X, Sigma, Mu, S, "pattern")$Y
  expect_lt(DMAR:::.bs_saturated_score(raw, Sigma, Mu, !is.na(X)), 1e-8)
  expect_gt(max(abs(raw - standardized), na.rm = TRUE), 0.1)
})

test_that("set aside replications are counted and too few of them leave the p-value NA", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[bs_test_items()]
  fit <- lavaan::cfa(bs_test_model(), data = d, std.lv = TRUE)
  original <- DMAR:::.bs_refit_one
  state <- new.env(parent = emptyenv()); state$n <- 0L
  fake <- function(rows, Y, group_var, group_labels, pt, model, options) {
    state$n <- state$n + 1L
    if (state$n %% 4L == 0L) return(c(statistic = NA_real_, status = 3))
    if (state$n %% 4L == 1L) return(c(statistic = 5, status = 2))
    original(rows, Y, group_var, group_labels, pt, model, options)
  }
  local_mocked_bindings(.bs_refit_one = fake)
  expect_warning(x <- sem_exact_fit_bootstrap(fit, B = 8, seed = 1),
                 "2 of 8 bootstrap samples were set aside")
  v <- function(term) x$value[x$term == term]
  expect_equal(v("B_used"), 6)
  expect_equal(v("B_improper"), 2)
  expect_equal(v("B_failed"), 2)
  expect_equal(sum(attr(x, "bootstrap")$status == "not_converged"), 2)
  state$n <- 0L
  expect_warning(y <- sem_exact_fit_bootstrap(fit, B = 8, seed = 1, improper = "drop"),
                 "gave an improper solution")
  expect_equal(y$value[y$term == "B_used"], 4)
  expect_output(print(y), "set aside")
  always_fail <- function(rows, Y, group_var, group_labels, pt, model, options) {
    c(statistic = NA_real_, status = 4)
  }
  local_mocked_bindings(.bs_refit_one = always_fail)
  # One warning says both that the samples were set aside and that no
  # p-value is formed (two warnings through the first versions of 1.1.0).
  w <- testthat::capture_warnings(
    z <- sem_exact_fit_bootstrap(fit, B = 4, seed = 1))
  expect_length(w, 1L)
  expect_match(w, "4 of 4 bootstrap samples were set aside")
  expect_match(w, "Fewer than half")
  expect_warning(sem_exact_fit_bootstrap(fit, B = 4, seed = 1),
                 class = "dmar_bootstrap_insufficient")
  expect_identical(attr(z, "bootstrap_status"), "insufficient_share")
  expect_true(is.na(z$value[z$term == "p_bootstrap"]))
  expect_true(is.na(z$value[z$term == "bootstrap_mean"]))
  expect_equal(z$value[z$term == "B_failed"], 4)
})

test_that("the bootstrap refuses what the transformation does not cover and restores the caller's RNG", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[bs_test_items()]
  fit <- lavaan::cfa(bs_test_model(), data = d, std.lv = TRUE)
  # B is held to the rule of ?dmar_bootstrap, whose refusal is classed
  # (the count it requires is lowered to 2 for this file).
  arg <- tryCatch(sem_exact_fit_bootstrap(fit, B = 1), error = function(e) e)
  expect_match(conditionMessage(arg), "at least 2")
  expect_s3_class(arg, "dmar_count_invalid")
  expect_error(sem_exact_fit_bootstrap(fit, B = 4, ncpus = 0), "ncpus")
  # A fit outside the scope of the transformation is refused with an error
  # of class dmar_exact_fit_bootstrap_scope whose reason names the case.
  scope <- function(expr, pattern, reason) {
    e <- tryCatch(expr, dmar_refusal = function(e) e)
    expect_s3_class(e, c("dmar_exact_fit_bootstrap_scope", "dmar_refusal", "error"))
    expect_identical(e$code, "exact_fit_bootstrap_scope")
    expect_match(conditionMessage(e), pattern)
    expect_identical(e$reason, reason)
    invisible(e)
  }
  summary_fit <- lavaan::cfa(bs_test_model(), sample.cov = cov(d), sample.nobs = nrow(d),
                             std.lv = TRUE)
  scope(sem_exact_fit_bootstrap(summary_fit, B = 3), "raw data", "summary_data")
  d2 <- d; d2[1:30, 2] <- NA
  two_stage <- lavaan::cfa(bs_test_model(), data = d2, missing = "two.stage", std.lv = TRUE)
  scope(sem_exact_fit_bootstrap(two_stage, B = 3), "two-stage", "two_stage")
  fixed_x <- lavaan::sem(paste("visual =~ t1_visual_perception + t2_cubes + t4_lozenges",
                               "visual ~ t9_word_meaning", sep = "\n"), data = d)
  e <- scope(sem_exact_fit_bootstrap(fixed_x, B = 3), "fixed.x = FALSE", "fixed_x")
  expect_identical(e$variables, "t9_word_meaning")
  expect_match(conditionMessage(e), "Refit with fixed.x = FALSE.", fixed = TRUE)
  ordered <- d
  for (j in seq_along(ordered)) ordered[[j]] <- as.integer(cut(ordered[[j]],
    unique(quantile(ordered[[j]], seq(0, 1, length.out = 5))), include.lowest = TRUE))
  categorical <- lavaan::cfa(bs_test_model(), data = ordered, ordered = TRUE, std.lv = TRUE)
  e <- scope(sem_exact_fit_bootstrap(categorical, B = 3), "ordered indicators", "ordered")
  expect_setequal(e$variables, bs_test_items())
  gls <- lavaan::cfa(bs_test_model(), data = d, std.lv = TRUE, estimator = "GLS")
  e <- scope(sem_exact_fit_bootstrap(gls, B = 3), "this fit used GLS", "estimator")
  expect_identical(e$estimator, "GLS")
  # Pairwise deletion leaves the missing values in the analysis data.
  pairwise <- lavaan::cfa(bs_test_model(), data = d2, std.lv = TRUE, missing = "pairwise")
  e <- scope(sem_exact_fit_bootstrap(pairwise, B = 3), "did not use FIML",
             "missing_without_fiml")
  expect_identical(e$group, 1L)
  listwise <- lavaan::cfa(bs_test_model(), data = d2, std.lv = TRUE)
  expect_warning(w <- sem_exact_fit_bootstrap(listwise, B = 3, seed = 1), "listwise deletion")
  expect_equal(w$value[w$term == "N"], nrow(d2) - 30)
  no_test <- lavaan::cfa(bs_test_model(), data = d, std.lv = TRUE, test = "none")
  scope(sem_exact_fit_bootstrap(no_test, B = 3), "no standard test statistic", "no_test")
  hs <- lavaan::HolzingerSwineford1939
  # A fixture for the refusal alone: lavaan remarks that the covariance
  # matrix of this clustered fit's estimates has an eigenvalue of -1e-16,
  # which is not under test here.
  clustered <- suppressWarnings(
    lavaan::cfa("visual =~ x1 + x2 + x3\ntextual =~ x4 + x5 + x6", data = hs,
                cluster = "school", estimator = "MLR"))
  e <- scope(sem_exact_fit_bootstrap(clustered, B = 3), "cluster variable", "cluster")
  expect_identical(e$cluster, "school")
  set.seed(7)
  hs$w <- runif(nrow(hs), 0.5, 1.5)
  weighted <- lavaan::cfa("visual =~ x1 + x2 + x3\ntextual =~ x4 + x5 + x6", data = hs,
                          sampling.weights = "w")
  scope(sem_exact_fit_bootstrap(weighted, B = 3), "sampling weights", "weights")
  conditional <- lavaan::sem("visual =~ x1 + x2 + x3\nvisual ~ x9", data = hs,
                             conditional.x = TRUE)
  e <- scope(sem_exact_fit_bootstrap(conditional, B = 3), "conditional.x", "fixed_x")
  expect_identical(e$variables, "x9")
  # lavaan does not take fixed.x = FALSE beside conditional.x = TRUE.
  expect_match(conditionMessage(e),
               "Refit with conditional.x = FALSE and fixed.x = FALSE.", fixed = TRUE)
  set.seed(42); expected <- runif(2)
  set.seed(42); invisible(sem_exact_fit_bootstrap(fit, B = 3, seed = 9))
  expect_equal(runif(2), expected)
  skip_on_cran()
  two_level <- lavaan::sem("level: 1
                              fw =~ y1 + y2 + y3
                            level: 2
                              fb =~ y1 + y2 + y3",
                           data = lavaan::Demo.twolevel, cluster = "cluster")
  scope(sem_exact_fit_bootstrap(two_level, B = 3), "two-level model", "two_level")
})

test_that("several groups are transformed and resampled group by group", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939
  fit <- lavaan::cfa("visual =~ x1 + x2 + x3\ntextual =~ x4 + x5 + x6", data = hs,
                     group = "school", std.lv = TRUE)
  x <- sem_exact_fit_bootstrap(fit, B = 10, seed = 2)
  expect_equal(x$value[x$term == "n_groups"], 2)
  expect_identical(names(attr(x, "transformation")), lavaan::lavInspect(fit, "group.label"))
  oracle <- suppressWarnings(lavaan::bootstrapLavaan(
    fit, R = 10, type = "bollen.stine", iseed = 2,
    FUN = function(f) lavaan::lavTech(f, "test")[[1L]]$stat))
  mine <- attr(x, "bootstrap")$chi_square
  ok <- !is.na(oracle)
  expect_equal(mine[ok], as.numeric(oracle)[ok], tolerance = 1e-8)
  # Several groups with incomplete data and a constraint across groups:
  # each group is transformed on its own and refitted by FIML.
  set.seed(5)
  hs2 <- hs
  hs2[sample(nrow(hs2), 50), "x2"] <- NA
  hs2[sample(nrow(hs2), 40), "x5"] <- NA
  fiml <- lavaan::cfa("visual =~ x1 + x2 + x3\ntextual =~ x4 + x5 + x6", data = hs2,
                      group = "school", group.equal = "loadings", missing = "ml")
  y <- sem_exact_fit_bootstrap(fiml, B = 8, seed = 5)
  rec <- attr(y, "transformation")
  expect_identical(unname(vapply(rec, `[[`, character(1), "method")), c("common", "common"))
  expect_true(all(vapply(rec, `[[`, numeric(1), "score") < 1e-9 * nrow(hs2)))
  expect_true(all(vapply(rec, `[[`, numeric(1), "mean_residual") < 1e-10))
  expect_equal(y$value[y$term == "N"], sum(lavaan::lavInspect(fiml, "nobs")))
  expect_equal(y$value[y$term == "B_used"] + y$value[y$term == "B_failed"], 8)
})

# Six of the tests with the cube and sentence scores missing at random given
# the general information score, the auxiliary variable of the fits below.
bs_aux_data <- function() {
  d <- holzinger_swineford[c(bs_test_items(), "t5_general_information")]
  set.seed(113)
  z <- as.numeric(scale(d$t5_general_information))
  drop <- runif(nrow(d)) < plogis(1.5 * z - 1)
  d[drop, c("t2_cubes", "t7_sentence")] <- NA
  d
}
bs_aux_fit <- function(data, ...) {
  # lavaan notes that it adds the auxiliary variable as a saturated
  # correlate, and its own post-estimation check, at the end of the fit,
  # warns that theta, the auxiliary variable's residual covariances
  # included, is not positive definite: the case these tests are about.
  suppressMessages(suppressWarnings(
    lavaan::cfa(bs_test_model(), data = data, missing = "ml", std.lv = TRUE,
                aux = "t5_general_information", ...)))
}

test_that("a FIML fit with an auxiliary variable is judged on the model's own parameters", {
  skip_if_not_installed("lavaan")
  fit <- bs_aux_fit(bs_aux_data())
  # lavaan's check fails on the auxiliary variable's residual covariances,
  # which amount to correlations above 1 with the verbal tests, although
  # the model's own solution is proper.
  expect_false(suppressWarnings(lavaan::lavInspect(fit, "post.check")))
  r <- stats::cov2cor(lavaan::lavInspect(fit, "theta"))["t5_general_information", ]
  expect_gt(min(r[c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")]), 1)
  expect_identical(DMAR:::.sem_aux_names(fit), "t5_general_information")
  expect_identical(DMAR:::.sem_inadmissible_conditions(fit)$reason, "theta")
  ok <- DMAR:::.sem_admissible(fit)
  expect_true(ok)
  expect_identical(nrow(attr(ok, "problems")), 0L)
  expect_identical(attr(ok, "auxiliary"), "t5_general_information")
  # Both functions accept it and report lavaan's statistic.
  stat <- lavaan::lavTech(fit, "test")[["standard"]]$stat
  x <- sem_exact_fit_bootstrap(fit, B = 10, seed = 113)
  expect_equal(x$value[x$term == "chi_square"], stat)
  expect_equal(x$value[x$term == "B_used"] + x$value[x$term == "B_failed"], 10)
  m <- sem_fit_measures(fit)
  expect_equal(m$value[m$term == "chi_square"], stat)
  expect_equal(m$value[m$term == "chi_square"], unname(lavaan::fitMeasures(fit, "chisq")))
  # A two-stage fit's auxiliary variables enter its first stage only and
  # are not among the model's variables, so it has none to set aside.
  two_stage <- suppressMessages(lavaan::cfa(bs_test_model(), data = bs_aux_data(),
                                            missing = "two.stage", std.lv = TRUE,
                                            aux = "t5_general_information"))
  expect_identical(DMAR:::.sem_aux_names(two_stage), character(0))
})

test_that("a Heywood case with an auxiliary variable is still refused, with its reason", {
  skip_if_not_installed("lavaan")
  # Three items whose correlations no factor reproduces with proper
  # variances (the first loading is 1.46 in the population), and an
  # auxiliary variable related to all three.
  R <- diag(4)
  R[1, 2:3] <- R[2:3, 1] <- 0.8
  R[2, 3] <- R[3, 2] <- 0.3
  R[4, 1:3] <- R[1:3, 4] <- c(0.5, 0.4, 0.4)
  set.seed(113)
  z <- scale(matrix(rnorm(300 * 4), 300, 4))
  z <- z %*% solve(chol(stats::cov(z))) %*% chol(R)
  d <- as.data.frame(10 * z + 50)
  names(d) <- c("y1", "y2", "y3", "a")
  d$y2[runif(300) < stats::plogis((d$a - 60) / 5)] <- NA
  fit <- suppressMessages(suppressWarnings(
    lavaan::cfa("f =~ y1 + y2 + y3", data = d, missing = "ml", aux = "a")))
  expect_true(lavaan::lavInspect(fit, "converged"))
  ok <- DMAR:::.sem_admissible(fit)
  expect_false(ok)
  expect_identical(attr(ok, "problems")$variable, "y1")
  e <- tryCatch(sem_fit_measures(fit), dmar_refusal = function(e) e)
  expect_s3_class(e, c("dmar_inadmissible", "dmar_refusal", "error"))
  expect_identical(e$reason, "negative_variance")
  expect_identical(e$variables, "y1")
  expect_identical(e$auxiliary, "a")
  expect_lt(e$problems$value, 0)
  expect_match(conditionMessage(e), paste0(
    "post-estimation check of the model's own parameters, with the auxiliary ",
    "variable a set aside, finds a negative variance estimate for y1"), fixed = TRUE)
  expect_error(sem_exact_fit_bootstrap(fit, B = 3), class = "dmar_inadmissible")
})

test_that("a multiple group fit with an auxiliary variable is judged group by group", {
  skip_if_not_installed("lavaan")
  d <- cbind(bs_aux_data(), school = holzinger_swineford$school)
  fit <- bs_aux_fit(d, group = "school")
  expect_false(suppressWarnings(lavaan::lavInspect(fit, "post.check")))
  kept <- DMAR:::.sem_inadmissible_conditions(fit)
  expect_identical(kept$reason, c("theta", "theta"))
  expect_identical(kept$group, lavaan::lavInspect(fit, "group.label"))
  expect_true(DMAR:::.sem_admissible(fit))
  x <- sem_exact_fit_bootstrap(fit, B = 6, seed = 113)
  expect_equal(x$value[x$term == "chi_square"], lavaan::lavTech(fit, "test")[["standard"]]$stat)
  expect_equal(x$value[x$term == "n_groups"], 2)
  expect_equal(x$value[x$term == "B_used"] + x$value[x$term == "B_failed"], 6)
})

test_that("with an auxiliary variable a resample is improper only when the model's own solution is", {
  skip_on_cran()  # 100 refits by DMAR and 100 more to check them
  skip_if_not_installed("lavaan")
  items <- c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning",
             "t10_addition", "t12_counting_groups_of_dots",
             "t13_straight_and_curved_capitals")
  d <- holzinger_swineford[c(items, "age_years")]
  set.seed(20261106)
  z <- as.numeric(scale(d$t6_paragraph_comprehension))
  drop <- runif(nrow(d)) < plogis(-1.4 + 0.9 * z)
  d[drop, c("t7_sentence", "t12_counting_groups_of_dots")] <- NA
  d[runif(nrow(d)) < 0.15, "t13_straight_and_curved_capitals"] <- NA
  model <- paste(
    "textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    "speed =~ t10_addition + t12_counting_groups_of_dots + t13_straight_and_curved_capitals",
    sep = "\n")
  # lavaan warns that some observed variances are more than 1000 times
  # others, which is not under test here.
  fit <- suppressMessages(suppressWarnings(
    lavaan::cfa(model, data = d, missing = "ml", std.lv = TRUE, aux = "age_years")))
  B <- 100
  x <- sem_exact_fit_bootstrap(fit, B = B, seed = 113)
  boot <- attr(x, "bootstrap")
  # The resamples rebuilt: the common transformation of the fit's data and
  # the draws that follow the seed.
  ov <- colnames(unclass(lavaan::lavInspect(fit, "cov.ov")))
  X <- as.matrix(lavaan::lavInspect(fit, "data"))
  colnames(X) <- lavaan::lavNames(fit, "ov")
  Y <- DMAR:::.bs_transform_group(
    X[, ov], unclass(lavaan::lavInspect(fit, "cov.ov"))[ov, ov],
    as.numeric(lavaan::lavInspect(fit, "mean.ov")[ov]),
    unclass(lavaan::lavInspect(fit, "sampstat")$cov)[ov, ov])$Y
  set.seed(113)
  rows <- lapply(seq_len(B), function(b) sample.int(nrow(Y), replace = TRUE))
  # Each resample refitted and judged here with eigen(): improper when a
  # variance of the model is negative, or the residual covariance matrix of
  # the six tests, psi, or the latent covariance matrix has a negative
  # eigenvalue; and lavaan's verdict beside it.
  smallest <- function(m) min(eigen(m, symmetric = TRUE, only.values = TRUE)$values)
  judged <- t(vapply(rows, function(r) {
    f <- suppressMessages(suppressWarnings(lavaan::cfa(
      model, data = as.data.frame(Y[r, ]), missing = "ml", std.lv = TRUE,
      aux = "age_years", se = "none", baseline = FALSE, start = fit)))
    est <- lavaan::lavInspect(f, "est")
    theta <- est$theta[items, items]
    c(statistic = lavaan::lavTech(f, "test")[["standard"]]$stat,
      own = any(diag(theta) < 0) || any(diag(est$psi) < 0) || smallest(theta) < 0 ||
        smallest(est$psi) < 0 || smallest(lavaan::lavInspect(f, "cov.lv")) < 0,
      lavaan = !isTRUE(suppressWarnings(lavaan::lavInspect(f, "post.check"))))
  }, numeric(3)))
  ok <- !is.na(boot$chi_square)
  expect_true(all(ok))
  expect_equal(boot$chi_square, unname(judged[, "statistic"]), tolerance = 1e-6)
  expect_identical(boot$status == "improper", unname(judged[, "own"] == 1))
  expect_equal(x$value[x$term == "B_improper"], sum(judged[, "own"]))
  # lavaan's own check would count many more, almost all of them for the
  # residual covariances of the auxiliary variable alone.
  expect_gt(sum(judged[, "lavaan"]), 5 * sum(judged[, "own"]))
  expect_true(all(judged[judged[, "own"] == 1, "lavaan"] == 1))
})
