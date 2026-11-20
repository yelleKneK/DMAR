# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# Tests for the progress callback and cooperative cancellation of the
# bootstrap engines, added in 1.1.0: .bootstrap_ci() in
# R/reliability_internals.R, and the bootstrap then impute loop of
# categorical omega (.boot_mi() in R/multiple_imputation_internals.R),
# through every exported function that runs them. A callback is called
# with the replications attempted so far and B, after each replication
# in this process and after each block on several processes; FALSE
# stops the bootstrap with an error of class dmar_bootstrap_cancelled, an
# error inside the callback stops it with one of class
# dmar_progress_error, and nothing it does changes a result. The engine
# tests use coefficient alpha from the covariance matrix as the point
# function, which is cheap.

skip_if_not_installed("boot")

set.seed(113)
N <- 60
J <- 6
eta <- rnorm(N)
items <- outer(eta, rep(0.7, J)) +
  matrix(rnorm(N * J), N, J) * sqrt(1 - 0.7^2)
colnames(items) <- paste0("y", seq_len(J))

alpha_fn <- function(d) DMAR:::.alpha_from_S(stats::cov(d))

# A callback that records every call in an environment of its own and
# returns FALSE on call number `stop_at` (TRUE before it).
recorder <- function(stop_at = Inf) {
  calls <- new.env(parent = emptyenv())
  calls$done <- integer()
  calls$total <- integer()
  calls$callback <- function(done, total) {
    calls$done <- c(calls$done, done)
    calls$total <- c(calls$total, total)
    length(calls$done) < stop_at
  }
  calls
}

engine <- function(...) DMAR:::.bootstrap_ci(...)


test_that("in this process the callback is called after every replication, with done counting up to B", {
  r <- recorder()
  engine(items, alpha_fn, B = 40, conf_level = 0.95, kind = "percentile",
         seed = 113, progress = r$callback)
  expect_identical(r$done, 1:40)
  expect_identical(r$total, rep(40L, 40))
  # B given as a double reaches the callback as an integer.
  r <- recorder()
  engine(items, alpha_fn, B = 5, conf_level = 0.5, kind = "percentile",
         seed = 113, progress = r$callback)
  expect_identical(r$done, 1:5)
  expect_type(r$total, "integer")
})

test_that("the callback changes no interval, count, draw, or warning", {
  # Errors and improper estimates together, under both regimes of the BCa
  # acceleration (regression with B >= n, the jackknife with B < n).
  fn <- function(d) {
    if (!any(d[, 1] == items[2, 1])) stop("row 2 is missing from this resample")
    v <- alpha_fn(d)
    if (!any(d[, 1] == items[1, 1])) attr(v, "dmar_reason") <- "improper"
    v
  }
  for (B in c(120L, 50L)) {
    for (kind in c("percentile", "bca")) {
      w0 <- testthat::capture_warnings(
        r0 <- engine(items, fn, B = B, conf_level = 0.9, kind = kind,
                     seed = 113))
      r <- recorder()
      w1 <- testthat::capture_warnings(
        r1 <- engine(items, fn, B = B, conf_level = 0.9, kind = kind,
                     seed = 113, progress = r$callback))
      expect_identical(r1, r0)
      expect_identical(w1, w0)
      expect_gt(r0$bootstrap$failed, 0L)
    }
  }
  expect_identical(r0$bootstrap$acceleration_method, "jackknife")
})

test_that("while the jackknife of a BCa interval runs, the callback is called after each refit with done and total both B", {
  r <- recorder()
  out <- engine(items, alpha_fn, B = 40, conf_level = 0.8, kind = "bca",
                seed = 113, progress = r$callback)
  expect_identical(out$bootstrap$acceleration_method, "jackknife")
  expect_identical(r$done, c(1:40, rep(40L, N)))
  expect_false(is.unsorted(r$done))
  expect_identical(r$total, rep(40L, 40 + N))
  # With B >= n the acceleration comes from the regression, and the
  # callback is called B times.
  r <- recorder()
  out <- engine(items, alpha_fn, B = 80, conf_level = 0.8, kind = "bca",
                seed = 113, progress = r$callback)
  expect_identical(out$bootstrap$acceleration_method, "regression")
  expect_identical(r$done, 1:80)
})

test_that("a FALSE stops the bootstrap before the next replication, with the class, done, and total", {
  fits <- new.env(parent = emptyenv())
  fits$n <- 0L
  counted <- function(d) {
    fits$n <- fits$n + 1L
    alpha_fn(d)
  }
  r <- recorder(stop_at = 7)
  err <- expect_error(
    engine(items, counted, B = 40, conf_level = 0.95, kind = "percentile",
           seed = 113, progress = r$callback),
    class = "dmar_bootstrap_cancelled")
  expect_s3_class(err, "error")
  expect_identical(err$done, 7L)
  expect_identical(err$total, 40L)
  expect_identical(conditionMessage(err), paste0(
    "The bootstrap was cancelled after 7 of 40 replications: the progress ",
    "callback returned FALSE."))
  expect_null(conditionCall(err))
  # The estimate on the data and the seven replications, and no more.
  expect_identical(fits$n, 8L)
  expect_identical(r$done, 1:7)
  # A FALSE on the last call stops it too, before the interval is formed.
  r <- recorder(stop_at = 40)
  err <- expect_error(
    engine(items, alpha_fn, B = 40, conf_level = 0.95, kind = "percentile",
           seed = 113, progress = r$callback),
    class = "dmar_bootstrap_cancelled")
  expect_identical(err$done, 40L)
  # And during the jackknife, after the B replications.
  r <- recorder(stop_at = 43)
  err <- expect_error(
    engine(items, alpha_fn, B = 40, conf_level = 0.95, kind = "bca",
           seed = 113, progress = r$callback),
    class = "dmar_bootstrap_cancelled")
  expect_identical(c(err$done, err$total), c(40L, 40L))
  expect_length(r$done, 43L)
})

test_that("only FALSE stops the bootstrap; any other value is ignored", {
  ref <- engine(items, alpha_fn, B = 30, conf_level = 0.8,
                kind = "percentile", seed = 113)
  for (value in list(TRUE, NULL, NA, "FALSE", 0L, c(FALSE, FALSE),
                     c(cancel = FALSE), list(FALSE))) {
    out <- engine(items, alpha_fn, B = 30, conf_level = 0.8,
                  kind = "percentile", seed = 113,
                  progress = function(done, total) value)
    expect_identical(out, ref)
  }
})

test_that("an error inside the callback stops the call with an error of class dmar_progress_error, the callback's condition its parent", {
  err <- expect_error(
    engine(items, alpha_fn, B = 20, conf_level = 0.95, kind = "percentile",
           seed = 113, progress = function(done, total) {
             if (done == 3L) stop("the display is gone")
             TRUE
           }),
    class = "dmar_progress_error")
  expect_identical(class(err), c("dmar_progress_error", "error", "condition"))
  expect_identical(conditionMessage(err), paste0(
    "The progress callback failed after 3 of 20 bootstrap replications: ",
    "the display is gone"))
  expect_null(conditionCall(err))
  expect_identical(c(err$done, err$total), c(3L, 20L))
  expect_identical(conditionMessage(err$parent), "the display is gone")
  # A condition of the application's own keeps its class behind it, as the
  # silent error a Shiny req() signals does.
  silent <- structure(
    class = c("shiny.silent.error", "validation", "error", "condition"),
    list(message = "", call = NULL))
  err <- expect_error(
    engine(items, alpha_fn, B = 20, conf_level = 0.95, kind = "percentile",
           seed = 113, progress = function(done, total) stop(silent)),
    class = "dmar_progress_error")
  expect_false(inherits(err, "shiny.silent.error"))
  expect_identical(err$parent, silent)
  expect_identical(c(err$done, err$total), c(1L, 20L))
  # During the jackknife of a BCa interval, done and total are both B.
  r <- recorder()
  err <- expect_error(
    engine(items, alpha_fn, B = 40, conf_level = 0.8, kind = "bca",
           seed = 113, progress = function(done, total) {
             r$callback(done, total)
             if (length(r$done) > 41L) stop("the session ended")
             TRUE
           }),
    class = "dmar_progress_error")
  expect_identical(c(err$done, err$total), c(40L, 40L))
  expect_length(r$done, 42L)
  # Through an exported function as well.
  err <- expect_error(
    reliability_alpha(items, ci_method = "percentile", B = 20, seed = 113,
                      progress = function(done, total) stop("no session")),
    class = "dmar_progress_error")
  expect_identical(conditionMessage(err), paste0(
    "The progress callback failed after 1 of 20 bootstrap replications: ",
    "no session"))
})

test_that("a callback that draws random numbers changes neither the result nor the caller's stream", {
  noisy <- function(done, total) {
    stats::runif(3)
    TRUE
  }
  a0 <- reliability_alpha(items, ci_method = "bca", B = 100, seed = 113,
                          conf_level = 0.8)
  a1 <- reliability_alpha(items, ci_method = "bca", B = 100, seed = 113,
                          conf_level = 0.8, progress = noisy)
  expect_identical(a1, a0)
  # With seed = NULL the resamples come from the caller's stream, which
  # continues after the call as it would without the callback.
  set.seed(9)
  b0 <- reliability_alpha(items, ci_method = "percentile", B = 100)
  u0 <- stats::runif(2)
  set.seed(9)
  b1 <- reliability_alpha(items, ci_method = "percentile", B = 100,
                          progress = noisy)
  u1 <- stats::runif(2)
  expect_identical(b1, b0)
  expect_identical(u1, u0)
})

test_that("progress is validated, and with NULL the bootstrap says nothing", {
  for (bad in list("yes", 1, TRUE, list(function(d, t) TRUE))) {
    expect_error(reliability_alpha(items, ci_method = "percentile", B = 20,
                                   progress = bad),
                 "'progress' must be NULL or a function", fixed = TRUE)
  }
  binary <- (items > 0) * 1
  expect_error(reliability_kr20(binary, progress = "yes"),
               "'progress' must be NULL or a function", fixed = TRUE)
  expect_error(reliability(items, type = "alpha", progress = 2),
               "'progress' must be NULL or a function", fixed = TRUE)
  expect_error(DMAR:::.bootstrap_ci(items, alpha_fn, B = 20,
                                    conf_level = 0.95, kind = "percentile",
                                    seed = 113, progress = "yes"),
               "'progress' must be NULL or a function", fixed = TRUE)
  expect_silent(reliability_alpha(items, ci_method = "percentile", B = 50,
                                  seed = 113, progress = NULL))
  # A callback is not called when no bootstrap runs.
  r <- recorder()
  reliability_alpha(items, progress = r$callback)
  expect_length(r$done, 0L)
})

test_that("a block on several processes is about a fiftieth of B, a multiple of ncpus", {
  block <- DMAR:::.boot_block_size
  expect_identical(block(100, 1L), 1L)
  expect_identical(block(10000, 1L), 1L)
  expect_identical(block(100, 2L), 2L)
  expect_identical(block(250, 2L), 6L)
  expect_identical(block(1000, 4L), 20L)
  expect_identical(block(10000, 4L), 200L)
  expect_identical(block(3, 4L), 4L)
})

test_that("on forked workers the callback is called after each block, and the result is the serial one", {
  skip_on_os("windows")   # forked workers; the socket path is tested below
  skip_on_cran()
  fn <- function(d) {
    if (!any(d[, 1] == items[2, 1])) stop("row 2 is missing from this resample")
    alpha_fn(d)
  }
  ref <- suppressWarnings(engine(items, fn, B = 250, conf_level = 0.9,
                                 kind = "bca", seed = 113))
  r <- recorder()
  out <- suppressWarnings(engine(items, fn, B = 250, conf_level = 0.9,
                                 kind = "bca", seed = 113, ncpus = 2L,
                                 progress = r$callback))
  expect_identical(out[setdiff(names(out), "bootstrap")],
                   ref[setdiff(names(ref), "bootstrap")])
  expect_identical(out$bootstrap[setdiff(names(out$bootstrap), "ncpus")],
                   ref$bootstrap[setdiff(names(ref$bootstrap), "ncpus")])
  expect_identical(out$bootstrap$ncpus, 2L)
  # Blocks of 2 * ceiling(250 / 100) = 6 replications, the last of 4; the
  # acceleration came from the regression (B >= n), so no jackknife ran.
  expect_identical(out$bootstrap$acceleration_method, "regression")
  expect_identical(r$done, c(seq(6L, 246L, by = 6L), 250L))
  expect_identical(r$total, rep(250L, 42))
  # With the jackknife: blocks of 2 * ceiling(60 / 100) = 2 refits.
  r <- recorder()
  out <- engine(items, alpha_fn, B = 40, conf_level = 0.9, kind = "bca",
                seed = 113, ncpus = 2L, progress = r$callback)
  expect_identical(out$bootstrap$acceleration_method, "jackknife")
  expect_identical(r$done, c(seq(2L, 40L, by = 2L), rep(40L, N / 2)))
  ref <- engine(items, alpha_fn, B = 40, conf_level = 0.9, kind = "bca",
                seed = 113, ncpus = 2L)
  expect_identical(out, ref)
})

test_that("on forked workers a FALSE stops the bootstrap at the end of a block", {
  skip_on_os("windows")
  skip_on_cran()
  # Blocks of 2 * ceiling(250 / 100) = 6 replications: a FALSE on the
  # third call stops after 18.
  r <- recorder(stop_at = 3)
  err <- expect_error(
    engine(items, alpha_fn, B = 250, conf_level = 0.95, kind = "percentile",
           seed = 113, ncpus = 2L, progress = r$callback),
    class = "dmar_bootstrap_cancelled")
  expect_identical(c(err$done, err$total), c(18L, 250L))
  expect_identical(r$done, c(6L, 12L, 18L))
  expect_match(conditionMessage(err), "cancelled after 18 of 250 replications",
               fixed = TRUE)
})

test_that("the functions sent to socket workers with a callback carry their arguments as values", {
  skip_on_os("windows")
  skip_on_cran()
  # A fork cluster sends each function to its workers serialized, as the
  # socket cluster the engine uses on Windows does, without needing an
  # installed DMAR.
  cl <- parallel::makeForkCluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  st <- DMAR:::.boot_statistic(alpha_fn, "keep")
  on_rows <- DMAR:::.boot_on_rows(st, items)
  expect_setequal(ls(environment(on_rows)), c("statistic", "data"))
  r <- recorder()
  b1 <- withr::with_seed(113, DMAR:::.boot_with_progress(
    items, st, B = 30L, m = 1L, par_type = "snow", ncpus = 2L, cl = cl,
    progress = r$callback))
  b0 <- withr::with_seed(113, boot::boot(items, st, R = 30, stype = "i"))
  expect_identical(b1$t, b0$t)
  expect_identical(b1$t0, b0$t0)
  expect_identical(b1$seed, b0$seed)
  expect_identical(r$done, seq(2L, 30L, by = 2L))
  r <- recorder()
  j <- DMAR:::.boot_jackknife(alpha_fn, items, "keep", "snow", 2L, cl,
                              progress = r$callback, B = 30L)
  expect_identical(j, DMAR:::.boot_jackknife(alpha_fn, items, "keep", "no",
                                             1L, NULL))
  expect_identical(r$done, rep(30L, N / 2))
})

test_that("on a socket cluster the function is sent to each worker once, and a block sends only its elements", {
  skip_on_os("windows")
  skip_on_cran()
  cl <- parallel::makeForkCluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # The function counts its calls in an environment of its own. Sent to a
  # worker once, a copy counts on across the blocks; sent again with every
  # block, each copy would start from zero and every count would be 1.
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  fun <- function(x) {
    counter$n <- counter$n + 1L
    c(x, counter$n, Sys.getpid())
  }
  r <- recorder()
  res <- DMAR:::.boot_lapply(seq_len(40L), fun, "snow", 2L, cl,
                             progress = r$callback)
  # Blocks of 2 * ceiling(40 / 100) = 2, one element to each worker.
  expect_identical(r$done, seq(2L, 40L, by = 2L))
  res <- do.call(rbind, res)
  expect_identical(res[, 1], seq_len(40L))
  by_worker <- split(res[, 2], res[, 3])
  expect_length(by_worker, 2L)
  for (counts in by_worker) expect_identical(counts, seq_len(20L))
  # Nothing ran, and nothing was staged, in this process.
  expect_identical(counter$n, 0L)
  expect_false(exists("fun", envir = DMAR:::.boot_worker, inherits = FALSE))
  # Without a callback the one parLapply() call is unchanged.
  res <- DMAR:::.boot_lapply(seq_len(40L), fun, "snow", 2L, cl)
  expect_identical(do.call(rbind, res)[, 1], seq_len(40L))
})

test_that("every exported bootstrap takes progress, and a seed gives the same result with or without it", {
  a0 <- reliability_alpha(items, ci_method = "percentile", B = 30, seed = 113,
                          conf_level = 0.8)
  r <- recorder()
  a1 <- reliability_alpha(items, ci_method = "percentile", B = 30, seed = 113,
                          conf_level = 0.8, progress = r$callback)
  expect_identical(a1, a0)
  expect_identical(r$done, 1:30)
  binary <- (items > 0) * 1
  k0 <- reliability_kr20(binary, ci_method = "bca", B = 30, seed = 113,
                         conf_level = 0.8)
  r <- recorder()
  k1 <- reliability_kr20(binary, ci_method = "bca", B = 30, seed = 113,
                         conf_level = 0.8, progress = r$callback)
  expect_identical(k1, k0)
  expect_identical(r$done, c(1:30, rep(30L, N)))
  d0 <- reliability(items, type = "alpha", ci_method = "percentile", B = 30,
                    conf_level = 0.8, seed = 113)
  r <- recorder()
  d1 <- reliability(items, type = "alpha", ci_method = "percentile", B = 30,
                    conf_level = 0.8, seed = 113, progress = r$callback)
  expect_identical(d1, d0)
  expect_identical(r$done, 1:30)
  r <- recorder(stop_at = 4)
  expect_error(reliability(binary, type = "kr20", ci_method = "percentile",
                           B = 30, seed = 113, progress = r$callback),
               class = "dmar_bootstrap_cancelled")
  skip_if_not_installed("lavaan")
  skip_on_cran()
  o0 <- reliability_omega(items, ci_method = "percentile", B = 20, seed = 113,
                          conf_level = 0.8)
  r <- recorder()
  o1 <- reliability_omega(items, ci_method = "percentile", B = 20, seed = 113,
                          conf_level = 0.8, progress = r$callback)
  expect_identical(o1, o0)
  expect_identical(r$done, 1:20)
  m0 <- reliability_alpha(items, estimator = "model_implied",
                          ci_method = "percentile", B = 20, seed = 113,
                          conf_level = 0.8)
  r <- recorder()
  m1 <- reliability_alpha(items, estimator = "model_implied",
                          ci_method = "percentile", B = 20, seed = 113,
                          conf_level = 0.8, progress = r$callback)
  expect_identical(m1, m0)
  expect_identical(r$done, 1:20)
  h0 <- reliability_H(data = items, ci_method = "percentile", B = 20,
                      conf_level = 0.8, seed = 113)
  r <- recorder()
  h1 <- reliability_H(data = items, ci_method = "percentile", B = 20,
                      conf_level = 0.8, seed = 113, progress = r$callback)
  expect_identical(h1, h0)
  expect_identical(r$done, 1:20)
  fit <- lavaan::cfa("f =~ y1 + y2 + y3 + y4 + y5 + y6",
                     data = as.data.frame(items), std.lv = TRUE)
  f0 <- reliability_H(fit = fit, ci_method = "percentile", B = 20, seed = 113,
                      conf_level = 0.8)
  r <- recorder()
  f1 <- reliability_H(fit = fit, ci_method = "percentile", B = 20, seed = 113,
                      conf_level = 0.8, progress = r$callback)
  expect_identical(f1, f0)
  expect_identical(r$done, 1:20)
  r <- recorder(stop_at = 5)
  err <- expect_error(
    reliability_omega(items, ci_method = "bca", B = 20, seed = 113,
                      progress = r$callback),
    class = "dmar_bootstrap_cancelled")
  expect_identical(c(err$done, err$total), c(5L, 20L))
})

test_that("categorical omega reports the progress of both of its bootstraps", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  set.seed(113)
  n <- 300
  e <- rnorm(n)
  ordinal <- as.data.frame(vapply(1:4, function(j) {
    as.integer(cut(0.7 * e + sqrt(1 - 0.49) * rnorm(n),
                   breaks = c(-Inf, -1, 0, 1, Inf), labels = FALSE))
  }, integer(n)))
  names(ordinal) <- paste0("u", 1:4)
  c0 <- reliability_omega_categorical(ordinal, ci_method = "percentile",
                                      B = 10, seed = 113, conf_level = 0.5)
  r <- recorder()
  c1 <- reliability_omega_categorical(ordinal, ci_method = "percentile",
                                      B = 10, seed = 113, conf_level = 0.5,
                                      progress = r$callback)
  expect_identical(c1, c0)
  expect_identical(r$done, 1:10)
  # The bootstrap then impute replications under missing = "mi".
  skip_if_not_installed("mice")
  testthat::local_mocked_bindings(
    .mi_imputations_floor = function(imputations) as.integer(imputations),
    .mi_imputations_rule = function(fmi_upper, cv = 0.05) NA_integer_,
    .package = "DMAR")
  incomplete <- ordinal
  incomplete$u4[seq(1, n, by = 5)] <- NA
  run <- function(...) {
    suppressWarnings(reliability_omega_categorical(
      incomplete, ci_method = "percentile", missing = "mi", B = 6,
      imputations = 2, seed = 113, ...))
  }
  i0 <- run()
  r <- recorder()
  i1 <- run(progress = r$callback)
  expect_identical(i1, i0)
  expect_identical(r$done, 1:6)
  expect_identical(r$total, rep(6L, 6))
  r <- recorder(stop_at = 2)
  err <- expect_error(run(progress = r$callback),
                      class = "dmar_bootstrap_cancelled")
  expect_identical(c(err$done, err$total), c(2L, 6L))
})

test_that("ci_bifactor_indices() takes progress, and a seed gives the same intervals with or without it", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  set.seed(113)
  n <- 300
  g <- rnorm(n)
  grp <- list(rnorm(n), rnorm(n), rnorm(n))
  X <- vapply(1:9, function(i) {
    0.6 * g + 0.55 * grp[[ceiling(i / 3)]] +
      sqrt(1 - 0.6^2 - 0.55^2) * rnorm(n)
  }, numeric(n))
  colnames(X) <- paste0("x", 1:9)
  fit <- lavaan::cfa("g  =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
                      f1 =~ x1 + x2 + x3
                      f2 =~ x4 + x5 + x6
                      f3 =~ x7 + x8 + x9",
                     data = as.data.frame(X), orthogonal = TRUE,
                     std.lv = TRUE)
  b0 <- suppressWarnings(ci_bifactor_indices(fit, B = 12, seed = 113))
  r <- recorder()
  b1 <- suppressWarnings(ci_bifactor_indices(fit, B = 12, seed = 113,
                                             progress = r$callback))
  expect_identical(b1, b0)
  expect_identical(r$done, 1:12)
  r <- recorder(stop_at = 3)
  err <- expect_error(ci_bifactor_indices(fit, B = 12, seed = 113,
                                          progress = r$callback),
                      class = "dmar_bootstrap_cancelled")
  expect_identical(c(err$done, err$total), c(3L, 12L))
  expect_error(ci_bifactor_indices(fit, B = 12, progress = "yes"),
               "'progress' must be NULL or a function", fixed = TRUE)
})
