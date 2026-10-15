# Tests for the counted bootstrap behind the reliability family's
# intervals (.bootstrap_ci in R/reliability_internals.R): every
# replication is counted, the interval is formed from the replications
# that returned a value, and the counts travel with the result. The
# point function is coefficient alpha from the covariance matrix, which
# is cheap, so the mechanics can be exercised with hundreds of
# replications; the constructs that make a replication fail are built
# by hand so the counts have known values.

skip_if_not_installed("boot")

set.seed(113)
N <- 60
J <- 6
lam <- 0.7
eta <- rnorm(N)
items <- outer(eta, rep(lam, J)) +
  matrix(rnorm(N * J), N, J) * sqrt(1 - lam^2)
colnames(items) <- paste0("y", seq_len(J))

alpha_fn <- function(d) DMAR:::.alpha_from_S(stats::cov(d))

# The same resampling boot::boot performs inside .bootstrap_ci, so the
# tests can recompute what the helper should have returned.
boot_alpha <- function(fn, data, B, seed) {
  withr::with_seed(seed, boot::boot(data, function(d, i) fn(d[i, , drop = FALSE]),
                                    R = B, stype = "i"))
}


test_that("every replication is counted, with the reason it was set aside", {
  # Row 1 of `marked` carries the only 1 in its first column and row 2 the
  # only 2: a resample without row 1 returns a non-finite value, a resample
  # without row 2 raises an error, and the two reasons are counted apart.
  marked <- items
  marked[, 1] <- 0
  marked[1, 1] <- 1
  marked[2, 1] <- 2
  fn <- function(d) {
    if (!any(d[, 1] == 2)) stop("row 2 is missing from this resample")
    if (!any(d[, 1] == 1)) return(NA_real_)
    alpha_fn(d)
  }
  expect_warning(
    r <- DMAR:::.bootstrap_ci(marked, fn, B = 100, conf_level = 0.95,
                              kind = "percentile", seed = 113),
    regexp = "of 100 bootstrap replications were set aside"
  )
  bs <- r$bootstrap
  expect_identical(bs$attempted, 100L)
  expect_identical(bs$returned + bs$failed, 100L)
  expect_gt(bs$failed, 0L)
  expect_identical(sum(bs$reasons), bs$failed)
  expect_named(bs$reasons, c("error", "non_finite"))
  expect_gt(bs$reasons[["error"]], 0L)
  expect_gt(bs$reasons[["non_finite"]], 0L)
  expect_identical(bs$error_messages,
                   c("row 2 is missing from this resample" = bs$reasons[["error"]]))

  # The counts are what boot saw: the same resampling under the same seed
  # returns a finite value exactly `returned` times, and the percentile
  # limits are the quantiles of those returned replicates, as boot.ci
  # forms them after dropping the non-finite ones.
  bo <- boot_alpha(function(d) tryCatch(fn(d), error = function(e) NA_real_),
                   marked, B = 100, seed = 113)
  expect_identical(sum(is.finite(bo$t)), bs$returned)
  perc <- boot::boot.ci(bo, conf = 0.95, type = "perc")$percent
  expect_equal(c(r$lower, r$upper), perc[4:5], tolerance = 1e-12)
  expect_equal(r$se, stats::sd(bo$t, na.rm = TRUE), tolerance = 1e-12)
})

test_that("a clean bootstrap counts every replication as returned and warns about nothing", {
  expect_no_warning(
    r <- DMAR:::.bootstrap_ci(items, alpha_fn, B = 200, conf_level = 0.95,
                              kind = "percentile", seed = 113)
  )
  expect_identical(r$bootstrap$attempted, 200L)
  expect_identical(r$bootstrap$returned, 200L)
  expect_identical(r$bootstrap$failed, 0L)
  expect_identical(r$bootstrap$reasons, c(error = 0L, non_finite = 0L))
  expect_length(r$bootstrap$error_messages, 0L)
  expect_identical(r$ci_method, "percentile")
})

test_that("when fewer than two replications return a value there is no interval, and a warning says so", {
  expect_warning(
    r <- DMAR:::.bootstrap_ci(items, function(d) NA_real_, B = 20,
                              conf_level = 0.95, kind = "percentile", seed = 113),
    regexp = "fewer than two returned a value, so no interval is reported"
  )
  expect_identical(r$bootstrap$returned, 0L)
  expect_true(is.na(r$se))
  expect_true(is.na(r$lower))
  expect_true(is.na(r$upper))
})

test_that("the exported functions carry the counts, print them in one line, and add nothing elsewhere", {
  r <- reliability_alpha(data = items, ci_method = "percentile", B = 200,
                         seed = 113)
  bs <- attr(r, "bootstrap")
  expect_type(bs, "list")
  expect_identical(bs$attempted, 200L)
  expect_identical(bs$returned, 200L)
  expect_identical(attr(r, "B"), 200)
  expect_output(print(r),
                "Bootstrap: 200 replications attempted, 200 returned a value")

  # The same attribute from every function that bootstraps.
  binary <- (items > 0) * 1
  r_kr <- reliability_kr20(data = binary, ci_method = "percentile", B = 100,
                           seed = 113)
  expect_identical(attr(r_kr, "bootstrap")$attempted, 100L)
  skip_if_not_installed("lavaan")
  r_om <- reliability_omega(data = items, denominator = "observed",
                            ci_method = "percentile", B = 40, seed = 113)
  expect_identical(attr(r_om, "bootstrap")$attempted, 40L)
  expect_identical(attr(r_om, "bootstrap")$returned +
                     attr(r_om, "bootstrap")$failed, 40L)

  # No existing return field changes: a closed-form interval carries no
  # bootstrap attribute, the table has the same rows, and the print has
  # no counts line.
  r_cf <- reliability_alpha(data = items, ci_method = "ml")
  expect_null(attr(r_cf, "bootstrap"))
  expect_identical(r_cf$term, r$term)
  expect_false(any(grepl("Bootstrap:", capture.output(print(r_cf)))))
  expect_setequal(names(attributes(r)),
                  c(names(attributes(r_cf)), "B", "bootstrap"))
})


# ---- The BCa acceleration ---------------------------------------------------
#
# boot::boot.ci chooses between two estimates of the empirical influence
# values (boot::empinf): a regression on the resampling frequencies when
# at least n replications are available, the jackknife otherwise. The
# helper keeps that rule and computes the values itself, so on a clean
# construct it must reproduce boot's acceleration and limits exactly.

test_that("with at least n replications the BCa acceleration and limits equal boot's (regression)", {
  expect_no_warning(
    r <- DMAR:::.bootstrap_ci(items, alpha_fn, B = 200, conf_level = 0.95,
                              kind = "bca", seed = 113)
  )
  bo <- boot_alpha(alpha_fn, items, B = 200, seed = 113)
  bca <- boot::boot.ci(bo, conf = 0.95, type = "bca")$bca
  expect_equal(c(r$lower, r$upper), bca[4:5], tolerance = 1e-10)
  L <- boot::empinf(bo, type = "reg")
  expect_equal(r$bootstrap$acceleration, sum(L^3) / (6 * sum(L^2)^1.5),
               tolerance = 1e-10)
  expect_identical(r$bootstrap$acceleration_method, "regression")
  expect_identical(r$bootstrap$jackknife$attempted, 0L)
  expect_identical(r$ci_method, "bca")
})

test_that("with fewer than n replications the BCa acceleration and limits equal boot's (jackknife)", {
  expect_no_warning(
    r <- DMAR:::.bootstrap_ci(items, alpha_fn, B = 50, conf_level = 0.80,
                              kind = "bca", seed = 113)
  )
  bo <- boot_alpha(alpha_fn, items, B = 50, seed = 113)
  bca <- boot::boot.ci(bo, conf = 0.80, type = "bca")$bca
  expect_equal(c(r$lower, r$upper), bca[4:5], tolerance = 1e-10)
  L <- boot::empinf(bo, type = "jack")
  expect_equal(r$bootstrap$acceleration, sum(L^3) / (6 * sum(L^2)^1.5),
               tolerance = 1e-10)
  expect_identical(r$bootstrap$acceleration_method, "jackknife")
  expect_identical(r$bootstrap$jackknife$attempted, as.integer(N))
  expect_identical(r$bootstrap$jackknife$returned, as.integer(N))
  expect_identical(r$bootstrap$jackknife$failed, 0L)
})

test_that("one failed leave-one-out refit no longer stops the BCa interval (the 2026-10-14 report)", {
  skip_if_not_installed("lavaan")
  # The first item varies in one case only. Without that case the item is
  # constant and lavaan refuses the fit, so exactly one jackknife refit
  # returns no value; through 1.0.0 that single NA made boot's whole
  # influence vector NA and the call stopped with "estimated adjustment
  # 'a' is NA" (the Psychological Methods task hit it on every BFI
  # construct with B = 1,000 below n). B = 50 stays below n = 60, which
  # forces the jackknife regime; the 80 percent level keeps boot's
  # extreme-order-statistic warning out of a small run.
  items_bad <- items
  items_bad[, 1] <- 0
  items_bad[1, 1] <- 1
  # lavaan prints its variable table when it refuses the constant item;
  # the capture keeps that print out of the test output.
  w <- testthat::capture_warnings(invisible(utils::capture.output(
    r <- reliability_omega(data = items_bad, denominator = "observed",
                           ci_method = "bca", B = 50, conf_level = 0.80,
                           seed = 113)
  )))
  expect_length(w, 1L)
  expect_match(w, "of 50 bootstrap replications were set aside")
  expect_match(w, "1 of 60 jackknife replicates behind the BCa acceleration failed")
  expect_match(w, "estimated from the 59 that returned a value")
  expect_identical(attr(r, "ci_method"), "bca")
  bs <- attr(r, "bootstrap")
  expect_identical(bs$jackknife$attempted, 60L)
  expect_identical(bs$jackknife$returned, 59L)
  expect_identical(bs$jackknife$failed, 1L)
  expect_identical(bs$acceleration_method, "jackknife")
  expect_identical(bs$returned + bs$failed, 50L)
  expect_gt(bs$failed, 0L)   # the resamples without case 1 fail the same way

  # The acceleration is the jackknife over the 59 refits that returned a
  # value, formed as boot forms it: (n - 1)(t0 - t_(i)).
  est <- r$value[r$term == "estimate"]
  t_jack <- vapply(seq_len(N), function(i) {
    v <- NA_real_
    invisible(utils::capture.output(
      v <- tryCatch(DMAR:::.omega_observed_from_data(items_bad[-i, , drop = FALSE]),
                    error = function(e) NA_real_)))
    v
  }, numeric(1))
  expect_identical(which(!is.finite(t_jack)), 1L)
  L <- (N - 1) * (est - t_jack[is.finite(t_jack)])
  expect_equal(bs$acceleration, sum(L^3) / (6 * sum(L^2)^1.5), tolerance = 1e-10)

  ll <- r$value[r$term == "lower_limit"]
  ul <- r$value[r$term == "upper_limit"]
  expect_lte(ll, est)
  expect_gte(ul, est)
  expect_output(print(r),
                "jackknife for the BCa acceleration: 60 attempted, 59 returned a value \\(1 set aside: 1 returned a non-finite value\\)")
})

test_that("when the acceleration cannot be estimated the percentile interval is reported, labeled, and warned about", {
  # A point function that fails on every leave-one-out set (n - 1 rows)
  # and on no resample (n rows): the bootstrap is clean, the jackknife
  # returns nothing, and no acceleration exists.
  fn <- function(d) {
    if (nrow(d) == N - 1L) stop("leave-one-out refit failed")
    alpha_fn(d)
  }
  expect_warning(
    r <- DMAR:::.bootstrap_ci(items, fn, B = 50, conf_level = 0.80,
                              kind = "bca", seed = 113),
    regexp = "the percentile interval is reported instead, and the ci_method attribute says so"
  )
  expect_identical(
    r$ci_method,
    "percentile (BCa acceleration could not be estimated: 60 jackknife replicates failed)"
  )
  expect_identical(r$bootstrap$jackknife$attempted, 60L)
  expect_identical(r$bootstrap$jackknife$returned, 0L)
  expect_identical(r$bootstrap$jackknife$failed, 60L)
  expect_identical(r$bootstrap$jackknife$reasons, c(error = 60L, non_finite = 0L))
  expect_identical(r$bootstrap$jackknife$error_messages,
                   c("leave-one-out refit failed" = 60L))
  expect_true(is.na(r$bootstrap$acceleration))
  expect_true(is.na(r$bootstrap$acceleration_method))
  expect_identical(r$bootstrap$failed, 0L)
  # The limits are the percentile limits of the same bootstrap.
  p <- DMAR:::.bootstrap_ci(items, alpha_fn, B = 50, conf_level = 0.80,
                            kind = "percentile", seed = 113)
  expect_equal(c(r$lower, r$upper), c(p$lower, p$upper), tolerance = 1e-12)
})

test_that("when the bias correction cannot be estimated the percentile interval is reported and labeled", {
  # Equal to 0.8 on the full data and above it on every resample, so no
  # replicate falls below the estimate and boot's bias correction is
  # infinite; the acceleration itself is estimable.
  S0 <- sum(items[, 1])
  fn <- function(d) 0.8 + (sum(d[, 1]) - S0)^2 / 1e4
  expect_warning(
    r <- DMAR:::.bootstrap_ci(items, fn, B = 200, conf_level = 0.95,
                              kind = "bca", seed = 113),
    regexp = "percentile interval is reported instead"
  )
  expect_match(r$ci_method,
               "^percentile \\(BCa bias correction could not be estimated: every returned replicate fell on one side of the estimate\\)$")
  expect_true(is.na(r$bootstrap$acceleration))
})

test_that("the fallback label and the jackknife counts reach the exported result, its print, and glance()", {
  fallback <- "percentile (BCa acceleration could not be estimated: 3 jackknife replicates failed)"
  testthat::local_mocked_bindings(
    .bootstrap_ci = function(data, point_fn, B, conf_level, kind, seed) {
      list(se = 0.01, lower = 0.5, upper = 0.9, ci_method = fallback,
           bootstrap = list(
             attempted = 10L, returned = 10L, failed = 0L,
             reasons = c(error = 0L, non_finite = 0L),
             error_messages = integer(),
             acceleration = NA_real_, acceleration_method = NA_character_,
             jackknife = list(attempted = 60L, returned = 57L, failed = 3L,
                              reasons = c(error = 3L, non_finite = 0L),
                              error_messages = c("refit failed" = 3L))))
    },
    .package = "DMAR"
  )
  r <- reliability_alpha(data = items, ci_method = "bca", B = 10)
  expect_identical(attr(r, "ci_method"), fallback)
  expect_identical(generics::glance(r)$ci_method, fallback)
  expect_identical(attr(r, "bootstrap")$jackknife$failed, 3L)
  expect_output(
    print(r),
    "Bootstrap: 10 replications attempted, 10 returned a value; jackknife for the BCa acceleration: 60 attempted, 57 returned a value \\(3 set aside: 3 raised an error \\(most often: refit failed\\)\\)",
    fixed = FALSE
  )
})
