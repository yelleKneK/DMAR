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
