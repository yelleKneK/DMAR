# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# Tests for the counted bootstrap behind the reliability family's
# intervals (.bootstrap_ci in R/reliability_internals.R): every
# replication is counted, the interval is formed from the replications
# that returned a value, and the counts travel with the result. The
# point function is coefficient alpha from the covariance matrix, which
# is cheap, so the mechanics can be exercised with hundreds of
# replications; the constructs that make a replication fail are built
# by hand so the counts have known values.

# The engine as every caller in the package runs it: no interval from
# fewer than half of the replications. (Through the first versions of
# 1.1.0 these tests ran it with a switch, require_half = FALSE, that kept
# the floor of two returned values of 1.0.1 and that no caller used; the
# switch is gone, and a test that needs an interval lets at least half of
# its replications return a value.)
bootstrap_for_mechanics <- function(...) DMAR:::.bootstrap_ci(...)

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
  # Rows 1 and 2 of `marked` carry the only 1s in its first column and
  # rows 3 and 4 the only 2s: a resample without a 1 returns a non-finite
  # value, a resample without a 2 raises an error, and the two reasons are
  # counted apart. Two rows each, so that about three resamples in four
  # return a value and the interval is formed (no interval is formed from
  # fewer than half of the replications).
  marked <- items
  marked[, 1] <- 0
  marked[1:2, 1] <- 1
  marked[3:4, 1] <- 2
  fn <- function(d) {
    if (!any(d[, 1] == 2)) stop("row 2 is missing from this resample")
    if (!any(d[, 1] == 1)) return(NA_real_)
    alpha_fn(d)
  }
  expect_warning(
    r <- bootstrap_for_mechanics(marked, fn, B = 100, conf_level = 0.95,
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

test_that("with fewer than two replications returned there is no caution about limits, and the printed line says no interval is reported", {
  # Regression test. The caution about the limits ran whenever the
  # replications set aside reached the tail area, and so followed "no
  # interval is reported" in the warning, and the printed line gave it
  # without saying that no interval was formed; it named an improper
  # solution for coefficients that fit no model; and the quoted refusal
  # kept its own period ("refused.))").
  for (kind in c("percentile", "bootstrap_se")) {
    w <- testthat::capture_warnings(
      r <- bootstrap_for_mechanics(items, function(d) NA_real_, B = 20,
                                conf_level = 0.95, kind = kind,
                                seed = 113))
    expect_match(w, "fewer than two returned a value, so no interval is reported",
                 fixed = TRUE, all = FALSE)
    expect_false(any(grepl("Caution", w, fixed = TRUE)))
    expect_true(r$bootstrap$caution)
    line <- DMAR:::.bootstrap_footer_line(
      structure(data.frame(term = "conf_level", value = 0.95),
                bootstrap = r$bootstrap))
    expect_match(line, "; fewer than two returned a value, so no interval is reported$")
    expect_no_match(line, "caution", fixed = TRUE)
  }
  # One replication returned (the first call is the estimate itself).
  calls <- new.env()
  calls$n <- 0L
  one <- function(d) {
    calls$n <- calls$n + 1L
    if (calls$n <= 2L) alpha_fn(d) else stop("refused.")
  }
  w <- testthat::capture_warnings(
    r <- bootstrap_for_mechanics(items, one, B = 20, conf_level = 0.95,
                              kind = "percentile", seed = 113))
  expect_identical(r$bootstrap$returned, 1L)
  expect_false(any(grepl("Caution", w, fixed = TRUE)))
  expect_match(w, "(most often: refused)", fixed = TRUE, all = FALSE)
  # With an interval formed and nothing improper kept, the caution names
  # only the replications set aside.
  # (About three resamples in four return a value: no interval is formed
  # from fewer than half.)
  half <- function(d) {
    if (mean(d[, 1]) > mean(items[, 1]) - 0.1) alpha_fn(d) else NA_real_
  }
  w <- testthat::capture_warnings(
    r <- bootstrap_for_mechanics(items, half, B = 100, conf_level = 0.95,
                              kind = "percentile", seed = 113))
  expect_true(r$bootstrap$caution && r$bootstrap$returned >= 2L)
  expect_match(w, "replications \\([0-9.]+%\\) were set aside, at least",
               all = FALSE)
  expect_false(any(grepl("improper solution", w, fixed = TRUE)))
  line <- DMAR:::.bootstrap_footer_line(
    structure(data.frame(term = "conf_level", value = 0.95),
              bootstrap = r$bootstrap))
  expect_match(line, "caution: [0-9.]+% of the replications were set aside, at least")
  # The bootstrap then impute warning of categorical omega follows the
  # same rule.
  counts <- list(failed = 19L, returned = 1L,
                 reasons = c(error = 0L, non_finite = 19L, imputation = 0L),
                 error_messages = integer(), undetermined = 0L,
                 improper = 0L, caution = TRUE, require_half = TRUE,
                 imputation_messages = integer())
  notes <- DMAR:::.boot_mi_notes(counts, B = 20L, conf_level = 0.95,
                                 imputations = 2L, fmi = NA_real_,
                                 widening = NA_real_)
  expect_false(any(grepl("Caution", notes, fixed = TRUE)))
})


test_that("a clean bootstrap counts every replication as returned and warns about nothing", {
  expect_no_warning(
    r <- bootstrap_for_mechanics(items, alpha_fn, B = 200, conf_level = 0.95,
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
    r <- bootstrap_for_mechanics(items, function(d) NA_real_, B = 20,
                              conf_level = 0.95, kind = "percentile", seed = 113),
    regexp = "fewer than two returned a value, so no interval is reported"
  )
  expect_identical(r$bootstrap$returned, 0L)
  expect_true(is.na(r$se))
  expect_true(is.na(r$lower))
  expect_true(is.na(r$upper))
  # With no replication set aside but only one attempted, the warning says
  # so too; through 1.0.1 it was silent there.
  expect_warning(
    r <- bootstrap_for_mechanics(items, alpha_fn, B = 1, conf_level = 0.95,
                              kind = "percentile", seed = 113),
    "Only 1 bootstrap replication was attempted", fixed = TRUE)
  expect_true(is.na(r$lower))
})

test_that("finite replicates that are all equal give the common value at both limits", {
  # With no replication set aside, boot.ci() would print "All values of t
  # are equal" and return NULL; with some set aside (NA), its test for
  # that case meets if (NA) and stops. Either way the limits are the
  # common value, for every kind, with a warning.
  # About half of the resamples hold more of the positive y1 scores than
  # the sample does, and those return no value.
  thr <- sum(items[, 1] > 0)
  half <- function(d) if (sum(d[, 1] > 0) > thr) NA_real_ else 0.5
  for (fn in list(function(d) 0.5, half)) {
    for (k in c("percentile", "bca", "bootstrap_se", "bootstrap_se_logistic")) {
      w <- testthat::capture_warnings(
        out <- utils::capture.output(
          r <- bootstrap_for_mechanics(items, fn, B = 50, conf_level = 0.95,
                                    kind = k, seed = 113)))
      expect_identical(c(r$lower, r$upper), c(0.5, 0.5))
      expect_match(w, "gave the same estimate", fixed = TRUE, all = FALSE)
      expect_false(any(grepl("All values of t", out, fixed = TRUE)))
    }
  }
  expect_true(DMAR:::.boot_constant(c(0.5, NA, 0.5 + 1e-12)))
  expect_false(DMAR:::.boot_constant(c(0.5, NA, 0.51)))
  expect_false(DMAR:::.boot_constant(c(NA_real_, NA_real_)))
})

test_that("the exported functions carry the counts, print them in one line, and add nothing elsewhere", {
  r <- reliability_alpha(data = items, ci_method = "percentile", B = 200,
                         seed = 113)
  bs <- attr(r, "bootstrap")
  expect_type(bs, "list")
  expect_identical(bs$attempted, 200L)
  expect_identical(bs$returned, 200L)
  expect_identical(attr(r, "B"), 200L)
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
  # no counts line. A bootstrap adds `B`, `bootstrap`, and (since 1.1.0)
  # `B_used`, the effective count every bootstrap in the package reports.
  r_cf <- reliability_alpha(data = items, ci_method = "ml")
  expect_null(attr(r_cf, "bootstrap"))
  expect_null(attr(r_cf, "B_used"))
  expect_identical(r_cf$term, r$term)
  expect_false(any(grepl("Bootstrap:", capture.output(print(r_cf)))))
  expect_setequal(names(attributes(r)),
                  c(names(attributes(r_cf)), "B", "bootstrap", "B_used",
                    "n_dropped", "B_required", "bootstrap_status"))
  expect_identical(attr(r, "B_used"), bs$returned)
  expect_identical(attr(r, "n_dropped"), bs$attempted - bs$returned)
  expect_identical(attr(r, "B_required"), bs$required)
  expect_identical(attr(r, "bootstrap_status"), bs$bootstrap_status)
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
    r <- bootstrap_for_mechanics(items, alpha_fn, B = 200, conf_level = 0.95,
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
    r <- bootstrap_for_mechanics(items, alpha_fn, B = 50, conf_level = 0.80,
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
  # constant and is refused before the fit (lavaan refused it until the
  # check reached the replications), so exactly one jackknife refit
  # returns no value; through 1.0.0 that single NA made boot's whole
  # influence vector NA and the call stopped with "estimated adjustment
  # 'a' is NA" (the Psychological Methods task hit it on every BFI
  # construct with B = 1,000 below n). B = 50 stays below n = 60, which
  # forces the jackknife regime; the 80 percent level keeps boot's
  # extreme-order-statistic warning out of a small run.
  items_bad <- items
  items_bad[, 1] <- 0
  items_bad[1, 1] <- 1
  # The jackknife values below are computed with the fit itself, which
  # reaches lavaan, and lavaan prints its variable table when it refuses
  # the constant item; the captures keep any such print out of the test
  # output.
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

  ll <- r$value[r$term == "ci_lower"]
  ul <- r$value[r$term == "ci_upper"]
  expect_lte(ll, est)
  expect_gte(ul, est)
  expect_output(print(r),
                "jackknife for the BCa acceleration: 60 attempted, 59 returned a value \\(1 set aside: 1 raised an error \\(most often: [^ ]+ takes a single value in this resample")
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
    r <- bootstrap_for_mechanics(items, fn, B = 50, conf_level = 0.80,
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
  p <- bootstrap_for_mechanics(items, alpha_fn, B = 50, conf_level = 0.80,
                            kind = "percentile", seed = 113)
  expect_equal(c(r$lower, r$upper), c(p$lower, p$upper), tolerance = 1e-12)

  # The floor is two returned values: with exactly one (only the refit
  # without row 1 succeeds; row 1 carries the only 1 in its column) boot's
  # centering on t0 would give the finite extreme of 1/6 in magnitude,
  # and the percentile interval is reported instead.
  marked <- items
  marked[, 1] <- 0
  marked[1, 1] <- 1
  fn_one <- function(d) {
    if (nrow(d) == N - 1L && any(d[, 1] == 1)) stop("leave-one-out refit failed")
    alpha_fn(d)
  }
  expect_warning(
    r1 <- bootstrap_for_mechanics(marked, fn_one, B = 50, conf_level = 0.80,
                               kind = "bca", seed = 113),
    regexp = "the percentile interval is reported instead"
  )
  expect_identical(r1$bootstrap$jackknife$returned, 1L)
  expect_identical(r1$bootstrap$jackknife$failed, 59L)
  expect_identical(
    r1$ci_method,
    "percentile (BCa acceleration could not be estimated: 59 jackknife replicates failed)"
  )
  expect_true(is.na(r1$bootstrap$acceleration))
})

test_that("when the bias correction cannot be estimated the percentile interval is reported and labeled", {
  # Equal to 0.8 on the full data and above it on every resample, so no
  # replicate falls below the estimate and boot's bias correction is
  # infinite; the acceleration itself is estimable.
  S0 <- sum(items[, 1])
  fn <- function(d) 0.8 + (sum(d[, 1]) - S0)^2 / 1e4
  expect_warning(
    r <- bootstrap_for_mechanics(items, fn, B = 200, conf_level = 0.95,
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
    .bootstrap_ci = function(data, point_fn, B, conf_level, kind, seed,
                             ncpus = 1L, improper = "keep",
                             two_rows = NULL, progress = NULL, fun = NULL) {
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

test_that("a replication's named cause is counted within the non-finite ones and named in the text", {
  # .omega_c_from_data() returns NA_real_ with a dmar_reason attribute
  # when a resample leaves a polychoric correlation undetermined or gives
  # an improper solution. The two causes are counted apart, as subsets
  # of the non-finite count, so `reasons` keeps its 1.0.1 shape.
  marked <- items
  marked[, 1] <- 0
  marked[1, 1] <- 1
  marked[2, 1] <- 2
  fn <- function(d) {
    if (!any(d[, 1] == 1)) {
      return(structure(NA_real_, dmar_reason = "undetermined_polychoric"))
    }
    if (!any(d[, 1] == 2)) return(structure(NA_real_, dmar_reason = "improper"))
    alpha_fn(d)
  }
  w <- testthat::capture_warnings(
    r <- bootstrap_for_mechanics(marked, fn, B = 100, conf_level = 0.95,
                              kind = "percentile", seed = 113))
  bs <- r$bootstrap
  expect_named(bs$reasons, c("error", "non_finite"))
  expect_identical(bs$reasons[["error"]], 0L)
  expect_gt(bs$undetermined, 0L)
  expect_gt(bs$improper, 0L)
  expect_identical(bs$undetermined + bs$improper, bs$reasons[["non_finite"]])
  expect_identical(bs$failed, bs$reasons[["non_finite"]])
  expect_length(w, 1L)
  expect_match(w, sprintf("(%d left a polychoric correlation undetermined, %d gave an improper (Heywood) solution)",
                          bs$undetermined, bs$improper), fixed = TRUE)
  expect_false(grepl("returned a non-finite value", w, fixed = TRUE))
  expect_match(w, 'denominator = "model_implied"', fixed = TRUE)
  expect_match(w, "trim the upper tail", fixed = TRUE)

  # A clean run counts neither, and every coefficient's bootstrap carries
  # the two counts.
  r0 <- bootstrap_for_mechanics(items, alpha_fn, B = 50, conf_level = 0.95,
                             kind = "percentile", seed = 113)
  expect_identical(r0$bootstrap$undetermined, 0L)
  expect_identical(r0$bootstrap$improper, 0L)
  expect_identical(r0$bootstrap$jackknife$undetermined, 0L)

  # A bootstrap attribute written before the counts existed still prints.
  expect_identical(DMAR:::.bootstrap_reason_phrase(c(error = 0L, non_finite = 2L)),
                   "2 returned a non-finite value")
  expect_identical(DMAR:::.bootstrap_reason_phrase(c(error = 0L, non_finite = 3L),
                                                   integer(), 2L, NULL),
                   "2 left a polychoric correlation undetermined, 1 returned a non-finite value")
  # A refit that failed or did not converge, where the bootstrap cannot
  # tell the two apart (lavaan's, behind average_variance_extracted()), is
  # named for both causes, in the line under the table as in the warning;
  # an attribute without the count prints as before.
  expect_identical(DMAR:::.bootstrap_reason_phrase(c(error = 0L, non_finite = 3L),
                                                   refit_failed = 2L),
                   "2 failed to refit or did not converge, 1 returned a non-finite value")
  x <- data.frame(term = c("estimate", "conf_level"), value = c(0.8, 0.95))
  attr(x, "bootstrap") <- list(
    attempted = 100L, returned = 97L, failed = 3L,
    reasons = c(error = 0L, non_finite = 3L), error_messages = integer(),
    refit_failed = 3L)
  expect_identical(DMAR:::.bootstrap_footer_line(x),
                   paste0("Bootstrap: 100 replications attempted, 97 ",
                          "returned a value (3 set aside: 3 failed to refit ",
                          "or did not converge)"))
  attr(x, "bootstrap")$refit_failed <- NULL
  expect_match(DMAR:::.bootstrap_footer_line(x),
               "(3 set aside: 3 returned a non-finite value)", fixed = TRUE)
})


# ---- One engine: each replication's fate, the improper rule, workers -------
#
# Since 1.1.0 the statistic boot runs returns, for every replication, the
# value the interval uses, a status code, and the estimate, so the counts
# are read from boot's own replicate matrix (and survive forked workers),
# and the draws and their fates are returned in the attribute.

# Row 1 carries the only 1 in the first column and row 2 the only 2, so a
# resample without row 1 or without row 2 is identifiable.
marked <- items
marked[, 1] <- 0
marked[1, 1] <- 1
marked[2, 1] <- 2
fn_mixed <- function(d) {
  if (!any(d[, 1] == 2)) stop("row 2 is missing from this resample")
  if (!any(d[, 1] == 1)) return(NA_real_)
  alpha_fn(d)
}

test_that("the interval is boot's own on the first column, bit for bit", {
  # The engine's interval for a clean run equals boot.ci() on the
  # one-column statistic the engine used through 1.0.1.
  bo <- boot_alpha(alpha_fn, items, B = 200, seed = 113)
  p <- bootstrap_for_mechanics(items, alpha_fn, B = 200, conf_level = 0.95,
                            kind = "percentile", seed = 113)
  expect_identical(c(p$lower, p$upper),
                   boot::boot.ci(bo, conf = 0.95, type = "perc")$percent[4:5])
  b <- bootstrap_for_mechanics(items, alpha_fn, B = 200, conf_level = 0.95,
                            kind = "bca", seed = 113)
  expect_identical(c(b$lower, b$upper),
                   boot::boot.ci(bo, conf = 0.95, type = "bca")$bca[4:5])
  expect_identical(p$se, stats::sd(as.numeric(bo$t), na.rm = TRUE))
})

test_that("every replication's estimate and fate travel in the attribute", {
  r <- suppressWarnings(
    bootstrap_for_mechanics(marked, fn_mixed, B = 100, conf_level = 0.95,
                         kind = "percentile", seed = 113))
  bs <- r$bootstrap
  expect_length(bs$t, 100L)
  expect_s3_class(bs$status, "factor")
  expect_identical(levels(bs$status),
                   c("proper", "improper", "error", "non_finite",
                     "undetermined", "not_converged"))
  expect_identical(sum(bs$status == "error"), bs$reasons[["error"]])
  expect_identical(sum(bs$status == "non_finite"), bs$reasons[["non_finite"]])
  expect_identical(sum(bs$status == "proper"), bs$returned)
  expect_true(all(is.na(bs$t[bs$status != "proper"])))
  # t holds the estimate boot saw, replication by replication.
  bo <- boot_alpha(function(d) tryCatch(fn_mixed(d), error = function(e) NA_real_),
                   marked, B = 100, seed = 113)
  expect_identical(bs$t, as.numeric(bo$t))
  expect_identical(bs$improper_rule, "keep")
  expect_identical(bs$improper_kept, 0L)
  expect_identical(bs$ncpus, 1L)
  # About 60 percent of the resamples omit row 1 or row 2, far beyond the
  # 2.5 percent beyond each limit, so the caution is set.
  expect_true(bs$caution)
})

test_that("improper replications are kept or set aside by the rule and counted either way", {
  # The replications whose resample omits row 1 are marked improper; the
  # estimate is alpha itself, so "keep" reproduces the clean run.
  fn_improper <- function(d) {
    v <- alpha_fn(d)
    if (!any(d[, 1] == 1)) attr(v, "dmar_reason") <- "improper"
    v
  }
  one <- marked
  one[2, 1] <- 0
  w_keep <- testthat::capture_warnings(
    k <- bootstrap_for_mechanics(one, fn_improper, B = 200, conf_level = 0.95,
                              kind = "percentile", seed = 113))
  w_drop <- testthat::capture_warnings(
    d <- bootstrap_for_mechanics(one, fn_improper, B = 200, conf_level = 0.95,
                              kind = "percentile", seed = 113,
                              improper = "drop"))
  bk <- k$bootstrap
  bd <- d$bootstrap
  n_imp <- sum(bk$status == "improper")
  expect_gt(n_imp, 0L)

  expect_identical(bk$improper_kept, n_imp)
  expect_identical(bk$improper, 0L)
  expect_identical(bk$returned, 200L)
  expect_identical(bk$failed, 0L)
  clean <- bootstrap_for_mechanics(one, alpha_fn, B = 200, conf_level = 0.95,
                                kind = "percentile", seed = 113)
  expect_identical(c(k$lower, k$upper), c(clean$lower, clean$upper))

  expect_identical(bd$improper, n_imp)
  expect_identical(bd$improper_kept, 0L)
  expect_identical(bd$returned, 200L - n_imp)
  expect_identical(bd$reasons, c(error = 0L, non_finite = n_imp))
  expect_identical(bd$improper_rule, "drop")
  bo <- boot_alpha(function(x) if (!any(x[, 1] == 1)) NA_real_ else alpha_fn(x),
                   one, B = 200, seed = 113)
  expect_equal(c(d$lower, d$upper),
               boot::boot.ci(bo, conf = 0.95, type = "perc")$percent[4:5],
               tolerance = 1e-12)

  # The same draws and estimates under either rule: the estimate of a
  # replication set aside under "drop" stays in t, so the other interval
  # can be formed without refitting.
  expect_identical(bk$t, bd$t)
  expect_identical(bk$status, bd$status)
  expect_true(all(is.finite(bd$t)))

  # Both warn once, with the caution (about 37 percent are improper).
  expect_length(w_keep, 1L)
  expect_match(w_keep, sprintf("Caution: %d of 200 replications", n_imp))
  expect_length(w_drop, 1L)
  expect_match(w_drop, sprintf("%d of 200 bootstrap replications were set aside (%d gave an improper (Heywood) solution)",
                               n_imp, n_imp), fixed = TRUE)
})

test_that("an improper solution that yields no coefficient is set aside under either rule", {
  fn_none <- function(d) {
    if (!any(d[, 1] == 1)) return(structure(NA_real_, dmar_reason = "improper"))
    alpha_fn(d)
  }
  one <- marked
  one[2, 1] <- 0
  w <- testthat::capture_warnings(
    r <- bootstrap_for_mechanics(one, fn_none, B = 100, conf_level = 0.95,
                              kind = "percentile", seed = 113))
  bs <- r$bootstrap
  expect_gt(bs$improper, 0L)
  expect_identical(bs$improper_kept, 0L)
  expect_identical(bs$improper, bs$failed)
  expect_true(all(is.na(bs$t[bs$status == "improper"])))
  expect_match(w, "set aside only when its solution yields no coefficient",
               fixed = TRUE)
})

test_that("the caution is set when the replications in question reach the tail area, ties included", {
  # 25 of 1000 at 95 percent is exactly the 2.5 percent beyond a limit,
  # although (1 - 0.95) / 2 is slightly above 0.025 in floating point.
  expect_true(DMAR:::.boot_caution(25L, 1000L, 0.95))
  expect_false(DMAR:::.boot_caution(24L, 1000L, 0.95))
  expect_true(DMAR:::.boot_caution(1L, 40L, 0.95))
  expect_false(DMAR:::.boot_caution(0L, 20L, 0.95))
  expect_true(DMAR:::.boot_caution(50L, 1000L, 0.90))
  expect_false(DMAR:::.boot_caution(49L, 1000L, 0.90))
})

test_that("the printed line names the improper replications kept and the caution", {
  x <- data.frame(term = c("estimate", "conf_level"), value = c(0.8, 0.95))
  attr(x, "bootstrap") <- list(
    attempted = 200L, returned = 190L, failed = 10L,
    reasons = c(error = 0L, non_finite = 10L), error_messages = integer(),
    undetermined = 0L, improper = 4L, improper_kept = 12L,
    improper_rule = "keep", caution = TRUE,
    jackknife = list(attempted = 0L))
  expect_identical(
    DMAR:::.bootstrap_footer_line(x),
    paste0("Bootstrap: 200 replications attempted, 190 returned a value ",
           "(10 set aside: 4 gave an improper (Heywood) solution, 6 ",
           "returned a non-finite value); 12 of the 190 came from an ",
           "improper (Heywood) solution and were kept; caution: 11% of ",
           "the replications came from an improper solution or were set ",
           "aside, at least the 2.5% beyond each limit"))
  # Replications with no single-factor solution are named apart from the
  # Heywood cases, and the caution says which way the lower limit errs.
  attr(x, "bootstrap") <- list(
    attempted = 200L, returned = 170L, failed = 30L,
    reasons = c(error = 0L, non_finite = 30L), error_messages = integer(),
    undetermined = 0L, improper = 30L, improper_kept = 5L,
    no_solution = 26L, improper_rule = "keep", caution = TRUE,
    jackknife = list(attempted = 0L))
  expect_identical(
    DMAR:::.bootstrap_footer_line(x),
    paste0("Bootstrap: 200 replications attempted, 170 returned a value ",
           "(30 set aside: 4 gave an improper (Heywood) solution, 26 had ",
           "no single-factor solution); 5 of the 170 came from an improper ",
           "(Heywood) solution and were kept; caution: 17.5% of the ",
           "replications came from an improper solution or were set aside, ",
           "at least the 2.5% beyond each limit; the lower limit is likely ",
           "too high"))
  # An attribute written before 1.1.0 prints as it did.
  attr(x, "bootstrap") <- list(attempted = 10L, returned = 10L, failed = 0L,
                               reasons = c(error = 0L, non_finite = 0L),
                               error_messages = integer())
  expect_identical(DMAR:::.bootstrap_footer_line(x),
                   "Bootstrap: 10 replications attempted, 10 returned a value")
  attr(x, "bootstrap") <- NULL
  expect_null(DMAR:::.bootstrap_footer_line(x))
})

test_that("the interval, counts, draws, and messages for a seed do not depend on ncpus", {
  skip_on_os("windows")   # forked workers; the socket path is tested below
  # Errors, non-finite returns, and improper estimates together, under
  # both regimes of the BCa acceleration: regression (B >= n) and the
  # jackknife (B < n), which also runs on the workers.
  fn <- function(d) {
    if (!any(d[, 1] == 2)) stop("row 2 is missing from this resample")
    v <- alpha_fn(d)
    if (!any(d[, 1] == 1)) attr(v, "dmar_reason") <- "improper"
    v
  }
  for (B in c(120L, 50L)) {
    for (rule in c("keep", "drop")) {
      w1 <- testthat::capture_warnings(
        r1 <- bootstrap_for_mechanics(marked, fn, B = B, conf_level = 0.9,
                                   kind = "bca", seed = 113,
                                   improper = rule))
      w2 <- testthat::capture_warnings(
        r2 <- bootstrap_for_mechanics(marked, fn, B = B, conf_level = 0.9,
                                   kind = "bca", seed = 113, ncpus = 2L,
                                   improper = rule))
      expect_identical(r1[c("se", "lower", "upper", "ci_method")],
                       r2[c("se", "lower", "upper", "ci_method")])
      b1 <- r1$bootstrap
      b2 <- r2$bootstrap
      expect_identical(b1[setdiff(names(b1), "ncpus")],
                       b2[setdiff(names(b2), "ncpus")])
      expect_identical(b2$ncpus, 2L)
      expect_identical(w1, w2)
      expect_identical(b1$error_messages,
                       c("row 2 is missing from this resample" = b1$reasons[["error"]]))
      # With B below n and at least half of the replications returning a
      # value (the improper ones kept), the acceleration is the
      # jackknife's.
      if (B == 50L && rule == "keep") {
        expect_gte(b1$returned, 25L)
        expect_identical(b1$acceleration_method, "jackknife")
      }
    }
  }
})

test_that("a session's boot.parallel option no longer loses the counts", {
  skip_on_os("windows")
  # Through 1.0.1 boot::boot() was called without `parallel`, so this
  # option ran the replications in forked workers, whose counters the
  # parent never saw: every failure went uncounted.
  r0 <- suppressWarnings(
    bootstrap_for_mechanics(marked, fn_mixed, B = 100, conf_level = 0.95,
                         kind = "percentile", seed = 113))
  withr::local_options(boot.parallel = "multicore", boot.ncpus = 2L)
  r <- suppressWarnings(
    bootstrap_for_mechanics(marked, fn_mixed, B = 100, conf_level = 0.95,
                         kind = "percentile", seed = 113))
  expect_gt(r$bootstrap$failed, 0L)
  expect_identical(r$bootstrap, r0$bootstrap)
})

test_that("ncpus is validated", {
  for (bad in list(0, 1.5, NA, "2", c(1, 2))) {
    expect_error(DMAR:::.check_ncpus(bad), "'ncpus' must be a single whole number")
  }
  expect_identical(DMAR:::.check_ncpus(2), 2L)
  expect_error(reliability_alpha(items, ci_method = "percentile", B = 20,
                                 ncpus = 0),
               "'ncpus' must be a single whole number")
})

test_that("the socket cluster path, Windows's, returns the serial interval", {
  skip_on_cran()
  # The socket workers load the installed DMAR, so the test runs only
  # where the installed version carries this engine (not against a
  # source tree loaded over an older installation).
  cl <- parallel::makePSOCKcluster(1L)
  has_engine <- tryCatch(
    isTRUE(parallel::clusterEvalQ(cl, exists(".boot_classify",
                                             envir = asNamespace("DMAR")))[[1]]),
    error = function(e) FALSE)
  parallel::stopCluster(cl)
  skip_if_not(has_engine, "the installed DMAR predates this engine")
  testthat::local_mocked_bindings(
    .boot_parallel_type = function(ncpus) if (ncpus > 1L) "snow" else "no",
    .package = "DMAR")
  r1 <- bootstrap_for_mechanics(items, alpha_fn, B = 50, conf_level = 0.9,
                             kind = "bca", seed = 113)
  r2 <- bootstrap_for_mechanics(items, alpha_fn, B = 50, conf_level = 0.9,
                             kind = "bca", seed = 113, ncpus = 2L)
  expect_identical(r1[c("se", "lower", "upper")], r2[c("se", "lower", "upper")])
  expect_identical(r1$bootstrap$t, r2$bootstrap$t)
  expect_identical(r1$bootstrap$acceleration, r2$bootstrap$acceleration)
})

test_that("the exported functions take ncpus, and a seed gives one interval for any ncpus", {
  skip_on_os("windows")
  a1 <- reliability_alpha(items, ci_method = "percentile", B = 200, seed = 113)
  a2 <- reliability_alpha(items, ci_method = "percentile", B = 200, seed = 113,
                          ncpus = 2)
  expect_identical(a1$value, a2$value)
  expect_identical(attr(a2, "bootstrap")$ncpus, 2L)
  binary <- (items > 0) * 1
  k1 <- reliability_kr20(binary, ci_method = "bca", B = 50, seed = 113,
                         conf_level = 0.8)
  k2 <- reliability_kr20(binary, ci_method = "bca", B = 50, seed = 113,
                         conf_level = 0.8, ncpus = 2)
  expect_identical(k1$value, k2$value)
  skip_on_cran()
  skip_if_not_installed("lavaan")
  o1 <- reliability_omega(items, ci_method = "percentile", B = 40, seed = 113)
  o2 <- reliability_omega(items, ci_method = "percentile", B = 40, seed = 113,
                          ncpus = 2)
  expect_identical(o1$value, o2$value)
  expect_identical(attr(o1, "bootstrap")$t, attr(o2, "bootstrap")$t)
})

test_that("the functions sent to workers carry their arguments as values", {
  skip_on_os("windows")
  skip_on_cran()
  # A fork cluster shares this session's DMAR but, like the socket cluster
  # the engine uses on Windows, sends each function to its workers
  # serialized, so it exercises that path without an installed DMAR. An
  # argument left as a promise would arrive unevaluable: before this was
  # caught, every leave-one-out refit of the jackknife failed on a
  # socket cluster.
  one <- DMAR:::.boot_leave_one_out(alpha_fn, items, "keep")
  expect_setequal(ls(environment(one)), c("point_fn", "data", "improper", "m"))
  st <- DMAR:::.boot_statistic(alpha_fn, "keep")
  expect_setequal(ls(environment(st)), c("point_fn", "improper", "m"))
  cl <- parallel::makeForkCluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  j1 <- DMAR:::.boot_jackknife(alpha_fn, items, "keep", "no", 1L, NULL)
  j2 <- DMAR:::.boot_jackknife(alpha_fn, items, "keep", "snow", 2L, cl)
  expect_identical(j1, j2)
  expect_true(all(j2[, 2] == 1))
  b1 <- withr::with_seed(113, boot::boot(items, st, R = 30, stype = "i"))
  b2 <- withr::with_seed(113, boot::boot(items, st, R = 30, stype = "i",
                                         parallel = "snow", ncpus = 2L,
                                         cl = cl))
  expect_identical(b1$t, b2$t)
})

test_that("the replications that raised an error are refit on the engine's processes, with the serial messages and counts", {
  skip_on_os("windows")   # forked workers and a fork cluster
  skip_on_cran()
  # Regression test. The refits that recover the error messages ran in
  # this process alone whatever ncpus was, and a replication that raises
  # an error is most often one whose fit ran off (no maximum), the slowest
  # kind, so on several processes they took most of the run. They now run
  # on the engine's processes; the messages, their counts, and the count
  # of refusals for too few cases are those of the serial path.
  d <- data.frame(x = seq_len(30))
  fn <- function(dd) {
    r <- sum(dd$x) %% 10
    if (r < 2) {
      stop(structure(class = c("dmar_too_few_cases", "error", "condition"),
                     list(message = "too few distinct cases", call = NULL)))
    }
    if (r < 4) stop("first message")
    if (r < 5) stop("second message")
    mean(dd$x)
  }
  sets <- withr::with_seed(113, lapply(seq_len(40), function(k)
    sample.int(30, replace = TRUE)))
  one <- DMAR:::.boot_refit_one(fn, d)
  expect_setequal(ls(environment(one)), c("point_fn", "data"))
  serial <- DMAR:::.boot_refit_errors(fn, d, sets)
  expect_gt(serial$too_few_cases, 0L)
  expect_setequal(names(serial$messages), c("first message",
                                            "second message",
                                            "too few distinct cases"))
  # Every refit runs in a worker: a count kept in this process stays at
  # zero.
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  counted <- function(dd) {
    calls$n <- calls$n + 1L
    fn(dd)
  }
  expect_identical(DMAR:::.boot_refit_errors(counted, d, sets, "multicore",
                                             2L), serial)
  expect_identical(calls$n, 0L)
  cl <- parallel::makeForkCluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  expect_identical(DMAR:::.boot_refit_errors(fn, d, sets, "snow", 2L, cl),
                   serial)
  expect_identical(DMAR:::.boot_refit_errors(fn, d, list()),
                   list(messages = integer(), too_few_cases = 0L))
  # A forked worker that dies returns no record, and its replications are
  # refit in this process, so nothing is lost from the counts.
  parent <- Sys.getpid()
  dies <- function(dd) {
    if (Sys.getpid() != parent) tools::pskill(Sys.getpid(), tools::SIGKILL)
    fn(dd)
  }
  expect_identical(suppressWarnings(
    DMAR:::.boot_refit_errors(dies, d, sets, "multicore", 2L)), serial)
  # Through the engine, on the bootstrap's replications and on the
  # jackknife's (fewer replications returned than cases, so the BCa
  # acceleration comes from the jackknife), a seed gives the same counts
  # and limits on one process and on two.
  quiet_boot <- function(ncpus) {
    suppressWarnings(bootstrap_for_mechanics(d, fn, B = 25L,
                                             conf_level = 0.95, kind = "bca",
                                             seed = 113, ncpus = ncpus))
  }
  r1 <- quiet_boot(1L)
  r2 <- quiet_boot(2L)
  expect_gt(r1$bootstrap$reasons[["error"]], 0L)
  expect_gt(r1$bootstrap$too_few_cases, 0L)
  expect_gt(r1$bootstrap$jackknife$reasons[["error"]], 0L)
  expect_identical(r1$bootstrap$acceleration_method, "jackknife")
  expect_identical(r2$bootstrap$ncpus, 2L)
  r2$bootstrap$ncpus <- 1L
  expect_identical(r2, r1)
})


# ---- Several statistics from one set of replications -----------------------
#
# Since 1.1.0 a point function may return several statistics, named by
# stat_names, as ci_bifactor_indices() does for the bifactor indices: each
# interval is boot's own on its column of the replicate matrix, and a
# replication is used or set aside whole, so every interval rests on the
# same replications.

two_fn <- function(d) {
  S <- stats::cov(d)
  c(alpha_fn(d), mean(S[upper.tri(S)]))
}

test_that("several statistics each get boot's own interval on their column", {
  bo <- withr::with_seed(113, boot::boot(
    items, function(d, i) two_fn(d[i, , drop = FALSE]), R = 200, stype = "i"))
  p <- bootstrap_for_mechanics(items, two_fn, B = 200, conf_level = 0.95,
                            kind = "percentile", seed = 113,
                            stat_names = c("alpha", "mean_cov"))
  b <- bootstrap_for_mechanics(items, two_fn, B = 200, conf_level = 0.95,
                            kind = "bca", seed = 113,
                            stat_names = c("alpha", "mean_cov"))
  for (j in 1:2) {
    expect_identical(c(p$lower[j], p$upper[j]),
                     boot::boot.ci(bo, conf = 0.95, type = "perc",
                                   index = j)$percent[4:5])
    expect_identical(c(b$lower[j], b$upper[j]),
                     boot::boot.ci(bo, conf = 0.95, type = "bca",
                                   index = j)$bca[4:5])
    expect_identical(p$se[j], stats::sd(bo$t[, j], na.rm = TRUE))
  }
  expect_identical(p$ci_method, c("percentile", "percentile"))
  expect_identical(b$ci_method, c("bca", "bca"))
  expect_length(b$bootstrap$acceleration, 2L)
  expect_identical(b$bootstrap$acceleration_method,
                   c("regression", "regression"))
  t_ref <- bo$t
  dimnames(t_ref) <- list(NULL, c("alpha", "mean_cov"))
  expect_identical(p$bootstrap$t, t_ref)
  # The first statistic alone is the single statistic run, bit for bit.
  one <- bootstrap_for_mechanics(items, alpha_fn, B = 200, conf_level = 0.95,
                              kind = "bca", seed = 113)
  expect_identical(c(one$lower, one$upper, one$se),
                   c(b$lower[1], b$upper[1], b$se[1]))
})

test_that("a replication with one statistic missing is set aside whole", {
  # A resample without row 1 has no second statistic, so the replication
  # is set aside for both, and both intervals rest on the rest.
  fn <- function(d) {
    S <- stats::cov(d)
    c(alpha_fn(d), if (any(d[, 1] == 1)) mean(S[upper.tri(S)]) else NA_real_)
  }
  w <- testthat::capture_warnings(
    r <- bootstrap_for_mechanics(marked, fn, B = 40, conf_level = 0.8,
                              kind = "percentile", seed = 113,
                              stat_names = c("alpha", "mean_cov")))
  bs <- r$bootstrap
  aside <- bs$status == "non_finite"
  expect_gt(sum(aside), 0L)
  expect_true(all(is.na(bs$t[aside, ])))
  expect_true(all(is.finite(bs$t[!aside, ])))
  expect_identical(bs$returned, sum(!aside))
  expect_match(w, "the eligible intervals are computed from the", all = FALSE)
})

test_that("the jackknife runs once for every statistic, and a BCa fallback names its statistic", {
  # The second statistic, the share of distinct rows, is one in the
  # sample and in every leave-one-out sample, so its influence values
  # have no spread and its BCa interval falls back to the percentile
  # interval; with fewer replications than cases the acceleration of the
  # first comes from the same, single jackknife.
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  fn <- function(d) {
    calls$n <- calls$n + 1L
    c(alpha_fn(d), nrow(unique(d)) / nrow(d))
  }
  w <- testthat::capture_warnings(
    r <- bootstrap_for_mechanics(items, fn, B = 40, conf_level = 0.8,
                              kind = "bca", seed = 113,
                              stat_names = c("alpha", "distinct")))
  expect_identical(calls$n, 1L + 40L + nrow(items))
  expect_identical(r$bootstrap$jackknife$attempted, nrow(items))
  expect_identical(r$ci_method[1], "bca")
  expect_identical(r$bootstrap$acceleration_method, c("jackknife", NA))
  expect_match(r$ci_method[2], "^percentile \\(BCa acceleration could not be estimated")
  expect_match(w, "The BCa interval for distinct could not be formed",
               all = FALSE)
})

test_that("a fit that did not converge is counted under that name", {
  fn <- function(d) {
    if (!any(d[, 1] == 1)) {
      out <- NA_real_
      attr(out, "dmar_reason") <- "not_converged"
      return(out)
    }
    alpha_fn(d)
  }
  w <- testthat::capture_warnings(
    r <- bootstrap_for_mechanics(marked, fn, B = 100, conf_level = 0.95,
                              kind = "percentile", seed = 113))
  bs <- r$bootstrap
  expect_gt(bs$not_converged, 0L)
  expect_identical(sum(bs$status == "not_converged"), bs$not_converged)
  expect_identical(bs$reasons[["non_finite"]], bs$not_converged)
  expect_match(w, sprintf("\\(%d did not converge\\)", bs$not_converged),
               all = FALSE)
  x <- data.frame(term = "estimate", value = 0.8)
  attr(x, "bootstrap") <- bs
  expect_match(DMAR:::.bootstrap_footer_line(x),
               sprintf("set aside: %d did not converge\\)", bs$not_converged))
  expect_identical(
    DMAR:::.bootstrap_reason_phrase(c(error = 0L, non_finite = 6L),
                                    integer(), 0L, 2L, 3L),
    "2 gave an improper (Heywood) solution, 3 did not converge, 1 returned a non-finite value")
})

test_that("no interval is formed from fewer than half of the replications", {
  skip_if_not_installed("boot")
  # A resample returns no value about seven times in ten here, as a
  # categorical omega resample does when it leaves a polychoric
  # correlation undetermined. Through 1.0.1 the engine formed an interval
  # from any two replications that returned a value; since 1.1.0 it asks
  # for half of them, for every caller, the criterion the point estimate
  # of categorical omega under multiple imputation applies to its
  # completed data sets. (The first versions of 1.1.0 kept a switch,
  # require_half = FALSE, that no caller used; it is gone, and the
  # documented element bootstrap$require_half is always TRUE.)
  d <- data.frame(x = seq_len(30))
  fn <- function(dd) {
    if (sum(dd$x) %% 10 < 7) {
      structure(NA_real_, dmar_reason = "undetermined_polychoric")
    } else {
      mean(dd$x)
    }
  }
  expect_false("require_half" %in% names(formals(DMAR:::.bootstrap_ci)))
  w <- testthat::capture_warnings(
    r <- bootstrap_for_mechanics(d, fn, B = 60L, conf_level = 0.95,
                       kind = "percentile", seed = 113))
  expect_lt(r$bootstrap$returned, 30L)
  expect_gte(r$bootstrap$returned, 2L)
  expect_identical(r$bootstrap$bootstrap_status, "insufficient_share")
  expect_true(is.na(r$se) && is.na(r$lower) && is.na(r$upper))
  expect_true(is.na(r$replicate_mean))
  expect_true(r$bootstrap$require_half)
  expect_match(w, "fewer than half returned a value, so no interval is reported",
               fixed = TRUE, all = FALSE)
  # No caution about limits that were not formed.
  expect_false(any(grepl("Caution:", w, fixed = TRUE)))
  line <- .bootstrap_footer_line(structure(data.frame(),
                                           bootstrap = r$bootstrap))
  expect_match(line, "fewer than half returned a value, so no interval is reported",
               fixed = TRUE)
  # With half or more returning a value the interval is formed from the
  # replications that did, and the caution is given.
  fn2 <- function(dd) {
    if (sum(dd$x) %% 10 < 3) {
      structure(NA_real_, dmar_reason = "undetermined_polychoric")
    } else {
      mean(dd$x)
    }
  }
  w2 <- testthat::capture_warnings(
    b <- bootstrap_for_mechanics(d, fn2, B = 60L, conf_level = 0.95,
                                 kind = "percentile", seed = 113))
  expect_gte(b$bootstrap$returned, 30L)
  expect_lt(b$bootstrap$returned, 60L)
  expect_true(b$bootstrap$require_half)
  expect_true(is.finite(b$se) && is.finite(b$lower) && is.finite(b$upper))
  expect_true(any(grepl("Caution:", w2, fixed = TRUE)))
  # The limits are the percentile limits of the replicates that returned
  # a value, as boot.ci() forms them.
  bo <- withr::with_seed(113, boot::boot(
    d, function(dd, i) as.numeric(fn2(dd[i, , drop = FALSE])), R = 60L,
    stype = "i"))
  expect_identical(sum(is.finite(bo$t)), b$bootstrap$returned)
  expect_equal(c(b$lower, b$upper),
               suppressWarnings(boot::boot.ci(bo, conf = 0.95,
                                              type = "perc"))$percent[4:5],
               tolerance = 1e-12)
  expect_equal(b$se, stats::sd(bo$t, na.rm = TRUE), tolerance = 1e-12)
})

test_that("more than half of the replications set aside for too few distinct cases leaves no interval, and the warning says that this is why", {
  skip_if_not_installed("boot")
  # Regression test. The rule on the number of cases keeps only the
  # resamples that repeat few enough rows, chosen by their composition,
  # and an interval was formed from any two of them (at the rule's bound,
  # from resamples that are the data reordered, an interval of no width).
  # The engine now counts these refusals (class dmar_too_few_cases) and
  # forms no interval when they are more than half of the B, with a
  # sentence of its own. A refusal of another kind in more than half of
  # the B leaves no interval either, by the rule that at least half must
  # return a value, with that rule's sentence.
  d <- data.frame(x = seq_len(30))
  refuse <- function(cls, below) {
    function(dd) {
      if (sum(dd$x) %% 10 < below) {
        stop(structure(class = c(cls, "error", "condition"),
                       list(message = "too few distinct cases", call = NULL)))
      }
      mean(dd$x)
    }
  }
  w0 <- testthat::capture_warnings(
    r0 <- bootstrap_for_mechanics(d, refuse("other_error", 7), B = 60L,
                        conf_level = 0.95, kind = "percentile", seed = 113))
  expect_identical(r0$bootstrap$too_few_cases, 0L)
  expect_true(is.na(r0$lower) && is.na(r0$upper))
  expect_true(is.na(r0$replicate_mean))
  expect_match(w0, "fewer than half returned a value, so no interval is reported",
               fixed = TRUE, all = FALSE)
  expect_false(any(grepl("too few distinct cases for the coefficient", w0,
                         fixed = TRUE)))
  w <- testthat::capture_warnings(
    r <- bootstrap_for_mechanics(d, refuse("dmar_too_few_cases", 7), B = 60L,
                       conf_level = 0.95, kind = "percentile", seed = 113))
  bs <- r$bootstrap
  expect_identical(bs$returned, r0$bootstrap$returned)
  expect_identical(bs$too_few_cases, bs$reasons[["error"]])
  expect_gt(bs$too_few_cases, 30L)
  expect_true(bs$require_half)
  expect_false(any(grepl("fewer than half returned a value", w,
                         fixed = TRUE)))
  expect_true(is.na(r$se) && is.na(r$lower) && is.na(r$upper))
  # No bias is estimated from the minority either.
  expect_true(is.na(r$replicate_mean))
  expect_match(w, sprintf(paste0(
    "; %d of them had too few distinct cases for the coefficient, more ",
    "than half of the 60, so no interval is reported: the %d that returned ",
    "a value are the resamples that repeat the fewest rows"),
    bs$too_few_cases, bs$returned), fixed = TRUE, all = FALSE)
  expect_false(any(grepl("Caution:", w, fixed = TRUE)))
  expect_false(any(grepl("Caution:", w0, fixed = TRUE)))
  line <- .bootstrap_footer_line(structure(data.frame(), bootstrap = bs))
  expect_match(line, paste0("; more than half had too few distinct cases, ",
                            "so no interval is reported$"))
  # The jackknife carries the count too, and fewer than half leave the
  # interval as it was.
  expect_identical(bs$jackknife$too_few_cases, 0L)
  a <- suppressWarnings(bootstrap_for_mechanics(d, refuse("other_error", 3), B = 60L,
                                      conf_level = 0.95, kind = "bca",
                                      seed = 113))
  b <- suppressWarnings(bootstrap_for_mechanics(d, refuse("dmar_too_few_cases", 3),
                                      B = 60L, conf_level = 0.95,
                                      kind = "bca", seed = 113))
  expect_gt(b$bootstrap$too_few_cases, 0L)
  expect_lte(b$bootstrap$too_few_cases, 30L)
  expect_identical(c(a$se, a$lower, a$upper), c(b$se, b$lower, b$upper))
})
