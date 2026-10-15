# Tests for reliability_omega_categorical(). Point estimate is checked for being in
# (0, 1) and close to the population value used in the simulation. CI
# methods are the bootstrap families. Computations involve a categorical
# CFA via lavaan which is heavy; bootstrap reps are kept small here.

skip_if_not_installed("lavaan")
skip_if_not_installed("mvtnorm")

set.seed(113)
N <- 400
J <- 5
loadings <- rep(0.7, J)
eta <- rnorm(N)
latent <- outer(eta, loadings) +
          matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
items_categorical <- apply(latent, 2, function(x)
  as.integer(cut(x, breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf),
                 labels = FALSE)))
colnames(items_categorical) <- paste0("y", seq_len(J))


test_that("reliability_omega_categorical() matches the MBESS reference implementation", {
  # One ordered-categorical CFA via lavaan with WLSMV, about a quarter of a
  # second, and the same fit runs unguarded further down this file, so there
  # is nothing to save by skipping it on CRAN.
  #
  # The assertion is an equality against MBESS, not a window. A window wide
  # enough to survive sampling variability (the old one spanned 0.6 to 0.95)
  # also admits at least three different coefficients in place of the
  # Green and Yang categorical omega: coefficient alpha on the raw integer
  # scores, linear omega on the raw scores, and alpha on the polychoric
  # matrix. Anchoring to the reference implementation is what makes this a
  # test of the coefficient rather than of its order of magnitude.
  r <- reliability_omega_categorical(data = items_categorical, ci_method = "none")
  est <- r$value[r$term == "estimate"]
  expect_gt(est, 0)
  expect_lt(est, 1)

  # Pinned from MBESS::ci.reliability (MBESS 4.9.3, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  expect_equal(est, 0.8106007039872396, tolerance = 1e-6)
})

test_that("reliability_omega_categorical() returns the expected tidy data.frame and attrs", {
  r <- reliability_omega_categorical(data = items_categorical, ci_method = "none")
  expect_s3_class(r, "data.frame")
  expect_named(r, c("term", "value"))
  expect_true(all(c("estimate", "se", "lower_limit", "upper_limit",
                    "conf_level", "N", "J") %in% r$term))
  expect_equal(attr(r, "coefficient"), "omega_categorical")
  expect_equal(attr(r, "ci_method"), "none")
  # The interval menu is bootstrap only, and a bootstrap standard error
  # is already on the coefficient scale, so no se_transformed row (or
  # se_transform_scale attribute) ever appears here.
  expect_false("se_transformed" %in% r$term)
  expect_null(attr(r, "se_transform_scale"))
})

test_that("reliability_omega_categorical() recommends bca but never bootstraps by default", {
  # "bca" leads the menu as the recommended method (KP 2016), but the
  # default call reports the point estimate with a message; no bootstrap
  # runs unless requested.
  expect_equal(eval(formals(reliability_omega_categorical)$ci_method)[1], "bca")
  expect_message(
    r <- reliability_omega_categorical(data = items_categorical),
    regexp = "bootstrap based"
  )
  expect_equal(attr(r, "ci_method"), "none")
  expect_true(is.na(r$value[r$term == "lower_limit"]))
})

test_that("reliability_omega_categorical() BCa bootstrap CI brackets the estimate", {
  # Skipped on CRAN: B = 60 ordered-categorical CFA refits via lavaan
  # WLSMV takes roughly 30-40 s on CRAN's check machines.
  skip_on_cran()
  skip_if_not_installed("boot")
  r <- reliability_omega_categorical(data = items_categorical, ci_method = "bca",
                           B = 60)
  est <- r$value[r$term == "estimate"]
  ll  <- r$value[r$term == "lower_limit"]
  ul  <- r$value[r$term == "upper_limit"]
  expect_lte(ll, est)
  expect_gte(ul, est)
  expect_gte(ll, 0)
  expect_lte(ul, 1)
  expect_equal(attr(r, "B"), 60)
})

test_that("reliability_omega_categorical() percentile bootstrap CI brackets the estimate", {
  # Same reason as the BCa block above: heavy lavaan WLSMV inside the
  # bootstrap loop. Skipped on CRAN.
  skip_on_cran()
  skip_if_not_installed("boot")
  r <- reliability_omega_categorical(data = items_categorical, ci_method = "percentile",
                           B = 60)
  est <- r$value[r$term == "estimate"]
  ll  <- r$value[r$term == "lower_limit"]
  ul  <- r$value[r$term == "upper_limit"]
  expect_lte(ll, est)
  expect_gte(ul, est)
  expect_gte(ll, 0)
  expect_lte(ul, 1)
})

test_that("reliability_omega_categorical() rejects non-integer (continuous) data", {
  expect_error(
    reliability_omega_categorical(data = items_categorical + 0.1),
    regexp = "integer category"
  )
})

test_that("reliability_omega_categorical() errors when data is missing", {
  expect_error(reliability_omega_categorical(), regexp = "'data'.*is required")
})

test_that("reliability_omega_categorical() rejects too few items", {
  expect_error(
    reliability_omega_categorical(data = items_categorical[, 1, drop = FALSE]),
    regexp = "At least two items"
  )
})

test_that("reliability_omega_categorical() rejects out-of-range conf_level", {
  expect_error(
    reliability_omega_categorical(data = items_categorical, conf_level = -0.1),
    regexp = "conf_level"
  )
})


test_that("a leave-one-out refit that fails no longer stops the BCa interval", {
  # The 2026-10-14 report: with B below n, boot estimates the BCa
  # acceleration by the jackknife, and one leave-one-out refit that
  # returned NA made the acceleration NA and stopped the call with
  # "estimated adjustment 'a' is NA". Here the first of six binary items
  # (n = 60) is endorsed by a single respondent, so the refit without
  # that case refuses a constant item: exactly one jackknife replicate
  # fails, and the resamples that omit the case fail the same way. The
  # acceleration is now estimated from the refits that returned a value,
  # the counts are on the result, and the interval is BCa. Skipped on
  # CRAN: n + B ordered-categorical refits.
  skip_on_cran()
  skip_if_not_installed("boot")
  set.seed(113)
  n_bad <- 60
  j_bad <- 6
  eta_bad <- rnorm(n_bad)
  items_bad <- (outer(eta_bad, rep(0.7, j_bad)) +
                  matrix(rnorm(n_bad * j_bad), n_bad, j_bad) * sqrt(1 - 0.7^2)) > 0
  items_bad <- apply(items_bad, 2, as.integer)
  colnames(items_bad) <- paste0("y", seq_len(j_bad))
  items_bad[, 1] <- 0L
  items_bad[1, 1] <- 1L
  # lavaan prints its variable table when it refuses the constant item
  # and may warn about an improper solution on a sparse resample; the
  # capture and the handler keep both out of the test record, which
  # looks only for the package's own warning.
  w <- testthat::capture_warnings(invisible(utils::capture.output(
    r <- withCallingHandlers(
      reliability_omega_categorical(data = items_bad, ci_method = "bca",
                                    B = 30, conf_level = 0.80, seed = 113),
      warning = function(w) {
        if (grepl("lavaan", conditionMessage(w))) invokeRestart("muffleWarning")
      })
  )))
  # Two thirds of the resamples omit case 1 or hit a degenerate table and
  # are set aside, so boot adds its own warning that extreme order
  # statistics served as endpoints; the package warning is the one that
  # names the jackknife count, and it is the only other one.
  own <- grepl("1 of 60 jackknife replicates behind the BCa acceleration failed", w)
  expect_identical(sum(own), 1L)
  expect_match(w[own], "bootstrap replications were set aside")
  expect_true(all(grepl("extreme order statistics used as endpoints", w[!own])))
  expect_identical(attr(r, "ci_method"), "bca")
  bs <- attr(r, "bootstrap")
  expect_identical(bs$attempted, 30L)
  expect_identical(bs$returned + bs$failed, 30L)
  expect_identical(bs$jackknife$attempted, 60L)
  expect_identical(bs$jackknife$returned, 59L)
  expect_identical(bs$jackknife$failed, 1L)
  expect_identical(bs$acceleration_method, "jackknife")
  expect_true(is.finite(bs$acceleration))
  est <- r$value[r$term == "estimate"]
  ll  <- r$value[r$term == "lower_limit"]
  ul  <- r$value[r$term == "upper_limit"]
  expect_lte(ll, est)
  expect_gte(ul, est)
  expect_gte(ll, 0)
  expect_lte(ul, 1)
})
