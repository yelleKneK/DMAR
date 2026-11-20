# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# Tests for reliability_kr20(). KR-20 on dichotomous items is algebraically
# equal to coefficient alpha; we exploit that to cross-check the point
# estimate against the alpha formula on the same data, and we check the
# CI invariants in the same way as for reliability_alpha().

set.seed(113)
N <- 300
J <- 10
ability <- rnorm(N)
loadings <- rep(0.6, J)
latent <- outer(ability, loadings) +
          matrix(rnorm(N * J, sd = sqrt(1 - 0.6^2)), N, J)
items_binary <- (latent > 0) * 1
colnames(items_binary) <- paste0("y", seq_len(J))


test_that("reliability_kr20() equals reliability_alpha() on dichotomous data", {
  r_kr20  <- reliability_kr20 (data = items_binary, ci_method = "none")
  r_alpha <- reliability_alpha(data = items_binary, ci_method = "none")
  expect_equal(
    r_kr20$value [r_kr20$term  == "estimate"],
    r_alpha$value[r_alpha$term == "estimate"],
    tolerance = 1e-10
  )
})

test_that("reliability_kr20() matches the KR-20 formula with N/(N-1) bias correction", {
  # Pair the 1/N-denominator p_j q_j numerator with the 1/N-denominator
  # composite variance so the conventions match; that is the form that
  # equals coefficient alpha computed from cov(data) with Bessel
  # correction. (See note above .kr20_from_data in reliability_internals.R.)
  p <- colMeans(items_binary)
  q <- 1 - p
  Nrow <- nrow(items_binary)
  s2_biased <- var(rowSums(items_binary)) * (Nrow - 1) / Nrow
  manual <- (J / (J - 1)) * (1 - sum(p * q) / s2_biased)
  r <- reliability_kr20(data = items_binary, ci_method = "none")
  expect_equal(r$value[r$term == "estimate"], manual, tolerance = 1e-10)
})

test_that("reliability_kr20() returns the expected tidy data.frame and attrs", {
  r <- reliability_kr20(data = items_binary)
  expect_s3_class(r, "data.frame")
  expect_named(r, c("term", "value"))
  expect_true(all(c("estimate", "se", "ci_lower", "ci_upper",
                    "conf_level", "N", "n_items") %in% r$term))
  expect_equal(attr(r, "coefficient"), "kr20")
  expect_equal(attr(r, "ci_method"), "feldt")
})

test_that("reliability_kr20() CIs bracket the estimate (closed-form methods)", {
  for (m in c("feldt", "fisher", "bonett", "hakstian_whalen", "ml", "ml_logistic")) {
    r <- reliability_kr20(data = items_binary, ci_method = m)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "ci_lower"]
    ul  <- r$value[r$term == "ci_upper"]
    expect_lte(ll, est, label = paste("lower limit <= estimate for", m))
    expect_gte(ul, est, label = paste("upper limit >= estimate for", m))
    expect_gte(ll, 0,   label = paste("lower limit >= 0 for", m))
    expect_lte(ul, 1,   label = paste("upper limit <= 1 for", m))
  }
})

test_that("reliability_kr20() reports both standard error scales for a transformation interval", {
  # The se row is on the coefficient scale (the delta method
  # back-transform of the log-scale standard error); the log-scale value
  # travels in the se_transformed row, named by the se_transform_scale
  # attribute (decided 2026-08-12).
  r <- reliability_kr20(data = items_binary, ci_method = "bonett")
  est  <- r$value[r$term == "estimate"]
  se_t <- r$value[r$term == "se_transformed"]
  expect_identical(attr(r, "se_transform_scale"), "log(1-alpha)")
  expect_equal(se_t, sqrt(2 * J / ((J - 1) * (N - 2))), tolerance = 1e-12)
  expect_equal(r$value[r$term == "se"], se_t * (1 - est),
               tolerance = 1e-12)
  # The Feldt default is an F pivot: no standard error on either scale,
  # and no se_transformed row.
  r_f <- reliability_kr20(data = items_binary)
  expect_true(is.na(r_f$value[r_f$term == "se"]))
  expect_false("se_transformed" %in% r_f$term)
  expect_null(attr(r_f, "se_transform_scale"))
})

test_that("reliability_kr20() bootstrap CIs bracket the estimate", {
  skip_on_cran()
  skip_if_not_installed("boot")
  for (m in c("bootstrap_se", "percentile", "bca")) {
    r <- reliability_kr20(data = items_binary, ci_method = m, B = 200)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "ci_lower"]
    ul  <- r$value[r$term == "ci_upper"]
    expect_lte(ll, est)
    expect_gte(ul, est)
    expect_gte(ll, 0)
    expect_lte(ul, 1)
    expect_equal(attr(r, "B"), 200)
  }
})

test_that("reliability_kr20() higher conf_level gives a wider CI (Feldt)", {
  r80 <- reliability_kr20(data = items_binary, ci_method = "feldt",
                          conf_level = 0.80)
  r95 <- reliability_kr20(data = items_binary, ci_method = "feldt",
                          conf_level = 0.95)
  width80 <- r80$value[r80$term == "ci_upper"] -
             r80$value[r80$term == "ci_lower"]
  width95 <- r95$value[r95$term == "ci_upper"] -
             r95$value[r95$term == "ci_lower"]
  expect_gt(width95, width80)
})

test_that("reliability_kr20() rejects non-binary data", {
  items_continuous <- items_binary + 0.1
  expect_error(reliability_kr20(data = items_continuous),
               regexp = "0/1")
})

test_that("reliability_kr20() rejects 0/1 data with too few items", {
  expect_error(reliability_kr20(data = items_binary[, 1, drop = FALSE]),
               regexp = "At least two items")
})

test_that("reliability_kr20() errors when data is missing", {
  expect_error(reliability_kr20(),
               regexp = "'data'.*is required")
})

test_that("reliability_kr20() rejects out-of-range conf_level", {
  expect_error(reliability_kr20(data = items_binary, conf_level = 1),
               regexp = "conf_level")
})

test_that("reliability_kr20() takes ncpus and reports the effective count", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  r <- reliability_kr20(data = items_binary, ci_method = "percentile", B = 50,
                        seed = 113)
  expect_identical(attr(r, "B_used"), 50L)
  expect_identical(attr(r, "bootstrap")$improper_kept, 0L)
  expect_error(reliability_kr20(data = items_binary, ncpus = 1.5),
               "'ncpus' must be a single whole number")
})

test_that("a negative KR-20 is warned about and its interval is not inverted", {
  # Regression test (the review's case): one item is nearly the reverse of
  # another, so KR-20 is -0.821 and Feldt's interval lies wholly below
  # zero. Through 1.0.1 the result was [0, -0.171], the lower limit above
  # the upper; each limit is now clamped into [0, 1].
  set.seed(3)
  z <- stats::rbinom(60, 1, .5)
  Xb <- cbind(i1 = z, i2 = 1 - z, i3 = stats::rbinom(60, 1, .5),
              i4 = stats::rbinom(60, 1, .5))
  Xb[1:6, 2] <- z[1:6]
  expect_warning(r <- reliability_kr20(Xb), class = "dmar_negative_reliability")
  v <- stats::setNames(r$value, r$term)
  expect_equal(v[["estimate"]], -0.8213478, tolerance = 1e-6)
  expect_identical(unname(v[c("ci_lower", "ci_upper")]), c(0, 0))
  expect_warning(reliability_kr20(Xb), "The estimate of KR-20 is negative",
                 fixed = TRUE)
})

test_that("identical binary items give a zero-width bootstrap interval, not an error", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # KR-20 is 1 in every resample in which the item varies, and a resample
  # in which it does not is refused, as the data would be, and set aside
  # (it returned NA before the check reached the replications). Through
  # 1.0.1 boot's test for equal replicates then met if (NA) and the call
  # stopped with "missing value where TRUE/FALSE needed".
  b <- c(rep(0, 5), rep(1, 3))
  Xb <- cbind(b, b, b)
  for (m in c("percentile", "bca")) {
    w <- testthat::capture_warnings(
      r <- reliability_kr20(Xb, ci_method = m, B = 200, seed = 113))
    expect_equal(r$value[r$term %in% c("ci_lower", "ci_upper")],
                 c(1, 1))
    expect_true(any(grepl("gave the same estimate", w, fixed = TRUE)))
  }
  expect_gt(attr(r, "bootstrap")$failed, 0L)
})

test_that("reliability_kr20() validates B when a bootstrap is asked for", {
  expect_error(reliability_kr20(items_binary, ci_method = "percentile",
                                B = 1),
               "'B' must be a single whole number of at least 2",
               fixed = TRUE)
  expect_error(reliability_kr20(items_binary, ci_method = "bca", B = 2.5),
               "'B' must be a single whole number of at least 2",
               fixed = TRUE)
})


test_that("reliability_kr20()'s ml standard error is van Zyl's variance on the covariance matrix", {
  # Binary items of unequal difficulty have unequal variances, so the
  # variance evaluated on cov2cor(S), as before 1.1.0, was not that of KR-20
  # of the raw items.
  b <- c(-1.2, -.8, -.4, 0, .4, .8, 1.2, 1.6, -1.6, .2)
  items <- (latent > matrix(b, N, J, byrow = TRUE)) * 1
  S <- stats::cov(items)
  g <- sum(S)
  vz <- J^2 / (J - 1)^2 * 2 / g^3 *
    (g * (sum(diag(S %*% S)) + sum(diag(S))^2) -
       2 * sum(diag(S)) * sum(S %*% S))
  r <- reliability_kr20(data = items, ci_method = "ml")
  expect_equal(r$value[r$term == "se"], sqrt(vz / (N - 1)), tolerance = 1e-12)
})


test_that("listwise deletion that leaves fewer than two rows stops KR-20 before anything is computed", {
  skip_on_cran()
  # Regression test. Through 1.0.1 a 0/1 matrix in which every row misses
  # an item (each respondent saw all but one, as in a planned missingness
  # design) returned NA for the estimate and every limit, with two
  # warnings that NaNs were produced. The call now stops with the family's
  # message (.stop_relia_too_few_complete_rows()), before any interval,
  # and, since reliability_kr20() has no missing argument, offers full
  # information maximum likelihood through reliability_alpha().
  e <- items_binary
  for (i in seq_len(nrow(e))) e[i, (i %% ncol(e)) + 1L] <- NA
  msg0 <- paste0("No row of 'data' is observed on every item, so listwise ",
                 "deletion (the only treatment reliability_kr20() offers) ",
                 "leaves nothing to analyze. Use reliability_alpha() with ",
                 "missing = \"fiml\", which computes the same coefficient ",
                 "(alpha of 0/1 items is KR-20) from the incomplete rows as ",
                 "well and assumes the values missing at random, where ",
                 "listwise deletion assumes them missing completely at ",
                 "random.")
  for (m in c("feldt", "ml", "percentile")) {
    err <- tryCatch(
      withCallingHandlers(reliability_kr20(e, ci_method = m, B = 20),
                          warning = function(w) stop(w)),
      error = conditionMessage)
    expect_identical(err, msg0)
  }
  # The wrapper reaches the same refusal.
  expect_error(reliability(e, type = "kr20"), msg0, fixed = TRUE)
  # One complete row is no better.
  e1 <- e
  e1[1, ] <- items_binary[1, ]
  expect_error(reliability_kr20(e1),
               "Only one row of 'data' is observed on every item",
               fixed = TRUE)
  # The function the message offers computes KR-20: on complete data the
  # two agree, and on these data it uses every row.
  expect_equal(
    suppressMessages(reliability_alpha(items_binary, missing = "fiml",
                                       ci_method = "none"))$value[1],
    reliability_kr20(items_binary, ci_method = "none")$value[1],
    tolerance = 1e-10)
  fiml <- suppressMessages(reliability_alpha(e, missing = "fiml",
                                             ci_method = "none"))
  expect_true(is.finite(fiml$value[fiml$term == "estimate"]))
  expect_identical(fiml$value[fiml$term == "N"], as.numeric(nrow(e)))
  # FIML is not offered when a pair of items is never observed together,
  # and an item no row observes is named first.
  g <- items_binary[1:60, ]
  g[1:30, "y1"] <- NA
  g[31:60, "y2"] <- NA
  expect_error(reliability_kr20(g),
               paste0("Full information maximum likelihood ",
                      "(reliability_alpha() with missing = \"fiml\") cannot ",
                      "help either: no row observes both y1 and y2"),
               fixed = TRUE)
  h <- items_binary[1:60, ]
  h[, "y4"] <- NA
  expect_error(reliability_kr20(h),
               "No row observes y4, so neither its variance", fixed = TRUE)
  # Two complete rows reach the checks of the items: an item the two rows
  # answer alike takes a single value in them and is named (through 1.0.1
  # KR-20 was computed with it). The binary and item count checks come
  # first, as before.
  e2 <- e1
  e2[2, ] <- items_binary[2, ]
  expect_error(reliability_kr20(e2), "a single value in the 2 complete rows",
               fixed = TRUE)
  expect_error(reliability_kr20(e1[, 1, drop = FALSE]),
               "At least two items are required.", fixed = TRUE)
  e1[1, 1] <- 2
  expect_error(reliability_kr20(e1), "non-binary values were found",
               fixed = TRUE)
})


test_that("a bootstrap replication that has lost a rare item's only passer is set aside with the refusal", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # Regression test. An item that takes a single value is refused in the
  # data, since KR-20 would count it among the J items (J(J - 2)/(J - 1)^2
  # times the KR-20 of the same total without it), and a resample that
  # leaves out the one respondent who passed an item is such data.
  # Through version 1.0.1 KR-20 of that resample was computed with the
  # item and kept: with these data 73 of 200 resamples lost the passer,
  # each kept at 0.96 = 6 * 4 / 25 times the KR-20 of the other five items.
  # They are now set aside as errors that quote the refusal, and the
  # caution reports the share.
  set.seed(1)
  n <- 60
  f <- stats::rnorm(n)
  b <- as.data.frame(sapply(c(0, .3, -.2, .1, .5), function(th) {
    as.numeric(f + stats::rnorm(n) > th)
  }))
  names(b) <- paste0("i", 1:5)
  b$i6 <- 0
  b$i6[which.max(f)] <- 1
  w <- testthat::capture_warnings(
    r <- reliability_kr20(b, ci_method = "percentile", B = 200, seed = 3))
  bt <- attr(r, "bootstrap")
  # boot draws the resamples in the order the engine does, so the same
  # seed finds the resamples without the passer.
  set.seed(3)
  bb <- boot::boot(b, function(d, i) sum(d$i6[i]), R = 200)
  lost <- bb$t[, 1] == 0
  expect_identical(sum(lost), 73L)
  expect_identical(bt$failed, 73L)
  expect_identical(bt$reasons[["error"]], 73L)
  expect_identical(as.character(bt$status[lost]), rep("error", 73L))
  expect_identical(
    names(bt$error_messages),
    paste0("i6 takes a single value in this resample, so its variance is ",
           "zero, and an item that does not vary is refused."))
  expect_true(bt$caution)
  expect_match(w, "73 of 200 bootstrap replications were set aside",
               all = FALSE, fixed = TRUE)
  # The replications kept are the KR-20 of their resamples.
  idx <- boot::boot.array(bb, indices = TRUE)
  keep <- which(!lost)[1:3]
  expect_equal(bt$t[keep], vapply(keep, function(k) {
    reliability_kr20(b[idx[k, ], ], ci_method = "none")$value[1]
  }, numeric(1)), tolerance = 1e-12)
  # The jackknife deletion of the passer is refused the same way.
  w <- testthat::capture_warnings(
    r <- reliability_kr20(b, ci_method = "bca", B = 50, seed = 3))
  jk <- attr(r, "bootstrap")$jackknife
  expect_identical(jk$attempted, 60L)
  expect_identical(jk$reasons[["error"]], 1L)
})

