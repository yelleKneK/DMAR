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
  expect_true(all(c("estimate", "se", "lower_limit", "upper_limit",
                    "conf_level", "N", "J") %in% r$term))
  expect_equal(attr(r, "coefficient"), "kr20")
  expect_equal(attr(r, "ci_method"), "feldt")
})

test_that("reliability_kr20() CIs bracket the estimate (closed-form methods)", {
  for (m in c("feldt", "fisher", "bonett", "hakstian_whalen", "ml", "ml_logistic")) {
    r <- reliability_kr20(data = items_binary, ci_method = m)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "lower_limit"]
    ul  <- r$value[r$term == "upper_limit"]
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
  skip_if_not_installed("boot")
  for (m in c("bootstrap_se", "percentile", "bca")) {
    r <- reliability_kr20(data = items_binary, ci_method = m, B = 200)
    est <- r$value[r$term == "estimate"]
    ll  <- r$value[r$term == "lower_limit"]
    ul  <- r$value[r$term == "upper_limit"]
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
  width80 <- r80$value[r80$term == "upper_limit"] -
             r80$value[r80$term == "lower_limit"]
  width95 <- r95$value[r95$term == "upper_limit"] -
             r95$value[r95$term == "lower_limit"]
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
