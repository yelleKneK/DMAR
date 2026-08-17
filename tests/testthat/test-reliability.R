# Tests for the general reliability() wrapper: dispatch correctness,
# auto-detection of type, and consistency with the family functions.

skip_if_not_installed("lavaan")

set.seed(113)
J <- 6
N <- 250
loadings <- seq(0.4, 0.8, length.out = J)
eta <- rnorm(N)
errors <- matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
items_continuous <- sweep(matrix(rep(eta, J), N, J), 2, loadings, `*`) + errors
colnames(items_continuous) <- paste0("y", seq_len(J))
S_continuous <- cov(items_continuous)

items_binary <- (items_continuous > 0) * 1
items_categorical <- apply(items_continuous, 2, function(x)
  as.integer(cut(x, breaks = quantile(x, probs = seq(0, 1, length.out = 6)),
                 include.lowest = TRUE)))
colnames(items_categorical) <- colnames(items_continuous)


test_that("reliability() dispatches to reliability_alpha when type = 'alpha'", {
  r_wrap   <- reliability      (data = items_continuous, type = "alpha",
                                ci_method = "bonett")
  r_direct <- reliability_alpha(data = items_continuous,
                                ci_method = "bonett")
  expect_equal(r_wrap$value, r_direct$value)
  expect_equal(attr(r_wrap, "coefficient"), "alpha")
  # The transformation-scale standard error row and its scale label pass
  # through the wrapper unchanged.
  expect_true("se_transformed" %in% r_wrap$term)
  expect_identical(attr(r_wrap, "se_transform_scale"),
                   attr(r_direct, "se_transform_scale"))
})

test_that("reliability() dispatches to reliability_kr20 when type = 'kr20'", {
  r_wrap   <- reliability     (data = items_binary, type = "kr20")
  r_direct <- reliability_kr20(data = items_binary)
  expect_equal(r_wrap$value, r_direct$value)
  expect_equal(attr(r_wrap, "coefficient"), "kr20")
})

test_that("reliability() dispatches to reliability_omega when type = 'omega'", {
  r_wrap   <- reliability      (data = items_continuous, type = "omega",
                                denominator = "model_implied",
                                ci_method = "ml")
  r_direct <- reliability_omega(data = items_continuous,
                                denominator = "model_implied",
                                ci_method = "ml")
  expect_equal(r_wrap$value, r_direct$value)
  expect_equal(attr(r_wrap, "coefficient"), "omega")
})

test_that("reliability() forwards 'denominator' to reliability_omega()", {
  r_wrap   <- reliability      (data = items_continuous, type = "omega",
                                denominator = "observed",
                                ci_method = "none")
  r_direct <- reliability_omega(data = items_continuous,
                                denominator = "observed",
                                ci_method = "none")
  expect_equal(r_wrap$value, r_direct$value)
  expect_equal(attr(r_wrap, "coefficient"), "omega")
  expect_equal(attr(r_wrap, "denominator"), "observed")
})

test_that("reliability() rejects 'denominator' for non-omega types", {
  expect_error(
    reliability(data = items_continuous, type = "alpha",
                denominator = "observed"),
    regexp = "applies only to type"
  )
})

test_that("reliability() rejects the retired type = 'omega_h'", {
  expect_error(
    reliability(data = items_continuous, type = "omega_h"),
    regexp = "arg"
  )
})

test_that("reliability() auto-detects robust omega for covariance matrix input", {
  msgs <- testthat::capture_messages(r <- reliability(S = S_continuous, N = N))
  expect_true(any(grepl("type = \"omega\"", msgs)))
  expect_true(any(grepl("covariance matrix", msgs)))
  expect_equal(attr(r, "coefficient"), "omega")
  expect_equal(attr(r, "denominator"), "observed")
  # Robust omega intervals are bootstrap based and need raw data, so the
  # covariance-input default is the point estimate.
  expect_equal(attr(r, "ci_method"), "none")
})

test_that("reliability() falls back to 'ml' for model implied omega with S only", {
  suppressMessages(
    r <- reliability(S = S_continuous, N = N, type = "omega",
                     denominator = "model_implied")
  )
  expect_equal(attr(r, "ci_method"), "ml")
})

test_that("reliability() auto-detects 'omega' for continuous data", {
  expect_message(
    r <- reliability(data = items_continuous, ci_method = "none"),
    regexp = "type = \"omega\""
  )
  expect_equal(attr(r, "coefficient"), "omega")
  expect_equal(attr(r, "denominator"), "observed")
})

test_that("reliability() auto-detects 'omega_categorical' for ordered categorical data", {
  expect_message(
    r <- reliability(data = items_categorical, ci_method = "none"),
    regexp = "type = \"omega_categorical\""
  )
  expect_equal(attr(r, "coefficient"), "omega_categorical")
})

test_that("reliability() accepts 'omega_c' as a type shorthand", {
  r_word  <- reliability(data = items_categorical,
                         type = "omega_categorical", ci_method = "none")
  r_short <- reliability(data = items_categorical,
                         type = "omega_c", ci_method = "none")
  expect_equal(r_word$value, r_short$value)
})

test_that("reliability() forwards ci_method to the chosen family function", {
  r <- reliability(data = items_continuous, type = "alpha",
                   ci_method = "feldt")
  expect_equal(attr(r, "ci_method"), "feldt")
})

test_that("reliability() respects the family-specific default when ci_method = NULL", {
  r <- reliability(data = items_continuous, type = "alpha")
  expect_equal(attr(r, "ci_method"), "bonett")
})

test_that("reliability() rejects S input for type = 'kr20' and 'omega_categorical'", {
  expect_error(reliability(S = S_continuous, N = N, type = "kr20"),
               regexp = "requires raw 'data'")
  expect_error(reliability(S = S_continuous, N = N, type = "omega_categorical"),
               regexp = "requires raw 'data'")
})

test_that("reliability() rejects unknown type", {
  expect_error(reliability(data = items_continuous, type = "lambda4"))
})

test_that("reliability() errors when neither data nor S is supplied (no type)", {
  expect_error(reliability(),
               regexp = "Either 'data' or 'S' must be supplied")
})
