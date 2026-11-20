## ss_seq_two_phase(): the stopping rule evaluated at the data in hand, its
## agreement with ci_two_phase(), the projection, and its refusals.

.v <- function(res) stats::setNames(res$value, res$term)

.pool <- function(seed = 113, M = 5000) {
  set.seed(seed)
  proxy  <- rnorm(M, mean = 50, sd = 10)
  target <- 20 + 0.6 * proxy + rnorm(M, sd = 5)
  list(proxy = proxy, target = target, order = sample(M))
}

test_that("the current half-width is the half-width of the two-phase interval", {
  p <- .pool()
  first <- p$order[1:60]
  s <- .v(ss_seq_two_phase(p$target[first], p$proxy[first], p$proxy[-first],
                           half_width = 1, alpha_level = 0.05))
  ci <- .v(ci_two_phase(p$target[first], p$proxy[first], p$proxy[-first],
                        conf_level = 0.90))
  expect_equal(s[["half_width_current"]],
               (ci[["ci_upper"]] - ci[["ci_lower"]]) / 2,
               tolerance = 1e-12)
  expect_equal(s[["half_width_target"]], 1)
  expect_equal(s[["n_labeled_current"]], 60)
  expect_equal(s[["n_pool"]], 5000)
})

test_that("the rule stops only at or after m0 and only when the interval is narrow enough", {
  p <- .pool()
  at <- function(k, h, m0 = 30) {
    first <- p$order[1:k]
    .v(ss_seq_two_phase(p$target[first], p$proxy[first], p$proxy[-first],
                        half_width = h, m0 = m0))
  }
  wide <- at(40, h = 0.5)
  expect_gt(wide[["half_width_current"]], 0.5)
  expect_equal(wide[["stop"]], 0)
  narrow <- at(100, h = 1)
  expect_lte(narrow[["half_width_current"]], 1)
  expect_equal(narrow[["stop"]], 1)
  # The same narrow interval does not stop before the minimum is reached.
  early <- at(100, h = 1, m0 = 150)
  expect_equal(early[["half_width_current"]], narrow[["half_width_current"]])
  expect_equal(early[["stop"]], 0)
  # The rule is monotone in the target: a wider target never un-stops.
  expect_equal(at(100, h = 5)[["stop"]], 1)
})

test_that("the projection is the normal-quantile size at the current estimates", {
  p <- .pool()
  first <- p$order[1:40]
  y <- p$target[first]; x <- p$proxy[first]; xu <- p$proxy[-first]
  s <- .v(ss_seq_two_phase(y, x, xu, half_width = 1))
  b <- cov(y, x) / var(x)
  s2_e <- summary(lm(y ~ x))$sigma^2
  M <- length(y) + length(xu)
  room <- (1 / qnorm(0.95))^2 - b^2 * var(c(x, xu)) / M
  expect_equal(s[["n_labeled_projected"]], max(40, ceiling(s2_e / room)))
  expect_gte(s[["n_labeled_projected"]], s[["n_labeled_current"]])
})

test_that("a pool too small for the target gives an NA projection with a warning", {
  p <- .pool(M = 300)
  first <- p$order[1:40]
  expect_warning(
    s <- .v(ss_seq_two_phase(p$target[first], p$proxy[first], p$proxy[-first],
                             half_width = 0.3)),
    "every unit of the pool")
  expect_true(is.na(s[["n_labeled_projected"]]))
  expect_equal(s[["stop"]], 0)
})

test_that("a target constant on the labeled units is an error, not a stop", {
  # A binary target that has been 1 on every unit labeled so far would give
  # an interval of zero width, on which the rule would stop.
  p <- .pool()
  first <- p$order[1:40]
  expect_error(ss_seq_two_phase(rep(1, 40), p$proxy[first], p$proxy[-first],
                                half_width = 0.1),
               "target takes a single value")
})

test_that("a target the proxy reproduces exactly on the labeled units is an error, not a stop", {
  # A binary judge (the proxy) and a binary human label (the target) that
  # agree on all 30 units labeled so far, 21 ones and 9 zeros. The residual
  # variance is zero, the half-width (0.011) is the pool's part alone, and
  # the rule would stop at a target of 0.03 on no evidence of disagreement.
  set.seed(113)
  y <- rep(c(1, 0), c(21, 9))
  x_unlabeled <- rbinom(4970, 1, 0.70)
  e <- tryCatch(ss_seq_two_phase(y, y, x_unlabeled, half_width = 0.03),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, c("y", "x"))
  expect_match(conditionMessage(e), "exact linear function of the proxy")
  expect_match(conditionMessage(e), "label more units")
  # The same when the residuals are rounding error and not exactly zero.
  p <- .pool()
  first <- p$order[1:40]
  linear <- 20 + 0.6 * p$proxy[first]
  expect_gt(sum(stats::resid(stats::lm(linear ~ p$proxy[first]))^2), 0)
  expect_error(ss_seq_two_phase(linear, p$proxy[first], p$proxy[-first],
                                half_width = 1),
               class = "dmar_argument_invalid")
  # Residuals just above rounding error are data, and the rule is evaluated:
  # with a residual variance of 4e-14 of the target's (the limit is
  # .Machine$double.eps, 2.2e-16) the half-width is the pool's part, 0.142,
  # and the rule stops at a target of 1.
  near <- linear + 1e-6 * rnorm(40)
  s_near <- .v(ss_seq_two_phase(near, p$proxy[first], p$proxy[-first],
                                half_width = 1))
  expect_equal(s_near[["stop"]], 1)
  expect_equal(s_near[["half_width_current"]], 0.1424, tolerance = 1e-3)
  # The refusal does not depend on whether the rule would have been met:
  # it is made below m0 (12 agreeing units) and with a pool of 200, whose
  # part of the half-width alone exceeds the target of 0.03.
  few <- c(1:8, 22:25)
  expect_error(ss_seq_two_phase(y[few], y[few], x_unlabeled,
                                half_width = 0.03),
               class = "dmar_argument_invalid")
  e_small <- tryCatch(ss_seq_two_phase(y, y, x_unlabeled[1:170],
                                       half_width = 0.03),
                      error = function(e) e)
  expect_s3_class(e_small, "dmar_argument_invalid")
  expect_match(conditionMessage(e_small), "size of the pool alone")
  # One disagreement among the 30 and the rule is evaluated: it does not
  # stop, and the interval is five times as wide.
  y_one <- y
  y_one[1] <- 0
  s <- .v(ss_seq_two_phase(y_one, y, x_unlabeled, half_width = 0.03))
  expect_equal(s[["stop"]], 0)
  expect_gt(s[["half_width_current"]], 0.05)
})

test_that("the projection leaves at least one unit with the proxy alone", {
  # A pool of 60 with 40 labeled: when the projection needs every unit of
  # the pool, it is past the M - 1 units the interval allows.
  p <- .pool(M = 60)
  first <- p$order[1:40]
  y <- p$target[first]; x <- p$proxy[first]; xu <- p$proxy[-first]
  b <- cov(y, x) / var(x)
  s2_e <- summary(lm(y ~ x))$sigma^2
  # The half-width at which the projection, s2_e / room rounded up, is 60.
  h <- qnorm(0.95) * sqrt(s2_e / 59.5 + b^2 * var(c(x, xu)) / 60)
  expect_warning(s <- .v(ss_seq_two_phase(y, x, xu, half_width = h)),
                 "but one")
  expect_true(is.na(s[["n_labeled_projected"]]))
})

test_that("missing values are set aside as in ci_two_phase()", {
  p <- .pool()
  first <- p$order[1:50]
  y <- p$target[first]; y[c(3, 4)] <- NA
  s <- .v(ss_seq_two_phase(y, p$proxy[first], p$proxy[-first], half_width = 1))
  expect_equal(s[["n_labeled_current"]], 48)
})

test_that("malformed inputs are refused", {
  p <- .pool()
  first <- p$order[1:40]
  call_with <- function(...) ss_seq_two_phase(p$target[first], p$proxy[first],
                                              p$proxy[-first], ...)
  expect_error(call_with(half_width = -1), "half_width")
  expect_error(call_with(half_width = 1, alpha_level = 0.6), "alpha_level")
  expect_error(call_with(half_width = 1, m0 = 2), "m0")
  expect_error(ss_seq_two_phase(p$target[first], rep(1, 40), p$proxy[-first],
                                half_width = 1), "single value")
})
