test_that("lin_ccc() reproduces the Lin (1989) formula and decomposition", {
  set.seed(113)
  x <- rnorm(40, 100, 15)
  y <- x + rnorm(40, 2, 5)
  res <- lin_ccc(x, y)
  expect_setequal(res$term,
                   c("ccc", "lower_limit", "upper_limit",
                     "pearson_r", "C_b", "u", "v"))
  ccc <- res$value[res$term == "ccc"]
  r   <- res$value[res$term == "pearson_r"]
  c_b <- res$value[res$term == "C_b"]
  expect_equal(ccc, r * c_b, tolerance = 1e-6)
})

test_that("lin_ccc() CCC = 1 only at perfect agreement (y = x)", {
  set.seed(113); x <- rnorm(40)
  expect_equal(suppressWarnings(lin_ccc(x, x))$value[1],
               1, tolerance = 1e-10)
})

test_that("lin_ccc() detects systematic offset (CCC < r)", {
  set.seed(113); x <- rnorm(40)
  # y = x + 5: Pearson r = 1, CCC < 1 because of mean shift.
  res <- lin_ccc(x, x + 5)
  expect_equal(res$value[res$term == "pearson_r"], 1, tolerance = 1e-10)
  expect_lt(res$value[res$term == "ccc"], 1)
})

test_that("lin_ccc() CI bounds stay in [-1, 1]", {
  set.seed(113); x <- rnorm(20)
  y <- x + rnorm(20, 0, 0.5)
  res <- suppressWarnings(lin_ccc(x, y))
  expect_gte(res$value[res$term == "lower_limit"], -1)
  expect_lte(res$value[res$term == "upper_limit"], 1)
})

test_that("lin_ccc() pins the help-page example on the Lin interval", {
  # Regression anchor for the 1.0.0 fix that removed the inconsistent
  # "king_chinchilli" variance (its interval here was [-0.977, 1.000]).
  # The pinned values agree to ten decimals with the independent
  # implementation DescTools::CCC(method_a, method_b, ci = "z-transform"),
  # which returns est = 0.9279709115, lwr.ci = 0.8706469494,
  # upr.ci = 0.9604286311 for these data (DescTools 0.99.x).
  set.seed(113)
  method_a <- rnorm(40, mean = 100, sd = 15)
  method_b <- method_a + rnorm(40, mean = 2, sd = 5)
  res <- lin_ccc(method_a, method_b)
  expect_equal(res$value[res$term == "ccc"],         0.9279709115,
               tolerance = 1e-8)
  expect_equal(res$value[res$term == "lower_limit"], 0.8706469494,
               tolerance = 1e-8)
  expect_equal(res$value[res$term == "upper_limit"], 0.9604286311,
               tolerance = 1e-8)
})

test_that("lin_ccc() rejects the removed king_chinchilli method", {
  set.seed(113); x <- rnorm(20)
  y <- x + rnorm(20, 0, 0.5)
  expect_error(lin_ccc(x, y, method = "king_chinchilli"))
})

test_that("lin_ccc() interval covers at the nominal rate", {
  # Monte Carlo confirmation of the 1.0.0 interval fix; the fast anchor
  # that stays on CRAN is the pinned help-page example test above.
  skip_on_cran()
  set.seed(4861)
  rho <- 0.5; n <- 50; G <- 2000
  cover <- logical(G); width <- numeric(G)
  for (g in seq_len(G)) {
    x <- rnorm(n)
    y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
    # Equal means and unit variances, so the population CCC equals rho.
    res <- lin_ccc(x, y)
    lo <- res$value[res$term == "lower_limit"]
    hi <- res$value[res$term == "upper_limit"]
    cover[g] <- lo <= rho && rho <= hi
    width[g] <- hi - lo
  }
  # Nominal .95; binomial SE with G = 2000 is about 0.005. The removed
  # default measured 0.9815 coverage and 0.667 mean width in this cell.
  expect_gt(mean(cover), 0.93)
  expect_lt(mean(cover), 0.97)
  expect_lt(mean(width), 0.5)
})
