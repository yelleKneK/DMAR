# Tests for vargha_delaney_A() against the Mann-Whitney U formulation and
# the DeLong et al. (1988) variance.

test_that("vargha_delaney_A() recovers A from the rank/U formulation", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(3, 4, 5, 6, 7)
  res <- vargha_delaney_A(x, y)
  # Hand calculation: sum of indicators (ties at half) is 4.5 over 25 pairs.
  expect_equal(res$A, 4.5 / 25, tolerance = 1e-10)
})

test_that("vargha_delaney_A() equals 0.5 for identical samples", {
  set.seed(113)
  x <- stats::rnorm(40)
  res <- vargha_delaney_A(x, x)
  expect_equal(res$A, 0.5, tolerance = 1e-10)
})

test_that("vargha_delaney_A() and 1 - A are symmetric across group order", {
  set.seed(113)
  x <- stats::rnorm(40, mean = 0.6)
  y <- stats::rnorm(40, mean = 0)
  res_xy <- vargha_delaney_A(x, y)
  res_yx <- vargha_delaney_A(y, x)
  expect_equal(res_xy$A, 1 - res_yx$A, tolerance = 1e-10)
  expect_equal(res_xy$se, res_yx$se,   tolerance = 1e-10)
})

test_that("vargha_delaney_A() matches the Mann-Whitney U / (n_1 n_2) computation", {
  set.seed(113)
  x <- stats::rnorm(50, mean = 0.5)
  y <- stats::rnorm(50)
  res <- vargha_delaney_A(x, y)
  # Independent computation: U_x via wilcox.test (R reports U for x as the
  # MWU "statistic" when 'x' is the first argument).
  U_x <- as.numeric(stats::wilcox.test(x, y, exact = FALSE)$statistic)
  expect_equal(res$A, U_x / (length(x) * length(y)), tolerance = 1e-8)
})

test_that("vargha_delaney_A() DeLong variance reproduces a hand calculation on a small sample", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(3, 4, 5, 6, 7)
  res <- vargha_delaney_A(x, y, ci_method = "wald")
  # V_10 placements per X_i:
  #   X=1: 0; X=2: 0; X=3: 0.1; X=4: 0.3; X=5: 0.5
  V_10 <- c(0, 0, 0.1, 0.3, 0.5)
  V_01 <- c(0.5, 0.3, 0.1, 0, 0)
  var_check <- stats::var(V_10) / 5 + stats::var(V_01) / 5
  expect_equal(res$se^2, var_check, tolerance = 1e-10)
})

test_that("vargha_delaney_A() formula interface matches two-vector interface", {
  set.seed(113)
  res_form <- vargha_delaney_A(len ~ supp, data = ToothGrowth)
  x_oj <- ToothGrowth$len[ToothGrowth$supp == "OJ"]
  x_vc <- ToothGrowth$len[ToothGrowth$supp == "VC"]
  res_vec <- vargha_delaney_A(x_oj, x_vc)
  expect_equal(res_form$A,  res_vec$A,  tolerance = 1e-10)
  expect_equal(res_form$se, res_vec$se, tolerance = 1e-10)
})

test_that("vargha_delaney_A() CI lies in [0, 1] and brackets the point estimate", {
  set.seed(113)
  x <- stats::rnorm(40, mean = 0.6)
  y <- stats::rnorm(40, mean = 0)
  res <- vargha_delaney_A(x, y)
  expect_gte(res$lower_limit, 0)
  expect_lte(res$upper_limit, 1)
  expect_lt(res$lower_limit, res$A)
  expect_gt(res$upper_limit, res$A)
})

test_that("vargha_delaney_A() wald CI agrees with normal approximation", {
  set.seed(113)
  x <- stats::rnorm(40, mean = 0.6)
  y <- stats::rnorm(40, mean = 0)
  res <- vargha_delaney_A(x, y, ci_method = "wald")
  z   <- stats::qnorm(0.975)
  expect_equal(res$lower_limit, max(0, res$A - z * res$se),
               tolerance = 1e-10)
  expect_equal(res$upper_limit, min(1, res$A + z * res$se),
               tolerance = 1e-10)
})
