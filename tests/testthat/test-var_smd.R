test_that("var_smd() returns exact and approximate rows", {
  res <- var_smd(true_smd = 0.5, n_1 = 20)
  expect_named(res, c("term", "value"))
  expect_setequal(res$term, c("var_smd_exact", "var_smd_approx"))
  expect_true(all(res$value > 0))
})

test_that("var_smd() approximation converges to exact as n grows", {
  small <- var_smd(0.5, n_1 = 10)
  large <- var_smd(0.5, n_1 = 1000)
  # Drift between exact and approx is bigger at small n:
  rel_small <- abs(small$value[1] - small$value[2]) / small$value[1]
  rel_large <- abs(large$value[1] - large$value[2]) / large$value[1]
  expect_gt(rel_small, rel_large)
})

test_that("var_smd() unbiased = TRUE rescales by J^2", {
  res_d <- var_smd(0.5, n_1 = 10)
  res_g <- var_smd(0.5, n_1 = 10, unbiased = TRUE)
  df <- 18
  J  <- gamma(df / 2) / (sqrt(df / 2) * gamma((df - 1) / 2))
  expect_equal(res_g$value[1] / res_d$value[1], J^2, tolerance = 1e-8)
})

test_that("var_smd() rejects bad inputs", {
  expect_error(var_smd("a", 10), "numeric")
  expect_error(var_smd(0.5, 2),  ">=")
})

test_that("var_smd() exact variance keeps its digits at very large samples", {
  # Reference values from a 60-digit evaluation of
  # {df (1 + lambda^2) / (df - 2) - lambda^2 / J^2} / design. The form used
  # through 1.0.x subtracted two numbers of about lambda^2 and returned a
  # negative variance at some sizes.
  v <- function(...) var_smd(...)$value[1]
  expect_equal(v(0.5, 1e8, 1e8), 2.0625000217968754e-8, tolerance = 1e-13)
  expect_equal(v(0.5, 1e8, 2e8), 1.5416666774652779e-8, tolerance = 1e-13)
  expect_equal(v(0.5, 1e5, 5e4), 3.0833765289532727e-5, tolerance = 1e-13)
  expect_equal(var_smd(0.5, 1e8, 1e8, unbiased = TRUE)$value[1],
               2.0625000063281251e-8, tolerance = 1e-13)
  expect_gt(v(3, 1e8, 2e8), 0)
  # Small samples are as before.
  expect_equal(v(0.5, 5, 5), 0.55987050908953842, tolerance = 1e-13)
  expect_equal(v(0.5, 20, 20), 0.10919487978477938, tolerance = 1e-13)
  expect_equal(v(0.5, 5, 2), 1.2296556820180104, tolerance = 1e-13)
  expect_equal(var_smd(0.5, 20, 10, unbiased = TRUE)$value[1],
               0.15776481868012051, tolerance = 1e-13)
})

test_that("var_smd() help-page ratios agree with a closed-form oracle", {
  # Oracle written independently of R/var_smd.R: the noncentral t variance
  # df (1 + lambda^2) / (df - 2) - lambda^2 c(df)^2 with
  # c(df) = sqrt(df / 2) Gamma((df - 1) / 2) / Gamma(df / 2)
  # (Johnson, Kotz, & Balakrishnan, 1995, Sec. 31.3), divided by
  # n_1 n_2 / (n_1 + n_2), and the Hedges-Olkin form
  # (n_1 + n_2) / (n_1 n_2) + delta^2 / (2 df).
  oracle <- function(delta, n1, n2 = n1) {
    df <- n1 + n2 - 2
    design <- n1 * n2 / (n1 + n2)
    lambda <- delta * sqrt(design)
    c_df <- sqrt(df / 2) * exp(lgamma((df - 1) / 2) - lgamma(df / 2))
    exact <- (df * (1 + lambda^2) / (df - 2) - lambda^2 * c_df^2) / design
    approx <- (n1 + n2) / (n1 * n2) + delta^2 / (2 * df)
    c(exact = exact, approx = approx, J2 = 1 / c_df^2)
  }
  ratio <- function(delta, n) {
    v <- var_smd(true_smd = delta, n_1 = n)$value
    v[2] / v[1]
  }
  # Agreement with the oracle, both signs of the effect.
  for (delta in c(-2, -1, 0, 0.5, 1, 2)) {
    for (n in c(5, 10, 30, 100)) {
      o <- oracle(delta, n)
      v <- var_smd(true_smd = delta, n_1 = n)$value
      expect_equal(v, unname(o[c("exact", "approx")]), tolerance = 1e-10)
      expect_lt(v[2], v[1])
    }
  }
  # The numbers the description quotes.
  expect_equal(ratio(0, 5), 0.75, tolerance = 1e-12)
  expect_equal(round(ratio(0, 30), 2), 0.97)
  expect_equal(round(ratio(1, 5), 2), 0.72)
  expect_equal(round(ratio(1, 30), 2), 0.96)
  expect_equal(round(ratio(-1, 30), 2), 0.96)
  expect_equal(round(unname(oracle(0, 5)["J2"]), 2), 0.81)
  expect_equal(round(unname(oracle(0, 20)["J2"]), 2), 0.96)
  g_over_d <- var_smd(0.5, n_1 = 5, unbiased = TRUE)$value[1] /
    var_smd(0.5, n_1 = 5)$value[1]
  expect_equal(g_over_d, unname(oracle(0.5, 5)["J2"]), tolerance = 1e-10)
  # The ratio falls as |delta| grows.
  r <- vapply(c(0, 0.5, 1, 2, 4), ratio, numeric(1), n = 10)
  expect_true(all(diff(r) < 0))
})
