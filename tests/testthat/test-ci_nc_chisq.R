test_that("ci_nc_chisq() returns a 2-row data frame with expected columns", {
  result <- ci_nc_chisq(chisq = 30, df = 15, conf_level = 0.95)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2L)
  expect_equal(result$term, c("ci_lower", "ci_upper"))
  expect_named(result, c("term", "value", "prob_less", "prob_greater"))
})

test_that("ci_nc_chisq() ci_lower < ci_upper for significant Chi", {
  result <- ci_nc_chisq(chisq = 30, df = 15, conf_level = 0.95)
  ll <- result$value[result$term == "ci_lower"]
  ul <- result$value[result$term == "ci_upper"]
  expect_lt(ll, ul)
})

test_that("ci_nc_chisq() achieves target tail probabilities to high precision", {
  result <- ci_nc_chisq(chisq = 30, df = 15, conf_level = 0.95)
  alpha <- 0.025
  expect_equal(result$prob_greater[result$term == "ci_lower"], alpha, tolerance = 1e-10)
  expect_equal(result$prob_less[result$term == "ci_upper"],    alpha, tolerance = 1e-10)
})

test_that("ci_nc_chisq() one-sided upper interval: ci_lower is 0", {
  result <- ci_nc_chisq(chisq = 30, df = 15, alpha_lower = 0, alpha_upper = 0.05, conf_level = NULL)
  expect_equal(result$value[result$term == "ci_lower"], 0)
})

test_that("ci_nc_chisq() one-sided lower interval: ci_upper is Inf", {
  result <- ci_nc_chisq(chisq = 30, df = 15, alpha_lower = 0.05, alpha_upper = 0, conf_level = NULL)
  expect_equal(result$value[result$term == "ci_upper"], Inf)
})

test_that("ci_nc_chisq() warns and clamps when chisq is below central critical value", {
  # chisq = 5 on 15 df is so small that BOTH tails degenerate: the lower
  # limit clamps to 0 and the upper limit is undefined (NA); assert both warnings.
  expect_warning(
    expect_warning(
      result <- ci_nc_chisq(chisq = 5, df = 15, conf_level = 0.95),
      "below the alpha_lower critical value"),
    "upper noncentrality limit is undefined"
  )
  expect_equal(result$value[result$term == "ci_lower"], 0)
  expect_true(is.na(result$value[result$term == "ci_upper"]))
})

test_that("ci_nc_chisq() errors on bad inputs", {
  expect_error(ci_nc_chisq(),                              "chisq")
  expect_error(ci_nc_chisq(chisq = -1, df = 15),      "non-negative")
  expect_error(ci_nc_chisq(chisq = 30),               "degrees of freedom")
  expect_error(ci_nc_chisq(chisq = 30, df = 15, conf_level = 0.95, alpha_lower = 0.025), "Specify either")
  expect_error(ci_nc_chisq(chisq = 30, df = 15, conf_level = NULL, alpha_lower = 0.5, alpha_upper = 0.6), "non-negative and sum to less than 1")
})

test_that("ci_nc_chisq(verbose = FALSE) drops probability columns", {
  result <- ci_nc_chisq(chisq = 30, df = 15, conf_level = 0.95, verbose = FALSE)
  expect_named(result, c("term", "value"))
})

test_that("ci_nc_chisq() asymmetric alphas use them in the correct direction", {
  # Chi = 50 with df = 15 is highly significant (central upper-tail ~ 1e-5),
  # so a positive lower limit exists at the stricter alpha_lower = 0.01.
  result <- ci_nc_chisq(chisq = 50, df = 15, alpha_lower = 0.01, alpha_upper = 0.04, conf_level = NULL)
  expect_equal(result$prob_greater[result$term == "ci_lower"], 0.01, tolerance = 1e-10)
  expect_equal(result$prob_less[result$term == "ci_upper"],    0.04, tolerance = 1e-10)
})

# An evaluation of the noncentral chi square distribution function that
# shares nothing with pchisq(ncp =) or the Poisson mixture:
# X = (Z + sqrt(ncp))^2 + a central chi square on df - 1, integrated over Z.
pnc_integral <- function(q, df, ncp) {
  s <- sqrt(ncp)
  f <- function(z) stats::pchisq(pmax(q - (z + s)^2, 0), df - 1) * stats::dnorm(z)
  stats::integrate(f, max(-sqrt(q) - s, -12), sqrt(q) - s, rel.tol = 1e-13,
                   abs.tol = 0, subdivisions = 2000)$value
}

test_that("a statistic in the millions gets limits with the requested tails", {
  # R's pchisq() stops its series after 10^6 terms, short of the Poisson
  # weights of a noncentrality near 3e6: both limits once collapsed onto
  # 2979734, below the statistic, with prob_less 1 and 0.
  x <- ci_nc_chisq(chisq = 3e6, df = 3000, conf_level = 0.90)
  lo <- x$value[1]; hi <- x$value[2]
  expect_lt(lo, 3e6 - 3000)
  expect_gt(hi, 3e6 - 3000)
  expect_equal(x$prob_greater[1], 0.05, tolerance = 1e-8)
  expect_equal(x$prob_less[2], 0.05, tolerance = 1e-8)
  expect_equal(1 - pnc_integral(3e6, 3000, lo), 0.05, tolerance = 1e-7)
  expect_equal(pnc_integral(3e6, 3000, hi), 0.05, tolerance = 1e-7)
  # At 2e6 the lower limit was the one far off (1982913 rather than about
  # 1992353).
  y <- ci_nc_chisq(chisq = 2e6, df = 3000, conf_level = 0.90)
  expect_equal(1 - pnc_integral(2e6, 3000, y$value[1]), 0.05, tolerance = 1e-7)
  expect_equal(pnc_integral(2e6, 3000, y$value[2]), 0.05, tolerance = 1e-7)
})

test_that("the distribution function is R's at and below a noncentrality of 1e5", {
  for (ncp in c(0, 10, 1e3, 1e5)) {
    expect_identical(.pchisq_nc(ncp + 20, 20, ncp), stats::pchisq(ncp + 20, 20, ncp = ncp))
    expect_identical(.pchisq_nc(ncp + 20, 20, ncp, lower.tail = FALSE),
                     stats::pchisq(ncp + 20, 20, ncp = ncp, lower.tail = FALSE))
  }
  # Above it the mixture agrees with the integral where R's is still usable.
  q <- 1e6 + 10 + 1.5 * sqrt(2 * (10 + 2e6))
  expect_equal(.pchisq_nc(q, 10, 1e6), pnc_integral(q, 10, 1e6), tolerance = 1e-12)
  expect_equal(.pchisq_nc(q, 10, 1e6, lower.tail = FALSE),
               1 - pnc_integral(q, 10, 1e6), tolerance = 1e-10)
})
