test_that("bessel_errors has the documented shape", {
  data(bessel_errors, package = "DMAR", envir = environment())
  expect_s3_class(bessel_errors, "data.frame")
  expect_equal(dim(bessel_errors), c(9L, 6L))
  expect_equal(names(bessel_errors),
               c("bin", "lower", "upper", "midpoint",
                 "observed", "expected"))
})

test_that("bin edges are consistent and the midpoints are correct", {
  data(bessel_errors, package = "DMAR", envir = environment())
  # Bins are abutting 0.0-0.9 in 0.1-wide intervals.
  expect_equal(bessel_errors$lower, seq(0.0, 0.8, by = 0.1))
  expect_equal(bessel_errors$upper, seq(0.1, 0.9, by = 0.1))
  expect_equal(bessel_errors$midpoint,
               (bessel_errors$lower + bessel_errors$upper) / 2)
  expect_equal(bessel_errors$bin, 1:9)
})

test_that("frequencies match Bessel (1818) / MDK4 Table 1.4", {
  data(bessel_errors, package = "DMAR", envir = environment())
  expect_equal(bessel_errors$observed,
               c(114L, 84L, 53L, 24L, 14L, 6L, 3L, 1L, 1L))
  expect_equal(bessel_errors$expected,
               c(107L, 87L, 57L, 30L, 13L, 5L, 1L, 0L, 0L))
  expect_equal(sum(bessel_errors$observed), 300L)
  expect_equal(sum(bessel_errors$expected), 300L)
})

test_that("frequency-weighted mean using bin midpoints recovers MDK4 commentary", {
  data(bessel_errors, package = "DMAR", envir = environment())
  wmean_obs <- with(bessel_errors,
                    sum(observed * midpoint) / sum(observed))
  wmean_exp <- with(bessel_errors,
                    sum(expected * midpoint) / sum(expected))
  # Both means should sit in the 0.15 to 0.20 seconds-of-arc range
  # given the bulk of the distribution is in the first two bins,
  # and the two should agree to roughly one decimal place.
  expect_true(wmean_obs > 0.10 && wmean_obs < 0.25)
  expect_true(wmean_exp > 0.10 && wmean_exp < 0.25)
  expect_lt(abs(wmean_obs - wmean_exp), 0.02)
})

test_that("the expected counts imply a normal sigma near 0.22, not 0.2", {
  # Anchors the ?bessel_errors @format prose: a least squares fit of the
  # half-normal bin probabilities to Bessel's expected frequencies gives
  # sigma = 0.216, documented as approximately 0.22. The 0.2 the page
  # previously reported misses the first bin by about eight observations
  # (114.9 expected at sigma = 0.2 against the tabled 107).
  data(bessel_errors, package = "DMAR", envir = environment())
  ss <- function(sigma) {
    pr <- 2 * (pnorm(bessel_errors$upper / sigma) -
               pnorm(bessel_errors$lower / sigma))
    sum((bessel_errors$expected - 300 * pr)^2)
  }
  sigma_hat <- optimize(ss, c(0.05, 1))$minimum
  expect_equal(sigma_hat, 0.216, tolerance = 0.005)
  expect_lt(ss(sigma_hat), 1)
  expect_gt(ss(0.2), 50)
})
