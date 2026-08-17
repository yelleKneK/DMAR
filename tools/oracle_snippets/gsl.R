## from tests/testthat/test-expected_r.R
local({
  # The expected_r() help page prints the Olkin-Pratt (Hotelling, 1953)
  # expectation; this evaluates that expression with gsl's independent
  # hypergeometric and checks the package agrees. Small n with large rho on
  # purpose: that is where a misstated parameter shows.
  documented <- function(rho, n) {
    rho * gsl::hyperg_2F1(0.5, 0.5, (n + 1) / 2, rho^2) *
      exp(2 * lgamma(n / 2) - lgamma((n - 1) / 2) - lgamma((n + 1) / 2))
  }
  grid <- list(c(0.9, 5), c(0.8, 6), c(0.5, 8), c(0.3, 10), c(0.95, 7))
  for (g in grid) {
    dmar   <- DMAR::expected_r(rho = g[1], n = g[2])$expected_r
    oracle <- documented(g[1], g[2])
    stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-10)))
  }
})

## from tests/testthat/test-hyperg_2F1.R
local({
  # DMAR's base-R 2F1 (.hyperg_2F1) against gsl over 60 random points away
  # from the x -> 1 boundary, where gsl is itself accurate.
  set.seed(113)
  for (i in 1:60) {
    a  <- sample(1:2, 1)
    N  <- sample(8:1000, 1)
    x  <- runif(1, 0, 0.9)
    cc <- (N + 1) / 2
    dmar   <- DMAR:::.hyperg_2F1(a, a, cc, x)
    oracle <- gsl::hyperg_2F1(a, a, cc, x)
    stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-9)))
  }
})
