## from tests/testthat/test-power_fisher_exact.R and R/power_fisher_exact.R
local({
  # DMAR's internal Fisher noncentral hypergeometric density against
  # BiasedUrn's, over the sweep used when the dependency was removed.
  worst <- 0
  for (m1 in c(5, 12, 40)) for (m2 in c(7, 25)) for (s in c(3, 10, 30)) {
    for (psi in c(0.2, 1, 2.7, 8)) {
      if (s > m1 + m2) next
      xs <- max(0, s - m2):min(m1, s)
      a <- DMAR:::.dfnc_hypergeo(xs, m1, m2, s, psi)
      b <- BiasedUrn::dFNCHypergeo(xs, m1, m2, s, psi)
      worst <- max(worst, max(abs(a - b)))
    }
  }
  stopifnot(worst < 1e-10)
  # And the Wallenius contrast: a different distribution, not a variant.
  n_1 <- 8; n_2 <- 6; s <- 7; psi <- 2.5
  x <- max(0, s - n_2):min(n_1, s)
  cond <- DMAR:::.dfnc_hypergeo(x, n_1, n_2, s, psi)
  stopifnot(max(abs(BiasedUrn::dWNCHypergeo(x, n_1, n_2, s, psi) - cond)) > 0.01)
})
