test_that("power_fisher_exact() returns the documented rows", {
  res <- power_fisher_exact(n_1 = 20, n_2 = 20, p_1 = 0.6, p_2 = 0.3)
  expect_setequal(res$term,
                  c("power", "odds_ratio_alt", "p_1", "p_2",
                    "expected_s", "n_1", "n_2", "alpha_level"))
})

test_that("power_fisher_exact() power is in [0, 1]", {
  res <- power_fisher_exact(n_1 = 30, n_2 = 30, p_1 = 0.6, p_2 = 0.3)
  pwr <- res$value[res$term == "power"]
  expect_gte(pwr, 0)
  expect_lte(pwr, 1)
})

test_that("power_fisher_exact() approaches alpha when p_1 = p_2", {
  res <- power_fisher_exact(n_1 = 30, n_2 = 30, p_1 = 0.5, p_2 = 0.5,
                            alpha_level = 0.05)
  pwr <- res$value[res$term == "power"]
  expect_lt(pwr, 0.10)
})

test_that("power_fisher_exact() rejects invalid probabilities", {
  expect_error(power_fisher_exact(20, 20, 1.2, 0.3),  "in \\(0, 1\\)")
  expect_error(power_fisher_exact(20, 20, 0.6, 1.5),  "in \\(0, 1\\)")
})

test_that("the alternative is Fisher's noncentral hypergeometric", {
  # Fisher's noncentral hypergeometric is the conditional distribution
  # of one binomial count given the total of two independent binomials;
  # the internal .dfnc_hypergeo() must reproduce that pmf to machine
  # precision. The same check against BiasedUrn::dFNCHypergeo (BiasedUrn
  # 2.0.12; max deviation 1.2e-13 over a 66-configuration sweep,
  # 2026-08-10) lives in tools/oracle_checks.R, along with the contrast
  # against Wallenius' distribution, which differs by 0.068 here.
  n_1 <- 8; n_2 <- 6; s <- 7; psi <- 2.5
  x <- max(0, s - n_2):min(n_1, s)
  cond <- choose(n_1, x) * choose(n_2, s - x) * psi^x
  cond <- cond / sum(cond)
  expect_equal(DMAR:::.dfnc_hypergeo(x, n_1, n_2, s, psi), cond,
               tolerance = 1e-12)
  # Vectorized input outside the support returns zero density there.
  expect_equal(DMAR:::.dfnc_hypergeo(c(-1, x, n_1 + 1), n_1, n_2, s, psi),
               c(0, cond, 0), tolerance = 1e-12)
})

test_that("power_fisher_exact() equals a direct two-binomial enumeration", {
  # Independent anchor: weighting each conditional rejection probability
  # by the distribution of the margin total must equal the
  # unconditional two-binomial probability of landing in the same
  # conditional rejection regions. The identity holds exactly for
  # Fisher's noncentral hypergeometric and not for Wallenius'.
  n_1 <- 12; n_2 <- 10; p_1 <- 0.7; p_2 <- 0.35; alpha <- 0.05
  res <- power_fisher_exact(n_1, n_2, p_1, p_2, alpha_level = alpha)
  pwr <- res$value[res$term == "power"]
  direct <- 0
  for (x_1 in 0:n_1) for (x_2 in 0:n_2) {
    s <- x_1 + x_2
    x_seq <- max(0, s - n_2):min(n_1, s)
    p0 <- stats::dhyper(x_seq, n_1, n_2, s)
    ord <- order(p0); cum <- cumsum(p0[ord])
    in_R <- logical(length(x_seq)); in_R[ord[cum <= alpha]] <- TRUE
    if (in_R[match(x_1, x_seq)])
      direct <- direct + stats::dbinom(x_1, n_1, p_1) *
                         stats::dbinom(x_2, n_2, p_2)
  }
  expect_equal(pwr, direct, tolerance = 1e-12)
})
