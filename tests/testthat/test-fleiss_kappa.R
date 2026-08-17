test_that("fleiss_kappa() returns a 1-row data frame with the expected columns", {
  ratings <- matrix(c(2, 1, 0, 1, 2, 0, 0, 3, 0, 1, 0, 2), nrow = 4, byrow = TRUE)
  result <- fleiss_kappa(ratings)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_true(all(c("kappa", "se", "lower_limit", "upper_limit", "z_value", "p_value",
                    "n_subjects", "n_raters", "n_categories") %in% names(result)))
})

test_that("fleiss_kappa() lower_limit <= kappa <= upper_limit", {
  ratings <- matrix(c(2, 1, 0, 1, 2, 0, 0, 3, 0, 1, 0, 2), nrow = 4, byrow = TRUE)
  result <- fleiss_kappa(ratings)
  expect_lte(result$lower_limit, result$kappa)
  expect_gte(result$upper_limit, result$kappa)
})

test_that("fleiss_kappa() perfect agreement across varied subjects gives kappa = 1", {
  # 6 subjects with 3 different true categories; 3 raters all agree on each.
  # Using only one category for every subject would make P_e = 1 and leave
  # kappa undefined (0/0); spreading across categories keeps P_e < 1.
  ratings <- rbind(
    c(3, 0, 0), c(3, 0, 0),
    c(0, 3, 0), c(0, 3, 0),
    c(0, 0, 3), c(0, 0, 3)
  )
  result <- fleiss_kappa(ratings)
  expect_equal(result$kappa, 1, tolerance = 1e-10)
})

# Fleiss (1971) Table 1: 30 subjects, 6 raters, 5 diagnostic categories.
# The published point estimate is kappa = 0.430; the corrected null
# variance of Fleiss, Nee, and Landis (1979) gives z = 17.65 on these data.
fleiss_1971_table1 <- matrix(c(
  0, 0, 0, 6, 0,
  0, 3, 0, 0, 3,
  0, 1, 4, 0, 1,
  0, 0, 0, 0, 6,
  0, 3, 0, 3, 0,
  2, 0, 4, 0, 0,
  0, 0, 4, 0, 2,
  2, 0, 3, 1, 0,
  2, 0, 0, 4, 0,
  0, 0, 0, 0, 6,
  1, 0, 0, 5, 0,
  1, 1, 0, 4, 0,
  0, 3, 3, 0, 0,
  1, 0, 0, 5, 0,
  0, 2, 0, 3, 1,
  0, 0, 5, 0, 1,
  3, 0, 0, 1, 2,
  5, 1, 0, 0, 0,
  0, 2, 0, 4, 0,
  1, 0, 2, 0, 3,
  0, 0, 0, 0, 6,
  0, 1, 0, 5, 0,
  0, 2, 0, 1, 3,
  2, 0, 0, 4, 0,
  1, 0, 0, 4, 1,
  0, 5, 0, 1, 0,
  4, 0, 0, 0, 2,
  0, 2, 0, 4, 0,
  1, 0, 5, 0, 0,
  0, 0, 0, 0, 6
), nrow = 30, byrow = TRUE)

test_that("fleiss_kappa() z on the Fleiss (1971) Table 1 example uses the corrected null variance", {
  result <- fleiss_kappa(fleiss_1971_table1)
  # Published point estimate (Fleiss, 1971): kappa = 0.430.
  expect_equal(result$kappa, 0.430, tolerance = 0.002)
  # The Fleiss, Nee, and Landis (1979) null variance gives z = 17.65;
  # the known-marginals variance the function once used gave 15.64.
  expect_equal(round(result$z_value, 2), 17.65)
  # Independent computation of the null variance in the p_j q_j form
  # printed by Fleiss, Nee, and Landis (1979).
  N <- nrow(fleiss_1971_table1)
  m <- sum(fleiss_1971_table1[1L, ])
  p_j <- colSums(fleiss_1971_table1) / (N * m)
  q_j <- 1 - p_j
  var_0 <- 2 * (sum(p_j * q_j)^2 - sum(p_j * q_j * (q_j - p_j))) /
           (N * m * (m - 1) * sum(p_j * q_j)^2)
  expect_equal(result$z_value, result$kappa / sqrt(var_0),
               tolerance = 1e-12)
  expect_equal(result$p_value, 2 * stats::pnorm(-abs(result$z_value)),
               tolerance = 1e-12)
})

test_that("fleiss_kappa() agrees with irr::kappam.fleiss on the Fleiss (1971) example", {
  # Pinned from irr::kappam.fleiss (irr 0.85, 2026-08-09); live comparison in
  # tools/oracle_checks.R.
  result <- fleiss_kappa(fleiss_1971_table1)
  expect_equal(result$kappa, 0.4302445200601406, tolerance = 1e-12)
  expect_equal(result$z_value, 17.65183058299136, tolerance = 1e-12)
})

test_that("fleiss_kappa() interval and se are untouched by the null variance correction", {
  # The se, lower_limit, and upper_limit come from the Gwet (2008)
  # linearization, not the null variance; these pins were computed on the
  # Table 1 example before the null variance fix and must not move.
  result <- fleiss_kappa(fleiss_1971_table1)
  expect_equal(result$se,          0.0541989355153328, tolerance = 1e-12)
  expect_equal(result$lower_limit, 0.3240165584496800, tolerance = 1e-12)
  expect_equal(result$upper_limit, 0.5364724816706020, tolerance = 1e-12)
})

test_that("fleiss_kappa() null variance matches the empirical null at skewed marginals", {
  # Monte Carlo confirmation of the corrected null variance at skewed
  # marginals (N = 1000, m = 5, p = (.85, .10, .05), the audit grid where
  # the known-marginals variance was 3.8 times too large in SD terms).
  # skip_on_cran: the fast anchor that stays on CRAN is the Fleiss (1971)
  # Table 1 pin z = 17.65 above.
  skip_on_cran()
  set.seed(113)
  N <- 1000; m <- 5; p <- c(0.85, 0.10, 0.05)
  G <- 2000
  sims <- vapply(seq_len(G), function(g) {
    ratings <- t(stats::rmultinom(N, m, p))
    result <- fleiss_kappa(ratings)
    c(result$kappa, result$z_value)
  }, numeric(2))
  # Under H_0 the z statistic is asymptotically standard normal. The
  # known-marginals variance gave sd(z) near 0.26 on this grid; the
  # corrected variance gives sd(z) near 1.
  expect_gt(stats::sd(sims[2L, ]), 0.9)
  expect_lt(stats::sd(sims[2L, ]), 1.1)
  # And the closed-form null SD evaluated at the true marginals matches
  # the empirical SD of kappa-hat (about .0080 here, versus the .0303 of
  # the known-marginals formula).
  P_e <- sum(p^2)
  sd_0 <- sqrt(2 * (P_e + P_e^2 - 2 * sum(p^3)) /
               (N * m * (m - 1) * (1 - P_e)^2))
  expect_gt(stats::sd(sims[1L, ]) / sd_0, 0.9)
  expect_lt(stats::sd(sims[1L, ]) / sd_0, 1.1)
})

test_that("bootstrap intervals: reproducible, ordered, near Wald at moderate N", {
  set.seed(113)
  N <- 120; m <- 6; k <- 4
  # Subjects with real agreement structure: a modal category per subject.
  modal <- sample(seq_len(k), N, replace = TRUE)
  ratings <- t(vapply(modal, function(mo) {
    pr <- rep(0.15, k); pr[mo] <- 1 - 0.15 * (k - 1)
    as.vector(stats::rmultinom(1, m, prob = pr))
  }, numeric(k)))

  w  <- fleiss_kappa(ratings)
  p1 <- fleiss_kappa(ratings, ci_method = "percentile", B = 2000, seed = 113)
  p2 <- fleiss_kappa(ratings, ci_method = "percentile", B = 2000, seed = 113)
  b1 <- fleiss_kappa(ratings, ci_method = "bca", B = 2000, seed = 113)

  expect_identical(as.data.frame(p1), as.data.frame(p2))
  expect_equal(p1$kappa, w$kappa)
  expect_equal(p1$se, w$se)
  for (r in list(p1, b1)) {
    expect_lt(r$lower_limit, r$kappa)
    expect_gt(r$upper_limit, r$kappa)
  }
  expect_equal(p1$lower_limit, w$lower_limit, tolerance = 0.06)
  expect_equal(p1$upper_limit, w$upper_limit, tolerance = 0.06)
  expect_identical(attr(p1, "B_used"), 2000L)
})
