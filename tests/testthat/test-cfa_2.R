cfa_2_hs_factors <- function() {
  list(
    factor_1 = c("t6_paragraph_comprehension", "t7_sentence",
                 "t9_word_meaning"),
    factor_2 = c("t20_deduction", "t22_problem_reasoning",
                 "t23_series_completion")
  )
}

test_that("cfa_2() is exactly the two factor cfa_k() call it wraps", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford)
  fs <- cfa_2_hs_factors()

  direct <- cfa_k(holzinger_swineford,
                  factors = list(f1 = fs$factor_1, f2 = fs$factor_2))
  wrapped <- cfa_2(holzinger_swineford,
                   factor_1 = fs$factor_1, factor_2 = fs$factor_2)
  expect_identical(as.data.frame(wrapped), as.data.frame(direct))
  expect_s3_class(wrapped, "dmar_cfa_k")

  expect_true(all(c("lambda_f1_1", "lambda_f2_1",
                    "phi_f1_f2", "omega_f1", "omega_f2") %in% wrapped$term))
})

test_that("cfa_2() passes through the two factor options", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford)
  fs <- cfa_2_hs_factors()

  # Uncorrelated factors: the factor covariance is fixed at zero and
  # its row is omitted from the table.
  orth <- cfa_2(holzinger_swineford,
                factor_1 = fs$factor_1, factor_2 = fs$factor_2,
                correlated_factors = FALSE)
  expect_false("phi_f1_f2" %in% orth$term)
  corr <- cfa_2(holzinger_swineford,
                factor_1 = fs$factor_1, factor_2 = fs$factor_2)
  expect_true("phi_f1_f2" %in% corr$term)

  # Per-factor constraints use the wrapper's factor names f1 and f2.
  mixed <- cfa_2(holzinger_swineford,
                 factor_1 = fs$factor_1, factor_2 = fs$factor_2,
                 equal_loading = c(f1 = TRUE, f2 = FALSE))
  # Constrained loadings collapse to the shared label lambda_f1 (no item
  # index); the free factor keeps indexed lambda_f2_1..3 terms.
  lam_f1 <- mixed$estimate[mixed$term == "lambda_f1"]
  lam_f2 <- mixed$estimate[grepl("^lambda_f2_", mixed$term)]
  expect_length(lam_f1, 3)
  expect_length(lam_f2, 3)
  expect_true(max(abs(diff(lam_f1))) < 1e-10)
  expect_true(max(abs(diff(lam_f2))) > 1e-4)

  # The measurement output reports omega, ave, and H per factor plus htmt.
  meas <- cfa_2(holzinger_swineford,
                factor_1 = fs$factor_1, factor_2 = fs$factor_2,
                output = "measurement")
  expect_true(all(c("omega_f1", "omega_f2", "ave_f1", "ave_f2",
                    "H_f1", "H_f2") %in% meas$term))
})

test_that("cfa_2() validates its item vectors", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford)
  fs <- cfa_2_hs_factors()

  expect_error(cfa_2(holzinger_swineford,
                     factor_1 = fs$factor_1[1], factor_2 = fs$factor_2),
               "two or more items")
  expect_error(cfa_2(holzinger_swineford,
                     factor_1 = fs$factor_1,
                     factor_2 = c(fs$factor_1[1], fs$factor_2[1])),
               "exactly one factor")

  # From a covariance matrix with N, as a paper reports it.
  items <- unlist(fs, use.names = FALSE)
  hs_c <- holzinger_swineford[stats::complete.cases(
    holzinger_swineford[, items]), items]
  from_S <- cfa_2(S = cov(hs_c), N = nrow(hs_c),
                  factor_1 = fs$factor_1, factor_2 = fs$factor_2)
  from_raw <- cfa_2(hs_c, factor_1 = fs$factor_1, factor_2 = fs$factor_2)
  expect_equal(from_S$estimate[grepl("^lambda_", from_S$term)],
               from_raw$estimate[grepl("^lambda_", from_raw$term)],
               tolerance = 1e-6)
})
