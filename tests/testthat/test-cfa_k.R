cfa_k_hs_factors <- function() {
  list(
    verbal = c("t6_paragraph_comprehension", "t7_sentence",
               "t9_word_meaning"),
    deduction = c("t20_deduction", "t22_problem_reasoning",
               "t23_series_completion")
  )
}

cfa_k_hs_data <- function() {
  data(holzinger_swineford, package = "DMAR", envir = environment())
  holzinger_swineford
}

cfa_k_example_cov <- function() {
  S <- matrix(
    c(1.384, 1.484, 1.988, 2.429, 3.031,
      1.484, 2.756, 2.874, 3.588, 4.390,
      1.988, 2.874, 4.845, 4.894, 6.080,
      2.429, 3.588, 4.894, 6.951, 7.476,
      3.031, 4.390, 6.080, 7.476, 10.313),
    nrow = 5
  )
  dimnames(S) <- list(paste0("y", 1:5), paste0("y", 1:5))
  S
}

test_that("cfa_k() default output has the documented shape and terms", {
  skip_if_not_installed("lavaan")

  res <- cfa_k(cfa_k_hs_data(), cfa_k_hs_factors())

  expect_s3_class(res, "dmar_cfa_k")
  expect_s3_class(res, "dmar_tbl")
  expect_named(res, c("syntax", "term", "estimate", "se", "z_value",
                      "p_value", "ci_lower", "ci_upper"))
  expect_true(all(c("lambda_verbal_1", "lambda_deduction_3",
                    "psi_verbal_1", "psi_deduction_3",
                    "phi_verbal", "phi_deduction", "phi_verbal_deduction",
                    "loading_sum_verbal", "error_sum_deduction",
                    "omega_verbal", "omega_deduction",
                    "ave_verbal", "H_deduction") %in% res$term))
  expect_true(all(c("chi_square", "df", "p_chi_square", "cfi", "tli",
                    "rmsea", "rmsea_ci_lower", "rmsea_ci_upper",
                    "rmsea_ci_level", "srmr", "AIC", "BIC")
                  %in% res$term))
  expect_equal(res$estimate[res$term == "rmsea_ci_level"], 0.90)
  expect_equal(res$estimate[res$term == "phi_verbal"], 1)

  ci <- res[res$term == "omega_verbal", ]
  expect_true(is.finite(ci$se))
  expect_lt(ci$ci_lower, ci$estimate)
  expect_gt(ci$ci_upper, ci$estimate)

  model <- attr(res, "model")
  expect_named(model, c("verbal", "deduction"))
  expect_match(model[["verbal"]], "^congeneric")
})

test_that("cfa_k() names the implied classical structures correctly", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # five fits over the descriptor grid; the check below runs on CRAN
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()

  ess_tau <- cfa_k(hs, fs, equal_loading = TRUE)
  expect_match(attr(ess_tau, "model")[["verbal"]],
               "^essentially tau-equivalent")
  expect_true("lambda_verbal" %in% ess_tau$term)
  expect_false("lambda_verbal_1" %in% ess_tau$term)

  tau <- cfa_k(hs, fs, equal_loading = TRUE, equal_intercept = TRUE)
  expect_match(attr(tau, "model")[["deduction"]], "^tau-equivalent")
  expect_true("nu_deduction" %in% tau$term)

  par <- cfa_k(hs, fs, equal_loading = TRUE, equal_intercept = TRUE,
               equal_error = TRUE)
  expect_match(attr(par, "model")[["verbal"]], "^parallel")

  ess_par <- cfa_k(hs, fs, equal_loading = TRUE, equal_error = TRUE)
  expect_match(attr(ess_par, "model")[["verbal"]], "^essentially parallel")

  odd <- cfa_k(hs, fs, equal_error = TRUE)
  expect_match(attr(odd, "model")[["verbal"]], "^no classical name")
})

test_that("cfa_k() descriptor flags can differ by factor", {
  skip_if_not_installed("lavaan")

  res <- cfa_k(cfa_k_hs_data(), cfa_k_hs_factors(),
               equal_loading = c(verbal = TRUE, deduction = FALSE))
  expect_match(attr(res, "model")[["verbal"]],
               "^essentially tau-equivalent")
  expect_match(attr(res, "model")[["deduction"]], "^congeneric")
  expect_true("lambda_verbal" %in% res$term)
  expect_true("lambda_deduction_1" %in% res$term)
})

test_that("cfa_k() omega agrees with cfa_1() on the same one-factor model", {
  skip_if_not_installed("lavaan")

  S <- cfa_k_example_cov()
  res_k <- cfa_k(S = S, factors = list(f = paste0("y", 1:5)), N = 300)
  res_1 <- cfa_1(S = unname(S), N = 300)

  expect_equal(res_k$estimate[res_k$term == "omega_f"],
               res_1$estimate[res_1$term == "omega_f1"], tolerance = 1e-5)
  expect_equal(res_k$se[res_k$term == "omega_f"],
               res_1$se[res_1$term == "omega_f1"], tolerance = 1e-5)

  # Agreement with cfa_1() is an internal consistency check: both routes
  # share DMAR's own omega machinery, so it would hold just as well if that
  # machinery were wrong. Anchor the coefficient and its standard error to
  # MBESS, the reference implementation, so a defect in the shared path
  # cannot pass unnoticed. Pinned from MBESS::ci.reliability (MBESS 4.9.3,
  # 2026-08-09); live comparison in tools/oracle_checks.R.
  expect_equal(res_k$estimate[res_k$term == "omega_f"], 0.9621012361618826,
               tolerance = 1e-6)
  expect_equal(res_k$se[res_k$term == "omega_f"], 0.003536944590571952,
               tolerance = 1e-6)
})

test_that("cfa_k() ave agrees with average_variance_extracted()", {
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()

  res <- cfa_k(hs, fs)
  fit <- cfa_k(hs, fs, output = "fit")
  ave <- average_variance_extracted(fit)

  expect_equal(res$estimate[res$term == "ave_verbal"],
               ave$ave[ave$factor == "verbal"], tolerance = 1e-5)
  expect_equal(res$estimate[res$term == "ave_deduction"],
               ave$ave[ave$factor == "deduction"], tolerance = 1e-5)
})

test_that("cfa_k(output='measurement') reports properties, phi, and htmt", {
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()

  res <- cfa_k(hs, fs, output = "measurement")
  expect_true(all(c("omega_verbal", "omega_deduction", "ave_verbal",
                    "ave_deduction", "H_verbal", "H_deduction",
                    "phi_verbal_deduction", "htmt_verbal_deduction")
                  %in% res$term))
  expect_false(any(grepl("^lambda_", res$term)))

  ht <- htmt(hs, blocks = fs)
  expect_equal(res$estimate[res$term == "htmt_verbal_deduction"],
               ht$htmt[1], tolerance = 1e-8)

  phi <- res[res$term == "phi_verbal_deduction", ]
  expect_true(is.finite(phi$ci_lower) && is.finite(phi$ci_upper))
})

test_that("cfa_k() supports a covariance matrix with means for intercepts", {
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()
  items <- unlist(fs)
  hs_c <- hs[stats::complete.cases(hs[, items]), items]
  S <- stats::cov(hs_c)
  M <- colMeans(hs_c)

  res <- cfa_k(S = S, factors = fs, N = nrow(hs_c), M = M, equal_loading = TRUE,
               equal_intercept = TRUE)
  expect_match(attr(res, "model")[["verbal"]], "^tau-equivalent")
  expect_true("nu_verbal" %in% res$term)
})

test_that("cfa_k() nested descriptor fits differ by the right df", {
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()

  fit_free <- cfa_k(hs, fs, output = "fit")
  fit_equal <- cfa_k(hs, fs, equal_loading = TRUE, output = "fit")
  df_free <- unname(lavaan::fitMeasures(fit_free, "df"))
  df_equal <- unname(lavaan::fitMeasures(fit_equal, "df"))
  # Two constraints per three-item factor, two factors.
  expect_equal(df_equal - df_free, 4)
})

test_that("cfa_k() validates its inputs", {
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()

  expect_error(cfa_k(hs, list(c("t10_addition", "t11_code"))),
               "named list")
  expect_error(cfa_k(hs, list(a = c("t10_addition"),
                              b = c("t11_code", "t12_counting_groups_of_dots"))),
               "two or more items")
  expect_error(cfa_k(hs, list(a = c("t10_addition", "t11_code"),
                              b = c("t11_code", "t12_counting_groups_of_dots"))),
               "exactly one factor")
  expect_error(cfa_k(hs, list(f = c("t10_addition", "t11_code"))),
               "three or more items")
  expect_error(cfa_k(hs, fs, estimator = "PML"), "'estimator'")
  expect_error(cfa_k(hs, fs, equal_loading = c(verbal = TRUE)),
               "named with exactly the factor names")
  expect_error(cfa_k(S = cfa_k_example_cov(),
                     factors = list(f = paste0("y", 1:5)), N = 300,
                     equal_intercept = TRUE),
               "requires the item means")
  expect_error(cfa_k(hs, fs, equal_intercept = TRUE,
                     meanstructure = FALSE),
               "meanstructure")
})

test_that("cfa_k() tidy() and glance() follow the broom conventions", {
  skip_if_not_installed("lavaan")

  res <- cfa_k(cfa_k_hs_data(), cfa_k_hs_factors())

  td <- generics::tidy(res)
  expect_named(td, c("term", "estimate", "se", "statistic",
                     "p_value", "ci_lower", "ci_upper"))
  expect_false(any(td$term %in% c("cfi", "rmsea", "AIC")))

  gl <- generics::glance(res)
  expect_equal(nrow(gl), 1L)
  expect_true(is.finite(gl$cfi))
  expect_true(is.finite(gl$srmr))
  expect_equal(gl$df, unname(
    lavaan::fitMeasures(cfa_k(cfa_k_hs_data(), cfa_k_hs_factors(),
                              output = "fit"), "df")))
})


# ---------------------------------------------------------------------------
# Ordered-categorical items and the Green and Yang per-factor omega
# ---------------------------------------------------------------------------

set.seed(113)
N_ord <- 400
lam_ord <- c(0.5, 0.6, 0.7, 0.75, 0.55, 0.65, 0.7, 0.6)
eta1_o <- rnorm(N_ord)
eta2_o <- 0.5 * eta1_o + sqrt(0.75) * rnorm(N_ord)
lat_o <- cbind(
  sweep(matrix(rep(eta1_o, 4), N_ord, 4), 2, lam_ord[1:4], `*`),
  sweep(matrix(rep(eta2_o, 4), N_ord, 4), 2, lam_ord[5:8], `*`)) +
  matrix(rnorm(N_ord * 8), N_ord, 8) %*% diag(sqrt(1 - lam_ord^2))
colnames(lat_o) <- paste0("V", 1:8)
items_ord <- apply(lat_o, 2, function(x)
  as.integer(cut(x, breaks = c(-Inf, -1, -0.2, 0.6, Inf))))
colnames(items_ord) <- paste0("o", 1:8)
factors_ord <- list(F1 = paste0("o", 1:4), F2 = paste0("o", 5:8))

test_that("ordered single-factor cfa_k matches reliability_omega_categorical", {
  msgs <- testthat::capture_messages(
    r <- cfa_k(items_ord[, 1:4], list(F1 = paste0("o", 1:4)),
               ordered = TRUE, output = "measurement"))
  expect_true(any(grepl("WLSMV", msgs)))
  expect_true(any(grepl("categorical sum score", msgs)))
  o_k <- as.data.frame(r)$estimate[as.data.frame(r)$term == "omega_F1"]
  ref <- suppressMessages(
    reliability_omega_categorical(data = items_ord[, 1:4],
                                  ci_method = "none"))
  expect_equal(o_k, ref$value[ref$term == "estimate"], tolerance = 1e-6)
})

test_that("ordered multi-factor cfa_k reports GY omegas with NA intervals", {
  skip_on_cran()  # WLSMV on eight ordered items; the ordered anchor above runs on CRAN
  r <- suppressMessages(
    cfa_k(items_ord, factors_ord, ordered = TRUE, output = "measurement"))
  d <- as.data.frame(r)
  for (f in c("F1", "F2")) {
    row_f <- d$term == paste0("omega_", f)
    expect_true(is.finite(d$estimate[row_f]))
    expect_gt(d$estimate[row_f], 0)
    expect_lt(d$estimate[row_f], 1)
    expect_true(is.na(d$se[row_f]))
    expect_true(is.na(d$ci_lower[row_f]))
  }
  # phi keeps its delta method interval.
  phi_row <- d$term == "phi_F1_F2"
  expect_false(is.na(d$ci_lower[phi_row]))
  expect_identical(unname(attr(r, "omega_metric")),
                   c("categorical_sum_score", "categorical_sum_score"))
  # Per-factor values sit near the single-block standalone estimates.
  ref1 <- suppressMessages(
    reliability_omega_categorical(data = items_ord[, 1:4],
                                  ci_method = "none"))
  expect_equal(d$estimate[d$term == "omega_F1"],
               ref1$value[ref1$term == "estimate"], tolerance = 0.05)
})

test_that("mixed model keeps the delta interval for the continuous factor", {
  items_mixed <- cbind(as.data.frame(items_ord[, 1:4]),
                       as.data.frame(lat_o[, 5:8]))
  names(items_mixed) <- c(paste0("o", 1:4), paste0("c", 1:4))
  fac_mixed <- list(F1 = paste0("o", 1:4), F2 = paste0("c", 1:4))
  r <- suppressMessages(
    cfa_k(items_mixed, fac_mixed, ordered = paste0("o", 1:4),
          output = "measurement"))
  d <- as.data.frame(r)
  expect_true(is.na(d$se[d$term == "omega_F1"]))
  expect_false(is.na(d$se[d$term == "omega_F2"]))
  expect_identical(unname(attr(r, "omega_metric")),
                   c("categorical_sum_score", "linear"))
})

test_that("ordered input validation fails loudly", {
  expect_error(
    cfa_k(items_ord, factors_ord, ordered = paste0("o", 1:3)),
    regexp = "mixes ordered and continuous")
  expect_error(
    cfa_k(items_ord, factors_ord, ordered = c("o1", "nope")),
    regexp = "not in 'factors'")
  expect_error(
    cfa_k(S = cov(items_ord), factors = factors_ord, N = N_ord, ordered = TRUE),
    regexp = "require raw data")
  expect_error(
    suppressMessages(cfa_k(items_ord, factors_ord, ordered = TRUE,
                           equal_intercept = TRUE)),
    regexp = "Thresholds replace intercepts")
  expect_error(
    suppressMessages(cfa_k(items_ord, factors_ord, ordered = TRUE,
                           missing = "fiml")),
    regexp = "pairwise")
})

test_that("continuous models record the linear omega metric", {
  r <- cfa_k(S = cov(lat_o), factors = list(F1 = paste0("V", 1:4)), N = N_ord,
             output = "measurement")
  expect_identical(unname(attr(r, "omega_metric")), "linear")
})

test_that("FIML estimates the mean structure with free intercepts", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two FIML fits on incomplete data; cfa_1()'s FIML tests run on CRAN
  # With missing = "ml", the bare lavaan() interface forces the mean
  # structure on but leaves int.ov.free = FALSE, so every intercept was
  # fixed at 0 and the loadings absorbed the item means: on data shifted
  # away from zero the reported omega was 0.995 against a true 0.802,
  # exactly the defect cfa_1() carried until its own FIML intercept fix.
  # cfa_k() now writes free nu lines under FIML, and the two functions
  # must agree exactly on the same data.
  set.seed(113)
  J <- 5; n <- 300
  lam <- seq(.5, .8, length.out = J)
  eta <- rnorm(n)
  X <- outer(eta, lam) + matrix(rnorm(n * J), n, J) %*% diag(sqrt(1 - lam^2))
  colnames(X) <- paste0("y", seq_len(J))
  set.seed(7)
  X[sample(length(X), 150)] <- NA
  X_shifted <- as.data.frame(X + 5)   # nonzero means expose the defect

  r1 <- cfa_1(data = X_shifted, missing = "ml")
  rk <- cfa_k(X_shifted, list(f = paste0("y", seq_len(J))), missing = "ml")
  expect_equal(rk$estimate[rk$term == "omega_f"],
               r1$estimate[r1$term == "omega_f1"])
  expect_equal(rk$se[rk$term == "omega_f"],
               r1$se[r1$term == "omega_f1"])

  # The impossible combination is refused rather than silently mis-fit.
  expect_error(
    cfa_k(X_shifted, list(f = paste0("y", seq_len(J))), missing = "ml",
          meanstructure = FALSE),
    "cannot be combined with", fixed = TRUE
  )
})

test_that("cfa_k() signals a Heywood case with the classed condition cfa_1() uses", {
  skip_if_not_installed("lavaan")
  # For a just-identified one factor model of three indicators the ML
  # solution is exact: lambda_1^2 = r12 r13 / r23. With r12 = r13 = 0.8
  # and r23 = 0.5 that is 1.28, so the first error variance is -0.28, a
  # Heywood case by construction rather than by a lucky draw. The matrix
  # is positive definite (leading minors 1, 0.36, 0.11).
  S <- matrix(c(1, 0.8, 0.8,
                0.8, 1, 0.5,
                0.8, 0.5, 1), 3, 3,
              dimnames = list(paste0("y", 1:3), paste0("y", 1:3)))
  classes <- character()
  withCallingHandlers(
    res <- cfa_k(S = S, N = 200, factors = list(f = paste0("y", 1:3))),
    warning = function(w) {
      classes <<- c(classes, class(w)[1L])
      invokeRestart("muffleWarning")
    }
  )
  expect_true("dmar_heywood_warning" %in% classes)
  expect_s3_class(res, "dmar_cfa_k")
})
