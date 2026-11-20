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

test_that("fa_confirmatory() default output has the documented shape and terms", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  res <- fa_confirmatory(cfa_k_hs_data(), cfa_k_hs_factors())

  expect_s3_class(res, "dmar_fa_confirmatory")
  expect_s3_class(res, "dmar_tbl")
  # The 1.0.0 class of cfa_k() trails, for code that tests it, in both
  # table outputs and after the 1.1.0 class.
  expect_identical(class(res), c("dmar_fa_confirmatory", "dmar_cfa_k",
                                 "dmar_tbl", "data.frame"))
  expect_s3_class(fa_confirmatory(cfa_k_hs_data(), cfa_k_hs_factors(),
                                  output = "measurement"), "dmar_cfa_k")
  expect_named(res, c("syntax", "term", "estimate", "se", "test_statistic",
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
  # The RMSEA interval is at rmsea_conf_level, 0.90 by default, apart from
  # the 0.95 of the parameter intervals.
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

test_that("fa_confirmatory() names the implied classical structures correctly", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # five fits over the descriptor grid; the check below runs on CRAN
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()

  ess_tau <- fa_confirmatory(hs, fs, equal_loading = TRUE)
  expect_match(attr(ess_tau, "model")[["verbal"]],
               "^essentially tau-equivalent")
  expect_true("lambda_verbal" %in% ess_tau$term)
  expect_false("lambda_verbal_1" %in% ess_tau$term)

  tau <- fa_confirmatory(hs, fs, equal_loading = TRUE, equal_intercept = TRUE)
  expect_match(attr(tau, "model")[["deduction"]], "^tau-equivalent")
  expect_true("nu_deduction" %in% tau$term)

  par <- fa_confirmatory(hs, fs, equal_loading = TRUE, equal_intercept = TRUE,
               equal_error = TRUE)
  expect_match(attr(par, "model")[["verbal"]], "^parallel")

  ess_par <- fa_confirmatory(hs, fs, equal_loading = TRUE, equal_error = TRUE)
  expect_match(attr(ess_par, "model")[["verbal"]], "^essentially parallel")

  odd <- fa_confirmatory(hs, fs, equal_error = TRUE)
  expect_match(attr(odd, "model")[["verbal"]], "^no classical name")
})

test_that("fa_confirmatory() descriptor flags can differ by factor", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")

  res <- fa_confirmatory(cfa_k_hs_data(), cfa_k_hs_factors(),
               equal_loading = c(verbal = TRUE, deduction = FALSE))
  expect_match(attr(res, "model")[["verbal"]],
               "^essentially tau-equivalent")
  expect_match(attr(res, "model")[["deduction"]], "^congeneric")
  expect_true("lambda_verbal" %in% res$term)
  expect_true("lambda_deduction_1" %in% res$term)
})

test_that("fa_confirmatory() omega agrees with cfa_1() on the same one-factor model", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  S <- cfa_k_example_cov()
  res_k <- fa_confirmatory(S = S, factors = list(f = paste0("y", 1:5)), N = 300)
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

test_that("fa_confirmatory() ave agrees with average_variance_extracted()", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()

  res <- fa_confirmatory(hs, fs)
  fit <- fa_confirmatory(hs, fs, output = "fit")
  ave <- average_variance_extracted(fit)

  expect_equal(res$estimate[res$term == "ave_verbal"],
               ave$ave[ave$factor == "verbal"], tolerance = 1e-5)
  expect_equal(res$estimate[res$term == "ave_deduction"],
               ave$ave[ave$factor == "deduction"], tolerance = 1e-5)
})

test_that("fa_confirmatory(output='measurement') reports properties, phi, and htmt", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()

  res <- fa_confirmatory(hs, fs, output = "measurement")
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

test_that("fa_confirmatory() supports a covariance matrix with means for intercepts", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()
  items <- unlist(fs)
  hs_c <- hs[stats::complete.cases(hs[, items]), items]
  S <- stats::cov(hs_c)
  M <- colMeans(hs_c)

  res <- fa_confirmatory(S = S, factors = fs, N = nrow(hs_c), M = M, equal_loading = TRUE,
               equal_intercept = TRUE)
  expect_match(attr(res, "model")[["verbal"]], "^tau-equivalent")
  expect_true("nu_verbal" %in% res$term)
})

test_that("fa_confirmatory() nested descriptor fits differ by the right df", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()

  fit_free <- fa_confirmatory(hs, fs, output = "fit")
  fit_equal <- fa_confirmatory(hs, fs, equal_loading = TRUE, output = "fit")
  df_free <- unname(lavaan::fitMeasures(fit_free, "df"))
  df_equal <- unname(lavaan::fitMeasures(fit_equal, "df"))
  # Two constraints per three-item factor, two factors.
  expect_equal(df_equal - df_free, 4)
})

test_that("fa_confirmatory() validates its inputs", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()

  expect_error(fa_confirmatory(hs, list(c("t10_addition", "t11_code"))),
               "named list")
  expect_error(fa_confirmatory(hs, list(a = c("t10_addition"),
                              b = c("t11_code", "t12_counting_groups_of_dots"))),
               "two or more items")
  expect_error(fa_confirmatory(hs, list(a = c("t10_addition", "t11_code"),
                              b = c("t11_code", "t12_counting_groups_of_dots"))),
               "exactly one factor")
  expect_error(fa_confirmatory(hs, list(f = c("t10_addition", "t11_code"))),
               "three or more items")
  expect_error(fa_confirmatory(hs, fs, estimator = "PML"), "'estimator'")
  expect_error(fa_confirmatory(hs, fs, equal_loading = c(verbal = TRUE)),
               "named with exactly the factor names")
  expect_error(fa_confirmatory(S = cfa_k_example_cov(),
                     factors = list(f = paste0("y", 1:5)), N = 300,
                     equal_intercept = TRUE),
               "requires the item means")
  expect_error(fa_confirmatory(hs, fs, equal_intercept = TRUE,
                     meanstructure = FALSE),
               "meanstructure")
})

test_that("an item the input lacks is refused with the class of the factor analysis items", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  fs <- cfa_k_hs_factors()
  misnamed <- fs
  misnamed$verbal[2] <- "t7_sentences"
  items_refusal <- function(expr, message) {
    e <- tryCatch(expr, dmar_items_invalid = function(e) e)
    expect_s3_class(e, c("dmar_items_invalid", "dmar_refusal", "error"))
    expect_identical(e$code, "items_invalid")
    expect_identical(e$reason, "unknown")
    expect_identical(e$items, "t7_sentences")
    expect_identical(conditionMessage(e), message)
  }
  # The messages are those DMAR 1.0.x gave without a class.
  items_refusal(fa_confirmatory(hs, misnamed),
                "Items not found in the data: t7_sentences.")
  S <- cov(hs[unlist(fs)])
  items_refusal(fa_confirmatory(S = S, factors = misnamed, N = nrow(hs)),
                "Items not found in the covariance matrix: t7_sentences.")
  # The wrappers reach the same check.
  items_refusal(cfa_1(hs, items = misnamed$verbal),
                "Items not found in the data: t7_sentences.")
  items_refusal(cfa_2(hs, misnamed$verbal, fs$deduction),
                "Items not found in the data: t7_sentences.")
  # A name that is NA is refused for its form, as the items of
  # fa_exploratory() are, with the message it had before the class.
  na_item <- fs
  na_item$verbal[2] <- NA
  e <- tryCatch(fa_confirmatory(hs, na_item), dmar_items_invalid = function(e) e)
  expect_s3_class(e, c("dmar_items_invalid", "dmar_refusal", "error"))
  expect_identical(e$reason, "form")
  expect_identical(e$items, character(0))
  expect_identical(conditionMessage(e), "Items not found in the data: NA.")
  # A factor with too few items is a question of the model rather than of
  # the names, and is not this class.
  e <- tryCatch(fa_confirmatory(hs, list(a = "t10_addition", b = fs$verbal)),
                error = function(e) e)
  expect_false(inherits(e, "dmar_items_invalid"))
})

test_that("fa_confirmatory() tidy() and glance() follow the broom conventions", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")

  res <- fa_confirmatory(cfa_k_hs_data(), cfa_k_hs_factors())

  td <- generics::tidy(res)
  expect_named(td, c("term", "estimate", "se", "statistic",
                     "p_value", "ci_lower", "ci_upper", "conf_level",
                     "wald_df", "h1_diagnostic", "wald_calibration_note"))
  expect_identical(unique(td$conf_level), attr(res, "conf_level"))
  expect_identical(unique(td$wald_df), Inf)
  expect_false(any(td$term %in% c("cfi", "rmsea", "AIC")))

  gl <- generics::glance(res)
  expect_equal(nrow(gl), 1L)
  expect_true(is.finite(gl$cfi))
  expect_true(is.finite(gl$srmr))
  expect_equal(gl$df, unname(
    lavaan::fitMeasures(fa_confirmatory(cfa_k_hs_data(), cfa_k_hs_factors(),
                              output = "fit"), "df")))
  expect_identical(gl[c("estimator", "missing", "test", "fit_indices")],
                   data.frame(estimator = "ML", missing = "fiml",
                              test = "standard", fit_indices = "robust"))
  expect_identical(gl$rmsea_ci_level, 0.90)
})

test_that("tidy() and glance() rows of different fits can be told apart when stacked", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # The rows of a two-stage, a FIML, and a listwise fit looked alike in
  # glance(), though AIC and BIC of the first are the two-stage criteria of
  # Lai (2021), comparable only with two-stage fits, and the tidy() rows of
  # a fit on the t reference looked like those on the normal reference.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  f <- list(visual = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
            verbal = c("t6_paragraph_comprehension", "t7_sentence",
                       "t9_word_meaning"))
  set.seed(1)
  d <- hs[c(unlist(f), "t5_general_information")]
  d$t2_cubes[sample(nrow(d), 70)] <- NA
  d$t7_sentence[sample(nrow(d), 60)] <- NA
  quiet <- function(x) suppressMessages(suppressWarnings(x))
  g <- rbind(
    generics::glance(quiet(fa_confirmatory(d, f, missing = "two_stage",
                                 aux = "t5_general_information"))),
    generics::glance(quiet(fa_confirmatory(d, f, missing = "robust_two_stage",
                                 aux = "t5_general_information"))),
    generics::glance(quiet(fa_confirmatory(d, f, missing = "fiml",
                                 aux = "t5_general_information"))),
    generics::glance(quiet(fa_confirmatory(d, f, missing = "fiml", estimator = "MLR"))),
    generics::glance(quiet(fa_confirmatory(d, f, missing = "listwise"))),
    generics::glance(fa_confirmatory(S = stats::cov(hs[unlist(f)]), N = nrow(hs),
                           factors = f)))
  expect_identical(g$missing, c("two_stage", "robust_two_stage", "fiml",
                                "fiml", "listwise", "none"))
  expect_identical(g$estimator, c("ML", "ML", "ML", "MLR", "ML", "ML"))
  expect_match(g$test[1:2], "^two-stage ")
  expect_identical(g$test[3:6], c("standard", "scaled", "standard",
                                  "standard"))
  expect_identical(g$fit_indices, c("robust", "robust", "robust", "robust",
                                    "standard", "standard"))
  expect_identical(g$rmsea_ci_level, rep(0.90, 6))
  # A saturated model's index rows are those of exact fit, read from the
  # standard keys.
  expect_identical(generics::glance(quiet(fa_confirmatory(d, f["verbal"],
                                                missing = "fiml")))$fit_indices,
                   "standard")
  # tidy(): the reference of each row, also where it differs within a
  # table (the keyword puts the loadings and omega rows alone on t).
  t_num <- generics::tidy(fa_confirmatory(hs[1:40, unlist(f)], f, wald_df = 20))
  t_inf <- generics::tidy(fa_confirmatory(hs[1:40, unlist(f)], f))
  both <- rbind(t_inf, t_num)
  expect_identical(both$wald_df, rep(c(Inf, 20), each = nrow(t_inf)))
  k <- quiet(fa_confirmatory(d, f, missing = "fiml", wald_df = "missing_information"))
  tk <- generics::tidy(k)
  on_t <- names(attr(k, "wald_df"))
  expect_equal(tk$wald_df[match(on_t, tk$term)],
               unname(attr(k, "wald_df")))
  expect_true(all(is.finite(tk$wald_df[tk$term %in% on_t])))
  expect_true(all(tk$wald_df[!tk$term %in% on_t] == Inf))
  expect_true(all(startsWith(on_t, "lambda_") | startsWith(on_t, "omega_")))
})


# ---------------------------------------------------------------------------
# Ordered-categorical items and the categorical per-factor omega
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

test_that("ordered single-factor fa_confirmatory matches reliability_omega_categorical", {
  skip_on_cran()
  msgs <- testthat::capture_messages(
    r <- fa_confirmatory(items_ord[, 1:4], list(F1 = paste0("o", 1:4)),
               ordered = TRUE, output = "measurement"))
  expect_true(any(grepl("WLSMV", msgs)))
  expect_true(any(grepl("categorical sum score", msgs)))
  o_k <- as.data.frame(r)$estimate[as.data.frame(r)$term == "omega_F1"]
  ref <- suppressMessages(
    reliability_omega_categorical(data = items_ord[, 1:4],
                                  ci_method = "none"))
  expect_equal(o_k, ref$value[ref$term == "estimate"], tolerance = 1e-6)
})

test_that("ordered multi-factor fa_confirmatory reports GY omegas with NA intervals", {
  skip_on_cran()  # WLSMV on eight ordered items; the ordered anchor above runs on CRAN
  r <- suppressMessages(
    fa_confirmatory(items_ord, factors_ord, ordered = TRUE, output = "measurement"))
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
  skip_on_cran()  # Additional model fits run in the full local suite.
  items_mixed <- cbind(as.data.frame(items_ord[, 1:4]),
                       as.data.frame(lat_o[, 5:8]))
  names(items_mixed) <- c(paste0("o", 1:4), paste0("c", 1:4))
  fac_mixed <- list(F1 = paste0("o", 1:4), F2 = paste0("c", 1:4))
  r <- suppressMessages(
    fa_confirmatory(items_mixed, fac_mixed, ordered = paste0("o", 1:4),
          output = "measurement"))
  d <- as.data.frame(r)
  expect_true(is.na(d$se[d$term == "omega_F1"]))
  expect_false(is.na(d$se[d$term == "omega_F2"]))
  expect_identical(unname(attr(r, "omega_metric")),
                   c("categorical_sum_score", "linear"))
})

test_that("ordered input validation fails loudly", {
  expect_error(
    fa_confirmatory(items_ord, factors_ord, ordered = paste0("o", 1:3)),
    regexp = "mixes ordered and continuous")
  expect_error(
    fa_confirmatory(items_ord, factors_ord, ordered = c("o1", "nope")),
    regexp = "not in 'factors'")
  expect_error(
    fa_confirmatory(S = cov(items_ord), factors = factors_ord, N = N_ord, ordered = TRUE),
    regexp = "require raw data")
  expect_error(
    suppressMessages(fa_confirmatory(items_ord, factors_ord, ordered = TRUE,
                           equal_intercept = TRUE)),
    regexp = "Thresholds replace intercepts")
  expect_error(
    suppressMessages(fa_confirmatory(items_ord, factors_ord, ordered = TRUE,
                           missing = "fiml")),
    regexp = "pairwise")
})

test_that("continuous models record the linear omega metric", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  r <- fa_confirmatory(S = cov(lat_o), factors = list(F1 = paste0("V", 1:4)), N = N_ord,
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
  # fa_confirmatory() now writes free nu lines under FIML, and the two functions
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
  rk <- fa_confirmatory(X_shifted, list(f = paste0("y", seq_len(J))), missing = "ml")
  expect_equal(rk$estimate[rk$term == "omega_f"],
               r1$estimate[r1$term == "omega_f1"])
  expect_equal(rk$se[rk$term == "omega_f"],
               r1$se[r1$term == "omega_f1"])

  # The impossible combination is refused rather than silently mis-fit.
  expect_error(
    fa_confirmatory(X_shifted, list(f = paste0("y", seq_len(J))), missing = "ml",
          meanstructure = FALSE),
    "cannot be combined with", fixed = TRUE
  )
})

test_that("fa_confirmatory() signals a Heywood case with the classed condition cfa_1() uses", {
  skip_on_cran()  # Additional model fits run in the full local suite.
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
    res <- fa_confirmatory(S = S, N = 200, factors = list(f = paste0("y", 1:3))),
    warning = function(w) {
      classes <<- c(classes, class(w)[1L])
      invokeRestart("muffleWarning")
    }
  )
  expect_true("dmar_heywood_warning" %in% classes)
  expect_s3_class(res, "dmar_fa_confirmatory")
})

test_that("an ordered factor with an undetermined polychoric reports omega as NA with a warning", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  set.seed(113)
  N <- 500
  p_top <- c(.011, .027, .047, .118, .797)
  cut_top <- c(-Inf, stats::qnorm(cumsum(p_top)[-5]), Inf)
  cut_bottom <- c(-Inf, stats::qnorm(cumsum(rev(p_top))[-5]), Inf)
  latent <- outer(stats::rnorm(N), rep(0.8, 4)) +
    matrix(stats::rnorm(N * 4, sd = 0.6), N, 4)
  it <- cbind(y1 = cut(latent[, 1], cut_top, labels = FALSE),
              y2 = cut(latent[, 2], cut_top, labels = FALSE),
              y3 = cut(latent[, 3], cut_bottom, labels = FALSE),
              y4 = cut(latent[, 4], cut_bottom, labels = FALSE))
  expect_warning(
    r <- suppressMessages(fa_confirmatory(it, list(G = paste0("y", 1:4)),
                                ordered = TRUE, output = "measurement")),
    "omega_G is NA: the data do not determine the polychoric correlation of y1 and y3")
  rr <- as.data.frame(r)
  expect_true(is.na(rr$estimate[rr$term == "omega_G"]))
})

test_that("an ordered factor with an improper solution reports omega, AVE, and H as NA with one warning", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # Four symmetric items whose first loading runs past one in this sample
  # (the improper sample of test-reliability_omega_categorical.R). The
  # categorical omega assembly evaluates the bivariate normal at each
  # communality, so through 1.0.1 fa_confirmatory() stopped with mvtnorm's
  # "'corr' is not a correlation matrix". Once omega became NA the same
  # table still reported H_G = 1.079 and an AVE averaging a communality
  # of 1.057, neither of them a proportion of variance.
  set.seed(1)
  lam <- c(.95, .95, .6, .6)
  f <- stats::rnorm(80)
  ys <- sapply(lam, function(l) l * f + sqrt(1 - l^2) * stats::rnorm(80))
  it <- apply(ys, 2, function(y) findInterval(y, c(-1, 0, 1)) + 1L)
  colnames(it) <- paste0("y", 1:4)
  w <- testthat::capture_warnings(
    r <- suppressMessages(fa_confirmatory(it, list(G = paste0("y", 1:4)),
                                ordered = TRUE, output = "measurement")))
  rr <- as.data.frame(r)
  for (term in c("omega_G", "ave_G", "H_G")) {
    row <- rr[rr$term == term, c("estimate", "se", "test_statistic", "p_value",
                                 "ci_lower", "ci_upper")]
    expect_true(all(is.na(unlist(row))), label = term)
  }
  improper_msg <- grepl("are NA: the solution is improper", w, fixed = TRUE)
  expect_identical(sum(improper_msg), 1L)
  expect_match(w[improper_msg],
    "^omega_G, ave_G, and H_G are NA: the solution is improper \\(a Heywood case\\): the latent response behind y1 has a negative residual variance \\(its communality is 1\\.[0-9]{3}\\)\\. Categorical omega is not defined")
  # The fit itself keeps lavaan's defined parameters, H above one among
  # them, and H is at least the communality, since y1 alone attains it.
  fit <- suppressWarnings(suppressMessages(
    fa_confirmatory(it, list(G = paste0("y", 1:4)), ordered = TRUE, output = "fit")))
  pe <- lavaan::parameterEstimates(fit)
  H_fit <- pe$est[pe$lhs == "H_G" & pe$op == ":="]
  expect_gt(H_fit, 1)
  # The helper itself returns NA with the reason and the communalities.
  v <- .omega_categorical_from_kfit(fit, paste0("y", 1:4))
  expect_true(is.na(v))
  expect_identical(attr(v, "dmar_reason"), "improper")
  expect_identical(names(attr(v, "improper")), "y1")
  expect_gt(attr(v, "improper")[["y1"]], 1)
  expect_gte(H_fit, attr(v, "improper")[["y1"]])
  # .cfa_k_improper() finds the same item and communality from the
  # residual variance.
  imp <- .cfa_k_improper(pe, list(G = paste0("y", 1:4)))
  expect_named(imp, "G")
  expect_equal(imp$G, attr(v, "improper"), tolerance = 1e-6)
})

test_that("a continuous factor with a negative error variance reports AVE and H as NA, and omega while the error variances sum above zero", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # The just-identified Heywood case of the classed-condition test above:
  # lambda_1^2 = 1.28 and psi_1 = -0.28, so y1's communality is 1.28.
  # Through 1.0.1 the table reported H_f = 1.636 and an AVE of .760 that
  # averaged that communality; omega, .900, rests on an error variance
  # sum of .72 and is still a proportion.
  S <- matrix(c(1, 0.8, 0.8,
                0.8, 1, 0.5,
                0.8, 0.5, 1), 3, 3,
              dimnames = list(paste0("y", 1:3), paste0("y", 1:3)))
  w <- testthat::capture_warnings(
    r <- fa_confirmatory(S = S, N = 200, factors = list(f = paste0("y", 1:3))))
  rr <- as.data.frame(r)
  expect_true(all(is.na(unlist(rr[rr$term %in% c("ave_f", "H_f"),
                                  c("estimate", "se", "ci_lower",
                                    "ci_upper")]))))
  expect_equal(rr$estimate[rr$term == "omega_f"], 0.9, tolerance = 1e-6)
  expect_true(is.finite(rr$ci_lower[rr$term == "omega_f"]))
  # The classed Heywood warning is kept, and one warning names the rows.
  expect_true(any(grepl("Negative error variance estimate (Heywood case) for: y1",
                        w, fixed = TRUE)))
  improper_msg <- grepl("are NA: the solution is improper", w, fixed = TRUE)
  expect_identical(sum(improper_msg), 1L)
  expect_match(w[improper_msg],
    "^ave_f and H_f are NA: the solution is improper \\(a Heywood case\\): y1 has a negative error variance \\(its communality is 1\\.280\\)\\.")
  expect_match(w[improper_msg], "omega_f is reported", fixed = TRUE)
  # lavaan's own defined parameter, the value the table no longer shows:
  # H = q / (1 + q) with q = sum(lambda^2 / psi) = 1.28 / -0.28 + 2, which
  # is 1.636, at least the communality 1.28.
  fit <- suppressWarnings(
    fa_confirmatory(S = S, N = 200, factors = list(f = paste0("y", 1:3)),
          output = "fit"))
  pe <- lavaan::parameterEstimates(fit)
  q <- 1.28 / -0.28 + 2
  expect_equal(pe$est[pe$lhs == "H_f" & pe$op == ":="], q / (1 + q),
               tolerance = 1e-6)
})

test_that("a two-item factor of holzinger_swineford with a Heywood case loses AVE and H, and omega when the error variances sum below zero", {
  # The two fixtures above carry the same checks on CRAN; these fits of
  # the package's data anchor the values NEWS reports.
  skip_on_cran()
  hs <- holzinger_swineford
  # Visual perception and lozenges beside three verbal tests: visual
  # perception's error variance is negative, and through 1.0.1 the
  # measurement table reported H_f1 = 1.321 [0.166, 2.476] and an AVE of
  # .726; omega_f1, .707, rests on a positive error variance sum.
  w <- testthat::capture_warnings(r <- fa_confirmatory(hs, list(
    f1 = c("t1_visual_perception", "t4_lozenges"),
    f2 = c("t5_general_information", "t6_paragraph_comprehension",
           "t7_sentence")), output = "measurement"))
  rr <- as.data.frame(r)
  est <- stats::setNames(rr$estimate, rr$term)
  expect_true(all(is.na(est[c("ave_f1", "H_f1")])))
  expect_equal(unname(est["omega_f1"]), 0.7067640, tolerance = 1e-5)
  # The proper factor keeps every row.
  expect_true(all(is.finite(est[c("omega_f2", "ave_f2", "H_f2")])))
  expect_true(any(grepl("^ave_f1 and H_f1 are NA: .* t1_visual_perception has a negative error variance",
                        w)))
  # Visual perception and cubes beside three speed tests: here the error
  # variances of the factor sum below zero, and 1.0.1 reported omega_f1 =
  # 1.050, H_f1 = 1.536, and an AVE of .788; all three are NA.
  w <- testthat::capture_warnings(r <- fa_confirmatory(hs, list(
    f1 = c("t1_visual_perception", "t2_cubes"),
    f2 = c("t10_addition", "t11_code", "t12_counting_groups_of_dots"))))
  rr <- as.data.frame(r)
  est <- stats::setNames(rr$estimate, rr$term)
  expect_lt(est[["error_sum_f1"]], 0)
  expect_true(all(is.na(est[c("omega_f1", "ave_f1", "H_f1")])))
  expect_true(all(is.finite(est[c("omega_f2", "ave_f2", "H_f2")])))
  expect_true(any(grepl("^omega_f1, ave_f1, and H_f1 are NA: .* The error variances sum to a negative value",
                        w)))
  # A proper solution is untouched, with no such warning.
  expect_silent(r <- fa_confirmatory(hs, list(
    f1 = c("t5_general_information", "t6_paragraph_comprehension",
           "t7_sentence")), output = "measurement"))
  expect_true(all(is.finite(r$estimate)))
})


test_that("fa_confirmatory() refuses a group passed through ..., whose labels lavaan would hold equal across the groups", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # Through 1.0.1, fa_confirmatory(..., group = "school") fit one set of labeled
  # parameters for both schools, with lavaan's warning about a single
  # label in a multiple group setting, and its table repeated each labeled
  # row once per school with the same estimate.
  hs <- holzinger_swineford
  f <- list(verbal = c("t6_paragraph_comprehension", "t7_sentence",
                       "t9_word_meaning"),
            deduction = c("t20_deduction", "t22_problem_reasoning",
                          "t23_series_completion"))
  msg <- "'group' is not available in fa_confirmatory()"
  expect_error(fa_confirmatory(hs, f, group = "school"), msg, fixed = TRUE)
  expect_error(fa_confirmatory(hs, f, group = "school", output = "fit"), msg,
               fixed = TRUE)
  # A partial name would reach lavaan's 'group' formal through do.call().
  expect_error(fa_confirmatory(hs, f, gr = "school"), msg, fixed = TRUE)
  expect_error(cfa_1(hs, items = f$verbal, group = "school"), msg,
               fixed = TRUE)
  expect_error(cfa_2(hs, f$verbal, f$deduction, group = "school"), msg,
               fixed = TRUE)
  expect_error(fa_confirmatory(hs, f, group = "school"), "measurement_invariance()",
               fixed = TRUE)
})

test_that("the route the group refusal names gives one group's H", {
  # The configural model of measurement_invariance() frees every parameter
  # in each group, so H for one group read from its fit is H from that
  # group's rows alone.
  skip_on_cran()
  hs <- holzinger_swineford
  f <- list(verbal = c("t6_paragraph_comprehension", "t7_sentence",
                       "t9_word_meaning"),
            deduction = c("t20_deduction", "t22_problem_reasoning",
                          "t23_series_completion"))
  mi <- measurement_invariance(hs, model = f, group = "school",
                               levels = "configural")
  h_group <- reliability_H(fit = attr(mi, "fits")$configural,
                           factor = "verbal", group = "Pasteur")
  alone <- fa_confirmatory(hs[hs$school == "Pasteur", ], f, output = "measurement")
  expect_equal(h_group$value[h_group$term == "estimate"],
               alone$estimate[alone$term == "H_verbal"], tolerance = 1e-5)
})


# ---- wald_df: a t reference distribution for the Wald intervals ------------

test_that("fa_confirmatory() with wald_df = Inf returns lavaan's intervals and p-values untouched", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  a <- fa_confirmatory(hs, cfa_k_hs_factors())
  expect_identical(a, fa_confirmatory(hs, cfa_k_hs_factors(), wald_df = Inf))
  expect_identical(attr(a, "wald_df"), Inf)
})

test_that("wald_df = Inf leaves the other tables and the wrappers untouched", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # six fits; the verbose table's identity above runs on CRAN
  hs <- cfa_k_hs_data()
  for (out in c("measurement", "standardized")) {
    a <- fa_confirmatory(hs, cfa_k_hs_factors(), output = out)
    b <- fa_confirmatory(hs, cfa_k_hs_factors(), output = out, wald_df = Inf)
    expect_identical(a, b)
    expect_identical(attr(a, "wald_df"), Inf)
  }
  # The wrappers forward the argument, and the default there too is the
  # normal reference.
  sp <- hs[, c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
               "t4_lozenges")]
  expect_identical(cfa_1(sp), cfa_1(sp, wald_df = Inf))
})

test_that("a finite wald_df refers fa_confirmatory()'s intervals and p-values to t", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # seven fits; the arithmetic is checked on CRAN in
  # test-wald_df_internals.R and test-reliability_alpha.R
  hs <- cfa_k_hs_data()
  df <- 20
  a <- fa_confirmatory(hs, cfa_k_hs_factors())
  r <- fa_confirmatory(hs, cfa_k_hs_factors(), wald_df = df)
  crit <- stats::qt(0.975, df)
  # Estimates, standard errors, the z ratios, and every fit row are
  # unchanged; the limits and p-values of the parameter rows move.
  expect_identical(r[c("syntax", "term", "estimate", "se", "test_statistic")],
                   a[c("syntax", "term", "estimate", "se", "test_statistic")])
  fit_rows <- is.na(a$se)
  for (col in names(a)) {
    expect_identical(r[[col]][fit_rows], a[[col]][fit_rows])
  }
  # The omega, AVE, and H rows form their limits on the logit scale (with
  # the same t quantile); every other parameter row is the estimate plus or
  # minus the t quantile times the standard error.
  h_rows <- grepl("^(omega|ave|H)_", a$term)
  free <- !is.na(a$test_statistic) & !h_rows
  expect_equal(r$ci_lower[free], a$estimate[free] - crit * a$se[free],
               tolerance = 1e-12)
  expect_equal(r$ci_upper[free], a$estimate[free] + crit * a$se[free],
               tolerance = 1e-12)
  for (i in which(h_rows)) {
    lim <- DMAR:::.logistic_ci(a$estimate[i], a$se[i], crit)
    expect_equal(r$ci_lower[i], lim$lower, tolerance = 1e-12)
    expect_equal(r$ci_upper[i], lim$upper, tolerance = 1e-12)
  }
  free <- !is.na(a$test_statistic)
  expect_identical(r$p_value[free], 2 * stats::pt(-abs(a$test_statistic[free]), df))
  expect_true(all(r$ci_upper[free] - r$ci_lower[free] >
                    a$ci_upper[free] - a$ci_lower[free]))
  expect_identical(attr(r, "wald_df"), df)
  printed <- capture.output(print(r))
  expect_true(any(grepl(paste0("Wald intervals and p-values: t reference ",
                               "on 20 degrees of freedom"), printed)))

  # The measurement table and the standardized solution follow suit.
  m <- fa_confirmatory(hs, cfa_k_hs_factors(), output = "measurement", wald_df = df)
  om <- m[m$term == "omega_verbal", ]
  expect_equal(om$ci_lower, DMAR:::.logistic_ci(om$estimate, om$se, crit)$lower,
               tolerance = 1e-12)
  expect_true(is.na(m$ci_lower[grepl("^htmt_", m$term)]))
  s0 <- fa_confirmatory(hs, cfa_k_hs_factors(), output = "standardized")
  s <- fa_confirmatory(hs, cfa_k_hs_factors(), output = "standardized", wald_df = df)
  ok <- !is.na(s0$z) & !(s0$op == ":=" & grepl("^(omega|ave|H)_", s0$lhs))
  expect_equal(s$ci.upper[ok], s0$est.std[ok] + crit * s0$se[ok],
               tolerance = 1e-12)
  hs_rows <- which(s0$op == ":=" & grepl("^(omega|ave|H)_", s0$lhs))
  for (i in hs_rows) {
    lim <- DMAR:::.logistic_ci(s0$est.std[i], s0$se[i], crit)
    expect_equal(s$ci.upper[i], lim$upper, tolerance = 1e-12)
  }
  ok <- !is.na(s0$z)
  expect_identical(s$pvalue[ok], 2 * stats::pt(-abs(s0$z[ok]), df))
  expect_identical(attr(s, "wald_df"), df)

  # The wrappers forward it.
  sp <- hs[, c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
               "t4_lozenges")]
  expect_identical(cfa_1(sp, wald_df = df),
                   fa_confirmatory(sp, list(f1 = names(sp)), wald_df = df))
})

test_that("fa_confirmatory() refuses a finite wald_df where nothing would change", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  expect_error(fa_confirmatory(hs, cfa_k_hs_factors(), output = "fit", wald_df = 20),
               "output = \"fit\" returns lavaan's own object")
  expect_error(fa_confirmatory(hs, cfa_k_hs_factors(), output = "summary",
                     wald_df = 20),
               "output = \"summary\" returns lavaan's own object")
  expect_error(fa_confirmatory(hs, cfa_k_hs_factors(), se = "none", wald_df = 20),
               "se = \"none\" gives no standard errors")
  expect_error(cfa_2(hs, cfa_k_hs_factors()$verbal,
                     cfa_k_hs_factors()$deduction, wald_df = -2),
               "single positive number")
  # The default stays accepted with every output.
  expect_s4_class(fa_confirmatory(hs, cfa_k_hs_factors(), output = "fit",
                        wald_df = Inf), "lavaan")
})

test_that("wald_df = 'missing_information' on complete data is t on N - 1 for the loadings", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # With no value missing every fraction of missing information is 0
  # exactly, so the calibrated rule max(3, N(1 - FMI) - 1) is N - 1, and
  # the loadings' intervals and p-values are those of wald_df = N - 1 to
  # the last bit. Every other row keeps the normal reference, omega, AVE,
  # and H with their logit-scale limits.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  f <- cfa_k_hs_factors()
  N <- nrow(hs)
  z <- fa_confirmatory(hs, f)
  r <- fa_confirmatory(hs, f, wald_df = "missing_information")
  n1 <- fa_confirmatory(hs, f, wald_df = N - 1)
  t_rows <- grepl("^lambda_", r$term)
  expect_identical(sum(t_rows), 6L)
  cols <- c("ci_lower", "ci_upper", "p_value")
  for (col in names(r)) {
    expect_identical(r[[col]][!t_rows], z[[col]][!t_rows], info = col)
    expect_identical(r[[col]][t_rows], n1[[col]][t_rows], info = col)
  }
  expect_false(identical(r$ci_lower[t_rows], z$ci_lower[t_rows]))
  wd <- attr(r, "wald_df")
  expect_identical(unname(wd), rep(N - 1, 6))
  expect_identical(names(wd), r$term[t_rows])
  tab <- attr(r, "wald_df_table")
  expect_s3_class(tab, "dmar_df_missing")
  expect_identical(tab$fmi, rep(0, 6))
  expect_identical(tab$k, rep(1, 6))
  expect_identical(attr(r, "wald_df_rows"), "loadings")
  printed <- capture.output(print(r))
  expect_true(any(grepl(paste0(
    "Wald intervals and p-values of the loadings: t reference on ",
    "300 degrees of freedom, from von Hippel's \\(2016\\)"), printed)))
  expect_true(any(grepl("N = 301 and FMI = 0 ", printed)))
  expect_true(any(grepl("the other rows keep the normal reference",
                        printed)))
  # Listwise deletion is complete data too, on the complete cases.
  sp <- c("t1_visual_perception", "t2_cubes", "t25_paper_form_board_r")
  lw <- fa_confirmatory(hs, list(s = sp), missing = "listwise",
                        wald_df = "missing_information")
  n_cc <- sum(stats::complete.cases(hs[sp]))
  expect_identical(unname(attr(lw, "wald_df")), rep(n_cc - 1, 3))
})

test_that("wald_df = 'missing_information' under FIML is max(3, N(1 - FMI) - 1) per loading", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # FIML fits with the FMI's two refits
  hs <- cfa_k_hs_data()
  # Flags and the second paper form board were given to the Grant-White
  # pupils only, so 156 of the 301 values of each are missing.
  f <- list(spatial = c("t1_visual_perception", "t2_cubes", "t4_lozenges",
                        "t25_paper_form_board_r", "t26_flags"))
  r <- fa_confirmatory(hs, f, missing = "fiml", wald_df = "missing_information")
  fit <- fa_confirmatory(hs, f, missing = "fiml", output = "fit")
  fm <- fraction_missing_information(fit)
  N <- lavaan::lavInspect(fit, "nobs")
  key <- paste("spatial =~", f$spatial)
  fmi <- fm$fmi[match(key, fm$term)]
  expect_false(anyNA(fmi))
  wd <- attr(r, "wald_df")
  expect_identical(names(wd), paste0("lambda_spatial_", 1:5))
  expect_equal(unname(wd), pmax(3, N * (1 - fmi) - 1), tolerance = 1e-10)
  # The two incomplete items lose about half their information, and their
  # degrees of freedom with it.
  expect_true(all(wd[4:5] < 160) && all(wd[1:3] > 240))
  # Each covered row's interval is its own t interval; the rest are
  # lavaan's.
  z <- fa_confirmatory(hs, f, missing = "fiml")
  rows <- match(names(wd), r$term)
  expect_equal(r$ci_upper[rows],
               z$estimate[rows] + stats::qt(0.975, unname(wd)) * z$se[rows],
               tolerance = 1e-12)
  expect_identical(r$p_value[rows],
                   2 * stats::pt(-abs(z$test_statistic[rows]), unname(wd)))
  others <- setdiff(seq_len(nrow(r)), rows)
  expect_identical(r[others, c("ci_lower", "ci_upper", "p_value")],
                   z[others, c("ci_lower", "ci_upper", "p_value")])
  printed <- capture.output(print(r))
  expect_true(any(grepl(paste0("t reference on [0-9.]+ to [0-9.]+ degrees ",
                               "of freedom, by parameter"), printed)))
  # The omega rows keep the normal reference with their logit-scale
  # limits, where the calibration found the t reference overcovered; the
  # measurement table, which has no loading, refuses the keyword.
  om <- r[r$term == "omega_spatial", ]
  lim <- DMAR:::.logistic_ci(om$estimate, om$se, stats::qnorm(0.975))
  expect_equal(c(om$ci_lower, om$ci_upper), c(lim$lower, lim$upper),
               tolerance = 1e-12)
  expect_error(fa_confirmatory(hs, f, missing = "fiml", output = "measurement",
                               wald_df = "missing_information"),
               "which that table does not show")
  # The wrappers forward it.
  expect_identical(cfa_1(hs, items = f$spatial, missing = "fiml",
                         wald_df = "missing_information")$ci_upper,
                   fa_confirmatory(hs, list(f1 = f$spatial), missing = "fiml",
                         wald_df = "missing_information")$ci_upper)
})

test_that("cfa_2() forwards wald_df = 'missing_information'", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  f <- cfa_k_hs_factors()
  a <- cfa_2(hs, f$verbal, f$deduction, wald_df = "missing_information")
  b <- fa_confirmatory(hs, list(f1 = f$verbal, f2 = f$deduction),
             wald_df = "missing_information")
  expect_identical(a, b)
})

test_that("fa_confirmatory() refuses wald_df = 'missing_information' outside its calibration", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  f <- cfa_k_hs_factors()
  K <- "missing_information"
  expect_error(fa_confirmatory(hs, f, output = "fit", wald_df = K),
               "the keyword applies to output = \"verbose\"")
  expect_error(fa_confirmatory(hs, f, output = "measurement", wald_df = K),
               "which that table does not show")
  expect_error(fa_confirmatory(hs, f, output = "standardized", wald_df = K),
               "evaluated the unstandardized loadings")
  expect_error(fa_confirmatory(hs, f, se = "none", wald_df = K),
               "gives no standard errors")
  expect_error(fa_confirmatory(hs, f, estimator = "WLS", wald_df = K),
               "estimator = \"WLS\": the rule takes")
  expect_error(fa_confirmatory(hs, f, ordered = TRUE, wald_df = K),
               "with ordered items")
  expect_error(fa_confirmatory(hs, f, missing = "pairwise", wald_df = K),
               "pairwise present data")
  expect_error(fa_confirmatory(hs, f, missing = "two_stage", wald_df = K),
               "DMAR does not compute for a two-stage fit")
  # Backstops on the fit itself: standard errors that came from an
  # argument passed on to lavaan.
  expect_error(suppressWarnings(
    fa_confirmatory(hs, f, estimator = "MLR", se = "robust.huber.white",
          cluster = "school", wald_df = K)),
    "its standard errors are lavaan's se = \"robust.cluster\"")
})

test_that("wald_df = 'missing_information' refuses robust.sem under FIML for the reason lavaan gives", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # lavaan computes the robust.sem sandwich for complete data only: asked
  # for it with missing = "ml", it sets missing to "listwise" (with a
  # warning) and fits the complete cases. The refusal said lavaan
  # evaluated the sandwich at the expected information under FIML, which
  # it never computes; it now says what lavaan does, whether or not a
  # value is missing, and points to the two fits that are accepted.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  K <- "missing_information"
  spatial <- list(s = c("t1_visual_perception", "t2_cubes", "t4_lozenges",
                        "t25_paper_form_board_r", "t26_flags"))
  msg <- 'se = "robust.sem" requires complete data and is incompatible with FIML'
  expect_error(fa_confirmatory(hs, spatial, estimator = "MLR", se = "robust.sem",
    missing = "fiml", output = "fit"), msg, fixed = TRUE)
  expect_error(fa_confirmatory(hs, spatial, estimator = "MLR", se = "robust.sem",
                     missing = "fiml", wald_df = K), msg)
  # Three items no pupil missed: lavaan's fallback drops nothing, and the
  # request is refused all the same, with the fit that is accepted named.
  complete3 <- list(s = spatial$s[1:3])
  expect_error(fa_confirmatory(hs, complete3, estimator = "MLR", se = "robust.sem",
                     missing = "fiml", wald_df = K),
               "robust.huber.white")
  lw <- fa_confirmatory(hs, complete3, estimator = "MLR", se = "robust.sem",
              wald_df = K)
  expect_identical(unname(attr(lw, "wald_df")), rep(nrow(hs) - 1, 3))
  # From a covariance matrix lavaan computes robust.sem only with the
  # asymptotic covariance matrix of the moments as its NACOV, which the
  # help now says (it said robust.sem was accepted from a covariance
  # matrix, which S alone is not, by lavaan). With NACOV the intervals are
  # those of the raw data listwise fit.
  cc <- hs[stats::complete.cases(hs[spatial$s]), spatial$s]
  S <- stats::cov(cc)
  expect_error(fa_confirmatory(S = S, N = nrow(cc), factors = spatial,
                     se = "robust.sem", wald_df = K), "NACOV")
  G <- lavaan::lavInspect(
    lavaan::cfa(paste("s =~", paste(spatial$s, collapse = " + ")), data = cc,
                se = "robust.sem", meanstructure = FALSE), "gamma")
  from_s <- fa_confirmatory(S = S, N = nrow(cc), factors = spatial, se = "robust.sem",
                  NACOV = G, wald_df = K)
  from_raw <- fa_confirmatory(cc, spatial, se = "robust.sem", wald_df = K)
  expect_identical(unname(attr(from_s, "wald_df")),
                   rep(nrow(cc) - 1, length(attr(from_s, "wald_df"))))
  keep <- startsWith(from_raw$term, "lambda")
  expect_equal(from_s$ci_lower[startsWith(from_s$term, "lambda")],
               from_raw$ci_lower[keep], tolerance = 1e-6)
})

test_that("an MLR fit uses matching robust errors unless an explicit method is requested", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # fa_confirmatory() always passes se to lavaan, so lavaan's own pairing of MLR with
  # the Huber-White sandwich does not apply; the help of 'se' and of
  # 'estimator' said the opposite, that MLR's standard errors were
  # corrected for nonnormality.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  spatial <- list(s = c("t1_visual_perception", "t2_cubes", "t4_lozenges",
                        "t25_paper_form_board_r", "t26_flags"))
  q <- function(...) suppressWarnings(fa_confirmatory(hs, spatial, estimator = "MLR",
                                            missing = "fiml", output = "fit",
                                            ...))
  expect_identical(lavaan::lavInspect(q(), "options")$se, "robust.huber.white")
  expect_identical(lavaan::lavInspect(q(se = "robust.huber.white"),
                                      "options")$se, "robust.huber.white")
  own <- lavaan::cfa(paste("s =~", paste(spatial$s, collapse = " + ")),
                     data = hs, estimator = "MLR", missing = "ml",
                     std.lv = TRUE)
  expect_identical(lavaan::lavInspect(own, "options")$se, "robust.huber.white")
  # se = "standard" is not normal theory whatever the estimator, as the
  # help of 'se' said after the first repair: with ordered items fa_confirmatory()
  # requests robust.sem in its place, and the standard errors of WLS are
  # the asymptotically distribution free ones, (D' Gamma^-1 D)^-1 / N,
  # which differ from the normal theory standard errors of the same
  # estimates.
  sp4 <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
           "t4_lozenges")
  h4 <- hs
  for (v in sp4) h4[[v]] <- as.integer(cut(h4[[v]], 4))
  ord <- suppressMessages(suppressWarnings(
    fa_confirmatory(h4, list(s = sp4), ordered = TRUE, output = "fit")))
  expect_identical(lavaan::lavInspect(ord, "options")$se, "robust.sem")
  wls <- suppressWarnings(fa_confirmatory(hs, list(s = sp4), estimator = "WLS",
                                output = "fit"))
  expect_identical(lavaan::lavInspect(wls, "options")$se, "standard")
  at <- lavaan::cfa(paste("s =~", paste(sp4, collapse = " + ")), data = hs,
                    std.lv = TRUE, estimator = "ML", start = wls,
                    optim.method = "none", se = "standard",
                    information = "expected")
  pw <- lavaan::parameterEstimates(wls)
  pa <- lavaan::parameterEstimates(at)
  key <- function(p) paste(p$lhs, p$op, p$rhs)
  load <- pw$op == "=~"
  expect_equal(pa$est[match(key(pw), key(pa))][load], pw$est[load],
               tolerance = 1e-10)
  expect_gt(max(abs(pw$se[load] / pa$se[match(key(pw), key(pa))][load] - 1)),
            0.1)
})

test_that("from a covariance matrix MLR and the Huber-White sandwich need raw data", {
  # lavaan 0.7-3 fits no MLR and computes no Huber-White sandwich from a
  # covariance matrix, with or without the asymptotic covariance matrix of
  # the moments as its NACOV: both need raw data, where the help of 'se'
  # after the first repair said that NACOV enabled MLR too. robust.sem with
  # ML runs with NACOV. The expectations pin lavaan's behavior that the
  # help states, so that a lavaan that changes it shows here.
  skip_if_not_installed("lavaan")
  skip_on_cran()
  hs <- cfa_k_hs_data()
  spatial <- list(s = c("t1_visual_perception", "t2_cubes", "t4_lozenges",
                        "t25_paper_form_board_r", "t26_flags"))
  cc <- hs[stats::complete.cases(hs[spatial$s]), spatial$s]
  S <- stats::cov(cc)
  G <- lavaan::lavInspect(
    lavaan::cfa(paste("s =~", paste(spatial$s, collapse = " + ")), data = cc,
                se = "robust.sem", meanstructure = FALSE), "gamma")
  for (se in c("standard", "robust.sem")) {
    expect_error(fa_confirmatory(S = S, N = nrow(cc), factors = spatial,
                       estimator = "MLR", se = se, NACOV = G,
                       output = "fit"), info = se)
  }
  expect_error(fa_confirmatory(S = S, N = nrow(cc), factors = spatial,
                     se = "robust.huber.white", NACOV = G, output = "fit"),
               "requires full data")
  ml <- fa_confirmatory(S = S, N = nrow(cc), factors = spatial, se = "robust.sem",
              NACOV = G, output = "fit")
  expect_identical(lavaan::lavInspect(ml, "options")$se, "robust.sem")
})

test_that("wald_df = 'missing_information' refuses every fit in which lavaan replaced FIML with listwise deletion", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # The 1.1.0 review: lavaan also fits the complete cases in place of the
  # FIML fit the call asks for when a test statistic of the
  # Satorra-Bentler kind reaches it through ..., and the keyword then
  # reported FMI = 0 and N - 1 of the complete cases for every loading
  # and omega (240 on these data, where the FIML fit gives 231.4 to
  # 285.5), and accepted the expected information it refuses under FIML.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  K <- "missing_information"
  sp <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
          "t4_lozenges")
  f <- list(spatial = sp)
  set.seed(1)
  d <- hs
  d$t1_visual_perception[sample(nrow(d), 60)] <- NA
  q <- function(expr) suppressWarnings(expr)
  # Refused before the fit, whatever the estimator or the output.
  for (test in c("satorra.bentler", "mean.var.adjusted", "scaled.shifted")) {
    expect_error(q(fa_confirmatory(d, f, missing = "fiml", test = test, wald_df = K)),
                 paste0("is not available with test = \"", test,
                        "\" and missing = \"fiml\": lavaan does not compute ",
                        "that test under full information maximum ",
                        "likelihood but falls back to listwise deletion"),
                 fixed = TRUE)
  }
  expect_error(q(fa_confirmatory(d, f, missing = "fiml", test = "mean.var.adjusted",
                       output = "measurement", wald_df = K)),
               "which that table does not show", fixed = TRUE)
  expect_error(q(fa_confirmatory(d, f, missing = "fiml", estimator = "MLR",
                       test = "satorra.bentler", wald_df = K)),
               "Use estimator = \"MLR\" with missing = \"fiml\"",
               fixed = TRUE)
  # With no value missing the request is refused all the same, as
  # robust.sem is.
  expect_error(q(fa_confirmatory(hs, f, missing = "fiml", test = "satorra.bentler",
                       wald_df = K)),
               "falls back to listwise deletion", fixed = TRUE)
  # Any other route is caught on the fit itself: lavaan's own spelling
  # of the test, which the check before the fit does not know.
  expect_error(q(fa_confirmatory(d, f, missing = "fiml", test = "SB", wald_df = K)),
               "lavaan changed missing-data handling to listwise",
               fixed = TRUE)
  # With no value missing that route is refused too, as the check before
  # the fit refuses it, and the stop says that no case was set aside; it
  # had said that lavaan fit the complete cases instead, "301 of the 301
  # cases read".
  msg <- tryCatch(q(fa_confirmatory(hs, f, missing = "fiml", test = "SB",
                          wald_df = K)), error = conditionMessage)
  expect_match(msg, "lavaan changed missing-data handling to listwise", fixed = TRUE)
  expect_match(msg, "No result is returned from the changed case set", fixed = TRUE)
  # Accepted: the test with listwise deletion, a covariance matrix (which
  # lavaan fits as complete data), and FIML with lavaan's default test.
  lw <- q(fa_confirmatory(d, f, missing = "listwise", test = "satorra.bentler",
                wald_df = K))
  expect_identical(unname(attr(lw, "wald_df")), rep(240, 4))
  cv <- fa_confirmatory(S = stats::cov(hs[sp]), N = nrow(hs), factors = f,
              missing = "fiml", wald_df = K)
  expect_identical(unname(attr(cv, "wald_df")), rep(300, 4))
  fiml <- fa_confirmatory(d, f, missing = "fiml", wald_df = K)
  expect_equal(unname(attr(fiml, "wald_df")),
               c(231.4, 285.5, 285.1, 278.3), tolerance = 1e-3)
  expect_output(print(fiml), "calibration study")
  expect_match(generics::tidy(fiml)$wald_calibration_note[1L], "affine-collapse screening")
  expect_match(generics::glance(fiml)$wald_calibration_note, "affine-collapse screening")
})

test_that("wald_df = 'missing_information' gives lavaan's coverage notice once and refuses a pair never observed together", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # The 1.1.0 review: the refit behind the keyword repeated lavaan's
  # notice of low coverage that the fit had already given, and with a
  # pair of items never observed together (a two-form planned missing
  # design) it repeated the zero coverage notice and returned degrees of
  # freedom resting on an arbitrary covariance.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  K <- "missing_information"
  f <- list(spatial = c("t1_visual_perception", "t2_cubes",
                        "t3_paper_form_board", "t4_lozenges"),
            verbal = c("t6_paragraph_comprehension", "t7_sentence",
                       "t9_word_meaning"))
  coverage_notices <- function(expr) {
    said <- new.env()
    said$w <- character(0)
    res <- withCallingHandlers(
      tryCatch(expr, error = identity),
      warning = function(w) {
        said$w <- c(said$w, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(result = res, n = sum(grepl("coverage", said$w, fixed = TRUE)))
  }
  set.seed(3)
  form <- sample(1:2, nrow(hs), replace = TRUE)
  # Ten pupils of the first form take both tests.
  both <- sample(which(form == 1), 10)
  low <- hs
  low$t1_visual_perception[form == 1 & !seq_len(nrow(hs)) %in% both] <- NA
  low$t6_paragraph_comprehension[form == 2] <- NA
  plain <- coverage_notices(fa_confirmatory(low, f, missing = "fiml"))
  keyword <- coverage_notices(fa_confirmatory(low, f, missing = "fiml", wald_df = K))
  expect_identical(plain$n, 1L)
  expect_identical(keyword$n, 1L)
  expect_false(inherits(keyword$result, "error"))
  # Never observed together: the fit stands, the keyword is refused.
  never <- hs
  never$t1_visual_perception[form == 1] <- NA
  never$t6_paragraph_comprehension[form == 2] <- NA
  plain <- coverage_notices(fa_confirmatory(never, f, missing = "fiml"))
  expect_false(inherits(plain$result, "error"))
  keyword <- coverage_notices(fa_confirmatory(never, f, missing = "fiml",
                                    wald_df = K))
  expect_identical(keyword$n, 1L)
  expect_s3_class(keyword$result, "dmar_fmi_never_together")
  expect_match(conditionMessage(keyword$result),
               "t1_visual_perception with t6_paragraph_comprehension",
               fixed = TRUE)
  expect_match(conditionMessage(keyword$result),
               "Pass a number of degrees of freedom", fixed = TRUE)
})

test_that("wald_df = 'missing_information' takes the information the calibration used, and names the one it refuses", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # Under FIML the study's standard errors came from the observed
  # information computed from the Hessian of the fitted model, lavaan's
  # default there. The check read only lavaan's 'information', so the
  # first order information was refused as "the expected information" and
  # the observed information of the unrestricted model
  # (observed.information = "h1") was accepted.
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  K <- "missing_information"
  spatial <- list(s = c("t1_visual_perception", "t2_cubes",
                        "t25_paper_form_board_r", "t26_flags"))
  expect_error(fa_confirmatory(hs, spatial, missing = "fiml",
                     information = "expected", wald_df = K),
               paste0("come from lavaan's information = \"expected\", and ",
                      "under full information maximum likelihood the ",
                      "calibration study used the observed information ",
                      "from the Hessian"))
  expect_error(fa_confirmatory(hs, spatial, missing = "fiml",
                     information = "first.order", wald_df = K),
               "come from lavaan's information = \"first.order\", and under")
  expect_error(fa_confirmatory(hs, spatial, missing = "fiml",
                     observed.information = "h1", wald_df = K),
               paste0("information = \"observed\" with ",
                      "observed.information = \"h1\", and under"))
  # With complete data the expected information, lavaan's default for ML
  # there, and the Hessian are taken; the first order information is not.
  complete3 <- list(s = c("t1_visual_perception", "t2_cubes",
                          "t4_lozenges"))
  N1 <- nrow(hs) - 1
  expect_identical(unname(attr(fa_confirmatory(hs, complete3, wald_df = K),
                               "wald_df")), rep(N1, 3))
  expect_identical(unname(attr(fa_confirmatory(hs, complete3, information = "observed",
                                     wald_df = K), "wald_df")), rep(N1, 3))
  expect_error(fa_confirmatory(hs, complete3, information = "first.order",
                     wald_df = K),
               paste0("information = \"first.order\", and with complete ",
                      "data the keyword takes the expected information"))
})

test_that("the measurement structure introduces the table in a knitted report too", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # print.dmar_fa_confirmatory() wrote the structure lines above the table itself,
  # so knit_print() and as_kable(), which never call print(), dropped
  # them. They are now the header block of the one footer every display
  # reads, and the console prints them exactly as before.
  skip_if_not_installed("lavaan")
  skip_if_not_installed("knitr")
  res <- fa_confirmatory(cfa_k_hs_data(), cfa_k_hs_factors(), equal_loading = TRUE)
  header <- c("Measurement structure, per factor:",
              "  verbal: essentially tau-equivalent (equal loadings)",
              "  deduction: essentially tau-equivalent (equal loadings)")
  printed <- capture.output(print(res))
  expect_identical(printed[1:4], c(header, ""))
  expect_match(printed[5], "^ syntax +term")

  footer <- DMAR:::.dmar_tbl_footer(res)
  expect_identical(unname(footer[names(footer) == "header"]), header)
  # The level of the parameter intervals, with the RMSEA interval's own.
  expect_true("Confidence level: 95% (90% for the RMSEA interval)" %in% footer)

  # knit_print() sets the lines as paragraphs before the table.
  knitted <- as.character(knit_print.dmar_tbl(res))
  expect_match(knitted, paste0(
    "^\n\nMeasurement structure, per factor:\n\n",
    "verbal: essentially tau-equivalent \\(equal loadings\\)\n\n",
    "deduction: essentially tau-equivalent \\(equal loadings\\)\n\n",
    "\\|syntax"))
  # as_kable() puts them first in the caption.
  kab <- as.character(as_kable(res, format = "pipe"))
  expect_identical(kab[1L], paste0(
    "Table: Measurement structure, per factor: verbal: essentially ",
    "tau-equivalent (equal loadings); deduction: essentially ",
    "tau-equivalent (equal loadings); Confidence level: 95% (90% for the ",
    "RMSEA interval)"))

  # A column subset loses the model attribute and has no header.
  sub <- res[, c("term", "estimate")]
  expect_false(any(grepl("Measurement structure",
                         capture.output(print(sub)), fixed = TRUE)))
})

test_that("the lines under the table name the version of the fit rows after FIML, and the srmr row is the covariance SRMR", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # three FIML fits and their corrected indices
  # Under FIML the cfi and rmsea rows are lavaan's robust indices, the
  # FIML corrected indices of Zhang and Savalei (2023), which are not
  # computed from the chi_square row, and the line under the table says
  # so. The srmr row is the SRMR of the covariance residuals alone, as
  # sem_fit_measures() reports it (through 1.0.x lavaan's srmr, which
  # averages in the zero residuals of the intercepts).
  hs <- cfa_k_hs_data()
  set.seed(1)
  hs$t7_sentence[sample(301, 90)] <- NA
  f <- cfa_k_hs_factors()
  x <- fa_confirmatory(hs, f, missing = "fiml")
  fit <- fa_confirmatory(hs, f, missing = "fiml", output = "fit")
  fm <- lavaan::fitMeasures(fit)
  val <- function(t) x$estimate[x$term == t]
  expect_equal(val("cfi"), fm[["cfi.robust"]])
  expect_equal(val("rmsea"), fm[["rmsea.robust"]])
  expect_equal(val("chi_square"), fm[["chisq"]])
  expect_equal(val("srmr"), fm[["srmr_bentler_nomean"]])
  expect_gt(abs(val("srmr") - fm[["srmr"]]), 1e-3)
  sf <- sem_fit_measures(fit)
  expect_equal(val("rmsea"), sf$value[sf$term == "rmsea"])
  expect_equal(val("srmr"), sf$value[sf$term == "srmr"])
  lines <- attr(x, "fit_rows")
  expect_length(lines, 1L)
  expect_match(lines[1], paste0(
    "^cfi, rmsea: the FIML corrected indices of Zhang and Savalei \\(2023\\), ",
    "version 6 \\(normal theory\\)"))
  footer <- DMAR:::.dmar_tbl_footer(x)
  expect_identical(unname(footer[names(footer) == "fit_rows"]), lines)
  printed <- paste(capture.output(print(x)), collapse = "\n")
  expect_match(printed, "FIML corrected indices of Zhang and Savalei",
               fixed = TRUE)
  # Complete data under the default FIML say nothing more: with no value
  # missing the corrected indices equal the standard ones.
  expect_null(attr(fa_confirmatory(cfa_k_hs_data(), f), "fit_rows"))
  expect_null(attr(fa_confirmatory(cfa_k_hs_data(), f, missing = "listwise"),
                   "fit_rows"))
})

test_that("the FIML line names version 3 under MLR, and a mean structure leaves the srmr row the covariance SRMR", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two more FIML and mean structure fits
  hs <- cfa_k_hs_data()
  set.seed(1)
  hs$t7_sentence[sample(301, 90)] <- NA
  f <- cfa_k_hs_factors()
  xr <- fa_confirmatory(hs, f, missing = "fiml", estimator = "MLR")
  expect_match(attr(xr, "fit_rows")[1],
               "version 3 \\(no normality assumption\\)")
  xm <- fa_confirmatory(cfa_k_hs_data(), f, missing = "listwise",
                        meanstructure = TRUE)
  expect_null(attr(xm, "fit_rows"))
  fm <- lavaan::fitMeasures(fa_confirmatory(cfa_k_hs_data(), f,
                                            missing = "listwise",
                                            meanstructure = TRUE,
                                            output = "fit"))
  expect_equal(xm$estimate[xm$term == "srmr"], fm[["srmr_bentler_nomean"]])
})

test_that("a requested FIML analysis never silently becomes listwise", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  hs$t7_sentence[1:40] <- NA
  f <- cfa_k_hs_factors()
  for (test in c("mean.var.adjusted", "satorra.bentler")) {
    expect_error(suppressWarnings(fa_confirmatory(hs, f, missing = "ml", test = test)),
                 "lavaan changed missing-data handling to listwise", fixed = TRUE)
  }
})

test_that("one version holds every index row, and the mean and variance adjusted statistic gets the robust indices", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # four fits
  # With the older mean and variance adjusted statistic lavaan computes no
  # robust indices, and through 1.0.1 each row fell to the scaled version
  # (a CFI of .961 here). The robust indices do not depend on which
  # corrected statistic accompanies them (Savalei, 2018), so DMAR computes
  # them, equal to lavaan's for the same model fitted by MLM, and names
  # them under the table.
  hs <- cfa_k_hs_data()
  f <- cfa_k_hs_factors()
  x <- fa_confirmatory(hs, f, test = "mean.var.adjusted")
  fit <- fa_confirmatory(hs, f, test = "mean.var.adjusted", output = "fit")
  mlm <- lavaan::cfa(paste("verbal =~", paste(f$verbal, collapse = " + "),
                           "\ndeduction =~", paste(f$deduction, collapse = " + ")),
                     data = hs, estimator = "MLM")
  val <- function(t) x$estimate[x$term == t]
  rows <- c("cfi", "tli", "nnfi", "rmsea", "rmsea_ci_lower", "rmsea_ci_upper")
  expect_equal(vapply(rows, val, numeric(1L), USE.NAMES = FALSE),
               unname(lavaan::fitMeasures(mlm, c(
                 "cfi.robust", "tli.robust", "nnfi.robust", "rmsea.robust",
                 "rmsea.ci.lower.robust", "rmsea.ci.upper.robust"),
                 fm.args = list(rmsea.ci.level = .90))),
               tolerance = 1e-7)
  expect_equal(val("cfi"), 0.9932502, tolerance = 1e-6)
  expect_equal(unname(lavaan::fitMeasures(fit, "cfi.scaled")), 0.9612618,
               tolerance = 1e-6)
  expect_equal(val("chi_square"), unname(lavaan::fitMeasures(fit, "chisq.scaled")))
  expect_match(attr(x, "fit_rows"), "computed by DMAR", fixed = TRUE)
  # Were the robust indices not to be had, a maximum likelihood fit of
  # continuous items would report its standard indices, which estimate the
  # same population values, in every row, and say so.
  xs <- with_mocked_bindings(fa_confirmatory(hs, f, test = "mean.var.adjusted"),
                             .sem_robust_mva = function(...) numeric(0))
  vs <- function(t) xs$estimate[xs$term == t]
  expect_equal(vapply(rows, vs, numeric(1L), USE.NAMES = FALSE),
               unname(lavaan::fitMeasures(fit, c("cfi", "tli", "nnfi", "rmsea",
                                                 "rmsea.ci.lower",
                                                 "rmsea.ci.upper"),
                 fm.args = list(rmsea.ci.level = .90))))
  expect_match(attr(xs, "fit_rows"),
               "so these rows hold its standard versions", fixed = TRUE)
  # A robust ML fit to the complete cases reports lavaan's robust indices
  # and says nothing more (under FIML, the default, a line names the
  # corrected version 3 indices, which differ from these).
  xr <- fa_confirmatory(hs, f, estimator = "MLR", missing = "listwise")
  fr <- fa_confirmatory(hs, f, estimator = "MLR", missing = "listwise",
                        output = "fit")
  expect_equal(xr$estimate[xr$term == "tli"],
               unname(lavaan::fitMeasures(fr, "tli.robust")))
  expect_null(attr(xr, "fit_rows"))
})

test_that("a fit of ordered items whose robust indices lavaan did not compute names the scaled version it reports", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two ordered fits
  # lavaan leaves the robust indices of this ordered fit NA, so every index
  # row holds its scaled version, which puts the scaled statistic into the
  # standard formulas, and the line under the table says so. The fit_rows
  # attribute is new in 1.1.0; 1.0.1 printed the same numbers with no line.
  hs <- holzinger_swineford
  v <- list(verbal = c("t6_paragraph_comprehension", "t7_sentence",
                       "t9_word_meaning"))
  for (it in v$verbal) {
    hs[[it]] <- as.integer(cut(hs[[it]], stats::quantile(
      hs[[it]], c(0, 0.25, 0.5, 0.75, 1)), include.lowest = TRUE))
  }
  x <- suppressMessages(suppressWarnings(fa_confirmatory(hs, v, ordered = TRUE)))
  fit <- suppressMessages(suppressWarnings(
    fa_confirmatory(hs, v, ordered = TRUE, output = "fit")))
  fm <- lavaan::fitMeasures(fit, fm.args = list(rmsea.ci.level = .90))
  expect_gt(fm[["df"]], 0)
  expect_true(is.na(fm[["cfi.robust"]]))
  rows <- c("cfi", "tli", "nnfi", "rmsea", "rmsea_ci_lower", "rmsea_ci_upper")
  expect_equal(vapply(rows, function(t) x$estimate[x$term == t], numeric(1L),
                      USE.NAMES = FALSE),
               unname(fm[c("cfi.scaled", "tli.scaled", "nnfi.scaled",
                           "rmsea.scaled", "rmsea.ci.lower.scaled",
                           "rmsea.ci.upper.scaled")]))
  expect_match(attr(x, "fit_rows"),
               "so these rows hold its scaled versions, which put the scaled",
               fixed = TRUE)
})

test_that("a model with no degrees of freedom reports exact fit in every version", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # four saturated fits
  # One factor of three items fits exactly: its statistic is zero in every
  # version, so its CFI is 1, and its TLI and RMSEA, 0/0, take lavaan's
  # conventions, 1 and 0. lavaan's robust indices of such a model come
  # out NA by rounding; through 1.0.1 those rows fell back to the scaled
  # or standard values one by one, which the version chosen for every row
  # must not lose.
  hs <- cfa_k_hs_data()
  f <- list(verbal = cfa_k_hs_factors()$verbal)
  rows <- c("cfi", "tli", "nnfi", "rmsea", "rmsea_ci_lower", "rmsea_ci_upper")
  # The rows of such a model do not depend on the version, so no line
  # under the table names one or a fallback from one (with the mean and
  # variance adjusted statistic the first version of the 1.1.0 fix named a
  # fallback to the standard indices).
  fits <- list(ML = fa_confirmatory(hs, f), MLR = fa_confirmatory(hs, f, estimator = "MLR"),
               MVA = fa_confirmatory(hs, f, test = "mean.var.adjusted"))
  for (est in names(fits)) {
    x <- fits[[est]]
    val <- vapply(rows, function(t) x$estimate[x$term == t], numeric(1L),
                  USE.NAMES = FALSE)
    expect_identical(x$estimate[x$term == "df"], 0)
    expect_equal(val, c(1, 1, 1, 0, 0, 0), info = est)
    expect_match(attr(x, "fit_rows"), "saturated", info = est)
  }
  fit <- fa_confirmatory(hs, f, estimator = "MLR", missing = "listwise",
                         output = "fit")
  expect_true(is.na(lavaan::fitMeasures(fit, "cfi.robust")))
  # The same holds under FIML and after two-stage estimation: the rows are
  # those of exact fit, and no line names the corrected indices of Zhang
  # and Savalei (2023), which through the first version of the rule above
  # the FIML line and the two-stage footer still named (lavaan's robust CFI
  # and TLI of the two-stage fit are NA). The two-stage lines on the test,
  # the criteria, and the cases stay.
  set.seed(1)
  hm <- hs
  hm$t7_sentence[sample(nrow(hm), 60)] <- NA
  quiet <- function(x) suppressMessages(suppressWarnings(x))
  incomplete <- list(
    FIML = quiet(fa_confirmatory(hm, f, missing = "fiml")),
    FIML_MLR = quiet(fa_confirmatory(hm, f, missing = "fiml", estimator = "MLR")),
    two_stage = quiet(fa_confirmatory(hm, f, missing = "two_stage")),
    robust_two_stage = quiet(fa_confirmatory(hm, f, missing = "robust_two_stage")))
  for (est in names(incomplete)) {
    x <- incomplete[[est]]
    val <- vapply(rows, function(t) x$estimate[x$term == t], numeric(1L),
                  USE.NAMES = FALSE)
    expect_identical(x$estimate[x$term == "df"], 0, info = est)
    expect_equal(val, c(1, 1, 1, 0, 0, 0), info = est)
    expect_length(attr(x, "fit_rows"), 1L)
    expect_match(attr(x, "fit_rows"), "saturated", info = est)
    printed <- paste(capture.output(print(x)), collapse = " ")
    expect_false(grepl("Zhang and Savalei", printed, fixed = TRUE),
                 info = est)
  }
  for (est in c("two_stage", "robust_two_stage")) {
    footer <- attr(incomplete[[est]], "two_stage")$footer
    expect_false("two_stage_indices" %in% names(footer), info = est)
    expect_true(all(c("two_stage", "two_stage_test", "two_stage_ic") %in%
                      names(footer)), info = est)
  }
  ts_fit <- quiet(fa_confirmatory(hm, f, missing = "two_stage", output = "fit"))
  expect_true(is.na(lavaan::fitMeasures(ts_fit, "cfi.robust")))
  # No degrees of freedom does not make a model fit exactly when an
  # inequality bound binds, which lavaan does not count among them: one
  # factor of three variables correlated .8, .8, and .5 with bounds =
  # "pos.var" has a statistic of 32.8 on 0 df, and its cfi is .925, as in
  # 1.0.1 and lavaan, where the first version of the rule above set it to 1.
  S <- matrix(c(1, .8, .8, .8, 1, .5, .8, .5, 1), 3, 3,
              dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  b <- fa_confirmatory(S = S, N = 200, factors = list(f = c("a", "b", "c")),
             bounds = "pos.var")
  bfit <- fa_confirmatory(S = S, N = 200, factors = list(f = c("a", "b", "c")),
                bounds = "pos.var", output = "fit")
  expect_identical(b$estimate[b$term == "df"], 0)
  expect_gt(b$estimate[b$term == "chi_square"], 30)
  expect_equal(b$estimate[b$term == "cfi"], 0.925, tolerance = 1e-3)
  expect_equal(b$estimate[b$term == "cfi"],
               unname(lavaan::fitMeasures(bfit, "cfi")), tolerance = 1e-12)
  # Its tli, nnfi, rmsea, and RMSEA limits divide by the degrees of
  # freedom and are undefined, so they are NA, as sem_fit_measures()
  # reports them for the same fit, with a line that says why. Through 1.0.1
  # (and the first version of the rule above) they held lavaan's values of
  # exact fit, TLI and NNFI 1 and RMSEA 0 [0, 0], beside a CFI of .925.
  undefined <- c("tli", "nnfi", "rmsea", "rmsea_ci_lower", "rmsea_ci_upper")
  expect_true(all(is.na(b$estimate[match(undefined, b$term)])))
  expect_identical(as.numeric(lavaan::fitMeasures(bfit, c("tli", "rmsea"))),
                   c(1, 0))
  sf <- suppressWarnings(sem_fit_measures(bfit))
  expect_true(all(is.na(sf$value[match(c("tli", "rmsea", "rmsea_ci_lower",
                                         "rmsea_ci_upper"), sf$term)])))
  expect_equal(sf$value[sf$term == "cfi"], b$estimate[b$term == "cfi"],
               tolerance = 1e-12)
  expect_length(attr(b, "fit_rows"), 1L)
  expect_match(attr(b, "fit_rows"),
               paste0("^tli, nnfi, rmsea, and the RMSEA interval: NA, since ",
                      "the model has zero degrees of freedom and a statistic ",
                      "that is not zero"))
  # Under FIML the same holds, and the line naming the corrected indices of
  # Zhang and Savalei (2023) names the cfi alone, where the first version of
  # the rule above attributed the tli, nnfi, and RMSEA interval rows, which
  # held lavaan's values of exact fit, to their correction. Three items
  # whose pairwise correlations have a negative product (-.336, .363, and
  # .558, N = 52) cannot be fitted exactly by one factor, and with no bound
  # the statistic is 4.25 on 0 df.
  set.seed(27)
  n <- sample(c(60, 100, 200), 1)
  r <- sample(c(.05, .15, .3), 1)
  pm <- sample(c(.3, .5, .6), 1)
  R <- matrix(r, 3, 3)
  diag(R) <- 1
  X <- matrix(rnorm(n * 3), n) %*% chol(R) * 10 + 50
  dz <- as.data.frame(X)
  names(dz) <- c("a", "b", "c")
  for (v in names(dz)) dz[[v]][runif(n) < pm] <- NA
  dz <- dz[rowSums(!is.na(dz)) > 0, ]
  z <- quiet(fa_confirmatory(dz, list(f = c("a", "b", "c")), missing = "fiml"))
  expect_identical(z$estimate[z$term == "df"], 0)
  expect_equal(z$estimate[z$term == "chi_square"], 4.25, tolerance = 1e-3)
  expect_equal(z$estimate[z$term == "cfi"], 0.156, tolerance = 1e-2)
  expect_true(all(is.na(z$estimate[match(undefined, z$term)])))
  expect_match(attr(z, "fit_rows")[1L],
               "^cfi: the FIML corrected index of Zhang and Savalei \\(2023\\)")
  expect_false(any(grepl("tli, nnfi, and the RMSEA interval: lavaan's",
                         attr(z, "fit_rows"), fixed = TRUE)))
  expect_match(attr(z, "fit_rows")[2L], "^tli, nnfi, rmsea, and the RMSEA ")
})

test_that("the measurement table keeps its rows when an htmt ratio is undefined", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # The 1.1.0 review: with raw data under FIML in which two items of
  # different factors are never observed together (a two-form planned
  # missing design), htmt() stopped, and fa_confirmatory(output = "measurement") lost
  # every phi, omega, AVE, and H row with it. Through 1.0.1 the table came
  # back with that htmt row NA, without a word.
  skip_if_not_installed("lavaan")
  set.seed(113)
  n <- 300
  f1 <- rnorm(n); f2 <- 0.4 * f1 + sqrt(.84) * rnorm(n)
  f3 <- 0.3 * f1 + sqrt(.91) * rnorm(n)
  d <- data.frame(
    a1 = .8 * f1 + rnorm(n, 0, .6), a2 = .7 * f1 + rnorm(n, 0, .7),
    a3 = .7 * f1 + rnorm(n, 0, .7),
    b1 = .8 * f2 + rnorm(n, 0, .6), b2 = .7 * f2 + rnorm(n, 0, .7),
    b3 = .7 * f2 + rnorm(n, 0, .7),
    c1 = .8 * f3 + rnorm(n, 0, .6), c2 = .7 * f3 + rnorm(n, 0, .7),
    c3 = .7 * f3 + rnorm(n, 0, .7))
  form <- rep(1:2, length.out = n)
  d$a1[form == 1] <- NA
  d$b1[form == 2] <- NA
  bl <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
  # lavaan warns that the pair's coverage is zero; the FIML fit itself is
  # identified through the factor structure.
  w <- testthat::capture_warnings(
    k <- fa_confirmatory(d, factors = bl, missing = "fiml", output = "measurement"))
  expect_identical(k$term, c("phi_A_B", "omega_A", "ave_A", "H_A",
                             "omega_B", "ave_B", "H_B", "htmt_A_B"))
  expect_true(all(is.finite(k$estimate[k$term != "htmt_A_B"])))
  expect_true(is.na(k$estimate[k$term == "htmt_A_B"]))
  expect_true(any(grepl(paste0("^HTMT is undefined for 'A' and 'B'.*",
                               "The htmt rows are NA\\.$"), w)))
  # The model based rows are those of the verbose table of the same fit.
  v <- suppressWarnings(fa_confirmatory(d, factors = bl, missing = "fiml"))
  expect_equal(k$estimate[1:7], v$estimate[match(k$term[1:7], v$term)],
               tolerance = 1e-12)
  # With a third factor, the defined pairs keep their ratios.
  bl3 <- c(bl, list(C = c("c1", "c2", "c3")))
  k3 <- suppressWarnings(fa_confirmatory(d, factors = bl3, missing = "fiml",
                               output = "measurement"))
  ht <- k3$estimate[grepl("^htmt_", k3$term)]
  expect_true(is.na(ht[1L]))
  expect_equal(ht[2:3], suppressWarnings(htmt(d, bl3))$htmt[2:3],
               tolerance = 0)
  # The pair's phi is identified here, through a2, a3, b2, and b3, and its
  # standard error is not zero.
  expect_gt(k$se[k$term == "phi_A_B"], 0)
})

test_that("two factors no case observes together are refused, by name", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # When no item of one factor is ever observed with an item of the other,
  # their correlation does not enter the likelihood and is not identified.
  # Through 1.0.1 (and the help of the 1.1.0 repairs, which said FIML
  # estimates it through the factor structure) the FIML fit reported phi
  # as lavaan's starting value, 0, and every omega, AVE, and H with a
  # standard error of 0 and an interval of no width; two-stage estimation
  # fitted it to an arbitrary covariance with a standard error of 0.
  skip_if_not_installed("lavaan")
  hs <- holzinger_swineford
  sp <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
          "t4_lozenges")
  vb <- c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  sd3 <- c("t10_addition", "t12_counting_groups_of_dots",
           "t13_straight_and_curved_capitals")
  set.seed(3)
  form <- sample(1:2, nrow(hs), replace = TRUE)
  d <- hs
  d[form == 1, sp] <- NA
  d[form == 2, vb] <- NA
  f <- list(spatial = sp, verbal = vb)
  for (m in c("fiml", "two_stage", "robust_two_stage")) {
    expect_error(fa_confirmatory(d, factors = f, missing = m),
                 paste0("No case observes an item of 'spatial' with an item ",
                        "of 'verbal', so the correlation of the two factors ",
                        "is not identified"), fixed = TRUE, info = m)
  }
  # A third factor observed with both does not identify it.
  expect_error(fa_confirmatory(d, factors = c(f, list(speed = sd3)), missing = "fiml"),
               "'spatial' with an item of 'verbal'", fixed = TRUE)
  # Orthogonal factors are not refused: their correlation is fixed at 0,
  # not estimated, and the likelihood separates into the two factors'
  # parts under direct FIML, so each factor's omega, AVE, and H, with their standard errors,
  # are those of its own one factor fit. The first version of the refusal
  # refused them too, saying that the correlation was not identified.
  props <- c("omega_", "ave_", "H_")
  quiet <- function(x) suppressMessages(suppressWarnings(x))
  for (m in "fiml") {
    o <- quiet(fa_confirmatory(d, factors = f, missing = m, correlated_factors = FALSE,
                     output = "measurement"))
    expect_false(any(startsWith(o$term, "phi_")), info = m)
    for (nm in names(f)) {
      one <- quiet(fa_confirmatory(d, factors = f[nm], missing = m,
                         output = "measurement"))
      terms <- paste0(props, nm)
      expect_equal(o$estimate[match(terms, o$term)],
                   one$estimate[match(terms, one$term)], tolerance = 1e-4,
                   info = paste(m, nm))
      expect_equal(o$se[match(terms, o$term)],
                   one$se[match(terms, one$term)], tolerance = 1e-3,
                   info = paste(m, nm))
    }
  }
  for (m in c("two_stage", "robust_two_stage")) {
    expect_error(quiet(fa_confirmatory(d, factors = f, missing = m,
      correlated_factors = FALSE, output = "measurement")),
      "requires unrestricted first-stage moments")
  }
  o <- quiet(fa_confirmatory(d, factors = f, missing = "fiml",
                   correlated_factors = FALSE, output = "measurement"))
  expect_equal(o$estimate[o$term == "omega_spatial"], 0.650, tolerance = 1e-3)
  expect_equal(o$se[o$term == "omega_spatial"], 0.0507, tolerance = 1e-2)
  expect_equal(o$estimate[o$term == "omega_verbal"], 0.845, tolerance = 1e-3)
  # One item of each given on both forms links them, and the fit runs.
  d2 <- hs
  d2[form == 1, sp[-1]] <- NA
  d2[form == 2, vb[-1]] <- NA
  k <- suppressWarnings(fa_confirmatory(d2, factors = f, missing = "fiml",
                              output = "measurement"))
  expect_true(all(k$se[k$term != "htmt_spatial_verbal"] > 0))
  expect_identical(DMAR:::.cfa_k_unlinked_factors(d2, f), character(0))
})

test_that("a fit whose standard errors lavaan cannot compute reports none", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # lavaan gives every defined quantity a standard error of 0 when it
  # cannot invert the information matrix; through 1.0.1 fa_confirmatory() reported
  # omega, AVE, and H with an interval of no width there. Here two items
  # are made uncorrelated with the other factor's items in the sample, so
  # phi is 0 and their loadings are not identified at the estimates.
  skip_if_not_installed("lavaan")
  hs <- holzinger_swineford
  A <- c("t1_visual_perception", "t2_cubes", "t4_lozenges")
  B <- c("t6_paragraph_comprehension", "t7_sentence")
  d <- hs[A]
  X <- cbind(1, as.matrix(d))
  for (v in B) d[[v]] <- as.numeric(hs[[v]] - X %*% qr.solve(X, hs[[v]])) + 10
  f <- list(a = A, b = B)
  # Listwise ML: lavaan cannot invert the expected information.
  w <- testthat::capture_warnings(
    k <- fa_confirmatory(d, factors = f, missing = "listwise",
                         output = "measurement"))
  expect_true(any(grepl("lavaan could not compute the standard errors of ",
                        w, fixed = TRUE)))
  # FIML, the default: lavaan inverts the observed information without
  # complaint (an H standard error of 534), and the singularity check of
  # the information restricted to the constrained directions stops it.
  w2 <- testthat::capture_warnings(
    k2 <- fa_confirmatory(d, factors = f, output = "measurement"))
  expect_true(any(grepl("The information matrix of this fit is singular",
                        w2, fixed = TRUE)))
  expect_true(all(is.na(c(k2$se, k2$ci_lower, k2$ci_upper, k2$p_value))))
  expect_true(all(is.finite(k$estimate[k$term != "htmt_a_b"])))
  expect_true(all(is.na(c(k$se, k$ci_lower, k$ci_upper, k$p_value))))
  v <- suppressWarnings(fa_confirmatory(d, factors = f))
  expect_true(all(is.na(v$se[!is.na(v$estimate) &
                               !v$term %in% c("chi_square", "df",
                                              "p_chi_square", "cfi", "tli",
                                              "nnfi", "rmsea",
                                              "rmsea_ci_lower",
                                              "rmsea_ci_upper",
                                              "rmsea_ci_level", "srmr",
                                              "AIC", "BIC", "H0", "H1")])))
})

test_that("the measurement table of a two-form design keeps the values 1.0.1 reported", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # The values NEWS quotes for the 1.1.0 repair: two three-item factors of
  # 300 cases, a1 missing on one form and b1 on the other. 1.0.1 returned
  # these rows with the htmt row NA and no warning; the first versions of
  # 1.1.0 stopped.
  skip_if_not_installed("lavaan")
  set.seed(113)
  n <- 300
  f1 <- rnorm(n); f2 <- 0.4 * f1 + sqrt(.84) * rnorm(n)
  d <- data.frame(a1 = .8 * f1 + rnorm(n, 0, .6),
                  a2 = .7 * f1 + rnorm(n, 0, .7),
                  a3 = .7 * f1 + rnorm(n, 0, .7),
                  b1 = .8 * f2 + rnorm(n, 0, .6),
                  b2 = .7 * f2 + rnorm(n, 0, .7),
                  b3 = .7 * f2 + rnorm(n, 0, .7))
  form <- rep(1:2, length.out = n)
  d$a1[form == 1] <- NA
  d$b1[form == 2] <- NA
  bl <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
  w <- testthat::capture_warnings(
    k <- fa_confirmatory(d, factors = bl, missing = "fiml", output = "measurement"))
  expect_true(any(grepl("^HTMT is undefined for 'A' and 'B'", w)))
  est <- setNames(k$estimate, k$term)
  expect_equal(unname(est[c("phi_A_B", "omega_A", "ave_A", "H_A",
                            "omega_B", "ave_B", "H_B")]),
               c(0.3801163, 0.7927398, 0.5578866, 0.7983507, 0.7838866,
                 0.5494796, 0.7919046), tolerance = 1e-5)
  expect_equal(c(k$ci_lower[k$term == "phi_A_B"],
                 k$ci_upper[k$term == "phi_A_B"]),
               c(0.2392005, 0.5210321), tolerance = 1e-5)
  expect_true(is.na(est[["htmt_A_B"]]))
})


test_that("the help describes released history rather than an interim implementation", {
  path <- file.path("..", "..", "man", "fa_confirmatory.Rd")
  text <- if (file.exists(path)) {
    paste(readLines(path, warn = FALSE), collapse = " ")
  } else {
    rd <- tools::Rd_db("DMAR")[["fa_confirmatory.Rd"]]
    paste(unlist(rd), collapse = " ")
  }
  text <- gsub("[[:space:]]+", " ", text)
  expect_false(grepl("Through the first versions", text, fixed = TRUE))
})

test_that("cfa_k() is fa_confirmatory() under its 1.0.0 name", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  a <- cfa_k(cfa_k_hs_data(), cfa_k_hs_factors())
  b <- fa_confirmatory(cfa_k_hs_data(), cfa_k_hs_factors())
  # Attributes may carry a lavaan object with wall-clock timings, so the
  # table is compared bare and through its tidiers.
  bare <- function(x) {
    attributes(x) <- attributes(x)[c("names", "row.names", "class")]
    x
  }
  expect_identical(bare(a), bare(b))
  expect_identical(tidy(a), tidy(b))
  expect_identical(glance(a), glance(b))
  expect_identical(formals(cfa_k), formals(fa_confirmatory))
})

test_that("omega, AVE, and H rows carry logit-scale limits, as reliability_H() and fa_exploratory() form H's", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  d <- holzinger_swineford[items]
  factors <- list(visual = items[1:3], verbal = items[4:6])
  res <- fa_confirmatory(d, factors = factors, missing = "fiml")
  fit <- fa_confirmatory(d, factors = factors, missing = "fiml", output = "fit")
  crit <- qnorm(0.975)
  for (f in names(factors)) {
    h <- res[res$term == paste0("H_", f), ]
    expect_true(h$ci_lower > 0 && h$ci_upper < 1)
    expect_true(h$ci_lower < h$estimate && h$estimate < h$ci_upper)
    lim <- DMAR:::.logistic_ci(h$estimate, h$se, crit)
    expect_equal(h$ci_lower, lim$lower)
    expect_equal(h$ci_upper, lim$upper)
    # The same fit through reliability_H(): the same estimate, standard
    # error, and limits.
    rh <- reliability_H(fit = fit, factor = f)
    expect_equal(h$estimate, rh$value[rh$term == "estimate"], tolerance = 1e-6)
    expect_equal(h$se, rh$value[rh$term == "se"], tolerance = 1e-5)
    expect_equal(h$ci_lower, rh$value[rh$term == "ci_lower"], tolerance = 1e-5)
    expect_equal(h$ci_upper, rh$value[rh$term == "ci_upper"], tolerance = 1e-5)
    # Omega and AVE take the same form from their own standard errors.
    for (q in c("omega_", "ave_")) {
      row <- res[res$term == paste0(q, f), ]
      lim <- DMAR:::.logistic_ci(row$estimate, row$se, crit)
      expect_equal(c(row$ci_lower, row$ci_upper), c(lim$lower, lim$upper))
      expect_true(row$ci_lower > 0 && row$ci_upper < 1)
    }
  }
  # A finite wald_df puts the t quantile on the logit scale.
  t_res <- fa_confirmatory(d, factors = factors, missing = "fiml", wald_df = 50)
  h <- t_res[t_res$term == "H_visual", ]
  lim <- DMAR:::.logistic_ci(h$estimate, h$se, qt(0.975, 50))
  expect_equal(h$ci_lower, lim$lower)
  expect_equal(h$ci_upper, lim$upper)
  # The t quantile widens the interval on the logit scale as well.
  expect_gt(res$ci_lower[res$term == "H_visual"] - h$ci_lower, 0)
  expect_lt(res$ci_upper[res$term == "H_visual"] - h$ci_upper, 0)
  # The standardized output's H rows take the same form.
  std <- fa_confirmatory(d, factors = factors, missing = "fiml", output = "standardized")
  hs <- std[std$op == ":=" & std$lhs == "H_visual", ]
  lim <- DMAR:::.logistic_ci(hs$est.std, hs$se, crit)
  expect_equal(hs$ci.lower, lim$lower)
  expect_equal(hs$ci.upper, lim$upper)
  expect_true(all(is.na(c(hs$ci.lower, hs$ci.upper))) || (hs$ci.lower > 0 && hs$ci.upper < 1))
  # With one factor the EFA's H interval is the same interval.
  one <- fa_exploratory(d, 1, seed = 113)
  cfa1 <- fa_confirmatory(d, factors = list(f1 = items), missing = "fiml")
  expect_equal(one$ci_lower[one$op == "H"], cfa1$ci_lower[cfa1$term == "H_f1"], tolerance = 1e-5)
  expect_equal(one$ci_upper[one$op == "H"], cfa1$ci_upper[cfa1$term == "H_f1"], tolerance = 1e-5)
})


test_that("the default 'missing' is FIML, resolved by rule for the paths FIML cannot serve", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  f <- cfa_k_hs_factors()
  d <- hs
  d$t7_sentence[1:40] <- NA
  quiet <- function(x) suppressMessages(suppressWarnings(x))
  # Complete data: the default fits by FIML, with the estimates of the
  # listwise fit.
  a <- fa_confirmatory(hs, f)
  b <- fa_confirmatory(hs, f, missing = "listwise")
  expect_identical(attr(a, "missing"), "fiml")
  ld <- grepl("^lambda_", a$term)
  expect_equal(a$estimate[ld], b$estimate[match(a$term[ld], b$term)],
               tolerance = 1e-6)
  # Incomplete data: every case with an observed item is used, without a
  # warning.
  expect_no_warning(fa_confirmatory(d, f))
  expect_identical(lavaan::lavInspect(fa_confirmatory(d, f, output = "fit"),
                                      "nobs"), nrow(d))
  expect_silent(fa_confirmatory(d, f, missing = "fiml", output = "measurement"))
  # The wrappers forward an unsupplied 'missing' unsupplied.
  expect_no_warning(cfa_1(d, items = f$verbal))
  expect_silent(cfa_1(d, items = f$verbal, missing = "fiml",
                      output = "measurement"))
  # se = "robust.sem": fitted as it stands with complete data, refused
  # with incomplete data unless FIML or listwise is asked for by name.
  r <- fa_confirmatory(hs, f, se = "robust.sem")
  expect_identical(attr(r, "missing"), "listwise")
  expect_error(fa_confirmatory(d, f, se = "robust.sem"),
               "robust.huber.white")
  expect_error(fa_confirmatory(d, f, se = "robust.sem", missing = "fiml"),
               "incompatible with FIML")
  expect_identical(attr(fa_confirmatory(d, f, se = "robust.sem",
                                        missing = "listwise"), "missing"),
                   "listwise")
  # Ordered items: as they stand with complete data, refused with missing
  # values, naming pairwise and listwise deletion.
  o <- hs[unlist(f)]
  o[] <- lapply(o, function(v) cut(v, 4, labels = FALSE))
  expect_identical(attr(quiet(fa_confirmatory(o, f, ordered = TRUE)),
                        "missing"), "listwise")
  o$t7_sentence[1:40] <- NA
  expect_error(quiet(fa_confirmatory(o, f, ordered = TRUE)),
               "missing = \"pairwise\" for ordered items")
  # A least squares estimator likewise.
  expect_error(fa_confirmatory(d, f, estimator = "GLS"),
               "estimator = \"GLS\"")
})

test_that("the default of se matches the estimator, and se = \"standard\" gives the normal theory standard errors under MLR", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # DMAR 1.0.0's cfa_k(), cfa_1(), and cfa_2() passed se = "standard" to
  # lavaan under every estimator, so their MLR standard errors were the
  # normal theory ones; the default now gives Huber-White ones under MLR.
  hs <- cfa_k_hs_data()
  f <- cfa_k_hs_factors()
  expect_no_warning(x <- fa_confirmatory(hs, f, estimator = "MLR",
                                         missing = "listwise"))
  # No warning with se left at its default or given, under MLR or ML, or
  # from a wrapper.
  expect_silent(fa_confirmatory(hs, f, estimator = "MLR",
                                missing = "listwise",
                                output = "measurement"))
  expect_silent(fa_confirmatory(hs, f, estimator = "MLR",
                                missing = "listwise", se = "standard",
                                output = "measurement"))
  expect_silent(fa_confirmatory(hs, f, missing = "listwise",
                                output = "measurement"))
  expect_silent(cfa_1(hs, items = f$verbal, estimator = "MLR",
                      missing = "listwise", output = "measurement"))
  # se = "standard" gives lavaan's standard errors under MLR, which DMAR
  # 1.0.0 reported (checked against DMAR 1.0.0's cfa_k() on 2026-11-05:
  # identical to the last digit).
  y <- fa_confirmatory(hs, f, estimator = "MLR", missing = "listwise",
                       se = "standard")
  fit <- suppressWarnings(lavaan::cfa(
    paste0("verbal =~ ", paste(f$verbal, collapse = " + "), "\n",
           "deduction =~ ", paste(f$deduction, collapse = " + ")),
    data = hs, estimator = "MLR", se = "standard", std.lv = TRUE))
  pe <- lavaan::parameterEstimates(fit)
  expect_equal(y$se[startsWith(y$term, "lambda_")],
               pe$se[pe$op == "=~"], tolerance = 1e-6)
  expect_false(isTRUE(all.equal(x$se[startsWith(x$term, "lambda_")],
                                y$se[startsWith(y$term, "lambda_")])))
})

test_that("the singularity check leaves a model identified by its equality constraints alone", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  f <- list(a = c("t1_visual_perception", "t2_cubes"),
            b = c("t6_paragraph_comprehension", "t7_sentence",
                  "t9_word_meaning"))
  expect_silent(x <- fa_confirmatory(hs, f, equal_loading = c(a = TRUE, b = FALSE),
                                     output = "measurement"))
  expect_true(all(is.finite(x$se[x$term %in% c("omega_a", "omega_b")])))
  # The expected information of a listwise fit and of a covariance matrix
  # is tested as well, and this model passes there too.
  expect_silent(y <- fa_confirmatory(hs, f, equal_loading = c(a = TRUE, b = FALSE),
                                     missing = "listwise",
                                     output = "measurement"))
  expect_true(all(is.finite(y$se[y$term %in% c("omega_a", "omega_b")])))
  items <- unlist(f, use.names = FALSE)
  expect_silent(s <- fa_confirmatory(S = stats::cov(hs[items]), N = nrow(hs),
                                     factors = f,
                                     equal_loading = c(a = TRUE, b = FALSE),
                                     output = "measurement"))
  expect_true(all(is.finite(s$se[s$term %in% c("omega_a", "omega_b")])))
})

test_that("a singular expected information that lavaan inverts withholds the standard errors of every table and of reliability_H()", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # eight more fits
  # The fixture of the test above in which lavaan cannot compute the
  # standard errors, with equal loadings on the spatial factor: lavaan then
  # inverts the singular expected information of the listwise and the
  # covariance matrix fits, with only its own warning that the result
  # does not appear to be positive definite. Through 1.1.0 development
  # fa_confirmatory() tested the observed information alone, and its
  # tables reported H_b = .853 with a standard error of .028 and the
  # interval [.790, .900], while reliability_H() of the same fit withheld
  # its standard error.
  hs <- holzinger_swineford
  A <- c("t1_visual_perception", "t2_cubes", "t4_lozenges")
  B <- c("t6_paragraph_comprehension", "t7_sentence")
  d <- hs[A]
  X <- cbind(1, as.matrix(d))
  for (v in B) d[[v]] <- as.numeric(hs[[v]] - X %*% qr.solve(X, hs[[v]])) + 10
  f <- list(a = A, b = B)
  eq <- c(a = TRUE, b = FALSE)
  qualified <- function(expr) {
    seen <- new.env(parent = emptyenv())
    seen$w <- list()
    value <- withCallingHandlers(expr, warning = function(w) {
      if (inherits(w, "dmar_result_qualified")) seen$w <- c(seen$w, list(w))
      invokeRestart("muffleWarning")
    })
    list(value = value, w = seen$w)
  }
  calls <- list(
    listwise = function(out) fa_confirmatory(d, factors = f, equal_loading = eq,
                                             missing = "listwise", output = out),
    covariance = function(out) fa_confirmatory(S = stats::cov(d), N = nrow(d),
                                               factors = f, equal_loading = eq,
                                               output = out))
  for (k in names(calls)) {
    fit <- suppressWarnings(calls[[k]]("fit"))
    expect_identical(lavaan::lavInspect(fit, "options")$information[1L],
                     "expected")
    # Another code path: the smallest eigenvalue of the expected
    # information scaled to a unit diagonal and restricted to the
    # directions the equality constraints leave free.
    I <- unclass(lavaan::lavInspect(fit, "information"))
    s <- 1 / sqrt(abs(diag(I)))
    J <- lavaan::lavInspect(fit, "con.jac")
    ce <- attr(J, "ceq.idx")
    N0 <- MASS::Null(t(unclass(J)[ce, , drop = FALSE] *
                          rep(s, each = length(ce))))
    R <- crossprod(N0, (I * outer(s, s)) %*% N0)
    R <- R / sqrt(outer(diag(R), diag(R)))
    expect_lt(min(eigen((R + t(R)) / 2, symmetric = TRUE)$values), 1e-6)
    m <- qualified(calls[[k]]("measurement"))
    expect_length(m$w, 1L)
    expect_true(all(is.na(m$value$se)), label = k)
    expect_true(all(is.na(c(m$value$ci_lower, m$value$ci_upper))), label = k)
    expect_equal(m$value$estimate[m$value$term == "H_b"], 0.85291,
                 tolerance = 1e-4)
    z <- qualified(calls[[k]]("standardized"))
    expect_length(z$w, 1L)
    expect_true(all(is.na(z$value$se)), label = k)
    h <- qualified(reliability_H(fit = fit, factor = "b"))
    expect_true(all(is.na(h$value$value[h$value$term %in%
                                          c("se", "ci_lower", "ci_upper")])))
    expect_identical(conditionMessage(h$w[[1L]]),
                     conditionMessage(m$w[[1L]]))
  }
})

test_that("the RMSEA interval is at rmsea_conf_level, apart from the parameter intervals", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # seven fits; the default level is checked on CRAN above
  # DMAR 1.0.0's cfa_k() reported the RMSEA interval at lavaan's 0.90
  # whatever conf_level was; the default reproduces it, and the level can
  # now be chosen without moving the parameter intervals.
  hs <- cfa_k_hs_data()
  f <- cfa_k_hs_factors()
  fit <- fa_confirmatory(hs, f, output = "fit")
  limits <- function(x) {
    x$estimate[match(c("rmsea_ci_lower", "rmsea_ci_upper", "rmsea_ci_level"),
                     x$term)]
  }
  lavaan_limits <- function(level) {
    unname(lavaan::fitMeasures(fit, c("rmsea.ci.lower", "rmsea.ci.upper"),
                               fm.args = list(rmsea.ci.level = level)))
  }
  params <- function(x) {
    p <- !x$term %in% c("chi_square", "df", "p_chi_square", "cfi", "tli",
                        "nnfi", "rmsea", "rmsea_ci_lower", "rmsea_ci_upper",
                        "rmsea_ci_level", "srmr", "AIC", "BIC", "H0", "H1")
    x[p, c("ci_lower", "ci_upper")]
  }
  expect_equal(eval(formals(fa_confirmatory)$rmsea_conf_level), 0.90)
  expect_identical(formals(cfa_k), formals(fa_confirmatory))
  x <- fa_confirmatory(hs, f)
  # The default matches lavaan's own (0.90) limits, as cfa_k() did in 1.0.0.
  expect_equal(limits(x)[1:2], lavaan_limits(0.90), tolerance = 1e-10)
  expect_equal(limits(cfa_k(hs, f)), limits(x))
  expect_equal(limits(x)[3], 0.90)
  expect_equal(attr(x, "conf_level"), 0.95)
  expect_equal(attr(x, "rmsea_conf_level"), 0.90)
  # rmsea_conf_level moves the RMSEA interval and nothing else.
  y <- fa_confirmatory(hs, f, rmsea_conf_level = 0.95)
  expect_equal(limits(y)[1:2], lavaan_limits(0.95), tolerance = 1e-10)
  expect_equal(limits(y)[3], 0.95)
  expect_identical(params(y), params(x))
  # conf_level moves the parameter intervals and not the RMSEA interval.
  z <- fa_confirmatory(hs, f, conf_level = 0.99)
  expect_identical(limits(z), limits(x))
  expect_false(isTRUE(all.equal(params(z), params(x))))
  # The line under the table states both levels when they differ, and one
  # when they agree.
  out <- utils::capture.output(print(x))
  expect_true("Confidence level: 95% (90% for the RMSEA interval)" %in% out)
  same <- utils::capture.output(print(fa_confirmatory(hs, f, conf_level = 0.90)))
  expect_true("Confidence level: 90%" %in% same)
  expect_false(any(grepl("for the RMSEA interval", same, fixed = TRUE)))
  expect_false(any(grepl("for the RMSEA interval", utils::capture.output(
    print(fa_confirmatory(hs, f, output = "measurement"))), fixed = TRUE)))
  # glance(), tidy(), results_sentence(), and dmar_metadata() keep the two
  # levels apart.
  expect_identical(generics::glance(y)$rmsea_ci_level, 0.95)
  expect_identical(unique(generics::tidy(y)$conf_level), 0.95)
  expect_identical(results_sentence(x, term = "rmsea", digits = 3),
                   sprintf("rmsea = %.3f, 90%% CI [%.3f, %.3f]",
                           x$estimate[x$term == "rmsea"], limits(x)[1],
                           limits(x)[2]))
  md <- dmar_metadata(y)
  expect_identical(c(md$conf_level, md$rmsea_conf_level), c(0.95, 0.95))
  expect_identical(dmar_metadata(z)$rmsea_conf_level, 0.90)
  expect_error(fa_confirmatory(hs, f, rmsea_conf_level = 1),
               "'rmsea_conf_level' must be a single number in (0, 1).",
               fixed = TRUE)
  expect_error(fa_confirmatory(hs, f, rmsea_conf_level = c(.9, .95)),
               "'rmsea_conf_level'")
})

test_that("cfa_1() and cfa_2() forward rmsea_conf_level apart from conf_level", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # six more fits
  hs <- cfa_k_hs_data()
  f <- cfa_k_hs_factors()
  four <- c(f$verbal, "t5_general_information")
  rows <- c("rmsea_ci_lower", "rmsea_ci_upper", "rmsea_ci_level")
  get <- function(x, t) x$estimate[match(t, x$term)]
  one <- cfa_1(hs, items = four)
  one_95 <- cfa_1(hs, items = four, rmsea_conf_level = 0.95)
  one_99 <- cfa_1(hs, items = four, conf_level = 0.99)
  ref <- fa_confirmatory(hs, list(f1 = four), rmsea_conf_level = 0.95)
  expect_equal(get(one, "rmsea_ci_level"), 0.90)
  expect_equal(get(one_95, rows), get(ref, rows))
  expect_identical(get(one_99, rows), get(one, rows))
  expect_identical(one_95$ci_lower[1:4], one$ci_lower[1:4])
  two <- cfa_2(hs, factor_1 = f$verbal, factor_2 = f$deduction)
  two_95 <- cfa_2(hs, factor_1 = f$verbal, factor_2 = f$deduction,
                  rmsea_conf_level = 0.95)
  two_99 <- cfa_2(hs, factor_1 = f$verbal, factor_2 = f$deduction,
                  conf_level = 0.99)
  expect_equal(get(two, "rmsea_ci_level"), 0.90)
  expect_equal(get(two_95, "rmsea_ci_level"), 0.95)
  expect_gt(get(two_95, "rmsea_ci_upper"), get(two, "rmsea_ci_upper"))
  expect_identical(get(two_99, rows), get(two, rows))
  expect_identical(two_95$ci_lower[1:6], two$ci_lower[1:6])
})

# The standardized defined rows under equality constraints ---------------
#
# An oracle on another code path: the standardized indices of every factor
# as functions of the free parameter vector, by matrix algebra on the
# lavInspect() matrices (Sigma = Lambda Phi Lambda' + Theta, each item
# standardized by its model implied variance), with the delta method
# standard error from a central difference Jacobian and lavaan's vcov().
cfa_std_oracle <- function(fit) {
  fr <- lavaan::lavInspect(fit, "free")
  es <- lavaan::lavInspect(fit, "est")
  at <- function(x) {
    L <- es$lambda; P <- es$psi; Th <- es$theta
    L[fr$lambda > 0] <- x[fr$lambda[fr$lambda > 0]]
    P[fr$psi > 0] <- x[fr$psi[fr$psi > 0]]
    Th[fr$theta > 0] <- x[fr$theta[fr$theta > 0]]
    P[upper.tri(P)] <- t(P)[upper.tri(P)]
    Th[upper.tri(Th)] <- t(Th)[upper.tri(Th)]
    Sig <- L %*% P %*% t(L) + Th
    d <- 1 / sqrt(diag(Sig))
    Ls <- d * L %*% diag(sqrt(diag(P)), ncol(P))
    e_all <- diag(Th) * d^2
    out <- numeric()
    for (k in seq_len(ncol(L))) {
      on <- which(L[, k] != 0)
      l <- Ls[on, k]
      e <- e_all[on]
      out[paste0(c("loading_sum_", "error_sum_", "omega_", "ave_", "H_"),
                 colnames(L)[k])] <-
        c(sum(l), sum(e), sum(l)^2 / (sum(l)^2 + sum(e)),
          mean(l^2 / (l^2 + e)), 1 / (1 + 1 / sum(l^2 / e)))
    }
    out
  }
  x0 <- lavaan::coef(fit)
  g0 <- at(x0)
  J <- vapply(seq_along(x0), function(k) {
    h <- 1e-5 * max(1, abs(x0[k]))
    xp <- x0
    xm <- x0
    xp[k] <- xp[k] + h
    xm[k] <- xm[k] - h
    (at(xp) - at(xm)) / (2 * h)
  }, numeric(length(g0)))
  list(est = g0, se = sqrt(diag(J %*% lavaan::vcov(fit) %*% t(J))))
}

test_that("standardized omega, AVE, and H use each item's own standardized values under every equality constraint", {
  skip_on_cran()  # numerical oracle check, run in the full local suite
  skip_if_not_installed("lavaan")
  # Four items with equal unstandardized loadings .7 and error variances
  # .1, .2, 1, and 2, so their standardized loadings differ. lavaan
  # evaluates a defined parameter at the first item's standardized value
  # of a shared label: through 1.1.0 development the standardized AVE,
  # H, and omega under equal_loading were .658, .909, and .873, and under
  # equal_error .316, .664, and .652, with lavaan's warning that they
  # were ambiguous.
  nms <- paste0("y", 1:4)
  S <- matrix(.49, 4, 4) + diag(c(.1, .2, 1, 2))
  dimnames(S) <- list(nms, nms)
  expected <- list(
    "FALSE.FALSE" = c(ave_f = .5165749, H_f = .8899285, omega_f = .7988381),
    "FALSE.TRUE" = c(ave_f = .3816192, H_f = .7709822, omega_f = .7007462),
    "TRUE.FALSE" = c(ave_f = .5165749, H_f = .8899285, omega_f = .7988381),
    "TRUE.TRUE" = c(ave_f = .3726236, H_f = .7037702, omega_f = .7037702))
  for (el in c(FALSE, TRUE)) for (ee in c(FALSE, TRUE)) {
    expect_no_warning(z <- fa_confirmatory(
      S = S, N = 500, factors = list(f = nms), equal_loading = el,
      equal_error = ee, output = "standardized"))
    fit <- fa_confirmatory(S = S, N = 500, factors = list(f = nms),
                           equal_loading = el, equal_error = ee,
                           output = "fit")
    # Direct algebra on the items' own rows of the same table.
    l <- z$est.std[z$op == "=~"]
    e <- z$est.std[z$op == "~~" & z$lhs %in% nms & z$lhs == z$rhs]
    direct <- c(loading_sum_f = sum(l), error_sum_f = sum(e),
                omega_f = sum(l)^2 / (sum(l)^2 + sum(e)),
                ave_f = mean(l^2), H_f = sum(l^2 / e) / (1 + sum(l^2 / e)))
    def <- z[z$op == ":=", ]
    got <- stats::setNames(def$est.std, def$lhs)[names(direct)]
    expect_equal(got, direct, tolerance = 1e-7)
    key <- paste(el, ee, sep = ".")
    expect_equal(got[names(expected[[key]])], expected[[key]],
                 tolerance = 1e-6)
    # The delta method standard errors match the oracle's.
    o <- cfa_std_oracle(fit)
    se <- stats::setNames(def$se, def$lhs)[names(direct)]
    expect_equal(se, o$se[names(direct)], tolerance = 1e-6)
    expect_equal(got, o$est[names(direct)], tolerance = 1e-8)
    # AVE and H do not depend on the metric of the items, so they are the
    # rows of the unstandardized table, standard errors included.
    m <- fa_confirmatory(S = S, N = 500, factors = list(f = nms),
                         equal_loading = el, equal_error = ee,
                         output = "measurement")
    for (t in c("ave_f", "H_f")) {
      expect_equal(unname(got[t]), m$estimate[m$term == t], tolerance = 1e-8)
      expect_equal(unname(se[t]), m$se[m$term == t], tolerance = 1e-6)
    }
    # The z and p-value are formed as lavaan forms them for every other
    # row.
    expect_equal(def$z, def$est.std / def$se)
    expect_equal(def$pvalue, 2 * (1 - stats::pnorm(abs(def$z))))
    # Omega and AVE and H take their logit-scale limits from these errors.
    for (t in c("omega_f", "ave_f", "H_f")) {
      r <- def[def$lhs == t, ]
      lim <- DMAR:::.logistic_ci(r$est.std, r$se, qnorm(0.975))
      expect_equal(c(r$ci.lower, r$ci.upper), c(lim$lower, lim$upper))
    }
  }
  # The raw-score omega of the unstandardized table stays apart from the
  # omega of the standardized items under equal loadings (.704 against
  # .799), and the defined rows of the unconstrained model are lavaan's.
  m <- fa_confirmatory(S = S, N = 500, factors = list(f = nms),
                       equal_loading = TRUE, output = "measurement")
  expect_equal(m$estimate[m$term == "omega_f"], .7037702, tolerance = 1e-6)
})

test_that("the recomputed standardized rows are lavaan's where lavaan's are unambiguous", {
  skip_if_not_installed("lavaan")
  hs <- cfa_k_hs_data()
  f <- cfa_k_hs_factors()
  fit <- fa_confirmatory(hs, f, output = "fit")
  std <- lavaan::standardizedSolution(fit)
  re <- DMAR:::.fa_cfa_std_defined(fit, std, f, 0.95)
  def <- std$op == ":="
  for (col in c("est.std", "se", "z", "pvalue", "ci.lower", "ci.upper")) {
    expect_equal(re[[col]][def], std[[col]][def], tolerance = 1e-7)
  }
  expect_identical(re[!def, ], std[!def, ])
})

test_that("the standardized rows read lavaan's covariance matrix in either ceq.simple layout", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # four more fits
  hs <- cfa_k_hs_data()
  f <- list(a = c("t1_visual_perception", "t2_cubes"),
            b = c("t6_paragraph_comprehension", "t7_sentence",
                  "t9_word_meaning"))
  # Two items with equal loadings identify their factor. Under
  # ceq.simple = TRUE the singularity check took the information of the
  # unconstrained parameters for that of the model, and through 1.1.0
  # development every standard error of this fit was withheld.
  run <- function(ceq, out) {
    fa_confirmatory(hs, f, equal_loading = c(a = TRUE, b = FALSE),
                    equal_error = c(a = FALSE, b = TRUE), output = out,
                    ceq.simple = ceq)
  }
  expect_no_warning(s1 <- run(TRUE, "standardized"))
  s0 <- run(FALSE, "standardized")
  def <- s0$op == ":="
  expect_equal(s1$est.std[def], s0$est.std[def], tolerance = 1e-5)
  expect_equal(s1$se[def], s0$se[def], tolerance = 1e-4)
  # The oracle indexes the parameters by free number, the layout of
  # lavaan's vcov() under ceq.simple = FALSE.
  o <- cfa_std_oracle(run(FALSE, "fit"))
  expect_equal(stats::setNames(s0$se[def], s0$lhs[def])[names(o$se)], o$se,
               tolerance = 1e-6)
  expect_no_warning(m1 <- run(TRUE, "measurement"))
  expect_true(all(is.finite(m1$se[m1$term != "htmt_a_b"])))
})

test_that("bootstrap standard errors and limits of the standardized rows come from the draws", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two bootstraps of 40 draws
  hs <- cfa_k_hs_data()
  f <- cfa_k_hs_factors()
  set.seed(113)
  fit <- fa_confirmatory(hs, f, se = "bootstrap", bootstrap = 40,
                         equal_loading = c(verbal = TRUE, deduction = FALSE),
                         output = "fit")
  std <- withCallingHandlers(lavaan::standardizedSolution(fit),
                             warning = function(w) {
                               if (DMAR:::.fa_cfa_std_ambiguous(w)) {
                                 invokeRestart("muffleWarning")
                               }
                             })
  re <- DMAR:::.fa_cfa_std_defined(fit, std, f["verbal"], 0.95)
  # The oracle: the standardized indices of the verbal factor at every
  # draw, by matrix algebra, and the standard deviation lavaan forms
  # (divisor the number of draws).
  draws <- lavaan::lavInspect(fit, "coef.boot")
  fr <- lavaan::lavInspect(fit, "free")
  es <- lavaan::lavInspect(fit, "est")
  vals <- t(apply(unclass(draws), 1L, function(x) {
    L <- es$lambda; Th <- es$theta
    L[fr$lambda > 0] <- x[fr$lambda[fr$lambda > 0]]
    Th[fr$theta > 0] <- x[fr$theta[fr$theta > 0]]
    on <- which(L[, "verbal"] != 0)
    v <- L[on, "verbal"]^2 + diag(Th)[on]
    l <- L[on, "verbal"] / sqrt(v)
    e <- diag(Th)[on] / v
    c(sum(l), sum(e), sum(l)^2 / (sum(l)^2 + sum(e)), mean(l^2),
      sum(l^2 / e) / (1 + sum(l^2 / e)))
  }))
  n <- nrow(vals)
  terms <- paste0(c("loading_sum_", "error_sum_", "omega_", "ave_", "H_"),
                  "verbal")
  got <- re[match(terms, re$lhs), ]
  expect_equal(got$se, unname(apply(vals, 2, sd)) * sqrt((n - 1) / n),
               tolerance = 1e-8)
  # The deduction factor shares no label, so lavaan's rows are kept, and
  # recomputed they are the same, standard errors and percentile limits.
  again <- DMAR:::.fa_cfa_std_defined(fit, std, f["deduction"], 0.95)
  d_rows <- std$op == ":=" & grepl("deduction$", std$lhs)
  for (col in c("est.std", "se", "ci.lower", "ci.upper")) {
    expect_equal(again[[col]][d_rows], std[[col]][d_rows], tolerance = 1e-8)
  }
})

test_that("ordered items with equal loadings give standardized rows from each item's latent response", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a WLSMV fit
  hs <- cfa_k_hs_data()
  items <- c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning",
             "t5_general_information")
  ord <- as.data.frame(lapply(hs[items], function(x) {
    cut(x, stats::quantile(x, c(0, .3, .7, 1)), include.lowest = TRUE,
        labels = FALSE)
  }))
  z <- suppressMessages(fa_confirmatory(ord, list(v = items), ordered = TRUE,
                                        equal_loading = TRUE,
                                        output = "standardized"))
  l <- z$est.std[z$op == "=~"]
  e <- z$est.std[z$op == "~~" & z$lhs %in% items & z$lhs == z$rhs]
  expect_gt(diff(range(l)), 0.01)
  def <- z[z$op == ":=", ]
  get <- function(t) def$est.std[def$lhs == t]
  expect_equal(get("ave_v"), mean(l^2), tolerance = 1e-7)
  expect_equal(get("omega_v"), sum(l)^2 / (sum(l)^2 + sum(e)), tolerance = 1e-7)
  expect_equal(get("H_v"), sum(l^2 / e) / (1 + sum(l^2 / e)), tolerance = 1e-7)
})

test_that("the standardized solution without standard errors gives its estimates and NA columns", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # standardizedSolution() has no standard error, test, or interval
  # columns under se = "none", and through 1.1.0 development the
  # logit-scale limits stopped such a call with the R error "missing
  # value where TRUE/FALSE needed", with or without equality constraints.
  hs <- cfa_k_hs_data()
  f <- list(vis = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
            verb = c("t6_paragraph_comprehension", "t7_sentence",
                     "t9_word_meaning"))
  for (el in c(FALSE, TRUE)) {
    expect_no_warning(n <- fa_confirmatory(hs, f, se = "none",
                                           equal_loading = el,
                                           output = "standardized"))
    z <- fa_confirmatory(hs, f, equal_loading = el, output = "standardized")
    expect_equal(n$est.std, z$est.std, tolerance = 1e-10)
    expect_true(all(is.na(unlist(n[c("se", "z", "pvalue", "ci.lower",
                                     "ci.upper")]))))
    # The defined rows of equal loadings are the standardized items' own.
    l <- n$est.std[n$op == "=~" & n$lhs == "vis"]
    e <- n$est.std[n$op == "~~" & n$lhs == n$rhs & n$lhs %in% f$vis]
    expect_equal(n$est.std[n$lhs == "omega_vis"],
                 sum(l)^2 / (sum(l)^2 + sum(e)), tolerance = 1e-10)
  }
})

test_that("with auxiliary variables the srmr row is over the items", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939
  set.seed(1)
  for (j in paste0("x", 1:9)) hs[sample(nrow(hs), 40), j] <- NA
  fac <- list(visual = c("x1", "x2", "x3"), textual = c("x4", "x5", "x6"),
              speed = c("x7", "x8", "x9"))
  d <- as.data.frame(suppressMessages(fa_confirmatory(hs, factors = fac, missing = "fiml",
                                                     aux = "ageyr")))
  fit <- suppressMessages(fa_confirmatory(hs, factors = fac, missing = "fiml",
                                         aux = "ageyr", output = "fit"))
  m <- sem_fit_measures(fit)
  expect_equal(d$estimate[d$term == "srmr"], m$value[m$term == "srmr"], tolerance = 1e-12)
  expect_gt(d$estimate[d$term == "srmr"],
            lavaan::fitMeasures(fit, "srmr_bentler_nomean")[[1L]] + 0.003)
})
