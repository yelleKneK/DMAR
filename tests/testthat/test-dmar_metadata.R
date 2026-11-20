# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

md_items <- function() c("t1_visual_perception", "t2_cubes", "t4_lozenges",
  "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
md_fields <- function() c("result", "coefficient", "estimator", "se_method", "test",
  "missing", "ci_method", "wald_df", "conf_level", "rmsea_conf_level",
  "tail_areas", "denominator",
  "interval_degenerate",
  "version", "version_reason", "correction", "N", "groups", "fit_rows",
  "auxiliary", "eligibility", "bootstrap", "simulation", "transformation",
  "methods",
  "rotation", "structure", "footer")

test_that("dmar_metadata() reads every family into one vocabulary", {
  d <- holzinger_swineford[md_items()]
  alpha <- dmar_metadata(reliability_alpha(d[1:3]))
  expect_identical(names(alpha), md_fields())
  expect_identical(alpha$result, "dmar_reliability")
  expect_identical(alpha$coefficient, "alpha")
  expect_identical(alpha$estimator, "analytic")
  expect_identical(alpha$ci_method, "bonett")
  expect_identical(alpha$missing, "listwise")
  expect_identical(alpha$eligibility, list(status = "eligible", code = NULL, note = NULL))
  expect_equal(alpha$N, nrow(d))
  expect_null(alpha$auxiliary)
  expect_null(alpha$bootstrap)
  expect_null(alpha$version)
  # A table without an RMSEA interval has no RMSEA level; one whose own
  # interval is the RMSEA's has its level there.
  expect_null(alpha$rmsea_conf_level)
  fm <- dmar_metadata(sem_fit_measures(chisq = 85.3, df = 24, N = 301,
                                       conf_level = 0.95))
  expect_identical(fm$rmsea_conf_level, 0.95)
  expect_identical(dmar_metadata(equivalence_rmsea(
    85.3, 24, 301, rmsea_bound = 0.08))$rmsea_conf_level, 0.90)
  # The auxiliary variables of a FIML fit are read from the family's record.
  d_aux <- d[1:3]
  d_aux$t5_general_information <- holzinger_swineford$t5_general_information
  d_aux[seq(1, nrow(d_aux), 7), 1] <- NA
  fiml <- dmar_metadata(reliability_alpha(d_aux, missing = "fiml",
                                          aux = "t5_general_information"))
  expect_identical(fiml$auxiliary, "t5_general_information")
  expect_identical(fiml$missing, "fiml")
  omega <- dmar_metadata(suppressWarnings(
    reliability_omega(d[1:3], ci_method = "percentile", B = 30, seed = 1)))
  expect_identical(omega$ci_method, "percentile")
  expect_identical(omega$denominator, "observed")
  expect_equal(omega$bootstrap$B, 30)
  expect_equal(omega$bootstrap$used + omega$bootstrap$failed, 30)
  expect_true(all(c("undetermined", "not_converged", "too_few_cases") %in% names(omega$bootstrap)))
  expect_error(dmar_metadata(data.frame(a = 1)), "dmar_tbl")
})

test_that("dmar_metadata() carries the version, correction, rows, and bootstrap counts of the SEM family", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[md_items()]
  model <- paste("visual  =~ t1_visual_perception + t2_cubes + t4_lozenges",
                 "textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
                 sep = "\n")
  fit <- lavaan::cfa(model, data = d, std.lv = TRUE, estimator = "MLR")
  f <- dmar_metadata(sem_fit_measures(fit))
  expect_identical(f$result, "dmar_sem_fit")
  expect_identical(f$version, "robust")
  expect_match(f$correction, "Brosseau-Liard")
  expect_identical(f$estimator, "ML")
  expect_match(f$test, "yuan.bentler")
  expect_equal(f$N, nrow(d))
  expect_equal(f$groups, 1)
  expect_identical(f$eligibility$status, "eligible")
  expect_true(is.character(f$footer) && length(f$footer) >= 2L)
  ic <- dmar_metadata(sem_information_criteria(fit))
  expect_identical(ic$result, "dmar_sem_information")
  expect_equal(ic$N, nrow(d))
  expect_identical(ic$estimator, "MLR")
  expect_equal(ic$groups, 1)
  expect_true(is.character(ic$footer))
  b <- dmar_metadata(sem_exact_fit_bootstrap(lavaan::cfa(model, data = d, std.lv = TRUE),
                                             B = 6, seed = 1))
  expect_equal(b$bootstrap$B, 6)
  expect_equal(b$bootstrap$used + b$bootstrap$failed, 6)
  expect_identical(b$transformation[[1L]]$method, "bollen_stine")
  expect_match(b$test, "likelihood ratio")
  e <- dmar_metadata(fa_exploratory(d, 2, estimator = "ULS", seed = 113))
  expect_identical(e$estimator, "ULS")
  expect_identical(e$se_method, "robust.sem.nt")
  expect_identical(e$rotation, "oblimin")
  expect_identical(e$structure, "correlation")
  expect_match(e$footer, "minimum residual")
  s <- dmar_metadata(factor_scores(fit))
  expect_identical(s$methods, c("regression", "Bartlett"))
  expect_null(names(s$methods))
  expect_identical(s$missing, "listwise")
  expect_equal(s$groups, 1)
  g <- generics::glance(fa_exploratory(d, 2, estimator = "ULS", seed = 113))
  expect_identical(g$structure, "correlation")
  set.seed(2)
  d2 <- d; d2[sample(nrow(d2), 60), "t2_cubes"] <- NA
  d2$t5_general_information <- holzinger_swineford$t5_general_information
  ts <- fa_confirmatory(d2, factors = list(visual = md_items()[1:3], textual = md_items()[4:6]),
              missing = "two_stage", aux = "t5_general_information")
  t <- dmar_metadata(ts)
  expect_identical(t$missing, "two_stage")
  expect_match(t$test, "residual based")
  expect_identical(t$auxiliary, "t5_general_information")
  expect_true(is.list(t$fit_rows))
  expect_equal(t$fit_rows$read, nrow(d2))
  expect_equal(t$fit_rows$nobs, nrow(d2))
  expect_identical(t$eligibility$status, "eligible")
  expect_true(is.character(t$footer) && length(t$footer) >= 1L)
})

test_that("a two-stage fit not computed again with every case is qualified, not eligible", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[md_items()]
  set.seed(2)
  d[sample(nrow(d), 60), "t2_cubes"] <- NA
  d$t5_general_information <- holzinger_swineford$t5_general_information
  local_mocked_bindings(.ts_all_rows = function(fit) {
    list(fit = fit, read = 301L, nobs = 290L, extra = 11L, corrected = FALSE,
         reestimated = FALSE)
  })
  expect_warning(
    ts <- fa_confirmatory(d, factors = list(visual = md_items()[1:3], textual = md_items()[4:6]),
                missing = "two_stage", aux = "t5_general_information"),
    class = "dmar_two_stage_n")
  m <- dmar_metadata(ts)
  expect_identical(m$eligibility$status, "qualified")
  expect_identical(m$eligibility$code, "two_stage_rows_unavailable")
  expect_match(m$eligibility$note, "could not compute the fit again with all 301")
  expect_match(m$eligibility$note, "^fa_confirmatory\\(\\)")
  expect_false(m$fit_rows$corrected)
})

test_that("a qualified result and a refusal are told apart by status and class", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  f <- rnorm(60)
  d <- as.data.frame(sapply(seq_len(6), function(j) .7 * f + rnorm(60, sd = sqrt(.51))))
  names(d) <- paste0("y", 1:6)
  d[cbind(5:60, rep(1:6, length.out = 56))] <- NA_real_
  expect_warning(x <- reliability_H(data = d, missing = "fiml", ci_method = "ml"),
                 class = "dmar_h1_unavailable")
  m <- dmar_metadata(x)
  expect_identical(m$eligibility$status, "qualified")
  expect_identical(m$eligibility$code, "h1_unavailable")
  expect_match(m$eligibility$note, "H0-based uncertainty")
  model <- paste0("f =~ ", paste0("y", 1:6, collapse = " + "))
  ref <- lavaan::cfa(model, data = d, std.lv = TRUE, missing = "ml")
  caught <- tryCatch(sem_fit_measures(ref), dmar_refusal = function(e) e)
  expect_s3_class(caught, c("dmar_h1_unavailable", "dmar_refusal", "error"))
  expect_identical(caught$code, "h1_unavailable")
  expect_match(caught$diagnostic, "H1|unrestricted|affine|collapse")
  expect_match(conditionMessage(caught), "Fit indices requiring the unrestricted H1")
  expect_error(sem_fit_measures(ref), class = "dmar_h1_unavailable")
  # fa_exploratory() keeps the constrained model when lavaan's H0 covariance matrix of
  # the estimates does not need the H1 information, and qualifies it; when
  # it does, the same class is a refusal.
  expect_warning(q <- fa_exploratory(d, 1, missing = "fiml"), class = "dmar_h1_unavailable")
  mq <- dmar_metadata(q)
  expect_identical(mq$eligibility$status, "qualified")
  expect_identical(mq$eligibility$code, "h1_unavailable")
  local_mocked_bindings(.fiml_h0_vcov_supported = function(fit) FALSE)
  caught_efa <- tryCatch(fa_exploratory(d, 1, missing = "fiml"), dmar_refusal = function(e) e)
  expect_s3_class(caught_efa, c("dmar_h1_unavailable", "dmar_refusal", "error"))
  expect_match(conditionMessage(caught_efa), "requires unavailable H1 information")
})

test_that("the two-stage all-row refusal is classed, with its row counts", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[md_items()]
  fit <- lavaan::cfa("f =~ t1_visual_perception + t2_cubes + t4_lozenges", data = d,
                     std.lv = TRUE)
  rows <- list(read = c(301L), nobs = c(290L), extra = 11L, corrected = FALSE,
               reestimated = FALSE, fit = fit)
  local_mocked_bindings(.ts_all_rows = function(fit) rows)
  caught <- tryCatch(sem_fit_measures(fit), dmar_refusal = function(e) e)
  expect_s3_class(caught, "dmar_two_stage_rows_unavailable")
  expect_identical(caught$code, "two_stage_rows_unavailable")
  expect_equal(caught$rows$extra, 11L)
  expect_match(conditionMessage(caught), "same all-row sample")
})

test_that(".dmar_refuse() builds the condition a caller catches", {
  e <- tryCatch(DMAR:::.dmar_refuse(c("a ", "b"), "example_code", extra = 1),
                error = function(e) e)
  expect_identical(class(e), c("dmar_example_code", "dmar_refusal", "error", "condition"))
  expect_identical(conditionMessage(e), "a b")
  expect_identical(e$code, "example_code")
  expect_equal(e$extra, 1)
  expect_null(conditionCall(e))
})

test_that("dmar_metadata() reads the tables whose method attribute is a label", {
  # randomization_test(), the meta-analysis functions, and ci_two_phase()
  # keep a character string in "method"; dmar_metadata() read it as the fit
  # families' list and stopped with "$ operator is invalid for atomic
  # vectors".
  r <- randomization_test(group_1 = c(12, 15, 11, 18, 14, 16, 13),
                          group_2 = c(9, 11, 10, 13, 8, 12, 10))
  m <- meta_es(c(0.21, 0.45, 0.10, 0.62, 0.33, 0.28),
               c(0.02, 0.05, 0.03, 0.08, 0.04, 0.025))
  set.seed(8)
  xx <- rnorm(40)
  yy <- 0.7 * xx + rnorm(40, sd = 0.7)
  tp <- ci_two_phase(yy, xx, rnorm(400))
  nn <- nnt_from_smd(smd = 0.5, n_1 = 50, n_2 = 50)
  for (res in list(r, m, tp, nn)) {
    md <- dmar_metadata(res)
    expect_type(md, "list")
    expect_identical(md$result, class(res)[1L])
    expect_null(md$estimator)
  }
  expect_identical(dmar_metadata(r)$conf_level, 0.95)
})

test_that("dmar_metadata() reports the replications mediate() dropped", {
  set.seed(113)
  n <- 40
  x <- c(rep(1, 8), rep(0, 32))
  w <- x
  w[c(1, 9)] <- 1 - w[c(1, 9)]
  m <- 0.5 * x + rnorm(n)
  y <- 0.4 * m + rnorm(n)
  res <- suppressWarnings(
    mediate(data.frame(x, m, y, w), "x", "m", "y", covariates = "w",
            B = 2000, seed = 113))
  boot <- dmar_metadata(res)$bootstrap
  expect_identical(boot$B, 2000L)
  expect_identical(boot$used, 1752L)
  expect_identical(boot$failed, 248L)
})
