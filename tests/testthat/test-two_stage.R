## Two-stage maximum likelihood through the CFA and SEM functions:
## fa_confirmatory() (and cfa_1(), cfa_2()), measurement_invariance(),
## sem_fit_measures(), sem_information_criteria(), and compare_sem() with
## missing = "two_stage" / "robust_two_stage" (lavaan's "two.stage" /
## "robust.two.stage").
##
## The oracle values were computed with lavaan 0.7-3 on fixture O1 (nine
## Holzinger and Swineford tests with 90 values of two deleted completely at
## random; see test-two_stage_internals.R); the residual statistic was
## checked against semTools 0.5-9 there.

.ts_o1 <- function() {
  hs <- holzinger_swineford[c(
    "t1_visual_perception", "t2_cubes", "t4_lozenges",
    "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning",
    "t10_addition", "t12_counting_groups_of_dots",
    "t13_straight_and_curved_capitals", "school", "age_years")]
  set.seed(1)
  hs$t7_sentence[sample(301, 90)] <- NA
  hs$t13_straight_and_curved_capitals[sample(301, 90)] <- NA
  for (v in c("t10_addition", "t12_counting_groups_of_dots",
              "t13_straight_and_curved_capitals")) hs[[v]] <- hs[[v]] / 10
  hs
}

.ts_o1_factors <- list(
  spatial = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
  verbal  = c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning"),
  speed   = c("t10_addition", "t12_counting_groups_of_dots",
              "t13_straight_and_curved_capitals"))

.ts_skip <- function() {
  skip_if_not_installed("lavaan")
  skip_if(utils::packageVersion("lavaan") < "0.7-3",
          "the two-stage oracles were computed with lavaan 0.7-3")
}

.row <- function(x, term, col = "estimate") x[[col]][x$term == term]

test_that("fa_confirmatory() reports the residual statistic and the TS V2 indices after two-stage estimation", {
  .ts_skip()
  x <- fa_confirmatory(.ts_o1(), .ts_o1_factors, missing = "two_stage", conf_level = .90)
  expect_equal(.row(x, "chi_square"), 65.2930632982, tolerance = 1e-6)
  expect_equal(.row(x, "df"), 24)
  expect_equal(.row(x, "p_chi_square"), 1.10240e-05, tolerance = 1e-4)
  expect_equal(.row(x, "rmsea"), 0.0869182753, tolerance = 1e-6)
  expect_equal(.row(x, "rmsea_ci_lower"), 0.0631545514, tolerance = 1e-6)
  expect_equal(.row(x, "rmsea_ci_upper"), 0.1114399248, tolerance = 1e-6)
  expect_equal(.row(x, "cfi"), 0.9374251252, tolerance = 1e-6)
  expect_equal(.row(x, "tli"), 0.9061376877, tolerance = 1e-6)
  # The SRMR moves with the optimizer's stopping point more than the
  # statistics do (5e-6 between this fit and lavaan's cfa() on the
  # unscaled tests).
  # The srmr row is the SRMR of the covariance residuals alone
  # (lavaan's srmr_bentler_nomean); through 1.0.x it was lavaan's srmr,
  # 0.0567490868 here, which averages in the mean residuals.
  expect_equal(.row(x, "srmr"), 0.0621652473, tolerance = 1e-5)
  expect_match(attr(x, "test"), "residual based statistic")
  ts <- attr(x, "two_stage")
  expect_identical(ts$lavaan_missing, "two.stage")
  expect_named(ts$footer, c("two_stage", "two_stage_test",
                             "two_stage_indices", "two_stage_ic"))
  out <- paste(utils::capture.output(print(x)), collapse = "\n")
  expect_match(out, "residual based statistic (Savalei & Bentler, 2009)",
               fixed = TRUE)
  expect_match(out, "Lai (2021)", fixed = TRUE)
  # Normal theory two-stage estimation gives the normal theory analog of
  # TS V2 (test-two_stage_internals.R builds both from the data), and
  # Zhang and Savalei define only the RMSEA and the CFI.
  expect_match(ts$footer[["two_stage_indices"]],
               paste0("^cfi, rmsea: the normal theory analog of the TS V2 ",
                      "indices of Zhang and Savalei \\(2023\\); tli, nnfi, ",
                      "and the RMSEA interval: lavaan's extension"))
  # The srmr row needs no line of its own.
  expect_false(any(grepl("srmr", attr(x, "fit_rows"), fixed = TRUE)))
  expect_false(grepl("srmr: lavaan's SRMR", out, fixed = TRUE))
  # The row set does not change.
  expect_identical(utils::tail(x$term, 15),
                   c("chi_square", "df", "p_chi_square", "cfi", "tli", "nnfi",
                     "rmsea", "rmsea_ci_lower", "rmsea_ci_upper",
                     "rmsea_ci_level", "srmr", "AIC", "BIC", "H0", "H1"))
  # lavaan's spelling is a synonym.
  y <- fa_confirmatory(.ts_o1(), .ts_o1_factors, missing = "two.stage", conf_level = .90)
  expect_equal(y$estimate, x$estimate)
  # Every other path names its statistic too.
  expect_identical(attr(fa_confirmatory(holzinger_swineford, .ts_o1_factors), "test"),
                   "standard")
  expect_identical(attr(fa_confirmatory(holzinger_swineford, .ts_o1_factors,
                              estimator = "MLR"), "test"), "scaled")
})

test_that("robust two-stage estimation reports the scaled statistic and its indices", {
  .ts_skip()
  x <- fa_confirmatory(.ts_o1(), .ts_o1_factors, missing = "robust_two_stage")
  expect_equal(.row(x, "chi_square"), 65.6657098130, tolerance = 1e-6)
  expect_equal(.row(x, "rmsea"), 0.0859051292, tolerance = 1e-6)
  expect_equal(.row(x, "cfi"), 0.9387864528, tolerance = 1e-6)
  expect_equal(.row(x, "tli"), 0.9081796791, tolerance = 1e-6)
  expect_match(attr(x, "test"), "scaled statistic")
  expect_match(attr(x, "two_stage")$footer[["two_stage_indices"]],
               "^cfi, rmsea: the TS V2 indices of Zhang and Savalei \\(2023\\)")
  # The estimates are the normal theory two-stage estimates; the standard
  # errors are lavaan's robust two-stage ones.
  n <- fa_confirmatory(.ts_o1(), .ts_o1_factors, missing = "two_stage")
  est <- !x$term %in% utils::tail(x$term, 15)
  expect_equal(x$estimate[est], n$estimate[est], tolerance = 1e-8)
  fit <- lavaan::cfa(paste(
    "spatial =~ t1_visual_perception + t2_cubes + t4_lozenges",
    "verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    "speed =~ t10_addition + t12_counting_groups_of_dots +",
    "t13_straight_and_curved_capitals", sep = "\n"), data = .ts_o1(),
    std.lv = TRUE, missing = "robust.two.stage")
  pe <- lavaan::parameterEstimates(fit)
  se_t7 <- pe$se[pe$op == "=~" & pe$rhs == "t7_sentence"]
  expect_equal(.row(x, "lambda_verbal_2", "se"), se_t7, tolerance = 1e-6)
})

test_that("fa_confirmatory() frees the intercepts under two-stage estimation (J = 3 identity with FIML)", {
  .ts_skip()
  # Through bare lavaan() a two-stage fit without intercept lines fixes
  # every intercept at 0, which on these uncentered tests gave omega
  # 0.9811142726 and a chi square of 1086 on 3 df. With three items the
  # single factor is saturated, so two-stage estimation gives the FIML
  # estimate and standard error.
  d <- holzinger_swineford[c("t1_visual_perception", "t2_cubes",
                             "t4_lozenges")]
  set.seed(1)
  d$t2_cubes[sample(301, 60)] <- NA
  ts <- cfa_1(d, missing = "two_stage")
  fiml <- cfa_1(d, missing = "fiml")
  expect_equal(.row(ts, "omega_f1"), 0.6646608062, tolerance = 1e-7)
  expect_equal(.row(ts, "omega_f1", "se"), 0.0390200280, tolerance = 1e-6)
  expect_equal(.row(ts, "omega_f1"), .row(fiml, "omega_f1"), tolerance = 1e-7)
  expect_equal(.row(ts, "omega_f1", "se"), .row(fiml, "omega_f1", "se"),
               tolerance = 1e-6)
  expect_gt(abs(.row(ts, "omega_f1") - 0.9811142726), 0.3)
  expect_true(attr(ts, "has_means"))
  expect_equal(.row(ts, "df"), 0)
  expect_equal(.row(ts, "chi_square"), 0, tolerance = 1e-8)
  # cfa_1() and cfa_2() pass missing through to fa_confirmatory().
  k <- fa_confirmatory(d, list(f1 = names(d)), missing = "two_stage")
  expect_equal(ts$estimate, k$estimate)
  hs <- .ts_o1()
  two <- cfa_2(hs, .ts_o1_factors$spatial, .ts_o1_factors$verbal,
               missing = "two_stage")
  k2 <- fa_confirmatory(hs, list(f1 = .ts_o1_factors$spatial,
                       f2 = .ts_o1_factors$verbal), missing = "two_stage")
  expect_equal(two$estimate, k2$estimate)
})

test_that("fa_confirmatory()'s two-stage AIC and BIC are Lai's criteria, as sem_information_criteria() reports them", {
  .ts_skip()
  hs <- .ts_o1()
  x <- fa_confirmatory(hs, .ts_o1_factors, missing = "two_stage")
  ic <- sem_information_criteria(fa_confirmatory(hs, .ts_o1_factors,
                                       missing = "two_stage",
                                       output = "fit"))
  expect_equal(.row(x, "AIC"), ic$value[ic$term == "AIC"])
  expect_equal(.row(x, "BIC"), ic$value[ic$term == "BIC"])
  expect_equal(.row(x, "H0"), ic$value[ic$term == "log_likelihood"])
  expect_match(attr(ic, "sem_information")$likelihood, "two-stage")
  # Lai (2021, Eq. 32): AIC_TS = n F_TS + 2q, where n F_TS is the naive
  # statistic -2(H0 - H1); the reported AIC differs from it by a constant
  # shared by every model of the same variables, -2 H1.
  q <- (.row(x, "AIC") + 2 * .row(x, "H0")) / 2
  expect_equal(.row(x, "AIC") + 2 * .row(x, "H1"),
               -2 * (.row(x, "H0") - .row(x, "H1")) + 2 * q)
})

test_that("fa_confirmatory() refuses two-stage requests it cannot honor", {
  .ts_skip()
  hs <- .ts_o1()
  f <- .ts_o1_factors
  expect_error(fa_confirmatory(hs, f, missing = "two_stage", estimator = "MLR"),
               "robust_two_stage")
  expect_error(fa_confirmatory(hs, f, missing = "two_stage", estimator = "WLS"),
               "is not available with it")
  expect_error(fa_confirmatory(hs, f, missing = "two_stage", se = "robust.sem"),
               "robust_two_stage")
  expect_error(fa_confirmatory(hs, f, missing = "two_stage", se = "bootstrap"),
               "must be \"standard\"")
  expect_error(fa_confirmatory(S = stats::cov(holzinger_swineford[unlist(f)]),
                     factors = f, N = 301, missing = "two_stage"),
               "requires raw 'data'")
  ord <- as.data.frame(lapply(holzinger_swineford[unlist(f)], function(v)
    as.integer(cut(v, stats::quantile(v, 0:4 / 4), include.lowest = TRUE))))
  expect_error(fa_confirmatory(ord, f, ordered = TRUE, missing = "two_stage"),
               "not available with ordered items")
  expect_error(fa_confirmatory(hs, f, missing = "two_stage", meanstructure = FALSE),
               "free item intercepts")
  # se = "none" is accepted and gives the two-stage standard errors.
  expect_s3_class(fa_confirmatory(hs, f, missing = "two_stage", se = "none"),
                  "dmar_fa_confirmatory")
})

test_that("fa_confirmatory() with an auxiliary variable in stage one", {
  .ts_skip()
  x <- fa_confirmatory(.ts_o1(), .ts_o1_factors, missing = "two_stage",
             aux = "age_years")
  # semTools 0.5-9: 64.6563565.
  expect_equal(.row(x, "chi_square"), 64.6563565, tolerance = 1e-6)
  m <- fa_confirmatory(.ts_o1(), .ts_o1_factors, missing = "two_stage",
             aux = "age_years", output = "measurement")
  expect_named(attr(m, "two_stage")$footer, "two_stage")
  # Robust two-stage estimation with an auxiliary reads lavaan's scaled
  # statistic and robust standard errors.
  fr <- lavaan::cfa(paste(
    "spatial =~ t1_visual_perception + t2_cubes + t4_lozenges",
    "verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    "speed =~ t10_addition + t12_counting_groups_of_dots +",
    "t13_straight_and_curved_capitals", sep = "\n"), data = .ts_o1(),
    std.lv = TRUE, missing = "robust.two.stage", aux = "age_years")
  r <- fa_confirmatory(.ts_o1(), .ts_o1_factors, missing = "robust_two_stage",
             aux = "age_years")
  expect_equal(.row(r, "chi_square"),
               lavaan::fitMeasures(fr, "chisq.scaled")[[1]], tolerance = 1e-6)
  pe <- lavaan::parameterEstimates(fr)
  expect_equal(.row(r, "lambda_verbal_2", "se"),
               pe$se[pe$op == "=~" & pe$rhs == "t7_sentence"],
               tolerance = 1e-6)
})

# Six tests with the whole block missing for 77 of 301 pupils, at random
# given general information, the auxiliary (26 at Pasteur, 51 at
# Grant-White).
.ts_block <- function() {
  hs <- holzinger_swineford
  set.seed(113)
  pr <- stats::plogis(-1.5 + 1.2 *
                        as.numeric(scale(hs$t5_general_information)))
  hs[sample(301, 77, prob = pr), .ts_block_items] <- NA
  hs
}
.ts_block_items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
                     "t6_paragraph_comprehension", "t7_sentence",
                     "t9_word_meaning")
.ts_block_factors <- list(spatial = .ts_block_items[1:3],
                          verbal = .ts_block_items[4:6])
.ts_block_model <- paste(
  "spatial =~ t1_visual_perception + t2_cubes + t4_lozenges",
  "verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
  sep = "\n")

test_that("a case observed only on an auxiliary variable is kept, and every value counts the cases the first stage read", {
  .ts_skip()
  # lavaan's two-stage estimator reads the 77 cases in its first stage and
  # leaves them out of its N (224), which makes every standard error
  # sqrt(301/224) times too large (test-two_stage_internals.R). fa_confirmatory()
  # keeps them and computes the fit again with N = 301. Through the first
  # versions of 1.1.0 it set them aside, which conditioned the estimates on
  # an observed item: the verbal loading was 2.796 and the intercept of
  # paragraph comprehension 8.754, against 2.840 and 9.210 with the cases
  # kept (2.972 and 9.183 with complete data).
  hs <- .ts_block()
  lav <- suppressWarnings(lavaan::cfa(
    .ts_block_model, data = hs, std.lv = TRUE, missing = "two.stage",
    aux = "t5_general_information", meanstructure = TRUE))
  expect_identical(lavaan::lavInspect(lav, "nobs"), 224L)
  # lavaan's warning that the 77 cases are empty and "will be ignored" is
  # wrong about them and is muffled; the message says what is done.
  m <- testthat::capture_messages(expect_no_warning(
    x <- fa_confirmatory(hs, .ts_block_factors, missing = "two_stage",
               aux = "t5_general_information")))
  expect_length(m, 1L)
  expect_match(m, paste0("fa_confirmatory(): 77 cases are observed on an auxiliary ",
                         "variable but on none of the model's variables. ",
                         "Two-stage estimation keeps them"), fixed = TRUE)
  expect_match(m, paste0("computed on all 301 cases it reads (lavaan's ",
                         "two-stage estimator counts only the 224"),
               fixed = TRUE)
  # The message ends its line, as base R's message() does, so the table
  # printed next starts on a line of its own (through the first versions
  # of 1.1.0 it ran into it).
  expect_true(endsWith(m, ".\n"))
  pe <- lavaan::parameterEstimates(lav)
  ld <- pe[pe$op == "=~", ]
  # The estimates are lavaan's two-stage estimates from every case.
  expect_equal(x$estimate[startsWith(x$term, "lambda_")], ld$est,
               tolerance = 1e-6)
  expect_equal(.row(x, "lambda_verbal_1"), 2.8396460577, tolerance = 1e-8)
  expect_equal(.row(x, "nu_verbal_1"), 9.2101744574, tolerance = 1e-8)
  # The standard errors are lavaan's times sqrt(224/301): the two-stage
  # sandwich over every case (test-two_stage_internals.R).
  expect_equal(x$se[startsWith(x$term, "lambda_")],
               ld$se * sqrt(224 / 301), tolerance = 1e-6)
  expect_equal(.row(x, "lambda_verbal_1", "se"), 0.1901966739,
               tolerance = 1e-7)
  # The fit rows, each from lavaan's naive statistic and scaling constants
  # with N = 301 and the statistics multiplied by 301/224, by the formulas
  # of Zhang and Savalei (2023) and lavaan's extension of them, written
  # out here.
  fm <- lavaan::fitMeasures(lav)
  r <- 301 / 224
  T <- fm[["chisq"]] * r
  c_m <- fm[["chisq.scaling.factor"]]
  df <- fm[["df"]]
  Tb <- fm[["baseline.chisq"]] * r
  c_b <- fm[["baseline.chisq.scaling.factor"]]
  df_b <- fm[["baseline.df"]]
  # The naive statistic over the 301 cases is the one semTools 0.5-9's
  # twostage(), whose first stage counts every case, prints (36.160).
  expect_equal(T, 36.160, tolerance = 1e-4)
  ncp <- function(p) stats::uniroot(function(l) {
    stats::pchisq(T / c_m, df, ncp = l) - p
  }, c(0, 4 * T), tol = 1e-12)$root
  expect_equal(.row(x, "chi_square"), .ts_residual_test(lav)$statistic)
  expect_equal(.row(x, "chi_square"), 23.7685737907, tolerance = 1e-8)
  expect_equal(.row(x, "rmsea"),
               sqrt(max(T - c_m * df, 0) / (301 * df)), tolerance = 1e-8)
  # The RMSEA interval is at rmsea_conf_level, 90 percent by default.
  expect_equal(.row(x, "rmsea_ci_lower"),
               sqrt(c_m * ncp(0.95) / (301 * df)), tolerance = 1e-6)
  expect_equal(.row(x, "rmsea_ci_upper"),
               sqrt(c_m * ncp(0.05) / (301 * df)), tolerance = 1e-6)
  expect_equal(.row(x, "cfi"), 1 - max(T - c_m * df, 0) /
                 max(T - c_m * df, Tb - c_b * df_b, 0), tolerance = 1e-8)
  expect_equal(.row(x, "tli"), 1 - ((T - c_m * df) * df_b) /
                 ((Tb - c_b * df_b) * df), tolerance = 1e-8)
  # The log likelihoods are 301 times the complete data log likelihood per
  # case at the EM estimates, and BIC takes log 301.
  ss <- lavaan::lavInspect(lav, "sampstat")
  imp <- lavaan::lavInspect(lav, "implied")
  e <- ss$mean - imp$mean
  per_case <- -0.5 * (length(e) * log(2 * pi) +
                        as.numeric(determinant(imp$cov)$modulus) +
                        sum(diag(solve(imp$cov, ss$cov))) +
                        sum(e * solve(imp$cov, e)))
  expect_equal(.row(x, "H0"), 301 * per_case, tolerance = 1e-10)
  expect_equal(.row(x, "H0"), -5429.08810301, tolerance = 1e-10)
  expect_equal(.row(x, "H1"), fm[["unrestricted.logl"]] * r,
               tolerance = 1e-10)
  expect_equal(.row(x, "AIC"), -2 * 301 * per_case + 2 * fm[["npar"]],
               tolerance = 1e-10)
  expect_equal(.row(x, "BIC"),
               -2 * 301 * per_case + fm[["npar"]] * log(301),
               tolerance = 1e-10)
  # Robust two-stage estimation: lavaan's scaled statistic times 301/224.
  fr <- suppressWarnings(lavaan::cfa(
    .ts_block_model, data = hs, std.lv = TRUE,
    missing = "robust.two.stage", aux = "t5_general_information",
    meanstructure = TRUE))
  xr <- suppressMessages(fa_confirmatory(hs, .ts_block_factors,
                               missing = "robust_two_stage",
                               aux = "t5_general_information"))
  expect_equal(.row(xr, "chi_square"),
               lavaan::fitMeasures(fr, "chisq.scaled")[[1]] * r,
               tolerance = 1e-8)
  per <- lavaan::parameterEstimates(fr)
  expect_equal(xr$se[startsWith(xr$term, "lambda_")],
               per$se[per$op == "=~"] * sqrt(224 / 301), tolerance = 1e-6)
  # A line under the table names the N, and the attribute keeps the counts.
  ts <- attr(x, "two_stage")
  expect_identical(ts$rows$read, 301L)
  expect_identical(ts$rows$nobs, 224L)
  expect_identical(ts$rows$extra, 77L)
  expect_true(ts$rows$corrected)
  expect_match(ts$footer[["two_stage_n"]],
               "^N: the 301 cases the first stage read, 77 of them")
  mm <- suppressMessages(fa_confirmatory(hs, .ts_block_factors, missing = "two_stage",
                               aux = "t5_general_information",
                               output = "measurement"))
  expect_named(attr(mm, "two_stage")$footer, c("two_stage", "two_stage_n"))
  expect_equal(.row(mm, "omega_verbal", "se"),
               .row(x, "omega_verbal", "se"))
  # The lavaan fit itself counts the 301 cases, with lavaan's estimates and
  # its record of their estimation.
  xf <- suppressMessages(fa_confirmatory(hs, .ts_block_factors, missing = "two_stage",
                               aux = "t5_general_information",
                               output = "fit"))
  expect_identical(lavaan::lavInspect(xf, "nobs"), 301L)
  expect_identical(nrow(lavaan::lavInspect(xf, "data")), 301L)
  expect_equal(unname(lavaan::coef(xf)[c(paste0("lambda_spatial_", 1:3),
                                         paste0("lambda_verbal_", 1:3))]),
               ld$est, tolerance = 1e-6)
  expect_gt(lavaan::lavInspect(xf, "iterations"), 0L)
  expect_identical(lavaan::lavInspect(xf, "options")$optim.method, "nlminb")
  sx <- suppressMessages(fa_confirmatory(hs, .ts_block_factors, missing = "two_stage",
                               aux = "t5_general_information",
                               output = "standardized"))
  sl <- lavaan::standardizedSolution(lav)
  expect_equal(sx$se[sx$op == "=~"], sl$se[sl$op == "=~"] * sqrt(224 / 301),
               tolerance = 1e-6)
  # cfa_1() and cfa_2() pass through fa_confirmatory(), and measurement_invariance()
  # gives the message once for its ladder.
  expect_message(cfa_1(hs, items = .ts_block_items[1:3],
                       missing = "two_stage",
                       aux = "t5_general_information", output = "fit"),
                 "fa_confirmatory\\(\\): 77 cases")
  mi_msg <- testthat::capture_messages(
    measurement_invariance(hs, .ts_block_factors, group = "school",
                           missing = "two_stage",
                           aux = "t5_general_information",
                           levels = c("configural", "metric")))
  expect_identical(sum(grepl("measurement_invariance(): 77 cases", mi_msg,
                             fixed = TRUE)), 1L)
  # A case with nothing observed, not even the auxiliary, is dropped by
  # lavaan from both stages, with its warning, and is not counted.
  he <- hs
  he[1, c(.ts_block_items, "t5_general_information")] <- NA
  expect_warning(
    xe <- suppressMessages(fa_confirmatory(he, .ts_block_factors,
                                 missing = "two_stage",
                                 aux = "t5_general_information")),
    "some cases are empty")
  expect_identical(attr(xe, "two_stage")$rows$read, 300L)
  # Without an auxiliary there is nothing to correct, and no line.
  x0 <- suppressWarnings(fa_confirmatory(hs, .ts_block_factors, missing = "two_stage"))
  expect_true(is.na(attr(x0, "two_stage")$rows$corrected))
  expect_false("two_stage_n" %in% names(attr(x0, "two_stage")$footer))
})

test_that("when the fit cannot be computed again, fa_confirmatory() says which values rest on lavaan's N", {
  .ts_skip()
  # The path for a lavaan that does not keep its sample statistics where
  # .ts_refit_n() sets them, simulated.
  hs <- .ts_block()
  y <- suppressMessages(fa_confirmatory(hs, .ts_block_factors, missing = "two_stage",
                              aux = "t5_general_information"))
  local_mocked_bindings(.ts_refit_n = function(fit, n, ...) stop("not here"))
  w <- testthat::capture_warnings(m <- testthat::capture_messages(
    x <- fa_confirmatory(hs, .ts_block_factors, missing = "two_stage",
               aux = "t5_general_information")))
  expect_length(w, 1L)
  # The message that every value counts the 301 cases is not given, since
  # none does; the warning alone says what rests on lavaan's N (through the
  # first versions of 1.1.0 both were given).
  expect_length(m, 0L)
  expect_match(w, paste0("fa_confirmatory(): the first stage of the two-stage fit ",
                         "read 301 cases, 77 of them observed on an ",
                         "auxiliary variable alone, which lavaan's N of 224 ",
                         "leaves out, and DMAR could not compute the fit ",
                         "again with all 301"), fixed = TRUE)
  expect_match(w, "about 1.159 times too large", fixed = TRUE)
  expect_false(attr(x, "two_stage")$rows$corrected)
  expect_match(attr(x, "two_stage")$footer[["two_stage_n"]],
               "^N: lavaan's 224, which leaves out 77 cases")
  expect_warning(
    xm <- suppressMessages(fa_confirmatory(hs, .ts_block_factors,
                                 missing = "two_stage",
                                 aux = "t5_general_information",
                                 output = "measurement")),
    class = "dmar_two_stage_n")
  # Every standard error, of a parameter or a defined quantity, is then
  # lavaan's, sqrt(301/224) times the one computed from every case.
  has <- is.finite(x$se) & x$se > 0
  expect_equal(x$se[has] / y$se[has], rep(sqrt(301 / 224), sum(has)),
               tolerance = 1e-6)
  # The chi_square row of a normal theory two-stage fit is the residual
  # based statistic, which sums the information over every case read, so
  # it is the same whether or not the fit was computed again, and neither
  # the line nor the warning says it rests on lavaan's N (they spoke of a
  # scaled statistic, which the table does not report). The indices and
  # criteria do rest on it.
  val <- function(t, term) t$estimate[t$term == term]
  expect_equal(val(x, "chi_square"), val(y, "chi_square"), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(val(x, "AIC"), val(y, "AIC"))))
  n_line <- attr(x, "two_stage")$footer[["two_stage_n"]]
  expect_false(grepl("scaled statistic", n_line, fixed = TRUE))
  expect_false(grepl("scaled statistic", w, fixed = TRUE))
  expect_match(n_line, paste0("the indices and criteria rest on 224 cases; ",
                              "the residual based statistic in chi_square ",
                              "uses all 301"), fixed = TRUE)
  expect_match(w, paste0("the fit indices and the information criteria rest ",
                         "on 224 cases; the residual based statistic in ",
                         "chi_square uses every case read"), fixed = TRUE)
  # After robust two-stage estimation the row is the scaled statistic,
  # which does rest on lavaan's N.
  wr <- testthat::capture_warnings(suppressMessages(
    xr <- fa_confirmatory(hs, .ts_block_factors, missing = "robust_two_stage",
                aux = "t5_general_information")))
  expect_match(attr(xr, "two_stage")$footer[["two_stage_n"]],
               "the scaled statistic, indices, and criteria rest on 224",
               fixed = TRUE)
  expect_match(wr, paste0("the scaled statistic, the fit indices, and the ",
                          "information criteria rest on 224 cases"),
               fixed = TRUE, all = FALSE)
})

test_that("measurement_invariance() under two-stage estimation", {
  .ts_skip()
  hs <- .ts_o1()
  mi <- measurement_invariance(hs, .ts_o1_factors, group = "school",
                               missing = "two_stage",
                               levels = c("configural", "metric"))
  expect_equal(mi$chi_square, c(81.4134204798, 89.8059318463),
               tolerance = 1e-6)
  expect_equal(mi$df, c(48, 54))
  # The rung test is lavaan's scaled difference (Satorra & Bentler, 2001).
  expect_equal(mi$delta_chi_square[2], 7.9180536574, tolerance = 1e-6)
  expect_equal(mi$delta_df[2], 6)
  expect_equal(mi$p_value[2], 0.2441686084, tolerance = 1e-6)
  expect_equal(mi$cfi, c(0.9379526571, 0.9353246187), tolerance = 1e-6)
  expect_match(attr(mi, "test"), "satorra.bentler.2001")
  expect_match(attr(mi, "model_test"), "residual based statistic")
  expect_identical(attr(mi, "fit_indices"), "robust")
  out <- paste(utils::capture.output(print(mi)), collapse = "\n")
  expect_match(out, "scaled difference test (Satorra & Bentler, 2001)",
               fixed = TRUE)
  # lavaan's spelling, and the robust variant, which reports the scaled
  # statistic of each rung.
  mi2 <- measurement_invariance(hs, .ts_o1_factors, group = "school",
                                missing = "two.stage", levels = "configural")
  expect_equal(mi2$chi_square, mi$chi_square[1])
  mr <- measurement_invariance(hs, .ts_o1_factors, group = "school",
                               missing = "robust_two_stage",
                               levels = c("configural", "metric"))
  fm <- lavaan::fitMeasures(attr(mr, "fits")$configural)
  expect_equal(mr$chi_square[1], fm[["chisq.scaled"]])
  expect_equal(mr$cfi[1], fm[["cfi.robust"]])
  expect_match(attr(mr, "model_test"), "scaled statistic")
  # The existing paths name their statistic as well.
  expect_identical(attr(measurement_invariance(
    holzinger_swineford, .ts_o1_factors[1:2], group = "school",
    levels = "configural"), "model_test"), "standard")
  expect_named(attr(mi, "two_stage")$footer,
               c("two_stage", "two_stage_test", "two_stage_indices",
                 "two_stage_diff"))
  expect_match(attr(mi, "two_stage")$footer[["two_stage_indices"]],
               "^cfi, rmsea: the normal theory analog of the TS V2")
})

test_that("measurement_invariance() with an auxiliary variable in each group's first stage", {
  .ts_skip()
  # The residual statistic of each group aligns the group's auxiliary rows
  # with its item rows, including rows with a missing auxiliary value and
  # rows whose group is missing, which lavaan drops. semTools 0.5-9's
  # twostage(), anova(), on the same data printed 80.275 and 78.355 on 48
  # df (unrounded 80.2750295703 and 78.3554788717; its stage one is a FIML
  # fit of the saturated model, which stops at a looser tolerance than
  # lavaan's EM). semTools cannot check the metric rung, whose constraints
  # it projects on the wrong parameters (see test-two_stage_internals.R).
  hs <- .ts_o1()
  hs$age_years[c(10, 11, 160)] <- NA
  mi <- measurement_invariance(hs, .ts_o1_factors, group = "school",
                               missing = "two_stage", aux = "age_years",
                               levels = c("configural", "metric"))
  expect_equal(mi$chi_square, c(80.2752102357, 88.2769888968),
               tolerance = 1e-6)
  expect_equal(mi$chi_square[1], 80.2750295703, tolerance = 1e-5)
  expect_equal(mi$df, c(48, 54))
  hn <- hs
  hn$school[c(3, 50, 200, 290)] <- NA
  mn <- suppressWarnings(measurement_invariance(
    hn, .ts_o1_factors, group = "school", missing = "two_stage",
    aux = "age_years", levels = "configural"))
  expect_equal(mn$chi_square, 78.3556875782, tolerance = 1e-6)
  expect_equal(mn$chi_square, 78.3554788717, tolerance = 1e-5)
  # Robust two-stage estimation with an auxiliary reads lavaan's scaled
  # statistic and robust indices.
  mr <- measurement_invariance(hs, .ts_o1_factors, group = "school",
                               missing = "robust_two_stage",
                               aux = "age_years", levels = "configural")
  fm <- lavaan::fitMeasures(attr(mr, "fits")$configural)
  expect_equal(mr$chi_square, fm[["chisq.scaled"]])
  expect_equal(mr$cfi, fm[["cfi.robust"]])
})

test_that("measurement_invariance() computes every rung again with the cases its first stage read", {
  .ts_skip()
  # The six tests by school, the block missing for 77 pupils at random
  # given general information: 26 such cases at Pasteur and 51 at
  # Grant-White, so the two groups' counts grow in different ratios (130 to
  # 156 and 94 to 145). Each rung's fit is computed again with the counts
  # of the cases read. The second stage weights each group's discrepancy
  # by its share of N, so a rung that constrains across groups is
  # estimated again: its estimates are the two-stage estimates with every
  # case read, those of an ordinary maximum likelihood fit of the stage one
  # moments with those counts (test-two_stage_internals.R checks the
  # standard errors and the scaling constants against the data). Through
  # the first versions of 1.1.0 the cases were set aside, and then
  # lavaan's estimates were kept with the new counts, which gave a scalar
  # difference test of 15.21 and a strict one of 12.66 (p = .049).
  hs <- .ts_block()
  d <- hs[c(.ts_block_items, "t5_general_information", "school")]
  run <- function(missing, levels = c("configural", "metric", "scalar",
                                      "strict")) {
    suppressMessages(measurement_invariance(
      d, .ts_block_factors, group = "school", missing = missing,
      aux = "t5_general_information", levels = levels))
  }
  m <- testthat::capture_messages(mi <- measurement_invariance(
    d, .ts_block_factors, group = "school", missing = "two_stage",
    aux = "t5_general_information"))
  expect_length(m, 1L)
  expect_match(m, paste0("and the estimates of the rungs that constrain ",
                         "across groups, the test statistics, fit indices, ",
                         "and difference tests of every rung are computed ",
                         "on all 301 cases it reads"), fixed = TRUE)
  expect_true(endsWith(m, ".\n"))
  fits <- attr(mi, "fits")
  ge <- list(character(0), "loadings", c("loadings", "intercepts"),
             c("loadings", "intercepts", "residuals"))
  lav <- lapply(ge, function(g) {
    args <- list(model = attr(mi, "model"), data = d, group = "school",
                 missing = "two.stage", aux = "t5_general_information",
                 meanstructure = TRUE)
    if (length(g)) args$group.equal <- g
    suppressWarnings(do.call(lavaan::cfa, args))
  })
  n_g <- lavaan::lavInspect(lav[[1L]], "nobs")
  read_g <- n_g + c(26L, 51L)
  ss <- lavaan::lavInspect(lav[[1L]], "sampstat")
  for (i in seq_along(fits)) {
    expect_identical(lavaan::lavInspect(fits[[i]], "nobs"), read_g)
    ml_args <- list(model = attr(mi, "model"),
                    sample.cov = lapply(ss, function(x) unclass(x$cov)),
                    sample.mean = lapply(ss, function(x) unclass(x$mean)),
                    sample.nobs = read_g, sample.cov.rescale = FALSE,
                    likelihood = "normal", meanstructure = TRUE)
    if (length(ge[[i]])) ml_args$group.equal <- ge[[i]]
    ml <- do.call(lavaan::cfa, ml_args)
    expect_equal(lavaan::coef(fits[[i]]), lavaan::coef(ml), tolerance = 1e-5)
    expect_equal(lavaan::fitMeasures(fits[[i]], "chisq")[[1]],
                 lavaan::fitMeasures(ml, "chisq")[[1]], tolerance = 1e-9)
    if (i == 1L) {
      # The configural rung constrains nothing across groups: lavaan's
      # estimates, and each group's naive statistic times its ratio.
      expect_identical(lavaan::coef(fits[[i]]), lavaan::coef(lav[[i]]))
      expect_equal(
        lavaan::fitMeasures(fits[[i]], "chisq")[[1]],
        sum(lavaan::lavInspect(lav[[i]], "test")[[1L]]$stat.group *
              read_g / n_g), tolerance = 1e-10)
    } else {
      expect_gt(max(abs(lavaan::coef(fits[[i]]) - lavaan::coef(lav[[i]]))),
                0.05)
    }
    # The residual based statistic reads every case, at the rung's own
    # estimates.
    expect_equal(mi$chi_square[i], .ts_residual_test(fits[[i]])$statistic)
  }
  # semTools 0.5-9, twostage(model, data, group = "school", aux =
  # "t5_general_information", fun = "cfa", group.equal = c("loadings",
  # "intercepts"), meanstructure = TRUE), 2026-11-02: naive statistic
  # 64.13918 and log likelihood -5393.338 on the scalar rung.
  expect_equal(lavaan::fitMeasures(fits$scalar, "chisq")[[1]], 64.13918,
               tolerance = 1e-6)
  expect_equal(lavaan::fitMeasures(fits$scalar, "logl")[[1]], -5393.338,
               tolerance = 1e-6)
  # The common intercept of lozenges at the scalar rung is 18.00 (SE .697),
  # where lavaan's is 18.21 (SE .782).
  pe <- lavaan::parameterEstimates(fits$scalar)
  lz <- pe[pe$lhs == "t4_lozenges" & pe$op == "~1", ]
  expect_equal(lz$est, rep(18.00168, 2), tolerance = 1e-6)
  expect_equal(lz$se, rep(0.6970, 2), tolerance = 1e-3)
  # The difference tests are the scaled differences of the fits returned,
  # from naive statistics and scaling constants that count every case.
  sb <- function(a, b) {
    fa <- lavaan::fitMeasures(a, c("chisq", "df", "chisq.scaling.factor"))
    fb <- lavaan::fitMeasures(b, c("chisq", "df", "chisq.scaling.factor"))
    cd <- (fb[[2]] * fb[[3]] - fa[[2]] * fa[[3]]) / (fb[[2]] - fa[[2]])
    (fb[[1]] - fa[[1]]) / cd
  }
  expect_equal(mi$delta_chi_square[2:4],
               c(sb(fits[[1]], fits[[2]]), sb(fits[[2]], fits[[3]]),
                 sb(fits[[3]], fits[[4]])),
               tolerance = 1e-8)
  expect_equal(mi$delta_chi_square[3:4], c(14.968078, 12.503334),
               tolerance = 1e-5)
  expect_equal(mi$p_value[4], 0.0516372, tolerance = 1e-4)
  # lavaan's own ladder on its N gives another scalar test.
  expect_equal(lavaan::lavTestLRT(lav[[2]], lav[[3]])[2L, "Chisq diff"],
               10.8638, tolerance = 1e-4)
  rows <- attr(mi, "two_stage")$rows
  expect_identical(rows$read, unname(read_g))
  expect_true(rows$corrected)
  expect_true(rows$reestimated)
  expect_match(attr(mi, "two_stage")$footer[["two_stage_n"]],
               paste0("the estimates of the rungs that constrain across ",
                      "groups and the statistics, indices, and difference ",
                      "tests of every rung use all 301"), fixed = TRUE)
  # The configural rung alone has nothing to estimate again.
  mc <- run("two_stage", "configural")
  expect_false(attr(mc, "two_stage")$rows$reestimated)
  expect_match(attr(mc, "two_stage")$footer[["two_stage_n"]],
               "; the statistics, indices, and difference tests of every rung")
  # When one rung cannot be computed again (simulated), none is, so the
  # rungs compared share one N, and the warning and the line say so.
  orig <- .ts_refit_n
  local_mocked_bindings(.ts_refit_n = function(fit, n, ...) {
    if ("intercepts" %in% lavaan::lavInspect(fit, "options")$group.equal) {
      stop("not here")
    }
    orig(fit, n, ...)
  })
  w <- testthat::capture_warnings(m <- testthat::capture_messages(
    mf <- measurement_invariance(d, .ts_block_factors, group = "school",
                                 missing = "two_stage",
                                 aux = "t5_general_information",
                                 levels = c("configural", "metric",
                                            "scalar"))))
  expect_length(m, 0L)
  expect_length(w, 1L)
  expect_match(w, paste0("The estimates of a model that constrains ",
                         "parameters across groups then weight the groups ",
                         "by lavaan's counts"), fixed = TRUE)
  expect_identical(lavaan::lavInspect(attr(mf, "fits")$metric, "nobs"), n_g)
  expect_equal(mf$delta_chi_square[3], 10.8638, tolerance = 1e-4)
  expect_false(attr(mf, "two_stage")$rows$corrected)
  expect_match(attr(mf, "two_stage")$footer[["two_stage_n"]],
               paste0("^N: lavaan's 224.*the estimates of the rungs that ",
                      "constrain across groups weight the groups by ",
                      "lavaan's counts, which moves the residual based ",
                      "statistics of those rungs as well$"))
  # Under normal theory the chi_square column holds the residual based
  # statistics, which use every case read: the configural rung's is the
  # one of the ladder computed again, and only the indices and difference
  # tests rest on lavaan's N (the line and the warning spoke of scaled
  # statistics, and the warning of information criteria the table does
  # not report).
  expect_equal(mf$chi_square[1], mi$chi_square[1], tolerance = 1e-10)
  expect_false(grepl("scaled statistic",
                     attr(mf, "two_stage")$footer[["two_stage_n"]],
                     fixed = TRUE))
  expect_match(attr(mf, "two_stage")$footer[["two_stage_n"]],
               paste0("the indices and difference tests rest on 224 cases ",
                      "(the residual based statistics in chi_square use ",
                      "every case read)"), fixed = TRUE)
  expect_false(grepl("scaled statistic|information criteria", w))
  expect_match(w, paste0("the fit indices, and the difference tests rest on ",
                         "224 cases; the residual based statistics in ",
                         "chi_square use every case read"), fixed = TRUE)
})

test_that("two-stage estimation reads lavaan's ceq.simple = TRUE", {
  .ts_skip()
  # lavaan's ceq.simple = TRUE keeps simple equality constraints without
  # == rows; the fit is the same model, and the residual statistic is the
  # one without it (it was refused for a rank it did not have).
  hs <- .ts_o1()
  lv <- c("configural", "metric")
  mi <- measurement_invariance(hs, .ts_o1_factors, group = "school",
                               missing = "two_stage", levels = lv)
  mc <- measurement_invariance(hs, .ts_o1_factors, group = "school",
                               missing = "two_stage", levels = lv,
                               ceq.simple = TRUE)
  expect_equal(mc$chi_square, mi$chi_square, tolerance = 1e-6)
  expect_equal(mc$df, mi$df)
  expect_equal(mc$delta_chi_square, mi$delta_chi_square, tolerance = 1e-4)
  f <- .ts_o1_factors["spatial"]
  k <- fa_confirmatory(hs, f, equal_loading = TRUE, missing = "two_stage")
  kc <- fa_confirmatory(hs, f, equal_loading = TRUE, missing = "two_stage",
              ceq.simple = TRUE)
  expect_equal(.row(kc, "chi_square"), .row(k, "chi_square"),
               tolerance = 1e-6)
  expect_equal(.row(kc, "df"), .row(k, "df"))
})

test_that("measurement_invariance() refuses two-stage requests it cannot honor", {
  .ts_skip()
  hs <- .ts_o1()
  f <- .ts_o1_factors
  expect_error(measurement_invariance(hs, f, group = "school",
                                      missing = "two_stage",
                                      estimator = "MLR"),
               "robust_two_stage")
  expect_error(measurement_invariance(hs, f, group = "school",
                                      missing = "two_stage",
                                      se = "robust.sem"),
               "robust_two_stage")
  expect_error(measurement_invariance(hs, f, group = "school",
                                      missing = "two_stage", ordered = TRUE),
               "not available with ordered items")
})

test_that("sem_fit_measures() reports a two-stage fit as fa_confirmatory() does", {
  .ts_skip()
  hs <- .ts_o1()
  model <- paste(
    "spatial =~ t1_visual_perception + t2_cubes + t4_lozenges",
    "verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    "speed =~ t10_addition + t12_counting_groups_of_dots +",
    "t13_straight_and_curved_capitals", sep = "\n")
  fit <- lavaan::cfa(model, data = hs, missing = "two.stage")
  x <- sem_fit_measures(fit, conf_level = 0.90)
  k <- fa_confirmatory(hs, .ts_o1_factors, missing = "two_stage", conf_level = .90)
  val <- function(t) x$value[x$term == t]
  # The two fits differ only in identification (a marker indicator here,
  # unit factor variances in fa_confirmatory()), so they agree to the optimizers'
  # tolerance.
  expect_equal(val("chi_square"), .row(k, "chi_square"), tolerance = 1e-6)
  expect_equal(val("chi_square"), 65.2930632982, tolerance = 1e-6)
  expect_equal(val("rmsea"), .row(k, "rmsea"), tolerance = 1e-6)
  expect_equal(val("rmsea_ci_lower"), .row(k, "rmsea_ci_lower"),
               tolerance = 1e-6)
  expect_equal(val("cfi"), .row(k, "cfi"), tolerance = 1e-6)
  expect_equal(val("tli"), .row(k, "tli"), tolerance = 1e-6)
  det <- attr(x, "sem_fit")
  expect_identical(det$version, "robust")
  expect_match(det$version_reason, "two-stage")
  expect_match(det$test, "residual based statistic")
  expect_match(det$correction, "TS V2")
  # After normal theory two-stage estimation the indices are the normal
  # theory analog of TS V2 (test-two_stage_internals.R builds both from
  # the data), and the footer says so, as fa_confirmatory()'s does.
  expect_identical(det$correction,
                   paste0("two-stage corrected indices of Zhang and Savalei ",
                          "(2023), the normal theory analog of TS V2"))
  expect_match(det$footer[2], "the normal theory analog of TS V2$")
  expect_error(sem_fit_measures(fit, version = "standard"),
               "not chi square distributed")
  expect_error(sem_fit_measures(fit, version = "standard"),
               "the corrected indices of Zhang and Savalei (2023).",
               fixed = TRUE)
  expect_error(sem_fit_measures(fit, version = "scaled"),
               "population values")
  fr <- lavaan::cfa(model, data = hs, missing = "robust.two.stage")
  xr <- sem_fit_measures(fr)
  expect_equal(xr$value[xr$term == "chi_square"],
               lavaan::fitMeasures(fr, "chisq.scaled")[[1]])
  expect_match(attr(xr, "sem_fit")$test, "scaled statistic")
  # After robust two-stage estimation they are TS V2 itself.
  expect_identical(attr(xr, "sem_fit")$correction,
                   "two-stage corrected indices of Zhang and Savalei (2023), TS V2")
  # Cases that observe nothing are dropped by lavaan with no claim of
  # listwise deletion.
  he <- hs
  he[1:3, unlist(.ts_o1_factors)] <- NA
  fe <- suppressWarnings(lavaan::cfa(model, data = he, missing = "two.stage"))
  expect_no_warning(sem_fit_measures(fe))
})

test_that("sem_fit_measures() corrects every endpoint to the cases its first stage read", {
  .ts_skip()
  # The data of the test of fa_confirmatory() above: six tests with the whole block
  # missing for 77 of 301 pupils, at random given general information, the
  # auxiliary. A lavaan fit of the user's own keeps those cases in its
  # first stage and out of its N (lavaan warns that they are empty).
  hs <- holzinger_swineford
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  set.seed(113)
  pr <- stats::plogis(-1.5 + 1.2 *
                        as.numeric(scale(hs$t5_general_information)))
  hs[sample(301, 77, prob = pr), items] <- NA
  model <- paste(
    "spatial =~ t1_visual_perception + t2_cubes + t4_lozenges",
    "verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    sep = "\n")
  two_stage_fit <- function(data, missing = "two.stage", ...) {
    suppressWarnings(lavaan::cfa(model, data = data, std.lv = TRUE,
                                 missing = missing, ...))
  }
  fit <- two_stage_fit(hs, aux = "t5_general_information")
  expect_identical(lavaan::lavInspect(fit, "nobs"), 224L)
  expect_identical(nrow(lavaan::lavInspect(fit, "data")), 301L)
  expect_no_warning(x <- sem_fit_measures(fit))
  val <- function(x, t) x$value[x$term == t]
  expect_equal(val(x, "N"), 301)
  rows <- attr(x, "sem_fit")$two_stage_rows
  expect_true(rows$corrected)
  expect_equal(sum(rows$read), 301)
  expect_equal(rows$extra, 77)
  expect_output(print(x), "auxiliary")
  corrected <- .ts_all_rows(fit)$fit
  expect_equal(val(x, "chi_square"), .ts_residual_test(corrected)$statistic)
  expect_equal(val(x, "rmsea"), lavaan::fitMeasures(corrected, "rmsea.robust")[[1]])
  robust <- two_stage_fit(hs, "robust.two.stage", aux = "t5_general_information")
  expect_no_warning(xr <- sem_fit_measures(robust))
  expect_equal(val(xr, "N"), 301)
  expect_equal(val(xr, "chi_square"),
    lavaan::fitMeasures(.ts_all_rows(robust)$fit, "chisq.scaled")[[1]])
  fg <- two_stage_fit(hs, aux = "t5_general_information", group = "school")
  expect_identical(sum(lavaan::lavInspect(fg, "nobs")), 224L)
  expect_no_warning(xg <- sem_fit_measures(fg))
  expect_equal(val(xg, "N"), 301)
  expect_true(attr(xg, "sem_fit")$two_stage_rows$corrected)
  h1 <- holzinger_swineford
  h1[5, items] <- NA
  expect_no_warning(single <- sem_fit_measures(two_stage_fit(h1, aux = "t5_general_information")))
  expect_equal(val(single, "N"), 301)
  expect_equal(attr(single, "sem_fit")$two_stage_rows$extra, 1)
  # No warning for the fit fa_confirmatory() returns, which keeps the cases and
  # counts them (its N row is 301, and its rows are those of lavaan with
  # that N), when the auxiliary is missing on them too (they are empty, and
  # lavaan drops them from both stages), without an auxiliary, or under
  # FIML.
  xf <- suppressMessages(fa_confirmatory(hs, list(spatial = items[1:3],
                                        verbal = items[4:6]),
                               missing = "two_stage",
                               aux = "t5_general_information",
                               output = "fit"))
  expect_no_warning(xs <- sem_fit_measures(xf))
  expect_equal(val(xs, "N"), 301)
  expect_equal(val(xs, "chi_square"), val(x, "chi_square"))
  expect_equal(val(xs, "rmsea"),
               lavaan::fitMeasures(xf, "rmsea.robust")[[1]])
  h2 <- hs
  h2$t5_general_information[rowSums(!is.na(hs[items])) == 0] <- NA
  expect_no_warning(
    sem_fit_measures(two_stage_fit(h2, aux = "t5_general_information")))
  expect_no_warning(sem_fit_measures(two_stage_fit(hs)))
  expect_no_warning(sem_fit_measures(two_stage_fit(hs, "ml")))
})

test_that("sem_information_criteria() and compare_sem() count every case the first stage read", {
  .ts_skip()
  # The data of the tests above: the six tests missing as a block for 77
  # of 301 pupils, at random given general information, the auxiliary.
  # lavaan forms the two-stage log likelihood as its N, 224, times the
  # complete data log likelihood per case at EM estimates that read all
  # 301 rows, and the BIC, SABIC, and HBIC penalties with the same N. The
  # criteria are computed instead with every case: the log likelihood is
  # 301 times the one per case, and N is 301. Through the first versions
  # of 1.1.0 they were lavaan's (BIC 8183.32), with a warning.
  hs <- .ts_block()
  items <- .ts_block_items
  one <- paste("g =~", paste(items, collapse = " + "))
  two_stage_fit <- function(model, missing = "two.stage", data = hs, ...) {
    suppressWarnings(lavaan::cfa(model, data = data, std.lv = TRUE,
                                 missing = missing, ...))
  }
  fit2 <- two_stage_fit(.ts_block_model, aux = "t5_general_information")
  fit1 <- two_stage_fit(one, aux = "t5_general_information")
  val <- function(x, t) x$value[x$term == t]
  expect_no_warning(ic <- sem_information_criteria(fit2))
  ss <- lavaan::lavInspect(fit2, "sampstat")
  imp <- lavaan::lavInspect(fit2, "implied")
  e <- ss$mean - imp$mean
  per_case <- -0.5 * (length(e) * log(2 * pi) +
                        as.numeric(determinant(imp$cov)$modulus) +
                        sum(diag(solve(imp$cov, ss$cov))) +
                        sum(e * solve(imp$cov, e)))
  # lavaan's own log likelihood is 224 times the one per case.
  expect_equal(as.numeric(lavaan::logLik(fit2)), 224 * per_case)
  expect_equal(val(ic, "log_likelihood"), 301 * per_case, tolerance = 1e-10)
  expect_equal(val(ic, "log_likelihood"), -5429.08810301, tolerance = 1e-10)
  expect_equal(val(ic, "N"), 301)
  d <- val(ic, "npar")
  expect_equal(val(ic, "AIC"), -2 * 301 * per_case + 2 * d)
  expect_equal(val(ic, "BIC"), -2 * 301 * per_case + d * log(301))
  expect_equal(val(ic, "BIC"), 10966.611301, tolerance = 1e-10)
  expect_equal(val(ic, "SABIC"), -2 * 301 * per_case + d * log(303 / 24))
  expect_equal(val(ic, "HBIC"), -2 * 301 * per_case + d * log(301 / (2 * pi)))
  rows <- attr(ic, "sem_information")$two_stage_rows
  expect_identical(c(rows$read, rows$nobs, rows$extra), c(301L, 224L, 77L))
  expect_match(attr(ic, "sem_information")$footer,
               "^Two-stage fit: the log likelihood and N count the 301 cases",
               all = FALSE)
  # In several groups each group's log likelihood is multiplied by the
  # ratio of its own counts (26 such cases at Pasteur, 51 at Grant-White).
  fg <- two_stage_fit(.ts_block_model, aux = "t5_general_information",
                      group = "school")
  icg <- sem_information_criteria(fg)
  n_g <- lavaan::lavInspect(fg, "nobs")
  read_g <- vapply(split(hs, hs$school)[lavaan::lavInspect(fg, "group.label")],
                   function(dg) sum(rowSums(!is.na(dg[c(
                     items, "t5_general_information")])) > 0), integer(1))
  expect_identical(unname(read_g - n_g), c(26L, 51L))
  expect_equal(val(icg, "log_likelihood"),
               sum(fg@loglik$loglik.group * read_g / n_g), tolerance = 1e-10)
  expect_equal(val(icg, "N"), 301)
  expect_false(attr(icg, "sem_information")$two_stage_rows$reestimated)
  # With loadings and intercepts equal across the groups the estimates
  # depend on the groups' counts, so the log likelihood is the one at the
  # two-stage estimates with every case read: an ordinary maximum
  # likelihood fit of the stage one moments with those counts, and
  # semTools 0.5-9's twostage() (-5393.338). At lavaan's estimates it would
  # be -5393.519, and lavaan's own is -4020.785.
  fs <- two_stage_fit(.ts_block_model, aux = "t5_general_information",
                      group = "school", meanstructure = TRUE,
                      group.equal = c("loadings", "intercepts"))
  ics <- sem_information_criteria(fs)
  ss <- lavaan::lavInspect(fs, "sampstat")
  ml <- lavaan::cfa(.ts_block_model, std.lv = TRUE, meanstructure = TRUE,
                    sample.cov = lapply(ss, function(x) unclass(x$cov)),
                    sample.mean = lapply(ss, function(x) unclass(x$mean)),
                    sample.nobs = unname(read_g), sample.cov.rescale = FALSE,
                    likelihood = "normal",
                    group.equal = c("loadings", "intercepts"))
  expect_equal(val(ics, "log_likelihood"),
               lavaan::fitMeasures(ml, "logl")[[1]], tolerance = 1e-11)
  expect_equal(val(ics, "log_likelihood"), -5393.338, tolerance = 1e-6)
  expect_equal(val(ics, "BIC"), 10957.8898558, tolerance = 1e-9)
  expect_equal(as.numeric(lavaan::logLik(fs)), -4020.78541387,
               tolerance = 1e-10)
  expect_true(attr(ics, "sem_information")$two_stage_rows$reestimated)
  expect_match(attr(ics, "sem_information")$footer,
               "so the estimates were computed again with every case read",
               fixed = TRUE, all = FALSE)
  # The scalar rung measurement_invariance() returns is the same fit.
  mi <- suppressMessages(measurement_invariance(
    hs[c(items, "t5_general_information", "school")], .ts_block_factors,
    group = "school", missing = "two_stage", std.lv = TRUE,
    aux = "t5_general_information",
    levels = c("configural", "metric", "scalar")))
  expect_equal(val(sem_information_criteria(attr(mi, "fits")$scalar), "BIC"),
               val(ics, "BIC"), tolerance = 1e-8)
  # The fit fa_confirmatory() returns already counts the 301 cases, so the criteria
  # are the same.
  f2 <- suppressMessages(fa_confirmatory(hs, .ts_block_factors, missing = "two_stage",
                               aux = "t5_general_information",
                               output = "fit"))
  f1 <- suppressMessages(fa_confirmatory(hs, list(g = items), missing = "two_stage",
                               aux = "t5_general_information",
                               output = "fit"))
  expect_equal(val(sem_information_criteria(f2), "BIC"), val(ic, "BIC"))
  # compare_sem(): no warning, the criteria of sem_information_criteria(),
  # N = 301, and a line that says why; a fit of the user's and the fit of
  # fa_confirmatory() of the same data compare.
  expect_no_warning(x <- compare_sem(list(two = fit2, one = fit1)))
  expect_equal(x$value[match(c("two", "one"), x$term)],
               c(val(ic, "BIC"), val(sem_information_criteria(fit1), "BIC")))
  expect_identical(attr(x, "comparison")$N, 301)
  expect_match(attr(x, "comparison")$footer,
               "count the 77 cases the first stage read", all = FALSE)
  expect_no_warning(xm <- compare_sem(list(two = f2, one = fit1)))
  expect_equal(xm$value, x$value)
  expect_no_warning(compare_sem(list(two = f2, one = f1), criterion = "AIC"))
  # When the fit cannot be computed again (simulated), the criteria are
  # lavaan's and the warning names what rests on its N; compare_sem()
  # gives it once for each fit, naming it.
  local_mocked_bindings(.ts_refit_n = function(fit, n, ...) stop("not here"))
  w <- testthat::capture_warnings(icl <- sem_information_criteria(fit2))
  expect_length(w, 1L)
  expect_match(w, paste0("The first stage of this two-stage fit read 77 ",
                         "cases observed on an auxiliary variable but on ",
                         "none of the model's variables, which lavaan ",
                         "leaves out of its N of 224."), fixed = TRUE)
  expect_match(w, "could not compute the fit again", fixed = TRUE)
  expect_equal(val(icl, "BIC"),
               lavaan::fitMeasures(fit2, "bic")[[1]])
  expect_equal(val(icl, "BIC"), 8183.32449806, tolerance = 1e-8)
  expect_equal(val(icl, "N"), 224)
  wc <- testthat::capture_warnings(compare_sem(list(two = fit2, one = fit1)))
  expect_length(wc, 2L)
  expect_match(wc[1], paste0("The first stage of the two-stage fit fits$two ",
                             "read 77 cases"), fixed = TRUE)
  expect_match(wc[2], "fits$one", fixed = TRUE)
  # No warning without an auxiliary, or under FIML.
  expect_no_warning(sem_information_criteria(two_stage_fit(.ts_block_model)))
  expect_no_warning(sem_information_criteria(
    two_stage_fit(.ts_block_model, "ml")))
})

test_that("compare_sem() compares two-stage fits by Lai's criteria and reports no evidence", {
  .ts_skip()
  hs <- .ts_o1()
  three <- fa_confirmatory(hs, .ts_o1_factors, missing = "two_stage", output = "fit")
  one <- fa_confirmatory(hs, list(g = unlist(.ts_o1_factors, use.names = FALSE)),
               missing = "two_stage", output = "fit")
  x <- compare_sem(list(three = three, one = one))
  bic <- c(three = sem_information_criteria(three)$value[2],
           one = sem_information_criteria(one)$value[2])
  expect_equal(x$value[match(names(bic), x$term)], unname(bic))
  expect_true(all(is.na(x$log_evidence_best_over_model)))
  expect_match(attr(x, "comparison")$footer[1], "Lai (2021)", fixed = TRUE)
  expect_match(attr(x, "comparison")$footer[3], "no evidence")
  # A FIML fit and a two-stage fit have different likelihoods.
  fiml <- fa_confirmatory(hs, .ts_o1_factors, missing = "fiml", output = "fit")
  expect_error(compare_sem(list(three = three, fiml = fiml)),
               "handle missing data differently")
  # Different auxiliary variables move the EM estimates.
  aux <- fa_confirmatory(hs, list(g = unlist(.ts_o1_factors, use.names = FALSE)),
               missing = "two_stage", aux = "age_years", output = "fit")
  expect_error(compare_sem(list(three = three, aux = aux)),
               "different auxiliary variables")
})
