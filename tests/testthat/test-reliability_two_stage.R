# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

## Two-stage maximum likelihood in the reliability family:
## reliability_omega(), reliability_alpha(), and reliability() with
## missing = "two_stage" (and aux).
##
## The oracles are lavaan 0.7-3 two-stage fits specified from scratch here,
## and the FIML results the two-stage ones must equal where theory says so
## (three items, where the model is saturated; analytic alpha, a function
## of the first stage alone).

# Five verbal tests of the Holzinger and Swineford battery with values of
# y3 and y5 deleted at random given y1, which is always observed (MAR):
# 67 and 74 values missing, 165 complete cases.
.ts_verbal <- function() {
  d <- holzinger_swineford[c("t5_general_information",
                             "t6_paragraph_comprehension", "t7_sentence",
                             "t8_word_classification", "t9_word_meaning")]
  names(d) <- paste0("y", 1:5)
  set.seed(20261101)
  hi <- d$y1 > stats::median(d$y1)
  d$y3[stats::runif(301) < ifelse(hi, 0.45, 0.10)] <- NA
  d$y5[stats::runif(301) < ifelse(hi, 0.10, 0.45)] <- NA
  d
}

# The single-factor model of the five items with labeled loadings and
# error variances, free intercepts, and omega as a defined parameter.
.ts_verbal_syntax <- paste(
  "f =~ NA*y1 + l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5",
  "f ~~ 1*f",
  paste0("y", 1:5, " ~~ p", 1:5, "*y", 1:5, collapse = "\n"),
  paste0("y", 1:5, " ~ 1", collapse = "\n"),
  "omega := (l1+l2+l3+l4+l5)^2/((l1+l2+l3+l4+l5)^2 + p1+p2+p3+p4+p5)",
  sep = "\n")

.ts_skip <- function() {
  skip_if_not_installed("lavaan")
  skip_if(utils::packageVersion("lavaan") < "0.7-3",
          "the two-stage oracles were computed with lavaan 0.7-3")
}

.val <- function(x, term) x$value[x$term == term]

test_that("model implied omega under two-stage estimation is lavaan's, with both standard errors", {
  .ts_skip()
  d <- .ts_verbal()
  r <- reliability_omega(d, missing = "two_stage",
                         denominator = "model_implied", ci_method = "ml")
  expect_equal(.val(r, "estimate"), 0.902130, tolerance = 1e-6)
  expect_equal(.val(r, "se"), 0.009928, tolerance = 1e-4)
  expect_equal(.val(r, "N"), 301)
  expect_equal(.val(r, "N_complete"), 165)
  expect_identical(attr(r, "missing"), "two_stage")
  fit <- lavaan::lavaan(.ts_verbal_syntax, data = d, missing = "two.stage")
  pe <- lavaan::parameterEstimates(fit)
  expect_equal(.val(r, "estimate"), pe$est[pe$lhs == "omega"],
               tolerance = 1e-7)
  expect_equal(.val(r, "se"), pe$se[pe$lhs == "omega"], tolerance = 1e-6)
  rr <- reliability_omega(d, missing = "two_stage",
                          denominator = "model_implied", ci_method = "mlr")
  fr <- lavaan::lavaan(.ts_verbal_syntax, data = d,
                       missing = "robust.two.stage")
  per <- lavaan::parameterEstimates(fr)
  expect_equal(.val(rr, "estimate"), .val(r, "estimate"), tolerance = 1e-8)
  expect_equal(.val(rr, "se"), 0.009611, tolerance = 1e-4)
  expect_equal(.val(rr, "se"), per$se[per$lhs == "omega"], tolerance = 1e-6)
  # The default interval of the model implied coefficient is the robust one.
  rd <- reliability_omega(d, missing = "two_stage",
                          denominator = "model_implied")
  expect_equal(.val(rd, "se"), .val(rr, "se"))
})

test_that("robust omega under two-stage estimation divides by the EM composite variance", {
  .ts_skip()
  d <- .ts_verbal()
  expect_message(reliability_omega(d, missing = "two_stage"),
                 "bootstrap based")
  r <- suppressMessages(reliability_omega(d, missing = "two_stage"))
  expect_equal(.val(r, "estimate"), 0.899805, tolerance = 1e-6)
  fit <- lavaan::lavaan(.ts_verbal_syntax, data = d, missing = "two.stage")
  pe <- lavaan::parameterEstimates(fit)
  em <- unclass(lavaan::lavInspect(fit, "sampstat")$cov)
  expect_equal(.val(r, "estimate"),
               sum(pe$est[pe$op == "=~"])^2 / sum(em), tolerance = 1e-7)
  expect_identical(attr(r, "denominator"), "observed")
  # Its interval is a bootstrap that refits both stages; a delta method
  # interval is not offered.
  expect_error(reliability_omega(d, missing = "two_stage", ci_method = "ml"),
               "model_implied")
})

test_that("the two-stage bootstrap refits both stages and is reproducible", {
  .ts_skip()
  skip_on_cran()
  d <- .ts_verbal()
  # B = 20 with an 80 percent interval keeps the percentile limits off the
  # extreme order statistics while the test stays quick.
  a <- reliability_omega(d, missing = "two_stage", ci_method = "percentile",
                         B = 20, seed = 113, conf_level = 0.80)
  b <- reliability_omega(d, missing = "two_stage", ci_method = "percentile",
                         B = 20, seed = 113, conf_level = 0.80)
  expect_equal(a$value, b$value)
  bs <- attr(a, "bootstrap")
  expect_equal(bs$attempted, 20)
  expect_equal(bs$returned, 20)
  # Every replication refit both stages and returned a proper estimate,
  # spread about the point estimate.
  expect_true(all(bs$status == "proper"))
  expect_lt(abs(stats::median(bs$t) - .val(a, "estimate")), 0.02)
  m <- reliability_omega(d, missing = "two_stage",
                         denominator = "model_implied",
                         ci_method = "percentile", B = 20, seed = 113,
                         conf_level = 0.80)
  expect_equal(attr(m, "bootstrap")$returned, 20)
  expect_lt(.val(m, "ci_lower"), .val(m, "estimate"))
  expect_gt(.val(m, "ci_upper"), .val(m, "estimate"))
})

test_that("a two-stage bootstrap replication fits the second stage alone and gives the two-stage coefficients", {
  .ts_skip()
  # lavaan's two-stage fit computes the two-stage standard errors and the
  # scaled statistic whatever is asked, and a replication reads neither,
  # so a replication fits the complete data model to the EM moments
  # (.omega_fit_cfa_ts_estimates()). Its coefficients are the two-stage
  # fit's to rounding error where the likelihood has a well defined
  # minimum, as on these resamples of 301 cases, with an auxiliary
  # variable, for the congeneric and the tau-equivalent model and for
  # robust omega (on a small, heavily incomplete resample the two
  # optimizers can stop at different points; see the comment on the
  # function); its fit carries no standard error.
  d <- .ts_verbal()
  d$z <- holzinger_swineford$t1_visual_perception
  d$z[c(4, 40, 140)] <- NA
  items <- paste0("y", 1:5)
  set.seed(113)
  for (r in 1:2) {
    b <- d[sample(nrow(d), replace = TRUE), ]
    for (aux in list(NULL, "z")) {
      full <- .omega_fit_cfa(data = b[, c(items, aux)], se = "none",
                             missing = "two.stage", aux = aux)
      fast <- .omega_fit_cfa_ts_estimates(b, items = items, aux = aux)
      expect_equal(fast$omega, full$omega, tolerance = 1e-10)
      expect_equal(fast$loadings, full$loadings, tolerance = 1e-8)
      expect_equal(fast$errors, full$errors, tolerance = 1e-8)
      tau_full <- .omega_fit_cfa(data = b[, c(items, aux)], se = "none",
                                 equal_loadings = TRUE,
                                 missing = "two.stage", aux = aux)
      tau_fast <- .omega_fit_cfa_ts_estimates(b, items = items, aux = aux,
                                              equal_loadings = TRUE)
      expect_equal(tau_fast$omega, tau_full$omega, tolerance = 1e-10)
      ro_full <- .omega_observed_ts_parts(b, items = items, aux = aux)
      ro_fast <- .omega_observed_ts_parts(b, items = items, aux = aux,
                                          estimates_only = TRUE)
      expect_equal(ro_fast$estimate, ro_full$estimate, tolerance = 1e-10)
    }
  }
  expect_identical(lavaan::lavInspect(fast$fit, "options")$se, "none")
  expect_identical(lavaan::lavInspect(full$fit, "options")$se, "two.stage")
  # Three items with no single-factor solution are refused before a fit,
  # and a replication of one case, in which every pair is observed
  # together in one row only, is refused before stage one.
  # y3 is built to correlate negatively with y1 and positively with y2,
  # which correlate positively with each other: the weight on y1 lies
  # between r12 and 1/r12.
  bad <- d[items[1:3]]
  r12 <- stats::cor(bad$y1, bad$y2)
  bad$y3 <- as.numeric(scale(bad$y2)) -
    (r12 + 1 / r12) / 2 * as.numeric(scale(bad$y1)) +
    stats::rnorm(301, 0, 0.3)
  expect_lt(prod(stats::cov(bad)[c(2, 3, 6)]), 0)
  none <- .omega_fit_cfa_ts_estimates(bad, items = items[1:3])
  expect_true(none$no_solution)
  # Before the refusal of a pair observed together in one row, lavaan
  # printed its table of the variables and refused the case, and the
  # replication returned no fit.
  expect_error(.omega_fit_cfa_ts_estimates(d[1, ], items = items),
               "y1 with y2 (one row observes both)", fixed = TRUE)
})

test_that("a stage one that cannot be estimated leaves the two-stage replication with no fit", {
  .ts_skip()
  # The test of this branch read a one-row replication, which the refusal
  # of a pair observed together in one row now stops before stage one, so
  # no test reached it. The refusals leave no data on which the EM fails,
  # so the failure is simulated.
  d <- .ts_verbal()
  items <- paste0("y", 1:5)
  d[c(3, 30, 90), "y2"] <- NA
  test_frame <- environment()
  ok <- .omega_fit_cfa_ts_estimates(d, items = items)
  expect_true(ok$converged)
  local_mocked_bindings(.ts_em_moments = function(...) stop("simulated"),
                        .env = test_frame)
  failed <- .omega_fit_cfa_ts_estimates(d, items = items)
  expect_false(failed$converged)
  expect_false(failed$no_solution)
  expect_identical(failed, .omega_fit_none(5L, items))
  local_mocked_bindings(.ts_em_moments = function(x, ...) {
    p <- ncol(x)
    list(cov = matrix(NA_real_, p, p), mean = rep(NA_real_, p))
  }, .env = test_frame)
  expect_identical(.omega_fit_cfa_ts_estimates(d, items = items),
                   .omega_fit_none(5L, items))
})

test_that("with three items two-stage omega is the FIML omega", {
  .ts_skip()
  d <- holzinger_swineford[c("t1_visual_perception", "t2_cubes",
                             "t4_lozenges")]
  set.seed(1)
  d$t2_cubes[sample(301, 60)] <- NA
  ts <- reliability_omega(d, missing = "two_stage",
                          denominator = "model_implied", ci_method = "ml")
  fi <- reliability_omega(d, missing = "fiml",
                          denominator = "model_implied", ci_method = "ml")
  expect_equal(.val(ts, "estimate"), 0.6646608062, tolerance = 1e-7)
  expect_equal(.val(ts, "estimate"), .val(fi, "estimate"), tolerance = 1e-7)
  expect_equal(.val(ts, "se"), .val(fi, "se"), tolerance = 1e-6)
  # The saturated model reproduces every residual correlation.
  expect_true(attr(ts, "residual_cor")$saturated)
})

test_that("auxiliary variables enter the first stage only", {
  .ts_skip()
  d <- .ts_verbal()
  d$z <- holzinger_swineford$t1_visual_perception
  r <- reliability_omega(d, missing = "two_stage", aux = "z",
                         denominator = "model_implied", ci_method = "ml")
  fit <- lavaan::lavaan(.ts_verbal_syntax, data = d, missing = "two.stage",
                        aux = "z")
  pe <- lavaan::parameterEstimates(fit)
  expect_equal(.val(r, "estimate"), pe$est[pe$lhs == "omega"],
               tolerance = 1e-7)
  expect_equal(.val(r, "se"), pe$se[pe$lhs == "omega"], tolerance = 1e-6)
  expect_identical(attr(r, "aux"), "z")
  expect_identical(attr(r, "missing"), "two_stage")
  ro <- suppressMessages(reliability_omega(d, missing = "two_stage",
                                           aux = "z"))
  em <- unclass(lavaan::lavInspect(fit, "sampstat")$cov)
  expect_equal(.val(ro, "estimate"),
               sum(pe$est[pe$op == "=~"])^2 / sum(em), tolerance = 1e-7)
  # aux without missing still means FIML.
  expect_identical(attr(suppressMessages(reliability_omega(d, aux = "z")),
                        "missing"), "fiml")
})

test_that("analytic alpha under two-stage estimation is FIML's analytic alpha", {
  .ts_skip()
  d <- .ts_verbal()
  ts <- reliability_alpha(d, missing = "two_stage")
  fi <- reliability_alpha(d, missing = "fiml")
  expect_identical(ts$term, fi$term)
  expect_equal(ts$value, fi$value)
  expect_identical(attr(ts, "ci_method"), "ml")
  expect_identical(attr(ts, "missing"), "two_stage")
  tl <- reliability_alpha(d, missing = "two_stage", ci_method = "ml_logistic")
  fl <- reliability_alpha(d, missing = "fiml", ci_method = "ml_logistic")
  expect_equal(tl$value, fl$value)
})

test_that("model implied alpha under two-stage estimation is the tau-equivalent two-stage fit", {
  .ts_skip()
  d <- .ts_verbal()
  r <- reliability_alpha(d, missing = "two_stage",
                         estimator = "model_implied", ci_method = "ml")
  k <- cfa_1(d, equal_loading = TRUE, missing = "two_stage")
  expect_equal(.val(r, "estimate"),
               k$estimate[k$term == "omega_f1"], tolerance = 1e-8)
  expect_equal(.val(r, "se"), k$se[k$term == "omega_f1"], tolerance = 1e-6)
  rr <- reliability_alpha(d, missing = "two_stage",
                          estimator = "model_implied")
  kr <- cfa_1(d, equal_loading = TRUE, missing = "robust_two_stage")
  expect_equal(.val(rr, "se"), kr$se[kr$term == "omega_f1"], tolerance = 1e-6)
})

test_that("reliability() forwards two-stage estimation to alpha and omega only", {
  .ts_skip()
  d <- .ts_verbal()
  w <- reliability(d, type = "omega", missing = "two_stage",
                   denominator = "model_implied", ci_method = "ml")
  o <- reliability_omega(d, missing = "two_stage",
                         denominator = "model_implied", ci_method = "ml")
  expect_equal(w$value, o$value)
  wa <- reliability(d, type = "alpha", missing = "two_stage")
  expect_equal(wa$value, reliability_alpha(d, missing = "two_stage")$value)
  expect_error(reliability(d, type = "kr20", missing = "two_stage"),
               "KR-20 analyzes the complete rows only")
  # Categorical omega checks the treatment against its own choices.
  expect_error(reliability(d, type = "omega_categorical",
                           missing = "two_stage"),
               "should be one of")
})

test_that("the two-stage refusals say what to use instead", {
  .ts_skip()
  d <- .ts_verbal()
  S <- stats::cov(stats::na.omit(d))
  expect_error(reliability_omega(S = S, N = 165, missing = "two_stage"),
               "requires raw 'data'")
  expect_error(reliability_alpha(S = S, N = 165, missing = "two_stage"),
               "two-stage maximum likelihood to use")
  for (m in c("feldt", "bonett")) {
    expect_error(reliability_omega(d, missing = "two_stage",
                                   denominator = "model_implied",
                                   ci_method = m),
                 "two-stage sandwich")
  }
  expect_error(reliability_alpha(d, missing = "two_stage",
                                 ci_method = "bonett"),
               "complete-data closed form")
  expect_error(reliability_omega(d, missing = "two_stage",
                                 denominator = "model_implied",
                                 ci_method = "adf"),
               "two-stage maximum likelihood")
  expect_error(reliability_omega(d, missing = "two_stage",
                                 denominator = "model_implied",
                                 ci_method = "likelihood"),
               "maximizes no likelihood")
  expect_error(reliability_alpha(d, missing = "two_stage",
                                 estimator = "model_implied",
                                 ci_method = "likelihood"),
               "maximizes no likelihood")
  d$z <- holzinger_swineford$t1_visual_perception
  expect_error(reliability_omega(d, missing = "listwise", aux = "z"),
               "incomplete data method")
  # The FIML refusal keeps its wording.
  expect_error(reliability_omega(d[1:5], missing = "fiml",
                                 denominator = "model_implied",
                                 ci_method = "likelihood"),
               "has not been implemented for missing = \"fiml\"")
})

test_that("a case observed only on an auxiliary variable is kept and counted, with a message once per call", {
  .ts_skip()
  # The first stage reads such a case, and the family's fits go through
  # fa_confirmatory(), which computes them again with every case read; through the
  # first versions of 1.1.0 the family set the case aside (silently at
  # first, then with a message), which conditioned the estimates on an
  # observed item.
  d <- .ts_verbal()
  d$z <- holzinger_swineford$t1_visual_perception
  d[c(1, 4, 9), paste0("y", 1:5)] <- NA
  msg <- paste0("3 cases are observed on an auxiliary variable but on none ",
                "of the model's variables. Two-stage estimation keeps them")
  m <- testthat::capture_messages(
    r <- reliability_omega(d, missing = "two_stage", aux = "z",
                           denominator = "model_implied", ci_method = "ml"))
  expect_identical(sum(grepl(msg, m, fixed = TRUE)), 1L)
  expect_match(m, "reliability_omega(): 3 cases", fixed = TRUE, all = FALSE)
  expect_match(m, "the coefficient is estimated from all 301 cases it reads",
               fixed = TRUE, all = FALSE)
  expect_true(all(endsWith(m, "\n")))
  expect_identical(.val(r, "N"), 301)
  # The estimate is lavaan's two-stage estimate from every case, and the
  # standard error lavaan's times sqrt(298/301), the sandwich over every
  # case.
  fit <- suppressWarnings(lavaan::lavaan(.ts_verbal_syntax, data = d,
                                         missing = "two.stage", aux = "z"))
  expect_identical(lavaan::lavInspect(fit, "nobs"), 298L)
  pe <- lavaan::parameterEstimates(fit)
  expect_equal(.val(r, "estimate"), pe$est[pe$lhs == "omega"],
               tolerance = 1e-7)
  expect_equal(.val(r, "se"), pe$se[pe$lhs == "omega"] * sqrt(298 / 301),
               tolerance = 1e-6)
  # Analytic alpha is then FIML's, value and interval, as its help says.
  m <- testthat::capture_messages(
    a <- reliability_alpha(d, missing = "two_stage", aux = "z"))
  expect_identical(sum(grepl(msg, m, fixed = TRUE)), 1L)
  expect_match(m, "reliability_alpha(): 3 cases", fixed = TRUE, all = FALSE)
  af <- suppressMessages(reliability_alpha(d, missing = "fiml", aux = "z"))
  expect_identical(a$term, af$term)
  expect_equal(a$value, af$value)
  # Model implied alpha: the tau-equivalent two-stage fit, from every case.
  ai <- suppressMessages(reliability_alpha(d, missing = "two_stage",
                                           aux = "z",
                                           estimator = "model_implied",
                                           ci_method = "ml"))
  kz <- suppressMessages(cfa_1(d, items = paste0("y", 1:5),
                               equal_loading = TRUE, missing = "two_stage",
                               aux = "z"))
  expect_equal(.val(ai, "se"), kz$se[kz$term == "omega_f1"], tolerance = 1e-8)
  expect_equal(.val(ai, "N"), 301)
  # A bootstrap gives the message once, not once per replication, and its
  # replications resample all 301 cases.
  skip_on_cran()
  m <- testthat::capture_messages(
    rb <- reliability_omega(d, missing = "two_stage", aux = "z",
                            denominator = "model_implied",
                            ci_method = "bootstrap_se", B = 5, seed = 113))
  expect_identical(sum(grepl(msg, m, fixed = TRUE)), 1L)
  expect_identical(.val(rb, "N"), 301)
  # FIML keeps the cases too, and gives no such message.
  m <- testthat::capture_messages(
    reliability_omega(d, missing = "fiml", aux = "z",
                      denominator = "model_implied", ci_method = "ml"))
  expect_false(any(grepl(msg, m, fixed = TRUE)))
})

test_that("when the two-stage fit cannot be computed again, the family says what that means for the coefficient", {
  .ts_skip()
  # The path for a lavaan that does not keep its sample statistics where
  # .ts_refit_n() sets them, simulated. fa_confirmatory()'s warning speaks of its own
  # table (its fit indices and criteria); the family gives its own, once,
  # in its place (through the first versions of 1.1.0 the message that the
  # values count every case came first, and then fa_confirmatory()'s warning said the
  # opposite).
  d <- .ts_verbal()
  d$z <- holzinger_swineford$t1_visual_perception
  d[c(1, 4, 9), paste0("y", 1:5)] <- NA
  ok <- suppressMessages(reliability_omega(d, missing = "two_stage",
                                           aux = "z",
                                           denominator = "model_implied",
                                           ci_method = "ml"))
  ok_a <- suppressMessages(reliability_alpha(d, missing = "two_stage",
                                             aux = "z",
                                             estimator = "model_implied",
                                             ci_method = "ml"))
  local_mocked_bindings(.ts_refit_n = function(fit, n, ...) stop("not here"))
  w <- testthat::capture_warnings(m <- testthat::capture_messages(
    r <- reliability_omega(d, missing = "two_stage", aux = "z",
                           denominator = "model_implied", ci_method = "ml")))
  expect_length(w, 1L)
  expect_match(w, paste0("^The first stage of the two-stage fit read 301 ",
                         "cases, 3 of them observed on an auxiliary variable ",
                         "alone, which lavaan's N of 298 leaves out"))
  expect_match(w, "The coefficient is right, but its analytic standard error",
               fixed = TRUE)
  expect_match(w, "A bootstrap interval counts every case", fixed = TRUE)
  expect_match(m, "the coefficient is estimated from all 301 cases it reads",
               fixed = TRUE, all = FALSE)
  expect_equal(.val(r, "estimate"), .val(ok, "estimate"), tolerance = 1e-10)
  expect_equal(.val(r, "se") / .val(ok, "se"), sqrt(301 / 298),
               tolerance = 1e-6)
  # The model implied alpha gives the family's warning too, not fa_confirmatory()'s,
  # which carries the same class and speaks of fit indices and criteria
  # alpha does not report: the class alone could not tell them apart.
  wa <- testthat::capture_warnings(a <- suppressMessages(
    reliability_alpha(d, missing = "two_stage", aux = "z",
                      estimator = "model_implied", ci_method = "ml")))
  expect_length(wa, 1L)
  expect_match(wa, "^The first stage of the two-stage fit read 301 cases")
  expect_match(wa, "The coefficient is right, but its analytic standard error",
               fixed = TRUE)
  expect_no_match(wa, "fa_confirmatory()", fixed = TRUE)
  expect_equal(.val(a, "estimate"), .val(ok_a, "estimate"), tolerance = 1e-10)
  expect_equal(.val(a, "se") / .val(ok_a, "se"), sqrt(301 / 298),
               tolerance = 1e-6)
  # Analytic alpha, from the first stage alone, has no fit to compute again.
  expect_no_warning(suppressMessages(
    reliability_alpha(d, missing = "two_stage", aux = "z")))
})
