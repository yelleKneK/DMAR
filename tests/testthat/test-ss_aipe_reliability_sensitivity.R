# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

## ss_aipe_reliability_sensitivity(): planning under missing data and the
## omega estimator. The family-wide schema is asserted in
## test-ss_aipe_sensitivity_family.R; these tests cover what this member
## adds. Monte Carlo, so skip_on_cran().

# The rows of this member's table, in order, with no assurance supplied.
.sens_terms <- c("mean_reliability", "median_reliability", "sd_reliability",
                 "mean_ci_width", "median_ci_width", "sd_ci_width",
                 "pct_ci_less_width", "pct_ci_miss_low", "pct_ci_miss_high",
                 "total_type_I_error", "total_N", "n_items",
                 "true_reliability", "supposed_reliability", "width",
                 "conf_level")

.sens_call <- function(...) {
  ss_aipe_reliability_sensitivity(true_reliability = .8,
                                  supposed_reliability = .8, n_items = 4,
                                  width = .15, ...)
}
.expect_sens_core <- function(res) {
  core <- DMAR:::.SS_AIPE_SENS_CORE_TERMS
  expect_true(all(core %in% res$term))
  expect_type(res$value, "double")
  v <- stats::setNames(res$value, res$term)
  props <- v[c("pct_ci_less_width", "pct_ci_miss_low", "pct_ci_miss_high",
               "total_type_I_error")]
  expect_false(anyNA(props))
  expect_true(all(props >= 0 & props <= 1))
  expect_equal(unname(v["total_type_I_error"]),
               unname(v["pct_ci_miss_low"] + v["pct_ci_miss_high"]),
               tolerance = 1e-12)
  expect_true(is.finite(v["mean_ci_width"]) && v["mean_ci_width"] > 0)
}

test_that("the complete data analysis keeps 1.0.x's numbers for the interval 1.0.x used", {
  skip_on_cran()
  # Values of the 1.0.x function at this call (set.seed(113), G = 20),
  # whose default interval for alpha was Bonett's; naming it reproduces
  # them, and fmi = 0 with missing_prop = 0 changes nothing.
  set.seed(113)
  res <- .sens_call(G = 20, ci_method = "bonett")
  set.seed(113)
  zero <- .sens_call(G = 20, ci_method = "bonett", fmi = 0,
                     missing_prop = 0)
  expect_identical(zero, res)
  expect_null(attr(res, "fmi"))
  expect_null(attr(res, "missing_prop"))
  expect_identical(res$term, .sens_terms)
  v <- stats::setNames(res$value, res$term)
  expect_identical(unname(v["total_N"]), 74)
  expect_equal(unname(v[c("mean_reliability", "mean_ci_width",
                          "pct_ci_less_width", "total_type_I_error")]),
               c(0.781799204193218, 0.168539782754444, 0.4, 0.05),
               tolerance = 1e-12)
  # The default interval is now "ml", the one the plan is made for: the
  # same samples, the same estimates, and narrower intervals.
  set.seed(113)
  ml <- .sens_call(G = 20)
  set.seed(113)
  expect_identical(.sens_call(G = 20, ci_method = "ml"), ml)
  w <- stats::setNames(ml$value, ml$term)
  expect_identical(unname(w["total_N"]), 74)
  expect_equal(unname(w[c("mean_reliability", "mean_ci_width",
                          "pct_ci_less_width", "total_type_I_error")]),
               c(0.781799204193218, 0.163036656351598, 0.45, 0.05),
               tolerance = 1e-12)
})

test_that("the replications evaluate the interval the plan was made for, which meets its assurance", {
  skip_on_cran()
  # The 1.1.0 review: the plan is ss_aipe_reliability()'s for the parallel
  # model with the normal theory interval, but the replications computed
  # Bonett's interval for alpha (and mlr for omega), so with the true
  # reliability equal to the planning value a correctly specified
  # assurance .80 plan fell short of .80. The default is now "ml", the
  # interval the plan describes.
  set.seed(113)
  res <- .sens_call(G = 2000, assurance = .8)
  v <- stats::setNames(res$value, res$term)
  expect_identical(unname(v["total_N"]), 99)
  expect_equal(unname(v["pct_ci_less_width"]), 0.813, tolerance = 1e-12)
  expect_gte(unname(v["pct_ci_less_width"]), .8)
  expect_equal(unname(v["mean_ci_width"]), 0.131160558183130,
               tolerance = 1e-12)
  set.seed(113)
  bonett <- .sens_call(G = 2000, assurance = .8, ci_method = "bonett")
  b <- stats::setNames(bonett$value, bonett$term)
  expect_identical(unname(b["total_N"]), 99)
  expect_equal(unname(b["pct_ci_less_width"]), 0.7675, tolerance = 1e-12)
  expect_lt(unname(b["pct_ci_less_width"]), .8)
  # At a population of parallel items the variance the plan takes, the
  # compound symmetry form of van Zyl, Neudecker, and Nel (2000), is the
  # variance of reliability_alpha(ci_method = "ml") at the population and
  # the delta method variance of the model implied omega, which
  # reliability_omega(denominator = "model_implied", ci_method = "ml")
  # uses.
  u <- .8 / (4 * .2 + .8)
  S <- matrix(u, 4, 4)
  diag(S) <- 1
  cs <- 2 * (1 - .8)^2 * 4 / 3
  expect_equal(DMAR:::.ci_alpha_ml(.8, S = S, N = 101, J = 4,
                                   conf_level = .95)$se^2 * 100, cs,
               tolerance = 1e-12)
  expect_equal(DMAR:::.ss_aipe_omega_avar(S, rep(sqrt(u), 4),
                                          rep(1 - u, 4)), cs,
               tolerance = 1e-12)
  # Omega's default is "ml" too.
  skip_if_not_installed("lavaan")
  set.seed(113)
  om <- .sens_call(G = 3, estimator = "omega")
  set.seed(113)
  expect_identical(.sens_call(G = 3, estimator = "omega", ci_method = "ml"),
                   om)
})

test_that("fmi enlarges the planned N exactly as the planner does", {
  skip_on_cran()
  u <- .8 / (4 * .2 + .8)
  plan <- ss_aipe_reliability(model = "parallel", type = "Normal Theory",
                              width = .15, n_items = 4, true_factor_loading = sqrt(u),
                              true_residual_var = 1 - u, fmi = .18)
  set.seed(113)
  res <- suppressWarnings(.sens_call(G = 3, fmi = .18))
  v <- stats::setNames(res$value, res$term)
  expect_identical(unname(v["total_N"]), plan$value)
  expect_gt(unname(v["total_N"]), 74)
  expect_identical(attr(res, "fmi"), .18)
  expect_identical(attr(res, "missing_prop"), 0)
  expect_identical(attr(res, "necessary_N_complete"), 74)
  expect_output(print(res),
                "Planned for a fraction of missing information of 0.18; simulated with complete data.",
                fixed = TRUE)
})

test_that("missing_prop simulates MCAR data analyzed by FIML", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(113)
  res <- suppressWarnings(.sens_call(G = 5, fmi = .18, missing_prop = .2))
  .expect_sens_core(res)
  expect_identical(res$term, .sens_terms)
  expect_identical(attr(res, "missing_prop"), .2)
  expect_output(print(res),
                "simulated with each item score missing completely at random with probability 0.2 and analyzed by full information maximum likelihood.",
                fixed = TRUE)
  # The default interval under missing data is alpha's FIML "ml", the
  # Wald interval the fmi adjustment plans for (it was "ml_logistic" before
  # 1.1.0): the same seed and an explicit ci_method give the same table.
  set.seed(113)
  explicit <- suppressWarnings(.sens_call(G = 5, fmi = .18, missing_prop = .2,
                                          ci_method = "ml"))
  expect_identical(explicit, res)
  # The omega estimator runs under FIML too.
  set.seed(113)
  om <- suppressWarnings(.sens_call(G = 3, missing_prop = .2,
                                    estimator = "omega"))
  .expect_sens_core(om)
  expect_identical(attr(om, "fmi"), 0)
})

test_that("missing values are deleted completely at random at missing_prop", {
  set.seed(113)
  Y <- matrix(stats::rnorm(4000), 1000, 4)
  set.seed(1)
  Ym <- DMAR:::.mcar_delete(Y, .2)
  # 4000 Bernoulli(.2) draws: the standard error of the proportion is .0063.
  expect_lt(abs(mean(is.na(Ym)) - .2), .025)
  expect_identical(Ym[!is.na(Ym)], Y[!is.na(Ym)])
  set.seed(1)
  expect_identical(is.na(Ym), matrix(stats::runif(4000) < .2, 1000, 4))
})

test_that("intervals that ignore missing data, and fmi without a plan, are refused", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  for (m in c("bonett", "feldt", "fisher", "hakstian_whalen", "adf")) {
    expect_error(.sens_call(G = 2, missing_prop = .2, ci_method = m),
                 "'ci_method' must be an interval whose standard error comes from the FIML information matrix",
                 fixed = TRUE)
  }
  expect_error(.sens_call(G = 2, missing_prop = .2, estimator = "omega",
                          ci_method = "adf"),
               "For estimator = \"omega\" use one of", fixed = TRUE)
  expect_error(ss_aipe_reliability_sensitivity(true_reliability = .8,
                                               specified_N = 60, n_items = 4,
                                               width = .15, G = 2, fmi = .2),
               "'fmi' is a planning value", fixed = TRUE)
  for (bad in list(1, -.1, c(.1, .2), NA_real_)) {
    expect_error(.sens_call(G = 2, missing_prop = bad),
                 "'missing_prop' must be a single number in [0, 1)",
                 fixed = TRUE)
    expect_error(.sens_call(G = 2, fmi = bad),
                 "'fmi' must be a single number in [0, 1)", fixed = TRUE)
  }
})

test_that("with specified_N and missing_prop the planning value is recorded as NA", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(113)
  res <- suppressWarnings(ss_aipe_reliability_sensitivity(
    true_reliability = .8, specified_N = 60, n_items = 4, width = .15, G = 3,
    missing_prop = .1))
  expect_true(is.na(attr(res, "fmi")))
  expect_output(print(res), "Evaluated at the specified N; simulated with",
                fixed = TRUE)
})

test_that("the omega estimator returns intervals (it returned none in 1.0.x)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # reliability_omega() defaults to the observed denominator, for which the
  # mlr interval is not defined, so every replication used to fail
  # silently and the table reported no widths and zero non-coverage. The
  # model implied denominator is now requested explicitly.
  set.seed(113)
  res <- .sens_call(G = 5, estimator = "omega")
  .expect_sens_core(res)
})

test_that("failed studies are left out of the proportions and counted, and invalid settings stop before simulation", {
  skip_on_cran()
  # Every replication fails: alpha's analytic estimator has no mlr
  # interval. The function stops with the reason instead of returning a
  # table of missing values.
  expect_error(.sens_call(G = 3, ci_method = "mlr"),
               "'ci_method' must name an interval available for estimator = 'alpha':",
               fixed = TRUE)
  # Every other replication fails: the proportions are over the ones that
  # returned an interval, and the per-replication file records NA for the
  # rest.
  real_alpha <- reliability_alpha
  state <- new.env(parent = emptyenv())
  state$n <- 0L
  local_mocked_bindings(reliability_alpha = function(...) {
    state$n <- state$n + 1L
    if (state$n %% 2L == 0L) stop("simulated failure")
    real_alpha(...)
  })
  out_file <- tempfile(fileext = ".csv")
  on.exit(unlink(out_file), add = TRUE)
  set.seed(113)
  expect_warning(res <- .sens_call(G = 6, filename = out_file),
                 "3 of 6 replications returned no confidence interval. The width, width-target, and coverage rows describe the 3 replications that returned one. The first failure: simulated failure",
                 fixed = TRUE)
  written <- utils::read.csv(out_file)
  expect_identical(which(is.na(written$type_I_lower)), c(2L, 4L, 6L))
  v <- stats::setNames(res$value, res$term)
  ok <- !is.na(written$type_I_lower)
  expect_equal(unname(v["pct_ci_miss_low"]), mean(written$type_I_lower[ok]))
  expect_equal(unname(v["pct_ci_less_width"]),
               mean(written$ci_width[ok] <= .15))
  expect_identical(attr(res, "simulation")[c("attempted", "returned",
                                             "failed", "denominator")],
                   list(attempted = 6L, returned = 3L, failed = 3L,
                        denominator = 3L))
})

test_that("missing_prop records the FMI the simulated data realized for the coefficient", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  out_file <- tempfile(fileext = ".csv")
  on.exit(unlink(out_file), add = TRUE)
  set.seed(113)
  res <- suppressWarnings(.sens_call(G = 4, fmi = .18, missing_prop = .2,
                                     filename = out_file))
  written <- utils::read.csv(out_file)
  expect_named(written, c("rel_hat", "ci_lower", "ci_upper", "ci_width",
                          "type_I_lower", "type_I_upper", "realized_fmi"))
  # Each replication's FMI lies in (0, 1): items missing completely at
  # random with probability .2 cost alpha some, but not all, information.
  expect_true(all(written$realized_fmi > 0 & written$realized_fmi < 1))
  expect_identical(attr(res, "n_realized_fmi"), 4L)
  expect_equal(attr(res, "mean_realized_fmi"), mean(written$realized_fmi),
               tolerance = 1e-12)
  expect_identical(res$term, .sens_terms)
  expect_output(print(res), paste0(
    "The fraction of missing information the simulated data realized for ",
    "the coefficient averaged ",
    formatC(attr(res, "mean_realized_fmi"), digits = 3, format = "fg"),
    " over the 4 replications in which it was computed."), fixed = TRUE)
  # The omega estimator records it too.
  set.seed(113)
  om <- suppressWarnings(.sens_call(G = 2, missing_prop = .2,
                                    estimator = "omega"))
  expect_identical(attr(om, "n_realized_fmi"), 2L)
  expect_gt(attr(om, "mean_realized_fmi"), 0)
  # With complete data there is nothing to realize, and nothing is added.
  set.seed(113)
  expect_null(attr(.sens_call(G = 2), "mean_realized_fmi"))
  set.seed(113)
  expect_null(attr(.sens_call(G = 2, fmi = .18), "mean_realized_fmi"))
})

test_that("the realized FMI is the coefficient's unstructured FMI, and computing it draws nothing", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(113)
  u <- .8 / (4 * .2 + .8)
  n <- 120
  Y <- matrix(sqrt(u) * stats::rnorm(n), n, 4) +
    matrix(stats::rnorm(n * 4, sd = sqrt(1 - u)), n, 4)
  Y <- as.data.frame(DMAR:::.mcar_delete(Y, .25))
  # Two cases that lose every item score. The fit sets them aside, and the
  # realized FMI folds their information back in: it is the FMI relative to
  # all nrow(Y) cases, 1 - (n_used / n)(1 - fmi_used), where fmi_used is the
  # FMI relative to the n_used cases the fit used.
  Y <- rbind(Y, Y[1:2, ] * NA)
  n <- nrow(Y)
  fold <- function(fmi_used, n_used) 1 - (n_used / n) * (1 - fmi_used)
  # Alpha: lavaan's own FMI for the same defined parameter of the saturated
  # model, computed by its own refit under the unstructured information.
  # (In a saturated model the three estimators coincide, Chen & Savalei,
  # 2021; lavaan's default numeric Hessian agrees to about 3e-7 here.) The
  # realized FMI is evaluated at the EM estimates of the saturated model and
  # this fit where lavaan's optimizer stopped, whose covariances differ by up
  # to 2e-6; the two FMIs agree to about 1e-8.
  sat <- paste(c(
    "V1 ~~ a11*V1", "V2 ~~ a22*V2", "V3 ~~ a33*V3", "V4 ~~ a44*V4",
    "V1 ~~ a12*V2", "V1 ~~ a13*V3", "V1 ~~ a14*V4", "V2 ~~ a23*V3",
    "V2 ~~ a24*V4", "V3 ~~ a34*V4",
    paste0("alpha := (4/3) * (1 - (a11 + a22 + a33 + a44) / (a11 + a22 + ",
           "a33 + a44 + 2*(a12 + a13 + a14 + a23 + a24 + a34)))")),
    collapse = "\n")
  # lavaan notes the two empty cases it sets aside. Its EM is run to the
  # tolerance of the refits behind the realized FMI (1e-10, where lavaan's
  # default is 1e-5), so that its own FMI comes from the same moments.
  expect_warning(
    fit <- lavaan::lavaan(sat, data = Y, missing = "ml", meanstructure = TRUE,
                          int.ov.free = TRUE, information = "observed",
                          observed.information = "h1",
                          h1.information = "unstructured",
                          em.h1.args = list(tol = 1e-10, max_iter = 10000L)),
    "some cases are empty", fixed = TRUE)
  pe <- lavaan::parameterEstimates(fit, fmi = TRUE)
  n_used <- lavaan::lavInspect(fit, "nobs")
  expect_lte(n_used, n - 2L)
  expect_equal(DMAR:::.aipe_realized_fmi(Y, "alpha"),
               fold(pe$fmi[pe$lhs == "alpha" & pe$op == ":="], n_used),
               tolerance = 1e-6)
  # Omega: the omega row of fraction_missing_information() for the
  # single-factor FIML fit.
  expect_warning(
    cfit <- fa_confirmatory(data = Y, factors = list(f1 = names(Y)), missing = "fiml",
                  output = "fit"),
    "some cases are empty", fixed = TRUE)
  tab <- fraction_missing_information(cfit)
  expect_equal(attr(tab, "N_excluded"), n - n_used)
  expect_equal(DMAR:::.aipe_realized_fmi(Y, "omega"),
               fold(tab$fmi[tab$term == "omega_f1"], n_used),
               tolerance = 1e-10)
  # A data set whose fit fails gives NA rather than an error (lavaan prints
  # its table of the variables before it stops).
  invisible(utils::capture.output(
    degenerate <- DMAR:::.aipe_realized_fmi(Y[1:2, ], "alpha")))
  expect_identical(degenerate, NA_real_)

  # The step draws no random numbers: the table equals that of a run in
  # which it is skipped, and a run in which no replication returns an FMI
  # says so.
  set.seed(113)
  with_fmi <- suppressWarnings(.sens_call(G = 3, missing_prop = .2))
  local_mocked_bindings(.aipe_realized_fmi = function(...) NA_real_)
  set.seed(113)
  without <- suppressWarnings(.sens_call(G = 3, missing_prop = .2))
  expect_identical(without$value, with_fmi$value)
  expect_identical(attr(without, "n_realized_fmi"), 0L)
  expect_identical(attr(without, "mean_realized_fmi"), NA_real_)
  expect_output(print(without),
                "could not be computed in any replication.", fixed = TRUE)
})

test_that("an assurance plan is searched at the planner's default numbers of iterations", {
  skip_on_cran()
  # 1.0.x searched at 50 and 200 iterations, whose noise
  # let a size below the one the assurance needs pass by chance: for this
  # call, 93 to 101 cases over 30 seeds (mean 96.9) against 99 to 101 at
  # the defaults, where the assurance of .80 is first reached near 99.
  u <- .8 / (4 * .2 + .8)
  set.seed(5)
  plan <- ss_aipe_reliability(model = "parallel", type = "Normal Theory",
                              width = .15, n_items = 4, true_factor_loading = sqrt(u),
                              true_residual_var = 1 - u, assurance = .8)
  set.seed(5)
  res <- .sens_call(G = 2, assurance = .8)
  expect_identical(res$value[res$term == "total_N"],
                   plan$value[plan$term == "necessary_N"])
  expect_identical(res$value[res$term == "assurance"], .8)
})

test_that("at a true reliability of zero a clamped interval covers zero, so no replication misses high", {
  skip_on_cran()
  # The reliability family clamps both limits into [0, 1] (since 1.1.0;
  # through 1.0.1 an alpha interval wholly below zero came back inverted, a
  # lower limit of zero above a negative upper one, and counted here as a
  # miss high). Coverage is counted on the reported interval, and [0, 0]
  # contains zero. With four items, N = 40, and Bonett's interval, 13 of
  # 400 replications lie wholly below zero: 1.0.1 reported a miss rate
  # of .0400, .0325 of it high, and it is now the .0075 missed low.
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f))
  set.seed(1)
  res <- ss_aipe_reliability_sensitivity(true_reliability = 0, n_items = 4,
                                         width = .3, specified_N = 40,
                                         estimator = "alpha",
                                         ci_method = "bonett", G = 400,
                                         filename = f)
  v <- stats::setNames(res$value, res$term)
  per <- utils::read.csv(f)
  expect_identical(sum(per$ci_lower == 0 & per$ci_upper == 0), 13L)
  expect_false(any(per$ci_upper < per$ci_lower))
  expect_identical(unname(v["pct_ci_miss_high"]), 0)
  expect_equal(unname(v[c("pct_ci_miss_low", "total_type_I_error")]),
               c(.0075, .0075), tolerance = 1e-12)
  # Above zero the same kind of interval misses high, as it always did.
  set.seed(1)
  res <- ss_aipe_reliability_sensitivity(true_reliability = .05, n_items = 4,
                                         width = .3, specified_N = 40,
                                         estimator = "alpha",
                                         ci_method = "bonett", G = 400)
  v <- stats::setNames(res$value, res$term)
  expect_equal(unname(v[c("pct_ci_miss_high", "total_type_I_error")]),
               c(.0325, .0450), tolerance = 1e-12)
})
