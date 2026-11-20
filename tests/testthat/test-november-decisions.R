# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

test_that("omitted levels yield to tail areas and explicit levels must agree", {
  for (call in list(
    function(...) ci_r(.35, N = 80, ...),
    function(...) ci_R(R = .5, N = 100, p = 3, ...),
    function(...) ci_cv(cv = .2, N = 40, ...)
  )) {
    actual <- call(alpha_lower = .01, alpha_upper = .09)
    expected <- call(conf_level = NULL, alpha_lower = .01, alpha_upper = .09)
    expect_equal(actual, expected, tolerance = 1e-12)
    expect_equal(attr(actual, "conf_level"), .90)
    expect_equal(call(conf_level = .90, alpha_lower = .01, alpha_upper = .09), expected)
    expect_error(call(conf_level = .95, alpha_lower = .01, alpha_upper = .09), "inconsistent")
    expect_error(call(alpha_lower = .01), "both")
    expect_error(call(alpha_lower = NA_real_, alpha_upper = .05), "finite")
  }
  actual <- ci_r(.35, 80, alpha_lower = 0, alpha_upper = .05)
  expect_identical(actual$value[actual$term == "ci_lower"], -1)
  expect_equal(actual$value[actual$term == "ci_upper"],
    tanh(atanh(.35) + qnorm(.95) / sqrt(77)), tolerance = 1e-12)
  values <- seq(2, 5, length.out = 40)
  for (call in list(function(...) ci_cv(data = values, ...),
                   function(...) ci_cv(mean = mean(values), sd = sd(values), N = 40, ...),
                   function(...) ci_cv(cv = sd(values) / mean(values), N = 40, ...))) {
    lower <- call(alpha_lower = 0, alpha_upper = .05)
    upper <- call(alpha_lower = .05, alpha_upper = 0)
    expect_identical(lower$value[lower$term == "ci_lower"], 0)
    expect_identical(upper$value[upper$term == "ci_upper"], Inf)
    expect_true(is.finite(lower$value[lower$term == "ci_upper"]))
    expect_true(is.finite(upper$value[upper$term == "ci_lower"]))
  }
})

test_that("the width proportion is over the studies that return an interval, and the others are counted", {
  counter <- new.env(parent = emptyenv())
  counter$i <- 0L
  local_mocked_bindings(reliability_alpha = function(...) {
    counter$i <- counter$i + 1L
    if (counter$i > 2L) stop("deliberate fitting failure")
    bounds <- if (counter$i == 1L) c(.6, .8) else c(.2, .9)
    data.frame(term = c("estimate", "ci_lower", "ci_upper"),
      value = c(.7, bounds))
  })
  expect_warning(out <- ss_aipe_reliability_sensitivity(true_reliability = .8,
    specified_N = 80, n_items = 4, width = .3, G = 4),
    "describe the 2 replications that returned one")
  expect_equal(out$value[out$term == "pct_ci_less_width"], 1 / 2)
  info <- attr(out, "simulation")
  expect_equal(info$attempted, 4)
  expect_equal(info$returned, 2)
  expect_equal(info$failed, 2)
  expect_equal(info$denominator, 2)
  expect_equal(info$width_target_probability, .5)
  expect_match(info$first_failure, "deliberate")
  expect_equal(out$value[out$term == "mean_ci_width"], .45)
  expect_equal(info$width_target_mc_ci, unname(binom.test(1, 2)$conf.int))
  expect_output(print(out),
    "Replications: 4 attempted, 2 with an interval, 2 without", fixed = TRUE)
})

test_that("no returned intervals means every width and coverage summary is unavailable", {
  local_mocked_bindings(reliability_alpha = function(...) stop("deliberate fitting failure"))
  expect_warning(out <- ss_aipe_reliability_sensitivity(true_reliability = .8,
    specified_N = 80, n_items = 4, width = .3, G = 3),
    "3 of 3 replications returned no confidence interval")
  expect_true(is.na(out$value[out$term == "pct_ci_less_width"]))
  expect_true(is.na(out$value[out$term == "mean_ci_width"]))
  expect_equal(attr(out, "simulation")$failed, 3)
  expect_identical(attr(out, "simulation")$returned, 0L)
  expect_true(is.na(attr(out, "simulation")$width_target_probability))
})

test_that("sentences pair simultaneous intervals with their adjusted tests", {
  x <- data.frame(estimate = .2, ci_lower = -.1, ci_upper = .5,
    p_value = .02, p_adjusted = .08)
  attr(x, "conf_level") <- .95
  attr(x, "simultaneous") <- list(method = "Bonferroni", n_comparisons = 4L)
  attr(x, "p_adjustment") <- "bonferroni"
  expect_match(results_sentence(x), "multiplicity-adjusted p = 0.0800", fixed = TRUE)
  expect_match(results_sentence(x), "simultaneous across the 4 comparisons", fixed = TRUE)
  x$p_adjusted <- NULL
  expect_match(results_sentence(x), "no matching multiplicity-adjusted p-value", fixed = TRUE)
  expect_false(grepl("0.0200", results_sentence(x), fixed = TRUE))
  attr(x, "simultaneous") <- NULL
  expect_match(results_sentence(x), "unadjusted p = 0.0200", fixed = TRUE)
})

test_that("equivalence sentences identify the TOST p-value", {
  x <- equivalence_r(r = .05, N = 100, cor_lower = .2, cor_upper = .2)
  expect_match(results_sentence(x), "TOST p", fixed = TRUE)
  expect_match(results_sentence(generics::glance(x)), "TOST p", fixed = TRUE)
})

test_that("positive composite variance does not license an indefinite covariance matrix", {
  S <- diag(3)
  S[1, 2] <- S[2, 1] <- 2
  expect_gt(sum(S), 0)
  expect_error(reliability_alpha(S = S, N = 100, ci_method = "none"), "negative variance")
  expect_error(reliability_omega(S = S, N = 100, ci_method = "none"), "negative variance")
  for (scales in list(rep(1e-5, 3), c(1e-5, 1e5, 1))) {
    scaled <- S * outer(scales, scales)
    expect_error(reliability_alpha(S = scaled, N = 100, ci_method = "none"), "negative variance")
  }
  psd <- matrix(1, 3, 3)
  expect_null(DMAR:::.relia_negative_eigenvalue(psd))
  valid <- matrix(.4, 3, 3)
  diag(valid) <- 1
  expect_equal(reliability_alpha(S = valid, N = 100, ci_method = "none")$value[1],
    3 / 2 * (1 - sum(diag(valid)) / sum(valid)), tolerance = 1e-12)
})

test_that("a fixed invalid sensitivity specification is refused before simulation", {
  base <- list(true_reliability = .8, specified_N = 80, n_items = 4, width = .3, G = 4)
  for (bad in list(list(G = 0), list(G = 2.5), list(G = NA_real_),
    list(conf_level = 1.2), list(width = 0), list(ci_method = "none"),
    list(ci_method = "mlr"), list(ci_method = "not_a_method"))) {
    expect_error(do.call(ss_aipe_reliability_sensitivity, utils::modifyList(base, bad)))
  }
})

test_that("Bonferroni intervals never borrow Holm or unadjusted p-values", {
  set.seed(113)
  data <- matrix(rnorm(90), 30, 3)
  colnames(data) <- c("A", "B", "C")
  for (adjust in c("holm", "none", "BH")) {
    x <- pairwise_within(data, adjust = adjust, bonferroni_ci = TRUE)
    expect_match(results_sentence(x[1, ]), "no matching multiplicity-adjusted", fixed = TRUE)
    expect_match(results_sentence(generics::tidy(x)[1, ]), "no matching multiplicity-adjusted", fixed = TRUE)
  }
  x <- pairwise_within(data, adjust = "bonferroni", bonferroni_ci = TRUE)
  for (row in list(x[1, ], generics::tidy(x)[1, ])) {
    expect_match(results_sentence(row), paste0("multiplicity-adjusted p = ", format_p(x$p_adjusted[1])), fixed = TRUE)
  }
})

test_that("CFA uncertainty follows the estimator and preserves FIML", {
  skip_if_not_installed("lavaan")
  data <- lavaan::HolzingerSwineford1939
  data$x1[seq_len(30)] <- NA_real_
  factors <- list(v = c("x1", "x2", "x3"), t = c("x4", "x5", "x6"), s = c("x7", "x8", "x9"))
  fit <- fa_confirmatory(data = data, factors = factors, estimator = "MLR", missing = "fiml", output = "fit")
  options <- lavaan::lavInspect(fit, "options")
  expect_identical(options$se, "robust.huber.white")
  expect_identical(options$missing, "ml")
  expect_equal(sum(lavaan::lavInspect(fit, "nobs")), nrow(data))
  ml <- fa_confirmatory(data = data, factors = factors, estimator = "ML", missing = "fiml", output = "fit")
  expect_identical(lavaan::lavInspect(ml, "options")$se, "standard")
  expect_equal(lavaan::coef(fit), lavaan::coef(ml), tolerance = 1e-6)
  explicit <- fa_confirmatory(data = data, factors = factors, estimator = "MLR", se = "standard", missing = "fiml", output = "fit")
  expect_identical(lavaan::lavInspect(explicit, "options")$se, "standard")
  expect_error(fa_confirmatory(data = data, factors = factors, estimator = "MLR", se = "robust.sem", missing = "fiml"), "robust.huber.white", fixed = TRUE)
  expect_error(suppressWarnings(fa_confirmatory(data = data, factors = factors,
    estimator = "ML", missing = "fiml", test = "satorra.bentler")),
    "changed missing-data handling to listwise", fixed = TRUE)
})

test_that("CFA RMSEA intervals use the requested rmsea_conf_level", {
  skip_if_not_installed("lavaan")
  # Ken's ruling of November 8, 2026: every RMSEA interval defaults to 90
  # percent, with an argument to change it. The CFA family's RMSEA interval
  # takes rmsea_conf_level (0.90 by default, the level DMAR 1.0.0 reported),
  # and conf_level, the level of the parameter intervals, no longer reaches
  # it, which withdraws the 1.1.0 development change that made it follow
  # conf_level.
  data <- lavaan::HolzingerSwineford1939
  factors <- list(f = paste0("x", 1:6))
  x <- fa_confirmatory(data, factors, rmsea_conf_level = .99)
  fit <- fa_confirmatory(data, factors, output = "fit")
  ref <- lavaan::fitMeasures(fit, c("rmsea.ci.lower", "rmsea.ci.upper"),
    fm.args = list(rmsea.ci.level = .99))
  expect_equal(x$estimate[x$term == "rmsea_ci_level"], .99)
  expect_equal(x$estimate[x$term %in% c("rmsea_ci_lower", "rmsea_ci_upper")],
               unname(ref), tolerance = 1e-12)
  y <- fa_confirmatory(data, factors, conf_level = .99)
  ref_90 <- lavaan::fitMeasures(fit, c("rmsea.ci.lower", "rmsea.ci.upper"))
  expect_equal(y$estimate[y$term == "rmsea_ci_level"], .90)
  expect_equal(y$estimate[y$term %in% c("rmsea_ci_lower", "rmsea_ci_upper")],
               unname(ref_90), tolerance = 1e-12)
})

test_that("saturated fit conventions are explicit and preserve parameter uncertainty", {
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa("f =~ x1 + x2 + x3", data = lavaan::HolzingerSwineford1939,
    std.lv = TRUE, estimator = "MLR")
  expect_warning(x <- sem_fit_measures(fit), "saturated")
  for (term in c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper")) {
    expect_equal(x$value[x$term == term], 0)
  }
  expect_equal(x$value[x$term == "tli"], 1)
  expect_true(any(lavaan::parameterEstimates(fit)$se > 0))
  expect_true(any(grepl("Parameter standard errors", utils::capture.output(print(x)), fixed = TRUE)))
})

test_that("the shared bootstrap default withholds minority-based intervals", {
  skip_if_not_installed("boot")
  data <- data.frame(id = seq_len(20), y = seq_len(20)^2)
  point <- function(d) {
    if (!all(1:4 %in% d$id)) return(NA_real_)
    mean(d$y)
  }
  expect_warning(x <- DMAR:::.bootstrap_ci(data, point, B = 200, conf_level = .95,
    kind = "percentile", seed = 113), "fewer than half")
  expect_gt(x$bootstrap$returned, 1)
  expect_lt(x$bootstrap$returned, 100)
  expect_true(is.na(x$lower))
  expect_true(is.na(x$se))
  expect_true(x$bootstrap$require_half)
})

test_that("categorical MI forwards and reports 20 inner imputations by default", {
  skip_if_not_installed("lavaan")
  d <- data.frame(y1 = rep(1:4, 5), y2 = rep(c(2, 4, 1, 3), 5),
                  y3 = rep(c(3, 1, 4, 2), 5))
  d$y3[1] <- NA
  seen <- NULL
  local_mocked_bindings(.omega_c_mi = function(data, items, impute,
      imputations, denominator, kind, B, conf_level, seed, ncpus,
      progress = NULL) {
    seen <<- imputations
    list(estimate = .8, se = NA_real_, lower = NA_real_, upper = NA_real_,
         ci_method = "none", bootstrap = NULL, notes = character(),
         undetermined = data.frame(), residual_cor = NULL,
         mi = list(imputations_replicate = imputations, imputations_point = imputations,
                   imputations_returned = imputations))
  })
  expect_message(x <- reliability_omega_categorical(d, missing = "mi"),
                 "each imputed 20 times")
  expect_identical(seen, 20L)
  expect_identical(attr(x, "imputation")$imputations_replicate, 20L)
})

test_that("Parallel population checks are invariant to a common change of units", {
  S <- matrix(.4, 4, 4)
  diag(S) <- 1
  args <- list(model = "Parallel", type = "Factor Analytic", width = .1)
  ref <- do.call(ss_aipe_reliability, c(args, list(true_cov_matrix = S)))
  expect_identical(do.call(ss_aipe_reliability, c(args, list(true_cov_matrix = S * 1e-16)))$value,
                   ref$value)
  unequal <- S
  diag(unequal) <- c(1, 2, 3, 4)
  expect_error(do.call(ss_aipe_reliability, c(args, list(true_cov_matrix = unequal))),
               "not Parallel.*True Score.*Congeneric")
  expect_error(do.call(ss_aipe_reliability, c(args, list(true_cov_matrix = unequal * 1e-16))),
               "not Parallel")
  unequal <- S
  unequal[1, 2] <- unequal[2, 1] <- .45
  expect_error(do.call(ss_aipe_reliability, c(args, list(true_cov_matrix = unequal))),
               "not Parallel")
  noisy <- S
  noisy[1, 2] <- noisy[2, 1] <- .4 + 1e-10
  expect_equal(do.call(ss_aipe_reliability, c(args, list(true_cov_matrix = noisy)))$value,
               ref$value)
})

test_that("congeneric feasibility floors precede missing-information inflation", {
  skip_if_not_installed("lavaan")
  args <- list(model = "Congeneric", type = "Factor Analytic", width = .2,
               n_items = 4, true_factor_loading = rep(.99, 4), true_residual_var = rep(1 - .99^2, 4))
  expect_warning(x <- do.call(ss_aipe_reliability, args), "only a feasibility")
  expect_identical(x$value, 5)
  expect_equal(attr(x, "feasibility_floor"), 5)
  expect_warning(y <- do.call(ss_aipe_reliability, c(args, list(fmi = .5))),
                 "only a feasibility")
  expect_identical(y$value, 9)
  expect_identical(attr(y, "necessary_N_complete"), 5)
})

test_that("Welch simple effects retain the test and withhold unsupported eta inference", {
  z <- as.numeric(scale(seq_len(1000)))
  d <- expand.grid(A = c("low", "high"), B = c("one", "two"), id = seq_len(1000))
  d$y <- ifelse(d$A == "low", z[d$id], 4 + 10 * z[d$id])
  expect_warning(x <- simple_effects_AB(lm(y ~ A * B, d), which = "A_at_B",
                                        error_term = "welch"), "unavailable")
  ref <- oneway.test(y ~ A, d[d$B == "one", ], var.equal = FALSE)
  expect_equal(x$F_value, rep(unname(ref$statistic), 2), tolerance = 1e-12)
  expect_equal(x$F_value, rep(16 / (101 / 1000), 2), tolerance = 1e-12)
  expect_equal(x$df_error, rep((101 / 1000)^2 /
    ((1 / 1000)^2 / 999 + (.1)^2 / 999), 2), tolerance = 1e-10)
  expect_equal(x$p_value, rep(ref$p.value, 2), tolerance = 1e-12)
  expect_true(all(is.na(x$eta_squared_partial)))
  expect_true(all(is.na(x$ci_lower)))
  expect_true(all(is.na(x$ci_upper)))
  expect_output(print(x), "does not establish the ordinary variance-ratio estimand")
})

test_that("likelihood collapse uses the full support and every overlapping row", {
  x <- rbind(c(-1, 0), c(1, 0), c(-1000, NA), c(1000, NA),
              c(NA, -sqrt(2)), c(NA, sqrt(2)))
  colnames(x) <- c("x", "y")
  # Independent profile: with q = 1 - rho^2, optimized variances are
  # a = (10^6 + 1/q)/2 and b = 1. The finite minimum is q = 10^-6;
  # the profile diverges upward as q approaches zero.
  profile <- function(q) 4 * log((1e6 + 1 / q) / 2) + 2 * log(q)
  expect_lt(profile(1e-6), profile(1))
  expect_lt(profile(1e-6), profile(1e-12))
  expect_null(DMAR:::.fmi_unbounded_subset(x))
  tiny <- x
  tiny[1:2, 1] <- c(-1e-15, 1e-15)
  tiny[3:4, 1] <- c(-1, 1)
  expect_null(DMAR:::.fmi_unbounded_subset(tiny))
  S <- matrix(c(1e6, sqrt(999999), sqrt(999999), 1), 2)
  expect_no_error(DMAR:::.fmi_check_saturated(list(S), list(x)))
  expect_null(DMAR:::.fmi_unbounded_subset(sweep(x[, 2:1], 2, c(1e-8, 1e8), "*")))
  y <- cbind(x = 1:5, y = 6 - (1:5), z = c(2, 5, 1, 3, 4))
  why <- DMAR:::.fmi_unbounded_subset(y)
  expect_setequal(why$vars, c("x", "y"))
  expect_identical(why$n, 5L)
  # Recollecting an additional partially observed row disproves that relation.
  y <- rbind(y, c(7, 4, NA))
  expect_null(DMAR:::.fmi_unbounded_subset(y))
  set.seed(113)
  y <- matrix(rnorm(24), 3, 8, dimnames = list(NULL, paste0("y", 1:8)))
  why <- DMAR:::.fmi_unbounded_subset(y)
  expect_true(why$certified_by_dimension)
  expect_error(DMAR:::.stop_fiml_unbounded(y), class = "dmar_fiml_unbounded")
})


test_that("a vector bootstrap applies the half rule separately to each target", {
  skip_if_not_installed("boot")
  d <- data.frame(id = 1:20, y = (1:20)^2)
  point <- function(d) {
    values <- c(all = mean(d$y),
                 rare = if (all(1:4 %in% d$id)) mean(d$y) else NA_real_)
    attr(values, "dmar_undefined") <- c(FALSE, !is.finite(values[2]))
    values
  }
  expect_warning(x <- DMAR:::.bootstrap_ci(d, point, B = 200,
    conf_level = .95, kind = "percentile", seed = 113,
    stat_names = c("all", "rare")), "Intervals withheld for rare")
  expect_true(is.finite(x$lower[[1L]]))
  expect_true(is.na(x$lower[[2L]]))
  expect_identical(unname(x$bootstrap$target_returned["all"]), 200)
  expect_lt(x$bootstrap$target_returned[["rare"]], 100)
})

test_that("a failed H1 diagnostic does not discard valid constrained FIML uncertainty", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  f <- rnorm(60)
  d <- as.data.frame(sapply(seq_len(6), function(j) .7 * f + rnorm(60, sd = sqrt(.51))))
  names(d) <- paste0("y", 1:6)
  d[cbind(5:60, rep(1:6, length.out = 56))] <- NA_real_
  expect_equal(sum(complete.cases(d)), 4)
  model <- paste0("f =~ ", paste0("l", 1:6, "*y", 1:6, collapse = " + "),
    "\n", paste0("y", 1:6, " ~~ p", 1:6, "*y", 1:6, collapse = "\n"),
    "\nomega := (", paste0("l", 1:6, collapse = "+"), ")^2 / ((",
    paste0("l", 1:6, collapse = "+"), ")^2 + ", paste0("p", 1:6, collapse = "+"), ")",
    "\nH := (", paste0("l", 1:6, "^2/p", 1:6, collapse = "+"), ") / (1 + ",
    paste0("l", 1:6, "^2/p", 1:6, collapse = "+"), ")")
  ref <- lavaan::cfa(model, data = d, std.lv = TRUE, missing = "ml")
  pe <- lavaan::parameterEstimates(ref, ci = TRUE, level = .95)
  for (target in c("omega", "H")) {
    expect_warning(x <- if (target == "omega") {
      reliability_omega(d, missing = "fiml", denominator = "model_implied", ci_method = "ml")
    } else reliability_H(data = d, missing = "fiml", ci_method = "ml"),
      class = "dmar_h1_unavailable")
    expected <- pe[pe$lhs == target & pe$op == ":=", ]
    expect_equal(x$value[x$term == "estimate"], expected$est, tolerance = 1e-6)
    expect_equal(x$value[x$term == "se"], expected$se, tolerance = 1e-6)
    expect_equal(x$value[x$term == "ci_lower"], expected$ci.lower, tolerance = 1e-6)
    expect_equal(x$value[x$term == "ci_upper"], expected$ci.upper, tolerance = 1e-6)
    expect_null(attr(x, "residual_cor"))
    expect_match(attr(x, "h1_diagnostic"), "constrained model")
  }
  expect_error(reliability_alpha(d, missing = "fiml", ci_method = "ml"),
               class = "dmar_fiml_unbounded")
  expect_error(reliability_omega(d, missing = "fiml", denominator = "observed"),
               class = "dmar_fiml_unbounded")
  expect_error(sem_fit_measures(ref), "Fit indices requiring the unrestricted H1")
  expect_warning(h <- reliability_H(fit = ref, factor = "f", ci_method = "delta"),
                 class = "dmar_h1_unavailable")
  expect_equal(h$value[h$term == "estimate"], pe$est[pe$lhs == "H" & pe$op == ":="],
               tolerance = 1e-6)
  expect_warning(cfa <- fa_confirmatory(d, list(f = names(d)), missing = "fiml"),
                 class = "dmar_h1_unavailable")
  expect_true(all(is.na(cfa$estimate[cfa$term %in% c("chi_square", "rmsea", "cfi", "H1")])))
  expect_true(any(is.finite(cfa$se)))
  expect_error(fa_confirmatory(d, list(f = names(d)), missing = "fiml", output = "summary"),
               "Use output = \"verbose\"", fixed = TRUE)
  expect_warning(returned_fit <- fa_confirmatory(d, list(f = names(d)), missing = "fiml",
                                       output = "fit"), class = "dmar_h1_unavailable")
  expect_match(attr(returned_fit, "h1_diagnostic"), "H0-based uncertainty")
  expect_match(generics::tidy(cfa)$h1_diagnostic[1L], "H0-based uncertainty")
  expect_match(generics::glance(cfa)$h1_diagnostic, "H0-based uncertainty")
  expect_match(generics::tidy(h)$h1_diagnostic[1L], "H0-based uncertainty")
  expect_match(generics::glance(h)$h1_diagnostic, "H0-based uncertainty")
})

test_that("linked forms retain identified H0 estimates while identifying unavailable H1 pairs", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  f <- rnorm(400)
  d <- as.data.frame(sapply(1:6, function(j) .7 * f + rnorm(400, sd = sqrt(.51))))
  names(d) <- paste0("y", 1:6)
  d[1:200, 1] <- NA; d[201:400, 6] <- NA
  ref <- suppressWarnings(lavaan::cfa("f =~ y1+y2+y3+y4+y5+y6", d,
                                     missing = "ml", std.lv = TRUE))
  pe <- lavaan::parameterEstimates(ref)
  warnings <- testthat::capture_warnings(x <- fa_confirmatory(d, list(f = names(d)), missing = "fiml"))
  expect_length(warnings, 2L)
  expect_true(any(grepl("zero coverage", warnings, fixed = TRUE)))
  expect_match(attr(x, "h1_diagnostic"), "y1 and y6")
  loadings <- x[x$term %in% paste0("loading_", names(d)), ]
  if (!nrow(loadings)) loadings <- x[grepl("=~", x$syntax, fixed = TRUE), ]
  expect_equal(loadings$estimate, pe$est[pe$op == "=~"], tolerance = 1e-6)
  expect_true(all(is.na(x$estimate[x$term %in% c("chi_square", "cfi", "rmsea")])) )
  expect_error(sem_fit_measures(ref), "no case observes both members")
  quiet_coverage <- function(expr) withCallingHandlers(expr,
    warning = function(w) {
      if (grepl("zero coverage|less than 10%.*coverage", conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
  expect_warning(raw_h <- quiet_coverage(reliability_H(data = d, missing = "fiml",
    ci_method = "ml")), class = "dmar_h1_unavailable")
  expect_warning(fit_h <- reliability_H(fit = ref, factor = "f",
    ci_method = "delta"), class = "dmar_h1_unavailable")
  expect_equal(raw_h$value[raw_h$term == "estimate"],
                fit_h$value[fit_h$term == "estimate"], tolerance = 1e-6)
  expect_warning(raw_omega <- quiet_coverage(reliability_omega(d, missing = "fiml",
    denominator = "model_implied", ci_method = "ml")),
    class = "dmar_h1_unavailable")
  expect_true(all(is.finite(raw_omega$value[raw_omega$term %in%
    c("estimate", "ci_lower", "ci_upper")])))
  expect_match(generics::glance(raw_h)$h1_diagnostic,
                "no case observes both members", fixed = TRUE)

  expect_no_error(suppressWarnings(sem_information_criteria(ref, measures = c("AIC", "BIC"))))
  expect_equal(x$estimate[x$term == "AIC"], -2 * as.numeric(lavaan::logLik(ref)) +
                 2 * attr(lavaan::logLik(ref), "df"))
  for (extra in list(list(h1.information.meat = "unstructured"),
                     list(information.bread = "expected"))) {
    expect_error(suppressWarnings(do.call(fa_confirmatory, c(list(data = d,
      factors = list(f = names(d)), missing = "fiml", estimator = "MLR"), extra))),
      "custom SE/information combination")
  }
})

test_that("two-stage H1 diagnostics include named auxiliary columns", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  f <- matrix(rnorm(800), 400, 2)
  d <- as.data.frame(sapply(1:6, function(j) .8 * f[, if (j <= 3) 1 else 2] +
                             .6 * rnorm(400)))
  names(d) <- paste0("y", 1:6)
  d$a <- d$y1 + d$y2
  d$y1[seq(2, 400, 3)] <- NA
  model <- "f1 =~ y1+y2+y3\nf2 =~ y4+y5+y6"
  fit <- suppressWarnings(lavaan::cfa(model, d, missing = "two.stage", aux = "a"))
  why <- DMAR:::.fiml_fit_h1_diagnostic(fit)
  expect_match(why, "y1, y2, and a", fixed = TRUE)
  expect_error(average_variance_extracted(fit), "requires eligible unrestricted moments")
  expect_error(reliability_H(fit = fit, factor = "f1", ci_method = "none"),
               "requires unrestricted first-stage moments")
  expect_error(suppressWarnings(sem_fit_measures(fit)), "unrestricted H1")
  expect_error(suppressWarnings(sem_information_criteria(fit, measures = "AIC")),
               "unrestricted first-stage moments")
  expect_error(suppressWarnings(fa_confirmatory(d, list(f1 = paste0("y", 1:3),
    f2 = paste0("y", 4:6)), missing = "two_stage", aux = "a")),
    "requires unrestricted first-stage moments")
})


test_that("interval availability follows actual bootstrap limits", {
  d <- data.frame(y = seq(-3, -1, length.out = 40))
  x <- .bootstrap_ci(d, function(x) mean(x$y), B = 40,
    conf_level = .95, kind = "bootstrap_se_logistic", seed = 113)
  expect_true(is.finite(x$se))
  expect_true(is.na(x$lower) && is.na(x$upper))
  expect_true(unname(x$bootstrap$target_withheld))
})


test_that("omega sensitivity rejects a structurally unidentified design", {
  expect_error(ss_aipe_reliability_sensitivity(true_reliability = .8,
    n_items = 2, width = .2, specified_N = 100, estimator = "omega", G = 1),
    "at least three items", fixed = TRUE)
})


test_that("bootstrap refits muffle repeated coverage notices but retain other warnings", {
  d <- data.frame(y = 1:20)
  fn <- function(d) {
    warning("lavaan: some pairwise combinations have zero coverage", call. = FALSE)
    mean(d$y)
  }
  expect_no_warning(.boot_replicate(fn, d, "keep"))
  other <- function(d) {
    fn(d)
    warning("an independent fitting diagnostic", call. = FALSE)
    mean(d$y)
  }
  expect_warning(.boot_replicate(other, d, "keep"), "independent fitting diagnostic")
})

test_that("the 1.0.0 confirmatory functions keep their argument order up to the dots", {
  # rmsea_conf_level (1.1.0) comes after ..., beside wald_df, so a 1.0.0
  # call that gives its arguments by position still reaches the same ones.
  before_dots <- function(f) {
    n <- names(formals(f))
    n[seq_len(match("...", n))]
  }
  expect_identical(before_dots(cfa_k),
                   c("data", "factors", "S", "N", "M", "equal_loading",
                     "equal_intercept", "equal_error", "correlated_factors",
                     "meanstructure", "estimator", "missing", "ordered", "se",
                     "conf_level", "output", "..."))
  expect_identical(before_dots(cfa_1),
                   c("data", "items", "S", "N", "equal_loading",
                     "equal_error", "estimator", "missing", "se",
                     "conf_level", "output", "..."))
  expect_identical(before_dots(cfa_2),
                   c("data", "factor_1", "factor_2", "S", "N",
                     "equal_loading", "equal_error", "correlated_factors",
                     "estimator", "missing", "se", "conf_level", "output",
                     "..."))
  for (f in list(cfa_k, cfa_1, cfa_2)) {
    expect_identical(formals(f)$rmsea_conf_level, 0.90)
  }
})
