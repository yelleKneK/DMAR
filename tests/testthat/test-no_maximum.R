# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# A likelihood with no maximum: the check (.lavaan_no_maximum(),
# R/no_maximum_internals.R), the refusal of class dmar_no_maximum it leads
# to on every fit path that runs it, the lavaan code the refusal prints,
# and the count of a bootstrap replication with no maximum among the
# replications that raised an error. The structural equation models come
# last: the SEM planners' fits to the population moments and
# mediation_mbco()'s lavaan fit, which are refused with the user's model
# printed whole; the replications of the planners' Monte Carlo, which are
# set aside; the functions that read only the chi square, which do not
# ask; and the observed variables lavaan carries as latent variables of
# their own, which the check leaves where they are.

skip_if_not_installed("lavaan")

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
caught <- function(expr) {
  tryCatch(quiet(expr), dmar_no_maximum = function(e) e)
}

# Records the warnings of a call, muffled, beside its value.
with_warnings <- function(expr) {
  seen <- new.env(parent = emptyenv())
  seen$msgs <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    seen$msgs <- c(seen$msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = seen$msgs)
}

# The lavaan code of a refusal, run as a user would run it: the model
# syntax without its ":=" lines, which do not change the fit (the message
# prints DMAR's own syntax without them, and a model the user wrote
# whole), and the call, with the placeholders (my_data, my_S, my_M) given
# in `...`. Returns the log likelihood of the fit.
run_code <- function(e, ...) {
  env <- list2env(list(...), parent = environment())
  if (!is.null(e$syntax)) {
    env$model <- paste(grep(":=", strsplit(e$syntax, "\n")[[1L]],
                            fixed = TRUE, invert = TRUE, value = TRUE),
                       collapse = "\n")
  }
  as.numeric(lavaan::logLik(quiet(eval(parse(text = e$lavaan_call),
                                       envir = env))))
}

# The model implied covariance matrix of block `b` of a lavaan parameter
# table `pt` at the values `st`, written out: loadings on the observed
# variables, regressions among the factors (a "=~" whose right side is a
# factor, as a second-order loading is), and the variances and
# covariances of the factors and of the errors.
implied_at <- function(pt, st, b = 1L) {
  in_b <- pt$block == b
  lv <- unique(pt$lhs[in_b & pt$op == "=~"])
  ov <- setdiff(unique(pt$rhs[in_b & pt$op == "=~"]), lv)
  L <- matrix(0, length(ov), length(lv), dimnames = list(ov, lv))
  B <- matrix(0, length(lv), length(lv), dimnames = list(lv, lv))
  P <- B
  Th <- matrix(0, length(ov), length(ov), dimnames = list(ov, ov))
  for (r in which(in_b & pt$op %in% c("=~", "~~"))) {
    l <- pt$lhs[r]
    h <- pt$rhs[r]
    if (pt$op[r] == "=~") {
      if (h %in% lv) B[h, l] <- st[r] else L[h, l] <- st[r]
    } else if (l %in% lv) {
      P[l, h] <- P[h, l] <- st[r]
    } else {
      Th[l, h] <- Th[h, l] <- st[r]
    }
  }
  A <- solve(diag(length(lv)) - B)
  L %*% A %*% P %*% t(A) %*% t(L) + Th
}

# The model of the 2026-11-05 note: spatial, verbal, and speed group
# factors and a fourth on t24, t25, and t26, whose residual covariances
# beyond the general factor are all negative.
f4x <- list(
  spatial = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
  verbal  = c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning"),
  speed   = c("t10_addition", "t12_counting_groups_of_dots",
              "t13_straight_and_curved_capitals"),
  extra   = c("t24_woody_mccall", "t25_paper_form_board_r", "t26_flags"))

# Three items whose correlations have a negative product (the sample of
# test-reliability_omega.R, population r12 = .4, r13 = .3, r23 = -.05):
# no single factor reproduces them, and lavaan stops on the ridge on
# which y1's loading grows without bound.
set.seed(3)
S_neg <- matrix(c(1, .4, .3, .4, 1, -.05, .3, -.05, 1), 3)
triad_neg <- as.data.frame(matrix(rnorm(300 * 3), 300) %*% chol(S_neg))
names(triad_neg) <- paste0("y", 1:3)

# Two groups: the triad in A, three items correlated .5 in B.
set.seed(4)
S_ok <- matrix(.5, 3, 3)
diag(S_ok) <- 1
other <- as.data.frame(matrix(rnorm(300 * 3), 300) %*% chol(S_ok))
names(other) <- paste0("y", 1:3)
two_groups <- rbind(cbind(triad_neg, grp = "A"), cbind(other, grp = "B"))

# Three items correlated .8, .8, and .5: the one factor model is
# saturated, and its unique solution gives the first an error variance of
# 1 - .8 * .8 / .5 = -.28, a finite maximum.
set.seed(113)
R_fin <- matrix(c(1, .8, .8, .8, 1, .5, .8, .5, 1), 3)
triad_fin <- as.data.frame(MASS::mvrnorm(200, rep(0, 3), R_fin,
                                         empirical = TRUE))
names(triad_fin) <- paste0("y", 1:3)

# A stand-in for the check, for the tests of what a caller does with its
# answer: a ridge on `variable` for every fit, or (with `unless`) for
# every fit whose data do not hold the row `unless`, so that the fit to
# the whole sample has a maximum and the resamples that leave the row out
# have none.
fake_ridge <- function(variable = "y2", unless = NULL) {
  function(fit, step = 10) {
    if (!is.null(unless)) {
      X <- lavaan::lavInspect(fit, "data")
      held <- rowSums(abs(sweep(unname(X), 2, unname(unless)))) < 1e-12
      if (any(held)) return(NULL)
    }
    list(variable = variable, group = NULL, estimate = -50, fixed_at = -500,
         loglik = c(fit = -100, further = -99.9), tol = 1e-4, N = 100)
  }
}

test_that("the bifactor model of the note stops under FIML, with its fields and code", {
  skip_on_cran()
  e <- caught(fa_bifactor(holzinger_swineford, f4x))
  expect_s3_class(e, c("dmar_no_maximum", "dmar_refusal", "error",
                       "condition"), exact = TRUE)
  expect_identical(e$code, "no_maximum")
  expect_identical(e$fn, "fa_bifactor")
  expect_identical(e$variables, "t25_paper_form_board_r")
  expect_null(e$group)
  # Where the optimizer stopped (-1163 with lavaan 0.7.3) is set by its
  # starting values and tolerance, not by the data, so only the side of
  # it the variance is on is pinned, with the fields.
  expect_lt(e$estimate, -1000)
  expect_equal(e$fixed_at, 10 * e$estimate)
  # The fit further out has a higher log likelihood: the likelihood is
  # still rising there.
  expect_named(e$loglik, c("fit", "further"))
  expect_gt(e$loglik[["further"]], e$loglik[["fit"]])
  expect_equal(e$N, 301)
  expect_equal(e$tolerance, 1e-6 * 301)
  msg <- conditionMessage(e)
  expect_match(msg, paste0("^fa_bifactor\\(\\) stopped: the fit failed ",
                           "because the likelihood has no maximum\\. The ",
                           "error variance of t25_paper_form_board_r runs ",
                           "off below zero: lavaan's optimizer stopped with ",
                           "it at -[0-9]+, and with it held at -[0-9]+ the ",
                           "log likelihood is higher, "))
  expect_match(msg, "the likelihood is still rising as the variance falls",
               fixed = TRUE)
  expect_match(msg, paste("DMAR fits the model with lavaan, the estimation",
                          "engine; see lavaan for more details."),
               fixed = TRUE)
  expect_false(grepl("\u2014", msg))
  # The choices in lavaan come before the code, which comes last, so that
  # nothing R would cut from a long message is lost but code; and the
  # message says where the code is kept.
  expect_match(msg, paste("There, bounds = \"pos.var\" keeps every error",
                          "variance at zero or above"), fixed = TRUE)
  expect_lt(regexpr("bounds = \"pos.var\"", msg, fixed = TRUE),
            regexpr("model <- '", msg, fixed = TRUE))
  expect_match(msg, "the error carries it as its syntax and lavaan_call elements",
               fixed = TRUE)
  # The code: DMAR's syntax (whole in the condition, without the lines
  # that define its indices in the message) and the call, the data a
  # placeholder.
  expect_match(e$syntax, "extra =~ NA*t24_woody_mccall", fixed = TRUE)
  expect_match(e$syntax, "omega_H_general :=", fixed = TRUE)
  expect_false(grepl(":=", msg, fixed = TRUE))
  expect_match(msg, "model <- '\nextra|model <- '\ngeneral")
  expect_match(e$lavaan_call,
               "^lavaan::lavaan\\(model, data = my_data, estimator = \"ML\"")
  expect_match(e$lavaan_call, "missing = \"fiml\")", fixed = TRUE)
  expect_true(endsWith(msg, e$lavaan_call))
})

test_that("listwise deletion sends two error variances off together, and output = \"fit\" stops", {
  skip_on_cran()
  e <- caught(fa_bifactor(holzinger_swineford, f4x, missing = "listwise",
                          output = "fit"))
  expect_s3_class(e, "dmar_no_maximum")
  expect_setequal(e$variables, c("t2_cubes", "t25_paper_form_board_r"))
  expect_length(e$estimate, 2L)
  expect_true(all(e$estimate < -1000))
  expect_match(conditionMessage(e), paste0(
    "The error variances of [a-z0-9_]+ and [a-z0-9_]+ run off below zero: ",
    "lavaan's optimizer stopped with them at"))
  # The code printed in the message reproduces the fit in lavaan.
  expect_equal(run_code(e, my_data = holzinger_swineford), e$loglik[["fit"]],
               tolerance = 1e-10)
  # bounds = "pos.var", passed on to lavaan, has a maximum on its boundary.
  b <- quiet(fa_bifactor(holzinger_swineford, f4x, missing = "listwise",
                         bounds = "pos.var", output = "measurement"))
  expect_s3_class(b, "dmar_fa_bifactor")
})

test_that("a finite negative error variance is not a ridge, and is reported", {
  skip_on_cran()
  # The nine tests of fa_bifactor()'s help: visual perception at -41, a
  # maximum of the likelihood, which falls on either side of it.
  fit <- quiet(fa_bifactor(holzinger_swineford, f4x[1:3], output = "fit"))
  pe <- lavaan::parameterEstimates(fit)
  psi <- pe$est[pe$op == "~~" & pe$lhs == "t1_visual_perception" &
                  pe$rhs == pe$lhs]
  expect_lt(psi, -40)
  expect_null(.lavaan_no_maximum(fit))
  # The table is reported, and the warning about the improper solution
  # names the possible causes of a negative error variance at a maximum,
  # not a weak group factor alone.
  r <- with_warnings(suppressMessages(fa_bifactor(holzinger_swineford,
                                                  f4x[1:3])))
  expect_s3_class(r$value, "dmar_fa_bifactor")
  expect_true(any(grepl(paste(
    "A negative error variance at a maximum of the likelihood can come",
    "from sampling error around a small error variance, from a model that",
    "does not describe these items, or from a group factor the data barely",
    "separate from the general factor."), r$warnings, fixed = TRUE)))
  expect_false(any(grepl("often gives an improper solution when the data",
                         r$warnings, fixed = TRUE)))
  # The saturated triad.
  x <- quiet(cfa_1(triad_fin))
  expect_s3_class(x, "dmar_fa_confirmatory")
  expect_equal(x$estimate[x$term == "psi_f1_1"],
               (1 - .8 * .8 / .5) * 199 / 200, tolerance = 1e-4)
  # Two-stage estimation, whose second stage has a likelihood: the check
  # is made and finds the maximum.
  set.seed(7)
  miss <- triad_fin
  miss[cbind(sample(200, 25), sample(3, 25, TRUE))] <- NA
  ts <- quiet(cfa_1(miss, missing = "two_stage", output = "fit"))
  expect_true(is.finite(.no_max_loglik(ts)))
  pe_ts <- lavaan::parameterEstimates(ts)
  expect_lt(pe_ts$est[pe_ts$op == "~~" & pe_ts$lhs == "y1" &
                        pe_ts$rhs == "y1"], 0)
  expect_null(.lavaan_no_maximum(ts))
  # A multiple group fit with the factor scaled by a marker item, y2, and
  # the finite negative variance in one group, which the refit's start
  # moves with the marker loading held at one: the ladder is reported.
  fin_groups <- rbind(cbind(triad_fin, grp = "A"), cbind(other, grp = "B"))
  mi <- quiet(measurement_invariance(fin_groups, items = c("y2", "y1", "y3"),
                                     group = "grp"))
  expect_s3_class(mi, "dmar_measurement_invariance")
  # A proper solution and a least squares fit are not checked at all.
  expect_null(.lavaan_no_maximum(quiet(cfa_1(triad_fin[1:2],
                                             equal_loading = TRUE,
                                             output = "fit"))))
  uls <- quiet(lavaan::cfa("f =~ y1 + y2 + y3", data = triad_fin,
                           estimator = "ULS"))
  expect_null(.lavaan_no_maximum(uls))
})

test_that("a one factor fit with no maximum stops under the name of the function called", {
  skip_on_cran()
  e <- caught(cfa_1(triad_neg, missing = "listwise"))
  expect_s3_class(e, "dmar_no_maximum")
  expect_identical(e$fn, "cfa_1")
  expect_identical(e$variables, "y1")
  expect_lt(e$estimate, -1000)
  expect_match(conditionMessage(e), "^cfa_1\\(\\) stopped")
  # Robust maximum likelihood has the likelihood of maximum likelihood.
  m <- caught(cfa_1(triad_neg, missing = "listwise", estimator = "MLR"))
  expect_s3_class(m, "dmar_no_maximum")
  expect_match(m$lavaan_call, "estimator = \"MLR\"", fixed = TRUE)
  # A covariance matrix with a variable more, and its means in another
  # order: lavaan reads the covariance matrix by name but pairs the means
  # with its rows by position, so the code selects the items from both,
  # and it reproduces the fit. output = "fit" stops too.
  big <- cbind(w = seq_len(300) %% 7 - 3, triad_neg)
  S_big <- stats::cov(big)
  M_big <- rev(colMeans(big))
  s <- caught(fa_confirmatory(S = S_big, M = M_big, N = 300,
                              factors = list(f = paste0("y", 1:3)),
                              output = "fit"))
  expect_s3_class(s, "dmar_no_maximum")
  expect_match(conditionMessage(s), "^fa_confirmatory\\(\\) stopped")
  expect_true(startsWith(s$lavaan_call, "lavaan::lavaan(model,"))
  expect_match(s$lavaan_call, paste0(
    "sample.cov = my_S[c(\"y1\", \"y2\", \"y3\"), ",
    "c(\"y1\", \"y2\", \"y3\")],"), fixed = TRUE)
  expect_match(s$lavaan_call, "sample.nobs = 300,", fixed = TRUE)
  expect_match(s$lavaan_call, "sample.mean = my_M[c(\"y1\", \"y2\", \"y3\")]",
               fixed = TRUE)
  expect_match(conditionMessage(s), paste(
    "your covariance matrix in place of my_S and its means in place of",
    "my_M"), fixed = TRUE)
  expect_equal(run_code(s, my_S = S_big, my_M = M_big), s$loglik[["fit"]],
               tolerance = 1e-8)
})

test_that("a rung of an invariance ladder with no maximum stops, naming the rung and the group", {
  skip_on_cran()
  e <- caught(measurement_invariance(two_groups, items = paste0("y", 1:3),
                                     group = "grp", std.lv = TRUE))
  expect_s3_class(e, "dmar_no_maximum")
  expect_identical(e$fn, "measurement_invariance")
  expect_identical(e$variables, "y1")
  expect_identical(e$group, "A")
  expect_identical(e$where, "the configural rung")
  expect_match(conditionMessage(e), paste0(
    "no maximum \\(the configural rung\\)\\. The error variance of y1 in ",
    "group A runs off below zero"))
  expect_match(e$lavaan_call, "^lavaan::cfa\\(model, data = my_data")
  expect_match(e$lavaan_call, "group = \"grp\"", fixed = TRUE)
  expect_match(e$lavaan_call, "std.lv = TRUE", fixed = TRUE)
  expect_equal(run_code(e, my_data = two_groups), e$loglik[["fit"]],
               tolerance = 1e-8)
})

test_that("a rung of an invariance ladder that does not converge stops, naming the rung", {
  skip_on_cran()
  # Under the default scaling by a marker item, lavaan's optimizer does not
  # converge on the triad's group; through DMAR 1.0.x lavaan's error from
  # its fit measures ended the call without naming the rung.
  e <- tryCatch(quiet(measurement_invariance(two_groups,
                                             items = paste0("y", 1:3),
                                             group = "grp")),
                dmar_not_converged = function(e) e)
  expect_s3_class(e, c("dmar_not_converged", "dmar_refusal", "error", "condition"),
                  exact = TRUE)
  expect_identical(e$level, "configural")
  expect_match(conditionMessage(e), paste0(
    "^measurement_invariance\\(\\) stopped: lavaan's optimizer did not ",
    "converge on the configural rung"))
  expect_identical(e$syntax, "f =~ y1 + y2 + y3")
  expect_match(e$lavaan_call, "^lavaan::cfa\\(model, data = my_data")
  expect_match(e$lavaan_call, "group = \"grp\"", fixed = TRUE)
  # The code reproduces the fit lavaan did not converge on.
  env <- list2env(list(model = e$syntax, my_data = two_groups))
  f <- quiet(eval(parse(text = e$lavaan_call), envir = env))
  expect_false(lavaan::lavInspect(f, "converged"))
})

test_that("measurement_alignment() refuses a configural fit with no maximum", {
  skip_on_cran()
  # DMAR 1.0.0 aligned the point where lavaan stopped (y1 in group A with a
  # loading of 21.8 and an error variance of -475) and gave group B a
  # factor variance of 1282.
  e <- caught(measurement_alignment(two_groups, items = paste0("y", 1:3),
                                    group = "grp", seed = 113))
  expect_s3_class(e, c("dmar_no_maximum", "dmar_refusal", "error",
                       "condition"), exact = TRUE)
  expect_identical(e$fn, "measurement_alignment")
  expect_identical(e$variables, "y1")
  expect_identical(e$group, "A")
  expect_identical(e$where, "the configural model")
  expect_lt(e$estimate, -100)
  expect_match(conditionMessage(e), paste0(
    "^measurement_alignment\\(\\) stopped: the fit failed because the ",
    "likelihood has no maximum \\(the configural model\\)"))
  expect_identical(e$syntax, "f =~ y1 + y2 + y3")
  expect_match(e$lavaan_call, paste0(
    "^lavaan::cfa\\(model, data = my_data, group = \"grp\", ",
    "estimator = \"ML\","))
  expect_match(e$lavaan_call, "std.lv = TRUE, meanstructure = TRUE)",
               fixed = TRUE)
})

test_that("an exploratory fit with no maximum is refused, rotated or not, with code that selects the items", {
  skip_on_cran()
  items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
             "t4_lozenges", "t5_general_information",
             "t6_paragraph_comprehension", "t7_sentence",
             "t8_word_classification", "t9_word_meaning", "t10_addition",
             "t11_code", "t12_counting_groups_of_dots",
             "t13_straight_and_curved_capitals")
  # The whole data set, whose school and sex are factors, with the items
  # named, and the default rotation: the check reads the unrotated fit.
  e <- caught(fa_exploratory(holzinger_swineford, items = items,
                             nfactors = 7))
  expect_s3_class(e, c("dmar_no_maximum", "dmar_efa_no_solution",
                       "dmar_refusal", "error", "condition"), exact = TRUE)
  expect_identical(e$reason, "no maximum")
  expect_identical(e$fn, "fa_exploratory")
  expect_null(e$syntax)
  expect_match(e$lavaan_call, paste0(
    "^lavaan::efa\\(data = my_data\\[c\\(\"t1_visual_perception\", ",
    "\"t2_cubes\","))
  expect_match(e$lavaan_call, "missing = \"ml\"", fixed = TRUE)
  expect_match(e$lavaan_call, "nfactors = 7L", fixed = TRUE)
  expect_match(e$lavaan_call, "rotation = \"none\"", fixed = TRUE)
  expect_match(e$lavaan_call, "bounds = \"none\"", fixed = TRUE)
  expect_false(grepl("  \"", e$lavaan_call))
  # Run on the whole data set, the code fits the items alone, as
  # fa_exploratory() did; lavaan::efa() given the whole data set would fit
  # every column.
  expect_equal(run_code(e, my_data = holzinger_swineford), e$loglik[["fit"]],
               tolerance = 1e-8)
  # Covariance input: lavaan::efa() reads the variables from the row names
  # of the matrix it is given, so the code selects the items' rows and
  # columns when the matrix holds others.
  s <- tryCatch(.efa_no_maximum(
    list(variable = "y1", group = NULL, estimate = -50, fixed_at = -500,
         loglik = c(fit = -100, further = -99.9), tol = 1e-4, N = 100),
    list(sample_cov = diag(3), sample_nobs = 100, nfactors = 1L),
    select = paste0("y", 1:3)), dmar_no_maximum = function(e) e)
  expect_true(startsWith(s$lavaan_call, "lavaan::efa(sample_cov = my_S["))
  expect_match(s$lavaan_call, paste0(
    "sample_cov = my_S[c(\"y1\", \"y2\", \"y3\"), ",
    "c(\"y1\", \"y2\", \"y3\")],"), fixed = TRUE)
  expect_s3_class(s, "dmar_efa_no_solution")
})

test_that("factor_count() reports a count with no maximum as a row", {
  skip_on_cran()
  orig <- fa_exploratory
  local_mocked_bindings(fa_exploratory = function(...) {
    args <- list(...)
    if (identical(as.integer(args$nfactors), 2L)) {
      .efa_no_maximum(list(variable = "t2_cubes", group = NULL,
                           estimate = -50, fixed_at = -500,
                           loglik = c(fit = -100, further = -99.9)),
                      list(data = args$data, nfactors = 2L))
    }
    orig(...)
  })
  items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
             "t4_lozenges", "t5_general_information",
             "t6_paragraph_comprehension")
  x <- quiet(factor_count(holzinger_swineford[items], max_factors = 2))
  expect_identical(x$solution, c("proper", "no maximum"))
  expect_true(is.na(x$chi_square[2L]))
  expect_identical(unname(attr(x, "reasons")[2L]), paste(
    "the likelihood has no maximum: the error variance of t2_cubes runs",
    "off below zero"))
  expect_true(any(grepl(paste("no maximum: a residual variance runs off",
                              "below zero while the likelihood is flat or",
                              "still rising"),
                        .dmar_tbl_footer(x), fixed = TRUE)))
})

test_that("a second-order fit with no maximum stops, and at the boundary the code is the one checked", {
  skip_on_cran()
  # Four first-order factors under one, the fourth on the triad: lavaan
  # stops on the ridge of y1, the second-order loading of that factor
  # near zero. Checked in fa_second_order()'s own parameterization, where
  # the factors regress on the second-order factor.
  set.seed(11)
  n <- 300
  g <- rnorm(n)
  fs <- lapply(c(.8, .7, .6), function(l) l * g + sqrt(1 - l^2) * rnorm(n))
  block <- function(f, pre) {
    x <- sapply(1:3, function(i) .7 * f + sqrt(.51) * rnorm(n))
    colnames(x) <- paste0(pre, 1:3)
    x
  }
  blocks <- Map(block, fs, c("x", "z", "w"))
  tri <- matrix(rnorm(n * 3), n) %*% chol(S_neg)
  colnames(tri) <- paste0("y", 1:3)
  d <- as.data.frame(do.call(cbind, c(blocks, list(tri))))
  fac <- list(fx = paste0("x", 1:3), fz = paste0("z", 1:3),
              fw = paste0("w", 1:3), fy = paste0("y", 1:3))
  e <- caught(fa_second_order(S = stats::cov(d), N = n, factors = fac,
                              second_order = list(g = names(fac))))
  expect_s3_class(e, "dmar_no_maximum")
  expect_identical(e$fn, "fa_second_order")
  expect_identical(e$variables, "y1")
  expect_lt(e$estimate, -1000)
  expect_match(e$syntax, "g =~", fixed = TRUE)
  expect_false(grepl("lower(0)", e$syntax, fixed = TRUE))
  expect_match(conditionMessage(e), "This lavaan code reproduces the fit",
               fixed = TRUE)
  # A standardized second-order loading at one: the check reads the fit
  # in the parameterization that reaches the boundary, and the code is
  # that fit's syntax and arguments, for output = "fit" as well. The
  # stand-in for the check finds a ridge in that fit alone, the one with
  # a bound.
  L <- matrix(0, 9L, 3L)
  for (j in 1:3) L[3L * j - 2:0, j] <- 0.7
  Phi <- matrix(c(1, .9, .9, .9, 1, .6, .9, .6, 1), 3L)
  R <- L %*% Phi %*% t(L)
  diag(R) <- 1
  dimnames(R) <- list(paste0("x", 1:9), paste0("x", 1:9))
  f <- list(a = paste0("x", 1:3), b = paste0("x", 4:6), c = paste0("x", 7:9))
  ridge <- fake_ridge("x2")
  local_mocked_bindings(.lavaan_no_maximum = function(fit, step = 10) {
    if (any(is.finite(lavaan::parTable(fit)$lower))) ridge(fit, step)
  })
  for (out in c("verbose", "fit")) {
    b <- caught(fa_second_order(S = R, N = 300, factors = f,
                                second_order = list(g = names(f)),
                                output = out))
    expect_s3_class(b, "dmar_no_maximum")
    expect_identical(b$fn, "fa_second_order")
    expect_match(b$syntax, "lower(0)*a", fixed = TRUE)
    expect_match(b$lavaan_call, "ceq.simple = TRUE", fixed = TRUE)
    expect_match(b$lavaan_call,
                 "^lavaan::lavaan\\(model, sample.cov = my_S, sample.nobs = 300,")
    expect_match(conditionMessage(b),
                 "This lavaan code fits the model as DMAR checked it",
                 fixed = TRUE)
  }
})

test_that("the refit's start keeps every variable's implied variance, with a marker held and with factors on factors", {
  skip_on_cran()
  # A marker item (y2) scales the factor in each group; y1's error
  # variance in group A is held at ten times its estimate.
  fin_groups <- rbind(cbind(triad_fin, grp = "A"), cbind(other, grp = "B"))
  fit <- quiet(lavaan::cfa("f =~ y2 + y1 + y3", data = fin_groups,
                           group = "grp", missing = "ml"))
  pt <- lavaan::parTable(fit)
  r <- which(pt$op == "~~" & pt$lhs == "y1" & pt$rhs == "y1" &
               pt$group == 1L)
  st <- .no_max_start(pt, r, 10 * pt$est[r], .no_max_matrices(fit))
  expect_equal(st[r], 10 * pt$est[r])
  marker <- pt$op == "=~" & pt$rhs == "y2"
  expect_equal(st[marker], c(1, 1))
  for (b in 1:2) {
    expect_equal(diag(implied_at(pt, st, b)), diag(implied_at(pt, pt$est, b)),
                 tolerance = 1e-10)
  }
  # y1's loading grows, and the factor's variance in group A shrinks to
  # keep the marker's loading at one.
  expect_gt(st[pt$op == "=~" & pt$rhs == "y1" & pt$group == 1L],
            pt$est[pt$op == "=~" & pt$rhs == "y1" & pt$group == 1L])
  # A second-order model, whose first-order factors regress on the second
  # order (lavaan's beta): t2_cubes's error variance held one below its
  # estimate.
  hs_four <- list(
    spatial = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
    verbal  = c("t6_paragraph_comprehension", "t7_sentence",
                "t9_word_meaning"),
    speed   = c("t10_addition", "t12_counting_groups_of_dots",
                "t13_straight_and_curved_capitals"),
    memory  = c("t14_word_recognition", "t15_number_recognition",
                "t17_object_number"))
  so <- quiet(fa_second_order(holzinger_swineford, hs_four,
                              second_order = list(g = names(hs_four)),
                              output = "fit"))
  mats <- .no_max_matrices(so)
  expect_false(is.null(mats[[1L]]$beta))
  pt <- lavaan::parTable(so)
  r <- which(pt$op == "~~" & pt$lhs == "t2_cubes" & pt$rhs == "t2_cubes")
  st <- .no_max_start(pt, r, pt$est[r] - 1, mats)
  expect_equal(st[r], pt$est[r] - 1)
  expect_equal(diag(implied_at(pt, st)), diag(implied_at(pt, pt$est)),
               tolerance = 1e-8)
})

test_that("the reliability family refuses under its own name, with the items' names", {
  local_mocked_bindings(.lavaan_no_maximum = fake_ridge("y2"))
  d <- holzinger_swineford[c("t1_visual_perception", "t2_cubes",
                             "t3_paper_form_board", "t4_lozenges")]
  e <- caught(reliability_omega(d, denominator = "model_implied",
                                ci_method = "ml"))
  expect_s3_class(e, "dmar_no_maximum")
  expect_identical(e$fn, "reliability_omega")
  expect_identical(e$variables, "t2_cubes")
  expect_match(conditionMessage(e), paste0(
    "^reliability_omega\\(\\) stopped: the fit failed because the ",
    "likelihood has no maximum\\. The error variance of t2_cubes runs off"))
  expect_match(e$syntax, "lambda_f1_2*t2_cubes", fixed = TRUE)
  expect_false(grepl("\\by[1-4]\\b", e$syntax, perl = TRUE))
  expect_equal(e$N, 100)
  expect_match(conditionMessage(caught(reliability_omega(d))),
               "^reliability_omega\\(\\) stopped")
  expect_match(conditionMessage(caught(reliability_alpha(
    d, estimator = "model_implied", ci_method = "ml"))),
    "^reliability_alpha\\(\\) stopped")
  expect_match(conditionMessage(caught(reliability_H(data = d))),
               "^reliability_H\\(\\) stopped")
  # FIML with an auxiliary variable fits through lavaan directly, and is
  # checked as well.
  a <- holzinger_swineford[c("t1_visual_perception", "t2_cubes",
                             "t3_paper_form_board", "t4_lozenges",
                             "t5_general_information")]
  a$t2_cubes[1:20] <- NA
  ea <- caught(reliability_omega(a, aux = "t5_general_information",
                                 missing = "fiml",
                                 denominator = "model_implied",
                                 ci_method = "ml"))
  expect_s3_class(ea, "dmar_no_maximum")
  expect_identical(ea$variables, "t2_cubes")
  expect_match(ea$syntax, "t5_general_information", fixed = TRUE)
})

test_that("a bootstrap replication with no maximum is counted among the errors", {
  d <- holzinger_swineford[c("t1_visual_perception", "t2_cubes",
                             "t3_paper_form_board", "t4_lozenges")]
  # The fit to every row has a maximum; a resample without the first case
  # has none.
  local_mocked_bindings(.lavaan_no_maximum = fake_ridge(
    "y2", unless = unlist(d[1L, ])))
  r <- quiet(reliability_omega(d, ci_method = "percentile", B = 20,
                               seed = 113))
  bs <- attr(r, "bootstrap")
  n_err <- sum(bs$status == "error")
  expect_gt(n_err, 0L)
  expect_lt(n_err, 20L)
  expect_identical(bs$reasons[["error"]], n_err)
  expect_identical(bs$attempted, 20L)
  expect_identical(bs$returned + bs$failed, 20L)
  expect_identical(names(bs$error_messages), paste(
    "the likelihood has no maximum: the error variance of t2_cubes runs off",
    "below zero"))
  expect_identical(unname(bs$error_messages), n_err)
  expect_match(.bootstrap_footer_line(r), paste0(
    "raised an error \\(most often: the likelihood has no maximum"))
})

test_that("the replications with no maximum are refit on the engine's processes, with the counts of one process", {
  skip_on_os("windows")   # forked workers, which see the stand-in check
  skip_on_cran()
  # This guards that the results are the same on one process and on two;
  # that the refits run in the workers is tested in
  # test-reliability_bootstrap.R.
  d <- holzinger_swineford[c("t1_visual_perception", "t2_cubes",
                             "t3_paper_form_board", "t4_lozenges")]
  local_mocked_bindings(.lavaan_no_maximum = fake_ridge(
    "y2", unless = unlist(d[1L, ])))
  r1 <- quiet(reliability_omega(d, ci_method = "percentile", B = 20,
                                seed = 113))
  r2 <- quiet(reliability_omega(d, ci_method = "percentile", B = 20,
                                seed = 113, ncpus = 2))
  b1 <- attr(r1, "bootstrap")
  b2 <- attr(r2, "bootstrap")
  expect_gt(b1$reasons[["error"]], 0L)
  expect_identical(b2$ncpus, 2L)
  b2$ncpus <- 1L
  expect_identical(b2, b1)
  expect_identical(r2$value, r1$value)
})

test_that("ci_bifactor_indices() and reliability_H() on a fit count a replication with no maximum", {
  skip_on_cran()
  hs_groups <- list(
    spatial = c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
                "t4_lozenges"),
    verbal  = c("t5_general_information", "t6_paragraph_comprehension",
                "t7_sentence", "t8_word_classification", "t9_word_meaning"),
    memory  = c("t14_word_recognition", "t15_number_recognition",
                "t16_figure_recognition", "t17_object_number",
                "t18_number_figure", "t19_figure_word"))
  fit <- quiet(fa_bifactor(holzinger_swineford, hs_groups, output = "fit"))
  one <- quiet(cfa_1(holzinger_swineford[hs_groups$spatial],
                     output = "fit"))
  pool <- lavaan::lavInspect(fit, "data")
  local_mocked_bindings(.lavaan_no_maximum = fake_ridge(
    "t2_cubes", unless = pool[1L, ]))
  ci <- quiet(ci_bifactor_indices(fit, B = 10, seed = 113))
  bs <- attr(ci, "bootstrap")
  n_err <- sum(bs$status == "error")
  expect_gt(n_err, 0L)
  expect_identical(bs$reasons[["error"]], n_err)
  expect_identical(names(bs$error_messages), paste(
    "the likelihood has no maximum: the error variance of t2_cubes runs off",
    "below zero"))
  # The refit of the data stored with the fit is held to the same rule:
  # with no row exempt, the call stops.
  local_mocked_bindings(.lavaan_no_maximum = fake_ridge("t2_cubes"))
  e <- caught(ci_bifactor_indices(fit, B = 10, seed = 113))
  expect_s3_class(e, "dmar_no_maximum")
  expect_identical(e$fn, "ci_bifactor_indices")
  expect_null(e$lavaan_call)
  # reliability_H() on a single factor fit.
  local_mocked_bindings(.lavaan_no_maximum = fake_ridge(
    "t2_cubes", unless = lavaan::lavInspect(one, "data")[1L, ]))
  h <- quiet(reliability_H(fit = one, ci_method = "percentile", B = 20,
                           seed = 113))
  hb <- attr(h, "bootstrap")
  expect_gt(sum(hb$status == "error"), 0L)
  expect_match(names(hb$error_messages)[1L], "the likelihood has no maximum")
})

test_that("the refusal helpers read the constraints, build the code, and word the message to the check", {
  pt <- lavaan::parTable(quiet(lavaan::lavaan(paste(
    "f =~ NA*y1 + l1*y1 + l2*y2 + l3*y3", "f ~~ 1*f",
    "y1 ~~ p*y1", "y2 ~~ p*y2", "y3 ~~ q*y3", sep = "\n"),
    data = triad_neg)))
  r1 <- which(pt$op == "~~" & pt$lhs == "y1" & pt$rhs == "y1")
  r2 <- which(pt$op == "~~" & pt$lhs == "y2" & pt$rhs == "y2")
  # Rows that share a label are held together.
  expect_identical(.no_max_tied_rows(pt, r1), sort(c(r1, r2)))
  # A constraint other than a simple equality leaves the check undefined.
  pt2 <- lavaan::parTable(quiet(lavaan::lavaan(paste(
    "f =~ NA*y1 + l1*y1 + l2*y2 + l3*y3", "f ~~ 1*f",
    "y1 ~~ p*y1", "y2 ~~ q*y2", "y3 ~~ s*y3", "p == 2*q", sep = "\n"),
    data = triad_neg)))
  expect_null(.no_max_tied_rows(pt2, which(pt2$op == "~~" &
                                              pt2$lhs == "y1" &
                                              pt2$rhs == "y1")))
  # The call text: the data a placeholder, broken between arguments only,
  # and a selection where lavaan would read every variable it is given.
  txt <- .no_max_call_text("lavaan::lavaan", list(
    model = "f =~ y1", estimator = "ML", se = "standard",
    meanstructure = TRUE, data = triad_neg, missing = "fiml",
    control = list(rel.tol = 1e-10)))
  expect_match(txt, "^lavaan::lavaan\\(model, data = my_data, estimator")
  expect_true(all(nchar(strsplit(txt, "\n")[[1L]]) <= 72L))
  expect_false(grepl("= *\n", txt))
  expect_match(txt, "control = list(rel.tol = 1e-10))", fixed = TRUE)
  sel <- .no_max_call_text("lavaan::efa", list(data = triad_neg,
                                               nfactors = 1L),
                           select = list(data = c("y1", "y3")))
  expect_identical(sel, paste0("lavaan::efa(data = my_data[c(\"y1\", ",
                               "\"y3\")], nfactors = 1L)"))
  # The message: a refit further out that is higher says the likelihood
  # is still rising; one lower by no more than the tolerance says it is
  # flat, with the tolerance and the cases.
  nm <- list(variable = "y1", group = NULL, estimate = -50, fixed_at = -500,
             loglik = c(fit = -100, further = -100.00002), tol = 3e-4,
             N = 300)
  flat <- .no_max_message("cfa_1", nm, syntax = NULL, lavaan_call = NULL)
  expect_match(flat, paste(
    "the log likelihood is lower by only 0.00002, within the tolerance of",
    "1e-6 per case (0.0003 for these 300 cases): the likelihood is flat as",
    "the variance falls."), fixed = TRUE)
  expect_true(endsWith(flat, "or a different model can be fit."))
  # The number of cases is written with commas, as a population fit's
  # million is.
  big <- .no_max_message("ss_aipe_sem_path",
                         utils::modifyList(nm, list(tol = 1, N = 1e6)),
                         syntax = NULL, lavaan_call = NULL)
  expect_match(big, "(1 for these 1,000,000 cases)", fixed = TRUE)
  # The ":=" lines of DMAR's syntax define its indices and are left out,
  # with a note; those of a model the user wrote are the user's own
  # quantities, and the model is printed whole.
  syn <- "f =~ y1 + y2 + y3\nomega := 1"
  call <- "lavaan::sem(model, data = my_data)"
  dropped <- .no_max_message("cfa_1", nm, syntax = syn, lavaan_call = call)
  expect_match(dropped, "model <- '\nf =~ y1 + y2 + y3\n'\n", fixed = TRUE)
  expect_match(dropped, paste("(the lines of DMAR's syntax that define its",
                              "indices are left out; they do not change the",
                              "fit)"), fixed = TRUE)
  kept <- .no_max_message("ss_aipe_sem_path", nm, syntax = syn,
                          lavaan_call = call, defined = "keep")
  expect_match(kept, paste0("model <- '\n", syn, "\n'\n"), fixed = TRUE)
  expect_false(grepl("left out", kept, fixed = TRUE))
  # A model whose syntax would carry the message past what R prints of an
  # error is printed as its call, the syntax left in the condition.
  long <- paste(sprintf("y%d ~~ psi_%d*y%d", 1:600, 1:600, 1:600),
                collapse = "\n")
  short <- .no_max_message("fa_confirmatory", nm, syntax = long,
                           lavaan_call = "lavaan::lavaan(model, data = my_data)")
  expect_lte(nchar(short, type = "bytes"), 8000L)
  expect_match(short, "The model syntax DMAR built, 600 lines, is too long",
               fixed = TRUE)
  expect_true(endsWith(short, "fit <- lavaan::lavaan(model, data = my_data)"))
  expect_false(grepl("model <- '", short, fixed = TRUE))
  users <- .no_max_message("ss_aipe_sem_path", nm, syntax = long,
                           lavaan_call = "lavaan::sem(model)",
                           defined = "keep")
  expect_match(users, "The model syntax, 600 lines, is too long", fixed = TRUE)
  # The condition carries the fields, however the message is cut, and the
  # limit on what R prints of a message is restored once the refusal has
  # been caught.
  before <- getOption("warning.length")
  e <- tryCatch(.no_max_refuse("fa_confirmatory", nm, syntax = long,
                               lavaan_call = "lavaan::lavaan(model)"),
                dmar_no_maximum = function(e) e)
  expect_identical(e$syntax, long)
  expect_equal(e$tolerance, 3e-4)
  expect_identical(.no_max_fields(e)$tol, 3e-4)
  expect_identical(getOption("warning.length"), before)
})

# ---- The structural equation models -----------------------------------------

# The covariance matrix of the triad above as a population: the single
# factor model with free loadings has no maximum there. lavaan stops with
# the error variance of y1 near -1900 at a million cases and reports
# convergence, after thousands of iterations (4 s); with y1 as a marker it
# reports no convergence, and the check, which reads only a converged fit,
# does not apply.
S_tri <- stats::cov(triad_neg)
one_factor <- "f =~ NA*y1 + b*y2 + y3\nf ~~ 1*f"

# The planning model of ss_aipe_sem_path()'s help, a proper population for
# the tests of what a planner does with the check's answer.
pop_2f <- paste(
  "f1 =~ 1*y1 + 0.8*y2 + 0.8*y3", "f2 =~ 1*y4 + 0.8*y5 + 0.8*y6",
  "f2 ~ 0.5*f1", "f1 ~~ 1*f1", "f2 ~~ 0.75*f2",
  "y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3",
  "y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6", sep = "\n")
Sigma_2f <- cov_sem(pop_2f)$sigma_theta
two_factor <- "f1 =~ y1 + y2 + y3\nf2 =~ y4 + y5 + y6\nf2 ~ b*f1"

# A latent mediator between x and y, and the saturated triad as its items:
# m1's error variance is negative, at a maximum of the likelihood.
sem_model <- "M =~ m1 + m2 + m3\nM ~ a*x\ny ~ b*M + x"
heywood_sem <- function() {
  d <- stats::setNames(triad_fin, paste0("m", 1:3))
  set.seed(5)
  d$x <- stats::rnorm(200)
  d$y <- 0.3 * d$x + stats::rnorm(200)
  d
}

# A stand-in for the check that finds a ridge in every fit for which
# `flag(fit)` is TRUE and records the log likelihood of that fit;
# `seen$calls` counts the fits it is given.
recording_ridge <- function(flag = function(fit) TRUE, variable = "y1") {
  seen <- new.env(parent = emptyenv())
  seen$calls <- 0L
  list(seen = seen, check = function(fit, step = 10) {
    seen$calls <- seen$calls + 1L
    if (!flag(fit)) return(NULL)
    seen$ll <- as.numeric(lavaan::logLik(fit))
    fake_ridge(variable)(fit, step)
  })
}
# The planners fit the population moments at 1e5 or 1e6 cases.
at_population <- function(fit) lavaan::lavInspect(fit, "ntotal") >= 1e5

# A stand-in for the check that finds a ridge in every other fit to a
# simulated sample, and counts those fits.
every_other_sample <- function() {
  state <- new.env(parent = emptyenv())
  state$samples <- 0L
  list(state = state, check = function(fit, step = 10) {
    if (at_population(fit)) return(NULL)
    state$samples <- state$samples + 1L
    if (state$samples %% 2L == 1L) fake_ridge("y1")(fit, step)
  })
}

# The model implied covariance matrix of the observed variables of a single
# group structural model's parameter table `pt` at the values `st`, in the
# reticular action form, which needs no latent variable for an observed
# one: every loading and regression a directed path, every variance and
# covariance a two-headed one, (I - A)^-1 S (I - A)^-T. lavaan's own form,
# in which an observed variable that enters a regression is a latent
# variable of its own, implies the same matrix.
implied_sem <- function(pt, st) {
  rows <- which(pt$op %in% c("=~", "~", "~~"))
  vars <- unique(c(pt$lhs[rows], pt$rhs[rows]))
  ov <- setdiff(vars, unique(pt$lhs[pt$op == "=~"]))
  A <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  S <- A
  for (r in rows) {
    l <- pt$lhs[r]
    h <- pt$rhs[r]
    if (pt$op[r] == "=~") {
      A[h, l] <- st[r]
    } else if (pt$op[r] == "~") {
      A[l, h] <- st[r]
    } else {
      S[l, h] <- S[h, l] <- st[r]
    }
  }
  IA <- solve(diag(length(vars)) - A)
  (IA %*% S %*% t(IA))[ov, ov]
}

# The three factor model of ss_power_composite_sem()'s help: a population
# and an analysis model whose paths a and b carry an indirect effect.
pop_3f <- paste(
  "f1 =~ 1*y1 + 0.8*y2 + 0.8*y3", "f2 =~ 1*y4 + 0.8*y5 + 0.8*y6",
  "f3 =~ 1*y7 + 0.8*y8 + 0.8*y9", "f2 ~ 0.4*f1", "f3 ~ 0.3*f2 + 0.25*f1",
  "f1 ~~ 1*f1", "f2 ~~ 0.84*f2", "f3 ~~ 0.8*f3",
  "y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3",
  "y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6",
  "y7 ~~ 0.5*y7; y8 ~~ 0.5*y8; y9 ~~ 0.5*y9", sep = "\n")
three_factor <- paste("f1 =~ y1 + y2 + y3", "f2 =~ y4 + y5 + y6",
                      "f3 =~ y7 + y8 + y9", "f2 ~ a*f1", "f3 ~ b*f2 + c*f1",
                      sep = "\n")

# Four weak indicators of one factor: the fit has a maximum, and resamples
# of its transformed data can fall on a ridge (seed 2 at B = 60 under seed
# 113 gives 7), which lavaan reports as converged and its post-estimation
# check flags.
weak <- function(seed) {
  set.seed(seed)
  f <- stats::rnorm(60)
  d <- as.data.frame(sapply(1:4, function(j) 0.35 * f + stats::rnorm(60)))
  names(d) <- paste0("v", 1:4)
  quiet(lavaan::cfa("f =~ v1 + v2 + v3 + v4", data = d, std.lv = TRUE))
}

test_that("ss_aipe_sem_path() refuses an analysis model with no maximum at true_cov_matrix, with code that reproduces the fit", {
  skip_on_cran()
  e <- caught(ss_aipe_sem_path(one_factor, true_cov_matrix = S_tri, width = 0.3,
                               which_path = "b"))
  expect_s3_class(e, c("dmar_no_maximum", "dmar_refusal", "error",
                       "condition"), exact = TRUE)
  expect_identical(e$fn, "ss_aipe_sem_path")
  expect_identical(e$where, "the analysis model fit to true_cov_matrix")
  expect_identical(e$variables, "y1")
  expect_lt(e$estimate, -1000)
  expect_gt(e$loglik[["further"]], e$loglik[["fit"]])
  expect_equal(e$N, 1e6)
  msg <- conditionMessage(e)
  expect_match(msg, paste0(
    "^ss_aipe_sem_path\\(\\) stopped: the fit failed because the ",
    "likelihood has no maximum \\(the analysis model fit to true_cov_matrix\\)\\. ",
    "The error variance of y1 runs off below zero"))
  expect_match(msg, paste("There, bounds = \"pos.var\" keeps every error",
                          "variance at zero or above"), fixed = TRUE)
  expect_match(msg, paste("This lavaan code reproduces the fit, with the",
                          "matrix given as true_cov_matrix in place of my_S;"),
               fixed = TRUE)
  expect_identical(e$syntax, one_factor)
  expect_identical(e$lavaan_call,
                   "lavaan::sem(model, sample.cov = my_S, sample.nobs = 1e+06)")
  expect_equal(run_code(e, my_S = S_tri), e$loglik[["fit"]], tolerance = 1e-8)
})

test_that("the other SEM planners refuse their fit to the population moments under their own names, before any replication", {
  skip_on_cran()
  # The arguments passed on to lavaan are in the code, which reproduces the
  # fit the planner made.
  ridge <- recording_ridge(at_population)
  local_mocked_bindings(.lavaan_no_maximum = ridge$check)
  p <- caught(ss_power_composite_sem(two_factor, true_cov_matrix = Sigma_2f, specified_N = 100,
                                     G = 10, std.lv = TRUE))
  expect_s3_class(p, "dmar_no_maximum")
  expect_identical(p$fn, "ss_power_composite_sem")
  expect_identical(p$where,
                   "the analysis model fit to the population covariance matrix")
  expect_identical(ridge$seen$calls, 1L)
  expect_identical(gsub("\\s+", " ", p$lavaan_call), paste(
    "lavaan::sem(model, sample.cov = my_S, sample.nobs = 1e+06,",
    "std.lv = TRUE)"))
  expect_match(conditionMessage(p),
               "with the matrix given as true_cov_matrix in place of my_S;", fixed = TRUE)
  expect_equal(run_code(p, my_S = Sigma_2f), ridge$seen$ll, tolerance = 1e-8)
  # From true_model, the matrix is the one cov_sem() derives.
  ridge <- recording_ridge(at_population)
  local_mocked_bindings(.lavaan_no_maximum = ridge$check)
  a <- caught(ss_aipe_composite_sem(two_factor, true_model = pop_2f,
                                    width = 0.3, specified_N = 100, G = 10))
  expect_identical(a$fn, "ss_aipe_composite_sem")
  expect_match(conditionMessage(a),
               "with cov_sem(true_model)$sigma_theta in place of my_S;",
               fixed = TRUE)
  expect_equal(run_code(a, my_S = cov_sem(pop_2f)$sigma_theta),
               ridge$seen$ll, tolerance = 1e-8)
  # With a mean structure, means named in another order than the rows of
  # Sigma are selected by name, since lavaan pairs them by position.
  mu <- stats::setNames(c(0.6, 0.5, 0.4, 0.3, 0.2, 0.1), paste0("y", 6:1))
  ridge <- recording_ridge(at_population)
  local_mocked_bindings(.lavaan_no_maximum = ridge$check)
  m <- caught(ss_aipe_composite_sem(paste(two_factor, "y1 ~ m1*1", sep = "\n"),
                                    true_cov_matrix = Sigma_2f, true_mean = mu,
                                    width = 0.3, specified_N = 100, G = 10))
  expect_identical(m$where, paste("the analysis model fit to the population",
                                  "means and covariance matrix"))
  expect_match(gsub("\\s+", " ", m$lavaan_call), paste0(
    "sample.mean = my_M[c(\"y1\", \"y2\", \"y3\", \"y4\", \"y5\", \"y6\")]"),
    fixed = TRUE)
  expect_match(conditionMessage(m), paste(
    "with the matrix given as true_cov_matrix in place of my_S and the means given",
    "as true_mean in place of my_M;"), fixed = TRUE)
  expect_equal(run_code(m, my_S = Sigma_2f, my_M = mu), ridge$seen$ll,
               tolerance = 1e-8)
  # ss_aipe_sem_path_sensitivity() raises the planner's refusal under its
  # own name, and fits true_cov_matrix before any replication: the check is made
  # twice, of the planning fit and of that one.
  ridge <- recording_ridge(at_population)
  local_mocked_bindings(.lavaan_no_maximum = ridge$check)
  s <- caught(ss_aipe_sem_path_sensitivity(two_factor, supposed_cov_matrix = Sigma_2f,
                                           which_path = "b",
                                           width = 0.3, G = 5))
  expect_identical(s$fn, "ss_aipe_sem_path_sensitivity")
  expect_identical(s$where, "the analysis model fit to supposed_cov_matrix")
  expect_match(conditionMessage(s),
               "^ss_aipe_sem_path_sensitivity\\(\\) stopped")
  expect_match(conditionMessage(s),
               "with the matrix given as supposed_cov_matrix in place of my_S;",
               fixed = TRUE)
  ridge <- recording_ridge(function(fit) {
    lavaan::lavInspect(fit, "ntotal") == 1e5
  })
  local_mocked_bindings(.lavaan_no_maximum = ridge$check)
  t <- caught(ss_aipe_sem_path_sensitivity(two_factor, supposed_cov_matrix = Sigma_2f,
                                           which_path = "b",
                                           width = 0.3, specified_N = 100,
                                           G = 5))
  expect_identical(t$where, "the analysis model fit to true_cov_matrix")
  expect_identical(ridge$seen$calls, 2L)
  expect_identical(t$lavaan_call,
                   "lavaan::sem(model, sample.cov = my_S, sample.nobs = 1e+05)")
  expect_equal(run_code(t, my_S = Sigma_2f), ridge$seen$ll, tolerance = 1e-8)
})

test_that("a planner prints the user's model whole, its := lines included", {
  skip_on_cran()
  # The := lines of a model the user wrote define the user's own
  # quantities, an indirect effect here, and are printed with the model;
  # only the lines of DMAR's own syntax that define its indices are left
  # out, and the message says so only then.
  ridge <- recording_ridge(at_population)
  local_mocked_bindings(.lavaan_no_maximum = ridge$check)
  model <- paste(three_factor, "ab := a*b", sep = "\n")
  e <- caught(ss_power_composite_sem(model, true_model = pop_3f, specified_N = 100,
                                     G = 10))
  msg <- conditionMessage(e)
  expect_identical(e$syntax, model)
  expect_match(msg, paste0("model <- '\n", model, "\n'\nfit <- "),
               fixed = TRUE)
  expect_false(grepl("left out", msg, fixed = TRUE))
  expect_equal(run_code(e, my_S = cov_sem(pop_3f)$sigma_theta),
               ridge$seen$ll, tolerance = 1e-8)
  # So does ss_aipe_sem_path_sensitivity(), raising the planner's refusal.
  ridge <- recording_ridge(at_population)
  local_mocked_bindings(.lavaan_no_maximum = ridge$check)
  s <- caught(ss_aipe_sem_path_sensitivity(model,
                                           supposed_cov_matrix = cov_sem(pop_3f)$sigma_theta,
                                           which_path = "b",
                                           width = 0.3, G = 5))
  expect_match(conditionMessage(s), "\nab := a*b\n'\n", fixed = TRUE)
  expect_false(grepl("left out", conditionMessage(s), fixed = TRUE))
})

test_that("a finite negative error variance at the population moments is planned for", {
  skip_on_cran()
  S_fin <- stats::cov(triad_fin)
  r <- quiet(ss_aipe_sem_path(one_factor, true_cov_matrix = S_fin, width = 0.3,
                              which_path = "b"))
  expect_s3_class(r, "dmar_ss_aipe")
  fit <- quiet(lavaan::sem(one_factor, sample.cov = S_fin, sample.nobs = 1e6))
  pe <- lavaan::parameterEstimates(fit)
  expect_lt(pe$est[pe$op == "~~" & pe$lhs == "y1" & pe$rhs == "y1"], -0.2)
  expect_null(.lavaan_no_maximum(fit))
})

test_that("a Monte Carlo replication with no maximum is set aside as one that did not converge", {
  skip_on_cran()
  # The composite planners draw a fresh sample in its place: twenty
  # samples give the ten replications asked for, and no warning.
  sigma_2 <- matrix(c(1, 0.5, 0.5, 1), 2, 2,
                    dimnames = list(c("y1", "y2"), c("y1", "y2")))
  ridge <- every_other_sample()
  local_mocked_bindings(.lavaan_no_maximum = ridge$check)
  p <- with_warnings(ss_power_composite_sem("y2 ~ b*y1", true_cov_matrix = sigma_2,
                                            specified_N = 40, G = 10, seed = 113))
  expect_equal(p$value$value[p$value$term == "converged_replications"], 10)
  expect_length(p$warnings, 0L)
  expect_identical(ridge$state$samples, 20L)
  # So does ss_aipe_sem_path_sensitivity().
  ridge <- every_other_sample()
  local_mocked_bindings(.lavaan_no_maximum = ridge$check)
  set.seed(113)
  s <- quiet(ss_aipe_sem_path_sensitivity(two_factor, supposed_cov_matrix = Sigma_2f,
                                          which_path = "b",
                                          width = 0.3, specified_N = 100,
                                          G = 4))
  expect_identical(attr(s, "simulation")$returned, 4L)
  expect_identical(ridge$state$samples, 8L)
  # Every other data set was set aside and drawn again: four of the eight.
  expect_identical(attr(s, "simulation")$n_redrawn, 4L)
})

test_that("the functions that read only the chi square do not ask whether a fit has a maximum", {
  skip_on_cran()
  # Along a ridge the likelihood levels off, so the chi square, and the
  # RMSEA formed from it, do not depend on where the optimizer stopped:
  # ss_aipe_rmsea_sensitivity() and sem_exact_fit_bootstrap() never call
  # the check, and a stand-in that finds a ridge in every fit changes
  # nothing.
  ridge <- recording_ridge()
  set.seed(113)
  before <- quiet(ss_aipe_rmsea_sensitivity(width = 0.05, model = two_factor,
                                            true_cov_matrix = Sigma_2f, specified_N = 200, G = 6))
  local_mocked_bindings(.lavaan_no_maximum = ridge$check)
  set.seed(113)
  after <- quiet(ss_aipe_rmsea_sensitivity(width = 0.05, model = two_factor,
                                           true_cov_matrix = Sigma_2f, specified_N = 200, G = 6))
  expect_identical(after$value, before$value)
  expect_identical(attr(after, "simulation")$returned, 6L)
  expect_identical(ridge$seen$calls, 0L)
  # A bootstrap sample with no maximum is improper like any other with a
  # negative variance, counted or set aside as `improper` says.
  fit <- weak(2)
  b_keep <- quiet(sem_exact_fit_bootstrap(fit, B = 8, seed = 113))
  b_drop <- quiet(sem_exact_fit_bootstrap(fit, B = 8, seed = 113,
                                          improper = "drop"))
  expect_identical(ridge$seen$calls, 0L)
  bs <- attr(b_keep, "bootstrap")
  v <- function(x, term) x$value[x$term == term]
  expect_gt(sum(bs$status == "improper"), 0L)
  expect_true(all(is.finite(bs$chi_square[bs$status == "improper"])))
  expect_equal(v(b_keep, "B_used"), sum(bs$status %in% c("proper", "improper")))
  expect_equal(v(b_drop, "B_used"), sum(bs$status == "proper"))
  expect_equal(v(b_keep, "B_improper"), v(b_drop, "B_improper"))
  expect_equal(v(b_keep, "B_failed"),
               sum(bs$status %in% c("not_converged", "error")))
})

test_that("a structural model's observed outcomes and covariates are not error variances", {
  skip_on_cran()
  # lavaan carries x and y as latent variables of their own, which its
  # loadings hold and its covariance matrix of the latent variables leaves
  # out: the check forms that matrix with them, moves the loadings of M
  # with m1's error variance held further out while every variable keeps
  # its implied variance, leaves y's disturbance variance alone, and finds
  # the maximum at which m1's error variance is negative.
  fit <- quiet(lavaan::sem(sem_model, data = heywood_sem(), std.lv = TRUE))
  pt <- lavaan::parTable(fit)
  r <- which(pt$op == "~~" & pt$lhs == "m1" & pt$rhs == "m1")
  expect_lt(pt$est[r], 0)
  mats <- .no_max_matrices(fit)
  expect_setequal(mats[[1L]]$dummy, c("x", "y"))
  theta <- unclass(lavaan::lavInspect(fit, "est")$theta)
  expect_equal(c(mats[[1L]]$lambda %*% mats[[1L]]$cov_lv %*%
                   t(mats[[1L]]$lambda) + theta),
               c(unclass(mats[[1L]]$implied)), tolerance = 1e-8)
  st <- .no_max_start(pt, r, 10 * pt$est[r], mats)
  y_row <- which(pt$op == "~~" & pt$lhs == "y" & pt$rhs == "y")
  expect_identical(st[y_row], pt$est[y_row])
  expect_equal(diag(implied_sem(pt, st)), diag(implied_sem(pt, pt$est)),
               tolerance = 1e-8)
  expect_null(.lavaan_no_maximum(fit))
})

test_that("an indicator with a direct effect on it is left where it is", {
  skip_on_cran()
  # A direct effect of x on y2, as in a MIMIC model with a direct effect,
  # makes lavaan carry y2 as a latent variable of its own, with its
  # loading among the regressions and a zero in its row of the loadings.
  # With y1's error variance held further out, y2's loading starts at its
  # estimate (not at the zero of the loadings), y1's grows, y3's shrinks,
  # and every observed variable keeps its implied variance, whether the
  # factor's variance is fixed or y2's loading fixes its scale.
  d <- triad_fin
  set.seed(5)
  d$x <- stats::rnorm(200)
  d$y2 <- d$y2 + 0.3 * d$x
  for (model in c("f =~ y1 + y2 + y3\ny2 ~ x", "f =~ y2 + y1 + y3\ny2 ~ x")) {
    fit <- quiet(lavaan::sem(model, data = d,
                             std.lv = startsWith(model, "f =~ y1")))
    pt <- lavaan::parTable(fit)
    mats <- .no_max_matrices(fit)
    expect_setequal(mats[[1L]]$dummy, c("y2", "x"))
    expect_identical(unname(mats[[1L]]$lambda["y2", "f"]), 0)
    r <- which(pt$op == "~~" & pt$lhs == "y1" & pt$rhs == "y1")
    expect_lt(pt$est[r], 0)
    st <- .no_max_start(pt, r, 10 * pt$est[r], mats)
    row <- function(op, lhs, rhs) which(pt$op == op & pt$lhs == lhs &
                                          pt$rhs == rhs)
    expect_identical(st[row("=~", "f", "y2")], pt$est[row("=~", "f", "y2")])
    expect_identical(st[row("~~", "y2", "y2")], pt$est[row("~~", "y2", "y2")])
    expect_identical(st[row("~~", "f", "f")], pt$est[row("~~", "f", "f")])
    expect_gt(st[row("=~", "f", "y1")], pt$est[row("=~", "f", "y1")])
    expect_lt(st[row("=~", "f", "y3")], pt$est[row("=~", "f", "y3")])
    expect_equal(diag(implied_sem(pt, st)), diag(implied_sem(pt, pt$est)),
                 tolerance = 1e-8)
    expect_null(.lavaan_no_maximum(fit))
  }
  # The variance of such an indicator is not a candidate: with the direct
  # effect on y1, whose error variance is the negative one, the check does
  # not refit, the limit its comments and NEWS state.
  d1 <- triad_fin
  d1$x <- d$x
  d1$y1 <- d1$y1 + 0.3 * d1$x
  fit1 <- quiet(lavaan::sem("f =~ y1 + y2 + y3\ny1 ~ x", data = d1,
                            std.lv = TRUE))
  pt1 <- lavaan::parTable(fit1)
  expect_lt(pt1$est[pt1$op == "~~" & pt1$lhs == "y1" & pt1$rhs == "y1"], 0)
  local_mocked_bindings(.no_max_refit = function(...) stop("refit"))
  expect_null(.lavaan_no_maximum(fit1))
})

test_that("mediation_mbco() refuses a lavaan fit with no maximum, printing the code of its own fit or none for one it is given", {
  skip_on_cran()
  skip_if_not_installed("OpenMx")
  set.seed(9)
  d <- data.frame(x = stats::rnorm(200))
  latent <- 0.5 * d$x + stats::rnorm(200)
  d$y <- 0.4 * latent + stats::rnorm(200)
  d$m1 <- latent + stats::rnorm(200, sd = 0.6)
  d$m2 <- 0.8 * latent + stats::rnorm(200, sd = 0.6)
  d$m3 <- 0.7 * latent + stats::rnorm(200, sd = 0.6)
  ridge <- recording_ridge(variable = "m1")
  local_mocked_bindings(.lavaan_no_maximum = ridge$check)
  # A lavaan fit given as the model is refused as it stands, before any
  # OpenMx fit, with no code and with a different model, not bounds, as
  # the choice, since the OpenMx model carries no bounds.
  fit <- quiet(lavaan::sem(sem_model, data = d, meanstructure = TRUE))
  e <- caught(mediation_mbco(fit, x = "x", y = "y", ci_method = "wald"))
  expect_s3_class(e, c("dmar_no_maximum", "dmar_refusal", "error",
                       "condition"), exact = TRUE)
  expect_identical(e$fn, "mediation_mbco")
  expect_identical(e$where, "the lavaan fit given as model")
  expect_null(e$syntax)
  expect_null(e$lavaan_call)
  msg <- conditionMessage(e)
  expect_match(msg, paste0(
    "^mediation_mbco\\(\\) stopped: the fit failed because the likelihood ",
    "has no maximum \\(the lavaan fit given as model\\)\\."))
  expect_match(msg, paste(
    "DMAR fits the model with OpenMx, the estimation engine, starting from",
    "the estimates of this lavaan fit of the same model and data; see",
    "lavaan for more details."), fixed = TRUE)
  expect_true(endsWith(msg, "so the choice is a different model."))
  expect_false(grepl("There, bounds", msg, fixed = TRUE))
  # DMAR's own lavaan fit, from raw data: the call with the data a
  # placeholder, which reproduces the fit.
  r <- caught(mediation_mbco(sem_model, data = d, x = "x", y = "y",
                             ci_method = "wald"))
  expect_identical(r$where, "the lavaan fit the OpenMx fits start from")
  expect_identical(r$syntax, sem_model)
  expect_identical(gsub("\\s+", " ", r$lavaan_call), paste(
    "lavaan::sem(model, data = my_data, meanstructure = TRUE,",
    "fixed.x = FALSE, missing = \"listwise\", se = \"none\", warn = FALSE)"))
  expect_equal(run_code(r, my_data = d), ridge$seen$ll, tolerance = 1e-8)
  # From summary statistics: the call that fits them, the same fit as the
  # one to the data mediation_mbco() builds from them, with means named in
  # another order selected by name.
  S_d <- stats::cov(d)
  M_d <- rev(colMeans(d))
  s <- caught(mediation_mbco(sem_model, S = S_d, M = M_d, N = 200, x = "x",
                             y = "y", ci_method = "wald"))
  expect_match(gsub("\\s+", " ", s$lavaan_call),
               "^lavaan::sem\\(model, sample.cov = my_S\\[c\\(")
  expect_match(s$lavaan_call, "sample.nobs = 200", fixed = TRUE)
  expect_match(conditionMessage(s), paste(
    "with your covariance matrix in place of my_S and its means in place of",
    "my_M;"), fixed = TRUE)
  expect_equal(run_code(s, my_S = S_d, my_M = M_d), ridge$seen$ll,
               tolerance = 1e-8)
  # Without means, the data it builds have zero means, and so has the code.
  z <- caught(mediation_mbco(sem_model, S = S_d, N = 200, x = "x", y = "y",
                             ci_method = "wald"))
  expect_match(conditionMessage(z), paste(
    "with your covariance matrix in place of my_S and a zero mean for each",
    "of its rows in place of my_M;"), fixed = TRUE)
  expect_equal(run_code(z, my_S = S_d,
                        my_M = stats::setNames(rep(0, ncol(S_d)),
                                               colnames(S_d))),
               ridge$seen$ll, tolerance = 1e-8)
  # A covariance matrix given as a data frame is one group's: the call
  # fits it directly, passed as a matrix, since lavaan reads a list given
  # as sample.cov as one matrix per group.
  f <- caught(mediation_mbco(sem_model, S = as.data.frame(S_d), M = M_d,
                             N = 200, x = "x", y = "y", ci_method = "wald"))
  expect_match(gsub("\\s+", " ", f$lavaan_call),
               "^lavaan::sem\\(model, sample.cov = my_S\\[c\\(")
  expect_match(conditionMessage(f), paste(
    "with as.matrix() of your covariance matrix in place of my_S and its",
    "means in place of my_M;"), fixed = TRUE)
  expect_equal(run_code(f, my_S = as.matrix(as.data.frame(S_d)), my_M = M_d),
               ridge$seen$ll, tolerance = 1e-8)
})

test_that("mediation_mbco()'s refusal describes the data it builds from summary statistics exactly", {
  skip_on_cran()
  skip_if_not_installed("OpenMx")
  set.seed(9)
  d <- data.frame(x = stats::rnorm(200))
  latent <- 0.5 * d$x + stats::rnorm(200)
  d$y <- 0.4 * latent + stats::rnorm(200)
  d$m1 <- latent + stats::rnorm(200, sd = 0.6)
  d$m2 <- 0.8 * latent + stats::rnorm(200, sd = 0.6)
  d$m3 <- 0.7 * latent + stats::rnorm(200, sd = 0.6)
  d$w <- stats::rnorm(200)
  ridge <- recording_ridge(variable = "m1")
  local_mocked_bindings(.lavaan_no_maximum = ridge$check)
  # The rows .mbco_surrogate() draws for one group's summary statistics.
  built <- function(S, M, N) {
    set.seed(1)
    as.data.frame(MASS::mvrnorm(N, mu = M[colnames(S)], Sigma = S,
                                empirical = TRUE))
  }
  # Summary statistics of two groups: the data stacked in the order of S,
  # with the group's name in the column group.
  vars <- c("x", "y", "m1", "m2", "m3")
  S_g <- list(A = stats::cov(d[1:100, vars]), B = stats::cov(d[101:200, vars]))
  M_g <- list(A = colMeans(d[1:100, vars]), B = colMeans(d[101:200, vars]))
  g <- caught(mediation_mbco(sem_model, S = S_g, M = M_g,
                             N = c(A = 100, B = 100), x = "x", y = "y",
                             ci_method = "wald"))
  expect_match(conditionMessage(g), paste(
    "with the data mediation_mbco() builds from the summary statistics in",
    "place of my_data (for each group, set.seed(1) and then",
    "MASS::mvrnorm(N, mu = M[colnames(S)], Sigma = S, empirical = TRUE) on",
    "that group's N, M, and S, with mu = 0 without M, the groups stacked in",
    "the order of S with the group's name in the column group);"),
    fixed = TRUE)
  expect_match(g$lavaan_call, "group = \"group\"", fixed = TRUE)
  stacked <- rbind(cbind(built(S_g$A, M_g$A, 100), group = "A"),
                   cbind(built(S_g$B, M_g$B, 100), group = "B"))
  expect_equal(run_code(g, my_data = stacked), ridge$seen$ll,
               tolerance = 1e-8)
  # A product term x:w, from one group's summary statistics: the moments of
  # the product column depend on the rows drawn, which the seed fixes, and
  # the code reproduces the fit only from those rows.
  pm <- "M =~ m1 + m2 + m3\nM ~ x + w + x:w\ny ~ M + x"
  S_w <- stats::cov(d)
  M_w <- colMeans(d)
  p <- caught(mediation_mbco(pm, S = S_w, M = M_w, N = 200, x = "x", y = "y",
                             moderator = "w", ci_method = "wald"))
  expect_match(conditionMessage(p), paste(
    "with the data mediation_mbco() builds from the summary statistics in",
    "place of my_data (set.seed(1) and then MASS::mvrnorm(N, mu =",
    "M[colnames(S)], Sigma = S, empirical = TRUE), with mu = 0 without M),",
    "once x_by_w = x * w is added to it;"), fixed = TRUE)
  expect_match(p$syntax, "M ~ x + w + x_by_w", fixed = TRUE)
  rows <- built(S_w, M_w, 200)
  rows$x_by_w <- rows$x * rows$w
  expect_equal(run_code(p, my_data = rows), ridge$seen$ll, tolerance = 1e-8)
  # From raw data, the product column is the user's data's.
  q <- caught(mediation_mbco(pm, data = d, x = "x", y = "y",
                             moderator = "w", ci_method = "wald"))
  expect_match(conditionMessage(q), paste(
    "with your data in place of my_data, once x_by_w = x * w is added to",
    "it;"), fixed = TRUE)
  d$x_by_w <- d$x * d$w
  expect_equal(run_code(q, my_data = d), ridge$seen$ll, tolerance = 1e-8)
})

test_that("a proper Monte Carlo is the same with the check as without it", {
  skip_on_cran()
  # The three factor model of ss_power_composite_sem()'s help, at its
  # example's size: no replication has an error variance below zero, so
  # the check reads each parameter table and nothing more.
  with_check <- ss_power_composite_sem(three_factor, true_model = pop_3f,
                                       specified_N = 200, G = 20, seed = 113)
  local_mocked_bindings(.lavaan_no_maximum = function(fit, step = 10) NULL)
  without <- ss_power_composite_sem(three_factor, true_model = pop_3f,
                                    specified_N = 200, G = 20, seed = 113)
  expect_identical(with_check$value, without$value)
})
