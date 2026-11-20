## Shared fixtures -----------------------------------------------------
# The three-factor Holzinger and Swineford measurement model, grouped by
# school, is the workhorse for the multi-factor checks; the three
# visual-perception items are the one-factor convenience path.

hs_factors <- function() {
  list(
    visual = c("t1_visual_perception", "t2_cubes", "t3_paper_form_board"),
    verbal = c("t6_paragraph_comprehension", "t7_sentence",
               "t9_word_meaning"),
    speed  = c("t20_deduction", "t22_problem_reasoning",
               "t23_series_completion"))
}

hs_syntax <- function() {
  factors <- hs_factors()
  paste(vapply(names(factors), function(f)
    paste(f, "=~", paste(factors[[f]], collapse = " + ")), character(1L)),
    collapse = "\n")
}

# Four ordered categories per item, cut at the sample quartiles. Fully
# deterministic, so no seed is involved.
hs_discretized <- function(data, items) {
  for (item in items) {
    data[[item]] <- as.integer(cut(
      data[[item]],
      breaks = stats::quantile(data[[item]], c(0, .25, .5, .75, 1),
                               na.rm = TRUE),
      include.lowest = TRUE))
  }
  data
}

# Two categories per item, split at the sample median. Also deterministic.
hs_dichotomized <- function(data, items) {
  for (item in items) {
    data[[item]] <- as.integer(
      data[[item]] > stats::median(data[[item]], na.rm = TRUE))
  }
  data
}

# Runs an expression with its messages and warnings collected rather than
# printed, so a test can assert on exactly what the user would have been
# told, including that a particular condition was not raised.
hs_conditions <- function(expr) {
  messages <- character(0)
  warnings <- character(0)
  value <- withCallingHandlers(
    expr,
    message = function(cond) {
      messages <<- c(messages, conditionMessage(cond))
      invokeRestart("muffleMessage")
    },
    warning = function(cond) {
      warnings <<- c(warnings, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  list(value = value, messages = messages, warnings = warnings)
}


test_that("the continuous multi-factor ladder matches a hand-written lavaan ladder", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # eight ladder fits plus a lavaan oracle; anchors below run on CRAN
  data(holzinger_swineford)
  res <- measurement_invariance(holzinger_swineford, hs_factors(),
                                group = "school")
  expect_s3_class(res, "dmar_tbl")
  expect_identical(res$level,
                   c("configural", "metric", "scalar", "strict"))

  # The oracle: the same four models fit directly with lavaan::cfa(), one
  # group.equal set per rung, and lavTestLRT() for the step comparisons.
  model <- hs_syntax()
  equal <- list(configural = character(0),
                metric     = "loadings",
                scalar     = c("loadings", "intercepts"),
                strict     = c("loadings", "intercepts", "residuals"))
  oracle <- lapply(equal, function(eq) {
    args <- list(model = model, data = holzinger_swineford,
                 group = "school", meanstructure = TRUE)
    if (length(eq)) args$group.equal <- eq
    do.call(lavaan::cfa, args)
  })

  fm <- function(fit) lavaan::fitMeasures(
    fit, c("chisq", "df", "pvalue", "cfi", "rmsea"))
  for (i in seq_along(oracle)) {
    target <- fm(oracle[[i]])
    expect_equal(res$chi_square[i], unname(target[["chisq"]]),
                 tolerance = 1e-6)
    expect_equal(res$df[i], unname(target[["df"]]))
    expect_equal(res$p_chi_square[i], unname(target[["pvalue"]]),
                 tolerance = 1e-6)
    expect_equal(res$cfi[i], unname(target[["cfi"]]), tolerance = 1e-6)
    expect_equal(res$rmsea[i], unname(target[["rmsea"]]), tolerance = 1e-6)
  }
  for (i in 2:length(oracle)) {
    lrt <- lavaan::lavTestLRT(oracle[[i - 1L]], oracle[[i]])
    expect_equal(res$delta_chi_square[i], lrt[2L, "Chisq diff"],
                 tolerance = 1e-6)
    expect_equal(res$delta_df[i], lrt[2L, "Df diff"])
    expect_equal(res$p_value[i], lrt[2L, "Pr(>Chisq)"], tolerance = 1e-6)
    expect_equal(res$delta_cfi[i],
                 unname(fm(oracle[[i]])[["cfi"]] -
                          fm(oracle[[i - 1L]])[["cfi"]]),
                 tolerance = 1e-6)
    expect_equal(res$delta_rmsea[i],
                 unname(fm(oracle[[i]])[["rmsea"]] -
                          fm(oracle[[i - 1L]])[["rmsea"]]),
                 tolerance = 1e-6)
  }

  # Constraints accumulate, and the first row has no step comparison.
  expect_true(all(diff(res$df) > 0))
  expect_true(all(is.na(c(res$delta_chi_square[1L], res$delta_df[1L],
                          res$p_value[1L], res$delta_cfi[1L],
                          res$delta_rmsea[1L]))))
})


test_that("model syntax and the named-list form give the same ladder", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # fits the ladder twice over to compare two input routes
  data(holzinger_swineford)
  from_list <- measurement_invariance(holzinger_swineford, hs_factors(),
                                      group = "school")
  from_syntax <- measurement_invariance(holzinger_swineford,
                                        model = hs_syntax(),
                                        group = "school")
  num <- vapply(as.data.frame(unclass(from_list)), is.numeric, logical(1L))
  expect_equal(as.data.frame(unclass(from_syntax))[num],
               as.data.frame(unclass(from_list))[num], tolerance = 1e-12)
  expect_identical(attr(from_list, "model"), attr(from_syntax, "model"))
})


test_that("the one-factor items path returns what it always returned", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford)
  items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board")
  res <- measurement_invariance(holzinger_swineford, items = items,
                                group = "school")
  expect_identical(res$level,
                   c("configural", "metric", "scalar", "strict"))

  # The oracle is the call the pre-generalization implementation made:
  # "f =~ ..." fit by ML with meanstructure and nothing else set.
  model <- paste("f =~", paste(items, collapse = " + "))
  equal <- list(character(0), "loadings", c("loadings", "intercepts"),
                c("loadings", "intercepts", "residuals"))
  oracle <- lapply(equal, function(eq) {
    args <- list(model = model, data = holzinger_swineford,
                 group = "school", estimator = "ML", meanstructure = TRUE)
    if (length(eq)) args$group.equal <- eq
    do.call(lavaan::cfa, args)
  })
  fm <- function(fit) lavaan::fitMeasures(
    fit, c("chisq", "df", "pvalue", "cfi", "rmsea"))
  target <- t(vapply(oracle, function(f) as.numeric(fm(f)), numeric(5L)))
  expect_equal(res$chi_square, target[, 1L], tolerance = 1e-12)
  expect_equal(res$df, target[, 2L])
  expect_equal(res$p_chi_square, target[, 3L], tolerance = 1e-12)
  expect_equal(res$cfi, target[, 4L], tolerance = 1e-12)
  expect_equal(res$rmsea, target[, 5L], tolerance = 1e-12)
  for (i in 2:4) {
    lrt <- lavaan::lavTestLRT(oracle[[i - 1L]], oracle[[i]])
    expect_equal(res$delta_chi_square[i], lrt[2L, "Chisq diff"],
                 tolerance = 1e-12)
    expect_equal(res$delta_df[i], lrt[2L, "Df diff"])
    expect_equal(res$p_value[i], lrt[2L, "Pr(>Chisq)"], tolerance = 1e-12)
  }

  # Two further checks that do not amount to running the same lavaan call
  # twice. First, the degrees of freedom follow from counting moments, so
  # they can be asserted without fitting anything. Three indicators in two
  # groups give 3 variances, 3 covariances, and 3 means per group, 18 in
  # all. The configural model spends 9 per group (2 free loadings, the
  # first indicator being the marker; 1 factor variance; 3 residual
  # variances; 3 intercepts; the latent mean fixed at 0), so it is
  # saturated. Metric adds 2 loading equalities. Scalar adds 3 intercept
  # equalities and frees the latent mean in the second group, a net 2.
  # Strict adds 3 residual variance equalities.
  expect_equal(res$df, c(0, 2, 4, 7))
  expect_equal(res$delta_df, c(NA, 2, 2, 3))

  # Second, the same restrictions specified a different way: equality
  # across groups imposed by explicit per-group labels instead of by
  # lavaan's group.equal machinery, which is the mechanism
  # measurement_invariance() uses and the mechanism the oracle above
  # shares with it. Two routes to one set of restrictions must give one
  # set of statistics.
  loadings <- paste("f =~ 1*t1_visual_perception + c(l2, l2)*t2_cubes +",
                    "c(l3, l3)*t3_paper_form_board")
  labeled_metric <- lavaan::cfa(loadings, data = holzinger_swineford,
                                group = "school", meanstructure = TRUE)
  labeled_scalar <- lavaan::cfa(
    paste(loadings,
          "t1_visual_perception ~ c(i1, i1)*1",
          "t2_cubes ~ c(i2, i2)*1",
          "t3_paper_form_board ~ c(i3, i3)*1",
          "f ~ c(0, NA)*1", sep = "\n"),
    data = holzinger_swineford, group = "school", meanstructure = TRUE)
  expect_equal(res$chi_square[2L],
               unname(lavaan::fitMeasures(labeled_metric, "chisq")),
               tolerance = 1e-8)
  expect_equal(res$df[2L], unname(lavaan::fitMeasures(labeled_metric, "df")))
  expect_equal(res$chi_square[3L],
               unname(lavaan::fitMeasures(labeled_scalar, "chisq")),
               tolerance = 1e-8)
  expect_equal(res$df[3L], unname(lavaan::fitMeasures(labeled_scalar, "df")))
  expect_equal(res$p_value[3L],
               lavaan::lavTestLRT(labeled_metric,
                                  labeled_scalar)[2L, "Pr(>Chisq)"],
               tolerance = 1e-8)

  # Fixed anchors recorded from this fixture. They catch numerical drift
  # that a comparison against a freshly fitted model cannot see, since both
  # sides of that comparison would move together. These are values of this
  # ladder on these data, not constants published by Holzinger and
  # Swineford (1939).
  expect_equal(res$chi_square[2L], 1.18816045755321, tolerance = 1e-8)
  expect_equal(res$p_value[3L], 0.262780928631595, tolerance = 1e-8)
})


test_that("the returned attributes carry the non-numeric information", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a full ladder fit to read its attributes; anchors above run on CRAN
  data(holzinger_swineford)
  res <- measurement_invariance(holzinger_swineford, hs_factors(),
                                group = "school")
  fits <- attr(res, "fits")
  expect_type(fits, "list")
  expect_identical(names(fits), res$level)
  expect_true(all(vapply(fits, inherits, logical(1L), what = "lavaan")))
  expect_identical(attr(res, "estimator"), "ML")
  expect_false(attr(res, "ordered"))
  expect_identical(attr(res, "test"), "standard chi square difference test")
  expect_identical(attr(res, "fit_indices"), "standard")
  expect_true(grepl("visual =~", attr(res, "model"), fixed = TRUE))
  # Every value column stays numeric; metadata never becomes a column.
  expect_true(all(vapply(res[setdiff(names(res), "level")], is.numeric,
                         logical(1L))))
  # The stored fits are usable without refitting the ladder.
  expect_equal(unname(lavaan::fitMeasures(fits$metric, "chisq")),
               res$chi_square[2L], tolerance = 1e-12)
})


test_that("a robust estimator triggers the scaled difference test", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a second full ladder under MLR; the ML anchors above run on CRAN
  data(holzinger_swineford)
  res <- measurement_invariance(holzinger_swineford, hs_factors(),
                                group = "school", estimator = "MLR")
  expect_identical(attr(res, "fit_indices"), "robust")
  expect_match(attr(res, "test"), "^scaled chi square difference test")

  fits <- attr(res, "fits")
  lrt <- lavaan::lavTestLRT(fits$configural, fits$metric)
  expect_equal(res$delta_chi_square[2L], lrt[2L, "Chisq diff"],
               tolerance = 1e-10)
  # The scaled difference is not the difference of the scaled statistics;
  # that is what "scaled" means here, and the reported value follows
  # lavTestLRT rather than a subtraction.
  expect_false(isTRUE(all.equal(res$delta_chi_square[2L],
                                res$chi_square[2L] - res$chi_square[1L],
                                tolerance = 1e-6)))
  # The tabled model-level statistics are lavaan's scaled and robust ones.
  expect_equal(res$chi_square[1L],
               unname(lavaan::fitMeasures(fits$configural, "chisq.scaled")),
               tolerance = 1e-12)
  expect_equal(res$cfi[1L],
               unname(lavaan::fitMeasures(fits$configural, "cfi.robust")),
               tolerance = 1e-12)
})


test_that("group_partial frees the named parameter at every constrained rung", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two full ladders to compare a released constraint
  data(holzinger_swineford)
  full <- measurement_invariance(holzinger_swineford, hs_factors(),
                                 group = "school")
  partial <- measurement_invariance(holzinger_swineford, hs_factors(),
                                    group = "school",
                                    group_partial = "visual =~ t2_cubes")
  # Two groups, one loading equality released: the configural model has no
  # constraint to release, and each of the three constrained rungs loses
  # exactly one degree of freedom.
  expect_equal(full$df - partial$df, c(0, 1, 1, 1))
  expect_equal(full$chi_square[1L], partial$chi_square[1L],
               tolerance = 1e-10)
  # The loading really is free by group in the metric fit.
  est <- lavaan::parameterEstimates(attr(partial, "fits")$metric)
  loading <- est[est$op == "=~" & est$lhs == "visual" &
                   est$rhs == "t2_cubes", "est"]
  expect_length(loading, 2L)
  expect_false(isTRUE(all.equal(loading[1L], loading[2L],
                                tolerance = 1e-8)))
  # Without group_partial the same loading is constrained equal.
  est_full <- lavaan::parameterEstimates(attr(full, "fits")$metric)
  loading_full <- est_full[est_full$op == "=~" & est_full$lhs == "visual" &
                             est_full$rhs == "t2_cubes", "est"]
  expect_equal(loading_full[1L], loading_full[2L], tolerance = 1e-10)
})


test_that("ordered indicators get the thresholds ladder and the scaled test", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  data(holzinger_swineford)
  items <- unlist(hs_factors(), use.names = FALSE)
  ordered_data <- hs_discretized(holzinger_swineford, items)
  res <- suppressMessages(
    measurement_invariance(ordered_data, hs_factors(), group = "school",
                           ordered = TRUE))

  # The thresholds rung is on the ladder, and it comes before metric.
  expect_identical(res$level, c("configural", "thresholds", "metric",
                                "scalar", "strict"))
  expect_true(attr(res, "ordered"))
  expect_identical(attr(res, "estimator"), "WLSMV")
  expect_match(attr(res, "test"), "^scaled chi square difference test")
  expect_match(attr(res, "test"), "satorra")

  fits <- attr(res, "fits")
  expect_true(all(vapply(fits, function(f)
    lavaan::lavInspect(f, "converged"), logical(1L))))
  # Every rung genuinely adds constraints, including the strict rung: that
  # holds only because the theta parameterization is in force, since under
  # delta the residual variances are not free parameters.
  expect_true(all(diff(res$df) > 0))
  expect_equal(res$delta_df[-1L], diff(res$df))

  # The metric rung constrains thresholds as well as loadings.
  expect_identical(
    sort(lavaan::lavInspect(fits$metric, "options")$group.equal),
    sort(c("loadings", "thresholds")))
  expect_identical(lavaan::lavInspect(fits$thresholds,
                                      "options")$group.equal, "thresholds")
  expect_identical(lavaan::lavInspect(fits$strict, "options")$parameterization,
                   "theta")
})


test_that("the ordered ladder matches the pinned reparameterization anchors", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # the ladder itself is five WLSMV fits

  data(holzinger_swineford)
  items <- unlist(hs_factors(), use.names = FALSE)
  ordered_data <- hs_discretized(holzinger_swineford, items)
  res <- suppressMessages(
    measurement_invariance(ordered_data, hs_factors(), group = "school",
                           ordered = TRUE))

  # The anchors are an external build of the same ladder with the Wu and
  # Estabrook (2016) identification written out explicitly and a different
  # factor identification (std.lv rather than lavaan's marker variable).
  # The models are equivalent reparameterizations, so the test statistics
  # and degrees of freedom must agree; agreement is to estimation-
  # convergence precision rather than to machine precision.
  # Pinned from semTools::measEq.syntax (semTools 0.5.7, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  oracle_df <- c(48, 57, 63, 69, 78)
  oracle_chi_square <- c(58.32879773268559, 70.51320723129015,
                         66.75343278406135, 86.67680449597691,
                         95.73938512700711)
  oracle_delta_chi_square <- c(NA, 13.6914002836, 1.6053579187,
                               14.8920154517, 10.3239423871)
  oracle_delta_df <- c(NA, 9, 6, 6, 9)
  oracle_p_value <- c(NA, 0.1337339440753325, 0.9521914355019727,
                      0.02111342558031126, 0.3249023590393668)

  for (i in seq_along(oracle_df)) {
    expect_equal(res$df[i], oracle_df[i])
    expect_equal(res$chi_square[i], oracle_chi_square[i], tolerance = 1e-4)
  }
  for (i in 2:length(oracle_df)) {
    expect_equal(res$delta_chi_square[i], oracle_delta_chi_square[i],
                 tolerance = 1e-4)
    expect_equal(res$delta_df[i], oracle_delta_df[i])
    expect_equal(res$p_value[i], oracle_p_value[i], tolerance = 1e-4)
  }
})


test_that("full information maximum likelihood handles incomplete cases", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two ladders refit under FIML; the anchors above run on CRAN
  data(holzinger_swineford)
  items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board")
  incomplete <- holzinger_swineford
  set.seed(113)
  incomplete$t2_cubes[sample(nrow(incomplete), 20L)] <- NA
  res <- measurement_invariance(incomplete, items = items, group = "school",
                                missing = "fiml")
  expect_identical(res$level,
                   c("configural", "metric", "scalar", "strict"))
  expect_true(all(is.finite(res$chi_square)))
  # Every case contributes, unlike listwise deletion.
  expect_equal(sum(lavaan::lavInspect(attr(res, "fits")$metric, "nobs")),
               nrow(incomplete))
  listwise <- measurement_invariance(incomplete, items = items,
                                     group = "school", missing = "listwise")
  expect_equal(sum(lavaan::lavInspect(attr(listwise, "fits")$metric,
                                      "nobs")),
               nrow(incomplete) - 20L)
})

test_that("under FIML an auxiliary variable's mean and variance stay free in every group", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # four ladders refit under FIML
  # lavaan adds an auxiliary variable (aux) to a FIML model as a saturated
  # correlate, its mean and variance ordinary parameters of each group, so
  # group.equal reached them: through version 1.0.1 the scalar rung equated
  # the means of general information across the schools (44.8 and 36.6)
  # and the strict rung its variances, and the scalar test was 57.59 on 5
  # df against 24.03 on 4 without the auxiliary. The rung differences are
  # now those of the same ladder with the auxiliary's moments free.
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  f <- list(spatial = items[1:3], verbal = items[4:6])
  d <- holzinger_swineford[c(items, "t5_general_information", "school")]
  set.seed(1)
  d$t7_sentence[sample(301, 60)] <- NA
  quiet <- function(x) suppressMessages(suppressWarnings(x))
  a <- quiet(measurement_invariance(d, f, group = "school", missing = "fiml",
                                    aux = "t5_general_information"))
  b <- quiet(measurement_invariance(d[c(items, "school")], f,
                                    group = "school", missing = "fiml"))
  expect_identical(a$df, b$df)
  expect_identical(a$delta_df, b$delta_df)
  expect_identical(a$df, c(16, 20, 24, 30))
  free <- c("t5_general_information ~ 1",
            "t5_general_information ~~ t5_general_information")
  lad <- lapply(list(character(0), "loadings", c("loadings", "intercepts"),
                     c("loadings", "intercepts", "residuals")), function(ge) {
    args <- list(model = attr(a, "model"), data = d, group = "school",
                 missing = "ml", aux = "t5_general_information",
                 meanstructure = TRUE, group.partial = free)
    if (length(ge)) args$group.equal <- ge
    quiet(do.call(lavaan::cfa, args))
  })
  lrt <- vapply(2:4, function(i) {
    unlist(lavaan::lavTestLRT(lad[[i - 1L]], lad[[i]])[2L, c(
      "Chisq diff", "Df diff", "Pr(>Chisq)")])
  }, numeric(3))
  expect_equal(a$delta_chi_square[2:4], lrt[1, ], tolerance = 1e-8)
  expect_equal(a$delta_df[2:4], lrt[2, ])
  expect_equal(a$p_value[2:4], lrt[3, ], tolerance = 1e-8)
  expect_equal(a$delta_chi_square[3], 23.875175, tolerance = 1e-6)
  fm <- lapply(lad, lavaan::fitMeasures, c("cfi.robust", "rmsea.robust"))
  expect_equal(a$cfi, vapply(fm, `[[`, numeric(1), 1L), tolerance = 1e-8)
  expect_equal(a$rmsea, vapply(fm, `[[`, numeric(1), 2L), tolerance = 1e-8)
  # The auxiliary's mean is estimated in each group, near its two sample
  # means (44.8 and 36.6).
  pt <- lavaan::parTable(attr(a, "fits")$scalar)
  aux_mean <- pt[pt$lhs == "t5_general_information" & pt$op == "~1", ]
  expect_identical(aux_mean$label, c("", ""))
  expect_gt(abs(diff(aux_mean$est)), 7)
  # A group_partial of the caller's is kept beside them.
  p <- quiet(measurement_invariance(d, f, group = "school", missing = "fiml",
                                    aux = "t5_general_information",
                                    group_partial = "verbal =~ t7_sentence",
                                    levels = c("configural", "metric")))
  expect_identical(p$delta_df[2], 3)
  expect_setequal(lavaan::lavInspect(attr(p, "fits")$metric,
                                     "options")$group.partial,
                  gsub(" ", "", c("verbal =~ t7_sentence", free)))
})

test_that("lavaan's group.partial passed through ... joins group_partial, and group.equal is refused there", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # several ladders refit under FIML
  # The help says that ... goes to lavaan::cfa() at every rung. A rung
  # passes one group.partial to lavaan, with the caller's terms and, under
  # FIML, the auxiliary's, so a group.partial in ... is merged into it:
  # lavaan stops on a second one with an error that a name is duplicated,
  # which it did through version 1.0.1 for a call that gave group_partial
  # as well.
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  f <- list(spatial = items[1:3], verbal = items[4:6])
  d <- holzinger_swineford[c(items, "t5_general_information", "school")]
  set.seed(1)
  d$t7_sentence[sample(301, 60)] <- NA
  quiet <- function(x) suppressMessages(suppressWarnings(x))
  free <- c("t5_general_information ~ 1",
            "t5_general_information ~~ t5_general_information")
  a <- quiet(measurement_invariance(d, f, group = "school", missing = "fiml",
                                    aux = "t5_general_information",
                                    group.partial = "t2_cubes ~ 1"))
  b <- quiet(measurement_invariance(d, f, group = "school", missing = "fiml",
                                    aux = "t5_general_information",
                                    group_partial = "t2_cubes ~ 1"))
  expect_equal(a$delta_chi_square, b$delta_chi_square, tolerance = 1e-10)
  expect_identical(a$df, c(16, 20, 23, 29))
  expect_equal(a$delta_chi_square[3], 15.093130, tolerance = 1e-6)
  expect_setequal(lavaan::lavInspect(attr(a, "fits")$scalar,
                                     "options")$group.partial,
                  gsub(" ", "", c("t2_cubes ~ 1", free)))
  # Both spellings at once free the terms of both.
  both <- quiet(measurement_invariance(d[c(items, "school")], f,
                                       group = "school", missing = "fiml",
                                       group_partial = "t2_cubes ~ 1",
                                       group.partial = "t4_lozenges ~ 1"))
  one <- quiet(measurement_invariance(d[c(items, "school")], f,
                                      group = "school", missing = "fiml",
                                      group_partial = c("t2_cubes ~ 1",
                                                        "t4_lozenges ~ 1")))
  expect_equal(both$delta_chi_square, one$delta_chi_square,
               tolerance = 1e-10)
  expect_identical(both$df, c(16, 20, 22, 28))
  expect_error(measurement_invariance(d, f, group = "school",
                                      group.partial = NA_character_),
               "'group.partial' must be a character vector")
  expect_error(measurement_invariance(d, f, group = "school",
                                      group.partial = 42),
               "'group.partial' must be a character vector")
  # An empty set of terms, in either spelling, frees nothing, as lavaan
  # reads it, so the ladder is the one fitted without the argument (and
  # the auxiliary's terms are still freed). Through the first version of
  # the merge, group.partial = NULL or character(0) in ... stopped with the
  # error above, and so did the idiom group.partial = if (cond) terms else
  # NULL, which DMAR 1.0.1 passed to lavaan.
  none <- quiet(measurement_invariance(d, f, group = "school",
                                       missing = "fiml",
                                       aux = "t5_general_information"))
  expect_identical(none$df, c(16, 20, 24, 30))
  free_nothing <- list(list(group.partial = NULL),
                       list(group.partial = character(0)),
                       list(group.partial = ""),
                       list(group.partial = if (FALSE) "t2_cubes ~ 1"),
                       list(group_partial = character(0)))
  for (arg in free_nothing) {
    e <- quiet(do.call(measurement_invariance,
                       c(list(d, f, group = "school", missing = "fiml",
                              aux = "t5_general_information"), arg)))
    expect_equal(e$delta_chi_square, none$delta_chi_square,
                 tolerance = 1e-10, info = deparse(arg))
    expect_identical(e$df, none$df, info = deparse(arg))
    expect_setequal(lavaan::lavInspect(attr(e, "fits")$scalar,
                                       "options")$group.partial,
                    gsub(" ", "", free))
  }
  # An empty one beside group_partial leaves group_partial's terms.
  e <- quiet(measurement_invariance(d, f, group = "school", missing = "fiml",
                                    aux = "t5_general_information",
                                    group_partial = "t2_cubes ~ 1",
                                    group.partial = NULL))
  expect_equal(e$delta_chi_square, b$delta_chi_square, tolerance = 1e-10)
  # group.equal is what the ladder sets; through version 1.0.1 it
  # constrained the configural rung (25.62 on 20 df, against 15.27 on 16)
  # and stopped the metric rung with lavaan's error.
  expect_error(measurement_invariance(d, f, group = "school",
                                      group.equal = "loadings",
                                      levels = "configural"),
               "'group.equal' is set by the ladder", fixed = TRUE)
  expect_error(measurement_invariance(d, f, group = "school",
                                      group.equal = "loadings"),
               "Choose the rungs with 'levels'", fixed = TRUE)
})


test_that("under FIML the fit indices are the corrected indices of Zhang and Savalei (2023)", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two ladders refit under FIML
  # The 1.1.0 review: measurement_invariance(missing = "fiml") reported
  # lavaan's standard cfi and rmsea, which move with the amount of missing
  # data, while fa_confirmatory() reports the corrected indices for the same kind of
  # fit. On the verbal and deduction tests with 90 sentence completion
  # scores deleted at random, the configural CFI was .9936 and the RMSEA
  # .0426, and the metric rung's delta_cfi was -.0106.
  hs <- holzinger_swineford
  set.seed(113)
  hs$t7_sentence[sample(nrow(hs), 90L)] <- NA
  f <- list(verbal = c("t6_paragraph_comprehension", "t7_sentence",
                       "t9_word_meaning"),
            deduction = c("t20_deduction", "t22_problem_reasoning",
                          "t23_series_completion"))
  mi <- measurement_invariance(hs, f, group = "school", missing = "fiml")
  fm <- lapply(attr(mi, "fits"), lavaan::fitMeasures)
  pick <- function(m) vapply(fm, function(x) x[[m]], numeric(1L),
                             USE.NAMES = FALSE)
  expect_equal(mi$cfi, pick("cfi.robust"), tolerance = 1e-12)
  expect_equal(mi$rmsea, pick("rmsea.robust"), tolerance = 1e-12)
  expect_equal(mi$cfi, c(0.9929797510, 0.9825999520, 0.9775428135,
                         0.9731335112), tolerance = 1e-8)
  expect_equal(mi$rmsea, c(0.04733930771, 0.06666008494, 0.06913177197,
                           0.06763172511), tolerance = 1e-8)
  expect_equal(mi$delta_cfi[-1L], diff(pick("cfi.robust")),
               tolerance = 1e-12)
  expect_equal(mi$delta_rmsea[-1L], diff(pick("rmsea.robust")),
               tolerance = 1e-12)
  # The standard indices it reported before, and the chi square column,
  # which stays the FIML likelihood ratio statistic.
  expect_equal(pick("cfi")[1:2], c(0.9936235117, 0.9830425652),
               tolerance = 1e-8)
  expect_equal(pick("rmsea")[1], 0.04261202471, tolerance = 1e-8)
  expect_equal(mi$chi_square, pick("chisq"), tolerance = 1e-12)
  expect_identical(attr(mi, "fit_indices"), "robust")
  expect_identical(generics::glance(mi)$fit_indices, "robust")
  line <- paste0(
    "cfi, rmsea: the FIML corrected indices of Zhang and Savalei (2023), ",
    "version 6 (normal theory), which approximate the values complete ",
    "data would give and are not computed from the chi_square column; ",
    "delta_cfi and delta_rmsea compare them")
  expect_identical(unname(attr(mi, "fit_columns")), line)
  expect_output(print(mi), line, fixed = TRUE)
  skip_if_not_installed("knitr")
  expect_match(as.character(knitr::knit_print(mi)),
               "version 6 (normal theory)", fixed = TRUE)

  # MLR under FIML already read the robust indices, version 3; the
  # numbers are unchanged, and the line now names them.
  mr <- measurement_invariance(hs, f, group = "school", missing = "fiml",
                               estimator = "MLR",
                               levels = c("configural", "metric"))
  fr <- lapply(attr(mr, "fits"), lavaan::fitMeasures)
  expect_equal(mr$cfi, vapply(fr, function(x) x[["cfi.robust"]],
                              numeric(1L), USE.NAMES = FALSE),
               tolerance = 1e-12)
  expect_identical(attr(mr, "fit_indices"), "robust")
  expect_match(attr(mr, "fit_columns"), "version 3 (no normality",
               fixed = TRUE)

  # Complete data and listwise deletion keep the standard indices and add
  # no line.
  lw <- measurement_invariance(hs, f, group = "school",
                               levels = "configural", missing = "listwise")
  expect_identical(attr(lw, "fit_indices"), "standard")
  expect_null(attr(lw, "fit_columns"))
})

test_that("under FIML every rung falls back together when a corrected index is missing", {
  # When lavaan cannot compute the corrected indices at some rung, every
  # rung reports one other version, so the differences compare like with
  # like, and the line says which.
  fms <- list(c(cfi = .99, rmsea = .04, cfi.robust = .985,
                rmsea.robust = .05, cfi.scaled = .98, rmsea.scaled = .06),
              c(cfi = .98, rmsea = .05, cfi.robust = NA,
                rmsea.robust = .06, cfi.scaled = .97, rmsea.scaled = .07))
  ok <- .measurement_invariance_indices(fms[c(1, 1)], scaled = FALSE,
                                        "fiml")
  expect_identical(ok$measures, c("cfi.robust", "rmsea.robust"))
  expect_identical(ok$version, "robust")
  std <- .measurement_invariance_indices(fms, scaled = FALSE, "fiml")
  expect_identical(std$measures, c("cfi", "rmsea"))
  expect_identical(std$version, "standard")
  expect_match(std$lines, paste0(
    "lavaan could not compute the FIML corrected indices of Zhang and ",
    "Savalei (2023) at every rung, so these columns hold its standard ",
    "versions, which change with the amount of missing data"), fixed = TRUE)
  sc <- .measurement_invariance_indices(fms, scaled = TRUE, "fiml")
  expect_identical(sc$measures, c("cfi.scaled", "rmsea.scaled"))
  expect_identical(sc$version, "scaled")
  expect_match(sc$lines, "its scaled versions, which put the scaled",
               fixed = TRUE)
})


test_that("a scaled ladder names the versions it reports at every rung", {
  # Through 1.0.1 the names were read from the first rung alone: when
  # lavaan had no robust cfi and rmsea for a fit with a scaled statistic,
  # the columns held the standard ones while fit_indices said "robust",
  # and when it had no scaled statistic, chi_square held the standard one
  # while model_test said "scaled". One version holds every rung, and the
  # attributes name it.
  rung <- function(df, robust = TRUE, scaled_test = TRUE, cfi_robust = .985,
                   chisq = if (df == 0) 0 else 30) {
    c(chisq = chisq, df = df, pvalue = .02,
      if (scaled_test) c(chisq.scaled = 25, df.scaled = df,
                         pvalue.scaled = .05),
      cfi = .98, rmsea = .05, cfi.scaled = .97, rmsea.scaled = .06,
      if (robust) c(cfi.robust = cfi_robust, rmsea.robust = .055))
  }
  ok <- .measurement_invariance_indices(list(rung(16), rung(20)),
                                        scaled = TRUE, "ml")
  expect_identical(ok$version, "robust")
  expect_null(ok$lines)
  # No robust indices at a rung with degrees of freedom: every rung of a
  # maximum likelihood ladder of continuous items falls back to the
  # standard indices, which estimate the same population values as the
  # robust ones, never to the scaled ones (Savalei, 2018).
  for (fms in list(list(rung(16, robust = FALSE), rung(20, robust = FALSE)),
                   list(rung(16), rung(20, robust = FALSE)),
                   list(rung(16), rung(20, cfi_robust = NA)))) {
    std <- .measurement_invariance_indices(fms, scaled = TRUE, "ml")
    expect_identical(std$measures, c("cfi", "rmsea"))
    expect_identical(std$version, "standard")
    expect_identical(unname(std$lines), paste0(
      "cfi, rmsea: lavaan did not compute its robust versions at every ",
      "rung, so these columns hold its standard versions, which estimate ",
      "the same population values with the bias corrected as if the data ",
      "were normal (Savalei, 2018); delta_cfi and delta_rmsea compare them"))
  }
  # A ladder of ordered items (or by least squares) falls back to the
  # scaled indices, since its standard statistic is not the basis of an
  # index.
  ord <- .measurement_invariance_indices(list(rung(16), rung(20, robust = FALSE)),
                                         scaled = TRUE, "other")
  expect_identical(ord$measures, c("cfi.scaled", "rmsea.scaled"))
  expect_identical(ord$version, "scaled")
  expect_identical(unname(ord$lines), paste0(
    "cfi, rmsea: lavaan did not compute its robust versions at every ",
    "rung, so these columns hold its scaled versions, which put the scaled ",
    "statistic into the standard formulas; delta_cfi and delta_rmsea ",
    "compare them"))
  # Without a scaled statistic the standard indices are the only version.
  plain <- .measurement_invariance_indices(list(rung(16)), scaled = FALSE,
                                           "ml")
  expect_identical(plain$version, "standard")
  expect_null(plain$lines)
  # A two-stage ladder falls back as fa_confirmatory() does, and its line replaces
  # the two-stage line.
  ts <- .measurement_invariance_indices(list(rung(16), rung(20, robust = FALSE)),
                                        scaled = TRUE, "two_stage")
  expect_identical(ts$version, "scaled")
  expect_match(ts$lines, paste0("lavaan could not compute the two-stage ",
                                "corrected indices of Zhang and Savalei ",
                                "(2023) at every rung"), fixed = TRUE)
  # A saturated rung does not decide: lavaan can leave its robust CFI NA
  # from a ratio of zeros, and the other rungs keep the robust indices
  # (under FIML the strict rule sent every rung to the standard ones).
  sat <- list(rung(0, cfi_robust = NA), rung(2), rung(4))
  expect_identical(.measurement_invariance_indices(sat, TRUE, "ml")$version,
                   "robust")
  expect_identical(.measurement_invariance_indices(sat, TRUE, "other")$version,
                   "robust")
  expect_identical(.measurement_invariance_indices(sat, FALSE, "fiml")$version,
                   "robust")
  expect_identical(
    .measurement_invariance_indices(list(rung(0), rung(2, cfi_robust = NA)),
                                    FALSE, "fiml")$version, "standard")
  # A rung with no degrees of freedom whose statistic is not zero (an
  # inequality bound that binds) does not fit exactly, and decides as any
  # other rung does.
  misfit <- list(rung(0, cfi_robust = NA, chisq = 30), rung(2), rung(4))
  expect_identical(.measurement_invariance_indices(misfit, TRUE, "ml")$version,
                   "standard")
  # It decides the indices through its CFI alone, since its RMSEA divides
  # by the degrees of freedom and is undefined in every version, and it
  # does not decide the test, which no version defines on zero degrees of
  # freedom (lavaan leaves every p-value NA there). Through the first
  # version of the rule above an RMSEA or p-value of NA at such a rung sent
  # every rung to the standard version.
  r0 <- rung(0, chisq = 30)
  r0[c("rmsea.robust", "rmsea.scaled", "pvalue", "pvalue.scaled")] <- NA
  expect_identical(
    .measurement_invariance_indices(list(r0, rung(2), rung(4)), TRUE,
                                    "ml")$version, "robust")
  expect_identical(
    .measurement_invariance_test_columns(list(r0, rung(2), rung(4)),
                                         TRUE)$version, "scaled")
  # The test columns: the scaled statistic when every rung has it.
  tc <- .measurement_invariance_test_columns(list(rung(16), rung(20)), TRUE)
  expect_identical(tc$version, "scaled")
  expect_identical(tc$measures, c("chisq.scaled", "df.scaled",
                                  "pvalue.scaled"))
  tc <- .measurement_invariance_test_columns(
    list(rung(16), rung(20, scaled_test = FALSE)), TRUE)
  expect_identical(tc$version, "standard")
  expect_identical(tc$measures, c("chisq", "df", "pvalue"))
  expect_identical(
    .measurement_invariance_test_columns(list(rung(16)), FALSE)$version,
    "standard")
})


test_that("a rung with no degrees of freedom that does not fit exactly has no rmsea", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two ladders of four rungs
  # One factor of three items in two groups of 200, correlated .8, .8, and
  # .5 in the first group, whose one factor solution is improper. With
  # bounds = "pos.var" the bound binds, which lavaan's count of the degrees
  # of freedom leaves out, so the configural rung has a statistic of 32.8 on
  # 0 df and does not fit exactly. Its cfi is formed from the statistic (the
  # 1.1.0 repairs first set it to 1), and its rmsea, which divides by the
  # degrees of freedom, is NA, as is the metric delta_rmsea; through 1.0.1
  # they were lavaan's 0, the value of exact fit, and .302.
  draw <- function(R, n, seed) {
    set.seed(seed)
    X <- matrix(rnorm(3 * n), n)
    X <- scale(X) %*% solve(chol(cor(X))) %*% chol(R)
    as.data.frame(X)
  }
  R1 <- matrix(c(1, .8, .8, .8, 1, .5, .8, .5, 1), 3, 3)
  R2 <- matrix(c(1, .6, .5, .6, 1, .55, .5, .55, 1), 3, 3)
  d <- rbind(draw(R1, 200, 1), draw(R2, 200, 2))
  names(d) <- c("a", "b", "c")
  d$g <- rep(c("x", "y"), each = 200)
  quiet <- function(x) suppressMessages(suppressWarnings(x))
  ml <- quiet(measurement_invariance(d, "f =~ a + b + c", group = "g",
                                     bounds = "pos.var"))
  expect_identical(ml$df[1L], 0)
  expect_equal(ml$chi_square[1L], 32.79, tolerance = 1e-3)
  cfi0 <- unname(lavaan::fitMeasures(attr(ml, "fits")[[1L]], "cfi"))
  expect_equal(ml$cfi[1L], cfi0, tolerance = 1e-12)
  expect_lt(ml$cfi[1L], 0.95)
  expect_true(is.na(ml$rmsea[1L]))
  expect_true(is.na(ml$delta_rmsea[2L]))
  expect_equal(ml$delta_cfi[2L], ml$cfi[2L] - ml$cfi[1L], tolerance = 1e-12)
  expect_true(all(is.finite(ml$rmsea[-1L])))
  line <- attr(ml, "fit_columns")
  expect_length(line, 1L)
  expect_match(line, paste0("^rmsea: NA at the configural rung, which has ",
                            "zero degrees of freedom and a statistic that is ",
                            "not zero"))
  printed <- paste(capture.output(print(ml)), collapse = " ")
  expect_match(printed, "rmsea: NA at the configural rung", fixed = TRUE)
  # Under MLR the configural rung does not decide the test, which no
  # version defines on zero degrees of freedom: the chi_square column is the
  # scaled statistic at every rung, as the scaled difference tests are.
  # Through the first version of the rule above it was the standard
  # statistic, under a model_test attribute of "standard". lavaan leaves the
  # robust CFI of that rung NA, and its CFI decides the indices, so every
  # rung reports the standard ones.
  mlr <- quiet(measurement_invariance(d, "f =~ a + b + c", group = "g",
                                      estimator = "MLR", bounds = "pos.var"))
  scaled <- vapply(attr(mlr, "fits")[-1L], function(x) {
    unname(lavaan::fitMeasures(x, "chisq.scaled"))
  }, numeric(1L))
  expect_equal(mlr$chi_square[-1L], unname(scaled), tolerance = 1e-12)
  expect_identical(attr(mlr, "model_test"), "scaled")
  expect_match(attr(mlr, "test"), "^scaled")
  expect_identical(attr(mlr, "fit_indices"), "standard")
  expect_equal(mlr$cfi[1L], cfi0, tolerance = 1e-12)
  expect_true(is.na(mlr$rmsea[1L]))
  expect_true(is.na(mlr$delta_rmsea[2L]))
  expect_true(any(startsWith(attr(mlr, "fit_columns"),
                             "rmsea: NA at the configural rung")))
})


test_that("a ladder with the mean and variance adjusted statistic reports the robust indices DMAR computes", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # three ladders of two rungs
  # With the older mean and variance adjusted statistic lavaan computes no
  # robust cfi and rmsea. Through 1.0.1 the ladder reported the standard
  # ones while fit_indices said "robust"; the 1.1.0 review relabeled them
  # "standard". The robust indices do not depend on which corrected
  # statistic accompanies them (Savalei, 2018), and DMAR computes them,
  # equal to those of the same ladder under MLM.
  f <- list(verbal = c("t6_paragraph_comprehension", "t7_sentence",
                       "t9_word_meaning"),
            deduction = c("t20_deduction", "t22_problem_reasoning",
                          "t23_series_completion"))
  mi <- measurement_invariance(holzinger_swineford, f, group = "school",
                               test = "mean.var.adjusted",
                               levels = c("configural", "metric"))
  fm <- lapply(attr(mi, "fits"), lavaan::fitMeasures)
  pick <- function(m) vapply(fm, function(x) x[[m]], numeric(1L),
                             USE.NAMES = FALSE)
  expect_false(any(c("cfi.robust", "rmsea.robust") %in% names(fm[[1L]])))
  expect_identical(attr(mi, "fit_indices"), "robust")
  expect_identical(generics::glance(mi)$fit_indices, "robust")
  expect_identical(attr(mi, "model_test"), "scaled")
  expect_equal(mi$chi_square, pick("chisq.scaled"), tolerance = 1e-12)
  expect_equal(mi$df, pick("df.scaled"), tolerance = 1e-12)
  mm <- measurement_invariance(holzinger_swineford, f, group = "school",
                               estimator = "MLM",
                               levels = c("configural", "metric"))
  expect_equal(mi$cfi, mm$cfi, tolerance = 1e-10)
  expect_equal(mi$rmsea, mm$rmsea, tolerance = 1e-10)
  expect_equal(mi$cfi, c(0.9927891929, 0.9794338307), tolerance = 1e-8)
  expect_equal(mi$delta_cfi[2L], -0.01335536224, tolerance = 1e-7)
  # 1.0.1 reported the standard indices, and the scaled ones differ again.
  expect_equal(pick("cfi"), c(0.9915882855, 0.9788618835), tolerance = 1e-8)
  expect_equal(pick("cfi.scaled")[1L], 0.9628364532, tolerance = 1e-8)
  expect_match(attr(mi, "fit_columns"), "computed by DMAR", fixed = TRUE)
  expect_output(print(mi), "computed by DMAR", fixed = TRUE)
  # Were the robust indices not to be had, the ladder would fall back to
  # the standard indices, never to the scaled ones (Savalei, 2018).
  ms <- with_mocked_bindings(
    measurement_invariance(holzinger_swineford, f, group = "school",
                           test = "mean.var.adjusted",
                           levels = c("configural", "metric")),
    .sem_robust_mva = function(...) numeric(0))
  expect_identical(attr(ms, "fit_indices"), "standard")
  expect_equal(ms$cfi, pick("cfi"), tolerance = 1e-12)
  expect_match(attr(ms, "fit_columns"),
               "lavaan did not compute its robust versions at every rung",
               fixed = TRUE)
  # A saturated configural rung does not decide the version: lavaan
  # leaves the robust CFI of this one NA, and the other rungs keep the
  # robust indices, as through 1.0.1.
  m <- suppressWarnings(measurement_invariance(
    holzinger_swineford, items = f$verbal, group = "school",
    estimator = "MLR", levels = c("configural", "metric", "scalar")))
  fr <- lapply(attr(m, "fits"), lavaan::fitMeasures,
               c("df", "cfi.robust", "rmsea.robust"))
  expect_identical(fr$configural[["df"]], 0)
  expect_true(is.na(fr$configural[["cfi.robust"]]))
  expect_identical(attr(m, "fit_indices"), "robust")
  expect_equal(m$cfi[2:3], vapply(fr[2:3], function(x) x[["cfi.robust"]],
                                  numeric(1L), USE.NAMES = FALSE),
               tolerance = 1e-12)
  expect_null(attr(m, "fit_columns"))
  # The saturated rung's statistic is zero in every version, so its CFI is
  # 1 whatever the scaling constant (0/0 there), and its RMSEA keeps
  # lavaan's 0; through 1.0.1 its CFI, and the metric rung's delta_cfi,
  # were NA.
  expect_identical(m$cfi[1L], 1)
  expect_identical(m$rmsea[1L], 0)
  expect_equal(m$delta_cfi[2L], fr$metric[["cfi.robust"]] - 1,
               tolerance = 1e-12)
  expect_equal(m$delta_cfi[2L], -0.009713357, tolerance = 1e-6)
  expect_equal(m$delta_rmsea[2L], fr$metric[["rmsea.robust"]],
               tolerance = 1e-12)
  # Nor does it decide the version of a ladder with the mean and variance
  # adjusted statistic, whose robust indices DMAR computes at every rung,
  # the saturated one included. That ladder equals the MLM ladder; the
  # first version of the 1.1.0 fix computed nothing at the saturated rung,
  # which sent every rung to the standard indices (a metric delta CFI of
  # -.0109 against the robust -.0100).
  mv <- suppressWarnings(measurement_invariance(
    holzinger_swineford, items = f$verbal, group = "school",
    test = "mean.var.adjusted", levels = c("configural", "metric", "scalar")))
  mm <- suppressWarnings(measurement_invariance(
    holzinger_swineford, items = f$verbal, group = "school",
    estimator = "MLM", levels = c("configural", "metric", "scalar")))
  expect_identical(attr(mv, "fit_indices"), "robust")
  for (col in c("cfi", "rmsea", "delta_cfi", "delta_rmsea")) {
    expect_equal(mv[[col]], mm[[col]], tolerance = 1e-8, info = col)
  }
  expect_equal(mv$delta_cfi[2L], -0.01000519, tolerance = 1e-6)
  expect_match(attr(mv, "fit_columns"), "computed by DMAR", fixed = TRUE)
  # A bootstrap test is listed in lavaan's test option as a corrected
  # statistic but gives no scaled one: the chi_square column and the
  # difference test are the standard ones, and the attributes say so
  # (through 1.0.1 model_test said "scaled", fit_indices "robust", and test
  # "scaled chi square difference test"). The bootstrap p-values are not
  # read, so five draws suffice.
  set.seed(113)
  b <- suppressWarnings(measurement_invariance(
    holzinger_swineford, f, group = "school", test = "bollen.stine",
    bootstrap = 5, levels = c("configural", "metric")))
  fb <- lapply(attr(b, "fits"), lavaan::fitMeasures, c("chisq", "cfi"))
  expect_identical(attr(b, "model_test"), "standard")
  expect_identical(attr(b, "fit_indices"), "standard")
  expect_identical(attr(b, "test"), "standard chi square difference test")
  expect_equal(b$chi_square, vapply(fb, function(x) x[["chisq"]], numeric(1L),
                                    USE.NAMES = FALSE), tolerance = 1e-12)
  expect_equal(b$cfi, vapply(fb, function(x) x[["cfi"]], numeric(1L),
                             USE.NAMES = FALSE), tolerance = 1e-12)
  expect_equal(b$delta_chi_square[2L], diff(b$chi_square), tolerance = 1e-10)
})


test_that("an ordered ladder without robust indices at a rung reports the scaled ones at every rung", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # an ordered ladder under WLSMV
  # The ordered example of the help page. Under the delta parameterization
  # lavaan computes no robust cfi and rmsea for the thresholds rung, and
  # through 1.0.1 the table reported the robust indices of the configural
  # rung (.980 and .084) beside NA at the thresholds rung, with no delta
  # CFI. Every rung now reports the scaled indices, the version lavaan
  # computes at both, and the line under the table says so.
  f <- list(verbal = c("t6_paragraph_comprehension", "t7_sentence",
                       "t9_word_meaning"),
            deduction = c("t20_deduction", "t22_problem_reasoning",
                          "t23_series_completion"))
  ordered_data <- hs_discretized(holzinger_swineford,
                                 unlist(f, use.names = FALSE))
  mi <- suppressMessages(measurement_invariance(
    ordered_data, f, group = "school", ordered = TRUE,
    levels = c("configural", "thresholds")))
  fm <- lapply(attr(mi, "fits"), lavaan::fitMeasures,
               c("cfi.robust", "rmsea.robust", "cfi.scaled", "rmsea.scaled"))
  pick <- function(m) vapply(fm, function(x) x[[m]], numeric(1L),
                             USE.NAMES = FALSE)
  expect_true(is.finite(pick("cfi.robust")[1L]))
  expect_true(is.na(pick("cfi.robust")[2L]))
  expect_identical(attr(mi, "fit_indices"), "scaled")
  expect_equal(mi$cfi, pick("cfi.scaled"), tolerance = 1e-12)
  expect_equal(mi$rmsea, pick("rmsea.scaled"), tolerance = 1e-12)
  expect_equal(mi$cfi, c(0.9985285478, 0.9962103896), tolerance = 1e-8)
  expect_equal(mi$rmsea, c(0.02900155231, 0.03969118812), tolerance = 1e-8)
  expect_equal(mi$delta_cfi[2L], -0.002318158225, tolerance = 1e-6)
  expect_equal(pick("cfi.robust")[1L], 0.9795737476, tolerance = 1e-8)
  expect_match(attr(mi, "fit_columns"), "hold its scaled versions",
               fixed = TRUE)
})


test_that("a dichotomous ladder drops the rungs that cannot be tested", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  data(holzinger_swineford)
  factors <- hs_factors()
  items <- unlist(factors, use.names = FALSE)
  binary <- hs_dichotomized(holzinger_swineford, items)
  expect_true(all(.measurement_invariance_n_categories(binary, items) == 2L))
  # The count is of realized categories, not declared factor levels, so an
  # ordered factor carrying a level the data never realize is still binary.
  padded <- binary
  for (item in items) {
    padded[[item]] <- factor(binary[[item]], levels = c(0, 1, 2),
                             ordered = TRUE)
  }
  expect_true(all(.measurement_invariance_n_categories(padded, items) == 2L))

  run <- hs_conditions(measurement_invariance(binary, factors,
                                              group = "school",
                                              ordered = TRUE))
  res <- run$value

  # A dichotomous item has one threshold, which is not separately identified
  # from the intercept of its underlying response, and no free residual
  # variance (Millsap & Yun-Tein, 2004; Wu & Estabrook, 2016). The thresholds
  # and strict rungs are therefore not on the ladder, and the user is told so.
  expect_identical(res$level, c("configural", "metric", "scalar"))
  expect_true(any(grepl("dichotomous", run$messages)))
  expect_true(any(grepl("not separately identified", run$messages)))

  # This is the reason the thresholds rung is dropped, asserted directly:
  # fitting it by hand costs nothing, so it can carry no test.
  hand <- function(equal) {
    args <- list(model = hs_syntax(), data = binary, group = "school",
                 ordered = items, estimator = "WLSMV", meanstructure = TRUE)
    if (length(equal)) args$group.equal <- equal
    suppressWarnings(do.call(lavaan::cfa, args))
  }
  expect_equal(unname(lavaan::fitMeasures(hand("thresholds"), "df.scaled")),
               unname(lavaan::fitMeasures(hand(character(0)), "df.scaled")))

  # Every rung that is fitted costs degrees of freedom, so every reported
  # difference test is a test of something. Before the ladder was corrected
  # the thresholds row carried delta_df of 0 and a small negative
  # delta_chi_square.
  expect_true(all(diff(res$df) > 0))
  expect_equal(res$delta_df[-1L], diff(res$df))
  expect_true(all(res$delta_chi_square[-1L] > 0))
  expect_true(all(is.finite(res$p_value[-1L])))
  expect_identical(attr(res, "estimator"), "WLSMV")
  expect_true(attr(res, "ordered"))

  # lavaan's own note about two models sharing their degrees of freedom is
  # muffled in favor of the package's explanation, and nothing else is.
  expect_false(any(grepl("same degrees of freedom", run$warnings)))

  # Asking for a rung that is not identified says why rather than fitting it.
  expect_error(
    suppressMessages(measurement_invariance(
      binary, factors, group = "school", ordered = TRUE,
      levels = c("configural", "metric", "scalar", "strict"))),
    "not on the dichotomous ordered ladder")
})


test_that("dichotomous items beside continuous ones keep the strict rung", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  data(holzinger_swineford)
  factors <- hs_factors()
  items <- unlist(factors, use.names = FALSE)
  # Three dichotomous indicators declared ordered; the other six stay
  # continuous. Their residual variances are free parameters, so the strict
  # rung still restricts something even though the thresholds rung does not.
  partly <- hs_dichotomized(holzinger_swineford, items[1:3])
  run <- hs_conditions(measurement_invariance(partly, factors,
                                              group = "school",
                                              ordered = items[1:3]))
  res <- run$value
  expect_identical(res$level, c("configural", "metric", "scalar", "strict"))
  expect_true(any(grepl("strict rung is kept", run$messages)))
  expect_true(all(diff(res$df) > 0))
  expect_true(all(res$delta_chi_square[-1L] > 0))
})


test_that("a partly dichotomous ladder keeps every rung and names the items", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  data(holzinger_swineford)
  factors <- hs_factors()
  items <- unlist(factors, use.names = FALSE)
  # Three dichotomous items among six with four categories: the polytomous
  # items still carry thresholds and residual variances to test, so the full
  # ordered ladder stays in force.
  mixed <- hs_dichotomized(hs_discretized(holzinger_swineford, items[4:9]),
                           items[1:3])
  run <- hs_conditions(measurement_invariance(mixed, factors,
                                              group = "school",
                                              ordered = TRUE))
  res <- run$value
  expect_identical(res$level, c("configural", "thresholds", "metric",
                                "scalar", "strict"))
  expect_true(any(grepl("3 of 9", run$messages, fixed = TRUE)))
  expect_true(any(grepl("t1_visual_perception", run$messages, fixed = TRUE)))
  expect_true(all(diff(res$df) > 0))
  expect_true(all(res$delta_chi_square[-1L] > 0))
})


test_that("a rung costing no degrees of freedom reports no difference test", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford)
  items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board")
  # Freeing both free loadings across groups leaves the metric rung with
  # nothing to constrain, which is the continuous analog of the vacuous
  # thresholds rung: the two models are one model, so their chi square
  # difference is numerical noise on zero degrees of freedom.
  run <- hs_conditions(measurement_invariance(
    holzinger_swineford, items = items, group = "school",
    levels = c("configural", "metric"),
    group_partial = c("f =~ t2_cubes", "f =~ t3_paper_form_board")))
  res <- run$value
  expect_equal(res$delta_df[2L], 0)
  expect_true(is.na(res$delta_chi_square[2L]))
  expect_true(is.na(res$p_value[2L]))
  expect_true(any(grepl("cost no degrees of freedom", run$messages)))
  expect_false(any(grepl("same degrees of freedom", run$warnings)))
})


test_that("ordered data promote DWLS and ULS to a scaled difference test", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  data(holzinger_swineford)
  factors <- hs_factors()
  items <- unlist(factors, use.names = FALSE)
  ordered_data <- hs_discretized(holzinger_swineford, items)
  two <- c("configural", "thresholds")

  fit_with <- function(estimator) {
    hs_conditions(measurement_invariance(ordered_data, factors,
                                         group = "school", ordered = TRUE,
                                         estimator = estimator,
                                         levels = two))
  }
  estimates <- function(res) {
    lavaan::parameterEstimates(attr(res, "fits")$thresholds)$est
  }

  # The promotion changes the test, not the fit: in lavaan "WLSMV" is "DWLS"
  # and "ULSMV" is "ULS" with the mean and variance adjusted statistic, so
  # every parameter estimate is identical to the unpromoted run.
  for (pair in list(c("DWLS", "WLSMV"), c("ULS", "ULSMV"))) {
    run <- fit_with(pair[1L])
    target <- fit_with(pair[2L])
    expect_identical(attr(run$value, "estimator"), pair[2L])
    expect_match(attr(run$value, "test"), "^scaled chi square difference test")
    # lavaan computes no robust indices for the thresholds rung under the
    # delta parameterization, so both rungs report the scaled ones.
    expect_identical(attr(run$value, "fit_indices"), "scaled")
    expect_true(any(grepl(paste0("switching to estimator = \"", pair[2L], "\""),
                          run$messages, fixed = TRUE)))
    expect_equal(estimates(run$value), estimates(target$value),
                 tolerance = 1e-12)
    expect_equal(run$value$chi_square, target$value$chi_square,
                 tolerance = 1e-12)
    expect_equal(run$value$p_value, target$value$p_value, tolerance = 1e-12)
  }

  # The full weight matrix estimator is left alone: its statistic is
  # asymptotically chi square already (Browne, 1984; Muthen, 1984), and a
  # substitution would change the weight matrix rather than only the test.
  wls <- fit_with("WLS")
  expect_identical(attr(wls$value, "estimator"), "WLS")
  expect_identical(attr(wls$value, "test"),
                   "standard chi square difference test")
  expect_false(any(grepl("switching to estimator", wls$messages)))
})


test_that("p_chi_square is displayed as a p-value, like p_value beside it", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two more ladders fit only to format their p-values
  data(holzinger_swineford)
  items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board")
  res <- measurement_invariance(holzinger_swineford, items = items,
                                group = "school")
  disp <- format(res)
  # Fixed decimals, not significant figures, in both columns. The configural
  # rung is saturated, so its exact-fit p-value is NA and is not formatted.
  expect_true(all(grepl("^[0-9]+\\.[0-9]{4}$", disp$p_chi_square[-1L])))
  expect_true(all(grepl("^[0-9]+\\.[0-9]{4}$", disp$p_value[-1L])))
  expect_identical(disp$p_chi_square[2L],
                   formatC(res$p_chi_square[2L], format = "f", digits = 4L))
  expect_identical(disp$p_value[2L],
                   formatC(res$p_value[2L], format = "f", digits = 4L))
  # Stored values keep full precision whatever the display does.
  expect_false(isTRUE(all.equal(res$p_chi_square[2L],
                                round(res$p_chi_square[2L], 4L),
                                tolerance = 0)))

  # A p-value below the displayable magnitude gets the floor rather than
  # scientific notation or a bare zero.
  nine <- unlist(hs_factors(), use.names = FALSE)
  poor <- measurement_invariance(holzinger_swineford, items = nine,
                                 group = "school")
  expect_identical(unique(format(poor)$p_chi_square), "< 0.0001")
})


test_that("levels must be a leading subset of the ladder in force", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford)
  items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board")
  res <- measurement_invariance(holzinger_swineford, items = items,
                                group = "school",
                                levels = c("configural", "metric"))
  expect_identical(res$level, c("configural", "metric"))
  # A single rung compares nothing, so no difference test is named.
  one <- measurement_invariance(holzinger_swineford, items = items,
                                group = "school", levels = "configural")
  expect_identical(one$level, "configural")
  expect_true(is.na(attr(one, "test")))
  expect_error(measurement_invariance(holzinger_swineford, items = items,
                                      group = "school",
                                      levels = c("configural", "scalar")),
               "leading subset")
  # "thresholds" is not a rung of the continuous ladder.
  expect_error(measurement_invariance(holzinger_swineford, items = items,
                                      group = "school",
                                      levels = c("configural",
                                                 "thresholds")),
               "not on the continuous ladder")
})


test_that("the ordered ladder accepts its own leading subsets only", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  data(holzinger_swineford)
  items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board")
  ordered_data <- hs_discretized(holzinger_swineford, items)
  res <- suppressMessages(
    measurement_invariance(ordered_data, items = items, group = "school",
                           ordered = TRUE,
                           levels = c("configural", "thresholds")))
  expect_identical(res$level, c("configural", "thresholds"))
  # Skipping the thresholds rung is exactly the methodological error the
  # ordered ladder exists to prevent.
  expect_error(
    suppressMessages(
      measurement_invariance(ordered_data, items = items, group = "school",
                             ordered = TRUE,
                             levels = c("configural", "metric"))),
    "leading subset")
})


test_that("measurement_invariance() validates its inputs and fails cleanly", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford)
  items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board")

  expect_error(measurement_invariance(holzinger_swineford, group = "school"),
               "exactly one")
  expect_error(measurement_invariance(holzinger_swineford, hs_factors(),
                                      group = "school", items = items),
               "exactly one")
  expect_error(measurement_invariance(as.matrix(holzinger_swineford[items]),
                                      items = items, group = "school"),
               "data.frame")
  expect_error(measurement_invariance(holzinger_swineford,
                                      items = items[1:2], group = "school"),
               "three or more")
  expect_error(measurement_invariance(holzinger_swineford, items = items,
                                      group = "nope"),
               "one column")
  one_group <- holzinger_swineford[holzinger_swineford$school ==
                                     holzinger_swineford$school[1L], ]
  expect_error(measurement_invariance(one_group, items = items,
                                      group = "school"),
               "at least two groups")
  expect_error(measurement_invariance(holzinger_swineford,
                                      model = "visual + verbal",
                                      group = "school"),
               "measurement relation")
  expect_error(measurement_invariance(holzinger_swineford,
                                      model = list(items),
                                      group = "school"),
               "must be named")
  expect_error(measurement_invariance(holzinger_swineford,
                                      model = list(visual = items[1:2]),
                                      group = "school"),
               "three or more items")
  expect_error(measurement_invariance(
    holzinger_swineford,
    model = list(visual = items, verbal = c("t7_sentence", items[1L],
                                            "t9_word_meaning")),
    group = "school"), "exactly one factor")
  expect_error(measurement_invariance(
    holzinger_swineford,
    model = list(visual = c(items[1:2], "no_such_item")),
    group = "school"), "not found in 'data'")
  expect_error(measurement_invariance(
    holzinger_swineford,
    model = "visual =~ t1_visual_perception + t2_cubes + no_such_item",
    group = "school"), "not found in 'data'")
  expect_error(measurement_invariance(holzinger_swineford, items = items,
                                      group = "school",
                                      estimator = "OLS"),
               "'estimator' must be one of")
  expect_error(measurement_invariance(holzinger_swineford, items = items,
                                      group = "school",
                                      missing = "impute"),
               "'missing' must be one of")
  expect_error(measurement_invariance(holzinger_swineford, items = items,
                                      group = "school", ordered = TRUE,
                                      missing = "fiml"),
               "Full information maximum likelihood is not available")
  expect_error(measurement_invariance(holzinger_swineford, items = items,
                                      group = "school",
                                      ordered = c("t1_visual_perception",
                                                  "not_an_item")),
               "not indicators of the model")
  expect_error(measurement_invariance(holzinger_swineford, items = items,
                                      group = "school", ordered = 1),
               "must be TRUE")
  expect_error(measurement_invariance(holzinger_swineford, items = items,
                                      group = "school",
                                      group_partial = 42),
               "character vector of parameter labels")
})

test_that("a misnamed item is refused with the class of the factor analysis items", {
  skip_if_not_installed("lavaan")
  items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board")
  hs <- holzinger_swineford
  # The reasons and their order are those of fa_exploratory()'s items
  # (.fa_items_check()); the messages are the ones the function gave
  # without a class.
  refused <- function(expr, reason, named, message) {
    e <- tryCatch(expr, dmar_items_invalid = function(e) e)
    expect_s3_class(e, c("dmar_items_invalid", "dmar_refusal", "error"))
    expect_identical(e$code, "items_invalid")
    expect_identical(e$reason, reason)
    expect_identical(e$items, named)
    expect_identical(conditionMessage(e), message)
  }
  three <- "'items' must name three or more columns of 'data'."
  # The items of the one-factor form.
  refused(measurement_invariance(hs, items = c(items[1:2], "no_such_item"),
                                 group = "school"),
          "unknown", "no_such_item", three)
  refused(measurement_invariance(hs, items = items[1:2], group = "school"),
          "count", character(0), three)
  refused(measurement_invariance(hs, items = c(items, NA), group = "school"),
          "form", character(0), three)
  refused(measurement_invariance(hs, items = c(items, items[1L]),
                                 group = "school"),
          "duplicate", items[1L], "'items' must not repeat an item.")
  # The factors of the list form, as fa_confirmatory() refuses its own.
  refused(measurement_invariance(
    hs, model = list(visual = c(items[1:2], "no_such_item")),
    group = "school"),
    "unknown", "no_such_item", "Items not found in 'data': no_such_item.")
  refused(measurement_invariance(
    hs, model = list(visual = c(items[1:2], NA)), group = "school"),
    "form", character(0), "Items not found in 'data': NA.")
  # The indicators of model syntax.
  refused(measurement_invariance(
    hs, model = "visual =~ t1_visual_perception + t2_cubes + no_such_item",
    group = "school"),
    "unknown", "no_such_item", "Indicators not found in 'data': no_such_item.")
  # A factor with too few items, or an item named twice among the factors
  # of the list form (on two factors or twice on one), is a question of
  # the model rather than of the names, and is not this class; only
  # `items`, from which DMAR writes the model, is refused for a repeat.
  e <- tryCatch(measurement_invariance(hs, model = list(visual = items[1:2]),
                                       group = "school"),
                error = function(e) e)
  expect_false(inherits(e, "dmar_items_invalid"))
  e <- tryCatch(measurement_invariance(
    hs, model = list(visual = items, verbal = c("t7_sentence", items[1L],
                                                "t9_word_meaning")),
    group = "school"), error = function(e) e)
  expect_s3_class(e, "error")
  expect_false(inherits(e, "dmar_items_invalid"))
  e <- tryCatch(measurement_invariance(
    hs, model = list(visual = c(items, items[1L])), group = "school"),
    error = function(e) e)
  expect_s3_class(e, "error")
  expect_false(inherits(e, "dmar_items_invalid"))
})


test_that("Browne's residual statistics do not make a fit a robust one", {
  # lavaan 0.7-3 lists Browne's normal theory residual statistic beside the
  # standard one for an ML fit, and the ADF version ahead of the
  # Yuan-Bentler statistic for a weighted MLR fit. Neither is a correction
  # of the standard statistic, so an ML ladder stays a standard one (the
  # attribute test above checks the ladder itself).
  expect_identical(.lavaan_corrected_tests("standard"), character(0))
  expect_identical(
    .lavaan_corrected_tests(c("standard", "browne.residual.nt.model")),
    character(0))
  expect_identical(.lavaan_corrected_tests("none"), character(0))
  expect_identical(
    .lavaan_corrected_tests(c("standard", "yuan.bentler.mplus")),
    "yuan.bentler.mplus")
  expect_identical(
    .lavaan_corrected_tests(c("standard", "browne.residual.adf",
                              "yuan.bentler.mplus")),
    "yuan.bentler.mplus")
  expect_identical(
    .lavaan_corrected_tests(c("standard", "scaled.shifted")),
    "scaled.shifted")
})


test_that("a saturated rung reports a CFI of 1 under FIML as well", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a FIML ladder of three rungs
  # lavaan's FIML corrected CFI of the saturated configural rung is NA on
  # these data, which left the metric rung's delta_cfi NA.
  d <- holzinger_swineford
  set.seed(1)
  d$t2_cubes[sample(nrow(d), 60)] <- NA
  m <- suppressWarnings(measurement_invariance(
    d, items = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
    group = "school", missing = "fiml",
    levels = c("configural", "metric", "scalar")))
  # lavaan warns that its corrected CFI of the saturated rung is NA.
  fr <- suppressWarnings(lapply(attr(m, "fits"), lavaan::fitMeasures,
                                c("df", "cfi.robust", "rmsea.robust")))
  expect_identical(fr$configural[["df"]], 0)
  expect_true(is.na(fr$configural[["cfi.robust"]]))
  expect_identical(attr(m, "fit_indices"), "robust")
  expect_identical(m$cfi[1L], 1)
  expect_identical(m$rmsea[1L], 0)
  expect_equal(m$cfi[2:3], vapply(fr[2:3], function(x) x[["cfi.robust"]],
                                  numeric(1L), USE.NAMES = FALSE),
               tolerance = 1e-12)
  expect_equal(m$delta_cfi[2L], m$cfi[2L] - 1, tolerance = 1e-12)
  expect_false(is.na(m$delta_cfi[2L]))
})


test_that("an explicit FIML request lavaan would not honor stops, and a listwise ladder is not labeled FIML", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # four ladders of two rungs
  # lavaan cannot combine FIML with the Satorra-Bentler or the mean and
  # variance adjusted statistic, or with se = "robust.sem", and sets
  # missing to "listwise" with a warning. Through 1.0.1 an explicit
  # missing = "fiml" (or "ml") with one of them returned that listwise
  # ladder; it now stops before any fit, whether or not a value is
  # missing, as fa_confirmatory() does, and any other route by which
  # lavaan replaces FIML is caught on the fits.
  f <- list(verbal = c("t6_paragraph_comprehension", "t7_sentence",
                       "t9_word_meaning"),
            deduction = c("t20_deduction", "t22_problem_reasoning",
                          "t23_series_completion"))
  d <- holzinger_swineford
  set.seed(1)
  d$t7_sentence[sample(301, 40)] <- NA
  d$t20_deduction[sample(301, 30)] <- NA
  lv <- c("configural", "metric")
  for (m in c("ml", "fiml")) {
    for (test in c("mean.var.adjusted", "satorra.bentler")) {
      for (x in list(d, holzinger_swineford)) {
        expect_error(measurement_invariance(x, f, group = "school",
                                            missing = m, test = test,
                                            levels = lv),
                     paste0("missing = \"", m, "\" asks for full ",
                            "information maximum likelihood, which lavaan ",
                            "does not compute with test = \"", test, "\""),
                     fixed = TRUE)
      }
    }
  }
  expect_error(measurement_invariance(d, f, group = "school",
                                      missing = "fiml",
                                      test = "satorra.bentler", levels = lv),
               paste0("Use missing = \"listwise\" to fit the complete cases ",
                      "with it, or estimator = \"MLR\" in place of the test"),
               fixed = TRUE)
  expect_error(measurement_invariance(d, f, group = "school",
                                      missing = "fiml", estimator = "MLR",
                                      test = "satorra.bentler", levels = lv),
               paste0("or leave out 'test', since estimator = \"MLR\" ",
                      "computes its own scaled statistic under FIML."),
               fixed = TRUE)
  expect_error(measurement_invariance(d, f, group = "school",
                                      missing = "fiml", se = "robust.sem",
                                      levels = lv),
               "or se = \"robust.huber.white\" in place of se = \"robust.sem\"",
               fixed = TRUE)
  # Caught on the fits: lavaan's own spelling of the test, and GLS, which
  # lavaan fits in two stages under missing = "ml" without a warning.
  expect_error(suppressWarnings(measurement_invariance(
    d, f, group = "school", missing = "fiml", test = "SB",
    levels = "configural")),
    "but lavaan fitted the ladder with missing = \"listwise\" instead",
    fixed = TRUE)
  expect_error(suppressWarnings(measurement_invariance(
    d, f, group = "school", missing = "fiml", estimator = "GLS",
    levels = "configural")),
    "but lavaan fitted the ladder with missing = \"two.stage\" instead",
    fixed = TRUE)
  # The tests with listwise deletion asked for by name: whether a ladder
  # used FIML is read from its fits, so the line under the table names
  # the robust indices DMAR computed, and no FIML version.
  mi <- suppressWarnings(measurement_invariance(
    d, f, group = "school", missing = "listwise", test = "mean.var.adjusted",
    levels = lv))
  expect_identical(lavaan::lavInspect(attr(mi, "fits")[[1L]],
                                      "options")$missing, "listwise")
  expect_false(any(grepl("FIML", attr(mi, "fit_columns"), fixed = TRUE)))
  expect_match(attr(mi, "fit_columns"), "computed by DMAR", fixed = TRUE)
  ml <- suppressWarnings(measurement_invariance(
    d, f, group = "school", estimator = "MLM", missing = "listwise",
    levels = lv))
  expect_equal(mi$cfi, ml$cfi, tolerance = 1e-8)
  sb <- suppressWarnings(measurement_invariance(
    d, f, group = "school", missing = "listwise", test = "satorra.bentler",
    levels = lv))
  expect_null(attr(sb, "fit_columns"))
  expect_equal(sb$cfi, ml$cfi, tolerance = 1e-8)
})


test_that("a ULS ladder of continuous items reports no index that depends on the units", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # six ULS ladders
  # The ULS statistic is in the units of the squared covariances. lavaan's
  # robust RMSEA of a ULSM ladder puts it into the correction of
  # Brosseau-Liard et al. (2012) for the ML fit function: through 1.0.1 the
  # rmsea column held it, 1.57 and 4.92 here, and a hundredth of each with
  # the scores divided by 10. At the scalar and strict rungs, which
  # restrict the means, the fit function adds the squared mean residuals to
  # the squared covariance residuals, in other units, so the CFI moves too:
  # a scalar delta CFI of -.00005 for the scores and -.0046 for the scores
  # divided by 10. Those columns are NA; the robust CFI of the rungs with
  # free means is also withheld under the ML-derivation rule. Eligible
  # scaled indices, the chi_square column and difference test remain.
  f <- list(verbal = c("t6_paragraph_comprehension", "t7_sentence",
                       "t9_word_meaning"),
            deduction = c("t20_deduction", "t22_problem_reasoning",
                          "t23_series_completion"))
  h10 <- holzinger_swineford
  h10[unlist(f)] <- h10[unlist(f)] / 10
  ladder <- function(d, ...) suppressWarnings(measurement_invariance(
    d, f, group = "school", ...))
  mi <- ladder(holzinger_swineford, estimator = "ULSM")
  mi10 <- ladder(h10, estimator = "ULSM")
  fm <- lapply(attr(mi, "fits"), lavaan::fitMeasures,
               c("rmsea.robust", "cfi.robust", "cfi.scaled", "rmsea.scaled"))
  fm10 <- lapply(attr(mi10, "fits"), lavaan::fitMeasures, "cfi.robust")
  expect_gt(fm$configural[["rmsea.robust"]], 1.5)
  expect_equal(fm$scalar[["cfi.robust"]] - fm$metric[["cfi.robust"]],
               -0.0000487, tolerance = 0.05)
  expect_equal(fm10$scalar[["cfi.robust"]] - fm10$metric[["cfi.robust"]],
               -0.0046, tolerance = 0.05)
  expect_identical(attr(mi, "fit_indices"), "scaled")
  expect_true(all(is.finite(mi$rmsea[1:2])))
  expect_true(all(is.na(mi$rmsea[3:4])))
  expect_true(all(is.na(mi$delta_rmsea[3:4])))
  expect_equal(mi$rmsea[1:2], vapply(fm[1:2], function(x) x[["rmsea.scaled"]],
                                    numeric(1L), USE.NAMES = FALSE), tolerance = 1e-12)
  expect_equal(mi$cfi[1:2], vapply(fm[1:2], function(x) x[["cfi.scaled"]],
                                   numeric(1L), USE.NAMES = FALSE),
               tolerance = 1e-12)
  expect_equal(mi$cfi[1:2], mi10$cfi[1:2], tolerance = 1e-6)
  expect_true(all(is.na(mi$cfi[3:4])))
  expect_true(all(is.na(mi$delta_cfi[3:4])))
  expect_false(anyNA(mi$chi_square))
  expect_false(anyNA(mi$p_value[2:4]))
  expect_match(attr(mi, "fit_columns"), "at the scalar and strict rungs",
               fixed = TRUE)
  expect_output(print(mi), "ML-derived robust indices are withheld", fixed = TRUE)
  expect_match(attr(mi, "fit_columns"), "ML cutoffs do not automatically apply", fixed = TRUE)
  # The standard indices of a ULS ladder (test = "standard") set the degrees
  # of freedom against the ULS statistic and are NA at every rung.
  ms <- ladder(holzinger_swineford, estimator = "ULS", test = "standard",
               levels = c("configural", "metric"))
  expect_identical(attr(ms, "fit_indices"), "standard")
  expect_true(all(is.na(c(ms$cfi, ms$rmsea))))
  expect_match(attr(ms, "fit_columns")[["fit_columns_uls"]],
               "(N - G) times the ULS fit function", fixed = TRUE)
  # ULS alone is tested by Browne's residual based statistic, which has no
  # units, so its cfi and rmsea are kept while the ULS estimates do not
  # change with the units, at the configural and metric rungs; but
  # lavaan's difference test refers the difference of the ULS statistics
  # to the chi square, and through 1.0.1 the metric rung had a
  # delta_chi_square of 91712 (p < .0001), and of 9.17 (p = .057) with
  # the scores divided by 10.
  mu <- ladder(holzinger_swineford, estimator = "ULS")
  mu10 <- ladder(h10, estimator = "ULS")
  lrt <- lavaan::lavTestLRT(attr(mu, "fits")[[1L]], attr(mu, "fits")[[2L]])
  expect_gt(lrt[2L, "Chisq diff"], 9e4)
  expect_equal(mu$cfi[1:2], mu10$cfi[1:2], tolerance = 1e-6)
  expect_equal(mu$rmsea[1:2], mu10$rmsea[1:2], tolerance = 1e-6)
  expect_false(anyNA(c(mu$cfi[1:2], mu$rmsea[1:2], mu$delta_cfi[2L])))
  expect_true(is.na(mu$delta_chi_square[2L]))
  expect_true(is.na(mu$p_value[2L]))
  expect_identical(mu$delta_df[2L], 4)
  expect_true(is.na(attr(mu, "test")))
  expect_match(attr(mu, "fit_columns")[["fit_columns_uls_test"]],
               "so there is none", fixed = TRUE)
  # At the scalar and strict rungs the means are restricted through the
  # loadings, which enter them through the latent means, so the ULS
  # estimates change with the units, and so does Browne's statistic
  # evaluated at them: through the first version of the ULS rule the
  # scalar CFI was .8874 for the scores and .8861 for the scores divided
  # by 10, and the scalar delta CFI -.0212 and -.0225. Those columns are
  # NA; the chi_square column is kept.
  fs <- lapply(list(attr(mu, "fits")$scalar, attr(mu10, "fits")$scalar),
               lavaan::fitMeasures, "cfi")
  expect_equal(unname(fs[[1L]]), 0.8874, tolerance = 1e-3)
  expect_equal(unname(fs[[2L]]), 0.8861, tolerance = 1e-3)
  expect_true(all(is.na(c(mu$cfi[3:4], mu$rmsea[3:4], mu$delta_cfi[3:4],
                          mu$delta_rmsea[3:4]))))
  expect_false(anyNA(mu$chi_square))
  expect_equal(mu$chi_square[3L],
               unname(lavaan::fitMeasures(attr(mu, "fits")$scalar, "chisq")))
  expect_match(attr(mu, "fit_columns")[["fit_columns_uls"]],
               paste("at the scalar and strict rungs, which restrict the means",
                     "through parameters the covariances share"),
               fixed = TRUE)
  expect_output(print(mu), "formed from Browne's residual based statistic",
                fixed = TRUE)
})


test_that("an all-saturated ULS ladder still excludes the robust candidate", {
  skip_if_not_installed("lavaan")
  items <- c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  x <- suppressWarnings(measurement_invariance(holzinger_swineford,
    items = items, group = "school", estimator = "ULSM", levels = "configural"))
  expect_identical(attr(x, "fit_indices"), "scaled")
  expect_identical(x$cfi, 1)
  expect_identical(x$rmsea, 0)
  expect_match(attr(x, "fit_columns"), "ML-derived robust indices are withheld", fixed = TRUE)
})


test_that("the default 'missing' fits incomplete items by FIML and complete data as they stand", {
  skip_if_not_installed("lavaan")
  # DMAR 1.0.x deleted incomplete cases listwise by default. Since 1.1.0
  # the default fits incomplete items by FIML, as fa_confirmatory() does,
  # and fits complete data as they stand, so that no complete data ladder
  # changes.
  m <- "verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning"
  lv <- c("configural", "metric")
  hs <- holzinger_swineford
  d <- hs
  d$t7_sentence[seq(5, 300, by = 10)] <- NA
  ladder <- function(x, levels = lv, ...) {
    measurement_invariance(x, model = m, group = "school", levels = levels,
                           ...)
  }
  columns <- function(x) {
    stats::setNames(lapply(names(x), function(n) x[[n]]), names(x))
  }
  lavaan_missing <- function(x) {
    lavaan::lavInspect(attr(x, "fits")[[1L]], "options")$missing
  }
  a <- ladder(hs)
  b <- ladder(hs, missing = "listwise")
  expect_identical(lavaan_missing(a), "listwise")
  expect_identical(attr(a, "fit_indices"), "standard")
  expect_identical(columns(a), columns(b))
  # Incomplete items: FIML with every case, without a warning.
  expect_no_warning(x <- ladder(d))
  expect_identical(lavaan_missing(x), "ml")
  expect_identical(lavaan::lavInspect(attr(x, "fits")[[1L]], "ntotal"),
                   nrow(d))
  expect_equal(columns(x), columns(ladder(d, missing = "fiml")))
  expect_no_warning(l <- ladder(d, missing = "listwise"))
  expect_identical(lavaan::lavInspect(attr(l, "fits")[[1L]], "ntotal"),
                   sum(!is.na(d$t7_sentence)))
  # What FIML cannot accompany is fitted as before with complete data and,
  # with values missing, stops and names the options.
  expect_identical(lavaan_missing(ladder(hs, estimator = "MLM")),
                   "listwise")
  expect_error(ladder(d, estimator = "MLM"),
               "not available with estimator = \"MLM\"", fixed = TRUE)
  expect_error(ladder(d, estimator = "MLM"),
               "or estimator = \"MLR\", whose robust standard errors",
               fixed = TRUE)
  expect_error(ladder(d, se = "robust.sem"),
               "se = \"robust.huber.white\" in place of", fixed = TRUE)
  expect_error(ladder(d, test = "satorra.bentler"),
               "not available with test = \"satorra.bentler\"", fixed = TRUE)
  # With MLR already the estimator, the remedy for the test is to leave it
  # out, since MLR's own scaled statistic is computed under FIML.
  expect_error(ladder(d, estimator = "MLR", test = "satorra.bentler"),
               paste0("Use missing = \"listwise\" to fit the complete ",
                      "cases, or leave out 'test', since estimator = ",
                      "\"MLR\" computes its own scaled statistic under ",
                      "FIML."), fixed = TRUE)
  expect_error(ladder(d, estimator = "MLM", test = "satorra.bentler"),
               "or estimator = \"MLR\" without 'test', whose robust",
               fixed = TRUE)
  o <- hs_discretized(hs, c("t6_paragraph_comprehension", "t7_sentence",
                            "t9_word_meaning"))
  quiet <- function(x) suppressWarnings(suppressMessages(x))
  expect_identical(lavaan_missing(quiet(ladder(o, ordered = TRUE,
                                              levels = "configural"))),
                   "listwise")
  o$t7_sentence[seq(5, 300, by = 10)] <- NA
  expect_error(quiet(ladder(o, ordered = TRUE, levels = "configural")),
               "or missing = \"pairwise\" for ordered items", fixed = TRUE)
  # A value missing only where the group is missing changes nothing,
  # since lavaan sets such a case aside under either treatment.
  g <- hs
  g$school <- as.character(g$school)
  g$school[1:3] <- NA
  g$t7_sentence[1:3] <- NA
  expect_identical(lavaan_missing(quiet(ladder(g))), "listwise")
  # Nor does a case with every model variable missing: listwise deletion
  # and FIML both set it aside, so complete data with a blank row keep the
  # complete data ladder, with no warning.
  blank <- hs
  blank[nrow(blank) + 1L, ] <- NA
  blank$school[nrow(blank)] <- "Pasteur"
  expect_no_warning(e <- ladder(blank))
  expect_identical(lavaan_missing(e), "listwise")
  expect_equal(columns(e), columns(a))
  # The default FIML ladder is lavaan's own FIML fit of each rung.
  fc <- lavaan::cfa(m, data = d, group = "school", missing = "ml")
  fm <- lavaan::cfa(m, data = d, group = "school", missing = "ml",
                    group.equal = "loadings")
  expect_equal(x$chi_square,
               unname(c(lavaan::fitMeasures(fc, "chisq"),
                        lavaan::fitMeasures(fm, "chisq"))),
               tolerance = 1e-8)
  expect_equal(x$cfi[2L], unname(lavaan::fitMeasures(fm, "cfi.robust")),
               tolerance = 1e-8)
  expect_equal(x$rmsea[2L], unname(lavaan::fitMeasures(fm, "rmsea.robust")),
               tolerance = 1e-8)
})


test_that("under FIML a group that never observes a variable, or a pair, is refused by name", {
  skip_if_not_installed("lavaan")
  # holzinger_swineford gives t25_paper_form_board_r and t26_flags in
  # Grant-White only. Under FIML, the default for these incomplete items,
  # lavaan stopped with "subscript out of bounds"; the call is now refused
  # before any fit, with what to do instead.
  hs <- holzinger_swineford
  items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
             "t25_paper_form_board_r")
  e <- tryCatch(measurement_invariance(hs, items = items, group = "school"),
                error = function(e) e)
  expect_s3_class(e, c("dmar_group_coverage", "dmar_refusal", "error"))
  expect_identical(e$code, "group_coverage")
  expect_equal(e$unobserved,
               data.frame(group = "Pasteur",
                          variable = "t25_paper_form_board_r",
                          n_group = 156L))
  expect_identical(nrow(e$pairs), 0L)
  msg <- conditionMessage(e)
  expect_match(msg, paste0("group Pasteur never observes ",
                           "t25_paper_form_board_r (none of its 156 cases ",
                           "has a value)"), fixed = TRUE)
  expect_match(msg, paste0("Leave t25_paper_form_board_r out of the model ",
                           "for this comparison, so that the ladder ",
                           "compares the groups on what every group ",
                           "observes."), fixed = TRUE)
  expect_match(msg, paste0("Listwise deletion (missing = \"listwise\") is ",
                           "no remedy: it keeps only the cases observed on ",
                           "every variable of the model, and group Pasteur ",
                           "has none"), fixed = TRUE)
  expect_no_match(msg, "out of 'data'", fixed = TRUE)
  # An explicit request is refused the same way.
  expect_error(measurement_invariance(hs, items = items, group = "school",
                                      missing = "fiml"),
               class = "dmar_group_coverage")
  # Listwise deletion indeed leaves no case of Pasteur, and lavaan stops.
  expect_error(measurement_invariance(hs, items = items, group = "school",
                                      missing = "listwise"),
               "Pasteur", fixed = TRUE)
  # The remedy the message names fits.
  four <- c(items[1:3], "t4_lozenges")
  ok <- measurement_invariance(hs, items = four, group = "school",
                               levels = "configural")
  expect_identical(as.character(ok$level), "configural")
  # A default FIML ladder of incomplete items fits without a warning.
  d <- hs
  d$t2_cubes[seq(4, 300, by = 10)] <- NA
  expect_no_warning(measurement_invariance(d, items = four, group = "school",
                                           levels = "configural"))
  skip_on_cran()  # the rest builds further data sets
  # A pair no case of a group observes together.
  p <- hs
  pas <- which(p$school == "Pasteur")
  p$t1_visual_perception[pas[1:80]] <- NA
  p$t2_cubes[pas[81:156]] <- NA
  e <- tryCatch(measurement_invariance(p, items = items[1:3],
                                       group = "school"),
                error = function(e) e)
  expect_s3_class(e, "dmar_group_coverage")
  expect_equal(e$pairs,
               data.frame(group = "Pasteur",
                          variable_1 = "t1_visual_perception",
                          variable_2 = "t2_cubes"))
  expect_match(conditionMessage(e),
               paste0("no case of group Pasteur observes ",
                      "t1_visual_perception and t2_cubes together. Leave ",
                      "t1_visual_perception or t2_cubes out of the model"),
               fixed = TRUE)
  # An auxiliary variable a group never observes is to be left out of
  # 'aux', and listwise deletion is not discussed, since it concerns the
  # model's variables.
  e <- tryCatch(measurement_invariance(hs, items = items[1:3],
                                       group = "school", missing = "fiml",
                                       aux = "t26_flags"),
                error = function(e) e)
  expect_s3_class(e, "dmar_group_coverage")
  expect_match(conditionMessage(e), "Leave t26_flags out of 'aux'.",
               fixed = TRUE)
  expect_no_match(conditionMessage(e), "Listwise", fixed = TRUE)
  # With more groups, leaving out the groups at issue is offered when at
  # least two would remain.
  g <- hs
  g$cell <- interaction(g$school, g$grade, sep = "_")
  g$t2_cubes[g$cell == "Pasteur_7"] <- NA
  e <- tryCatch(measurement_invariance(g, items = items[1:3],
                                       group = "cell"),
                error = function(e) e)
  expect_s3_class(e, "dmar_group_coverage")
  expect_match(conditionMessage(e),
               paste0("; or leave group Pasteur_7 out of 'data', so that ",
                      "the ladder compares the other groups on every ",
                      "variable."), fixed = TRUE)
})
