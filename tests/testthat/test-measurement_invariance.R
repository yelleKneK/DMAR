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
                                     group = "school")
  expect_equal(sum(lavaan::lavInspect(attr(listwise, "fits")$metric,
                                      "nobs")),
               nrow(incomplete) - 20L)
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
  # nothing to constrain, which is the continuous analogue of the vacuous
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
    expect_identical(attr(run$value, "fit_indices"), "robust")
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
