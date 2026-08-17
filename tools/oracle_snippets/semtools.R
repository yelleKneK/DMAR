## from tests/testthat/test-measurement_invariance.R
local({
  # The ordered ladder against semTools::measEq.syntax. semTools builds the
  # Wu and Estabrook (2016) identification explicitly and with a different
  # factor identification (std.lv rather than lavaan's marker variable).
  # The models are equivalent reparameterizations, so the test statistics
  # and degrees of freedom must agree; agreement is to estimation-
  # convergence precision rather than to machine precision.
  hs_factors <- list(
    visual = c("t1_visual_perception", "t2_cubes", "t3_paper_form_board"),
    verbal = c("t6_paragraph_comprehension", "t7_sentence",
               "t9_word_meaning"),
    speed  = c("t20_deduction", "t22_problem_reasoning",
               "t23_series_completion"))
  hs_syntax <- paste(vapply(names(hs_factors), function(f)
    paste(f, "=~", paste(hs_factors[[f]], collapse = " + ")), character(1L)),
    collapse = "\n")

  data(holzinger_swineford, package = "DMAR")
  items <- unlist(hs_factors, use.names = FALSE)
  ordered_data <- holzinger_swineford
  for (item in items) {
    ordered_data[[item]] <- as.integer(cut(
      ordered_data[[item]],
      breaks = stats::quantile(ordered_data[[item]], c(0, .25, .5, .75, 1),
                               na.rm = TRUE),
      include.lowest = TRUE))
  }

  res <- suppressMessages(
    DMAR::measurement_invariance(ordered_data, hs_factors, group = "school",
                                 ordered = TRUE))

  equal <- list(configural = character(0),
                thresholds = "thresholds",
                metric     = c("thresholds", "loadings"),
                scalar     = c("thresholds", "loadings", "intercepts"),
                strict     = c("thresholds", "loadings", "intercepts",
                               "residuals"))
  external <- lapply(equal, function(eq) {
    syntax <- semTools::measEq.syntax(
      configural.model = hs_syntax, data = ordered_data, ordered = items,
      parameterization = "theta", ID.fac = "std.lv",
      ID.cat = "Wu.Estabrook.2016", group = "school", group.equal = eq,
      return.fit = FALSE)
    lavaan::cfa(as.character(syntax), data = ordered_data, group = "school",
                ordered = items, parameterization = "theta",
                estimator = "WLSMV")
  })

  # as.numeric() strips the lavaan.vector class so all.equal compares the
  # values rather than the attributes.
  for (i in seq_along(external)) {
    stopifnot(isTRUE(all.equal(
      res$df[i],
      as.numeric(lavaan::fitMeasures(external[[i]], "df.scaled")))))
    stopifnot(isTRUE(all.equal(
      res$chi_square[i],
      as.numeric(lavaan::fitMeasures(external[[i]], "chisq.scaled")),
      tolerance = 1e-4)))
  }
  for (i in 2:length(external)) {
    lrt <- lavaan::lavTestLRT(external[[i - 1L]], external[[i]])
    stopifnot(isTRUE(all.equal(res$delta_chi_square[i],
                               lrt[2L, "Chisq diff"], tolerance = 1e-4)))
    stopifnot(isTRUE(all.equal(res$delta_df[i], lrt[2L, "Df diff"])))
    stopifnot(isTRUE(all.equal(res$p_value[i], lrt[2L, "Pr(>Chisq)"],
                               tolerance = 1e-4)))
  }
})
