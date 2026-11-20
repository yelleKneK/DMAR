# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

fc_items <- function() c("t1_visual_perception", "t2_cubes", "t4_lozenges",
                         "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")

test_that("factor_count() sets fit, difference tests, criteria, eigenvalues, and MAP beside each count", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[fc_items()]
  fc <- factor_count(d, B = 30, seed = 113)
  expect_s3_class(fc, "dmar_factor_count")
  expect_s3_class(fc, "dmar_tbl")
  # Six items: the two-factor model keeps four degrees of freedom and the
  # three-factor model none, so the default stops at two.
  expect_identical(fc$factors, 1:2)
  expect_identical(fc$solution, c("proper", "proper"))
  expect_identical(fc$term, c("factors_1", "factors_2"))
  # The fit columns are fa_exploratory()'s fit measures at each count.
  fits <- lapply(1:2, function(k) fa_exploratory(d, k, rotation = "none"))
  for (k in 1:2) {
    fm <- attr(fits[[k]], "fit_measures")
    value <- function(term) fm$value[fm$term == term]
    expect_equal(fc$chi_square[k], value("chi_square"))
    expect_equal(fc$df[k], value("df"))
    expect_equal(fc$p_value[k], value("p_chi_square"))
    expect_equal(fc$rmsea[k], value("rmsea"))
    expect_equal(fc$rmsea_ci_upper[k], value("rmsea_ci_upper"))
    expect_equal(fc$cfi[k], value("cfi"))
    expect_equal(fc$tli[k], value("tli"))
    expect_equal(fc$srmr[k], value("srmr"))
    ic <- sem_information_criteria(attr(fits[[k]], "fit"))
    expect_equal(fc$aic[k], ic$value[ic$term == "AIC"])
    expect_equal(fc$bic[k], ic$value[ic$term == "BIC"])
    expect_equal(fc$sabic[k], ic$value[ic$term == "SABIC"])
    expect_equal(fc$hbic[k], ic$value[ic$term == "HBIC"])
  }
  expect_identical(fc$df, c(9, 4))
  # The criterion comparison is compare_sem()'s over the same fits.
  cs <- compare_sem(list(factors_1 = attr(fits[[1]], "fit"), factors_2 = attr(fits[[2]], "fit")))
  at <- match(fc$term, cs$term)
  expect_equal(fc$criterion_delta, cs$delta[at])
  expect_equal(fc$criterion_weight, cs$relative_weight[at])
  expect_equal(fc$criterion_log_evidence, cs$log_evidence_best_over_model[at])
  expect_equal(fc$criterion_delta, fc$bic - min(fc$bic))
  expect_identical(attr(fc, "method")$criterion, "BIC")
  expect_s3_class(attr(fc, "comparison"), "dmar_tbl")
  # The difference test is lavaan's test of the nested pair.
  lrt <- lavaan::lavTestLRT(attr(fits[[1]], "fit"), attr(fits[[2]], "fit"))
  expect_true(is.na(fc$delta_chi_square[1]) && is.na(fc$p_delta[1]))
  expect_equal(fc$delta_chi_square[2], lrt[2, "Chisq diff"])
  expect_equal(fc$delta_df[2], 5)
  expect_equal(fc$p_delta[2], lrt[2, "Pr(>Chisq)"])
  expect_equal(fc$delta_chi_square[2], fc$chi_square[1] - fc$chi_square[2])
  # Eigenvalues and references are parallel_analysis()'s with the same seed.
  pa <- parallel_analysis(d, B = 30, seed = 113)
  expect_equal(fc$eigenvalue, pa$estimate[1:2])
  expect_equal(fc$reference_quantile, pa$reference_quantile[1:2])
  expect_identical(fc$exceeds_reference, pa$exceeds_reference[1:2])
  expect_identical(attr(fc, "method")$suggested_factors, attr(pa, "suggested_factors"))
  expect_s3_class(attr(fc, "parallel_analysis"), "dmar_parallel_analysis")
  # The MAP: the average squared partial correlation after m components.
  R <- cor(d)
  map <- attr(fc, "map")
  expect_identical(names(map), paste0("m_", 0:4))
  expect_equal(unname(map[1]), mean(R[upper.tri(R)]^2))
  e <- eigen(R)
  A <- e$vectors[, 1, drop = FALSE] * sqrt(e$values[1])
  C <- R - A %*% t(A)
  Rs <- C / sqrt(outer(diag(C), diag(C)))
  expect_equal(unname(map[2]), mean(Rs[upper.tri(Rs)]^2))
  expect_equal(fc$map, unname(map[2:3]))
  expect_identical(attr(fc, "method")$map_minimum, unname(which.min(map)) - 1L)
  # (The comparison with psych's VSS()$map is tools/oracle_snippets/psych_map.R,
  # since psych is not a declared dependency.)
  # Attributes and tidy().
  expect_length(attr(fc, "fits"), 2)
  expect_s4_class(attr(fc, "fits")[[2]], "lavaan")
  expect_true(all(c("delta_df", "sabic", "exceeds_reference") %in% names(fc)))
  expect_s3_class(generics::tidy(fc), "data.frame")
  # The display.
  local_reproducible_output(width = 200)
  out <- capture.output(print(fc))
  expect_match(out, "^Factor counts 1 to 2 for 6 items, N = 301$", all = FALSE)
  expect_match(out, "^ *factors +solution +chi_square +df +p_value +rmsea +cfi +tli +srmr +bic +delta +weight +p_delta +eigenvalue +reference +map *$", all = FALSE)
  expect_match(out, "^BIC comparison over the proper solutions \\(compare_sem\\(\\)\\): lowest at 2 factors, weight 1\\.000; the nearest other count, 1 factor, is 72\\.6 higher \\(weight 0\\.000, log Bayes factor 36\\.3 for 2 factors over it\\)\\.$", all = FALSE)
  expect_match(out, "with delta \\(its gap from the lowest\\) and weight", all = FALSE)
  # An RMSEA of 0 with its interval (at the default 90 percent), and a CFI
  # of 1, keep their decimals.
  expect_match(out, "^ *2 +proper +2\\.109 +4 +0\\.7157 +0\\.000 \\[0\\.000, 0\\.064\\] +1\\.000 +1\\.011", all = FALSE)
  expect_match(out, "^Information criteria over the proper solutions: lowest AIC at 2 factors", all = FALSE)
  expect_match(out, "^Parallel analysis \\(normal reference, B = 30, 95th percentile, reference N = 301\\)", all = FALSE)
  expect_match(out, "^Velicer's MAP: the average squared partial correlation is smallest", all = FALSE)
  expect_match(out, "No cutoff is applied", all = FALSE)
  expect_match(out, "Confidence level: 90%", all = FALSE)
  expect_false(any(grepl("Without a proper solution", out)))
  # A subset prints as an ordinary table.
  sub <- capture.output(print(fc[, c("term", "factors", "bic")]))
  expect_false(any(grepl("Velicer", sub)))
  skip_if_not_installed("knitr")
  expect_match(as.character(knitr::knit_print(fc)), "Velicer's MAP", all = FALSE)
})

test_that("factor_count() reports counts without a proper solution and keeps going", {
  skip_if_not_installed("lavaan")
  nine <- grep("^t[1-9]_", names(holzinger_swineford), value = TRUE)
  fc <- factor_count(holzinger_swineford, items = nine, B = 0)
  # Nine items: five factors keep one degree of freedom.
  expect_identical(fc$factors, 1:5)
  expect_identical(fc$solution[1:2], c("proper", "proper"))
  expect_true(all(fc$solution[3:5] %in% c("improper", "not converged")))
  expect_true(all(is.na(fc$chi_square[3:5])))
  expect_true(all(is.finite(fc$chi_square[1:2])))
  expect_true(all(is.na(fc$reference_quantile)))
  expect_true(all(is.finite(fc$criterion_delta[1:2])) && all(is.na(fc$criterion_delta[3:5])))
  expect_identical(nrow(attr(fc, "comparison")), 2L)
  expect_true(all(is.finite(fc$eigenvalue)))
  expect_true(all(is.finite(fc$map)))
  expect_true(is.null(attr(fc, "fits")[[3]]))
  reasons <- attr(fc, "reasons")
  expect_true(all(is.na(reasons[1:2])) && all(!is.na(reasons[3:5])))
  expect_match(reasons[["factors_3"]], "improper solution|did not converge")
  local_reproducible_output(width = 200)
  out <- capture.output(print(fc))
  expect_match(out, "^Without a proper solution, so their fit is not reported: 3 factors \\((the optimizer did not converge|improper: )", all = FALSE)
  expect_match(out, "4 factors \\((the optimizer did not converge|improper: )", all = FALSE)
  expect_false(any(grepl("^Parallel analysis", out)))
  expect_false(any(grepl(" reference ", out[grepl("^ *factors", out)])))
  # An explicit count beyond the default is accepted up to zero degrees of
  # freedom and refused past it.
  d <- holzinger_swineford[fc_items()]
  three <- factor_count(d, max_factors = 3, B = 0)
  expect_identical(three$factors, 1:3)
  expect_error(factor_count(d, max_factors = 4, B = 0), "largest is 3")
  expect_error(factor_count(d, max_factors = 0, B = 0), "whole number")
  expect_error(factor_count(d, B = -1), "'B'")
  expect_error(factor_count(d[1:3], B = 0), "at least four")
  # An input problem stops the comparison rather than filling the rows.
  expect_error(factor_count(d, estimator = "WLSMV", B = 0), "ML or MLR")
  expect_error(factor_count(d, criterion = "SPBIC", B = 0), "arg")
  # Another criterion: AIC has weights but no log evidence.
  aic <- factor_count(d, criterion = "aic", B = 0)
  expect_identical(attr(aic, "method")$criterion, "AIC")
  expect_equal(aic$criterion_delta, aic$aic - min(aic$aic))
  expect_true(all(is.na(aic$criterion_log_evidence)))
  expect_true(all(is.finite(aic$criterion_weight)))
  local_reproducible_output(width = 200)
  out_aic <- capture.output(print(aic))
  expect_match(out_aic, "^ *factors +solution +chi_square +df +p_value +rmsea +cfi +tli +srmr +aic +delta +weight +p_delta +eigenvalue +map *$", all = FALSE)
  expect_match(out_aic, "^AIC comparison over the proper solutions .*higher \\(weight [0-9.e-]+\\)\\.$", all = FALSE)
})

test_that("factor_count() takes a covariance matrix and the ULS and MLR estimators", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[fc_items()]
  # Complete data fitted listwise and its covariance matrix agree in every
  # column; the FIML fit adds the mean structure, which moves the
  # information criteria but nothing else.
  from_data <- factor_count(d, missing = "listwise", B = 0)
  from_S <- factor_count(S = cov(d), N = nrow(d), B = 0)
  for (col in c("chi_square", "df", "rmsea", "cfi", "srmr", "aic", "bic", "eigenvalue", "map")) {
    expect_equal(from_S[[col]], from_data[[col]], tolerance = 1e-6, info = col)
  }
  fiml <- factor_count(d, B = 0)
  for (col in c("chi_square", "df", "rmsea", "cfi", "srmr", "eigenvalue", "map")) {
    expect_equal(fiml[[col]], from_data[[col]], tolerance = 1e-6, info = col)
  }
  expect_false(isTRUE(all.equal(fiml$aic, from_data$aic)))
  expect_identical(attr(from_S, "method")$missing, "covariance input")
  uls <- factor_count(d, estimator = "ULS", B = 0)
  expect_true(all(is.finite(uls$chi_square)) && all(is.finite(uls$srmr)))
  expect_true(all(is.na(c(uls$p_value, uls$rmsea, uls$cfi, uls$aic, uls$p_delta,
                          uls$criterion_delta, uls$criterion_weight))))
  expect_null(attr(uls, "comparison"))
  local_reproducible_output(width = 200)
  out <- capture.output(print(uls))
  expect_match(out, "^ *factors +solution +chi_square +df +srmr +eigenvalue +map *$", all = FALSE)
  expect_match(out, "Under ULS the statistic has no reference distribution", all = FALSE)
  mlr <- factor_count(d, estimator = "MLR", B = 0)
  fits <- lapply(1:2, function(k) attr(fa_exploratory(d, k, rotation = "none", estimator = "MLR"), "fit"))
  lrt <- suppressMessages(lavaan::lavTestLRT(fits[[1]], fits[[2]]))
  expect_equal(mlr$delta_chi_square[2], lrt[2, "Chisq diff"])
  expect_equal(mlr$p_delta[2], lrt[2, "Pr(>Chisq)"])
  expect_identical(attr(mlr, "method")$estimator, "MLR")
})

test_that("factor_count() reduces a covariance matrix to the named items", {
  skip_if_not_installed("lavaan")
  nine <- grep("^t[1-9]_", names(holzinger_swineford), value = TRUE)
  S9 <- cov(holzinger_swineford[nine])
  N <- nrow(holzinger_swineford)
  # Six of the nine items, in an order other than that of S. Naming them
  # analyzes the same matrix as reducing S to them first.
  items <- rev(fc_items())
  named <- factor_count(S = S9, N = N, items = items, B = 20, seed = 113)
  reduced <- factor_count(S = S9[items, items], N = N, B = 20, seed = 113)
  expect_identical(attr(named, "method")$p, 6L)
  expect_identical(named$factors, 1:2)
  for (col in names(named)) expect_equal(named[[col]], reduced[[col]], info = col)
  expect_equal(attr(named, "map"), attr(reduced, "map"))
  expect_equal(attr(named, "method"), attr(reduced, "method"))
  expect_equal(attr(named, "parallel_analysis"), attr(reduced, "parallel_analysis"))
  expect_identical(nrow(attr(named, "parallel_analysis")), 6L)
  expect_identical(lavaan::lavNames(attr(named, "fits")[[2]], "ov"), items)
  # Without items every item of S is analyzed.
  expect_identical(attr(factor_count(S = S9, N = N, max_factors = 1, B = 0), "method")$p, 9L)
  # The same items from the raw data give the same analysis (complete data
  # fitted listwise, as in the covariance matrix).
  from_data <- factor_count(holzinger_swineford, items = items, missing = "listwise",
                            B = 20, seed = 113)
  for (col in c("chi_square", "df", "rmsea", "cfi", "srmr", "aic", "bic",
                "eigenvalue", "reference_quantile", "map")) {
    expect_equal(named[[col]], from_data[[col]], tolerance = 1e-6, info = col)
  }
  # A name that is not an item of S is refused, as are too few items and a
  # repeated one; fa_exploratory() reads the names the same way, and raw
  # data refuse a name that is not a column. Each refusal carries the
  # class dmar_items_invalid.
  expect_error(factor_count(S = S9, N = N, items = c(items[1:5], "t99_unknown"), B = 0),
               "at least three distinct items of 'S'. Not items of 'S': 't99_unknown'.",
               fixed = TRUE, class = "dmar_items_invalid")
  expect_error(factor_count(S = S9, N = N, items = c(items[1:5], ""), B = 0),
               "Not items of 'S': ''.", fixed = TRUE, class = "dmar_items_invalid")
  expect_error(factor_count(S = S9, N = N, items = items[1:2], B = 0),
               "at least three distinct items of 'S'.", fixed = TRUE,
               class = "dmar_items_invalid")
  expect_error(factor_count(S = S9, N = N, items = c(items, items[1]), B = 0),
               "at least three distinct items of 'S'.", fixed = TRUE,
               class = "dmar_items_invalid")
  expect_error(fa_exploratory(S = S9, N = N, items = c(items[1:5], "t99_unknown")),
               "Not items of 'S': 't99_unknown'.", fixed = TRUE,
               class = "dmar_items_invalid")
  expect_error(factor_count(holzinger_swineford, items = c(items[1:5], "t99_unknown"), B = 0),
               "at least three distinct columns of data", fixed = TRUE,
               class = "dmar_items_invalid")
})

test_that("a set of items that cannot be read is one refusal class, with the check that failed", {
  skip_if_not_installed("lavaan")
  nine <- grep("^t[1-9]_", names(holzinger_swineford), value = TRUE)
  S9 <- cov(holzinger_swineford[nine])
  N <- nrow(holzinger_swineford)
  items <- fc_items()
  d <- holzinger_swineford[items]
  items_refusal <- function(expr, reason, named = character(0)) {
    e <- tryCatch(expr, dmar_items_invalid = function(e) e)
    expect_s3_class(e, c("dmar_items_invalid", "dmar_refusal", "error"))
    expect_identical(e$code, "items_invalid")
    expect_identical(e$reason, reason)
    expect_identical(e$items, named)
    invisible(e)
  }
  # The covariance path and the raw data path give the same reason for the
  # same names: the first check that fails, in the order form, unknown,
  # duplicate, count.
  for (input in c("S", "data")) {
    call_fc <- function(items) {
      if (input == "S") factor_count(S = S9, N = N, items = items, B = 0)
      else factor_count(d, items = items, B = 0)
    }
    items_refusal(call_fc(c(items[1:5], "t99_unknown", "t99_unknown")),
                  "unknown", "t99_unknown")
    items_refusal(call_fc(c(items[1:2], "")), "unknown", "")
    items_refusal(call_fc(c(items[1:2], "t99_unknown")), "unknown", "t99_unknown")
    items_refusal(call_fc(c(items, items[2], items[2])), "duplicate", items[2])
    items_refusal(call_fc(items[1:2]), "count")
    items_refusal(call_fc(c(items[1:3], NA)), "form")
    items_refusal(call_fc(1:3), "form")
  }
  # The message of each path is the one it gave before the class.
  e <- items_refusal(factor_count(S = S9, N = N, items = items[1:2], B = 0), "count")
  expect_identical(conditionMessage(e), "'items' must name at least three distinct items of 'S'.")
  e <- items_refusal(factor_count(d, items = items[1:2], B = 0), "count")
  expect_identical(conditionMessage(e), paste0("'items' must name at least three ",
    "distinct columns of data, whose column names must be unique."))
  # Raw data whose columns repeat a name are refused for their names,
  # whether or not the repeated column is among the items; with items =
  # NULL every column is an item, and an empty name is refused the same way.
  twice <- d
  names(twice)[2] <- names(twice)[1]
  items_refusal(factor_count(twice, B = 0), "data_names", items[1])
  items_refusal(factor_count(twice, items = items[3:6], B = 0), "data_names", items[1])
  unnamed <- d
  names(unnamed)[3] <- ""
  items_refusal(factor_count(unnamed, B = 0), "data_names", "")
  expect_s3_class(factor_count(unnamed, items = items[-3], max_factors = 1, B = 0),
                  "dmar_tbl")
  # parallel_analysis() and fa_exploratory() read the names through the
  # same checks.
  items_refusal(parallel_analysis(d, items = c(items[1:2], items[1]), B = 2), "duplicate", items[1])
  items_refusal(parallel_analysis(S = S9, N = N, items = items[1:2], B = 2), "count")
  items_refusal(fa_exploratory(d, 1, items = c(items[1:5], "t99_unknown")),
                "unknown", "t99_unknown")
  items_refusal(fa_exploratory(S = S9, N = N, nfactors = 1, items = 1:3), "form")
})

test_that("factor_count() checks only the named items of a covariance matrix", {
  skip_if_not_installed("lavaan")
  nine <- grep("^t[1-9]_", names(holzinger_swineford), value = TRUE)
  S9 <- cov(holzinger_swineford[nine])
  N <- nrow(holzinger_swineford)
  items <- fc_items()
  other <- setdiff(nine, items)
  same_analysis <- function(x, y) {
    for (col in names(y)) expect_equal(x[[col]], y[[col]], info = col)
    for (a in c("map", "method", "parallel_analysis")) {
      expect_equal(attr(x, a), attr(y, a), info = a)
    }
  }
  reduced <- factor_count(S = S9[items, items], N = N, B = 20, seed = 113)
  # Entries outside the named items are not read, as the columns of raw
  # data not named are not: a missing covariance, a variance of zero, and
  # an asymmetric pair, which a matrix holding other variables can have.
  holed <- S9
  holed[other[1], other[2]] <- holed[other[2], other[1]] <- NA
  holed[other[3], other[3]] <- 0
  holed[other[1], other[3]] <- holed[other[1], other[3]] + 1
  expect_error(factor_count(S = holed, N = N, B = 0),
               "'S' must be a finite symmetric covariance matrix", fixed = TRUE)
  same_analysis(factor_count(S = holed, N = N, items = items, B = 20, seed = 113),
                reduced)
  # A whole matrix that is indefinite, as one from pairwise deletion can
  # be, while the six-item submatrix is positive definite.
  indefinite <- S9
  indefinite[other[1], other[2]] <- indefinite[other[2], other[1]] <-
    1.5 * sqrt(S9[other[1], other[1]] * S9[other[2], other[2]])
  expect_lt(min(eigen(indefinite, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_error(factor_count(S = indefinite, N = N, B = 0),
               "'S' must be positive definite", fixed = TRUE)
  same_analysis(factor_count(S = indefinite, N = N, items = items, B = 20, seed = 113),
                reduced)
  # A problem among the named items is still refused, and the message says
  # it was found in the reduced matrix.
  inside <- S9
  inside[items[1], items[2]] <- inside[items[2], items[1]] <- NA
  expect_error(factor_count(S = inside, N = N, items = items, B = 0),
               "'S' reduced to 'items' must be a finite symmetric", fixed = TRUE)
  inside <- S9
  inside[items[1], items[2]] <- inside[items[2], items[1]] <-
    1.5 * sqrt(S9[items[1], items[1]] * S9[items[2], items[2]])
  expect_error(factor_count(S = inside, N = N, items = items, B = 0),
               "'S' reduced to 'items' must be positive definite", fixed = TRUE)
  # The names are read from the whole matrix, so they must still be sound.
  unnamed <- holed
  dimnames(unnamed) <- NULL
  expect_error(factor_count(S = unnamed, N = N, items = items, B = 0),
               "matching, unique row and column names", fixed = TRUE)
  # fa_exploratory() reads covariance input the same way.
  efa <- fa_exploratory(S = holed, N = N, nfactors = 2, items = items, seed = 113)
  efa_reduced <- fa_exploratory(S = S9[items, items], N = N, nfactors = 2, seed = 113)
  for (col in names(efa_reduced)) expect_equal(efa[[col]], efa_reduced[[col]], info = col)
  for (a in c("method", "fit_measures", "variance_accounted", "factor_quality")) {
    expect_equal(attr(efa, a), attr(efa_reduced, a), info = a)
  }
})

test_that("the RMSEA intervals are at conf_level, 0.90 by default, which reaches nothing else", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two comparisons
  # Every RMSEA interval defaults to 90 percent, whose limits decide the
  # tests of close and not close fit at .05. conf_level sets the RMSEA
  # intervals of the counts alone; the Monte Carlo limits of the parallel
  # analysis reference keep parallel_analysis()'s 95 percent.
  d <- holzinger_swineford[fc_items()]
  expect_equal(eval(formals(factor_count)$conf_level), 0.90)
  x <- factor_count(d, B = 10, seed = 113)
  y <- factor_count(d, B = 10, seed = 113, conf_level = 0.95)
  for (k in 1:2) {
    fm <- sem_fit_measures(attr(fa_exploratory(d, k, rotation = "none"), "fit"))
    fm95 <- sem_fit_measures(attr(fa_exploratory(d, k, rotation = "none"), "fit"),
                             conf_level = 0.95)
    lim <- c("rmsea_ci_lower", "rmsea_ci_upper")
    expect_equal(c(x$rmsea_ci_lower[k], x$rmsea_ci_upper[k]),
                 fm$value[match(lim, fm$term)])
    expect_equal(c(y$rmsea_ci_lower[k], y$rmsea_ci_upper[k]),
                 fm95$value[match(lim, fm95$term)])
  }
  expect_equal(attr(x, "conf_level"), 0.90)
  expect_identical(dmar_metadata(x)$rmsea_conf_level, 0.90)
  expect_identical(attr(attr(x, "parallel_analysis"), "method")$mc_conf_level,
                   0.95)
  expect_identical(attr(x, "parallel_analysis"), attr(y, "parallel_analysis"))
  out <- capture.output(print(x))
  expect_match(out, "RMSEA with 90% interval", all = FALSE, fixed = TRUE)
  expect_match(out, "Confidence level: 90%", all = FALSE, fixed = TRUE)
  # tidy() states the level of the RMSEA intervals beside them, with the
  # count kept as the estimate, and glance() gives it too.
  tx <- generics::tidy(x)
  ty <- generics::tidy(y)
  expect_identical(tx$estimate, c(1, 2))
  expect_equal(tx$rmsea_ci_upper, x$rmsea_ci_upper)
  expect_identical(tx$conf_level, c(0.90, 0.90))
  expect_identical(ty$conf_level, c(0.95, 0.95))
  expect_identical(generics::glance(x)$conf_level, 0.90)
})

# Six items cut into four ordered categories, n = 400: the fixture of the
# 2026-11-06 deep review (finding F12), whose items are the first six of
# the nine it generates (the same draws).
fc_ordered_fixture <- function() {
  withr::with_seed(113, {
    n <- 400L
    g <- rnorm(n)
    groups <- matrix(rnorm(n * 3), n, 3)
    X <- sapply(seq_len(6), function(j) .5 * g + .4 * groups[, ceiling(j / 3)] +
                  sqrt(.59) * rnorm(n))
    colnames(X) <- paste0("x", seq_len(6))
    as.data.frame(lapply(as.data.frame(X), function(v)
      as.integer(cut(v, breaks = c(-Inf, -.7, 0, .7, Inf)))))
  })
}
# The same items with 40 values removed (39 cases incomplete).
fc_ordered_incomplete <- function() {
  ord <- fc_ordered_fixture()
  withr::with_seed(7, {
    ord[cbind(sample.int(nrow(ord), 40, replace = TRUE),
              sample.int(6, 40, replace = TRUE))] <- NA
  })
  ord
}
# The columns of a table without its attributes (the lavaan fits carry
# their own timings).
fc_columns <- function(x) {
  y <- as.data.frame(x)
  attributes(y) <- attributes(y)[c("names", "row.names")]
  y
}

test_that("complete ordered items left to the default are fitted as fa_exploratory() fits them, listwise", {
  skip_if_not_installed("lavaan")
  # The 2026-11-06 deep review (F12): fa_exploratory() fits complete
  # ordered items with 'missing' left out as they stand (listwise deletion
  # with nothing to delete), but factor_count() matched its own formal
  # default, "fiml", passed it on explicitly, and was refused ("Ordered
  # EFA requires ... an explicit missing = ..."), while missing =
  # "listwise" ran. Both functions now resolve the default once, in the
  # same way.
  ord <- fc_ordered_fixture()
  omitted <- factor_count(ord, estimator = "WLSMV", ordered = TRUE, B = 0)
  listwise <- factor_count(ord, estimator = "WLSMV", ordered = TRUE,
                           missing = "listwise", B = 0)
  expect_identical(fc_columns(omitted), fc_columns(listwise))
  expect_identical(attr(omitted, "method"), attr(listwise, "method"))
  expect_identical(attr(omitted, "method")$missing, "listwise")
  expect_identical(omitted$solution, c("proper", "proper"))
  for (k in 1:2) {
    expect_identical(lavaan::coef(attr(omitted, "fits")[[k]]),
                     lavaan::coef(attr(listwise, "fits")[[k]]))
  }
  # The rotation's random starts take the seed, so that the two calls
  # rotate from the same starts.
  efa <- fa_exploratory(ord, 2, estimator = "WLSMV", ordered = TRUE, seed = 113)
  expect_identical(attr(efa, "method")$missing, "listwise")
  expect_identical(fc_columns(efa), fc_columns(fa_exploratory(
    ord, 2, estimator = "WLSMV", ordered = TRUE, missing = "listwise", seed = 113)))
  # The oracle is lavaan called directly: the one-factor model as lavaan's
  # CFA with one factor (a route other than lavaan::efa()), the
  # two-factor model through lavaan::efa(). The chi_square column is the
  # scaled (mean and variance adjusted) statistic, the RMSEA and CFI the
  # robust ones.
  items <- names(ord)
  cfa1 <- lavaan::cfa(paste("f =~", paste(items, collapse = " + ")), data = ord,
                      ordered = items, estimator = "WLSMV", std.lv = TRUE)
  efa2 <- lavaan::efa(data = ord, nfactors = 2, ordered = items,
                      estimator = "WLSMV", rotation = "none")
  measures <- c("chisq.scaled", "df", "rmsea.robust", "cfi.robust", "srmr")
  fm1 <- lavaan::fitMeasures(cfa1, measures)
  fm2 <- lavaan::fitMeasures(efa2, measures)
  expect_equal(omitted$chi_square, c(fm1[["chisq.scaled"]], fm2["chisq.scaled", 1]),
               tolerance = 1e-6)
  expect_equal(omitted$df, c(fm1[["df"]], fm2["df", 1]))
  expect_equal(omitted$rmsea, c(fm1[["rmsea.robust"]], fm2["rmsea.robust", 1]),
               tolerance = 1e-6)
  expect_equal(omitted$cfi, c(fm1[["cfi.robust"]], fm2["cfi.robust", 1]),
               tolerance = 1e-6)
  expect_equal(omitted$srmr, c(fm1[["srmr"]], fm2["srmr", 1]), tolerance = 1e-6)
  # The eigenvalues are those of the polychoric correlations with the
  # communalities of the one-factor maximum likelihood solution on the
  # diagonal, formed here from lavaan's polychoric matrix and
  # stats::factanal().
  R <- as.matrix(lavaan::lavCor(ord, ordered = items, output = "cor"))
  f1 <- stats::factanal(covmat = list(cov = R, n.obs = nrow(ord)), factors = 1)
  diag(R) <- 1 - f1$uniquenesses
  expect_equal(omitted$eigenvalue,
               eigen(R, symmetric = TRUE, only.values = TRUE)$values[1:2],
               tolerance = 1e-6)
})

test_that("complete ordered items left to the default reach the parallel analysis listwise", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # reference data sets are simulated
  # The parallel analysis branch of F12 (DR27): parallel_analysis() takes
  # ordered items only with an explicit pairwise or listwise choice, and
  # factor_count() passes it the method resolved for the fits. A B = 0
  # call does not reach it.
  ord <- fc_ordered_fixture()
  omitted <- factor_count(ord, estimator = "WLSMV", ordered = TRUE, B = 5, seed = 113)
  listwise <- factor_count(ord, estimator = "WLSMV", ordered = TRUE,
                           missing = "listwise", B = 5, seed = 113)
  expect_identical(fc_columns(omitted), fc_columns(listwise))
  expect_identical(attr(omitted, "method"), attr(listwise, "method"))
  pa <- attr(omitted, "parallel_analysis")
  expect_identical(fc_columns(pa), fc_columns(attr(listwise, "parallel_analysis")))
  expect_identical(attr(omitted, "method")$missing, "listwise")
  expect_identical(attr(pa, "method")$missing, "listwise")
  expect_identical(attr(pa, "counts")[c("attempted", "usable", "failed")],
                   list(attempted = 5L, usable = 5L, failed = 0L))
  expect_true(all(is.finite(omitted$reference_quantile)))
  # The oracle is the ordinal normal reference written from its
  # description: after the seed, B blocks of N x p standard normal
  # deviates, each item cut at the normal quantiles of its observed
  # cumulative category proportions, polychoric correlations, the reduced
  # eigenvalues of the one-factor maximum likelihood solution, and the
  # 95th percentile of each rank (quantile type 7).
  reference <- withr::with_seed(113, {
    N <- nrow(ord)
    p <- ncol(ord)
    draws <- t(vapply(seq_len(5), function(b) {
      x <- matrix(stats::rnorm(N * p), N, p)
      for (j in seq_len(p)) {
        prop <- as.numeric(table(ord[[j]])) / N
        cuts <- stats::qnorm(c(0, cumsum(prop)[-length(prop)], 1))
        x[, j] <- sort(unique(ord[[j]]))[findInterval(x[, j], cuts, all.inside = TRUE)]
      }
      x <- as.data.frame(x)
      names(x) <- names(ord)
      R <- as.matrix(lavaan::lavCor(x, ordered = names(x), output = "cor"))
      f1 <- stats::factanal(covmat = list(cov = R, n.obs = N), factors = 1)
      diag(R) <- 1 - f1$uniquenesses
      eigen(R, symmetric = TRUE, only.values = TRUE)$values
    }, numeric(p)))
    apply(draws, 2L, stats::quantile, probs = 0.95, names = FALSE, type = 7)
  })
  expect_equal(pa$reference_quantile, reference, tolerance = 1e-6)
  expect_equal(omitted$reference_quantile, reference[1:2], tolerance = 1e-6)
})

test_that("incomplete ordered items left to the default are refused before anything is computed", {
  skip_if_not_installed("lavaan")
  ord_na <- fc_ordered_incomplete()
  expect_identical(sum(!stats::complete.cases(ord_na)), 39L)
  # fa_exploratory() refuses these items with 'missing' left out, and
  # factor_count() refuses them with the same message.
  efa_refusal <- tryCatch(fa_exploratory(ord_na, 1, estimator = "WLSMV", ordered = TRUE),
                          error = function(e) e)
  expect_s3_class(efa_refusal, "dmar_argument_invalid")
  # Nothing may be reached: the correlations, the fits, the parallel
  # analysis, and the random number stream are all left alone.
  reached <- function(...) stop("reached before the refusal")
  local_mocked_bindings(.pa_cor = reached, cov_fiml = reached,
                        fa_exploratory = reached, parallel_analysis = reached)
  for (B in c(0, 200)) {
    set.seed(99)
    before <- get(".Random.seed", envir = globalenv())
    e <- tryCatch(factor_count(ord_na, estimator = "WLSMV", ordered = TRUE, B = B),
                  error = function(e) e)
    expect_s3_class(e, c("dmar_argument_invalid", "dmar_refusal", "error"))
    expect_match(conditionMessage(e),
                 "explicit missing = \"pairwise\" or \"listwise\" choice", fixed = TRUE)
    expect_identical(conditionMessage(e), conditionMessage(efa_refusal))
    expect_identical(get(".Random.seed", envir = globalenv()), before)
  }
})

test_that("incomplete ordered items are fitted under the explicit pairwise or listwise choice", {
  skip_if_not_installed("lavaan")
  ord_na <- fc_ordered_incomplete()
  items <- names(ord_na)
  model <- paste("f =~", paste(items, collapse = " + "))
  pairwise <- factor_count(ord_na, estimator = "WLSMV", ordered = TRUE,
                           missing = "pairwise", max_factors = 1, B = 0)
  listwise <- factor_count(ord_na, estimator = "WLSMV", ordered = TRUE,
                           missing = "listwise", max_factors = 1, B = 0)
  expect_identical(attr(pairwise, "method")$missing, "pairwise")
  expect_identical(attr(pairwise, "method")$N, 400L)
  expect_identical(attr(listwise, "method")$missing, "listwise")
  expect_identical(attr(listwise, "method")$N, 361L)
  # Each equals lavaan called directly on the same cases.
  fit_pw <- lavaan::cfa(model, data = ord_na, ordered = items, estimator = "WLSMV",
                        missing = "pairwise", std.lv = TRUE)
  fit_lw <- lavaan::cfa(model, data = ord_na[stats::complete.cases(ord_na), ],
                        ordered = items, estimator = "WLSMV", std.lv = TRUE)
  expect_equal(pairwise$chi_square, lavaan::fitMeasures(fit_pw, "chisq.scaled")[[1]],
               tolerance = 1e-6)
  expect_equal(listwise$chi_square, lavaan::fitMeasures(fit_lw, "chisq.scaled")[[1]],
               tolerance = 1e-6)
  # An explicit FIML is refused for ordered items, complete or not, as
  # fa_exploratory() refuses it.
  ord <- fc_ordered_fixture()
  expect_error(factor_count(ord, estimator = "WLSMV", ordered = TRUE, missing = "fiml",
                            B = 0),
               "Ordered EFA requires", class = "dmar_argument_invalid")
  expect_error(fa_exploratory(ord, 1, estimator = "WLSMV", ordered = TRUE,
                              missing = "fiml"),
               "Ordered EFA requires", class = "dmar_argument_invalid")
})
