# fa_exploratory(point_estimates_only = TRUE): the estimates of the full
# analysis, with no standard errors, tests, or intervals computed. The
# other tests of fa_exploratory() are in test-factor_analysis.R and
# test-fa_exploratory_target.R.

po_items <- function() c("t1_visual_perception", "t2_cubes", "t4_lozenges",
  "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")

# The table and its attributes without the lavaan object, which carries
# timings of its own.
po_strip <- function(x) {
  attr(x, "fit") <- NULL
  x
}

# The full table `full` and the point estimates table `point` of one
# analysis: the same rows and estimates (to 1e-10), every standard error,
# statistic, p-value, and limit NA, the same derived tables, the same fit
# measures, and the setting recorded.
po_expect_point <- function(full, point, info = NULL) {
  expect_identical(point$term, full$term, info = info)
  expect_identical(point$metric, full$metric, info = info)
  expect_equal(point$estimate, full$estimate, tolerance = 1e-10, info = info)
  inference <- c("se", "statistic", "p_value", "ci_lower", "ci_upper")
  expect_true(all(is.na(as.matrix(as.data.frame(point)[inference]))), info = info)
  expect_identical(point$conf_level, full$conf_level, info = info)
  expect_equal(attr(point, "variance_accounted"), attr(full, "variance_accounted"),
               tolerance = 1e-10, info = info)
  qf <- attr(full, "factor_quality")
  qp <- attr(point, "factor_quality")
  for (col in c("determinacy", "min_score_correlation", "H")) {
    expect_equal(qp[[col]], qf[[col]], tolerance = 1e-10, info = paste(info, col))
  }
  expect_true(all(is.na(as.matrix(qp[grepl("_(se|lower|upper)$", names(qp))]))),
              info = info)
  fm_full <- attr(full, "fit_measures")
  fm_point <- attr(point, "fit_measures")
  expect_identical(fm_point$term, fm_full$term, info = info)
  expect_equal(fm_point$value, fm_full$value, tolerance = 1e-10, info = info)
  expect_identical(attr(fm_point, "sem_fit")$footer, attr(fm_full, "sem_fit")$footer,
                   info = info)
  expect_false(attr(full, "point_estimates_only"), info = info)
  expect_true(attr(point, "point_estimates_only"), info = info)
  expect_identical(attr(point, "method")$se, "none", info = info)
  expect_identical(lavaan::lavInspect(attr(point, "fit"), "options")$se, "none",
                   info = info)
  method_full <- attr(full, "method")
  method_point <- attr(point, "method")
  method_full$se <- method_point$se <- NULL
  expect_identical(method_point, method_full, info = info)
  expect_true(generics::glance(point)$point_estimates_only, info = info)
  expect_false(generics::glance(full)$point_estimates_only, info = info)
  expect_true(all(generics::tidy(point)$point_estimates_only), info = info)
  invisible(point)
}

test_that("point_estimates_only = FALSE is the default analysis", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[po_items()]
  x <- fa_exploratory(d, 2, seed = 113)
  y <- fa_exploratory(d, 2, seed = 113, point_estimates_only = FALSE)
  expect_identical(po_strip(y), po_strip(x))
  expect_false(attr(x, "point_estimates_only"))
  expect_true(all(is.finite(x$se[x$op == "=~"])))
  expect_false(any(grepl("point_estimates_only", capture.output(print(x)))))
  # efa_k() is the same function, with the same argument.
  expect_identical(names(formals(efa_k)), names(formals(fa_exploratory)))
  expect_false(formals(fa_exploratory)$point_estimates_only)
})

test_that("point estimates alone equal the full analysis under ML, MLR with FIML, and S input", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[po_items()]
  full <- fa_exploratory(d, 2, seed = 113)
  point <- fa_exploratory(d, 2, seed = 113, point_estimates_only = TRUE)
  po_expect_point(full, point, "ML")
  # The display: the method line and the note say what was not computed,
  # and the cells hold the estimates alone.
  out <- capture.output(print(point))
  expect_true(any(grepl("SE: not computed (point_estimates_only = TRUE)", out,
                        fixed = TRUE)))
  expect_true(any(grepl("Standard errors, tests, and intervals were not computed",
                        out, fixed = TRUE)))
  expect_false(any(grepl("[", out[grepl("^ *t1_visual", out)], fixed = TRUE)))
  expect_true(any(grepl("^Fit: chi square", out)))
  # The only interval left is the RMSEA's, whose level the fit line names;
  # the table's level (conf_level) sets nothing and is not stated.
  expect_true(any(grepl("^Fit: .*; RMSEA with 90% interval ", out)))
  expect_false(any(grepl("Confidence level", out, fixed = TRUE)))
  expect_match(attr(point, "analysis_note"), "(point_estimates_only = TRUE)",
               fixed = TRUE)
  # An orthogonal rotation, and the unrotated solution, whose fit is the
  # rotated fit's.
  po_expect_point(fa_exploratory(d, 2, rotation = "varimax", seed = 113),
                  fa_exploratory(d, 2, rotation = "varimax", seed = 113,
                                 point_estimates_only = TRUE), "varimax")
  po_expect_point(fa_exploratory(d, 2, rotation = "none"),
                  fa_exploratory(d, 2, rotation = "none",
                                 point_estimates_only = TRUE), "none")
  # MLR with FIML on incomplete items: the robust test and the FIML
  # corrected indices are unchanged.
  dm <- d
  dm[seq(2, 300, 3), 2] <- NA
  dm[seq(5, 300, 7), 5] <- NA
  full <- fa_exploratory(dm, 2, estimator = "MLR", seed = 113)
  point <- fa_exploratory(dm, 2, estimator = "MLR", seed = 113,
                          point_estimates_only = TRUE)
  po_expect_point(full, point, "MLR FIML")
  expect_identical(attr(point, "method")$missing, "fiml")
  fm <- attr(point, "fit_measures")
  expect_true(all(is.finite(fm$value[fm$term %in% c("chi_square", "rmsea", "cfi",
                                                     "tli", "srmr")])))
  # Covariance input under ML, with rmsea_conf_level setting the RMSEA
  # interval.
  S <- cov(d)
  full <- fa_exploratory(S = S, N = nrow(d), nfactors = 2, seed = 113,
                         rmsea_conf_level = 0.95)
  point <- fa_exploratory(S = S, N = nrow(d), nfactors = 2, seed = 113,
                          rmsea_conf_level = 0.95, point_estimates_only = TRUE)
  po_expect_point(full, point, "S")
  fm <- attr(point, "fit_measures")
  expect_equal(fm$value[fm$term == "rmsea_ci_upper"],
               unname(lavaan::fitMeasures(
                 lavaan::efa(sample.cov = S, sample.nobs = nrow(d), nfactors = 2,
                             rotation = "none", output = "lavaan"),
                 "rmsea.ci.upper", fm.args = list(rmsea.ci.level = 0.95))),
               tolerance = 1e-8)
  expect_true(any(grepl("; RMSEA with 95% interval ",
                        capture.output(print(point)), fixed = TRUE)))
})

test_that("point estimates alone equal the full analysis for ordered items and ULS", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[po_items()]
  o <- d
  for (j in seq_along(o)) o[[j]] <- as.integer(cut(o[[j]],
    unique(quantile(o[[j]], seq(0, 1, length.out = 5))), include.lowest = TRUE))
  full <- fa_exploratory(o, 2, estimator = "WLSMV", ordered = TRUE, seed = 113)
  point <- fa_exploratory(o, 2, estimator = "WLSMV", ordered = TRUE, seed = 113,
                          point_estimates_only = TRUE)
  po_expect_point(full, point, "WLSMV")
  expect_true(any(grepl("of the latent responses", capture.output(print(point)))))
  # ULS from raw data, whose sandwich errors are skipped, and from a
  # supplied matrix, which has none either way.
  po_expect_point(fa_exploratory(d, 2, estimator = "ULS", seed = 113),
                  fa_exploratory(d, 2, estimator = "ULS", seed = 113,
                                 point_estimates_only = TRUE), "ULS raw")
  R <- cor(d)
  po_expect_point(fa_exploratory(S = R, N = nrow(d), nfactors = 2, estimator = "ULS",
                                 seed = 113),
                  fa_exploratory(S = R, N = nrow(d), nfactors = 2, estimator = "ULS",
                                 seed = 113, point_estimates_only = TRUE), "ULS S")
  # ULS has no RMSEA, so with no interval at all no confidence level is
  # stated.
  out <- capture.output(print(fa_exploratory(d, 2, estimator = "ULS", seed = 113,
                                             point_estimates_only = TRUE)))
  expect_false(any(grepl("Confidence level", out, fixed = TRUE)))
})

test_that("point estimates alone equal the full analysis under target rotation", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  items <- po_items()
  d <- holzinger_swineford[items]
  # Overidentified: every loading estimated.
  keys <- list(spatial = items[1:3], verbal = items[4:6])
  po_expect_point(fa_exploratory(d, rotation = "target", target = keys, seed = 113),
                  fa_exploratory(d, rotation = "target", target = keys, seed = 113,
                                 point_estimates_only = TRUE), "target keys")
  # Just identified: both zeros attained exactly and shown at zero in both
  # metrics, as in the full analysis.
  B <- matrix(NA_real_, 6, 2, dimnames = list(items, c("spatial", "verbal")))
  B["t7_sentence", "spatial"] <- 0
  B["t4_lozenges", "verbal"] <- 0
  full <- fa_exploratory(d, rotation = "target", target = B, seed = 113)
  point <- fa_exploratory(d, rotation = "target", target = B, seed = 113,
                          point_estimates_only = TRUE)
  po_expect_point(full, point, "target just identified")
  fixed <- point[point$op == "=~" & paste(point$lhs, point$rhs) %in%
                   c("f1 t7_sentence", "f2 t4_lozenges"), ]
  expect_identical(nrow(fixed), 4L)
  expect_identical(fixed$estimate, rep(0, 4))
  expect_match(attr(point, "analysis_note"),
               "2 of the 2 specified loadings are attained exactly", fixed = TRUE)
  # Orthogonal, with one zero attained.
  Bo <- B
  Bo["t7_sentence", "spatial"] <- NA
  full <- fa_exploratory(d, rotation = "target_orthogonal", target = Bo, seed = 113)
  point <- fa_exploratory(d, rotation = "target_orthogonal", target = Bo, seed = 113,
                          point_estimates_only = TRUE)
  po_expect_point(full, point, "target orthogonal")
  expect_match(attr(point, "analysis_note"),
               "1 of the 1 specified loadings is attained exactly", fixed = TRUE)
})

test_that("point_estimates_only keeps every refusal that does not concern the errors", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[po_items()]
  for (bad in list(NA, "yes", c(TRUE, FALSE), 1, NULL)) {
    expect_error(fa_exploratory(d, 2, point_estimates_only = bad),
                 "'point_estimates_only' must be TRUE or FALSE", fixed = TRUE)
  }
  # The target checks run before any fit.
  B <- matrix(NA_real_, 6, 2, dimnames = list(po_items(), c("spatial", "verbal")))
  B["t7_sentence", "spatial"] <- 0
  e <- tryCatch(fa_exploratory(d, rotation = "target", target = B,
                               point_estimates_only = TRUE),
                error = function(e) e)
  expect_s3_class(e, "dmar_target_unidentified")
  expect_identical(e$reason, "counts")
  e <- tryCatch(fa_exploratory(d, 2, target = B, point_estimates_only = TRUE),
                error = function(e) e)
  expect_s3_class(e, "dmar_target_invalid")
  # The input checks are those of the full analysis.
  expect_error(fa_exploratory(d, 4, point_estimates_only = TRUE), "identified EFA")
  expect_error(fa_exploratory(S = cov(d), N = nrow(d), estimator = "MLR",
                              point_estimates_only = TRUE), "raw data")
})

test_that("a fit without a maximum is refused with point_estimates_only", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  items <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
             "t4_lozenges", "t5_general_information",
             "t6_paragraph_comprehension", "t7_sentence",
             "t8_word_classification", "t9_word_meaning", "t10_addition",
             "t11_code", "t12_counting_groups_of_dots",
             "t13_straight_and_curved_capitals")
  # lavaan warns that the solution ran away and is improper, which the
  # refusal reports; only the refusal is examined here.
  caught <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) e)
  full <- caught(fa_exploratory(holzinger_swineford, items = items, nfactors = 7))
  point <- caught(fa_exploratory(holzinger_swineford, items = items, nfactors = 7,
                                 point_estimates_only = TRUE))
  expect_s3_class(point, c("dmar_no_maximum", "dmar_efa_no_solution",
                           "dmar_refusal", "error", "condition"), exact = TRUE)
  expect_identical(point$reason, "no maximum")
  expect_identical(point[["variables"]], full[["variables"]])
  expect_equal(point$loglik, full$loglik, tolerance = 1e-8)
  expect_identical(point$lavaan_call, full$lavaan_call)
})

test_that("point_estimates_only refuses a solution whose information cannot be inverted", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[po_items()]
  # An inverse lavaan could not form (solve() failed), or one with a
  # diagonal entry that is not finite and positive, is unusable.
  expect_true(.efa_inverse_unusable(structure("Error in solve", class = "try-error")))
  expect_true(.efa_inverse_unusable(NULL))
  expect_true(.efa_inverse_unusable(diag(c(1, -1))))
  expect_true(.efa_inverse_unusable(diag(c(1, NA))))
  expect_false(.efa_inverse_unusable(diag(2)))
  fit <- lavaan::efa(d, nfactors = 2, rotation = "none", se = "none",
                     output = "lavaan")
  expect_false(.efa_information_singular(fit))
  # A singular information is refused with the class and reason of the
  # full analysis's refusal of unavailable standard errors, which
  # factor_count() reads.
  local_mocked_bindings(.efa_information_singular = function(fit) TRUE)
  e <- tryCatch(fa_exploratory(d, 2, seed = 113, point_estimates_only = TRUE),
                error = function(e) e)
  expect_s3_class(e, c("dmar_efa_no_solution", "dmar_refusal", "error", "condition"), exact = TRUE)
  expect_identical(e$reason, "standard errors unavailable")
  expect_match(conditionMessage(e), "expected information", fixed = TRUE)
  # The full analysis reads its standard errors instead, and ULS from a
  # supplied matrix, which has none, is not checked in either mode.
  expect_s3_class(fa_exploratory(d, 2, seed = 113), "dmar_fa_exploratory")
  expect_s3_class(fa_exploratory(S = cor(d), N = nrow(d), nfactors = 2,
                                 estimator = "ULS", seed = 113,
                                 point_estimates_only = TRUE),
                  "dmar_fa_exploratory")
})

test_that("point_estimates_only on a FIML fit without unrestricted (H1) information", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[po_items()]
  n <- nrow(d)
  # The H1 warning is kept and the other warnings of these thin data
  # (lavaan's on pairwise coverage, among them) are muffled.
  run <- function(y, po) {
    h1 <- NULL
    x <- withCallingHandlers(
      fa_exploratory(y, 1, point_estimates_only = po),
      dmar_h1_unavailable = function(w) {
        h1 <<- w
        invokeRestart("muffleWarning")
      },
      warning = function(w) invokeRestart("muffleWarning"))
    list(x = x, h1 = h1)
  }
  # Two blocks of three items observed together in three rows only.
  y <- d
  h <- (n - 3L) %/% 2L
  y[seq_len(h), 4:6] <- NA
  y[h + seq_len(n - 3L - h), 1:3] <- NA
  full <- run(y, FALSE)
  point <- run(y, TRUE)
  expect_s3_class(point$h1, "dmar_h1_unavailable")
  expect_match(conditionMessage(full$h1), "H0-based uncertainty", fixed = TRUE)
  # No uncertainty was computed, and the note says so.
  note <- attr(point$x, "h1_diagnostic")
  expect_identical(conditionMessage(point$h1), note)
  expect_true(startsWith(note, paste0("The constrained model's estimates (no ",
    "uncertainty was computed: point_estimates_only = TRUE) are retained")))
  expect_false(grepl("uncertainty are", note, fixed = TRUE))
  # The fit measures are withheld in both modes, and the estimates agree.
  expect_null(attr(full$x, "fit_measures"))
  expect_null(attr(point$x, "fit_measures"))
  expect_equal(point$x$estimate, full$x$estimate, tolerance = 1e-10)
  inference <- c("se", "statistic", "p_value", "ci_lower", "ci_upper")
  expect_true(all(is.na(as.matrix(as.data.frame(point$x)[inference]))))
  out <- capture.output(print(point$x))
  expect_true(any(grepl("Fit measures: not reported", out, fixed = TRUE)))
  expect_false(any(grepl("Confidence level", out, fixed = TRUE)))
  # No row observes both blocks. lavaan 0.7-3 starts one block's loadings
  # at zero, a stationary point at which the expected information is
  # singular and the observed one, which the full analysis's standard
  # errors use, is not; so neither mode refuses the fit, and the
  # estimates agree.
  y <- d
  h <- n %/% 2L
  y[seq_len(h), 4:6] <- NA
  y[h + seq_len(n - h), 1:3] <- NA
  full <- run(y, FALSE)
  point <- run(y, TRUE)
  expect_equal(point$x$estimate, full$x$estimate, tolerance = 1e-10)
  expect_s3_class(point$h1, "dmar_h1_unavailable")
})

test_that("the RMSEA interval is at rmsea_conf_level, apart from the intervals of the solution", {
  skip_if_not_installed("lavaan")
  # Every RMSEA interval defaults to 90 percent, whose limits decide the
  # tests of close and not close fit at .05; the loadings, communalities,
  # and factor correlations keep conf_level, 0.95 by default, and neither
  # level reaches the other.
  d <- holzinger_swineford[po_items()]
  rows <- c("rmsea_ci_lower", "rmsea_ci_upper")
  ci <- function(x) attr(x, "fit_measures")$value[
    match(rows, attr(x, "fit_measures")$term)]
  x <- fa_exploratory(d, 2, seed = 113)
  y <- fa_exploratory(d, 2, seed = 113, rmsea_conf_level = 0.95)
  z <- fa_exploratory(d, 2, seed = 113, conf_level = 0.99)
  expect_equal(eval(formals(fa_exploratory)$rmsea_conf_level), 0.90)
  expect_identical(formals(efa_k)$rmsea_conf_level, 0.90)
  expect_equal(attr(attr(x, "fit_measures"), "conf_level"), 0.90)
  expect_equal(attr(attr(y, "fit_measures"), "conf_level"), 0.95)
  sf <- sem_fit_measures(attr(fa_exploratory(d, 2, rotation = "none"), "fit"))
  expect_equal(ci(x), sf$value[match(rows, sf$term)], tolerance = 1e-10)
  expect_gt(ci(y)[2], ci(x)[2])
  expect_identical(ci(z), ci(x))
  expect_identical(y$ci_lower, x$ci_lower)
  expect_false(isTRUE(all.equal(z$ci_lower, x$ci_lower)))
  expect_identical(attr(x, "rmsea_conf_level"), 0.90)
  expect_identical(dmar_metadata(y)[c("conf_level", "rmsea_conf_level")],
                   list(conf_level = 0.95, rmsea_conf_level = 0.95))
  # The fit line names the RMSEA level, the line beneath the table the
  # level of the intervals in it.
  out <- capture.output(print(x))
  expect_true(any(grepl("^Fit: .*; RMSEA with 90% interval [0-9.]+ \\[", out)))
  expect_true("Confidence level: 95%" %in% out)
  expect_true(any(grepl("; RMSEA with 95% interval ", capture.output(print(y)),
                        fixed = TRUE)))
  expect_error(fa_exploratory(d, 2, rmsea_conf_level = 1.2),
               "'rmsea_conf_level' must be one number strictly between 0 and 1.",
               fixed = TRUE)
})
