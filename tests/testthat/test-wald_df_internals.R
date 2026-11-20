## The t reference distribution of a Wald interval: the wald_df helpers in
## R/missing_information_internals.R, shared by reliability_omega(),
## reliability_alpha(), fa_confirmatory(), mlmr(), and mlmr_mv().

test_that("the critical value and p-value at Inf are the normal ones to the last bit", {
  # R evaluates qt() and pt() at infinite degrees of freedom by calling
  # qnorm() and pnorm(), which is what lets wald_df = Inf reproduce every
  # interval and p-value of DMAR 1.0.x exactly.
  for (cl in c(0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.999)) {
    expect_identical(.wald_crit(cl), stats::qnorm(1 - (1 - cl) / 2))
    expect_identical(.wald_crit(cl, Inf), stats::qnorm(1 - (1 - cl) / 2))
  }
  z <- c(b0 = -40, b_1 = -8.3, b_2 = -1.96, b_3 = 0, b_4 = 0.5,
         b_5 = 15.49, b_6 = NA, b_7 = Inf)
  expect_identical(.wald_p(z), 2 * stats::pnorm(-abs(z)))
  expect_identical(.wald_p(z, Inf), 2 * stats::pnorm(-abs(z)))
  zm <- matrix(c(-2, 0.3, 5, NA), 2, dimnames = list(c("a", "b"),
                                                     c("y1", "y2")))
  expect_identical(.wald_p(zm, Inf), 2 * stats::pnorm(-abs(zm)))
})

test_that("a finite wald_df gives the t quantile and p-value, vectorized", {
  expect_identical(.wald_crit(0.95, 10), stats::qt(0.975, 10))
  expect_identical(.wald_crit(0.90, c(5, 50)), stats::qt(0.95, c(5, 50)))
  expect_identical(.wald_p(c(-2, 2), 7), 2 * stats::pt(-2, 7) * c(1, 1))
  expect_true(.wald_crit(0.95, 10) > .wald_crit(0.95))
  expect_true(.wald_p(2, 10) > .wald_p(2))
})

test_that("wald_df is a single positive number, or the keyword where a rule exists", {
  skip_on_cran()
  for (fn in c("reliability_omega", "reliability_alpha", "fa_confirmatory", "mlmr",
               "mlmr_mv")) {
    expect_identical(.wald_df_check(Inf, fn), Inf)
    expect_identical(.wald_df_check(25, fn), 25)
    expect_identical(.wald_df_check(25L, fn), 25)
    expect_identical(.wald_df_check(10.5, fn), 10.5)
    for (bad in list(0, -3, NA_real_, c(10, 20), TRUE, "inf", NULL,
                     numeric(0))) {
      expect_error(.wald_df_check(bad, fn), "single positive number")
    }
  }
  for (fn in c("mlmr", "mlmr_mv", "fa_confirmatory", "reliability_omega")) {
    expect_identical(.wald_df_check("missing_information", fn),
                     "missing_information")
  }
  # Alpha has no rule: no variant passed the calibration's decision rule,
  # and the error gives its basis as it is. The t reference came barely
  # closer to .95 on average and exceeded .96 in one cell, the cell where
  # the study's saturated fit converged in 53 percent of the replications
  # (the message said it gave "no gain", which criterion 1 contradicts).
  msg <- tryCatch(.wald_df_check("missing_information", "reliability_alpha"),
                  error = conditionMessage)
  expect_match(msg, paste0("not available for reliability_alpha\\(\\): in ",
                           "the calibration study behind the keyword no rule ",
                           "passed for alpha"))
  expect_match(msg, "barely closer to \\.95.*\\.0042 against \\.0043")
  expect_match(msg, paste0("\\(12 items, N = 50, 30 percent missing ",
                           "completely at random, where the study's fit of ",
                           "the saturated model converged in 53 percent of ",
                           "the replications\\) its coverage reached \\.966"))
  expect_false(grepl("no gain", msg))
  # A function the study never saw gets the general refusal.
  expect_error(.wald_df_check("missing_information", "reliability_kr20"),
               "no rule for its degrees of freedom has been calibrated")
  # The error names the keyword only where it is accepted.
  expect_error(.wald_df_check(-1, "mlmr"), "or \"missing_information\"")
  expect_error(.wald_df_check(-1, "fa_confirmatory"), "or \"missing_information\"")
  msg <- tryCatch(.wald_df_check(-1, "reliability_alpha"),
                  error = conditionMessage)
  expect_false(grepl("missing_information", msg))
})

test_that("the keyword is granted in one place, .wald_df_rule()", {
  for (fn in c("mlmr", "mlmr_mv")) {
    rule <- .wald_df_rule(fn)
    expect_identical(rule$method, "effective_n")
    rows <- data.frame(lhs = c("y", "y"), op = c("~1", "~"),
                       rhs = c("", "x_1"))
    expect_identical(rule$k(rows, list(n_coef = 2L)), c(2L, 2L))
    # Every coefficient is covered, and no calibration scope applies.
    expect_null(rule$rows)
    expect_null(rule$calibrated)
  }
  # The calibration study's rule: von Hippel's effective sample size rule
  # with k = 1, for the loadings and the omega rows only.
  rows <- data.frame(lhs = c("f", "f", "f", "y1", "y1", "omega_f",
                             "ave_f", "H_f", "f"),
                     op = c("=~", "=~", "~~", "~~", "~1", ":=", ":=", ":=",
                            "~~"),
                     rhs = c("y1", "y2", "f", "y1", "", "", "", "", "g"))
  for (fn in c("fa_confirmatory", "reliability_omega")) {
    rule <- .wald_df_rule(fn)
    expect_identical(rule$method, "effective_n")
    expect_identical(rule$k(rows, list()), rep(1, nrow(rows)))
    expect_identical(rule$rows(rows, list(omega = "omega_f")),
                     c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,
                       FALSE))
    expect_true(rule$calibrated)
  }
  for (fn in c("reliability_alpha", "reliability_kr20")) {
    expect_null(.wald_df_rule(fn))
  }
})

test_that("a rule that covers some rows leaves the others on the normal reference", {
  # With no fit (a closed form from complete data) every FMI is 0, so the
  # covered rows get N - 1 exactly and the others Inf; the table holds the
  # covered rows only.
  rows <- data.frame(lhs = c("f", "f", "y1", "omega_f"),
                     op = c("=~", "=~", "~~", ":="),
                     rhs = c("y1", "y2", "y1", ""))
  r <- .wald_df_resolve("missing_information", "fa_confirmatory", fit = NULL,
                        rows = rows,
                        context = list(N = 40, omega = "omega_f",
                                       terms = c("l1", "l2", "psi1",
                                                 "omega_f")))
  expect_identical(r$df, c(l1 = 39, l2 = 39, psi1 = Inf, omega_f = 39))
  expect_identical(r$table$term, c("l1", "l2", "omega_f"))
  expect_identical(r$table$fmi, c(0, 0, 0))
  expect_identical(r$table$k, c(1, 1, 1))
  # A rule that covers no row returns no table.
  none <- .wald_df_resolve("missing_information", "fa_confirmatory", fit = NULL,
                           rows = rows[3, ],
                           context = list(N = 40, terms = "psi1"))
  expect_identical(none$df, c(psi1 = Inf))
  expect_null(none$table)
})

test_that("a t reference with no Wald interval to change is refused", {
  expect_invisible(.wald_df_require_wald(Inf, FALSE, "x."))
  expect_invisible(.wald_df_require_wald(10, TRUE, "x."))
  expect_error(.wald_df_require_wald(10, FALSE, "and the call has none."),
               paste0("'wald_df' sets the reference distribution of a ",
                      "Wald interval, and the call has none. Leave ",
                      "wald_df at Inf"))
  expect_error(.wald_df_require_wald("missing_information", FALSE, "y."),
               "Leave wald_df at Inf")
})

test_that("a number is resolved once per parameter, named by term", {
  rows <- data.frame(lhs = c("f", "f"), op = c("=~", "=~"),
                     rhs = c("a", "b"))
  r <- .wald_df_resolve(12, "fa_confirmatory", rows = rows,
                        context = list(terms = c("l1", "l2")))
  expect_identical(r$df, c(l1 = 12, l2 = 12))
  expect_null(r$table)
  r_inf <- .wald_df_resolve(Inf, "reliability_omega")
  expect_identical(unname(r_inf$df), Inf)
})

test_that("the keyword on complete data gives the complete data N - k exactly", {
  # A closed form computed from complete data has no missing information,
  # so the effective sample size rule returns N - k to the last bit.
  rows <- data.frame(lhs = c("y", "y", "y"), op = c("~1", "~", "~"),
                     rhs = c("", "x_1", "x_2"))
  r <- .wald_df_resolve("missing_information", "mlmr", fit = NULL,
                        rows = rows,
                        context = list(N = 40, n_coef = 3L,
                                       terms = c("(Intercept)", "a", "b")))
  expect_identical(r$df, c("(Intercept)" = 37, a = 37, b = 37))
  expect_s3_class(r$table, "dmar_df_missing")
  expect_identical(r$table$fmi, c(0, 0, 0))
  expect_identical(r$table$k, c(3, 3, 3))
  expect_error(.wald_df_resolve("missing_information", "reliability_alpha",
                                rows = rows),
               "has no rule for reliability_alpha")
})

test_that("the recomputed table keeps every row without a standard error", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  tbl <- data.frame(estimate = c(1, 2, 3, 0.5),
                    se = c(0.2, NA, 0, 0.1),
                    test_statistic = c(5, NA, NA, 5),
                    p_value = c(1e-6, NA, NA, 1e-6),
                    ci_lower = c(0.6, NA, 3, 0.3),
                    ci_upper = c(1.4, NA, 3, 0.7))
  out <- .wald_df_apply(tbl, 8, 0.95)
  crit <- stats::qt(0.975, 8)
  expect_equal(out$ci_lower[c(1, 4)], c(1, 0.5) - crit * c(0.2, 0.1),
               tolerance = 1e-15)
  expect_equal(out$ci_upper[c(1, 4)], c(1, 0.5) + crit * c(0.2, 0.1),
               tolerance = 1e-15)
  expect_identical(out$p_value[c(1, 4)], 2 * stats::pt(-5, 8) * c(1, 1))
  # A fixed parameter (se 0, no z ratio) keeps its point interval, and a
  # row with no standard error keeps its NAs.
  expect_identical(out[2:3, ], tbl[2:3, ])
  # One value per row is honored.
  out2 <- .wald_df_apply(tbl, c(8, 8, 8, 30), 0.95)
  expect_equal(out2$ci_upper[4], 0.5 + stats::qt(0.975, 30) * 0.1,
               tolerance = 1e-15)
  # A row on Inf keeps what it had to the last bit, whatever the arithmetic
  # would give: the keyword's rows outside its rule stay lavaan's own.
  tbl$ci_lower[1] <- 0.6000000000000001
  out3 <- .wald_df_apply(tbl, c(Inf, Inf, Inf, 30), 0.95)
  expect_identical(out3[1:3, ], tbl[1:3, ])
  expect_identical(out3$p_value[4], 2 * stats::pt(-5, 30))
})

test_that("the footer line appears only for a finite wald_df", {
  long <- data.frame(term = c("estimate", "se"), value = c(0.8, 0.02))
  attr(long, "wald_df") <- Inf
  expect_null(.wald_df_footer_line(long))
  attr(long, "wald_df") <- NULL
  expect_null(.wald_df_footer_line(long))
  attr(long, "wald_df") <- 30
  expect_identical(.wald_df_footer_line(long), paste0(
    "Wald interval: t reference on 30 degrees of freedom (wald_df), in ",
    "place of the normal"))
  wide <- data.frame(term = "a", estimate = 1, p_value = 0.01)
  attr(wide, "wald_df") <- 136.52114
  expect_match(.wald_df_footer_line(wide),
               "^Wald intervals and p-values: t reference on 136.5 ")
  attr(wide, "wald_df") <- c(118.0076, 136.5211)
  expect_match(.wald_df_footer_line(wide), "118 to 136.5 degrees of freedom")
  # .dmar_tbl_footer() carries the line, so a knitted report shows it too.
  tbl <- .as_dmar_tbl(long)
  expect_identical(.dmar_tbl_footer(tbl)[["wald_df"]],
                   .wald_df_footer_line(long))
  expect_true(any(grepl("t reference on 30 degrees",
                        capture.output(print(tbl)))))
})

test_that("under the keyword the line names the rule, N, the FMI, and the rows", {
  # A reliability table: one coefficient.
  long <- .as_dmar_tbl(data.frame(term = c("estimate", "se"),
                                  value = c(0.8, 0.02)))
  tab <- df_missing(c(omega = 0.2836), N = 301, k = 1)
  attr(long, "wald_df") <- tab$df
  attr(long, "wald_df_table") <- tab
  expect_identical(.wald_df_footer_line(long), paste0(
    "Wald interval: t reference on 214.6 degrees of freedom, from von ",
    "Hippel's (2016) max(3, N(1 - FMI) - 1) with N = 301 and FMI = 0.2836 ",
    "(wald_df = \"missing_information\"), in place of the normal"))
  expect_identical(.dmar_tbl_footer(long)[["wald_df"]],
                   .wald_df_footer_line(long))
  # A fa_confirmatory() table: the loadings and omega took the t reference, by
  # parameter, and the other rows did not.
  wide <- data.frame(term = c("lambda_f_1", "lambda_f_2", "psi_f_1",
                              "omega_f"),
                     estimate = 1, p_value = 0.01)
  tab <- df_missing(c(lambda_f_1 = 0.1, lambda_f_2 = 0.5, omega_f = 0),
                    N = 101, k = 1)
  x <- .wald_df_attach(wide, "missing_information", tab,
                       c("loadings", "omega"))
  expect_identical(attr(x, "wald_df"),
                   c(lambda_f_1 = 89.9, lambda_f_2 = 49.5, omega_f = 100))
  expect_identical(.wald_df_footer_line(x), paste0(
    "Wald intervals and p-values of the loadings and omega: t reference on ",
    "49.5 to 100 degrees of freedom, by parameter, from von Hippel's (2016) ",
    "max(3, N(1 - FMI) - 1) with N = 101 (wald_df = ",
    "\"missing_information\"); the other rows keep the normal reference"))
  # Only the parameters the table reports, once each (an equality
  # constraint repeats a label), and the measurement table's omega alone.
  dup <- df_missing(c(lambda_f = 0.1, lambda_f = 0.1, omega_f = 0),
                    N = 101, k = 1)
  m <- .wald_df_attach(data.frame(term = c("omega_f", "ave_f"),
                                  estimate = 1, p_value = 0.01),
                       "missing_information", dup, "omega")
  expect_identical(attr(m, "wald_df"), c(omega_f = 100))
  expect_match(.wald_df_footer_line(m),
               "^Wald intervals and p-values of omega: t reference on 100 ")
  expect_match(.wald_df_footer_line(m), "N = 101 and FMI = 0 ")
  w2 <- .wald_df_attach(wide, "missing_information", dup)
  expect_identical(names(attr(w2, "wald_df")), "omega_f")
  # A number is recorded as it was, with no table.
  n <- .wald_df_attach(wide, 25)
  expect_identical(attr(n, "wald_df"), 25)
  expect_null(attr(n, "wald_df_table"))
})

test_that("the keyword is refused, with its reason, outside the calibration", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  omega <- function(...) {
    spec <- utils::modifyList(list(denominator = "model_implied",
                                   ci_method = "ml", missing = "fiml"),
                              list(...))
    .wald_df_scope("missing_information", "reliability_omega", spec)
  }
  for (m in c("ml", "mlr")) {
    for (miss in c("listwise", "fiml")) {
      expect_invisible(omega(ci_method = m, missing = miss))
    }
  }
  expect_error(omega(denominator = "observed"),
               "robust omega \\(denominator = \"observed\"\\) has a bootstrap")
  for (m in c("ml_logistic", "mlr_logistic")) {
    expect_error(omega(ci_method = m),
                 paste0("ci_method = \"", m, "\": in the calibration study ",
                        "the t reference on the logit scale covered up to ",
                        "\\.9645"))
  }
  expect_error(omega(ci_method = "adf_logistic"),
               "not the asymptotically distribution free one")
  for (m in c("bca", "percentile", "likelihood", "bonett", "none")) {
    expect_error(omega(ci_method = m),
                 paste0("ci_method = \"", m, "\" is not one; the keyword ",
                        "applies to ci_method = \"ml\" and \"mlr\""))
  }
  expect_error(omega(missing = "two_stage"),
               "DMAR does not compute for a two-stage fit")

  cfa <- function(...) {
    spec <- utils::modifyList(list(output = "verbose", se = "standard",
                                   estimator = "ML", missing = "ml",
                                   two_stage = FALSE, complete = FALSE),
                              list(...))
    .wald_df_scope("missing_information", "fa_confirmatory", spec)
  }
  expect_invisible(cfa())
  expect_invisible(cfa(estimator = "MLR", se = "robust.huber.white"))
  expect_error(cfa(output = "measurement"), "which that table does not show")
  expect_invisible(cfa(estimator = "MLR", se = "robust", missing = "fiml"))
  expect_invisible(cfa(se = "robust.sem", missing = "listwise",
                       complete = TRUE))
  expect_error(cfa(output = "fit"), "output = \"fit\" returns lavaan's own")
  expect_error(cfa(output = "standardized"),
               "evaluated the unstandardized loadings")
  expect_error(cfa(se = "none"), "gives no standard errors")
  expect_error(cfa(se = "bootstrap"), "percentile bootstrap intervals")
  expect_error(cfa(ordered = TRUE), "with ordered items")
  expect_error(cfa(ordered = "y1"), "with ordered items")
  for (e in c("WLS", "WLSMV", "GLS")) {
    expect_error(cfa(estimator = e),
                 paste0("estimator = \"", e, "\": the rule takes"))
  }
  expect_error(cfa(missing = "pairwise"), "a fit to pairwise present data")
  expect_error(cfa(missing = "two_stage", two_stage = TRUE),
               "DMAR does not compute for a two-stage fit")
  # lavaan computes robust.sem for complete data only and falls back to
  # listwise deletion under FIML; the refusal says so for either spelling
  # of FIML.
  for (m in c("ml", "fiml")) {
    expect_error(cfa(se = "robust.sem", missing = m),
                 paste0("se = \"robust.sem\" and missing = \"", m, "\": ",
                        "lavaan does not compute the robust.sem sandwich ",
                        "under full information maximum likelihood but ",
                        "falls back to listwise deletion"))
  }
  expect_invisible(cfa(se = "robust.sem", complete = TRUE))
  expect_error(cfa(se = "first.order"),
               "se = \"first.order\": the calibration study used")
  # A number is never checked against the calibration.
  expect_invisible(.wald_df_scope(20, "fa_confirmatory", list(output = "fit")))
})

test_that("the summary lines name the rule and each coefficient", {
  expect_identical(.wald_df_summary_lines(Inf), character(0))
  expect_identical(.wald_df_summary_lines(NULL), character(0))
  expect_match(.wald_df_summary_lines(20),
               "t reference on 20 degrees of freedom")
  tab <- df_missing(c("(Intercept)" = 0.5, slope = 1e-14, rare = 0.99),
                    N = 50, k = 2)
  lines <- .wald_df_summary_lines(setNames(tab$df, tab$term), tab)
  expect_true(any(grepl("von Hippel's \\(2016\\)", lines)))
  expect_true(any(grepl("N = 50 and k = 2", paste(lines, collapse = " "))))
  # An FMI that differs from zero only by rounding shows as 0, and a
  # floored value says so.
  slope_line <- lines[grepl("slope", lines)]
  expect_match(slope_line, "slope +0[.0]* ")
  expect_match(lines[grepl("rare", lines)], "\\(floor\\)")
})
