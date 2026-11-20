## df_missing() -- von Hippel's (2016) small-sample degrees of freedom for
## FIML Wald intervals, from the fraction of missing information.

test_that("df_missing() reproduces von Hippel's (2016, p. 432) worked example", {
  # A simple regression with n = 25, k = 2, and gamma = .5: the general
  # rule (Eq. 32) on df_complete = n - k = 23 gives 23(.5)(24/26) = 10.6, and
  # the effective sample size rule (Eq. 31) gives 25(.5) - 2 = 10.5.
  general <- df_missing(0.5, df_complete = 23, method = "complete_df")
  expect_equal(general$df, 23 * 0.5 * 24 / 26, tolerance = 1e-12)
  expect_equal(round(general$df, 1), 10.6)
  effective <- df_missing(0.5, N = 25, k = 2)
  expect_equal(effective$df, 10.5, tolerance = 1e-12)
  expect_identical(attr(effective, "method"), "effective_n")
  expect_identical(attr(general, "method"), "complete_df")
})

test_that("df_missing() returns one fixed shape whatever the method", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  cols <- c("term", "df", "df_unbounded", "fmi", "N", "k", "df_complete",
            "at_floor")
  a <- df_missing(c(0.2, 0.4), N = 50, k = 3)
  b <- df_missing(c(0.2, 0.4), df_complete = 47, method = "complete_df")
  for (res in list(a, b)) {
    expect_named(res, cols)
    expect_s3_class(res, c("dmar_df_missing", "dmar_tbl", "data.frame"),
                    exact = TRUE)
    expect_type(res$df, "double")
    expect_type(res$at_floor, "logical")
    expect_identical(attr(res, "floor"), 3)
  }
  # A column the method does not use is NA.
  expect_true(all(is.na(a$df_complete)))
  expect_true(all(is.na(b$k)))
  expect_true(all(is.na(b$N)))
  expect_identical(a$N, c(50, 50))
  expect_identical(b$df_complete, c(47, 47))
})

test_that("with no missing information each rule reduces as documented", {
  # The effective sample size rule is the complete data N - k exactly.
  expect_identical(df_missing(0, N = 40, k = 3)$df, 37)
  # The general rule keeps Barnard and Rubin's observed data factor.
  nu <- 23
  expect_equal(df_missing(0, df_complete = nu, method = "complete_df")$df,
               nu * (nu + 1) / (nu + 3), tolerance = 1e-12)
  # An infinite df_complete is the normal reference.
  expect_identical(
    df_missing(0.3, df_complete = Inf, method = "complete_df")$df, Inf)
})

test_that("the floor bounds the degrees of freedom below and says so", {
  # fmi = .95, N = 25, k = 2: 25(.05) - 2 = -0.75.
  res <- df_missing(0.95, N = 25, k = 2)
  expect_equal(res$df_unbounded, -0.75, tolerance = 1e-12)
  expect_identical(res$df, 3)
  expect_true(res$at_floor)
  expect_identical(df_missing(0.95, N = 25, k = 2, floor = 2)$df, 2)
  expect_identical(df_missing(0.95, N = 25, k = 2, floor = 0)$df, 0)
  # The general rule near an FMI of 1.
  near_one <- df_missing(0.99, df_complete = 23, method = "complete_df")
  expect_equal(near_one$df_unbounded, 23 * 0.01 * 24 / 26,
               tolerance = 1e-12)
  expect_identical(near_one$df, 3)
  expect_true(near_one$at_floor)
  # floor = 0 leaves the general rule unbounded.
  expect_equal(
    df_missing(0.99, df_complete = 23, method = "complete_df", floor = 0)$df,
    23 * 0.01 * 24 / 26, tolerance = 1e-12)
  # A value above the floor is untouched, and a value equal to it is not
  # marked.
  expect_false(df_missing(0.5, N = 10, k = 2)$at_floor)
  expect_identical(df_missing(0.5, N = 10, k = 2)$df, 3)
})

test_that("df_missing() is vectorized over parameters, N, k, and df_complete", {
  g <- c(slope = 0.5, intercept = 0.1, complete = 0)
  res <- df_missing(g, N = 25, k = 2)
  expect_identical(res$term, c("slope", "intercept", "complete"))
  expect_equal(res$df, 25 * (1 - g) - 2, tolerance = 1e-12,
               ignore_attr = TRUE)
  # Per-row N and k, as a multigroup or multivariate fit would need.
  res2 <- df_missing(c(0.2, 0.3), N = c(156, 145), k = c(2, 3))
  expect_equal(res2$df, c(156 * 0.8 - 2, 145 * 0.7 - 3), tolerance = 1e-12)
  # Per-row df_complete.
  res3 <- df_missing(c(0.2, 0.3), df_complete = c(30, Inf),
                     method = "complete_df")
  expect_equal(res3$df[1], 30 * 0.8 * 31 / 33, tolerance = 1e-12)
  expect_identical(res3$df[2], Inf)
  # Unnamed input is labeled parameter_1, parameter_2, ...
  expect_identical(df_missing(c(0.1, 0.2), N = 30, k = 1)$term,
                   c("parameter_1", "parameter_2"))
  expect_error(df_missing(c(0.1, 0.2, 0.3), N = c(30, 40), k = 1),
               "'N' must have length 1 or the length")
})

test_that("a missing FMI gives missing degrees of freedom for that row only", {
  res <- df_missing(c(0.2, NA), N = 30, k = 2)
  expect_equal(res$df[1], 30 * 0.8 - 2, tolerance = 1e-12)
  expect_true(is.na(res$df[2]))
  expect_true(is.na(res$at_floor[2]))
})

test_that("df_missing() has no default for k or df_complete and refuses unused arguments", {
  expect_error(df_missing(0.3, N = 50), "needs 'k'")
  expect_error(df_missing(0.3, k = 2), "needs 'N'")
  expect_error(df_missing(0.3, method = "complete_df"), "needs 'df_complete'")
  expect_error(df_missing(0.3, N = 50, k = 2, df_complete = 48),
               "'df_complete' is used by method = \"complete_df\"")
  expect_error(df_missing(0.3, df_complete = 48, k = 2, method = "complete_df"),
               "'k' is used by method = \"effective_n\"")
  expect_error(df_missing(0.3, df_complete = 48, N = 50, method = "complete_df"),
               "'N' is used by method = \"effective_n\"")
})

test_that("df_missing() validates its inputs", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  expect_error(df_missing(0.3, N = 50, k = 1.5), "whole number")
  expect_error(df_missing(0.3, N = 50, k = -1), "whole number")
  expect_error(df_missing(0.3, N = -5, k = 1), "finite and positive")
  expect_error(df_missing(0.3, N = 2, k = 2), "'N' must exceed 'k'")
  expect_error(df_missing(0.3, df_complete = 0, method = "complete_df"),
               "'df_complete' must be positive")
  expect_error(df_missing(0.3, N = 50, k = 1, floor = -1), "'floor'")
  expect_error(df_missing(0.3, N = 50, k = 1, floor = c(2, 3)), "'floor'")
  expect_error(df_missing("a", N = 50, k = 1), "must be numeric")
  expect_error(df_missing(numeric(0), N = 50, k = 1), "has no values")
})

test_that("an FMI at or above 1, or clearly negative, is refused", {
  expect_error(df_missing(1, N = 50, k = 1), "must be below 1")
  expect_error(df_missing(1.2, N = 50, k = 1), "must be below 1")
  # A negative FMI from the structured or Hessian estimators.
  expect_error(df_missing(-0.01, N = 50, k = 1),
               "information = \"unstructured\"")
  # Rounding error around a complete data parameter is set to 0.
  res <- df_missing(-1e-12, N = 50, k = 1)
  expect_identical(res$fmi, 0)
  expect_identical(res$df, 49)
})

test_that("a table with fmi, N, and term columns supplies all three", {
  # The shape fraction_missing_information() returns: one row per
  # parameter, with the N basis of each row.
  tab <- data.frame(term = c("visual =~ x2", "visual =~ x2"),
                    op = "=~", group = c("Pasteur", "Grant-White"),
                    fmi = c(0.3, 0.4), N = c(156, 145),
                    stringsAsFactors = FALSE)
  class(tab) <- c("dmar_fmi", "dmar_tbl", "data.frame")
  res <- df_missing(tab, k = 1)
  expect_identical(res$term, tab$term)
  expect_identical(res$N, c(156, 145))
  expect_equal(res$df, c(156 * 0.7 - 1, 145 * 0.6 - 1), tolerance = 1e-12)
  expect_error(df_missing(tab, N = 301, k = 1), "leave 'N' unset")
  # The general rule echoes the table's N without using it.
  res_cd <- df_missing(tab, df_complete = 100, method = "complete_df")
  expect_identical(res_cd$N, c(156, 145))
  expect_equal(res_cd$df, 100 * c(0.7, 0.6) * 101 / 103, tolerance = 1e-12)
  # A table without an fmi column is refused, and one without N needs it.
  expect_error(df_missing(data.frame(x = 0.3), N = 50, k = 1),
               "must have an fmi column")
  no_n <- data.frame(fmi = c(0.1, 0.2))
  expect_error(df_missing(no_n, k = 1), "needs 'N'")
  expect_identical(df_missing(no_n, N = 40, k = 1)$term,
                   c("parameter_1", "parameter_2"))
})

test_that("the general rule is Barnard and Rubin's observed data degrees of freedom", {
  # df_barnard_rubin()'s nu_obs is the complete_df rule at the
  # imputations' own gamma.
  br <- df_barnard_rubin(between = c(0.04, 0.01), within = c(0.10, 0.20),
                         m = 5, df_complete = 23)
  cd <- df_missing(br$gamma, df_complete = 23, method = "complete_df")
  expect_equal(cd$df_unbounded, br$nu_obs, tolerance = 1e-12)
})

test_that("tidy() and glance() give one row per parameter and one summary row", {
  res <- df_missing(c(a = 0.2, b = 0.97), N = 30, k = 2)
  td <- generics::tidy(res)
  expect_named(td, c("term", "estimate", "df_unbounded", "fmi", "N", "k",
                     "df_complete", "at_floor"))
  expect_identical(td$term, c("a", "b"))
  expect_identical(td$estimate, res$df)
  gl <- generics::glance(res)
  expect_named(gl, c("n_terms", "n_at_floor", "method", "floor"))
  expect_identical(gl$n_at_floor, 1L)
  expect_identical(gl$method, "effective_n")
})

test_that("the printed table names its rule and floor, at the console and knitted", {
  # The df column was the only record of which rule made it, read from
  # the N and k (or df_complete) columns the method filled. One line under the
  # table now names the rule and the floor, the same line in a knitted
  # report.
  res <- df_missing(c(a = 0.2, b = 0.97), N = 30, k = 2)
  out <- capture.output(print(res))
  rule <- c(paste0("Degrees of freedom: von Hippel's (2016, Eq. 31) ",
                   "effective sample size rule,"),
            "N(1 - fmi) - k, bounded below at 3")
  expect_identical(utils::tail(out, 3L), c("", rule))
  footer <- DMAR:::.dmar_tbl_footer(res)
  expect_identical(unname(footer), paste(rule, collapse = "\n"))
  expect_identical(names(footer), "df_rule")

  cd <- capture.output(print(df_missing(0.5, df_complete = 23,
                                        method = "complete_df",
                                        floor = 0)))
  expect_identical(utils::tail(cd, 2L), c(
    "Degrees of freedom: von Hippel's (2016, Eq. 32) complete data rule,",
    "df_complete (1 - fmi) (df_complete + 1) / (df_complete + 3), bounded below at 0"))

  skip_if_not_installed("knitr")
  knitted <- as.character(knit_print.dmar_tbl(res))
  expect_match(knitted, paste0(
    "|\n\nDegrees of freedom: von Hippel's (2016, Eq. 31) effective ",
    "sample size rule,\nN(1 - fmi) - k, bounded below at 3"),
    fixed = TRUE)
  kab <- as.character(as_kable(res, format = "pipe"))
  expect_identical(kab[1L], paste0(
    "Table: Degrees of freedom: von Hippel's (2016, Eq. 31) effective ",
    "sample size rule, N(1 - fmi) - k, bounded below at 3"))

  # A column subset loses the attributes and prints no line.
  sub <- capture.output(print(res[, c("term", "df")]))
  expect_false(any(grepl("Degrees of freedom", sub, fixed = TRUE)))
})

test_that("df_missing() takes lavaan's FMI for a FIML regression (real data)", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  # t26_flags is observed only for the 145 Grant-White students of 301. In
  # this saturated regression lavaan's default (Hessian) FMI and the
  # unstructured FMI agree to about 1e-9, so the effective sample size
  # degrees of freedom match those computed from the unstructured FMI:
  # 301(1 - .5364747) - 3 = 136.5211 and 301(1 - .5979813) - 3 = 118.0076.
  fit <- lavaan::sem("t26_flags ~ t1_visual_perception + t4_lozenges",
                     data = holzinger_swineford, missing = "ml",
                     fixed.x = FALSE)
  pe <- lavaan::parameterEstimates(fit, fmi = TRUE)
  slopes <- pe[pe$op == "~", ]
  N <- lavaan::lavInspect(fit, "nobs")
  expect_identical(as.numeric(N), 301)
  res <- df_missing(setNames(slopes$fmi, slopes$rhs), N = N, k = 3)
  expect_equal(res$df, c(136.521123, 118.0076156), tolerance = 1e-6)
  expect_equal(res$df, N * (1 - slopes$fmi) - 3, tolerance = 1e-12)
  expect_identical(res$term, c("t1_visual_perception", "t4_lozenges"))
})

test_that("df_missing() takes the table fraction_missing_information() returns", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  # The six spatial tests of the bundled battery: t25 and t26 are missing
  # for the 156 Pasteur students, so their loadings carry about half the
  # information complete data would. The table supplies each row's FMI,
  # its N basis (the 301 cases of the single group), and its term, and k
  # is the caller's.
  fit <- fa_confirmatory(data = holzinger_swineford,
               factors = list(spatial = c("t1_visual_perception",
                                          "t2_cubes", "t3_paper_form_board",
                                          "t4_lozenges",
                                          "t25_paper_form_board_r",
                                          "t26_flags")),
               missing = "fiml", output = "fit")
  fmi <- fraction_missing_information(fit)
  res <- df_missing(fmi, k = 1)
  expect_identical(res$term, fmi$term)
  expect_identical(res$N, fmi$N)
  expect_true(all(res$N == 301))
  # Each FMI is used as estimated; the unstructured estimator cannot be
  # negative beyond rounding, which df_missing() sets to 0.
  expect_identical(res$fmi, pmax(fmi$fmi, 0))
  expect_equal(res$df_unbounded, 301 * (1 - res$fmi) - 1,
               tolerance = 1e-12)
  expect_false(any(res$at_floor))
  # omega_spatial has an FMI of .2442136 (Chen and Savalei's third
  # estimator), so 301(1 - .2442136) - 1 = 226.4917 degrees of freedom; the
  # loading of t25 (FMI .5116) has about half of the complete data 300.
  om <- res[res$term == "omega_spatial", ]
  expect_equal(om$fmi, 0.2442136260, tolerance = 1e-6)
  expect_equal(om$df, 301 * (1 - 0.2442136260) - 1, tolerance = 1e-4)
  l25 <- res[res$term == "spatial =~ t25_paper_form_board_r", ]
  expect_equal(l25$df, 301 * (1 - 0.5116258498) - 1, tolerance = 1e-4)
  # The general rule echoes the table's N without using it.
  gen <- df_missing(fmi, df_complete = 300, method = "complete_df")
  expect_identical(gen$N, fmi$N)
  expect_equal(gen$df, 300 * (1 - res$fmi) * 301 / 303, tolerance = 1e-12)
})

test_that("df_missing() takes the table of fraction_missing_information() for an mlmr() fit", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  skip_if_not_installed("lavaan")
  # The help page's example: the rows of the FMI table carry the N behind
  # each slope, and the degrees of freedom are those computed above from
  # lavaan's FMI for the same regression, since in this saturated model the
  # unstructured and Hessian FMIs agree to about 1e-9.
  fit <- mlmr(t26_flags ~ t1_visual_perception + t4_lozenges,
              data = holzinger_swineford, ci_method = "wald")
  tab <- fraction_missing_information(fit)
  slopes <- tab[tab$op == "~", ]
  res <- df_missing(slopes, k = 3)
  expect_identical(res$term, c("t26_flags ~ t1_visual_perception",
                               "t26_flags ~ t4_lozenges"))
  expect_identical(res$N, c(301, 301))
  expect_equal(res$df, c(136.521123, 118.0076156), tolerance = 1e-6)
  expect_equal(res$df, 301 * (1 - slopes$fmi) - 3, tolerance = 1e-12)
  expect_false(any(res$at_floor))
})
