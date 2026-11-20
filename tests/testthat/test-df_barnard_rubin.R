## df_barnard_rubin() -- Barnard and Rubin's (1999) small-sample degrees of
## freedom for an estimate pooled over m imputations.

test_that("df_barnard_rubin() matches mice's barnard.rubin()", {
  # mice 3.19.0, mice:::barnard.rubin(m = 5, b = 0.04, t = 0.148,
  # dfcom = 23) returns 10.4159384445, and with dfcom = Inf it returns
  # Rubin's (m - 1) / lambda^2 = 38.0277777778. The values are hard-coded
  # because mice is not a dependency of DMAR.
  res <- df_barnard_rubin(between = 0.04, within = 0.10, m = 5,
                          df_complete = 23, floor = 0)
  expect_equal(res$df, 10.4159384445, tolerance = 1e-10)
  expect_equal(res$gamma, 0.3243243243, tolerance = 1e-9)
  expect_equal(res$total_variance, 0.10 + 1.2 * 0.04, tolerance = 1e-14)
  expect_equal(res$relative_increase, 1.2 * 0.04 / 0.10, tolerance = 1e-14)
  expect_equal(res$nu_m, 38.0277777778, tolerance = 1e-10)
  expect_equal(res$nu_obs, 14.3451143451, tolerance = 1e-10)
  rubin <- df_barnard_rubin(between = 0.04, within = 0.10, m = 5,
                            floor = 0)
  expect_equal(rubin$df, 38.0277777778, tolerance = 1e-10)
  expect_identical(rubin$nu_obs, Inf)
})

test_that("the formula is von Hippel's Eqs. 23 to 25 and mice's algebraic form", {
  between <- c(0.04, 0.002, 0.3)
  within  <- c(0.10, 0.05, 0.02)
  m <- 7
  df_complete <- 40
  res <- df_barnard_rubin(between, within, m = m, df_complete = df_complete,
                          floor = 0)
  lambda <- (1 + 1 / m) * between / (within + (1 + 1 / m) * between)
  # mice's form: (m - 1) tmp / ((df_complete + 3)(m - 1) + lambda^2 tmp), with
  # tmp = (1 - lambda)(1 + df_complete) df_complete.
  tmp <- (1 - lambda) * (1 + df_complete) * df_complete
  mice_form <- (m - 1) * tmp / ((df_complete + 3) * (m - 1) + lambda^2 * tmp)
  expect_equal(res$df, mice_form, tolerance = 1e-12)
  expect_equal(res$gamma, lambda, tolerance = 1e-14)
  # Barnard and Rubin's df never exceed the observed data df, which never
  # exceed the complete data df.
  expect_true(all(res$df <= res$nu_obs))
  expect_true(all(res$nu_obs < df_complete))
})

test_that("the default floor of 3 applies and is marked", {
  # A large between-imputation variance with few imputations.
  res <- df_barnard_rubin(between = 0.5, within = 0.05, m = 3,
                          df_complete = 20)
  expect_lt(res$df_unbounded, 3)
  expect_identical(res$df, 3)
  expect_true(res$at_floor)
  expect_identical(
    df_barnard_rubin(between = 0.5, within = 0.05, m = 3, df_complete = 20,
                     floor = 0)$df,
    res$df_unbounded)
  expect_identical(attr(res, "floor"), 3)
})

test_that("the boundary cases follow the formula's limits", {
  # Imputations that agree exactly: gamma = 0, Rubin's df are infinite,
  # and Barnard and Rubin's df are the observed data df alone.
  res <- df_barnard_rubin(between = 0, within = 0.1, m = 5, df_complete = 23)
  expect_identical(res$gamma, 0)
  expect_identical(res$nu_m, Inf)
  expect_equal(res$df, 23 * 24 / 26, tolerance = 1e-12)
  expect_identical(
    df_barnard_rubin(between = 0, within = 0.1, m = 5)$df, Inf)
  # No within-imputation variance: gamma = 1 and no information.
  none <- df_barnard_rubin(between = 0.1, within = 0, m = 5, df_complete = Inf,
                           floor = 0)
  expect_identical(none$gamma, 1)
  expect_identical(none$df, 0)
})

test_that("df_barnard_rubin() is vectorized and labels rows by name", {
  est <- cbind(slope = c(0.42, 0.51, 0.38, 0.47, 0.55),
               intercept = c(1.90, 1.94, 1.88, 1.93, 1.91))
  u <- cbind(slope = c(0.010, 0.011, 0.009, 0.010, 0.012),
             intercept = c(0.020, 0.021, 0.019, 0.020, 0.022))
  res <- df_barnard_rubin(between = apply(est, 2, var),
                          within = colMeans(u), m = 5, df_complete = 23)
  expect_identical(res$term, c("slope", "intercept"))
  one <- df_barnard_rubin(between = var(est[, "slope"]),
                          within = mean(u[, "slope"]), m = 5, df_complete = 23)
  expect_equal(res$df[1], one$df, tolerance = 1e-12)
  # Unnamed input: parameter_1, parameter_2; a per-row df_complete.
  res2 <- df_barnard_rubin(c(0.04, 0.02), c(0.1, 0.1), m = 5,
                           df_complete = c(23, 50))
  expect_identical(res2$term, c("parameter_1", "parameter_2"))
  expect_identical(attr(res2, "df_complete"), c(23, 50))
  # A missing value stays in its row.
  res3 <- df_barnard_rubin(c(0.04, NA), c(0.1, 0.1), m = 5, df_complete = 23)
  expect_true(is.na(res3$df[2]))
  expect_false(is.na(res3$df[1]))
})

test_that("df_barnard_rubin() returns one fixed shape", {
  res <- df_barnard_rubin(0.04, 0.10, m = 5, df_complete = 23)
  expect_named(res, c("term", "df", "df_unbounded", "gamma",
                      "relative_increase", "total_variance", "nu_m",
                      "nu_obs", "at_floor"))
  expect_s3_class(res, c("dmar_df_barnard_rubin", "dmar_tbl", "data.frame"),
                  exact = TRUE)
  expect_identical(attr(res, "m"), 5)
  td <- generics::tidy(res)
  expect_identical(td$estimate, res$df)
  expect_identical(names(td)[1:2], c("term", "estimate"))
  gl <- generics::glance(res)
  expect_named(gl, c("n_terms", "n_at_floor", "m", "floor"))
})

test_that("the printed table names its rule, m, df_complete, and floor", {
  skip_on_cran()
  # One line under the table, the same line in a knitted report.
  res <- df_barnard_rubin(between = 0.04, within = 0.10, m = 5,
                          df_complete = 23)
  rule <- c(paste0("Degrees of freedom: Barnard and Rubin's (1999) rule ",
                   "for m = 5 imputations,"),
            paste0("1 / (1 / nu_m + 1 / nu_obs) with df_complete = 23, ",
                   "bounded below at 3"))
  out <- capture.output(print(res))
  expect_identical(utils::tail(out, 3L), c("", rule))
  footer <- DMAR:::.dmar_tbl_footer(res)
  expect_identical(unname(footer), paste(rule, collapse = "\n"))
  expect_identical(names(footer), "df_rule")

  # With every df_complete infinite the degrees of freedom are Rubin's.
  ru <- capture.output(print(df_barnard_rubin(0.04, 0.10, m = 5)))
  expect_identical(utils::tail(ru, 2L), c(
    paste0("Degrees of freedom: Rubin's (1987) large-sample rule for ",
           "m = 5 imputations,"),
    "(m - 1) / gamma^2 (df_complete = Inf), bounded below at 3"))
  # Several values of df_complete are listed, and the floor is the one used.
  several <- df_barnard_rubin(c(a = 0.04, b = 0.02, c = 0.01),
                              c(0.10, 0.10, 0.10), m = 5,
                              df_complete = c(23, Inf, 40.5), floor = 2.5)
  expect_match(DMAR:::.dmar_tbl_footer(several),
               "with df_complete = 23, Inf, and 40.5, bounded below at 2.5$")
  two <- df_barnard_rubin(c(0.04, 0.02), c(0.10, 0.10), m = 5,
                          df_complete = c(23, Inf), floor = 0)
  expect_match(DMAR:::.dmar_tbl_footer(two),
               "with df_complete = 23 and Inf, bounded below at 0$")

  skip_if_not_installed("knitr")
  kab <- as.character(as_kable(res, format = "pipe"))
  expect_identical(kab[1L], paste("Table:", paste(rule, collapse = " ")))
  expect_match(as.character(knit_print.dmar_tbl(res)),
               paste(rule, collapse = "\n"), fixed = TRUE)

  # A column subset loses the attributes and prints no line.
  sub <- capture.output(print(res[, c("term", "df")]))
  expect_false(any(grepl("Degrees of freedom", sub, fixed = TRUE)))
})

test_that("df_barnard_rubin() validates its inputs", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  expect_error(df_barnard_rubin(0.04, 0.10, m = 1), "'m' must be")
  expect_error(df_barnard_rubin(0.04, 0.10, m = 2.5), "'m' must be")
  expect_error(df_barnard_rubin(-0.04, 0.10, m = 5), "cannot be negative")
  expect_error(df_barnard_rubin(0.04, c(0.1, 0.2), m = 5),
               "one per value of 'between'")
  expect_error(df_barnard_rubin(0, 0, m = 5), "total variance is 0")
  expect_error(df_barnard_rubin(0.04, 0.10, m = 5, df_complete = 0),
               "'df_complete' must be positive")
  expect_error(df_barnard_rubin(c(0.04, 0.02), c(0.1, 0.1), m = 5,
                                df_complete = c(10, 20, 30)),
               "'df_complete' must have length 1")
  expect_error(df_barnard_rubin("a", 0.10, m = 5), "'between' must be")
  expect_error(df_barnard_rubin(0.04, 0.10, m = 5, floor = NA), "'floor'")
})
