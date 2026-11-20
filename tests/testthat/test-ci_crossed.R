## ci_crossed(): a published anchor for the mean squares (Shrout & Fleiss,
## 1979, Table 3); the interval against an independent implementation,
## lmerTest's Satterthwaite interval for the intercept of the crossed
## random effects model, both live (when lmerTest is installed) and as
## pinned values; and a re-derivation from the two-way ANOVA table of
## stats::anova(stats::lm()), which checks the arithmetic but shares the
## method, so it is not an oracle for it.

# Shrout and Fleiss (1979), Table 2: four judges each rated six targets.
shrout_fleiss <- matrix(c( 9, 2, 5, 8,
                           6, 1, 3, 2,
                           8, 4, 6, 8,
                           7, 1, 2, 6,
                          10, 5, 6, 9,
                           6, 2, 4, 7),
                        nrow = 6, byrow = TRUE)

val <- function(x, term) x$value[x$term == term]

# The interval re-derived: the mean squares from lm() on the long form of
# the array, the truncated ANOVA components, the plug-in variance of the
# grand mean, and Satterthwaite's degrees of freedom. The truncation and
# the Satterthwaite formula are the implementation's own, so agreement
# here checks the arithmetic, not the method; the lmerTest tests below
# check the method.
anova_interval <- function(x, conf_level = 0.95) {
  n_u <- nrow(x); n_r <- ncol(x)
  d <- data.frame(y = c(x), unit = factor(rep(seq_len(n_u), n_r)),
                  rater = factor(rep(seq_len(n_r), each = n_u)))
  tab <- stats::anova(stats::lm(y ~ unit + rater, data = d))
  ms_u <- tab["unit", "Mean Sq"]; ms_r <- tab["rater", "Mean Sq"]
  ms_e <- tab["Residuals", "Mean Sq"]
  vu <- max(0, (ms_u - ms_e) / n_r)
  vr <- max(0, (ms_r - ms_e) / n_u)
  V <- vu / n_u + vr / n_r + ms_e / (n_u * n_r)
  # The same V as a combination of mean squares, for the Satterthwaite df.
  a <- c(u = as.numeric(ms_u > ms_e), r = as.numeric(ms_r > ms_e))
  a <- c(a, e = 1 - a[["u"]] - a[["r"]])
  ms <- c(ms_u, ms_r, ms_e)
  dfs <- c(n_u - 1, n_r - 1, (n_u - 1) * (n_r - 1))
  nu <- sum(a * ms)^2 / sum((a * ms)^2 / dfs)
  m <- mean(x)
  q <- qt(1 - (1 - conf_level) / 2, nu)
  list(ms = c(ms_u, ms_r, ms_e), mean = m, V = V, nu = nu,
       lower = m - q * sqrt(V), upper = m + q * sqrt(V),
       vu = vu, vr = vr, ve = ms_e)
}

test_that("the mean squares reproduce Shrout and Fleiss (1979, Table 3)", {
  ms <- DMAR:::.crossed_ms(shrout_fleiss)
  # Published to two decimals: between targets 11.24, between judges 32.49,
  # residual 1.02.
  expect_equal(round(unname(ms[c("ms_unit", "ms_rater", "ms_residual")]), 2),
               c(11.24, 32.49, 1.02))
  expect_equal(unname(ms[c("ms_unit", "ms_rater", "ms_residual")]),
               anova_interval(shrout_fleiss)$ms, tolerance = 1e-12)
  expect_equal(unname(ms["grand_mean"]), 127 / 24)
})

# The Satterthwaite interval for the intercept of
# lmer(y ~ 1 + (1 | unit) + (1 | rater)), fitted by REML. In a balanced
# crossed design whose ANOVA component estimates are positive, REML
# returns those estimates and lmerTest's Satterthwaite degrees of freedom
# for the intercept are those of the same combination of mean squares,
# reached through the asymptotic covariance of the REML variance
# parameters and a numerical gradient rather than through the ANOVA
# table.
lmer_interval <- function(x, conf_level = 0.95) {
  n_u <- nrow(x); n_r <- ncol(x)
  d <- data.frame(y = c(x), unit = factor(rep(seq_len(n_u), n_r)),
                  rater = factor(rep(seq_len(n_r), each = n_u)))
  fit <- lmerTest::lmer(y ~ 1 + (1 | unit) + (1 | rater), data = d)
  co <- summary(fit)$coefficients
  q <- qt(1 - (1 - conf_level) / 2, co[1, "df"])
  vc <- as.data.frame(lme4::VarCorr(fit))
  list(mean = co[1, "Estimate"], se = co[1, "Std. Error"], df = co[1, "df"],
       lower = co[1, "Estimate"] - q * co[1, "Std. Error"],
       upper = co[1, "Estimate"] + q * co[1, "Std. Error"],
       vu = vc$vcov[vc$grp == "unit"], vr = vc$vcov[vc$grp == "rater"],
       ve = vc$vcov[vc$grp == "Residual"])
}

test_that("the interval matches lmerTest's Satterthwaite interval for the intercept", {
  skip_on_cran()
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("lme4")
  res <- ci_crossed(shrout_fleiss)
  ref <- lmer_interval(shrout_fleiss)
  expect_equal(val(res, "grand_mean"), ref$mean, tolerance = 1e-10)
  expect_equal(val(res, "se"), ref$se, tolerance = 1e-5)
  expect_equal(val(res, "df"), ref$df, tolerance = 1e-5)
  expect_equal(val(res, "ci_lower"), ref$lower, tolerance = 1e-5)
  expect_equal(val(res, "ci_upper"), ref$upper, tolerance = 1e-5)
  expect_equal(c(val(res, "var_unit"), val(res, "var_rater"),
                 val(res, "var_residual")),
               c(ref$vu, ref$vr, ref$ve), tolerance = 1e-5)

  # Arrays drawn with both components clearly positive, so that neither
  # estimate is truncated (REML would then set it to zero and part ways
  # with the ANOVA estimate on the degrees of freedom). lme4's optimizer
  # stops at about 1e-5 of the optimum, hence the looser tolerance.
  set.seed(20261105)
  for (i in 1:3) {
    x <- 50 + outer(rnorm(15, 0, 1), rnorm(5, 0, 0.8), "+") +
      matrix(rnorm(75, 0, 0.6), 15, 5)
    res <- ci_crossed(x)
    expect_gt(val(res, "var_unit"), 0)
    expect_gt(val(res, "var_rater"), 0)
    ref <- lmer_interval(x)
    expect_equal(val(res, "df"), ref$df, tolerance = 1e-4, info = i)
    expect_equal(val(res, "se"), ref$se, tolerance = 1e-4, info = i)
    expect_equal(c(val(res, "ci_lower"), val(res, "ci_upper")),
                 c(ref$lower, ref$upper), tolerance = 1e-4, info = i)
  }
})

test_that("the interval reproduces values pinned from lmerTest", {
  # lmerTest 3.2.1 with lme4 2.0.6 on the Shrout and Fleiss data gives an
  # intercept of 5.291667 with standard error 1.333985 on 4.836569
  # Satterthwaite degrees of freedom, and the 95% limits 1.827402 and
  # 8.755932. Pinned here so the check runs where lmerTest is absent.
  res <- ci_crossed(shrout_fleiss)
  expect_equal(val(res, "grand_mean"), 5.291667, tolerance = 1e-6)
  expect_equal(val(res, "se"), 1.333985, tolerance = 1e-5)
  expect_equal(val(res, "df"), 4.836569, tolerance = 1e-5)
  expect_equal(val(res, "ci_lower"), 1.827402, tolerance = 1e-5)
  expect_equal(val(res, "ci_upper"), 8.755932, tolerance = 1e-5)
})

test_that("the interval agrees with a re-derivation from the ANOVA mean squares", {
  ref <- anova_interval(shrout_fleiss)
  res <- ci_crossed(shrout_fleiss)
  expect_equal(val(res, "grand_mean"), ref$mean, tolerance = 1e-12)
  expect_equal(val(res, "se"), sqrt(ref$V), tolerance = 1e-12)
  expect_equal(val(res, "df"), ref$nu, tolerance = 1e-12)
  expect_equal(val(res, "ci_lower"), ref$lower, tolerance = 1e-12)
  expect_equal(val(res, "ci_upper"), ref$upper, tolerance = 1e-12)
  expect_equal(val(res, "var_unit"), ref$vu, tolerance = 1e-12)
  expect_equal(val(res, "var_rater"), ref$vr, tolerance = 1e-12)
  expect_equal(val(res, "var_residual"), ref$ve, tolerance = 1e-12)
  # Neither component is truncated here, so V is (MS_u + MS_r - MS_e) / 24.
  expect_equal(val(res, "se")^2, sum(c(1, 1, -1) * ref$ms) / 24,
               tolerance = 1e-12)
  # At 90% as well.
  ref90 <- anova_interval(shrout_fleiss, conf_level = 0.90)
  res90 <- ci_crossed(shrout_fleiss, conf_level = 0.90)
  expect_equal(c(val(res90, "ci_lower"), val(res90, "ci_upper")),
               c(ref90$lower, ref90$upper), tolerance = 1e-12)
  expect_equal(attr(res90, "conf_level"), 0.90)
})

test_that("a negative component estimate is truncated, and the degrees of freedom follow", {
  set.seed(20261105)
  n_u <- 25; n_r <- 6
  # Raters that hardly differ: MS_rater falls below MS_residual.
  x <- outer(rnorm(n_u, 0, 1), rep(0, n_r), "+") +
    matrix(rnorm(n_u * n_r, 0, 0.7), n_u, n_r)
  x <- sweep(x, 2, colMeans(x) - mean(x))   # equal rater means
  ref <- anova_interval(x)
  expect_lt(ref$ms[2], ref$ms[3])
  res <- ci_crossed(x)
  expect_equal(val(res, "var_rater"), 0)
  expect_equal(val(res, "df"), ref$nu, tolerance = 1e-12)
  # With only the unit component kept, the combination is MS_u alone.
  expect_equal(val(res, "df"), n_u - 1, tolerance = 1e-12)
  expect_equal(val(res, "se"), sqrt(ref$ms[1] / (n_u * n_r)), tolerance = 1e-12)

  # Units and raters that both hardly differ: only MS_e is left, on
  # (n_u - 1)(n_r - 1) degrees of freedom, and V-hat is MS_e / (n_u n_r).
  y <- matrix(rnorm(n_u * n_r), n_u, n_r)
  y <- y - outer(rowMeans(y), colMeans(y), "+") + mean(y)
  ref <- anova_interval(y)
  res <- ci_crossed(y)
  expect_equal(c(val(res, "var_unit"), val(res, "var_rater")), c(0, 0))
  expect_equal(val(res, "df"), (n_u - 1) * (n_r - 1), tolerance = 1e-12)
  expect_equal(val(res, "se"), sqrt(ref$ms[3] / (n_u * n_r)), tolerance = 1e-12)
})

test_that("the interval is equivariant under a change of location and scale", {
  a <- ci_crossed(shrout_fleiss)
  b <- ci_crossed(10 + 2.5 * shrout_fleiss)
  expect_equal(val(b, "ci_lower"), 10 + 2.5 * val(a, "ci_lower"))
  expect_equal(val(b, "ci_upper"), 10 + 2.5 * val(a, "ci_upper"))
  expect_equal(val(b, "df"), val(a, "df"))
  expect_equal(val(b, "var_rater"), 2.5^2 * val(a, "var_rater"))
})

test_that("a data frame gives the matrix's result, and the table has the documented rows", {
  res <- ci_crossed(shrout_fleiss)
  expect_equal(ci_crossed(as.data.frame(shrout_fleiss))$value, res$value)
  expect_s3_class(res, c("dmar_tbl", "data.frame"))
  expect_identical(res$term,
                   c("ci_lower", "grand_mean", "ci_upper", "se", "df",
                     "var_unit", "var_rater", "var_residual", "n_units",
                     "n_raters"))
  expect_type(res$value, "double")
  expect_equal(c(val(res, "n_units"), val(res, "n_raters")), c(6, 4))
  td <- generics::tidy(res)
  expect_equal(td$conf_level[td$term %in% c("ci_lower", "grand_mean",
                                            "ci_upper")],
               rep(0.95, 3))
  expect_true(all(is.na(td$conf_level[td$term == "se"])))
})

test_that("the components plan the next study through ss_aipe_crossed()", {
  res <- ci_crossed(shrout_fleiss)
  plan <- ss_aipe_crossed(true_var_unit = val(res, "var_unit"),
                          true_var_rater = val(res, "var_rater"),
                          true_var_residual = val(res, "var_residual"), width = 2)
  # The planning values come back as the plan's echoes.
  expect_equal(val(plan, "true_var_rater"), val(res, "var_rater"))
  expect_lte(val(plan, "expected_width"), 2)
})

test_that("arrays the interval cannot use are refused", {
  x <- shrout_fleiss
  x[2, 3] <- NA
  expect_error(ci_crossed(x), "missing values")
  expect_error(ci_crossed(shrout_fleiss[, 1, drop = FALSE]), "at least 2")
  expect_error(ci_crossed(shrout_fleiss[1, , drop = FALSE]), "at least 2")
  expect_error(ci_crossed(matrix(3, 4, 3)), "do not vary")
  y <- shrout_fleiss
  y[1, 1] <- Inf
  expect_error(ci_crossed(y), "finite")
  d <- as.data.frame(shrout_fleiss)
  d$V1 <- as.character(d$V1)
  expect_error(ci_crossed(d), "must be numeric")
  expect_error(ci_crossed(letters[1:4]), "numeric matrix")
  expect_error(ci_crossed(shrout_fleiss, conf_level = 95), "'conf_level'")
})
