make_mlmr_missing_data <- function(seed = 113, n_miss_y = 5, n_miss_x = 8) {
  set.seed(seed)
  d <- mtcars
  d$hp[sample.int(nrow(d), n_miss_x)] <- NA
  d$mpg[sample.int(nrow(d), n_miss_y)] <- NA
  d
}


test_that("mlmr() returns an mlmr object with the expected slots", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a lavaan refit; the lm() agreement anchors run on CRAN
  fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
  expect_s3_class(fit, "mlmr")
  expect_named(fit$coefficients, c("(Intercept)", "wt", "hp"))
  expect_s3_class(fit$coef_table, "data.frame")
  expect_named(fit$coef_table, c("term", "estimate", "std_estimate",
                                  "se", "z_value", "p_value",
                                  "ci_lower", "ci_upper"))
  expect_true(is.matrix(fit$vcov))
  expect_equal(dim(fit$vcov), c(3L, 3L))
})


test_that("mlmr() coefficients match lm() to working precision on complete data", {
  skip_if_not_installed("lavaan")
  fit_m <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
  fit_l <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(unname(coef(fit_m)), unname(coef(fit_l)), tolerance = 1e-8)
})


test_that("mlmr() fits an intercept-only (null) model and composes with anova()", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # four lavaan refits; the lm() anchors run on CRAN
  null_fit <- mlmr(mpg ~ 1, data = mtcars, ci_method = "wald",
                   effect_sizes = FALSE)
  expect_s3_class(null_fit, "mlmr")
  expect_named(null_fit$coefficients, "(Intercept)")
  # the lone coefficient is the ML mean; an intercept-only model explains
  # none of the variance, so R^2 is 0
  expect_equal(unname(coef(null_fit)), mean(mtcars$mpg), tolerance = 1e-8)
  expect_equal(null_fit$R2, 0, tolerance = 1e-8)
  expect_equal(dim(null_fit$vcov), c(1L, 1L))
  expect_null(null_fit$effect_sizes)
  expect_equal(unname(predict(null_fit)[1]), mean(mtcars$mpg), tolerance = 1e-8)

  # the null model is the natural restricted model in a likelihood ratio test
  full_fit <- mlmr(mpg ~ wt, data = mtcars, ci_method = "wald",
                   effect_sizes = FALSE)
  lrt <- anova(null_fit, full_fit)
  expect_s3_class(lrt, "anova")
  expect_equal(nrow(lrt), 2L)

  # effect_sizes = TRUE on a null model is a no-op, not an error
  expect_error(mlmr(mpg ~ 1, data = mtcars), NA)
  # a formula with neither predictors nor an intercept still errors
  expect_error(mlmr(mpg ~ 0, data = mtcars), "nothing for mlmr")
})


test_that("mlmr() interactions and factors match lm()", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan refits; the lm() anchors run on CRAN
  fit_int_m <- mlmr(mpg ~ wt * hp, data = mtcars, ci_method = "wald")
  fit_int_l <- lm(mpg ~ wt * hp, data = mtcars)
  expect_equal(unname(coef(fit_int_m)), unname(coef(fit_int_l)),
               tolerance = 1e-6)

  fit_fac_m <- mlmr(mpg ~ wt + factor(cyl), data = mtcars,
                    ci_method = "wald")
  fit_fac_l <- lm(mpg ~ wt + factor(cyl), data = mtcars)
  expect_equal(unname(coef(fit_fac_m)), unname(coef(fit_fac_l)),
               tolerance = 1e-6)
})


test_that("mlmr() Wald and listwise reproduce lm() coefficients with missing X", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a lavaan refit; the lm() agreement anchors run on CRAN
  d <- make_mlmr_missing_data()
  fit_m <- mlmr(mpg ~ wt + hp, data = d, missing = "listwise",
                ci_method = "wald")
  fit_l <- lm(mpg ~ wt + hp, data = d)
  expect_equal(unname(coef(fit_m)), unname(coef(fit_l)), tolerance = 1e-8)
  expect_equal(nobs(fit_m), nobs(fit_l))
})


test_that("mlmr() FIML uses more data than listwise when X has missing values", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan refits; the lm() anchors run on CRAN
  d <- make_mlmr_missing_data()
  fit_fiml <- mlmr(mpg ~ wt + hp, data = d, ci_method = "wald")
  fit_lwd  <- mlmr(mpg ~ wt + hp, data = d, missing = "listwise",
                   ci_method = "wald")
  expect_gt(nobs(fit_fiml), nobs(fit_lwd))
  # The two estimators should produce different point estimates
  # (listwise discards rows under MAR; FIML uses them).
  expect_false(isTRUE(all.equal(unname(coef(fit_fiml)),
                                unname(coef(fit_lwd)),
                                tolerance = 1e-4)))
})


test_that("mlmr() profile likelihood CIs are produced and bracket the estimate", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- mlmr(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$ci_method, "profile")
  ci <- confint(fit)
  expect_true(all(is.finite(ci)))
  expect_true(all(ci[, 1L] < coef(fit)))
  expect_true(all(ci[, 2L] > coef(fit)))
})


test_that("mlmr() bootstrap CIs are produced and bracket the estimate", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(113)
  fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "boot",
              B = 100)
  expect_equal(fit$ci_method, "boot")
  ci <- confint(fit)
  expect_true(all(is.finite(ci)))
  expect_true(all(ci[, 1L] < coef(fit)))
  expect_true(all(ci[, 2L] > coef(fit)))
})


test_that("mlmr() predict() matches lm() on new data", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a lavaan refit; the lm() agreement anchors run on CRAN
  fit_m <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
  fit_l <- lm(mpg ~ wt + hp, data = mtcars)
  new_d <- data.frame(wt = c(2, 3, 4), hp = c(80, 150, 220))
  expect_equal(predict(fit_m, newdata = new_d),
               unname(predict(fit_l, newdata = new_d)),
               tolerance = 1e-6, ignore_attr = TRUE)
})


test_that("mlmr() residuals length matches data nrow and equals lm residuals on complete data", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a lavaan refit; the lm() agreement anchors run on CRAN
  fit_m <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
  fit_l <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(length(residuals(fit_m)), nrow(mtcars))
  expect_equal(unname(residuals(fit_m)), unname(residuals(fit_l)),
               tolerance = 1e-6)
})


test_that("mlmr() summary() prints without error and returns expected fields", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a lavaan refit; the lm() agreement anchors run on CRAN
  fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
  s <- summary(fit)
  expect_s3_class(s, "summary.mlmr")
  expect_true(all(c("coef_table", "R2", "adj_R2", "logLik", "N") %in%
                  names(s)))
  expect_output(print(s), "Coefficients")
  expect_output(print(s), "Confidence intervals")
})


test_that("mlmr() no-intercept formula fixes intercept to zero", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a lavaan refit; the lm() agreement anchors run on CRAN
  fit_m <- mlmr(mpg ~ 0 + wt + hp, data = mtcars, ci_method = "wald")
  fit_l <- lm(mpg ~ 0 + wt + hp, data = mtcars)
  # Iterative maximum likelihood against closed-form least squares:
  # convergence differs at the seventh digit across platforms (the
  # r-devel Windows and Debian pretests sit 2.5e-6 apart relatively).
  expect_equal(unname(coef(fit_m)), unname(coef(fit_l)), tolerance = 1e-5)
  expect_false("(Intercept)" %in% names(coef(fit_m)))
})


test_that("mlmr() errors on bad input", {
  skip_if_not_installed("lavaan")
  expect_error(mlmr(mpg ~ wt, data = matrix(0, 4, 2)),
               "must be a data.frame")
  expect_error(mlmr("mpg ~ wt", data = mtcars),
               "must be a formula")
  expect_error(mlmr(~ wt + hp, data = mtcars),
               "two-sided")
  expect_error(mlmr(mpg ~ wt, data = mtcars, conf_level = 1.1),
               "conf_level")
})


test_that("mlmr() model implied R-squared is between 0 and 1 and matches lm to working precision", {
  skip_if_not_installed("lavaan")
  fit_m <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
  fit_l <- lm(mpg ~ wt + hp, data = mtcars)
  expect_gt(fit_m$R2, 0)
  expect_lt(fit_m$R2, 1)
  expect_equal(fit_m$R2, summary(fit_l)$r.squared, tolerance = 1e-6)
})


test_that("mlmr() effect sizes (sr2, f2) are produced and finite", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # the effect sizes cost a refit per predictor
  fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
  es <- fit$effect_sizes
  expect_s3_class(es, "data.frame")
  expect_named(es, c("term", "sr2", "f2"))
  expect_equal(nrow(es), 2L)
  expect_true(all(is.finite(es$sr2)))
  expect_true(all(es$sr2 >= 0))
  expect_true(all(is.finite(es$f2)))
  expect_true(all(es$f2 >= 0))
  # The two slopes' sr^2 should not exceed the full model R^2.
  expect_lte(max(es$sr2), fit$R2 + 1e-8)
})


test_that("mlmr() omnibus LR test recovers the lm() F-test conclusion", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a lavaan refit; the lm() agreement anchors run on CRAN
  fit_m <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
  fit_l <- lm(mpg ~ wt + hp, data = mtcars)
  expect_lt(fit_m$omnibus_test$p_value, 0.01)
  # F vs chi^2: both should be highly significant on mtcars.
  expect_lt(summary(fit_l)$fstatistic[1L], Inf)
})


test_that("tidy.mlmr returns broom-style columns", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a lavaan refit; the lm() agreement anchors run on CRAN
  fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
  tdy <- tidy.mlmr(fit)
  expect_named(tdy, c("term", "estimate", "se", "statistic",
                      "p_value"))
  tdy_ci <- tidy.mlmr(fit, conf.int = TRUE)
  expect_true(all(c("ci_lower", "ci_upper") %in% colnames(tdy_ci)))
  tdy_std <- tidy.mlmr(fit, standardized = TRUE)
  expect_true("std_estimate" %in% colnames(tdy_std))
})


test_that("glance.mlmr returns a one-row data.frame with broom-style columns", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a lavaan refit; the lm() agreement anchors run on CRAN
  fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
  gl <- glance.mlmr(fit)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("R2", "adj_R2", "sigma", "logLik",
                    "AIC", "BIC", "nobs", "df_residual", "statistic",
                    "p_value") %in% colnames(gl)))
})


test_that("anova.mlmr compares nested fits with different formulas", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan refits; the lm() anchors run on CRAN
  fit1 <- mlmr(mpg ~ wt, data = mtcars, ci_method = "wald",
               effect_sizes = FALSE)
  fit2 <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
               effect_sizes = FALSE)
  res <- anova(fit1, fit2)
  expect_s3_class(res, "anova")
  expect_equal(nrow(res), 2L)
  # The LR p-value for adding hp should be highly significant on mtcars.
  expect_lt(res[["Pr(>Chisq)"]][2], 0.01)
})


test_that("anova.mlmr errors on non-nested formulas", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan refits; the lm() anchors run on CRAN
  fit_a <- mlmr(mpg ~ wt, data = mtcars, ci_method = "wald",
                effect_sizes = FALSE)
  fit_b <- mlmr(mpg ~ hp, data = mtcars, ci_method = "wald",
                effect_sizes = FALSE)
  expect_error(anova(fit_a, fit_b), "nested")
})


test_that("update.mlmr works via the default update() generic", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan refits; the lm() anchors run on CRAN
  fit_full <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
                   effect_sizes = FALSE)
  fit_red <- update(fit_full, . ~ . - hp)
  expect_s3_class(fit_red, "mlmr")
  expect_equal(names(coef(fit_red)), c("(Intercept)", "wt"))
})


test_that("MLR estimator produces robust SEs that differ from ML", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan refits; the lm() anchors run on CRAN
  fit_ml  <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
                  effect_sizes = FALSE)
  fit_mlr <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
                  estimator = "MLR", effect_sizes = FALSE)
  # Point estimates should match; SEs should differ for robust ones.
  expect_equal(unname(coef(fit_ml)), unname(coef(fit_mlr)),
               tolerance = 1e-6)
  ml_se  <- sqrt(diag(vcov(fit_ml)))
  mlr_se <- sqrt(diag(vcov(fit_mlr)))
  expect_false(isTRUE(all.equal(unname(ml_se), unname(mlr_se),
                                tolerance = 1e-4)))
})


test_that("confint(fit) does not warn when level argument is omitted", {
  skip_if_not_installed("lavaan")
  fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
              conf_level = 0.90, effect_sizes = FALSE)
  expect_silent(confint(fit))
})


test_that("confint(fit, level = X) warns only when level differs from fit-time", {
  skip_if_not_installed("lavaan")
  fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
              conf_level = 0.90, effect_sizes = FALSE)
  expect_warning(confint(fit, level = 0.95), "ignores 'level'")
  expect_silent(confint(fit, level = 0.90))
})


test_that("broom::tidy dispatches via generics", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("generics")
  fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
              effect_sizes = FALSE)
  out <- generics::tidy(fit)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("term", "estimate", "se", "statistic",
                      "p_value"))
})


test_that("broom::glance dispatches via generics", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("generics")
  fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
              effect_sizes = FALSE)
  out <- generics::glance(fit)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
})


test_that("anova.mlmr errors when fits use different missing-data settings", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan refits; the lm() anchors run on CRAN
  fit_fiml <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
                   effect_sizes = FALSE)
  fit_lwd  <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
                   missing = "listwise", effect_sizes = FALSE)
  expect_error(anova(fit_fiml, fit_lwd), "missing")
})


test_that("anova.mlmr errors when fits have different numbers of observations", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan refits; the lm() anchors run on CRAN
  fit1 <- mlmr(mpg ~ wt + hp, data = mtcars,             ci_method = "wald",
               effect_sizes = FALSE)
  d2 <- mtcars[1:25, ]
  fit2 <- mlmr(mpg ~ wt + hp, data = d2, ci_method = "wald",
               effect_sizes = FALSE)
  expect_error(anova(fit1, fit2), "same data")
  expect_error(anova(fit1, fit2), "32 observations")
})


test_that("anova.mlmr errors on two fits of equal N run on different data", {
  # HIGH-06 regression. Equal complete-case counts do not make two fits
  # comparable: if the response or a shared predictor holds different
  # values, the likelihoods are on different data and the LR test is
  # meaningless. Oracle: the two datasets are independently generated, so
  # by construction they are not the same data and anova() must refuse.
  skip_if_not_installed("lavaan")
  skip_on_cran()  # four lavaan refits on simulated data
  set.seed(113)
  n <- 120L
  d1 <- data.frame(x1 = rnorm(n))
  d1$y <- 10 + 3 * d1$x1 + rnorm(n)
  d2 <- data.frame(x1 = rnorm(n))
  d2$y <- -5 - 0.6 * d2$x1 + rnorm(n)
  f1 <- mlmr(y ~ x1, data = d1, ci_method = "wald", effect_sizes = FALSE)
  f2 <- mlmr(y ~ 1,  data = d2, ci_method = "wald", effect_sizes = FALSE)
  expect_identical(f1$N_complete, f2$N_complete)   # equal N; the trap
  expect_error(anova(f1, f2), "same data")
  # A subset of d1 with the same N as a differently generated set must
  # also error even when both fits carry a slope on x1.
  d3 <- data.frame(x1 = rnorm(n))
  d3$y <- rnorm(n)
  g1 <- mlmr(y ~ x1, data = d1, ci_method = "wald", effect_sizes = FALSE)
  g2 <- mlmr(y ~ 1,  data = d3, ci_method = "wald", effect_sizes = FALSE)
  expect_error(anova(g1, g2), "same data")
})


test_that("anova.mlmr same-data LRT matches the ML nested-model oracle", {
  # HIGH-06 companion: the same-data fix must not change the correct LR
  # statistic. Independent oracle: for a joint-normal FIML regression with
  # fixed.x = FALSE, the saturated predictor distribution is common to the
  # two nested fits, so the LR statistic reduces to the conditional
  # Y | X part, G^2 = N * log(RSS_reduced / RSS_full), with RSS taken from
  # ordinary least squares on the complete data.
  skip_if_not_installed("lavaan")
  skip_on_cran()  # three refits plus the constrained refits anova() makes
  fit_small <- mlmr(mpg ~ wt,      data = mtcars, ci_method = "wald",
                    effect_sizes = FALSE)
  fit_large <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
                    effect_sizes = FALSE)
  res <- anova(fit_small, fit_large)
  chisq_diff <- res[["Chisq diff"]][2L]

  lm_small <- lm(mpg ~ wt,      data = mtcars)
  lm_large <- lm(mpg ~ wt + hp, data = mtcars)
  N <- nrow(mtcars)
  rss_small <- sum(residuals(lm_small)^2)
  rss_large <- sum(residuals(lm_large)^2)
  oracle <- N * log(rss_small / rss_large)

  expect_equal(chisq_diff, oracle, tolerance = 1e-6)

  # Nesting with the shared predictor in a different position must still
  # be accepted (the same-data check compares by display name, not by the
  # internal positional column name).
  fit_hp <- mlmr(mpg ~ hp, data = mtcars, ci_method = "wald",
                 effect_sizes = FALSE)
  expect_s3_class(anova(fit_hp, fit_large), "anova")
})


test_that("anova.mlmr accepts nested FIML fits whose complete-case counts differ", {
  # The same-data guard must compare the analysis data, not a
  # complete-case count. Adding a predictor that is itself incompletely
  # observed lowers the number of rows complete on all modeled variables
  # without changing the observations, and comparing y ~ x1 with
  # y ~ x1 + x2 in that situation is exactly what FIML is for. A guard
  # keyed on N_complete refused it.
  #
  # Oracle: the observed data normal log-likelihood written out directly
  # (each row contributes the density of the subvector it observes) and
  # maximized outside DMAR and outside lavaan.
  #
  #   Larger model. With free means, free variances, and a free residual
  #   variance, y ~ x1 + x2 with x1 and x2 modeled as random is the
  #   saturated trivariate normal. Under this monotone pattern its MLE is
  #   closed form: the mean vector and divisor-N covariance of the
  #   complete (y, x1) block from all 200 rows, and the regression of x2
  #   on (y, x1) from the 160 rows where x2 is observed. That gives
  #   logLik = -795.3425413646, confirmed to ten decimal places by an
  #   unconstrained BFGS maximization over a Cholesky parameterization.
  #
  #   Smaller model. The same likelihood maximized subject to the partial
  #   slope of x2 being zero, parameterized by the joint normal of
  #   (x1, x2) plus the regression of y on x1 alone (Nelder-Mead followed
  #   by BFGS restarts): logLik = -803.0930756027.
  #
  # LR = 2 * (-795.3425413646 - (-803.0930756027)) = 15.5010684763 on
  # 1 degree of freedom, p = 8.24586e-05.
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two FIML refits on simulated data plus a constrained refit
  set.seed(113)
  n  <- 200L
  x1 <- rnorm(n)
  x2 <- 0.4 * x1 + rnorm(n)
  y  <- 0.5 * x1 + 0.3 * x2 + rnorm(n)
  d  <- data.frame(y = y, x1 = x1, x2 = x2)
  d$x2[sample(n, 40L)] <- NA

  m1 <- mlmr(y ~ x1,      data = d, missing = "fiml", ci_method = "wald",
             effect_sizes = FALSE)
  m2 <- mlmr(y ~ x1 + x2, data = d, missing = "fiml", ci_method = "wald",
             effect_sizes = FALSE)

  # The condition the old guard tripped on: same data, different
  # complete-case counts.
  expect_identical(m1$N_complete, 200L)
  expect_identical(m2$N_complete, 160L)

  res <- anova(m1, m2)
  expect_s3_class(res, "anova")
  expect_equal(res[["Chisq diff"]][2L], 15.5010684763, tolerance = 1e-7)
  expect_identical(as.integer(res[["Df diff"]][2L]), 1L)
  expect_equal(res[["Pr(>Chisq)"]][2L], 8.24586e-05, tolerance = 1e-5)
})


test_that("anova.mlmr still refuses nested FIML fits on genuinely different data", {
  # Companion to the test above: relaxing the complete-case requirement
  # must not let two different data frames through. Both data sets here
  # have 200 rows and 40 missing values on x2, so every count matches;
  # only the values differ. Oracle: the two data frames are independently
  # generated, so by construction they are not the same data.
  skip_if_not_installed("lavaan")
  skip_on_cran()  # three FIML refits on simulated data
  make_d <- function(seed) {
    set.seed(seed)
    n  <- 200L
    x1 <- rnorm(n)
    x2 <- 0.4 * x1 + rnorm(n)
    d  <- data.frame(y = 0.5 * x1 + 0.3 * x2 + rnorm(n), x1 = x1, x2 = x2)
    d$x2[sample(n, 40L)] <- NA
    d
  }
  d1 <- make_d(113)
  d2 <- make_d(114)
  m1 <- mlmr(y ~ x1,      data = d1, missing = "fiml", ci_method = "wald",
             effect_sizes = FALSE)
  m2 <- mlmr(y ~ x1 + x2, data = d2, missing = "fiml", ci_method = "wald",
             effect_sizes = FALSE)
  expect_identical(m1$N_complete, 200L)
  expect_identical(m2$N_complete, 160L)
  expect_error(anova(m1, m2), "same data")

  # Equal predictor sets on different data, where no constrained refit
  # intervenes, must also still error (the commit that introduced the
  # analysis-data comparison).
  m3 <- mlmr(y ~ x1 + x2, data = d1, missing = "fiml", ci_method = "wald",
             effect_sizes = FALSE)
  expect_error(anova(m3, m2), "same data")
})


test_that("with only the outcome missing, mlmr slopes and SEs match lm (MEDIUM-10)", {
  # Supports the corrected documentation. When missingness is confined to
  # the outcome, the rows with an observed predictor but a missing outcome
  # inform the marginal distribution of X, not the conditional
  # distribution of Y given X that identifies the slopes, so FIML gives no
  # slope-SE advantage over listwise deletion. Oracle: lm() on the
  # complete cases, which is independent of the FIML machinery.
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a FIML fit on 300 simulated rows
  set.seed(113)
  N <- 300L
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  y  <- 1 + 0.5 * x1 - 0.3 * x2 + rnorm(N)
  d  <- data.frame(y = y, x1 = x1, x2 = x2)
  d$y[sample(N, 120L)] <- NA          # ~120 outcomes missing, X fully observed

  fit  <- mlmr(y ~ x1 + x2, data = d, ci_method = "wald")
  l_cc <- lm(y ~ x1 + x2, data = d)   # listwise deletion == complete cases

  # Slopes coincide with listwise to working precision: the extra X-only
  # rows do not move the conditional-model estimates. (Measured max
  # absolute coefficient difference here is on the order of 1e-6.)
  expect_equal(unname(coef(fit)), unname(coef(l_cc)), tolerance = 1e-4)

  # The raw standard errors differ, but the difference is entirely the
  # maximum likelihood N vs (N - K - 1) variance divisor, not information
  # recovered from the X-only rows. Rescaling the lm standard errors onto
  # the ML N divisor removes the whole discrepancy (measured max absolute
  # difference on the order of 1e-7), which is the point of the corrected
  # documentation: there is no slope-SE advantage when only the outcome is
  # missing.
  n_obs    <- sum(!is.na(d$y))
  p        <- length(coef(l_cc))
  se_fiml  <- sqrt(diag(vcov(fit)))
  se_lm_ml <- sqrt(diag(vcov(l_cc))) * sqrt((n_obs - p) / n_obs)
  expect_equal(unname(se_fiml), unname(se_lm_ml), tolerance = 1e-4)
})


test_that("bootstrap CIs are reproducible when boot_seed is supplied", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit_a <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "boot",
                B = 50, boot_seed = 113, effect_sizes = FALSE)
  fit_b <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "boot",
                B = 50, boot_seed = 113, effect_sizes = FALSE)
  expect_equal(fit_a$ci, fit_b$ci, tolerance = 1e-10)
})


test_that("boot_seed restores the user RNG state", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(42)
  before <- .Random.seed
  fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "boot",
              B = 50, boot_seed = 113, effect_sizes = FALSE)
  after <- .Random.seed
  expect_identical(before, after)
})


test_that("enforce_es_bounds clamps negative sr2/f2 to zero", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two fits, each with a refit per predictor
  # Use a small/noisy design where one slope is essentially zero so
  # the constrained and full models have nearly identical sigma2_e
  # and the raw sr2 may dip below zero in some samples.
  set.seed(113)
  d <- data.frame(y = rnorm(40), x1 = rnorm(40), x2 = rnorm(40))
  fit_raw <- mlmr(y ~ x1 + x2, data = d, ci_method = "wald")
  fit_clamped <- mlmr(y ~ x1 + x2, data = d, ci_method = "wald",
                      enforce_es_bounds = TRUE)
  expect_true(all(fit_clamped$effect_sizes$sr2 >= 0))
  expect_true(all(fit_clamped$effect_sizes$f2 >= 0))
})


test_that("auxiliary variables leave complete-data coefficients unchanged", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan refits; auxiliary validation runs on CRAN
  set.seed(113)
  d <- mtcars
  d$aux <- d$qsec + rnorm(nrow(d))
  f0 <- mlmr(mpg ~ wt + hp, data = d, ci_method = "wald",
             effect_sizes = FALSE)
  fa <- mlmr(mpg ~ wt + hp, data = d, ci_method = "wald",
             effect_sizes = FALSE, auxiliary = "aux")
  # Saturated correlates do not change the focal regression on
  # complete data; the coefficients must match to working precision.
  expect_equal(unname(coef(f0)), unname(coef(fa)), tolerance = 1e-6)
  expect_identical(fa$auxiliary, "aux")
  expect_null(f0$auxiliary)
})


test_that("auxiliary variables validate their inputs", {
  skip_if_not_installed("lavaan")
  expect_error(mlmr(mpg ~ wt, data = mtcars, auxiliary = "nope"),
               "not found in 'data'")
  expect_error(mlmr(mpg ~ wt, data = mtcars, auxiliary = "wt"),
               "also appear in 'formula'")
  expect_error(mlmr(mpg ~ wt, data = mtcars, auxiliary = "hp",
                    fixed_x = TRUE),
               "fixed_x = FALSE")
  d <- mtcars
  d$grp <- factor(rep(c("a", "b"), length.out = nrow(d)))
  expect_error(mlmr(mpg ~ wt, data = d, auxiliary = "grp"),
               "must be numeric")
})


test_that("an auxiliary recovers information lost to MAR dropout", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a FIML fit on 4000 simulated rows with MAR dropout
  set.seed(113)
  n <- 4000
  X <- rnorm(n)
  Z <- rnorm(n)
  Tt <- rbinom(n, 1, 0.5)
  Y <- 0.3 * X + 0.6 * Z + 0.5 * Tt + rnorm(n)
  d <- data.frame(T = Tt, X = X, Z = Z, Y = Y)
  # Outcome drops out as a function of the auxiliary Z (MAR given Z).
  d$Y[runif(n) < plogis(1.5 * (d$Z - 0.3))] <- NA
  truth_mean <- 0.5 * mean(Tt) + 0.3 * mean(X)
  fa <- mlmr(Y ~ T + X, data = d, ci_method = "wald",
             effect_sizes = FALSE, auxiliary = "Z")
  # The complete-case mean is badly biased; conditioning on the
  # dropout-related auxiliary recovers the population mean.
  expect_gt(abs(mean(d$Y, na.rm = TRUE) - truth_mean), 0.15)
  expect_lt(abs(mean(predict(fa)) - truth_mean), 0.1)
  # The treatment contrast is recovered near its population value.
  expect_equal(unname(coef(fa)["T"]), 0.5, tolerance = 0.08)
})


test_that(".mlmr_require_lavaan guards the detectCores() == NA case (MEDIUM-09)", {
  # lavaan 0.6-19 builds its default 'ncpus' from
  # max(1L, parallel::detectCores() - 1L); an undetectable core count
  # makes that NA and lavaan then dies in option validation with the
  # opaque message "missing value where TRUE/FALSE needed". The guard
  # converts that into an actionable upgrade message and lets a fixed
  # lavaan through. Oracle: hand-supplied core counts and a version
  # comparison against the 0.7-2 floor, independent of the installed
  # lavaan and of the fit machinery.
  skip_if_not_installed("lavaan")

  # Undetectable cores with a pre-0.7-2 lavaan is the known-bad case.
  expect_error(
    .mlmr_require_lavaan("mlmr", n_cores = NA_integer_,
                         lavaan_version = "0.6-19"),
    "Upgrade lavaan"
  )
  # A detectable core count is fine on any version.
  expect_true(.mlmr_require_lavaan("mlmr", n_cores = 8L,
                                   lavaan_version = "0.6-19"))
  # The 0.7-2 floor, and anything newer, is trusted even with NA cores.
  expect_true(.mlmr_require_lavaan("mlmr", n_cores = NA_integer_,
                                   lavaan_version = "0.7-2"))
  # Version ordering must be numeric, not lexicographic: 0.10-1 > 0.7-2,
  # so an undetectable core count on 0.10-1 must not error.
  expect_true(.mlmr_require_lavaan("mlmr", n_cores = NA_integer_,
                                   lavaan_version = "0.10-1"))
  # mlmr_mv() shares the guard and names itself in the message.
  expect_error(
    .mlmr_require_lavaan("mlmr_mv", n_cores = NA_integer_,
                         lavaan_version = "0.6-19"),
    "mlmr_mv\\(\\) fails inside lavaan"
  )
})
