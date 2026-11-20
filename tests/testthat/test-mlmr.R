# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

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
                                  "se", "test_statistic", "p_value",
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
                      "p_value", "wald_df"))
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
                      "p_value", "wald_df"))
  expect_identical(out$wald_df, rep(Inf, 3))
  # The t reference travels in the rows, so that a stack of tidy tables
  # still says which rows are on it (it was only an attribute, which
  # rbind() keeps from the first table).
  fit_t <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
                effect_sizes = FALSE, wald_df = 29)
  both <- rbind(out, generics::tidy(fit_t))
  expect_identical(both$wald_df, rep(c(Inf, 29), each = 3))
  expect_equal(both$p_value[4:6],
               2 * stats::pt(-abs(both$statistic[4:6]), 29))
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


test_that("an auxiliary is refused under listwise deletion", {
  skip_if_not_installed("lavaan")
  # lavaan's listwise deletion drops every row missing any variable of the
  # model, the auxiliary included. With t26_flags, observed for the 145
  # Grant-White pupils, as the auxiliary, mlmr() fitted 145 of the 301
  # rows and returned coefficients of 0.971, 0.313, and 0.178 where lm()
  # gives 1.117, 0.317, and 0.167; an auxiliary carries information only
  # through the full information likelihood.
  f <- t6_paragraph_comprehension ~ t7_sentence + t9_word_meaning
  expect_error(mlmr(f, data = holzinger_swineford, missing = "listwise",
                    auxiliary = "t26_flags"),
               "cannot be combined with missing = \"listwise\"")
  expect_error(mlmr(f, data = holzinger_swineford, missing = "listwise",
                    auxiliary = "t13_straight_and_curved_capitals"),
               "set missing = \"fiml\"")
})

test_that("an auxiliary is refused under pairwise deletion, where it changes nothing", {
  skip_if_not_installed("lavaan")
  # Under missing = "pairwise" (and "available.cases") the moments are
  # computed from each pair of variables on its own, and with the model
  # saturated the coefficients depend on the moments of the outcome and
  # the predictors alone. Before 1.1.0 the call was accepted: with t13 as
  # the auxiliary, every coefficient and standard error of t6 on t7, t9,
  # and t26 agreed with the fit without it to 1e-7, while the result
  # recorded the auxiliary as used.
  f <- t6_paragraph_comprehension ~ t7_sentence + t9_word_meaning +
    t26_flags
  for (m in c("pairwise", "available.cases")) {
    expect_error(mlmr(f, data = holzinger_swineford, missing = m,
                      auxiliary = "t13_straight_and_curved_capitals"),
                 paste0("'auxiliary' cannot be combined with missing = \"",
                        m, "\""), fixed = TRUE)
    expect_error(mlmr(f, data = holzinger_swineford, missing = m,
                      auxiliary = "t13_straight_and_curved_capitals"),
                 "leave every estimate unchanged. Drop 'auxiliary' or set missing = \"fiml\".",
                 fixed = TRUE)
  }
  # The refusal comes before anything is fit, and the same call without the
  # auxiliary, or with it under FIML, still runs.
  expect_s3_class(mlmr(f, data = holzinger_swineford, missing = "pairwise",
                       ci_method = "wald", effect_sizes = FALSE), "mlmr")
  fit <- mlmr(f, data = holzinger_swineford, ci_method = "wald",
              effect_sizes = FALSE,
              auxiliary = "t13_straight_and_curved_capitals")
  expect_s3_class(fit, "mlmr")
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


# Paragraph comprehension (t6) regressed on sentence completion (t7) and
# word meaning (t9), with t9 made missing at random given t7 (110 of the
# 301 scores), and straight and curved capitals (t13) as an auxiliary.
# This is the example on which the defect of 1.0.x was found.
.mlmr_aux_example <- function() {
  v <- c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning",
         "t13_straight_and_curved_capitals")
  d <- holzinger_swineford[v]
  set.seed(113)
  z <- as.numeric(scale(d$t7_sentence))
  d$t9_word_meaning[stats::runif(nrow(d)) < stats::plogis(-1 + 2 * z)] <- NA
  d
}

test_that("an auxiliary leaves the covariances among the predictors free", {
  # Regression test for a defect of 1.0.x. Writing "aux ~~ x_k" made the
  # predictors endogenous in lavaan's sense, lavaan stopped adding their
  # covariances, and every covariance among them was fixed at zero. On
  # this example the model that should be saturated had chi square 129.5
  # on 1 df, an implied covariance of 0 between t7 and t9 (observed
  # correlation .69), and a t7 slope of 0.449, against 0.362 without the
  # auxiliary and 0.317 on the complete data.
  skip_if_not_installed("lavaan")
  d <- .mlmr_aux_example()
  expect_identical(sum(is.na(d$t9_word_meaning)), 110L)
  f <- t6_paragraph_comprehension ~ t7_sentence + t9_word_meaning
  aux <- "t13_straight_and_curved_capitals"
  fa <- mlmr(f, data = d, ci_method = "wald", effect_sizes = FALSE,
             auxiliary = aux)
  f0 <- mlmr(f, data = d, ci_method = "wald", effect_sizes = FALSE)

  # Saturated: no degrees of freedom and no misfit.
  fm <- lavaan::fitMeasures(fa$lavaan_fit, c("chisq", "df", "npar"))
  expect_equal(unname(fm[["df"]]), 0)
  expect_lt(unname(fm[["chisq"]]), 1e-6)
  expect_equal(unname(fm[["npar"]]), 14)

  # The predictors' covariance, variances, and means are free parameters.
  pt <- lavaan::parTable(fa$lavaan_fit)
  free_row <- function(lhs, op, rhs) {
    pt$free[(pt$lhs == lhs & pt$op == op & pt$rhs == rhs) |
              (pt$lhs == rhs & pt$op == op & pt$rhs == lhs)]
  }
  expect_gt(free_row("x_1", "~~", "x_2"), 0)
  expect_gt(free_row("x_1", "~~", "x_1"), 0)
  expect_gt(free_row("x_2", "~~", "x_2"), 0)
  expect_gt(free_row("x_1", "~1", ""), 0)
  expect_gt(free_row("x_2", "~1", ""), 0)
  expect_gt(free_row("aux_1", "~1", ""), 0)

  # Oracle: in a saturated model the maximum likelihood estimates are
  # those of the unrestricted means and covariances, which cov_fiml()
  # finds by the EM algorithm, a different route from lavaan's
  # quasi-Newton fit of the regression. The regression they imply is the
  # mlmr() fit, and so is their covariance between the predictors.
  em <- cov_fiml(d, aux = aux)
  x <- c("t7_sentence", "t9_word_meaning")
  y <- "t6_paragraph_comprehension"
  b <- solve(em$S[x, x], em$S[x, y])
  b0 <- em$mu[[y]] - sum(b * em$mu[x])
  expect_equal(unname(coef(fa)), unname(c(b0, b)), tolerance = 1e-6)
  implied <- lavaan::lavInspect(fa$lavaan_fit, "cov.ov")
  expect_equal(implied["x_1", "x_2"], em$S[x[1], x[2]], tolerance = 1e-6)
  expect_equal(unname(lavaan::lavInspect(fa$lavaan_fit, "mean.ov")[
    c("x_1", "x_2")]), unname(em$mu[x]), tolerance = 1e-6)
  expect_equal(fa$R2,
               1 - (em$S[y, y] - sum(em$S[y, x] * b)) / em$S[y, y],
               tolerance = 1e-6)

  # The values themselves (semTools::sem.auxiliary() gives the same to
  # 6e-9), and the complete data and no-auxiliary fits beside them.
  expect_equal(unname(coef(fa)), c(0.7864647905, 0.3619103682, 0.1430313507),
               tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(fa)))),
               c(0.4741596925, 0.0435295850, 0.0367545038), tolerance = 1e-5)
  expect_equal(unname(coef(f0)), c(0.7841720544, 0.3621366904, 0.1429605205),
               tolerance = 1e-6)
  complete <- coef(stats::lm(f, data = holzinger_swineford))
  expect_equal(unname(complete), c(1.1169740321, 0.3174143443, 0.1669888063),
               tolerance = 1e-8)
  # t9 is missing at random given t7, which is in the model, so both fits
  # are consistent and t13, which correlates about .2 with the three
  # tests, moves them little; each coefficient lies nearer the complete
  # data value with the auxiliary than without it.
  expect_true(all(abs(coef(fa) - complete) < abs(coef(f0) - complete)))
})

test_that("an auxiliary that predicts the missingness moves the slopes toward the complete data", {
  # Word meaning (t9) missing at random given general information (t5),
  # which is not in the model: the fit without t5 is biased, and t5 as an
  # auxiliary recovers most of the bias. Under 1.0.x the auxiliary fit
  # was worse than either (t7 slope 0.438 and t9 slope 0.090, against
  # 0.317 and 0.167 on the complete data), since it treated t7 and t9 as
  # uncorrelated.
  skip_if_not_installed("lavaan")
  v <- c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning",
         "t5_general_information")
  d <- holzinger_swineford[v]
  set.seed(113)
  z <- as.numeric(scale(d$t5_general_information))
  d$t9_word_meaning[stats::runif(nrow(d)) < stats::plogis(-1 + 2 * z)] <- NA
  expect_identical(sum(is.na(d$t9_word_meaning)), 100L)
  f <- t6_paragraph_comprehension ~ t7_sentence + t9_word_meaning
  fa <- mlmr(f, data = d, ci_method = "wald", effect_sizes = FALSE,
             auxiliary = "t5_general_information")
  f0 <- mlmr(f, data = d, ci_method = "wald", effect_sizes = FALSE)
  complete <- coef(stats::lm(f, data = holzinger_swineford))
  expect_equal(unname(coef(fa)), c(0.8056606, 0.3376791, 0.1713900),
               tolerance = 1e-5)
  expect_true(all(abs(coef(fa) - complete) < abs(coef(f0) - complete)))
  expect_lt(unname(lavaan::fitMeasures(fa$lavaan_fit, "chisq")), 1e-6)
})

test_that("with complete data an auxiliary changes no regression result", {
  # Under 1.0.x the slopes agreed with lm() on complete data (their
  # maximum likelihood estimates do not depend on the predictors'
  # covariance there), but the implied variance of the outcome left the
  # predictors' covariance out: R^2 was .471 rather than lm()'s .602 on
  # this model, and the standardized coefficients and semi-partial R^2
  # were wrong with it.
  skip_if_not_installed("lavaan")
  aux <- "t13_straight_and_curved_capitals"
  f <- t6_paragraph_comprehension ~ t7_sentence + t9_word_meaning + sex
  fa <- mlmr(f, data = holzinger_swineford, ci_method = "wald",
             auxiliary = aux)
  f0 <- mlmr(f, data = holzinger_swineford, ci_method = "wald")
  ref <- stats::lm(f, data = holzinger_swineford)
  fm <- lavaan::fitMeasures(fa$lavaan_fit, c("chisq", "df"))
  expect_equal(unname(fm[["df"]]), 0)
  expect_lt(unname(fm[["chisq"]]), 1e-6)
  expect_equal(unname(coef(fa)), unname(coef(ref)), tolerance = 1e-6)
  expect_equal(unname(coef(f0)), unname(coef(ref)), tolerance = 1e-6)
  expect_equal(fa$R2, summary(ref)$r.squared, tolerance = 1e-6)
  expect_equal(fa$adj_R2, summary(ref)$adj.r.squared, tolerance = 1e-6)
  expect_equal(fa$coef_table$std_estimate, f0$coef_table$std_estimate,
               tolerance = 1e-6)
  expect_equal(fa$effect_sizes$sr2, f0$effect_sizes$sr2, tolerance = 1e-5)
  expect_equal(fa$effect_sizes$f2, f0$effect_sizes$f2, tolerance = 1e-5)
  expect_equal(fa$omnibus_test$statistic, f0$omnibus_test$statistic,
               tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(fa)))), unname(sqrt(diag(vcov(f0)))),
               tolerance = 1e-6)
})

test_that("the three FMI estimators agree on a fit with an auxiliary", {
  # The model is saturated with or without auxiliaries, so the Hessian,
  # structured, and unstructured estimators of Chen and Savalei (2021)
  # coincide. Under 1.0.x they did not on this example (t7 slope:
  # unstructured .016, Hessian -.559), and wald_df =
  # "missing_information" gave its slope 293.1 degrees of freedom.
  skip_if_not_installed("lavaan")
  d <- .mlmr_aux_example()
  f <- t6_paragraph_comprehension ~ t7_sentence + t9_word_meaning
  fa <- mlmr(f, data = d, ci_method = "wald", effect_sizes = FALSE,
             auxiliary = "t13_straight_and_curved_capitals")
  u <- fraction_missing_information(fa, information = "unstructured")
  s <- fraction_missing_information(fa, information = "structured")
  h <- fraction_missing_information(fa, information = "hessian")
  expect_identical(u$term, h$term)
  expect_true("t7_sentence ~~ t9_word_meaning" %in% u$term)
  expect_equal(u$fmi, s$fmi, tolerance = 1e-6)
  expect_equal(u$fmi, h$fmi, tolerance = 1e-6)
  coefs <- u$op %in% c("~1", "~") & grepl("^t6_", u$term)
  expect_equal(u$fmi[coefs], c(0.0399233, 0.2542425, 0.3593729),
               tolerance = 1e-5)
  fw <- mlmr(f, data = d, ci_method = "wald", effect_sizes = FALSE,
             auxiliary = "t13_straight_and_curved_capitals",
             wald_df = "missing_information")
  expect_equal(unname(attr(fw, "wald_df")),
               c(285.9830751, 221.4729933, 189.8287479), tolerance = 1e-6)
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


# ---- wald_df: a t reference distribution for the Wald interval -------------

test_that("mlmr() with wald_df = Inf reproduces the normal intervals and p-values exactly", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # four lavaan fits; the df_missing() anchors run on CRAN
  strip <- function(f) {
    f$call <- NULL
    f$lavaan_fit <- NULL
    f
  }
  for (ci in c("wald", "profile")) {
    a <- mlmr(t6_paragraph_comprehension ~ t7_sentence + t26_flags,
              data = holzinger_swineford, ci_method = ci,
              effect_sizes = FALSE)
    b <- mlmr(t6_paragraph_comprehension ~ t7_sentence + t26_flags,
              data = holzinger_swineford, ci_method = ci,
              effect_sizes = FALSE, wald_df = Inf)
    expect_identical(strip(a), strip(b))
    expect_identical(attr(a, "wald_df"), Inf)
    expect_null(attr(a, "wald_df_table"))
  }
  tab <- a$coef_table
  expect_identical(tab$p_value, 2 * stats::pnorm(-abs(tab$test_statistic)))
})

test_that("a finite wald_df refers mlmr()'s Wald interval and p-values to t", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan fits
  a <- mlmr(t6_paragraph_comprehension ~ t7_sentence + t26_flags,
            data = holzinger_swineford, ci_method = "wald",
            effect_sizes = FALSE)
  r <- mlmr(t6_paragraph_comprehension ~ t7_sentence + t26_flags,
            data = holzinger_swineford, ci_method = "wald",
            effect_sizes = FALSE, wald_df = 25)
  ta <- a$coef_table
  tr <- r$coef_table
  expect_identical(tr[c("term", "estimate", "se", "test_statistic")],
                   ta[c("term", "estimate", "se", "test_statistic")])
  crit <- stats::qt(0.975, 25)
  expect_equal(tr$ci_lower, ta$estimate - crit * ta$se, tolerance = 1e-14)
  expect_equal(tr$ci_upper, ta$estimate + crit * ta$se, tolerance = 1e-14)
  expect_identical(tr$p_value, 2 * stats::pt(-abs(ta$test_statistic), 25))
  expect_equal(unname(confint(r)), unname(as.matrix(tr[c("ci_lower",
                                                         "ci_upper")])))
  expect_identical(attr(r, "wald_df"), 25)
  expect_true(any(grepl("t reference on 25 degrees of freedom",
                        capture.output(print(summary(r))))))
  expect_false(any(grepl("t reference",
                         capture.output(print(summary(a))))))
})

test_that("summary() heads the p-values with their reference distribution", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # three lavaan fits
  # With a finite wald_df the p-values are 2 P(T > |z|) on wald_df degrees
  # of freedom (here .000571, .009346, .008414 on 5, against normal tail
  # probabilities of 9e-15, 4e-5, and 3e-5), and the printed table was
  # still headed Pr(>|z|). It is headed as summary.lm() heads it, and the
  # line naming the degrees of freedom sits directly beneath it.
  f <- t26_flags ~ t1_visual_perception + t4_lozenges
  t5 <- mlmr(f, data = holzinger_swineford, ci_method = "wald",
             effect_sizes = FALSE, wald_df = 5)
  expect_equal(t5$coef_table$p_value,
               2 * stats::pt(-abs(t5$coef_table$test_statistic), 5))
  printed <- capture.output(print(summary(t5)))
  head_row <- grep("Std. Error", printed, fixed = TRUE)
  expect_length(head_row, 1L)
  expect_match(printed[head_row], "t value +Pr\\(>\\|t\\|\\)")
  expect_false(any(grepl("Pr(>|z|)", printed, fixed = TRUE)))
  expect_lt(grep("t reference on 5 degrees of freedom", printed),
            grep("Confidence intervals", printed, fixed = TRUE))
  mi <- mlmr(f, data = holzinger_swineford, ci_method = "wald",
             effect_sizes = FALSE, wald_df = "missing_information")
  printed <- capture.output(print(summary(mi)))
  expect_match(printed[grep("Std. Error", printed, fixed = TRUE)],
               "t value +Pr\\(>\\|t\\|\\)")
  expect_lt(grep("von Hippel's (2016)", printed, fixed = TRUE),
            grep("Confidence intervals", printed, fixed = TRUE))
  normal <- mlmr(f, data = holzinger_swineford, ci_method = "wald",
                 effect_sizes = FALSE)
  printed <- capture.output(print(summary(normal)))
  expect_match(printed[grep("Std. Error", printed, fixed = TRUE)],
               "z value +Pr\\(>\\|z\\|\\)")
  # The stored column keeps its name.
  expect_true("test_statistic" %in% names(t5$coef_table))
})

test_that("wald_df = 'missing_information' applies von Hippel's rule per coefficient", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # FIML fits and the two refits behind the FMI
  # t26_flags is observed for the 145 Grant-White students of 301, and the
  # predictors are complete. Each coefficient's degrees of freedom are
  # max(3, N(1 - FMI) - k) with N = 301, k = 3 (intercept and two slopes),
  # and the FMI that fraction_missing_information() returns for the fit:
  # 301(1 - .5364747) - 3 = 136.5211 and 301(1 - .5979813) - 3 = 118.0076.
  fit <- mlmr(t26_flags ~ t1_visual_perception + t4_lozenges,
              data = holzinger_swineford, ci_method = "wald",
              effect_sizes = FALSE, wald_df = "missing_information")
  wd <- attr(fit, "wald_df")
  expect_named(wd, c("(Intercept)", "t1_visual_perception", "t4_lozenges"))
  expect_equal(unname(wd[2:3]), c(136.521123, 118.0076156),
               tolerance = 1e-6)
  fmi <- fraction_missing_information(fit)
  fmi_coef <- fmi$fmi[fmi$op %in% c("~1", "~") &
                        grepl("^t26_flags", fmi$term)]
  expect_equal(unname(wd), pmax(3, 301 * (1 - fmi_coef) - 3),
               tolerance = 1e-12)
  # The same numbers df_missing() gives from the FMI table.
  dm <- df_missing(fmi[fmi$op %in% c("~1", "~") &
                         grepl("^t26_flags", fmi$term), ], k = 3)
  expect_equal(unname(wd), dm$df, tolerance = 1e-12)
  tbl <- attr(fit, "wald_df_table")
  expect_s3_class(tbl, "dmar_df_missing")
  expect_identical(tbl$term, names(wd))
  expect_identical(tbl$N, c(301, 301, 301))
  expect_identical(tbl$k, c(3, 3, 3))
  # The intervals and p-values use them.
  tab <- fit$coef_table
  crit <- stats::qt(0.975, wd)
  expect_equal(tab$ci_upper, tab$estimate + unname(crit) * tab$se,
               tolerance = 1e-14)
  expect_equal(tab$p_value, 2 * stats::pt(-abs(tab$test_statistic), unname(wd)),
               tolerance = 1e-15)
  printed <- capture.output(print(summary(fit)))
  expect_true(any(grepl("von Hippel's \\(2016\\)", printed)))
  expect_true(any(grepl("t4_lozenges +0.598 +118.0", printed)))
})

test_that("wald_df = 'missing_information' stops when a coefficient's FMI cannot be computed", {
  skip_if_not_installed("lavaan")
  # The FMI behind each coefficient's degrees of freedom comes from two
  # refits of the fit, and .fmi_core() stops with the reason when one
  # cannot be computed. An FMI that came back NA all the same gave NA
  # degrees of freedom, p-values, and limits without a word (with an
  # outcome whose standard deviation was near 1e8, before the refits were
  # rescaled). The stand-in returns NA for every FMI.
  local_mocked_bindings(.fmi_core = function(fit, information) {
    pt <- lavaan::parTable(fit)
    pt$fmi <- NA_real_
    pt
  })
  expect_error(
    mlmr(t26_flags ~ t1_visual_perception + t4_lozenges,
         data = holzinger_swineford, ci_method = "wald",
         effect_sizes = FALSE, wald_df = "missing_information"),
    paste0("could not compute the fraction of missing information of ",
           "(Intercept), t1_visual_perception, t4_lozenges from the fit, ",
           "so mlmr() has no degrees of freedom for their intervals"),
    fixed = TRUE)
})

test_that("wald_df = 'missing_information' is refused where lavaan replaced the FIML fit with listwise deletion", {
  skip_if_not_installed("lavaan")
  # lavaan fits the complete cases in place of the FIML fit missing =
  # "fiml" asks for when the robust.sem sandwich or a test statistic of the
  # Satorra-Bentler kind reaches it, and the keyword then gave every
  # coefficient FMI = 0 and the N - k of the complete cases (142 of the
  # 145 Grant-White students here, where the FIML fit gives 136.5 and
  # 118.0 to the slopes).
  f <- t26_flags ~ t1_visual_perception + t4_lozenges
  msg <- paste0("and lavaan fit the complete cases instead (145 of the 301 ",
                "cases read), as it does with ")
  expect_error(suppressWarnings(
    mlmr(f, data = holzinger_swineford, ci_method = "wald",
         effect_sizes = FALSE, wald_df = "missing_information",
         test = "satorra.bentler")),
    paste0(msg, "test = \"satorra.bentler\""), fixed = TRUE)
  expect_error(suppressWarnings(
    mlmr(f, data = holzinger_swineford, ci_method = "wald",
         effect_sizes = FALSE, wald_df = "missing_information",
         se = "robust.sem")),
    paste0(msg, "se = \"robust.sem\""), fixed = TRUE)
  # Listwise deletion asked for is the fit the call asked for.
  lw <- suppressWarnings(
    mlmr(f, data = holzinger_swineford, ci_method = "wald",
         effect_sizes = FALSE, wald_df = "missing_information",
         missing = "listwise", test = "satorra.bentler"))
  expect_identical(unname(attr(lw, "wald_df")), rep(142, 3))
  # With no value missing the listwise fit is the FIML fit, and every
  # coefficient has the N - k of lm(), as the help says. The refusal above
  # had refused these calls as well, saying that lavaan fit the complete
  # cases instead, "145 of the 145 cases read".
  cc <- holzinger_swineford[stats::complete.cases(
    holzinger_swineford[all.vars(f)]), ]
  for (arg in list(list(se = "robust.sem"), list(test = "satorra.bentler"))) {
    fit <- suppressWarnings(do.call(mlmr, c(list(
      f, data = cc, ci_method = "wald", effect_sizes = FALSE,
      wald_df = "missing_information"), arg)))
    expect_identical(unname(attr(fit, "wald_df")), rep(142, 3))
  }
})

test_that("with complete data the missing information rule gives N - k exactly", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan fits
  # No value is missing, so every FMI is zero by definition and each
  # coefficient gets the complete data degrees of freedom of lm().
  d <- holzinger_swineford[, c("t6_paragraph_comprehension", "t7_sentence",
                               "t9_word_meaning")]
  fit <- mlmr(t6_paragraph_comprehension ~ t7_sentence + t9_word_meaning,
              data = d, ci_method = "wald", effect_sizes = FALSE,
              wald_df = "missing_information")
  expect_identical(unname(attr(fit, "wald_df")), rep(301 - 3, 3))
  expect_identical(attr(fit, "wald_df_table")$fmi, c(0, 0, 0))
  # The interval is lm()'s t interval on the same degrees of freedom with
  # the maximum likelihood standard error, whose variance divides by N
  # rather than N - k, so it is shorter by sqrt((N - k)/N); von Hippel's
  # construction, without a variance correction.
  ols <- lm(t6_paragraph_comprehension ~ t7_sentence + t9_word_meaning,
            data = d)
  width <- function(ci) ci[, 2] - ci[, 1]
  expect_equal(unname(width(confint(fit)) / width(confint(ols))),
               rep(sqrt(298 / 301), 3), tolerance = 1e-6)
  # A listwise fit analyzes complete cases only: N_complete - k.
  lw <- mlmr(t26_flags ~ t1_visual_perception + t4_lozenges,
             data = holzinger_swineford, ci_method = "wald",
             missing = "listwise", effect_sizes = FALSE,
             wald_df = "missing_information")
  expect_identical(unname(attr(lw, "wald_df")), rep(145 - 3, 3))
})

test_that("mlmr() refuses wald_df with intervals other than Wald", {
  skip_if_not_installed("lavaan")
  for (ci in c("profile", "boot")) {
    expect_error(mlmr(t6_paragraph_comprehension ~ t7_sentence,
                      data = holzinger_swineford, ci_method = ci,
                      wald_df = 25),
                 paste0("ci_method = \"", ci, "\" is not one"))
  }
  # The default interval is the profile likelihood one.
  expect_error(mlmr(t6_paragraph_comprehension ~ t7_sentence,
                    data = holzinger_swineford,
                    wald_df = "missing_information"),
               "ci_method = \"profile\" is not one")
  expect_error(mlmr(t6_paragraph_comprehension ~ t7_sentence,
                    data = holzinger_swineford, ci_method = "wald",
                    wald_df = "missing"),
               "or \"missing_information\"")
})

test_that("the missing information rule needs a FIML fit or complete data", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a lavaan fit
  expect_error(mlmr(t26_flags ~ t1_visual_perception + t4_lozenges,
                    data = holzinger_swineford, ci_method = "wald",
                    missing = "pairwise", effect_sizes = FALSE,
                    wald_df = "missing_information"),
               "takes each coefficient's fraction of missing information")
})


test_that("summary() and confint() write the level and the columns whatever getOption(\"digits\") is", {
  skip_if_not_installed("lavaan")
  # format() followed getOption("digits"), so under options(digits = 2) a
  # .975 fit printed "level = 0.97" over columns "1.3 %" and "99 %".
  fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
              conf_level = 0.975, effect_sizes = FALSE)
  withr::local_options(digits = 2)
  out <- utils::capture.output(print(summary(fit)))
  expect_true(any(grepl("Confidence intervals (method: wald, level = 0.975):",
                        out, fixed = TRUE)))
  expect_true(any(grepl("1.25 %", out, fixed = TRUE)))
  expect_true(any(grepl("98.75 %", out, fixed = TRUE)))
  expect_identical(colnames(confint(fit)), c("1.25 %", "98.75 %"))
  # The default level keeps its labels, as stats::confint() writes them.
  fit95 <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald",
                effect_sizes = FALSE)
  expect_identical(colnames(confint(fit95)), c("2.5 %", "97.5 %"))
  expect_identical(colnames(confint(fit95)),
                   colnames(stats::confint(stats::lm(mpg ~ wt + hp,
                                                     data = mtcars))))
  out95 <- utils::capture.output(print(summary(fit95)))
  expect_true(any(grepl("level = 0.95):", out95, fixed = TRUE)))
  # A level below 1 never reads 1.
  expect_identical(.format_conf_prop(c(0.95, 0.975, 1 - 1e-9)),
                   c("0.95", "0.975", "0.999999999"))
})
