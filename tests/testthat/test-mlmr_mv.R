test_that("mlmr_mv() returns an mlmr_mv object with the expected slots", {
  skip_if_not_installed("lavaan")
  fit <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = mtcars,
                 ci_method = "wald", effect_sizes = FALSE)
  expect_s3_class(fit, "mlmr_mv")
  expect_true(is.matrix(fit$coefficients))
  expect_equal(dim(fit$coefficients), c(3L, 2L))
  expect_equal(colnames(fit$coefficients), c("mpg", "disp"))
  expect_equal(rownames(fit$coefficients), c("(Intercept)", "wt", "hp"))
  expect_true(is.matrix(fit$residual_cov))
  expect_equal(dim(fit$residual_cov), c(2L, 2L))
  expect_s3_class(fit$coef_table, "data.frame")
  expect_true(all(c("outcome", "term", "estimate", "se", "z_value",
                    "p_value", "ci_lower", "ci_upper") %in%
                  colnames(fit$coef_table)))
})


test_that("vcov.mlmr_mv is coefficient-only and aligned with coef()", {
  skip_if_not_installed("lavaan")
  fit <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = mtcars,
                 ci_method = "wald", effect_sizes = FALSE)
  V <- vcov(fit)
  cf <- coef(fit)

  # Coefficient-only: 3 terms times 2 outcomes, not the full lavaan
  # parameter covariance, which also carries the residual variances
  # and covariances and, with fixed_x = FALSE, the predictor moments.
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(6L, 6L))
  expect_lt(nrow(V), length(lavaan::coef(fit$lavaan_fit)))

  # Dimnames are outcome-major "outcome:term", the column-major
  # flattening of the coef() matrix, so a Wald quadratic form built
  # from as.vector(coef(fit)) and vcov(fit) lines up element by
  # element. The order is also coef_table's row order.
  coef_names <- paste0(rep(colnames(cf), each = nrow(cf)), ":",
                       rownames(cf))
  expect_identical(rownames(V), coef_names)
  expect_identical(colnames(V), coef_names)
  expect_identical(rownames(V),
                   paste0(fit$coef_table$outcome, ":",
                          fit$coef_table$term))

  # The diagonal is the squared SEs the fit table reports.
  expect_equal(unname(diag(V)), fit$coef_table$se^2, tolerance = 1e-10)

  # Same naming scheme as base R's vcov on the OLS multivariate fit.
  V_lm <- vcov(lm(cbind(mpg, disp) ~ wt + hp, data = mtcars))
  expect_identical(rownames(V), rownames(V_lm))
})


test_that("mlmr_mv() point estimates match lm() multivariate to working precision", {
  skip_if_not_installed("lavaan")
  fit_m <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = mtcars,
                   ci_method = "wald", effect_sizes = FALSE)
  fit_l <- lm(cbind(mpg, disp) ~ wt + hp, data = mtcars)
  expect_equal(unname(coef(fit_m)), unname(coef(fit_l)),
               tolerance = 1e-4)
})


test_that("mlmr_mv() errors on a single-outcome formula", {
  skip_if_not_installed("lavaan")
  expect_error(mlmr_mv(mpg ~ wt + hp, data = mtcars),
               "cbind\\(\\) on the left-hand side")
})


test_that("mlmr_mv() per-outcome R^2 matches lm() per-outcome R^2", {
  skip_if_not_installed("lavaan")
  fit_m <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = mtcars,
                   ci_method = "wald", effect_sizes = FALSE)
  fit_lm_mpg  <- lm(mpg ~ wt + hp, data = mtcars)
  fit_lm_disp <- lm(disp ~ wt + hp, data = mtcars)
  expect_equal(unname(fit_m$R2["mpg"]),
               summary(fit_lm_mpg)$r.squared, tolerance = 1e-6)
  expect_equal(unname(fit_m$R2["disp"]),
               summary(fit_lm_disp)$r.squared, tolerance = 1e-6)
})


test_that("mlmr_mv() FIML uses more rows than listwise when one outcome is partly missing", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two lavaan refits; the lm() anchors above run on CRAN
  set.seed(113)
  d <- mtcars
  d$disp[sample.int(nrow(d), 6)] <- NA
  fit_fiml <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = d,
                      ci_method = "wald", effect_sizes = FALSE)
  fit_lwd  <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = d,
                      missing = "listwise", ci_method = "wald",
                      effect_sizes = FALSE)
  expect_gt(nobs(fit_fiml), nobs(fit_lwd))
})


test_that("mlmr_mv() profile LR CIs are produced and bracket the estimate", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = mtcars,
                 effect_sizes = FALSE)
  ci_df <- confint(fit)
  expect_s3_class(ci_df, "data.frame")
  expect_equal(nrow(ci_df), 6L)  # 3 terms * 2 outcomes
  lo <- ci_df[[3L]]
  hi <- ci_df[[4L]]
  expect_true(all(is.finite(lo)))
  expect_true(all(is.finite(hi)))
  expect_true(all(lo < hi))
})


test_that("mlmr_mv() bootstrap CIs are reproducible with boot_seed", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # suppressWarnings: with B = 50 on n = 32, a few bootstrap resamples trip
  # lavaan's post-fit checks; reproducibility, not fit quality, is under test.
  fit_a <- suppressWarnings(mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = mtcars,
                   ci_method = "boot", B = 50, boot_seed = 113,
                   effect_sizes = FALSE))
  fit_b <- suppressWarnings(mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = mtcars,
                   ci_method = "boot", B = 50, boot_seed = 113,
                   effect_sizes = FALSE))
  expect_equal(fit_a$ci, fit_b$ci, tolerance = 1e-8)
})


test_that("mlmr_mv() predict() returns a matrix shaped (n, J)", {
  skip_if_not_installed("lavaan")
  fit <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = mtcars,
                 ci_method = "wald", effect_sizes = FALSE)
  new_d <- data.frame(wt = c(2, 3, 4), hp = c(80, 150, 220))
  pred <- predict(fit, newdata = new_d)
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(3L, 2L))
  expect_equal(colnames(pred), c("mpg", "disp"))
})


test_that("mlmr_mv() effect sizes are produced when requested", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # the effect sizes cost a refit per predictor and outcome
  fit <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = mtcars,
                 ci_method = "wald", effect_sizes = TRUE)
  es <- fit$effect_sizes
  expect_s3_class(es, "data.frame")
  expect_named(es, c("outcome", "term", "sr2", "f2"))
  # 2 predictors * 2 outcomes = 4 rows
  expect_equal(nrow(es), 4L)
})

test_that("tidy() and glance() work on mlmr_mv fits as documented", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a fit plus a refit per predictor and outcome
  fit <- mlmr_mv(cbind(mpg, hp) ~ wt, data = mtcars, ci_method = "wald")
  td <- generics::tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_identical(names(td), c("response", "term", "estimate", "se",
                                "statistic", "p_value"))
  expect_setequal(unique(td$response), unique(fit$coef_table$outcome))
  expect_equal(nrow(td), nrow(fit$coef_table))
  td_ci <- generics::tidy(fit, conf.int = TRUE, standardized = TRUE)
  expect_true(all(c("ci_lower", "ci_upper", "std_estimate") %in%
                    names(td_ci)))
  gl <- generics::glance(fit)
  expect_equal(nrow(gl), 1)
  expect_equal(gl$n.responses, 2)
  expect_equal(gl$R2, mean(fit$R2))
})


test_that("mlmr_mv auxiliary leaves complete-data coefficients unchanged", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  d <- mtcars
  d$aux <- d$qsec + rnorm(nrow(d))
  f0 <- mlmr_mv(cbind(mpg, drat) ~ wt + hp, data = d, ci_method = "wald",
                effect_sizes = FALSE)
  fa <- mlmr_mv(cbind(mpg, drat) ~ wt + hp, data = d, ci_method = "wald",
                effect_sizes = FALSE, auxiliary = "aux")
  expect_equal(unname(coef(f0)), unname(coef(fa)), tolerance = 1e-4)
  expect_identical(fa$auxiliary, "aux")
  expect_null(f0$auxiliary)
})


test_that("mlmr_mv auxiliary validates inputs", {
  skip_if_not_installed("lavaan")
  expect_error(mlmr_mv(cbind(mpg, drat) ~ wt, data = mtcars,
                       auxiliary = "nope"), "not found in 'data'")
  expect_error(mlmr_mv(cbind(mpg, drat) ~ wt, data = mtcars,
                       auxiliary = "mpg"), "also appear in 'formula'")
})


test_that("anova.mlmr_mv errors on two fits of equal N run on different data", {
  # HIGH-06 regression for the multivariate path. Two independently
  # generated datasets of the same row count are, by construction, not the
  # same data; anova() must refuse them rather than return a bogus
  # chi square, even though their complete-case counts match.
  skip_if_not_installed("lavaan")
  set.seed(113)
  n <- 120L
  d1 <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d1$y <- 10 + 3 * d1$x1 + rnorm(n)
  d1$z <- 2 - d1$x2 + rnorm(n)
  d2 <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d2$y <- -5 - 0.6 * d2$x1 + rnorm(n)
  d2$z <- 1 + d2$x2 + rnorm(n)
  m1 <- mlmr_mv(cbind(y, z) ~ x1 + x2, data = d1, ci_method = "wald",
                effect_sizes = FALSE)
  m2 <- mlmr_mv(cbind(y, z) ~ x1,      data = d2, ci_method = "wald",
                effect_sizes = FALSE)
  expect_identical(m1$N_complete, m2$N_complete)   # equal N; the trap
  expect_error(anova(m1, m2), "same data")
})


test_that("anova.mlmr_mv same-data LRT matches the multivariate ML oracle", {
  # HIGH-06 companion for mlmr_mv: the same-data fix must leave the correct
  # LR statistic intact. Independent oracle: for a joint-normal
  # multivariate FIML regression with fixed.x = FALSE, the LR statistic for
  # nested fits is the Wilks statistic G^2 = N * log(det(E_reduced) /
  # det(E_full)), with E the residual sum-of-squares-and-cross-products
  # matrix from ordinary least squares on the complete data.
  skip_if_not_installed("lavaan")
  fit_small <- mlmr_mv(cbind(mpg, disp) ~ wt,      data = mtcars,
                       ci_method = "wald", effect_sizes = FALSE)
  fit_large <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = mtcars,
                       ci_method = "wald", effect_sizes = FALSE)
  res <- anova(fit_small, fit_large)
  chisq_diff <- res[["Chisq diff"]][2L]

  lm_small <- lm(cbind(mpg, disp) ~ wt,      data = mtcars)
  lm_large <- lm(cbind(mpg, disp) ~ wt + hp, data = mtcars)
  N <- nrow(mtcars)
  E_small <- crossprod(residuals(lm_small))
  E_large <- crossprod(residuals(lm_large))
  oracle <- N * log(det(E_small) / det(E_large))

  expect_equal(chisq_diff, oracle, tolerance = 1e-6)
})


test_that("anova.mlmr_mv accepts nested FIML fits whose complete-case counts differ", {
  # Multivariate counterpart of the mlmr test of the same name. The
  # same-data guard must compare the analysis data, not a complete-case
  # count: adding a partially observed predictor changes the number of
  # rows complete on all modeled variables without changing which
  # observations were analyzed.
  #
  # Oracle: the observed data normal log-likelihood for the four
  # variables (y1, y2, x1, x2), written out directly and maximized
  # outside DMAR and outside lavaan.
  #
  #   Larger model. cbind(y1, y2) ~ x1 + x2 with random predictors and a
  #   free residual covariance is the saturated 4-variate normal, whose
  #   MLE is closed form under this monotone pattern (mean vector and
  #   divisor-N covariance of the complete (y1, y2, x1) block from all
  #   200 rows, regression of x2 on that block from the 160 rows where x2
  #   is observed): logLik = -1088.0941228199.
  #
  #   Smaller model. The same likelihood with both x2 slopes fixed at
  #   zero, parameterized by the joint normal of (x1, x2), the intercepts
  #   and x1 slopes, and a free residual covariance (Nelder-Mead followed
  #   by BFGS restarts): logLik = -1116.8440745476.
  #
  # LR = 2 * (-1088.0941228199 - (-1116.8440745476)) = 57.4999034553 on
  # 2 degrees of freedom, one per outcome.
  skip_if_not_installed("lavaan")
  set.seed(113)
  n  <- 200L
  x1 <- rnorm(n)
  x2 <- 0.4 * x1 + rnorm(n)
  y1 <- 0.5 * x1 + 0.3 * x2 + rnorm(n)
  y2 <- -0.2 * x1 + 0.6 * x2 + rnorm(n)
  d  <- data.frame(y1 = y1, y2 = y2, x1 = x1, x2 = x2)
  d$x2[sample(n, 40L)] <- NA

  m1 <- mlmr_mv(cbind(y1, y2) ~ x1,      data = d, missing = "fiml",
                ci_method = "wald", effect_sizes = FALSE)
  m2 <- mlmr_mv(cbind(y1, y2) ~ x1 + x2, data = d, missing = "fiml",
                ci_method = "wald", effect_sizes = FALSE)

  # The condition the old guard tripped on: same data, different
  # complete-case counts.
  expect_identical(m1$N_complete, 200L)
  expect_identical(m2$N_complete, 160L)

  res <- anova(m1, m2)
  expect_s3_class(res, "anova")
  expect_equal(res[["Chisq diff"]][2L], 57.4999034553, tolerance = 1e-7)
  expect_identical(as.integer(res[["Df diff"]][2L]), 2L)
})


test_that("anova.mlmr_mv still refuses nested FIML fits on genuinely different data", {
  # Companion to the test above. Both data sets have 200 rows and 40
  # missing values on x2, so every count matches and only the values
  # differ. Oracle: the two data frames are independently generated, so
  # by construction they are not the same data.
  skip_if_not_installed("lavaan")
  skip_on_cran()  # four FIML refits on simulated data
  make_d <- function(seed) {
    set.seed(seed)
    n  <- 200L
    x1 <- rnorm(n)
    x2 <- 0.4 * x1 + rnorm(n)
    d  <- data.frame(y1 = 0.5 * x1 + 0.3 * x2 + rnorm(n),
                     y2 = -0.2 * x1 + 0.6 * x2 + rnorm(n),
                     x1 = x1, x2 = x2)
    d$x2[sample(n, 40L)] <- NA
    d
  }
  d1 <- make_d(113)
  d2 <- make_d(114)
  m1 <- mlmr_mv(cbind(y1, y2) ~ x1,      data = d1, missing = "fiml",
                ci_method = "wald", effect_sizes = FALSE)
  m2 <- mlmr_mv(cbind(y1, y2) ~ x1 + x2, data = d2, missing = "fiml",
                ci_method = "wald", effect_sizes = FALSE)
  expect_identical(m1$N_complete, 200L)
  expect_identical(m2$N_complete, 160L)
  expect_error(anova(m1, m2), "same data")

  # A row-count mismatch is reported as such.
  m3 <- mlmr_mv(cbind(y1, y2) ~ x1, data = d1[1:150, ], missing = "fiml",
                ci_method = "wald", effect_sizes = FALSE)
  m4 <- mlmr_mv(cbind(y1, y2) ~ x1 + x2, data = d1, missing = "fiml",
                ci_method = "wald", effect_sizes = FALSE)
  expect_error(anova(m3, m4), "150 observations")
})
