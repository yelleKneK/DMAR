test_that("contrast_adjusted() returns the documented five-row schema", {
  set.seed(113)
  d <- data.frame(
    A = factor(rep(c("a1", "a2"), each = 40)),
    B = factor(rep(rep(c("b1", "b2"), each = 20), 2)),
    x = rnorm(80)
  )
  d$y <- 5 + 2 * (d$A == "a2") + 1.5 * (d$B == "b2") + 0.8 * d$x + rnorm(80)
  fit <- lm(y ~ A * B + x, data = d)

  out <- contrast_adjusted(fit, contrast = c(-0.5, 0.5, -0.5, 0.5))
  expect_s3_class(out, "dmar_tbl")
  expect_identical(out$term, c("contrast", "lower_limit", "upper_limit", "t", "p"))
  expect_true(is.numeric(out$value))
  expect_identical(attr(out, "conf_level"), 0.95)
  expect_identical(attr(out, "p_terms"), "p")
  # The interval brackets the point estimate.
  expect_lt(out$value[2], out$value[1])
  expect_gt(out$value[3], out$value[1])
})

test_that("contrast_adjusted() matches pinned emmeans values to 1e-6", {
  set.seed(113)
  d <- data.frame(
    A = factor(rep(c("a1", "a2"), each = 40)),
    B = factor(rep(rep(c("b1", "b2"), each = 20), 2)),
    x = rnorm(80)
  )
  d$y <- 5 + 2 * (d$A == "a2") + 1.5 * (d$B == "b2") + 0.8 * d$x + rnorm(80)
  fit <- lm(y ~ A * B + x, data = d)

  # Reference-grid cell order (first factor fastest): (a1,b1),(a2,b1),(a1,b2),(a2,b2).
  # This is also the order emmeans::emmeans(fit, ~ A * B) enumerates its grid.
  mycontrast <- c(-0.5, 0.5, -0.5, 0.5)
  out <- contrast_adjusted(fit, contrast = mycontrast)

  # Pinned from emmeans::contrast and confint (emmeans 2.0.3, 2026-08-09);
  # live comparison in tools/oracle_checks.R.
  expect_equal(out$value[out$term == "contrast"],    1.901554725395017, tolerance = 1e-6)
  expect_equal(out$value[out$term == "lower_limit"], 1.488955687103666, tolerance = 1e-6)
  expect_equal(out$value[out$term == "upper_limit"], 2.314153763686368, tolerance = 1e-6)
})

test_that("a more elaborate contrast also matches pinned emmeans values", {
  set.seed(113)
  # 2 x 3 factorial ANCOVA with two covariates.
  d <- expand.grid(A = c("a1", "a2"), C = c("c1", "c2", "c3"))
  d <- d[rep(seq_len(nrow(d)), each = 25), ]
  d$A <- factor(d$A); d$C <- factor(d$C)
  d$x1 <- rnorm(nrow(d))
  d$x2 <- rnorm(nrow(d))
  d$y <- 10 + 1.2 * (d$A == "a2") + as.numeric(d$C) +
    0.5 * d$x1 - 0.3 * d$x2 + rnorm(nrow(d))
  fit <- lm(y ~ A * C + x1 + x2, data = d)

  # Cells in reference-grid order (A fastest within C):
  # (a1,c1),(a2,c1),(a1,c2),(a2,c2),(a1,c3),(a2,c3).
  # Compare mean(c1, c2) versus c3, marginal over A.
  w <- c(0.25, 0.25, 0.25, 0.25, -0.5, -0.5)
  out <- contrast_adjusted(fit, contrast = w)

  # Pinned from emmeans::contrast and confint (emmeans 2.0.3, 2026-08-09);
  # live comparison in tools/oracle_checks.R.
  expect_equal(out$value[out$term == "contrast"],    -1.621321974083438, tolerance = 1e-6)
  expect_equal(out$value[out$term == "lower_limit"], -1.95610771191861,  tolerance = 1e-6)
  expect_equal(out$value[out$term == "upper_limit"], -1.286536236248265, tolerance = 1e-6)
})

test_that("cell-means and crossed parameterizations agree", {
  set.seed(113)
  d <- data.frame(
    A = factor(rep(c("a1", "a2"), each = 40)),
    B = factor(rep(rep(c("b1", "b2"), each = 20), 2)),
    x = rnorm(80)
  )
  d$y <- 5 + 2 * (d$A == "a2") + 1.5 * (d$B == "b2") + 0.8 * d$x + rnorm(80)
  w <- c(-0.5, 0.5, -0.5, 0.5)

  fit_crossed <- lm(y ~ A * B + x, data = d)
  d$cell <- interaction(d$A, d$B)          # levels: a1.b1, a2.b1, a1.b2, a2.b2
  fit_cm <- lm(y ~ 0 + cell + x, data = d)

  r1 <- contrast_adjusted(fit_crossed, contrast = w)
  r2 <- contrast_adjusted(fit_cm, contrast = w)
  expect_equal(r1$value, r2$value, tolerance = 1e-8)
})

test_that("input validation errors are informative", {
  set.seed(113)
  d <- data.frame(g = factor(rep(c("g1", "g2"), each = 20)), x = rnorm(40))
  d$y <- d$x + as.numeric(d$g) + rnorm(40)
  fit <- lm(y ~ g + x, data = d)

  expect_error(contrast_adjusted(fit, contrast = c(1, -1, 0)), "one weight per cell")
  expect_error(contrast_adjusted(fit, contrast = c(1, -1), conf_level = 2),
               "strictly between 0 and 1")
  # A model with no covariate is out of scope for an adjusted-means contrast.
  fit_no_cov <- lm(y ~ g, data = d)
  expect_error(contrast_adjusted(fit_no_cov, contrast = c(1, -1)), "covariate")
})

test_that("contrast_adjusted() rejects a non-estimable empty-cell contrast", {
  # CRITICAL-04 regression. A 2x2 ANCOVA with the (a2, b2) cell empty makes lm
  # alias the Aa2:Bb2 interaction (coefficient NA). emmeans reports the (a2, b2)
  # adjusted mean as non-estimable; contrast_adjusted() previously dropped the
  # aliased column and returned a fabricated, highly significant value. Any
  # contrast that places weight on the empty cell must now error, while a
  # contrast confined to the observed cells still matches emmeans exactly.
  set.seed(113)
  d <- do.call(rbind, lapply(list(c("a1", "b1"), c("a2", "b1"), c("a1", "b2")),
    function(cell) data.frame(A = cell[1], B = cell[2])[rep(1, 25), ]))
  d$A <- factor(d$A, levels = c("a1", "a2"))
  d$B <- factor(d$B, levels = c("b1", "b2"))
  d$x <- rnorm(75)
  d$y <- 5 + 2 * (d$A == "a2") + 1.5 * (d$B == "b2") + 0.8 * d$x + rnorm(75)
  fit <- lm(y ~ A * B + x, data = d)
  expect_true(anyNA(stats::coef(fit)))   # confirm the empty cell aliased a term

  # Reference-grid order (A fastest): (a1,b1), (a2,b1), (a1,b2), (a2,b2).
  # Isolating the empty (a2,b2) cell mean, and the 2x2 interaction that involves
  # it, are both non-estimable.
  expect_error(contrast_adjusted(fit, c(0, 0, 0, 1)), "not estimable")
  expect_error(contrast_adjusted(fit, c(1, -1, -1, 1)), "not estimable")

  # A contrast confined to observed cells (B within a1) is estimable and must
  # still agree with emmeans. Pinned from emmeans::contrast and confint
  # (emmeans 2.0.3, 2026-08-09); live comparison in tools/oracle_checks.R.
  w   <- c(-1, 0, 1, 0)
  out <- contrast_adjusted(fit, w)
  expect_equal(out$value[out$term == "contrast"],    1.105756219364711,  tolerance = 1e-6)
  expect_equal(out$value[out$term == "lower_limit"], 0.5797769453566087, tolerance = 1e-6)
  expect_equal(out$value[out$term == "upper_limit"], 1.631735493372813,  tolerance = 1e-6)
})
