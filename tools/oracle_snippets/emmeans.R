# Live oracle comparisons against emmeans, removed from the test suite when the
# anchors were pinned (emmeans 2.0.3, 2026-08-09). Run at release time when the
# oracle packages are installed locally.

## from tests/testthat/test-equivalence_c.R
local({
  set.seed(113)
  d <- data.frame(
    g = factor(rep(c("a", "b", "c"), times = c(20, 25, 30))),
    y = rnorm(75, mean = rep(c(10, 10.5, 12), times = c(20, 25, 30)), sd = 4)
  )
  fit <- lm(y ~ g, data = d)
  emm <- emmeans::emmeans(fit, ~g)
  chk_e <- summary(emmeans::test(emmeans::contrast(emm, list(ba = c(-1, 1, 0))),
                                 delta = 2, side = "equivalence",
                                 adjust = "none"))
  chk_n <- summary(emmeans::test(emmeans::contrast(emm, list(ba = c(-1, 1, 0))),
                                 delta = 2, side = "noninferiority",
                                 adjust = "none"))
  res <- DMAR::equivalence_c(means = as.data.frame(emm)$emmean,
                      s_anova = summary(fit)$sigma,
                      c_weights = c(-1, 1, 0), n = c(20, 25, 30),
                      delta_upper = 2)
  dmar   <- c(res$value[res$term == "p_tost"],
              res$value[res$term == "p_noninferiority"])
  oracle <- c(chk_e$p.value, chk_n$p.value)
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-12)))
})

## from tests/testthat/test-contrast_adjusted.R
local({
  # 2x2 ANCOVA; reference-grid cell order (first factor fastest):
  # (a1,b1),(a2,b1),(a1,b2),(a2,b2).
  set.seed(113)
  d <- data.frame(
    A = factor(rep(c("a1", "a2"), each = 40)),
    B = factor(rep(rep(c("b1", "b2"), each = 20), 2)),
    x = rnorm(80)
  )
  d$y <- 5 + 2 * (d$A == "a2") + 1.5 * (d$B == "b2") + 0.8 * d$x + rnorm(80)
  fit <- lm(y ~ A * B + x, data = d)
  mycontrast <- c(-0.5, 0.5, -0.5, 0.5)
  out <- DMAR::contrast_adjusted(fit, contrast = mycontrast)
  emm <- emmeans::emmeans(fit, ~ A * B)
  ref <- as.data.frame(
    confint(emmeans::contrast(emm, list(mycontrast = mycontrast)))
  )
  dmar   <- c(out$value[out$term == "contrast"],
              out$value[out$term == "lower_limit"],
              out$value[out$term == "upper_limit"])
  oracle <- c(ref$estimate, ref$lower.CL, ref$upper.CL)
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-6)))
})

## from tests/testthat/test-contrast_adjusted.R
local({
  # 2x3 factorial ANCOVA with two covariates; cells in reference-grid order
  # (A fastest within C). Compare mean(c1, c2) versus c3, marginal over A.
  set.seed(113)
  d <- expand.grid(A = c("a1", "a2"), C = c("c1", "c2", "c3"))
  d <- d[rep(seq_len(nrow(d)), each = 25), ]
  d$A <- factor(d$A); d$C <- factor(d$C)
  d$x1 <- rnorm(nrow(d))
  d$x2 <- rnorm(nrow(d))
  d$y <- 10 + 1.2 * (d$A == "a2") + as.numeric(d$C) +
    0.5 * d$x1 - 0.3 * d$x2 + rnorm(nrow(d))
  fit <- lm(y ~ A * C + x1 + x2, data = d)
  w <- c(0.25, 0.25, 0.25, 0.25, -0.5, -0.5)
  out <- DMAR::contrast_adjusted(fit, contrast = w)
  emm <- emmeans::emmeans(fit, ~ A * C)
  ref <- as.data.frame(confint(emmeans::contrast(emm, list(w = w))))
  dmar   <- c(out$value[out$term == "contrast"],
              out$value[out$term == "lower_limit"],
              out$value[out$term == "upper_limit"])
  oracle <- c(ref$estimate, ref$lower.CL, ref$upper.CL)
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-6)))
})

## from tests/testthat/test-contrast_adjusted.R
local({
  # Empty-cell 2x2 ANCOVA (the (a2, b2) cell is unobserved); the contrast
  # confined to the observed cells (B within a1) is estimable and must agree
  # with emmeans.
  set.seed(113)
  d <- do.call(rbind, lapply(list(c("a1", "b1"), c("a2", "b1"), c("a1", "b2")),
    function(cell) data.frame(A = cell[1], B = cell[2])[rep(1, 25), ]))
  d$A <- factor(d$A, levels = c("a1", "a2"))
  d$B <- factor(d$B, levels = c("b1", "b2"))
  d$x <- rnorm(75)
  d$y <- 5 + 2 * (d$A == "a2") + 1.5 * (d$B == "b2") + 0.8 * d$x + rnorm(75)
  fit <- lm(y ~ A * B + x, data = d)
  w   <- c(-1, 0, 1, 0)
  out <- DMAR::contrast_adjusted(fit, w)
  emm <- emmeans::emmeans(fit, ~ A * B)
  ref <- as.data.frame(confint(emmeans::contrast(emm, list(w = w))))
  dmar   <- c(out$value[out$term == "contrast"],
              out$value[out$term == "lower_limit"],
              out$value[out$term == "upper_limit"])
  oracle <- c(ref$estimate, ref$lower.CL, ref$upper.CL)
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-6)))
})
