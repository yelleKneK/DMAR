## from tests/testthat/test-R2_mixed_effects.R
local({
  # Nakagawa marginal and conditional R2 for a random-slope model.
  fit <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                    data = lme4::sleepstudy)
  res <- DMAR::R2_mixed_effects(fit)
  ref <- performance::r2_nakagawa(fit)
  dmar   <- c(res$value[res$term == "R2_marginal"],
              res$value[res$term == "R2_conditional"])
  oracle <- c(as.numeric(ref$R2_marginal), as.numeric(ref$R2_conditional))
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-3)))
})

## from tests/testthat/test-R2_mixed_effects.R
local({
  # Nakagawa marginal and conditional R2 for a random-intercept model.
  fit <- lme4::lmer(Reaction ~ Days + (1 | Subject),
                    data = lme4::sleepstudy)
  res <- DMAR::R2_mixed_effects(fit)
  ref <- performance::r2_nakagawa(fit)
  dmar   <- c(res$value[res$term == "R2_marginal"],
              res$value[res$term == "R2_conditional"])
  oracle <- c(as.numeric(ref$R2_marginal), as.numeric(ref$R2_conditional))
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-3)))
})

## from tests/testthat/test-R2_mixed_effects.R
local({
  # nlme::lme fit against performance on the equivalent lme4 REML fit.
  fit_lme  <- nlme::lme(Reaction ~ Days, random = ~ 1 | Subject,
                        data = lme4::sleepstudy)
  fit_lmer <- lme4::lmer(Reaction ~ Days + (1 | Subject),
                         data = lme4::sleepstudy, REML = TRUE)
  res_lme <- DMAR::R2_mixed_effects(fit_lme)
  ref <- performance::r2_nakagawa(fit_lmer)
  dmar   <- c(res_lme$value[res_lme$term == "R2_marginal"],
              res_lme$value[res_lme$term == "R2_conditional"])
  oracle <- c(as.numeric(ref$R2_marginal), as.numeric(ref$R2_conditional))
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-2)))
})

## from tests/testthat/test-R2_mixed_effects.R
local({
  # Independent implementation cross-check for separate random intercept and
  # slope terms: r2mlm's f (marginal) and fvm (conditional) totals.
  ss  <- lme4::sleepstudy
  fit <- lme4::lmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
                    data = ss, REML = TRUE)
  res <- DMAR::R2_mixed_effects(fit)
  # as.numeric() drops the noquote class r2mlm puts on $R2s.
  rr  <- r2mlm::r2mlm(fit, bargraph = FALSE)$R2s
  dmar   <- c(res$value[res$term == "R2_marginal"],
              res$value[res$term == "R2_conditional"])
  oracle <- as.numeric(c(rr[1], rr[length(rr)]))
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 5e-3)))
})

## from tests/testthat/test-R2_mixed_effects_decomposition.R
local({
  # Rights and Sterba (2019) totals against r2mlm on a
  # non-cluster-mean-centered model.
  fit <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                    data = lme4::sleepstudy)
  d <- DMAR::R2_mixed_effects_decomposition(fit)
  v <- d$value; names(v) <- d$term
  # as.numeric() drops the noquote class r2mlm puts on $R2s.
  ref <- suppressWarnings(r2mlm::r2mlm(fit, bargraph = FALSE))$R2s[, "total"]
  dmar   <- unname(v[c("total_f", "total_v", "total_m", "total_fv", "total_fvm")])
  oracle <- as.numeric(ref[c("f", "v", "m", "fv", "fvm")])
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-8)))
})

## from tests/testthat/test-R2_mixed_effects_decomposition.R
local({
  # All 12 measures against r2mlm for a cluster-mean-centered model.
  set.seed(113)
  J <- 60; nj <- 15; N <- J * nj; cl <- rep(1:J, each = nj)
  x <- rnorm(N); xbar <- ave(x, cl); xw <- x - xbar; w <- xbar
  u0 <- rnorm(J, 0, 1.2)[cl]; u1 <- rnorm(J, 0, .5)[cl]
  y <- 2 + 0.5 * xw + 0.8 * w + u0 + u1 * xw + rnorm(N, 0, 1)
  dd <- data.frame(y, xw, w, cl = factor(cl))
  fit <- suppressWarnings(
    lme4::lmer(y ~ xw + w + (xw | cl), data = dd))
  d <- DMAR::R2_mixed_effects_decomposition(fit)
  v <- d$value; names(v) <- d$term
  O <- suppressWarnings(r2mlm::r2mlm(fit, bargraph = FALSE))$R2s
  ref <- c(total_f1 = O["f1", "total"], total_f2 = O["f2", "total"],
           total_v = O["v", "total"], total_m = O["m", "total"],
           total_f = O["f", "total"], total_fv = O["fv", "total"],
           total_fvm = O["fvm", "total"], within_f1 = O["f1", "within"],
           within_v = O["v", "within"], within_fv = O["fv", "within"],
           between_f2 = O["f2", "between"], between_m = O["m", "between"])
  dmar   <- unname(v[names(ref)])
  # as.numeric() drops the noquote class r2mlm puts on $R2s.
  oracle <- as.numeric(ref)
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-8)))
})

## from tests/testthat/test-R2_mixed_effects_decomposition.R
local({
  # total_f / total_fvm coincide with the Nakagawa marginal / conditional R2.
  fit <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
  d <- DMAR::R2_mixed_effects_decomposition(fit)
  v <- d$value; names(v) <- d$term
  nk <- performance::r2_nakagawa(fit)
  dmar   <- unname(v[c("total_f", "total_fvm")])
  oracle <- c(as.numeric(nk$R2_marginal), as.numeric(nk$R2_conditional))
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-6)))
})

## from tests/testthat/test-R2_mixed_effects_decomposition.R
local({
  # Separate random intercept and slope terms reproduce the r2mlm totals.
  fit <- lme4::lmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
                    data = lme4::sleepstudy, REML = TRUE)
  d <- DMAR::R2_mixed_effects_decomposition(fit)
  v <- d$value; names(v) <- d$term
  # as.numeric() drops the noquote class r2mlm puts on $R2s.
  ref <- suppressWarnings(r2mlm::r2mlm(fit, bargraph = FALSE))$R2s[, "total"]
  dmar   <- unname(v[c("total_f", "total_v", "total_m", "total_fv", "total_fvm")])
  oracle <- as.numeric(ref[c("f", "v", "m", "fv", "fvm")])
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-8)))
})
