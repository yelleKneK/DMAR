## from tests/testthat/test-ci_dunnett.R
## ci_dunnett() matches multcomp::glht() on PlantGrowth (two-sided limits).
local({
  fit <- lm(weight ~ group, data = PlantGrowth)
  res <- DMAR::ci_dunnett(fit, control = "ctrl")
  g  <- multcomp::glht(fit,
                       linfct = multcomp::mcp(group = "Dunnett"))
  set.seed(113)  # multcomp's quantile is stochastic; fix its RNG state
  ci <- confint(g)$confint
  # The estimated mean differences are the ordinary contrasts of group
  # means and must reproduce multcomp's estimates exactly.
  stopifnot(isTRUE(all.equal(res$mean_difference, unname(ci[, "Estimate"]),
                             tolerance = 1e-8)))
  # The simultaneous Dunnett limits use the same multivariate-t critical
  # value, so the interval endpoints agree with glht()'s confint().
  stopifnot(isTRUE(all.equal(res$lower_limit, unname(ci[, "lwr"]),
                             tolerance = 1e-4)))
  stopifnot(isTRUE(all.equal(res$upper_limit, unname(ci[, "upr"]),
                             tolerance = 1e-4)))
})

## from tests/testthat/test-ci_dunnett.R
## ci_dunnett() adjusted p-values match multcomp::glht().
local({
  fit <- lm(weight ~ group, data = PlantGrowth)
  res <- DMAR::ci_dunnett(fit, control = "ctrl")
  g   <- multcomp::glht(fit, linfct = multcomp::mcp(group = "Dunnett"))
  set.seed(113)  # multcomp's adjusted p-values are stochastic; fix its RNG state
  p_glht <- summary(g)$test$pvalues
  stopifnot(isTRUE(all.equal(res$p_adjusted, as.numeric(p_glht),
                             tolerance = 1e-3)))
})

## from tests/testthat/test-ci_dunnett.R
## ci_dunnett() one-sided limits match multcomp::glht().
local({
  fit <- lm(weight ~ group, data = PlantGrowth)
  r_l <- DMAR::ci_dunnett(fit, control = "ctrl", alternative = "less")
  r_g <- DMAR::ci_dunnett(fit, control = "ctrl", alternative = "greater")
  set.seed(113)  # multcomp's quantile is stochastic; fix its RNG state
  ci_l <- confint(multcomp::glht(fit,
                                 linfct = multcomp::mcp(group = "Dunnett"),
                                 alternative = "less"))$confint
  ci_g <- confint(multcomp::glht(fit,
                                 linfct = multcomp::mcp(group = "Dunnett"),
                                 alternative = "greater"))$confint
  # Compare on the absolute scale: multcomp's simulated critical value carries
  # Monte Carlo error of a few units in the fourth decimal, which a relative
  # tolerance would magnify for the near-zero trt2 "greater" limit.
  stopifnot(max(abs(r_l$upper_limit - unname(ci_l[, "upr"]))) < 5e-4)
  stopifnot(max(abs(r_g$lower_limit - unname(ci_g[, "lwr"]))) < 5e-4)
})
