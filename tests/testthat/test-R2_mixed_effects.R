test_that("R2_mixed_effects matches the Nakagawa reference values for a random-slope model", {
  skip_if_not_installed("lme4")

  fit <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                    data = lme4::sleepstudy)
  res <- R2_mixed_effects(fit)

  # Pinned from performance::r2_nakagawa (performance 0.16.0, 2026-08-09); live comparison in tools/oracle_checks.R.
  expect_equal(res$value[res$term == "R2_marginal"],
               0.2786510581631936, tolerance = 1e-3)
  expect_equal(res$value[res$term == "R2_conditional"],
               0.7992199227797347, tolerance = 1e-3)
})

test_that("R2_mixed_effects matches the Nakagawa reference values for a random-intercept model", {
  skip_if_not_installed("lme4")

  fit <- lme4::lmer(Reaction ~ Days + (1 | Subject),
                    data = lme4::sleepstudy)
  res <- R2_mixed_effects(fit)

  # Pinned from performance::r2_nakagawa (performance 0.16.0, 2026-08-09); live comparison in tools/oracle_checks.R.
  expect_equal(res$value[res$term == "R2_marginal"],
               0.2798856367631553, tolerance = 1e-3)
  expect_equal(res$value[res$term == "R2_conditional"],
               0.7042554523353618, tolerance = 1e-3)
})

test_that("R2_mixed_effects reproduces a hand-derived variance decomposition", {
  skip_if_not_installed("lme4")

  # Random-intercept model: hand-build the three variance components and the
  # two R2 values from VarCorr, sigma, and the fixed-effect linear predictor,
  # then confirm R2_mixed_effects returns the same numbers.
  fit <- lme4::lmer(Reaction ~ Days + (1 | Subject),
                    data = lme4::sleepstudy)

  var_fixed <- stats::var(as.vector(stats::model.matrix(fit) %*% lme4::fixef(fit)))
  vc <- as.data.frame(lme4::VarCorr(fit))
  var_random <- vc$vcov[vc$grp == "Subject" & vc$var1 == "(Intercept)" & is.na(vc$var2)]
  var_resid <- stats::sigma(fit)^2
  total <- var_fixed + var_random + var_resid

  R2_marginal_hand    <- var_fixed / total
  R2_conditional_hand <- (var_fixed + var_random) / total

  res <- R2_mixed_effects(fit)
  expect_equal(res$value[res$term == "R2_marginal"],
               R2_marginal_hand, tolerance = 1e-10)
  expect_equal(res$value[res$term == "R2_conditional"],
               R2_conditional_hand, tolerance = 1e-10)

  # Marginal is never larger than conditional, both in [0, 1].
  expect_lte(R2_marginal_hand, R2_conditional_hand)
  expect_gte(R2_marginal_hand, 0)
  expect_lte(R2_conditional_hand, 1)
})

test_that("R2_mixed_effects returns a dmar_tbl with broom methods", {
  skip_if_not_installed("lme4")

  # Under devtools::load_all() the NAMESPACE in the worktree may not yet carry
  # the S3method() line for these methods (it is regenerated centrally by
  # devtools::document()), so register them here to exercise dispatch.
  registerS3method("tidy", "dmar_R2_mixed_effects", tidy.dmar_R2_mixed_effects,
                   envir = asNamespace("generics"))
  registerS3method("glance", "dmar_R2_mixed_effects", glance.dmar_R2_mixed_effects,
                   envir = asNamespace("generics"))

  fit <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                    data = lme4::sleepstudy)
  res <- R2_mixed_effects(fit)

  expect_s3_class(res, "dmar_tbl")
  expect_s3_class(res, "dmar_R2_mixed_effects")
  expect_true(is.numeric(res$value))
  expect_setequal(res$term, c("R2_marginal", "R2_conditional"))

  td <- generics::tidy(res)
  expect_equal(td$term, c("R2_marginal", "R2_conditional"))
  expect_true(all(c("term", "estimate") %in% names(td)))
  expect_equal(td$estimate[td$term == "R2_marginal"],
               res$value[res$term == "R2_marginal"])

  gl <- generics::glance(res)
  expect_equal(nrow(gl), 2L)
})

test_that("R2_mixed_effects matches the Nakagawa reference values for an nlme::lme fit", {
  skip_if_not_installed("nlme")
  skip_if_not_installed("lme4")

  fit_lme <- nlme::lme(Reaction ~ Days, random = ~ 1 | Subject,
                       data = lme4::sleepstudy)

  res_lme <- R2_mixed_effects(fit_lme)

  # Pinned from performance::r2_nakagawa on the equivalent lme4 REML fit of
  # Reaction ~ Days + (1 | Subject) (performance 0.16.0, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  expect_equal(res_lme$value[res_lme$term == "R2_marginal"],
               0.2798856367631553, tolerance = 1e-2)
  expect_equal(res_lme$value[res_lme$term == "R2_conditional"],
               0.7042554523353618, tolerance = 1e-2)
})

test_that("R2_mixed_effects handles separate random intercept and slope terms", {
  skip_if_not_installed("lme4")

  # CRITICAL-02 regression. For (1 | g) + (0 + x | g), VarCorr names the two
  # components 'Subject' and 'Subject.1' while getME(., "cnms") repeats
  # 'Subject'; indexing VarCorr by the repeated name used the intercept
  # covariance for both terms and dropped the slope entirely. The random-effect
  # design matrices and covariance components must instead be aligned by
  # position. Oracle is the direct Johnson (2014) mean quadratic form built by
  # hand: for uncorrelated terms, var_random = tau0^2 + tau1^2 * mean(Days^2).
  # performance/insight is deliberately not used as the oracle here: it silently
  # omits the random-slope variance for this split specification.
  ss  <- lme4::sleepstudy
  fit <- lme4::lmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
                    data = ss, REML = TRUE)

  vcd  <- as.data.frame(lme4::VarCorr(fit))
  # which() so the NA var1 on the Residual row does not leak an NA into the
  # logical subset.
  tau0 <- vcd$vcov[which(vcd$var1 == "(Intercept)" & is.na(vcd$var2))]
  tau1 <- vcd$vcov[which(vcd$var1 == "Days" & is.na(vcd$var2))]
  var_fixed  <- stats::var(as.vector(stats::model.matrix(fit) %*% lme4::fixef(fit)))
  var_resid  <- stats::sigma(fit)^2
  var_random <- tau0 * mean(1^2) + tau1 * mean(ss$Days^2)
  total      <- var_fixed + var_random + var_resid
  marg_hand  <- var_fixed / total
  cond_hand  <- (var_fixed + var_random) / total

  res  <- R2_mixed_effects(fit)
  marg <- res$value[res$term == "R2_marginal"]
  cond <- res$value[res$term == "R2_conditional"]
  expect_equal(marg, marg_hand, tolerance = 1e-10)
  expect_equal(cond, cond_hand, tolerance = 1e-10)

  # Guard against the specific defect: the intercept-twice bug returned a
  # marginal near 0.045 and a conditional near 0.967.
  expect_gt(marg, 0.20)
  expect_lt(cond, 0.90)

  # Independent implementation cross-check.
  # Pinned from r2mlm::r2mlm (r2mlm 0.3.8, 2026-08-09); live comparison in tools/oracle_checks.R.
  expect_equal(marg, 0.2828350562331848, tolerance = 5e-3)
  expect_equal(cond, 0.7966272707001041, tolerance = 5e-3)
})

test_that("R2_mixed_effects treats the double-bar form as the split form", {
  skip_if_not_installed("lme4")

  # (Days || Subject) is lme4 shorthand for (1 | Subject) + (0 + Days | Subject),
  # so it hits the same duplicate-grouping-name path and must give the same R2.
  ss <- lme4::sleepstudy
  fit_bar   <- lme4::lmer(Reaction ~ Days + (Days || Subject),
                          data = ss, REML = TRUE)
  fit_split <- lme4::lmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
                          data = ss, REML = TRUE)
  expect_equal(R2_mixed_effects(fit_bar)$value, R2_mixed_effects(fit_split)$value,
               tolerance = 1e-8)
})
