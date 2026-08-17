# Tests for the unified dmar_tidiers topic (R/dmar_tidiers.R), which holds the
# tidy() / glance() methods for the confidence interval family, the post hoc
# family, and the power based sample size planners.
#
# The consolidation must not change which methods are registered or what any of
# them returns, so every expected value below comes from an external oracle
# computed with base R, never from DMAR's own output.

# The @exportS3Method tags register these methods in NAMESPACE once
# devtools::document() has run. Under a bare load_all() in an isolated worktree
# the NAMESPACE may lag the source, so register them here so generics::tidy() /
# generics::glance() dispatch during the test. This is a no-op once the
# regenerated NAMESPACE carries the S3method lines.
.dmar_tidier_classes <- c("dmar_ci_long", "dmar_ci_anova", "dmar_post_hoc_ci",
                          "dmar_ss_power", "dmar_ss_power_sensitivity")
for (.gen in c("tidy", "glance")) {
  for (.cls in .dmar_tidier_classes) {
    .nm <- paste0(.gen, ".", .cls)
    if (exists(.nm))
      registerS3method(.gen, .cls, get(.nm), envir = asNamespace("generics"))
  }
}


test_that("every tidier the consolidation moved is still registered", {
  is_registered <- function(generic, cls) {
    !is.null(utils::getS3method(generic, cls, optional = TRUE,
                                envir = asNamespace("generics")))
  }
  for (cls in .dmar_tidier_classes) {
    expect_true(is_registered("tidy", cls),
                info = paste("tidy method missing for", cls))
    expect_true(is_registered("glance", cls),
                info = paste("glance method missing for", cls))
  }
})


test_that("tidy() on a long CI matches the Fisher's Z interval from base R", {
  # Oracle: the Fisher's Z transformed interval for a correlation, in closed
  # form from base R, r = 0.5 with n = 50.
  #   z  = atanh(r), se = 1 / sqrt(n - 3)
  #   limits = tanh(z -/+ qnorm(0.975) * se)
  r <- 0.5; n <- 50
  se <- 1 / sqrt(n - 3)
  oracle_low  <- tanh(atanh(r) - stats::qnorm(0.975) * se)
  oracle_high <- tanh(atanh(r) + stats::qnorm(0.975) * se)

  res <- ci_r(r = r, n = n, conf_level = 0.95)
  expect_s3_class(res, "dmar_ci_long")

  td <- generics::tidy(res)
  expect_equal(nrow(td), 1L)
  expect_named(td, c("term", "estimate", "ci_lower", "ci_upper", "conf_level"))
  expect_equal(td$estimate, r, tolerance = 1e-12)
  expect_equal(td$ci_lower, oracle_low, tolerance = 1e-10)
  expect_equal(td$ci_upper, oracle_high, tolerance = 1e-10)
  expect_equal(td$conf_level, 0.95)

  # A single estimand has nothing extra at the model level, so glance()
  # coincides with tidy(). That equivalence is part of the documented contract.
  expect_equal(generics::glance(res), td)
})


test_that("tidy() on ci_pvaf() names the estimate row and carries a real number", {
  # Bargman (1970): a 5-group ANOVA with F = 11.221 on 4 and 50 degrees of
  # freedom, N = 55. Oracle: the sample proportion of variance accounted
  # for in closed form, df_1 * F / (df_1 * F + df_2).
  F_value <- 11.221; df_1 <- 4; df_2 <- 50; N <- 55
  oracle_estimate <- df_1 * F_value / (df_1 * F_value + df_2)

  res <- ci_pvaf(F_value = F_value, df_1 = df_1, df_2 = df_2, N = N)
  expect_s3_class(res, "dmar_ci_long")

  td <- generics::tidy(res)
  expect_equal(nrow(td), 1L)
  expect_named(td, c("term", "estimate", "ci_lower", "ci_upper", "conf_level"))
  expect_equal(td$term, "pvaf")
  expect_equal(td$estimate, oracle_estimate, tolerance = 1e-12)
  expect_equal(td$ci_lower, res$value[res$term == "lower_limit"])
  expect_equal(td$ci_upper, res$value[res$term == "upper_limit"])
  expect_equal(td$conf_level, 0.95)
  expect_equal(generics::glance(res), td)

  # The same contract on a dmar_ci_long sibling: the shared tidier names the
  # estimate row for the quantity and reports its value, never an NA row.
  sib <- generics::tidy(ci_r(r = 0.35, n = 100))
  expect_named(sib, c("term", "estimate", "ci_lower", "ci_upper", "conf_level"))
  expect_equal(sib$term, "r")
  expect_equal(sib$estimate, 0.35, tolerance = 1e-12)
  expect_false(anyNA(sib$estimate))
})


test_that("tidy() on an ANOVA effect size CI matches the noncentral F pivot", {
  F_value <- 5.2; df_effect <- 2; df_error <- 60; N <- 63

  # Oracle 1: the eta squared point estimate in closed form.
  oracle_estimate <- df_effect * F_value / (df_effect * F_value + df_error)

  # Oracle 2: the confidence limits from the noncentral F pivot, solved with
  # base R's pf() and uniroot(). The limits on the noncentrality lambda satisfy
  # P(F(df_effect, df_error, lambda) > F_value) = 1 - alpha/2 and alpha/2, and
  # eta squared = lambda / (lambda + N).
  lambda_at <- function(tail_prob) {
    stats::uniroot(function(L)
      stats::pf(F_value, df_effect, df_error, ncp = L,
                lower.tail = FALSE) - tail_prob,
      interval = c(0, 500), tol = 1e-10)$root
  }
  oracle_low  <- lambda_at(0.025) / (lambda_at(0.025) + N)
  oracle_high <- lambda_at(0.975) / (lambda_at(0.975) + N)

  res <- ci_eta_squared(F_value = F_value, df_effect = df_effect,
                        df_error = df_error, N = N, conf_level = 0.95)
  expect_s3_class(res, "dmar_ci_anova")

  td <- generics::tidy(res)
  expect_equal(nrow(td), 1L)
  expect_equal(td$term, "eta_squared")
  expect_equal(td$estimate, oracle_estimate, tolerance = 1e-10)
  expect_equal(td$ci_lower, oracle_low, tolerance = 1e-6)
  expect_equal(td$ci_upper, oracle_high, tolerance = 1e-6)

  expect_equal(generics::glance(res), td)
})


test_that("tidy() on a post hoc family matches stats::TukeyHSD", {
  # Oracle: base R's own Tukey HSD implementation, an independent computation
  # of the same simultaneous intervals and adjusted p-values.
  fit    <- stats::aov(weight ~ group, data = PlantGrowth)
  oracle <- stats::TukeyHSD(fit)$group

  res <- ci_tukey_kramer(fit, conf_level = 0.95)
  expect_s3_class(res, "dmar_post_hoc_ci")

  td <- generics::tidy(res)
  expect_equal(nrow(td), nrow(oracle))
  expect_named(td, c("term", "estimate", "ci_lower", "ci_upper",
                     "p_adjusted", "conf_level"))
  # DMAR spaces the contrast label ("trt1 - ctrl"); TukeyHSD does not.
  expect_equal(gsub(" ", "", td$term), rownames(oracle))
  expect_equal(td$estimate,   unname(oracle[, "diff"]),  tolerance = 1e-10)
  expect_equal(td$ci_lower,   unname(oracle[, "lwr"]),   tolerance = 1e-10)
  expect_equal(td$ci_upper,  unname(oracle[, "upr"]),   tolerance = 1e-10)
  expect_equal(td$p_adjusted, unname(oracle[, "p adj"]), tolerance = 1e-8)

  # glance() summarizes the family rather than repeating it.
  gl <- generics::glance(res)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$n_contrasts, nrow(oracle))
  expect_equal(gl$conf_level, 0.95)
})


test_that("tidy() on a sample size planner matches the noncentral t power", {
  # Oracle: the exact two-tailed power of the two-group t-test at the planned
  # per-group n, from base R's noncentral t distribution. With n per group,
  # nu = 2(n - 1) and the noncentrality is sqrt(n / 2) * smd.
  smd <- 0.5
  power_at <- function(n) {
    nu   <- 2 * (n - 1)
    ncp  <- sqrt(n / 2) * smd
    crit <- stats::qt(0.025, nu, lower.tail = FALSE)
    stats::pt(crit, nu, ncp, lower.tail = FALSE) +
      stats::pt(-crit, nu, ncp, lower.tail = TRUE)
  }
  # The smallest per-group n reaching 0.80 is the planner's answer.
  oracle_n <- Position(function(n) power_at(n) >= 0.80, 2:400) + 1L

  plan <- ss_power_smd(smd = smd, desired_power = 0.80, alpha_level = 0.05)
  expect_s3_class(plan, "dmar_ss_power")

  td <- generics::tidy(plan)
  expect_equal(nrow(td), 1L)
  expect_named(td, c("term", "estimate", "power"))
  expect_equal(td$term, "sample_size")
  expect_equal(td$estimate, as.numeric(oracle_n))
  expect_equal(td$power, power_at(oracle_n), tolerance = 1e-10)

  # glance() keeps the same head and appends the echoed planning inputs.
  gl <- generics::glance(plan)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl[, c("term", "estimate", "power")], td)
  expect_true(all(c("supposed_smd", "desired_power", "alpha_level") %in%
                    names(gl)))
  expect_equal(gl$supposed_smd, smd)
})


test_that("the tidiers return full precision, not the printed rounding", {
  # The dmar_tbl layer rounds for display only; tidy() must not.
  res <- ci_r(r = 0.5, n = 50, conf_level = 0.95)
  td  <- generics::tidy(res)
  expect_gt(abs(td$ci_lower - round(td$ci_lower, 3)), 0)
  expect_equal(td$ci_lower, res$value[res$term == "lower_limit"])
})
