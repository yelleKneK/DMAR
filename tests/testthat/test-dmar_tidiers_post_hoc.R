# The @exportS3Method tags register these methods in NAMESPACE once
# devtools::document() has run. Under a bare load_all() in an isolated
# worktree the NAMESPACE may lag the source, so register them here so
# generics::tidy() / generics::glance() dispatch during the test. This
# is a no-op once the regenerated NAMESPACE carries the S3method lines.
if (exists("tidy.dmar_post_hoc_ci"))
  registerS3method("tidy", "dmar_post_hoc_ci", tidy.dmar_post_hoc_ci,
                   envir = asNamespace("generics"))
if (exists("glance.dmar_post_hoc_ci"))
  registerS3method("glance", "dmar_post_hoc_ci", glance.dmar_post_hoc_ci,
                   envir = asNamespace("generics"))

test_that("ci_dunnett returns a dmar_tbl / dmar_post_hoc_ci object", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  out <- ci_dunnett(fit, control = "ctrl")
  expect_s3_class(out, "dmar_tbl")
  expect_s3_class(out, "dmar_post_hoc_ci")
  expect_s3_class(out, "data.frame")
})

test_that("tidy.dmar_post_hoc_ci gives one broom row per Dunnett contrast", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  out <- ci_dunnett(fit, control = "ctrl")
  td  <- generics::tidy(out)

  expect_setequal(names(td),
    c("term", "estimate", "ci_lower", "ci_upper", "p_adjusted", "conf_level"))
  expect_equal(nrow(td), nrow(out))
  # broom columns match the source object's columns exactly.
  expect_equal(td$term,      out$contrast)
  expect_equal(td$estimate,  out$mean_difference)
  expect_equal(td$ci_lower,  out$lower_limit)
  expect_equal(td$ci_upper, out$upper_limit)
  expect_equal(td$p_adjusted, out$p_adjusted)
  expect_equal(unique(td$conf_level), 0.95)
})

test_that("glance.dmar_post_hoc_ci returns a one-row summary", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  out <- ci_dunnett(fit, control = "ctrl")
  gl  <- generics::glance(out)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$n_contrasts, nrow(out))
  expect_equal(gl$conf_level, 0.95)
})

test_that("ci_tukey_kramer tidies to broom columns (mean_difference estimate)", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  out <- ci_tukey_kramer(fit)
  expect_s3_class(out, "dmar_tbl")
  td  <- generics::tidy(out)
  expect_equal(nrow(td), nrow(out))
  expect_equal(td$estimate,  out$mean_difference)
  expect_equal(td$ci_lower,  out$lower_limit)
  expect_equal(td$ci_upper, out$upper_limit)
  expect_equal(nrow(generics::glance(out)), 1L)
})

test_that("ci_scheffe tidies to broom columns (contrast_value estimate)", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  out <- ci_scheffe(fit)
  expect_s3_class(out, "dmar_tbl")
  td  <- generics::tidy(out)
  expect_equal(nrow(td), nrow(out))
  # Scheffe's point estimate column is contrast_value, not mean_difference.
  expect_equal(td$estimate,  out$contrast_value)
  expect_equal(td$ci_lower,  out$lower_limit)
  expect_equal(td$ci_upper, out$upper_limit)
  expect_equal(nrow(generics::glance(out)), 1L)
})
