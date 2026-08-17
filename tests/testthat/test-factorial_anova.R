# factorial_anova() computes Type I/II/III sums of squares by model comparison
# in base R. These tests pin the numbers to base anova() (Type I) and, where
# available, to car::Anova() (Types II and III), and check the balanced-design
# coincidence, the effect sizes, and the argument handling.

test_that("factorial_anova() reproduces the Type I sums of squares of base anova()", {
  d <- mtcars; d$cyl <- factor(d$cyl); d$am <- factor(d$am)
  # Weak effects legitimately clamp some effect-size CI lower limits to 0;
  # the clamp warning is asserted wherever these designs are fit below.
  expect_warning(fa <- factorial_anova(mpg ~ cyl * am, data = d, ss_type = 1),
                 "lower-limit clamp")
  a  <- stats::anova(stats::lm(mpg ~ cyl * am, data = d))
  expect_equal(fa$SS[1:3], a[["Sum Sq"]][1:3], tolerance = 1e-8)
  expect_equal(fa$df[1:3], a[["Df"]][1:3])
  expect_equal(fa$F_value[1:3], a[["F value"]][1:3], tolerance = 1e-8)
  # residual row
  expect_equal(fa$SS[4], a[["Sum Sq"]][4], tolerance = 1e-8)
  expect_equal(fa$df[4], a[["Df"]][4])
})

test_that("factorial_anova() matches car::Anova() for Types II and III", {
  skip_if_not_installed("car")
  d <- mtcars; d$cyl <- factor(d$cyl); d$am <- factor(d$am)
  fit <- stats::lm(mpg ~ cyl * am, data = d,
                   contrasts = list(cyl = contr.sum, am = contr.sum))

  # suppressWarnings on car::Anova: its internals partially match summary(corr =)
  # and warn under warnPartialMatchArgs; not a property under test.
  expect_warning(fa2 <- factorial_anova(mpg ~ cyl * am, data = d, ss_type = 2),
                 "lower-limit clamp")
  c2  <- suppressWarnings(car::Anova(fit, type = 2))
  keep <- rownames(c2) %in% c("cyl", "am", "cyl:am")
  expect_equal(fa2$SS[1:3], unname(c2[["Sum Sq"]][keep]), tolerance = 1e-6)
  expect_equal(fa2$F_value[1:3], unname(c2[["F value"]][keep]), tolerance = 1e-6)

  expect_warning(fa3 <- factorial_anova(mpg ~ cyl * am, data = d, ss_type = 3),
                 "lower-limit clamp")
  c3  <- suppressWarnings(car::Anova(fit, type = 3))
  keep3 <- rownames(c3) %in% c("cyl", "am", "cyl:am")
  expect_equal(fa3$SS[1:3], unname(c3[["Sum Sq"]][keep3]), tolerance = 1e-6)
  expect_equal(fa3$F_value[1:3], unname(c3[["F value"]][keep3]), tolerance = 1e-6)
})

test_that("factorial_anova() Type III matches the validated mtcars values", {
  d <- mtcars; d$cyl <- factor(d$cyl); d$am <- factor(d$am)
  expect_warning(fa <- factorial_anova(mpg ~ cyl * am, data = d, ss_type = 3),
                 "lower-limit clamp")
  expect_equal(fa$SS[1:3], c(410.4639, 29.86735, 25.43651), tolerance = 1e-3)
  expect_equal(fa$df, c(2, 1, 2, 26))
  expect_equal(fa$effect, c("cyl", "am", "cyl:am", "Residuals"))
})

test_that("the three types coincide for a balanced design", {
  d <- warpbreaks  # wool (2) x tension (3), 9 per cell: balanced
  expect_warning(s1 <- factorial_anova(breaks ~ wool * tension, data = d, ss_type = 1)$SS,
                 "lower-limit clamp")
  expect_warning(s2 <- factorial_anova(breaks ~ wool * tension, data = d, ss_type = 2)$SS,
                 "lower-limit clamp")
  expect_warning(s3 <- factorial_anova(breaks ~ wool * tension, data = d, ss_type = 3)$SS,
                 "lower-limit clamp")
  expect_equal(s1, s2, tolerance = 1e-8)
  expect_equal(s2, s3, tolerance = 1e-8)
})

test_that("factorial_anova() returns a dmar_tbl with effect sizes and CIs, and records the SS type", {
  d <- mtcars; d$cyl <- factor(d$cyl); d$am <- factor(d$am)
  expect_warning(fa <- factorial_anova(mpg ~ cyl * am, data = d, ss_type = 3),
                 "lower-limit clamp")
  expect_s3_class(fa, "dmar_tbl")
  expect_true(all(c("effect", "SS", "df", "F_value", "p_value",
                    "eta_squared_partial", "eta_squared_partial_lower",
                    "eta_squared_partial_upper", "omega_squared_partial",
                    "omega_squared_partial_lower", "omega_squared_partial_upper")
                  %in% names(fa)))
  expect_equal(attr(fa, "ss_type"), 3L)
  # every effect (not the residual row) has a partial eta^2 and its CI
  eff <- fa[fa$effect != "Residuals", ]
  expect_false(any(is.na(eff$eta_squared_partial)))
  expect_true(all(eff$eta_squared_partial_lower <= eff$eta_squared_partial &
                    eff$eta_squared_partial <= eff$eta_squared_partial_upper))
  # partial eta^2 equals SS_effect / (SS_effect + SS_error)
  sse <- fa$SS[fa$effect == "Residuals"]
  expect_equal(eff$eta_squared_partial, eff$SS / (eff$SS + sse), tolerance = 1e-8)
})

test_that("ss_type accepts integers and roman-numeral strings equivalently", {
  d <- mtcars; d$cyl <- factor(d$cyl); d$am <- factor(d$am)
  expect_warning(s3i <- factorial_anova(mpg ~ cyl * am, data = d, ss_type = 3)$SS,
                 "lower-limit clamp")
  expect_warning(s3r <- factorial_anova(mpg ~ cyl * am, data = d, ss_type = "III")$SS,
                 "lower-limit clamp")
  expect_equal(s3i, s3r)
  expect_warning(s2i <- factorial_anova(mpg ~ cyl * am, data = d, ss_type = 2)$SS,
                 "lower-limit clamp")
  expect_warning(s2r <- factorial_anova(mpg ~ cyl * am, data = d, ss_type = "II")$SS,
                 "lower-limit clamp")
  expect_equal(s2i, s2r)
})

test_that("three-way designs work and Type I sums equal the base anova partition", {
  set.seed(113)
  d <- expand.grid(A = factor(1:2), B = factor(1:2), C = factor(1:2))
  d <- d[rep(seq_len(nrow(d)), c(6, 4, 5, 7, 3, 8, 6, 5)), ]  # unbalanced
  d$y <- rnorm(nrow(d)) + as.numeric(d$A) + 0.5 * as.numeric(d$B)
  expect_warning(fa <- factorial_anova(y ~ A * B * C, data = d, ss_type = 1),
                 "lower-limit clamp")
  a  <- stats::anova(stats::lm(y ~ A * B * C, data = d))
  expect_equal(fa$SS[1:7], a[["Sum Sq"]][1:7], tolerance = 1e-8)
})

test_that("Types II and III are general: a 2x2x3 unbalanced design matches car::Anova", {
  # The marginality (Type II) and drop-one-with-contr.sum (Type III) rules must
  # hold for every effect of a three-way design, not just a two-way one. This
  # exercises the two-way interactions and the three-way interaction.
  skip_if_not_installed("car")
  set.seed(113)
  d <- expand.grid(A = factor(1:2), B = factor(1:2), C = factor(1:3))
  d <- d[rep(seq_len(nrow(d)), sample(3:9, 12, replace = TRUE)), ]
  d$y <- rnorm(nrow(d)) + as.numeric(d$A) + 0.5 * as.numeric(d$B) +
    0.3 * as.numeric(d$C)
  cs <- list(A = contr.sum, B = contr.sum, C = contr.sum)
  fit_t <- stats::lm(y ~ A * B * C, data = d)
  fit_s <- stats::lm(y ~ A * B * C, data = d, contrasts = cs)
  for (ty in 2:3) {
    expect_warning(fa <- factorial_anova(y ~ A * B * C, data = d, ss_type = ty),
                   "lower-limit clamp")
    fa <- fa[fa$effect != "Residuals", ]
    # suppressWarnings: car's internal partial argument matches, as above.
    ca <- suppressWarnings(car::Anova(if (ty == 3) fit_s else fit_t, type = ty))
    keep <- !rownames(ca) %in% c("Residuals", "(Intercept)")
    css <- stats::setNames(ca[["Sum Sq"]][keep], rownames(ca)[keep])
    cf  <- stats::setNames(ca[["F value"]][keep], rownames(ca)[keep])
    expect_equal(fa$SS, unname(css[fa$effect]), tolerance = 1e-6,
                 info = paste("SS, Type", ty))
    expect_equal(fa$F_value, unname(cf[fa$effect]), tolerance = 1e-6,
                 info = paste("F, Type", ty))
  }
})

test_that("inline factor() in the formula gives the same result as pre-made factors", {
  mt <- mtcars; mt$cyl <- factor(mt$cyl); mt$am <- factor(mt$am)
  for (ty in 1:3) {
    expect_warning(
      s_inline <- factorial_anova(mpg ~ factor(cyl) * factor(am), data = mtcars,
                                  ss_type = ty)$SS,
      "lower-limit clamp")
    expect_warning(
      s_pre <- factorial_anova(mpg ~ cyl * am, data = mt, ss_type = ty)$SS,
      "lower-limit clamp")
    expect_equal(s_inline, s_pre, tolerance = 1e-9, info = paste("Type", ty))
  }
})

test_that("factorial_anova() rejects bad input clearly", {
  d <- mtcars; d$cyl <- factor(d$cyl); d$am <- factor(d$am)
  expect_error(factorial_anova("mpg ~ cyl", data = d), "formula")
  expect_error(factorial_anova(mpg ~ cyl * am, data = as.list(d)), "data.frame")
  expect_error(factorial_anova(mpg ~ cyl * am, data = d, ss_type = 4), "ss_type")
  expect_error(factorial_anova(mpg ~ cyl * am, data = d, conf_level = 1.5), "conf_level")
  # empty cell -> rank deficient -> not estimable
  d2 <- subset(d, !(cyl == "8" & am == "1"))
  expect_error(factorial_anova(mpg ~ cyl * am, data = d2), "rank deficient|estimable")
})
