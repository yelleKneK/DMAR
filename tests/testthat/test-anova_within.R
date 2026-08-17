# Built-in within-subjects example used throughout these tests.
get_orthodont <- function() {
  skip_if_not_installed("nlme")
  nlme::Orthodont
}

test_that("anova_within() reproduces base R aov(... + Error()) on Orthodont", {
  ortho <- get_orthodont()
  res <- anova_within(ortho, id = "Subject", time = "age", outcome = "distance")

  ortho2 <- ortho; ortho2$age_f <- factor(ortho2$age)
  fit <- aov(distance ~ age_f + Error(Subject / age_f), data = ortho2)
  ref <- summary(fit)[["Error: Subject:age_f"]][[1]]

  none_row <- res[res$adjustment == "none", ]
  expect_equal(none_row$F_value, ref[["F value"]][1], tolerance = 1e-6)
  expect_equal(none_row$df_1,    ref[["Df"]][1])
  expect_equal(none_row$df_2,    ref[["Df"]][2])
  expect_equal(none_row$p_value, ref[["Pr(>F)"]][1], tolerance = 1e-12)
})

test_that("mauchly_test() returns 1 with k = 2 (sphericity is trivial)", {
  set.seed(113)
  Y <- matrix(rnorm(40), nrow = 20, ncol = 2)
  res <- mauchly_test(Y)
  expect_equal(res$W, 1)
  expect_equal(res$p_value, 1)
})

test_that("mauchly_test() chi square approximation is non-negative and df is correct", {
  set.seed(113)
  Y <- matrix(rnorm(30 * 5), nrow = 30)
  res <- mauchly_test(Y)
  expect_gte(res$statistic, 0)
  expect_equal(res$df, 5 * 4 / 2 - 1)  # = 9
  expect_true(res$W > 0 && res$W <= 1)
})

test_that("epsilon_corrections() returns three rows in a stable order", {
  set.seed(113)
  Y <- matrix(rnorm(20 * 4), nrow = 20)
  res <- epsilon_corrections(Y)
  expect_equal(res$epsilon_method,
               c("Greenhouse-Geisser", "Huynh-Feldt", "lower_bound"))
  # Bounds: lower_bound = 1/(k-1); GG and HF in [lower_bound, 1].
  expect_equal(res$epsilon[3], 1 / 3)
  expect_lte(res$epsilon[1], 1 + 1e-12)
  expect_lte(res$epsilon[2], 1 + 1e-12)
  expect_gte(res$epsilon[1], res$epsilon[3])
})

test_that("epsilon_corrections() returns 1 for all three when sphericity holds (perfectly spherical data)", {
  # Construct y_ij = subject_i + epsilon_ij with no condition effect: spherical.
  set.seed(113)
  n <- 100; k <- 4
  subj <- rnorm(n, sd = 2)
  Y <- matrix(NA, nrow = n, ncol = k)
  for (j in 1:k) Y[, j] <- subj + rnorm(n)
  res <- epsilon_corrections(Y)
  # Sample variability means GG/HF won't be exactly 1, but should be close.
  expect_gt(res$epsilon[1], 0.85)
  expect_gt(res$epsilon[2], 0.85)
})

test_that("anova_within() partial eta squared matches SS_effect / (SS_effect + SS_error)", {
  set.seed(113)
  n <- 25; k <- 4
  Y <- matrix(rnorm(n * k), nrow = n) + matrix(rep(c(0, .5, 1, 1.5), each = n),
                                                nrow = n)
  res  <- anova_within(Y)
  none <- res[res$adjustment == "none", ]
  ss_eff <- none$F_value * none$df_1 / none$df_2 / (1 + none$F_value * none$df_1 / none$df_2)
  expect_equal(attr(res, "partial_eta_squared"), unname(ss_eff), tolerance = 1e-12)
})

test_that("anova_within() epsilon-adjusted df scale linearly with the unadjusted df", {
  ortho <- get_orthodont()
  res <- anova_within(ortho, id = "Subject", time = "age", outcome = "distance")
  none <- res[res$adjustment == "none", ]
  for (adj in c("Greenhouse-Geisser", "Huynh-Feldt", "lower_bound")) {
    row <- res[res$adjustment == adj, ]
    expect_equal(row$df_1, none$df_1 * row$epsilon, tolerance = 1e-12)
    expect_equal(row$df_2, none$df_2 * row$epsilon, tolerance = 1e-12)
  }
})

test_that("anova_within() attaches the Mauchly test result as an attribute", {
  set.seed(113)
  Y <- matrix(rnorm(20 * 4), nrow = 20)
  res <- anova_within(Y)
  m <- attr(res, "mauchly")
  expect_s3_class(m, "data.frame")
  expect_equal(m$n_subjects, 20)
  expect_equal(m$n_levels, 4)
})

test_that("anova_within() rejects bad inputs", {
  expect_error(anova_within(matrix(1, nrow = 1, ncol = 4)),
               "At least 2 subjects")
  expect_error(anova_within(matrix(1:4, nrow = 4, ncol = 1)),
               "At least 2 within-subjects levels")
  bad <- matrix(c(1, NA), nrow = 2, ncol = 2)
  expect_error(anova_within(bad), "Missing")
})

test_that("long-format input gives the same result as wide-format", {
  ortho <- get_orthodont()
  long_res <- anova_within(ortho, id = "Subject",
                           time = "age", outcome = "distance")
  Y <- stats::reshape(ortho, timevar = "age", idvar = "Subject",
                      direction = "wide")
  Y <- as.matrix(Y[, grep("distance", names(Y))])
  storage.mode(Y) <- "double"
  wide_res <- anova_within(Y)
  expect_equal(long_res$F_value, wide_res$F_value, tolerance = 1e-12)
  expect_equal(long_res$p_value, wide_res$p_value, tolerance = 1e-12)
})
