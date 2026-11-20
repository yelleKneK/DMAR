test_that("ci_tukey_kramer() returns documented columns", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  res <- ci_tukey_kramer(fit)
  expect_setequal(colnames(res),
                  c("contrast", "mean_difference", "se", "q_statistic",
                    "ci_lower", "ci_upper", "p_adjusted"))
  expect_equal(nrow(res), 3)  # choose(3, 2)
})

test_that("ci_tukey_kramer() agrees with stats::TukeyHSD() on a balanced design", {
  fit  <- aov(weight ~ group, data = PlantGrowth)
  res  <- ci_tukey_kramer(fit)
  ref  <- as.data.frame(stats::TukeyHSD(fit)$group)
  # Match row ordering by contrast label structure.
  expect_equal(sort(res$mean_difference), sort(ref$diff), tolerance = 1e-6)
})

test_that("ci_tukey_kramer() vector / group interface matches lm interface", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  r1  <- ci_tukey_kramer(fit)
  r2  <- ci_tukey_kramer(PlantGrowth$weight, group = PlantGrowth$group)
  expect_equal(r1$mean_difference, r2$mean_difference, tolerance = 1e-10)
})

test_that("ci_tukey_kramer() uses the exact studentized range at 2 error df", {
  skip_on_cran()  # a spline of the range distribution for twenty means
  # Twenty groups, two of them with a second observation, so 2 error df.
  # Oracles, computed without package code (the 1.1.0 review's F3
  # evidence, 11_test_oracles.R): MS_E = 0.04031269637356, the "2 - 1"
  # standard error 0.1419730544391, its q = 11.29198760883, the upper tail of
  # the studentized range there, 0.1064739333761, and
  # q_.99;20,2 = 37.943462287637 (qtukey()
  # gives 44.526266).
  set.seed(113)
  g <- factor(c(rep(1:20, each = 1), 1, 2))
  y <- rnorm(length(g)) + as.integer(g) / 4
  out <- ci_tukey_kramer(y, g, conf_level = 0.99)
  r <- out[out$contrast == "2 - 1", ]
  expect_equal(r$se, 0.1419730544391, tolerance = 1e-10)
  expect_equal(r$q_statistic, 11.29198760883, tolerance = 1e-10)
  expect_equal(r$p_adjusted, 0.1064739333761, tolerance = 1e-7)
  expect_equal((r$ci_upper - r$ci_lower) / 2, 37.943462287637 * r$se,
               tolerance = 1e-8)
})
