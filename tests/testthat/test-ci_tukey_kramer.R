test_that("ci_tukey_kramer() returns documented columns", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  res <- ci_tukey_kramer(fit)
  expect_setequal(colnames(res),
                  c("contrast", "mean_difference", "se", "q_statistic",
                    "lower_limit", "upper_limit", "p_adjusted"))
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
