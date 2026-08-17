test_that("mixed_anova() one-way returns expected columns", {
  set.seed(113)
  d <- data.frame(A = factor(rep(1:3, each = 10)),
                  y = rnorm(30) + 0.5 * rep(1:3, each = 10))
  res <- mixed_anova(d, "y", "A", A_type = "fixed")
  expect_setequal(colnames(res),
                  c("effect", "ss", "df", "ms", "denominator",
                    "F_value", "p_value"))
  expect_equal(res$effect, "A")
})

test_that("mixed_anova() two-way fixed-fixed denominator is MS_within", {
  set.seed(113)
  d <- expand.grid(A = factor(1:3), B = factor(1:2), rep = 1:5)
  d$y <- rnorm(nrow(d))
  res <- mixed_anova(d, "y", "A", "B", A_type = "fixed", B_type = "fixed")
  expect_true(all(res$denominator == "MS_within"))
})

test_that("mixed_anova() A fixed / B random uses MS_AB for A", {
  set.seed(113)
  d <- expand.grid(A = factor(1:3), B = factor(1:4), rep = 1:5)
  d$y <- rnorm(nrow(d))
  res <- mixed_anova(d, "y", "A", "B", A_type = "fixed", B_type = "random")
  expect_equal(res$denominator[res$effect == "A"],   "MS_AB")
  expect_equal(res$denominator[res$effect == "B"],   "MS_within")
  expect_equal(res$denominator[res$effect == "A:B"], "MS_within")
})

test_that("mixed_anova() both random uses MS_AB for both main effects", {
  set.seed(113)
  d <- expand.grid(A = factor(1:3), B = factor(1:4), rep = 1:4)
  d$y <- rnorm(nrow(d))
  res <- mixed_anova(d, "y", "A", "B", A_type = "random", B_type = "random")
  expect_equal(res$denominator[res$effect == "A"], "MS_AB")
  expect_equal(res$denominator[res$effect == "B"], "MS_AB")
})
