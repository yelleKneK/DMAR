test_that("ci_scheffe() default returns all pairwise contrasts", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  res <- ci_scheffe(fit)
  expect_equal(nrow(res), 3)
  expect_setequal(colnames(res),
                  c("contrast", "contrast_value", "se", "F_statistic",
                    "lower_limit", "upper_limit", "p_adjusted"))
})

test_that("ci_scheffe() supports custom contrasts", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  cmat <- matrix(c(-1, 0.5, 0.5), 3, 1,
                 dimnames = list(c("ctrl", "trt1", "trt2"),
                                 "avg(trt) - ctrl"))
  res <- ci_scheffe(fit, contrasts = cmat)
  expect_equal(nrow(res), 1)
  expect_equal(res$contrast, "avg(trt) - ctrl")
})

test_that("ci_scheffe() is wider than Tukey-Kramer on pairwise contrasts", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  s <- ci_scheffe(fit)
  t <- ci_tukey_kramer(fit)
  s_w <- s$upper_limit - s$lower_limit
  t_w <- t$upper_limit - t$lower_limit
  expect_true(all(s_w >= t_w - 1e-10))
})

test_that("ci_scheffe() rejects contrasts that do not sum to zero", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  expect_error(ci_scheffe(fit, contrasts = matrix(c(1, 1, 0), 3, 1)),
               "sum to zero")
})

test_that("ci_scheffe() default returns all a(a-1)/2 pairwise contrasts", {
  # Four groups give choose(4, 2) = 6 pairwise contrasts; the help page
  # previously said a - 1 (which would be 3 here).
  set.seed(113)
  d <- data.frame(y = rnorm(40), g = factor(rep(letters[1:4], each = 10)))
  res <- ci_scheffe(lm(y ~ g, data = d))
  expect_equal(nrow(res), choose(4, 2))
})
