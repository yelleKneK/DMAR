test_that("obrien_test() returns a tidy one-row data.frame with expected columns", {
  res <- obrien_test(weight ~ group, data = PlantGrowth)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("statistic", "df_1", "df_2", "p_value",
                      "n_groups", "n_total", "method"))
})

test_that("obrien_test() formula and (vector, group) interfaces give the same result", {
  res_formula <- obrien_test(weight ~ group, data = PlantGrowth)
  res_vec     <- obrien_test(PlantGrowth$weight, PlantGrowth$group)
  expect_equal(res_formula$statistic, res_vec$statistic)
  expect_equal(res_formula$p_value,   res_vec$p_value)
})

test_that("obrien_test() degrees of freedom are k-1 and N-k", {
  res <- obrien_test(weight ~ group, data = PlantGrowth)
  expect_equal(res$df_1, length(unique(PlantGrowth$group)) - 1L)
  expect_equal(res$df_2, nrow(PlantGrowth) - length(unique(PlantGrowth$group)))
})

test_that("obrien_test() within-group mean of transformed scores equals sample variance", {
  # Recompute the transformation by hand on the Abdi example.
  D <- data.frame(
    group = factor(rep(1:3, each = 6)),
    score = c(5,6,7,8,9,10, 5,11,9,3,2,2, 0,1,2,4,9,16)
  )
  ybar <- tapply(D$score, D$group, mean)
  s2   <- tapply(D$score, D$group, var)
  n_g  <- as.integer(table(D$group))
  idx  <- as.integer(D$group)
  r_ij <- ((n_g[idx] - 1.5) * n_g[idx] * (D$score - ybar[idx])^2 -
           0.5 * s2[idx] * (n_g[idx] - 1)) /
          ((n_g[idx] - 1) * (n_g[idx] - 2))
  expect_equal(as.numeric(tapply(r_ij, D$group, mean)), as.numeric(s2))
})

test_that("obrien_test() flags grossly unequal variances with a small p-value", {
  set.seed(113)
  D <- data.frame(
    g = factor(rep(c("low", "high"), each = 50)),
    y = c(rnorm(50, sd = 1), rnorm(50, sd = 5))
  )
  res <- obrien_test(y ~ g, data = D)
  expect_lt(res$p_value, 0.001)
})

test_that("obrien_test() does not reject under equal variances at typical alpha", {
  set.seed(113)
  D <- data.frame(
    g = factor(rep(letters[1:3], each = 40)),
    y = rnorm(120)
  )
  res <- obrien_test(y ~ g, data = D)
  expect_gt(res$p_value, 0.05)
})

test_that("obrien_test() errors on too-small groups", {
  D <- data.frame(g = factor(c("a","a","b","b")), y = c(1, 2, 3, 4))
  expect_error(obrien_test(y ~ g, data = D), "at least 3 observations")
})

test_that("obrien_test() errors when only one group is present", {
  D <- data.frame(g = factor(rep("a", 10)), y = rnorm(10))
  expect_error(obrien_test(y ~ g, data = D), "two groups")
})

test_that("obrien_test() handles missing values via na.omit by default", {
  D <- data.frame(
    g = factor(rep(c("a","b","c"), each = 5)),
    y = c(1,2,3,4,NA, 6,7,8,9,10, 11,12,13,14,15)
  )
  res <- obrien_test(y ~ g, data = D)
  expect_equal(res$n_total, 14L)
})

test_that("the Abdi (2007) Hunter example reproduces F = 1.29 with p = .260", {
  # Anchors the ?obrien_test example prose. Abdi (2007, Table 6) prints
  # F = 1.29 on df = 1, 62; the p-value computed here is 0.2598, quoted as
  # .260 on the help page (Abdi's printed .2595 matches neither his own F
  # nor the exact statistic). Recomputed independently by feeding the
  # O'Brien transformed scores to a plain one-way ANOVA via lm()/anova().
  hunter_1964 <- data.frame(
    group = factor(
      c(rep("Control", 32), rep("Experimental", 32)),
      levels = c("Control", "Experimental")
    ),
    recall = c(
      5, 5, 5, 5, 5,
      6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
      7, 7, 7, 7, 7, 7, 7, 7, 7,
      8, 8, 8,
      9, 9,
      10, 10,
      6,
      7, 7,
      8, 8, 8, 8,
      9, 9, 9, 9, 9, 9, 9, 9, 9,
      10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
    )
  )
  res <- obrien_test(recall ~ group, data = hunter_1964)
  ybar <- tapply(hunter_1964$recall, hunter_1964$group, mean)
  s2   <- tapply(hunter_1964$recall, hunter_1964$group, var)
  n_g  <- as.integer(table(hunter_1964$group))
  idx  <- as.integer(hunter_1964$group)
  r_ij <- ((n_g[idx] - 1.5) * n_g[idx] * (hunter_1964$recall - ybar[idx])^2 -
           0.5 * s2[idx] * (n_g[idx] - 1)) /
          ((n_g[idx] - 1) * (n_g[idx] - 2))
  a <- anova(lm(r_ij ~ hunter_1964$group))
  expect_equal(res$statistic, a[1, "F value"], tolerance = 1e-10)
  expect_equal(res$p_value,   a[1, "Pr(>F)"],  tolerance = 1e-10)
  expect_equal(round(res$statistic, 2), 1.29)
  expect_equal(round(res$p_value, 3), 0.260)
})
