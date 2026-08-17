# manova_split_plot() is a thin tidy wrapper around car::Anova for a mixed-
# design MANOVA. We verify that (a) the function runs, (b) the output
# has the documented structure, and (c) the numeric values match car's
# printed table exactly. The car package is in Suggests, so the tests
# are gated on requireNamespace().

mm_fixture <- function() {
  set.seed(113)
  n_per <- 10
  data.frame(
    subject = factor(1:(2 * n_per)),
    group   = factor(rep(c("A", "B"), each = n_per)),
    t1      = c(rnorm(n_per, 0,   1), rnorm(n_per, 0,   1)),
    t2      = c(rnorm(n_per, 0.4, 1), rnorm(n_per, 0.8, 1)),
    t3      = c(rnorm(n_per, 0.8, 1), rnorm(n_per, 1.6, 1))
  )
}

test_that("manova_split_plot() returns the documented columns and effects", {
  skip_if_not_installed("car")
  d <- mm_fixture()
  res <- manova_split_plot(d, within = c("t1", "t2", "t3"), between = "group")

  expect_named(res, c("effect", "statistic_name", "statistic_value",
                      "F_approx", "df_1", "df_2", "p_value"))
  expect_setequal(res$effect,
                  c("A", "B", "interaction", "sum_of_squares_type"))

  tests <- res[res$effect != "sum_of_squares_type", ]
  expect_setequal(tests$statistic_name,
                  c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"))
  # 3 effects x 4 statistics = 12 test rows, plus 1 metadata row.
  expect_equal(nrow(tests), 12L)
  expect_equal(nrow(res), 13L)
})

test_that("manova_split_plot() reports the sum-of-squares type numerically", {
  skip_if_not_installed("car")
  d <- mm_fixture()

  # Default is Type III.
  res <- manova_split_plot(d, within = c("t1", "t2", "t3"), between = "group")
  type_row <- res[res$effect == "sum_of_squares_type", ]
  expect_equal(nrow(type_row), 1L)
  expect_equal(type_row$statistic_value, 3)
  # The value column stays numeric across the whole table.
  expect_type(res$statistic_value, "double")

  # ss_type = 2 changes the car call and is reported as 2 (smoke test
  # that the alternative type runs end to end).
  res2 <- manova_split_plot(d, within = c("t1", "t2", "t3"),
                       between = "group", ss_type = 2)
  expect_equal(
    res2$statistic_value[res2$effect == "sum_of_squares_type"], 2)
  expect_type(res2$statistic_value, "double")

  # The Roman-numeral string form is accepted and normalized.
  res3 <- manova_split_plot(d, within = c("t1", "t2", "t3"),
                       between = "group", ss_type = "III")
  expect_equal(
    res3$statistic_value[res3$effect == "sum_of_squares_type"], 3)

  # Type I is not available for this multivariate path: car computes
  # only Type II and Type III, so requesting it errors clearly.
  expect_error(
    manova_split_plot(d, within = c("t1", "t2", "t3"),
                 between = "group", ss_type = 1),
    "Type I")
  expect_error(
    manova_split_plot(d, within = c("t1", "t2", "t3"),
                 between = "group", ss_type = "I"),
    "Type I")

  # An out-of-range type is rejected before reaching car.
  expect_error(
    manova_split_plot(d, within = c("t1", "t2", "t3"),
                 between = "group", ss_type = 4),
    "ss_type")
})

test_that("manova_split_plot() values match car's printed table", {
  skip_if_not_installed("car")
  d <- mm_fixture()
  res <- manova_split_plot(d, within = c("t1", "t2", "t3"), between = "group")

  # Reference values are taken from print() of summary(car::Anova(...))
  # under R 4.5.2 with car 3.x; they only depend on the seed-113 data.
  # When the four statistics agree on s = 1 (the between effect with
  # 1 df), Pillai + Wilks = 1 and Hotelling-Lawley = Roy.
  between <- res[res$effect == "A", ]
  expect_equal(between$F_approx, rep(1.680826, 4L), tolerance = 1e-5)
  expect_equal(between$df_1,     rep(1, 4L))
  expect_equal(between$df_2,     rep(18, 4L))
  expect_equal(round(between$statistic_value[1] +     # Pillai
                     between$statistic_value[2], 6),  # Wilks
               1)
  expect_equal(between$statistic_value[3],            # HL
               between$statistic_value[4],            # Roy
               tolerance = 1e-12)

  within_eff <- res[res$effect == "B", ]
  expect_equal(within_eff$F_approx, rep(5.012657, 4L), tolerance = 1e-5)
  expect_equal(within_eff$df_1,     rep(2, 4L))
  expect_equal(within_eff$df_2,     rep(17, 4L))
})

test_that("manova_split_plot() rejects bad input", {
  skip_if_not_installed("car")
  d <- mm_fixture()
  expect_error(manova_split_plot(d, within = "t1", between = "group"),
               "2\\+ column names")
  expect_error(manova_split_plot(d, within = c("t1", "t2"), between = "not_there"),
               "single column name")
  # Between factor with 1 level.
  d_bad <- d
  d_bad$group <- factor("A", levels = "A")
  expect_error(manova_split_plot(d_bad, within = c("t1", "t2"), between = "group"),
               "at least 2 levels")
})
