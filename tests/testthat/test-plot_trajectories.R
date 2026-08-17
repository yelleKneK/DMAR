test_that("plot_trajectories() returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("nlme")
  d <- nlme::Orthodont
  p <- plot_trajectories(d, id = "Subject", time = "age",
                         outcome = "distance", group = "Sex")
  expect_s3_class(p, "ggplot")
})

test_that("plot_trajectories() facet=TRUE adds a FacetWrap layer", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("nlme")
  d <- nlme::Orthodont
  p <- plot_trajectories(d, id = "Subject", time = "age",
                         outcome = "distance", facet = TRUE)
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("plot_trajectories() rejects unknown columns", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(id = 1, t = 1, y = 1)
  expect_error(plot_trajectories(d, id = "id", time = "t", outcome = "missing"),
               "not found in 'data'")
  expect_error(plot_trajectories(d, id = "id", time = "missing", outcome = "y"),
               "not found in 'data'")
})

test_that("plot_trajectories() ids subset is enforced", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(
    id = rep(1:3, each = 4),
    t  = rep(1:4, 3),
    y  = rnorm(12)
  )
  p <- plot_trajectories(d, id = "id", time = "t", outcome = "y", ids = c(1, 3))
  used <- unique(p$data$id)
  expect_setequal(used, c(1, 3))
  expect_error(plot_trajectories(d, id = "id", time = "t", outcome = "y", ids = c(1, 99)),
               "not present in the data")
})

test_that("plot_trajectories() complains when more than one subset arg is supplied", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(id = rep(1:5, each = 2), t = rep(1:2, 5), y = 1:10)
  expect_error(
    plot_trajectories(d, id = "id", time = "t", outcome = "y",
                      ids = 1:2, n_random = 3),
    "at most one"
  )
})

test_that("plot_trajectories() pct_random accepts proportions and percentages", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(id = rep(1:10, each = 2), t = rep(1:2, 10), y = 1:20)
  set.seed(113)
  p1 <- plot_trajectories(d, id = "id", time = "t", outcome = "y",
                          pct_random = 0.5)
  set.seed(113)
  p2 <- plot_trajectories(d, id = "id", time = "t", outcome = "y",
                          pct_random = 50)
  expect_equal(nrow(p1$data), nrow(p2$data))
})

test_that("plot_trajectories_fitted() returns a ggplot with a quality_of_fit attribute (lme)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("nlme")
  fm <- nlme::lme(distance ~ age, random = ~ age | Subject,
                  data = nlme::Orthodont)
  p <- plot_trajectories_fitted(fm)
  expect_s3_class(p, "ggplot")
  qof <- attr(p, "quality_of_fit")
  expect_s3_class(qof, "data.frame")
  expect_named(qof, c("Subject", "r_squared", "rmse"))
  expect_equal(nrow(qof), length(unique(nlme::Orthodont$Subject)))
})

test_that("plot_trajectories_fitted() works on an lmerMod fit (sleepstudy)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("lme4")
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                   data = lme4::sleepstudy)
  p <- plot_trajectories_fitted(fm)
  expect_s3_class(p, "ggplot")
  qof <- attr(p, "quality_of_fit")
  expect_equal(nrow(qof), length(unique(lme4::sleepstudy$Subject)))
  # All R^2 values should be in [0, 1].
  expect_true(all(qof$r_squared >= 0 & qof$r_squared <= 1))
  # RMSE is non-negative.
  expect_true(all(qof$rmse >= 0))
})

test_that("plot_trajectories_fitted() per-subject R^2 matches a hand calculation", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("lme4")
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                   data = lme4::sleepstudy)
  p <- plot_trajectories_fitted(fm)
  qof <- attr(p, "quality_of_fit")

  # Hand-compute for one subject from fitted() and observed.
  d <- lme4::sleepstudy
  yhat <- fitted(fm)
  s1 <- d$Subject == levels(d$Subject)[1]
  expected_R2 <- cor(d$Reaction[s1], yhat[s1])^2
  expected_rmse <- sqrt(mean((d$Reaction[s1] - yhat[s1])^2))
  row <- qof[as.character(qof$Subject) == levels(d$Subject)[1], ]
  expect_equal(row$r_squared, expected_R2,   tolerance = 1e-8)
  expect_equal(row$rmse,      expected_rmse, tolerance = 1e-8)
})

test_that("plot_trajectories_fitted() rejects unsupported model classes", {
  skip_if_not_installed("ggplot2")
  fit <- lm(mpg ~ wt, data = mtcars)
  expect_error(plot_trajectories_fitted(fit), "lme.*nlme|lmer")
})

test_that("plot_trajectories_fitted() respects the n_random subset", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("lme4")
  fm <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                   data = lme4::sleepstudy)
  p <- plot_trajectories_fitted(fm, n_random = 5)
  qof <- attr(p, "quality_of_fit")
  # quality_of_fit is computed for ALL subjects; only the plot is subset.
  # The plot data should have 5 unique subjects.
  obs_layer <- p$layers[[1]]$data  # geom_point layer
  expect_equal(length(unique(obs_layer$Subject)), 5L)
})
