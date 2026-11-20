# The directional test family's canonical vocabulary is snake_case
# ("two_sided"), with base R's "two.sided" accepted as an alias. Each
# function must return identical results under either spelling, and
# anything stored on the object must carry the canonical form.

test_that("two.sided and two_sided are the same alternative everywhere", {
  set.seed(113)
  x <- rnorm(12); y <- rnorm(12, 0.4)

  a <- welch_t(x, y, alternative = "two_sided")
  b <- welch_t(x, y, alternative = "two.sided")
  expect_equal(a, b)

  a <- summary_t_test(mean_1 = 3, sd_1 = 1, n_1 = 20,
                      mean_2 = 2.5, sd_2 = 1.2, n_2 = 18,
                      alternative = "two_sided")
  b <- summary_t_test(mean_1 = 3, sd_1 = 1, n_1 = 20,
                      mean_2 = 2.5, sd_2 = 1.2, n_2 = 18,
                      alternative = "two.sided")
  expect_equal(a, b)

  set.seed(113)
  a <- randomization_test(group_1 = x, group_2 = y, seed = 113,
                          alternative = "two_sided")
  set.seed(113)
  b <- randomization_test(group_1 = x, group_2 = y, seed = 113,
                          alternative = "two.sided")
  expect_equal(a, b)
  expect_identical(attr(a, "alternative"), "two_sided")
  expect_identical(attr(b, "alternative"), "two_sided")

  set.seed(113)
  a <- randomization_test_paired(x, y, seed = 113, alternative = "two_sided")
  set.seed(113)
  b <- randomization_test_paired(x, y, seed = 113, alternative = "two.sided")
  expect_equal(a, b)

  set.seed(113)
  dd <- data.frame(y = rnorm(45, rep(c(0, 0.5, 1), each = 15)),
                   g = factor(rep(c("control", "low", "high"), each = 15),
                              levels = c("control", "low", "high")))
  a <- ci_dunnett(dd$y, group = dd$g, alternative = "two_sided")
  b <- ci_dunnett(dd$y, group = dd$g, alternative = "two.sided")
  expect_equal(a, b)

  a <- power_fisher_exact(n_1 = 40, n_2 = 40, p_1 = 0.3, p_2 = 0.6,
                          alternative = "two_sided")
  b <- power_fisher_exact(n_1 = 40, n_2 = 40, p_1 = 0.3, p_2 = 0.6,
                          alternative = "two.sided")
  expect_equal(a, b)
})
