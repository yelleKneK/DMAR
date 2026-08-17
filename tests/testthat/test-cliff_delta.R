test_that("cliff_delta() returns point estimate, CI, variance, and tie/win/lose counts", {
  set.seed(113)
  res <- cliff_delta(rnorm(30, 0, 1), rnorm(30, 0.5, 1))
  expect_setequal(res$term,
                   c("cliff_delta", "lower_limit", "upper_limit",
                     "var_cliff_delta", "p_y1_greater", "p_y1_less",
                     "p_tied"))
  expect_lte(res$value[res$term == "cliff_delta"], 1)
  expect_gte(res$value[res$term == "cliff_delta"], -1)
})

test_that("cliff_delta() = 2A - 1 (Vargha-Delaney equivalence)", {
  set.seed(113)
  g1 <- rnorm(40); g2 <- rnorm(40, 0.4)
  d  <- cliff_delta(g1, g2)$value[1]
  # A = p_y1_greater + 0.5 * p_tied
  cd_obj <- cliff_delta(g1, g2)
  A <- cd_obj$value[cd_obj$term == "p_y1_greater"] + 0.5 *
       cd_obj$value[cd_obj$term == "p_tied"]
  expect_equal(d, 2 * A - 1, tolerance = 1e-10)
})

test_that("cliff_delta() CI bounds stay within [-1, 1]", {
  set.seed(113)
  res <- cliff_delta(rnorm(20), rnorm(20, 1))
  expect_gte(res$value[res$term == "lower_limit"], -1)
  expect_lte(res$value[res$term == "upper_limit"], 1)
})

test_that("cliff_delta() handles ties", {
  o1 <- c(1, 2, 2, 3, 3, 3, 4, 4, 5)
  o2 <- c(2, 3, 3, 4, 4, 5, 5, 5)
  res <- cliff_delta(o1, o2)
  expect_gt(res$value[res$term == "p_tied"], 0)
})
