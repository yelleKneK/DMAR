test_that("equivalence_smd() returns documented rows", {
  res <- equivalence_smd(smd = 0.05, n_1 = 40, n_2 = 40, delta_upper = 0.4)
  expect_setequal(res$term,
                  c("smd", "t_lower", "t_upper", "df",
                    "p_lower", "p_upper", "p_tost",
                    "lower_limit", "upper_limit",
                    "delta_lower", "delta_upper", "equivalent"))
})

test_that("equivalence_smd() declares equivalence for a tiny d and large n", {
  res <- equivalence_smd(smd = 0.02, n_1 = 400, n_2 = 400, delta_upper = 0.3)
  expect_equal(res$value[res$term == "equivalent"], 1)
})

test_that("equivalence_smd() does not declare equivalence for a large d", {
  res <- equivalence_smd(smd = 0.7, n_1 = 50, n_2 = 50, delta_upper = 0.3)
  expect_equal(res$value[res$term == "equivalent"], 0)
})

test_that("equivalence_smd() raw-data interface agrees with summary interface", {
  set.seed(113)
  x <- rnorm(40, 100, 15); y <- rnorm(40, 101, 15)
  s_p <- sqrt(((40 - 1) * var(x) + (40 - 1) * var(y)) / (40 + 40 - 2))
  d_obs <- (mean(x) - mean(y)) / s_p
  r1 <- equivalence_smd(x = x, y = y, delta_upper = 0.3)
  r2 <- equivalence_smd(smd = d_obs, n_1 = 40, n_2 = 40, delta_upper = 0.3)
  expect_equal(r1$value[r1$term == "smd"],
               r2$value[r2$term == "smd"], tolerance = 1e-10)
  expect_equal(r1$value[r1$term == "p_tost"],
               r2$value[r2$term == "p_tost"], tolerance = 1e-10)
})

test_that("equivalence_smd() errors on bad bounds", {
  expect_error(equivalence_smd(smd = 0, n_1 = 30, n_2 = 30, delta_upper = -1),
               "positive number")
})
