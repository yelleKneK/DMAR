test_that("equivalence_r() returns documented rows", {
  res <- equivalence_r(r = 0.05, N = 200, cor_upper = 0.10)
  expect_setequal(res$term,
                  c("r", "z_lower_test", "z_upper_test",
                    "p_lower", "p_upper", "p_tost",
                    "ci_lower", "ci_upper",
                    "cor_lower", "cor_upper",
                    "equivalent", "N"))
  expect_identical(res$value[res$term == "N"], 200)
})

test_that("equivalence_r() declares equivalence for r near 0 with large n", {
  res <- equivalence_r(r = 0.01, N = 1000, cor_upper = 0.10)
  expect_equal(res$value[res$term == "equivalent"], 1)
})

test_that("equivalence_r() does not declare equivalence for r near the bound", {
  res <- equivalence_r(r = 0.18, N = 100, cor_upper = 0.20)
  expect_equal(res$value[res$term == "equivalent"], 0)
})

test_that("equivalence_r() raw-data interface matches summary interface", {
  set.seed(113)
  x <- rnorm(150); y <- 0.04 * x + rnorm(150)
  r_obs <- cor(x, y)
  r1 <- equivalence_r(x = x, y = y, cor_upper = 0.15)
  r2 <- equivalence_r(r = r_obs, N = 150, cor_upper = 0.15)
  expect_equal(r1$value[r1$term == "p_tost"],
               r2$value[r2$term == "p_tost"], tolerance = 1e-10)
})

test_that("results_sentence() reports the TOST p-value, on the table and on its tidy()", {
  # The 1.1.0 review: the sentence took the first p_terms row, p_lower, one
  # of the two one-sided tests, so a result that is not equivalent (TOST
  # p = .187) was written "p = < 0.0001". Through 1.0.1 the table's own
  # sentence did this; the first versions of 1.1.0 spread it to tidy().
  e <- equivalence_r(r = .15, N = 300, cor_lower = .2, cor_upper = .2)
  p <- setNames(e$value, e$term)
  expect_lt(p[["p_lower"]], 1e-4)
  expect_identical(p[["equivalent"]], 0)
  want <- paste0("r = 0.15, 90% CI [0.06, 0.24], TOST p = ",
                 format_p(p[["p_tost"]]))
  expect_identical(want, "r = 0.15, 90% CI [0.06, 0.24], TOST p = 0.1870")
  expect_identical(results_sentence(e), want)
  expect_identical(results_sentence(generics::tidy(e)), want)
  expect_identical(attr(e, "p_sentence"), "p_tost")
})
