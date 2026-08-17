test_that("ss_aipe_equivalence_r() returns documented rows", {
  res <- ss_aipe_equivalence_r(population_r = 0, width = 0.20)
  expect_setequal(res$term,
                  c("necessary_N", "width", "population_r",
                    "ci_width_expected"))
})

test_that("ss_aipe_equivalence_r() at rho = 0 matches the closed form", {
  # At rho = 0 the smallest N with width <= omega solves
  # z_{1-alpha} / sqrt(N - 3) <= atanh(omega / 2) exactly, so the search
  # must return ceiling(3 + (z / atanh(omega / 2))^2), computed here
  # independently of the function.
  for (w in c(0.10, 0.20, 0.30)) {
    oracle <- ceiling(3 + (qnorm(0.95) / atanh(w / 2))^2)
    res <- ss_aipe_equivalence_r(population_r = 0, width = w)
    expect_equal(res$value[res$term == "necessary_N"], oracle)
  }
})

test_that("ss_aipe_equivalence_r() is tight: N works and N - 1 does not", {
  res <- ss_aipe_equivalence_r(population_r = 0.30, width = 0.15)
  N <- res$value[res$term == "necessary_N"]
  width_at <- function(N, rho = 0.30, alpha = 0.05) {
    h <- qnorm(1 - alpha) / sqrt(N - 3)
    tanh(atanh(rho) + h) - tanh(atanh(rho) - h)
  }
  expect_lte(width_at(N), 0.15)
  expect_gt(width_at(N - 1), 0.15)
  expect_equal(res$value[res$term == "ci_width_expected"], width_at(N))
})

test_that("ss_aipe_equivalence_r() smaller width => larger N; larger |rho| => smaller N", {
  n_wide  <- ss_aipe_equivalence_r(population_r = 0, width = 0.30)$value[1]
  n_tight <- ss_aipe_equivalence_r(population_r = 0, width = 0.10)$value[1]
  expect_lt(n_wide, n_tight)
  # The interval narrows on the correlation scale away from zero, so the
  # widest-case default rho = 0 is conservative for any other value.
  n_zero <- ss_aipe_equivalence_r(population_r = 0,   width = 0.20)$value[1]
  n_half <- ss_aipe_equivalence_r(population_r = 0.5, width = 0.20)$value[1]
  expect_lt(n_half, n_zero)
})

test_that("ss_aipe_equivalence_r() assurance can only increase N and is seeded by the caller", {
  n_plain <- ss_aipe_equivalence_r(population_r = 0.30, width = 0.20)$value[1]
  set.seed(113)
  n_assur <- ss_aipe_equivalence_r(population_r = 0.30, width = 0.20,
                            assurance = 0.80)$value[1]
  expect_gte(n_assur, n_plain)
  set.seed(113)
  n_again <- ss_aipe_equivalence_r(population_r = 0.30, width = 0.20,
                            assurance = 0.80)$value[1]
  expect_identical(n_assur, n_again)
})

test_that("ss_aipe_equivalence_r() rejects bad inputs", {
  expect_error(ss_aipe_equivalence_r(width = 0), "in \\(0, 2\\)")
  expect_error(ss_aipe_equivalence_r(width = 2), "in \\(0, 2\\)")
  expect_error(ss_aipe_equivalence_r(population_r = 1, width = 0.2), "\\|r\\| < 1")
  expect_error(ss_aipe_equivalence_r(width = 0.2, alpha_level = 0.6),
               "in \\(0, 0.5\\)")
  expect_error(ss_aipe_equivalence_r(width = 0.2, assurance = 1.2),
               "in \\(0, 1\\)")
})
