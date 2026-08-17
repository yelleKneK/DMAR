test_that("ss_seq_c_sensitivity() exhibits first-order efficiency and
           near-nominal coverage", {
  skip_on_cran()
  res <- ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 2.5,
                              true_sigma = 15.67, G = 400, seed = 113)
  ratio <- res$value[res$term == "ratio_mean_N_n_star"]
  cover <- res$value[res$term == "coverage"]
  expect_gt(ratio, 0.9)
  expect_lt(ratio, 1.1)
  expect_gt(cover, 0.85)
  expect_lt(cover, 0.95)
})

test_that("ss_seq_c_sensitivity() oracle matches the known-sigma formula", {
  skip_on_cran()
  res <- ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 2.5,
                              true_sigma = 15.67, G = 10, seed = 113)
  n_star <- res$value[res$term == "n_star"]
  expected <- 2 * qnorm(0.95)^2 * 15.67^2 * 2 / 2.5^2
  expect_equal(n_star, expected, tolerance = 1e-10)
})

test_that("ss_seq_c_sensitivity() covers the true contrast under nonzero
           means", {
  skip_on_cran()
  res <- ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 3,
                              true_sigma = 10,
                              true_means = c(70, 65),
                              G = 300, seed = 113)
  cover <- res$value[res$term == "coverage"]
  expect_gt(cover, 0.83)
})

test_that("ss_seq_c_sensitivity() restores the caller's RNG state when
           seeded", {
  skip_on_cran()
  set.seed(20)
  before <- .Random.seed
  invisible(ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 3,
                                 true_sigma = 10, G = 10, seed = 113))
  expect_identical(.Random.seed, before)
})

test_that("ss_seq_c_sensitivity() leaves no RNG state when none existed", {
  skip_on_cran()
  # Simulate a fresh session that has not drawn a random number: the seeded
  # call must not leave a .Random.seed behind (the else branch of the restore).
  if (exists(".Random.seed", envir = globalenv()))
    rm(".Random.seed", envir = globalenv())
  invisible(ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 3,
                                 true_sigma = 10, G = 10, seed = 113))
  expect_false(exists(".Random.seed", envir = globalenv()))
})

test_that("ss_seq_c_sensitivity() is reproducible under a seed", {
  skip_on_cran()
  r1 <- ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 3,
                             true_sigma = 10, G = 50, seed = 113)
  r2 <- ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 3,
                             true_sigma = 10, G = 50, seed = 113)
  expect_equal(r1$value, r2$value)
})

test_that("ss_seq_c_sensitivity() validates its inputs", {
  expect_error(ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 0,
                                    true_sigma = 10),
               "half_width")
  expect_error(ss_seq_c_sensitivity(c_weights = c(1, -1), half_width = 2,
                                    true_sigma = 10,
                                    true_means = c(1, 2, 3)),
               "one mean per group")
})
