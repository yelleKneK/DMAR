## compare_cov_structures() -- fits competing within-subject covariance structures
## and returns a model-fit comparison data.frame.

test_that("compare_cov_structures() returns the full eight-structure comparison frame", {
  skip_if_not_installed("nlme")
  set.seed(113)
  n <- 30; k <- 4
  dat <- data.frame(
    id   = rep(1:n, each = k),
    time = rep(1:k, n),
    y    = as.vector(t(matrix(rnorm(n * k), n, k) +
                         rep(rnorm(n, 0, 1), each = k)))
  )
  res <- compare_cov_structures(dat, outcome = "y", subject = "id", time = "time")
  expect_s3_class(res, "data.frame")
  expect_setequal(res$structure,
                  c("IND", "CS", "CSH", "AR1", "ARH1", "TOEP", "TOEPH", "UN"))
  expect_true(all(c("log_lik", "AIC", "BIC", "n_par") %in% names(res)))
  expect_true(all(is.finite(res$log_lik)))
  expect_true(all(res$n_par >= 4))
})

test_that("compare_cov_structures() fits each new structure without error", {
  skip_if_not_installed("nlme")
  skip_on_cran()  # refits UN beside each structure; the full frame above runs on CRAN
  set.seed(113)
  n <- 30; k <- 4
  dat <- data.frame(
    id   = rep(1:n, each = k),
    time = rep(1:k, n),
    y    = as.vector(t(matrix(rnorm(n * k), n, k) +
                         rep(rnorm(n, 0, 1), each = k)))
  )
  for (s in c("CSH", "ARH1", "TOEP", "TOEPH")) {
    res <- compare_cov_structures(dat, outcome = "y", subject = "id",
                                  time = "time", structures = c(s, "UN"))
    expect_true(is.finite(res$log_lik[res$structure == s]),
                info = paste("structure", s, "should fit"))
  }
})

test_that("compare_cov_structures() honors the structures argument", {
  skip_if_not_installed("nlme")
  set.seed(113)
  n <- 12
  dat <- data.frame(
    id   = rep(1:n, each = 3),
    time = rep(1:3, n),
    y    = rnorm(n * 3)
  )
  res <- compare_cov_structures(dat, outcome = "y", subject = "id", time = "time",
                                structures = c("IND", "UN"))
  expect_setequal(res$structure, c("IND", "UN"))
})

test_that("compare_cov_structures() accepts lowercase aliases", {
  skip_if_not_installed("nlme")
  set.seed(113)
  n <- 30; k <- 4
  dat <- data.frame(
    id   = rep(1:n, each = k),
    time = rep(1:k, n),
    y    = as.vector(t(matrix(rnorm(n * k), n, k) +
                         rep(rnorm(n, 0, 1), each = k)))
  )
  lower <- compare_cov_structures(dat, outcome = "y", subject = "id", time = "time",
                                  structures = c("cs", "csh", "ar1", "arh1",
                                                 "toep", "toeph", "un", "ind"))
  expect_setequal(lower$structure,
                  c("CS", "CSH", "AR1", "ARH1", "TOEP", "TOEPH", "UN", "IND"))
  # Lowercase and uppercase requests yield identical results.
  upper <- compare_cov_structures(dat, outcome = "y", subject = "id", time = "time",
                                  structures = c("CS", "CSH"))
  mixed <- compare_cov_structures(dat, outcome = "y", subject = "id", time = "time",
                                  structures = c("cs", "Csh"))
  expect_equal(mixed$log_lik, upper$log_lik)
})

test_that("compare_cov_structures() reports LRTs against UN for restricted structures", {
  skip_if_not_installed("nlme")
  set.seed(113)
  n <- 30; k <- 4
  dat <- data.frame(
    id   = rep(1:n, each = k),
    time = rep(1:k, n),
    y    = as.vector(t(matrix(rnorm(n * k), n, k) +
                         rep(rnorm(n, 0, 1), each = k)))
  )
  res <- compare_cov_structures(dat, outcome = "y", subject = "id", time = "time")
  # UN has the most parameters; every other structure is nested in it.
  expect_equal(res$n_par[res$structure == "UN"], max(res$n_par))
  restricted <- res[res$structure != "UN", ]
  expect_true(all(restricted$LRT_vs_UN_df > 0))
  expect_true(all(restricted$LRT_vs_UN_p >= 0 & restricted$LRT_vs_UN_p <= 1))
})
