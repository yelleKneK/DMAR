test_that("procrustes_phi() = 1 for identical loadings", {
  l <- c(0.6, 0.7, 0.8, 0.5, 0.9)
  res <- procrustes_phi(l, l, n_perm = 0)
  expect_equal(res$value[res$term == "tucker_phi"], 1, tolerance = 1e-12)
})

test_that("procrustes_phi() is scale-invariant", {
  l1 <- c(0.6, 0.7, 0.8)
  l2 <- 2 * l1
  expect_equal(procrustes_phi(l1, l2, n_perm = 0)$value[1], 1,
               tolerance = 1e-12)
})

test_that("procrustes_phi() permutation test runs and returns a p-value in [0, 1]", {
  set.seed(113)
  l1 <- c(0.72, 0.65, 0.81, 0.55, 0.69, 0.62)
  l2 <- c(0.70, 0.62, 0.83, 0.58, 0.66, 0.65)
  res <- procrustes_phi(l1, l2, n_perm = 2000)
  expect_true("p_value_perm" %in% res$term)
  p <- res$value[res$term == "p_value_perm"]
  expect_gte(p, 0); expect_lte(p, 1)
})

test_that("procrustes_phi() permutation p-value carries the add-one correction", {
  # Identical vectors with eight distinct loadings: the observed phi is 1,
  # so a sampled permutation counts only when it reproduces phi exactly.
  # Before the (r + 1) / (m + 1) correction this call reported a p-value
  # of exactly 0, a value a sampled permutation test cannot support
  # (Phipson & Smyth, 2010).
  set.seed(113)
  a <- c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)
  res <- procrustes_phi(a, a, n_perm = 10000)
  p <- res$value[res$term == "p_value_perm"]
  expect_gt(p, 0)
  expect_gte(p, 1 / 10001)
  expect_lte(p, 1)
  # The (r + 1) / (m + 1) form means (m + 1) * p recovers the whole
  # number r + 1; the uncorrected proportion r / m fails this whenever
  # 0 < r < m.
  expect_equal(10001 * p, round(10001 * p), tolerance = 1e-8)
})

test_that("procrustes_phi() returns the documented rows and class", {
  l1 <- c(0.72, 0.65, 0.81, 0.55, 0.69)
  l2 <- c(0.70, 0.62, 0.83, 0.58, 0.66)
  res <- procrustes_phi(l1, l2, n_perm = 200)
  expect_identical(res$term, c("tucker_phi", "p_value_perm", "n_perm"))
  expect_s3_class(res, "dmar_tbl")
  no_perm <- procrustes_phi(l1, l2, n_perm = 0)
  expect_identical(no_perm$term, "tucker_phi")
})

test_that("procrustes_phi() rejects mismatched lengths", {
  expect_error(procrustes_phi(c(0.5, 0.6), c(0.5, 0.6, 0.7)),
               "same length")
})
