test_that("procrustes_phi() = 1 for identical loadings", {
  l <- c(0.6, 0.7, 0.8, 0.5, 0.9)
  res <- procrustes_phi(l, l, n_resamples = 0)
  expect_equal(res$value[res$term == "tucker_phi"], 1, tolerance = 1e-12)
})

test_that("procrustes_phi() is scale-invariant", {
  l1 <- c(0.6, 0.7, 0.8)
  l2 <- 2 * l1
  expect_equal(procrustes_phi(l1, l2, n_resamples = 0)$value[1], 1,
               tolerance = 1e-12)
})

test_that("procrustes_phi() permutation test runs and returns a p-value in [0, 1]", {
  set.seed(113)
  l1 <- c(0.72, 0.65, 0.81, 0.55, 0.69, 0.62)
  l2 <- c(0.70, 0.62, 0.83, 0.58, 0.66, 0.65)
  res <- procrustes_phi(l1, l2, n_resamples = 2000)
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
  res <- procrustes_phi(a, a, n_resamples = 10000)
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
  res <- procrustes_phi(l1, l2, n_resamples = 200)
  expect_identical(res$term, c("tucker_phi", "p_value_perm", "n_resamples"))
  expect_s3_class(res, "dmar_tbl")
  no_perm <- procrustes_phi(l1, l2, n_resamples = 0)
  expect_identical(no_perm$term, "tucker_phi")
})

test_that("procrustes_phi() rejects mismatched lengths", {
  expect_error(procrustes_phi(c(0.5, 0.6), c(0.5, 0.6, 0.7)),
               "same length")
})

test_that("procrustes_phi() requires n_resamples to be a whole number of at least 0", {
  l1 <- c(0.72, 0.65, 0.81, 0.55, 0.69)
  l2 <- c(0.70, 0.62, 0.83, 0.58, 0.66)
  for (bad in list(2.5, -1, NA_real_, Inf, "10", c(10, 20))) {
    expect_error(procrustes_phi(l1, l2, n_resamples = bad), "'n_resamples'",
                 class = "dmar_count_invalid")
  }
  set.seed(113)
  a <- procrustes_phi(l1, l2, n_resamples = 200)
  set.seed(113)
  b <- procrustes_phi(l1, l2, n_resamples = 200L)
  expect_identical(a$value, b$value)
  expect_identical(a$value[a$term == "n_resamples"], 200)
  # 0 skips the test and leaves the point estimate alone.
  expect_identical(procrustes_phi(l1, l2, n_resamples = 0)$term, "tucker_phi")
})

test_that("procrustes_phi() draws and counts as DMAR 1.1.0 did when no loadings are tied", {
  # The p-value recorded from the reviewed 1.1.0 build (df120a0d).
  l1 <- c(0.72, 0.65, 0.81, 0.55, 0.69)
  l2 <- c(0.70, 0.62, 0.83, 0.58, 0.66)
  set.seed(113)
  res <- procrustes_phi(l1, l2, n_resamples = 2000)
  expect_equal(res$value[res$term == "p_value_perm"], 0.0079960019990004995,
               tolerance = 1e-14)
})

# Exact share of the permutations of loadings_2 whose |phi| is at least the
# observed one, for loadings given in hundredths: the numerator of phi is
# then a whole number and the denominator does not change, so ties are
# exact.
exact_phi_share <- function(l1, l2) {
  k1 <- round(100 * l1)
  k2 <- round(100 * l2)
  p <- length(k2)
  perms <- matrix(seq_len(p), nrow = 1L)
  for (j in seq_len(p - 1L)) {
    perms <- do.call(rbind, lapply(seq_len(nrow(perms)), function(r) {
      row <- perms[r, ]
      t(vapply(j:p, function(s) { row[c(j, s)] <- row[c(s, j)]; row },
               numeric(p)))
    }))
  }
  numerators <- apply(perms, 1L, function(i) sum(k1 * k2[i]))
  mean(abs(numerators) >= abs(sum(k1 * k2)))
}

test_that("procrustes_phi() counts a permutation that exchanges tied loadings as at least as extreme", {
  # loadings_1 has two items at 0.72. Exchanging their partners leaves phi
  # unchanged, so 2 of the 120 permutations are at least as extreme as the
  # observed arrangement.
  l1 <- c(0.6, 0.72, 0.65, 0.72, 0.81)
  l2 <- c(0.61, 0.76, 0.64, 0.68, 0.85)
  expect_equal(exact_phi_share(l1, l2), 2 / 120)
  set.seed(113)
  res <- procrustes_phi(l1, l2, n_resamples = 60000)
  expect_lt(abs(res$value[res$term == "p_value_perm"] - 2 / 120), 0.002)
})

test_that("procrustes_phi() agrees with whole-number enumeration when loadings are tied", {
  skip_on_cran()
  # Two-decimal loadings with ties in both vectors. The p-value is
  # (r + 1) / (m + 1) with r binomial, so it is compared with the exact
  # share within four Monte Carlo standard errors.
  set.seed(20261109)
  m <- 20000
  for (i in 1:8) {
    l1 <- sample(c(0.45, 0.55, 0.60, 0.72, 0.81), 5, replace = TRUE)
    l2 <- sample(c(0.40, 0.55, 0.61, 0.72, 0.85), 5, replace = TRUE)
    share <- exact_phi_share(l1, l2)
    res <- procrustes_phi(l1, l2, n_resamples = m)
    p <- res$value[res$term == "p_value_perm"]
    expect_lt(abs(p - (m * share + 1) / (m + 1)),
              4 * sqrt(share * (1 - share) / m) + 1e-12)
  }
})
