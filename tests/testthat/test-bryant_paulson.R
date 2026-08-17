# Tests for the Bryant-Paulson generalized studentized range and the
# simultaneous ANCOVA confidence intervals it produces. The published
# critical values of Bryant & Paulson (1976, Table 1) and Bryant & Bruvold
# (1980, Tables 1-2) provide exact targets.

test_that("qbryant_paulson reproduces the published critical value q_.05;1,6,14 = 4.83", {
  expect_equal(qbryant_paulson(0.95, num_covariates = 1, num_groups = 6, df = 14),
               4.83, tolerance = 5e-3)
})

test_that("p = 0 reduces exactly to the ordinary studentized range (Tukey)", {
  for (k in c(2, 4, 6)) for (nu in c(10, 20, 60)) {
    expect_equal(qbryant_paulson(0.95, num_covariates = 0, num_groups = k, df = nu),
                 qtukey(0.95, nmeans = k, df = nu), tolerance = 1e-6)
    expect_equal(pbryant_paulson(3.5, num_covariates = 0, num_groups = k, df = nu),
                 ptukey(3.5, nmeans = k, df = nu), tolerance = 1e-8)
  }
})

test_that("Bryant-Paulson critical values exceed Tukey's (random-covariate inflation)", {
  skip_on_cran()  # repeated numeric integration; fast anchors remain on CRAN
  for (p in 1:3) {
    bp <- qbryant_paulson(0.95, num_covariates = p, num_groups = 5, df = 25)
    tk <- qtukey(0.95, nmeans = 5, df = 25)
    expect_gt(bp, tk)
  }
  # More covariates -> larger critical value.
  q123 <- qbryant_paulson(0.95, num_covariates = 1:3, num_groups = 4, df = 20)
  expect_true(all(diff(q123) > 0))
})

test_that("k = 2 entries of Bryant & Bruvold (1980) Table 2 are reproduced", {
  skip_on_cran()  # repeated numeric integration; fast anchors remain on CRAN
  # Table 2, k = 2 column (where the Duncan value equals the BP q value).
  expect_equal(qbryant_paulson(0.95, 1, 2, 20),  3.03, tolerance = 5e-3) # p=1 a=.05
  expect_equal(qbryant_paulson(0.99, 1, 2, 20),  4.14, tolerance = 5e-3) # p=1 a=.01
  expect_equal(qbryant_paulson(0.95, 2, 2, 20),  3.10, tolerance = 5e-3) # p=2 a=.05
  expect_equal(qbryant_paulson(0.99, 3, 2, 120), 3.75, tolerance = 5e-3) # p=3 a=.01
})

test_that("p and q functions are mutual inverses", {
  skip_on_cran()  # repeated numeric integration; fast anchors remain on CRAN
  q <- qbryant_paulson(0.90, num_covariates = 2, num_groups = 4, df = 30)
  expect_equal(pbryant_paulson(q, 2, 4, 30), 0.90, tolerance = 1e-5)
})

test_that("lower_tail toggles the complementary probability", {
  expect_equal(pbryant_paulson(4, 1, 6, 14, lower_tail = FALSE),
               1 - pbryant_paulson(4, 1, 6, 14), tolerance = 1e-10)
})

test_that("invalid arguments are rejected", {
  expect_error(qbryant_paulson(0.95, num_covariates = -1, num_groups = 4, df = 10))
  expect_error(qbryant_paulson(0.95, num_covariates = 1,  num_groups = 1, df = 10))
  expect_error(qbryant_paulson(0.95, num_covariates = 1,  num_groups = 4, df = 0))
  expect_error(qbryant_paulson(1.5,  num_covariates = 1,  num_groups = 4, df = 10))
})

test_that("dbryant_paulson is a valid density: non-negative and integrates to one", {
  skip_on_cran()
  qs <- seq(0.5, 8, by = 0.5)
  dens <- dbryant_paulson(qs, num_covariates = 1, num_groups = 4, df = 20)
  expect_true(all(dens >= 0))
  area <- integrate(function(q) dbryant_paulson(q, 1, 4, 20),
                    lower = 0, upper = 30)$value
  expect_equal(area, 1, tolerance = 1e-4)
})

test_that("dbryant_paulson at p = 0 matches the derivative of the studentized range", {
  q0 <- 3.5
  h  <- 1e-4
  fd <- (ptukey(q0 + h, nmeans = 4, df = 20) -
         ptukey(q0 - h, nmeans = 4, df = 20)) / (2 * h)
  expect_equal(dbryant_paulson(q0, num_covariates = 0, num_groups = 4, df = 20),
               fd, tolerance = 1e-6)
})

test_that("ci_c_ancova_bp reproduces the paper's 0.278 critical difference", {
  adj <- c(3.595, 3.619, 4.102, 4.515, 4.618, 4.876)
  ci <- ci_c_ancova_bp(adj_means = adj, s_ancova = sqrt(0.01326), n = 4,
                       num_covariates = 1, df = 14)
  expect_equal(nrow(ci), choose(6, 2))
  half <- (ci$upper_limit - ci$lower_limit) / 2
  expect_true(all(abs(half - 0.278) < 1e-3))
  expect_equal(attr(ci, "critical_value"),
               qbryant_paulson(0.95, 1, 6, 14), tolerance = 1e-8)
})

test_that("ci_c_ancova_bp point estimates equal the contrasts of adjusted means", {
  skip_on_cran()  # a root find for an arithmetic check; the 0.278 anchor runs on CRAN
  adj <- c(3.595, 3.619, 4.102, 4.515, 4.618, 4.876)
  ci <- ci_c_ancova_bp(adj_means = adj, s_ancova = sqrt(0.01326), n = 4, df = 14,
                       c_weights = c(1, -1, 0, 0, 0, 0))
  expect_equal(ci$estimate, adj[1] - adj[2], tolerance = 1e-10)
})

test_that("pairwise and allowance widths agree for a pairwise contrast", {
  skip_on_cran()  # two root finds for one width; the 0.278 anchor above runs on CRAN
  adj <- c(1, 2, 3, 4)
  a <- ci_c_ancova_bp(adj_means = adj, s_ancova = 2, n = 10, df = 35,
                      c_weights = c(1, -1, 0, 0), contrast_type = "pairwise")
  b <- ci_c_ancova_bp(adj_means = adj, s_ancova = 2, n = 10, df = 35,
                      c_weights = c(1, -1, 0, 0), contrast_type = "allowance")
  expect_equal(a$lower_limit, b$lower_limit, tolerance = 1e-8)
})

test_that("default df for a one-way ANCOVA is N - k - p", {
  skip_on_cran()  # two root finds for a default check; the 0.278 anchor above runs on CRAN
  adj <- c(1, 2, 3)
  ci <- ci_c_ancova_bp(adj_means = adj, s_ancova = 1, n = 10, num_covariates = 1)
  # N = 30, k = 3, p = 1 -> df = 26; check the critical value matches.
  expect_equal(attr(ci, "critical_value"),
               qbryant_paulson(0.95, 1, 3, 26), tolerance = 1e-8)
})
