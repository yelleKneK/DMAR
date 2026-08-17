## design_effect() -- Kish (1965) design effect, DEFT, and effective sample size.

test_that("design_effect() reproduces the textbook formula for balanced clusters", {
  # K = 30 clusters of size 20 each, ICC = 0.10. The textbook formula
  # is DEFF = 1 + (m - 1) * rho = 1 + 19 * 0.10 = 2.9.
  res <- design_effect(cluster_sizes = rep(20, 30), icc = 0.10)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_type(res$value, "double")

  expect_equal(res$value[res$term == "design_effect"], 2.9, tolerance = 1e-10)
  expect_equal(res$value[res$term == "deft"], sqrt(2.9), tolerance = 1e-10)
  expect_equal(res$value[res$term == "n_total"], 600)
  expect_equal(res$value[res$term == "effective_n"], 600 / 2.9, tolerance = 1e-10)
  expect_equal(res$value[res$term == "m_bar"], 20)
  expect_equal(res$value[res$term == "m_kish"], 20)
})

test_that("design_effect() returns DEFF = 1 when icc = 0 (no clustering effect)", {
  res <- design_effect(cluster_sizes = c(5, 10, 15, 20), icc = 0)
  expect_equal(res$value[res$term == "design_effect"], 1)
  expect_equal(res$value[res$term == "deft"], 1)
  expect_equal(res$value[res$term == "effective_n"], sum(c(5, 10, 15, 20)))
})

test_that("design_effect() handles empty clusters correctly (they drop out cleanly)", {
  res <- design_effect(cluster_sizes = c(0, 0, 20, 20, 20), icc = 0.10)
  # Empty clusters contribute 0 to N and to sum(m^2); they should leave
  # the design effect equal to the balanced-of-three-twenties case.
  expect_equal(res$value[res$term == "n_clusters_empty"], 2)
  expect_equal(res$value[res$term == "n_clusters_with_data"], 3)
  expect_equal(res$value[res$term == "n_clusters_total"], 5)
  expect_equal(res$value[res$term == "n_total"], 60)
  # m_kish = (3 * 20^2) / 60 = 20. DEFF = 1 + (20 - 1) * 0.10 = 2.9.
  expect_equal(res$value[res$term == "m_kish"], 20)
  expect_equal(res$value[res$term == "design_effect"], 2.9, tolerance = 1e-10)
})

test_that("design_effect() pulls DEFF toward 1 when many clusters are singletons", {
  # Many singletons: m_kish is small, so DEFF stays near 1 even at high ICC.
  many_singles <- c(rep(1, 50), 10, 10)
  res <- design_effect(cluster_sizes = many_singles, icc = 0.30)
  # m_kish = (50 * 1 + 2 * 100) / 70 = 250 / 70 = 3.5714...
  expect_equal(res$value[res$term == "m_kish"], 250 / 70, tolerance = 1e-10)
  # DEFF = 1 + (3.5714 - 1) * 0.30 = 1 + 0.7714 = 1.7714
  expect_equal(res$value[res$term == "design_effect"],
               1 + (250 / 70 - 1) * 0.30, tolerance = 1e-10)
  expect_equal(res$value[res$term == "n_clusters_singletons"], 50)
})

test_that("design_effect() effective_n equals N / DEFT^2", {
  res <- design_effect(cluster_sizes = c(15, 18, 22, 25, 30), icc = 0.20)
  N    <- res$value[res$term == "n_total"]
  deft_val <- res$value[res$term == "deft"]
  eff <- res$value[res$term == "effective_n"]
  expect_equal(eff, N / deft_val^2, tolerance = 1e-10)
})

test_that("design_effect() Kish's m_kish equals m_bar for equal cluster sizes; m_kish > m_bar otherwise", {
  # Equal sizes: m_kish = m_bar.
  res_eq <- design_effect(cluster_sizes = rep(15, 8), icc = 0.10)
  expect_equal(res_eq$value[res_eq$term == "m_kish"],
               res_eq$value[res_eq$term == "m_bar"])

  # Unequal sizes: m_kish > m_bar (size-weighted average pulls toward large
  # clusters because they carry more observations).
  res_uneq <- design_effect(cluster_sizes = c(2, 4, 6, 50), icc = 0.10)
  expect_gt(res_uneq$value[res_uneq$term == "m_kish"],
            res_uneq$value[res_uneq$term == "m_bar"])
})

test_that("design_effect() validates inputs", {
  expect_error(design_effect(cluster_sizes = "x", icc = 0.1),
               "must be a non-empty numeric vector")
  expect_error(design_effect(cluster_sizes = c(5, NA), icc = 0.1),
               "must not contain NA")
  expect_error(design_effect(cluster_sizes = c(5, -1), icc = 0.1),
               "non-negative integers")
  expect_error(design_effect(cluster_sizes = c(5, 2.5), icc = 0.1),
               "non-negative integers")
  expect_error(design_effect(cluster_sizes = c(5, 10), icc = -0.1),
               "in \\[0, 1\\)")
  expect_error(design_effect(cluster_sizes = c(5, 10), icc = 1.0),
               "in \\[0, 1\\)")
  expect_error(design_effect(cluster_sizes = c(0, 0, 0), icc = 0.1),
               "total sample size is 0")
})

test_that("design_effect() handles a single cluster (degenerate but valid)", {
  # 1 cluster of 20 observations; m_kish = 20, DEFF = 1 + 19 * rho.
  res <- design_effect(cluster_sizes = 20, icc = 0.10)
  expect_equal(res$value[res$term == "design_effect"], 2.9, tolerance = 1e-10)
  expect_equal(res$value[res$term == "n_clusters_total"], 1)
})

test_that("the deft row is the square root of the design_effect row", {
  res <- design_effect(cluster_sizes = c(5, 12, 20, 20, 30), icc = 0.15)
  expect_equal(res$value[res$term == "deft"],
               sqrt(res$value[res$term == "design_effect"]))
})
