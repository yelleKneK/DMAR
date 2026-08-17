test_that("ss_aipe_icc() returns documented rows and width <= target", {
  res <- ss_aipe_icc(rho = 0.7, k = 3, width = 0.20)
  expect_true(all(c("necessary_N", "expected_width", "rho", "k",
                    "width_target", "conf_level") %in% res$term))
  expect_lte(res$value[res$term == "expected_width"], 0.20 + 1e-6)
})

test_that("ss_aipe_icc() assurance inflates the sample size", {
  res_50 <- ss_aipe_icc(rho = 0.7, k = 3, width = 0.20)
  res_80 <- ss_aipe_icc(rho = 0.7, k = 3, width = 0.20, assurance = 0.80)
  expect_gt(res_80$value[1], res_50$value[1])
})

test_that("ss_aipe_icc() rejects bad inputs", {
  expect_error(ss_aipe_icc(rho = 1.5, k = 3, width = 0.2), "in \\[0, 1\\)")
  expect_error(ss_aipe_icc(rho = 0.7, k = 1, width = 0.2), ">= 2")
})

test_that("ss_aipe_icc() validates 'type' and records it on the result", {
  expect_error(ss_aipe_icc(rho = 0.7, k = 3, width = 0.2, type = "garbage"))
  res <- ss_aipe_icc(rho = 0.7, k = 3, width = 0.2, type = "ICC(1,k)")
  expect_identical(attr(res, "icc_type"), "ICC(1,k)")
})

test_that("ss_aipe_icc() average-of-k planning matches the closed-form inversion", {
  # Independent anchor. On the average-of-k scale the Bonett interval is
  # rho_k = 1 - exp(-2L) with L = -log(1 - rho_k) / 2 and
  # se_L = sqrt(k / (2 (k - 1) (n - 2))), so the full width is
  # 2 (1 - rho_k) sinh(2 z se_L) and the smallest n with width <= w
  # inverts in closed form. This is a different route from the function's
  # grid search over back-transformed widths.
  cf_n_avg <- function(rho_k, k, w, conf_level = 0.95) {
    z  <- qnorm(1 - (1 - conf_level) / 2)
    se <- asinh(w / (2 * (1 - rho_k))) / (2 * z)
    max(4, ceiling(2 + k / (2 * (k - 1) * se^2)))
  }
  cells <- list(c(0.70, 3, 0.20), c(0.70, 3, 0.10), c(0.90, 3, 0.10))
  for (cell in cells) {
    res <- ss_aipe_icc(rho = cell[1], k = cell[2], width = cell[3],
                       type = "ICC(1,k)")
    expect_identical(res$value[res$term == "necessary_N"],
                     as.numeric(cf_n_avg(cell[1], cell[2], cell[3])))
    expect_lte(res$value[res$term == "expected_width"], cell[3] + 1e-6)
  }
  # The three cells pin the audit's Monte Carlo confirmed values.
  n_at <- function(rho_k, w) {
    r <- ss_aipe_icc(rho = rho_k, k = 3, width = w, type = "ICC(1,k)")
    r$value[r$term == "necessary_N"]
  }
  expect_identical(n_at(0.70, 0.20), 110)
  expect_identical(n_at(0.70, 0.10), 421)
  expect_identical(n_at(0.90, 0.10), 52)
})

test_that("ss_aipe_icc() planning is shared within the single-rater and average-of-k families", {
  # The three single-rater forms share one planning variance, as do the
  # three average-of-k forms; the two families differ from each other.
  n_of <- function(type) {
    r <- ss_aipe_icc(rho = 0.7, k = 3, width = 0.20, type = type)
    r$value[r$term == "necessary_N"]
  }
  expect_identical(n_of("ICC(1,1)"), 69)   # unchanged single-rater plan
  expect_identical(n_of("ICC(2,1)"), 69)
  expect_identical(n_of("ICC(3,1)"), 69)
  expect_identical(n_of("ICC(1,k)"), 110)
  expect_identical(n_of("ICC(2,k)"), 110)
  expect_identical(n_of("ICC(3,k)"), 110)
})

test_that("ss_aipe_icc() average-of-k realized widths land on target (Monte Carlo)", {
  # skip_on_cran: the fast anchor that stays on CRAN is the closed-form
  # inversion test above ("average-of-k planning matches the closed-form
  # inversion"), which pins the same three cells deterministically.
  skip_on_cran()
  set.seed(2026)
  cells <- list(c(0.70, 3, 0.20), c(0.70, 3, 0.10), c(0.90, 3, 0.10))
  for (cell in cells) {
    res <- ss_aipe_icc_sensitivity(
      true_rho = cell[1], estimated_rho = cell[1], k = cell[2],
      width = cell[3], type = "ICC(1,k)", G = 2000
    )
    mean_width <- res$value[res$term == "mean_ci_width"]
    mean_icc   <- res$value[res$term == "mean_icc"]
    miss       <- res$value[res$term == "total_type_I_error"]
    # Realized F-based widths within noise (and Bonett-approximation
    # slack) of the target; the estimator tracks true rho on the
    # average-of-k scale; coverage near the nominal .95.
    expect_lt(abs(mean_width - cell[3]), 0.01)
    expect_lt(abs(mean_icc - cell[1]), 0.02)
    expect_gt(miss, 0.02)
    expect_lt(miss, 0.08)
  }
})
