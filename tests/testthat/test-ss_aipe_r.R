# The exact full width of the back-transformed Fisher's Z interval at sample
# size n, the interval ss_aipe_r() plans for and correlations_test() reports
# for method = "pearson" (SE(z) = 1 / sqrt(n - 3); Bonett & Wright, 2000).
.fisher_z_width <- function(n, rho, conf_level) {
  z    <- atanh(rho)
  crit <- qnorm(1 - (1 - conf_level) / 2)
  se   <- 1 / sqrt(n - 3)
  tanh(z + crit * se) - tanh(z - crit * se)
}

test_that("ss_aipe_r() returns documented rows with the family classes", {
  res <- ss_aipe_r(rho = 0.30, width = 0.20)
  expect_s3_class(res, "dmar_ss_aipe")
  expect_s3_class(res, "dmar_tbl")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("necessary_N", "expected_width", "rho",
                    "width_target", "conf_level") %in% res$term))
  expect_true(is.numeric(res$value))
  n_val <- res$value[res$term == "necessary_N"]
  expect_gte(n_val, 4)
  expect_equal(n_val, round(n_val))
})

test_that("ss_aipe_r() expected_width row is the exact width at necessary_N", {
  res <- ss_aipe_r(rho = 0.30, width = 0.20)
  n   <- res$value[res$term == "necessary_N"]
  expect_equal(res$value[res$term == "expected_width"],
               .fisher_z_width(n, rho = 0.30, conf_level = 0.95))
  expect_lte(res$value[res$term == "expected_width"], 0.20)
})

test_that("ss_aipe_r() matches a brute-force search over a grid", {
  brute_force_n <- function(rho, width, conf_level) {
    n <- 4L
    while (.fisher_z_width(n, rho, conf_level) > width) n <- n + 1L
    n
  }
  for (rho in c(0, 0.30, -0.50, 0.70, 0.90)) {
    for (width in c(0.10, 0.20, 0.30)) {
      for (conf_level in c(0.90, 0.95, 0.99)) {
        planned <- ss_aipe_r(rho = rho, width = width,
                             conf_level = conf_level)
        n <- planned$value[planned$term == "necessary_N"]
        expect_equal(n, brute_force_n(rho, width, conf_level),
                     info = sprintf("rho = %.2f, width = %.2f, conf = %.2f",
                                    rho, width, conf_level))
        # Minimality restated directly: the width at n meets the target
        # and the width at n - 1 does not.
        expect_lte(.fisher_z_width(n, rho, conf_level), width)
        if (n > 4) expect_gt(.fisher_z_width(n - 1, rho, conf_level), width)
      }
    }
  }
})

test_that("ss_aipe_r() planning value closer to zero is more conservative", {
  n_0  <- ss_aipe_r(rho = 0.00, width = 0.20)$value[1]
  n_5  <- ss_aipe_r(rho = 0.50, width = 0.20)$value[1]
  n_9  <- ss_aipe_r(rho = 0.90, width = 0.20)$value[1]
  expect_gte(n_0, n_5)
  expect_gte(n_5, n_9)
})

test_that("ss_aipe_r() is symmetric in the sign of rho", {
  expect_equal(ss_aipe_r(rho = 0.40, width = 0.20)$value[1],
               ss_aipe_r(rho = -0.40, width = 0.20)$value[1])
})

test_that("ss_aipe_r() assurance inflates the sample size", {
  res_50 <- ss_aipe_r(rho = 0.30, width = 0.20)
  res_80 <- ss_aipe_r(rho = 0.30, width = 0.20, assurance = 0.80)
  res_99 <- ss_aipe_r(rho = 0.30, width = 0.20, assurance = 0.99)
  expect_gt(res_80$value[1], res_50$value[1])
  expect_gt(res_99$value[1], res_80$value[1])
})

test_that("ss_aipe_r() rejects bad inputs", {
  expect_error(ss_aipe_r(rho = 1.5, width = 0.2), "in \\(-1, 1\\)")
  expect_error(ss_aipe_r(rho = c(0.2, 0.3), width = 0.2), "in \\(-1, 1\\)")
  expect_error(ss_aipe_r(rho = 0.3, width = 0), "positive")
  expect_error(ss_aipe_r(rho = 0.3, width = -0.1), "positive")
  expect_error(ss_aipe_r(rho = 0.3, width = 0.2, conf_level = 1.5),
               "in \\(0, 1\\)")
  expect_error(ss_aipe_r(rho = 0.3, width = 0.2, assurance = 0.4),
               "in \\(0.5, 1\\)")
  expect_error(ss_aipe_r(rho = 0.3, width = 0.2, assurance = 1),
               "in \\(0.5, 1\\)")
})

test_that("ss_aipe_r() tidies through the dmar_ss_aipe registry", {
  res <- ss_aipe_r(rho = 0.30, width = 0.20)
  td  <- generics::tidy(res)
  expect_equal(td$estimate, res$value[res$term == "necessary_N"])
  expect_equal(td$width, 0.20)
  gl <- generics::glance(res)
  expect_equal(gl$estimate, td$estimate)
  expect_true("rho" %in% names(gl))
})
