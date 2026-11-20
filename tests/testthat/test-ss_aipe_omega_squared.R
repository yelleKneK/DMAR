test_that("ss_aipe_omega_squared() returns documented rows and width <= target", {
  skip_on_cran()
  # suppressMessages: the iterative search summarizes its noncentral F
  # lower-limit clamps in a message that is not what this test checks.
  res <- suppressMessages(
    ss_aipe_omega_squared(true_omega_squared = 0.10,
                          df_effect = 2, width = 0.10))
  expect_true(all(c("necessary_N", "expected_width",
                    "true_omega_squared", "df_effect",
                    "width", "conf_level") %in% res$term))
  expect_lte(res$value[res$term == "expected_width"], 0.10 + 1e-6)
})

test_that("ss_aipe_omega_squared() assurance inflates the sample size", {
  skip_on_cran()
  res_50 <- suppressMessages(
    ss_aipe_omega_squared(true_omega_squared = 0.10,
                          df_effect = 2, width = 0.10))
  res_80 <- suppressMessages(
    ss_aipe_omega_squared(true_omega_squared = 0.10,
                          df_effect = 2, width = 0.10,
                          assurance = 0.80))
  expect_gt(res_80$value[1], res_50$value[1])
})

test_that("ss_aipe_omega_squared() rejects bad inputs", {
  expect_error(ss_aipe_omega_squared(true_omega_squared = 1.2,
                                      df_effect = 2, width = 0.1),
               "in \\[0, 1\\)")
  expect_error(ss_aipe_omega_squared(true_omega_squared = 0.1,
                                      df_effect = 0, width = 0.1),
               ">=")
})
