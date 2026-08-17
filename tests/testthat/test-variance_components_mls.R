test_that("variance_components_mls() returns documented rows", {
  res <- variance_components_mls(ms_between = 6.0, ms_within = 1.5,
                                 df_between = 9, df_within = 40, n = 5)
  expect_setequal(res$term,
                   c("sigma2_between", "sigma2_b_lower", "sigma2_b_upper",
                     "sigma2_within",  "sigma2_w_lower", "sigma2_w_upper",
                     "icc", "icc_lower", "icc_upper"))
})

test_that("variance_components_mls() point estimates: sigma_b = (MS_b-MS_w)/n; sigma_w = MS_w", {
  res <- variance_components_mls(6.0, 1.5, 9, 40, 5)
  expect_equal(res$value[res$term == "sigma2_between"], (6 - 1.5) / 5,
               tolerance = 1e-12)
  expect_equal(res$value[res$term == "sigma2_within"], 1.5,
               tolerance = 1e-12)
  rho <- (6 - 1.5) / 5 / ((6 - 1.5) / 5 + 1.5)
  expect_equal(res$value[res$term == "icc"], rho, tolerance = 1e-12)
})

test_that("variance_components_mls() lower bound is non-negative", {
  res <- variance_components_mls(2.0, 1.5, 9, 40, 5)
  expect_gte(res$value[res$term == "sigma2_b_lower"], 0)
})

test_that("variance_components_mls() interval follows the documented Burdick-Graybill V_L and V_U", {
  # Recomputes the documented MLS form (squared constants and the
  # G_12 / H_12 cross terms; Burdick & Graybill, 1992, equations
  # 2.4.1-2.4.5) independently of the function. The symmetric form the
  # page previously printed, with unsquared constants and a
  # -MS_b MS_w / n cross term, gives [0.037, 1.763] on these inputs
  # instead of the function's [0.236, 3.692], so this recomputation
  # discriminates the two.
  msb <- 6; msw <- 1.5; dfb <- 9; dfw <- 40; n <- 5; alpha <- 0.05
  res <- variance_components_mls(msb, msw, dfb, dfw, n)
  v <- function(t) res$value[res$term == t]
  G1 <- 1 - 1 / stats::qf(1 - alpha / 2, dfb, Inf)
  G2 <- 1 - 1 / stats::qf(1 - alpha / 2, dfw, Inf)
  H1 <- 1 / stats::qf(alpha / 2, dfb, Inf) - 1
  H2 <- 1 / stats::qf(alpha / 2, dfw, Inf) - 1
  F_u <- stats::qf(1 - alpha / 2, dfb, dfw)
  F_l <- stats::qf(alpha / 2, dfb, dfw)
  G12 <- ((F_u - 1)^2 - (G1 * F_u)^2 - H2^2) / F_u
  H12 <- ((1 - F_l)^2 - (H1 * F_l)^2 - G2^2) / F_l
  V_L <- G1^2 * msb^2 + H2^2 * msw^2 + G12 * msb * msw
  V_U <- H1^2 * msb^2 + G2^2 * msw^2 + H12 * msb * msw
  expect_equal(v("sigma2_b_lower"), max(0, (msb - msw - sqrt(V_L)) / n),
               tolerance = 1e-12)
  expect_equal(v("sigma2_b_upper"), (msb - msw + sqrt(V_U)) / n,
               tolerance = 1e-12)
})
