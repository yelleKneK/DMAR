test_that("smd_trimmed() returns documented rows", {
  set.seed(113)
  x <- rnorm(40, 0, 1); y <- rnorm(40, 0.5, 1)
  res <- smd_trimmed(x, y)
  expect_setequal(res$term,
                  c("smd_trimmed", "lower_limit", "upper_limit",
                    "trimmed_mean_x", "trimmed_mean_y",
                    "winsorized_sd_x", "winsorized_sd_y",
                    "winsorized_sd_pooled",
                    "h_x", "h_y", "trim", "df_yuen"))
})

test_that("smd_trimmed() reproduces the Algina-Keselman-Penfield d at a fixed seed", {
  set.seed(113)
  x <- rnorm(40, 0, 1); y <- rnorm(40, 0.5, 1)
  res <- smd_trimmed(x, y)

  # Hand computation of d_R = c * (m_t,x - m_t,y) / s_W,pooled at
  # trim = 0.20 on this fixed sample. The scaling constant c is
  # SD(X_W) for X ~ N(0, 1), obtained here by numerical integration
  # (an independent route from the closed form the function uses):
  # c = sqrt(0.41208673) = 0.6419398, the value Algina, Keselman,
  # and Penfield (2005) report rounded as 0.642.
  trim <- 0.20
  wsd <- function(z, tr) {
    z <- sort(z); m <- length(z); g <- floor(tr * m)
    if (g > 0) { z[seq_len(g)] <- z[g + 1L]; z[(m - g + 1L):m] <- z[m - g] }
    sd(z)
  }
  h_x <- length(x) - 2L * floor(trim * length(x))
  h_y <- length(y) - 2L * floor(trim * length(y))
  swp <- sqrt(((h_x - 1) * wsd(x, trim)^2 + (h_y - 1) * wsd(y, trim)^2) /
              (h_x + h_y - 2))
  z_c <- qnorm(1 - trim)
  win_var <- integrate(function(v) v^2 * dnorm(v), -z_c, z_c)$value +
    2 * trim * z_c^2
  d_hand <- sqrt(win_var) *
    (mean(x, trim = trim) - mean(y, trim = trim)) / swp

  expect_equal(res$value[res$term == "smd_trimmed"], d_hand,
               tolerance = 1e-10)
  expect_equal(res$value[res$term == "smd_trimmed"], -0.404725392337871,
               tolerance = 1e-9)
})

test_that("smd_trimmed()'s CI construction reproduces the Keselman et al. (2008) worked example", {
  # Keselman, Algina, Lix, Wilcox, and Deering (2008, Psychological
  # Methods, 13, 110-129) analyze a two-group example with 20%
  # trimming (their Tables 1 and 3): control n = 23, trimmed mean
  # 58.8, Winsorized SD 3.77; experimental n = 22, trimmed mean 69.8,
  # Winsorized SD 11.51. The paper reports the trimmed Welch-James
  # statistic t = 2.68 on 15.52 Yuen-Welch degrees of freedom (their
  # Equations 8 and 9) and noncentral t confidence intervals
  # [0.31, 3.37] and [0.10, 1.11] for the robust effect sizes
  # delta_R1 = 1.87 and delta_R2 = 0.61 (their Equation 19,
  # pp. 118-119). Feeding the published summary statistics through
  # the construction smd_trimmed() uses reproduces every printed
  # value to the precision the rounded inputs support.
  n_1 <- 23; n_2 <- 22
  h_1 <- n_1 - 2L * floor(0.20 * n_1)   # 15
  h_2 <- n_2 - 2L * floor(0.20 * n_2)   # 14
  m_t1 <- 58.8;  m_t2 <- 69.8
  s_W1 <- 3.77;  s_W2 <- 11.51

  v_1 <- s_W1^2 * (n_1 - 1) / (h_1 * (h_1 - 1))
  v_2 <- s_W2^2 * (n_2 - 1) / (h_2 * (h_2 - 1))
  t_yuen <- (m_t2 - m_t1) / sqrt(v_1 + v_2)
  df_y   <- (v_1 + v_2)^2 / (v_1^2 / (h_1 - 1) + v_2^2 / (h_2 - 1))
  expect_equal(t_yuen, 2.68, tolerance = 0.005)
  expect_equal(df_y, 15.52, tolerance = 0.0005)

  d_R1 <- 0.642 * (m_t2 - m_t1) / s_W1
  d_R2 <- 0.642 * (m_t2 - m_t1) / s_W2
  expect_equal(d_R1, 1.87, tolerance = 0.005)
  expect_equal(d_R2, 0.61, tolerance = 0.01)

  lam <- conf_limits_nct(ncp = t_yuen, df = df_y, conf_level = 0.95,
                         verbose = FALSE)$value
  expect_lt(max(abs(lam * d_R1 / t_yuen - c(0.31, 3.37))), 0.03)
  expect_lt(max(abs(lam * d_R2 / t_yuen - c(0.10, 1.11))), 0.015)

  # The pooled degrees of freedom the interval formerly used cannot
  # have produced the printed limits: at df = h_1 + h_2 - 2 = 27 the
  # lower limit for delta_R1 is 0.40, three times farther from the
  # printed 0.31 than the Yuen-Welch reproduction above.
  lam_pooled <- conf_limits_nct(ncp = t_yuen, df = h_1 + h_2 - 2,
                                conf_level = 0.95,
                                verbose = FALSE)$value
  expect_gt(abs(lam_pooled[1] * d_R1 / t_yuen - 0.31), 0.06)
})

test_that("smd_trimmed()'s CI inverts Yuen's t at the Yuen-Welch df", {
  # Independent recomputation of the interval on the help-page
  # example data, plus pins on the corrected values. Before the
  # interval moved to the Keselman et al. (2008) construction it
  # inverted the noncentral t at h_1 + h_2 - 2 degrees of freedom
  # with noncentrality d_R * sqrt(h_1 h_2 / (h_1 + h_2)), which gave
  # [-0.9743, 0.1692] on these data.
  set.seed(113)
  x <- rnorm(40, 0, 1); y <- rnorm(40, 0.5, 1)
  res <- smd_trimmed(x, y)

  trim <- 0.20
  wsd <- function(z, tr) {
    z <- sort(z); m <- length(z); g <- floor(tr * m)
    if (g > 0) { z[seq_len(g)] <- z[g + 1L]; z[(m - g + 1L):m] <- z[m - g] }
    sd(z)
  }
  n_1 <- length(x); n_2 <- length(y)
  h_1 <- n_1 - 2L * floor(trim * n_1)
  h_2 <- n_2 - 2L * floor(trim * n_2)
  v_1 <- wsd(x, trim)^2 * (n_1 - 1) / (h_1 * (h_1 - 1))
  v_2 <- wsd(y, trim)^2 * (n_2 - 1) / (h_2 * (h_2 - 1))
  t_yuen <- (mean(x, trim = trim) - mean(y, trim = trim)) /
    sqrt(v_1 + v_2)
  df_y <- (v_1 + v_2)^2 / (v_1^2 / (h_1 - 1) + v_2^2 / (h_2 - 1))
  expect_equal(res$value[res$term == "df_yuen"], df_y,
               tolerance = 1e-12)

  # d_R = scale_to_d * t_yuen, so the noncentrality limits transfer
  # to the d_R metric by the same factor.
  d_R <- res$value[res$term == "smd_trimmed"]
  lam <- conf_limits_nct(ncp = t_yuen, df = df_y, conf_level = 0.95,
                         verbose = FALSE)$value
  expect_equal(res$value[res$term == "lower_limit"],
               lam[1] * d_R / t_yuen, tolerance = 1e-9)
  expect_equal(res$value[res$term == "upper_limit"],
               lam[2] * d_R / t_yuen, tolerance = 1e-9)

  # Regression pins on the corrected interval.
  expect_equal(res$value[res$term == "lower_limit"],
               -0.882686927279091, tolerance = 1e-9)
  expect_equal(res$value[res$term == "upper_limit"],
               0.077499049060744, tolerance = 1e-9)
  expect_equal(res$value[res$term == "df_yuen"],
               45.996690063811, tolerance = 1e-9)
})

test_that("smd_trimmed() with trim = 0 reduces to standard d-like estimate", {
  set.seed(113)
  x <- rnorm(60, 0, 1); y <- rnorm(60, 0.5, 1)
  d_R <- smd_trimmed(x, y, trim = 0)$value[1]
  d_C <- smd(x, y)$value[smd(x, y)$term == "smd"]
  expect_equal(d_R, d_C, tolerance = 0.05)
})

test_that("smd_trimmed() is more stable than Cohen d under contamination", {
  set.seed(113)
  x <- rnorm(50, 0, 1)
  y <- c(rnorm(48, 0.5, 1), 50, -50)  # heavy contamination
  d_R <- abs(smd_trimmed(x, y)$value[1])
  d_C <- abs(smd(x, y)$value[smd(x, y)$term == "smd"])
  # d_R should be roughly comparable while d_C is suppressed by outliers,
  # or vice versa; check that they differ meaningfully.
  expect_gt(abs(d_R - d_C), 0.1)
})

test_that("smd_trimmed() rejects bad trim", {
  expect_error(smd_trimmed(1:10, 1:10, trim = 0.6),  "in \\[0, 0.5\\)")
  expect_error(smd_trimmed(1:10, 1:10, trim = -0.1), "in \\[0, 0.5\\)")
})

test_that(".akp_constant() is the square root of the Winsorized variance", {
  # Closed-form anchor that runs on CRAN: for X ~ N(0, 1) Winsorized
  # at gamma in each tail, Var(X_W) equals the integral of x^2 phi(x)
  # from -z to z plus 2 * gamma * z^2, with z = qnorm(1 - gamma), and
  # the Algina-Keselman-Penfield scaling constant is sqrt(Var(X_W)).
  for (tr in c(0.10, 0.30)) {
    z <- qnorm(1 - tr)
    win_var <- integrate(function(v) v^2 * dnorm(v), -z, z)$value +
      2 * tr * z^2
    expect_equal(DMAR:::.akp_constant(tr), sqrt(win_var),
                 tolerance = 1e-10)
  }
  # Independent integration gives Var(X_W) = 0.41208673 at trim = 0.20,
  # so the constant at the default is sqrt(0.41208673) = 0.6419398.
  expect_equal(DMAR:::.akp_constant(0.20), sqrt(0.41208673),
               tolerance = 1e-7)
})

test_that("smd_trimmed() is continuous in trim at the 0.20 default", {
  # 41 observations per group keeps floor(trim * n) = 8 on both sides
  # of 0.20, so a jump at the default could come only from the scaling
  # constant.
  set.seed(113)
  x <- rnorm(41, 0, 1); y <- rnorm(41, 0.5, 1)
  d_lo  <- smd_trimmed(x, y, trim = 0.1999)$value[1]
  d_mid <- smd_trimmed(x, y, trim = 0.20)$value[1]
  d_hi  <- smd_trimmed(x, y, trim = 0.2001)$value[1]
  expect_equal(d_lo, d_mid, tolerance = 1e-3)
  expect_equal(d_hi, d_mid, tolerance = 1e-3)
})

test_that("smd_trimmed() recovers Cohen's delta under normality at non-default trims", {
  # Monte Carlo confirmation of the Algina-Keselman-Penfield scaling
  # at trim = 0.10 and 0.30: under normal data with population
  # delta = 0.5, the average d_R must estimate 0.5, not 0.5 divided
  # by the Winsorized variance (0.74 and 2.50, the pre-fix values).
  # The fast anchor that stays on CRAN is the closed-form check above
  # that .akp_constant() equals the square root of the numerically
  # integrated Winsorized variance at these same trims.
  skip_on_cran()
  for (tr in c(0.10, 0.30)) {
    set.seed(113)
    d_bar <- mean(vapply(seq_len(200), function(i) {
      x <- rnorm(1000, 0.5, 1); y <- rnorm(1000, 0, 1)
      smd_trimmed(x, y, trim = tr)$value[1]
    }, numeric(1)))
    expect_equal(d_bar, 0.5, tolerance = 0.02)
  }
})
