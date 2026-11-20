# Data and independent oracles shared by the tests of eta_squared(),
# omega_squared(), ci_eta_squared(), and ci_omega_squared().

# Maxwell, Delaney, and Kelley (2027), Table 7.5 (p. 323 of the 7 July
# 2026 draft): blood pressure in a 2 (biofeedback, A) by 3 (drug, B)
# design with 5 per cell, N = 30, typed from the table. Table 7.7 (p. 324)
# gives the analysis: SS 1,080 (A, 1 df), 3,420 (B, 2 df), 780 (AB, 2 df),
# 3,738 (within cells, 24 df), and Table 7.8 (p. 329) the effect sizes to
# two decimals. An exact reproduction of the raw data.
.mdk_blood_pressure <- function() {
  data.frame(
    A = factor(rep(c("present", "absent"), each = 15),
               levels = c("present", "absent")),
    B = factor(rep(rep(c("X", "Y", "Z"), each = 5), 2)),
    y = c(170, 175, 165, 180, 160,   # biofeedback and drug X
          186, 194, 201, 215, 219,   # biofeedback and drug Y
          180, 187, 199, 170, 204,   # biofeedback and drug Z
          173, 194, 197, 190, 176,   # drug X alone
          189, 194, 217, 206, 199,   # drug Y alone
          202, 228, 190, 206, 224))  # drug Z alone
}

# Kim (2016), Table 3: the 2 (LIGHT) by 4 (MATERIAL) design with 10 per
# cell, N = 80, sums of squares as printed (Type III, which equal the
# sequential ones in this balanced design). The corrected total, 5863.715,
# is the sum of the three effects and the error; the table's uncorrected
# "Total" is 53758.357.
.kim_table3 <- list(
  SS_effect = c(LIGHT = 34.716, MATERIAL = 1996.998,
                "LIGHT:MATERIAL" = 1570.389),
  df_effect = c(1, 3, 3), SS_error = 2261.612, df_error = 72,
  SS_total = 5863.715, N = 80)

# The noncentrality limits of an F statistic, by uniroot() on pf(),
# independently of ci_nc_F().
.ncp_by_pf <- function(F_value, df_1, df_2, alpha_lower = 0.025,
                       alpha_upper = 0.025) {
  gap <- function(ncp, p) stats::pf(F_value, df_1, df_2, ncp = ncp) - p
  lower <- if (stats::pf(F_value, df_1, df_2) < 1 - alpha_lower) 0 else
    stats::uniroot(gap, c(0, 1e4), p = 1 - alpha_lower, tol = 1e-12)$root
  upper <- stats::uniroot(gap, c(0, 1e4), p = alpha_upper, tol = 1e-12)$root
  c(lower, upper)
}
