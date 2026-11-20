v <- function(tab, t) tab$value[tab$term == t]

test_that("design_consequences() returns the documented tidy schema", {
  d <- design_consequences(true_effect = 0.1, se = 0.3)
  expect_s3_class(d, "dmar_tbl")
  expect_identical(d$term, c("power", "type_s_error", "exaggeration_ratio",
                             "expected_half_width", "mean_ci_width",
                             "median_ci_width", "sd_ci_width",
                             "pct_ci_less_width", "width",
                             "true_effect", "se", "df", "alpha_level"))
  expect_identical(attr(d, "conf_level"), 0.95)
  expect_type(d$value, "double")
  expect_equal(v(d, "true_effect"), 0.1)
  expect_equal(v(d, "se"), 0.3)
  expect_equal(v(d, "alpha_level"), 0.05)
})

test_that("design_consequences() power matches the two-sided normal power formula", {
  th <- 0.1; s <- 0.3; a <- 0.05
  z  <- qnorm(1 - a / 2); lam <- th / s
  expect_equal(v(design_consequences(th, s, alpha_level = a), "power"),
               (1 - pnorm(z - lam)) + pnorm(-z - lam))
})

test_that("design_consequences() limits behave correctly", {
  # Precise design: power near 1, no sign errors, no exaggeration.
  big <- design_consequences(true_effect = 1, se = 0.05)
  expect_gt(v(big, "power"), 0.999)
  expect_lt(v(big, "type_s_error"), 1e-12)
  expect_equal(v(big, "exaggeration_ratio"), 1, tolerance = 1e-6)
  # Null effect: power = alpha, Type S = 1/2 by symmetry, Type M undefined.
  null <- design_consequences(true_effect = 0, se = 1, alpha_level = 0.05)
  expect_equal(v(null, "power"), 0.05, tolerance = 1e-12)
  expect_equal(v(null, "type_s_error"), 0.5)
  expect_true(is.na(v(null, "exaggeration_ratio")))
  # The sign of the true effect does not change any of the three quantities.
  expect_equal(design_consequences(0.2, 0.5)$value[1:3],
               design_consequences(-0.2, 0.5)$value[1:3])
})

test_that("design_consequences() underpowered designs exaggerate, as in Gelman and Carlin (2014)", {
  d <- design_consequences(true_effect = 0.1, se = 0.3)
  expect_lt(v(d, "power"), 0.10)             # badly underpowered
  expect_gt(v(d, "type_s_error"), 0.01)      # real chance the sign is wrong
  expect_gt(v(d, "exaggeration_ratio"), 3)   # significant estimates inflate severely
})

test_that("design_consequences() finite df reduces toward the normal case as df grows", {
  norm  <- design_consequences(0.4, 0.26)$value[1:3]
  t_big <- design_consequences(0.4, 0.26, df = 1e5)$value[1:3]
  expect_equal(t_big, norm, tolerance = 1e-3)
  # And finite df costs power relative to the normal reference.
  t_small <- design_consequences(0.4, 0.26, df = 10)
  expect_lt(v(t_small, "power"), norm[1])
})

test_that("finite df uses the noncentral t, not retrodesign's shifted central t", {
  # n = 5 per group, d = 1, sd = 1 (se = sqrt(0.4), df = 8). A Monte
  # Carlo of the two-group design (2e6 replications: draw both samples,
  # pool the variance, run the t test) gives power 0.2860 (MCSE 0.0003),
  # Type S 0.00129, and Type M 1.6576 (MCSE 0.0005). The noncentral t
  # (and the chi-mixture Type M) reproduce those: 0.2862955 / 0.0012850 /
  # 1.6575870. The location-shifted central t of retrodesign(), which
  # this function used before the 2026-08 fix, gives 0.2469 / 0.00937 /
  # 1.9150, which the Monte Carlo rules out by over 50 MCSEs on power.
  d <- design_consequences(true_effect = 1, sd = 1, n_1 = 5, n_2 = 5)
  expect_equal(v(d, "power"), 0.2862954934, tolerance = 1e-8)
  expect_equal(v(d, "type_s_error"), 0.0012850209, tolerance = 1e-6)
  expect_equal(v(d, "exaggeration_ratio"), 1.6575869641, tolerance = 1e-7)
  # Independent anchor: base R's power.t.test() computes two-group
  # noncentral t power, and strict = TRUE adds the wrong-sign rejection
  # region, exactly the power reported here.
  expect_equal(v(d, "power"),
               power.t.test(n = 5, delta = 1, sd = 1, sig.level = 0.05,
                            strict = TRUE)$power,
               tolerance = 1e-12)
})

test_that("the noncentral t agrees with retrodesign's regime at large df", {
  # retrodesign()'s location-shifted central t is the known-se
  # approximation; at df = 5000 the two sampling models nearly coincide
  # (power 0.35244 vs 0.35250), while at df = 8 they do not (0.2469 vs
  # 0.2863; the regression test above). This ties the implementation
  # back to Gelman and Carlin (2014) where their approximation holds.
  L <- 0.4 / 0.26; df <- 5000; crit <- qt(0.975, df)
  retro_power  <- (1 - pt(crit - L, df)) + pt(-crit - L, df)
  retro_type_s <- pt(-crit - L, df) / retro_power
  retro_e_abs  <- integrate(function(t) (L + t) * dt(t, df),
                            lower = crit - L, upper = Inf,
                            rel.tol = 1e-10)$value +
                  integrate(function(t) -(L + t) * dt(t, df),
                            lower = -Inf, upper = -crit - L,
                            rel.tol = 1e-10)$value
  retro_type_m <- retro_e_abs / (retro_power * L)
  d <- design_consequences(true_effect = 0.4, se = 0.26, df = df)
  expect_equal(v(d, "power"), retro_power, tolerance = 5e-4)
  expect_equal(v(d, "exaggeration_ratio"), retro_type_m, tolerance = 5e-4)
  # Type S is a ratio of small tails and converges more slowly in
  # relative terms; check the relative difference explicitly.
  expect_lt(abs(v(d, "type_s_error") - retro_type_s) / retro_type_s, 0.01)
})

test_that("design_consequences() validates its arguments", {
  expect_error(design_consequences("a", 1), "single number")
  expect_error(design_consequences(0.2, -1), "positive")
  expect_error(design_consequences(0.2, 1, alpha_level = 1.2), "\\(0, 1\\)")
  expect_error(design_consequences(0.2, 1, df = -3), "positive")
})

test_that("design_consequences() matches a Monte Carlo simulation (normal and t)", {
  skip_on_cran()
  # Fast anchors that stay on CRAN: the pinned noncentral t values and the
  # power.t.test(strict = TRUE) agreement in the regression test above.
  set.seed(113)
  check_mc <- function(th, s, df, G = 4e5) {
    # Simulate the design the function models: estimate ~ N(theta, se^2)
    # with an independently estimated standard error se_hat =
    # se * sqrt(W / df), W ~ chi-square_df, and the t test
    # |estimate / se_hat| > crit. (This test previously drew
    # est = th + s * rt(G, df), the location-shifted central t of
    # retrodesign(); that is not the distribution a t test's estimate
    # and estimated se jointly follow, and the function now models the
    # estimated-se design exactly, so the simulation was corrected with
    # the code.)
    est    <- th + s * rnorm(G)
    crit   <- if (is.finite(df)) qt(0.975, df) else qnorm(0.975)
    se_hat <- if (is.finite(df)) s * sqrt(rchisq(G, df) / df) else s
    sig  <- abs(est) > crit * se_hat
    ana  <- design_consequences(th, s, df = df)
    expect_equal(v(ana, "power"), mean(sig), tolerance = 0.01)
    expect_equal(v(ana, "type_s_error"),
                 mean(sign(est[sig]) != sign(th)), tolerance = 0.02)
    expect_equal(v(ana, "exaggeration_ratio"),
                 mean(abs(est[sig])) / abs(th), tolerance = 0.02)
  }
  check_mc(0.1, 0.3, Inf)       # the Gelman-Carlin normal case
  check_mc(0.4, 0.26, 58)       # a two-group t-test design
  check_mc(1, sqrt(0.4), 8)     # n = 5 per group, d = 1: small df
})

test_that("the n = 5 / d = 1 pins match a raw-data two-group simulation", {
  skip_on_cran()
  # Fast anchor that stays on CRAN: the pinned noncentral t values and
  # the power.t.test(strict = TRUE) agreement in the regression test
  # above. This simulation draws the actual data (no model shortcut), so
  # it also arbitrates the (sd, n_1, n_2) -> (se, df) mapping.
  set.seed(113)
  n1 <- 5; n2 <- 5; theta <- 1; G <- 2e5
  df <- n1 + n2 - 2; crit <- qt(0.975, df)
  g1 <- matrix(rnorm(G * n1, mean = theta), nrow = G)
  g2 <- matrix(rnorm(G * n2, mean = 0), nrow = G)
  m1 <- rowMeans(g1); m2 <- rowMeans(g2)
  ss1 <- rowSums(g1^2) - n1 * m1^2
  ss2 <- rowSums(g2^2) - n2 * m2^2
  est    <- m1 - m2
  se_hat <- sqrt((ss1 + ss2) / df) * sqrt(1 / n1 + 1 / n2)
  sig <- abs(est / se_hat) > crit
  d <- design_consequences(true_effect = theta, sd = 1, n_1 = n1, n_2 = n2)
  expect_equal(v(d, "power"), mean(sig), tolerance = 0.01)
  expect_equal(v(d, "type_s_error"),
               mean(sign(est[sig]) != sign(theta)), tolerance = 0.02)
  expect_equal(v(d, "exaggeration_ratio"),
               mean(abs(est[sig])) / theta, tolerance = 0.01)
})


test_that("the precision lens has the right closed forms and identities", {
  d <- design_consequences(true_effect = 0.4, sd = 1, n_1 = 60, n_2 = 60,
                           width = 0.7)
  v <- function(t) d$value[d$term == t]
  df <- 118; se <- sqrt(1/60 + 1/60); tci <- qt(.975, df)
  c_df <- sqrt(2/df) * exp(lgamma((df+1)/2) - lgamma(df/2))
  expect_equal(v("se"), se)
  expect_equal(v("df"), df)
  expect_equal(v("mean_ci_width"), 2 * tci * se * c_df)
  expect_equal(v("expected_half_width"), v("mean_ci_width") / 2)
  expect_equal(v("median_ci_width"), 2 * tci * se * sqrt(qchisq(.5, df)/df))
  expect_equal(v("sd_ci_width"), 2 * tci * se * sqrt(1 - c_df^2))
  expect_equal(v("pct_ci_less_width"), pchisq(df * (0.7 / (2*tci*se))^2, df))
  expect_equal(v("width"), 0.7)
})

test_that("ordinary designs keep their precision rows to 1e-12", {
  # Through 1.0.x the mean factor E(S) / sigma was the difference of two
  # lgamma() values below, and sd_ci_width took the square root of one
  # minus its square. Where those expressions are accurate (the mean rows
  # through hundreds of df, the sd row below about 60 df) the rows are
  # unchanged to 1e-12. No other row involves the factor.
  old_c <- function(df) sqrt(2 / df) * exp(lgamma((df + 1) / 2) - lgamma(df / 2))
  for (df in c(0.5, 1, 2.5, 8, 30, 39, 48, 58, 118, 500)) {
    d <- design_consequences(true_effect = 0.4, se = 0.26, df = df, width = 1)
    fixed <- 2 * qt(0.975, df) * 0.26
    expect_equal(v(d, "mean_ci_width"), fixed * old_c(df), tolerance = 1e-12,
                 info = paste("df =", df))
    expect_equal(v(d, "expected_half_width"), fixed * old_c(df) / 2,
                 tolerance = 1e-12, info = paste("df =", df))
    if (df < 60) {
      expect_equal(v(d, "sd_ci_width"), fixed * sqrt(1 - old_c(df)^2),
                   tolerance = 1e-12, info = paste("df =", df))
    }
  }
  # From about 100 df on, the former sd row lost digits to the subtraction:
  # it was off by 1e-12 at 118 df (60 per group, as on the help page), 2e-10
  # at 500, 1e-8 at 5000, and 8e-6 at 1e5. It is now right to rounding,
  # here against Var(S) / sigma^2 evaluated with Rmpfr at 320 bits.
  df_big <- c(118, 500, 5000, 1e5)
  var_big <- c(0.0042282729961787635, 0.00099949950062787167,
               9.9994999500062527e-05, 4.9999874999375004e-06)
  for (i in seq_along(df_big)) {
    d <- design_consequences(true_effect = 0.4, se = 0.26, df = df_big[i])
    fixed <- 2 * qt(0.975, df_big[i]) * 0.26
    expect_equal(v(d, "sd_ci_width"), fixed * sqrt(var_big[i]),
                 tolerance = 1e-14, info = paste("df =", df_big[i]))
  }
})

test_that("the precision lens stays right at a very large df", {
  # Through 1.0.x the mean factor was 1 or more at about half of all df
  # above 1e8 (first at 10810066 df), so mean_ci_width could exceed the
  # width with the standard error known and sd_ci_width was 0; elsewhere
  # sd_ci_width was 57 times its value (1e9 df). E(S) / sigma and
  # Var(S) / sigma^2 below were evaluated with Rmpfr at 320 bits.
  df_big <- c(1e6, 1e8, 1e12)
  mean_big <- c(0.99999975000003127, 0.99999999750000002, 0.99999999999974998)
  var_big <- c(4.9999987499993748e-07, 4.9999999874999999e-09,
               4.9999999999987498e-13)
  for (i in seq_along(df_big)) {
    d <- design_consequences(true_effect = NULL, se = 1, df = df_big[i])
    fixed <- 2 * qt(0.975, df_big[i])
    expect_lt(v(d, "mean_ci_width"), fixed)
    expect_gt(v(d, "sd_ci_width"), 0)
    expect_equal(v(d, "mean_ci_width"), fixed * mean_big[i],
                 tolerance = 4e-15, info = paste("df =", df_big[i]))
    expect_equal(v(d, "sd_ci_width"), fixed * sqrt(var_big[i]),
                 tolerance = 1e-14, info = paste("df =", df_big[i]))
  }
  # Formed through lbeta(), the mean factor came out a few units in the
  # last place above 1 at some df from about 1.7e14 on, which put
  # mean_ci_width above the fixed width (by 2e-15 at 1e15 df). It is now
  # never above 1, and exactly 1 only where the true factor rounds to 1.
  for (df in c(1e14, 1.7e14, 1e15, 2.25e15, 1e16, 1e20)) {
    d <- design_consequences(true_effect = NULL, se = 1, df = df)
    fixed <- 2 * qt(0.975, df)
    expect_lte(v(d, "mean_ci_width"), fixed)
    expect_lte(v(d, "expected_half_width"), fixed / 2)
    expect_gt(v(d, "sd_ci_width"), 0)
  }
  d <- design_consequences(true_effect = NULL, se = 1, df = 1e15)
  expect_equal(v(d, "mean_ci_width"), 2 * qt(0.975, 1e15) * 0.99999999999999975,
               tolerance = 1e-15)
})

test_that("the precision lens is unchanged at ordinary df by the stable factors", {
  # The mean factor through lbeta() was right below 100 df to a few units
  # in the last place, so the mean rows there move by less than 1e-14.
  # The variance factor of 1 minus its square was off by up to 3e-13 below
  # 100 df, so sd_ci_width moves by less than 1e-12.
  lbeta_mean <- function(df) sqrt(2 / df) * exp(lgamma(0.5) - lbeta(df / 2, 0.5))
  for (df in c(2, 5, 18, 48, 99)) {
    d <- design_consequences(true_effect = 0.4, se = 0.26, df = df)
    fixed <- 2 * qt(0.975, df) * 0.26
    expect_equal(v(d, "mean_ci_width"), fixed * lbeta_mean(df),
                 tolerance = 1e-14, info = paste("df =", df))
    expect_equal(v(d, "sd_ci_width"), fixed * sqrt(1 - lbeta_mean(df)^2),
                 tolerance = 1e-12, info = paste("df =", df))
  }
})

test_that("the precision lens matches a Monte Carlo simulation of realized widths", {
  skip_on_cran()
  set.seed(113)
  n1 <- 25; n2 <- 25; sdev <- 2; G <- 40000
  d <- design_consequences(true_effect = 0.3, sd = sdev, n_1 = n1, n_2 = n2,
                           width = 1.5)
  v <- function(t) d$value[d$term == t]
  df <- n1 + n2 - 2
  widths <- replicate(G, {
    s2p <- sdev^2 * rchisq(1, df) / df
    2 * qt(.975, df) * sqrt(s2p) * sqrt(1/n1 + 1/n2)
  })
  expect_equal(v("mean_ci_width"),   mean(widths),   tolerance = 0.005)
  expect_equal(v("median_ci_width"), median(widths), tolerance = 0.005)
  expect_equal(v("sd_ci_width"),     sd(widths),     tolerance = 0.02)
  expect_equal(v("pct_ci_less_width"),   mean(widths <= 1.5), tolerance = 0.01)
})

test_that("df = Inf treats the width as known: zero spread, step pct", {
  v <- function(d, t) d$value[d$term == t]
  hi <- design_consequences(true_effect = .2, se = .1, width = 1)
  lo <- design_consequences(true_effect = .2, se = .1, width = .3)
  expect_equal(v(hi, "sd_ci_width"), 0)
  expect_equal(v(hi, "mean_ci_width"), 2 * qnorm(.975) * .1)
  expect_equal(v(hi, "pct_ci_less_width"), 1)
  expect_equal(v(lo, "pct_ci_less_width"), 0)
})

test_that("true_effect = NULL keeps the schema and NAs the significance lens", {
  d <- design_consequences(true_effect = NULL, sd = 1, n_1 = 40, width = 0.6)
  v <- function(t) d$value[d$term == t]
  expect_true(all(is.na(c(v("power"), v("type_s_error"),
                          v("exaggeration_ratio"), v("true_effect")))))
  expect_false(is.na(v("mean_ci_width")))
  expect_equal(v("df"), 39)            # one-sample derivation: n - 1
  expect_equal(v("se"), 1 / sqrt(40))
  expect_equal(nrow(d), 13)            # the schema never changes shape
})

test_that("the (sd, n) route and the se route agree, and cannot be mixed", {
  a <- design_consequences(0.4, se = sqrt(2/60), df = 118)
  b <- design_consequences(0.4, sd = 1, n_1 = 60, n_2 = 60)
  expect_equal(a$value[1:9], b$value[1:9], tolerance = 1e-12)
  expect_error(design_consequences(0.4, se = 0.2, sd = 1, n_1 = 60),
               "not both")
  expect_error(design_consequences(0.4), "Supply 'se'")
  expect_error(design_consequences(0.4, sd = 1, n_1 = 60, width = -1),
               "positive target width")
})

test_that("the significance lens is exact past a noncentrality of 37.62", {
  # With an estimated standard error the statistic is noncentral t. At
  # lambda = 40, df = 2, and alpha .001 the power is .7981439583; R's pt()
  # (a normal approximation past 37.62; DMAR through 3338679) gave
  # .7823613465, a Type S error of 1.2e-5 (it is below 1e-300), and an
  # exaggeration ratio of 1.0207 (1.0005).
  d <- expect_no_warning(design_consequences(true_effect = 40, se = 1, df = 2,
                                             alpha_level = .001))
  v <- function(t) d$value[d$term == t]
  crit <- qt(1 - .001 / 2, 2)
  expect_equal(v("power"), .oracle_nct(crit, 2, 40, lower_tail = FALSE) +
                 .oracle_nct(-crit, 2, 40), tolerance = 1e-12)
  expect_equal(v("power"), 0.7981439583, tolerance = 1e-9)
  expect_lt(v("type_s_error"), 1e-300)
  expect_equal(v("exaggeration_ratio"), 1.000505561, tolerance = 1e-8)
  # A negative effect: the same power, by symmetry.
  e <- design_consequences(true_effect = -40, se = 1, df = 2, alpha_level = .001)
  expect_equal(e$value[1], v("power"), tolerance = 1e-14)
  # Inside the range the power is pt()'s.
  d <- design_consequences(true_effect = 3, se = 1, df = 5)
  crit <- qt(.975, 5)
  expect_equal(d$value[1], pt(crit, 5, 3, lower.tail = FALSE) + pt(-crit, 5, 3),
               tolerance = 1e-14)
})
