v <- function(tab, t) tab$value[tab$term == t]

test_that("design_consequences() returns the documented tidy schema", {
  d <- design_consequences(true_effect = 0.1, se = 0.3)
  expect_s3_class(d, "dmar_tbl")
  expect_identical(d$term, c("power", "type_s_error", "exaggeration_ratio",
                             "expected_half_width", "mean_ci_width",
                             "median_ci_width", "sd_ci_width",
                             "pct_ci_less_w", "target_width",
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
                           w = 0.7)
  v <- function(t) d$value[d$term == t]
  df <- 118; se <- sqrt(1/60 + 1/60); tci <- qt(.975, df)
  c_df <- sqrt(2/df) * exp(lgamma((df+1)/2) - lgamma(df/2))
  expect_equal(v("se"), se)
  expect_equal(v("df"), df)
  expect_equal(v("mean_ci_width"), 2 * tci * se * c_df)
  expect_equal(v("expected_half_width"), v("mean_ci_width") / 2)
  expect_equal(v("median_ci_width"), 2 * tci * se * sqrt(qchisq(.5, df)/df))
  expect_equal(v("sd_ci_width"), 2 * tci * se * sqrt(1 - c_df^2))
  expect_equal(v("pct_ci_less_w"), pchisq(df * (0.7 / (2*tci*se))^2, df))
  expect_equal(v("target_width"), 0.7)
})

test_that("the precision lens matches a Monte Carlo simulation of realized widths", {
  skip_on_cran()
  set.seed(113)
  n1 <- 25; n2 <- 25; sdev <- 2; G <- 40000
  d <- design_consequences(true_effect = 0.3, sd = sdev, n_1 = n1, n_2 = n2,
                           w = 1.5)
  v <- function(t) d$value[d$term == t]
  df <- n1 + n2 - 2
  widths <- replicate(G, {
    s2p <- sdev^2 * rchisq(1, df) / df
    2 * qt(.975, df) * sqrt(s2p) * sqrt(1/n1 + 1/n2)
  })
  expect_equal(v("mean_ci_width"),   mean(widths),   tolerance = 0.005)
  expect_equal(v("median_ci_width"), median(widths), tolerance = 0.005)
  expect_equal(v("sd_ci_width"),     sd(widths),     tolerance = 0.02)
  expect_equal(v("pct_ci_less_w"),   mean(widths <= 1.5), tolerance = 0.01)
})

test_that("df = Inf treats the width as known: zero spread, step pct", {
  v <- function(d, t) d$value[d$term == t]
  hi <- design_consequences(true_effect = .2, se = .1, w = 1)
  lo <- design_consequences(true_effect = .2, se = .1, w = .3)
  expect_equal(v(hi, "sd_ci_width"), 0)
  expect_equal(v(hi, "mean_ci_width"), 2 * qnorm(.975) * .1)
  expect_equal(v(hi, "pct_ci_less_w"), 1)
  expect_equal(v(lo, "pct_ci_less_w"), 0)
})

test_that("true_effect = NULL keeps the schema and NAs the significance lens", {
  d <- design_consequences(true_effect = NULL, sd = 1, n_1 = 40, w = 0.6)
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
  expect_error(design_consequences(0.4, sd = 1, n_1 = 60, w = -1),
               "positive target width")
})
