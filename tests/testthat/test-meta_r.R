v <- function(tab, t) tab$value[tab$term == t]

# Studies that differ in correlation, both reliabilities, and size (F01).
.rel_studies <- function() {
  list(r   = c(.21, .35, .42, .18, .30, .47),
       n   = c(48, 120, 75, 300, 60, 150),
       rxx = c(.72, .85, .78, .90, .65, .80),
       ryy = c(.80, .70, .88, .75, .82, .90))
}

# The corrected Fisher Z values and their sampling variances, computed
# here on another path than the one meta_r() takes: the derivative of
# atanh(r / sqrt(a)) with respect to atanh(r) by a central difference in
# atanh(r), times the variance 1 / (n - 3) of atanh(r).
.corrected_yivi_oracle <- function(r, n, a, h = 1e-6) {
  z  <- atanh(r)
  dz <- (atanh(tanh(z + h) / sqrt(a)) - atanh(tanh(z - h) / sqrt(a))) /
    (2 * h)
  list(yi = atanh(r / sqrt(a)), vi = dz^2 / (n - 3))
}

test_that("meta_r() matches metafor's Fisher's Z pipeline, back-transformed", {
  r <- c(.28, .35, .22, .40, .31)
  n <- c(120, 85, 200, 60, 150)
  ours <- meta_r(r, n, hartung_knapp = FALSE)
  # Pinned from metafor::rma(yi = atanh(r), vi = 1 / (n - 3), method = "REML"),
  # back-transformed with tanh (metafor 5.0.1, 2026-08-09); live comparison in
  # tools/oracle_checks.R.
  expect_equal(v(ours, "estimate"), 0.2897158161550648, tolerance = 1e-6)
  expect_equal(v(ours, "ci_lower"), 0.2148407135896157, tolerance = 1e-6)
  expect_equal(v(ours, "ci_upper"), 0.3612051691843595, tolerance = 1e-6)
  expect_equal(v(ours, "tau2"), 0, tolerance = 1e-5)
})

test_that("meta_r() reliability corrections disattenuate before pooling", {
  r <- c(.28, .35, .22, .40, .31)
  n <- c(120, 85, 200, 60, 150)
  plain <- meta_r(r, n, hartung_knapp = FALSE)
  fixed <- meta_r(r, n, reliability_y = 0.80, hartung_knapp = FALSE)
  # Correcting for criterion unreliability raises the pooled validity.
  expect_gt(v(fixed, "estimate"), v(plain, "estimate"))
  # This test also compared the corrected analysis with an analysis of the
  # hand-corrected correlations as if they were observed ones, which is
  # the earlier variance, 1 / (n - 3), that finding F01 showed to be too
  # small. The comparison is now with the corrected Fisher Z values and
  # their delta method variances (the tests that follow); per-study
  # reliabilities are honored elementwise.
  rel <- c(.7, .8, .75, .9, .85)
  el <- meta_r(r, n, reliability_x = rel, method = "fe",
               hartung_knapp = FALSE)
  o <- .corrected_yivi_oracle(r, n, rel)
  w <- 1 / o$vi
  expect_equal(v(el, "estimate"), tanh(sum(w * o$yi) / sum(w)),
               tolerance = 1e-10)
  expect_equal(v(el, "se"), sqrt(1 / sum(w)), tolerance = 1e-8)
})

test_that("meta_r() corrects the sampling variance of each corrected correlation, not only the correlation", {
  # Five studies with r = .40 and n = 100, both reliabilities .80: each
  # corrected correlation is .50, and the derivative of atanh(r / .8) with
  # respect to atanh(r) is (1 - .16) / (.8 * (1 - .25)) = 1.4, so each
  # sampling variance is 1.96 / 97. The variance of an observed
  # correlation, 1 / 97, was used before, which gave a standard error of
  # 0.0454 and the interval .4303 to .5637 in place of .0636 and .4009 to
  # .5875.
  fit <- meta_r(rep(.4, 5), rep(100, 5), reliability_x = .8,
                reliability_y = .8, method = "fe", hartung_knapp = FALSE)
  se <- sqrt(1.96 / 97 / 5)
  expect_equal(v(fit, "estimate"), .5, tolerance = 1e-12)
  expect_equal(v(fit, "se"), se, tolerance = 1e-12)
  expect_equal(v(fit, "se"), 0.0635707252861, tolerance = 1e-10)
  expect_equal(c(v(fit, "ci_lower"), v(fit, "ci_upper")),
               tanh(atanh(.5) + c(-1, 1) * qnorm(.975) * se),
               tolerance = 1e-12)
  expect_equal(c(v(fit, "ci_lower"), v(fit, "ci_upper")),
               c(0.400890762639, 0.587541066273), tolerance = 1e-10)
})

test_that("the corrected effect sizes and variances are the delta method values, and exactly the uncorrected ones without reliabilities", {
  s <- .rel_studies()
  a <- s$rxx * s$ryy
  ours <- .meta_r_yivi(s$r, s$n, s$rxx, s$ryy)
  oracle <- .corrected_yivi_oracle(s$r, s$n, a)
  expect_equal(ours$yi, oracle$yi, tolerance = 1e-14)
  expect_equal(ours$vi, oracle$vi, tolerance = 1e-8)
  # The closed form of the variance.
  expect_equal(ours$vi,
               (1 - s$r^2)^2 / (a * (1 - s$r^2 / a)^2 * (s$n - 3)),
               tolerance = 1e-14)
  # Every corrected variance exceeds 1 / (n - 3).
  expect_true(all(ours$vi > 1 / (s$n - 3)))
  # Reliabilities of 1 change nothing, to the last bit.
  plain <- .meta_r_yivi(s$r, s$n, rep(1, 6), rep(1, 6))
  expect_identical(plain$yi, atanh(s$r))
  expect_identical(plain$vi, 1 / (s$n - 3))
  for (m in c("reml", "pm", "dl", "fe")) {
    expect_identical(
      meta_r(s$r, s$n, reliability_x = 1, reliability_y = 1, method = m),
      meta_r(s$r, s$n, method = m), label = m)
  }
})

test_that("meta_r() with reliabilities agrees with metafor given the corrected effect sizes and their variances", {
  skip_on_cran()
  # Every row of the table, under each estimator and with and without the
  # Hartung-Knapp adjustment: the weights, the pooled estimate with its
  # standard error and interval, the heterogeneity quantities with their
  # intervals, and the prediction interval all use the corrected variances.
  # Pinned from metafor::rma(yi, vi, method = toupper(m), test = "z" or
  # "knha", control = list(threshold = 1e-12, tol = 1e-14)) and confint()
  # of the same fit (control = list(tol = 1e-14)), with yi and vi from
  # .corrected_yivi_oracle() above, not from DMAR (metafor 5.0.1,
  # 2026-11-11); live comparison in tools/oracle_snippets/metafor.R. The
  # estimate and the confidence and prediction limits are back-transformed
  # with tanh; the prediction limits are those of Higgins, Thompson, and
  # Spiegelhalter (2009), t on k - 2 df, from metafor's tau2 and SE. Q and
  # the Q-profile limits do not depend on the estimator.
  common <- c(Q = 12.04502945983407, Q_p = 0.03417595524031483,
              tau2_ci_lower = 0, tau2_ci_upper = 0.130951255283132,
              I2_ci_lower = 0, I2_ci_upper = 88.87527909601774)
  pinned <- list(
    fe_z = c(estimate = 0.3677201263850035, se = 0.04852893982218589,
             stat = 7.94956947320182, p_value = 1.871608001504996e-15,
             ci_lower = 0.2827506319264864, ci_upper = 0.4469634227478088,
             tau2 = 0),
    dl_z = c(estimate = 0.4064241598836317, se = 0.08424125908295793,
             stat = 5.120060670335484, p_value = 3.054373767644048e-07,
             ci_lower = 0.260095211140328, ci_upper = 0.5345045019890109,
             tau2 = 0.02309567298856977, I2 = 58.48910111284292,
             H2 = 2.409005891966813,
             prediction_lower = -0.05106808400641398,
             prediction_upper = 0.7229288301035528),
    dl_knha = c(estimate = 0.4064241598836317, se = 0.06903749996983385,
                stat = 6.247624227972644, p_value = 0.001539365890601472,
                ci_lower = 0.2485378708520765, ci_upper = 0.5432726035689651,
                tau2 = 0.02309567298856977, I2 = 58.48910111284292,
                H2 = 2.409005891966813,
                prediction_lower = -0.03210918333869466,
                prediction_upper = 0.71373729517804),
    reml_z = c(estimate = 0.4051365270071033, se = 0.0799989899065425,
               stat = 5.372304251820107, p_value = 7.773679830890336e-08,
               ci_lower = 0.2663990531197422, ci_upper = 0.5274271581594645,
               tau2 = 0.01928990108841982, I2 = 54.0615518662791,
               H2 = 2.17682581938582,
               prediction_lower = -0.01522905933172964,
               prediction_upper = 0.7037986627195194),
    reml_knha = c(estimate = 0.4051365270071033, se = 0.06960917416176411,
                  stat = 6.174170557137668, p_value = 0.001623631347270321,
                  ci_lower = 0.2457107760454, ci_upper = 0.5432219165035099,
                  tau2 = 0.01928990108841982, I2 = 54.0615518662791,
                  H2 = 2.17682581938582,
                  prediction_lower = -0.001557163747622605,
                  prediction_upper = 0.6968316989781346),
    pm_z = c(estimate = 0.4009094267417121, se = 0.07108206670633549,
             stat = 5.975234967648932, p_value = 2.297587842786971e-09,
             ci_lower = 0.2779081981470288, ci_upper = 0.5109764329583364,
             tau2 = 0.01217215731468523, I2 = 42.61419125162129,
             H2 = 1.742591107123244,
             prediction_lower = 0.06026895544691285,
             prediction_upper = 0.6579113712754049),
    pm_knha = c(estimate = 0.4009094267417121, se = 0.07108206670633549,
                stat = 5.975234967648932, p_value = 0.001880612675841696,
                ci_lower = 0.2373931710500061, ci_upper = 0.5423326475770627,
                tau2 = 0.01217215731468523, I2 = 42.61419125162129,
                H2 = 1.742591107123244,
                prediction_lower = 0.06026895544691285,
                prediction_upper = 0.6579113712754049)
  )
  s <- .rel_studies()
  o <- .corrected_yivi_oracle(s$r, s$n, s$rxx * s$ryy)
  for (lab in names(pinned)) {
    m  <- sub("_.*", "", lab)
    hk <- grepl("knha", lab)
    ref <- pinned[[lab]]
    ours <- meta_r(s$r, s$n, reliability_x = s$rxx, reliability_y = s$ryy,
                   method = m, hartung_knapp = hk)
    # The weights: the estimate is the weighted mean of the corrected
    # Fisher Z values with weights 1 / (v_i + tau2), v_i the oracle's.
    w <- 1 / (o$vi + ref[["tau2"]])
    expect_equal(v(ours, "estimate"), tanh(sum(w * o$yi) / sum(w)),
                 tolerance = 1e-7, label = lab)
    expect_equal(v(ours, "estimate"), ref[["estimate"]], tolerance = 1e-7,
                 label = lab)
    expect_equal(v(ours, "se"), ref[["se"]], tolerance = 1e-6, label = lab)
    expect_equal(v(ours, "test_statistic"), ref[["stat"]],
                 tolerance = 1e-6, label = lab)
    expect_equal(v(ours, "p_value"), ref[["p_value"]], tolerance = 1e-5,
                 label = lab)
    expect_equal(c(v(ours, "ci_lower"), v(ours, "ci_upper")),
                 unname(ref[c("ci_lower", "ci_upper")]), tolerance = 1e-6,
                 label = lab)
    expect_equal(v(ours, "tau2"), ref[["tau2"]], tolerance = 1e-6,
                 label = lab)
    expect_equal(v(ours, "Q"), common[["Q"]], tolerance = 1e-10, label = lab)
    expect_equal(v(ours, "Q_p"), common[["Q_p"]], tolerance = 1e-8,
                 label = lab)
    if (m == "fe") next
    for (q in c("I2", "H2", "prediction_lower", "prediction_upper")) {
      expect_equal(v(ours, q), ref[[q]], tolerance = 1e-6,
                   label = paste(lab, q))
    }
    for (q in c("tau2_ci_lower", "tau2_ci_upper", "I2_ci_lower", "I2_ci_upper")) {
      expect_equal(v(ours, q), common[[q]], tolerance = 1e-9,
                   label = paste(lab, q))
    }
  }
})

test_that("intervals for the correlation between constructs cover at about the nominal rate when the reliabilities are known", {
  # Studies of observed scores X = sqrt(rxx) T_x + sqrt(1 - rxx) E_x (and
  # likewise Y) with true scores correlated .5, so the observed correlation
  # is .5 sqrt(rxx ryy). Five studies of 100 with both reliabilities .8,
  # as in the first test of the corrected variance. The interval covers .5
  # in 385 of these 400 samples; with the earlier variance, 1 / (n - 3), it
  # covered .5 in 328 (in theory 84 percent, the interval being 1 / 1.4 of
  # the width it needs).
  skip_on_cran()
  set.seed(20261111)
  G <- 400
  covered <- logical(G)
  rel <- .8
  for (g in seq_len(G)) {
    r <- vapply(1:5, function(i) {
      tx <- rnorm(100)
      ty <- .5 * tx + sqrt(1 - .5^2) * rnorm(100)
      stats::cor(sqrt(rel) * tx + sqrt(1 - rel) * rnorm(100),
                 sqrt(rel) * ty + sqrt(1 - rel) * rnorm(100))
    }, numeric(1))
    fit <- meta_r(r, rep(100, 5), reliability_x = rel, reliability_y = rel,
                  method = "fe", hartung_knapp = FALSE)
    covered[g] <- v(fit, "ci_lower") < .5 && .5 < v(fit, "ci_upper")
  }
  # 400 samples: the Monte Carlo standard error at .95 is about .011.
  expect_gt(mean(covered), .92)
  expect_lt(mean(covered), .98)
})


test_that("meta_r() validates its arguments", {
  r <- c(.28, .35, .22); n <- c(50, 60, 70)
  expect_error(meta_r(c(.5, 1), c(50, 50)), "\\(-1, 1\\)")
  expect_error(meta_r(r, c(50, 60)), "for each study")
  expect_error(meta_r(r, n, reliability_x = c(.8, .9)), "length 1 or one per")
  expect_error(meta_r(c(.9, .95), c(50, 50), reliability_x = 0.5),
               "cannot support")
})

test_that("meta_r() warns when the correlations are identical under Hartung-Knapp", {
  expect_warning(res <- meta_r(r = c(.3, .3, .3), n = c(40, 60, 90)),
                 class = "dmar_meta_hk_zero_se")
  expect_equal(v(res, "estimate"), 0.3)
  expect_equal(c(v(res, "ci_lower"), v(res, "ci_upper")), c(0.3, 0.3))
  expect_no_warning(
    plain <- meta_r(r = c(.3, .3, .3), n = c(40, 60, 90),
                    hartung_knapp = FALSE))
  expect_lt(v(plain, "ci_lower"), 0.3)
  expect_gt(v(plain, "ci_upper"), 0.3)
})

test_that("meta_r() refuses an infinite sample size", {
  # It stopped with "missing value where TRUE/FALSE needed" or, under the
  # fixed effect model, returned an estimate of NaN with no message.
  for (m in c("reml", "fe")) {
    expect_error(meta_r(c(0.1, 0.2, 0.3), c(20, Inf, 20), method = m),
                 "'n' must give an integer sample size")
  }
})

test_that("meta_r() states the conventional interval as correlations and labels the scale of the standard error", {
  r <- c(.3, .3, .3)
  n <- c(40, 60, 90)
  w <- expect_warning(res <- meta_r(r = r, n = n),
                      class = "dmar_meta_hk_zero_se")
  se_z <- sqrt(1 / sum(n - 3))
  expect_identical(w$n_studies, 3L)
  expect_equal(w$estimate, 0.3)
  expect_equal(w$se_conventional, se_z)
  expect_equal(c(w$lower_conventional, w$upper_conventional),
               tanh(atanh(0.3) + c(-1, 1) * qnorm(0.975) * se_z))
  plain <- meta_r(r = r, n = n, hartung_knapp = FALSE)
  expect_equal(c(w$lower_conventional, w$upper_conventional),
               c(v(plain, "ci_lower"), v(plain, "ci_upper")))
  expect_match(conditionMessage(w), "has no width (0.3, 0.3)", fixed = TRUE)
  expect_match(conditionMessage(w), "is 0.07433 (on the Fisher Z scale), with",
               fixed = TRUE)
  expect_match(conditionMessage(w), "normal interval 0.1624 to 0.4262",
               fixed = TRUE)
  # The attributes are on the Fisher Z scale, as the se row is.
  expect_identical(attr(res, "hartung_knapp_q"), 0)
  expect_equal(attr(res, "se_conventional"), se_z)
  expect_identical(nrow(res), 20L)
  fit <- meta_r(c(.28, .35, .22, .40, .31), c(120, 85, 200, 60, 150))
  expect_equal(v(fit, "se"),
               sqrt(attr(fit, "hartung_knapp_q")) * attr(fit, "se_conventional"))
})
