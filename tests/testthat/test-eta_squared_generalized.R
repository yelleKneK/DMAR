test_that("SS interface matches the closed-form formula (scalar SS_observed)", {
  res <- eta_squared_generalized(SS_effect = 100, SS_observed = 70,
                                 SS_error = 200)
  expect_equal(res$eta_squared_generalized, 100 / (100 + 70 + 200),
               tolerance = 1e-12)
  expect_equal(res$effect, "overall")
})

test_that("SS interface accepts a vector of observed SS", {
  res <- eta_squared_generalized(SS_effect = 100, SS_observed = c(40, 30),
                                 SS_error = 200)
  expect_equal(res$eta_squared_generalized, 100 / (100 + 70 + 200),
               tolerance = 1e-12)
})

test_that("F/df interface agrees with SS interface", {
  # Choose SS values, derive F and dfs, check both interfaces agree.
  SS_eff <- 100; df_eff <- 2
  SS_err <- 200; df_err <- 50
  SS_obs_vec <- c(40, 30); df_obs_vec <- c(1, 1)
  F_eff <- (SS_eff / df_eff) / (SS_err / df_err)
  F_obs_vec <- (SS_obs_vec / df_obs_vec) / (SS_err / df_err)

  res_ss  <- eta_squared_generalized(SS_effect = SS_eff,
                                     SS_observed = SS_obs_vec,
                                     SS_error = SS_err)
  res_fdf <- eta_squared_generalized(F_effect = F_eff, df_effect = df_eff,
                                     F_observed = F_obs_vec,
                                     df_observed = df_obs_vec,
                                     df_error = df_err)
  expect_equal(res_ss$eta_squared_generalized,
               res_fdf$eta_squared_generalized, tolerance = 1e-12)
})

test_that("dual interface: consistent inputs proceed silently", {
  SS_eff <- 100; df_eff <- 2
  SS_err <- 200; df_err <- 50
  SS_obs_vec <- c(40, 30); df_obs_vec <- c(1, 1)
  F_eff <- (SS_eff / df_eff) / (SS_err / df_err)
  F_obs_vec <- (SS_obs_vec / df_obs_vec) / (SS_err / df_err)

  res <- eta_squared_generalized(
    SS_effect   = SS_eff,    SS_observed = SS_obs_vec, SS_error = SS_err,
    F_effect    = F_eff,     df_effect   = df_eff,
    F_observed  = F_obs_vec, df_observed = df_obs_vec, df_error = df_err
  )
  expect_equal(res$eta_squared_generalized,
               SS_eff / (SS_eff + sum(SS_obs_vec) + SS_err),
               tolerance = 1e-12)
})

test_that("dual interface: inconsistent inputs stop with a detailed message", {
  expect_error(
    eta_squared_generalized(
      SS_effect = 100, SS_observed = 70, SS_error = 200,
      F_effect  = 1.0, df_effect = 2,
      F_observed = 1, df_observed = 1, df_error = 50
    ),
    "Inconsistent inputs"
  )
})

test_that("model interface (no observed): equals total/partial eta^2", {
  # With no observed factors specified, eta_g = SS_eff / (SS_eff + SS_err)
  # which equals partial eta^2 (and total in one-way).
  fit <- aov(weight ~ group, data = PlantGrowth)
  res_g <- eta_squared_generalized(fit)
  res_p <- eta_squared(fit)
  expect_equal(res_g$eta_squared_generalized, res_p$eta_squared,
               tolerance = 1e-12)
})

test_that("model interface, factorial with observed factor", {
  fit <- aov(len ~ supp * dose, data = ToothGrowth)
  res <- eta_squared_generalized(fit, observed = "supp")
  expect_equal(nrow(res), 3L)
  expect_setequal(res$effect, c("supp", "dose", "supp:dose"))

  tbl <- anova(fit)
  ss_err <- tbl["Residuals", "Sum Sq"]
  ss_supp <- tbl["supp", "Sum Sq"]
  ss_dose <- tbl["dose", "Sum Sq"]
  ss_int  <- tbl["supp:dose", "Sum Sq"]
  # Olejnik and Algina (2003, Eq. 5): listing supp as measured also makes
  # supp:dose a measured source, so both join dose's denominator.
  expected_dose <- ss_dose / (ss_dose + ss_supp + ss_int + ss_err)
  expect_equal(res$eta_squared_generalized[res$effect == "dose"],
               expected_dose, tolerance = 1e-12)

  # For 'supp' as focal: its own SS leads the denominator, and the
  # supp:dose interaction (a measured source) joins it.
  expected_supp <- ss_supp / (ss_supp + ss_int + ss_err)
  expect_equal(res$eta_squared_generalized[res$effect == "supp"],
               expected_supp, tolerance = 1e-12)
})

test_that("model with observed = all effects equals total eta^2 per row", {
  fit <- aov(breaks ~ wool * tension, data = warpbreaks)
  res <- eta_squared_generalized(fit,
                                 observed = c("wool", "tension", "wool:tension"))
  tbl <- anova(fit)
  ss_total <- sum(tbl[, "Sum Sq"])
  for (eff in res$effect) {
    expect_equal(res$eta_squared_generalized[res$effect == eff],
                 tbl[eff, "Sum Sq"] / ss_total,
                 tolerance = 1e-12)
  }
})

test_that("eta_g is bounded in [0, 1]", {
  set.seed(113)
  for (i in seq_len(50)) {
    ss_eff <- runif(1, 0, 500)
    ss_obs <- runif(sample(0:4, 1), 0, 200)
    ss_err <- runif(1, 50, 500)
    val <- eta_squared_generalized(SS_effect = ss_eff,
                                   SS_observed = if (length(ss_obs)) ss_obs else 0,
                                   SS_error = ss_err)$eta_squared_generalized
    expect_gte(val, 0)
    expect_lte(val, 1)
  }
})

test_that("errors on missing arguments / bad input", {
  expect_error(eta_squared_generalized(), "Provide a fitted model")
  expect_error(eta_squared_generalized(SS_effect = 100, SS_error = 200),
               "SS interface requires")
  expect_error(eta_squared_generalized(F_effect = 1, df_effect = 2),
               "F/df interface requires")
  expect_error(eta_squared_generalized(SS_effect = -1, SS_observed = 0,
                                       SS_error = 100),
               "'SS_effect' must be")
  expect_error(eta_squared_generalized(object = list(foo = 1)), "aov, lm, or aovlist")
  fit <- aov(len ~ supp * dose, data = ToothGrowth)
  expect_error(eta_squared_generalized(fit, observed = "not_a_factor"),
               "not in anova")
})

test_that("warn when both object and raw args are supplied", {
  fit <- aov(weight ~ group, data = PlantGrowth)
  expect_warning(
    eta_squared_generalized(object = fit, SS_effect = 100,
                            SS_observed = 0, SS_error = 200),
    "raw SS or F/df arguments are ignored"
  )
})

test_that("the Olejnik and Algina (2003) worked examples reproduce", {
  # Keppel (1991) ABc example, O&A p. 439: eta2_G(A) = 33.63 / 170.03.
  keppel <- eta_squared_generalized(SS_effect = 33.63, SS_observed = 52.40,
                                    SS_error = 84.00)
  expect_equal(keppel$eta_squared_generalized[1],
               33.63 / (33.63 + 52.40 + 84.00), tolerance = 1e-12)
  expect_equal(round(keppel$eta_squared_generalized[1], 3), 0.198)
  # Kirk mixed-design example, O&A p. 441: eta2_G = 3.125 / 21.625 = .145.
  kirk <- eta_squared_generalized(SS_effect = 3.125,
                                  SS_observed = 21.625 - 3.125 - 10.0,
                                  SS_error = 10.0)
  expect_equal(kirk$eta_squared_generalized[1],
               3.125 / 21.625, tolerance = 1e-12)
})

test_that("listing a measured factor pulls its interactions into the denominator", {
  # Olejnik and Algina (2003, Eq. 5): an effect is a measured source when
  # any factor in its term is measured. Listing the bare factor must give
  # the same result as listing the factor and all its interactions.
  set.seed(113)
  d <- expand.grid(A = factor(1:3), B = factor(1:2), c = factor(1:2),
                   rep = 1:5)
  d$y <- rnorm(nrow(d)) + as.numeric(d$A) + 0.5 * as.numeric(d$c)
  fit <- aov(y ~ A * B * c, data = d)
  via_factor   <- eta_squared_generalized(fit, observed = "c")
  via_explicit <- eta_squared_generalized(fit,
    observed = c("c", "A:c", "B:c", "A:B:c"))
  expect_equal(via_factor$eta_squared_generalized,
               via_explicit$eta_squared_generalized, tolerance = 1e-12)
  # And the interaction SS genuinely lands in a manipulated effect's
  # denominator: eta2_G(A) must be smaller than the listed-only variant
  # that ignored the mixed interactions.
  tbl <- anova(fit)
  ss  <- setNames(tbl[, "Sum Sq"], rownames(tbl))
  denom_paper <- ss[["A"]] + ss[["c"]] + ss[["A:c"]] + ss[["B:c"]] +
    ss[["A:B:c"]] + ss[["Residuals"]]
  a_row <- via_factor$eta_squared_generalized[via_factor$effect == "A"]
  expect_equal(a_row, unname(ss[["A"]] / denom_paper), tolerance = 1e-12)
})
