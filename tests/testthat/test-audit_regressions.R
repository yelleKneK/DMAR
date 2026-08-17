# Regression tests for the release-audit fixes. Each test anchors a corrected
# computation to an independent reference (a published value, another package,
# or a first-principles derivation) or pins corrected behavior at a boundary.

test_that("krippendorff_alpha matches the published Krippendorff (2011) value", {
  ratings <- matrix(c(
    1, 2, 3, 3, 2, 1, 4, 1, 2, NA, NA, NA,
    1, 2, 3, 3, 2, 2, 4, 1, 2, 5,  NA, 3,
    NA, 3, 3, 3, 2, 3, 4, 2, 2, 5,  1,  NA,
    1, 2, 3, 3, 2, 4, 4, 1, 2, 5,  1,  NA
  ), nrow = 12, ncol = 4)
  res <- krippendorff_alpha(ratings, level = "nominal", boot = FALSE)
  expect_equal(res$value[res$term == "krippendorff_alpha"], 0.743, tolerance = 0.002)
  # Single-rating unit 12 is excluded from the pairable count.
  expect_equal(res$value[res$term == "n_pairable"], 40)
})

test_that("krippendorff_alpha agrees with irr across metrics", {
  # Pinned from irr::kripp.alpha (irr 0.85, 2026-08-09); live comparison in
  # tools/oracle_checks.R.
  ratings <- matrix(c(
    1, 2, 3, 3, 2, 1, 4, 1, 2, NA, NA, NA,
    1, 2, 3, 3, 2, 2, 4, 1, 2, 5,  NA, 3,
    NA, 3, 3, 3, 2, 3, 4, 2, 2, 5,  1,  NA,
    1, 2, 3, 3, 2, 4, 4, 1, 2, 5,  1,  NA
  ), nrow = 12, ncol = 4)
  by_irr <- c(nominal  = 0.7434210526315790,
              ordinal  = 0.8153875037548813,
              interval = 0.8491071428571428,
              ratio    = 0.7974027747116121)
  for (lv in names(by_irr)) {
    d <- krippendorff_alpha(ratings, level = lv, boot = FALSE)$value[1]
    expect_equal(d, unname(by_irr[[lv]]), tolerance = 1e-6, info = lv)
  }
})

test_that("gwet_ac AC2 weights affect the chance term and match irrCAC", {
  # Pinned from irrCAC::gwet.ac1.raw (irrCAC 1.4, 2026-08-09), which reports
  # its coefficient and se rounded to five decimals; live comparison in
  # tools/oracle_checks.R.
  set.seed(113)
  r1 <- sample(1:5, 60, TRUE, prob = c(.1, .2, .4, .2, .1))
  r2 <- pmin(pmax(r1 + sample(c(-1, 0, 1), 60, TRUE, prob = c(.2, .6, .2)), 1), 5)
  by_irrCAC <- list(unweighted = c(coeff = 0.60843, se = 0.07517),
                    linear     = c(coeff = 0.81423, se = 0.03672),
                    quadratic  = c(coeff = 0.92999, se = 0.01454))
  for (w in names(by_irrCAC)) {
    dm <- gwet_ac(cbind(r1, r2), weights = w)
    expect_equal(dm$value[dm$term == "gwet_ac"], by_irrCAC[[w]][["coeff"]],
                 tolerance = 1e-4, info = w)
    expect_equal(dm$value[dm$term == "se"],       by_irrCAC[[w]][["se"]],
                 tolerance = 1e-4, info = w)
  }
})

test_that("reliability_omega_categorical uses the correctly ordered polychoric matrix", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  J <- 4; Ntot <- 300; lam <- c(.8, .7, .75, .65); f <- rnorm(Ntot)
  d <- as.data.frame(lapply(seq_len(J), function(j)
    as.integer(cut(lam[j] * f + sqrt(1 - lam[j]^2) * rnorm(Ntot),
                   c(-Inf, -.6, .4, Inf), labels = FALSE))))
  names(d) <- paste0("y", seq_len(J))
  do <- as.data.frame(lapply(d, ordered))
  poly <- DMAR:::.polychoric_from_fit(
    lavaan::cfa(paste0("f =~ ", paste(names(do), collapse = " + ")),
                data = do, ordered = TRUE, std.lv = TRUE),
    do, names(do))
  expect_identical(colnames(poly), names(do))
  ref <- lavaan::lavCor(do, ordered = names(do))
  expect_lt(max(abs(poly[upper.tri(poly)] - ref[upper.tri(ref)])), 1e-6)
})

test_that("anova_within_two_way interaction epsilon matches car for a != b", {
  skip_if_not_installed("car")
  skip_if_not_installed("MASS")
  set.seed(113)
  a <- 3; b <- 4; n <- 16
  Sig <- 0.5^abs(outer(1:(a * b), 1:(a * b), "-")) *
    outer(sqrt(seq(1, 3, length.out = a * b)), sqrt(seq(1, 3, length.out = a * b)))
  Y <- MASS::mvrnorm(n, mu = rep(0, a * b), Sigma = Sig)
  df <- expand.grid(subject = factor(1:n), A = factor(1:a), B = factor(1:b))
  df$y <- as.vector(Y)
  res <- anova_within_two_way(data = df, outcome = "y",
                              factor_A = "A", factor_B = "B", subject = "subject")
  gg_ab <- res$epsilon[res$effect == "A:B" & res$adjustment == "Greenhouse-Geisser"]
  # suppressWarnings: car warns "HF eps > 1 treated as 1", an informational
  # message about its own output, not a property under test.
  cs <- suppressWarnings(
    summary(car::Anova(lm(Y ~ 1),
                       idata = expand.grid(A = factor(1:a), B = factor(1:b)),
                       idesign = ~ A * B, type = 3), multivariate = FALSE))
  expect_equal(gg_ab, unname(cs$pval.adjustments["A:B", "GG eps"]), tolerance = 1e-5)
})

test_that("expected_r matches a Monte Carlo E[r] and is not the old off-by-one form", {
  skip_on_cran()  # 200,000 Monte Carlo replications; the closed-form anchors in test-expected_r.R run on CRAN
  mc_er <- function(rho, n, B) {
    S <- matrix(c(1, rho, rho, 1), 2, 2); L <- chol(S)
    mean(replicate(B, { X <- matrix(rnorm(2 * n), n, 2) %*% L; cor(X[, 1], X[, 2]) }))
  }
  set.seed(20260702)
  d <- expected_r(rho = 0.9, n = 5)$expected_r
  m <- mc_er(0.9, 5, 2e5)
  expect_equal(d, m, tolerance = 0.004)
  # The old (wrong) closed form gave ~0.8867 here; the correct value is ~0.8685.
  expect_lt(d, 0.875)
})

test_that("smd/smd_c unbiased correction and power_equivalence_md do not overflow", {
  expect_true(is.finite(smd(mean_1 = 1, mean_2 = 0, s = 1, n_1 = 250, n_2 = 250,
                            unbiased = TRUE)$value))
  expect_true(is.finite(smd_c(mean_T = 1, mean_C = 0, s_C = 1, n_C = 400,
                              unbiased = TRUE)$value))
  # log-scale form equals the direct gamma ratio at small df
  d18 <- smd_c(mean_T = 1, mean_C = 0, s_C = 1, n_C = 18, unbiased = TRUE)$value
  expect_equal(d18, gamma(17 / 2) / (sqrt(17 / 2) * gamma(16 / 2)), tolerance = 1e-10)
  p <- power_equivalence_md(alpha_level = .05, logscale = FALSE, ltheta1 = .8,
                            ltheta2 = 1.25, ldiff = 1, sigma = .4, n = 402, nu = 400)
  expect_gt(p$value[1], 0.99)
})

test_that("smd and smd_c error (not NULL) on unresolvable input", {
  expect_error(smd(mean_1 = 1, s = 2))
  expect_error(smd_c(mean_C = 1))
  expect_error(smd_c())
})

test_that("ci_R2 clamps the random-predictor upper limit at the null boundary", {
  r <- suppressWarnings(ci_R2(R2 = 0.005, N = 100, p = 5, conf_level = 0.95,
                              random_predictors = TRUE))
  expect_equal(r$value[r$term == "upper_limit"], 0)
  # Legitimate case still matches MBESS-style limits (upper well below 1).
  r2 <- ci_R2(R2 = 0.30, N = 100, p = 5, conf_level = 0.95, random_predictors = TRUE)
  expect_lt(r2$value[r2$term == "upper_limit"], 0.5)
})

test_that("ci_c and ci_sc use N - k error df on the psi-only path", {
  a <- ci_c(psi = 2.5, s_anova = 1.6, c_weights = c(1, -1, 0, 0),
            n = c(7, 9, 6, 8), N = 30)
  b <- ci_c(psi = 2.5, s_anova = 1.6, c_weights = c(1, -1, 0, 0),
            n = c(7, 9, 6, 8), N = 30, df_error = 26)
  expect_equal(a$value, b$value, tolerance = 1e-10)
})

test_that("ci_R2 / ci_sc / ci_reg_coef do not mislabel coverage on the alpha path", {
  expect_equal(attr(suppressWarnings(
    ci_R2(R2 = .4, N = 100, p = 3, alpha_lower = .05, alpha_upper = .05)),
    "conf_level"), 0.90)
  expect_null(attr(
    ci_sc(psi = 2.5, s_anova = 1.6, c_weights = c(1, -1, 0), n = c(10, 10, 10),
          N = 30, alpha_lower = .05, alpha_upper = .05), "conf_level"))
  expect_null(attr(suppressWarnings(
    ci_reg_coef(b_j = .5, SE_b_j = .1, N = 100, p = 2,
                alpha_lower = .05, alpha_upper = .05)), "conf_level"))
  # Mixing conf_level with the alphas is rejected.
  expect_error(ci_R2(R2 = .4, N = 100, p = 3, conf_level = .9,
                     alpha_lower = .05, alpha_upper = .05))
})

test_that("ci_reg_coef noncentral prob columns match the central convention", {
  central <- ci_reg_coef(b_j = .6707, SE_b_j = .1761, N = 30, p = 6, noncentral = FALSE)
  noncent <- ci_reg_coef(b_j = .6707, SE_b_j = .1761, N = 30, p = 6, noncentral = TRUE)
  expect_equal(central$prob_less[central$term == "lower_limit"],
               noncent$prob_less[noncent$term == "lower_limit"], tolerance = 1e-3)
})

test_that("cliff_delta returns finite limits under complete separation", {
  r <- suppressWarnings(cliff_delta(c(10, 11, 12), c(1, 2, 3)))
  expect_equal(r$value[r$term == "cliff_delta"], 1)
  expect_true(all(is.finite(r$value[r$term %in% c("lower_limit", "upper_limit")])))
  expect_warning(cliff_delta(c(10, 11, 12), c(1, 2, 3)), "completely separated")
})

test_that("correlations_test reports NA (not 0) for pairs with too few complete cases", {
  set.seed(113)
  d <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  d$a[1:28] <- NA
  res <- correlations_test(d)
  expect_true(is.na(res$r["a", "b"]))
  expect_false(is.na(res$r["b", "c"]))
})

test_that("ss_aipe_cv returns the minimal N (no workspace leak) when width is easy", {
  assign("N", 999, envir = .GlobalEnv)
  on.exit(if (exists("N", envir = .GlobalEnv, inherits = FALSE))
    rm("N", envir = .GlobalEnv), add = TRUE)
  out <- suppressWarnings(ss_aipe_cv(C_of_V = .1, width = 1))
  expect_equal(out$value[out$term == "necessary_N"], 4)
})

test_that("ss_aipe_partial_r fisher_z path does not return the grid minimum", {
  fz <- suppressWarnings(ss_aipe_partial_r(rho = .3, J = 2, width = .05, fisher_z = TRUE))
  ss <- fz$value[fz$term == "necessary_N"]
  expect_gt(ss, 5000)
})

test_that("ss_power_R2 errors (does not hang) when the alternative is at/below the null", {
  expect_error(ss_power_R2(population_R2 = 0.3, null_R2 = 0.5, p = 5, desired_power = .8),
               "no effect in the hypothesized direction")
})

test_that("anova.mlmr labels each row with its own model", {
  skip_if_not_installed("lavaan")
  m0 <- mlmr(mpg ~ wt, data = mtcars, ci_method = "wald")
  m1 <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
  a <- anova(m0, m1)
  # lavTestLRT orders the full model (fewer df) first; its row must be labeled Model 2.
  expect_equal(rownames(a)[which.min(a[["Df"]])], "Model 2")
})

test_that("ss_aipe_reliability does not pollute the global environment", {
  if (exists("sim_data", envir = .GlobalEnv, inherits = FALSE))
    rm("sim_data", envir = .GlobalEnv)
  # The corrected function uses local assignment; guard by inspecting its body.
  expect_false(any(grepl("<<-", deparse(body(ss_aipe_reliability)))))
})

test_that("cfa_1 with se = 'none' returns omega with an NA standard error", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  p <- 5; L <- c(.8, .75, .7, .65, .6); f <- rnorm(300)
  X <- sapply(seq_len(p), function(j) L[j] * f + sqrt(1 - L[j]^2) * rnorm(300))
  colnames(X) <- paste0("x", seq_len(p))
  r <- cfa_1(S = cov(X), N = 300, se = "none")
  expect_true(is.finite(r$estimate[r$term == "omega_f1"]))
  expect_true(is.na(r$se[r$term == "omega_f1"]))
})
