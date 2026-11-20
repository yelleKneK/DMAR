## ss_aipe_rc() / ss_aipe_src() -- sample size for unstandardized / standardized
## regression coefficient AIPE planning.

.src_R <- matrix(c(1, .3, .3, 1), 2)
.src_r <- c(.4, .3)
.src_R5 <- local({ R <- matrix(.4, 5, 5); diag(R) <- 1; R })
.src_c <- function(...) ss_aipe_src(width = .3, true_cor_XX = .src_R,
                                    true_cor_YX = .src_r, which_predictor = 1,
                                    noncentral = FALSE, ...)

test_that("ss_aipe_src(noncentral = FALSE) plans the expected width, and the chi square plan, of Kelley and Maxwell (2003)", {
  # Without assurance the plan is Equation 2 with the t quantile at the
  # expected squared multiple correlations, the central plan of
  # ss_aipe_reg_coef() at unit standard deviations.
  expect_identical(.src_c()$value,
                   ss_aipe_reg_coef(width = .3, true_cor_XX = .src_R,
                                    true_cor_YX = .src_r, which_predictor = 1)$value)
  expect_equal(.src_c()$value, 156)
  # With assurance and the default assurance_method, "analytic", their
  # Equation 3 (the chi square factor): 181 here, and 486 for their five
  # exchangeable predictors (pp. 312--313; 485 in the article, which uses
  # the normal quantile in place of t, and 487 in DMAR 1.0.x and MBESS,
  # which took the expected squared multiple correlation of X_j on p
  # predictors). These are the plans of the first 1.1.0 builds (79a03795)
  # with their published step.
  pub <- .src_c(assurance = .9)
  expect_identical(pub$term, c("necessary_N", "assurance", "assurance_attained"))
  expect_equal(pub$value, c(181, .9, NA))
  expect_identical(pub, .src_c(assurance = .9, assurance_method = "analytic"))
  expect_equal(ss_aipe_src(width = .2, true_cor_XX = .src_R5, true_cor_YX = rep(.3, 5),
                           which_predictor = 1, noncentral = FALSE,
                           assurance = .8)$value[1], 486)
  # The squared multiple correlations path gives the same plans.
  rr <- c(R2 = drop(.src_r %*% solve(.src_R, .src_r)), R2_j = .09)
  beta <- solve(.src_R, .src_r)[1]
  expect_equal(ss_aipe_src(true_R2_Y_X = rr[["R2"]], true_R2_j_X_without_j = rr[["R2_j"]],
                           p = 2, true_std_reg_coef_j = beta, width = .3,
                           noncentral = FALSE)$value, 156)
  # noncentral must be TRUE or FALSE.
  for (bad in list(NA, "no", c(TRUE, FALSE), 0)) {
    e <- tryCatch(ss_aipe_src(width = .3, true_cor_XX = .src_R, true_cor_YX = .src_r,
                              which_predictor = 1, noncentral = bad),
                  error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, "noncentral")
  }
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_aipe_rc() returns a tidy data.frame with necessary_N", {
  res <- ss_aipe_rc(true_R2_Y_X = 0.5, true_R2_j_X_without_j = 0.3,
                    p = 3, true_reg_coef_j = 0.4, width = 0.20)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_N" %in% res$term)
  expect_gt(res$value[res$term == "necessary_N"], 0)
})

test_that("ss_aipe_rc() requires more N for tighter widths", {
  wide  <- ss_aipe_rc(0.5, 0.3, p = 3, true_reg_coef_j = 0.4, width = 0.30)$value[1]
  tight <- ss_aipe_rc(0.5, 0.3, p = 3, true_reg_coef_j = 0.4, width = 0.10)$value[1]
  expect_gt(tight, wide)
})

test_that("ss_aipe_src() returns a tidy data.frame with necessary_N", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  res <- ss_aipe_src(true_R2_Y_X = 0.5, true_R2_j_X_without_j = 0.3,
                     p = 3, true_std_reg_coef_j = 0.4, width = 0.20)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_N" %in% res$term)
  expect_gt(res$value[res$term == "necessary_N"], 0)
})

test_that("ss_aipe_src() requires more N for tighter widths", {
  wide  <- ss_aipe_src(0.5, 0.3, p = 3, true_std_reg_coef_j = 0.4, width = 0.30)$value[1]
  tight <- ss_aipe_src(0.5, 0.3, p = 3, true_std_reg_coef_j = 0.4, width = 0.10)$value[1]
  expect_gt(tight, wide)
})

test_that("ss_aipe_rc() and ss_aipe_src() report the assurance rows of their plans", {
  R <- matrix(c(1, .3, .3, 1), 2)
  r <- c(.4, .3)
  rc <- ss_aipe_rc(width = .3, true_cor_XX = R, true_cor_YX = r,
                   which_predictor = 1, assurance = .9)
  expect_identical(rc$term, c("necessary_N", "assurance", "assurance_attained"))
  expect_equal(rc$value[1:2], c(187, .9))
  # The central plan is that of ss_aipe_reg_coef(); the exact probability
  # needs no simulation, so there is no G row and no standard error.
  expect_identical(rc, ss_aipe_reg_coef(
    width = .3, true_cor_XX = R, true_cor_YX = r, which_predictor = 1,
    assurance = .9))
  expect_null(attr(rc, "assurance_mc_se"))
  expect_identical(attr(rc, "assurance_method"), "analytic")
  expect_identical(attr(rc, "assurance_probability"), "exact")
  expect_identical(ss_aipe_rc(width = .3, true_cor_XX = R, true_cor_YX = r,
                              which_predictor = 1, assurance = .9,
                              assurance_method = "analytic"), rc)
  # The approximations of ss_aipe_src() report NA, the probability they do
  # not compute.
  for (nc in c(TRUE, FALSE)) {
    src <- ss_aipe_src(width = .3, true_cor_XX = R, true_cor_YX = r,
                       which_predictor = 1, assurance = .9, noncentral = nc)
    expect_identical(src$term, c("necessary_N", "assurance", "assurance_attained"))
    expect_true(is.na(src$value[3]))
    expect_identical(attr(src, "assurance_probability"), "not computed")
  }
})

test_that("ss_aipe_src() reports the rows of its simulated noncentral step", {
  # The one simulated plan on CRAN is the smoke plan of the central step
  # below; the two steps share .reg_coef_std_assurance_n().
  skip_on_cran()
  R <- matrix(c(1, .3, .3, 1), 2)
  r <- c(.4, .3)
  src <- ss_aipe_src(width = .3, true_cor_XX = R, true_cor_YX = r,
                     which_predictor = 1, assurance = .9,
                     assurance_method = "simulation",
                     G = 200, seed = 21)
  expect_identical(src$term, c("necessary_N", "assurance", "assurance_attained", "G"))
  expect_equal(src$value[4], 200)
  a <- src$value[3]
  expect_gte(a, .9)
  expect_equal(attr(src, "assurance_mc_se"), sqrt(a * (1 - a) / 200))
})

test_that("ss_aipe_src() reproduces its simulated plan with a seed, and it is that of ss_aipe_reg_coef(noncentral = TRUE)", {
  skip_on_cran()
  R <- matrix(c(1, .3, .3, 1), 2)
  r <- c(.4, .3)
  src <- ss_aipe_src(width = .3, true_cor_XX = R, true_cor_YX = r,
                     which_predictor = 1, assurance = .9,
                     assurance_method = "simulation",
                     G = 1000, seed = 21)
  expect_identical(src, ss_aipe_src(
    width = .3, true_cor_XX = R, true_cor_YX = r, which_predictor = 1,
    assurance = .9, assurance_method = "simulation", G = 1000, seed = 21))
  expect_identical(src, ss_aipe_reg_coef(
    width = .3, true_cor_XX = R, true_cor_YX = r, which_predictor = 1,
    noncentral = TRUE, assurance = .9, assurance_method = "simulation",
    G = 1000, seed = 21))
})

# ss_aipe_src(noncentral = FALSE): the central interval for the standardized
# coefficient (Kelley & Maxwell, 2003, Equation 1) ---------------------------

test_that("ss_aipe_src(noncentral = FALSE) with assurance runs one small seeded smoke plan", {
  # A small seeded run of the simulated assurance step (G = 200); the
  # block below checks the plans at the default G against raw data.
  plan <- .src_c(assurance = .9, assurance_method = "simulation", G = 200, seed = 4)
  expect_identical(plan$term, c("necessary_N", "assurance", "assurance_attained", "G"))
  expect_gte(plan$value[3], .9)
  expect_gte(plan$value[1], 156)
  expect_equal(attr(plan, "assurance_mc_se"),
               sqrt(plan$value[3] * (1 - plan$value[3]) / 200))
})

test_that("ss_aipe_src(noncentral = FALSE) attains its assurance, pinned to raw multivariate normal data", {
  skip_on_cran()
  # The plans, with seed 1 and the default G = 20000, and the share of
  # intervals no wider than the target among samples of raw multivariate
  # normal data standardized and fit by least squares without DMAR
  # (work/round6/evidence/g1_final_fixes/03_src_central_oracle.R):
  #   the probe, .3 and .90: 171, .9136 at 171 and .8951 at 170 (50,000
  #     samples each, Monte Carlo standard errors .0013 and .0014); the F
  #     law of ss_aipe_reg_coef() gives 187 (.997), the chi square 181
  #     (.986);
  #   Kelley and Maxwell (2003, pp. 312--313), .2 and .80: 478, .8069 at
  #     478 and .7969 at 477 (100,000 each, .0012 and .0013); the chi square
  #     gives 486 (.877).
  probe <- .src_c(assurance = .9, assurance_method = "simulation", seed = 1)
  expect_equal(probe$value[c(1, 2, 4)], c(171, .9, 20000))
  km <- ss_aipe_src(width = .2, true_cor_XX = .src_R5, true_cor_YX = rep(.3, 5),
                    which_predictor = 1, noncentral = FALSE, assurance = .8,
                    assurance_method = "simulation", seed = 1)
  expect_equal(km$value[1], 478)
  # The plan's own estimate agrees with the raw-data share at the plan.
  for (x in list(list(plan = probe, raw = .9136, se = .0013),
                 list(plan = km, raw = .8069, se = .0012))) {
    est <- x$plan$value[3]
    mc_se <- attr(x$plan, "assurance_mc_se")
    expect_lt(abs(est - x$raw), 4 * sqrt(mc_se^2 + x$se^2))
  }
  # A raw-data check in the test itself, at the probe's plan and one less.
  raw_share <- function(N, G, seed) {
    set.seed(seed)
    Rc <- chol(rbind(c(1, .src_r), cbind(.src_r, .src_R)))
    w <- vapply(seq_len(G), function(g) {
      D <- scale(matrix(rnorm(N * 3), N) %*% Rc)
      X <- cbind(1, D[, -1])
      XtXi <- solve(crossprod(X))
      b <- XtXi %*% crossprod(X, D[, 1])
      2 * qt(.975, N - 3) * sqrt(sum((D[, 1] - X %*% b)^2) / (N - 3) * XtXi[2, 2])
    }, 0)
    c(usable = sum(is.finite(w)), share = mean(w <= .3))
  }
  at_N <- raw_share(171, 4000, 171)
  expect_equal(at_N[["usable"]], 4000)  # attempted 4000, failed 0
  expect_lt(abs(at_N[["share"]] - probe$value[3]),
            4 * sqrt(probe$value[3] * (1 - probe$value[3]) / 4000 +
                       attr(probe, "assurance_mc_se")^2))
})

test_that("the central step's samples are the same at every N, its share rises with N, and its plan is the smallest N it allows", {
  skip_on_cran()
  sets <- list(
    list(R = .src_R, r = .src_r, j = 1, width = .3, a = .9),
    list(R = .src_R5, r = rep(.3, 5), j = 2, width = .25, a = .8),
    list(R = rbind(c(1, .4, .6), c(.4, 1, .05), c(.6, .05, 1)), r = c(.5, .3, .1),
         j = 3, width = .3, a = .95),
    list(R = matrix(1), r = .4, j = 1, width = .2, a = .6))
  for (s in sets) {
    p <- ncol(s$R)
    beta <- solve(s$R, s$r)[s$j]
    R2 <- sum(s$r * solve(s$R, s$r))
    R2_j <- if (p > 1) 1 - 1 / solve(s$R)[s$j, s$j] else 0
    ad <- DMAR:::.reg_coef_std_ad_values(R2, R2_j, beta, p, s$r[s$j])
    # A negative coefficient is planned as its magnitude (reverse the sign
    # of Y), as the step does.
    if (beta < 0) {
      beta <- -beta
      ad <- -ad
    }
    N_0 <- ss_aipe_src(width = s$width, true_cor_XX = s$R, true_cor_YX = s$r,
                       which_predictor = s$j, noncentral = FALSE)$value
    for (seed in 1:2) {
      plan <- ss_aipe_src(width = s$width, true_cor_XX = s$R, true_cor_YX = s$r,
                          which_predictor = s$j, noncentral = FALSE,
                          assurance = s$a, assurance_method = "simulation",
                          G = 500, seed = seed)
      N <- plan$value[1]
      set.seed(seed)
      base <- DMAR:::.reg_coef_std_base_draws(500, p)
      share <- function(n) min(DMAR:::.reg_coef_std_share(
        DMAR:::.reg_coef_std_draws_at(base, n, beta, R2, R2_j), ad, s$width,
        .025, .025, interval = "central"))
      expect_identical(share(N), plan$value[3])
      expect_gte(plan$value[3], s$a)
      if (N > N_0) expect_lt(share(N - 1), s$a)
      expect_true(all(diff(vapply(seq(N_0, N + 15), share, 0)) >= 0))
    }
  }
  # The central width of the draws is the t quantiles times the standard
  # error of the standardized coefficient formed from the two sample
  # squared multiple correlations, the interval of ci_reg_coef() with
  # noncentral = FALSE (see the noncentral counterpart in
  # test-ss_aipe_reg_coef.R for the standard error itself).
  set.seed(3)
  d <- DMAR:::.reg_coef_std_draws(5, 50, 2, .3, .2, .09)
  one <- function(t) rep(1, length(t))
  se <- DMAR:::.reg_coef_std_width(d, .1, one)
  delta <- DMAR:::.reg_coef_central_delta(47, .01, .04)
  expect_equal(DMAR:::.reg_coef_std_width(d, .1, delta),
               (qt(.99, 47) + qt(.96, 47)) * se)
  for (g in 1:5) {
    ci <- ci_reg_coef(reg_coef_j = .3, se_reg_coef_j = se[g], sd_Y = 1, sd_X = 1,
                      N = 50, p = 2, conf_level = NULL, alpha_lower = .01,
                      alpha_upper = .04)
    expect_equal(diff(ci$value[ci$term %in% c("ci_lower", "ci_upper")]),
                 DMAR:::.reg_coef_std_width(d, .1, delta)[g])
  }
})

test_that("ss_aipe_src(noncentral = FALSE) gives one plan for a coefficient and its negative, and refuses what it cannot plan", {
  skip_on_cran()
  plan <- function(r, ...) ss_aipe_src(width = .3, true_cor_XX = .src_R, true_cor_YX = r,
                                       which_predictor = 1, noncentral = FALSE,
                                       assurance = .85, assurance_method = "simulation",
                                       G = 1000, seed = 6, ...)
  expect_identical(plan(-.src_r), plan(.src_r))
  # With unequal tails the central width is the same whichever tail is which.
  expect_identical(plan(-.src_r, conf_level = NULL, alpha_lower = .01, alpha_upper = .04),
                   plan(.src_r, conf_level = NULL, alpha_lower = .04, alpha_upper = .01))
  # A seed leaves the caller's random number stream as it was.
  set.seed(9); before <- .Random.seed
  invisible(plan(.src_r))
  expect_identical(.Random.seed, before)
  # A limit at infinity is never narrow enough.
  e <- tryCatch(.src_c(conf_level = NULL, alpha_lower = 0, alpha_upper = .05,
                       assurance = .9), error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, c("alpha_lower", "alpha_upper"))
  # Squared multiple correlations no population allows with the coefficient,
  # refused in this function's names.
  e <- tryCatch(ss_aipe_src(true_R2_Y_X = .1, true_R2_j_X_without_j = .2, p = 3,
                            true_std_reg_coef_j = .5, width = .3, noncentral = FALSE,
                            assurance = .9, assurance_method = "simulation",
                            G = 200), error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_true("true_std_reg_coef_j" %in% e$arg)
  expect_false(grepl("(^|[^_])true_reg_coef_j", conditionMessage(e)))
  # G is checked when the central step simulates.
  expect_error(.src_c(assurance = .9, assurance_method = "simulation", G = 50),
               class = "dmar_count_invalid")
  # The squared multiple correlations path plans for the least favorable
  # correlation of Y with the predictor, at least the covariance path's plan
  # up to the Monte Carlo error of the two searches.
  R3 <- rbind(c(1, .4, .6), c(.4, 1, .05), c(.6, .05, 1)); r3 <- c(.5, .3, .1)
  beta <- solve(R3, r3)[1]
  r2_path <- ss_aipe_src(true_R2_Y_X = sum(r3 * solve(R3, r3)),
                         true_R2_j_X_without_j = 1 - 1 / solve(R3)[1, 1], p = 3,
                         true_std_reg_coef_j = beta, width = .3, noncentral = FALSE,
                         assurance = .9, assurance_method = "simulation",
                         G = 5000, seed = 2)
  cov_path <- ss_aipe_src(width = .3, true_cor_XX = R3, true_cor_YX = r3,
                          which_predictor = 1, noncentral = FALSE, assurance = .9,
                          assurance_method = "simulation", G = 5000, seed = 2)
  expect_gte(r2_path$value[1], cov_path$value[1] - 2)
})

test_that("ss_aipe_rc() with unequal tails plans with both quantiles, without a warning", {
  skip_on_cran()
  R <- matrix(c(1, .3, .3, 1), 2)
  r <- c(.4, .3)
  plan <- function(...) ss_aipe_rc(
    width = .3, true_cor_XX = R, true_cor_YX = r, which_predictor = 1,
    conf_level = NULL, alpha_lower = .04, alpha_upper = .01, ...)
  # The exact plan, 201 (the chi square step with both quantiles is 194),
  # and 169 without assurance.
  expect_no_warning(expect_identical(plan(assurance = .9)$value[1], 201))
  expect_no_warning(expect_identical(plan()$value, 169))
  # Width .4, tails .2 and .1, and assurance .6 give 31, where the chi
  # square step with the tails given is 32 (test-ss_aipe_reg_coef.R
  # derives both).
  wide <- function(...) ss_aipe_rc(width = .4, true_cor_XX = R, true_cor_YX = r,
                                   which_predictor = 1, conf_level = NULL, alpha_lower = .2,
                                   alpha_upper = .1, assurance = .6, ...)
  expect_no_warning(expect_identical(wide()$value[1], 31))
})
