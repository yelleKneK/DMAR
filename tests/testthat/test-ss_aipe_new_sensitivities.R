## Tests for the 10 new AIPE sensitivity functions added in this commit.
## All Monte Carlo tests behind skip_on_cran() to preserve CRAN's runtime
## budget; the planners they pair with are covered on the CRAN-side check
## by the closed-form planner tests.

# ---------- ss_aipe_c_sensitivity ----------

test_that("ss_aipe_c_sensitivity() returns the tidy schema and respects monotonicity in true_var_error", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_c_sensitivity(
    true_var_error = 4, supposed_var_error = 4,
    c_weights = c(-1, 0, 1), width = 1,
    G = 60, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_type(res$value, "double")
  expect_true(all(c("mean_contrast", "mean_ci_width", "pct_ci_less_width",
                    "n_per_group") %in% res$term))

  # Higher true error variance should inflate the mean realized CI width.
  set.seed(113)
  res_hi <- ss_aipe_c_sensitivity(
    true_var_error = 16, supposed_var_error = 4,
    c_weights = c(-1, 0, 1), width = 1,
    G = 60, print_iter = FALSE
  )
  w_lo <- res$value[res$term == "mean_ci_width"]
  w_hi <- res_hi$value[res_hi$term == "mean_ci_width"]
  expect_gt(w_hi, w_lo)
})

# ---------- ss_aipe_omega_squared_sensitivity ----------

test_that("ss_aipe_omega_squared_sensitivity() returns the tidy schema and recovers the planning value", {
  skip_on_cran()
  set.seed(113)
  # suppressMessages: the planner's iterative search summarizes its
  # noncentral F lower-limit clamps in a message that is not under test.
  res <- suppressMessages(ss_aipe_omega_squared_sensitivity(
    true_omega_squared = 0.10, supposed_omega_squared = 0.10,
    df_effect = 2, width = 0.15, G = 60, print_iter = FALSE
  ))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true(all(c("mean_omega_squared", "mean_ci_width", "total_N") %in% res$term))
  # The mean realized estimate should be roughly near the truth at moderate G.
  expect_lt(abs(res$value[res$term == "mean_omega_squared"] - 0.10), 0.10)
})

# ---------- ss_aipe_partial_r_sensitivity ----------

test_that("ss_aipe_partial_r_sensitivity() recovers the planning rho on average", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_partial_r_sensitivity(
    true_cor = 0.40, supposed_cor = 0.40, p = 4, width = 0.30,
    G = 60, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean_partial_r", "mean_ci_width", "total_N") %in% res$term))
  expect_lt(abs(res$value[res$term == "mean_partial_r"] - 0.40), 0.10)
})

# ---------- ss_aipe_semipartial_r_sensitivity ----------

test_that("ss_aipe_semipartial_r_sensitivity() recovers the planning value on average", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_semipartial_r_sensitivity(
    true_r_sp = 0.30, supposed_r_sp = 0.30, p = 4, width = 0.30,
    G = 60, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean_r_sp", "mean_ci_width") %in% res$term))
  expect_lt(abs(res$value[res$term == "mean_r_sp"] - 0.30), 0.10)
})

# ---------- Both functions simulate p - 1 variables partialed out ----------
# In the tests below J is the number of variables partialed out, J = p - 1.

test_that("ss_aipe_partial_r_sensitivity() partials p - 1 variables and uses the Fisher's Z interval for them", {
  skip_if_not_installed("MASS")
  n <- 30; J <- 2; G <- 3; rho <- 0.3
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  # fisher_z = TRUE: the interval this function formed for every call
  # before it took the argument, so the widths below are those it gave.
  set.seed(7)
  res <- ss_aipe_partial_r_sensitivity(true_cor = rho, specified_N = n, p = J + 1,
                                       width = 0.2, G = G, filename = f,
                                       fisher_z = TRUE)
  expect_identical(attr(res, "interval"), "fisher_z")
  expect_output(print(res), "Interval: Fisher's Z scale", fixed = TRUE)
  per_rep <- utils::read.csv(f)

  # The same draws, made here: Y, X, and two variables to partial out.
  Sigma <- diag(J + 2)
  Sigma[1, 2] <- Sigma[2, 1] <- rho
  set.seed(7)
  own <- vapply(seq_len(G), function(g) {
    dat <- MASS::mvrnorm(n, rep(0, J + 2), Sigma)
    Z <- dat[, 3:4]
    stats::cor(stats::resid(stats::lm(dat[, 1] ~ Z)),
               stats::resid(stats::lm(dat[, 2] ~ Z)))
  }, numeric(1))
  expect_equal(per_rep$partial_r, own, tolerance = 1e-10)

  z <- atanh(per_rep$partial_r)
  half <- stats::qnorm(0.975) / sqrt(n - J - 3)
  expect_equal(per_rep$ci_width, tanh(z + half) - tanh(z - half),
               tolerance = 1e-10)
  expect_equal(res$value[res$term == "p"], J + 1)
  expect_false("J" %in% res$term)
  expect_equal(res$value[res$term == "total_N"], n)
})

test_that("ss_aipe_partial_r_sensitivity() simulates the raw-scale interval it plans by default", {
  skip_if_not_installed("MASS")
  n <- 30; J <- 2; G <- 4L; rho <- 0.3
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  res <- ss_aipe_partial_r_sensitivity(true_cor = rho, specified_N = n, p = J + 1,
                                       width = 0.2, G = G, filename = f,
                                       seed = 7)
  per_rep <- utils::read.csv(f)
  expect_identical(nrow(per_rep), G)
  # The same sample partial correlations as the Fisher run with this seed,
  # with the interval r +/- z (1 - r^2) / sqrt(n - J - 1), not truncated.
  half <- stats::qnorm(0.975) * (1 - per_rep$partial_r^2) / sqrt(n - J - 1)
  expect_equal(per_rep$ci_width, 2 * half, tolerance = 1e-10)
  expect_equal(per_rep$ci_lower, per_rep$partial_r - half, tolerance = 1e-10)
  expect_equal(per_rep$ci_upper, per_rep$partial_r + half, tolerance = 1e-10)
  expect_equal(res$value[res$term == "mean_ci_width"], mean(2 * half),
               tolerance = 1e-10)
  expect_identical(attr(res, "interval"), "raw")
  expect_output(print(res), "Interval: raw scale", fixed = TRUE)
  expect_s3_class(res, "dmar_tbl")
  # The raw-scale interval exists from N = p + 2; Fisher's from p + 3.
  expect_no_error(
    ss_aipe_partial_r_sensitivity(true_cor = 0.3, specified_N = 5, p = 3,
                                  width = 0.2, G = 2, seed = 1))
  expect_error(
    ss_aipe_partial_r_sensitivity(true_cor = 0.3, specified_N = 4, p = 3,
                                  width = 0.2, G = 2),
    "specified_N")
  expect_error(
    ss_aipe_partial_r_sensitivity(true_cor = 0.3, specified_N = 30, p = 3,
                                  width = 0.2, G = 2, fisher_z = NA),
    "'fisher_z' must be TRUE or FALSE")
})

test_that("ss_aipe_partial_r_sensitivity() plans the interval it simulates", {
  skip_if_not_installed("MASS")
  for (fz in c(FALSE, TRUE)) {
    res <- ss_aipe_partial_r_sensitivity(true_cor = 0.5, supposed_cor = 0.5,
                                         p = 3, width = 0.3, G = 2,
                                         fisher_z = fz, seed = 3)
    expect_equal(res$value[res$term == "total_N"],
                 ss_aipe_partial_r(true_cor = 0.5, p = 3, width = 0.3,
                                   fisher_z = fz)$value[1])
  }
  # With assurance the plan is the planner's, and the echo is appended.
  res <- ss_aipe_partial_r_sensitivity(true_cor = 0.5, supposed_cor = 0.5,
                                       p = 3, width = 0.3, G = 2,
                                       assurance = 0.9, seed = 3)
  expect_equal(res$value[res$term == "total_N"],
               ss_aipe_partial_r(true_cor = 0.5, p = 3, width = 0.3,
                                 assurance = 0.9)$value[1])
  expect_equal(res$value[res$term == "assurance"], 0.9)
})

test_that("the sensitivity functions of the lane take a seed and leave the caller's stream alone", {
  skip_if_not_installed("MASS")
  calls <- list(
    partial = function(seed) ss_aipe_partial_r_sensitivity(
      true_cor = 0.4, specified_N = 30, p = 3, width = 0.3, G = 4,
      seed = seed),
    semipartial = function(seed) ss_aipe_semipartial_r_sensitivity(
      true_r_sp = 0.4, specified_N = 30, p = 3, width = 0.3, G = 4,
      seed = seed),
    semipartial_assurance = function(seed) ss_aipe_semipartial_r_sensitivity(
      true_r_sp = 0.4, supposed_r_sp = 0.4, p = 2, width = 0.4, G = 4,
      assurance = 0.8, seed = seed),
    cv = function(seed) ss_aipe_cv_sensitivity(
      true_cv = 0.2, specified_N = 20, width = 0.2, G = 4, seed = seed),
    R2 = function(seed) suppressWarnings(ss_aipe_R2_sensitivity(
      true_R2 = 0.4, specified_N = 60, width = 0.3, p = 2, G = 4,
      print_iter = FALSE, seed = seed))
  )
  for (nm in names(calls)) {
    fn <- calls[[nm]]
    set.seed(2026)
    before <- .Random.seed
    a <- fn(11)
    expect_identical(.Random.seed, before, info = nm)
    b <- fn(11)
    expect_identical(a, b, info = nm)
    expect_false(identical(a$value, fn(12)$value), info = nm)
    # Without a seed the draws come from, and advance, the caller's stream.
    set.seed(2026)
    c1 <- fn(NULL)
    expect_false(identical(.Random.seed, before), info = nm)
    set.seed(2026)
    expect_identical(fn(NULL), c1, info = nm)
  }
  # seed is the last formal of each.
  for (nm in c("ss_aipe_partial_r_sensitivity",
               "ss_aipe_semipartial_r_sensitivity", "ss_aipe_cv_sensitivity",
               "ss_aipe_R2_sensitivity", "ss_aipe_semipartial_r")) {
    current <- names(formals(get(nm, envir = asNamespace("DMAR"))))
    expect_identical(current[length(current)], "seed", info = nm)
    expect_null(formals(get(nm, envir = asNamespace("DMAR")))$seed)
  }
  # fisher_z was appended after the released formals, before seed.
  expect_identical(
    utils::tail(names(formals(ss_aipe_partial_r_sensitivity)), 3L),
    c("filename", "fisher_z", "seed"))
})

test_that("ss_aipe_semipartial_r_sensitivity() partials p - 1 variables out of X and forms the planned Wald interval", {
  skip_if_not_installed("MASS")
  n <- 40; J <- 2; G <- 3; r_sp <- 0.4; R2 <- 0.5
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  set.seed(11)
  res <- ss_aipe_semipartial_r_sensitivity(true_r_sp = r_sp, specified_N = n,
                                           p = J + 1, width = 0.2, G = G,
                                           filename = f, true_R2_full = R2)
  per_rep <- utils::read.csv(f)

  Sigma <- diag(J + 2)
  Sigma[1, 2] <- Sigma[2, 1] <- r_sp
  Sigma[1, 3] <- Sigma[3, 1] <- sqrt(R2 - r_sp^2)
  set.seed(11)
  own <- vapply(seq_len(G), function(g) {
    dat <- MASS::mvrnorm(n, rep(0, J + 2), Sigma)
    Y <- dat[, 1]; X <- dat[, 2]; Z <- dat[, 3:4]
    r <- stats::cor(Y, stats::resid(stats::lm(X ~ Z)))
    R2_hat <- summary(stats::lm(Y ~ X + Z))$r.squared
    # The interval of var_semipartial_r(under_null = FALSE), from the
    # public function.
    v <- var_semipartial_r(r, n, p = J + 1, R2_full = R2_hat, under_null = FALSE)$value
    c(r = r, width = 2 * stats::qnorm(0.975) * sqrt(v))
  }, numeric(2))
  expect_equal(per_rep$r_sp, unname(own["r", ]), tolerance = 1e-10)
  expect_equal(per_rep$ci_width, unname(own["width", ]), tolerance = 1e-8)
  expect_equal(per_rep$ci_upper - per_rep$r_sp, per_rep$r_sp - per_rep$ci_lower,
               tolerance = 1e-10)

  # The population has the stated full-model R2: Sigma gives it directly.
  expect_equal(drop(Sigma[1, -1] %*% solve(Sigma[-1, -1]) %*% Sigma[-1, 1]), R2)
  expect_equal(res$value[res$term == "true_R2_full"], R2)
  expect_true(is.na(res$value[res$term == "supposed_R2_full"]))
  expect_equal(res$term[seq(nrow(res) - 1, nrow(res))],
               c("true_R2_full", "supposed_R2_full"))
})

test_that("ss_aipe_semipartial_r_sensitivity() echoes the R2_full values and passes the planning value on", {
  skip_if_not_installed("MASS")
  # The plan with assurance is a simulation, so the two calls share a seed:
  # the sensitivity function sets it before it calls the planner.
  res <- ss_aipe_semipartial_r_sensitivity(
    true_r_sp = 0.8, supposed_r_sp = 0.8, p = 2, width = 0.10, G = 3,
    true_R2_full = 0.94, supposed_R2_full = 0.94, assurance = 0.80,
    seed = 113)
  plan <- ss_aipe_semipartial_r(true_r_sp = 0.8, p = 2, width = 0.10,
                                true_R2_full = 0.94, assurance = 0.80, seed = 113)
  expect_equal(res$value[res$term == "total_N"], plan$value[1])
  expect_equal(res$value[res$term == "supposed_R2_full"], 0.94)
  # The optional assurance echo stays last.
  expect_equal(res$term[seq(nrow(res) - 2, nrow(res))],
               c("true_R2_full", "supposed_R2_full", "assurance"))
  # Without true_R2_full the population has R2_full = true_r_sp^2.
  set.seed(113)
  res0 <- ss_aipe_semipartial_r_sensitivity(
    true_r_sp = 0.3, specified_N = 50, p = 3, width = 0.3, G = 3)
  expect_equal(res0$value[res0$term == "true_R2_full"], 0.09)
})

test_that("the partial and semipartial sensitivity functions check p, G, specified_N, and R2_full", {
  for (bad_p in list(2.5, 1, 0, NA_real_)) {
    expect_error(
      ss_aipe_partial_r_sensitivity(true_cor = 0.3, supposed_cor = 0.3,
                                    p = bad_p, width = 0.2, G = 2),
      "whole number of at least 2")
    expect_error(
      ss_aipe_semipartial_r_sensitivity(true_r_sp = 0.3, supposed_r_sp = 0.3,
                                        p = bad_p, width = 0.2, G = 2),
      "whole number of at least 2")
  }
  expect_error(
    ss_aipe_partial_r_sensitivity(true_cor = 0.3, supposed_cor = 0.3,
                                  width = 0.2, G = 2),
    "'p' (the number of predictors) is required.", fixed = TRUE)
  expect_error(
    ss_aipe_semipartial_r_sensitivity(true_r_sp = 0.3, supposed_r_sp = 0.3,
                                      width = 0.2, G = 2),
    "'p' (the number of predictors) is required.", fixed = TRUE)
  expect_error(
    ss_aipe_partial_r_sensitivity(true_cor = 0.3, specified_N = 40, p = 2,
                                  width = 0.2, G = 2.5),
    class = "dmar_count_invalid")
  expect_error(
    ss_aipe_semipartial_r_sensitivity(true_r_sp = 0.3, specified_N = 40,
                                      p = 2, width = 0.2, G = 0),
    "'G'")
  # The Fisher's Z interval needs N >= p + 3; a fractional size is refused.
  expect_error(
    ss_aipe_partial_r_sensitivity(true_cor = 0.3, specified_N = 5, p = 3,
                                  width = 0.2, G = 2, fisher_z = TRUE),
    "specified_N")
  expect_error(
    ss_aipe_semipartial_r_sensitivity(true_r_sp = 0.3, specified_N = 40.5,
                                      p = 3, width = 0.2, G = 2),
    "specified_N")

  # An R2_full below the squared semipartial correlation is refused by class,
  # for the population and for the planning value.
  expect_error(
    ss_aipe_semipartial_r_sensitivity(true_r_sp = 0.9, specified_N = 40,
                                      p = 2, width = 0.2, G = 2,
                                      true_R2_full = 0.1),
    class = "dmar_semipartial_r_incoherent")
  expect_error(
    ss_aipe_semipartial_r_sensitivity(true_r_sp = 0.9, supposed_r_sp = 0.9,
                                      p = 2, width = 0.2, G = 2,
                                      supposed_R2_full = 0.1),
    class = "dmar_refusal")
  expect_error(
    ss_aipe_semipartial_r_sensitivity(true_r_sp = 0.6, specified_N = 40,
                                      p = 2, width = 0.2, G = 2,
                                      true_R2_full = 1),
    "below 1")
  expect_error(
    ss_aipe_semipartial_r_sensitivity(true_r_sp = 0.3, specified_N = 40,
                                      p = 2, width = 0.2, G = 2,
                                      supposed_R2_full = 0.5),
    "specified_N")
})

test_that("the partial and semipartial sensitivity functions give no warning", {
  skip_if_not_installed("MASS")
  set.seed(113)
  expect_no_warning(
    ss_aipe_partial_r_sensitivity(true_cor = 0.3, specified_N = 30, p = 2,
                                  width = 0.2, G = 2))
  expect_no_warning(
    ss_aipe_semipartial_r_sensitivity(true_r_sp = 0.3, specified_N = 30,
                                      p = 2, width = 0.2, G = 2))
})

test_that("the corrected intervals of the two sensitivity functions cover at the nominal rate", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  get <- function(res, term) res$value[res$term == term]

  # Partial: five variables partialed out, N = 30 (Monte Carlo SE of a
  # 5 percent rate at G = 4000 is .0034).
  set.seed(20261109)
  res <- ss_aipe_partial_r_sensitivity(true_cor = 0.4, specified_N = 30,
                                       p = 6, width = 0.2, G = 4000,
                                       fisher_z = TRUE)
  expect_lt(abs(get(res, "total_type_I_error") - 0.05), 0.012)

  # Semipartial, controls that explain much of the outcome, above and
  # below r_sp^2 = 1/2. DMAR 1.0.x's two interval forms missed these
  # populations at rates near .13 and .001.
  set.seed(20261109)
  res <- ss_aipe_semipartial_r_sensitivity(
    true_r_sp = 0.8, supposed_r_sp = 0.8, p = 2, width = 0.10, G = 4000,
    true_R2_full = 0.94, supposed_R2_full = 0.94)
  expect_equal(get(res, "total_N"), 331)
  expect_lt(abs(get(res, "total_type_I_error") - 0.05), 0.015)
  expect_lt(abs(get(res, "mean_ci_width") - 0.10), 0.003)

  res <- ss_aipe_semipartial_r_sensitivity(
    true_r_sp = 0, supposed_r_sp = 0, p = 2, width = 0.10, G = 4000,
    true_R2_full = 0.64, supposed_R2_full = 0.64)
  expect_lt(abs(get(res, "total_type_I_error") - 0.05), 0.015)
  expect_lt(abs(get(res, "mean_ci_width") - 0.10), 0.003)
})

# ---------- ss_aipe_indirect_effect_sensitivity ----------

test_that("ss_aipe_indirect_effect_sensitivity() returns a tidy summary; closed_form and monte_carlo both work", {
  skip_on_cran()
  set.seed(113)
  res_s <- ss_aipe_indirect_effect_sensitivity(
    true_std_reg_coef_a = 0.4, true_std_reg_coef_b = 0.3, supposed_std_reg_coef_a = 0.4, supposed_std_reg_coef_b = 0.3,
    width = 0.30, method = "closed_form", G = 50, print_iter = FALSE
  )
  expect_s3_class(res_s, "data.frame")
  expect_true(all(c("mean_std_indirect_effect", "mean_ci_width", "true_std_indirect_effect") %in% res_s$term))
  expect_equal(res_s$value[res_s$term == "true_std_indirect_effect"], 0.12, tolerance = 1e-10)

  set.seed(113)
  res_m <- ss_aipe_indirect_effect_sensitivity(
    true_std_reg_coef_a = 0.4, true_std_reg_coef_b = 0.3, supposed_std_reg_coef_a = 0.4, supposed_std_reg_coef_b = 0.3,
    width = 0.30, method = "monte_carlo", B = 400, G = 30, print_iter = FALSE
  )
  expect_s3_class(res_m, "data.frame")
})

# ---------- ss_aipe_equivalence_smd_sensitivity ----------

test_that("ss_aipe_equivalence_smd_sensitivity() reports pct_equivalent inside the bounds", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_equivalence_smd_sensitivity(
    true_smd = 0.0, supposed_smd = 0.0, width = 0.40,
    smd_upper = 0.5,
    G = 60, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean_smd", "pct_equivalent",
                    "smd_lower", "smd_upper") %in% res$term))
  expect_equal(res$value[res$term == "smd_lower"], -0.5)
  expect_equal(res$value[res$term == "smd_upper"],  0.5)
})

# ---------- ss_aipe_cliff_delta_sensitivity ----------

test_that("ss_aipe_cliff_delta_sensitivity() recovers Cliff's delta on average", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_cliff_delta_sensitivity(
    true_delta = 0.30, supposed_delta = 0.30, width = 0.40,
    G = 60, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean_cliff_delta", "mean_ci_width", "n_1", "n_2") %in% res$term))
  expect_lt(abs(res$value[res$term == "mean_cliff_delta"] - 0.30), 0.15)
})

# ---------- ss_aipe_pcm_sensitivity ----------

test_that("ss_aipe_pcm_sensitivity() returns the tidy schema", {
  skip_on_cran()
  set.seed(113)
  res <- ss_aipe_pcm_sensitivity(
    true_var_trend = 0.003, true_var_error = 0.0262,
    supposed_var_trend = 0.003, supposed_var_error = 0.0262,
    duration = 4, frequency = 1, width = 0.025,
    G = 30, print_iter = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_type(res$value, "double")
  expect_true(all(c("mean_slope_diff", "mean_ci_width", "n_per_group",
                    "n_occasions") %in% res$term))
  expect_equal(res$value[res$term == "n_occasions"], 5)  # f*D + 1 = 1*4 + 1
})

test_that("ss_aipe_pcm_sensitivity() evaluates the planner's two-group estimand", {
  skip_on_cran()
  # ss_aipe_pcm() sizes n PER GROUP so the expected width of the confidence
  # interval on the BETWEEN-group difference in change coefficients is at most
  # `width` (a two-sample pooled-variance t interval on per-subject estimated
  # slopes, with 2n - 2 degrees of freedom). The sensitivity simulator must
  # evaluate that same two-group estimand, so it has to (a) resolve the per-
  # group n the planner returns and (b) realize a mean CI width at or just
  # under the target with coverage near conf_level. A regression to a one-group
  # quantity (a single group, df = n - 1, standard error sd / sqrt(n)) would
  # land near half the target width and fail the width band below, which is the
  # load-bearing assertion here.

  # frequency = 1 (the Kelley & Rausch tables): f^(2p) factor is 1.
  n1 <- ss_aipe_pcm(true_var_trend = 0.003, true_var_error = 0.0262,
                    duration = 4, frequency = 1, width = 0.025)$value
  set.seed(113)
  r1 <- ss_aipe_pcm_sensitivity(
    true_var_trend = 0.003, true_var_error = 0.0262,
    supposed_var_trend = 0.003, supposed_var_error = 0.0262,
    duration = 4, frequency = 1, width = 0.025, G = 100, print_iter = FALSE
  )
  v1 <- function(t) r1$value[r1$term == t]
  expect_equal(v1("n_per_group"), n1)        # simulator uses the planned n
  expect_lt(abs(v1("mean_ci_width") - 0.025), 0.20 * 0.025)  # width on target
  expect_lt(v1("total_type_I_error"), 0.15)      # coverage ~ conf_level

  # frequency = 2 (f^(2p) = 4 at p = 1): the planner's V, and therefore n, must
  # scale with frequency^(2p). If that factor were dropped from the planner the
  # returned n would be too small and the simulator (which fits OLS on the true
  # time grid) would realize an interval wider than the target, failing the
  # width band. So this case guards the f^(2p) factor end to end. The pure
  # frequency effect is isolated by comparing the planner at f = 1 and f = 2 at
  # a common width: more frequent sampling inflates the per-unit-time slope
  # variance V by f^(2p), so the f = 2 design needs the larger n.
  n2 <- ss_aipe_pcm(true_var_trend = 0.003, true_var_error = 0.0262,
                    duration = 2, frequency = 2, width = 0.05)$value
  n1_w05 <- ss_aipe_pcm(true_var_trend = 0.003, true_var_error = 0.0262,
                        duration = 4, frequency = 1, width = 0.05)$value
  expect_gt(n2, n1_w05)                          # V scales by f^(2p) = 4
  set.seed(113)
  r2 <- ss_aipe_pcm_sensitivity(
    true_var_trend = 0.003, true_var_error = 0.0262,
    supposed_var_trend = 0.003, supposed_var_error = 0.0262,
    duration = 2, frequency = 2, width = 0.05, G = 100, print_iter = FALSE
  )
  v2 <- function(t) r2$value[r2$term == t]
  expect_equal(v2("n_per_group"), n2)
  expect_lt(abs(v2("mean_ci_width") - 0.05), 0.20 * 0.05)
  expect_lt(v2("total_type_I_error"), 0.15)
})

# ---------- ss_aipe_mixed_effects_sensitivity ----------

test_that("ss_aipe_mixed_effects_sensitivity() returns the tidy schema", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  set.seed(113)
  # suppressMessages: some replications land on lme4's boundary (singular)
  # fit, which lme4 reports as a message.
  res <- suppressMessages(suppressWarnings(ss_aipe_mixed_effects_sensitivity(
    true_total_var = 1, true_var_X = 1, true_icc = 0.10, true_reg_coef = 0.30,
    supposed_total_var = 1, supposed_var_X = 1, supposed_icc = 0.10,
    width = 0.30, cluster_size = 20,
    G = 25, print_iter = FALSE
  )))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean_reg_coef", "mean_ci_width", "n_clusters",
                    "cluster_size") %in% res$term))
})

# ---------- ss_aipe_reliability_sensitivity ----------

test_that("ss_aipe_reliability_sensitivity() recovers the planning reliability for alpha", {
  skip_on_cran()
  set.seed(113)
  res <- suppressWarnings(ss_aipe_reliability_sensitivity(
    true_reliability = 0.80, supposed_reliability = 0.80, n_items = 6,
    width = 0.15, estimator = "alpha", G = 30, print_iter = FALSE
  ))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mean_reliability", "mean_ci_width", "total_N",
                    "n_items") %in% res$term))
  expect_lt(abs(res$value[res$term == "mean_reliability"] - 0.80), 0.10)
})

# ---------- Input-validation sanity check (sample) ----------

test_that("the new sensitivity functions reject conflicting estimated_* and specified_N", {
  expect_error(
    ss_aipe_c_sensitivity(true_var_error = 4, supposed_var_error = 4,
                          specified_n_per_group = 50, c_weights = c(-1, 0, 1), width = 1, G = 5),
    "but not both"
  )
  expect_error(
    ss_aipe_omega_squared_sensitivity(true_omega_squared = 0.1,
                                      supposed_omega_squared = 0.1,
                                      specified_N = 50, df_effect = 2, width = 0.1, G = 5),
    "but not both"
  )
  expect_error(
    ss_aipe_partial_r_sensitivity(true_cor = 0.4, supposed_cor = 0.4,
                                  specified_N = 30, p = 4, width = 0.2, G = 5),
    "but not both"
  )
})

test_that("the partial and semipartial sensitivity functions check specified_N and true_cor", {
  # With p = 3 the interval needs N - p > 0 and the fit N >= p + 2.
  expect_error(
    ss_aipe_semipartial_r_sensitivity(true_r_sp = .3, specified_N = 4, p = 3,
                                      width = .2, G = 2),
    "specified_N")
  # A vector true_cor stopped with R's "'length = 2' in coercion to
  # 'logical(1)'".
  expect_error(
    ss_aipe_partial_r_sensitivity(true_cor = c(.4, .5), specified_N = 40,
                                  p = 3, width = .3, G = 50),
    "'true_cor' must be a single value")
})

test_that("the partial and semipartial sensitivity functions take p in its place and do not accept J", {
  skip_if_not_installed("MASS")
  # The count holds the third place and is the number of predictors.
  expect_identical(
    ss_aipe_partial_r_sensitivity(.3, NULL, 3, .3, specified_N = 30, G = 3,
                                  seed = 5),
    ss_aipe_partial_r_sensitivity(true_cor = .3, specified_N = 30, p = 3,
                                  width = .3, G = 3, seed = 5))
  expect_identical(
    ss_aipe_semipartial_r_sensitivity(.3, NULL, 3, .3, specified_N = 30,
                                      G = 3, seed = 5),
    ss_aipe_semipartial_r_sensitivity(true_r_sp = .3, specified_N = 30,
                                      p = 3, width = .3, G = 3, seed = 5))
  expect_error(
    ss_aipe_partial_r_sensitivity(true_cor = .3, specified_N = 30, J = 2,
                                  width = .3, G = 3), "unused argument")
  expect_error(
    ss_aipe_semipartial_r_sensitivity(true_r_sp = .3, specified_N = 30,
                                      J = 2, width = .3, G = 3),
    "unused argument")
  # The planned sample size message and the smallest specified_N are stated
  # in p: p + 2 on the raw scale and for the semipartial, p + 3 on Fisher's.
  expect_no_error(
    ss_aipe_semipartial_r_sensitivity(true_r_sp = .3, specified_N = 6, p = 4,
                                      width = .3, G = 2, seed = 1))
  expect_error(
    ss_aipe_semipartial_r_sensitivity(true_r_sp = .3, specified_N = 5, p = 4,
                                      width = .3, G = 2),
    "specified_N")
  expect_no_error(
    ss_aipe_partial_r_sensitivity(true_cor = .3, specified_N = 7, p = 4,
                                  width = .3, G = 2, fisher_z = TRUE,
                                  seed = 1))
  expect_error(
    ss_aipe_partial_r_sensitivity(true_cor = .3, specified_N = 6, p = 4,
                                  width = .3, G = 2, fisher_z = TRUE),
    "specified_N")
})
