# The correlation and regression functions name a population value with
# true_ and an English quantity (true_cor, true_cor_XX, true_reg_coef_j,
# true_std_reg_coef_j, true_R2_Y_X, true_mean_Y), a planning value of a
# sensitivity function with supposed_ (supposed_cor, supposed_R2,
# supposed_cov_XX), an observed coefficient and its standard error with
# reg_coef_j, std_reg_coef_j, and se_ (reg_coef_j, se_reg_coef_j,
# sd_Y), and an equivalence bound with the quantity it bounds (cor_lower,
# cor_upper). The names of DMAR 1.0.x and of the first builds of 1.1.0
# are not accepted: each stops with R's "unused argument" error. The calls
# are built with do.call() so that an old name is data here and not a call
# that a later rename sweep would change.

.old_names_refused <- function(fn, olds) {
  f <- get(fn, envir = asNamespace("DMAR"))
  expect_false("..." %in% names(formals(f)), info = fn)
  for (old in olds) {
    expect_error(do.call(f, stats::setNames(list(1), old)),
                 "unused argument", info = paste(fn, old))
  }
}

test_that("the correlation functions refuse rho and the other earlier names", {
  .old_names_refused("var_r", "rho")
  .old_names_refused("expected_r", "rho")
  .old_names_refused("expected_partial_r", "rho")
  .old_names_refused("equivalence_r", c("rho_lower", "rho_upper"))
  .old_names_refused("ss_aipe_r", c("rho", "true_rho"))
  .old_names_refused("ss_aipe_r_sensitivity", c("true_rho", "estimated_rho"))
  .old_names_refused("ss_aipe_partial_r", c("rho", "true_rho"))
  .old_names_refused("ss_aipe_partial_r_sensitivity",
                     c("true_rho", "estimated_rho"))
  .old_names_refused("ss_aipe_semipartial_r_sensitivity",
                     c("estimated_r_sp", "estimated_R2_full"))
  .old_names_refused("ss_aipe_equivalence_r", c("population_r", "true_rho"))
  .old_names_refused("ss_aipe_equivalence_r_sensitivity",
                     c("true_r", "estimated_r", "true_rho", "estimated_rho",
                       "rho_lower", "rho_upper"))
  .old_names_refused("ss_power_r", c("rho", "rho_0", "true_rho", "null_rho"))
  .old_names_refused("ss_aipe_two_phase", "true_rho")
  .old_names_refused("ss_aipe_two_phase_cost", "true_rho")
  .old_names_refused("ss_aipe_two_phase_sensitivity",
                     c("true_rho", "estimated_rho", "estimated_sd_Y"))
})

test_that("the regression functions refuse b_j, beta_j, rho2_, and rho_ names", {
  ci_old <- c("b_j", "SE_b_j", "s_Y", "s_X")
  .old_names_refused("ci_reg_coef", ci_old)
  .old_names_refused("ci_rc", ci_old)
  .old_names_refused("ci_src", c("beta_j", "SE_beta_j", ci_old))
  .old_names_refused("simulate_regression_data",
                     c("rho_YX", "rho_XX", "mu_Y", "mu_X"))
  later <- c("true_rho2_Y_X", "true_rho2_j_X_without_j", "true_rho_XX",
             "true_rho_YX")
  .old_names_refused("ss_aipe_reg_coef",
                     c("rho2_Y_X", "rho2_j_X_without_j", "b_j", "rho_XX",
                       "rho_YX", "true_b_j", later))
  .old_names_refused("ss_aipe_rc",
                     c("rho2_Y_X", "Rho2_j_X_without_j", "b_j", "rho_XX",
                       "rho_YX", "true_b_j", later))
  .old_names_refused("ss_aipe_src",
                     c("rho2_Y_X", "Rho2_j_X_without_j", "beta_j", "rho_XX",
                       "rho_YX", "true_beta_j", later))
  .old_names_refused("ss_power_reg_coef",
                     c("rho2_Y_X", "rho2_Y_X_without_j", "rho2_j_X_without_j",
                       "beta_j", "rho_XX", "rho_YX", "true_b_j",
                       "true_rho2_Y_X_without_j", later))
  .old_names_refused("ss_power_rc",
                     c("rho2_Y_X", "rho2_Y_X_without_j", "Rho2_j_X_without_j",
                       "beta_j", "rho_XX", "rho_YX", "true_b_j",
                       "true_rho2_Y_X_without_j", later))
  for (fn in c("ss_aipe_reg_coef_sensitivity", "ss_aipe_rc_sensitivity",
               "ss_aipe_src_sensitivity", "ss_power_reg_coef_sensitivity")) {
    .old_names_refused(fn, c("estimated_var_Y", "estimated_cov_YX",
                             "estimated_cov_XX"))
  }
  for (fn in c("ss_aipe_R2_sensitivity", "ss_power_R2_sensitivity")) {
    .old_names_refused(fn, c("estimated_R2", "rho_yx", "rho_xx",
                             "true_rho_YX", "true_rho_XX"))
  }
})

test_that("each function takes the final names", {
  has <- function(fn, new) {
    expect_true(all(new %in% names(formals(get(fn, envir = asNamespace("DMAR"))))),
                info = fn)
  }
  for (fn in c("var_r", "expected_r", "expected_partial_r", "ss_aipe_r",
               "ss_aipe_partial_r", "ss_aipe_equivalence_r",
               "ss_aipe_two_phase", "ss_aipe_two_phase_cost")) {
    has(fn, "true_cor")
  }
  for (fn in c("ss_aipe_r_sensitivity", "ss_aipe_partial_r_sensitivity",
               "ss_aipe_equivalence_r_sensitivity",
               "ss_aipe_two_phase_sensitivity")) {
    has(fn, c("true_cor", "supposed_cor"))
  }
  has("ss_aipe_two_phase_sensitivity", "supposed_sd_Y")
  has("ss_power_r", c("true_cor", "null_cor"))
  has("equivalence_r", c("cor_lower", "cor_upper"))
  has("ss_aipe_equivalence_r_sensitivity", c("cor_lower", "cor_upper"))
  has("ss_aipe_semipartial_r_sensitivity",
      c("supposed_r_sp", "supposed_R2_full"))
  has("ci_reg_coef", c("reg_coef_j", "se_reg_coef_j", "sd_Y", "sd_X"))
  has("ci_rc", c("reg_coef_j", "se_reg_coef_j", "sd_Y", "sd_X"))
  has("ci_src", c("std_reg_coef_j", "se_std_reg_coef_j", "reg_coef_j",
                  "se_reg_coef_j", "sd_Y", "sd_X"))
  has("simulate_regression_data",
      c("true_cor_YX", "true_cor_XX", "true_mean_Y", "true_mean_X"))
  for (fn in c("ss_aipe_reg_coef", "ss_aipe_rc", "ss_power_reg_coef",
               "ss_power_rc")) {
    has(fn, c("true_R2_Y_X", "true_R2_j_X_without_j", "true_reg_coef_j",
              "true_cor_XX", "true_cor_YX"))
  }
  has("ss_aipe_src", c("true_R2_Y_X", "true_R2_j_X_without_j",
                       "true_std_reg_coef_j", "true_cor_XX", "true_cor_YX"))
  has("ss_power_reg_coef", "true_R2_Y_X_without_j")
  has("ss_power_rc", "true_R2_Y_X_without_j")
  for (fn in c("ss_aipe_reg_coef_sensitivity", "ss_aipe_rc_sensitivity",
               "ss_aipe_src_sensitivity", "ss_power_reg_coef_sensitivity")) {
    has(fn, c("supposed_var_Y", "supposed_cov_YX", "supposed_cov_XX"))
  }
  for (fn in c("ss_aipe_R2_sensitivity", "ss_power_R2_sensitivity")) {
    has(fn, c("supposed_R2", "true_cor_YX", "true_cor_XX"))
  }
})

test_that("the echo rows carry the final names", {
  expect_true("true_cor" %in% ss_aipe_r(true_cor = .3, width = .2)$term)
  pr <- ss_power_r(true_cor = .4, null_cor = .2, desired_power = .8)
  expect_true(all(c("true_cor", "null_cor") %in% pr$term))
  expect_false(any(c("true_rho", "null_rho", "rho", "rho_0") %in% pr$term))
  eq <- equivalence_r(r = .05, N = 200, cor_upper = .1)
  expect_true(all(c("cor_lower", "cor_upper") %in% eq$term))
  expect_equal(eq$value[eq$term == "cor_lower"], -.1)
  expect_named(expected_r(true_cor = .5, N = 10),
               c("true_cor", "N", "expected_r", "bias", "relative_bias"))
  expect_true("true_cor" %in% names(expected_partial_r(.4, 30, 3)))
  tp <- ss_aipe_two_phase(true_sd_Y = .43, true_cor = .55, width = .1)
  expect_true("true_cor" %in% tp$term)
})

test_that("the regression sensitivity rows are named for the coefficient", {
  skip_on_cran()
  cov_XX <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX <- c(0.4, 0.3)
  pw <- ss_power_reg_coef_sensitivity(
    true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = cov_XX,
    which_predictor = 1, desired_power = 0.80, G = 20, print_iter = FALSE,
    seed = 113)
  rows <- c("mean_reg_coef_j", "median_reg_coef_j", "sd_reg_coef_j",
            "mean_se_reg_coef_j", "median_se_reg_coef_j", "sd_se_reg_coef_j",
            "true_reg_coef_j", "supposed_reg_coef_j")
  expect_true(all(rows %in% pw$term))
  expect_false(any(grepl("b_j$", pw$term)))
  expect_false(any(grepl("^estimated_", pw$term)))
  # The simulated t statistics are summarized under the t_value name.
  expect_true(all(c("mean_t_value", "median_t_value", "sd_t_value") %in%
                    pw$term))
  expect_false(any(c("mean_t", "median_t", "sd_t") %in% pw$term))
  r2 <- ss_power_R2_sensitivity(true_R2 = 0.3, specified_N = 40, p = 2,
                                G = 20, print_iter = FALSE, seed = 113)
  expect_true(all(c("mean_F_value", "median_F_value", "sd_F_value") %in%
                    r2$term))
  expect_false(any(c("mean_F", "median_F", "sd_F") %in% r2$term))
})

test_that("the per-replication files name each column for what it holds", {
  skip_on_cran()
  cov_XX <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  cov_YX <- c(0.4, 0.3)
  columns_written <- function(fn, ...) {
    path <- tempfile(fileext = ".csv")
    on.exit(unlink(path))
    suppressMessages(fn(true_var_Y = 1, true_cov_YX = cov_YX,
                        true_cov_XX = cov_XX, specified_N = 40,
                        which_predictor = 1, G = 3, print_iter = FALSE,
                        filename = path, seed = 113, ...))
    names(utils::read.csv(path))
  }
  aipe_raw <- c("reg_coef_j", "ci_lower", "ci_upper", "R2", "se_reg_coef_j",
                "t_value")
  aipe_std <- c("std_reg_coef_j", "ci_lower", "ci_upper", "R2",
                "se_std_reg_coef_j", "t_value")
  expect_identical(columns_written(ss_aipe_rc_sensitivity, width = 0.4),
                   aipe_raw)
  expect_identical(columns_written(ss_aipe_reg_coef_sensitivity, width = 0.4),
                   aipe_raw)
  # ss_aipe_src_sensitivity() reports the standardized coefficient; its
  # file and its rows both name it.
  expect_identical(columns_written(ss_aipe_src_sensitivity, width = 0.4),
                   aipe_std)
  expect_identical(columns_written(ss_power_reg_coef_sensitivity),
                   c("reg_coef_j", "se_reg_coef_j", "t_value", "R2"))
})

test_that("the two-phase rows are named in English", {
  skip_on_cran()
  cost <- ss_aipe_two_phase_cost(true_sd_Y = .43, true_cor = .55,
                                 width = .1, cost_ratio = .05)
  expect_true("threshold_cor" %in% cost$term)
  expect_false("threshold_rho" %in% cost$term)
  sens <- ss_aipe_two_phase_sensitivity(true_sd_Y = 1, true_cor = .7,
                                        supposed_sd_Y = 1, supposed_cor = .7,
                                        width = .3, n_pool = 2000, G = 20,
                                        seed = 4)
  expect_true(all(c("mean_mean_Y", "median_mean_Y", "sd_mean_Y",
                    "true_mean_Y") %in% sens$term))
  expect_false(any(grepl("theta", sens$term)))
})
