## The filename contract of the sensitivity functions that write their
## per-replication results directly with utils::write.csv() or
## utils::write.table(): nothing is written unless a path is supplied, a
## path that is not a single string is refused before any replication runs,
## and a supplied path receives one row per replication with the documented
## columns. The 2026-09-04 CRAN review retired the save argument and the
## default path these functions inherited from MBESS. Every call below
## bypasses the planner (specified_N or n_per_group) and uses a tiny G, so
## the file is the object under test, not the Monte Carlo precision.

Sigma_X <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
cov_YX  <- c(0.4, 0.3)

.direct_writer_calls <- function(G, path) {
  list(
    ss_aipe_R2_sensitivity = bquote(ss_aipe_R2_sensitivity(
      true_R2 = 0.3, specified_N = 40, w = 0.3, p = 2,
      G = .(G), print_iter = FALSE, filename = .(path))),
    ss_power_R2_sensitivity = bquote(ss_power_R2_sensitivity(
      true_R2 = 0.3, specified_N = 40, p = 2,
      G = .(G), print_iter = FALSE, filename = .(path))),
    ss_aipe_reg_coef_sensitivity = bquote(ss_aipe_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = .(cov_YX), true_cov_XX = .(Sigma_X),
      specified_N = 40, which_predictor = 1, w = 0.4,
      G = .(G), print_iter = FALSE, filename = .(path))),
    ss_power_reg_coef_sensitivity = bquote(ss_power_reg_coef_sensitivity(
      true_var_Y = 1, true_cov_YX = .(cov_YX), true_cov_XX = .(Sigma_X),
      specified_N = 40, which_predictor = 1,
      G = .(G), print_iter = FALSE, filename = .(path))),
    ss_aipe_sc_ancova_sensitivity = bquote(ss_aipe_sc_ancova_sensitivity(
      true_psi = 0.5, n_per_group = 10, c_weights = c(-1, 0, 1),
      desired_width = 0.5, rho = 0.4,
      G = .(G), print_iter = FALSE, filename = .(path))),
    ss_aipe_sc_ancova_sensitivity_s_anova = bquote(ss_aipe_sc_ancova_sensitivity(
      true_psi = 0.5, n_per_group = 10, c_weights = c(-1, 0, 1),
      desired_width = 0.5, rho = 0.4, divisor = "s_anova",
      G = .(G), print_iter = FALSE, filename = .(path)))
  )
}

.expected_columns <- list(
  ss_aipe_R2_sensitivity = c("lower_limit", "r2", "upper_limit",
                             "lower_width_ci", "upper_width_ci", "width_ci"),
  ss_power_R2_sensitivity = c("r2", "f_stat"),
  ss_aipe_reg_coef_sensitivity = c("b_j", "ll_ci_beta_j", "ul_ci_beta_j",
                                   "r_2", "se_b_j", "t_for_b_j"),
  ss_power_reg_coef_sensitivity = c("b_j", "se_b_j", "t_stat", "r_2"),
  ss_aipe_sc_ancova_sensitivity = c("psi_obs", "full_width", "width_lower",
                                    "width_upper", "type_I_error_upper",
                                    "type_I_error_lower", "type_I_error",
                                    "lower_limit", "upper_limit"),
  ss_aipe_sc_ancova_sensitivity_s_anova = c("psi_obs", "full_width", "width_lower",
                                            "width_upper", "type_I_error_upper",
                                            "type_I_error_lower", "type_I_error",
                                            "lower_limit", "upper_limit")
)

# The write.csv() members (the two R2 members and the power member for the
# regression coefficient) overwrite; the write.table() members append.
.appends <- c(ss_aipe_R2_sensitivity = FALSE, ss_power_R2_sensitivity = FALSE,
              ss_aipe_reg_coef_sensitivity = TRUE, ss_power_reg_coef_sensitivity = FALSE,
              ss_aipe_sc_ancova_sensitivity = TRUE,
              ss_aipe_sc_ancova_sensitivity_s_anova = TRUE)

test_that("a supplied filename receives one row per replication with the documented columns", {
  for (member in names(.expected_columns)) {
    path <- tempfile(fileext = ".csv")
    call <- .direct_writer_calls(G = 4, path = path)[[member]]
    set.seed(113)
    # The write.csv() members announce the path in a message; the
    # write.table() members print a line only when they append.
    res <- suppressMessages(eval(call))
    expect_s3_class(res, "data.frame")
    expect_true(file.exists(path), info = member)
    per_rep <- utils::read.csv(path)
    expect_named(per_rep, .expected_columns[[member]])
    expect_equal(nrow(per_rep), 4L, info = member)

    set.seed(113)
    out <- utils::capture.output(res2 <- suppressMessages(eval(call)))
    expect_equal(nrow(utils::read.csv(path)),
                 if (.appends[[member]]) 8L else 4L, info = member)
    if (.appends[[member]]) {
      expect_match(paste(out, collapse = "\n"), "already exists", info = member)
    }
    unlink(path)
    expect_false(file.exists(path))
  }
})

test_that("nothing is written by default", {
  before <- list.files(getwd(), all.files = TRUE)
  for (call in .direct_writer_calls(G = 2, path = NULL)) {
    set.seed(113)
    res <- suppressMessages(eval(call))
    expect_s3_class(res, "data.frame")
  }
  expect_identical(list.files(getwd(), all.files = TRUE), before)
})

test_that("a filename that is not a single string is refused before any replication runs", {
  fns <- list(ss_aipe_R2_sensitivity, ss_power_R2_sensitivity,
              ss_aipe_reg_coef_sensitivity, ss_power_reg_coef_sensitivity,
              ss_aipe_sc_ancova_sensitivity, ss_aipe_c_sensitivity)
  for (fn in fns) {
    expect_null(formals(fn)$filename)
    expect_false("save" %in% names(formals(fn)))
  }
  for (bad in list(TRUE, 1, c("a.csv", "b.csv"), NA_character_, "")) {
    calls <- .direct_writer_calls(G = 10000, path = bad)
    for (member in names(calls)) {
      expect_error(eval(calls[[member]]),
                   "'filename' must be NULL or a single string",
                   info = member)
    }
  }
})
