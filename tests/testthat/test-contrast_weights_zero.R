# Contrast weights with no nonzero value define no contrast. Every exported
# function that takes contrast weights refuses them with a condition of
# class dmar_argument_invalid whose `arg` names the argument as that
# function spells it (.dmar_check_contrast_nonzero(), R/dmar_conditions.R).
# Such weights pass the rules that the weights sum to zero and that the
# positive ones sum to at most 1. Before the check, ci_c() and the other
# interval functions returned the interval [0, 0], ss_aipe_c() a plan of 2
# per group, ss_power_c() a power of 1 with an infinite noncentrality
# parameter, and ss_aipe_sc_ancova_sensitivity() with a specified sample
# size stopped with "object 'k' not found"; the functions that refused
# them already did so with a message about the positive weights.

# One expectation: a refusal of class dmar_argument_invalid, with `arg`
# and a message matching `regexp`.
expect_no_contrast <- function(expr, arg = "c_weights",
                               regexp = "no nonzero weight") {
  e <- tryCatch({
    expr
    NULL
  }, error = function(e) e)
  got <- if (is.null(e)) "no error" else paste0(
    "class ", paste(class(e), collapse = "/"), ", arg ", deparse(e$arg),
    ", message: ", conditionMessage(e))
  expect(inherits(e, "dmar_argument_invalid") && identical(e$arg, arg) &&
           grepl(regexp, conditionMessage(e), fixed = TRUE),
         paste0("Expected a dmar_argument_invalid refusal with arg ",
                deparse(arg), " and a message with \"", regexp, "\"; got ",
                got, "."))
  invisible(e)
}

z <- c(0, 0)
one_way <- data.frame(
  g = factor(rep(c("a", "b", "c"), each = 4)),
  x = c(1.2, 0.4, 2.2, 1.5, 0.9, 1.8, 0.3, 1.1, 2.0, 1.4, 0.7, 1.6),
  y = c(3.1, 2.4, 4.0, 3.3, 4.2, 5.0, 3.6, 4.4, 2.9, 2.5, 2.0, 3.0))
fit_anova <- stats::lm(y ~ g, data = one_way)
fit_ancova <- stats::lm(y ~ g + x, data = one_way)

test_that("the confidence interval functions refuse all-zero weights", {
  expect_no_contrast(ci_c(means = c(1, 2), s_anova = 1, c_weights = z,
                          n_per_group = c(20, 20)))
  expect_no_contrast(ci_sc(means = c(1, 2), s_anova = 1, c_weights = z,
                           n_per_group = c(20, 20)))
  expect_no_contrast(ci_c_ancova(adj_means = c(1, 2), s_ancova = .8,
                                 c_weights = z, n_per_group = c(20, 20),
                                 cov_means = c(0, .1), SSwithin_x = 40))
  expect_no_contrast(ci_c_ancova(adj_means = c(1, 2), s_ancova = .8,
                                 c_weights = z, n_per_group = c(20, 20),
                                 cov_means = c(0, .1), SSwithin_x = 40,
                                 procedure = "bryant_paulson"))
  expect_no_contrast(ci_c_ancova_bp(adj_means = c(1, 2), s_ancova = .8,
                                    c_weights = z, n_per_group = 20))
  expect_no_contrast(ci_sc_ancova(adj_means = c(1, 2), s_anova = 1,
                                  s_ancova = .8, c_weights = z,
                                  n_per_group = c(20, 20),
                                  cov_means = c(0, .1), SSwithin_x = 40))
  expect_no_contrast(ci_scheffe(fit_anova, contrasts = c(0, 0, 0)),
                     arg = "contrasts")
  expect_no_contrast(contrast_test(fit_anova, contrasts = c(0, 0, 0)),
                     arg = "contrasts")
  expect_no_contrast(contrast_adjusted(fit_ancova, contrast = c(0, 0, 0)),
                     arg = "contrast")
})

test_that("the accuracy planners and their sensitivity functions refuse all-zero weights", {
  expect_no_contrast(ss_aipe_c(true_var_error = 1, c_weights = z,
                               width = .5))
  expect_no_contrast(ss_aipe_c_ancova(true_var_error_anova = 1,
                                      true_cor = .22, c_weights = z,
                                      width = .5))
  expect_no_contrast(ss_aipe_c_sensitivity(true_var_error = 1,
                                           supposed_var_error = 1,
                                           c_weights = z, width = .5,
                                           G = 10, seed = 1))
  expect_no_contrast(ss_aipe_c_ancova_sensitivity(
    true_var_error_ancova = 1, supposed_var_error_ancova = 1,
    true_cor = .2, true_mean_Y = c(10, 12), true_mean_X = 2,
    true_sd_X = 1.3, true_sd_Y = 2, c_weights = z, width = .5, G = 10,
    seed = 1))
  expect_no_contrast(ss_aipe_sc(true_std_contrast = .5, c_weights = z,
                                width = .5))
  expect_no_contrast(ss_aipe_sc_sensitivity(true_std_contrast = .5,
                                            supposed_std_contrast = .5,
                                            c_weights = z, width = .5,
                                            G = 10, print_iter = FALSE,
                                            seed = 1))
  expect_no_contrast(ss_aipe_sc_ancova(true_std_contrast = .5,
                                       true_cor = .5, c_weights = z,
                                       width = .5))
  # With a specified sample size no planner is called, so the function's
  # own checks of the weights, one for each divisor, refuse them.
  expect_no_contrast(ss_aipe_sc_ancova_sensitivity(
    true_std_contrast = .5, c_weights = z, specified_n_per_group = 20,
    width = .5, true_cor = .5, G = 10, print_iter = FALSE, seed = 1))
  expect_no_contrast(ss_aipe_sc_ancova_sensitivity(
    true_std_contrast = .5, c_weights = c(0, 0, 0),
    specified_n_per_group = 20, width = .5, true_cor = .5,
    divisor = "s_anova", G = 10, print_iter = FALSE, seed = 1))
  expect_no_contrast(ss_seq_c(c_weights = z, half_width = 2.5, pilot = TRUE))
  expect_no_contrast(ss_seq_c_sensitivity(c_weights = z, half_width = 2.5,
                                          true_sd = 15, G = 10, seed = 1))
})

test_that("the power and equivalence functions refuse all-zero weights", {
  expect_no_contrast(ss_power_c(true_contrast = .5, c_weights = z,
                                true_sd = 1))
  expect_no_contrast(ss_power_c_ancova(true_contrast = .5, c_weights = z,
                                       true_sd = 1, true_cor = .5))
  expect_no_contrast(ss_power_sc(true_std_contrast = .5, c_weights = z))
  expect_no_contrast(ss_power_contrast(c_weights = z, true_contrast = .5,
                                       true_var_error = 1))
  expect_no_contrast(ss_power_equivalence_c(c_weights = z, true_sd = 1,
                                            contrast_upper = .5))
  expect_no_contrast(power_equivalence_c(c_weights = z,
                                         n_per_group = c(20, 20),
                                         true_sd = 1, contrast_upper = .5))
  expect_no_contrast(equivalence_c(means = c(1, 2), s_anova = 1,
                                   c_weights = z, n_per_group = c(20, 20),
                                   contrast_upper = .5))
})

test_that("meta_contrast() and is_orthogonal_set() refuse all-zero weights", {
  expect_no_contrast(meta_contrast(yi = c(.2, .4, .6), vi = c(.04, .05, .06),
                                   weights = c(0, 0, 0)), arg = "weights")
  # Weights that are all equal are zero once centered (with a message
  # that they were centered); that refusal now names the argument too.
  expect_no_contrast(suppressMessages(
    meta_contrast(yi = c(.2, .4, .6), vi = c(.04, .05, .06),
                  weights = c(1, 1, 1))),
    arg = "weights", regexp = "zero after centering")
  expect_no_contrast(is_orthogonal_set(cbind(c(1, -1, 0), c(0, 0, 0))),
                     arg = "contrasts")
})

test_that("a set of contrasts with an all-zero one is refused and the message names it", {
  # Rows of ci_c_ancova_bp()'s c_weights and of contrast_test()'s
  # resolved contrasts are contrasts; columns of ci_scheffe()'s and
  # is_orthogonal_set()'s are.
  expect_no_contrast(ci_c_ancova_bp(adj_means = c(1, 2, 3), s_ancova = .8,
                                    c_weights = rbind(c(1, -1, 0), c(0, 0, 0)),
                                    n_per_group = 20),
                     regexp = "Row 2 of 'c_weights' has no nonzero weight")
  expect_no_contrast(ci_scheffe(fit_anova,
                                contrasts = cbind(c(1, -1, 0), c(0, 0, 0))),
                     arg = "contrasts",
                     regexp = "Column 2 of 'contrasts' has no nonzero weight")
  expect_no_contrast(contrast_test(fit_anova,
                                   contrasts = list(a_b = c(1, -1, 0),
                                                    none = c(0, 0, 0))),
                     arg = "contrasts",
                     regexp = "The contrast 'none' in 'contrasts' has no nonzero weight")
  expect_no_contrast(is_orthogonal_set(cbind(c(0, 0, 0), c(1, -1, 0),
                                             c(0, 0, 0))),
                     arg = "contrasts",
                     regexp = "Columns 1 and 3 of 'contrasts' have no nonzero weight")
})

test_that("empty weights are refused too", {
  # numeric(0) has no nonzero weight. On the contrast path, and on the
  # means path with empty means, nothing else checked its length: ci_c()
  # and ci_c_ancova() returned the limits NaN (or [.5, .5] with N given),
  # and the accuracy planners stopped with "missing value where TRUE/FALSE
  # needed".
  e <- numeric(0)
  expect_no_contrast(ci_c(contrast = .5, s_anova = 1, c_weights = e,
                          n_per_group = 20))
  expect_no_contrast(ci_c(contrast = .5, s_anova = 1, c_weights = e,
                          n_per_group = 20, N = 40))
  expect_no_contrast(ci_c(means = e, s_anova = 1, c_weights = e,
                          n_per_group = 20))
  expect_no_contrast(ci_c_ancova(contrast = .5, s_ancova = .8,
                                 c_weights = e, n_per_group = 20,
                                 cov_means = e, SSwithin_x = 40))
  expect_no_contrast(ss_aipe_c(true_var_error = 1, c_weights = e,
                               width = .5))
  expect_no_contrast(ss_aipe_c_ancova(true_var_error_anova = 1,
                                      true_cor = .22, c_weights = e,
                                      width = .5))
  expect_no_contrast(ss_aipe_c_sensitivity(true_var_error = 1,
                                           supposed_var_error = 1,
                                           c_weights = e, width = .5,
                                           G = 10, seed = 1))
  expect_no_contrast(ss_aipe_sc_ancova_sensitivity(
    true_std_contrast = .5, c_weights = e, specified_n_per_group = 20,
    width = .5, true_cor = .5, G = 10, print_iter = FALSE, seed = 1))
  # With empty means as well; with two means the length of the weights is
  # refused first, under arg "true_mean_Y".
  expect_no_contrast(ss_aipe_c_ancova_sensitivity(
    true_var_error_ancova = 1, supposed_var_error_ancova = 1,
    true_cor = .2, true_mean_Y = e, true_mean_X = 2, true_sd_X = 1.3,
    true_sd_Y = 2, c_weights = e, width = .5, G = 10, seed = 1))
  # A matrix or list with no contrast in it stopped with a message from R
  # ("arguments imply differing number of rows: 1, 0", "length of
  # 'dimnames' [2] not equal to array extent", "length of 'dimnames' [1]
  # not equal to array extent", and "argument must be coercible to
  # non-negative integer").
  expect_no_contrast(ci_c_ancova_bp(adj_means = c(1, 2), s_ancova = .8,
                                    c_weights = matrix(numeric(0), 0, 2),
                                    n_per_group = 20))
  expect_no_contrast(ci_scheffe(fit_anova,
                                contrasts = matrix(numeric(0), 3, 0)),
                     arg = "contrasts")
  expect_no_contrast(contrast_test(fit_anova,
                                   contrasts = matrix(numeric(0), 0, 3)),
                     arg = "contrasts")
  expect_no_contrast(contrast_test(fit_anova, contrasts = list()),
                     arg = "contrasts")
})

test_that("logical and complex weights with no nonzero value are refused", {
  # R computes with FALSE as 0, so c(FALSE, FALSE) passes the rule that the
  # weights sum to zero. ci_c() returned [0, 0], ss_aipe_c() a plan of 2,
  # and ci_scheffe() a row with a standard error of 0. ci_sc_ancova() and
  # ss_aipe_sc_ancova() refused them before their checks moved to the
  # shared helper, and must still do so (for complex zeros too).
  f2 <- c(FALSE, FALSE)
  expect_no_contrast(ci_c(means = c(1, 2), s_anova = 1, c_weights = f2,
                          n_per_group = 20))
  expect_no_contrast(ss_aipe_c(true_var_error = 1, c_weights = f2,
                               width = .5))
  expect_no_contrast(ci_scheffe(fit_anova,
                                contrasts = c(FALSE, FALSE, FALSE)),
                     arg = "contrasts")
  expect_no_contrast(ci_sc_ancova(adj_means = c(1, 2), s_anova = 1,
                                  s_ancova = .8, c_weights = f2,
                                  n_per_group = 20, cov_means = c(0, .1),
                                  SSwithin_x = 40))
  expect_no_contrast(ss_aipe_sc_ancova(true_std_contrast = .5,
                                       true_cor = .5, c_weights = f2,
                                       width = .5))
  expect_no_contrast(ss_aipe_sc_ancova(true_std_contrast = .5,
                                       true_cor = .5,
                                       c_weights = complex(2), width = .5))
  # Logical weights with a TRUE do not sum to zero and are refused, as
  # before, by that rule.
  expect_no_contrast(ci_c(means = c(1, 2), s_anova = 1,
                          c_weights = c(TRUE, FALSE), n_per_group = 20),
                     regexp = "should equal zero")
})

test_that("the rule that the positive weights sum to 1 names 'c_weights' in the equivalence functions", {
  # c(2, -2) sums to zero, so it reaches the rule; the refusal carried no
  # arg before.
  expect_no_contrast(ss_power_equivalence_c(c_weights = c(2, -2), true_sd = 1,
                                            contrast_upper = .5),
                     regexp = "The positive weights must sum to 1")
  expect_no_contrast(power_equivalence_c(c_weights = c(2, -2),
                                         n_per_group = c(20, 20),
                                         true_sd = 1, contrast_upper = .5),
                     regexp = "The positive weights must sum to 1")
  expect_no_contrast(equivalence_c(means = c(1, 1.1), s_anova = 1,
                                   c_weights = c(2, -2),
                                   n_per_group = c(20, 20),
                                   contrast_upper = .5),
                     regexp = "The positive weights must sum to 1")
})

test_that("nonzero weights are accepted as before", {
  # A spot check that the new rule leaves weights with a nonzero value
  # alone, including a zero weight among nonzero ones.
  expect_s3_class(ci_c(means = c(1, 2, 3), s_anova = 1,
                       c_weights = c(1, 0, -1), n_per_group = 20),
                  "data.frame")
  expect_s3_class(ss_aipe_c(true_var_error = 1, c_weights = c(1, 0, -1),
                            width = .5), "data.frame")
  expect_s3_class(ss_power_c(true_contrast = .5, c_weights = c(1, 0, -1),
                             true_sd = 1), "data.frame")
  expect_s3_class(ss_aipe_sc_ancova_sensitivity(
    true_std_contrast = .5, c_weights = c(1, -1),
    specified_n_per_group = 20, width = .5, true_cor = .5, G = 10,
    print_iter = FALSE, seed = 1), "data.frame")
})
