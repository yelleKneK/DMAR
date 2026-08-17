v <- function(tab, t) tab$value[tab$term == t]

# Study-level Raudenbush (1984) inputs: the 18 experiments, with the
# Pellegrini & Hicks study at its study-level values (d = .52, p = .010,
# total n = 44; see ?teacher_expectancy).
.r84 <- function() {
  s <- teacher_expectancy[-c(4, 5), ]
  list(p  = append(s$p_one_tailed, .010, after = 3),
       df = append(s$n_experimental + s$n_control - 2, 42, after = 3),
       d  = append(s$d, .52, after = 3),
       wk = append(s$weeks, 0, after = 3),
       ne = append(s$n_experimental, 22, after = 3),
       nc = append(s$n_control, 22, after = 3))
}

test_that("combine_p() reproduces Raudenbush (1984) Table 2", {
  r <- .r84()
  res <- combine_p(r$p, weights = r$df)
  # Fisher: chi-square(36) = 62.17; the paper prints p = .0025, though
  # 62.17 on 36 df implies p = .0043 (the recomputed side pinned below).
  expect_equal(v(res, "fisher_chi_square"), 62.17, tolerance = 0.01)
  expect_equal(v(res, "fisher_df"), 36)
  expect_lt(v(res, "fisher_p"), 0.005)
  # Edgington: the paper prints sum p = 6.84 with p = .04 (Table 2, p. 90);
  # the p-values as tabled sum to 7.00, just above the .05 level (p = .0509).
  # ?combine_p describes both, so pin the computed side from both
  # directions.
  expect_equal(v(res, "edgington_sum_p"), 7.00, tolerance = 0.01)
  expect_gt(v(res, "edgington_p"), 0.05)
  expect_lt(v(res, "edgington_p"), 0.06)
  # Mosteller-Bush adding Zs: z = 2.12 in the paper.
  expect_equal(v(res, "stouffer_z"), 2.12, tolerance = 0.15)
  expect_lt(v(res, "stouffer_p"), 0.05)
  # Weighted by degrees of freedom: z = 0.834, NOT significant; the large
  # studies found smaller effects.
  expect_equal(v(res, "stouffer_weighted_z"), 0.834, tolerance = 0.05)
  expect_gt(v(res, "stouffer_weighted_p"), 0.05)
  expect_equal(v(res, "k"), 18)
})

test_that("combine_p() reproduces the Raudenbush (1984) Table 3 subgroups", {
  r <- .r84()
  lo <- r$wk <= 2
  res_lo <- combine_p(r$p[lo], weights = r$df[lo])
  res_hi <- combine_p(r$p[!lo], weights = r$df[!lo])
  # Low prior contact (10 studies): all four tests significant.
  expect_equal(v(res_lo, "fisher_chi_square"), 54.18, tolerance = 0.05)
  expect_equal(v(res_lo, "stouffer_z"), 4.16, tolerance = 0.05)
  expect_equal(v(res_lo, "stouffer_weighted_z"), 2.52, tolerance = 0.02)
  # High prior contact (8 studies): nothing is significant.
  expect_equal(v(res_hi, "fisher_chi_square"), 7.98, tolerance = 0.02)
  expect_equal(v(res_hi, "edgington_sum_p"), 5.26, tolerance = 0.01)
  expect_lt(v(res_hi, "stouffer_z"), 0)
  expect_gt(v(res_hi, "stouffer_p"), 0.5)
})

test_that("combine_p() methods agree with first principles on small cases", {
  p <- c(0.02, 0.20, 0.40)
  res <- combine_p(p)
  expect_equal(v(res, "fisher_chi_square"), -2 * sum(log(p)))
  expect_equal(v(res, "fisher_p"),
               pchisq(-2 * sum(log(p)), df = 6, lower.tail = FALSE))
  expect_equal(v(res, "stouffer_z"), sum(qnorm(1 - p)) / sqrt(3))
  expect_equal(v(res, "edgington_sum_p"), sum(p))
  # A single-method request returns only that method's rows (plus k).
  fish <- combine_p(p, method = "fisher")
  expect_identical(fish$term, c("fisher_chi_square", "fisher_df",
                                "fisher_p", "k"))
})

test_that("combine_p() validates its arguments", {
  expect_error(combine_p(c(0.5, 1)), "inside \\(0, 1\\)")
  expect_error(combine_p(c(0, 0.5)), "inside \\(0, 1\\)")
  expect_error(combine_p(0.5), "two or more")
  expect_error(combine_p(c(.2, .3), method = "stouffer_weighted"),
               "needs 'weights'")
  expect_error(combine_p(c(.2, .3), weights = c(1, -1)), "non-negative")
  expect_error(combine_p(c(.2, .3), weights = 1), "one weight per")
})
