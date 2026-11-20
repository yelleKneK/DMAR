v <- function(tab, t) tab$value[tab$term == t]

test_that("meta_contrast() reproduces the Raudenbush (1984) prior-contact contrast", {
  s  <- teacher_expectancy[-c(4, 5), ]
  d  <- append(s$d, .52, after = 3)
  wk <- append(s$weeks, 0, after = 3)
  ne <- append(s$n_experimental, 22, after = 3)
  nc <- append(s$n_control, 22, after = 3)
  vi <- (ne + nc) / (ne * nc) + d^2 / (2 * (ne + nc))
  # Weights inversely proportional to weeks of prior contact (constant 2,
  # the same constant as his reciprocal transformation), mean-centered.
  expect_message(res <- meta_contrast(d, vi, weights = 1 / (wk + 2)),
                 "mean-centered")
  expect_equal(v(res, "test_statistic"), 2.75, tolerance = 0.05)   # paper: z = 2.75
  expect_lt(v(res, "p_value") / 2, 0.005)             # one-tailed p = .003
  expect_equal(v(res, "n_studies"), 18)
})

test_that("meta_contrast() matches first principles with preset zero-sum weights", {
  yi <- c(0.5, 0.1, -0.2); vi <- c(0.04, 0.02, 0.05)
  lam <- c(1, 0, -1)
  res <- meta_contrast(yi, vi, weights = lam)
  expect_equal(v(res, "estimate"), sum(lam * yi))
  expect_equal(v(res, "se"), sqrt(sum(lam^2 * vi)))
  expect_equal(v(res, "test_statistic"), sum(lam * yi) / sqrt(sum(lam^2 * vi)))
  expect_equal(v(res, "p_value"), 2 * pnorm(-abs(v(res, "test_statistic"))))
})

test_that("the statistic is the row test_statistic, referred to the standard normal", {
  # Ruling 6 of 2026-11-11: one row name, test_statistic, for the
  # statistic of meta_es(), meta_r(), meta_smd(), and meta_contrast(),
  # with the reference distribution in attributes. DMAR 1.0.x named this
  # row z (and the count k); a lookup by the old name now finds nothing.
  yi <- c(0.5, 0.1, -0.2); vi <- c(0.04, 0.02, 0.05)
  res <- meta_contrast(yi, vi, weights = c(1, 0, -1))
  expect_identical(class(res), c("dmar_meta_contrast", "dmar_tbl", "data.frame"))
  expect_identical(res$term,
                   c("estimate", "se", "test_statistic", "p_value", "n_studies"))
  expect_false(any(c("z", "k") %in% res$term))
  expect_length(res$value[res$term == "z"], 0L)
  expect_identical(attr(res, "test_distribution"), "normal")
  expect_identical(attr(res, "test_df"), NA_real_)
  # The reference is the one meta_es() records without the Hartung-Knapp
  # adjustment, and the printed line is the same.
  fe <- meta_es(yi, vi, method = "fe")
  expect_identical(attr(res, "test_distribution"), attr(fe, "test_distribution"))
  expect_output(print(res), "test_statistic: referred to the standard normal",
                fixed = TRUE)
  # The two-sided p-value is that of the standard normal at the statistic
  # 0.7 / 0.3.
  expect_equal(v(res, "test_statistic"), 7 / 3, tolerance = 1e-14)
  expect_equal(v(res, "p_value"), 2 * pnorm(-7 / 3), tolerance = 1e-14)
})

test_that("meta_contrast() validates its arguments", {
  yi <- c(0.5, 0.1, -0.2); vi <- c(0.04, 0.02, 0.05)
  expect_error(meta_contrast(yi, c(0.1, -0.1, 0.2), weights = c(1, 0, -1)),
               "positive sampling variance")
  expect_error(meta_contrast(yi, vi, weights = c(1, 1)), "for each effect")
  expect_error(meta_contrast(yi, vi, weights = c(1, 1, 1), center = FALSE),
               "sum to zero")
  expect_error(suppressMessages(
    meta_contrast(yi, vi, weights = c(2, 2, 2))), "weights are zero")
})
