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
  expect_equal(v(res, "z"), 2.75, tolerance = 0.05)   # paper: z = 2.75
  expect_lt(v(res, "p_value") / 2, 0.005)             # one-tailed p = .003
  expect_equal(v(res, "k"), 18)
})

test_that("meta_contrast() matches first principles with preset zero-sum weights", {
  yi <- c(0.5, 0.1, -0.2); vi <- c(0.04, 0.02, 0.05)
  lam <- c(1, 0, -1)
  res <- meta_contrast(yi, vi, weights = lam)
  expect_equal(v(res, "estimate"), sum(lam * yi))
  expect_equal(v(res, "se"), sqrt(sum(lam^2 * vi)))
  expect_equal(v(res, "z"), sum(lam * yi) / sqrt(sum(lam^2 * vi)))
  expect_equal(v(res, "p_value"), 2 * pnorm(-abs(v(res, "z"))))
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
