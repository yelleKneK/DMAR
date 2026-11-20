test_that("anova_within_two_way() returns documented rows", {
  set.seed(113)
  n_sub <- 12
  grid <- expand.grid(subject = factor(1:n_sub),
                      A = factor(c("a1", "a2", "a3")),
                      B = factor(c("b1", "b2")))
  grid$y <- rnorm(nrow(grid))
  res <- anova_within_two_way(grid, "y", "A", "B", "subject")
  expect_equal(nrow(res), 12)  # 3 effects x 4 adjustments
  expect_setequal(unique(res$effect), c("A", "B", "A:B"))
  expect_setequal(unique(res$adjustment),
                  c("none", "Greenhouse-Geisser", "Huynh-Feldt", "lower_bound"))
})

test_that("anova_within_two_way() eta_squared_partial is constant across adjustments", {
  set.seed(113)
  n_sub <- 8
  grid <- expand.grid(subject = factor(1:n_sub),
                      A = factor(c("a1", "a2", "a3")),
                      B = factor(c("b1", "b2", "b3")))
  grid$y <- rnorm(nrow(grid)) + 0.6 * as.integer(grid$A)
  res <- anova_within_two_way(grid, "y", "A", "B", "subject")
  eta_A <- res$eta_squared_partial[res$effect == "A"]
  expect_true(all(diff(eta_A) == 0))
})

test_that("anova_within_two_way() A:B F and GG epsilon match car::Anova", {
  skip_if_not_installed("car")
  # Fixed, balanced 3-by-2 within-subjects design (mirrors the one-way
  # sibling test that anchors against aov()/car::Anova).
  set.seed(113)
  n_sub <- 8
  grid <- expand.grid(subject = factor(1:n_sub),
                      A = factor(c("a1", "a2", "a3")),
                      B = factor(c("b1", "b2")))
  grid$y <- rnorm(nrow(grid)) + 0.6 * as.integer(grid$A) +
            0.3 * as.integer(grid$B)
  res <- anova_within_two_way(grid, "y", "A", "B", "subject")

  # Reference: car::Anova on the wide multivariate within-subjects fit.
  cells <- c("a1", "a2", "a3")
  wm <- matrix(NA_real_, n_sub, 6L)
  col <- 1L
  for (bb in c("b1", "b2")) for (aa in cells) {
    sel <- grid$A == aa & grid$B == bb
    wm[as.integer(grid$subject[sel]), col] <- grid$y[sel]
    col <- col + 1L
  }
  idata <- expand.grid(A = factor(cells), B = factor(c("b1", "b2")))
  av <- car::Anova(lm(wm ~ 1), idata = idata, idesign = ~ A * B,
                   type = "III")
  s <- summary(av, multivariate = FALSE)
  f_ref  <- s$univariate.tests["A:B", "F value"]
  gg_ref <- s$pval.adjustments["A:B", "GG eps"]

  ab <- res[res$effect == "A:B", ]
  expect_equal(ab$F_value[ab$adjustment == "none"], unname(f_ref),
               tolerance = 1e-6)
  expect_equal(ab$epsilon[ab$adjustment == "Greenhouse-Geisser"],
               unname(gg_ref), tolerance = 1e-6)
})

test_that("anova_within_two_way() errors on unbalanced design", {
  set.seed(113)
  d <- expand.grid(subject = factor(1:5),
                   A = factor(c("a1", "a2")),
                   B = factor(c("b1", "b2")))
  d$y <- rnorm(nrow(d))
  d <- d[-1, ]  # drop one row
  expect_error(anova_within_two_way(d, "y", "A", "B", "subject"),
               "not balanced")
})

test_that("epsilon rows are NA with a warning when n - 1 < df_effect", {
  # A 3 x 3 design with n = 4 subjects: the interaction has df = 4,
  # but the covariance matrix of its 4 orthonormal contrasts has rank
  # at most n - 1 = 3, so it is singular and the Greenhouse-Geisser
  # and Huynh-Feldt epsilons are artifacts bounded above by
  # (n - 1)/df = 3/4 regardless of the population epsilon. The main
  # effects (df = 2) are estimable and unaffected.
  set.seed(113)
  g <- expand.grid(subject = factor(1:4),
                   A = factor(paste0("a", 1:3)),
                   B = factor(paste0("b", 1:3)))
  g$y <- rnorm(nrow(g))
  expect_warning(res <- anova_within_two_way(g, "y", "A", "B", "subject"),
                 "cannot be estimated")

  # The output shape is unchanged: 3 effects x 4 adjustments.
  expect_equal(nrow(res), 12L)
  expect_setequal(unique(res$adjustment),
                  c("none", "Greenhouse-Geisser", "Huynh-Feldt",
                    "lower_bound"))

  ab <- res[res$effect == "A:B", ]
  for (adj in c("Greenhouse-Geisser", "Huynh-Feldt")) {
    row <- ab[ab$adjustment == adj, ]
    expect_true(is.na(row$epsilon))
    expect_true(is.na(row$df_1))
    expect_true(is.na(row$df_2))
    expect_true(is.na(row$p_value))
    # The observed F is never adjusted, so it stays reportable.
    expect_false(is.na(row$F_value))
    expect_false(is.na(row$eta_squared_partial))
  }

  # The unadjusted and lower-bound rows for A:B are intact; the lower
  # bound 1/df needs no estimate of the covariance matrix.
  un <- ab[ab$adjustment == "none", ]
  expect_false(anyNA(un[, c("F_value", "df_1", "df_2", "p_value")]))
  lb <- ab[ab$adjustment == "lower_bound", ]
  expect_equal(lb$epsilon, 1 / 4)
  expect_false(anyNA(lb[, c("F_value", "df_1", "df_2", "p_value",
                            "epsilon")]))

  # Main effects with df = 2 <= n - 1 = 3 still carry estimated
  # epsilons and adjusted tests.
  for (eff in c("A", "B")) {
    est <- res[res$effect == eff &
               res$adjustment %in% c("Greenhouse-Geisser", "Huynh-Feldt"), ]
    expect_false(anyNA(est[, c("epsilon", "df_1", "df_2", "p_value")]))
  }
})

test_that("exactly one warning, naming every affected effect", {
  # A 4 x 4 design with n = 3 subjects: every effect's df (3, 3, and 9)
  # exceeds n - 1 = 2, so all three effects are named in one warning.
  set.seed(113)
  g <- expand.grid(subject = factor(1:3),
                   A = factor(paste0("a", 1:4)),
                   B = factor(paste0("b", 1:4)))
  g$y <- rnorm(nrow(g))
  w <- capture_warnings(res <- anova_within_two_way(g, "y", "A", "B",
                                                    "subject"))
  expect_length(w, 1L)
  expect_match(w, "A \\(df = 3\\)")
  expect_match(w, "B \\(df = 3\\)")
  expect_match(w, "A:B \\(df = 9\\)")
  expect_match(w, "n = 3 subjects")
  gg <- res[res$adjustment == "Greenhouse-Geisser", ]
  hf <- res[res$adjustment == "Huynh-Feldt", ]
  expect_true(all(is.na(gg$epsilon)))
  expect_true(all(is.na(hf$epsilon)))
  lb <- res[res$adjustment == "lower_bound", ]
  expect_equal(lb$epsilon[match(c("A", "B", "A:B"), lb$effect)],
               c(1 / 3, 1 / 3, 1 / 9))
})

test_that("no warning and estimated epsilons when n - 1 >= df_effect", {
  # n = 5 subjects in a 3 x 2 design: the largest effect df is 2,
  # within n - 1 = 4, so every epsilon is estimable. Includes the
  # boundary n - 1 = df via the 3 x 3 design with n = 5 (df_AB = 4).
  set.seed(113)
  g <- expand.grid(subject = factor(1:5),
                   A = factor(paste0("a", 1:3)),
                   B = factor(paste0("b", 1:2)))
  g$y <- rnorm(nrow(g))
  expect_no_warning(res <- anova_within_two_way(g, "y", "A", "B",
                                                "subject"))
  est <- res[res$adjustment != "none", ]
  expect_false(anyNA(est[, c("epsilon", "df_1", "df_2", "p_value")]))

  g3 <- expand.grid(subject = factor(1:5),
                    A = factor(paste0("a", 1:3)),
                    B = factor(paste0("b", 1:3)))
  g3$y <- rnorm(nrow(g3))
  expect_no_warning(res3 <- anova_within_two_way(g3, "y", "A", "B",
                                                 "subject"))
  est3 <- res3[res3$adjustment != "none", ]
  expect_false(anyNA(est3[, c("epsilon", "df_1", "df_2", "p_value")]))
})

test_that("the lower bound epsilon is 1/df, the attainable infimum, per effect", {
  # The help page previously printed 1 / (df - 1); the code computes
  # 1 / df with df the effect's numerator degrees of freedom: 1/2 for A
  # (df = 2), 1/2 for B (df = 2), 1/4 for A:B (df = 4) in a 3 x 3
  # design. The old expression would give 1, 1, and 1/3 instead.
  set.seed(113)
  n_sub <- 10
  g3 <- expand.grid(subject = factor(1:n_sub),
                    A = factor(paste0("a", 1:3)),
                    B = factor(paste0("b", 1:3)))
  g3$y <- rnorm(nrow(g3)) + rep(rnorm(n_sub), times = 9)
  res <- anova_within_two_way(g3, "y", "A", "B", "subject")
  lb <- res[res$adjustment == "lower_bound", ]
  un <- res[res$adjustment == "none", ]
  df_e <- un$df_1[match(lb$effect, un$effect)]
  expect_equal(lb$epsilon, 1 / df_e)
  expect_equal(lb$epsilon[match(c("A", "B", "A:B"), lb$effect)],
               c(1 / 2, 1 / 2, 1 / 4))
})

# The design of the help example: a balanced 3 x 2 within-subjects design
# with 12 subjects.
within_example <- function() {
  set.seed(113)
  n_sub <- 12
  grid <- expand.grid(subject = factor(1:n_sub),
                      A = factor(c("a1", "a2", "a3")),
                      B = factor(c("b1", "b2")))
  grid$y <- with(grid,
    3 * (A == "a2") + 1 * (A == "a3") +
    2 * (B == "b2") + 1 * ((A == "a3") & (B == "b2")) +
    rnorm(nrow(grid), 0, 1) +
    rep(rnorm(n_sub, 0, 1.5), times = 6))
  grid
}

test_that("anova_within_two_way() returns a DMAR table that holds the 1.0.x data frame", {
  res <- anova_within_two_way(within_example(), "y", "A", "B", "subject")
  # DMAR 1.0.x returned a plain data.frame. The table is now a dmar_tbl,
  # and as.data.frame() gives that data frame back: the same columns,
  # values, and row names, with no attribute beside them.
  expect_identical(class(res), c("dmar_tbl", "data.frame"))
  plain <- as.data.frame(res)
  expect_identical(class(plain), "data.frame")
  expect_identical(unclass(res), unclass(plain))
  expect_setequal(names(attributes(plain)), c("names", "row.names", "class"))
})

test_that("anova_within_two_way() prints its p-values to four decimals", {
  res <- anova_within_two_way(within_example(), "y", "A", "B", "subject")
  # A width that keeps each row on one line.
  withr::local_options(width = 200)

  # The plain data frame of DMAR 1.0.x printed the p_value column in
  # scientific notation (1.183465e-09 for the unadjusted test of A).
  expect_true(any(grepl("1.183465e-09",
                        capture.output(print(as.data.frame(res))),
                        fixed = TRUE)))

  out <- capture.output(print(res))
  expect_false(any(grepl("[0-9]e[-+][0-9]", out)))
  # Whole-number degrees of freedom without a decimal part, adjusted ones
  # and the other values to 3 decimal places, the p-values to four
  # decimals, and NA where the unadjusted row has no epsilon.
  expect_true(any(grepl(
    "^ *A +none +60\\.273 +2 +22 +< 0\\.0001 +<NA> +0\\.846 *$", out)))
  expect_true(any(grepl(paste0(
    "^ *A +Greenhouse-Geisser +60\\.273 +1\\.830 +20\\.131 +< 0\\.0001 ",
    "+0\\.915 +0\\.846 *$"), out)))
  expect_true(any(grepl(
    "^ *A:B +none +0\\.973 +2 +22 +0\\.3936 +<NA> +0\\.081 *$", out)))
})

test_that("anova_within_two_way() answers tidy() and glance() through the default methods", {
  res <- anova_within_two_way(within_example(), "y", "A", "B", "subject")
  # One row per row of the table: the term joins the effect and the
  # adjustment, and the estimate is the F statistic.
  td <- generics::tidy(res)
  expect_identical(names(td), c("term", "estimate", "df_1", "df_2",
                                "p_value", "epsilon",
                                "eta_squared_partial"))
  expect_identical(td$term, paste(res$effect, res$adjustment, sep = ":"))
  expect_identical(td$term[1:4], c("A:none", "A:Greenhouse-Geisser",
                                   "A:Huynh-Feldt", "A:lower_bound"))
  expect_identical(td$estimate, res$F_value)
  for (nm in c("df_1", "df_2", "p_value", "epsilon", "eta_squared_partial"))
    expect_identical(td[[nm]], res[[nm]])
  # glance() is one row: the number of rows, with no level and no
  # bootstrap count, since no interval or bootstrap is formed.
  gl <- generics::glance(res)
  expect_identical(names(gl), c("n_terms", "conf_level", "B_used"))
  expect_identical(gl$n_terms, nrow(res))
  expect_true(is.na(gl$conf_level))
  expect_true(is.na(gl$B_used))
  expect_identical(dmar_metadata(res)$result, "dmar_tbl")
  expect_error(results_sentence(res), "confidence interval")
})
