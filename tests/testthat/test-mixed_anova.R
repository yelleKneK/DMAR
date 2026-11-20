test_that("mixed_anova() one-way returns expected columns", {
  set.seed(113)
  d <- data.frame(A = factor(rep(1:3, each = 10)),
                  y = rnorm(30) + 0.5 * rep(1:3, each = 10))
  res <- mixed_anova(d, "y", "A", A_type = "fixed")
  expect_setequal(colnames(res),
                  c("effect", "ss", "df", "ms", "denominator",
                    "F_value", "p_value"))
  expect_equal(res$effect, "A")
})

test_that("mixed_anova() two-way fixed-fixed denominator is MS_within", {
  set.seed(113)
  d <- expand.grid(A = factor(1:3), B = factor(1:2), rep = 1:5)
  d$y <- rnorm(nrow(d))
  res <- mixed_anova(d, "y", "A", "B", A_type = "fixed", B_type = "fixed")
  expect_true(all(res$denominator == "MS_within"))
})

test_that("mixed_anova() A fixed / B random uses MS_AB for A", {
  set.seed(113)
  d <- expand.grid(A = factor(1:3), B = factor(1:4), rep = 1:5)
  d$y <- rnorm(nrow(d))
  res <- mixed_anova(d, "y", "A", "B", A_type = "fixed", B_type = "random")
  expect_equal(res$denominator[res$effect == "A"],   "MS_AB")
  expect_equal(res$denominator[res$effect == "B"],   "MS_within")
  expect_equal(res$denominator[res$effect == "A:B"], "MS_within")
})

test_that("mixed_anova() both random uses MS_AB for both main effects", {
  set.seed(113)
  d <- expand.grid(A = factor(1:3), B = factor(1:4), rep = 1:4)
  d$y <- rnorm(nrow(d))
  res <- mixed_anova(d, "y", "A", "B", A_type = "random", B_type = "random")
  expect_equal(res$denominator[res$effect == "A"], "MS_AB")
  expect_equal(res$denominator[res$effect == "B"], "MS_AB")
})

# The design of the help example: A fixed with three levels, B random with
# five, and four replicates per cell.
mixed_example <- function() {
  set.seed(113)
  grid <- expand.grid(A = factor(1:3), B = factor(1:5), rep = 1:4)
  grid$y <- with(grid, 0.8 * as.integer(A) +
                   rep(rnorm(5, 0, 1.5), each = 1)[as.integer(B)] +
                   rnorm(nrow(grid), 0, 1))
  grid
}

test_that("mixed_anova() returns a DMAR table that holds the 1.0.x data frame", {
  grid <- mixed_example()
  two_way <- mixed_anova(grid, "y", "A", "B", A_type = "fixed",
                         B_type = "random")
  one_way <- mixed_anova(grid, "y", "A")
  for (res in list(two_way, one_way)) {
    # DMAR 1.0.x returned a plain data.frame carrying one attribute of its
    # own, sum_of_squares_type. The table is now a dmar_tbl, and
    # as.data.frame() gives that data frame back unchanged.
    expect_identical(class(res), c("dmar_tbl", "data.frame"))
    plain <- as.data.frame(res)
    expect_identical(class(plain), "data.frame")
    expect_identical(unclass(res), unclass(plain))
    expect_setequal(names(attributes(plain)),
                    c("names", "row.names", "class", "sum_of_squares_type"))
    expect_identical(attr(res, "sum_of_squares_type", exact = TRUE), 3)
  }
})

test_that("mixed_anova() prints its p-values to four decimals", {
  res <- mixed_anova(mixed_example(), "y", "A", "B", A_type = "fixed",
                     B_type = "random")
  # A width that keeps each row on one line.
  withr::local_options(width = 200)

  # The plain data frame of DMAR 1.0.x printed the p_value column in
  # scientific notation: 2.148782e-05, 1.817787e-12, and 9.565106e-01.
  expect_true(any(grepl("2.148782e-05",
                        capture.output(print(as.data.frame(res))),
                        fixed = TRUE)))

  out <- capture.output(print(res))
  expect_false(any(grepl("[0-9]e[-+][0-9]", out)))
  # The sums of squares, mean squares, and F to 3 decimal places, the
  # degrees of freedom as whole numbers, and the p-values to four decimals.
  expect_true(any(grepl(
    "^ *A +42\\.288 +2 +21\\.144 +MS_AB +54\\.751 +< 0\\.0001 *$", out)))
  expect_true(any(grepl(
    "^ *B +153\\.418 +4 +38\\.354 +MS_within +31\\.236 +< 0\\.0001 *$",
    out)))
  expect_true(any(grepl(
    "^ *A:B +3\\.089 +8 +0\\.386 +MS_within +0\\.315 +0\\.9565 *$", out)))
})

test_that("mixed_anova() answers tidy() and glance() through the default methods", {
  res <- mixed_anova(mixed_example(), "y", "A", "B", A_type = "fixed",
                     B_type = "random")
  # One row per effect: the term joins the effect and the mean square it
  # is tested against, and the estimate is the sum of squares.
  td <- generics::tidy(res)
  expect_identical(names(td), c("term", "estimate", "df", "ms", "F_value",
                                "p_value"))
  expect_identical(td$term, c("A:MS_AB", "B:MS_within", "A:B:MS_within"))
  expect_identical(td$estimate, res$ss)
  for (nm in c("df", "ms", "F_value", "p_value"))
    expect_identical(td[[nm]], res[[nm]])
  # glance() is one row: the number of rows, with no level and no
  # bootstrap count, since no interval or bootstrap is formed.
  gl <- generics::glance(res)
  expect_identical(names(gl), c("n_terms", "conf_level", "B_used"))
  expect_identical(gl$n_terms, nrow(res))
  expect_true(is.na(gl$conf_level))
  expect_true(is.na(gl$B_used))
  expect_error(results_sentence(res), "confidence interval")
})

test_that("mixed_anova() tables reach dmar_metadata() and as_kable()", {
  res <- mixed_anova(mixed_example(), "y", "A", "B", A_type = "fixed",
                     B_type = "random")
  # Both take only a DMAR table, so the plain data frame is refused; DMAR
  # 1.0.x returned that data frame, and its as_kable() had no method for it.
  expect_identical(dmar_metadata(res)$result, "dmar_tbl")
  expect_error(dmar_metadata(as.data.frame(res)), "dmar_tbl")
  skip_if_not_installed("knitr")
  k <- as_kable(res, format = "pipe")
  expect_s3_class(k, "knitr_kable")
  # The kable holds the display strings, the p-values to four decimals.
  txt <- paste(as.character(k), collapse = "\n")
  expect_match(txt, "< 0.0001", fixed = TRUE)
  expect_match(txt, "0.9565", fixed = TRUE)
  expect_false(grepl("[0-9]e[-+][0-9]", txt))
  expect_error(as_kable(as.data.frame(res)), "no applicable method")
})
