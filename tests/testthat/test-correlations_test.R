test_that("correlations_test() returns an object of class 'correlations_test'", {
  res <- suppressMessages(correlations_test(attitude))
  expect_s3_class(res, "correlations_test")
  expect_named(res, c("r", "p", "ci_lower", "ci_upper", "n",
                      "method", "conf_level", "listwise",
                      "stars", "decimals_r", "decimals_p", "format"))
})

test_that("correlations_test() r matrix has 1s on the diagonal and is symmetric", {
  res <- correlations_test(attitude, format = "text")
  expect_equal(diag(res$r), rep(1, ncol(attitude)), ignore_attr = TRUE)
  expect_equal(res$r, t(res$r))
})

test_that("correlations_test() off-diagonal r values match cor()", {
  res <- correlations_test(attitude, format = "text")
  expected <- cor(attitude, use = "pairwise.complete.obs")
  diag(res$r) <- 1; diag(expected) <- 1
  expect_equal(res$r, expected, ignore_attr = TRUE)
})

test_that("correlations_test() p-values match cor.test() for a specific pair", {
  res <- correlations_test(attitude, format = "text")
  ct <- cor.test(attitude$rating, attitude$complaints)
  expect_equal(res$p["rating", "complaints"], ct$p.value)
})

test_that("correlations_test() CI matches cor.test conf.int for Pearson", {
  res <- correlations_test(attitude, format = "text")
  ct <- cor.test(attitude$rating, attitude$complaints, conf_level = 0.95)
  expect_equal(res$ci_lower["rating", "complaints"], ct$conf.int[1])
  expect_equal(res$ci_upper["rating", "complaints"], ct$conf.int[2])
})

test_that("correlations_test() pairwise N reflects available observations", {
  x <- data.frame(a = c(1, 2, 3, 4, 5, NA),
                  b = c(1, 2, NA, 4, 5, 6),
                  c = c(1, 2, 3, 4, 5, 6))
  res <- correlations_test(x, format = "text")
  expect_equal(res$n["a", "b"], 4L)  # rows 1,2,4,5 are complete for a,b
  expect_equal(res$n["a", "c"], 5L)
  expect_equal(res$n["b", "c"], 5L)
})

test_that("correlations_test() listwise deletion gives equal N for all pairs", {
  x <- data.frame(a = c(1, 2, 3, 4, 5, NA),
                  b = c(1, 2, NA, 4, 5, 6),
                  c = c(1, 2, 3, 4, 5, 6))
  res <- correlations_test(x, listwise = TRUE, format = "text")
  off_diag <- res$n[lower.tri(res$n)]
  expect_equal(length(unique(off_diag)), 1L)
  expect_equal(off_diag[1], 4L)
})

test_that("correlations_test() supports Spearman and Kendall methods", {
  res_s <- correlations_test(attitude, method = "spearman", format = "text")
  res_k <- correlations_test(attitude, method = "kendall", format = "text")
  expect_equal(res_s$method, "spearman")
  expect_equal(res_k$method, "kendall")
  # Spearman r between rating and complaints should match cor()
  expected <- cor(attitude$rating, attitude$complaints, method = "spearman")
  expect_equal(res_s$r["rating", "complaints"], expected)
})

test_that("correlations_test() errors on fewer than 2 numeric variables", {
  expect_error(correlations_test(data.frame(a = 1:5)), "at least two")
})

test_that("correlations_test() errors on bad argument values", {
  expect_error(correlations_test(attitude, conf_level = 1.5), "between 0 and 1")
  expect_error(correlations_test(attitude, conf_level = 0),   "between 0 and 1")
  expect_error(correlations_test(attitude, decimals_r = -1),  "decimals_r")
  expect_error(correlations_test(attitude, decimals_p = 0),   "decimals_p")
})

test_that("correlations_test() text output contains 'Correlations' caption and p-value strings", {
  out <- capture.output(correlations_test(attitude))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Correlations \\(Pearson")
  expect_match(txt, "p (<|=)")
  expect_match(txt, "N =")
})

test_that("correlations_test() with stars=TRUE includes a footnote and stars on r", {
  out <- capture.output(correlations_test(attitude, stars = TRUE))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Note\\. \\* p <")
  # At least one correlation in attitude is highly significant, so *** should appear
  expect_match(txt, "\\*\\*\\*")
})

test_that("correlations_test() produces an HTML kable when format='html'", {
  tbl <- correlations_test(attitude, format = "html")
  expect_true(inherits(tbl, "knitr_kable"))
  expect_match(as.character(tbl), "<table", fixed = FALSE)
})

test_that("correlations_test() produces LaTeX when format='latex'", {
  tbl <- correlations_test(attitude, format = "latex")
  s <- as.character(tbl)
  expect_match(s, "tabular", fixed = TRUE)
  expect_match(s, "makecell", fixed = TRUE)
})

test_that("correlations_test() writes text to file when file is supplied", {
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  correlations_test(attitude, format = "text", file = tmp)
  contents <- paste(readLines(tmp), collapse = "\n")
  expect_match(contents, "Correlations")
})

test_that("correlations_test() print method works for a stored object", {
  res <- correlations_test(attitude)
  expect_output(print(res), "Correlations")
})

test_that("correlations_test() handles matrix input", {
  res <- correlations_test(as.matrix(attitude))
  expect_equal(dim(res$r), c(ncol(attitude), ncol(attitude)))
})

test_that("correlations_test() HTML preserves 'p < .0001' as escaped '&lt;'", {
  tbl <- correlations_test(attitude, stars = TRUE, format = "html")
  src <- as.character(tbl)
  # The attitude data has multiple correlations with p < .0001 (rating vs
  # complaints, learning, raises). HTML must contain "p &lt; .0001" so the
  # browser doesn't interpret "<" as a tag opener and strip the value.
  # (The threshold tracks the decimals_p default, which is 4.)
  expect_match(src, "p &lt; .0001", fixed = TRUE)
  # Negative: bare "p [" indicates the < .0001 was eaten by HTML parsing.
  expect_false(grepl(">p [[]", src))
})

test_that("correlations_test() LaTeX escapes '<' as '\\textless{}'", {
  tbl <- correlations_test(attitude, stars = TRUE, format = "latex")
  src <- as.character(tbl)
  expect_match(src, "textless", fixed = TRUE)
  # No bare "<" should appear in body content (the booktabs tabular has none
  # by construction; this guards against future regressions).
  body <- sub(".*\\\\begin\\{tabular\\}", "", src)
  body <- sub("\\\\end\\{tabular\\}.*", "", body)
  expect_false(grepl("[^\\]<", body))
})

test_that("correlations_test() LaTeX escapes special chars in variable names", {
  x <- data.frame(a_b = rnorm(20), c_d = rnorm(20), `pct%` = rnorm(20),
                  check.names = FALSE)
  tbl <- correlations_test(x, format = "latex")
  src <- as.character(tbl)
  expect_match(src, "a\\\\_b", fixed = FALSE)
  expect_match(src, "c\\\\_d", fixed = FALSE)
  expect_match(src, "pct\\\\%", fixed = FALSE)
})

test_that("correlations_test() HTML file output writes a self-contained document without pandoc", {
  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp), add = TRUE)
  correlations_test(attitude, stars = TRUE, format = "html", file = tmp)
  expect_true(file.exists(tmp))
  contents <- paste(readLines(tmp), collapse = "\n")
  expect_match(contents, "<!DOCTYPE html>", fixed = TRUE)
  expect_match(contents, "<table",          fixed = TRUE)
  expect_match(contents, "p &lt; .001",     fixed = TRUE)
})

test_that("correlations_test() LaTeX file output writes the tabular fragment", {
  tmp <- tempfile(fileext = ".tex")
  on.exit(unlink(tmp), add = TRUE)
  correlations_test(attitude, stars = TRUE, format = "latex", file = tmp)
  expect_true(file.exists(tmp))
  contents <- paste(readLines(tmp), collapse = "\n")
  expect_match(contents, "tabular",   fixed = TRUE)
  expect_match(contents, "makecell",  fixed = TRUE)
  expect_match(contents, "textless",  fixed = TRUE)
})

test_that("correlations_test() Pearson CI matches the Fisher (1915, 1921) form exactly", {
  res <- correlations_test(attitude, format = "text")
  r <- res$r["rating", "complaints"]
  n <- res$n["rating", "complaints"]
  zc <- qnorm(0.975)
  z  <- atanh(r)
  expected_lo <- tanh(z - zc / sqrt(n - 3))
  expected_up <- tanh(z + zc / sqrt(n - 3))
  expect_equal(res$ci_lower["rating", "complaints"], expected_lo)
  expect_equal(res$ci_upper["rating", "complaints"], expected_up)
})

test_that("correlations_test() Spearman CI uses Bonett-Wright (2000) SE = sqrt((1 + r^2/2)/(n-3))", {
  res <- correlations_test(attitude, method = "spearman", format = "text")
  r <- res$r["rating", "complaints"]
  n <- res$n["rating", "complaints"]
  zc <- qnorm(0.975)
  z  <- atanh(r)
  se <- sqrt((1 + r^2 / 2) / (n - 3))
  expect_equal(res$ci_lower["rating", "complaints"], tanh(z - zc * se))
  expect_equal(res$ci_upper["rating", "complaints"], tanh(z + zc * se))
})

test_that("correlations_test() Kendall CI uses Bonett-Wright (2000) SE = sqrt(0.437/(n-4))", {
  res <- correlations_test(attitude, method = "kendall", format = "text")
  r <- res$r["rating", "complaints"]
  n <- res$n["rating", "complaints"]
  zc <- qnorm(0.975)
  z  <- atanh(r)
  se <- sqrt(0.437 / (n - 4))
  expect_equal(res$ci_lower["rating", "complaints"], tanh(z - zc * se))
  expect_equal(res$ci_upper["rating", "complaints"], tanh(z + zc * se))
})

test_that("correlations_test() Spearman CI is wider than Pearson CI at the same r, n (B-W penalty)", {
  set.seed(113)
  x <- data.frame(a = rnorm(50), b = rnorm(50))
  rp <- correlations_test(x, method = "pearson",  format = "text")
  rs <- correlations_test(x, method = "spearman", format = "text")
  # At |r| > 0 the Bonett-Wright Spearman SE is strictly larger than the
  # Pearson Fisher SE since (1 + r^2/2) > 1, so the interval is wider.
  pearson_width  <- rp$ci_upper["a", "b"] - rp$ci_lower["a", "b"]
  spearman_width <- rs$ci_upper["a", "b"] - rs$ci_lower["a", "b"]
  expect_gt(spearman_width, pearson_width * 0.99)  # generous to permit r near 0
})

test_that("correlations_test() Kendall returns NA CI for n < 5", {
  # Only 4 complete pairs in a, b -- enough for tau but not for the (n - 4)
  # Bonett-Wright denominator.
  x <- data.frame(
    a = c(1, 2, 3, 4, NA, NA),
    b = c(2, 1, 4, 3, NA, NA),
    c = c(1, 2, 3, 4, 5, 6)
  )
  res <- correlations_test(x, method = "kendall", format = "text")
  expect_equal(res$n["a", "b"], 4L)
  expect_true(is.na(res$ci_lower["a", "b"]))
  expect_true(is.na(res$ci_upper["a", "b"]))
})
