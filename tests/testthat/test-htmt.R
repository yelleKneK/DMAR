# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

.htmt_data <- function(n = 300) {
  set.seed(113)
  f1 <- rnorm(n); f2 <- 0.3 * f1 + sqrt(1 - 0.09) * rnorm(n)
  data.frame(
    a1 = .8 * f1 + rnorm(n, 0, .6), a2 = .7 * f1 + rnorm(n, 0, .7),
    a3 = .6 * f1 + rnorm(n, 0, .8),
    b1 = .8 * f2 + rnorm(n, 0, .6), b2 = .7 * f2 + rnorm(n, 0, .7),
    b3 = .6 * f2 + rnorm(n, 0, .8))
}

test_that("htmt() matches the Henseler definition by hand", {
  d <- .htmt_data()
  blocks <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
  res <- htmt(d, blocks)
  expect_s3_class(res, "dmar_tbl")
  expect_equal(nrow(res), 1)
  # By hand: the ratio of the mean absolute between-construct correlation
  # to the geometric mean of the two mean absolute within-construct
  # correlations. It is the definition, exactly, so no second
  # implementation is needed to confirm it. (semTools::htmt() with
  # htmt2 = FALSE and absolute = TRUE gives 0.2691979 on these data; the
  # check was run once by hand and is not repeated here, so the test does
  # not depend on semTools.)
  R <- abs(cor(d))
  hetero <- mean(R[1:3, 4:6])
  mono_a <- mean(R[1:3, 1:3][lower.tri(diag(3))])
  mono_b <- mean(R[4:6, 4:6][lower.tri(diag(3))])
  expect_equal(res$htmt, hetero / sqrt(mono_a * mono_b))
  expect_equal(res$htmt, 0.2691979, tolerance = 1e-6)
})

test_that("an unrecoded reverse keyed item leaves htmt() where recoding it would", {
  # The 1.1.0 review: htmt() took absolute values of the heterotrait
  # correlations but kept the monotrait ones signed, so one reverse keyed
  # item in a block of four (a4 below) shrank the denominator and not the
  # numerator, and on these data the ratio was 1.346 for two constructs
  # whose latent correlation is about .47, a discriminant validity failure
  # that is not there; with a block of three such items the call stopped.
  # Every correlation is now absolute, so reversing an item changes
  # nothing.
  set.seed(113)
  n <- 300
  f1 <- rnorm(n)
  f2 <- 0.5 * f1 + sqrt(0.75) * rnorm(n)
  d <- data.frame(
    a1 = .8 * f1 + rnorm(n, 0, .6), a2 = .7 * f1 + rnorm(n, 0, .7),
    a3 = .8 * f1 + rnorm(n, 0, .6), a4 = -.6 * f1 + rnorm(n, 0, .8),
    b1 = .8 * f2 + rnorm(n, 0, .6), b2 = .7 * f2 + rnorm(n, 0, .7),
    b3 = .6 * f2 + rnorm(n, 0, .8))
  blocks <- list(A = c("a1", "a2", "a3", "a4"), B = c("b1", "b2", "b3"))
  recoded <- d
  recoded$a4 <- -recoded$a4
  as_keyed <- htmt(d, blocks)$htmt
  expect_equal(as_keyed, htmt(recoded, blocks)$htmt)
  # By hand, on the signed correlations of the recoded data.
  R <- cor(recoded)
  a <- blocks$A; b <- blocks$B
  by_hand <- mean(R[a, b]) /
    sqrt(mean(R[a, a][lower.tri(diag(4))]) * mean(R[b, b][lower.tri(diag(3))]))
  expect_equal(as_keyed, by_hand)
  expect_equal(as_keyed, 0.4630144, tolerance = 1e-6)
  expect_lt(as_keyed, 1)
  # The block of three that stopped through 1.0.1 now has a ratio, the
  # one its recoded items give.
  three <- list(A = c("a1", "a2", "a4"), B = c("b1", "b2", "b3"))
  expect_equal(htmt(d, three)$htmt, htmt(recoded, three)$htmt)
  # Reversing every item of a block, or the other block, changes nothing
  # either.
  flipped <- d
  flipped[b] <- -flipped[b]
  expect_equal(htmt(flipped, blocks)$htmt, as_keyed)
})

test_that("htmt() bootstrap upper bound behaves and is seeded", {
  skip_on_cran()
  d <- .htmt_data(150)
  blocks <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
  r1 <- htmt(d, blocks, B = 200, seed = 7)
  r2 <- htmt(d, blocks, B = 200, seed = 7)
  expect_identical(r1$ci_upper, r2$ci_upper)
  expect_gt(r1$ci_upper, r1$htmt)
  expect_identical(attr(r1, "conf_level"), 0.95)
})

test_that("htmt() validates and handles three constructs", {
  d <- .htmt_data()
  set.seed(7)
  f3 <- rnorm(nrow(d))
  d$c1 <- .7 * f3 + rnorm(nrow(d), 0, .7)
  d$c2 <- .6 * f3 + rnorm(nrow(d), 0, .8)
  res <- htmt(d, list(A = c("a1", "a2"), B = c("b1", "b2"),
                      C = c("c1", "c2")))
  expect_equal(nrow(res), 3)             # all pairs
  expect_error(htmt(d, list(A = "a1", B = c("b1", "b2"))), "two or more")
  expect_error(htmt(d, list(c("a1", "a2"), B = c("b1", "b2"))), "named")
  expect_error(htmt(d, list(A = c("a1", "a2"), B = c("b1", "b2")),
                    B = 1), "single whole number",
               class = "dmar_count_invalid")
})

test_that("a degenerate resample is dropped, not fatal", {
  skip_on_cran()
  # Item b2 takes the value 1 for 3 of 40 cases and 0 otherwise, so about
  # one resample in twenty leaves it constant and its correlations
  # undefined; those resamples must be dropped with one warning that
  # gives the counts, not abort the bootstrap.
  set.seed(113)
  n <- 40
  f <- rnorm(n)
  g <- rnorm(n)
  d <- data.frame(
    a1 = f + rnorm(n, sd = 0.4), a2 = f + rnorm(n, sd = 0.4),
    b1 = g, b2 = as.numeric(rank(g) > n - 3)
  )
  blocks <- list(A = c("a1", "a2"), B = c("b1", "b2"))
  w <- testthat::capture_warnings(
    res <- htmt(d, blocks = blocks, B = 300, seed = 113))
  expect_length(w, 1L)
  expect_match(w, paste0("^htmt\\(\\): [0-9]+ of 300 bootstrap resamples ",
                         "\\([0-9.]+%\\) left HTMT undefined"))
  expect_match(w, "that remained", fixed = TRUE)
  expect_true(all(is.finite(res$ci_upper)))
  # The resamples the bound used are the ones that remained, and the
  # warning's count is the difference.
  used <- attr(res, "B_used")
  expect_lt(used, 300L)
  expect_gte(used, 100L)
  expect_match(w, paste0("^htmt\\(\\): ", 300L - used, " of 300 "))
  expect_match(w, paste0("from the ", used, " that remained"), fixed = TRUE)
  expect_identical(attr(res, "n_dropped"), 300L - used)
  # (DMAR 1.0.x stopped when fewer than 100 resamples remained; the bound
  # is now withheld under a classed warning, tested in
  # test-dmar_bootstrap.R.)
  # A constant item in the observed data has no ratio, and the error says
  # why (cor()'s own warning is not repeated).
  d$b2 <- 0
  expect_error(expect_no_warning(htmt(d, blocks)),
               "a correlation among their items cannot be computed")
})

test_that("htmt() returns the same columns whatever B, with NA where nothing is computed", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # The 1.1.0 review: ci_upper existed only when B > 0, so the shape of
  # the table depended on an argument. It is now always present.
  d <- .htmt_data(150)
  blocks <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
  cols <- c("construct_1", "construct_2", "htmt", "ci_upper")
  point <- htmt(d, blocks)
  boot <- htmt(d, blocks, B = 200, seed = 7)
  expect_identical(names(point), cols)
  expect_identical(names(boot), cols)
  expect_identical(vapply(point, typeof, ""), vapply(boot, typeof, ""))
  expect_true(all(is.na(point$ci_upper)))
  expect_identical(point$htmt, boot$htmt)
  expect_true(all(is.finite(boot$ci_upper)))
  # The confidence level and the resamples used travel with a bound only.
  expect_null(attr(point, "conf_level"))
  expect_null(attr(point, "B_used"))
  expect_identical(attr(boot, "B_used"), 200L)
  # The broom view has the same shape either way.
  expect_identical(names(generics::tidy(point)), names(generics::tidy(boot)))
})

test_that("a pair whose ratio is undefined is NA, with a warning, beside the pairs that are defined", {
  skip_on_cran()
  # The 1.1.0 review: an undefined pair stopped the call, so with three
  # blocks the defined pairs were lost, and fa_confirmatory(output = "measurement")
  # lost its whole table. Through 1.0.1 a pair whose items were never
  # observed together across its blocks was NA without a word.
  set.seed(113)
  n <- 300
  f1 <- rnorm(n); f2 <- 0.4 * f1 + sqrt(.84) * rnorm(n)
  f3 <- 0.3 * f1 + sqrt(.91) * rnorm(n)
  d <- data.frame(
    a1 = .8 * f1 + rnorm(n, 0, .6), a2 = .7 * f1 + rnorm(n, 0, .7),
    a3 = .7 * f1 + rnorm(n, 0, .7),
    b1 = .8 * f2 + rnorm(n, 0, .6), b2 = .7 * f2 + rnorm(n, 0, .7),
    b3 = .7 * f2 + rnorm(n, 0, .7),
    c1 = .8 * f3 + rnorm(n, 0, .6), c2 = .7 * f3 + rnorm(n, 0, .7),
    c3 = .7 * f3 + rnorm(n, 0, .7))
  form <- rep(1:2, length.out = n)
  d$a1[form == 1] <- NA
  d$b1[form == 2] <- NA   # a1 and b1 are never observed together
  bl <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"),
             C = c("c1", "c2", "c3"))
  w <- testthat::capture_warnings(res <- htmt(d, bl))
  expect_identical(w, paste0(
    "HTMT is undefined for 'A' and 'B': a correlation among their items ",
    "cannot be computed (an item with no variance, or two items never ",
    "observed together). Its row is NA."))
  expect_true(is.na(res$htmt[1L]))
  # The defined pairs are those of a call on their two blocks alone.
  expect_equal(res$htmt[2L], htmt(d, bl[c("A", "C")])$htmt, tolerance = 0)
  expect_equal(res$htmt[3L], htmt(d, bl[c("B", "C")])$htmt, tolerance = 0)
  # The bootstrap bounds the defined pairs and leaves the other NA.
  wb <- testthat::capture_warnings(
    rb <- htmt(d, bl, B = 100, seed = 113))
  expect_length(wb, 1L)
  expect_identical(is.na(rb$ci_upper), c(TRUE, FALSE, FALSE))
  expect_identical(attr(rb, "B_used"), 100L)
  # With nothing defined the call stops, and says why.
  expect_error(htmt(d, bl[c("A", "B")]),
               "HTMT is undefined for 'A' and 'B': a correlation",
               class = "dmar_htmt_undefined")
})

test_that("the bootstrap resamples the absolute ratio it reports", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # The 1.1.0 review: every observed within-block correlation of block B is
  # positive (.04, .11, .20), so the point estimate is that of 1.0.x, but a
  # resample can draw a negative one. Through 1.0.1 a resample's
  # within-block correlations entered with their signs, and a resample
  # whose average was not positive was dropped: the bound was .869 from
  # 198 resamples. Each resample now gives the ratio that is reported.
  set.seed(4)
  n <- 80
  f <- rnorm(n); g <- rnorm(n)
  d <- data.frame(a1 = f + rnorm(n, sd = .5), a2 = f + rnorm(n, sd = .5),
                  b1 = .45 * g + rnorm(n), b2 = .45 * g + rnorm(n),
                  b3 = .45 * g + rnorm(n))
  bl <- list(A = c("a1", "a2"), B = c("b1", "b2", "b3"))
  R <- cor(d)
  expect_true(all(R[3:5, 3:5][lower.tri(diag(3))] > 0))
  res <- htmt(d, bl, B = 200, seed = 113)
  ratio <- function(dd) {
    R <- abs(cor(dd))
    mean(R[bl$A, bl$B]) /
      sqrt(mean(R[bl$A, bl$A][lower.tri(diag(2))]) *
             mean(R[bl$B, bl$B][lower.tri(diag(3))]))
  }
  boots <- withr::with_seed(113, replicate(200, ratio(
    d[sample.int(n, n, replace = TRUE), ])))
  expect_equal(res$htmt, ratio(d), tolerance = 1e-12)
  expect_equal(res$ci_upper, unname(stats::quantile(boots, 0.95)),
               tolerance = 1e-12)
  expect_equal(res$ci_upper, 0.7278264, tolerance = 1e-6)
  expect_identical(attr(res, "B_used"), 200L)
})

# Nine items of three blocks, n = 400: the fixture of the 2026-11-06 deep
# review (finding F13).
.htmt_review_data <- function() {
  withr::with_seed(113, {
    n <- 400L
    g <- rnorm(n)
    groups <- matrix(rnorm(n * 3), n, 3)
    X <- sapply(seq_len(9), function(j) .5 * g + .4 * groups[, ceiling(j / 3)] +
                  sqrt(.59) * rnorm(n))
    colnames(X) <- paste0("x", seq_len(9))
    as.data.frame(X)
  })
}

test_that("blocks that repeat, share, or rename items are refused before any correlation or resample", {
  skip_on_cran()
  # The 2026-11-06 deep review (F13): on these data the blocks
  # A = x1-x3 and B = x4-x6 give HTMT .615; x1 listed twice in A gave
  # .544 and x3 placed in both blocks gave .944, each without a word.
  # The blocks are now checked to be disjoint sets of distinct items under
  # distinct names before anything is computed, whatever B.
  d <- .htmt_review_data()
  bl <- list(A = c("x1", "x2", "x3"), B = c("x4", "x5", "x6"))
  # Disjoint blocks are unchanged: the definition by direct algebra.
  R <- abs(cor(d))
  by_hand <- mean(R[1:3, 4:6]) /
    sqrt(mean(R[1:3, 1:3][lower.tri(diag(3))]) * mean(R[4:6, 4:6][lower.tri(diag(3))]))
  expect_equal(htmt(d, bl)$htmt, by_hand, tolerance = 1e-12)
  expect_equal(by_hand, 0.6154196, tolerance = 1e-6)
  refused <- function(blocks, reason, at, items, data = d) {
    for (B in c(0, 200)) {
      set.seed(99)
      before <- get(".Random.seed", envir = globalenv())
      e <- tryCatch(htmt(data, blocks, B = B), error = function(e) e)
      expect_s3_class(e, c("dmar_htmt_blocks_invalid", "dmar_refusal", "error"))
      expect_identical(e$code, "htmt_blocks_invalid")
      expect_identical(e$reason, reason)
      expect_identical(e$blocks, at)
      expect_identical(e$items, items)
      for (nm in c(at, items)) expect_match(conditionMessage(e), nm, fixed = TRUE)
      # No resample was drawn: the random number stream is where it was.
      expect_identical(get(".Random.seed", envir = globalenv()), before)
    }
    # A column of the blocks that cor() cannot read would stop the call
    # with cor()'s own error if a correlation were computed first.
    unreadable <- data
    unreadable[[if (length(items)) items[1L] else blocks[[1L]][1L]]] <- "not a number"
    expect_error(htmt(unreadable, blocks), class = "dmar_htmt_blocks_invalid")
  }
  refused(list(A = c("x1", "x1", "x2", "x3"), B = bl$B), "repeated", "A", "x1")
  refused(list(A = c("x1", "x1", "x2", "x3"), B = c("x4", "x5", "x5", "x6", "x6")),
          "repeated", c("A", "B"), c("x1", "x5", "x6"))
  refused(list(A = bl$A, B = c("x3", "x4", "x5", "x6")), "shared", c("A", "B"), "x3")
  refused(list(A = bl$A, B = c("x4", "x5", "x6"), C = c("x3", "x7", "x8")),
          "shared", c("A", "C"), "x3")
  refused(list(A = bl$A, B = bl$A), "shared", c("A", "B"), c("x1", "x2", "x3"))
  refused(list(A = bl$A, A = bl$B), "construct_names", "A", character(0))
  dup <- d
  names(dup)[7] <- "x1"
  refused(bl, "data_names", "A", "x1", data = dup)
  # The messages say what is wrong and what to do.
  expect_error(htmt(d, list(A = c("x1", "x1", "x2", "x3"), B = bl$B)),
               "x1 is listed more than once in block 'A'", fixed = TRUE)
  expect_error(htmt(d, list(A = bl$A, B = c("x3", "x4", "x5", "x6"))),
               "x3 is in blocks 'A' and 'B'", fixed = TRUE)
  expect_error(htmt(d, list(A = bl$A, A = bl$B)),
               "'A' names more than one block", fixed = TRUE)
  # A block of two entries that are one item is a repeated item, not a
  # block of two.
  expect_error(htmt(d, list(A = c("x1", "x1"), B = bl$B)),
               class = "dmar_htmt_blocks_invalid")
  # The checks of form are unchanged: a block that is not two or more
  # columns of data, or a list that is not named.
  expect_error(htmt(d, list(A = c("x1", "x99"), B = bl$B)),
               "Block 'A' must name two or more columns", class = "dmar_argument_invalid")
  expect_error(htmt(d, stats::setNames(list(bl$A, bl$B), c("A", NA))),
               "named list", class = "dmar_argument_invalid")
})

test_that("disjoint blocks keep their ratios and bootstrap bounds", {
  skip_on_cran()  # 200 resamples, recomputed here by the oracle
  # Three disjoint blocks: each pair is its definition by direct algebra,
  # and its bound is the 95th percentile of the ratios of the same
  # resamples, drawn here by sample.int() after the seed.
  d <- .htmt_review_data()
  bl <- list(A = c("x1", "x2", "x3"), B = c("x4", "x5", "x6"), C = c("x7", "x8", "x9"))
  ratio <- function(dd, a, b) {
    R <- abs(cor(dd))
    mean(R[a, b]) / sqrt(mean(R[a, a][lower.tri(diag(3))]) *
                           mean(R[b, b][lower.tri(diag(3))]))
  }
  pairs <- utils::combn(names(bl), 2, simplify = FALSE)
  res <- htmt(d, bl, B = 200, seed = 113)
  expect_identical(res$construct_1, c("A", "A", "B"))
  expect_identical(res$construct_2, c("B", "C", "C"))
  expect_equal(res$htmt, vapply(pairs, function(p) ratio(d, bl[[p[1]]], bl[[p[2]]]),
                                numeric(1)), tolerance = 1e-12)
  boots <- withr::with_seed(113, vapply(seq_len(200), function(i) {
    dd <- d[sample.int(nrow(d), nrow(d), replace = TRUE), ]
    vapply(pairs, function(p) ratio(dd, bl[[p[1]]], bl[[p[2]]]), numeric(1))
  }, numeric(3)))
  expect_equal(res$ci_upper, apply(boots, 1, stats::quantile, probs = 0.95, names = FALSE),
               tolerance = 1e-12)
  expect_identical(attr(res, "B_used"), 200L)
})
