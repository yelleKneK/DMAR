.mp_spatial <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
                 "t4_lozenges", "t25_paper_form_board_r", "t26_flags")

# lavaan's nine Holzinger and Swineford tests with values missing at random
# on six of them (the fixture of test-fraction_missing_information.R).
.mp_hs_mar <- function(seed = 113) {
  d <- lavaan::HolzingerSwineford1939[, c(paste0("x", 1:9), "ageyr",
                                          "agemo")]
  set.seed(seed)
  n <- nrow(d)
  age <- d$ageyr + d$agemo / 12
  zs <- function(v) as.numeric(scale(v))
  drive <- list(x2 = zs(d$x1), x5 = zs(d$x4), x8 = zs(d$x7),
                x3 = zs(age), x6 = zs(age), x9 = zs(age))
  for (v in names(drive)) {
    p <- stats::plogis(-1.3 + 1.2 * drive[[v]])
    d[[v]][stats::runif(n) < p] <- NA
  }
  d[paste0("x", 1:9)]
}

test_that("the spatial tests show the two patterns the design created", {
  res <- missing_patterns(holzinger_swineford, items = .mp_spatial)
  expect_s3_class(res, "dmar_tbl")
  expect_identical(names(res), c("pattern", "missing_items", "n",
                                 "proportion", "n_observed", "n_missing"))
  expect_identical(res$pattern, c("111100", "111111"))
  expect_identical(res$missing_items,
                   c("t25_paper_form_board_r, t26_flags", ""))
  expect_identical(res$n, c(156L, 145L))
  expect_equal(res$proportion, c(156, 145) / 301)
  expect_identical(res$n_observed, c(4L, 6L))
  expect_identical(res$n_missing, c(2L, 0L))
  expect_identical(sum(res$n), nrow(holzinger_swineford))
  expect_identical(attr(res, "N"), 301L)
  expect_identical(attr(res, "N_complete"), 145L)
  expect_identical(attr(res, "N_all_missing"), 0L)
  expect_identical(attr(res, "n_patterns"), 2L)
  expect_identical(attr(res, "items"), .mp_spatial)
  m <- attr(res, "matrix")
  expect_identical(dim(m), c(2L, 6L))
  expect_identical(colnames(m), .mp_spatial)
  expect_identical(m[1, ], c(rep(TRUE, 4), FALSE, FALSE),
                   ignore_attr = TRUE)
})

test_that("the pattern with nothing observed is kept, and the rest match lavaan", {
  skip_if_not_installed("lavaan")
  d <- .mp_hs_mar()
  d[1:3, ] <- NA
  res <- missing_patterns(d)
  empty <- res[res$n_observed == 0L, ]
  expect_identical(nrow(empty), 1L)
  expect_identical(empty$n, 3L)
  expect_identical(empty$pattern, "000000000")
  expect_identical(attr(res, "N_all_missing"), 3L)
  expect_identical(sum(res$n), nrow(d))

  fit <- suppressWarnings(lavaan::cfa(
    "visual =~ x1 + x2 + x3
     textual =~ x4 + x5 + x6
     speed =~ x7 + x8 + x9",
    data = d, std.lv = TRUE, missing = "ml"))
  lav <- lavaan::lavInspect(fit, "patterns")
  expect_identical(attr(res, "n_patterns"), nrow(lav) + 1L)
  ours <- attr(res, "matrix")[res$n_observed > 0L, , drop = FALSE]
  key <- function(m) sort(apply(m, 1L, function(r)
    paste(as.integer(r), collapse = "")))
  expect_identical(key(ours), key(unclass(lav)))
  expect_identical(colnames(ours), colnames(lav))
})

test_that("the rows are sorted by frequency and then by missing items, whatever the row order", {
  d <- data.frame(a = c(1, NA, 3, NA, 5, 6, NA),
                  b = c(1, 2, NA, NA, 5, 6, NA),
                  c = c(1, 2, 3, 4, 5, 6, NA))
  res <- missing_patterns(d)
  # Three complete rows lead; the four single rows follow by the number of
  # items missing, and the two with one missing by the pattern string
  # (more items observed earlier in the item order first).
  expect_identical(res$pattern, c("111", "101", "011", "001", "000"))
  expect_identical(res$n, c(3L, 1L, 1L, 1L, 1L))
  expect_identical(res$n_missing, c(0L, 1L, 1L, 2L, 3L))
  expect_identical(res$missing_items, c("", "b", "a", "a, b", "a, b, c"))
  shuffled <- missing_patterns(d[c(7, 3, 1, 5, 2, 6, 4), ])
  expect_identical(shuffled$pattern, res$pattern)
  expect_identical(shuffled$n, res$n)
  expect_identical(attr(shuffled, "matrix"), attr(res, "matrix"))
})

test_that("any column type works and one item keeps the shape", {
  d <- data.frame(num = c(1, NA, 3),
                  chr = c("a", "b", NA),
                  fac = factor(c(NA, "x", "y")),
                  lgl = c(TRUE, FALSE, NA),
                  stringsAsFactors = FALSE)
  res <- missing_patterns(d)
  expect_identical(attr(res, "n_patterns"), 3L)
  expect_identical(sum(res$n), 3L)
  one <- missing_patterns(d, items = "num")
  expect_identical(names(one), names(res))
  expect_identical(one$pattern, c("1", "0"))
  expect_identical(one$missing_items, c("", "num"))
  single_row <- missing_patterns(d[2, ])
  expect_identical(single_row$pattern, "0111")
})

test_that("bad input is refused", {
  expect_error(missing_patterns(as.matrix(holzinger_swineford)),
               "'data' must be a data.frame")
  expect_error(missing_patterns(holzinger_swineford, items = "no_such_item"),
               "not in 'data': no_such_item")
  expect_error(missing_patterns(holzinger_swineford,
                                items = c("t2_cubes", "t2_cubes")),
               "more than once")
  expect_error(missing_patterns(holzinger_swineford, items = 3),
               "character vector")
  expect_error(missing_patterns(holzinger_swineford[0, ]), "no rows")
})
