.mc_spatial <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
                 "t4_lozenges", "t25_paper_form_board_r", "t26_flags")

# lavaan's nine Holzinger and Swineford tests with values missing at random
# on six of them (the fixture of test-fraction_missing_information.R).
.mc_hs_mar <- function(seed = 113) {
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

test_that("the spatial tests are covered by Grant-White wherever t25 or t26 enters", {
  res <- missing_coverage(holzinger_swineford, items = .mc_spatial)
  expect_s3_class(res, "dmar_tbl")
  expect_identical(names(res), c("var_1", "var_2", "n_both", "coverage"))
  expect_identical(nrow(res), 21L)
  pair <- res$var_1 == "t25_paper_form_board_r" & res$var_2 == "t26_flags"
  expect_identical(res$n_both[pair], 145L)
  expect_equal(res$coverage[pair], 145 / 301)
  expect_equal(res$coverage[pair], 0.4817276, tolerance = 1e-7)
  # The first pair is the first item with itself, then the item order.
  expect_identical(res$var_1[1:3], rep("t1_visual_perception", 3))
  expect_identical(res$var_2[1:3], .mc_spatial[1:3])
  m <- attr(res, "matrix")
  expect_identical(dimnames(m), list(.mc_spatial, .mc_spatial))
  expect_true(isSymmetric(m))
  expect_equal(unname(diag(m)),
               colMeans(!is.na(holzinger_swineford[.mc_spatial])),
               ignore_attr = TRUE)
  expect_identical(attr(res, "N_used"), 301L)
  expect_identical(attr(res, "N_all_missing"), 0L)
  expect_equal(attr(res, "min_coverage"), 145 / 301)
})

test_that("the matrix is lavaan's coverage, whose denominator leaves out empty rows", {
  skip_if_not_installed("lavaan")
  d <- .mc_hs_mar()
  d[1:3, ] <- NA
  res <- missing_coverage(d)
  fit <- suppressWarnings(lavaan::cfa(
    "visual =~ x1 + x2 + x3
     textual =~ x4 + x5 + x6
     speed =~ x7 + x8 + x9",
    data = d, std.lv = TRUE, missing = "ml"))
  lav <- unclass(lavaan::lavInspect(fit, "coverage"))
  expect_identical(attr(res, "N_used"), 298L)
  expect_identical(attr(res, "N_all_missing"), 3L)
  expect_equal(attr(res, "matrix"), lav, tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_identical(dimnames(attr(res, "matrix")), dimnames(lav))
  # Dividing by every row instead would not reproduce lavaan.
  all_rows <- crossprod(!is.na(as.matrix(d))) / nrow(d)
  expect_gt(max(abs(all_rows - lav)), 0.005)
})

test_that("complete data are fully covered and every case missing everything is refused", {
  res <- missing_coverage(holzinger_swineford, items = .mc_spatial[1:4])
  expect_true(all(res$coverage == 1))
  expect_identical(nrow(res), 10L)
  expect_identical(attr(res, "min_coverage"), 1)

  d <- data.frame(a = c(NA, NA), b = c(NA, NA))
  expect_error(missing_coverage(d), "Every case is missing every item")
  expect_error(missing_coverage(list(a = 1)), "'data' must be a data.frame")
  one <- missing_coverage(data.frame(a = c(1, NA, 3)))
  expect_identical(nrow(one), 1L)
  expect_equal(one$coverage, 1)
  expect_identical(attr(one, "N_all_missing"), 1L)
})
