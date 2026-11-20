.cmar_columns <- c("marker_r", "variable_1", "variable_2", "r", "estimate",
                   "ci_lower", "ci_upper", "t_value", "df", "p_value")

.cmar_R <- matrix(c(1, .50, .40, .30,
                    .50, 1, .45, .20,
                    .40, .45, 1, .35,
                    .30, .20, .35, 1), 4, 4,
                  dimnames = list(c("a", "b", "c", "d"),
                                  c("a", "b", "c", "d")))

test_that("the table has one row per marker value and pair, in a fixed order", {
  res <- common_method_adjusted_r(.cmar_R, N = 120,
                                  marker_r = c(0.05, 0.10, 0.20))
  expect_s3_class(res, "dmar_common_method_adjusted_r")
  expect_s3_class(res, "dmar_tbl")
  expect_identical(names(res), .cmar_columns)
  expect_identical(nrow(res), 18L)
  expect_identical(res$marker_r, rep(c(0.05, 0.10, 0.20), each = 6))
  expect_identical(res$variable_1[1:6], c("a", "a", "a", "b", "b", "c"))
  expect_identical(res$variable_2[1:6], c("b", "c", "d", "c", "d", "d"))
  expect_identical(res$r[1:6], c(.50, .40, .30, .45, .20, .35))
  expect_true(all(res$df == 117))
  expect_identical(attr(res, "N"), 120)
  expect_identical(attr(res, "conf_level"), 0.95)
  expect_identical(attr(res, "marker_rule"), "supplied")
  adj <- attr(res, "adjusted")
  expect_identical(names(adj), c("0.05", "0.1", "0.2"))
  expect_output(print(res), "Confidence level: 95%")
})

test_that("the adjustment is common_method_marker()'s, proxy included", {
  for (m in c(0.10, 0.25)) {
    res <- common_method_adjusted_r(.cmar_R, N = 120, marker_r = m)
    marker <- common_method_marker(.cmar_R, marker_r = m)
    expect_equal(attr(res, "adjusted")[[1L]], attr(marker, "adjusted"),
                 tolerance = 1e-15)
    expect_equal(res$estimate, (res$r - m) / (1 - m), tolerance = 1e-15)
  }
  proxy <- common_method_adjusted_r(.cmar_R, N = 120)
  expect_identical(attr(proxy, "marker_rule"), "smallest_positive")
  expect_identical(unique(proxy$marker_r),
                   common_method_marker(.cmar_R)$value[1L])
  expect_identical(unique(proxy$marker_r), 0.20)
  # The pair that supplied the proxy is adjusted to exactly zero.
  supplier <- proxy$variable_1 == "b" & proxy$variable_2 == "d"
  expect_identical(proxy$estimate[supplier], 0)
  expect_identical(proxy$t_value[supplier], 0)
  expect_identical(proxy$p_value[supplier], 1)
})

test_that("with no adjustment the interval is ci_r()'s and the test is on N - 3", {
  N <- 120
  res <- common_method_adjusted_r(.cmar_R, N = N, marker_r = 0)
  for (k in seq_len(nrow(res))) {
    ci <- ci_r(res$r[k], N = N)
    expect_equal(res$ci_lower[k], ci$value[ci$term == "ci_lower"],
                 tolerance = 1e-12)
    expect_equal(res$ci_upper[k], ci$value[ci$term == "ci_upper"],
                 tolerance = 1e-12)
  }
  expect_identical(res$estimate, res$r)
  t_value <- res$r / sqrt((1 - res$r^2) / (N - 3))
  expect_equal(res$t_value, t_value, tolerance = 1e-14)
  expect_equal(res$p_value, 2 * stats::pt(-abs(t_value), N - 3),
               tolerance = 1e-14)
})

test_that("each limit is the interval for rho carried through the adjustment", {
  skip_on_cran()
  N <- 75
  for (m in c(0.05, 0.15, 0.30)) {
    res <- common_method_adjusted_r(.cmar_R, N = N, marker_r = m,
                                    conf_level = 0.90)
    expect_identical(attr(res, "conf_level"), 0.90)
    for (k in seq_len(nrow(res))) {
      ci <- ci_r(res$r[k], N = N, conf_level = 0.90)
      lo <- ci$value[ci$term == "ci_lower"]
      hi <- ci$value[ci$term == "ci_upper"]
      expect_equal(res$ci_lower[k], (lo - m) / (1 - m), tolerance = 1e-12)
      expect_equal(res$ci_upper[k], (hi - m) / (1 - m), tolerance = 1e-12)
      expect_lt(res$ci_lower[k], res$estimate[k])
      expect_gt(res$ci_upper[k], res$estimate[k])
    }
    # The interval excludes zero exactly when the interval for rho
    # excludes the marker correlation.
    rho_lo <- (res$ci_lower * (1 - m)) + m
    expect_identical(res$ci_lower > 0, rho_lo > m)
  }
})

test_that("visual with textual at a marker of .10 reproduces the published check", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939
  scores <- data.frame(visual = rowMeans(hs[c("x1", "x2", "x3")]),
                       textual = rowMeans(hs[c("x4", "x5", "x6")]),
                       speed = rowMeans(hs[c("x7", "x8", "x9")]))
  res <- common_method_adjusted_r(stats::cor(scores), N = nrow(scores),
                                  marker_r = 0.10)
  row <- res[res$variable_1 == "visual" & res$variable_2 == "textual", ]
  expect_equal(row$estimate, 0.237, tolerance = 5e-4 / 0.237)
  expect_equal(row$ci_lower, 0.120, tolerance = 5e-4 / 0.120)
  expect_equal(row$ci_upper, 0.347, tolerance = 5e-4 / 0.347)
  expect_equal(row$t_value, 4.22, tolerance = 5e-3 / 4.22)
  expect_identical(row$df, 298)
})

test_that("the bundled composites give the values on the help page", {
  hs <- holzinger_swineford
  scores <- data.frame(
    visual = rowMeans(hs[c("t1_visual_perception", "t2_cubes",
                           "t4_lozenges")]),
    textual = rowMeans(hs[c("t6_paragraph_comprehension", "t7_sentence",
                            "t9_word_meaning")]),
    speed = rowMeans(hs[c("t10_addition", "t12_counting_groups_of_dots",
                          "t13_straight_and_curved_capitals")]))
  res <- common_method_adjusted_r(stats::cor(scores), N = nrow(scores),
                                  marker_r = c(0, 0.1, 0.2, 0.3))
  vt <- res[res$variable_1 == "visual" & res$variable_2 == "textual", ]
  expect_equal(vt$estimate,
               c(0.31206603, 0.23562892, 0.14008254, 0.01723719),
               tolerance = 1e-7)
  expect_equal(vt$ci_lower,
               c(0.206291652, 0.118101836, 0.007864565, -0.133869068),
               tolerance = 1e-7)
  expect_equal(vt$p_value[3], 0.01517582, tolerance = 1e-6)
})

test_that("an adjusted correlation beyond the range of a correlation gets no inference", {
  R <- matrix(c(1, -0.5, 0.4, -0.5, 1, -0.2, 0.4, -0.2, 1), 3, 3)
  res <- common_method_adjusted_r(R, N = 20, marker_r = 0.3)
  expect_identical(res$variable_1, c("V1", "V1", "V2"))
  expect_identical(res$variable_2, c("V2", "V3", "V3"))
  # -0.5 adjusts to (-0.5 - 0.3) / 0.7, below -1.
  expect_equal(res$estimate[1], -0.8 / 0.7)
  expect_true(all(is.na(unlist(res[1, c("ci_lower", "ci_upper",
                                        "t_value", "p_value")]))))
  # -0.2 adjusts to -5/7, inside; its lower limit is cut at -1.
  expect_equal(res$estimate[3], -0.5 / 0.7)
  expect_identical(res$ci_lower[3], -1)
  expect_gt(res$ci_upper[3], res$estimate[3])
  expect_false(is.na(res$t_value[3]))
})

test_that("a lower limit is cut at -1 exactly when the limit for rho is below 2 r_M - 1", {
  skip_on_cran()  # Extended shape and boundary checks run in the full local suite.
  # The 1.1.0 review: the help page said the cut happens only for a
  # negative correlation, but it depends on the Fisher lower limit for rho
  # alone. A positive correlation with a positive adjusted correlation is
  # cut in a small sample (r = .35, N = 5, r_M = .3 gives .071 with a
  # lower limit of -1), and with r_M of .5 or more a positive lower limit
  # for rho is cut as well.
  one <- function(r, N, m) {
    common_method_adjusted_r(matrix(c(1, r, r, 1), 2), N = N, marker_r = m)
  }
  res <- one(0.35, 5, 0.3)
  expect_equal(res$estimate, 0.05 / 0.7)
  expect_gt(res$estimate, 0)
  expect_identical(res$ci_lower, -1)
  # Over a grid of correlations, sample sizes, and markers, the lower limit
  # is -1 exactly when tanh(atanh(r) - z / sqrt(N - 3)) < 2 r_M - 1.
  grid <- expand.grid(r = c(-0.3, 0.05, 0.35, 0.71, 0.92), N = c(5, 15, 50),
                      m = c(-0.2, 0.1, 0.3, 0.5, 0.7, 0.9))
  grid <- grid[abs((grid$r - grid$m) / (1 - grid$m)) < 1, ]
  crit <- stats::qnorm(0.975)
  for (k in seq_len(nrow(grid))) {
    g <- grid[k, ]
    rho_lower <- tanh(atanh(g$r) - crit / sqrt(g$N - 3))
    lo <- one(g$r, g$N, g$m)$ci_lower
    if (rho_lower < 2 * g$m - 1) {
      expect_identical(lo, -1, info = paste(g, collapse = " "))
    } else {
      expect_equal(lo, (rho_lower - g$m) / (1 - g$m),
                   info = paste(g, collapse = " "))
      expect_gt(lo, -1)
    }
  }
  # A positive lower limit for rho is cut when r_M is above .5.
  rho_lower <- tanh(atanh(0.71) - crit / sqrt(17))
  expect_gt(rho_lower, 0)
  expect_identical(one(0.71, 20, 0.7)$ci_lower, -1)
})

test_that("tidy() and glance() name the pairs and keep the marker values apart", {
  res <- common_method_adjusted_r(.cmar_R, N = 120, marker_r = c(0, 0.1))
  td <- generics::tidy(res)
  expect_identical(names(td), c("term", "marker_r", "r", "estimate",
                                "ci_lower", "ci_upper", "statistic", "df",
                                "p_value", "conf_level"))
  expect_identical(td$term[1:3], c("a:b", "a:c", "a:d"))
  expect_identical(td$marker_r, res$marker_r)
  expect_identical(td$estimate, res$estimate)
  expect_identical(td$statistic, res$t_value)
  expect_true(all(td$conf_level == 0.95))
  gl <- generics::glance(res)
  expect_identical(nrow(gl), 1L)
  expect_identical(gl$n_pairs, 6L)
  expect_identical(gl$n_markers, 2L)
  expect_identical(gl$N, 120)
  expect_identical(gl$conf_level, 0.95)
})

test_that("common_method_adjusted_r() validates its input", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  R <- .cmar_R
  expect_error(common_method_adjusted_r(matrix(2, 2, 2), N = 50),
               "correlation matrix")
  expect_error(common_method_adjusted_r(R * 2, N = 50), "correlation matrix")
  asym <- R
  asym[1, 2] <- 0.1
  expect_error(common_method_adjusted_r(asym, N = 50), "symmetric")
  with_na <- R
  with_na[1, 2] <- with_na[2, 1] <- NA
  expect_error(common_method_adjusted_r(with_na, N = 50), "missing values")
  expect_error(common_method_adjusted_r(as.data.frame(R), N = 50),
               "correlation matrix")
  expect_error(common_method_adjusted_r(R, N = 3), "at least 4")
  expect_error(common_method_adjusted_r(R, N = c(50, 60)), "single")
  expect_error(common_method_adjusted_r(R, N = 50, marker_r = 1),
               "\\(-1, 1\\)")
  expect_error(common_method_adjusted_r(R, N = 50, marker_r = c(0.1, NA)),
               "\\(-1, 1\\)")
  expect_error(common_method_adjusted_r(R, N = 50, marker_r = numeric(0)),
               "\\(-1, 1\\)")
  expect_error(common_method_adjusted_r(R, N = 50, conf_level = 95),
               "conf_level")
  negative <- matrix(c(1, -0.3, -0.3, 1), 2, 2)
  expect_error(common_method_adjusted_r(negative, N = 50),
               "No positive correlation")
})
