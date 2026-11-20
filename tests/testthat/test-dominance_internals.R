# The contract shared by the order-based effect sizes vargha_delaney_A(),
# cliff_delta(), and probability_of_superiority_paired(): what they accept
# as observations, how they compute the pairwise comparisons, and what they
# say when the estimated variance is 0. Everything is tested through the
# three public functions.

# The pairwise formulas the functions used through DMAR 1.1.0's first
# candidate, which build the n_1 by n_2 matrix of comparisons. They are kept
# here as the reference for the rank-based computation.
pairwise_reference <- function(x, y) {
  n_1 <- length(x)
  n_2 <- length(y)
  psi <- outer(x, y, function(a, b) (a > b) + 0.5 * (a == b))
  signs <- sign(outer(x, y, "-"))
  list(A = mean(psi),
       se = sqrt(max(0, stats::var(rowMeans(psi)) / n_1 +
                        stats::var(colMeans(psi)) / n_2)),
       delta = mean(signs),
       var_delta = ((n_2 - 1) * stats::var(rowMeans(signs)) +
                    (n_1 - 1) * stats::var(colMeans(signs)) +
                    stats::var(as.vector(signs))) / (n_1 * n_2),
       p_greater = mean(signs == 1),
       p_less = mean(signs == -1),
       p_tied = mean(signs == 0))
}

cliff_value <- function(res, term) res$value[res$term == term]

expect_matches_reference <- function(x, y, tolerance = 1e-12) {
  ref <- pairwise_reference(x, y)
  a <- suppressWarnings(vargha_delaney_A(x, y))
  d <- suppressWarnings(cliff_delta(x, y))
  expect_equal(a$A, ref$A, tolerance = tolerance)
  expect_equal(a$se, ref$se, tolerance = tolerance)
  expect_equal(cliff_value(d, "cliff_delta"), ref$delta, tolerance = tolerance)
  expect_equal(cliff_value(d, "var_cliff_delta"), ref$var_delta,
               tolerance = tolerance)
  expect_equal(cliff_value(d, "p_y1_greater"), ref$p_greater,
               tolerance = tolerance)
  expect_equal(cliff_value(d, "p_y1_less"), ref$p_less, tolerance = tolerance)
  expect_equal(cliff_value(d, "p_tied"), ref$p_tied, tolerance = tolerance)
}

# Every sample of n_1 and n_2 observations on the points 1, 2, 3, up to the
# order of the observations within a group.
three_point_samples <- function(n_1, n_2) {
  one_group <- function(n) {
    g <- expand.grid(rep(list(1:3), n))
    g <- g[apply(g, 1L, function(r) !is.unsorted(r)), , drop = FALSE]
    lapply(seq_len(nrow(g)), function(i) as.numeric(g[i, ]))
  }
  g_1 <- one_group(n_1)
  g_2 <- one_group(n_2)
  out <- vector("list", length(g_1) * length(g_2))
  k <- 0L
  for (x in g_1) for (y in g_2) {
    k <- k + 1L
    out[[k]] <- list(x = x, y = y)
  }
  out
}

# ---- The rank-based computation ------------------------------------------

test_that("the rank-based estimates and variances equal the pairwise formulas", {
  set.seed(20261109)
  # Continuous, unequal n.
  expect_matches_reference(stats::rnorm(23, 0.4), stats::rnorm(41))
  expect_matches_reference(stats::rt(60, 2), stats::rt(7, 2) + 1)
  # Heavily tied.
  expect_matches_reference(stats::rpois(35, 2), stats::rpois(28, 3))
  expect_matches_reference(sample(1:7, 50, TRUE, prob = (7:1)^2),
                           sample(1:7, 19, TRUE))
  # Binary.
  expect_matches_reference(stats::rbinom(30, 1, 0.3), stats::rbinom(12, 1, 0.6))
  # One constant group, and the smallest design beside a large group.
  expect_matches_reference(rep(2, 6), c(1, 2, 2, 3, 5))
  expect_matches_reference(c(0.2, 1.7), stats::rnorm(50))
  # Values that differ in the last place, and a signed zero.
  expect_matches_reference(c(0.1 + 0.2, 1, -0), c(0.3, 1, 0))
})

test_that("the rank-based computation equals the pairwise formulas on every 3-point sample with n in 2..4", {
  for (n_1 in 2:4) for (n_2 in 2:4) {
    for (s in three_point_samples(n_1, n_2)) {
      expect_matches_reference(s$x, s$y)
    }
  }
})

test_that("infinite values are ordered as the largest and smallest values", {
  with_inf <- cliff_delta(c(1, 2, Inf), c(0, 1, Inf))
  finite   <- cliff_delta(c(1, 2, 10), c(0, 1, 10))
  expect_true(all(is.finite(with_inf$value)))
  expect_equal(with_inf$value, finite$value, tolerance = 1e-12)
  expect_equal(vargha_delaney_A(c(1, 2, Inf), c(0, 1, Inf))$A,
               vargha_delaney_A(c(1, 2, 10), c(0, 1, 10))$A, tolerance = 1e-12)
  expect_equal(vargha_delaney_A(c(-Inf, 2, 3), c(0, 1, 4))$se,
               vargha_delaney_A(c(-9, 2, 3), c(0, 1, 4))$se, tolerance = 1e-12)
  # A pair of equal infinite values is a tie.
  paired_inf <- probability_of_superiority_paired(c(1, 2, Inf, 0),
                                                  c(2, 1, Inf, 3))
  paired_fin <- probability_of_superiority_paired(c(1, 2, 10, 0),
                                                  c(2, 1, 10, 3))
  expect_equal(paired_inf$value, paired_fin$value, tolerance = 1e-12)
})

test_that("20,000 observations per group need no n_1 by n_2 matrix", {
  skip_on_cran()
  set.seed(113)
  x <- round(stats::rnorm(20000, 0.2), 1)
  y <- round(stats::rnorm(20000), 1)
  a <- vargha_delaney_A(x, y)
  d <- cliff_delta(x, y)
  expect_equal(a$A, as.numeric(stats::wilcox.test(x, y)$statistic) / 20000^2,
               tolerance = 1e-12)
  expect_equal(cliff_value(d, "cliff_delta"), 2 * a$A - 1, tolerance = 1e-12)
  expect_equal(sum(d$value[d$term %in% c("p_y1_greater", "p_y1_less",
                                         "p_tied")]), 1, tolerance = 1e-12)
  # The counts of pairs pass the largest integer (4e8 pairs here do not,
  # 60,000 by 40,000 do) without overflow.
  big <- cliff_delta(rep(c(1, 2, 3), 20000), rep(c(1, 2), 20000))
  expect_equal(big$value[big$term %in% c("p_y1_greater", "p_y1_less",
                                         "p_tied")],
               c(1 / 2, 1 / 6, 1 / 3), tolerance = 1e-12)
})

# ---- The bias of the DeLong variance estimator ---------------------------

test_that("with 3 and 3 continuous observations from one distribution, Var(A) is 7/108 and the estimator averages 8/108", {
  # Under one continuous distribution every allocation of the ranks 1..6 to
  # the two groups is equally likely.
  alloc <- utils::combn(6, 3)
  res <- lapply(seq_len(ncol(alloc)), function(k) {
    suppressWarnings(vargha_delaney_A(alloc[, k], setdiff(1:6, alloc[, k])))
  })
  A  <- vapply(res, function(r) r$A, numeric(1))
  se <- vapply(res, function(r) r$se, numeric(1))
  expect_length(A, 20L)
  expect_equal(mean(A), 0.5, tolerance = 1e-12)
  expect_equal(mean((A - 0.5)^2), 7 / 108, tolerance = 1e-12)
  expect_equal(mean(se^2), 8 / 108, tolerance = 1e-12)
})

# ---- What is accepted as observations ------------------------------------

test_that("input with no common order is refused with a classed condition", {
  refused <- function(expr) {
    expect_error(expr, class = "dmar_scores_not_ordered")
  }
  lev <- c("low", "middle", "high")
  o <- factor(c("low", "high", "middle", "high", "low"), levels = lev,
              ordered = TRUE)
  # Labels that are numbers: the codes of "10", "9", "8" are 1, 3, 2.
  refused(vargha_delaney_A(factor(c("10", "9", "8")), c(10, 9, 8)))
  refused(vargha_delaney_A(c("1", "2", "3"), c(1, 2, 3)))
  refused(vargha_delaney_A(list(1, 2, 3), c(1, 2, 3)))
  refused(vargha_delaney_A(data.frame(a = 1:3), c(1, 2, 3)))
  refused(vargha_delaney_A(matrix(1:6, 2), c(1, 2, 3)))
  refused(vargha_delaney_A(c(1i, 2, 3), c(1, 2, 3)))
  refused(vargha_delaney_A(o, c(1, 2, 3)))
  refused(vargha_delaney_A(o, factor(o, levels = rev(lev), ordered = TRUE)))
  # After droplevels() the same label has another code.
  refused(vargha_delaney_A(droplevels(o[o != "low"]), o))
  # 1, 2, 3 hours against 30, 90, 150 minutes.
  refused(vargha_delaney_A(as.difftime(c(1, 2, 3), units = "hours"),
                           as.difftime(c(30, 90, 150), units = "mins")))
  expect_error(vargha_delaney_A(as.difftime(c(1, 2, 3), units = "hours"),
                                c(1, 2, 3)),
               "as.numeric(x, units = ", fixed = TRUE)
  # Days against seconds.
  day <- as.Date("2026-01-01") + c(0, 5, 9, 12)
  refused(vargha_delaney_A(day, as.POSIXct(day)))

  # The formula interface: an unordered factor, a character, and a matrix
  # outcome (which a row index was recycled over).
  d <- data.frame(g = rep(c("a", "b"), each = 5),
                  s = rep(c("agree", "neutral", "disagree", "agree",
                            "neutral"), 2),
                  t = 1:10, stringsAsFactors = FALSE)
  d$f <- factor(d$s)
  refused(vargha_delaney_A(f ~ g, data = d))
  refused(vargha_delaney_A(s ~ g, data = d))
  refused(vargha_delaney_A(cbind(t, t) ~ g, data = d))

  # The other two functions share the contract.
  refused(cliff_delta(factor(c("10", "9", "8")), c(10, 9, 8)))
  refused(cliff_delta(c("1", "2", "3"), c(1, 2, 3)))
  refused(cliff_delta(matrix(1:6, 2), c(1, 2, 3)))
  refused(cliff_delta(o, c(1, 2, 3)))
  refused(probability_of_superiority_paired(factor(c("10", "9", "8")),
                                            c(10, 9, 8)))
  refused(probability_of_superiority_paired(c("1", "2", "3"), c(1, 2, 3)))
  refused(probability_of_superiority_paired(
    as.difftime(c(1, 2, 3), units = "hours"),
    as.difftime(c(30, 90, 150), units = "mins")))

  # The refusal is a dmar_refusal with its code.
  e <- tryCatch(cliff_delta(factor(1:3), 1:3), error = identity)
  expect_s3_class(e, "dmar_refusal")
  expect_identical(e$code, "scores_not_ordered")
})

test_that("ordered factors with the same levels are read through their order", {
  set.seed(4)
  lev <- c("strongly disagree", "disagree", "neutral", "agree",
           "strongly agree")
  a <- sample(lev, 30, TRUE, prob = c(.30, .30, .20, .10, .10))
  b <- sample(lev, 30, TRUE, prob = c(.10, .20, .20, .25, .25))
  o_a <- factor(a, levels = lev, ordered = TRUE)
  o_b <- factor(b, levels = lev, ordered = TRUE)
  from_factor <- vargha_delaney_A(o_a, o_b)
  from_codes  <- vargha_delaney_A(match(a, lev), match(b, lev))
  rescored    <- vargha_delaney_A(c(1, 10, 100, 1000, 10000)[match(a, lev)],
                                  c(1, 10, 100, 1000, 10000)[match(b, lev)])
  expect_equal(from_factor, from_codes, tolerance = 1e-12)
  expect_equal(from_factor, rescored, tolerance = 1e-12)
  # The alphabetical codes of the labels give another value.
  alphabetical <- vargha_delaney_A(as.integer(factor(a)), as.integer(factor(b)))
  expect_gt(abs(alphabetical$A - from_factor$A), 0.05)

  d <- data.frame(s = factor(c(a, b), levels = lev, ordered = TRUE),
                  g = rep(c("first", "second"), each = 30))
  expect_equal(vargha_delaney_A(s ~ g, data = d), from_factor,
               tolerance = 1e-12)

  expect_equal(cliff_delta(o_a, o_b)$value,
               cliff_delta(match(a, lev), match(b, lev))$value,
               tolerance = 1e-12)
  expect_equal(probability_of_superiority_paired(o_a, o_b)$value,
               probability_of_superiority_paired(match(a, lev),
                                                 match(b, lev))$value,
               tolerance = 1e-12)
})

test_that("logical vectors, one-column matrices, 1-d arrays, and dates are read as vectors", {
  x <- c(1, 2, 5, 7)
  y <- c(1, 3, 4)
  plain <- vargha_delaney_A(x, y)
  expect_equal(vargha_delaney_A(matrix(x), y), plain, tolerance = 1e-12)
  expect_equal(vargha_delaney_A(tapply(x, 1:4, sum), matrix(y)), plain,
               tolerance = 1e-12)
  expect_equal(cliff_delta(matrix(x), y)$value, cliff_delta(x, y)$value,
               tolerance = 1e-12)

  l_1 <- c(TRUE, FALSE, TRUE, TRUE)
  l_2 <- c(FALSE, FALSE, TRUE, FALSE)
  expect_equal(vargha_delaney_A(l_1, l_2),
               vargha_delaney_A(as.numeric(l_1), as.numeric(l_2)),
               tolerance = 1e-12)

  # Dates and date-times keep the answers they had.
  day_1 <- as.Date("2026-01-01") + c(0, 5, 9, 12)
  day_2 <- day_1 + c(1, -2, 3, 7)
  from_days <- vargha_delaney_A(as.numeric(day_1), as.numeric(day_2))
  expect_equal(vargha_delaney_A(day_1, day_2), from_days, tolerance = 1e-12)
  expect_equal(vargha_delaney_A(as.POSIXct(day_1), as.POSIXct(day_2)),
               from_days, tolerance = 1e-12)
})

test_that("cliff_delta() is 2A - 1 for each accepted kind of input", {
  lev <- c("low", "middle", "high")
  pairs <- list(
    numeric = list(c(1, 2, 2, 5, 7), c(1, 3, 4, 4)),
    integer = list(c(1L, 4L, 4L, 6L), c(2L, 4L, 5L)),
    logical = list(c(TRUE, FALSE, TRUE, TRUE), c(FALSE, FALSE, TRUE, FALSE)),
    ordered = list(factor(c("low", "high", "middle", "high"), levels = lev,
                          ordered = TRUE),
                   factor(c("low", "low", "middle", "high", "middle"),
                          levels = lev, ordered = TRUE)),
    infinite = list(c(1, 2, Inf), c(0, 1, Inf)),
    matrix = list(matrix(c(1, 2, 5, 7)), c(1, 3, 4)))
  for (p in pairs) {
    expect_equal(cliff_delta(p[[1L]], p[[2L]])$value[1L],
                 2 * vargha_delaney_A(p[[1L]], p[[2L]])$A - 1,
                 tolerance = 1e-12)
  }
})

test_that("missing values are removed and counted in n_missing", {
  complete <- vargha_delaney_A(c(1, 2, 5, 7), c(1, 3, 4))
  with_na  <- vargha_delaney_A(c(1, 2, NA, 5, 7), c(1, 3, 4))
  expect_identical(attr(with_na, "n_missing"), c(group_1 = 1L, group_2 = 0L))
  expect_identical(attr(complete, "n_missing"), c(group_1 = 0L, group_2 = 0L))
  expect_equal(with_na, complete, ignore_attr = "n_missing", tolerance = 1e-12)

  # The formula interface counts the missing outcomes of each group. A row
  # with a missing group label belongs to neither group and is not counted,
  # and a level with no observation is not a group.
  d <- data.frame(y = c(1, 2, NA, 5, 7, 1, 3, 4, 9),
                  g = factor(c(rep("a", 5), rep("b", 3), NA),
                             levels = c("a", "b", "c")))
  from_formula <- vargha_delaney_A(y ~ g, data = d)
  expect_identical(attr(from_formula, "n_missing"),
                   c(group_1 = 1L, group_2 = 0L))
  expect_equal(from_formula, complete, ignore_attr = "n_missing",
               tolerance = 1e-12)
  expect_error(vargha_delaney_A(y ~ g, data = d[d$g %in% "a", ]),
               "exactly two observed levels")

  cd <- cliff_delta(c(1, 2, NA), c(1, 2, 3))
  expect_identical(attr(cd, "n_missing"), c(group_1 = 1L, group_2 = 0L))
  expect_equal(cd$value, cliff_delta(c(1, 2), c(1, 2, 3))$value,
               tolerance = 1e-12)
  expect_error(cliff_delta(c(1, NA), c(1, 2, 3)),
               "at least 2 non-missing observations")

  # A pair with a missing member leaves as a pair.
  ps <- probability_of_superiority_paired(c(1, NA, 3, 5, 2), c(0, 1, 2, 6, NaN))
  expect_identical(attr(ps, "n_missing"), 2L)
  expect_equal(ps$value,
               probability_of_superiority_paired(c(1, 3, 5), c(0, 2, 6))$value,
               tolerance = 1e-12)
  expect_error(probability_of_superiority_paired(c(1, NA, 3), c(0, 1, 2)),
               "at least 3 complete pairs")
})

test_that("a missing or vector conf_level gets the worded error in all three functions", {
  msg <- "'conf_level' must be a single number in (0, 1)."
  for (bad in list(NA_real_, c(.90, .95), "0.95", 1, 0)) {
    expect_error(vargha_delaney_A(1:4, 2:5, conf_level = bad), msg,
                 fixed = TRUE)
    expect_error(cliff_delta(1:4, 2:5, conf_level = bad), msg, fixed = TRUE)
    expect_error(probability_of_superiority_paired(1:4, c(2, 1, 5, 6),
                                                   conf_level = bad),
                 msg, fixed = TRUE)
  }
})

# ---- An estimated variance of 0 ------------------------------------------

# The classed warning of a call, or NULL when the call gives none.
degenerate_warning <- function(expr) {
  tryCatch({
    expr
    NULL
  }, dmar_interval_degenerate = function(w) w)
}

test_that("complete separation gives degenerate limits and the classed warning in all three functions", {
  for (method in c("logit", "wald")) {
    expect_warning(high <- vargha_delaney_A(c(4, 5, 6), c(1, 2, 3),
                                            ci_method = method),
                   class = "dmar_interval_degenerate")
    expect_warning(low <- vargha_delaney_A(c(1, 2, 3), c(4, 5, 6),
                                           ci_method = method),
                   class = "dmar_interval_degenerate")
    expect_identical(c(high$A, high$se, high$ci_lower, high$ci_upper),
                     c(1, 0, 1, 1))
    expect_identical(c(low$A, low$se, low$ci_lower, low$ci_upper),
                     c(0, 0, 0, 0))
    expect_true(is.na(high$test_statistic) && is.na(high$p_value))
    expect_identical(high$ci_method, method)
  }
  w <- degenerate_warning(vargha_delaney_A(c(4, 5, 6), c(1, 2, 3)))
  expect_identical(w$reason, "complete_separation")
  expect_match(conditionMessage(w), "completely separated (A = 1)",
               fixed = TRUE)
  expect_match(conditionMessage(w),
               "not a confidence interval with the stated coverage",
               fixed = TRUE)

  # cliff_delta() keeps its 1.0.x limits and wording, now with a class.
  expect_warning(cd <- cliff_delta(c(10, 11, 12), c(1, 2, 3)),
                 class = "dmar_interval_degenerate")
  expect_identical(cd$value[1:4], c(1, 1, 1, 0))
  w <- degenerate_warning(cliff_delta(c(1, 2, 3), c(10, 11, 12)))
  expect_identical(w$reason, "complete_separation")
  # The opening words are those of 1.0.x through "completely separated";
  # the rest of the message was rewritten in 1.1.0.
  expect_match(conditionMessage(w),
               "cliff_delta(): the groups are completely separated (|delta| = 1), so the estimated variance is 0. Both confidence limits are reported equal to the estimate (-1, -1).",
               fixed = TRUE)
  expect_match(conditionMessage(w),
               "not a confidence interval with the stated coverage",
               fixed = TRUE)

  # The paired function returned NaN limits here, with no warning.
  expect_warning(ps <- probability_of_superiority_paired(c(1, 2, 3, 4),
                                                         c(4, 5, 6, 7)),
                 class = "dmar_interval_degenerate")
  expect_identical(ps$value[1:4], c(1, 1, 1, 0))
  expect_warning(ps <- probability_of_superiority_paired(c(4, 5, 6, 7),
                                                         c(1, 2, 3, 4)),
                 class = "dmar_interval_degenerate")
  expect_identical(ps$value[1:4], c(0, 0, 0, 0))
  w <- degenerate_warning(probability_of_superiority_paired(c(1, 2, 3, 4),
                                                            c(4, 5, 6, 7)))
  expect_identical(w$reason, "complete_separation")
  expect_match(conditionMessage(w), "every pair favors the same member",
               fixed = TRUE)
})

test_that("observations that are all equal give the warning with reason no_variation", {
  w <- degenerate_warning(a <- vargha_delaney_A(c(2, 2, 2), c(2, 2, 2)))
  expect_identical(w$reason, "no_variation")
  a <- suppressWarnings(vargha_delaney_A(c(2, 2, 2), c(2, 2, 2)))
  expect_identical(c(a$A, a$se, a$ci_lower, a$ci_upper),
                   c(0.5, 0, 0.5, 0.5))

  w <- degenerate_warning(cliff_delta(c(2, 2, 2), c(2, 2, 2)))
  expect_identical(w$reason, "no_variation")
  cd <- suppressWarnings(cliff_delta(c(2, 2, 2), c(2, 2, 2)))
  expect_identical(cd$value[1:4], c(0, 0, 0, 0))

  w <- degenerate_warning(probability_of_superiority_paired(c(2, 2, 2, 2),
                                                            c(2, 2, 2, 2)))
  expect_identical(w$reason, "no_variation")
  ps <- suppressWarnings(probability_of_superiority_paired(c(2, 2, 2, 2),
                                                           c(2, 2, 2, 2)))
  expect_identical(ps$value[1:4], c(0.5, 0.5, 0.5, 0))
})

test_that("groups that touch or overlap give no warning and an interval", {
  expect_no_warning(a <- vargha_delaney_A(c(2, 2, 2), c(1, 2, 3)))
  expect_gt(a$ci_upper - a$ci_lower, 0)
  expect_no_warning(a <- vargha_delaney_A(c(3, 4, 5), c(1, 2, 3)))
  expect_lt(a$ci_upper, 1)
  expect_no_warning(cd <- cliff_delta(c(3, 4, 5), c(1, 2, 3)))
  expect_lt(cd$value[3L], 1)
  expect_no_warning(ps <- probability_of_superiority_paired(c(1, 2, 3, 4),
                                                            c(2, 3, 4, 4)))
  expect_lt(ps$value[2L], ps$value[3L])
})

test_that("on every 3-point sample with n in 2..4 the warning is given exactly when the variance is 0", {
  for (n_1 in 2:4) for (n_2 in 2:4) {
    for (s in three_point_samples(n_1, n_2)) {
      separated <- min(s$x) > max(s$y) || max(s$x) < min(s$y)
      flat <- min(c(s$x, s$y)) == max(c(s$x, s$y))
      w_a <- degenerate_warning(vargha_delaney_A(s$x, s$y))
      w_d <- degenerate_warning(cliff_delta(s$x, s$y))
      a <- suppressWarnings(vargha_delaney_A(s$x, s$y))
      d <- suppressWarnings(cliff_delta(s$x, s$y))
      expect_identical(!is.null(w_a), separated || flat)
      expect_identical(!is.null(w_d), separated || flat)
      expect_identical(a$se == 0, separated || flat)
      expect_identical(cliff_value(d, "var_cliff_delta") == 0,
                       separated || flat)
      expect_identical(a$A %in% c(0, 1), separated)
      if (separated || flat) {
        expect_identical(w_a$reason,
                         if (flat) "no_variation" else "complete_separation")
        expect_identical(w_d$reason, w_a$reason)
        expect_identical(c(a$ci_lower, a$ci_upper), c(a$A, a$A))
      }
    }
  }
})

test_that("randomization_test() passes on the one classed warning of separated groups", {
  d <- data.frame(y = c(1, 2, 3, 4, 11, 12, 13, 14),
                  g = rep(c("a", "b"), each = 4))
  seen <- new.env(parent = emptyenv())
  seen$n <- 0L
  res <- withCallingHandlers(
    randomization_test(y ~ g, data = d),
    dmar_interval_degenerate = function(w) {
      seen$n <- seen$n + 1L
      invokeRestart("muffleWarning")
    })
  expect_identical(seen$n, 1L)
  expect_identical(res$value[res$term %in% c("cliff_delta",
                                             "cliff_delta_ci_lower",
                                             "cliff_delta_ci_upper")],
                   c(-1, -1, -1))
})

test_that("time differences in the same units are ordered, and in different units refused", {
  h1 <- as.difftime(c(1, 2, 3), units = "hours")
  h2 <- as.difftime(c(0.5, 1.5, 2.5), units = "hours")
  expect_equal(vargha_delaney_A(h1, h2),
               vargha_delaney_A(c(1, 2, 3), c(0.5, 1.5, 2.5)))
  expect_equal(vargha_delaney_A(h1, h2)$A, 2 / 3)
  expect_equal(cliff_delta(h1, h2),
               cliff_delta(c(1, 2, 3), c(0.5, 1.5, 2.5)))
  expect_equal(suppressWarnings(probability_of_superiority_paired(h1, h2)),
               suppressWarnings(
                 probability_of_superiority_paired(c(1, 2, 3),
                                                   c(0.5, 1.5, 2.5))))
  # One column has one unit, so the formula interface accepts it.
  d <- data.frame(s = as.difftime(c(3, 1, 4, 1, 5, 9, 2, 6, 2, 7, 1, 8, 2, 8),
                                  units = "mins"),
                  g = rep(c("b", "a"), c(8, 6)))
  d$n <- as.numeric(d$s)
  res <- vargha_delaney_A(s ~ g, data = d)
  expect_equal(res, vargha_delaney_A(n ~ g, data = d))
  expect_equal(res$A, 0.5625)
  # Different units, or a time difference beside plain numbers.
  expect_error(vargha_delaney_A(h1, as.difftime(c(30, 90, 150),
                                                units = "mins")),
               class = "dmar_scores_not_ordered")
  expect_error(cliff_delta(h1, c(1, 2, 3)),
               class = "dmar_scores_not_ordered")
})

test_that("a level of the grouping variable with no outcome is not a group", {
  # Level "c" has only missing outcomes, so the two groups are "a" and "b",
  # whether the grouping variable is a character vector or a factor.
  d <- data.frame(s = c(1, 4, 2, 8, 5, 7, NA, NA),
                  g = c("a", "a", "a", "b", "b", "b", "c", "c"),
                  stringsAsFactors = FALSE)
  two <- suppressWarnings(vargha_delaney_A(c(1, 4, 2), c(8, 5, 7)))
  for (grouping in list(d$g, factor(d$g))) {
    d$grp <- grouping
    res <- suppressWarnings(vargha_delaney_A(s ~ grp, data = d))
    expect_identical(c(res$n_1, res$n_2), c(3L, 3L))
    expect_equal(as.data.frame(res), as.data.frame(two),
                 ignore_attr = TRUE)
  }
})

test_that("a POSIXlt vector is read as the date-time it holds", {
  t1 <- as.POSIXct("2026-03-01 12:00:00", tz = "UTC") +
    3600 * c(0, 5, 9, 2)
  t2 <- as.POSIXct("2026-03-01 12:00:00", tz = "UTC") +
    3600 * c(1, 5, 20, 3)
  expect_equal(vargha_delaney_A(as.POSIXlt(t1), t2),
               vargha_delaney_A(t1, t2))
  expect_equal(cliff_delta(as.POSIXlt(t1), t2), cliff_delta(t1, t2))
})

test_that("an ordered factor beside a vector that is not one is refused in those words", {
  lev <- c("low", "mid", "high")
  o <- factor(c("low", "mid", "high", "mid"), levels = lev, ordered = TRUE)
  expect_error(vargha_delaney_A(o, c(1, 2, 3)),
               "'x' is an ordered factor and 'y' is a numeric vector",
               fixed = TRUE)
})

# ---- What each refusal says (1.1.0) --------------------------------------

test_that("each refusal says why, what would have been computed, and what to supply", {
  refusal <- function(expr) {
    e <- tryCatch(expr, dmar_scores_not_ordered = function(e) e)
    expect_s3_class(e, "dmar_scores_not_ordered")
    expect_s3_class(e, "dmar_refusal")
    expect_match(conditionMessage(e), "Nothing was computed.", fixed = TRUE)
    e
  }
  lev <- c("low", "mid", "high")
  f1 <- factor(c("low", "mid", "high", "high"))
  f2 <- factor(c("low", "low", "mid", "high"))

  # An unordered factor: the levels and codes of the object itself.
  e <- refusal(vargha_delaney_A(f1, f2))
  expect_identical(e$arg, "x")
  expect_match(conditionMessage(e), "labels with no order", fixed = TRUE)
  expect_match(conditionMessage(e),
               "here the levels are high, low, mid, stored as 1, 2, 3",
               fixed = TRUE)
  expect_match(conditionMessage(e),
               "would rank high below low below mid", fixed = TRUE)
  expect_match(conditionMessage(e),
               "factor(x, levels = c(...), ordered = TRUE)", fixed = TRUE)
  # The codes give A = 0.40625 and the ordered levels 0.65625, as the help
  # says.
  expect_equal(vargha_delaney_A(as.integer(f1), as.integer(f2))$A, 0.40625)
  expect_equal(vargha_delaney_A(factor(f1, levels = lev, ordered = TRUE),
                                factor(f2, levels = lev, ordered = TRUE))$A,
               0.65625)
  # The argument at fault is named, in each function and on the formula
  # path, and a long list of levels is cut at six.
  e <- refusal(cliff_delta(c(1, 2, 3), f2[1:3]))
  expect_identical(e$arg, "group_2")
  expect_match(conditionMessage(e),
               "cliff_delta() cannot use an unordered factor ('group_2')",
               fixed = TRUE)
  d <- data.frame(g = rep(c("a", "b"), each = 4), f = c(f1, f2))
  e <- refusal(vargha_delaney_A(f ~ g, data = d))
  expect_identical(e$arg, "f")
  e <- refusal(probability_of_superiority_paired(factor(letters[1:9]), 1:9))
  expect_match(conditionMessage(e),
               "the levels are a, b, c, d, e, f, and so on, stored as 1, 2, 3, 4, 5, 6, and so on. Computed",
               fixed = TRUE)
  expect_no_match(conditionMessage(e), "f, ...", fixed = TRUE)
  expect_no_match(conditionMessage(e), "6, ...", fixed = TRUE)
  # The paired function speaks of the members of a pair, not of groups.
  expect_match(conditionMessage(e),
               "supply both members of the pairs as ordered factors",
               fixed = TRUE)
  expect_no_match(conditionMessage(e), "groups", fixed = TRUE)
  e <- refusal(probability_of_superiority_paired(
    factor(c("lo", "hi", "lo"), levels = c("lo", "hi"), ordered = TRUE),
    c(1, 2, 3)))
  expect_match(conditionMessage(e),
               "probability_of_superiority_paired() needs both members of the pairs on one scale",
               fixed = TRUE)
  expect_no_match(conditionMessage(e), "group", fixed = TRUE)
  e <- refusal(vargha_delaney_A(
    factor(c("lo", "hi", "lo"), levels = c("lo", "hi"), ordered = TRUE),
    c(1, 2, 3)))
  expect_match(conditionMessage(e), "needs both groups on one scale",
               fixed = TRUE)
  # A factor with no levels has none to list.
  e <- refusal(vargha_delaney_A(factor(c(NA, NA)), 1:2))
  expect_match(conditionMessage(e),
               "which is alphabetical unless it was set. Computed from those codes, the statistic would change",
               fixed = TRUE)
  expect_no_match(conditionMessage(e), "here the levels are", fixed = TRUE)

  # A character vector has its own message.
  e <- refusal(vargha_delaney_A(c(1, 2, 3), c("1", "2", "3")))
  expect_identical(e$arg, "y")
  expect_match(conditionMessage(e), "alphabetical order, which is rarely",
               fixed = TRUE)
  expect_match(conditionMessage(e), "as.numeric(y)", fixed = TRUE)
  e <- refusal(cliff_delta(c("a", "b"), c(1, 2)))
  expect_match(conditionMessage(e),
               "cliff_delta() cannot use a character vector ('group_1')",
               fixed = TRUE)

  # An ordered factor beside anything else.
  o <- factor(c("low", "mid", "high", "mid"), levels = lev, ordered = TRUE)
  e <- refusal(vargha_delaney_A(c(TRUE, FALSE, TRUE), o))
  expect_match(conditionMessage(e), "positions in a list of levels",
               fixed = TRUE)
  expect_match(conditionMessage(e),
               "'y' is an ordered factor and 'x' is a logical vector",
               fixed = TRUE)
  e <- refusal(cliff_delta(o, as.Date("2026-01-01") + 0:3))
  expect_match(conditionMessage(e), "is a vector of dates (class Date)",
               fixed = TRUE)

  # Two ordered factors with different levels: a shared label with two
  # codes, no shared label, and shared labels with the same codes.
  e <- refusal(vargha_delaney_A(o, factor(c("low", "high"),
                                          levels = c("low", "high"),
                                          ordered = TRUE)))
  expect_match(conditionMessage(e), "numbers its own levels", fixed = TRUE)
  expect_match(conditionMessage(e),
               "'x' has low < mid < high and 'y' has low < high",
               fixed = TRUE)
  expect_match(conditionMessage(e),
               "(high is 3 in one and 2 in the other)", fixed = TRUE)
  e <- refusal(vargha_delaney_A(
    factor(c("a", "b"), levels = c("a", "b"), ordered = TRUE),
    factor(c("c", "d"), levels = c("c", "d"), ordered = TRUE)))
  expect_match(conditionMessage(e), "numbers its own levels", fixed = TRUE)
  expect_match(conditionMessage(e), "(1 is a in one and c in the other)",
               fixed = TRUE)
  e <- refusal(vargha_delaney_A(o, droplevels(o[o != "high"])))
  expect_match(conditionMessage(e), "numbers its own levels", fixed = TRUE)
  expect_match(conditionMessage(e),
               "the labels the two share have the same codes here",
               fixed = TRUE)

  # Other classes, shapes, time differences, and a Date beside a POSIXct
  # keep their sentences and gain the closing one.
  e <- refusal(vargha_delaney_A(c(1i, 2, 3), c(1, 2, 3)))
  expect_match(conditionMessage(e),
               "an object of this class has no order the function can read",
               fixed = TRUE)
  e <- refusal(vargha_delaney_A(matrix(1:6, 2), c(1, 2, 3)))
  expect_match(conditionMessage(e), "it is a matrix with 3 columns",
               fixed = TRUE)
  refusal(vargha_delaney_A(as.difftime(c(1, 2, 3), units = "hours"),
                           as.difftime(c(30, 90, 150), units = "mins")))
  day <- as.Date("2026-01-01") + c(0, 5, 9, 12)
  refusal(vargha_delaney_A(day, as.POSIXct(day)))
})

# ---- What the zero-variance warning carries (1.1.0) ----------------------

test_that("the warning carries the estimate, the sample sizes, and the Feng and Cliff limit", {
  # Feng and Cliff (2004, Equation 8) at the 95 percent level:
  # (n_b - z^2) / (n_b + z^2) for n_b = 3, 4, 5, 10, 20.
  published <- c(`3` = -0.123, `4` = 0.020, `5` = 0.131, `10` = 0.445,
                 `20` = 0.678)
  for (n_b in c(3, 4, 5, 10, 20)) {
    hi <- seq_len(n_b) + 100
    lo <- seq_len(n_b + 2)
    w <- degenerate_warning(cliff_delta(hi, lo))
    expect_identical(w$reason, "complete_separation")
    expect_identical(w$estimate, 1)
    expect_identical(w$n, c(as.integer(n_b), as.integer(n_b + 2)))
    expect_equal(round(w$limit_feng_cliff, 3), published[[as.character(n_b)]])
    # Its negative when delta-hat is -1, with the smaller group second.
    w <- degenerate_warning(cliff_delta(lo, hi))
    expect_identical(w$estimate, -1)
    expect_equal(round(w$limit_feng_cliff, 3),
                 -published[[as.character(n_b)]])
  }
  z2 <- stats::qnorm(0.975)^2
  w <- degenerate_warning(cliff_delta(6:10, 1:5))
  expect_match(conditionMessage(w),
               "give the limit (n_b - z^2) / (n_b + z^2) for delta, with n_b the smaller group size and z the normal quantile of the level: 0.131 here at the 95% level, with the other limit 1.",
               fixed = TRUE)
  expect_match(conditionMessage(w), "(here 5 and 5 observations)",
               fixed = TRUE)
  expect_match(conditionMessage(w),
               "not a confidence interval with the stated coverage",
               fixed = TRUE)
  expect_match(conditionMessage(w), "randomization_test()", fixed = TRUE)
  # The level of the call sets z.
  w <- degenerate_warning(cliff_delta(6:10, 1:5, conf_level = 0.90))
  z2_90 <- stats::qnorm(0.95)^2
  expect_equal(w$limit_feng_cliff, (5 - z2_90) / (5 + z2_90),
               tolerance = 1e-12)
  expect_match(conditionMessage(w), "at the 90% level", fixed = TRUE)

  # On the scale of A: n_b / (n_b + z^2) when A is 1, z^2 / (n_b + z^2)
  # when A is 0; 0.566 with five in the smaller group.
  w <- degenerate_warning(vargha_delaney_A(6:10, 1:5))
  expect_identical(w$estimate, 1)
  expect_identical(w$n, c(5L, 5L))
  expect_equal(w$limit_feng_cliff, 5 / (5 + z2), tolerance = 1e-12)
  expect_equal(round(w$limit_feng_cliff, 3), 0.566)
  expect_match(conditionMessage(w), "on the scale of A, 0.566", fixed = TRUE)
  w <- degenerate_warning(vargha_delaney_A(1:5, 6:12))
  expect_identical(w$estimate, 0)
  expect_equal(w$limit_feng_cliff, z2 / (5 + z2), tolerance = 1e-12)

  # The paired function, and no variation: no such limit.
  w <- degenerate_warning(probability_of_superiority_paired(c(1, 2, 3, 4),
                                                            c(4, 5, 6, 7)))
  expect_identical(w$estimate, 1)
  expect_identical(w$n, 4L)
  expect_true(is.na(w$limit_feng_cliff))
  expect_match(conditionMessage(w), "(here 4 pairs)", fixed = TRUE)
  expect_match(conditionMessage(w), "randomization_test_paired()",
               fixed = TRUE)
  expect_false(grepl("Feng and Cliff", conditionMessage(w), fixed = TRUE))
  for (w in list(
    degenerate_warning(vargha_delaney_A(c(2, 2, 2), c(2, 2, 2))),
    degenerate_warning(cliff_delta(c(2, 2, 2), c(2, 2, 2))),
    degenerate_warning(probability_of_superiority_paired(c(2, 2, 2, 2),
                                                         c(2, 2, 2, 2))))) {
    expect_identical(w$reason, "no_variation")
    expect_true(is.na(w$limit_feng_cliff))
    expect_match(conditionMessage(w),
                 "not a confidence interval with the stated coverage",
                 fixed = TRUE)
    expect_match(conditionMessage(w), "Check that the scores were read",
                 fixed = TRUE)
  }
})

test_that("infinite values, constant groups, and groups that touch behave as the help says", {
  # A group of Inf beside finite scores is completely separated.
  w <- degenerate_warning(a <- vargha_delaney_A(c(Inf, Inf), c(1, 2, 3)))
  expect_identical(w$reason, "complete_separation")
  w <- degenerate_warning(cliff_delta(c(1, 2, 3), c(Inf, Inf)))
  expect_identical(c(w$reason, w$estimate), c("complete_separation", "-1"))
  # The same infinite value throughout is no variation.
  w <- degenerate_warning(vargha_delaney_A(c(Inf, Inf), c(Inf, Inf, Inf)))
  expect_identical(w$reason, "no_variation")
  w <- degenerate_warning(probability_of_superiority_paired(rep(Inf, 3),
                                                            rep(Inf, 3)))
  expect_identical(w$reason, "no_variation")
  # Two constant groups at different values are completely separated.
  w <- degenerate_warning(vargha_delaney_A(c(2, 2), c(1, 1)))
  expect_identical(w$reason, "complete_separation")
  # Groups that touch have an ordinary, wide interval (the help's numbers).
  expect_no_warning(a <- vargha_delaney_A(c(3, 4, 5), c(1, 2, 3)))
  expect_equal(round(c(a$A, a$ci_lower, a$ci_upper), 3),
               c(0.944, 0.475, 0.997))
  expect_no_warning(cd <- cliff_delta(c(3, 4, 5), c(1, 2, 3)))
  expect_equal(round(cd$value[1:3], 3), c(0.889, -0.167, 0.995))
  # Three wins and a tie.
  expect_no_warning(ps <- probability_of_superiority_paired(c(1, 2, 3, 4),
                                                            c(2, 3, 4, 4)))
  expect_equal(round(ps$value[1:3], 3), c(0.875, 0.502, 0.980))
  # Fewer than three complete pairs is an error, not the degenerate case.
  expect_error(probability_of_superiority_paired(c(1, 2), c(2, 3)),
               "at least 3 complete pairs")
})

# ---- The mark a degenerate table carries (1.1.0) -------------------------

test_that("a degenerate table says so in its attribute, sentence, footer, tidy(), and metadata", {
  a  <- suppressWarnings(vargha_delaney_A(c(4, 5, 6), c(1, 2, 3)))
  cd <- suppressWarnings(cliff_delta(group_1 = c(4, 5, 6),
                                     group_2 = c(1, 2, 3)))
  ps <- suppressWarnings(probability_of_superiority_paired(c(1, 2, 3, 4),
                                                           c(4, 5, 6, 7)))
  for (res in list(a, cd, ps)) {
    expect_identical(attr(res, "interval_degenerate"), "complete_separation")
    expect_identical(dmar_metadata(res)$interval_degenerate,
                     "complete_separation")
    expect_true(
      "Limits are not a confidence interval: the estimated variance is 0." %in%
        utils::capture.output(print(res)))
    expect_true(all(generics::tidy(res)$interval_degenerate))
    expect_true(generics::glance(res)$interval_degenerate)
  }
  # The sentences. Through the first versions of 1.1.0 each read
  # "... = 1.00, 95% CI [1.00, 1.00]".
  expect_identical(
    results_sentence(a),
    "A = 1.00 (no confidence interval: the groups are completely separated, so the estimated variance is 0)")
  expect_identical(
    results_sentence(cd),
    "cliff_delta = 1.00 (no confidence interval: the groups are completely separated, so the estimated variance is 0)")
  expect_identical(
    results_sentence(ps),
    "probability_of_superiority = 1.00 (no confidence interval: every pair favors the same member, so the estimated variance is 0)")
  # The label and the digits follow the usual rules.
  expect_identical(
    results_sentence(cd, label = "Cliff's delta", digits = 3),
    "Cliff's delta = 1.000 (no confidence interval: the groups are completely separated, so the estimated variance is 0)")
  # No variation.
  flat_a <- suppressWarnings(vargha_delaney_A(c(2, 2, 2), c(2, 2, 2)))
  expect_identical(attr(flat_a, "interval_degenerate"), "no_variation")
  expect_identical(
    results_sentence(flat_a),
    "A = 0.50 (no confidence interval: every score is equal, so the estimated variance is 0)")
  flat_cd <- suppressWarnings(cliff_delta(group_1 = c(2, 2, 2),
                                          group_2 = c(2, 2, 2)))
  expect_identical(attr(flat_cd, "interval_degenerate"), "no_variation")
  expect_identical(
    results_sentence(flat_cd),
    "cliff_delta = 0.00 (no confidence interval: every score is equal, so the estimated variance is 0)")
  expect_true(all(generics::tidy(flat_cd)$interval_degenerate))
  expect_true(
    "Limits are not a confidence interval: the estimated variance is 0." %in%
      utils::capture.output(print(flat_cd)))
  flat_ps <- suppressWarnings(
    probability_of_superiority_paired(c(2, 2, 2, 2), c(2, 2, 2, 2)))
  expect_identical(
    results_sentence(flat_ps),
    "probability_of_superiority = 0.50 (no confidence interval: every pair is tied, so the estimated variance is 0)")
  # A p-value is appended only when the table has one: the permutation
  # test of completely separated groups.
  perm <- suppressWarnings(vargha_delaney_A(6:10, 1:5, test = "permutation"))
  expect_identical(
    results_sentence(perm),
    "A = 1.00 (no confidence interval: the groups are completely separated, so the estimated variance is 0), unadjusted p = 0.0079")
  # A tidy() or glance() row carries the mark as a column, which says that
  # the variance is 0 and not which case gave it; a stack keeps it by row.
  expect_identical(
    results_sentence(generics::tidy(a)),
    "A = 1.00 (no confidence interval: the estimated variance is 0)")
  expect_identical(
    results_sentence(generics::tidy(cd)),
    "cliff_delta = 1.00 (no confidence interval: the estimated variance is 0)")
  expect_identical(
    results_sentence(generics::glance(ps)),
    "probability_of_superiority = 1.00 (no confidence interval: the estimated variance is 0)")
  ordinary <- vargha_delaney_A(c(1, 2, 5, 7), c(1, 3, 4))
  stack <- rbind(generics::tidy(a), generics::tidy(ordinary))
  expect_identical(stack$interval_degenerate, c(TRUE, FALSE))
  expect_match(results_sentence(stack[2L, ]), "95% CI [", fixed = TRUE)
  expect_match(results_sentence(stack[1L, ]), "no confidence interval",
               fixed = TRUE)

  # An ordinary table has no attribute, no footer line, and its sentence.
  cd_ordinary <- cliff_delta(c(1, 2, 5, 7), c(1, 3, 4))
  ps_ordinary <- probability_of_superiority_paired(c(1, 2, 3, 4),
                                                   c(2, 1, 5, 6))
  for (res in list(ordinary, cd_ordinary, ps_ordinary)) {
    expect_null(attr(res, "interval_degenerate"))
    expect_null(dmar_metadata(res)$interval_degenerate)
    expect_false(any(grepl("Limits are not a confidence interval",
                           utils::capture.output(print(res)), fixed = TRUE)))
    expect_false(any(generics::tidy(res)$interval_degenerate))
    expect_match(results_sentence(res), "95% CI [", fixed = TRUE)
  }
  expect_s3_class(cd_ordinary, c("dmar_cliff_delta", "dmar_tbl",
                                 "data.frame"), exact = TRUE)
  expect_s3_class(ps_ordinary, c("dmar_probability_of_superiority_paired",
                                 "dmar_tbl", "data.frame"), exact = TRUE)
  # The verbs of the long tables are the default ones plus the column.
  plain <- cd_ordinary
  class(plain) <- c("dmar_tbl", "data.frame")
  td <- generics::tidy(cd_ordinary)
  expect_identical(td[names(td) != "interval_degenerate"],
                   generics::tidy(plain))
})
