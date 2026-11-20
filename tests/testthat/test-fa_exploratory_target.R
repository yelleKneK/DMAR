# Rotation to a partially specified target in fa_exploratory() (Browne,
# 1972a, 1972b, 2001): the two forms of the target, the identification
# checks lavaan does not make, loadings the rotation attains exactly, and
# the display.

target_items <- function() c("t1_visual_perception", "t2_cubes", "t4_lozenges",
  "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")

# The key list's target as a matrix: each item at zero on the factor it is
# not keyed to.
target_keys_matrix <- function() {
  items <- target_items()
  B <- matrix(NA_real_, 6, 2, dimnames = list(items, c("spatial", "verbal")))
  B[4:6, "spatial"] <- 0
  B[1:3, "verbal"] <- 0
  B
}

# The condition an expression signals, for its class and fields.
target_condition <- function(expr) tryCatch(expr, error = function(e) e)

# Another oblique rotation of a two-factor solution L with factor
# correlation phi: with T0'T0 = Phi (unit columns), A = L T0' is an
# orthogonal unrotated form of the solution, and turning the jth column of
# T0 by `angle` gives another oblique rotation of it, A T1'^-1.
target_turn <- function(L, phi, angle, j) {
  T0 <- chol(matrix(c(1, phi, phi, 1), 2))
  a <- atan2(T0[2, j], T0[1, j]) + angle
  T1 <- T0
  T1[, j] <- c(cos(a), sin(a))
  L %*% t(T0) %*% t(solve(T1))
}

# The standardized loadings of factor f, in the order of `items`, and the
# correlation of f1 and f2.
target_loadings <- function(x, f, items) {
  s <- x[x$op == "=~" & x$metric == "standardized" & x$lhs == f, ]
  s$estimate[match(items, s$rhs)]
}
target_phi <- function(x) {
  x$estimate[x$op == "~~" & x$metric == "standardized" & x$lhs == "f1" & x$rhs == "f2"]
}

expect_target_refusal <- function(expr, class, reason) {
  e <- target_condition(expr)
  expect_s3_class(e, c(paste0("dmar_", class), "dmar_refusal", "error"))
  expect_identical(e$reason, reason)
  invisible(e)
}

test_that("the key list and the matrix give the same target rotation", {
  skip_if_not_installed("lavaan")
  items <- target_items()
  keys <- list(spatial = items[1:3], verbal = items[4:6])
  B <- target_keys_matrix()
  # The items and the factor count come from the target.
  x <- fa_exploratory(holzinger_swineford, rotation = "target", target = keys, seed = 113)
  m <- attr(x, "method")
  expect_identical(m$rotation, "target")
  expect_false(m$orthogonal)
  expect_identical(m$row_weights, "none")
  expect_identical(m$target, B)
  expect_identical(m$n_specified, 6L)
  expect_identical(unique(x$rhs[x$op == "=~"]), items)
  expect_identical(attr(x, "variance_accounted")$factor, c("f1", "f2", "total"))
  # The matrix gives the same table, and so does naming the items and the
  # count. A matrix in another row order is put in the order of the items
  # named; with none named, its rows are the items, in its order.
  strip <- function(x) { attr(x, "fit") <- NULL; x }
  y <- fa_exploratory(holzinger_swineford, rotation = "target", target = B, seed = 113)
  expect_identical(strip(y), strip(x))
  z <- fa_exploratory(holzinger_swineford, nfactors = 2, items = items,
                      rotation = "target", target = keys, seed = 113)
  expect_identical(strip(z), strip(x))
  shuffled <- B[c(6, 2, 4, 1, 5, 3), ]
  w <- fa_exploratory(holzinger_swineford, items = items, rotation = "target",
                      target = shuffled, seed = 113)
  expect_identical(strip(w), strip(x))
  expect_identical(attr(w, "method")$target, B)
  v <- fa_exploratory(holzinger_swineford, rotation = "target", target = shuffled,
                      seed = 113)
  expect_identical(unique(v$rhs[v$op == "=~"]), rownames(shuffled))
  expect_identical(attr(v, "method")$target, shuffled)
  # lavaan received the target, the mask that marks the specified entries,
  # no row weights, and the columns in the target's order.
  rot <- lavaan::lavInspect(attr(x, "fit"), "options")$rotation.args
  expect_identical(lavaan::lavInspect(attr(x, "fit"), "options")$rotation, "pst")
  expect_equal(unname(rot$target_mask), unname(1 * !is.na(B)))
  expect_identical(rot$order_lv_by, "none")
  expect_true(rot$reflect)
  # The rotated loadings are the least squares rotation to the target: no
  # other rotation brings the specified loadings closer to zero. A small
  # oblique turn of the solution (Lambda T^-1' with T's columns of unit
  # length, Phi = T'T) raises the criterion.
  s <- x[x$op == "=~" & x$metric == "standardized", ]
  L <- matrix(s$estimate, 6, 2, dimnames = list(items, c("f1", "f2")))
  phi <- x$estimate[x$op == "~~" & x$metric == "standardized" & x$lhs == "f1" &
                      x$rhs == "f2"]
  criterion <- function(L) sum(L[!is.na(B)]^2)
  base <- criterion(L)
  expect_equal(target_turn(L, phi, 0, 1), L, ignore_attr = TRUE)
  for (angle in c(-0.02, 0.02)) for (j in 1:2) {
    expect_gt(criterion(target_turn(L, phi, angle, j)), base)
  }
})

test_that("the key list targets zero off the keys and leaves the keys free", {
  items <- target_items()
  # An item keyed to two factors is unspecified on both; an analyzed item
  # keyed to none is unspecified on every factor.
  B <- DMAR:::.efa_target_matrix(list(a = items[1:4], b = items[3:6]), items)
  expect_identical(dimnames(B), list(items, c("a", "b")))
  expect_true(all(is.na(B[3:4, ])))
  expect_identical(unname(B[c(1, 2, 5, 6), "a"]), c(NA, NA, 0, 0))
  expect_identical(unname(B[c(1, 2, 5, 6), "b"]), c(0, 0, NA, NA))
  B <- DMAR:::.efa_target_matrix(list(a = items[1:2], b = items[4:5]), items)
  expect_true(all(is.na(B[c(3, 6), ])))
  expect_identical(unname(B[1:2, "b"]), c(0, 0))
  # Without the items, the keyed items in order of first appearance.
  expect_identical(rownames(DMAR:::.efa_target_matrix(list(a = items[4:6],
                                                           b = items[1:3]))),
                   items[c(4:6, 1:3)])
  # A matrix without column names names its factors f1, f2, ...
  B <- unname(target_keys_matrix())
  rownames(B) <- items
  expect_identical(colnames(DMAR:::.efa_target_matrix(B, items)), c("f1", "f2"))
})

test_that("a target DMAR cannot read is refused by class, with its reason", {
  skip_if_not_installed("lavaan")
  items <- target_items()
  d <- holzinger_swineford[items]
  keys <- list(spatial = items[1:3], verbal = items[4:6])
  B <- target_keys_matrix()
  expect_target_refusal(fa_exploratory(d, rotation = "target"),
                        "target_invalid", "rotation")
  expect_target_refusal(fa_exploratory(d, 2, target = keys),
                        "target_invalid", "rotation")
  expect_target_refusal(fa_exploratory(d, rotation = "target", target = "spatial"),
                        "target_invalid", "form")
  expect_target_refusal(fa_exploratory(d, rotation = "target",
                                       target = as.data.frame(B)),
                        "target_invalid", "form")
  # Fewer than two factors.
  expect_target_refusal(fa_exploratory(d, rotation = "target",
                                       target = list(spatial = items[1:3])),
                        "target_invalid", "factors")
  expect_target_refusal(fa_exploratory(d, rotation = "target",
                                       target = B[, 1, drop = FALSE]),
                        "target_invalid", "factors")
  # Names and dimensions that do not match the items.
  expect_target_refusal(fa_exploratory(d, rotation = "target",
                                       target = list(spatial = c(items[1:3], "t3_paper_form_board"),
                                                     verbal = items[4:6]),
                                       items = items),
                        "target_invalid", "names")
  expect_target_refusal(fa_exploratory(d, rotation = "target",
                                       target = list(spatial = c(items[1:3], "not_an_item"),
                                                     verbal = items[4:6])),
                        "target_invalid", "names")
  expect_target_refusal(fa_exploratory(d, rotation = "target", target = B[1:5, ],
                                       items = items),
                        "target_invalid", "names")
  expect_target_refusal(fa_exploratory(d, rotation = "target", target = unname(B)),
                        "target_invalid", "names")
  expect_target_refusal(fa_exploratory(d, rotation = "target",
                                       target = `colnames<-`(B, c("a", "a"))),
                        "target_invalid", "names")
  expect_target_refusal(fa_exploratory(d, rotation = "target",
                                       target = list(spatial = items[1:3], items[4:6])),
                        "target_invalid", "names")
  # Values that are not finite, or not standardized loadings.
  for (v in c(Inf, NaN, 1, -1.5)) {
    bad <- B
    bad[1, 2] <- v
    expect_target_refusal(fa_exploratory(d, rotation = "target", target = bad),
                          "target_invalid", "values")
  }
  # No specified entry.
  none <- B
  none[] <- NA
  expect_target_refusal(fa_exploratory(d, 2, rotation = "target", target = none),
                        "target_invalid", "none_specified")
  # A factor count that disagrees with the target.
  expect_target_refusal(fa_exploratory(d, 3, rotation = "target", target = keys),
                        "target_invalid", "nfactors")
})

test_that("a target that does not identify the rotation is refused before fitting", {
  skip_if_not_installed("lavaan")
  items <- target_items()
  d <- holzinger_swineford[items]
  B <- target_keys_matrix()
  # Oblique: every column needs m - 1 = 1 specified zero. lavaan alone
  # would return this rotation with finite standard errors and no warning.
  one <- B
  one[, "verbal"] <- NA
  e <- expect_target_refusal(fa_exploratory(d, rotation = "target", target = one),
                             "target_unidentified", "counts")
  expect_identical(e$factors, "verbal")
  expect_match(conditionMessage(e), "verbal has 0, needs 1", fixed = TRUE)
  # Orthogonal rotation needs fewer: 0 and 1, so the same target is
  # identified there, and two empty columns are not.
  x <- fa_exploratory(d, rotation = "target_orthogonal", target = one, seed = 113)
  expect_true(all(is.finite(x$se[x$op == "=~"])))
  expect_identical(attr(x, "method")$n_specified, 3L)
  sparse <- B
  sparse[] <- NA
  sparse[4, 1] <- 0
  expect_no_error(DMAR:::.efa_target_counts(sparse, orthogonal = TRUE))
  expect_target_refusal(DMAR:::.efa_target_counts(sparse, orthogonal = FALSE),
                        "target_unidentified", "counts")
  # A column with a nonzero target value needs one entry more (Algina,
  # 1980, Theorems 1 and 3): m = 2 in oblique rotation.
  nonzero <- B
  nonzero[, "verbal"] <- NA
  nonzero[1, "verbal"] <- 0.3
  expect_target_refusal(fa_exploratory(d, rotation = "target", target = nonzero),
                        "target_unidentified", "counts")
  nonzero[2, "verbal"] <- 0
  expect_no_error(DMAR:::.efa_target_counts(nonzero, orthogonal = FALSE))
  # The orthogonal counts, three factors: sorted, at least 0, 1, 2.
  B3 <- matrix(NA_real_, 9, 3, dimnames = list(paste0("y", 1:9), c("a", "b", "c")))
  B3[1:2, "c"] <- 0
  expect_target_refusal(DMAR:::.efa_target_counts(B3, orthogonal = TRUE),
                        "target_unidentified", "counts")
  B3[3, "b"] <- 0
  expect_no_error(DMAR:::.efa_target_counts(B3, orthogonal = TRUE))
  # Two columns specified alike cannot be told apart (lavaan alone stops
  # with a message about degenerate random starts).
  twin <- B
  twin[] <- NA
  twin[4:6, ] <- 0
  e <- expect_target_refusal(fa_exploratory(d, rotation = "target", target = twin),
                             "target_unidentified", "duplicate")
  expect_identical(e$factors, c("spatial", "verbal"))
})

test_that("a target whose specified rows are rank deficient is refused", {
  skip_if_not_installed("lavaan")
  # A covariance matrix built exactly from three correlated clusters of
  # three items (interleaved, so that lavaan's unrotated solution is
  # regular). Items y2 and y5 load on the second cluster's factor alone, so
  # their unrotated loadings are proportional: rank 1, below the m - 1 = 2
  # an oblique column needs.
  L <- matrix(0, 9, 3)
  L[1:3, 1] <- c(.7, .6, .8)
  L[4:6, 2] <- c(.75, .65, .7)
  L[7:9, 3] <- c(.6, .7, .8)
  L <- L[c(1, 4, 7, 2, 5, 8, 3, 6, 9), ]
  S <- L %*% matrix(c(1, .3, .3, .3, 1, .3, .3, .3, 1), 3) %*% t(L)
  diag(S) <- 1
  dimnames(S) <- list(paste0("y", 1:9), paste0("y", 1:9))
  B <- matrix(NA_real_, 9, 3, dimnames = list(paste0("y", 1:9), c("a", "b", "c")))
  B[c(2, 5), "a"] <- 0
  B[c(1, 3), "b"] <- 0
  B[c(1, 2), "c"] <- 0
  # The rank check itself, on the population loadings, where the rank of
  # rows y2 and y5 is exactly 1, and on a turn of them (the rank does not
  # depend on the rotation). Orthogonal rotation asks ranks 0, 1, and 2 of
  # the columns in increasing order, which these rows meet.
  e <- expect_target_refusal(DMAR:::.efa_target_rank(B, L, FALSE),
                             "target_unidentified", "rank")
  expect_identical(e$factors, "a")
  turned <- L %*% qr.Q(qr(matrix(c(2, 1, 0, -1, 3, 1, 0, 2, 4), 3)))
  expect_target_refusal(DMAR:::.efa_target_rank(B, turned, FALSE),
                        "target_unidentified", "rank")
  expect_no_error(DMAR:::.efa_target_rank(B, turned, TRUE))
  # The same refusal through fa_exploratory(), which reads the rank from
  # the unrotated fit, whose loadings meet the deficiency only to the
  # optimizer's precision (about 3e-9 of the largest singular value,
  # against the threshold of 1e-6).
  e <- expect_target_refusal(fa_exploratory(S = S, N = 300, rotation = "target",
                                            target = B),
                             "target_unidentified", "rank")
  expect_identical(e$factors, "a")
  expect_match(conditionMessage(e), "column a have rank 1, below the 2", fixed = TRUE)
  # Items from two clusters separate the factors, and the rotation runs.
  # This target is just identified, with its zeros where the population
  # loadings are zero, so the rotation recovers the population structure:
  # the loadings, the factor correlations of 0.3, and the six specified
  # zeros attained exactly, fixed at zero in both metrics.
  B[c(2, 5), "a"] <- NA
  B[c(2, 3), "a"] <- 0
  x <- fa_exploratory(S = S, N = 300, rotation = "target", target = B, seed = 113)
  expect_identical(attr(x, "method")$rotation, "target")
  fit <- attr(x, "fit")
  expect_lt(max(abs(unname(unclass(lavaan::lavInspect(fit, "std")$lambda)) - L)), 1e-6)
  factors <- c("f1", "f2", "f3")
  phi <- x[x$op == "~~" & x$metric == "standardized" & x$lhs %in% factors &
             x$rhs %in% factors & x$lhs != x$rhs, ]
  expect_identical(nrow(phi), 3L)
  expect_equal(phi$estimate, rep(0.3, 3), tolerance = 1e-6)
  fixed <- x[x$op == "=~" & is.na(x$se), ]
  expect_identical(nrow(fixed), 12L)
  expect_setequal(paste(fixed$lhs, fixed$rhs),
                  c("f1 y2", "f1 y3", "f2 y1", "f2 y3", "f3 y1", "f3 y2"))
  expect_identical(table(fixed$metric)[["original"]], 6L)
  expect_true(all(fixed$estimate == 0 & fixed$ci_lower == 0 & fixed$ci_upper == 0))
  expect_match(attr(x, "analysis_note"), "Every specified loading is attained exactly",
               fixed = TRUE)
})

test_that("a just identified oblique target is the confirmatory model with those zeros", {
  skip_if_not_installed("lavaan")
  items <- target_items()
  d <- holzinger_swineford[items]
  # m - 1 = 1 zero per column: the rotation attains both, and the solution
  # is the CFA with those two loadings fixed at zero and the factor
  # variances at one.
  B <- matrix(NA_real_, 6, 2, dimnames = list(items, c("spatial", "verbal")))
  B["t7_sentence", "spatial"] <- 0
  B["t4_lozenges", "verbal"] <- 0
  x <- fa_exploratory(d, rotation = "target", target = B, seed = 113)
  model <- paste0(
    "spatial =~ t1_visual_perception + t2_cubes + t4_lozenges + ",
    "t6_paragraph_comprehension + 0*t7_sentence + t9_word_meaning\n",
    "verbal =~ t1_visual_perception + t2_cubes + 0*t4_lozenges + ",
    "t6_paragraph_comprehension + t7_sentence + t9_word_meaning")
  cfa <- lavaan::cfa(model, data = d, std.lv = TRUE, missing = "ml")
  key <- function(lhs, rhs) paste(lhs, rhs)
  lv <- c(f1 = "spatial", f2 = "verbal")
  for (metric in c("standardized", "original")) {
    tab <- if (metric == "standardized") lavaan::standardizedSolution(cfa) else
      lavaan::parameterEstimates(cfa)
    est <- if (metric == "standardized") tab$est.std else tab$est
    load <- x[x$op == "=~" & x$metric == metric, ]
    i <- match(key(lv[load$lhs], load$rhs), key(tab$lhs, tab$rhs))
    expect_equal(load$estimate, est[i], tolerance = 1e-5, info = metric)
    free <- tab$se[i] > 0
    expect_equal(load$se[free], tab$se[i][free], tolerance = 1e-4, info = metric)
    # The two zeros are fixed by the target: at zero, with no standard
    # error, statistic, or p-value, and both limits at zero.
    fixed <- load[!free, ]
    expect_identical(key(fixed$lhs, fixed$rhs),
                     c("f1 t7_sentence", "f2 t4_lozenges"), info = metric)
    expect_identical(fixed$estimate, c(0, 0), info = metric)
    expect_true(all(is.na(c(fixed$se, fixed$statistic, fixed$p_value))), info = metric)
    expect_identical(c(fixed$ci_lower, fixed$ci_upper), rep(0, 4), info = metric)
  }
  phi <- x[x$op == "~~" & x$metric == "standardized" & x$lhs == "f1" & x$rhs == "f2", ]
  sc <- lavaan::standardizedSolution(cfa)
  phi_cfa <- sc[sc$lhs == "spatial" & sc$op == "~~" & sc$rhs == "verbal", ]
  expect_equal(phi$estimate, phi_cfa$est.std, tolerance = 1e-5)
  expect_equal(phi$se, phi_cfa$se, tolerance = 1e-4)
  fm <- attr(x, "fit_measures")
  expect_equal(fm$value[fm$term == "chi_square"], unname(lavaan::fitMeasures(cfa, "chisq")),
               tolerance = 1e-6)
  expect_equal(fm$value[fm$term == "df"], unname(lavaan::fitMeasures(cfa, "df")))
  expect_match(attr(x, "analysis_note"), "Every specified loading is attained exactly",
               fixed = TRUE)
  expect_match(attr(x, "analysis_note"), "confirmatory model with those loadings fixed",
               fixed = TRUE)
  # The same holds for an orthogonal rotation with one zero, against the
  # CFA with uncorrelated factors.
  Bo <- B
  Bo["t7_sentence", "spatial"] <- NA
  xo <- fa_exploratory(d, rotation = "target_orthogonal", target = Bo, seed = 113)
  cfa_o <- lavaan::cfa(sub("0\\*t7_sentence", "t7_sentence", model), data = d,
                       std.lv = TRUE, missing = "ml", orthogonal = TRUE)
  so <- lavaan::standardizedSolution(cfa_o)
  load <- xo[xo$op == "=~" & xo$metric == "standardized", ]
  i <- match(key(lv[load$lhs], load$rhs), key(so$lhs, so$rhs))
  expect_equal(load$estimate, so$est.std[i], tolerance = 1e-5)
  expect_identical(sum(is.na(load$se)), 1L)
  phi <- xo[xo$op == "~~" & xo$lhs == "f1" & xo$rhs == "f2", ]
  expect_true(nrow(phi) == 0L || isTRUE(all.equal(phi$estimate, rep(0, nrow(phi)))))
})

test_that("an overidentified target estimates every loading, and a negative target keeps its sign", {
  skip_if_not_installed("lavaan")
  items <- target_items()
  d <- holzinger_swineford[items]
  x <- fa_exploratory(d, rotation = "target",
                      target = list(spatial = items[1:3], verbal = items[4:6]), seed = 113)
  expect_true(all(is.finite(x$se[x$op == "=~"])))
  expect_match(attr(x, "analysis_note"), "^Intervals condition on the selected factor count and on a target")
  expect_false(grepl("attained exactly", attr(x, "analysis_note")))
  # A target of -0.8 for the verbal items on the second factor: with a
  # nonzero value specified no column is reflected, so the verbal loadings
  # stay negative (lavaan's reflection would turn the column to a positive
  # sum).
  B <- target_keys_matrix()
  B[4:6, "verbal"] <- -0.8
  y <- fa_exploratory(d, rotation = "target", target = B, seed = 113)
  rot <- lavaan::lavInspect(attr(y, "fit"), "options")$rotation.args
  expect_false(rot$reflect)
  expect_equal(unname(rot$target_mask), unname(1 * !is.na(B)))
  expect_true(all(target_loadings(y, "f2", items)[4:6] < -0.5))
  # The nonzero values are specified, not dropped from the mask: the
  # rotated loadings minimize the squared differences from every specified
  # value, the -0.8s among them, and a small turn of the solution raises
  # that criterion.
  L <- cbind(target_loadings(y, "f1", items), target_loadings(y, "f2", items))
  criterion <- function(L) sum((L - B)[!is.na(B)]^2)
  for (angle in c(-0.02, 0.02)) for (j in 1:2) {
    expect_gt(criterion(target_turn(L, target_phi(y), angle, j)), criterion(L))
  }
  # The same communalities and fit as any other rotation.
  expect_equal(y$estimate[y$op == "communality"], x$estimate[x$op == "communality"],
               tolerance = 1e-6)
  expect_output(print(y), "6 at 0, 3 nonzero", fixed = TRUE)
})

test_that("with a nonzero value no factor is reflected, and a factor of zeros takes its start's sign", {
  skip_if_not_installed("lavaan")
  items <- target_items()
  d <- holzinger_swineford[items]
  strip <- function(x) { attr(x, "fit") <- NULL; x }
  # lavaan reflects all the columns or none. With a -0.8 specified on the
  # second factor, none is: that factor keeps the target's sign in every
  # run, while the first, whose specified values are all zero, takes the
  # sign of the random start that reached the minimum. The solution is the
  # same up to that sign, which also turns the factor correlation.
  B <- target_keys_matrix()
  B[4:6, "verbal"] <- -0.8
  runs <- lapply(1:3, function(s) fa_exploratory(d, rotation = "target", target = B,
                                                 seed = s))
  first <- runs[[1]]
  for (x in runs) {
    expect_false(lavaan::lavInspect(attr(x, "fit"), "options")$rotation.args$reflect)
    expect_equal(target_loadings(x, "f2", items), target_loadings(first, "f2", items),
                 tolerance = 1e-6)
    turn <- sign(sum(target_loadings(x, "f1", items))) *
      sign(sum(target_loadings(first, "f1", items)))
    expect_equal(turn * target_loadings(x, "f1", items),
                 target_loadings(first, "f1", items), tolerance = 1e-6)
    expect_equal(turn * target_phi(x), target_phi(first), tolerance = 1e-6)
  }
  # A seed fixes the sign.
  expect_identical(strip(fa_exploratory(d, rotation = "target", target = B, seed = 2)),
                   strip(runs[[2]]))
  # With every specified value zero, every factor is reflected to a
  # positive sum of loadings, whatever the seed.
  for (s in 1:2) {
    z <- fa_exploratory(d, rotation = "target", target = target_keys_matrix(), seed = s)
    expect_true(lavaan::lavInspect(attr(z, "fit"), "options")$rotation.args$reflect)
    expect_gt(sum(target_loadings(z, "f1", items)), 0)
    expect_gt(sum(target_loadings(z, "f2", items)), 0)
  }
})

test_that("raw data and a covariance matrix analyze the same items", {
  skip_if_not_installed("lavaan")
  items <- target_items()
  nine <- c(items, "t10_addition", "t12_counting_groups_of_dots",
            "t13_straight_and_curved_capitals")
  d <- holzinger_swineford[nine]
  expect_false(anyNA(d))
  S <- stats::cov(d)
  keys <- list(spatial = items[1:3], verbal = items[4:6])
  # With a target and no items, the items the target names, from either
  # input; the other columns of S are not analyzed.
  x <- fa_exploratory(d, rotation = "target", target = keys, seed = 113)
  y <- fa_exploratory(S = S, N = nrow(d), rotation = "target", target = keys, seed = 113)
  expect_identical(unique(y$rhs[y$op == "=~"]), items)
  expect_identical(unique(x$rhs[x$op == "=~"]), unique(y$rhs[y$op == "=~"]))
  expect_identical(attr(y, "method")$target, attr(x, "method")$target)
  # Complete data: the same estimates, in both metrics.
  key <- function(t) paste(t$term, t$lhs, t$op, t$rhs, t$metric)
  expect_false(anyDuplicated(key(x)) > 0L)
  i <- match(key(y), key(x))
  expect_false(anyNA(i))
  expect_equal(y$estimate, x$estimate[i], tolerance = 1e-5)
  # Named items reduce S under any rotation, as they select columns of data.
  z <- fa_exploratory(S = S, N = nrow(d), nfactors = 2, items = items, seed = 113)
  w <- fa_exploratory(d, nfactors = 2, items = items, seed = 113)
  expect_identical(unique(z$rhs[z$op == "=~"]), items)
  zl <- z[z$op == "=~" & z$metric == "standardized", ]
  wl <- w[w$op == "=~" & w$metric == "standardized", ]
  expect_equal(zl$estimate, wl$estimate[match(key(zl), key(wl))], tolerance = 1e-5)
  expect_error(fa_exploratory(S = S, N = nrow(d), nfactors = 2,
                              items = c(items, "t3_paper_form_board")),
               "at least three distinct items of 'S'", fixed = TRUE,
               class = "dmar_items_invalid")
  expect_target_refusal(fa_exploratory(S = S[items[-1], items[-1]], N = nrow(d),
                                       rotation = "target", target = keys),
                        "target_invalid", "names")
})

test_that("the target rotation display states the target under the table", {
  skip_if_not_installed("lavaan")
  items <- target_items()
  d <- holzinger_swineford[items]
  x <- fa_exploratory(d, rotation = "target",
                      target = list(spatial = items[1:3], verbal = items[4:6]), seed = 113)
  local_reproducible_output(width = 200)
  out <- capture.output(print(x))
  expect_match(out, paste0("rotation: target \\(oblique\\), 6 of 12 loadings specified ",
                           "\\(all 0\\); SE: standard"), all = FALSE)
  head <- which(out == "Target (the value where a loading is specified, ? where it is not):")
  expect_length(head, 1L)
  expect_match(out[head + 1L], "^ +f1 f2$")
  expect_match(out[head + 2L], "^t1_visual_perception +\\? +0$")
  expect_match(out[head + 5L], "^t6_paragraph_comprehension +0 +\\?$")
  expect_identical(out[head + 8L],
                   "Factors: f1 = spatial, f2 = verbal (the table keeps lavaan's names).")
  # Beneath the loadings, before the variance accounted for.
  expect_true(head > which(grepl("^ *item +f1 +f2 +communality", out)))
  expect_true(head < which(grepl("^Variance accounted for", out)))
  expect_match(out, "condition on the selected factor count and on a target fixed",
               all = FALSE)
  expect_false(any(grepl("departure", out, ignore.case = TRUE)))
  expect_identical(generics::glance(x)$rotation, "target")
  # The orthogonal form, and a matrix without column names (no mapping line).
  B <- unname(target_keys_matrix())
  rownames(B) <- items
  o <- fa_exploratory(d, rotation = "target_orthogonal", target = B, seed = 113)
  oo <- capture.output(print(o))
  expect_match(oo, "rotation: target \\(orthogonal\\), 6 of 12 loadings specified", all = FALSE)
  expect_match(oo, "Factor correlations: fixed at zero \\(orthogonal rotation\\)", all = FALSE)
  expect_false(any(grepl("^Factors: ", oo)))
  expect_true(any(grepl("^Target \\(the value", oo)))
  # knit_print() sets the target off as fixed-width text.
  skip_if_not_installed("knitr")
  kp <- as.character(knitr::knit_print(x))
  expect_match(kp, "Target (the value where a loading is specified", fixed = TRUE)
  expect_match(kp, "```\n +f1 f2\nt1_visual_perception", all = FALSE)
})
