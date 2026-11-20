# factor_congruence(): Tucker's congruence coefficients between the factors
# of two solutions, matched one to one and reflected. The oracle is
# psych::fa.congruence(), whose values are pinned below because psych is
# not a dependency; tools/oracle_snippets/psych_congruence.R reruns the
# comparison live. This file is in the CRAN subset (tests/testthat.R) for
# the pinned values and the assignment solver's brute force check; the
# blocks that fit lavaan models skip on CRAN to keep it near a second.

# Two fixed loading matrices on eight items, x with three factors and y
# with three, whose best pairing is X1-Y2 (negative), X2-Y1, X3-Y3.
cong_Lx <- matrix(c(.71, .64, .58, .12, .05, -.08, .20, .11,
                    .09, .14, -.03, .68, .74, .61, .02, .19,
                    .05, -.10, .21, .07, .12, .03, .66, .59), 8, 3,
                  dimnames = list(paste0("v", 1:8), c("X1", "X2", "X3")))
cong_Ly <- matrix(c(.10, .18, .02, .70, .69, .55, .08, .25,
                    -.62, -.70, -.49, -.15, .02, .11, -.30, -.04,
                    .12, -.05, .30, .02, .20, -.06, .58, .71), 8, 3,
                  dimnames = list(paste0("v", 1:8), c("Y1", "Y2", "Y3")))

# All permutations of 1:n, one per row.
cong_perms <- function(n) {
  if (n == 1L) return(matrix(1L))
  p <- cong_perms(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(k)
    cbind(k, ifelse(p >= k, p + 1L, p))))
}

# The brute force answer to .assignment_max(): the largest total over every
# one-to-one assignment of the zero-padded square matrix, and among the
# assignments that reach it the lexicographically first, with a padding
# column read as "unassigned" (NA).
cong_brute <- function(score) {
  nr <- nrow(score)
  nc <- ncol(score)
  n <- max(nr, nc)
  S <- matrix(0, n, n)
  S[seq_len(nr), seq_len(nc)] <- score
  P <- cong_perms(n)
  total <- apply(P, 1L, function(a) sum(S[cbind(seq_len(n), a)]))
  best <- max(total)
  tol <- sqrt(.Machine$double.eps) * max(1, abs(best))
  ok <- P[total >= best - tol, seq_len(nr), drop = FALSE]
  first <- ok[do.call(order, as.data.frame(ok)), , drop = FALSE][1L, ]
  first[first > nc] <- NA_integer_
  list(best = best, assignment = as.integer(first))
}

test_that("the congruence matrix is the uncentered cosine and matches psych", {
  res <- factor_congruence(cong_Lx, cong_Ly)
  C <- attr(res, "congruence_matrix")
  hand <- crossprod(cong_Lx, cong_Ly) /
    sqrt(outer(colSums(cong_Lx^2), colSums(cong_Ly^2)))
  expect_equal(C, hand, tolerance = 1e-14)
  expect_identical(dimnames(C), list(c("X1", "X2", "X3"), c("Y1", "Y2", "Y3")))
  # psych::fa.congruence(cong_Lx, cong_Ly, digits = 15), psych 2.6.5,
  # by rows of x.
  psych_C <- matrix(c(0.23347238276659499, -0.98440209920224497,
                      0.38252099289133701, 0.99394757310223403,
                      -0.12893655829710499, 0.22232199567199501,
                      0.31178929012075501, -0.280910653898859,
                      0.973654104075802), 3, 3, byrow = TRUE)
  expect_equal(unname(C), psych_C, tolerance = 1e-12)
  # Not the Pearson correlation of the columns: loadings that differ by a
  # constant correlate perfectly but are not congruent.
  a <- c(.4, .4, 0)
  b <- c(.1, .1, -.3)
  phi <- factor_congruence(a, b)$congruence
  expect_equal(phi, .08 / sqrt(.32 * .11), tolerance = 1e-14)
  expect_equal(stats::cor(a, b), 1)
})

test_that("the fixed matrices are matched, reflected, and ordered as x", {
  res <- factor_congruence(cong_Lx, cong_Ly)
  expect_s3_class(res, c("dmar_factor_congruence", "dmar_tbl", "data.frame"),
                  exact = TRUE)
  expect_identical(names(res), c("factor_x", "factor_y", "congruence",
                                 "reflected", "next_largest"))
  expect_identical(res$factor_x, c("X1", "X2", "X3"))
  expect_identical(res$factor_y, c("Y2", "Y1", "Y3"))
  C <- attr(res, "congruence_matrix")
  expect_equal(res$congruence, abs(c(C[1, 2], C[2, 1], C[3, 3])))
  expect_identical(res$reflected, c(TRUE, FALSE, FALSE))
  expect_equal(res$next_largest,
               c(max(abs(C[1, -2])), max(abs(C[2, -1])), max(abs(C[3, -3]))))
  aligned <- attr(res, "loadings_aligned")
  expect_identical(colnames(aligned), c("Y2", "Y1", "Y3"))
  expect_equal(aligned[, "Y2"], -cong_Ly[, "Y2"])
  expect_equal(aligned[, c("Y1", "Y3")], cong_Ly[, c("Y1", "Y3")])
  expect_null(attr(res, "Phi_aligned"))
  expect_identical(attr(res, "items"), paste0("v", 1:8))
  expect_identical(attr(res, "match"), "congruence")
  expect_identical(attr(res, "loadings"), "pattern")
  expect_identical(attr(res, "inputs"),
                   list(x = "matrix", y = "matrix", items = "name"))
})

test_that("a solution compared with itself has coefficients of 1", {
  res <- factor_congruence(cong_Lx, cong_Lx)
  expect_equal(res$congruence, rep(1, 3), tolerance = 1e-14)
  expect_identical(res$factor_y, res$factor_x)
  expect_identical(res$reflected, rep(FALSE, 3))
  one <- factor_congruence(cong_Lx[, 1], cong_Lx[, 1])
  expect_equal(one$congruence, 1, tolerance = 1e-14)
  expect_identical(one$next_largest, NA_real_)
})

test_that("a permuted, reflected, reordered copy is recovered exactly", {
  copy <- cong_Lx[c(8, 3, 5, 1, 7, 2, 6, 4), c("X3", "X1", "X2")]
  copy[, "X1"] <- -copy[, "X1"]
  colnames(copy) <- c("c", "a", "b")
  res <- factor_congruence(cong_Lx, copy)
  expect_identical(res$factor_y, c("a", "b", "c"))
  expect_identical(res$reflected, c(TRUE, FALSE, FALSE))
  expect_equal(res$congruence, rep(1, 3), tolerance = 1e-14)
  aligned <- attr(res, "loadings_aligned")
  expect_equal(unname(aligned), unname(cong_Lx), tolerance = 1e-15)
  expect_identical(rownames(aligned), rownames(cong_Lx))
  # Rescaling a column by a positive constant changes nothing.
  scaled <- sweep(copy, 2L, c(2, 0.5, 10), `*`)
  expect_equal(factor_congruence(cong_Lx, scaled)$congruence,
               res$congruence, tolerance = 1e-14)
})

test_that("unnamed items are aligned by position", {
  res <- factor_congruence(unname(cong_Lx), unname(cong_Ly))
  expect_identical(attr(res, "items"), 1:8)
  expect_identical(attr(res, "inputs")$items, "position")
  expect_identical(res$factor_x, c("F1", "F2", "F3"))
  expect_equal(unname(attr(res, "congruence_matrix")),
               unname(attr(factor_congruence(cong_Lx, cong_Ly),
                           "congruence_matrix")))
  # One named and one unnamed: by position, with the names that exist.
  half <- factor_congruence(cong_Lx, unname(cong_Ly))
  expect_identical(attr(half, "items"), rownames(cong_Lx))
  expect_identical(attr(half, "inputs")$items, "position")
  expect_error(factor_congruence(unname(cong_Lx), unname(cong_Ly)[-1, ]),
               "matched by position")
})

test_that(".assignment_max() agrees with brute force, square and rectangular", {
  skip_on_cran()
  # The 240 comparisons are gathered and checked by one expectation per
  # property, which keeps the file fast enough for the CRAN subset; a
  # failure names the matrix by its position in the list.
  set.seed(20261105)
  n_it <- 240L
  got <- want <- vector("list", n_it)
  total_got <- total_want <- numeric(n_it)
  one_to_one <- logical(n_it)
  for (it in seq_len(n_it)) {
    nr <- sample(1:6, 1)
    nc <- sample(1:6, 1)
    # Every third matrix has many ties, which the lexicographic rule settles.
    score <- if (it %% 3 == 0) {
      matrix(sample(0:2, nr * nc, replace = TRUE), nr, nc)
    } else {
      matrix(stats::runif(nr * nc), nr, nc)
    }
    got[[it]] <- DMAR:::.assignment_max(score)
    brute <- cong_brute(score)
    want[[it]] <- brute$assignment
    used <- !is.na(got[[it]])
    total_got[it] <- sum(score[cbind(which(used), got[[it]][used])])
    total_want[it] <- brute$best
    one_to_one[it] <- !anyDuplicated(got[[it]][used]) &&
      sum(used) == min(nr, nc)
  }
  expect_identical(got, want)
  expect_equal(total_got, total_want, tolerance = 1e-12)
  expect_identical(which(!one_to_one), integer(0))
})

test_that(".assignment_max() is optimal where greedy matching is not", {
  # Row-wise best partners collide (both rows prefer column 1), and taking
  # the largest entry first gives .90 + .10 = 1.00; the optimal pairing
  # gives .80 + .85 = 1.65.
  score <- matrix(c(.90, .85, .80, .10), 2, 2)
  expect_identical(DMAR:::.assignment_max(score), c(2L, 1L))
  # Ties keep the listed order.
  expect_identical(DMAR:::.assignment_max(matrix(1, 4, 4)), 1:4)
  expect_identical(DMAR:::.assignment_max(matrix(0, 2, 3)), 1:2)
  expect_identical(DMAR:::.assignment_max(matrix(0, 3, 2)), c(1L, 2L, NA))
  expect_error(DMAR:::.assignment_max(matrix(c(1, NA), 1)), "finite")
})

test_that("unequal factor counts leave the surplus factors without a partner", {
  # x has two factors, y three: two pairs and one row for y's surplus.
  x2 <- cong_Lx[, 1:2]
  res <- factor_congruence(x2, cong_Ly)
  expect_identical(nrow(res), 3L)
  expect_identical(res$factor_x, c("X1", "X2", NA))
  expect_identical(res$factor_y, c("Y2", "Y1", "Y3"))
  expect_identical(res$reflected, c(TRUE, FALSE, NA))
  expect_true(is.na(res$congruence[3]))
  C <- attr(res, "congruence_matrix")
  expect_identical(dim(C), c(2L, 3L))
  expect_equal(res$next_largest[3], max(abs(C[, "Y3"])))
  expect_identical(colnames(attr(res, "loadings_aligned")), c("Y2", "Y1", "Y3"))
  # Reversed: x has three, y two, and x's surplus has the row.
  rev <- factor_congruence(cong_Ly, x2)
  expect_identical(rev$factor_x, c("Y1", "Y2", "Y3"))
  expect_identical(rev$factor_y, c("X2", "X1", NA))
  expect_identical(rev$reflected, c(FALSE, TRUE, NA))
  expect_equal(rev$next_largest[3], max(abs(attr(rev, "congruence_matrix")["Y3", ])))
  expect_identical(generics::glance(rev)$n_matched, 2L)
})

test_that("match = \"order\" pairs by position, signed and unreflected", {
  res <- factor_congruence(cong_Lx, cong_Ly, match = "order")
  C <- attr(res, "congruence_matrix")
  expect_identical(res$factor_y, c("Y1", "Y2", "Y3"))
  expect_equal(res$congruence, diag(C))
  expect_true(any(res$congruence < 0))
  expect_identical(res$reflected, rep(FALSE, 3))
  expect_equal(attr(res, "loadings_aligned"), cong_Ly)
  expect_error(factor_congruence(cong_Lx, cong_Ly[, 1:2], match = "order"),
               "same number of factors")
})

test_that("a \"loadings\" object is read as its matrix", {
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  fa <- stats::factanal(holzinger_swineford[items], factors = 2)
  expect_s3_class(fa$loadings, "loadings")
  L <- unclass(fa$loadings)
  copy <- L[6:1, 2:1]
  copy[, 1] <- -copy[, 1]
  res <- factor_congruence(fa$loadings, copy)
  expect_equal(res, factor_congruence(L, copy))
  expect_equal(res$congruence, c(1, 1), tolerance = 1e-14)
  expect_identical(res$reflected, c(FALSE, TRUE))
  expect_identical(attr(res, "items"), items)
  expect_identical(attr(res, "inputs")$x, "matrix")
})

test_that("structure loadings need factor correlations", {
  expect_error(factor_congruence(cong_Lx, cong_Ly, loadings = "structure"),
               "factor correlations")
  expect_error(factor_congruence(cong_Lx, cong_Ly, match = "best"))
})

test_that("the refusals name what is wrong", {
  na <- cong_Lx
  na[3, 2] <- NA
  expect_error(factor_congruence(na, cong_Ly), "items v3")
  zero <- cong_Lx
  zero[, "X2"] <- 0
  expect_error(factor_congruence(cong_Lx, zero), "Factor X2 of 'y'")
  expect_error(factor_congruence(cong_Lx[1, , drop = FALSE],
                                 cong_Ly[1, , drop = FALSE]),
               "at least 2 items")
  other <- cong_Ly
  rownames(other)[c(2, 5)] <- c("w2", "w5")
  expect_error(factor_congruence(cong_Lx, other),
               "In 'x' and not in 'y': v2, v5. In 'y' and not in 'x': w2, w5.")
  expect_error(factor_congruence(as.data.frame(cong_Lx), cong_Ly),
               "data frame")
  expect_error(factor_congruence(list(cong_Lx), cong_Ly), "must be")
  expect_error(factor_congruence(cong_Lx, letters[1:8]), "must be")
  dup <- cong_Lx
  colnames(dup) <- c("A", "A", "B")
  expect_error(factor_congruence(dup, cong_Ly), "factor \\(column\\) names")
})

test_that("tidy() and glance() give the pairs and the summary", {
  res <- factor_congruence(cong_Lx, cong_Ly)
  td <- generics::tidy(res)
  expect_identical(names(td), c("term", "estimate", "reflected", "next_largest"))
  expect_identical(td$term, c("X1 | Y2", "X2 | Y1", "X3 | Y3"))
  expect_identical(td$estimate, res$congruence)
  expect_identical(class(td), "data.frame")
  gl <- generics::glance(res)
  expect_identical(nrow(gl), 1L)
  expect_identical(names(gl), c("n_items", "n_factors_x", "n_factors_y",
                                "n_matched", "match", "loadings",
                                "mean_congruence"))
  expect_identical(gl$n_items, 8L)
  expect_identical(gl$n_matched, 3L)
  expect_equal(gl$mean_congruence, mean(res$congruence))
  expect_identical(gl$match, "congruence")
  # A factor without a partner is named in its term, and the missing
  # partner is "(none)", not the text "NA".
  surplus <- generics::tidy(factor_congruence(cong_Lx[, 1:2], cong_Ly))
  expect_identical(surplus$term[3], "(none) | Y3")
  expect_true(is.na(surplus$estimate[3]))
  expect_true(is.na(surplus$reflected[3]))
  surplus_x <- generics::tidy(factor_congruence(cong_Lx, cong_Ly[, 1:2]))
  expect_identical(surplus_x$term[3], "X3 | (none)")
  expect_false(any(grepl("NA", c(surplus$term, surplus_x$term),
                         fixed = TRUE)))
})

test_that("the print lines state the comparison and apply no threshold", {
  res <- factor_congruence(cong_Lx[, 1:2], cong_Ly)
  foot <- DMAR:::.dmar_tbl_footer(res)
  expect_match(foot[["header"]][1],
               "Tucker's congruence coefficients: 8 items; x has 2 factors, y has 3",
               fixed = TRUE)
  expect_match(paste(foot[names(foot) == "matching"], collapse = " "),
               "maximizes the sum of the absolute coefficients; 2 pairs, 1 factor of y reflected")
  expect_match(paste(foot[names(foot) == "matching"], collapse = " "),
               "Without a partner: Y3 of y.", fixed = TRUE)
  table_lines <- foot[names(foot) == "congruence_matrix_table"]
  expect_length(table_lines, 3L)
  expect_match(table_lines[1], "Y1.*Y2.*Y3")
  text <- paste(foot, collapse = " ")
  expect_false(grepl("excellent|good|fair|poor|borderline|equal|similar|replicat|threshold|cutoff",
                     text, ignore.case = TRUE))
  out <- utils::capture.output(print(res))
  expect_true(any(grepl("Congruence matrix", out, fixed = TRUE)))
  ordered <- DMAR:::.dmar_tbl_footer(factor_congruence(cong_Lx, cong_Ly,
                                                       match = "order"))
  expect_match(paste(ordered, collapse = " "), "by position, column j of y")
  # A subset of the columns keeps the class but not the attributes, and
  # prints as an ordinary table.
  sub <- res[, c("factor_x", "factor_y", "congruence")]
  expect_s3_class(sub, "dmar_factor_congruence")
  expect_false(any(grepl("Congruence matrix", DMAR:::.dmar_tbl_footer(sub))))
  expect_no_error(utils::capture.output(print(sub)))
})

test_that("the table and the matrix beneath it print the coefficients alike", {
  # The help page's example: a copy with its items reversed, its factors
  # swapped, and one factor reflected. The matched coefficients are 1 and
  # 1 up to rounding error, and both are measured values, which print as
  # 1.000 in the table as they do in the matrix beneath it.
  L <- matrix(c(0.7, 0.6, 0.8, 0.1, 0.0, 0.2,
                0.1, 0.2, 0.0, 0.7, 0.8, 0.6), ncol = 2,
              dimnames = list(paste0("item_", 1:6), c("A", "B")))
  L_copy <- L[6:1, c("B", "A")]
  L_copy[, "A"] <- -L_copy[, "A"]
  matched <- factor_congruence(L, L_copy)
  expect_identical(format(matched)$congruence, c("1.000", "1.000"))
  out <- utils::capture.output(print(matched))
  expect_true(any(grepl("^ A +A +1\\.000 +TRUE +0\\.247 *$", out)))
  expect_true(any(grepl("^ B +B +1\\.000 +FALSE +0\\.247 *$", out)))
  expect_true(any(grepl("^B +1\\.000 +-0\\.247$", out)))
  # digits reaches the matrix as well as the table.
  out5 <- utils::capture.output(print(matched, digits = 5))
  expect_true(any(grepl("^ B +B +1\\.00000 +FALSE +0\\.24675 *$", out5)))
  expect_true(any(grepl("^A +0\\.24675 +-1\\.00000$", out5)))
  expect_false(any(grepl("0\\.247( |$)", out5)))
})

test_that("fa_exploratory(), lavaan, and matrix routes agree", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  efa <- fa_exploratory(holzinger_swineford, nfactors = 2, items = items,
                        seed = 113)
  cfa <- fa_confirmatory(holzinger_swineford,
                         list(spatial = items[1:3], verbal = items[4:6]),
                         output = "fit")
  std_efa <- lavaan::lavInspect(attr(efa, "fit"), "std")
  std_cfa <- lavaan::lavInspect(cfa, "std")
  by_fit <- factor_congruence(efa, cfa)
  by_matrix <- factor_congruence(unclass(std_efa$lambda),
                                 unclass(std_cfa$lambda))
  expect_equal(attr(by_fit, "congruence_matrix"),
               attr(by_matrix, "congruence_matrix"), tolerance = 1e-10)
  expect_identical(by_fit$factor_y, c("spatial", "verbal"))
  expect_identical(attr(by_fit, "inputs"),
                   list(x = "fa_exploratory", y = "lavaan", items = "name"))
  # The table and its lavaan fit are the same solution.
  same <- factor_congruence(efa, attr(efa, "fit"))
  expect_equal(same$congruence, c(1, 1), tolerance = 1e-10)
  # Phi_aligned is D Phi D in the matched order: y's factors are listed
  # the other way around and the first of them is reflected.
  L <- unclass(std_efa$lambda)
  res <- factor_congruence(cbind(a = -L[, "f2"], b = L[, "f1"]), efa)
  expect_identical(res$factor_y, c("f2", "f1"))
  expect_identical(res$reflected, c(TRUE, FALSE))
  Phi <- unclass(std_efa$psi)
  expect_equal(attr(res, "Phi_aligned"),
               Phi[2:1, 2:1] * outer(c(-1, 1), c(-1, 1)), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(attr(res, "loadings_aligned"), cbind(-L[, "f2"], L[, "f1"]),
               tolerance = 1e-10, ignore_attr = TRUE)
  # Structure loadings are the pattern times the factor correlations.
  st <- factor_congruence(efa, cfa, loadings = "structure")
  S_efa <- unclass(std_efa$lambda) %*% unclass(std_efa$psi)
  S_cfa <- unclass(std_cfa$lambda) %*% unclass(lavaan::lavInspect(cfa, "cor.lv"))
  expect_equal(unname(attr(st, "congruence_matrix")),
               unname(crossprod(S_efa, S_cfa) /
                        sqrt(outer(colSums(S_efa^2), colSums(S_cfa^2)))),
               tolerance = 1e-10)
  st_foot <- DMAR:::.dmar_tbl_footer(st)
  expect_match(st_foot[names(st_foot) == "header"][2],
               "standardized structure loadings", fixed = TRUE)
  # A subset of an fa_exploratory() table without its standardized rows.
  orig <- efa[efa$metric == "original", ]
  expect_error(factor_congruence(orig, cfa), "Pass the whole table")
})

test_that("fit inputs are refused where the comparison is not defined", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  factors <- list(spatial = items[1:3], verbal = items[4:6])
  expect_error(factor_congruence(fa_confirmatory(holzinger_swineford, factors),
                                 cong_Lx[1:6, 1:2]),
               "output = \"fit\"")
  model <- paste("spatial =~ t1_visual_perception + t2_cubes + t4_lozenges",
                 "verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
                 sep = "\n")
  multi <- lavaan::cfa(model, data = holzinger_swineford, group = "school")
  expect_error(factor_congruence(multi, multi), "multiple-group")
  # A higher-order factor has no item loadings.
  higher <- lavaan::cfa(paste(model,
    "speed =~ t10_addition + t12_counting_groups_of_dots + t13_straight_and_curved_capitals",
    "g =~ spatial + verbal + speed", sep = "\n"), data = holzinger_swineford)
  expect_error(factor_congruence(higher, higher), "higher-order factor")
  # An observed outcome of the structural model is not an item, and the
  # factor correlations are the model implied ones, not psi.
  structural <- lavaan::sem(paste(model, "t10_addition ~ spatial + verbal",
                                  sep = "\n"), data = holzinger_swineford)
  plain <- lavaan::cfa(model, data = holzinger_swineford)
  res <- factor_congruence(structural, plain)
  expect_identical(attr(res, "items"), items)
  expect_identical(res$factor_x, c("spatial", "verbal"))
  efa_list <- lavaan::efa(data = holzinger_swineford[items], nfactors = 1:2)
  expect_error(factor_congruence(efa_list, plain), "efaList")
  # A multilevel fit has loadings at each level; the refusal comes before
  # the fit is used, so an unfitted model shows it.
  two_level <- holzinger_swineford[items]
  two_level$cluster <- rep(1:30, length.out = nrow(two_level))
  multilevel <- lavaan::cfa(paste("level: 1", model, "level: 2", model,
                                  sep = "\n"),
                            data = two_level, cluster = "cluster",
                            do.fit = FALSE)
  expect_error(factor_congruence(plain, multilevel),
               "'y' is a multilevel lavaan fit")
  # Two iterations do not reach a solution; lavaan's own warning that the
  # optimizer stopped short is expected here and is not under test.
  stopped <- suppressWarnings(lavaan::cfa(model, data = holzinger_swineford,
                                          control = list(iter.max = 2)))
  expect_false(lavaan::lavInspect(stopped, "converged"))
  expect_error(factor_congruence(stopped, plain),
               "The lavaan fit 'x' did not converge")
})

test_that("a CFA results table is refused, naming the three functions", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  one <- cfa_1(holzinger_swineford, items = items[1:3])
  expect_error(factor_congruence(one, cong_Lx[1:3, 1, drop = FALSE]),
               "fa_confirmatory\\(\\), cfa_1\\(\\), or cfa_2\\(\\)")
  expect_error(factor_congruence(cong_Lx[1:3, 1, drop = FALSE], one),
               "'y' is a table of results.*output = \"fit\"")
  std <- cfa_2(holzinger_swineford, factor_1 = items[1:3],
               factor_2 = items[4:6], output = "standardized")
  expect_s3_class(std, "lavaan.data.frame")
  expect_error(factor_congruence(std, cong_Lx[1:6, 1:2]),
               "lavaan table of parameter estimates.*output = \"fit\"")
})

test_that("Phi_aligned is D Phi D under structure loadings with a reflection", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  model <- paste("spatial =~ t1_visual_perception + t2_cubes + t4_lozenges",
                 "verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
                 sep = "\n")
  plain <- lavaan::cfa(model, data = holzinger_swineford)
  # Fixing the marker loading at -1 defines the reflected spatial factor:
  # its loadings and its correlation with verbal change sign, and nothing
  # else in the model changes.
  flipped <- lavaan::cfa(sub("t1_visual_perception",
                             "-1*t1_visual_perception", model, fixed = TRUE),
                         data = holzinger_swineford)
  Phi_y <- unclass(lavaan::lavInspect(flipped, "cor.lv"))
  expect_true(Phi_y["verbal", "spatial"] < 0)
  st <- factor_congruence(plain, flipped, loadings = "structure")
  expect_identical(st$factor_y, c("spatial", "verbal"))
  expect_identical(st$reflected, c(TRUE, FALSE))
  expect_equal(st$congruence, c(1, 1), tolerance = 1e-6)
  D <- diag(c(-1, 1))
  expect_equal(attr(st, "Phi_aligned"), D %*% Phi_y %*% D,
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(attr(st, "Phi_aligned"),
               unclass(lavaan::lavInspect(plain, "cor.lv")),
               tolerance = 1e-5, ignore_attr = TRUE)
  # The aligned loadings are y's structure loadings, reflected to match.
  S_y <- unclass(lavaan::lavInspect(flipped, "std")$lambda) %*% Phi_y
  expect_equal(attr(st, "loadings_aligned"), S_y %*% D,
               tolerance = 1e-14, ignore_attr = TRUE)
})
