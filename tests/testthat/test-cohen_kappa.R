test_that("cohen_kappa() reproduces Cohen (1960) Table 1 (kappa ~ .49)", {
  levels_3 <- c("Personality", "Neurosis", "Brain")
  M <- matrix(c(88, 14, 18,
                10, 40, 10,
                 2,  6, 12),
              nrow = 3, byrow = TRUE,
              dimnames = list(levels_3, levels_3))
  df <- as.data.frame(as.table(M))
  r1 <- factor(rep(df$Var1, df$Freq), levels = levels_3)
  r2 <- factor(rep(df$Var2, df$Freq), levels = levels_3)
  res <- cohen_kappa(r1, r2)
  expect_equal(round(res$kappa, 2), 0.49)
  expect_equal(res$n, 200L)
  expect_equal(res$n_categories, 3L)
})

test_that("cohen_kappa() with quadratic weights matches a hand calculation", {
  # 4 categories on the diagonal -> perfect agreement -> kappa = 1.
  perfect <- factor(c("A", "B", "C", "D", "A", "B", "C", "D"),
                    levels = c("A", "B", "C", "D"))
  res <- cohen_kappa(perfect, perfect, weights = "quadratic")
  expect_equal(res$kappa, 1)
})

test_that("cohen_kappa() unweighted equals Cohen formulation; weighted shifts kappa", {
  # Off-by-one ordinal disagreements -> unweighted underestimates agreement.
  set.seed(113)
  x <- sample(1:5, 100, replace = TRUE)
  y <- pmin(pmax(x + sample(-1:1, 100, replace = TRUE), 1), 5)
  un <- cohen_kappa(x, y, weights = "unweighted")$kappa
  qd <- cohen_kappa(x, y, weights = "quadratic")$kappa
  ln <- cohen_kappa(x, y, weights = "linear")$kappa
  expect_lt(un, qd)   # weighting credits near-misses, raising kappa
  expect_lt(un, ln)
  expect_lt(ln, qd)   # quadratic > linear (more lenient on near-misses)
})

test_that("cohen_kappa() handles a custom weight matrix", {
  W <- matrix(c(1.0, 0.5, 0.0,
                0.5, 1.0, 0.5,
                0.0, 0.5, 1.0),
              nrow = 3, byrow = TRUE)
  set.seed(113)
  r1 <- sample(1:3, 50, replace = TRUE)
  r2 <- sample(1:3, 50, replace = TRUE)
  res <- cohen_kappa(r1, r2, weights = W)
  expect_equal(res$weights, "custom")
})

test_that("cohen_kappa() drops NA pairs before computing", {
  x <- c(1, 2, 3, 4, NA, 1, 2)
  y <- c(1, 2, 3, NA, 5, 1, 2)
  res <- cohen_kappa(x, y)
  expect_equal(res$n, 5L)  # 5 complete pairs
})

test_that("cohen_kappa() rejects bad inputs", {
  expect_error(cohen_kappa(1:5, 1:6),                 "same length")
  expect_error(cohen_kappa(1:5, 1:5, conf_level = 1.5), "conf_level")
  expect_error(cohen_kappa(1:5, 1:5, weights = matrix(0, 2, 2)),
               "must be 5 x 5")
})

test_that("cohen_kappa() requires a complete 'categories' set (HIGH-08)", {
  # Independent oracle: irr::kappa2 on the full data. A 'categories' set that
  # omits an observed rating used to silently drop those ratings yet still
  # divide by the original N, returning a corrupted kappa instead of erroring.
  r1 <- c("a", "a", "b", "b", "c", "c", "a", "b")
  r2 <- c("a", "b", "b", "b", "c", "a", "a", "c")

  # Pinned from irr::kappa2 (irr 0.85, 2026-08-09); live comparison in
  # tools/oracle_checks.R.
  oracle <- 0.4285714285714285

  # Complete set (default, and an explicit reordered set) matches the oracle.
  expect_equal(cohen_kappa(r1, r2)$kappa, oracle, tolerance = 1e-7)
  expect_equal(cohen_kappa(r1, r2, categories = c("c", "b", "a"))$kappa,
               oracle, tolerance = 1e-7)

  # Omitting an observed category now errors (listing the omitted label)
  # instead of returning a corrupted value.
  expect_error(cohen_kappa(r1, r2, categories = c("a", "b")),
               "omits observed rating")
  # A duplicated or missing category is also rejected.
  expect_error(cohen_kappa(r1, r2, categories = c("a", "b", "c", "c")),
               "must be unique")
  expect_error(cohen_kappa(r1, r2, categories = c("a", "b", NA)),
               "must not contain missing")
})

test_that("fleiss_kappa() reproduces Fleiss (1971) Table 1 (kappa ~ .43)", {
  fleiss_1971 <- matrix(c(
    0,0,0,0,6, 0,3,0,0,3, 0,1,4,0,1, 0,0,0,0,6, 0,3,0,3,0,
    2,0,4,0,0, 0,0,4,0,2, 2,0,3,1,0, 2,0,0,4,0, 0,0,0,0,6,
    1,0,0,5,0, 1,1,0,1,3, 0,3,3,0,0, 1,0,0,5,0, 0,2,0,3,1,
    0,0,5,0,1, 3,0,0,1,2, 5,1,0,0,0, 0,2,0,4,0, 1,0,2,0,3,
    0,0,0,0,6, 0,1,0,5,0, 0,2,0,1,3, 2,0,0,4,0, 1,0,0,4,1,
    0,5,0,1,0, 4,0,0,0,2, 0,2,0,4,0, 1,0,5,0,0, 0,0,0,0,6
  ), nrow = 30, byrow = TRUE)
  res <- fleiss_kappa(fleiss_1971)
  # Fleiss (1971) Table 2 reports P_bar = 0.5475 and P_e = 0.2179, giving
  # kappa = (0.5475 - 0.2179) / (1 - 0.2179) = 0.4214; the paper itself
  # rounds this to 0.43, but the precise value is 0.42 to two decimals.
  expect_equal(round(res$kappa, 2), 0.42)
  expect_equal(res$n_subjects, 30L)
  expect_equal(res$n_raters, 6L)
  expect_equal(res$n_categories, 5L)
})

test_that("fleiss_kappa() = 1 when raters fully agree on every subject", {
  N <- 10; m <- 4; k <- 3
  perfect <- matrix(0, nrow = N, ncol = k)
  for (i in seq_len(N)) {
    perfect[i, sample(k, 1)] <- m
  }
  res <- fleiss_kappa(perfect)
  expect_equal(res$kappa, 1)
})

test_that("fleiss_kappa() rejects unequal row sums", {
  bad <- rbind(c(2, 2, 0), c(1, 2, 0))   # row sums 4 and 3
  expect_error(fleiss_kappa(bad), "same value")
})

test_that("fleiss_kappa() rejects fewer than 2 raters", {
  one_rater <- matrix(c(1, 0, 0,  0, 1, 0,  0, 0, 1), nrow = 3, byrow = TRUE)
  expect_error(fleiss_kappa(one_rater), "at least 2 raters")
})


# ---------------------------------------------------------------------------
# Cohen (1968) Table 1 replication: table input and disagreement scaling
# ---------------------------------------------------------------------------

cohen_1968 <- matrix(
  c(88, 10,  2,
    14, 40,  6,
    18, 10, 12),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("Personality disorder", "Neurosis", "Psychosis"),
                  c("Personality disorder", "Neurosis", "Psychosis")))
v_1968 <- matrix(c(0, 1, 3,
                   1, 0, 6,
                   3, 6, 0), nrow = 3, byrow = TRUE)

test_that("cohen_kappa(table =) reproduces Cohen (1968) unweighted kappa = .492", {
  res <- cohen_kappa(table = cohen_1968)
  expect_equal(res$kappa, 1 - 0.30 / 0.59, tolerance = 1e-12)
  expect_equal(round(res$kappa, 3), 0.492)
  expect_equal(res$n, 200)
  expect_equal(res$n_categories, 3L)
})

test_that("cohen_kappa() reproduces Cohen (1968) weighted kappa = .348", {
  res <- cohen_kappa(table = cohen_1968, weights = v_1968,
                     weight_scaling = "disagreement")
  # Cohen's Formula 8: kappa_w = 1 - q'_o / q'_c = 1 - .90 / 1.38.
  expect_equal(res$kappa, 1 - 0.90 / 1.38, tolerance = 1e-12)
  expect_equal(round(res$kappa, 3), 0.348)
})

test_that("cohen_kappa() reproduces Cohen's 6 <-> 1 weight interchange, kappa_w = .574", {
  v_swap <- matrix(c(0, 6, 3,
                     6, 0, 1,
                     3, 1, 0), nrow = 3, byrow = TRUE)
  res <- cohen_kappa(table = cohen_1968, weights = v_swap,
                     weight_scaling = "disagreement")
  expect_equal(res$kappa, 1 - 1.10 / 2.58, tolerance = 1e-12)
  expect_equal(round(res$kappa, 3), 0.574)
})

test_that("cohen_kappa() reproduces Cohen's (1968) asymmetric validity example, kappa_w = .353", {
  v_validity <- matrix(c(0, 1, 2,
                         1, 0, 2,
                         4, 6, 0), nrow = 3, byrow = TRUE)
  res <- cohen_kappa(table = cohen_1968, weights = v_validity,
                     weight_scaling = "disagreement")
  # Cohen reports sum(v * P) = .86 and sum(v * E) = 1.33.
  expect_equal(res$kappa, 1 - 0.86 / 1.33, tolerance = 1e-12)
  expect_equal(round(res$kappa, 3), 0.353)
})

test_that("Cohen's (1968) Formulas 10 and 13 arithmetic matches his printed values", {
  # These verify the Table 1 quantities used in the help file example;
  # the function itself reports the Fleiss-Cohen-Everitt (1969) SE.
  P <- cohen_1968 / sum(cohen_1968)
  E <- outer(rowSums(P), colSums(P))
  q_o <- sum(v_1968 * P)
  q_c <- sum(v_1968 * E)
  expect_equal(q_o, 0.90, tolerance = 1e-12)
  expect_equal(q_c, 1.38, tolerance = 1e-12)
  expect_equal(sum(v_1968^2 * P), 3.90, tolerance = 1e-12)
  expect_equal(sum(v_1968^2 * E), 5.10, tolerance = 1e-12)
  se_10 <- sqrt((3.90 - 0.90^2) / (200 * 1.38^2))
  se_13 <- sqrt((5.10 - 1.38^2) / (200 * 1.38^2))
  expect_equal(round(se_10, 4), 0.0901)
  expect_equal(round(se_13, 4), 0.0916)
  expect_equal(round((1 - q_o / q_c) / se_13, 2), 3.80)
})

test_that("disagreement scaling equals the manual agreement conversion", {
  w_agree <- 1 - v_1968 / max(v_1968)
  res_d <- cohen_kappa(table = cohen_1968, weights = v_1968,
                       weight_scaling = "disagreement")
  res_a <- cohen_kappa(table = cohen_1968, weights = w_agree)
  expect_equal(res_d$kappa, res_a$kappa, tolerance = 1e-12)
  expect_equal(res_d$se, res_a$se, tolerance = 1e-12)
  expect_equal(res_d$weights, "custom_disagreement")
  expect_equal(res_a$weights, "custom")
})

test_that("kappa_w is invariant to positive rescaling of disagreement weights", {
  res_1 <- cohen_kappa(table = cohen_1968, weights = v_1968,
                       weight_scaling = "disagreement")
  res_2 <- cohen_kappa(table = cohen_1968, weights = 10 * v_1968,
                       weight_scaling = "disagreement")
  expect_equal(res_1$kappa, res_2$kappa, tolerance = 1e-12)
  expect_equal(res_1$se, res_2$se, tolerance = 1e-12)
})

test_that("table input equals expanded rater-vector input", {
  labs <- rownames(cohen_1968)
  df <- as.data.frame(as.table(cohen_1968))
  r1 <- factor(rep(df$Var1, df$Freq), levels = labs)
  r2 <- factor(rep(df$Var2, df$Freq), levels = labs)
  res_tab <- cohen_kappa(table = cohen_1968, weights = v_1968,
                         weight_scaling = "disagreement")
  res_vec <- cohen_kappa(r1, r2, weights = v_1968,
                         weight_scaling = "disagreement")
  expect_equal(res_tab$kappa, res_vec$kappa, tolerance = 1e-12)
  expect_equal(res_tab$se, res_vec$se, tolerance = 1e-12)
  expect_equal(res_tab$n, res_vec$n)
})

test_that("cohen_kappa() agrees with psych::cohen.kappa on Cohen (1968) Table 1", {
  # Pinned from psych::cohen.kappa (psych 2.5.3, 2026-08-09); psych reports
  # the FCE interval, so each pinned se backs out the CI half-width divided by
  # qnorm(0.975). Live comparison in tools/oracle_checks.R.
  psych_kappa_un <- 0.4915254237288135
  psych_se_un    <- 0.05100181557607786
  psych_kappa_w  <- 0.3478260869565211
  psych_se_w     <- 0.07550401525482291

  res_un <- cohen_kappa(table = cohen_1968)
  expect_equal(res_un$kappa, psych_kappa_un, tolerance = 1e-8)
  expect_equal(res_un$se, psych_se_un, tolerance = 1e-6)

  res_w <- cohen_kappa(table = cohen_1968, weights = v_1968,
                       weight_scaling = "disagreement")
  expect_equal(res_w$kappa, psych_kappa_w, tolerance = 1e-8)
  expect_equal(res_w$se, psych_se_w, tolerance = 1e-6)
})

test_that("cohen_kappa() input validation for tables and disagreement weights", {
  expect_error(cohen_kappa(1:5, 1:5, table = cohen_1968),
               regexp = "not both")
  expect_error(cohen_kappa(), regexp = "Supply either")
  expect_error(cohen_kappa(table = cohen_1968[, 1:2]), regexp = "square")
  bad_names <- cohen_1968
  colnames(bad_names) <- c("a", "b", "c")
  expect_error(cohen_kappa(table = bad_names), regexp = "identical")
  expect_error(
    cohen_kappa(table = cohen_1968, categories = c("x", "y", "z")),
    regexp = "conflicts")
  # Disagreement weights must be a zero-diagonal, non-negative matrix.
  expect_error(
    cohen_kappa(table = cohen_1968, weights = "linear",
                weight_scaling = "disagreement"),
    regexp = "requires a custom")
  v_bad <- v_1968; diag(v_bad) <- 1
  expect_error(
    cohen_kappa(table = cohen_1968, weights = v_bad,
                weight_scaling = "disagreement"),
    regexp = "0 on the diagonal")
  v_neg <- v_1968; v_neg[1, 2] <- -1
  expect_error(
    cohen_kappa(table = cohen_1968, weights = v_neg,
                weight_scaling = "disagreement"),
    regexp = "non-negative")
  # An agreement-scaled call with Cohen-style weights warns.
  expect_warning(
    cohen_kappa(table = cohen_1968, weights = v_1968),
    regexp = "disagreement")
})

test_that("unnamed table input uses integer category labels", {
  tab <- unname(cohen_1968)
  res <- cohen_kappa(table = tab)
  expect_equal(res$n_categories, 3L)
  expect_equal(res$kappa, 1 - 0.30 / 0.59, tolerance = 1e-12)
})


test_that("factor inputs keep their own level order for weighting", {
  # Levels in substantive (non-alphabetical) order; a custom weight
  # matrix must align with that order, not with sort() order.
  labs <- c("Personality disorder", "Neurosis", "Psychosis")
  df <- as.data.frame(as.table(cohen_1968))
  r1 <- factor(rep(df$Var1, df$Freq), levels = labs)
  r2 <- factor(rep(df$Var2, df$Freq), levels = labs)
  res_vec <- cohen_kappa(r1, r2, weights = v_1968,
                         weight_scaling = "disagreement")
  expect_equal(res_vec$kappa, 1 - 0.90 / 1.38, tolerance = 1e-12)
})


# ---------------------------------------------------------------------------
# The diagnosis_agreement data set and the cells attribute
# ---------------------------------------------------------------------------

test_that("diagnosis_agreement reproduces Cohen's Table 1 exactly", {
  data(diagnosis_agreement, envir = environment())
  expect_equal(nrow(diagnosis_agreement), 9L)
  expect_equal(sum(diagnosis_agreement$frequency), 200)
  tab <- xtabs(frequency ~ judge_b + judge_a, data = diagnosis_agreement)
  # Judge B margins .50/.30/.20 (rows), Judge A margins .60/.30/.10 (columns).
  expect_equal(unname(rowSums(tab)), c(100, 60, 40))
  expect_equal(unname(colSums(tab)), c(120, 60, 20))
  expect_equal(diagnosis_agreement$observed_proportion,
               diagnosis_agreement$frequency / 200)
  # Chance-expected proportions, Cohen's parenthetical values.
  expect_equal(diagnosis_agreement$expected_proportion,
               c(.30, .15, .05, .18, .09, .03, .12, .06, .02))
  v <- unclass(xtabs(disagreement_weight ~ judge_b + judge_a,
                     data = diagnosis_agreement))
  expect_true(all(v == t(v)))
  expect_true(all(diag(v) == 0))
  expect_equal(cohen_kappa(table = tab)$kappa, 1 - 0.30 / 0.59,
               tolerance = 1e-12)
  expect_equal(
    cohen_kappa(table = tab, weights = v,
                weight_scaling = "disagreement")$kappa,
    1 - 0.90 / 1.38, tolerance = 1e-12)
})

test_that("cohen_kappa() attaches the per-cell detail as the cells attribute", {
  data(diagnosis_agreement, envir = environment())
  tab <- xtabs(frequency ~ judge_b + judge_a, data = diagnosis_agreement)
  v <- unclass(xtabs(disagreement_weight ~ judge_b + judge_a,
                     data = diagnosis_agreement))
  res <- cohen_kappa(table = tab, weights = v,
                     weight_scaling = "disagreement")
  cells <- attr(res, "cells")
  expect_s3_class(cells, "data.frame")
  expect_named(cells, c("rater_1", "rater_2", "observed_proportion",
                        "expected_proportion", "weight",
                        "disagreement_weight"))
  # Row-major over the table, so the rows align with the data set.
  expect_equal(cells$observed_proportion,
               diagnosis_agreement$observed_proportion)
  expect_equal(cells$expected_proportion,
               diagnosis_agreement$expected_proportion)
  expect_equal(cells$disagreement_weight,
               diagnosis_agreement$disagreement_weight)
  expect_equal(cells$weight,
               1 - diagnosis_agreement$disagreement_weight / 6)
  expect_equal(as.character(cells$rater_1),
               as.character(diagnosis_agreement$judge_b))
  expect_equal(as.character(cells$rater_2),
               as.character(diagnosis_agreement$judge_a))

  # Unweighted: identity weights and no disagreement_weight column.
  cells0 <- attr(cohen_kappa(table = tab), "cells")
  expect_false("disagreement_weight" %in% names(cells0))
  expect_equal(cells0$weight,
               as.numeric(cells0$rater_1 == cells0$rater_2))

  # Vector input carries the attribute too.
  labs <- levels(diagnosis_agreement$judge_b)
  df <- as.data.frame(as.table(tab))
  r1 <- factor(rep(df$judge_b, df$Freq), levels = labs)
  r2 <- factor(rep(df$judge_a, df$Freq), levels = labs)
  cells_vec <- attr(cohen_kappa(r1, r2), "cells")
  expect_equal(cells_vec$observed_proportion,
               cells0$observed_proportion)
})

test_that("validity orientation reproduces Cohen's published .353; the printed transpose does not", {
  data(diagnosis_agreement, envir = environment())
  tab <- xtabs(frequency ~ judge_b + judge_a, data = diagnosis_agreement)
  v_validity <- matrix(c(0, 1, 2,
                         1, 0, 2,
                         4, 6, 0), nrow = 3, byrow = TRUE)
  expect_equal(
    cohen_kappa(table = tab, weights = v_validity,
                weight_scaling = "disagreement")$kappa,
    1 - 0.86 / 1.33, tolerance = 1e-12)
  # The transpose corresponds to the paper's printed weight display read
  # with its own labels; it yields a different value than the published
  # ones. The weighted kappa vignette works both computations in full.
  expect_equal(
    cohen_kappa(table = tab, weights = t(v_validity),
                weight_scaling = "disagreement")$kappa,
    1 - 0.62 / 1.07, tolerance = 1e-12)
})

test_that("bootstrap intervals: reproducible, ordered, near Wald at moderate N", {
  set.seed(113)
  n <- 200
  r1 <- sample(letters[1:3], n, replace = TRUE)
  r2 <- ifelse(runif(n) < 0.7, r1, sample(letters[1:3], n, replace = TRUE))

  w  <- cohen_kappa(r1, r2)
  p1 <- cohen_kappa(r1, r2, ci_method = "percentile", B = 2000, seed = 113)
  p2 <- cohen_kappa(r1, r2, ci_method = "percentile", B = 2000, seed = 113)
  b1 <- cohen_kappa(r1, r2, ci_method = "bca", B = 2000, seed = 113)

  # Same seed, same interval; the point estimate and the asymptotic se
  # are untouched by the choice of interval.
  expect_identical(as.data.frame(p1), as.data.frame(p2))
  expect_equal(p1$kappa, w$kappa)
  expect_equal(p1$se, w$se)

  # Ordered and bracketing the estimate.
  for (r in list(p1, b1)) {
    expect_lt(r$lower_limit, r$kappa)
    expect_gt(r$upper_limit, r$kappa)
  }

  # At this N the bootstrap and Wald intervals agree closely.
  expect_equal(p1$lower_limit, w$lower_limit, tolerance = 0.05)
  expect_equal(p1$upper_limit, w$upper_limit, tolerance = 0.05)

  expect_identical(attr(p1, "ci_method"), "percentile")
  expect_identical(attr(p1, "B_used"), 2000L)
  expect_identical(attr(w, "ci_method"), "wald")
  expect_null(attr(w, "B_used"))
})

test_that("table input bootstraps as the equivalent paired ratings", {
  set.seed(113)
  n <- 150
  r1 <- sample(1:3, n, replace = TRUE)
  r2 <- ifelse(runif(n) < 0.6, r1, sample(1:3, n, replace = TRUE))
  raw <- cohen_kappa(r1, r2, ci_method = "percentile", B = 2000, seed = 113)
  tab <- cohen_kappa(table = table(r1, r2), ci_method = "percentile",
                     B = 2000, seed = 113)
  # The expansion reorders subjects, so the same seed draws different
  # resamples; the intervals agree up to Monte Carlo error only.
  expect_equal(raw$kappa, tab$kappa)
  expect_equal(raw$lower_limit, tab$lower_limit, tolerance = 0.03)
  expect_equal(raw$upper_limit, tab$upper_limit, tolerance = 0.03)
})

test_that("bootstrap input validation", {
  expect_error(cohen_kappa(c("a", "b"), c("a", "b"),
                           ci_method = "percentile", B = 50),
               "at least 100")
})
