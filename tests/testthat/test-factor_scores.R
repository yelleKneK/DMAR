# Tests for factor_scores(). The oracle throughout is lavaan::lavPredict()
# itself, aligned by lavaan's case indices, together with the closed forms
# the scores and standard errors must satisfy pattern by pattern and the
# determinacies lavaan reports for complete data. The fixtures are the
# three-factor model of nine Holzinger and Swineford tests, with missing
# values made to depend on age (a mechanism outside the model), two rows
# with every item missing, and two rows missing every textual item.

fs_items <- list(
  visual  = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
  textual = c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning"),
  speed   = c("t10_addition", "t12_counting_groups_of_dots",
              "t13_straight_and_curved_capitals"))
fs_model <- paste(vapply(names(fs_items), function(f) {
  paste(f, "=~", paste(fs_items[[f]], collapse = " + "))
}, character(1)), collapse = "\n")

# The columns of the one canonical table, with their types.
fs_shape <- c(row = "integer", group = "character", factor = "character",
              .derived_regression = "numeric",
              .derived_regression_se = "numeric",
              .derived_regression_ci_lower = "numeric",
              .derived_regression_ci_upper = "numeric",
              .derived_bartlett = "numeric", .derived_bartlett_se = "numeric",
              .derived_bartlett_ci_lower = "numeric",
              .derived_bartlett_ci_upper = "numeric",
              items_observed = "integer", items = "integer",
              pattern = "integer", determinacy_regression = "numeric",
              determinacy_bartlett = "numeric", reason = "character")
expect_fs_shape <- function(x) {
  expect_identical(names(x), names(fs_shape))
  expect_identical(vapply(x, function(v) class(v)[1L], character(1)),
                   fs_shape)
}

fs_mar_data <- function(extra = character(0)) {
  hs <- holzinger_swineford
  d <- hs[, c(unlist(fs_items), extra)]
  set.seed(113)
  pm <- stats::plogis(-2 + 0.5 * scale(hs$age_years)[, 1])
  for (v in unlist(fs_items)[c(2, 3, 5, 6, 8, 9)]) {
    d[[v]][stats::runif(nrow(d)) < pm] <- NA
  }
  d[c(5, 77), unlist(fs_items)] <- NA
  d
}

fs_fit_fiml <- function(d, ...) {
  # lavaan warns that the two rows with every item missing are ignored.
  suppressWarnings(lavaan::cfa(fs_model, data = d, std.lv = TRUE,
                               missing = "ml", ...))
}

test_that("under FIML both scores equal lavPredict() row for row", {
  skip_if_not_installed("lavaan")
  d <- fs_mar_data()
  d[c(10, 11), fs_items$textual] <- NA
  fit <- fs_fit_fiml(d)
  fs <- factor_scores(fit, data = d)
  expect_s3_class(fs, "dmar_factor_scores")
  expect_s3_class(fs, "dmar_tbl")
  expect_fs_shape(fs)
  expect_identical(nrow(fs), 3L * nrow(d))
  expect_identical(fs$row, rep(seq_len(nrow(d)), 3L))
  expect_identical(fs$factor, rep(names(fs_items), each = nrow(d)))
  idx <- lavaan::lavInspect(fit, "case.idx")
  for (m in c("regression", "Bartlett")) {
    lp <- lavaan::lavPredict(fit, method = m, se = "standard")
    se <- attr(lp, "se")[[1]]
    col <- paste0(".derived_", tolower(m))
    for (f in names(fs_items)) {
      mine <- fs[[col]][fs$factor == f][idx]
      expect_identical(mine, unname(lp[, f]))
      expect_identical(fs[[paste0(col, "_se")]][fs$factor == f][idx],
                       unname(se[, f]))
    }
  }
  expect_identical(attr(fs, "methods"),
                   c(regression = "regression", bartlett = "Bartlett"))
  expect_identical(attr(fs, "scored_from"), "fit")
  expect_identical(attr(fs, "missing"), "ml")
  expect_identical(attr(fs, "conf_level"), 0.95)

  # The interval is the score plus and minus the normal quantile times
  # the standard error.
  z <- stats::qnorm(0.975)
  expect_equal(fs$.derived_regression_ci_lower, fs$.derived_regression - z * fs$.derived_regression_se)
  expect_equal(fs$.derived_bartlett_ci_upper, fs$.derived_bartlett + z * fs$.derived_bartlett_se)
  fs90 <- factor_scores(fit, data = d, conf_level = 0.90)
  expect_equal(fs90$.derived_bartlett_ci_upper - fs90$.derived_bartlett,
               stats::qnorm(0.95) * fs$.derived_bartlett_se)
  expect_identical(fs90$.derived_bartlett, fs$.derived_bartlett)
})

test_that("each case is scored from its own observed items", {
  skip_if_not_installed("lavaan")
  d <- fs_mar_data()
  fit <- fs_fit_fiml(d)
  fs <- factor_scores(fit, data = d)
  est <- lavaan::lavInspect(fit, "est")
  lam <- est$lambda
  th <- est$theta
  phi <- est$psi
  imp <- lavaan::lavInspect(fit, "implied")
  sig <- imp$cov
  mu <- imp$mean
  x <- as.matrix(d[, rownames(lam)])
  lv <- colnames(lam)
  got <- function(col) {
    sapply(lv, function(f) fs[[col]][fs$factor == f])
  }
  rows <- which(rowSums(!is.na(x)) > 0L)
  want <- lapply(rows, function(i) {
    o <- !is.na(x[i, ])
    yc <- x[i, o] - mu[o]
    lo <- lam[o, , drop = FALSE]
    s_oo <- sig[o, o]
    post <- phi - phi %*% t(lo) %*% solve(s_oo) %*% lo %*% phi
    info <- t(lo) %*% solve(th[o, o]) %*% lo
    list(reg = drop(phi %*% t(lo) %*% solve(s_oo, yc)),
         reg_se = sqrt(diag(post)),
         bar = drop(solve(info, t(lo) %*% solve(th[o, o], yc))),
         bar_se = sqrt(diag(solve(info))))
  })
  stack <- function(el) unname(do.call(rbind, lapply(want, `[[`, el)))
  expect_equal(unname(got(".derived_regression")[rows, ]), stack("reg"),
               tolerance = 1e-10)
  expect_equal(unname(got(".derived_regression_se")[rows, ]), stack("reg_se"),
               tolerance = 1e-10)
  expect_equal(unname(got(".derived_bartlett")[rows, ]), stack("bar"),
               tolerance = 1e-10)
  expect_equal(unname(got(".derived_bartlett_se")[rows, ]), stack("bar_se"),
               tolerance = 1e-10)
})

test_that("per-pattern determinacy agrees with lavaan and the errors", {
  skip_if_not_installed("lavaan")
  d <- fs_mar_data()
  fit <- fs_fit_fiml(d)
  fs <- factor_scores(fit, data = d)
  det <- attr(fs, "determinacy")
  expect_equal(det$determinacy_regression,
               unname(lavaan::lavInspect(fit, "fs.determinacy")),
               tolerance = 1e-10)
  expect_equal(det$determinacy_bartlett,
               unname(lavaan::lavInspect(fit, "fs.determinacy.bartlett")),
               tolerance = 1e-10)
  # Pattern 1 is the complete one, and its rows carry those values.
  pt <- attr(fs, "patterns")
  expect_identical(unique(pt$items_observed[pt$pattern == 1L]), 3L)
  for (f in names(fs_items)) {
    one <- fs$pattern == 1L & fs$factor == f
    expect_equal(unique(fs$determinacy_regression[one]),
                 det$determinacy_regression[det$factor == f])
  }
  # With every factor variance fixed at 1, the regression determinacy is
  # sqrt(1 - posterior variance) and the Bartlett determinacy is
  # sqrt(1 / (1 + error variance)): an independent check of the
  # per-pattern values against lavaan's per-case standard errors.
  ok <- is.finite(fs$.derived_regression_se) & is.finite(fs$determinacy_regression)
  expect_equal(fs$determinacy_regression[ok],
               sqrt(1 - fs$.derived_regression_se[ok]^2), tolerance = 1e-10)
  ok <- is.finite(fs$.derived_bartlett_se) & is.finite(fs$determinacy_bartlett)
  expect_equal(fs$determinacy_bartlett[ok],
               sqrt(1 / (1 + fs$.derived_bartlett_se[ok]^2)), tolerance = 1e-10)
  # The regression score is the more determinate, and a missing item
  # costs determinacy.
  expect_true(all(fs$determinacy_regression >= fs$determinacy_bartlett - 1e-12,
                  na.rm = TRUE))
  tex <- pt[pt$factor == "textual" & pt$items_observed > 0L, ]
  expect_true(all(diff(tapply(tex$determinacy_bartlett, tex$items_observed,
                              mean)) > 0))
  # The pattern table's standard errors are those in the rows.
  m <- merge(fs[fs$items_observed > 0, c("pattern", "factor", 
                  ".derived_bartlett_se")],
             pt[, c("pattern", "factor", "bartlett_se")],
             by = c("pattern", "factor"))
  expect_equal(m$.derived_bartlett_se, m$bartlett_se)
  expect_identical(sum(pt$n[pt$factor == "visual"]), nrow(d))
})

test_that("for one factor both determinacies squared equal coefficient H", {
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa("verbal =~ t5_general_information +
                      t6_paragraph_comprehension + t7_sentence +
                      t8_word_classification + t9_word_meaning",
                     data = holzinger_swineford, std.lv = TRUE)
  fs <- factor_scores(fit)
  # Coefficient H from the standardized loadings, by its closed form
  # (Hancock & Mueller, 2001).
  std <- lavaan::standardizedSolution(fit)
  l <- std$est.std[std$op == "=~"]
  q <- sum(l^2 / (1 - l^2))
  h <- q / (1 + q)
  det <- attr(fs, "determinacy")
  expect_equal(det$determinacy_regression^2, h, tolerance = 1e-10)
  expect_equal(det$determinacy_bartlett^2, h, tolerance = 1e-10)
  # Complete data: both scores are lavPredict()'s, in one pattern, with
  # one standard error for every row.
  for (m in c("regression", "Bartlett")) {
    expect_identical(fs[[paste0(".derived_", tolower(m))]],
                     unname(lavaan::lavPredict(fit, method = m)[, "verbal"]))
  }
  expect_identical(unique(fs$pattern), 1L)
  expect_length(unique(fs$.derived_bartlett_se), 1L)
  expect_true(all(is.na(fs$group)))
  expect_true(all(fs$reason == ""))
})

test_that("a row with no observed item and a factor with none of its items", {
  skip_if_not_installed("lavaan")
  d <- fs_mar_data()
  d[c(10, 11), fs_items$textual] <- NA
  fit <- fs_fit_fiml(d)
  fs <- factor_scores(fit, data = d)
  empty <- fs$row %in% c(5, 77)
  expect_true(all(is.na(fs$.derived_regression[empty])))
  expect_true(all(is.na(fs$.derived_bartlett[empty])))
  expect_true(all(is.na(fs$determinacy_regression[empty])))
  expect_true(all(fs$reason[empty] == "no observed item"))
  expect_true(all(fs$items_observed[empty] == 0L))
  tx <- fs$row %in% c(10, 11) & fs$factor == "textual"
  expect_true(all(is.finite(fs$.derived_regression[tx])))
  expect_true(all(is.finite(fs$.derived_regression_se[tx])))
  expect_true(all(is.na(fs$.derived_bartlett[tx])))
  expect_true(all(is.na(fs$.derived_bartlett_se[tx])))
  expect_true(all(is.na(fs$determinacy_bartlett[tx])))
  expect_true(all(grepl("no item of this factor observed", fs$reason[tx])))
  # Borrowed from correlated factors, the score is far less determinate.
  expect_true(all(fs$determinacy_regression[tx] < 0.5))
  expect_true(all(fs$items_observed[tx] == 0L))
  # The other factors of those rows are scored as usual.
  other <- fs$row %in% c(10, 11) & fs$factor != "textual"
  expect_true(all(is.finite(fs$.derived_bartlett[other])))
  expect_true(all(fs$reason[other] == ""))
})

test_that("a loading fixed or constrained to zero does not make an item", {
  skip_if_not_installed("lavaan")
  # Paragraph comprehension listed under the visual factor with its
  # loading held at zero, written 0*x or a*x with a == 0: either is the
  # model written without that row, and the visual factor has three
  # items, not four. Rows 10 and 11 lack every visual item but have
  # paragraph comprehension, so they have no visual item observed, and no
  # Bartlett score for it. A rule that read only rows fixed at zero
  # counted the constrained row as a fourth visual item, observed in rows
  # 10 and 11.
  d <- fs_mar_data()
  d[c(10, 11), fs_items$visual] <- NA
  expect_false(anyNA(d[c(10, 11), "t6_paragraph_comprehension"]))
  fit_with <- function(extra) {
    suppressWarnings(lavaan::cfa(paste(fs_model, extra), data = d,
                                 std.lv = TRUE, missing = "ml"))
  }
  plain <- factor_scores(fit_with(""), data = d)
  fixed <- factor_scores(
    fit_with("\nvisual =~ 0*t6_paragraph_comprehension"), data = d)
  constrained <- factor_scores(
    fit_with("\nvisual =~ a*t6_paragraph_comprehension\na == 0"), data = d)
  # Held at zero through a chain of equalities, the loadings are returned
  # by lavaan at its optimizer's rounding (about 1e-15), not at zero; read
  # by value, they gave rows 10 and 11 the false reason that their observed
  # items do not separate the visual factor from another.
  fit_chain <- fit_with(paste0(
    "\nvisual =~ a*t6_paragraph_comprehension",
    "\nspeed =~ b*t6_paragraph_comprehension\na == b\nb == 0"))
  pt <- lavaan::parTable(fit_chain)
  expect_true(any(pt$label %in% c("a", "b") & pt$est != 0))
  chain <- factor_scores(fit_chain, data = d)
  for (fs in list(fixed, constrained, chain)) {
    expect_identical(fs$items, plain$items)
    expect_identical(fs$items_observed, plain$items_observed)
    expect_identical(fs$reason, plain$reason)
    expect_equal(fs$.derived_regression, plain$.derived_regression, tolerance = 1e-5)
    expect_equal(fs$.derived_bartlett, plain$.derived_bartlett, tolerance = 1e-5)
  }
  vis <- constrained$factor == "visual"
  expect_true(all(constrained$items[vis] == 3L))
  gone <- vis & constrained$row %in% c(10, 11)
  expect_true(all(constrained$items_observed[gone] == 0L))
  expect_true(all(is.na(constrained$.derived_bartlett[gone])))
  expect_true(all(grepl("no item of this factor observed",
                        constrained$reason[gone])))
  # A higher order factor with direct loadings held at zero through a
  # chain has no items of its own, as with the loadings written 0*x or
  # left out: lavaan's Bartlett score of it is formed through the loadings
  # of the factors it measures, and it is the same in every spelling. Read
  # by value, the chain gave the factor items of its own, and its Bartlett
  # score was NA in every row (DMAR's check) or, from lavaan, off by 5.
  ho <- paste(fs_model, "\ng =~ visual + textual + speed")
  fit_ho <- function(extra) {
    suppressWarnings(lavaan::cfa(paste(ho, extra), data = d,
                                 std.lv = TRUE, missing = "ml"))
  }
  plain_ho <- factor_scores(fit_ho(""), data = d)
  fit_ho_chain <- fit_ho(paste0(
    "\ng =~ a*t1_visual_perception + b*t6_paragraph_comprehension",
    "\na == b\nb == 0"))
  expect_equal(lavaan::fitMeasures(fit_ho_chain, c("chisq", "df")),
               lavaan::fitMeasures(fit_ho(""), c("chisq", "df")),
               tolerance = 1e-6)
  for (fit in list(fit_ho_chain, fit_ho("\ng =~ 0*t1_visual_perception"))) {
    fs <- factor_scores(fit, data = d)
    expect_identical(fs$reason, plain_ho$reason)
    expect_identical(is.na(fs$.derived_bartlett), is.na(plain_ho$.derived_bartlett))
    expect_equal(fs$.derived_bartlett, plain_ho$.derived_bartlett, tolerance = 1e-4)
    expect_equal(fs$.derived_regression, plain_ho$.derived_regression, tolerance = 1e-4)
  }
  expect_identical(sum(!is.na(plain_ho$.derived_bartlett[plain_ho$factor == "g"])),
                   299L)
})

test_that("an entry held at zero is placed in its cell in the fit's own order", {
  skip_if_not_installed("lavaan")
  # The data's columns in the reverse of the model's order, fitted with
  # ov_order = "data", so that the fit's matrices order the tests as the
  # data do. Placed by the parameter table alone (in the order of the
  # model's syntax), the cross-loading of t1 on speed, written 0*x, was
  # taken for the cell of t13 on speed, a real loading of 24.1, which was
  # set to zero in the matrices and in the fit given to lavPredict(): the
  # regression scores moved by up to 2.7 and the Bartlett scores by 3.8.
  # Held through a chain of equalities, the cross-loading of paragraph
  # comprehension gave rows 10 and 11 (no visual item observed) the false
  # reason that their items do not separate the visual factor from
  # another, as in model order.
  d <- fs_mar_data()
  d[c(10, 11), fs_items$visual] <- NA
  d <- d[, rev(names(d))]
  fit_with <- function(extra) {
    suppressWarnings(lavaan::cfa(paste(fs_model, extra), data = d,
                                 std.lv = TRUE, missing = "ml",
                                 ov_order = "data"))
  }
  fit_p <- fit_with("")
  expect_identical(rownames(lavaan::lavInspect(fit_p, "est")$lambda),
                   names(d))
  plain <- factor_scores(fit_p, data = d)
  fit_c <- fit_with(paste0(
    "\nvisual =~ a*t6_paragraph_comprehension",
    "\nspeed =~ b*t6_paragraph_comprehension\na == b\nb == 0"))
  pt <- lavaan::parTable(fit_c)
  expect_true(any(pt$label %in% c("a", "b") & pt$est != 0))
  cells <- .fs_held_cells(fit_c)
  lam <- lavaan::lavInspect(fit_c, "est")$lambda
  expect_setequal(paste(rownames(lam)[cells$row], colnames(lam)[cells$col]),
                  c("t6_paragraph_comprehension visual",
                    "t6_paragraph_comprehension speed"))
  for (fit in list(fit_with("\nspeed =~ 0*t1_visual_perception"), fit_c)) {
    fs <- factor_scores(fit, data = d)
    expect_identical(fs$items, plain$items)
    expect_identical(fs$reason, plain$reason)
    expect_equal(fs$.derived_regression, plain$.derived_regression, tolerance = 1e-5)
    expect_equal(fs$.derived_bartlett, plain$.derived_bartlett, tolerance = 1e-5)
  }
  # A residual variance fixed at zero (the single indicator remedy's
  # z ~~ 0*z) is zero in the matrices already and is left there; placed
  # in model order it was taken for the residual variance of word meaning,
  # which was set to zero, and every determinacy was NA. The scores and
  # determinacies are lavaan's.
  dc <- holzinger_swineford[, rev(unlist(fs_items))]
  fit_z <- lavaan::cfa(paste(fs_model, "\nt6_paragraph_comprehension ~~",
                             "0*t6_paragraph_comprehension"),
                       data = dc, std.lv = TRUE, ov_order = "data")
  expect_null(.fs_held_cells(fit_z))
  fz <- factor_scores(fit_z)
  lp <- lavaan::lavPredict(fit_z)
  for (f in names(fs_items)) {
    expect_equal(fz$.derived_regression[fz$factor == f], unname(lp[, f]),
                 tolerance = 1e-8)
  }
  expect_equal(attr(fz, "determinacy")$determinacy_regression,
               unname(lavaan::lavInspect(fit_z, "fs.determinacy")),
               tolerance = 1e-6)
  # A cell that does not hold the row's estimate is left as it is.
  wrong <- cells
  wrong$row[] <- match("t13_straight_and_curved_capitals", rownames(lam))
  wrong$col[] <- match("speed", colnames(lam))
  est <- .fs_as_list(lavaan::lavInspect(fit_c, "est"), inner = "lambda")
  expect_identical(.fs_zero_held(est, wrong), est)
  expect_identical(.fs_zero_held_fit(fit_c, wrong)@Model@GLIST,
                   fit_c@Model@GLIST)
})

test_that("ordered items: a cross-loading held at zero through a chain is no loading", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # EBM and ML scores optimize case by case
  # Two factors of three Holzinger and Swineford tests cut at their
  # quartiles, with cross-loadings held at zero through a == b, b == 0:
  # the model of simple structure (the same chi square and degrees of
  # freedom). Read by value, lavaan's rounding of the held loadings (about
  # 1e-16) made an item in a middle category bound the other factor's
  # likelihood, and most rows whose items of a factor all sit at one end
  # were told the likelihood rises as the factor moves with others that
  # share its items, which no item does.
  its <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
           "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  db <- as.data.frame(lapply(holzinger_swineford[, its], function(x) {
    cut(x, stats::quantile(x, 0:4 / 4), include.lowest = TRUE,
        labels = FALSE)
  }))
  names(db) <- paste0("y", 1:6)
  m <- "f1 =~ y1 + y2 + y3\nf2 =~ y4 + y5 + y6"
  fit_with <- function(extra) {
    lavaan::sem(paste(m, extra), data = db, std.lv = TRUE, ordered = TRUE)
  }
  fit_p <- fit_with("")
  fit_c <- fit_with("\nf1 =~ a*y4\nf2 =~ b*y1\na == b\nb == 0")
  expect_equal(lavaan::fitMeasures(fit_c, c("chisq", "df")),
               lavaan::fitMeasures(fit_p, c("chisq", "df")),
               tolerance = 1e-6)
  pt <- lavaan::parTable(fit_c)
  expect_true(any(pt$label %in% c("a", "b") & pt$est != 0))
  fp <- factor_scores(fit_p, data = db)
  fc <- factor_scores(fit_c, data = db)
  ext <- grepl("every observed item of this factor is at an end", fp$reason)
  expect_identical(sum(ext), 107L)
  expect_identical(fc$reason, fp$reason)
  expect_identical(is.na(fc$.derived_bartlett), is.na(fp$.derived_bartlett))
  expect_equal(fc$.derived_bartlett, fp$.derived_bartlett, tolerance = 1e-5)
  expect_equal(fc$.derived_regression, fp$.derived_regression, tolerance = 1e-5)
})

test_that("SIEVER's factor score records are reproduced", {
  skip_if_not_installed("lavaan")
  # Values recorded from SIEVER's sieve_factor_score_record() on this
  # fixture (computed through lavaan 0.7-3); the Bartlett and regression
  # estimates, standard errors, and limits agreed with this function
  # exactly.
  d <- fs_mar_data()
  fit <- fs_fit_fiml(d)
  fs <- factor_scores(fit, data = d)
  pick <- function(col, f, r) fs[[col]][fs$factor == f & fs$row == r]
  expect_equal(pick(".derived_bartlett", "visual", 1), -1.139120747510906,
               tolerance = 1e-6)
  expect_equal(pick(".derived_bartlett_se", "visual", 1), 0.748523476486752,
               tolerance = 1e-6)
  expect_equal(pick(".derived_bartlett", "textual", 8), 0.287002678812326,
               tolerance = 1e-6)
  expect_equal(pick(".derived_bartlett_se", "textual", 8), 0.436587070979264,
               tolerance = 1e-6)
  expect_equal(pick(".derived_bartlett_ci_upper", "speed", 10), 0.8957240234923352,
               tolerance = 1e-6)
  expect_equal(pick(".derived_regression", "visual", 10), -1.1315366372845874,
               tolerance = 1e-6)
  expect_equal(pick(".derived_regression_se", "speed", 8), 0.527062699525430,
               tolerance = 1e-6)
  expect_equal(pick(".derived_regression_ci_lower", "textual", 1),
               -0.825243522359024, tolerance = 1e-6)
})

test_that("two groups: rows without a group are kept with data", {
  skip_if_not_installed("lavaan")
  d <- fs_mar_data(extra = "school")
  d$school <- as.character(d$school)
  d$school[c(3, 300, 301)] <- NA
  # lavaan warns about the missing group values and the empty rows.
  fit <- suppressWarnings(lavaan::cfa(
    fs_model, data = d, std.lv = TRUE, missing = "ml", group = "school",
    group.equal = c("loadings", "intercepts")))
  fs <- factor_scores(fit, data = d)
  expect_fs_shape(fs)
  expect_identical(nrow(fs), 3L * nrow(d))
  idx <- lavaan::lavInspect(fit, "case.idx")
  lp <- lavaan::lavPredict(fit, method = "Bartlett", assemble = FALSE)
  for (g in seq_along(idx)) {
    expect_identical(fs$.derived_bartlett[fs$factor == "speed"][idx[[g]]],
                     unname(lp[[g]][, "speed"]))
    expect_true(all(fs$group[fs$factor == "speed"][idx[[g]]] ==
                      names(idx)[g]))
  }
  nog <- fs$row %in% c(3, 300, 301)
  expect_true(all(is.na(fs$.derived_regression[nog])))
  expect_true(all(is.na(fs$group[nog])))
  expect_true(all(grepl("group variable school is missing", fs$reason[nog])))
  expect_true(all(fs$items[nog] == 3L))
  # Group-specific determinacy, as lavaan reports it.
  det <- attr(fs, "determinacy")
  lav <- lavaan::lavInspect(fit, "fs.determinacy")
  for (g in names(lav)) {
    expect_equal(det$determinacy_regression[det$group == g],
                 unname(lav[[g]]), tolerance = 1e-10)
  }
  expect_identical(attr(fs, "ngroups"), 2L)
  # Without data the table stops at the last row the fit holds.
  fs0 <- factor_scores(fit)
  expect_identical(max(fs0$row), max(unlist(idx)))
  expect_identical(fs0$.derived_regression, fs$.derived_regression[fs$row <= max(unlist(idx))])
})

test_that("data that are not the rows given to the fit are refused", {
  skip_if_not_installed("lavaan")
  d <- fs_mar_data()
  fit <- fs_fit_fiml(d)
  expect_error(factor_scores(fit, data = d[nrow(d):1, ]),
               "values of .* differ")
  expect_error(factor_scores(fit, data = d[-1, ]), "has 300 rows")
  expect_error(factor_scores(fit, data = d[, -1]),
               "does not hold the model's observed variable")
  expect_error(factor_scores(fit, data = "d"), "data frame")
  dg <- fs_mar_data(extra = "school")
  fitg <- suppressWarnings(lavaan::cfa(fs_model, data = dg, std.lv = TRUE,
                                       missing = "ml", group = "school"))
  expect_error(factor_scores(fitg, data = dg[, unlist(fs_items)]),
               "group variable school")
  dg2 <- dg
  dg2$school <- as.character(dg2$school)
  dg2$school[1] <- "Elsewhere"
  expect_error(factor_scores(fitg, data = dg2), "not fitted to")
})

test_that("a group the fit leaves out (group.label) is kept with a reason", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[, c(unlist(fs_items), "school")]
  fit <- lavaan::cfa(fs_model, data = d, std.lv = TRUE, group = "school",
                     group.label = "Pasteur")
  fs <- factor_scores(fit, data = d)
  gw <- d$school[fs$row] == "Grant-White"
  expect_true(all(is.na(fs$.derived_regression[gw])))
  expect_true(all(grepl("group.label", fs$reason[gw])))
  expect_true(all(fs$group[gw] == "Grant-White"))
  expect_true(all(is.finite(fs$.derived_bartlett[!gw])))
})

test_that("a listwise fit scores its complete rows and names the others", {
  skip_if_not_installed("lavaan")
  d <- fs_mar_data()
  fit <- lavaan::cfa(fs_model, data = d, std.lv = TRUE)
  fs <- factor_scores(fit, data = d)
  cc <- stats::complete.cases(d)
  expect_identical(sum(is.finite(fs$.derived_bartlett[fs$factor == "visual"])), sum(cc))
  gone <- !cc[fs$row]
  expect_true(all(is.na(fs$.derived_regression[gone])))
  expect_true(all(grepl("deleted listwise", fs$reason[gone])))
  # The items answered by a deleted row are counted from the data.
  expect_identical(fs$items_observed[fs$factor == "speed"],
                   as.integer(rowSums(!is.na(d[, fs_items$speed]))))
  expect_true(all(fs$reason[!gone] == ""))
  fs0 <- factor_scores(fit)
  expect_identical(nrow(fs0), nrow(fs))
  expect_identical(fs0$.derived_bartlett, fs$.derived_bartlett)
  expect_true(all(grepl("deleted listwise", fs0$reason[gone])))
  expect_identical(attr(fs, "missing"), "listwise")
})

test_that("auxiliary variables are set aside; the indicators are scored", {
  skip_if_not_installed("lavaan")
  # A factor with five items, missingness that depends on an auxiliary
  # variable correlated with the factor, and the auxiliary carried as a
  # saturated correlate.
  set.seed(113)
  n <- 300
  lam <- c(.8, .7, .6, .6, .5)
  eta <- stats::rnorm(n)
  x <- sapply(lam, function(l) l * eta + stats::rnorm(n, sd = sqrt(1 - l^2)))
  colnames(x) <- paste0("y", 1:5)
  z <- stats::rnorm(n, mean = 0.6 * eta, sd = 0.8)
  for (j in 2:5) x[stats::runif(n) < stats::plogis(-1.4 + 1.2 * z), j] <- NA
  d <- data.frame(x, z = z)
  d[3, 1:5] <- NA
  fit <- suppressWarnings(lavaan::cfa(
    "f =~ y1 + y2 + y3 + y4 + y5
     z ~~ y1 + y2 + y3 + y4 + y5", data = d, std.lv = TRUE,
    missing = "ml"))
  expect_message(fs <- factor_scores(fit, data = d), "saturated correlates")
  expect_identical(attr(fs, "scored_from"), "indicators")
  expect_identical(attr(fs, "auxiliary"), "z")
  expect_fs_shape(fs)
  # The oracle: lavPredict() on the same model without the auxiliary,
  # every parameter fixed at the auxiliary fit's estimate.
  pt <- lavaan::parTable(fit)
  pt <- pt[!(pt$lhs == "z" | pt$rhs == "z"), ]
  pt$free <- 0L
  pt$ustart <- pt$est
  pt$start <- pt$est
  pt$se <- NULL
  fixed <- suppressWarnings(lavaan::lavaan(pt, data = d[, 1:5],
                                           missing = "ml"))
  idx <- lavaan::lavInspect(fixed, "case.idx")
  for (m in c("regression", "Bartlett")) {
    lp <- lavaan::lavPredict(fixed, method = m, se = "standard")
    col <- paste0(".derived_", tolower(m))
    expect_equal(fs[[col]][idx], unname(lp[, "f"]), tolerance = 1e-10)
    expect_equal(fs[[paste0(col, "_se")]][idx],
                 unname(attr(lp, "se")[[1]][, "f"]), tolerance = 1e-10)
  }
  # A row whose indicators are all missing is not scored from the
  # auxiliary alone.
  expect_true(is.na(fs$.derived_regression[3]))
  expect_identical(fs$reason[3], "no observed item")
  # The indicator-only scores track the factor far better than scores
  # that condition on the auxiliary.
  r_fit <- stats::cor(lavaan::lavPredict(fit)[, "f"], eta,
                      use = "complete.obs")^2
  r_ind <- stats::cor(fs$.derived_regression, eta, use = "complete.obs")^2
  expect_gt(r_ind, r_fit + 0.2)
  # The auxiliary written with a loading held at zero, 0*z or a*z with
  # a == 0, is the same model (the same log likelihood and degrees of
  # freedom), and it is set aside the same way. Read as an indicator, it
  # was conditioned on, and the scores differed from these.
  for (extra in c("f =~ 0*z", "f =~ a*z\na == 0")) {
    fit_0 <- suppressWarnings(lavaan::cfa(
      paste("f =~ y1 + y2 + y3 + y4 + y5
       z ~~ y1 + y2 + y3 + y4 + y5\n", extra), data = d, std.lv = TRUE,
      missing = "ml"))
    expect_equal(lavaan::fitMeasures(fit_0, c("logl", "df")),
                 lavaan::fitMeasures(fit, c("logl", "df")), tolerance = 1e-8)
    expect_message(fs_0 <- factor_scores(fit_0, data = d),
                   "saturated correlates")
    expect_identical(attr(fs_0, "scored_from"), "indicators")
    expect_identical(attr(fs_0, "auxiliary"), "z")
    expect_identical(fs_0$reason, fs$reason)
    expect_equal(fs_0$.derived_regression, fs$.derived_regression, tolerance = 1e-5)
    expect_equal(fs_0$.derived_bartlett, fs$.derived_bartlett, tolerance = 1e-5)
  }
})

test_that("a two-stage fit is scored from each case's observed items, as under FIML", {
  skip_if_not_installed("lavaan")
  # lavPredict() scores a two-stage fit as complete data, so every case
  # with a missing item, which the fit used, had no score ("score
  # unavailable from lavaan") beside the standard error of the complete
  # pattern. The scores are now computed here with the two-stage
  # estimates, pattern by pattern.
  d <- fs_mar_data()
  # lavaan warns that the two rows with every item missing are ignored.
  fit <- suppressWarnings(lavaan::cfa(fs_model, data = d, std.lv = TRUE,
                                      missing = "two.stage"))
  fs <- factor_scores(fit, data = d)
  expect_fs_shape(fs)
  expect_identical(attr(fs, "missing"), "two.stage")
  expect_identical(attr(fs, "scored_from"), "fit")
  # The oracle: lavaan's own FIML scoring, pattern by pattern, of the same
  # model with every parameter fixed at its two-stage estimate.
  pt <- lavaan::parTable(fit)
  pt$free <- 0L
  pt$ustart <- pt$est
  pt$start <- pt$est
  pt$se <- NULL
  fixed <- suppressWarnings(lavaan::lavaan(pt, data = d, missing = "ml"))
  idx <- lavaan::lavInspect(fixed, "case.idx")
  for (m in c("regression", "Bartlett")) {
    lp <- lavaan::lavPredict(fixed, method = m, se = "standard")
    col <- paste0(".derived_", tolower(m))
    for (f in names(fs_items)) {
      rows <- fs$factor == f
      expect_equal(fs[[col]][rows][idx], unname(lp[, f]), tolerance = 1e-10)
      expect_equal(fs[[paste0(col, "_se")]][rows][idx],
                   unname(attr(lp, "se")[[1]][, f]), tolerance = 1e-10)
    }
  }
  # Every row with an observed item is scored, and nothing is attributed
  # to lavaan; the complete rows are lavaan's own scores of the two-stage
  # fit.
  empty <- fs$row %in% c(5, 77)
  expect_true(all(is.finite(fs$.derived_regression[!empty])))
  expect_true(all(fs$reason[empty] == "no observed item"))
  expect_false(any(grepl("lavaan", fs$reason)))
  idx_ts <- lavaan::lavInspect(fit, "case.idx")
  cc <- stats::complete.cases(d[idx_ts, ])
  lp_ts <- lavaan::lavPredict(fit)
  for (f in names(fs_items)) {
    expect_equal(fs$.derived_regression[fs$factor == f][idx_ts][cc],
                 unname(lp_ts[cc, f]), tolerance = 1e-10)
  }
  # The robust variant has the same estimates, so the same scores; so does
  # fa_confirmatory(missing = "two_stage"), up to the optimizer's tolerance.
  fit_r <- suppressWarnings(lavaan::cfa(fs_model, data = d, std.lv = TRUE,
                                        missing = "robust.two.stage"))
  expect_equal(factor_scores(fit_r, data = d)$.derived_regression, fs$.derived_regression,
               tolerance = 1e-10)
  fit_k <- suppressWarnings(fa_confirmatory(d, fs_items, missing = "two_stage",
                                  output = "fit"))
  expect_equal(factor_scores(fit_k, data = d)$.derived_bartlett, fs$.derived_bartlett,
               tolerance = 1e-4)
  # A higher order factor has no Bartlett score when the scores are
  # computed here, and its reason does not credit lavaan.
  fit_g <- suppressWarnings(lavaan::cfa(
    paste(fs_model, "\ng =~ visual + textual + speed"), data = d,
    std.lv = TRUE, missing = "two.stage"))
  fs_g <- factor_scores(fit_g, data = d)
  g <- fs_g$factor == "g" & !(fs_g$row %in% c(5, 77))
  expect_true(all(is.na(fs_g$.derived_bartlett[g])))
  expect_true(all(is.finite(fs_g$.derived_regression[g])))
  expect_identical(unique(fs_g$reason[g]), paste0(
    "no Bartlett score for a factor without items of its own when the ",
    "scores are computed from the two-stage estimates rather than by ",
    "lavaan"))
  # Conditional on exogenous covariates, the computation here does not
  # apply, and lavaan's scores leave out the incomplete cases: refused.
  dx <- fs_mar_data("age_years")
  fit_x <- suppressWarnings(lavaan::sem(
    paste(fs_model, "\nvisual ~ age_years"), data = dx, std.lv = TRUE,
    missing = "two.stage", conditional.x = TRUE))
  expect_error(factor_scores(fit_x, data = dx),
               "A two-stage fit conditional on exogenous covariates",
               fixed = TRUE)
})

test_that("a fit by pairwise deletion is scored from each case's observed items", {
  skip_if_not_installed("lavaan")
  # lavPredict() scores a fit by pairwise deletion of continuous items as
  # complete data, so every case with a missing item had no score ("score
  # unavailable from lavaan") beside the standard error of the complete
  # pattern, 459 of the 903 cells here.
  d <- fs_mar_data()
  # lavaan warns that the two rows with every item missing are ignored.
  fit <- suppressWarnings(lavaan::cfa(fs_model, data = d, std.lv = TRUE,
                                      missing = "pairwise"))
  fs <- factor_scores(fit, data = d)
  expect_fs_shape(fs)
  expect_identical(attr(fs, "missing"), "pairwise")
  expect_identical(attr(fs, "scored_from"), "fit")
  empty <- fs$row %in% c(5, 77)
  expect_true(all(is.finite(fs$.derived_regression[!empty])))
  expect_true(all(is.finite(fs$.derived_bartlett[!empty])))
  expect_true(all(fs$reason[empty] == "no observed item"))
  expect_false(any(grepl("lavaan", fs$reason)))
  expect_true(all(is.na(fs$.derived_regression_se[is.na(fs$.derived_regression)])))
  # The complete rows are lavaan's own scores of the fit, which has no mean
  # structure: they also check the means the scores are centered at.
  idx <- lavaan::lavInspect(fit, "case.idx")
  cc <- stats::complete.cases(d[idx, ])
  for (m in c("regression", "Bartlett")) {
    lp <- lavaan::lavPredict(fit, method = m, se = "standard")
    col <- paste0(".derived_", tolower(m))
    for (f in names(fs_items)) {
      expect_equal(fs[[col]][fs$factor == f][idx][cc], unname(lp[cc, f]),
                   tolerance = 1e-10)
      # lavaan gives one row of standard errors, that of complete data.
      expect_equal(fs[[paste0(col, "_se")]][fs$factor == f][idx][cc],
                   rep(unname(attr(lp, "se")[[1]][1L, f]), sum(cc)),
                   tolerance = 1e-10)
    }
  }
  # Every incomplete case: the closed forms with the estimates from
  # pairwise deletion, centered at the means of the observed values.
  est <- lavaan::lavInspect(fit, "est")
  lam <- est$lambda
  th <- est$theta
  phi <- est$psi
  x <- as.matrix(d[, rownames(lam)])
  mu <- colMeans(x, na.rm = TRUE)
  sig <- lam %*% phi %*% t(lam) + th
  rows <- which(rowSums(!is.na(x)) > 0L & !stats::complete.cases(x))
  want <- lapply(rows, function(i) {
    o <- !is.na(x[i, ])
    yc <- x[i, o] - mu[o]
    lo <- lam[o, , drop = FALSE]
    info <- t(lo) %*% solve(th[o, o]) %*% lo
    list(reg = drop(phi %*% t(lo) %*% solve(sig[o, o], yc)),
         bar = drop(solve(info, t(lo) %*% solve(th[o, o], yc))),
         bar_se = sqrt(diag(solve(info))))
  })
  got <- function(col) sapply(colnames(lam), function(f) {
    fs[[col]][fs$factor == f]
  })[rows, , drop = FALSE]
  stack <- function(el) unname(do.call(rbind, lapply(want, `[[`, el)))
  expect_equal(unname(got(".derived_regression")), stack("reg"), tolerance = 1e-10)
  expect_equal(unname(got(".derived_bartlett")), stack("bar"), tolerance = 1e-10)
  expect_equal(unname(got(".derived_bartlett_se")), stack("bar_se"),
               tolerance = 1e-10)
  # The available cases give lavaan the same estimates, and so the same
  # scores.
  fit_a <- suppressWarnings(lavaan::cfa(fs_model, data = d, std.lv = TRUE,
                                        missing = "available.cases"))
  expect_equal(factor_scores(fit_a, data = d)$.derived_bartlett, fs$.derived_bartlett,
               tolerance = 1e-10)
  # Conditional on exogenous covariates, the computation here does not
  # apply: refused, as for a two-stage fit.
  dx <- fs_mar_data("age_years")
  fit_x <- suppressWarnings(lavaan::sem(
    paste(fs_model, "\nvisual ~ age_years"), data = dx, std.lv = TRUE,
    missing = "pairwise", conditional.x = TRUE))
  expect_error(factor_scores(fit_x, data = dx),
               "A fit by pairwise deletion conditional on exogenous", fixed = TRUE)
  # A higher order factor's reason names the estimates the scores use.
  fit_g <- suppressWarnings(lavaan::cfa(
    paste(fs_model, "\ng =~ visual + textual + speed"), data = d,
    std.lv = TRUE, missing = "pairwise"))
  fs_g <- factor_scores(fit_g, data = d)
  g <- fs_g$factor == "g" & !(fs_g$row %in% c(5, 77))
  expect_true(all(is.na(fs_g$.derived_bartlett[g])))
  expect_identical(unique(fs_g$reason[g]), paste0(
    "no Bartlett score for a factor without items of its own when the ",
    "scores are computed from the estimates from pairwise deletion rather ",
    "than by lavaan"))
})

# Three spatial tests and age, with the third test missing for 40 pupils:
# the fixture of the models with a covariate below.
fs_mimic_data <- function() {
  d <- holzinger_swineford[, c("t1_visual_perception", "t2_cubes",
                               "t4_lozenges", "age_years")]
  names(d) <- c("y1", "y2", "y3", "x")
  set.seed(4)
  d$y3[sample(nrow(d), 40)] <- NA
  d
}

# The oracle for scores computed here: lavPredict() under FIML of the
# same model with every parameter fixed at the fit's estimate (the fit
# must have a mean structure), aligned to the rows of the data.
fs_fixed_oracle <- function(fit, d, method) {
  pt <- lavaan::parTable(fit)
  pt$free <- 0L
  pt$ustart <- pt$est
  pt$start <- pt$est
  pt$se <- NULL
  fixed <- lavaan::lavaan(pt, data = d, missing = "ml")
  lp <- lavaan::lavPredict(fixed, method = method, se = "standard")
  idx <- lavaan::lavInspect(fixed, "case.idx")
  score <- se <- matrix(NA_real_, nrow(d), ncol(lp),
                        dimnames = list(NULL, colnames(lp)))
  score[idx, ] <- lp
  se[idx, ] <- attr(lp, "se")[[1]]
  list(score = score, se = se)
}

test_that("without a mean structure a covariate's factor is centered as lavaan centers it", {
  skip_if_not_installed("lavaan")
  # A fit by pairwise deletion without a mean structure (the default for one
  # group) took every latent intercept as zero, the covariate's included,
  # while the covariate itself was centered at its mean. Every regression
  # and Bartlett score of a factor the covariate predicts, complete rows
  # included, was shifted by -beta times the covariate's mean (-0.26 here),
  # against lavPredict() of the same fit, the fit with a mean structure, and
  # FIML.
  d <- fs_mimic_data()
  model <- "f =~ y1 + y2 + y3\nf ~ x"
  cc <- stats::complete.cases(d)
  for (miss in c("pairwise", "available.cases")) {
    # With the available cases, lavaan warns that its Browne residual test
    # cannot be computed, which has no bearing on the estimates.
    fit <- suppressWarnings(lavaan::sem(model, data = d, std.lv = TRUE,
                                        missing = miss))
    expect_false(lavaan::lavInspect(fit, "meanstructure"))
    fit_ms <- suppressWarnings(lavaan::sem(model, data = d, std.lv = TRUE,
                                           missing = miss,
                                           meanstructure = TRUE))
    fs <- factor_scores(fit, data = d)
    fs_ms <- factor_scores(fit_ms, data = d)
    # The two fits' estimates agree to about 1e-6, and so do the scores.
    expect_equal(fs$.derived_regression, fs_ms$.derived_regression, tolerance = 1e-5)
    expect_equal(fs$.derived_bartlett, fs_ms$.derived_bartlett, tolerance = 1e-5)
    # The complete rows are lavaan's own scores of the fit.
    for (m in c("regression", "Bartlett")) {
      lp <- lavaan::lavPredict(fit, method = m)
      expect_equal(fs[[paste0(".derived_", tolower(m))]][cc], unname(lp[cc, "f"]),
                   tolerance = 1e-10)
    }
    # Every row is FIML's with the estimates held fixed.
    for (m in c("regression", "Bartlett")) {
      o <- fs_fixed_oracle(fit_ms, d, m)
      expect_equal(fs_ms[[paste0(".derived_", tolower(m))]], unname(o$score[, "f"]),
                   tolerance = 1e-10)
    }
    # The mean the scores are centered at is lavaan's.
    expect_equal(mean(fs$.derived_regression[cc] - lavaan::lavPredict(fit)[cc, "f"]),
                 0, tolerance = 1e-10)
  }
  # The scores from the indicators with an auxiliary variable set aside had
  # the same shift in a fit without a mean structure (-0.67 here).
  da <- holzinger_swineford[, c("t1_visual_perception", "t2_cubes",
                                "t4_lozenges", "age_years",
                                "t5_general_information")]
  names(da) <- c("y1", "y2", "y3", "x", "a")
  model_a <- paste(model, "\na ~~ y1 + y2 + y3")
  fit_a <- lavaan::sem(model_a, data = da, std.lv = TRUE)
  fit_am <- lavaan::sem(model_a, data = da, std.lv = TRUE,
                        meanstructure = TRUE)
  fa <- suppressMessages(factor_scores(fit_a))
  fam <- suppressMessages(factor_scores(fit_am))
  expect_equal(fa$.derived_regression, fam$.derived_regression, tolerance = 1e-5)
  expect_equal(fa$.derived_bartlett, fam$.derived_bartlett, tolerance = 1e-5)
  # An outcome the factor predicts, with and without a covariate.
  dz <- d
  dz$z <- holzinger_swineford$t5_general_information
  for (mz in c("f =~ y1 + y2 + y3\nz ~ f", paste(model, "\nz ~ f"))) {
    fz <- lavaan::sem(mz, data = dz, std.lv = TRUE, missing = "pairwise")
    fzm <- lavaan::sem(mz, data = dz, std.lv = TRUE, missing = "pairwise",
                       meanstructure = TRUE)
    expect_equal(factor_scores(fz, data = dz)$.derived_regression,
                 factor_scores(fzm, data = dz)$.derived_regression, tolerance = 1e-5)
  }
})

test_that("without a mean structure an outcome's mean is counted once", {
  skip_if_not_installed("lavaan")
  # Without a mean structure or an exogenous covariate, lavPredict() sets
  # the intercept of an outcome the factor predicts (z ~ f) to its sample
  # mean, which the latent variable lavaan holds it with already carries,
  # so z's mean is counted twice and every score of f is shifted: by -0.38
  # (regression) and -0.57 (Bartlett) here, against the same fit with a
  # mean structure. The scores of such a fit are computed here.
  d <- holzinger_swineford[, c("t1_visual_perception", "t2_cubes",
                               "t4_lozenges", "t5_general_information")]
  names(d) <- c("y1", "y2", "y3", "z")
  model <- "f =~ y1 + y2 + y3\nz ~ f"
  fit <- lavaan::sem(model, data = d, std.lv = TRUE)
  fit_ms <- lavaan::sem(model, data = d, std.lv = TRUE, meanstructure = TRUE)
  shift <- lavaan::lavPredict(fit)[, "f"] - lavaan::lavPredict(fit_ms)[, "f"]
  expect_lt(max(shift), -0.3)
  fs <- factor_scores(fit)
  fs_ms <- factor_scores(fit_ms)
  for (col in c(".derived_regression", ".derived_regression_se",
                ".derived_bartlett", ".derived_bartlett_se",
                "determinacy_regression", "determinacy_bartlett")) {
    expect_equal(fs[[col]], fs_ms[[col]], tolerance = 1e-5)
  }
  # The scores of the fit with a mean structure are lavaan's.
  for (m in c("regression", "Bartlett")) {
    lp <- lavaan::lavPredict(fit_ms, method = m, se = "standard")
    expect_equal(fs_ms[[paste0(".derived_", tolower(m))]], unname(lp[, "f"]), tolerance = 1e-10)
    # lavaan gives one row of standard errors, that of complete data.
    expect_equal(fs[[paste0(paste0(".derived_", tolower(m)), "_se")]],
                 rep(unname(attr(lp, "se")[[1]][1L, "f"]), nrow(d)),
                 tolerance = 1e-5)
  }
  expect_identical(attr(fs, "scored_from"), "fit")
  expect_true(all(fs$reason == ""))
  # The double count shifts every score whose weights involve z: the
  # regression score of a second factor correlated with f, which does not
  # predict z, as well (by -0.07 here), though not its Bartlett score.
  d2 <- d
  d2[c("y4", "y5", "y6")] <- holzinger_swineford[
    , c("t10_addition", "t12_counting_groups_of_dots",
        "t13_straight_and_curved_capitals")]
  model2 <- "f =~ y1 + y2 + y3\ng =~ y4 + y5 + y6\nz ~ f"
  fit2 <- lavaan::sem(model2, data = d2, std.lv = TRUE)
  fit2_ms <- lavaan::sem(model2, data = d2, std.lv = TRUE,
                         meanstructure = TRUE)
  for (m in c("regression", "Bartlett")) {
    shift <- colMeans(lavaan::lavPredict(fit2, method = m) -
                        lavaan::lavPredict(fit2_ms, method = m))
    expect_equal(unname(shift), if (m == "regression") {
      c(-0.390, -0.0697)
    } else {
      c(-0.613, 0)
    }, tolerance = 1e-3)
  }
  fs2 <- factor_scores(fit2)
  fs2_ms <- factor_scores(fit2_ms)
  for (col in c(".derived_regression", ".derived_regression_se",
                ".derived_bartlett", ".derived_bartlett_se")) {
    expect_equal(fs2[[col]], fs2_ms[[col]], tolerance = 1e-5)
  }
  # With an exogenous covariate in the model lavaan's means are right, and
  # its scores are kept.
  dx <- d
  dx$x <- holzinger_swineford$age_years
  fit_x <- lavaan::sem("f =~ y1 + y2 + y3\nz ~ f\nf ~ x", data = dx,
                       std.lv = TRUE)
  expect_identical(factor_scores(fit_x)$.derived_regression,
                   unname(lavaan::lavPredict(fit_x)[, "f"]))
})

test_that("the model implied covariance matrix is formed from the raw matrices", {
  skip_if_not_installed("lavaan")
  # The scores' Sigma is lavaan's model implied covariance matrix. The
  # matrices .fs_dummies() builds for lavPredict()'s loadings and lavaan's
  # residual covariance matrix (lambda_eff, theta_eff) are not a
  # decomposition of it: a variable that keeps its unit row has its block
  # of Psi in both (f ~~ z, f ~ x), and a chain z1 ~ f, z2 ~ z1 is not
  # reproduced. Pinned so that the eff form is never swapped in.
  hs <- holzinger_swineford
  d <- data.frame(y1 = hs$t1_visual_perception, y2 = hs$t2_cubes,
                  y3 = hs$t4_lozenges, z = hs$t3_paper_form_board,
                  z2 = hs$t5_general_information, x = hs$age_years)
  for (m in c("f ~~ z", "f ~ x", "z ~ f\nz2 ~ z")) {
    fit <- lavaan::sem(paste("f =~ y1 + y2 + y3\n", m), data = d,
                       std.lv = TRUE)
    e <- .fs_dummies(lavaan::lavInspect(fit, "est"),
                     lavaan::lavNames(fit, "lv"),
                     lavaan::lavNames(fit, "ov.x"))
    mom <- .fs_moments(e, as.matrix(d))
    v <- rownames(mom$sigma)
    expect_equal(mom$sigma, lavaan::lavInspect(fit, "cov.ov")[v, v],
                 tolerance = 1e-12, ignore_attr = TRUE)
    eff <- e$lambda_eff %*% mom$veta %*% t(e$lambda_eff) + e$theta_eff
    expect_gt(max(abs(eff - mom$sigma)), 0.5)
  }
})

test_that("an item with a direct effect of a covariate enters its factor's Bartlett score", {
  skip_if_not_installed("lavaan")
  # lavaan holds an item regressed on a covariate (y1 ~ x) as a dummy latent
  # variable, with its loading on f in B rather than in Lambda. The closed
  # forms took the loadings from Lambda, so y1 was a factor of its own and
  # no item of f: the Bartlett scores of a two-stage or pairwise fit came
  # from y2 and y3 alone (off by up to 2.57 from FIML's), the pattern table
  # and the determinacy were wrong for every fit, and with an auxiliary
  # variable set aside no row had a Bartlett score.
  d <- fs_mimic_data()
  model <- "f =~ y1 + y2 + y3\nf ~ x\ny1 ~ x"
  for (miss in c("two.stage", "pairwise")) {
    fit <- lavaan::sem(model, data = d, std.lv = TRUE, missing = miss,
                       meanstructure = TRUE)
    fs <- factor_scores(fit, data = d)
    for (m in c("regression", "Bartlett")) {
      o <- fs_fixed_oracle(fit, d, m)
      col <- paste0(".derived_", tolower(m))
      expect_equal(fs[[col]], unname(o$score[, "f"]), tolerance = 1e-10)
      expect_equal(fs[[paste0(col, "_se")]], unname(o$se[, "f"]),
                   tolerance = 1e-10)
    }
  }
  # Under FIML the rows are lavaan's, and the pattern table and the
  # determinacies now agree with them and with lavaan.
  d2 <- d
  d2[1:5, c("y2", "y3")] <- NA
  fit_ml <- lavaan::sem(model, data = d2, std.lv = TRUE, missing = "ml")
  fm <- factor_scores(fit_ml, data = d2)
  lb <- lavaan::lavPredict(fit_ml, method = "Bartlett", se = "standard")
  expect_identical(fm$.derived_bartlett, unname(lb[, "f"]))
  pt <- attr(fm, "patterns")
  expect_equal(pt$bartlett_se[fm$pattern],
               unname(attr(lb, "se")[[1]][, "f"]), tolerance = 1e-10)
  # The pattern with y1 alone of f's items has a Bartlett determinacy.
  expect_identical(pt$items_observed[3], 1L)
  expect_true(all(is.finite(fm$determinacy_bartlett[1:5])))
  expect_true(all(fm$reason[1:5] == ""))
  # So does the same pattern of the two-stage fit, and a Bartlett score.
  fit_2s <- lavaan::sem(model, data = d2, std.lv = TRUE,
                        missing = "two.stage")
  f2s <- factor_scores(fit_2s, data = d2)
  expect_true(all(is.finite(f2s$.derived_bartlett[1:5])))
  expect_true(all(f2s$reason[1:5] == ""))
  # With complete data the determinacies are lavaan's.
  dc <- holzinger_swineford[, c("t1_visual_perception", "t2_cubes",
                                "t4_lozenges", "age_years")]
  names(dc) <- c("y1", "y2", "y3", "x")
  fit_c <- lavaan::sem(model, data = dc, std.lv = TRUE)
  det <- attr(factor_scores(fit_c), "determinacy")
  expect_equal(det$determinacy_regression,
               unname(lavaan::lavInspect(fit_c, "fs.determinacy")["f"]),
               tolerance = 1e-10)
  expect_equal(det$determinacy_bartlett,
               unname(lavaan::lavInspect(fit_c,
                                         "fs.determinacy.bartlett")["f"]),
               tolerance = 1e-10)
  # With an auxiliary variable set aside lavaan holds every indicator as a
  # dummy latent variable; the scores from the indicators are those of the
  # model without the auxiliary, every parameter fixed at the fit's
  # estimate.
  da <- d
  da$a <- holzinger_swineford$t5_general_information
  fit_a <- lavaan::sem(paste(model, "\na ~~ y1 + y2 + y3"), data = da,
                       std.lv = TRUE, missing = "ml")
  fa <- suppressMessages(factor_scores(fit_a, data = da))
  expect_true(all(is.finite(fa$.derived_bartlett)))
  expect_true(all(fa$reason == ""))
  pa <- lavaan::parTable(fit_a)
  pa <- pa[!(pa$lhs == "a" | pa$rhs == "a"), ]
  pa$free <- 0L
  pa$ustart <- pa$est
  pa$start <- pa$est
  pa$se <- NULL
  fixed <- lavaan::lavaan(pa, data = d, missing = "ml")
  for (m in c("regression", "Bartlett")) {
    lp <- lavaan::lavPredict(fixed, method = m, se = "standard")
    col <- paste0(".derived_", tolower(m))
    expect_equal(fa[[col]], unname(lp[, "f"]), tolerance = 1e-10)
    expect_equal(fa[[paste0(col, "_se")]], unname(attr(lp, "se")[[1]][, "f"]),
                 tolerance = 1e-10)
  }
})

test_that("an observed variable regressed on others that predicts a factor is refused", {
  skip_if_not_installed("lavaan")
  # lavPredict() scores an observed variable regressed on others as a
  # measure of the latent variables it depends on. When it predicts a
  # factor (an observed mediator) that is not its measurement equation, and
  # lavaan's regression score is not the factor's conditional mean.
  d <- holzinger_swineford[, c("t1_visual_perception", "t2_cubes",
                               "t4_lozenges", "age_years",
                               "t5_general_information")]
  names(d) <- c("y1", "y2", "y3", "x", "a")
  fit <- lavaan::sem("f =~ y1 + y2 + y3\nf ~ a\na ~ x", data = d,
                     std.lv = TRUE, meanstructure = TRUE)
  est <- lavaan::lavInspect(fit, "est")
  lam <- est$lambda
  ib <- solve(diag(ncol(lam)) - est$beta)
  v <- ib %*% est$psi %*% t(ib)
  sig <- lam %*% v %*% t(lam) + est$theta
  y <- as.matrix(lavaan::lavInspect(fit, "data"))[, rownames(lam)]
  eeta <- drop(ib %*% est$alpha)
  ey <- drop(est$nu) + drop(lam %*% eeta)
  cond_mean <- drop(sweep(y, 2L, ey) %*% solve(sig, lam %*% v[, "f"])) +
    eeta[1L]
  lp <- lavaan::lavPredict(fit)[, "f"]
  expect_gt(max(abs(lp - cond_mean)), 1)
  expect_error(factor_scores(fit), paste0(
    "The disturbance of a, an observed variable regressed on other ",
    "variables, covaries with f"), fixed = TRUE)
  # Carried as a single indicator latent variable, the mediator is scored,
  # and the scores are the conditional means (up to the factor's origin,
  # which the two parameterizations place differently).
  fit_r <- lavaan::sem("f =~ NA*y1 + y2 + y3\nf ~~ 1*f\nA =~ 1*a
                        a ~~ 0*a\nf ~ A\nA ~ x", data = d,
                       meanstructure = TRUE)
  fr <- factor_scores(fit_r)
  dif <- fr$.derived_regression[fr$factor == "f"] - cond_mean
  expect_lt(diff(range(dif)), 1e-5)
})

test_that("the single indicator remedy a refusal names is the refused model under std.lv = TRUE", {
  skip_if_not_installed("lavaan")
  # A verbal factor predicts the counting test (z ~ f), whose disturbance
  # covaries with age (z ~~ a). The remedy, refitted exactly as the
  # refusal words it and with the refused fit's std.lv = TRUE, must be the
  # same model: lavaan fixes a latent variable's variance at 1 under
  # std.lv = TRUE even when a loading is fixed, so the remedy written
  # z_lv =~ 1*z fixed both and gave a different model (log likelihood
  # -4502.0 against -4262.4) whose factor scores missed the conditional
  # means by up to 5.6.
  hs <- holzinger_swineford
  d <- data.frame(y1 = hs$t6_paragraph_comprehension, y2 = hs$t7_sentence,
                  y3 = hs$t9_word_meaning,
                  z = hs$t12_counting_groups_of_dots, a = hs$age_years)
  fit <- lavaan::sem("f =~ y1 + y2 + y3\nz ~ f\nz ~~ a", data = d,
                     std.lv = TRUE, meanstructure = TRUE)
  scores <- factor_scores(fit, data = d)
  expect_true(all(is.finite(scores$.derived_regression)))
  expect_true(all(is.na(scores$.derived_bartlett)))
  expect_true(all(grepl("Bartlett scores withheld", scores$reason)))
  # An explicit single-indicator representation remains an equivalent
  # fitting alternative, with exactly one scale constraint.
  remedy <- c("", "z_lv =~ z", "z ~~ 0*z")
  fit_r <- lavaan::sem(paste("f =~ y1 + y2 + y3", remedy[2L], remedy[3L],
                             "z_lv ~ f", "z_lv ~~ a", sep = "\n"),
                       data = d, std.lv = TRUE, meanstructure = TRUE)
  expect_equal(as.numeric(lavaan::logLik(fit_r)),
               as.numeric(lavaan::logLik(fit)), tolerance = 1e-8)
  expect_identical(lavaan::fitMeasures(fit_r, "df"),
                   lavaan::fitMeasures(fit, "df"))
  # The factor's scores are its conditional means under the refused model.
  ca <- lavaan::lavInspect(fit, "cov.all")
  ov <- lavaan::lavNames(fit, "ov")
  mu <- lavaan::lavInspect(fit, "mean.ov")[ov]
  cond_mean <- drop(sweep(as.matrix(d[, ov]), 2L, mu) %*%
                      solve(ca[ov, ov], ca[ov, "f"])) +
    lavaan::lavInspect(fit, "mean.lv")[["f"]]
  fr <- factor_scores(fit_r, data = d)
  expect_equal(fr$.derived_regression[fr$factor == "f"], unname(cond_mean),
               tolerance = 1e-5)
  # Two schools with equal loadings and intercepts: lavaan frees the new
  # latent variable's mean in the second school, so the remedy as first
  # worded had one more parameter (log likelihood -4246.42 on 15 degrees
  # of freedom, against -4246.51 on 16). The refusal now names the
  # intercept fixed at zero, and the remedy so written is the refused
  # model.
  d$school <- hs$school
  ge <- c("loadings", "intercepts")
  fit_g <- lavaan::sem("f =~ y1 + y2 + y3\nz ~ f\nz ~~ a", data = d,
                       std.lv = TRUE, group = "school", group.equal = ge)
  scores_g <- factor_scores(fit_g, data = d)
  expect_true(all(is.finite(scores_g$.derived_regression)))
  expect_true(all(is.na(scores_g$.derived_bartlett)))
  mean0 <- c("", "z_lv ~ 0*1")
  refit <- function(extra, group.equal, ...) {
    lavaan::sem(paste("f =~ y1 + y2 + y3", remedy[2L], remedy[3L],
                      "z_lv ~ f", extra, sep = "\n"),
                data = d, std.lv = TRUE, group = "school",
                group.equal = group.equal, ...)
  }
  fit_gr <- refit(paste("z_lv ~~ a", mean0[2L], sep = "\n"), ge)
  expect_equal(as.numeric(lavaan::logLik(fit_gr)),
               as.numeric(lavaan::logLik(fit_g)), tolerance = 1e-8)
  expect_identical(lavaan::fitMeasures(fit_gr, "df"),
                   lavaan::fitMeasures(fit_g, "df"))
  fit_g1 <- refit("z_lv ~~ a", ge)
  expect_equal(as.numeric(lavaan::fitMeasures(fit_g1, "df")),
               as.numeric(lavaan::fitMeasures(fit_g, "df")) - 1)
  # With equal residual variances, the variable's residual variance is held
  # equal across the schools, the latent variable's is not, and the
  # refusal asks for that equality: written with a label, the remedy is
  # the refused model.
  ge <- c("loadings", "residuals")
  fit_r2 <- lavaan::sem("f =~ y1 + y2 + y3\nz ~ f\nz ~~ a", data = d,
                        std.lv = TRUE, group = "school", group.equal = ge)
  scores_r2 <- factor_scores(fit_r2, data = d)
  expect_true(all(is.finite(scores_r2$.derived_regression)))
  expect_true(all(is.na(scores_r2$.derived_bartlett)))
  fit_rr <- refit("z_lv ~~ a\nz_lv ~~ c(r, r)*z_lv", ge)
  expect_equal(as.numeric(lavaan::logLik(fit_rr)),
               as.numeric(lavaan::logLik(fit_r2)), tolerance = 1e-8)
  expect_identical(lavaan::fitMeasures(fit_rr, "df"),
                   lavaan::fitMeasures(fit_r2, "df"))
})

# The oracle for a fit whose scores condition on observed variables: from
# a fit of the same model with a mean structure, the factor's conditional
# mean given every observed variable, E(f | y), and the properties the
# Bartlett score must have given the factor and the variables it conditions
# on (`cond`): with its weights A recovered from the scores by least
# squares, a slope of 1 on the factor and 0 on each conditioning variable,
# and a standard error of sqrt(A Var(y | f, cond) A').
fs_cond_oracle <- function(fit_ms, d, f, cond) {
  ca <- lavaan::lavInspect(fit_ms, "cov.all")
  ov <- lavaan::lavNames(fit_ms, "ov")
  y <- as.matrix(d[, ov])
  mu <- lavaan::lavInspect(fit_ms, "mean.ov")[ov]
  reg <- drop(sweep(y, 2L, mu) %*% solve(ca[ov, ov], ca[ov, f])) +
    lavaan::lavInspect(fit_ms, "mean.lv")[[f]]
  w <- c(f, cond)
  bartlett <- function(score) {
    a <- stats::coef(stats::lm(score ~ y))[-1L]
    v_res <- ca[ov, ov] - ca[ov, w] %*% solve(ca[w, w], ca[w, ov])
    list(slope = drop(a %*% ca[ov, w] %*% solve(ca[w, w])),
         se = sqrt(drop(a %*% v_res %*% a)))
  }
  list(regression = unname(reg), bartlett = bartlett,
       target = stats::setNames(c(1, rep(0, length(cond))), w))
}

test_that("an observed variable regressed on nothing is conditioned on, as an exogenous covariate is", {
  skip_if_not_installed("lavaan")
  # lavaan holds an observed variable that covaries with a factor (f ~~ z),
  # or a covariate that also covaries with another variable (f ~ x with
  # x ~~ a), as a latent variable regressed on nothing. The scores condition
  # on it, as on an exogenous covariate. Such fits were refused as if the
  # variable were an observed mediator; in a model with regressions,
  # lavPredict() gives the variable a row of zeros among the loadings and
  # misses the factor's conditional mean (by 0.03 and 1.9 below), and
  # without a mean structure it counts the variable's mean twice (1.06).
  hs <- holzinger_swineford
  d <- data.frame(y1 = hs$t1_visual_perception, y2 = hs$t2_cubes,
                  y3 = hs$t4_lozenges, y4 = hs$t6_paragraph_comprehension,
                  y5 = hs$t7_sentence, y6 = hs$t9_word_meaning,
                  z = hs$t3_paper_form_board, x = hs$age_years,
                  a = hs$t5_general_information)
  cases <- list(
    list(model = "f =~ y1 + y2 + y3\nf ~~ z", f = "f", cond = "z",
         lavaan_off = c(1.0, 0)),
    list(model = "f =~ y1 + y2 + y3\nf ~ x\nf ~~ z", f = "f",
         cond = c("x", "z"), lavaan_off = c(1.8, 1.8)),
    list(model = "f =~ y1 + y2 + y3\nf ~ x\nx ~~ a", f = "f",
         cond = c("x", "a"), lavaan_off = c(0.03, 0.03)),
    list(model = "f1 =~ y1 + y2 + y3\nf2 =~ y4 + y5 + y6\nf2 ~ f1\nf1 ~~ z",
         f = "f2", cond = c("f1", "z"), lavaan_off = c(2, 0.9)))
  for (cs in cases) {
    fit <- lavaan::sem(cs$model, data = d, std.lv = TRUE)
    fit_ms <- lavaan::sem(cs$model, data = d, std.lv = TRUE,
                          meanstructure = TRUE)
    o <- fs_cond_oracle(fit_ms, d, cs$f, cs$cond)
    for (k in 1:2) {
      ft <- list(fit, fit_ms)[[k]]
      off <- max(abs(lavaan::lavPredict(ft)[, cs$f] - o$regression))
      if (cs$lavaan_off[k] > 0) expect_gt(off, cs$lavaan_off[k])
      fs <- factor_scores(ft, data = d)
      k_f <- fs$factor == cs$f
      # The two fits' estimates differ by about 1e-6.
      expect_equal(fs$.derived_regression[k_f], o$regression,
                   tolerance = if (k == 1L) 1e-5 else 1e-10)
      b <- o$bartlett(fs$.derived_bartlett[k_f])
      expect_equal(b$slope, o$target, tolerance = 1e-6)
      expect_equal(fs$.derived_bartlett_se[k_f], rep(b$se, nrow(d)),
                   tolerance = 1e-6)
      expect_true(all(fs$reason == ""))
    }
  }
  # Under FIML each case is scored from its own observed variables, z among
  # them where it is observed: the regression score is the factor's mean
  # given them, and its standard error the posterior standard deviation
  # (lavPredict() misses the mean by 1.8 here).
  dm <- d[, c("y1", "y2", "y3", "z", "x")]
  set.seed(9)
  dm$y1[sample(nrow(dm), 40)] <- NA
  dm$z[sample(nrow(dm), 50)] <- NA
  dm$y3[sample(nrow(dm), 30)] <- NA
  fit_m <- lavaan::sem("f =~ y1 + y2 + y3\nf ~ x\nf ~~ z", data = dm,
                       std.lv = TRUE, missing = "ml")
  fm <- factor_scores(fit_m, data = dm)
  ca <- lavaan::lavInspect(fit_m, "cov.all")
  ov <- lavaan::lavNames(fit_m, "ov")
  mu <- lavaan::lavInspect(fit_m, "mean.ov")[ov]
  y <- as.matrix(dm[, ov])
  post <- t(vapply(seq_len(nrow(y)), function(i) {
    v <- ov[!is.na(y[i, ])]
    w <- solve(ca[v, v], ca[v, "f"])
    c(lavaan::lavInspect(fit_m, "mean.lv")[["f"]] +
        sum((y[i, v] - mu[v]) * w),
      sqrt(ca["f", "f"] - sum(ca["f", v] * w)))
  }, numeric(2)))
  expect_gt(max(abs(lavaan::lavPredict(fit_m)[, "f"] - post[, 1])), 1.5)
  expect_equal(fm$.derived_regression, post[, 1], tolerance = 1e-10)
  expect_equal(fm$.derived_regression_se, post[, 2], tolerance = 1e-10)
  # Graham's saturated correlate in a MIMIC model, which also covaries with
  # the covariate: lavaan then holds the covariate and the auxiliary among
  # the latent variables. The scores from the indicators are those of the
  # model without the auxiliary, every parameter fixed at the fit's
  # estimate.
  dm <- d[, c("y1", "y2", "y3", "z", "x")]
  set.seed(4)
  dm$y3[sample(nrow(dm), 60)] <- NA
  dm$y1[sample(nrow(dm), 50)] <- NA
  fit_a <- lavaan::sem("f =~ y1 + y2 + y3\nf ~ x\nz ~~ y1 + y2 + y3 + x",
                       data = dm, std.lv = TRUE, missing = "ml")
  expect_identical(lavaan::lavNames(fit_a, "ov.x"), character(0))
  fa <- suppressMessages(factor_scores(fit_a, data = dm))
  expect_identical(attr(fa, "auxiliary"), "z")
  expect_true(all(fa$reason == ""))
  pa <- lavaan::parTable(fit_a)
  pa <- pa[!(pa$lhs == "z" | pa$rhs == "z"), ]
  pa$free <- 0L
  pa$ustart <- pa$est
  pa$start <- pa$est
  pa$se <- NULL
  fixed <- lavaan::lavaan(pa, data = dm[, c("y1", "y2", "y3", "x")],
                          missing = "ml")
  for (m in c("regression", "Bartlett")) {
    lp <- lavaan::lavPredict(fixed, method = m, se = "standard")
    col <- paste0(".derived_", tolower(m))
    expect_equal(fa[[col]], unname(lp[, "f"]), tolerance = 1e-10)
    expect_equal(fa[[paste0(col, "_se")]], unname(attr(lp, "se")[[1]][, "f"]),
                 tolerance = 1e-10)
  }
  # Conditional on the covariates, with regressions, lavPredict() misses
  # the conditional means (by 1.8 here), and the scores computed here
  # assume conditional.x = FALSE: the fit is refused, and refitted so it is
  # scored exactly.
  m_c <- "f1 =~ y1 + y2 + y3\nf2 =~ y4 + y5 + y6\nf2 ~ f1 + x\nf1 ~~ z"
  fit_c <- lavaan::sem(m_c, data = d, std.lv = TRUE, conditional.x = TRUE)
  fit_u <- lavaan::sem(m_c, data = d, std.lv = TRUE, meanstructure = TRUE)
  o <- fs_cond_oracle(fit_u, d, "f1", c("f2", "x", "z"))
  expect_gt(max(abs(lavaan::lavPredict(fit_c)[, "f1"] - o$regression)), 1.5)
  expect_error(factor_scores(fit_c, data = d), paste0(
    "regresses on nothing (z): lavaan::lavPredict() then gives such a ",
    "variable a row of zeros among the loadings"), fixed = TRUE)
  fu <- factor_scores(fit_u, data = d)
  expect_equal(fu$.derived_regression[fu$factor == "f1"], o$regression,
               tolerance = 1e-10)
  # Without regressions lavaan keeps the variable's row, and its scores of a
  # conditional fit are kept.
  fit_c0 <- lavaan::sem("f =~ y1 + y2 + y3\nf ~~ z\nf ~ x", data = d,
                        std.lv = TRUE, conditional.x = TRUE)
  fc0 <- factor_scores(fit_c0, data = d)
  expect_identical(fc0$.derived_regression, unname(lavaan::lavPredict(fit_c0)[, "f"]))
  # Only Bartlett's stronger conditional-unbiasedness condition fails.
  fit_a <- lavaan::sem("f =~ y1 + y2 + y3\nz ~ f\nz ~~ a",
                        data = d, std.lv = TRUE, meanstructure = TRUE)
  fa <- factor_scores(fit_a, data = d)
  oracle_a <- fs_cond_oracle(fit_a, d, "f", "a")
  expect_equal(fa$.derived_regression, oracle_a$regression, tolerance = 1e-10)
  expect_true(all(is.na(fa$.derived_bartlett)))
  expect_true(all(grepl("Bartlett scores withheld", fa$reason)))
  # The score columns are named by role in the derived_columns attribute,
  # and the unprefixed names are not columns.
  map <- attr(fa, "derived_columns")
  expect_identical(names(map),
                   c("regression", "regression_se", "regression_ci_lower",
                     "regression_ci_upper", "bartlett", "bartlett_se",
                     "bartlett_ci_lower", "bartlett_ci_upper"))
  expect_identical(unname(map), paste0(".derived_", names(map)))
  expect_true(all(map %in% names(fa)))
  expect_false(any(names(map) %in% names(fa)))

})

# The oracle for a fit conditional on the covariates without a B matrix,
# from its own estimates: per case and its observed variables o, with Psi
# the latent covariance matrix given the covariates (every observed
# variable lavaan holds as latent keeping its unit row among the loadings),
# Sigma_oo = Lambda_o Psi Lambda_o' + Theta_oo, and E(eta | x) = alpha +
# Gamma x, the regression score E(f | y_o, x), and the Bartlett score with
# weights A = (Lambda_h' Sigma_oo^{-1} Lambda_h)^{-1} Lambda_h'
# Sigma_oo^{-1} over the latent variables with an observed indicator (h),
# whose row of A Lambda_h for f is its slope on (f, z) given x.
fs_gamma_oracle <- function(fit, d, f, group = NULL) {
  est <- lavaan::lavInspect(fit, "est")
  if (!is.null(est$lambda)) est <- list(est)
  gi <- if (is.null(group)) rep(1L, nrow(d)) else {
    match(d[[group]], lavaan::lavInspect(fit, "group.label"))
  }
  out <- t(vapply(seq_len(nrow(d)), function(i) {
    e <- est[[gi[i]]]
    lam <- e$lambda
    v <- e$psi
    s <- lam %*% v %*% t(lam) + e$theta
    eeta <- drop(e$alpha) + drop(e$gamma %*% d$x[i])
    ey <- drop(e$nu) + drop(lam %*% eeta)
    y <- unlist(d[i, rownames(lam)])
    o <- !is.na(y)
    lo <- lam[o, , drop = FALSE]
    so <- s[o, o, drop = FALSE]
    lh <- lo[, colSums(lo != 0) > 0, drop = FALSE]
    a <- solve(t(lh) %*% solve(so, lh), t(lh) %*% solve(so))
    sl <- c(f = 0, z = 0)
    sl[colnames(lh)] <- (a %*% lh)[f, ]
    c(eeta[[f]] + sum((y[o] - ey[o]) * solve(so, lo %*% v[, f])),
      eeta[[f]] + sum(a[f, ] * (y[o] - ey[o])), sl)
  }, numeric(4)))
  list(regression = out[, 1L], bartlett = out[, 2L], slope = out[, 3:4])
}

test_that("conditional on the covariates and without B, a variable regressed only on them is conditioned on", {
  skip_if_not_installed("lavaan")
  # In a fit conditional on the covariates with no regression among the
  # variables lavaan holds as latent (no B matrix), lavaan keeps the unit
  # row among the loadings of an observed variable regressed only on a
  # covariate (z ~ x, through Gamma), as of one regressed on nothing, so
  # lavPredict()'s scores condition on it: its regression score is
  # E(f | y, z, x) and its Bartlett score has slope 1 on f and 0 on z. Such
  # a fit was refused as if lavaan measured the latent variables by z,
  # whenever z covaried with the factor, which lavaan::sem() frees by
  # default when the covariate also predicts the factor.
  hs <- holzinger_swineford
  d <- data.frame(y1 = hs$t1_visual_perception, y2 = hs$t2_cubes,
                  y3 = hs$t4_lozenges, y4 = hs$t6_paragraph_comprehension,
                  y5 = hs$t7_sentence, y6 = hs$t9_word_meaning,
                  z = hs$t3_paper_form_board, x = hs$age_years,
                  school = hs$school)
  m <- "f =~ y1 + y2 + y3\nf ~ x\nz ~ x"
  fit <- lavaan::sem(m, data = d, std.lv = TRUE, conditional.x = TRUE)
  pt <- lavaan::parTable(fit)
  expect_true(any(pt$lhs == "f" & pt$op == "~~" & pt$rhs == "z" &
                    pt$free > 0 & pt$est != 0))
  expect_null(lavaan::lavInspect(fit, "est")$beta)
  one <- matrix(c(1, 0), nrow(d), 2L, byrow = TRUE,
                dimnames = list(NULL, c("f", "z")))
  o <- fs_gamma_oracle(fit, d, "f")
  expect_equal(o$slope, one, tolerance = 1e-10)
  fs <- factor_scores(fit, data = d)
  k <- fs$factor == "f"
  expect_equal(fs$.derived_regression[k], o$regression, tolerance = 1e-10)
  expect_equal(fs$.derived_bartlett[k], o$bartlett, tolerance = 1e-10)
  expect_true(all(is.finite(fs$.derived_bartlett_se[k])))
  # The same model with z carried as a single indicator latent variable
  # (the same log likelihood) gives the same scores, to the optimizer's
  # tolerance.
  fit_r <- lavaan::sem("f =~ y1 + y2 + y3\nf ~ x\nZ =~ z\nz ~~ 0*z\nZ ~ x
                        f ~~ Z", data = d, std.lv = TRUE,
                       conditional.x = TRUE)
  expect_equal(lavaan::logLik(fit_r), lavaan::logLik(fit), tolerance = 1e-8)
  fr <- factor_scores(fit_r, data = d)
  expect_equal(fr$.derived_regression[fr$factor == "f"], fs$.derived_regression[k],
               tolerance = 1e-5)
  expect_equal(fr$.derived_bartlett[fr$factor == "f"], fs$.derived_bartlett[k],
               tolerance = 1e-5)
  # Under FIML, each case from its own observed variables (z among them
  # where it is observed), and with two groups.
  dm <- d
  set.seed(9)
  dm$y1[sample(nrow(dm), 40)] <- NA
  dm$z[sample(nrow(dm), 50)] <- NA
  dm$y3[sample(nrow(dm), 30)] <- NA
  # suppressWarnings: lavaan's FIML with conditional.x partially matches
  # `cov` to `cov.x` and warns under warnPartialMatchDollar; not a property
  # under test.
  fit_m <- suppressWarnings(lavaan::sem(m, data = dm, std.lv = TRUE,
                                        conditional.x = TRUE, missing = "ml"))
  om <- fs_gamma_oracle(fit_m, dm, "f")
  fm <- factor_scores(fit_m, data = dm)
  km <- fm$factor == "f"
  expect_equal(fm$.derived_regression[km], om$regression, tolerance = 1e-10)
  expect_equal(fm$.derived_bartlett[km], om$bartlett, tolerance = 1e-10)
  expect_equal(om$slope[, "f"], rep(1, nrow(dm)), tolerance = 1e-10)
  fit_g <- lavaan::sem(m, data = d, std.lv = TRUE, conditional.x = TRUE,
                       group = "school")
  og <- fs_gamma_oracle(fit_g, d, "f", group = "school")
  fg <- factor_scores(fit_g, data = d)
  kg <- fg$factor == "f"
  expect_equal(fg$.derived_regression[kg], og$regression, tolerance = 1e-10)
  expect_equal(fg$.derived_bartlett[kg], og$bartlett, tolerance = 1e-10)
  # With a B matrix lavaan measures the latent variables by z, through its
  # row of B, and the fit is still refused: not conditional on the
  # covariates (x then a latent variable of its own), or with a regression
  # among the latent variables.
  for (fit_b in list(
    lavaan::sem(m, data = d, std.lv = TRUE),
    lavaan::sem(paste(m, "\nf2 =~ y4 + y5 + y6\nf2 ~ f"), data = d,
                std.lv = TRUE, conditional.x = TRUE))) {
    expect_false(is.null(lavaan::lavInspect(fit_b, "est")$beta))
    expect_error(factor_scores(fit_b, data = d), paste0(
      "The disturbance of z, an observed variable regressed on other ",
      "variables, covaries with f"), fixed = TRUE)
  }
})

test_that("missingness through an auxiliary variable: the designs the help quotes", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two fits of 5,000 cases
  # The help quoted two cases as one design while their auxiliaries
  # differed (one correlating .7 with the factor, the other .4 with the
  # factor and .5 with the first item's residual). Both auxiliaries now
  # correlate .4 with the factor: four items loading .7, the last three
  # missing whenever the auxiliary falls below -0.3, the auxiliary carried
  # as a saturated correlate under FIML, and the cases left with one item.
  aux_design <- function(rcor) {
    set.seed(113)
    n <- 5000
    eta <- stats::rnorm(n)
    e <- matrix(stats::rnorm(4 * n, sd = sqrt(0.51)), n, 4)
    d <- as.data.frame(0.7 * eta + e)
    names(d) <- paste0("y", 1:4)
    d$a <- 0.4 * eta + rcor * e[, 1] / sqrt(0.51) +
      stats::rnorm(n, sd = sqrt(1 - 0.16 - rcor^2))
    one <- d$a < -0.3
    d[one, c("y2", "y3", "y4")] <- NA
    # lavaan's post-estimation check finds the residual covariance matrix
    # with the auxiliary not positive definite, as a saturated correlate
    # related to the factor makes it; the indicators' block, which the
    # scores use, is (factor_scores() gives no warning).
    fit <- suppressWarnings(lavaan::cfa(
      "f =~ y1 + y2 + y3 + y4\na ~~ y1 + y2 + y3 + y4", data = d,
      std.lv = TRUE, missing = "ml"))
    expect_no_warning(fs <- suppressMessages(factor_scores(fit)))
    err <- (fs$.derived_bartlett - eta)[one]
    c(bartlett = mean(err), mc_se = stats::sd(err) / sqrt(sum(one)),
      cover = mean((fs$.derived_bartlett_ci_lower <= eta &
                      eta <= fs$.derived_bartlett_ci_upper)[one]),
      regression = mean((fs$.derived_regression - eta)[one]))
  }
  through <- aux_design(0)
  expect_equal(unname(through[c("bartlett", "cover", "regression")]),
               c(-0.046, 0.956, 0.170), tolerance = 0.01)
  expect_lt(abs(through[["bartlett"]]), 2 * through[["mc_se"]])
  shared <- aux_design(0.5)
  expect_equal(unname(shared[c("bartlett", "cover")]), c(-0.562, 0.931),
               tolerance = 0.01)
})

test_that("auxiliary variables are set aside in every group", {
  skip_if_not_installed("lavaan")
  d <- fs_mar_data(c("school", "age_years"))
  d$age_years <- as.numeric(scale(d$age_years))
  aux <- paste("age_years ~~", paste(unlist(fs_items), collapse = " + "))
  # lavaan warns that the two rows with every item missing are ignored.
  fit <- suppressWarnings(lavaan::cfa(paste(fs_model, "\n", aux), data = d,
                                      std.lv = TRUE, missing = "ml",
                                      group = "school"))
  expect_message(fs <- factor_scores(fit, data = d), "saturated correlates")
  # The oracle: lavaan's own scores of the model without the auxiliary,
  # every parameter of each group fixed at the auxiliary fit's estimate.
  pt <- lavaan::parTable(fit)
  pt <- pt[!(pt$lhs == "age_years" | pt$rhs == "age_years"), ]
  pt$free <- 0L
  pt$ustart <- pt$est
  pt$start <- pt$est
  pt$se <- NULL
  fixed <- suppressWarnings(lavaan::lavaan(
    pt, data = d[, c(unlist(fs_items), "school")], missing = "ml",
    group = "school"))
  idx <- lavaan::lavInspect(fixed, "case.idx")
  for (m in c("regression", "Bartlett")) {
    lp <- lavaan::lavPredict(fixed, method = m, se = "standard",
                             assemble = FALSE)
    col <- paste0(".derived_", tolower(m))
    for (g in seq_along(idx)) {
      for (f in names(fs_items)) {
        rows <- fs$factor == f
        expect_equal(fs[[col]][rows][idx[[g]]], unname(lp[[g]][, f]),
                     tolerance = 1e-10)
        expect_equal(fs[[paste0(col, "_se")]][rows][idx[[g]]],
                     unname(attr(lp, "se")[[g]][, f]), tolerance = 1e-10)
      }
    }
  }
})

test_that("a pattern whose items do not separate two factors has no Bartlett score for them", {
  skip_if_not_installed("lavaan")
  # Under FIML a case that answered only y3, which loads on both factors,
  # has items that do not identify either factor's value. lavaan forms the
  # Bartlett weights from a generalized inverse there, and its minimum
  # norm "score" and "standard error" (0.18, against 1.44 for y3 as a
  # measure of f2 alone) were passed through with an interval.
  set.seed(113)
  n <- 500
  L <- matrix(c(.7, .7, .5, .4, 0, 0, 0, 0, .4, .5, .7, .7), 6, 2)
  eta <- matrix(stats::rnorm(2 * n), n, 2) %*%
    chol(matrix(c(1, .4, .4, 1), 2))
  d <- as.data.frame(eta %*% t(L) + matrix(stats::rnorm(6 * n, sd = .6), n, 6))
  names(d) <- paste0("y", 1:6)
  d[1:5, c(1, 2, 4, 5, 6)] <- NA
  model <- "f1 =~ y1 + y2 + y3 + y4\nf2 =~ y3 + y4 + y5 + y6"
  fit <- lavaan::cfa(model, data = d, std.lv = TRUE, missing = "ml")
  fs <- factor_scores(fit, data = d)
  lp <- lavaan::lavPredict(fit, method = "Bartlett")
  expect_true(all(is.finite(lp[1:5, ])))
  one <- fs$row %in% 1:5
  expect_true(all(is.na(fs$.derived_bartlett[one])))
  expect_true(all(is.na(fs$.derived_bartlett_se[one])))
  expect_true(all(is.na(fs$.derived_bartlett_ci_lower[one])))
  expect_true(all(is.na(fs$determinacy_bartlett[one])))
  expect_true(all(fs$reason[one] == paste0(
    "the observed items do not separate this factor from another: no ",
    "Bartlett score")))
  expect_true(all(is.finite(fs$.derived_regression[one])))
  # Every other row is lavaan's, and the pattern table says what the rows
  # say.
  for (f in c("f1", "f2")) {
    expect_identical(fs$.derived_bartlett[fs$factor == f][-(1:5)],
                     unname(lp[-(1:5), f]))
  }
  pt <- attr(fs, "patterns")
  m <- merge(fs[, c("pattern", "factor", 
                  ".derived_bartlett_se")],
             pt[, c("pattern", "factor", "bartlett_se")],
             by = c("pattern", "factor"))
  expect_equal(m$.derived_bartlett_se, m$bartlett_se)
  # Scores computed here (an auxiliary variable set aside) say the same.
  # lavaan's post-estimation check warns that the residual covariance
  # matrix with z in it is not positive definite: z is mostly y1, which
  # also measures f1. The indicators' own is, so factor_scores(), which sets
  # z aside, does not warn.
  d$z <- d$y1 + stats::rnorm(n)
  fit_z <- suppressWarnings(lavaan::cfa(
    paste(model, "\nz ~~ y1 + y2 + y3 + y4 + y5 + y6"), data = d,
    std.lv = TRUE, missing = "ml"))
  expect_false(suppressWarnings(lavaan::lavInspect(fit_z, "post.check")))
  expect_message(fs_z <- factor_scores(fit_z, data = d),
                 "saturated correlates")
  expect_true(all(is.na(fs_z$.derived_bartlett[one])))
  expect_identical(fs_z$reason[one], fs$reason[one])
  # A factor the pattern does identify keeps its score. With an item of f1
  # alone and an item shared by f2 and f3, f1 is identified and the other
  # two are not; f1's score is then that of its one item, with that item's
  # error, which lavaan's generalized inverse gives uniquely.
  L3 <- cbind(c(.7, .7, .6, 0, 0, 0, 0, 0), c(0, 0, 0, .7, .7, .5, 0, 0),
              c(0, 0, 0, 0, 0, .5, .7, .7))
  eta3 <- matrix(stats::rnorm(3 * n), n, 3)
  d3 <- as.data.frame(eta3 %*% t(L3) + matrix(stats::rnorm(8 * n, sd = .6), n, 8))
  names(d3) <- paste0("x", 1:8)
  d3[1:3, c(2, 3, 4, 5, 7, 8)] <- NA
  fit3 <- lavaan::cfa("a =~ x1 + x2 + x3\nb =~ x4 + x5 + x6\nc =~ x6 + x7 + x8",
                      data = d3, std.lv = TRUE, missing = "ml")
  fs3 <- factor_scores(fit3, data = d3)
  est <- lavaan::lavInspect(fit3, "est")
  mu <- lavaan::lavInspect(fit3, "implied")$mean
  a <- fs3$factor == "a" & fs3$row %in% 1:3
  lam1 <- est$lambda["x1", "a"]
  expect_equal(fs3$.derived_bartlett[a], unname((d3$x1[1:3] - mu[["x1"]]) / lam1),
               tolerance = 1e-8)
  expect_equal(fs3$.derived_bartlett_se[a],
               rep(sqrt(est$theta["x1", "x1"]) / lam1, 3), tolerance = 1e-8)
  expect_equal(fs3$.derived_bartlett[a],
               unname(lavaan::lavPredict(fit3, method = "Bartlett")[1:3, "a"]),
               tolerance = 1e-8)
  bc <- fs3$factor %in% c("b", "c") & fs3$row %in% 1:3
  expect_true(all(is.na(fs3$.derived_bartlett[bc])))
  expect_true(all(grepl("do not separate", fs3$reason[bc])))
})

test_that("ordered items give EBM and ML scores; extreme patterns set aside", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # EBM and ML scores optimize case by case
  its <- c(fs_items$visual, fs_items$textual)
  db <- as.data.frame(lapply(holzinger_swineford[1:150, its], function(x) {
    cut(x, stats::quantile(x, c(0, .25, .5, .75, 1)),
        include.lowest = TRUE, labels = FALSE)
  }))
  db$t2_cubes[1:10] <- NA
  db[20, fs_items$textual] <- NA
  db[21, ] <- NA
  model <- "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
            textual =~ t6_paragraph_comprehension + t7_sentence +
                        t9_word_meaning"
  # lavaan warns about the empty row.
  fit <- suppressWarnings(lavaan::cfa(model, data = db, std.lv = TRUE,
                                      ordered = TRUE, missing = "pairwise"))
  fs <- factor_scores(fit, data = db)
  expect_fs_shape(fs)
  expect_identical(attr(fs, "methods"), c(regression = "EBM", bartlett = "ML"))
  ebm <- lavaan::lavPredict(fit, method = "EBM")
  idx <- lavaan::lavInspect(fit, "case.idx")
  expect_equal(fs$.derived_regression[fs$factor == "visual"][idx],
               unname(ebm[, "visual"]))
  # Every observed visual item at the same end of its scale: no finite ML
  # score; the EBM score stays.
  v <- db[, fs_items$visual]
  ends <- apply(v, 1L, function(r) {
    r <- r[!is.na(r)]
    length(r) > 0L && (all(r == 1) || all(r == 4))
  })
  vis <- fs$factor == "visual"
  expect_true(any(ends))
  expect_true(all(is.na(fs$.derived_bartlett[vis][ends])))
  expect_true(all(grepl("end of its scale", fs$reason[vis][ends])))
  expect_true(all(is.finite(fs$.derived_regression[vis][ends])))
  ml <- lavaan::lavPredict(fit, method = "ML")
  bv <- fs$.derived_bartlett[vis][idx]
  keep <- !ends[idx] & is.finite(bv)
  expect_equal(bv[keep], unname(ml[keep, "visual"]))
  # The standard errors are lavaan's too.
  ebm_se <- attr(lavaan::lavPredict(fit, method = "EBM", se = "standard"),
                 "se")[[1]]
  ml_se <- attr(lavaan::lavPredict(fit, method = "ML", se = "standard"),
                "se")[[1]]
  expect_equal(fs$.derived_regression_se[vis][idx], unname(ebm_se[, 1]))
  expect_equal(fs$.derived_bartlett_se[vis][idx][keep], unname(ml_se[keep, 1]))
  # Determinacy is not computed for ordered items, and every row says so.
  expect_true(all(is.na(fs$determinacy_regression)))
  scored <- fs$items_observed > 0L & !(fs$row == 21)
  expect_true(all(grepl("not computed for ordered items",
                        fs$reason[scored])))
  # No textual item: no ML score, the EBM score from the visual items.
  tx <- fs$row == 20 & fs$factor == "textual"
  expect_true(is.na(fs$.derived_bartlett[tx]))
  expect_true(is.finite(fs$.derived_regression[tx]))
  expect_true(all(fs$reason[fs$row == 21] == "no observed item"))
})

test_that("the ML check finds every direction in which the likelihood rises without bound", {
  # Two factors, y1 (two categories) on f1, y2 (three) on f2, and y3
  # (three) on both. A row's ML score of a factor is finite and unique when
  # no direction of the factor values that lowers no observed item's
  # probability moves the factor.
  lam <- matrix(c(1, 0, .8, 0, 1, .6), 3, 2,
                dimnames = list(c("y1", "y2", "y3"), c("f1", "f2")))
  tau <- matrix(c(-.5, -.4, .6, 0, 1), ncol = 1,
                dimnames = list(c("y1|t1", "y2|t1", "y2|t2", "y3|t1",
                                  "y3|t2"), "threshold"))
  x <- rbind(c(2, 1, 2),    # f1 up and f2 down by .6 to .8 keep y3's mean
             c(2, 3, 2),    # no such direction: both finite
             c(2, NA, 2),   # the same direction as the first row
             c(NA, NA, 2),  # y3 alone cannot separate them
             c(1, 1, 1),    # every item at its lowest category
             c(2, 1, NA),   # simple structure in effect: both at an end
             c(2, 2, 3),    # f1 at an end; f2 finite, though it shares y3
             c(2, 1, 3),    # f1 at an end; f2 down with f1 up raises all
             c(NA, 2, NA))  # f1 unobserved; f2 finite
  colnames(x) <- rownames(lam)
  code <- .fs_ml_finite(x, list(lambda = lam, tau = tau), c("f1", "f2"),
                        rownames(lam), c(TRUE, TRUE))
  expect_identical(unname(code), rbind(
    c("joint", "joint"), c("", ""), c("joint", "joint"), c("flat", "flat"),
    c("extreme", "extreme"), c("extreme", "extreme"), c("extreme", ""),
    c("extreme", "joint"), c("", "")))
  # A factor without items of its own is not judged here.
  expect_identical(unname(.fs_ml_finite(x, list(lambda = lam, tau = tau),
                                        c("f1", "f2"), rownames(lam),
                                        c(TRUE, FALSE))[1, ]),
                   c("joint", ""))
  # The cone membership behind it (Farkas' lemma).
  b <- rbind(c(1, 0), c(0, 1))
  expect_true(.fs_in_cone(b, c(1, 2)))
  expect_false(.fs_in_cone(b, c(-1, 0)))
  expect_true(.fs_in_cone(rbind(c(1, 0), c(-1, 0), c(0, 1)), c(-2, 3)))
  expect_true(.fs_in_cone(rbind(c(1, 1), c(1, 1), c(2, 2)), c(1, 1)))
  expect_false(.fs_in_cone(rbind(c(1, 1), c(1, 1), c(2, 2)), c(1, 0)))
  expect_false(.fs_in_cone(matrix(0, 0, 2), c(1, 0)))
  # Identification of a factor by a pattern's items.
  expect_identical(.fs_identified(matrix(c(.4, .5), 1)), c(FALSE, FALSE))
  expect_identical(.fs_identified(rbind(c(1, 0, 0), c(0, 1, 1))),
                   c(TRUE, FALSE, FALSE))
  expect_identical(.fs_identified(rbind(c(1, 0, 0), c(1, 1, 0))),
                   c(TRUE, TRUE, FALSE))
})

test_that("ordered items with a cross-loading: no ML score where the likelihood has no finite maximum", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # EBM and ML scores optimize case by case
  # y4 loads on both factors. The check of extreme patterns used to skip
  # any factor with such an item, so rows 71 and 123, every item in its
  # lowest and its highest category, kept lavaan's ML scores (where its
  # ridge stopped, with standard errors of about 71 and 316) and intervals.
  its <- c(fs_items$visual, fs_items$textual)
  db <- as.data.frame(lapply(holzinger_swineford[, its], function(x) {
    cut(x, stats::quantile(x, 0:4 / 4), include.lowest = TRUE,
        labels = FALSE)
  }))
  names(db) <- paste0("y", 1:6)
  fit <- lavaan::cfa("f1 =~ y1 + y2 + y3 + y4\nf2 =~ y4 + y5 + y6",
                     data = db, std.lv = TRUE, ordered = TRUE)
  fs <- factor_scores(fit, data = db)
  ends <- fs$row %in% c(71, 123)
  expect_true(all(is.na(fs$.derived_bartlett[ends])))
  expect_true(all(grepl("end of its scale", fs$reason[ends])))
  expect_true(all(is.finite(fs$.derived_regression[ends])))
  # No ML score kept has a standard error the size of lavaan's ridge.
  expect_lt(max(fs$.derived_bartlett_se, na.rm = TRUE), 10)
  # Row 47 answered y1 to y3 in their highest category and y4 to y6 in
  # their lowest: raising f1 while lowering f2 raises every probability,
  # though y4, which loads on both, keeps f1 from being at an end alone.
  r47 <- fs$row == 47
  expect_identical(fs$reason[r47 & fs$factor == "f1"], paste0(
    "no finite ML score: the likelihood rises without bound as this ",
    "factor moves with others that share its items; determinacy is not ",
    "computed for ordered items"))
  # Every ML score lavaan gives and the check keeps is lavaan's.
  ml <- lavaan::lavPredict(fit, method = "ML")
  for (f in c("f1", "f2")) {
    b <- fs$.derived_bartlett[fs$factor == f]
    keep <- is.finite(b) & is.finite(ml[, f])
    expect_identical(b[keep], unname(ml[keep, f]))
  }
})

# The ordered probit likelihood of each case of a single group fit of
# ordered items (delta parameterization, scale factors of 1), maximized
# directly: the brute force oracle of the ML scores. Each item's latent
# response mean is lambda_eff'eta, with an item regressed on the factors
# (held by lavaan as a dummy latent variable) loading through its row of B.
# A ridge of 1e-12 keeps a factor without a finite maximum from running
# off, so the factors with one reach it; the standard errors are from the
# Hessian of the unridged log likelihood there.
fs_ml_brute <- function(fit, lv, rows) {
  est <- lavaan::lavInspect(fit, "est")
  lam <- est$lambda
  dum <- intersect(setdiff(colnames(lam), lv), rownames(lam))
  if (length(dum)) lam[dum, ] <- est$beta[dum, colnames(lam)]
  lam <- lam[, lv, drop = FALSE]
  sdv <- sqrt(diag(lavaan::lavInspect(fit, "theta")))[rownames(lam)]
  y <- lavaan::lavInspect(fit, "data")[, rownames(lam), drop = FALSE]
  cuts <- lapply(rownames(lam), function(j) {
    c(-Inf, est$tau[startsWith(rownames(est$tau), paste0(j, "|")), 1], Inf)
  })
  nll <- function(e, yi, r) {
    mu <- drop(lam %*% e)
    s <- 0
    for (j in which(!is.na(yi))) {
      a <- (cuts[[j]][yi[j]] - mu[j]) / sdv[j]
      b <- (cuts[[j]][yi[j] + 1L] - mu[j]) / sdv[j]
      s <- s + if (is.infinite(a)) {
        stats::pnorm(b, log.p = TRUE)
      } else if (is.infinite(b)) {
        stats::pnorm(-a, log.p = TRUE)
      } else {
        log(stats::pnorm(b) - stats::pnorm(a))
      }
    }
    r / 2 * sum(e^2) - s
  }
  out <- t(vapply(rows, function(i) {
    best <- NULL
    for (st in c(0, 1, -1)) {
      o <- stats::nlminb(rep(st, length(lv)), nll, yi = y[i, ], r = 1e-12,
                         control = list(rel.tol = 1e-14, x.tol = 1e-12,
                                        iter.max = 5000, eval.max = 1e4))
      if (is.null(best) || o$objective < best$objective) best <- o
    }
    h <- stats::optimHess(best$par, nll, yi = y[i, ], r = 0)
    se <- suppressWarnings(sqrt(diag(MASS::ginv(h))))
    c(best$par, se)
  }, numeric(2L * length(lv))))
  list(score = out[, seq_along(lv), drop = FALSE],
       se = out[, length(lv) + seq_along(lv), drop = FALSE])
}

test_that("ordered items: a factor sharing items with one without a finite ML score keeps its own", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # EBM and ML scores optimize case by case
  # f1 and f2 share y4. When f2's items all sit at one end, f2's likelihood
  # rises without bound, but f1's maximum is finite and unique: y4's
  # probability tends to 1 as f2 moves out, whatever f1's value. The ML
  # scores of f1 in those 65 rows were withheld ("it shares observed items
  # with a factor that has no finite ML score"), although lavaan's values
  # agree with the maximum to 1.2e-4, and lavaan's default optimizer stops
  # on its iteration limit in the 6 rows of f2 sharing y4 with an f1 at an
  # end.
  its <- c(fs_items$visual, fs_items$textual)
  db <- as.data.frame(lapply(holzinger_swineford[, its], function(x) {
    cut(x, stats::quantile(x, 0:4 / 4), include.lowest = TRUE,
        labels = FALSE)
  }))
  names(db) <- paste0("y", 1:6)
  fit <- lavaan::cfa("f1 =~ y1 + y2 + y3 + y4\nf2 =~ y4 + y5 + y6",
                     data = db, std.lv = TRUE, ordered = TRUE)
  fs <- factor_scores(fit, data = db)
  div <- grepl("end of its scale|rises without bound", fs$reason)
  b <- matrix(fs$.derived_bartlett, ncol = 2L)
  se <- matrix(fs$.derived_bartlett_se, ncol = 2L)
  d <- matrix(div, ncol = 2L)
  expect_true(all(is.na(b[d])))
  shared <- cbind(d[, 2L] & !d[, 1L], d[, 1L] & !d[, 2L])
  expect_identical(colSums(shared), c(65, 6))
  expect_true(all(is.finite(b[shared])))
  expect_true(all(is.finite(se[shared])))
  expect_false(any(grepl("shares observed items", fs$reason)))
  rows <- which(shared[, 1L] | shared[, 2L])
  bf <- fs_ml_brute(fit, c("f1", "f2"), rows)
  expect_lt(max(abs(b[rows, ][shared[rows, ]] - bf$score[shared[rows, ]])),
            2.5e-4)
  expect_lt(max(abs(se[rows, ][shared[rows, ]] - bf$se[shared[rows, ]])),
            1e-4)
  # f1's scores are lavaan's; f2's, where lavaan's default optimizer gave
  # none, come from its nlminb optimizer.
  ml <- lavaan::lavPredict(fit, method = "ML")
  expect_identical(b[shared[, 1L], 1L], unname(ml[shared[, 1L], "f1"]))
  expect_true(all(is.na(ml[shared[, 2L], "f2"])))
  nl <- lavaan::lavPredict(fit, method = "ML", optim_method = "nlminb")
  expect_identical(b[shared[, 2L], 2L], unname(nl[shared[, 2L], "f2"]))
})

test_that("ordered items: an item regressed on the factor bounds the ML score", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # EBM and ML scores optimize case by case
  # The same three items, y1 written as an outcome the factor predicts
  # rather than as its indicator: the same likelihood, which lavaan holds
  # with y1's loading in B. The check of the ML score read the loadings from
  # Lambda alone, so the 22 rows with y2 and y3 at one end and y1 in a
  # middle category were told that every item of the factor sat at an end,
  # and had no ML score.
  its <- c("t1_visual_perception", "t2_cubes", "t4_lozenges")
  db <- as.data.frame(lapply(holzinger_swineford[, its], function(x) {
    cut(x, stats::quantile(x, 0:4 / 4), include.lowest = TRUE,
        labels = FALSE)
  }))
  names(db) <- c("y1", "y2", "y3")
  fit_i <- lavaan::cfa("f =~ y1 + y2 + y3", data = db, std.lv = TRUE,
                       ordered = TRUE)
  fit_r <- lavaan::sem("f =~ y2 + y3\ny1 ~ f", data = db, std.lv = TRUE,
                       ordered = TRUE)
  expect_identical(unname(lavaan::lavInspect(fit_r, "est")$lambda["y1", "f"]),
                   0)
  fs_i <- factor_scores(fit_i, data = db)
  fs_r <- factor_scores(fit_r, data = db)
  expect_identical(is.na(fs_r$.derived_bartlett), is.na(fs_i$.derived_bartlett))
  expect_identical(fs_r$reason, fs_i$reason)
  expect_equal(fs_r$.derived_bartlett, fs_i$.derived_bartlett, tolerance = 1e-5)
  mid <- db$y2 == db$y3 & db$y2 %in% c(1, 4) & db$y1 %in% 2:3
  expect_identical(sum(mid), 22L)
  expect_true(all(is.finite(fs_r$.derived_bartlett[mid])))
  bf <- fs_ml_brute(fit_r, "f", which(mid))
  expect_equal(fs_r$.derived_bartlett[mid], bf$score[, 1L], tolerance = 1e-3)
  expect_equal(fs_r$.derived_bartlett_se[mid], bf$se[, 1L], tolerance = 1e-3)
})

test_that("ordered items: a direct effect of a covariate on an item is refused", {
  skip_if_not_installed("lavaan")
  # lavaan's EBM and ML scores leave the direct effect (y1 ~ x) out of each
  # case's likelihood, while the thresholds are estimated given x, so they
  # are the scores of a case with x at zero: the same for every age within
  # a response pattern.
  its <- c("t1_visual_perception", "t2_cubes", "t4_lozenges")
  db <- as.data.frame(lapply(holzinger_swineford[, its], function(x) {
    cut(x, stats::quantile(x, 0:4 / 4), include.lowest = TRUE,
        labels = FALSE)
  }))
  names(db) <- c("y1", "y2", "y3")
  db$x <- holzinger_swineford$age_years
  fit <- lavaan::sem("f =~ y1 + y2 + y3\ny1 ~ x", data = db, std.lv = TRUE,
                     ordered = c("y1", "y2", "y3"))
  expect_true(lavaan::lavInspect(fit, "est")$gamma["y1", "x"] != 0)
  ebm <- lavaan::lavPredict(fit, method = "EBM")[, "f"]
  key <- apply(db[, c("y1", "y2", "y3")], 1L, paste, collapse = "")
  spread <- tapply(ebm, key, function(v) diff(range(v)))
  ages <- tapply(db$x, key, function(v) diff(range(v)))
  expect_lt(max(spread[ages > 0]), 1e-6)
  expect_error(factor_scores(fit, data = db),
               paste0("y1 is regressed on an exogenous covariate (x) and ",
                      "also depends on a factor directly."), fixed = TRUE)
  # A covariate that predicts the factor is not a direct effect.
  fit_f <- lavaan::sem("f =~ y1 + y2 + y3\nf ~ x", data = db, std.lv = TRUE,
                       ordered = c("y1", "y2", "y3"))
  expect_s3_class(factor_scores(fit_f, data = db), "dmar_factor_scores")
  # Not conditional on the covariates, lavPredict() stops on the fit.
  fit_u <- lavaan::sem("f =~ y1 + y2 + y3\nf ~ x", data = db, std.lv = TRUE,
                       ordered = c("y1", "y2", "y3"), conditional.x = FALSE)
  expect_error(lavaan::lavPredict(fit_u), "replacement has length zero")
  expect_error(factor_scores(fit_u, data = db),
               "that is not conditional on them (conditional.x = FALSE)",
               fixed = TRUE)
})

test_that("ordered items: an observed variable that covaries with a factor is refused", {
  skip_if_not_installed("lavaan")
  # lavaan's EBM score takes the factor's prior from its model implied mean,
  # which leaves out the value of a variable that covaries with the factor
  # (f ~~ x) without being an exogenous covariate: the scores are the same
  # for every x within a response pattern, up to the optimizer's tolerance.
  # The refusal named the variable an observed variable regressed on
  # others.
  its <- c("t1_visual_perception", "t2_cubes", "t4_lozenges")
  db <- as.data.frame(lapply(holzinger_swineford[, its], function(x) {
    cut(x, stats::quantile(x, 0:4 / 4), include.lowest = TRUE,
        labels = FALSE)
  }))
  names(db) <- c("y1", "y2", "y3")
  db$x <- holzinger_swineford$t5_general_information
  fit <- lavaan::sem("f =~ y1 + y2 + y3\nf ~~ x", data = db, std.lv = TRUE,
                     ordered = c("y1", "y2", "y3"))
  expect_gt(abs(lavaan::lavInspect(fit, "cor.all")["f", "x"]), 0.3)
  ebm <- lavaan::lavPredict(fit, method = "EBM")[, "f"]
  key <- apply(db[, c("y1", "y2", "y3")], 1L, paste, collapse = "")
  expect_lt(max(tapply(ebm, key, function(v) diff(range(v)))), 1e-3)
  expect_error(factor_scores(fit, data = db), paste0(
    "The observed variable x, regressed on nothing, covaries with the ",
    "factor f, and lavaan holds it among the latent variables rather than ",
    "as an exogenous covariate."), fixed = TRUE)
  # Written as a predictor of the factor, x is an exogenous covariate the
  # EBM score conditions on.
  fit_x <- lavaan::sem("f =~ y1 + y2 + y3\nf ~ x", data = db, std.lv = TRUE,
                       ordered = c("y1", "y2", "y3"))
  fx <- factor_scores(fit_x, data = db)
  expect_gt(max(tapply(fx$.derived_regression, key, function(v) diff(range(v)))), 0.1)
})

# The ML score of one case of ordered items, one factor, whose items are
# independent given the factor but for one pair, by direct maximization of
# the model's likelihood: `l` and `s` are each item's loading and residual
# standard deviation (of its latent response given the factor), `pair` the
# two items whose residuals correlate `rho`, `cuts` the thresholds of each
# item with -Inf and Inf, and `y` the case's categories. The probability of
# the pair's two categories is a one-dimensional integral of the bivariate
# normal density.
fs_ml_pair <- function(cuts, l, s, pair, rho, y) {
  r <- sqrt(1 - rho^2)
  nll <- function(f) {
    mu <- l * f
    ll <- 0
    for (j in setdiff(seq_along(l), pair)) {
      ll <- ll + log(stats::pnorm((cuts[[j]][y[j] + 1L] - mu[j]) / s[j]) -
                       stats::pnorm((cuts[[j]][y[j]] - mu[j]) / s[j]))
    }
    a <- (cuts[[pair[1L]]][y[pair[1L]] + 0:1] - mu[pair[1L]]) / s[pair[1L]]
    b <- (cuts[[pair[2L]]][y[pair[2L]] + 0:1] - mu[pair[2L]]) / s[pair[2L]]
    p <- stats::integrate(function(u) {
      stats::dnorm(u) * (stats::pnorm((b[2L] - rho * u) / r) -
                           stats::pnorm((b[1L] - rho * u) / r))
    }, a[1L], a[2L], rel.tol = 1e-10)$value
    -(ll + log(p))
  }
  stats::optimize(nll, c(-8, 8), tol = 1e-10)$minimum
}

test_that("ordered items: a regression on an observed variable is refused, a residual covariance warned of", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # EBM and ML scores optimize case by case
  # lavaan's EBM and ML scores take each case's likelihood as the product of
  # its items' probabilities given the factor, with the latent response
  # means of lav_lisrel_eyetax(), which takes an observed variable held
  # among the latent variables as zero where it predicts another. With
  # local dependence written y4 ~ y3, y4's latent response then has the
  # wrong mean and variance given the factor, and the fit was scored
  # without a word; with y3 ~~ y4 the pair's residual correlation is left
  # out. Against the maximizers of the model's own likelihood over the 160
  # response patterns with a finite score, lavaan's ML scores are off by up
  # to 0.18 in the first and 0.11 (15% of the standard error) in the
  # second.
  its <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
           "t4_lozenges")
  db <- as.data.frame(lapply(holzinger_swineford[, its], function(x) {
    cut(x, stats::quantile(x, 0:4 / 4), include.lowest = TRUE,
        labels = FALSE)
  }))
  names(db) <- paste0("y", 1:4)
  key <- apply(db, 1L, paste, collapse = "")
  for (spec in c("y4 ~ y3", "y3 ~~ y4")) {
    fit <- lavaan::sem(paste("f =~ y1 + y2 + y3 + y4\n", spec), data = db,
                       std.lv = TRUE, ordered = TRUE)
    e <- lavaan::lavInspect(fit, "est")
    th <- lavaan::lavInspect(fit, "theta")
    cuts <- lapply(names(db), function(j) {
      c(-Inf, e$tau[startsWith(rownames(e$tau), paste0(j, "|")), 1], Inf)
    })
    if (spec == "y4 ~ y3") {
      b <- e$beta["y4", "y3"]
      p3 <- e$psi["y3", "y3"]
      l <- c(e$lambda[c("y1", "y2"), "f"], e$beta["y3", "f"],
             e$beta["y4", "f"] + b * e$beta["y3", "f"])
      s <- c(sqrt(diag(th))[c("y1", "y2")], sqrt(p3),
             sqrt(b^2 * p3 + e$psi["y4", "y4"]))
      rho <- b * p3 / (s[3L] * s[4L])
    } else {
      l <- e$lambda[, "f"]
      s <- sqrt(diag(th))
      rho <- th["y3", "y4"] / (s[3L] * s[4L])
    }
    ml <- lavaan::lavPredict(fit, method = "ML", se = "standard")
    rows <- which(!duplicated(key) & is.finite(ml[, "f"]) &
                    abs(ml[, "f"]) < 5)
    expect_length(rows, 160L)
    model <- vapply(rows, function(i) {
      fs_ml_pair(cuts, l, s, 3:4, rho, unlist(db[i, ]))
    }, numeric(1))
    off <- abs(ml[rows, "f"] - model)
    # lavaan's scores maximize the product of the items' probabilities,
    # with y4's loading and residual variance read from its rows of B and
    # Psi alone.
    brute <- fs_ml_brute(fit, "f", rows[1:20])
    expect_lt(max(abs(ml[rows[1:20], "f"] - brute$score[, 1])), 2.5e-4)
    if (spec == "y4 ~ y3") {
      expect_gt(max(off), 0.18)
      expect_error(factor_scores(fit, data = db), paste0(
        "The regression y4 ~ y3 is on an observed variable that lavaan ",
        "holds among the latent variables rather than as an exogenous ",
        "covariate."), fixed = TRUE)
    } else {
      expect_gt(max(off), 0.11)
      expect_lt(max(off), 0.12)
      expect_gt(max(off / attr(ml, "se")[[1]][rows, "f"]), 0.15)
      expect_equal(stats::median(off), 0.03, tolerance = 0.05)
      expect_warning(fs <- factor_scores(fit, data = db),
                     "leaves out the residual covariance y3 ~~ y4: ",
                     fixed = TRUE)
      fin <- is.finite(fs$.derived_bartlett)
      expect_identical(sum(fin), 280L)
      expect_identical(fs$.derived_bartlett[fin], unname(ml[fin, "f"]))
    }
  }
  # An observed outcome of the factor (z ~ f) is not such a regression.
  db$z <- holzinger_swineford$t5_general_information
  fit_z <- lavaan::sem("f =~ y1 + y2 + y3 + y4\nz ~ f", data = db,
                       std.lv = TRUE, ordered = paste0("y", 1:4))
  expect_no_warning(fz <- factor_scores(fit_z, data = db))
  expect_s3_class(fz, "dmar_factor_scores")
})

test_that("ordered items: a term lavaan leaves out that is constant in the factors does not refuse the fit", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # EBM and ML scores optimize case by case
  # lavaan's EBM and ML scores take an observed variable held among the
  # latent variables as zero where it predicts another, and leave a
  # covariate's direct effect out. In a fully mediated chain of continuous
  # variables (m ~ f, w ~ m) and for a covariate's effect on an outcome
  # unrelated to the factor (w ~ x), the term left out is constant in the
  # factor, so the scores are the model's; such fits were refused.
  its <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
           "t4_lozenges")
  db <- as.data.frame(lapply(holzinger_swineford[, its], function(x) {
    cut(x, stats::quantile(x, 0:4 / 4), include.lowest = TRUE,
        labels = FALSE)
  }))
  ord <- paste0("y", 1:4)
  names(db) <- ord
  db$m <- holzinger_swineford$t6_paragraph_comprehension
  db$w <- holzinger_swineford$t7_sentence
  db$x <- holzinger_swineford$age_years
  fit_with <- function(extra, data = db, ...) {
    lavaan::sem(paste("f =~ y1 + y2 + y3 + y4\n", extra), data = data,
                std.lv = TRUE, ordered = ord, ...)
  }
  # The model's ML and EBM scores by direct maximization of each case's
  # likelihood, the items' ordered probit probabilities times `extra`
  # (the density of the observed variables that depend on the factor), and
  # for EBM the factor's standard normal prior.
  direct <- function(fit, rows, extra) {
    e <- lavaan::lavInspect(fit, "est")
    lam <- e$lambda[ord, "f"]
    s <- sqrt(diag(lavaan::lavInspect(fit, "theta"))[ord])
    cuts <- lapply(ord, function(j) {
      c(-Inf, e$tau[startsWith(rownames(e$tau), paste0(j, "|")), 1], Inf)
    })
    t(vapply(rows, function(i) {
      y <- unlist(db[i, ord])
      nll <- function(f) {
        -sum(log(stats::pnorm((mapply(`[`, cuts, y + 1L) - lam * f) / s) -
                   stats::pnorm((mapply(`[`, cuts, y) - lam * f) / s))) -
          extra(f, i, e)
      }
      c(ml = stats::optimize(nll, c(-8, 8), tol = 1e-10)$minimum,
        ebm = stats::optimize(function(f) {
          nll(f) - stats::dnorm(f, log = TRUE)
        }, c(-8, 8), tol = 1e-10)$minimum)
    }, numeric(2)))
  }
  rows <- 1:60
  fit_a <- fit_with("m ~ f\nw ~ m")
  fa <- factor_scores(fit_a, data = db)
  expect_true(all(fa$reason[rows] ==
                    "determinacy is not computed for ordered items"))
  o <- direct(fit_a, rows, function(f, i, e) {
    mu <- e$nu["m", 1] + e$alpha["m", 1] + e$beta["m", "f"] * f
    stats::dnorm(db$m[i], mu, sqrt(e$psi["m", "m"]), log = TRUE)
  })
  expect_lt(max(abs(fa$.derived_bartlett[rows] - o[, "ml"])), 2e-3)
  expect_lt(max(abs(fa$.derived_regression[rows] - o[, "ebm"])), 1e-3)
  # m's density enters: the scores are not those of the items alone.
  o_items <- direct(fit_a, rows, function(f, i, e) 0)
  expect_gt(max(abs(fa$.derived_regression[rows] - o_items[, "ebm"])), 0.3)
  fit_x <- fit_with("w ~ x")
  expect_true(lavaan::lavInspect(fit_x, "est")$gamma["w", "x"] != 0)
  fx <- factor_scores(fit_x, data = db)
  ox <- direct(fit_x, rows, function(f, i, e) 0)
  # A row whose items all sit at one end has no finite ML score.
  fin <- is.finite(fx$.derived_bartlett[rows])
  expect_true(all(grepl("every observed item of this factor is at an end",
                        fx$reason[rows][!fin])))
  expect_lt(max(abs(fx$.derived_bartlett[rows][fin] - ox[fin, "ml"])), 2e-3)
  expect_lt(max(abs(fx$.derived_regression[rows] - ox[, "ebm"])), 1e-3)
  # The refusals concern a variable that depends on a factor directly: w,
  # which depends on f only through the observed m, may have a direct
  # effect of a covariate, and the scores are those of the items and m.
  # (Age is centered: the conditional estimator fits m's intercept given
  # the covariate, which with m not regressed on it is m's mean only for a
  # centered covariate.)
  dx <- db
  dx$x <- dx$x - mean(dx$x)
  fit_mx <- fit_with("m ~ f\nw ~ m + x", data = dx)
  expect_true(lavaan::lavInspect(fit_mx, "est")$gamma["w", "x"] != 0)
  fmx <- factor_scores(fit_mx, data = dx)
  omx <- direct(fit_mx, rows, function(f, i, e) {
    mu <- e$nu["m", 1] + e$alpha["m", 1] + e$beta["m", "f"] * f
    stats::dnorm(db$m[i], mu, sqrt(e$psi["m", "m"]), log = TRUE)
  })
  expect_lt(max(abs(fmx$.derived_bartlett[rows] - omx[, "ml"])), 2e-3)
  expect_lt(max(abs(fmx$.derived_regression[rows] - omx[, "ebm"])), 1e-3)
  # Where the term left out involves the factor, the fit is still refused:
  # the variable predicted also depends on the factor, the predictor is an
  # ordered item, or a row has the predictor missing and the variable it
  # predicts observed (pairwise deletion), and an outcome of a covariate
  # whose disturbance covaries with the factor (lavaan frees f ~~ w by
  # default when x also predicts f).
  expect_error(factor_scores(fit_with("m ~ f\nw ~ m + f"), data = db),
               "here the variable predicted depends on a factor, so",
               fixed = TRUE)
  expect_error(factor_scores(fit_with("w ~ y3"), data = db), paste0(
    "only when the variable predicted is an observed variable that depends ",
    "on no factor directly and the predictor is a continuous variable ",
    "(where a case lacks it, one that reaches no factor except through a ",
    "variable the case has observed); here the predictor is an ordered ",
    "item, so"), fixed = TRUE)
  dp <- db
  set.seed(3)
  dp$m[sample(nrow(dp), 20)] <- NA
  expect_error(factor_scores(fit_with("m ~ f\nw ~ m", data = dp,
                                      missing = "pairwise"), data = dp),
               paste0("and m is missing in 20 rows the fit keeps where w is ",
                      "observed and where m depends on a factor through no ",
                      "variable the row has observed"),
               fixed = TRUE)
  # A variable the row has observed on the predictor's path to the factor
  # cuts the dependence: with m ~ f, w ~ m, and s ~ w, and w missing where
  # m and s are observed, the likelihood of s given m does not involve f,
  # so the EBM and ML scores of those rows are the model's (those of the
  # items and m). Such a fit was refused, with the claim that the
  # likelihood of s depends on the factor through w. With m missing as
  # well, the path from w to f runs through no observed variable, and the
  # fit is refused.
  ds <- db
  ds$s <- holzinger_swineford$t9_word_meaning
  set.seed(3)
  gap <- sample(nrow(ds), 20)
  ds$w[gap] <- NA
  fit_s <- fit_with("m ~ f\nw ~ m\ns ~ w", data = ds, missing = "pairwise")
  fs_s <- factor_scores(fit_s, data = ds)
  os <- direct(fit_s, gap, function(f, i, e) {
    mu <- e$nu["m", 1] + e$alpha["m", 1] + e$beta["m", "f"] * f
    stats::dnorm(db$m[i], mu, sqrt(e$psi["m", "m"]), log = TRUE)
  })
  fin <- is.finite(fs_s$.derived_bartlett[gap])
  expect_gt(sum(fin), 10L)
  expect_lt(max(abs(fs_s$.derived_bartlett[gap][fin] - os[fin, "ml"])), 2e-3)
  expect_lt(max(abs(fs_s$.derived_regression[gap] - os[, "ebm"])), 1e-3)
  ds$m[gap[1:5]] <- NA
  expect_error(factor_scores(fit_with("m ~ f\nw ~ m\ns ~ w", data = ds,
                                      missing = "pairwise"), data = ds),
               paste0("and w is missing in 5 rows the fit keeps where s is ",
                      "observed and where w depends on a factor through no ",
                      "variable the row has observed"),
               fixed = TRUE)
  # A predictor that depends on no factor may be missing.
  expect_s3_class(factor_scores(fit_with("m ~ x\nw ~ m", data = dp,
                                         missing = "pairwise"), data = dp),
                  "dmar_factor_scores")
  expect_error(factor_scores(fit_with("f ~ x\nw ~ x"), data = db), paste0(
    "The disturbance of the observed variable w, regressed on exogenous ",
    "covariates only, covaries with the factor f."), fixed = TRUE)
  expect_s3_class(factor_scores(fit_with("f ~ x\nw ~ x\nf ~~ 0*w"),
                                data = db), "dmar_factor_scores")
})

test_that("ordered items: another factor's extreme pattern does not cost a finite ML score", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # EBM and ML scores optimize case by case
  # lavaan optimizes the ML scores of all factors of a case jointly, and a
  # factor whose items all sit at one end can stop its optimizer on the
  # iteration limit, which leaves every factor of the case without a
  # score. The textual score, which the likelihood separates from the
  # visual one under simple structure, was reported as "score unavailable
  # from lavaan".
  its <- c(fs_items$visual, fs_items$textual)
  db <- as.data.frame(lapply(holzinger_swineford[1:150, its], function(x) {
    as.integer(x > stats::median(x))
  }))
  model <- paste("visual =~", paste(fs_items$visual, collapse = " + "),
                 "\ntextual =~", paste(fs_items$textual, collapse = " + "))
  fit <- lavaan::cfa(model, data = db, std.lv = TRUE, ordered = TRUE)
  fs <- factor_scores(fit, data = db)
  ml <- lavaan::lavPredict(fit, method = "ML")
  tx <- fs$factor == "textual"
  lost <- which(is.na(ml[, "textual"]) &
                  !grepl("end of its scale", fs$reason[tx]))
  expect_gt(length(lost), 10L)
  expect_true(all(is.finite(fs$.derived_bartlett[tx][lost])))
  expect_true(all(is.finite(fs$.derived_bartlett_se[tx][lost])))
  expect_false(any(grepl("unavailable", fs$reason)))
  # The recovered score is that of the textual items alone: the score of
  # every other case with the same textual responses.
  key <- apply(db[, fs_items$textual], 1L, paste, collapse = "")
  for (i in lost) {
    same <- setdiff(which(key == key[i]), lost)
    expect_equal(fs$.derived_bartlett[tx][i], stats::median(fs$.derived_bartlett[tx][same]),
                 tolerance = 1e-3)
  }
  # It comes from lavaan's nlminb optimizer.
  nl <- lavaan::lavPredict(fit, method = "ML", optim_method = "nlminb")
  expect_equal(fs$.derived_bartlett[tx][lost], unname(nl[lost, "textual"]))
})

test_that("ordered items: a higher order factor has no ML score", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # EBM and ML scores optimize case by case
  # lavaan's ML criterion depends on a factor without items of its own
  # only through the small ridge its optimizer adds, so its ML score of g
  # was 0 for every case, with a standard error of 316, an interval of
  # [-620, 620], and no reason.
  db <- as.data.frame(lapply(holzinger_swineford[1:150, unlist(fs_items)],
                             function(x) {
    cut(x, stats::quantile(x, c(0, 1 / 3, 2 / 3, 1)), include.lowest = TRUE,
        labels = FALSE)
  }))
  fit <- lavaan::cfa(paste(fs_model, "\ng =~ visual + textual + speed"),
                     data = db, ordered = TRUE)
  expect_true(lavaan::lavInspect(fit, "post.check"))
  ml <- lavaan::lavPredict(fit, method = "ML")
  expect_true(all(abs(ml[, "g"]) < 1e-8, na.rm = TRUE))
  fs <- factor_scores(fit, data = db)
  g <- fs$factor == "g"
  expect_true(all(is.na(fs$.derived_bartlett[g])))
  expect_true(all(is.na(fs$.derived_bartlett_se[g])))
  expect_true(all(grepl(paste0("no ML score for a factor without items of ",
                               "its own"), fs$reason[g])))
  expect_false(any(grepl("end of its scale", fs$reason[g])))
  # Its EBM score, which the model informs, is kept, and is lavaan's.
  ebm <- lavaan::lavPredict(fit, method = "EBM")
  expect_equal(fs$.derived_regression[g], unname(ebm[, "g"]))
})

test_that("ordered items: the reasons name what lavaan can and cannot do", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # EBM and ML scores optimize case by case
  its <- fs_items$visual
  db <- as.data.frame(lapply(holzinger_swineford[, its], function(x) {
    cut(x, stats::quantile(x, 0:4 / 4), include.lowest = TRUE,
        labels = FALSE)
  }))
  db$t2_cubes[1:10] <- NA
  model <- paste("v =~", paste(its, collapse = " + "))
  # Listwise deletion: the remedy is pairwise deletion, since lavaan has no
  # FIML for ordered items. The reason used to recommend missing = "ml".
  fit <- lavaan::cfa(model, data = db, std.lv = TRUE, ordered = TRUE)
  expect_error(lavaan::cfa(model, data = db, std.lv = TRUE, ordered = TRUE,
                           missing = "ml"), "categorical")
  fs <- factor_scores(fit, data = db)
  why <- paste0("deleted listwise by the fit (with ordered items, a fit ",
                "with missing = \"pairwise\" scores it)")
  expect_true(all(fs$reason[fs$row <= 10] == why))
  expect_identical(unique(factor_scores(fit)$reason[1:10]), why)
  # The summary counts it as one reason.
  out <- capture.output(print(fs))
  expect_true(any(grepl(why, out, fixed = TRUE)))
  fit_pw <- lavaan::cfa(model, data = db, std.lv = TRUE, ordered = TRUE,
                        missing = "pairwise")
  expect_true(all(is.finite(factor_scores(fit_pw,
                                          data = db)$.derived_regression[1:10])))
  # A case that answered none of the items, only a covariate, in a fit
  # conditional on it: lavaan skips the case, and the reason says there is
  # no EBM score rather than that one was borrowed.
  db$age <- holzinger_swineford$age_years
  db[5, its] <- NA
  # lavaan warns that the case is empty.
  fit_x <- suppressWarnings(lavaan::sem(paste(model, "\nv ~ age"), data = db,
                                        std.lv = TRUE, ordered = TRUE,
                                        missing = "pairwise"))
  fs_x <- factor_scores(fit_x, data = db)
  expect_true(is.na(fs_x$.derived_regression[5]))
  expect_identical(fs_x$reason[5], paste0(
    "no item of this factor observed: no ML score, and no EBM score from ",
    "lavaan; determinacy is not computed for ordered items"))
})

test_that("ordered items: a solution with a negative residual variance is not scored", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # a WLSMV fit of a bifactor model
  # In this orthogonal bifactor model t4_lozenges' residual variance is far
  # below zero, and lavaan computes no EBM or ML score. The warning said
  # the scores inherit the improper solution and every reason read "score
  # unavailable from lavaan".
  dq <- as.data.frame(lapply(holzinger_swineford[, unlist(fs_items)],
                             function(x) {
    cut(x, stats::quantile(x, 0:4 / 4), include.lowest = TRUE,
        labels = FALSE)
  }))
  model <- paste("g =~", paste(unlist(fs_items), collapse = " + "), "\n",
                 fs_model)
  # lavaan warns about the fit.
  fit <- suppressWarnings(lavaan::cfa(model, data = dq, ordered = TRUE,
                                      std.lv = TRUE, orthogonal = TRUE))
  expect_lt(lavaan::lavInspect(fit, "est")$theta["t4_lozenges",
                                                  "t4_lozenges"], 0)
  # DMAR's warning, and lavaan's own (twice) that it computed no scores.
  seen <- new.env(parent = emptyenv())
  seen$warnings <- character(0)
  fs <- withCallingHandlers(factor_scores(fit, data = dq),
                            warning = function(w) {
    seen$warnings <- c(seen$warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_true(any(grepl(paste0("negative residual variance for t4_lozenges)",
                               ". lavaan computes no EBM or ML score"),
                        seen$warnings, fixed = TRUE)))
  expect_true(all(is.na(fs$.derived_regression)))
  expect_true(all(grepl(paste0("no EBM or ML score: lavaan does not score a ",
                               "solution with a negative residual variance"),
                        fs$reason, fixed = TRUE)))
})

test_that("a covariate that predicts a factor is left to lavPredict()", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[, c(unlist(fs_items), "age_years")]
  d$age_years[c(10, 20)] <- NA
  # lavaan warns that the two rows with a missing covariate are deleted.
  fit <- suppressWarnings(lavaan::cfa(paste(fs_model, "\nvisual ~ age_years"),
                                      data = d, std.lv = TRUE,
                                      missing = "ml"))
  fs <- factor_scores(fit, data = d)
  expect_identical(attr(fs, "scored_from"), "fit")
  idx <- lavaan::lavInspect(fit, "case.idx")
  lp <- lavaan::lavPredict(fit)
  expect_identical(fs$.derived_regression[fs$factor == "visual"][idx],
                   unname(lp[, "visual"]))
  det <- attr(fs, "determinacy")
  expect_equal(det$determinacy_regression,
               unname(lavaan::lavInspect(fit, "fs.determinacy"))[1:3],
               tolerance = 1e-10)
  expect_equal(det$determinacy_bartlett,
               unname(lavaan::lavInspect(fit, "fs.determinacy.bartlett"))[1:3],
               tolerance = 1e-10)
  gone <- fs$row %in% c(10, 20)
  expect_true(all(is.na(fs$.derived_regression[gone])))
  expect_true(all(grepl("exogenous covariate", fs$reason[gone])))
})

test_that("a row deleted for a missing covariate says what keeps it, whatever the fit", {
  skip_if_not_installed("lavaan")
  # The reason for a row missing an exogenous covariate was given only under
  # FIML. Every other fit told it that missing = "ml" scores it, which
  # deletes it too (fixed.x = TRUE), and an ordered fit, pairwise ones
  # included, that missing = "pairwise" does; without data, a listwise fit
  # still told it so.
  d <- holzinger_swineford[, c("t1_visual_perception", "t2_cubes",
                               "t4_lozenges", "age_years")]
  names(d) <- c("y1", "y2", "y3", "x")
  d$x[1:5] <- NA
  d$y1[6:8] <- NA
  model <- "f =~ y1 + y2 + y3\nf ~ x"
  why_x <- paste0("deleted by the fit for a missing exogenous covariate ",
                  "(fixed.x = TRUE), which missing = \"ml.x\" scores")
  for (miss in c("listwise", "pairwise", "available.cases", "two.stage",
                 "ml")) {
    # lavaan warns that the rows with a missing covariate are deleted.
    fit <- suppressWarnings(lavaan::sem(model, data = d, std.lv = TRUE,
                                        missing = miss))
    expect_identical(lavaan::lavInspect(fit, "nobs"),
                     if (miss == "listwise") 293L else 296L)
    fs <- factor_scores(fit, data = d)
    expect_true(all(fs$reason[fs$row <= 5] == why_x))
    # Without data the reason cannot tell a missing covariate from a
    # missing item, and it named only what scores the second ("ml").
    fs_nd <- factor_scores(fit)
    expect_true(all(fs_nd$reason[fs_nd$row <= 5] == if (miss == "listwise") {
      paste0("deleted listwise by the fit (a fit by FIML, missing = ",
             "\"ml\", scores a row with a missing item, missing = ",
             "\"ml.x\" one with a missing exogenous covariate, and ",
             "supplying data identifies the reason)")
    } else {
      paste0("not used by the fit, which deletes a row with a missing ",
             "exogenous covariate (missing = \"ml.x\" scores such a row, ",
             "and supplying data identifies the reason)")
    }))
  }
  # missing = "ml.x" keeps and scores them.
  fit_mx <- lavaan::sem(model, data = d, std.lv = TRUE, missing = "ml.x")
  expect_true(all(is.finite(factor_scores(fit_mx, data = d)$.derived_regression)))
  # Conditional on the covariate, missing = "ml.x" needs conditional.x =
  # FALSE as well.
  fit_cx <- suppressWarnings(lavaan::sem(model, data = d, std.lv = TRUE,
                                         missing = "ml",
                                         conditional.x = TRUE))
  expect_error(lavaan::sem(model, data = d, std.lv = TRUE, missing = "ml.x",
                           conditional.x = TRUE), "ml.x")
  expect_true(all(factor_scores(fit_cx, data = d)$reason[1:5] == paste0(
    "deleted by the fit for a missing exogenous covariate (conditional.x ",
    "= TRUE), which missing = \"ml.x\" with conditional.x = FALSE scores")))
  expect_true(all(factor_scores(fit_cx)$reason[1:5] == paste0(
    "not used by the fit, which deletes a row with a missing exogenous ",
    "covariate (missing = \"ml.x\" with conditional.x = FALSE scores such ",
    "a row, and supplying data identifies the reason)")))
  # Ordered items: listwise and pairwise fits delete the row alike.
  db <- as.data.frame(lapply(d[, 1:3], function(v) {
    cut(v, stats::quantile(v, 0:4 / 4, na.rm = TRUE), include.lowest = TRUE,
        labels = FALSE)
  }))
  db$x <- d$x
  why_o <- paste0("deleted by the fit for a missing exogenous covariate ",
                  "(with ordered items, no fit that lavaan can score keeps ",
                  "it)")
  for (miss in c("listwise", "pairwise")) {
    fit_o <- suppressWarnings(lavaan::sem(model, data = db, std.lv = TRUE,
                                          ordered = c("y1", "y2", "y3"),
                                          missing = miss))
    fs_o <- factor_scores(fit_o, data = db)
    expect_true(all(fs_o$reason[1:5] == why_o))
    fs_ond <- factor_scores(fit_o)
    expect_true(all(fs_ond$reason[1:5] == if (miss == "listwise") {
      paste0("deleted listwise by the fit (with ordered items, a fit with ",
             "missing = \"pairwise\" scores a row with a missing item, no ",
             "fit that lavaan can score keeps one with a missing exogenous ",
             "covariate, and supplying data identifies the reason)")
    } else {
      paste0("not used by the fit, which deletes a row with a missing ",
             "exogenous covariate (with ordered items no fit that lavaan can ",
             "score keeps such a row, and supplying data identifies the ",
             "reason)")
    }))
    # A row missing an item only is deleted listwise, or scored under
    # pairwise deletion.
    if (miss == "listwise") {
      expect_true(all(grepl("deleted listwise", fs_o$reason[6:8])))
    } else {
      expect_true(all(is.finite(fs_o$.derived_regression[6:8])))
    }
  }
})

test_that("a higher order factor's items are those of its factors", {
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa(paste(fs_model, "\ng =~ visual + textual + speed"),
                     data = holzinger_swineford, std.lv = TRUE)
  fs <- factor_scores(fit)
  expect_identical(unique(fs$items[fs$factor == "g"]), 9L)
  lp <- lavaan::lavPredict(fit, method = "Bartlett")
  expect_identical(fs$.derived_bartlett[fs$factor == "g"], unname(lp[, "g"]))
  det <- attr(fs, "determinacy")
  expect_equal(det$determinacy_regression,
               unname(lavaan::lavInspect(fit, "fs.determinacy")),
               tolerance = 1e-10)
  expect_true(is.na(det$determinacy_bartlett[det$factor == "g"]))
  expect_true(all(grepl("without items of its own",
                        fs$reason[fs$factor == "g"])))
  expect_true(all(fs$reason[fs$factor != "g"] == ""))
})

test_that("the pattern table carries lavaan's Bartlett standard error of a higher order factor", {
  skip_if_not_installed("lavaan")
  # The rows carried lavaan's Bartlett score and standard error of g, which
  # lavaan forms through the loadings of the factors g measures, and the
  # pattern table reported NA for that standard error.
  d <- holzinger_swineford[, unlist(fs_items)]
  d$t2_cubes[1:20] <- NA
  fit <- lavaan::cfa(paste(fs_model, "\ng =~ visual + textual + speed"),
                     data = d, std.lv = TRUE, missing = "ml")
  fs <- factor_scores(fit, data = d)
  pt <- attr(fs, "patterns")
  m <- merge(fs[, c("pattern", "factor", 
                  ".derived_bartlett_se")],
             pt[, c("pattern", "factor", "bartlett_se")],
             by = c("pattern", "factor"))
  expect_equal(m$.derived_bartlett_se, m$bartlett_se)
  g <- pt$factor == "g"
  expect_true(all(is.finite(pt$bartlett_se[g])))
  expect_true(all(is.na(pt$determinacy_bartlett[g])))
  lp <- lavaan::lavPredict(fit, method = "Bartlett", se = "standard")
  expect_setequal(round(pt$bartlett_se[g], 12),
                  round(unique(attr(lp, "se")[[1]][, "g"]), 12))
})

test_that("a case without any of a factor's items is scored from its other observed variables", {
  skip_if_not_installed("lavaan")
  # The reason said the regression score "comes only from correlated
  # factors"; with one factor and a covariate it comes from the covariate.
  its <- fs_items$visual
  d <- holzinger_swineford[, c(its, "age_years")]
  d[5, its] <- NA
  fit <- lavaan::sem(paste("v =~", paste(its, collapse = " + "),
                           "\nv ~ age_years"),
                     data = d, std.lv = TRUE, missing = "ml")
  fs <- factor_scores(fit, data = d)
  b <- lavaan::coef(fit)[["v~age_years"]]
  expect_equal(fs$.derived_regression[5], b * d$age_years[5], tolerance = 1e-8)
  expect_identical(fs$reason[5], paste0(
    "no item of this factor observed: no Bartlett score, and the ",
    "regression score comes only from the case's other observed variables"))
  # Conditional on the covariate, lavaan sets the case aside, and the
  # reason says there is no regression score.
  # lavaan warns that the case is empty.
  fit_x <- suppressWarnings(lavaan::sem(
    paste("v =~", paste(its, collapse = " + "), "\nv ~ age_years"),
    data = d, std.lv = TRUE, missing = "ml", conditional.x = TRUE))
  fs_x <- factor_scores(fit_x, data = d)
  expect_true(is.na(fs_x$.derived_regression[5]))
  expect_true(is.na(fs_x$.derived_regression_se[5]))
  expect_match(fs_x$reason[5], paste0(
    "^no item of this factor observed: no Bartlett score, and no ",
    "regression score from lavaan"))
})

test_that("a higher order factor scored from the indicators says why it has no Bartlett score", {
  skip_if_not_installed("lavaan")
  # With an auxiliary variable carried as a saturated correlate the scores
  # are computed here from the indicators, and a higher order factor gets
  # no Bartlett score. Its reason used to read "score unavailable from
  # lavaan", although lavaan was not asked for it.
  d <- fs_mar_data("age_years")
  d$age_years <- as.numeric(scale(d$age_years))
  aux <- paste("age_years ~~", paste(unlist(fs_items), collapse = " + "))
  # lavaan warns that the two rows with every item missing are ignored.
  fit <- suppressWarnings(lavaan::cfa(
    paste(fs_model, "\ng =~ visual + textual + speed\n", aux), data = d,
    std.lv = TRUE, missing = "ml"))
  expect_message(fs <- factor_scores(fit, data = d), "saturated correlates")
  expect_identical(attr(fs, "scored_from"), "indicators")
  g <- fs$factor == "g"
  empty <- fs$row %in% c(5, 77)
  expect_true(all(is.na(fs$.derived_bartlett[g])))
  expect_true(all(is.finite(fs$.derived_regression[g & !empty])))
  expect_identical(unique(fs$reason[g & !empty]), paste0(
    "no Bartlett score for a factor without items of its own when the ",
    "scores are computed from the indicators (auxiliary variables set ",
    "aside)"))
  expect_false(any(grepl("lavaan", fs$reason)))
  expect_true(all(fs$reason[!g & !empty] == ""))
  # The regression score of g is lavaan's for the same model without the
  # auxiliary, every parameter fixed at the auxiliary fit's estimate.
  pt <- lavaan::parTable(fit)
  pt <- pt[!(pt$lhs == "age_years" | pt$rhs == "age_years"), ]
  pt$free <- 0L
  pt$ustart <- pt$est
  pt$start <- pt$est
  pt$se <- NULL
  fixed <- suppressWarnings(lavaan::lavaan(pt, data = d[, unlist(fs_items)],
                                           missing = "ml"))
  idx <- lavaan::lavInspect(fixed, "case.idx")
  lp <- lavaan::lavPredict(fixed)
  expect_equal(fs$.derived_regression[g][idx], unname(lp[, "g"]), tolerance = 1e-8)
})

test_that("a higher order factor keeps its Bartlett score when lavaan counts an outcome's mean twice", {
  skip_if_not_installed("lavaan")
  # Without a mean structure lavaan counts an outcome's mean twice, and the
  # scores of such a fit are computed here (215fce2); a higher order factor
  # then had no Bartlett score, which lavaan forms through the loadings of
  # the factors it measures. lavaan's score is off by a constant only (by
  # -1.09 when the spatial factor predicts the paper form board test, and
  # -0.36 when the verbal factor predicts general information), and less
  # that constant it is the score of the fit with a mean structure.
  hs <- holzinger_swineford
  d <- data.frame(y1 = hs$t1_visual_perception, y2 = hs$t2_cubes,
                  y3 = hs$t4_lozenges, y4 = hs$t6_paragraph_comprehension,
                  y5 = hs$t7_sentence, y6 = hs$t9_word_meaning,
                  y7 = hs$t10_addition / 10,
                  y8 = hs$t12_counting_groups_of_dots / 5,
                  y9 = hs$t13_straight_and_curved_capitals / 10,
                  z = hs$t5_general_information / 3,
                  w = hs$t3_paper_form_board)
  base <- paste("f1 =~ y1 + y2 + y3\nf2 =~ y4 + y5 + y6",
                "f3 =~ y7 + y8 + y9\ng =~ f1 + f2 + f3", sep = "\n")
  for (out in c("w ~ f1", "z ~ f2")) {
    model <- paste(base, out, sep = "\n")
    fit <- lavaan::sem(model, data = d, std.lv = TRUE)
    fit_ms <- lavaan::sem(model, data = d, std.lv = TRUE,
                          meanstructure = TRUE)
    shift <- lavaan::lavPredict(fit, method = "Bartlett")[, "g"] -
      lavaan::lavPredict(fit_ms, method = "Bartlett")[, "g"]
    expect_lt(diff(range(shift)), 1e-8)
    expect_equal(mean(shift), if (out == "w ~ f1") -1.086 else -0.361,
                 tolerance = 1e-3)
    fs <- factor_scores(fit, data = d)
    fs_ms <- factor_scores(fit_ms, data = d)
    for (col in c(".derived_regression", ".derived_regression_se",
                ".derived_bartlett", ".derived_bartlett_se",
                  "determinacy_regression", "determinacy_bartlett",
                  "reason")) {
      expect_equal(fs[[col]], fs_ms[[col]], tolerance = 1e-8)
    }
    g <- fs$factor == "g"
    expect_true(all(is.finite(fs$.derived_bartlett[g])))
    expect_identical(unique(fs$reason[g]), paste0(
      "Bartlett determinacy is not computed for a factor without items of ",
      "its own"))
    pt <- attr(fs, "patterns")
    expect_equal(pt$bartlett_se, attr(fs_ms, "patterns")$bartlett_se,
                 tolerance = 1e-8)
  }
  # Where lavaan gives an observed variable regressed on nothing a row of
  # zeros among the loadings (g ~~ w), its higher order Bartlett score is
  # not of this model, and the reason names the variable; the regression
  # score is the factor's conditional mean.
  model_w <- paste(base, "g ~~ w", sep = "\n")
  fit_w <- lavaan::sem(model_w, data = d, std.lv = TRUE,
                       meanstructure = TRUE)
  fw <- factor_scores(fit_w, data = d)
  g <- fw$factor == "g"
  expect_true(all(is.na(fw$.derived_bartlett[g])))
  expect_identical(unique(fw$reason[g]), paste0(
    "no Bartlett score for a factor without items of its own when the ",
    "scores are computed from the fit's estimates rather than by lavaan, ",
    "whose loadings give w, regressed on nothing, a row of zeros"))
  o <- fs_cond_oracle(fit_w, d, "g", character(0))
  expect_equal(fw$.derived_regression[g], o$regression, tolerance = 1e-10)
})

test_that("ordered items require lavaan 0.7-3, whose EBM and ML scores count the intercepts once", {
  skip_if_not_installed("lavaan")
  expect_error(.fs_require_lavaan_ordered(lavaan_version = "0.7-2"),
               "require lavaan 0.7-3 or newer (installed: 0.7-2)",
               fixed = TRUE)
  expect_error(.fs_require_lavaan_ordered(lavaan_version = "0.7-2"),
               "counted the intercepts of the ordered indicators twice",
               fixed = TRUE)
  expect_true(.fs_require_lavaan_ordered(lavaan_version = "0.7-3"))
  expect_true(.fs_require_lavaan_ordered(lavaan_version = "0.10-1"))
  # factor_scores() asks for the floor exactly when the fit has ordered
  # items; the check stands in here for an installed lavaan 0.7-2.
  real <- .fs_require_lavaan_ordered
  local_mocked_bindings(.fs_require_lavaan_ordered = function(...) {
    real(lavaan_version = "0.7-2")
  })
  db <- as.data.frame(lapply(holzinger_swineford[, fs_items$visual],
                             function(x) {
    cut(x, stats::quantile(x, c(0, .25, .5, .75, 1)),
        include.lowest = TRUE, labels = FALSE)
  }))
  fit_o <- lavaan::cfa("visual =~ t1_visual_perception + t2_cubes +
                        t4_lozenges", data = db, std.lv = TRUE,
                       ordered = TRUE)
  expect_error(factor_scores(fit_o), "require lavaan 0.7-3", fixed = TRUE)
  fit_c <- lavaan::cfa(fs_model, data = holzinger_swineford, std.lv = TRUE)
  expect_no_error(factor_scores(fit_c))
})

test_that("a fit conditional on exogenous covariates is scored by lavaan", {
  skip_if_not_installed("lavaan")
  fit <- lavaan::sem(paste(fs_model, "\nvisual ~ age_years"),
                     data = holzinger_swineford, std.lv = TRUE,
                     conditional.x = TRUE)
  fs <- factor_scores(fit)
  expect_identical(fs$.derived_regression[fs$factor == "visual"],
                   unname(lavaan::lavPredict(fit)[, "visual"]))
  expect_true(all(is.na(fs$determinacy_regression)))
  expect_true(all(grepl("conditional on exogenous covariates", fs$reason)))
})

test_that("an uncorrelated factor without items and an error free indicator", {
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[, c(unlist(fs_items), "age_years")]
  d[c(10, 11), fs_items$textual] <- NA
  fit <- lavaan::cfa(fs_model, data = d, std.lv = TRUE, missing = "ml",
                     orthogonal = TRUE)
  fs <- factor_scores(fit, data = d)
  # With no correlated factor to borrow from, the regression score is the
  # factor mean, its standard error the factor's standard deviation, and
  # its determinacy undefined.
  tx <- fs$row %in% c(10, 11) & fs$factor == "textual"
  expect_equal(fs$.derived_regression[tx], c(0, 0))
  expect_equal(fs$.derived_regression_se[tx], c(1, 1))
  expect_true(all(is.na(fs$determinacy_regression[tx])))
  expect_true(all(grepl("no item of this factor observed", fs$reason[tx])))
  # A factor measured by one indicator without error is known exactly:
  # lavaan reports no standard error, and the reason says why.
  # lavaan warns about the scale of age_years against the test scores.
  fit1 <- suppressWarnings(lavaan::cfa(
    paste(fs_model, "\nage =~ 1*age_years\nage_years ~~ 0*age_years"),
    data = d, missing = "ml"))
  fs1 <- factor_scores(fit1, data = d)
  age <- fs1$factor == "age"
  expect_equal(fs1$.derived_regression[age], d$age_years - mean(d$age_years),
               tolerance = 1e-6)
  expect_true(all(is.na(fs1$.derived_regression_se[age])))
  expect_true(all(grepl("no regression standard error", fs1$reason[age])))
  expect_equal(unique(fs1$determinacy_regression[age]), 1)
})

test_that("what cannot be scored is refused, and an improper solution warned", {
  skip_if_not_installed("lavaan")
  expect_error(factor_scores(stats::lm(t1_visual_perception ~ t2_cubes,
                                       data = holzinger_swineford)),
               "lavaan fit object")
  s <- stats::cov(holzinger_swineford[, unlist(fs_items)])
  fit_s <- lavaan::cfa(fs_model, sample.cov = s, sample.nobs = 301)
  expect_error(factor_scores(fit_s), "summary statistics")
  fit_nc <- suppressWarnings(lavaan::cfa(fs_model, data = holzinger_swineford,
                                         control = list(iter.max = 2)))
  expect_error(factor_scores(fit_nc), "did not converge")
  fit <- lavaan::cfa(fs_model, data = holzinger_swineford, std.lv = TRUE)
  expect_error(factor_scores(fit, conf_level = 1), "conf_level")
  expect_error(factor_scores(fit, conf_level = c(.9, .95)), "conf_level")

  # A residual variance estimated below zero (a Heywood case).
  set.seed(19)
  n <- 60
  f <- stats::rnorm(n)
  d <- data.frame(y1 = .95 * f + stats::rnorm(n, 0, .3),
                  y2 = .6 * f + stats::rnorm(n, 0, .8),
                  y3 = .5 * f + stats::rnorm(n, 0, .87))
  # lavaan warns about the negative variance itself.
  fit_h <- suppressWarnings(lavaan::cfa("f =~ y1 + y2 + y3", data = d,
                                        std.lv = TRUE))
  expect_lt(lavaan::lavInspect(fit_h, "est")$theta["y1", "y1"], 0)
  expect_warning(factor_scores(fit_h), "inadmissible.*y1")

  # A correlation of the factors beyond 1, which lavaan's post-estimation
  # check flags, was scored without a warning, although the regression
  # weights and both determinacies are built from it.
  set.seed(2)
  n <- 60
  e <- stats::rnorm(n)
  d2 <- as.data.frame(sapply(rep(.6, 6), function(l) {
    l * e + stats::rnorm(n, sd = .8)
  }))
  names(d2) <- paste0("y", 1:6)
  # lavaan warns that the covariance matrix of the factors is not positive
  # definite.
  fit_c <- suppressWarnings(lavaan::cfa("a =~ y1 + y2 + y3\nb =~ y4 + y5 + y6",
                                        data = d2, std.lv = TRUE))
  expect_gt(lavaan::lavInspect(fit_c, "est")$psi["a", "b"], 1)
  expect_warning(factor_scores(fit_c),
                 "a correlation of 1.07 between the factors a and b",
                 fixed = TRUE)
  # The residual covariance matrix of the continuous items, judged without
  # an auxiliary variable (a saturated correlate that shares much with the
  # factor makes the whole matrix indefinite while the items' is not).
  th <- matrix(c(.5, .6, 0, .6, .5, 0, 0, 0, 1), 3, 3,
               dimnames = list(c("y1", "y2", "z"), c("y1", "y2", "z")))
  est <- list(list(lambda = matrix(c(.7, .7, 0), 3, 1,
                                   dimnames = list(rownames(th), "f")),
                   theta = th,
                   psi = matrix(1, 1, 1, dimnames = list("f", "f"))))
  expect_identical(.fs_inadmissible(est, NA, "f"), paste0(
    "a residual covariance matrix of the continuous observed variables ",
    "that is not positive definite"))
  th_z <- matrix(c(.5, 0, .6, 0, .5, .6, .6, .6, 1), 3, 3,
                 dimnames = dimnames(th))
  est[[1]]$theta <- th_z
  expect_match(.fs_inadmissible(est, NA, "f"), "not positive definite")
  expect_identical(.fs_inadmissible(est, NA, "f", aux = "z"), character(0))
  expect_identical(.fs_inadmissible(est, NA, "f", ordered = c("y1", "y2"),
                                    aux = "z"), character(0))

  # An item regressed on a covariate is held by lavaan as a dummy latent
  # variable: its residual variance is in Psi, and its Theta entry is zero.
  # A negative one was scored without a warning, with a regression
  # determinacy of 1.56, though lavaan's post-estimation check flags it.
  set.seed(1)
  n <- 60
  x <- stats::rnorm(n)
  f <- 0.3 * x + stats::rnorm(n)
  dh <- data.frame(y1 = .95 * f + stats::rnorm(n, sd = .2),
                   y2 = .6 * f + stats::rnorm(n, sd = .8),
                   y3 = .6 * f + stats::rnorm(n, sd = .8), x = x)
  # lavaan warns about the negative variance itself.
  fit_d <- suppressWarnings(lavaan::sem("f =~ y1 + y2 + y3\nf ~ x\ny1 ~ x",
                                        data = dh, std.lv = TRUE))
  est_d <- lavaan::lavInspect(fit_d, "est")
  expect_identical(unname(est_d$theta["y1", "y1"]), 0)
  expect_lt(est_d$psi["y1", "y1"], 0)
  expect_false(suppressWarnings(lavaan::lavInspect(fit_d, "post.check")))
  expect_warning(factor_scores(fit_d, data = dh),
                 "a negative residual variance for y1", fixed = TRUE)
  # The residual covariance matrix judged is lavaan's, with the dummies'
  # block of Psi folded in; and Psi as a whole is judged when everything
  # else passes, as lavaan's check does.
  lam <- matrix(c(.7, .7, 0, 0, 0, 1), 3, 2,
                dimnames = list(c("y1", "y2", "y3"), c("f", "y3")))
  th0 <- diag(c(.5, .5, 0))
  dimnames(th0) <- list(rownames(lam), rownames(lam))
  ps <- matrix(c(1, 0, 0, .5), 2, 2, dimnames = list(colnames(lam),
                                                      colnames(lam)))
  e1 <- .fs_dummies(list(lambda = lam, theta = th0, psi = ps,
                         beta = matrix(c(0, .6, 0, 0), 2, 2,
                                       dimnames = dimnames(ps))), "f")
  expect_identical(unname(e1$theta_eff["y3", "y3"]), .5)
  expect_identical(unname(e1$lambda_eff["y3", ]), c(.6, 0))
  expect_identical(.fs_inadmissible(list(e1), NA, "f"), character(0))
  e1$psi["f", "y3"] <- e1$psi["y3", "f"] <- .9
  expect_identical(.fs_inadmissible(list(e1), NA, "f"), paste0(
    "a residual covariance matrix of the structural part of the model ",
    "(Psi) that is not positive definite"))
  e2 <- e1
  e2$psi["f", "y3"] <- e2$psi["y3", "f"] <- 0
  e2$psi["y3", "y3"] <- -.2
  e2 <- .fs_dummies(e2, "f")
  expect_identical(.fs_inadmissible(list(e2), NA, "f"),
                   "a negative residual variance for y3")
})

test_that("the broom verbs and the print method", {
  skip_if_not_installed("lavaan")
  d <- fs_mar_data()
  fit <- fs_fit_fiml(d)
  fs <- factor_scores(fit, data = d)
  td <- generics::tidy(fs)
  expect_identical(class(td), "data.frame")
  expect_identical(names(td), names(fs_shape))
  expect_null(attr(td, "patterns"))
  expect_identical(td$.derived_bartlett, fs$.derived_bartlett)
  gl <- generics::glance(fs)
  expect_identical(nrow(gl), 1L)
  expect_identical(gl$n_rows, nrow(d))
  expect_identical(gl$n_scored, nrow(d) - 2L)
  expect_identical(gl$n_factors, 3L)
  expect_identical(gl$n_patterns,
                   length(unique(fs$pattern[!is.na(fs$pattern)])))
  expect_identical(gl$method_bartlett, "Bartlett")
  expect_identical(gl$scored_from, "fit")
  out <- capture.output(print(fs))
  expect_true(any(grepl("regression and Bartlett scores, 3 factors", out)))
  expect_true(any(grepl("Determinacy", out)))
  expect_true(any(grepl("Skrondal & Laake", out)))
  expect_true(any(grepl("no observed item", out)))
  expect_true(any(grepl("Confidence level: 95%", out)))
  expect_true("(the first 3 rows of each factor)" %in% out)
  # One row of each factor is "the first row", not "the first 1 rows",
  # at the console and in a knitted report.
  out1 <- capture.output(print(fs, n = 1))
  expect_true("(the first row of each factor)" %in% out1)
  expect_false(any(grepl("first 1 rows", out1, fixed = TRUE)))
  expect_identical(sum(grepl("^ +[0-9]+ +(visual|textual|speed) ", out1)),
                   3L)
  foot1 <- DMAR:::.dmar_tbl_footer(fs, n = 1)
  expect_identical(unname(foot1[names(foot1) == "table_note"]),
                   "(the first row of each factor)")
  out2 <- capture.output(print(fs, n = 2))
  expect_true("(the first 2 rows of each factor)" %in% out2)
  # A subset keeps the class but not the attributes, and prints as an
  # ordinary table.
  sub <- fs[fs$row == 1, c("row", "factor", ".derived_bartlett")]
  expect_output(print(sub), ".derived_bartlett", fixed = TRUE)
  expect_identical(generics::glance(sub)$n_rows, 1L)
})

test_that("the interval columns name the level as it is", {
  skip_if_not_installed("lavaan")
  # sprintf("%g") wrote six significant digits, so a level of .9999999
  # named its columns "100% CI", a level only an unbounded interval has.
  d <- fs_mar_data()
  fit <- fs_fit_fiml(d)
  fs <- factor_scores(fit, data = d, conf_level = 0.9999999)
  out <- capture.output(print(fs, n = 1))
  expect_true(any(grepl("99.99999% CI", out, fixed = TRUE)))
  expect_false(any(grepl("100% CI", out, fixed = TRUE)))
})

test_that("the factor_scores() summary follows digits and dmar.digits", {
  skip_if_not_installed("lavaan")
  # Through the first versions of 1.1.0 the summary printed three decimal
  # places whatever print()'s digits or options(dmar.digits = ) asked for.
  d <- fs_mar_data()
  fs <- factor_scores(fs_fit_fiml(d), data = d)
  score_rows <- function(out) {
    grep("^ +[0-9]+ +(visual|textual|speed) ", out, value = TRUE)
  }
  three <- capture.output(print(fs, n = 1))
  five <- capture.output(print(fs, n = 1, digits = 5))
  expect_length(score_rows(three), 3L)
  expect_length(score_rows(five), 3L)
  expect_true(all(grepl("[0-9]\\.[0-9]{3} ", score_rows(three))))
  expect_false(any(grepl("[0-9]\\.[0-9]{4}", score_rows(three))))
  expect_true(all(grepl("[0-9]\\.[0-9]{5} ", score_rows(five))))
  # The determinacy table beneath it follows as well.
  det_line <- function(out) {
    out[grep("^ *pattern ", out)[1L] + 1L]
  }
  expect_match(det_line(three), "[0-9]\\.[0-9]{3}( |$)")
  expect_match(det_line(five), "[0-9]\\.[0-9]{5}( |$)")
  withr::local_options(dmar.digits = 5)
  expect_identical(capture.output(print(fs, n = 1)), five)
})

test_that("factor_scores() reports 2 rho^2 - 1 beside the regression determinacy", {
  skip_if_not_installed("lavaan")
  d <- fs_mar_data()
  fs <- factor_scores(fs_fit_fiml(d), data = d)
  det <- attr(fs, "determinacy")
  expect_equal(det$min_score_correlation, 2 * det$determinacy_regression^2 - 1)
  pt <- attr(fs, "patterns")
  ok <- is.finite(pt$determinacy_regression)
  expect_true(any(ok))
  expect_equal(pt$min_score_correlation[ok], 2 * pt$determinacy_regression[ok]^2 - 1)
  expect_true(all(is.na(pt$min_score_correlation[!ok])))
  out <- capture.output(print(fs))
  expect_true(any(grepl("min_r, the lowest correlation two sets of", out, fixed = TRUE)))
  expect_true(any(grepl("min_r", out[grepl("pattern", out)], fixed = TRUE)))
})
