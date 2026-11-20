# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# The @export / @exportS3Method tags on as_kable.dmar_tbl and
# knit_print.dmar_tbl register these methods once NAMESPACE is regenerated with
# devtools::document(). Register them here as well so the suite passes under a
# bare devtools::load_all() before that regeneration; registerS3method is
# idempotent, so this is a harmless no-op afterward.
if (exists("as_kable.dmar_tbl"))
  registerS3method("as_kable", "dmar_tbl", as_kable.dmar_tbl)
if (exists("knit_print.dmar_tbl") && requireNamespace("knitr", quietly = TRUE))
  registerS3method("knit_print", "dmar_tbl", knit_print.dmar_tbl,
                   envir = asNamespace("knitr"))

test_that("knit_print.dmar_tbl renders formatted values, not raw doubles", {
  skip_if_not_installed("knitr")
  x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
  kp <- knit_print.dmar_tbl(x)
  # A knit_asis object (what knitr places in the document).
  expect_true(inherits(kp, "knit_asis") || is.character(kp))
  txt <- as.character(kp)
  # The formatted display values appear; the raw full-precision double does not.
  expect_match(txt, "0.101", fixed = TRUE)
  expect_match(txt, "0.897", fixed = TRUE)
  expect_false(grepl("0.1005857", txt, fixed = TRUE))
  expect_false(grepl("0.8969414", txt, fixed = TRUE))
  # The confidence-level footer travels with the table.
  expect_match(txt, "Confidence level: 95%", fixed = TRUE)
})

test_that("as_kable returns a knitr_kable retaining all columns", {
  skip_if_not_installed("knitr")
  # Wide table: several typed columns must all survive.
  x <- ci_eta_squared(aov(len ~ supp * factor(dose), data = ToothGrowth))
  k <- as_kable(x)
  expect_s3_class(k, "knitr_kable")
  # Without kableExtra the kable is one string per line, not one string;
  # collapse so the content assertions hold under either rendering.
  txt <- paste(as.character(k), collapse = "\n")
  for (nm in names(x)) expect_match(txt, nm, fixed = TRUE)

  # Long table: the single value column formats to the display precision.
  k2 <- as_kable(ci_smd(smd = 0.5, n_1 = 50, n_2 = 50))
  expect_s3_class(k2, "knitr_kable")
  expect_match(paste(as.character(k2), collapse = "\n"), "0.101", fixed = TRUE)
})

test_that("knit_print() and as_kable() keep the decimal places of the console", {
  skip_if_not_installed("knitr")
  # Three decimal places with their trailing zeros (0.5 reads 0.500), and a
  # lower limit truncated at 0 reads 0.000, in a knitted report as at the
  # console; options(dmar.digits = ) moves all three alike.
  x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
  kp <- as.character(knit_print.dmar_tbl(x))
  pipe <- paste(as.character(as_kable(x, format = "pipe")), collapse = "\n")
  for (txt in c(kp, pipe)) {
    expect_match(txt, "|smd      |0.500 |", fixed = TRUE)
    expect_match(txt, "|ci_lower |0.101 |", fixed = TRUE)
  }
  r2 <- ci_R2(R2 = 0.05, N = 50, p = 5)
  expect_match(as.character(knit_print.dmar_tbl(r2)),
               "|ci_lower |0.000 |", fixed = TRUE)
  withr::local_options(dmar.digits = 5)
  kp5 <- as.character(knit_print.dmar_tbl(x))
  pipe5 <- paste(as.character(as_kable(x, format = "pipe")), collapse = "\n")
  for (txt in c(kp5, pipe5)) {
    expect_match(txt, "|smd      |0.50000 |", fixed = TRUE)
    expect_match(txt, "|ci_lower |0.10059 |", fixed = TRUE)
  }
  expect_true(any(grepl("smd +0\\.50000", utils::capture.output(print(x)))))
})

# A LaTeX special character that is not escaped: "%", "_", "&", "#", or
# "$" without a backslash before it.
bare_latex_special <- function(txt) {
  grepl("(?<!\\\\)[%_&#$]", txt, perl = TRUE)
}

test_that("as_kable escapes the caption for LaTeX and HTML", {
  # The 1.1.0 review: the caption went to knitr::kable() as it was, and
  # knitr does not escape a caption, so with format = "latex" (or the
  # default format in an Rnw document) the "%" of "Confidence level: 95%"
  # started a LaTeX comment that swallowed the closing brace, and the
  # document did not compile; an underscore in an item name stopped it as
  # well.
  skip_if_not_installed("knitr")
  x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
  tex <- paste(as.character(as_kable(x, format = "latex")), collapse = "\n")
  expect_match(tex, "\\caption{Confidence level: 95\\%}", fixed = TRUE)
  cap <- regmatches(tex, regexpr("\\\\caption\\{[^\n]*\\}", tex))
  expect_false(bare_latex_special(cap))
  # A footer with an underscored item name, a pairwise deletion line with
  # an ampersand, and the escapes LaTeX needs for each.
  lines <- c(conf_level = "Confidence level: 95%",
             residual_cor = "Largest residual correlation: 0.05 (item_4 with item_5)",
             missing = "Missing data: pairwise deletion (Asparouhov & Muthen, 2010) <")
  y <- DMAR:::.as_dmar_tbl(data.frame(term = c("a", "b"), value = c(1, 2)),
                           subclass = "dmar_caption_test")
  registerS3method(".dmar_tbl_footer", "dmar_caption_test",
                   function(x, ...) lines, envir = asNamespace("DMAR"))
  tex <- paste(as.character(as_kable(y, format = "latex")), collapse = "\n")
  cap <- regmatches(tex, regexpr("\\\\caption\\{[^\n]*\\}", tex))
  expect_false(bare_latex_special(cap))
  expect_match(cap, "95\\%; Largest", fixed = TRUE)
  expect_match(cap, "item\\_4 with item\\_5", fixed = TRUE)
  expect_match(cap, "Asparouhov \\& Muthen", fixed = TRUE)
  expect_match(cap, "2010) \\textless{}}", fixed = TRUE)
  # HTML: the caption's markup characters are entities.
  html <- paste(as.character(as_kable(y, format = "html")), collapse = "\n")
  expect_match(html, "<caption>Confidence level: 95%; ", fixed = TRUE)
  expect_match(html, "Asparouhov &amp; Muthen, 2010) &lt;</caption>",
               fixed = TRUE)
  # Markdown keeps the lines as they are, as before.
  pipe <- paste(as.character(as_kable(y, format = "pipe")), collapse = "\n")
  expect_match(pipe, "Table: Confidence level: 95%; Largest residual correlation: 0.05 (item_4 with item_5); Missing data: pairwise deletion (Asparouhov & Muthen, 2010) <",
               fixed = TRUE)
  # A table without footer lines has no caption.
  z <- DMAR:::.as_dmar_tbl(data.frame(term = "a", value = 1))
  expect_false(grepl("caption", paste(as.character(
    as_kable(z, format = "latex")), collapse = "\n"), fixed = TRUE))
  # A caption of the caller's own is refused in words, with or without
  # lines of the table's own.
  for (tab in list(x, z)) {
    expect_error(as_kable(tab, format = "latex", caption = "Mine"),
                 "writes the caption itself")
  }

  # knit_print() in an Rnw document (a LaTeX kable) escapes its paragraphs
  # and sets a small table of lines as verbatim text, which it leaves as
  # it is.
  lines <- c(lines, stats::setNames(c(" n  v_1", " 1  2%"),
                                    rep("own_table", 2L)))
  registerS3method(".dmar_tbl_footer", "dmar_caption_test",
                   function(x, ...) lines, envir = asNamespace("DMAR"))
  kp <- as.character(knit_print.dmar_tbl(y, format = "latex"))
  expect_match(kp, "\\end{tabular}\n\nConfidence level: 95\\%\n\n",
               fixed = TRUE)
  expect_match(kp, "(item\\_4 with item\\_5)", fixed = TRUE)
  expect_match(kp, "\\begin{verbatim}\n n  v_1\n 1  2%\n\\end{verbatim}",
               fixed = TRUE)
  paras <- strsplit(sub("\\\\begin\\{verbatim\\}.*$", "", kp), "\n\n")[[1]]
  paras <- paras[nzchar(paras)]
  expect_false(any(bare_latex_special(paras[-1L])))
  # In HTML, paragraphs and a preformatted block, escaped.
  kh <- as.character(knit_print.dmar_tbl(y, format = "html"))
  expect_match(kh, "<p>Confidence level: 95%</p>", fixed = TRUE)
  expect_match(kh, "(Asparouhov &amp; Muthen, 2010) &lt;</p>", fixed = TRUE)
  expect_match(kh, "<pre>\n n  v_1\n 1  2%\n</pre>", fixed = TRUE)
  # In Markdown, unchanged.
  km <- as.character(knit_print.dmar_tbl(y))
  expect_match(km, "\n\nConfidence level: 95%\n\n", fixed = TRUE)
  expect_match(km, "```\n n  v_1\n 1  2%\n```", fixed = TRUE)
})

test_that("results_sentence(ci_R2(...)) is the exact publication string", {
  x <- ci_R2(R2 = 0.25, N = 100, p = 5)
  expect_equal(results_sentence(x), "R2 = 0.25, 95% CI [0.08, 0.37]")
})

test_that("results_sentence honors label, digits, and conf_level", {
  x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
  expect_equal(results_sentence(x, label = "Cohen's d"),
               "Cohen's d = 0.50, 95% CI [0.10, 0.90]")
  # digits controls the reported precision.
  expect_equal(results_sentence(x, label = "d", digits = 3),
               "d = 0.500, 95% CI [0.101, 0.897]")
  # A non-default confidence level is reflected in the coverage.
  x90 <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.90)
  expect_match(results_sentence(x90), "90% CI", fixed = TRUE)
})

test_that("results_sentence writes a value that rounds to zero without a minus sign", {
  # As the table does: a lower limit of -0.0041, which the table shows to
  # three decimal places as -0.004, rounds to zero at the sentence's two
  # and reads 0.00, where through 1.0.x the sentence read -0.00. A limit
  # that rounds away from zero keeps its sign, and so does an infinite one.
  x <- DMAR:::.as_dmar_tbl(data.frame(
    term = c("smd", "ci_lower", "ci_upper"),
    value = c(0.27, -0.0041, 0.5512), stringsAsFactors = FALSE),
    conf_level = 0.95)
  expect_identical(results_sentence(x), "smd = 0.27, 95% CI [0.00, 0.55]")
  expect_identical(format(x)$value[2L], "-0.004")
  expect_identical(format(x, digits = 2)$value[2L], "0.00")
  x$value[2L] <- -0.0061
  expect_identical(results_sentence(x), "smd = 0.27, 95% CI [-0.01, 0.55]")
  expect_identical(results_sentence(x, digits = 3),
                   "smd = 0.270, 95% CI [-0.006, 0.551]")
  x$value[2:3] <- c(-Inf, 0.5512)
  expect_match(results_sentence(x), "-Inf", fixed = TRUE)
})

test_that("results_sentence works on a single row of a wide table", {
  x <- ci_eta_squared(aov(len ~ supp * factor(dose), data = ToothGrowth))
  s <- results_sentence(x[1, ], label = "eta squared")
  expect_match(s, "^eta squared = 0\\.[0-9]{2}, 95% CI \\[0\\.[0-9]{2}, 0\\.[0-9]{2}\\]$")
})

test_that("results_sentence appends a p-value when the table carries one", {
  # A p-value below format_p()'s floor is written as the bound it is.
  # Through 1.0.1 the sentence read "p = < 0.0001".
  w <- data.frame(effect = "a", estimate = 0.5, ci_lower = 0.1,
                  ci_upper = 0.9, p_value = 1e-7)
  w <- DMAR:::.as_dmar_tbl(w, conf_level = 0.95)
  expect_equal(results_sentence(w),
               "estimate = 0.50, 95% CI [0.10, 0.90], unadjusted p < 0.0001")
  w$p_value <- 0.0123
  expect_equal(results_sentence(w),
               "estimate = 0.50, 95% CI [0.10, 0.90], unadjusted p = 0.0123")
  # The same in a long table, whose p-value is a p_terms row: the t test
  # of two groups of 30 whose means differ by two standard deviations.
  tt <- summary_t_test(mean_1 = 2, sd_1 = 1, n_1 = 30, mean_2 = 0,
                       sd_2 = 1, n_2 = 30)
  expect_identical(results_sentence(tt),
                   "mean_difference = 2.00, 95% CI [1.48, 2.52], unadjusted p < 0.0001")
})

test_that("results_sentence errors clearly when there is no interval", {
  d <- DMAR:::.as_dmar_tbl(data.frame(term = c("mean", "sd"),
                                      value = c(1, 2)))
  # The error names what it looked for, and does not claim the table has
  # no interval (it may hold one under other names).
  expect_error(results_sentence(d),
               "reads a confidence interval from ci_lower and ci_upper")
})

test_that("results_sentence reads an interval held as ci_lower and ci_upper", {
  # The 1.1.0 review: a row of common_method_adjusted_r(), fa_confirmatory(), or any
  # tidy() result, whose interval is ci_lower / ci_upper, was told "This
  # table carries no confidence interval".
  R <- matrix(c(1, .5, .3, .5, 1, .4, .3, .4, 1), 3,
              dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  cm <- common_method_adjusted_r(R, N = 200, marker_r = .1)
  # The estimate is the estimate column, not marker_r, the first numeric
  # one.
  expect_identical(
    results_sentence(cm[1, ]),
    sprintf("estimate = %s, 95%% CI [%s, %s], unadjusted p %s",
            formatC(cm$estimate[1], format = "f", digits = 2),
            formatC(cm$ci_lower[1], format = "f", digits = 2),
            formatC(cm$ci_upper[1], format = "f", digits = 2),
            format_p(cm$p_value[1])))
  expect_match(results_sentence(cm[1, ], label = "adjusted r"),
               "^adjusted r = 0\\.44, 95% CI \\[0\\.32, 0\\.55\\], unadjusted p < 0\\.0001$")
  # A tidy() row gives the sentence the table itself gives, its term the
  # label and its conf_level column the level.
  x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
  expect_identical(results_sentence(generics::tidy(x)),
                   results_sentence(x))
  x90 <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.90)
  td90 <- generics::tidy(x90)
  attr(td90, "conf_level") <- NULL
  expect_match(results_sentence(td90), "^smd = 0\\.50, 90% CI ")
  # A coefficient row with a term column.
  w <- DMAR:::.as_dmar_tbl(data.frame(
    term = "b1", estimate = 0.3, se = 0.1, z_value = 3, p_value = 0.0027,
    ci_lower = 0.104, ci_upper = 0.496), conf_level = 0.95)
  expect_identical(results_sentence(w),
                   "b1 = 0.30, 95% CI [0.10, 0.50], unadjusted p = 0.0027")
  # Several rows still need one to be chosen.
  expect_error(results_sentence(cm), "single-row table")
})

test_that("results_sentence takes the estimate, not a count or a statistic, beside the limits", {
  # Through 1.0.1 the estimate of a wide table was its first numeric column
  # other than the limits, N, and the p-value, so a row that leads with a
  # count (responder_analysis()) or an F statistic (simple_effects_AB())
  # named that count or statistic as the estimate of the interval.
  w <- DMAR:::.as_dmar_tbl(data.frame(
    group = "a", n = 40, responders = 12, estimate = 0.3,
    ci_lower = 0.17, ci_upper = 0.46, chi_square = 1, df = 1,
    p_value = 0.3), conf_level = 0.95)
  expect_match(results_sentence(w), "^estimate = 0\\.30, 95% CI")
  w <- DMAR:::.as_dmar_tbl(data.frame(
    focal_factor = "A", F_value = 6.2, df_effect = 1, df_error = 36,
    p_value = 0.018, eta_squared_partial = 0.147, ci_lower = 0.006,
    ci_upper = 0.34), conf_level = 0.95)
  expect_identical(results_sentence(w),
                   "eta_squared_partial = 0.15, 95% CI [0.01, 0.34], unadjusted p = 0.0180")
})

test_that("results_sentence reports the level of a reliability coefficient's interval", {
  # Through 1.0.1 the confidence level was read from the attribute only,
  # and the reliability family carries it as a row, so a 90% interval for
  # omega was written as a 95% one.
  skip_if_not_installed("lavaan")
  set.seed(113)
  f <- rnorm(200)
  X <- sapply(c(.8, .7, .6, .7, .5),
              function(l) l * f + rnorm(200, sd = sqrt(1 - l^2)))
  colnames(X) <- paste0("y", 1:5)
  o <- suppressMessages(reliability_omega(data = X, conf_level = 0.90,
                                          ci_method = "ml",
                                          denominator = "model_implied"))
  expect_null(attr(o, "conf_level"))
  expect_match(results_sentence(o), "^estimate = 0\\.[0-9]{2}, 90% CI \\[")
})

test_that("results_sentence refuses limits that were never formed", {
  # A table of fixed shape keeps its interval NA when none was formed;
  # "95% CI [NA, NA]" is not a sentence.
  long <- DMAR:::.as_dmar_tbl(data.frame(
    term = c("est", "ci_lower", "ci_upper"), value = c(0.4, NA, NA)))
  expect_error(results_sentence(long), "limits are NA")
  wide <- DMAR:::.as_dmar_tbl(data.frame(
    label = "a", est = 0.4, ci_lower = 0.1, ci_upper = NA_real_))
  expect_error(results_sentence(wide), "limits are NA")
})

test_that("results_sentence reads an interval held in <estimate>_ci_lower and _ci_upper rows", {
  # The 1.1.0 review: a sem_fit_measures() table holds the RMSEA's interval
  # in rmsea_ci_lower and rmsea_ci_upper rows, which results_sentence() did
  # not read, so the table had no sentence.
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa(paste(
    "visual =~ t1_visual_perception + t2_cubes + t4_lozenges",
    "textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    "speed =~ t10_addition + t12_counting_groups_of_dots +",
    "  t13_straight_and_curved_capitals", sep = "\n"),
    data = holzinger_swineford)
  fm <- sem_fit_measures(fit, conf_level = 0.90)
  v <- stats::setNames(fm$value, fm$term)
  f3 <- function(z) formatC(z, format = "f", digits = 3)
  # The RMSEA, the first estimate with such rows, at the level the table
  # records (90 percent here), with no p-value: the table's p-values are
  # the tests of exact, close, and not close fit.
  expect_identical(
    results_sentence(fm, label = "RMSEA", digits = 3),
    sprintf("RMSEA = %s, 90%% CI [%s, %s]", f3(v[["rmsea"]]),
            f3(v[["rmsea_ci_lower"]]), f3(v[["rmsea_ci_upper"]])))
  # lavaan's own RMSEA and 90 percent interval for the fit, .092 [.071,
  # .114].
  lv <- lavaan::fitMeasures(fit, c("rmsea", "rmsea.ci.lower",
                                   "rmsea.ci.upper"))
  expect_identical(results_sentence(fm, label = "RMSEA", digits = 3),
                   sprintf("RMSEA = %s, 90%% CI [%s, %s]", f3(lv[[1]]),
                           f3(lv[[2]]), f3(lv[[3]])))
  expect_identical(results_sentence(fm, label = "RMSEA", digits = 3),
                   "RMSEA = 0.092, 90% CI [0.071, 0.114]")
  expect_match(results_sentence(fm), "^rmsea = 0\\.09, 90% CI ")
  # Gamma hat by name, and the level of another conf_level.
  expect_identical(
    results_sentence(fm, term = "gamma_hat", digits = 3),
    sprintf("gamma_hat = %s, 90%% CI [%s, %s]", f3(v[["gamma_hat"]]),
            f3(v[["gamma_hat_ci_lower"]]), f3(v[["gamma_hat_ci_upper"]])))
  expect_match(results_sentence(sem_fit_measures(fit, conf_level = 0.95)),
               "95% CI", fixed = TRUE)
  # A term without such rows is refused with the ones that have them.
  expect_error(results_sentence(fm, term = "cfi"),
               "the estimates that have such rows are \"rmsea\", \"gamma_hat\", \"adjusted_gamma_hat\"",
               fixed = TRUE)
  expect_error(results_sentence(ci_smd(smd = 0.5, n_1 = 50, n_2 = 50),
                                term = "smd"),
               "the table has no such rows", fixed = TRUE)
  expect_error(results_sentence(fm, term = c("rmsea", "gamma_hat")),
               "'term' must be NULL or a single term")

  # A fa_confirmatory() table records its RMSEA interval's level in
  # rmsea_ci_level (rmsea_conf_level, .90 by default), apart from the .95
  # of its parameter intervals; term = "rmsea" reads the interval at its own
  # level. Without term the table still needs one row chosen, and the
  # message names the term.
  k <- fa_confirmatory(holzinger_swineford, list(
    verbal = c("t6_paragraph_comprehension", "t7_sentence",
               "t9_word_meaning"),
    speed = c("t10_addition", "t12_counting_groups_of_dots",
              "t13_straight_and_curved_capitals")))
  expect_identical(attr(k, "conf_level"), 0.95)
  kv <- stats::setNames(k$estimate, k$term)
  expect_identical(
    results_sentence(k, term = "rmsea", digits = 3),
    sprintf("rmsea = %s, 90%% CI [%s, %s]", f3(kv[["rmsea"]]),
            f3(kv[["rmsea_ci_lower"]]), f3(kv[["rmsea_ci_upper"]])))
  expect_error(results_sentence(k), "e.g. term = \"rmsea\"", fixed = TRUE)

  # A table with an NA limit has no sentence, as for the other forms.
  na_rows <- DMAR:::.as_dmar_tbl(data.frame(
    term = c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper"),
    value = c(0.05, NA, 0.08)))
  expect_error(results_sentence(na_rows), "limits are NA")
})

test_that("results_sentence gives a selection of the RMSEA rows the RMSEA interval's own level", {
  # The 2026-11-09 assessment (R08): selecting the rmsea, rmsea_ci_lower,
  # and rmsea_ci_upper rows of a fa_confirmatory() table leaves the
  # rmsea_ci_level row out, and the sentence then took the level of the
  # parameter intervals (95 percent, beside the 90 percent limits it
  # printed) and refused the right level given as conf_level. The table's
  # rmsea_conf_level attribute, which the selection keeps, is now read.
  skip_if_not_installed("lavaan")
  factors <- list(
    visual  = c("t1_visual_perception", "t2_cubes", "t3_paper_form_board"),
    textual = c("t5_general_information", "t6_paragraph_comprehension",
                "t7_sentence"),
    speed   = c("t10_addition", "t11_code", "t12_counting_groups_of_dots"))
  hs <- holzinger_swineford[, unlist(factors)]
  keep <- c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper")
  sentence <- function(x, ...) {
    results_sentence(x, term = "rmsea", digits = 3, ...)
  }

  # Defaults: parameters at 95 percent, the RMSEA at 90 percent. The
  # limits are lavaan's 90 percent limits, .053 and .097.
  k <- fa_confirmatory(hs, factors)
  s <- k[k$term %in% keep, ]
  expect_false("rmsea_ci_level" %in% s$term)
  expect_identical(attr(s, "rmsea_conf_level"), 0.90)
  expect_identical(sentence(s), sentence(k))
  expect_match(sentence(s), "90% CI", fixed = TRUE)
  expect_identical(sentence(s), "rmsea = 0.075, 90% CI [0.053, 0.097]")
  fit <- fa_confirmatory(hs, factors, output = "fit")
  lv <- lavaan::fitMeasures(fit, c("rmsea", "rmsea.ci.lower",
                                   "rmsea.ci.upper"))
  expect_identical(
    sentence(s),
    sprintf("rmsea = %.3f, 90%% CI [%.3f, %.3f]", lv[[1]], lv[[2]], lv[[3]]))

  # The level given as conf_level is accepted when it is the RMSEA
  # interval's and refused, naming that level, when it is not.
  expect_identical(sentence(s, conf_level = 0.90), sentence(s))
  expect_error(
    sentence(s, conf_level = 0.95),
    "not the level this table records for its interval, 90%", fixed = TRUE)

  # Parameters at 99 percent and the RMSEA at 95 percent.
  k_99 <- fa_confirmatory(hs, factors, conf_level = 0.99,
                          rmsea_conf_level = 0.95)
  s_99 <- k_99[k_99$term %in% keep, ]
  expect_match(sentence(s_99), "95% CI", fixed = TRUE)
  expect_identical(sentence(s_99), sentence(k_99))

  # Both levels equal.
  k_90 <- fa_confirmatory(hs, factors, conf_level = 0.90)
  s_90 <- k_90[k_90$term %in% keep, ]
  expect_match(sentence(s_90), "90% CI", fixed = TRUE)
  expect_identical(sentence(s_90), sentence(k_90))

  # The bifactor and second-order tables record the level the same way.
  b <- fa_bifactor(hs, factors)
  s_b <- b[b$term %in% keep, ]
  expect_identical(sentence(s_b), sentence(b))
  expect_match(sentence(s_b), "90% CI", fixed = TRUE)
  so <- fa_second_order(hs, factors,
                        second_order = list(g = names(factors)))
  s_so <- so[so$term %in% keep, ]
  expect_identical(sentence(s_so), sentence(so))
  expect_match(sentence(s_so), "90% CI", fixed = TRUE)

  # A conf_level column on the interval's rows still comes before any
  # attribute, as for stacked tidy() rows, whose attributes are those of
  # the first frame of the stack.
  d <- s
  d$conf_level <- 0.80
  expect_identical(attr(d, "rmsea_conf_level"), 0.90)
  expect_match(sentence(d), "80% CI", fixed = TRUE)

  # Unchanged: a table with one level, and a selection that drops the
  # attributes, which has no level to state.
  fm <- sem_fit_measures(fit, conf_level = 0.95)
  expect_match(results_sentence(fm), "95% CI", fixed = TRUE)
  expect_error(sentence(subset(k, term %in% keep)),
               "records no confidence level")
})

test_that("a knitted table carries every line the console prints under it", {
  # The 1.1.0 review: the knitr methods wrote only the sum-of-squares and
  # confidence level lines, so a knitted omega lost its denominator and
  # residual lines and a knitted bootstrap interval its replication
  # counts. The console and the knitr methods now read one footer.
  skip_if_not_installed("knitr")
  skip_if_not_installed("lavaan")
  set.seed(113)
  f <- rnorm(200)
  X <- sapply(c(.8, .7, .6, .7, .5),
              function(l) l * f + rnorm(200, sd = sqrt(1 - l^2)))
  colnames(X) <- paste0("y", 1:5)
  footer_of_print <- function(x) {
    out <- capture.output(print(x))
    out[grepl("^(Sum of squares|Confidence level|Denominator|Largest residual|Bootstrap)", out)]
  }
  o <- suppressMessages(reliability_omega(data = X))
  a <- reliability_alpha(data = X, ci_method = "percentile", B = 50,
                         seed = 113)
  for (x in list(o, a, ci_smd(smd = 0.5, n_1 = 50, n_2 = 50))) {
    console <- footer_of_print(x)
    expect_gt(length(console), 0L)
    knitted <- as.character(knit_print.dmar_tbl(x))
    for (line in console) expect_match(knitted, line, fixed = TRUE)
    kab <- paste(as.character(as_kable(x, format = "pipe")), collapse = "\n")
    for (line in console) expect_match(kab, line, fixed = TRUE)
  }
  ko <- as.character(knit_print.dmar_tbl(o))
  expect_match(ko, "Denominator: observed composite variance (robust omega)",
               fixed = TRUE)
  expect_match(ko, "Largest residual correlation under the single-factor model:",
               fixed = TRUE)
  # One paragraph per line, so Markdown does not run them together.
  expect_match(ko, "(robust omega)\n\nLargest residual", fixed = TRUE)
  expect_match(as.character(knit_print.dmar_tbl(a)),
               "Bootstrap: 50 replications attempted, 50 returned a value",
               fixed = TRUE)
})

# The lines a display writes beside its table: every console line that is not
# a line of the table itself. The table is printed the way the console prints
# it (left aligned, or right aligned for a factor score summary), so the lines
# left over are the ones the result adds.
lines_beside_table <- function(x, right = FALSE) {
  console <- utils::capture.output(print(x))
  shown <- format(DMAR:::.dmar_tbl_display(x))
  tab <- utils::capture.output(print.data.frame(shown, row.names = FALSE,
                                                right = right))
  squish(setdiff(console, tab))
}
squish <- function(s) {
  s <- trimws(gsub("\\s+", " ", s))
  s[nzchar(s)]
}

test_that("every class that writes lines of its own carries them into a report", {
  # The 1.1.0 review, again: the stage 1 classes printed lines of their own
  # (which version of a fit measure was reported, what a comparison checked,
  # how an FMI was estimated, a factor score summary) from print methods that
  # knit_print() and as_kable() never called, and the AIPE line that names a
  # plan's missing data was written by print.dmar_tbl() alone. Every line now
  # comes from .dmar_tbl_footer(), and every display shows the table the
  # console shows (.dmar_tbl_display()).
  skip_if_not_installed("knitr")
  skip_if_not_installed("lavaan")
  hs <- holzinger_swineford
  fit <- lavaan::cfa(
    "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
     textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
     speed   =~ t10_addition + t12_counting_groups_of_dots +
                t13_straight_and_curved_capitals",
    data = hs, std.lv = TRUE)
  one <- lavaan::cfa(
    "g =~ t1_visual_perception + t2_cubes + t4_lozenges +
          t6_paragraph_comprehension + t7_sentence + t9_word_meaning +
          t10_addition + t12_counting_groups_of_dots +
          t13_straight_and_curved_capitals",
    data = hs, std.lv = TRUE)
  spatial <- lavaan::cfa(
    "spatial =~ t1_visual_perception + t2_cubes + t3_paper_form_board +
                t4_lozenges + t25_paper_form_board_r + t26_flags",
    data = hs, std.lv = TRUE, missing = "ml")
  set.seed(113)
  objs <- list(
    sem_fit_measures = sem_fit_measures(fit),
    equivalence_rmsea = equivalence_rmsea(85.3055, 24, 301, rmsea_bound = 0.08),
    sem_residual_fit = sem_residual_fit(fit, conf_level = 0.90,
                                        close_fit = 0.05),
    sem_information_criteria = sem_information_criteria(fit),
    compare_sem = compare_sem(list(three = fit, one = one)),
    fraction_missing_information = fraction_missing_information(spatial),
    factor_scores = factor_scores(spatial, data = hs),
    ss_aipe_reliability = ss_aipe_reliability(
      model = "Parallel", type = "Normal Theory", width = .1, n_items = 6,
      true_item_covariance = .3, true_residual_var = .2, fmi = .30),
    ss_aipe_reliability_sensitivity = ss_aipe_reliability_sensitivity(
      true_reliability = 0.80, specified_N = 60, n_items = 4, width = 0.15,
      estimator = "alpha", G = 2, missing_prop = 0.2),
    ss_aipe_c = ss_aipe_c(true_var_error = 1, c_weights = c(1, -1),
                          width = 0.8, assurance = 0.80, true_kurtosis = -1.2),
    df_missing = df_missing(c(slope = 0.5, sparse = 0.95), N = 25, k = 2),
    df_barnard_rubin = df_barnard_rubin(between = 0.04, within = 0.1, m = 5,
                                        df_complete = 23))
  # A line each class must show, so the check below cannot pass vacuously.
  expected <- c(
    sem_fit_measures = "Version: standard",
    equivalence_rmsea = "Not close fit (RMSEA >= 0.08) is not rejected",
    sem_residual_fit = "Residual covariance: normal theory",
    sem_information_criteria = "Smaller is preferred",
    compare_sem = "Checked: the same cases and data values",
    fraction_missing_information = "the third estimator of Chen & Savalei, 2021",
    factor_scores = "Determinacy, the correlation of each score",
    ss_aipe_reliability = "with complete data, N = ",
    ss_aipe_reliability_sensitivity = "Evaluated at the specified N",
    ss_aipe_c = "allows for errors with an excess kurtosis of -1.2",
    df_missing = "von Hippel's (2016, Eq. 31) effective sample size rule",
    df_barnard_rubin = "Barnard and Rubin's (1999) rule for m = 5")
  expect_setequal(names(expected), names(objs))
  for (nm in names(objs)) {
    x <- objs[[nm]]
    beside <- lines_beside_table(x, right = nm == "factor_scores")
    knitted <- squish(paste(as.character(knit_print.dmar_tbl(x)),
                            collapse = "\n"))
    kab <- squish(paste(as.character(as_kable(x, format = "pipe")),
                        collapse = "\n"))
    footer <- squish(paste(DMAR:::.dmar_tbl_footer(x), collapse = "\n"))
    for (line in beside) {
      # Each line comes from the one source, and reaches both knitr methods.
      expect_true(grepl(line, footer, fixed = TRUE), info = paste(nm, line))
      expect_true(grepl(line, knitted, fixed = TRUE), info = paste(nm, line))
      expect_true(grepl(line, kab, fixed = TRUE), info = paste(nm, line))
    }
    expect_true(any(grepl(expected[[nm]], beside, fixed = TRUE)),
                info = nm)
  }

  # The knitted table is the table the console shows: a factor score
  # summary shows the first three rows of each factor, not all 301 rows,
  # and the FMI table leaves out the operator column.
  fs <- objs$factor_scores
  k <- as.character(as_kable(fs, format = "pipe"))
  expect_identical(sum(grepl("^\\|[0-9]", k)), 3L)
  expect_match(as.character(knit_print.dmar_tbl(fs)),
               "^\n\nFactor scores: regression and Bartlett scores")
  kf <- paste(as.character(as_kable(objs$fraction_missing_information,
                                    format = "pipe")), collapse = "\n")
  expect_false(grepl("|op ", kf, fixed = TRUE))
  # The determinacy table keeps its columns in a knitted report.
  expect_match(as.character(knit_print.dmar_tbl(fs)),
               "```\n pattern +n +factor items", perl = TRUE)
})

test_that("footer blocks are spaced alike at the console and set apart when knitted", {
  skip_if_not_installed("knitr")
  x <- DMAR:::.as_dmar_tbl(data.frame(term = c("a", "b"), value = c(1, 2.5)),
                           conf_level = 0.95, subclass = "dmar_footer_test")
  # A class adds its lines through a method of the internal generic; the
  # header and the note on the table go where they belong wherever the
  # method puts them.
  local_footer <- function(x, ...) {
    out <- DMAR:::.dmar_footer_add(NextMethod(), "Two rows", "header")
    out <- DMAR:::.dmar_footer_add(out, "(both rows)", "table_note")
    out <- DMAR:::.dmar_footer_add(out, c("First line", "Second line"), "own")
    out <- DMAR:::.dmar_footer_add(out, c(" n  v", " 1  2"), "own_table")
    out
  }
  registerS3method(".dmar_tbl_footer", "dmar_footer_test", local_footer,
                   envir = asNamespace("DMAR"))
  out <- utils::capture.output(print(x))
  expect_identical(out, c("Two rows", "", " term value", " a    1    ",
                          " b    2.500", "(both rows)", "",
                          "Confidence level: 95%", "", "First line",
                          "Second line", " n  v", " 1  2"))
  kp <- as.character(knit_print.dmar_tbl(x))
  expect_match(kp, "^\n\nTwo rows\n\n\\|term")
  expect_match(kp, "(both rows)\n\nConfidence level: 95%\n\nFirst line\n\nSecond line\n\n```\n n  v\n 1  2\n```",
               fixed = TRUE)
  cap <- DMAR:::.dmar_footer_caption(DMAR:::.dmar_tbl_footer(x))
  expect_identical(cap, paste0("Two rows; (both rows); Confidence level: 95%; ",
                               "First line; Second line; n v; 1 2"))
})

test_that("results_sentence states no confidence level a table does not record", {
  # The 1.1.0 review: a table that recorded no level was taken to be at
  # 95%, so a 90% interval was written "95% CI". A table that records none
  # is now refused unless the caller gives the level, and a level that
  # differs from the recorded one is refused.
  hand <- data.frame(term = "b1", estimate = 0.3, ci_lower = 0.15,
                     ci_upper = 0.45)
  expect_error(results_sentence(hand), "records no confidence level")
  expect_error(results_sentence(hand), "conf_level = 0.95", fixed = TRUE)
  expect_identical(results_sentence(hand, conf_level = 0.90),
                   "b1 = 0.30, 90% CI [0.15, 0.45]")
  rec <- DMAR:::.as_dmar_tbl(hand, conf_level = 0.95)
  expect_identical(results_sentence(rec), "b1 = 0.30, 95% CI [0.15, 0.45]")
  expect_identical(results_sentence(rec, conf_level = 0.95),
                   results_sentence(rec))
  expect_error(results_sentence(rec, conf_level = 0.90),
               "is not the level this table records for its interval, 95%")
  expect_error(results_sentence(rec, conf_level = 1.5),
               "'conf_level' must be NULL or a single number in \\(0, 1\\)")
  # A tidy() row whose conf_level is NA records none.
  na_row <- hand
  na_row$conf_level <- NA_real_
  expect_error(results_sentence(na_row), "records no confidence level")
})

test_that("a 90% interval is written as one through every route", {
  # The 1.1.0 review, at conf_level = .90: the coefficient table of mlmr()
  # and the tidy() results of the reliability family, fa_confirmatory(),
  # content_validity_index(), adjusted_means(), icc(), cohen_kappa(), and
  # every long table recorded no level, and each sentence said "95% CI"
  # around the 90% limits.
  skip_if_not_installed("lavaan")
  hs <- holzinger_swineford
  f <- t6_paragraph_comprehension ~ t7_sentence + t9_word_meaning
  m90 <- suppressWarnings(mlmr(f, data = hs, conf_level = 0.90))
  expect_identical(attr(m90$coef_table, "conf_level"), 0.90)
  expect_identical(results_sentence(m90$coef_table[2, ]),
                   "t7_sentence = 0.32, 90% CI [0.26, 0.38], unadjusted p < 0.0001")
  td <- generics::tidy(m90, conf.int = TRUE)
  expect_identical(td$conf_level, rep(0.90, nrow(td)))
  expect_identical(results_sentence(td[2, ]),
                   results_sentence(m90$coef_table[2, ]))
  # The reliability family: the tidy row says what the table says.
  o <- reliability_omega(hs[, c("t6_paragraph_comprehension", "t7_sentence",
                                "t9_word_meaning")], ci_method = "ml",
                         denominator = "model_implied", conf_level = 0.90)
  expect_identical(results_sentence(generics::tidy(o)),
                   sub("^estimate", "omega", results_sentence(o)))
  expect_match(results_sentence(generics::tidy(o)), "^omega = 0\\.87, 90% CI")
  # fa_confirmatory(): the phi row of the table and of its tidy() result.
  k <- suppressWarnings(fa_confirmatory(hs, list(
    verbal = c("t5_general_information", "t6_paragraph_comprehension",
               "t7_sentence"),
    spatial = c("t1_visual_perception", "t2_cubes", "t4_lozenges")),
    conf_level = 0.90))
  tk <- generics::tidy(k)
  expect_identical(results_sentence(tk[tk$term == "phi_verbal_spatial", ]),
                   results_sentence(k[k$term == "phi_verbal_spatial", ]))
  expect_match(results_sentence(tk[tk$term == "phi_verbal_spatial", ]),
               "90% CI", fixed = TRUE)
  # A long table through the default tidy(), which keeps the level.
  p90 <- ci_proportion(successes = 30, N = 100, conf_level = 0.90)
  expect_identical(results_sentence(generics::tidy(p90)),
                   results_sentence(p90))
  expect_match(results_sentence(p90), "90% CI", fixed = TRUE)
  # A wide table through the default tidy(), which adds the level beside
  # the interval.
  ck <- cohen_kappa(rater_1 = rep(1:3, 20), rater_2 = rep(c(1:3, 3:1), 10),
                    conf_level = 0.90)
  tck <- generics::tidy(ck)
  expect_identical(tck$conf_level, rep(0.90, nrow(tck)))
  expect_match(results_sentence(tck[1, ]), "90% CI", fixed = TRUE)
  ratings <- matrix(c(9, 2, 5, 8, 6, 1, 3, 2, 8, 4, 6, 8, 7, 1, 2, 6, 10,
                      5, 6, 9, 6, 2, 4, 7), ncol = 4, byrow = TRUE)
  ic <- icc(ratings, conf_level = 0.90)
  expect_match(results_sentence(generics::tidy(ic)[1, ]), "90% CI",
               fixed = TRUE)
  cvi <- content_validity_index(matrix(c(4, 4, 3, 4, 3, 1, 4, 3), nrow = 2),
                                conf_level = 0.90)
  expect_match(results_sentence(generics::tidy(cvi)[1, ]), "90% CI",
               fixed = TRUE)
  # The tidy() row names its estimate as the table does ("i_cvi = 1.00"),
  # not by the item ("item_1 = 1.00").
  expect_identical(results_sentence(generics::tidy(cvi)[1, ]),
                   results_sentence(cvi[1, ]))
  expect_match(results_sentence(generics::tidy(cvi)[1, ]), "^i_cvi = ")
})

test_that("the confidence level is written in full whatever getOption(\"digits\") is", {
  # The level was formatted with getOption("digits") significant digits,
  # so under options(digits = 2) a 97.5% interval was labeled "98%" at the
  # console, in a knitted caption, and in its sentence.
  x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.975)
  withr::local_options(digits = 2)
  expect_true("Confidence level: 97.5%" %in% utils::capture.output(print(x)))
  expect_match(results_sentence(x), "^smd = 0\\.50, 97\\.5% CI \\[")
  cap <- utils::capture.output(print(as_kable(x, format = "pipe")))
  expect_true(any(grepl("Confidence level: 97.5%", cap, fixed = TRUE)))
  # A level below 1 never reads 100, at any setting.
  y <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.9995)
  expect_true("Confidence level: 99.95%" %in%
                utils::capture.output(print(y)))
  withr::local_options(digits = 7)
  z <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.999999995)
  expect_true("Confidence level: 99.9999995%" %in%
                utils::capture.output(print(z)))
})

test_that("the confidence level is written under a comma decimal mark", {
  # The test of whether the display read 100 parsed the formatted text with
  # as.numeric(), which reads "97,5" as NA under options(OutDec = ","), so
  # every level that is not a whole percent stopped with "missing value
  # where TRUE/FALSE needed": the printed table, its caption, its sentence,
  # a plot label, and every Bonferroni line (whose per-interval level is
  # never whole). The level follows OutDec, as the table's numbers do.
  withr::local_options(OutDec = ",")
  x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.975)
  expect_no_warning(out <- utils::capture.output(print(x)))
  expect_true("Confidence level: 97,5%" %in% out)
  expect_match(results_sentence(x), "^smd = 0,50, 97,5% CI \\[")
  cap <- utils::capture.output(print(as_kable(x, format = "pipe")))
  expect_true(any(grepl("Confidence level: 97,5%", cap, fixed = TRUE)))
  expect_identical(
    plot_ci(estimate = 0.4, lower = 0.1, upper = 0.7,
            conf_level = 0.975)$labels$title,
    "97,5% Confidence Interval")
  expect_identical(.format_conf_pct(c(0.95, 0.975, 1 - 0.001 / 28)),
                   c("95", "97,5", "99,99643"))
  expect_identical(.format_conf_pct(1 - 0.001 / 28, digits = 4L), "99,996")
  withr::local_seed(2)
  g <- factor(rep(letters[1:4], each = 5))
  v <- stats::rnorm(20)
  ct <- contrast_test(stats::aov(v ~ g), adjust = "bonferroni")
  expect_no_warning(out <- utils::capture.output(print(ct)))
  expect_true(any(grepl(paste0(
    "Confidence level: 95%, simultaneous across the 6 comparisons ",
    "(Bonferroni: each interval at 99,17%)"), out, fixed = TRUE)))
})

test_that("every interval formed from alpha_lower and alpha_upper records its level", {
  # Through 1.0.1 an interval formed from tail areas recorded no level (or,
  # for ci_eta_squared_generalized(), the default 0.95 beside a 90%
  # interval), so the console printed none and results_sentence() wrote
  # "95% CI". The level is now 1 - alpha_lower - alpha_upper.
  a <- list(conf_level = NULL, alpha_lower = 0.05, alpha_upper = 0.05)
  calls <- list(
    ci_c = list(means = c(1, 2, 3), s_anova = 1, c_weights = c(1, -1, 0),
                n_per_group = c(10, 10, 10)),
    ci_cv = list(cv = 0.3, N = 30),
    ci_mahalanobis = list(D2 = 1, n_1 = 20, n_2 = 20, p = 3),
    ci_r = list(r = 0.3, N = 50),
    ci_R = list(R = 0.5, df_1 = 3, df_2 = 50),
    ci_rc = list(reg_coef_j = 0.3, se_reg_coef_j = 0.1, sd_Y = 1, sd_X = 1,
                 N = 100, p = 3),
    ci_reg_coef = list(reg_coef_j = 0.3, se_reg_coef_j = 0.1, sd_Y = 1,
                       sd_X = 1, N = 100, p = 3),
    ci_rmsea = list(rmsea = 0.06, df = 24, N = 300),
    ci_sc = list(means = c(1, 2, 3), s_anova = 1, c_weights = c(1, -1, 0),
                 n_per_group = c(10, 10, 10)),
    ci_sm = list(sm = 0.5, N = 30),
    ci_smd = list(smd = 0.5, n_1 = 30, n_2 = 30),
    ci_smd_c = list(smd_c = 0.5, n_C = 30, n_E = 30),
    ci_src = list(std_reg_coef_j = 0.3, se_std_reg_coef_j = 0.1, N = 100,
                  p = 3),
    ci_eta_squared = list(F_value = 5, df_effect = 2, df_error = 27, N = 30),
    ci_omega_squared = list(F_value = 5, df_effect = 2, df_error = 27,
                            N = 30))
  for (fn in names(calls)) {
    res <- suppressWarnings(do.call(fn, c(calls[[fn]], a)))
    expect_equal(attr(res, "conf_level"), 0.90, tolerance = 1e-12,
                 info = fn)
    expect_true(any(utils::capture.output(print(res)) ==
                      "Confidence level: 90%"), info = fn)
  }
  expect_match(results_sentence(ci_r(r = 0.3, N = 50, conf_level = NULL,
                                     alpha_lower = 0.05,
                                     alpha_upper = 0.05)),
               "^r = 0\\.30, 90% CI \\[")
  # The generalized eta squared recorded the default 0.95 beside the
  # interval the tail areas define.
  fit <- aov(len ~ supp * factor(dose), data = ToothGrowth)
  g <- suppressWarnings(ci_eta_squared_generalized(
    fit, method = "parametric", alpha_lower = 0.05, alpha_upper = 0.05))
  expect_equal(attr(g, "conf_level"), 0.90, tolerance = 1e-12)
  g90 <- suppressWarnings(ci_eta_squared_generalized(
    fit, method = "parametric", conf_level = 0.90))
  expect_equal(g$ci_lower, g90$ci_lower, tolerance = 1e-12)
})

test_that("an interval with unequal tail areas says so on every route", {
  # The level of an interval formed from alpha_lower and alpha_upper is its
  # coverage, 1 - alpha_lower - alpha_upper, and that is all the table
  # recorded: a .01/.04 interval printed "Confidence level: 95%" and was
  # written "95% CI [L, U]", as the equal-tailed interval is, and a
  # one-sided interval was written with its uninformative limit,
  # "95% CI [-Inf, 0.83]" or, for a bounded parameter, "[0.00, 0.08]".
  lines_of <- function(x) utils::capture.output(print(x))
  x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, alpha_lower = 0.01,
              alpha_upper = 0.04)
  expect_identical(attr(x, "conf_level"), 1 - 0.01 - 0.04)
  expect_identical(attr(x, "tail_areas"), c(lower = 0.01, upper = 0.04))
  expect_true(paste0("Confidence level: 95% (unequal tail areas: 0.01 ",
                     "below the lower limit, 0.04 above the upper limit)")
              %in% lines_of(x))
  s <- "smd = 0.50, 95% CI [0.03, 0.85] (unequal tail areas: 0.01 below, 0.04 above)"
  expect_identical(results_sentence(x), s)
  td <- generics::tidy(x)
  expect_named(td, c("term", "estimate", "ci_lower", "ci_upper",
                     "conf_level", "alpha_lower", "alpha_upper"))
  expect_identical(c(td$alpha_lower, td$alpha_upper), c(0.01, 0.04))
  expect_identical(results_sentence(td), s)
  expect_identical(results_sentence(generics::glance(x)), s)
  cap <- utils::capture.output(print(as_kable(x, format = "pipe")))
  expect_true(any(grepl("unequal tail areas: 0.01 below the lower limit",
                        cap, fixed = TRUE)))
  # One-sided, from either tail.
  up <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, alpha_lower = 0,
               alpha_upper = 0.05)
  expect_true(paste0("Confidence level: 95% (one-sided: the upper limit is ",
                     "the 95% upper confidence bound, with tail area 0.05 ",
                     "above it and none below the lower limit)")
              %in% lines_of(up))
  expect_identical(results_sentence(up),
                   "smd = 0.50, one-sided 95% upper confidence bound 0.83")
  expect_identical(results_sentence(generics::tidy(up)),
                   results_sentence(up))
  lo <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, alpha_lower = 0.05,
               alpha_upper = 0)
  expect_match(results_sentence(lo),
               "^smd = 0\\.50, one-sided 95% lower confidence bound 0\\.1")
  expect_true(paste0("Confidence level: 95% (one-sided: the lower limit is ",
                     "the 95% lower confidence bound, with tail area 0.05 ",
                     "below it and none above the upper limit)")
              %in% lines_of(lo))
  # A bounded parameter's one-sided interval has a limit at the bound (0
  # for the RMSEA), which looked exactly like a two-sided interval.
  rm <- suppressMessages(ci_rmsea(rmsea = 0.06, df = 24, N = 301,
                                  alpha_lower = 0, alpha_upper = 0.05))
  expect_identical(rm$value[rm$term == "ci_lower"], 0)
  expect_identical(results_sentence(rm),
                   "rmsea = 0.06, one-sided 95% upper confidence bound 0.08")
  # A table with achieved tail areas in prob_less and prob_greater columns
  # (ci_R2) says the same beneath it and in its sentence, and its tidy()
  # keeps the split.
  r2 <- ci_R2(R2 = 0.25, N = 100, p = 5, conf_level = NULL,
              alpha_lower = 0.01, alpha_upper = 0.04)
  expect_true(any(grepl("(unequal tail areas: 0.01 below the lower limit",
                        lines_of(r2), fixed = TRUE)))
  expect_match(results_sentence(r2),
               "95% CI \\[.*\\] \\(unequal tail areas: 0.01 below, 0.04 above\\)$")
  expect_identical(results_sentence(generics::tidy(r2)),
                   results_sentence(r2))
  # The ANOVA intervals, through the table, a row of it, and tidy().
  fit <- aov(len ~ supp * factor(dose), data = ToothGrowth)
  es <- suppressWarnings(ci_eta_squared(fit, alpha_lower = 0.01,
                                        alpha_upper = 0.04))
  tes <- generics::tidy(es)
  # Three effects, each with a total (no interval) and a partial value.
  expect_identical(tes$alpha_lower, rep(0.01, 6))
  expect_identical(results_sentence(es[1, ]), results_sentence(tes[4, ]))
  expect_match(results_sentence(es[1, ]), "(unequal tail areas: 0.01 below, 0.04 above)",
               fixed = TRUE)
  # Rows of stacked tidy() results give the split of their own interval.
  a <- generics::tidy(ci_c(means = c(1, 2, 3), s_anova = 1,
                           c_weights = c(1, -1, 0), n_per_group = 10,
                           alpha_lower = 0.01, alpha_upper = 0.04))
  b <- generics::tidy(ci_c(means = c(1, 2, 3), s_anova = 1,
                           c_weights = c(1, -1, 0), n_per_group = 10,
                           alpha_lower = 0.04, alpha_upper = 0.01))
  expect_identical(a$alpha_lower, rep(0.01, 3))
  ab <- rbind(a, b)
  expect_match(results_sentence(ab[4:6, ]),
               "(unequal tail areas: 0.04 below, 0.01 above)", fixed = TRUE)
  expect_match(results_sentence(ab[1:3, ]),
               "(unequal tail areas: 0.01 below, 0.04 above)", fixed = TRUE)
  # Equal tails, given either way, add nothing, and no column: the line
  # and the sentence of 1.0.1, at the level the tails define.
  eq <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, alpha_lower = 0.05,
               alpha_upper = 0.05)
  expect_null(attr(eq, "tail_areas"))
  expect_true("Confidence level: 90%" %in% lines_of(eq))
  expect_identical(results_sentence(eq), "smd = 0.50, 90% CI [0.16, 0.83]")
  expect_named(generics::tidy(eq), c("term", "estimate", "ci_lower",
                                     "ci_upper", "conf_level"))
  # One-sidedness is read from a tail area of 0, not from a limit: an
  # equal-tailed interval whose lower limit is 0 because F fell below its
  # critical value (Steiger, 2004) is the two-sided interval it is.
  cl <- suppressWarnings(ci_eta_squared(F_value = 0.5, df_effect = 2,
                                        df_error = 27, N = 30))
  expect_identical(cl$eta_squared_partial_ci_lower, 0)
  expect_match(results_sentence(cl),
               "^eta_squared_partial = 0\\.04, 95% CI \\[0\\.00, ")
  expect_true("Confidence level: 95%" %in% lines_of(cl))
})

test_that("results_sentence reports the estimate a long table's interval belongs to", {
  # The 1.1.0 review: a long table's estimate was its first row that is
  # not a limit, and two tables list another quantity first, so their
  # sentences paired it with the interval of the estimate after it.
  x <- correction_for_attenuation(r = .3, reliability_x = .8,
                                  reliability_y = .7, N = 200)
  expect_identical(attr(x, "ci_estimate"), "correlation_corrected")
  expect_identical(results_sentence(x),
                   "correlation_corrected = 0.40, 95% CI [0.22, 0.56]")
  # The limits are the Fisher's Z interval for r divided by the
  # attenuation factor, the interval of the corrected correlation.
  lims <- tanh(atanh(.3) + c(-1, 1) * qnorm(.975) / sqrt(197)) / sqrt(.56)
  expect_equal(x$value[x$term %in% c("ci_lower", "ci_upper")], lims,
               tolerance = 1e-12)
  expect_identical(results_sentence(x, label = "corrected r"),
                   "corrected r = 0.40, 95% CI [0.22, 0.56]")
  ps <- proportion_of_superiority(smd = .5, n_1 = 50, n_2 = 50)
  expect_identical(results_sentence(ps),
                   "proportion_of_superiority = 0.69, 95% CI [0.54, 0.82]")
  # tidy() names the estimate and the level in its rows, so its sentence
  # is the table's.
  expect_identical(results_sentence(generics::tidy(x)), results_sentence(x))
  expect_identical(results_sentence(generics::tidy(ps)),
                   results_sentence(ps))
  # Without an interval there is nothing to name.
  expect_null(attr(proportion_of_superiority(smd = .5), "ci_estimate"))
  # A name that is not a row of the table falls back to the first row.
  odd <- DMAR:::.as_dmar_tbl(data.frame(
    term = c("est", "ci_lower", "ci_upper"), value = c(.5, .2, .8)),
    conf_level = 0.95)
  attr(odd, "ci_estimate") <- "not_a_row"
  expect_identical(results_sentence(odd), "est = 0.50, 95% CI [0.20, 0.80]")
})

test_that("rows of stacked tidy() results give the sentence of their own table", {
  # The 1.1.0 review: the tidy() of a long table carried its level and the
  # name of its estimate as attributes, which rbind() keeps from the first
  # frame of a stack, so the rows of a 99% interval stacked under a 90%
  # one were written "90% CI", and the rows of the corrected correlation
  # under another table were written with the observed correlation.
  # Through 1.0.1 the tidy() rows recorded no level at all.
  a <- generics::tidy(ci_proportion(12, 40, conf_level = .90))
  b <- generics::tidy(ci_proportion(30, 60, conf_level = .99))
  # The level sits on the rows of the interval, NA on the others.
  expect_identical(a$term, c("ci_lower", "proportion", "ci_upper",
                             "successes", "N"))
  expect_identical(a$conf_level, c(.90, .90, .90, NA, NA))
  expect_null(attr(a, "conf_level"))
  ab <- rbind(a, b)
  expect_identical(results_sentence(ab[6:10, ]),
                   "proportion = 0.50, 99% CI [0.34, 0.66]")
  expect_identical(results_sentence(ab[6:10, ]),
                   results_sentence(ci_proportion(30, 60, conf_level = .99)))
  expect_identical(results_sentence(ab[1:5, ]), results_sentence(a))
  expect_match(results_sentence(a), "90% CI", fixed = TRUE)
  x <- generics::tidy(correction_for_attenuation(
    r = .3, reliability_x = .8, reliability_y = .7, N = 200))
  expect_identical(x$conf_level[x$term == "correlation_observed"],
                   NA_real_)
  ax <- rbind(a, x)
  expect_identical(results_sentence(ax[-(1:5), ]),
                   "correlation_corrected = 0.40, 95% CI [0.22, 0.56]")
  # Rows that keep the limits without their estimate are refused rather
  # than paired with another row.
  expect_error(results_sentence(x[x$term != "correlation_corrected", ]),
               "not its estimate")
  # An interval held in rows named for its estimate is marked at its own
  # level, so the RMSEA rows of a 95% table stacked under a 90% one give
  # 95%.
  f90 <- generics::tidy(sem_fit_measures(chisq = 85.3, df = 24, N = 301,
                                         p = 9, conf_level = 0.90))
  f95 <- generics::tidy(sem_fit_measures(chisq = 85.3, df = 24, N = 301,
                                         p = 9, conf_level = 0.95))
  expect_identical(f90$conf_level[f90$term %in% c("rmsea", "rmsea_ci_lower",
                                                  "rmsea_ci_upper")],
                   rep(0.90, 3L))
  expect_true(is.na(f90$conf_level[f90$term == "chi_square"]))
  expect_identical(results_sentence(f90),
                   "rmsea = 0.09, 90% CI [0.07, 0.11]")
  ff <- rbind(f90, f95)
  expect_identical(results_sentence(ff[-seq_len(nrow(f90)), ]),
                   results_sentence(sem_fit_measures(chisq = 85.3, df = 24,
                                                     N = 301, p = 9,
                                                     conf_level = 0.95)))
  expect_match(results_sentence(ff[-seq_len(nrow(f90)), ]), "95% CI",
               fixed = TRUE)
  # A table without an interval has no level to carry.
  m <- DMAR:::.as_dmar_tbl(data.frame(term = c("mean", "sd"),
                                      value = c(1, 2)))
  expect_identical(names(generics::tidy(m)), c("term", "estimate"))
})

test_that("tidy() of a wide table takes as estimate the column its interval belongs to", {
  # The 1.1.0 review: the default tidy() took a wide table's first numeric
  # column as the estimate, so a row of simple_effects_AB() gave its F
  # statistic beside the interval for partial eta squared, and
  # results_sentence() on the tidy row wrote F as the estimate.
  fit <- aov(len ~ supp * factor(dose), data = ToothGrowth)
  se <- suppressWarnings(simple_effects_AB(fit))
  td <- generics::tidy(se)
  expect_equal(td$estimate, se$eta_squared_partial, tolerance = 0)
  expect_equal(td$F_value, se$F_value, tolerance = 0)
  expect_false("eta_squared_partial" %in% names(td))
  expect_identical(results_sentence(td[2, ]),
                   paste0("eta_squared_partial = 0.20, 95% CI [0.04, 0.35], ",
                          "unadjusted p = 0.0006"))
  expect_identical(results_sentence(td[2, ]), results_sentence(se[2, ]))
  # A table with an estimate column keeps every other column: the count
  # before the estimate of responder_analysis() was lost.
  w <- DMAR:::.as_dmar_tbl(data.frame(
    group = "a", n = 40, responders = 12, estimate = 0.3,
    ci_lower = 0.17, ci_upper = 0.46), conf_level = 0.95)
  tw <- generics::tidy(w)
  expect_identical(names(tw), c("term", "estimate", "n", "responders",
                                "ci_lower", "ci_upper", "conf_level"))
  expect_identical(tw$n, 40)
  # A table without an interval keeps its first numeric column.
  plain <- DMAR:::.as_dmar_tbl(data.frame(effect = "a", F_value = 3,
                                          df_1 = 1, epsilon = 0.8))
  expect_identical(generics::tidy(plain)$estimate, 3)
})

test_that("the sentence of a tidy() row names the quantity, not the row", {
  # The 1.1.0 review: a sentence from a row the default tidy() built from
  # a wide table was labeled with the row's term, which is the row number
  # when the table has no label column ("1 = -0.07"), a value of a label
  # column ("two-sample = 0.50"), or several pasted together. Through
  # 1.0.1 it read "estimate = ...". The tidy() now keeps the name of the
  # column the estimate came from in an estimand column, and the sentence
  # is the table's own.
  set.seed(113)
  counts <- t(apply(matrix(sample(1:3, 90, TRUE), 30), 1,
                    function(r) tabulate(r, 3)))
  fk <- fleiss_kappa(counts)
  tf <- generics::tidy(fk)
  expect_identical(tf$term, "1")
  expect_identical(tf$estimand, "kappa")
  expect_identical(names(tf)[1:3], c("term", "estimate", "estimand"))
  expect_identical(results_sentence(tf), results_sentence(fk))
  expect_match(results_sentence(tf), "^kappa = ")
  m <- ci_mahalanobis(D2 = .5, n_1 = 40, n_2 = 40, p = 3, conf_level = .9)
  tm <- generics::tidy(m)
  expect_identical(tm$term, "two-sample")
  expect_identical(results_sentence(tm), "D2 = 0.50, 90% CI [0.03, 1.03]")
  expect_identical(results_sentence(tm), results_sentence(m))
  # A row of a generic estimate column is still labeled by its term, which
  # names it (a group, a pair of factors, a coefficient).
  w <- DMAR:::.as_dmar_tbl(data.frame(
    group = "control", estimate = 0.3, ci_lower = 0.17,
    ci_upper = 0.46), conf_level = 0.95)
  tw <- generics::tidy(w)
  expect_false("estimand" %in% names(tw))
  expect_identical(results_sentence(tw),
                   "control = 0.30, 95% CI [0.17, 0.46]")
  # The column stays with its rows when tidy() results are stacked.
  st <- rbind(tf[, c("term", "estimate", "estimand", "ci_lower",
                     "ci_upper", "conf_level")],
              tm[, c("term", "estimate", "estimand", "ci_lower",
                     "ci_upper", "conf_level")])
  expect_match(results_sentence(st[2, ]), "^D2 = 0\\.50, 90% CI")
  expect_match(results_sentence(st[1, ]), "^kappa = ")
})

test_that("the glance() row of a long table gives the table's own sentence", {
  # The 1.1.0 review: glance() of a long table dropped the level and the
  # name of the estimate, so its sentence, which 1.0.1 wrote at an assumed
  # 95%, was refused, and with conf_level given it named the first column:
  # the observed correlation beside the interval of the corrected one.
  ca <- correction_for_attenuation(r = .3, reliability_x = .8,
                                   reliability_y = .7, N = 200,
                                   conf_level = .9)
  g <- generics::glance(ca)
  expect_identical(g$estimand, "correlation_corrected")
  expect_identical(g$conf_level, 0.9)
  expect_identical(names(g)[seq_len(nrow(ca))], ca$term)
  expect_identical(results_sentence(g),
                   "correlation_corrected = 0.40, 90% CI [0.25, 0.54]")
  expect_identical(results_sentence(g), results_sentence(ca))
  ps <- proportion_of_superiority(smd = .5, n_1 = 50, n_2 = 50,
                                  conf_level = .9)
  expect_identical(results_sentence(generics::glance(ps)),
                   results_sentence(ps))
  expect_match(results_sentence(ps), "^proportion_of_superiority = 0\\.69, ")
  p90 <- ci_proportion(30, 100, conf_level = .90)
  expect_identical(results_sentence(generics::glance(p90)),
                   "proportion = 0.30, 90% CI [0.23, 0.38]")
  expect_identical(results_sentence(generics::glance(p90)),
                   results_sentence(p90))
  # The TOST p-value of an equivalence table, as on the table, and the
  # p-value a table names through p_terms alone (contrast_adjusted()'s
  # row "p").
  e <- equivalence_r(r = .15, N = 300, cor_lower = .2, cor_upper = .2)
  expect_identical(results_sentence(generics::glance(e)),
                   results_sentence(e))
  set.seed(113)
  d <- data.frame(A = factor(rep(c("a1", "a2"), each = 40)),
                  B = factor(rep(rep(c("b1", "b2"), each = 20), 2)),
                  x = rnorm(80))
  d$y <- 5 + 2 * (d$A == "a2") + 1.5 * (d$B == "b2") + 0.8 * d$x + rnorm(80)
  ca <- contrast_adjusted(lm(y ~ A * B + x, data = d),
                          contrast = c(-0.5, 0.5, -0.5, 0.5))
  expect_match(results_sentence(ca), ", unadjusted p < 0\\.0001$")
  expect_identical(results_sentence(generics::glance(ca)),
                   results_sentence(ca))
  # The level and the estimate are columns, so they stay with their row
  # when glance() results of the same kind are stacked.
  p99 <- ci_proportion(30, 60, conf_level = .99)
  st <- rbind(generics::glance(p90), generics::glance(p99))
  expect_identical(results_sentence(st[2, ]), results_sentence(p99))
  expect_match(results_sentence(st[2, ]), "99% CI", fixed = TRUE)
  # A long table without an interval glances as before.
  m <- DMAR:::.as_dmar_tbl(data.frame(term = c("mean", "sd"),
                                      value = c(1, 2)))
  expect_identical(names(generics::glance(m)), c("mean", "sd"))
  # A table that reports its limits alone has no estimate: its glance()
  # stopped with "replacement has 0 rows" (1.0.1 returned the row of its
  # limits), and now names no estimand, so its sentence is refused as the
  # table's is.
  for (lim_only in list(
    ci_snr(F_value = 5, df_1 = 2, df_2 = 27, N = 30),
    ci_srsnr(F_value = 5, df_1 = 2, df_2 = 27, N = 30))) {
    g <- generics::glance(lim_only)
    expect_identical(names(g), c("ci_lower", "ci_upper", "estimand",
                                 "conf_level"))
    expect_identical(g$ci_lower,
                     lim_only$value[lim_only$term == "ci_lower"])
    expect_identical(g$estimand, NA_character_)
    expect_identical(g$conf_level, 0.95)
    expect_error(results_sentence(g), "point-estimate")
  }
})

test_that("two tables knitted in one chunk stay apart", {
  # The 1.1.0 review: knit_print() wrote no blank line before or after its
  # output, and knitr writes a chunk's outputs with nothing between them,
  # so the second table's header row ran onto the first table's last
  # footer line, which pandoc then read as a header cell.
  skip_if_not_installed("knitr")
  md <- knitr::knit(text = c("```{r, echo = FALSE}",
                             "ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)",
                             "ci_smd(smd = 0.3, n_1 = 50, n_2 = 50)",
                             "```"),
                    quiet = TRUE, envir = new.env())
  expect_false(grepl("95%|term", md, fixed = TRUE))
  expect_match(md, "Confidence level: 95%\n\n+\\|term")
  expect_identical(lengths(regmatches(md, gregexpr("\\|term", md))), 2L)
})

test_that("results_sentence reads the rmsea_ci_level row before the attribute, and the attribute for the RMSEA only", {
  # Two properties of the lookup that no table DMAR returns today can show,
  # because the row and the attribute agree in every one of them. A
  # hand-built long table in which they differ: the row is the level of the
  # interval in the rows beside it, so it comes first.
  tbl <- data.frame(
    term = c("rmsea", "rmsea_ci_lower", "rmsea_ci_upper", "rmsea_ci_level",
             "gamma_hat", "gamma_hat_ci_lower", "gamma_hat_ci_upper"),
    value = c(0.072, 0.030, 0.110, 0.99, 0.95, 0.92, 0.98),
    stringsAsFactors = FALSE)
  tbl <- DMAR:::.as_dmar_tbl(tbl, conf_level = 0.95)
  attr(tbl, "rmsea_conf_level") <- 0.90
  expect_match(results_sentence(tbl, term = "rmsea", digits = 3),
               "99% CI", fixed = TRUE)
  # Without the row, the attribute gives the level of the RMSEA interval.
  no_row <- tbl[tbl$term != "rmsea_ci_level", ]
  expect_match(results_sentence(no_row, term = "rmsea", digits = 3),
               "90% CI", fixed = TRUE)
  # The attribute is the level of the RMSEA interval only: another
  # estimate's interval in the same table is at the table's level.
  expect_match(results_sentence(no_row, term = "gamma_hat", digits = 3),
               "95% CI", fixed = TRUE)
})

test_that("results_sentence reads an interval held in columns named for its estimate", {
  # 1.1.0: ci_eta_squared() and ci_omega_squared() hold each estimand's
  # limits in <estimate>_ci_lower and <estimate>_ci_upper columns. A row
  # with one interval formed gives its sentence; a row with two needs
  # `term` to name one.
  w <- DMAR:::.as_dmar_tbl(data.frame(
    effect = "a", eta_squared_total = 0.3, eta_squared_total_ci_lower = NA,
    eta_squared_total_ci_upper = NA, eta_squared_partial = 0.4,
    eta_squared_partial_ci_lower = 0.1, eta_squared_partial_ci_upper = 0.6,
    F_value = 5, df_effect = 2, df_error = 30, N = 40),
    conf_level = 0.90)
  expect_identical(results_sentence(w),
                   "eta_squared_partial = 0.40, 90% CI [0.10, 0.60]")
  expect_identical(results_sentence(w, term = "eta_squared_partial",
                                    label = "partial eta squared"),
                   "partial eta squared = 0.40, 90% CI [0.10, 0.60]")
  expect_error(results_sentence(w, term = "eta_squared_total"),
               "limits are NA")
  w$eta_squared_total_ci_lower <- 0.05
  w$eta_squared_total_ci_upper <- 0.5
  expect_error(results_sentence(w), "intervals of several estimates",
               class = "dmar_argument_invalid")
  expect_identical(results_sentence(w, term = "eta_squared_total"),
                   "eta_squared_total = 0.30, 90% CI [0.05, 0.50]")
  expect_error(results_sentence(w, term = "omega_squared_total"),
               "names no estimate")
})
