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

test_that("results_sentence works on a single row of a wide table", {
  x <- ci_eta_squared(aov(len ~ supp * factor(dose), data = ToothGrowth))
  s <- results_sentence(x[1, ], label = "eta squared")
  expect_match(s, "^eta squared = 0\\.[0-9]{2}, 95% CI \\[0\\.[0-9]{2}, 0\\.[0-9]{2}\\]$")
})

test_that("results_sentence appends a p-value when the table carries one", {
  w <- data.frame(effect = "a", estimate = 0.5, lower_limit = 0.1,
                  upper_limit = 0.9, p_value = 1e-7)
  w <- DMAR:::.as_dmar_tbl(w, conf_level = 0.95)
  expect_equal(results_sentence(w),
               "estimate = 0.50, 95% CI [0.10, 0.90], p = < 0.0001")
})

test_that("results_sentence errors clearly when there is no interval", {
  d <- DMAR:::.as_dmar_tbl(data.frame(term = c("mean", "sd"),
                                      value = c(1, 2)))
  expect_error(results_sentence(d), "no confidence interval")
})
