# equivalence_proportion_paired(): the score tests of noninferiority,
# nonsuperiority, and equivalence for two paired proportions.

val <- function(res, term) res$value[res$term == term]

test_that("equivalence_proportion_paired() returns the documented rows", {
  res <- equivalence_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
                                       prop_diff_lower = 0.05,
                                       prop_diff_upper = 0.05)
  expect_s3_class(res, "dmar_tbl")
  # The flags close the table, as in equivalence_c(), with nonsuperior
  # last.
  expect_identical(res$term,
                   c("difference", "n_10", "n_01", "N",
                     "z_lower", "z_upper",
                     "p_lower", "p_upper", "p_tost", "p_noninferiority",
                     "ci_lower", "ci_upper",
                     "prop_diff_lower", "prop_diff_upper",
                     "equivalent", "noninferior", "superior", "inferior",
                     "nonsuperior"))
  expect_type(res$value, "double")
  flags <- c("equivalent", "noninferior", "superior", "inferior",
             "nonsuperior")
  expect_true(all(res$value[res$term %in% flags] %in% c(0, 1)))
  expect_equal(attr(res, "conf_level"), 0.90)
  expect_identical(attr(res, "p_sentence"), "p_tost")
  expect_identical(attr(res, "p_terms"),
                   c("p_lower", "p_upper", "p_tost", "p_noninferiority"))
  expect_true(attr(res, "verdict") %in%
                c("Equivalent", "Superior", "Inferior", "Noninferior only",
                  "Not relevantly superior only", "Inconclusive"))
})

test_that("Nam's (1997) noninferiority test of the liver scan data is reproduced", {
  # Nam (1997, p. 1427): n = 344 with 32 pairs positive by scan only and
  # 27 by pathology only; H0: psi = -.05 against psi > -.05 gives the
  # score statistic 2.795 with p = .0026.
  res <- equivalence_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
                                       prop_diff_lower = 0.05)
  expect_equal(round(val(res, "z_lower"), 3), 2.795)
  expect_equal(round(val(res, "p_lower"), 4), 0.0026)
  expect_identical(val(res, "p_noninferiority"), val(res, "p_lower"))
  expect_equal(val(res, "noninferior"), 1)
})

test_that("the decisions agree with the one-sided tests and the interval", {
  skip_on_cran()
  # Duality over a grid of tables, bounds, and levels: each test rejects
  # at alpha exactly when the 1 - 2 alpha score interval clears its bound.
  for (cs in list(c(32, 27, 344), c(38, 52, 400), c(12, 5, 100),
                  c(3, 0, 20), c(0, 7, 30), c(5, 5, 30))) {
    for (d in c(0.03, 0.05, 0.10, 0.20)) {
      for (a in c(0.025, 0.05, 0.10)) {
        res <- equivalence_proportion_paired(n_10 = cs[1], n_01 = cs[2],
                                             N = cs[3], prop_diff_lower = d,
                                             prop_diff_upper = d,
                                             alpha_level = a)
        ci <- ci_proportion_paired(n_10 = cs[1], n_01 = cs[2], N = cs[3],
                                   conf_level = 1 - 2 * a)
        lo <- ci$value[ci$term == "ci_lower"]
        hi <- ci$value[ci$term == "ci_upper"]
        expect_equal(c(val(res, "ci_lower"), val(res, "ci_upper")),
                     c(lo, hi))
        expect_identical(val(res, "p_lower") < a, lo > -d)
        expect_identical(val(res, "p_upper") < a, hi < d)
        expect_identical(val(res, "noninferior") == 1, lo > -d)
        expect_identical(val(res, "nonsuperior") == 1, hi < d)
        expect_identical(val(res, "equivalent") == 1,
                         val(res, "p_tost") < a)
      }
    }
  }
})

test_that("the TOST p-value is the larger one-sided p-value", {
  res <- equivalence_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
                                       prop_diff_lower = 0.05,
                                       prop_diff_upper = 0.05)
  expect_identical(val(res, "p_tost"),
                   max(val(res, "p_lower"), val(res, "p_upper")))
  expect_equal(val(res, "p_lower"),
               stats::pnorm(val(res, "z_lower"), lower.tail = FALSE))
  expect_equal(val(res, "p_upper"), stats::pnorm(val(res, "z_upper")))
  # The upper side's anchor, from Nam's (1997) equations 1 and 3 written
  # out here rather than from the package's functions: the restricted
  # estimate of p_01 at psi0 = .05 is the larger root of
  # 2n q^2 + {(2n + n_01 - n_10) psi0 - (n_10 + n_01)} q
  #   - n_01 psi0 (1 - psi0) = 0,
  # and Z(.05) for the liver scan data is -1.5625180.
  qa <- 2 * 344
  qb <- (2 * 344 + 27 - 32) * 0.05 - (32 + 27)
  qc <- -27 * 0.05 * (1 - 0.05)
  p01 <- (-qb + sqrt(qb^2 - 4 * qa * qc)) / (2 * qa)
  z <- (32 - 27 - 344 * 0.05) / sqrt(344 * (2 * p01 + 0.05 * (1 - 0.05)))
  expect_equal(z, -1.5625180323, tolerance = 1e-9)
  expect_equal(val(res, "z_upper"), z, tolerance = 1e-10)
})

test_that("each of the six verdicts is reached where it should be", {
  verdict <- function(n_10, n_01, n, d) {
    attr(equivalence_proportion_paired(n_10 = n_10, n_01 = n_01, N = n,
                                       prop_diff_lower = d, prop_diff_upper = d),
         "verdict")
  }
  # 90% limits [-.074, .004] against bounds of .10: inside.
  expect_identical(verdict(38, 52, 400, 0.10), "Equivalent")
  # The same limits against .05: below .05, but containing -.05.
  expect_identical(verdict(38, 52, 400, 0.05), "Not relevantly superior only")
  # 90% limits [-.023, .052] against .05: above -.05, containing .05.
  expect_identical(verdict(32, 27, 344, 0.05), "Noninferior only")
  # An estimate of .275 clears .05 entirely, and its mirror image falls
  # entirely below -.05.
  expect_identical(verdict(60, 5, 200, 0.05), "Superior")
  expect_identical(verdict(5, 60, 200, 0.05), "Inferior")
  # Thirty units with ten discordant pairs: the interval contains both
  # bounds.
  expect_identical(verdict(5, 5, 30, 0.05), "Inconclusive")
  # The flags behind the two new outcomes.
  ns <- equivalence_proportion_paired(n_10 = 38, n_01 = 52, N = 400,
                                      prop_diff_lower = 0.05, prop_diff_upper = 0.05)
  expect_equal(val(ns, "noninferior"), 0)
  expect_equal(val(ns, "nonsuperior"), 1)
  inc <- equivalence_proportion_paired(n_10 = 5, n_01 = 5, N = 30,
                                       prop_diff_lower = 0.05,
                                       prop_diff_upper = 0.05)
  expect_equal(inc$value[inc$term %in% c("equivalent", "noninferior",
                                         "superior", "inferior",
                                         "nonsuperior")], rep(0, 5))
})

test_that("prop_diff_upper alone makes the bounds symmetric", {
  e2 <- equivalence_proportion_paired(n_10 = 38, n_01 = 52, N = 400,
                                      prop_diff_lower = 0.10,
                                      prop_diff_upper = 0.10)
  e_up <- equivalence_proportion_paired(n_10 = 38, n_01 = 52, N = 400,
                                        prop_diff_upper = 0.10)
  expect_equal(e_up$value, e2$value)
  expect_identical(attr(e_up, "verdict"), attr(e2, "verdict"))
  expect_identical(attr(e_up, "p_sentence"), "p_tost")
  # Asymmetric bounds are stored as given, the lower one signed.
  e_as <- equivalence_proportion_paired(n_10 = 38, n_01 = 52, N = 400,
                                        prop_diff_lower = 0.08,
                                        prop_diff_upper = 0.03)
  expect_equal(val(e_as, "prop_diff_lower"), -0.08)
  expect_equal(val(e_as, "prop_diff_upper"), 0.03)
})

test_that("prop_diff_lower alone is a noninferiority test with no upper bound assumed", {
  # The liver scan question of Nam (1997, p. 1427) is noninferiority
  # alone. The table keeps its shape, the rows that would need an upper
  # bound are NA rather than read against a mirrored one, and every other
  # row agrees with a call that also gives an upper bound.
  ni <- equivalence_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
                                      prop_diff_lower = 0.05)
  both <- equivalence_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
                                        prop_diff_lower = 0.05,
                                        prop_diff_upper = 0.05)
  expect_identical(ni$term, both$term)
  expect_type(ni$value, "double")
  upper <- ni$term %in% c("z_upper", "p_upper", "p_tost", "prop_diff_upper",
                          "equivalent", "superior", "nonsuperior")
  expect_true(all(is.na(ni$value[upper])))
  expect_false(anyNA(ni$value[!upper]))
  expect_equal(ni$value[!upper], both$value[!upper])
  expect_identical(attr(ni, "verdict"), "Noninferior")
  expect_identical(attr(ni, "p_sentence"), "p_noninferiority")
  # Where the two-bound call says "Equivalent" against a mirrored bound,
  # the noninferiority call says only what was asked.
  expect_identical(attr(equivalence_proportion_paired(
    n_10 = 38, n_01 = 52, N = 400, prop_diff_lower = 0.10,
    prop_diff_upper = 0.10), "verdict"), "Equivalent")
  expect_identical(attr(equivalence_proportion_paired(
    n_10 = 38, n_01 = 52, N = 400, prop_diff_lower = 0.10), "verdict"),
    "Noninferior")
  # The other two outcomes of a noninferiority test.
  inf <- equivalence_proportion_paired(n_10 = 5, n_01 = 60, N = 200,
                                       prop_diff_lower = 0.05)
  expect_identical(attr(inf, "verdict"), "Inferior")
  expect_equal(val(inf, "inferior"), 1)
  expect_identical(attr(equivalence_proportion_paired(
    n_10 = 5, n_01 = 5, N = 30, prop_diff_lower = 0.05), "verdict"),
    "Inconclusive")
  # The table prints with the NA rows in place.
  expect_output(print(ni), "p_noninferiority")
})

test_that("raw vectors and counts give the same result", {
  set.seed(113)
  y <- stats::rbinom(200, 1, 0.70)
  x <- ifelse(stats::runif(200) < 0.15, 1 - y, y)
  a <- equivalence_proportion_paired(x = x, y = y, prop_diff_lower = 0.05)
  b <- equivalence_proportion_paired(n_10 = sum(x == 1 & y == 0),
                                     n_01 = sum(x == 0 & y == 1), N = 200,
                                     prop_diff_lower = 0.05)
  expect_equal(a$value, b$value)
  expect_identical(attr(a, "verdict"), attr(b, "verdict"))
})

test_that("results_sentence() reports the p-value of the question asked", {
  # Equivalence: the TOST p-value beside the 90% interval.
  res <- equivalence_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
                                       prop_diff_lower = 0.05,
                                       prop_diff_upper = 0.05)
  s <- results_sentence(res)
  expect_match(s, "^difference = ")
  expect_match(s, "90% CI", fixed = TRUE)
  expect_match(s, "TOST p", fixed = TRUE)
  # Noninferiority alone: its own p-value, never a TOST p-value against
  # a bound the user did not set.
  ni <- equivalence_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
                                      prop_diff_lower = 0.05)
  s <- results_sentence(ni)
  expect_match(s, "90% CI", fixed = TRUE)
  expect_match(s, "p_noninferiority = 0.0026", fixed = TRUE)
  expect_no_match(s, "TOST", fixed = TRUE)
})

test_that("equivalence_proportion_paired() refuses malformed input", {
  expect_error(equivalence_proportion_paired(n_10 = 1, n_01 = 2, N = 10),
               "Supply 'prop_diff_lower'")
  expect_error(equivalence_proportion_paired(n_10 = 1, n_01 = 2, N = 10,
                                             prop_diff_lower = 5),
               "in \\(0, 1\\)")
  expect_error(equivalence_proportion_paired(n_10 = 1, n_01 = 2, N = 10,
                                             prop_diff_upper = -0.1),
               "'prop_diff_upper' must be a single number in \\(0, 1\\)")
  expect_error(equivalence_proportion_paired(n_10 = 1, n_01 = 2, N = 10,
                                             prop_diff_lower = 0.1,
                                             alpha_level = 0.5),
               "alpha_level")
  expect_error(equivalence_proportion_paired(n_10 = 6, n_01 = 7, N = 10,
                                             prop_diff_lower = 0.1),
               "at least")
})
