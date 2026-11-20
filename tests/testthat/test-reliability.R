# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# Tests for the general reliability() wrapper: dispatch correctness,
# auto-detection of type, and consistency with the family functions.

skip_if_not_installed("lavaan")

set.seed(113)
J <- 6
N <- 250
loadings <- seq(0.4, 0.8, length.out = J)
eta <- rnorm(N)
errors <- matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
items_continuous <- sweep(matrix(rep(eta, J), N, J), 2, loadings, `*`) + errors
colnames(items_continuous) <- paste0("y", seq_len(J))
S_continuous <- cov(items_continuous)

items_binary <- (items_continuous > 0) * 1
items_categorical <- apply(items_continuous, 2, function(x)
  as.integer(cut(x, breaks = quantile(x, probs = seq(0, 1, length.out = 6)),
                 include.lowest = TRUE)))
colnames(items_categorical) <- colnames(items_continuous)


test_that("reliability() dispatches to reliability_alpha when type = 'alpha'", {
  r_wrap   <- reliability      (data = items_continuous, type = "alpha",
                                ci_method = "bonett")
  r_direct <- reliability_alpha(data = items_continuous,
                                ci_method = "bonett")
  expect_equal(r_wrap$value, r_direct$value)
  expect_equal(attr(r_wrap, "coefficient"), "alpha")
  # The transformation-scale standard error row and its scale label pass
  # through the wrapper unchanged.
  expect_true("se_transformed" %in% r_wrap$term)
  expect_identical(attr(r_wrap, "se_transform_scale"),
                   attr(r_direct, "se_transform_scale"))
})

test_that("reliability() dispatches to reliability_kr20 when type = 'kr20'", {
  r_wrap   <- reliability     (data = items_binary, type = "kr20")
  r_direct <- reliability_kr20(data = items_binary)
  expect_equal(r_wrap$value, r_direct$value)
  expect_equal(attr(r_wrap, "coefficient"), "kr20")
})

test_that("reliability() dispatches to reliability_omega when type = 'omega'", {
  r_wrap   <- reliability      (data = items_continuous, type = "omega",
                                denominator = "model_implied",
                                ci_method = "ml")
  r_direct <- reliability_omega(data = items_continuous,
                                denominator = "model_implied",
                                ci_method = "ml")
  expect_equal(r_wrap$value, r_direct$value)
  expect_equal(attr(r_wrap, "coefficient"), "omega")
})

test_that("reliability() forwards 'denominator' to reliability_omega()", {
  r_wrap   <- reliability      (data = items_continuous, type = "omega",
                                denominator = "observed",
                                ci_method = "none")
  r_direct <- reliability_omega(data = items_continuous,
                                denominator = "observed",
                                ci_method = "none")
  expect_equal(r_wrap$value, r_direct$value)
  expect_equal(attr(r_wrap, "coefficient"), "omega")
  expect_equal(attr(r_wrap, "denominator"), "observed")
})

test_that("reliability() rejects 'denominator' for non-omega types", {
  expect_error(
    reliability(data = items_continuous, type = "alpha",
                denominator = "observed"),
    regexp = "applies only to type"
  )
})

test_that("reliability() forwards 'denominator' to reliability_omega_categorical()", {
  r_wrap <- reliability(data = items_categorical, type = "omega_categorical",
                        denominator = "model_implied", ci_method = "none")
  r_direct <- reliability_omega_categorical(data = items_categorical,
                                            denominator = "model_implied",
                                            ci_method = "none")
  expect_equal(r_wrap$value, r_direct$value)
  expect_equal(attr(r_wrap, "denominator"), "model_implied")
  expect_error(
    reliability(data = items_binary, type = "kr20", denominator = "observed"),
    regexp = "applies only to type"
  )
})

test_that("reliability() rejects the retired type = 'omega_h'", {
  expect_error(
    reliability(data = items_continuous, type = "omega_h"),
    regexp = "arg"
  )
})

test_that("reliability() auto-detects robust omega for covariance matrix input", {
  msgs <- testthat::capture_messages(r <- reliability(S = S_continuous, N = N))
  expect_true(any(grepl("type = \"omega\"", msgs)))
  expect_true(any(grepl("covariance matrix", msgs)))
  expect_equal(attr(r, "coefficient"), "omega")
  expect_equal(attr(r, "denominator"), "observed")
  # Robust omega intervals are bootstrap based and need raw data, so the
  # covariance-input default is the point estimate.
  expect_equal(attr(r, "ci_method"), "none")
})

test_that("reliability() falls back to 'ml_logistic' for model implied omega with S only", {
  # The logit-scale form of the robust default, which needs raw data; DMAR
  # 1.0.x fell back to "ml" here.
  suppressMessages(
    r <- reliability(S = S_continuous, N = N, type = "omega",
                     denominator = "model_implied")
  )
  expect_equal(attr(r, "ci_method"), "ml_logistic")
})

test_that("reliability() auto-detects 'omega' for continuous data", {
  expect_message(
    r <- reliability(data = items_continuous, ci_method = "none"),
    regexp = "type = \"omega\""
  )
  expect_equal(attr(r, "coefficient"), "omega")
  expect_equal(attr(r, "denominator"), "observed")
})

test_that("reliability() auto-detects 'omega_categorical' for ordered categorical data", {
  expect_message(
    r <- reliability(data = items_categorical, ci_method = "none"),
    regexp = "type = \"omega_categorical\""
  )
  expect_equal(attr(r, "coefficient"), "omega_categorical")
})

test_that("reliability() accepts 'omega_c' as a type shorthand", {
  r_word  <- reliability(data = items_categorical,
                         type = "omega_categorical", ci_method = "none")
  r_short <- reliability(data = items_categorical,
                         type = "omega_c", ci_method = "none")
  expect_equal(r_word$value, r_short$value)
})

test_that("reliability() forwards ci_method to the chosen family function", {
  r <- reliability(data = items_continuous, type = "alpha",
                   ci_method = "feldt")
  expect_equal(attr(r, "ci_method"), "feldt")
})

test_that("reliability() respects the family-specific default when ci_method = NULL", {
  r <- reliability(data = items_continuous, type = "alpha")
  expect_equal(attr(r, "ci_method"), "bonett")
})

test_that("reliability() rejects S input for type = 'kr20' and 'omega_categorical'", {
  expect_error(reliability(S = S_continuous, N = N, type = "kr20"),
               regexp = "requires raw 'data'")
  expect_error(reliability(S = S_continuous, N = N, type = "omega_categorical"),
               regexp = "requires raw 'data'")
})

test_that("reliability() rejects unknown type", {
  expect_error(reliability(data = items_continuous, type = "lambda4"))
})

test_that("reliability() errors when neither data nor S is supplied (no type)", {
  expect_error(reliability(),
               regexp = "Either 'data' or 'S' must be supplied")
})

test_that("reliability() lets the two-item refusal of categorical omega through", {
  # Two integer items route to categorical omega, which is not identified
  # with two items; the wrapper reports that rather than substituting
  # another coefficient.
  skip_if_not_installed("lavaan")
  set.seed(113)
  two <- data.frame(a = sample(1:4, 100, replace = TRUE),
                    b = sample(1:4, 100, replace = TRUE))
  expect_error(suppressMessages(reliability(data = two)),
               "Categorical omega needs three or more items", fixed = TRUE)
  expect_error(reliability(data = two, type = "omega_categorical",
                           denominator = "model_implied"),
               "Categorical omega needs three or more items", fixed = TRUE)
})

test_that("reliability() forwards ncpus and improper, and refuses improper where it does not apply", {
  r <- reliability(data = items_continuous, type = "omega",
                   denominator = "model_implied", ci_method = "percentile",
                   B = 40, seed = 113, improper = "drop")
  expect_identical(attr(r, "bootstrap")$improper_rule, "drop")
  r_default <- reliability(data = items_continuous, type = "alpha",
                           ci_method = "percentile", B = 40, seed = 113)
  expect_identical(attr(r_default, "bootstrap")$improper_rule, "keep")
  expect_identical(attr(r_default, "bootstrap")$ncpus, 1L)
  expect_error(reliability(data = items_binary, type = "kr20",
                           improper = "drop"),
               "'improper' applies only to type = \"alpha\" and type = \"omega\": KR-20 fits no model")
  expect_error(reliability(data = items_categorical,
                           type = "omega_categorical", improper = "keep"),
               "categorical omega is not defined for an improper solution")
  expect_error(reliability(data = items_continuous, type = "alpha",
                           ncpus = -1),
               "'ncpus' must be a single whole number")
})

test_that("reliability() forwards categorical omega's missing data treatment, its imputation arguments, and B", {
  d <- as.data.frame(items_categorical)
  d$y5[1:40] <- NA
  d$y6[41:80] <- NA
  pw <- reliability(data = d, type = "omega_categorical",
                    missing = "pairwise", ci_method = "none")
  expect_identical(pw$value, reliability_omega_categorical(
    d, ci_method = "none", missing = "pairwise")$value)
  expect_identical(attr(pw, "missing"), "pairwise")
  expect_error(reliability(data = d, type = "omega", impute = "polr"),
               "apply only to type = \"omega_categorical\"", fixed = TRUE)
  expect_error(reliability(data = items_binary, type = "kr20",
                           missing = "pairwise"),
               "apply only to")
  seen <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(.omega_c_mi = function(...) {
    seen$args <- list(...)
    list(estimate = 0.8, se = 0.01, lower = 0.78, upper = 0.82,
         ci_method = "percentile (bootstrap then impute)", bootstrap = NULL,
         mi = list(), residual_cor = NULL,
         undetermined = matrix(character(0), ncol = 2L), notes = character())
  })
  r <- reliability(data = d, type = "omega_categorical", missing = "mi",
                   impute = "polr", imputations = 7, ci_method = "percentile")
  expect_identical(seen$args$impute, "polr")
  expect_identical(seen$args$imputations, 7L)
  # B left at NULL is the family function's default: 1000 under
  # multiple imputation.
  expect_identical(seen$args$B, 1000)
  expect_identical(attr(r, "B"), 1000L)
  reliability(data = d, type = "omega_categorical", missing = "mi",
              ci_method = "percentile", B = 50)
  expect_identical(seen$args$B, 50)
})

test_that("reliability() with B left at NULL keeps every family function's default of 10000", {
  seen <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    .bootstrap_ci = function(data, point_fn, B, ...) {
      seen$B <- c(seen$B, B)
      list(se = NA_real_, lower = NA_real_, upper = NA_real_,
           ci_method = "percentile", bootstrap = NULL)
    })
  expect_null(formals(reliability)$B)
  reliability(data = items_continuous, type = "alpha",
              ci_method = "percentile")
  reliability(data = items_binary, type = "kr20", ci_method = "percentile")
  reliability(data = items_continuous, type = "omega",
              ci_method = "percentile")
  reliability(data = items_categorical, type = "omega_categorical",
              ci_method = "percentile")
  expect_identical(seen$B, rep(10000, 4))
})

test_that("reliability() forwards wald_df to alpha and omega only", {
  expect_identical(
    reliability(data = items_continuous, type = "alpha", ci_method = "ml",
                wald_df = 20),
    reliability_alpha(data = items_continuous, ci_method = "ml",
                      wald_df = 20))
  expect_identical(
    reliability(S = S_continuous, N = N, type = "omega",
                denominator = "model_implied", ci_method = "ml",
                wald_df = 20),
    reliability_omega(S = S_continuous, N = N,
                      denominator = "model_implied", ci_method = "ml",
                      wald_df = 20))
  # The default leaves the family function's own default in place.
  expect_identical(
    reliability(data = items_continuous, type = "alpha", ci_method = "ml"),
    reliability_alpha(data = items_continuous, ci_method = "ml"))
  # KR-20, categorical omega, and H take no wald_df: a finite value is an
  # error there, which names the interval each builds instead (the error
  # said they kept "the normal reference", which neither Feldt's F
  # interval nor a bootstrap interval refers to), and Inf is accepted and
  # not forwarded.
  expect_error(reliability(data = items_binary, type = "kr20",
                           ci_method = "ml", wald_df = 20),
               paste0("applies only to type = \"alpha\" and type = ",
                      "\"omega\".*type = \"kr20\" takes no wald_df, and ",
                      "its interval is its own \\(Feldt's F interval by ",
                      "default\\)"))
  expect_error(reliability(data = items_categorical,
                           type = "omega_categorical", wald_df = 20),
               paste0("type = \"omega_categorical\" takes no wald_df, and ",
                      "its interval is its own \\(by default a ",
                      "bias-corrected and accelerated bootstrap interval\\)"))
  expect_error(reliability(data = items_continuous, type = "H",
                           wald_df = 20),
               paste0("type = \"H\" takes no wald_df, and its interval is ",
                      "its own \\(by default a Wald interval on the logit ",
                      "scale with the normal reference\\)"))
  expect_identical(
    reliability(data = items_binary, type = "kr20", ci_method = "ml",
                wald_df = Inf),
    reliability(data = items_binary, type = "kr20", ci_method = "ml"))
})

test_that("reliability() forwards wald_df = 'missing_information' to omega and says where it applies", {
  K <- "missing_information"
  expect_identical(
    reliability(S = S_continuous, N = N, type = "omega",
                denominator = "model_implied", wald_df = K),
    reliability_omega(S = S_continuous, N = N,
                      denominator = "model_implied", ci_method = "ml",
                      wald_df = K))
  # Alpha refuses it with the calibration's verdict, and the other
  # coefficients name the one it applies to.
  expect_error(reliability(data = items_continuous, type = "alpha",
                           ci_method = "ml", wald_df = K),
               "not available for reliability_alpha\\(\\): in the calibration")
  expect_error(reliability(data = items_binary, type = "kr20",
                           ci_method = "ml", wald_df = K),
               paste0("applies only to type = \"omega\" \\(with denominator ",
                      "= \"model_implied\"\\).*type = \"kr20\" takes no ",
                      "wald_df, and its interval is its own \\(Feldt's F"))
  expect_error(reliability(data = items_categorical,
                           type = "omega_categorical", wald_df = K),
               "type = \"omega_categorical\" takes no wald_df")
  expect_error(reliability(data = items_continuous, type = "H",
                           wald_df = K),
               "type = \"H\" takes no wald_df")
})

test_that("auto-detection reads the items, not the auxiliary variables", {
  # Through 1.0.1 every column of data entered the rule, so a continuous
  # auxiliary variable (or an integer one with more than 10 values) beside
  # ordered items sent them to continuous omega under FIML, and a factor
  # auxiliary stopped the call with a base R error.
  d <- as.data.frame(items_categorical)
  d$z <- eta + stats::rnorm(N)
  expect_identical(
    suppressMessages(DMAR:::.resolve_relia_type(NULL, d[1:J], NULL)),
    "omega_categorical")
  for (z in list(d$z, round(10 * d$z), factor(d$z > 0))) {
    d$z <- z
    expect_message(
      type <- DMAR:::.resolve_relia_type(NULL, d, NULL, aux = "z"),
      "type = \"omega_categorical\"")
    expect_identical(type, "omega_categorical")
  }
  # Without aux the auxiliary is read as an item, as the rule says.
  d$z <- eta
  expect_identical(suppressMessages(DMAR:::.resolve_relia_type(NULL, d, NULL)),
                   "omega")
})

test_that("reliability() with an auxiliary variable and ordered items computes categorical omega", {
  skip_on_cran()
  d <- as.data.frame(items_categorical[, 1:4])
  d$z <- eta + stats::rnorm(N)
  d[1:20, 2] <- NA
  r <- suppressMessages(reliability(data = d, aux = "z"))
  expect_identical(attr(r, "coefficient"), "omega_categorical")
  expect_identical(attr(r, "missing"), "mi")
})

test_that("SPSS user-defined missing values are missing throughout the family", {
  # haven::read_sav(user_na = TRUE) keeps a code declared missing as a
  # value of a haven_labelled_spss column (built here without haven).
  # Through 1.0.1 alpha, omega, and KR-20 analyzed it as a score: a 9 for
  # "refused" in a column of standardized scores, or beside 0/1 items.
  spss <- function(x, na = 9) {
    structure(as.double(ifelse(is.na(x), na, x)), labels = c(Refused = na),
              na_values = na,
              class = c("haven_labelled_spss", "haven_labelled",
                        "vctrs_vctr", "double"))
  }
  d <- as.data.frame(items_continuous)
  d$y2[1:30] <- NA
  sav <- d
  sav$y2 <- spss(d$y2)
  # Listwise by name: the test is about the codes, not the missing data
  # method.
  msgs <- testthat::capture_messages(
    a <- reliability_alpha(sav, missing = "listwise", ci_method = "none"))
  expect_match(msgs, "30 in y2", fixed = TRUE, all = FALSE)
  expect_identical(a$value, reliability_alpha(d, missing = "listwise",
                                              ci_method = "none")$value)
  o <- suppressMessages(reliability_omega(sav, missing = "listwise",
                                          ci_method = "none"))
  expect_identical(o$value, reliability_omega(d, missing = "listwise",
                                              ci_method = "none")$value)
  w <- suppressMessages(reliability(data = sav, type = "omega",
                                    missing = "listwise", ci_method = "none"))
  expect_identical(w$value, o$value)
  b <- as.data.frame(items_binary)
  b$y1[1:25] <- NA
  bsav <- b
  bsav$y1 <- spss(b$y1)
  k <- suppressMessages(reliability_kr20(bsav, ci_method = "none"))
  expect_identical(k$value, reliability_kr20(b, ci_method = "none")$value)
  # A column that declares nothing missing, or a matrix, is untouched.
  expect_identical(.relia_declared_missing(d), d)
  expect_identical(.relia_declared_missing(items_continuous), items_continuous)
})

test_that("SPSS user-defined missing values are missing under FIML and two-stage estimation, in the items and the auxiliary variables, and in a bootstrap", {
  skip_if_not_installed("lavaan")
  # Test of combinations no test reached: the declared values are set to NA
  # once, over every column, before the auxiliary variables are selected,
  # the cases observed on an auxiliary alone are counted, and the rule on
  # the number of cases is applied, so each call gives exactly what the
  # same data coded NA give.
  spss <- function(x, na_values = NULL, na_range = NULL, code) {
    structure(as.double(ifelse(is.na(x), code, x)),
              na_values = na_values, na_range = na_range,
              class = c("haven_labelled_spss", "haven_labelled",
                        "vctrs_vctr", "double"))
  }
  hs <- holzinger_swineford
  v <- c("t5_general_information", "t6_paragraph_comprehension",
         "t7_sentence", "t8_word_classification", "t9_word_meaning")
  d <- hs[v]
  names(d) <- paste0("v", 1:5)
  d$z <- hs$t1_visual_perception
  set.seed(1)
  d$v3[d$z > stats::median(d$z) & stats::runif(301) < .5] <- NA
  d$v5[sample(301, 35)] <- NA
  d$z[sample(301, 25)] <- NA
  d[c(7, 70, 170, 270), paste0("v", 1:5)] <- NA
  d$z[c(7, 70, 170, 270)] <- c(20, 25, 30, 35)
  sav <- d
  sav$v3 <- spss(d$v3, na_values = c(98, 99), code = 99)
  sav$v5 <- spss(d$v5, na_range = c(900, 999), code = 950)
  sav$z <- spss(d$z, na_values = -9, code = -9)
  same <- function(f) {
    a <- suppressWarnings(suppressMessages(f(sav)))
    b <- suppressWarnings(suppressMessages(f(d)))
    expect_identical(a$value, b$value)
  }
  same(function(x) reliability_alpha(x, aux = "z", missing = "fiml",
                                     ci_method = "ml"))
  same(function(x) reliability_alpha(x, aux = "z", missing = "two_stage",
                                     ci_method = "ml"))
  same(function(x) reliability_omega(x, aux = "z", missing = "fiml",
                                     denominator = "model_implied",
                                     ci_method = "ml"))
  same(function(x) reliability_omega(x, aux = "z", missing = "two_stage",
                                     denominator = "model_implied",
                                     ci_method = "ml"))
  same(function(x) reliability_H(data = x, aux = "z", missing = "fiml",
                                 ci_method = "ml"))
  same(function(x) reliability(data = x, aux = "z"))
  msgs <- testthat::capture_messages(suppressWarnings(
    reliability_alpha(sav, aux = "z", missing = "two_stage",
                      ci_method = "none")))
  expect_match(msgs, "in v3", fixed = TRUE, all = FALSE)
  expect_match(msgs, "in z", fixed = TRUE, all = FALSE)
  expect_match(msgs, "4 cases are observed on an auxiliary variable",
               fixed = TRUE, all = FALSE)
  skip_on_cran()
  same(function(x) reliability_omega(x, aux = "z", missing = "fiml",
                                     ci_method = "percentile", B = 20,
                                     seed = 1))
  # Declared missing values must not manufacture complete observations:
  # this selected small sample has a genuine affine-collapse obstacle,
  # including the auxiliary-only likelihood rows.
  small <- sav[c(1:6, 7, 70), ]
  for (m in c("fiml", "two_stage")) {
    expect_error(suppressMessages(reliability_alpha(small, aux = "z",
                                                    missing = m)),
                 class = "dmar_fiml_unbounded")
  }
})


test_that("SPSS user-defined missing values are missing in every column of a repeated name", {
  # The declared values were set to NA by column name, so with two
  # columns of one name both passes acted on the first, and the second
  # column's 9s were analyzed as scores (alpha 0.518 against 0.797).
  spss <- function(x, na = 9) {
    structure(as.double(ifelse(is.na(x), na, x)), na_values = na,
              class = c("haven_labelled_spss", "haven_labelled",
                        "vctrs_vctr", "double"))
  }
  d <- as.data.frame(items_continuous)
  names(d)[5] <- "y2"
  d[[2]][1:30] <- NA
  d[[5]][31:60] <- NA
  sav <- d
  sav[[2]] <- spss(d[[2]])
  sav[[5]] <- spss(d[[5]])
  msgs <- testthat::capture_messages(r <- .relia_declared_missing(sav))
  expect_identical(c(sum(is.na(r[[2]])), sum(is.na(r[[5]]))), c(30L, 30L))
  expect_false(inherits(r[[5]], "haven_labelled_spss"))
  expect_match(msgs, "30 in y2 (column 2), 30 in y2 (column 5)",
               fixed = TRUE)
  # Listwise deletion reads the columns by position, so the two y2 columns
  # are both analyzed (FIML selects the items by name and refuses them).
  a <- suppressMessages(reliability_alpha(sav, ci_method = "none",
                                          missing = "listwise"))
  expect_identical(a$value, reliability_alpha(d, ci_method = "none",
                                              missing = "listwise")$value)
})

test_that("SPSS user-defined missing values in holzinger_swineford give the values of NA", {
  # The NEWS example: five spatial tests, with the 156 missing scores of
  # t25_paper_form_board_r coded 999 and declared missing, as
  # haven::read_sav(user_na = TRUE) would return them. Version 1.0.1
  # analyzed 999 as a score; that is the analysis of the plain codes.
  d <- holzinger_swineford[c("t1_visual_perception", "t2_cubes",
                             "t3_paper_form_board", "t4_lozenges",
                             "t25_paper_form_board_r")]
  x <- d$t25_paper_form_board_r
  expect_identical(sum(is.na(x)), 156L)
  codes <- as.double(ifelse(is.na(x), 999, x))
  sav <- d
  sav$t25_paper_form_board_r <- structure(
    codes, na_values = 999,
    class = c("haven_labelled_spss", "haven_labelled", "vctrs_vctr",
              "double"))
  plain <- d
  plain$t25_paper_form_board_r <- codes
  # The NEWS values are those of listwise deletion, the default through
  # 1.0.x; the default analyzes the declared values as missing too.
  lw <- "listwise"
  a <- suppressMessages(reliability_alpha(sav, ci_method = "none",
                                          missing = lw))
  o <- suppressMessages(reliability_omega(sav, ci_method = "none",
                                          missing = lw))
  expect_equal(a$value[1], 0.6872627085, tolerance = 1e-9)
  expect_equal(o$value[1], 0.7296656975, tolerance = 1e-7)
  expect_identical(a$value, reliability_alpha(d, ci_method = "none",
                                              missing = lw)$value)
  expect_identical(o$value, reliability_omega(d, ci_method = "none",
                                              missing = lw)$value)
  quiet <- function(x) suppressWarnings(suppressMessages(x))
  expect_identical(quiet(reliability_omega(sav, ci_method = "none"))$value,
                   quiet(reliability_omega(d, ci_method = "none"))$value)
  # The values 1.0.1 returned for the declared codes.
  expect_equal(reliability_alpha(plain, ci_method = "none")$value[1],
               0.0083954724, tolerance = 1e-8)
  expect_equal(suppressWarnings(
    reliability_omega(plain, ci_method = "none"))$value[1],
    0.0146112183, tolerance = 1e-6)
})
