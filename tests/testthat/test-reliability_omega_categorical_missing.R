# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# Tests for the missing data treatments of reliability_omega_categorical():
# pairwise deletion and multiple imputation with a bootstrap then impute
# interval (missing = "pairwise", "mi"), added in 1.1.0. The fixture is
# the one of test-reliability_omega_categorical.R (five items, N = 400,
# loadings 0.7, four cut points) with a missing at random copy in which
# y4 and y5 go missing more often for respondents high on y1 to y3,
# which are always observed (130 and 119 missing; 219 complete rows).
# The pinned values were computed with lavaan 0.7-3 and mice 3.19.0.

skip_if_not_installed("lavaan")

set.seed(113)
N <- 400
J <- 5
loadings <- rep(0.7, J)
eta <- rnorm(N)
latent <- outer(eta, loadings) +
          matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
full <- as.data.frame(apply(latent, 2, function(x)
  as.integer(cut(x, breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf),
                 labels = FALSE))))
colnames(full) <- paste0("y", seq_len(J))
set.seed(2026)
z_mar <- as.numeric(scale(rowMeans(full[, 1:3])))
mar <- full
mar$y4[runif(N) < plogis(-1.2 + 1.5 * z_mar)] <- NA
mar$y5[runif(N) < plogis(-1.2 + 1.5 * z_mar)] <- NA

# The number of imputations of the observed data is at least 20; most
# tests below need neither that many nor von Hippel's rule, so they
# replace the floor by `imputations` itself and the rule by a fixed
# number. Mocks are registered on the frame of the test that calls this.
few_imputations <- function(rule = NA_integer_, env = parent.frame()) {
  testthat::local_mocked_bindings(
    .mi_imputations_floor = function(imputations) as.integer(imputations),
    .mi_imputations_rule = function(fmi_upper, cv = 0.05) rule,
    .package = "DMAR", .env = env)
}

est_of <- function(r) r$value[r$term == "estimate"]


test_that("the fixture has the pinned missingness", {
  expect_identical(sum(is.na(mar$y4)), 130L)
  expect_identical(sum(is.na(mar$y5)), 119L)
  expect_identical(sum(stats::complete.cases(mar)), 219L)
})

test_that("pairwise deletion of complete items is listwise deletion, to the last bit", {
  for (den in c("observed", "model_implied")) {
    expect_identical(.omega_c_from_data(full, den, missing = "pairwise"),
                     .omega_c_from_data(full, den))
  }
  expect_equal(.omega_c_from_data(full), 0.8106007095, tolerance = 1e-9)
})

test_that("with no item value missing, pairwise deletion and imputation return the complete data result", {
  ref <- reliability_omega_categorical(full, ci_method = "none")
  expect_message(
    r_pw <- reliability_omega_categorical(full, ci_method = "none",
                                          missing = "pairwise"),
    "No item value is missing")
  # No imputation runs: an imputation would fail here.
  testthat::local_mocked_bindings(
    .mi_complete = function(...) stop("an imputation ran"))
  expect_message(
    r_mi <- reliability_omega_categorical(full, ci_method = "none",
                                          missing = "mi"),
    "has nothing to do")
  for (r in list(r_pw, r_mi)) {
    expect_identical(r$value, ref$value)
    expect_identical(attr(r, "missing"), "listwise")
    expect_null(attr(r, "imputation", exact = TRUE))
  }
  # A row with no observed item is set aside first, so the rest decide.
  with_empty <- rbind(full, NA)
  expect_message(
    r_e <- reliability_omega_categorical(with_empty, ci_method = "none",
                                         missing = "pairwise"),
    "in the rows with an observed item")
  expect_identical(r_e$value, ref$value)
})

test_that("on the missing at random fixture listwise < pairwise < imputation < complete", {
  few_imputations()
  complete <- est_of(reliability_omega_categorical(full, ci_method = "none"))
  lw <- reliability_omega_categorical(mar, ci_method = "none")
  pw <- reliability_omega_categorical(mar, ci_method = "none",
                                      missing = "pairwise")
  expect_equal(est_of(lw), 0.7204162981, tolerance = 1e-8)
  expect_equal(est_of(pw), 0.7710249021, tolerance = 1e-8)
  expect_identical(attr(lw, "missing"), "listwise")
  expect_identical(attr(pw, "missing"), "pairwise")
  # N is every row with an observed item; N_complete the complete rows.
  expect_identical(pw$value[pw$term %in% c("N", "N_complete")], c(400, 219))
  expect_identical(lw$value[lw$term %in% c("N", "N_complete")], c(219, 219))
  mi <- suppressMessages(reliability_omega_categorical(
    mar, ci_method = "none", missing = "mi", imputations = 5, seed = 113))
  expect_identical(attr(mi, "missing"), "mi")
  expect_identical(mi$value[mi$term %in% c("N", "N_complete")], c(400, 219))
  expect_lt(est_of(lw), est_of(pw))
  expect_lt(est_of(pw), est_of(mi))
  expect_lt(est_of(mi), complete)
  # The version-free form of the prototype's finding: imputation lands
  # nearer the complete data value than pairwise deletion does.
  expect_lt(abs(est_of(mi) - complete), abs(est_of(pw) - complete))
  expect_lt(abs(est_of(mi) - complete), 0.03)
})

test_that("the imputations of the observed data reproduce the prototype (mice 3.19.0)", {
  skip_if_not_installed("mice")
  skip_if_not(packageVersion("mice") == "3.19.0")
  few_imputations()
  # The prototype called set.seed(113) and then mice(m = 5, maxit = 10);
  # the function's stream must match it draw for draw.
  r <- suppressMessages(reliability_omega_categorical(
    mar, ci_method = "none", missing = "mi", imputations = 5, seed = 113))
  pmm <- c(0.8204585985, 0.7858649497, 0.7964735150, 0.7907035809,
           0.7830097807)
  expect_equal(attr(r, "imputation")$estimates, pmm, tolerance = 1e-8)
  expect_equal(est_of(r), 0.7953020850, tolerance = 1e-8)
  expect_equal(est_of(r), mean(attr(r, "imputation")$estimates), tolerance = 1e-15)
  skip_on_cran()
  r <- suppressMessages(reliability_omega_categorical(
    mar, ci_method = "none", missing = "mi", impute = "polr",
    imputations = 5, seed = 113))
  polr <- c(0.8043311202, 0.8008970958, 0.8012855819, 0.8035216011,
            0.8068955979)
  expect_equal(attr(r, "imputation")$estimates, polr, tolerance = 1e-8)
  expect_equal(est_of(r), 0.8033861994, tolerance = 1e-8)
  expect_identical(attr(r, "imputation")$impute, "polr")
})

test_that("both denominators work under every treatment", {
  few_imputations()
  pw <- reliability_omega_categorical(mar, ci_method = "none",
                                      missing = "pairwise",
                                      denominator = "model_implied")
  expect_equal(est_of(pw), 0.7663908063, tolerance = 1e-8)
  expect_identical(attr(pw, "denominator"), "model_implied")
  for (den in c("observed", "model_implied")) {
    mi <- suppressMessages(reliability_omega_categorical(
      mar, ci_method = "none", missing = "mi", imputations = 3,
      denominator = den, seed = 113))
    expect_identical(attr(mi, "denominator"), den)
    expect_true(is.finite(est_of(mi)))
    # The residual polychoric correlations are averaged over imputations.
    expect_true(is.finite(attr(mi, "residual_cor")$largest))
  }
})

test_that("without a bootstrap the point estimate averages max(imputations, 20) imputations, and says so", {
  skip_if_not_installed("mice")
  skip_on_cran()
  expect_message(
    r <- reliability_omega_categorical(mar, missing = "mi", seed = 113),
    "request it with ci_method = \"percentile\"", fixed = TRUE)
  mi <- attr(r, "imputation")
  expect_identical(mi$imputations_point, 20L)
  expect_length(mi$estimates, 20L)
  expect_identical(mi$imputations_floor, 20L)
  expect_true(is.na(mi$imputations_rule))
  expect_true(is.na(mi$fmi))
  expect_null(mi$replicate_estimates)
  expect_identical(mi$maxit, 10L)
  expect_true(is.na(r$value[r$term == "se"]))
  expect_identical(attr(r, "ci_method"), "none")
  expect_null(attr(r, "B"))
  expect_message(
    reliability_omega_categorical(mar, missing = "mi", ci_method = "none",
                                  imputations = 25, seed = 113),
    "averages 25 imputations")
})

test_that("a seed reproduces an imputation and restores the caller's stream", {
  skip_if_not_installed("mice")
  few_imputations()
  run <- function() suppressMessages(reliability_omega_categorical(
    mar, ci_method = "none", missing = "mi", imputations = 3, seed = 113))
  expect_identical(run(), run())
  set.seed(1)
  x <- runif(1)
  set.seed(1)
  run()
  expect_identical(runif(1), x)
})

test_that("the refusals say what to do instead", {
  expect_error(reliability_omega_categorical(mar, missing = "mi",
                                             ci_method = "bca"),
               "use ci_method = \"percentile\"", fixed = TRUE)
  # The refusal does not depend on the data.
  expect_error(reliability_omega_categorical(full, missing = "mi",
                                             ci_method = "bca"),
               "has not been studied with bootstrap then impute")
  d <- cbind(mar, x = z_mar)
  for (m in c("listwise", "pairwise")) {
    expect_error(reliability_omega_categorical(d, missing = m, aux = "x"),
                 "use missing = \"mi\"", fixed = TRUE)
  }
  expect_error(reliability_omega_categorical(mar, impute = "polr"),
               "apply only to missing = \"mi\"", fixed = TRUE)
  expect_error(reliability_omega_categorical(mar, missing = "pairwise",
                                             imputations = 10),
               "apply only to missing = \"mi\"", fixed = TRUE)
  for (bad in list(1, 2.5, NA, c(3, 4), "5")) {
    expect_error(reliability_omega_categorical(mar, missing = "mi",
                                               imputations = bad),
                 "'imputations' must be a single whole number of at least 2")
  }
  expect_error(reliability_omega_categorical(mar, B = 0),
               "'B' must be NULL or a single whole number")
  expect_error(reliability_omega_categorical(d, missing = "mi",
                                             aux = "absent"),
               "not found in 'data'")
  # Non-integer items are refused under every treatment.
  bad <- mar
  bad$y1 <- bad$y1 + 0.5
  for (m in c("listwise", "pairwise", "mi")) {
    expect_error(reliability_omega_categorical(bad, missing = m),
                 "integer category values")
  }
  expect_error(reliability_omega_categorical(mar[, 1:2],
                                             missing = "pairwise"),
               "three or more items")
  # Two items never observed together have no pairwise polychoric.
  apart <- mar
  apart$y4[!is.na(apart$y5)] <- NA
  apart$y4[1:50] <- 1:5
  apart$y5[1:50] <- NA
  expect_error(reliability_omega_categorical(apart, missing = "pairwise"),
               "No row observes both y4 and y5")
  # An item no row observes is named before anything is fitted or
  # imputed.
  none <- mar
  none$y5 <- NA_integer_
  for (m in c("pairwise", "mi")) {
    expect_error(reliability_omega_categorical(none, missing = m),
                 "No row observes y5, so neither its thresholds")
  }
})

test_that("without mice, multiple imputation asks for it by name", {
  testthat::local_mocked_bindings(.mi_have_mice = function() FALSE)
  expect_error(reliability_omega_categorical(mar, missing = "mi"),
               "install.packages(\"mice\")", fixed = TRUE)
  # Complete items need no imputation, so no mice either.
  expect_message(reliability_omega_categorical(full, missing = "mi",
                                               ci_method = "none"),
                 "has nothing to do")
})

test_that("auxiliary variables enter the imputation model and nothing else", {
  skip_if_not_installed("mice")
  few_imputations()
  d <- cbind(mar, x = z_mar + rnorm(N, sd = 0.5))
  d$x[1:10] <- NA
  # aux implies missing = "mi", as it implies FIML in reliability_omega().
  r <- suppressMessages(reliability_omega_categorical(
    d, ci_method = "none", aux = "x", imputations = 3, seed = 113))
  expect_identical(attr(r, "missing"), "mi")
  expect_identical(attr(r, "aux"), "x")
  expect_identical(r$value[r$term %in% c("N", "N_complete", "n_items")],
                   c(400, 219, 5))
  expect_output(print(r), "given the items and the auxiliary variable x",
                fixed = TRUE)
  # mice receives the auxiliary column as a predictor.
  seen <- new.env(parent = emptyenv())
  orig <- .mi_complete
  testthat::local_mocked_bindings(.mi_complete = function(d, items, ...) {
    seen$cols <- names(d)
    seen$items <- items
    orig(d, items, ...)
  })
  suppressMessages(reliability_omega_categorical(
    d, ci_method = "none", aux = "x", imputations = 2, seed = 113))
  expect_identical(seen$cols, c(paste0("y", 1:5), "x"))
  expect_identical(seen$items, paste0("y", 1:5))
})

test_that("a row observed only on an auxiliary variable is imputed, not set aside", {
  skip_if_not_installed("mice")
  few_imputations()
  set.seed(7)
  d <- cbind(mar, x = z_mar + rnorm(N, sd = 0.5))
  items <- paste0("y", 1:5)
  d[1:12, items] <- NA
  d$x[12] <- NA                         # nothing observed on row 12
  seen <- new.env(parent = emptyenv())
  orig <- .mi_complete
  testthat::local_mocked_bindings(.mi_complete = function(d, items, ...) {
    seen$n <- nrow(d)
    orig(d, items, ...)
  })
  # The eleven rows with no observed item but an observed auxiliary are
  # imputed from it, as full information maximum likelihood elsewhere in
  # the family keeps them; the row with nothing observed is set aside.
  r <- suppressMessages(reliability_omega_categorical(
    d, ci_method = "none", aux = "x", imputations = 2, seed = 113))
  expect_identical(seen$n, 399L)
  expect_identical(r$value[r$term %in% c("N", "N_complete")],
                   c(399, sum(stats::complete.cases(d[items]))))
  expect_true(is.finite(est_of(r)))
  # Without an auxiliary, a row with no observed item carries nothing.
  rp <- suppressMessages(reliability_omega_categorical(
    d[items], missing = "pairwise", ci_method = "none"))
  expect_identical(rp$value[rp$term == "N"], 388)
})

test_that("the variance components recover a known fraction of missing information", {
  set.seed(9)
  B <- 2000
  m <- 5
  s2_boot <- 4e-4
  s2_imp <- 1e-4
  Tm <- 0.8 + matrix(rep(rnorm(B, 0, sqrt(s2_boot)), m), B, m) +
    matrix(rnorm(B * m, 0, sqrt(s2_imp)), B, m)
  vc <- .mi_variance_components(Tm)
  expect_equal(vc$s2_boot, s2_boot, tolerance = 0.1)
  expect_equal(vc$s2_imp, s2_imp, tolerance = 0.05)
  expect_equal(vc$fmi, 0.25, tolerance = 0.1)
  expect_gt(vc$fmi_upper, vc$fmi)
  expect_lt(vc$fmi_upper, 0.3)
  expect_equal(vc$sd_replicates, stats::sd(rowMeans(Tm)))
  # The one-way decomposition by hand.
  means <- rowMeans(Tm)
  msw <- mean(apply(Tm, 1, stats::var))
  expect_equal(vc$s2_imp, msw, tolerance = 1e-12)
  expect_equal(vc$s2_boot, (m * stats::var(means) - msw) / m,
               tolerance = 1e-12)
  rho_lower <- (m * stats::var(means) / msw /
                  stats::qf(0.975, B - 1, B * (m - 1)) - 1) / m
  expect_equal(vc$fmi_upper, 1 / rho_lower, tolerance = 1e-12)
  # Degenerate inputs.
  expect_true(is.na(.mi_variance_components(Tm[1, , drop = FALSE])$fmi))
  expect_identical(.mi_variance_components(matrix(0.8, 10, 3))$fmi, 0)
  # von Hippel's rule: 1 + (1/2)(FMI / .05)^2.
  expect_identical(.mi_imputations_rule(0.5), 51L)
  expect_identical(.mi_imputations_rule(1), 201L)
  expect_identical(.mi_imputations_rule(0.1), 3L)
  expect_true(is.na(.mi_imputations_rule(NA_real_)))
  expect_identical(.mi_imputations_floor(5L), 20L)
  expect_identical(.mi_imputations_floor(30L), 30L)
  # The imputations per replication that keep the widening at 5%.
  for (fmi in c(0.1, 0.39, 0.5, 0.8)) {
    k <- .mi_widening_imputations(fmi)
    expect_lte(sqrt(1 + fmi / k), 1.05 + 1e-12)
    if (k > 2L) expect_gt(sqrt(1 + fmi / (k - 1)), 1.05)
  }
})

test_that("a bootstrap then impute run brackets its estimate and records its pieces", {
  skip_if_not_installed("mice")
  skip_if_not_installed("boot")
  skip_on_cran()
  few_imputations(rule = 4L)
  w <- testthat::capture_warnings(
    r <- reliability_omega_categorical(mar, ci_method = "percentile",
                                       missing = "mi", B = 20,
                                       imputations = 2, seed = 113))
  est <- est_of(r)
  lo <- r$value[r$term == "ci_lower"]
  hi <- r$value[r$term == "ci_upper"]
  expect_true(lo < est && est < hi)
  expect_identical(attr(r, "ci_method"), "percentile (bootstrap then impute)")
  expect_identical(attr(r, "B"), 20L)
  bs <- attr(r, "bootstrap")
  expect_identical(bs$attempted, 20L)
  expect_identical(attr(r, "B_used"), bs$returned)
  expect_named(bs$reasons, c("error", "non_finite", "imputation"))
  # The engine's six levels, then the imputation failures of this loop.
  expect_identical(levels(bs$status), c("proper", "improper", "error",
                                        "non_finite", "undetermined",
                                        "not_converged", "imputation"))
  expect_identical(bs$improper_rule, "drop")
  expect_identical(bs$jackknife$attempted, 0L)
  mi <- attr(r, "imputation")
  expect_true(mi$fmi >= 0 && mi$fmi <= 1)
  expect_true(mi$fmi_ci_upper >= mi$fmi && mi$fmi_ci_upper <= 1)
  # The rule asked for 4 and the floor gave 2, so two more were made.
  expect_identical(mi$imputations_rule, 4L)
  expect_identical(mi$imputations_point, 4L)
  expect_length(mi$estimates, 4L)
  expect_equal(est, mean(mi$estimates), tolerance = 1e-15)
  # A replicate's value is the mean of its imputations.
  expect_identical(dim(mi$replicate_estimates), c(20L, 2L))
  used <- bs$status == "proper"
  expect_equal(bs$t[used], rowMeans(mi$replicate_estimates[used, ]),
               tolerance = 1e-15)
  # The percentile interval is that of the replicate means.
  # (boot warns that 20 replications put the limits at extreme order
  # statistics, as it does for the function's own call.)
  bci <- suppressWarnings(
    boot::boot.ci(list(R = 20), conf = 0.95, type = "perc", t0 = est,
                  t = ifelse(used, bs$t, NA)))
  expect_equal(c(lo, hi), bci$percent[4:5], tolerance = 1e-15)
  # The standard error is that of the point estimate's own M_0.
  expect_equal(r$value[r$term == "se"],
               sqrt(mi$s2_boot + mi$s2_imp / mi$imputations_point),
               tolerance = 1e-12)
  expect_equal(mi$widening, sqrt(1 + mi$fmi / 2), tolerance = 1e-15)
  out <- capture.output(print(r))
  expect_true(any(grepl("^Bootstrap: 20 replications attempted", out)))
  expect_true(any(grepl(paste0("^Missing data: multiple imputation by ",
                               "predictive mean matching"), out)))
  expect_true(any(grepl(paste0("^Imputations: 4 of the observed data ",
                               "\\(von Hippel's rule asks for 4, the floor ",
                               "is 2\\); 2 in each bootstrap replication"),
                        out)))
  # The Wald variants use the same standard error, centered at the
  # estimate.
  r_se <- suppressWarnings(reliability_omega_categorical(
    mar, ci_method = "bootstrap_se", missing = "mi", B = 20,
    imputations = 2, seed = 113))
  expect_identical(est_of(r_se), est)
  se <- r_se$value[r_se$term == "se"]
  expect_equal(r_se$value[r_se$term == "ci_upper"],
               est + stats::qnorm(0.975) * se, tolerance = 1e-12)
})

test_that("von Hippel's rule raises the imputations of the observed data after the bootstrap", {
  skip_if_not_installed("mice")
  skip_if_not_installed("boot")
  skip_on_cran()
  # The rule itself is checked above; here it records the upper limit it
  # is handed and asks for five imputations, so the run stays short
  # (with B = 12 the upper limit is often 1, for which the rule asks
  # for 201).
  seen <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    .mi_imputations_floor = function(imputations) as.integer(imputations),
    .mi_imputations_rule = function(fmi_upper, cv = 0.05) {
      seen$fmi_upper <- fmi_upper
      5L
    })
  r <- suppressWarnings(reliability_omega_categorical(
    mar, ci_method = "percentile", missing = "mi", B = 12,
    imputations = 2, seed = 7))
  mi <- attr(r, "imputation")
  expect_identical(seen$fmi_upper, mi$fmi_ci_upper)
  expect_identical(mi$imputations_rule, 5L)
  expect_identical(mi$imputations_point, 5L)
  expect_length(mi$estimates, 5L)
  # The first two imputations are those a call without the bootstrap
  # makes, so the bootstrap does not move them.
  r0 <- suppressMessages(reliability_omega_categorical(
    mar, ci_method = "none", missing = "mi", imputations = 2, seed = 7))
  expect_identical(mi$estimates[1:2], attr(r0, "imputation")$estimates)
})

test_that("a replication whose imputation fails is set aside and counted", {
  skip_if_not_installed("mice")
  skip_if_not_installed("boot")
  skip_on_cran()
  few_imputations()
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  orig <- .mi_complete
  testthat::local_mocked_bindings(.mi_complete = function(...) {
    calls$n <- calls$n + 1L
    # The first call imputes the observed data; calls 3 and 5 are the
    # second and fourth replications.
    if (calls$n %in% c(3L, 5L)) stop("the chains diverged")
    orig(...)
  })
  w <- testthat::capture_warnings(
    r <- reliability_omega_categorical(mar, ci_method = "percentile",
                                       missing = "mi", B = 8,
                                       imputations = 2, seed = 113))
  bs <- attr(r, "bootstrap")
  expect_identical(bs$reasons[["imputation"]], 2L)
  # An imputation failure is not a fit that did not converge (the
  # engine's code 6 must not catch it).
  expect_identical(bs$not_converged, 0L)
  expect_identical(bs$failed, 2L)
  expect_identical(bs$returned, 6L)
  expect_identical(as.character(bs$status[c(2, 4)]),
                   c("imputation", "imputation"))
  expect_identical(bs$imputation_messages, c("the chains diverged" = 2L))
  expect_true(all(is.na(attr(r, "imputation")$replicate_estimates[c(2, 4), ])))
  expect_match(w, "2 of 8 bootstrap replications were set aside (2 could not be imputed)",
               fixed = TRUE, all = FALSE)
  expect_match(w, "failed most often with: the chains diverged",
               fixed = TRUE, all = FALSE)
  expect_output(print(r), "(2 set aside: 2 could not be imputed)",
                fixed = TRUE)
})

test_that("bootstrap then impute reports no interval when fewer than half of the replications return a value", {
  skip_if_not_installed("mice")
  skip_if_not_installed("boot")
  skip_on_cran()
  few_imputations()
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  orig <- .mi_complete
  testthat::local_mocked_bindings(.mi_complete = function(...) {
    calls$n <- calls$n + 1L
    # Call 1 imputes the observed data; calls 3 to 7 are replications 2
    # to 6 of 8, so 3 return a value, fewer than half.
    if (calls$n %in% 3:7) stop("the chains diverged")
    orig(...)
  })
  w <- testthat::capture_warnings(
    r <- reliability_omega_categorical(mar, ci_method = "percentile",
                                       missing = "mi", B = 8,
                                       imputations = 2, seed = 113))
  bs <- attr(r, "bootstrap")
  expect_identical(bs$returned, 3L)
  expect_true(bs$require_half)
  expect_true(is.finite(est_of(r)))
  expect_true(all(is.na(r$value[r$term %in% c("se", "ci_lower",
                                               "ci_upper")])))
  expect_match(w, "fewer than half returned a value, so no interval is reported",
               fixed = TRUE, all = FALSE)
  expect_false(any(grepl("Caution:", w, fixed = TRUE)))
  expect_output(print(r), "fewer than half returned a value", fixed = TRUE)
})

test_that("a completed data set of the observed data without an estimate is set aside and counted", {
  skip_if_not_installed("mice")
  few_imputations()
  r0 <- suppressMessages(reliability_omega_categorical(
    mar, ci_method = "none", missing = "mi", imputations = 4, seed = 113))
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  orig <- .omega_c_parts
  testthat::local_mocked_bindings(.omega_c_parts = function(...) {
    calls$n <- calls$n + 1L
    p <- orig(...)
    if (calls$n == 3L) p$improper <- "y4"
    p
  })
  w <- testthat::capture_warnings(
    r <- suppressMessages(reliability_omega_categorical(
      mar, ci_method = "none", missing = "mi", imputations = 4,
      seed = 113)))
  mi <- attr(r, "imputation")
  expect_identical(mi$imputations_point, 4L)
  expect_identical(mi$imputations_returned, 3L)
  expect_identical(mi$imputations_set_aside[["improper"]], 1L)
  expect_identical(sum(mi$imputations_set_aside), 1L)
  expect_identical(mi$status, c("proper", "proper", "improper", "proper"))
  expect_true(is.na(mi$estimates[3L]))
  # The estimate averages the three that returned a value, and the
  # imputations themselves are those of the same call without the mock.
  expect_equal(est_of(r), mean(mi$estimates[-3L]), tolerance = 1e-15)
  expect_equal(mi$between_sd, stats::sd(mi$estimates[-3L]),
               tolerance = 1e-15)
  expect_identical(mi$estimates[-3L], attr(r0, "imputation")$estimates[-3L])
  expect_true(is.finite(attr(r, "residual_cor")$largest))
  expect_match(w, paste0("1 of the 4 completed data sets of the observed ",
                         "data was set aside (1 gave an improper (Heywood) ",
                         "solution), and the point estimate averages the 3"),
               fixed = TRUE, all = FALSE)
  expect_match(w, paste0("An improper set is one in which an error ",
                         "variance fell below zero (a standardized loading ",
                         "past one; here for y4)"), fixed = TRUE,
               all = FALSE)
  expect_output(print(r), paste0("Imputations: 4 of the observed data, 1 ",
                                 "of them set aside (1 gave an improper ",
                                 "(Heywood) solution)"), fixed = TRUE)
})

test_that("the call stops when fewer than half of the completed data sets of the observed data return a value", {
  skip_if_not_installed("mice")
  few_imputations()
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  orig <- .omega_c_parts
  testthat::local_mocked_bindings(.omega_c_parts = function(...) {
    calls$n <- calls$n + 1L
    if (calls$n %in% c(2L, 4L)) return(NULL)
    p <- orig(...)
    if (calls$n == 3L) p$improper <- "y4"
    p
  })
  expect_error(
    reliability_omega_categorical(mar, ci_method = "none", missing = "mi",
                                  imputations = 4, seed = 113),
    paste0("could not be computed on 3 of the 4 completed data sets of the ",
           "observed data (1 gave an improper (Heywood) solution, 2 did ",
           "not converge), more than half"),
    fixed = TRUE)
  # Exactly half is enough.
  calls$n <- 0L
  testthat::local_mocked_bindings(.omega_c_parts = function(...) {
    calls$n <- calls$n + 1L
    if (calls$n %in% c(2L, 4L)) return(NULL)
    orig(...)
  })
  r <- suppressWarnings(suppressMessages(reliability_omega_categorical(
    mar, ci_method = "none", missing = "mi", imputations = 4, seed = 113)))
  expect_identical(attr(r, "imputation")$imputations_returned, 2L)
  expect_identical(
    attr(r, "imputation")$imputations_set_aside[["non_converged"]], 2L)
})

test_that("a completed data set whose fit lavaan stops with an error is counted as an error, with its message", {
  skip_if_not_installed("mice")
  few_imputations()
  # .omega_c_parts() returned NULL for a lavaan error as for a fit that
  # did not converge, so such a set was counted as not converged and its
  # message was lost. It now carries the message, and the set is counted
  # under "error", whose first message the warning quotes.
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  orig <- .omega_c_parts
  testthat::local_mocked_bindings(.omega_c_parts = function(...) {
    calls$n <- calls$n + 1L
    if (calls$n == 2L) {
      return(list(converged = FALSE,
                  fit_error = "lavaan->lav_x():\n   something went wrong"))
    }
    if (calls$n == 4L) return(list(converged = FALSE))
    orig(...)
  })
  w <- testthat::capture_warnings(
    r <- suppressMessages(reliability_omega_categorical(
      mar, ci_method = "none", missing = "mi", imputations = 4,
      seed = 113)))
  mi <- attr(r, "imputation")
  expect_identical(mi$status, c("proper", "error", "proper", "non_converged"))
  expect_identical(mi$imputations_set_aside[["error"]], 1L)
  expect_identical(mi$imputations_set_aside[["non_converged"]], 1L)
  expect_match(w, "1 raised an error, 1 did not converge", fixed = TRUE,
               all = FALSE)
  expect_match(w, "The first error was: lavaan->lav_x(): something went wrong",
               fixed = TRUE, all = FALSE)
})

test_that("every completed data set of the observed data is fitted before the bootstrap, so none can discard it", {
  skip_if_not_installed("mice")
  skip_if_not_installed("boot")
  skip_on_cran()
  # The floor is 2 imputations and the rule's most is 4 (mocked), so the
  # observed data are fitted on calls 1 to 4 of .omega_c_parts() and the
  # bootstrap's replications come after.
  few_imputations(rule = 4L)
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  calls$at_boot <- NA_integer_
  orig_parts <- .omega_c_parts
  orig_boot <- .boot_mi
  testthat::local_mocked_bindings(
    .omega_c_parts = function(...) {
      calls$n <- calls$n + 1L
      p <- orig_parts(...)
      if (calls$n == 4L) p$improper <- "y5"
      p
    },
    .boot_mi = function(...) {
      calls$at_boot <- calls$n
      orig_boot(...)
    })
  w <- testthat::capture_warnings(
    r <- reliability_omega_categorical(mar, ci_method = "percentile",
                                       missing = "mi", B = 6,
                                       imputations = 2, seed = 113))
  expect_identical(calls$at_boot, 4L)
  mi <- attr(r, "imputation")
  # The rule asked for 4, so the reserve's improper set is among those
  # the estimate would average; it is set aside, and the finished
  # bootstrap is kept.
  expect_identical(mi$imputations_point, 4L)
  expect_identical(mi$imputations_returned, 3L)
  expect_identical(attr(r, "bootstrap")$attempted, 6L)
  expect_true(is.finite(r$value[r$term == "ci_lower"]))
  expect_equal(est_of(r), mean(mi$estimates[1:3]), tolerance = 1e-15)
  expect_equal(r$value[r$term == "se"],
               sqrt(mi$s2_boot + mi$s2_imp / 3), tolerance = 1e-12)
  expect_match(w, paste0("1 of the 4 completed data sets of the observed ",
                         "data was set aside"), fixed = TRUE, all = FALSE)
  # When more than half of all the completed data sets made fail, the
  # call stops before the bootstrap: floor 1 of 2 failed, all 3 of 4.
  calls$n <- 0L
  calls$at_boot <- NA_integer_
  testthat::local_mocked_bindings(.omega_c_parts = function(...) {
    calls$n <- calls$n + 1L
    p <- orig_parts(...)
    if (calls$n %in% 2:4) p$improper <- "y5"
    p
  })
  expect_error(
    reliability_omega_categorical(mar, ci_method = "percentile",
                                  missing = "mi", B = 6, imputations = 2,
                                  seed = 113),
    "could not be computed on 3 of the 4 completed data sets", fixed = TRUE)
  expect_true(is.na(calls$at_boot))
})

test_that("on the help page's example the improper completed data sets are set aside, not fatal (mice 3.19.0)", {
  skip_if_not_installed("mice")
  skip_if_not(packageVersion("mice") == "3.19.0")
  skip_on_cran()
  # The example of ?reliability_omega_categorical: Green and Yang's
  # mixed-skew design, with y3 and y4 missing at random given y1 + y2.
  set.seed(113)
  n <- 1000
  p_top <- c(.011, .027, .047, .118, .797)
  cut_top <- c(-Inf, qnorm(cumsum(p_top)[-5]), Inf)
  cut_bottom <- c(-Inf, qnorm(cumsum(rev(p_top))[-5]), Inf)
  latent <- outer(rnorm(n), rep(0.8, 4)) +
    matrix(rnorm(n * 4, sd = 0.6), n, 4)
  items <- cbind(y1 = cut(latent[, 1], cut_top, labels = FALSE),
                 y2 = cut(latent[, 2], cut_top, labels = FALSE),
                 y3 = cut(latent[, 3], cut_bottom, labels = FALSE),
                 y4 = cut(latent[, 4], cut_bottom, labels = FALSE))
  z <- as.numeric(scale(items[, "y1"] + items[, "y2"]))
  items_mar <- items
  items_mar[runif(n) < plogis(-1.5 - 1.2 * z), "y3"] <- NA
  items_mar[runif(n) < plogis(-1.5 - 1.2 * z), "y4"] <- NA
  # Neither the observed data nor the complete data give an improper fit.
  expect_length(.omega_c_parts(items_mar, missing = "pairwise")$improper, 0L)
  expect_length(.omega_c_parts(items)$improper, 0L)
  # Stopping on any improper completed data set, this call failed, as
  # it did for every seed from 1 to 5; 4 of its 20 are improper.
  w <- testthat::capture_warnings(
    r <- suppressMessages(reliability_omega_categorical(
      items_mar, denominator = "model_implied", missing = "mi", seed = 1)))
  mi <- attr(r, "imputation")
  expect_identical(mi$imputations_set_aside[["improper"]], 4L)
  expect_identical(mi$imputations_returned, 16L)
  expect_equal(est_of(r), 0.6961702045, tolerance = 1e-8)
  expect_match(w, "here for y3", fixed = TRUE, all = FALSE)
})

test_that("the seed gives the same bootstrap then impute run whatever ncpus is", {
  skip_if_not_installed("mice")
  skip_if_not_installed("boot")
  skip_on_cran()
  skip_on_os("windows")   # forked workers
  few_imputations(rule = 3L)
  run <- function(ncpus) suppressWarnings(reliability_omega_categorical(
    mar, ci_method = "percentile", missing = "mi", B = 6, imputations = 2,
    seed = 113, ncpus = ncpus))
  r1 <- run(1L)
  r2 <- run(2L)
  expect_identical(r1$value, r2$value)
  expect_identical(attr(r1, "imputation"), attr(r2, "imputation"))
  b1 <- attr(r1, "bootstrap")
  b2 <- attr(r2, "bootstrap")
  expect_identical(b1[setdiff(names(b1), "ncpus")],
                   b2[setdiff(names(b2), "ncpus")])
  expect_identical(b2$ncpus, 2L)
})

test_that("a pairwise bootstrap resamples the rows and refits by pairwise deletion", {
  skip_if_not_installed("boot")
  skip_on_cran()
  seen <- new.env(parent = emptyenv())
  seen$missing <- character()
  orig <- .omega_c_from_data
  testthat::local_mocked_bindings(.omega_c_from_data = function(data, denominator = "observed", missing = "listwise") {
    seen$missing <- c(seen$missing, missing)
    orig(data, denominator, missing)
  })
  r <- suppressWarnings(
    reliability_omega_categorical(mar, ci_method = "percentile",
                                  missing = "pairwise", B = 20, seed = 113))
  expect_true(all(seen$missing == "pairwise"))
  expect_identical(attr(r, "bootstrap")$attempted, 20L)
  est <- est_of(r)
  expect_true(r$value[r$term == "ci_lower"] < est &&
                est < r$value[r$term == "ci_upper"])
})

test_that("the printed lines name the treatment and its assumption, and listwise deletion prints as before", {
  pw <- reliability_omega_categorical(mar, ci_method = "none",
                                      missing = "pairwise")
  expect_output(print(pw), paste0("Missing data: pairwise deletion, ",
                                  "consistent only when values are missing ",
                                  "completely at random"), fixed = TRUE)
  lw <- reliability_omega_categorical(mar, ci_method = "none")
  expect_false(any(grepl("Missing data", capture.output(print(lw)))))
  expect_identical(.missing_treatment_footer_lines(lw), character())
  # Full information maximum likelihood, elsewhere in the family, gets
  # no line here either.
  x <- lw
  attr(x, "missing") <- "fiml"
  expect_identical(.missing_treatment_footer_lines(x), character())
  # The knitted table carries the same lines.
  expect_identical(.dmar_tbl_footer(pw)[["missing"]],
                   .missing_treatment_footer_lines(pw)[["missing"]])
  few_imputations()
  mi <- suppressMessages(reliability_omega_categorical(
    mar, ci_method = "none", missing = "mi", imputations = 2, seed = 113))
  lines <- .missing_treatment_footer_lines(mi)
  expect_named(lines, c("missing", "imputation"))
  expect_match(lines[["imputation"]],
               "Imputations: 2 of the observed data; the number von Hippel's (2020) rule asks for",
               fixed = TRUE)
  expect_false(grepl("bootstrap then impute", lines[["missing"]]))
})

test_that("an improper replication is described by its error variance, as the engine describes it", {
  # The bootstrap then impute warning said a loading "ran past one",
  # while the engine's own warning names the criterion, an error
  # variance below zero.
  counts <- list(failed = 3L, returned = 17L,
                 reasons = c(error = 0L, non_finite = 3L),
                 error_messages = integer(), undetermined = 0L,
                 improper = 3L, caution = FALSE,
                 imputation_messages = integer())
  notes <- .boot_mi_notes(counts, B = 20L, conf_level = 0.95,
                          imputations = 2L, fmi = NA_real_,
                          widening = NA_real_)
  expect_match(notes, paste0("An improper replication is one in which an ",
                             "error variance fell below zero (a ",
                             "standardized loading past one), so setting ",
                             "these aside can trim the upper tail"),
               fixed = TRUE, all = FALSE)
  expect_false(any(grepl("loading ran past one", notes, fixed = TRUE)))
})

test_that("the reason phrase names imputation failures apart from the others", {
  expect_identical(
    .bootstrap_reason_phrase(c(error = 1L, non_finite = 1L,
                               imputation = 2L)),
    "1 raised an error, 1 returned a non-finite value, 2 could not be imputed")
  # Every other result's reasons have no imputation element.
  expect_identical(.bootstrap_reason_phrase(c(error = 0L, non_finite = 2L)),
                   "2 returned a non-finite value")
})

test_that("labelled and factor columns are imputed as their plain codes", {
  skip_if_not_installed("mice")
  few_imputations()
  # A column read with haven::read_sav() is a double of class
  # haven_labelled carrying value labels; the structure below is that
  # object without loading haven. Under predictive mean matching mice
  # transposes residuals that inherit the class, which the vctrs method
  # behind haven refuses, so left to mice such items stopped the call.
  # An ordered factor whose labels are its codes reached mice as a factor
  # and was imputed by another model (polynomial contrasts for the other
  # items, optimal scaling for the item itself).
  labelled <- function(x, cls = "haven_labelled") {
    structure(as.double(x), labels = c(Never = 1, Always = 5),
              label = "An item", class = c(cls, "vctrs_vctr", "double"))
  }
  d <- cbind(mar, x = z_mar)
  plain <- .mi_plain_columns(d, paste0("y", 1:5))
  sav <- d
  for (j in names(sav)) sav[[j]] <- labelled(sav[[j]])
  sav$y2 <- labelled(d$y2, "haven_labelled_spss")
  fac <- d
  for (j in paste0("y", 1:5)) {
    fac[[j]] <- factor(fac[[j]], levels = 1:5, ordered = TRUE)
  }
  for (input in list(sav, fac)) {
    out <- .mi_plain_columns(input, paste0("y", 1:5))
    expect_identical(out, plain)
    expect_true(all(vapply(out[paste0("y", 1:5)], is.integer, logical(1L))))
    expect_true(is.double(out$x) && is.null(attributes(out$x)))
  }
  ref <- suppressMessages(reliability_omega_categorical(
    d, ci_method = "none", aux = "x", imputations = 2, seed = 1))
  for (input in list(sav, fac)) {
    r <- suppressMessages(reliability_omega_categorical(
      input, ci_method = "none", aux = "x", imputations = 2, seed = 1))
    expect_identical(r$value, ref$value)
    expect_identical(attr(r, "imputation")$estimates,
                     attr(ref, "imputation")$estimates)
  }
})

test_that("an item observed together with another on two rows is still imputed", {
  skip_if_not_installed("mice")
  few_imputations()
  # A split form design: y4 observed on rows 1 to 200, y5 on rows 199 to
  # 400, so the two share two rows and their pairwise correlation is
  # one. mice's collinearity check would stop imputing y5 (or y4), and
  # every fit would then use only the rows that observe it.
  set.seed(2)
  n <- 400
  e <- rnorm(n)
  s <- as.data.frame(sapply(1:6, function(j) {
    as.integer(cut(0.7 * e + rnorm(n, sd = sqrt(0.51)),
                   c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf)))
  }))
  names(s) <- paste0("y", 1:6)
  s$y4[201:400] <- NA
  s$y5[1:198] <- NA
  expect_identical(sum(!is.na(s$y4) & !is.na(s$y5)), 2L)
  set.seed(1)
  comp <- .mi_complete(s, names(s), "pmm", m = 2)
  expect_false(any(vapply(comp$data, anyNA, logical(1L))))
  r <- suppressMessages(reliability_omega_categorical(
    s, ci_method = "none", missing = "mi", imputations = 2, seed = 1))
  expect_true(is.finite(est_of(r)))
  # An item mice still leaves unimputed stops the imputation, by name.
  cst <- s
  cst$y5 <- ifelse(is.na(s$y5), NA_integer_, 2L)
  expect_error(.mi_complete(cst, names(cst), "pmm", m = 2),
               "mice did not impute y5 (set aside as constant)",
               fixed = TRUE)
  # The function refuses such an item before imputing it.
  expect_error(reliability_omega_categorical(cst, missing = "mi"),
               paste0("y5 takes a single observed value, so its ",
                      "thresholds cannot be estimated, and imputation ",
                      "could only repeat it"), fixed = TRUE)
})

test_that("complete items under missing = \"mi\" say which B their bootstrap uses", {
  seen <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(.bootstrap_ci = function(data, point_fn,
                                                           B, ...) {
    seen$B <- B
    list(se = NA_real_, lower = NA_real_, upper = NA_real_,
         ci_method = "percentile", bootstrap = NULL)
  })
  expect_message(
    r <- reliability_omega_categorical(full, missing = "mi",
                                       ci_method = "percentile"),
    paste0("(missing = \"listwise\"), whose bootstrap has B = 10000 ",
           "replications (its default; each is one refit, with nothing ",
           "to impute)"), fixed = TRUE)
  expect_identical(seen$B, 10000)
  expect_message(
    reliability_omega_categorical(full, missing = "mi",
                                  ci_method = "percentile", B = 500),
    "whose bootstrap has B = 500 replications.", fixed = TRUE)
  expect_identical(seen$B, 500)
  # Without a bootstrap, B is not mentioned.
  expect_message(
    reliability_omega_categorical(full, missing = "mi", ci_method = "none"),
    "(missing = \"listwise\").", fixed = TRUE)
})


# Identification of the loadings under multiple imputation --------------

mixed_skew_mar <- function(N, seed) {
  set.seed(seed)
  p_top <- c(.011, .027, .047, .118, .797)
  cut_top <- c(-Inf, stats::qnorm(cumsum(p_top)[-5]), Inf)
  cut_bottom <- c(-Inf, stats::qnorm(cumsum(rev(p_top))[-5]), Inf)
  latent <- outer(stats::rnorm(N), rep(0.8, 4)) +
    matrix(stats::rnorm(N * 4, sd = 0.6), N, 4)
  items <- cbind(y1 = cut(latent[, 1], cut_top, labels = FALSE),
                 y2 = cut(latent[, 2], cut_top, labels = FALSE),
                 y3 = cut(latent[, 3], cut_bottom, labels = FALSE),
                 y4 = cut(latent[, 4], cut_bottom, labels = FALSE))
  z <- as.numeric(scale(items[, "y1"] + items[, "y2"]))
  items_mar <- items
  items_mar[stats::runif(N) < stats::plogis(-1.5 - 1.2 * z), "y3"] <- NA
  items_mar[stats::runif(N) < stats::plogis(-1.5 - 1.2 * z), "y4"] <- NA
  list(items = items, items_mar = items_mar)
}

test_that("completed data sets whose determined correlations do not identify the loadings are set aside, not averaged", {
  skip_if_not_installed("mice")
  skip_on_cran()
  # The seed 10 sample of the help page's design at N = 300, with the
  # help page's missing at random mechanism. In most completed data sets
  # every correlation of y1 or y2 with y3 or y4 is undetermined, so the
  # loadings are not identified; they were counted as proper under the
  # model implied denominator, and 12 of the 17 averaged sat between 0.04
  # and 0.07, which made the estimate 0.243 (population value 0.716).
  d <- mixed_skew_mar(300, 10)$items_mar
  err <- tryCatch(suppressMessages(reliability_omega_categorical(
    d, denominator = "model_implied", missing = "mi", seed = 1,
    ci_method = "none")), error = function(e) conditionMessage(e))
  expect_type(err, "character")
  expect_match(err, "left a polychoric correlation undetermined",
               fixed = TRUE)
  expect_match(err, "leave the single-factor loadings unidentified",
               fixed = TRUE)
  expect_match(err, "links {y1, y2} to {y3, y4}", fixed = TRUE)
  # Under the observed denominator the same imputations leave one set
  # more undetermined, one whose loadings are identified. The stop
  # advised the model implied denominator for it, although under that
  # denominator 19 of the 20 sets still return no value (the stop
  # above); it now says that the switch is not enough.
  err_obs <- tryCatch(suppressMessages(reliability_omega_categorical(
    d, denominator = "observed", missing = "mi", seed = 1,
    ci_method = "none")), error = function(e) conditionMessage(e))
  expect_match(err_obs, "20 of the 20 completed data sets", fixed = TRUE)
  expect_match(err_obs, "on only 1 of the 20 in all, still fewer than half",
               fixed = TRUE)
  expect_false(grepl('Use denominator = "model_implied"', err_obs,
                     fixed = TRUE))
  # Where the switch is enough, it is offered (a record built by hand:
  # 8 proper, 4 identified but undetermined, 8 improper of 20).
  pt <- list(status = rep(c("proper", "undetermined", "improper"),
                          c(8L, 4L, 8L)),
             improper = rep(list("y4"), 20L),
             undetermined = rep(list(matrix(c("y1", "y3"), 1L)), 20L),
             unidentified = logical(20L), messages = rep(NA_character_, 20L),
             items = c("y1", "y2", "y3", "y4"))
  adv <- .omega_c_mi_advice(pt, .omega_c_mi_tally(pt$status))
  expect_match(adv, 'Use denominator = "model_implied"', fixed = TRUE,
               all = FALSE)
  expect_match(adv, "12 of the 20", fixed = TRUE, all = FALSE)
  pt$status[7:10] <- "improper"
  adv <- .omega_c_mi_advice(pt, .omega_c_mi_tally(pt$status))
  expect_false(any(grepl("Use denominator", adv, fixed = TRUE)))
  expect_match(adv, "on only 8 of the 20 in all", fixed = TRUE,
               all = FALSE)
  # The bookkeeping, set by set, with a mocked fit: a set marked
  # unidentified is "undetermined" under either denominator, before the
  # improper verdict, and the notes say why.
  few_imputations()
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  orig <- .omega_c_parts
  testthat::local_mocked_bindings(.omega_c_parts = function(...) {
    calls$n <- calls$n + 1L
    p <- orig(...)
    if (calls$n == 2L) {
      p$identified <- FALSE
      p$improper <- "y4"
      p$undetermined <- matrix(c("y1", "y4", "y2", "y5"), 2L, byrow = TRUE)
    }
    p
  })
  w <- testthat::capture_warnings(r <- suppressMessages(
    reliability_omega_categorical(mar, ci_method = "none", missing = "mi",
                                  imputations = 4, seed = 113,
                                  denominator = "model_implied")))
  mi <- attr(r, "imputation")
  expect_identical(mi$status[2L], "undetermined")
  expect_identical(mi$imputations_set_aside[["undetermined"]], 1L)
  expect_identical(mi$imputations_set_aside[["improper"]], 0L)
  expect_match(w, "do not identify the single-factor loadings",
               fixed = TRUE, all = FALSE)
  expect_false(any(grepl("keeps every completed data set", w, fixed = TRUE)))
})

test_that("a replication averages the completed data sets that return a value, as the point estimate does", {
  # Through the first versions of 1.1.0 a replication was set aside when
  # any one of its completed data sets returned no value, and its other
  # sets were discarded.
  testthat::local_mocked_bindings(.mi_complete = function(d, items, impute,
                                                          m, maxit) {
    list(data = lapply(seq_len(m), function(k) data.frame(k = k)),
         polr_fallbacks = 0L)
  })
  est_fn <- function(dk) {
    if (dk$k %in% c(2L, 4L)) {
      structure(NA_real_, dmar_reason = "improper")
    } else {
      0.7 + dk$k / 100
    }
  }
  # Three of five return a value: the replication is their mean.
  r <- .mi_replicate(data.frame(), "y", est_fn, "pmm", 5L, 1L)
  expect_identical(r$row[2L], 1)
  expect_equal(r$row[1L], mean(c(0.71, 0.73, 0.75)), tolerance = 1e-15)
  expect_equal(r$row[3L + 1:5], c(0.71, NA, 0.73, NA, 0.75))
  expect_identical(r$row[3L + 5L + 1:5], c(1, 2, 1, 2, 1))
  # Two of four is half, still enough, as for the point estimate.
  r <- .mi_replicate(data.frame(), "y", est_fn, "pmm", 4L, 1L)
  expect_identical(r$row[2L], 1)
  expect_equal(r$row[1L], mean(c(0.71, 0.73)), tolerance = 1e-15)
  # Fewer than half: set aside, under the first set's reason.
  est_fn2 <- function(dk) {
    if (dk$k == 1L) 0.7 else structure(NA_real_, dmar_reason = "improper")
  }
  r <- .mi_replicate(data.frame(), "y", est_fn2, "pmm", 3L, 1L)
  expect_true(is.na(r$row[1L]))
  expect_identical(r$row[2L], 2)
  # The tally of the sets set aside in the replications kept.
  codes <- rbind(c(1, 2, 1, 5), c(1, 1, 3, 1))
  expect_identical(.boot_mi_sets_set_aside(codes),
                   c(improper = 1L, undetermined = 1L, error = 1L,
                     non_converged = 0L, non_finite = 0L))
})

test_that("the share of replications set aside falls with imputations of one parity and rises from an even number to the next odd one", {
  # With a share p of completed data sets returning no value, the old
  # rule set a replication aside with probability 1 - (1 - p)^M_b: at
  # p = .2, 36% at M_b = 2 and 89% at M_b = 10. Under the point
  # estimate's rule it is the chance that fewer than half return a
  # value, P(Bin(M_b, 0.8) < M_b / 2): 4.0, 10.4, 2.7, 5.8% at M_b = 2 to
  # 5 and 0.6% at 10. An odd M_b needs one more set to return a value
  # than the even number below it, so it sets more aside; the texts
  # first said that the share falls as M_b grows.
  testthat::local_mocked_bindings(.mi_complete = function(d, items, impute,
                                                          m, maxit) {
    list(data = lapply(seq_len(m), function(k) data.frame(k = k)),
         polr_fallbacks = 0L)
  })
  est_fn <- function(dk) {
    if (stats::runif(1L) < 0.2) {
      structure(NA_real_, dmar_reason = "improper")
    } else {
      0.7 + stats::rnorm(1L, sd = 0.01)
    }
  }
  share <- function(m) {
    set.seed(113)
    b <- .boot_mi(data.frame(y = 1:30), "y", est_fn, "pmm", m, 1L, 2000L, 1L)
    mean(!is.finite(b$t[, 1L]))
  }
  s <- vapply(c(2L, 3L, 4L, 5L, 10L), share, numeric(1L))
  exact <- stats::pbinom(ceiling(c(2, 3, 4, 5, 10) / 2) - 1,
                         c(2, 3, 4, 5, 10), 0.8)
  # Within three Monte Carlo standard errors at B = 2000 (0.004 at 4%).
  expect_lt(max(abs(s - exact)), 0.015)
  expect_gt(s[2L], s[1L])  # 2 to 3 rises
  expect_gt(s[4L], s[3L])  # 4 to 5 rises
  expect_lt(s[3L], s[1L])  # 2 to 4 falls
  expect_lt(s[4L], s[2L])  # 3 to 5 falls
  expect_lt(s[5L], s[1L])
  expect_lt(s[5L], 0.03)
  # The number of imputations the widening warning names is raised to an
  # even number when completed data sets were set aside, never otherwise.
  expect_identical(.mi_widening_imputations(0.5, kept = 1), 5L)
  expect_identical(.mi_widening_imputations(0.5, kept = 1, even = TRUE), 6L)
  expect_identical(.mi_widening_imputations(0.8, kept = 0.8, even = TRUE),
                   10L)
  expect_identical(.mi_widening_imputations(1, kept = 0.8, even = TRUE), 14L)
})

test_that("the variance components allow for replications that average fewer imputations", {
  set.seed(4)
  B <- 400
  m <- 5
  Tm <- 0.8 + matrix(rep(rnorm(B, 0, 0.02), m), B, m) +
    matrix(rnorm(B * m, 0, 0.01), B, m)
  full_vc <- .mi_variance_components(Tm)
  # Balanced: the formulas of the balanced design, unchanged.
  means <- rowMeans(Tm)
  expect_equal(full_vc$s2_imp, mean(apply(Tm, 1, stats::var)),
               tolerance = 1e-12)
  expect_equal(full_vc$inv_n, 1 / m)
  # Set aside one or two sets in a third of the replications.
  Tu <- Tm
  Tu[1:130, 2] <- NA
  Tu[1:60, 4] <- NA
  vc <- .mi_variance_components(Tu)
  n_b <- rowSums(!is.na(Tu))
  N <- sum(n_b)
  mb <- rowSums(Tu, na.rm = TRUE) / n_b
  grand <- sum(n_b * mb) / N
  msb <- sum(n_b * (mb - grand)^2) / (B - 1)
  msw <- sum((Tu - mb)^2, na.rm = TRUE) / (N - B)
  n0 <- (N - sum(n_b^2) / N) / (B - 1)
  expect_equal(vc$s2_imp, msw, tolerance = 1e-12)
  expect_equal(vc$s2_boot, (msb - msw) / n0, tolerance = 1e-12)
  expect_equal(vc$fmi, min(msw / ((msb - msw) / n0), 1), tolerance = 1e-12)
  rho_lower <- (msb / msw / stats::qf(0.975, B - 1, N - B) - 1) / n0
  expect_equal(vc$fmi_upper, min(1 / rho_lower, 1), tolerance = 1e-12)
  expect_equal(vc$inv_n, mean(1 / n_b))
  expect_equal(vc$sd_replicates, stats::sd(mb))
  # Both recover the components the data were drawn with.
  expect_equal(vc$s2_imp, 1e-4, tolerance = 0.1)
  expect_equal(vc$s2_boot, 4e-4, tolerance = 0.15)
  # No replication with two sets: nothing to estimate the within part.
  expect_true(is.na(.mi_variance_components(cbind(Tm[, 1], NA))$fmi))
  # The imputations named for a 5% widening allow for the sets set aside.
  expect_identical(.mi_widening_imputations(0.5, kept = 0.5),
                   .mi_widening_imputations(1))
})

test_that("the warning of a bootstrap then impute run counts the sets set aside in the replications kept", {
  counts <- list(failed = 2L, returned = 18L,
                 reasons = c(error = 0L, non_finite = 2L, imputation = 0L),
                 error_messages = integer(), undetermined = 0L,
                 improper = 2L, caution = FALSE,
                 imputation_messages = integer())
  sets <- c(improper = 5L, undetermined = 1L, error = 0L,
            non_converged = 0L, non_finite = 0L)
  notes <- .boot_mi_notes(counts, B = 20L, conf_level = 0.95,
                          imputations = 5L, fmi = NA_real_,
                          widening = NA_real_, denominator = "model_implied",
                          sets_aside = sets)
  expect_match(notes, paste0("a replication averages those of its 5 ",
                             "completed data sets that return a value, as ",
                             "the point estimate does, and is set aside ",
                             "when fewer than half do"),
               fixed = TRUE, all = FALSE)
  expect_match(notes, paste0("In the 18 replications that returned a value, ",
                             "6 of their 90 completed data sets returned none ",
                             "(5 gave an improper (Heywood) solution, 1 left a ",
                             "polychoric correlation undetermined)"),
               fixed = TRUE, all = FALSE)
  expect_match(notes, "lowers its replication's mean", fixed = TRUE,
               all = FALSE)
  expect_match(notes, "identify them only when", fixed = TRUE, all = FALSE)
  expect_false(any(grepl("keeps every", notes, fixed = TRUE)))
  # The number of imputations named for the widening is even when sets
  # were set aside (13 raised to 14 here), and as computed otherwise.
  notes <- .boot_mi_notes(counts, B = 20L, conf_level = 0.95,
                          imputations = 5L, fmi = 1, widening = 1.08,
                          denominator = "model_implied", sets_aside = sets,
                          kept = 0.8)
  expect_match(notes, "imputations = 14 would bring the widening under 5% (an even number",
               fixed = TRUE, all = FALSE)
  clean <- counts
  clean$failed <- 0L
  notes <- .boot_mi_notes(clean, B = 20L, conf_level = 0.95,
                          imputations = 5L, fmi = 1, widening = 1.08,
                          denominator = "model_implied",
                          sets_aside = sets * 0L, kept = 1)
  expect_match(notes, "imputations = 10 would bring the widening under 5%",
               fixed = TRUE, all = FALSE)
  expect_false(any(grepl("an even number", notes, fixed = TRUE)))
  # Regression test. Replications set aside only because their imputation
  # failed have no completed data set, and the half rule behind the even
  # number never applies to them: the number was raised all the same (5
  # to 6 here), with the half rule given as the reason.
  imp_only <- list(failed = 3L, returned = 17L,
                   reasons = c(error = 0L, non_finite = 0L, imputation = 3L),
                   error_messages = integer(), undetermined = 0L,
                   improper = 0L, caution = TRUE, require_half = TRUE,
                   imputation_messages = c("the chains diverged" = 3L))
  notes <- .boot_mi_notes(imp_only, B = 20L, conf_level = 0.95,
                          imputations = 2L, fmi = 0.5,
                          widening = sqrt(1 + 0.5 / 2),
                          sets_aside = sets * 0L, kept = 1)
  expect_match(notes, "imputations = 5 would bring the widening under 5%",
               fixed = TRUE, all = FALSE)
  expect_false(any(grepl("an even number", notes, fixed = TRUE)))
  # Replications set aside for a reason of their sets keep the even number.
  set_level <- imp_only
  set_level$reasons <- c(error = 0L, non_finite = 3L, imputation = 0L)
  set_level$improper <- 3L
  notes <- .boot_mi_notes(set_level, B = 20L, conf_level = 0.95,
                          imputations = 2L, fmi = 0.5,
                          widening = sqrt(1 + 0.5 / 2),
                          sets_aside = sets * 0L, kept = 1)
  expect_match(notes, "imputations = 6 would bring the widening under 5% (an even number",
               fixed = TRUE, all = FALSE)
})


test_that("with no bootstrap then impute interval there is no widening of one, in the warning or the printed line", {
  # Regression test. When fewer than half of the replications returned a
  # value, the warning said no interval was reported and then how much
  # imputation noise widened it, naming the imputations that would bring
  # the widening under 5%; the printed line did the same, under the line
  # that said no interval was reported.
  counts <- list(failed = 6L, returned = 4L,
                 reasons = c(error = 0L, non_finite = 0L, imputation = 6L),
                 error_messages = integer(), undetermined = 0L,
                 improper = 0L, caution = TRUE, require_half = TRUE,
                 imputation_messages = c("the chains diverged" = 6L))
  notes <- .boot_mi_notes(counts, B = 10L, conf_level = 0.95,
                          imputations = 2L, fmi = 0.98,
                          widening = sqrt(1 + 0.98 / 2))
  expect_match(notes, "fewer than half returned a value, so no interval is reported",
               fixed = TRUE, all = FALSE)
  expect_false(any(grepl("widens the interval", notes, fixed = TRUE)))
  expect_false(any(grepl("Caution", notes, fixed = TRUE)))
  x <- structure(
    data.frame(term = c("estimate", "conf_level"), value = c(.8, .95)),
    missing = "mi",
    imputation = list(impute = "pmm", imputations_point = 20L,
                      imputations_rule = 25L, imputations_floor = 20L,
                      imputations_replicate = 2L, fmi = 0.98,
                      fmi_ci_upper = 1, widening = sqrt(1 + 0.98 / 2),
                      status = rep("proper", 20L), polr_fallbacks = 0L),
    bootstrap = list(attempted = 10L, returned = 4L, failed = 6L,
                     require_half = TRUE))
  line <- .missing_treatment_footer_lines(x)[["imputation"]]
  expect_no_match(line, "widens the interval", fixed = TRUE)
  expect_match(line, paste0("fraction of missing information 0.98 (95% ",
                            "upper limit 1.00), estimated from the 4 ",
                            "replications that returned a value"),
               fixed = TRUE)
  # With an interval the widening is given as before.
  attr(x, "bootstrap")$returned <- 8L
  expect_match(.missing_treatment_footer_lines(x)[["imputation"]],
               "imputation noise widens the interval by about", fixed = TRUE)
})


# Auxiliary variables that duplicate other predictors -------------------

test_that("an auxiliary variable that is a linear function of other predictors is left out, with a message", {
  skip_if_not_installed("mice")
  few_imputations()
  # x is the standardized mean of the complete items y1 to y3, x2 = 2x +
  # 1. With mice's collinearity check turned off both reached every
  # regression, mice's remove.lindep() removed only one of the two exact
  # dependencies, and its ridge fallback made the imputations unstable:
  # 0.768 (between-imputation SD 0.062) against 0.800 (0.008) without
  # them, at 20 imputations. Neither carries information the items do
  # not, so the result is now the one without auxiliary variables.
  d <- cbind(mar, x = z_mar, x2 = 2 * z_mar + 1)
  ref <- suppressMessages(reliability_omega_categorical(
    mar, ci_method = "none", missing = "mi", imputations = 4, seed = 1))
  msgs <- testthat::capture_messages(
    r <- reliability_omega_categorical(d, ci_method = "none",
                                       aux = c("x", "x2"), imputations = 4,
                                       seed = 1))
  expect_match(msgs,
               paste0("Among the auxiliary variables, x is a linear function ",
                      "of y1, y2, and y3; x2 is a linear function of y1, y2, ",
                      "and y3"),
               fixed = TRUE, all = FALSE)
  expect_identical(est_of(r), est_of(ref))
  mi <- attr(r, "imputation")
  expect_identical(mi$aux_set_aside, c("x", "x2"))
  expect_identical(mi$ridge, c(point = 0L, replicates = 0L))
  # The result still names the auxiliary variables the user gave.
  expect_identical(attr(r, "aux"), c("x", "x2"))
})

test_that("the auxiliary screen keeps what carries information and finds what does not", {
  set.seed(5)
  n <- 200
  base <- data.frame(y1 = sample(1:5, n, TRUE), y2 = sample(1:5, n, TRUE),
                     y3 = sample(1:5, n, TRUE))
  base$y3[1:40] <- NA
  base$age <- round(stats::rnorm(n, 40, 10))
  base$months <- 12 * base$age
  base$noise <- stats::rnorm(n)
  # An auxiliary computed from an incomplete item and observed only where
  # that item is: a linear function of it on every row that observes it.
  base$copy <- 2 * base$y3
  # Observed on two rows only: too few to tell, so kept.
  base$rare <- NA_real_
  base$rare[1:2] <- c(1, 2)
  s <- .mi_redundant_aux(base, c("y1", "y2", "y3"),
                         c("age", "months", "noise", "copy", "rare"))
  expect_identical(s$aux, c("months", "copy"))
  expect_identical(s$why$months, "age")
  expect_identical(s$why$copy, "y3")
  # The order given decides which of a pair is kept.
  s2 <- .mi_redundant_aux(base, c("y1", "y2", "y3"), c("months", "age"))
  expect_identical(s2$aux, "age")
  # The reason given: no information beyond the predictors already in the
  # model. The instability of more than one exact dependency is named only
  # when more than one variable is left out; with one, mice removes it
  # itself and imputes as without it, so the message first gave a false
  # reason for the common case.
  m1 <- .mi_redundant_aux_message(s2)
  expect_match(m1, paste0("The auxiliary variable age is a linear function ",
                          "of months"), fixed = TRUE)
  expect_match(m1, "carries no information beyond the predictors already",
               fixed = TRUE)
  expect_false(grepl("more than one exact", m1, fixed = TRUE))
  expect_false(grepl("unstable", m1, fixed = TRUE))
  m2 <- .mi_redundant_aux_message(s)
  expect_match(m2, "these carry no information beyond the predictors",
               fixed = TRUE)
  expect_match(m2, "more than one exact linear dependency", fixed = TRUE)
  expect_identical(.mi_redundant_aux(base, c("y1", "y2", "y3"),
                                     character())$aux, character())
})

test_that("mice's ridge penalty and polr fallbacks are counted by imputation and reported", {
  events <- data.frame(
    it = c(1, 1, 2, 3), im = c(1, 2, 2, 3), dep = "y4", meth = "pmm",
    out = c("mice detected that your data are (nearly) multi-collinear.\nIt applied a ridge penalty",
            "polr falls back to multinom",
            "mice detected that your data are (nearly) multi-collinear.\nIt applied a ridge penalty",
            "x, x2"))
  got <- .mi_logged_counts(events, 3L)
  expect_identical(got$ridge, c(1L, 1L, 0L))
  expect_identical(got$polr, c(0L, 1L, 0L))
  expect_identical(.mi_logged_counts(NULL, 2L)$ridge, c(0L, 0L))
  skip_if_not_installed("mice")
  few_imputations()
  orig <- .mi_complete
  testthat::local_mocked_bindings(.mi_complete = function(...) {
    out <- orig(...)
    out$ridge_each <- c(2L, rep(0L, length(out$data) - 1L))
    out$ridge <- 2L
    out
  })
  w <- testthat::capture_warnings(r <- suppressMessages(
    reliability_omega_categorical(mar, ci_method = "none", missing = "mi",
                                  imputations = 2, seed = 1)))
  expect_identical(attr(r, "imputation")$ridge,
                   c(point = 2L, replicates = 0L))
  expect_match(w, "mice applied its ridge penalty 2 times in imputing the observed data",
               fixed = TRUE, all = FALSE)
})


# Factors and SPSS user-defined missing values --------------------------

test_that("a factor whose labels are whole numbers is the item those numbers code, under every treatment", {
  # A factor passed only when its labels were 1, 2, ..., K in order; a
  # Likert factor labeled 0 to 4 was refused as "non-integer values".
  f0 <- mar
  i0 <- mar
  for (j in names(mar)) {
    f0[[j]] <- factor(mar[[j]] - 1L, levels = 0:4, ordered = TRUE)
    i0[[j]] <- mar[[j]] - 1L
  }
  for (mm in c("listwise", "pairwise")) {
    a <- suppressMessages(reliability_omega_categorical(
      f0, ci_method = "none", missing = mm))
    b <- suppressMessages(reliability_omega_categorical(
      i0, ci_method = "none", missing = mm))
    expect_identical(a$value, b$value)
  }
  # A factor of the observed codes whose lowest category nobody chose.
  sk <- full
  sk$y1[sk$y1 == 1L] <- 2L
  fk <- sk
  fk$y1 <- factor(sk$y1)
  expect_identical(levels(fk$y1), c("2", "3", "4", "5"))
  expect_identical(
    suppressMessages(reliability_omega_categorical(fk, ci_method = "none"))$value,
    suppressMessages(reliability_omega_categorical(sk, ci_method = "none"))$value)
  # Words, or levels out of the order of their numbers, are refused.
  tf <- full
  tf$y2 <- factor(c("SD", "D", "N", "A", "SA")[tf$y2],
                  levels = c("SD", "D", "N", "A", "SA"), ordered = TRUE)
  expect_error(reliability_omega_categorical(tf, ci_method = "none"),
               "y2 (column 2) is a factor whose labels are not all whole numbers (\"SD\", \"D\", \"N\", \"A\", \"SA\")",
               fixed = TRUE)
  rf <- full
  rf$y1 <- factor(rf$y1, levels = 5:1, ordered = TRUE)
  expect_error(reliability_omega_categorical(rf, ci_method = "none"),
               "y1 (column 1) is an ordered factor whose levels",
               fixed = TRUE)
  skip_if_not_installed("mice")
  few_imputations()
  a <- suppressMessages(reliability_omega_categorical(
    f0, ci_method = "none", missing = "mi", imputations = 2, seed = 1))
  b <- suppressMessages(reliability_omega_categorical(
    i0, ci_method = "none", missing = "mi", imputations = 2, seed = 1))
  expect_identical(a$value, b$value)
})

test_that("SPSS user-defined missing values are missing under every treatment", {
  # haven::read_sav(user_na = TRUE) keeps a code declared missing (9,
  # "refused") as a value of a haven_labelled_spss column, built here
  # without haven. It was analyzed as a sixth category.
  spss <- function(x, na = 9) {
    structure(as.double(ifelse(is.na(x), na, x)), labels = c(Refused = na),
              na_values = na,
              class = c("haven_labelled_spss", "haven_labelled",
                        "vctrs_vctr", "double"))
  }
  sav <- mar
  sav$y4 <- spss(mar$y4)
  msgs <- testthat::capture_messages(
    r <- reliability_omega_categorical(sav, ci_method = "none"))
  expect_match(msgs,
               "treated as missing, as SPSS and haven's is.na() treat them: 130 in y4",
               fixed = TRUE, all = FALSE)
  expect_identical(r$value, suppressMessages(
    reliability_omega_categorical(mar, ci_method = "none"))$value)
  rp <- suppressMessages(reliability_omega_categorical(
    sav, ci_method = "none", missing = "pairwise"))
  expect_identical(rp$value, suppressMessages(reliability_omega_categorical(
    mar, ci_method = "none", missing = "pairwise"))$value)
  # A declared range is read the same way.
  rng <- mar
  rng$y5 <- structure(as.double(ifelse(is.na(mar$y5), 8, mar$y5)),
                      na_range = c(8, 9),
                      class = c("haven_labelled_spss", "haven_labelled",
                                "vctrs_vctr", "double"))
  expect_identical(sum(is.na(suppressMessages(
    .relia_declared_missing(rng))$y5)), 119L)
  skip_if_not_installed("mice")
  few_imputations()
  a <- suppressMessages(reliability_omega_categorical(
    sav, ci_method = "none", missing = "mi", imputations = 2, seed = 1))
  b <- suppressMessages(reliability_omega_categorical(
    mar, ci_method = "none", missing = "mi", imputations = 2, seed = 1))
  expect_identical(a$value, b$value)
})


test_that("completed categorical resamples name a constant item before lavaan fits it", {
  flat <- full
  flat$y2 <- 1L
  for (denominator in c("observed", "model_implied")) {
    point <- .omega_c_point_fn(denominator)
    output <- utils::capture.output(
      err <- tryCatch(point(flat), error = identity))
    expect_length(output, 0L)
    expect_s3_class(err, "error")
    if (inherits(err, "error")) {
      expect_match(conditionMessage(err), "y2 takes a single value", fixed = TRUE)
    }
  }
})
