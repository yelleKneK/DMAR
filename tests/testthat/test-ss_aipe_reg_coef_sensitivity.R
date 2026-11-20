## The scale of the coefficient in ss_aipe_reg_coef_sensitivity(),
## ss_aipe_rc_sensitivity(), ss_aipe_src_sensitivity(), and
## ss_power_reg_coef_sensitivity(). The name of the function fixes the scale
## (Ken's ruling of 2026-11-11): the first two and the power function report
## the unstandardized coefficient, ss_aipe_src_sensitivity() the
## standardized one, and none takes the DMAR 1.0.x argument `standardize`.
## ss_aipe_src_sensitivity() standardizes every replication, so its estimate
## and interval are for the standardized coefficient; the true value they are
## judged against, the supposed value echoed beside it, and the planned
## width are on that scale too. Through the first 1.1.0 builds (and in DMAR
## 1.0.x) the true and supposed values stayed unstandardized and the plan was
## made with the supposed standard deviations, so for var(Y) = 4,
## cov(Y, X) = (.8, .6), and unit predictor variances correlated .3, the
## intervals for the standardized coefficient .341 were judged against .681
## (every interval missed, total_type_I_error 1) and the plan for width .3
## was 610, the plan for the unstandardized coefficient, where the
## standardized one needs 157.

.probe <- list(var_Y = 4, cov_YX = c(.8, .6),
               cov_XX = matrix(c(1, .3, .3, 1), 2))

# The oracle, by hand from the covariances (Kelley & Maxwell, 2003, p. 316):
# b = solve(cov_XX, cov_YX) and beta_j = b_j sd_Xj / sd_Y.
.hand <- function(var_Y, cov_YX, cov_XX, j = 1) {
  b <- solve(cov_XX, cov_YX)[j]
  c(b = b, beta = b * sqrt(cov_XX[j, j] / var_Y))
}

.sens <- function(fn = ss_aipe_reg_coef_sensitivity, s = .probe, ...) {
  r <- fn(true_var_Y = s$var_Y, true_cov_YX = s$cov_YX,
          true_cov_XX = s$cov_XX, print_iter = FALSE, ...)
  r
}
.v <- function(r) stats::setNames(r$value, r$term)

test_that("the standardize argument is gone: a call that gives it stops", {
  for (fn in c("ss_aipe_reg_coef_sensitivity", "ss_aipe_rc_sensitivity",
               "ss_aipe_src_sensitivity", "ss_power_reg_coef_sensitivity")) {
    f <- get(fn, envir = asNamespace("DMAR"))
    expect_false("standardize" %in% names(formals(f)), info = fn)
    for (value in c(TRUE, FALSE)) {
      expect_error(f(true_var_Y = 1, true_cov_YX = c(.4, .3),
                     true_cov_XX = matrix(c(1, .3, .3, 1), 2),
                     specified_N = 40, G = 2, print_iter = FALSE,
                     standardize = value),
                   "unused argument", info = paste(fn, value))
    }
  }
})

test_that("ss_aipe_src_sensitivity() reports the standardized coefficient and ss_aipe_reg_coef_sensitivity() the unstandardized one", {
  skip_on_cran()
  h <- .hand(.probe$var_Y, .probe$cov_YX, .probe$cov_XX)
  expect_equal(unname(h[["beta"]]), .62 / .91 / 2)
  std <- .v(.sens(ss_aipe_src_sensitivity, specified_N = 60, width = .3,
                  G = 20, seed = 1))
  raw <- .v(.sens(specified_N = 60, width = .3, G = 20, seed = 1))
  rc  <- .v(.sens(ss_aipe_rc_sensitivity, specified_N = 60, width = .3,
                  G = 20, seed = 1))
  expect_equal(std[["true_std_reg_coef_j"]], h[["beta"]], tolerance = 1e-12)
  expect_equal(raw[["true_reg_coef_j"]], h[["b"]], tolerance = 1e-12)
  # ss_aipe_rc_sensitivity() is the same analysis: the same draws give the
  # same table.
  expect_identical(rc, raw)
  # The supposed covariances, here different from the true ones, are put on
  # the scale of the function with the supposed standard deviations.
  sup <- list(var_Y = 9, cov_YX = c(1.2, .3), cov_XX = matrix(c(4, .6, .6, 1), 2))
  hs <- .hand(sup$var_Y, sup$cov_YX, sup$cov_XX)
  args <- list(supposed_var_Y = sup$var_Y, supposed_cov_YX = sup$cov_YX,
               supposed_cov_XX = sup$cov_XX, specified_N = 60, width = .3,
               G = 5, seed = 1)
  both_std <- .v(do.call(.sens, c(list(fn = ss_aipe_src_sensitivity), args)))
  both_raw <- .v(do.call(.sens, args))
  expect_equal(both_std[["supposed_std_reg_coef_j"]], hs[["beta"]], tolerance = 1e-12)
  expect_equal(both_raw[["supposed_reg_coef_j"]], hs[["b"]], tolerance = 1e-12)
  # The rows are named for the coefficient computed, and none for the
  # other scale.
  expect_true(all(paste0(c("mean_", "median_", "sd_", "true_", "supposed_"),
                         "std_reg_coef_j") %in% names(std)))
  expect_false(any(grepl("(^|_)reg_coef_j$", sub("std_reg_coef_j$", "", names(std)))))
  expect_true(all(paste0(c("mean_", "median_", "sd_", "true_", "supposed_"),
                         "reg_coef_j") %in% names(raw)))
  expect_false(any(grepl("std_reg_coef_j", names(raw))))
})

test_that("the misses are counted against the standardized truth, replication by replication", {
  skip_on_cran()
  # The per-replication file holds every interval; counting its misses of
  # the hand-computed standardized coefficient reproduces the tail rows
  # exactly.
  h <- .hand(.probe$var_Y, .probe$cov_YX, .probe$cov_XX)
  for (noncentral in c(FALSE, TRUE)) {
    path <- tempfile(fileext = ".csv")
    r <- .v(.sens(ss_aipe_src_sensitivity, specified_N = 200, width = .3,
                  noncentral = noncentral, G = 300, seed = 2026,
                  filename = path))
    f <- utils::read.csv(path)
    unlink(path)
    expect_named(f, c("std_reg_coef_j", "ci_lower", "ci_upper", "R2",
                      "se_std_reg_coef_j", "t_value"))
    expect_equal(r[["pct_ci_miss_low"]], mean(h[["beta"]] < f$ci_lower))
    expect_equal(r[["pct_ci_miss_high"]], mean(h[["beta"]] > f$ci_upper))
    expect_equal(r[["mean_std_reg_coef_j"]], mean(f$std_reg_coef_j))
    # The probe of the earlier repair: the 1.0.x truth .681 was missed by
    # every interval; the standardized truth .341 is covered at about the
    # nominal rate (the oracle simulation in the c1_regcoef evidence gives
    # .969 and .973 for these two intervals with 4000 replications).
    expect_lt(r[["total_type_I_error"]], .10)
    expect_equal(mean(h[["b"]] < f$ci_lower | h[["b"]] > f$ci_upper), 1)
  }
  # The unstandardized file names the unstandardized coefficient.
  path <- tempfile(fileext = ".csv")
  .sens(specified_N = 40, width = .3, G = 3, seed = 1, filename = path)
  expect_named(utils::read.csv(path), c("reg_coef_j", "ci_lower", "ci_upper",
                                        "R2", "se_reg_coef_j", "t_value"))
  unlink(path)
})

test_that("coverage of each coefficient matches an independent simulation", {
  skip_on_cran()
  # The oracle: 4000 replications of the probe structure at N = 200,
  # simulated without DMAR (Cholesky draws, standardization, least
  # squares, the central t interval), gave coverage .9688 (Monte Carlo SE
  # .0028) of the standardized coefficient and .9550 (.0033) of the
  # unstandardized one (evidence/c1_regcoef, oracle_coverage.log). The
  # package's coverage over 1000 replications must agree within four
  # combined Monte Carlo standard errors.
  for (case in list(list(fn = ss_aipe_src_sensitivity, oracle = .9688, se = .0028),
                    list(fn = ss_aipe_reg_coef_sensitivity, oracle = .9550, se = .0033))) {
    r <- .sens(case$fn, specified_N = 200, width = .3, noncentral = FALSE,
               G = 1000, seed = 99)
    sim <- attr(r, "simulation")
    expect_equal(sim$attempted, 1000L)
    expect_equal(sim$failed, 0L)
    cover <- 1 - .v(r)[["total_type_I_error"]]
    mcse <- sqrt(cover * (1 - cover) / sim$returned)
    expect_lt(abs(cover - case$oracle), 4 * sqrt(mcse^2 + case$se^2))
  }
})

test_that("the standardized truth follows a reversed sign and is unchanged by the units", {
  skip_on_cran()
  h <- .hand(.probe$var_Y, .probe$cov_YX, .probe$cov_XX)
  # X1 reflected: the coefficient changes sign, the size does not.
  flipped <- list(var_Y = 4, cov_YX = c(-.8, .6),
                  cov_XX = matrix(c(1, -.3, -.3, 1), 2))
  rf <- .v(.sens(ss_aipe_src_sensitivity, s = flipped, specified_N = 60,
                 width = .3, G = 5, seed = 1))
  expect_equal(rf[["true_std_reg_coef_j"]], -h[["beta"]], tolerance = 1e-12)
  rfr <- .v(.sens(s = flipped, specified_N = 60, width = .3, G = 5, seed = 1))
  expect_equal(rfr[["true_reg_coef_j"]], -h[["b"]], tolerance = 1e-12)
  # Y in units 10 times smaller and X1 in units 3 times smaller: the
  # unstandardized coefficient is multiplied by 10/3, the standardized one
  # is unchanged.
  units <- list(var_Y = 400, cov_YX = c(24, 6), cov_XX = matrix(c(9, .9, .9, 1), 2))
  ru_std <- .v(.sens(ss_aipe_src_sensitivity, s = units, specified_N = 60,
                     width = .3, G = 5, seed = 1))
  ru_raw <- .v(.sens(s = units, specified_N = 60, width = .3, G = 5, seed = 1))
  expect_equal(ru_std[["true_std_reg_coef_j"]], h[["beta"]], tolerance = 1e-12)
  expect_equal(ru_raw[["true_reg_coef_j"]], h[["b"]] * 10 / 3, tolerance = 1e-12)
  # A third predictor, the second one targeted, unequal variances.
  three <- list(var_Y = 2.5, cov_YX = c(.5, -.9, .3),
                cov_XX = matrix(c(1, .4, .2, .4, 2, -.3, .2, -.3, .5), 3))
  h3 <- .hand(three$var_Y, three$cov_YX, three$cov_XX, j = 2)
  r3 <- .v(.sens(ss_aipe_src_sensitivity, s = three, specified_N = 60,
                 which_predictor = 2, width = .3, G = 5, seed = 1))
  expect_equal(r3[["true_std_reg_coef_j"]], h3[["beta"]], tolerance = 1e-12)
  r3r <- .v(.sens(ss_aipe_rc_sensitivity, s = three, specified_N = 60,
                  which_predictor = 2, width = .3, G = 5, seed = 1))
  expect_equal(r3r[["true_reg_coef_j"]], h3[["b"]], tolerance = 1e-12)
})

test_that("the plan for a supposed structure is made on the scale of the function", {
  skip_on_cran()
  cor_all <- cov2cor(rbind(c(.probe$var_Y, .probe$cov_YX),
                           cbind(.probe$cov_YX, .probe$cov_XX)))
  plan_std <- function(...) ss_aipe_reg_coef(
    width = .3, true_cor_XX = cor_all[-1, -1], true_cor_YX = cor_all[1, -1],
    which_predictor = 1, ...)$value[1]
  sup <- function(fn = ss_aipe_reg_coef_sensitivity, ...)
    .v(.sens(fn, supposed_var_Y = .probe$var_Y, supposed_cov_YX = .probe$cov_YX,
             supposed_cov_XX = .probe$cov_XX, width = .3, G = 2, seed = 1, ...))
  expect_equal(plan_std(), 156)
  # The standardized coefficient is planned with unit standard deviations.
  expect_equal(sup(ss_aipe_src_sensitivity, noncentral = FALSE)[["total_N"]],
               plan_std())
  # With assurance the central interval for the standardized coefficient
  # is planned by ss_aipe_src(noncentral = FALSE), not by the F law of the
  # interval in the units of the data (187 here, which the first 1.1.0
  # builds returned): by simulation, 171, from the draws the seed fixes,
  # and by default the chi square approximation, 181.
  expect_equal(sup(ss_aipe_src_sensitivity, noncentral = FALSE, assurance = .9,
                   assurance_method = "simulation")[["total_N"]],
               ss_aipe_src(width = .3, true_cor_XX = cor_all[-1, -1],
                           true_cor_YX = cor_all[1, -1], which_predictor = 1,
                           noncentral = FALSE, assurance = .9,
                           assurance_method = "simulation", seed = 1)$value[1])
  expect_equal(sup(ss_aipe_src_sensitivity, noncentral = FALSE, assurance = .9,
                   assurance_method = "simulation")[["total_N"]], 171)
  expect_equal(sup(ss_aipe_src_sensitivity, noncentral = FALSE,
                   assurance = .9)[["total_N"]],
               ss_aipe_src(width = .3, true_cor_XX = cor_all[-1, -1],
                           true_cor_YX = cor_all[1, -1], which_predictor = 1,
                           noncentral = FALSE, assurance = .9)$value[1])
  expect_equal(sup(ss_aipe_src_sensitivity, noncentral = FALSE,
                   assurance = .9)[["total_N"]], 181)
  # The unstandardized one in the units of the data: sd_Y / sd_X = 2
  # doubles the standard error.
  expect_equal(sup()[["total_N"]], plan_std(true_sd_Y = 2, true_sd_X = 1))
  expect_equal(sup(ss_aipe_rc_sensitivity)[["total_N"]],
               plan_std(true_sd_Y = 2, true_sd_X = 1))
  # The default of ss_aipe_src_sensitivity() (noncentral) plans the
  # noncentral standardized coefficient; it was refused through the first
  # 1.1.0 builds whenever the supposed variance of Y was not 1.
  expect_equal(sup(ss_aipe_src_sensitivity)[["total_N"]],
               plan_std(noncentral = TRUE))
  # Noncentral planning for the unstandardized coefficient is refused by
  # name: the noncentral plan is for a standardized coefficient.
  for (fn in list(ss_aipe_reg_coef_sensitivity, ss_aipe_rc_sensitivity)) {
    e <- tryCatch(sup(fn, noncentral = TRUE), error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_true(all(c("noncentral", "supposed_var_Y") %in% e$arg))
    expect_match(conditionMessage(e), "ss_aipe_src_sensitivity()", fixed = TRUE)
  }
  # With unit supposed variances the two scales coincide and the call runs.
  cor_s <- .v(.sens(supposed_var_Y = 1, supposed_cov_YX = cor_all[1, -1],
                    supposed_cov_XX = cor_all[-1, -1], width = .3,
                    noncentral = TRUE, G = 2, seed = 1))
  expect_equal(cor_s[["total_N"]], plan_std(noncentral = TRUE))
})

test_that("ss_aipe_src_sensitivity(noncentral = FALSE) with assurance simulates at the plan that attains it", {
  skip_on_cran()
  # Unit variances, the probe of test-ss_aipe_reg_coef.R: the plan is 171,
  # where raw multivariate normal data standardized and fit without DMAR
  # give .9136 (Monte Carlo standard error .0013; 50,000 samples;
  # work/round6/evidence/g1_final_fixes/03_src_central_oracle.R). Through
  # the first 1.1.0 builds the plan was 187, where the proportion is .997.
  R <- matrix(c(1, .3, .3, 1), 2)
  r <- ss_aipe_src_sensitivity(true_var_Y = 1, true_cov_YX = c(.4, .3), true_cov_XX = R,
                               width = .3, noncentral = FALSE, assurance = .9,
                               assurance_method = "simulation",
                               G = 2000, print_iter = FALSE, seed = 5)
  v <- .v(r)
  expect_equal(v[["total_N"]], 171)
  # Every replication returns an interval (attempted 2000, failed 0).
  expect_lt(abs(v[["pct_ci_less_width"]] - .9136),
            4 * sqrt(.9136 * (1 - .9136) / 2000 + .0013^2))
})

test_that("each function names its coefficient rows for its own scale", {
  skip_on_cran()
  rows <- function(r) grep("reg_coef_j$", r$term, value = TRUE)
  std_rows <- paste0(c("mean_", "median_", "sd_", "true_", "supposed_"),
                     "std_reg_coef_j")
  raw_rows <- sub("std_", "", std_rows)
  expect_setequal(rows(.sens(ss_aipe_src_sensitivity, specified_N = 40, width = .3,
                             G = 3, seed = 1)), std_rows)
  expect_setequal(rows(.sens(ss_aipe_src_sensitivity, specified_N = 40, width = .3,
                             noncentral = FALSE, G = 3, seed = 1)), std_rows)
  expect_setequal(rows(.sens(ss_aipe_rc_sensitivity, specified_N = 40, width = .3,
                             G = 3, seed = 1)), raw_rows)
  expect_setequal(rows(.sens(ss_aipe_reg_coef_sensitivity, specified_N = 40,
                             width = .3, noncentral = TRUE, G = 3, seed = 1)),
                  raw_rows)
})

test_that("each help page names the per-replication columns its own function writes", {
  skip_on_cran()
  # The oracle is the file each function writes. The filename entry of each
  # page must name those six columns and no coefficient column of the other
  # scale. The rc and src pages once said the file was "forwarded to"
  # ss_aipe_reg_coef_sensitivity(), whose entry lists the unstandardized
  # columns, which are not the ones ss_aipe_src_sensitivity() writes.
  filename_words <- function(topic) {
    path <- file.path("..", "..", "man", paste0(topic, ".Rd"))
    rd <- if (file.exists(path)) tools::parse_Rd(path) else
      tools::Rd_db("DMAR")[[paste0(topic, ".Rd")]]
    txt <- paste(utils::capture.output(tools::Rd2txt(rd)), collapse = " ")
    from <- regexpr("filename: ", txt, fixed = TRUE, useBytes = TRUE)
    to <- regexpr(" seed: ", txt, fixed = TRUE, useBytes = TRUE)
    expect_true(from > 0 && to > from, info = topic)
    entry <- substr(txt, from + 10L, to)
    strsplit(entry, "[^A-Za-z0-9_]+", useBytes = TRUE)[[1]]
  }
  for (fn in c("ss_aipe_reg_coef_sensitivity", "ss_aipe_rc_sensitivity",
               "ss_aipe_src_sensitivity")) {
    path <- tempfile(fileext = ".csv")
    .sens(get(fn), specified_N = 40, width = .3, noncentral = FALSE, G = 3,
          seed = 1, filename = path)
    written <- names(utils::read.csv(path))
    unlink(path)
    words <- filename_words(fn)
    expect_length(written, 6L)
    for (col in written) expect_true(col %in% words, info = paste(fn, col))
    other <- if (fn == "ss_aipe_src_sensitivity") {
      c("reg_coef_j", "se_reg_coef_j")
    } else {
      c("std_reg_coef_j", "se_std_reg_coef_j")
    }
    expect_false(any(other %in% words), info = fn)
    expect_false("forwarded" %in% words, info = fn)
    expect_true("appended" %in% words, info = fn)
  }
})

test_that("ss_aipe_reg_coef_sensitivity() refuses a missing variance, a predictor out of range, and a size without residual degrees of freedom", {
  e <- tryCatch(ss_aipe_reg_coef_sensitivity(
    true_cov_YX = .probe$cov_YX, true_cov_XX = .probe$cov_XX, specified_N = 40,
    width = .3, G = 2, print_iter = FALSE), error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "true_var_Y")
  e <- tryCatch(.sens(specified_N = 40, which_predictor = 3, width = .3, G = 2),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "which_predictor")
  expect_error(.sens(specified_N = 3, width = .3, G = 2),
               class = "dmar_argument_invalid")
  expect_error(.sens(specified_N = 40.5, width = .3, G = 2),
               class = "dmar_argument_invalid")
  expect_error(.sens(ss_aipe_src_sensitivity, specified_N = 40.5, width = .3,
                     G = 2), class = "dmar_argument_invalid")
})

test_that("ss_power_reg_coef_sensitivity() reports the unstandardized coefficient, and its t statistics are those of the standardized data", {
  skip_on_cran()
  h <- .hand(.probe$var_Y, .probe$cov_YX, .probe$cov_XX)
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  G <- 200; N <- 60
  r <- suppressMessages(ss_power_reg_coef_sensitivity(
    true_var_Y = .probe$var_Y, true_cov_YX = .probe$cov_YX,
    true_cov_XX = .probe$cov_XX, specified_N = N, G = G,
    print_iter = FALSE, seed = 5, filename = path))
  raw <- .v(r)
  f <- utils::read.csv(path)
  expect_named(f, c("reg_coef_j", "se_reg_coef_j", "t_value", "R2"))
  expect_equal(raw[["true_reg_coef_j"]], h[["b"]], tolerance = 1e-12)
  expect_equal(raw[["supposed_reg_coef_j"]], h[["b"]], tolerance = 1e-12)
  expect_true(all(paste0(c("mean_", "median_", "sd_", "mean_se_", "median_se_",
                           "sd_se_", "true_", "supposed_"), "reg_coef_j")
                  %in% names(raw)))
  expect_false(any(grepl("std_reg_coef_j", names(raw))))
  expect_lt(abs(raw[["mean_reg_coef_j"]] - h[["b"]]), 4 * raw[["sd_reg_coef_j"]] / sqrt(G))
  # The oracle: the same draws (MASS::mvrnorm under the same seed), each
  # standardized and fit by lm() here, give the same t statistics, so the
  # power is that of the test of the standardized coefficient too.
  t_std <- withr::with_seed(5, vapply(seq_len(G), function(i) {
    d <- scale(MASS::mvrnorm(N, mu = rep(0, 3),
                             Sigma = rbind(c(.probe$var_Y, .probe$cov_YX),
                                           cbind(.probe$cov_YX, .probe$cov_XX))))
    summary(lm(d[, 1] ~ d[, -1]))$coefficients[2, 3]
  }, numeric(1)))
  expect_equal(f$t_value, t_std, tolerance = 1e-10)
  crit <- stats::qt(.975, N - 3)
  expect_equal(raw[["empirical_power"]], mean(abs(t_std) > crit))
})

test_that("ss_aipe_src_sensitivity() plans and judges on the standardized scale whatever the variances, without a warning", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  R <- matrix(c(1, .3, .3, 1), 2)
  # Each call's total_N and true standardized coefficient.
  run <- function(..., which_predictor = 1) {
    out <- expect_no_warning(ss_aipe_src_sensitivity(
      ..., which_predictor = which_predictor, width = .3, G = 2,
      print_iter = FALSE, seed = 1))
    list(N = out$value[out$term == "total_N"],
         true = out$value[out$term == "true_std_reg_coef_j"])
  }
  # The variances of Y and of the predictor differ.
  var4 <- list(true_var_Y = 4, true_cov_YX = c(.8, .6), true_cov_XX = R)
  r <- do.call(run, c(var4, noncentral = FALSE))
  expect_identical(r$N, 156)
  expect_equal(r$true, .3407, tolerance = 1e-3)
  for (m in list(NULL, "analytic")) {
    expect_identical(do.call(run, c(var4, noncentral = FALSE, assurance = .9,
                                    assurance_method = m))$N, 181)
  }
  expect_identical(do.call(run, c(var4, noncentral = FALSE, assurance = .9,
                                  assurance_method = "simulation"))$N, 171)
  # Only the predictor's variance differs.
  varX <- list(true_var_Y = 1, true_cov_YX = c(.8, .3), true_cov_XX = matrix(c(4, .6, .6, 1), 2))
  expect_identical(do.call(run, c(varX, noncentral = FALSE, assurance = .9))$N, 181)
  # The true variance of Y differs and the supposed values are unit ones.
  mixed <- c(var4, list(supposed_var_Y = 1, supposed_cov_YX = c(.4, .3), supposed_cov_XX = R))
  expect_identical(do.call(run, c(mixed, assurance = .9))$N, 179)
  # Every variance 4: the results of unit variances.
  all4 <- list(true_var_Y = 4, true_cov_YX = 4 * c(.4, .3), true_cov_XX = 4 * R)
  expect_identical(do.call(run, c(all4, noncentral = FALSE))$N, 156)
  expect_identical(do.call(run, c(all4, noncentral = FALSE, assurance = .9))$N, 181)
  # A predictor that is not targeted has another variance (true, supposed,
  # or both): the results of unit variances.
  X2var4 <- matrix(c(1, .6, .6, 4), 2)
  X1var4 <- matrix(c(4, .6, .6, 1), 2)
  cases <- list(
    both = list(s = list(true_var_Y = 1, true_cov_YX = c(.4, .6), true_cov_XX = X2var4),
                N = c(156, 181, 164, 179), true = .3407),
    supposed = list(s = list(true_var_Y = 1, true_cov_YX = c(.4, .3), true_cov_XX = R,
                             supposed_var_Y = 1, supposed_cov_YX = c(.4, .6),
                             supposed_cov_XX = X2var4),
                    N = c(156, 181, 164, 179), true = .3407),
    true = list(s = list(true_var_Y = 1, true_cov_YX = c(.4, .6), true_cov_XX = X2var4,
                         supposed_var_Y = 1, supposed_cov_YX = c(.4, .3),
                         supposed_cov_XX = R),
                N = c(156, 181, 164, 179), true = .3407),
    second = list(s = list(true_var_Y = 1, true_cov_YX = c(.8, .3), true_cov_XX = X1var4,
                           which_predictor = 2),
                  N = c(156, 181, 157, 173), true = .1978))
  for (nm in names(cases)) {
    s <- cases[[nm]]$s
    N <- cases[[nm]]$N
    r <- do.call(run, c(s, noncentral = FALSE))
    expect_identical(r$N, N[1], info = nm)
    expect_equal(r$true, cases[[nm]]$true, tolerance = 1e-3, info = nm)
    expect_identical(do.call(run, c(s, noncentral = FALSE, assurance = .9))$N,
                     N[2], info = nm)
    expect_identical(do.call(run, s)$N, N[3], info = nm)
    expect_identical(do.call(run, c(s, assurance = .9))$N, N[4], info = nm)
  }
  # Y and the targeted predictor of variance 4 and the other predictor of
  # variance 1.
  s <- list(true_var_Y = 4, true_cov_YX = c(1.6, .6), true_cov_XX = X1var4)
  expect_identical(do.call(run, c(s, noncentral = FALSE))$N, 156)
  expect_identical(do.call(run, c(s, noncentral = FALSE, assurance = .9))$N, 181)
})
