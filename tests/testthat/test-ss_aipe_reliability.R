## ss_aipe_reliability() -- sample size for the alpha / omega reliability CI.

test_that("ss_aipe_reliability(parallel/Normal Theory) returns a tidy necessary_N", {
  skip_on_cran()
  set.seed(113)
  res <- suppressMessages(suppressWarnings(
    ss_aipe_reliability(model = "parallel", type = "Normal Theory",
                        width = 0.10, n_items = 4, true_factor_loading = 0.7, true_residual_var = 0.51,
                        G_initial = 50, G = 200)
  ))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_N" %in% res$term)
  expect_gt(res$value[res$term == "necessary_N"], 0)
})

test_that("ss_aipe_reliability(Congeneric / Factor Analytic) plans a sample size", {
  # Regression test: the Factor Analytic path read the legacy cfa_1() list
  # layout ($factor_loadings / $parameter_cov), which no longer exists, so the
  # path errored. It now routes through the maintained omega internals.
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(113)
  res <- suppressMessages(suppressWarnings(
    ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
                        width = 0.10, n_items = 5,
                        true_factor_loading = c(.4, .4, .3, .3, .5),
                        true_residual_var = c(.2, .4, .3, .3, .2),
                        conf_level = 0.95)
  ))
  expect_s3_class(res, "data.frame")
  expect_true("necessary_N" %in% res$term)
  expect_gt(res$value[res$term == "necessary_N"], 0)
})

test_that("ss_aipe_reliability() requires more N for a tighter width target", {
  skip_on_cran()
  set.seed(113)
  wide <- suppressMessages(suppressWarnings(
    ss_aipe_reliability(model = "parallel", type = "Normal Theory",
                        width = 0.20, n_items = 4, true_factor_loading = 0.7, true_residual_var = 0.51,
                        G_initial = 50, G = 200)
  ))$value[1]
  set.seed(113)
  tight <- suppressMessages(suppressWarnings(
    ss_aipe_reliability(model = "parallel", type = "Normal Theory",
                        width = 0.05, n_items = 4, true_factor_loading = 0.7, true_residual_var = 0.51,
                        G_initial = 50, G = 200)
  ))$value[1]
  expect_gt(tight, wide)
})

## Planning under missing data: fmi and seed (1.1.0).

# The complete data closed forms, written out independently of the function
# so the fmi oracles do not rest on the code they test. Parallel model,
# normal theory: N - 1 = z^2 8 (1 - alpha)^2 p / (w^2 (p - 1)). Congeneric
# model (1.1.0): N - 1 = 4 z^2 V / w^2, with V / N the squared standard
# error of reliability_omega(denominator = "model_implied", ci_method =
# "ml"): V = 0.0739065022 read from that function (lavaan 0.7-3's delta
# method) at N = 1e5 and 1e6, which agree to 1e-13; the test of the
# congeneric closed form below checks it against lavaan again.
.closed_form_parallel <- function() {
  S <- covmat_from_cfa(factor_loading = matrix(sqrt(.3), 1, 6),
                       residual_var = matrix(.2, 1, 6),
                       tol_det = 1e-5)$model_implied_cov
  p <- 6
  a <- p / (p - 1) * (1 - sum(diag(S)) / sum(S))
  z <- stats::qnorm(.975)
  z^2 * 8 * (1 - a)^2 * p / (.1^2 * (p - 1)) + 1
}
.closed_form_congeneric <- function() {
  V <- 0.0739065022
  z <- stats::qnorm(.975)
  4 * z^2 * V / .15^2 + 1
}
.plan_parallel <- function(...) {
  ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
                      n_items = 6, true_item_covariance = .3, true_residual_var = .2, ...)
}
.plan_congeneric <- function(...) {
  ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
                      width = .15, n_items = 4, true_factor_loading = c(.8, .7, .7, .8),
                      true_residual_var = c(.4, .5, .5, .4), ...)
}

test_that("fmi divides the complete data N - 1 by 1 - fmi (closed form oracles)", {
  # Unrounded complete data answers 37.878005 and 51.472673 (52.549566,
  # alpha's Eq. 19 on the correlation matrix, before 1.1.0).
  Nc_par <- .closed_form_parallel()
  Nc_con <- .closed_form_congeneric()
  expect_equal(Nc_par, 37.878005, tolerance = 1e-7)
  expect_equal(Nc_con, 51.472673, tolerance = 1e-6)
  fmis <- c(0, .2, .3, .5)
  par <- vapply(fmis, function(f) .plan_parallel(fmi = f)$value, numeric(1))
  con <- vapply(fmis, function(f) .plan_congeneric(fmi = f)$value, numeric(1))
  expect_identical(par, c(38, 48, 54, 75))
  expect_identical(con, c(52, 65, 74, 102))
  expect_identical(par, ceiling((Nc_par - 1) / (1 - fmis) + 1))
  expect_identical(con, ceiling((Nc_con - 1) / (1 - fmis) + 1))
  # The factor scales N - 1, not N: the naive ceiling(N_c / (1 - fmi)) is
  # larger at .3 and .5 for both plans, by a case or two.
  expect_identical(ceiling(ceiling(Nc_par) / (1 - c(.3, .5))), c(55, 76))
  expect_identical(ceiling(ceiling(Nc_con) / (1 - c(.3, .5))), c(75, 104))
  # The plan maps the unrounded N_c, so the documented formula applied to
  # the printed whole number can exceed the plan, as the help says: 103
  # against the plan's 102 at an FMI of .5, and 43 against 42 at .1.
  out <- .plan_congeneric(fmi = .5)
  expect_identical(attr(out, "necessary_N_complete"), 52)
  expect_identical(ceiling((52 - 1) / (1 - .5) + 1), 103)
  out <- .plan_parallel(fmi = .1)
  expect_identical(out$value, 42)
  expect_identical(attr(out, "necessary_N_complete"), 38)
  expect_identical(ceiling((38 - 1) / (1 - .1) + 1), 43)
})

test_that("fmi = 0 returns the complete data plan unchanged, with no new attribute", {
  # The 1.0.x object, built the way the function built it before fmi.
  expected <- DMAR:::.as_dmar_tbl(data.frame(term = "necessary_N", value = 38),
                                  conf_level = 0.95, subclass = "dmar_ss_aipe")
  expect_identical(.plan_parallel(), expected)
  expect_identical(.plan_parallel(fmi = 0), expected)
  expect_null(attr(.plan_parallel(), "fmi"))
  expect_null(attr(.plan_parallel(), "necessary_N_complete"))
  expect_identical(.plan_congeneric(fmi = 0), .plan_congeneric())
})

test_that("a plan with fmi > 0 keeps its rows and records the complete data N", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  out <- .plan_parallel(fmi = .3)
  expect_named(out, c("term", "value"))
  expect_identical(out$term, "necessary_N")
  expect_identical(class(out), class(.plan_parallel()))
  expect_identical(attr(out, "fmi"), .3)
  expect_identical(attr(out, "necessary_N_complete"), 38)
  expect_output(print(out),
                "Planned for a fraction of missing information of 0.3 for the coefficient; with complete data, N = 38.",
                fixed = TRUE)
  expect_false(any(grepl("missing information",
                         utils::capture.output(print(.plan_parallel())))))
})

test_that("fmi must be a single number in [0, 1)", {
  for (bad in list(1, -0.1, c(.1, .2), NA_real_, Inf, "0.3")) {
    expect_error(.plan_parallel(fmi = bad), "'fmi' must be a single number in [0, 1)",
                 fixed = TRUE)
  }
  expect_error(.plan_parallel(seed = "a"), "'seed' must be a single number or NULL",
               fixed = TRUE)
})

test_that("a seed makes the assurance search reproducible and restores the caller's stream", {
  skip_on_cran()
  args <- list(assurance = .80, G_initial = 20, G = 40)
  # seed = 113 draws what set.seed(113) before the call draws.
  set.seed(113)
  unseeded <- do.call(.plan_parallel, args)
  set.seed(42)
  before <- get(".Random.seed", envir = globalenv())
  seeded <- do.call(.plan_parallel, c(args, seed = 113))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_identical(seeded, unseeded)
  expect_identical(do.call(.plan_parallel, c(args, seed = 113)), seeded)
  # Without a seed the search draws from, and advances, the current stream.
  set.seed(42)
  invisible(do.call(.plan_parallel, args))
  expect_false(identical(get(".Random.seed", envir = globalenv()), before))
})

test_that("the assurance search runs at the complete data size and maps back", {
  skip_on_cran()
  args <- list(assurance = .80, G_initial = 20, G = 40, seed = 7)
  complete <- do.call(.plan_parallel, args)
  N_e <- complete$value[complete$term == "necessary_N"]
  # The planned size in exact integer arithmetic: the smallest N with
  # (N - 1)(1 - fmi) >= N_e - 1, for fmi = .2 (4/5) and .5 (1/2).
  exact <- c(`0.2` = ceiling((N_e - 1) * 5 / 4) + 1, `0.5` = (N_e - 1) * 2 + 1)
  for (f in c(.2, .5)) {
    out <- do.call(.plan_parallel, c(args, fmi = f))
    expect_identical(out$term, complete$term)
    expect_identical(out$value[out$term == "necessary_N"],
                     unname(exact[as.character(f)]))
    # Every other row is the complete data search's, which ran at N_e.
    expect_identical(out$value[-1], complete$value[-1])
    expect_identical(attr(out, "necessary_N_complete"), N_e)
    expect_identical(attr(out, "fmi"), f)
    expect_output(print(out), "the size at which the assurance was simulated",
                  fixed = TRUE)
  }
  # A start_N is a size under missing data, converted to its effective
  # size, ceiling((start_N - 1)(1 - fmi) + 1), before the search.
  out <- do.call(.plan_parallel, c(args, fmi = .3, start_N = 70))
  ref <- do.call(.plan_parallel, c(args, start_N = ceiling(69 * 7 / 10) + 1))
  N_ref <- ref$value[ref$term == "necessary_N"]
  expect_identical(out$value[out$term == "necessary_N"],
                   ceiling((N_ref - 1) * 10 / 7) + 1)
  expect_identical(attr(out, "necessary_N_complete"), N_ref)
})

test_that("a planned size that is an integer in exact arithmetic is not rounded up a case", {
  # 1 - .3, 1 - .8, and 1 - .9 are inexact in binary, so (N - 1) / (1 - fmi)
  # + 1 can land a rounding error above an exact integer (50 at .3 is 71).
  # The exact answers in integer arithmetic: the smallest N with
  # (N - 1)(1 - fmi) >= N_c - 1.
  n <- 2:400
  map <- function(f) DMAR:::.aipe_ceiling(DMAR:::.aipe_fmi_inflate(n, f))
  expect_identical(map(.3), ceiling((n - 1) * 10 / 7) + 1)
  expect_identical(map(.8), (n - 1) * 5 + 1)
  expect_identical(map(.9), (n - 1) * 10 + 1)
  expect_identical(map(.25), ceiling((n - 1) * 4 / 3) + 1)
  expect_identical(DMAR:::.aipe_ceiling(DMAR:::.aipe_fmi_inflate(50, .3)), 71)
  # A genuine fraction of a case is still rounded up.
  expect_identical(DMAR:::.aipe_ceiling(37.878005), 38)
  expect_identical(DMAR:::.aipe_ceiling(71 + 1e-4), 72)
})

test_that("a true_residual_var vector for the Parallel model stops with the function's own message", {
  # The length was checked after the comparison true_residual_var <= 0, which
  # stops on a vector with R's "the condition has length > 1".
  expect_error(ss_aipe_reliability(model = "Parallel", type = "Normal Theory",
                                   width = .1, n_items = 4, true_item_covariance = .3,
                                   true_residual_var = c(.2, .3, .2, .3)),
               "'true_residual_var' must be one number for the Parallel model",
               fixed = TRUE)
  expect_error(ss_aipe_reliability(model = "Parallel", type = "Normal Theory",
                                   width = .1, n_items = 4, true_factor_loading = .6,
                                   true_residual_var = c(.2, .3, .2, .3)),
               "'true_residual_var' must be one number for the Parallel model",
               fixed = TRUE)
  # A single nonpositive value still gets its own message.
  expect_error(ss_aipe_reliability(model = "Parallel", type = "Normal Theory",
                                   width = .1, n_items = 4, true_item_covariance = .3,
                                   true_residual_var = 0),
               "must be greater than zero", fixed = TRUE)
  # The same order for n_items in the True Score model planned from
  # true_factor_loading: it was compared with the length of true_residual_var before it
  # was checked for presence, and a missing value stopped with "argument is
  # of length zero".
  expect_error(ss_aipe_reliability(model = "True Score", type = "Normal Theory",
                                   width = .1, true_factor_loading = .6,
                                   true_residual_var = c(.2, .3, .2, .3)),
               "You need to enter a quantity for n_items", fixed = TRUE)
})

test_that("a true_factor_loading vector for the Parallel or True Score model stops with a message naming the Congeneric model", {
  # The Parallel model planned from true_factor_loading passed the vector to
  # matrix(true_factor_loading, nrow = 1, ncol = n_items), so four loadings planned the
  # congeneric pattern they describe under the Parallel label: here
  # N = 96 from the alpha of that pattern, without a word.
  for (m in c("Parallel", "True Score")) {
    expect_error(ss_aipe_reliability(model = m, type = "Normal Theory",
                                     width = .1, n_items = 4,
                                     true_factor_loading = c(.8, .7, .7, .8),
                                     true_residual_var = if (m == "Parallel") .4
                                     else c(.4, .5, .5, .4)),
                 paste0("'true_factor_loading' must be one number for the ", m,
                        " model. If you want to enter multiple 'true_factor_loading' ",
                        "values, you must use the Congeneric model."),
                 fixed = TRUE)
  }
  # An empty vector is not one number either.
  expect_error(ss_aipe_reliability(model = "Parallel", type = "Normal Theory",
                                   width = .1, n_items = 4, true_factor_loading = numeric(0),
                                   true_residual_var = .4),
               "'true_factor_loading' must be one number for the Parallel model",
               fixed = TRUE)
  # One loading still plans the Parallel model, as before, and gives the
  # plan of true_item_covariance = .49, from which the function takes the loading
  # sqrt(.49) = .7 and so the same covariance matrix.
  from_lambda <- ss_aipe_reliability(model = "Parallel",
                                     type = "Normal Theory", width = .1,
                                     n_items = 4, true_factor_loading = .7, true_residual_var = .51)
  from_cor <- ss_aipe_reliability(model = "Parallel", type = "Normal Theory",
                                  width = .1, n_items = 4, true_item_covariance = .49,
                                  true_residual_var = .51)
  expect_identical(from_lambda$value[from_lambda$term == "necessary_N"],
                   from_cor$value[from_cor$term == "necessary_N"])
})

test_that("the assurance path carries the planner class, so tidy() and glance() agree with the closed form", {
  skip_on_cran()
  # Before 1.1.0 only the closed form returned class dmar_ss_aipe; the
  # assurance table fell through to the generic dmar_tbl methods.
  closed <- .plan_parallel()
  assured <- .plan_parallel(assurance = .80, G_initial = 20,
                            G = 40, seed = 7)
  expect_identical(class(assured), class(closed))
  expect_s3_class(assured, "dmar_ss_aipe")
  N <- assured$value[assured$term == "necessary_N"]
  td <- generics::tidy(assured)
  expect_identical(td, data.frame(term = "sample_size", estimate = N,
                                  width = .1, stringsAsFactors = FALSE))
  expect_named(generics::tidy(closed), names(td))
  # The closed form's table holds only the size, so its width is NA, as in
  # 1.0.1 (the help said both forms gave it).
  expect_identical(generics::tidy(closed)$width, NA_real_)
  gl <- generics::glance(assured)
  expect_named(gl, c("term", "estimate", "width", "assurance",
                     "assurance_attained", "G", "no_interval"))
  expect_identical(gl$no_interval, 0L)
  expect_identical(gl$estimate, N)
  expect_identical(gl$assurance, .80)
  expect_identical(gl$G, 40)
  # A plan under missing data keeps the class on both paths.
  expect_s3_class(.plan_parallel(fmi = .3, assurance = .80, G_initial = 20,
                                 G = 40, seed = 7), "dmar_ss_aipe")
  # A plan under missing data says so in its one-row summary too.
  m <- .plan_parallel(fmi = .2)
  gm <- generics::glance(m)
  expect_identical(gm$fmi, .2)
  expect_identical(gm$necessary_N_complete,
                   attr(m, "necessary_N_complete"))
  expect_false("no_interval" %in% names(gm))
})

## The coefficient planned for, the closed forms, and the assurance search
## (1.1.0).

# The van Zyl, Neudecker, and Nel (2000) variance of alpha times N - 1,
# written out from Terry and Kelley (2012, Eqs. 7 and 8) on whatever matrix
# it is given.
.vz <- function(S) {
  J <- ncol(S)
  g <- sum(S)
  J^2 / (J - 1)^2 * 2 / g^3 *
    (g * (sum(diag(S %*% S)) + sum(diag(S))^2) -
       2 * sum(diag(S)) * sum(S %*% S))
}
.pop <- function(lambda, psi) tcrossprod(lambda) + diag(psi, length(lambda))

test_that("the congeneric model with type = \"Normal Theory\" is refused, with or without an assurance", {
  # The closed form returned 53 without reading type, and with an assurance
  # every iteration's refusal was swallowed and the search stopped with
  # "missing value where TRUE/FALSE needed".
  msg <- "The congeneric model has no normal theory interval"
  expect_error(ss_aipe_reliability(model = "Congeneric", type = "Normal Theory",
                                   width = .15, n_items = 4,
                                   true_factor_loading = c(.8, .7, .7, .8),
                                   true_residual_var = c(.4, .5, .5, .4)),
               msg, fixed = TRUE)
  expect_error(ss_aipe_reliability(model = "omega", type = "nt", width = .15,
                                   n_items = 4, true_factor_loading = c(.8, .7, .7, .8),
                                   true_residual_var = c(.4, .5, .5, .4),
                                   assurance = .8, G_initial = 20,
                                   G = 40, seed = 1),
               msg, fixed = TRUE)
  # Two items do not identify the congeneric model.
  expect_error(ss_aipe_reliability(model = "Congeneric",
                                   type = "Factor Analytic", width = .1,
                                   n_items = 2, true_factor_loading = c(.8, .7),
                                   true_residual_var = c(.4, .5)),
               "The congeneric model needs at least three items", fixed = TRUE)
})

test_that("the tau equivalent plan is for alpha of the covariance matrix, not of its correlation matrix", {
  # Finding of the 1.1.0 review: equal loadings .7 and error variances .1,
  # .2, 1, 2 (item standard deviations .77 to 1.58). Eq. 19 on the
  # covariance matrix needs 367 cases for a width of .10; 1.0.x evaluated
  # it on the correlation matrix and returned 206, at which the interval
  # for raw alpha covers about 85 percent of the time.
  Sig <- .pop(rep(.7, 4), c(.1, .2, 1, 2))
  z <- stats::qnorm(.975)
  expected <- ceiling(.vz(Sig) / (.10 / (2 * z))^2 + 1)
  expect_identical(expected, 367)
  plan <- function(...) {
    ss_aipe_reliability(model = "True Score", type = "Normal Theory",
                        width = .10, ...)$value
  }
  expect_identical(plan(true_cov_matrix = Sig), expected)
  expect_identical(plan(n_items = 4, true_factor_loading = .7, true_residual_var = c(.1, .2, 1, 2)),
                   expected)
  expect_identical(plan(n_items = 4, true_item_covariance = .49, true_residual_var = c(.1, .2, 1, 2)),
                   expected)
  # A common rescaling of the items leaves alpha and the plan as they are;
  # a correlation matrix is the population of standardized items and plans
  # for them.
  expect_identical(plan(true_cov_matrix = 9 * Sig), expected)
  expect_identical(plan(true_cov_matrix = stats::cov2cor(Sig)), 206)
  # With equal item variances nothing changes from 1.0.x.
  expect_identical(plan(n_items = 6, true_item_covariance = .3, true_residual_var = rep(.7, 6)), 291)
})

test_that("the tau equivalent closed form for the factor analytic interval is that interval's own", {
  # Finding of the 1.1.0 review: with type = "Factor Analytic" the tau
  # equivalent model returned the normal theory size, 367 for loadings .7
  # and error variances .1, .2, 1, 2 at a width of .10, a plan for a wider
  # interval than the one asked for. 1.0.x, which evaluated that form on
  # the correlation matrix, returned 206, too few for the factor analytic
  # interval of these items (mean width .108). The variance of the model
  # implied coefficient with equal loadings is V = 0.156501199 from
  # reliability_alpha(estimator = "model_implied", ci_method = "ml")
  # (lavaan 0.7-3) at N = 1e6, and the plan is 242.
  lam <- rep(.7, 4)
  psi <- c(.1, .2, 1, 2)
  expect_equal(DMAR:::.ss_aipe_omega_avar(NULL, lam, psi,
                                          equal_loadings = TRUE),
               0.156501199, tolerance = 1e-7)
  plan <- function(type, ...) {
    ss_aipe_reliability(model = "True Score", type = type, width = .10,
                        ...)$value
  }
  expect_identical(plan("Factor Analytic", n_items = 4, true_factor_loading = .7,
                        true_residual_var = psi), 242)
  expect_identical(plan("Factor Analytic", n_items = 4, true_item_covariance = .49,
                        true_residual_var = psi), 242)
  expect_identical(plan("Normal Theory", n_items = 4, true_factor_loading = .7,
                        true_residual_var = psi), 367)
  # With equal error variances the two variances coincide, so the parallel
  # model keeps one closed form for both types: the compound symmetry form
  # of van Zyl, Neudecker, and Nel (2000) is the variance of the model
  # implied coefficient with equal loadings at a population of parallel
  # items.
  for (J in c(3, 6)) for (l in c(.5, .7)) for (ps in c(.3, 1)) {
    a <- (J * l)^2 / ((J * l)^2 + J * ps)
    expect_equal(DMAR:::.ss_aipe_omega_avar(NULL, rep(l, J), rep(ps, J),
                                            equal_loadings = TRUE),
                 2 * (1 - a)^2 * J / (J - 1), tolerance = 1e-12)
  }
  for (type in c("Factor Analytic", "Normal Theory")) {
    expect_identical(ss_aipe_reliability(model = "Parallel", type = type,
                                         width = .1, n_items = 6, true_item_covariance = .6,
                                         true_residual_var = .4)$value, 38)
  }
})

test_that("the tau equivalent factor analytic plan agrees with lavaan and meets its width by simulation", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  Sig <- .pop(rep(.7, 4), c(.1, .2, 1, 2))
  plan <- function(...) {
    ss_aipe_reliability(model = "True Score", type = "Factor Analytic",
                        width = .10, ...)$value
  }
  # A population given as a matrix is fit by the equal loadings model
  # first, and gives the same plan, as does any common rescaling of it.
  expect_identical(plan(true_cov_matrix = Sig), 242)
  expect_identical(plan(true_cov_matrix = 9 * Sig), 242)
  # Where the items do not load equally, the plan is still for the
  # interval asked for, with the variance lavaan reports for the model
  # implied coefficient of the equal loadings model fit to the population.
  raw <- .pop(c(.3, .5, 1, 2), c(.1, .3, 1, 3))
  r <- suppressWarnings(reliability_alpha(S = raw, N = 1e6,
                                          estimator = "model_implied",
                                          ci_method = "ml"))
  expect_equal(DMAR:::.ss_aipe_omega_avar(raw, equal_loadings = TRUE),
               r$value[r$term == "se"]^2 * 1e6, tolerance = 1e-6)
  # At the planned 242 cases the factor analytic interval averages the
  # target width; at the 367 of the normal theory form it averaged .081 in
  # the review's simulation.
  U <- chol(Sig)
  set.seed(2026)
  w <- replicate(200, {
    S <- stats::var(matrix(stats::rnorm(242 * 4), 242, 4) %*% U)
    ci <- DMAR:::.ss_aipe_reliability_ci(S, 242, "True Score",
                                         "Factor Analytic", .95)
    ci$upper - ci$lower
  })
  expect_lt(abs(mean(w) - .10), .004)
})

test_that("the congeneric closed form is the delta method variance of the model implied omega", {
  # V / N is the squared standard error of reliability_omega(denominator =
  # "model_implied", ci_method = "ml"). Loadings .9, .8, .3, .2 of
  # standardized items: V = 0.256274867 from that function (lavaan 0.7-3)
  # at N = 1e6, and 177 cases for a width of .15, where 1.0.x returned
  # alpha's Eq. 19, 297. Loadings .8, .7, .7, .8: the V of the oracle above.
  lam <- c(.9, .8, .3, .2)
  expect_equal(DMAR:::.ss_aipe_omega_avar(NULL, lam, 1 - lam^2),
               0.256274867, tolerance = 1e-7)
  expect_equal(DMAR:::.ss_aipe_omega_avar(NULL, c(.8, .7, .7, .8),
                                          c(.4, .5, .5, .4)),
               0.0739065022, tolerance = 1e-7)
  plan <- function(...) {
    ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
                        width = .15, ...)$value
  }
  expect_identical(plan(n_items = 4, true_factor_loading = lam, true_residual_var = 1 - lam^2), 177)
  # Unequal item variances are planned for as they are: loadings .3, .5, 1,
  # 2 with error variances .1, .3, 1, 3 need 151 cases; 1.0.x returned 75,
  # the plan for alpha of their correlation matrix.
  expect_identical(plan(n_items = 4, true_factor_loading = c(.3, .5, 1, 2),
                        true_residual_var = c(.1, .3, 1, 3)), 151)
})

test_that("the congeneric closed form agrees with lavaan and gives one plan however the population is given", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  lam <- c(.9, .8, .3, .2)
  Sig <- .pop(lam, 1 - lam^2)
  N <- 1e6
  r <- reliability_omega(S = Sig, N = N, denominator = "model_implied",
                         ci_method = "ml")
  V <- r$value[r$term == "se"]^2 * N
  expect_equal(DMAR:::.ss_aipe_omega_avar(Sig, lam, 1 - lam^2), V,
               tolerance = 1e-6)
  z <- stats::qnorm(.975)
  expect_identical(ceiling(4 * z^2 * V / .15^2 + 1), 177)
  plan <- function(...) {
    ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
                        width = .15, ...)$value
  }
  # The same population given as a matrix is fit by the single-factor model
  # first, and gives the same plan; so does any common rescaling of it.
  expect_identical(plan(true_cov_matrix = Sig), 177)
  expect_identical(plan(true_cov_matrix = 4 * Sig), 177)
  # Three items, the exactly identified case, against lavaan.
  lam3 <- c(.7, .6, .5)
  Sig3 <- .pop(lam3, 1 - lam3^2)
  r3 <- reliability_omega(S = Sig3, N = N, denominator = "model_implied",
                          ci_method = "ml")
  expect_equal(DMAR:::.ss_aipe_omega_avar(Sig3),
               r3$value[r3$term == "se"]^2 * N, tolerance = 1e-6)
  # The correlation matrix of a population with unequal item variances is
  # another population, of standardized items, with another plan.
  raw <- .pop(c(.3, .5, 1, 2), c(.1, .3, 1, 3))
  expect_identical(plan(true_cov_matrix = raw), 151)
  expect_false(plan(true_cov_matrix = stats::cov2cor(raw)) == 151)
  # A population with no proper single-factor solution is refused.
  R3 <- matrix(c(1, .5, .5, .5, 1, -.2, .5, -.2, 1), 3)
  expect_error(plan(true_cov_matrix = R3), "is not a proper solution", fixed = TRUE)
  R4 <- matrix(c(1, .8, .8, .5, .8, 1, .4, .3, .8, .4, 1, .3,
                 .5, .3, .3, 1), 4)
  expect_error(plan(true_cov_matrix = R4), "is not a proper solution", fixed = TRUE)
})

test_that("the congeneric plan meets its expected width by simulation", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # Loadings .3, .5, 1, 2 and error variances .1, .3, 1, 3, width .15. At
  # the planned 151 cases the factor analytic interval averaged .1495 wide
  # in 1,500 replications; at the 75 of 1.0.x it averaged .214.
  Sig <- .pop(c(.3, .5, 1, 2), c(.1, .3, 1, 3))
  N <- ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
                           width = .15, true_cov_matrix = Sig)$value
  U <- chol(Sig)
  set.seed(2026)
  w <- replicate(300, {
    S <- stats::var(matrix(stats::rnorm(N * 4), N, 4) %*% U)
    ci <- DMAR:::.ss_aipe_reliability_ci(S, N, "Congeneric",
                                         "Factor Analytic", .95)
    ci$upper - ci$lower
  })
  # Monte Carlo standard error of the mean width about .0013.
  expect_lt(abs(mean(w) - .15), .006)
})

test_that("a population given by data plans with an assurance, as from its covariance matrix", {
  skip_on_cran()
  # 1.0.x never set the number of items from data, and every assurance plan from data stopped
  # in rep(0, i) with "invalid 'times' argument".
  set.seed(1)
  df <- as.data.frame(matrix(stats::rnorm(800), 200, 4) %*%
                        chol(.pop(rep(.7, 4), rep(.51, 4))))
  # The data route computes var(na.omit(data), na.rm = TRUE), which differs
  # from var(data) in the last bits; the Monte Carlo draws, through the
  # eigenvectors of a matrix with nearly equal eigenvalues, do not forgive
  # even that.
  S <- stats::var(df, na.rm = TRUE)
  for (m in "True Score") {
    from_data <- ss_aipe_reliability(model = m, type = "Normal Theory",
                                     width = .1, data = df, assurance = .8,
                                     G_initial = 20, G = 40,
                                     seed = 1)
    from_S <- ss_aipe_reliability(model = m, type = "Normal Theory",
                                  width = .1, true_cov_matrix = S,
                                  assurance = .8, G_initial = 20,
                                  G = 40, seed = 1)
    expect_identical(from_data, from_S)
  }
  expect_identical(
    ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
                        width = .1, data = df),
    ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
                        width = .1, true_cov_matrix = S))
})

test_that("the assurance search moves below a closed form size that already meets the target", {
  skip_on_cran()
  # With an assurance of .30 the expected width size (291) is more than
  # enough. 1.0.x returned the starting size, 292, whatever the assurance.
  args <- list(model = "Parallel", type = "Normal Theory", width = .1, n_items = 6,
               true_item_covariance = .3, true_residual_var = .7)
  closed <- do.call(ss_aipe_reliability, args)$value
  expect_identical(closed, 291)
  low <- do.call(ss_aipe_reliability,
                 c(args, assurance = .30, G_initial = 200,
                   G = 1000, seed = 1))
  v <- stats::setNames(low$value, low$term)
  expect_lt(v[["necessary_N"]], closed)
  expect_gt(v[["necessary_N"]], 240)
  expect_gte(v[["assurance_attained"]], .30)
  # A search whose starting size falls short draws exactly what 1.0.x drew
  # and returns its answer (values of the 1.0.x function at this call).
  kept <- .plan_parallel(assurance = .85, G_initial = 50, G = 200,
                         seed = 113)
  expect_identical(kept$value[kept$term == "necessary_N"], 57)
  expect_identical(kept$value[kept$term == "assurance_attained"], .875)
})

test_that("a population with an error variance of zero or below, or a true_cov_matrix that is not a covariance matrix, is refused on every route", {
  # Finding of the 1.1.0 review: a congeneric population given as S whose
  # single-factor fit is improper was refused, while the same population
  # given by its loadings was planned without a word (1.0.x: 81 cases for
  # the congeneric model; and 46 for the tau equivalent one with one
  # loading of .7 on all four items and these error variances), and the
  # refusal pointed to that route.
  lam <- c(.7, .6, .5, .5)
  psi <- c(.51, -.1, .75, .75)
  msg <- "Every value of 'true_residual_var' must be greater than zero"
  expect_error(ss_aipe_reliability(model = "Congeneric",
                                   type = "Factor Analytic", width = .15,
                                   n_items = 4, true_factor_loading = lam, true_residual_var = psi),
               msg, fixed = TRUE)
  expect_error(ss_aipe_reliability(model = "Congeneric",
                                   type = "Factor Analytic", width = .15,
                                   n_items = 4, true_item_covariance = .4,
                                   true_residual_var = c(.6, -.1, .6, .6)),
               msg, fixed = TRUE)
  expect_error(ss_aipe_reliability(model = "Congeneric",
                                   type = "Factor Analytic", width = .15,
                                   n_items = 4, true_factor_loading = lam,
                                   true_residual_var = c(.51, 0, .75, .75)),
               msg, fixed = TRUE)
  expect_error(ss_aipe_reliability(model = "True Score",
                                   type = "Factor Analytic", width = .15,
                                   n_items = 4, true_factor_loading = .7, true_residual_var = psi),
               msg, fixed = TRUE)
  expect_error(ss_aipe_reliability(model = "True Score",
                                   type = "Factor Analytic", width = .15,
                                   n_items = 4, true_item_covariance = .49, true_residual_var = psi),
               msg, fixed = TRUE)
  # With the normal theory interval of the tau equivalent model nothing is
  # fit, and the loading and error variances only build the covariance
  # matrix, so the plan is the one the same matrix gives as S (a review of
  # the first repair found the loadings routes refusing what the S route
  # planned). Both matrices here are positive definite.
  for (ps in list(psi, c(0, .5, .5, .5))) {
    from_S <- ss_aipe_reliability(model = "True Score", type = "Normal Theory",
                                  width = .15, true_cov_matrix = .pop(rep(.7, 4), ps))$value
    expect_identical(ss_aipe_reliability(model = "True Score",
                                         type = "Normal Theory", width = .15,
                                         n_items = 4, true_factor_loading = .7,
                                         true_residual_var = ps)$value, from_S)
    expect_identical(ss_aipe_reliability(model = "True Score",
                                         type = "Normal Theory", width = .15,
                                         n_items = 4, true_item_covariance = .49,
                                         true_residual_var = ps)$value, from_S)
  }
  expect_identical(ss_aipe_reliability(model = "True Score",
                                       type = "Normal Theory", width = .15,
                                       n_items = 4, true_factor_loading = .7,
                                       true_residual_var = psi)$value, 74)
  # A matrix that is not a covariance matrix is refused by both routes,
  # before covmat_from_cfa() warns about its determinant.
  expect_error(ss_aipe_reliability(model = "True Score",
                                   type = "Normal Theory", width = .15, n_items = 4,
                                   true_factor_loading = .7,
                                   true_residual_var = c(.51, -.3, .75, .75)),
               paste0("The matrix that 'true_factor_loading' (or 'true_item_covariance') and ",
                      "'true_residual_var' imply is not a covariance matrix: its ",
                      "smallest eigenvalue after standardizing item units is -0.442"), fixed = TRUE)
  expect_error(ss_aipe_reliability(model = "True Score",
                                   type = "Normal Theory", width = .15,
                                   true_cov_matrix = .pop(rep(.7, 4), c(.51, -.3, .75, .75))),
               "'true_cov_matrix' is not a covariance matrix: its smallest eigenvalue after standardizing item units is -0.442",
               fixed = TRUE)
  expect_error(ss_aipe_reliability(model = "True Score",
                                   type = "Normal Theory", width = .15, n_items = 4,
                                   true_factor_loading = .7,
                                   true_residual_var = c(.51, NA, .75, .75)),
               "Every value of 'true_residual_var' must be a finite number.",
               fixed = TRUE)
  # A matrix with a negative eigenvalue (-0.144 here) is no population's
  # covariance matrix; the alpha routes planned from it, 1.0.x included.
  S_bad <- .pop(lam, c(.51, -.3, .75, .75))
  for (m in c("Parallel", "True Score")) {
    expect_error(ss_aipe_reliability(model = m, type = "Normal Theory",
                                     width = .15, true_cov_matrix = S_bad),
                 "'true_cov_matrix' is not a covariance matrix: its smallest eigenvalue after standardizing item units is -1.15",
                 fixed = TRUE)
  }
  expect_error(ss_aipe_reliability(model = "Congeneric",
                                   type = "Factor Analytic", width = .15,
                                   true_cov_matrix = S_bad),
               "'true_cov_matrix' is not a covariance matrix", fixed = TRUE)
  # Fewer than three items with a nonzero loading do not identify the
  # congeneric model; the information matrix is singular, and the plan
  # stopped with R's "system is computationally singular" (1.0.x: 937 and
  # 1078, alpha's closed form).
  ident <- "The congeneric model needs at least three items with a nonzero loading"
  expect_error(ss_aipe_reliability(model = "Congeneric",
                                   type = "Factor Analytic", width = .15,
                                   n_items = 3, true_factor_loading = c(.7, .6, 0),
                                   true_residual_var = c(.51, .64, 1)),
               ident, fixed = TRUE)
  expect_error(ss_aipe_reliability(model = "Congeneric",
                                   type = "Factor Analytic", width = .15,
                                   n_items = 4, true_factor_loading = c(.7, .6, 0, 0),
                                   true_residual_var = c(.51, .64, 1, 1)),
               ident, fixed = TRUE)
  expect_error(ss_aipe_reliability(model = "Congeneric",
                                   type = "Factor Analytic", width = .15,
                                   n_items = 3, true_factor_loading = c(.7, .6, 1e-9),
                                   true_residual_var = c(.51, .64, 1)),
               "The single-factor model is not identified at this population",
               fixed = TRUE)
  # The same populations given as S are refused too, and the refusal no
  # longer points to the lambda route.
  skip_on_cran()
  skip_if_not_installed("lavaan")
  expect_error(ss_aipe_reliability(model = "Congeneric",
                                   type = "Factor Analytic", width = .15,
                                   true_cov_matrix = .pop(lam, psi)),
               "A population to plan for is one the model describes with positive error variances",
               fixed = TRUE)
  expect_error(ss_aipe_reliability(model = "Congeneric",
                                   type = "Factor Analytic", width = .15,
                                   true_cov_matrix = .pop(c(.7, .6, 0, 0),
                                            c(.51, .64, 1, 1))),
               "The single-factor model is not identified at this population",
               fixed = TRUE)
})

test_that("the congeneric closed form is refused only where its delta method is degenerate", {
  # Findings of the 1.1.0 reviews. Loadings that cancel in the unweighted
  # sum, as reverse keyed items left unrecoded do, make omega zero, and the
  # delta method closed form returned N = 1; at an omega of .0038 it
  # returned 49, at which the interval averages .22 wide for a target of .10
  # (1.0.x, with alpha's closed form: 14078 and 12606). A first repair
  # refused every factor analytic plan for a coefficient below half the
  # width, closed form and assurance alike, which withdrew plans that met
  # their target. The closed form is now refused where the reported
  # interval would average wider than the target at the N it gives, under
  # a normal distribution for the sum of the loadings.
  plan <- function(lam, ...) {
    ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
                        width = .10, n_items = length(lam), true_factor_loading = lam,
                        true_residual_var = 1 - lam^2, ...)
  }
  advice <- "reverse score such items before planning"
  expect_error(plan(c(.7, -.7, .5, -.5)), "The sum of the loadings is zero",
               fixed = TRUE)
  expect_error(plan(c(.7, -.7, .5, -.5)), advice, fixed = TRUE)
  expect_error(plan(c(.7, -.7, .5, -.4)),
               paste0("At the N it gives, 49 cases, the sum of the loadings ",
                      "(0.1) would be estimated with a standard error 3.35 ",
                      "times its size"), fixed = TRUE)
  expect_error(plan(c(.7, -.7, .5, -.4), fmi = .2),
               "At the N it gives, 49 cases with complete data,", fixed = TRUE)
  expect_error(plan(c(.7, -.7, .5, -.37)), "standard error 2 times",
               fixed = TRUE)
  # Below the point at which the clamped interval averages the target, the
  # plan is made, with the warning that the loadings differ in sign. At a
  # standard error 1.08 times the sum, the interval averaged .082 wide at
  # these 157 cases in simulation.
  opp <- "items 2 and 4 load in the direction opposite to the other items"
  expect_warning(n <- plan(c(.7, -.7, .5, -.32))$value, opp, fixed = TRUE)
  expect_identical(n, 157)
  expect_warning(n <- plan(c(.7, -.7, .5, -.2))$value, opp, fixed = TRUE)
  expect_identical(n, 419)
  # The point moves with the confidence level: 1.52 times the sum is
  # refused at 95 percent and planned at 99.
  expect_error(plan(c(.7, -.7, .5, -.35)), "standard error 1.52 times",
               fixed = TRUE)
  expect_warning(n <- plan(c(.7, -.7, .5, -.35), conf_level = .99)$value,
                 opp, fixed = TRUE)
  expect_identical(n, 188)
  expect_warning(n <- plan(c(.7, -.7, .5, .5))$value,
                 "item 2 loads in the direction opposite", fixed = TRUE)
  expect_identical(n, 1419)
  # Weak loadings of one sign are planned as in 1.0.x, without a warning,
  # until the sum of the loadings is itself lost in its standard error,
  # where the refusal gives no advice about keying.
  expect_silent(n <- plan(rep(.1, 4))$value)
  expect_identical(n, 3787)
  msg <- tryCatch(plan(rep(.04, 4)), error = conditionMessage)
  expect_match(msg, "standard error 2.01 times", fixed = TRUE)
  expect_no_match(msg, advice, fixed = TRUE)
  # The ratio of the mean clamped width to the planned one, from the normal
  # approximation: .838 at r = 1, one at r = 1.4925, and 1.863 at r = 3.354
  # (95 percent).
  z <- stats::qnorm(.975)
  expect_equal(DMAR:::.ss_aipe_clamped_width_ratio(1, z), .838, tolerance = 1e-3)
  expect_equal(DMAR:::.ss_aipe_clamped_width_ratio(1.4925, z), 1,
               tolerance = 1e-3)
  expect_equal(DMAR:::.ss_aipe_clamped_width_ratio(3.354, z), 1.863,
               tolerance = 1e-3)
  expect_identical(DMAR:::.ss_aipe_clamped_width_ratio(Inf, z), Inf)
  # The parallel and tau equivalent models need no such check, and plan a
  # small coefficient as 1.0.x did: their clamped interval only makes the
  # plan conservative (362 cases averaged .236 wide for a target of .30).
  expect_identical(ss_aipe_reliability(model = "Parallel",
                                       type = "Factor Analytic", width = .3,
                                       n_items = 4, true_item_covariance = .03,
                                       true_residual_var = .97)$value, 362)
  expect_identical(ss_aipe_reliability(model = "True Score",
                                       type = "Factor Analytic", width = .1,
                                       n_items = 4, true_item_covariance = .01,
                                       true_residual_var = rep(.99, 4))$value, 3787)
  skip_on_cran()
  skip_if_not_installed("lavaan")
  lam <- c(.7, -.7, .5, -.5)
  expect_error(ss_aipe_reliability(model = "Congeneric",
                                   type = "Factor Analytic", width = .10,
                                   true_cov_matrix = .pop(lam, 1 - lam^2)),
               "The closed form cannot plan this interval.", fixed = TRUE)
  # Regression test. r was taken at the unrounded closed form, whose N - 1
  # is a fraction of one case for a highly reliable population and a wide
  # target, so four loadings of .99 at a width of .2 were refused (r 4.4,
  # "about 2.39 times as wide") while four of .95 were not, and six
  # loadings near .97 were planned at a width of .1 and refused at .2.
  # It is taken at a realizable N, at least i + 1, and the closed form is
  # floored at the feasibility minimum of i + 1 cases.
  for (l in c(.99, .95)) {
    expect_warning(plan <- ss_aipe_reliability(model = "Congeneric",
                                         type = "Factor Analytic", width = .2,
                                         n_items = 4, true_factor_loading = rep(l, 4),
                                         true_residual_var = rep(1 - l^2, 4)),
                   "only a feasibility")
    expect_identical(plan$value, 5)
  }
  l6 <- c(.97, .95, .99, .96, .98, .97)
  for (w in c(.1, .2)) {
    expect_identical(suppressWarnings(ss_aipe_reliability(
      model = "Congeneric", type = "Factor Analytic", width = w, n_items = 6,
      true_factor_loading = l6, true_residual_var = 1 - l6^2))$value, 7)
  }
  # Loadings that nearly cancel are refused at i + 1 cases, where the
  # interval exists, when the closed form gives fewer.
  l <- c(.7, -.7, .5, -.49)
  expect_error(suppressWarnings(ss_aipe_reliability(
    model = "Congeneric", type = "Factor Analytic", width = .1, n_items = 4,
    true_factor_loading = l, true_residual_var = 1 - l^2)),
    "At 5 cases, the fewest at which the interval can be computed (the closed form gives 2)",
    fixed = TRUE)
  S <- .pop(c(.7, -.7, .5, -.2), 1 - c(.7, -.7, .5, -.2)^2)
  dimnames(S) <- list(paste0("t", 1:4), paste0("t", 1:4))
  expect_warning(n <- ss_aipe_reliability(model = "Congeneric",
                                          type = "Factor Analytic",
                                          width = .10, true_cov_matrix = S)$value,
                 "items t2 and t4 load in the direction opposite",
                 fixed = TRUE)
  expect_identical(n, 419)
})

test_that("an assurance plan is not refused where the congeneric closed form is", {
  # The search counts the widths of the intervals themselves. The interval
  # is mocked, failing below i + 1 cases as the fit does, so the search runs
  # quickly; it starts at 5, from a closed form of 1, and stops at the first
  # size whose widths, 2 / N, are below .10.
  local_mocked_bindings(.ss_aipe_reliability_ci = function(S, N, ...) {
    if (N <= ncol(S)) stop("singular", call. = FALSE)
    list(estimate = .01, se = 1 / N, lower = 0, upper = 2 / N)
  })
  lam <- c(.7, -.7, .5, -.5)
  expect_warning(
    res <- ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
                               width = .10, n_items = 4, true_factor_loading = lam,
                               true_residual_var = 1 - lam^2, assurance = .8,
                               G_initial = 5, G = 5, seed = 1),
    "items 2 and 4 load in the direction opposite", fixed = TRUE)
  expect_identical(res$value[res$term == "necessary_N"], 21)
})

test_that("items that share no common factor are refused with the factor analytic interval of the equal loadings models", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # Four items that covary -.2 have an alpha of -2. 1.0.1 planned 36880
  # cases for the parallel model at a width of .10 and 4099 at .30, and
  # true_item_covariance = 0 gave 4099 at .10 for both models. A first repair refused
  # them with a width of at most "  -4" to plan for.
  S <- matrix(-.2, 4, 4)
  diag(S) <- 1
  msg <- "The items of this population covary zero or less on average (alpha of the population matrix is -2), so"
  expect_error(ss_aipe_reliability(model = "Parallel", type = "Factor Analytic",
                                   width = .1, true_cov_matrix = S), msg, fixed = TRUE)
  expect_error(ss_aipe_reliability(model = "True Score",
                                   type = "Factor Analytic", width = .3, true_cov_matrix = S),
               msg, fixed = TRUE)
  for (m in c("Parallel", "True Score")) {
    psi <- if (m == "Parallel") 1 else rep(1, 4)
    expect_error(ss_aipe_reliability(model = m, type = "Factor Analytic",
                                     width = .1, n_items = 4, true_item_covariance = 0,
                                     true_residual_var = psi),
                 "(alpha of the population matrix is 0)", fixed = TRUE)
  }
  # The normal theory interval is for alpha itself and still plans.
  expect_identical(ss_aipe_reliability(model = "Parallel",
                                       type = "Normal Theory", width = .3,
                                       true_cov_matrix = S)$value, 4099)
  skip_if_not_installed("lavaan")
  # Regression test. Alpha decided for a population given as S, which is
  # right only for items of equal variance: the fit weighs each covariance
  # by the inverse of its items' variances. Three items of variance 1 that
  # correlate .3 and a fourth of variance 4 that correlates -.25 have alpha
  # -0.276 and a fitted coefficient of .315, and were refused (1.0.1
  # planned 786 from alpha of the correlation matrix); the tau equivalent
  # model plans them from the fitted coefficient, and the parallel model,
  # whose closed form is in alpha, refuses them and names it.
  R <- matrix(.3, 4, 4)
  R[4, 1:3] <- R[1:3, 4] <- -.25
  diag(R) <- 1
  S3 <- R * tcrossprod(c(1, 1, 1, 2))
  expect_identical(ss_aipe_reliability(model = "True Score",
                                       type = "Factor Analytic", width = .2,
                                       true_cov_matrix = S3)$value, 244)
  expect_error(ss_aipe_reliability(model = "Parallel",
                                   type = "Factor Analytic", width = .2,
                                   true_cov_matrix = S3),
               "not Parallel.*True Score.*Congeneric")
  # Two items of variance 10 that covary 4 and two of variance 1 that
  # covary -.3 have alpha .336 and a fitted loading of zero: the tau
  # equivalent model stopped with the message for an unidentified model
  # and the parallel model planned 202 cases.
  S2 <- matrix(0, 4, 4)
  S2[1:2, 1:2] <- matrix(c(10, 4, 4, 10), 2)
  S2[3:4, 3:4] <- matrix(c(1, -.3, -.3, 1), 2)
  for (m in "True Score") {
    err <- tryCatch(ss_aipe_reliability(model = m, type = "Factor Analytic",
                                        width = .3, true_cov_matrix = S2),
                    error = conditionMessage)
    expect_match(err, paste0("The single-factor model with equal loadings, ",
                             "fit to this population by maximum likelihood, ",
                             "gives its items a loading of zero"),
                 fixed = TRUE)
    expect_no_match(err, "not identified", fixed = TRUE)
  }
})

test_that("the factor analytic search starts at no fewer than i + 1 cases", {
  # Loadings .95 on four items and a width of .2 need two cases by the
  # closed form, so the search started at three, where the sample
  # covariance matrix of four items is singular, no single-factor model can
  # be fit, and the search stopped. The interval is mocked so that a sample
  # of n <= n_items fails, as the fit does.
  sizes <- new.env()
  sizes$n <- integer(0)
  local_mocked_bindings(.ss_aipe_reliability_ci = function(S, N, ...) {
    sizes$n <- c(sizes$n, N)
    if (N <= ncol(S)) stop("singular", call. = FALSE)
    list(estimate = .97, se = 1 / N, lower = .97 - 1 / N, upper = .97 + 1 / N)
  })
  plan <- ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic",
                              width = .2, n_items = 4, true_factor_loading = rep(.95, 4),
                              true_residual_var = rep(1 - .95^2, 4),
                              assurance = .8, G_initial = 5,
                              G = 5, seed = 1)
  expect_identical(min(sizes$n), 5)
  expect_identical(plan$value[plan$term == "necessary_N"], 11)
})

test_that("the search's intervals are the ones the estimation functions return", {
  # Finding of the 1.1.0 review: the interval behind the assurance search
  # raised the lower limit to 0 and lowered the upper limit to 1, each on
  # its own, so a negative alpha whose interval lies wholly below zero came
  # back inverted, [0, -0.750] here, where reliability_alpha() reports
  # [0, 0].
  S <- matrix(-.2, 4, 4)
  diag(S) <- 1
  ra <- suppressWarnings(reliability_alpha(S = S, N = 60, ci_method = "ml"))
  expect_identical(ra$value[ra$term %in% c("ci_lower", "ci_upper")],
                   c(0, 0))
  for (m in "True Score") {
    ci <- DMAR:::.ss_aipe_reliability_ci(S, 60, m, "Normal Theory", .95)
    expect_lt(ci$estimate, 0)
    expect_identical(c(ci$lower, ci$upper), c(0, 0))
  }
  # A sample whose single-factor solution has run away (error variances
  # summing below zero, a model implied omega above one) is refused by
  # reliability_omega(), and the search now sets it aside too. It was kept,
  # with an estimate of 104 and limits 104 and 1, a width of -103 counted as
  # a narrow interval. Its likelihood has no maximum (the error variance of
  # y1 falls without bound, -565 where lavaan stopped), so the refusal is
  # the one for that (.lavaan_no_maximum()), which comes before the one for
  # error variances that sum below zero.
  skip_on_cran()
  skip_if_not_installed("lavaan")
  lam <- c(.8, .3, .2, .2)
  N <- 288
  set.seed(2)
  X <- matrix(stats::rnorm(N * 4), N, 4) %*% chol(.pop(lam, 1 - lam^2))
  expect_error(suppressWarnings(
    reliability_omega(S = stats::var(X), N = N, denominator = "model_implied",
                      ci_method = "ml")), class = "dmar_no_maximum")
  expect_error(suppressWarnings(
    DMAR:::.ss_aipe_reliability_ci(stats::var(X), N, "Congeneric",
                                   "Factor Analytic", .95)),
    "The single-factor solution has run away", fixed = TRUE)
  # An improper solution whose error variances still sum above zero is
  # returned by reliability_omega() with a warning, and the search keeps it
  # with the same interval.
  set.seed(1)
  X <- matrix(stats::rnorm(78 * 4), 78, 4) %*% chol(.pop(lam, 1 - lam^2))
  ci <- DMAR:::.ss_aipe_reliability_ci(stats::var(X), 78, "Congeneric",
                                       "Factor Analytic", .95)
  expect_warning(
    ro <- reliability_omega(S = stats::var(X), N = 78,
                            denominator = "model_implied", ci_method = "ml"),
    "Improper solution (Heywood case)", fixed = TRUE)
  expect_equal(c(ci$estimate, ci$lower, ci$upper),
               ro$value[match(c("estimate", "ci_lower", "ci_upper"),
                              ro$term)], tolerance = 1e-12)
})

test_that("a Monte Carlo sample that returns no interval counts as one wider than the target", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # Finding of the 1.1.0 review: the assurance was the proportion of narrow
  # intervals among the samples that returned one, so a plan at which a
  # third of the studies would yield no coefficient could meet it (1.0.x
  # set a sample whose fit did not converge aside in the same way). The
  # assurance is the probability of an interval no wider than the target
  # (Terry & Kelley, 2012, p. 373), and a study with no interval has none.
  # The interval is mocked: below 45 cases every other sample returns
  # none, and every interval returned is narrow. 1.0.x and the earlier
  # 1.1.0 search, which took the proportion over the rest, met the target
  # at the starting 39 and stepped down to 7.
  calls <- new.env()
  calls$k <- 0L
  local_mocked_bindings(.ss_aipe_reliability_ci = function(S, N, ...) {
    calls$k <- calls$k + 1L
    if (N < 45 && calls$k %% 2L == 0L) stop("no interval", call. = FALSE)
    list(estimate = .5, se = .001, lower = .5, upper = .51)
  })
  res <- .plan_parallel(assurance = .8, G_initial = 10, G = 20,
                        seed = 1)
  expect_identical(res$value[res$term == "necessary_N"], 45)
  expect_identical(res$value[res$term == "assurance_attained"], 1)
  expect_identical(attr(res, "no_interval"), 0L)
  expect_false(any(grepl("returned no interval",
                         utils::capture.output(print(res)), fixed = TRUE)))
  # One sample in ten returns none at every size: the target of .8 is met
  # at the smallest size, with an assurance of .9 rather than 1, and the
  # two samples of the final pass with no interval are named under the
  # table.
  calls$k <- 0L
  local_mocked_bindings(.ss_aipe_reliability_ci = function(S, N, ...) {
    calls$k <- calls$k + 1L
    if (calls$k %% 10L == 0L) stop("no interval", call. = FALSE)
    list(estimate = .5, se = .001, lower = .5, upper = .51)
  })
  res <- .plan_parallel(assurance = .8, G_initial = 10, G = 20,
                        seed = 1)
  expect_identical(res$value[res$term == "necessary_N"], 7)
  expect_identical(res$value[res$term == "assurance_attained"], .9)
  expect_identical(attr(res, "no_interval"), 2L)
  # The one-row summary carries the count, which glance() left out.
  expect_identical(generics::glance(res)$no_interval, 2L)
  expect_true(any(grepl(paste0(
    "Of the 20 samples of the final Monte Carlo pass, 2 returned no ",
    "interval"), utils::capture.output(print(res)), fixed = TRUE)))
})

test_that("a size at which no Monte Carlo sample returns an interval stops with the first error", {
  local_mocked_bindings(.ss_aipe_reliability_ci = function(...) {
    stop("no interval at this size", call. = FALSE)
  })
  expect_error(.plan_parallel(assurance = .8, G_initial = 5, G = 5,
                              seed = 1),
               paste0("None of the Monte Carlo samples of size 39 returned a ",
                      "confidence interval, so the assurance cannot be ",
                      "estimated. The first failure: no interval at this ",
                      "size"),
               fixed = TRUE)
})
