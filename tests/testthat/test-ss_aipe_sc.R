## assurance_method (Ken's decision of 2026-11-13): only "analytic" is
## offered, the smallest plan, at or above the plan without assurance,
## whose exact probability attains the request.
test_that("ss_aipe_sc() plans the smallest sample size that attains the request and checks assurance_method", {
  plan <- function(...) ss_aipe_sc(true_std_contrast = 0, c_weights = c(1, -1),
                                   width = .8, assurance = .95, ...)
  expect_identical(formals(ss_aipe_sc)$assurance_method, "analytic")
  # The published step of DMAR 1.0.x gave 49, which attains .949878; 50
  # attain .994636 (the exact law coded apart from the package, evidence
  # j1_sc/05; also the values of ss_aipe_smd() for this two-group design).
  def <- plan()
  expect_identical(def$term, c("necessary_n_per_group", "assurance",
                               "assurance_attained"))
  expect_identical(def$value[1:2], c(50, .95))
  expect_equal(def$value[3], .994636, tolerance = 1e-5)
  expect_identical(attr(def, "assurance_method"), "analytic")
  expect_identical(attr(def, "assurance_probability"), "exact")
  expect_identical(plan(assurance_method = "analytic"), def)
  # "simulation" is refused, with or without assurance, in this function's
  # name, and so is anything else.
  for (a in list(.95, NULL)) {
    err <- tryCatch(ss_aipe_sc(true_std_contrast = 0, c_weights = c(1, -1),
                               width = .8, assurance = a,
                               assurance_method = "simulation"),
                    error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "assurance_method")
    expect_match(conditionMessage(err),
                 "is not offered by ss_aipe_sc(): the analytic probability",
                 fixed = TRUE)
  }
  for (bad in list(NA, "anal", c("analytic", "x"), 1,
                   c("analytic", "simulation"))) {
    err <- tryCatch(plan(assurance_method = bad), error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "assurance_method")
    expect_match(conditionMessage(err), "the only method ss_aipe_sc() offers",
                 fixed = TRUE)
  }
  # Where the published plan already attains the request and nothing
  # smaller does, the plan is the same, with its probability.
  keep <- ss_aipe_sc(true_std_contrast = -1, c_weights = c(.5, .5, -1), width = .6,
                     conf_level = NULL, alpha_lower = .005, alpha_upper = .095,
                     assurance = .8)
  expect_identical(keep$value[1], 73)
  expect_equal(keep$value[3], .8865, tolerance = 1e-3)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_aipe_sc() returns a tidy single-row data frame with the planned sample size", {
  res <- ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("term", "value"))
  expect_equal(res$term, "necessary_n_per_group")
  expect_true(is.numeric(res$value) && res$value > 0)
})

test_that("ss_aipe_sc() default conf_level = .95 produces a sensible per-group N", {
  res <- ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4)
  # Documented value: a moderately large per-group sample is needed for a
  # standardized contrast of .6 to be estimated to within .4 at 95%.
  expect_gt(res$value, 50)
  expect_lt(res$value, 500)
})

test_that("ss_aipe_sc() assurance path no longer errors on inner ci_nc_t calls", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # Regression test: the documented example used to fail with
  # "Specify either 'conf_level' or both of 'alpha_lower' and 'alpha_upper'"
  # because the inner ci_nc_t() calls did not pass conf_level = NULL.
  res <- ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0),
                    width = .4, assurance = .90)
  expect_s3_class(res, "data.frame")
  expect_equal(res$term, c("necessary_n_per_group", "assurance",
                           "assurance_attained"))
  # Assurance must (weakly) increase the planned N versus the expected-width
  # plan, because the constraint is more demanding.
  res_noassur <- ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0),
                            width = .4)
  expect_gte(res$value[1], res_noassur$value)
})

test_that("ss_aipe_sc() rejects illegal contrast weights", {
  expect_error(
    ss_aipe_sc(true_std_contrast = .5, c_weights = c(.5, .5, -.5), width = .3),
    "should equal zero"
  )
  expect_error(
    ss_aipe_sc(true_std_contrast = .5, c_weights = c(2, 1, -3), width = .3),
    "fractions"
  )
})

test_that("ss_aipe_sc() accepts explicit alpha_lower/alpha_upper for symmetric CI", {
  # Explicit symmetric alphas must reproduce the default conf_level = .95 plan.
  ref <- ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4)
  alt <- ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
                    conf_level = NULL,
                    alpha_lower = .025, alpha_upper = .025)
  expect_equal(alt$value, ref$value)
})

test_that("ss_aipe_sc() accepts asymmetric alpha_lower/alpha_upper", {
  # Asymmetric split with the same total alpha (here .05) generally requires
  # a different N. Just check the function runs and returns a sensible value.
  res <- ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
                    conf_level = NULL,
                    alpha_lower = .01, alpha_upper = .04)
  expect_s3_class(res, "data.frame")
  expect_equal(res$term, "necessary_n_per_group")
  expect_gt(res$value, 0)
})

test_that("ss_aipe_sc() rejects mixing conf_level and alphas", {
  expect_error(
    ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
               conf_level = .95, alpha_lower = .025, alpha_upper = .025),
    "cannot mix them"
  )
})

test_that("ss_aipe_sc() rejects single-tail alpha (need both)", {
  expect_error(
    ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
               conf_level = NULL, alpha_lower = .025),
    "Supply both"
  )
})

test_that("ss_aipe_sc() rejects out-of-range alpha bounds", {
  expect_error(
    ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
               conf_level = NULL, alpha_lower = -.01, alpha_upper = .05),
    "strictly positive"
  )
  expect_error(
    ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0), width = .4,
               conf_level = NULL, alpha_lower = .5, alpha_upper = .5),
    "sum to less than 1"
  )
})

test_that("ss_aipe_sc() rejects out-of-range assurance values", {
  expect_error(
    ss_aipe_sc(true_std_contrast = .5, c_weights = c(.5, .5, -.5, -.5), width = .3,
               assurance = 1.1),
    "assurance"
  )
  expect_error(
    ss_aipe_sc(true_std_contrast = .5, c_weights = c(.5, .5, -.5, -.5), width = .3,
               assurance = .3),
    "larger than 0\\.5"
  )
})

## The weights rule of ci_sc() (register finding 9).
test_that("ss_aipe_sc() refuses the weights ci_sc() refuses, with its message", {
  for (cw in list(c(.25, .25, -.5), c(.5, .5, -.5), c(1, 1, -2), c(0, 0, 0))) {
    caught <- tryCatch(ss_aipe_sc(true_std_contrast = .5, c_weights = cw,
                                  width = .4), error = function(e) e)
    by_ci_sc <- tryCatch(ci_sc(means = rep(0, length(cw)), s_anova = 1,
                               c_weights = cw, n_per_group = 10),
                         error = function(e) e)
    expect_s3_class(caught, "dmar_argument_invalid")
    expect_identical(caught$arg, "c_weights")
    expect_identical(conditionMessage(caught), conditionMessage(by_ci_sc))
  }
  # The weights ci_sc() accepts are accepted, and integer weights are
  # judged by their values.
  expect_identical(
    ss_aipe_sc(true_std_contrast = .5, c_weights = -1:1, width = .5)$value,
    ss_aipe_sc(true_std_contrast = .5, c_weights = c(-1, 0, 1), width = .5)$value)
  expect_silent(ss_aipe_sc(true_std_contrast = .5, c_weights = c(1/3, 1/3, 1/3, -1),
                           width = .8))
})

## The sign with unequal tails (register finding 2): a negative contrast
## gets the positive contrast's plan with the tails exchanged. DMAR 1.0.x
## gave 522 and 342 at the two cases below (evidence i1_sc_contrasts/06).
test_that("ss_aipe_sc() with assurance plans a negative contrast as its magnitude with the tails exchanged", {
  plan <- function(d, al, au, w = .4, a = .9, cw = c(1, -1))
    ss_aipe_sc(true_std_contrast = d, c_weights = cw, width = w,
               conf_level = NULL, alpha_lower = al, alpha_upper = au,
               assurance = a)$value[1]
  # A small case (6418c4f gave 72, the plan for the tails as given; the
  # exact probabilities are in the next test).
  three <- c(.5, .5, -1)
  expect_identical(plan(-1, .005, .095, w = .6, a = .8, cw = three), 73)
  expect_identical(plan(-1, .005, .095, w = .6, a = .8, cw = three),
                   plan(1, .095, .005, w = .6, a = .8, cw = three))
  expect_identical(plan(1, .005, .095, w = .6, a = .8, cw = three), 72)
  skip_on_cran()
  # The smallest plans that attain the request (the published step of the
  # first 1.1.0 builds, raised where it fell short, gave 524 and 347).
  expect_identical(plan(-2.5, .001, .049), 523)
  expect_identical(plan(-2.5, .001, .049), plan(2.5, .049, .001))
  expect_identical(plan(-3, .001, .099, w = .5), 346)
  expect_identical(plan(-3, .001, .099, w = .5), plan(3, .099, .001, w = .5))
  # A positive contrast, and either sign with equal tails, plan as for the
  # magnitude with the tails as given.
  expect_identical(plan(2.5, .001, .049), 520)
  expect_identical(plan(-.8, .025, .025, w = .6), plan(.8, .025, .025, w = .6))
  # Without assurance the signed noncentrality already gave this.
  no_a <- function(d, al, au)
    ss_aipe_sc(true_std_contrast = d, c_weights = c(1, -1), width = .4,
               conf_level = NULL, alpha_lower = al, alpha_upper = au)$value
  expect_identical(no_a(-2.5, .001, .049), no_a(2.5, .049, .001))
})

# The exact P(ci_sc() width <= w): {t : width <= w} = [t_lo, t_hi] (the
# width falls then rises in t), and P(t_lo <= T <= t_hi) is integrated
# over the chi square, E_V[Phi(t_hi S - lambda) - Phi(t_lo S - lambda)].
.sc_exact_width_probability <- function(n, d, cw, al, au, w) {
  df <- length(cw) * (n - 1)
  lam <- d / sqrt(sum(cw^2) / n)
  width_at <- function(t) {
    ci <- ci_sc(ncp = t, c_weights = cw, n_per_group = n, s_anova = 1,
                conf_level = NULL, alpha_lower = al, alpha_upper = au)
    ci$value[3] - ci$value[1]
  }
  m <- stats::optimize(width_at, c(-abs(lam) - 20, abs(lam) + 20), tol = 1e-11)
  hi <- m$minimum + 1
  while (width_at(hi) <= w) hi <- hi + 2 * abs(hi) + 1
  lo <- m$minimum - 1
  while (width_at(lo) <= w) lo <- lo - 2 * abs(lo) - 1
  t_lo <- stats::uniroot(function(t) width_at(t) - w, c(lo, m$minimum), tol = 1e-12)$root
  t_hi <- stats::uniroot(function(t) width_at(t) - w, c(m$minimum, hi), tol = 1e-12)$root
  g <- function(u) {
    s <- sqrt(stats::qchisq(u, df) / df)
    stats::pnorm(t_hi * s - lam) - stats::pnorm(t_lo * s - lam)
  }
  stats::integrate(g, 0, 1, rel.tol = 1e-12, subdivisions = 5000L)$value
}

test_that("ss_aipe_sc() plans with unequal tails attain the assurance for both signs (exact probability)", {
  skip_on_cran()  # each probability inverts ci_sc() widths and integrates
  cw <- c(1, -1)
  # Values of the exact probability recorded in the lane evidence
  # (work/round6/evidence/i1_sc_contrasts/07_exact_p_old_plans.log): .9123
  # at 524, .8899 at 522 (DMAR 1.0.x), and .8636 at 520 (the first 1.1.0
  # builds). The plan is 523, where it is .9016.
  p_524 <- .sc_exact_width_probability(524, -2.5, cw, .001, .049, .4)
  expect_equal(p_524, .9123, tolerance = 1e-3)
  p_523 <- .sc_exact_width_probability(523, -2.5, cw, .001, .049, .4)
  expect_equal(p_523, .9016, tolerance = 1e-3)
  expect_gte(p_523, .90)
  expect_equal(.sc_exact_width_probability(523, 2.5, cw, .049, .001, .4), p_523,
               tolerance = 1e-8)
  expect_lt(.sc_exact_width_probability(522, -2.5, cw, .001, .049, .4), .90)
  expect_lt(.sc_exact_width_probability(520, -2.5, cw, .001, .049, .4), .87)
  # The small case of the test above: .8865 at 73, .7928 at 72.
  three <- c(.5, .5, -1)
  expect_equal(.sc_exact_width_probability(73, -1, three, .005, .095, .6), .8865,
               tolerance = 1e-3)
  expect_lt(.sc_exact_width_probability(72, -1, three, .005, .095, .6), .80)
  # Both signs with unequal tails at assurance .80 and .99 meet the request.
  for (k in list(c(1.5, .01, .04, .8), c(.8, .001, .049, .99), c(2.5, .001, .049, .8))) {
    for (sgn in c(-1, 1)) {
      d <- sgn * k[1]
      n <- ss_aipe_sc(true_std_contrast = d, c_weights = cw, width = .6,
                      conf_level = NULL, alpha_lower = k[2], alpha_upper = k[3],
                      assurance = k[4])$value[1]
      expect_gte(.sc_exact_width_probability(n, d, cw, k[2], k[3], .6), k[4])
    }
  }
})

test_that("ss_aipe_sc() with unequal tails and a negative contrast attains its assurance in raw two-group data", {
  skip_on_cran()  # Monte Carlo
  # Two groups of 523 at -2.5, tails .001/.049: the share of 2000 data sets
  # (all usable) no wider than .4 against the exact .9016.
  set.seed(5240)
  n <- 523
  R <- 2000
  width <- vapply(seq_len(R), function(r) {
    y1 <- stats::rnorm(n, -2.5, 1)
    y2 <- stats::rnorm(n, 0, 1)
    s <- sqrt((stats::var(y1) + stats::var(y2)) / 2)
    ci <- ci_sc(means = c(mean(y1), mean(y2)), s_anova = s, c_weights = c(1, -1),
                n_per_group = n, conf_level = NULL, alpha_lower = .001,
                alpha_upper = .049)
    ci$value[3] - ci$value[1]
  }, numeric(1))
  expect_true(all(is.finite(width)))
  share <- mean(width <= .4)
  mcse <- sqrt(share * (1 - share) / R)
  expect_lt(abs(share - .9016) / mcse, 4)
  expect_gt(share, .9 - 3 * mcse)
})

test_that("ss_aipe_sc() and ss_aipe_sc_sensitivity() plan with assurance without a warning", {
  skip_on_cran()  # plans with assurance and a short simulation
  expect_no_warning(expect_identical(
    ss_aipe_sc(true_std_contrast = -2.5, c_weights = c(1, -1), width = .4,
               conf_level = NULL, alpha_lower = .001, alpha_upper = .049,
               assurance = .9)$value[1], 523))
  # 102 per group here, the published plan and the smallest that attains
  # the request.
  expect_no_warning(expect_identical(
    ss_aipe_sc(true_std_contrast = .6, c_weights = c(.5, .5, -.5, -.5, 0),
               width = .4, assurance = .9)$value[1], 102))
  expect_no_warning(ss_aipe_sc(true_std_contrast = -2.5, c_weights = c(1, -1),
                               width = .4, conf_level = NULL, alpha_lower = .001,
                               alpha_upper = .049))
  expect_no_warning(ss_aipe_sc_sensitivity(
    true_std_contrast = 0, supposed_std_contrast = 0, c_weights = c(1, -1),
    width = .8, assurance = .95, G = 2, print_iter = FALSE, seed = 1))
})

test_that("ss_aipe_sc() reports the exact probability (separate integral) at the smallest plan", {
  skip_on_cran()  # each probability inverts ci_sc() widths and integrates
  cases <- list(
    list(d = 0, cw = c(1, -1), w = 1.6, al = .025, au = .025, a = .95),
    list(d = -2.5, cw = c(1, -1), w = .6, al = .01, au = .04, a = .9),
    list(d = 2.5, cw = c(1, -1), w = .6, al = .04, au = .01, a = .9),
    list(d = .6, cw = c(.5, .5, -.5, -.5, 0), w = .4, al = .025, au = .025, a = .9))
  for (k in cases) {
    n0 <- ss_aipe_sc(true_std_contrast = k$d, c_weights = k$cw, width = k$w,
                     conf_level = NULL, alpha_lower = k$al,
                     alpha_upper = k$au)$value
    ex <- ss_aipe_sc(true_std_contrast = k$d, c_weights = k$cw, width = k$w,
                     conf_level = NULL, alpha_lower = k$al, alpha_upper = k$au,
                     assurance = k$a)$value
    p <- .sc_exact_width_probability(ex[1], k$d, k$cw, k$al, k$au, k$w)
    expect_equal(ex[3], p, tolerance = 1e-7)
    expect_gte(ex[3], k$a)
    expect_gte(ex[1], n0)
    if (ex[1] > n0)
      expect_lt(.sc_exact_width_probability(ex[1] - 1, k$d, k$cw, k$al, k$au, k$w), k$a)
  }
  # The two reflected cases get one plan and one probability.
  a <- ss_aipe_sc(-2.5, c(1, -1), .6, conf_level = NULL, alpha_lower = .01,
                  alpha_upper = .04, assurance = .9)
  b <- ss_aipe_sc(2.5, c(1, -1), .6, conf_level = NULL, alpha_lower = .04,
                  alpha_upper = .01, assurance = .9)
  expect_equal(a$value, b$value, tolerance = 1e-9)
})

test_that("ss_aipe_sc() returns the smallest plan at or above the plan without assurance, by an exhaustive scan", {
  skip_on_cran()  # each probability inverts ci_sc() widths and integrates
  # The scan raises n one at a time from the plan without assurance and
  # stops at the first n whose probability, by the integral above (coded
  # apart from the package), reaches the request.
  cases <- list(
    list(d = 0, cw = c(1, -1), w = .8, al = .025, au = .025, a = .95, n = 50),
    list(d = .6, cw = c(.5, .5, -.5, -.5, 0), w = .4, al = .025, au = .025,
         a = .9, n = 102),
    list(d = -1, cw = c(.5, .5, -1), w = .6, al = .005, au = .095, a = .8,
         n = 73),
    list(d = 2.5, cw = c(1, -1), w = .6, al = .001, au = .049, a = .8,
         n = 231))
  for (k in cases) {
    n0 <- ss_aipe_sc(true_std_contrast = k$d, c_weights = k$cw, width = k$w,
                     conf_level = NULL, alpha_lower = k$al,
                     alpha_upper = k$au)$value
    m <- n0
    while (.sc_exact_width_probability(m, k$d, k$cw, k$al, k$au, k$w) < k$a)
      m <- m + 1
    expect_identical(m, k$n)
    expect_identical(ss_aipe_sc(true_std_contrast = k$d, c_weights = k$cw,
                                width = k$w, conf_level = NULL,
                                alpha_lower = k$al, alpha_upper = k$au,
                                assurance = k$a)$value[1], k$n)
  }
})

test_that("the plan of ss_aipe_sc() attains its probability in raw two-group data", {
  skip_on_cran()  # Monte Carlo
  # 4000 data sets (all usable) of two groups in test-score units at the
  # published plan of DMAR 1.0.x (49) and the plan (50): the shares against
  # the exact .949878 and .994636.
  share <- function(n, R = 4000) {
    w <- vapply(seq_len(R), function(r) {
      y1 <- stats::rnorm(n, 100, 15)
      y2 <- stats::rnorm(n, 100, 15)
      s <- sqrt((stats::var(y1) + stats::var(y2)) / 2)
      ci <- ci_sc(means = c(mean(y1), mean(y2)), s_anova = s, c_weights = c(1, -1),
                  n_per_group = n)
      ci$value[3] - ci$value[1]
    }, numeric(1))
    expect_true(all(is.finite(w)))
    p <- mean(w <= .8)
    c(p = p, se = sqrt(p * (1 - p) / R))
  }
  set.seed(4950)
  at_49 <- share(49)
  at_50 <- share(50)
  expect_lt(abs(at_49[["p"]] - .949878) / at_49[["se"]], 4)
  expect_lt(abs(at_50[["p"]] - .994636) / at_50[["se"]], 4)
  expect_gt(at_50[["p"]], .95)
})

test_that("ss_aipe_sc_sensitivity() passes assurance_method to the planner", {
  sens <- function(...) {
    res <- ss_aipe_sc_sensitivity(true_std_contrast = 0, supposed_std_contrast = 0,
                                  c_weights = c(1, -1), width = .8, assurance = .95,
                                  G = 2, print_iter = FALSE, seed = 1, ...)
    res$value[res$term == "n_per_group"]
  }
  expect_identical(formals(ss_aipe_sc_sensitivity)$assurance_method, "analytic")
  expect_identical(sens(), 50)
  expect_identical(sens(assurance_method = "analytic"), 50)
  expect_identical(sens(), ss_aipe_sc(true_std_contrast = 0, c_weights = c(1, -1),
                                      width = .8, assurance = .95)$value[1])
  for (bad in list("simulation", NA)) {
    err <- tryCatch(sens(assurance_method = bad), error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "assurance_method")
  }
  err <- tryCatch(sens(assurance_method = "simulation"), error = function(e) e)
  expect_match(conditionMessage(err), "is not offered by ss_aipe_sc_sensitivity()",
               fixed = TRUE)
})

# The width of the ci_sc() interval at the planning value with n per group:
# the criterion of the plan without assurance.
.sc_planned_width <- function(n, d, cw, al = .025, au = .025) {
  ci <- ci_sc(means = cw * d / sum(cw^2), s_anova = 1, c_weights = cw,
              n_per_group = rep(n, length(cw)), conf_level = NULL,
              alpha_lower = al, alpha_upper = au)
  ci$value[ci$term == "ci_upper"] - ci$value[ci$term == "ci_lower"]
}

test_that("ss_aipe_sc() plans the smallest n of at least 4 whose interval meets the width", {
  skip_on_cran()
  # The search doubles a step and bisects; the plan must be what raising n
  # by one from the documented start, max(4, ceiling(n_0) - 5) with n_0 =
  # sum(c^2) (z_l + z_u)^2 / width^2, gives, checked here through ci_sc().
  cases <- list(list(d = 2, cw = c(1, -1), w = .5, al = .025, au = .025),
                list(d = -1, cw = c(.5, .5, -1), w = .4, al = .01, au = .04),
                list(d = .6, cw = c(.5, .5, -.5, -.5, 0), w = .4, al = .025, au = .025))
  for (k in cases) {
    plan <- ss_aipe_sc(true_std_contrast = k$d, c_weights = k$cw, width = k$w,
                       conf_level = NULL, alpha_lower = k$al,
                       alpha_upper = k$au)$value
    n <- max(4, ceiling(sum(k$cw^2) * (stats::qnorm(1 - k$al) +
                                         stats::qnorm(1 - k$au))^2 / k$w^2) - 5)
    while (.sc_planned_width(n, k$d, k$cw, k$al, k$au) > k$w) n <- n + 1
    expect_identical(plan, n)
    expect_gt(.sc_planned_width(plan - 1, k$d, k$cw, k$al, k$au), k$w)
  }
  # A narrow target, about 8 million per group (raising n by one from the
  # start took minutes): the interval meets the width at the plan and not
  # at one fewer.
  plan <- ss_aipe_sc(true_std_contrast = .5, c_weights = c(1, -1), width = .002)$value
  expect_identical(plan, 7923009)
  expect_lte(.sc_planned_width(plan, .5, c(1, -1)), .002)
  expect_gt(.sc_planned_width(plan - 1, .5, c(1, -1)), .002)
  # The floor of 4: at a contrast of 0 the interval at t = 0 is
  # +/- qnorm(.975) sqrt(2 / n) exactly, 3.92 wide at 2 per group, yet the
  # plan for a width of 5 is 4.
  expect_identical(ss_aipe_sc(true_std_contrast = 0, c_weights = c(1, -1), width = 5)$value, 4)
  expect_equal(.sc_planned_width(2, 0, c(1, -1)), 2 * stats::qnorm(.975), tolerance = 1e-9)
})

test_that("ss_aipe_sc() refuses a width that needs more than 2^53 per group, on every path", {
  skip_on_cran()
  for (args in list(list(), list(assurance = .9),
                    list(assurance = .9, assurance_method = "analytic"))) {
    err <- tryCatch(do.call(ss_aipe_sc, c(list(true_std_contrast = .5,
                                               c_weights = c(1, -1),
                                               width = 1e-9), args)),
                    error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "width")
    expect_match(conditionMessage(err), "more than 2^53 per group", fixed = TRUE)
  }
  # Large plans are still returned: DMAR 1.0.0 returned 122926683,
  # 30731670565548, and 3073167056555294 for these calls at once. Past
  # about 1e14 per group the last digits rest on the precision of doubles
  # (see Details of ?ss_aipe_sc), so they are compared to 1e-12.
  expect_equal(ss_aipe_sc(true_std_contrast = 0, c_weights = c(1, -1),
                          width = 1e-6)$value, 30731670565548, tolerance = 1e-12)
  expect_equal(ss_aipe_sc(true_std_contrast = 0, c_weights = c(1, -1),
                          width = 1e-7)$value, 3073167056555294, tolerance = 1e-12)
  expect_identical(ss_aipe_sc(true_std_contrast = 0, c_weights = c(1, -1),
                              width = .0005)$value, 122926683)
  expect_identical(ss_aipe_sc(true_std_contrast = .5, c_weights = c(1, -1),
                              width = .0006)$value, 88033432)
})

test_that("ss_aipe_sc() refuses from the search a plan whose start is below 2^53 and whose plan is above it", {
  skip_on_cran()
  # The start, ceiling(n_0) - 5, is 9,006,298,534,815,513, below 2^53, and
  # the plan at a contrast of .1 is above it, so the refusal comes from
  # the search (.sc_smallest_n() at n_max). It stopped with
  # dmar_computation_failed before the noncentral t was computed without
  # rounding S (see R/nc_t_internals.R).
  orig <- DMAR:::.sc_smallest_n
  reached <- FALSE
  local_mocked_bindings(
    .sc_smallest_n = function(n, meets, n_max, unreachable)
      orig(n, meets, n_max, function() {
        reached <<- TRUE
        unreachable()
      }),
    .package = "DMAR")
  width <- sqrt(2 * (2 * stats::qnorm(.975))^2 / (.9999 * 2^53))
  err <- tryCatch(ss_aipe_sc(true_std_contrast = .1, c_weights = c(1, -1),
                             width = width), error = function(e) e)
  expect_true(reached)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "width")
  expect_identical(err$n_max, 2^53)
  expect_match(conditionMessage(err),
               "more than 2^53 per group (about 9.0e15, the largest number up to which R holds every whole number exactly)",
               fixed = TRUE)
})

test_that("ss_aipe_sc() returns a plan above 3e13 per group at a nonzero contrast", {
  skip_on_cran()
  # It stopped with dmar_computation_failed ("roundoff error was detected")
  # from about 3e13 per group (6.1e13 degrees of freedom) at a contrast of
  # 2, and the same at .1 near 2^53. The smallest n is 46,097,505,848,330
  # by a 512-bit computation of the interval (work/round6/evidence/m4_nct,
  # f6_oracle.R); the plan rests on the difference of two limits near
  # 9.6e6, which doubles hold to 1.9e-9, and is 1,850 below it (see
  # Details of ?ss_aipe_sc).
  plan <- ss_aipe_sc(true_std_contrast = 2, c_weights = c(1, -1), width = 1e-6)$value
  expect_equal(plan, 46097505848330, tolerance = 1e-9)
  plan <- ss_aipe_sc(true_std_contrast = .1, c_weights = c(1, -1), width = 5.9e-8)$value
  expect_lte(plan, 2^53)
  expect_equal(plan, 2 * (2 * stats::qnorm(.975))^2 * (1 + .1^2 / 8) / 5.9e-8^2,
               tolerance = 1e-8)
})

test_that("at a contrast of 0 the plan and the exact probability take the width at t = 0 from the normal quantiles", {
  # At t = 0 the limits are -/+ qnorm(.975) at any df. ci_nc_t() found them
  # only to its tolerance, which past about 1e14 per group decided the
  # plan: at width 6e-7 the smallest n is 85,365,751,570,981 (n* is
  # 85,365,751,570,980.55 by Rmpfr, for the width and the tail areas
  # (1 - .95) / 2 as R holds them) and the plan was one below it.
  expect_identical(ss_aipe_sc(true_std_contrast = 0, c_weights = c(1, -1),
                              width = 6e-7)$value, 85365751570981)
  # The example of ?ss_aipe_sc: at width 1e-6 n* is
  # 30,731,670,565,552.9985 by Rmpfr, so the plan is the smallest n
  # (work/round6/evidence/m4_nct/fix2, contrast0_help_examples.log).
  expect_identical(ss_aipe_sc(true_std_contrast = 0, c_weights = c(1, -1),
                              width = 1e-6)$value, 30731670565553)
  # At 3,073,167,056,555,294 per group the interval at t = 0 is
  # 9.5e-16 of the width too wide, so the exact probability that it is no
  # wider than 1e-7 is 0; it came out 1.
  a <- (1 - .95) / 2
  expect_identical(DMAR:::.sc_width_probability(3073167056555294, 0, c(1, -1),
                                                1e-7, a, a), 0)
})
