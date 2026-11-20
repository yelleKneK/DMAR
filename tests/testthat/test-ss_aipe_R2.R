## ss_aipe_R2() -- sample size planning for the squared multiple correlation.

test_that("ss_aipe_R2() returns a tidy data.frame with a positive necessary_N", {
  skip_on_cran()  # Shape and monotonicity grids run locally; numerical anchors stay on CRAN.
  res <- suppressWarnings(ss_aipe_R2(true_R2 = 0.30, width = 0.20, p = 5))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("term", "value"))
  expect_true("necessary_N" %in% res$term)
  expect_gt(res$value[res$term == "necessary_N"], 0)
})

test_that("ss_aipe_R2() requires more N for tighter target widths", {
  skip_on_cran()  # Shape and monotonicity grids run locally; numerical anchors stay on CRAN.
  wide  <- suppressWarnings(ss_aipe_R2(0.30, width = 0.30, p = 5)$value[1])
  tight <- suppressWarnings(ss_aipe_R2(0.30, width = 0.10, p = 5)$value[1])
  expect_gt(tight, wide)
})

test_that("ss_aipe_R2() requires more N for more predictors at the same R^2", {
  skip_on_cran()  # Extended shape and boundary checks run in the full local suite.
  few  <- suppressWarnings(ss_aipe_R2(0.30, width = 0.20, p = 2)$value[1])
  many <- suppressWarnings(ss_aipe_R2(0.30, width = 0.20, p = 10)$value[1])
  expect_gt(many, few)
})

test_that("ss_aipe_R2() accepts assurance > 0.5 and bumps N", {
  skip_on_cran()  # Extended shape and boundary checks run in the full local suite.
  no_assurance <- suppressWarnings(ss_aipe_R2(0.30, width = 0.20, p = 5)$value[1])
  with_85      <- suppressWarnings(ss_aipe_R2(0.30, width = 0.20, p = 5, assurance = 0.85)$value[1])
  expect_gte(with_85, no_assurance)
})

test_that("ss_aipe_R2() with assurance returns where the search stopped with \"invalid 'xmin' value\"", {
  # A small population R2: the one-sided lower limit the procedure uses is
  # 0, no grid value has an expected R2 that low, and optimize() was handed
  # a missing end point. Eight conditions that stopped through DMAR 1.0.x,
  # with the sample sizes now returned.
  #
  # What these values check, and what they do not. They show that the
  # search returns where it stopped. They are not a check of the value the
  # repair substitutes for the missing end point: the interval handed to
  # optimize() matters only when the population value with the widest
  # interval lies strictly inside it, and in these conditions it does not,
  # so the same sample sizes are returned when the missing lower end is
  # set to the other end of the grid, to .5, or to .999, and when the upper
  # end is always the nearest grid value or always .999 (two reviews of
  # round 2 tried seven such substitutions; none changed a result). A later
  # reader should not take a pinned value here as evidence about the
  # substituted value itself.
  N <- function(R2, p, w, a) suppressWarnings(
    ss_aipe_R2(true_R2 = R2, p = p, width = w, assurance = a))$value[1]
  expect_equal(N(.03, 1, .2, .90), 113)
  expect_equal(N(.05, 4, .2, .90), 131)
  expect_equal(N(.03, 1, .2, .99), 148)
  expect_equal(N(.03, 4, .2, .90), 111)
})

test_that("the repaired search of ss_aipe_R2() returns the remaining conditions that stopped", {
  skip_on_cran()
  N <- function(R2, p, w, a) suppressWarnings(
    ss_aipe_R2(true_R2 = R2, p = p, width = w, assurance = a))$value[1]
  expect_equal(N(.10, 10, .2, .90), 173)
  expect_equal(N(.05, 4, .2, .99), 169)
  expect_equal(N(.05, 1, .2, .99), 166)
  expect_equal(N(.05, 10, .2, .90), 107)
  expect_equal(N(.03, 4, .2, .99), 152)
})

test_that("ss_aipe_R2() assurance plans that returned before the repair are unchanged", {
  skip_on_cran()
  # The 22 conditions of a random draw of 30 that returned in DMAR 1.0.x
  # and at the earlier 1.1.0 builds, with the sample sizes they returned.
  pinned <- rbind(
    c(.20, 10, .1, .99, 880), c(.70, 4, .2, .90, 139), c(.03, 1, .1, .90, 327),
    c(.70, 4, .2, .99, 160), c(.50, 4, .1, .90, 821), c(.10, 4, .1, .99, 680),
    c(.10, 1, .1, .99, 679), c(.30, 4, .2, .99, 230), c(.30, 10, .1, .90, 921),
    c(.50, 1, .1, .99, 846), c(.70, 1, .1, .90, 457), c(.70, 1, .2, .99, 155),
    c(.30, 4, .2, .90, 229), c(.90, 1, .2, .99, 52), c(.30, 10, .2, .99, 239),
    c(.20, 10, .2, .90, 226), c(.10, 4, .2, .90, 172), c(.50, 10, .1, .99, 862),
    c(.50, 4, .2, .90, 218), c(.20, 10, .1, .90, 846), c(.20, 1, .2, .90, 217),
    c(.20, 1, .2, .99, 225))
  for (i in seq_len(nrow(pinned))) {
    a <- pinned[i, ]
    n <- suppressWarnings(ss_aipe_R2(true_R2 = a[1], p = a[2],
                                     width = a[3], assurance = a[4]))$value[1]
    expect_equal(n, a[5], info = paste(a[1:4], collapse = " "))
  }
})

test_that("an abbreviated which_width plans for the width it stands for", {
  skip_on_cran()
  # DMAR 1.0.x checked which_width with char.expand() and kept the value as
  # given, so "F", "L", "Up", and "" matched none of the three plans and
  # the function returned NULL.
  plan <- function(w) suppressWarnings(ss_aipe_R2(
    true_R2 = .5, width = .15, which_width = w, p = 5,
    random_predictors = FALSE))$value
  expect_identical(plan("Full"), 264)
  expect_identical(plan("Lower"), 117)
  expect_identical(plan("Upper"), 13)
  for (w in c("F", "Fu")) expect_identical(plan(w), 264)
  for (w in c("L", "Lo", "Low")) expect_identical(plan(w), 117)
  for (w in c("U", "Up", "Upp")) expect_identical(plan(w), 13)
  for (w in list("", NA_character_, "full", "x", 1, c("Full", "Lower"))) {
    e <- tryCatch(plan(w), error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, "which_width")
  }
})

test_that("assurance with a half-width is refused, naming both arguments", {
  skip_on_cran()
  # The assurance step is for the full width. For the distance to one limit
  # it does not deliver the probability asked for, so the combination is
  # refused before the search starts, with either kind of predictors and
  # for an abbreviated which_width.
  for (w in c("Lower", "L", "Upper", "Up")) for (random in c(TRUE, FALSE)) {
    e <- tryCatch(
      ss_aipe_R2(true_R2 = .30, width = .10, p = 3, which_width = w,
                 assurance = .80, random_predictors = random),
      error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, c("assurance", "which_width"))
    expect_match(conditionMessage(e), "can be used only with which_width = \"Full\"")
    expect_match(conditionMessage(e),
                 if (startsWith(w, "L")) "which_width = \"Lower\"" else "which_width = \"Upper\"")
    # The message does not send the reader to the half-width plan.
    expect_no_match(conditionMessage(e), "without 'assurance'")
  }
  # Each half-width is planned without assurance, and the full width with it.
  plan <- function(...) suppressWarnings(
    ss_aipe_R2(true_R2 = .30, width = .10, p = 3, ...))$value
  expect_identical(plan(which_width = "Lower"), 253)
  expect_identical(plan(which_width = "Upper"), 198)
  expect_identical(plan(which_width = "Full", assurance = .80), 910)
})

test_that("a half-width plan is the smallest N from which the distance stays within the width", {
  skip_on_cran()
  # The case on the help page. The distance from the upper limit to the
  # expected sample R2 is .043 at N = 23, where the search begins, is
  # larger than .05 at every N from 27 to 201, and the plan is 202.
  plan <- suppressWarnings(
    ss_aipe_R2(true_R2 = .70, width = .05, p = 10, which_width = "Upper"))
  expect_identical(plan$value, 202)
  distance <- function(N) {
    e <- expected_R2(true_R2 = .70, N = N, p = 10)$value[1]
    ci <- suppressWarnings(ci_R2(R2 = e, N = N, p = 10))
    ci$value[ci$term == "ci_upper"] - e
  }
  expect_identical(round(distance(23), 3), .043)
  expect_lte(distance(26), .05)
  for (N in c(27:60, 100, 150, 201)) expect_gt(distance(N), .05)
  expect_gt(distance(56), .06)
  for (N in c(202:230, 300, 1000, 5000)) expect_lte(distance(N), .05)

  # The same distance from the distribution of the sample R2 with random
  # predictors itself, in base R alone: given the predictors the F statistic
  # is noncentral F with noncentrality rho2 / (1 - rho2) times a chi square
  # on N - 1 degrees of freedom. (ci_R2() uses Lee's approximation to it.)
  expected <- function(rho2, N, p) {
    k <- 0:3000
    1 - (N - p - 1) / (N - 1) * (1 - rho2) *
      sum(exp(lgamma(1 + k) + lgamma((N + 1) / 2) - lgamma((N + 1) / 2 + k) +
                k * log(rho2)))
  }
  cdf_random <- function(R2, rho2, N, p) {
    F_obs <- (R2 / p) / ((1 - R2) / (N - p - 1))
    stats::integrate(
      function(u) stats::pf(F_obs, p, N - p - 1, ncp = rho2 / (1 - rho2) * u) *
        stats::dchisq(u, N - 1),
      stats::qchisq(1e-12, N - 1), stats::qchisq(1 - 1e-12, N - 1),
      rel.tol = 1e-10, subdivisions = 500L)$value
  }
  exact_distance <- function(rho2, N, p) {
    e <- expected(rho2, N, p)
    stats::uniroot(function(r) cdf_random(e, r, N, p) - .025,
                   c(e, 1 - 1e-9), tol = 1e-10)$root - e
  }
  expect_gt(exact_distance(.70, 201, 10), .05)
  expect_lte(exact_distance(.70, 202, 10), .05)
  expect_gt(exact_distance(.70, 35, 5), .10)
  expect_lte(exact_distance(.70, 36, 5), .10)
  # Two more plans with random predictors that lie above a first sample
  # size that meets the target. (For the second the exact distribution
  # gives 34: at N = 33 its distance is .1503 and that of ci_R2() is under
  # .15.)
  plan_upper <- function(R2, p, w, ...) suppressWarnings(ss_aipe_R2(
    true_R2 = R2, width = w, p = p, which_width = "Upper", ...))$value
  expect_identical(plan_upper(.70, 5, .10), 36)
  expect_identical(plan_upper(.50, 5, .15), 33)
})

test_that("half-width plans with fixed predictors are those of the noncentral F in base R", {
  skip_on_cran()
  # Base R alone: the expected sample R2 from the series of
  # 2F1(1, 1; (N + 1) / 2; rho2), and the 95 percent limits from the
  # noncentral F distribution by uniroot(). The plan is one more than the
  # largest N from 2p + 3 to 500 at which the distance from the expected
  # sample R2 to the limit is larger than the width, and 2p + 3 when there
  # is none; the distance is also within the width at four larger N.
  expected <- function(rho2, N, p) {
    k <- 0:3000
    1 - (N - p - 1) / (N - 1) * (1 - rho2) *
      sum(exp(lgamma(1 + k) + lgamma((N + 1) / 2) - lgamma((N + 1) / 2 + k) +
                k * log(rho2)))
  }
  limit <- function(R2, N, p, prob) {
    F_obs <- (R2 / p) / ((1 - R2) / (N - p - 1))
    f <- function(ncp) stats::pf(F_obs, p, N - p - 1, ncp = ncp) - prob
    if (f(0) <= 0) return(0)
    ncp <- stats::uniroot(f, c(0, 20 * N), tol = 1e-10)$root
    ncp / (ncp + N)
  }
  distance <- function(side, rho2, N, p) {
    e <- expected(rho2, N, p)
    if (side == "Upper") limit(e, N, p, .025) - e else e - limit(e, N, p, .975)
  }
  # which_width, true_R2, p, width, the plan, and the first N from 2p + 3
  # at which the distance is within the width. In the first six that first
  # N is 2p + 3 and larger sample sizes miss the target, so the plan is
  # above it; in the others no larger N misses and it is the plan.
  cases <- data.frame(
    side  = c(rep("Upper", 9), rep("Lower", 3)),
    rho2  = c(.30, .50, .70, .10, .50, .30, .10, .50, .70, .30, .50, .70),
    p     = c(  5,   5,   5,  10,  10,  10,   3,   5,  10,   5,  10,   3),
    width = c(.10, .10, .05, .05, .05, .05, .10, .15, .05, .10, .15, .05),
    plan  = c(107,  44, 123, 216, 279, 449, 117,  13,  23, 265, 148, 351),
    first = c( 13,  13,  13,  23,  23,  23, 117,  13,  23, 265, 148, 351))
  for (i in seq_len(nrow(cases))) {
    g <- cases[i, ]
    info <- paste(g$side, g$rho2, g$p, g$width)
    Ns <- as.numeric((2 * g$p + 3):500)
    d <- vapply(Ns, function(N) distance(g$side, g$rho2, N, g$p), numeric(1))
    far <- vapply(c(1000, 2500, 10000, 50000),
                  function(N) distance(g$side, g$rho2, N, g$p), numeric(1))
    expect_true(all(far <= g$width), info = info)
    expect_identical(Ns[which(d <= g$width)[1L]], g$first, info = info)
    expect_identical(
      if (any(d > g$width)) max(Ns[d > g$width]) + 1 else 2 * g$p + 3,
      g$plan, info = info)
    expect_identical(suppressWarnings(ss_aipe_R2(
      true_R2 = g$rho2, width = g$width, p = g$p, which_width = g$side,
      random_predictors = FALSE))$value, g$plan, info = info)
  }
})

test_that("a half-width plan is the first N that meets the target when no larger N misses it", {
  skip_on_cran()
  # Plans at which the distance is within the width from the first N of
  # an upward search from 2p + 3 that meets the target: which_width, random
  # predictors (1) or fixed (0), conf_level, true_R2, p, width, and the
  # plan.
  pinned <- list(
    list("Lower", 1, .95, .10,  3, .05,  486), list("Lower", 1, .95, .30,  3, .10, 253),
    list("Lower", 1, .95, .50,  5, .15,  129), list("Lower", 1, .95, .70, 10, .10, 180),
    list("Lower", 1, .95, .30, 10, .05, 1150), list("Lower", 0, .95, .10,  3, .05, 479),
    list("Lower", 1, .99, .50,  5, .10,  410), list("Lower", 0, .90, .90,  2, .05,  54),
    list("Lower", 1, .95, .05, 20, .20,  128), list("Upper", 1, .95, .10,  3, .05, 526),
    list("Upper", 1, .95, .30,  5, .10,  170), list("Upper", 1, .95, .50,  3, .15,  54),
    list("Upper", 1, .95, .70,  3, .05,  306), list("Upper", 1, .95, .30, 10, .15,  23),
    list("Upper", 1, .95, .10, 10, .05,  289), list("Upper", 0, .95, .30,  3, .05, 670),
    list("Upper", 0, .95, .50,  3, .15,    9), list("Upper", 1, .90, .30,  8, .10,  19),
    list("Upper", 0, .99, .50,  1, .10,  188), list("Upper", 0, .90, .90, 15, .03,  33))
  for (a in pinned) {
    n <- suppressWarnings(ss_aipe_R2(
      true_R2 = a[[4]], width = a[[6]], p = a[[5]], which_width = a[[1]],
      conf_level = a[[3]], random_predictors = a[[2]] == 1))$value
    expect_identical(n, a[[7]], info = paste(unlist(a[1:6]), collapse = " "))
  }
})

test_that("the distance is within the width at a half-width plan and beyond it, and not just below it", {
  skip_on_cran()
  # For each plan: the distance from the expected sample R2 to the limit,
  # from expected_R2() and ci_R2(), is no larger than the width at the
  # plan, at the 10 sample sizes after it, and at a spread of larger ones,
  # and is larger than the width at plan - 1 unless the plan is 2p + 3,
  # where the search begins.
  conditions <- expand.grid(
    side = c("Lower", "Upper"), random = c(TRUE, FALSE),
    true_R2 = c(.10, .50, .80), p = c(2, 6, 12), stringsAsFactors = FALSE)
  conditions$width <- rep_len(c(.04, .08, .12, .20), nrow(conditions))
  conditions$conf_level <- rep_len(c(.95, .90, .99), nrow(conditions))
  at_floor <- 0
  for (i in seq_len(nrow(conditions))) {
    g <- conditions[i, ]
    info <- paste(unlist(g), collapse = " ")
    distance <- function(N) {
      e <- expected_R2(true_R2 = g$true_R2, N = N, p = g$p)$value[1]
      ci <- suppressWarnings(ci_R2(R2 = e, N = N, p = g$p,
                                   conf_level = g$conf_level,
                                   random_predictors = g$random))
      if (g$side == "Lower") e - ci$value[ci$term == "ci_lower"]
      else ci$value[ci$term == "ci_upper"] - e
    }
    plan <- suppressWarnings(ss_aipe_R2(
      true_R2 = g$true_R2, width = g$width, p = g$p, which_width = g$side,
      conf_level = g$conf_level, random_predictors = g$random))$value
    expect_gte(plan, 2 * g$p + 3)
    later <- unique(c(plan + 0:10, round(plan * c(1.2, 1.5, 2, 3, 5, 10, 30))))
    expect_true(all(vapply(later, distance, numeric(1)) <= g$width + 1e-9),
                info = info)
    if (plan > 2 * g$p + 3) {
      expect_gt(distance(plan - 1), g$width)
    } else {
      at_floor <- at_floor + 1
    }
  }
  # Both kinds of plan are among the 36: at 2p + 3 and above it.
  expect_gt(at_floor, 0)
  expect_lt(at_floor, nrow(conditions))
})

test_that("the half-width search returns one more than the largest N that misses", {
  skip_on_cran()
  # .R2_half_width_search_N() on profiles whose answer is known by
  # evaluating every N: a distance A / sqrt(N) - B / N, which rises to one
  # maximum (of height h at N = N_peak) and then falls, beside a covering
  # width that falls, with and without a rounding term of the size the
  # limits carry.
  profile <- function(h, N_peak, rounding) {
    A <- 2 * h * sqrt(N_peak)
    B <- h * N_peak
    upper <- function(N) A / sqrt(N) - B / N + rounding * sin(12.9898 * N)
    list(upper = upper, distances = function(N) {
      up <- upper(N)
      low <- A / sqrt(N) + B / N
      c(up, max(up, low, up + low))
    })
  }
  brute <- function(upper, N_min, width, N_top) {
    Ns <- N_min:N_top
    miss <- upper(Ns) - width > 1e-9
    if (any(miss)) max(Ns[miss]) + 1 else N_min
  }
  set.seed(20261114)
  kinds <- c(first = 0, later = 0, above_first = 0)
  for (i in 1:200) {
    h <- runif(1, .03, .3)
    N_peak <- exp(runif(1, log(2), log(400)))
    N_min <- sample(5:40, 1)
    width <- h * runif(1, .5, 1.5)
    pr <- profile(h, N_peak, rounding = sample(c(0, 1e-6, 1e-5), 1))
    got <- .R2_half_width_search_N(pr$distances, N_min, width, tol = 1e-9)
    expect_identical(got, brute(pr$upper, N_min, width, N_top = 20000),
                     info = paste(h, N_peak, N_min, width))
    # The plan is the first sample size; or the first sample size misses
    # and the plan is later; or the first sample size meets the target and
    # larger ones miss it.
    kind <- if (got == N_min) "first"
            else if (pr$upper(N_min) - width > 1e-9) "later" else "above_first"
    kinds[kind] <- kinds[kind] + 1
  }
  expect_true(all(kinds >= 20))
  # One of the third kind: the distance is within .05 at N = 10 and to
  # N = 30, larger from 31 to 209, and within it from 210 on.
  pr <- profile(h = .0625, N_peak = 64, rounding = 0)
  expect_lte(pr$upper(10), .05)
  expect_lte(pr$upper(30), .05)
  expect_gt(pr$upper(31), .05)
  expect_gt(pr$upper(209), .05)
  expect_identical(.R2_half_width_search_N(pr$distances, 10, .05, tol = 1e-9), 210)
  # A sample size at which the distance cannot be computed misses.
  d <- function(N) { x <- if (N < 40) NA_real_ else 3 / sqrt(N); c(x, 2 * x) }
  expect_identical(.R2_half_width_search_N(d, 10, 1, tol = 1e-9), 40)
  expect_identical(.R2_half_width_search_N(d, 10, .1, tol = 1e-9), 900)
  # Within the width from the first sample size on: that size is the plan.
  expect_identical(
    .R2_half_width_search_N(function(N) c(.01, .02), 7, .1, tol = 1e-9), 7)
  # A maximum one sample size wide. With the width 2e-6 below the height of
  # the maximum only N_peak misses, the first sample size meets the width,
  # and the search has to close on the maximum to find the miss.
  peaks <- expand.grid(N_peak = 20:80, N_min = c(7, 10, 13, 19))
  peaks$misses <- peaks$plan <- NA_real_
  for (i in seq_len(nrow(peaks))) {
    pr <- profile(h = .0625, N_peak = peaks$N_peak[i], rounding = 0)
    Ns <- peaks$N_min[i]:5000
    miss <- Ns[pr$upper(Ns) - (.0625 - 2e-6) > 1e-9]
    peaks$misses[i] <- if (length(miss) == 1L) miss else NA_real_
    peaks$plan[i] <- .R2_half_width_search_N(pr$distances, peaks$N_min[i],
                                             .0625 - 2e-6, tol = 1e-9)
  }
  expect_identical(peaks$misses, as.numeric(peaks$N_peak))
  expect_identical(peaks$plan, peaks$N_peak + 1)
  # Sample sizes that miss beyond the one the bisection ends on, as the
  # rounding of the limits produces. The distance is scale / N, within .01
  # from N = 100 * scale, and is raised above .01 at the sizes in `at`.
  # With misses at 102, 108, and 112 the bisection ends at 103 and the plan
  # is 113: not 103 (no step past the end of the bisection), 109 (the
  # search ended at the first sizes within the width after a miss), or 112
  # (a size that misses).
  bumps <- function(at, scale = 1) {
    function(N) c(scale / N + .002 * (N %in% at), 4 * scale / N + .002)
  }
  search_bumps <- function(at, scale = 1) {
    .R2_half_width_search_N(bumps(at, scale), 10, .01, tol = 1e-9)
  }
  expect_identical(search_bumps(c(102, 108, 112)), 113)
  # The search ends only after max(10, N / 50) consecutive sizes within the
  # width. Nine of them, 103 to 111 and 113 to 121, are not enough at
  # N = 103; nor are nineteen, 1000 to 1018, at N = 1000.
  expect_identical(search_bumps(c(102, 112, 122)), 123)
  expect_identical(search_bumps(1019, scale = 10), 1020)
  # A sample size the search cannot bound by N_max is refused, and the
  # message says what the search did not find.
  e <- tryCatch(.R2_half_width_search_N(function(N) c(1, 1), 7, .1,
                                        tol = 1e-9, N_max = 1e4),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "width")
  expect_match(conditionMessage(e), "narrower than can be planned")
  expect_match(conditionMessage(e), "above 10,000, the largest sample size")
  # So is a search that has its bound and reaches N_max while it steps
  # through the sample sizes after the plan: with the distance 10 / N the
  # plan for .01 is 1000, and 20 sizes within the width are needed.
  ten_over_N <- function(N) c(10 / N, 10 / N)
  e <- tryCatch(.R2_half_width_search_N(ten_over_N, 10, .01, tol = 1e-9,
                                        N_max = 1005),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "width")
  expect_match(conditionMessage(e),
               "reached 1,005, the largest sample size it tries")
  expect_identical(
    .R2_half_width_search_N(ten_over_N, 10, .01, tol = 1e-9, N_max = 1100),
    1000)
  # Through ss_aipe_R2(), where the search stops at 10,000,000.
  e <- tryCatch(
    suppressWarnings(ss_aipe_R2(true_R2 = .30, width = .0001, p = 5,
                                which_width = "Upper")),
    error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "width")
  expect_match(conditionMessage(e), "narrower than can be planned")
  expect_match(conditionMessage(e), "10,000,000")
  # The refusal is about the bound and does not mean that no sample size
  # meets the width. With fixed predictors, true_R2 = .98, and p = 5 the
  # distance to the upper limit is within .0005 from N = 11,351, and the
  # call is refused: the bound is a sample size at which the full width is
  # within .0005, and before that the noncentrality parameter passes two
  # million, where pf() is no longer accurate and the lower limit goes
  # wrong (the distance to it is .00035 at N = 26,624 and .00093 at
  # N = 53,248, where it should be smaller).
  distances <- function(N) {
    e <- expected_R2(true_R2 = .98, N = N, p = 5)$value[1]
    ci <- ci_R2(R2 = e, N = N, p = 5, random_predictors = FALSE)
    c(lower = e - ci$value[ci$term == "ci_lower"],
      upper = ci$value[ci$term == "ci_upper"] - e)
  }
  expect_gt(distances(11350)[["upper"]], .0005)
  for (N in c(11351, 11400, 13312, 26624)) {
    expect_lte(distances(N)[["upper"]], .0005)
  }
  expect_gt(sum(distances(26624)), .0005)
  expect_gt(distances(53248)[["lower"]], 2 * distances(26624)[["lower"]])
  e <- tryCatch(
    ss_aipe_R2(true_R2 = .98, width = .0005, p = 5, which_width = "Upper",
               random_predictors = FALSE),
    error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "width")
  expect_match(conditionMessage(e), "narrower than can be planned")
  expect_match(conditionMessage(e), "(pf() and pbeta()) are accurate",
               fixed = TRUE)
})

test_that("a half-width plan steps past the rounding of the limits for random predictors", {
  skip_on_cran()
  # The limits for random predictors are computed to a tolerance, so where
  # the distance passes the width on its way down it can pass it more than
  # once. With true_R2 = .30, p = 15, and a 99 percent interval the
  # distance to the upper limit is within .03 at N = 3687 to 3689, above
  # .03 by less than .000001 at N = 3690, and within .03 from N = 3691 on.
  # The plan is 3691, the smallest N from which the distance stays within
  # the width, and not 3687, the first N of the run that meets it.
  distance <- function(N) {
    e <- expected_R2(true_R2 = .30, N = N, p = 15)$value[1]
    ci <- suppressWarnings(ci_R2(R2 = e, N = N, p = 15, conf_level = .99))
    ci$value[ci$term == "ci_upper"] - e
  }
  expect_gt(distance(3686), .03)
  for (N in 3687:3689) expect_lte(distance(N), .03)
  expect_gt(distance(3690), .03)
  expect_lt(distance(3690), .03 + 1e-6)
  for (N in c(3691:3770, 4000, 6000, 20000)) expect_lte(distance(N), .03)
  expect_identical(suppressWarnings(ss_aipe_R2(
    true_R2 = .30, width = .03, p = 15, conf_level = .99,
    which_width = "Upper"))$value, 3691)
})

test_that("a lower limit within 1e-6 of 0 is taken as 0 in a half-width search", {
  skip_on_cran()
  # With true_R2 = .10, p = 3, and a 99 percent interval the lower limit at
  # the expected sample R2 is 0 up to N = 98 and .0003 at N = 100. At
  # N = 99 it lies between 0 and 1e-6: the distribution function of the
  # sample R2 at the expected value is above .995 for a population value of
  # 0 and below .995 for one of 1e-6. The bisection for the limit begins at
  # 1e-6 and does not find it, and ci_R2() warns there that the lower limit
  # is above the upper limit. A half-width search that passes through
  # N = 99 takes the lower limit as 0 and does not pass the warning on.
  expected <- function(N) expected_R2(true_R2 = .10, N = N, p = 3)$value[1]
  limits <- function(N) {
    ci <- suppressWarnings(ci_R2(R2 = expected(N), N = N, p = 3,
                                 conf_level = .99))
    c(lower = ci$value[ci$term == "ci_lower"],
      upper = ci$value[ci$term == "ci_upper"])
  }
  expect_gt(.lee_random_R2_cdf(expected(99), 0, 99, 3), .995)
  expect_lt(.lee_random_R2_cdf(expected(99), 1e-6, 99, 3), .995)
  expect_identical(limits(98)[["lower"]], 0)
  expect_identical(round(limits(100)[["lower"]], 4), .0003)
  plan <- function(which_width, width) {
    seen <- new.env()
    seen$classes <- character()
    value <- withCallingHandlers(
      ss_aipe_R2(true_R2 = .10, width = width, p = 3, conf_level = .99,
                 which_width = which_width)$value,
      warning = function(w) {
        seen$classes <- c(seen$classes, class(w)[1L])
        invokeRestart("muffleWarning")
      })
    list(value = value, classes = seen$classes)
  }
  # "Upper", width .175: the plan is 94, and the search evaluates the ten
  # sample sizes after it.
  expect_gt(limits(93)[["upper"]] - expected(93), .175)
  expect_lte(limits(94)[["upper"]] - expected(94), .175)
  upper <- plan("Upper", .175)
  expect_identical(upper$value, 94)
  expect_false("dmar_result_qualified" %in% upper$classes)
  # "Lower", width .1255: the distance at N = 99 is the expected sample R2
  # less a limit below 1e-6, which is larger than the width, and the
  # distance at N = 100 is within it, so the plan is 100.
  expect_gt(expected(99) - 1e-6, .1255)
  expect_lte(expected(100) - limits(100)[["lower"]], .1255)
  lower <- plan("Lower", .1255)
  expect_identical(lower$value, 100)
  expect_false("dmar_result_qualified" %in% lower$classes)
  # "Lower", width .127: the plan, 95, is below N = 99, where the distance
  # is within the width as well.
  expect_gt(expected(94) - limits(94)[["lower"]], .127)
  expect_lte(expected(95), .127)
  lower <- plan("Lower", .127)
  expect_identical(lower$value, 95)
  expect_false("dmar_result_qualified" %in% lower$classes)
})

test_that("verify_ss = TRUE is refused with a half-width", {
  skip_on_cran()
  for (w in c("Lower", "Upper")) {
    e <- tryCatch(
      suppressWarnings(ss_aipe_R2(true_R2 = .30, width = .10, p = 3,
                                  which_width = w, verify_ss = TRUE)),
      error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_match(conditionMessage(e), "verify_ss can only be used with")
  }
})

test_that("a width that is left out is refused by name", {
  skip_on_cran()
  # width has no default to plan with. Left out, or given as NULL, it is
  # refused in its own name before any search starts, for each which_width
  # and with assurance.
  for (w in c("Full", "Lower", "Upper")) {
    e <- tryCatch(ss_aipe_R2(true_R2 = .30, p = 3, which_width = w),
                  error = function(e) e)
    expect_s3_class(e, "dmar_argument_invalid")
    expect_identical(e$arg, "width")
    expect_match(conditionMessage(e), "You need to specify 'width'", fixed = TRUE)
  }
  e <- tryCatch(ss_aipe_R2(true_R2 = .30, p = 3, width = NULL, assurance = .80),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "width")
})
