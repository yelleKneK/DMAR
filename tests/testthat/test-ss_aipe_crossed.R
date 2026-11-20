## ss_aipe_crossed(): closed-form anchors and independent re-derivations of
## the planner's criterion. Every expected value below is computed in this
## file from the variance of the grand mean of a units by raters array,
##   V = var_unit / n_u + var_rater / n_r + var_residual / (n_u n_r),
## and, for the t quantile, from the Satterthwaite degrees of freedom of
## (MS_u + MS_r - MS_e) / (n_u n_r) written out from the expected mean
## squares, never by calling the package's internal helpers.

# The planned half-width at (n_u, n_r), written out independently. With
# an assurance the chi square factor is floored at 1, so that an assurance
# plan is never smaller than the plan without one; floor = FALSE gives the
# unfloored factor, to show what the floor prevents.
crossed_half_width <- function(n_u, n_r, vu, vr, ve, conf_level = 0.95,
                               quantile = "t", assurance = NULL,
                               floor = TRUE) {
  V <- vu / n_u + vr / n_r + ve / (n_u * n_r)
  p <- 1 - (1 - conf_level) / 2
  if (quantile == "normal") return(qnorm(p) * sqrt(V))
  emu <- ve + n_r * vu
  emr <- ve + n_u * vr
  eme <- ve
  nu <- (emu + emr - eme)^2 /
    (emu^2 / (n_u - 1) + emr^2 / (n_r - 1) + eme^2 / ((n_u - 1) * (n_r - 1)))
  k <- if (is.null(assurance)) 1 else qchisq(assurance, nu) / nu
  if (floor) k <- pmax(1, k)
  qt(p, nu) * sqrt(V * k)
}

# The smallest n_u meeting the half-width h with n_r raters, by a plain
# linear scan (no bisection, no shape assumption); NA if none up to max_u.
scan_units <- function(n_r, h, max_u, ...) {
  hw <- crossed_half_width(2:max_u, n_r, ...)
  hit <- which(hw <= h)
  if (length(hit)) hit[1] + 1 else NA_real_
}

val <- function(x, term) x$value[x$term == term]

test_that("with the normal quantile and a fixed panel, the number of units is the closed form", {
  # z^2 V <= h^2 solves to n_u >= (vu + ve / n_r) / ((h / z)^2 - vr / n_r).
  vu <- 0.50; vr <- 0.10; ve <- 0.40; w <- 0.40
  z <- qnorm(0.975); h <- w / 2
  for (n_r in c(12, 20, 30, 75)) {
    res <- ss_aipe_crossed(vu, vr, ve, width = w, n_raters = n_r,
                           quantile = "normal")
    closed <- ceiling((vu + ve / n_r) / ((h / z)^2 - vr / n_r))
    n_u <- val(res, "necessary_n_units")
    expect_equal(n_u, max(2, closed), info = n_r)
    expect_equal(val(res, "specified_n_raters"), n_r)
    # The design meets the width and the next smaller one does not.
    expect_lte(crossed_half_width(n_u, n_r, vu, vr, ve, quantile = "normal"), h)
    expect_gt(crossed_half_width(n_u - 1, n_r, vu, vr, ve, quantile = "normal"), h)
  }
})

# On CRAN the file stops here, after the anchor above; the blocks below run
# in the full local suite (see cran_files in tests/testthat.R).
skip_on_cran()

test_that("with judgments as the only cost, the plan is the cheapest whole-number design near the Marcoulides-Goldstein ratio", {
  vu <- 0.50; vr <- 0.10; ve <- 0.40; w <- 0.40
  z <- qnorm(0.975); v <- (w / 2 / z)^2
  res <- ss_aipe_crossed(vu, vr, ve, width = w, quantile = "normal")
  n_u <- val(res, "necessary_n_units"); n_r <- val(res, "necessary_n_raters")
  expect_equal(val(res, "total_judgments"), n_u * n_r)
  expect_equal(val(res, "cost"), n_u * n_r)

  # Brute force over every number of raters with the closed-form units.
  R <- 2:200
  room <- v - vr / R
  U <- ifelse(room > 0, pmax(2, ceiling((vu + ve / R) / room)), NA)
  best <- which.min(U * R)
  expect_equal(c(n_u, n_r), c(U[best], R[best]))

  # The continuous optimum: for a fixed number of judgments K the variance
  # vu / n_u + vr / n_r + ve / K is least at n_u / n_r = vu / vr
  # (Marcoulides & Goldstein, 1990, Equations 5 and 6), where it equals
  # 2 sqrt(vu vr / K) + ve / K. Setting that to v and solving the quadratic
  # in x = 1 / sqrt(K) gives the fewest judgments any design could use.
  x <- (-2 * sqrt(vu * vr) + sqrt(4 * vu * vr + 4 * ve * v)) / (2 * ve)
  K <- 1 / x^2
  expect_gte(n_u * n_r, K)
  expect_lt(n_u * n_r, 1.02 * K)
  expect_equal(n_u / n_r, vu / vr, tolerance = 0.2)
})

test_that("costs on units and raters move the plan to the cheapest priced design", {
  vu <- 0.50; vr <- 0.10; ve <- 0.40; w <- 0.40
  z <- qnorm(0.975); v <- (w / 2 / z)^2
  for (costs in list(c(0, 40, 1), c(5, 0, 1), c(3, 25, 0.5), c(0, 1, 0))) {
    res <- ss_aipe_crossed(vu, vr, ve, width = w, quantile = "normal",
                           unit_cost = costs[1], rater_cost = costs[2],
                           judgment_cost = costs[3])
    R <- 2:500
    room <- v - vr / R
    U <- ifelse(room > 0, pmax(2, ceiling((vu + ve / R) / room)), NA)
    cost <- costs[1] * U + costs[2] * R + costs[3] * U * R
    best <- which.min(cost)
    expect_equal(c(val(res, "necessary_n_units"), val(res, "necessary_n_raters")),
                 c(U[best], R[best]), info = paste(costs, collapse = "/"))
    expect_equal(val(res, "cost"), cost[best])
  }
})

test_that("with no rater variance the cheapest design uses two raters", {
  # V = vu / n_u + ve / (n_u n_r): for a fixed number of judgments it falls
  # as n_u grows, so the fewest raters the variance estimate allows is best.
  # This pins the normal plan's design, not its performance: the analysis
  # still estimates the rater component, which with two raters is kept in
  # about a third of samples (MS_r / MS_e is F on 1 and n_u - 1 df), and
  # the Satterthwaite df then fall to about 2, so this plan's intervals
  # are often far wider than the target. test-ss_aipe_crossed_sensitivity.R
  # shows that, and that the t plan (three raters here) does not suffer it.
  vu <- 0.30; ve <- 0.60; w <- 0.25
  z <- qnorm(0.975); v <- (w / 2 / z)^2
  res <- ss_aipe_crossed(vu, 0, ve, width = w, quantile = "normal")
  expect_equal(val(res, "necessary_n_raters"), 2)
  expect_equal(val(res, "necessary_n_units"), ceiling((vu + ve / 2) / v))
  expect_equal(val(res, "rater_floor"), 0)
})

test_that("with the t quantile the design meets the width and the next smaller does not", {
  skip_on_cran()
  vu <- 0.50; vr <- 0.10; ve <- 0.40; w <- 0.40; h <- w / 2
  res <- ss_aipe_crossed(vu, vr, ve, width = w)
  n_u <- val(res, "necessary_n_units"); n_r <- val(res, "necessary_n_raters")
  expect_lte(crossed_half_width(n_u, n_r, vu, vr, ve), h)
  expect_gt(crossed_half_width(n_u - 1, n_r, vu, vr, ve), h)
  # The t plan costs at least as much as the normal one.
  nrm <- ss_aipe_crossed(vu, vr, ve, width = w, quantile = "normal")
  expect_gte(val(res, "cost"), val(nrm, "cost"))
  # It is the cheapest over the raters, with each number of raters given
  # its smallest number of units by a linear scan.
  R <- 2:60
  U <- vapply(R, scan_units, numeric(1), h = h, max_u = 2000,
              vu = vu, vr = vr, ve = ve)
  best <- which.min(U * R)
  expect_equal(c(n_u, n_r), c(U[best], R[best]))
  # The df and expected width rows describe the planned design.
  V <- vu / n_u + vr / n_r + ve / (n_u * n_r)
  emu <- ve + n_r * vu; emr <- ve + n_u * vr
  nu <- (emu + emr - ve)^2 /
    (emu^2 / (n_u - 1) + emr^2 / (n_r - 1) + ve^2 / ((n_u - 1) * (n_r - 1)))
  expect_equal(val(res, "df"), nu, tolerance = 1e-12)
  expect_equal(val(res, "expected_width"), 2 * qt(0.975, nu) * sqrt(V),
               tolerance = 1e-12)
  expect_lte(val(res, "expected_width"), w)
})

test_that("an assurance is met at the plan and missed one unit short of it", {
  vu <- 0.50; vr <- 0.10; ve <- 0.40; w <- 0.40; h <- w / 2
  for (gam in c(0.70, 0.80, 0.95)) {
    res <- ss_aipe_crossed(vu, vr, ve, width = w, n_raters = 30,
                           assurance = gam)
    n_u <- val(res, "necessary_n_units")
    # The chi square approximation: nu V-hat / V ~ chi square on nu df, so
    # P(width <= w) = P(chi square <= nu h^2 / (t^2 V)).
    prob <- function(n) {
      V <- vu / n + vr / 30 + ve / (n * 30)
      emu <- ve + 30 * vu; emr <- ve + n * vr
      nu <- (emu + emr - ve)^2 /
        (emu^2 / (n - 1) + emr^2 / 29 + ve^2 / ((n - 1) * 29))
      pchisq(nu * h^2 / (qt(0.975, nu)^2 * V), nu)
    }
    expect_gte(prob(n_u), gam)
    expect_lt(prob(n_u - 1), gam)
    expect_equal(val(res, "assurance"), gam)
  }
  expected <- ss_aipe_crossed(vu, vr, ve, width = w, n_raters = 30)
  assured <- ss_aipe_crossed(vu, vr, ve, width = w, n_raters = 30,
                             assurance = 0.80)
  expect_gt(val(assured, "necessary_n_units"),
            val(expected, "necessary_n_units"))
  expect_lt(val(assured, "expected_width"), w)
})

test_that("an assurance plan is never smaller than the plan without one, even where the chi square quantile is below its df", {
  skip_on_cran()
  # Few raters and a rater component that dominates hold the Satterthwaite
  # df near n_r - 1 = 5. There the plan without an assurance already meets
  # the width with probability P(chi square_5 <= 5) = .58 under the chi
  # square approximation, so an assurance of .55 asks for less, and the
  # unfloored factor qchisq(.55, nu) / nu is below 1.
  vu <- 0.10; vr <- 0.20; ve <- 2.00; w <- 1.2; h <- w / 2; n_r <- 6
  expected <- ss_aipe_crossed(vu, vr, ve, width = w, n_raters = n_r)
  n_0 <- val(expected, "necessary_n_units")
  nu <- val(expected, "df")
  expect_true(nu > 4.5 && nu < 6)
  expect_gt(pchisq(nu, nu), 0.55)
  expect_lt(qchisq(0.55, nu) / nu, 1)
  # Unfloored, the factor would shrink the plan below n_0 units.
  shrunk <- scan_units(n_r, h, 2000, vu = vu, vr = vr, ve = ve,
                       assurance = 0.55, floor = FALSE)
  expect_lt(shrunk, n_0)

  assured <- ss_aipe_crossed(vu, vr, ve, width = w, n_raters = n_r,
                             assurance = 0.55)
  expect_equal(val(assured, "necessary_n_units"), n_0)
  expect_equal(val(assured, "necessary_n_units"),
               scan_units(n_r, h, 2000, vu = vu, vr = vr, ve = ve,
                          assurance = 0.55))
  expect_lte(val(assured, "expected_width"), w)
  expect_equal(val(assured, "expected_width"), val(expected, "expected_width"))
  expect_equal(val(assured, "assurance"), 0.55)

  # With the raters planned as well, the plan costs no less.
  free_0 <- ss_aipe_crossed(vu, vr, ve, width = w)
  free_a <- ss_aipe_crossed(vu, vr, ve, width = w, assurance = 0.55)
  expect_gte(val(free_a, "cost"), val(free_0, "cost"))
  expect_lte(val(free_a, "expected_width"), w)

  # Over a range of assurances the units never fall below n_0, and they
  # rise once the assurance asks for more than the plan delivers. (Six
  # raters cannot reach an assurance of .90 here at all: the panel, not
  # the units, then binds.)
  n_gam <- vapply(c(0.05, 0.3, 0.5, 0.52, 0.55, 0.58, 0.6, 0.7, 0.8),
                  function(g) val(ss_aipe_crossed(vu, vr, ve, width = w,
                                                  n_raters = n_r,
                                                  assurance = g),
                                  "necessary_n_units"),
                  numeric(1))
  expect_true(all(n_gam >= n_0))
  expect_false(is.unsorted(n_gam))
  expect_gt(n_gam[length(n_gam)], n_0)
  expect_error(ss_aipe_crossed(vu, vr, ve, width = w, n_raters = n_r,
                               assurance = 0.9),
               "More raters are needed")

  # An assurance at or below one half returns the plan without one,
  # every row the same apart from the assurance echo.
  low <- ss_aipe_crossed(vu, vr, ve, width = w, n_raters = n_r,
                         assurance = 0.3)
  expect_equal(low$value[low$term != "assurance"], expected$value)
})

test_that("past a point more units lengthen the t interval, and the plan takes the first feasible design", {
  # A large rater component and few raters: as units are added, V falls
  # toward vr / n_r while the Satterthwaite df fall toward n_r - 1, so the
  # planned half-width falls and then rises.
  vu <- 0.90; vr <- 0.35; ve <- 1.80; n_r <- 4
  hw <- crossed_half_width(2:5000, n_r, vu, vr, ve)
  turn <- which.min(hw) + 1
  expect_true(turn > 10 && turn < 5000)
  expect_gt(hw[length(hw)], min(hw))
  # A target between the minimum and the large-n limit is met only on a
  # window of units; the plan is the window's first design.
  h <- (min(hw) + hw[length(hw)]) / 2
  res <- ss_aipe_crossed(vu, vr, ve, width = 2 * h, n_raters = n_r)
  n_u <- val(res, "necessary_n_units")
  expect_equal(n_u, scan_units(n_r, h, 5000, vu = vu, vr = vr, ve = ve))
  expect_lt(n_u, turn)
  # A target below the minimum is unattainable with four raters, although
  # the normal quantile's rater floor lies below four.
  h_low <- 0.99 * min(hw)
  expect_lt(qnorm(0.975)^2 * vr / h_low^2, n_r)
  expect_error(ss_aipe_crossed(vu, vr, ve, width = 2 * h_low, n_raters = n_r),
               "No number of units")
})

test_that("the planned width falls then rises in the units, and the bisection matches a linear scan", {
  # The search assumes the t half-width, as a function of the number of
  # units, never rises and then falls again. Check that shape, and the
  # search's answer against a plain scan, over a spread of planning
  # values, panels, confidence levels, and assurances, with targets
  # chosen so that some are met on a window, some on a ray, and some not
  # at all.
  skip_on_cran()
  set.seed(1990)
  max_u <- 3000
  for (i in 1:200) {
    vu <- runif(1) * (runif(1) > 0.05)
    vr <- runif(1) * runif(1) * (runif(1) > 0.10)
    ve <- runif(1, 0.01, 2)
    n_r <- sample(2:100, 1)
    cl <- sample(c(0.80, 0.90, 0.95, 0.99), 1)
    gam <- if (runif(1) < 0.5) NULL else runif(1, 0.55, 0.99)
    hw <- crossed_half_width(2:max_u, n_r, vu, vr, ve, conf_level = cl,
                             assurance = gam)
    s <- sign(diff(hw))
    s <- s[s != 0]
    expect_false(any(diff(s) < 0), info = i)
    h <- runif(1, 0.98 * min(hw), hw[1])
    hit <- which(hw <= h)
    expected <- if (length(hit)) hit[1] + 1 else NA_real_
    got <- DMAR:::.crossed_n_units(n_r, h, max_u, function(n_u, n_r)
      DMAR:::.crossed_half_width(n_u, n_r, vu, vr, ve, cl, "t", gam))
    expect_identical(got, expected, info = i)
  }

  # The same near one half, where the floor of the assurance factor at 1
  # binds at some designs and not at others, with small panels, whose
  # df are small enough for the chi square quantile to fall below them.
  for (i in 1:100) {
    vu <- runif(1)
    vr <- runif(1) * runif(1)
    ve <- runif(1, 0.01, 2)
    n_r <- sample(2:12, 1)
    gam <- runif(1, 0.50, 0.62)
    hw <- crossed_half_width(2:max_u, n_r, vu, vr, ve, assurance = gam)
    s <- sign(diff(hw))
    s <- s[s != 0]
    expect_false(any(diff(s) < 0), info = paste("near half", i))
    expect_true(all(hw >= crossed_half_width(2:max_u, n_r, vu, vr, ve)),
                info = paste("near half", i))
    h <- runif(1, 0.98 * min(hw), hw[1])
    hit <- which(hw <= h)
    expected <- if (length(hit)) hit[1] + 1 else NA_real_
    got <- DMAR:::.crossed_n_units(n_r, h, max_u, function(n_u, n_r)
      DMAR:::.crossed_half_width(n_u, n_r, vu, vr, ve, 0.95, "t", gam))
    expect_identical(got, expected, info = paste("near half", i))
  }
})

test_that("the rater floor is reported, and a panel at or below it is refused", {
  vu <- 0.50; vr <- 0.10; ve <- 0.40; w <- 0.40
  floor_r <- qnorm(0.975)^2 * vr / (w / 2)^2
  res <- ss_aipe_crossed(vu, vr, ve, width = w)
  expect_equal(val(res, "rater_floor"), floor_r, tolerance = 1e-12)
  expect_gt(val(res, "necessary_n_raters"), floor_r)
  expect_error(ss_aipe_crossed(vu, vr, ve, width = w, n_raters = floor(floor_r),
                               quantile = "normal"),
               "rater floor")
  expect_error(ss_aipe_crossed(vu, vr, ve, width = 0.05, max_raters = 50),
               "No design with at most 50 raters")
})

test_that("an infeasible target names the limit that binds", {
  skip_on_cran()
  vu <- 0.50; vr <- 0.10; ve <- 0.40; w <- 0.40
  # Raters: the floor is 3.84 raters, and three are too few however many
  # units are added.
  msg <- tryCatch(ss_aipe_crossed(vu, vr, ve, width = w, n_raters = 3),
                  error = conditionMessage)
  expect_match(msg, "No number of units up to 100,000", fixed = TRUE)
  expect_match(msg, "(the rater floor)", fixed = TRUE)
  expect_match(msg, "More raters are needed", fixed = TRUE)
  expect_no_match(msg, "a larger 'max_units'", fixed = TRUE)
  msg <- tryCatch(ss_aipe_crossed(vu, vr, ve, width = 0.05, max_raters = 50),
                  error = conditionMessage)
  expect_match(msg, "(the rater floor)", fixed = TRUE)
  expect_match(msg, "More raters than 'max_raters' are needed", fixed = TRUE)

  # Units: with no rater variance there is no rater floor, and the review's
  # case of a tight target with a small max_units names max_units and the
  # number of units that would do.
  msg <- tryCatch(ss_aipe_crossed(0.5, 0, 0.4, width = 0.05, max_units = 50,
                                  n_raters = 4),
                  error = conditionMessage)
  need <- val(ss_aipe_crossed(0.5, 0, 0.4, width = 0.05, n_raters = 4),
              "necessary_n_units")
  expect_gt(need, 50)
  expect_match(msg, sprintf("With 4 raters the target needs %s units",
                            format(need, big.mark = ",")), fixed = TRUE)
  expect_match(msg, "a larger 'max_units' allows them", fixed = TRUE)
  expect_no_match(msg, "rater floor", fixed = TRUE)
  # The same with the raters planned: the design given is the plan the
  # default limits find.
  msg <- tryCatch(ss_aipe_crossed(vu, vr, ve, width = w, max_units = 20),
                  error = conditionMessage)
  plan <- ss_aipe_crossed(vu, vr, ve, width = w)
  expect_match(msg, sprintf("has %s units and %s raters",
                            val(plan, "necessary_n_units"),
                            val(plan, "necessary_n_raters")), fixed = TRUE)
  expect_match(msg, "a larger 'max_units' allows it", fixed = TRUE)

  # The t quantile: four raters are above the floor, but the planned
  # interval stops narrowing before it reaches the target.
  vu <- 0.90; vr <- 0.35; ve <- 1.80; n_r <- 4
  h_low <- 0.99 * min(crossed_half_width(2:5000, n_r, vu, vr, ve))
  msg <- tryCatch(ss_aipe_crossed(vu, vr, ve, width = 2 * h_low,
                                  n_raters = n_r),
                  error = conditionMessage)
  expect_match(msg, "above the rater floor", fixed = TRUE)
  expect_match(msg, "with the t quantile the planned interval stops narrowing",
               fixed = TRUE)
})

test_that("the table has the documented rows, class, and tidy summary", {
  skip_on_cran()
  res <- ss_aipe_crossed(0.50, 0.10, 0.40, width = 0.40)
  expect_s3_class(res, c("dmar_ss_aipe", "dmar_tbl", "data.frame"))
  expect_named(res, c("term", "value"))
  expect_type(res$value, "double")
  expect_identical(res$term,
                   c("necessary_n_units", "necessary_n_raters",
                     "total_judgments", "cost", "expected_width", "df",
                     "rater_floor", "true_var_unit", "true_var_rater", "true_var_residual",
                     "width", "conf_level"))
  expect_identical(attr(res, "quantile"), "t")
  expect_identical(attr(ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4,
                                        quantile = "normal"), "quantile"),
                   "normal")
  expect_equal(attr(res, "conf_level"), 0.95)
  # The registry reads the number of units as the size and the target as
  # the width; glance() keeps the raters.
  td <- generics::tidy(res)
  expect_equal(td$estimate, val(res, "necessary_n_units"))
  expect_equal(td$width, 0.40)
  gl <- generics::glance(res)
  expect_equal(gl$necessary_n_raters, val(res, "necessary_n_raters"))
  fixed <- ss_aipe_crossed(0.50, 0.10, 0.40, width = 0.40, n_raters = 25)
  expect_equal(generics::tidy(fixed)$estimate, val(fixed, "necessary_n_units"))
  expect_equal(generics::glance(fixed)$specified_n_raters, 25)
  expect_false("assurance" %in% res$term)
})

test_that("a cheapest design at the rater limit is flagged", {
  # Raters are free and units are costly, so the cheapest design uses as
  # many raters as the search allows.
  expect_warning(
    res <- ss_aipe_crossed(0.50, 0.10, 0.40, width = 0.40, unit_cost = 1,
                           judgment_cost = 0, max_raters = 40),
    "search limit 'max_raters'")
  expect_equal(val(res, "necessary_n_raters"), 40)
})

test_that("bad arguments are refused with a message naming them", {
  skip_on_cran()
  expect_error(ss_aipe_crossed(-0.1, 0.1, 0.4, width = 0.4), "'true_var_unit'")
  expect_error(ss_aipe_crossed(0.5, -0.1, 0.4, width = 0.4), "'true_var_rater'")
  expect_error(ss_aipe_crossed(0.5, 0.1, 0, width = 0.4),
               "'true_var_residual' must be a single positive number")
  expect_error(ss_aipe_crossed(0.5, 0.1, 0.4, width = 0), "'width'")
  expect_error(ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4, conf_level = 1),
               "'conf_level'")
  expect_error(ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4, assurance = 1),
               "'assurance' must be a single number strictly between 0 and 1")
  expect_error(ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4, assurance = 0),
               "'assurance' must be a single number strictly between 0 and 1")
  expect_error(ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4, assurance = "0.8"),
               "'assurance'")
  expect_error(ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4, assurance = 0.8,
                               quantile = "normal"),
               "normal quantile")
  expect_error(ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4, n_raters = 1),
               "'n_raters' must be a single whole number of at least 2")
  expect_error(ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4, n_raters = 12.5),
               "'n_raters'")
  expect_error(ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4, judgment_cost = 0),
               "must be positive")
  expect_error(ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4, rater_cost = -1),
               "'rater_cost'")
  expect_error(ss_aipe_crossed(0.5, 0.1, 0.4, width = 0.4, quantile = "z"))
})
