v <- function(tab, t) tab$value[tab$term == t]

# Oracles written here, sharing no code with the function.
# .ie_r_power(rho, n, alpha): the power of the t test of a correlation from
# n cases, from the exact density of r in its integral form. With X random
# the a test is that of r_XM on N cases, and the b test (either design) that
# of r_YM.X, a correlation from N - 1 cases.
.ie_r_power <- function(rho, n, alpha = .05) {
  dens <- function(r) vapply(r, function(ri) {
    lg <- log(n - 2) - log(pi) + (n - 1) / 2 * log1p(-rho^2) +
      (n - 4) / 2 * log1p(-ri^2)
    f <- function(w) exp(lg - (n - 1) *
                           (w + log((1 + exp(-2 * w)) / 2 - rho * ri * exp(-w))))
    stats::integrate(f, 0, Inf, rel.tol = 1e-12, abs.tol = 0)$value
  }, 0)
  tc <- stats::qt(1 - alpha / 2, n - 2)
  rc <- tc / sqrt(tc^2 + n - 2)
  stats::integrate(dens, rc, 1, rel.tol = 1e-11)$value +
    stats::integrate(dens, -1, -rc, rel.tol = 1e-11)$value
}
# .ie_power(): the power by integrate() over W = SS_M.X / (1 - a^2), chi
# square on N - 2, given which t_b is noncentral t with ncp delta_b sqrt(W)
# (pt()) and, independent of it, t_a = (Z + delta_a sqrt(SS_X)) /
# sqrt(W / (N - 2)), SS_X chi square on N - 1 or, with X fixed, N. For the
# Sobel test (X fixed), the integral over y = |t_a| > z of its density
# times P(|t_b| > y z / sqrt(y^2 - z^2) | W).
.ie_power <- function(N, a, b, cp = 0, alpha = .05, random = TRUE, sobel = FALSE) {
  da <- abs(a) / sqrt(1 - a^2)
  db <- abs(b) * sqrt(1 - a^2) / sqrt(1 - (b^2 + cp^2 + 2 * a * b * cp))
  ca <- stats::qt(1 - alpha / 2, N - 2); cb <- stats::qt(1 - alpha / 2, N - 3)
  z <- stats::qnorm(1 - alpha / 2)
  pb <- function(k, w) stats::pt(k, N - 3, db * sqrt(w), lower.tail = FALSE) + stats::pt(-k, N - 3, db * sqrt(w))
  given <- function(w) {
    s <- sqrt((N - 2) / w)
    if (sobel) {
      m <- da * sqrt(N) * s
      f <- function(y) (stats::dnorm((y - m) / s) + stats::dnorm((y + m) / s)) / s * pb(y * z / sqrt(y^2 - z^2), w)
      br <- sort(c(z, z + 2^(-12:2), m + s * c(-6, -3, 0, 3, 6, 12)))
      br <- br[br >= z & br <= m + 12 * s]
      return(sum(vapply(seq_along(br[-1]), function(j)
        stats::integrate(f, br[j], br[j + 1], rel.tol = 1e-10, abs.tol = 1e-14)$value, 0)))
    }
    pa <- function(v) stats::pnorm(da * sqrt(v) - ca / s) + stats::pnorm(-da * sqrt(v) - ca / s)
    pb(cb, w) * if (random) stats::integrate(function(u) pa(stats::qchisq(u, N - 1)), 0, 1,
                                             rel.tol = 1e-11)$value else pa(N)
  }
  stats::integrate(function(u) vapply(stats::qchisq(u, N - 2), given, 0), 0, 1, rel.tol = 1e-10)$value
}
# .ie_simulate(N, a, b, cp, R, seed, random): R raw-data studies fitted by
# least squares written here (lm()'s t values), with X ~ N(5, 2^2) when
# random and two groups coded 0/1 when fixed, M = 10 + 3 z_M, and
# Y = -4 + 0.5 z_Y: the rejection counts of both tests, and the attempted,
# usable, and failed counts.
.ie_simulate <- function(N, a, b, cp, R, seed, random = TRUE, alpha = .05) {
  set.seed(seed)
  sY <- sqrt(1 - (b^2 + cp^2 + 2 * a * b * cp))
  ca <- stats::qt(1 - alpha / 2, N - 2); cb <- stats::qt(1 - alpha / 2, N - 3)
  z <- stats::qnorm(1 - alpha / 2)
  out <- c(joint = 0, sobel = 0, usable = 0)
  done <- 0
  while (done < R) {
    n <- min(10000, R - done)
    if (random) {
      X <- matrix(stats::rnorm(n * N, 5, 2), n); zx <- (X - 5) / 2
    } else {
      X <- matrix(rep(0:1, length.out = N), n, N, byrow = TRUE)
      zx <- X - rowMeans(X); zx <- zx / sqrt(rowMeans(zx^2))
    }
    zm <- a * zx + sqrt(1 - a^2) * matrix(stats::rnorm(n * N), n)
    M <- 10 + 3 * zm
    Y <- -4 + .5 * (cp * zx + b * zm + sY * matrix(stats::rnorm(n * N), n))
    X <- X - rowMeans(X); M <- M - rowMeans(M); Y <- Y - rowMeans(Y)
    sxx <- rowSums(X^2); smm <- rowSums(M^2); sxm <- rowSums(X * M)
    sxy <- rowSums(X * Y); smy <- rowSums(M * Y); syy <- rowSums(Y^2)
    ah <- sxm / sxx
    ta <- ah / sqrt((smm - ah * sxm) / (N - 2) / sxx)
    det <- sxx * smm - sxm^2
    bm <- (sxx * smy - sxm * sxy) / det
    bx <- (smm * sxy - sxm * smy) / det
    tb <- bm / sqrt((syy - bx * sxy - bm * smy) / (N - 3) * sxx / det)
    ok <- is.finite(ta) & is.finite(tb)
    ta <- ta[ok]; tb <- tb[ok]
    out <- out + c(sum(abs(ta) > ca & abs(tb) > cb),
                   sum(abs(ta * tb) / sqrt(ta^2 + tb^2) > z), sum(ok))
    done <- done + n
  }
  c(out, attempted = R, failed = R - out[["usable"]])
}

test_that("power_a with X fixed uses the exact noncentral t where pt() is not accurate", {
  # Past a noncentrality of 37.62, checked against the oracle that
  # conditions on the normal variate (helper-nc_t_oracle.R).
  a <- .5; N <- 6000
  ncp <- a / sqrt((1 - a^2) / N)
  expect_gt(ncp, 37.62)
  cv <- stats::qt(.975, N - 2)
  oracle <- .oracle_nct(cv, N - 2, ncp, lower_tail = FALSE) + .oracle_nct(-cv, N - 2, ncp)
  r <- ss_power_indirect_effect(true_std_reg_coef_a = a, true_std_reg_coef_b = .05,
                                specified_N = N, random_predictors = FALSE)
  expect_equal(v(r, "power_a"), oracle, tolerance = 1e-12)
  calls <- all.names(body(ss_power_indirect_effect))
  expect_false("pt" %in% calls)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_power_indirect_effect() solves and evaluates coherently", {
  plan <- ss_power_indirect_effect(true_std_reg_coef_a = .39, true_std_reg_coef_b = .39, desired_power = .80)
  expect_s3_class(plan, "dmar_tbl")
  n <- v(plan, "necessary_N")
  expect_identical(n, 69)
  expect_gte(v(plan, "power_ab"), 0.80)
  # One fewer case dips below the target (minimality of the search).
  at <- ss_power_indirect_effect(true_std_reg_coef_a = .39, true_std_reg_coef_b = .39, specified_N = n - 1)
  expect_lt(v(at, "power_ab"), 0.80)
  # The power of the indirect effect is reported first, then the paths';
  # it is below their product, since the tests are dependent.
  expect_identical(plan$term[2:4], c("power_ab", "power_a", "power_b"))
  expect_lt(v(plan, "power_ab"), v(plan, "power_a") * v(plan, "power_b"))
  expect_false("actual_power" %in% plan$term)
  expect_equal(v(plan, "true_std_indirect_effect"), .39 * .39)
  expect_true(attr(plan, "random_predictors"))
  expect_identical(attr(plan, "method"), "joint_significance")
  fixed <- ss_power_indirect_effect(.39, .39, desired_power = .80, random_predictors = FALSE)
  expect_identical(v(fixed, "necessary_N"), 68)
  expect_false(attr(fixed, "random_predictors"))
})

test_that("the weak link drives the requirement, and Sobel needs more", {
  small_a <- ss_power_indirect_effect(true_std_reg_coef_a = .14, true_std_reg_coef_b = .39,
                                      desired_power = .80)
  both_med <- ss_power_indirect_effect(true_std_reg_coef_a = .39, true_std_reg_coef_b = .39,
                                       desired_power = .80)
  expect_gt(v(small_a, "necessary_N"), v(both_med, "necessary_N"))
  sobel <- ss_power_indirect_effect(true_std_reg_coef_a = .39, true_std_reg_coef_b = .39, desired_power = .80,
                                    method = "sobel")
  expect_gte(v(sobel, "necessary_N"), v(both_med, "necessary_N"))
  expect_true(is.na(v(sobel, "power_a")))
  expect_true(is.na(v(sobel, "power_b")))
})

test_that("ss_power_indirect_effect() validates its arguments", {
  expect_error(ss_power_indirect_effect(true_std_reg_coef_a = 0, true_std_reg_coef_b = .3,
                                        desired_power = .8), "nonzero")
  expect_error(ss_power_indirect_effect(true_std_reg_coef_a = 1.1, true_std_reg_coef_b = .3,
                                        desired_power = .8), "\\(-1, 1\\)")
  expect_error(ss_power_indirect_effect(true_std_reg_coef_a = .3, true_std_reg_coef_b = .3), "exactly one")
  expect_error(ss_power_indirect_effect(true_std_reg_coef_a = .3, true_std_reg_coef_b = .3, specified_N = 100,
                                        desired_power = .8), "exactly one")
  expect_error(ss_power_indirect_effect(true_std_reg_coef_a = .7, true_std_reg_coef_b = .7, true_std_reg_coef_direct = .7,
                                        desired_power = .8),
               "not admissible")
  for (bad in list(NA, "yes", c(TRUE, FALSE), 1)) {
    expect_error(ss_power_indirect_effect(.3, .3, specified_N = 100, random_predictors = bad),
                 "'random_predictors' must be TRUE or FALSE", class = "dmar_argument_invalid")
  }
})

test_that("the path powers follow the exact law of the sample correlation", {
  # With X fixed, power_a is the noncentral t power of DMAR 1.0.x.
  two_sided <- function(ncp, df) {
    cv <- stats::qt(.975, df)
    stats::pt(cv, df, ncp, lower.tail = FALSE) + stats::pt(-cv, df, ncp)
  }
  for (x in list(c(120, .3, -.25, .1), c(12, .6, .5, .2), c(65, .39, .39, 0))) {
    N <- x[1]; a <- x[2]; b <- x[3]; cp <- x[4]
    var_eY <- 1 - (b^2 + cp^2 + 2 * a * b * cp)
    rho_p <- b * sqrt(1 - a^2) / sqrt(b^2 * (1 - a^2) + var_eY)
    r <- ss_power_indirect_effect(a, b, cp, specified_N = N)
    f <- ss_power_indirect_effect(a, b, cp, specified_N = N, random_predictors = FALSE)
    expect_equal(v(r, "power_a"), .ie_r_power(a, N), tolerance = 1e-9)
    expect_equal(v(r, "power_b"), .ie_r_power(rho_p, N - 1), tolerance = 1e-9)
    expect_equal(v(f, "power_b"), v(r, "power_b"), tolerance = 1e-12)
    expect_equal(v(f, "power_a"), two_sided(a * sqrt(N / (1 - a^2)), N - 2), tolerance = 1e-12)
  }
})

test_that("power_ab is the probability that the tests reject, by an independent integral", {
  for (x in list(c(65, .39, .39, 0), c(12, .6, .5, .2), c(120, .3, -.25, .1))) {
    for (rp in c(TRUE, FALSE)) {
      r <- ss_power_indirect_effect(x[2], x[3], x[4], specified_N = x[1], random_predictors = rp)
      expect_equal(v(r, "power_ab"), .ie_power(x[1], x[2], x[3], x[4], random = rp), tolerance = 1e-9)
    }
    s <- ss_power_indirect_effect(x[2], x[3], x[4], specified_N = x[1], method = "sobel",
                                  random_predictors = FALSE)
    expect_equal(v(s, "power_ab"), .ie_power(x[1], x[2], x[3], x[4], random = FALSE, sobel = TRUE),
                 tolerance = 1e-9)
  }
  # The product of the path powers, the 1.0.x combination, overstates it.
  r <- ss_power_indirect_effect(.6, .5, .2, specified_N = 12)
  expect_gt(v(r, "power_a") * v(r, "power_b") - v(r, "power_ab"), .017)
})

test_that("reversing the sign of X, M, or Y leaves every power unchanged", {
  a <- .3; b <- -.25; cp <- .1; N <- 120
  for (m in c("joint_significance", "sobel")) for (rp in c(TRUE, FALSE)) {
    ref <- ss_power_indirect_effect(a, b, cp, specified_N = N, method = m, random_predictors = rp)
    # X reversed: a and c' change sign; M: a and b; Y: b and c'.
    for (s in list(c(-1, 1, -1), c(-1, -1, 1), c(1, -1, -1))) {
      r <- ss_power_indirect_effect(s[1] * a, s[2] * b, s[3] * cp, specified_N = N,
                                    method = m, random_predictors = rp)
      expect_equal(r$value[2:4], ref$value[2:4], tolerance = 1e-12)
    }
  }
})

test_that("the powers and the plans hold in raw data", {
  # At the plan for .80 and one case fewer, for either test and either
  # predictor model, the share of raw-data studies that reject is within 4
  # Monte Carlo standard errors of power_ab; at N = 12 the raw data agree
  # with the joint power and not with the product of the path powers.
  skip_on_cran()
  R <- 40000
  seed <- 1130
  for (m in c("joint_significance", "sobel")) for (rp in c(TRUE, FALSE)) {
    plan <- ss_power_indirect_effect(.39, .39, desired_power = .8, method = m, random_predictors = rp)
    N <- v(plan, "necessary_N")
    for (n in c(N, N - 1)) {
      seed <- seed + 1
      pk <- v(ss_power_indirect_effect(.39, .39, specified_N = n, method = m, random_predictors = rp), "power_ab")
      s <- .ie_simulate(n, .39, .39, 0, R, seed, random = rp)
      expect_identical(s[["attempted"]], R)
      expect_identical(s[["failed"]], 0)
      hits <- if (m == "sobel") s[["sobel"]] else s[["joint"]]
      p <- hits / s[["usable"]]
      mcse <- sqrt(pk * (1 - pk) / s[["usable"]])
      expect_lt(abs(p - pk), 4 * mcse)
    }
    expect_gte(v(plan, "power_ab"), .8)
  }
  s <- .ie_simulate(12, .6, .5, .2, R, 1129)
  r <- ss_power_indirect_effect(.6, .5, .2, specified_N = 12)
  p <- s[["joint"]] / s[["usable"]]
  mcse <- sqrt(p * (1 - p) / s[["usable"]])
  expect_lt(abs(p - v(r, "power_ab")), 4 * mcse)
  expect_gt(v(r, "power_a") * v(r, "power_b") - p, 6 * mcse)
  sobel <- v(ss_power_indirect_effect(.6, .5, .2, specified_N = 12, method = "sobel"), "power_ab")
  expect_lt(abs(s[["sobel"]] / s[["usable"]] - sobel), 4 * mcse)
})

test_that("ss_power_indirect_effect errors instead of hanging on a tiny effect (HIGH-04)", {
  expect_error(ss_power_indirect_effect(true_std_reg_coef_a = 1e-12, true_std_reg_coef_b = 0.39, desired_power = 0.80),
               "too small|Could not reach|did not converge")
  expect_error(ss_power_indirect_effect(1e-12, 0.39, desired_power = 0.80, method = "sobel"),
               "Could not reach")
  r <- ss_power_indirect_effect(true_std_reg_coef_a = 0.3, true_std_reg_coef_b = 0.39, desired_power = 0.80)
  expect_gte(r$value[r$term == "power_ab"], 0.80)
})

test_that("the search finds the N a scan one case at a time finds", {
  # Plans (X random) that a search one case at a time also found; the power
  # at N - 1 is below the target.
  pinned <- list(list(c(.39, .39, .8), c(69, 84)), list(c(.14, .39, .8), c(398, 412)),
                 list(c(.1, .59, .95), c(1293, 1296)), list(c(-.3, .2, .5), c(119, 157)))
  for (x in pinned) for (k in 1:2) {
    m <- c("joint_significance", "sobel")[k]
    plan <- ss_power_indirect_effect(x[[1]][1], x[[1]][2], desired_power = x[[1]][3], method = m)
    n <- v(plan, "necessary_N")
    expect_identical(n, x[[2]][k])
    expect_gte(v(plan, "power_ab"), x[[1]][3])
    expect_lt(v(ss_power_indirect_effect(x[[1]][1], x[[1]][2], specified_N = n - 1, method = m),
                "power_ab"), x[[1]][3])
  }
  # The Sobel power for paths .9 and .05 falls from N = 10 to N = 59 and
  # then rises. A target the power at 10 reaches gives N = 10; one above it
  # gives a plan past the fall.
  sob <- function(n) v(ss_power_indirect_effect(.9, .05, specified_N = n, method = "sobel"), "power_ab")
  p <- vapply(c(10, 11, 58, 59, 60), sob, 0)
  expect_true(all(diff(p[1:4]) < 0) && p[5] > p[4])
  expect_identical(v(ss_power_indirect_effect(.9, .05, desired_power = p[1], method = "sobel"),
                     "necessary_N"), 10)
  n <- v(ss_power_indirect_effect(.9, .05, desired_power = .08, method = "sobel"), "necessary_N")
  expect_gte(sob(n), .08)
  expect_lt(sob(n - 1), .08)
  # A plan in the hundreds of thousands, quickly: 405,414 for paths .0044
  # and .5 (405,411 with X fixed, the 1.0.x plan, whose power_a past 400,000
  # degrees of freedom is the exact integral).
  t0 <- proc.time()[["elapsed"]]
  plan <- ss_power_indirect_effect(.0044, .5, desired_power = .8)
  expect_lt(proc.time()[["elapsed"]] - t0, 5)
  n <- v(plan, "necessary_N")
  expect_identical(n, 405414)
  expect_lt(v(ss_power_indirect_effect(.0044, .5, specified_N = n - 1), "power_ab"), .8)
  fixed <- ss_power_indirect_effect(.0044, .5, desired_power = .8, random_predictors = FALSE)
  n <- v(fixed, "necessary_N")
  expect_identical(n, 405411)
  ncp <- .0044 / sqrt((1 - .0044^2) / n)
  cv <- stats::qt(.975, n - 2)
  expect_equal(v(fixed, "power_a"),
               .oracle_nct(cv, n - 2, ncp, lower_tail = FALSE) + .oracle_nct(-cv, n - 2, ncp),
               tolerance = 1e-10)
})

# .ie_first_n(p, N, target): the first N of a scan one case at a time, the
# smallest N whose power p reaches the target (NA if none in N).
.ie_first_n <- function(p, N, target) N[which(p >= target)[1]]

test_that("the plan is the smallest N whose power reaches the target where the Sobel power falls", {
  # Paths .8 and .02 at alpha .5: the Sobel power rises from N = 10 to 13,
  # falls to N = 59, and does not reach .507 again until N = 216, so the
  # plan for .507 is 12, not 216.
  sob <- function(n) v(ss_power_indirect_effect(.8, .02, specified_N = n, alpha_level = .5,
                                                method = "sobel"), "power_ab")
  N <- c(10:14, 59, 215, 216)
  p <- vapply(N, sob, 0)
  expect_identical(.ie_first_n(p, N, .507), 12)
  expect_true(all(diff(p[1:4]) > 0) && p[4] > p[5] && p[5] > p[6])
  expect_true(p[6] < .507 && p[7] < .507 && p[8] >= .507)
  plan <- ss_power_indirect_effect(.8, .02, desired_power = .507, alpha_level = .5, method = "sobel")
  expect_identical(v(plan, "necessary_N"), 12)
  expect_identical(v(plan, "power_ab"), p[3])
})

test_that("on a grid, the plan is the first N of a scan one case at a time from 10", {
  # Paths .8 and .05 at alpha .5 (Sobel): the power rises to .509076 at
  # N = 14, falls to .508639 at 24, and passes .509076 again at 33. Targets
  # below the rise, within it, and just past its peak; joint significance
  # in both designs. Every N below the plan falls short and the plan does
  # not.
  grid <- list(list(.8, .05, 0, .5, "sobel", TRUE, 10:40, c(.5087, .5090, .50908)),
               list(.39, .39, 0, .05, "joint_significance", TRUE, 10:200, c(.2, .5, .8, .95)),
               list(.3, -.25, .1, .2, "joint_significance", FALSE, 10:400, c(.3, .8, .95)),
               list(.6, .5, .2, .01, "joint_significance", TRUE, 10:200, c(.05, .5, .9)),
               list(-.14, .39, 0, .5, "joint_significance", FALSE, 10:800, c(.25, .8, .99)),
               list(.9, .05, 0, .8, "joint_significance", TRUE, 10:400, c(.1, .5)))
  for (g in grid) {
    p <- vapply(g[[7]], function(n) v(ss_power_indirect_effect(
      g[[1]], g[[2]], g[[3]], specified_N = n, alpha_level = g[[4]], method = g[[5]],
      random_predictors = g[[6]]), "power_ab"), 0)
    for (target in g[[8]]) {
      first <- .ie_first_n(p, g[[7]], target)
      expect_false(is.na(first))
      plan <- v(ss_power_indirect_effect(g[[1]], g[[2]], g[[3]], desired_power = target,
                                         alpha_level = g[[4]], method = g[[5]],
                                         random_predictors = g[[6]]), "necessary_N")
      expect_identical(plan, as.numeric(first))
      expect_true(all(p[g[[7]] < plan] < target))
    }
  }
})

test_that("the bound the search passes over N with is at least the power at each N", {
  bound <- function(a, b, lo, hi, al, m) {
    DMAR:::.indirect_power(hi, abs(a) / sqrt(1 - a^2), abs(b) * sqrt(1 - a^2) / sqrt(1 - b^2),
                           al, m, TRUE, lo)[1]
  }
  # Every N of each block; of the Sobel block of 601 N, longer than the
  # bound's 256 runs, every 20th N.
  for (x in list(list(.8, .02, 11, 14, .5, "sobel", 1), list(.8, .02, 13, 30, .5, "sobel", 1),
                 list(.9, .05, 20, 620, .05, "sobel", 20),
                 list(.3, -.25, 50, 80, .05, "joint_significance", 1),
                 list(.3, -.25, 100, 400, .05, "joint_significance", 1))) {
    p <- vapply(seq(x[[3]], x[[4]], by = x[[7]]), function(n) v(ss_power_indirect_effect(
      x[[1]], x[[2]], specified_N = n, alpha_level = x[[5]], method = x[[6]]), "power_ab"), 0)
    expect_gte(bound(x[[1]], x[[2]], x[[3]], x[[4]], x[[5]], x[[6]]), max(p) - 1e-9)
    expect_equal(bound(x[[1]], x[[2]], x[[4]], x[[4]], x[[5]], x[[6]]), p[length(p)],
                 tolerance = 1e-14)
  }
})

test_that("the bound for a near-zero path is at least the power, and the search uses it", {
  # The bound from its formula (Details): the smaller, over the two paths, of
  # the rate at which the path's test rejects when the path is zero
  # (alpha_level; for Sobel, P(|t| > z) on N - 2 or N - 3 df) plus dnorm(1)
  # times the expected squared noncentrality.
  nb <- function(a, b, N, al, m, rp) {
    z <- stats::qnorm(1 - al / 2)
    ra <- if (m == "sobel") 2 * stats::pt(-z, N - 2) else al
    rb <- if (m == "sobel") 2 * stats::pt(-z, N - 3) else al
    min(ra + stats::dnorm(1) * a^2 / (1 - a^2) * (if (rp) N - 1 else N),
        rb + stats::dnorm(1) * b^2 * (1 - a^2) / (1 - b^2) * (N - 2))
  }
  for (x in list(list(.01, .5, .05, "joint_significance", TRUE),
                 list(.6, -.003, .3, "joint_significance", FALSE),
                 list(.02, .4, .05, "sobel", TRUE), list(.5, .01, .2, "sobel", FALSE))) {
    for (N in c(10, 40, 400, 4000)) {
      p <- v(ss_power_indirect_effect(x[[1]], x[[2]], specified_N = N, alpha_level = x[[3]],
                                      method = x[[4]], random_predictors = x[[5]]), "power_ab")
      expect_lte(p, nb(x[[1]], x[[2]], N, x[[3]], x[[4]], x[[5]]))
    }
  }
  # Paths .02 and .5: the bound shows that no N up to 206 reaches .07, and
  # the plan is the first N of a scan, 435.
  expect_lt(nb(.02, .5, 206, .05, "joint_significance", TRUE), .07)
  N <- 10:450
  p <- vapply(N, function(n) v(ss_power_indirect_effect(.02, .5, specified_N = n), "power_ab"), 0)
  expect_identical(.ie_first_n(p, N, .07), 435L)
  expect_identical(v(ss_power_indirect_effect(.02, .5, desired_power = .07), "necessary_N"), 435)
  # A path of 1e-12 cannot give power .050001 below ten million cases: joint
  # significance refuses without computing a power, and Sobel (at .0501)
  # after a few.
  calls <- 0
  power <- DMAR:::.indirect_power
  local_mocked_bindings(.indirect_power = function(...) {
    calls <<- calls + 1
    power(...)
  })
  expect_error(ss_power_indirect_effect(1e-12, .39, desired_power = .050001),
               class = "dmar_computation_failed")
  expect_error(ss_power_indirect_effect(.39, 1e-12, desired_power = .050001),
               class = "dmar_computation_failed")
  expect_identical(calls, 0)
  expect_error(ss_power_indirect_effect(1e-12, .39, desired_power = .0501, method = "sobel"),
               class = "dmar_computation_failed")
  expect_lt(calls, 100)
})

test_that("a failed Sobel integral names the N, or the range of N a bound was for", {
  local_mocked_bindings(integrate = function(...) stop("no convergence"), .package = "stats")
  expect_error(ss_power_indirect_effect(.39, .39, specified_N = 40, method = "sobel"),
               "The power of the Sobel test at N = 40 could not", class = "dmar_computation_failed")
  expect_error(DMAR:::.indirect_power(20, .4, .4, .05, "sobel", TRUE, 12),
               "A bound of the power of the Sobel test for N from 12 to 20 could not",
               class = "dmar_computation_failed")
})

test_that("the Sobel power at alpha_level = .999 is computed, not refused", {
  # Pieces of the integral next to z, worth less than 1e-45, come back from
  # integrate() with "roundoff error is detected in the extrapolation
  # table" and an estimated error near 1e-47, and the call was refused
  # with that status until such a piece was accepted within its error
  # budget. X fixed, against the independent integral above; X random,
  # against an independent nested integral (integrate() over both chi
  # squares, pt() for t_b; too slow to run here), which gives
  # .998264864866.
  f <- ss_power_indirect_effect(.2, .1, specified_N = 10, alpha_level = .999,
                                method = "sobel", random_predictors = FALSE)
  expect_equal(v(f, "power_ab"),
               .ie_power(10, .2, .1, alpha = .999, random = FALSE, sobel = TRUE),
               tolerance = 1e-9)
  r <- ss_power_indirect_effect(.2, .1, specified_N = 10, alpha_level = .999, method = "sobel")
  expect_equal(v(r, "power_ab"), .998264864866, tolerance = 1e-9)
  for (s in list(c(-1, 1), c(1, -1), c(-1, -1))) {
    expect_equal(v(ss_power_indirect_effect(s[1] * .2, s[2] * .1, specified_N = 10,
                                            alpha_level = .999, method = "sobel"), "power_ab"),
                 v(r, "power_ab"), tolerance = 1e-12)
  }
})

test_that("the Sobel power is computed where the integral starts within rounding of z", {
  # At alpha_level = .99999, z = 1.25e-5, and at N = 10 k_top is above
  # 3,000 at every W node, so at two of them y_lo = k(k_top) rounds just
  # below z. y^2 - z^2 was then negative and the integrand NaN, with the
  # warning "NaNs produced"; once roundoff pieces were accepted, the call
  # was refused with "non-finite function value" (before, a roundoff
  # piece refused it first). The test cannot reject at |t_a| <= z,
  # so k is infinite there and the integrand 0. X fixed, against the
  # independent integral above; X random, against the independent nested
  # integral named in the test at alpha_level = .999, which gives
  # .999989909203 here.
  f <- expect_no_warning(ss_power_indirect_effect(
    .39, .39, specified_N = 10, alpha_level = .99999, method = "sobel",
    random_predictors = FALSE))
  expect_equal(v(f, "power_ab"),
               .ie_power(10, .39, .39, alpha = .99999, random = FALSE, sobel = TRUE),
               tolerance = 1e-10)
  r <- expect_no_warning(ss_power_indirect_effect(
    .39, .39, specified_N = 10, alpha_level = .99999, method = "sobel"))
  expect_equal(v(r, "power_ab"), .999989909203, tolerance = 1e-10)
  p <- expect_no_warning(ss_power_indirect_effect(
    .39, .39, desired_power = .999, alpha_level = .99999, method = "sobel"))
  expect_identical(v(p, "necessary_N"), 10)
})

test_that("with X random and a near 1, the Sobel integral finds each peak of the density of t_a", {
  # Given W, the density of |t_a| is a mixture over the SS_X nodes of
  # normal densities with a common sd s1, and with a = .999999 their means
  # stand up to about 700 s1 apart. integrate() stepped over some of the
  # peaks: the power was .352 at alpha_level = .5, so the plan for .52
  # stopped at ten million, and .377 at alpha_level = .999. Against an
  # independent integral over W and SS_X by their quantiles, given which
  # t_a is normal with a single peak, with pt() for t_b (too slow to run
  # here), which gives .5216460213764 and .9990349721533.
  r <- ss_power_indirect_effect(.999999, .5, specified_N = 10, alpha_level = .5,
                                method = "sobel")
  expect_equal(v(r, "power_ab"), .5216460213764, tolerance = 1e-10)
  r <- ss_power_indirect_effect(.999999, .5, specified_N = 10, alpha_level = .999,
                                method = "sobel")
  expect_equal(v(r, "power_ab"), .9990349721533, tolerance = 1e-10)
  p <- ss_power_indirect_effect(.999999, .5, desired_power = .52, alpha_level = .5,
                                method = "sobel")
  expect_identical(v(p, "necessary_N"), 10)
})

test_that("a piece of the Sobel integral returned OK is accepted whatever its error estimate", {
  # .99931434686885845 is the power before roundoff statuses were accepted;
  # a rule that also refused OK pieces over a budget refused this call.
  r <- ss_power_indirect_effect(.39, .39, specified_N = 200, alpha_level = .05, method = "sobel")
  expect_equal(v(r, "power_ab"), .99931434686885845, tolerance = 1e-12)
  real <- stats::integrate
  ref <- v(ss_power_indirect_effect(.39, .39, specified_N = 40, method = "sobel"), "power_ab")
  big_error <- function(...) {
    q <- real(...)
    q$abs.error <- 1
    q
  }
  local_mocked_bindings(integrate = big_error, .package = "stats")
  expect_identical(v(ss_power_indirect_effect(.39, .39, specified_N = 40, method = "sobel"),
                     "power_ab"), ref)
})

test_that("a roundoff status is accepted within its error budget, and every other failure refused", {
  real <- stats::integrate
  ref <- v(ss_power_indirect_effect(.39, .39, specified_N = 40, method = "sobel"), "power_ab")
  # integrate() as it is, but with the status, error estimate, or value of
  # each piece replaced.
  as_if <- function(message, abs.error = NULL, value = NULL) function(...) {
    q <- real(...)
    q$message <- message
    if (!is.null(abs.error)) q$abs.error <- abs.error
    if (!is.null(value)) q$value <- value
    q
  }
  sobel_40 <- function(fake) with_mocked_bindings(
    v(ss_power_indirect_effect(.39, .39, specified_N = 40, method = "sobel"), "power_ab"),
    integrate = fake, .package = "stats")
  for (status in c("roundoff error was detected",
                   "roundoff error is detected in the extrapolation table")) {
    expect_identical(sobel_40(as_if(status, abs.error = 0)), ref)
    expect_error(sobel_40(as_if(status, abs.error = 1e-6)),
                 paste0("^The power of the Sobel test at N = 40 could not be computed reliably: ",
                        "the numerical integral over the distribution of the path t statistics ",
                        "failed \\(", status, ", with an estimated absolute error of 1e-06 ",
                        "against a budget of [0-9.]+e-[0-9]+ for this piece\\)\\. ",
                        "The joint significance test \\(method = \"joint_significance\"\\) ",
                        "does not need this integral\\.$"),
                 class = "dmar_computation_failed")
    expect_error(sobel_40(as_if(status, abs.error = NaN)),
                 "estimated absolute error of NaN", class = "dmar_computation_failed")
    expect_error(sobel_40(as_if(status, abs.error = 0, value = -1)),
                 "and a value of -1 for the integral of a nonnegative function\\)\\. The joint",
                 class = "dmar_computation_failed")
  }
  for (status in c("maximum number of subdivisions reached",
                   "the integral is probably divergent")) {
    expect_error(sobel_40(as_if(status, abs.error = 0)),
                 paste0("path t statistics failed \\(", status, "\\)\\.$"),
                 class = "dmar_computation_failed")
  }
})

test_that("tidy() and glance() read power_ab as the power of the design", {
  r <- ss_power_indirect_effect(true_std_reg_coef_a = .39,
                                true_std_reg_coef_b = .39, specified_N = 75)
  expect_equal(generics::tidy(r)$power, v(r, "power_ab"))
  g <- generics::glance(r)
  expect_equal(g$power, v(r, "power_ab"))
  expect_equal(c(g$power_a, g$power_b), c(v(r, "power_a"), v(r, "power_b")))
  expect_false("power_ab" %in% names(g))
  # The Sobel method has no component powers; tidy() still reads power_ab.
  s <- ss_power_indirect_effect(true_std_reg_coef_a = .39,
                                true_std_reg_coef_b = .39, specified_N = 75,
                                method = "sobel")
  expect_equal(generics::tidy(s)$power, v(s, "power_ab"))
  expect_true(is.finite(generics::tidy(s)$power))
})

test_that("a call gives no warning, and the power at N = 65 for paths of .39 is .7675", {
  expect_no_warning(r <- ss_power_indirect_effect(.39, .39, specified_N = 65))
  expect_identical(sprintf("%.4f", v(r, "power_ab")), "0.7675")
  expect_no_warning(ss_power_indirect_effect(.39, .39, specified_N = 65, method = "sobel",
                                             random_predictors = FALSE))
  expect_error(ss_power_indirect_effect(0, .3, specified_N = 50), "nonzero")
})
