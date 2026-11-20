## .pbinorm() -- DMAR's bivariate normal distribution function, which
## replaced mvtnorm::pmvnorm() in categorical omega -- and the vectorized
## Green and Yang (2009) assembly, .gy_omega_from_components(), that
## calls it. Nothing here needs lavaan or mvtnorm: the mvtnorm values are
## stored, and tools/oracle_snippets/mvtnorm.R makes the live comparison.

test_that(".gauss_legendre() gives the Gauss-Legendre rules", {
  # The six-point rule as tabulated, to 16 digits.
  six <- .gauss_legendre(6L)
  expect_equal(sort(six$x),
               c(-0.9324695142031521, -0.6612093864662645,
                 -0.2386191860831969, 0.2386191860831969,
                 0.6612093864662645, 0.9324695142031521),
               tolerance = 1e-14)
  expect_equal(six$w[order(six$x)],
               c(0.1713244923791704, 0.3607615730481386,
                 0.4679139345726910, 0.4679139345726910,
                 0.3607615730481386, 0.1713244923791704),
               tolerance = 1e-14)
  # An n-point rule integrates x^m over [-1, 1] exactly for m up to
  # 2n - 1: 2 / (m + 1) for even m and 0 for odd m.
  for (n in c(6L, 12L, 20L)) {
    rule <- .gauss_legendre(n)
    m <- 0:(2L * n - 1L)
    exact <- ifelse(m %% 2L == 0L, 2 / (m + 1), 0)
    quad <- vapply(m, function(p) sum(rule$w * rule$x^p), numeric(1))
    expect_lt(max(abs(quad - exact)), 1e-14)
  }
  # The rules .pbinorm() reads are these, computed at installation.
  expect_identical(.GAUSS_LEGENDRE,
                   lapply(c(6L, 12L, 20L), .gauss_legendre))
})

test_that(".pbinorm() reproduces mvtnorm::pmvnorm()", {
  # Pinned from mvtnorm::pmvnorm(upper = c(h, k), corr = ...) (mvtnorm
  # 1.4.2, 2026-11-01); the live comparison, over these points and a
  # grid of 2,873, is in tools/oracle_snippets/mvtnorm.R. Twenty-four
  # points in each of the four ranges of |r| in which the algorithm
  # changes (6, 12, or 20 nodes below .925, the series above it), then
  # the six boundaries themselves; only the reference values are
  # stored. The largest difference here is 1.1e-16 (4.0e-16 on the
  # grid); the tolerance leaves room for another platform's exp(),
  # sin(), and pnorm(), not for any change in the method.
  oracle <- c(
    0.00042464521235648843, 0.24019082896897503, 0.10058048057338051,
    0.97448826552351908, 0.005715580789979623, 0.82926698124385501,
    0.014796615157071667, 0.56306609772883176, 0.15344120558816407,
    0.27757324240767051, 0.25366766976966243, 0.0066527650383515526,
    0.029226762665221166, 0.038003687556804763, 0.13407332213039874,
    0.8771563760542358, 0.042721181638332127, 0.00044966118401549887,
    0.030308177176576943, 0.029721138281604772, 0.0019811474683126551,
    0.3203176330585818, 0.10615994487292044, 0.096205337523475501,
    0.0032392653863809301, 1.8892548398304443e-06, 0.15708492539903615,
    0.79076367336801512, 0.013036449302394849, 2.3051808829194516e-06,
    0.50339839446721313, 0.0054266375514223416, 0.073448141652885393,
    0.68416139468186687, 0.57094039752495251, 0.25289980116827854,
    0.30221175570295838, 0.0064650444557325694, 0.55701510704041124,
    0.84793501383004044, 0.021227226356511614, 0.12104721080053733,
    0.98713375218751365, 0.2090887228852904, 0.016200676811855635,
    0.80434509394195675, 0.49894246044035451, 0.1926403324482345,
    0.013639550706805456, 0.035962450316314259, 0.031736185073252224,
    0.012320956356470547, 0.0018093393055732946, 0.48899576266795858,
    0.0045446903322128758, 0.58168747643413377, 0.34885898628312989,
    0.12525756013344341, 0.84632453780550898, 0.020728058648403817,
    1.0345566949316967e-18, 0.80427843449824787, 1.5789243729067481e-07,
    0.14814660901131724, 0.75774565939432936, 0.23093473631443429,
    0.004338487208774537, -6.434096135825073e-20, 0.47507586228921289,
    0.078437628024841041, 5.686301240686318e-22, 0.0053426157476266587,
    0.045885034591996519, 0.97416569413165244, 0.72774154170712535,
    0.05575533023599568, 0.15810103039067031, 3.9319109200571753e-15,
    0.0048704560785367045, 0.45786793357443195, 0.39457936688784173,
    0.0037495844183093709, 1.5712216484777855e-19, 7.3710681159310816e-06,
    0.26110223904704344, 0.97982965946671974, 3.3867186919656293e-68,
    0.65708034809006954, 0.030063724909587187, 0.14133972122330615,
    0.0078672153144958151, 0.031308629787635475, 0.59076874216583031,
    0.52088046285387368, 3.4020039061281034e-05, 0.67266914383231047,
    0.0016444659152474079, 0.018393198975041148, 0.063272380865661629,
    0.11126216320754038, 0.13397222073588291, 0.13566422338921663)
  set.seed(113)
  sgn <- function(n) sample(c(-1, 1), n, replace = TRUE)
  r <- c(runif(24, -0.3, 0.3), sgn(24) * runif(24, 0.3, 0.75),
         sgn(24) * runif(24, 0.75, 0.925), sgn(24) * runif(24, 0.925, 0.9999),
         c(-0.925, -0.75, -0.3, 0.3, 0.75, 0.925))
  h <- c(runif(96, -3, 3), rep(0.4, 6))
  k <- c(runif(96, -3, 3), rep(-1.1, 6))
  expect_lt(max(abs(.pbinorm(h, k, r) - oracle)), 2e-15)
  # One call over all the points equals a call per point: the branches
  # do not leak into each other.
  one_at_a_time <- mapply(.pbinorm, h, k, r)
  expect_identical(.pbinorm(h, k, r), one_at_a_time)
})

test_that(".pbinorm() is exact at its limits and where a closed form exists", {
  h <- c(-2.2, -0.7, 0, 0.35, 1.6)
  k <- c(0.9, -1.3, 0.2, 0.35, -0.4)
  # r = 0: the product of the margins, exactly.
  expect_identical(.pbinorm(h, k, 0), stats::pnorm(h) * stats::pnorm(k))
  # r = 1: the smaller margin; r = -1: max(0, Phi(h) + Phi(k) - 1).
  expect_identical(.pbinorm(h, k, 1), stats::pnorm(pmin(h, k)))
  expect_lt(max(abs(.pbinorm(h, k, -1) -
                      pmax(0, stats::pnorm(h) + stats::pnorm(k) - 1))),
            1e-16)
  # At h = k = 0 the orthant probability is 1/4 + asin(r) / (2 pi) for
  # every r. Near |r| = 1 this is where mvtnorm 1.4.2 falls back to the
  # singular limit (off by 1.3e-6 at 1 - |r| = 3e-11, and by 4e-14 at
  # 1 - |r| = 1e-7); the port stays exact.
  rr <- c(-1, 1) * rep(1 - 10^-c(1, 4, 7, 9, 10.5, 12, 15), each = 2)
  rr <- c(rr, seq(-0.95, 0.95, by = 0.05))
  expect_lt(max(abs(.pbinorm(0, 0, rr) - (0.25 + asin(rr) / (2 * pi)))),
            2e-16)
  # Symmetric in the two limits, and the four orthants sum to one:
  # P(X <= h, Y <= k) + P(X <= h, Y > k) + P(X > h, Y <= k) +
  # P(X > h, Y > k), with P(X <= h, Y > k) = P(X <= h, -Y < -k) at -r.
  r <- c(-0.97, -0.6, -0.1, 0.2, 0.8, 0.95)
  for (rho in r) {
    expect_lt(max(abs(.pbinorm(h, k, rho) - .pbinorm(k, h, rho))), 1e-16)
    total <- .pbinorm(h, k, rho) + .pbinorm(h, -k, -rho) +
      .pbinorm(-h, k, -rho) + .pbinorm(-h, -k, rho)
    expect_lt(max(abs(total - 1)), 1e-15)
  }
  # An infinite limit gives the margin or the bound.
  expect_identical(.pbinorm(c(Inf, Inf, -Inf, 0.3, Inf),
                            c(0.3, Inf, 0.3, -Inf, -Inf), 0.5),
                   c(stats::pnorm(0.3), 1, 0, 0, 0))
  # Near the singular end, against the defining integral
  # int_{-Inf}^{h} phi(x) Phi((k - r x) / sqrt(1 - r^2)) dx, taken by
  # stats::integrate() split where the inner distribution function
  # steps (the values are stored; mvtnorm 1.4.2 is off by 2e-11 at the
  # first point).
  expect_lt(abs(.pbinorm(0.6665284, 0.7633714, -0.99999999991364019) -
                  0.52484232570225375), 1e-15)
  expect_lt(abs(.pbinorm(-0.6900522, -0.6895479, 0.99999995549354481) -
                  0.24507892381075119), 1e-15)
})

test_that(".pbinorm() refuses what mvtnorm refused, and takes a bound within rounding", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  expect_error(.pbinorm(0.2, NA, 0.5), "missing limit or correlation")
  expect_error(.pbinorm(0.2, 0.1, NA_real_), "missing limit or correlation")
  # A correlation above one is what an improper (Heywood) solution
  # would send; categorical omega refuses such a solution before the
  # assembly, so this is a backstop, as mvtnorm's error was.
  expect_error(.pbinorm(0.2, 0.1, 1.2), "needs correlations in \\[-1, 1\\]")
  expect_error(.pbinorm(0.2, 0.1, -1 - 1e-6),
               "needs correlations in \\[-1, 1\\]")
  # Within sqrt(.Machine$double.eps) of a bound the correlation is the
  # bound (mvtnorm 1.4.2 accepts it too, but returns 0 at 1 + 1e-10).
  expect_identical(.pbinorm(0.3, -0.2, 1 + 1e-10), stats::pnorm(-0.2))
  expect_identical(.pbinorm(numeric(0), 1, 0.5), numeric(0))
  # More probabilities than one block (65,536) are computed block by
  # block with the same result.
  n <- 70000L
  h <- seq(-3, 3, length.out = n)
  k <- rev(h)
  r <- rep(c(-0.95, -0.5, 0.1, 0.6, 0.93), length.out = n)
  whole <- .pbinorm(h, k, r)
  parts <- c(.pbinorm(h[1:40000], k[1:40000], r[1:40000]),
             .pbinorm(h[40001:n], k[40001:n], r[40001:n]))
  expect_identical(whole, parts)
})

# The assembly ----------------------------------------------------------

# The assembly as it was written through version 1.0.1, four nested
# loops with one bivariate normal call per cell, here with .pbinorm() as
# the call and with item weights, so the vectorized assembly can be
# checked against it cell for cell.
gy_loop <- function(thresh, rho_num, rho_den,
                    weights = rep(1, length(thresh))) {
  J <- length(thresh)
  num <- 0
  den <- 0
  for (j in seq_len(J)) {
    for (jp in seq_len(J)) {
      t1 <- thresh[[j]]
      t2 <- thresh[[jp]]
      p2_num <- 0
      p2_den <- 0
      for (c in seq_along(t1)) {
        for (cp in seq_along(t2)) {
          p2_num <- p2_num + .pbinorm(t1[c], t2[cp], rho_num[j, jp])
          p2_den <- p2_den + .pbinorm(t1[c], t2[cp], rho_den[j, jp])
        }
      }
      product <- sum(stats::pnorm(t1)) * sum(stats::pnorm(t2))
      num <- num + weights[j] * weights[jp] * (p2_num - product)
      den <- den + weights[j] * weights[jp] * (p2_den - product)
    }
  }
  num / den
}

# A random single-factor configuration: J items with between one and
# six thresholds each (two to seven categories), loadings that reach
# the steep range of the bivariate normal (|r| above .925) and include
# a negative one, and a denominator matrix that departs from the model
# (a doublet on the first two items), with a unit diagonal.
gy_config <- function(J) {
  thresh <- lapply(seq_len(J), function(j) {
    sort(stats::runif(sample(1:6, 1), -2.2, 2.2))
  })
  lambda <- stats::runif(J, 0.35, 0.98) * sample(c(1, 1, 1, -1), J, TRUE)
  rho_num <- outer(lambda, lambda)
  diag(rho_num) <- lambda^2
  rho_den <- outer(lambda, lambda)
  doublet <- sign(rho_den[1, 2]) * min(abs(rho_den[1, 2]) + 0.1, 0.995)
  rho_den[1, 2] <- rho_den[2, 1] <- doublet
  diag(rho_den) <- 1
  list(thresh = thresh, rho_num = rho_num, rho_den = rho_den)
}

test_that("the vectorized assembly equals the loops it replaced", {
  skip_on_cran()
  set.seed(113)
  for (i in 1:12) {
    cfg <- gy_config(sample(3:7, 1))
    expect_equal(.gy_omega_from_components(cfg$thresh, cfg$rho_num,
                                           cfg$rho_den),
                 gy_loop(cfg$thresh, cfg$rho_num, cfg$rho_den),
                 tolerance = 1e-14)
  }
  # Steep correlations off the diagonal, where the series branch of
  # .pbinorm() does the work.
  thresh <- list(c(-1, 0.4), c(-0.5, 0.6, 1.5), c(0.1), c(-1.8, -0.2, 0.9))
  lambda <- c(0.99, 0.97, -0.98, 0.965)
  rho_num <- outer(lambda, lambda)
  rho_den <- rho_num
  diag(rho_den) <- 1
  expect_equal(.gy_omega_from_components(thresh, rho_num, rho_den),
               gy_loop(thresh, rho_num, rho_den), tolerance = 1e-14)
})

test_that("the categorical omega of Green and Yang's (2009) designs is unchanged from the mvtnorm loops", {
  # The ten population configurations of Green and Yang (2009, Table 1)
  # that test-reliability_omega_categorical.R checks against the
  # published values, here pinned from the loops of version 1.0.1,
  # which called mvtnorm::pmvnorm() (mvtnorm 1.4.2) once per cell. The
  # inputs are exact (normal quantiles of fixed category proportions,
  # products of fixed loadings), so the comparison does not depend on a
  # model fit; the largest difference is 8.9e-16.
  thr <- function(P) stats::qnorm(cumsum(P)[-length(P)])
  gy <- function(P_list, lambda) {
    rho_num <- outer(lambda, lambda)
    rho_den <- rho_num
    diag(rho_den) <- 1
    .gy_omega_from_components(lapply(P_list, thr), rho_num, rho_den)
  }
  u2 <- c(.5, .5); neg2 <- c(.10, .90); pos2 <- rev(neg2)
  u5 <- rep(.2, 5); bell <- c(.09, .15, .51, .15, .09)
  matched <- c(.485, .01, .01, .01, .485)
  smooth_neg <- c(.011, .027, .047, .118, .797); smooth_pos <- rev(smooth_neg)
  l8 <- rep(.8, 4); l85 <- c(.8, .8, .5, .5)
  computed <- c(
    gy(rep(list(u2), 4), l8), gy(rep(list(u5), 4), l8),
    gy(rep(list(bell), 4), l8), gy(rep(list(matched), 4), l8),
    gy(rep(list(neg2), 4), l8), gy(list(neg2, neg2, pos2, pos2), l8),
    gy(rep(list(smooth_neg), 4), l8),
    gy(list(smooth_neg, smooth_neg, smooth_pos, smooth_pos), l8),
    gy(rep(list(u2), 4), l85), gy(rep(list(u5), 4), l85))
  mvtnorm_loops <- c(0.76020026004737684, 0.85023756061651168,
                     0.83970425333128162, 0.76815156508193849,
                     0.68778218573069116, 0.59042241893746183,
                     0.78646814803241005, 0.71631179221373398,
                     0.61735380858811773, 0.71976245758355861)
  expect_lt(max(abs(computed - mvtnorm_loops)), 2e-15)
})

test_that("weights give the reliability of a weighted composite of the category scores", {
  skip_on_cran()
  set.seed(113)
  cfg <- gy_config(5)
  unit <- .gy_omega_from_components(cfg$thresh, cfg$rho_num, cfg$rho_den)
  # Unit weights are the unweighted sum score, to the bit; doubling
  # every weight multiplies numerator and denominator by four, exactly.
  expect_identical(.gy_omega_from_components(cfg$thresh, cfg$rho_num,
                                             cfg$rho_den, rep(1, 5)), unit)
  w <- c(2, 1, 0.5, 1.5, 1)
  weighted <- .gy_omega_from_components(cfg$thresh, cfg$rho_num,
                                        cfg$rho_den, w)
  expect_identical(.gy_omega_from_components(cfg$thresh, cfg$rho_num,
                                             cfg$rho_den, 2 * w), weighted)
  expect_equal(weighted, gy_loop(cfg$thresh, cfg$rho_num, cfg$rho_den, w),
               tolerance = 1e-14)
  expect_false(isTRUE(all.equal(weighted, unit)))
  # A weight of zero leaves the item out of the composite.
  expect_equal(.gy_omega_from_components(cfg$thresh, cfg$rho_num,
                                         cfg$rho_den, c(1, 1, 0, 1, 1)),
               .gy_omega_from_components(cfg$thresh[-3],
                                         cfg$rho_num[-3, -3],
                                         cfg$rho_den[-3, -3]),
               tolerance = 1e-14)
  expect_error(.gy_omega_from_components(cfg$thresh, cfg$rho_num,
                                         cfg$rho_den, c(1, 1, 1)),
               "must be 5 finite numbers, one per item")
  expect_error(.gy_omega_from_components(cfg$thresh, cfg$rho_num,
                                         cfg$rho_den, c(1, NA, 1, 1, 1)),
               "must be 5 finite numbers, one per item")
})

test_that("the assembly refuses a communality above one through .pbinorm()", {
  # An improper solution never reaches the assembly (both callers refuse
  # it first); were one to, the bivariate normal at a correlation above
  # one stops it rather than returning a number.
  thresh <- list(c(-0.5, 0.5), c(0), c(-1, 1))
  rho_num <- matrix(0.5, 3, 3)
  diag(rho_num) <- c(1.1, 0.5, 0.5)
  rho_den <- matrix(0.5, 3, 3)
  diag(rho_den) <- 1
  expect_error(.gy_omega_from_components(thresh, rho_num, rho_den),
               "needs correlations in \\[-1, 1\\]; it was given 1.1")
})
