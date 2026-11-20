# Tests for the Bryant-Paulson generalized studentized range and the
# simultaneous ANCOVA confidence intervals it produces. The published
# critical values of Bryant & Paulson (1976, Table 1) and Bryant & Bruvold
# (1980, Tables 1-2) provide exact targets.

test_that("qbryant_paulson reproduces the published critical value q_.05;1,6,14 = 4.83", {
  expect_equal(qbryant_paulson(0.95, n_covariates = 1, n_groups = 6, df = 14),
               4.83, tolerance = 5e-3)
})

test_that("p = 0 reduces exactly to the ordinary studentized range (Tukey)", {
  for (k in c(2, 4, 6)) for (nu in c(10, 20, 60)) {
    expect_equal(qbryant_paulson(0.95, n_covariates = 0, n_groups = k, df = nu),
                 qtukey(0.95, nmeans = k, df = nu), tolerance = 1e-6)
    expect_equal(pbryant_paulson(3.5, n_covariates = 0, n_groups = k, df = nu),
                 ptukey(3.5, nmeans = k, df = nu), tolerance = 1e-8)
  }
})

test_that("Bryant-Paulson critical values exceed Tukey's (random-covariate inflation)", {
  skip_on_cran()  # repeated numeric integration; fast anchors remain on CRAN
  for (p in 1:3) {
    bp <- qbryant_paulson(0.95, n_covariates = p, n_groups = 5, df = 25)
    tk <- qtukey(0.95, nmeans = 5, df = 25)
    expect_gt(bp, tk)
  }
  # More covariates -> larger critical value.
  q123 <- qbryant_paulson(0.95, n_covariates = 1:3, n_groups = 4, df = 20)
  expect_true(all(diff(q123) > 0))
})

test_that("k = 2 entries of Bryant & Bruvold (1980) Table 2 are reproduced", {
  skip_on_cran()  # repeated numeric integration; fast anchors remain on CRAN
  # Table 2, k = 2 column (where the Duncan value equals the BP q value).
  expect_equal(qbryant_paulson(0.95, 1, 2, 20),  3.03, tolerance = 5e-3) # p=1 a=.05
  expect_equal(qbryant_paulson(0.99, 1, 2, 20),  4.14, tolerance = 5e-3) # p=1 a=.01
  expect_equal(qbryant_paulson(0.95, 2, 2, 20),  3.10, tolerance = 5e-3) # p=2 a=.05
  expect_equal(qbryant_paulson(0.99, 3, 2, 120), 3.75, tolerance = 5e-3) # p=3 a=.01
})

test_that("p and q functions are mutual inverses", {
  skip_on_cran()  # repeated numeric integration; fast anchors remain on CRAN
  q <- qbryant_paulson(0.90, n_covariates = 2, n_groups = 4, df = 30)
  expect_equal(pbryant_paulson(q, 2, 4, 30), 0.90, tolerance = 1e-5)
})

test_that("lower_tail toggles the complementary probability", {
  expect_equal(pbryant_paulson(4, 1, 6, 14, lower_tail = FALSE),
               1 - pbryant_paulson(4, 1, 6, 14), tolerance = 1e-10)
})

test_that("invalid arguments are rejected", {
  expect_error(qbryant_paulson(0.95, n_covariates = -1, n_groups = 4, df = 10))
  expect_error(qbryant_paulson(0.95, n_covariates = 1,  n_groups = 1, df = 10))
  expect_error(qbryant_paulson(0.95, n_covariates = 1,  n_groups = 4, df = 0))
  expect_error(qbryant_paulson(1.5,  n_covariates = 1,  n_groups = 4, df = 10))
})

test_that("dbryant_paulson is a valid density: non-negative and integrates to one", {
  skip_on_cran()
  qs <- seq(0.5, 8, by = 0.5)
  dens <- dbryant_paulson(qs, n_covariates = 1, n_groups = 4, df = 20)
  expect_true(all(dens >= 0))
  area <- integrate(function(q) dbryant_paulson(q, 1, 4, 20),
                    lower = 0, upper = 30)$value
  expect_equal(area, 1, tolerance = 1e-4)
})

test_that("dbryant_paulson at p = 0 matches the derivative of the studentized range", {
  q0 <- 3.5
  h  <- 1e-4
  fd <- (ptukey(q0 + h, nmeans = 4, df = 20) -
         ptukey(q0 - h, nmeans = 4, df = 20)) / (2 * h)
  expect_equal(dbryant_paulson(q0, n_covariates = 0, n_groups = 4, df = 20),
               fd, tolerance = 1e-6)
})

test_that("ci_c_ancova_bp reproduces the paper's 0.278 critical difference", {
  adj <- c(3.595, 3.619, 4.102, 4.515, 4.618, 4.876)
  ci <- ci_c_ancova_bp(adj_means = adj, s_ancova = sqrt(0.01326), n_per_group = 4,
                       n_covariates = 1, df = 14)
  expect_equal(nrow(ci), choose(6, 2))
  half <- (ci$ci_upper - ci$ci_lower) / 2
  expect_true(all(abs(half - 0.278) < 1e-3))
  expect_equal(attr(ci, "critical_value"),
               qbryant_paulson(0.95, 1, 6, 14), tolerance = 1e-8)
})

test_that("ci_c_ancova_bp point estimates equal the contrasts of adjusted means", {
  skip_on_cran()  # a root find for an arithmetic check; the 0.278 anchor runs on CRAN
  adj <- c(3.595, 3.619, 4.102, 4.515, 4.618, 4.876)
  ci <- ci_c_ancova_bp(adj_means = adj, s_ancova = sqrt(0.01326), n_per_group = 4, df = 14,
                       c_weights = c(1, -1, 0, 0, 0, 0))
  expect_equal(ci$estimate, adj[1] - adj[2], tolerance = 1e-10)
})

test_that("pairwise and allowance widths agree for a pairwise contrast", {
  skip_on_cran()  # two root finds for one width; the 0.278 anchor above runs on CRAN
  adj <- c(1, 2, 3, 4)
  a <- ci_c_ancova_bp(adj_means = adj, s_ancova = 2, n_per_group = 10, df = 35,
                      c_weights = c(1, -1, 0, 0), contrast_type = "pairwise")
  b <- ci_c_ancova_bp(adj_means = adj, s_ancova = 2, n_per_group = 10, df = 35,
                      c_weights = c(1, -1, 0, 0), contrast_type = "allowance")
  expect_equal(a$ci_lower, b$ci_lower, tolerance = 1e-8)
})

test_that("default df for a one-way ANCOVA is N - k - p", {
  skip_on_cran()  # two root finds for a default check; the 0.278 anchor above runs on CRAN
  adj <- c(1, 2, 3)
  ci <- ci_c_ancova_bp(adj_means = adj, s_ancova = 1, n_per_group = 10, n_covariates = 1)
  # N = 30, k = 3, p = 1 -> df = 26; check the critical value matches.
  expect_equal(attr(ci, "critical_value"),
               qbryant_paulson(0.95, 1, 3, 26), tolerance = 1e-8)
})

test_that("ci_c_ancova_bp refuses an s_ancova that is not a positive number in words", {
  # An NA stopped with "missing value where TRUE/FALSE needed".
  for (bad in list(NULL, 0, -1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(ci_c_ancova_bp(adj_means = c(1, 2, 3), s_ancova = bad,
                                n_per_group = 10),
                 "s_ancova", info = deparse(bad))
  }
  expect_error(ci_c_ancova_bp(adj_means = c(1, 2, 3), s_ancova = NA_real_,
                              n_per_group = 10),
               "'s_ancova' must be a single positive, finite number")
  expect_error(ci_c_ancova_bp(adj_means = c(1, 2, 3), n_per_group = 10),
               "standard deviation of the errors of the ANCOVA model")
})

test_that("the size of the Bryant-Paulson excess over Tukey stated on ?ci_c_ancova_bp", {
  # The help page says that over Bryant and Paulson's (1976) Table 1 the
  # Bryant-Paulson critical value exceeds Tukey's studentized range value by
  # less than 5% for every number of means and both tabled alpha levels
  # when nu is at least 14 with one covariate, 30 with two, or 40 with
  # three; that the excess grows as nu falls and as the numbers of
  # covariates and means rise; and that it reaches 34% at alpha = .05 with
  # one covariate and nu = 2, 84% with three covariates, and 92% at
  # alpha = .01.
  # Oracles, with no package code: the values printed in Table 1 of Bryant
  # and Paulson (1976), read from the scan to the two decimals printed; the
  # studentized range quantile from stats::qtukey() at nu >= 12, where
  # R's ptukey is accurate to about 1e-7; and, at nu = 2, where it is not
  # (it is off by up to 0.1% at alpha = .05 and by 17% at alpha = .01),
  # the quantile by direct quadrature of the range distribution against the
  # chi distribution, checked against the closed form for two means,
  # sqrt(2) * qt(1 - alpha / 2, nu). The ratio rises with k in every row of
  # the table, so k = 20 is the binding cell for each bound.
  tk <- function(conf, k, nu) qtukey(conf, nmeans = k, df = nu)
  # Below 1.05 at the thresholds, k = 20, both alpha levels.
  # Columns: conf, covariates, nu, printed value at k = 20.
  below <- rbind(c(.95, 1, 14, 6.29), c(.99, 1, 14, 7.75),
                 c(.95, 2, 30, 5.68), c(.99, 2, 30, 6.66),
                 c(.95, 3, 40, 5.57), c(.99, 3, 40, 6.47))
  for (i in seq_len(nrow(below))) {
    r <- below[i, 4] / tk(below[i, 1], 20, below[i, 3])
    expect_gt(r, 1)
    expect_lt(r, 1.05)
  }
  # Above 1.05 one tabled nu earlier for one and three covariates, and at
  # nu = 20 for two (at nu = 24 the k = 20 value, 1.0506, is too near 1.05
  # for two printed decimals to settle).
  above <- rbind(c(.99, 1, 12, 8.18), c(.99, 2, 20, 7.25), c(.99, 3, 30, 6.78))
  for (i in seq_len(nrow(above))) {
    expect_gt(above[i, 4] / tk(above[i, 1], 20, above[i, 3]), 1.05)
  }
  skip_on_cran()  # numeric quadrature for the nu = 2 quantiles and the package checks
  q_sr <- function(conf, k, nu) {
    range_upper <- function(w) 1 - integrate(function(z)
      k * dnorm(z) * (pnorm(z) - pnorm(z - w))^(k - 1), -Inf, Inf,
      rel.tol = 1e-12)$value
    log_g <- (nu / 2) * log(nu) - lgamma(nu / 2) - (nu / 2 - 1) * log(2)
    upper <- function(q) integrate(function(s)
      sapply(q * s, range_upper) * exp(log_g + (nu - 1) * log(s) - nu * s^2 / 2),
      0, Inf, rel.tol = 1e-10)$value
    uniroot(function(q) upper(q) - (1 - conf), c(1, 100), tol = 1e-10)$root
  }
  expect_equal(q_sr(0.95, 2, 2), sqrt(2) * qt(0.975, 2), tolerance = 1e-8)
  q20_05 <- q_sr(0.95, 20, 2)
  q20_01 <- q_sr(0.99, 20, 2)
  expect_equal(round(22.40 / q20_05 - 1, 2), 0.34)  # one covariate
  expect_equal(round(30.83 / q20_05 - 1, 2), 0.84)  # three covariates
  expect_equal(round(73.01 / q20_01 - 1, 2), 0.92)  # three covariates, alpha = .01
  # The package reproduces the printed values in the binding cells.
  for (i in seq_len(nrow(below))) {
    expect_equal(round(qbryant_paulson(below[i, 1], n_covariates = below[i, 2],
                                       n_groups = 20, df = below[i, 3]), 2),
                 below[i, 4])
  }
})

test_that("Bryant-Paulson intervals hold their familywise coverage where Tukey's do not (raw data)", {
  skip_on_cran()  # a seeded simulation of 4000 one-way ANCOVA data sets
  # Oracle: one-way balanced ANCOVA data with three random normal covariates,
  # k = 4 groups of n = 3 (nu = 12 - 4 - 3 = 5), equal adjusted means. The
  # adjusted means and the error standard deviation are computed here by
  # within-group regression, without package code. All pairwise intervals
  # cover their contrasts exactly when the range of the adjusted means is at
  # most q * s / sqrt(n), so the familywise coverage of the Bryant-Paulson
  # and of the Tukey intervals is the proportion of data sets with the
  # studentized range of the adjusted means below each critical value.
  k <- 4; n <- 3; p <- 3; R <- 4000
  nu <- k * n - k - p
  g <- rep(seq_len(k), each = n)
  stat <- withr::with_seed(20261111, vapply(seq_len(R), function(r) {
    X <- matrix(rnorm(k * n * p), k * n, p)
    y <- drop(X %*% rep(0.6, p)) + rnorm(k * n)
    xm <- rowsum(X, g) / n
    ym <- rowsum(y, g) / n
    Xw <- X - xm[g, , drop = FALSE]
    yw <- y - ym[g]
    b <- solve(crossprod(Xw), crossprod(Xw, yw))
    s <- sqrt(sum((yw - Xw %*% b)^2) / nu)
    adj <- drop(ym - (xm - matrix(colMeans(X), k, p, byrow = TRUE)) %*% b)
    diff(range(adj)) / (s / sqrt(n))
  }, numeric(1)))
  mcse <- sqrt(0.95 * 0.05 / R)
  cover_bp <- mean(stat <= qbryant_paulson(0.95, n_covariates = p,
                                           n_groups = k, df = nu))
  cover_tukey <- mean(stat <= qtukey(0.95, nmeans = k, df = nu))
  expect_lt(abs(cover_bp - 0.95), 4 * mcse)
  expect_lt(cover_tukey, 0.95 - 4 * mcse)
})

test_that("p = 0 at small df is the exact studentized range, not qtukey()", {
  skip_on_cran()  # a spline of the range distribution for twenty means
  # Oracle: an independent quadrature of the studentized range (the 1.1.0
  # review's F3 evidence), q_.99;20,2 = 37.943462287637; stats::qtukey()
  # gives 44.526266.
  expect_equal(qbryant_paulson(0.99, n_covariates = 0, n_groups = 20, df = 2),
               37.943462287637, tolerance = 1e-8)
  expect_equal(pbryant_paulson(37.943462287637, n_covariates = 0, n_groups = 20,
                               df = 2), 0.99, tolerance = 1e-8)
  # The reduction to the Tukey critical value holds there too.
  expect_equal(cv_bryant_paulson(.01, df = 2, n_groups = 20, n_covariates = 0,
                                 verbose = FALSE)$value,
               sqrt(2) * cv_tukey_hsd(.01, df = 2, n_groups = 20,
                                      verbose = FALSE)$value,
               tolerance = 1e-10)
})

test_that("with covariates the far-tail quantile integrates the exact studentized range", {
  skip_on_cran()
  # Through ptukey() at every node the quantile was 12.977072; an
  # independent Gauss-Legendre oracle of the Bryant-Paulson integral with
  # the exact studentized range gives 12.977142 (review evidence,
  # work/round6/evidence/g2_qtukey_extend/review2).
  q <- qbryant_paulson(.9999, n_covariates = 1, n_groups = 10, df = 11)
  expect_equal(q, 12.977142217, tolerance = 1e-8)
  expect_equal(pbryant_paulson(q, 1, 10, 11), .9999, tolerance = 1e-9)
})

test_that("with covariates the distribution reaches its infinite-df limit at very large df", {
  skip_on_cran()  # the exact studentized range at every node up to 1e8 df
  # From about 4e7 df DMAR 1.0.x integrated the shrinkage weight where it had
  # no mass: pbryant_paulson(3.7, 1, 4, 1e8) was 5.6e-25 (it is .956),
  # qbryant_paulson() stopped or, at 4.09e7 df, returned 4.53; df = Inf
  # stopped. Oracle, no package code: Q_p = R_k / V, V = S sqrt(delta), nu S^2
  # chi square on nu df, 1 - delta Beta(p/2, (nu + 1)/2), independent (see
  # R/bryant_paulson.R), so P(Q_p <= w) = E[H(w V)], H the CDF of R_k; to
  # second order in V - 1, with E[S] and E[sqrt(delta)] from the series of
  # Gamma(x + 1/2) / Gamma(x) and E[V^2] = a / (a + b). Remainder O(nu^-2).
  k <- 4
  H <- function(w, d = 0) {
    f <- switch(d + 1,
      function(z) dnorm(z) * (pnorm(z) - pnorm(z - w))^(k - 1),
      function(z) (k - 1) * dnorm(z) * (pnorm(z) - pnorm(z - w))^(k - 2) *
        dnorm(z - w),
      function(z) {
        D <- pnorm(z) - pnorm(z - w)
        g <- dnorm(z - w)
        (k - 1) * dnorm(z) * ((k - 2) * D^(k - 3) * g^2 + D^(k - 2) * (z - w) * g)
      })
    k * integrate(f, -Inf, Inf, rel.tol = 1e-13, abs.tol = 0)$value
  }
  lratio <- function(x) log1p(-1 / (8 * x) + 1 / (128 * x^2) + 5 / (1024 * x^3))
  oracle <- function(w, p, nu) {
    if (is.infinite(nu)) return(H(w))
    a <- (nu + 1) / 2
    b <- p / 2
    m1 <- expm1(lratio(nu / 2) - 0.5 * log1p(b / a) + lratio(a) - lratio(a + b))
    m2 <- -b / (a + b) - 2 * m1
    H(w) + w * H(w, 1) * m1 + w^2 * H(w, 2) * m2 / 2
  }
  oracle_q <- function(prob, p, nu, range) uniroot(
    function(w) oracle(w, p, nu) - prob, range, tol = 1e-13)$root
  for (nu in c(4.09e7, 5e7, 1e8, 1e10, Inf)) for (p in c(1, 3)) {
    expect_lt(abs(pbryant_paulson(3.7, p, k, nu) - oracle(3.7, p, nu)), 1e-9)
    expect_lt(abs(qbryant_paulson(.95, p, k, nu) - oracle_q(.95, p, nu, c(3.6, 3.7))),
              2e-8)
  }
  # The far tails, where every node uses the exact studentized range; at
  # 1e10 df the search for the root had no bracket.
  for (nu in c(1e8, 1e10)) {
    expect_lt(abs(qbryant_paulson(.9999, 1, k, nu) - oracle_q(.9999, 1, nu, c(6, 6.2))),
              2e-8)
  }
  # At 1e8 df the values are not yet the infinite-df ones.
  expect_gt(oracle(3.7, 1, Inf) - pbryant_paulson(3.7, 1, k, 1e8), 5e-9)
  # At infinite df the number of covariates no longer matters.
  expect_identical(pbryant_paulson(3.7, 2, k, Inf), pbryant_paulson(3.7, 0, k, Inf))
  expect_identical(qbryant_paulson(.95, 2, k, Inf), qbryant_paulson(.95, 0, k, Inf))
  # The critical value and the intervals built on the quantile follow.
  expect_equal(cv_bryant_paulson(.05, df = 1e9, n_groups = k, n_covariates = 1,
                                 verbose = FALSE)$value,
               oracle_q(.95, 1, 1e9, c(3.6, 3.7)), tolerance = 1e-8)
  ci <- ci_c_ancova_bp(adj_means = c(1, 2, 3, 4), s_ancova = 2, n_per_group = 10,
                       n_covariates = 1, df = 1e9)
  expect_equal(attr(ci, "critical_value"), oracle_q(.95, 1, 1e9, c(3.6, 3.7)),
               tolerance = 1e-8)
  # Many covariates, where the factors of the weight overflow and underflow
  # apart (1.0.x returned 0 with 150 covariates at 1e6 df). Oracle for two
  # means: 2 pt(w sqrt(delta) / sqrt 2, nu) - 1 averaged over the quantiles
  # of 1 - delta.
  oracle_2 <- function(w, p, nu) integrate(function(t)
    2 * pt(w * sqrt(1 - qbeta(t, p / 2, (nu + 1) / 2, lower.tail = FALSE)) /
             sqrt(2), nu) - 1, 0, 1, rel.tol = 1e-12)$value
  for (cell in list(c(150, 1e6), c(200, 1e5), c(40, 1e20)))
    expect_lt(abs(pbryant_paulson(3, cell[1], 2, cell[2]) -
                    oracle_2(3, cell[1], cell[2])), 1e-9)
  expect_lt(abs(qbryant_paulson(.95, 200, 2, 1e5) - uniroot(
    function(w) oracle_2(w, 200, 1e5) - .95, c(2.7, 2.9), tol = 1e-12)$root), 1e-8)
})
