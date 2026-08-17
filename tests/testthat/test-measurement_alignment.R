# Simulate a multiple group single factor data set with known group factor
# means and standard deviations. `Lambda` and `Nu` are G x I matrices of the
# generating loadings and intercepts, so planting noninvariance is a matter
# of changing one cell.
.alignment_sim <- function(n_g, factor_mean, factor_sd, Lambda, Nu,
                           resid_sd = 0.6) {
  G <- length(n_g)
  I <- ncol(Lambda)
  do.call(rbind, lapply(seq_len(G), function(g) {
    eta <- stats::rnorm(n_g[g], factor_mean[g], factor_sd[g])
    x <- vapply(seq_len(I), function(i) {
      Nu[g, i] + Lambda[g, i] * eta + stats::rnorm(n_g[g], 0, resid_sd)
    }, numeric(n_g[g]))
    out <- as.data.frame(x)
    names(out) <- paste0("x", seq_len(I))
    out$cohort <- paste0("cohort_", g)
    out
  }))
}

# Data whose maximum likelihood covariance matrix (divisor N, the convention
# lavaan's normal likelihood uses) and mean vector equal the supplied
# population values exactly. When Sigma is exactly a one factor structure the
# configural solution then equals the population solution to optimizer
# precision, which turns the exact invariance case into an analytic anchor.
.alignment_exact <- function(n, mu, Sigma) {
  p <- length(mu)
  z <- scale(matrix(stats::rnorm(n * p), n, p), center = TRUE, scale = FALSE)
  z <- z %*% solve(chol(crossprod(z) / n)) %*% chol(Sigma)
  sweep(z, 2L, mu, "+")
}

test_that("measurement_alignment() agrees with sirt's invariance alignment", {
  skip_if_not_installed("lavaan")
  set.seed(20)
  Lambda <- matrix(0.8, nrow = 4, ncol = 6)
  Nu <- matrix(1.0, nrow = 4, ncol = 6)
  Lambda[2, 3] <- 0.3
  Nu[3, 5] <- 1.8
  Nu[4, 1] <- 0.4
  d <- .alignment_sim(c(200, 250, 180, 300), c(0, 0.3, -0.5, 0.8),
                      c(1, 1.2, 0.8, 1.1), Lambda, Nu)
  out <- measurement_alignment(d, items = paste0("x", 1:6),
                               group = "cohort", seed = 113)

  # sirt takes the configural loading and intercept matrices directly, so it
  # was fed exactly the ones DMAR aligned when the reference solution below
  # was computed. Its weights are normalized by their column sum and it
  # counts each group pair in both orders, so its criterion is DMAR's
  # divided by sum(sqrt(N_g))^2 / 2; align.pow = 0.5 is its current coding
  # of the fourth root loss of Asparouhov and Muthen. The pinned solution is
  # tied to the seed 20 data above and to the configural fit on it.
  L0 <- attr(out, "configural_loadings")
  N0 <- attr(out, "configural_intercepts")
  n_g <- out$n

  # Pinned from sirt::invariance.alignment (sirt 4.2.133, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  oracle_alpha0 <- c(0, 0.3212615479342715,
                     -0.5096655006535324, 0.7510122420205143)
  oracle_psi0 <- c(1, 1.270279738625676,
                   0.8249961967097325, 1.204505111842333)
  oracle_lambda_aligned <- rbind(
    c(0.7356074076277439, 0.8319186422757626, 0.7534652660756466,
      0.7962795592449338, 0.8020121752464777, 0.7556140049305852),
    c(0.7701837934505398, 0.7742611528195155, 0.2952081097520632,
      0.7716890354033354, 0.7798405054949993, 0.7841414684011249),
    c(0.7882693604563553, 0.8836430952394124, 0.7511412093859433,
      0.717001664958853, 0.6613707669597664, 0.7432010724665805),
    c(0.7504886599221401, 0.775725710096056, 0.767718563001213,
      0.7960344406695111, 0.8116065583558054, 0.7044932369016907)
  )
  oracle_nu_aligned <- rbind(
    c(1.007523845355167, 1.009291699334871, 1.044792584335339,
      1.016263912833932, 1.008401161557097, 0.9890085035845576),
    c(0.9777226844269117, 1.021940353254991, 1.001117379170691,
      0.9910834139464593, 1.084300681053957, 1.017297594680787),
    c(1.016726897273215, 1.040704858281373, 0.9697977685963808,
      0.9829786068045168, 1.739240362580484, 0.9769897532132121),
    c(0.4625915992242624, 1.01041605881392, 0.9937780606175273,
      1.015919505568633, 1.0557243755891, 1.093186617434623)
  )
  oracle_fopt <- 3.471791690212806
  oracle_R2_total <- c(0.9833842071339252, 0.9800452721430867)

  expect_lt(max(abs(out$factor_mean - oracle_alpha0)), 5e-4)
  expect_lt(max(abs(sqrt(out$factor_variance) - oracle_psi0)), 5e-4)
  expect_lt(max(abs(attr(out, "aligned_loadings") - oracle_lambda_aligned)),
            5e-4)
  expect_lt(max(abs(attr(out, "aligned_intercepts") - oracle_nu_aligned)),
            5e-4)

  # The criterion itself, once sirt's weight normalization is undone. The
  # sharper version evaluates DMAR's criterion at sirt's own solution: the
  # two implementations compute the same function, so the agreement there is
  # to machine precision rather than to optimizer precision. DMAR's minimum
  # must in turn be no larger than the value at sirt's solution.
  rescaled <- oracle_fopt * sum(sqrt(n_g))^2 / 2
  machinery <- .alignment_machinery(L0, N0, n_g, 0.01, "fixed")
  at_oracle <- machinery$fn(c(oracle_alpha0[-1],
                              log(oracle_psi0[-1])))
  expect_equal(at_oracle, rescaled, tolerance = 1e-12)
  expect_equal(attr(out, "simplicity_function"), rescaled, tolerance = 1e-5)
  expect_lte(attr(out, "simplicity_function"), at_oracle + 1e-8)

  # The per item R2 measures pool to the overall effect sizes of approximate
  # invariance that sirt reports.
  expect_equal(unname(attr(out, "R2_total")), oracle_R2_total,
               tolerance = 1e-5)
})

test_that("measurement_alignment() recovers known factor means and variances", {
  skip_if_not_installed("lavaan")
  # The reference group's factor mean is 0 and its factor standard
  # deviation is 1, which is the metric "fixed" alignment reports on, so
  # the estimates compare directly with the generating values.
  set.seed(113)
  factor_mean <- c(0, 0.30, -0.50, 0.80, 0.20)
  factor_sd <- c(1, 1.20, 0.80, 1.10, 0.90)
  Lambda <- matrix(0.8, nrow = 5, ncol = 6)
  Nu <- matrix(1.0, nrow = 5, ncol = 6)
  Lambda[2, 3] <- 0.3                       # planted loading noninvariance
  Nu[4, 5] <- 1.8                           # planted intercept noninvariance
  d <- .alignment_sim(c(400, 450, 380, 500, 420), factor_mean, factor_sd,
                      Lambda, Nu)
  out <- measurement_alignment(d, items = paste0("x", 1:6),
                               group = "cohort", seed = 113)

  expect_s3_class(out, "dmar_tbl")
  expect_equal(nrow(out), 5L)
  expect_true(attr(out, "converged"))
  expect_lt(max(abs(out$factor_mean - factor_mean)), 0.08)
  expect_lt(max(abs(sqrt(out$factor_variance) - factor_sd)), 0.08)

  # The planted noninvariance is where the criterion says it is.
  R2_lambda <- attr(out, "R2_loadings")
  R2_nu <- attr(out, "R2_intercepts")
  expect_equal(names(which.min(R2_lambda)), "x3")
  expect_equal(names(which.min(R2_nu)), "x5")
  loss <- attr(out, "item_loss")
  expect_equal(rownames(loss)[which.max(loss[, "loadings"])], "x3")
  expect_equal(rownames(loss)[which.max(loss[, "intercepts"])], "x5")
  expect_true(all(R2_lambda[c("x1", "x2", "x4", "x5", "x6")] > 0.99))
  expect_true(all(R2_nu[c("x1", "x2", "x3", "x4", "x6")] > 0.99))

  # The same design at the population moments, which removes sampling error
  # and shows that planted noninvariance does not bias what alignment
  # reports: the estimates return the generating values, not values pulled
  # toward the two noninvariant parameters.
  resid_var <- rep(0.36, 6)
  exact <- do.call(rbind, lapply(1:5, function(g) {
    Sigma <- factor_sd[g]^2 * outer(Lambda[g, ], Lambda[g, ]) + diag(resid_var)
    mu <- Nu[g, ] + Lambda[g, ] * factor_mean[g]
    z <- as.data.frame(.alignment_exact(400, mu, Sigma))
    names(z) <- paste0("x", 1:6)
    z$cohort <- paste0("cohort_", g)
    z
  }))
  pop <- measurement_alignment(exact, items = paste0("x", 1:6),
                               group = "cohort", seed = 113)
  expect_lt(max(abs(pop$factor_mean - factor_mean)), 0.02)
  expect_lt(max(abs(pop$factor_variance - factor_sd^2)), 0.03)
  expect_equal(names(which.min(attr(pop, "R2_loadings"))), "x3")
  expect_equal(names(which.min(attr(pop, "R2_intercepts"))), "x5")
})

test_that("measurement_alignment() is exact under exact invariance", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  loading <- c(0.80, 0.70, 0.90, 0.60, 0.75)
  intercept <- c(1.00, 0.50, 1.20, 0.80, 0.30)
  resid_var <- rep(0.36, 5)
  # Neither the first group's factor standard deviation nor the geometric
  # mean of the standard deviations is 1, so the two identification rules
  # put the answer on two different metrics and the test can pin both.
  factor_mean <- c(0.5, 0.9, 0.2, 1.4, -0.1)
  factor_sd <- c(1.25, 1.3, 0.7, 1.1, 0.9)
  n_g <- c(300, 350, 400, 320, 290)
  d <- do.call(rbind, lapply(seq_along(n_g), function(g) {
    Sigma <- factor_sd[g]^2 * outer(loading, loading) + diag(resid_var)
    mu <- intercept + loading * factor_mean[g]
    out <- as.data.frame(.alignment_exact(n_g[g], mu, Sigma))
    names(out) <- paste0("x", 1:5)
    out$cohort <- paste0("cohort_", g)
    out
  }))
  out <- measurement_alignment(d, items = paste0("x", 1:5),
                               group = "cohort", seed = 113)

  # Fixing the reference group's factor variance at 1 puts the solution on
  # the reference group's latent metric, so the factor means come back in
  # units of its factor standard deviation and the variances as ratios.
  expect_lt(max(abs(out$factor_mean -
                      (factor_mean - factor_mean[1]) / factor_sd[1])), 1e-5)
  expect_lt(max(abs(out$factor_variance -
                      (factor_sd / factor_sd[1])^2)), 1e-5)
  expect_lt(max(abs(sweep(attr(out, "aligned_loadings"), 2L,
                          loading * factor_sd[1]))), 1e-5)
  expect_lt(max(abs(sweep(attr(out, "aligned_intercepts"), 2L,
                          intercept + loading * factor_mean[1]))), 1e-5)

  # With every difference exactly zero each of the 2 * I terms in a group
  # pair contributes epsilon^(1/4), so the criterion sits exactly on its
  # analytic floor.
  pairs <- utils::combn(length(n_g), 2L)
  floor_value <- 2 * 5 * 0.01^0.25 *
    sum(sqrt(n_g[pairs[1L, ]] * n_g[pairs[2L, ]]))
  expect_equal(attr(out, "simplicity_function"), floor_value,
               tolerance = 1e-8)
  expect_equal(unname(attr(out, "R2_loadings")), rep(1, 5), tolerance = 1e-8)
  expect_equal(unname(attr(out, "R2_intercepts")), rep(1, 5), tolerance = 1e-8)
  expect_equal(unname(attr(out, "R2_total")), c(1, 1), tolerance = 1e-8)

  # The free rule puts the same solution on the metric whose group factor
  # standard deviations have a geometric mean of 1, and centers the means.
  free <- measurement_alignment(d, items = paste0("x", 1:5),
                                group = "cohort", alignment = "free",
                                seed = 113)
  gm <- exp(mean(log(factor_sd)))
  expect_lt(max(abs(free$factor_mean -
                      (factor_mean - mean(factor_mean)) / gm)), 1e-5)
  expect_lt(max(abs(sqrt(free$factor_variance) - factor_sd / gm)), 1e-5)
  expect_equal(attr(free, "simplicity_function"), floor_value,
               tolerance = 1e-8)
})

test_that("measurement_alignment() discards the degenerate branch", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # On real data a real share of random starts slides down the direction in
  # which every factor variance but the reference group's runs off and the
  # aligned loadings collapse. Those runs stop with the optimizer reporting
  # convergence, so they have to be recognized and thrown away.
  data(holzinger_swineford, package = "DMAR", envir = environment())
  hs <- holzinger_swineford
  hs$school_sex <- interaction(hs$school, hs$sex, sep = ", ")
  verbal <- c("t5_general_information", "t6_paragraph_comprehension",
              "t7_sentence", "t8_word_classification", "t9_word_meaning")
  hs[verbal] <- scale(hs[verbal])
  out <- measurement_alignment(hs, items = verbal, group = "school_sex",
                               n_starts = 40, seed = 7)
  starts <- attr(out, "simplicity_starts")
  expect_length(starts, 40L)
  expect_true(any(!is.finite(starts)))           # some starts went degenerate
  expect_true(all(is.finite(out$factor_variance)))
  expect_true(all(out$factor_variance > 1e-4 & out$factor_variance < 1e4))
  expect_equal(attr(out, "simplicity_function"), min(starts))

  # Eight groups whose factor means are far apart. From the null start alone
  # the optimizer goes degenerate, which the function refuses to report; the
  # second start finds the real minimum, and it is the one 20 starts find.
  set.seed(1)
  G <- 8L
  n_g <- sample(120:300, G)
  Lambda <- matrix(0.8, nrow = G, ncol = 6)
  Nu <- matrix(1.0, nrow = G, ncol = 6)
  Lambda[2, 3] <- 0.2
  Lambda[5, 1] <- 1.5
  Nu[4, 5] <- 2.5
  Nu[7, 2] <- -0.5
  d <- .alignment_sim(n_g, seq(-2.5, 2.5, length.out = G),
                      exp(seq(-0.5, 0.5, length.out = G)), Lambda, Nu)
  expect_error(measurement_alignment(d, items = paste0("x", 1:6),
                                     group = "cohort", n_starts = 1L),
               "usable alignment")
  two <- measurement_alignment(d, items = paste0("x", 1:6), group = "cohort",
                               n_starts = 2L)
  many <- measurement_alignment(d, items = paste0("x", 1:6), group = "cohort",
                                n_starts = 20L, seed = 113)
  expect_identical(attr(two, "simplicity_starts")[1], Inf)   # the null start
  expect_equal(attr(two, "simplicity_function"),
               attr(many, "simplicity_function"), tolerance = 1e-6)
  expect_true(all(is.finite(many$factor_variance)))
  expect_true(all(many$factor_variance > 0))
})

test_that("measurement_alignment() searches past genuine local minima", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # Six groups with loadings and intercepts noninvariant everywhere, which
  # is when the criterion really is multimodal. Several distinct minima are
  # reached; the default number of starts still lands on the lowest.
  set.seed(24)
  G <- 6L
  n_g <- sample(100:250, G)
  Lambda <- matrix(runif(G * 6, 0.3, 1.1), nrow = G)
  Nu <- matrix(runif(G * 6, -1, 2), nrow = G)
  d <- .alignment_sim(n_g, runif(G, -1, 1), runif(G, 0.7, 1.4), Lambda, Nu)
  many <- measurement_alignment(d, items = paste0("x", 1:6), group = "cohort",
                                n_starts = 40L, seed = 3)
  starts <- attr(many, "simplicity_starts")
  finite_starts <- starts[is.finite(starts)]
  expect_gte(attr(many, "n_optima"), 2L)
  expect_gt(max(finite_starts), min(finite_starts))
  expect_equal(attr(many, "simplicity_function"), min(starts))
  default <- measurement_alignment(d, items = paste0("x", 1:6),
                                   group = "cohort", seed = 113)
  expect_equal(attr(default, "simplicity_function"),
               attr(many, "simplicity_function"), tolerance = 1e-6)
})

test_that("measurement_alignment() aligns Holzinger and Swineford's schools", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford, package = "DMAR", envir = environment())
  hs <- holzinger_swineford
  hs$school_sex <- interaction(hs$school, hs$sex, sep = ", ")
  verbal <- c("t5_general_information", "t6_paragraph_comprehension",
              "t7_sentence", "t8_word_classification", "t9_word_meaning")
  hs[verbal] <- scale(hs[verbal])
  out <- measurement_alignment(hs, items = verbal, group = "school_sex",
                               seed = 113)

  expect_equal(nrow(out), 4L)
  expect_equal(sum(out$n), sum(!is.na(rowSums(hs[verbal]))))
  expect_equal(out$factor_mean[1], 0)          # "fixed" reference group
  expect_equal(out$factor_variance[1], 1)
  expect_true(all(is.finite(out$factor_mean)))
  expect_true(all(out$factor_variance > 0))
  expect_true(attr(out, "converged"))
  expect_equal(dim(attr(out, "aligned_loadings")), c(4L, 5L))
  expect_equal(colnames(attr(out, "aligned_intercepts")), verbal)
  expect_true(inherits(attr(out, "fit"), "lavaan"))

  # Both Grant-White groups score above both Pasteur groups on the aligned
  # verbal factor, which is the comparison the exact invariance ladder
  # cannot license here.
  expect_true(all(out$factor_mean[3:4] > out$factor_mean[1:2]))

  # The aligned parameters are what the definition says they are.
  L0 <- attr(out, "configural_loadings")
  N0 <- attr(out, "configural_intercepts")
  sd_g <- sqrt(out$factor_variance)
  expect_equal(attr(out, "aligned_loadings"), L0 / sd_g)
  expect_equal(attr(out, "aligned_intercepts"),
               N0 - out$factor_mean * (L0 / sd_g))
})

test_that("the alignment gradient matches a numerical derivative", {
  set.seed(5)
  L0 <- matrix(runif(20, 0.4, 1.2), nrow = 4)
  N0 <- matrix(runif(20, -1, 2), nrow = 4)
  for (rule in c("fixed", "free")) {
    m <- .alignment_machinery(L0, N0, c(120, 200, 90, 160), 0.01, rule)
    par <- c(0.4, -0.3, 0.5, 0.2, -0.1, 0.35)
    numeric_gradient <- vapply(seq_along(par), function(j) {
      up <- par; up[j] <- up[j] + 1e-6
      dn <- par; dn[j] <- dn[j] - 1e-6
      (m$fn(up) - m$fn(dn)) / 2e-6
    }, numeric(1L))
    expect_equal(m$gr(par), numeric_gradient, tolerance = 1e-5)
  }
  # The criterion is the definition, computed the slow way.
  m <- .alignment_machinery(L0, N0, c(120, 200, 90, 160), 0.01, "fixed")
  par <- c(0.4, -0.3, 0.5, 0.2, -0.1, 0.35)
  p <- m$build(par)
  a <- m$aligned(p)
  n_g <- c(120, 200, 90, 160)
  by_hand <- 0
  for (g1 in 1:3) for (g2 in (g1 + 1):4) {
    w <- sqrt(n_g[g1] * n_g[g2])
    by_hand <- by_hand + w * sum(
      ((a$lambda[g1, ] - a$lambda[g2, ])^2 + 0.01)^0.25 +
        ((a$nu[g1, ] - a$nu[g2, ])^2 + 0.01)^0.25)
  }
  expect_equal(m$fn(par), by_hand)
  expect_equal(.alignment_n_optima(c(1, 1, 1 + 1e-12, 2, 2, 5)), 3L)
})

test_that("the free identification rule satisfies its own constraints", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  Lambda <- matrix(0.8, nrow = 4, ncol = 5)
  Nu <- matrix(1.0, nrow = 4, ncol = 5)
  Lambda[3, 2] <- 0.35
  d <- .alignment_sim(c(220, 260, 200, 240), c(0, 0.4, -0.3, 0.6),
                      c(1, 1.15, 0.85, 1.05), Lambda, Nu)
  free <- measurement_alignment(d, items = paste0("x", 1:5),
                                group = "cohort", alignment = "free",
                                seed = 113)
  fixed <- measurement_alignment(d, items = paste0("x", 1:5),
                                 group = "cohort", alignment = "fixed",
                                 seed = 113)
  expect_equal(mean(free$factor_mean), 0)
  expect_equal(prod(sqrt(free$factor_variance)), 1)
  expect_equal(attr(free, "alignment"), "free")
  expect_equal(fixed$factor_mean[1], 0)
  expect_equal(fixed$factor_variance[1], 1)

  # The two rules index different slices of the same family, so the factor
  # mean differences they report agree closely without being identical.
  expect_lt(max(abs(diff(free$factor_mean) - diff(fixed$factor_mean))), 0.05)
})

test_that("measurement_alignment() honors seed and leaves the RNG alone", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  Lambda <- matrix(0.8, nrow = 3, ncol = 4)
  Nu <- matrix(1.0, nrow = 3, ncol = 4)
  d <- .alignment_sim(c(150, 180, 160), c(0, 0.4, -0.2), c(1, 1.2, 0.9),
                      Lambda, Nu)
  set.seed(999)
  before <- .Random.seed
  a <- measurement_alignment(d, items = paste0("x", 1:4), group = "cohort",
                             seed = 42)
  expect_identical(.Random.seed, before)
  b <- measurement_alignment(d, items = paste0("x", 1:4), group = "cohort",
                             seed = 42)
  expect_equal(a$factor_mean, b$factor_mean)
  expect_equal(attr(a, "simplicity_starts"), attr(b, "simplicity_starts"))
  expect_equal(attr(a, "n_starts"), 10L)
  expect_equal(attr(a, "epsilon"), 0.01)
})

test_that("measurement_alignment() fails clearly on unusable input", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  Lambda <- matrix(0.8, nrow = 3, ncol = 5)
  Nu <- matrix(1.0, nrow = 3, ncol = 5)
  d <- .alignment_sim(c(150, 180, 160), c(0, 0.4, -0.2), c(1, 1.2, 0.9),
                      Lambda, Nu)
  items <- paste0("x", 1:5)

  expect_error(measurement_alignment(as.matrix(d), items, "cohort"),
               "must be a data.frame")
  expect_error(measurement_alignment(d, items[1:2], "cohort"),
               "three or more columns")
  expect_error(measurement_alignment(d, c(items, "not_a_column"), "cohort"),
               "three or more columns")
  expect_error(measurement_alignment(d, c(items[1], items), "cohort"),
               "must not repeat")
  expect_error(measurement_alignment(d, c(items[1:4], "cohort"), "cohort"),
               "must be numeric")
  expect_error(measurement_alignment(d, items, "nope"),
               "must name one column")
  expect_error(measurement_alignment(d[d$cohort == "cohort_1", ], items,
                                     "cohort"),
               "at least two groups")
  expect_error(measurement_alignment(d, items, "cohort", n_starts = 0),
               "positive integer")
  expect_error(measurement_alignment(d, items, "cohort", n_starts = 2.5),
               "positive integer")
  expect_error(measurement_alignment(d, items, "cohort", epsilon = 0),
               "positive number")
  expect_error(measurement_alignment(d, items, "cohort", seed = "abc"),
               "single number")
  expect_error(measurement_alignment(d, items, "cohort", model = 3),
               "lavaan model syntax")
  expect_error(measurement_alignment(d, items, "cohort", alignment = "loose"),
               "should be one of")
  expect_error(
    suppressWarnings(
      measurement_alignment(d, items, "cohort",
                            model = "f1 =~ x1 + x2 + x3\nf2 =~ x3 + x4 + x5")),
    "exactly one latent variable")
  expect_error(
    measurement_alignment(d, items[1:4], "cohort",
                          model = paste("f =~", paste(items, collapse = " + "))),
    "exactly the columns named in 'items'")
  expect_error(
    suppressWarnings(measurement_alignment(
      d, items, "cohort",
      model = paste("f =~", paste(items, collapse = " + "), "\nf ~~ 4*f"))),
    "standardize the factor")

  # A model that cannot be fitted at all fails loudly rather than returning
  # NA factor means.
  narrow <- d[d$cohort %in% c("cohort_1", "cohort_2"), ]
  narrow[narrow$cohort == "cohort_2", items] <- NA
  expect_error(suppressWarnings(
    measurement_alignment(narrow, items, "cohort")))
})

test_that("measurement_alignment() refuses a model that sets the scale through the loadings", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  Lambda <- matrix(0.8, nrow = 3, ncol = 4)
  Nu <- matrix(1.0, nrow = 3, ncol = 4)
  dat <- .alignment_sim(c(150, 150, 150), c(0, 0.3, -0.2), c(1, 1.2, 0.9),
                        Lambda, Nu)
  items <- paste0("x", 1:4)

  # A marker constraint pins x1's loading to 1 in every group while leaving
  # the factor variance at 1 and the mean at 0, so the psi/alpha check
  # cannot see it; the parameter-table guard must.
  expect_error(
    measurement_alignment(dat, items = items, group = "cohort",
                          model = "f =~ 1*x1 + x2 + x3 + x4"),
    "fixes or constrains a loading"
  )

  # An equality label ties x1's loading across groups; same guard.
  expect_error(
    measurement_alignment(dat, items = items, group = "cohort",
                          model = "f =~ c(a1, a1, a1)*x1 + x2 + x3 + x4"),
    "fixes or constrains a loading"
  )

  # The harmless explicit standardization is still accepted and returns
  # the same aligned solution as the default model (same seed, so the
  # multi-start search is identical; attributes such as the per-start
  # criterion values are start-order noise and not compared).
  base <- measurement_alignment(dat, items = items, group = "cohort",
                                seed = 113)
  std <- measurement_alignment(dat, items = items, group = "cohort",
                               model = "f =~ x1 + x2 + x3 + x4\nf ~~ 1*f",
                               seed = 113)
  expect_equal(as.data.frame(std), as.data.frame(base),
               tolerance = 1e-6, ignore_attr = TRUE)
})
