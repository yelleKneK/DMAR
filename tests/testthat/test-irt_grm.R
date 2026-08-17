# Generate ordered responses directly through the normal ogive boundary
# response functions P*_ik(theta) = Phi(a_i (theta - b_ik)). Generating this
# way (rather than by cutting a latent response variate at thresholds) keeps
# the check on the lambda-to-a and tau-to-b conversions independent of the
# algebra those conversions use.
.grm_simulate <- function(n, a, b, seed) {
  set.seed(seed)
  theta <- stats::rnorm(n)
  responses <- vapply(seq_along(a), function(i) {
    p_star <- outer(theta, b[i, ], function(z, bk) stats::pnorm(a[i] * (z - bk)))
    as.integer(1 + rowSums(stats::runif(n) < p_star))
  }, integer(n))
  colnames(responses) <- paste0("item", seq_along(a))
  as.data.frame(responses)
}

# Key one five-category item against the rest of the scale, so a respondent
# high on theta answers low on it.
.grm_reverse_key <- function(d, j) {
  d[[j]] <- 6L - d[[j]]
  names(d)[j] <- "rev_item"
  d
}

.grm_a_pop <- c(1.2, 0.8, 1.6, 1.0, 1.4, 1.1)
.grm_b_pop <- rbind(c(-1.6, -0.6, 0.3, 1.2),
                    c(-1.4, -0.4, 0.5, 1.5),
                    c(-1.8, -0.7, 0.2, 1.1),
                    c(-1.2, -0.2, 0.7, 1.6),
                    c(-1.5, -0.5, 0.4, 1.3),
                    c(-1.3, -0.3, 0.6, 1.4))


test_that("irt_grm() recovers the generating graded response parameters", {
  skip_if_not_installed("lavaan")
  skip_on_cran()                       # n = 5000, six items, one WLSMV fit
  d <- .grm_simulate(5000, .grm_a_pop, .grm_b_pop, seed = 113)
  res <- irt_grm(d)

  expect_s3_class(res, "dmar_tbl")
  expect_equal(nrow(res), nrow(.grm_b_pop) * ncol(.grm_b_pop))

  a_hat <- res$a[!duplicated(res$item)]
  b_hat <- matrix(res$b, nrow = nrow(.grm_b_pop), byrow = TRUE)

  # The tolerances below are calibrated to this design rather than chosen as
  # round numbers. The test is deterministic (the responses come from
  # set.seed(113)), and the clean-build deviation at that seed is 0.035 on a
  # and 0.090 on b. Rerunning the same n = 5000 design over the twelve seeds
  # 113, 7, 19, 42, 101, 202, 303, 404, 505, 606, 707, and 808 puts the whole
  # clean spread in [0.031, 0.084] on a and [0.039, 0.115] on b, so 0.11 and
  # 0.13 sit above everything measured while cutting the earlier round 0.15
  # to something a conversion error has to fit inside.
  #
  # The limit of what this check can do was measured too, and it is worth
  # recording. Scaling the standardized loadings by 1.015, a 1.5 percent
  # error in the metric that feeds both a and b, moves the deviation on a
  # into [0.061, 0.170] across those same twelve seeds. That range overlaps
  # the clean range, so no threshold on this comparison separates a clean
  # build from a systematically biased one: at n = 5000 the sampling error
  # and the bias are the same size. Resolving an error that small is the job
  # of the just identified analytic anchor below, which compares against a
  # closed form instead of against a population value and so is not limited
  # by sampling error at all.
  expect_lt(max(abs(a_hat - .grm_a_pop)), 0.11)   # 3.2x the seed 113 deviation
  expect_lt(max(abs(b_hat - .grm_b_pop)), 0.13)   # 1.4x the seed 113 deviation
  expect_gt(stats::cor(a_hat, .grm_a_pop), 0.99)
  expect_gt(stats::cor(as.vector(b_hat), as.vector(.grm_b_pop)), 0.99)

  # The boundary locations of a positively keyed item are ordered by
  # construction, and every item here is positively keyed.
  expect_true(all(res$a > 0))
  for (nm in unique(res$item)) {
    bb <- res$b[res$item == nm]
    expect_equal(bb, sort(bb))
  }
})


test_that("irt_grm() matches the closed form of a just identified model", {
  skip_if_not_installed("lavaan")
  # With three items the single-factor model is just identified: three
  # polychoric correlations and three loadings, so df = 0 and the weighted
  # least squares solution reproduces the sample statistics exactly rather
  # than approximately. Both the loadings and the thresholds therefore have a
  # closed form, and it is computed here from quantities the conversion never
  # touches: the thresholds straight from the observed cumulative
  # proportions, the loadings from the polychoric correlations through
  # lambda_i lambda_j = R_ij. Checking the returned columns against that
  # closed form pins the conversion at optimizer precision, so a systematic
  # error in the standardized loading metric of even a fraction of a percent
  # fails here. The identity check further down cannot do this, because it
  # recomputes a and b from the returned lambda and so is satisfied by any
  # lambda the function chooses to report.
  d <- .grm_simulate(1500, .grm_a_pop[1:3], .grm_b_pop[1:3, , drop = FALSE],
                     seed = 113)
  res <- irt_grm(d)
  fit <- attr(res, "fit")

  # Just identified, hence an exactly reproduced sample polychoric matrix.
  expect_equal(lavaan::fitMeasures(fit)[["df"]], 0)
  expect_lt(lavaan::fitMeasures(fit)[["chisq"]], 1e-8)

  # lambda_1 lambda_2 = R_12, lambda_1 lambda_3 = R_13, lambda_2 lambda_3 =
  # R_23 solve to lambda_1 = sqrt(R_12 R_13 / R_23) and cyclically. The
  # positive root is the right one because the sign convention points the
  # factor in the direction the scale as a whole measures.
  R <- lavaan::lavInspect(fit, "sampstat")$cov
  lambda_closed <- c(sqrt(R[1, 2] * R[1, 3] / R[2, 3]),
                     sqrt(R[1, 2] * R[2, 3] / R[1, 3]),
                     sqrt(R[1, 3] * R[2, 3] / R[1, 2]))
  # A threshold is the normal quantile of an observed cumulative proportion.
  tau_closed <- unlist(lapply(d, function(x) {
    p <- cumsum(prop.table(table(x)))
    stats::qnorm(p[-length(p)])
  }), use.names = FALSE)

  n_boundaries <- as.integer(table(factor(res$item, levels = unique(res$item))))
  lambda_long <- rep(lambda_closed, times = n_boundaries)
  a_closed <- lambda_long / sqrt(1 - lambda_long^2)
  b_closed <- tau_closed / lambda_long

  expect_equal(res$lambda, lambda_long, tolerance = 1e-6)
  expect_equal(res$tau, tau_closed, tolerance = 1e-6)
  expect_equal(res$a, a_closed, tolerance = 1e-6)
  expect_equal(res$b, b_closed, tolerance = 1e-6)
})


test_that("irt_grm() agrees with mirt's marginal maximum likelihood fit", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # Exact agreement is not expected and is not asserted. lavaan fits the
  # model by limited information weighted least squares on the polychoric
  # correlations; mirt fits it by marginal maximum likelihood on the full
  # response patterns. The two are consistent for the same population
  # parameters but are different estimators, so they differ on a finite
  # sample. mirt reports slopes in the logistic metric, so the slopes
  # pinned below were divided by 1.702 to put them on the normal ogive
  # metric. The pinned values are tied to the seed 113 data above.
  d <- .grm_simulate(5000, .grm_a_pop, .grm_b_pop, seed = 113)
  res <- irt_grm(d)
  a_lavaan <- res$a[!duplicated(res$item)]
  b_lavaan <- matrix(res$b, nrow = nrow(.grm_b_pop), byrow = TRUE)

  # Pinned from mirt::mirt (mirt 1.46.1, 2026-08-09); live comparison in
  # tools/oracle_checks.R.
  a_mirt <- c(1.240697009039211, 0.8053765558454936, 1.661579227276265,
              1.00372256208225, 1.472715025230759, 1.093234183592852)
  b_mirt <- rbind(
    c(-1.643647417339596, -0.6197973522731184,
      0.2602490049858169, 1.193977102164867),
    c(-1.400208569169318, -0.4063540574952019,
      0.4296832018598173, 1.49013130747088),
    c(-1.902713372100115, -0.7556234880945267,
      0.1672152003788583, 1.078322538257959),
    c(-1.247186061207261, -0.2604712558518215,
      0.6458177246924814, 1.553404541368955),
    c(-1.501530328537965, -0.5388726992687777,
      0.3616992512529445, 1.238273429659636),
    c(-1.393367005853113, -0.3668837991809947,
      0.5618250892563006, 1.384803776851262)
  )

  expect_gt(stats::cor(a_lavaan, a_mirt), 0.99)
  expect_gt(stats::cor(as.vector(b_lavaan), as.vector(as.matrix(b_mirt))), 0.99)
  expect_lt(mean(abs(a_lavaan - a_mirt)), 0.10)
  expect_lt(mean(abs(b_lavaan - as.matrix(b_mirt))), 0.10)
})


test_that("irt_grm() satisfies the loading-to-discrimination identity", {
  skip_if_not_installed("lavaan")
  d <- .grm_simulate(1200, .grm_a_pop, .grm_b_pop, seed = 7)
  res <- irt_grm(d)

  # Recomputed from the returned columns, so the identity is checked on the
  # object the user receives, not on an intermediate. This is a consistency
  # check on the returned table and nothing more: it is satisfied by any
  # lambda the function reports, which is why the closed-form anchor above
  # exists to pin the metric of lambda itself.
  expect_equal(res$a, res$lambda / sqrt(1 - res$lambda^2), tolerance = 1e-12)
  expect_equal(res$b, res$tau / res$lambda, tolerance = 1e-12)

  # Every reported column is numeric except the two label columns, and the
  # metadata lives on attributes.
  expect_type(res$item, "character")
  expect_type(res$factor, "character")
  for (nm in c("category", "lambda", "tau", "a", "b")) {
    expect_true(is.numeric(res[[nm]]))
  }
  expect_identical(attr(res, "metric"), "normal_ogive")
  expect_identical(attr(res, "estimator"), "WLSMV")
  expect_identical(unname(attr(res, "n_categories")), rep(5L, 6L))
  expect_identical(attr(res, "N"), 1200L)
  expect_identical(attr(res, "factor_sign_flipped"), FALSE)
  expect_s4_class(attr(res, "fit"), "lavaan")
  expect_true(is.numeric(attr(res, "fit_measures")))
  expect_true(all(c("cfi", "rmsea") %in% names(attr(res, "fit_measures"))))
})


test_that("irt_grm() fixes the direction of the factor across column orders", {
  skip_if_not_installed("lavaan")
  # A single-factor model is identified only up to the direction of theta, so
  # lavaan returns whichever direction its starting values point toward. On a
  # scale containing a reverse-keyed item those starting values depend on the
  # order of the columns, and before the sign convention was imposed the same
  # six items fitted in two column orders came back with every discrimination
  # and every boundary location flipped, on an identical fit, with no
  # warning. The ordering b_i1 < b_i2 < ... then failed for the five items
  # that happened to disagree with the sign of the first column.
  d <- .grm_reverse_key(.grm_simulate(1200, .grm_a_pop, .grm_b_pop, seed = 7),
                        2L)
  d_moved <- d[, c(2, 1, 3, 4, 5, 6)]

  res <- suppressWarnings(irt_grm(d))
  res_moved <- suppressWarnings(irt_grm(d_moved))

  # Same model either way, so the fit must be identical.
  expect_equal(lavaan::fitMeasures(attr(res, "fit"))[["chisq"]],
               lavaan::fitMeasures(attr(res_moved, "fit"))[["chisq"]],
               tolerance = 1e-8)

  # And now so must every reported parameter. The tolerance here is optimizer
  # noise between two runs from different starting values, not sampling
  # error; the observed disagreement is on the order of 1e-7.
  m <- match(paste(res$item, res$category),
             paste(res_moved$item, res_moved$category))
  expect_false(anyNA(m))
  for (nm in c("lambda", "tau", "a", "b")) {
    expect_equal(res[[nm]], res_moved[[nm]][m], tolerance = 1e-5)
  }

  # This fixture is the one that used to flip: exactly one of the two column
  # orders needs the correction, and the attribute reports which.
  expect_type(attr(res, "factor_sign_flipped"), "logical")
  expect_length(attr(res, "factor_sign_flipped"), 1L)
  expect_true(xor(attr(res, "factor_sign_flipped"),
                  attr(res_moved, "factor_sign_flipped")))

  # The convention is that the factor runs in the direction the scale as a
  # whole measures, so the loadings sum to a positive number either way and
  # the five positively keyed items keep ascending boundary locations.
  for (r in list(res, res_moved)) {
    expect_gt(sum(r$lambda[!duplicated(r$item)]), 0)
    for (nm in setdiff(unique(r$item), "rev_item")) {
      bb <- r$b[r$item == nm]
      expect_equal(bb, sort(bb))
      expect_true(all(r$a[r$item == nm] > 0))
    }
  }
})


test_that("irt_grm() names an item keyed against the rest of the scale", {
  skip_if_not_installed("lavaan")
  d_ok <- .grm_simulate(1200, .grm_a_pop, .grm_b_pop, seed = 7)
  d_rev <- .grm_reverse_key(d_ok, 2L)

  # The sign indeterminacy is handled silently because it is a relabeling; a
  # genuinely reverse-keyed item is a property of the data and is reported.
  expect_warning(irt_grm(d_rev), "rev_item")
  expect_warning(irt_grm(d_rev), "Reverse score")
  res_rev <- suppressWarnings(irt_grm(d_rev))
  expect_true(all(res_rev$a[res_rev$item == "rev_item"] < 0))
  # Descending boundary locations, which is exactly what the warning says.
  b_rev <- res_rev$b[res_rev$item == "rev_item"]
  expect_equal(b_rev, sort(b_rev, decreasing = TRUE))

  # A scale with no reverse-keyed item needs no flip and says nothing.
  expect_no_warning(res_ok <- irt_grm(d_ok))
  expect_false(attr(res_ok, "factor_sign_flipped"))
  expect_true(all(res_ok$a > 0))

  # Reverse scoring the offending item is the documented remedy, and it is
  # what makes the result usable by irt_information(), which requires every
  # discrimination to be positive.
  d_fixed <- d_rev
  d_fixed$rev_item <- 6L - d_fixed$rev_item
  expect_no_warning(res_fixed <- irt_grm(d_fixed))
  expect_true(all(res_fixed$a > 0))
  for (nm in unique(res_fixed$item)) {
    bb <- res_fixed$b[res_fixed$item == nm]
    expect_equal(bb, sort(bb))
  }
  expect_error(irt_information(grm = res_rev), "greater than zero")
  expect_s3_class(irt_information(grm = res_fixed), "data.frame")
})


test_that("irt_grm() reports both discrimination metrics, 1.702 apart", {
  skip_if_not_installed("lavaan")
  d <- .grm_simulate(1200, .grm_a_pop, .grm_b_pop, seed = 7)
  normal <- irt_grm(d, metric = "normal_ogive")
  logistic <- irt_grm(d, metric = "logistic")

  expect_equal(logistic$a, 1.702 * normal$a, tolerance = 1e-12)
  # The scaling multiplies the slope and leaves the location alone.
  expect_equal(logistic$b, normal$b, tolerance = 1e-12)
  expect_equal(logistic$lambda, normal$lambda, tolerance = 1e-12)
  # The metric not reported in the a column is on an attribute, so both are
  # always available and the table keeps one shape.
  expect_equal(unname(attr(normal, "a_logistic")),
               1.702 * normal$a[!duplicated(normal$item)], tolerance = 1e-12)
  expect_equal(unname(attr(logistic, "a_normal_ogive")),
               normal$a[!duplicated(normal$item)], tolerance = 1e-12)
  expect_null(attr(normal, "a_normal_ogive"))
  expect_null(attr(logistic, "a_logistic"))
  expect_identical(dim(normal), dim(logistic))
})


test_that("irt_grm() reduces to the two parameter normal ogive on binary items", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  a_pop <- c(1.3, 0.9, 1.7, 1.1, 1.5)
  b_pop <- matrix(c(-0.8, -0.2, 0.3, 0.9, 0.0), ncol = 1L)
  d <- .grm_simulate(6000, a_pop, b_pop, seed = 113)
  expect_true(all(vapply(d, function(x) length(unique(x)), integer(1)) == 2L))

  res <- irt_grm(d)
  expect_equal(nrow(res), length(a_pop))       # one boundary per item
  expect_true(all(res$category == 1L))
  expect_true(all(is.finite(res$a)))
  expect_true(all(is.finite(res$b)))
  expect_true(all(res$a > 0))

  # Calibrated the same way as the graded design above. At seed 113 the
  # clean-build deviation is 0.054 on a and 0.062 on b. A binary item carries
  # one bit, so the spread on a across the same twelve seeds is wide, up to
  # 0.170 at seed 505, while b stays at 0.062. The thresholds below are set
  # for the seed this test actually uses, and changing that seed calls for
  # remeasuring them. At this seed the tightened bound on a does catch the
  # 1.015 loading-scale mutation described above, which lands at 0.169.
  expect_lt(max(abs(res$a - a_pop)), 0.15)     # 2.8x the seed 113 deviation
  expect_lt(max(abs(res$b - as.vector(b_pop))), 0.10)  # 1.6x that deviation

  # The item response function of a binary item is Phi(a (theta - b)), so
  # the model implied probability of the higher category at theta = b is
  # one half.
  expect_equal(stats::pnorm(res$a * (res$b - res$b)), rep(0.5, nrow(res)))
})


test_that("irt_grm() rejects continuous, too-few, and degenerate items", {
  skip_if_not_installed("lavaan")
  d <- .grm_simulate(200, .grm_a_pop, .grm_b_pop, seed = 7)

  # A continuous column handed to a graded response model.
  d_cont <- d
  set.seed(113)
  d_cont$item2 <- stats::rnorm(nrow(d_cont))
  expect_error(irt_grm(d_cont), "looks continuous")

  # Integer coded but with far too many distinct values to be categories.
  d_many <- d
  d_many$item3 <- seq_len(nrow(d_many))
  expect_error(irt_grm(d_many), "looks continuous")

  # Fewer than three items cannot identify a single-factor model.
  expect_error(irt_grm(d[, 1, drop = FALSE]), "At least 3 items")
  expect_error(irt_grm(d, items = c("item1", "item2")), "At least 3 items")

  # An item with a single observed category carries no information.
  d_const <- d
  d_const$item4 <- 3L
  expect_error(irt_grm(d_const), "at least 2 categories")

  # Argument validation.
  expect_error(irt_grm(d, estimator = "ML"), "'estimator' must be one of")
  expect_error(irt_grm(d, items = c("item1", "nope", "item3")),
               "not columns of 'data'")
  expect_error(irt_grm(d, metric = "probit"), "should be one of")
  expect_error(irt_grm("not data"), "must be a data.frame or a matrix")
  names(d)[1] <- "bad name"
  expect_error(irt_grm(d), "syntactically valid R names")
})
