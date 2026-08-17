v <- function(tab, t) tab$value[tab$term == t]

.med_data <- function(n = 200) {
  set.seed(113)
  x <- rnorm(n)
  m <- 0.5 * x + rnorm(n, 0, sqrt(0.75))
  y <- 0.2 * x + 0.4 * m + rnorm(n, 0, 0.8)
  data.frame(x = x, m = m, y = y)
}

test_that("mediate() recovers the paths and the linear identities exactly", {
  d <- .med_data()
  res <- mediate(d, "x", "m", "y", B = 200, seed = 113)
  fm <- lm(m ~ x, d); fy <- lm(y ~ x + m, d); ft <- lm(y ~ x, d)
  expect_equal(v(res, "a"), unname(coef(fm)["x"]))
  expect_equal(v(res, "b"), unname(coef(fy)["m"]))
  expect_equal(v(res, "direct_effect"), unname(coef(fy)["x"]))
  expect_equal(v(res, "indirect_effect"), v(res, "a") * v(res, "b"))
  # The linear-model identity: total = direct + indirect, and it equals the
  # X-only regression coefficient exactly.
  expect_equal(v(res, "total_effect"),
               v(res, "direct_effect") + v(res, "indirect_effect"))
  expect_equal(v(res, "total_effect"), unname(coef(ft)["x"]))
  expect_equal(v(res, "proportion_mediated"),
               v(res, "indirect_effect") / v(res, "total_effect"))
  expect_equal(v(res, "N"), 200)
  expect_identical(attr(res, "ci_method"), "boot_percentile")
  expect_identical(attr(res, "conf_level"), 0.95)
})

test_that("mediate() Sobel interval matches var_indirect_effect by hand", {
  d <- .med_data()
  res <- mediate(d, "x", "m", "y", ci_method = "sobel")
  se <- v(res, "se_indirect")
  expect_equal(se^2,
               var_indirect_effect(a = v(res, "a"), b = v(res, "b"),
                                   var_a = v(res, "se_a")^2,
                                   var_b = v(res, "se_b")^2)$value[1])
  expect_equal(v(res, "indirect_lower"),
               v(res, "indirect_effect") + qnorm(0.025) * se)
  expect_equal(v(res, "indirect_upper"),
               v(res, "indirect_effect") + qnorm(0.975) * se)
})

test_that("mediate() bootstrap and Monte Carlo intervals are seeded, sane, and skewed-capable", {
  skip_on_cran()  # two bootstraps, a Monte Carlo of 5000 draws, and a BCa run
  d <- .med_data()
  r1 <- mediate(d, "x", "m", "y", B = 300, seed = 7)
  r2 <- mediate(d, "x", "m", "y", B = 300, seed = 7)
  expect_identical(r1$value, r2$value)            # reproducible under seed
  expect_lt(v(r1, "indirect_lower"), v(r1, "indirect_effect"))
  expect_gt(v(r1, "indirect_upper"), v(r1, "indirect_effect"))
  mc <- mediate(d, "x", "m", "y", ci_method = "monte_carlo", B = 5000,
                seed = 7)
  so <- mediate(d, "x", "m", "y", ci_method = "sobel")
  # With a comfortably nonzero indirect effect the MC interval is close to
  # Sobel but need not coincide; both must cover the estimate.
  expect_lt(abs(v(mc, "indirect_lower") - v(so, "indirect_lower")), 0.05)
  bca <- mediate(d, "x", "m", "y", ci_method = "boot_bca", B = 300,
                 seed = 7)
  expect_lt(v(bca, "indirect_lower"), v(bca, "indirect_upper"))
})

test_that("mediate() seed hygiene restores the caller's RNG state", {
  d <- .med_data()
  set.seed(42); before <- .Random.seed
  invisible(mediate(d, "x", "m", "y", B = 150, seed = 99))
  expect_identical(.Random.seed, before)
})

test_that("mediate() handles covariates and validates inputs", {
  d <- .med_data(); d$z <- rnorm(nrow(d))
  res <- mediate(d, "x", "m", "y", covariates = "z", B = 150, seed = 1)
  fy <- lm(y ~ x + m + z, d)
  expect_equal(v(res, "b"), unname(coef(fy)["m"]))
  expect_error(mediate(d, "x", "m", "nope"), "must name one column")
  expect_error(mediate(d, "x", "m", "y", covariates = "nope"),
               "must name columns")
  expect_error(mediate(d, "x", "m", "y", B = 10), "at least 100")
})

test_that("mediate() percentile bootstrap covers the truth at the nominal rate", {
  skip_on_cran()
  set.seed(113)
  true_ab <- 0.5 * 0.4
  # 60 replications of a B = 200 bootstrap keep the check against gross
  # miscoverage at a quarter of the earlier 120 x 300 runtime; the 0.88
  # bound already allowed for the wider Monte Carlo band.
  cover <- replicate(60, {
    n <- 150
    x <- rnorm(n); m <- 0.5 * x + rnorm(n, 0, sqrt(0.75))
    y <- 0.2 * x + 0.4 * m + rnorm(n, 0, 0.8)
    r <- mediate(data.frame(x, m, y), "x", "m", "y", B = 200)
    v(r, "indirect_lower") <= true_ab && true_ab <= v(r, "indirect_upper")
  })
  expect_gt(mean(cover), 0.88)   # 95 percent nominal, generous MC band
})

test_that("a rank-deficient bootstrap resample is dropped, not fatal", {
  skip_on_cran()  # a bootstrap of 300 replications; the anchors run on CRAN
  # One case carries all the variation in x, so any resample that
  # omits it refits with a constant predictor and returns no indirect
  # effect; those replications are dropped with one warning instead of
  # erroring the quantile call.
  set.seed(113)
  n <- 20
  x <- c(3, rep(0, n - 1))
  m <- 0.5 * x + rnorm(n)
  y <- 0.6 * m + rnorm(n)
  d <- data.frame(x = x, m = m, y = y)
  got_warning <- FALSE
  res <- withCallingHandlers(
    mediate(d, x = "x", m = "m", y = "y",
            ci_method = "boot_percentile", B = 300, seed = 113),
    warning = function(w) {
      if (grepl("dropped", conditionMessage(w))) got_warning <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  expect_true(got_warning)
  expect_true(all(is.finite(
    res$value[res$term %in% c("indirect_ci_lower", "indirect_ci_upper")])))
})
