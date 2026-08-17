# R2_mixed_effects_decomposition() implements the Rights & Sterba (2019) framework; validated against
# the authors' reference implementation (the r2mlm package) and, for the
# omnibus/marginal measures, against the Nakagawa marginal and conditional R2.
# The reference values below are pinned constants; the live oracle comparisons
# are collected in tools/oracle_checks.R.

test_that("R2_mixed_effects_decomposition() matches r2mlm on a non-cluster-mean-centered model", {
  skip_if_not_installed("lme4")
  fit <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                    data = lme4::sleepstudy)
  d <- R2_mixed_effects_decomposition(fit)
  v <- d$value; names(v) <- d$term
  expect_s3_class(d, "dmar_tbl")
  # non-CMC returns the five total-variance measures
  expect_setequal(d$term, c("total_f", "total_v", "total_m", "total_fv", "total_fvm"))
  # Pinned from r2mlm::r2mlm (r2mlm 0.3.8, 2026-08-09); live comparison in tools/oracle_checks.R.
  expect_equal(unname(v["total_f"]),   0.2785130443514864,  tolerance = 1e-8)
  expect_equal(unname(v["total_v"]),   0.08915267098225658, tolerance = 1e-8)
  expect_equal(unname(v["total_m"]),   0.431653652335793,   tolerance = 1e-8)
  expect_equal(unname(v["total_fv"]),  0.367665715333743,   tolerance = 1e-8)
  expect_equal(unname(v["total_fvm"]), 0.799319367669536,   tolerance = 1e-8)
})

test_that("R2_mixed_effects_decomposition() matches r2mlm on all 12 measures for a cluster-mean-centered model", {
  skip_if_not_installed("lme4")
  set.seed(113)
  J <- 60; nj <- 15; N <- J * nj; cl <- rep(1:J, each = nj)
  x <- rnorm(N); xbar <- ave(x, cl); xw <- x - xbar; w <- xbar
  u0 <- rnorm(J, 0, 1.2)[cl]; u1 <- rnorm(J, 0, .5)[cl]
  y <- 2 + 0.5 * xw + 0.8 * w + u0 + u1 * xw + rnorm(N, 0, 1)
  dd <- data.frame(y, xw, w, cl = factor(cl))
  fit <- suppressWarnings(
    lme4::lmer(y ~ xw + w + (xw | cl), data = dd))
  d <- R2_mixed_effects_decomposition(fit)
  v <- d$value; names(v) <- d$term
  # Pinned from r2mlm::r2mlm (r2mlm 0.3.8, 2026-08-09); live comparison in tools/oracle_checks.R.
  ref <- c(total_f1 = 0.0597275798372667, total_f2 = 0.00684771825130456,
           total_v = 0.0926289859730671, total_m = 0.474518079504468,
           total_f = 0.06657529808857129, total_fv = 0.159204284061638,
           total_fvm = 0.633722363566106, within_f1 = 0.115163210561151,
           within_v = 0.178601768977526, within_fv = 0.293764979538677,
           between_f2 = 0.0142256019917286, between_m = 0.985774398008271)
  expect_setequal(d$term, names(ref))
  expect_equal(max(abs(v[names(ref)] - ref)), 0, tolerance = 1e-8)
  # the variance decomposition is carried as an attribute and sums to 1
  expect_equal(sum(attr(d, "decomposition")[1, ]), 1, tolerance = 1e-8)
})

test_that("R2_mixed_effects_decomposition() total_f / total_fvm coincide with Nakagawa marginal / conditional", {
  skip_if_not_installed("lme4")
  fit <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
  d <- R2_mixed_effects_decomposition(fit); v <- d$value; names(v) <- d$term
  # Pinned from performance::r2_nakagawa (performance 0.16.0, 2026-08-09); live comparison in tools/oracle_checks.R.
  expect_equal(unname(v["total_f"]),   0.2798856367631553, tolerance = 1e-6)
  expect_equal(unname(v["total_fvm"]), 0.7042554523353618, tolerance = 1e-6)
})

test_that("R2_mixed_effects_decomposition() rejects unsupported input", {
  expect_error(R2_mixed_effects_decomposition(lm(mpg ~ wt, data = mtcars)),
               "merMod|lme")
})

test_that("R2_mixed_effects_decomposition() handles separate random intercept and slope terms", {
  skip_if_not_installed("lme4")

  # CRITICAL-03 regression. For (1 | g) + (0 + x | g) the random effects span
  # two VarCorr components (named 'Subject' and 'Subject.1'), while ranef()[[1]]
  # merges both columns; reading a single VarCorr component made the variance
  # decomposition non-conformable and R2_mixed_effects_decomposition() aborted. The reassembled
  # covariance must reproduce the r2mlm reference for the same split model.
  fit <- lme4::lmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
                    data = lme4::sleepstudy, REML = TRUE)
  d <- R2_mixed_effects_decomposition(fit)
  v <- d$value; names(v) <- d$term
  # Pinned from r2mlm::r2mlm (r2mlm 0.3.8, 2026-08-09); live comparison in tools/oracle_checks.R.
  expect_equal(unname(v["total_f"]),   0.2828350562331848,  tolerance = 1e-8)
  expect_equal(unname(v["total_v"]),   0.09256690012171188, tolerance = 1e-8)
  expect_equal(unname(v["total_m"]),   0.4212253143452074,  tolerance = 1e-8)
  expect_equal(unname(v["total_fv"]),  0.3754019563548966,  tolerance = 1e-8)
  expect_equal(unname(v["total_fvm"]), 0.7966272707001041,  tolerance = 1e-8)
})

test_that("R2_mixed_effects_decomposition() treats the double-bar form as the split form", {
  skip_if_not_installed("lme4")

  # (x || g) is lme4 shorthand for (1 | g) + (0 + x | g); both reach the
  # multi-component reassembly path and must give identical measures.
  ss <- lme4::sleepstudy
  fit_bar   <- lme4::lmer(Reaction ~ Days + (Days || Subject),
                          data = ss, REML = TRUE)
  fit_split <- lme4::lmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
                          data = ss, REML = TRUE)
  expect_equal(R2_mixed_effects_decomposition(fit_bar)$value, R2_mixed_effects_decomposition(fit_split)$value, tolerance = 1e-8)
})
