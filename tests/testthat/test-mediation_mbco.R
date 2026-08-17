# Tests for mediation_mbco(): structure and agreement with lm/lavaan on
# a small simulated model, then exact replication of the empirical
# example in Tofighi and Kelley (2020), Psychological Methods, whose
# data are the memory experiment of MacKinnon, Valente, and Wurpts
# (2018). The replication runs from the exact sample moments (complete
# data maximum likelihood depends on the data only through the sample
# means and covariance matrix), so the published values are recovered
# without shipping the raw data.

# Exact moments of the memory study data (N = 369), variables
# instruction (0 = repetition, 1 = imagery), imagery, repetition,
# recall. Computed once from the study data; divisor N - 1.
.tk_vars <- c("instruction", "imagery", "repetition", "recall")
.tk_M <- c(instruction = 0.5094850949, imagery = 5.6612466125,
           repetition = 6.0785907859, recall = 12.0650406504)
.tk_S <- matrix(c(
   0.2505891363,  0.9121892306, -0.9531931189,  0.5428596677,
   0.9121892306,  8.7626517026, -4.7368917167,  5.1470926122,
  -0.9531931189, -4.7368917167,  8.0454371392, -2.7143646165,
   0.5428596677,  5.1470926122, -2.7143646165, 11.5283669141),
  nrow = 4L, dimnames = list(.tk_vars, .tk_vars))
.tk_N <- 369L

.tk_single <- "
  imagery ~ b1*instruction
  recall  ~ b2*imagery + b3*instruction
"
.tk_parallel <- "
  imagery    ~ b1*instruction
  repetition ~ b3*instruction
  recall     ~ b2*imagery + b4*repetition + b5*instruction
  imagery ~~ repetition
"

test_that("mediation_mbco matches lm and lavaan on a simple model", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  set.seed(113)
  n <- 150
  x <- rnorm(n)
  m <- 0.5 * x + rnorm(n)
  y <- 0.2 * x + 0.4 * m + rnorm(n)
  d <- data.frame(x = x, m = m, y = y)

  res <- mediation_mbco("m ~ x \n y ~ m + x", data = d,
                        x = "x", y = "y", ci_method = "wald")

  expect_s3_class(res, "dmar_mediation_mbco")
  expect_s3_class(res, "dmar_tbl")
  expect_identical(res$term, c("total_effect", "direct_effect",
                               "indirect_via_m"))
  expect_true(all(c("pathway", "estimate", "se", "ci_lower", "ci_upper",
                    "lrt", "df", "p_value", "delta_aic", "delta_bic")
                  %in% names(res)))

  # Maximum likelihood path coefficients equal the OLS coefficients.
  a_hat  <- unname(coef(lm(m ~ x, d))["x"])
  b_hat  <- unname(coef(lm(y ~ m + x, d))["m"])
  cp_hat <- unname(coef(lm(y ~ m + x, d))["x"])
  expect_equal(res$estimate[res$term == "indirect_via_m"],
               a_hat * b_hat, tolerance = 1e-4)
  expect_equal(res$estimate[res$term == "direct_effect"], cp_hat,
               tolerance = 1e-4)
  expect_equal(res$estimate[res$term == "total_effect"],
               cp_hat + a_hat * b_hat, tolerance = 1e-4)

  # The full-model deviance equals lavaan's for the same likelihood.
  lav <- lavaan::sem("m ~ x \n y ~ m + x", data = d,
                     meanstructure = TRUE, fixed.x = FALSE, se = "none")
  expect_equal(attr(res, "deviance"),
               -2 * as.numeric(lavaan::fitMeasures(lav, "logl")),
               tolerance = 1e-3)

  # Wald interval identity, test bookkeeping, information criteria.
  z <- qnorm(0.975)
  expect_equal(res$ci_lower, res$estimate - z * res$se,
               tolerance = 1e-10)
  expect_true(all(res$lrt >= 0))
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
  expect_identical(res$df, rep(1L, 3L))
  expect_equal(res$delta_aic, res$lrt - 2, tolerance = 1e-10)
  expect_equal(res$delta_bic, res$lrt - log(n), tolerance = 1e-10)

  # Broom verbs.
  td <- generics::tidy(res)
  expect_identical(names(td), c("term", "estimate", "se",
                                "statistic", "p_value", "ci_lower",
                                "ci_upper"))
  gl <- generics::glance(res)
  expect_identical(nrow(gl), 1L)
  expect_equal(gl$nobs, n)
  expect_equal(gl$AIC, gl$deviance + 2 * gl$npar, tolerance = 1e-8)

  # R2 attributes cover the endogenous variables.
  expect_named(attr(res, "R2"), c("m", "y"))
  expect_identical(dim(attr(res, "delta_R2")), c(2L, 3L))
})

test_that("single-mediator replication of Tofighi and Kelley (2020)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  res <- mediation_mbco(.tk_single, S = .tk_S, M = .tk_M, N = .tk_N,
                        x = "instruction", y = "recall",
                        ci_method = "monte_carlo", B = 1e5, seed = 113)

  ind <- res[res$term == "indirect_via_imagery", ]
  # Published: estimate 2.121, SE 0.276, 95% Monte Carlo CI
  # [1.600, 2.682].
  expect_equal(ind$estimate, 2.121041, tolerance = 1e-3)
  expect_equal(ind$se, 0.276, tolerance = 2e-3)
  expect_equal(ind$ci_lower, 1.600, tolerance = 0.03)
  expect_equal(ind$ci_upper, 2.686, tolerance = 0.03)

  # The likelihood ratio statistic is defined by the best-fitting null
  # model. The null set {b1 b2 = 0} is the union of two branches; the
  # better one sets b2 = 0 and gives LRT = 72.544 (the value 175.766 in
  # the published example is the other, worse-fitting branch, which the
  # NPSOL run stopped on). The exact value here comes from the
  # regression identity LRT = N log(RSS_null / RSS_full) for the
  # recall equation.
  expect_equal(ind$lrt, 72.54357, tolerance = 0.02)
  expect_equal(ind$df, 1L)
  expect_true(ind$p_value < 1e-15)

  # Total = direct + indirect (an identity in this just-identified
  # model).
  expect_equal(res$estimate[res$term == "total_effect"],
               res$estimate[res$term == "direct_effect"] + ind$estimate,
               tolerance = 1e-6)

  # On the retained branch the imagery equation is untouched, so its
  # R2 is essentially unchanged while the recall R2 drops.
  dr2 <- attr(res, "delta_R2")
  expect_lt(abs(dr2["imagery", "indirect_via_imagery"]), 0.01)
  expect_gt(dr2["recall", "indirect_via_imagery"], 0.10)

  # Full-model R2 as published (.38 imagery, .26 recall).
  r2 <- attr(res, "R2")
  expect_equal(unname(r2["imagery"]), 0.38, tolerance = 0.01)
  expect_equal(unname(r2["recall"]), 0.26, tolerance = 0.01)
})

test_that("profile likelihood interval matches the published range", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  res <- mediation_mbco(.tk_single, S = .tk_S, M = .tk_M, N = .tk_N,
                        x = "instruction", y = "recall",
                        ci_method = "profile_likelihood")
  ind <- res[res$term == "indirect_via_imagery", ]
  expect_equal(ind$ci_lower, 1.604, tolerance = 0.02)
  expect_equal(ind$ci_upper, 2.690, tolerance = 0.02)
})

test_that("parallel two-mediator replication (RQ2 and RQ3)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  res <- mediation_mbco(
    .tk_parallel, S = .tk_S, M = .tk_M, N = .tk_N,
    x = "instruction", y = "recall",
    hypotheses = c(imagery_minus_repetition =
                     "indirect_via_imagery - indirect_via_repetition"),
    ci_method = "monte_carlo", B = 1e5, seed = 113)

  expect_identical(res$term, c(
    "total_effect", "direct_effect", "total_indirect",
    "indirect_via_imagery", "indirect_via_repetition",
    "imagery_minus_repetition"))

  # Published RQ2: specific indirect through repetition -0.082
  # (SE 0.285), LRT = 0.083, p = .773; through imagery 2.139
  # (SE 0.284).
  rep_row <- res[res$term == "indirect_via_repetition", ]
  expect_equal(rep_row$estimate, -0.082, tolerance = 2e-3)
  expect_equal(rep_row$se, 0.285, tolerance = 3e-3)
  expect_equal(rep_row$lrt, 0.0831, tolerance = 0.01)
  expect_equal(rep_row$p_value, 0.773, tolerance = 0.01)
  expect_equal(res$estimate[res$term == "indirect_via_imagery"],
               2.139, tolerance = 2e-3)

  # Published RQ3: contrast 2.222 (SE 0.445), LRT = 25.828,
  # p = 3.7e-07.
  con <- res[res$term == "imagery_minus_repetition", ]
  expect_equal(con$estimate, 2.222, tolerance = 2e-3)
  expect_equal(con$se, 0.445, tolerance = 2e-3)
  expect_equal(con$lrt, 25.828, tolerance = 0.05)
  expect_lt(con$p_value, 1e-5)
  expect_gt(con$p_value, 1e-8)

  # Additivity of the enumerated effects.
  expect_equal(res$estimate[res$term == "total_indirect"],
               sum(res$estimate[res$term %in%
                     c("indirect_via_imagery",
                       "indirect_via_repetition")]),
               tolerance = 1e-8)
  expect_equal(res$estimate[res$term == "total_effect"],
               res$estimate[res$term == "direct_effect"] +
                 res$estimate[res$term == "total_indirect"],
               tolerance = 1e-6)
})

test_that("latent-mediator models fit and agree with lavaan", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  set.seed(113)
  n <- 400
  x <- rnorm(n)
  eta <- 0.6 * x + rnorm(n)
  d <- data.frame(
    x = x,
    m1 = 0.8 * eta + rnorm(n, 0, 0.6),
    m2 = 0.7 * eta + rnorm(n, 0, 0.6),
    m3 = 0.9 * eta + rnorm(n, 0, 0.6),
    yy = 0.4 * eta + 0.1 * x + rnorm(n))
  model <- "
    M  =~ m1 + m2 + m3
    M  ~ a*x
    yy ~ b*M + cp*x
    ind := a*b
  "
  res <- mediation_mbco(model, data = d, x = "x", y = "yy",
                        ci_method = "wald")

  # The user-defined quantity is carried along with the enumerated
  # effects, and duplicates the enumerated indirect pathway.
  expect_true(all(c("indirect_via_M", "ind") %in% res$term))
  expect_equal(res$estimate[res$term == "ind"],
               res$estimate[res$term == "indirect_via_M"],
               tolerance = 1e-8)

  lav <- lavaan::sem(model, data = d, meanstructure = TRUE,
                     fixed.x = FALSE, se = "none")
  expect_equal(attr(res, "deviance"),
               -2 * as.numeric(lavaan::fitMeasures(lav, "logl")),
               tolerance = 1e-2)
  est <- lavaan::parameterEstimates(lav)
  lav_ind <- est$est[est$label == "ind"]
  expect_equal(res$estimate[res$term == "ind"], lav_ind,
               tolerance = 1e-3)
  expect_true(res$lrt[res$term == "indirect_via_M"] >= 0)
})

test_that("input validation errors are informative", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  d <- data.frame(x = rnorm(30), m = rnorm(30), y = rnorm(30))
  expect_error(
    mediation_mbco("m ~ x \n y ~ m + x", data = d, S = diag(3)),
    "not both")
  expect_error(
    mediation_mbco("m ~ x \n y ~ m + x",
                   S = cov(d)),
    "'N'")
  expect_error(
    mediation_mbco("m ~ x \n y ~ m + x", data = d, x = "w", y = "y"),
    "no regression path")
  # A feedback loop makes pathway enumeration ill defined.
  d_cyc <- data.frame(x = rnorm(60))
  d_cyc$m <- 0.5 * d_cyc$x + rnorm(60)
  d_cyc$y <- 0.4 * d_cyc$m + rnorm(60)
  expect_error(
    mediation_mbco("m ~ x + y \n y ~ m", data = d_cyc,
                   x = "x", y = "y"),
    "nonrecursive")
  expect_error(
    mediation_mbco("y ~ x", data = d, x = "x", y = "y"),
    "no indirect pathway")
  expect_error(
    mediation_mbco("m ~ x \n y ~ m + x", data = d, x = "x", y = "y",
                   hypotheses = c(h = "b1 * nope")),
    "unknown parameter")
  # Ambiguous source variables require explicit x and y.
  d2 <- data.frame(x1 = rnorm(30), x2 = rnorm(30))
  d2$m <- d2$x1 + rnorm(30)
  d2$y <- d2$m + d2$x2 + rnorm(30)
  expect_error(
    mediation_mbco("m ~ x1 \n y ~ m + x2", data = d2),
    "ambiguous")
})

test_that("multiple groups: per-group effects and difference tests", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  set.seed(113)
  n <- 250
  make_grp <- function(a, b, cp, label) {
    x <- rnorm(n)
    m <- a * x + rnorm(n)
    y <- cp * x + b * m + rnorm(n)
    data.frame(x = x, m = m, y = y, g = label)
  }
  d <- rbind(make_grp(0.6, 0.5, 0.2, "one"),
             make_grp(0.2, 0.1, 0.2, "two"))

  res <- mediation_mbco("m ~ x \n y ~ m + x", data = d, group = "g",
                        x = "x", y = "y", ci_method = "wald")
  expect_identical(attr(res, "groups"), c("one", "two"))
  expect_true(all(c("indirect_via_m_one", "indirect_via_m_two",
                    "indirect_via_m_two_minus_one",
                    "direct_effect_one", "direct_effect_two",
                    "direct_effect_two_minus_one") %in% res$term))

  lav <- lavaan::sem("m ~ x \n y ~ m + x", data = d, group = "g",
                     meanstructure = TRUE, fixed.x = FALSE,
                     se = "none")
  expect_equal(attr(res, "deviance"),
               -2 * as.numeric(lavaan::fitMeasures(lav, "logl")),
               tolerance = 1e-3)
  pe <- lavaan::parameterEstimates(lav)
  ab <- function(g) {
    a <- pe$est[pe$lhs == "m" & pe$op == "~" & pe$rhs == "x" &
                  pe$group == g]
    b <- pe$est[pe$lhs == "y" & pe$op == "~" & pe$rhs == "m" &
                  pe$group == g]
    a * b
  }
  expect_equal(res$estimate[res$term == "indirect_via_m_one"], ab(1L),
               tolerance = 1e-4)
  expect_equal(res$estimate[res$term == "indirect_via_m_two"], ab(2L),
               tolerance = 1e-4)
  diff_row <- res[res$term == "indirect_via_m_two_minus_one", ]
  expect_equal(diff_row$estimate, ab(2L) - ab(1L), tolerance = 1e-4)
  # The groups were generated with truly different indirect effects.
  expect_gt(diff_row$lrt, 10)
  expect_lt(diff_row$p_value, 0.001)
  expect_true(all(c("m (one)", "y (one)", "m (two)", "y (two)") %in%
                    names(attr(res, "R2"))))

  # Grouped summary statistics reproduce the raw-data analysis.
  d1 <- d[d$g == "one", c("x", "m", "y")]
  d2 <- d[d$g == "two", c("x", "m", "y")]
  res_s <- mediation_mbco(
    "m ~ x \n y ~ m + x",
    S = list(one = cov(d1), two = cov(d2)),
    M = list(one = colMeans(d1), two = colMeans(d2)),
    N = c(one = nrow(d1), two = nrow(d2)),
    x = "x", y = "y", ci_method = "wald")
  expect_equal(res_s$estimate, res$estimate, tolerance = 1e-5)
  expect_equal(attr(res_s, "deviance"), attr(res, "deviance"),
               tolerance = 1e-4)

  # A label shared across groups is a cross-group equality constraint;
  # the degenerate difference row is dropped, the others remain.
  # (lavaan issues an advisory about the shared-label semantics.)
  res_eq <- suppressWarnings(
    mediation_mbco("m ~ b1*x \n y ~ b2*m + x", data = d,
                   group = "g", x = "x", y = "y",
                   ci_method = "wald"))
  expect_false("indirect_via_m_two_minus_one" %in% res_eq$term)
  expect_true("direct_effect_two_minus_one" %in% res_eq$term)
})

test_that("guardrails warn on binary endogenous and interaction terms", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  set.seed(113)
  n <- 150
  db <- data.frame(x = rnorm(n))
  db$m <- as.numeric(db$x + rnorm(n) > 0)
  db$y <- 0.6 * db$m + 0.2 * db$x + rnorm(n)
  w_bin <- testthat::capture_warnings(
    mediation_mbco("m ~ x \n y ~ m + x", data = db, x = "x", y = "y",
                   ci_method = "wald"))
  expect_true(any(grepl("two distinct values", w_bin)))

  dp <- data.frame(x = rnorm(n), w = rnorm(n))
  dp$xw <- dp$x * dp$w
  dp$m <- 0.5 * dp$x + 0.3 * dp$w + 0.4 * dp$xw + rnorm(n)
  dp$y <- 0.5 * dp$m + 0.2 * dp$x + 0.1 * dp$w + rnorm(n)
  w_int <- testthat::capture_warnings(
    mediation_mbco("m ~ x + w + xw \n y ~ m + x + w", data = dp,
                   x = "x", y = "y", ci_method = "wald"))
  expect_true(any(grepl("Interaction term", w_int)))

  # Declaring the moderator resolves the interaction warning.
  res_mod <- mediation_mbco("m ~ x + w + xw \n y ~ m + x + w",
                            data = dp, x = "x", y = "y",
                            moderator = "w", ci_method = "wald")
  expect_true("indirect_via_m_moderation" %in% res_mod$term)

  # A declared moderator with no interaction in the model is ignored,
  # with a warning.
  d0 <- data.frame(x = rnorm(n), w = rnorm(n))
  d0$m <- 0.5 * d0$x + rnorm(n)
  d0$y <- 0.5 * d0$m + rnorm(n)
  expect_warning(
    mediation_mbco("m ~ x \n y ~ m + x", data = d0, x = "x", y = "y",
                   moderator = "w", ci_method = "wald"),
    "is ignored")

  # Moderator and group cannot yet be combined.
  d0$g <- rep(c("a", "b"), length.out = n)
  expect_error(
    mediation_mbco("m ~ x \n y ~ m + x", data = d0, group = "g",
                   moderator = "w"),
    "one moderator at a time")
})

test_that("moderator probing: first-stage moderation and the index", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  set.seed(113)
  n <- 400
  x <- rnorm(n)
  w <- rnorm(n)
  m <- 0.5 * x + 0.3 * w + 0.4 * x * w + rnorm(n)
  y <- 0.5 * m + 0.2 * x + 0.1 * w + rnorm(n)
  d <- data.frame(x = x, w = w, m = m, y = y)

  # The 'x:w' interaction syntax is rewritten into a product column.
  res <- mediation_mbco("m ~ x + w + x:w \n y ~ m + x + w", data = d,
                        x = "x", y = "y", moderator = "w",
                        ci_method = "wald")
  expect_true(all(c("indirect_via_m_at_low", "indirect_via_m_at_mean",
                    "indirect_via_m_at_high",
                    "indirect_via_m_moderation",
                    "total_effect_moderation") %in% res$term))
  # The direct effect is unmoderated here: one plain row, no probes.
  expect_true("direct_effect" %in% res$term)
  expect_false(any(grepl("^direct_effect_at_", res$term)))

  # Regression algebra anchors for the conditional effects and the
  # index of moderated mediation.
  d$xw <- d$x * d$w
  fm <- lm(m ~ x + w + xw, d)
  fy <- lm(y ~ m + x + w, d)
  a1 <- unname(coef(fm)["x"])
  a3 <- unname(coef(fm)["xw"])
  b  <- unname(coef(fy)["m"])
  mw <- mean(d$w)
  sw <- sd(d$w)
  expect_equal(res$estimate[res$term == "indirect_via_m_at_low"],
               (a1 + a3 * (mw - sw)) * b, tolerance = 1e-4)
  expect_equal(res$estimate[res$term == "indirect_via_m_at_mean"],
               (a1 + a3 * mw) * b, tolerance = 1e-4)
  expect_equal(res$estimate[res$term == "indirect_via_m_at_high"],
               (a1 + a3 * (mw + sw)) * b, tolerance = 1e-4)
  expect_equal(res$estimate[res$term == "indirect_via_m_moderation"],
               a3 * b, tolerance = 1e-4)
  # The data were generated with strong first-stage moderation.
  expect_lt(res$p_value[res$term == "indirect_via_m_moderation"],
            1e-4)

  # Named probe values are honored and labeled.
  res2 <- mediation_mbco("m ~ x + w + x:w \n y ~ m + x + w", data = d,
                         x = "x", y = "y", moderator = "w",
                         probe_values = c(minus1 = -1, plus1 = 1),
                         ci_method = "wald")
  expect_true(all(c("indirect_via_m_at_minus1",
                    "indirect_via_m_at_plus1") %in% res2$term))
  expect_equal(res2$estimate[res2$term == "indirect_via_m_at_plus1"],
               (a1 + a3) * b, tolerance = 1e-4)
})

test_that("moderator probing: two-stage moderation, joint constancy", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  set.seed(113)
  n <- 400
  x <- rnorm(n)
  w <- rnorm(n)
  m <- 0.5 * x + 0.3 * w + 0.4 * x * w + rnorm(n)
  y <- 0.5 * m + 0.3 * m * w + 0.2 * x + 0.1 * w + rnorm(n)
  d <- data.frame(x = x, w = w, m = m, y = y)

  res <- mediation_mbco("m ~ x + w + x:w \n y ~ m + x + w + m:w",
                        data = d, x = "x", y = "y", moderator = "w",
                        ci_method = "wald")
  # The conditional indirect effect is quadratic in w: linear and
  # quadratic coefficient rows plus a joint constancy test.
  expect_true(all(c("indirect_via_m_moderation_w",
                    "indirect_via_m_moderation_w2",
                    "indirect_via_m_moderation") %in% res$term))
  joint <- res[res$term == "indirect_via_m_moderation", ]
  expect_identical(joint$df, 2L)
  expect_true(is.na(joint$estimate))
  expect_lt(joint$p_value, 1e-6)

  # Polynomial coefficient anchors: w carries a3 b1 + a1 b3, w^2
  # carries a3 b3.
  d$xw <- d$x * d$w
  d$mw <- d$m * d$w
  fm <- lm(m ~ x + w + xw, d)
  fy <- lm(y ~ m + x + w + mw, d)
  a1 <- unname(coef(fm)["x"])
  a3 <- unname(coef(fm)["xw"])
  b1 <- unname(coef(fy)["m"])
  b3 <- unname(coef(fy)["mw"])
  expect_equal(res$estimate[res$term == "indirect_via_m_moderation_w"],
               a3 * b1 + a1 * b3, tolerance = 1e-4)
  expect_equal(
    res$estimate[res$term == "indirect_via_m_moderation_w2"],
    a3 * b3, tolerance = 1e-4)
})

test_that("joint hypotheses match the nested-model likelihood ratio", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  set.seed(113)
  n <- 300
  x <- rnorm(n)
  w <- rnorm(n)
  m <- 0.5 * x + rnorm(n)          # w truly irrelevant
  y <- 0.5 * m + 0.2 * x + rnorm(n)
  d <- data.frame(x = x, w = w, m = m, y = y)
  model <- "m ~ a1*x + a2*w \n y ~ b*m + cp*x + cw*w"

  res <- mediation_mbco(model, data = d, x = "x", y = "y",
                        hypotheses = list(w_irrelevant = c("a2", "cw")),
                        ci_method = "wald")
  row <- res[res$term == "w_irrelevant", ]
  expect_identical(row$df, 2L)
  expect_true(is.na(row$estimate))
  expect_equal(row$delta_aic, row$lrt - 4, tolerance = 1e-10)

  # For linear constraints the joint MBCO test must equal the ordinary
  # nested-model likelihood ratio, here computed by lavaan directly.
  lav_full <- lavaan::sem(model, data = d, meanstructure = TRUE,
                          fixed.x = FALSE, se = "none")
  lav_null <- lavaan::sem("m ~ a1*x + 0*w \n y ~ b*m + cp*x + 0*w",
                          data = d, meanstructure = TRUE,
                          fixed.x = FALSE, se = "none")
  lrt_lav <- -2 * (as.numeric(lavaan::fitMeasures(lav_null, "logl")) -
                     as.numeric(lavaan::fitMeasures(lav_full, "logl")))
  expect_equal(row$lrt, lrt_lav, tolerance = 0.01)
})

test_that("moderated mediation with a continuous moderator via ':='", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  set.seed(113)
  n <- 400
  x <- rnorm(n)
  w <- rnorm(n)
  xw <- x * w
  m <- 0.5 * x + 0.3 * w + 0.4 * xw + rnorm(n)
  y <- 0.5 * m + 0.2 * x + 0.1 * w + rnorm(n)
  d <- data.frame(x = x, w = w, xw = xw, m = m, y = y)
  mod <- "
    m ~ a1*x + a2*w + a3*xw
    y ~ b*m + cp*x + cw*w
    index    := a3*b
    ind_low  := (a1 - a3)*b
    ind_mean := a1*b
    ind_high := (a1 + a3)*b
  "
  # Omitting x and y skips the (here inappropriate) automatic
  # enumeration; only the defined conditional quantities are tested.
  res <- mediation_mbco(mod, data = d, ci_method = "wald")
  expect_identical(res$term,
                   c("index", "ind_low", "ind_mean", "ind_high"))

  fm <- lm(m ~ x + w + xw, d)
  fy <- lm(y ~ m + x + w, d)
  a1h <- unname(coef(fm)["x"])
  a3h <- unname(coef(fm)["xw"])
  bh  <- unname(coef(fy)["m"])
  expect_equal(res$estimate[res$term == "index"], a3h * bh,
               tolerance = 1e-4)
  expect_equal(res$estimate[res$term == "ind_low"], (a1h - a3h) * bh,
               tolerance = 1e-4)
  expect_equal(res$estimate[res$term == "ind_high"], (a1h + a3h) * bh,
               tolerance = 1e-4)
  expect_true(all(res$lrt >= 0))
  # The data were generated with a strong index of moderated mediation.
  expect_lt(res$p_value[res$term == "index"], 0.001)
})

test_that("the caller's RNG state is preserved", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")

  set.seed(20)
  before <- .Random.seed
  invisible(mediation_mbco(.tk_single, S = .tk_S, M = .tk_M, N = .tk_N,
                           x = "instruction", y = "recall",
                           ci_method = "monte_carlo", B = 1000,
                           seed = 7))
  expect_identical(before, .Random.seed)
})
