# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

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
                                "ci_upper", "conf_level"))
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

test_that("the lines above an MBCO table reach a knitted report", {
  # print.dmar_mediation_mbco() wrote these lines above the table itself,
  # so knit_print() and as_kable(), which never call print(), dropped
  # them. They are now the header block of the one footer every display
  # reads, and the console prints them exactly as before.
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  skip_if_not_installed("knitr")
  set.seed(113)
  n <- 150
  x <- rnorm(n)
  m <- 0.5 * x + rnorm(n)
  y <- 0.2 * x + 0.4 * m + rnorm(n)
  res <- mediation_mbco("m ~ x \n y ~ m + x", data = data.frame(x, m, y),
                        x = "x", y = "y", ci_method = "wald")
  fmt <- function(v) format(v, nsmall = 3, digits = 10)
  r2 <- attr(res, "R2")
  header <- c(
    "Mediation tests via model-based constrained optimization (MBCO)",
    "  Effects of x on y; N = 150; optimizer: SLSQP",
    paste0("  Full model: deviance = ", fmt(attr(res, "deviance")),
           ", AIC = ", fmt(attr(res, "aic")),
           ", BIC = ", fmt(attr(res, "bic"))),
    paste0("  R2 (full model): m = ", formatC(r2[["m"]], digits = 3,
                                              format = "f"),
           ", y = ", formatC(r2[["y"]], digits = 3, format = "f")),
    "  Each p_value is from the MBCO likelihood ratio test of the",
    "  null model constraining that effect to zero. Causal readings",
    "  rest on the no omitted confounder assumption.")
  printed <- capture.output(print(res))
  expect_identical(printed[seq_len(length(header) + 1L)], c(header, ""))
  expect_match(printed[length(header) + 2L], "^ pathway +term")

  footer <- DMAR:::.dmar_tbl_footer(res)
  expect_identical(
    strsplit(paste(footer[names(footer) == "header"], collapse = "\n"),
             "\n", fixed = TRUE)[[1L]],
    header)

  # knit_print() sets the lines as paragraphs before the table, and
  # as_kable() puts them first in the caption.
  squish <- function(s) trimws(gsub("\\s+", " ", s))
  knitted <- as.character(knit_print.dmar_tbl(res))
  table_at <- regexpr("|pathway", knitted, fixed = TRUE)
  expect_gt(table_at, 0L)
  before <- squish(substr(knitted, 1L, table_at - 1L))
  expect_identical(before, squish(paste(header, collapse = " ")))
  expect_match(knitted, "^\n\nMediation tests via model-based constrained ")
  expect_match(knitted, "rest on the no omitted confounder assumption.\n\n|pathway",
               fixed = TRUE)
  kab <- as.character(as_kable(res, format = "pipe"))
  expect_match(kab[1L], paste0(
    "Table: Mediation tests via model-based constrained optimization ",
    "(MBCO); Effects of x on y; N = 150; optimizer: SLSQP; Full model: ",
    "deviance = "), fixed = TRUE)
  expect_match(kab[1L], paste0(
    "Causal readings rest on the no omitted confounder assumption. ",
    "Confidence level: 95%$"))

  # A column subset loses the fit attributes; only what exists is written.
  sub <- capture.output(print(res[, c("term", "estimate")]))
  expect_identical(sub[1:5], c(header[1L], header[5:7], ""))
})

# Least squares oracle for the single mediator model "m ~ a*x; y ~ b*m +
# cp*x" with observed variables: the likelihood factors into x, m | x, and
# y | m, x, so the likelihood ratio statistic for a = 0 is
# N log(SSE(m ~ 1) / SSE(m ~ x)), that for b = 0 is
# N log(SSE(y ~ x) / SSE(y ~ m + x)), and that for a = 0 and b = 0 jointly
# is their sum. Coded apart from the package (lm.fit only).
.mbco_ls_branches <- function(d) {
  N <- nrow(d)
  sse <- function(X, v) sum(stats::lm.fit(X, v)$residuals^2)
  one <- matrix(1, N, 1L)
  c(a = N * log(sse(one, d$m) / sse(cbind(one, d$x), d$m)),
    b = N * log(sse(cbind(one, d$x), d$y) / sse(cbind(one, d$x, d$m), d$y)))
}
.mbco_dr38_data <- function() {
  set.seed(731)
  d <- data.frame(x = stats::rnorm(100))
  d$m <- .2 * d$x + stats::rnorm(100)
  d$y <- .2 * d$x + .3 * d$m + stats::rnorm(100)
  d
}

test_that("redundant joint hypotheses are refused before any fit (DR38)", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  # Through the first builds of 1.1.0 c("a", "a") and c("a", "2*a") gave
  # the statistic of "a" alone (4.980) with 2 df, p = .0829 where "a"
  # gave .0256; the refusal now names the redundant expression. The cases
  # go beyond a repeated string: a repeat after normalization, a scalar
  # multiple, a sign change, and a sum of two others.
  d <- .mbco_dr38_data()
  model <- "m ~ a*x\ny ~ b*m + cp*x"
  cases <- list(
    list(h = c("a", "a"), redundant = "a", rank = 1L),
    list(h = c("a", " a "), redundant = "a", rank = 1L),
    list(h = c("a", "2*a"), redundant = "2*a", rank = 1L),
    list(h = c("a - b", "b - a"), redundant = "b - a", rank = 1L),
    list(h = c("a", "b", "a + b"), redundant = "a + b", rank = 2L),
    list(h = c("a*b", "2*(a*b)"), redundant = "2*(a*b)", rank = 1L),
    list(h = c("a + cp", "b", "a + cp + 3*b"), redundant = "a + cp + 3*b",
         rank = 2L))
  for (cs in cases) {
    err <- tryCatch(
      mediation_mbco(model, data = d, hypotheses = list(joint = cs$h),
                     ci_method = "wald"),
      error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "hypotheses")
    expect_identical(err$redundant, cs$redundant)
    expect_identical(err$rank, cs$rank)
    expect_identical(err$n_constraints, length(cs$h))
    expect_match(conditionMessage(err), "Nothing was fitted.", fixed = TRUE)
  }
  # A single expression that restricts nothing.
  err <- tryCatch(
    mediation_mbco(model, data = d, hypotheses = c(none = "a - a"),
                   ci_method = "wald"),
    error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$rank, 0L)
})

test_that("equivalent independent hypotheses give the same test, as the oracle (DR38)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  d <- .mbco_dr38_data()
  N <- nrow(d)
  res <- mediation_mbco("m ~ a*x\ny ~ b*m + cp*x", data = d,
                        hypotheses = list(single = "a", scaled = "2*a",
                                          pair = c("a", "b"),
                                          rotated = c("a + b", "a - b")),
                        ci_method = "wald")
  row <- function(tm) res[res$term == tm, c("lrt", "df", "p_value",
                                            "delta_aic", "delta_bic")]
  orc <- .mbco_ls_branches(d)
  # One constraint: a and 2*a are the same restriction.
  expect_equal(row("scaled"), row("single"), tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_equal(row("single")$lrt, unname(orc["a"]), tolerance = 1e-4)
  expect_identical(row("single")$df, 1L)
  # Two constraints: c(a, b) and c(a + b, a - b) state the same null set.
  expect_equal(row("rotated"), row("pair"), tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_identical(row("pair")$df, 2L)
  expect_equal(row("pair")$lrt, sum(orc), tolerance = 1e-4)
  expect_equal(row("pair")$p_value,
               stats::pchisq(sum(orc), 2, lower.tail = FALSE), tolerance = 1e-4)
  expect_equal(row("pair")$delta_aic, sum(orc) - 4, tolerance = 1e-4)
  expect_equal(row("pair")$delta_bic, sum(orc) - 2 * log(N), tolerance = 1e-4)
})

# Exact null-set algebra for constraints that are products of parameter
# labels (times nonzero constants, with positive powers), coded apart from
# the package: the null set of such constraints is the union of the
# coordinate subspaces theta_H = 0 over the hitting sets H of their factor
# sets, so a constraint is 0 wherever the others are exactly when its
# factor set meets every hitting set of theirs. Returns the constraints
# found redundant in order, each against the ones still kept, as the
# package names them.
.mbco_exact_redundant <- function(exprs) {
  fac <- lapply(exprs, function(e) sort(unique(all.vars(parse(text = e)[[1]]))))
  labels <- sort(unique(unlist(fac)))
  hits <- function(F) {
    out <- list()
    for (m in 0:length(labels)) {
      for (H in utils::combn(labels, m, simplify = FALSE)) {
        if (all(vapply(F, function(f) any(f %in% H), logical(1)))) {
          out[[length(out) + 1]] <- H
        }
      }
    }
    out
  }
  active <- seq_along(exprs)
  red <- integer(0)
  for (k in seq_along(exprs)) {
    o <- setdiff(active, k)
    if (!length(o)) next
    if (all(vapply(hits(fac[o]), function(H) any(fac[[k]] %in% H),
                   logical(1)))) {
      red <- c(red, k)
      active <- o
    }
  }
  red
}

test_that("a constraint that is 0 wherever the others are is refused (DR38)", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  # These sets pass the Jacobian check (their constraints are independent
  # off the null set) but their null sets impose one restriction, not two
  # or three: c("a*b", "b") is the null set b = 0. Through the first round
  # of the DR38 repair c("a*b", "b") gave the statistic of "b" (6.908) with
  # 2 df, p = .0316 where the 1 df test gives .00858, and c("a", "a*cp")
  # gave .0829 where "a" gives .0256. The expected redundant constraints
  # come from the exact algebra above, not from the package.
  d <- .mbco_dr38_data()
  model <- "m ~ a*x\ny ~ b*m + cp*x"
  for (h in list(c("a*b", "b"), c("b", "a*b"), c("a*b", "a"),
                 c("a", "a*cp"), c("a*b", "a*b*cp"), c("a", "a*b", "a*cp"),
                 c("a^2", "a*cp"))) {
    expected <- .mbco_exact_redundant(h)
    expect_true(length(expected) > 0L)
    err <- tryCatch(
      mediation_mbco(model, data = d, hypotheses = list(joint = h),
                     ci_method = "wald"),
      error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "hypotheses")
    expect_identical(err$found_by, "null_set")
    expect_identical(err$redundant, h[expected])
    expect_identical(err$kept, h[-expected])
    expect_identical(err$n_constraints, length(h))
    expect_match(conditionMessage(err), "Nothing was fitted.", fixed = TRUE)
  }
})

test_that("the null-set check agrees with the exact algebra of products (DR38)", {
  # Random sets of two or three products of labels, with constant
  # multiples and squares, at parameter values from 1e-6 to 1e3 in size
  # (badly scaled units included): the refusal decision of the two checks
  # mediation_mbco() runs before fitting, and the constraints the null-set
  # check names, agree with the exact algebra.
  withr::local_seed(3804)
  labs <- c("a", "b", "cp", "d")
  for (i in seq_len(150)) {
    m <- sample(2:3, 1)
    h <- vapply(seq_len(m), function(j) {
      f <- sample(labs, sample(1:3, 1))
      paste0(if (stats::runif(1) < .3) "-2*" else "",
             paste0(f, ifelse(stats::runif(length(f)) < .15, "^2", ""),
                    collapse = "*"))
    }, character(1))
    theta <- stats::setNames(sample(c(-1, 1), 4, TRUE) *
                               10^stats::runif(4, -6, 3), labs)
    expected <- .mbco_exact_redundant(h)
    red <- .mbco_redundant(h, theta)
    if (length(red$redundant) == 0L) {
      red <- .mbco_vanishing(h, theta)
      expect_identical(red$redundant, expected, info = paste(h, collapse = ", "))
    } else {
      # The Jacobian check names a different member of a dependent set (the
      # one that does not raise the rank of those before it); the decision
      # is the same.
      expect_true(length(expected) > 0L, info = paste(h, collapse = ", "))
    }
  }
})

test_that("a joint test whose branches differ in restrictions is kept, as the oracle (DR38)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  # c("a*b", "a*cp") has the null set a = 0 (one restriction) or b = cp = 0
  # (two); no constraint is 0 wherever the other is, so it is tested with
  # 2 df, and its statistic is the smaller of the two branch statistics.
  # Least squares oracle: LR(a = 0) and LR(b = cp = 0) = N log(SSE(y ~ 1) /
  # SSE(y ~ m + x)).
  d <- .mbco_dr38_data()
  N <- nrow(d)
  res <- mediation_mbco("m ~ a*x\ny ~ b*m + cp*x", data = d,
                        hypotheses = list(mixed = c("a*b", "a*cp")),
                        ci_method = "wald")
  sse <- function(X, v) sum(stats::lm.fit(X, v)$residuals^2)
  one <- matrix(1, N, 1L)
  lr_bcp <- N * log(sse(one, d$y) / sse(cbind(one, d$m, d$x), d$y))
  orc <- min(.mbco_ls_branches(d)[["a"]], lr_bcp)
  row <- res[res$term == "mixed", ]
  expect_identical(row$df, 2L)
  expect_equal(row$lrt, orc, tolerance = 1e-4)
  expect_equal(row$p_value, stats::pchisq(orc, 2, lower.tail = FALSE),
               tolerance = 1e-4)
})

test_that("a redundant constraint in a sequential model is refused (DR38)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  # Review probe: c("a*d*b", "a*d") has the null set ad = 0, that of "a*d"
  # alone; it gave LRT 8.708 with 2 df, p = .0129, where "a*d" gives .0032.
  set.seed(2611)
  n <- 150
  x <- stats::rnorm(n)
  m1 <- .25 * x + stats::rnorm(n)
  m2 <- .3 * m1 + stats::rnorm(n)
  y <- .3 * m2 + .1 * x + stats::rnorm(n)
  ds <- data.frame(x, m1, m2, y)
  err <- tryCatch(
    mediation_mbco("m1 ~ a*x\nm2 ~ d*m1 + e*x\ny ~ b*m2 + f*m1 + cp*x",
                   data = ds, x = "x", y = "y", ci_method = "wald",
                   hypotheses = list(joint = c("a*d*b", "a*d"))),
    error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$redundant, "a*d*b")
  expect_identical(err$kept, "a*d")
})

test_that("a common label is taken out of a sum of products exactly (review)", {
  # The branch fits of a conditional effect rest on this factoring; the
  # product of the factors must equal the sum at any parameter values.
  withr::local_seed(3805)
  for (e in c("(a0*b) + (a1*b)*(10)", "a1*b1 + a1*d*b2", "-(a*b) + 2*a*cp",
              "a*b - (a*c - a*d)", "(2*a*b*c) + (a*b)*(-3)")) {
    f <- .mbco_common_factors(parse(text = e)[[1L]])
    expect_gt(length(f), 1L)
    for (i in 1:5) {
      env <- as.list(stats::setNames(stats::rnorm(8),
                                     c("a", "a0", "a1", "b", "b1", "b2",
                                       "c", "cp")))
      env$d <- stats::rnorm(1)
      lhs <- eval(parse(text = e)[[1L]], env)
      rhs <- prod(vapply(f, function(x) eval(parse(text = x)[[1L]], env),
                         numeric(1)))
      expect_equal(rhs, lhs, tolerance = 1e-12, info = e)
    }
  }
  # No label common to every term, or a term that is not a product: the
  # sum is one factor.
  expect_identical(.mbco_common_factors(quote(a*b + cp)), "a * b + cp")
  expect_length(.mbco_common_factors(quote(a*b + a*exp(b))), 1L)
})

test_that("the statistics do not depend on the units of the variables (review)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  # Through the first builds of 1.1.0, with x in thousands and m in
  # thousandths, "a" gave LRT 1316.2, "b" 96.49, c("a", "b") 21.98, and the
  # indirect effect 7.56. The oracle is least squares on the same data: the
  # likelihood factors into x, m | x, y | m, x, and the total effect
  # constraint cp + ab = 0 is cov(x, y) = 0, a restriction of the (x, y)
  # margin, whose statistic is -N log(1 - r_xy^2).
  d <- .mbco_dr38_data()
  N <- nrow(d)
  d2 <- d
  d2$x <- d$x * 1000
  d2$m <- d$m / 1000
  # lavaan, whose fit gives the starting values, notes that it rescaled the
  # badly scaled variables itself; the note is a message, set aside here.
  res <- suppressMessages(mediation_mbco(
    "m ~ a*x\ny ~ b*m + cp*x", data = d2,
    hypotheses = list(ha = "a", hb = "b", pair = c("a", "b"),
                      rotated = c("a + b", "a - b")),
    ci_method = "wald"))
  orc <- .mbco_ls_branches(d)
  lrt <- stats::setNames(res$lrt, res$term)
  expect_equal(lrt[["ha"]], orc[["a"]], tolerance = 1e-5)
  expect_equal(lrt[["hb"]], orc[["b"]], tolerance = 1e-5)
  expect_equal(lrt[["pair"]], sum(orc), tolerance = 1e-5)
  expect_equal(lrt[["rotated"]], sum(orc), tolerance = 1e-5)
  expect_equal(lrt[["indirect_via_m"]], min(orc), tolerance = 1e-5)
  expect_equal(lrt[["total_effect"]], -N * log(1 - stats::cor(d$x, d$y)^2),
               tolerance = 1e-5)
  # The deviance is that of the variables as given: the unscaled deviance
  # plus 2 N log(1000) - 2 N log(1000) (x times 1000, m divided by 1000).
  res0 <- mediation_mbco("m ~ a*x\ny ~ b*m + cp*x", data = d,
                         ci_method = "wald")
  expect_equal(attr(res, "deviance"), attr(res0, "deviance"),
               tolerance = 1e-8)
  expect_false(is.null(attr(res, "rescale")))
  # Estimates and standard errors are in the units of the data: the
  # indirect effect of x on y is divided by 1000.
  expect_equal(res$estimate[res$term == "indirect_via_m"] * 1000,
               res0$estimate[res0$term == "indirect_via_m"], tolerance = 1e-6)
  expect_equal(res$se[res$term == "indirect_via_m"] * 1000,
               res0$se[res0$term == "indirect_via_m"], tolerance = 1e-4)
})

test_that("a conditional indirect effect is tested on its best branch, as least squares (review)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  # First-stage moderated mediation: the conditional indirect effect at a
  # probe value v is (a0 + a1 v) b, zero on the branch b = 0 or on the
  # branch a0 + a1 v = 0. Through the first builds of 1.1.0 only products
  # of labels had branch fits, and at w = 1 the constrained runs missed
  # b = 0: 37.294 where the smaller branch statistic is 33.631. Least
  # squares oracle: LR(b = 0) = N log(SSE(y ~ x) / SSE(y ~ m + x)),
  # LR(a0 + a1 v = 0) = N log(SSE(m ~ w + I(xw - v x)) / SSE(m ~ x + w + xw)).
  set.seed(3812)
  n <- 300
  w <- stats::rnorm(n)
  x <- stats::rnorm(n)
  m <- .4 * x + .2 * x * w + stats::rnorm(n)
  y <- .3 * m + .1 * x + stats::rnorm(n)
  sse <- function(f) sum(stats::resid(f)^2)
  xw <- x * w
  lr_b <- n * log(sse(stats::lm(y ~ x)) / sse(stats::lm(y ~ m + x)))
  full <- sse(stats::lm(m ~ x + w + xw))
  lr_av <- function(v) n * log(sse(stats::lm(m ~ w + I(xw - v * x))) / full)
  orc <- c(min(lr_b, lr_av(-1)), min(lr_b, lr_av(0)), min(lr_b, lr_av(1)))
  for (k in c(1, 100)) {
    res <- mediation_mbco("m ~ x + w + x:w\ny ~ m + x",
                          data = data.frame(x, w = w * k, m, y),
                          x = "x", y = "y", moderator = "w",
                          probe_values = c(-k, 0, k), ci_method = "wald")
    got <- res$lrt[grepl("^indirect_via_m_at", res$term)]
    expect_equal(got, orc, tolerance = 5e-3, info = paste("w times", k))
  }
})

test_that("the indirect effect statistic is the smaller branch statistic (DR44)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  # With a = b = 0 the null set of ab = 0 is the crossing of the branches
  # a = 0 and b = 0, the best-fitting null model lies on one of them, and
  # the statistic is the smaller of the two branch statistics, which is
  # asymptotically min(Za^2, Zb^2) and not chi square on 1 df (see the help
  # page). The oracle is the least squares form of each branch statistic.
  for (seed in c(44101, 44102, 44103)) {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(200), m = stats::rnorm(200),
                    y = stats::rnorm(200))
    res <- mediation_mbco("m ~ a*x\ny ~ b*m + cp*x", data = d,
                          hypotheses = c(indirect = "a*b"),
                          ci_method = "wald")
    orc <- .mbco_ls_branches(d)
    expect_equal(res$lrt[res$term == "indirect"], min(orc), tolerance = 1e-4)
    expect_identical(res$df[res$term == "indirect"], 1L)
  }
})

test_that("at a = b = 0 the test rejects near .0025 in the design of Tofighi and Kelley (2020)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  # Ruling 18 of 2026-11-11: the help page states the large sample rate
  # .05^2 = .0025 at the crossing and the published simulated rates (.025
  # to .075, about .03 to .04 in their Figure 3, at N = 50 to 200), and
  # that this function's statistic rejects near the limit at those sizes.
  # The design is theirs: X binary with probability .5, M and Y standard
  # normal, a = b = 0, the direct path fitted. First, the statistic the
  # function returns is the smaller branch statistic for this design too
  # (least squares oracle, .mbco_ls_branches() above).
  set.seed(18508)
  for (i in 1:2) {
    d <- data.frame(x = stats::rbinom(50, 1, .5), m = stats::rnorm(50),
                    y = stats::rnorm(50))
    res <- mediation_mbco("m ~ a*x\ny ~ b*m + cp*x", data = d, x = "x",
                          y = "y", ci_method = "wald")
    expect_equal(res$lrt[res$term == "indirect_via_m"],
                 min(.mbco_ls_branches(d)), tolerance = 1e-5)
  }
  # Second, the rate of that statistic at N = 50, from the least squares
  # form of the two branch statistics in 40,000 data sets (the help page
  # reports .0032 from 400,000, Monte Carlo SE .00009). The SE here is
  # about .00028; the rate is far below the published .025 to .075 and
  # within 5 SEs of .0032.
  set.seed(18509)
  N <- 50; R <- 40000
  X <- matrix(stats::rbinom(N * R, 1, .5), N, R)
  M <- matrix(stats::rnorm(N * R), N, R)
  Y <- matrix(stats::rnorm(N * R), N, R)
  n1 <- colSums(X)
  keep <- n1 > 0 & n1 < N
  within <- function(V) {
    m1 <- colSums(V * X) / n1
    m0 <- colSums(V * (1 - X)) / (N - n1)
    V - (X * rep(m1, each = N) + (1 - X) * rep(m0, each = N))
  }
  rm <- within(M); ry <- within(Y)
  lrt_a <- N * log(colSums(sweep(M, 2, colMeans(M))^2) / colSums(rm^2))
  sse_y <- colSums(ry^2)
  lrt_b <- N * log(sse_y / (sse_y - colSums(ry * rm)^2 / colSums(rm^2)))
  rate <- mean((pmin(lrt_a, lrt_b) > stats::qchisq(.95, 1))[keep])
  expect_identical(sum(!keep), 0L)
  expect_lt(abs(rate - .0032), 5 * sqrt(.0032 * (1 - .0032) / R))
  expect_lt(rate, .025)
})
