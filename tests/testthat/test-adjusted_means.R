test_that("adjusted_means() returns the documented wide schema", {
  set.seed(113)
  d <- data.frame(
    A = factor(rep(c("a1", "a2"), times = c(30, 50))),
    B = factor(c(rep(c("b1", "b2"), times = c(10, 20)),
                 rep(c("b1", "b2"), times = c(35, 15)))),
    x = rnorm(80)
  )
  d$y <- 5 + 2 * (d$A == "a2") + 1.5 * (d$B == "b2") + 0.8 * d$x + rnorm(80)
  fit <- lm(y ~ A * B + x, data = d)

  out <- adjusted_means(fit)
  expect_s3_class(out, "dmar_tbl")
  expect_identical(names(out), c("A", "B", "estimate", "se",
                                 "ci_lower", "ci_upper"))
  # One row per cell of the reference grid, first factor varying fastest.
  expect_identical(out$A, c("a1", "a2", "a1", "a2"))
  expect_identical(out$B, c("b1", "b1", "b2", "b2"))
  expect_true(all(vapply(out[c("estimate", "se", "ci_lower", "ci_upper")],
                         is.numeric, logical(1L))))
  expect_identical(attr(out, "conf_level"), 0.95)
  expect_identical(attr(out, "df_residual"), df.residual(fit))
  expect_null(attr(out, "weights"))
  expect_true(all(out$ci_lower < out$estimate & out$estimate < out$ci_upper))

  marg <- adjusted_means(fit, by = "A")
  expect_identical(names(marg), c("A", "estimate", "se",
                                  "ci_lower", "ci_upper"))
  expect_identical(marg$A, c("a1", "a2"))
  expect_identical(attr(marg, "weights"), "equal")

  # Without by there is no averaging, so the weights argument changes nothing.
  expect_equal(adjusted_means(fit, weights = "proportional"), out)

  # The default wide-table broom verbs apply.
  td <- generics::tidy(out)
  expect_s3_class(td, "data.frame")
  expect_identical(nrow(td), nrow(out))
})

test_that("one-way cell means match pinned emmeans values", {
  fit <- lm(weight ~ group, data = PlantGrowth)
  out <- adjusted_means(fit)

  # Pinned from emmeans::emmeans() (emmeans 2.0.3, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  expect_identical(out$group, c("ctrl", "trt1", "trt2"))
  expect_equal(out$estimate, c(5.032, 4.661, 5.526), tolerance = 1e-8)
  expect_equal(out$se, rep(0.197128365773657, 3), tolerance = 1e-8)
  expect_equal(out$ci_lower,
               c(4.62752600344172, 4.25652600344172, 5.12152600344172),
               tolerance = 1e-8)
  expect_equal(out$ci_upper,
               c(5.43647399655828, 5.06547399655828, 5.93047399655828),
               tolerance = 1e-8)
  expect_identical(attr(out, "df_residual"), 27L)

  # Balanced data: the two weightings coincide, and the marginal means of the
  # only factor are the cell means themselves.
  eq <- adjusted_means(fit, by = "group")
  pr <- adjusted_means(fit, by = "group", weights = "proportional")
  expect_equal(eq$estimate, out$estimate, tolerance = 1e-12)
  expect_equal(pr$estimate, out$estimate, tolerance = 1e-12)
  expect_equal(pr$se, eq$se, tolerance = 1e-12)
})

test_that("unbalanced two-way marginal means match pinned emmeans values for both weightings", {
  set.seed(407)
  n_cell <- c(7, 12, 9, 15)  # a1b1, a2b1, a1b2, a2b2
  d <- do.call(rbind, lapply(1:4, function(i) {
    g <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2"))
    g[rep(i, n_cell[i]), ]
  }))
  d$A <- factor(d$A); d$B <- factor(d$B)
  d$y <- 2 + (d$A == "a2") + 2 * (d$B == "b2") +
    1.5 * (d$A == "a2") * (d$B == "b2") + rnorm(nrow(d))
  fit <- lm(y ~ A * B, data = d)

  # With no covariate the cell means are the observed cell means.
  cells <- adjusted_means(fit)
  obs <- as.numeric(t(tapply(d$y, list(d$B, d$A), mean)))
  expect_equal(cells$estimate, obs, tolerance = 1e-10)

  # Pinned from emmeans::emmeans(fit, ~ A, weights = ) (emmeans 2.0.3,
  # 2026-08-09); live comparison in tools/oracle_checks.R.
  eq <- adjusted_means(fit, by = "A")
  expect_equal(eq$estimate, c(2.83453697890022, 4.51887731317308),
               tolerance = 1e-8)
  expect_equal(eq$se, c(0.272874275362897, 0.209709694870538),
               tolerance = 1e-8)
  expect_equal(eq$ci_lower, c(2.28259665981208, 4.09469941751476),
               tolerance = 1e-8)
  expect_equal(eq$ci_upper, c(3.38647729798836, 4.9430552088314),
               tolerance = 1e-8)

  pr <- adjusted_means(fit, by = "A", weights = "proportional")
  expect_identical(attr(pr, "weights"), "proportional")
  expect_equal(pr$estimate, c(2.93374252479812, 4.74850523919718),
               tolerance = 1e-8)
  expect_equal(pr$se, c(0.27074451042749, 0.208413988214084),
               tolerance = 1e-8)
  expect_equal(pr$ci_lower, c(2.38611006190664, 4.3269481576279),
               tolerance = 1e-8)
  expect_equal(pr$ci_upper, c(3.48137498768961, 5.17006232076646),
               tolerance = 1e-8)

  # In unbalanced data the two weightings genuinely differ.
  expect_gt(max(abs(eq$estimate - pr$estimate)), 0.05)
})

test_that("three-way proportional averaging matches pinned emmeans values", {
  set.seed(431)
  n_cell <- c(4, 7, 5, 9, 6, 3, 8, 5, 7, 4, 6, 8)
  d <- do.call(rbind, lapply(seq_along(n_cell), function(i) {
    g <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2", "b3"),
                     C = c("c1", "c2"))
    g[rep(i, n_cell[i]), ]
  }))
  d$A <- factor(d$A); d$B <- factor(d$B); d$C <- factor(d$C)
  d$y <- 1 + as.numeric(d$A) + 0.5 * as.numeric(d$B) +
    0.3 * as.numeric(d$C) + rnorm(nrow(d))
  fit <- lm(y ~ A * B * C, data = d)

  # Averaging over A and B weights each (A, B) combination by its joint
  # frequency in the data. Pinned from emmeans::emmeans(fit, ~ C,
  # weights = "proportional") (emmeans 2.0.3, 2026-08-09); live comparison
  # in tools/oracle_checks.R.
  pr <- adjusted_means(fit, by = "C", weights = "proportional")
  expect_identical(pr$C, c("c1", "c2"))
  expect_equal(pr$estimate, c(3.43509820943712, 4.19380393870638),
               tolerance = 1e-8)
  expect_equal(pr$se, c(0.198907010552974, 0.186894238613107),
               tolerance = 1e-8)
  expect_equal(pr$ci_lower, c(3.03722494944464, 3.81995980026157),
               tolerance = 1e-8)
  expect_equal(pr$ci_upper, c(3.8329714694296, 4.5676480771512),
               tolerance = 1e-8)

  # A two-factor margin works the same way and keeps the requested order,
  # first named factor varying fastest.
  m2 <- adjusted_means(fit, by = c("B", "A"))
  expect_identical(names(m2)[1:2], c("B", "A"))
  expect_identical(m2$B, rep(c("b1", "b2", "b3"), 2))
  expect_identical(m2$A, rep(c("a1", "a2"), each = 3))
})

test_that("a one-way ANCOVA matches pinned emmeans values and ancova()", {
  set.seed(113)
  d <- data.frame(g = factor(rep(c("g1", "g2", "g3"),
                                 times = c(12, 18, 15))))
  d$x <- rnorm(45, 10, 2)
  d$y <- 3 + 2 * (d$g == "g2") + 4 * (d$g == "g3") + 0.6 * d$x + rnorm(45)
  fit <- lm(y ~ g + x, data = d)
  out <- adjusted_means(fit)

  # Pinned from emmeans::emmeans(fit, ~ g) (emmeans 2.0.3, 2026-08-09);
  # live comparison in tools/oracle_checks.R.
  expect_equal(out$estimate,
               c(9.50834662038919, 10.9218615796478, 13.3882062892459),
               tolerance = 1e-8)
  expect_equal(out$se,
               c(0.335700089962901, 0.273036903087181, 0.29587928175921),
               tolerance = 1e-8)
  expect_equal(out$ci_lower,
               c(8.83038653492825, 10.3704523674208, 12.7906659574284),
               tolerance = 1e-8)
  expect_equal(out$ci_upper,
               c(10.1863067058501, 11.4732707918748, 13.9857466210634),
               tolerance = 1e-8)
  expect_identical(attr(out, "df_residual"), 41L)

  # The same adjusted means and standard errors that ancova() reports.
  av <- ancova(data = d, outcome = "y", treatment = "g", covariates = "x")
  expect_equal(av$value[grep("^adjusted_mean\\[", av$term)], out$estimate,
               tolerance = 1e-10)
  expect_equal(av$value[grep("^se_adjusted_mean\\[", av$term)], out$se,
               tolerance = 1e-10)
})

test_that("the cell table and contrast_adjusted() agree", {
  set.seed(509)
  d <- data.frame(
    A = factor(rep(c("a1", "a2"), times = c(28, 36))),
    B = factor(rep(rep(c("b1", "b2"), times = c(9, 7)), 4)),
    x1 = rnorm(64), x2 = runif(64, -2, 2)
  )
  d$y <- 1 + (d$A == "a2") + 0.5 * (d$B == "b2") + 0.7 * d$x1 -
    0.4 * d$x2 + rnorm(64)
  fit <- lm(y ~ A * B + x1 + x2, data = d)
  cells <- adjusted_means(fit)

  # A unit contrast picks out one adjusted cell mean, so contrast_adjusted()
  # must reproduce that row's estimate and standard error exactly (its se is
  # recovered as estimate / t).
  for (j in seq_len(nrow(cells))) {
    ca <- contrast_adjusted(fit,
                            contrast = as.numeric(seq_len(nrow(cells)) == j))
    expect_equal(ca$value[ca$term == "contrast"], cells$estimate[j],
                 tolerance = 1e-10)
    expect_equal(ca$value[ca$term == "contrast"] / ca$value[ca$term == "t"],
                 cells$se[j], tolerance = 1e-10)
  }
  # And a genuine contrast of the cells is the same weighted sum of the
  # table's estimates.
  w <- c(-0.5, 0.5, -0.5, 0.5)
  ca <- contrast_adjusted(fit, contrast = w)
  expect_equal(ca$value[ca$term == "contrast"], sum(w * cells$estimate),
               tolerance = 1e-10)

  # Pinned from emmeans::emmeans(fit, ~ B, weights = "proportional")
  # (emmeans 2.0.3, 2026-08-09); live comparison in tools/oracle_checks.R.
  pr <- adjusted_means(fit, by = "B", weights = "proportional")
  expect_equal(pr$estimate, c(1.30638243555216, 2.06200767538313),
               tolerance = 1e-8)
  expect_equal(pr$se, c(0.15243236286902, 0.174274041142233),
               tolerance = 1e-8)
})

test_that("a covariate-by-factor interaction adjusts at the covariate mean", {
  set.seed(617)
  d <- data.frame(A = factor(rep(c("lo", "mid", "hi"),
                                 times = c(14, 11, 17)),
                             levels = c("lo", "mid", "hi")))
  d$x <- rnorm(42, 5, 1.5)
  d$y <- 2 + (d$A == "mid") + 3 * (d$A == "hi") +
    (0.5 + 0.4 * (d$A == "hi")) * d$x + rnorm(42)
  fit <- lm(y ~ A * x, data = d)
  out <- adjusted_means(fit)

  # Pinned from emmeans::emmeans(fit, ~ A) (emmeans 2.0.3, 2026-08-09);
  # live comparison in tools/oracle_checks.R.
  expect_identical(out$A, c("lo", "mid", "hi"))
  expect_equal(out$estimate,
               c(4.82201472270881, 5.43542523896586, 9.0712075924237),
               tolerance = 1e-8)
  expect_equal(out$se,
               c(0.315543488172012, 0.411183395716631, 0.25873013864845),
               tolerance = 1e-8)
  expect_equal(out$ci_lower,
               c(4.18206286729871, 4.60150666081019, 8.54647855035794),
               tolerance = 1e-8)
  expect_equal(out$ci_upper,
               c(5.46196657811892, 6.26934381712153, 9.59593663448946),
               tolerance = 1e-8)
})

test_that("transformed covariates are evaluated at the raw covariate mean", {
  # poly(x, 2): the basis is evaluated at mean(x) with the coefficients
  # stored at fit time. Pinned from emmeans::emmeans(fit, ~ A)
  # (emmeans 2.0.3, 2026-08-09); live comparison in tools/oracle_checks.R.
  set.seed(701)
  d1 <- data.frame(A = factor(rep(c("p", "q"), each = 24)))
  d1$x <- runif(48, 0, 4)
  d1$y <- 1 + 2 * (d1$A == "q") + 0.8 * d1$x - 0.15 * d1$x^2 +
    rnorm(48, 0, 0.7)
  fit1 <- lm(y ~ A + poly(x, 2), data = d1)
  out1 <- adjusted_means(fit1)
  expect_equal(out1$estimate, c(1.90635642250777, 3.98251037096841),
               tolerance = 1e-8)
  expect_equal(out1$se, c(0.170778740501233, 0.220123421059232),
               tolerance = 1e-8)
  expect_equal(out1$ci_lower, c(1.56217448649724, 3.53888076579),
               tolerance = 1e-8)
  expect_equal(out1$ci_upper, c(2.2505383585183, 4.42613997614681),
               tolerance = 1e-8)

  # log(x): the grid carries mean(x), so the cell mean is the prediction at
  # the raw covariate mean, exactly what predict() returns there. Pinned
  # from emmeans::emmeans(fit, ~ A) (emmeans 2.0.3, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  set.seed(809)
  d2 <- data.frame(A = factor(rep(c("u", "v"), times = c(21, 17))))
  d2$x <- rexp(38) + 0.5
  d2$y <- 2 + (d2$A == "v") + 1.4 * log(d2$x) + rnorm(38, 0, 0.6)
  fit2 <- lm(y ~ A + log(x), data = d2)
  out2 <- adjusted_means(fit2)
  expect_equal(out2$estimate, c(2.64819433120467, 3.78648250403861),
               tolerance = 1e-8)
  expect_equal(out2$se, c(0.125632923818136, 0.131938528213672),
               tolerance = 1e-8)
  expect_equal(out2$ci_lower, c(2.3931459365122, 3.51863305187035),
               tolerance = 1e-8)
  expect_equal(out2$ci_upper, c(2.90324272589714, 4.05433195620686),
               tolerance = 1e-8)
  pred <- predict(fit2, newdata = data.frame(
    A = factor(c("u", "v"), levels = c("u", "v")), x = mean(d2$x)),
    se.fit = TRUE)
  expect_equal(out2$estimate, unname(pred$fit), tolerance = 1e-10)
  expect_equal(out2$se, unname(pred$se.fit), tolerance = 1e-10)
})

test_that("a character predictor uses the fit's sorted levels", {
  set.seed(907)
  d <- data.frame(g = rep(c("zeta", "alpha", "mid"), times = c(13, 16, 11)),
                  stringsAsFactors = FALSE)
  d$x <- rnorm(40)
  d$y <- 1 + (d$g == "mid") + 2 * (d$g == "zeta") + 0.5 * d$x + rnorm(40)
  fit <- lm(y ~ g + x, data = d)
  out <- adjusted_means(fit)

  # Pinned from emmeans::emmeans(fit, ~ g) (emmeans 2.0.3, 2026-08-09);
  # live comparison in tools/oracle_checks.R.
  expect_identical(out$g, c("alpha", "mid", "zeta"))
  expect_equal(out$estimate,
               c(1.08986150477403, 2.22465303326002, 2.6194456292884),
               tolerance = 1e-8)
  expect_equal(out$se,
               c(0.251513169304767, 0.307947967578857, 0.280472857364745),
               tolerance = 1e-8)
})

test_that("an ordered factor keeps its polynomial contrasts", {
  set.seed(1013)
  d <- data.frame(
    dose = factor(rep(c("low", "med", "high"), times = c(15, 12, 18)),
                  levels = c("low", "med", "high"), ordered = TRUE),
    B = factor(rep(c("b1", "b2"), length.out = 45))
  )
  d$x <- rnorm(45)
  d$y <- 1 + 0.8 * as.numeric(d$dose) + 0.5 * (d$B == "b2") +
    0.3 * d$x + rnorm(45)
  fit <- lm(y ~ dose * B + x, data = d)

  cells <- adjusted_means(fit)
  expect_identical(nrow(cells), 6L)

  # Pinned from emmeans::emmeans(fit, ~ dose, weights = "proportional")
  # (emmeans 2.0.3, 2026-08-09); live comparison in tools/oracle_checks.R.
  pr <- adjusted_means(fit, by = "dose", weights = "proportional")
  expect_identical(pr$dose, c("low", "med", "high"))
  expect_equal(pr$estimate,
               c(2.06839076702721, 2.74939521846345, 3.3487910138557),
               tolerance = 1e-8)
  expect_equal(pr$se,
               c(0.237794138008398, 0.263441212529475, 0.217640205280305),
               tolerance = 1e-8)
  expect_equal(pr$ci_lower,
               c(1.58700170183053, 2.21608636528488, 2.90820145245364),
               tolerance = 1e-8)
  expect_equal(pr$ci_upper,
               c(2.54977983222389, 3.28270407164201, 3.78938057525775),
               tolerance = 1e-8)
})

test_that("a weighted fit averages by total prior weight", {
  set.seed(407)
  n_cell <- c(7, 12, 9, 15)
  d <- do.call(rbind, lapply(1:4, function(i) {
    g <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2"))
    g[rep(i, n_cell[i]), ]
  }))
  d$A <- factor(d$A); d$B <- factor(d$B)
  d$y <- 2 + (d$A == "a2") + 2 * (d$B == "b2") +
    1.5 * (d$A == "a2") * (d$B == "b2") + rnorm(nrow(d))
  set.seed(1201)
  d$w <- runif(nrow(d), 0.5, 2)
  fit <- lm(y ~ A * B, data = d, weights = w)

  # Pinned from emmeans::emmeans(fit, ~ A, weights = "proportional")
  # (emmeans 2.0.3, 2026-08-09); live comparison in tools/oracle_checks.R.
  pr <- adjusted_means(fit, by = "A", weights = "proportional")
  expect_equal(pr$estimate, c(2.92170737045355, 4.63703805293535),
               tolerance = 1e-8)
  expect_equal(pr$se, c(0.26703789214755, 0.19791956370645),
               tolerance = 1e-8)
  expect_equal(pr$ci_lower, c(2.38157225070094, 4.23670794852868),
               tolerance = 1e-8)
  expect_equal(pr$ci_upper, c(3.46184249020615, 5.03736815734203),
               tolerance = 1e-8)
})

test_that("aov fits and other confidence levels work", {
  set.seed(509)
  d <- data.frame(
    A = factor(rep(c("a1", "a2"), times = c(28, 36))),
    B = factor(rep(rep(c("b1", "b2"), times = c(9, 7)), 4)),
    x1 = rnorm(64), x2 = runif(64, -2, 2)
  )
  d$y <- 1 + (d$A == "a2") + 0.5 * (d$B == "b2") + 0.7 * d$x1 -
    0.4 * d$x2 + rnorm(64)
  fit <- aov(y ~ A * B + x1 + x2, data = d)

  # Pinned from emmeans::emmeans(fit, ~ A, options = list(level = 0.90))
  # (emmeans 2.0.3, 2026-08-09); live comparison in tools/oracle_checks.R.
  out <- adjusted_means(fit, by = "A", conf_level = 0.90)
  expect_identical(attr(out, "conf_level"), 0.90)
  expect_equal(out$estimate, c(1.34205401612113, 1.95030475273715),
               tolerance = 1e-8)
  expect_equal(out$se, c(0.182163275194002, 0.152646287370058),
               tolerance = 1e-8)
  expect_equal(out$ci_lower, c(1.03755849025277, 1.69514842940525),
               tolerance = 1e-8)
  expect_equal(out$ci_upper, c(1.64654954198949, 2.20546107606905),
               tolerance = 1e-8)
})

test_that("cell-means and crossed parameterizations agree", {
  set.seed(407)
  n_cell <- c(7, 12, 9, 15)
  d <- do.call(rbind, lapply(1:4, function(i) {
    g <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2"))
    g[rep(i, n_cell[i]), ]
  }))
  d$A <- factor(d$A); d$B <- factor(d$B)
  d$y <- 2 + (d$A == "a2") + 2 * (d$B == "b2") +
    1.5 * (d$A == "a2") * (d$B == "b2") + rnorm(nrow(d))
  d$cell <- interaction(d$A, d$B)  # levels: a1.b1, a2.b1, a1.b2, a2.b2

  am_crossed <- adjusted_means(lm(y ~ A * B, data = d))
  am_cells   <- adjusted_means(lm(y ~ 0 + cell, data = d))
  expect_equal(am_cells$estimate, am_crossed$estimate, tolerance = 1e-10)
  expect_equal(am_cells$se, am_crossed$se, tolerance = 1e-10)
})

test_that("rank-deficient designs are refused with the offending rows named", {
  # A 2 x 2 ANCOVA with the (a2, b2) cell empty aliases the interaction
  # coefficient; emmeans reports the same rows as nonEst. The affected rows
  # must be named, and only those rows.
  set.seed(1103)
  d <- do.call(rbind,
               lapply(list(c("a1", "b1"), c("a2", "b1"), c("a1", "b2")),
                      function(cell) {
                        data.frame(A = cell[1], B = cell[2])[rep(1, 20), ]
                      }))
  d$A <- factor(d$A, levels = c("a1", "a2"))
  d$B <- factor(d$B, levels = c("b1", "b2"))
  d$x <- rnorm(60)
  d$y <- 5 + 2 * (d$A == "a2") + 1.5 * (d$B == "b2") + 0.8 * d$x + rnorm(60)
  fit <- lm(y ~ A * B + x, data = d)
  expect_true(anyNA(coef(fit)))  # the empty cell aliased a coefficient

  expect_error(adjusted_means(fit), "not estimable")
  msg <- tryCatch(adjusted_means(fit), error = conditionMessage)
  expect_match(msg, "A = a2, B = b2", fixed = TRUE)
  expect_false(grepl("a1", msg))
  expect_false(grepl("A = a2, B = b1", msg, fixed = TRUE))

  # A marginal mean that averages over the empty cell is refused too, and
  # only for the level whose average needs that cell.
  msg_by <- tryCatch(adjusted_means(fit, by = "A"), error = conditionMessage)
  expect_match(msg_by, "not estimable")
  expect_match(msg_by, "A = a2", fixed = TRUE)
  expect_false(grepl("a1", msg_by))
})

test_that("rows dropped by na.action or subset= do not shift the covariate mean", {
  set.seed(29)
  d <- data.frame(g = factor(rep(c("g1", "g2", "g3"), 14)), x = rnorm(42))
  d$y <- as.numeric(d$g) + 0.7 * d$x + rnorm(42)
  d$x[c(4, 9)] <- NA
  fit <- lm(y ~ g + log(x + 10), data = d, subset = -3)

  out <- adjusted_means(fit)
  used <- as.integer(rownames(model.frame(fit)))
  pred <- predict(fit, newdata = data.frame(
    g = factor(levels(d$g), levels = levels(d$g)),
    x = mean(d$x[used])), se.fit = TRUE)
  expect_equal(out$estimate, unname(pred$fit), tolerance = 1e-10)
  expect_equal(out$se, unname(pred$se.fit), tolerance = 1e-10)
})

test_that("logical predictors define cells like factors", {
  set.seed(23)
  d <- data.frame(trt = rep(c(TRUE, FALSE), times = c(18, 22)),
                  x = rnorm(40))
  d$y <- 1 + 2 * d$trt + 0.5 * d$x + rnorm(40)
  fit <- lm(y ~ trt + x, data = d)
  out <- adjusted_means(fit)
  expect_identical(out$trt, c("FALSE", "TRUE"))
  pred <- predict(fit, newdata = data.frame(trt = c(FALSE, TRUE),
                                            x = mean(d$x)), se.fit = TRUE)
  expect_equal(out$estimate, unname(pred$fit), tolerance = 1e-10)
  expect_equal(out$se, unname(pred$se.fit), tolerance = 1e-10)
})

test_that("input validation errors are informative", {
  set.seed(113)
  d <- data.frame(g = factor(rep(c("g1", "g2"), each = 20)), x = rnorm(40))
  d$y <- d$x + as.numeric(d$g) + rnorm(40)
  fit <- lm(y ~ g + x, data = d)

  # Within-subjects (multi-stratum) fits are out of scope.
  dw <- data.frame(subj = factor(rep(1:8, each = 3)),
                   cond = factor(rep(c("c1", "c2", "c3"), 8)),
                   y = rnorm(24))
  fit_ws <- aov(y ~ cond + Error(subj), data = dw)
  expect_error(adjusted_means(fit_ws), "multi-stratum")
  expect_error(adjusted_means(fit_ws), "Error\\(\\) term")

  # Other unsupported fits.
  expect_error(adjusted_means(glm(I(y > 3) ~ g + x, family = binomial,
                                  data = d)),
               "must be a fitted lm or aov")
  expect_error(adjusted_means(lm(cbind(mpg, hp) ~ wt, data = mtcars)),
               "multivariate")
  expect_error(adjusted_means(lm(mpg ~ wt, data = mtcars)),
               "at least one factor")
  expect_error(adjusted_means(lm(mpg ~ factor(cyl), data = mtcars)),
               "created inside the model formula")
  expect_error(adjusted_means(lm(y ~ g + x + offset(x), data = d)),
               "offset")

  # Argument validation.
  expect_error(adjusted_means(fit, conf_level = 2),
               "strictly between 0 and 1")
  expect_error(adjusted_means(fit, by = 1),
               "character vector of factor names")
  expect_error(adjusted_means(fit, by = "x"), "'g'")
  expect_error(adjusted_means(fit, weights = "cells"), "arg")
})
