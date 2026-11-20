## from tests/testthat/test-adjusted_means.R
local({
  # adjusted_means() one-way cells and marginals vs emmeans (emmeans 2.0.3,
  # 2026-08-09): estimate, SE, df, and both CI limits at 1e-10.
  agree <- function(am, emm, facs, tol = 1e-10) {
    s <- summary(emm)
    key <- function(d) do.call(paste, c(unname(lapply(d[facs], as.character)),
                                        list(sep = "|")))
    i <- match(key(am), key(s))
    max(abs(am$estimate - s$emmean[i]), abs(am$se - s$SE[i]),
        abs(attr(am, "df_residual") - s$df[i]),
        abs(am$ci_lower - s$lower.CL[i]),
        abs(am$ci_upper - s$upper.CL[i])) < tol
  }
  fit <- lm(weight ~ group, data = PlantGrowth)
  stopifnot(
    agree(DMAR::adjusted_means(fit), emmeans::emmeans(fit, ~ group), "group"),
    agree(DMAR::adjusted_means(fit, by = "group", weights = "proportional"),
          emmeans::emmeans(fit, ~ group, weights = "proportional"), "group"))
})

## from tests/testthat/test-adjusted_means.R
local({
  # Unbalanced 2 x 2 with interaction: cell means and both marginal-mean
  # weightings vs emmeans (emmeans 2.0.3, 2026-08-09). "equal" is emmeans's
  # default; "proportional" is emmeans's weights = "proportional".
  agree <- function(am, emm, facs, tol = 1e-10) {
    s <- summary(emm)
    key <- function(d) do.call(paste, c(unname(lapply(d[facs], as.character)),
                                        list(sep = "|")))
    i <- match(key(am), key(s))
    max(abs(am$estimate - s$emmean[i]), abs(am$se - s$SE[i]),
        abs(attr(am, "df_residual") - s$df[i]),
        abs(am$ci_lower - s$lower.CL[i]),
        abs(am$ci_upper - s$upper.CL[i])) < tol
  }
  set.seed(407)
  n_cell <- c(7, 12, 9, 15)
  d <- do.call(rbind, lapply(1:4, function(i) {
    g <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2"))
    g[rep(i, n_cell[i]), ]
  }))
  d$A <- factor(d$A); d$B <- factor(d$B)
  d$y <- 2 + (d$A == "a2") + 2 * (d$B == "b2") +
    1.5 * (d$A == "a2") * (d$B == "b2") + rnorm(nrow(d))
  fit <- lm(y ~ A * B, data = d)
  stopifnot(
    agree(DMAR::adjusted_means(fit), emmeans::emmeans(fit, ~ A * B),
          c("A", "B")),
    agree(DMAR::adjusted_means(fit, by = "A"), emmeans::emmeans(fit, ~ A),
          "A"),
    agree(DMAR::adjusted_means(fit, by = "A", weights = "proportional"),
          emmeans::emmeans(fit, ~ A, weights = "proportional"), "A"),
    agree(DMAR::adjusted_means(fit, by = "B"), emmeans::emmeans(fit, ~ B),
          "B"),
    agree(DMAR::adjusted_means(fit, by = "B", weights = "proportional"),
          emmeans::emmeans(fit, ~ B, weights = "proportional"), "B"))
})

## from tests/testthat/test-adjusted_means.R
local({
  # Unbalanced three-way: proportional averaging weights each averaged-over
  # factor combination by its joint frequency, matching emmeans
  # (emmeans 2.0.3, 2026-08-09).
  agree <- function(am, emm, facs, tol = 1e-10) {
    s <- summary(emm)
    key <- function(d) do.call(paste, c(unname(lapply(d[facs], as.character)),
                                        list(sep = "|")))
    i <- match(key(am), key(s))
    max(abs(am$estimate - s$emmean[i]), abs(am$se - s$SE[i]),
        abs(attr(am, "df_residual") - s$df[i]),
        abs(am$ci_lower - s$lower.CL[i]),
        abs(am$ci_upper - s$upper.CL[i])) < tol
  }
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
  stopifnot(
    agree(DMAR::adjusted_means(fit), emmeans::emmeans(fit, ~ A * B * C),
          c("A", "B", "C")),
    agree(DMAR::adjusted_means(fit, by = c("A", "B")),
          emmeans::emmeans(fit, ~ A * B), c("A", "B")),
    agree(DMAR::adjusted_means(fit, by = c("A", "B"),
                               weights = "proportional"),
          emmeans::emmeans(fit, ~ A * B, weights = "proportional"),
          c("A", "B")),
    agree(DMAR::adjusted_means(fit, by = "C", weights = "proportional"),
          emmeans::emmeans(fit, ~ C, weights = "proportional"), "C"))
})

## from tests/testthat/test-adjusted_means.R
local({
  # One-way ANCOVA: adjusted means at the covariate mean vs emmeans
  # (emmeans 2.0.3, 2026-08-09).
  agree <- function(am, emm, facs, tol = 1e-10) {
    s <- summary(emm)
    key <- function(d) do.call(paste, c(unname(lapply(d[facs], as.character)),
                                        list(sep = "|")))
    i <- match(key(am), key(s))
    max(abs(am$estimate - s$emmean[i]), abs(am$se - s$SE[i]),
        abs(attr(am, "df_residual") - s$df[i]),
        abs(am$ci_lower - s$lower.CL[i]),
        abs(am$ci_upper - s$upper.CL[i])) < tol
  }
  set.seed(113)
  d <- data.frame(g = factor(rep(c("g1", "g2", "g3"),
                                 times = c(12, 18, 15))))
  d$x <- rnorm(45, 10, 2)
  d$y <- 3 + 2 * (d$g == "g2") + 4 * (d$g == "g3") + 0.6 * d$x + rnorm(45)
  fit <- lm(y ~ g + x, data = d)
  stopifnot(agree(DMAR::adjusted_means(fit), emmeans::emmeans(fit, ~ g), "g"))
})

## from tests/testthat/test-adjusted_means.R
local({
  # Two-factor ANCOVA with two covariates: cells and proportional B margins
  # vs emmeans (emmeans 2.0.3, 2026-08-09).
  agree <- function(am, emm, facs, tol = 1e-10) {
    s <- summary(emm)
    key <- function(d) do.call(paste, c(unname(lapply(d[facs], as.character)),
                                        list(sep = "|")))
    i <- match(key(am), key(s))
    max(abs(am$estimate - s$emmean[i]), abs(am$se - s$SE[i]),
        abs(attr(am, "df_residual") - s$df[i]),
        abs(am$ci_lower - s$lower.CL[i]),
        abs(am$ci_upper - s$upper.CL[i])) < tol
  }
  set.seed(509)
  d <- data.frame(
    A = factor(rep(c("a1", "a2"), times = c(28, 36))),
    B = factor(rep(rep(c("b1", "b2"), times = c(9, 7)), 4)),
    x1 = rnorm(64), x2 = runif(64, -2, 2)
  )
  d$y <- 1 + (d$A == "a2") + 0.5 * (d$B == "b2") + 0.7 * d$x1 -
    0.4 * d$x2 + rnorm(64)
  fit <- lm(y ~ A * B + x1 + x2, data = d)
  stopifnot(
    agree(DMAR::adjusted_means(fit), emmeans::emmeans(fit, ~ A * B),
          c("A", "B")),
    agree(DMAR::adjusted_means(fit, by = "B", weights = "proportional"),
          emmeans::emmeans(fit, ~ B, weights = "proportional"), "B"))
  # The same fit as an aov, at conf_level 0.90.
  fit_aov <- aov(y ~ A * B + x1 + x2, data = d)
  s90 <- summary(emmeans::emmeans(fit_aov, ~ A,
                                  options = list(level = 0.90)))
  am90 <- DMAR::adjusted_means(fit_aov, by = "A", conf_level = 0.90)
  i <- match(am90$A, as.character(s90$A))
  stopifnot(max(abs(am90$estimate - s90$emmean[i]),
                abs(am90$se - s90$SE[i]),
                abs(am90$ci_lower - s90$lower.CL[i]),
                abs(am90$ci_upper - s90$upper.CL[i])) < 1e-10)
})

## from tests/testthat/test-adjusted_means.R
local({
  # Covariate-by-factor interaction: each cell mean is the prediction at the
  # covariate mean under that cell's own slope, matching emmeans
  # (emmeans 2.0.3, 2026-08-09).
  set.seed(617)
  d <- data.frame(A = factor(rep(c("lo", "mid", "hi"),
                                 times = c(14, 11, 17)),
                             levels = c("lo", "mid", "hi")))
  d$x <- rnorm(42, 5, 1.5)
  d$y <- 2 + (d$A == "mid") + 3 * (d$A == "hi") +
    (0.5 + 0.4 * (d$A == "hi")) * d$x + rnorm(42)
  fit <- lm(y ~ A * x, data = d)
  am <- DMAR::adjusted_means(fit)
  s <- summary(emmeans::emmeans(fit, ~ A))
  i <- match(am$A, as.character(s$A))
  stopifnot(max(abs(am$estimate - s$emmean[i]), abs(am$se - s$SE[i]),
                abs(am$ci_lower - s$lower.CL[i]),
                abs(am$ci_upper - s$upper.CL[i])) < 1e-10)
})

## from tests/testthat/test-adjusted_means.R
local({
  # Transformed covariates: the poly(x, 2) basis and log(x) are evaluated at
  # the mean of the raw covariate, as emmeans does (emmeans 2.0.3,
  # 2026-08-09).
  set.seed(701)
  d1 <- data.frame(A = factor(rep(c("p", "q"), each = 24)))
  d1$x <- runif(48, 0, 4)
  d1$y <- 1 + 2 * (d1$A == "q") + 0.8 * d1$x - 0.15 * d1$x^2 +
    rnorm(48, 0, 0.7)
  fit1 <- lm(y ~ A + poly(x, 2), data = d1)
  am1 <- DMAR::adjusted_means(fit1)
  s1 <- summary(emmeans::emmeans(fit1, ~ A))
  i1 <- match(am1$A, as.character(s1$A))
  stopifnot(max(abs(am1$estimate - s1$emmean[i1]),
                abs(am1$se - s1$SE[i1]),
                abs(am1$ci_lower - s1$lower.CL[i1]),
                abs(am1$ci_upper - s1$upper.CL[i1])) < 1e-10)

  set.seed(809)
  d2 <- data.frame(A = factor(rep(c("u", "v"), times = c(21, 17))))
  d2$x <- rexp(38) + 0.5
  d2$y <- 2 + (d2$A == "v") + 1.4 * log(d2$x) + rnorm(38, 0, 0.6)
  fit2 <- lm(y ~ A + log(x), data = d2)
  am2 <- DMAR::adjusted_means(fit2)
  s2 <- summary(emmeans::emmeans(fit2, ~ A))
  i2 <- match(am2$A, as.character(s2$A))
  stopifnot(max(abs(am2$estimate - s2$emmean[i2]),
                abs(am2$se - s2$SE[i2]),
                abs(am2$ci_lower - s2$lower.CL[i2]),
                abs(am2$ci_upper - s2$upper.CL[i2])) < 1e-10)
})

## from tests/testthat/test-adjusted_means.R
local({
  # A character predictor uses the fit's sorted level set, matching emmeans
  # (emmeans 2.0.3, 2026-08-09).
  set.seed(907)
  d <- data.frame(g = rep(c("zeta", "alpha", "mid"), times = c(13, 16, 11)),
                  stringsAsFactors = FALSE)
  d$x <- rnorm(40)
  d$y <- 1 + (d$g == "mid") + 2 * (d$g == "zeta") + 0.5 * d$x + rnorm(40)
  fit <- lm(y ~ g + x, data = d)
  am <- DMAR::adjusted_means(fit)
  s <- summary(emmeans::emmeans(fit, ~ g))
  i <- match(am$g, as.character(s$g))
  stopifnot(identical(am$g, c("alpha", "mid", "zeta")),
            max(abs(am$estimate - s$emmean[i]), abs(am$se - s$SE[i]),
                abs(am$ci_lower - s$lower.CL[i]),
                abs(am$ci_upper - s$upper.CL[i])) < 1e-10)
})

## from tests/testthat/test-adjusted_means.R
local({
  # An ordered factor keeps its polynomial contrasts; cells and proportional
  # margins match emmeans (emmeans 2.0.3, 2026-08-09).
  agree <- function(am, emm, facs, tol = 1e-10) {
    s <- summary(emm)
    key <- function(d) do.call(paste, c(unname(lapply(d[facs], as.character)),
                                        list(sep = "|")))
    i <- match(key(am), key(s))
    max(abs(am$estimate - s$emmean[i]), abs(am$se - s$SE[i]),
        abs(attr(am, "df_residual") - s$df[i]),
        abs(am$ci_lower - s$lower.CL[i]),
        abs(am$ci_upper - s$upper.CL[i])) < tol
  }
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
  stopifnot(
    agree(DMAR::adjusted_means(fit), emmeans::emmeans(fit, ~ dose * B),
          c("dose", "B")),
    agree(DMAR::adjusted_means(fit, by = "dose", weights = "proportional"),
          emmeans::emmeans(fit, ~ dose, weights = "proportional"), "dose"))
})

## from tests/testthat/test-adjusted_means.R
local({
  # Rank-deficient design (empty (a2, b2) cell): emmeans reports the cell
  # and the a2 margin as nonEst (NA); adjusted_means() refuses, naming
  # exactly those rows (emmeans 2.0.3, 2026-08-09).
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
  s_cells <- summary(emmeans::emmeans(fit, ~ A * B))
  s_margA <- summary(emmeans::emmeans(fit, ~ A))
  msg <- tryCatch(DMAR::adjusted_means(fit), error = conditionMessage)
  msg_by <- tryCatch(DMAR::adjusted_means(fit, by = "A"),
                     error = conditionMessage)
  stopifnot(
    sum(is.na(s_cells$emmean)) == 1L,
    is.na(s_cells$emmean[s_cells$A == "a2" & s_cells$B == "b2"]),
    is.na(s_margA$emmean[s_margA$A == "a2"]),
    !is.na(s_margA$emmean[s_margA$A == "a1"]),
    is.character(msg), grepl("A = a2, B = b2", msg, fixed = TRUE),
    !grepl("a1", msg),
    is.character(msg_by), grepl("A = a2", msg_by, fixed = TRUE),
    !grepl("a1", msg_by))
})

## from tests/testthat/test-adjusted_means.R
local({
  # A weighted lm: proportional averaging uses the total prior weight of
  # each averaged-over combination, matching emmeans's .wgt.
  # (emmeans 2.0.3, 2026-08-09).
  agree <- function(am, emm, facs, tol = 1e-10) {
    s <- summary(emm)
    key <- function(d) do.call(paste, c(unname(lapply(d[facs], as.character)),
                                        list(sep = "|")))
    i <- match(key(am), key(s))
    max(abs(am$estimate - s$emmean[i]), abs(am$se - s$SE[i]),
        abs(attr(am, "df_residual") - s$df[i]),
        abs(am$ci_lower - s$lower.CL[i]),
        abs(am$ci_upper - s$upper.CL[i])) < tol
  }
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
  stopifnot(
    agree(DMAR::adjusted_means(fit), emmeans::emmeans(fit, ~ A * B),
          c("A", "B")),
    agree(DMAR::adjusted_means(fit, by = "A", weights = "proportional"),
          emmeans::emmeans(fit, ~ A, weights = "proportional"), "A"))
})

## from tests/testthat/test-adjusted_means.R
local({
  # Rows dropped by na.action or subset= are dropped from the covariate mean
  # too, matching emmeans (emmeans 2.0.3, 2026-08-09).
  set.seed(29)
  d <- data.frame(g = factor(rep(c("g1", "g2", "g3"), 14)), x = rnorm(42))
  d$y <- as.numeric(d$g) + 0.7 * d$x + rnorm(42)
  d$x[c(4, 9)] <- NA
  fit <- lm(y ~ g + log(x + 10), data = d, subset = -3)
  am <- DMAR::adjusted_means(fit)
  s <- summary(emmeans::emmeans(fit, ~ g))
  i <- match(am$g, as.character(s$g))
  stopifnot(max(abs(am$estimate - s$emmean[i]), abs(am$se - s$SE[i]),
                abs(am$ci_lower - s$lower.CL[i]),
                abs(am$ci_upper - s$upper.CL[i])) < 1e-10)
})

## from tests/testthat/test-adjusted_means.R
local({
  # A logical predictor defines cells like a two-level factor, matching
  # emmeans (emmeans 2.0.3, 2026-08-09).
  set.seed(23)
  d <- data.frame(trt = rep(c(TRUE, FALSE), times = c(18, 22)),
                  x = rnorm(40))
  d$y <- 1 + 2 * d$trt + 0.5 * d$x + rnorm(40)
  fit <- lm(y ~ trt + x, data = d)
  am <- DMAR::adjusted_means(fit)
  s <- summary(emmeans::emmeans(fit, ~ trt))
  i <- match(am$trt, as.character(s$trt))
  stopifnot(max(abs(am$estimate - s$emmean[i]), abs(am$se - s$SE[i]),
                abs(am$ci_lower - s$lower.CL[i]),
                abs(am$ci_upper - s$upper.CL[i])) < 1e-10)
})
