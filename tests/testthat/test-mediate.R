# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

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
               var_indirect_effect(reg_coef_a = v(res, "a"),
                                   reg_coef_b = v(res, "b"),
                                   var_reg_coef_a = v(res, "se_a")^2,
                                   var_reg_coef_b = v(res, "se_b")^2)$value[1])
  expect_equal(v(res, "indirect_effect_ci_lower"),
               v(res, "indirect_effect") + qnorm(0.025) * se)
  expect_equal(v(res, "indirect_effect_ci_upper"),
               v(res, "indirect_effect") + qnorm(0.975) * se)
})

test_that("mediate() bootstrap and Monte Carlo intervals are seeded, sane, and skewed-capable", {
  skip_on_cran()  # two bootstraps, a Monte Carlo of 5000 draws, and a BCa run
  d <- .med_data()
  r1 <- mediate(d, "x", "m", "y", B = 300, seed = 7)
  r2 <- mediate(d, "x", "m", "y", B = 300, seed = 7)
  expect_identical(r1$value, r2$value)            # reproducible under seed
  expect_lt(v(r1, "indirect_effect_ci_lower"), v(r1, "indirect_effect"))
  expect_gt(v(r1, "indirect_effect_ci_upper"), v(r1, "indirect_effect"))
  mc <- mediate(d, "x", "m", "y", ci_method = "monte_carlo", B = 5000,
                seed = 7)
  so <- mediate(d, "x", "m", "y", ci_method = "sobel")
  # With a comfortably nonzero indirect effect the MC interval is close to
  # Sobel but need not coincide; both must cover the estimate.
  expect_lt(abs(v(mc, "indirect_effect_ci_lower") - v(so, "indirect_effect_ci_lower")), 0.05)
  bca <- mediate(d, "x", "m", "y", ci_method = "boot_bca", B = 300,
                 seed = 7)
  expect_lt(v(bca, "indirect_effect_ci_lower"), v(bca, "indirect_effect_ci_upper"))
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
  expect_error(mediate(d, "x", "m", "y", B = 1), "at least 2",
               class = "dmar_count_invalid")
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
    v(r, "indirect_effect_ci_lower") <= true_ab && true_ab <= v(r, "indirect_effect_ci_upper")
  })
  expect_gt(mean(cover), 0.88)   # 95 percent nominal, generous MC band
})

# The data of the rare predictor: one case carries all the variation in x,
# so a resample that omits it has a constant predictor and no indirect
# effect.
.med_rare_x <- function() {
  set.seed(113)
  n <- 20
  x <- c(3, rep(0, n - 1))
  m <- 0.5 * x + rnorm(n)
  y <- 0.6 * m + rnorm(n)
  data.frame(x = x, m = m, y = y)
}

# The data of the rare covariate category: 2 of 60 cases in category "b",
# so about one resample in eight holds no case of it.
.med_rare_category <- function() {
  set.seed(2026)
  n <- 60
  x <- rnorm(n)
  m <- 0.4 * x + rnorm(n)
  y <- 0.3 * m + 0.2 * x + rnorm(n)
  data.frame(x = x, m = m, y = y,
             g = factor(c("b", "b", rep("a", n - 2))))
}

# mediate() recomputed with lm() on formulas, for data with plain names:
# the 17 rows the table had before B_used was added.
.med_reference <- function(d, covariates = NULL, ci_method, B, seed,
                           conf_level = 0.95) {
  rhs_m <- paste(c("x", covariates), collapse = " + ")
  f_m <- as.formula(paste("m ~", rhs_m))
  f_y <- as.formula(paste("y ~ m +", rhs_m))
  f_t <- as.formula(paste("y ~", rhs_m))
  ab <- function(dd) {
    unname(coef(lm(f_m, dd))["x"] * coef(lm(f_y, dd))["m"])
  }
  fm <- lm(f_m, d); fy <- lm(f_y, d); ft <- lm(f_t, d)
  a <- unname(coef(fm)["x"]); b <- unname(coef(fy)["m"])
  cp <- unname(coef(fy)["x"])
  se_a <- sqrt(vcov(fm)["x", "x"]); se_b <- sqrt(vcov(fy)["m", "m"])
  se_ind <- sqrt(a^2 * se_b^2 + b^2 * se_a^2)
  probs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
  n <- nrow(d)
  lims <- withr::with_seed(seed, {
    if (ci_method %in% c("boot_percentile", "boot_bca")) {
      boots <- vapply(seq_len(B), function(i) {
        ab(d[sample.int(n, n, replace = TRUE), , drop = FALSE])
      }, numeric(1))
      if (ci_method == "boot_bca") {
        z0 <- qnorm(mean(boots < a * b))
        jack <- vapply(seq_len(n), function(i) ab(d[-i, , drop = FALSE]),
                       numeric(1))
        acc <- sum((mean(jack) - jack)^3) /
          (6 * sum((mean(jack) - jack)^2)^1.5)
        z <- z0 + qnorm(probs)
        probs <- pnorm(z0 + z / (1 - acc * z))
      }
      quantile(boots, probs, names = FALSE)
    } else if (ci_method == "monte_carlo") {
      a_sim <- rnorm(B, a, se_a)
      b_sim <- rnorm(B, b, se_b)
      quantile(a_sim * b_sim, probs, names = FALSE)
    } else {
      a * b + qnorm(probs) * se_ind
    }
  })
  c(indirect_effect = a * b, indirect_effect_ci_lower = lims[1],
    indirect_effect_ci_upper = lims[2], direct_effect = cp,
    direct_effect_ci_lower = confint(fy, "x", level = conf_level)[1],
    direct_effect_ci_upper = confint(fy, "x", level = conf_level)[2],
    total_effect = unname(coef(ft)["x"]),
    total_effect_ci_lower = confint(ft, "x", level = conf_level)[1],
    total_effect_ci_upper = confint(ft, "x", level = conf_level)[2],
    a = a, b = b, se_a = se_a, se_b = se_b, se_indirect = se_ind,
    proportion_mediated = a * b / unname(coef(ft)["x"]), N = n,
    B = if (ci_method == "sobel") NA_real_ else B)
}

test_that("a resample with a constant predictor is dropped, counted, and not fatal", {
  skip_on_cran()  # a bootstrap of 300 replications; the anchors run on CRAN
  d <- .med_rare_x()
  w <- NULL
  res <- withCallingHandlers(
    mediate(d, x = "x", m = "m", y = "y",
            ci_method = "boot_percentile", B = 300, seed = 113),
    warning = function(cnd) {
      w[[length(w) + 1L]] <<- cnd
      invokeRestart("muffleWarning")
    }
  )
  # One warning, classed, carrying the counts.
  expect_length(w, 1L)
  expect_s3_class(w[[1L]], "dmar_bootstrap_replications_dropped")
  expect_match(conditionMessage(w[[1L]]), "120 of 300")
  expect_identical(w[[1L]]$n_dropped, 120L)
  expect_identical(w[[1L]]$B_used, 180L)
  expect_identical(w[[1L]]$B, 300L)
  lims <- res$value[res$term %in% c("indirect_effect_ci_lower", "indirect_effect_ci_upper")]
  expect_length(lims, 2L)
  expect_true(all(is.finite(lims)))
  expect_equal(v(res, "B"), 300)
  expect_equal(v(res, "B_used"), 180)
  expect_identical(attr(res, "B"), 300L)
  expect_identical(attr(res, "B_used"), 180L)
  expect_identical(attr(res, "n_dropped"), 120L)
  expect_identical(attr(res, "n_covariate_redundant"), 0L)
  # The limits are the percentiles of the 180 resamples that hold case 1.
  ref <- withr::with_seed(113, vapply(seq_len(300), function(i) {
    idx <- sample.int(20, 20, replace = TRUE)
    if (!any(idx == 1L)) return(NA_real_)
    dd <- d[idx, ]
    unname(coef(lm(m ~ x, dd))["x"] * coef(lm(y ~ x + m, dd))["m"])
  }, numeric(1)))
  expect_equal(sum(is.na(ref)), 120L)
  expect_equal(lims, quantile(ref[!is.na(ref)], c(0.025, 0.975),
                              names = FALSE), tolerance = 1e-10)
})

test_that("replications dropped at the required B leave the interval reported, and a run settled below the count gives NA limits and a classed warning, not an error", {
  skip_on_cran()  # two bootstraps of 400 replications
  d <- .med_rare_x()
  collect <- function(expr) {
    w <- list()
    res <- withCallingHandlers(expr, warning = function(cnd) {
      w[[length(w) + 1L]] <<- cnd
      invokeRestart("muffleWarning")
    })
    list(res = res, w = w)
  }
  # The rule itself, so the requirement is not lowered here: of B = 400,
  # about four in ten resamples lose case 1, and fewer than the 400 a 95
  # percent interval must request remain. The count is held to the
  # request (ruled 2026-11-10), so the interval is reported from those
  # that remain, and the warning says that they are fewer than 400 and to
  # request more. (The first versions of 1.1.0 reported no interval here;
  # DMAR 1.0.x stopped when fewer than 100 remained.)
  testthat::local_mocked_bindings(
    .dmar_boot_required = function(kind = "quantile", tail = NULL,
                                   floor = 100L) 400L, .package = "DMAR")
  got <- collect(mediate(d, "x", "m", "y", B = 400, seed = 113))
  res <- got$res
  w <- got$w
  expect_length(w, 1L)
  expect_s3_class(w[[1L]], "dmar_bootstrap_replications_dropped")
  expect_identical(w[[1L]]$status, "dropped")
  expect_identical(w[[1L]]$B, 400L)
  expect_identical(w[[1L]]$B_required, 400L)
  expect_identical(w[[1L]]$B_used + w[[1L]]$n_dropped, 400L)
  expect_lt(w[[1L]]$B_used, 400L)
  expect_gte(w[[1L]]$B_used, 200L)
  expect_true(w[[1L]]$caution)
  expect_match(conditionMessage(w[[1L]]),
               "that remained are fewer, so request more bootstrap replications where some can fail, enough that at least 400 remain",
               fixed = TRUE)
  expect_true(all(is.finite(
    res$value[res$term %in% c("indirect_effect_ci_lower", "indirect_effect_ci_upper")])))
  expect_identical(attr(res, "bootstrap_status"), "dropped")
  expect_identical(attr(res, "B_used"), w[[1L]]$B_used)
  expect_identical(attr(res, "B_required"), 400L)
  # The same run settled as one whose request was below the count: NA
  # limits, every other row as the Sobel call gives it, a classed warning.
  local_count_unmet()
  got <- collect(mediate(d, "x", "m", "y", B = 400, seed = 113))
  res <- got$res
  w <- got$w
  expect_length(w, 1L)
  expect_s3_class(w[[1L]], "dmar_bootstrap_insufficient")
  expect_identical(w[[1L]]$status, "insufficient_count")
  expect_match(conditionMessage(w[[1L]]), "no interval is reported")
  expect_true(all(is.na(
    res$value[res$term %in% c("indirect_effect_ci_lower", "indirect_effect_ci_upper")])))
  sobel <- mediate(d, "x", "m", "y", ci_method = "sobel")
  same <- !(res$term %in% c("indirect_effect_ci_lower", "indirect_effect_ci_upper", "B",
                            "B_used"))
  expect_equal(res$value[same], sobel$value[same], tolerance = 1e-12)
  expect_identical(attr(res, "bootstrap_status"), "insufficient_count")
  expect_identical(attr(res, "B_used"), w[[1L]]$B_used)
  expect_identical(attr(res, "B_required"), 400L)
})

test_that("the BCa interval is refused when a leave-one-out effect is not estimable", {
  skip_on_cran()  # a bootstrap of 300 replications before the refusal
  d <- .med_rare_x()
  err <- tryCatch(
    suppressWarnings(mediate(d, "x", "m", "y", ci_method = "boot_bca",
                             B = 300, seed = 113)),
    error = function(e) e)
  expect_s3_class(err, "dmar_bca_undefined")
  expect_s3_class(err, "dmar_refusal")
  expect_match(conditionMessage(err), "with case 1 left out")
  expect_match(conditionMessage(err), "boot_percentile")
})

test_that("a rare covariate category neither aborts the call nor drops a resample", {
  skip_on_cran()  # five bootstraps of 400 replications
  d <- .med_rare_category()
  res <- expect_no_warning(
    mediate(d, "x", "m", "y", covariates = "g", B = 400, seed = 5))
  expect_equal(v(res, "B_used"), 400)
  expect_identical(attr(res, "n_dropped"), 0L)
  # The resamples with no case of the rare category are kept and counted.
  no_b <- withr::with_seed(5, sum(vapply(seq_len(400), function(i) {
    !any(sample.int(60, 60, replace = TRUE) <= 2L)
  }, logical(1))))
  expect_gt(no_b, 0L)
  expect_identical(attr(res, "n_covariate_redundant"), no_b)
  # The coding of the category does not matter: character, another
  # reference level, and a numeric 0/1 column give the same table.
  d_chr <- d; d_chr$g <- as.character(d$g)
  d_rel <- d; d_rel$g <- relevel(d$g, ref = "b")
  d_num <- d; d_num$g <- as.numeric(d$g == "b")
  for (alt in list(d_chr, d_rel, d_num)) {
    res_alt <- expect_no_warning(
      mediate(alt, "x", "m", "y", covariates = "g", B = 400, seed = 5))
    expect_equal(res_alt$value, res$value, tolerance = 1e-10)
  }
  fy <- lm(y ~ x + m + g, d)
  expect_equal(v(res, "b"), unname(coef(fy)["m"]), tolerance = 1e-10)
  expect_equal(v(res, "se_b"), sqrt(vcov(fy)["m", "m"]), tolerance = 1e-10)
})

test_that("a resample that aliases the predictor with a covariate is dropped", {
  skip_on_cran()  # a bootstrap of 2000 replications and its lm() recomputation
  set.seed(113)
  n <- 40
  x <- rep(c(1, 0), c(8, 32))
  g <- x
  g[c(1, 9)] <- 1 - g[c(1, 9)]
  m <- 0.5 * x + 0.8 * g + rnorm(n)
  y <- 0.4 * m + 0.3 * g + rnorm(n)
  d <- data.frame(x = x, m = m, y = y, g = g)
  w <- NULL
  res <- withCallingHandlers(
    mediate(d, "x", "m", "y", covariates = "g", B = 2000, seed = 113),
    warning = function(cnd) {
      w[[length(w) + 1L]] <<- cnd
      invokeRestart("muffleWarning")
    }
  )
  expect_length(w, 1L)
  expect_s3_class(w[[1L]], "dmar_bootstrap_replications_dropped")
  expect_identical(w[[1L]]$n_dropped, 248L)
  expect_equal(v(res, "B_used"), 1752)
  expect_identical(attr(res, "B_used"), 1752L)
  # a is estimable in a resample exactly when some level of g holds both
  # values of x (cell counts, no matrix decomposition). lm() returns a
  # finite product in the other resamples too, which is why they have to
  # be identified and not merely screened for NA.
  ref <- withr::with_seed(113, vapply(seq_len(2000), function(i) {
    dd <- d[sample.int(n, n, replace = TRUE), ]
    cells <- table(factor(dd$g, 0:1), factor(dd$x, 0:1)) > 0
    estimable <- any(cells[, 1L] & cells[, 2L])
    ab <- unname(coef(lm(m ~ x + g, dd))["x"] *
                   coef(lm(y ~ x + m + g, dd))["m"])
    c(estimable, ab)
  }, numeric(2)))
  expect_equal(sum(ref[1L, ] == 0), 248)
  expect_equal(
    c(v(res, "indirect_effect_ci_lower"), v(res, "indirect_effect_ci_upper")),
    quantile(ref[2L, ref[1L, ] == 1], c(0.025, 0.975), names = FALSE),
    tolerance = 1e-10)
})

test_that("the units and the origin of the predictor do not decide which resamples are kept", {
  skip_on_cran()  # four bootstraps of 400 replications
  d <- .med_rare_category()
  r0 <- mediate(d, "x", "m", "y", covariates = "g", B = 400, seed = 5)
  rows <- c("indirect_effect", "indirect_effect_ci_lower", "indirect_effect_ci_upper")
  for (k in c(1e-8, 1e8)) {
    d_k <- d; d_k$x <- d$x * k
    r_k <- expect_no_warning(
      mediate(d_k, "x", "m", "y", covariates = "g", B = 400, seed = 5))
    expect_equal(v(r_k, "B_used"), 400)
    # ab is a per unit of x times b, so it scales by the reciprocal.
    expect_equal(r_k$value[match(rows, r_k$term)] * k,
                 r0$value[match(rows, r0$term)], tolerance = 1e-10)
  }
  d_s <- d; d_s$x <- d$x + 1e9
  r_s <- expect_no_warning(
    mediate(d_s, "x", "m", "y", covariates = "g", B = 400, seed = 5))
  expect_equal(v(r_s, "B_used"), 400)
  expect_equal(r_s$value[match(rows, r_s$term)],
               r0$value[match(rows, r0$term)], tolerance = 1e-6)
})

test_that("mediate() takes column names as they are", {
  set.seed(9)
  n <- 120
  w1 <- rnorm(n); w2 <- rnorm(n); w3 <- rexp(n); w4 <- w1^2
  x <- rnorm(n) + 0.5 * w4
  m <- 0.5 * x + 0.6 * w4 + rnorm(n)
  y <- 0.3 * m + 0.5 * w4 + rnorm(n)
  plain <- data.frame(x, m, y, c1 = w1, c2 = w2, c3 = w3, c4 = w4)
  odd <- stats::setNames(plain, c("predictor score", "mediator`score",
                                  "outcome $ score", "w:x", ".", "log(w)",
                                  "w^2"))
  # An object named 'w' in the calling frame must not be reached by the
  # names "w:x" and "log(w)".
  w <- rep(c(-50, 50), length.out = n)
  r_plain <- mediate(plain, "x", "m", "y",
                     covariates = c("c1", "c2", "c3", "c4"), B = 150,
                     seed = 3)
  r_odd <- mediate(odd, "predictor score", "mediator`score",
                   "outcome $ score",
                   covariates = c("w:x", ".", "log(w)", "w^2"), B = 150,
                   seed = 3)
  expect_identical(r_odd$value, r_plain$value)
  # A column named w^2 beside w is its own covariate, not w a second time.
  sq <- data.frame(x, m, y, w = w1, check.names = FALSE)
  sq[["w^2"]] <- w4
  sq$w_sq <- w4
  r_sq <- mediate(sq, "x", "m", "y", covariates = c("w", "w^2"),
                  ci_method = "sobel")
  expect_equal(v(r_sq, "a"), unname(coef(lm(m ~ x + w + w_sq, sq))["x"]),
               tolerance = 1e-10)
  expect_equal(v(r_sq, "b"),
               unname(coef(lm(y ~ x + m + w + w_sq, sq))["m"]),
               tolerance = 1e-10)
  expect_gt(abs(v(r_sq, "a") - unname(coef(lm(m ~ x + w, sq))["x"])), 0.1)
})

test_that("mediate() refuses columns it cannot use, each with its reason", {
  d <- .med_data(60)
  d$z <- rnorm(60)
  d$x_factor <- factor(d$x > 0)
  d$x_logical <- d$x > 0
  d$x_character <- ifelse(d$x > 0, "high", "low")
  d$m_factor <- factor(d$m > 0)
  d$z_inf <- d$z; d$z_inf[4] <- Inf
  d$x_constant <- 1
  for (bad in c("x_factor", "x_logical", "x_character")) {
    expect_error(mediate(d, bad, "m", "y", ci_method = "sobel"),
                 "'x' must name a numeric column.*coded 0 and 1")
  }
  expect_error(mediate(d, "x", "m_factor", "y", ci_method = "sobel"),
               "'m' must name a numeric column")
  expect_error(mediate(d, "x", "x", "y", ci_method = "sobel"),
               "must name different columns; 'x' is named more than once")
  expect_error(mediate(d, "x", "m", "y", covariates = "y",
                       ci_method = "sobel"),
               "must name different columns; 'y' is named more than once")
  expect_error(mediate(d, "x", "m", "y", covariates = c("z", "z"),
                       ci_method = "sobel"),
               "must name different columns; 'z' is named more than once")
  expect_error(mediate(d, "x", "m", "y", covariates = c("z", NA),
                       ci_method = "sobel"),
               "'covariates' must name columns of 'data'")
  expect_error(mediate(d, NA_character_, "m", "y", ci_method = "sobel"),
               "'x' must name one column of 'data'")
  expect_error(mediate(d, "x", "m", "y", covariates = "z_inf",
                       ci_method = "sobel"),
               "must be finite")
  expect_error(mediate(d, "x_constant", "m", "y", ci_method = "sobel"),
               "not estimable in these data: 'x' is constant")
  d$m_alias <- 2 * d$x + 3 * d$z
  expect_error(mediate(d, "x", "m_alias", "y", covariates = "z",
                       ci_method = "sobel"),
               "'m' or 'x' is a linear combination")
})

test_that("mediate() keeps the inputs it accepted: time classes, integers, logical covariates", {
  d <- .med_data(80)
  d$x_date <- as.Date("2026-01-01") + round(10 * d$x)
  d$x_days <- as.numeric(d$x_date)
  d$x_int <- as.integer(round(3 * d$x))
  d$flag <- d$y + rnorm(80) > 0
  r_date <- mediate(d, "x_date", "m", "y", B = 150, seed = 2)
  r_days <- mediate(d, "x_days", "m", "y", B = 150, seed = 2)
  expect_equal(r_date$value, r_days$value, tolerance = 1e-8)
  expect_equal(v(r_date, "a"), unname(coef(lm(m ~ x_days, d))["x_days"]),
               tolerance = 1e-8)
  r_int <- mediate(d, "x_int", "m", "y", covariates = "flag",
                   ci_method = "sobel")
  fm <- lm(m ~ x_int + flag, d); fy <- lm(y ~ x_int + m + flag, d)
  expect_equal(v(r_int, "a"), unname(coef(fm)["x_int"]), tolerance = 1e-10)
  expect_equal(v(r_int, "b"), unname(coef(fy)["m"]), tolerance = 1e-10)
  expect_equal(v(r_int, "se_a"), sqrt(vcov(fm)["x_int", "x_int"]),
               tolerance = 1e-10)
  expect_equal(v(r_int, "direct_effect"), unname(coef(fy)["x_int"]),
               tolerance = 1e-10)
})

test_that("mediate() reads the residual degrees of freedom from the fit, not from a column count", {
  set.seed(2)
  n <- 16
  d <- data.frame(x = rnorm(n), m = rnorm(n), y = rnorm(n),
                  f2 = factor(rep(1:8, 2)),
                  f3 = factor(rep(1:8, each = 2)),
                  f14 = factor(1:16 %% 14))
  # The two factors give more coded columns than cases but a design of
  # lower rank, which lm() fits with residual degrees of freedom to spare.
  res <- mediate(d, "x", "m", "y", covariates = c("f2", "f3"),
                 ci_method = "sobel")
  fm <- lm(m ~ x + f2 + f3, d); fy <- lm(y ~ x + m + f2 + f3, d)
  expect_gt(fy$df.residual, 0)
  expect_gt(length(coef(fy)), n)
  expect_equal(v(res, "a"), unname(coef(fm)["x"]), tolerance = 1e-10)
  expect_equal(v(res, "b"), unname(coef(fy)["m"]), tolerance = 1e-10)
  expect_equal(v(res, "se_a"), sqrt(vcov(fm)["x", "x"]), tolerance = 1e-10)
  expect_equal(v(res, "se_b"), sqrt(vcov(fy)["m", "m"]), tolerance = 1e-10)
  expect_equal(c(v(res, "direct_effect_ci_lower"), v(res, "direct_effect_ci_upper")),
               unname(confint(fy, "x")[1, ]), tolerance = 1e-10)
  expect_error(mediate(d, "x", "m", "y", covariates = "f14",
                       ci_method = "sobel"),
               "no residual degrees of freedom")
})

test_that("a covariate with a single value is left out with a message", {
  d <- .med_data(60)
  d$site <- factor("a")
  d$wave <- "first"
  d$dose <- 2
  ref <- mediate(d, "x", "m", "y", B = 150, seed = 4)
  expect_message(
    res <- mediate(d, "x", "m", "y", covariates = "site", B = 150,
                   seed = 4),
    "Covariate 'site' has a single value")
  expect_identical(res$value, ref$value)
  expect_identical(attr(res, "n_covariate_redundant"), 0L)
  expect_message(
    res3 <- mediate(d, "x", "m", "y",
                    covariates = c("site", "wave", "dose"), B = 150,
                    seed = 4),
    "Covariates 'site', 'wave', 'dose' have a single value.*are left out")
  expect_identical(res3$value, ref$value)
  # A category emptied by listwise deletion is the same case.
  d$group <- factor(c("b", rep("a", 59)))
  d$y[1] <- NA
  ref_59 <- mediate(d, "x", "m", "y", B = 150, seed = 4)
  expect_message(
    res_59 <- mediate(d, "x", "m", "y", covariates = "group", B = 150,
                      seed = 4),
    "single value")
  expect_identical(res_59$value, ref_59$value)
})

test_that("mediate() agrees with lm() on every earlier row, four methods, with and without covariates", {
  skip_on_cran()  # eight calls, each recomputed with lm() in every resample
  set.seed(113)
  n <- 120
  x <- rnorm(n); z <- rnorm(n)
  g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  site <- sample(c("north", "south"), n, replace = TRUE)
  flag <- rnorm(n) > 0
  m <- 0.5 * x + 0.3 * z + 0.4 * (g == "b") + rnorm(n)
  y <- 0.2 * x + 0.4 * m + 0.2 * z + 0.3 * flag + rnorm(n)
  d <- data.frame(x, m, y, z, g, site, flag)
  for (covs in list(NULL, c("z", "g", "site", "flag"))) {
    for (method in c("boot_percentile", "boot_bca", "monte_carlo",
                     "sobel")) {
      res <- mediate(d, "x", "m", "y", covariates = covs,
                     ci_method = method, B = 150, seed = 21,
                     conf_level = 0.90)
      ref <- .med_reference(d, covs, method, B = 150, seed = 21,
                            conf_level = 0.90)
      expect_identical(nrow(res), 18L)
      expect_identical(res$term[1:17], names(ref))
      expect_equal(res$value[1:17], unname(ref), tolerance = 1e-10,
                   info = paste(method, length(covs)))
    }
  }
})

test_that("mediate() records the replications where dmar_metadata() reads them", {
  d <- .med_data()
  so <- mediate(d, "x", "m", "y", ci_method = "sobel")
  expect_true(is.na(v(so, "B")))
  expect_true(is.na(v(so, "B_used")))
  expect_null(attr(so, "B"))
  expect_null(attr(so, "B_used"))
  expect_null(attr(so, "n_dropped"))
  expect_null(dmar_metadata(so)$bootstrap)
  bo <- mediate(d, "x", "m", "y", B = 150, seed = 8)
  expect_identical(bo$term[17:18], c("B", "B_used"))
  expect_identical(dmar_metadata(bo)$bootstrap$B, 150L)
  expect_identical(dmar_metadata(bo)$bootstrap$used, 150L)
  expect_identical(attr(bo, "n_dropped"), 0L)
  mc <- mediate(d, "x", "m", "y", ci_method = "monte_carlo", B = 500,
                seed = 8)
  expect_equal(v(mc, "B_used"), 500)
  expect_identical(attr(mc, "B_used"), 500L)
  expect_identical(attr(mc, "n_dropped"), 0L)
  expect_identical(attr(mc, "bootstrap_status"), "complete")
})

test_that("mediate() reads a one-column matrix as its column", {
  # scale() leaves a one-column matrix in a data frame; the earlier build
  # read it, and a check on dim() had refused it.
  set.seed(31337)
  N <- 90
  w <- rnorm(N)
  x <- rnorm(N) + 0.4 * w
  m <- 0.5 * x + 0.3 * w + rnorm(N)
  y <- 0.4 * m + 0.2 * x - 0.3 * w + rnorm(N)
  z <- data.frame(x, m, y, w)
  z$xs <- scale(z$x)
  z$xn <- as.numeric(z$xs)
  z$ms <- scale(z$m)
  z$mn <- as.numeric(z$ms)
  z$ya <- array(z$y, dim = N)
  ref <- mediate(z, "xn", "mn", "y", "w", B = 200, seed = 6)
  res <- mediate(z, "xs", "ms", "ya", "w", B = 200, seed = 6)
  expect_equal(res$value, ref$value)
  # A matrix of two columns is not one variable.
  z$two <- cbind(z$x, z$w)
  expect_error(mediate(z, "two", "m", "y", B = 200, seed = 6),
               "'x' must name a numeric column")
  # A logical outcome is refused, as a logical predictor or mediator is.
  z$yl <- z$y > 0
  expect_error(mediate(z, "x", "m", "yl", B = 200, seed = 6),
               "'y' must name a numeric column")
})

test_that("mediate() counts covariate redundancy against the rank in the data", {
  # A covariate entered beside twice itself adds no column to the model in
  # the data, so no resample loses rank among the covariates relative to
  # the data, and the table is that of the covariate entered once.
  set.seed(31337)
  N <- 90
  w <- rnorm(N)
  x <- rnorm(N) + 0.4 * w
  m <- 0.5 * x + 0.3 * w + rnorm(N)
  y <- 0.4 * m + 0.2 * x - 0.3 * w + rnorm(N)
  z <- data.frame(x, m, y, w, w2 = 2 * w)
  both <- mediate(z, "x", "m", "y", c("w", "w2"), B = 200, seed = 6)
  one  <- mediate(z, "x", "m", "y", "w", B = 200, seed = 6)
  expect_identical(attr(both, "n_covariate_redundant"), 0L)
  expect_equal(both$value, one$value, tolerance = 1e-8)
})

test_that("mediate() requires a finite B and refuses a BCa interval without the dropped-replications warning", {
  set.seed(1)
  d <- data.frame(x = c(3, rep(0, 19)))
  d$m <- d$x + rnorm(20)
  d$y <- d$m + rnorm(20)
  expect_error(mediate(d, "x", "m", "y", B = Inf), "single whole number",
               class = "dmar_count_invalid")
  # The Sobel interval draws nothing, so B is not read for it.
  expect_true(is.na(v(mediate(d, "x", "m", "y", B = Inf,
                              ci_method = "sobel"), "B")))
  # The interval is refused, so no warning describes an interval computed
  # from the replications that remained.
  expect_no_warning(
    err <- tryCatch(mediate(d, "x", "m", "y", ci_method = "boot_bca",
                            B = 300, seed = 113),
                    dmar_bca_undefined = function(e) e))
  expect_s3_class(err, "dmar_bca_undefined")
})
