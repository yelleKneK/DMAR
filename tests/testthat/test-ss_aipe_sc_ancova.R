## assurance_method = "analytic", the default: the published step of Lai
## and Kelley (2012), from the plan without assurance, which does not
## compute the probability it attains. At these values the plans are those
## of DMAR 1.0.x, whose plan without assurance, 127, was already the
## smallest, and those of the published step of the first 1.1.0 builds
## (79a03795).
test_that("ss_aipe_sc_ancova() returns the published plan by default", {
  plan <- function(div, a, ...)
    ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5, divisor = div,
                      c_weights = c(1, -1), width = .5, assurance = a, ...)
  expect_identical(eval(formals(ss_aipe_sc_ancova)$assurance_method),
                   c("analytic", "simulation"))
  res <- plan("s_ancova", .8)
  expect_identical(res$term, c("necessary_n_per_group", "assurance",
                               "assurance_attained"))
  expect_identical(res$value, c(129, .8, NA))
  expect_identical(attr(res, "assurance_method"), "analytic")
  expect_identical(attr(res, "assurance_probability"), "not computed")
  expect_null(attr(res, "assurance_mc_se"))
  expect_identical(plan("s_ancova", .8, assurance_method = "analytic"), res)
  footer <- paste(utils::capture.output(print(res)), collapse = " ")
  expect_match(footer, paste0(
    "Assurance method: analytic; the plan is the published approximation, ",
    "which does not compute the probability"), fixed = TRUE)
  expect_match(footer, paste0(
    "assurance_method = \"simulation\" computes it and returns the smallest ",
    "plan that attains it"), fixed = TRUE)
  expect_identical(plan("s_ancova", .9)$value[1], 130)
  expect_identical(plan("s_anova", .8)$value[1], 98)
  expect_identical(plan("s_anova", .9)$value[1], 100)
  expect_identical(
    ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5, c_weights = c(1, -1),
                      width = .5)$value, 127)
  # The step starts from the smallest plan without assurance: at a contrast
  # of 0 and width .5 that is 123, and the published plan is 124 (DMAR 1.0.0
  # gave 125 for both, from a start above the smallest).
  expect_identical(
    ss_aipe_sc_ancova(true_std_contrast = 0, c_weights = c(1, -1), width = .5,
                      assurance = .9)$value[1], 124)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_aipe_sc_ancova() returns a tidy single-row data frame", {
  res <- ss_aipe_sc_ancova(true_std_contrast = .8, width = .5, c_weights = c(.5, .5, 0, -1))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("term", "value"))
  expect_equal(res$term, "necessary_n_per_group")
  expect_true(is.numeric(res$value) && res$value > 0)
})

test_that("ss_aipe_sc_ancova() ratio path returns a smaller N than the unconditioned plan", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # ratio < 1 represents an effective reduction in residual variance from a
  # covariate, so fewer subjects are needed.
  base   <- ss_aipe_sc_ancova(true_std_contrast = .8, width = .5, c_weights = c(.5, .5, 0, -1))
  ratio  <- ss_aipe_sc_ancova(true_std_contrast = .8, true_sd_ratio = .6, width = .5,
                              c_weights = c(.5, .5, 0, -1), divisor = "s_anova")
  expect_lt(ratio$value, base$value)
})

test_that("ss_aipe_sc_ancova() accepts explicit alpha_lower/alpha_upper (s_ancova)", {
  ref <- ss_aipe_sc_ancova(true_std_contrast = .8, width = .5, c_weights = c(.5, .5, 0, -1))
  alt <- ss_aipe_sc_ancova(true_std_contrast = .8, width = .5, c_weights = c(.5, .5, 0, -1),
                           conf_level = NULL,
                           alpha_lower = .025, alpha_upper = .025)
  expect_equal(alt$value, ref$value)
})

test_that("ss_aipe_sc_ancova() refuses unequal tails, which ci_sc_ancova() cannot report", {
  # ci_sc_ancova() takes conf_level only (register finding 9).
  for (div in c("s_ancova", "s_anova")) {
    err <- tryCatch(
      ss_aipe_sc_ancova(true_std_contrast = .8, true_sd_ratio = if (div == "s_anova") .6,
                        width = .5, c_weights = c(.5, .5, 0, -1), divisor = div,
                        conf_level = NULL, alpha_lower = .01, alpha_upper = .04),
      error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, c("alpha_lower", "alpha_upper"))
    expect_match(conditionMessage(err), "must be equal")
  }
  # Equal tails give the plan of the conf_level they define.
  expect_identical(
    ss_aipe_sc_ancova(true_std_contrast = .8, true_sd_ratio = .6, width = .5,
                      c_weights = c(.5, .5, 0, -1), divisor = "s_anova",
                      conf_level = NULL, alpha_lower = .05, alpha_upper = .05)$value,
    ss_aipe_sc_ancova(true_std_contrast = .8, true_sd_ratio = .6, width = .5,
                      c_weights = c(.5, .5, 0, -1), divisor = "s_anova",
                      conf_level = .90)$value)
})

test_that("ss_aipe_sc_ancova() rejects mixing conf_level and alphas", {
  expect_error(
    ss_aipe_sc_ancova(true_std_contrast = .8, width = .5, c_weights = c(.5, .5, 0, -1),
                      conf_level = .95, alpha_lower = .025, alpha_upper = .025),
    "cannot mix them"
  )
})

test_that("ss_aipe_sc_ancova() s_anova assurance path runs end-to-end", {
  # The s_anova + assurance branch previously also forwarded conf_level alongside
  # alphas to ci_nc_t, which is the same bug pattern as s_ancova.
  res <- ss_aipe_sc_ancova(true_std_contrast = .8, true_sd_ratio = .6, width = .5,
                           c_weights = c(.5, .5, 0, -1), divisor = "s_anova",
                           assurance = .90)
  expect_s3_class(res, "data.frame")
  expect_equal(res$term, c("necessary_n_per_group", "assurance",
                           "assurance_attained"))
})

test_that("ss_aipe_sc_ancova() assurance path no longer errors on inner ci_nc_t calls", {
  # Regression test: same bug as ss_aipe_sc() -- the inner ci_nc_t
  # calls forwarded the default conf_level = .95 alongside alpha_lower /
  # alpha_upper, which ci_nc_t refuses. Those calls are in the published
  # step, the default (the rows of the simulated step are checked below).
  res <- ss_aipe_sc_ancova(true_std_contrast = .8, width = .5,
                           c_weights = c(.5, .5, 0, -1), assurance = .90)
  expect_s3_class(res, "data.frame")
  expect_equal(res$term, c("necessary_n_per_group", "assurance",
                           "assurance_attained"))
})

test_that("ss_aipe_sc_ancova() checks assurance_method, G, and the ratio of the standard deviations", {
  call <- function(...) tryCatch(
    ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5, c_weights = c(1, -1),
                      width = .5, assurance = .8, ...), error = function(e) e)
  for (bad in list(NA, "yes", "sim", TRUE, c("simulation", "analytic"))) {
    err <- call(assurance_method = bad)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "assurance_method")
  }
  expect_match(conditionMessage(call(assurance_method = "sim")),
               "'assurance_method' must be \"analytic\" (the default) or \"simulation\"",
               fixed = TRUE)
  expect_s3_class(call(G = 50, assurance_method = "simulation"),
                  "dmar_count_invalid")
  # G is not used, and not checked, by the published plan, nor without
  # assurance.
  expect_identical(call(G = 50)$value[1], 129)
  expect_identical(ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5,
                                     c_weights = c(1, -1), width = .5, G = 50,
                                     assurance_method = "simulation")$value, 127)
  # With divisor = "s_anova" the ratio sqrt(1 - rho^2) is in (0, 1].
  for (args in list(list(true_std_contrast = .5, true_cor = 1),
                    list(true_std_contrast = .5, true_sd_ratio = 0),
                    list(true_contrast = 1, true_sd_anova = 1, true_sd_ancova = 2))) {
    err <- tryCatch(do.call(ss_aipe_sc_ancova,
                            c(args, list(divisor = "s_anova", c_weights = c(1, -1),
                                         width = .5))),
                    error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
  }
})

# One raw-data study: the ci_sc_ancova() width from summaries computed here,
# and the five statistics of the simulated step (R/ss_aipe_sc_internals.R).
.sc_ancova_raw_study <- function(n, cw, delta, rho, divisor) {
  J <- length(cw)
  sigma_e <- 2
  sigma_x <- 3
  sigma_anova <- sigma_e / sqrt(1 - rho^2)
  beta <- rho * sigma_anova / sigma_x
  psi <- delta * if (divisor == "s_ancova") sigma_e else sigma_anova
  mu <- 5 + cw * psi / sum(cw^2)
  g <- rep(seq_len(J), each = n)
  x <- stats::rnorm(n * J, 10, sigma_x)
  e <- stats::rnorm(n * J, 0, sigma_e)
  y <- mu[g] + beta * x + e
  xm <- as.vector(tapply(x, g, mean))
  ym <- as.vector(tapply(y, g, mean))
  xc <- x - xm[g]
  yc <- y - ym[g]
  ssw_x <- sum(xc^2)
  b_w <- sum(xc * yc) / ssw_x
  sse <- sum((yc - b_w * xc)^2)
  ci <- ci_sc_ancova(adj_means = ym - b_w * (xm - mean(x)),
                     s_anova = sqrt(sum(yc^2) / (n * J - J)),
                     s_ancova = sqrt(sse / (n * J - J - 1)),
                     standardizer = divisor, c_weights = cw, n_per_group = n,
                     cov_means = xm, SSwithin_x = ssw_x)
  C <- sum(cw^2)
  list(width = ci$value[3] - ci$value[1],
       z_x = sum(cw * xm) / (sigma_x * sqrt(C / n)),
       u_w = stats::pchisq(ssw_x / sigma_x^2, J * (n - 1)),
       z_e = sum(cw * as.vector(tapply(e, g, mean))) / (sigma_e * sqrt(C / n)),
       z_b = (b_w - beta) * sqrt(ssw_x) / sigma_e,
       u_q = stats::pchisq(sse / sigma_e^2, J * (n - 1) - 1))
}

test_that("the simulated assurance step rebuilds the ci_sc_ancova() width from its five statistics", {
  skip_on_cran()
  set.seed(1112)
  for (d in list(list(n = 12, cw = c(1, -1), delta = .5, rho = .5),
                 list(n = 9, cw = c(.5, .5, 0, -1), delta = -.8, rho = .7))) {
    for (div in c("s_ancova", "s_anova")) {
      studies <- lapply(1:25, function(r)
        .sc_ancova_raw_study(d$n, d$cw, d$delta, d$rho, div))
      base <- lapply(c(z_x = "z_x", u_w = "u_w", z_e = "z_e", z_b = "z_b",
                       u_q = "u_q"), function(k) vapply(studies, `[[`, 0, k))
      ratio <- sqrt(1 - d$rho^2)
      rebuilt <- DMAR:::.sc_ancova_widths(
        base, d$n, d$cw, if (div == "s_ancova") d$delta else d$delta / ratio,
        div, if (div == "s_ancova") 1 else ratio, .025, .025)
      expect_equal(rebuilt, vapply(studies, `[[`, 0, "width"), tolerance = 1e-8)
    }
  }
})

test_that("the noncentrality distance the simulated step interpolates is that of ci_nc_t()", {
  set.seed(77)
  t <- stats::rnorm(500, 4, 1.5)
  f <- DMAR:::.sc_nc_distance(t, 120, .025, .025)
  check <- sample(t, 25)
  direct <- vapply(check, function(v) {
    lim <- ci_nc_t(ncp = v, df = 120, conf_level = .95, verbose = FALSE)
    lim$value[2] - lim$value[1]
  }, 0)
  expect_equal(f(check), direct, tolerance = 1e-8)
})

test_that("ss_aipe_sc_ancova() returns the smallest n at which its simulated share meets the request", {
  skip_on_cran()
  plan <- ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5, divisor = "s_anova",
                            c_weights = c(1, -1), width = .5, assurance = .8,
                            assurance_method = "simulation", G = 500, seed = 3)
  n <- plan$value[1]
  expect_identical(plan$term, c("necessary_n_per_group", "assurance",
                                "assurance_attained", "G"))
  expect_gte(plan$value[3], .8)
  expect_identical(plan$value[4], 500)
  expect_equal(attr(plan, "assurance_mc_se"),
               sqrt(plan$value[3] * (1 - plan$value[3]) / 500))
  expect_identical(attr(plan, "assurance_method"), "simulation")
  expect_identical(attr(plan, "assurance_probability"), "simulation")
  expect_match(paste(utils::capture.output(print(plan)), collapse = " "),
               "Assurance method: simulation; assurance_attained is estimated from G = 500 simulated samples (Monte Carlo SE",
               fixed = TRUE)
  # The same draws, one case fewer, fall short (the planner draws them first).
  base <- withr::with_seed(3, DMAR:::.sc_ancova_base_draws(500))
  ratio <- sqrt(1 - .5^2)
  share <- function(m) mean(DMAR:::.sc_ancova_widths(
    base, m, c(1, -1), .5 / ratio, "s_anova", ratio, .025, .025) <= .5)
  expect_equal(share(n), plan$value[3])
  expect_lt(share(n - 1), .8)
  # Not below the plan without assurance; the published plan is 98.
  expect_gte(n, 97)
})

test_that("ss_aipe_sc_ancova() simulated step searches up from the smallest plan without assurance", {
  skip_on_cran()
  # The plan without assurance is 213. Before the search started at 2 per
  # group it was 215, and the simulated step returned 215
  # (assurance_attained .9651). Under the same draws 213 falls short and 214
  # attains .8861.
  args <- list(true_std_contrast = 0, c_weights = c(1, -1), width = .5,
               conf_level = .99)
  expect_identical(do.call(ss_aipe_sc_ancova, args)$value, 213)
  plan <- do.call(ss_aipe_sc_ancova, c(args, list(
    assurance = .8, assurance_method = "simulation", seed = 1)))
  expect_identical(plan$value[1], 214)
  expect_equal(plan$value[3], .8861)
  base <- withr::with_seed(1, DMAR:::.sc_ancova_base_draws(10000))
  share <- function(m) mean(DMAR:::.sc_ancova_widths(
    base, m, c(1, -1), 0, "s_ancova", 1, .005, .005) <= .5)
  expect_lt(share(213), .8)
  expect_equal(share(214), .8861)
})

## The plan of the simulated step, and its Monte Carlo standard error, are the
## same under a change of the sign of the contrast or of the weights, the
## order of the groups, the sign of the correlation, the units of the
## inputs, or the route to the ratio of the standard deviations.
test_that("ss_aipe_sc_ancova() simulated plans are invariant to sign, group order, units, and the route to the ratio", {
  skip_on_cran()  # fourteen plans
  plan <- function(...) {
    p <- ss_aipe_sc_ancova(..., width = .6, assurance = .9,
                           assurance_method = "simulation", G = 500, seed = 5)
    list(p$value, attr(p, "assurance_mc_se"))
  }
  four <- c(.5, .5, 0, -1)
  ref <- plan(true_std_contrast = .8, true_cor = .6, divisor = "s_anova",
              c_weights = four)
  same <- list(
    plan(true_std_contrast = -.8, true_cor = .6, divisor = "s_anova", c_weights = four),
    plan(true_std_contrast = -.8, true_cor = .6, divisor = "s_anova", c_weights = -four),
    plan(true_std_contrast = .8, true_cor = .6, divisor = "s_anova",
         c_weights = c(-1, 0, .5, .5)),
    plan(true_std_contrast = .8, true_cor = -.6, divisor = "s_anova", c_weights = four),
    plan(true_std_contrast = .8, true_sd_ratio = .8, divisor = "s_anova", c_weights = four),
    plan(true_contrast = .8, true_sd_anova = 1, true_sd_ancova = .8,
         divisor = "s_anova", c_weights = four),
    plan(true_contrast = 8, true_sd_anova = 10, true_sd_ancova = 8,
         divisor = "s_anova", c_weights = four))
  for (s in same) expect_identical(s, ref)
  ref <- plan(true_std_contrast = .5, c_weights = c(1, -1))
  for (s in list(plan(true_std_contrast = -.5, c_weights = c(1, -1)),
                 plan(true_std_contrast = .5, c_weights = c(-1, 1)),
                 plan(true_contrast = 5, true_sd_ancova = 10, c_weights = c(1, -1)),
                 plan(true_contrast = -1.5, true_sd_ancova = 3, c_weights = c(1, -1))))
    expect_identical(s, ref)
  # The seed makes the plan reproducible and leaves the caller's stream.
  set.seed(10)
  before <- stats::runif(1)
  set.seed(10)
  again <- plan(true_std_contrast = .5, c_weights = c(1, -1))
  expect_identical(stats::runif(1), before)
  expect_identical(again, ref)
})

test_that("ss_aipe_sc_ancova() plans attain the assurance in raw-data ANCOVA studies, both divisors", {
  skip_on_cran()  # Monte Carlo
  # 2000 raw-data studies at the simulated plan (seed 1): the share agrees
  # with assurance_attained; the s_anova published plan falls far short
  # (the s_ancova shortfall, .772 for .80, is in 13_precise_probabilities.log
  # by an exact integral, and in 03_ancova_exact_oracle.log with lm() fits).
  # All studies are usable.
  cases <- list(
    list(div = "s_ancova", delta = .5, rho = .5, cw = c(1, -1), w = .5, a = .8,
         short = FALSE),
    list(div = "s_anova", delta = .8, rho = .6, cw = c(.5, .5, 0, -1), w = .6, a = .9,
         short = TRUE))
  for (k in cases) {
    plan <- ss_aipe_sc_ancova(true_std_contrast = k$delta, true_cor = k$rho,
                              divisor = k$div, c_weights = k$cw, width = k$w,
                              assurance = k$a, assurance_method = "simulation",
                              seed = 1)
    published <- ss_aipe_sc_ancova(true_std_contrast = k$delta, true_cor = k$rho,
                                   divisor = k$div, c_weights = k$cw, width = k$w,
                                   assurance = k$a)$value[1]
    raw_share <- function(n, R = 2000) {
      w <- vapply(seq_len(R), function(r)
        .sc_ancova_raw_study(n, k$cw, k$delta, k$rho, k$div)$width, 0)
      expect_true(all(is.finite(w)))
      p <- mean(w <= k$w)
      c(p = p, se = sqrt(p * (1 - p) / R))
    }
    set.seed(2026)
    at_plan <- raw_share(plan$value[1])
    expect_lt(abs(at_plan[["p"]] - plan$value[3]) /
                sqrt(at_plan[["se"]]^2 + attr(plan, "assurance_mc_se")^2), 4)
    expect_gt(at_plan[["p"]], k$a - 3 * at_plan[["se"]])
    if (k$short) {
      at_published <- raw_share(published)
      expect_lt(at_published[["p"]], k$a - 5 * at_published[["se"]])
    }
  }
})

test_that("ss_aipe_sc_ancova() gives the plans of the first 1.1.0 builds, without a warning", {
  skip_on_cran()  # four plans, two of them simulated
  # Pinned from DMAR at 79a03795 (2026-11-13), whose logical option gave the
  # published step as 129 (divisor "s_ancova") and 98 ("s_anova") at
  # true_std_contrast = .5, true_cor = .5, c_weights = c(1, -1), width = .5,
  # and assurance = .80, and the simulated step, with G = 2000 and seed = 1,
  # as 130 (assurance_attained .860) and 103 (.802). "analytic" is the first and
  # "simulation" the second.
  plan <- function(div, ...) ss_aipe_sc_ancova(
    true_std_contrast = .5, true_cor = .5, divisor = div, c_weights = c(1, -1),
    width = .5, assurance = .8, ...)
  for (k in list(list(div = "s_ancova", a = 129, s = c(130, .86)),
                 list(div = "s_anova", a = 98, s = c(103, .802)))) {
    expect_no_warning(an <- plan(k$div))
    expect_identical(an$value[1], k$a)
    expect_no_warning(sim <- plan(k$div, assurance_method = "simulation",
                                  G = 2000, seed = 1))
    expect_identical(sim$value[c(1, 4)], c(k$s[1], 2000))
    expect_equal(sim$value[3], k$s[2], tolerance = 1e-12)
  }
})

test_that("ss_aipe_sc_ancova() refuses an assurance outside (.5, 1) by class, before any step", {
  skip_on_cran()
  for (a in list(1, 1.2, Inf, 0, -.2, .5, .3))
    for (m in c("analytic", "simulation")) for (div in c("s_ancova", "s_anova")) {
      err <- tryCatch(
        ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = if (div == "s_anova") .5,
                          divisor = div, c_weights = c(1, -1), width = .5,
                          assurance = a, assurance_method = m, G = 100, seed = 1),
        error = function(e) e)
      expect_s3_class(err, "dmar_argument_invalid")
      expect_identical(err$arg, "assurance")
      expect_match(conditionMessage(err),
                   if (a >= 1 || a <= 0) "less than unity" else "larger than 0.5")
    }
  # Through the sensitivity function as well.
  err <- tryCatch(ss_aipe_sc_ancova_sensitivity(
    true_std_contrast = .5, supposed_std_contrast = .5, c_weights = c(1, -1),
    width = .5, true_cor = .5, assurance = 1, G = 2, print_iter = FALSE, seed = 1),
    error = function(e) e)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "assurance")
})

# The ci_sc_ancova() interval at the planning value, with equal covariate
# means (the term the plan leaves out is then 0): the interval the plan
# without assurance is for.
.sc_ancova_planned_width <- function(n, delta, ratio, cw, divisor,
                                     conf_level = .95) {
  ci <- ci_sc_ancova(adj_means = cw * delta / sum(cw^2), s_anova = 1,
                     s_ancova = ratio, standardizer = divisor, c_weights = cw,
                     n_per_group = n, cov_means = rep(0, length(cw)),
                     SSwithin_x = 100, conf_level = conf_level)
  ci$value[3] - ci$value[1]
}

test_that("ss_aipe_sc_ancova() plans wide targets and near-perfect correlations with at least 2 per group", {
  skip_on_cran()
  # DMAR 1.0.x passed qt() an error df at or below 0 for these: the first
  # did not return, the others stopped with "missing value where TRUE/FALSE
  # needed" (evidence j1_sc/02). The plan is at least 2 per group (error df
  # J(n - 1) - 1 at least 1), and the ci_sc_ancova() interval at the
  # planning value is no wider than the target there (and it is the
  # smallest such n; see the next test).
  cases <- list(
    list(d = .5, cor = .5, div = "s_anova", cw = c(1, -1), w = 3.6, n = 2),
    list(d = .5, cor = NULL, div = "s_ancova", cw = c(1, -1), w = 5, n = 2),
    list(d = .5, cor = .999, div = "s_anova", cw = c(1, -1), w = .5, n = 6),
    list(d = .5, cor = .99, div = "s_anova", cw = c(.5, .5, 0, -1), w = 2, n = 2),
    list(d = 0, cor = NULL, div = "s_ancova", cw = c(1, -1), w = 50, n = 2))
  for (k in cases) {
    plan <- ss_aipe_sc_ancova(true_std_contrast = k$d, true_cor = k$cor,
                              divisor = k$div, c_weights = k$cw, width = k$w)$value
    expect_identical(plan, k$n)
    ratio <- if (k$div == "s_anova") sqrt(1 - k$cor^2) else 1
    expect_lte(.sc_ancova_planned_width(plan, k$d, ratio, k$cw, k$div), k$w)
  }
  # With assurance, both steps start there, from 2 per group, and give 4.
  p <- ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5, divisor = "s_anova",
                         c_weights = c(1, -1), width = 3.6, assurance = .8,
                         assurance_method = "simulation", G = 500, seed = 2)
  expect_identical(p$value[1], 4)
  expect_gte(p$value[3], .8)
  expect_identical(ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5,
                                     divisor = "s_anova", c_weights = c(1, -1),
                                     width = 3.6, assurance = .8)$value[1], 4)
})

test_that("ss_aipe_sc_ancova() plans the smallest n that meets the width", {
  skip_on_cran()
  # The width of the ci_sc_ancova() interval at the planning value meets
  # the target at the plan and not at one fewer, unless the plan is 2, the
  # smallest n with positive error df. DMAR 1.0.0 started its search above
  # the smallest n where the central t quantiles exceed the noncentral
  # limits, and gave 125 at a contrast of 0 and width .5, and 56 for the
  # case at conf_level = .99; at 3.6 and 5 it did not return.
  cases <- list(
    list(d = .5, cor = .5, div = "s_anova", cw = c(1, -1), w = 3.6, n = 2),
    list(d = .5, cor = NULL, div = "s_ancova", cw = c(1, -1), w = 5, n = 2),
    list(d = .5, cor = .999, div = "s_anova", cw = c(1, -1), w = .5, n = 6),
    list(d = 0, cor = NULL, div = "s_ancova", cw = c(1, -1), w = .5, n = 123),
    list(d = 2, cor = NULL, div = "s_ancova", cw = c(1, -1), w = .5, n = 185),
    list(d = -1.5, cor = .6, div = "s_anova", cw = c(.5, .5, 0, -1), w = .6, n = 54),
    list(d = .2, cor = NULL, div = "s_ancova", cw = c(1, -1), w = 1, cl = .99, n = 54))
  for (k in cases) {
    r <- if (k$div == "s_anova") sqrt(1 - k$cor^2) else 1
    cl <- if (is.null(k$cl)) .95 else k$cl
    plan <- ss_aipe_sc_ancova(true_std_contrast = k$d, true_cor = k$cor,
                              divisor = k$div, c_weights = k$cw, width = k$w,
                              conf_level = cl)$value
    expect_identical(plan, k$n)
    expect_lte(.sc_ancova_planned_width(plan, k$d, r, k$cw, k$div, cl), k$w)
    if (plan > 2)
      expect_gt(.sc_ancova_planned_width(plan - 1, k$d, r, k$cw, k$div, cl), k$w)
  }
  # The width at 2 per group of the first case from the noncentral t oracle
  # of helper-nc_t_oracle.R (error df 1): the limits for the noncentrality
  # at the observed t, times r sqrt(2 / n).
  r <- sqrt(1 - .5^2)
  t_obs <- .5 / (r * sqrt(2 / 2))
  lo <- stats::uniroot(function(m) .oracle_nct(t_obs, 1, m, lower_tail = FALSE) - .025,
                       c(-50, t_obs), tol = 1e-12)$root
  hi <- stats::uniroot(function(m) .oracle_nct(t_obs, 1, m) - .025,
                       c(t_obs, 100), tol = 1e-12)$root
  w_2 <- r * sqrt(2 / 2) * (hi - lo)
  expect_lt(w_2, 3.6)
  expect_equal(.sc_ancova_planned_width(2, .5, r, c(1, -1), "s_anova"), w_2, tolerance = 1e-7)
  expect_equal(round(w_2, 3), 3.596)
})

test_that("ss_aipe_sc_ancova() plans narrow targets promptly, the smallest n meeting the width", {
  skip_on_cran()
  # About 8 million per group (raising n by one took minutes): the
  # interval meets the width at the plan and not at one fewer, with and
  # without the published assurance step started from it.
  plan <- ss_aipe_sc_ancova(true_std_contrast = .5, c_weights = c(1, -1), width = .002)$value
  expect_identical(plan, 7923009)
  expect_lte(.sc_ancova_planned_width(plan, .5, 1, c(1, -1), "s_ancova"), .002)
  expect_gt(.sc_ancova_planned_width(plan - 1, .5, 1, c(1, -1), "s_ancova"), .002)
  pub <- ss_aipe_sc_ancova(true_std_contrast = .5, c_weights = c(1, -1), width = .002,
                           assurance = .8)$value[1]
  expect_gt(pub, plan)
})

test_that("ss_aipe_sc_ancova() refuses a width it cannot reach within 2^53 per group", {
  skip_on_cran()
  for (a in list(NULL, .8)) for (m in c("analytic", "simulation")) {
    err <- tryCatch(ss_aipe_sc_ancova(true_std_contrast = .5, c_weights = c(1, -1),
                                      width = 1e-9, assurance = a,
                                      assurance_method = m, G = 100, seed = 1),
                    error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "width")
    expect_match(conditionMessage(err), "more than 2^53 per group", fixed = TRUE)
  }
  # Plans DMAR 1.0.0 returned at once are still returned: here the smallest
  # n, 3,073,167,057 by a 512-bit computation of the interval (1.0.0 gave
  # one more, from a start above it).
  expect_identical(ss_aipe_sc_ancova(true_std_contrast = 0, true_cor = .5,
                                     c_weights = c(1, -1), width = .0001)$value[1],
                   3073167057)
  # A width whose normal bound C (z_l + z_u)^2 / width^2 is below 2^53 and
  # whose plan is above it is refused from the search (.sc_smallest_n() at
  # n_max); at a contrast of .1 it stopped with dmar_computation_failed
  # before the noncentral t was computed without rounding S.
  orig <- DMAR:::.sc_smallest_n
  reached <- FALSE
  local_mocked_bindings(
    .sc_smallest_n = function(n, meets, n_max, unreachable)
      orig(n, meets, n_max, function() {
        reached <<- TRUE
        unreachable()
      }),
    .package = "DMAR")
  width <- sqrt(2 * (2 * stats::qnorm(.975))^2 / (.9999 * 2^53))
  err <- tryCatch(ss_aipe_sc_ancova(true_std_contrast = .1, c_weights = c(1, -1),
                                    width = width), error = function(e) e)
  expect_true(reached)
  expect_s3_class(err, "dmar_argument_invalid")
  expect_identical(err$arg, "width")
  expect_match(conditionMessage(err), "more than 2^53 per group", fixed = TRUE)
})

test_that("ss_aipe_sc_ancova() returns a plan near 3e15 per group at a nonzero contrast", {
  skip_on_cran()
  # It stopped with dmar_computation_failed at 6.2e15 degrees of freedom.
  # The plan of ss_aipe_sc() for the same width has n* =
  # 3,169,203,527,072,653.0 by a 512-bit computation of the interval
  # (work/round6/evidence/m4_nct, f6_oracle.R; one degree of freedom
  # fewer moves it by less than 1); both plans rest on the difference of
  # two limits near 2e7, which doubles hold to 3.7e-9, and are within
  # 1.6e-9 of it (see Details of ?ss_aipe_sc_ancova).
  plan <- ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5,
                            c_weights = c(1, -1), width = 1e-7)$value
  expect_lte(plan, 2^53)
  expect_equal(plan, 3169203527072654, tolerance = 1e-8)
})

test_that("ss_aipe_sc_ancova() simulated plans at a wide target attain the assurance in raw-data studies", {
  skip_on_cran()  # Monte Carlo
  # At 3.6 the plan without assurance is 4 per group (error df 5). The
  # simulated plan and the share of 4000 raw-data ANCOVA studies (all
  # usable) at it, and at one fewer when that is not below 2.
  plan <- ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5, divisor = "s_anova",
                            c_weights = c(1, -1), width = 3.6, assurance = .8,
                            assurance_method = "simulation", G = 4000, seed = 36)
  n <- plan$value[1]
  share <- function(m, R = 4000) {
    w <- vapply(seq_len(R), function(r)
      .sc_ancova_raw_study(m, c(1, -1), .5, .5, "s_anova")$width, 0)
    expect_true(all(is.finite(w)))
    p <- mean(w <= 3.6)
    c(p = p, se = sqrt(p * (1 - p) / R))
  }
  set.seed(3600)
  at <- share(n)
  expect_lt(abs(at[["p"]] - plan$value[3]) /
              sqrt(at[["se"]]^2 + attr(plan, "assurance_mc_se")^2), 4)
  expect_gt(at[["p"]], .8 - 3 * at[["se"]])
})

test_that("all-zero contrast weights are refused by ss_aipe_sc_ancova() and ci_sc_ancova()", {
  # They pass the sum and positive-sum rules and define no contrast; they
  # stopped with "missing value where TRUE/FALSE needed" and with a message
  # about 'ncp'.
  e <- tryCatch(ss_aipe_sc_ancova(true_std_contrast = .5, true_cor = .5,
                                  c_weights = c(0, 0), width = .5),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "c_weights")
  expect_match(conditionMessage(e), "no nonzero weight", fixed = TRUE)
  e <- tryCatch(ci_sc_ancova(adj_means = c(1, 2), s_anova = 1, s_ancova = .8,
                             c_weights = c(0, 0), n_per_group = c(20, 20),
                             cov_means = c(0, .1), SSwithin_x = 40),
                error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "c_weights")
})

test_that("ss_aipe_sc_ancova() at a contrast of 0 takes the width from the normal quantiles past 1e15 per group", {
  # At a contrast of 0 the interval's limits are -z and z, and the width is
  # now taken from qnorm() rather than from the root of the noncentral t
  # limits, whose tolerance decided the last digits past about 4e15 per
  # group: these plans were 7,827,306,117,430,336 and 7,768,059,475,289,048.
  expect_identical(ss_aipe_sc_ancova(true_std_contrast = 0, c_weights = c(1, -1),
                                     width = 6.26595e-8)$value[1], 7827306117430337)
  expect_identical(ss_aipe_sc_ancova(true_std_contrast = 0,
                                     c_weights = c(.5, -.5, .5, -.5),
                                     width = 4.44756e-8)$value[1], 7768059475289050)
})
