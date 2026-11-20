# Tests for within-subjects (aovlist / multi-stratum) support across the
# eta_squared family. Data is simulated inline with set.seed(113) so the
# package needs no external dataset dependency.

# Helper: build a balanced one-way within-subjects design (subject x time)
.make_rm_data <- function(n = 20, sd_subject = 1.5, sd_within = 1.2,
                          time_effect = 0.7) {
  set.seed(113)
  subject <- factor(rep(seq_len(n), each = 3))
  time    <- factor(rep(c("Pre", "Mid", "Post"), n), levels = c("Pre", "Mid", "Post"))
  y <- rnorm(n, sd = sd_subject)[as.integer(subject)] +
       time_effect * as.integer(time) +
       rnorm(n * 3, sd = sd_within)
  data.frame(subject = subject, time = time, y = y)
}

# Helper: build a 2 (group: between) x 3 (time: within) mixed design
.make_mixed_data <- function(n_per_group = 10) {
  set.seed(113)
  n <- n_per_group * 2
  subject <- factor(rep(seq_len(n), each = 3))
  group   <- factor(rep(c("Treatment", "Control"), each = 3 * n_per_group))
  time    <- factor(rep(c("Pre", "Mid", "Post"), n), levels = c("Pre", "Mid", "Post"))
  # Subject random intercept + main effects + interaction + noise
  y <- rnorm(n, sd = 1)[as.integer(subject)] +
       0.4 * (group == "Treatment") +
       0.5 * as.integer(time) +
       0.3 * (group == "Treatment") * as.integer(time) +
       rnorm(n * 3, sd = 1)
  data.frame(subject = subject, group = group, time = time, y = y)
}


# ---------- eta_squared on aovlist ----------

test_that("eta_squared(aovlist) returns the within-subjects effect with its stratum", {
  fit <- aov(y ~ time + Error(subject/time), data = .make_rm_data())
  res <- eta_squared(fit)
  expect_true(inherits(fit, "aovlist"))
  expect_true("stratum" %in% names(res))
  expect_equal(nrow(res), 1L)
  expect_equal(res$effect, "time")
  # Within-subjects effect's stratum should include "time"
  expect_match(res$stratum[1], "time")
})

test_that("eta_squared(aovlist) point estimate matches the F/df formula per stratum", {
  rm_data <- .make_rm_data()
  fit <- aov(y ~ time + Error(subject/time), data = rm_data)
  res <- eta_squared(fit)
  s <- summary(fit)
  # Find the within-subjects stratum
  inner <- s[[grep("time", names(s))]]
  if (is.list(inner) && !is.data.frame(inner)) inner <- inner[[1L]]
  F_val <- inner["time", "F value"]
  df_e  <- inner["time", "Df"]
  df_r  <- inner["Residuals", "Df"]
  expected <- (df_e * F_val) / (df_e * F_val + df_r)
  expect_equal(res$eta_squared_partial, expected, tolerance = 1e-9)
})

test_that("estimand = \"partial\" on an aovlist gives the partial column of estimand = \"all\"", {
  fit <- aov(y ~ time + Error(subject/time), data = .make_rm_data())
  expect_identical(eta_squared(fit)$eta_squared_partial,
                   eta_squared(fit, estimand = "partial")$eta_squared_partial)
})

test_that("eta_squared(aovlist) on a mixed design returns effects from both strata", {
  fit <- aov(y ~ group * time + Error(subject/time),
             data = .make_mixed_data())
  res <- eta_squared(fit)
  # Effects should include 'group' (between) and 'time' and 'group:time' (within)
  expect_setequal(res$effect, c("group", "time", "group:time"))
  # group's stratum should NOT contain "time"
  stratum_group <- res$stratum[res$effect == "group"]
  expect_false(grepl("time", stratum_group))
  # time and group:time should be in the within stratum
  expect_match(res$stratum[res$effect == "time"], "time")
  expect_match(res$stratum[res$effect == "group:time"], "time")
})


# ---------- ci_eta_squared on aovlist ----------

test_that("ci_eta_squared(aovlist) brackets the point estimate", {
  fit <- aov(y ~ time + Error(subject/time), data = .make_rm_data())
  res <- ci_eta_squared(fit)
  expect_true(all(res$eta_squared_partial_ci_lower <= res$eta_squared_partial))
  expect_true(all(res$eta_squared_partial_ci_upper >= res$eta_squared_partial))
  expect_true(all(res$eta_squared_partial_ci_lower >= 0))
  expect_true(all(res$eta_squared_partial_ci_upper <= 1))
})

test_that("ci_eta_squared(aovlist) reports the correct total N", {
  rm_data <- .make_rm_data()
  fit <- aov(y ~ time + Error(subject/time), data = rm_data)
  res <- ci_eta_squared(fit)
  expect_equal(res$N[1], nrow(rm_data))   # 20 subjects * 3 times = 60
})


# The noncentral F inverted by uniroot() on pf(), independently of
# ci_nc_F(), and its limits converted with M.
.limits_by_pf <- function(F_value, df_1, df_2, M, a = 0.025) {
  gap <- function(ncp, p) stats::pf(F_value, df_1, df_2, ncp = ncp) - p
  lower <- if (stats::pf(F_value, df_1, df_2) < 1 - a) 0 else
    stats::uniroot(gap, c(0, 1e4), p = 1 - a, tol = 1e-12)$root
  upper <- stats::uniroot(gap, c(0, 1e4), p = a, tol = 1e-12)$root
  c(lower, upper) / (c(lower, upper) + M)
}

test_that("ci_eta_squared(aovlist) converts each effect with its stratum's size, not the rows", {
  # Through 1.0.x every effect was converted with the 60 rows. The group
  # effect is tested among the 20 subjects, one more than the subject
  # stratum's 1 + 18 degrees of freedom; time and group:time are tested in
  # the subject:time stratum, whose 2 + 2 + 36 = 40 degrees of freedom are
  # the 20 subjects' 2 within-subject contrasts each.
  fit <- aov(y ~ group * time + Error(subject/time),
             data = .make_mixed_data())
  expect_identical(DMAR:::.aovlist_stratum_sizes(fit),
                   c(subject = 20, "subject:time" = 40))
  res <- suppressWarnings(ci_eta_squared(fit))
  expect_equal(res$N, rep(60, 3))
  M <- c(group = 20, time = 40, "group:time" = 40)[res$effect]
  for (k in seq_len(nrow(res))) {
    expect_equal(c(res$eta_squared_partial_ci_lower[k], res$eta_squared_partial_ci_upper[k]),
                 .limits_by_pf(res$F_value[k], res$df_effect[k],
                               res$df_error[k], M[[k]]),
                 tolerance = 1e-6, info = res$effect[k])
    expect_gte(res$eta_squared_partial[k], res$eta_squared_partial_ci_lower[k])
    expect_lte(res$eta_squared_partial[k], res$eta_squared_partial_ci_upper[k])
  }
  # The group effect's upper limit, converted with the 60 rows in 1.0.x
  # and with the 20 subjects now.
  g <- res[res$effect == "group", ]
  expect_equal(round(.limits_by_pf(g$F_value, 1, 18, 60)[2], 3), 0.159)
  expect_equal(round(g$eta_squared_partial_ci_upper, 3), 0.362)
})

test_that("two within-subjects factors: one stratum per effect, or one stratum for all", {
  set.seed(113)
  d <- expand.grid(A = factor(1:2), B = factor(1:3), id = factor(1:12))
  d$y <- rnorm(12)[d$id] + 0.5 * as.integer(d$A) +
    0.3 * as.integer(d$B) + rnorm(nrow(d))
  # Error(id / (A * B)): each effect has its own stratum of 12, 24, and 24
  # degrees of freedom (12 subjects times 1, 2, and 2 contrasts).
  sep <- aov(y ~ A * B + Error(id / (A * B)), data = d)
  expect_identical(DMAR:::.aovlist_stratum_sizes(sep),
                   c(id = 12, "id:A" = 12, "id:B" = 24, "id:A:B" = 24))
  # Error(id): the three effects share the Within stratum, whose
  # 1 + 2 + 2 + 55 = 60 degrees of freedom are 12 subjects times 5.
  pooled <- aov(y ~ A * B + Error(id), data = d)
  expect_identical(DMAR:::.aovlist_stratum_sizes(pooled),
                   c(id = 12, Within = 60))
  for (fit in list(sep, pooled)) {
    res <- suppressWarnings(ci_eta_squared(fit, estimand = "partial"))
    M <- DMAR:::.aovlist_stratum_sizes(fit)[res$stratum]
    for (k in seq_len(nrow(res))) {
      expect_equal(c(res$eta_squared_partial_ci_lower[k], res$eta_squared_partial_ci_upper[k]),
                   .limits_by_pf(res$F_value[k], res$df_effect[k],
                                 res$df_error[k], M[[k]]),
                   tolerance = 1e-6, info = res$effect[k])
    }
  }
})

test_that("the converted interval of a within-subjects effect is for the value its estimate converges to", {
  skip_on_cran()
  # 3000 subjects, three occasions, occasion effects -0.25, 0, 0.25, and
  # error variance 1 after subject intercepts of SD 1.5. The sample partial
  # eta squared converges to 0.125 / (0.125 + 2 * 1) = 0.0588; converted
  # with the 9000 rows the interval is for 0.125 / (0.125 + 3 * 1), 0.04.
  set.seed(113)
  n <- 3000
  alpha <- c(-0.25, 0, 0.25)
  d <- data.frame(id = factor(rep(seq_len(n), times = 3)),
                  occasion = factor(rep(1:3, each = n)))
  d$y <- rnorm(n, sd = 1.5)[as.integer(d$id)] +
    alpha[as.integer(d$occasion)] + rnorm(3 * n)
  res <- ci_eta_squared(estimand = "partial", aov(y ~ occasion + Error(id / occasion),
                                    data = d))
  target <- sum(alpha^2) / (sum(alpha^2) + 2)
  expect_gt(target, res$eta_squared_partial_ci_lower)
  expect_lt(target, res$eta_squared_partial_ci_upper)
  rows <- ci_eta_squared(estimand = "partial", F_value = res$F_value, df_effect = 2,
                                 df_error = res$df_error, N = 3 * n)
  expect_lt(rows$eta_squared_partial_ci_upper, target)
})

# ---------- eta_squared_generalized on aovlist ----------

test_that("eta_squared_generalized(aovlist) sums all stratum residuals in denominator", {
  rm_data <- .make_rm_data()
  fit <- aov(y ~ time + Error(subject/time), data = rm_data)
  res <- eta_squared_generalized(fit)
  s <- summary(fit)
  # Compute expected: SS_time / (SS_time + SS_subject_residual + SS_subject:time_residual)
  inner_subject <- s[["Error: subject"]]
  if (is.list(inner_subject) && !is.data.frame(inner_subject))
    inner_subject <- inner_subject[[1L]]
  inner_within <- s[["Error: subject:time"]]
  if (is.list(inner_within) && !is.data.frame(inner_within))
    inner_within <- inner_within[[1L]]
  ss_subject_resid <- inner_subject["Residuals", "Sum Sq"]
  ss_within_resid  <- inner_within["Residuals", "Sum Sq"]
  ss_time          <- inner_within["time", "Sum Sq"]
  expected <- ss_time / (ss_time + ss_subject_resid + ss_within_resid)
  expect_equal(res$eta_squared_generalized, expected, tolerance = 1e-9)
})

test_that("eta_squared_generalized(aovlist) honors `observed` argument", {
  fit <- aov(y ~ group * time + Error(subject/time),
             data = .make_mixed_data())
  res_no_obs <- eta_squared_generalized(fit)
  res_obs    <- eta_squared_generalized(fit, observed = "group")
  # Treating 'group' as measured should put SS_group in the denominator for
  # 'time' and 'group:time', shrinking their eta_g values.
  for (eff in c("time", "group:time")) {
    expect_lt(
      res_obs$eta_squared_generalized[res_obs$effect == eff],
      res_no_obs$eta_squared_generalized[res_no_obs$effect == eff]
    )
  }
  # For the focal 'group' effect the denominator is the measured sources
  # plus the error strata (Olejnik & Algina, 2003, Eq. 5, delta = 0):
  # group's own SS leads it, and group:time joins it as a measured source
  # because listing 'group' marks every interaction containing it.
  eff_tbl <- DMAR:::.aovlist_effects_table(fit)
  ss <- setNames(eff_tbl$SS_effect, eff_tbl$effect)
  ss_err_total <- DMAR:::.aovlist_total_residual_ss(fit)
  expect_equal(
    res_obs$eta_squared_generalized[res_obs$effect == "group"],
    unname(ss[["group"]] / (ss[["group"]] + ss[["group:time"]] + ss_err_total)),
    tolerance = 1e-12
  )
})

test_that("eta_squared_generalized(aovlist) errors on unknown observed factor", {
  fit <- aov(y ~ time + Error(subject/time), data = .make_rm_data())
  expect_error(eta_squared_generalized(fit, observed = "not_a_factor"),
               "not in the aovlist")
})


# ---------- ci_eta_squared_generalized on aovlist ----------

test_that("ci_eta_squared_generalized(aovlist, parametric) returns a valid CI", {
  fit <- aov(y ~ time + Error(subject/time), data = .make_rm_data())
  res <- suppressWarnings(
    ci_eta_squared_generalized(fit, method = "parametric")
  )
  for (i in seq_len(nrow(res))) {
    expect_gte(res$eta_squared_generalized[i], res$ci_lower[i])
    expect_lte(res$eta_squared_generalized[i], res$ci_upper[i])
    expect_gte(res$ci_lower[i], 0)
    expect_lte(res$ci_upper[i], 1)
  }
  expect_equal(res$method[1], "parametric")
})

test_that("ci_eta_squared_generalized(aovlist, bootstrap) errors with informative message", {
  fit <- aov(y ~ time + Error(subject/time), data = .make_rm_data())
  expect_error(
    suppressWarnings(ci_eta_squared_generalized(
      fit, method = "bootstrap", B = 1000
    )),
    "not yet supported for aovlist"
  )
})

test_that("ci_eta_squared_generalized(aovlist) preserves stratum column from point estimate", {
  fit <- aov(y ~ group * time + Error(subject/time),
             data = .make_mixed_data())
  res <- suppressWarnings(
    ci_eta_squared_generalized(fit, observed = "group", method = "parametric")
  )
  expect_true("stratum" %in% names(res))
  expect_setequal(res$effect, c("group", "time", "group:time"))
})
