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


# ---------- eta_squared / eta_squared_partial on aovlist ----------

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
  expect_equal(res$eta_squared, expected, tolerance = 1e-9)
})

test_that("eta_squared_partial(aovlist) gives the same numbers as eta_squared(aovlist)", {
  fit <- aov(y ~ time + Error(subject/time), data = .make_rm_data())
  expect_equal(eta_squared(fit)$eta_squared,
               eta_squared_partial(fit)$eta_squared_partial, tolerance = 1e-12)
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
  expect_true(all(res$lower_limit <= res$eta_squared))
  expect_true(all(res$upper_limit >= res$eta_squared))
  expect_true(all(res$lower_limit >= 0))
  expect_true(all(res$upper_limit <= 1))
})

test_that("ci_eta_squared(aovlist) reports the correct total N", {
  rm_data <- .make_rm_data()
  fit <- aov(y ~ time + Error(subject/time), data = rm_data)
  res <- ci_eta_squared(fit)
  expect_equal(res$N[1], nrow(rm_data))   # 20 subjects * 3 times = 60
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
    expect_gte(res$eta_squared_generalized[i], res$lower_limit[i])
    expect_lte(res$eta_squared_generalized[i], res$upper_limit[i])
    expect_gte(res$lower_limit[i], 0)
    expect_lte(res$upper_limit[i], 1)
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
