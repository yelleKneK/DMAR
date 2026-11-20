test_that("default method = 'none' returns NA limits and a message", {
  expect_message(
    res <- ci_eta_squared_generalized(SS_effect = 100, SS_observed = 70,
                                      SS_error = 200),
    "No CI computed"
  )
  expect_equal(res$method, "none")
  expect_true(is.na(res$ci_lower))
  expect_true(is.na(res$ci_upper))
  expect_equal(res$eta_squared_generalized, 100 / (100 + 70 + 200),
               tolerance = 1e-12)
})

test_that("parametric method returns a CI and warns about approximation", {
  fit <- aov(len ~ supp * dose, data = ToothGrowth)
  expect_warning(
    res <- ci_eta_squared_generalized(fit, observed = "supp",
                                      method = "parametric"),
    "Parametric CI uses an approximate"
  )
  expect_equal(res$method[1], "parametric")
  # CI brackets the point estimate for each row
  for (i in seq_len(nrow(res))) {
    expect_gte(res$eta_squared_generalized[i], res$ci_lower[i])
    expect_lte(res$eta_squared_generalized[i], res$ci_upper[i])
  }
  # All limits in [0, 1]
  expect_true(all(res$ci_lower >= 0))
  expect_true(all(res$ci_upper <= 1))
})

test_that("parametric method with F/df interface requires N", {
  expect_error(
    suppressWarnings(ci_eta_squared_generalized(
      F_effect = 6, df_effect = 2,
      F_observed = 2.5, df_observed = 1, df_error = 50,
      method = "parametric"
    )),
    "requires 'N'"
  )
})

test_that("parametric F/df interface returns a valid CI when N is supplied", {
  res <- suppressWarnings(ci_eta_squared_generalized(
    F_effect = 6, df_effect = 2,
    F_observed = 2.5, df_observed = 1, df_error = 50,
    N = 54, method = "parametric"
  ))
  expect_gte(res$eta_squared_generalized, res$ci_lower)
  expect_lte(res$eta_squared_generalized, res$ci_upper)
  expect_gte(res$ci_lower, 0)
  expect_lte(res$ci_upper, 1)
})

test_that("bootstrap method requires a fitted model", {
  expect_error(
    suppressWarnings(ci_eta_squared_generalized(
      SS_effect = 100, SS_observed = 70, SS_error = 200,
      method = "bootstrap"
    )),
    "requires a fitted model"
  )
})

test_that("bootstrap method returns CI; reproducible with a seed", {
  # Two B = 1000 bootstraps back-to-back; refits the aov() inside the
  # boot() loop, slow on CRAN. The parametric CI test above already
  # checks the same bracketing and "method" attribute.
  skip_on_cran()
  fit <- aov(len ~ supp * dose, data = ToothGrowth)
  res1 <- suppressWarnings(ci_eta_squared_generalized(
    fit, observed = "supp", method = "bootstrap", B = 1000, seed = 113
  ))
  res2 <- suppressWarnings(ci_eta_squared_generalized(
    fit, observed = "supp", method = "bootstrap", B = 1000, seed = 113
  ))
  expect_equal(res1$ci_lower, res2$ci_lower)
  expect_equal(res1$ci_upper, res2$ci_upper)
  for (i in seq_len(nrow(res1))) {
    expect_gte(res1$eta_squared_generalized[i], res1$ci_lower[i])
    expect_lte(res1$eta_squared_generalized[i], res1$ci_upper[i])
    expect_gte(res1$ci_lower[i], 0)
    expect_lte(res1$ci_upper[i], 1)
  }
  expect_equal(res1$method[1], "bootstrap")
})

test_that("invalid R is rejected (minimum 1000)", {
  fit <- aov(len ~ supp * dose, data = ToothGrowth)
  expect_error(
    suppressWarnings(ci_eta_squared_generalized(
      fit, observed = "supp", method = "bootstrap", B = 999
    )),
    "'B' must be"
  )
})

test_that("method argument is matched with match.arg", {
  expect_error(
    ci_eta_squared_generalized(SS_effect = 100, SS_observed = 70,
                               SS_error = 200, method = "garbage"),
    "should be one of"
  )
})

test_that("output columns are consistent across methods", {
  # The bootstrap arm of this column-consistency check uses B = 1000.
  # Skipped on CRAN; the default and "parametric" arms still run there.
  skip_on_cran()
  cols <- c("effect", "eta_squared_generalized", "ci_lower",
            "ci_upper", "method")
  fit  <- aov(len ~ supp * dose, data = ToothGrowth)
  r1   <- suppressMessages(ci_eta_squared_generalized(fit, observed = "supp"))
  r2   <- suppressWarnings(ci_eta_squared_generalized(
    fit, observed = "supp", method = "parametric"
  ))
  r3   <- suppressWarnings(ci_eta_squared_generalized(
    fit, observed = "supp", method = "bootstrap", B = 1000, seed = 113
  ))
  for (r in list(r1, r2, r3)) {
    expect_true(all(cols %in% names(r)))
  }
})

test_that("parametric CI reports the lower-limit clamp once, in terms of generalized eta squared", {
  set.seed(113)
  d <- data.frame(
    y = rnorm(60),
    a = factor(rep(1:3, each = 20)),
    b = factor(rep(1:2, 30))
  )
  fit <- aov(y ~ a * b, data = d)
  msgs <- capture_warnings(result <- ci_eta_squared_generalized(fit, method = "parametric"))
  clamp_msgs <- msgs[grepl("lower confidence limit on generalized eta squared is 0", msgs)]
  expect_length(clamp_msgs, 1L)
  expect_false(any(grepl("prob_greater", msgs)))
  expect_true(all(result$ci_lower == 0))
  expect_true(all(result$ci_upper > 0))
})


# ---------- The size with which the limits are converted (1.1.0) ----------
#
# An effect of an aovlist fit is converted with the size of its stratum, as
# ci_eta_squared() converts it; a fit without an Error() term with N, as in
# DMAR 1.0.0. The checks below invert the noncentral F by uniroot() on pf(),
# independently of ci_nc_F(), and map the limits for lambda through
# r = lambda / size and r / (r + r_obs + 1).

.ncp_limits_pf <- function(F_value, df_1, df_2, conf_level = 0.95) {
  a <- (1 - conf_level) / 2
  gap <- function(ncp, p) stats::pf(F_value, df_1, df_2, ncp = ncp) - p
  lower <- if (stats::pf(F_value, df_1, df_2) < 1 - a) 0 else
    stats::uniroot(gap, c(0, 1e4), p = 1 - a, tol = 1e-12)$root
  upper <- stats::uniroot(gap, c(0, 1e4), p = a, tol = 1e-12)$root
  c(lower, upper)
}

.eta_g_limits_pf <- function(F_value, df_1, df_2, size, r_obs) {
  r <- .ncp_limits_pf(F_value, df_1, df_2) / size
  r / (r + r_obs + 1)
}

# The three IQ occasions of pygmalion in long form, as vignette("DMAR")
# analyzes them, with the expectancy treatment and the grade beside each
# child's rows.
.pygmalion_long <- function() {
  n <- nrow(pygmalion)
  occasions <- c("pretest", "intermediate", "end_of_study")
  data.frame(
    id        = factor(rep(seq_len(n), times = 3)),
    occasion  = factor(rep(occasions, each = n), levels = occasions),
    treatment = rep(pygmalion$treatment, times = 3),
    grade     = factor(rep(pygmalion$grade, times = 3)),
    iq        = c(pygmalion$iq_pre, pygmalion$iq_4, pygmalion$iq_8)
  )
}

# The residual sum of squares of each stratum of an aovlist, by name.
.stratum_residual_ss <- function(fit) {
  s <- summary(fit)
  vapply(s, function(x) {
    tab <- x[[1L]]
    tab[trimws(rownames(tab)) == "Residuals", "Sum Sq"]
  }, numeric(1))
}

test_that("the parametric interval of the pygmalion occasions is converted with the stratum's 620 and contains its estimate", {
  # F(2, 618) for 310 children: the occasion stratum holds 2 + 618 = 620
  # degrees of freedom. DMAR 1.0.x converted with the 930 rows and gave
  # [.0193, .0366], below the estimate of .0408.
  fit <- aov(iq ~ occasion + Error(id / occasion), data = .pygmalion_long())
  res <- suppressWarnings(ci_eta_squared_generalized(fit, method = "parametric"))
  expect_identical(res$stratum, "id:occasion")
  ss <- .stratum_residual_ss(fit)
  r_obs <- ss[["Error: id"]] / ss[["Error: id:occasion"]]
  F_occ <- summary(fit)[["Error: id:occasion"]][[1L]][1L, "F value"]
  expect_equal(c(res$ci_lower, res$ci_upper),
               .eta_g_limits_pf(F_occ, 2, 618, 620, r_obs), tolerance = 1e-6)
  expect_gt(res$eta_squared_generalized, res$ci_lower)
  expect_lt(res$eta_squared_generalized, res$ci_upper)
  expect_equal(round(c(res$eta_squared_generalized, res$ci_lower,
                       res$ci_upper), 4), c(0.0408, 0.0287, 0.0539))
  # The 1.0.x conversion with the rows, an interval below the estimate.
  rows <- .eta_g_limits_pf(F_occ, 2, 618, 930, r_obs)
  expect_equal(round(rows, 4), c(0.0193, 0.0366))
  expect_lt(rows[2], res$eta_squared_generalized)
  # The F interface reproduces the row with the stratum's size as N and the
  # children's residual given as an observed source, its F the ratio of its
  # mean square to the occasions' error mean square on 309 df.
  F_id <- (ss[["Error: id"]] / 309) / (ss[["Error: id:occasion"]] / 618)
  raw <- suppressWarnings(ci_eta_squared_generalized(
    F_effect = F_occ, df_effect = 2, F_observed = F_id, df_observed = 309,
    df_error = 618, N = 620, method = "parametric"))
  expect_equal(raw$eta_squared_generalized, res$eta_squared_generalized,
               tolerance = 1e-12)
  expect_equal(c(raw$ci_lower, raw$ci_upper),
               c(res$ci_lower, res$ci_upper), tolerance = 1e-8)
})

test_that("an effect in the subjects stratum of a mixed design is converted with the number of subjects", {
  # The expectancy treatment of pygmalion between subjects, the occasions
  # within: treatment is tested among the 310 children (1 + 308 degrees of
  # freedom), and occasion and the interaction in a stratum of
  # 2 + 2 + 616 = 620.
  fit <- aov(iq ~ treatment * occasion + Error(id / occasion),
             data = .pygmalion_long())
  expect_identical(DMAR:::.aovlist_stratum_sizes(fit),
                   c(id = 310, "id:occasion" = 620))
  res <- suppressWarnings(ci_eta_squared_generalized(fit, method = "parametric"))
  expect_identical(res$effect, c("treatment", "occasion", "treatment:occasion"))
  ss <- .stratum_residual_ss(fit)
  s <- summary(fit)
  between <- s[["Error: id"]][[1L]]
  within  <- s[["Error: id:occasion"]][[1L]]
  # Treatment: its error is the children's residual; the within-subjects
  # residual joins the denominator. Occasion and the interaction: the
  # reverse.
  lim_t <- .eta_g_limits_pf(between[1L, "F value"], 1, 308, 310,
                            ss[["Error: id:occasion"]] / ss[["Error: id"]])
  lim_o <- .eta_g_limits_pf(within[1L, "F value"], 2, 616, 620,
                            ss[["Error: id"]] / ss[["Error: id:occasion"]])
  expect_equal(c(res$ci_lower[1], res$ci_upper[1]), lim_t,
               tolerance = 1e-6)
  expect_equal(c(res$ci_lower[2], res$ci_upper[2]), lim_o,
               tolerance = 1e-6)
  for (k in seq_len(nrow(res))) {
    expect_gte(res$eta_squared_generalized[k], res$ci_lower[k])
    expect_lte(res$eta_squared_generalized[k], res$ci_upper[k])
  }
  # Treatment's upper limit, converted with the 930 rows in 1.0.x and with
  # the 310 children now.
  rows_t <- .eta_g_limits_pf(between[1L, "F value"], 1, 308, 930,
                             ss[["Error: id:occasion"]] / ss[["Error: id"]])
  expect_equal(round(c(rows_t[2], res$ci_upper[1]), 4), c(0.0152, 0.0442))
})

test_that("a measured between-subjects factor of a mixed design keeps the within-subjects residual and the interaction in its denominator", {
  set.seed(113)
  d <- data.frame(subject = factor(rep(1:20, each = 3)),
                  group   = factor(rep(c("a", "b"), each = 30)),
                  time    = factor(rep(1:3, 20)))
  d$y <- rnorm(20)[d$subject] + 0.4 * (d$group == "b") +
    0.5 * as.integer(d$time) + rnorm(60)
  fit <- aov(y ~ group * time + Error(subject / time), data = d)
  res <- suppressWarnings(ci_eta_squared_generalized(fit, observed = "group",
                                                     method = "parametric"))
  ss <- .stratum_residual_ss(fit)
  within <- summary(fit)[["Error: subject:time"]][[1L]]
  ss_gt <- within[trimws(rownames(within)) == "group:time", "Sum Sq"]
  g <- res[res$effect == "group", ]
  F_g <- summary(fit)[["Error: subject"]][[1L]][1L, "F value"]
  r_obs <- (ss_gt + ss[["Error: subject:time"]]) / ss[["Error: subject"]]
  expect_equal(c(g$ci_lower, g$ci_upper),
               .eta_g_limits_pf(F_g, 1, 18, 20, r_obs), tolerance = 1e-6)
  expect_gte(g$eta_squared_generalized, g$ci_lower)
  expect_lte(g$eta_squared_generalized, g$ci_upper)
})

test_that("fits without an Error() term give the parametric limits of DMAR 1.0.0", {
  # DMAR 1.0.0 (bcc0624), lambda / (lambda + N) with N = nobs(fit) = 310.
  pyg <- pygmalion
  pyg$grade <- factor(pyg$grade)
  one <- suppressWarnings(ci_eta_squared_generalized(
    aov(iq_gain ~ grade, data = pyg), method = "parametric"))
  expect_equal(c(one$eta_squared_generalized, one$ci_lower, one$ci_upper),
               c(0.11831301182860803, 0.048071990418499709,
                 0.17578896675246219), tolerance = 1e-10)
  fac <- suppressWarnings(ci_eta_squared_generalized(
    aov(iq_8 ~ treatment * grade, data = pyg), observed = "grade",
    method = "parametric"))
  expect_identical(fac$effect, c("treatment", "grade", "treatment:grade"))
  expect_equal(fac$eta_squared_generalized,
               c(0.020153905615505580, 0.043961812490175479,
                 0.018725635081059330), tolerance = 1e-10)
  expect_equal(fac$ci_lower,
               c(0.00097395816332176348, 0.00115921848045894395, 0),
               tolerance = 1e-10)
  expect_equal(fac$ci_upper,
               c(0.058319106180302539, 0.080679867677746195,
                 0.041291225984204323), tolerance = 1e-10)
})

test_that("the F interface accepts the size of a one-effect stratum and nothing smaller", {
  at_size <- suppressWarnings(ci_eta_squared_generalized(
    F_effect = 83.11, df_effect = 2, F_observed = 3, df_observed = 1,
    df_error = 618, N = 620, method = "parametric"))
  expect_true(all(is.finite(c(at_size$ci_lower, at_size$ci_upper))))
  expect_error(suppressWarnings(ci_eta_squared_generalized(
    F_effect = 83.11, df_effect = 2, F_observed = 3, df_observed = 1,
    df_error = 618, N = 619, method = "parametric")),
    "at least df_effect \\+ df_error")
})

test_that("an effect tested in more than one stratum stops the estimate and the interval", {
  # Twenty subjects in two groups observed three times, with four
  # observations missing: aov() divides time and group:time between the
  # subject and subject:time strata and tests each in both. DMAR 1.0.x
  # returned nine rows for the five tests, each split effect four times
  # with its share in the subject stratum, and its parametric interval
  # stopped inside the noncentral F computation.
  set.seed(11)
  d <- data.frame(subject = factor(rep(1:20, each = 3)),
                  group   = factor(rep(c("a", "b"), each = 30)),
                  time    = factor(rep(1:3, 20)))
  d$y <- rnorm(20)[d$subject] + 0.4 * (d$group == "b") +
    0.5 * as.integer(d$time) + rnorm(60)
  fit <- suppressWarnings(aov(y ~ group * time + Error(subject / time),
                              data = d[-c(2, 5, 31, 44), ]))
  tab <- DMAR:::.aovlist_effects_table(fit)
  expect_identical(tab$effect,
                   c("group", "time", "group:time", "time", "group:time"))
  msg <- paste0("tested in more than one.*: time, group:time ",
                "\\(strata subject, subject:time\\)")
  expect_error(eta_squared_generalized(fit), msg)
  expect_error(suppressMessages(ci_eta_squared_generalized(fit)), msg)
  expect_error(suppressWarnings(ci_eta_squared_generalized(
    fit, method = "parametric")), msg)
  # A covariate that varies within subjects is divided the same way.
  set.seed(5)
  dc <- data.frame(subject = factor(rep(1:12, each = 3)),
                   time    = factor(rep(1:3, 12)))
  dc$x <- rnorm(36)
  dc$y <- rnorm(12)[dc$subject] + 0.5 * dc$x + 0.3 * as.integer(dc$time) +
    rnorm(36)
  expect_error(eta_squared_generalized(
    aov(y ~ x + time + Error(subject / time), data = dc), observed = "x"),
    "varies within them: x \\(strata subject, subject:time\\)")
  # Unequal groups, with every subject observed at every time, leave each
  # effect in one stratum: dropping subject 20 leaves 10 and 9.
  fit_unequal <- aov(y ~ group * time + Error(subject / time),
                     data = droplevels(d[d$subject != "20", ]))
  res <- suppressWarnings(ci_eta_squared_generalized(fit_unequal,
                                                     method = "parametric"))
  expect_identical(res$effect, c("group", "time", "group:time"))
  expect_identical(res$stratum, c("subject", "subject:time", "subject:time"))
  ss <- .stratum_residual_ss(fit_unequal)
  F_t <- summary(fit_unequal)[["Error: subject:time"]][[1L]][1L, "F value"]
  # 19 subjects: the subject:time stratum holds 2 + 2 + 34 = 38 df.
  expect_equal(c(res$ci_lower[2], res$ci_upper[2]),
               .eta_g_limits_pf(F_t, 2, 34, 38,
                                ss[["Error: subject"]] /
                                  ss[["Error: subject:time"]]),
               tolerance = 1e-6)
})


# ---------- The summary-statistics interface (1.1.0) ----------
#
# DMAR 1.0.x stopped with "'F_effect' must be a single non-negative number"
# whenever the sums of squares came with df_effect and df_error, which the
# interval needs, because eta_squared_generalized() read the call as the F
# interface as well. The interval is built on the F the sums of squares
# imply, with SS_error the effect's error term.

.lims <- function(x) {
  c(x$eta_squared_generalized, x$ci_lower, x$ci_upper)
}

test_that("the SS interface gives the interval of the F interface and of the fitted model in a one-way design", {
  pyg <- pygmalion
  pyg$grade <- factor(pyg$grade)
  fit <- aov(iq_gain ~ grade, data = pyg)
  tab <- anova(fit)
  ss <- suppressWarnings(ci_eta_squared_generalized(
    SS_effect = tab["grade", "Sum Sq"], SS_observed = 0,
    SS_error = tab["Residuals", "Sum Sq"], df_effect = tab["grade", "Df"],
    df_error = tab["Residuals", "Df"], N = nobs(fit), method = "parametric"))
  fi <- suppressWarnings(ci_eta_squared_generalized(
    F_effect = tab["grade", "F value"], df_effect = tab["grade", "Df"],
    F_observed = 0, df_observed = 1, df_error = tab["Residuals", "Df"],
    N = nobs(fit), method = "parametric"))
  expect_equal(.lims(ss), .lims(fi), tolerance = 1e-10)
  # The values of the fitted model, as in DMAR 1.0.0.
  expect_equal(.lims(ss), c(0.11831301182860803, 0.048071990418499709,
                            0.17578896675246219), tolerance = 1e-10)
  expect_identical(ss$effect, "overall")
  expect_identical(ss$method, "parametric")
  expect_s3_class(ss, "dmar_ci_anova")
  expect_equal(attr(ss, "conf_level"), 0.95)
})

test_that("the SS interface gives the interval of the F interface and of the fitted model with a measured factor", {
  pyg <- pygmalion
  pyg$grade <- factor(pyg$grade)
  fit <- aov(iq_8 ~ treatment * grade, data = pyg)
  tab <- anova(fit)
  model <- suppressWarnings(ci_eta_squared_generalized(
    fit, observed = "grade", method = "parametric"))
  # With grade measured, the other sources in each effect's denominator;
  # every effect is tested against the residual.
  others <- list(treatment = c("grade", "treatment:grade"),
                 grade = "treatment:grade",
                 "treatment:grade" = "grade")
  for (eff in names(others)) {
    o <- others[[eff]]
    ss <- suppressWarnings(ci_eta_squared_generalized(
      SS_effect = tab[eff, "Sum Sq"], SS_observed = tab[o, "Sum Sq"],
      SS_error = tab["Residuals", "Sum Sq"], df_effect = tab[eff, "Df"],
      df_error = tab["Residuals", "Df"], N = nobs(fit), method = "parametric"))
    fi <- suppressWarnings(ci_eta_squared_generalized(
      F_effect = tab[eff, "F value"], df_effect = tab[eff, "Df"],
      F_observed = tab[o, "F value"], df_observed = tab[o, "Df"],
      df_error = tab["Residuals", "Df"], N = nobs(fit), method = "parametric"))
    expect_equal(.lims(ss), .lims(fi), tolerance = 1e-10, info = eff)
    expect_equal(.lims(ss), .lims(model[model$effect == eff, ]),
                 tolerance = 1e-10, info = eff)
    # The estimate is the one eta_squared_generalized() computes.
    expect_identical(ss$eta_squared_generalized,
                     eta_squared_generalized(
                       SS_effect = tab[eff, "Sum Sq"],
                       SS_observed = tab[o, "Sum Sq"],
                       SS_error = tab["Residuals", "Sum Sq"]
                     )$eta_squared_generalized, info = eff)
  }
})

test_that("for an effect tested within subjects, SS_error is the residual of its stratum", {
  # The pygmalion occasions: the residual of the occasion stratum is the
  # error term, and the children's residual joins the denominator through
  # SS_observed, as the fitted model computes it.
  fit <- aov(iq ~ occasion + Error(id / occasion), data = .pygmalion_long())
  model <- suppressWarnings(ci_eta_squared_generalized(fit, method = "parametric"))
  res <- .stratum_residual_ss(fit)
  ss_occ <- summary(fit)[["Error: id:occasion"]][[1L]][1L, "Sum Sq"]
  ss <- suppressWarnings(ci_eta_squared_generalized(
    SS_effect = ss_occ, SS_observed = res[["Error: id"]],
    SS_error = res[["Error: id:occasion"]], df_effect = 2, df_error = 618,
    N = 620, method = "parametric"))
  expect_equal(.lims(ss), .lims(model), tolerance = 1e-10)
  expect_equal(round(.lims(ss), 4), c(0.0408, 0.0287, 0.0539))
  # The two residuals together as SS_error give the same estimate but not
  # the occasions' F, and so another interval.
  lumped <- suppressWarnings(ci_eta_squared_generalized(
    SS_effect = ss_occ, SS_observed = 0,
    SS_error = res[["Error: id"]] + res[["Error: id:occasion"]],
    df_effect = 2, df_error = 618, N = 620, method = "parametric"))
  expect_equal(lumped$eta_squared_generalized, ss$eta_squared_generalized,
               tolerance = 1e-12)
  expect_gt(abs(lumped$ci_lower - ss$ci_lower), 0.01)
})

test_that("the SS interface reproduces each effect of a mixed design with a measured between-subjects factor", {
  # Grade, measured, between subjects and the occasions within. Grade is
  # tested against the children within grades, so its N is the 310
  # children; occasion and the interaction share the occasion stratum,
  # whose size is its total of 2 + 10 + 608 = 620 degrees of freedom. The
  # other sources of each denominator, the other stratum's residual among
  # them, go in SS_observed.
  fit <- aov(iq ~ grade * occasion + Error(id / occasion),
             data = .pygmalion_long())
  model <- suppressWarnings(ci_eta_squared_generalized(
    fit, observed = "grade", method = "parametric"))
  s <- summary(fit)
  between <- s[["Error: id"]][[1L]]
  within  <- s[["Error: id:occasion"]][[1L]]
  rownames(between) <- trimws(rownames(between))
  rownames(within)  <- trimws(rownames(within))
  expect_equal(sum(within[, "Df"]), 620)
  rows <- list(
    grade = list(tab = between, N = 310,
                 others = c(within["grade:occasion", "Sum Sq"],
                            within["Residuals", "Sum Sq"])),
    occasion = list(tab = within, N = 620,
                    others = c(between["grade", "Sum Sq"],
                               within["grade:occasion", "Sum Sq"],
                               between["Residuals", "Sum Sq"])),
    "grade:occasion" = list(tab = within, N = 620,
                            others = c(between["grade", "Sum Sq"],
                                       between["Residuals", "Sum Sq"])))
  for (eff in names(rows)) {
    r <- rows[[eff]]
    ss <- suppressWarnings(ci_eta_squared_generalized(
      SS_effect = r$tab[eff, "Sum Sq"], SS_observed = r$others,
      SS_error = r$tab["Residuals", "Sum Sq"], df_effect = r$tab[eff, "Df"],
      df_error = r$tab["Residuals", "Df"], N = r$N, method = "parametric"))
    expect_equal(.lims(ss), .lims(model[model$effect == eff, ]),
                 tolerance = 1e-10, info = eff)
  }
})

test_that("the SS interface no longer stops on the degrees of freedom its interval needs", {
  args <- list(SS_effect = 100, SS_observed = c(40, 30), SS_error = 200,
               df_effect = 2, df_error = 50, N = 60)
  point <- eta_squared_generalized(SS_effect = 100, SS_observed = c(40, 30),
                                   SS_error = 200)$eta_squared_generalized
  expect_message(none <- do.call(ci_eta_squared_generalized, args),
                 "No CI computed")
  expect_identical(none$eta_squared_generalized, point)
  expect_true(is.na(none$ci_lower) && is.na(none$ci_upper))
  msgs <- capture_warnings(par <- do.call(ci_eta_squared_generalized,
                                          c(args, method = "parametric")))
  expect_false(any(grepl("F_effect", msgs)))
  expect_identical(par$eta_squared_generalized, point)
  # The F the sums of squares imply, (100 / 2) / (200 / 50) = 12.5, with
  # F_observed = c(40, 30) / (200 / 50) on 1 df each.
  fi <- suppressWarnings(ci_eta_squared_generalized(
    F_effect = 12.5, df_effect = 2, F_observed = c(10, 7.5),
    df_observed = c(1, 1), df_error = 50, N = 60, method = "parametric"))
  expect_equal(.lims(par), .lims(fi), tolerance = 1e-10)
  expect_true(par$ci_lower > 0 && par$ci_upper < 1)
  # The bootstrap still needs a fitted model.
  expect_error(suppressWarnings(do.call(ci_eta_squared_generalized,
                                        c(args, method = "bootstrap"))),
               "requires a fitted model")
})

test_that("the SS interface refuses arguments that describe no interval", {
  ss <- function(...) {
    suppressWarnings(ci_eta_squared_generalized(
      SS_effect = 100, SS_observed = 70, SS_error = 200, ...,
      method = "parametric"))
  }
  needs <- "additionally requires 'df_effect', 'df_error', and 'N'"
  expect_error(ss(), needs)
  expect_error(ss(df_effect = 2, df_error = 50), needs)
  expect_error(ss(df_effect = 2, N = 60), needs)
  at_least <- "'N' must be a single finite number of at least df_effect \\+ df_error"
  expect_error(ss(df_effect = 2, df_error = 50, N = 51), at_least)
  expect_error(ss(df_effect = 2, df_error = 50, N = NA), at_least)
  expect_error(ss(df_effect = 2, df_error = 50, N = c(60, 70)), at_least)
  expect_error(ss(df_effect = 2, df_error = 50, N = "60"), at_least)
  # An infinite N would convert both limits to 0.
  expect_error(ss(df_effect = 2, df_error = 50, N = Inf), at_least)
  # N = df_effect + df_error, the size of a stratum that tests one effect,
  # is accepted, as the F interface accepts it.
  edge <- ss(df_effect = 2, df_error = 50, N = 52)
  expect_true(all(is.finite(c(edge$ci_lower, edge$ci_upper))))
  expect_error(ss(df_effect = 0, df_error = 50, N = 60),
               "'df_effect' must be a single finite positive number")
  expect_error(ss(df_effect = NA_real_, df_error = 50, N = 60),
               "'df_effect' must be a single finite positive number")
  expect_error(ss(df_effect = Inf, df_error = 50, N = 60),
               "'df_effect' must be a single finite positive number")
  expect_error(ss(df_effect = 2, df_error = -1, N = 60),
               "'df_error' must be a single finite positive number")
  expect_error(ss(df_effect = 2, df_error = c(50, 51), N = 60),
               "'df_error' must be a single finite positive number")
  # The checks of the estimate come first.
  expect_error(suppressWarnings(ci_eta_squared_generalized(
    SS_effect = -1, SS_observed = 70, SS_error = 200, df_effect = 2,
    df_error = 50, N = 60, method = "parametric")), "'SS_effect' must be")
  expect_error(suppressWarnings(ci_eta_squared_generalized(
    SS_effect = 100, SS_error = 200, df_effect = 2, df_error = 50, N = 60,
    method = "parametric")), "SS interface requires all of")
})

test_that("missing and infinite summary statistics are refused in both interfaces", {
  # DMAR 1.0.x stopped on most missing values with R's "missing value where
  # TRUE/FALSE needed", returned an estimate of 0 or NaN for an infinite
  # one, and gave the F interface the interval [0, 0] for an infinite N.
  ss <- function(...) {
    args <- utils::modifyList(
      list(SS_effect = 100, SS_observed = 70, SS_error = 200, df_effect = 2,
           df_error = 50, N = 53, method = "parametric"), list(...))
    suppressWarnings(do.call(ci_eta_squared_generalized, args))
  }
  expect_error(ss(SS_effect = NA_real_),
               "'SS_effect' must be a single finite non-negative number")
  expect_error(ss(SS_effect = Inf),
               "'SS_effect' must be a single finite non-negative number")
  expect_error(ss(SS_observed = c(40, NA)),
               "'SS_observed' must be a finite non-negative")
  expect_error(ss(SS_error = Inf),
               "'SS_error' must be a single finite positive number")
  expect_error(ss(SS_error = Inf, method = "none"),
               "'SS_error' must be a single finite positive number")
  fi <- function(...) {
    args <- utils::modifyList(
      list(F_effect = 25, df_effect = 2, F_observed = 85, df_observed = 1,
           df_error = 50, N = 53, method = "parametric"), list(...))
    suppressWarnings(do.call(ci_eta_squared_generalized, args))
  }
  expect_error(fi(F_effect = Inf),
               "'F_effect' must be a single finite non-negative number")
  # DMAR 1.0.x accepted an infinite df_effect for the estimate and then
  # reported it as an N too small.
  expect_error(fi(df_effect = Inf),
               "'df_effect' must be a single finite positive number")
  expect_error(fi(F_observed = NA_real_),
               "'F_observed' must be a finite non-negative")
  expect_error(fi(N = Inf),
               "'N' must be a single finite number of at least df_effect \\+ df_error")
})

test_that("the sums of squares and the F interface together are answered when they agree and refused when they do not", {
  # The same estimate, 100 / 370, from the sums of squares, whose F is
  # (100 / 2) / (200 / 50) = 12.5, and from F = 25 with F_observed = 85:
  # DMAR 1.0.x accepted both and built the interval from F = 25.
  both <- list(SS_effect = 100, SS_observed = 70, SS_error = 200,
               F_effect = 25, df_effect = 2, F_observed = 85, df_observed = 1,
               df_error = 50, N = 53)
  expect_equal(eta_squared_generalized(SS_effect = 100, SS_observed = 70,
                                       SS_error = 200, F_effect = 25,
                                       df_effect = 2, F_observed = 85,
                                       df_observed = 1, df_error = 50
                                       )$eta_squared_generalized,
               100 / 370, tolerance = 1e-12)
  from_F <- suppressWarnings(ci_eta_squared_generalized(
    F_effect = 25, df_effect = 2, F_observed = 85, df_observed = 1,
    df_error = 50, N = 53, method = "parametric"))
  from_SS <- suppressWarnings(ci_eta_squared_generalized(
    SS_effect = 100, SS_observed = 70, SS_error = 200, df_effect = 2,
    df_error = 50, N = 53, method = "parametric"))
  expect_equal(from_F$eta_squared_generalized, from_SS$eta_squared_generalized,
               tolerance = 1e-12)
  # The intervals NEWS reports for the two.
  expect_equal(round(c(from_SS$ci_lower, from_SS$ci_upper), 3),
               c(0.089, 0.413))
  expect_equal(round(c(from_F$ci_lower, from_F$ci_upper), 3),
               c(0.127, 0.381))
  refusal <- "describe different tests of the effect"
  expect_error(suppressWarnings(do.call(ci_eta_squared_generalized,
                                        c(both, method = "parametric"))),
               refusal)
  expect_error(do.call(ci_eta_squared_generalized, both), refusal)
  # Part of the F interface with the sums of squares cannot be checked.
  expect_error(ci_eta_squared_generalized(SS_effect = 100, SS_observed = 70,
                                          SS_error = 200, df_observed = 1),
               "in full, with df_effect and df_error")
  # The F-values these sums of squares imply, F = 12.5 for the effect and
  # 17.5 for the measured source on 1 df (ratio 70 / 200 = 17.5 / 50),
  # describe the same test: the call is answered, as in DMAR 1.0.x, with
  # the sums of squares' interval.
  same <- list(SS_effect = 100, SS_observed = 70, SS_error = 200,
               F_effect = 12.5, df_effect = 2, F_observed = 17.5,
               df_observed = 1, df_error = 50, N = 53, method = "parametric")
  agreed <- suppressWarnings(do.call(ci_eta_squared_generalized, same))
  expect_equal(c(agreed$ci_lower, agreed$ci_upper),
               c(from_SS$ci_lower, from_SS$ci_upper), tolerance = 1e-10)
  same$F_effect <- 12.6
  expect_error(suppressWarnings(do.call(ci_eta_squared_generalized, same)))
  # The pygmalion one-way from one anova() table, both interfaces given:
  # DMAR 1.0.0 answered it, 0.118 [0.048, 0.176].
  tab <- anova(lm(iq_gain ~ factor(grade), data = pygmalion))
  pyg <- suppressWarnings(ci_eta_squared_generalized(
    SS_effect = tab[1, "Sum Sq"], SS_observed = 0,
    SS_error = tab[2, "Sum Sq"],
    F_effect = tab[1, "F value"], F_observed = 0, df_observed = 1,
    df_effect = tab[1, "Df"],
    df_error = tab[2, "Df"], N = 310, method = "parametric"))
  expect_equal(round(c(pyg$eta_squared_generalized, pyg$ci_lower,
                       pyg$ci_upper), 3), c(0.118, 0.048, 0.176))
  # A fitted model sets the summary statistics aside, as before.
  fit <- aov(iq_gain ~ factor(grade), data = pygmalion)
  expect_warning(
    with_fit <- suppressMessages(do.call(ci_eta_squared_generalized,
                                         c(list(fit), both))),
    "arguments are ignored")
  expect_identical(with_fit$eta_squared_generalized,
                   eta_squared_generalized(fit)$eta_squared_generalized)
})
