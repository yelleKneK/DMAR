# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# Fixtures. The spatial tests of DMAR's Holzinger and Swineford data carry
# real missing values (t25 and t26 were given only at Grant-White). The
# lavaan copy of the nine classic tests gets values missing at random:
# x1, x4, and x7 stay complete, and each other test is missing with a
# probability logistic in the marker of its own factor (x2, x5, x8) or in
# age (x3, x6, x9).
.fmi_spatial <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
                  "t4_lozenges", "t25_paper_form_board_r", "t26_flags")

.fmi_spatial_fit <- function(items = .fmi_spatial, ...) {
  fa_confirmatory(data = holzinger_swineford, factors = list(spatial = items),
        missing = "fiml", output = "fit", ...)
}

.fmi_hs_mar <- function(seed = 113) {
  d <- lavaan::HolzingerSwineford1939[, c(paste0("x", 1:9), "ageyr",
                                          "agemo", "school")]
  set.seed(seed)
  n <- nrow(d)
  age <- d$ageyr + d$agemo / 12
  zs <- function(v) as.numeric(scale(v))
  drive <- list(x2 = zs(d$x1), x5 = zs(d$x4), x8 = zs(d$x7),
                x3 = zs(age), x6 = zs(age), x9 = zs(age))
  for (v in names(drive)) {
    p <- stats::plogis(-1.3 + 1.2 * drive[[v]])
    d[[v]][stats::runif(n) < p] <- NA
  }
  d
}

.fmi_m3 <- "visual =~ x1 + x2 + x3
            textual =~ x4 + x5 + x6
            speed =~ x7 + x8 + x9"

# lavaan's EM settings for the unrestricted moments under which its own
# FMI is computed as precisely as DMAR's: the refits behind
# fraction_missing_information() run the EM to an absolute tolerance of
# 1e-10, and lavaan's parameterEstimates(fmi = TRUE) runs it to the fit's
# own (1e-5 by default).
.fmi_em_tight <- list(tol = 1e-10, max_iter = 10000L)

# lavaan's own FMI for the rows fraction_missing_information() reports
# (free and defined parameters, in parameter table order), refitting with
# the options that select each of Chen and Savalei's estimators. The refit
# call is evaluated where the fit was made, so its data are found.
.fmi_lavaan <- function(fit, information, envir = parent.frame()) {
  extra <- switch(information,
    hessian = list(),
    structured = list(observed.information = "h1"),
    unstructured = list(observed.information = "h1",
                        h1.information = "unstructured"))
  refit <- if (length(extra)) {
    eval(do.call(lavaan::update, c(list(fit), extra, evaluate = FALSE)),
         envir)
  } else {
    fit
  }
  pe <- lavaan::parameterEstimates(refit, fmi = TRUE, remove.def = FALSE,
                                   remove.eq = FALSE, remove.ineq = FALSE,
                                   remove.system.eq = FALSE,
                                   remove.nonfree = FALSE)
  pt <- lavaan::parTable(fit)
  pe$fmi[pt$free > 0L | pt$op == ":="]
}


test_that("Chen and Savalei's (2021) Table 6 is reproduced by all three estimators", {
  skip_if_not_installed("lavaan")
  # Their example draws 61 values per column with sample() under R's
  # sampler before 3.6.0, although they report R 3.6.1; with the current
  # sampler the draw, and so the table, differs by up to .15.
  dm <- withr::with_rng_version("3.5.0", {
    set.seed(202104)
    as.data.frame(lapply(lavaan::HolzingerSwineford1939[paste0("x", 1:9)],
                         function(x) {
                           x[sample(seq_along(x), 61)] <- NA
                           x
                         }))
  })
  fit <- lavaan::cfa("visual =~ x1 + x2 + x3
                      verbal =~ x4 + x5 + x6
                      speed =~ x7 + x8 + x9",
                     data = dm, std.lv = TRUE, missing = "ml")
  # Table 6: nine loadings, then the visual-verbal, visual-speed, and
  # verbal-speed factor correlations; delta 1, 2, 3 are the Hessian,
  # structured, and unstructured estimators.
  published <- list(
    hessian      = c(.29, .26, .28, .24, .17, .19, .27, .23, .19, .18, .23, .18),
    structured   = c(.27, .25, .28, .22, .17, .18, .27, .24, .22, .15, .24, .18),
    unstructured = c(.29, .28, .31, .24, .19, .19, .28, .26, .22, .18, .29, .20))
  for (inf in names(published)) {
    res <- fraction_missing_information(fit, information = inf)
    sel <- res$op == "=~" |
      (res$op == "~~" & grepl("^(visual|verbal|speed) ~~ (visual|verbal|speed)$",
                              res$term))
    expect_identical(sum(sel), 12L)
    expect_equal(round(res$fmi[sel], 2), published[[inf]], tolerance = 1e-12,
                 info = inf)
  }
  res <- fraction_missing_information(fit)
  sel <- res$op == "=~" |
    (res$op == "~~" & grepl("^(visual|verbal|speed) ~~ (visual|verbal|speed)$",
                            res$term))
  expect_equal(res$fmi[sel],
               c(.2920103, .2815117, .3106551, .2413521, .1863283, .1863431,
                 .2838877, .2611273, .2249327, .1815048, .2862098, .1978106),
               tolerance = 1e-6)
})

test_that("the FMI equals lavaan's for every information, one group and two", {
  skip_if_not_installed("lavaan")
  # The fits run lavaan's EM as far as the refits run theirs, so that
  # lavaan's own FMI is computed from the same unrestricted moments.
  # lavaan's nine tests have standard deviations from 1.0 to 1.3, so the
  # refits keep their metric and the FMIs are lavaan's to rounding error.
  # DMAR's spatial tests have standard deviations from 2.8 to 9.0, so the
  # refits divide them by 4 or 8, and lavaan computes its FMI in the raw
  # metric: the two agree to about 1e-9.
  fit_spatial <- .fmi_spatial_fit(em.h1.args = .fmi_em_tight)
  d <- .fmi_hs_mar()
  fit_configural <- lavaan::cfa(.fmi_m3, data = d, std.lv = TRUE,
                                missing = "ml", group = "school",
                                em.h1.args = .fmi_em_tight)
  fit_metric <- lavaan::cfa(.fmi_m3, data = d, std.lv = TRUE,
                            missing = "ml", group = "school",
                            group.equal = "loadings",
                            em.h1.args = .fmi_em_tight)
  for (fit in list(fit_configural, fit_metric)) {
    for (inf in c("unstructured", "structured", "hessian")) {
      res <- fraction_missing_information(fit, information = inf)
      expect_equal(res$fmi, .fmi_lavaan(fit, inf), tolerance = 1e-10,
                   info = inf)
    }
  }
  for (inf in c("unstructured", "structured", "hessian")) {
    res <- fraction_missing_information(fit_spatial, information = inf)
    expect_lt(max(abs(res$fmi - .fmi_lavaan(fit_spatial, inf))), 1e-8)
  }
  # information = "hessian" is lavaan's default FMI, with no refit.
  pe <- lavaan::parameterEstimates(fit_spatial, fmi = TRUE,
                                   remove.nonfree = TRUE)
  res <- fraction_missing_information(fit_spatial, information = "hessian")
  expect_lt(max(abs(res$fmi - pe$fmi)), 1e-8)
})

test_that("the unrestricted moments are estimated to a tolerance of 1e-10", {
  skip_if_not_installed("lavaan")
  # lavaan's EM stops when no mean or covariance changes by more than an
  # absolute 1e-5, its default, and EM converges slowly where the
  # unrestricted model's fraction of missing information is high. The
  # refits ran it to the source fit's tolerance, and the FMIs moved with
  # it: with x5 missing at random on x4 and a quarter of x9 missing, the
  # Hessian FMIs of the nine tests moved by 1.9e-5, and the unstructured
  # ones by 5.9e-6, when the EM was run to 1e-14 instead. They now come
  # from moments run to 1e-10 whatever the fit's tolerance, and agree with
  # lavaan's own FMI when lavaan runs its EM as far.
  hs <- lavaan::HolzingerSwineford1939
  set.seed(11)
  hs$x5[stats::runif(301) < stats::plogis(-1 + hs$x4)] <- NA
  hs$x9[stats::runif(301) < 0.25] <- NA
  d <- hs[paste0("x", 1:9)]
  fit <- lavaan::cfa(.fmi_m3, data = d, missing = "ml")
  fit_tight <- lavaan::cfa(.fmi_m3, data = d, missing = "ml",
                           em.h1.args = list(tol = 1e-14, max_iter = 10000L))
  expect_equal(lavaan::coef(fit), lavaan::coef(fit_tight), tolerance = 1e-6)
  for (inf in c("unstructured", "hessian")) {
    res <- fraction_missing_information(fit, information = inf)
    ref <- fraction_missing_information(fit_tight, information = inf)
    expect_lt(max(abs(res$fmi - ref$fmi)), 1e-8)
  }
  # lavaan's own FMI, with its EM run to 1e-10.
  res <- fraction_missing_information(fit, information = "hessian")
  fit_10 <- lavaan::cfa(.fmi_m3, data = d, missing = "ml",
                        em.h1.args = .fmi_em_tight)
  own_10 <- lavaan::parameterEstimates(fit_10, fmi = TRUE,
                                       remove.nonfree = TRUE)
  expect_lt(max(abs(res$fmi - own_10$fmi)), 1e-8)
})

test_that("the spatial factor of DMAR's Holzinger and Swineford data matches its oracle values", {
  skip_if_not_installed("lavaan")
  fit <- .fmi_spatial_fit()
  oracle <- data.frame(
    term = c("spatial =~ t1_visual_perception",
             "spatial =~ t25_paper_form_board_r", "spatial =~ t26_flags",
             "t25_paper_form_board_r ~~ t25_paper_form_board_r",
             "t25_paper_form_board_r ~1", "omega_spatial"),
    unstructured = c(0.1125059265, 0.5116258498, 0.5061975406, 0.5165527874,
                     0.4436848778, 0.2442136260),
    structured = c(0.1114128510, 0.5230155033, 0.5221441497, 0.5288508466,
                   0.4586661598, 0.2629622705),
    hessian = c(0.1279294122, 0.5127424158, 0.5128900313, 0.5238948182,
                0.4579232555, 0.2462885730),
    stringsAsFactors = FALSE)
  for (inf in c("unstructured", "structured", "hessian")) {
    res <- fraction_missing_information(fit, information = inf)
    expect_equal(res$fmi[match(oracle$term, res$term)], oracle[[inf]],
                 tolerance = 1e-6, info = inf)
  }
  res <- fraction_missing_information(fit)
  defined <- res[res$op == ":=", ]
  expect_identical(defined$term, c("loading_sum_spatial", "error_sum_spatial",
                                   "omega_spatial", "ave_spatial",
                                   "H_spatial"))
  expect_equal(defined$fmi, c(0.1837593448, 0.3479011082, 0.2442136260,
                              0.2213079569, 0.2163195458), tolerance = 1e-6)
  omega <- res[res$term == "omega_spatial", ]
  expect_equal(omega$estimate, 0.7607995657, tolerance = 1e-6)
  expect_equal(omega$se_observed, 0.0255261246, tolerance = 1e-6)
  expect_equal(omega$se_complete, 0.0221913853, tolerance = 1e-6)
  expect_equal(omega$N_effective, 227.491699, tolerance = 1e-6)
  expect_equal(omega$N_target, 398.260686, tolerance = 1e-6)
  expect_equal(omega$width_inflation, 1.15027180, tolerance = 1e-6)
  # The four tests every pupil took have fully observed intercepts.
  full <- res$term %in% paste(.fmi_spatial[1:4], "~1")
  expect_lt(max(abs(res$fmi[full])), 1e-12)
  expect_equal(res$N, rep(301, nrow(res)))
})

test_that("with complete data every FMI is zero under all three estimators", {
  skip_if_not_installed("lavaan")
  fit <- .fmi_spatial_fit(items = .fmi_spatial[1:4])
  for (inf in c("unstructured", "structured", "hessian")) {
    res <- fraction_missing_information(fit, information = inf)
    expect_lt(max(abs(res$fmi)), 1e-8)
  }
})

test_that("the unstructured estimator is not negative across small MAR samples", {
  skip_if_not_installed("lavaan")
  set.seed(20261101)
  lam <- c(.8, .7, .7, .6, .6, .5)
  Sigma <- tcrossprod(lam) + diag(1 - lam^2)
  model <- paste("f =~", paste0("y", 1:6, collapse = " + "))
  mins <- numeric(0)
  for (r in 1:20) {
    N <- c(30, 50, 80)[(r - 1) %% 3 + 1]
    Y <- as.data.frame(matrix(stats::rnorm(N * 6), N) %*% chol(Sigma))
    names(Y) <- paste0("y", 1:6)
    p_miss <- stats::plogis(stats::qlogis(.35) + 1.5 * as.numeric(scale(Y$y1)))
    for (j in 3:6) Y[[j]][stats::runif(N) < p_miss] <- NA
    # Small samples can produce a Heywood case or a nonconvergence; such a
    # draw has no FMI to check and is skipped.
    f <- tryCatch(suppressWarnings(
      lavaan::cfa(model, data = Y, std.lv = TRUE, missing = "ml")),
      error = function(e) NULL)
    if (is.null(f) || !lavaan::lavInspect(f, "converged")) next
    res <- tryCatch(suppressWarnings(fraction_missing_information(f)),
      dmar_fmi_saturated_singular = identity)
    if (inherits(res, "dmar_fmi_saturated_singular")) {
      expect_false(is.null(DMAR:::.fmi_unbounded_subset(Y)))
      next
    }
    mins <- c(mins, min(res$fmi, na.rm = TRUE))
  }
  expect_gt(length(mins), 10L)
  expect_gte(min(mins), -1e-10)
})

test_that("under MAR the FMI can exceed the proportion of values missing", {
  skip_if_not_installed("lavaan")
  d <- .fmi_hs_mar()
  fit <- lavaan::cfa(.fmi_m3, data = d, std.lv = TRUE, missing = "ml")
  res <- fraction_missing_information(fit)
  x2 <- res$fmi[res$term == "visual =~ x2"]
  expect_equal(x2, 0.3851317344, tolerance = 1e-6)
  expect_gt(x2, mean(is.na(d$x2)))
  # The Hessian estimator breaks down for this misspecified model.
  hess <- fraction_missing_information(fit, information = "hessian")
  expect_equal(hess$fmi[hess$term == "speed =~ x7"], -2.200404,
               tolerance = 1e-5)
})

test_that("an MLR fit gets the normal theory FMI of the ML fit, and says so", {
  skip_if_not_installed("lavaan")
  d <- .fmi_hs_mar()
  fit_ml <- lavaan::cfa(.fmi_m3, data = d, std.lv = TRUE, missing = "ml")
  fit_mlr <- lavaan::cfa(.fmi_m3, data = d, std.lv = TRUE, missing = "ml",
                         estimator = "MLR")
  ml <- fraction_missing_information(fit_ml)
  mlr <- fraction_missing_information(fit_mlr)
  expect_equal(mlr$fmi, ml$fmi, tolerance = 1e-12)
  expect_equal(mlr$se_observed, ml$se_observed, tolerance = 1e-12)
  expect_identical(attr(mlr, "source_estimator"), "MLR")
  expect_identical(attr(mlr, "source_se"), "robust.huber.white")
  expect_match(attr(mlr, "basis"), "normal theory")
  expect_output(print(mlr), "MLR with robust.huber.white standard errors")
  expect_identical(attr(ml, "source_se"), "standard")
})

test_that("a source fit made without the unrestricted model gets the same FMI", {
  skip_if_not_installed("lavaan")
  # A fit made with h1 = FALSE carries no EM moments, and lavInspect()
  # returns other sample statistics in their place; the complete data side
  # must use the EM moments all the same (a first version read them from
  # the source fit and was off by .21 here).
  d <- .fmi_hs_mar()
  fit <- lavaan::cfa(.fmi_m3, data = d, std.lv = TRUE, missing = "ml")
  fit_no_h1 <- lavaan::cfa(.fmi_m3, data = d, std.lv = TRUE, missing = "ml",
                           h1 = FALSE, baseline = FALSE)
  for (inf in c("unstructured", "hessian")) {
    expect_equal(fraction_missing_information(fit_no_h1, information = inf)$fmi,
                 fraction_missing_information(fit, information = inf)$fmi,
                 tolerance = 1e-12, info = inf)
  }
})

test_that("a fit with standard errors robust to clustering gets the ML fit's FMI", {
  skip_if_not_installed("lavaan")
  d <- .fmi_hs_mar()
  d$clus <- rep(1:40, length.out = nrow(d))
  fit_cl <- lavaan::cfa(.fmi_m3, data = d, std.lv = TRUE, missing = "ml",
                        cluster = "clus")
  fit_ml <- lavaan::cfa(.fmi_m3, data = d, std.lv = TRUE, missing = "ml")
  cl <- fraction_missing_information(fit_cl)
  expect_identical(attr(cl, "source_se"), "robust.cluster")
  expect_equal(cl$fmi, fraction_missing_information(fit_ml)$fmi,
               tolerance = 1e-12)
})

test_that("a regression keeps the cases missing a predictor only under ml.x", {
  skip_if_not_installed("lavaan")
  reg <- .fmi_hs_mar()
  reg$age <- reg$ageyr + reg$agemo / 12
  set.seed(2)
  reg$age[sample.int(nrow(reg), 30)] <- NA
  # Under missing = "ml" with fixed.x = TRUE lavaan drops the 30 cases
  # missing age and says so; that notice is part of the fixture.
  fit_ml <- suppressWarnings(lavaan::sem("x5 ~ x4 + age", data = reg,
                                         missing = "ml", fixed.x = TRUE,
                                         em.h1.args = .fmi_em_tight))
  fit_mlx <- lavaan::sem("x5 ~ x4 + age", data = reg, missing = "ml.x",
                         fixed.x = TRUE, em.h1.args = .fmi_em_tight)
  res_ml <- fraction_missing_information(fit_ml)
  res_mlx <- fraction_missing_information(fit_mlx)
  expect_equal(res_ml$N, rep(271, nrow(res_ml)))
  expect_identical(attr(res_ml, "N_excluded"), 30)
  expect_equal(res_mlx$N, rep(301, nrow(res_mlx)))
  expect_identical(attr(res_mlx, "source_missing"), "ml.x")
  expect_equal(res_ml$fmi,
               suppressWarnings(.fmi_lavaan(fit_ml, "unstructured")),
               tolerance = 1e-10)
  expect_equal(res_mlx$fmi, .fmi_lavaan(fit_mlx, "unstructured"),
               tolerance = 1e-10)
})

test_that("the N behind each row follows the groups and the equality constraints", {
  skip_if_not_installed("lavaan")
  d <- .fmi_hs_mar()

  # Three rows with nothing observed are set aside by the fit.
  d3 <- d
  d3[1:3, paste0("x", 1:9)] <- NA
  fit3 <- suppressWarnings(
    lavaan::cfa(.fmi_m3, data = d3, std.lv = TRUE, missing = "ml"))
  res3 <- fraction_missing_information(fit3)
  expect_identical(attr(res3, "N_excluded"), 3)
  expect_equal(res3$N, rep(298, nrow(res3)))
  expect_equal(unname(attr(res3, "N_used")),
               as.numeric(lavaan::lavInspect(fit3, "nobs")))
  expect_output(print(res3), "set aside 3 cases")

  # Configural: each group's rows carry that group's N.
  fit_c <- lavaan::cfa(.fmi_m3, data = d, std.lv = TRUE, missing = "ml",
                       group = "school")
  res_c <- fraction_missing_information(fit_c)
  expect_setequal(unique(res_c$group), c("Pasteur", "Grant-White"))
  expect_true(all(res_c$N[res_c$group == "Pasteur"] == 156))
  expect_true(all(res_c$N[res_c$group == "Grant-White"] == 145))
  expect_equal(attr(res_c, "N_used"),
               c(Pasteur = 156, `Grant-White` = 145))
  expect_output(print(res_c),
                "301 cases \\(156 in Pasteur, 145 in Grant-White\\)")

  # Metric: the loadings constrained equal across the schools rest on all
  # 301 cases and have the same FMI in both groups.
  fit_m <- lavaan::cfa(.fmi_m3, data = d, std.lv = TRUE, missing = "ml",
                       group = "school", group.equal = "loadings")
  res_m <- fraction_missing_information(fit_m)
  load <- res_m[res_m$op == "=~", ]
  expect_true(all(load$N == 301))
  expect_equal(load$fmi[load$group == "Pasteur"],
               load$fmi[load$group == "Grant-White"], tolerance = 1e-8)
  expect_true(all(res_m$N[res_m$op == "~1" & res_m$group == "Pasteur"] == 156))

  # A defined parameter reaches the groups whose labels it uses.
  fit_d <- lavaan::cfa("visual =~ c(a1, b1)*x1 + c(a2, b2)*x2 + c(a3, b3)*x3
                        textual =~ x4 + x5 + x6
                        speed =~ x7 + x8 + x9
                        def_pasteur := a1 + a2 + a3
                        def_both := a1 + b1",
                       data = d, std.lv = TRUE, missing = "ml",
                       group = "school")
  res_d <- fraction_missing_information(fit_d)
  expect_equal(res_d$N[res_d$term == "def_pasteur"], 156)
  expect_equal(res_d$N[res_d$term == "def_both"], 301)
  expect_true(is.na(res_d$group[res_d$term == "def_pasteur"]))
})

test_that("the planning columns are the stated functions of the FMI", {
  skip_if_not_installed("lavaan")
  res <- fraction_missing_information(.fmi_spatial_fit())
  expect_identical(res$N_effective, res$N * (1 - res$fmi))
  expect_identical(res$N_target, res$N / (1 - res$fmi))
  expect_identical(res$width_inflation, 1 / sqrt(1 - res$fmi))
  expect_equal(res$width_inflation, res$se_observed / res$se_complete,
               tolerance = 1e-10)
  # Savalei and Rhemtulla's (2012, p. 488) example: N = 494, FMI = .31506.
  sr <- .fmi_precision_terms(0.31506, 494)
  expect_equal(sr$N_effective, 338.3604, tolerance = 1e-6)
  expect_equal(sr$N_target, 721.2311, tolerance = 1e-6)
})

test_that("the table has one fixed shape, prints its footer, and tidies", {
  skip_if_not_installed("lavaan")
  res <- fraction_missing_information(.fmi_spatial_fit())
  cols <- c("term", "op", "group", "fmi", "fmi_ci_lower", "fmi_ci_upper",
            "width_inflation", "N", "N_effective", "N_target", "estimate",
            "se_observed", "se_complete")
  expect_identical(names(res), cols)
  expect_s3_class(res, c("dmar_fmi", "dmar_tbl", "data.frame"), exact = TRUE)
  expect_true(all(is.na(res$fmi_ci_lower)) && all(is.na(res$fmi_ci_upper)))
  expect_true(all(is.na(res$group)))
  expect_identical(attr(res, "information"), "unstructured")
  expect_identical(attr(res, "ci_method"), "none")
  expect_null(attr(res, "conf_level"))
  expect_identical(attr(res, "n_patterns"), 2)
  expect_equal(attr(res, "fmi_max"), max(res$fmi))
  expect_false(is.null(attr(res, "information_description")))

  out <- capture.output(print(res))
  expect_false(any(grepl("fmi_ci_lower|^ *op ", out)))
  expect_true(any(grepl("third estimator of Chen & Savalei, 2021\\), from 301 cases", out)))
  expect_true(any(grepl("imprecise unless N is several hundred", out)))
  # A subset no longer carries the attributes and prints as a plain table.
  expect_output(print(res[res$op == ":=", c("term", "fmi")]), "omega_spatial")

  td <- generics::tidy(res)
  expect_identical(names(td), c("term", "group", "estimate", "estimand",
                                "ci_lower", "ci_upper", "conf_level",
                                "N_effective"))
  expect_identical(td$estimate, res$fmi)
  # The estimate is the FMI, which the estimand column names, so a row's
  # sentence reports it as the FMI and not as the parameter's value: it
  # read "omega_spatial = 0.24, 90% CI [...]" for an omega of .76.
  expect_identical(td$estimand, rep("fmi", nrow(td)))
  row <- td[td$term == "omega_spatial", ]
  row$ci_lower <- 0.18
  row$ci_upper <- 0.32
  row$conf_level <- 0.90
  expect_match(results_sentence(row), "^fmi = 0\\.[0-9]{2}, 90% CI")
  # No percentile interval, so no level beside the NA limits.
  expect_true(all(is.na(td$conf_level)))
  gl <- generics::glance(res)
  expect_identical(names(gl), c("n_parameters", "fmi_max", "fmi_median", "N",
                                "information"))
  expect_identical(gl$N, 301)
  expect_identical(gl$n_parameters, nrow(res))
  # A subset keeps the class; the verbs keep their shape.
  part <- res[res$op == ":=", c("term", "fmi")]
  expect_identical(names(generics::tidy(part)), names(td))
  expect_true(all(is.na(generics::tidy(part)$N_effective)))
  expect_identical(names(generics::glance(part)), names(gl))
  expect_true(is.na(generics::glance(part)$N))
})

test_that("negative estimates from the Hessian estimator are noted under the table", {
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa(.fmi_m3, data = .fmi_hs_mar(), std.lv = TRUE,
                     missing = "ml")
  expect_output(print(fraction_missing_information(fit,
                                                   information = "hessian")),
                "Some estimates are negative")
})

test_that("mlmr() and mlmr_mv() fits are accepted and read in the user's names", {
  skip_if_not_installed("lavaan")
  fit <- mlmr(t26_flags ~ t1_visual_perception + t4_lozenges,
              data = holzinger_swineford, ci_method = "wald")
  res <- fraction_missing_information(fit)
  slopes <- res[res$op == "~", ]
  expect_identical(slopes$term, c("t26_flags ~ t1_visual_perception",
                                  "t26_flags ~ t4_lozenges"))
  expect_equal(slopes$fmi, c(0.5364746748, 0.5979813433), tolerance = 1e-6)
  # The three tests have standard deviations of 7 to 9, and the refits
  # divide them by 8 (see the test of lavaan's FMI above).
  expect_lt(max(abs(res$fmi - .fmi_lavaan(fit$lavaan_fit, "unstructured"))),
            5e-7)

  fit_mv <- mlmr_mv(cbind(t25_paper_form_board_r, t26_flags) ~
                      t1_visual_perception,
                    data = holzinger_swineford, ci_method = "wald",
                    auxiliary = "t2_cubes")
  res_mv <- fraction_missing_information(fit_mv)
  expect_true("t26_flags ~ t1_visual_perception" %in% res_mv$term)
  expect_true("t26_flags ~~ t2_cubes" %in% res_mv$term)
})

test_that("fits outside the definition are refused with the reason", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939
  d <- .fmi_hs_mar()
  m1 <- "visual =~ x1 + x2 + x3"

  expect_error(fraction_missing_information(lavaan::cfa(m1, data = d)),
               "listwise deletion")
  dd <- hs
  dd[c("x1", "x2", "x3")] <- lapply(dd[c("x1", "x2", "x3")], function(x)
    cut(x, 4, ordered_result = TRUE))
  expect_error(
    fraction_missing_information(lavaan::cfa(m1, data = dd, ordered = TRUE)),
    "ordered")
  expect_error(
    fraction_missing_information(lavaan::cfa(m1, data = hs,
                                             estimator = "GLS")),
    "estimator GLS")
  expect_error(
    fraction_missing_information(lavaan::cfa(m1, data = d,
                                             missing = "two.stage")),
    "two-stage")
  expect_error(
    fraction_missing_information(lavaan::cfa(m1, data = d,
                                             missing = "pairwise")),
    "pairwise")

  dl <- d
  dl$clus <- rep(1:30, length.out = nrow(dl))
  fit_ml2 <- suppressWarnings(lavaan::cfa(
    "level: 1
       f =~ x1 + x2 + x3
     level: 2
       x1 ~~ x2 + x3
       x2 ~~ x3",
    data = dl, cluster = "clus", missing = "ml"))
  expect_error(fraction_missing_information(fit_ml2), "multilevel")

  dw <- d
  dw$w <- stats::runif(nrow(dw), .5, 1.5)
  fit_w <- lavaan::cfa(m1, data = dw, missing = "ml",
                       sampling.weights = "w")
  expect_error(fraction_missing_information(fit_w), "sampling weights")

  fit_nc <- suppressWarnings(lavaan::cfa(.fmi_m3, data = d, missing = "ml",
                                         control = list(iter.max = 2)))
  expect_false(lavaan::lavInspect(fit_nc, "converged"))
  expect_error(fraction_missing_information(fit_nc), "did not converge")

  dr <- d
  dr$age <- dr$ageyr + dr$agemo / 12
  fit_cx <- suppressWarnings(lavaan::sem("x5 ~ x4 + age", data = dr,
                                         missing = "ml", fixed.x = TRUE,
                                         conditional.x = TRUE))
  expect_error(fraction_missing_information(fit_cx), "conditional.x")

  fit_efa <- lavaan::efa(data = d[paste0("x", 1:9)], nfactors = 3,
                         missing = "ml", rotation = "oblimin",
                         output = "lavaan")
  expect_error(fraction_missing_information(fit_efa), "rotated exploratory")
  # The guard behind the refusal: a refit at the estimates of a rotated
  # solution rotates it again, and the computation stops.
  expect_error(.fmi_core(fit_efa, "unstructured"), "moved them")

  expect_error(fraction_missing_information(
    fa_confirmatory(holzinger_swineford, list(spatial = .fmi_spatial),
          missing = "fiml")),
    "output = \"fit\"")
  expect_error(fraction_missing_information(lm(x1 ~ x2, data = hs)),
               "must be a fitted lavaan model")
  fit_list <- lavaan::lavaanList(.fmi_m3,
                                 dataList = list(d[1:150, ], d[151:301, ]),
                                 cmd = "cfa", missing = "ml")
  expect_error(fraction_missing_information(fit_list), "lavaanList")
})

test_that("an improper solution lavaan will not refit at its estimates stops the computation and is named", {
  skip_if_not_installed("lavaan")
  # Three items whose correlations (.7, .7, .2) make the first loading
  # larger than the item's standard deviation: the solution's residual
  # variance for y1 is about -1.1 against an observed variance near 1,
  # which lavaan 0.7-3 reads as a solution that has run away. A refit with
  # the optimizer off then restarts from lavaan's simple starting values
  # (loadings .7, unit variances, zero intercepts) instead of holding the
  # estimates, and the guard stops the computation.
  R <- matrix(c(1, .7, .7, .7, 1, .2, .7, .2, 1), 3)
  set.seed(113)
  d <- as.data.frame(matrix(stats::rnorm(600), 200, 3) %*% chol(R))
  names(d) <- c("y1", "y2", "y3")
  d$y2[1:40] <- NA
  d$y3[41:80] <- NA
  # lavaan warns that the solution has run away and that a variance is
  # negative.
  fit <- suppressWarnings(lavaan::cfa("f =~ y1 + y2 + y3", data = d,
                                      missing = "ml", std.lv = TRUE))
  expect_true(lavaan::lavInspect(fit, "converged"))
  theta_y1 <- lavaan::parameterEstimates(fit)
  theta_y1 <- theta_y1$est[theta_y1$lhs == "y1" & theta_y1$op == "~~"]
  expect_lt(theta_y1, -1)
  err <- tryCatch(fraction_missing_information(fit), error = identity)
  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  expect_match(msg, "moved them", fixed = TRUE)
  expect_match(msg, "a rotated exploratory one", fixed = TRUE)
  expect_match(msg, "an improper solution, such as a Heywood case",
               fixed = TRUE)
  expect_match(msg, "This fit's solution is improper", fixed = TRUE)
  # lavaan's own FMI makes the same refit and is silently wrong here: the
  # intercept of y1, which is never missing, gets a negative FMI.
  pe <- suppressWarnings(lavaan::parameterEstimates(fit, fmi = TRUE))
  expect_lt(pe$fmi[pe$lhs == "y1" & pe$op == "~1"], -0.1)
  # A Heywood case that is not that far below zero (-0.07 here) is refit at
  # its estimates, and its FMI is computed.
  R2 <- matrix(c(1, .6, .6, .6, 1, .3, .6, .3, 1), 3)
  set.seed(113)
  d2 <- as.data.frame(matrix(stats::rnorm(600), 200, 3) %*% chol(R2))
  names(d2) <- c("y1", "y2", "y3")
  d2$y2[1:40] <- NA
  d2$y3[41:80] <- NA
  fit2 <- suppressWarnings(lavaan::cfa("f =~ y1 + y2 + y3", data = d2,
                                       missing = "ml", std.lv = TRUE))
  pe2 <- lavaan::parameterEstimates(fit2)
  expect_lt(pe2$est[pe2$lhs == "y1" & pe2$op == "~~"], 0)
  tab <- fraction_missing_information(fit2)
  expect_true(all(is.finite(tab$fmi)))
})

test_that("a solution with factors correlated above 1 is refit at its estimates", {
  skip_if_not_installed("lavaan")
  # Two factors of three indicators each, drawn from one factor, n = 60:
  # the fit's factor correlation is 1.077 and no variance is negative.
  # lavaan's check of starting values would shrink the covariance in the
  # refit to a correlation of .91, and its own parameterEstimates(fmi =
  # TRUE), which runs that check, reports an FMI of -1.62 for it. Every
  # refit holds the estimates, so the FMIs are those of this solution.
  set.seed(21)
  n <- 60
  f <- stats::rnorm(n)
  mk <- function(l) l * f + stats::rnorm(n, sd = sqrt(1 - l^2))
  d <- data.frame(a1 = mk(0.7), a2 = mk(0.6), a3 = mk(0.65),
                  b1 = mk(0.7), b2 = mk(0.6), b3 = mk(0.62))
  d$b3[stats::runif(n) < stats::plogis(-1 + d$a1)] <- NA
  d$a2[stats::runif(n) < 0.2] <- NA
  set.seed(1)
  d[] <- lapply(d, function(v) v + stats::rnorm(n, sd = 0.3))
  # lavaan warns that the covariance matrix of the latent variables is not
  # positive definite.
  fit <- suppressWarnings(lavaan::cfa("A =~ a1 + a2 + a3\n B =~ b1 + b2 + b3",
                                      data = d, missing = "ml",
                                      std.lv = TRUE))
  expect_gt(lavaan::lavInspect(fit, "cor.lv")[1, 2], 1)
  expect_no_warning(res <- fraction_missing_information(fit))
  ab <- res[res$term == "A ~~ B", ]
  expect_equal(ab$estimate, 1.077003, tolerance = 1e-6)
  expect_equal(ab$fmi, 0.2492837, tolerance = 1e-6)
  pe <- suppressWarnings(lavaan::parameterEstimates(fit, fmi = TRUE))
  expect_lt(pe$fmi[pe$lhs == "A" & pe$op == "~~" & pe$rhs == "B"], -1)
})

# The reviewer's regression: y on x1, x2, and x3, with x3 missing at
# random given x1 for about 65 percent of the cases.
.fmi_units_data <- function() {
  set.seed(3)
  n <- 200
  R <- matrix(0.5, 4, 4)
  diag(R) <- 1
  x <- as.data.frame(matrix(stats::rnorm(n * 4), n) %*% chol(R))
  names(x) <- c("y", "x1", "x2", "x3")
  x$x3[stats::runif(n) < stats::plogis(1 + 2 * x$x1)] <- NA
  x$x2[sample(n, 3)] <- NA
  x
}

.fmi_units_fit <- function(d, ...) {
  suppressWarnings(suppressMessages(
    lavaan::sem("y ~ x1 + x2 + x3", data = d, missing = "ml",
                fixed.x = FALSE, ...)))
}

test_that("the FMI does not depend on the units of the variables", {
  skip_if_not_installed("lavaan")
  # The FMI was first computed from refits in the variables' own metric:
  # with x3 multiplied by .001 the FMIs moved by up to .024, with every
  # variable so multiplied the slope FMIs were .092, .292, and .740 in
  # place of .059, .188, and .700, and the mlmr() degrees of freedom of the
  # slope of x2 fell from 158.5 to 137.7. Then the refits took the source
  # fit's moments, which left the Hessian estimator off by .5 at a
  # standard deviation of 100,000 and every FMI NA, and every mlmr()
  # degree of freedom with it, once the outcome's standard deviation
  # reached about 7e7, where lavaan could not invert the refits'
  # information. The refits are now made with the standard deviations
  # near 1.
  x <- .fmi_units_data()
  infs <- c(unstructured = "unstructured", structured = "structured",
            hessian = "hessian")
  base <- lapply(infs, function(inf) {
    fraction_missing_information(.fmi_units_fit(x), inf)
  })
  expect_equal(base$unstructured$fmi[1:3], c(0.0592, 0.1877, 0.6999),
               tolerance = 1e-3)
  for (k in c(1e-3, 1e4, 1e5)) {
    for (v in list("x3", "y", names(x))) {
      d <- x
      d[v] <- d[v] * k
      fit <- .fmi_units_fit(d)
      for (inf in names(infs)) {
        expect_no_warning(res <- fraction_missing_information(fit, inf))
        expect_lt(max(abs(res$fmi - base[[inf]]$fmi)), 1e-5)
      }
    }
  }
  # An outcome in dollars, with a standard deviation of 100,000,000.
  d <- x
  d$y <- d$y * 1e8
  fit <- .fmi_units_fit(d)
  for (inf in names(infs)) {
    res <- fraction_missing_information(fit, inf)
    expect_true(all(is.finite(res$fmi)))
    expect_lt(max(abs(res$fmi - base[[inf]]$fmi)), 1e-5)
  }

  # The missing information degrees of freedom of mlmr() and mlmr_mv().
  df_of <- function(f, d) {
    attr(suppressWarnings(suppressMessages(f(d))), "wald_df")
  }
  one <- function(d) {
    mlmr(y ~ x1 + x2 + x3, data = d, ci_method = "wald",
         wald_df = "missing_information", effect_sizes = FALSE)
  }
  two <- function(d) {
    mlmr_mv(cbind(y, x3) ~ x1 + x2, data = d, ci_method = "wald",
            wald_df = "missing_information", effect_sizes = FALSE)
  }
  df_one <- df_of(one, x)
  df_two <- df_of(two, x)
  expect_equal(unname(df_one), c(124.3, 184.2, 158.5, 56.0),
               tolerance = 1e-3)
  for (k in c(1e-3, 1e5, 1e8)) {
    expect_equal(df_of(one, x * k), df_one, tolerance = 1e-6, info = k)
    expect_equal(df_of(two, x * k), df_two, tolerance = 1e-6, info = k)
  }
  d <- x
  d$x3 <- d$x3 * 1e-3
  expect_equal(df_of(one, d), df_one, tolerance = 1e-6)
  d <- x
  d$y <- d$y * 1e8
  m <- suppressWarnings(suppressMessages(one(d)))
  expect_equal(attr(m, "wald_df"), df_one, tolerance = 1e-6)
  expect_true(all(is.finite(m$coef_table$p_value)))
  expect_true(all(is.finite(m$coef_table$ci_lower)))
})

test_that("the FMI does not depend on the units where lavaan does not rescale", {
  skip_if_not_installed("lavaan")
  # Two groups of the reviewer's regression, refit with every variable
  # multiplied by .01 and by 2^-8 at the estimates of the fit to the
  # original data, so that only the computation of the FMI is compared
  # (lavaan's own fit is less precise at such scales, which it does not
  # rescale). The refits estimated the unrestricted moments in the raw
  # metric, where lavaan's EM, which stops on an absolute change of 1e-5,
  # stops short for variances near 1e-4: group b's slope FMIs were off by
  # .0034 at .01 and by .022 at 2^-8 (.0039).
  set.seed(3)
  n <- 200
  R <- matrix(0.5, 4, 4)
  diag(R) <- 1
  x <- do.call(rbind, lapply(c("a", "b"), function(g) {
    d <- as.data.frame(matrix(stats::rnorm(n * 4), n) %*% chol(R))
    names(d) <- c("y", "x1", "x2", "x3")
    d$x3[stats::runif(n) < stats::plogis(1 + 2 * d$x1)] <- NA
    d$x2[sample(n, 3)] <- NA
    d$g <- g
    d
  }))
  v <- c("y", "x1", "x2", "x3")
  fit1 <- .fmi_units_fit(x, group = "g")
  base <- fraction_missing_information(fit1)
  at_scale <- function(k) {
    d <- x
    d[v] <- d[v] * k
    start <- lavaan::parTable(fit1)
    start$est <- start$est * ifelse(start$op == "~~", k^2,
                                     ifelse(start$op == "~1", k, 1))
    .fmi_units_fit(d, group = "g", start = start, optim.method = "none")
  }
  res <- fraction_missing_information(at_scale(0.01))
  expect_lt(max(abs(res$fmi - base$fmi)), 1e-5)
  # A power of two leaves the refits' data exactly as they were, so the
  # FMIs are the same to rounding error.
  res <- fraction_missing_information(at_scale(2^-8))
  expect_equal(res$fmi, base$fmi, tolerance = 1e-12)
  expect_equal(res$se_observed / res$se_complete,
               base$se_observed / base$se_complete, tolerance = 1e-12)
})

test_that("a fit that orders its variables as the data do gets the same FMI", {
  skip_if_not_installed("lavaan")
  # lavaan's ov_order = "data" orders the variables as the data columns,
  # here the reverse of the model's order. The refit was handed the source
  # fit's unrestricted moments by position and ordered its variables as
  # the model does, so it read a permuted covariance matrix: FMIs of -1.58
  # for x2 ~~ x2 and .88 for the intercept of x1, which is never missing,
  # under the unstructured information, and up to 43 off under the
  # structured one.
  hs <- lavaan::HolzingerSwineford1939
  set.seed(11)
  hs$x5[stats::runif(301) < stats::plogis(-1 + hs$x4)] <- NA
  hs$x9[stats::runif(301) < 0.25] <- NA
  d <- hs[c(paste0("x", 9:1), "school")]
  key <- function(r) paste(r$term, r$group)
  for (group in list(NULL, "school")) {
    f_model <- suppressWarnings(lavaan::cfa(.fmi_m3, data = d, missing = "ml",
                                            group = group))
    f_data <- suppressWarnings(lavaan::cfa(.fmi_m3, data = d, missing = "ml",
                                           group = group, ov.order = "data"))
    expect_identical(lavaan::lavNames(f_data, "ov")[1:2], c("x9", "x8"))
    for (inf in c("unstructured", "structured", "hessian")) {
      a <- fraction_missing_information(f_model, inf)
      b <- fraction_missing_information(f_data, inf)
      b <- b[match(key(a), key(b)), ]
      expect_equal(b$fmi, a$fmi, tolerance = 1e-8,
                   info = paste(inf, length(group)))
      expect_equal(b$se_observed, a$se_observed, tolerance = 1e-8)
    }
  }
  a <- fraction_missing_information(f_model)
  expect_lt(max(abs(a$fmi[a$op == "~1" & a$term == "x1 ~1"])), 1e-10)
})

test_that("the refits carry the fit's estimates and constraints across units exactly", {
  skip_if_not_installed("lavaan")
  # With no missing values the unrestricted moments are the sample
  # moments, so the observed standard errors of the refits, made with the
  # variables divided by powers of two, must be lavaan's own for the fit
  # in its units: every estimate, fixed value, defined parameter, and
  # constraint goes through the map. The model ties a loading of a test in
  # one unit to that of a test in another (a constraint lavaan itself does
  # not rescale), defines a parameter through both, and adds an equality
  # constraint on intercepts across the schools.
  hs <- lavaan::HolzingerSwineford1939
  d <- hs[c(paste0("x", 1:9), "school")]
  scale <- c(x1 = 20, x3 = 0.05, x4 = 8, x6 = 0.1, x7 = 30, x9 = 3)
  for (v in names(scale)) d[[v]] <- d[[v]] * scale[[v]]
  model <- "visual =~ x1 + c(a, a)*x2 + x3
            textual =~ x5 + c(b, b)*x4 + x6
            speed =~ x7 + x8 + x9
            ab := a*b + a
            a*8 == b"
  # lavaan notes that the variances differ by a factor of 1000. Made in
  # the variables' own metric (rescale.data = FALSE), lavaan's fit stops
  # short of the maximum (a Newton step moves the loading of x8 by one
  # percent), which the refits find; with no value missing every FMI is 0
  # there and at the maximum alike, so the warning, which says that the
  # FMIs are not those of the maximum likelihood estimates, is not given.
  # The standard errors are evaluated at the fit's estimates, as lavaan's
  # are.
  fit <- suppressWarnings(lavaan::cfa(
    model, data = d, missing = "ml", group = "school",
    group.equal = "intercepts", information = "observed",
    observed.information = "h1", h1.information = "unstructured",
    rescale.data = FALSE))
  expect_true(lavaan::lavInspect(fit, "converged"))
  expect_no_warning(res <- fraction_missing_information(fit))
  pt <- lavaan::parTable(fit)
  keep <- pt$free > 0L | pt$op == ":="
  expect_equal(res$se_observed, pt$se[keep], tolerance = 1e-8)
  expect_lt(max(abs(res$fmi)), 1e-10)
  expect_true("ab" %in% res$term)
  # The refits were made in another metric: their standard deviations
  # were near 1.
  sd_range <- attr(suppressWarnings(.fmi_core(fit, "unstructured")),
                   "sd_range")
  expect_gt(sd_range[1L], 0.5)
  expect_lt(sd_range[2L], 2)
})

test_that("the Hessian estimator warns when the refits keep a badly scaled metric", {
  skip_if_not_installed("lavaan")
  # A model with a product term cannot be carried to another metric (nor
  # does lavaan rescale it), so the refits keep the variables' own, and
  # the Hessian estimator's numeric step does not suit an outcome with a
  # standard deviation near 10,000.
  d <- .fmi_units_data()
  d$y <- d$y * 1e4
  fit <- suppressWarnings(lavaan::sem("y ~ x1 + x2 + x1:x2", data = d,
                                      missing = "ml.x", fixed.x = FALSE))
  expect_identical(attr(.fmi_core(fit, "hessian"), "scale_reason"),
                   "product terms")
  expect_warning(fraction_missing_information(fit, information = "hessian"),
                 paste0("standard deviations run from [0-9.e+]+ to [0-9.e+]+, ",
                        "and the refits behind the FMI could not be moved ",
                        "to a metric in which they are near 1 because the ",
                        "model has product terms\\. The Hessian"))
  expect_no_warning(fraction_missing_information(fit))
  expect_no_warning(fraction_missing_information(.fmi_spatial_fit(),
                                                 information = "hessian"))
})

test_that("a refit whose information cannot be inverted stops with the reason", {
  skip_if_not_installed("lavaan")
  # The product term model with the outcome's standard deviation near
  # 1e8: in its own metric neither lavaan's fit nor the refits can invert
  # the information. The FMIs were NA, with lavaan's warning that the
  # model may not be identified.
  d <- .fmi_units_data()
  d$y <- d$y * 1e8
  fit <- suppressWarnings(lavaan::sem("y ~ x1 + x2 + x1:x2", data = d,
                                      missing = "ml.x", fixed.x = FALSE))
  err <- tryCatch(fraction_missing_information(fit), error = identity)
  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  expect_match(msg, "information matrix of a refit at the fit's estimates ",
               fixed = TRUE)
  expect_match(msg, "because the model has product terms, so the failure is ",
               fixed = TRUE)
  expect_match(msg, "which leaves the fraction unchanged", fixed = TRUE)
  expect_no_match(msg, "not identified", fixed = TRUE)
  # A factor with two indicators is not identified. lavaan returns
  # standard errors for it all the same (580 for the free loading), and
  # the refits cannot invert the information.
  hs <- .fmi_hs_mar()
  fit_ni <- suppressWarnings(lavaan::cfa("f =~ x1 + x2", data = hs,
                                         missing = "ml"))
  expect_error(fraction_missing_information(fit_ni),
               "singular only for a model that is not identified",
               fixed = TRUE)
  # The other wordings.
  pt <- lavaan::parTable(lavaan::cfa(.fmi_m3, data = hs, missing = "ml"))
  opt <- list(se = "standard")
  expect_error(.fmi_stop_not_inverted(pt, opt, "hessian", c(1, 2), TRUE, NULL),
               "The hessian information can be singular", fixed = TRUE)
  pt$se[pt$free > 0L] <- NA_real_
  expect_error(.fmi_stop_not_inverted(pt, opt, "unstructured", c(1, 2), TRUE,
                                      NULL),
               paste0("standard errors either, which happens when the ",
                      "model is not identified"),
               fixed = TRUE)
})

test_that("a variable far from zero is named when the information cannot be inverted", {
  skip_if_not_installed("lavaan")
  # A predictor with a standard deviation of 1 and a mean of 10,000:
  # dividing by a power of two leaves it as far from zero, and neither
  # lavaan's fit nor the refits can invert the information. The stop said
  # that this happens when the model is not identified.
  d <- .fmi_units_data()
  d$x2 <- d$x2 + 1e4
  fit <- suppressWarnings(lavaan::sem("y ~ x1 + x2 + x3", data = d,
                                      missing = "ml", fixed.x = FALSE))
  err <- tryCatch(fraction_missing_information(fit), error = identity)
  expect_s3_class(err, "dmar_fmi_not_inverted")
  msg <- conditionMessage(err)
  expect_match(msg, "x2 lies 10010 standard deviations from zero",
               fixed = TRUE)
  expect_match(msg, "the coefficients of that variable are nearly collinear",
               fixed = TRUE)
  expect_match(msg, "(center it)", fixed = TRUE)
  expect_no_match(msg, "not identified", fixed = TRUE)
  # Centered, the slopes' FMIs are those of the data as generated.
  d$x2 <- d$x2 - 1e4
  centered <- fraction_missing_information(.fmi_units_fit(d))
  base <- fraction_missing_information(.fmi_units_fit(.fmi_units_data()))
  expect_equal(centered$fmi[1:3], base$fmi[1:3], tolerance = 1e-8)
  # Below zero the distance is the same; the stop said the mean was 10010
  # times the standard deviation, of a mean of -10,000.
  d$x2 <- d$x2 - 1e4
  msg <- tryCatch(fraction_missing_information(.fmi_units_fit(d)),
                  error = conditionMessage)
  expect_match(msg, "x2 lies 10010 standard deviations from zero",
               fixed = TRUE)
  # Two such variables are named together, in the plural; the stop said
  # "The mean of x1, x2 is between ... times the variable's standard
  # deviation", "the coefficients of the variable", and "(center it)".
  d <- .fmi_units_data()
  d$x2 <- d$x2 + 1e4
  d$x1 <- d$x1 + 2e4
  msg <- tryCatch(fraction_missing_information(.fmi_units_fit(d)),
                  error = conditionMessage)
  expect_match(msg, paste0("x1 and x2 lie between 10010 and 20184 standard ",
                           "deviations from zero"), fixed = TRUE)
  expect_match(msg, "the coefficients of those variables are nearly",
               fixed = TRUE)
  expect_match(msg, "from each of them (center them)", fixed = TRUE)
  expect_no_match(msg, "the variable's", fixed = TRUE)
  expect_no_match(msg, "(center it)", fixed = TRUE)
  # mlmr() fits its model to renamed variables (y, x_1, x_2, ...); its
  # missing information degrees of freedom stop with the variable named as
  # the formula names it.
  d <- stats::setNames(.fmi_units_data(), c("score", "age", "income", "hours"))
  d$income <- d$income + 1e4
  err <- tryCatch(
    suppressWarnings(mlmr(score ~ age + income + hours, data = d,
                          ci_method = "wald", effect_sizes = FALSE,
                          wald_df = "missing_information")),
    error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "income lies 10010 standard deviations",
               fixed = TRUE)
  expect_no_match(conditionMessage(err), "x_[0-9]")
})

test_that("a variable far from zero is not named when the model is not identified", {
  skip_if_not_installed("lavaan")
  # The 1.1.0 review: the stop named any variable 1000 or more standard
  # deviations from zero, before it considered identification, and said
  # to center it. A model that is not identified is not identified
  # whatever its variables' means, and an identified regression gave the
  # unshifted FMIs with a predictor 7000 standard deviations from zero.
  d <- .fmi_units_data()
  d$x2 <- d$x2 + 2000
  ni <- "y ~ x1 + x2 + x3
         y ~~ x1"
  fit <- suppressWarnings(lavaan::sem(ni, data = d, missing = "ml",
                                      fixed.x = FALSE))
  msg <- tryCatch(fraction_missing_information(fit), error = conditionMessage)
  expect_match(msg, "singular only for a model that is not identified",
               fixed = TRUE)
  expect_no_match(msg, "center", fixed = TRUE)
  expect_no_match(msg, "standard deviations from zero", fixed = TRUE)
  # The same regression identified computes the unshifted FMIs.
  base <- fraction_missing_information(.fmi_units_fit(.fmi_units_data()))
  shifted <- fraction_missing_information(.fmi_units_fit(d))
  expect_equal(shifted$fmi[1:3], base$fmi[1:3], tolerance = 1e-6)
  # An indicator's location enters only its own intercept and is never
  # the cause: a factor with two indicators, one of them 2000 standard
  # deviations from zero.
  hs <- .fmi_hs_mar()
  hs$x1 <- hs$x1 + 2000 * stats::sd(hs$x1)
  fit <- suppressWarnings(lavaan::cfa("f =~ x1 + x2", data = hs,
                                      missing = "ml"))
  msg <- tryCatch(fraction_missing_information(fit), error = conditionMessage)
  expect_match(msg, "singular only for a model that is not identified",
               fixed = TRUE)
  expect_no_match(msg, "center", fixed = TRUE)
  # Where centering would change an intercept the model constrains (two
  # outcomes of x2 whose intercepts are held equal), the centered model is
  # another model, the check is not made, and the stop names both causes.
  u <- .fmi_units_data()
  fit_eq <- suppressWarnings(lavaan::sem("y ~ x1 + x2 + x3
                                          x3 ~ x2
                                          y ~ a*1
                                          x3 ~ a*1", data = u,
                                         missing = "ml", fixed.x = FALSE))
  pt <- lavaan::parTable(fit_eq)
  opt <- lavaan::lavInspect(fit_eq, "options")
  cause <- .fmi_location_confirms(fit_eq, pt, opt,
                                  list(lavaan::lavInspect(fit_eq, "data")),
                                  "x2", "unstructured")
  expect_true(is.na(cause$confirmed))
  expect_match(cause$why, "fixes or constrains", fixed = TRUE)
  msg <- tryCatch(.fmi_stop_not_inverted(pt, opt, "unstructured", c(1, 2),
                                         TRUE, NULL, location = c(x2 = 1e4),
                                         location_cause = cause),
                  error = conditionMessage)
  expect_match(msg, paste0("x2 lies 10000 standard deviations from zero, ",
                           "and so far from zero an intercept and the ",
                           "coefficients of that variable can be nearly ",
                           "collinear"), fixed = TRUE)
  expect_match(msg, "could not be checked, because centering would",
               fixed = TRUE)
  expect_match(msg, "centering cannot repair a model that is not identified",
               fixed = TRUE)
})

test_that("a variable far from zero is named in a multiple group fit", {
  skip_if_not_installed("lavaan")
  # The 1.1.0 review: lavaan names the groups' moments by group, so the
  # locations came out as "a.x2" and "b.x2", which no regressor's name
  # matched; the stop never named x2 and said only that the model may not
  # be identified, while before the check of the cause it had named
  # "a.x2, b.x2".
  d <- .fmi_units_data()
  d$g <- rep(c("a", "b"), length.out = nrow(d))
  d$x2 <- d$x2 + 1e4
  fit <- .fmi_units_fit(d, group = "g")
  err <- tryCatch(fraction_missing_information(fit), error = identity)
  expect_s3_class(err, "dmar_fmi_not_inverted")
  msg <- conditionMessage(err)
  expect_match(msg, "x2 lies 10179 standard deviations from zero",
               fixed = TRUE)
  expect_match(msg, "the same model fit with it centered has an information ",
               fixed = TRUE)
  expect_match(msg, "(center it)", fixed = TRUE)
  expect_no_match(msg, "a.x2", fixed = TRUE)
  expect_no_match(msg, "not identified", fixed = TRUE)
  # Centered, as the stop says, the fraction is computed.
  d$x2 <- d$x2 - mean(d$x2, na.rm = TRUE)
  centered <- fraction_missing_information(.fmi_units_fit(d, group = "g"))
  expect_true(all(is.finite(centered$fmi[centered$term == "y ~ x2"])))
  # A latent outcome of x2 in two groups: centering would change its
  # intercept, which the model fixes at 0, so the check is not made and
  # both causes are named.
  set.seed(5)
  n <- 400
  x2 <- stats::rnorm(n)
  x1 <- stats::rnorm(n)
  f <- 0.5 * x1 + 0.4 * x2 + stats::rnorm(n)
  lat <- data.frame(y1 = f + stats::rnorm(n), y2 = 0.8 * f + stats::rnorm(n),
                    y3 = 0.7 * f + stats::rnorm(n), x1 = x1, x2 = x2 + 5000,
                    g = rep(c("a", "b"), length.out = n))
  lat$y3[stats::runif(n) < stats::plogis(-1 + x1)] <- NA
  fit_lat <- suppressWarnings(lavaan::sem("f =~ y1 + y2 + y3
                                           f ~ x1 + x2", data = lat,
                                          missing = "ml", fixed.x = FALSE,
                                          group = "g"))
  msg <- tryCatch(fraction_missing_information(fit_lat),
                  error = conditionMessage)
  expect_match(msg, paste0("x2 lies [0-9]+ standard deviations from zero, and ",
                           "so far from zero an intercept and the ",
                           "coefficients of that variable can be nearly ",
                           "collinear"))
  expect_match(msg, "could not be checked, because centering would",
               fixed = TRUE)
  expect_match(msg, "centering cannot repair a model that is not identified",
               fixed = TRUE)
  # The outcome's intercept held equal across the groups is such a
  # constraint as well.
  u <- .fmi_units_data()
  u$g <- rep(c("a", "b"), length.out = nrow(u))
  fit_eq <- .fmi_units_fit(u, group = "g", group.equal = "intercepts")
  cause <- .fmi_location_confirms(fit_eq, lavaan::parTable(fit_eq),
                                  lavaan::lavInspect(fit_eq, "options"),
                                  lavaan::lavInspect(fit_eq, "data"), "x2",
                                  "unstructured",
                                  lavaan::lavInspect(fit_eq, "group.label"))
  expect_true(is.na(cause$confirmed))
  expect_match(cause$why, "fixes or constrains", fixed = TRUE)
})

test_that("a fit that is not at the maximum of the likelihood is named", {
  skip_if_not_installed("lavaan")
  # lavaan reports these fits as converged, and its check of the gradient,
  # made in the variables' own metric, passes them; the FMIs were computed
  # at their estimates without a word. The factor model of the nine tests,
  # identified by marker indicators, with every variable multiplied by 100
  # (standard deviations from 100 to 160, which lavaan does not rescale):
  # its log likelihood is .11 below the maximum and its FMIs were off by
  # up to .011.
  hs <- lavaan::HolzingerSwineford1939
  set.seed(11)
  hs$x5[stats::runif(301) < stats::plogis(-1 + hs$x4)] <- NA
  hs$x9[stats::runif(301) < 0.25] <- NA
  v <- paste0("x", 1:9)
  base <- fraction_missing_information(lavaan::cfa(.fmi_m3, data = hs[v],
                                                   missing = "ml"))
  d <- hs[v] * 100
  fit <- lavaan::cfa(.fmi_m3, data = d, missing = "ml")
  expect_true(lavaan::lavInspect(fit, "converged"))
  expect_warning(res <- fraction_missing_information(fit),
                 "not at the maximum of the likelihood", fixed = TRUE)
  expect_gt(max(abs(res$fmi - base$fmi)), 0.005)
  # Fit with lavaan's rescaling, it is at the maximum, and its FMIs are
  # those of the variables in their own units.
  fit_r <- suppressMessages(lavaan::cfa(.fmi_m3, data = d, missing = "ml",
                                        rescale_data = TRUE))
  expect_no_warning(res_r <- fraction_missing_information(fit_r))
  expect_lt(max(abs(res_r$fmi - base$fmi)), 1e-5)
  # A regression whose columns are not in the model's order, made with
  # ov_order = "data" and rescaled by lavaan: in lavaan 0.7-3 it is not the
  # maximum likelihood fit (its log likelihood is 50 below that of the
  # same fit in the model's order), and its slopes' FMIs were .278, .221,
  # and .699 in place of .059, .188, and .700.
  x <- .fmi_units_data()
  x$x3 <- x$x3 * 0.004
  x$y <- x$y * 0.03
  x$x1 <- x$x1 * 0.5
  x <- x[c("x3", "x2", "x1", "y")]
  fit_data <- .fmi_units_fit(x, rescale.data = TRUE, ov.order = "data")
  fit_model <- .fmi_units_fit(x, rescale.data = TRUE)
  expect_true(lavaan::lavInspect(fit_data, "converged"))
  expect_warning(fraction_missing_information(fit_data),
                 "ov_order = \"data\"", fixed = TRUE)
  expect_no_warning(fraction_missing_information(fit_model))
  # A regression with every standard deviation near .003 stops short by a
  # step of .14 percent, with FMIs off by .001; near .01 the step is 6e-5
  # of an estimate, and the FMIs are off by 6e-5, which passes.
  u <- .fmi_units_data()
  expect_warning(fraction_missing_information(.fmi_units_fit(u * 0.003)),
                 "not at the maximum of the likelihood", fixed = TRUE)
  expect_no_warning(fraction_missing_information(.fmi_units_fit(u * 0.01)))
  # mlmr() takes the missing information degrees of freedom from the same
  # computation and gives the same warning, among lavaan's own.
  said <- new.env()
  said$w <- character(0)
  withCallingHandlers(
    mlmr(y ~ x1 + x2 + x3, data = u * 0.003, ci_method = "wald",
         wald_df = "missing_information", effect_sizes = FALSE),
    warning = function(w) {
      said$w <- c(said$w, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("not at the maximum of the likelihood", said$w,
                        fixed = TRUE)))
  # It names the estimate as the formula does, not as the fit does.
  expect_false(any(grepl("x_[0-9]", said$w)))
  said$w <- character(0)
  withCallingHandlers(
    fraction_missing_information(
      mlmr(y ~ x1 + x2 + x3, data = u * 0.003, ci_method = "wald",
           effect_sizes = FALSE)),
    warning = function(w) {
      said$w <- c(said$w, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("not at the maximum of the likelihood", said$w,
                        fixed = TRUE)))
  expect_false(any(grepl("x_[0-9]", said$w)))
  # A fit at a maximum with equality constraints across groups is not
  # named, nor is one whose inequality constraint is active, where the
  # gradient is not zero at the constrained maximum.
  fit_eq <- lavaan::cfa(.fmi_m3, data = .fmi_hs_mar(), std.lv = TRUE,
                        missing = "ml", group = "school",
                        group.equal = c("loadings", "intercepts"))
  expect_no_warning(fraction_missing_information(fit_eq))
  fit_in <- suppressWarnings(lavaan::sem("y ~ b1*x1 + x2 + x3
                                          b1 < -0.1", data = u,
                                         missing = "ml", fixed.x = FALSE))
  expect_no_warning(fraction_missing_information(fit_in))
})

test_that("a fit short of the maximum on a flat likelihood is named only when its FMIs move", {
  skip_if_not_installed("lavaan")
  # The 1.1.0 review: in a small sample with most of the outcome missing,
  # lavaan's optimizer stops about 5e-4 standard errors short of the
  # maximum (a log likelihood 5e-7 below it), and the Newton step, scaled
  # by the estimate rather than by how flat the likelihood is, passed its
  # threshold. The warning blamed units and ov_order, neither of which
  # applied, and prescribed rescale_data = TRUE, which did not remove it.
  # The FMIs one step on are now compared with those at the estimates.
  flat <- function(seed) {
    set.seed(seed)
    n <- 60
    R <- matrix(0.4, 5, 5)
    diag(R) <- 1
    x <- as.data.frame(matrix(stats::rnorm(n * 5), n) %*% chol(R))
    names(x) <- c("y", paste0("x", 1:4))
    x$y[stats::runif(n) < stats::plogis(0.8 + 1.5 * x$x1)] <- NA
    x$x4[stats::runif(n) < 0.5] <- NA
    suppressWarnings(lavaan::sem("y ~ x1 + x2 + x3 + x4", data = x,
                                 missing = "ml", fixed.x = FALSE))
  }
  # An interior EM point on this draw is not a global H1 maximum: its
  # three complete rows permit collapse in five dimensions.
  expect_error(fraction_missing_information(flat(7)),
    class = "dmar_fmi_saturated_singular")
  # FMIs 6e-4 off (y ~~ y, .9325 against .9332): the warning, with the
  # change and the cause that applies, and no remedy that does not help.
  expect_warning(fraction_missing_information(flat(3)),
                 paste0("changes the fraction of missing information of ",
                        "y ~~ y from 0.9325 to 0.9332.*flat in some ",
                        "direction"))
  w <- tryCatch(fraction_missing_information(flat(3)),
                warning = conditionMessage)
  expect_no_match(w, "Fit the model again with rescale_data", fixed = TRUE)
  expect_no_match(w, "small or large units", fixed = TRUE)
})

test_that("the maximum warning names units and their remedy only for a fit lavaan did not rescale", {
  skip_if_not_installed("lavaan")
  # The 1.1.0 review: the same small sample with every variable multiplied
  # by 20 (standard deviations of 11 to 22), fit with rescale_data = TRUE,
  # stops short on its flat likelihood as it does in its own units, and
  # the warning blamed variables "measured in small or large units that it
  # does not rescale" and said to fit the model again with rescale_data =
  # TRUE, the remedy already taken; fit without it, the warning is not
  # given at all.
  set.seed(3)
  n <- 60
  R <- matrix(0.4, 5, 5)
  diag(R) <- 1
  x <- as.data.frame(matrix(stats::rnorm(n * 5), n) %*% chol(R))
  names(x) <- c("y", paste0("x", 1:4))
  x$y[stats::runif(n) < stats::plogis(0.8 + 1.5 * x$x1)] <- NA
  x$x4[stats::runif(n) < 0.5] <- NA
  m <- "y ~ x1 + x2 + x3 + x4"
  fit_r <- suppressMessages(lavaan::sem(m, data = x * 20, missing = "ml",
                                        fixed.x = FALSE, rescale_data = TRUE))
  w <- tryCatch(fraction_missing_information(fit_r),
                warning = conditionMessage)
  expect_match(w, "changes the fraction of missing information of y ~~ y",
               fixed = TRUE)
  expect_match(w, "flat in some direction", fixed = TRUE)
  expect_no_match(w, "small or large units", fixed = TRUE)
  expect_no_match(w, "Fit the model again with rescale_data", fixed = TRUE)
  expect_no_warning(fraction_missing_information(
    lavaan::sem(m, data = x * 20, missing = "ml", fixed.x = FALSE)))
  # The rule itself, on one change of an FMI larger than 1e-4: units are
  # named for a fit lavaan did not rescale whose standard deviations reach
  # 100 (not 20), and the order of the variables only for a fit lavaan
  # rescaled.
  said <- function(fit) {
    pt <- lavaan::parTable(fit)
    here <- rep(0.5, length(pt$lhs))
    step <- here
    step[1L] <- 0.51
    h1 <- lavaan::lavInspect(fit, "sampstat.h1")
    rows <- which(pt$free > 0L)
    screen <- list(term = "y ~ x1", est = 1, est_step = 1.01, gain = 0.1,
                   rows = rows, rel = rep(1, length(rows)))
    tryCatch(.fmi_warn_maximum(screen, fit, pt, here, step, h1$cov, NULL),
             warning = function(w) {
               expect_s3_class(w, "dmar_fmi_not_at_maximum")
               expect_identical(w$lhs, pt$lhs[rows[1L]])
               expect_identical(w$from, screen$est)
               expect_identical(w$gain, screen$gain)
               conditionMessage(w)
             })
  }
  u <- .fmi_units_data()
  fit_20 <- .fmi_units_fit(u * 20, rescale.data = FALSE)
  fit_200 <- .fmi_units_fit(u * 200, rescale.data = FALSE)
  fit_200r <- .fmi_units_fit(u * 200, rescale.data = TRUE)
  expect_match(said(fit_20), "flat in some direction", fixed = TRUE)
  expect_match(said(fit_200), "small or large units that it does not rescale",
               fixed = TRUE)
  expect_match(said(fit_200), "Fit the model again with rescale_data = TRUE",
               fixed = TRUE)
  expect_match(said(fit_200r), "flat in some direction", fixed = TRUE)
  expect_no_match(said(fit_200r), "rescale_data = TRUE, and", fixed = TRUE)
  p <- u[c("x3", "x2", "x1", "y")]
  fit_order <- .fmi_units_fit(p, ov.order = "data", rescale.data = FALSE)
  fit_order_r <- .fmi_units_fit(p, ov.order = "data", rescale.data = TRUE)
  expect_no_match(said(fit_order), "ov_order", fixed = TRUE)
  expect_match(said(fit_order_r), "ov_order = \"data\"", fixed = TRUE)
  expect_match(said(fit_order_r),
               "columns of the data in the model's order", fixed = TRUE)
})

# The calibration study's hardest cell as the 1.1.0 review drew it: one
# factor, 12 items with loadings from .8 to .5, N = 50, and 30 percent of
# the values of every second item missing completely at random; draw r
# after set.seed(1).
.fmi_sparse_draw <- function(r) {
  set.seed(1)
  for (i in seq_len(r)) {
    lam <- seq(0.8, 0.5, length.out = 12)
    eta <- stats::rnorm(50)
    d <- as.data.frame(sapply(lam, function(l) {
      l * eta + stats::rnorm(50, sd = sqrt(1 - l^2))
    }))
    names(d) <- paste0("y", 1:12)
    for (j in seq(2, 12, by = 2)) d[[j]][stats::runif(50) < 0.3] <- NA
  }
  d
}

test_that("an unrestricted estimate that does not exist is refused at once, in DMAR's words", {
  skip_if_not_installed("lavaan")
  # The 1.1.0 review: with five complete cases of twelve items the
  # likelihood of the unrestricted model has no maximum, and lavaan's EM
  # drifts toward a singular covariance matrix. Run to 1e-10, the refit's
  # EM spent its 10,000 iterations following the drift (7 to 14 seconds,
  # against .4 where the estimate exists), passed on lavaan's warning to
  # raise max_iter, which cannot help, and returned FMIs that were those
  # of the point where it stopped (.168 for the loading of y2 at 10,000
  # iterations, .192 where it stalled near 30,000).
  d <- .fmi_sparse_draw(4)
  expect_identical(sum(stats::complete.cases(d)), 5L)
  model <- paste("f =~", paste(names(d), collapse = " + "))
  fit <- suppressWarnings(lavaan::cfa(model, data = d, missing = "ml",
                                      std.lv = TRUE))
  said <- new.env()
  said$w <- character(0)
  collect <- function(w) {
    said$w <- c(said$w, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  took <- system.time(err <- withCallingHandlers(
    tryCatch(fraction_missing_information(fit), error = identity),
    warning = collect))[["elapsed"]]
  expect_s3_class(err, "dmar_fmi_saturated_singular")
  msg <- conditionMessage(err)
  expect_match(msg, "Only 5 of the 50 cases observe all 12 variables",
               fixed = TRUE)
  expect_match(msg, "has no finite maximum", fixed = TRUE)
  expect_match(msg, "so none is computed", fixed = TRUE)
  expect_false(any(grepl("Maximum number of iterations", said$w,
                         fixed = TRUE)))
  expect_lt(took, 5)
  # The missing information degrees of freedom of fa_confirmatory() stop with the
  # same cause and the keyword's alternatives, without lavaan's notice
  # about the internal refit.
  said$w <- character(0)
  err <- withCallingHandlers(
    tryCatch(fa_confirmatory(d, list(f = names(d)), missing = "fiml",
                   wald_df = "missing_information"), error = identity),
    warning = collect)
  expect_s3_class(err, "dmar_fmi_saturated_singular")
  expect_match(conditionMessage(err), "has no finite maximum", fixed = TRUE)
  expect_match(conditionMessage(err), "Pass a number of degrees of freedom",
               fixed = TRUE)
  expect_false(any(grepl("Maximum number of iterations", said$w,
                         fixed = TRUE)))
  # An interior EM stationary point is not a global-maximum certificate.
  # This draw's eight complete rows admit a full-support collapse in twelve
  # dimensions, even though the local covariance is well conditioned.
  d2 <- .fmi_sparse_draw(2)
  fit2 <- suppressWarnings(lavaan::cfa(model, data = d2, missing = "ml",
                                       std.lv = TRUE))
  expect_error(fraction_missing_information(fit2),
    class = "dmar_fmi_saturated_singular")
  certificate <- DMAR:::.fmi_unbounded_subset(d2)
  expect_true(certificate$certified_by_dimension)
  expect_equal(certificate$n, 8)

})

test_that("EM iteration limits do not turn an unbounded H1 into a valid estimate", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # One data set of the calibration study's generator (loadings .8, .7, .7,
  # .6, .6, .5 twice, N = 50, 30 percent of the values of items 7 to 12
  # missing completely at random): the unrestricted estimate exists, and
  # lavaan's EM needs about 12,500 iterations to reach it. lavaan's notice
  # said to raise max_iter in em_h1_args without saying of what fit.
  set.seed(20261102)
  lam <- rep(c(0.8, 0.7, 0.7, 0.6, 0.6, 0.5), 2)
  for (r in 1:27) {
    f <- stats::rnorm(50)
    e <- matrix(stats::rnorm(50 * 12), 50, 12)
    y <- outer(f, lam) + sweep(e, 2L, sqrt(1 - lam^2), `*`)
    colnames(y) <- paste0("y", 1:12)
    for (j in 7:12) y[stats::runif(50) < 0.3, j] <- NA
  }
  d <- as.data.frame(y)
  fit <- lavaan::cfa(paste("f =~", paste(names(d), collapse = " + ")),
                     data = d, missing = "ml", std.lv = TRUE)
  said <- new.env()
  said$w <- character(0)
  # The old fixture was mistakenly called a finite maximum. The eight
  # complete rows permit a global affine collapse; raising the iteration
  # limit can only improve a local stationary point, not establish H1.
  expect_error(fraction_missing_information(fit),
    class = "dmar_fmi_saturated_singular")
  expect_true(DMAR:::.fmi_unbounded_subset(d)$certified_by_dimension)
  # A capped near-singular iterate without a certificate remains a
  # distinct fitting failure, not proof of an unbounded likelihood.
  z <- rbind(c(-1, 0), c(1, 0), c(-1000, NA), c(1000, NA),
             c(NA, -sqrt(2)), c(NA, sqrt(2)))
  colnames(z) <- c("x", "y")
  singular <- matrix(c(1, 1 - 1e-8, 1 - 1e-8, 1), 2)
  expect_error(DMAR:::.fmi_check_saturated(list(singular), list(z), capped = TRUE),
    "fitting failure, not a proof")

})

test_that("a pair of variables never observed together is refused by name", {
  skip_if_not_installed("lavaan")
  # A two-form design: t1 is missing on one form and t6 on the other. The
  # unrestricted estimate has no information about their covariance, and
  # lavaan's EM leaves it at 0, its starting value; FMIs computed from it
  # moved by up to .12 across values of it that the data could not tell
  # apart, with only lavaan's warning (given twice under fa_confirmatory()'s
  # keyword, once by the fit and once by the refit).
  d <- holzinger_swineford
  set.seed(3)
  form <- sample(1:2, nrow(d), replace = TRUE)
  d$t1_visual_perception[form == 1] <- NA
  d$t6_paragraph_comprehension[form == 2] <- NA
  f <- list(spatial = .fmi_spatial[1:4],
            verbal = c("t6_paragraph_comprehension", "t7_sentence",
                       "t9_word_meaning"))
  fit <- suppressWarnings(fa_confirmatory(d, f, missing = "fiml", output = "fit"))
  err <- tryCatch(fraction_missing_information(fit), error = identity)
  expect_s3_class(err, "dmar_fmi_never_together")
  expect_match(conditionMessage(err),
               "t1_visual_perception with t6_paragraph_comprehension",
               fixed = TRUE)
  expect_match(conditionMessage(err), "as in a planned missing design",
               fixed = TRUE)
})

test_that("a fit with auxiliary variables keeps its estimates in the refits", {
  skip_if_not_installed("lavaan")
  # The auxiliary fit of vignette("missing_data"): sentence completion
  # missing at random given word classification, which enters as a
  # saturated correlate through lavaan's aux argument. Word
  # classification shares the verbal factor, so its covariances with the
  # verbal tests run through their residuals, the augmented residual
  # covariance matrix is not positive definite, and the starting values
  # imply correlations above 1. lavaan's start check then shrank the
  # covariances in the refit, the estimates moved, and the computation
  # was refused as improper although no variance is negative.
  three <- list(
    visual = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
    verbal = c("t6_paragraph_comprehension", "t7_sentence",
               "t9_word_meaning"),
    speed  = c("t10_addition", "t12_counting_groups_of_dots",
               "t13_straight_and_curved_capitals"))
  nine <- unlist(three, use.names = FALSE)
  set.seed(113)
  d <- holzinger_swineford[c(nine, "t8_word_classification")]
  z_word <- as.numeric(scale(d$t8_word_classification))
  d$t7_sentence[stats::runif(nrow(d)) < stats::plogis(-0.7 - 1.5 * z_word)] <-
    NA
  fit <- suppressWarnings(suppressMessages(
    fa_confirmatory(data = d, factors = three, missing = "fiml",
          aux = "t8_word_classification", output = "fit")))
  pt <- lavaan::parTable(fit)
  expect_true(all(pt$est[pt$op == "~~" & pt$lhs == pt$rhs] > 0))
  expect_no_warning(res <- fraction_missing_information(fit))
  keep <- pt$free > 0L | pt$op == ":="
  expect_identical(res$estimate, pt$est[keep])
  expect_true(all(is.finite(res$fmi)))
  expect_true(all(res$fmi > -1e-8 & res$fmi < 1))
  expect_true("t8_word_classification ~~ t7_sentence" %in% res$term)
  # The FIML side is the fit's own: refit with the unstructured
  # information, lavaan's standard errors are the observed ones here.
  fit_h1 <- suppressWarnings(suppressMessages(
    fa_confirmatory(data = d, factors = three, missing = "fiml",
          aux = "t8_word_classification", output = "fit",
          observed.information = "h1", h1.information = "unstructured")))
  expect_equal(res$se_observed, lavaan::parTable(fit_h1)$se[keep],
               tolerance = 1e-6)
  expect_equal(res$fmi[res$term == "t7_sentence ~1"], 0.2357, tolerance = 1e-3)
  # With complete data the same specification has no missing
  # information.
  full <- suppressWarnings(suppressMessages(
    fa_confirmatory(data = holzinger_swineford[names(d)], factors = three,
          missing = "fiml", aux = "t8_word_classification",
          output = "fit")))
  res_full <- fraction_missing_information(full)
  expect_lt(max(abs(res_full$fmi)), 1e-8)
})

test_that("an equality constraint written through defined parameters links its rows", {
  skip_if_not_installed("lavaan")
  set.seed(4)
  d <- lavaan::HolzingerSwineford1939[c("x1", "x4", "x7", "school")]
  d$x4[stats::runif(301) < 0.3] <- NA
  d$x7[stats::runif(301) < 0.2] <- NA
  base <- "x4 ~ c(a1, a2)*x1
           x7 ~ c(b1, b2)*x4 + x1
           ind1 := a1*b1
           ind2 := a2*b2
           d1 := a1
           d2 := a2
           "
  fit_n <- function(con) {
    fit <- suppressWarnings(lavaan::sem(paste0(base, con), data = d,
                                        group = "school", missing = "ml",
                                        fixed.x = FALSE))
    res <- fraction_missing_information(fit)
    res[res$op %in% c("~", ":="), c("term", "group", "fmi", "N")]
  }
  # The same constraint written through the defined names and through
  # the labels links the same rows: the slopes it ties rest on both
  # schools.
  by_def <- fit_n("ind1 == ind2")
  by_lab <- fit_n("a1*b1 == a2*b2")
  expect_identical(by_def$N, by_lab$N)
  expect_equal(by_def$fmi, by_lab$fmi, tolerance = 1e-8)
  tied <- by_def$term %in% c("x4 ~ x1", "x7 ~ x4", "ind1", "ind2")
  expect_true(all(by_def$N[tied] == 301))
  expect_true(all(by_def$N[by_def$term == "x7 ~ x1"] %in% c(156, 145)))
  # A linear equality through defined names, d1 == d2, is a1 == a2.
  lin_def <- fit_n("d1 == d2")
  lin_lab <- fit_n("a1 == a2")
  expect_identical(lin_def$N, lin_lab$N)
  expect_true(all(lin_def$N[lin_def$term %in% c("x4 ~ x1", "d1", "d2")] ==
                    301))
  expect_true(all(lin_def$N[lin_def$term == "x7 ~ x4"] %in% c(156, 145)))
})

test_that("arguments are validated", {
  skip_if_not_installed("lavaan")
  fit <- .fmi_spatial_fit()
  expect_error(fraction_missing_information(fit, information = "expected"),
               "should be one of")
  expect_error(fraction_missing_information(fit, conf_level = 1.2),
               "'conf_level'")
  expect_error(fraction_missing_information(fit, ci_method = "percentile",
                                            B = 1),
               "single whole number", class = "dmar_count_invalid")
})

test_that("failed bootstrap replications are dropped with one warning", {
  skip_if_not_installed("lavaan")
  fit <- .fmi_spatial_fit(items = .fmi_spatial[c(1, 2, 5, 6)])
  point <- fraction_missing_information(fit)
  # lavBootstrap() is replaced by a stand-in that returns the point
  # estimates with small perturbations and marks three replications as
  # failed (NA) and one as non-finite, as lavaan does for a failed refit.
  local_mocked_bindings(
    lavBootstrap = function(object, r, fun, ...) {
      v <- fun(object)
      out <- matrix(rep(v, each = r), nrow = r) +
        matrix(seq_len(r) / 1e4, nrow = r, ncol = length(v))
      out[1:3, ] <- NA
      out[4, 2] <- Inf
      out
    },
    .package = "lavaan")
  expect_warning(
    res <- fraction_missing_information(fit, ci_method = "percentile",
                                        B = 100),
    "4 of 100 bootstrap replications")
  expect_identical(attr(res, "B"), 100L)
  expect_identical(attr(res, "B_used"), 96L)
  expect_identical(attr(res, "conf_level"), 0.95)
  expect_identical(names(res), names(point))
  expect_true(all(res$fmi_ci_lower < res$fmi_ci_upper))
  expect_output(print(res), "100 replications attempted, 96 returned")
})

test_that("the percentile interval is reproducible and leaves the caller's random stream alone", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(113)
  d <- holzinger_swineford[.fmi_spatial[1:4]]
  d$t3_paper_form_board[stats::runif(nrow(d)) < .25] <- NA
  d$t4_lozenges[stats::runif(nrow(d)) < .25] <- NA
  fit <- lavaan::cfa(paste("f =~", paste(names(d), collapse = " + ")),
                     data = d, std.lv = TRUE, missing = "ml")
  set.seed(1)
  before <- .Random.seed
  r1 <- fraction_missing_information(fit, ci_method = "percentile",
                                     B = 100, seed = 7)
  expect_identical(.Random.seed, before)
  r2 <- fraction_missing_information(fit, ci_method = "percentile",
                                     B = 100, seed = 7)
  expect_identical(r1$fmi_ci_lower, r2$fmi_ci_lower)
  expect_identical(r1$fmi_ci_upper, r2$fmi_ci_upper)
  expect_lte(attr(r1, "B_used"), 100L)
  expect_true(all(r1$fmi_ci_lower <= r1$fmi_ci_upper))
  load <- r1$op == "=~"
  expect_true(all(r1$fmi_ci_lower[load] <= r1$fmi[load] &
                    r1$fmi[load] <= r1$fmi_ci_upper[load]))
  expect_output(print(r1), "Confidence level: 95%")
  sentence <- results_sentence(generics::tidy(r1)[1, ])
  expect_match(sentence, "^fmi = [0-9.]+, 95% CI \\[")
})
