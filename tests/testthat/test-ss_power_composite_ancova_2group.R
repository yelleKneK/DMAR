test_that("a single-term composite is the ordinary noncentral t power", {
  # The composite machinery integrates the conditional rejection probability
  # over the error distribution. With one test that integral must reproduce
  # pt() exactly, which is what keeps the marginals and the composite from
  # drifting apart.
  for (d in c(0.5, 1.5, 3)) {
    for (df in c(10, 50, 200)) {
      got <- .composite_power_shared_sigma(d, df, 0.05, FALSE)$composite
      want <- (1 - pt(qt(0.975, df), df, ncp = d)) +
        pt(qt(0.025, df), df, ncp = d)
      expect_equal(got, want, tolerance = 1e-10)
    }
  }
})

test_that("the group test reproduces ss_power_c_ancova", {
  # With the interaction dropped, the model is the one-way ANCOVA that
  # ss_power_c_ancova plans for, and the group effect is the c(1, -1) contrast.
  # The two functions must agree to numerical precision.
  grid <- list(c(0.5, 0.3, 30), c(0.2, 0.1, 100), c(0.8, 0.6, 15))
  for (g in grid) {
    mine <- ss_power_composite_ancova_2group(smd = g[1], rho = g[2], n = g[3],
                                      composite_terms = "group",
                                      include_interaction = FALSE)
    theirs <- ss_power_c_ancova(psi = g[1], c_weights = c(1, -1), sigma = 1,
                                rho = g[2], n = g[3])
    expect_equal(mine$value[mine$term == "composite_power"],
                 theirs$value[theirs$term == "actual_power"],
                 tolerance = 1e-10)
  }
})

test_that("composite power lies between the product and the minimum of the marginals", {
  # Sharing an error estimate makes the tests positively dependent, so the
  # composite exceeds the product of the marginals; it can never exceed the
  # least powerful test on its own.
  set.seed(113)
  for (i in 1:50) {
    delta <- runif(3, 0, 3)
    df <- sample(5:300, 1)
    res <- .composite_power_shared_sigma(delta, df, 0.05, FALSE)
    expect_gte(res$composite, prod(res$marginal) - 1e-10)
    expect_lte(res$composite, min(res$marginal) + 1e-10)
  }
})

test_that("the composite approaches the product of the marginals as df grows", {
  # The dependence among the tests comes entirely from the shared error
  # estimate. As the residual degrees of freedom grow that estimate stabilizes,
  # the tests decouple, and the composite converges to the product.
  excess <- vapply(c(8, 20, 60, 500, 5000), function(df) {
    res <- .composite_power_shared_sigma(c(1.9, 2.1), df, 0.05, FALSE)
    res$composite - prod(res$marginal)
  }, numeric(1))
  expect_true(all(excess > 0))
  expect_true(all(diff(excess) < 0))
  expect_lt(excess[length(excess)], 1e-3)
})

test_that("equal correlations give equal slopes and no interaction", {
  out <- ss_power_composite_ancova_2group(smd = 0.2, rho = 0.3, n = 50)
  slopes <- attr(out, "slopes")
  expect_equal(unname(slopes[1]), unname(slopes[2]))
  expect_equal(unname(attr(out, "coefficients")["group_by_covariate"]), 0)
  # A null interaction rejects at the Type I error rate.
  expect_equal(out$value[out$term == "power_group_by_covariate"], 0.05,
               tolerance = 1e-6)
})

test_that("unequal correlations put a nonzero interaction in the model", {
  out <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.2, rho = c(0.1, 0.4), n = 100))
  slopes <- attr(out, "slopes")
  expect_gt(slopes[["group_2"]], slopes[["group_1"]])
  expect_equal(unname(attr(out, "coefficients")["group_by_covariate"]),
               unname(slopes[["group_2"]] - slopes[["group_1"]]))
  # The correlations also differ in absolute value, so the residual variances
  # differ and the powers are reported under their approximate names.
  expect_gt(out$value[out$term == "approximate_power_group_by_covariate"], 0.05)
})

test_that("equal correlations recover sigma_adj = sigma * sqrt(1 - rho^2)", {
  # The general parameterization averages the within-group outcome variances.
  # When the slopes are equal it must collapse to the familiar ANCOVA form.
  for (rho in c(0, 0.2, 0.5, 0.9)) {
    out <- ss_power_composite_ancova_2group(smd = 0.2, rho = rho, sigma = 2, n = 40)
    expect_equal(attr(out, "sigma_adj"), 2 * sqrt(1 - rho^2), tolerance = 1e-12)
  }
})

test_that("sigma is the within-group SD of the outcome in both groups", {
  # sigma is the anchor: it is the same in both groups whether or not the slopes
  # differ. Within group g the slope is rho_g * sigma / sd_cov and the residual
  # variance is sigma^2 * (1 - rho_g^2), so unequal correlations give unequal
  # slopes without sigma moving.
  rho <- c(0.1, 0.6)
  sigma <- 1.5
  d <- .ancova_design(smd = 0.2, rho = rho, sigma = sigma, sd_cov = 2,
                      n = 50, include_interaction = TRUE)
  expect_equal(d$slope_1, rho[1] * sigma / 2, tolerance = 1e-12)
  expect_equal(d$slope_2, rho[2] * sigma / 2, tolerance = 1e-12)
  # The pooled error variance the ANCOVA estimates is the average of the two
  # within-group residual variances.
  expect_equal(d$sigma_adj^2,
               mean(c(sigma^2 * (1 - rho[1]^2), sigma^2 * (1 - rho[2]^2))),
               tolerance = 1e-12)
})

test_that("sd_cov does not affect power", {
  # A correlation is scale free, so the spread of the covariate cannot change
  # the noncentralities. It sets the units the slopes are reported in.
  a <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.3, rho = c(0.1, 0.4), n = 60,
                                     sd_cov = 1))
  b <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.3, rho = c(0.1, 0.4), n = 60,
                                     sd_cov = 17))
  # These correlations differ in absolute value, so the power rows carry the
  # approximate names. Naming the row explicitly, rather than reading whichever
  # row happens to be there, is what keeps the comparison from succeeding on a
  # pair of zero-length vectors.
  expect_length(a$value[a$term == "approximate_composite_power"], 1L)
  expect_equal(a$value[a$term == "approximate_composite_power"],
               b$value[b$term == "approximate_composite_power"],
               tolerance = 1e-12)
  # but the slopes are in the covariate's units, so they do move
  expect_equal(attr(b, "slopes")[["group_1"]],
               attr(a, "slopes")[["group_1"]] / 17, tolerance = 1e-12)
})

test_that("the interaction costs one residual degree of freedom", {
  with_int <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.2, rho = c(0.1, 0.3), n = 50))
  no_int <- ss_power_composite_ancova_2group(smd = 0.2, rho = 0.1, n = 50,
                                      composite_terms = c("group", "covariate"),
                                      include_interaction = FALSE)
  expect_equal(with_int$value[with_int$term == "residual_df"], 2 * 50 - 4)
  expect_equal(no_int$value[no_int$term == "residual_df"], 2 * 50 - 3)
})

test_that("planning a sample size attains the desired composite power", {
  # Equal correlations are the exact case: the size is necessary_n_per_group and
  # the power it attains is composite_power.
  out <- ss_power_composite_ancova_2group(smd = 0.5, rho = 0.4,
                                   composite_terms = c("group", "covariate"),
                                   desired_power = 0.80)
  n <- out$value[out$term == "necessary_n_per_group"]
  expect_gte(out$value[out$term == "composite_power"], 0.80)
  # and the sample size is the smallest one that does
  below <- ss_power_composite_ancova_2group(smd = 0.5, rho = 0.4, n = n - 1,
                                     composite_terms = c("group", "covariate"))
  expect_lt(below$value[below$term == "composite_power"], 0.80)
})

test_that("an unequal-correlation plan reports an approximate size, not a necessary one", {
  # The same search runs, but the target is reached by an approximation, so the
  # size is not presented as one that attains it.
  out <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.5, rho = c(0.1, 0.5),
                                     composite_terms = c("group",
                                                         "group_by_covariate"),
                                     desired_power = 0.80))
  expect_false("necessary_n_per_group" %in% out$term)
  expect_false("necessary_N" %in% out$term)
  n <- out$value[out$term == "approximate_n_per_group"]
  expect_length(n, 1L)
  expect_gte(out$value[out$term == "approximate_composite_power"], 0.80)
  below <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.5, rho = c(0.1, 0.5), n = n - 1,
                                     composite_terms = c("group",
                                                         "group_by_covariate")))
  expect_lt(below$value[below$term == "approximate_composite_power"], 0.80)
})

test_that("a supposed effect of zero stops the sample size search", {
  # Power for a null effect is the Type I error rate at every sample size, so
  # no sample size attains the desired power and the search would not end.
  expect_error(
    ss_power_composite_ancova_2group(smd = 0.2, rho = 0.1, desired_power = 0.8),
    "supposed effect of zero")
  # The unequal |rho| pair fires the shared-error approximation warning
  # before the zero-effect error stops the search; assert both.
  expect_warning(
    expect_error(
      ss_power_composite_ancova_2group(smd = 0, rho = c(0.1, 0.4), desired_power = 0.8,
                                composite_terms = "group"),
      "supposed effect of zero"),
    "differ in absolute value")
})

test_that("an effect cannot be tested unless it is in the model", {
  expect_error(
    ss_power_composite_ancova_2group(composite_terms = "group_by_covariate",
                              include_interaction = FALSE, n = 50),
    "cannot be tested unless it is in the model")
})

test_that("the defaults describe a population with no effects", {
  out <- ss_power_composite_ancova_2group(n = 100)
  expect_equal(out$value[out$term == "supposed_smd"], 0)
  expect_equal(out$value[out$term == "supposed_rho_group_1"], 0)
  expect_equal(out$value[out$term == "power_group"], 0.05, tolerance = 1e-6)
  expect_equal(out$value[out$term == "noncentral_t_parm_group"], 0)
})

test_that("directional tests are more powerful in the supposed direction", {
  two <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.3, rho = c(0.1, 0.4), n = 60))
  one <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.3, rho = c(0.1, 0.4), n = 60,
                                     directional = TRUE))
  expect_gt(one$value[one$term == "approximate_composite_power"],
            two$value[two$term == "approximate_composite_power"])
  expect_equal(one$value[one$term == "tails"], 1)
  expect_equal(two$value[two$term == "tails"], 2)
})

test_that("the value column is numeric and the term names are attributes", {
  out <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.2, rho = c(0.1, 0.4), n = 50))
  expect_true(is.numeric(out$value))
  expect_identical(attr(out, "composite_terms"),
                   c("group", "covariate", "group_by_covariate"))
  expect_s3_class(out, "dmar_tbl")
  expect_s3_class(out, "dmar_ss_power")
})

test_that("the reported order of the composite terms is fixed", {
  # The order the terms are typed in must not change the returned table.
  a <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.2, rho = c(0.1, 0.4), n = 50,
                                     composite_terms = c("group_by_covariate",
                                                         "group")))
  b <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.2, rho = c(0.1, 0.4), n = 50,
                                     composite_terms = c("group",
                                                         "group_by_covariate")))
  expect_identical(a$term, b$term)
  expect_equal(a$value, b$value)
})

test_that("planning values are rejected when they are out of range", {
  expect_error(ss_power_composite_ancova_2group(rho = 1, n = 50), "must lie in")
  expect_error(ss_power_composite_ancova_2group(rho = c(0.1, 0.2, 0.3), n = 50),
               "length 1")
  expect_error(ss_power_composite_ancova_2group(sigma = 0, n = 50), "'sigma'")
  expect_error(ss_power_composite_ancova_2group(sd_cov = -1, n = 50), "'sd_cov'")
  expect_error(ss_power_composite_ancova_2group(alpha_level = 0, n = 50),
               "'alpha_level'")
  expect_error(ss_power_composite_ancova_2group(n = 1), "'n'")
  expect_error(ss_power_composite_ancova_2group(n = 50.5), "'n'")
  expect_error(ss_power_composite_ancova_2group(smd = 0.2, rho = 0.1, n = 50,
                                         composite_terms = "bogus"),
               "not in the model")
  expect_error(ss_power_composite_ancova_2group(smd = 0.2, rho = 0.1, n = 50,
                                         composite_terms = c("group", "group")),
               "duplicate")
})

test_that("dropping the interaction while the slopes differ is refused", {
  # Correlations that differ put a nonzero interaction in the population. A
  # model that omits it absorbs it into the error, so planning with the error
  # variance of the model that keeps the term overstates every power in the
  # table. The combination is refused rather than silently mis-planned.
  expect_error(
    ss_power_composite_ancova_2group(smd = 0.5, rho = c(0.1, 0.9), n = 50,
                              composite_terms = c("group", "covariate"),
                              include_interaction = FALSE),
    "differs across the groups")
  # Equal correlations are the case ss_power_c_ancova plans for and stay allowed,
  # including when given as a length-two vector.
  expect_silent(ss_power_composite_ancova_2group(smd = 0.5, rho = 0.3, n = 50,
                                          composite_terms = "group",
                                          include_interaction = FALSE))
  expect_silent(ss_power_composite_ancova_2group(smd = 0.5, rho = c(0.3, 0.3), n = 50,
                                          composite_terms = "group",
                                          include_interaction = FALSE))
})

test_that("a length-two rho of equal values matches the length-one form", {
  a <- ss_power_composite_ancova_2group(smd = 0.3, rho = 0.25, n = 40)
  b <- ss_power_composite_ancova_2group(smd = 0.3, rho = c(0.25, 0.25), n = 40)
  expect_equal(a$value, b$value)
})

test_that("plot() draws the population effects from a result in hand", {
  skip_if_not_installed("ggplot2")
  # The figure is reached from a result, for both an n-specified result and a
  # desired-power plan, and shape does not depend on an argument.
  # The plot method reads the powers and the size out of the table, so it has to
  # work whichever set of row names the table carries. The first two calls are
  # the approximate case and the last the exact one.
  expect_s3_class(
    plot(suppressWarnings(
      ss_power_composite_ancova_2group(smd = 0.2, rho = c(0.1, 0.4), n = 100))),
    "ggplot")
  expect_s3_class(
    plot(suppressWarnings(
      ss_power_composite_ancova_2group(smd = 0.5, rho = c(0.1, 0.5),
                                       desired_power = 0.80))),
    "ggplot")
  expect_s3_class(
    plot(ss_power_composite_ancova_2group(smd = 0.5, rho = 0.4,
                                composite_terms = c("group", "covariate"),
                                desired_power = 0.80)),
    "ggplot")
  # show_power = FALSE draws the effects without the power annotation.
  expect_s3_class(
    plot(suppressWarnings(
      ss_power_composite_ancova_2group(smd = 0.2, rho = c(0.1, 0.4), n = 100)),
         show_power = FALSE),
    "ggplot")
})

test_that("unequal absolute correlations warn that the composite power is approximate (HIGH-01)", {
  # Unequal residual variances make the shared-error composite an approximation
  # whose error grows with the correlation gap; the planner must not report it
  # as exact. Equal absolute correlations are exact and stay silent.
  expect_warning(
    ss_power_composite_ancova_2group(smd = 0.3, rho = c(0, 0.9), n = 20,
                                     composite_terms = c("group", "covariate")),
    "approximation")
  expect_no_warning(
    ss_power_composite_ancova_2group(smd = 0.3, rho = 0.3, n = 20,
                                     composite_terms = c("group", "covariate")))
})

test_that("equal absolute correlations keep the exact term names", {
  # The common case must not move. Equal correlations, and correlations of the
  # same magnitude and opposite sign, both leave the residual variances equal,
  # so both are exact and neither is relabeled.
  for (rho in list(0.3, c(0.25, 0.25), c(0.4, -0.4))) {
    at_n <- ss_power_composite_ancova_2group(smd = 0.3, rho = rho, n = 50)
    expect_true(all(c("specified_n_per_group", "specified_N",
                      "composite_power", "power_group", "power_covariate",
                      "power_group_by_covariate") %in% at_n$term))
    expect_false(any(grepl("^approximate_", at_n$term)))
    expect_false(attr(at_n, "approximate"))
  }
  planned <- ss_power_composite_ancova_2group(smd = 0.5, rho = 0.4,
                                composite_terms = c("group", "covariate"),
                                desired_power = 0.80)
  expect_true(all(c("necessary_n_per_group", "necessary_N", "composite_power",
                    "power_group", "power_covariate") %in% planned$term))
  expect_false(any(grepl("^approximate_", planned$term)))
})

test_that("unequal absolute correlations relabel the powers and the planned size", {
  # The output itself has to say the number is an approximation, because a
  # reader pulling a row out of the table never sees the warning.
  at_n <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.2, rho = c(0.1, 0.4), n = 100))
  expect_true(all(c("approximate_composite_power", "approximate_power_group",
                    "approximate_power_covariate",
                    "approximate_power_group_by_covariate") %in% at_n$term))
  expect_false(any(c("composite_power", "power_group", "power_covariate",
                     "power_group_by_covariate") %in% at_n$term))
  expect_true(attr(at_n, "approximate"))
  # The rows that are not powers keep their names: a noncentrality, a residual
  # degrees of freedom, and an echoed planning value are what they always were.
  expect_true(all(c("specified_n_per_group", "specified_N", "residual_df",
                    "noncentral_t_parm_group", "supposed_smd", "tails")
                  %in% at_n$term))
})

test_that("relabeling is a renaming and not a recomputation", {
  # The values are the ones the planner reported before the labels existed,
  # recorded here to full precision so a change in the arithmetic cannot pass
  # as a change in the names.
  out <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.5, rho = c(0.1, 0.5),
                                     composite_terms = c("group",
                                                         "group_by_covariate"),
                                     desired_power = 0.80))
  v <- stats::setNames(out$value, out$term)
  expect_equal(unname(v[["approximate_n_per_group"]]), 95)
  expect_equal(unname(v[["approximate_N"]]), 190)
  expect_equal(unname(v[["approximate_composite_power"]]), 0.8006971221,
               tolerance = 1e-9)
  expect_equal(unname(v[["approximate_power_group"]]), 0.9568646626,
               tolerance = 1e-9)
  expect_equal(unname(v[["approximate_power_group_by_covariate"]]),
               0.8365443417, tolerance = 1e-9)
  expect_equal(unname(v[["noncentral_t_parm_group"]]), 3.6945128620,
               tolerance = 1e-9)
})

test_that("the approximate and exact term maps are inverses", {
  # The plot methods read a table by its exact names, so the mapping has to be
  # reversible on every row a composite table can carry.
  exact <- c("specified_n_per_group", "specified_N", "necessary_n_per_group",
             "necessary_n_per_cell", "necessary_N", "composite_power",
             "residual_df", "power_group", "power_covariate",
             "power_group_by_covariate", "power_cov_x_1",
             "noncentral_t_parm_group", "f_1", "df_1", "cells",
             "desired_power", "supposed_smd", "alpha_level", "tails")
  approx <- DMAR:::.composite_approximate_terms(exact)
  expect_identical(DMAR:::.composite_exact_terms(approx), exact)
  # and only the powers and the resolved sizes move
  moved <- exact[approx != exact]
  expect_identical(moved,
                   c("necessary_n_per_group", "necessary_n_per_cell",
                     "necessary_N", "composite_power", "power_group",
                     "power_covariate", "power_group_by_covariate",
                     "power_cov_x_1"))
})

test_that("tidy() and glance() still find the size and the power when relabeled", {
  # The shared ss_power tidiers match fixed row names, so a relabeled table
  # would otherwise tidy to a pair of missing values without complaint.
  out <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.5, rho = c(0.1, 0.5),
                                     composite_terms = c("group",
                                                         "group_by_covariate"),
                                     desired_power = 0.80))
  td <- generics::tidy(out)
  expect_equal(td$estimate, out$value[out$term == "approximate_n_per_group"])
  expect_equal(td$power,
               out$value[out$term == "approximate_composite_power"])
  gl <- generics::glance(out)
  expect_equal(gl$estimate, td$estimate)
  expect_equal(gl$power, td$power)
  # The per-group size wins over the total, the same way it does for a table
  # that carries the exact names.
  expect_equal(td$estimate, 95)
  # and the exact case is still summarized the same way
  ok <- ss_power_composite_ancova_2group(smd = 0.5, rho = 0.4,
                                composite_terms = c("group", "covariate"),
                                desired_power = 0.80)
  expect_equal(generics::tidy(ok)$estimate,
               ok$value[ok$term == "necessary_n_per_group"])
  expect_equal(generics::tidy(ok)$power,
               ok$value[ok$term == "composite_power"])
})

test_that("the label is earned: simulation matches the exact case and not the approximate one", {
  # An external check that owes nothing to the composite machinery. The
  # covariate is held fixed at the population moments, so the design columns are
  # exactly orthogonal and the conditioning approximation plays no part; the
  # only thing separating the two calls below is whether the two groups share a
  # residual variance. Equal correlations reproduce the simulation to Monte
  # Carlo error. Unequal correlations do not, which is what the approximate_
  # labels are there to say.
  skip_on_cran()
  set.seed(113)
  simulate_composite <- function(smd, rho, n, G, sigma = 1, sd_cov = 1,
                                 alpha = 0.05) {
    x0 <- stats::qnorm((seq_len(n) - 0.5) / n)
    x0 <- x0 - mean(x0)
    x0 <- x0 / sqrt(mean(x0^2)) * sd_cov          # sum(x0^2) = n * sd_cov^2
    grp <- rep(c(-0.5, 0.5), each = n)
    cov <- rep(x0, 2)
    design <- cbind(1, grp, cov, grp * cov)
    N <- nrow(design); df <- N - ncol(design)
    xtxi <- solve(crossprod(design))
    se_mult <- sqrt(diag(xtxi))
    slope <- rho * sigma / sd_cov
    beta <- c(0, smd * sigma, mean(slope), slope[2] - slope[1])
    mu <- as.vector(design %*% beta)
    e_sd <- rep(sigma * sqrt(1 - rho^2), each = n)
    crit <- stats::qt(1 - alpha / 2, df)
    y <- mu + matrix(stats::rnorm(N * G, 0, e_sd), nrow = N)
    b <- xtxi %*% t(design) %*% y
    s <- sqrt(colSums((y - design %*% b)^2) / df)
    mean(abs(b[2, ]) / (s * se_mult[2]) > crit &
         abs(b[3, ]) / (s * se_mult[3]) > crit)
  }
  G <- 60000
  terms <- c("group", "covariate")

  exact <- ss_power_composite_ancova_2group(smd = 0.3, rho = 0.5, n = 8,
                                            composite_terms = terms)
  got <- exact$value[exact$term == "composite_power"]
  sim <- simulate_composite(0.3, c(0.5, 0.5), 8, G)
  se <- sqrt(sim * (1 - sim) / G)
  expect_lt(abs(got - sim), 4 * se)

  approx <- suppressWarnings(
    ss_power_composite_ancova_2group(smd = 0.3, rho = c(0, 0.9), n = 8,
                                     composite_terms = terms))
  got <- approx$value[approx$term == "approximate_composite_power"]
  sim <- simulate_composite(0.3, c(0, 0.9), 8, G)
  se <- sqrt(sim * (1 - sim) / G)
  expect_gt(abs(got - sim), 5 * se)
})
