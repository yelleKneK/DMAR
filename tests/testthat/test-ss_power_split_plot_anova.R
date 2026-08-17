## Tests for ss_power_split_plot_anova(), the between x within (split-plot)
## power planner. Its between-subjects test is the compound-symmetry special
## case of the two-level model, so it is checked against ss_power_mixed_effects
## and a split-plot aov() simulation.

test_that("ss_power_split_plot_anova() between-subjects test loses power as rho increases", {
  low_rho  <- ss_power_split_plot_anova(a = 2, b = 4, effect = "between", f = 0.25, rho = 0.1, desired_power = 0.80)
  high_rho <- ss_power_split_plot_anova(a = 2, b = 4, effect = "between", f = 0.25, rho = 0.6, desired_power = 0.80)
  expect_lt(low_rho[low_rho$term == "necessary_n_per_group", "value"],
            high_rho[high_rho$term == "necessary_n_per_group", "value"])
})

test_that("ss_power_split_plot_anova() within-subjects test gains power as rho increases", {
  low_rho  <- ss_power_split_plot_anova(a = 2, b = 4, effect = "within", f = 0.25, rho = 0.1, desired_power = 0.80)
  high_rho <- ss_power_split_plot_anova(a = 2, b = 4, effect = "within", f = 0.25, rho = 0.6, desired_power = 0.80)
  expect_gt(low_rho[low_rho$term  == "necessary_n_per_group", "value"],
            high_rho[high_rho$term == "necessary_n_per_group", "value"])
})

test_that("ss_power_split_plot_anova() interaction test gains power as rho increases", {
  low_rho  <- ss_power_split_plot_anova(a = 2, b = 4, effect = "interaction", f = 0.25, rho = 0.1, desired_power = 0.80)
  high_rho <- ss_power_split_plot_anova(a = 2, b = 4, effect = "interaction", f = 0.25, rho = 0.6, desired_power = 0.80)
  expect_gt(low_rho[low_rho$term  == "necessary_n_per_group", "value"],
            high_rho[high_rho$term == "necessary_n_per_group", "value"])
})

test_that("ss_power_split_plot_anova() within-subjects noncentrality uses the total N", {
  # When J = 2 (smallest allowed), the within-subjects power calculation should still produce
  # something sensible and consistent with rm_anova for the same nominal n*J subjects.
  res <- ss_power_split_plot_anova(a = 2, b = 4, effect = "within", f = 0.25, rho = 0.5, n = 20)
  ncp_split <- res[res$term == "noncentrality", "value"]
  expect_equal(ncp_split, 2 * 20 * 4 * 0.25^2 / (1 - 0.5))
})

test_that("ss_power_split_plot_anova() Greenhouse-Geisser epsilon < 1 increases necessary n on within and interaction tests", {
  no_eps  <- ss_power_split_plot_anova(a = 2, b = 4, effect = "within", f = 0.25, rho = 0.5, epsilon = 1.0,  desired_power = 0.80)
  eps_low <- ss_power_split_plot_anova(a = 2, b = 4, effect = "within", f = 0.25, rho = 0.5, epsilon = 0.6, desired_power = 0.80)
  expect_gte(eps_low[eps_low$term == "necessary_n_per_group", "value"],
             no_eps[no_eps$term   == "necessary_n_per_group", "value"])
})

test_that("ss_power_split_plot_anova() validates arguments", {
  expect_error(ss_power_split_plot_anova(a = 1, b = 4, effect = "within", f = 0.25, rho = 0.5),
               "'a' must be a single integer >= 2")
  expect_error(ss_power_split_plot_anova(a = 2, b = 4, effect = "bogus", f = 0.25, rho = 0.5),
               "'effect' must be one of")
  expect_error(ss_power_split_plot_anova(a = 2, b = 4, effect = "within", f = 0.25, rho = 1.5),
               "'rho' must be a single numeric value in")
  expect_error(ss_power_split_plot_anova(a = 2, b = 4, effect = "within",
                                   f = 0.25, partial_eta_squared = 0.06, rho = 0.5),
               "Specify either 'f' or 'partial_eta_squared'")
})

test_that("ss_power_split_plot_anova() between-subjects noncentrality uses the total N", {
  # Regression guard. The between-subjects noncentrality is
  # lambda_B = N b f^2 / (1 + (b - 1) rho) with N = n a. An earlier version used the
  # per-group n in place of N, so the noncentrality was too small by the factor a
  # (halved for a = 2), understating between-subjects power. This mirrors the
  # within-subjects anchor above.
  a <- 3; b <- 4; n <- 15; f <- 0.25; rho <- 0.3
  res <- ss_power_split_plot_anova(a = a, b = b, effect = "between", f = f, rho = rho, n = n)
  expect_equal(res[res$term == "noncentrality", "value"],
               (n * a) * b * f^2 / (1 + (b - 1) * rho))
})

test_that("ss_power_split_plot_anova() between-subjects (a = 2) equals the mixed-effects treatment test", {
  # A between x within split-plot is the compound-symmetry (random-intercept) special
  # case of the two-level mixed model: rho is the intraclass correlation, the b
  # occasions are the n level-1 units, and the between-subjects factor is the
  # between-cluster treatment. For a = 2 the between-subjects F(1, .) test is therefore
  # the two-level mixed-effects two-sided treatment t test squared (F = t^2), the
  # identity that ties ss_power_split_plot_anova to ss_power_mixed_effects.
  d <- 0.5; b <- 4; n <- 20; rho <- 0.3
  sp <- ss_power_split_plot_anova(a = 2, b = b, effect = "between", f = d / 2, rho = rho, n = n)
  me <- ss_power_mixed_effects(d = d, n = b, rho = rho, J = n)
  expect_equal(sp[sp$term == "noncentrality", "value"],
               me[me$term == "noncentrality", "value"]^2)
  expect_equal(sp[sp$term == "actual_power", "value"],
               me[me$term == "actual_power", "value"], tolerance = 1e-6)
})

test_that("ss_power_split_plot_anova() matches a split-plot aov() simulation for all three effects", {
  skip_on_cran()
  set.seed(113)
  a <- 3; b <- 3; n <- 18; rho <- 0.3; fval <- 0.35; reps <- 1000
  sim_pow <- function(effect) {
    g <- rep(0, a); o <- rep(0, b); AB <- matrix(0, a, b)
    ca <- scale(seq_len(a), TRUE, FALSE)[, 1]; cb <- scale(seq_len(b), TRUE, FALSE)[, 1]
    if (effect == "between") g <- ca / sqrt(mean(ca^2)) * fval
    if (effect == "within")  o <- cb / sqrt(mean(cb^2)) * fval
    if (effect == "interaction") { M <- outer(ca, cb); AB <- M / sqrt(mean(M^2)) * fval }
    subj <- factor(rep(seq_len(a * n), each = b))
    grp  <- factor(rep(rep(seq_len(a), each = n), each = b))
    occ  <- factor(rep(seq_len(b), times = a * n))
    gi <- as.integer(grp); oi <- as.integer(occ); si <- as.integer(subj)
    key <- switch(effect, between = "grp", within = "occ", interaction = "grp:occ")
    rej <- logical(reps)
    for (r in seq_len(reps)) {
      y <- g[gi] + o[oi] + AB[cbind(gi, oi)] +
        rnorm(a * n, 0, sqrt(rho))[si] + rnorm(a * n * b, 0, sqrt(1 - rho))
      ss <- summary(aov(y ~ grp * occ + Error(subj / occ)))
      p <- if (effect == "between") ss[["Error: subj"]][[1]]["grp", "Pr(>F)"]
           else ss[["Error: subj:occ"]][[1]][key, "Pr(>F)"]
      rej[r] <- isTRUE(p < 0.05)
    }
    mean(rej)
  }
  for (eff in c("between", "within", "interaction")) {
    res <- ss_power_split_plot_anova(a = a, b = b, effect = eff, f = fval, rho = rho, n = n)
    ana <- res[res$term == "actual_power", "value"]
    emp <- sim_pow(eff)
    expect_equal(ana, emp, tolerance = 4 * sqrt(emp * (1 - emp) / reps) + 0.01)
  }
})
