## ss_aipe_sc_sensitivity() applies the weights rule of ci_sc() and
## ss_aipe_sc() (finding 9 of the 2026-11-11 planner register), whether it
## plans the sample size or is given one.

test_that("ss_aipe_sc_sensitivity() refuses the weights ci_sc() refuses", {
  for (args in list(list(supposed_std_contrast = .5, width = .4),
                    list(specified_n_per_group = 20, width = .4))) {
    err <- tryCatch(do.call(ss_aipe_sc_sensitivity,
                            c(list(true_std_contrast = .5,
                                   c_weights = c(.25, .25, -.5), G = 3,
                                   print_iter = FALSE), args)),
                    error = function(e) e)
    expect_s3_class(err, "dmar_argument_invalid")
    expect_identical(err$arg, "c_weights")
    expect_match(conditionMessage(err), "must sum to 1", fixed = TRUE)
  }
  res <- ss_aipe_sc_sensitivity(true_std_contrast = .5, specified_n_per_group = 20,
                                c_weights = c(.5, .5, -1), width = .8, G = 3,
                                print_iter = FALSE, seed = 1)
  expect_identical(res$value[res$term == "n_per_group"], 20)
  # Integer weights are judged by their values, as in ss_aipe_sc().
  int <- ss_aipe_sc_sensitivity(true_std_contrast = .5, specified_n_per_group = 20,
                                c_weights = -1:1, width = .8, G = 3,
                                print_iter = FALSE, seed = 1)
  dbl <- ss_aipe_sc_sensitivity(true_std_contrast = .5, specified_n_per_group = 20,
                                c_weights = c(-1, 0, 1), width = .8, G = 3,
                                print_iter = FALSE, seed = 1)
  expect_identical(int$value, dbl$value)
})

## Each replication is the ci_sc() interval of its simulated data, with the
## pooled standard deviation (the root mean square error) as s_anova. DMAR
## 1.0.x divided by the mean of the group standard deviations.
test_that("ss_aipe_sc_sensitivity() records the ci_sc() interval of the data it simulates", {
  for (k in list(list(cw = c(.5, .5, -1), n = 10, d = .5, cl = .95),
                 list(cw = -1:1, n = 6, d = -.8, cl = .90),
                 list(cw = c(1, -1), n = 25, d = 1.2, cl = .99))) {
    f <- tempfile(fileext = ".csv")
    withr::defer(unlink(f))
    G <- 15
    res <- ss_aipe_sc_sensitivity(true_std_contrast = k$d, specified_n_per_group = k$n,
                                  c_weights = k$cw, width = 1, conf_level = k$cl,
                                  G = G, print_iter = FALSE, filename = f, seed = 21)
    rec <- utils::read.csv(f)
    expect_identical(nrow(rec), as.integer(G))
    # The same draws: J groups of n, SD 1, the contrast placed on the first
    # positive weight.
    J <- length(k$cw)
    mu <- rep(0, J)
    first <- which(k$cw > 0)[1]
    mu[first] <- k$d / k$cw[first]
    by_ci_sc <- withr::with_seed(21, t(vapply(seq_len(G), function(i) {
      y <- lapply(seq_len(J), function(j) stats::rnorm(k$n, mu[j], 1))
      mse <- sum(vapply(y, function(v) sum((v - mean(v))^2), 0)) / (J * k$n - J)
      ci_sc(means = vapply(y, mean, 0), s_anova = sqrt(mse), c_weights = k$cw,
            n_per_group = k$n, conf_level = k$cl)$value
    }, numeric(3))))
    expect_equal(rec$ci_lower, by_ci_sc[, 1], tolerance = 1e-10)
    expect_equal(rec$std_contrast, by_ci_sc[, 2], tolerance = 1e-10)
    expect_equal(rec$ci_upper, by_ci_sc[, 3], tolerance = 1e-10)
    v <- stats::setNames(res$value, res$term)
    expect_equal(unname(v["mean_std_contrast"]), mean(by_ci_sc[, 2]), tolerance = 1e-10)
    expect_equal(unname(v["mean_ci_width"]), mean(by_ci_sc[, 3] - by_ci_sc[, 1]),
                 tolerance = 1e-10)
  }
})

test_that("ss_aipe_sc_sensitivity() gives no warning, whatever the assurance method", {
  # The only method offered, "analytic", given or left at its default, is
  # the one every call uses.
  for (args in list(list(specified_n_per_group = 10),
                    list(specified_n_per_group = 10, assurance_method = "analytic"))) {
    expect_no_warning(do.call(ss_aipe_sc_sensitivity, c(list(
      true_std_contrast = .5, c_weights = c(.5, .5, -1), width = 1, G = 2,
      print_iter = FALSE, seed = 1), args)))
  }
})
