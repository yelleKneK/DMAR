# Informed priors for the Bayes t family: the Cauchy takes a location, and
# a normal prior takes the mean and standard deviation directly. The two
# families are linked by an exact identity: Cauchy(mu, r) is a normal
# N(mu, r^2 / z^2) whose z is standard normal (equivalently, the normal's
# variance is inverse-gamma(1/2, r^2/2)). The identity cross-validates the
# two code paths against each other with no external oracle.

test_that("the Cauchy marginal equals the normal-mixture marginal", {
  t_obs <- 2.3; df <- 28; n_eff <- 29
  for (loc in c(0, 0.4)) {
    r <- 0.6
    prior_c <- DMAR:::.bayes_t_prior(prior_location = loc, prior_scale = r)
    m_c <- DMAR:::.bayes_t_marginal(t_obs, df, n_eff, prior_c)
    mix_integrand <- function(z) {
      vapply(z, function(zz) {
        pr <- DMAR:::.bayes_t_prior(prior_mean = loc, prior_sd = r / zz)
        DMAR:::.bayes_t_marginal(t_obs, df, n_eff, pr)
      }, numeric(1)) * 2 * stats::dnorm(z)
    }
    m_mix <- stats::integrate(mix_integrand, 0, Inf, rel.tol = 1e-7)$value
    expect_equal(m_c, m_mix, tolerance = 1e-5)
  }
})

test_that("the default prior reproduces the original JZS results", {
  set.seed(113)
  x <- rnorm(25, 0.4)
  default_call <- bayes_one_sample_t(x)
  explicit     <- bayes_one_sample_t(x, prior_location = 0,
                                     prior_scale = sqrt(2) / 2)
  expect_equal(default_call$value, explicit$value, tolerance = 1e-12)
  expect_identical(attr(default_call, "prior_family"), "cauchy")
})

test_that("a diffuse normal prior recovers the data", {
  set.seed(113)
  x <- rnorm(40, 0.5)
  res <- bayes_one_sample_t(x, prior_mean = 0, prior_sd = 50)
  d_hat <- mean(x) / sd(x)
  med <- res$value[res$term == "delta_posterior_median"]
  expect_equal(med, d_hat, tolerance = 0.02)
  expect_identical(attr(res, "prior_family"), "normal")
})

test_that("an informed prior location moves the Bayes factor sensibly", {
  set.seed(113)
  x <- rnorm(20, 0.6)
  d_hat <- mean(x) / sd(x)
  bf_with    <- bayes_one_sample_t(x, prior_location = d_hat,
                                   prior_scale = 0.3)
  bf_against <- bayes_one_sample_t(x, prior_location = -d_hat,
                                   prior_scale = 0.3)
  expect_gt(bf_with$value[bf_with$term == "bf_10"],
            bf_against$value[bf_against$term == "bf_10"])
})

test_that("the posterior rides on the object and is a proper density", {
  set.seed(113)
  x <- rnorm(25, 0.4)
  res <- bayes_one_sample_t(x, prior_mean = 0.3, prior_sd = 0.5)
  post <- attr(res, "posterior")
  expect_s3_class(post, "data.frame")
  step <- diff(post$delta[1:2])
  expect_equal(sum(post$density) * step, 1, tolerance = 1e-6)
  # P(delta > 0) recomputed from the returned density matches the row.
  p_pos_row <- res$value[res$term == "p_delta_positive"]
  p_pos_att <- sum(post$density[post$delta > 0]) * step
  expect_equal(p_pos_att, p_pos_row, tolerance = 1e-3)
})

test_that("prior families are mutually exclusive and validated", {
  set.seed(113)
  x <- rnorm(10)
  expect_error(bayes_one_sample_t(x, prior_mean = 0.3),
               "both 'prior_mean' and 'prior_sd'")
  expect_error(bayes_one_sample_t(x, prior_mean = 0.3, prior_sd = 0.5,
                                  prior_scale = 1),
               "one prior family")
  expect_error(bayes_one_sample_t(x, prior_mean = 0.3, prior_sd = 0.5,
                                  prior_location = 0.2),
               "one prior family")
  expect_error(bayes_one_sample_t(x, prior_mean = 0, prior_sd = -1),
               "'prior_sd' must be positive")
  y <- rnorm(10)
  paired_raw  <- bayes_paired_t(x, y, prior_mean = 0.2, prior_sd = 0.4)
  paired_diff <- bayes_one_sample_t(x - y, prior_mean = 0.2, prior_sd = 0.4)
  expect_equal(paired_raw$value, paired_diff$value, tolerance = 1e-12)
})
