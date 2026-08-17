# Tests for ss_power_composite_sem(). Validation and setup errors run
# everywhere; the Monte Carlo paths are off CRAN for runtime. Small
# observed-variable path models keep each lavaan fit fast, and their
# population covariance matrices are known in closed form, so the empirical
# powers can be checked against the analytic Wald powers they estimate.

test_that("ss_power_composite_sem() validates its inputs", {
  skip_if_not_installed("lavaan")
  sigma_2 <- matrix(c(1, 0.5, 0.5, 1), 2, 2,
                    dimnames = list(c("y1", "y2"), c("y1", "y2")))
  pop_2 <- "y1 ~~ 1*y1
            y2 ~ 0.5*y1
            y2 ~~ 0.75*y2"

  # Exactly one way of stating the population.
  expect_error(ss_power_composite_sem("y2 ~ b*y1"), "exactly one")
  expect_error(ss_power_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                      pop_model = pop_2), "exactly one")

  # Sigma must carry the observed-variable names.
  expect_error(ss_power_composite_sem("y2 ~ b*y1",
                                      Sigma = matrix(c(1, .5, .5, 1), 2, 2)),
               "row and column names")

  # The parameters of interest must be labeled in the model.
  expect_error(ss_power_composite_sem("y2 ~ y1", Sigma = sigma_2),
               "No labeled parameters")
  expect_error(ss_power_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                      parameters = "not_a_label"),
               "not_a_label")

  # Scalar planning inputs.
  expect_error(ss_power_composite_sem("y2 ~ b*y1", Sigma = sigma_2, G = 5),
               "'G'")
  expect_error(ss_power_composite_sem("y2 ~ b*y1", Sigma = sigma_2, N = 3),
               "'N'")
  expect_error(ss_power_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                      desired_power = 1.2),
               "'desired_power'")
  expect_error(ss_power_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                      alpha_level = 0),
               "'alpha_level'")
})


test_that("planning refuses a parameter whose population value is zero", {
  skip_if_not_installed("lavaan")
  # y1 and y2 are independent, so the labeled path is 0 in the population and
  # no sample size gives it power above the Type I error rate. The refusal
  # comes from the analytic start value, before any simulation is spent.
  sigma_null <- diag(2)
  dimnames(sigma_null) <- list(c("y1", "y2"), c("y1", "y2"))
  expect_error(ss_power_composite_sem("y2 ~ b*y1", Sigma = sigma_null,
                                      desired_power = 0.80, G = 10),
               "zero or nearly zero")
})


test_that("empirical power tracks the analytic Wald power", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # Population: y2 = 0.2 y1 + e with unit variances, so the path estimate has
  # asymptotic variance 0.96 / N and analytic two-sided Wald power
  # pnorm(delta - z) + pnorm(-delta - z) with delta = 0.2 / sqrt(0.96 / N).
  sigma_2 <- matrix(c(1, 0.2, 0.2, 1), 2, 2,
                    dimnames = list(c("y1", "y2"), c("y1", "y2")))
  res <- ss_power_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                N = 100, G = 200, seed = 113)
  delta <- 0.2 / sqrt(0.96 / 100)
  analytic <- pnorm(delta - qnorm(0.975)) + pnorm(-delta - qnorm(0.975))
  emp <- res$value[res$term == "power_b"]
  expect_lt(abs(emp - analytic), 0.08)

  # With a single parameter of interest the composite event is the marginal
  # event, so the two proportions are identical in-sample.
  expect_identical(res$value[res$term == "composite_power"], emp)

  # The population value the analysis model implies is the generating value.
  expect_equal(res$value[res$term == "population_b"], 0.2, tolerance = 1e-6)

  # Echoes and classes.
  expect_true("specified_N" %in% res$term)
  expect_false("desired_power" %in% res$term)
  expect_s3_class(res, "dmar_ss_power")
  expect_s3_class(res, "dmar_tbl")
  expect_type(res$value, "double")
})


test_that("the composite power is at most the smallest marginal power", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  pop_3 <- "y1 ~~ 1*y1
            y2 ~~ 1*y2
            y1 ~~ 0.3*y2
            y3 ~ 0.25*y1 + 0.2*y2
            y3 ~~ 0.8*y3"
  res <- ss_power_composite_sem("y3 ~ b1*y1 + b2*y2", pop_model = pop_3,
                                N = 150, G = 100, seed = 113)
  comp <- res$value[res$term == "composite_power"]
  marg <- res$value[res$term %in% c("power_b1", "power_b2")]
  # The joint event is contained in each marginal event, exactly, in-sample.
  expect_true(all(comp <= marg))
})


test_that("Sigma and pop_model state the same population", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  pop_3 <- "y1 ~~ 1*y1
            y2 ~~ 1*y2
            y1 ~~ 0.3*y2
            y3 ~ 0.25*y1 + 0.2*y2
            y3 ~~ 0.8*y3"
  Sigma_3 <- cov_sem(pop_3)$sigma_theta
  res_sigma <- ss_power_composite_sem("y3 ~ b1*y1 + b2*y2", Sigma = Sigma_3,
                                      N = 120, G = 100, seed = 7)
  res_pop <- ss_power_composite_sem("y3 ~ b1*y1 + b2*y2", pop_model = pop_3,
                                    N = 120, G = 100, seed = 7)
  expect_equal(res_sigma$term, res_pop$term)
  expect_equal(res_sigma$value, res_pop$value)
})


test_that("a supplied seed reproduces the result and restores the RNG", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  sigma_2 <- matrix(c(1, 0.5, 0.5, 1), 2, 2,
                    dimnames = list(c("y1", "y2"), c("y1", "y2")))
  set.seed(999)
  state_before <- .Random.seed
  res_1 <- ss_power_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                  N = 60, G = 50, seed = 113)
  expect_identical(.Random.seed, state_before)
  res_2 <- ss_power_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                  N = 60, G = 50, seed = 113)
  expect_equal(res_1$value, res_2$value)
})


test_that("planning finds a sample size that meets the target", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # A generous path keeps the necessary N small: the analytic Wald size for
  # power 0.80 on b = 0.5 (residual variance 0.75) is about 24.
  sigma_2 <- matrix(c(1, 0.5, 0.5, 1), 2, 2,
                    dimnames = list(c("y1", "y2"), c("y1", "y2")))
  res <- ss_power_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                desired_power = 0.80, G = 100, seed = 113)
  n_planned <- res$value[res$term == "necessary_N"]
  expect_true(n_planned >= 10 && n_planned <= 80)
  # The reported composite power is the Monte Carlo evaluation at the
  # returned N, which met the target by construction of the search.
  expect_gte(res$value[res$term == "composite_power"], 0.80)
  expect_equal(res$value[res$term == "desired_power"], 0.80)

  # The broom verbs read the planned size and the composite power.
  td <- generics::tidy(res)
  expect_equal(td$estimate, n_planned)
  expect_equal(td$power, res$value[res$term == "composite_power"])
  gl <- generics::glance(res)
  expect_equal(gl$estimate, n_planned)
  expect_true("alpha_level" %in% names(gl))
})


test_that("population means are validated and matched by name", {
  skip_if_not_installed("lavaan")
  sigma_1 <- matrix(1, dimnames = list("y1", "y1"))
  # mu with pop_model is refused; the mean structure owns the means.
  expect_error(ss_power_composite_sem("y1 ~ m*1",
                                      pop_model = "y1 ~~ 1*y1",
                                      mu = c(y1 = 0.3)),
               "do not also give 'mu'")
  expect_error(ss_power_composite_sem("y1 ~ m*1", Sigma = sigma_1,
                                      mu = c(wrong = 0.3)),
               "one entry per observed variable")
  expect_error(ss_power_composite_sem("y1 ~ m*1", Sigma = sigma_1,
                                      mu = c(0.3, 0.4)),
               "one entry per row")
})


test_that("a mean structure target tracks the analytic one-sample power", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # The model estimates one mean from one variable with unit variance, so
  # the labeled intercept has asymptotic variance 1 / N and the analytic
  # two-sided Wald power at mu = 0.3, N = 100 is
  # pnorm(3 - 1.96) + pnorm(-3 - 1.96), about 0.851.
  sigma_1 <- matrix(1, dimnames = list("y1", "y1"))
  res <- ss_power_composite_sem("y1 ~ m*1", Sigma = sigma_1,
                                mu = c(y1 = 0.3), N = 100, G = 200,
                                seed = 113)
  analytic <- pnorm(3 - qnorm(0.975)) + pnorm(-3 - qnorm(0.975))
  expect_lt(abs(res$value[res$term == "power_m"] - analytic), 0.08)
  expect_equal(res$value[res$term == "population_m"], 0.3, tolerance = 1e-6)
})


test_that("a latent growth curve composite runs off a pop_model with means", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  pop_lgm <- "
    i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
    s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
    i ~~ 1*i
    s ~~ 0.2*s
    i ~~ -0.15*s
    t1 ~~ 0.5*t1; t2 ~~ 0.5*t2; t3 ~~ 0.5*t3; t4 ~~ 0.5*t4
    t1 ~ 0*1; t2 ~ 0*1; t3 ~ 0*1; t4 ~ 0*1
    i ~ 5*1
    s ~ 0.3*1
  "
  analysis_lgm <- "
    i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
    s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
    i ~~ cov_is*s
    t1 ~ 0*1; t2 ~ 0*1; t3 ~ 0*1; t4 ~ 0*1
    i ~ 1
    s ~ mu_s*1
  "
  res <- ss_power_composite_sem(model = analysis_lgm, pop_model = pop_lgm,
                                parameters = c("mu_s", "cov_is"),
                                N = 150, G = 50, seed = 113)
  # The analysis model recovers the generating values of both targets.
  expect_equal(res$value[res$term == "population_mu_s"], 0.3,
               tolerance = 1e-4)
  expect_equal(res$value[res$term == "population_cov_is"], -0.15,
               tolerance = 1e-4)
  comp <- res$value[res$term == "composite_power"]
  marg <- res$value[res$term %in% c("power_mu_s", "power_cov_is")]
  expect_true(all(comp <= marg))
})

test_that("a target finer than the Monte Carlo can resolve is refused", {
  skip_if_not_installed("lavaan")
  sem_Sigma <- matrix(c(1, 0.30, 0.35,
                        0.30, 1, 0.30,
                        0.35, 0.30, 1), 3, 3,
                      dimnames = list(c("y1", "y2", "y3"),
                                      c("y1", "y2", "y3")))
  # The composite power is a proportion of G replications, so it lives on
  # the grid 0, 1/G, ..., 1. A target above 1 - 1/G is met the moment every
  # replication happens to succeed, which the search reaches well below the
  # sample size the target requires (at G = 10 the planner returned N = 185
  # for desired_power = 0.999, where the true power is 0.938 and the oracle
  # size is 370), and the table then reported power 1.000 with a simulation
  # standard error of 0.
  expect_error(
    ss_power_composite_sem(model = "y3 ~ b1*y1 + b2*y2", Sigma = sem_Sigma,
                           parameters = c("b1", "b2"),
                           desired_power = 0.999, G = 200),
    "can resolve", fixed = TRUE
  )
  # The boundary is attainability, not roundness: 0.999 is exactly 999/1000,
  # so it is allowed at the default G, while 0.9995 is not.
  expect_silent(.composite_sem_check_resolution(0.999, 1000, "desired_power"))
  expect_error(.composite_sem_check_resolution(0.9995, 1000, "desired_power"),
               "Raise 'G' to at least 2001", fixed = TRUE)
  # The same ceiling governs the accuracy planner's assurance.
  expect_error(
    ss_aipe_composite_sem(model = "y3 ~ b1*y1 + b2*y2", Sigma = sem_Sigma,
                          parameters = c("b1", "b2"), width = 0.3,
                          assurance = 0.99, G = 50),
    "can resolve", fixed = TRUE
  )
})
