# Tests for ss_aipe_composite_sem(). Validation errors run everywhere; the
# Monte Carlo paths are off CRAN for runtime. The observed-variable path
# models have closed-form population covariance matrices, so the simulated
# interval widths can be checked against the asymptotic widths they estimate.

test_that("ss_aipe_composite_sem() validates its inputs", {
  skip_if_not_installed("lavaan")
  sigma_2 <- matrix(c(1, 0.5, 0.5, 1), 2, 2,
                    dimnames = list(c("y1", "y2"), c("y1", "y2")))

  expect_error(ss_aipe_composite_sem("y2 ~ b*y1", desired_width = 0.3),
               "exactly one")
  expect_error(ss_aipe_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                     desired_width = -0.3),
               "'desired_width'")
  # Several widths must be named so none can attach to the wrong parameter.
  expect_error(ss_aipe_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                     desired_width = c(0.3, 0.4)),
               "named")
  expect_error(ss_aipe_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                     desired_width = c(wrong = 0.3)),
               "one entry per parameter")
  expect_error(ss_aipe_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                     desired_width = 0.3, assurance = 0.3),
               "'assurance'")
  expect_error(ss_aipe_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                     desired_width = 0.3, conf_level = 1.2),
               "'conf_level'")
  expect_error(ss_aipe_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                     desired_width = 0.3, G = 5),
               "'G'")
})


test_that("simulated widths track the asymptotic Wald width", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # Population: y2 = 0.2 y1 + e with unit variances, so the path estimate has
  # asymptotic variance 0.96 / N and asymptotic 95% interval width
  # 2 * 1.96 * sqrt(0.96 / 100) at N = 100, about 0.384.
  sigma_2 <- matrix(c(1, 0.2, 0.2, 1), 2, 2,
                    dimnames = list(c("y1", "y2"), c("y1", "y2")))
  res <- ss_aipe_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                               desired_width = 0.40, N = 100, G = 300,
                               seed = 113)
  asymptotic_width <- 2 * qnorm(0.975) * sqrt(0.96 / 100)
  expect_lt(abs(res$value[res$term == "mean_width_b"] - asymptotic_width),
            0.03)

  # With one parameter of interest the joint proportion is the marginal one.
  expect_identical(res$value[res$term == "composite_assurance"],
                   res$value[res$term == "width_within_desired_b"])

  # Echoes and classes.
  expect_true("specified_N" %in% res$term)
  expect_false("assurance" %in% res$term)
  expect_equal(res$value[res$term == "desired_width_b"], 0.40)
  expect_equal(res$value[res$term == "conf_level"], 0.95)
  expect_s3_class(res, "dmar_tbl")
  expect_type(res$value, "double")
})


test_that("the joint proportion is at most each marginal proportion", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  pop_3 <- "y1 ~~ 1*y1
            y2 ~~ 1*y2
            y1 ~~ 0.3*y2
            y3 ~ 0.25*y1 + 0.2*y2
            y3 ~~ 0.8*y3"
  res <- ss_aipe_composite_sem("y3 ~ b1*y1 + b2*y2", pop_model = pop_3,
                               desired_width = 0.40, N = 120, G = 100,
                               seed = 113)
  joint <- res$value[res$term == "composite_assurance"]
  marg <- res$value[res$term %in% c("width_within_desired_b1",
                                    "width_within_desired_b2")]
  expect_true(all(joint <= marg))
})


test_that("a named desired_width vector matches widths to labels", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  pop_3 <- "y1 ~~ 1*y1
            y2 ~~ 1*y2
            y1 ~~ 0.3*y2
            y3 ~ 0.25*y1 + 0.2*y2
            y3 ~~ 0.8*y3"
  # Names in scrambled order still attach to the right parameters.
  res <- ss_aipe_composite_sem("y3 ~ b1*y1 + b2*y2", pop_model = pop_3,
                               desired_width = c(b2 = 0.50, b1 = 0.35),
                               N = 120, G = 50, seed = 113)
  expect_equal(res$value[res$term == "desired_width_b1"], 0.35)
  expect_equal(res$value[res$term == "desired_width_b2"], 0.50)
})


test_that("expected-width and assurance planning agree with the closed form and each other", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # One expected-width search serves both halves of this test; an identical
  # second search used to run in a separate test for the assurance
  # comparison alone.
  #
  # The closed-form no-assurance size for width 0.35 on a path with residual
  # variance 0.64 is 4 * 1.96^2 * 0.64 / 0.35^2, about 81; the Monte Carlo
  # criterion (mean simulated width within the desired width) should resolve
  # to its neighborhood.
  sigma_2 <- matrix(c(1, 0.6, 0.6, 1), 2, 2,
                    dimnames = list(c("y1", "y2"), c("y1", "y2")))
  res_mean <- ss_aipe_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                                    desired_width = 0.35, G = 100,
                                    seed = 113)
  n_planned <- res_mean$value[res_mean$term == "necessary_N"]
  expect_true(n_planned >= 55 && n_planned <= 115)
  # The mean width at the returned N met the criterion by construction.
  expect_lte(res_mean$value[res_mean$term == "mean_width_b"], 0.35)
  expect_false("assurance" %in% res_mean$term)

  res <- ss_aipe_composite_sem("y2 ~ b*y1", Sigma = sigma_2,
                               desired_width = 0.35, assurance = 0.80,
                               G = 100, seed = 113)
  expect_gte(res$value[res$term == "composite_assurance"], 0.80)
  expect_equal(res$value[res$term == "assurance"], 0.80)

  # Guaranteeing the joint event with 80% assurance cannot need a smaller
  # sample size than meeting the expected widths, up to Monte Carlo wobble.
  expect_gte(res$value[res$term == "necessary_N"], n_planned - 8)
})
