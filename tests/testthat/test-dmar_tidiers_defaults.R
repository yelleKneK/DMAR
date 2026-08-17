# The default tidy()/glance() on dmar_tbl, and the AIPE family verbs. The
# defaults are the floor under the whole package: any tidy-returning DMAR
# function answers the broom verbs, and family classes override with richer
# views where the structure warrants one. This is the contract the
# dmar_output vignette states as "the rule is the same everywhere."

test_that("every dmar_tbl answers tidy() and glance() through the defaults", {
  v <- var_ete(sigma2 = 4, sigma2_Z = 1, n_1 = 20, n_2 = 20,
               beta_1 = .5, beta_2 = .5, mu_Z = 0, type = "population")
  td <- generics::tidy(v)
  expect_named(td, c("term", "estimate"))
  expect_identical(td$term, as.character(v$term))
  expect_identical(td$estimate, as.numeric(v$value))

  gl <- generics::glance(v)
  expect_equal(nrow(gl), 1L)
  expect_identical(unname(unlist(gl)), as.numeric(v$value))

  set.seed(113)
  rt <- randomization_test(group_1 = rnorm(6), group_2 = rnorm(6) + 1)
  expect_equal(nrow(generics::tidy(rt)), nrow(rt))
  expect_equal(nrow(generics::glance(rt)), 1L)
})

test_that("glance() disambiguates repeated terms instead of dropping them", {
  # A regions-of-significance table reports one block per group pair, so its
  # terms repeat; the wide view suffixes them the way make.unique does.
  r <- regions_of_significance(lm(mpg ~ wt * factor(am), mtcars))
  g <- generics::glance(r)
  expect_equal(nrow(g), 1L)
  expect_false(anyDuplicated(names(g)) > 0)
  expect_equal(ncol(g), nrow(r))
})

test_that("the AIPE planners tidy to size and width through the registries", {
  a <- ss_aipe_smd(delta = .5, width = .3)
  expect_s3_class(a, "dmar_ss_aipe")
  td <- generics::tidy(a)
  expect_named(td, c("term", "estimate", "width"))
  expect_identical(td$term, "sample_size")
  expect_equal(td$estimate,
               a$value[a$term == "necessary_n_per_group"])
  expect_equal(td$width, a$value[a$term == "width"])

  gl <- generics::glance(a)
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("estimate", "width", "supposed_smd") %in% names(gl)))

  # A planner whose width row uses the family's other name resolves too.
  i <- ss_aipe_icc(rho = .7, k = 3, width = .2)
  expect_equal(generics::tidy(i)$width,
               i$value[i$term == "width_target"])
})

test_that("family verbs still win over the defaults", {
  p <- ss_power_smd(smd = .5, desired_power = .8)
  expect_true("power" %in% names(generics::tidy(p)))
  # ss_power_pcm, the one planner that used to error on tidy(), now
  # dispatches with its family.
  q <- ss_power_pcm(beta = .3, tau = .05, level_1_variance = 1,
                    frequency = 5, duration = 4, desired_power = .8)
  expect_s3_class(q, "dmar_ss_power")
  expect_true(is.finite(generics::tidy(q)$estimate))
  # ci_eta_squared_generalized rejoins its family.
  g <- suppressMessages(
    ci_eta_squared_generalized(aov(len ~ supp * factor(dose), ToothGrowth))
  )
  expect_s3_class(g, "dmar_ci_anova")
})
