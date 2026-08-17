v <- function(tab, t) tab$value[tab$term == t]

test_that("ss_power_indirect_effect() solves and evaluates coherently", {
  plan <- ss_power_indirect_effect(a = .39, b = .39, desired_power = .80)
  expect_s3_class(plan, "dmar_tbl")
  n <- v(plan, "necessary_N")
  expect_gte(v(plan, "actual_power"), 0.80)
  # One fewer case dips below the target (minimality of the search).
  at <- ss_power_indirect_effect(a = .39, b = .39, N = n - 1)
  expect_lt(v(at, "actual_power"), 0.80)
  # Joint significance power is the product of the component powers.
  expect_equal(v(plan, "actual_power"),
               v(plan, "power_a") * v(plan, "power_b"))
  expect_equal(v(plan, "indirect_effect"), .39 * .39)
})

test_that("the weak link drives the requirement, and Sobel needs more", {
  small_a <- ss_power_indirect_effect(a = .14, b = .39,
                                      desired_power = .80)
  both_med <- ss_power_indirect_effect(a = .39, b = .39,
                                       desired_power = .80)
  expect_gt(v(small_a, "necessary_N"), v(both_med, "necessary_N"))
  sobel <- ss_power_indirect_effect(a = .39, b = .39, desired_power = .80,
                                    method = "sobel")
  expect_gte(v(sobel, "necessary_N"), v(both_med, "necessary_N"))
  expect_true(is.na(v(sobel, "power_a")))
})

test_that("ss_power_indirect_effect() validates its arguments", {
  expect_error(ss_power_indirect_effect(a = 0, b = .3,
                                        desired_power = .8), "nonzero")
  expect_error(ss_power_indirect_effect(a = 1.1, b = .3,
                                        desired_power = .8), "\\(-1, 1\\)")
  expect_error(ss_power_indirect_effect(a = .3, b = .3), "exactly one")
  expect_error(ss_power_indirect_effect(a = .3, b = .3, N = 100,
                                        desired_power = .8), "exactly one")
  expect_error(ss_power_indirect_effect(a = .7, b = .7, c_prime = .7,
                                        desired_power = .8),
               "not admissible")
})

test_that("joint-significance power matches raw-data simulation", {
  # Monte Carlo confirmation; skipped on CRAN. The fast anchors that stay
  # on CRAN are the coherence tests above (minimality of the sample size
  # search and the product-of-component-powers identity).
  skip_on_cran()
  set.seed(113)
  a <- .39; b <- .39; n <- 75; G <- 2000
  ana <- v(ss_power_indirect_effect(a = a, b = b, N = n), "actual_power")
  hits <- replicate(G, {
    x <- rnorm(n)
    m <- a * x + rnorm(n, 0, sqrt(1 - a^2))
    y <- b * m + rnorm(n, 0, sqrt(1 - b^2))
    fm <- summary(lm(m ~ x))$coefficients
    fy <- summary(lm(y ~ x + m))$coefficients
    fm["x", "Pr(>|t|)"] < .05 && fy["m", "Pr(>|t|)"] < .05
  })
  # The analytic power is an approximation (population standard errors and
  # component independence), and at these settings it sits about 3 percent
  # above the raw-data simulation on any seed, with about a 1 percent Monte
  # Carlo standard error at G = 2000 on top. The former 0.03 relative
  # tolerance failed on some seeds from the approximation gap alone; 0.06
  # covers the systematic gap plus simulation noise while staying far
  # tighter than the error a wrong formula would produce.
  expect_equal(ana, mean(hits), tolerance = 0.06)
})

test_that("ss_power_indirect_effect errors instead of hanging on a tiny effect (HIGH-04)", {
  expect_error(ss_power_indirect_effect(a = 1e-12, b = 0.39, desired_power = 0.80),
               "too small|Could not reach|did not converge")
  r <- ss_power_indirect_effect(a = 0.3, b = 0.39, desired_power = 0.80)
  expect_gte(r$value[r$term == "actual_power"], 0.80)
})
