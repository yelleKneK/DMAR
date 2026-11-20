test_that("gwet_ac() returns documented rows", {
  set.seed(113)
  r1 <- sample(c("A","B","C"), 30, replace = TRUE)
  r2 <- ifelse(runif(30) < 0.7, r1, sample(c("A","B","C"), 30, TRUE))
  res <- gwet_ac(cbind(r1, r2))
  expect_setequal(res$term,
                  c("gwet_ac", "se", "ci_lower", "ci_upper",
                    "percent_agreement", "chance_agreement", "n_units"))
})

test_that("gwet_ac() is 1.0 on perfect agreement", {
  # 3 raters give identical scores on every unit.
  R <- matrix(rep(1:4, times = 3), 4, 3)
  res <- gwet_ac(R)
  expect_equal(res$value[res$term == "gwet_ac"], 1, tolerance = 1e-10)
})

test_that("gwet_ac() rejects too-few raters / categories", {
  expect_error(gwet_ac(matrix(1:5, 5, 1)), "at least 2 raters")
  expect_error(gwet_ac(matrix(rep(1, 10), 5, 2)), "at least 2 categories")
})

test_that("gwet_ac() linear weights yield AC2 <= AC1 on noisy ordinal", {
  set.seed(113)
  r1 <- sample(1:5, 50, replace = TRUE)
  r2 <- pmin(5, pmax(1, r1 + sample(-1:1, 50, replace = TRUE)))
  ac1 <- gwet_ac(cbind(r1, r2))$value[1]
  ac2 <- gwet_ac(cbind(r1, r2), weights = "linear")$value[1]
  expect_true(ac2 >= ac1 - 1e-6)  # linear weights are more forgiving
})

test_that("gwet_ac() is invariant to appended all-missing units (HIGH-07)", {
  # Information invariance oracle: a unit with no ratings carries no agreement
  # or chance information, so appending any number of all-NA rows must leave
  # every numeric result identical. An all-missing row used to divide the
  # category proportions by zero raters and turn AC and its SE into NaN.
  base <- rbind(c(1, 1, 1), c(2, 2, 3), c(1, 2, 1), c(3, 3, 3), c(2, 2, 2))
  ref  <- gwet_ac(base)
  one  <- gwet_ac(rbind(base, c(NA, NA, NA)))
  two  <- gwet_ac(rbind(base, c(NA, NA, NA), c(NA, NA, NA)))

  expect_false(any(is.nan(one$value)))
  # AC1 and SE specifically, then the whole numeric vector.
  expect_equal(one$value[one$term == "gwet_ac"], ref$value[ref$term == "gwet_ac"])
  expect_equal(one$value[one$term == "se"],      ref$value[ref$term == "se"])
  expect_equal(one$value, ref$value)
  expect_equal(two$value, ref$value)
})

test_that("gwet_ac() stays near the observed agreement where kappa collapses (help claim)", {
  # The help page says that when one category is very common, Cohen's kappa
  # can be near zero although the raters agree on almost every unit, and
  # that the chance-agreement term of AC1 shrinks in that case. Oracle: the
  # closed forms worked by hand for a 2 x 2 table with 90 (A, A), 5 (A, B),
  # 4 (B, A), and 1 (B, B) of 100 units. Observed agreement p_a = .91.
  # AC1: pi_A = (95 + 94) / 200 = .945, p_e = 2 * .945 * .055 = .10395,
  #      AC1 = (.91 - .10395) / (1 - .10395).
  # Kappa: p_e = .95 * .94 + .05 * .06 = .896, kappa = (.91 - .896) / .104.
  r1 <- c(rep("A", 95), rep("B", 5))
  r2 <- c(rep("A", 90), rep("B", 5), rep("A", 4), "B")
  ac <- gwet_ac(cbind(r1, r2))
  kp <- cohen_kappa(r1, r2)
  pe_ac1 <- 2 * 0.945 * 0.055
  expect_equal(ac$value[ac$term == "percent_agreement"], 0.91, tolerance = 1e-12)
  expect_equal(ac$value[ac$term == "chance_agreement"], pe_ac1, tolerance = 1e-12)
  expect_equal(ac$value[ac$term == "gwet_ac"], (0.91 - pe_ac1) / (1 - pe_ac1),
               tolerance = 1e-12)
  expect_equal(kp$kappa, (0.91 - 0.896) / (1 - 0.896), tolerance = 1e-12)
  expect_lt(kp$kappa, 0.2)
  expect_gt(ac$value[ac$term == "gwet_ac"], 0.85)
})
