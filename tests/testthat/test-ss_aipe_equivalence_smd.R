# The oracles below are coded here, apart from the package: the width of
# d +/- qt(1 - alpha, 2n - 2) sqrt(2/n + d^2/(4n)) by algebra, and the
# probability that it is no wider than the target from R's pt() (AS 243)
# or, beyond pt()'s accurate range, from an integral over the normal
# variate of the noncentral t.
eq_smd_width_oracle <- function(d, n, a)
  2 * stats::qt(1 - a, 2 * n - 2) * sqrt(2 / n + d^2 / (4 * n))
eq_smd_analysis_width <- function(d, n, a) {
  f <- equivalence_smd(smd = d, n_1 = n, n_2 = n, smd_upper = 50,
                       alpha_level = a)
  f$value[f$term == "ci_upper"] - f$value[f$term == "ci_lower"]
}
# P(|T| <= k) is the same for noncentrality lambda and -lambda; pt() is
# evaluated at |lambda|, where its series converges without the precision
# warning it gives for some negative noncentrality values.
eq_smd_p_pt <- function(n, d, w, a) {
  c2 <- n * w^2 / stats::qt(1 - a, 2 * n - 2)^2 - 8
  if (c2 <= 0) return(0)
  k <- sqrt(c2 * n / 2)
  stats::pt(k, 2 * n - 2, abs(d) * sqrt(n / 2)) -
    stats::pt(-k, 2 * n - 2, abs(d) * sqrt(n / 2))
}
eq_smd_p_z <- function(n, d, w, a) {
  c2 <- n * w^2 / stats::qt(1 - a, 2 * n - 2)^2 - 8
  if (c2 <= 0) return(0)
  k <- sqrt(c2 * n / 2); lam <- abs(d) * sqrt(n / 2); nu <- 2 * n - 2
  f <- function(z) stats::pchisq(nu * (z + lam)^2 / k^2, nu,
                                 lower.tail = FALSE) * stats::dnorm(z)
  s_hi <- sqrt(stats::qchisq(1e-15, nu, lower.tail = FALSE) / nu)
  brk <- c(-lam - k * s_hi, -lam - k, -lam, -lam + k, -lam + k * s_hi)
  sum(vapply(1:4, function(i) stats::integrate(
    f, brk[i], brk[i + 1], rel.tol = 1e-12, abs.tol = 1e-15)$value,
    numeric(1)))
}

test_that("the assurance step returns the smallest n whose exact probability meets the request", {
  cases <- list(c(0.5, 0.5, 0.05, 0.80), c(0, 0.5, 0.05, 0.80),
                c(0, 0.2, 0.05, 0.99), c(1, 0.2, 0.05, 0.90),
                c(-0.8, 0.9, 0.025, 0.95), c(1.5, 2, 0.05, 0.95))
  for (cs in cases) {
    info <- paste(cs, collapse = " ")
    p0 <- ss_aipe_equivalence_smd(true_smd = cs[1], width = cs[2],
                                  alpha_level = cs[3])
    p <- ss_aipe_equivalence_smd(true_smd = cs[1], width = cs[2],
                                 alpha_level = cs[3], assurance = cs[4])
    expect_identical(p$term, c("necessary_n_per_group", "total_N", "width",
                               "true_smd", "expected_width", "assurance",
                               "assurance_attained"))
    n0 <- p0$value[p0$term == "necessary_n_per_group"]
    n <- p$value[p$term == "necessary_n_per_group"]
    attained <- p$value[p$term == "assurance_attained"]
    expect_gte(n, n0)
    expect_equal(attained, eq_smd_p_pt(n, cs[1], cs[2], cs[3]),
                 tolerance = 1e-9, info = info)
    expect_gte(attained, cs[4])
    if (n > n0) expect_lt(eq_smd_p_pt(n - 1, cs[1], cs[2], cs[3]), cs[4])
    expect_lte(p$value[p$term == "expected_width"], cs[2])
  }
  # The values the help page quotes.
  p <- ss_aipe_equivalence_smd(true_smd = 0.5, width = 0.5, assurance = 0.80)
  expect_identical(p$value[p$term == "necessary_n_per_group"], 92)
  expect_equal(eq_smd_p_pt(91, 0.5, 0.5, 0.05), 0.671, tolerance = 5e-4)
  expect_equal(eq_smd_p_pt(87, 0, 0.5, 0.05), 0)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_aipe_equivalence_smd() returns documented rows", {
  res <- ss_aipe_equivalence_smd(true_smd = 0, width = 0.20)
  expect_setequal(res$term,
                  c("necessary_n_per_group", "total_N", "width",
                    "true_smd", "expected_width"))
})

test_that("ss_aipe_equivalence_smd() recommends a finite n", {
  res <- ss_aipe_equivalence_smd(true_smd = 0, width = 0.20)
  n   <- res$value[res$term == "necessary_n_per_group"]
  expect_gt(n, 2)
  expect_lt(n, 1e5)
})

test_that("ss_aipe_equivalence_smd() expected width is at most target", {
  res <- ss_aipe_equivalence_smd(true_smd = 0, width = 0.20)
  w   <- res$value[res$term == "expected_width"]
  expect_lte(w, 0.20 + 1e-8)
})

test_that("ss_aipe_equivalence_smd() smaller width => larger n", {
  n1 <- ss_aipe_equivalence_smd(true_smd = 0, width = 0.30)$value[1]
  n2 <- ss_aipe_equivalence_smd(true_smd = 0, width = 0.10)$value[1]
  expect_lt(n1, n2)
})

# ---- The plan is for the interval equivalence_smd() reports (DR16, F09) ----
test_that("the analysis interval at the planned n is no wider than the target, and at n - 1 it is wider", {
  skip_on_cran()
  grid <- expand.grid(true_smd = c(0, 0.5, -1), width = c(6, 2.5, 0.5, 0.1),
                      alpha = c(0.05, 0.01))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    p <- ss_aipe_equivalence_smd(true_smd = g$true_smd, width = g$width,
                                 alpha_level = g$alpha)
    n <- p$value[p$term == "necessary_n_per_group"]
    info <- paste(g$true_smd, g$width, g$alpha)
    aw <- eq_smd_analysis_width(g$true_smd, n, g$alpha)
    expect_lte(aw, g$width)
    expect_equal(aw, eq_smd_width_oracle(g$true_smd, n, g$alpha),
                 tolerance = 1e-13, info = info)
    expect_identical(p$value[p$term == "expected_width"], aw)
    if (n > 2) {
      expect_gt(eq_smd_analysis_width(g$true_smd, n - 1, g$alpha), g$width)
      expect_gt(eq_smd_width_oracle(g$true_smd, n - 1, g$alpha), g$width)
    } else {
      expect_identical(n, 2)
    }
    expect_equal(p$value[p$term == "total_N"], 2 * n)
    expect_equal(attr(p, "conf_level"), 1 - 2 * g$alpha)
  }
})

test_that("the ledger cases now plan for the reported interval", {
  # DMAR 1.0.x: true_smd = 0, width = 2.5 gave n = 4 and an expected width
  # of 2.326174307353 (normal quantile), where equivalence_smd() at n = 4
  # reports 2.748071906841; width = .5 gave 87 (.498784126626 against the
  # reported .501485176103).
  p <- ss_aipe_equivalence_smd(true_smd = 0, width = 2.5)
  expect_identical(p$value[p$term == "necessary_n_per_group"], 5)
  expect_equal(p$value[p$term == "expected_width"],
               eq_smd_width_oracle(0, 5, 0.05), tolerance = 1e-14)
  expect_equal(eq_smd_width_oracle(0, 4, 0.05), 2.748071906841,
               tolerance = 1e-12)
  p <- ss_aipe_equivalence_smd(true_smd = 0, width = 0.5)
  expect_identical(p$value[p$term == "necessary_n_per_group"], 88)
  expect_equal(eq_smd_width_oracle(0, 87, 0.05), 0.501485176103,
               tolerance = 1e-12)
})

test_that("the assurance probability holds beyond the range in which pt() is accurate", {
  # true_smd = 3 with a width of .1 needs a noncentrality near 145.
  p <- ss_aipe_equivalence_smd(true_smd = 3, width = 0.1, assurance = 0.9)
  n <- p$value[p$term == "necessary_n_per_group"]
  expect_gt(3 * sqrt(n / 2), 37.62)
  expect_equal(p$value[p$term == "assurance_attained"],
               eq_smd_p_z(n, 3, 0.1, 0.05), tolerance = 1e-9)
  expect_lt(eq_smd_p_z(n - 1, 3, 0.1, 0.05), 0.9)
})

test_that("the planner draws no random numbers and takes no seed", {
  set.seed(8); before <- .Random.seed
  a <- ss_aipe_equivalence_smd(true_smd = 0.5, width = 0.4, assurance = 0.9)
  expect_identical(.Random.seed, before)
  expect_identical(ss_aipe_equivalence_smd(true_smd = 0.5, width = 0.4,
                                           assurance = 0.9), a)
  expect_error(ss_aipe_equivalence_smd(true_smd = 0.5, width = 0.4,
                                       seed = 1), "unused argument")
})

test_that("the assurance probability matches simulation from normal scores", {
  skip_on_cran()
  set.seed(2026)
  p <- ss_aipe_equivalence_smd(true_smd = 0.3, width = 0.9, assurance = 0.95)
  n <- p$value[p$term == "necessary_n_per_group"]
  attained <- p$value[p$term == "assurance_attained"]
  R <- 2000L
  x <- matrix(stats::rnorm(R * n, mean = 0.3), R)
  y <- matrix(stats::rnorm(R * n), R)
  s_p <- sqrt((apply(x, 1, stats::var) + apply(y, 1, stats::var)) / 2)
  d <- (rowMeans(x) - rowMeans(y)) / s_p
  phat <- mean(eq_smd_width_oracle(d, n, 0.05) <= 0.9)
  expect_lt(abs(phat - attained), 4 * sqrt(attained * (1 - attained) / R))
})

# ---- Planning values that are not finite, and a width left out ----
# A missing or nonfinite true_smd had stopped with R's own "missing value
# where TRUE/FALSE needed" (NA, NaN) or reached the search and was refused
# as "Required n exceeded 1e6 per group" (Inf); a width left out stopped
# with R's "argument is missing". Each is now a classed refusal that names
# the argument.
test_that("a true_smd that is not a single finite number is refused by name", {
  skip_on_cran()
  for (v in list(NA_real_, NaN, Inf, -Inf, NA, c(0.1, 0.2), "0.2")) {
    e <- expect_error(ss_aipe_equivalence_smd(true_smd = v, width = 0.5),
                      class = "dmar_argument_invalid")
    expect_identical(e$arg, "true_smd")
    expect_false(grepl("1e6", conditionMessage(e)))
  }
})

test_that("a width left out is refused by name with the classed message", {
  e <- expect_error(ss_aipe_equivalence_smd(true_smd = 0),
                    class = "dmar_argument_invalid")
  expect_identical(e$arg, "width")
  expect_match(conditionMessage(e), "none was given")
  e <- expect_error(ss_aipe_equivalence_smd(true_smd = 0, width = NA_real_),
                    class = "dmar_argument_invalid")
  expect_identical(e$arg, "width")
})

test_that("ss_aipe_equivalence_smd() gives no warning", {
  expect_no_warning(ss_aipe_equivalence_smd(true_smd = 0, width = 0.5))
})
