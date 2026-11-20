# The exact full width of the back-transformed Fisher's Z interval at sample
# size n, the interval ss_aipe_r() plans for and correlations_test() reports
# for method = "pearson" (SE(z) = 1 / sqrt(n - 3); Bonett & Wright, 2000).
.fisher_z_width <- function(n, rho, conf_level) {
  z    <- atanh(rho)
  crit <- qnorm(1 - (1 - conf_level) / 2)
  se   <- 1 / sqrt(n - 3)
  tanh(z + crit * se) - tanh(z - crit * se)
}

test_that("ss_aipe_r() expected_width row is the exact width at necessary_N", {
  res <- ss_aipe_r(true_cor = 0.30, width = 0.20)
  n   <- res$value[res$term == "necessary_N"]
  expect_equal(res$value[res$term == "expected_width"],
               .fisher_z_width(n, rho = 0.30, conf_level = 0.95))
  # The smallest such N: the width is met at N and missed at N - 1.
  expect_true(.fisher_z_width(n, 0.30, 0.95) <= 0.20 &&
                .fisher_z_width(n - 1, 0.30, 0.95) > 0.20)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_aipe_r() returns documented rows with the family classes", {
  res <- ss_aipe_r(true_cor = 0.30, width = 0.20)
  expect_s3_class(res, "dmar_ss_aipe")
  expect_s3_class(res, "dmar_tbl")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("necessary_N", "expected_width", "true_cor",
                    "width", "conf_level") %in% res$term))
  expect_false(any(c("rho", "true_rho") %in% res$term))
  expect_true(is.numeric(res$value))
  n_val <- res$value[res$term == "necessary_N"]
  expect_gte(n_val, 4)
  expect_equal(n_val, round(n_val))
})

test_that("ss_aipe_r() matches a brute-force search over a grid", {
  skip_on_cran()
  brute_force_n <- function(rho, width, conf_level) {
    n <- 4L
    while (.fisher_z_width(n, rho, conf_level) > width) n <- n + 1L
    n
  }
  for (rho in c(0, 0.30, -0.50, 0.70, 0.90)) {
    for (width in c(0.10, 0.20, 0.30)) {
      for (conf_level in c(0.90, 0.95, 0.99)) {
        planned <- ss_aipe_r(true_cor = rho, width = width,
                             conf_level = conf_level)
        n <- planned$value[planned$term == "necessary_N"]
        expect_equal(n, brute_force_n(rho, width, conf_level),
                     info = sprintf("rho = %.2f, width = %.2f, conf = %.2f",
                                    rho, width, conf_level))
        # Minimality restated directly: the width at n meets the target
        # and the width at n - 1 does not.
        expect_lte(.fisher_z_width(n, rho, conf_level), width)
        if (n > 4) expect_gt(.fisher_z_width(n - 1, rho, conf_level), width)
      }
    }
  }
})

test_that("ss_aipe_r() planning value closer to zero is more conservative", {
  n_0  <- ss_aipe_r(true_cor = 0.00, width = 0.20)$value[1]
  n_5  <- ss_aipe_r(true_cor = 0.50, width = 0.20)$value[1]
  n_9  <- ss_aipe_r(true_cor = 0.90, width = 0.20)$value[1]
  expect_gte(n_0, n_5)
  expect_gte(n_5, n_9)
})

test_that("ss_aipe_r() is symmetric in the sign of rho", {
  expect_equal(ss_aipe_r(true_cor = 0.40, width = 0.20)$value[1],
               ss_aipe_r(true_cor = -0.40, width = 0.20)$value[1])
})

test_that("ss_aipe_r() assurance inflates the sample size", {
  res_50 <- ss_aipe_r(true_cor = 0.30, width = 0.20)
  res_80 <- ss_aipe_r(true_cor = 0.30, width = 0.20, assurance = 0.80)
  res_99 <- ss_aipe_r(true_cor = 0.30, width = 0.20, assurance = 0.99)
  expect_gt(res_80$value[1], res_50$value[1])
  expect_gt(res_99$value[1], res_80$value[1])
  expect_identical(res_80$term,
                   c("necessary_N", "expected_width", "true_cor", "width",
                     "conf_level", "assurance", "assurance_attained"))
  expect_false("assurance_attained" %in% res_50$term)
})

# ---- The assurance step against the exact law of r -------------------------
#
# An oracle coded apart from the package: the density of the sample
# correlation r of n bivariate normal observations in Hotelling's form,
#   f(r) = (n - 2) Gamma(n - 1) / (sqrt(2 pi) Gamma(n - 1/2))
#          (1 - rho^2)^((n - 1)/2) (1 - r^2)^((n - 4)/2)
#          (1 - rho r)^(-(n - 3/2)) 2F1(1/2, 1/2; n - 1/2; (1 + rho r)/2),
# with the hypergeometric series summed here, integrated on Fisher's Z
# scale. The package computes the distribution function from a central t
# integrated against a beta density (.partial_r_cdf()); nothing is shared.
.oracle_h2f1 <- function(x, n) {
  term <- 1
  total <- 1
  k <- 0
  repeat {
    term <- term * (k + 0.5)^2 / ((n - 0.5 + k) * (k + 1)) * x
    total <- total + term
    k <- k + 1
    if (abs(term) < 1e-17 * total || k > 1e5) break
  }
  total
}
.oracle_r_density <- function(r, rho, n) {
  log_c <- log(n - 2) + lgamma(n - 1) - 0.5 * log(2 * pi) - lgamma(n - 0.5)
  vapply(r, function(x) {
    exp(log_c + (n - 1) / 2 * log1p(-rho^2) + (n - 4) / 2 * log1p(-x^2) -
          (n - 1.5) * log1p(-rho * x)) * .oracle_h2f1((1 + rho * x) / 2, n)
  }, numeric(1))
}
# P(-c <= r <= c), integrated over Fisher's Z in pieces that bracket the
# peak of the density (standard deviation about 1 / sqrt(n - 3)).
.oracle_r_between <- function(c, rho, n) {
  f <- function(z) .oracle_r_density(tanh(z), rho, n) / cosh(z)^2
  cuts <- atanh(rho) + c(-12, -6, -3, -1, 0, 1, 3, 6, 12) / sqrt(n - 3)
  edges <- sort(unique(c(-atanh(c), cuts[abs(cuts) < atanh(c)], atanh(c))))
  sum(vapply(seq_len(length(edges) - 1L), function(i)
    stats::integrate(f, edges[i], edges[i + 1L], rel.tol = 1e-12,
                     subdivisions = 2000L)$value, numeric(1)))
}
# The probability that the realized interval is no wider than `width` at n.
.oracle_width_probability <- function(n, rho, width, conf_level = 0.95) {
  if (.fisher_z_width(n, 0, conf_level) <= width) return(1)
  r_star <- stats::uniroot(
    function(r) .fisher_z_width(n, r, conf_level) - width,
    c(0, 1 - 1e-13), tol = 1e-14)$root
  1 - .oracle_r_between(r_star, rho, n)
}

test_that("the oracle density of r integrates to 1 and gives the t law at rho = 0", {
  expect_equal(.oracle_r_between(1 - 1e-12, 0.6, 12), 1, tolerance = 1e-10)
  expect_equal(.oracle_r_between(1 - 1e-14, -0.9, 72), 1, tolerance = 1e-10)
  # At rho = 0, sqrt(n - 2) r / sqrt(1 - r^2) is Student's t on n - 2 df.
  c <- 0.3
  expect_equal(.oracle_r_between(c, 0, 20),
               1 - 2 * stats::pt(-sqrt(18) * c / sqrt(1 - c^2), 18),
               tolerance = 1e-10)
})

test_that("ss_aipe_r() with assurance attains the request under the exact law of r", {
  skip_on_cran()
  # The 1.0.x step multiplied the plan by qchisq(assurance, N - 2) / (N - 2):
  # at true_cor = .90 and width .10 it returned 72 for .80 and 92 for .99,
  # where the oracle gives .672 and .887.
  expect_equal(.oracle_width_probability(72, 0.90, 0.10), 0.672086002,
               tolerance = 1e-8)
  expect_equal(.oracle_width_probability(92, 0.90, 0.10), 0.887067986,
               tolerance = 1e-8)
  for (rho in c(0.90, 0.50, 0.10, 0, -0.50)) {
    for (width in c(0.10, 0.20)) {
      for (a in c(0.80, 0.99)) {
        info <- sprintf("rho = %.2f, width = %.2f, assurance = %.2f",
                        rho, width, a)
        res <- ss_aipe_r(true_cor = rho, width = width, assurance = a)
        n <- res$value[res$term == "necessary_N"]
        n_without <- ss_aipe_r(true_cor = rho, width = width)$value[1]
        attained <- res$value[res$term == "assurance_attained"]
        oracle <- .oracle_width_probability(n, rho, width)
        expect_equal(attained, oracle, tolerance = 1e-8, info = info)
        expect_gte(oracle, a)
        expect_gte(n, n_without)
        # Minimality: N - 1 falls short whenever N is above the plan
        # without assurance.
        if (n > n_without)
          expect_lt(.oracle_width_probability(n - 1, rho, width), a,
                    label = info)
      }
    }
  }
  # The reproduction cases of the review.
  res_80 <- ss_aipe_r(true_cor = 0.90, width = 0.10, assurance = 0.80)
  res_99 <- ss_aipe_r(true_cor = 0.90, width = 0.10, assurance = 0.99)
  expect_identical(res_80$value[1], 83)
  expect_identical(res_99$value[1], 122)
})

test_that("ss_aipe_r() with assurance is symmetric in the sign of rho and needs nothing extra at 0", {
  for (rho in c(0.10, 0.50, 0.90, 0.99)) {
    pos <- ss_aipe_r(true_cor = rho, width = 0.10, assurance = 0.90)
    neg <- ss_aipe_r(true_cor = -rho, width = 0.10, assurance = 0.90)
    expect_identical(neg$value[-3], pos$value[-3])
    expect_identical(neg$value[3], -rho)
  }
  zero <- ss_aipe_r(true_cor = 0, width = 0.10, assurance = 0.99)
  expect_identical(zero$value[1], ss_aipe_r(true_cor = 0, width = 0.10)$value[1])
  expect_identical(zero$value[zero$term == "assurance_attained"], 1)
})

test_that("the width probability of ss_aipe_r() does not decrease above the plan without assurance", {
  skip_on_cran()
  for (rho in c(0.90, -0.50)) {
    n_without <- ss_aipe_r(true_cor = rho, width = 0.20)$value[1]
    n_99 <- ss_aipe_r(true_cor = rho, width = 0.20, assurance = 0.99)$value[1]
    p <- vapply(n_without:(n_99 + 10), .oracle_width_probability, numeric(1),
                rho = rho, width = 0.20)
    expect_true(all(diff(p) >= -1e-12), label = paste("rho", rho))
  }
})

test_that("the attained assurance of ss_aipe_r() agrees with simulated bivariate normal data", {
  skip_on_cran()
  # Correlations of n pairs drawn from a bivariate normal population; the
  # interval is the Fisher's Z interval of correlations_test(). 4,000
  # samples; the standard error of a proportion near .8 is about .0063.
  check <- function(rho, width, a, G = 4000) {
    res <- ss_aipe_r(true_cor = rho, width = width, assurance = a)
    n <- res$value[1]
    attained <- res$value[res$term == "assurance_attained"]
    x <- matrix(stats::rnorm(n * G), n)
    y <- rho * x + sqrt(1 - rho^2) * matrix(stats::rnorm(n * G), n)
    r <- vapply(seq_len(G), function(g) stats::cor(x[, g], y[, g]), 0)
    share <- mean(.fisher_z_width(n, r, 0.95) <= width)
    expect_lt(abs(share - attained), 4 * sqrt(attained * (1 - attained) / G))
  }
  set.seed(20261111)
  check(0.90, 0.10, 0.80)
  check(-0.50, 0.20, 0.90)
})

test_that("ss_aipe_r() rejects bad inputs", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  expect_error(ss_aipe_r(true_cor = 1.5, width = 0.2), "in \\(-1, 1\\)")
  expect_error(ss_aipe_r(true_cor = c(0.2, 0.3), width = 0.2), "in \\(-1, 1\\)")
  expect_error(ss_aipe_r(true_cor = 0.3, width = 0), "positive")
  expect_error(ss_aipe_r(true_cor = 0.3, width = -0.1), "positive")
  expect_error(ss_aipe_r(true_cor = 0.3, width = 0.2, conf_level = 1.5),
               "in \\(0, 1\\)")
  expect_error(ss_aipe_r(true_cor = 0.3, width = 0.2, assurance = 0.4),
               "in \\(0.5, 1\\)")
  expect_error(ss_aipe_r(true_cor = 0.3, width = 0.2, assurance = 1),
               "in \\(0.5, 1\\)")
})

test_that("ss_aipe_r() tidies through the dmar_ss_aipe registry", {
  res <- ss_aipe_r(true_cor = 0.30, width = 0.20)
  td  <- generics::tidy(res)
  expect_equal(td$estimate, res$value[res$term == "necessary_N"])
  expect_equal(td$width, 0.20)
  gl <- generics::glance(res)
  expect_equal(gl$estimate, td$estimate)
  expect_true("true_cor" %in% names(gl))
})

test_that("ss_aipe_r() gives no warning, with or without assurance", {
  expect_no_warning(ss_aipe_r(true_cor = .9, width = .1, assurance = .8))
  expect_no_warning(ss_aipe_r(true_cor = .9, width = .1))
  expect_no_warning(ss_aipe_r(true_cor = .99, width = .05, assurance = .8))
})
