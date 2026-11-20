## ss_aipe_mean(): the AIPE planner for a single mean (one sample, or the
## mean of paired differences). The anchors are closed form: the expected
## width under normality is 2 t c4(N) sigma / sqrt(N) with c4(N) = E(S) /
## sigma, and the probability that the realized width is at most the target
## is a chi square probability. Each test recomputes its quantity here, from
## its definition, rather than through the function's own helpers.

# E(S) / sigma for a normal sample of size n, from the gamma function ratio
# directly (accurate to about 1e-12 for the moderate n used below).
c4_gamma <- function(n) sqrt(2 / (n - 1)) * exp(lgamma(n / 2) - lgamma((n - 1) / 2))

expected_width <- function(n, sd, conf_level) {
  2 * stats::qt(1 - (1 - conf_level) / 2, n - 1) * c4_gamma(n) * sd / sqrt(n)
}

# P(W <= width) under normality: W <= width exactly when
# (n - 1) S^2 / sigma^2 <= (n - 1) n width^2 / (4 t^2 sigma^2).
prob_narrow <- function(n, sd, width, conf_level) {
  t <- stats::qt(1 - (1 - conf_level) / 2, n - 1)
  stats::pchisq((n - 1) * n * width^2 / (4 * t^2 * sd^2), df = n - 1)
}

# The same probability under the scaled chi square of the kurtosis
# adjustment: S^2 / sigma^2 is taken as chi square on nu df over nu, with nu
# matched to Var(S^2) = sigma^4 {2 / (n - 1) + kappa / n}.
prob_narrow_k <- function(n, sd, width, conf_level, kappa) {
  t <- stats::qt(1 - (1 - conf_level) / 2, n - 1)
  nu <- 2 / (2 / (n - 1) + kappa / n)
  stats::pchisq(nu * n * width^2 / (4 * t^2 * sd^2), df = nu)
}

# The bias-corrected sample excess kurtosis, written from the moments with
# divisor n, G2 = {(n + 1) g2 + 6} (n - 1) / {(n - 2)(n - 3)} with
# g2 = m4 / m2^2 - 3, rather than through DMAR's kurtosis().
excess_kurtosis_G2 <- function(x) {
  n <- length(x)
  d <- x - mean(x)
  g2 <- mean(d^4) / mean(d^2)^2 - 3
  ((n + 1) * g2 + 6) * (n - 1) / ((n - 2) * (n - 3))
}

# The definition of necessary_N by brute force: every N from 2 to an n_max
# is evaluated, and the answer is one more than the largest N at which the
# criterion fails (2 when it fails nowhere). n_max is doubled until the
# criterion is certain to hold beyond it: for the expected width because
# E(W) falls in N, and with an assurance by Cantelli's inequality, under
# which qchisq(g, nu) / nu is at most 1 + sqrt(2 g / ((1 - g) nu)), a bound
# that falls as nu rises.
scan_N <- function(sd, w, cl, g = NULL, k = 0) {
  q <- 1 - (1 - cl) / 2
  ew <- function(n) 2 * stats::qt(q, n - 1) * c4_gamma(n) * sd / sqrt(n)
  bnd <- function(n) n * (w / 2)^2 / (stats::qt(q, n - 1)^2 * sd^2)
  nu_of <- function(n) 2 / (2 / (n - 1) + k / n)
  beyond <- if (is.null(g)) {
    function(n) ew(n) <= w
  } else {
    function(n) bnd(n) >= 1 + sqrt(2 * g / ((1 - g) * nu_of(n)))
  }
  n_max <- 512
  while (!beyond(n_max)) n_max <- 2 * n_max
  n <- 2:n_max
  ok <- if (is.null(g)) {
    ew(n) <= w
  } else {
    stats::pchisq(nu_of(n) * bnd(n), nu_of(n)) >= g
  }
  fails <- which(!ok)
  if (length(fails)) n[max(fails)] + 1 else 2
}

planned_N <- function(...) {
  r <- ss_aipe_mean(...)
  r$value[r$term == "necessary_N"]
}

test_that("the expected width plan matches the S = sigma plan of ss_aipe_c() or is one smaller", {
  # ss_aipe_c() sets S equal to sigma in the width. Applied to a single
  # mean that is the smallest N with 2 t sigma / sqrt(N) <= width; since
  # E(S) < sigma, the exact expected width plan needs the same N or one less.
  sigma_plan <- function(w, cl) {
    q <- 1 - (1 - cl) / 2
    n <- max(2, floor((stats::qnorm(q) / (w / 2))^2))
    while (2 * stats::qt(q, n - 1) / sqrt(n) > w) n <- n + 1
    n
  }
  d <- integer(0)
  for (cl in c(0.80, 0.90, 0.95, 0.99, 0.999)) {
    for (w in exp(seq(log(0.01), log(5), length.out = 25))) {
      d <- c(d, sigma_plan(w, cl) - planned_N(true_sd = 1, width = w, conf_level = cl))
    }
  }
  expect_true(all(d %in% c(0L, 1L)))
  expect_true(any(d == 0L))
  expect_true(any(d == 1L))
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

test_that("ss_aipe_mean() returns the planner table with the answer row necessary_N", {
  r <- ss_aipe_mean(true_sd = 1, width = 0.5)
  expect_s3_class(r, "dmar_ss_aipe")
  expect_s3_class(r, "dmar_tbl")
  expect_s3_class(r, "data.frame")
  expect_named(r, c("term", "value"))
  expect_type(r$value, "double")
  expect_identical(r$term, c("necessary_N", "expected_width", "true_sd", "width",
                             "conf_level"))
  expect_equal(attr(r, "conf_level"), 0.95)

  v <- stats::setNames(r$value, r$term)
  expect_equal(unname(v["true_sd"]), 1)
  expect_equal(unname(v["width"]), 0.5)
  expect_equal(unname(v["expected_width"]),
               expected_width(v[["necessary_N"]], 1, 0.95), tolerance = 1e-10)

  # The assurance and the kurtosis it was planned with appear only with an
  # assurance, the planner convention.
  ra <- ss_aipe_mean(true_sd = 1, width = 0.5, assurance = 0.80, true_kurtosis = -1)
  expect_identical(ra$term, c("necessary_N", "expected_width", "true_sd", "width",
                              "conf_level", "assurance", "true_kurtosis"))
  expect_equal(ra$value[ra$term == "assurance"], 0.80)
  expect_equal(ra$value[ra$term == "true_kurtosis"], -1)

  # The broom verbs read the planned size and the target width.
  td <- generics::tidy(ra)
  expect_equal(td$estimate, ra$value[ra$term == "necessary_N"])
  expect_equal(td$width, 0.5)
  gl <- generics::glance(ra)
  expect_equal(gl$estimate, td$estimate)
  expect_equal(gl$assurance, 0.80)
})

test_that("the expected width at necessary_N is at most the target and at necessary_N - 1 above it", {
  skip_on_cran()
  grid <- expand.grid(sd = c(1, 13.29708, 0.09587),
                      ratio = c(0.1, 0.3, 1.5),
                      conf_level = c(0.90, 0.95, 0.99))
  for (i in seq_len(nrow(grid))) {
    s  <- grid$sd[i]
    w  <- grid$ratio[i] * s
    cl <- grid$conf_level[i]
    n  <- planned_N(true_sd = s, width = w, conf_level = cl)
    info <- sprintf("sd = %g, width = %g, conf_level = %g", s, w, cl)
    expect_lte(expected_width(n, s, cl), w,
               label = paste0("E(W) at N (", info, ")"))
    if (n > 2) expect_gt(expected_width(n - 1, s, cl), w,
                         label = paste0("E(W) at N - 1 (", info, ")"))
  }

  # The pilot of the help page: the control children's IQ gains in
  # pygmalion, planned for an expected width of 4 points. The bracket is
  # computed here; the pin keeps the number the help page and NEWS quote.
  gain <- pygmalion$iq_gain[pygmalion$treatment == "Control"]
  expect_lte(expected_width(172, sd(gain), 0.95), 4)
  expect_gt(expected_width(171, sd(gain), 0.95), 4)
  expect_equal(planned_N(true_sd = sd(gain), width = 4), 172)
})

test_that("c4(N) agrees with E(S) / sigma by quadrature and with its asymptotic expansion", {
  # E(S) = sigma E{sqrt(X / (n - 1))} for X chi square on n - 1 df, by
  # numerical integration over the bulk of that distribution.
  for (n in c(2, 3, 5, 10, 30, 200)) {
    df <- n - 1
    lim <- stats::qchisq(c(1e-15, 1 - 1e-15), df)
    e_s <- stats::integrate(function(x) sqrt(x / df) * stats::dchisq(x, df),
                            lim[1], lim[2], rel.tol = 1e-12)$value
    expect_equal(DMAR:::.c4(n), e_s, tolerance = 1e-9, info = paste("n =", n))
  }
  # For a large n the expansion 1 - 1/(4n) - 7/(32n^2) - 19/(128n^3) is
  # accurate to O(n^-4); the difference of two lgamma() values is not (it
  # returns a c4 above 1 at n = 1e8).
  for (n in c(1e4, 1e6, 1e8, 1e12)) {
    expect_equal(DMAR:::.c4(n), 1 - 1 / (4 * n) - 7 / (32 * n^2) - 19 / (128 * n^3),
                 tolerance = 1e-14, info = paste("n =", n))
  }
})

test_that("with an assurance, P(W <= width) is at least the assurance at necessary_N and below it at necessary_N - 1", {
  skip_on_cran()
  grid <- expand.grid(ratio = c(0.05, 0.3, 1.5),
                      conf_level = c(0.90, 0.95, 0.99),
                      assurance = c(0.55, 0.80, 0.99))
  for (i in seq_len(nrow(grid))) {
    w  <- grid$ratio[i]
    cl <- grid$conf_level[i]
    g  <- grid$assurance[i]
    n  <- planned_N(true_sd = 1, width = w, conf_level = cl, assurance = g)
    info <- sprintf("width = %g, conf_level = %g, assurance = %g", w, cl, g)
    expect_gte(prob_narrow(n, 1, w, cl), g,
               label = paste0("P(W <= width) at N (", info, ")"))
    if (n > 2) expect_lt(prob_narrow(n - 1, 1, w, cl), g,
                         label = paste0("P(W <= width) at N - 1 (", info, ")"))
  }

  # The help page's assurance plans for the pygmalion pilot, which the
  # sensitivity page's example compares (six children apart). Each is
  # bracketed here from its own probability, with the pilot's kurtosis
  # computed by the moment formula above rather than by kurtosis(), and the
  # pins keep the numbers the help pages and NEWS quote.
  gain <- pygmalion$iq_gain[pygmalion$treatment == "Control"]
  s_gain <- sd(gain)
  expect_gte(prob_narrow(195, s_gain, 4, 0.95), 0.90)
  expect_lt(prob_narrow(194, s_gain, 4, 0.95), 0.90)
  expect_equal(planned_N(true_sd = s_gain, width = 4, assurance = 0.90), 195)

  k_gain <- excess_kurtosis_G2(gain)
  expect_equal(kurtosis(gain), k_gain, tolerance = 1e-12)
  expect_equal(k_gain, 1.2, tolerance = 0.05)
  expect_gte(prob_narrow_k(201, s_gain, 4, 0.95, k_gain), 0.90)
  expect_lt(prob_narrow_k(200, s_gain, 4, 0.95, k_gain), 0.90)
  expect_true(all(prob_narrow_k(201:5000, s_gain, 4, 0.95, k_gain) >= 0.90))
  expect_equal(scan_N(s_gain, 4, 0.95, 0.90, k_gain), 201)
  expect_equal(planned_N(true_sd = s_gain, width = 4, assurance = 0.90,
                         true_kurtosis = k_gain), 201)
})

test_that("necessary_N is the smallest N from which the criterion holds at every larger N", {
  skip_on_cran()
  # With a large kurtosis the scaled chi square can meet the assurance at the
  # smallest N, miss it over a stretch above, and meet it again. Here it
  # holds at N = 2 and 3, fails at 4 to 6, and holds from 7 on, so the
  # answer is 7, not 2.
  p <- prob_narrow_k(2:7, 1, 1.6, 0.80, 100)
  expect_true(all(p[1:2] >= 0.90))
  expect_true(all(p[3:5] < 0.90))
  expect_true(all(prob_narrow_k(7:5000, 1, 1.6, 0.80, 100) >= 0.90))
  expect_equal(planned_N(true_sd = 1, width = 1.6, conf_level = 0.80,
                         assurance = 0.90, true_kurtosis = 100), 7)
  # A larger kurtosis lifts the whole stretch above the assurance, and the
  # plan is then N = 2: in this region a larger kurtosis can give a smaller
  # plan, as the help page says.
  expect_true(all(prob_narrow_k(2:5000, 1, 1.6, 0.80, 150) >= 0.90))
  expect_equal(planned_N(true_sd = 1, width = 1.6, conf_level = 0.80,
                         assurance = 0.90, true_kurtosis = 150), 2)

  # The brute-force scan of the definition, over cases that include an
  # assurance close to one half and kurtosis values up to 1e4 (a lognormal
  # whose log has standard deviation 1.5 has an excess kurtosis near that).
  cases <- list(
    list(w = 0.30, cl = 0.95, g = NULL, k = 0),
    list(w = 3.00, cl = 0.99, g = NULL, k = 0),
    list(w = 0.40, cl = 0.90, g = 0.501, k = 0),
    list(w = 0.20, cl = 0.95, g = 0.95, k = -1.9),
    list(w = 0.50, cl = 0.95, g = 0.80, k = 6),
    list(w = 0.2760482, cl = 0.90, g = 0.5353618, k = 50),
    list(w = 0.60, cl = 0.80, g = 0.51, k = 100),
    list(w = 1.20, cl = 0.95, g = 0.90, k = 20),
    list(w = 1.60, cl = 0.80, g = 0.90, k = 100),
    list(w = 1.00, cl = 0.95, g = 0.95, k = 1500),
    list(w = 0.3126, cl = 0.90, g = 0.95, k = 1e4)
  )
  for (cs in cases) {
    info <- sprintf("width = %g, conf_level = %g, assurance = %s, true_kurtosis = %g",
                    cs$w, cs$cl, format(cs$g), cs$k)
    expect_equal(planned_N(true_sd = 1, width = cs$w, conf_level = cs$cl,
                           assurance = cs$g, true_kurtosis = cs$k),
                 scan_N(1, cs$w, cs$cl, cs$g, cs$k), info = info)
  }
  # In the last two the stretch of failures sits well above N = 2.
  expect_equal(scan_N(1, 1.00, 0.95, 0.95, 1500), 89)
  expect_equal(scan_N(1, 0.3126, 0.90, 0.95, 1e4), 629)
})

test_that("a brute-force scan agrees over a grid of standard deviations, widths, levels, assurances, and kurtosis values", {
  skip_on_cran()  # The full grid runs in the local suite.
  # About one case in nine has a criterion that meets the assurance at some
  # N below a stretch where it fails, so the grid exercises the definition
  # and not only the bracket at N - 1.
  grid <- expand.grid(sd = c(0.4, 7), ratio = c(0.15, 0.5, 1.6),
                      conf_level = c(0.50, 0.80, 0.95, 0.99),
                      assurance = c(0.501, 0.6, 0.8, 0.9, 0.95, 0.99),
                      kurtosis = c(-1.99, -1, 0, 1.2, 6, 30, 100, 500, 1500,
                                   1e4))
  got <- want <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    x <- grid[i, ]
    w <- x$ratio * x$sd
    got[i] <- planned_N(true_sd = x$sd, width = w, conf_level = x$conf_level,
                        assurance = x$assurance, true_kurtosis = x$kurtosis)
    want[i] <- scan_N(x$sd, w, x$conf_level, x$assurance, x$kurtosis)
  }
  expect_equal(got, want)

  # The expected width plan over the same standard deviations, widths, and
  # levels.
  ew_grid <- unique(grid[, c("sd", "ratio", "conf_level")])
  got <- want <- numeric(nrow(ew_grid))
  for (i in seq_len(nrow(ew_grid))) {
    x <- ew_grid[i, ]
    got[i] <- planned_N(true_sd = x$sd, width = x$ratio * x$sd,
                        conf_level = x$conf_level)
    want[i] <- scan_N(x$sd, x$ratio * x$sd, x$conf_level)
  }
  expect_equal(got, want)
})

test_that("under normality the assurance criterion, once met, holds at every larger N", {
  skip_on_cran()  # A dense check behind a statement in Details.
  # The criterion is (width / (2 sd))^2 a(N) >= Q(N - 1), with
  # a(N) = N / t^2 and Q(nu) = qchisq(assurance, nu) / nu. It holds at every
  # N beyond the first for every width and sd exactly when a(N) / Q(N - 1)
  # never falls. Q peaks at nu below 1415 for every assurance here and falls
  # after, while a(N) rises, so the scan reaches past every peak.
  n <- 2:3000
  for (g in c(0.51, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999)) {
    for (cl in c(0.5, 0.8, 0.9, 0.95, 0.99, 0.999)) {
      t <- stats::qt(1 - (1 - cl) / 2, n - 1)
      r <- log(n / t^2) - log(stats::qchisq(g, n - 1) / (n - 1))
      expect_true(all(diff(r) > 0),
                  info = sprintf("assurance = %g, conf_level = %g", g, cl))
    }
  }
})

test_that("the kurtosis adjusts the assurance plan in the direction of the variance of S^2", {
  base <- planned_N(true_sd = 1, width = 0.3, assurance = 0.90)
  expect_equal(planned_N(true_sd = 1, width = 0.3, assurance = 0.90, true_kurtosis = 0),
               base)
  expect_lt(planned_N(true_sd = 1, width = 0.3, assurance = 0.90, true_kurtosis = -1.2),
            base)
  expect_gt(planned_N(true_sd = 1, width = 0.3, assurance = 0.90, true_kurtosis = 3),
            base)
  # The expected width plan does not use it, and says so.
  expect_warning(n_k <- planned_N(true_sd = 1, width = 0.3, true_kurtosis = 3),
                 "unused when 'assurance' is NULL")
  expect_equal(n_k, planned_N(true_sd = 1, width = 0.3))
  # With an assurance close to one half the criterion turns on the median
  # of S^2, and a larger kurtosis gives a slightly smaller plan.
  expect_lt(planned_N(true_sd = 1, width = 1, assurance = 0.51, true_kurtosis = 10),
            planned_N(true_sd = 1, width = 1, assurance = 0.51))
})

test_that("with an assurance of .80 or more, no smaller kurtosis gives a larger plan than one of at least kurtosis + 3", {
  skip_on_cran()  # A grid behind two statements in the help page.
  ks <- c(seq(-1.99, -0.2, by = 0.3), 0:10, seq(15, 100, by = 5),
          seq(200, 3000, by = 200))
  for (g in c(0.80, 0.90, 0.95, 0.99)) {
    for (cl in c(0.5, 0.95, 0.99)) {
      for (w in c(0.1, 0.4, 1.2, 4)) {
        info <- sprintf("assurance = %g, conf_level = %g, width = %g", g, cl, w)
        N_k <- vapply(ks, function(k) {
          planned_N(true_sd = 1, width = w, conf_level = cl, assurance = g,
                    true_kurtosis = k)
        }, numeric(1))
        # The kurtosis values whose plan is at least kurtosis + 3 but below
        # the plan of some smaller kurtosis; there must be none.
        exceeded <- vapply(seq_along(ks), function(j) {
          N_k[j] >= ks[j] + 3 && any(N_k[seq_len(j)] > N_k[j])
        }, logical(1))
        expect_equal(ks[exceeded], numeric(0), info = info)
        expect_true(any(N_k >= ks + 3), info = info)
        # A planning value of 0, in place of a pilot's value at or below -2,
        # never gives a smaller plan than -1.9, the alternative the help
        # page names.
        expect_gte(planned_N(true_sd = 1, width = w, conf_level = cl,
                             assurance = g),
                   planned_N(true_sd = 1, width = w, conf_level = cl,
                             assurance = g, true_kurtosis = -1.9), label = info)
      }
    }
  }
})

test_that("the plan grows with the standard deviation, the level, and the assurance and shrinks with the width", {
  expect_lt(planned_N(true_sd = 1, width = 0.5), planned_N(true_sd = 2, width = 0.5))
  expect_gt(planned_N(true_sd = 1, width = 0.5), planned_N(true_sd = 1, width = 1))
  expect_lt(planned_N(true_sd = 1, width = 0.5, conf_level = 0.90),
            planned_N(true_sd = 1, width = 0.5, conf_level = 0.99))
  expect_lt(planned_N(true_sd = 1, width = 0.5, assurance = 0.80),
            planned_N(true_sd = 1, width = 0.5, assurance = 0.99))
  # The plan depends on sd and width only through their ratio.
  expect_equal(planned_N(true_sd = 3, width = 1.5), planned_N(true_sd = 1, width = 0.5))
})

test_that("a very wide target is met at N = 2 and a very narrow one at a very large N", {
  expect_equal(planned_N(true_sd = 1, width = 100), 2)
  expect_equal(planned_N(true_sd = 1, width = 100, assurance = 0.99), 2)
  # At a size in the billions the gamma ratio is evaluated through its
  # asymptotic expansion here, since the lgamma() difference has lost the
  # digits that separate N from N - 1.
  ew_large <- function(n) {
    2 * stats::qt(0.975, n - 1) * (1 - 1 / (4 * n) - 7 / (32 * n^2)) / sqrt(n)
  }
  n <- planned_N(true_sd = 1, width = 1e-4)
  expect_gt(n, 1e9)
  expect_lte(ew_large(n), 1e-4)
  expect_gt(ew_large(n - 1), 1e-4)
  expect_error(ss_aipe_mean(true_sd = 1, width = 1e-9), "exceed 1e15")
})

test_that("ss_aipe_mean() validates its inputs", {
  skip_on_cran()
  expect_error(ss_aipe_mean(true_sd = -1, width = 1), "'true_sd'")
  expect_error(ss_aipe_mean(true_sd = c(1, 2), width = 1), "'true_sd'")
  expect_error(ss_aipe_mean(true_sd = NA_real_, width = 1), "'true_sd'")
  expect_error(ss_aipe_mean(true_sd = 1, width = 0), "'width'")
  expect_error(ss_aipe_mean(true_sd = 1, width = Inf), "'width'")
  expect_error(ss_aipe_mean(true_sd = 1, width = 1, conf_level = 1), "'conf_level'")
  expect_error(ss_aipe_mean(true_sd = 1, width = 1, conf_level = 0), "'conf_level'")
  expect_error(ss_aipe_mean(true_sd = 1, width = 1, assurance = 0.5), "'assurance'")
  expect_error(ss_aipe_mean(true_sd = 1, width = 1, assurance = 1), "'assurance'")
  expect_error(ss_aipe_mean(true_sd = 1, width = 1, assurance = 0.9, true_kurtosis = -2),
               "'true_kurtosis'")
  # kurtosis() corrects for small-sample bias, so a pilot of two equally
  # frequent values gives an estimate below -2, which is refused with a
  # pointer to the help page.
  k_two <- kurtosis(rep(c(0, 1), 123))
  expect_lt(k_two, -2)
  expect_equal(k_two, excess_kurtosis_G2(rep(c(0, 1), 123)), tolerance = 1e-12)
  expect_error(ss_aipe_mean(true_sd = 1, width = 1, assurance = 0.9,
                            true_kurtosis = k_two), "-2 or less")
  expect_error(ss_aipe_mean(true_sd = 1, width = 1, assurance = 0.9,
                            true_kurtosis = NA_real_), "'true_kurtosis'")
})
