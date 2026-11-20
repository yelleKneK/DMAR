## ss_aipe_two_phase(): the labeled units for the regression interval of
## ci_two_phase(), checked against a brute-force search written here from the
## criterion itself, the target-only row against an independent computation of
## E(S) and against ss_aipe_mean(), the pool cap, and the refusals.

# The criterion, restated independently of the package: TRUE when n labeled
# units meet the goal. df_lost is 2 for the two-phase interval and 1 for the
# target alone (where rho = 0 and the pool term vanishes). Without an
# assurance the two-phase plan sets its estimates at their planning values,
# while the target-only plan uses the exact expected width of the t interval,
# E(S) = c4(n) sigma, written here through lgamma() (the package uses
# lbeta()).
.meets_two_phase <- function(n, sd_y, rho, width, M, conf_level, assurance,
                             kurtosis, df_lost = 2) {
  df <- n - df_lost
  t <- qt(1 - (1 - conf_level) / 2, df)
  fixed <- if (is.finite(M)) rho^2 * sd_y^2 / M else 0
  if (is.null(assurance)) {
    if (df_lost == 1) {
      c4 <- sqrt(2 / (n - 1)) * exp(lgamma(n / 2) - lgamma((n - 1) / 2))
      return(2 * t * c4 * sd_y / sqrt(n) <= width)
    }
    return(2 * t * sqrt(sd_y^2 * (1 - rho^2) / n + fixed) <= width)
  }
  room <- (width / (2 * t))^2 - fixed
  if (room <= 0) return(FALSE)
  nu <- 2 / (2 / df + kurtosis / n)
  pchisq(nu * n * room / (sd_y^2 * (1 - rho^2)), nu) >= assurance
}

# The smallest n from which the criterion holds at every larger n, found by
# evaluating it at every n from the smallest the interval allows (one degree
# of freedom) up to n_hi, which must lie where the criterion has settled.
.brute_n <- function(..., df_lost = 2, n_hi) {
  n <- (df_lost + 1):n_hi
  ok <- vapply(n, function(k) .meets_two_phase(k, ..., df_lost = df_lost),
               logical(1))
  stopifnot(all(ok[n > n_hi / 2]))
  if (all(ok)) df_lost + 1 else max(n[!ok]) + 1
}

.v <- function(res) stats::setNames(res$value, res$term)

test_that("the target-only expected width uses E(S) for a normal sample", {
  # E(S) / sigma computed by integrating over the chi square distribution of
  # (n - 1) S^2 / sigma^2, independently of the closed form. The range runs
  # between quantiles that leave out less than 1e-14 of the distribution.
  for (n in c(2, 3, 5, 10, 40, 250)) {
    df <- n - 1
    e_s <- integrate(function(u) sqrt(u / df) * dchisq(u, df),
                     qchisq(1e-14, df), qchisq(1e-14, df, lower.tail = FALSE),
                     rel.tol = 1e-10)$value
    expect_equal(DMAR:::.c4(n), e_s, tolerance = 1e-8,
                 info = paste("n =", n))
  }
  # At a large n the closed form stays accurate (no cancellation) and below 1.
  big <- DMAR:::.c4(1e8)
  expect_lt(big, 1)
  expect_equal(big, 1 - 1 / (4 * 1e8), tolerance = 1e-12)
  # Setting S equal to sigma in the width, as the plug-in rule does, gives
  # 123 here; the exact expected width gives 122.
  plan <- .v(ss_aipe_two_phase(true_sd_Y = 1, true_cor = 0.5, width = 0.3,
                               conf_level = 0.90))
  expect_equal(plan[["necessary_n_target_only"]], 122)
  w <- function(n) 2 * qt(0.95, n - 1) * DMAR:::.c4(n) / sqrt(n)
  expect_lte(w(122), 0.3)
  expect_gt(w(121), 0.3)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

# Planning values of an applied evaluation (target standard deviations and
# target-proxy correlations of two systems, and the residual excess kurtosis
# of the first), with the sizes found when the method was developed. These
# are regression anchors: they pin the planner's own earlier output, so that
# any change in it shows, and are not independent values. The independent
# checks are the brute-force search, the integral for E(S), and the
# comparison with ss_aipe_mean() below.
vic <- c(sd_y = 0.427741256538497, rho = 0.549109594316365,
         kurt_e = -0.27048375367956)
g4  <- c(sd_y = 0.377486619420778, rho = 0.269917091129373)

test_that("the planned sizes reproduce the anchors", {
  f <- function(p, M, g = 0.80, k = 0) {
    .v(ss_aipe_two_phase(p[["sd_y"]], p[["rho"]], width = 0.10, n_pool = M,
                         conf_level = 0.90, assurance = g, true_kurtosis = k))
  }
  expect_equal(f(vic, 20000)[["necessary_n_labeled"]], 154)
  expect_equal(f(vic, 20000)[["necessary_n_target_only"]], 216)
  expect_equal(f(vic, 2000)[["necessary_n_labeled"]], 158)
  expect_equal(f(vic, 20000, k = vic[["kurt_e"]])[["necessary_n_labeled"]], 153)
  expect_equal(f(g4, 20000)[["necessary_n_labeled"]], 159)
  expect_equal(f(g4, 2000)[["necessary_n_labeled"]], 160)
  expect_equal(f(g4, 20000)[["necessary_n_target_only"]], 170)
  expect_equal(f(vic, 20000, g = NULL)[["necessary_n_labeled"]], 141)
})

test_that("the planned size is the smallest from which the criterion holds, over a grid", {
  skip_on_cran()
  # Assurances below one half and a large kurtosis are included: there the
  # normal-quantile solution is not a lower bound, and with kurtosis 20 the
  # two-moment approximation meets a low assurance at a few small n before
  # failing over a stretch of larger ones.
  settings <- list(list(NULL, 0), list(0.1, 0), list(0.1, 20), list(0.5, 0),
                   list(0.5, 20), list(0.9, -1), list(0.9, 0), list(0.9, 2),
                   list(0.9, 20))
  for (rho in c(-0.7, 0.3, 0.9)) {
    for (M in c(500, 5000, Inf)) {
      for (st in settings) {
        assurance <- st[[1]]; kurtosis <- st[[2]]
        args <- list(sd_y = 1.5, rho = rho, width = 0.5, M = M,
                     conf_level = 0.95, assurance = assurance,
                     kurtosis = kurtosis)
        info <- sprintf("rho = %.1f, M = %s, assurance = %s, kurtosis = %g",
                        rho, M, format(assurance), kurtosis)
        plan <- .v(ss_aipe_two_phase(true_sd_Y = 1.5, true_cor = rho, width = 0.5,
                                     n_pool = M, assurance = assurance,
                                     true_kurtosis = kurtosis))
        n <- plan[["necessary_n_labeled"]]
        expect_equal(n, do.call(.brute_n, c(args, n_hi = 2 * n + 40)),
                     info = info)
        expect_true(do.call(.meets_two_phase, c(list(n = n), args)),
                    info = info)
        # Below the plan the criterion fails, unless the plan is already the
        # smallest size the interval allows (3, one degree of freedom).
        if (n > 3)
          expect_false(do.call(.meets_two_phase, c(list(n = n - 1), args)),
                       info = info)
        expect_equal(plan[["ratio"]], n / plan[["necessary_n_target_only"]])
      }
    }
  }
  # The target-only row, which depends on neither rho nor the pool.
  for (st in settings) {
    for (cl in c(0.90, 0.99)) {
      assurance <- st[[1]]; kurtosis <- st[[2]]
      plan <- .v(ss_aipe_two_phase(true_sd_Y = 1.5, true_cor = 0.5, width = 0.5,
                                   conf_level = cl, assurance = assurance,
                                   true_kurtosis = kurtosis))
      n_t <- plan[["necessary_n_target_only"]]
      expect_equal(n_t, .brute_n(sd_y = 1.5, rho = 0, width = 0.5, M = Inf,
                                 conf_level = cl, assurance = assurance,
                                 kurtosis = kurtosis, df_lost = 1,
                                 n_hi = 2 * n_t + 40),
                   info = sprintf("target only, %s, %g, %g", format(assurance),
                                  kurtosis, cl))
    }
  }
})

test_that("the target-only row is the plan ss_aipe_mean() makes", {
  skip_on_cran()
  # ss_aipe_mean() plans the same single-mean interval on its own; the two
  # must agree in every cell, with and without an assurance.
  ss_aipe_mean <- get0("ss_aipe_mean", envir = asNamespace("DMAR"),
                       inherits = FALSE)
  skip_if(is.null(ss_aipe_mean), "ss_aipe_mean() is not in this build")
  for (sd in c(0.43, 1, 2)) {
    for (w in c(0.1, 0.3, 0.5, 1, 20)) {
      for (cl in c(0.90, 0.95)) {
        for (g in list(NULL, 0.8, 0.95)) {
          for (k in if (is.null(g)) 0 else c(-1.5, 0, 3, 20)) {
            a <- .v(ss_aipe_two_phase(sd, 0.5, w, conf_level = cl,
                                      assurance = g, true_kurtosis = k))
            b <- .v(ss_aipe_mean(sd, w, conf_level = cl, assurance = g,
                                 true_kurtosis = k))
            expect_equal(a[["necessary_n_target_only"]], b[["necessary_N"]],
                         info = sprintf("sd %g, width %g, conf %g, %s, %g",
                                        sd, w, cl, format(g), k))
          }
        }
      }
    }
  }
})

test_that("the expected width and the floor are reported at the plan", {
  plan <- .v(ss_aipe_two_phase(true_sd_Y = 2, true_cor = 0.7, width = 0.6, n_pool = 3000,
                               conf_level = 0.90))
  n <- plan[["necessary_n_labeled"]]
  expect_equal(plan[["expected_width"]],
               2 * qt(0.95, n - 2) * sqrt(4 * (1 - 0.49) / n + 0.49 * 4 / 3000))
  expect_lte(plan[["expected_width"]], 0.6)
  # The floor: the normal-quantile width with all but one unit labeled.
  expect_equal(plan[["minimum_width"]],
               2 * qnorm(0.95) * sqrt(4 * (1 - 0.49) / 2999 + 0.49 * 4 / 3000))
  expect_equal(.v(ss_aipe_two_phase(2, 0.7, 0.6))[["minimum_width"]], 0)
  expect_false(any(c("assurance", "true_kurtosis", "kurtosis") %in% names(plan)))
  with_a <- .v(ss_aipe_two_phase(2, 0.7, 0.6, assurance = 0.9, true_kurtosis = -1))
  expect_equal(with_a[["assurance"]], 0.9)
  expect_equal(with_a[["true_kurtosis"]], -1)
})

test_that("a pool too small for the width is an error, not a missing value", {
  expect_error(ss_aipe_two_phase(0.4, 0.9, width = 0.01, n_pool = 100,
                                 conf_level = 0.90), "pool is too small")
  # Just above the floor the normal quantile could reach the width, but with
  # the t quantile and an assurance the search passes the pool.
  floor_w <- 2 * qnorm(0.95) * 0.4 / sqrt(200)
  expect_error(ss_aipe_two_phase(0.4, 0.5, width = 1.01 * floor_w,
                                 n_pool = 200, conf_level = 0.90,
                                 assurance = 0.90),
               "passes the pool")
  # A wider target fits in the same pool, with units left unlabeled.
  ok <- .v(ss_aipe_two_phase(0.4, 0.5, width = 2 * floor_w, n_pool = 200,
                             conf_level = 0.90))
  expect_lt(ok[["necessary_n_labeled"]], 200)
})

test_that("a plan always leaves at least one unit of the pool unlabeled", {
  skip_on_cran()
  # Labeling all 50 units of a pool of 50 would meet this width (the t
  # quantile on 48 degrees of freedom is what the width was built from), but
  # ci_two_phase() needs a unit with the proxy alone, and 49 labeled units
  # do not meet it.
  w <- 1.0001 * 2 * qt(0.975, 48) / sqrt(50)
  expect_error(ss_aipe_two_phase(1, 0.5, w, n_pool = 50), "passes the pool")
  expect_error(ss_aipe_two_phase_sensitivity(true_sd_Y = 1, true_cor = 0.5,
                                             supposed_sd_Y = 1,
                                             supposed_cor = 0.5, n_pool = 50,
                                             width = w, G = 2),
               "passes the pool")
  # Over widths from just above the floor upward, a plan is either an error
  # or a design that ci_two_phase() can analyze.
  for (M in c(10, 50, 200)) {
    floor_w <- .v(ss_aipe_two_phase(1, 0.5, 10, n_pool = M))[["minimum_width"]]
    for (f in c(1.0001, 1.001, 1.01, 1.05, 1.2, 2)) {
      for (g in list(NULL, 0.8)) {
        plan <- tryCatch(.v(ss_aipe_two_phase(1, 0.5, f * floor_w, n_pool = M,
                                              assurance = g)),
                         error = function(e) conditionMessage(e))
        info <- sprintf("M = %d, f = %g, %s", M, f, format(g))
        if (is.character(plan)) {
          expect_match(plan, "passes the pool", info = info)
        } else {
          n <- plan[["necessary_n_labeled"]]
          expect_lte(n, M - 1, label = info)
          set.seed(M)
          x <- rnorm(M)
          y <- 0.5 * x + rnorm(M)
          expect_s3_class(ci_two_phase(y[1:n], x[1:n], x[-(1:n)]), "dmar_tbl")
        }
      }
    }
  }
})

test_that("heavier residual tails need more labeled units for the same assurance", {
  n_k <- vapply(c(-1.2, 0, 3), function(k)
    .v(ss_aipe_two_phase(1, 0.5, 0.3, n_pool = 10000, assurance = 0.9,
                         true_kurtosis = k))[["necessary_n_labeled"]], numeric(1))
  expect_true(n_k[1] < n_k[2] && n_k[2] < n_k[3])
  expect_warning(ss_aipe_two_phase(1, 0.5, 0.3, true_kurtosis = 1),
                 "only a plan with an 'assurance'")
})

test_that("the planner joins the AIPE tidiers through the size registry", {
  res <- ss_aipe_two_phase(true_sd_Y = 0.43, true_cor = 0.55, width = 0.10,
                           n_pool = 20000, conf_level = 0.90)
  expect_s3_class(res, "dmar_ss_aipe")
  td <- generics::tidy(res)
  expect_equal(td$estimate, res$value[res$term == "necessary_n_labeled"])
  expect_equal(td$width, 0.10)
  gl <- generics::glance(res)
  expect_equal(gl$necessary_n_target_only,
               res$value[res$term == "necessary_n_target_only"])
})

test_that("malformed planning inputs are refused", {
  expect_error(ss_aipe_two_phase(-1, 0.5, 0.2), "'true_sd_Y'")
  expect_error(ss_aipe_two_phase(1, 1, 0.2), "'true_cor'")
  expect_error(ss_aipe_two_phase(1, 0.5, 0), "width")
  expect_error(ss_aipe_two_phase(1, 0.5, 0.2, conf_level = 95), "conf_level")
  expect_error(ss_aipe_two_phase(1, 0.5, 0.2, assurance = 1), "assurance")
  expect_error(ss_aipe_two_phase(1, 0.5, 0.2, assurance = 0.8, true_kurtosis = -3),
               "'true_kurtosis'")
  expect_error(ss_aipe_two_phase(1, 0.5, 0.2, n_pool = 1000.5), "n_pool")
  expect_error(ss_aipe_two_phase(1, 0.5, 0.2, n_pool = 3), "n_pool")
})
