test_that("ss_aipe_smd() returns a tidy data frame led by necessary_n_per_group", {
  result <- ss_aipe_smd(true_smd = 0.5, conf_level = 0.95, width = 0.4)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term[1], "necessary_n_per_group")
  expect_true(all(c("true_smd", "width") %in% result$term))
  expect_false("delta" %in% result$term)   # denoted as supposed, not bare delta
})

test_that("ss_aipe_smd() returns a positive integer-valued sample size", {
  r <- ss_aipe_smd(true_smd = 0.5, conf_level = 0.95, width = 0.4)
  n <- r$value[r$term == "necessary_n_per_group"]
  expect_gt(n, 0)
  expect_equal(n, as.integer(n))
})

test_that("ss_aipe_smd() requires a larger n for a narrower width", {
  n_of <- function(w) {
    r <- ss_aipe_smd(true_smd = 0.5, conf_level = 0.95, width = w)
    r$value[r$term == "necessary_n_per_group"]
  }
  expect_lt(n_of(0.5), n_of(0.3))
})

test_that("ss_aipe_smd() echoes the planning inputs as numeric rows", {
  r <- ss_aipe_smd(true_smd = 0.5, conf_level = 0.95, width = 0.4)
  expect_equal(r$value[r$term == "true_smd"], 0.5)
  expect_equal(r$value[r$term == "width"], 0.4)
  expect_false("assurance" %in% r$term)          # not supplied -> no row
  expect_type(r$value, "double")                 # value column stays numeric
  ra <- ss_aipe_smd(true_smd = 0.5, conf_level = 0.95, width = 0.4, assurance = 0.9)
  expect_equal(ra$value[ra$term == "assurance"], 0.9)
})

test_that("ss_aipe_smd() floors an impossibly wide target at the admissible minimum", {
  # A huge target width previously crashed the fixed-point search ("missing
  # value where TRUE/FALSE needed") because the central-t error df 2n - 2 went
  # negative. The smallest admissible design this planner returns is 4 per
  # group.
  r <- ss_aipe_smd(true_smd = 0.5, width = 100)
  expect_equal(r$value[r$term == "necessary_n_per_group"], 4)
})

test_that("ss_aipe_smd() validates its boundary inputs", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  expect_error(ss_aipe_smd(true_smd = 0.5, width = 0), "width")
  expect_error(ss_aipe_smd(true_smd = 0.5, width = -1), "width")
  expect_error(ss_aipe_smd(true_smd = 0.5, width = 0.3, conf_level = 1.2),
               "conf_level")
  expect_error(ss_aipe_smd(true_smd = NA_real_, width = 0.3), "true_smd")
})

# ---- The assurance step: group reversal and the attained probability --------
#
# Oracles coded apart from the package. The width event: at n per group the
# interval ci_smd() reports widens as |d| grows, so it is no wider than the
# target exactly when |d| <= d*, found here from ci_smd() widths with a
# tighter tolerance. Its probability is computed two ways, neither of which
# the package uses: (a) R's noncentral t with a nonnegative noncentrality,
# each tail taken as a small probability (accurate to about 1e-12 at the
# noncentralities used here, below 30 with at most 2,000 degrees of
# freedom; see ?ci_nc_t for where it is not); (b) conditioning on the
# normal variate,
#   P(|Z + lambda| <= k S) = E_Z[ P(V >= df (Z + lambda)^2 / k^2) ],
# with V chi square on df = 2n - 2 and S = sqrt(V / df). The package
# integrates over the chi square instead.
.oracle_smd_width <- function(d, n, conf_level = 0.95) {
  ci <- ci_smd(smd = d, n_1 = n, n_2 = n, conf_level = conf_level, tol = 1e-12)
  ci$value[ci$term == "ci_upper"] - ci$value[ci$term == "ci_lower"]
}
.oracle_smd_cutoff <- function(n, width) {
  if (.oracle_smd_width(0, n) > width) return(0)
  hi <- 1
  while (.oracle_smd_width(hi, n) <= width) hi <- 2 * hi
  stats::uniroot(function(d) .oracle_smd_width(d, n) - width, c(0, hi),
                 tol = 1e-12)$root
}
.oracle_smd_pt <- function(cut, n, delta) {
  s <- sqrt(n / 2)
  df <- 2 * n - 2
  lambda <- abs(delta) * s
  1 - stats::pt(cut * s, df, ncp = lambda, lower.tail = FALSE) -
    stats::pt(-cut * s, df, ncp = lambda)
}
.oracle_smd_z <- function(cut, n, delta) {
  if (cut <= 0) return(0)
  df <- 2 * n - 2
  k <- cut * sqrt(n / 2)
  lambda <- abs(delta) * sqrt(n / 2)
  f <- function(z) stats::dnorm(z) *
    stats::pchisq(df * (z + lambda)^2 / k^2, df, lower.tail = FALSE)
  edges <- sort(unique(c(-12, 12, (c(-lambda - k, -lambda, -lambda + k))[
    abs(c(-lambda - k, -lambda, -lambda + k)) < 12])))
  sum(vapply(seq_len(length(edges) - 1L), function(i)
    stats::integrate(f, edges[i], edges[i + 1L], rel.tol = 1e-12,
                     subdivisions = 2000L)$value, numeric(1)))
}

test_that("the two oracle forms agree with each other where R's noncentral t is exact", {
  for (n in c(10, 130, 500)) {
    for (delta in c(0, 0.5, -1)) {
      cut <- .oracle_smd_cutoff(n, 2 * 1.96 * sqrt(2 / n) * 1.02)
      # Absolute agreement: the pt() form subtracts two tails from 1, and
      # its series stops at an error of 1e-12.
      expect_lt(abs(.oracle_smd_pt(cut, n, delta) -
                      .oracle_smd_z(cut, n, delta)), 1e-10,
                label = paste(n, delta))
    }
  }
})

test_that("the width probability of ss_aipe_smd() is the noncentral t probability, beyond 37.62 as well", {
  skip_on_cran()
  # ss_aipe_smd() computes P(|t| <= k) at n per group with the shared
  # .nc_t_abs_probability() of R/nc_t_internals.R (it had a copy of its
  # own, .smd_abs_t_probability(), before).
  prob <- function(k, n, delta)
    DMAR:::.nc_t_abs_probability(k, 2 * n - 2, abs(delta) * sqrt(n / 2))
  for (n in c(4, 30, 400)) {
    for (delta in c(0, 0.4, -1.2)) {
      for (cut in c(0.2, 0.9)) {
        k <- cut * sqrt(n / 2)
        expect_lt(abs(prob(k, n, delta) - .oracle_smd_pt(cut, n, delta)),
                  1e-10, label = paste(n, delta, cut))
        expect_lt(abs(prob(k, n, delta) - .oracle_smd_z(cut, n, delta)),
                  1e-12, label = paste(n, delta, cut))
      }
    }
  }
  # n = 800 and delta = 2 put the noncentrality at 40, where pt() is a
  # normal approximation; the conditioning form is the oracle there.
  expect_lt(abs(prob(42, 800, 2) - .oracle_smd_z(42 / sqrt(400), 800, 2)),
            1e-12)
  expect_identical(prob(0, 50, 1), 0)
  # The one-sided form, against the conditioning form, at 100,000 per group
  # and a noncentrality of 36, where R's series stops early.
  cdf_z <- function(t, df, ncp) {
    f <- function(z) stats::dnorm(z) * ifelse(z + ncp <= 0, 1,
      stats::pchisq(df * (z + ncp)^2 / t^2, df, lower.tail = FALSE))
    pts <- c(-12, t - ncp + c(-8, 0, 8) * t / sqrt(2 * df), 12)
    edges <- sort(unique(pts[abs(pts) <= 12]))
    sum(vapply(seq_len(length(edges) - 1L), function(i)
      stats::integrate(f, edges[i], edges[i + 1L], rel.tol = 1e-12,
                       subdivisions = 2000L)$value, numeric(1)))
  }
  for (t in c(34, 36, 38)) {
    expect_lt(abs(DMAR:::.nc_t_cdf(t, 2e5 - 2, 36) - cdf_z(t, 2e5 - 2, 36)),
              1e-10, label = paste("t", t))
  }
})

test_that("ss_aipe_smd() with assurance gives the same plan when the groups are reversed", {
  skip_on_cran()
  # DMAR 1.0.x gave 125 and 133 per group for -.5 and -1 at width .5 and
  # assurance .90 (below the plans without assurance, 127 and 139), and
  # 130 and 144 for .5 and 1.
  for (delta in c(0.5, 1)) {
    pos <- ss_aipe_smd(true_smd = delta, width = 0.5, assurance = 0.9)
    neg <- ss_aipe_smd(true_smd = -delta, width = 0.5, assurance = 0.9)
    expect_identical(neg$value[-2], pos$value[-2])
    expect_identical(neg$value[neg$term == "true_smd"], -delta)
  }
  expect_identical(
    ss_aipe_smd(true_smd = -0.5, width = 0.5, assurance = 0.9)$value[1], 130)
  expect_identical(
    ss_aipe_smd(true_smd = -1, width = 0.5, assurance = 0.9)$value[1], 144)
  for (delta in c(0.2, 2)) {
    for (a in c(0.8, 0.99)) {
      expect_identical(
        ss_aipe_smd(true_smd = -delta, width = 1, assurance = a)$value[-2],
        ss_aipe_smd(true_smd = delta, width = 1, assurance = a)$value[-2])
    }
  }
})

test_that("ss_aipe_smd() with assurance attains the request by the oracle and never plans below the plan without it", {
  skip_on_cran()  # numerical oracle check, run in the full local suite
  # Exact width probabilities at the 1.0.x plans and at the repaired ones,
  # computed in the review by the noncentral t: .146146166043 (125 per
  # group, |true_smd| = .5), .064400167775 (133, 1), .915626559175 (130,
  # .5), and .904165403522 (144, 1), all at width .5.
  expect_equal(.oracle_smd_pt(.oracle_smd_cutoff(125, 0.5), 125, -0.5),
               0.146146166043, tolerance = 1e-8)
  expect_equal(.oracle_smd_pt(.oracle_smd_cutoff(133, 0.5), 133, -1),
               0.064400167775, tolerance = 1e-8)
  for (delta in c(0, 0.5, -1)) {
    for (width in c(0.5, 1)) {
      for (a in c(0.8, 0.99)) {
        info <- sprintf("true_smd = %.1f, width = %.2f, assurance = %.2f",
                        delta, width, a)
        res <- ss_aipe_smd(true_smd = delta, width = width, assurance = a)
        n <- res$value[res$term == "necessary_n_per_group"]
        n_without <- ss_aipe_smd(true_smd = delta, width = width)$value[1]
        attained <- res$value[res$term == "assurance_attained"]
        cut <- .oracle_smd_cutoff(n, width)
        expect_equal(attained, .oracle_smd_pt(cut, n, delta),
                     tolerance = 1e-7, info = info)
        expect_equal(attained, .oracle_smd_z(cut, n, delta),
                     tolerance = 1e-7, info = info)
        expect_gte(attained, a)
        expect_gte(n, n_without)
        # The plan is the smallest: one fewer, when not below the plan
        # without assurance, falls short.
        if (n > n_without)
          expect_lt(.oracle_smd_z(.oracle_smd_cutoff(n - 1, width), n - 1,
                                  delta), a, label = info)
      }
    }
  }
  res <- ss_aipe_smd(true_smd = -0.5, width = 0.5, assurance = 0.9)
  expect_equal(res$value[res$term == "assurance_attained"], 0.915626559175,
               tolerance = 1e-8)
})

test_that("ss_aipe_smd() returns the smallest plan, not below the plan without assurance, whose exact probability attains the request", {
  skip_on_cran()  # an exhaustive scan by the oracle, run in the full local suite
  # The scan raises n one at a time from the plan without assurance and
  # stops at the first n whose probability, by the oracle that conditions
  # on the normal variate, reaches the request. Through 1.1.0's first
  # builds the default was the published procedure of Kelley and Rausch
  # (2006), whose plans MBESS::ss.aipe.smd() (MBESS 5.0.1, 2026-11-11) gives
  # as 513, 150, 53, 125, and 34 for these five sets: it attained the
  # request in each, and in the second and third it is one above the
  # smallest plan that does (149 attains .9903, and 52 attains .9008).
  cases <- list(c(0.5, 0.25, 0.90, 513), c(1, 0.5, 0.99, 149),
                c(2, 1, 0.90, 52), c(0, 0.5, 0.99, 125), c(0.2, 1, 0.99, 34))
  for (x in cases) {
    info <- paste(x[1], x[2], x[3])
    n_without <- ss_aipe_smd(true_smd = x[1], width = x[2])$value[1]
    m <- n_without
    while (.oracle_smd_z(.oracle_smd_cutoff(m, x[2]), m, x[1]) < x[3])
      m <- m + 1
    expect_identical(m, x[4], info = info)
    for (sgn in c(1, -1)) {
      res <- ss_aipe_smd(true_smd = sgn * x[1], width = x[2], assurance = x[3])
      expect_identical(res$value[1], x[4], info = paste(sgn, info))
    }
  }
  expect_equal(.oracle_smd_pt(.oracle_smd_cutoff(149, 0.5), 149, 1), 0.99026,
               tolerance = 1e-5)
  expect_equal(.oracle_smd_pt(.oracle_smd_cutoff(150, 0.5), 150, 1), 0.99430,
               tolerance = 1e-5)
})

test_that("ss_aipe_smd() plans 50 where the published procedure of DMAR 1.0.x fell short with 49", {
  # At true_smd = 0, width .8, and assurance .95 the published procedure
  # (MBESS::ss.aipe.smd() gives 49 per group) attains slightly less than
  # .95 by the oracle, and 50 attain .9946.
  expect_lt(.oracle_smd_z(.oracle_smd_cutoff(49, 0.8), 49, 0), 0.95)
  res <- ss_aipe_smd(true_smd = 0, width = 0.8, assurance = 0.95)
  expect_identical(res$value[1], 50)
  attained <- res$value[res$term == "assurance_attained"]
  expect_gte(attained, 0.95)
  expect_equal(attained, .oracle_smd_z(.oracle_smd_cutoff(50, 0.8), 50, 0),
               tolerance = 1e-7)
})

test_that("ss_aipe_smd() plans its benchmarks without noncentral t warnings", {
  skip_on_cran()
  # DMAR 1.0.x emitted 15 and 2 warnings of pt() ("full precision may not
  # have been achieved in 'pnt{final}'") at true_smd = -1 and -.5, from one
  # minus a probability near 1 in the assurance step.
  for (delta in c(-1, -0.5, 0.5, 1)) {
    expect_no_warning(ss_aipe_smd(true_smd = delta, width = 0.5,
                                  assurance = 0.9))
  }
})

test_that("ss_aipe_smd() plans exactly past R's noncentral t limit, without a warning", {
  skip_on_cran()  # numerical oracle check, run in the full local suite
  # true_smd = 2 with width .25 plans about 740 per group, where the
  # noncentrality is near 39. ci_nc_t() warned once per interval (31 times
  # in 1.0.x), and until the limits were exact the probability of an
  # interval no wider than the target could not be certified there
  # (assurance_attained was NA). The plans are 738, and 759 with an
  # assurance of .90.
  res <- expect_no_warning(ss_aipe_smd(true_smd = 2, width = 0.25))
  expect_identical(res$value[1], 738)
  res <- expect_no_warning(ss_aipe_smd(true_smd = 2, width = 0.25,
                                       assurance = 0.9))
  expect_identical(res$value[1], 759)
  attained <- res$value[res$term == "assurance_attained"]
  expect_gte(attained, 0.9)
  expect_equal(attained, 0.905763, tolerance = 1e-6)
  expect_equal(attained,
               .oracle_smd_z(.oracle_smd_cutoff(759, 0.25), 759, 2),
               tolerance = 1e-7)
  # true_smd = 2 with width .5 stays in the accurate range (192 per group,
  # an observed noncentrality near 21 at the cutoff).
  res <- expect_no_warning(ss_aipe_smd(true_smd = 2, width = 0.5,
                                       assurance = 0.8))
  expect_identical(res$value[1], 192)
  expect_equal(res$value[res$term == "assurance_attained"],
               .oracle_smd_z(.oracle_smd_cutoff(192, 0.5), 192, 2),
               tolerance = 1e-7)
})

test_that("the attained assurance of ss_aipe_smd() agrees with simulated two-group data", {
  skip_on_cran()
  # Two groups of n normal observations with standardized difference
  # true_smd; d from the pooled standard deviation and the interval from
  # ci_smd(). 1,000 samples; the standard error of a proportion near .92 is
  # about .0088. Each sample's interval is no wider than the target exactly
  # when |d| is at most the cutoff.
  set.seed(20261111)
  res <- ss_aipe_smd(true_smd = -0.5, width = 0.5, assurance = 0.9)
  n <- res$value[1]
  attained <- res$value[res$term == "assurance_attained"]
  cut <- .oracle_smd_cutoff(n, 0.5)
  G <- 1000
  d <- vapply(seq_len(G), function(g) {
    y1 <- stats::rnorm(n, -0.5)
    y2 <- stats::rnorm(n)
    (mean(y1) - mean(y2)) / sqrt((stats::var(y1) + stats::var(y2)) / 2)
  }, numeric(1))
  narrow <- vapply(d, function(x) .oracle_smd_width(x, n) <= 0.5, logical(1))
  expect_identical(narrow, abs(d) <= cut)
  expect_lt(abs(mean(narrow) - attained),
            4 * sqrt(attained * (1 - attained) / G))
})

test_that("ss_aipe_smd() offers only the analytic assurance method and records it", {
  # The default is "analytic", and with assurance the table reports the
  # exact probability, with the method recorded and named in the footer.
  expect_identical(formals(ss_aipe_smd)$assurance_method, "analytic")
  res <- ss_aipe_smd(true_smd = 0.5, width = 0.5, assurance = 0.9)
  expect_identical(res$term, c("necessary_n_per_group", "true_smd", "width",
                               "assurance", "assurance_attained"))
  expect_identical(res$value[1], 130)
  expect_gte(res$value[5], 0.9)
  expect_identical(attr(res, "assurance_method"), "analytic")
  expect_identical(attr(res, "assurance_probability"), "exact")
  expect_identical(ss_aipe_smd(true_smd = 0.5, width = 0.5, assurance = 0.9,
                               assurance_method = "analytic"), res)
  expect_match(paste(utils::capture.output(print(res)), collapse = " "),
               "Assurance method: analytic; assurance_attained is the exact",
               fixed = TRUE)
  # Without assurance there are no assurance rows and no record.
  plain <- ss_aipe_smd(true_smd = 0.5, width = 0.5)
  expect_null(attr(plain, "assurance_method"))
  expect_false("assurance_attained" %in% plain$term)
  # "simulation" is refused, with or without assurance, and so is anything
  # other than the two methods, partial matches included.
  for (a in list(0.9, NULL)) {
    e <- expect_error(ss_aipe_smd(true_smd = 0.5, width = 0.5, assurance = a,
                                  assurance_method = "simulation"),
                      class = "dmar_argument_invalid")
    expect_identical(e$arg, "assurance_method")
    expect_match(conditionMessage(e), paste0(
      "is not offered by ss_aipe_smd(): the analytic probability that the ",
      "interval is no wider than 'width' is exact for this interval"),
      fixed = TRUE)
    expect_match(conditionMessage(e), paste0(
      "Simulation is offered by ss_aipe_sc_ancova(), by ss_aipe_reg_coef() ",
      "with noncentral = TRUE, and by ss_aipe_src()"), fixed = TRUE)
  }
  # The two methods together are not a default here, since the formal is
  # "analytic", and the message names the one method offered.
  for (bad in list("anal", "sim", NA_character_, TRUE, "Analytic", "exact",
                   c("analytic", "simulation"))) {
    e <- expect_error(ss_aipe_smd(true_smd = 0.5, width = 0.5, assurance = 0.9,
                                  assurance_method = bad),
                      class = "dmar_argument_invalid")
    expect_identical(e$arg, "assurance_method")
    expect_match(conditionMessage(e), paste0(
      "'assurance_method' must be \"analytic\", the only method ",
      "ss_aipe_smd() offers;"), fixed = TRUE)
  }
})

test_that("ss_aipe_smd_sensitivity() plans through ss_aipe_smd(), for the magnitude of supposed_smd, without a warning", {
  sens <- function(...) ss_aipe_smd_sensitivity(
    true_smd = -.5, supposed_smd = -.5, width = .5, assurance = .9, G = 10,
    seed = 1, ...)
  expect_no_warning(out <- sens())
  expect_identical(out$value[out$term == "n_per_group"], 130)
  expect_no_warning(ss_aipe_smd_sensitivity(true_smd = .5, supposed_smd = .5, width = .5,
                                            assurance = .9, G = 10, seed = 1))
  expect_no_warning(ss_aipe_smd_sensitivity(true_smd = -.5, specified_n_per_group = 125,
                                            width = .5, assurance = .9, G = 10, seed = 1))
  # assurance_method is passed to the planner, and "simulation" is refused
  # in this function's name.
  expect_identical(formals(ss_aipe_smd_sensitivity)$assurance_method, "analytic")
  e <- tryCatch(sens(assurance_method = "simulation"), error = function(e) e)
  expect_s3_class(e, "dmar_argument_invalid")
  expect_identical(e$arg, "assurance_method")
  expect_match(conditionMessage(e), "is not offered by ss_aipe_smd_sensitivity()",
               fixed = TRUE)
  out <- sens(assurance_method = "analytic")
  expect_identical(out$value[out$term == "n_per_group"],
                   ss_aipe_smd(true_smd = -.5, width = .5, assurance = .9)$value[1])
})
