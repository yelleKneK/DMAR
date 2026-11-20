# The rule for bootstrap and Monte Carlo replications (R/dmar_bootstrap.R).
# These tests run with the rule as shipped: no local_small_B() here.

required <- function(...) DMAR:::.dmar_boot_required(...)
check_B  <- function(...) DMAR:::.dmar_check_B(...)
settle   <- function(...) DMAR:::.dmar_boot_settle(...)

# Collects the warnings of a call, muffled, with its value.
with_warnings <- function(expr) {
  w <- list()
  value <- withCallingHandlers(expr, warning = function(cnd) {
    w[[length(w) + 1L]] <<- cnd
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

test_that("the required counts are 10 beyond each limit, never fewer than 100", {
  expect_identical(required("quantile", (1 - 0.90) / 2), 200L)
  expect_identical(required("quantile", (1 - 0.95) / 2), 400L)
  expect_identical(required("quantile", (1 - 0.99) / 2), 2000L)
  expect_identical(required("quantile", (1 - 0.999) / 2), 20000L)
  for (tail in c(0.10, 0.16, 0.25, 0.5)) {
    expect_identical(required("quantile", tail), 100L)
  }
  # A one-sided bound at the 95th percentile.
  expect_identical(required("quantile", 0.05), 200L)
  # BCa: 50 beyond each nominal limit.
  expect_identical(required("bca", (1 - 0.90) / 2), 1000L)
  expect_identical(required("bca", (1 - 0.95) / 2), 2000L)
  expect_identical(required("bca", (1 - 0.99) / 2), 10000L)
  expect_identical(required("se"), 100L)
  expect_identical(required("p"), 100L)
  # A function's own floor is respected and does not lower the count.
  expect_identical(required("quantile", 0.025, floor = 1000L), 1000L)
  expect_identical(required("quantile", 0.005, floor = 1000L), 2000L)
  # A tail so small that no integer count serves.
  expect_identical(required("quantile", 1e-12), NA_integer_)
  expect_error(required("quantile", 0), "tail")
  expect_error(required("quantile", 0.6), "tail")
})

test_that("a B below the requirement is refused with the class, the minimum, and the reason", {
  err <- tryCatch(check_B(399, tail = 0.025, conf_level = 0.95),
                  error = function(e) e)
  expect_s3_class(err, "dmar_count_invalid")
  expect_s3_class(err, "dmar_refusal")
  expect_identical(err$min, 400L)
  expect_identical(err$arg, "B")
  expect_match(err$reason, "10 / 0.025 = 400", fixed = TRUE)
  expect_match(conditionMessage(err),
               "'B' must be a single whole number of at least 400. A 95% interval",
               fixed = TRUE)
  expect_match(conditionMessage(err), "see ?dmar_bootstrap.", fixed = TRUE)
  expect_identical(check_B(400, tail = 0.025, conf_level = 0.95), 400L)
  expect_identical(check_B(400L, tail = 0.025, conf_level = 0.95), 400L)
  # The other reasons.
  expect_error(check_B(50, tail = 0.16, conf_level = 0.68),
               "never from fewer than 100 replications",
               class = "dmar_count_invalid")
  expect_error(check_B(150, tail = 0.05, what = "bound", position = 0.95),
               "A bound at the 95th percentile", class = "dmar_count_invalid")
  expect_error(check_B(1999, kind = "bca", tail = 0.025, conf_level = 0.95),
               "50 / 0.025 = 2000", class = "dmar_count_invalid", fixed = TRUE)
  expect_error(check_B(99, kind = "se", conf_level = 0.95),
               "standard error taken from the replications",
               class = "dmar_count_invalid")
  expect_error(check_B(99, kind = "p"),
               "Monte Carlo p-value is reported only from at least 100",
               class = "dmar_count_invalid")
  # Not a count at all.
  for (bad in list(2.5, NA, "400", c(400, 500), Inf, NULL, -400)) {
    expect_error(check_B(bad, tail = 0.025, conf_level = 0.95),
                 class = "dmar_count_invalid")
  }
  # B = 0 where it means "no bootstrap".
  expect_identical(check_B(0, tail = 0.05, zero_ok = TRUE), 0L)
  expect_error(check_B(0, tail = 0.05), class = "dmar_count_invalid")
  expect_error(check_B(50, tail = 0.05, zero_ok = TRUE), "B = 0 runs none",
               class = "dmar_count_invalid")
})

test_that("the count check keeps its message and fields without a reason", {
  err <- tryCatch(DMAR:::.dmar_check_count(1, "B", min = 2L),
                  error = function(e) e)
  expect_identical(conditionMessage(err),
                   "'B' must be a single whole number of at least 2.")
  expect_null(err$reason)
  err <- tryCatch(DMAR:::.dmar_check_count(1, "B", min = 2L, reason = "Why."),
                  error = function(e) e)
  expect_identical(conditionMessage(err),
                   "'B' must be a single whole number of at least 2. Why.")
  expect_identical(err$reason, "Why.")
  expect_identical(err$min, 2L)
})

test_that("the exploratory option accepts a small B, marks the run, and still refuses a non-count", {
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  expect_identical(check_B(20, tail = 0.025, conf_level = 0.95), 20L)
  expect_error(check_B(1, tail = 0.025, conf_level = 0.95),
               class = "dmar_count_invalid")
  expect_error(check_B(2.5, tail = 0.025, conf_level = 0.95),
               class = "dmar_count_invalid")
  got <- with_warnings(settle(20, 20, tail = 0.025, conf_level = 0.95,
                              fun = "f", quantity = "the statistic",
                              dropped_because = "had no value"))
  expect_length(got$warnings, 1L)
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_exploratory")
  expect_identical(w$status, "exploratory")
  expect_identical(w$B_required, 400L)
  expect_identical(w$B, 20L)
  expect_match(conditionMessage(w),
               "f(): B = 20 is below the 400 replications a 95% interval",
               fixed = TRUE)
  expect_match(conditionMessage(w),
               "options(dmar.bootstrap_exploratory = TRUE)", fixed = TRUE)
  expect_true(got$value$formed)
  expect_identical(got$value$status, "exploratory")
  # Dropped replications are still reported, before the mark.
  got <- with_warnings(settle(20, 19, tail = 0.025, conf_level = 0.95,
                              fun = "f", quantity = "the statistic",
                              dropped_because = "had no value"))
  expect_identical(vapply(got$warnings, function(w) class(w)[1L], ""),
                   c("dmar_bootstrap_replications_dropped",
                     "dmar_bootstrap_exploratory"))
  # The share rule still applies.
  got <- with_warnings(settle(20, 9, tail = 0.025, conf_level = 0.95,
                              fun = "f", quantity = "the statistic",
                              dropped_because = "had no value"))
  expect_false(got$value$formed)
  expect_identical(got$value$status, "insufficient_share")
})

test_that("the option is off by default and only TRUE turns it on", {
  withr::local_options(dmar.bootstrap_exploratory = NULL)
  expect_false(DMAR:::.dmar_boot_exploratory())
  withr::local_options(dmar.bootstrap_exploratory = "yes")
  expect_false(DMAR:::.dmar_boot_exploratory())
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  expect_true(DMAR:::.dmar_boot_exploratory())
})

test_that("the verdict follows the share and the count rules at their boundaries", {
  s <- function(B, used, ...) {
    with_warnings(settle(B, used, tail = 0.025, conf_level = 0.95,
                         fun = "f", quantity = "kappa",
                         dropped_because = "left kappa undefined", ...))
  }
  # Nothing dropped: no condition.
  got <- s(1000, 1000)
  expect_length(got$warnings, 0L)
  expect_identical(got$value[c("formed", "status", "B", "B_used", "n_dropped",
                               "B_required", "caution")],
                   list(formed = TRUE, status = "complete", B = 1000L,
                        B_used = 1000L, n_dropped = 0L, B_required = 400L,
                        caution = FALSE))
  # One dropped: formed, dropped, no caution.
  got <- s(1000, 999)
  expect_length(got$warnings, 1L)
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_replications_dropped")
  expect_identical(w[c("B", "B_used", "n_dropped", "B_required", "caution",
                       "status")],
                   list(B = 1000L, B_used = 999L, n_dropped = 1L,
                        B_required = 400L, caution = FALSE,
                        status = "dropped"))
  expect_equal(w$share_dropped, 0.001)
  expect_match(conditionMessage(w),
               "f(): 1 of 1000 bootstrap replications (0.1%) left kappa undefined and were dropped.",
               fixed = TRUE)
  expect_no_match(conditionMessage(w), "Caution")
  # Exactly tail * B dropped: the caution (24 is below it, 25 reaches it).
  expect_false(s(1000, 976)$value$caution)
  got <- s(1000, 975)
  expect_true(got$value$caution)
  expect_true(got$warnings[[1L]]$caution)
  expect_match(conditionMessage(got$warnings[[1L]]),
               "Caution: the dropped replications are at least the 2.5%",
               fixed = TRUE)
  # Exactly half dropped is formed; one more is not.
  got <- s(1000, 500)
  expect_true(got$value$formed)
  got <- s(1000, 499)
  expect_false(got$value$formed)
  expect_identical(got$value$status, "insufficient_share")
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_insufficient")
  expect_match(conditionMessage(got$warnings[[1L]]),
               "More than half is too many for an interval at any B")
  # The count is held to the number requested, not to the number that
  # remained (ruled 2026-11-10): a request that meets the count and loses
  # replications still reports its interval, and the warning that reports
  # the drops says when those that remained are fewer than the number
  # required, and to request more.
  got <- s(700, 400)
  expect_true(got$value$formed)
  expect_match(conditionMessage(got$warnings[[1L]]),
               paste0("The 400 required are counted on the number requested ",
                      "(700), not on the number that remained; request more ",
                      "than 400 where bootstrap replications can fail."),
               fixed = TRUE)
  got <- s(700, 399)
  expect_true(got$value$formed)
  expect_identical(got$value$status, "dropped")
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_replications_dropped")
  expect_identical(w[c("B", "B_used", "n_dropped", "B_required")],
                   list(B = 700L, B_used = 399L, n_dropped = 301L,
                        B_required = 400L))
  expect_match(conditionMessage(w),
               "f(): 301 of 700 bootstrap replications (43%) left kappa undefined and were dropped.",
               fixed = TRUE)
  expect_match(conditionMessage(w),
               paste0("The 400 required are counted on the number requested ",
                      "(700), not on the number that remained; the 399 that ",
                      "remained are fewer, so request more bootstrap ",
                      "replications where some can fail, enough that at ",
                      "least 400 remain."),
               fixed = TRUE)
  # The case the ruling was made for: one replication dropped at exactly
  # the required B. Through the first versions of 1.1.0 no interval was
  # reported.
  got <- s(400, 399)
  expect_true(got$value$formed)
  expect_identical(got$value$status, "dropped")
  expect_false(got$value$caution)
  expect_length(got$warnings, 1L)
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_replications_dropped")
  # The share rules then decide: the caution at the tail area, no interval
  # past one half.
  expect_true(s(400, 390)$value$caution)
  expect_true(s(400, 200)$value$formed)
  expect_identical(s(400, 199)$value$status, "insufficient_share")
  # A request below the count that no check refused (an internal caller)
  # is not formed, whatever was dropped.
  got <- s(399, 399)
  expect_false(got$value$formed)
  expect_identical(got$value$status, "insufficient_count")
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_insufficient")
  expect_identical(w$B_required, 400L)
  expect_match(conditionMessage(w),
               paste0("f(): no interval is reported. B = 399 is below the ",
                      "400 that must be requested. A 95% interval read from ",
                      "the quantiles of the replications has 2.5% of them ",
                      "beyond each limit"),
               fixed = TRUE)
  expect_match(conditionMessage(w), "(10 / 0.025 = 400)", fixed = TRUE)
  expect_match(conditionMessage(s(399, 390)$warnings[[1L]]),
               "9 of 399 bootstrap replications left kappa undefined and were dropped, and 390 remained. B = 399 is below",
               fixed = TRUE)
  # The one statement of the rule.
  expect_true(DMAR:::.dmar_boot_count_met(400L, 400L))
  expect_false(DMAR:::.dmar_boot_count_met(399L, 400L))
  expect_true(DMAR:::.dmar_boot_count_met(5L, NA_integer_))
  # signal = FALSE returns the verdict in silence.
  expect_no_warning(v <- settle(399, 399, tail = 0.025, conf_level = 0.95,
                                fun = "f", quantity = "kappa",
                                dropped_because = "left kappa undefined",
                                signal = FALSE))
  expect_false(v$formed)
  # Nothing usable: no division by zero in the message.
  got <- with_warnings(settle(400, 0, kind = "se", conf_level = 0.95,
                              fun = "f", quantity = "kappa",
                              dropped_because = "left kappa undefined"))
  expect_identical(got$value$status, "insufficient_share")
  v <- settle(400, 0, tail = 0.025, conf_level = 0.95, fun = "f",
              quantity = "kappa", dropped_because = "left kappa undefined",
              signal = FALSE)
  expect_match(DMAR:::.dmar_boot_message(v, "insufficient_count"),
               "An interval needs at least two usable bootstrap replications.",
               fixed = TRUE)
  # The standard error and the p-value kinds: the same rule on the
  # request, and each states its own count when the request is below it.
  got <- with_warnings(settle(150, 99, kind = "se", conf_level = 0.95,
                              fun = "f", quantity = "the coefficient",
                              dropped_because = "returned no value"))
  expect_identical(got$value$status, "dropped")
  expect_match(conditionMessage(got$warnings[[1L]]),
               "The standard error is computed from the 99 that remained",
               fixed = TRUE)
  expect_match(conditionMessage(got$warnings[[1L]]),
               "the 99 that remained are fewer", fixed = TRUE)
  got <- with_warnings(settle(99, 99, kind = "se", conf_level = 0.95,
                              fun = "f", quantity = "the coefficient",
                              dropped_because = "returned no value"))
  expect_identical(got$value$status, "insufficient_count")
  expect_match(conditionMessage(got$warnings[[1L]]),
               "A standard error taken from the replications is reported only from a request of at least 100 of them.",
               fixed = TRUE)
  got <- with_warnings(settle(99, 99, kind = "p", conf_level = 0.95,
                              fun = "f", quantity = "the statistic",
                              dropped_because = "returned no value"))
  expect_identical(got$value$status, "insufficient_count")
  expect_match(conditionMessage(got$warnings[[1L]]),
               "A Monte Carlo p-value is reported only from at least 100 replications",
               fixed = TRUE)
  # In exploratory mode the mark follows the request, not what remained.
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  expect_identical(s(700, 399)$value$status, "dropped")
  expect_identical(s(399, 399)$value$status, "exploratory")
  expect_identical(s(3, 1)$value$status, "insufficient_share")
  expect_identical(s(2, 1)$value$status, "insufficient_count")
})

test_that("a BCa interval is withheld when fewer than 10 replications lie beyond an adjusted limit", {
  v <- settle(2000, 2000, kind = "bca", tail = 0.025, conf_level = 0.95,
              fun = "f", quantity = "the statistic",
              dropped_because = "had no value", signal = FALSE)
  expect_true(v$formed)
  # Positions that leave 10 or more beyond each limit: kept.
  kept <- DMAR:::.dmar_boot_bca_positions(v, c(0.005, 0.99))
  expect_true(kept$formed)
  # 0.4 percent of 2000 is 8 replications: withheld.
  held <- DMAR:::.dmar_boot_bca_positions(v, c(0.004, 0.99))
  expect_false(held$formed)
  expect_identical(held$status, "insufficient_count")
  expect_identical(held$B_required, 2500L)
  expect_identical(held$B_suggest, 3800L)
  expect_equal(held$positions, c(0.004, 0.99))
  got <- with_warnings(DMAR:::.dmar_boot_signal(held))
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_insufficient")
  expect_equal(w$positions, c(0.004, 0.99))
  expect_identical(w$B_suggest, 3800L)
  expect_match(conditionMessage(w), "f(): no BCa interval is reported.",
               fixed = TRUE)
  expect_match(conditionMessage(w),
               "at the 0.4% and 99% points of the bootstrap distribution, not at 2.5% and 97.5%",
               fixed = TRUE)
  expect_match(conditionMessage(w), "raise B to about 3800", fixed = TRUE)
  expect_match(conditionMessage(w),
               "which takes a B of 2500 at those points (2000 were requested)",
               fixed = TRUE)
  # The count at the positions is held to the number requested, as the
  # count of the request is: 2600 requested and 2400 usable is kept.
  v2 <- settle(2600, 2400, kind = "bca", tail = 0.025, conf_level = 0.95,
               fun = "f", quantity = "the statistic",
               dropped_because = "had no value", signal = FALSE)
  expect_true(DMAR:::.dmar_boot_bca_positions(v2, c(0.004, 0.99))$formed)
  # A limit beyond every replication.
  gone <- DMAR:::.dmar_boot_bca_positions(v, c(0, 0.99))
  expect_false(gone$formed)
  expect_match(DMAR:::.dmar_boot_message(gone),
               "place a limit beyond every replication", fixed = TRUE)
  # Under the exploratory option the interval is not withheld.
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  expect_true(DMAR:::.dmar_boot_bca_positions(v, c(0.004, 0.99))$formed)
})

test_that("the footer line states the counts, the reason for no interval, and the exploratory mark", {
  foot <- DMAR:::.dmar_boot_footer
  x <- structure(list(), B = 2000L, B_used = 1988L, n_dropped = 12L,
                 B_required = 400L, bootstrap_status = "dropped")
  expect_identical(foot(x), paste0("Bootstrap: 2000 replications, 1988 used ",
                                   "(12 dropped; at least 400 must be ",
                                   "requested)"))
  attr(x, "bootstrap_status") <- "insufficient_share"
  expect_match(foot(x), "No interval: more than half", fixed = TRUE)
  attr(x, "bootstrap_status") <- "insufficient_count"
  expect_match(foot(x), "No interval: fewer than the 400 required were requested.",
               fixed = TRUE)
  attr(x, "B_used") <- 1L
  expect_match(foot(x), "No interval: fewer than two replications were usable.",
               fixed = TRUE)
  attr(x, "bootstrap_status") <- "exploratory"
  expect_match(foot(x),
               "Exploratory: below the number required for a reported interval.",
               fixed = TRUE)
  expect_null(foot(structure(list(), B = 10L)))
})

# ---- the rule through the public functions --------------------------------

test_that("mediate() refuses a small B, runs it under the option, and records the counts", {
  set.seed(113)
  n <- 80
  x <- rnorm(n); m <- 0.5 * x + rnorm(n); y <- 0.4 * m + rnorm(n)
  d <- data.frame(x = x, m = m, y = y)
  err <- tryCatch(mediate(d, "x", "m", "y", B = 399), error = function(e) e)
  expect_s3_class(err, "dmar_count_invalid")
  expect_identical(err$min, 400L)
  expect_error(mediate(d, "x", "m", "y", ci_method = "boot_bca", B = 1999),
               "at least 2000", class = "dmar_count_invalid")
  expect_error(mediate(d, "x", "m", "y", conf_level = 0.99, B = 1999),
               "at least 2000", class = "dmar_count_invalid")
  expect_error(mediate(d, "x", "m", "y", ci_method = "monte_carlo", B = 100),
               "at least 400", class = "dmar_count_invalid")
  res <- mediate(d, "x", "m", "y", B = 400, seed = 113)
  expect_identical(attr(res, "B"), 400L)
  expect_identical(attr(res, "B_used"), 400L)
  expect_identical(attr(res, "n_dropped"), 0L)
  expect_identical(attr(res, "B_required"), 400L)
  expect_identical(attr(res, "bootstrap_status"), "complete")
  meta <- dmar_metadata(res)$bootstrap
  expect_identical(meta$required, 400L)
  expect_identical(meta$status, "complete")
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  got <- with_warnings(mediate(d, "x", "m", "y", B = 20, seed = 113))
  expect_length(got$warnings, 1L)
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_exploratory")
  expect_identical(attr(got$value, "bootstrap_status"), "exploratory")
  expect_identical(attr(got$value, "B_required"), 400L)
  expect_true(all(is.finite(got$value$value[got$value$term %in%
                                            c("indirect_effect_ci_lower",
                                              "indirect_effect_ci_upper")])))
})

test_that("an interval that is still reported has the limits DMAR 1.0.x read from the same draws", {
  # The percentile limits are the type 7 quantiles of the replicated
  # indirect effects, as before the rule.
  set.seed(7)
  n <- 60
  x <- rnorm(n); m <- 0.5 * x + rnorm(n); y <- 0.4 * m + rnorm(n)
  d <- data.frame(x = x, m = m, y = y)
  res <- mediate(d, "x", "m", "y", B = 400, seed = 31)
  ref <- withr::with_seed(31, vapply(seq_len(400), function(i) {
    idx <- sample.int(n, n, replace = TRUE)
    dd <- d[idx, ]
    unname(coef(lm(m ~ x, dd))["x"] * coef(lm(y ~ x + m, dd))["m"])
  }, numeric(1)))
  expect_equal(res$value[res$term %in% c("indirect_effect_ci_lower", "indirect_effect_ci_upper")],
               quantile(ref, c(0.025, 0.975), names = FALSE),
               tolerance = 1e-10)
})

# Two raters whose resamples sometimes fall in one category, so that kappa
# is undefined there (chance agreement 1).
.rare_kappa <- function() {
  list(r1 = c("b", "a", rep("a", 6)), r2 = c("b", "b", rep("a", 6)))
}

test_that("cohen_kappa() drops undefined resamples, withholds when too few remain, and refuses a fractional B", {
  k <- .rare_kappa()
  expect_error(cohen_kappa(k$r1, k$r2, ci_method = "percentile", B = 2.5),
               class = "dmar_count_invalid")
  expect_error(cohen_kappa(k$r1, k$r2, ci_method = "percentile", B = 399),
               "at least 400", class = "dmar_count_invalid")
  expect_error(cohen_kappa(k$r1, k$r2, ci_method = "bca", B = 400),
               "at least 2000", class = "dmar_count_invalid")
  # B is not read for the Wald interval.
  expect_no_error(cohen_kappa(k$r1, k$r2, B = 2.5))
  got <- with_warnings(cohen_kappa(k$r1, k$r2, ci_method = "percentile",
                                   B = 2000, seed = 113))
  res <- got$value
  expect_length(got$warnings, 1L)
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_replications_dropped")
  expect_identical(w$B, 2000L)
  expect_identical(w$B_used + w$n_dropped, 2000L)
  expect_gt(w$n_dropped, 0L)
  expect_identical(attr(res, "B"), 2000L)
  expect_identical(attr(res, "B_used"), w$B_used)
  expect_identical(attr(res, "n_dropped"), w$n_dropped)
  expect_identical(attr(res, "B_required"), 400L)
  expect_identical(attr(res, "bootstrap_status"), "dropped")
  # The limits are the quantiles of the resamples in which kappa is defined.
  ref <- withr::with_seed(113, vapply(seq_len(2000), function(b) {
    idx <- sample.int(8, 8, replace = TRUE)
    tab <- table(factor(k$r1[idx], c("a", "b")), factor(k$r2[idx], c("a", "b"))) / 8
    pe <- sum(rowSums(tab) * colSums(tab))
    if (1 - pe <= .Machine$double.eps) NA_real_ else (sum(diag(tab)) - pe) / (1 - pe)
  }, numeric(1)))
  expect_identical(sum(is.finite(ref)), w$B_used)
  expect_equal(c(res$ci_lower, res$ci_upper),
               quantile(ref[is.finite(ref)], c(0.025, 0.975), names = FALSE),
               tolerance = 1e-10)
  # At exactly the required B with resamples dropped the interval is
  # reported (the count is held to the request), and the warning says
  # that fewer than the number required remained.
  got <- with_warnings(cohen_kappa(k$r1, k$r2, ci_method = "percentile",
                                   B = 400, seed = 113))
  expect_length(got$warnings, 1L)
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_replications_dropped")
  expect_identical(w$status, "dropped")
  expect_lt(w$B_used, 400L)
  expect_gte(w$B_used, 200L)
  expect_match(conditionMessage(w), "^cohen_kappa\\(\\): ")
  expect_match(conditionMessage(w),
               "request more bootstrap replications where some can fail",
               fixed = TRUE)
  expect_true(is.finite(got$value$ci_lower) &&
                is.finite(got$value$ci_upper))
  expect_identical(attr(got$value, "bootstrap_status"), "dropped")
  # A request that is settled as below the count: NA limits, the estimate
  # untouched, a classed warning.
  local_count_unmet()
  got <- with_warnings(cohen_kappa(k$r1, k$r2, ci_method = "percentile",
                                   B = 400, seed = 113))
  expect_length(got$warnings, 1L)
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_insufficient")
  expect_identical(got$warnings[[1L]]$status, "insufficient_count")
  expect_true(is.na(got$value$ci_lower) && is.na(got$value$ci_upper))
  expect_equal(got$value$kappa, cohen_kappa(k$r1, k$r2)$kappa)
  expect_equal(got$value$se, cohen_kappa(k$r1, k$r2)$se)
  expect_identical(attr(got$value, "bootstrap_status"), "insufficient_count")
})

test_that("the kappa functions refuse a BCa interval whose jackknife is undefined", {
  k <- .rare_kappa()
  # Leaving out case 1 leaves rater 1 with one category and rater 2 with
  # two, so kappa is defined; a table whose leave-one-out is degenerate:
  r1 <- c("b", rep("a", 7)); r2 <- c("b", rep("a", 7))
  err <- tryCatch(
    suppressWarnings(cohen_kappa(r1, r2, ci_method = "bca", B = 4000,
                                 seed = 113)),
    error = function(e) e, warning = function(w) w)
  # Either the BCa interval is refused for the undefined jackknife value,
  # or so many resamples were dropped that no interval is formed before
  # the jackknife is reached; both are the rule, never a silent drop.
  expect_true(inherits(err, "dmar_bca_undefined") ||
                identical(attr(err, "bootstrap_status"),
                          "insufficient_share") ||
                identical(attr(err, "bootstrap_status"),
                          "insufficient_count"))
  if (inherits(err, "dmar_bca_undefined")) {
    expect_s3_class(err, "dmar_refusal")
    expect_match(conditionMessage(err), "with case 1 left out")
    expect_match(conditionMessage(err),
                 "ci_method = \"percentile\" does not need these quantities.",
                 fixed = TRUE)
  }
  # fleiss_kappa(): two raters, all subjects but one in category 1.
  ratings <- cbind(c(0, rep(2, 7)), c(2, rep(0, 7)))
  err <- tryCatch(
    suppressWarnings(fleiss_kappa(ratings, ci_method = "bca", B = 4000,
                                  seed = 113)),
    error = function(e) e)
  expect_true(inherits(err, "dmar_bca_undefined") ||
                startsWith(attr(err, "bootstrap_status") %||% "", "insufficient"))
})

test_that("fleiss_kappa() and krippendorff_alpha() follow the rule", {
  ratings <- cbind(c(0, 1, rep(2, 6)), c(2, 1, rep(0, 6)))
  expect_error(fleiss_kappa(ratings, ci_method = "percentile", B = 2.5),
               class = "dmar_count_invalid")
  expect_error(fleiss_kappa(ratings, ci_method = "percentile", B = 399),
               "at least 400", class = "dmar_count_invalid")
  got <- with_warnings(fleiss_kappa(ratings, ci_method = "percentile",
                                    B = 2000, seed = 113))
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_replications_dropped")
  expect_identical(attr(got$value, "B"), 2000L)
  expect_identical(attr(got$value, "B_used") + attr(got$value, "n_dropped"),
                   2000L)
  expect_identical(attr(got$value, "B_required"), 400L)
  # At exactly the required B the interval is reported with the drops.
  got <- with_warnings(fleiss_kappa(ratings, ci_method = "percentile",
                                    B = 400, seed = 113))
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_replications_dropped")
  expect_lt(attr(got$value, "B_used"), 400L)
  expect_true(is.finite(got$value$ci_lower) &&
                is.finite(got$value$ci_upper))
  expect_identical(attr(got$value, "bootstrap_status"), "dropped")
  local({
    local_count_unmet()
    got <- with_warnings(fleiss_kappa(ratings, ci_method = "percentile",
                                      B = 400, seed = 113))
    expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_insufficient")
    expect_true(is.na(got$value$ci_lower) && is.na(got$value$ci_upper))
    expect_equal(got$value$kappa, fleiss_kappa(ratings)$kappa)
  })

  # Krippendorff's alpha: five units, four alike, so many resamples have
  # values all alike and alpha is undefined there.
  R <- cbind(c(1, 1, 1, 1, 2), c(1, 1, 1, 1, 2))
  expect_error(krippendorff_alpha(R, boot = TRUE, B = 2.5),
               class = "dmar_count_invalid")
  expect_error(krippendorff_alpha(R, boot = TRUE, B = 100), "at least 400",
               class = "dmar_count_invalid")
  expect_no_error(krippendorff_alpha(R, B = 2.5))  # B is not read
  got <- with_warnings(krippendorff_alpha(R, boot = TRUE, B = 2000,
                                          seed = 113))
  res <- got$value
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_replications_dropped")
  expect_identical(attr(res, "B"), 2000L)
  expect_equal(res$value[res$term == "B_used"], attr(res, "B_used"))
  expect_identical(attr(res, "B_required"), 400L)
  expect_true(all(is.finite(res$value[res$term %in%
                                      c("ci_lower", "ci_upper")])))
  got <- with_warnings(krippendorff_alpha(R, boot = TRUE, B = 400,
                                          seed = 113))
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_replications_dropped")
  res <- got$value
  expect_lt(attr(res, "B_used"), 400L)
  expect_true(all(is.finite(res$value[res$term %in%
                                      c("ci_lower", "ci_upper")])))
  local_count_unmet()
  got <- with_warnings(krippendorff_alpha(R, boot = TRUE, B = 400,
                                          seed = 113))
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_insufficient")
  res <- got$value
  expect_true(all(is.na(res$value[res$term %in%
                                  c("ci_lower", "ci_upper")])))
  expect_equal(res$value[res$term == "krippendorff_alpha"],
               krippendorff_alpha(R)$value[1L])
  expect_null(attr(res, "conf_level"))
})

test_that("htmt() needs 200 resamples for a 95 percent bound, reports it when some are dropped, and withholds it when the count is not met", {
  set.seed(113)
  n <- 12
  f <- rnorm(n)
  d <- data.frame(a1 = f + rnorm(n), a2 = f + rnorm(n),
                  b1 = f + rnorm(n), b2 = c(1, rep(0, n - 1)))
  blocks <- list(A = c("a1", "a2"), B = c("b1", "b2"))
  expect_error(htmt(d, blocks, B = 199), "at least 200",
               class = "dmar_count_invalid")
  expect_error(htmt(d, blocks, B = 199), "B = 0 runs none",
               class = "dmar_count_invalid")
  expect_error(htmt(d, blocks, B = 2.5), class = "dmar_count_invalid")
  expect_error(htmt(d, blocks, B = 999, conf_level = 0.99), "at least 1000",
               class = "dmar_count_invalid")
  expect_null(attr(suppressWarnings(htmt(d, blocks)), "B"))
  # Item b2 is constant in the resamples that miss case 1 (about a third).
  got <- with_warnings(htmt(d, blocks, B = 1000, seed = 113))
  res <- got$value
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_true("dmar_bootstrap_replications_dropped" %in% cls)
  w <- got$warnings[[which(cls == "dmar_bootstrap_replications_dropped")]]
  expect_identical(w$B, 1000L)
  expect_identical(w$B_required, 200L)
  expect_true(w$caution)
  expect_identical(attr(res, "B_used"), w$B_used)
  expect_identical(attr(res, "bootstrap_status"), "dropped")
  expect_true(is.finite(res$ci_upper))
  got <- with_warnings(htmt(d, blocks, B = 200, seed = 113))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_true("dmar_bootstrap_replications_dropped" %in% cls)
  expect_false("dmar_bootstrap_insufficient" %in% cls)
  expect_lt(attr(got$value, "B_used"), 200L)
  expect_true(is.finite(got$value$ci_upper))
  local_count_unmet()
  got <- with_warnings(htmt(d, blocks, B = 200, seed = 113))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_true("dmar_bootstrap_insufficient" %in% cls)
  expect_true(is.na(got$value$ci_upper))
  expect_equal(got$value$htmt, suppressWarnings(htmt(d, blocks))$htmt)
})

test_that("the ratio functions follow the rule at their tail areas", {
  set.seed(113)
  x <- rpois(12, 0.4) + 0; y <- rpois(12, 0.25) + 0
  y[1] <- 1; y[-1] <- 0  # the denominator mean is zero unless case 1 is drawn
  expect_error(ci_ratio_of_means(x, y, ci_method = "percentile", B = 399),
               "at least 400", class = "dmar_count_invalid")
  expect_error(ci_ratio_of_means(x, y, ci_method = "bca", B = 1999),
               "at least 2000", class = "dmar_count_invalid")
  expect_error(ci_ratio_of_means(x, y, ci_method = "percentile", B = 2.5),
               class = "dmar_count_invalid")
  # Unequal tail areas: the smaller one decides (10 / 0.01 = 1000).
  expect_error(ci_ratio_of_means(x, y, ci_method = "percentile", B = 999,
                                 alpha_lower = 0.01, alpha_upper = 0.04),
               "at least 1000", class = "dmar_count_invalid")
  # A one-sided bound: 10 beyond its one limit (10 / 0.05 = 200).
  expect_error(ci_ratio_of_means(x, y, ci_method = "percentile", B = 199,
                                 alpha_lower = 0, alpha_upper = 0.05),
               "at least 200", class = "dmar_count_invalid")
  got <- with_warnings(ci_ratio_of_means(x, y, ci_method = "percentile",
                                         B = 2000, seed = 113))
  res <- got$value
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_true("dmar_bootstrap_replications_dropped" %in% cls)
  expect_identical(attr(res, "B"), 2000L)
  expect_identical(attr(res, "B_used") + attr(res, "n_dropped"), 2000L)
  expect_identical(attr(res, "B_required"), 400L)
  expect_output(print(res), "at least 400 must be requested")
  expect_output(print(res), "denominator of zero")
  got <- with_warnings(ci_ratio_of_means(x, y, ci_method = "percentile",
                                         B = 400, seed = 113))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_true("dmar_bootstrap_replications_dropped" %in% cls)
  expect_false("dmar_bootstrap_insufficient" %in% cls)
  res <- got$value
  expect_lt(attr(res, "B_used"), 400L)
  expect_true(all(is.finite(res$value[res$term %in%
                                      c("ci_lower", "ci_upper")])))
  local({
    local_count_unmet()
    got <<- with_warnings(ci_ratio_of_means(x, y, ci_method = "percentile",
                                            B = 400, seed = 113))
  })
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_true("dmar_bootstrap_insufficient" %in% cls)
  res <- got$value
  expect_true(all(is.na(res$value[res$term %in%
                                  c("ci_lower", "ci_upper")])))
  expect_equal(res$value[res$term == "ratio"], mean(x) / mean(y))
  expect_output(print(res), "No interval")
  # Paired.
  expect_error(ci_ratio_of_means_paired(x, y, ci_method = "percentile",
                                        B = 399),
               "at least 400", class = "dmar_count_invalid")
  got <- with_warnings(ci_ratio_of_means_paired(x, y, ci_method = "percentile",
                                                B = 2000, seed = 113))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_true("dmar_bootstrap_replications_dropped" %in% cls)
  expect_identical(attr(got$value, "B_required"), 400L)
})

test_that("the reliability family follows the rule by interval method", {
  set.seed(113)
  X <- matrix(rnorm(60 * 4), 60, 4) + rnorm(60)
  colnames(X) <- paste0("x", 1:4)
  expect_error(reliability_alpha(X, ci_method = "percentile", B = 399),
               "at least 400", class = "dmar_count_invalid")
  expect_error(reliability_alpha(X, ci_method = "bca", B = 1999),
               "at least 2000", class = "dmar_count_invalid")
  expect_error(reliability_alpha(X, ci_method = "bootstrap_se", B = 99),
               "at least 100", class = "dmar_count_invalid")
  expect_error(reliability_alpha(X, ci_method = "percentile", B = 2.5),
               class = "dmar_count_invalid")
  expect_error(reliability_alpha(X, ci_method = "percentile", B = 1000,
                                 conf_level = 0.99),
               "at least 2000", class = "dmar_count_invalid")
  expect_error(reliability_kr20(1 * (X > 0), ci_method = "percentile",
                                B = 399),
               "at least 400", class = "dmar_count_invalid")
  expect_error(reliability_omega(X, ci_method = "percentile", B = 399),
               "at least 400", class = "dmar_count_invalid")
  # B is not read when no bootstrap is asked for.
  expect_no_error(reliability_alpha(X, B = 2.5))
  a <- reliability_alpha(X, ci_method = "bootstrap_se", B = 100, seed = 113)
  expect_identical(attr(a, "B_required"), 100L)
  expect_identical(attr(a, "bootstrap_status"), "complete")
  expect_identical(attr(a, "n_dropped"), 0L)
  bs <- attr(a, "bootstrap")
  expect_identical(bs$required, 100L)
  expect_identical(bs$bootstrap_status, "complete")
  expect_true(bs$require_half)
  expect_identical(dmar_metadata(a)$bootstrap$required, 100L)
  expect_identical(dmar_metadata(a)$bootstrap$status, "complete")
  # The exploratory option: one classed warning, the mark in the footer.
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  got <- with_warnings(reliability_alpha(X, ci_method = "percentile", B = 40,
                                         seed = 113))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_identical(cls, "dmar_bootstrap_exploratory")
  expect_identical(attr(got$value, "bootstrap_status"), "exploratory")
  expect_identical(attr(got$value, "B_required"), 400L)
  expect_output(print(got$value), "Exploratory: below the number required")
})

test_that("the engine holds the count to the request and withholds an interval when fewer than half return a value, with the class", {
  skip_if_not_installed("boot")
  set.seed(113)
  X <- as.data.frame(matrix(rnorm(40 * 3), 40, 3))
  # A statistic that returns no value in about 30 percent of resamples.
  point_fn <- function(d) {
    if (d[1L, 1L] > stats::quantile(X[[1L]], 0.7)) stop("no value")
    mean(d[[1L]])
  }
  # 500 requested, fewer than the 400 required return a value: the
  # interval is reported (the count is held to the request), and the
  # warning says that fewer than the number required returned a value.
  got <- with_warnings(DMAR:::.bootstrap_ci(
    X, point_fn, B = 500, conf_level = 0.95, kind = "percentile", seed = 113))
  r <- got$value
  expect_length(got$warnings, 1L)
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_replications_dropped")
  expect_identical(w$status, "dropped")
  expect_identical(w$B_required, 400L)
  expect_identical(w$B, 500L)
  expect_lt(w$B_used, 400L)
  expect_gte(w$B_used, 250L)
  expect_match(conditionMessage(w),
               paste0("the 400 required are counted on the number requested ",
                      "(500), not on the number that remained; the ",
                      w$B_used, " that remained are fewer, so request more ",
                      "replications where some can fail"), fixed = TRUE)
  expect_true(is.finite(r$lower) && is.finite(r$upper) && is.finite(r$se))
  expect_identical(r$bootstrap$bootstrap_status, "dropped")
  expect_true(w$caution)
  # With more replications those that return a value reach the count and
  # the note is the short one.
  got <- with_warnings(DMAR:::.bootstrap_ci(
    X, point_fn, B = 800, conf_level = 0.95, kind = "percentile", seed = 113))
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_replications_dropped")
  expect_match(conditionMessage(w),
               "not on the number that remained; request more than 400 where replications can fail",
               fixed = TRUE)
  expect_true(is.finite(got$value$lower))
  expect_identical(got$value$bootstrap$bootstrap_status, "dropped")
  expect_true(w$caution)
  # Fewer than half return a value: no interval, by the share rule.
  rare_fn <- function(d) {
    if (d[1L, 1L] > stats::quantile(X[[1L]], 0.3)) stop("no value")
    mean(d[[1L]])
  }
  got <- with_warnings(DMAR:::.bootstrap_ci(
    X, rare_fn, B = 500, conf_level = 0.95, kind = "percentile", seed = 113))
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_insufficient")
  expect_identical(w$status, "insufficient_share")
  expect_match(conditionMessage(w), "fewer than half returned a value")
  expect_true(is.na(got$value$lower) && is.na(got$value$upper) &&
                is.na(got$value$se))
  # A request below the count that no check refused: no interval, and the
  # message says that the request is what fell short.
  got <- with_warnings(DMAR:::.bootstrap_ci(
    X, point_fn, B = 300, conf_level = 0.95, kind = "percentile", seed = 113))
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_insufficient")
  expect_identical(w$status, "insufficient_count")
  expect_match(conditionMessage(w),
               "the 300 requested are fewer than the 400 required at this level",
               fixed = TRUE)
  expect_true(is.na(got$value$lower))
})

test_that("R2_mixed_effects() checks B, drops a replicate for both quantities, and records the counts", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  fit <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
  expect_error(R2_mixed_effects(fit, ci_method = "boot", B = 399),
               "at least 400", class = "dmar_count_invalid")
  expect_error(R2_mixed_effects(fit, ci_method = "boot", B = 2.5),
               class = "dmar_count_invalid")
  expect_no_error(R2_mixed_effects(fit, B = 2.5))  # B is not read
  # bootMer() replaced by a function that returns fixed replicates, with a
  # missing value for one quantity in some of them.
  t_mat <- cbind(seq(0.10, 0.30, length.out = 500),
                 seq(0.60, 0.80, length.out = 500))
  t_mat[c(3, 40, 77), 1L] <- NA
  t_mat[c(120, 121), 2L] <- NaN
  testthat::local_mocked_bindings(
    bootMer = function(x, FUN, nsim, ...) list(t = t_mat[seq_len(nsim), ]),
    .package = "lme4")
  got <- with_warnings(R2_mixed_effects(fit, ci_method = "boot", B = 500))
  res <- got$value
  expect_length(got$warnings, 1L)
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_replications_dropped")
  expect_identical(w$n_dropped, 5L)
  expect_identical(attr(res, "B"), 500L)
  expect_identical(attr(res, "B_used"), 495L)
  expect_identical(attr(res, "n_dropped"), 5L)
  expect_identical(attr(res, "B_required"), 400L)
  keep <- stats::complete.cases(t_mat)
  val <- function(term) res$value[res$term == term]
  expect_equal(c(val("R2_marginal_ci_lower"), val("R2_marginal_ci_upper")),
               quantile(t_mat[keep, 1L], c(0.025, 0.975), names = FALSE))
  expect_equal(c(val("R2_conditional_ci_lower"), val("R2_conditional_ci_upper")),
               quantile(t_mat[keep, 2L], c(0.025, 0.975), names = FALSE))
  # DMAR 1.0.x removed a missing value for one quantity only; the limits
  # of the other quantity then rested on other replicates.
  old_marginal <- quantile(t_mat[, 1L], c(0.025, 0.975), na.rm = TRUE,
                           names = FALSE)
  expect_false(isTRUE(all.equal(c(val("R2_marginal_ci_lower"),
                                  val("R2_marginal_ci_upper")), old_marginal)))
  # Fewer than the 400 required remain of the 500 requested: the limits
  # are reported, since the count is held to the request.
  t_mat[1:150, 1L] <- NA
  got <- with_warnings(R2_mixed_effects(fit, ci_method = "boot", B = 500))
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_replications_dropped")
  expect_lt(attr(got$value, "B_used"), 400L)
  expect_match(conditionMessage(got$warnings[[1L]]),
               "that remained are fewer, so request more")
  expect_true(all(is.finite(got$value$value[grepl("lower|upper",
                                                   got$value$term)])))
  # More than half dropped: NA limits and the estimates untouched.
  t_mat[1:260, 1L] <- NA
  got <- with_warnings(R2_mixed_effects(fit, ci_method = "boot", B = 500))
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_insufficient")
  expect_identical(attr(got$value, "bootstrap_status"), "insufficient_share")
  res <- got$value
  expect_true(all(is.na(res$value[grepl("lower|upper", res$term)])))
  expect_equal(res$value[res$term == "R2_marginal"],
               R2_mixed_effects(fit)$value[1L])
})

test_that("mlmr() and mlmr_mv() check B and drop a whole resample that returns no value", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford
  f <- t6_paragraph_comprehension ~ t5_general_information + t9_word_meaning
  expect_error(mlmr(f, data = d, ci_method = "boot", B = 399),
               "at least 400", class = "dmar_count_invalid")
  expect_error(mlmr(f, data = d, ci_method = "boot", B = 2.5),
               class = "dmar_count_invalid")
  expect_error(mlmr_mv(cbind(t6_paragraph_comprehension, t9_word_meaning) ~
                         t5_general_information, data = d,
                       ci_method = "boot", B = 399),
               "at least 400", class = "dmar_count_invalid")
  # B is not read for the other intervals.
  expect_no_error(mlmr(f, data = d, ci_method = "wald", B = 2.5,
                       effect_sizes = FALSE))
  # lavBootstrap() replaced by fixed replicates with missing coefficients.
  set.seed(113)
  boot_mat <- matrix(rnorm(500 * 3), 500, 3)
  boot_mat[c(5, 50), 2L] <- NA
  boot_mat[400, 3L] <- NA
  testthat::local_mocked_bindings(
    lavBootstrap = function(object, R = 1000L, r = R, ...) {
      boot_mat[seq_len(r), , drop = FALSE]
    }, .package = "lavaan")
  got <- with_warnings(mlmr(f, data = d, ci_method = "boot", B = 500,
                            effect_sizes = FALSE))
  fit <- got$value
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_true("dmar_bootstrap_replications_dropped" %in% cls)
  w <- got$warnings[[which(cls == "dmar_bootstrap_replications_dropped")]]
  expect_identical(w$n_dropped, 3L)
  expect_identical(fit$B, 500L)
  expect_identical(fit$B_used, 497L)
  expect_identical(fit$n_dropped, 3L)
  expect_identical(fit$B_required, 400L)
  expect_identical(fit$bootstrap_status, "dropped")
  keep <- stats::complete.cases(boot_mat)
  expect_equal(unname(fit$ci[, "lower"]),
               apply(boot_mat[keep, ], 2L, quantile, probs = 0.025,
                     names = FALSE))
  expect_equal(unname(fit$ci[, "upper"]),
               apply(boot_mat[keep, ], 2L, quantile, probs = 0.975,
                     names = FALSE))
  # DMAR 1.0.x: the first coefficient's limits used all 500 resamples.
  expect_false(isTRUE(all.equal(
    unname(fit$ci[1L, ]),
    quantile(boot_mat[, 1L], c(0.025, 0.975), names = FALSE))))
  expect_output(print(fit), "Bootstrap: 500 replications, 497 used")
  expect_output(print(summary(fit)), "at least 400 must be requested")
  # Fewer than the 400 required remain of the 500 requested: the limits
  # are reported, since the count is held to the request.
  boot_mat[1:120, 1L] <- NA
  got <- with_warnings(mlmr(f, data = d, ci_method = "boot", B = 500,
                            effect_sizes = FALSE))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_true("dmar_bootstrap_replications_dropped" %in% cls)
  expect_false("dmar_bootstrap_insufficient" %in% cls)
  expect_lt(got$value$B_used, 400L)
  expect_true(all(is.finite(got$value$ci)))
  expect_identical(got$value$bootstrap_status, "dropped")
  # More than half dropped: NA limits, coefficients untouched.
  boot_mat[1:260, 1L] <- NA
  got <- with_warnings(mlmr(f, data = d, ci_method = "boot", B = 500,
                            effect_sizes = FALSE))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_true("dmar_bootstrap_insufficient" %in% cls)
  expect_true(all(is.na(got$value$ci)))
  expect_equal(coef(got$value),
               coef(mlmr(f, data = d, ci_method = "wald",
                         effect_sizes = FALSE)))
  expect_identical(got$value$bootstrap_status, "insufficient_share")
})

test_that("parallel_analysis() and factor_count() require the reference data sets the percentile needs", {
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  d <- holzinger_swineford[items]
  expect_error(parallel_analysis(d, B = 199), "at least 200",
               class = "dmar_count_invalid")
  expect_error(parallel_analysis(d, B = 199),
               "A reference quantile at the 95th percentile",
               class = "dmar_count_invalid")
  expect_error(parallel_analysis(d, B = 2.5), class = "dmar_count_invalid")
  expect_error(parallel_analysis(d, B = 99, percentile = 0.5),
               "at least 100", class = "dmar_count_invalid")
  expect_error(parallel_analysis(d, B = 999, percentile = 0.99),
               "at least 1000", class = "dmar_count_invalid")
  pa <- parallel_analysis(d, B = 200, seed = 113)
  expect_identical(attr(pa, "B"), 200L)
  expect_identical(attr(pa, "B_used"), 200L)
  expect_identical(attr(pa, "n_dropped"), 0L)
  expect_identical(attr(pa, "B_required"), 200L)
  expect_identical(attr(pa, "bootstrap_status"), "complete")
  expect_identical(attr(pa, "counts")$usable, 200L)
  expect_error(factor_count(d, B = 199), "at least 200",
               class = "dmar_count_invalid")
  expect_error(factor_count(d, B = 2.5), class = "dmar_count_invalid")
  expect_error(factor_count(d, B = -1), class = "dmar_count_invalid")
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  got <- with_warnings(parallel_analysis(d, B = 20, seed = 113))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_identical(cls, "dmar_bootstrap_exploratory")
  expect_identical(attr(got$value, "bootstrap_status"), "exploratory")
  expect_output(print(got$value), "Exploratory: 20 reference data sets")
})

test_that("factor_count() carries the counts of its parallel analysis", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  fc <- factor_count(holzinger_swineford, items = items, B = 200, seed = 113)
  expect_identical(attr(fc, "B"), 200L)
  expect_identical(attr(fc, "B_required"), 200L)
  expect_identical(attr(fc, "bootstrap_status"), "complete")
  none <- factor_count(holzinger_swineford, items = items, B = 0)
  expect_null(attr(none, "B"))
  expect_null(attr(none, "bootstrap_status"))
})

test_that("the functions that refit a lavaan model refuse a small or fractional B before refitting", {
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa(
    "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
     textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    data = holzinger_swineford, std.lv = TRUE)
  for (B in list(2.5, 399)) {
    expect_error(average_variance_extracted(fit, ci_method = "percentile",
                                            B = B),
                 class = "dmar_count_invalid")
    expect_error(fraction_missing_information(fit, ci_method = "percentile",
                                              B = B),
                 class = "dmar_count_invalid")
  }
  expect_error(sem_exact_fit_bootstrap(fit, B = 99), "at least 100",
               class = "dmar_count_invalid")
  expect_error(sem_exact_fit_bootstrap(fit, B = 99),
               "Monte Carlo p-value is reported only from at least 100",
               class = "dmar_count_invalid")
  expect_error(sem_exact_fit_bootstrap(fit, B = 2.5),
               class = "dmar_count_invalid")
  expect_error(reliability_H(fit, ci_method = "percentile", B = 399),
               "at least 400", class = "dmar_count_invalid")
  expect_error(average_variance_extracted(fit, ci_method = "percentile",
                                          B = 1000, conf_level = 0.99),
               "at least 2000", class = "dmar_count_invalid")
  # B is not read without a bootstrap interval.
  expect_no_error(average_variance_extracted(fit, B = 2.5))
})

test_that("ci_eta_squared_generalized() keeps its floor of 1000 and counts failed refits", {
  set.seed(113)
  d <- data.frame(g = gl(3, 8), y = rnorm(24))
  fit <- stats::aov(y ~ g, data = d)
  expect_error(suppressWarnings(
    ci_eta_squared_generalized(fit, method = "bootstrap", B = 999)),
    "at least 1000", class = "dmar_count_invalid")
  expect_error(suppressWarnings(
    ci_eta_squared_generalized(fit, method = "bootstrap", B = 999)),
    "never from fewer than 1000 replications", class = "dmar_count_invalid")
  expect_error(suppressWarnings(
    ci_eta_squared_generalized(fit, method = "bootstrap", B = 1999,
                               conf_level = 0.99)),
    "at least 2000", class = "dmar_count_invalid")
  expect_error(suppressWarnings(
    ci_eta_squared_generalized(fit, method = "bootstrap", B = 1000.5)),
    class = "dmar_count_invalid")
  res <- suppressWarnings(
    ci_eta_squared_generalized(fit, method = "bootstrap", B = 1000,
                               seed = 113))
  expect_identical(attr(res, "B"), 1000L)
  expect_identical(attr(res, "B_used"), 1000L)
  expect_identical(attr(res, "n_dropped"), 0L)
  expect_identical(attr(res, "B_required"), 1000L)
  expect_identical(attr(res, "bootstrap_status"), "complete")
  expect_true(is.finite(res$ci_lower) && is.finite(res$ci_upper))
})

test_that("mediation_mbco() checks B only for the Monte Carlo interval", {
  skip_on_cran()
  skip_if_not_installed("OpenMx")
  skip_if_not_installed("lavaan")
  set.seed(113)
  n <- 120
  x <- rnorm(n); m <- 0.5 * x + rnorm(n); y <- 0.4 * m + rnorm(n)
  d <- data.frame(x = x, m = m, y = y)
  model <- "m ~ a*x\n y ~ b*m + cp*x"
  expect_error(mediation_mbco(model, d, x = "x", y = "y",
                              ci_method = "monte_carlo", B = 399),
               "at least 400", class = "dmar_count_invalid")
  expect_error(mediation_mbco(model, d, x = "x", y = "y",
                              ci_method = "monte_carlo", B = 2.5),
               class = "dmar_count_invalid")
  # The Wald interval does not read B; through 1.0.x a small B was an
  # error whatever the interval.
  res <- mediation_mbco(model, d, x = "x", y = "y", ci_method = "wald",
                        B = 10)
  expect_null(attr(res, "B"))
  mc <- mediation_mbco(model, d, x = "x", y = "y",
                       ci_method = "monte_carlo", B = 400, seed = 113)
  expect_identical(attr(mc, "B"), 400L)
  expect_identical(attr(mc, "B_used"), 400L)
  expect_identical(attr(mc, "n_dropped"), 0L)
  expect_identical(attr(mc, "B_required"), 400L)
  expect_identical(attr(mc, "bootstrap_status"), "complete")
  expect_output(print(mc), "Monte Carlo interval: 400 replications, 400 used")
})

test_that("every exported function of the lane with a B formal refuses B = 2.5", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("lme4")
  set.seed(113)
  X <- matrix(rnorm(40 * 4), 40, 4) + rnorm(40)
  colnames(X) <- paste0("x", 1:4)
  d3 <- data.frame(x = X[, 1], m = X[, 2], y = X[, 3])
  fit <- lavaan::cfa(
    "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
     textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    data = holzinger_swineford, std.lv = TRUE)
  bif <- lavaan::cfa(
    "g =~ t1_visual_perception + t2_cubes + t4_lozenges +
          t6_paragraph_comprehension + t7_sentence + t9_word_meaning
     s1 =~ t1_visual_perception + t2_cubes + t4_lozenges
     s2 =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    data = holzinger_swineford, std.lv = TRUE, orthogonal = TRUE,
    check.post = FALSE, warn = FALSE)
  lmm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
  aovfit <- stats::aov(y ~ g, data = data.frame(g = gl(3, 8), y = rnorm(24)))
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  B <- 2.5
  calls <- list(
    mediate = function() mediate(d3, "x", "m", "y", B = B),
    reliability_alpha = function() reliability_alpha(X, ci_method = "percentile", B = B),
    reliability_kr20 = function() reliability_kr20(1 * (X > 0), ci_method = "percentile", B = B),
    reliability_omega = function() reliability_omega(X, ci_method = "percentile", B = B),
    reliability_omega_categorical = function() reliability_omega_categorical(
      as.data.frame(lapply(as.data.frame(X), function(v) as.integer(cut(v, 3)))),
      ci_method = "percentile", B = B),
    reliability_H = function() reliability_H(fit, ci_method = "percentile", B = B),
    ci_bifactor_indices = function() ci_bifactor_indices(bif, B = B),
    average_variance_extracted = function() average_variance_extracted(fit, ci_method = "percentile", B = B),
    fraction_missing_information = function() fraction_missing_information(fit, ci_method = "percentile", B = B),
    cohen_kappa = function() cohen_kappa(c("a", "b", "a", "b"), c("a", "b", "b", "b"), ci_method = "percentile", B = B),
    fleiss_kappa = function() fleiss_kappa(cbind(c(2, 0, 1, 2), c(0, 2, 1, 0)), ci_method = "percentile", B = B),
    krippendorff_alpha = function() krippendorff_alpha(cbind(c(1, 2, 3, 1), c(1, 2, 3, 2)), boot = TRUE, B = B),
    htmt = function() htmt(as.data.frame(X), list(A = c("x1", "x2"), B = c("x3", "x4")), B = B),
    ci_ratio_of_means = function() ci_ratio_of_means(abs(X[, 1]) + 1, abs(X[, 2]) + 1, ci_method = "percentile", B = B),
    ci_ratio_of_means_paired = function() ci_ratio_of_means_paired(abs(X[, 1]) + 1, abs(X[, 2]) + 1, ci_method = "percentile", B = B),
    R2_mixed_effects = function() R2_mixed_effects(lmm, ci_method = "boot", B = B),
    ci_eta_squared_generalized = function() suppressWarnings(ci_eta_squared_generalized(aovfit, method = "bootstrap", B = B)),
    mlmr = function() mlmr(y ~ x, data = d3, ci_method = "boot", B = B),
    mlmr_mv = function() mlmr_mv(cbind(y, m) ~ x, data = d3, ci_method = "boot", B = B),
    parallel_analysis = function() parallel_analysis(holzinger_swineford, items = items, B = B),
    factor_count = function() factor_count(holzinger_swineford, items = items, B = B),
    sem_exact_fit_bootstrap = function() sem_exact_fit_bootstrap(fit, B = B),
    ss_aipe_indirect_effect = function() ss_aipe_indirect_effect(
      true_std_reg_coef_a = 0.4, true_std_reg_coef_b = 0.4, width = 0.3, method = "monte_carlo", B = B, G = 10),
    ss_aipe_indirect_effect_sensitivity = function()
      ss_aipe_indirect_effect_sensitivity(
        true_std_reg_coef_a = 0.4, true_std_reg_coef_b = 0.3, specified_N = 50, width = 0.3,
        method = "monte_carlo", B = B, G = 2),
    ss_aipe_reliability_sensitivity = function()
      ss_aipe_reliability_sensitivity(
        true_reliability = 0.8, specified_N = 40, n_items = 4, width = 0.2,
        estimator = "alpha", ci_method = "percentile", B = B, G = 2))
  for (nm in names(calls)) {
    expect_error(calls[[nm]](), class = "dmar_count_invalid",
                 label = paste0(nm, "(B = 2.5)"))
  }
  # Every exported function with a B formal is either tested above or
  # named here with the reason it is not: reliability() passes B on to the
  # functions above; mediation_mbco() and plot_mediation_mbco() need
  # OpenMx and are tested where it is available.
  with_B <- Filter(function(nm) {
    f <- get(nm, envir = asNamespace("DMAR"))
    is.function(f) && "B" %in% names(formals(f))
  }, getNamespaceExports("DMAR"))
  elsewhere <- c("reliability", "mediation_mbco", "plot_mediation_mbco")
  expect_setequal(setdiff(with_B, elsewhere), names(calls))
})

test_that("average_variance_extracted() and fraction_missing_information() hold the count to the request and withhold when more than half are dropped", {
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa("f =~ t6_paragraph_comprehension+t7_sentence+t9_word_meaning",
                     data = holzinger_swineford)
  state <- new.env(parent = emptyenv())
  # lavBootstrap() replaced by fixed replicates, of which only the first
  # `returned` have a value.
  testthat::local_mocked_bindings(lavBootstrap = function(object, r, fun, ...) {
    v <- fun(object)
    x <- matrix(NA_real_, r, length(v))
    keep <- seq_len(state$returned)
    x[keep, ] <- rep(v, each = length(keep)) +
      seq(-0.05, 0.05, length.out = length(keep)) * (seq_along(v) <= 1L)
    x
  }, .package = "lavaan")
  point <- average_variance_extracted(fit)
  # 399 of 600: more than half returned, fewer than the 400 required. The
  # count is held to the 600 requested, so the interval is reported.
  state$returned <- 399L
  got <- with_warnings(average_variance_extracted(fit, ci_method = "percentile",
                                                  B = 600))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_identical(cls, "dmar_bootstrap_replications_dropped")
  w <- got$warnings[[1L]]
  expect_identical(w$status, "dropped")
  expect_identical(w$B_required, 400L)
  expect_match(conditionMessage(w),
               "the 399 that remained are fewer, so request more replications where some can fail",
               fixed = TRUE)
  x <- got$value
  expect_equal(x$ave, point$ave)
  expect_true(is.finite(x$ci_lower) && is.finite(x$ci_upper))
  expect_identical(attr(x, "B"), 600L)
  expect_identical(attr(x, "B_used"), 399L)
  expect_identical(attr(x, "n_dropped"), 201L)
  expect_identical(attr(x, "bootstrap_status"), "dropped")
  expect_identical(attr(x, "bootstrap")$minimum_returned, 300L)
  # 299 of 600: fewer than half returned, so no interval.
  state$returned <- 299L
  got <- with_warnings(average_variance_extracted(fit, ci_method = "percentile",
                                                  B = 600))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_identical(cls, "dmar_bootstrap_insufficient")
  expect_identical(got$warnings[[1L]]$status, "insufficient_share")
  expect_match(conditionMessage(got$warnings[[1L]]),
               "at least 300 are required (half of those requested",
               fixed = TRUE)
  expect_true(is.na(got$value$ci_lower) && is.na(got$value$ci_upper))
  expect_equal(got$value$ave, point$ave)
  # 400 of 600: formed, with the class of a run that dropped replications.
  state$returned <- 400L
  got <- with_warnings(average_variance_extracted(fit, ci_method = "percentile",
                                                  B = 600))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_identical(cls, "dmar_bootstrap_replications_dropped")
  expect_true(is.finite(got$value$ci_lower))
  expect_identical(attr(got$value, "bootstrap_status"), "dropped")

  # The same for the fraction of missing information, which stopped below
  # two usable replications through the first versions of 1.1.0.
  d <- holzinger_swineford[c("t6_paragraph_comprehension", "t7_sentence",
                             "t9_word_meaning")]
  d[c(3, 30, 77, 120, 200), 2L] <- NA
  fitm <- lavaan::cfa("f =~ t6_paragraph_comprehension+t7_sentence+t9_word_meaning",
                      data = d, missing = "ml")
  point <- fraction_missing_information(fitm)
  state$returned <- 399L
  got <- with_warnings(fraction_missing_information(fitm,
                                                    ci_method = "percentile",
                                                    B = 600))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_identical(cls, "dmar_bootstrap_replications_dropped")
  res <- got$value
  expect_true(all(is.finite(res$fmi_ci_lower)) &&
                all(is.finite(res$fmi_ci_upper)))
  expect_equal(res$fmi, point$fmi)
  expect_identical(attr(res, "B_used"), 399L)
  expect_identical(attr(res, "bootstrap_status"), "dropped")
  state$returned <- 299L
  got <- with_warnings(fraction_missing_information(fitm,
                                                    ci_method = "percentile",
                                                    B = 600))
  expect_identical(vapply(got$warnings, function(w) class(w)[1L], ""),
                   "dmar_bootstrap_insufficient")
  res <- got$value
  expect_true(all(is.na(res$fmi_ci_lower)) && all(is.na(res$fmi_ci_upper)))
  expect_equal(res$fmi, point$fmi)
  expect_identical(attr(res, "bootstrap_status"), "insufficient_share")
  expect_output(print(res), "more than half were dropped, so no interval is reported")
  state$returned <- 0L
  got <- with_warnings(fraction_missing_information(fitm,
                                                    ci_method = "percentile",
                                                    B = 600))
  expect_identical(vapply(got$warnings, function(w) class(w)[1L], ""),
                   "dmar_bootstrap_insufficient")
  expect_identical(attr(got$value, "bootstrap_status"), "insufficient_share")
})

test_that("sem_exact_fit_bootstrap() holds the count of 100 to the request and withholds a p-value when more than half are set aside", {
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa(
    "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
     textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    data = holzinger_swineford, std.lv = TRUE)
  state <- new.env(parent = emptyenv()); state$n <- 0L; state$fail_every <- 3L
  # Every sample returns a fixed statistic; every third "does not converge".
  testthat::local_mocked_bindings(.bs_refit_one = function(...) {
    state$n <- state$n + 1L
    if (state$n %% state$fail_every == 0L) {
      return(c(statistic = NA_real_, status = 3))
    }
    c(statistic = 5 + state$n / 100, status = 1)
  }, .package = "DMAR")
  # 120 samples, 40 set aside: 80 usable, more than half but fewer than
  # 100. The count is held to the 120 requested, so the p-value is
  # reported, and the warning says that fewer than 100 remained.
  got <- with_warnings(sem_exact_fit_bootstrap(fit, B = 120, seed = 1))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_identical(cls, "dmar_bootstrap_replications_dropped")
  w <- got$warnings[[1L]]
  expect_identical(w$status, "dropped")
  expect_identical(w$B_required, 100L)
  expect_match(conditionMessage(w), "40 of 120 bootstrap samples were set aside")
  expect_match(conditionMessage(w),
               "the 80 that remained are fewer, so request more bootstrap samples where some can fail, enough that at least 100 remain",
               fixed = TRUE)
  x <- got$value
  v <- function(term) x$value[x$term == term]
  expect_true(is.finite(v("p_bootstrap")))
  expect_equal(v("B_used"), 80)
  expect_identical(attr(x, "B_required"), 100L)
  expect_identical(attr(x, "bootstrap_status"), "dropped")
  expect_identical(attr(x, "sem_fit")$bootstrap_status, "dropped")
  expect_identical(dmar_metadata(x)$bootstrap$required, 100L)
  # Two of every three set aside: fewer than half usable, no p-value.
  state$n <- 0L; state$fail_every <- 1L
  testthat::local_mocked_bindings(.bs_refit_one = function(...) {
    state$n <- state$n + 1L
    if (state$n %% 3L != 0L) return(c(statistic = NA_real_, status = 3))
    c(statistic = 5 + state$n / 100, status = 1)
  }, .package = "DMAR")
  got <- with_warnings(sem_exact_fit_bootstrap(fit, B = 120, seed = 1))
  expect_identical(vapply(got$warnings, function(w) class(w)[1L], ""),
                   "dmar_bootstrap_insufficient")
  expect_identical(got$warnings[[1L]]$status, "insufficient_share")
  expect_true(is.na(got$value$value[got$value$term == "p_bootstrap"]))
  testthat::local_mocked_bindings(.bs_refit_one = function(...) {
    state$n <- state$n + 1L
    if (state$n %% 3L == 0L) return(c(statistic = NA_real_, status = 3))
    c(statistic = 5 + state$n / 100, status = 1)
  }, .package = "DMAR")
  # 150 samples, 50 set aside: 100 usable, and the p-value is reported.
  state$n <- 0L
  got <- with_warnings(sem_exact_fit_bootstrap(fit, B = 150, seed = 1))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_identical(cls, "dmar_bootstrap_replications_dropped")
  expect_match(conditionMessage(got$warnings[[1L]]),
               "the p-value uses the 100 that remain. The 100 required are counted on the number requested (150), not on the number that remained; request more than 100 where bootstrap samples can fail.",
               fixed = TRUE)
  expect_true(is.finite(got$value$value[got$value$term == "p_bootstrap"]))
  expect_identical(attr(got$value, "bootstrap_status"), "dropped")
})

test_that("ss_aipe_reliability_sensitivity() takes its default B from the rule and passes B to the reliability functions", {
  # A stand-in for reliability_alpha() records the count it is given and
  # returns an interval at once, so no bootstrap runs here.
  state <- new.env(parent = emptyenv())
  state$B <- integer(0)
  testthat::local_mocked_bindings(reliability_alpha = function(..., B) {
    state$B <- c(state$B, B)
    data.frame(term = c("estimate", "ci_lower", "ci_upper"),
               value = c(.8, .75, .85))
  }, .package = "DMAR")
  sens <- function(...) {
    ss_aipe_reliability_sensitivity(true_reliability = .8,
                                    supposed_reliability = .8, n_items = 4,
                                    width = .15, ...)
  }
  expect_null(formals(ss_aipe_reliability_sensitivity)$B)
  # The default is the larger of 200 and one and a quarter times the number
  # the method and the level require, which leaves room for replications
  # that return no value.
  defaults <- list(percentile = 500L, bca = 2500L, bootstrap_se = 200L,
                   bootstrap_se_logistic = 200L)
  for (method in names(defaults)) {
    state$B <- integer(0)
    set.seed(113)
    invisible(sens(G = 2, ci_method = method))
    expect_identical(state$B, rep(defaults[[method]], 2L), label = method)
  }
  state$B <- integer(0)
  set.seed(113)
  invisible(sens(G = 2, ci_method = "percentile", conf_level = .99))
  expect_identical(state$B, c(2500L, 2500L))
  state$B <- integer(0)
  set.seed(113)
  invisible(sens(G = 2, ci_method = "percentile", B = 600))
  expect_identical(state$B, c(600L, 600L))
  # A count below the requirement, or not a whole number, is refused
  # before any replication runs; a method that does not resample ignores B.
  for (bad in list(50, 399, 450.5, c(400, 500), NA_real_, "400")) {
    err <- expect_error(sens(G = 2, ci_method = "percentile", B = bad),
                        class = "dmar_count_invalid")
    expect_identical(err$min, 400L)
  }
  expect_error(sens(G = 2, ci_method = "bca", B = 1999),
               class = "dmar_count_invalid")
  expect_length(state$B, 2L)
  set.seed(113)
  expect_no_error(sens(G = 2, ci_method = "bonett", B = 50))
  # With the exploratory option a small B runs, with one warning for the
  # call and not one for each simulated study.
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  set.seed(113)
  got <- with_warnings(sens(G = 3, ci_method = "percentile", B = 150))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_identical(cls, "dmar_bootstrap_exploratory")
  expect_identical(got$warnings[[1L]]$B_required, 400L)
  # The request each simulated study ran is recorded on the table, and the
  # printed table carries the mark.
  expect_identical(attr(got$value, "B"), 150L)
  expect_identical(attr(got$value, "B_used"), 150L)
  expect_identical(attr(got$value, "n_dropped"), 0L)
  expect_identical(attr(got$value, "B_required"), 400L)
  expect_identical(attr(got$value, "bootstrap_status"), "exploratory")
  expect_match(paste(capture.output(print(got$value)), collapse = " "),
               "Each simulated study ran B = 150 replications, below the 400 required",
               fixed = TRUE)
})

test_that("the default B of ss_aipe_reliability_sensitivity() leaves simulated studies their intervals when replications fail", {
  skip_on_cran()
  # Omega of three items in 25 cases: some resamples have no one-factor
  # solution and return no value. With a default B equal to the 400
  # required, every simulated study lost its interval (12 of 12).
  res <- suppressWarnings(ss_aipe_reliability_sensitivity(
    true_reliability = .6, n_items = 3, width = .3, specified_N = 25,
    estimator = "omega", ci_method = "percentile", G = 4, seed = 4))
  expect_identical(attr(res, "B"), 500L)
  expect_identical(attr(res, "B_required"), 400L)
  expect_identical(attr(res, "bootstrap_status"), "complete")
  expect_gte(attr(res, "simulation")$returned, 1L)
  # A complete run prints no line about the replications.
  expect_no_match(paste(capture.output(print(res)), collapse = " "),
                  "Each simulated study ran", fixed = TRUE)
})

test_that("the two indirect effect planners mark an exploratory run once and record the request", {
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  plan <- with_warnings(ss_aipe_indirect_effect(
    .4, .4, width = .3, method = "monte_carlo", G = 20, B = 50, seed = 3))
  sens <- with_warnings(ss_aipe_indirect_effect_sensitivity(
    true_std_reg_coef_a = .4, true_std_reg_coef_b = .4, supposed_std_reg_coef_a = .4, supposed_std_reg_coef_b = .4,
    width = .3, method = "monte_carlo", G = 20, B = 50, seed = 3))
  for (got in list(plan, sens)) {
    cls <- vapply(got$warnings, function(w) class(w)[1L], "")
    expect_identical(cls, "dmar_bootstrap_exploratory")
    expect_identical(got$warnings[[1L]]$B_required, 400L)
    expect_identical(attr(got$value, "B"), 50L)
    expect_identical(attr(got$value, "B_used"), 50L)
    expect_identical(attr(got$value, "B_required"), 400L)
    expect_identical(attr(got$value, "bootstrap_status"), "exploratory")
    expect_match(paste(capture.output(print(got$value)), collapse = " "),
                 "Each simulated study ran B = 50 replications", fixed = TRUE)
  }
  expect_match(conditionMessage(plan$warnings[[1L]]),
               "^ss_aipe_indirect_effect\\(\\): ")
  expect_match(conditionMessage(sens$warnings[[1L]]),
               "^ss_aipe_indirect_effect_sensitivity\\(\\): ")
  # The closed form draws nothing and records nothing.
  closed <- ss_aipe_indirect_effect(.4, .4, width = .3)
  expect_null(attr(closed, "B"))
  expect_null(attr(closed, "bootstrap_status"))
})

test_that("a BCa interval is withheld by the count at its adjusted positions, in each family", {
  skip_on_cran()
  catch <- function(expr) {
    w <- NULL
    value <- withCallingHandlers(
      expr,
      dmar_bootstrap_insufficient = function(c) {
        w <<- c
        invokeRestart("muffleWarning")
      },
      warning = function(c) invokeRestart("muffleWarning"))
    list(value = value, w = w)
  }
  # mediate(): the lower position, 12.7 percent, has enough replications
  # beyond it; the upper one, 99.8 percent, does not.
  set.seed(31337)
  for (g in 1:16) {
    x <- rnorm(50)
    m <- .14 * x + rnorm(50)
    y <- .14 * m + rnorm(50)
  }
  med <- catch(mediate(data.frame(x, m, y), "x", "m", "y",
                       ci_method = "boot_bca", B = 2000, seed = 16))
  expect_s3_class(med$w, "dmar_bootstrap_insufficient")
  expect_identical(med$w$status, "insufficient_count")
  expect_identical(med$w$B_required, 5914L)
  expect_identical(med$w$B_suggest, 8900L)
  expect_equal(med$w$positions, c(0.127, 0.998), tolerance = 5e-3)
  expect_gt(med$w$positions[1L], 10 / 2000)
  expect_lt(1 - med$w$positions[2L], 10 / 2000)
  expect_identical(attr(med$value, "bootstrap_status"), "insufficient_count")
  expect_identical(attr(med$value, "B_required"), 5914L)
  expect_true(all(is.na(med$value$value[med$value$term %in%
                                          c("ci_lower", "ci_upper")])))

  # fleiss_kappa() at the 99 percent level with the 10,000 the request is
  # held to: withheld by the upper position.
  fl <- matrix(c(
    3L, 1L, 2L, 1L, 1L, 2L, 3L, 4L, 3L, 4L, 3L, 1L, 2L, 3L, 4L, 3L, 1L, 2L,
    2L, 2L, 3L, 1L, 2L, 1L, 2L, 0L, 0L, 2L, 2L, 3L, 2L, 1L, 2L, 4L, 2L, 0L,
    0L, 5L, 2L, 4L, 0L, 2L, 1L, 3L, 2L, 1L, 2L, 1L, 0L, 1L, 1L, 2L, 1L, 1L,
    1L, 2L, 1L, 1L, 2L, 0L, 0L, 2L, 2L, 3L, 1L, 3L, 2L, 3L, 2L, 0L, 2L, 3L,
    3L, 1L, 2L, 3L, 2L, 0L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 0L, 0L, 2L, 0L,
    1L, 2L, 2L, 1L, 0L, 0L, 3L, 2L, 1L, 3L, 2L, 2L, 1L, 1L, 2L, 2L, 3L, 0L,
    1L, 2L, 1L, 1L, 0L, 0L, 1L, 2L, 3L, 0L, 1L, 0L), nrow = 40L, ncol = 3L)
  fk <- catch(fleiss_kappa(fl, ci_method = "bca", B = 10000,
                           conf_level = .99, seed = 907))
  expect_s3_class(fk$w, "dmar_bootstrap_insufficient")
  expect_identical(fk$w$B_required, 11807L)
  expect_gt(fk$w$positions[1L], 10 / 10000)
  expect_lt(1 - fk$w$positions[2L], 10 / 10000)
  expect_true(is.na(fk$value$ci_lower) && is.na(fk$value$ci_upper))
  # The request for a BCa interval is held to 50 beyond each nominal limit.
  err <- expect_error(fleiss_kappa(fl, ci_method = "bca", B = 1999),
                      class = "dmar_count_invalid")
  expect_identical(err$min, 2000L)

  # The reliability family records the count at the positions, the number
  # its message states, in the condition, the attribute, the bootstrap
  # record, and the printed line.
  set.seed(607)
  f <- rnorm(8)
  X <- as.data.frame(sapply(1:3, function(j) .5 * f + rexp(8)^2))
  f <- rnorm(8)
  X <- as.data.frame(sapply(1:3, function(j) .5 * f + rexp(8)^2))
  rel <- catch(reliability_alpha(X, ci_method = "bca", B = 2000, seed = 2))
  expect_s3_class(rel$w, "dmar_bootstrap_insufficient")
  expect_match(conditionMessage(rel$w), "^reliability_alpha\\(\\): ")
  expect_match(conditionMessage(rel$w),
               "which takes a B of 5105 at those points (2000 were requested)",
               fixed = TRUE)
  expect_identical(rel$w$B_required, 5105L)
  expect_identical(rel$w$B_suggest, 7700L)
  expect_identical(attr(rel$value, "B_required"), 5105L)
  expect_identical(attr(rel$value, "bootstrap_status"), "insufficient_count")
  expect_identical(dmar_metadata(rel$value)$bootstrap$required, 5105L)
  expect_true(all(is.na(rel$value$value[rel$value$term %in%
                                          c("ci_lower", "ci_upper")])))
  expect_match(paste(capture.output(print(rel$value)), collapse = " "),
               paste0("5105 are required where the bias correction and the ",
                      "acceleration place the limits, so no BCa interval is ",
                      "reported"), fixed = TRUE)
})

test_that("plot_mediation_mbco() holds B to the rule and its warnings do not promise a record", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  skip_if_not_installed("ggplot2")
  set.seed(8)
  n <- 150
  d <- data.frame(x = rnorm(n), w = rnorm(n))
  d$m <- .4 * d$x + .2 * d$x * d$w + rnorm(n)
  d$y <- .4 * d$m + rnorm(n)
  res <- suppressWarnings(
    mediation_mbco("m ~ x + w + x:w \n y ~ m + x + w", data = d,
                   x = "x", y = "y", moderator = "w", ci_method = "wald"))
  err <- expect_error(plot_mediation_mbco(res, B = 399),
                      class = "dmar_count_invalid")
  expect_identical(err$min, 400L)
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  got <- with_warnings(plot_mediation_mbco(res, B = 50, seed = 1))
  expl <- Filter(function(w) inherits(w, "dmar_bootstrap_exploratory"),
                 got$warnings)
  expect_length(expl, 1L)
  expect_no_match(conditionMessage(expl[[1L]]), "recorded in the result",
                  fixed = TRUE)
})

test_that("messages name a bound, a reference quantile, and unequal BCa tails for what they are", {
  # A reference quantile of a parallel analysis.
  err <- expect_error(
    parallel_analysis(holzinger_swineford[7:15], B = 99, percentile = .5),
    class = "dmar_count_invalid")
  expect_match(conditionMessage(err),
               "A reference quantile at the 50th percentile of the reference data sets",
               fixed = TRUE)
  expect_error(parallel_analysis(holzinger_swineford[7:15], B = 200,
                                 percentile = 50),
               "'percentile' must be one number strictly between 0 and 1",
               fixed = TRUE)
  expect_error(factor_count(holzinger_swineford[7:15], percentile = 50),
               "'percentile' must be one number strictly between 0 and 1",
               fixed = TRUE)
  # A BCa interval with unequal tails: the nominal positions are 1 and 96
  # percent, and the message names the smaller tail area only.
  set.seed(1)
  x <- rnorm(60) + 5
  y <- rnorm(60) + 5
  err <- expect_error(
    ci_ratio_of_means(x, y, ci_method = "bca", alpha_lower = .01,
                      alpha_upper = .04, B = 4999, seed = 1),
    class = "dmar_count_invalid")
  expect_identical(err$min, 5000L)
  expect_match(conditionMessage(err),
               "the smaller nominal tail area is 1%", fixed = TRUE)
  expect_no_match(conditionMessage(err), "99%", fixed = TRUE)
  # No level passed: the sentence names the tail area and no level.
  err <- expect_error(.dmar_check_B(150, kind = "quantile", tail = .05),
                      class = "dmar_count_invalid")
  expect_no_match(conditionMessage(err), "A % interval", fixed = TRUE)
  expect_match(conditionMessage(err),
               "The limit that cuts off the smaller tail area has 5%",
               fixed = TRUE)
  # Under the exploratory option a single bound is called a bound, and
  # factor_count() gives the mark under its own name.
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  items <- c("t1_visual_perception", "t3_paper_form_board", "t4_lozenges",
             "t26_flags", "t10_addition", "t11_code",
             "t12_counting_groups_of_dots",
             "t13_straight_and_curved_capitals")
  got <- with_warnings(htmt(holzinger_swineford[items],
                            blocks = list(v = items[1:4], s = items[5:8]),
                            B = 30, seed = 1))
  expl <- Filter(function(w) inherits(w, "dmar_bootstrap_exploratory"),
                 got$warnings)
  expect_length(expl, 1L)
  expect_match(conditionMessage(expl[[1L]]),
               "With 30 replications the bound rests on about 1.5 of them",
               fixed = TRUE)
  got <- with_warnings(factor_count(holzinger_swineford[items], B = 30,
                                    seed = 1))
  expl <- Filter(function(w) inherits(w, "dmar_bootstrap_exploratory"),
                 got$warnings)
  expect_length(expl, 1L)
  expect_match(conditionMessage(expl[[1L]]), "^factor_count\\(\\): ")
  expect_match(conditionMessage(expl[[1L]]),
               "the reference quantile rests on about 1.5 of them",
               fixed = TRUE)
})

test_that("a withheld ratio interval prints that no interval is reported and has no form", {
  # A resample without the fourth score has a denominator of zero, about
  # 32 percent of them, so the share rule does not withhold; the run is
  # settled as one whose request was below the count.
  local_count_unmet()
  res <- suppressWarnings(ci_ratio_of_means(
    c(2, 3, 4, 5), c(0, 0, 0, 1), ci_method = "percentile", B = 400,
    seed = 12))
  expect_identical(attr(res, "bootstrap_status"), "insufficient_count")
  expect_true(is.na(attr(res, "confidence_set")))
  expect_true(is.na(generics::tidy(res)$confidence_set))
  shown <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(shown, "Percentile bootstrap interval: no interval is reported.",
               fixed = TRUE)
  expect_no_match(shown, "the whole line", fixed = TRUE)
  expect_no_match(shown, "Fieller set", fixed = TRUE)
  expect_no_match(shown, "..", fixed = TRUE)
  expect_match(shown, "were requested. A dropped replication had a denominator of zero",
               fixed = TRUE)
  paired <- suppressWarnings(ci_ratio_of_means_paired(
    c(2, 3, 4, 5), c(0, 0, 0, 1), ci_method = "percentile", B = 400,
    seed = 12))
  if (identical(attr(paired, "bootstrap_status"), "insufficient_count")) {
    expect_true(is.na(attr(paired, "confidence_set")))
    expect_no_match(paste(capture.output(print(paired)), collapse = "\n"),
                    "the whole line", fixed = TRUE)
  }
})

test_that("a note about a run with nothing dropped is a classed warning, and B is an integer", {
  set.seed(1)
  z <- rnorm(60)
  same <- data.frame(a = z, b = z, c = z)
  got <- with_warnings(reliability_alpha(same, ci_method = "percentile",
                                         B = 400, seed = 1))
  expect_true(any(vapply(got$warnings, inherits, logical(1),
                         what = "dmar_bootstrap_note")))
  expect_false(any(vapply(got$warnings, function(w)
    identical(class(w)[1L], "simpleWarning"), logical(1))))
  for (B in list(400, 400L)) {
    res <- reliability_alpha(holzinger_swineford[10:13],
                             ci_method = "percentile", B = B, seed = 1)
    expect_identical(attr(res, "B"), 400L)
  }
  # An exploratory run below the required number gives DMAR's classed
  # warning and does not also pass on boot.ci()'s unclassed warning that
  # extreme order statistics were used as endpoints.
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  got <- with_warnings(reliability_alpha(holzinger_swineford[10:13],
                                         ci_method = "percentile", B = 30,
                                         seed = 1))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_true("dmar_bootstrap_exploratory" %in% cls)
  expect_false(any(grepl("extreme order statistics",
                         vapply(got$warnings, conditionMessage, ""),
                         fixed = TRUE)))
})

test_that("the default footer prints the bootstrap counts once, for a table with no bootstrap line of its own", {
  set.seed(1)
  r1 <- sample(c("a", "b", "c"), 60, TRUE)
  r2 <- ifelse(stats::runif(60) < .7, r1, sample(c("a", "b", "c"), 60, TRUE))
  k <- cohen_kappa(r1, r2, ci_method = "percentile", B = 400, seed = 113)
  lines <- DMAR:::.dmar_tbl_footer(k)
  expect_identical(
    unname(lines["bootstrap_counts"]),
    "Bootstrap: 400 replications, 400 used (0 dropped; at least 400 must be requested)")
  expect_output(print(k), "Bootstrap: 400 replications, 400 used", fixed = TRUE)
  # A table without replications has no such line, and a family that
  # prints its own line prints one.
  expect_false("bootstrap_counts" %in% names(DMAR:::.dmar_tbl_footer(cohen_kappa(r1, r2))))
  x <- abs(stats::rnorm(80)) + 1
  y <- abs(stats::rnorm(80)) + 1
  ratio <- ci_ratio_of_means(x, y, ci_method = "percentile", B = 400, seed = 113)
  printed <- utils::capture.output(print(ratio))
  expect_identical(sum(grepl("Bootstrap: 400 replications", printed, fixed = TRUE)), 1L)
})

test_that("the Monte Carlo intervals of a mediation_mbco() call give at most one warning", {
  # The verdicts are those mediation_mbco() forms effect by effect; the
  # step that turns them into conditions is tested on its own, so that
  # the cases (drops in several effects, a withheld interval beside
  # reported ones) do not need a model in which they happen.
  settle <- function(B, used, term) {
    DMAR:::.dmar_boot_settle(
      B, used, tail = 0.025, conf_level = 0.95, fun = "mediation_mbco",
      quantity = paste0("the effect '", term, "'"),
      dropped_because = paste0("gave a value of '", term,
                               "' that is not a finite number"),
      replications = "Monte Carlo draws", signal = FALSE)
  }
  collect <- function(verdicts) {
    seen <- list()
    out <- withCallingHandlers(
      DMAR:::.mbco_mc_settle(verdicts),
      warning = function(w) {
        seen[[length(seen) + 1L]] <<- w
        invokeRestart("muffleWarning")
      })
    list(out = out, warnings = seen)
  }

  # Nothing dropped: no warning, status complete.
  clean <- collect(list(ab = settle(2000, 2000, "ab"),
                        total = settle(2000, 2000, "total")))
  expect_length(clean$warnings, 0L)
  expect_identical(clean$out$status, "complete")

  # One affected effect: its own message, as before.
  one <- collect(list(ab = settle(2000, 1990, "ab"),
                      total = settle(2000, 2000, "total")))
  expect_length(one$warnings, 1L)
  w <- one$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_replications_dropped")
  expect_match(conditionMessage(w), "^mediation_mbco\\(\\): 10 of 2000 ")
  expect_identical(w$effects$term, "ab")
  expect_identical(c(w$B, w$B_used, w$n_dropped), c(2000L, 1990L, 10L))
  expect_identical(one$out$status, "dropped")

  # Drops in two effects, every interval formed: one warning that names
  # both, with the counts of the effect that has the fewest usable draws.
  two <- collect(list(ab = settle(2000, 1990, "ab"),
                      ratio = settle(2000, 1900, "ratio"),
                      total = settle(2000, 2000, "total")))
  expect_length(two$warnings, 1L)
  w <- two$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_replications_dropped")
  msg <- conditionMessage(w)
  expect_match(msg, "2 of the 3 effects", fixed = TRUE)
  expect_match(msg, "'ab': 10 of 2000 dropped, 1990 used, interval reported",
               fixed = TRUE)
  expect_match(msg, "'ratio': 100 of 2000 dropped, 1900 used, interval reported",
               fixed = TRUE)
  expect_identical(lengths(regmatches(msg, gregexpr("mediation_mbco\\(\\)", msg))),
                   1L)
  expect_identical(w$effects$term, c("ab", "ratio"))
  expect_identical(w$effects$n_dropped, c(10L, 100L))
  expect_identical(w$effects$formed, c(TRUE, TRUE))
  expect_identical(c(w$B_used, w$n_dropped), c(1900L, 100L))
  expect_identical(w$status, "dropped")
  expect_identical(two$out$B_used, 1900L)

  # One interval withheld (more than half dropped) beside one with a few
  # drops: one warning, of the less favorable class.
  mixed <- collect(list(ab = settle(2000, 1990, "ab"),
                        ratio = settle(2000, 900, "ratio")))
  expect_length(mixed$warnings, 1L)
  w <- mixed$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_insufficient")
  expect_match(conditionMessage(w), "'ratio': 1100 of 2000 dropped, 900 used, no interval reported",
               fixed = TRUE)
  expect_identical(w$effects$formed, c(TRUE, FALSE))
  expect_identical(w$status, "insufficient_share")
  expect_identical(mixed$out$status, "insufficient_share")

  # An exploratory run below the required number of draws: one warning.
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  few <- collect(list(ab = settle(200, 200, "ab"),
                      total = settle(200, 200, "total")))
  expect_length(few$warnings, 1L)
  expect_s3_class(few$warnings[[1L]], "dmar_bootstrap_exploratory")
  few_dropped <- collect(list(ab = settle(200, 190, "ab"),
                              total = settle(200, 195, "total")))
  expect_length(few_dropped$warnings, 1L)
  expect_s3_class(few_dropped$warnings[[1L]],
                  "dmar_bootstrap_replications_dropped")
  expect_match(conditionMessage(few_dropped$warnings[[1L]]),
               "dmar.bootstrap_exploratory", fixed = TRUE)
})

# ---- paths the rule governs that had no test of their own ------------------
#
# Each function below drops and counts a replication that returns no value,
# reports the interval from the rest when no more than half were dropped,
# and withholds it past one half. The replications are replaced by fixed
# ones, so the counts are known.

# Runs `expr` with MASS::mvrnorm() returning its draws with the first
# `n_bad` rows missing (the draws of a Monte Carlo interval; the call with
# empirical = TRUE that builds data from moments is left alone). The draws
# are kept in `seen$draws`.
with_missing_draws <- function(n_bad, seen, expr) {
  real <- MASS::mvrnorm
  testthat::local_mocked_bindings(
    mvrnorm = function(n = 1, mu, Sigma, ..., empirical = FALSE) {
      out <- real(n, mu, Sigma, ..., empirical = empirical)
      if (!empirical) {
        if (n_bad > 0L) out[seq_len(n_bad), ] <- NA_real_
        seen$draws <- out
      }
      out
    }, .package = "MASS")
  expr
}

test_that("mediation_mbco() drops and counts Monte Carlo draws without a value and withholds the interval past one half", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  # The memory study moments of test-mediation_mbco.R (Tofighi & Kelley,
  # 2020; N = 369).
  vars <- c("instruction", "imagery", "repetition", "recall")
  M <- c(instruction = 0.5094850949, imagery = 5.6612466125,
         repetition = 6.0785907859, recall = 12.0650406504)
  S <- matrix(c(
     0.2505891363,  0.9121892306, -0.9531931189,  0.5428596677,
     0.9121892306,  8.7626517026, -4.7368917167,  5.1470926122,
    -0.9531931189, -4.7368917167,  8.0454371392, -2.7143646165,
     0.5428596677,  5.1470926122, -2.7143646165, 11.5283669141),
    nrow = 4L, dimnames = list(vars, vars))
  model <- "
    imagery ~ b1*instruction
    recall  ~ b2*imagery + b3*instruction
  "
  run <- function(n_bad, seen = new.env()) {
    with_missing_draws(n_bad, seen, with_warnings(
      mediation_mbco(model, S = S, M = M, N = 369L, x = "instruction",
                     y = "recall", ci_method = "monte_carlo", B = 400,
                     seed = 7)))
  }
  clean <- run(0L)
  expect_length(clean$warnings, 0L)
  expect_identical(attr(clean$value, "bootstrap_status"), "complete")

  # Seven draws without a value: dropped, counted, and reported, and the
  # limits are the quantiles of the 393 that remain, although 393 is
  # below the 400 a request must ask for.
  seen <- new.env()
  got <- run(7L, seen)
  res <- got$value
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_gte(length(cls), 1L)
  expect_true(all(cls == "dmar_bootstrap_replications_dropped"))
  w <- got$warnings[[1L]]
  expect_identical(w$B, 400L)
  expect_identical(w$B_required, 400L)
  expect_match(conditionMessage(w), "^mediation_mbco\\(\\): ")
  expect_match(conditionMessage(w), "of 400 Monte Carlo draws", fixed = TRUE)
  expect_identical(attr(res, "B"), 400L)
  expect_identical(attr(res, "B_used"), 393L)
  expect_identical(attr(res, "n_dropped"), 7L)
  expect_identical(attr(res, "B_required"), 400L)
  expect_identical(attr(res, "bootstrap_status"), "dropped")
  expect_true(all(c("b1", "b2") %in% colnames(seen$draws)))
  ab <- seen$draws[, "b1"] * seen$draws[, "b2"]
  expect_identical(sum(is.finite(ab)), 393L)
  ind <- res$term == "indirect_via_imagery"
  expect_equal(c(res$ci_lower[ind], res$ci_upper[ind]),
               quantile(ab[is.finite(ab)], c(0.025, 0.975), names = FALSE),
               tolerance = 1e-10)
  expect_equal(res$estimate, clean$value$estimate)
  expect_equal(res$lrt, clean$value$lrt)
  expect_output(print(res), "Monte Carlo interval: 400 replications, 393 used",
                fixed = TRUE)

  # More than half without a value: no interval for any effect, and the
  # estimates and the tests are untouched.
  got <- run(201L)
  res <- got$value
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_gte(length(cls), 1L)
  expect_true(all(cls == "dmar_bootstrap_insufficient"))
  expect_identical(got$warnings[[1L]]$status, "insufficient_share")
  expect_true(all(is.na(res$ci_lower)) && all(is.na(res$ci_upper)))
  expect_equal(res$estimate, clean$value$estimate)
  expect_equal(res$lrt, clean$value$lrt)
  expect_equal(res$p_value, clean$value$p_value)
  expect_identical(attr(res, "B_used"), 199L)
  expect_identical(attr(res, "bootstrap_status"), "insufficient_share")
  expect_output(print(res), "No interval: more than half", fixed = TRUE)
})

test_that("plot_mediation_mbco() draws the band from the draws that have a value and omits it past one half", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("OpenMx")
  skip_if_not_installed("ggplot2")
  set.seed(113)
  n <- 250
  x <- rnorm(n)
  w <- rnorm(n)
  m <- 0.5 * x + 0.3 * w + 0.4 * x * w + rnorm(n)
  y <- 0.5 * m + 0.2 * x + 0.1 * w + rnorm(n)
  d <- data.frame(x = x, w = w, m = m, y = y)
  res <- mediation_mbco("m ~ x + w + x:w \n y ~ m + x + w", data = d,
                        x = "x", y = "y", moderator = "w",
                        ci_method = "wald")
  draw <- function(n_bad) {
    with_missing_draws(n_bad, new.env(), with_warnings(
      plot_mediation_mbco(res, effects = "indirect_via_m", n_grid = 5,
                          B = 400, seed = 113)))
  }
  clean <- draw(0L)
  expect_length(clean$warnings, 0L)
  # Seven draws without a value: one warning, and the band is drawn.
  got <- draw(7L)
  expect_length(got$warnings, 1L)
  wn <- got$warnings[[1L]]
  expect_s3_class(wn, "dmar_bootstrap_replications_dropped")
  expect_identical(wn[c("B", "B_used", "n_dropped", "B_required")],
                   list(B = 400L, B_used = 393L, n_dropped = 7L,
                        B_required = 400L))
  expect_match(conditionMessage(wn), "^plot_mediation_mbco\\(\\): 7 of 400 Monte Carlo draws")
  # A plot records nothing, and its warning does not say that it does.
  expect_no_match(conditionMessage(wn), "are recorded in the result",
                  fixed = TRUE)
  band <- got$value$data
  expect_true(all(is.finite(band$band_lower)) &&
                all(is.finite(band$band_upper)))
  expect_true(all(band$band_lower <= band$estimate + 1e-8 &
                    band$band_upper >= band$estimate - 1e-8))
  expect_equal(band$estimate, clean$value$data$estimate)
  # More than half without a value: one warning, the curve without a band.
  got <- draw(201L)
  expect_length(got$warnings, 1L)
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_insufficient")
  expect_identical(got$warnings[[1L]]$status, "insufficient_share")
  band <- got$value$data
  expect_true(all(is.na(band$band_lower)) && all(is.na(band$band_upper)))
  expect_equal(band$estimate, clean$value$data$estimate)
})

test_that("ci_eta_squared_generalized() drops and counts a replication whose refit or effect size fails", {
  skip_on_cran()  # three runs of the replications
  fit <- stats::aov(len ~ supp * dose,
                    data = transform(ToothGrowth, dose = factor(dose)))
  # The effect size of chosen replications is made to fail. The function
  # computes the effect size of the data `state$before` times before the
  # first replication (counted in the clean run below), so call
  # state$before + b is replication b.
  state <- new.env(parent = emptyenv())
  state$before <- 0L
  real <- eta_squared_generalized
  failing <- function(bad) {
    state$n <- 0L
    testthat::local_mocked_bindings(
      eta_squared_generalized = function(...) {
        state$n <- state$n + 1L
        if ((state$n - state$before) %in% bad) stop("no effect size")
        real(...)
      }, .package = "DMAR", .env = parent.frame())
  }
  boot_class <- function(ws) {
    cls <- vapply(ws, function(w) class(w)[1L], "")
    cls[startsWith(cls, "dmar_bootstrap")]
  }
  call <- function() with_warnings(ci_eta_squared_generalized(
    fit, observed = "supp", method = "bootstrap", B = 1000, seed = 113))
  failing(integer())
  clean <- call()
  state$before <- state$n - 1000L
  expect_gte(state$before, 1L)
  expect_length(boot_class(clean$warnings), 0L)
  expect_identical(attr(clean$value, "bootstrap_status"), "complete")
  # 12 replications fail of the 1000 requested, this function's own
  # minimum: the limits are reported from the 988 that remain.
  failing(c(3, 40, 77, 120, 121, 300, 301, 302, 500, 777, 900, 1000))
  got <- call()
  res <- got$value
  expect_identical(boot_class(got$warnings),
                   "dmar_bootstrap_replications_dropped")
  w <- Filter(function(w) inherits(w, "dmar_bootstrap_replications_dropped"),
              got$warnings)[[1L]]
  expect_identical(w[c("B", "B_used", "n_dropped", "B_required")],
                   list(B = 1000L, B_used = 988L, n_dropped = 12L,
                        B_required = 1000L))
  expect_match(conditionMessage(w),
               "ci_eta_squared_generalized(): 12 of 1000 bootstrap replications (1.2%) gave a refit that failed or an effect size that could not be computed and were dropped.",
               fixed = TRUE)
  expect_match(conditionMessage(w),
               "the 988 that remained are fewer, so request more bootstrap replications where some can fail, enough that at least 1000 remain",
               fixed = TRUE)
  expect_false(w$caution)
  expect_identical(attr(res, "B_used"), 988L)
  expect_identical(attr(res, "n_dropped"), 12L)
  expect_identical(attr(res, "bootstrap_status"), "dropped")
  expect_true(all(is.finite(res$ci_lower)) &&
                all(is.finite(res$ci_upper)))
  expect_equal(res$eta_squared_generalized,
               clean$value$eta_squared_generalized)
  # The residual draws are those of the clean run, so the limits differ
  # from it only by the 12 replications left out.
  expect_false(isTRUE(all.equal(res$ci_lower, clean$value$ci_lower)))
  expect_equal(res$ci_lower, clean$value$ci_lower, tolerance = 0.05)
  # More than half fail: no limits, the estimates untouched.
  failing(seq_len(501))
  got <- call()
  expect_identical(boot_class(got$warnings), "dmar_bootstrap_insufficient")
  res <- got$value
  expect_true(all(is.na(res$ci_lower)) && all(is.na(res$ci_upper)))
  expect_equal(res$eta_squared_generalized,
               clean$value$eta_squared_generalized)
  expect_identical(attr(res, "B_used"), 499L)
  expect_identical(attr(res, "bootstrap_status"), "insufficient_share")
})

test_that("parallel_analysis() drops and counts a reference data set whose fit fails", {
  skip_on_cran()  # three runs of the replications
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  d <- holzinger_swineford[items]
  # Call 1 of the reduction is the observed matrix; call b + 1 is
  # reference data set b.
  state <- new.env(parent = emptyenv())
  real <- DMAR:::.pa_reduced_eigen
  failing <- function(bad) {
    state$n <- 0L
    testthat::local_mocked_bindings(
      .pa_reduced_eigen = function(R, N) {
        state$n <- state$n + 1L
        if ((state$n - 1L) %in% bad) stop("no reference fit")
        real(R, N)
      }, .package = "DMAR", .env = parent.frame())
  }
  clean <- with_warnings(parallel_analysis(d, B = 200, seed = 113))
  expect_length(clean$warnings, 0L)
  # Five reference data sets fail of the 200 requested: thresholds from
  # the 195 that remain, with the counts.
  failing(c(2, 50, 51, 120, 200))
  got <- with_warnings(parallel_analysis(d, B = 200, seed = 113))
  res <- got$value
  expect_length(got$warnings, 1L)
  w <- got$warnings[[1L]]
  expect_s3_class(w, "dmar_bootstrap_replications_dropped")
  expect_identical(w[c("B", "B_used", "n_dropped", "B_required", "status")],
                   list(B = 200L, B_used = 195L, n_dropped = 5L,
                        B_required = 200L, status = "dropped"))
  expect_match(conditionMessage(w),
               "195 of 200 parallel-analysis reference datasets returned defined eigenvalues.",
               fixed = TRUE)
  expect_match(conditionMessage(w),
               "the 195 that remained are fewer, so request more reference data sets where some can fail, enough that at least 200 remain",
               fixed = TRUE)
  counts <- attr(res, "counts")
  expect_identical(counts$attempted, 200L)
  expect_identical(counts$usable, 195L)
  expect_identical(counts$failed, 5L)
  expect_identical(unname(c(counts$failure_reasons)), 5L)
  expect_identical(names(counts$failure_reasons), "no reference fit")
  expect_identical(attr(res, "B_used"), 195L)
  expect_identical(attr(res, "bootstrap_status"), "dropped")
  expect_true(all(is.finite(res$reference_quantile)))
  expect_false(is.na(attr(res, "suggested_factors")))
  # The thresholds are the quantiles of the reference draws that returned
  # eigenvalues.
  draws <- attr(res, "reference_draws")
  usable <- stats::complete.cases(draws)
  expect_identical(which(!usable), c(2L, 50L, 51L, 120L, 200L))
  expect_equal(res$reference_quantile,
               unname(apply(draws[usable, ], 2L, stats::quantile,
                            probs = 0.95, names = FALSE)))
  expect_equal(res$estimate, clean$value$estimate)
  # More than half fail: no thresholds and no factor count; the observed
  # eigenvalues are returned.
  failing(seq_len(101))
  got <- with_warnings(parallel_analysis(d, B = 200, seed = 113))
  res <- got$value
  expect_length(got$warnings, 1L)
  expect_s3_class(got$warnings[[1L]], "dmar_bootstrap_insufficient")
  expect_identical(got$warnings[[1L]]$status, "insufficient_share")
  expect_match(conditionMessage(got$warnings[[1L]]),
               "Thresholds and the factor count are withheld: more than half failed.",
               fixed = TRUE)
  expect_true(all(is.na(res$reference_quantile)))
  expect_true(all(is.na(res$exceeds_reference)))
  expect_true(is.na(attr(res, "suggested_factors")))
  expect_equal(res$estimate, clean$value$estimate)
  expect_identical(attr(res, "counts")$usable, 99L)
})

test_that("mlmr_mv() drops a whole resample that returns no value and withholds the intervals past one half", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford
  f <- cbind(t6_paragraph_comprehension, t9_word_meaning) ~
    t5_general_information
  point <- mlmr_mv(f, data = d, ci_method = "wald", effect_sizes = FALSE)
  k <- length(coef(point))
  # lavBootstrap() replaced by fixed replicates with missing coefficients.
  set.seed(113)
  boot_mat <- matrix(rnorm(400 * k), 400, k)
  boot_mat[c(5, 50), 1L] <- NA
  boot_mat[399, k] <- NaN
  testthat::local_mocked_bindings(
    lavBootstrap = function(object, R = 1000L, r = R, ...) {
      boot_mat[seq_len(r), , drop = FALSE]
    }, .package = "lavaan")
  got <- with_warnings(mlmr_mv(f, data = d, ci_method = "boot", B = 400,
                               effect_sizes = FALSE))
  fit <- got$value
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_identical(sum(cls == "dmar_bootstrap_replications_dropped"), 1L)
  expect_false("dmar_bootstrap_insufficient" %in% cls)
  w <- got$warnings[[which(cls == "dmar_bootstrap_replications_dropped")]]
  expect_match(conditionMessage(w), "^mlmr_mv\\(\\): 3 of 400 bootstrap replications")
  expect_match(conditionMessage(w),
               "describes resamples in which every coefficient is defined",
               fixed = TRUE)
  expect_identical(w$n_dropped, 3L)
  expect_identical(fit$B, 400L)
  expect_identical(fit$B_used, 397L)
  expect_identical(fit$n_dropped, 3L)
  expect_identical(fit$B_required, 400L)
  expect_identical(fit$bootstrap_status, "dropped")
  # A resample is dropped for every coefficient, so all limits rest on
  # the same 397 resamples.
  keep <- stats::complete.cases(boot_mat)
  # (The limits are an array of coefficients by limit by outcome, and
  # the replicate columns run outcome by outcome.)
  expect_identical(dim(fit$ci)[2L], 2L)
  expect_equal(as.numeric(fit$ci[, 1L, ]),
               apply(boot_mat[keep, ], 2L, quantile, probs = 0.025,
                     names = FALSE))
  expect_equal(as.numeric(fit$ci[, 2L, ]),
               apply(boot_mat[keep, ], 2L, quantile, probs = 0.975,
                     names = FALSE))
  expect_equal(coef(fit), coef(point))
  expect_output(print(fit), "Bootstrap: 400 replications, 397 used",
                fixed = TRUE)
  expect_identical(summary(fit)$B_used, 397L)
  # More than half dropped: NA limits, coefficients untouched.
  boot_mat[1:201, 2L] <- NA
  got <- with_warnings(mlmr_mv(f, data = d, ci_method = "boot", B = 400,
                               effect_sizes = FALSE))
  cls <- vapply(got$warnings, function(w) class(w)[1L], "")
  expect_true("dmar_bootstrap_insufficient" %in% cls)
  expect_true(all(is.na(got$value$ci)))
  expect_equal(coef(got$value), coef(point))
  expect_identical(got$value$bootstrap_status, "insufficient_share")
  expect_output(print(got$value), "No interval: more than half", fixed = TRUE)
})

test_that("mlmr() does not repeat lavaan's count of failed bootstrap runs beside its own, and passes lavaan's other warnings", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford
  f <- t6_paragraph_comprehension ~ t5_general_information + t9_word_meaning
  set.seed(113)
  boot_mat <- matrix(rnorm(400 * 3), 400, 3)
  boot_mat[c(5, 50, 300), ] <- NA
  # lavBootstrap() gives the warning lavaan gives for runs that failed, in
  # its singular and its plural form, and one about something else.
  state <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    lavBootstrap = function(object, R = 1000L, r = R, ...) {
      warning(state$count_warning, call. = FALSE)
      warning("lavaan->lavBootstrap(): something else about this run.",
              call. = FALSE)
      boot_mat[seq_len(r), , drop = FALSE]
    }, .package = "lavaan")
  for (text in c(
    "lavaan->lav_bootstrap_internal():  \n   3 bootstrap runs failed or did not converge.",
    "lavaan->lav_bootstrap_internal():  \n   1 bootstrap run failed or did not converge.")) {
    state$count_warning <- text
    got <- with_warnings(mlmr(f, data = d, ci_method = "boot", B = 400,
                              effect_sizes = FALSE))
    msgs <- vapply(got$warnings, conditionMessage, "")
    cls <- vapply(got$warnings, function(w) class(w)[1L], "")
    # lavaan's count is not passed on.
    expect_false(any(grepl("failed or did not converge.", msgs,
                           fixed = TRUE) & !startsWith(cls, "dmar_")))
    # Its other warning is.
    expect_identical(sum(grepl("something else about this run", msgs,
                               fixed = TRUE)), 1L)
    # DMAR's classed warning states the count, once.
    expect_identical(sum(cls == "dmar_bootstrap_replications_dropped"), 1L)
    w <- got$warnings[[which(cls == "dmar_bootstrap_replications_dropped")]]
    expect_identical(w$n_dropped, 3L)
    expect_match(conditionMessage(w),
                 "mlmr(): 3 of 400 bootstrap replications (0.75%) gave a refit that failed, did not converge, or returned a coefficient that is not a finite number and were dropped.",
                 fixed = TRUE)
    expect_identical(got$value$B_used, 397L)
    expect_true(all(is.finite(got$value$ci)))
  }
})

test_that("categorical omega under multiple imputation follows the rule: the request, the replications set aside, and the share", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("mice")
  skip_if_not_installed("boot")
  # The fixture of test-reliability_omega_categorical_missing.R: five
  # ordinal items, N = 400, y4 and y5 missing at random.
  set.seed(113)
  N <- 400
  J <- 5
  eta <- rnorm(N)
  latent <- outer(eta, rep(0.7, J)) +
    matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - 0.7^2), J)
  full <- as.data.frame(apply(latent, 2, function(x)
    as.integer(cut(x, breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf),
                   labels = FALSE))))
  colnames(full) <- paste0("y", seq_len(J))
  set.seed(2026)
  z <- as.numeric(scale(rowMeans(full[, 1:3])))
  mar <- full
  mar$y4[runif(N) < plogis(-1.2 + 1.5 * z)] <- NA
  mar$y5[runif(N) < plogis(-1.2 + 1.5 * z)] <- NA
  # Two imputations of the observed data, and no rule that asks for more.
  testthat::local_mocked_bindings(
    .mi_imputations_floor = function(imputations) as.integer(imputations),
    .mi_imputations_rule = function(fmi_upper, cv = 0.05) NA_integer_,
    .package = "DMAR")
  call <- function(B, ...) with_warnings(reliability_omega_categorical(
    mar, ci_method = "percentile", missing = "mi", B = B, imputations = 2,
    seed = 113, ...))

  # The request is held to the count before anything is imputed.
  err <- expect_error(call(399), class = "dmar_count_invalid")
  expect_identical(err$min, 400L)
  expect_error(call(1000, conf_level = 0.99), "at least 2000",
               class = "dmar_count_invalid")
  expect_error(reliability_omega_categorical(
    mar, ci_method = "bca", missing = "mi", B = 1999, imputations = 2),
    "at least 2000", class = "dmar_count_invalid")

  # The bootstrap then impute loop replaced by fixed replications, of
  # which the first `n_bad` could not be imputed (status code 7).
  fixed_boot <- function(n_bad) {
    testthat::local_mocked_bindings(
      .boot_mi = function(data, items, estimate_fn, impute, m, maxit, B,
                          ncpus, progress = NULL) {
        v <- 0.80 + seq(-0.05, 0.05, length.out = B)
        sets <- cbind(v - 0.004, v + 0.004)
        t <- cbind(v, 1, v, sets, matrix(1, B, m))
        bad <- seq_len(n_bad)
        t[bad, c(1L, 3L, 3L + seq_len(m))] <- NA_real_
        t[bad, 2L] <- 7
        t[bad, 3L + m + seq_len(m)] <- 3
        list(t = unname(t),
             messages = c(rep("the chains diverged", n_bad),
                          rep(NA_character_, B - n_bad)),
             polr_fallbacks = 0L, ridge = 0L, ncpus = 1L)
      }, .package = "DMAR", .env = parent.frame())
  }
  boot_class <- function(ws) {
    cls <- vapply(ws, function(w) class(w)[1L], "")
    cls[startsWith(cls, "dmar_bootstrap")]
  }
  val <- function(r, term) r$value[r$term == term]

  # Nothing set aside at exactly the required B.
  fixed_boot(0L)
  got <- call(400)
  r <- got$value
  expect_false(any(c("dmar_bootstrap_replications_dropped",
                     "dmar_bootstrap_insufficient",
                     "dmar_bootstrap_exploratory") %in%
                     boot_class(got$warnings)))
  expect_identical(attr(r, "bootstrap")$bootstrap_status, "complete")
  expect_identical(attr(r, "bootstrap")$required, 400L)
  expect_identical(attr(r, "B_required"), 400L)
  expect_true(is.finite(val(r, "ci_lower")) &&
                is.finite(val(r, "ci_upper")))
  point <- val(r, "estimate")

  # Three replications set aside at exactly the required B: the interval
  # is reported from the 397 that returned a value (through the first
  # versions of 1.1.0 it was withheld), and the one warning names the
  # function, counts them, and says to request more.
  fixed_boot(3L)
  got <- call(400)
  r <- got$value
  expect_identical(boot_class(got$warnings),
                   "dmar_bootstrap_replications_dropped")
  w <- Filter(function(w) inherits(w, "dmar_bootstrap_replications_dropped"),
              got$warnings)[[1L]]
  expect_identical(w$fun, "reliability_omega_categorical")
  expect_identical(w[c("B", "B_used", "n_dropped", "B_required", "status")],
                   list(B = 400L, B_used = 397L, n_dropped = 3L,
                        B_required = 400L, status = "dropped"))
  expect_match(conditionMessage(w),
               "reliability_omega_categorical(): 3 of 400 bootstrap replications were set aside (3 could not be imputed)",
               fixed = TRUE)
  expect_match(conditionMessage(w),
               "the interval is computed from the 397 that returned a value; the 400 required are counted on the number requested (400), not on the number that remained; the 397 that remained are fewer, so request more replications where some can fail, enough that at least 400 remain",
               fixed = TRUE)
  bs <- attr(r, "bootstrap")
  expect_identical(bs$returned, 397L)
  expect_identical(bs$reasons[["imputation"]], 3L)
  expect_identical(bs$bootstrap_status, "dropped")
  expect_identical(bs$B_used, 397L)
  expect_true(bs$require_half)
  expect_identical(attr(r, "bootstrap_status"), "dropped")
  expect_identical(attr(r, "n_dropped"), 3L)
  expect_true(is.finite(val(r, "se")) && is.finite(val(r, "ci_lower")) &&
                is.finite(val(r, "ci_upper")))
  expect_lt(val(r, "ci_lower"), val(r, "ci_upper"))
  expect_equal(val(r, "estimate"), point)
  expect_output(print(r), "(3 set aside: 3 could not be imputed)",
                fixed = TRUE)

  # Fewer than half return a value: no interval and no standard error.
  fixed_boot(201L)
  got <- call(400)
  r <- got$value
  expect_identical(boot_class(got$warnings), "dmar_bootstrap_insufficient")
  w <- Filter(function(w) inherits(w, "dmar_bootstrap_insufficient"),
              got$warnings)[[1L]]
  expect_identical(w$status, "insufficient_share")
  expect_match(conditionMessage(w), "^reliability_omega_categorical\\(\\): ")
  expect_match(conditionMessage(w),
               "fewer than half returned a value, so no interval is reported",
               fixed = TRUE)
  expect_true(all(is.na(c(val(r, "se"), val(r, "ci_lower"),
                          val(r, "ci_upper")))))
  expect_equal(val(r, "estimate"), point)
  expect_identical(attr(r, "bootstrap")$bootstrap_status,
                   "insufficient_share")

  # Under the option a smaller request runs and is marked, once, with the
  # function's name; replications set aside are still reported.
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  fixed_boot(2L)
  got <- call(40)
  expect_identical(boot_class(got$warnings),
                   c("dmar_bootstrap_replications_dropped",
                     "dmar_bootstrap_exploratory"))
  w <- Filter(function(w) inherits(w, "dmar_bootstrap_exploratory"),
              got$warnings)[[1L]]
  expect_match(conditionMessage(w),
               "reliability_omega_categorical(): B = 40 is below the 400 replications a 95% interval is reported from",
               fixed = TRUE)
  expect_identical(attr(got$value, "bootstrap")$bootstrap_status,
                   "exploratory")
  expect_true(is.finite(val(got$value, "ci_lower")))
})

# ---- the reliability engine names the function that was called -------------

test_that("every warning of the reliability engine begins with the name of the function that was called and carries it as a field", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("boot")
  # The engine itself: without `fun` the text is as it was, and with it
  # the text begins with the name and the condition has the field.
  set.seed(113)
  Z <- as.data.frame(matrix(rnorm(40 * 3), 40, 3))
  point_fn <- function(d) {
    if (d[1L, 1L] > stats::quantile(Z[[1L]], 0.8)) stop("no value")
    mean(d[[1L]])
  }
  bare <- with_warnings(DMAR:::.bootstrap_ci(
    Z, point_fn, B = 400, conf_level = 0.95, kind = "percentile",
    seed = 113))$warnings[[1L]]
  named <- with_warnings(DMAR:::.bootstrap_ci(
    Z, point_fn, B = 400, conf_level = 0.95, kind = "percentile",
    seed = 113, fun = "some_function"))$warnings[[1L]]
  expect_match(conditionMessage(bare), "^[0-9]+ of 400 bootstrap replications were set aside")
  expect_null(bare$fun)
  expect_identical(conditionMessage(named),
                   paste0("some_function(): ", conditionMessage(bare)))
  expect_identical(named$fun, "some_function")
  expect_s3_class(named, "dmar_bootstrap_replications_dropped")
  # The ten calls of the engine in the package all pass the name of the
  # function they are in.
  engine_calls <- function(e) {
    if (is.call(e)) {
      here <- if (identical(e[[1L]], as.name(".bootstrap_ci"))) list(e)
      c(here, unlist(lapply(as.list(e), engine_calls), recursive = FALSE))
    } else if (is.pairlist(e)) {
      unlist(lapply(as.list(e), engine_calls), recursive = FALSE)
    } else {
      NULL
    }
  }
  n_calls <- 0L
  for (f in c("reliability_alpha", "reliability_kr20", "reliability_omega",
              "reliability_omega_categorical", "reliability_H",
              "ci_bifactor_indices")) {
    found <- engine_calls(body(get(f, envir = asNamespace("DMAR"))))
    expect_gte(length(found), 1L, label = paste0("engine calls in ", f))
    for (cl in found) {
      expect_identical(cl[["fun"]], f,
                       label = paste0("the fun argument of the engine call in ", f))
    }
    n_calls <- n_calls + length(found)
  }
  expect_identical(n_calls, 10L)

  # Through the exported functions: a run below the required number under
  # the option gives the engine's exploratory warning, and it names the
  # function.
  set.seed(113)
  X <- matrix(rnorm(60 * 4), 60, 4) + rnorm(60)
  colnames(X) <- paste0("x", 1:4)
  ord <- as.data.frame(lapply(as.data.frame(X), function(v)
    as.integer(cut(v, 3))))
  fit <- lavaan::cfa("f =~ x1 + x2 + x3 + x4", data = as.data.frame(X),
                     std.lv = TRUE)
  bif <- lavaan::cfa(
    "g =~ t1_visual_perception + t2_cubes + t4_lozenges +
          t6_paragraph_comprehension + t7_sentence + t9_word_meaning
     s1 =~ t1_visual_perception + t2_cubes + t4_lozenges
     s2 =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    data = holzinger_swineford, std.lv = TRUE, orthogonal = TRUE,
    check.post = FALSE, warn = FALSE)
  withr::local_options(dmar.bootstrap_exploratory = TRUE)
  calls <- list(
    reliability_alpha = function()
      reliability_alpha(X, ci_method = "percentile", B = 12, seed = 1),
    reliability_kr20 = function()
      reliability_kr20(1 * (X > 0), ci_method = "percentile", B = 12,
                       seed = 1),
    reliability_omega = function()
      reliability_omega(X, ci_method = "percentile", B = 12, seed = 1),
    reliability_omega_categorical = function()
      reliability_omega_categorical(ord, ci_method = "percentile", B = 6,
                                    seed = 1),
    reliability_H = function()
      reliability_H(fit, ci_method = "percentile", B = 12, seed = 1),
    ci_bifactor_indices = function()
      ci_bifactor_indices(bif, B = 6, seed = 1))
  for (nm in names(calls)) {
    got <- with_warnings(suppressMessages(calls[[nm]]()))
    ours <- Filter(function(w) startsWith(class(w)[1L], "dmar_bootstrap"),
                   got$warnings)
    expect_gte(length(ours), 1L, label = paste0(nm, "(): classed warnings"))
    for (w in ours) {
      expect_true(startsWith(conditionMessage(w), paste0(nm, "(): ")),
                  label = paste0(nm, "(): the warning begins with the name"))
      expect_identical(w$fun, nm, label = paste0(nm, "(): the fun field"))
    }
    expl <- Filter(function(w) inherits(w, "dmar_bootstrap_exploratory"),
                   ours)
    expect_length(expl, 1L)
    expect_match(conditionMessage(expl[[1L]]),
                 paste0("^", nm, "\\(\\): B = [0-9]+ is below the 400 ",
                        "replications a 95% interval is reported from"),
                 label = paste0(nm, "(): the exploratory warning"))
  }
  # reliability() passes the call on, and the warning names the function
  # that ran the bootstrap.
  got <- with_warnings(reliability(X, type = "alpha", ci_method = "percentile",
                                   B = 12, seed = 1))
  expl <- Filter(function(w) inherits(w, "dmar_bootstrap_exploratory"),
                 got$warnings)
  expect_length(expl, 1L)
  expect_match(conditionMessage(expl[[1L]]), "^reliability_alpha\\(\\): ")
  # A note about a run with nothing dropped is named too.
  withr::local_options(dmar.bootstrap_exploratory = FALSE)
  set.seed(1)
  z <- rnorm(60)
  got <- with_warnings(reliability_alpha(data.frame(a = z, b = z, c = z),
                                         ci_method = "percentile", B = 400,
                                         seed = 1))
  note <- Filter(function(w) inherits(w, "dmar_bootstrap_note"),
                 got$warnings)
  expect_length(note, 1L)
  expect_match(conditionMessage(note[[1L]]),
               "^reliability_alpha\\(\\): Every replication that returned a value gave the same estimate")
  expect_identical(note[[1L]]$fun, "reliability_alpha")
})
