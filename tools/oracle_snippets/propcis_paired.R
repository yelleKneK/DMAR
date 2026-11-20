# Live PropCIs comparisons for the paired proportion family, written
# against PropCIs 0.3.0 (2026-11-05). PropCIs::scoreci.mp() is Tango's
# (1998) score interval for the difference of two paired proportions; its
# arguments are the discordant counts b = n_01 and c = n_10. It bisects for
# each limit and stops when the step falls below 1e-7 or the statistic is
# within 1e-6 of its target, so in small tables, where the statistic is
# flat, its limits can sit a few times 1e-7 from the root (3.6e-7 at worst
# over the sweep below, where DMAR's limits solve Z = +/- q to 2e-11 and
# PropCIs's to 1.4e-5). That sets the tolerance at 1e-6.

## from tests/testthat/test-ci_proportion_paired.R
local({
  # The pinned cases of the test, then a sweep over small and moderate
  # tables, both levels the equivalence family reads, and every boundary
  # (no discordant pair, all pairs discordant one way, n_10 = n).
  cases <- rbind(
    c(38, 52, 400, 0.90), c(38, 52, 400, 0.95), c(12, 5, 100, 0.90),
    c(0, 0, 50, 0.95), c(3, 0, 20, 0.90), c(0, 7, 30, 0.90),
    c(1, 9, 15, 0.90), c(40, 0, 40, 0.95), c(0, 40, 40, 0.95),
    c(32, 27, 344, 0.95), c(20, 40, 160, 0.95))
  for (n in c(1, 2, 5, 12, 30)) {
    for (n_10 in 0:n) for (n_01 in 0:(n - n_10)) {
      cases <- rbind(cases, c(n_10, n_01, n, 0.90), c(n_10, n_01, n, 0.95))
    }
  }
  worst <- 0
  for (i in seq_len(nrow(cases))) {
    r <- cases[i, ]
    ours <- DMAR::ci_proportion_paired(n_10 = r[1], n_01 = r[2], N = r[3],
                                       conf_level = r[4])
    lims <- ours$value[match(c("ci_lower", "ci_upper"), ours$term)]
    theirs <- PropCIs::scoreci.mp(b = r[2], c = r[1], n = r[3],
                                  conf.level = r[4])$conf.int
    worst <- max(worst, abs(lims - as.numeric(theirs)))
  }
  stopifnot(worst < 1e-6)
})

## from tests/testthat/test-ci_proportion_paired.R
local({
  # The exact coverage the help page quotes, .986 for the 95% score
  # interval at p_10 = .03, p_01 = .01, n = 50, summed over every possible
  # table with PropCIs's limits in place of DMAR's (.98612 on 2026-11-05).
  # The test computes it from an independent inversion of the score test.
  p10 <- 0.03; p01 <- 0.01; n <- 50; psi <- p10 - p01
  g <- expand.grid(n_10 = 0:n, n_01 = 0:n)
  g <- g[g$n_10 + g$n_01 <= n, ]
  pr <- mapply(function(a, b) {
    stats::dmultinom(c(a, b, n - a - b), prob = c(p10, p01, 1 - p10 - p01))
  }, g$n_10, g$n_01)
  covers <- mapply(function(a, b) {
    l <- PropCIs::scoreci.mp(b = b, c = a, n = n, conf.level = 0.95)$conf.int
    l[1] <= psi && psi <= l[2]
  }, g$n_10, g$n_01)
  stopifnot(round(sum(pr * covers), 3) == 0.986)
})

## from tests/testthat/test-equivalence_proportion_paired.R
local({
  # The equivalence decisions are read from the 1 - 2 alpha score
  # interval, so its limits must match PropCIs at that level as well.
  for (a in c(0.025, 0.05, 0.10)) {
    res <- DMAR::equivalence_proportion_paired(n_10 = 32, n_01 = 27,
                                               N = 344, prop_diff_lower = 0.05,
                                               alpha_level = a)
    lims <- res$value[match(c("ci_lower", "ci_upper"), res$term)]
    theirs <- PropCIs::scoreci.mp(b = 27, c = 32, n = 344,
                                  conf.level = 1 - 2 * a)$conf.int
    stopifnot(max(abs(lims - as.numeric(theirs))) < 1e-6)
  }
})
