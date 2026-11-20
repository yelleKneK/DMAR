# Live mratios comparisons for ci_ratio_of_means(). The pinned constants in
# tests/testthat/test-ci_ratio_of_means.R were captured from these calls
# (mratios 1.4.4, 2026-11-05). mratios is installed for this check only and
# is not a dependency of DMAR.

## from tests/testthat/test-ci_ratio_of_means.R
local({
  x <- c(37.4, 56.9, 53.1, 44.7, 56.1, 52.1, 52.5, 62.1, 58.9, 40.2, 64.6,
         48.2, 58, 54.7, 49.3, 52.7, 72.4, 56.9, 51.9, 65.9, 55.2, 36.1,
         48.2, 43.6)
  y <- c(31.4, 43.2, 47.7, 44.7, 56, 65.8, 66, 46.3, 41, 34.5, 43.2, 32.6,
         44.6, 46.1, 40.2, 54.5, 39, 27.5, 20, 50.5, -0.5, 40.1, 37, 41.9,
         61.3, 21, 35.5, 40.6, 52.8, 50.8)
  limits <- function(out) {
    out$value[out$term %in% c("ci_lower", "ci_upper")]
  }
  # Pooled variances: Fieller's interval, which ttestratio() computes in
  # closed form.
  for (level in c(0.95, 0.90)) {
    ours <- DMAR::ci_ratio_of_means(x, y, var_equal = TRUE,
                                    conf_level = level)
    theirs <- mratios::ttestratio(x, y, var.equal = TRUE, conf.level = level)
    stopifnot(isTRUE(all.equal(limits(ours), as.numeric(theirs$conf.int),
                               tolerance = 1e-10)),
              isTRUE(all.equal(ours$value[ours$term == "ratio"],
                               as.numeric(theirs$estimate[3]),
                               tolerance = 1e-12)))
  }
  # Unequal variances: ttestratio(iterativeCI = TRUE) inverts the same
  # Welch-type test, searching with uniroot() at its default tolerance,
  # which leaves its limits about 1e-5 from the exact inversion. DMAR's
  # limits are the exact ones: t.test(x, L * y) has p-value .05 at each.
  ours <- DMAR::ci_ratio_of_means(x, y)
  theirs <- mratios::ttestratio(x, y, var.equal = FALSE, iterativeCI = TRUE)
  stopifnot(
    max(abs(limits(ours) - as.numeric(theirs$conf.int))) < 1e-4,
    isTRUE(all.equal(vapply(limits(ours),
                            function(L) stats::t.test(x, L * y)$p.value, 0),
                     c(0.05, 0.05), tolerance = 1e-9))
  )
  one <- DMAR::ci_ratio_of_means(x, y, alpha_lower = 0.05, alpha_upper = 0)
  theirs <- mratios::ttestratio(x, y, var.equal = FALSE, iterativeCI = TRUE,
                                alternative = "greater")
  stopifnot(abs(limits(one)[1] - as.numeric(theirs$conf.int[1])) < 1e-4,
            is.infinite(limits(one)[2]))
})
