# Live sandwich comparison for the cluster-robust covariance of
# ci_ratio_of_means_paired(). The pinned constants in
# tests/testthat/test-ci_ratio_of_means_paired.R were captured from these
# calls (sandwich 3.1.3, 2026-11-05). The cluster sum covariance with the
# factor Q / (Q - 1) is the variance of an intercept that vcovCL() gives
# with type = "HC1", whose (n - 1) / (n - k) factor is 1 for an
# intercept-only model; the covariance of the two means follows from
# var(a - b) = var(a) + var(b) - 2 cov(a, b). (An mlm fit of both at once
# would count two coefficients in that factor, so separate fits are used.)

## from tests/testthat/test-ci_ratio_of_means_paired.R
local({
  # The data of the test: units in twelve clusters with a shared effect,
  # generated and rounded to two places as they were for the pins.
  set.seed(19540512)
  cluster <- rep(1:12, times = c(4, 5, 3, 6, 4, 5, 4, 3, 6, 5, 4, 4))
  n <- length(cluster)
  shared <- stats::rnorm(12)[cluster]
  z <- round(10 + shared + stats::rnorm(n), 2)
  y <- round(z + 3 + 0.5 * shared + stats::rnorm(n), 2)
  x <- round(z + 2.4 + 0.5 * shared + stats::rnorm(n), 2)
  ours <- DMAR::ci_ratio_of_means_paired(x, y, baseline = z,
                                         cluster = cluster)
  v <- function(t) ours$value[ours$term == t]
  a <- x - z
  b <- y - z
  vcl <- function(w) {
    sandwich::vcovCL(stats::lm(w ~ 1), cluster = cluster, type = "HC1")[1, 1]
  }
  va <- vcl(a)
  vb <- vcl(b)
  cov_ab <- (va + vb - vcl(a - b)) / 2
  stopifnot(
    isTRUE(all.equal(v("se_numerator")^2, va, tolerance = 1e-10)),
    isTRUE(all.equal(v("se_denominator")^2, vb, tolerance = 1e-10)),
    isTRUE(all.equal(v("r") * v("se_numerator") * v("se_denominator"),
                     cov_ab, tolerance = 1e-10)),
    v("df") == 11,
    # The constants the test pins.
    isTRUE(all.equal(c(va, vb, cov_ab), c(0.0675418874878823,
                                          0.0692726405168686,
                                          0.0469585608693910),
                     tolerance = 1e-12))
  )
})
