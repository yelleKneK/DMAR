# The row of Table 3 of Kelley (2007, p. 762) that shows a misprint: 99
# percent confidence, assurance .80, width .010, coefficients of variation
# of .30 to .50. The row prints 14,279; 28,321; 50,209; 38,120; and 82,817,
# which do not increase. ss_aipe_cv() gives 14,277; 20,473; 28,320; 38,118;
# and 50,207 with the exact noncentral t distribution, and the article's
# procedure gives 14,279; 20,475; 28,321; 38,120; and 50,209 with R's
# approximation to that distribution, with which the tables were computed,
# so the printed 28,321 and 50,209 are the values for .40 and .50, one and
# two columns to the left of where they belong (?ss_aipe_cv, "The
# published tables").
#
# The five plans take about six minutes of processor time (about two and
# a half minutes when run on five cores, as here), which is why they are
# checked here and not in the test suite; the rest of the tables is in
# tests/testthat/test-ss_aipe_cv.R.
#
# Usage, from the package root, after any change to ss_aipe_cv(), ci_cv(),
# or R/nc_t_internals.R, and before a release:
#   Rscript tools/kelley2007_table3_row.R
# No step of tools/release_gate.R runs this script; if it stops,
# ?ss_aipe_cv ("The published tables") is out of date.
# The source tree in the working directory is loaded, so the row is
# computed by the code as it stands and not by an installed build. The
# script stops with an error if a plan differs from the value recorded
# here.
stopifnot(file.exists("DESCRIPTION"))
suppressMessages(devtools::load_all(".", quiet = TRUE))

kappa <- c(.30, .35, .40, .45, .50)
printed <- c(14279, 28321, 50209, 38120, 82817)
recorded <- c(14277, 20473, 28320, 38118, 50207)
approximate <- c(14279, 20475, 28321, 38120, 50209)

## 1. The function, on the five cells at once.
t0 <- proc.time()[["elapsed"]]
plan <- unlist(parallel::mclapply(kappa, function(k)
  ss_aipe_cv(true_cv = k, width = .010, conf_level = .99, assurance = .80)$value,
  mc.cores = min(5L, parallel::detectCores())))
cat(sprintf("ss_aipe_cv(): %s (%.0f seconds)\n",
            paste(format(plan, big.mark = ","), collapse = "; "),
            proc.time()[["elapsed"]] - t0))

## 2. The article's procedure in base R with pt() and qt(), as in 2007: the
##    expected-width search, the quantile of the assurance step, and the
##    search again for the value of Equation 3 at that quantile.
eq3 <- function(k, N) k * (1 - 1 / (4 * (N - 1)) + k^2 / N + 1 / (2 * (N - 1)^2))
width_at <- function(k, N, level) {
  a <- (1 - level) / 2; nu <- N - 1; th <- sqrt(N) / k
  s <- 10 * sqrt(1 + th^2 / (2 * nu)) + 10
  upper_ncp <- uniroot(function(l) pt(th, nu, l) - a, c(th - s, th + s),
                       tol = 1e-12)$root
  lower_ncp <- if (pt(th, nu, 0) <= 1 - a) 0 else
    uniroot(function(l) pt(th, nu, l) - (1 - a), c(0, th + s), tol = 1e-12)$root
  (if (lower_ncp <= 0) Inf else sqrt(N) / lower_ncp) - sqrt(N) / upper_ncp
}
# The first N from 4 whose width meets w, by doubling, bisection, and a
# look at the sizes just below.
search <- function(k, w, level) {
  wid <- function(N) width_at(eq3(k, N), N, level)
  if (wid(4) <= w) return(4)
  lo <- 4; hi <- 8
  while (wid(hi) > w) { lo <- hi; hi <- 2 * hi }
  while (hi - lo > 1) {
    mid <- (lo + hi) %/% 2
    if (wid(mid) > w) lo <- mid else hi <- mid
  }
  for (n in max(4, hi - 15):(hi - 1)) if (wid(n) <= w) return(n)
  hi
}
# Past a noncentrality of 37.62 pt() and qt() warn that full precision may
# not have been achieved. That approximation is what is being reproduced,
# so the warnings are muffled here and nowhere else.
by_approximation <- suppressWarnings(vapply(kappa, function(k) {
  N0 <- search(k, .010, .99)
  k_gamma <- sqrt(N0) / qt(1 - .80, N0 - 1, sqrt(N0) / k)
  search(eq3(k_gamma, N0), .010, .99)
}, numeric(1)))
cat("With pt() and qt():", paste(format(by_approximation, big.mark = ","),
                                  collapse = "; "), "\n")
cat("Printed:           ", paste(format(printed, big.mark = ","),
                                  collapse = "; "), "\n")

stopifnot(identical(as.numeric(plan), recorded),
          identical(by_approximation, approximate),
          identical(approximate[c(3, 5)], printed[c(2, 3)]),
          identical(approximate[c(1, 4)], printed[c(1, 4)]))
cat("OK: the function returns the recorded plans, and the printed 28,321 and",
    "50,209 are the approximation's values for .40 and .50.\n")
