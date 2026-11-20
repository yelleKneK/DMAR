## from tests/testthat/test-bayes_t.R
local({
  v <- function(tab, t) tab$value[tab$term == t]
  set.seed(113)
  x <- rnorm(35, 0.3, 1); y <- rnorm(40, 0, 1)

  bf_one <- BayesFactor::extractBF(BayesFactor::ttestBF(x))$bf
  stopifnot(isTRUE(all.equal(v(DMAR::bayes_one_sample_t(x), "bf_10"),
                             bf_one, tolerance = 1e-4)))

  bf_two <- BayesFactor::extractBF(BayesFactor::ttestBF(x = x, y = y))$bf
  stopifnot(isTRUE(all.equal(v(DMAR::bayes_independent_t(x, y), "bf_10"),
                             bf_two, tolerance = 1e-4)))

  y2 <- x + rnorm(35, 0.2, 0.8)
  bf_pair <- BayesFactor::extractBF(
    BayesFactor::ttestBF(x = y2, y = x, paired = TRUE))$bf
  stopifnot(isTRUE(all.equal(v(DMAR::bayes_paired_t(y2, x), "bf_10"),
                             bf_pair, tolerance = 1e-4)))
})

## from tests/testthat/test-bayes_t.R
local({
  v <- function(tab, t) tab$value[tab$term == t]
  set.seed(113)
  x <- rnorm(50, 0.5, 1)
  res <- DMAR::bayes_one_sample_t(x)
  ch <- suppressMessages(
    BayesFactor::posterior(BayesFactor::ttestBF(x), iterations = 40000,
                           progress = FALSE))
  delta <- as.numeric(ch[, "delta"])
  stopifnot(isTRUE(all.equal(v(res, "delta_posterior_median"), median(delta),
                             tolerance = 0.02)))
  stopifnot(isTRUE(all.equal(v(res, "delta_posterior_mean"), mean(delta),
                             tolerance = 0.02)))
  stopifnot(isTRUE(all.equal(v(res, "delta_lower"),
                             unname(quantile(delta, 0.025)),
                             tolerance = 0.04)))
  stopifnot(isTRUE(all.equal(v(res, "delta_upper"),
                             unname(quantile(delta, 0.975)),
                             tolerance = 0.04)))
  stopifnot(isTRUE(all.equal(v(res, "p_delta_positive"), mean(delta > 0),
                             tolerance = 0.01)))
})
