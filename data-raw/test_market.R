## data-raw/test_market.R
##
## Construction script for the `test_market` data set: the controlled
## test-market experiment of Bryant and Bruvold (1980, Table 1), used to
## illustrate multiple-comparison procedures in the analysis of covariance
## with a random covariate.
##
## Design: k = 6 marketing panels (treatments) crossed with s = 4 blocks
## of retail outlets (a randomized complete block). Within each panel-by-
## block cell there is one observation of two quantities:
##   brand_movement    (y) -- test-brand movement during the test period,
##                            in hundreds of statistical cases (the outcome).
##   category_movement (x) -- remaining category movement, the random
##                            concomitant (covariate) that becomes available
##                            only during the experiment.
##
## The values below are transcribed directly from Table 1 of Bryant and
## Bruvold (1980). This script verifies that the data reproduce every
## published summary quantity (panel/block marginals and the full ANCOVA:
## slope 0.4079, error MS 0.01326 on 14 df, and the adjusted panel means
## 3.595, 3.619, 4.102, 4.515, 4.618, 4.876).
##
## To rebuild from the package root:
##   source("data-raw/test_market.R")

y <- matrix(c(
  2.98, 3.29, 4.33, 3.48,
  2.99, 3.36, 3.99, 4.18,
  4.49, 4.25, 4.19, 3.89,
  4.30, 4.40, 3.83, 5.03,
  4.33, 5.24, 4.25, 5.24,
  3.91, 4.83, 4.91, 5.61), nrow = 6, byrow = TRUE)
x <- matrix(c(
   9.19, 13.78, 11.20, 13.88,
   9.61, 13.59,  9.99, 15.70,
  11.74, 14.21,  9.54, 14.31,
  10.43, 13.58,  7.88, 15.67,
  10.18, 15.73,  8.28, 16.04,
   8.72, 14.05,  8.91, 16.51), nrow = 6, byrow = TRUE)

test_market <- data.frame(
  panel             = factor(rep(1:6, each = 4)),
  block             = factor(rep(1:4, times = 6)),
  brand_movement    = as.vector(t(y)),
  category_movement = as.vector(t(x)),
  stringsAsFactors  = FALSE
)

# ---- Structural checks ----------------------------------------------
stopifnot(
  identical(dim(test_market), c(24L, 4L)),
  nlevels(test_market$panel) == 6L,
  nlevels(test_market$block) == 4L
)

# ---- Published marginal means (Table 1) -----------------------------
panel_y <- round(tapply(test_market$brand_movement, test_market$panel, mean), 3)
panel_x <- round(tapply(test_market$category_movement, test_market$panel, mean), 3)
stopifnot(
  all(abs(panel_y - c(3.520, 3.630, 4.205, 4.390, 4.765, 4.815)) < 5e-3),
  all(abs(panel_x - c(12.013, 12.223, 12.450, 11.890, 12.558, 12.048)) < 5e-3),
  all(abs(tapply(test_market$brand_movement, test_market$block, sum) -
            c(23.00, 25.37, 25.50, 27.43)) < 5e-3)
)

# ---- Reproduce the published ANCOVA ---------------------------------
fit <- lm(brand_movement ~ panel + block + category_movement, data = test_market)
slope     <- unname(coef(fit)["category_movement"])
error_ms  <- sum(resid(fit)^2) / fit$df.residual
adj_means <- vapply(levels(test_market$panel), function(p) {
  nd <- data.frame(panel = factor(p, levels = levels(test_market$panel)),
                   block = factor(1:4, levels = levels(test_market$block)),
                   category_movement = mean(test_market$category_movement))
  mean(predict(fit, nd))
}, numeric(1))

stopifnot(
  fit$df.residual == 14L,
  abs(slope    - 0.4079)  < 5e-4,
  abs(error_ms - 0.01326) < 5e-5,
  all(abs(adj_means - c(3.595, 3.619, 4.102, 4.515, 4.618, 4.876)) < 1e-3)
)

save(test_market,
     file = "data/test_market.rda",
     version = 2L, compress = "bzip2")

message("test_market: built and verified (24 x 4); slope ", round(slope, 4),
        ", error MS ", round(error_ms, 5), " on ", fit$df.residual, " df.")
