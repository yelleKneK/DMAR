## ci_two_phase(): the regression estimator of double sampling and the
## prediction-powered estimate, each checked against its closed form computed
## here independently (lm() and predict() for the slope, the residual
## variance, and the estimate, base R for the rest), plus the fast path,
## missing values, clusters (an exact identity with one unit per cluster and
## Monte Carlo checks of the linearization), and refusals.

.two_phase_data <- function(n = 150, N = 1000, seed = 11) {
  set.seed(seed)
  x_all <- rnorm(n + N, mean = 50, sd = 10)
  y_all <- 20 + 0.6 * x_all + rnorm(n + N, sd = 5)
  list(y = y_all[1:n], x = x_all[1:n], xu = x_all[-(1:n)])
}

.v <- function(res) stats::setNames(res$value, res$term)

test_that("the regression interval matches the closed-form double sampling variance", {
  d <- .two_phase_data()
  n <- length(d$y); N <- length(d$xu); M <- n + N
  fit <- lm(d$y ~ d$x)
  b <- unname(coef(fit)[2])
  s2_e <- summary(fit)$sigma^2                 # residual mean square, n - 2 df
  x_all <- c(d$x, d$xu)
  theta <- mean(d$y) + b * (mean(x_all) - mean(d$x))
  se <- sqrt(s2_e / n + b^2 * var(x_all) / M)
  q <- qt(0.975, n - 2)

  v <- .v(ci_two_phase(d$y, d$x, d$xu))
  expect_equal(v[["mean_Y"]], theta, tolerance = 1e-12)
  expect_equal(v[["se"]], se, tolerance = 1e-12)
  expect_equal(v[["df"]], n - 2)
  expect_equal(v[["ci_lower"]], theta - q * se, tolerance = 1e-12)
  expect_equal(v[["ci_upper"]], theta + q * se, tolerance = 1e-12)
  expect_equal(v[["slope"]], b, tolerance = 1e-12)
  expect_equal(v[["correlation"]], cor(d$y, d$x), tolerance = 1e-12)
  expect_equal(v[["proxy_bias"]], mean(d$x) - mean(d$y), tolerance = 1e-12)
  expect_equal(v[["se_ratio_target_only"]], se / (sd(d$y) / sqrt(n)),
               tolerance = 1e-12)
  expect_equal(c(v[["n_labeled"]], v[["n_pool"]]), c(n, M))
})

test_that("the regression estimate is the mean over the pool of the line fitted to the labeled units", {
  # Cochran's form of the regression estimator, ybar + b (xbar_pool - xbar),
  # is the average over the whole pool of the least squares line fitted to
  # the labeled units; lm() and predict() compute it here.
  d <- .two_phase_data(seed = 23)
  fit <- lm(y ~ x, data = data.frame(y = d$y, x = d$x))
  pool <- data.frame(x = c(d$x, d$xu))
  expect_equal(.v(ci_two_phase(d$y, d$x, d$xu))[["mean_Y"]],
               mean(predict(fit, newdata = pool)), tolerance = 1e-12)
})

test_that("the regression estimate is the prediction-powered form with proxy weight b N / M", {
  d <- .two_phase_data(seed = 3)
  n <- length(d$y); N <- length(d$xu)
  v <- .v(ci_two_phase(d$y, d$x, d$xu, conf_level = 0.90))
  lam <- v[["slope"]] * N / (n + N)
  expect_equal(v[["proxy_weight"]], lam, tolerance = 1e-12)
  expect_equal(v[["mean_Y"]],
               lam * mean(d$xu) + (mean(d$y) - lam * mean(d$x)),
               tolerance = 1e-12)
})

test_that("the prediction-powered interval has its closed form", {
  d <- .two_phase_data(seed = 5)
  n <- length(d$y); N <- length(d$xu)
  theta <- mean(d$xu) + mean(d$y) - mean(d$x)
  se <- sqrt(var(d$y - d$x) / n + var(d$xu) / N)
  q <- qt(0.95, n - 1)
  res <- ci_two_phase(d$y, d$x, d$xu, conf_level = 0.90,
                      method = "prediction_powered")
  v <- .v(res)
  expect_equal(v[["mean_Y"]], theta, tolerance = 1e-12)
  expect_equal(v[["se"]], se, tolerance = 1e-12)
  expect_equal(v[["df"]], n - 1)
  expect_equal(v[["ci_lower"]], theta - q * se, tolerance = 1e-12)
  expect_equal(v[["ci_upper"]], theta + q * se, tolerance = 1e-12)
  expect_equal(v[["proxy_weight"]], 1)
  expect_identical(attr(res, "method"), "prediction_powered")
})

test_that("the regression estimate is invariant to the proxy's location and scale", {
  d <- .two_phase_data(seed = 7)
  a <- .v(ci_two_phase(d$y, d$x, d$xu))
  b <- .v(ci_two_phase(d$y, 3 - 0.2 * d$x, 3 - 0.2 * d$xu))
  for (k in c("mean_Y", "se", "ci_lower", "ci_upper",
              "se_ratio_target_only"))
    expect_equal(a[[k]], b[[k]], tolerance = 1e-10, info = k)
  # The slope carries the proxy's units: it is rescaled by -1 / 0.2.
  expect_equal(b[["slope"]], -5 * a[["slope"]], tolerance = 1e-10)
})

test_that("the fast path returns the public interval", {
  skip_on_cran()
  for (seed in 1:4) {
    d <- .two_phase_data(n = 40 + 20 * seed, N = 300 * seed, seed = seed)
    for (cl in c(0.80, 0.95, 0.99)) {
      pub <- .v(ci_two_phase(d$y, d$x, d$xu, conf_level = cl))
      fast <- DMAR:::.ci_two_phase_fast(d$y, d$x, length(d$xu), mean(d$xu),
                                        sum((d$xu - mean(d$xu))^2), cl)
      for (k in names(fast))
        expect_equal(fast[[k]], pub[[k]], tolerance = 1e-12,
                     info = paste(seed, cl, k))
    }
  }
})

test_that("the table has the documented rows and classes", {
  d <- .two_phase_data()
  res <- ci_two_phase(d$y, d$x, d$xu)
  expect_s3_class(res, "dmar_tbl")
  expect_named(res, c("term", "value"))
  expect_type(res$value, "double")
  expect_identical(res$term,
                   c("mean_Y", "se", "df", "ci_lower", "ci_upper",
                     "slope", "proxy_weight", "correlation", "proxy_bias",
                     "se_ratio_target_only", "n_labeled", "n_pool"))
  # The estimate is named for its estimand, the mean of the target, and the
  # weight for what it is (Ken's ruling of 2026-11-11); the first 1.1.0
  # builds called them theta_hat and lambda.
  expect_false(any(c("theta_hat", "lambda") %in% res$term))
  pp <- ci_two_phase(d$y, d$x, d$xu, method = "prediction_powered")
  expect_identical(pp$term, res$term)
  expect_identical(attr(res, "method"), "regression")
  expect_equal(attr(res, "conf_level"), 0.95)
  td <- generics::tidy(res)
  expect_equal(td$conf_level[td$term == "mean_Y"], 0.95)
})

test_that("missing values are set aside unit by unit and counted out of the size rows", {
  d <- .two_phase_data(seed = 9)
  y <- d$y; x <- d$x; xu <- d$xu
  y[c(2, 9)] <- NA
  x[c(5, 9)] <- NA
  xu[c(1, 50, 51)] <- NA
  keep <- !is.na(y) & !is.na(x)
  with_na <- .v(ci_two_phase(y, x, xu))
  clean <- .v(ci_two_phase(y[keep], x[keep], xu[!is.na(xu)]))
  expect_equal(with_na, clean)
  expect_equal(with_na[["n_labeled"]], length(y) - 3)
  expect_equal(with_na[["n_pool"]], length(y) - 3 + length(xu) - 3)
})

test_that("logical measurements are taken as 0 and 1", {
  set.seed(21)
  x <- rbinom(200, 1, 0.6) == 1
  y <- ifelse(runif(200) < 0.8, x, rbinom(200, 1, 0.6) == 1)
  xu <- rbinom(2000, 1, 0.6) == 1
  expect_equal(.v(ci_two_phase(y, x, xu)),
               .v(ci_two_phase(as.numeric(y), as.numeric(x), as.numeric(xu))))
})

test_that("the cluster standard error sums each unit's contribution within clusters", {
  d <- .two_phase_data(n = 120, N = 600, seed = 13)
  n <- length(d$y); N <- length(d$xu); M <- n + N
  cl_lab <- rep(1:30, each = 4)
  cl_unl <- rep(1:30, each = 20)
  res <- .v(ci_two_phase(d$y, d$x, d$xu, cluster = cl_lab,
                         cluster_unlabeled = cl_unl))
  # The linearization, written out here: a labeled unit contributes its
  # residual over n and every unit its share of the pool mean.
  b <- cov(d$y, d$x) / var(d$x)
  xbar_M <- mean(c(d$x, d$xu))
  e <- resid(lm(d$y ~ d$x))
  u <- c(e / n + b * (d$x - xbar_M) / M, b * (d$xu - xbar_M) / M)
  tot <- tapply(u, c(cl_lab, cl_unl), sum)
  se <- sqrt(30 / 29 * sum(tot^2))
  expect_equal(res[["se"]], se, tolerance = 1e-12)
  expect_equal(res[["df"]], 29)
  t_y <- tapply((d$y - mean(d$y)) / n, cl_lab, sum)
  expect_equal(res[["se_ratio_target_only"]], se / sqrt(30 / 29 * sum(t_y^2)),
               tolerance = 1e-12)
  # The estimate itself does not depend on the clustering.
  expect_equal(res[["mean_Y"]], .v(ci_two_phase(d$y, d$x, d$xu))[["mean_Y"]])
})

test_that("with one unit per cluster the cluster standard error is the unclustered sandwich", {
  # An identity that does not go through the influence values: least
  # squares residuals sum to zero and are orthogonal to the proxy, so with
  # singleton clusters the sum of squared contributions is
  # sum(e^2) / n^2 + b^2 SS_pool / M^2, and the factor is n / (n - 1), the
  # labeled units being the labeled clusters.
  d <- .two_phase_data(n = 80, N = 400, seed = 17)
  n <- length(d$y); N <- length(d$xu); M <- n + N
  fit <- lm(d$y ~ d$x)
  b <- unname(coef(fit)[2])
  x_all <- c(d$x, d$xu)
  ss_pool <- sum((x_all - mean(x_all))^2)
  res <- .v(ci_two_phase(d$y, d$x, d$xu, cluster = seq_len(n),
                         cluster_unlabeled = n + seq_len(N)))
  expect_equal(res[["se"]],
               sqrt(n / (n - 1) * (sum(resid(fit)^2) / n^2 + b^2 * ss_pool / M^2)),
               tolerance = 1e-12)
  expect_equal(res[["df"]], n - 1)
  # For the prediction-powered estimate the labeled part is then exactly the
  # unclustered s^2_{Y - X} / n.
  pp <- .v(ci_two_phase(d$y, d$x, d$xu, method = "prediction_powered",
                        cluster = seq_len(n), cluster_unlabeled = n + seq_len(N)))
  expect_equal(pp[["se"]],
               sqrt(var(d$y - d$x) / n +
                      n / (n - 1) * (N - 1) * var(d$xu) / N^2),
               tolerance = 1e-12)
  expect_equal(pp[["df"]], n - 1)
})

test_that("the cluster degrees of freedom count the clusters among the labeled units", {
  # Labels from 6 clusters of a pool whose other units fall in 94 clusters
  # that hold only unlabeled units: the interval rests on the 6.
  d <- .two_phase_data(n = 60, N = 940, seed = 19)
  cl_lab <- rep(1:6, each = 10)
  cl_unl <- c(rep(1:6, each = 40), rep(7:100, length.out = 700))
  for (method in c("regression", "prediction_powered")) {
    res <- .v(ci_two_phase(d$y, d$x, d$xu, method = method, cluster = cl_lab,
                           cluster_unlabeled = cl_unl))
    expect_equal(res[["df"]], 5, info = method)
    expect_equal(res[["ci_upper"]] - res[["mean_Y"]],
                 qt(0.975, 5) * res[["se"]], tolerance = 1e-12, info = method)
    # Splitting the unlabeled-only units into more clusters changes neither
    # the degrees of freedom nor the factor.
    alt <- .v(ci_two_phase(d$y, d$x, d$xu, method = method, cluster = cl_lab,
                           cluster_unlabeled = c(cl_unl[1:240],
                                                 1000 + seq_len(700))))
    expect_equal(alt[["df"]], 5, info = method)
  }
})

test_that("the cluster standard error tracks the sampling variability of the estimate", {
  skip_on_cran()
  # Pools of 60 clusters of 25 units, with cluster effects in the proxy and
  # in the target's residual, from a population whose target mean is 0; the
  # labeled units are a simple random subsample of 150. Whatever the
  # linearization, the mean squared standard error must match the variance
  # of the estimate over replications; the unclustered standard error, which
  # ignores the clusters, falls well short.
  set.seed(42)
  R <- 1500
  q <- rep(1:60, each = 25)
  est <- se <- se_plain <- numeric(R)
  for (r in seq_len(R)) {
    a <- rnorm(60); u <- rnorm(60)
    x <- 0.7 * a[q] + rnorm(1500)
    y <- 0.6 * x + 0.5 * u[q] + rnorm(1500)
    lab <- sample(1500, 150)
    v <- .v(ci_two_phase(y[lab], x[lab], x[-lab], cluster = q[lab],
                         cluster_unlabeled = q[-lab]))
    est[r] <- v[["mean_Y"]]
    se[r] <- v[["se"]]
    se_plain[r] <- .v(ci_two_phase(y[lab], x[lab], x[-lab]))[["se"]]
  }
  # The relative Monte Carlo SE of var(est) at R = 1500 is about .037.
  expect_lt(abs(mean(se^2) / var(est) - 1), 0.12)
  expect_lt(mean(se_plain^2) / var(est), 0.75)
})

test_that("labels from a few clusters of a large pool keep the interval near nominal", {
  skip_on_cran()
  # 100 clusters of 20 units, every unit of 6 clusters labeled. With the
  # degrees of freedom of all 100 clusters the 95% interval covered about
  # .85 of the time; on the labeled clusters' 5 it is near .95.
  set.seed(5)
  R <- 1000
  q <- rep(1:100, each = 20)
  hit <- logical(R)
  for (r in seq_len(R)) {
    a <- rnorm(100); u <- rnorm(100)
    x <- 0.7 * a[q] + rnorm(2000)
    y <- 0.6 * x + 0.5 * u[q] + rnorm(2000)
    lab <- which(q %in% sample(100, 6))
    v <- .v(ci_two_phase(y[lab], x[lab], x[-lab], cluster = q[lab],
                         cluster_unlabeled = q[-lab]))
    hit[r] <- v[["ci_lower"]] <= 0 && 0 <= v[["ci_upper"]]
  }
  # Monte Carlo SE of a .95 coverage at R = 1000 is about .007.
  expect_gt(mean(hit), 0.91)
  expect_lt(mean(hit), 0.98)
})

test_that("a target constant on the labeled units is refused rather than given a zero-width interval", {
  d <- .two_phase_data()
  y1 <- rep(1, length(d$y))
  expect_error(ci_two_phase(y1, d$x, d$xu), "target takes a single value")
  expect_error(ci_two_phase(y1 == 1, d$x, d$xu), "target takes a single value")
  expect_error(ci_two_phase(y1, d$x, d$xu, cluster = rep(1:30, 5),
                            cluster_unlabeled = rep(1:30, length.out = 1000)),
               "target takes a single value")
  # The prediction-powered interval does not rest on the residual variation,
  # and stays defined.
  v <- .v(ci_two_phase(y1, d$x, d$xu, method = "prediction_powered"))
  expect_gt(v[["se"]], 0)
})

test_that("malformed inputs are refused with a message naming the problem", {
  d <- .two_phase_data()
  expect_error(ci_two_phase(d$y, d$x[-1], d$xu), "same length")
  expect_error(ci_two_phase(as.character(d$y), d$x, d$xu), "numeric")
  expect_error(ci_two_phase(d$y, d$x, d$xu, conf_level = 1.2), "conf_level")
  expect_error(ci_two_phase(d$y[1:2], d$x[1:2], d$xu), "At least 3 labeled")
  expect_error(ci_two_phase(d$y, d$x, numeric(0)), "unlabeled")
  expect_error(ci_two_phase(d$y, d$x, d$xu[1], method = "prediction_powered"),
               "At least 2 unlabeled")
  expect_error(ci_two_phase(d$y, rep(1, length(d$y)), d$xu), "single value")
  expect_error(ci_two_phase(d$y, d$x, c(d$xu, Inf)), "infinite")
  expect_error(ci_two_phase(d$y, d$x, d$xu, cluster = seq_along(d$y)),
               "both")
  expect_error(ci_two_phase(d$y, d$x, d$xu, cluster = 1:3,
                            cluster_unlabeled = seq_along(d$xu)),
               "one entry per labeled unit")
  expect_error(ci_two_phase(d$y, d$x, d$xu, cluster = rep(1, length(d$y)),
                            cluster_unlabeled = seq_along(d$xu)),
               "fewer than two clusters")
})

test_that("a constant proxy leaves the prediction-powered interval defined", {
  d <- .two_phase_data()
  v <- .v(ci_two_phase(d$y, rep(50, length(d$y)), d$xu,
                       method = "prediction_powered"))
  expect_true(is.na(v[["slope"]]))
  expect_true(is.na(v[["correlation"]]))
  expect_true(is.finite(v[["mean_Y"]]))
})
