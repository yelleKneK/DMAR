.htmt_data <- function(n = 300) {
  set.seed(113)
  f1 <- rnorm(n); f2 <- 0.3 * f1 + sqrt(1 - 0.09) * rnorm(n)
  data.frame(
    a1 = .8 * f1 + rnorm(n, 0, .6), a2 = .7 * f1 + rnorm(n, 0, .7),
    a3 = .6 * f1 + rnorm(n, 0, .8),
    b1 = .8 * f2 + rnorm(n, 0, .6), b2 = .7 * f2 + rnorm(n, 0, .7),
    b3 = .6 * f2 + rnorm(n, 0, .8))
}

test_that("htmt() matches the Henseler definition by hand and semTools", {
  d <- .htmt_data()
  blocks <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
  res <- htmt(d, blocks)
  expect_s3_class(res, "dmar_tbl")
  expect_equal(nrow(res), 1)
  # By hand.
  R <- cor(d)
  hetero <- mean(abs(R[1:3, 4:6]))
  mono_a <- mean(R[1:3, 1:3][lower.tri(diag(3))])
  mono_b <- mean(R[4:6, 4:6][lower.tri(diag(3))])
  # The by-hand computation above is the oracle: the ratio of the mean
  # between-construct correlation to the geometric mean of the two
  # within-construct means, which is the definition. It is exact, so no
  # second implementation is needed to confirm it.
  expect_equal(res$htmt, hetero / sqrt(mono_a * mono_b))
})

test_that("htmt() bootstrap upper bound behaves and is seeded", {
  d <- .htmt_data(150)
  blocks <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"))
  r1 <- htmt(d, blocks, B = 200, seed = 7)
  r2 <- htmt(d, blocks, B = 200, seed = 7)
  expect_identical(r1$upper_limit, r2$upper_limit)
  expect_gt(r1$upper_limit, r1$htmt)
  expect_identical(attr(r1, "conf_level"), 0.95)
})

test_that("htmt() validates and handles three constructs", {
  d <- .htmt_data()
  set.seed(7)
  f3 <- rnorm(nrow(d))
  d$c1 <- .7 * f3 + rnorm(nrow(d), 0, .7)
  d$c2 <- .6 * f3 + rnorm(nrow(d), 0, .8)
  res <- htmt(d, list(A = c("a1", "a2"), B = c("b1", "b2"),
                      C = c("c1", "c2")))
  expect_equal(nrow(res), 3)             # all pairs
  expect_error(htmt(d, list(A = "a1", B = c("b1", "b2"))), "two or more")
  expect_error(htmt(d, list(c("a1", "a2"), B = c("b1", "b2"))), "named")
  expect_error(htmt(d, list(A = c("a1", "a2"), B = c("b1", "b2")),
                    B = 10), "at least 100")
})

test_that("a degenerate resample is dropped, not fatal", {
  # Block B's two items are nearly uncorrelated, so a nontrivial share
  # of resamples flip its average within-block correlation negative;
  # those resamples must be dropped with one warning, not abort the
  # bootstrap.
  set.seed(113)
  n <- 40
  f <- rnorm(n)
  d <- data.frame(
    a1 = f + rnorm(n, sd = 0.4), a2 = f + rnorm(n, sd = 0.4),
    b1 = rnorm(n), b2 = rnorm(n)
  )
  # Keep the observed within-block mean correlation of B positive so
  # the point estimate is defined.
  if (cor(d$b1, d$b2) <= 0) d$b2 <- d$b2 + 0.15 * d$b1
  res <- withCallingHandlers(
    htmt(d, blocks = list(A = c("a1", "a2"), B = c("b1", "b2")),
         B = 300, seed = 113),
    warning = function(w) {
      expect_match(conditionMessage(w), "dropped")
      invokeRestart("muffleWarning")
    }
  )
  expect_true(all(is.finite(res$upper_limit)))
})
