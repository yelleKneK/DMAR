test_that("pairwise_within() wide-format returns documented columns", {
  set.seed(113)
  Y <- matrix(rnorm(40), 10, 4) +
       matrix(rep(c(0, 0.3, 0.6, 0.9), each = 10), 10, 4)
  colnames(Y) <- paste0("T", 1:4)
  res <- pairwise_within(Y)
  expect_setequal(colnames(res),
                  c("contrast", "mean_difference", "sd_difference",
                    "t_statistic", "df", "p_value", "p_adjusted",
                    "lower_limit", "upper_limit", "n_pairs"))
  expect_equal(nrow(res), 6)  # choose(4, 2)
})

test_that("pairwise_within() Holm adjusts p-values upward", {
  set.seed(113)
  Y <- matrix(rnorm(50), 10, 5)
  colnames(Y) <- paste0("T", 1:5)
  res <- pairwise_within(Y, adjust = "holm")
  expect_true(all(res$p_adjusted >= res$p_value - 1e-12))
})

test_that("pairwise_within() long-format equals wide-format", {
  set.seed(113)
  n <- 8; k <- 3
  Y <- matrix(rnorm(n * k), n, k); colnames(Y) <- paste0("T", 1:k)
  long <- data.frame(subject = factor(rep(1:n, k)),
                     time    = factor(rep(paste0("T", 1:k), each = n)),
                     y       = as.vector(Y))
  wide_res <- pairwise_within(Y)
  long_res <- pairwise_within(long, subject = "subject",
                              condition = "time", outcome = "y")
  expect_equal(sort(wide_res$mean_difference),
               sort(long_res$mean_difference), tolerance = 1e-10)
})
