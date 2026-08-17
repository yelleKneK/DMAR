test_that("ss_seq_c() pilot mode returns the pilot size and an allocation
           that sums to one", {
  res <- ss_seq_c(c_weights = c(1, -1), half_width = 2.5, pilot = TRUE)
  expect_true("pilot_n_per_group" %in% res$term)
  expect_equal(res$value[res$term == "pilot_n_per_group"], 10)
  alloc <- res$value[grepl("^allocation_", res$term)]
  expect_equal(sum(alloc), 1, tolerance = 1e-12)
  expect_equal(alloc, c(0.5, 0.5), tolerance = 1e-12)
})

test_that("ss_seq_c() stops when the criterion is met and continues when it
           is not", {
  # Precise: q * sqrt(sum(c^2 s^2 / n)) well under the target.
  precise <- ss_seq_c(c_weights = c(1, -1), half_width = 2.5,
                      s = 10, n = c(500, 500))
  expect_equal(precise$value[precise$term == "stop"], 1)

  imprecise <- ss_seq_c(c_weights = c(1, -1), half_width = 2.5,
                        s = 15.4, n = c(60, 60))
  expect_equal(imprecise$value[imprecise$term == "stop"], 0)
  expect_gt(imprecise$value[imprecise$term == "N_projected"],
            imprecise$value[imprecise$term == "N_current"])
})

test_that("ss_seq_c() half_width_current matches the direct computation", {
  res <- ss_seq_c(c_weights = c(1, -1), half_width = 2.5,
                  s = 15.4, n = c(60, 60))
  nu <- 120 - 2
  expected <- qt(0.95, nu) * 15.4 * sqrt(1 / 60 + 1 / 60)
  expect_equal(res$value[res$term == "half_width_current"], expected,
               tolerance = 1e-12)
})

test_that("ss_seq_c() normal quantile stops earlier than the t quantile", {
  # At a boundary case the t criterion is stricter (larger quantile).
  s <- 15.4; n <- c(60, 60)
  ht <- ss_seq_c(c_weights = c(1, -1), half_width = 2.5, s = s, n = n,
                 quantile = "t")
  hz <- ss_seq_c(c_weights = c(1, -1), half_width = 2.5, s = s, n = n,
                 quantile = "normal")
  expect_gt(ht$value[ht$term == "half_width_current"],
            hz$value[hz$term == "half_width_current"])
})

test_that("ss_seq_c() cost-optimal allocation follows |c| * s / sqrt(cost)", {
  res <- ss_seq_c(c_weights = c(1, -1), half_width = 2.5,
                  s = c(15.4, 15.4), n = c(60, 60), cost = c(1, 4))
  alloc <- res$value[grepl("^allocation_", res$term)]
  raw <- abs(c(1, -1)) * 15.4 / sqrt(c(1, 4))
  expect_equal(alloc, raw / sum(raw), tolerance = 1e-12)
  # The expensive group receives the smaller share.
  expect_lt(alloc[2], alloc[1])
})

test_that("ss_seq_c() zero-weight groups get no allocation", {
  res <- ss_seq_c(c_weights = c(0.5, 0.5, -1, 0), half_width = 2,
                  pilot = TRUE)
  alloc <- res$value[grepl("^allocation_", res$term)]
  expect_equal(alloc[4], 0)
  expect_equal(sum(alloc), 1, tolerance = 1e-12)
})

test_that("ss_seq_c() validates its inputs", {
  expect_error(ss_seq_c(c_weights = c(1, -1), half_width = 2.5),
               "supply the current 's' and 'n'")
  expect_error(ss_seq_c(c_weights = c(2, -2), half_width = 2.5,
                        s = 10, n = c(50, 50)),
               "positive weights must sum to 1")
  expect_error(ss_seq_c(c_weights = c(1, -1), half_width = -1,
                        s = 10, n = c(50, 50)),
               "half_width")
})
