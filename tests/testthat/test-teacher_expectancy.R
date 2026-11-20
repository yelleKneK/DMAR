test_that("teacher_expectancy matches the Raudenbush (1984) Table 1 summaries", {
  data(teacher_expectancy)
  te <- teacher_expectancy
  expect_equal(nrow(te), 19)
  expect_equal(sum(grepl("Pellegrini", te$author)), 2)
  expect_equal(range(te$weeks), c(0, 24))
  expect_equal(sum(te$testing == "individual"), 3)   # paper: n = 3
  expect_equal(sum(te$tester == "aware"), 10)        # paper: n = 10 aware
  expect_equal(sum(te$tester == "blind"), 9)         # paper: n = 9 blind

  # Study-level reconstruction (Pellegrini & Hicks merged at d = .52):
  # the paper reports M = 0.11, SD = 0.20, range 0.55 to -0.13, five
  # studies significant (three at .05, two at .01), and r(d, weeks) = -.55.
  d18 <- append(te$d[-c(4, 5)], 0.52, after = 3)
  p18 <- append(te$p_one_tailed[-c(4, 5)], .010, after = 3)
  wk18 <- append(te$weeks[-c(4, 5)], 0, after = 3)
  expect_equal(round(mean(d18), 2), 0.11)
  expect_equal(round(sd(d18), 2), 0.20)
  expect_equal(max(d18), 0.55)
  expect_equal(min(d18), -0.13)
  expect_equal(sum(p18 <= .05), 5)
  expect_equal(sum(p18 <= .011), 2)   # .002 and .010; the paper counts "two at .01"
  expect_equal(cor(d18, wk18), -0.55, tolerance = 0.005)
  # And his linearized correlation: r = -.77 after the reciprocal
  # transformations x -> -1/(2 + x), y -> -1/(1 + y).
  expect_equal(cor(-1 / (1 + d18), -1 / (2 + wk18)), -0.77,
               tolerance = 0.015)
})

test_that("the prior-contact contrast stated on ?teacher_expectancy", {
  # The help page says that, over the 18 studies, the contrast with weights
  # 1 / (weeks + 2), centered to sum to zero, gives z = 2.76, one-tailed
  # p = .003. Oracle: the Rosenthal-Rubin contrast written out here,
  # z = sum(lambda d) / sqrt(sum(lambda^2 v)), with the large-sample variance
  # of d, v = (n_e + n_c) / (n_e n_c) + d^2 / (2 (n_e + n_c)), and the
  # Pellegrini and Hicks conditions merged as on the page (d = .52, 0 weeks,
  # 22 and 22 pupils), as in the teacher expectancy vignette.
  te <- teacher_expectancy[-c(4, 5), ]
  d18 <- append(te$d, 0.52, after = 3)
  wk18 <- append(te$weeks, 0, after = 3)
  ne18 <- append(te$n_experimental, 22, after = 3)
  nc18 <- append(te$n_control, 22, after = 3)
  v18 <- (ne18 + nc18) / (ne18 * nc18) + d18^2 / (2 * (ne18 + nc18))
  lambda <- 1 / (wk18 + 2)
  lambda <- lambda - mean(lambda)
  z <- sum(lambda * d18) / sqrt(sum(lambda^2 * v18))
  expect_equal(round(z, 2), 2.76)
  expect_equal(round(pnorm(-z), 3), 0.003)
  expect_equal(round(cor(d18, wk18), 2), -0.55)
  out <- suppressMessages(meta_contrast(d18, v18, weights = 1 / (wk18 + 2)))
  expect_equal(out$value[out$term == "test_statistic"], z, tolerance = 1e-12)
})
