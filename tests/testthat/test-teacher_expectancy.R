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
