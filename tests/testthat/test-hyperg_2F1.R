## .hyperg_2F1 -- base-R Gauss hypergeometric 2F1, the replacement for
## gsl::hyperg_2F1 across the R^2 moment functions.

test_that(".hyperg_2F1 reproduces the closed form 2F1(1,1;2;x) = -log(1-x)/x", {
  # 2F1(1,1;2;x) has the elementary closed form -log(1-x)/x, an exact check
  # that needs no external special-function library.
  for (x in c(0.1, 0.3, 0.5, 0.8, 0.95, 0.999)) {
    expect_equal(.hyperg_2F1(1, 1, 2, x), -log(1 - x) / x, tolerance = 1e-11)
  }
})

test_that(".hyperg_2F1 returns 1 at x = 0", {
  expect_identical(.hyperg_2F1(1, 1, 5.5, 0), 1)
  expect_identical(.hyperg_2F1(2, 2, 6.5, 0), 1)
})

test_that(".hyperg_2F1 matches independently verified values, including x -> 1", {
  # Reference values from the independent Euler integral representation (and a
  # high-precision series); these are the TRUE values. gsl is correct for the
  # first three but loses ~5e-3 of accuracy on the last one as x -> 1, which is
  # one reason DMAR computes the 2F1 itself.
  expect_equal(.hyperg_2F1(1, 1, 5.5, 0.30),   1.060286213596, tolerance = 1e-9)
  expect_equal(.hyperg_2F1(2, 2, 6.5, 0.30),   1.225659361603, tolerance = 1e-9)
  expect_equal(.hyperg_2F1(1, 1, 1.5, 0.999), 48.697129242100, tolerance = 1e-6)
  expect_equal(.hyperg_2F1(2, 2, 9.0, 0.999),  1.864802792578, tolerance = 1e-9)
})

test_that(".hyperg_2F1 agrees with gsl away from the x -> 1 boundary", {
  # Pinned from gsl::hyperg_2F1 (gsl 2.1.8, 2026-08-09); live comparison in
  # tools/oracle_checks.R. The 60 evaluation points are the draws the loop
  # below produces under set.seed(113), so only the reference values are
  # stored.
  oracle <- c(1.002579971343044, 1.005651033417876, 1.00237474688911,
              1.000715244297069, 1.00796045886712,  1.08897172769139,
              1.008671576648846, 1.008715759272523, 1.002943319541831,
              1.001163781448927, 1.002426084222118, 1.008410284351776,
              1.011431110411731, 1.001832846321943, 1.00909255189971,
              1.002397467167829, 1.00138818418789,  1.000694126916134,
              1.002507886286649, 1.000708506332314, 1.00152847294778,
              1.004569504280686, 1.000468501959024, 1.000172632104875,
              1.023804704372723, 1.000933453428324, 1.001440913844438,
              1.007839678260801, 1.001074247462407, 1.010632919914832,
              1.004166722349747, 1.002264168965575, 1.00058910046866,
              1.003868301147627, 1.001481501643864, 1.030752265237992,
              1.00810375136601,  1.002765427501132, 1.000935769191928,
              1.003996224017043, 1.014770773250152, 1.005842012513908,
              1.004246839404026, 1.000018992370147, 1.00356620277097,
              1.00364954658357,  1.002864858881864, 1.003675340021205,
              1.003742653640696, 1.005787885789253, 1.000167600523812,
              1.019969364694074, 1.000346681287388, 1.000517408682314,
              1.00502957033638,  1.002082461563719, 1.028166576799452,
              1.001998742060819, 1.006017108617989, 1.001449888147296)
  set.seed(113)
  for (i in 1:60) {
    a  <- sample(1:2, 1)
    N  <- sample(8:1000, 1)
    x  <- runif(1, 0, 0.9)            # away from 1, where gsl is itself accurate
    cc <- (N + 1) / 2
    expect_equal(.hyperg_2F1(a, a, cc, x), oracle[i], tolerance = 1e-9)
  }
})
