# Integrity anchors for the shipped datasets that lack a dedicated
# test file. Dimensions, column names, and simple aggregates computed
# from the shipped objects; any change to data/ that alters these is a
# substantive change and should be deliberate.

test_that("diagnosis_agreement matches its documented structure", {
  data(diagnosis_agreement, envir = environment())
  expect_identical(dim(diagnosis_agreement), c(9L, 6L))
  expect_identical(names(diagnosis_agreement),
                   c("judge_b", "judge_a", "frequency",
                     "disagreement_weight", "observed_proportion",
                     "expected_proportion"))
  # Cohen (1968), Table 1: 200 cases in the 3 x 3 agreement table.
  expect_equal(sum(diagnosis_agreement$frequency), 200)
  expect_equal(sum(diagnosis_agreement$observed_proportion), 1,
               tolerance = 1e-8)
})

test_that("pygmalion matches its documented structure", {
  data(pygmalion, envir = environment())
  expect_identical(dim(pygmalion), c(310L, 6L))
  expect_identical(names(pygmalion),
                   c("grade", "treatment", "iq_pre", "iq_4", "iq_8",
                     "iq_gain"))
  expect_equal(mean(pygmalion$iq_pre), 98.464516, tolerance = 1e-6)
  expect_equal(mean(pygmalion$iq_gain), 9.432258, tolerance = 1e-6)
  expect_equal(pygmalion$iq_gain, pygmalion$iq_8 - pygmalion$iq_pre,
               tolerance = 1e-8)
})

test_that("test_market matches its documented structure", {
  data(test_market, envir = environment())
  expect_identical(dim(test_market), c(24L, 4L))
  expect_identical(names(test_market),
                   c("panel", "block", "brand_movement",
                     "category_movement"))
  expect_equal(mean(test_market$brand_movement), 4.220833,
               tolerance = 1e-6)
  expect_equal(mean(test_market$category_movement), 12.196667,
               tolerance = 1e-6)
})
