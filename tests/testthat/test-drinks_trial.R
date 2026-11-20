test_that("drinks_trial has the documented shape", {
  data(drinks_trial, package = "DMAR", envir = environment())
  expect_s3_class(drinks_trial, "data.frame")
  expect_equal(dim(drinks_trial), c(88L, 5L))
  expect_true(all(c("id", "cohort", "treatment",
                    "drinks_per_week", "log_drinks") %in%
                  names(drinks_trial)))
})

test_that("cohort and treatment factor levels are as documented", {
  data(drinks_trial, package = "DMAR", envir = environment())
  expect_equal(levels(drinks_trial$cohort),
               c("1", "2"))
  expect_equal(levels(drinks_trial$treatment),
               c("Standard", "CRA", "CRA + Disulfiram"))
})

test_that("per-cell sample sizes reconcile with the published paper", {
  data(drinks_trial, package = "DMAR", envir = environment())
  ct <- with(drinks_trial, table(cohort, treatment))
  expect_equal(as.integer(ct["1", "Standard"]),         17L)
  expect_equal(as.integer(ct["1", "CRA"]),              15L)
  expect_equal(as.integer(ct["1", "CRA + Disulfiram"]), 19L)
  expect_equal(as.integer(ct["2", "Standard"]),         20L)
  expect_equal(as.integer(ct["2", "CRA"]),              17L)
  # Cohort 2 had no disulfiram cell by design.
  expect_equal(as.integer(ct["2", "CRA + Disulfiram"]), 0L)
})

test_that("marginal sample sizes match the documentation", {
  data(drinks_trial, package = "DMAR", envir = environment())
  expect_equal(sum(drinks_trial$cohort == "1"), 51L)
  expect_equal(sum(drinks_trial$cohort == "2"), 37L)
  expect_equal(sum(drinks_trial$treatment == "Standard"),         37L)
  expect_equal(sum(drinks_trial$treatment == "CRA"),              32L)
  expect_equal(sum(drinks_trial$treatment == "CRA + Disulfiram"), 19L)
})

test_that("drinks_per_week and log_drinks are consistent and on the expected scales", {
  data(drinks_trial, package = "DMAR", envir = environment())
  expect_equal(min(drinks_trial$drinks_per_week), 0)
  expect_true(max(drinks_trial$drinks_per_week) > 600)
  expect_equal(drinks_trial$log_drinks,
               log10(drinks_trial$drinks_per_week + 1),
               tolerance = 1e-4)
})
