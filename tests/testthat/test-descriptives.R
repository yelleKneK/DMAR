test_that("descriptives() returns a list with 'descriptives' and 'correlations'", {
  result <- descriptives(attitude)
  expect_type(result, "list")
  expect_named(result, c("descriptives", "correlations"))
})

test_that("descriptives() produces one row per variable with expected columns", {
  result <- descriptives(mtcars)
  expect_s3_class(result$descriptives, "data.frame")
  expect_equal(nrow(result$descriptives), ncol(mtcars))
  expected_cols <- c("variable", "type", "n", "n_missing", "prop_missing",
                     "mean", "median", "sd", "min", "max", "q25", "q75",
                     "skewness", "kurtosis")
  expect_equal(names(result$descriptives), expected_cols)
})

test_that("descriptives() correlations is NULL by default", {
  result <- descriptives(attitude)
  expect_null(result$correlations)
})

test_that("descriptives() correlations is a matrix matching cor() when requested", {
  result <- descriptives(attitude, correlations = TRUE)
  expect_true(is.matrix(result$correlations))
  expect_equal(dim(result$correlations), c(ncol(attitude), ncol(attitude)))
  expect_equal(result$correlations,
               stats::cor(attitude, use = "pairwise.complete.obs"))
})

test_that("descriptives() numeric summaries match base R", {
  result <- descriptives(mtcars)
  expect_equal(result$descriptives$mean, unname(colMeans(mtcars)))
  expect_equal(result$descriptives$sd,   unname(apply(mtcars, 2, stats::sd)))
  expect_equal(result$descriptives$min,  unname(apply(mtcars, 2, min)))
  expect_equal(result$descriptives$max,  unname(apply(mtcars, 2, max)))
})

test_that("descriptives() counts missing values correctly", {
  x <- data.frame(a = c(1, 2, NA, 4), b = c(1, 2, 3, 4))
  result <- descriptives(x)
  expect_equal(result$descriptives$n,            c(3, 4))
  expect_equal(result$descriptives$n_missing,    c(1, 0))
  expect_equal(result$descriptives$prop_missing, c(0.25, 0))
})

test_that("descriptives() listwise deletion uses complete cases", {
  x <- data.frame(a = c(1, 2, NA, 4), b = c(1, NA, 3, 4))
  result <- descriptives(x, listwise = TRUE)
  # Only rows 1 and 4 are complete, so N = 2 for both variables.
  expect_equal(result$descriptives$n, c(2L, 2L))
  expect_equal(result$descriptives$n_missing, c(0L, 0L))
})

test_that("descriptives() factor variables get type='factor' with NA numeric stats", {
  result <- descriptives(iris)
  row_species <- result$descriptives[result$descriptives$variable == "Species", ]
  expect_equal(row_species$type, "factor")
  expect_true(is.na(row_species$mean))
  expect_true(is.na(row_species$sd))
  expect_equal(row_species$n, 150L)
})

test_that("descriptives() skewness is zero for a symmetric vector", {
  result <- descriptives(data.frame(v = c(1, 2, 3, 4, 5)))
  expect_equal(result$descriptives$skewness, 0)
})

test_that("descriptives() kurtosis matches SAS/SPSS Type 2 formula by hand", {
  # For x = 1:5, the SAS/SPSS Type-2 excess kurtosis is -1.2.
  result <- descriptives(data.frame(v = c(1, 2, 3, 4, 5)))
  expect_equal(result$descriptives$kurtosis, -1.2)
})

test_that("descriptives() returns NA for skewness/kurtosis when sd = 0", {
  result <- descriptives(data.frame(v = rep(3, 10)))
  expect_true(is.na(result$descriptives$skewness))
  expect_true(is.na(result$descriptives$kurtosis))
})

test_that("descriptives() handles matrix input", {
  result <- descriptives(as.matrix(mtcars))
  expect_equal(nrow(result$descriptives), ncol(mtcars))
  expect_equal(result$descriptives$variable, colnames(mtcars))
})

test_that("descriptives() warns and returns NULL correlations when <2 numeric vars", {
  x <- data.frame(a = 1:5, f = factor(letters[1:5]))
  expect_warning(res <- descriptives(x, correlations = TRUE),
                 "at least two numeric")
  expect_null(res$correlations)
})

test_that("descriptives() errors informatively on bad input", {
  expect_error(descriptives("hello"),     "data frame")
  expect_error(descriptives(list(a = 1)), "data frame")
  expect_error(descriptives(data.frame()),"no columns")
})

test_that("descriptives() logicals are summarized as proportions via mean", {
  result <- descriptives(data.frame(l = c(TRUE, TRUE, FALSE, FALSE, TRUE)))
  expect_equal(result$descriptives$type, "logical")
  expect_equal(result$descriptives$mean, 0.6)
})
