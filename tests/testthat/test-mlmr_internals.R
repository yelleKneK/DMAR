## The saturated correlates syntax of mlmr_auxiliary_syntax() in
## R/mlmr_internals.R, shared by mlmr() and mlmr_mv(). The fits it leads
## to are tested in test-mlmr.R and test-mlmr_mv.R.

# The covariance lines of a syntax string, each as a sorted pair, so that
# "a ~~ b" and "b ~~ a" count as the same line.
.mlmr_cov_pairs <- function(syntax) {
  lines <- strsplit(syntax, "\n", fixed = TRUE)[[1L]]
  parts <- strsplit(lines, " ~~ ", fixed = TRUE)
  vapply(parts, function(p) paste(sort(p), collapse = " ~~ "),
         character(1L))
}

test_that("every predictor pair covaries once auxiliaries are present", {
  # A predictor that covaries with an auxiliary is no longer exogenous to
  # lavaan, whose automatic covariances among exogenous variables then
  # stop reaching it; before 1.1.0 the syntax omitted the predictor pairs,
  # and every covariance among the predictors was fixed at zero.
  s <- mlmr_auxiliary_syntax("y", c("x_1", "x_2"), "aux_1")
  expect_setequal(.mlmr_cov_pairs(s),
                  c("aux_1 ~~ y", "aux_1 ~~ x_1", "aux_1 ~~ x_2",
                    "x_1 ~~ x_2"))

  # J outcomes, K predictors, A auxiliaries: J * A residual covariances,
  # A * K auxiliary-predictor covariances, choose(A, 2) among the
  # auxiliaries, and choose(K, 2) among the predictors, each once.
  y <- c("y_1", "y_2")
  x <- paste0("x_", 1:4)
  a <- paste0("aux_", 1:3)
  pairs <- .mlmr_cov_pairs(mlmr_auxiliary_syntax(y, x, a))
  expect_false(anyDuplicated(pairs) > 0L)
  expect_length(pairs, 2 * 3 + 3 * 4 + choose(3, 2) + choose(4, 2))
  expected_x <- apply(utils::combn(x, 2L), 2L, paste, collapse = " ~~ ")
  expect_true(all(expected_x %in% pairs))
  expected_a <- apply(utils::combn(a, 2L), 2L, paste, collapse = " ~~ ")
  expect_true(all(expected_a %in% pairs))
  expected_y <- .mlmr_cov_pairs(paste(outer(y, a, paste, sep = " ~~ "),
                                      collapse = "\n"))
  expect_true(all(expected_y %in% pairs))
  # No line names an outcome with a predictor (that would be a residual
  # covariance the regression does not have) or two outcomes (mlmr_mv()
  # writes the residual covariances of the outcomes itself).
  expect_false(any(grepl("^x_.* ~~ y_|^y_.* ~~ x_|^y_.* ~~ y_", pairs)))
})

test_that("one predictor, or none, adds no predictor covariance", {
  # The syntax of 1.0.x for these models, which the fix leaves alone.
  expect_identical(mlmr_auxiliary_syntax("y", "x_1", "aux_1"),
                   "y ~~ aux_1\naux_1 ~~ x_1")
  expect_identical(mlmr_auxiliary_syntax("y", character(0),
                                         c("aux_1", "aux_2")),
                   "y ~~ aux_1\ny ~~ aux_2\naux_1 ~~ aux_2")
})
