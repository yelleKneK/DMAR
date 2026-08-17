test_that("ss_aipe_smd() returns a tidy data frame led by necessary_n_per_group", {
  result <- ss_aipe_smd(delta = 0.5, conf_level = 0.95, width = 0.4)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("term", "value"))
  expect_equal(result$term[1], "necessary_n_per_group")
  expect_true(all(c("supposed_smd", "width") %in% result$term))
  expect_false("delta" %in% result$term)   # denoted as supposed, not bare delta
})

test_that("ss_aipe_smd() returns a positive integer-valued sample size", {
  r <- ss_aipe_smd(delta = 0.5, conf_level = 0.95, width = 0.4)
  n <- r$value[r$term == "necessary_n_per_group"]
  expect_gt(n, 0)
  expect_equal(n, as.integer(n))
})

test_that("ss_aipe_smd() requires a larger n for a narrower width", {
  n_of <- function(w) {
    r <- ss_aipe_smd(delta = 0.5, conf_level = 0.95, width = w)
    r$value[r$term == "necessary_n_per_group"]
  }
  expect_lt(n_of(0.5), n_of(0.3))
})

test_that("ss_aipe_smd() echoes the planning inputs as numeric rows", {
  r <- ss_aipe_smd(delta = 0.5, conf_level = 0.95, width = 0.4)
  expect_equal(r$value[r$term == "supposed_smd"], 0.5)
  expect_equal(r$value[r$term == "width"], 0.4)
  expect_false("assurance" %in% r$term)          # not supplied -> no row
  expect_type(r$value, "double")                 # value column stays numeric
  ra <- ss_aipe_smd(delta = 0.5, conf_level = 0.95, width = 0.4, assurance = 0.9)
  expect_equal(ra$value[ra$term == "assurance"], 0.9)
})

test_that("ss_aipe_smd() floors an impossibly wide target at the admissible minimum", {
  # A huge target width previously crashed the fixed-point search ("missing
  # value where TRUE/FALSE needed") because the central-t error df 2n - 2 went
  # negative. The smallest admissible design this planner returns is 4 per
  # group.
  r <- ss_aipe_smd(delta = 0.5, width = 100)
  expect_equal(r$value[r$term == "necessary_n_per_group"], 4)
})

test_that("ss_aipe_smd() validates its boundary inputs", {
  expect_error(ss_aipe_smd(delta = 0.5, width = 0), "width")
  expect_error(ss_aipe_smd(delta = 0.5, width = -1), "width")
  expect_error(ss_aipe_smd(delta = 0.5, width = 0.3, conf_level = 1.2),
               "conf_level")
  expect_error(ss_aipe_smd(delta = NA_real_, width = 0.3), "delta")
})
