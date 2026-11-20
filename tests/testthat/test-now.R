## now() -- friendly current-date utility and lightweight stopwatch.

test_that("now() returns a dmar_now object that prints with month / day / year and a time block", {
  out <- now()
  expect_s3_class(out, "dmar_now")
  expect_s3_class(out, "POSIXct")
  printed <- format(out)
  expect_type(printed, "character")
  expect_match(printed, format(Sys.Date(), "%Y"))
  expect_match(printed, "\\(\\d+:\\d+ [AP]M\\)$")
})


test_that("now(time = FALSE) drops the time block", {
  out <- now(time = FALSE)
  expect_s3_class(out, "dmar_now")
  printed <- format(out)
  expect_match(printed, format(Sys.Date(), "%Y"))
  expect_false(grepl("[AP]M", printed))
})


test_that("subtraction of two now() captures returns a difftime", {
  start <- now()
  Sys.sleep(0.05)
  end <- now()
  elapsed <- end - start
  expect_s3_class(elapsed, "difftime")
  expect_true(as.numeric(elapsed, units = "secs") >= 0)
  expect_true(as.numeric(elapsed, units = "secs") < 5)
})


test_that("now(tidy = TRUE) returns a numeric value column with non-numeric metadata on attributes", {
  out <- now(tidy = TRUE)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("term", "value"))
  expect_true(is.numeric(out$value))
  expect_setequal(out$term, c("day", "year", "hour", "minute"))
  expect_true(is.character(attr(out, "month")))
  expect_true(is.character(attr(out, "am_pm")))
})


test_that("now(time = FALSE, tidy = TRUE) returns only the date rows and drops am_pm attribute", {
  out <- now(time = FALSE, tidy = TRUE)
  expect_s3_class(out, "data.frame")
  expect_true(is.numeric(out$value))
  expect_setequal(out$term, c("day", "year"))
  expect_null(attr(out, "am_pm"))
})
