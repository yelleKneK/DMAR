test_that("holzinger_swineford has the documented shape", {
  data(holzinger_swineford, package = "DMAR",
       envir = environment())
  expect_s3_class(holzinger_swineford, "data.frame")
  expect_equal(dim(holzinger_swineford), c(301L, 34L))
  expect_true(all(c(
    "id", "sex", "grade", "age", "month_since_birthday",
    "age_months", "age_years", "school",
    "t1_visual_perception", "t26_flags"
  ) %in% names(holzinger_swineford)))
})


test_that("school and grade breakdown matches the 1939 monograph", {
  data(holzinger_swineford, package = "DMAR",
       envir = environment())
  schools <- table(holzinger_swineford$school)
  expect_equal(as.integer(schools["Pasteur"]),     156L)
  expect_equal(as.integer(schools["Grant-White"]), 145L)
  grades <- table(holzinger_swineford$grade)
  expect_equal(as.integer(grades["7"]), 157L)
  expect_equal(as.integer(grades["8"]), 144L)
})


test_that("t25 and t26 are Grant-White only", {
  data(holzinger_swineford, package = "DMAR",
       envir = environment())
  d <- holzinger_swineford
  pas <- d$school == "Pasteur"
  gw  <- d$school == "Grant-White"
  expect_true(all(is.na(d$t25_paper_form_board_r[pas])))
  expect_true(all(is.na(d$t26_flags[pas])))
  expect_true(all(!is.na(d$t25_paper_form_board_r[gw])))
  expect_true(all(!is.na(d$t26_flags[gw])))
})


test_that("age_months is internally consistent with age and month_since_birthday", {
  data(holzinger_swineford, package = "DMAR",
       envir = environment())
  d <- holzinger_swineford
  expect_equal(d$age_months, 12L * d$age + d$month_since_birthday)
})


test_that("holzinger_swineford carries the post-4.6.0 MBESS corrections", {
  # Spot-check several cells that distinguish the corrected data from
  # the older snapshot still resident in the sem package's HS.data and
  # OpenMx HS.ability.data.
  data(holzinger_swineford, package = "DMAR",
       envir = environment())
  d <- holzinger_swineford
  expect_equal(d$t20_deduction[d$id == 2L],         -3L)
  expect_equal(d$t22_problem_reasoning[d$id == 2L], 21L)
  expect_equal(d$t24_woody_mccall[d$id == 2L],      12L)
  expect_equal(d$t15_number_recognition[d$id == 1L], 86L)
})


test_that("no alias binding shadows the canonical dataset name", {
  # The HS_Data active-binding alias was removed 2026-08-15; the
  # canonical snake_case name is the only route to the data.
  expect_false("HS_Data" %in% getNamespaceExports("DMAR"))
})
