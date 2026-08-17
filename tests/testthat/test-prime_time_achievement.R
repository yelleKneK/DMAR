test_that("prime_time_achievement has the documented shape", {
  data(prime_time_achievement, package = "DMAR", envir = environment())
  expect_s3_class(prime_time_achievement, "data.frame")
  expect_equal(dim(prime_time_achievement), c(10927L, 113L))
  expect_true(all(c("id", "region", "corp", "school", "class",
                    "corp_id", "school_id", "class_id",
                    "nctotal", "ncread", "ncmath", "nclang",
                    "ptia", "ptratio", "classize", "clenroll",
                    "ses", "geog", "race", "gender") %in%
                  names(prime_time_achievement)))
})

test_that("derived cluster IDs reconcile with the published counts", {
  data(prime_time_achievement, package = "DMAR", envir = environment())
  # Lapsley et al. (2002, ERIC ED466679):
  expect_equal(length(unique(prime_time_achievement$corp_id)),    61L)
  expect_equal(length(unique(prime_time_achievement$school_id)), 163L)
  # 586 by paste(corp, school, class); paper reports 573 by a
  # different counting rule (see Details).
  expect_equal(length(unique(prime_time_achievement$class_id)),  586L)
  expect_equal(length(unique(prime_time_achievement$region)),      9L)
})

test_that("corporation cluster sizes match the documented range and median", {
  # Anchors the ?prime_time_achievement Details prose. The summaries are
  # computed from the derived corp_id key (61 corporations); the bare corp
  # column merges the two corporations that share code 2400 and gave the
  # 15 to 808 (median 124) figures the page previously reported.
  data(prime_time_achievement, package = "DMAR", envir = environment())
  sizes <- as.numeric(table(prime_time_achievement$corp_id))
  expect_length(sizes, 61L)
  expect_equal(min(sizes), 15)
  expect_equal(max(sizes), 756)
  expect_equal(median(sizes), 117)
})

test_that("nesting structure is consistent", {
  data(prime_time_achievement, package = "DMAR", envir = environment())
  # school_id is unique within corp_id
  s <- unique(prime_time_achievement[, c("corp_id", "school_id")])
  expect_equal(length(unique(s$school_id)), nrow(s))
  # class_id is unique within school_id
  c <- unique(prime_time_achievement[, c("school_id", "class_id")])
  expect_equal(length(unique(c$class_id)), nrow(c))
})

test_that("demographic distributions match Lapsley et al. (2002)", {
  data(prime_time_achievement, package = "DMAR", envir = environment())
  expect_equal(sum(prime_time_achievement$gender == 1L, na.rm = TRUE), 5425L)
  expect_equal(sum(prime_time_achievement$gender == 2L, na.rm = TRUE), 5457L)
  expect_equal(sum(prime_time_achievement$race == 5L, na.rm = TRUE), 9207L)  # Caucasian
  expect_equal(sum(prime_time_achievement$race == 2L, na.rm = TRUE), 995L)   # African Am.
  expect_equal(sum(prime_time_achievement$race == 4L, na.rm = TRUE), 348L)   # Hispanic
  expect_equal(sum(prime_time_achievement$race == 1L, na.rm = TRUE), 16L)    # Am. Indian
  expect_equal(sum(prime_time_achievement$race == 3L, na.rm = TRUE), 62L)    # Asian
  expect_equal(sum(prime_time_achievement$race == 6L, na.rm = TRUE), 188L)   # Multi-racial
})

test_that("treatment indicators match the manuscript subgroup sizes", {
  data(prime_time_achievement, package = "DMAR", envir = environment())
  expect_equal(sum(prime_time_achievement$ptia == 1L), 4021L)  # aide
  expect_equal(sum(prime_time_achievement$ptia == 2L), 6789L)  # no aide
  expect_equal(sum(prime_time_achievement$ptia == 3L), 117L)   # other
  expect_setequal(unique(prime_time_achievement$classize), 1:4)
})

test_that("missing-data codes have been cleaned", {
  data(prime_time_achievement, package = "DMAR", envir = environment())
  for (v in c("age", "gender", "race",
              "geread", "gevocab", "gemath", "getotal",
              "ncread", "nclang", "ncmath", "nctotal",
              "aatotal", "csi", "npatotal")) {
    expect_false(any(prime_time_achievement[[v]] == 999, na.rm = TRUE),
                 info = paste(v, "still contains 999"))
  }
  expect_false(any(prime_time_achievement$typmulti == 888, na.rm = TRUE))
  expect_false(any(prime_time_achievement$ptstatus == 888, na.rm = TRUE))
})

test_that("SPSS variable labels are preserved", {
  data(prime_time_achievement, package = "DMAR", envir = environment())
  expect_equal(attr(prime_time_achievement$nctotal, "label"), "NCE TOTAL")
  expect_equal(attr(prime_time_achievement$ptia,    "label"),
               "PRESENCE OF A PRIME TIME IA?")
  expect_equal(attr(prime_time_achievement$ses,     "label"),
               "SOCIOECONOMIC STATUS")
  # Derived IDs also carry labels.
  expect_match(attr(prime_time_achievement$corp_id,   "label"),
               "corporation identifier")
  expect_match(attr(prime_time_achievement$school_id, "label"),
               "school identifier")
  expect_match(attr(prime_time_achievement$class_id,  "label"),
               "classroom identifier")
})

test_that("FILTER_$ has been dropped from the SPSS source", {
  data(prime_time_achievement, package = "DMAR", envir = environment())
  expect_false("filter_$"  %in% names(prime_time_achievement))
  expect_false("filter_dol" %in% names(prime_time_achievement))
})

test_that("corp 2400 spans regions 2 and 3 (data integrity note)", {
  data(prime_time_achievement, package = "DMAR", envir = environment())
  rs <- with(subset(prime_time_achievement, corp == 2400L),
             sort(unique(region)))
  expect_equal(rs, c(2L, 3L))
})

test_that("three-level null model recovers documented ICCs", {
  # Slow (a three-level fit on all 10,927 students), so it skips on CRAN;
  # the corporation cluster-size anchor above is the fast check that stays.
  # The expected values anchor the Details prose, computed from corp_id;
  # the looser 0.059 the page previously reported came from grouping on
  # the bare corp column.
  skip_on_cran()
  skip_if_not_installed("lme4")
  data(prime_time_achievement, package = "DMAR", envir = environment())
  m <- lme4::lmer(nctotal ~ 1 + (1 | corp_id/school_id),
                  data = prime_time_achievement,
                  control = lme4::lmerControl(check.conv.singular = "ignore"))
  vc  <- as.data.frame(lme4::VarCorr(m))
  vC  <- vc$vcov[vc$grp == "corp_id"]
  vS  <- vc$vcov[vc$grp == "school_id:corp_id"]
  vE  <- vc$vcov[vc$grp == "Residual"]
  tot <- vC + vS + vE
  expect_equal(vC, 16.294, tolerance = 1e-3)
  expect_equal(vS, 22.716, tolerance = 1e-3)
  expect_equal(vE, 240.434, tolerance = 1e-3)
  expect_equal(vC / tot, 0.0583, tolerance = 2e-3)
  expect_equal(vS / tot, 0.0813, tolerance = 2e-3)
})

test_that("the canonical name is the only route to the data", {
  # The Prime_Time active-binding alias was removed 2026-08-15; the
  # canonical snake_case name stands alone, and its Rd page carries
  # no alias entries beyond the canonical name.
  expect_false("Prime_Time" %in% getNamespaceExports("DMAR"))
  rd_path <- testthat::test_path("..", "..", "man",
                                 "prime_time_achievement.Rd")
  skip_if_not(file.exists(rd_path), "Rd source not on disk")
  rd_lines <- readLines(rd_path, warn = FALSE)
  expect_true(any(grepl("^\\\\alias\\{prime_time_achievement\\}$",
                        rd_lines)))
  expect_false(any(grepl("^\\\\alias\\{Prime_Time\\}$", rd_lines)))
  expect_false(any(grepl("^\\\\alias\\{prime_time\\}$", rd_lines)))
})
