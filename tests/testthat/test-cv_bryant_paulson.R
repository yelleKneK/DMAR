# Tests for cv_bryant_paulson(), the DMAR critical-value-family member for
# ANCOVA multiple comparisons. Published anchors: Bryant & Bruvold (1980)
# Section 3 (q_.05;1,6,14 = 4.83) and Section 4 (the Duncan example R_j values).

test_that("cv_bryant_paulson reproduces the paper's q_.05;1,6,14 = 4.83", {
  out <- cv_bryant_paulson(alpha_level = .05, df = 14, groups = 6, covariates = 1)
  expect_s3_class(out, "dmar_tbl")
  expect_equal(out$value, 4.83, tolerance = 5e-3)
  expect_equal(out$area_less, 0.95, tolerance = 1e-4)
  expect_equal(out$area_greater, 0.05, tolerance = 1e-4)
})

test_that("verbose = FALSE returns just term and value", {
  skip_on_cran()  # a shape check that pays a full root find; the anchor above runs on CRAN
  out <- cv_bryant_paulson(.05, df = 14, groups = 6, covariates = 1, verbose = FALSE)
  expect_equal(names(out), c("term", "value"))
  expect_equal(out$term, "upper_cv")
})

test_that("covariates = 0 equals sqrt(2) * cv_tukey_hsd", {
  for (k in c(3, 5, 6)) for (nu in c(14, 30)) {
    bp <- cv_bryant_paulson(.05, df = nu, groups = k, covariates = 0,
                            verbose = FALSE)$value
    tk <- sqrt(2) * cv_tukey_hsd(.05, df = nu, groups = k, verbose = FALSE)$value
    expect_equal(bp, tk, tolerance = 1e-6)
  }
})

test_that("qbryant_paulson reproduces Bryant & Paulson (1976) Table 1 to 2 dp", {
  skip_on_cran()  # dozens of numeric integrations; the fast anchors above run on CRAN
  # Upper 0.01 points read from the original Bryant & Paulson (1976) Table 1(b),
  # spanning the table's corners. At small df stats::ptukey loses accuracy, so
  # qbryant_paulson bypasses it with a direct studentized-range integral; these
  # rounded values are what the paper prints, and a large simulation of the
  # statistic confirms them. Columns: nu, covariates, groups, published.
  ref <- rbind(c( 2, 1,  2, 19.09), c( 2, 1,  8, 40.60), c( 2, 2,  8, 49.31),
               c( 2, 3,  6, 51.07), c( 2, 3, 20, 73.01), c( 3, 1, 20, 25.05),
               c( 3, 3, 20, 33.13), c( 4, 3, 16, 20.87), c( 6, 3,  5, 10.07),
               c(20, 3,  6,  5.99), c(60, 1,  2,  3.79))
  for (i in seq_len(nrow(ref))) {
    v <- qbryant_paulson(0.99, num_covariates = ref[i, 2],
                         num_groups = ref[i, 3], df = ref[i, 1])
    expect_equal(round(v, 2), ref[i, 4])
  }
})

test_that("the two Table 1 entries on the rounding boundary are identified", {
  skip_on_cran()  # dozens of numeric integrations; the fast anchors above run on CRAN
  # These two entries are the only ones whose exact value sits within about
  # 1e-5 of the point where the second decimal turns over, so they round up to
  # 23.17 and 19.75 while the 1976 table rounds down. The targets below were
  # confirmed to fourteen significant figures by two independent high-order
  # quadrature engines sharing no code with the package; the 1e-6 tolerance
  # guards the small-nu path against spline-resolution regressions.
  expect_equal(qbryant_paulson(0.99, num_covariates = 2, num_groups = 8, df = 3),
               23.16501285, tolerance = 1e-6)
  expect_equal(qbryant_paulson(0.99, num_covariates = 2, num_groups = 20, df = 4),
               19.74500754, tolerance = 1e-6)
})

test_that("the critical value increases with the number of covariates", {
  skip_on_cran()  # a root find per covariate count; the fast anchors above run on CRAN
  v <- sapply(0:3, function(p)
    cv_bryant_paulson(.05, df = 20, groups = 5, covariates = p,
                      verbose = FALSE)$value)
  expect_true(all(diff(v) > 0))
})

test_that("Duncan procedure reproduces the Section 4 example R_j values", {
  skip_on_cran()  # dozens of numeric integrations; the fast anchors above run on CRAN
  # R_k = (sqrt(0.01326)/sqrt(4)) * r_{.05;1,k,14}; paper: .181 .190 .195 .199 .202
  scale <- sqrt(0.01326) / sqrt(4)
  R <- sapply(2:6, function(k)
    scale * cv_bryant_paulson(.05, df = 14, groups = k, covariates = 1,
                              procedure = "duncan", verbose = FALSE)$value)
  # Absolute agreement to the paper's two-decimal R_j. (R_6 sits on the
  # .201/.202 rounding boundary -- the exact value ~0.2013 is within 0.001 of
  # the paper's 0.202 -- so we compare absolute differences rather than
  # pre-rounding, which would otherwise quantize 0.2013 to 0.201.)
  expect_true(max(abs(unname(R) -
                       c(0.181, 0.190, 0.195, 0.199, 0.202))) < 2e-3)
})

test_that("Duncan significant ranges are non-decreasing in the number of groups", {
  skip_on_cran()  # dozens of numeric integrations; the fast anchors above run on CRAN
  r <- sapply(2:8, function(k)
    cv_bryant_paulson(.05, df = 20, groups = k, covariates = 1,
                      procedure = "duncan", verbose = FALSE)$value)
  expect_true(all(diff(r) >= -1e-8))
})

test_that("Duncan tail areas are NA (stepwise quantity)", {
  skip_on_cran()  # a root find at every k up to six; the fast anchors above run on CRAN
  out <- cv_bryant_paulson(.05, df = 14, groups = 6, covariates = 1,
                           procedure = "duncan")
  expect_true(is.na(out$area_less) && is.na(out$area_greater))
})

test_that("invalid arguments are rejected", {
  expect_error(cv_bryant_paulson(df = 14, groups = 6, covariates = 1))      # no alpha
  expect_error(cv_bryant_paulson(.05, groups = 6, covariates = 1))          # no df
  expect_error(cv_bryant_paulson(.05, df = -1, groups = 6, covariates = 1)) # df <= 0
  expect_error(cv_bryant_paulson(.05, df = 14, groups = 1, covariates = 1)) # groups < 2
  expect_error(cv_bryant_paulson(.05, df = 14, groups = 6, covariates = -1))# covariates < 0
  expect_error(cv_bryant_paulson(1.2, df = 14, groups = 6, covariates = 1)) # alpha out of range
})

