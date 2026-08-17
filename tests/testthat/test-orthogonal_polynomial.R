test_that("orthonormal form matches stats::contr.poly()", {
  for (a in 3:8) {
    M <- orthogonal_polynomial(a, type = "orthonormal")
    expect_equal(matrix(as.numeric(M), nrow = a),
                 matrix(as.numeric(stats::contr.poly(a)), nrow = a),
                 tolerance = 1e-10)
  }
})

test_that("orthonormal columns sum to zero, have unit length, and are orthogonal", {
  M <- orthogonal_polynomial(6)
  expect_true(all(abs(colSums(M)) < 1e-10))            # mean zero
  expect_true(all(abs(colSums(M^2) - 1) < 1e-10))      # unit length
  cp <- crossprod(M)
  expect_true(max(abs(cp[upper.tri(cp)])) < 1e-10)     # mutually orthogonal
})

test_that("integer form reproduces the classic coefficient table", {
  M4 <- orthogonal_polynomial(4, type = "integer")
  expect_equal(M4[, "linear"],    c(L1 = -3, L2 = -1, L3 = 1, L4 = 3))
  expect_equal(M4[, "quadratic"], c(L1 =  1, L2 = -1, L3 = -1, L4 = 1))
  expect_equal(M4[, "cubic"],     c(L1 = -1, L2 =  3, L3 = -3, L4 = 1))

  M5 <- orthogonal_polynomial(5, type = "integer")
  expect_equal(unname(M5[, "linear"]),    c(-2, -1, 0, 1, 2))
  expect_equal(unname(M5[, "quadratic"]), c(2, -1, -2, -1, 2))
  expect_equal(unname(M5[, "quartic"]),   c(1, -4, 6, -4, 1))
})

test_that("integer columns are whole numbers, sum to zero, and stay orthogonal", {
  for (a in 3:10) {
    M <- orthogonal_polynomial(a, type = "integer")
    expect_true(all(abs(M - round(M)) < 1e-9))
    expect_true(all(abs(colSums(M)) < 1e-9))
    cp <- crossprod(M)
    expect_true(max(abs(cp[upper.tri(cp)])) < 1e-9)
  }
})

test_that("sum_sq attribute records the per-trend sum of squared coefficients", {
  M_on <- orthogonal_polynomial(4, type = "orthonormal")
  expect_equal(unname(attr(M_on, "sum_sq")), c(1, 1, 1), tolerance = 1e-10)

  M_int <- orthogonal_polynomial(4, type = "integer")
  expect_equal(attr(M_int, "sum_sq"),
               c(linear = 20, quadratic = 4, cubic = 20))
})

test_that("quantitative levels are used as labels and spacing", {
  M <- orthogonal_polynomial(c(1, 2, 3, 4), type = "integer")
  expect_equal(rownames(M), c("1", "2", "3", "4"))
  expect_equal(unname(M[, "linear"]), c(-3, -1, 1, 3))
})

test_that("degree truncation returns the requested number of trends", {
  M <- orthogonal_polynomial(5, degree = 2)
  expect_equal(ncol(M), 2L)
  expect_equal(colnames(M), c("linear", "quadratic"))
})

test_that("unequal spacing yields a valid orthonormal set", {
  M <- orthogonal_polynomial(c("a", "b", "c", "d"), scores = c(0, 1, 2, 4))
  expect_true(all(abs(colSums(M)) < 1e-10))
  expect_true(all(abs(colSums(M^2) - 1) < 1e-10))
  cp <- crossprod(M)
  expect_true(max(abs(cp[upper.tri(cp)])) < 1e-10)
})

test_that("integer form rejects unequal spacing", {
  expect_error(
    orthogonal_polynomial(c(0, 1, 2, 4), type = "integer"),
    "equally spaced"
  )
})

test_that("result works as a contrasts() assignment", {
  f <- factor(rep(1:4, each = 3))
  expect_silent(contrasts(f) <- orthogonal_polynomial(levels(f)))
  y <- as.numeric(f)
  fit <- lm(y ~ f)
  expect_equal(length(coef(fit)), 4L)
})

test_that("invalid arguments error informatively", {
  expect_error(orthogonal_polynomial(2), "at least 3 levels")
  expect_error(orthogonal_polynomial(4, degree = 0), "between 1 and 3")
  expect_error(orthogonal_polynomial(4, degree = 4), "between 1 and 3")
  expect_error(orthogonal_polynomial(c("a", "b", "c"), scores = c(1, 2)),
               "length 3")
})

test_that("print method shows book-style layout and returns invisibly", {
  M <- orthogonal_polynomial(4, type = "integer")
  expect_output(print(M), "Table A.10")          # book-style header
  expect_output(print(M), "sum c\\^2")            # sum c^2 as a column
  utils::capture.output(vis <- withVisible(print(M)))
  expect_identical(vis$visible, FALSE)
})

test_that("print lays trends in rows and levels in columns (transpose of storage)", {
  M  <- orthogonal_polynomial(4, type = "integer")
  out <- utils::capture.output(print(M))
  # A row label for a trend (linear) must appear; the stored matrix has the
  # trends as COLUMNS, so the printed orientation is the transpose.
  expect_true(any(grepl("^linear", out)))
  expect_true(any(grepl("quadratic", out)))
})
