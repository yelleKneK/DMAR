test_that("simple_structure() matches the Hofmann complexity by hand", {
  Lambda <- rbind(
    i1 = c(0.80, 0.05), i2 = c(0.75, 0.10), i3 = c(0.70, -0.05),
    i4 = c(0.08, 0.78), i5 = c(-0.04, 0.72), i6 = c(0.30, 0.60))
  res <- simple_structure(Lambda)
  expect_s3_class(res, "dmar_tbl")

  cx_hand <- rowSums(Lambda^2)^2 / rowSums(Lambda^4)
  expect_equal(unname(attr(res, "complexity")), unname(cx_hand))
  expect_equal(res$value[res$term == "mean_complexity"], mean(cx_hand))

  # Pure / complex at salient = .30: i6 loads .30 and .60 -> complex.
  expect_equal(res$value[res$term == "n_complex"], 1)
  expect_equal(res$value[res$term == "n_pure"], 5)
  expect_equal(res$value[res$term == "items"], 6)
  expect_equal(res$value[res$term == "factors"], 2)

  # Hyperplane proportion: |loading| < .10.
  expect_equal(res$value[res$term == "hyperplane_proportion"],
               mean(abs(Lambda) < 0.10))
})

test_that("simple_structure() agrees with psych::fa complexity", {
  # Pinned from psych::fa (psych 2.5.3, 2026-08-09): the oblimin loadings and
  # per-item Hofmann complexities for a two-factor solution on
  # psych::bfi[, 1:10] with set.seed(7). Live comparison in
  # tools/oracle_checks.R.
  L <- matrix(
    c( 0.0795515078441172,  0.0070031916378957, -0.0282778823739929,
       0.1449343431711100,  0.0274581076243546,  0.5707308291482110,
       0.6368103209561290,  0.5415625900515870, -0.6492041315662060,
      -0.5617797510344260,
      -0.4054913957089600,  0.6773122523760470,  0.7595180599883370,
       0.4387135304696980,  0.6023711085906150, -0.0607046910077694,
      -0.0132766740435319,  0.0316143727528687, -0.00366023075383061,
      -0.0579471993068160),
    nrow = 10,
    dimnames = list(c("A1", "A2", "A3", "A4", "A5",
                      "C1", "C2", "C3", "C4", "C5"),
                    c("MR1", "MR2")))
  complexity_psych <- c(1.07686358931520, 1.00021381794352, 1.00277234145027,
                        1.21570890925556, 1.00415567005727, 1.02262331571498,
                        1.00086933729720, 1.00681548518651, 1.00006357470836,
                        1.02127715575166)
  res <- simple_structure(L)
  expect_equal(unname(attr(res, "complexity")), complexity_psych,
               tolerance = 1e-8)
})

test_that("simple_structure() matches Hofmann (1978) computed from the paper", {
  # Hofmann (1978, p. 248) prints no worked numbers, so the anchors are the
  # printed equations. Equation 1 defines the complexity of item i as
  # c_i = (sum_j a_ij^2)^2 / sum_j a_ij^4; it is implemented here element by
  # element, straight from the paper, as an independent check on the
  # vectorized computation.
  set.seed(113)
  Lambda <- matrix(stats::runif(24, -0.2, 0.9), nrow = 8, ncol = 3)
  res <- simple_structure(Lambda)
  c_fn <- unname(attr(res, "complexity"))

  c_paper <- numeric(nrow(Lambda))
  for (i in seq_len(nrow(Lambda))) {
    sum_sq <- 0
    sum_quart <- 0
    for (j in seq_len(ncol(Lambda))) {
      sum_sq <- sum_sq + Lambda[i, j]^2
      sum_quart <- sum_quart + Lambda[i, j]^4
    }
    c_paper[i] <- sum_sq^2 / sum_quart
  }
  expect_equal(c_fn, c_paper)

  # The complexity ranges from one (a perfect unifactor item) to the number
  # of factors (Hofmann, 1978, p. 247).
  expect_true(all(c_fn >= 1 & c_fn <= ncol(Lambda)))

  # The mean complexity row is Hofmann's total matrix complexity, the
  # arithmetic average of the c_i (Hofmann, 1978, p. 248).
  expect_equal(res$value[res$term == "mean_complexity"], mean(c_paper))

  # Equation 3 (Hofmann, 1978, p. 248) expresses Kaiser's item simplicity
  # through the complexity: s_i = [1/(m - 1)][(m / c_i) - 1]. Pushing the
  # returned complexities through Equation 3 must recover Kaiser's (1974)
  # simplicity computed directly from its own definition (Hofmann's
  # Equation 2).
  m <- ncol(Lambda)
  s_kaiser <- (m * rowSums(Lambda^4) / rowSums(Lambda^2)^2 - 1) / (m - 1)
  expect_equal((1 / (m - 1)) * ((m / c_fn) - 1), s_kaiser)
})

test_that("simple_structure() validates input", {
  expect_error(simple_structure(1:5), "loading matrix")
  expect_error(simple_structure(matrix(NA_real_, 3, 2)), "loading matrix")
  expect_error(simple_structure(matrix(0.5, 3, 2), salient = 2), "\\[0, 1\\]")
})
