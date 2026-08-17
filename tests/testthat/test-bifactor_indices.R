.bifactor_fit <- function() {
  skip_if_not_installed("lavaan")
  # Simulate clean orthogonal bifactor data so the solution is proper.
  set.seed(11)
  n <- 600
  g <- rnorm(n); f <- list(rnorm(n), rnorm(n), rnorm(n))
  x <- vapply(1:9, function(i) {
    grp <- ceiling(i / 3)
    0.5 * g + 0.5 * f[[grp]] + sqrt(1 - 0.25 - 0.25) * rnorm(n)
  }, numeric(n))
  colnames(x) <- paste0("x", 1:9)
  d <- as.data.frame(x)
  model <- paste(
    "g  =~", paste(colnames(x), collapse = " + "), "\n",
    "f1 =~ x1 + x2 + x3\n",
    "f2 =~ x4 + x5 + x6\n",
    "f3 =~ x7 + x8 + x9")
  lavaan::cfa(model, data = d, orthogonal = TRUE, std.lv = TRUE)
}

test_that("bifactor_indices() reproduces ECV and omega_H by hand", {
  fit <- .bifactor_fit()
  res <- bifactor_indices(fit)
  expect_s3_class(res, "dmar_tbl")
  expect_equal(nrow(res), 4L)                       # general + 3 group factors
  expect_equal(res$factor[1], "g")
  expect_false(isTRUE(attr(res, "improper")))

  # Hand computation from the standardized loadings.
  std <- lavaan::standardizedSolution(fit)
  lam <- std[std$op == "=~", ]
  items <- unique(lam$rhs)
  lg <- setNames(lam$est.std[lam$lhs == "g"][match(items, lam$rhs[lam$lhs == "g"])], items)
  groups <- c("f1", "f2", "f3"); ls <- setNames(rep(NA_real_, length(items)), items)
  grp <- ls
  for (g in groups) {
    gi <- lam$rhs[lam$lhs == g]
    ls[gi] <- lam$est.std[lam$lhs == g][match(gi, lam$rhs[lam$lhs == g])]
    grp[gi] <- g
  }
  theta <- 1 - lg^2 - ls^2
  rv <- std[std$op == "~~" & std$lhs == std$rhs & std$lhs %in% items, ]
  theta[rv$lhs] <- rv$est.std

  ecv <- sum(lg^2) / (sum(lg^2) + sum(ls^2))
  total_common <- sum(lg)^2 + sum(vapply(groups, function(g) sum(ls[grp == g])^2, numeric(1)))
  total_var <- total_common + sum(theta)
  omega_h <- sum(lg)^2 / total_var

  expect_equal(res$ECV[1], ecv, tolerance = 1e-8)
  expect_equal(res$omega_H[1], omega_h, tolerance = 1e-8)

  # PUC: 9 items, 3 groups of 3 -> within pairs 3*3 = 9, total C(9,2) = 36.
  expect_equal(res$PUC[1], 1 - 9 / 36, tolerance = 1e-10)

  # Group rows carry omega_HS, not omega_H/PUC.
  expect_true(all(is.na(res$omega_H[-1])))
  expect_true(all(is.na(res$PUC[-1])))
  expect_true(all(is.finite(res$omega_HS[-1])))
})

test_that("bifactor_indices() validates and detects the general factor", {
  skip_if_not_installed("lavaan")
  fit1 <- lavaan::cfa(
    "a =~ t1_visual_perception + t2_cubes + t4_lozenges
     b =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    data = holzinger_swineford, std.lv = TRUE)
  expect_error(bifactor_indices(fit1), "general factor")  # no factor loads all
  expect_error(bifactor_indices(lm(1 ~ 1)), "lavaan fit")
})
