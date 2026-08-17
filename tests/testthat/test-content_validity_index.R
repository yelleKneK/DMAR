test_that("content_validity_index() reproduces Polit, Beck, and Owen (2007)", {
  # Their worked example: 6 experts, 5 of whom rate the item relevant.
  # I-CVI = 0.83, p_c = 0.094, modified kappa = 0.81.
  res <- content_validity_index(matrix(c(4, 4, 3, 4, 3, 1), nrow = 1))
  expect_s3_class(res, "dmar_tbl")
  expect_equal(nrow(res), 1L)
  expect_equal(res$n_experts, 6)
  expect_equal(res$n_relevant, 5)
  expect_equal(round(res$i_cvi, 2), 0.83)
  expect_equal(res$i_cvi, 5 / 6, tolerance = 1e-12)

  # p_c re-derived from the published formula with factorial() at this small N.
  p_c <- (factorial(6) / (factorial(5) * factorial(1))) * 0.5^6
  expect_equal(round(p_c, 3), 0.094)
  expect_equal(res$kappa, (5 / 6 - p_c) / (1 - p_c), tolerance = 1e-12)

  # The published kappa of 0.81 is what the paper's own display precision
  # gives: carrying their rounded I-CVI (0.83) and p_c (0.094) through the
  # formula returns 0.81, while the unrounded computation returns 0.8161.
  expect_equal(round((0.83 - 0.094) / (1 - 0.094), 2), 0.81)
  expect_equal(round(res$kappa, 4), 0.8161)
})

test_that("modified kappa stays stable for a large panel", {
  # 200 experts, 150 relevant: factorial(200) overflows to Inf, lchoose does not.
  ratings <- matrix(c(rep(4, 150), rep(1, 50)), nrow = 1)
  res <- content_validity_index(ratings)
  p_c <- exp(lchoose(200, 150) + 200 * log(0.5))
  expect_true(is.finite(res$kappa))
  expect_equal(res$kappa, (0.75 - p_c) / (1 - p_c), tolerance = 1e-12)
  expect_equal(res$i_cvi, 0.75, tolerance = 1e-12)
})

test_that("Lawshe's CVR matches hand computation at the four anchors", {
  # N = 10 experts; n_e essential varies by row.
  mk <- function(n_e) c(rep(4, n_e), rep(1, 10 - n_e))
  ratings <- rbind(mk(9), mk(5), mk(10), mk(0))
  res <- content_validity_index(ratings)
  expect_equal(res$cvr, c(0.8, 0, 1, -1), tolerance = 1e-12)

  # A separate essential dichotomization changes cvr but not i_cvi.
  res2 <- content_validity_index(ratings, essential = 4)
  expect_equal(res2$cvr, res$cvr, tolerance = 1e-12)
  res3 <- content_validity_index(ratings, essential = c(2, 3, 4))
  expect_equal(res3$cvr, res$cvr, tolerance = 1e-12)   # no 2s or 3s present
  expect_equal(res3$i_cvi, res$i_cvi, tolerance = 1e-12)

  # essential = 1 counts the "not relevant" ratings, so the CVR flips sign.
  res4 <- content_validity_index(ratings, essential = 1)
  expect_equal(res4$cvr, -res$cvr, tolerance = 1e-12)
  expect_identical(attr(res4, "essential"), 1)
  expect_identical(attr(res, "essential"), c(3, 4))
})

test_that("the interval is exactly stats::binom.test's exact interval", {
  ratings <- rbind(
    a = c(4, 4, 3, 4, 4, 3),
    b = c(4, 3, 4, 4, 3, 2),
    c = c(2, 3, 1, 2, 3, 2),
    d = c(4, 4, 4, 4, 4, 4))
  for (cl in c(0.90, 0.95, 0.99)) {
    res <- content_validity_index(ratings, conf_level = cl)
    for (i in seq_len(nrow(res))) {
      ci <- stats::binom.test(res$n_relevant[i], res$n_experts[i],
                              conf.level = cl)$conf.int
      expect_equal(res$ci_lower[i], ci[1], tolerance = 0)
      expect_equal(res$ci_upper[i], ci[2], tolerance = 0)
    }
    expect_equal(attr(res, "conf_level"), cl)
  }
  # A perfect item has an upper limit of exactly 1 and a lower limit below 1.
  res <- content_validity_index(ratings)
  expect_equal(res$ci_upper[4], 1)
  expect_true(res$ci_lower[4] < 1)
  expect_true(all(res$ci_lower <= res$i_cvi & res$i_cvi <= res$ci_upper))
})

test_that("scale level summaries are the mean and the universal agreement rate", {
  ratings <- rbind(
    a = c(4, 4, 3, 4, 4, 3),   # 6/6 = 1
    b = c(4, 3, 4, 4, 3, 2),   # 5/6
    c = c(2, 3, 1, 2, 3, 2),   # 2/6
    d = c(4, 4, 4, 4, 4, 4))   # 6/6 = 1
  res <- content_validity_index(ratings)
  expect_equal(res$i_cvi, c(1, 5 / 6, 2 / 6, 1), tolerance = 1e-12)
  expect_equal(attr(res, "s_cvi_ave"), mean(res$i_cvi), tolerance = 1e-12)
  expect_equal(attr(res, "s_cvi_ave"), mean(c(1, 5 / 6, 2 / 6, 1)),
               tolerance = 1e-12)
  expect_equal(attr(res, "s_cvi_ua"), mean(res$i_cvi == 1), tolerance = 1e-12)
  expect_equal(attr(res, "s_cvi_ua"), 0.5, tolerance = 1e-12)
  expect_identical(res$item, c("a", "b", "c", "d"))
  expect_identical(attr(res, "relevant"), c(3, 4))
})

test_that("a missing rating reduces only its own item's denominator", {
  full <- rbind(
    a = c(4, 4, 3, 4, 4, 3),
    b = c(4, 3, 4, 4, 3, 2))
  gapped <- full
  gapped[2, 6] <- NA                       # expert 6 did not rate item b
  res_full <- content_validity_index(full)
  res_gap  <- content_validity_index(gapped)

  expect_equal(res_gap$n_experts, c(6, 5))
  expect_equal(res_gap$n_relevant, c(6, 5))
  expect_equal(res_gap$i_cvi, c(1, 1), tolerance = 1e-12)   # 5/5, not 5/6
  expect_equal(res_full$i_cvi[2], 5 / 6, tolerance = 1e-12)
  # Item a is untouched by the gap in item b's row.
  expect_equal(res_gap$i_cvi[1], res_full$i_cvi[1], tolerance = 1e-12)
  expect_equal(res_gap$kappa[1], res_full$kappa[1], tolerance = 1e-12)
  expect_equal(res_gap$cvr[1], res_full$cvr[1], tolerance = 1e-12)
  # The reduced denominator propagates to the CVR and the interval too.
  expect_equal(res_gap$cvr[2], (5 - 2.5) / 2.5, tolerance = 1e-12)
  expect_equal(res_gap$ci_lower[2], stats::binom.test(5, 5)$conf.int[1],
               tolerance = 0)
})

test_that("value columns are numeric and nothing is pre-rounded", {
  res <- content_validity_index(rbind(a = c(4, 4, 3, 1, 1, 1)))
  num <- c("n_experts", "n_relevant", "i_cvi", "ci_lower", "ci_upper",
           "kappa", "cvr")
  expect_true(all(vapply(res[num], is.numeric, logical(1))))
  expect_type(res$item, "character")
  expect_equal(res$i_cvi, 0.5, tolerance = 1e-15)
  # Full precision is stored, not the 3 significant figures that print.
  expect_true(abs(res$ci_lower - 0.1181163) < 1e-6)
})

test_that("content_validity_index() fails clearly on bad input", {
  ok <- rbind(a = c(4, 4, 3, 2))
  expect_error(content_validity_index(matrix(letters[1:4], nrow = 1)),
               "numeric")
  expect_error(content_validity_index(matrix(c(4, 3), nrow = 2)),
               "at least 2 experts")
  expect_error(content_validity_index(matrix(c(4, 3, 5, 2), nrow = 1)),
               "4 point relevance scale")
  expect_error(content_validity_index(matrix(c(4, 3, 0, 2), nrow = 1)),
               "4 point relevance scale")
  expect_error(content_validity_index(ok, relevant = 5),
               "'relevant' must be")
  expect_error(content_validity_index(ok, essential = 0),
               "'essential' must be")
  expect_error(content_validity_index(ok, conf_level = 1),
               "conf_level")
  expect_error(content_validity_index(ok, conf_level = 0),
               "conf_level")
  expect_error(content_validity_index(rbind(a = c(4, 3, 2, 1),
                                            b = c(NA, NA, NA, NA))),
               "entirely missing")
  expect_error(content_validity_index(list(1, 2)), "matrix or data.frame")
  expect_error(content_validity_index(matrix(numeric(0), nrow = 0, ncol = 3)),
               "at least 1 item")
  # No silent NA: a legitimate call returns finite numbers throughout.
  res <- content_validity_index(ok)
  expect_false(anyNA(res$i_cvi))
  expect_false(anyNA(res$kappa))
  expect_false(anyNA(res$cvr))
})

test_that("a data.frame input is accepted and row names name the items", {
  d <- data.frame(e1 = c(4, 2), e2 = c(3, 1), e3 = c(4, 3),
                  row.names = c("clarity", "scope"))
  res <- content_validity_index(d)
  expect_identical(res$item, c("clarity", "scope"))
  expect_equal(res$i_cvi, c(1, 1 / 3), tolerance = 1e-12)
  # Without row names the items are labeled positionally.
  d2 <- data.frame(e1 = c(4, 2), e2 = c(3, 1), e3 = c(4, 3))
  expect_identical(content_validity_index(d2)$item, c("item_1", "item_2"))
})
