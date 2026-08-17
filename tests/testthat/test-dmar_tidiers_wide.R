# tidy() and glance() on the wide tables: the four bespoke pairs
# (content_validity_index, dmacs, measurement_invariance,
# measurement_alignment) and the generic wide fallback (htmt,
# average_variance_extracted). All columns are DMAR-native names;
# no dots.

test_that("content_validity_index() answers tidy() and glance()", {
  set.seed(113)
  r <- matrix(sample(2:4, 36, replace = TRUE), 12, 3)
  colnames(r) <- paste0("e", 1:3)
  res <- content_validity_index(r)

  td <- generics::tidy(res)
  expect_identical(names(td)[1:2], c("term", "estimate"))
  expect_identical(td$term, as.character(res$item))
  expect_identical(td$estimate, res$i_cvi)
  expect_true(all(c("ci_lower", "ci_upper", "cvr") %in% names(td)))
  expect_false(any(grepl(".", names(td), fixed = TRUE)))

  gl <- generics::glance(res)
  expect_identical(nrow(gl), 1L)
  expect_identical(gl$s_cvi_ave, attr(res, "s_cvi_ave"))
  expect_identical(gl$s_cvi_ua, attr(res, "s_cvi_ua"))
  expect_identical(gl$n_items, nrow(res))
})

test_that("dmacs() answers tidy() and glance()", {
  res <- dmacs(lambda_reference = c(.8, .7, .6),
               lambda_focal = c(.75, .5, .62),
               nu_reference = c(0, 0, 0), nu_focal = c(.2, 0, -.1),
               mean_focal = .3, sd_focal = 1.2, sd_pooled = c(1, 1, 1))
  td <- generics::tidy(res)
  expect_identical(td$term, as.character(res$item))
  expect_identical(td$estimate, res$dmacs)
  expect_true("lambda_focal" %in% names(td))
  gl <- generics::glance(res)
  expect_identical(gl$n_items, 3L)
  expect_identical(gl$mean_focal, 0.3)
  expect_identical(gl$sd_focal, 1.2)
})

test_that("measurement_invariance() answers tidy() and glance()", {
  skip_if_not_installed("lavaan")
  data(holzinger_swineford)
  items <- c("t6_paragraph_comprehension", "t7_sentence",
             "t9_word_meaning")
  res <- measurement_invariance(holzinger_swineford, items = items,
                                group = "grade")
  td <- generics::tidy(res)
  expect_identical(names(td)[1], "term")
  expect_identical(td$term, as.character(res$level))
  expect_identical(td$chi_square, res$chi_square)
  expect_identical(td$delta_cfi, res$delta_cfi)
  expect_null(attr(td, "fits"))
  gl <- generics::glance(res)
  expect_identical(gl$n_levels, nrow(res))
  expect_identical(gl$estimator, "ML")
  expect_identical(gl$fit_indices, "standard")
})

test_that("measurement_alignment() answers tidy() and glance()", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  sim <- do.call(rbind, lapply(1:3, function(g) {
    eta <- rnorm(120, c(0, .3, -.2)[g], c(1, 1.2, .9)[g])
    d <- as.data.frame(sapply(1:4, function(i)
      1 + .8 * eta + rnorm(120, 0, .6)))
    names(d) <- paste0("x", 1:4)
    d$cohort <- paste0("c", g)
    d
  }))
  res <- measurement_alignment(sim, items = paste0("x", 1:4),
                               group = "cohort", seed = 113)
  td <- generics::tidy(res)
  expect_identical(td$term, as.character(res$group))
  expect_identical(td$factor_mean, res$factor_mean)
  expect_identical(td$factor_variance, res$factor_variance)
  gl <- generics::glance(res)
  expect_identical(gl$n_groups, nrow(res))
  expect_equal(gl$R2_loadings_mean,
               mean(unlist(attr(res, "R2_loadings"))))
  expect_true(gl$converged)
})

test_that("the generic wide fallback serves htmt() and AVE", {
  skip_if_not_installed("lavaan")
  set.seed(113)
  f1 <- rnorm(80); f2 <- 0.4 * f1 + rnorm(80, 0, sqrt(1 - 0.16))
  d <- data.frame(a1 = f1 + rnorm(80, 0, .5), a2 = f1 + rnorm(80, 0, .5),
                  b1 = f2 + rnorm(80, 0, .5), b2 = f2 + rnorm(80, 0, .5))
  h <- htmt(d, blocks = list(A = c("a1", "a2"), B = c("b1", "b2")),
            B = 200, seed = 113)
  td <- generics::tidy(h)
  # Two label columns join into one term; the first numeric column is
  # the estimate; the rest pass through under their own names.
  expect_identical(td$term, "A:B")
  expect_identical(td$estimate, h$htmt)
  expect_identical(td$upper_limit, h$upper_limit)
  gl <- generics::glance(h)
  expect_identical(gl$n_terms, 1L)
  expect_identical(gl$conf_level, 0.95)

  # suppressWarnings: the deliberately toy 3-indicator model trips lavaan's
  # post-fit checks; the tidier, not the fit quality, is under test.
  fit <- suppressWarnings(lavaan::cfa("f =~ a1 + a2 + b1", data = d, std.lv = TRUE))
  a <- average_variance_extracted(fit)
  ta <- generics::tidy(a)
  expect_identical(ta$term, as.character(a$factor))
  expect_identical(ta$estimate, a$ave)
  ga <- generics::glance(a)
  expect_identical(ga$n_terms, nrow(a))
})

test_that("no tidy() or glance() column anywhere carries a dot", {
  # The package-wide decision of 2026-07-30: the broom verbs speak
  # DMAR's native underscore names.
  set.seed(113)
  objs <- list(
    ci_smd(smd = .5, n_1 = 40, n_2 = 40),
    welch_t(x = rnorm(20), y = rnorm(20) + .3)
  )
  for (o in objs) {
    expect_false(any(grepl(".", names(generics::tidy(o)), fixed = TRUE)))
    expect_false(any(grepl(".", names(generics::glance(o)), fixed = TRUE)))
  }
})
