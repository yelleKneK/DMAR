# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

test_that("average_variance_extracted() computes the Fornell-Larcker mean", {
  res <- average_variance_extracted(loadings = c(.8, .7, .6))
  expect_s3_class(res, "dmar_tbl")
  expect_equal(res$ave, mean(c(.8, .7, .6)^2))
})

test_that("average_variance_extracted() on a fit is the mean squared standardized loading", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(113)
  f <- rnorm(250)
  d <- data.frame(y1 = .8 * f + rnorm(250, 0, .6),
                  y2 = .7 * f + rnorm(250, 0, .7),
                  y3 = .6 * f + rnorm(250, 0, .8))
  fit <- lavaan::cfa("f =~ y1 + y2 + y3", data = d, std.lv = TRUE)
  res <- average_variance_extracted(fit)
  expect_equal(res$factor, "f")
  # The interval columns are present, NA, when no interval is formed.
  expect_named(res, c("factor", "ave", "ci_lower", "ci_upper"))
  expect_true(is.na(res$ci_lower) && is.na(res$ci_upper))
  expect_null(attr(res, "conf_level"))

  # The oracle is the Fornell and Larcker (1981) definition computed from
  # scratch off the fitted solution: the mean of the squared completely
  # standardized loadings. Deriving it here rather than calling another
  # package's AVE() keeps the assertion a test of the definition, and it
  # holds to machine precision rather than to a third decimal.
  std <- lavaan::standardizedSolution(fit)
  lam <- std$est.std[std$op == "=~" & std$lhs == "f"]
  expect_equal(res$ave, mean(lam^2))
})

test_that("average_variance_extracted() validates", {
  expect_error(average_variance_extracted(), "exactly one")
  expect_error(average_variance_extracted(loadings = .8), "two or more")
  expect_error(average_variance_extracted(fit = lm(mpg ~ wt, mtcars)),
               "lavaan fit")
})

test_that("AVE retains its point and counts while applying both bootstrap floors", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa("f =~ t6_paragraph_comprehension+t7_sentence+t9_word_meaning",
                     data = holzinger_swineford)
  state <- new.env(parent = emptyenv()); state$returned <- 120L
  local_mocked_bindings(lavBootstrap = function(object, r, fun, ...) {
    x <- matrix(NA_real_, r, 2L)
    x[seq_len(state$returned), 1L] <- seq(.3, .7, length.out = state$returned)
    x[seq_len(state$returned), 2L] <- 0
    x[1L, 2L] <- 1
    x
  }, .package = "lavaan")
  expect_warning(x <- average_variance_extracted(fit, ci_method = "percentile",
    B = 300, seed = 113), "at least 150 are required")
  expect_true(is.finite(x$ave))
  expect_true(is.na(x$ci_lower) && is.na(x$ci_upper))
  expect_equal(attr(x, "B_used"), 120)
  expect_equal(attr(x, "bootstrap")$improper_kept, 1)
  state$returned <- 100L
  expect_warning(y <- average_variance_extracted(fit, ci_method = "percentile",
    B = 200, seed = 113), "were set aside")
  expect_equal(y$ci_lower, unname(quantile(seq(.3, .7, length.out = 100), .025)))
  expect_equal(y$ci_upper, unname(quantile(seq(.3, .7, length.out = 100), .975)))
  expect_equal(attr(y, "bootstrap")$minimum_returned, 100)
})

test_that("a loading fixed or constrained to zero is not a loading, so the model written with it gets the AVE of the model without it", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # Through 1.0.1 a row held at zero entered the average as a loading of
  # zero: a visual factor of holzinger_swineford written with
  # 0*t7_sentence, or with a*t7_sentence and a == 0, the same model as the
  # one without that row, had an AVE of .279 for its three tests rather
  # than .371.
  skip_if_not_installed("lavaan")
  hs <- holzinger_swineford
  plain <- lavaan::cfa("
    vis  =~ t1_visual_perception + t2_cubes + t4_lozenges
    verb =~ t5_general_information + t6_paragraph_comprehension + t7_sentence",
    data = hs, std.lv = TRUE)
  zero <- lavaan::cfa("
    vis  =~ t1_visual_perception + t2_cubes + t4_lozenges + 0*t7_sentence
    verb =~ t5_general_information + t6_paragraph_comprehension + t7_sentence",
    data = hs, std.lv = TRUE)
  constrained <- lavaan::cfa("
    vis  =~ t1_visual_perception + t2_cubes + t4_lozenges + a*t7_sentence
    verb =~ t5_general_information + t6_paragraph_comprehension + t7_sentence
    a == 0", data = hs, std.lv = TRUE)
  a <- average_variance_extracted(plain)
  b <- average_variance_extracted(zero)
  e <- average_variance_extracted(constrained)
  expect_equal(b$ave, a$ave, tolerance = 1e-6)
  expect_equal(e$ave, a$ave, tolerance = 1e-6)
  expect_equal(b$ave[b$factor == "vis"], 0.3715, tolerance = 1e-3)
  # The verbal factor written with 0*t1_visual_perception beside the same
  # visual factor: .529 through 1.0.1 (the zero averaged in as a fourth
  # loading, .705 times 3/4), and now the .705 of the model without it.
  verb0 <- average_variance_extracted(lavaan::cfa("
    vis  =~ t1_visual_perception + t2_cubes + t4_lozenges
    verb =~ t5_general_information + t6_paragraph_comprehension +
            t7_sentence + 0*t1_visual_perception", data = hs, std.lv = TRUE))
  expect_equal(verb0$ave, a$ave, tolerance = 1e-6)
  expect_equal(verb0$ave[verb0$factor == "verb"], 0.7048, tolerance = 1e-3)
  # Nor does a loading held at zero make a factor improper. With a visual
  # factor of visual perception and lozenges alone, visual perception's
  # error variance is negative, so the visual factor's AVE is NA; written
  # with 0*t1_visual_perception, the verbal factor keeps the .705 of the
  # model without that row, where the check for an improper solution
  # counted visual perception as one of its indicators.
  het <- "vis  =~ t1_visual_perception + t4_lozenges
          verb =~ t5_general_information + t6_paragraph_comprehension +
                  t7_sentence"
  het0 <- paste(het, "+ 0*t1_visual_perception")
  h <- suppressWarnings(lavaan::cfa(het, data = hs, std.lv = TRUE))
  h0 <- suppressWarnings(lavaan::cfa(het0, data = hs, std.lv = TRUE))
  w <- testthat::capture_warnings(r0 <- average_variance_extracted(h0))
  r <- suppressWarnings(average_variance_extracted(h))
  expect_length(w, 1L)
  expect_match(w, "^The AVE of vis is NA")
  expect_equal(r0$ave, r$ave, tolerance = 1e-6)
  expect_equal(r0$ave[r0$factor == "verb"], 0.705, tolerance = 1e-3)
})

test_that("the percentile bootstrap brackets the estimate and reproduces", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # 2 x 300 CFA refits; the cheap-path tests cover CRAN
  set.seed(113)
  f <- rnorm(250)
  d <- data.frame(y1 = .8 * f + rnorm(250, 0, .6),
                  y2 = .7 * f + rnorm(250, 0, .7),
                  y3 = .6 * f + rnorm(250, 0, .8))
  fit <- lavaan::cfa("f =~ y1 + y2 + y3", data = d, std.lv = TRUE)

  plain <- average_variance_extracted(fit)
  expect_true(all(is.na(plain$ci_lower)) && all(is.na(plain$ci_upper)))

  b1 <- average_variance_extracted(fit, ci_method = "percentile",
                                   B = 300, seed = 113)
  b2 <- average_variance_extracted(fit, ci_method = "percentile",
                                   B = 300, seed = 113)
  expect_identical(as.data.frame(b1), as.data.frame(b2))
  expect_equal(b1$ave, plain$ave)
  expect_lt(b1$ci_lower, b1$ave)
  expect_gt(b1$ci_upper, b1$ave)
  expect_identical(attr(b1, "B_used"), 300L)
  expect_identical(attr(b1, "conf_level"), 0.95)

  # One shape whatever ci_method: the same columns, and the broom verbs
  # answer with the same columns, glance() reporting the resamples used.
  expect_identical(names(b1), names(plain))
  expect_identical(names(generics::tidy(b1)), names(generics::tidy(plain)))
  expect_identical(names(generics::glance(b1)),
                   names(generics::glance(plain)))
  expect_true(is.na(generics::glance(plain)$B_used))
  expect_identical(generics::glance(b1)$B_used, 300L)
})

test_that("loadings alone cannot carry an interval", {
  expect_error(
    average_variance_extracted(loadings = c(.8, .7, .6),
                               ci_method = "percentile"),
    "loadings' alone")
  # The loadings path keeps the canonical four columns.
  r <- average_variance_extracted(loadings = c(.8, .7, .6))
  expect_named(r, c("factor", "ave", "ci_lower", "ci_upper"))
})

test_that("a factor with a negative error variance has no AVE, as in fa_confirmatory()", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # The stage 2b review: average_variance_extracted(fit) reported an AVE of
  # .726 for a factor of visual perception and lozenges beside a verbal
  # factor of holzinger_swineford, where visual perception's error
  # variance is negative and its communality 1.304, while fa_confirmatory() gives
  # the factor's ave row as NA with a warning. It now does the same.
  skip_if_not_installed("lavaan")
  model <- "f1 =~ t1_visual_perception + t4_lozenges
            f2 =~ t5_general_information + t6_paragraph_comprehension +
                  t7_sentence"
  fit <- suppressWarnings(lavaan::cfa(model, data = holzinger_swineford))
  w <- testthat::capture_warnings(res <- average_variance_extracted(fit))
  expect_length(w, 1L)
  expect_identical(w, paste0(
    "The AVE of f1 is NA: the solution is improper (a Heywood case): ",
    "t1_visual_perception has a negative error variance (its communality ",
    "is 1.304). The average variance extracted averages the communalities, ",
    "so it is not a proportion of variance here. A larger sample, or ",
    "reconsidering the model for these items, are the usual remedies."))
  expect_identical(res$factor, c("f1", "f2"))
  expect_true(is.na(res$ave[1]))
  # The proper factor keeps its AVE, the mean of its squared standardized
  # loadings.
  std <- lavaan::standardizedSolution(fit)
  lam <- std$est.std[std$op == "=~" & std$lhs == "f2"]
  expect_equal(res$ave[2], mean(lam^2))
  # The same factors, rows, and warning as the measurement table of
  # fa_confirmatory(), whatever the identification of the fit.
  wk <- testthat::capture_warnings(k <- fa_confirmatory(holzinger_swineford, list(
    f1 = c("t1_visual_perception", "t4_lozenges"),
    f2 = c("t5_general_information", "t6_paragraph_comprehension",
           "t7_sentence")), output = "measurement"))
  kk <- as.data.frame(k)
  expect_true(is.na(kk$estimate[kk$term == "ave_f1"]))
  expect_equal(res$ave[2], kk$estimate[kk$term == "ave_f2"],
               tolerance = 1e-6)
  expect_true(any(grepl("its communality is 1.304", wk, fixed = TRUE)))
  fit_std <- suppressWarnings(lavaan::cfa(model, data = holzinger_swineford,
                                          std.lv = TRUE))
  w2 <- testthat::capture_warnings(res2 <- average_variance_extracted(fit_std))
  expect_identical(w2, w)
  expect_equal(as.data.frame(res2), as.data.frame(res), tolerance = 1e-6)
})

test_that("an ordered factor whose latent response has a negative residual variance has no AVE", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # The improper ordered sample of test-fa_confirmatory.R: y1's latent response has
  # a residual variance of -0.057 under the delta parameterization, so its
  # communality is 1.057.
  skip_if_not_installed("lavaan")
  set.seed(1)
  lam <- c(.95, .95, .6, .6)
  f <- stats::rnorm(80)
  ys <- sapply(lam, function(l) l * f + sqrt(1 - l^2) * stats::rnorm(80))
  it <- apply(ys, 2, function(y) findInterval(y, c(-1, 0, 1)) + 1L)
  colnames(it) <- paste0("y", 1:4)
  fit <- suppressWarnings(lavaan::cfa("G =~ y1 + y2 + y3 + y4",
                                      data = as.data.frame(it),
                                      ordered = TRUE, std.lv = TRUE))
  w <- testthat::capture_warnings(res <- average_variance_extracted(fit))
  expect_true(is.na(res$ave))
  expect_match(w, paste0("^The AVE of G is NA: the solution is improper ",
                         "\\(a Heywood case\\): the latent response behind ",
                         "y1 has a negative residual variance \\(its ",
                         "communality is 1\\.057\\)"))
  # The helper reports the indicator and its communality.
  imp <- .ave_improper(fit)
  expect_named(imp, "G")
  expect_equal(as.numeric(imp$G), 1.057, tolerance = 1e-3)
  expect_identical(attr(imp$G, "kind"), "ordered")
  # A proper solution is untouched, with no warning.
  ok <- lavaan::cfa("G =~ y2 + y3 + y4", data = as.data.frame(it),
                    ordered = TRUE, std.lv = TRUE)
  expect_identical(.ave_improper(ok), list())
  expect_no_warning(average_variance_extracted(ok))
})

test_that("the bootstrap gives no interval for an improper factor", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # 150 CFA refits
  fit <- suppressWarnings(lavaan::cfa(
    "f1 =~ t1_visual_perception + t4_lozenges
     f2 =~ t5_general_information + t6_paragraph_comprehension +
           t7_sentence", data = holzinger_swineford))
  # A few refits of this unstable two-indicator factor do not converge
  # and are dropped with their own warning.
  w <- testthat::capture_warnings(
    res <- average_variance_extracted(fit, ci_method = "percentile",
                                      B = 150, seed = 113))
  expect_true(any(grepl("^The AVE of f1 is NA", w)))
  expect_true(all(is.na(unlist(res[1, c("ave", "ci_lower", "ci_upper")]))))
  expect_true(all(is.finite(unlist(res[2, c("ave", "ci_lower", "ci_upper")]))))
  expect_lt(res$ci_lower[2], res$ave[2])
  expect_gt(res$ci_upper[2], res$ave[2])
})

test_that("the bootstrap counts the improper replications it keeps, per factor, and prints them", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # 2 x 200 CFA refits
  # The 1.1.0 review: a replication whose resample solution is improper for
  # a factor that is proper in the observed fit was kept, as the reliability
  # family's improper = "keep" keeps one, but neither counted nor reported.
  # On 80 pupils of holzinger_swineford the spatial factor is proper in the
  # sample, yet about half its replications are improper, and its upper
  # limit, 1.19, is not a proportion of variance.
  set.seed(113)
  d <- holzinger_swineford[sample(301, 80), ]
  fit <- suppressWarnings(lavaan::cfa(
    "verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
     spatial =~ t1_visual_perception + t2_cubes + t4_lozenges", data = d))
  expect_identical(.ave_improper(fit), list())
  w <- testthat::capture_warnings(
    res <- average_variance_extracted(fit, ci_method = "percentile",
                                      B = 200, seed = 113))
  # No interval changes: the limits of 1.0.x, to the last digit.
  expect_identical(res$factor, c("spatial", "verbal"))
  expect_equal(res$ci_lower, c(0.371596507609394, 0.652123670128383),
               tolerance = 1e-12)
  expect_equal(res$ci_upper, c(1.186703359303767, 0.806398145268667),
               tolerance = 1e-12)
  # One warning for the bootstrap, as the reliability family gives: the
  # replications set aside, for the reason the line under the table names,
  # and the caution, naming the factor it holds for (see the test below).
  expect_length(w, 1L)
  expect_match(w, paste0("^average_variance_extracted\\(\\): 2 of 200 ",
                         "bootstrap replications were set aside ",
                         "\\(2 failed to refit or did not converge\\)"))
  bs <- attr(res, "bootstrap")
  expect_identical(bs[c("attempted", "returned", "failed", "improper",
                        "improper_kept", "refit_failed", "improper_rule",
                        "caution")],
                   list(attempted = 200L, returned = 198L, failed = 2L,
                        improper = 0L, improper_kept = 101L,
                        refit_failed = 2L, improper_rule = "keep",
                        caution = TRUE))
  expect_identical(bs$reasons, c(error = 0L, non_finite = 2L))
  expect_identical(bs$improper_kept_by_factor,
                   c(spatial = 101L, verbal = 0L))
  expect_identical(attr(res, "B_used"), 198L)
  # Regression test. The caution, printed under the table only in the
  # first versions of 1.1.0, is given at call time as well, naming the
  # factor and its count: the 101 improper replications kept for spatial
  # and the 2 set aside are 103 of 200, against the 5 beyond each limit.
  # Through 1.0.1 the warning said that 2 replications "failed to refit
  # (or did not converge) and were dropped", and nothing more.
  expect_identical(w, paste0(
    "average_variance_extracted(): 2 of 200 bootstrap replications ",
    "were set aside (2 failed to refit or ",
    "did not converge); the intervals are computed from the 198 that ",
    "returned a value, 101 of them from an improper (Heywood) solution and ",
    # The count rule is held to the request (2 here, under
    # local_small_B()), and the warning says so.
    "kept; the 2 required are counted on the number requested (200), not ",
    "on the number that remained; request more than 2 where replications ",
    "can fail. Caution: for spatial, 103 of 200 replications (51.5%) came from ",
    "a solution improper for the factor or were set aside, at least the ",
    "2.5% of the bootstrap distribution that lies beyond each limit of a ",
    "95% interval, so the limits of its interval depend on how those ",
    "replications were treated; report the interval with these counts."))
  # An independent recount over the same resamples: lavaan's bootstrap at
  # the same seed, asking each refit whether an indicator of the factor has
  # a negative error variance.
  neg <- function(f) {
    pe <- lavaan::parameterEstimates(f, se = FALSE, zstat = FALSE,
                                     pvalue = FALSE, ci = FALSE)
    v <- pe[pe$op == "~~" & pe$lhs == pe$rhs, ]
    c(spatial = any(v$est[v$lhs %in% c("t1_visual_perception", "t2_cubes",
                                       "t4_lozenges")] < 0),
      verbal = any(v$est[v$lhs %in% c("t6_paragraph_comprehension",
                                      "t7_sentence", "t9_word_meaning")] < 0))
  }
  oracle <- withr::with_seed(113, suppressWarnings(
    lavaan::lavBootstrap(fit, r = 200L, fun = neg)))
  oracle <- oracle[stats::complete.cases(oracle), , drop = FALSE]
  expect_identical(nrow(oracle), 198L)
  expect_identical(as.integer(colSums(oracle == 1)), c(101L, 0L))
  # The family's footer line, then the count for each factor, at the
  # console and in a knitted report alike.
  lines <- DMAR:::.dmar_tbl_footer(res)
  expect_identical(unname(lines[names(lines) == "bootstrap"]), c(
    paste0("Bootstrap: 200 replications attempted, 198 returned a value ",
           "(2 set aside: 2 failed to refit or did not converge); 101 of ",
           "the 198 came from an improper (Heywood) solution and were kept; ",
           "caution: for spatial, 51.5% of the replications came from a ",
           "solution improper for the factor or were set aside, at least ",
           "the 2.5% beyond each limit"),
    paste0("Improper replications kept, by factor (a resample in which an ",
           "indicator of the factor had a negative error variance, so ",
           "that its AVE averages a communality above one): spatial 101, ",
           "verbal 0")))
  expect_output(print(res), "spatial 101, verbal 0", fixed = TRUE)
  # The line and the warning give the replications set aside the same
  # reason, which the line once called "a non-finite value".
  expect_match(lines[["bootstrap"]],
               "(2 set aside: 2 failed to refit or did not converge)",
               fixed = TRUE)
  expect_match(w, "(2 failed to refit or did not converge)", fixed = TRUE)
  skip_if_not_installed("knitr")
  expect_match(as.character(knitr::knit_print(res)),
               "101 of the 198 came from an improper", fixed = TRUE)
})

test_that("the bootstrap counts follow the reliability family's attribute", {
  # The counts behind the footer line, from a matrix of improper flags (the
  # replications used by the factors): a factor improper in the observed
  # fit has no interval but remains counted, a replication improper for two
  # factors counts once in improper_kept, and the caution holds when the
  # replications set aside and the improper ones kept for some factor reach
  # the tail area beyond a limit, here 2.5 of 100.
  flags <- cbind(a = c(TRUE, TRUE, FALSE, FALSE, rep(FALSE, 94)),
                 b = c(TRUE, FALSE, TRUE, FALSE, rep(FALSE, 94)),
                 c = rep(TRUE, 98))
  bs <- .ave_bootstrap_counts(flags, c("a", "b", "c"), c("a", "b"),
                              B = 100, conf_level = 0.95)
  expect_identical(bs$improper_kept_by_factor, c(a = 2L, b = 2L, c = 98L))
  expect_identical(bs$improper_kept, 98L)
  expect_identical(bs[c("attempted", "returned", "failed")],
                   list(attempted = 100L, returned = 98L, failed = 2L))
  expect_identical(bs$reasons, c(error = 0L, non_finite = 2L))
  # 2 set aside + 2 kept for a reach 2.5: the caution holds.
  expect_true(bs$caution)
  expect_false(.ave_bootstrap_counts(flags[, c("a", "b")] & FALSE,
                                     c("a", "b"), c("a", "b"), B = 98,
                                     conf_level = 0.95)$caution)
  # The by factor line needs two or more factors with an interval and a
  # replication to report; the family's line alone covers one factor.
  tbl <- function(by) {
    structure(data.frame(factor = names(by)),
              bootstrap = list(improper_kept_by_factor = by))
  }
  expect_match(.ave_bootstrap_footer(tbl(c(a = 2L, b = 0L))),
               ": a 2, b 0$")
  expect_null(.ave_bootstrap_footer(tbl(c(a = 2L))))
  expect_null(.ave_bootstrap_footer(tbl(c(a = 0L, b = 0L))))
  # The class carries the footer method on every path, so the shape of
  # the result does not depend on ci_method.
  expect_s3_class(average_variance_extracted(loadings = c(.8, .7, .6)),
                  c("dmar_ave", "dmar_tbl", "data.frame"), exact = TRUE)
  expect_null(attr(average_variance_extracted(loadings = c(.8, .7, .6)),
                   "bootstrap"))
})

test_that("the bootstrap warning gives the caution for each factor it holds for", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # The warning of a percentile interval, from the counts, without a refit:
  # the replications set aside with the reason the line under the table
  # gives, and the caution for each factor whose replications set aside and
  # improper ones kept reach the tail area beyond a limit (2.5 of 100 at
  # 95 percent), named with its count.
  flags <- cbind(a = c(TRUE, TRUE, FALSE, FALSE, rep(FALSE, 94)),
                 b = c(TRUE, FALSE, TRUE, FALSE, rep(FALSE, 94)),
                 c = c(rep(FALSE, 3), TRUE, rep(FALSE, 94)))
  bs <- .ave_bootstrap_counts(flags, c("a", "b", "c"), c("a", "b", "c"),
                              B = 100, conf_level = 0.95)
  expect_identical(bs$refit_failed, 2L)
  expect_identical(.ave_bootstrap_warning(bs, 0.95), paste0(
    "2 of 100 bootstrap replications were set aside (2 failed to refit or ",
    "did not converge); the intervals are computed from the 98 that ",
    "returned a value, 4 of them from an improper (Heywood) solution and ",
    "kept. Caution: for a, 4 of 100 replications (4%), for b, 4 of 100 ",
    "(4%), and for c, 3 of 100 (3%) came from a solution improper for the ",
    "factor or were set aside, at least the 2.5% of the bootstrap ",
    "distribution that lies beyond each limit of a 95% interval, so the ",
    "limits of their intervals depend on how those replications were ",
    "treated; report the intervals with these counts."))
  # Each share is formatted on its own (20 and 15 of 200 are 10% and
  # 7.5%, where the vector once gave "10.0%").
  bs2 <- bs
  bs2$attempted <- 200L
  bs2$failed <- 0L
  bs2$improper_kept_by_factor <- c(a = 20L, b = 15L)
  expect_match(.ave_bootstrap_warning(bs2, 0.95),
               "for a, 20 of 200 replications (10%) and for b, 15 of 200 (7.5%)",
               fixed = TRUE)
  # At 90 percent the tail area is 5 of 100, which no factor reaches.
  expect_identical(
    .ave_bootstrap_warning(bs, 0.90),
    paste0("2 of 100 bootstrap replications were set aside (2 failed to ",
           "refit or did not converge); the intervals are computed from ",
           "the 98 that returned a value, 4 of them from an improper ",
           "(Heywood) solution and kept."))
  # With nothing set aside, the caution alone, for the one factor it holds
  # for; with no caution either, no warning.
  flags0 <- cbind(a = c(TRUE, TRUE, TRUE, rep(FALSE, 97)),
                  b = rep(TRUE, 100),
                  c = c(TRUE, rep(FALSE, 99)))
  bs0 <- .ave_bootstrap_counts(flags0, c("a", "b", "c"), c("a", "c"),
                               B = 100, conf_level = 0.95)
  expect_true(bs0$caution)
  expect_identical(.ave_bootstrap_warning(bs0, 0.95), paste0(
    "Caution: for a, 3 of 100 replications (3%) came from a solution ",
    "improper for the factor or were set aside, at least the 2.5% of the ",
    "bootstrap distribution that lies beyond each limit of a 95% interval, ",
    "so the limits of its interval depend on how those replications were ",
    "treated; report the interval with these counts."))
  bs1 <- .ave_bootstrap_counts(flags0 & FALSE, c("a", "b", "c"),
                               c("a", "b", "c"), B = 100, conf_level = 0.95)
  expect_false(bs1$caution)
  expect_null(.ave_bootstrap_warning(bs1, 0.95))
})

test_that("no bootstrap is run when every factor is improper", {
  skip_on_cran()  # Additional model fits run in the full local suite.
  # The 1.1.0 review: with every factor improper in the observed fit, the
  # percentile interval still refit the model B times and recorded a
  # confidence level, B_used, and the counts, printing "Confidence level:
  # 95%" and a bootstrap line beneath a table that was all NA. No interval
  # is formed for an improper factor, so nothing is run or recorded, as
  # krippendorff_alpha() forms no interval for an undefined coefficient.
  skip_if_not_installed("lavaan")
  set.seed(1)
  lam <- c(.95, .95, .6, .6)
  f <- stats::rnorm(80)
  ys <- sapply(lam, function(l) l * f + sqrt(1 - l^2) * stats::rnorm(80))
  it <- apply(ys, 2, function(y) findInterval(y, c(-1, 0, 1)) + 1L)
  colnames(it) <- paste0("y", 1:4)
  fit <- suppressWarnings(lavaan::cfa("G =~ y1 + y2 + y3 + y4",
                                      data = as.data.frame(it),
                                      ordered = TRUE, std.lv = TRUE))
  counter <- new.env()
  counter$n <- 0L
  local_mocked_bindings(
    lavBootstrap = function(...) {
      counter$n <- counter$n + 1L
      stop("lavBootstrap() should not be called")
    },
    .package = "lavaan")
  w <- testthat::capture_warnings(
    res <- average_variance_extracted(fit, ci_method = "percentile",
                                      B = 100, seed = 1))
  expect_identical(counter$n, 0L)
  expect_length(w, 1L)
  expect_match(w, "^The AVE of G is NA")
  expect_true(all(is.na(unlist(res[, c("ave", "ci_lower", "ci_upper")]))))
  expect_null(attr(res, "conf_level"))
  expect_null(attr(res, "B_used"))
  expect_null(attr(res, "bootstrap"))
  out <- utils::capture.output(print(res))
  expect_false(any(grepl("Confidence level|Bootstrap", out)))
  # The arguments are still checked.
  expect_error(suppressWarnings(
    average_variance_extracted(fit, ci_method = "percentile", B = 1)),
    "single whole number", class = "dmar_count_invalid")
})

test_that("the caution under the table gives each factor its own share", {
  # The open item of the stage 5 review: the caution's percentage counted
  # the replications set aside with those improper for any factor, while
  # the caution itself, and the warning, count them for each factor. Here
  # 30 of 200 replications were improper for some factor, 20 for a and 15
  # for b: the old line said 15%, which is neither factor's share.
  tbl <- structure(
    data.frame(factor = c("a", "b"), ave = c(.5, .6),
               ci_lower = c(.4, .5), ci_upper = c(.9, .7)),
    class = c("dmar_ave", "dmar_tbl", "data.frame"),
    conf_level = 0.95,
    bootstrap = list(attempted = 200L, returned = 200L, failed = 0L,
                     reasons = c(error = 0L, non_finite = 0L),
                     error_messages = integer(), undetermined = 0L,
                     improper = 0L, improper_kept = 30L, no_solution = 0L,
                     not_converged = 0L, refit_failed = 0L,
                     improper_rule = "keep", caution = TRUE,
                     improper_kept_by_factor = c(a = 20L, b = 15L)))
  line <- DMAR:::.bootstrap_footer_line(tbl)
  expect_identical(line, paste0(
    "Bootstrap: 200 replications attempted, 200 returned a value; 30 of ",
    "the 200 came from an improper (Heywood) solution and were kept; ",
    "caution: for a, 10% of the replications, and for b, 7.5%, came from ",
    "a solution improper for the factor or were set aside, at least the ",
    "2.5% beyond each limit"))
  # A factor below the tail area is not named.
  attr(tbl, "bootstrap")$improper_kept_by_factor <- c(a = 20L, b = 3L)
  expect_match(DMAR:::.bootstrap_footer_line(tbl),
               "caution: for a, 10% of the replications came from a solution improper for the factor",
               fixed = TRUE)
  # One factor with an interval keeps the family's wording.
  attr(tbl, "bootstrap")$improper_kept_by_factor <- c(a = 30L)
  expect_match(DMAR:::.bootstrap_footer_line(tbl),
               "caution: 15% of the replications came from an improper solution or were set aside",
               fixed = TRUE)
})


test_that("AVE footer explains its separate minimum and limits its interval caution", {
  flags <- cbind(a = rep(FALSE, 80), b = rep(TRUE, 80))
  bs <- .ave_bootstrap_counts(flags, c("a", "b"), character(),
                              B = 100, conf_level = .95)
  bs$interval_withheld <- TRUE
  bs$minimum_returned <- 100L
  tbl <- .as_dmar_tbl(data.frame(factor = c("a", "b"), ave = c(.5, NA)))
  attr(tbl, "bootstrap") <- bs
  expect_match(.bootstrap_footer_line(tbl), "at least 100 returned estimates are required, so no interval", fixed = TRUE)
  bs <- .ave_bootstrap_counts(cbind(a = c(rep(TRUE, 10), rep(FALSE, 190)),
                                   b = rep(TRUE, 200)), c("a", "b"), "a",
                              B = 200, conf_level = .95)
  attr(tbl, "bootstrap") <- bs
  attr(tbl, "conf_level") <- .95
  line <- .bootstrap_footer_line(tbl)
  expect_match(line, "caution: 5% of the replications", fixed = TRUE)
  expect_no_match(line, "for b", fixed = TRUE)
  expect_identical(bs$improper_kept_by_factor, c(a = 10L, b = 200L))
})

# The correlation matrix of three items that measure one factor with
# standardized loadings sqrt(r): every correlation is r, so the AVE of a
# fit that reaches the solution is r exactly (each squared standardized
# loading is the correlation of two items).
.ave_equal_r <- function(r, items = paste0("x", 1:3)) {
  S <- matrix(r, length(items), length(items))
  diag(S) <- 1
  dimnames(S) <- list(items, items)
  S
}

test_that("a fit whose optimizer did not converge is refused before anything is read from it", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  S <- .ave_equal_r(.81)
  # Through DMAR 1.0.1 and the first builds of 1.1.0 the fit stopped after
  # one iteration had an AVE of .412, reported as any other. lavaan's own
  # warning that the optimizer stopped is the point of this fit.
  bad <- suppressWarnings(lavaan::cfa("f =~ x1 + x2 + x3", sample.cov = S,
    sample.nobs = 200, sample.cov.rescale = FALSE, std.lv = TRUE,
    control = list(iter.max = 1L)))
  expect_false(lavaan::lavInspect(bad, "converged"))
  e <- expect_error(average_variance_extracted(bad),
                    class = "dmar_not_converged")
  expect_s3_class(e, c("dmar_not_converged", "dmar_refusal", "error"))
  expect_identical(e$iterations, 1L)
  expect_identical(e$fn, "average_variance_extracted()")
  expect_match(conditionMessage(e), "did not converge (lavaan stopped its optimizer after 1 iteration)",
               fixed = TRUE)
  # A bootstrap request is refused the same way, before any resample is
  # drawn.
  local_mocked_bindings(lavBootstrap = function(...) stop("resampled"),
                        .package = "lavaan")
  expect_error(average_variance_extracted(bad, ci_method = "percentile",
                                          B = 400, seed = 113),
               class = "dmar_not_converged")
  # A fit that was never optimized (do.fit = FALSE) is refused as well.
  unfit <- lavaan::cfa("f =~ x1 + x2 + x3", sample.cov = S,
    sample.nobs = 200, sample.cov.rescale = FALSE, std.lv = TRUE,
    do.fit = FALSE)
  e0 <- expect_error(average_variance_extracted(unfit),
                     class = "dmar_not_converged")
  expect_identical(e0$iterations, 0L)
  # The same model fitted to convergence: the oracle is the algebra of
  # equal correlations, an AVE of .81.
  good <- lavaan::cfa("f =~ x1 + x2 + x3", sample.cov = S, sample.nobs = 200,
    sample.cov.rescale = FALSE, std.lv = TRUE)
  expect_equal(average_variance_extracted(good)$ave, .81, tolerance = 1e-6)
})

test_that("a converged model with equal raw loadings keeps its AVE", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # Four items with raw loadings of .7 held equal and error variances .1,
  # .2, 1, and 2, fitted to the population covariance matrix: the
  # standardized loadings are .7 / sqrt(.49 + theta), so the oracle is the
  # mean of .49 / (.49 + theta), .5165749.
  items <- paste0("y", 1:4)
  theta <- c(.1, .2, 1, 2)
  S <- matrix(.49, 4, 4) + diag(theta)
  dimnames(S) <- list(items, items)
  fit <- lavaan::cfa("f =~ a*y1 + a*y2 + a*y3 + a*y4", sample.cov = S,
    sample.nobs = 500, sample.cov.rescale = FALSE, std.lv = TRUE)
  expect_true(lavaan::lavInspect(fit, "converged"))
  expect_equal(average_variance_extracted(fit)$ave, mean(.49 / (.49 + theta)),
               tolerance = 1e-6)
})

test_that("a multiple-group fit is refused, and each group's AVE is that of its own loadings", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  S_A <- .ave_equal_r(.81)
  S_B <- .ave_equal_r(.09)
  mg <- function(n) lavaan::cfa("f =~ x1 + x2 + x3",
    sample.cov = list(A = S_A, B = S_B), sample.nobs = n,
    sample.cov.rescale = FALSE, std.lv = TRUE)
  # Through DMAR 1.0.1 and the first builds of 1.1.0 both fits gave one
  # unlabeled row of .450, the mean of the six squared standardized
  # loadings, for a factor whose AVE is .81 in group A and .09 in group B.
  for (n in list(c(1000, 100), c(100, 1000))) {
    fit <- mg(n)
    expect_true(lavaan::lavInspect(fit, "converged"))
    e <- expect_error(average_variance_extracted(fit),
                      class = "dmar_argument_invalid")
    expect_identical(e$arg, "fit")
    expect_identical(e$groups, c("A", "B"))
    expect_match(conditionMessage(e), "multiple-group fit (2 groups: A, B)",
                 fixed = TRUE)
  }
  # A bootstrap request meets the same refusal, before any resample.
  local_mocked_bindings(lavBootstrap = function(...) stop("resampled"),
                        .package = "lavaan")
  expect_error(average_variance_extracted(fit, ci_method = "percentile",
                                          B = 400, seed = 113),
               class = "dmar_argument_invalid")
  # The route the refusal names: one group's standardized loadings given
  # as `loadings`. The oracle is the algebra of equal correlations.
  std <- lavaan::standardizedSolution(fit)
  per_group <- vapply(1:2, function(g) {
    average_variance_extracted(
      loadings = std$est.std[std$op == "=~" & std$group == g])$ave
  }, numeric(1))
  expect_equal(per_group, c(.81, .09), tolerance = 1e-6)
})

test_that("a multilevel fit is refused, whatever its factors are named", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two two-level fits
  set.seed(113)
  J <- 150L; m <- 8L
  cl <- rep(seq_len(J), each = m)
  fb <- rnorm(J)[cl]
  fw <- rnorm(J * m)
  lw <- c(.7, .6, .5); lb <- c(.5, .45, .4); ub <- c(.15, .2, .25)
  Y <- vapply(1:3, function(j) lb[j] * fb + rnorm(J, 0, ub[j])[cl] +
                lw[j] * fw + rnorm(J * m, 0, sqrt(1 - lw[j]^2)), numeric(J * m))
  colnames(Y) <- paste0("y", 1:3)
  d <- data.frame(Y, cluster = cl)
  same <- lavaan::cfa("level: 1
                         f =~ y1 + y2 + y3
                       level: 2
                         f =~ y1 + y2 + y3", data = d, cluster = "cluster",
                      std.lv = TRUE)
  apart <- lavaan::cfa("level: 1
                          fw =~ y1 + y2 + y3
                        level: 2
                          fb =~ y1 + y2 + y3", data = d, cluster = "cluster",
                       std.lv = TRUE)
  for (fit in list(same, apart)) {
    expect_true(lavaan::lavInspect(fit, "converged"))
    e <- expect_error(average_variance_extracted(fit),
                      class = "dmar_argument_invalid")
    expect_identical(e$levels, c("within", "cluster"))
    expect_match(conditionMessage(e), "multilevel fit (2 levels: within, cluster)",
                 fixed = TRUE)
  }
  # Through DMAR 1.0.1 and the first builds of 1.1.0 the fit with the same
  # factor name at both levels gave one row of .591, the mean of the six
  # squared standardized loadings, where the AVE is .356 within clusters
  # and .826 between them, each level's own lambda matrix being the oracle.
  lam <- lapply(lavaan::lavInspect(same, "std"), function(b) b$lambda[, "f"])
  expect_equal(unname(vapply(lam, function(l) mean(l^2), numeric(1))),
               c(0.3557764, 0.8256386), tolerance = 1e-4)
  expect_equal(average_variance_extracted(loadings = lam$within)$ave,
               mean(lam$within^2))
})
