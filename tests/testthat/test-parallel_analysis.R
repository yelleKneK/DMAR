# The reference data sets in this file are mechanics tests run with a small
# B; the rule on the number of replications is tested in
# test-dmar_bootstrap.R.
local_small_B()

# Six items cut into four ordered categories, n = 400 (the fixture of
# test-factor_count.R, generated here from the same draws).
pa_ordered <- function() {
  withr::with_seed(113, {
    n <- 400L
    g <- rnorm(n)
    groups <- matrix(rnorm(n * 3), n, 3)
    X <- sapply(seq_len(6), function(j) .5 * g + .4 * groups[, ceiling(j / 3)] +
                  sqrt(.59) * rnorm(n))
    colnames(X) <- paste0("x", seq_len(6))
    as.data.frame(lapply(as.data.frame(X), function(v)
      as.integer(cut(v, breaks = c(-Inf, -.7, 0, .7, Inf)))))
  })
}
# The same items with 40 values removed (39 cases incomplete).
pa_ordered_incomplete <- function() {
  ord <- pa_ordered()
  withr::with_seed(7, {
    ord[cbind(sample.int(nrow(ord), 40, replace = TRUE),
              sample.int(6, 40, replace = TRUE))] <- NA
  })
  ord
}
# The columns of a table without its attributes.
pa_columns <- function(x) {
  y <- as.data.frame(x)
  attributes(y) <- attributes(y)[c("names", "row.names")]
  y
}

test_that("complete ordered items with 'missing' left out are analyzed listwise", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # reference data sets are simulated
  # fa_exploratory() and factor_count() resolve the default this way
  # (.efa_resolve_missing()); through the first builds of 1.1.0 a direct
  # call refused these items ("requires explicit pairwise or listwise
  # handling"), so the same items ran in factor_count() and stopped here.
  ord <- pa_ordered()
  expect_false(anyNA(ord))
  omitted <- parallel_analysis(ord, ordered = TRUE, B = 5, seed = 113)
  listwise <- parallel_analysis(ord, ordered = TRUE, missing = "listwise",
                                B = 5, seed = 113)
  expect_identical(pa_columns(omitted), pa_columns(listwise))
  expect_identical(attr(omitted, "method"), attr(listwise, "method"))
  expect_identical(attr(omitted, "method")$missing, "listwise")
  expect_identical(attr(omitted, "counts")[c("attempted", "usable", "failed")],
                   list(attempted = 5L, usable = 5L, failed = 0L))
  # The observed eigenvalues: lavaan's polychoric correlations with the
  # communalities of the one-factor maximum likelihood solution of
  # stats::factanal() on the diagonal.
  R <- as.matrix(lavaan::lavCor(ord, ordered = names(ord), output = "cor"))
  f1 <- stats::factanal(covmat = list(cov = R, n.obs = nrow(ord)), factors = 1)
  diag(R) <- 1 - f1$uniquenesses
  expect_equal(omitted$estimate,
               eigen(R, symmetric = TRUE, only.values = TRUE)$values,
               tolerance = 1e-6)
  # The reference quantiles, written from the description of the ordinal
  # normal reference: after the seed, B blocks of N x p standard normal
  # deviates, each item cut at the normal quantiles of its observed
  # cumulative category proportions, polychoric correlations, the reduced
  # eigenvalues, and the 95th percentile of each rank (quantile type 7).
  reference <- withr::with_seed(113, {
    N <- nrow(ord)
    p <- ncol(ord)
    draws <- t(vapply(seq_len(5), function(b) {
      x <- matrix(stats::rnorm(N * p), N, p)
      for (j in seq_len(p)) {
        prop <- as.numeric(table(ord[[j]])) / N
        cuts <- stats::qnorm(c(0, cumsum(prop)[-length(prop)], 1))
        x[, j] <- sort(unique(ord[[j]]))[findInterval(x[, j], cuts, all.inside = TRUE)]
      }
      x <- as.data.frame(x)
      names(x) <- names(ord)
      Rb <- as.matrix(lavaan::lavCor(x, ordered = names(x), output = "cor"))
      fb <- stats::factanal(covmat = list(cov = Rb, n.obs = N), factors = 1)
      diag(Rb) <- 1 - fb$uniquenesses
      eigen(Rb, symmetric = TRUE, only.values = TRUE)$values
    }, numeric(p)))
    apply(draws, 2L, stats::quantile, probs = 0.95, names = FALSE, type = 7)
  })
  expect_equal(omitted$reference_quantile, reference, tolerance = 1e-6)
  # The parallel analysis factor_count() runs on the same items with the
  # same seed is this one.
  fc <- factor_count(ord, estimator = "WLSMV", ordered = TRUE, B = 5,
                     seed = 113)
  expect_identical(pa_columns(attr(fc, "parallel_analysis")),
                   pa_columns(omitted))
})

test_that("incomplete ordered items with 'missing' left out are refused as fa_exploratory() refuses them", {
  skip_if_not_installed("lavaan")
  ord_na <- pa_ordered_incomplete()
  ord <- pa_ordered()
  expect_identical(sum(!stats::complete.cases(ord_na)), 39L)
  # The rule is fa_exploratory()'s (.efa_resolve_missing()): every call
  # refused below is refused by fa_exploratory(estimator = "WLSMV") on the
  # same items with the same 'missing', with the same class.
  efa_refuses <- function(d, ...) {
    e <- tryCatch(fa_exploratory(d, 1, estimator = "WLSMV", ordered = TRUE, ...),
                  error = function(e) e)
    inherits(e, "dmar_argument_invalid")
  }
  expect_true(efa_refuses(ord_na))
  expect_true(efa_refuses(ord_na, missing = "fiml"))
  expect_true(efa_refuses(ord_na, missing = "ml"))
  expect_true(efa_refuses(ord, missing = "fiml"))
  # Nothing may be computed before the refusal: no correlation, no
  # reference data set, and no draw from the random number stream.
  reached <- function(...) stop("reached before the refusal")
  local_mocked_bindings(.pa_cor = reached, cov_fiml = reached,
                        .pa_reduced_eigen = reached)
  for (call in list(
    quote(parallel_analysis(ord_na, ordered = TRUE, B = 200)),
    quote(parallel_analysis(ord_na, ordered = TRUE, missing = "fiml", B = 200)),
    quote(parallel_analysis(ord_na, ordered = TRUE, missing = "ml", B = 200)),
    # An explicit FIML request is refused for complete ordered items too,
    # as fa_exploratory() refuses it.
    quote(parallel_analysis(ord, ordered = TRUE, missing = "fiml", B = 200)))) {
    set.seed(99)
    before <- get(".Random.seed", envir = globalenv())
    e <- tryCatch(eval(call), error = function(e) e)
    expect_s3_class(e, c("dmar_argument_invalid", "dmar_refusal", "error"))
    expect_identical(e$arg, "missing")
    # The sentence is this function's: parallel_analysis() has no
    # 'estimator' argument, so the message does not name one.
    expect_match(conditionMessage(e), "^Ordered parallel analysis requires")
    expect_match(conditionMessage(e),
                 "explicit missing = \"pairwise\" or \"listwise\" choice",
                 fixed = TRUE)
    expect_no_match(conditionMessage(e), "estimator", fixed = TRUE)
    expect_no_match(conditionMessage(e), "EFA", fixed = TRUE)
    expect_identical(get(".Random.seed", envir = globalenv()), before)
  }
})

test_that("incomplete ordered items run under an explicit pairwise or listwise choice", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # reference data sets are simulated
  ord_na <- pa_ordered_incomplete()
  pw <- parallel_analysis(ord_na, ordered = TRUE, missing = "pairwise",
                          B = 5, seed = 113)
  expect_identical(attr(pw, "method")$missing, "pairwise")
  # Pairwise: the polychoric correlations of every available pair.
  R <- as.matrix(lavaan::lavCor(ord_na, ordered = names(ord_na),
                                missing = "pairwise", output = "cor"))
  f1 <- stats::factanal(covmat = list(cov = R, n.obs = 361), factors = 1)
  diag(R) <- 1 - f1$uniquenesses
  expect_equal(pw$estimate, eigen(R, symmetric = TRUE, only.values = TRUE)$values,
               tolerance = 1e-6)
  lw <- parallel_analysis(ord_na, ordered = TRUE, missing = "listwise",
                          B = 5, seed = 113)
  expect_identical(attr(lw, "method")$missing, "listwise")
  expect_identical(attr(lw, "method")$N, 361L)
  # Continuous items keep their default, FIML, whether or not it was
  # given.
  cont <- holzinger_swineford[c("t1_visual_perception", "t2_cubes",
                                "t4_lozenges", "t6_paragraph_comprehension")]
  expect_identical(
    attr(parallel_analysis(cont, B = 5, seed = 113), "method")$missing, "fiml")
})
