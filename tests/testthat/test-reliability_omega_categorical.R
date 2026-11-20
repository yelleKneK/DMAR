# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# Tests for reliability_omega_categorical(). Point estimate is checked for being in
# (0, 1) and close to the population value used in the simulation. CI
# methods are the bootstrap families. Computations involve a categorical
# CFA via lavaan which is heavy; bootstrap reps are kept small here.

skip_if_not_installed("lavaan")

set.seed(113)
N <- 400
J <- 5
loadings <- rep(0.7, J)
eta <- rnorm(N)
latent <- outer(eta, loadings) +
          matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
items_categorical <- apply(latent, 2, function(x)
  as.integer(cut(x, breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf),
                 labels = FALSE)))
colnames(items_categorical) <- paste0("y", seq_len(J))


test_that("reliability_omega_categorical() matches the MBESS reference implementation", {
  # One ordered-categorical CFA via lavaan with WLSMV, about a quarter of a
  # second, and the same fit runs unguarded further down this file, so there
  # is nothing to save by skipping it on CRAN.
  #
  # The assertion is an equality against MBESS, not a window. A window wide
  # enough to survive sampling variability (the old one spanned 0.6 to 0.95)
  # also admits at least three different coefficients in place of the
  # Green and Yang categorical omega: coefficient alpha on the raw integer
  # scores, linear omega on the raw scores, and alpha on the polychoric
  # matrix. Anchoring to the reference implementation is what makes this a
  # test of the coefficient rather than of its order of magnitude.
  r <- reliability_omega_categorical(data = items_categorical, ci_method = "none")
  est <- r$value[r$term == "estimate"]
  expect_gt(est, 0)
  expect_lt(est, 1)

  # Pinned from MBESS::ci.reliability (MBESS 4.9.3, 2026-08-09); live
  # comparison in tools/oracle_checks.R.
  expect_equal(est, 0.8106007039872396, tolerance = 1e-6)
})

test_that("reliability_omega_categorical() returns the expected tidy data.frame and attrs", {
  r <- reliability_omega_categorical(data = items_categorical, ci_method = "none")
  expect_s3_class(r, "data.frame")
  expect_named(r, c("term", "value"))
  expect_true(all(c("estimate", "se", "ci_lower", "ci_upper",
                    "conf_level", "N", "n_items") %in% r$term))
  expect_equal(attr(r, "coefficient"), "omega_categorical")
  expect_equal(attr(r, "ci_method"), "none")
  # The interval menu is bootstrap only, and a bootstrap standard error
  # is already on the coefficient scale, so no se_transformed row (or
  # se_transform_scale attribute) ever appears here.
  expect_false("se_transformed" %in% r$term)
  expect_null(attr(r, "se_transform_scale"))
})

test_that("reliability_omega_categorical() recommends bca but never bootstraps by default", {
  # "bca" leads the menu as the recommended method (KP 2016), but the
  # default call reports the point estimate with a message; no bootstrap
  # runs unless requested.
  expect_equal(eval(formals(reliability_omega_categorical)$ci_method)[1], "bca")
  expect_message(
    r <- reliability_omega_categorical(data = items_categorical),
    regexp = "bootstrap based"
  )
  expect_equal(attr(r, "ci_method"), "none")
  expect_true(is.na(r$value[r$term == "ci_lower"]))
})

test_that("reliability_omega_categorical() BCa bootstrap CI brackets the estimate", {
  # Skipped on CRAN: B = 60 ordered-categorical CFA refits via lavaan
  # WLSMV takes roughly 30-40 s on CRAN's check machines.
  skip_on_cran()
  skip_if_not_installed("boot")
  r <- reliability_omega_categorical(data = items_categorical, ci_method = "bca",
                           B = 60)
  est <- r$value[r$term == "estimate"]
  ll  <- r$value[r$term == "ci_lower"]
  ul  <- r$value[r$term == "ci_upper"]
  expect_lte(ll, est)
  expect_gte(ul, est)
  expect_gte(ll, 0)
  expect_lte(ul, 1)
  expect_equal(attr(r, "B"), 60)
})

test_that("reliability_omega_categorical() percentile bootstrap CI brackets the estimate", {
  # Same reason as the BCa block above: heavy lavaan WLSMV inside the
  # bootstrap loop. Skipped on CRAN.
  skip_on_cran()
  skip_if_not_installed("boot")
  r <- reliability_omega_categorical(data = items_categorical, ci_method = "percentile",
                           B = 60)
  est <- r$value[r$term == "estimate"]
  ll  <- r$value[r$term == "ci_lower"]
  ul  <- r$value[r$term == "ci_upper"]
  expect_lte(ll, est)
  expect_gte(ul, est)
  expect_gte(ll, 0)
  expect_lte(ul, 1)
})

test_that("reliability_omega_categorical() rejects non-integer (continuous) data", {
  expect_error(
    reliability_omega_categorical(data = items_categorical + 0.1),
    regexp = "integer category"
  )
})

test_that("reliability_omega_categorical() errors when data is missing", {
  expect_error(reliability_omega_categorical(), regexp = "'data'.*is required")
})

test_that("reliability_omega_categorical() rejects too few items", {
  expect_error(
    reliability_omega_categorical(data = items_categorical[, 1, drop = FALSE]),
    regexp = "three or more items"
  )
})

test_that("two ordered items are refused, since the single-factor model is not identified", {
  # With two items the one polychoric correlation determines only the
  # product of the two loadings, while categorical omega depends on each,
  # so the value a fit returned depended on the optimizer's start values
  # and on the order of the columns (0.748 one way and 0.758 the other
  # on one data set in the 1.1.0 review). An error now replaces that
  # number, under either denominator, from reliability() as well (whose
  # automatic routing sends integer items here and lets the error
  # through rather than substituting another coefficient), and from the
  # internal helper, so no caller can fit a two-item model.
  two <- items_categorical[, 1:2]
  for (den in c("observed", "model_implied")) {
    expect_error(
      reliability_omega_categorical(data = two, ci_method = "none",
                                    denominator = den),
      "Categorical omega needs three or more items: the single-factor model for ordered items is not identified with two",
      fixed = TRUE)
  }
  expect_error(suppressMessages(reliability(data = two)),
               "Categorical omega needs three or more items", fixed = TRUE)
  expect_error(.omega_c_parts(two), "three or more items", fixed = TRUE)
})

test_that("reliability_omega_categorical() rejects out-of-range conf_level", {
  expect_error(
    reliability_omega_categorical(data = items_categorical, conf_level = -0.1),
    regexp = "conf_level"
  )
})


test_that("a leave-one-out refit that fails no longer stops the BCa interval", {
  # The 2026-10-14 report: with B below n, boot estimates the BCa
  # acceleration by the jackknife, and one leave-one-out refit that
  # returned NA made the acceleration NA and stopped the call with
  # "estimated adjustment 'a' is NA". Here the first of six binary items
  # (n = 60) is endorsed by a single respondent, so the refit without
  # that case refuses a constant item: exactly one jackknife replicate
  # fails, and the resamples that omit the case fail the same way. The
  # acceleration is now estimated from the refits that returned a value,
  # the counts are on the result, and the interval is BCa. Skipped on
  # CRAN: n + B ordered-categorical refits.
  skip_on_cran()
  skip_if_not_installed("boot")
  set.seed(113)
  n_bad <- 60
  j_bad <- 6
  eta_bad <- rnorm(n_bad)
  items_bad <- (outer(eta_bad, rep(0.7, j_bad)) +
                  matrix(rnorm(n_bad * j_bad), n_bad, j_bad) * sqrt(1 - 0.7^2)) > 0
  items_bad <- apply(items_bad, 2, as.integer)
  colnames(items_bad) <- paste0("y", seq_len(j_bad))
  items_bad[, 1] <- 0L
  items_bad[1, 1] <- 1L
  # A sparse resample can still give a lavaan warning. The constant
  # item is refused before fitting, so no variable table is printed.
  quiet <- function(expr) {
    warnings <- testthat::capture_warnings(
      output <- utils::capture.output(withCallingHandlers(
        expr,
        warning = function(w) {
          if (grepl("lavaan", conditionMessage(w))) invokeRestart("muffleWarning")
        })))
    expect_length(output, 0L)
    warnings
  }
  # Two thirds of the resamples omit case 1 or hit a degenerate table and
  # are set aside, so categorical omega, which forms no interval from
  # fewer than half of the replications, reports none here. The engine's
  # jackknife is checked through the engine itself, with a point function
  # that fails as the refit does, when the single respondent who endorses
  # the first item is left out, and otherwise returns categorical omega of
  # the other five items, so that at least half of the replications
  # return a value and the BCa interval is formed.
  w <- quiet(r <- reliability_omega_categorical(
    data = items_bad, ci_method = "bca", B = 30, conf_level = 0.80,
    seed = 113))
  expect_lt(attr(r, "bootstrap")$returned, 15L)
  expect_true(is.na(r$value[r$term == "ci_lower"]))
  expect_match(w, "fewer than half returned a value", fixed = TRUE,
               all = FALSE)
  omega_five <- .omega_c_point_fn("observed")
  point_fn <- function(d) {
    if (!any(d[, 1L] == 1L)) stop("y1 takes a single value")
    omega_five(d[, -1L, drop = FALSE])
  }
  # Seed 1 and not 113: of the 30 resamples of seed 113, 16 leave that
  # respondent out, more than half, and no interval is formed from those.
  w <- quiet(e <- .bootstrap_ci(items_bad, point_fn, B = 30,
                                conf_level = 0.80, kind = "bca", seed = 1,
                                improper = "drop"))
  # boot adds its own warning that extreme order statistics served as
  # endpoints; the package warning is the one that names the jackknife
  # count, and it is the only other one.
  own <- grepl("1 of 60 jackknife replicates behind the BCa acceleration failed", w)
  expect_identical(sum(own), 1L)
  expect_match(w[own], "bootstrap replications were set aside")
  expect_true(all(grepl("extreme order statistics used as endpoints", w[!own])))
  expect_identical(e$ci_method, "bca")
  bs <- e$bootstrap
  expect_identical(bs$attempted, 30L)
  expect_identical(bs$returned + bs$failed, 30L)
  expect_gte(bs$returned, 15L)
  expect_lt(bs$returned, 30L)
  expect_identical(bs$jackknife$attempted, 60L)
  expect_identical(bs$jackknife$returned, 59L)
  expect_identical(bs$jackknife$failed, 1L)
  expect_identical(bs$acceleration_method, "jackknife")
  expect_true(is.finite(bs$acceleration))
  est <- as.numeric(omega_five(items_bad[, -1L]))
  expect_lte(e$lower, est)
  expect_gte(e$upper, est)
  expect_gte(e$lower, 0)
  expect_lte(e$upper, 1)
})


# The denominator argument -------------------------------------------------
#
# Green and Yang (2009) built the denominator of categorical omega from
# the latent-response correlations the fitted model reproduces (their
# appendix program takes "the reproduced (model-based) correlation
# matrix" as its input), and Kelley and Pornprasertmanit (2016, Eq. 18)
# from the unrestricted polychoric correlations. The default stays the
# unrestricted one, which the MBESS anchor above pins.

test_that("the categorical omega assembly reproduces Green and Yang's (2009) population values", {
  # Population values from Green and Yang (2009, Table 1, p. 162, the
  # nonlinear SEM column): four items, loadings of .8 (or .8, .8, .5,
  # .5), and the category probabilities of their design (p. 161). In the
  # population the model implied and the unrestricted polychoric
  # correlations coincide, so these values check the Eq. 21 assembly
  # that both denominators share. Nine rows agree to the third decimal;
  # the binary row with two negatively and two positively skewed items
  # computes as 0.59042 against the published .591 (a difference of
  # 0.00058, larger than rounding allows), so that row alone has a
  # tolerance of 0.001.
  thr <- function(P) stats::qnorm(cumsum(P)[-length(P)])
  gy <- function(P_list, lambda) {
    rho_num <- outer(lambda, lambda)
    rho_den <- rho_num
    diag(rho_den) <- 1
    .gy_omega_from_components(lapply(P_list, thr), rho_num, rho_den)
  }
  u2 <- c(.5, .5); neg2 <- c(.10, .90); pos2 <- rev(neg2)
  u5 <- rep(.2, 5); bell <- c(.09, .15, .51, .15, .09)
  matched <- c(.485, .01, .01, .01, .485)
  smooth_neg <- c(.011, .027, .047, .118, .797); smooth_pos <- rev(smooth_neg)
  l8 <- rep(.8, 4); l85 <- c(.8, .8, .5, .5)
  computed <- c(
    gy(rep(list(u2), 4), l8), gy(rep(list(u5), 4), l8),
    gy(rep(list(bell), 4), l8), gy(rep(list(matched), 4), l8),
    gy(rep(list(neg2), 4), l8), gy(list(neg2, neg2, pos2, pos2), l8),
    gy(rep(list(smooth_neg), 4), l8),
    gy(list(smooth_neg, smooth_neg, smooth_pos, smooth_pos), l8),
    gy(rep(list(u2), 4), l85), gy(rep(list(u5), 4), l85))
  published <- c(.760, .850, .840, .768, .688, .591, .786, .716, .617, .720)
  expect_true(all(abs(computed[-6] - published[-6]) < 0.0005))
  expect_lt(abs(computed[6] - published[6]), 0.001)
  expect_equal(computed[6], 0.59042, tolerance = 1e-5)
})

test_that("reliability_omega_categorical() records its denominator, observed by default", {
  r_def <- reliability_omega_categorical(data = items_categorical, ci_method = "none")
  r_obs <- reliability_omega_categorical(data = items_categorical, ci_method = "none",
                                         denominator = "observed")
  expect_equal(attr(r_def, "denominator"), "observed")
  expect_equal(r_def$value, r_obs$value)
  expect_equal(eval(formals(reliability_omega_categorical)$denominator)[1],
               "observed")
  expect_error(reliability_omega_categorical(data = items_categorical,
                                             ci_method = "none",
                                             denominator = "saturated"))
})

test_that("denominator = 'model_implied' is Green and Yang's Eq. 21 with the reproduced correlations", {
  # An assembly written apart from the package's own: the WLSMV fit's
  # loadings, thresholds, and implied correlations go into Eq. 21, with
  # each bivariate normal probability taken by one-dimensional numerical
  # integration of its defining integral rather than by the package's
  # .pbinorm(). (Through version 1.0.1 this test called mvtnorm here; the
  # same comparison through mvtnorm is in tools/oracle_snippets/mvtnorm.R.)
  r <- reliability_omega_categorical(data = items_categorical, ci_method = "none",
                                     denominator = "model_implied")
  est <- r$value[r$term == "estimate"]
  expect_equal(attr(r, "denominator"), "model_implied")

  d <- as.data.frame(items_categorical)
  d[] <- lapply(d, ordered)
  fit <- lavaan::cfa(paste("f =~", paste(colnames(d), collapse = " + ")),
                     data = d, ordered = colnames(d), std.lv = TRUE,
                     se = "none")
  lam <- lavaan::lavInspect(fit, "est")$lambda[, 1]
  implied <- stats::cov2cor(lavaan::lavInspect(fit, "implied")$cov)
  tau <- lavaan::lavInspect(fit, "est")$tau
  th <- split(as.numeric(tau), sub("\\|.*$", "", rownames(tau)))[colnames(d)]
  # P(Z1 <= a, Z2 <= b) = int_{-Inf}^{a} phi(x) Phi((b - r x) / s) dx,
  # with s = sqrt(1 - r^2), and Phi(min(a, b)) at r = 1 (the diagonal of
  # the denominator, which cov2cor() sets to one exactly).
  p2_cell <- function(a, b, r) {
    if (r == 1) return(stats::pnorm(min(a, b)))
    s <- sqrt(1 - r^2)
    stats::integrate(function(x) {
      stats::dnorm(x) * stats::pnorm((b - r * x) / s)
    }, -Inf, a, rel.tol = 1e-10)$value
  }
  cov_cat <- function(t1, t2, r) {
    p2 <- 0
    for (a in t1) for (b in t2) p2 <- p2 + p2_cell(a, b, r)
    p2 - sum(stats::pnorm(t1)) * sum(stats::pnorm(t2))
  }
  num <- den <- 0
  for (j in seq_along(lam)) for (k in seq_along(lam)) {
    num <- num + cov_cat(th[[j]], th[[k]], lam[j] * lam[k])
    den <- den + cov_cat(th[[j]], th[[k]], implied[j, k])
  }
  expect_equal(est, num / den, tolerance = 1e-6)
})

test_that("the two denominators nearly agree when one factor holds and part ways under a doublet", {
  r_obs <- reliability_omega_categorical(data = items_categorical, ci_method = "none")
  r_mod <- reliability_omega_categorical(data = items_categorical, ci_method = "none",
                                         denominator = "model_implied")
  expect_lt(abs(r_obs$value[1] - r_mod$value[1]), 0.01)

  # A shared specific factor on the first two items: the unrestricted
  # polychorics carry the doublet into the total variance while the
  # reproduced correlations spread it through the inflated loadings, so
  # the two coefficients part ways (the direction depends on the
  # configuration, so only the size of the gap is asserted), and the
  # residual correlations point at the pair.
  set.seed(113)
  n <- 600
  eta <- rnorm(n); s <- rnorm(n)
  lat <- outer(eta, rep(0.6, 6)) + matrix(rnorm(n * 6, sd = 0.7), n, 6)
  lat[, 1:2] <- lat[, 1:2] + 0.6 * s
  it <- apply(lat, 2, function(x)
    as.integer(cut(x, c(-Inf, -1, 0, 1, Inf))))
  colnames(it) <- paste0("q", 1:6)
  d_obs <- reliability_omega_categorical(data = it, ci_method = "none")
  d_mod <- reliability_omega_categorical(data = it, ci_method = "none",
                                         denominator = "model_implied")
  expect_gt(abs(d_mod$value[1] - d_obs$value[1]), 0.005)
  rc <- attr(d_obs, "residual_cor")
  expect_setequal(rc$items, c("q1", "q2"))
  expect_gt(rc$largest, 0)
})

test_that("reliability_omega_categorical() carries the residual polychoric correlations", {
  r <- reliability_omega_categorical(data = items_categorical, ci_method = "none")
  rc <- attr(r, "residual_cor")
  expect_type(rc, "list")
  expect_equal(dim(rc$matrix), c(J, J))
  expect_equal(rownames(rc$matrix), colnames(items_categorical))
  expect_equal(rc$matrix, t(rc$matrix))
  expect_equal(unname(diag(rc$matrix)), rep(0, J), tolerance = 1e-8)
  up <- rc$matrix[upper.tri(rc$matrix)]
  expect_equal(abs(rc$largest), max(abs(up)))
  expect_true(all(rc$items %in% colnames(items_categorical)))
  expect_output(print(r), "Denominator: variance of the summed categories")
  expect_output(print(r), "Largest residual correlation")
})

test_that("a bootstrap interval recomputes the coefficient with the chosen denominator", {
  skip_on_cran()
  r <- reliability_omega_categorical(data = items_categorical,
                                     ci_method = "percentile", B = 60,
                                     seed = 113,
                                     denominator = "model_implied")
  est <- r$value[r$term == "estimate"]
  expect_equal(attr(r, "denominator"), "model_implied")
  expect_lte(r$value[r$term == "ci_lower"], est)
  expect_gte(r$value[r$term == "ci_upper"], est)
})


# An undetermined polychoric correlation ------------------------------------
#
# In the mixed-skew design of Green and Yang (2009), two items piled into
# the top category and two into the bottom, samples of 500 often have a
# pair whose bivariate table leaves the polychoric correlation without a
# finite standard error. Before DMAR 1.1.0 the saturated refit returned
# its starting value, 0, for that correlation and the observed
# denominator silently dropped the pair.

mixed_skew_items <- function(N, seed) {
  set.seed(seed)
  p_top <- c(.011, .027, .047, .118, .797)
  cut_top <- c(-Inf, stats::qnorm(cumsum(p_top)[-5]), Inf)
  cut_bottom <- c(-Inf, stats::qnorm(cumsum(rev(p_top))[-5]), Inf)
  latent <- outer(stats::rnorm(N), rep(0.8, 4)) +
    matrix(stats::rnorm(N * 4, sd = 0.6), N, 4)
  cbind(y1 = cut(latent[, 1], cut_top, labels = FALSE),
        y2 = cut(latent[, 2], cut_top, labels = FALSE),
        y3 = cut(latent[, 3], cut_bottom, labels = FALSE),
        y4 = cut(latent[, 4], cut_bottom, labels = FALSE))
}

test_that("the unrestricted polychorics come from the single-factor fit", {
  # Equal to lavaan's pairwise polychorics, so equal to what a converged
  # saturated fit reproduces.
  parts <- .omega_c_parts(items_categorical)
  d <- as.data.frame(lapply(as.data.frame(items_categorical), ordered))
  ref <- unclass(lavaan::lavCor(d, ordered = names(d)))
  expect_equal(unname(parts$rho_sat), unname(ref), tolerance = 1e-8)
  expect_equal(nrow(parts$undetermined), 0L)
})

test_that("an undetermined polychoric stops the observed denominator and not the model implied one", {
  it <- mixed_skew_items(500, 113)
  # The pair: every respondent below the top category of y1 is in the
  # bottom category of y3.
  expect_true(all(it[it[, "y1"] < 5, "y3"] == 1))
  parts <- .omega_c_parts(it)
  expect_equal(sort(as.vector(parts$undetermined)), c("y1", "y3"))
  expect_error(
    suppressMessages(reliability_omega_categorical(data = it,
                                                   ci_method = "none")),
    "do not determine the polychoric correlation of y1 and y3")
  expect_true(is.na(.omega_c_from_data(it, denominator = "observed")))
  expect_message(
    r <- reliability_omega_categorical(data = it, ci_method = "none",
                                       denominator = "model_implied"),
    "do not determine the polychoric correlation of y1 and y3")
  est <- r$value[r$term == "estimate"]
  expect_gt(est, 0.68)
  expect_lt(est, 0.76)
  rc <- attr(r, "residual_cor")
  expect_true(is.na(rc$matrix["y1", "y3"]))
  expect_false(setequal(rc$items, c("y1", "y3")))
})

test_that("the help page's items as factors labeled 0 to 4 give the result of their codes", {
  # The NEWS example of the change: version 1.0.1 refused these factors
  # ("column 1 contains non-integer values"), admitting a factor only
  # when its labels were 1, 2, ..., K.
  it <- mixed_skew_items(1000, 113)
  fac <- as.data.frame(lapply(as.data.frame(it - 1L), factor))
  expect_identical(levels(fac$y1), as.character(0:4))
  r <- reliability_omega_categorical(data = fac, ci_method = "none")
  expect_equal(r$value[1], 0.7201391554, tolerance = 1e-9)
  expect_identical(r$value,
                   reliability_omega_categorical(data = it,
                                                 ci_method = "none")$value)
})

test_that("the undetermined messages agree in number with the pairs they name", {
  # The seed 17 sample at N = 300 leaves two pairs undetermined while its
  # loadings are identified; the messages spoke of one table, one
  # standard error, and "its residual correlation".
  it <- mixed_skew_items(300, 17)
  expect_identical(nrow(.omega_c_parts(it)$undetermined), 2L)
  err <- tryCatch(suppressMessages(reliability_omega_categorical(
    data = it, ci_method = "none")), error = function(e) conditionMessage(e))
  expect_match(err, "the bivariate tables have empty cells", fixed = TRUE)
  expect_match(err, "so the standard errors are unbounded", fixed = TRUE)
  expect_match(err, "the correlations the single-factor model implies for those pairs",
               fixed = TRUE)
  msgs <- testthat::capture_messages(reliability_omega_categorical(
    data = it, ci_method = "none", denominator = "model_implied"))
  expect_match(msgs, "does not need them", fixed = TRUE, all = FALSE)
  expect_match(msgs, "their residual correlations are reported as NA",
               fixed = TRUE, all = FALSE)
  # One pair keeps the singular.
  it1 <- mixed_skew_items(500, 113)
  msgs <- testthat::capture_messages(reliability_omega_categorical(
    data = it1, ci_method = "none", denominator = "model_implied"))
  expect_match(msgs, "(the bivariate table has empty cells). The model implied denominator does not need it,",
               fixed = TRUE, all = FALSE)
  expect_match(msgs, "its residual correlation is reported as NA",
               fixed = TRUE, all = FALSE)
})

test_that("the determined correlations identify the loadings only when they connect every item and hold an odd cycle", {
  items <- paste0("y", 1:4)
  und <- function(...) matrix(as.character(c(...)), ncol = 2L, byrow = TRUE)
  # Every correlation determined: the complete graph, with triangles.
  expect_true(.omega_c_identification(items, und())$identified)
  # Two undetermined pairs sharing y3 leave the triangle y1, y2, y4, with
  # y3 joined to it by its one determined correlation, with y4. y3 lies
  # on no cycle of an odd number of items and is identified all the same,
  # which is why the condition is a connected graph holding at least one
  # such cycle, not one through every item; the notes say so.
  expect_true(.omega_c_identification(items, und("y1", "y3", "y2", "y3"))$identified)
  note <- .omega_c_undetermined_note("model_implied")
  expect_match(note, "connect every item to the others, directly or through other items, and include at least one cycle of an odd number",
               fixed = TRUE)
  expect_false(grepl("through a cycle", note, fixed = TRUE))
  # Two disjoint undetermined pairs leave the 4-cycle y1-y2-y3-y4-y1,
  # which is bipartite: the sides' loadings trade a common factor.
  g <- .omega_c_identification(items, und("y1", "y3", "y2", "y4"))
  expect_false(g$identified)
  expect_length(g$blocks, 1L)
  expect_identical(g$sides[[1L]], list(c("y1", "y3"), c("y2", "y4")))
  # Every cross pair undetermined: two blocks, each a single edge.
  cross <- und("y1", "y3", "y2", "y3", "y1", "y4", "y2", "y4")
  g <- .omega_c_identification(items, cross)
  expect_false(g$identified)
  expect_identical(g$blocks, list(c("y1", "y2"), c("y3", "y4")))
  msg <- .omega_c_unidentified_message(g, cross)
  expect_match(msg, "links {y1, y2} to {y3, y4}", fixed = TRUE)
  expect_match(msg, "relative sign", fixed = TRUE)
  expect_match(msg, "under either denominator", fixed = TRUE)
  # An item with no determined correlation is a block of its own.
  lone <- und("y1", "y2", "y1", "y3", "y1", "y4")
  g <- .omega_c_identification(items, lone)
  expect_false(g$identified)
  expect_match(.omega_c_unidentified_message(g, lone),
               "the loading of y1, which shares a determined correlation with no item",
               fixed = TRUE)
  # Two blocks that each hold a triangle: the scales are fixed, the
  # relative sign is not.
  six <- paste0("y", 1:6)
  between <- as.matrix(expand.grid(c("y1", "y2", "y3"), c("y4", "y5", "y6"),
                                   stringsAsFactors = FALSE))
  g <- .omega_c_identification(six, between)
  expect_false(g$identified)
  expect_null(g$sides[[1L]])
  expect_false(grepl("scale", .omega_c_unidentified_message(g, between)))
})

test_that("a sample whose determined correlations do not identify the loadings is refused under either denominator", {
  # Green and Yang's mixed-skew design at N = 300 (seed 10): every
  # correlation of y1 or y2 with y3 or y4 is undetermined, so the
  # determined pairs y1-y2 and y3-y4 fix only two products of four
  # loadings. lavaan's default start values gave loadings of opposite
  # signs, and the model implied denominator returned categorical omega
  # 0.0485, against a population value of 0.716, with a message that it
  # did not need those correlations.
  it <- mixed_skew_items(300, 10)
  parts <- .omega_c_parts(it)
  expect_identical(nrow(parts$undetermined), 4L)
  expect_false(parts$identified)
  for (den in c("model_implied", "observed")) {
    err <- tryCatch(reliability_omega_categorical(data = it, ci_method = "none",
                                                  denominator = den),
                    error = function(e) conditionMessage(e))
    expect_type(err, "character")
    expect_match(err, "The single-factor loadings are not identified",
                 fixed = TRUE)
    expect_match(err, "links {y1, y2} to {y3, y4}", fixed = TRUE)
    # No advice to use the other denominator, which fails the same way.
    expect_false(grepl('Use denominator = "model_implied"', err, fixed = TRUE))
  }
  # Every replication of a bootstrap runs through this point function,
  # which now sets such a resample aside, counted as undetermined.
  for (den in c("model_implied", "observed")) {
    v <- .omega_c_from_data(it, denominator = den)
    expect_true(is.na(v))
    expect_identical(attr(v, "dmar_reason"), "undetermined_polychoric")
  }
  # A sample whose determined correlations do identify the loadings (one
  # undetermined pair) still gets the model implied value, and the
  # observed denominator's error offers it.
  it <- mixed_skew_items(500, 113)
  expect_true(.omega_c_parts(it)$identified)
  expect_error(suppressMessages(reliability_omega_categorical(data = it,
                                                              ci_method = "none")),
               "Use denominator = \"model_implied\".*identify the loadings here")
  expect_message(reliability_omega_categorical(data = it, ci_method = "none",
                                               denominator = "model_implied"),
                 "the correlations the data do determine identify")
})

test_that("on the help page's incomplete example, listwise deletion leaves the loadings unidentified", {
  # The example of ?reliability_omega_categorical. Its complete rows keep
  # none of the few respondents who determine the correlations of y2
  # with y3 and of y1 with y4, so the determined pairs form the 4-cycle
  # y1-y2-y4-y3-y1, and loadings that trade a factor between {y1, y4}
  # and {y2, y3} fit equally well: the model implied value was 0.6016 at
  # lavaan's start values and 0.598 to 0.651 at others. The example
  # called it.
  it <- mixed_skew_items(1000, 113)
  z <- as.numeric(scale(it[, "y1"] + it[, "y2"]))
  it[stats::runif(1000) < stats::plogis(-1.5 - 1.2 * z), "y3"] <- NA
  it[stats::runif(1000) < stats::plogis(-1.5 - 1.2 * z), "y4"] <- NA
  expect_error(reliability_omega_categorical(data = it, ci_method = "none",
                                             denominator = "model_implied"),
               "joins an item of {y1, y4} to an item of {y2, y3}",
               fixed = TRUE)
  # Pairwise deletion keeps them, and only y2 with y3 is undetermined.
  parts <- .omega_c_parts(it, missing = "pairwise")
  expect_identical(nrow(parts$undetermined), 1L)
  expect_true(parts$identified)
})

test_that("no interval is reported when fewer than half of the replications return a value", {
  # The seed 17 sample of the mixed-skew design at N = 300: its determined
  # correlations identify the loadings and its solution is proper, but
  # most resamples omit the few respondents who determine two of its
  # cross correlations. Under the model implied denominator 28 of 40
  # replications were unidentified and 7 improper, and the interval
  # reported through the first versions of 1.1.0 came from the other 5.
  skip_on_cran()
  skip_if_not_installed("boot")
  it <- mixed_skew_items(300, 17)
  expect_true(.omega_c_parts(it)$identified)
  w <- testthat::capture_warnings(r <- suppressMessages(
    reliability_omega_categorical(data = it, denominator = "model_implied",
                                  ci_method = "percentile", B = 40,
                                  seed = 1)))
  bs <- attr(r, "bootstrap")
  expect_identical(c(bs$returned, bs$undetermined, bs$improper),
                   c(5L, 28L, 7L))
  expect_true(bs$require_half)
  val <- function(t) r$value[r$term == t]
  expect_true(is.finite(val("estimate")))
  expect_true(is.na(val("se")) && is.na(val("ci_lower")) &&
                is.na(val("ci_upper")))
  expect_match(w, "fewer than half returned a value, so no interval is reported",
               fixed = TRUE, all = FALSE)
  expect_output(print(r), "fewer than half returned a value", fixed = TRUE)
})

test_that("a bootstrap counts the replications an undetermined polychoric sets aside, and says what they do to the interval", {
  # The 1.1.0 review: every polychoric of this sample is determined, but
  # the few respondents in the sparse cells of an oppositely skewed pair
  # are often left out of a resample, so the observed denominator sets
  # many replications aside. They were counted only as "returned a
  # non-finite value", the cause of a non-converged fit. Skipped on CRAN:
  # two bootstraps of ordered-categorical refits.
  skip_on_cran()
  skip_if_not_installed("boot")
  it <- mixed_skew_items(300, 3)
  parts <- .omega_c_parts(it)
  expect_equal(nrow(parts$undetermined), 0L)
  expect_length(parts$improper, 0L)
  w <- testthat::capture_warnings(
    r <- reliability_omega_categorical(data = it, ci_method = "percentile",
                                       B = 20, seed = 113))
  bs <- attr(r, "bootstrap")
  # The reasons keep the shape they have had since 1.0.1; the causes are
  # counted within the non-finite replications.
  expect_named(bs$reasons, c("error", "non_finite"))
  expect_gt(bs$undetermined, 0L)
  expect_lte(bs$undetermined + bs$improper, bs$reasons[["non_finite"]])
  expect_identical(bs$returned + bs$failed, 20L)
  own <- grepl("bootstrap replications were set aside", w)
  expect_identical(sum(own), 1L)
  expect_match(w[own], sprintf("%d left a polychoric correlation undetermined",
                               bs$undetermined), fixed = TRUE)
  expect_match(w[own], "redrew the respondents who determine it", fixed = TRUE)
  expect_match(w[own], 'denominator = "model_implied"', fixed = TRUE)
  expect_output(print(r), sprintf("%d left a polychoric correlation undetermined",
                                  bs$undetermined), fixed = TRUE)
  expect_false(any(grepl(sprintf("%d returned a non-finite value",
                                 bs$undetermined), w, fixed = TRUE)))

  # The model implied denominator needs no unrestricted polychoric and
  # sets aside only the replications whose determined correlations leave
  # the loadings unidentified (5 of these 20), a subset of those the
  # observed denominator sets aside on the same resamples. It kept every
  # one of them through the first versions of 1.1.0.
  w_mi <- testthat::capture_warnings(
    r_mi <- reliability_omega_categorical(data = it, ci_method = "percentile",
                                          B = 20, seed = 113,
                                          denominator = "model_implied"))
  bs_mi <- attr(r_mi, "bootstrap")
  expect_gt(bs_mi$undetermined, 0L)
  expect_lt(bs_mi$undetermined, bs$undetermined)
  expect_match(w_mi, "identify them only when they connect every item",
               fixed = TRUE, all = FALSE)
  expect_false(any(grepl("keeps every", c(w, w_mi), fixed = TRUE)))
})


# An improper single-factor solution ----------------------------------------
#
# A standardized loading above one puts a communality above one on the
# diagonal of the numerator's correlations, where the assembly evaluates
# the bivariate normal; that failed with mvtnorm's bare "'corr' is not a
# correlation matrix", and under the model implied denominator only
# after a message advising that denominator.

heywood_items <- function(N, lam, tau, seed) {
  set.seed(seed)
  f <- stats::rnorm(N); J <- length(lam)
  ys <- sapply(seq_len(J), function(j) lam[j] * f + sqrt(1 - lam[j]^2) * stats::rnorm(N))
  out <- sapply(seq_len(J), function(j) findInterval(ys[, j], tau[[j]]) + 1L)
  colnames(out) <- paste0("y", seq_len(J))
  as.data.frame(out)
}

test_that("an improper solution stops with an error naming the item, under either denominator", {
  # Symmetric items with no undetermined pair: the first item's loading
  # runs past one.
  d <- heywood_items(80, c(.95, .95, .6, .6), rep(list(c(-1, 0, 1)), 4), 1)
  parts <- .omega_c_parts(d)
  expect_identical(parts$improper, "y1")
  expect_gt(parts$loadings[["y1"]], 1)
  expect_equal(nrow(parts$undetermined), 0L)
  for (den in c("observed", "model_implied")) {
    expect_error(
      reliability_omega_categorical(data = d, ci_method = "none",
                                    denominator = den),
      "improper \\(a Heywood case\\): the standardized loading of y1 is 1\\.[0-9]{3}")
  }
  expect_true(is.na(.omega_c_from_data(d)))
  expect_identical(attr(.omega_c_from_data(d), "dmar_reason"), "improper")
  # lavaan's own warning about the negative variance does not repeat the
  # error once per fit.
  expect_no_warning(.omega_c_parts(d))
})

test_that("an improper solution is reported before an undetermined pair, with no advice to switch denominators", {
  tau <- list(c(0.8, 1.5, 2.0), c(0.9, 1.6, 2.2),
              c(-2.2, -1.6, -0.9), c(-2.0, -1.5, -0.8))
  d <- heywood_items(200, c(.85, .9, .9, .95), tau, 63)
  parts <- .omega_c_parts(d)
  expect_identical(parts$improper, "y2")
  expect_gt(nrow(parts$undetermined), 0L)
  expect_true(parts$identified)
  expect_no_message(expect_error(
    reliability_omega_categorical(data = d, ci_method = "none",
                                  denominator = "model_implied"),
    "the standardized loading of y2 is .*The data also do not determine"))
  # When the undetermined pairs leave the loadings unidentified, the
  # improper verdict reads arbitrary loadings, and the error says that the
  # loadings are not identified (seed 16: every cross pair undetermined;
  # lavaan's start values gave y4 a loading past one).
  d <- heywood_items(200, c(.85, .9, .9, .95), tau, 16)
  parts <- .omega_c_parts(d)
  expect_identical(parts$improper, "y4")
  expect_false(parts$identified)
  expect_no_message(expect_error(
    reliability_omega_categorical(data = d, ci_method = "none",
                                  denominator = "model_implied"),
    "The single-factor loadings are not identified"))
})

test_that("a bootstrap counts improper replications by cause, and lavaan's warning does not leak", {
  # Skipped on CRAN: a bootstrap of ordered-categorical refits.
  skip_on_cran()
  skip_if_not_installed("boot")
  # A proper sample whose largest communality is about .95, so that some
  # resamples push a loading past one.
  d <- heywood_items(80, c(.95, .95, .6, .6), rep(list(c(-1, 0, 1)), 4), 2)
  expect_length(.omega_c_parts(d)$improper, 0L)
  w <- testthat::capture_warnings(
    r <- reliability_omega_categorical(data = d, ci_method = "percentile",
                                       B = 20, seed = 113))
  bs <- attr(r, "bootstrap")
  expect_gt(bs$improper, 0L)
  expect_identical(bs$reasons[["error"]], 0L)
  # Categorical omega is not defined for an improper solution, so such a
  # replication has no estimate to keep: it is set aside, under the rule
  # the attribute records.
  expect_identical(bs$improper_kept, 0L)
  expect_identical(bs$improper_rule, "drop")
  expect_identical(sum(bs$status == "improper"), bs$improper)
  # The count the help documents as always zero for categorical omega.
  expect_identical(bs$no_solution, 0L)
  expect_true(all(is.na(bs$t[bs$status == "improper"])))
  expect_false(any(grepl("ov variances are negative", w, fixed = TRUE)))
  own <- grepl("bootstrap replications were set aside", w)
  expect_match(w[own], sprintf("%d gave an improper (Heywood) solution",
                               bs$improper), fixed = TRUE)
  expect_match(w[own], "trim the upper tail", fixed = TRUE)
})

test_that("with three items the saturated model leaves no residual to name", {
  # Three items fit three loadings to three polychoric correlations, so
  # the residuals are zero up to rounding and the pair carrying the
  # largest rounding error says nothing about local dependence.
  r3 <- reliability_omega_categorical(data = items_categorical[, 1:3],
                                      ci_method = "none")
  rc <- attr(r3, "residual_cor")
  expect_true(rc$saturated)
  expect_true(is.na(rc$largest))
  expect_null(rc$items)
  expect_lt(max(abs(rc$matrix)), 1e-6)
  expect_output(print(r3), "With three items the single-factor model is saturated (df = 0)",
                fixed = TRUE)
  expect_false(any(grepl("Largest residual", capture.output(print(r3)))))
  # Four items keep the largest residual and its pair.
  r4 <- reliability_omega_categorical(data = items_categorical[, 1:4],
                                      ci_method = "none")
  expect_false(attr(r4, "residual_cor")$saturated)
  expect_true(is.finite(attr(r4, "residual_cor")$largest))
})

test_that("each item's codes are scored as consecutive ranks, and a message says when the coded sum differs", {
  # Four three-category items; y1 and y2 recoded 1, 3, 5 beside y3 and y4
  # coded 1, 2, 3. Categorical omega (Green & Yang, 2009, Eq. 11; Kelley &
  # Pornprasertmanit, 2016, Eq. 16) scores each item 0, 1, 2, so the
  # coefficient is unchanged, and a message names the different steps.
  set.seed(1)
  n <- 400
  lat <- outer(rnorm(n), rep(0.7, 4)) +
    matrix(rnorm(n * 4), n, 4) * sqrt(1 - 0.49)
  d <- as.data.frame(apply(lat, 2, function(x) {
    as.integer(cut(x, c(-Inf, -0.6, 0.6, Inf)))
  }))
  names(d) <- paste0("y", 1:4)
  ref <- expect_silent(reliability_omega_categorical(d, ci_method = "none"))
  gap <- d
  gap$y1 <- c(1L, 3L, 5L)[d$y1]
  gap$y2 <- c(1L, 3L, 5L)[d$y2]
  expect_message(
    r <- reliability_omega_categorical(gap, ci_method = "none"),
    "(y1 and y2 step by 2, y3 and y4 step by 1)", fixed = TRUE)
  expect_identical(r$value, ref$value)
  # The coded sum is a different composite: its reliability, from the same
  # fit with each threshold of y1 and y2 counted twice (a step of 2), is
  # what the message warns of (0.684 here against 0.718).
  p <- .omega_c_parts(gap)
  th <- p$thresh
  th[1:2] <- lapply(th[1:2], rep, each = 2)
  coded <- .gy_omega_from_components(th, p$rho_num, p$rho_sat)
  expect_equal(coded, .gy_omega_from_components(p$thresh, p$rho_num,
                                                p$rho_sat,
                                                weights = c(2, 2, 1, 1)),
               tolerance = 1e-12)
  expect_gt(r$value[r$term == "estimate"] - coded, 0.03)
  # Codes equally spaced with one step on every item are a linear
  # rescaling of the ranks: same coefficient, no message.
  even <- as.data.frame(lapply(d, function(x) 2L * x - 2L))
  r_even <- expect_silent(reliability_omega_categorical(even,
                                                        ci_method = "none"))
  expect_identical(r_even$value, ref$value)
  # Codes that skip a value: bipolar coding, or a category no one chose.
  bip <- d
  bip$y1 <- c(-1L, 1L, 2L)[d$y1]
  expect_message(reliability_omega_categorical(bip, ci_method = "none"),
                 "y1 is coded -1, 1, 2", fixed = TRUE)
})

test_that("listwise deletion that leaves too few rows says so and names the remedy", {
  e <- as.data.frame(items_categorical)
  e[1:200, "y1"] <- NA
  e[201:400, "y2"] <- NA
  # No complete row, and y1 and y2 are never observed together, so only
  # imputation can use these data. lavaan is not reached, so nothing is
  # printed.
  out <- utils::capture.output(
    err <- tryCatch(reliability_omega_categorical(e), error = identity))
  expect_length(out, 0L)
  expect_match(conditionMessage(err), paste0(
    "No row of 'data' is observed on every item, so listwise deletion ",
    "(the default) leaves nothing to analyze. Use missing = \"mi\": no row ",
    "observes both y1 and y2"), fixed = TRUE)
  # With every pair observed together, pairwise deletion is offered too.
  f <- as.data.frame(items_categorical)
  f[1:100, "y1"] <- NA
  f[101:200, "y2"] <- NA
  f[201:300, "y3"] <- NA
  f[301:400, "y4"] <- NA
  expect_error(reliability_omega_categorical(f),
               "Use missing = \"pairwise\", which assumes values missing",
               fixed = TRUE)
  # One complete row is no better.
  f[301, "y4"] <- items_categorical[301, "y4"]
  expect_error(reliability_omega_categorical(f),
               "Only one row of 'data' is observed on every item",
               fixed = TRUE)
  # An item no row observes is named, whatever the treatment.
  g <- as.data.frame(items_categorical)
  g$y5 <- NA_integer_
  expect_error(reliability_omega_categorical(g),
               "No row observes y5, so neither its thresholds", fixed = TRUE)
  # An item that takes one value in the complete rows, though not in
  # all of them.
  h <- as.data.frame(items_categorical)
  h$y3[h$y3 != 3L] <- NA
  h$y3[1:2] <- c(1L, 2L)
  h$y1[1:2] <- NA
  expect_error(reliability_omega_categorical(h),
               paste0("y3 takes a single value in the ",
                      sum(stats::complete.cases(h)), " complete rows, so ",
                      "its thresholds cannot be estimated. Listwise ",
                      "deletion (the default) analyzes only the complete ",
                      "rows"), fixed = TRUE)
  # A constant item in complete data.
  k <- as.data.frame(items_categorical)
  k$y4 <- 2L
  expect_error(reliability_omega_categorical(k),
               paste0("y4 takes a single value in the 400 complete rows, ",
                      "so its thresholds cannot be estimated. Leave it ",
                      "out of 'data'."), fixed = TRUE)
})

test_that("a fit lavaan stops with an error is reported with lavaan's message, not as one that did not converge", {
  # .omega_c_parts() swallowed lavaan's error, so every failure of the
  # categorical fit reached the user as "did not converge". An item
  # with one category (which the function itself refuses before any fit)
  # makes lavaan stop, and the helper now carries its message back.
  one_level <- as.data.frame(items_categorical[1:50, 1:3])
  one_level$y3 <- 2L
  utils::capture.output(parts <- .omega_c_parts(one_level))
  expect_false(parts$converged)
  expect_match(parts$fit_error, "only 1 level", fixed = TRUE)
  expect_true(.omega_c_parts(items_categorical)$converged)
  # A bootstrap replication without a fit has no value, as before.
  utils::capture.output(value <- .omega_c_from_data(one_level))
  expect_true(is.na(value))
  expect_null(attr(value, "dmar_reason"))
  # The point estimate quotes the error; a fit that ran and did not
  # converge keeps the message it had.
  testthat::local_mocked_bindings(.omega_c_parts = function(...) {
    list(converged = FALSE,
         fit_error = "lavaan->lav_x():\n   something went wrong")
  })
  expect_error(
    reliability_omega_categorical(items_categorical, ci_method = "none"),
    paste0("The single-factor categorical CFA model could not be fit; ",
           "lavaan stopped with this error: lavaan->lav_x(): something ",
           "went wrong"), fixed = TRUE)
  testthat::local_mocked_bindings(
    .omega_c_parts = function(...) list(converged = FALSE))
  expect_error(
    reliability_omega_categorical(items_categorical, ci_method = "none"),
    "^The single-factor categorical CFA model did not converge\\.$")
})

test_that("a bootstrap needs B of at least 2, the family's rule", {
  # Through version 1.0.1 B was not checked, and B = 1 returned an
  # interval with NA limits (since 1.1.0 with a warning that one
  # replication cannot make an interval); B = 1 is now refused before
  # anything is fitted, with the message reliability_omega() gives. (The
  # count the rule of ?dmar_bootstrap requires is lowered to 2 for this
  # file; the counts themselves are tested in test-dmar_bootstrap.R.)
  msg <- "'B' must be a single whole number of at least 2."
  for (b in list(1, 0, 2.5, -3, NA, c(10, 20), "100")) {
    expect_error(reliability_omega_categorical(items_categorical,
                                               ci_method = "percentile",
                                               B = b),
                 msg, fixed = TRUE, class = "dmar_count_invalid")
  }
  expect_error(reliability_omega(items_categorical, ci_method = "percentile",
                                 B = 1), msg, fixed = TRUE,
               class = "dmar_count_invalid")
  # Without a bootstrap B is not used, and the earlier check stands.
  expect_no_error(reliability_omega_categorical(items_categorical,
                                                ci_method = "none", B = 1))
  expect_error(reliability_omega_categorical(items_categorical,
                                             ci_method = "none", B = 0),
               "'B' must be NULL or a single whole number of at least 1")
})

test_that("lavaan's error is quoted with the caller's item names", {
  # The categorical fit enters the items into lavaan as y1, ..., yJ, with
  # a1, ..., aJ as the labels of their loadings, so lavaan names an item
  # y1 when it stops; through the first versions of 1.1.0 its error was
  # quoted that way. lavaan 0.7-3 has no error for ordered items that
  # names one and that DMAR does not refuse first (an item with a single
  # category is refused by name before the fit), so lavaan's error is
  # simulated here, naming the tenth item and the first; the item names
  # are mapped back as whole words, and a loading label is left alone.
  set.seed(113)
  ten <- as.data.frame(apply(
    outer(stats::rnorm(300), rep(0.7, 10)) +
      matrix(stats::rnorm(3000, sd = sqrt(1 - 0.49)), 300, 10),
    2, function(x) as.integer(cut(x, c(-Inf, -1, 0, 1, Inf)))))
  names(ten) <- paste0("item_", c(LETTERS[1:9], "J"))
  testthat::local_mocked_bindings(
    cfa = function(...) {
      stop("lavaan->lav_x():\n   ordered variable y10 is wrong, and so is ",
           "y1 (loading a1)", call. = FALSE)
    },
    .package = "lavaan")
  parts <- .omega_c_parts(ten)
  expect_false(parts$converged)
  expect_identical(parts$fit_error,
                   paste0("lavaan->lav_x():\n   ordered variable item_J is ",
                          "wrong, and so is item_A (loading a1)"))
  expect_error(
    reliability_omega_categorical(ten, ci_method = "none"),
    paste0("The single-factor categorical CFA model could not be fit; ",
           "lavaan stopped with this error: lavaan->lav_x(): ordered ",
           "variable item_J is wrong, and so is item_A (loading a1)"),
    fixed = TRUE)
})


test_that("the help describes released history rather than an interim implementation", {
  path <- file.path("..", "..", "man", "reliability_omega_categorical.Rd")
  text <- if (file.exists(path)) {
    paste(readLines(path, warn = FALSE), collapse = " ")
  } else {
    rd <- tools::Rd_db("DMAR")[["reliability_omega_categorical.Rd"]]
    paste(unlist(rd), collapse = " ")
  }
  text <- gsub("[[:space:]]+", " ", text)
  expect_false(grepl("Through the first versions", text, fixed = TRUE))
})


test_that("categorical bootstrap point functions refuse constant items by name before fitting", {
  # A constant resample has no item thresholds. Checking before the fit
  # prevents lavaan's printed variable table and keeps the reason in the
  # engine's error counts, on listwise and pairwise paths alike.
  testthat::local_mocked_bindings(.bootstrap_ci = function(data, point_fn, ...) {
    flat <- data
    flat[, 3] <- 1L
    output <- utils::capture.output(
      err <- tryCatch(point_fn(flat), error = identity))
    expect_length(output, 0L)
    expect_s3_class(err, "error")
    if (inherits(err, "error")) {
      expect_match(conditionMessage(err), "y3 takes a single value", fixed = TRUE)
    }
    list(se = NA_real_, lower = NA_real_, upper = NA_real_)
  })
  for (method in c("listwise", "pairwise")) {
    d <- items_categorical
    if (method == "pairwise") d[1, 1] <- NA
    r <- reliability_omega_categorical(d, ci_method = "percentile", B = 10,
                                        missing = method)
    expect_true(is.finite(r$value[r$term == "estimate"]))
  }
})
