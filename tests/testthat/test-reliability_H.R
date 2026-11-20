# The bootstraps in this file are mechanics tests run with a small B; the
# rule on the number of replications is tested in test-dmar_bootstrap.R.
local_small_B()

# Helpers ----------------------------------------------------------------

# The value of `expr` and the classes of the warnings it signaled, in order.
with_warning_classes <- function(expr) {
  seen <- new.env(parent = emptyenv())
  seen$classes <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    seen$classes <- c(seen$classes, class(w)[1L])
    invokeRestart("muffleWarning")
  })
  list(value = value, classes = seen$classes)
}

verbal <- c("t5_general_information", "t6_paragraph_comprehension",
            "t7_sentence", "t8_word_classification", "t9_word_meaning")
verbal_model <- paste("verbal =~", paste(verbal, collapse = " + "))
three_factor_model <- "
  visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
  textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
  speed   =~ t10_addition + t12_counting_groups_of_dots +
             t13_straight_and_curved_capitals"

# Missing at random on four verbal tests, the chance rising with the visual
# perception score (an observed variable outside the scale): 167 of the 301
# rows stay complete.
mar_fixture <- function() {
  d <- holzinger_swineford[, c(verbal, "t1_visual_perception")]
  withr::with_seed(113, {
    p_miss <- stats::plogis(-2.2 + 0.08 * (d$t1_visual_perception -
                                             mean(d$t1_visual_perception)))
    for (v in verbal[2:5]) d[[v]][stats::runif(nrow(d)) < p_miss * 1.4] <- NA
  })
  d
}

# Three items whose one-factor solution is improper: the correlations imply
# a squared standardized loading of .8 * .6 / .3 = 1.6 for y1.
heywood_S <- function() {
  R <- matrix(c(1, .8, .6,
                .8, 1, .3,
                .6, .3, 1), 3, 3,
              dimnames = list(paste0("y", 1:3), paste0("y", 1:3)))
  R
}

# Agreement with SIEVER ----------------------------------------------------

test_that("the SIEVER complete-case fit reproduces SIEVER's recorded numbers", {
  skip_if_not_installed("lavaan")
  # SIEVER's continuous complete-case fit: ML with robust.sem standard
  # errors and the Satorra-Bentler test. sieve_H_delta() and its recorded
  # h_interval() give H = .918404655996, SE = .00707993242539, and the 95%
  # logit interval [.903409265314, .931249225144]. They read the standard
  # error from lavaan's vcov.std.all, whose Jacobian is numerical; DMAR's
  # gradient agrees with lavaan's own delta method for the H_<factor>
  # defined parameter to about 1e-11, and the two routes differ by about
  # 1e-10 in the standard error and 3e-10 in the limits.
  hs <- holzinger_swineford
  fit <- lavaan::cfa(verbal_model, data = hs[stats::complete.cases(hs[, verbal]), verbal],
                     std.lv = TRUE, estimator = "ML", se = "robust.sem",
                     test = "satorra.bentler", missing = "listwise")
  res <- reliability_H(fit, ci_method = "delta_logistic")
  absolute_gap <- function(term, value) abs(res$value[res$term == term] - value)
  expect_lt(absolute_gap("estimate", 0.918404655996), 1e-10)
  expect_lt(absolute_gap("se", 0.00707993242539), 5e-10)
  expect_lt(absolute_gap("ci_lower", 0.903409265314), 1e-9)
  expect_lt(absolute_gap("ci_upper", 0.931249225144), 1e-9)
  expect_identical(attr(res, "se_type"), "robust.sem")
  expect_equal(reliability_H(fit, ci_method = "mlr_logistic")$value,
               res$value)
})

# On CRAN the file stops here, after the anchor above, which comes first for
# that reason; the blocks below run in the full local suite (see cran_files
# in tests/testthat.R).
skip_on_cran()

# From standardized loadings ---------------------------------------------

test_that("reliability_H() reproduces the Hancock-Mueller (2001) formula", {
  loadings <- c(0.6, 0.7, 0.8)
  theta <- loadings^2 / (1 - loadings^2)
  H_expected <- sum(theta) / (1 + sum(theta))
  res <- reliability_H(loadings)
  expect_equal(res$value[res$term == "estimate"], H_expected,
               tolerance = 1e-12)
})

test_that("reliability_H() H monotone increasing in number of indicators", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  # With identical loadings, adding more indicators must increase H.
  H_two   <- reliability_H(rep(0.7, 2))$value[1]
  H_three <- reliability_H(rep(0.7, 3))$value[1]
  H_six   <- reliability_H(rep(0.7, 6))$value[1]
  expect_gt(H_three, H_two)
  expect_gt(H_six,   H_three)
})

test_that("reliability_H() rejects |loading| >= 1", {
  expect_error(reliability_H(c(0.5, 1.0)), "in \\(-1, 1\\)")
})

test_that("loadings alone return the family table with the estimate and no interval", {
  res <- with_warning_classes(reliability_H(c(0.6, 0.7, 0.8)))
  expect_identical(res$classes, character(0))
  out <- res$value
  expect_s3_class(out, "dmar_reliability")
  expect_identical(out$term, c("estimate", "se", "ci_lower",
                               "ci_upper", "conf_level", "N",
                               "N_complete", "n_items"))
  theta <- c(0.6, 0.7, 0.8)^2 / (1 - c(0.6, 0.7, 0.8)^2)
  expect_equal(out$value[1], sum(theta) / (1 + sum(theta)), tolerance = 1e-14)
  expect_true(all(is.na(out$value[out$term %in% c("se", "ci_lower",
                                                  "ci_upper")])))
  expect_identical(attr(out, "coefficient"), "H")
  expect_identical(attr(out, "ci_method"), "none")
  expect_identical(attr(out, "missing"), "unknown")
  expect_identical(attr(out, "vcov_basis"), "none")
  expect_match(DMAR:::.dmar_tbl_footer(out)[["H"]],
               "No interval .*needs 'vcov_loadings'")
  # The same table by name, by position, and with ci_method = "none".
  expect_identical(reliability_H(loadings = c(0.6, 0.7, 0.8),
                                 ci_method = "none"), out)
})

test_that("se_loadings give the logit interval from the marginal standard errors, with its warning every call", {
  # The stored values are those of reliability_H() as DMAR 1.0.0 computed
  # them (the CRAN commit bcc0624), written in hexadecimal:
  # H 0.76749928768164122, limits 0.71614658505680862 and
  # 0.81200049376189842, variance 0.00059899465096126714, and at
  # conf_level = 0.90 the limits 0.72484946113673421 and
  # 0.80531378590893765. The interval is unchanged; the table is the
  # family table, whose row se is the square root of that variance.
  l <- c(0.6, 0.7, 0.8)
  se <- c(0.05, 0.04, 0.03)
  first <- with_warning_classes(reliability_H(l, se_loadings = se))
  expect_identical(first$classes, "dmar_H_marginal_se")
  again <- with_warning_classes(reliability_H(l, se_loadings = se))
  expect_identical(again$classes, "dmar_H_marginal_se")
  out <- first$value
  val <- function(x, term) x$value[x$term == term]
  expect_s3_class(out, "dmar_reliability")
  expect_identical(attr(out, "ci_method"), "delta_logistic")
  expect_identical(attr(out, "vcov_basis"), "diagonal")
  expect_equal(val(out, "estimate"), 0x1.88f5aaa897861p-1, tolerance = 1e-14)
  expect_equal(val(out, "ci_lower"), 0x1.6eaac3e3ec121p-1,
               tolerance = 1e-14)
  expect_equal(val(out, "ci_upper"), 0x1.9fbe875a161f1p-1,
               tolerance = 1e-14)
  expect_equal(val(out, "se")^2, 0x1.3a0bb37da90aap-11, tolerance = 1e-14)
  expect_false(any(c("reliability_H", "var_H") %in% out$term))

  ninety <- suppressWarnings(reliability_H(l, se, 0.90))
  expect_equal(val(ninety, "ci_lower"), 0x1.731f77f43603ap-1,
               tolerance = 1e-14)
  expect_equal(val(ninety, "ci_upper"), 0x1.9c5216aafe60ep-1,
               tolerance = 1e-14)

  with_N <- suppressWarnings(reliability_H(l, se_loadings = se, N = 200,
                                           ci_method = "delta_logistic"))
  expect_equal(val(with_N, "N"), 200)
  expect_equal(val(with_N, "ci_lower"), val(out, "ci_lower"),
               tolerance = 1e-14)
})

test_that("a bootstrap replication of H in which a variable takes a single value is set aside by name", {
  # Regression test. A resample that left out the one case on which an
  # item differs from the rest reached lavaan, which printed its table of
  # the variables for every such resample; it is now refused before the
  # fit, by name, from data as from a fitted model.
  skip_on_cran()
  skip_if_not_installed("lavaan")
  set.seed(1)
  n <- 60
  f <- stats::rnorm(n)
  d <- as.data.frame(sapply(c(.8, .7, .6, .5), function(l) {
    l * f + stats::rnorm(n, sd = sqrt(1 - l^2))
  }))
  names(d) <- paste0("x", 1:4)
  d$x5 <- 0
  d$x5[which.max(f)] <- 1
  set.seed(3)
  bb <- boot::boot(d, function(x, i) sum(x$x5[i]), R = 20)
  lost <- sum(bb$t[, 1] == 0)
  expect_gt(lost, 0L)
  out <- utils::capture.output(
    r <- suppressWarnings(reliability_H(data = d, ci_method = "percentile",
                                        B = 20, seed = 3)))
  expect_identical(out, character(0))
  bt <- attr(r, "bootstrap")
  expect_identical(bt$reasons[["error"]], lost)
  expect_identical(
    names(bt$error_messages),
    paste0("x5 takes a single value in this resample, so its variance is ",
           "zero, and an item that does not vary is refused."))
  fit <- suppressWarnings(lavaan::cfa("f =~ x1 + x2 + x3 + x4 + x5",
                                      data = d, std.lv = TRUE))
  out <- utils::capture.output(
    r <- suppressWarnings(reliability_H(fit, ci_method = "percentile",
                                        B = 20, seed = 3)))
  expect_identical(out, character(0))
  bt <- attr(r, "bootstrap")
  expect_identical(bt$reasons[["error"]], lost)
  # The fit is listwise (complete data), and its replications are checked
  # as such, the variables named as the model's.
  expect_identical(
    names(bt$error_messages),
    paste0("x5 takes a single value in this resample, so its variance is ",
           "zero, and a variable that does not vary is refused."))
})

test_that("the table says what an interval needs only when nothing came with the loadings", {
  # ci_method = "none" leaves vcov_basis "none" whatever was supplied, so
  # the line was printed, naming vcov_loadings as the remedy, under a
  # call that had supplied it, or se_loadings, and asked for no interval.
  l <- c(0.6, 0.7, 0.8)
  needs <- "needs 'vcov_loadings'"
  foot <- function(x) paste(DMAR:::.dmar_tbl_footer(x), collapse = "\n")
  with_vcov <- reliability_H(loadings = l,
                             vcov_loadings = diag(c(0.05, 0.04, 0.03)^2),
                             ci_method = "none")
  with_se <- reliability_H(loadings = l, se_loadings = c(0.05, 0.04, 0.03),
                           ci_method = "none")
  expect_identical(attr(with_vcov, "vcov_supplied"), "full")
  expect_identical(attr(with_se, "vcov_supplied"), "diagonal")
  expect_identical(attr(with_vcov, "vcov_basis"), "none")
  expect_identical(attr(with_se, "vcov_basis"), "none")
  expect_false(grepl(needs, foot(with_vcov), fixed = TRUE))
  expect_false(grepl(needs, foot(with_se), fixed = TRUE))
  # Loadings alone, by default or by request, still get the line.
  alone <- reliability_H(loadings = l, N = 200)
  asked <- reliability_H(loadings = l, ci_method = "none")
  expect_identical(attr(alone, "vcov_supplied"), "none")
  expect_match(foot(alone), needs, fixed = TRUE)
  expect_match(foot(asked), needs, fixed = TRUE)
  # The interval from what was supplied records it too.
  formed <- reliability_H(loadings = l,
                          vcov_loadings = diag(c(0.05, 0.04, 0.03)^2))
  expect_identical(attr(formed, "vcov_supplied"), "full")
  expect_identical(attr(formed, "vcov_basis"), "full")
  expect_false(grepl(needs, foot(formed), fixed = TRUE))
})

# The fit path -------------------------------------------------------------

test_that("the fit path equals fa_confirmatory()'s H defined parameter (ML, MLR; listwise, FIML)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- mar_fixture()[, verbal]
  for (estimator in c("ML", "MLR")) {
    for (missing in c("listwise", "ml")) {
      fit <- suppressWarnings(fa_confirmatory(
        data = d, factors = list(verbal = verbal), estimator = estimator,
        se = if (estimator == "MLR") "robust.huber.white" else "standard",
        missing = missing, output = "fit"))
      pe <- lavaan::parameterEstimates(fit)
      ref <- pe[pe$lhs == "H_verbal" & pe$op == ":=", ]
      res <- reliability_H(fit, ci_method = "delta")
      expect_equal(res$value[res$term == "estimate"], ref$est,
                   tolerance = 1e-8)
      expect_equal(res$value[res$term == "se"], ref$se, tolerance = 1e-8)
      expect_identical(attr(res, "missing"),
                       if (missing == "ml") "fiml" else "listwise")
      expect_equal(res$value[res$term == "N"],
                   if (missing == "ml") 301 else 167)
      expect_equal(res$value[res$term == "N_complete"], 167)
    }
  }
})

test_that("H is fa_confirmatory()'s H of a factor of a larger model and the squared determinacy of its score", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  factors <- list(
    visual = c("t1_visual_perception", "t2_cubes", "t4_lozenges"),
    textual = c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning"),
    speed = c("t10_addition", "t12_counting_groups_of_dots",
              "t13_straight_and_curved_capitals"))
  fit <- fa_confirmatory(data = holzinger_swineford, factors = factors, output = "fit")
  pe <- lavaan::parameterEstimates(fit)
  for (f in names(factors)) {
    ref <- pe[pe$lhs == paste0("H_", f) & pe$op == ":=", ]
    res <- reliability_H(fit, factor = f, ci_method = "delta")
    expect_equal(res$value[1:2], c(ref$est, ref$se), tolerance = 1e-8)
  }
  one <- lavaan::cfa(verbal_model, data = holzinger_swineford, std.lv = TRUE)
  det <- attr(factor_scores(one), "determinacy")
  expect_equal(det$determinacy_regression^2, reliability_H(one)$value[1],
               tolerance = 1e-10)
})

test_that("H does not depend on how the factor is identified", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  for (estimator in c("ML", "MLR")) {
    unit <- lavaan::cfa(verbal_model, data = holzinger_swineford,
                        std.lv = TRUE, estimator = estimator)
    marker <- lavaan::cfa(verbal_model, data = holzinger_swineford,
                          estimator = estimator)
    a <- reliability_H(unit)
    b <- reliability_H(marker)
    expect_equal(a$value[a$term == "estimate"], b$value[b$term == "estimate"],
                 tolerance = 1e-7)
    expect_equal(a$value[a$term == "se"], b$value[b$term == "se"],
                 tolerance = 1e-7)
    expect_equal(attr(a, "loadings"), attr(b, "loadings"), tolerance = 1e-7)
  }
})

test_that("the full covariance standard error differs from the diagonal one", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  visual <- c("t1_visual_perception", "t2_cubes", "t4_lozenges")
  ratio <- function(items) {
    fit <- lavaan::cfa(paste("f =~", paste(items, collapse = " + ")),
                       data = holzinger_swineford, std.lv = TRUE)
    full <- reliability_H(fit)
    std <- lavaan::standardizedSolution(fit)
    lam <- std[std$op == "=~", ]
    expect_warning(
      diag_tbl <- reliability_H(loadings = lam$est.std, se_loadings = lam$se,
                                ci_method = "delta_logistic"),
      class = "dmar_H_marginal_se")
    expect_equal(full$value[1], diag_tbl$value[1], tolerance = 1e-10)
    full$value[full$term == "se"] / diag_tbl$value[diag_tbl$term == "se"]
  }
  # The verbal tests' standardized loadings are positively correlated, so
  # the diagonal understates the variance; in the saturated three-item
  # visual block it overstates it.
  expect_gt(ratio(verbal), 1.2)
  expect_lt(ratio(visual), 0.8)
})

test_that("loadings with their full covariance matrix reproduce the fit path", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa(three_factor_model, data = holzinger_swineford,
                     std.lv = TRUE)
  items <- c("t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  nm <- paste0("textual=~", items)
  V <- lavaan::lavInspect(fit, "vcov.std.all")[nm, nm]
  std <- lavaan::standardizedSolution(fit)
  lam <- std$est.std[std$op == "=~" & std$lhs == "textual"]
  from_loadings <- reliability_H(loadings = stats::setNames(lam, items),
                                 vcov_loadings = V, N = 301)
  from_fit <- reliability_H(fit, factor = "textual")
  expect_equal(from_loadings$value[1:4], from_fit$value[1:4],
               tolerance = 1e-8)
  expect_identical(attr(from_loadings, "vcov_basis"), "full")
  expect_equal(attr(from_fit, "loadings"), stats::setNames(lam, items),
               tolerance = 1e-10)
  expect_identical(attr(from_fit, "factor"), "textual")
})

test_that("residual covariances inside the block give Gabler and Raykov's maximal reliability", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa(paste(verbal_model,
                           "\n t6_paragraph_comprehension ~~ t7_sentence"),
                     data = holzinger_swineford, std.lv = TRUE)
  res <- reliability_H(fit, ci_method = "delta")
  est <- lavaan::lavInspect(fit, "est")
  lam <- est$lambda[, 1]
  w <- solve(est$theta, lam)
  Sigma <- lavaan::lavInspect(fit, "implied")$cov
  composite_reliability <- drop((t(w) %*% lam)^2 * est$psi[1, 1] /
                                  (t(w) %*% Sigma %*% w))
  expect_equal(res$value[res$term == "estimate"], composite_reliability,
               tolerance = 1e-10)
  # An independent numerical delta method over the unstandardized
  # parameters (factor variance fixed at 1).
  x0 <- lavaan::coef(fit)
  V <- lavaan::vcov(fit)
  H_at <- function(x) {
    l <- x[paste0("verbal=~", verbal)]
    Th <- diag(x[paste0(verbal, "~~", verbal)])
    Th[2, 3] <- Th[3, 2] <- x["t6_paragraph_comprehension~~t7_sentence"]
    q <- drop(t(l) %*% solve(Th, l))
    q / (1 + q)
  }
  g <- vapply(seq_along(x0), function(k) {
    xp <- x0
    xm <- x0
    xp[k] <- xp[k] + 1e-6
    xm[k] <- xm[k] - 1e-6
    (H_at(xp) - H_at(xm)) / 2e-6
  }, numeric(1))
  expect_equal(res$value[res$term == "se"], sqrt(drop(t(g) %*% V %*% g)),
               tolerance = 1e-6)
  expect_identical(attr(res, "residual_covariances"), 1L)
  expect_match(DMAR:::.dmar_tbl_footer(res)[["H"]], "Gabler & Raykov")
  # The Hancock and Mueller formula that ignores the covariance differs.
  l <- attr(res, "loadings")
  s <- sum(l^2 / (1 - l^2))
  expect_gt(abs(s / (1 + s) - composite_reliability), 1e-3)
})

test_that("ordered items give the same H under the delta and theta parameterizations", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  ordinal <- as.data.frame(lapply(holzinger_swineford[, verbal], function(x) {
    cut(x, stats::quantile(x, c(0, .25, .5, .75, 1)), include.lowest = TRUE,
        labels = FALSE)
  }))
  delta <- lavaan::cfa(verbal_model, data = ordinal, std.lv = TRUE,
                       ordered = verbal)
  theta <- lavaan::cfa(verbal_model, data = ordinal, std.lv = TRUE,
                       ordered = verbal, parameterization = "theta")
  a <- reliability_H(delta)
  b <- reliability_H(theta)
  expect_equal(a$value[1:4], b$value[1:4], tolerance = 1e-6)
  expect_identical(attr(a, "metric"), "latent_response")
  expect_identical(attr(b, "metric"), "latent_response")
  expect_match(DMAR:::.dmar_tbl_footer(a)[["H"]], "latent responses")
  # The standardized loadings route of vcov.std.all agrees.
  nm <- paste0("verbal=~", verbal)
  V <- lavaan::lavInspect(delta, "vcov.std.all")[nm, nm]
  l <- lavaan::standardizedSolution(delta)
  l <- l$est.std[l$op == "=~"]
  s <- sum(l^2 / (1 - l^2))
  g <- 2 * l / ((1 - l^2)^2 * (1 + s)^2)
  expect_equal(a$value[2], sqrt(drop(t(g) %*% V %*% g)), tolerance = 1e-6)
})

test_that("an endogenous factor of a structural model has its CFA value", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  sem_fit <- lavaan::sem(paste(three_factor_model,
                               "\n textual ~ visual \n speed ~ visual + textual"),
                         data = holzinger_swineford, std.lv = TRUE)
  cfa_fit <- lavaan::cfa(three_factor_model, data = holzinger_swineford,
                         std.lv = TRUE)
  a <- reliability_H(sem_fit, factor = "speed")
  b <- reliability_H(cfa_fit, factor = "speed")
  expect_equal(a$value[1:2], b$value[1:2], tolerance = 1e-5)
  # A first-order factor of a higher-order model, likewise.
  ho_fit <- lavaan::cfa(paste(three_factor_model,
                              "\n g =~ visual + textual + speed"),
                        data = holzinger_swineford, std.lv = TRUE)
  expect_equal(reliability_H(ho_fit, factor = "textual")$value[1],
               reliability_H(cfa_fit, factor = "textual")$value[1],
               tolerance = 1e-5)
})

test_that("a multiple group fit gives one H per named group", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fit <- lavaan::cfa(verbal_model, data = holzinger_swineford,
                     group = "school", group.equal = "loadings")
  labels <- lavaan::lavInspect(fit, "group.label")
  nobs <- lavaan::lavInspect(fit, "nobs")
  for (g in seq_along(labels)) {
    res <- reliability_H(fit, group = labels[g])
    expect_identical(attr(res, "group"), labels[g])
    expect_equal(res$value[res$term == "N"], nobs[g])
    one <- lavaan::cfa(verbal_model,
                       data = holzinger_swineford[
                         holzinger_swineford$school == labels[g], ])
    expect_gt(res$value[1], 0.85)
    expect_lt(abs(res$value[1] - reliability_H(one)$value[1]), 0.01)
  }
  expect_match(DMAR:::.dmar_tbl_footer(res)[["H"]], "group")
  expect_error(reliability_H(fit, ci_method = "percentile", B = 10,
                             group = labels[1]), "single-group")
})

test_that("labelled equalities give lavaan's own delta method SE under either ceq.simple setting", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # Under ceq.simple = TRUE parameters that share a label share one free
  # number, while vcov() keeps a row for each of them, so the variance of
  # H must read vcov() by row. Indexed by free number it read the wrong
  # rows: 0.00674 for the one-group model below against lavaan's 0.01274,
  # and 0.00808 for the second group of the two-group model against
  # 0.01275.
  hs <- holzinger_swineford
  tau <- paste0(
    "verbal =~ ", paste0("a*", verbal, collapse = " + "), "\n",
    paste0(verbal, " ~~ e", 1:5, "*", verbal, collapse = "\n"), "\n",
    "q := a^2 * (1/e1 + 1/e2 + 1/e3 + 1/e4 + 1/e5)\n",
    "H := q / (1 + q)")
  two <- paste0(
    "verbal =~ ", paste0("c(a", 1:5, ", a", 1:5, ")*", verbal,
                         collapse = " + "), "\n",
    paste0(verbal, " ~~ c(u", 1:5, ", v", 1:5, ")*", verbal,
           collapse = "\n"), "\n",
    "verbal ~~ c(1, NA)*verbal + c(p1, p2)*verbal\n",
    "q1 := a1^2/u1 + a2^2/u2 + a3^2/u3 + a4^2/u4 + a5^2/u5\n",
    "H1 := q1 / (1 + q1)\n",
    "q2 := p2 * (a1^2/v1 + a2^2/v2 + a3^2/v3 + a4^2/v4 + a5^2/v5)\n",
    "H2 := q2 / (1 + q2)")
  lavaan_H <- function(fit, name) {
    pe <- lavaan::parameterEstimates(fit)
    unlist(pe[pe$lhs == name, c("est", "se")], use.names = FALSE)
  }
  for (ceq in c(FALSE, TRUE)) {
    f1 <- lavaan::cfa(tau, data = hs, std.lv = TRUE, ceq.simple = ceq)
    expect_equal(reliability_H(f1, ci_method = "delta")$value[1:2],
                 lavaan_H(f1, "H"), tolerance = 1e-7)
    f2 <- lavaan::cfa(two, data = hs, group = "school", std.lv = TRUE,
                      ceq.simple = ceq)
    labels <- lavaan::lavInspect(f2, "group.label")
    for (g in 1:2) {
      expect_equal(
        reliability_H(f2, group = labels[g], ci_method = "delta")$value[1:2],
        lavaan_H(f2, paste0("H", g)), tolerance = 1e-7)
    }
  }
  # group.equal under either setting: the same model, so the same SE to
  # the precision the two optimizations share (it was 0.652 against
  # 0.0129 for the second group of this scalar model).
  m <- paste("verbal =~", paste(verbal, collapse = " + "))
  se <- vapply(c(FALSE, TRUE), function(ceq) {
    f <- lavaan::cfa(m, data = hs, group = "school", std.lv = TRUE,
                     group.equal = c("loadings", "intercepts"),
                     ceq.simple = ceq)
    vapply(lavaan::lavInspect(f, "group.label"), function(g) {
      reliability_H(f, group = g, ci_method = "delta")$value[2]
    }, numeric(1))
  }, numeric(2))
  expect_equal(se[, 1], se[, 2], tolerance = 1e-5)
})

# Refusals -----------------------------------------------------------------

test_that("fits H is not defined for are refused with the remedy named", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  hs <- holzinger_swineford
  cross <- lavaan::cfa("
    visual  =~ t1_visual_perception + t2_cubes + t4_lozenges + t9_word_meaning
    textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
    data = hs, std.lv = TRUE)
  expect_error(reliability_H(cross, factor = "visual"),
               "also load on another factor.*bifactor_indices")
  higher <- lavaan::cfa(paste(three_factor_model,
                              "\n g =~ visual + textual + speed"),
                        data = hs, std.lv = TRUE)
  expect_error(reliability_H(higher, factor = "g"), "higher-order")
  mimic <- lavaan::sem(paste(verbal_model, "\n verbal ~ age_years"), data = hs)
  expect_error(reliability_H(mimic), "observed covariate")
  regressed <- lavaan::sem(paste(verbal_model, "\n t5_general_information ~ age_years"),
                           data = hs)
  expect_error(reliability_H(regressed), "regressed on other variables")
  three <- lavaan::cfa(three_factor_model, data = hs, std.lv = TRUE)
  expect_error(reliability_H(three), "3 factors.*name one with 'factor'")
  expect_error(reliability_H(three, factor = "nope"), "one latent variable")
  groups <- lavaan::cfa(verbal_model, data = hs, group = "school")
  expect_error(reliability_H(groups), "2 groups.*name one with 'group'")
  expect_error(reliability_H(groups, group = "Nowhere"), "group labels")
  expect_error(reliability_H(three, factor = "textual", group = "Pasteur"),
               "several groups")
  no_se <- lavaan::cfa(verbal_model, data = hs, std.lv = TRUE, se = "none")
  expect_error(reliability_H(no_se), "no standard errors")
  expect_equal(reliability_H(no_se, ci_method = "none")$value[1],
               reliability_H(lavaan::cfa(verbal_model, data = hs))$value[1],
               tolerance = 1e-7)
  ml_fit <- lavaan::cfa(verbal_model, data = hs, std.lv = TRUE)
  expect_error(reliability_H(ml_fit, ci_method = "mlr"),
               "se = \"standard\".*delta_logistic")
  expect_equal(reliability_H(ml_fit, ci_method = "ml")$value[1:4],
               reliability_H(ml_fit, ci_method = "delta")$value[1:4])
  expect_error(reliability_H(fit = data.frame(a = 1)), "fitted lavaan model")
  expect_error(reliability_H(fit = ml_fit, data = hs[, verbal]),
               "exactly one")
})

test_that("a loading or regression fixed at zero is not one, so the 0* model gets the H of the model without it", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  hs <- holzinger_swineford
  same_H <- function(a, b, factor) {
    ha <- reliability_H(a, factor = factor)
    hb <- reliability_H(b, factor = factor)
    expect_equal(ha$value, hb$value, tolerance = 1e-5)
    expect_identical(attr(ha, "residual_covariances"),
                     attr(hb, "residual_covariances"))
    expect_identical(names(attr(ha, "loadings")), names(attr(hb, "loadings")))
    ha
  }
  # A cross-loading fixed at zero. Its row made t7 an item of both
  # factors, and the fit was refused for either one as an item that
  # "also loads on another factor".
  plain <- lavaan::cfa("
    vis  =~ t1_visual_perception + t2_cubes + t4_lozenges
    verb =~ t5_general_information + t6_paragraph_comprehension + t7_sentence",
    data = hs, std.lv = TRUE)
  zero <- lavaan::cfa("
    vis  =~ t1_visual_perception + t2_cubes + t4_lozenges + 0*t7_sentence
    verb =~ t5_general_information + t6_paragraph_comprehension + t7_sentence",
    data = hs, std.lv = TRUE)
  expect_equal(lavaan::fitMeasures(zero, "chisq"),
               lavaan::fitMeasures(plain, "chisq"), tolerance = 1e-6)
  h <- same_H(zero, plain, "verb")
  expect_equal(h$value[1:2], c(0.88325, 0.01260), tolerance = 1e-4)
  h <- same_H(zero, plain, "vis")
  expect_identical(h$value[h$term == "n_items"], 3)
  # An item listed under the factor alone with its loading fixed at zero
  # is not one of its items; with a residual covariance to an item of the
  # block it entered the optimal weights and gave 0.641 (J = 4), where the
  # same model written without it gives 0.637.
  listed <- lavaan::cfa("
    vis =~ t1_visual_perception + t2_cubes + t4_lozenges + 0*t3_paper_form_board
    t3_paper_form_board ~~ t1_visual_perception", data = hs, std.lv = TRUE)
  unlisted <- lavaan::cfa("
    vis =~ t1_visual_perception + t2_cubes + t4_lozenges
    t3_paper_form_board ~~ t1_visual_perception", data = hs, std.lv = TRUE)
  expect_equal(lavaan::fitMeasures(listed, "chisq"),
               lavaan::fitMeasures(unlisted, "chisq"), tolerance = 1e-6)
  h <- same_H(listed, unlisted, "vis")
  expect_identical(h$value[h$term == "n_items"], 3)
  expect_identical(attr(h, "residual_covariances"), 0L)
  expect_equal(h$value[1], 0.6369, tolerance = 1e-3)
  expect_false(any(grepl("residual covariance",
                         DMAR:::.dmar_tbl_footer(reliability_H(listed)))))
  # A regression of an item fixed at zero is not a regression.
  fixed_reg <- lavaan::cfa("
    vis =~ t1_visual_perception + t2_cubes + t4_lozenges
    t1_visual_perception ~ 0*t3_paper_form_board", data = hs, std.lv = TRUE)
  no_reg <- lavaan::cfa("
    vis =~ t1_visual_perception + t2_cubes + t4_lozenges
    t3_paper_form_board ~~ t3_paper_form_board", data = hs, std.lv = TRUE)
  same_H(fixed_reg, no_reg, "vis")
  # A free cross-loading is still refused, and the message says that a
  # loading fixed or constrained to zero would not count.
  cross <- lavaan::cfa("
    vis  =~ t1_visual_perception + t2_cubes + t4_lozenges + t7_sentence
    verb =~ t5_general_information + t6_paragraph_comprehension + t7_sentence",
    data = hs, std.lv = TRUE)
  expect_error(reliability_H(cross, factor = "verb"),
               "also load on another factor.*constrained to zero does not count")
})

test_that("a loading or regression constrained to zero is one fixed at zero, so a*x with a == 0 gets the H of the model without it", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  hs <- holzinger_swineford
  same_H <- function(a, b, factor) {
    ha <- reliability_H(a, factor = factor)
    hb <- reliability_H(b, factor = factor)
    expect_equal(ha$value, hb$value, tolerance = 1e-5)
    expect_identical(names(attr(ha, "loadings")), names(attr(hb, "loadings")))
    ha
  }
  plain <- lavaan::cfa("
    vis  =~ t1_visual_perception + t2_cubes + t4_lozenges
    verb =~ t5_general_information + t6_paragraph_comprehension + t7_sentence",
    data = hs, std.lv = TRUE)
  # A labeled cross-loading constrained to zero is the model written with
  # 0*t7_sentence and the model written without it (the same chi square
  # and degrees of freedom), but lavaan estimates it as a free parameter,
  # so a rule that read only rows fixed at zero refused the fit for
  # either factor as one whose item "also loads on another factor".
  eq <- lavaan::cfa("
    vis  =~ t1_visual_perception + t2_cubes + t4_lozenges + a*t7_sentence
    verb =~ t5_general_information + t6_paragraph_comprehension + t7_sentence
    a == 0", data = hs, std.lv = TRUE)
  expect_equal(lavaan::fitMeasures(eq, c("chisq", "df")),
               lavaan::fitMeasures(plain, c("chisq", "df")), tolerance = 1e-6)
  h <- same_H(eq, plain, "verb")
  expect_equal(h$value[1:2], c(0.88325, 0.01260), tolerance = 1e-4)
  h <- same_H(eq, plain, "vis")
  expect_identical(h$value[h$term == "n_items"], 3)
  # Through a chain of equalities: the loading equals a regression
  # constrained to zero.
  chain <- lavaan::cfa("
    vis  =~ t1_visual_perception + t2_cubes + t4_lozenges + a*t7_sentence
    verb =~ t5_general_information + t6_paragraph_comprehension + t7_sentence
    t6_paragraph_comprehension ~ b*t13_straight_and_curved_capitals
    a == b
    b == 0", data = hs, std.lv = TRUE)
  chain_plain <- lavaan::cfa("
    vis  =~ t1_visual_perception + t2_cubes + t4_lozenges
    verb =~ t5_general_information + t6_paragraph_comprehension + t7_sentence
    t13_straight_and_curved_capitals ~~ t13_straight_and_curved_capitals",
    data = hs, std.lv = TRUE)
  expect_equal(lavaan::fitMeasures(chain, c("chisq", "df")),
               lavaan::fitMeasures(chain_plain, c("chisq", "df")),
               tolerance = 1e-6)
  same_H(chain, chain_plain, "verb")
  same_H(chain, chain_plain, "vis")
  # A constraint to a nonzero value, or one that is not a plain equality
  # with a constant, leaves the loading a loading.
  half <- lavaan::cfa("
    vis  =~ t1_visual_perception + t2_cubes + t4_lozenges + a*t7_sentence
    verb =~ t5_general_information + t6_paragraph_comprehension + t7_sentence
    a == 0.5", data = hs, std.lv = TRUE)
  expect_error(reliability_H(half, factor = "verb"),
               "also load on another factor")
})

test_that(".pt_zero_rows() reads rows fixed at zero and labels constrained to zero, and nothing else", {
  pt <- data.frame(
    lhs = c("f", "f", "f", "f", "f", "f", "f", "y",
            "a", "0", "c", ".p6.", "e", "d", "h"),
    op = c("=~", "=~", "=~", "=~", "=~", "=~", "=~", "~",
           "==", "==", "==", "==", "==", "==", ":="),
    rhs = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8",
            "0", "b", "a", "0", "1", "2*a", "a*b"),
    free = c(0L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, rep(0L, 7)),
    est = c(1, 0, 0, 0, 0, 0, 1, 0, rep(0, 7)),
    label = c("", "", "a", "b", "c", "", "e", "d", rep("", 6), "h"),
    plabel = c(paste0(".p", 1:8, "."), rep("", 6), ".p15."),
    stringsAsFactors = FALSE)
  # x1 is fixed at one (a marker loading) and x2 at zero; x3 is
  # constrained to zero (a == 0), x4 through 0 == b written the other way
  # round, x5 through the chain c == a, and x6, which has no label, through
  # its plabel (.p6. == 0); x7 is constrained to one, and the regression
  # on x8 to 2*a, which is not a plain equality and is not followed. The
  # constraint and defined parameter rows are not parameters held at zero.
  expect_identical(DMAR:::.pt_zero_rows(pt),
                   c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
                     rep(FALSE, 7)))
  # Without constraints, only the row fixed at zero.
  expect_identical(DMAR:::.pt_zero_rows(pt[pt$op != "==", ]),
                   c(FALSE, TRUE, rep(FALSE, 7)))
})

test_that("an improper solution is refused on every path", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # lavaan warns about the negative variance before the refusal.
  R <- heywood_S()
  expect_error(suppressWarnings(reliability_H(S = R, N = 200)), "improper")
  fit <- suppressWarnings(lavaan::cfa("f =~ y1 + y2 + y3", sample.cov = R,
                                      sample.nobs = 200))
  expect_error(reliability_H(fit), class = "dmar_H_improper")
  # Raw data with the same correlations.
  z <- withr::with_seed(113, matrix(stats::rnorm(600), 200, 3))
  z <- scale(z, scale = FALSE)
  z <- z %*% solve(chol(stats::cov(z))) %*% chol(R)
  colnames(z) <- colnames(R)
  expect_error(suppressWarnings(reliability_H(data = z)), "improper")
})

test_that("a fit lavaan stops with an error is reported with lavaan's message, not as one that did not converge", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  skip_if_not_installed("lavaan")
  # An error lavaan stops the fit with, which the data path reported as
  # "The single-factor model did not converge.", is quoted, as
  # reliability_omega() quotes it. (A constant item, the example this test
  # used, is now refused by name before any fit; lavaan's own error named
  # no item and came after its table of the variables.)
  d <- holzinger_swineford[, verbal]
  d$t7_sentence <- 5
  out <- utils::capture.output(
    err <- tryCatch(reliability_H(data = d), error = conditionMessage))
  expect_length(out, 0L)
  expect_match(err, "^t7_sentence takes a single value in the 301 complete rows")
  ch <- holzinger_swineford[, verbal]
  ch$t7_sentence <- as.character(ch$t7_sentence)
  expect_match(
    .relia_fit_failure(.omega_fit_cfa(data = ch, se = "none"),
                       model = "single-factor"),
    paste0("^The single-factor model could not be fit; lavaan stopped ",
           "with this error: .*of type character: \"t7_sentence\""))
  # A fit that ran and did not converge keeps its message (through the
  # first versions of 1.1.0 fa_confirmatory()'s own stop on it was quoted as
  # lavaan's error), with lavaan's convergence flag set to FALSE.
  d <- holzinger_swineford[, verbal]
  d$t6_paragraph_comprehension[1:30] <- NA
  real <- lavaan::lavInspect
  testthat::local_mocked_bindings(
    lavInspect = function(object, what = "free", ...) {
      if (identical(what, "converged")) FALSE else real(object, what, ...)
    }, .package = "lavaan")
  expect_error(suppressWarnings(reliability_H(data = d,
                                              missing = "listwise")),
               "^The single-factor model did not converge\\.$")
  expect_error(suppressWarnings(reliability_H(data = d, missing = "fiml")),
               paste0("The single-factor model did not converge.",
                      .relia_convergence_note("fiml")), fixed = TRUE)
})

test_that("interval methods and arguments that do not apply are refused", {
  skip_on_cran()  # The extended validation grid runs in the full local suite.
  skip_if_not_installed("lavaan")
  items <- holzinger_swineford[, verbal]
  expect_error(reliability_H(data = items, ci_method = "delta"),
               "does not apply to 'data'")
  expect_error(reliability_H(S = stats::cov(items), N = 301,
                             ci_method = "mlr"), "does not apply to 'S'")
  expect_error(reliability_H(S = stats::cov(items), N = 301,
                             ci_method = "bca"), "need raw 'data'")
  expect_error(reliability_H(loadings = c(.6, .7, .8), ci_method = "ml"),
               "does not apply to 'loadings'")
  expect_error(reliability_H(loadings = c(.6, .7, .8), N = 100,
                             ci_method = "delta_logistic"),
               "vcov_loadings")
  expect_error(reliability_H(loadings = c(.6, .7, .8), vcov_loadings = diag(2)),
               "3 x 3")
  expect_error(reliability_H(loadings = c(.6, .7, .8), vcov_loadings = diag(3),
                             se_loadings = rep(.1, 3)), "not both")
  expect_error(reliability_H(data = items, N = 301), "'N' does not apply")
  expect_error(reliability_H(data = items, factor = "f1"),
               "'factor' does not apply")
  expect_error(reliability_H(data = items[, 1:2]), "three or more items")
  expect_error(reliability_H(S = stats::cov(items), N = 301,
                             missing = "fiml"), "requires raw 'data'")
  expect_error(reliability_H(data = mar_fixture(), aux = "t1_visual_perception",
                             missing = "listwise"), "cannot be combined")
  expect_error(reliability_H(data = items, conf_level = 95), "conf_level")
  expect_error(reliability_H(data = items, ci_method = "percentile", B = 0.5),
               "'B'")
  expect_error(reliability_H(), "exactly one")
})

# The data path ------------------------------------------------------------

test_that("the data path rests on the same fit as omega, and FIML uses the incomplete rows", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- mar_fixture()
  listwise <- reliability_H(data = d[, verbal], missing = "listwise")
  fiml <- reliability_H(data = d[, verbal], missing = "fiml")
  # The default is FIML for these incomplete items.
  expect_equal(reliability_H(data = d[, verbal])$value, fiml$value)
  aux <- reliability_H(data = d, aux = "t1_visual_perception")
  expect_identical(attr(listwise, "ci_method"), "mlr_logistic")
  expect_identical(attr(listwise, "se_type"), "robust.sem")
  expect_identical(attr(fiml, "se_type"), "robust.huber.white")
  expect_equal(listwise$value[listwise$term %in% c("N", "N_complete")],
               c(167, 167))
  expect_equal(fiml$value[fiml$term %in% c("N", "N_complete")], c(301, 167))
  expect_equal(aux$value[aux$term %in% c("N", "N_complete")], c(301, 167))
  expect_gt(abs(fiml$value[1] - listwise$value[1]), 0.005)
  expect_lt(fiml$value[2], listwise$value[2])
  expect_identical(attr(fiml, "missing"), "fiml")
  expect_identical(attr(aux, "aux"), "t1_visual_perception")
  expect_identical(names(attr(aux, "loadings")), verbal)
  # The same single-factor fit, through the fit path.
  same <- fa_confirmatory(data = d[, verbal], factors = list(f1 = verbal),
                estimator = "MLR", se = "robust.huber.white",
                missing = "ml", output = "fit")
  expect_equal(fiml$value[1:4],
               reliability_H(same, ci_method = "delta_logistic")$value[1:4],
               tolerance = 1e-8)
  # Covariance input: the maximum likelihood fit of the same model.
  hs <- holzinger_swineford
  from_S <- reliability_H(S = stats::cov(hs[, verbal]), N = 301)
  from_fit <- reliability_H(lavaan::cfa(verbal_model, data = hs,
                                        std.lv = TRUE))
  expect_equal(from_S$value[1:4], from_fit$value[1:4], tolerance = 1e-6)
  expect_identical(attr(from_S, "ci_method"), "ml_logistic")
})

# Bootstrap ------------------------------------------------------------------

test_that("the data path bootstrap is reproducible, restores the RNG, and estimates bias", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("boot")
  d <- mar_fixture()[, verbal]
  set.seed(20261101)
  before <- get(".Random.seed", envir = globalenv())
  a <- reliability_H(data = d, missing = "fiml", ci_method = "bootstrap_se",
                     B = 20, seed = 113)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  b <- reliability_H(data = d, missing = "fiml", ci_method = "bootstrap_se",
                     B = 20, seed = 113)
  expect_identical(a$value, b$value)
  bs <- attr(a, "bootstrap")
  expect_identical(bs$attempted, 20L)
  expect_identical(bs$returned + bs$failed, 20L)
  expect_gt(a$value[a$term == "se"], 0)
  expect_true(is.finite(bs$bias))
  expect_identical(attr(a, "B"), 20L)
  expect_identical(attr(a, "vcov_basis"), "none")
  expect_lt(a$value[a$term == "ci_lower"], a$value[1])
  expect_gt(a$value[a$term == "ci_upper"], a$value[1])
  expect_match(DMAR:::.dmar_tbl_footer(a)[["H_bias"]], "bias of H")
})

test_that("the fit path bootstrap refits the user's model", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("boot")
  fit <- lavaan::cfa(verbal_model, data = mar_fixture(), std.lv = TRUE,
                     missing = "ml", estimator = "MLR")
  # A refit on the fitted rows reproduces the estimate.
  h <- DMAR:::.H_from_fit(fit, se = FALSE)
  refit <- DMAR:::.H_refit_fun(fit, h$setup)
  expect_equal(refit(DMAR:::.H_fit_data(fit)), h$H, tolerance = 1e-8)
  # The refit starts from the fit's estimates, which lavaan reads from the
  # est column of the parameter table, so on the fitted rows it has
  # nothing left to iterate; from lavaan's own start values (a table
  # without est or start) the same refit takes many iterations.
  lavaan_fit <- lavaan::lavaan
  seen <- new.env(parent = emptyenv())
  local_mocked_bindings(lavaan = function(...) {
    r <- lavaan_fit(...)
    seen$iterations <- lavaan::lavInspect(r, "iterations")
    r
  }, .package = "lavaan")
  refit(DMAR:::.H_fit_data(fit))
  pt <- lavaan::parTable(fit)
  pt <- pt[setdiff(names(pt), c("est", "se", "start"))]
  opts <- lavaan::lavInspect(fit, "options")
  opts[c("se", "test", "baseline", "h1")] <- list("none", "none", FALSE,
                                                   FALSE)
  from_defaults <- lavaan::lavInspect(lavaan_fit(
    model = pt, data = DMAR:::.H_fit_data(fit), slotOptions = opts),
    "iterations")
  expect_lte(seen$iterations, 2L)
  expect_gt(from_defaults, 5L)
  a <- reliability_H(fit, ci_method = "bootstrap_se_logistic", B = 20,
                     seed = 113)
  b <- reliability_H(fit, ci_method = "bootstrap_se_logistic", B = 20,
                     seed = 113)
  expect_identical(a$value, b$value)
  expect_gt(a$value[a$term == "se"], 0)
  expect_identical(attr(a, "bootstrap")$attempted, 20L)
  expect_true(is.finite(attr(a, "bootstrap")$bias))
  expect_identical(attr(a, "factor"), "verbal")
  summary_fit <- lavaan::cfa(verbal_model,
                             sample.cov = stats::cov(holzinger_swineford[, verbal]),
                             sample.nobs = 301)
  expect_error(reliability_H(summary_fit, ci_method = "percentile", B = 10),
               "raw data")
})

test_that("a fit path replication with too few distinct rows is refused under the fit's own treatment of incomplete rows", {
  skip_on_cran()  # One hundred refits run in the full local suite.
  skip_if_not_installed("lavaan")
  skip_if_not_installed("boot")
  # Regression test. Every fit was checked as a FIML fit, so a replication
  # of a fit to complete data was refused in the words of the saturated
  # model that full information maximum likelihood estimates first, with
  # "the data would be refused for the same reason"; no EM step runs for
  # such a fit, and reliability_H(fit = ) does not refuse its data on this
  # rule. Nine pupils, five verbal tests, complete data (listwise, ML).
  d <- holzinger_swineford[1:9, verbal]
  fit <- lavaan::cfa(verbal_model, data = d, std.lv = TRUE)
  r <- suppressWarnings(reliability_H(fit = fit, ci_method = "percentile",
                                      B = 100, seed = 1))
  msgs <- names(attr(r, "bootstrap")$error_messages)
  listwise <- function(n) {
    paste0("Only ", n, " distinct rows of this resample, no more than the ",
           "5 variables of the model, so their sample covariance matrix is ",
           "singular (its rank is at most ", n - 1L, ").")
  }
  expect_true(listwise(4L) %in% msgs)
  expect_true(listwise(5L) %in% msgs)
  expect_false(any(grepl("full information|would be refused", msgs)))
})

test_that("the fit path checks distinct rows under its own missing data treatment", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[1:9, verbal]
  fit <- lavaan::cfa(verbal_model, data = d, std.lv = TRUE)
  listwise <- function(n) {
    paste0("Only ", n, " distinct rows of this resample, no more than the ",
           "5 variables of the model, so their sample covariance matrix is ",
           "singular (its rank is at most ", n - 1L, ").")
  }
  # A FIML fit is held to the rule of the saturated model, in words that
  # describe the resample alone.
  fiml <- lavaan::cfa(verbal_model, data = mar_fixture()[, verbal],
                      std.lv = TRUE, missing = "ml")
  chk <- function(fit, rows) {
    h <- DMAR:::.H_from_fit(fit, se = FALSE)
    refit <- DMAR:::.H_refit_fun(fit, h$setup)
    tryCatch(refit(DMAR:::.H_fit_data(fit)[rows, ]),
             error = conditionMessage)
  }
  expect_identical(
    chk(fiml, c(1:5, 5)),
    paste0("Only 5 distinct cases of this resample have an observed model ",
           "variable, no more than the 5 model variables. It does not meet ",
           "DMAR's distinct-case estimability screen for incomplete-data ",
           "resamples. This conservative screen does not establish ",
           "likelihood nonexistence."))
  expect_identical(chk(fit, c(1:5, 5)), listwise(5L))
  # A fit by pairwise deletion, whose covariances come from different
  # rows, is held to neither rule.
  expect_true(DMAR:::.relia_check_resample(
    d[c(1:5, 5), ], missing = "pairwise", items = verbal, fitted = TRUE,
    on_fit = "variables"))
})

test_that("the fit path bootstrap refuses weights and clusters and checks that its refit reproduces H", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("boot")
  d <- holzinger_swineford[, verbal]
  d$w <- ifelse(d$t5_general_information >
                  stats::median(d$t5_general_information), 5, 0.2)
  # The weights did not reach the refits, so every replication was the
  # unweighted H: the percentile interval [0.902, 0.932] lay above the
  # weighted estimate 0.889, and the reported bias, 0.030, was the gap
  # between the two.
  weighted <- lavaan::cfa(verbal_model, data = d, std.lv = TRUE,
                          sampling.weights = "w")
  expect_error(reliability_H(weighted, ci_method = "percentile", B = 10,
                             seed = 113),
               "sampling weights.*ci_method = \"delta_logistic\"")
  expect_equal(reliability_H(weighted)$value[1],
               DMAR:::.H_from_fit(weighted, se = FALSE)$H)
  # A cluster fit's rows were resampled as independent cases.
  d$cl <- rep(1:30, length.out = nrow(d))
  clustered <- lavaan::cfa(verbal_model, data = d, std.lv = TRUE,
                           cluster = "cl")
  expect_error(reliability_H(clustered, ci_method = "bca", B = 10),
               "cluster variable \\(cl\\)")
  # Before any resample is drawn, the refit of the fit's own rows must
  # reproduce the estimate.
  fit <- lavaan::cfa(verbal_model, data = holzinger_swineford,
                     std.lv = TRUE)
  local_mocked_bindings(.H_refit_fun = function(fit, s) function(d) 0.5)
  expect_error(reliability_H(fit, ci_method = "percentile", B = 10),
               "did not reproduce its H \\(0\\.500000 against 0\\.918")
  local_mocked_bindings(
    .H_refit_fun = function(fit, s) function(d) stop("no refit"))
  expect_error(reliability_H(fit, ci_method = "percentile", B = 10),
               "could not be refit to the data stored with 'fit': no refit")
})

test_that("a resample with no single-factor solution is counted as improper, as omega counts it", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("boot")
  # Three items, one of them weak: a resample whose three correlations
  # have a product that is not positive admits no single-factor
  # solution. H counted 19 of these 60 as non-finite while
  # reliability_omega() counted them as improper, and then counted them
  # as improper without the `no_solution` count omega gives them.
  d <- withr::with_seed(1, {
    n <- 120
    f <- stats::rnorm(n)
    data.frame(y1 = 0.8 * f + stats::rnorm(n), y2 = 0.8 * f + stats::rnorm(n),
               y3 = 0.12 * f + stats::rnorm(n))
  })
  w <- testthat::capture_warnings(
    h <- reliability_H(data = d, ci_method = "percentile", B = 60, seed = 1))
  o <- suppressWarnings(
    reliability_omega(data = d, ci_method = "percentile", B = 60, seed = 1))
  bh <- attr(h, "bootstrap")
  expect_identical(bh$status, attr(o, "bootstrap")$status)
  expect_identical(sum(bh$status == "non_finite"), 0L)
  expect_gt(sum(bh$status == "improper"), 0L)
  expect_identical(bh$improper, sum(bh$status == "improper"))
  expect_identical(bh$improper_kept, 0L)
  expect_gt(bh$no_solution, 0L)
  expect_identical(bh$no_solution, attr(o, "bootstrap")$no_solution)
  # H is not defined for an improper solution, and the record says so.
  expect_identical(bh$improper_rule, "drop")
  expect_match(w, "admit no single-factor solution. H is not defined",
               fixed = TRUE, all = FALSE)
  expect_match(w, "the second belongs on neither side of it", fixed = TRUE,
               all = FALSE)
  expect_match(w, sprintf("%d had no single-factor solution", bh$no_solution),
               fixed = TRUE, all = FALSE)
  expect_false(any(grepl("improper = \"keep\"", w)))
})

test_that("percentile and BCa intervals rest on the same replications on the data and fit paths", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("boot")
  # 80 rows, so the BCa acceleration comes from the jackknife (40
  # replications are fewer than the cases) and the test stays quick; an
  # 80 percent interval keeps its limits off the extreme replicates.
  d <- holzinger_swineford[1:80, verbal]
  fit <- lavaan::cfa(verbal_model, data = d, std.lv = TRUE)
  for (m in c("percentile", "bca")) {
    a <- reliability_H(data = d, ci_method = m, B = 40, seed = 3,
                       conf_level = 0.8)
    b <- reliability_H(fit, ci_method = m, B = 40, seed = 3,
                       conf_level = 0.8)
    expect_equal(b$value[1:4], a$value[1:4], tolerance = 1e-6)
    expect_identical(attr(b, "ci_method"), m)
    bs <- attr(b, "bootstrap")
    expect_identical(bs$attempted, 40L)
    expect_equal(bs$bias, mean(bs$t[bs$status == "proper"]) - b$value[1],
                 tolerance = 1e-12)
    expect_lt(b$value[b$term == "ci_lower"], b$value[1])
    expect_gt(b$value[b$term == "ci_upper"], b$value[1])
  }
  expect_identical(attr(b, "bootstrap")$acceleration_method, "jackknife")
})

# The dispatcher and the broom verbs -----------------------------------------

test_that("reliability(type = \"H\") dispatches to reliability_H()", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  items <- holzinger_swineford[, verbal]
  expect_identical(reliability(data = items, type = "H"),
                   reliability_H(data = items))
  expect_identical(reliability(S = stats::cov(items), N = 301, type = "H"),
                   reliability_H(S = stats::cov(items), N = 301))
  expect_identical(
    reliability(data = mar_fixture()[, verbal], type = "H", missing = "fiml",
                ci_method = "ml"),
    reliability_H(data = mar_fixture()[, verbal], missing = "fiml",
                  ci_method = "ml"))
})

test_that("tidy() and glance() read the H table", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  res <- reliability_H(lavaan::cfa(three_factor_model,
                                   data = holzinger_swineford, std.lv = TRUE),
                       factor = "visual")
  td <- generics::tidy(res)
  gl <- generics::glance(res)
  expect_identical(td$term, "H")
  expect_equal(td$estimate, res$value[res$term == "estimate"])
  expect_equal(td$ci_upper, res$value[res$term == "ci_upper"])
  expect_identical(gl$coefficient, "H")
  expect_equal(gl$nobs, 301)
  expect_identical(gl$n_items, 3)
  expect_identical(gl$ci_method, "delta_logistic")
  expect_identical(DMAR:::.dmar_tbl_footer(res)[["H"]], "Factor: visual")
  # H takes no wald_df: its row is on no t reference.
  expect_identical(td$wald_df, Inf)
  expect_identical(gl$wald_df, Inf)
})

test_that("the tidy() and glance() rows of a reliability coefficient carry its t reference", {
  skip_on_cran()  # Extended shape and boundary checks run in the full local suite.
  # The t reference of a Wald interval was only an attribute, so the
  # rows of an interval on t(30) and one on the normal reference looked
  # alike once stacked, differing only in their limits.
  skip_if_not_installed("lavaan")
  items <- holzinger_swineford[1:40, c("t6_paragraph_comprehension",
                                       "t7_sentence", "t9_word_meaning")]
  z <- reliability_omega(items, denominator = "model_implied",
                         ci_method = "ml")
  t30 <- reliability_omega(items, denominator = "model_implied",
                           ci_method = "ml", wald_df = 30)
  both <- rbind(generics::tidy(z), generics::tidy(t30))
  expect_identical(both$wald_df, c(Inf, 30))
  expect_identical(generics::glance(t30)$wald_df, 30)
  a <- reliability_alpha(items, ci_method = "ml", wald_df = 30)
  expect_identical(generics::tidy(a)$wald_df, 30)
})

test_that("from data the default 'missing' is FIML for incomplete items, as in reliability_omega(), without a warning", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- holzinger_swineford[, verbal]
  d_na <- d
  d_na$t7_sentence[seq(2, 120, by = 3)] <- NA
  # Complete data and a covariance matrix are fitted as they stand.
  a <- reliability_H(data = d)
  expect_identical(attr(a, "missing"), "listwise")
  expect_identical(a$value, reliability_H(data = d, missing = "listwise")$value)
  expect_identical(attr(reliability_H(S = cov(d), N = nrow(d)), "missing"),
                   "listwise")
  # Incomplete items are fitted by FIML.
  expect_no_warning(x <- reliability_H(data = d_na))
  expect_identical(attr(x, "missing"), "fiml")
  expect_equal(x$value[x$term == "N"], nrow(d_na))
  expect_equal(x$value, reliability_H(data = d_na, missing = "fiml")$value)
  l <- reliability_H(data = d_na, missing = "listwise")
  expect_equal(l$value[l$term == "N"], sum(stats::complete.cases(d_na)))
  # The single-factor model is the one reliability_omega() fits by
  # default, so the two use the same cases.
  o <- suppressWarnings(reliability_omega(d_na,
                                          denominator = "model_implied"))
  expect_identical(attr(o, "missing"), attr(x, "missing"))
  # reliability() forwards an unsupplied 'missing' unsupplied.
  expect_equal(reliability(d_na, type = "H")$value, x$value)
})

# The covariance matrix of the estimates, used only when it can be --------

# Two verbal tests made uncorrelated in the sample with the three spatial
# tests of the other factor, so that their factor's correlation with it is
# 0 and their two loadings are not identified at the estimates: the
# fixture of the test of the singularity check in test-fa_confirmatory.R.
singular_fixture <- function() {
  hs <- holzinger_swineford
  A <- c("t1_visual_perception", "t2_cubes", "t4_lozenges")
  B <- c("t6_paragraph_comprehension", "t7_sentence")
  d <- hs[A]
  X <- cbind(1, as.matrix(d))
  for (v in B) d[[v]] <- as.numeric(hs[[v]] - X %*% qr.solve(X, hs[[v]])) + 10
  list(data = d, factors = list(a = A, b = B))
}

test_that("H from a fit whose information is singular keeps its estimate and withholds its interval, as fa_confirmatory() does", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  fx <- singular_fixture()
  # FIML (the default): lavaan inverts the singular observed information
  # without complaint. Through 1.1.0 development reliability_H() reported
  # H = .853 for the two-item factor with a standard error of 534 and the
  # interval [0, 1], with no warning, beside a fa_confirmatory() table of
  # the same fit that withheld every standard error.
  seen <- new.env(parent = emptyenv())
  keep <- function(slot) {
    function(w) {
      assign(slot, w, envir = seen)
      invokeRestart("muffleWarning")
    }
  }
  tab <- withCallingHandlers(
    fa_confirmatory(fx$data, factors = fx$factors, output = "measurement"),
    dmar_result_qualified = keep("tab_w"))
  tab_w <- seen$tab_w
  fit <- suppressWarnings(fa_confirmatory(fx$data, factors = fx$factors,
                                          output = "fit"))
  expect_true(isTRUE(lavaan::lavInspect(fit, "converged")))
  for (f in c("b", "a")) {
    seen$h_w <- NULL
    res <- withCallingHandlers(
      reliability_H(fit = fit, factor = f),
      dmar_result_qualified = keep("h_w"))
    h_w <- seen$h_w
    h_val <- function(t) res$value[res$term == t]
    # The estimate is the H row of the table, lavaan's defined parameter.
    expect_equal(h_val("estimate"), tab$estimate[tab$term == paste0("H_", f)],
                 tolerance = 1e-7)
    expect_true(all(is.na(c(h_val("se"), h_val("ci_lower"), h_val("ci_upper")))))
    # The same class and the same words as the table's diagnostic.
    expect_s3_class(h_w, "dmar_result_qualified")
    expect_identical(conditionMessage(h_w), conditionMessage(tab_w))
    expect_identical(attr(res, "vcov_basis"), "none")
  }
  # The last factor read was "a", whose H the fit also reports.
  expect_equal(h_val("estimate"), 0.6598382, tolerance = 1e-6)
  b <- suppressWarnings(reliability_H(fit = fit, factor = "b"))
  expect_equal(b$value[b$term == "estimate"], 0.8530609, tolerance = 1e-6)
  # A point estimate alone asks nothing of the covariance matrix.
  expect_no_warning(n <- reliability_H(fit = fit, factor = "b",
                                       ci_method = "none"))
  expect_equal(n$value[1], b$value[1])
})

test_that("H from a fit whose standard errors lavaan could not compute withholds its interval rather than failing", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # two more fits
  fx <- singular_fixture()
  # Listwise ML: lavaan cannot invert the expected information. Through
  # 1.1.0 development reliability_H() stopped with an R error ("argument
  # is of length zero") reading the empty covariance matrix.
  fit <- suppressWarnings(fa_confirmatory(fx$data, factors = fx$factors,
                                          missing = "listwise",
                                          output = "fit"))
  seen <- with_warning_classes(reliability_H(fit = fit, factor = "b"))
  res <- seen$value
  # The class of the tables' warning; its last words name the standard
  # error and interval of H rather than the tables' columns.
  expect_true("dmar_result_qualified" %in% seen$classes)
  expect_true(is.finite(res$value[res$term == "estimate"]))
  expect_true(all(is.na(res$value[res$term %in% c("se", "ci_lower",
                                                   "ci_upper")])))
})

test_that("a two-item factor its equal loadings identify keeps the delta method interval", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  hs <- holzinger_swineford
  f <- list(a = c("t1_visual_perception", "t2_cubes"),
            b = c("t6_paragraph_comprehension", "t7_sentence",
                  "t9_word_meaning"))
  fit <- fa_confirmatory(hs, f, equal_loading = c(a = TRUE, b = FALSE),
                         output = "fit")
  tab <- fa_confirmatory(hs, f, equal_loading = c(a = TRUE, b = FALSE),
                         output = "measurement")
  expect_no_warning(res <- reliability_H(fit = fit, factor = "a"))
  # lavaan's delta method standard error of its defined H, another code
  # path, agrees with the numerical gradient over the free parameters.
  expect_equal(res$value[res$term == "se"], tab$se[tab$term == "H_a"],
               tolerance = 1e-6)
  expect_equal(res$value[res$term == "estimate"],
               tab$estimate[tab$term == "H_a"], tolerance = 1e-8)
  expect_identical(attr(res, "vcov_basis"), "full")
})
