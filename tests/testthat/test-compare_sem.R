criterion_of <- function(fit, criterion, ...) {
  x <- sem_information_criteria(fit, measures = criterion, ...)
  x$value[x$term == criterion]
}

test_that("values, labels, weights, and evidence agree with lavaan", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  one <- lavaan::cfa("f =~ x1 + x2 + x3 + x4 + x5 + x6", data = d)
  two <- lavaan::cfa("f1 =~ x1 + x2 + x3
                      f2 =~ x4 + x5 + x6", data = d)
  x <- compare_sem(list(one = one, two = two))
  expect_s3_class(x, "dmar_compare_sem")
  expect_named(x, c("term", "value", "delta", "relative_weight",
                    "log_evidence_best_over_model", "npar"))
  bic <- c(one = lavaan::fitMeasures(one, "bic")[[1]],
           two = lavaan::fitMeasures(two, "bic")[[1]])
  expect_identical(x$term, names(sort(bic)))
  expect_equal(x$value, unname(sort(bic)))
  expect_equal(x$delta, unname(sort(bic) - min(bic)))
  expect_equal(x$relative_weight[x$term == "two"],
               1 / (1 + exp(-(bic[["one"]] - bic[["two"]]) / 2)))
  expect_equal(sum(x$relative_weight), 1)
  expect_equal(x$log_evidence_best_over_model, x$delta / 2)
  expect_equal(x$npar[x$term == "one"], lavaan::fitMeasures(one, "npar")[[1]])
  expect_identical(attr(x, "criterion"), "BIC")
  expect_equal(attr(x, "comparison")$N, 301)
  # Ties keep the order of fits and split the weight.
  tie <- compare_sem(list(a = one, b = one))
  expect_identical(tie$term, c("a", "b"))
  expect_equal(tie$relative_weight, c(0.5, 0.5))
  expect_output(print(x), "Criterion: BIC")
  expect_output(print(x), "approximate posterior probability")
})

test_that("every criterion is the value sem_information_criteria() gives", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  one <- lavaan::cfa("f =~ x1 + x2 + x3 + x4 + x5 + x6", data = d)
  two <- lavaan::cfa("f1 =~ x1 + x2 + x3
                      f2 =~ x4 + x5 + x6", data = d)
  for (cr in c("AIC", "BIC", "SABIC", "HBIC", "IBIC", "SPBIC")) {
    x <- compare_sem(list(one = one, two = two), criterion = tolower(cr))
    expect_identical(attr(x, "criterion"), cr)
    expect_equal(x$value[match(c("one", "two"), x$term)],
                 c(criterion_of(one, cr), criterion_of(two, cr)), label = cr)
    if (cr %in% c("AIC", "SABIC")) {
      expect_true(all(is.na(x$log_evidence_best_over_model)), label = cr)
    } else {
      expect_equal(x$log_evidence_best_over_model, x$delta / 2, label = cr)
    }
  }
  expect_output(print(compare_sem(list(one = one, two = two), criterion = "AIC")),
                "Akaike weight")
})

test_that("explicit SPBIC prior means follow the order of the fits", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  one <- lavaan::cfa("f =~ x1 + x2 + x3 + x4 + x5 + x6", data = d)
  two <- lavaan::cfa("f1 =~ x1 + x2 + x3
                      f2 =~ x4 + x5 + x6", data = d)
  pm <- list(one = lavaan::coef(one) * 0, two = lavaan::coef(two) * 0)
  x <- compare_sem(list(one = one, two = two), criterion = "SPBIC", prior_means = pm)
  expect_equal(x$value[match(c("one", "two"), x$term)],
               c(criterion_of(one, "SPBIC", prior_mean = pm$one),
                 criterion_of(two, "SPBIC", prior_mean = pm$two)))
  expect_error(compare_sem(list(one = one, two = two), criterion = "SPBIC",
                           prior_means = rev(pm)), "names and order")
  expect_error(compare_sem(list(one = one, two = two), prior_means = pm),
               "only with criterion = \"SPBIC\"")
  expect_error(compare_sem(list(one = one, two = two), criterion = "SPBIC",
                           prior_means = pm[1]), "one prior mean vector per model")
})

test_that("models conditioned on different exogenous variables are refused", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  # Two saturated models of the same three variables. Under fixed.x = TRUE
  # (lavaan::sem()'s default) the first conditions on {x1, x2} and the
  # second on {x2}, and their BIC values differ by about 930.
  a <- lavaan::sem("x3 ~ x1 + x2", data = d)
  b <- lavaan::sem("x3 ~ x1 + x2
                    x1 ~ x2", data = d)
  expect_error(compare_sem(list(a = a, b = b)), "different exogenous observed variables")
  expect_error(compare_sem(list(a = a, b = b)), "fixed.x = FALSE")
  joint <- compare_sem(list(a = lavaan::sem("x3 ~ x1 + x2", data = d, fixed.x = FALSE),
                            b = lavaan::sem("x3 ~ x1 + x2
                                             x1 ~ x2", data = d, fixed.x = FALSE)))
  expect_equal(joint$delta, c(0, 0), tolerance = 1e-6)
  expect_equal(joint$relative_weight, c(0.5, 0.5), tolerance = 1e-6)
  # With the same exogenous variables the conditional comparison is the
  # joint one.
  restricted <- lavaan::sem("x3 ~ x1 + 0*x2", data = d)
  same <- compare_sem(list(a = a, r = restricted))
  joint_same <- compare_sem(list(
    a = lavaan::sem("x3 ~ x1 + x2", data = d, fixed.x = FALSE),
    r = lavaan::sem("x3 ~ x1 + 0*x2", data = d, fixed.x = FALSE)))
  expect_equal(same$delta, joint_same$delta, tolerance = 1e-6)
  expect_identical(attr(same, "comparison")$conditioning, c("x1", "x2"))
  expect_error(compare_sem(list(
    a = lavaan::sem("x3 ~ x1 + x2", data = d, conditional.x = TRUE),
    b = lavaan::sem("x3 ~ x1 + x2
                     x1 ~ x2", data = d, conditional.x = TRUE))),
    "different exogenous observed variables")
})

test_that("weighted fits are refused", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  d$w <- rep(c(0.5, 1.5), length.out = nrow(d))
  weighted <- lavaan::cfa("f =~ x1 + x2 + x3 + x4 + x5 + x6", data = d,
                          sampling.weights = "w")
  unweighted <- lavaan::cfa("f =~ x1 + x2 + x3 + x4 + x5 + x6", data = d)
  expect_error(compare_sem(list(w = weighted, u = unweighted)),
               "fits\\$w was fitted with sampling weights")
  expect_error(compare_sem(list(u = unweighted, w = weighted)),
               "fits\\$w was fitted with sampling weights")
})

test_that("models fitted to different data or under different conventions are refused", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  m <- "f =~ x1 + x2 + x3 + x4 + x5 + x6"
  base <- lavaan::cfa(m, data = d)
  altered <- d
  altered$x1 <- altered$x1 + 0.001
  expect_error(compare_sem(list(base = base, altered = lavaan::cfa(m, data = altered))),
               "different data values")
  expect_error(compare_sem(list(base = base, part = lavaan::cfa(m, data = d[1:200, ]))),
               "different numbers of cases \\(301 and 200\\)")
  expect_error(compare_sem(list(base = base,
                                more = lavaan::cfa(paste(m, "+ x7"), data = d))),
               "different observed variables \\(none only in the first; x7 only in the second\\)")
  incomplete <- d
  incomplete$x1[1:10] <- NA
  expect_error(compare_sem(list(listwise = lavaan::cfa(m, data = incomplete),
                                fiml = lavaan::cfa(m, data = incomplete, missing = "ml"))),
               "missing = \"listwise\" and \"ml\"")
  expect_error(compare_sem(list(base = base,
                                means = lavaan::cfa(m, data = d, meanstructure = TRUE))),
               "mean structure")
  expect_error(compare_sem(list(base = base,
                                groups = lavaan::cfa(m, data = d, group = "school"))),
               "different groups")
  moments <- lavaan::cfa(m, sample.cov = stats::cov(d[, paste0("x", 1:6)]),
                         sample.nobs = 301)
  expect_error(compare_sem(list(raw = base, moments = moments)), "summary statistics")
  # The order of the variables does not matter.
  two <- lavaan::cfa("f1 =~ x1 + x2 + x3
                      f2 =~ x4 + x5 + x6", data = d)
  reversed <- lavaan::cfa("f2 =~ x4 + x5 + x6
                           f1 =~ x1 + x2 + x3", data = d)
  x <- compare_sem(list(two = two, reversed = reversed))
  expect_equal(x$relative_weight, c(0.5, 0.5), tolerance = 1e-6)
})

test_that("models fitted to the same summary statistics are compared", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  S <- stats::cov(d[, paste0("x", 1:6)])
  one <- lavaan::cfa("f =~ x1 + x2 + x3 + x4 + x5 + x6", sample.cov = S, sample.nobs = 301)
  two <- lavaan::cfa("f1 =~ x1 + x2 + x3
                      f2 =~ x4 + x5 + x6", sample.cov = S, sample.nobs = 301)
  x <- compare_sem(list(one = one, two = two))
  expect_equal(x$value[match(c("one", "two"), x$term)],
               c(lavaan::fitMeasures(one, "bic")[[1]], lavaan::fitMeasures(two, "bic")[[1]]))
  expect_false(attr(x, "comparison")$raw_data)
  other <- lavaan::cfa("f =~ x1 + x2 + x3 + x4 + x5 + x6", sample.cov = S * 1.01,
                       sample.nobs = 301)
  expect_error(compare_sem(list(one = one, other = other)), "different data values")
})

test_that("IBIC and default prior SPBIC need one identification", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  m <- "visual =~ x1 + x2 + x3
        textual =~ x4 + x5 + x6"
  marker <- lavaan::cfa(m, data = d)
  variance <- lavaan::cfa(paste(m, "\n textual =~ x1"), data = d, std.lv = TRUE)
  expect_error(compare_sem(list(marker = marker, variance = variance), criterion = "IBIC"),
               "fits\\$marker: marker indicator; fits\\$variance: variance fixed")
  expect_error(compare_sem(list(marker = marker, variance = variance), criterion = "SPBIC"),
               "more than one identification")
  expect_s3_class(compare_sem(list(marker = marker, variance = variance)), "dmar_compare_sem")
  pm <- list(lavaan::coef(marker) * 0, lavaan::coef(variance) * 0)
  expect_s3_class(compare_sem(list(marker = marker, variance = variance),
                              criterion = "SPBIC", prior_means = pm), "dmar_compare_sem")
  same <- compare_sem(list(a = lavaan::cfa(m, data = d, std.lv = TRUE), variance = variance),
                      criterion = "IBIC")
  expect_identical(unique(unlist(attr(same, "comparison")$identification)), "variance fixed")
})

test_that("multiple group models are compared by the likelihood criteria", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  m <- "visual =~ x1 + x2 + x3
        textual =~ x4 + x5 + x6"
  configural <- lavaan::cfa(m, data = d, group = "school")
  metric <- lavaan::cfa(m, data = d, group = "school", group.equal = "loadings")
  x <- compare_sem(list(configural = configural, metric = metric))
  expect_equal(x$value[match(c("configural", "metric"), x$term)],
               c(stats::BIC(configural), stats::BIC(metric)))
  expect_equal(attr(x, "comparison")$groups, 2L)
  expect_error(compare_sem(list(configural = configural, metric = metric), criterion = "IBIC"),
               "more than one group")
})

test_that("ML and its robust variants share one likelihood", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  ml <- lavaan::cfa("f =~ x1 + x2 + x3 + x4 + x5 + x6", data = d)
  mlr <- lavaan::cfa("f1 =~ x1 + x2 + x3
                      f2 =~ x4 + x5 + x6", data = d, estimator = "MLR")
  x <- compare_sem(list(ml = ml, mlr = mlr))
  expect_equal(x$value[match("mlr", x$term)], lavaan::fitMeasures(mlr, "bic")[[1]])
})

test_that("refusals name the model and the argument checks are clear", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  one <- lavaan::cfa("f =~ x1 + x2 + x3", data = d)
  expect_error(compare_sem(list(one = one, b = lm(x1 ~ x2, data = d))),
               "fits\\$b must be a lavaan model")
  expect_error(compare_sem(list(one, one)), "uniquely named")
  expect_error(compare_sem(list(one = one)), "at least two")
  expect_error(compare_sem(one), "uniquely named")
  expect_error(compare_sem(list(a = one, b = one), criterion = c("BIC", "AIC")),
               "one name")
  expect_error(compare_sem(list(a = one, b = one), criterion = "AICc"), "one name")
  ordered <- d
  for (v in c("x1", "x2", "x3")) ordered[[v]] <- cut(ordered[[v]], 3, labels = FALSE)
  wlsmv <- lavaan::cfa("f =~ x1 + x2 + x3", data = ordered, ordered = c("x1", "x2", "x3"))
  expect_error(compare_sem(list(a = wlsmv, b = wlsmv)), "fits\\$a has ordered categorical")
})

test_that("tidy() keeps every column and glance() summarizes the set", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  one <- lavaan::cfa("f =~ x1 + x2 + x3 + x4 + x5 + x6", data = d)
  two <- lavaan::cfa("f1 =~ x1 + x2 + x3
                      f2 =~ x4 + x5 + x6", data = d)
  orth <- lavaan::cfa("f1 =~ x1 + x2 + x3
                       f2 =~ x4 + x5 + x6", data = d, orthogonal = TRUE)
  x <- compare_sem(list(one = one, two = two, orth = orth), criterion = "AIC")
  td <- generics::tidy(x)
  expect_named(td, c("term", "estimate", "delta", "relative_weight",
                     "log_evidence_best_over_model", "npar", "criterion"))
  expect_equal(td$estimate, x$value)
  expect_equal(td$relative_weight, x$relative_weight)
  expect_true(all(is.na(td$log_evidence_best_over_model)))
  expect_identical(unique(td$criterion), "AIC")
  gl <- generics::glance(x)
  expect_identical(nrow(gl), 1L)
  expect_identical(gl$best_model, x$term[1])
  expect_equal(gl$best_weight, x$relative_weight[1])
  expect_equal(gl$delta_second, x$delta[2])
  expect_identical(gl$n_models, 3L)
  expect_identical(gl$criterion, "AIC")
})

test_that("models that each mix the two identifications are refused with the rule, not a false difference", {
  skip_if_not_installed("lavaan")
  d <- lavaan::HolzingerSwineford1939
  a <- "visual =~ x1 + x2 + x3
        textual =~ NA*x4 + x5 + x6
        textual ~~ 1*textual"
  fa <- lavaan::cfa(a, data = d)
  fb <- lavaan::cfa(paste(a, "\n visual =~ x5"), data = d)
  for (crit in c("IBIC", "SPBIC")) {
    err <- tryCatch(compare_sem(list(a = fa, b = fb), criterion = crit),
                    error = conditionMessage)
    expect_match(err, "every latent variable of every model")
    expect_match(err, "more than one identification")
    expect_false(grepl("identify them differently", err))
  }
})

test_that("a two-stage fit and a robust two-stage fit are compared as one method", {
  skip_if_not_installed("lavaan")
  skip_on_cran()
  # The two settings fit the same EM moments by the same likelihood and
  # differ only in their standard errors and test, which no criterion uses.
  d <- lavaan::HolzingerSwineford1939[paste0("x", 1:9)]
  set.seed(1)
  d$x1[sample(301, 60)] <- NA
  three <- "vis =~ x1 + x2 + x3\n text =~ x4 + x5 + x6\n speed =~ x7 + x8 + x9"
  one <- "g =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"
  t3 <- lavaan::cfa(three, d, missing = "two.stage")
  t1 <- lavaan::cfa(one, d, missing = "two.stage")
  r1 <- lavaan::cfa(one, d, missing = "robust.two.stage")
  expect_equal(as.numeric(lavaan::logLik(t1)), as.numeric(lavaan::logLik(r1)))
  mixed <- compare_sem(list(three = t3, one = r1))
  same <- compare_sem(list(three = t3, one = t1))
  expect_equal(mixed$value, same$value)
  expect_equal(mixed$delta, same$delta)
  expect_identical(attr(mixed, "comparison")$missing, c("two.stage", "robust.two.stage"))
  expect_match(attr(mixed, "comparison")$footer[2],
               "missing = \"two.stage\" and \"robust.two.stage\"", fixed = TRUE)
  expect_true(all(is.na(mixed$log_evidence_best_over_model)))
  # FIML is still a different likelihood.
  fiml <- lavaan::cfa(one, d, missing = "ml")
  expect_error(compare_sem(list(three = t3, one = fiml)),
               "handle missing data differently")
})
