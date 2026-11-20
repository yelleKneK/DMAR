## The numerics of two-stage maximum likelihood (R/two_stage_internals.R):
## the observed information of the saturated model under missing data, the
## stage one moments, the residual based statistic of Savalei and Bentler
## (2009, Eqs. 8 and 14), the corrected indices of Zhang and Savalei (2023),
## and the rows lavaan's two-stage N leaves out.
##
## Oracles. The residual statistic was reproduced against semTools 0.5-9's
## twostage() (anova(), "Browne (1984) residual-based test statistic"),
## whose printed values are hard-coded below with their provenance rather
## than computed, since semTools is not a dependency. For models with
## equality constraints semTools differentiates with respect to the
## unconstrained parameters and reports a statistic on the wrong degrees of
## freedom (81.55 for the metric model below); the oracle there is lavaan's
## own browne.residual.nt statistic on complete data, which the two-stage
## statistic equals when nothing is missing. All values were computed with
## lavaan 0.7-3.

# Fixture O1: nine Holzinger and Swineford tests (they correspond exactly
# to lavaan's x1--x9) with 90 values of two of them deleted completely at
# random. The three speed tests are divided by 10, which leaves every test
# statistic and fit index unchanged and spares lavaan's warning about
# variances of very different sizes.
.ts_o1 <- function() {
  hs <- holzinger_swineford[c(
    "t1_visual_perception", "t2_cubes", "t4_lozenges",
    "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning",
    "t10_addition", "t12_counting_groups_of_dots",
    "t13_straight_and_curved_capitals", "school", "age_years")]
  set.seed(1)
  hs$t7_sentence[sample(301, 90)] <- NA
  hs$t13_straight_and_curved_capitals[sample(301, 90)] <- NA
  for (v in c("t10_addition", "t12_counting_groups_of_dots",
              "t13_straight_and_curved_capitals")) hs[[v]] <- hs[[v]] / 10
  hs
}

.ts_o1_model <- "
  spatial =~ t1_visual_perception + t2_cubes + t4_lozenges
  verbal  =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
  speed   =~ t10_addition + t12_counting_groups_of_dots +
             t13_straight_and_curved_capitals"

.ts_skip <- function() {
  skip_if_not_installed("lavaan")
  skip_if(utils::packageVersion("lavaan") < "0.7-3",
          "the two-stage oracles were computed with lavaan 0.7-3")
}

test_that("the duplication matrix maps vech to vec in lavaan's order", {
  set.seed(113)
  S <- crossprod(matrix(rnorm(16), 4, 4))
  vech <- S[lower.tri(S, diag = TRUE)]
  expect_equal(as.vector(.ts_dup_matrix(4) %*% vech), as.vector(S))
  expect_equal(dim(.ts_dup_matrix(3)), c(9L, 6L))
})

test_that("the saturated information is lavaan's FIML h1 observed information", {
  .ts_skip()
  d <- .ts_o1()
  fml <- lavaan::cfa(.ts_o1_model, data = d, missing = "ml",
                     h1.information = "unstructured",
                     observed.information = "h1", information = "observed")
  ov <- lavaan::lavNames(fml, "ov")
  ss <- lavaan::lavInspect(fml, "sampstat")
  A <- .ts_sat_information(as.matrix(d[ov]), as.numeric(ss$mean[ov]),
                           unclass(ss$cov)[ov, ov])
  A_lav <- unclass(lavaan::lavInspect(fml, "h1.information.observed"))
  # .ts_sat_information() sums over the cases; lavaan reports per case.
  expect_lt(max(abs(A / nrow(d) - A_lav)) / max(abs(A_lav)), 1e-12)
})

test_that("the stage one moments are lavaan's, with and without auxiliaries", {
  .ts_skip()
  d <- .ts_o1()
  f <- lavaan::cfa(.ts_o1_model, data = d, missing = "two.stage")
  ov <- lavaan::lavNames(f, "ov")
  ss <- lavaan::lavInspect(f, "sampstat")
  em <- .ts_em_moments(d[ov])
  expect_equal(em$cov, unname(unclass(ss$cov)[ov, ov]), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(em$mean, as.numeric(ss$mean[ov]), tolerance = 1e-10)
  # With an auxiliary the item block of the augmented EM estimate is the
  # moments lavaan's second stage used.
  fa <- lavaan::cfa(.ts_o1_model, data = d, missing = "two.stage",
                    aux = "age_years")
  ssa <- lavaan::lavInspect(fa, "sampstat")
  ema <- .ts_em_moments(d[c(ov, "age_years")])
  expect_equal(ema$cov[1:9, 1:9], unname(unclass(ssa$cov)[ov, ov]),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_identical(.ts_fit_aux(fa)$names, "age_years")
  expect_null(.ts_fit_aux(f))
})

test_that("the two-stage sandwich with this information is lavaan's vcov", {
  .ts_skip()
  # This pins the convention the residual statistic shares with lavaan's
  # standard errors: acov = (D'HD)^-1 D'H A^-1 H D (D'HD)^-1 (Savalei &
  # Bentler, 2009, Eqs. 6 and 13), with H the complete data information at
  # the EM moments (lavaan's wls.v) and A the information computed here,
  # summed over the cases; with an auxiliary, its Schur complement.
  d <- .ts_o1()
  sandwich <- function(fit, A) {
    D <- unclass(lavaan::lavInspect(fit, "delta"))
    H <- unclass(lavaan::lavInspect(fit, "wls.v"))
    B <- solve(crossprod(D, H %*% D))
    B %*% crossprod(D, H) %*% solve(A, H %*% D) %*% B
  }
  f <- lavaan::cfa(.ts_o1_model, data = d, missing = "two.stage")
  ov <- lavaan::lavNames(f, "ov")
  ss <- lavaan::lavInspect(f, "sampstat")
  A <- .ts_sat_information(as.matrix(d[ov]), as.numeric(ss$mean[ov]),
                           unclass(ss$cov)[ov, ov])
  V <- unclass(lavaan::vcov(f))
  expect_lt(max(abs(sandwich(f, A) - V)) / max(abs(V)), 1e-8)

  fa <- lavaan::cfa(.ts_o1_model, data = d, missing = "two.stage",
                    aux = "age_years")
  Y <- as.matrix(d[c(ov, "age_years")])
  em <- .ts_em_moments(Y)
  Aa <- .ts_schur_items(.ts_sat_information(Y, em$mean, em$cov),
                        p_all = 10L, p_items = 9L)
  Va <- unclass(lavaan::vcov(fa))
  expect_lt(max(abs(sandwich(fa, Aa) - Va)) / max(abs(Va)), 1e-8)
})

test_that("the residual statistic reproduces semTools on O1, without and with an auxiliary", {
  .ts_skip()
  d <- .ts_o1()
  f <- lavaan::cfa(.ts_o1_model, data = d, missing = "two.stage")
  r <- .ts_residual_test(f)
  # semTools 0.5-9: 65.293 (printed); the unrounded value of the same
  # computation on lavaan's moments is 65.2930632982, p = 1.10240e-05.
  expect_equal(r$statistic, 65.2930632982, tolerance = 1e-6)
  expect_equal(r$df, 24)
  expect_equal(r$p_value, 1.10240e-05, tolerance = 1e-4)

  fa <- lavaan::cfa(.ts_o1_model, data = d, missing = "two.stage",
                    aux = "age_years")
  ra <- .ts_residual_test(fa)
  # semTools 0.5-9: 64.6563565. The value built on lavInspect(fit, "gamma"),
  # which with an auxiliary is not the matrix lavaan's own standard errors
  # use, would be 64.2028232808.
  expect_equal(ra$statistic, 64.6563565, tolerance = 1e-6)
  expect_gt(abs(ra$statistic - 64.2028232808), 0.4)
})

test_that("multiple groups and equality constraints are projected correctly", {
  .ts_skip()
  d <- .ts_o1()
  conf <- lavaan::cfa(.ts_o1_model, data = d, missing = "two.stage",
                      group = "school")
  metric <- lavaan::cfa(.ts_o1_model, data = d, missing = "two.stage",
                        group = "school", group.equal = "loadings")
  # semTools 0.5-9 prints 81.41336 for the configural model.
  expect_equal(.ts_residual_test(conf)$statistic, 81.4134204798,
               tolerance = 1e-6)
  # For the metric model semTools prints 81.55182 (on 54 df), the value of
  # a projection on the unconstrained parameters; respecting the six
  # equality constraints gives 89.8059318463, 8.39 above the configural
  # model, beside a scaled difference of 7.92.
  rm <- .ts_residual_test(metric)
  expect_equal(rm$statistic, 89.8059318463, tolerance = 1e-6)
  expect_equal(rm$df, 54)
  expect_gt(abs(rm$statistic - 81.5517846845), 8)

  # On complete data the two-stage statistic is Browne's normal theory
  # residual statistic, which lavaan computes itself (browne.residual.nt),
  # constraints included.
  dc <- d[stats::complete.cases(d), ]
  browne <- function(...) {
    fit <- lavaan::cfa(.ts_o1_model, data = dc, group = "school",
                       meanstructure = TRUE,
                       test = c("standard", "browne.residual.nt"), ...)
    lavaan::lavInspect(fit, "test")$browne.residual.nt$stat
  }
  ts <- function(...) {
    .ts_residual_test(lavaan::cfa(.ts_o1_model, data = dc,
                                  group = "school", missing = "two.stage",
                                  ...))$statistic
  }
  expect_equal(ts(), browne(), tolerance = 1e-6)
  expect_equal(ts(group.equal = "loadings"),
               browne(group.equal = "loadings"), tolerance = 1e-6)
  # A general equality constraint, through the constraint Jacobian.
  m2 <- "f =~ a*t6_paragraph_comprehension + b*t7_sentence +
               t9_word_meaning + t1_visual_perception
         b == 1.5*a"
  fb <- lavaan::cfa(m2, data = dc, std.lv = TRUE, meanstructure = TRUE,
                    test = c("standard", "browne.residual.nt"))
  ft <- lavaan::cfa(m2, data = dc, std.lv = TRUE, missing = "two.stage")
  expect_equal(.ts_residual_test(ft)$statistic,
               lavaan::lavInspect(fb, "test")$browne.residual.nt$stat,
               tolerance = 1e-6)
  expect_equal(.ts_residual_test(ft)$df, 3)
})

test_that("lavaan's robust indices are TS V2 after robust two-stage estimation and its normal theory analog after two-stage estimation", {
  .ts_skip()
  # Zhang and Savalei (2023, Table 2, p. 268): RMSEA and CFI from the naive
  # statistic T less an estimate c = tr(U Gamma) of its expectation, with U
  # the residual weight matrix of the complete data information at the EM
  # estimates (W_c, the saturated estimates: TS V2) and Gamma the
  # asymptotic covariance matrix of the stage one moments. TS V2 takes
  # Gamma = n A^-1 V A^-1, the sandwich of the observed (A) and first order
  # (V) information summed over the cases; under normality Gamma = n A^-1.
  # Both are built here from the data, independently of lavaan's scaling
  # factors, and the baseline model (free means and variances) gets its own
  # U. lavaan's robust.two.stage indices are TS V2; its two.stage indices
  # are the normal theory analog, which differs from TS V2.
  d <- .ts_o1()
  f <- lavaan::cfa(.ts_o1_model, data = d, missing = "two.stage")
  fr <- lavaan::cfa(.ts_o1_model, data = d, missing = "robust.two.stage")
  ov <- lavaan::lavNames(f, "ov")
  X <- as.matrix(d[ov])
  n <- nrow(X)
  p <- ncol(X)
  ss <- lavaan::lavInspect(f, "sampstat")
  mu <- as.numeric(ss$mean[ov])
  S <- unclass(ss$cov)[ov, ov]
  A <- .ts_sat_information(X, mu, S)
  Dp <- .ts_dup_matrix(p)
  # The first order information: the sum over the cases of the outer
  # product of each case's score for the means and vech(Sigma).
  V <- matrix(0, nrow(A), ncol(A))
  for (i in seq_len(n)) {
    o <- which(!is.na(X[i, ]))
    Si <- solve(S[o, o, drop = FALSE])
    r <- X[i, o] - mu[o]
    g_mu <- numeric(p)
    g_mu[o] <- Si %*% r
    G <- matrix(0, p, p)
    G[o, o] <- (Si %*% tcrossprod(r) %*% Si - Si) / 2
    V <- V + tcrossprod(c(g_mu, crossprod(Dp, as.vector(G))))
  }
  Si <- solve(S)
  k <- p * (p + 1) / 2
  Wc <- rbind(cbind(Si, matrix(0, p, k)),
              cbind(matrix(0, k, p),
                    crossprod(Dp, kronecker(Si, Si) %*% Dp) / 2))
  U <- function(D) Wc - Wc %*% D %*% solve(crossprod(D, Wc %*% D),
                                           crossprod(D, Wc))
  D_model <- unclass(lavaan::lavInspect(f, "delta"))
  vi <- which(lower.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  D_base <- matrix(0, p + k, 2 * p)
  D_base[cbind(seq_len(p), seq_len(p))] <- 1
  D_base[cbind(p + which(vi[, 1] == vi[, 2]), p + seq_len(p))] <- 1
  A_inv <- solve(A)
  gamma_nt <- n * A_inv
  gamma_v2 <- n * A_inv %*% V %*% A_inv
  fm <- lavaan::fitMeasures(f)
  fmr <- lavaan::fitMeasures(fr)
  indices <- function(gamma) {
    c_m <- sum(diag(U(D_model) %*% gamma))
    c_b <- sum(diag(U(D_base) %*% gamma))
    t_m <- fm[["chisq"]] - c_m
    c(rmsea = sqrt(max(t_m, 0) / (n * fm[["df"]])),
      cfi = 1 - max(t_m, 0) / max(fm[["baseline.chisq"]] - c_b, t_m, 0),
      c = c_m)
  }
  nt <- indices(gamma_nt)
  v2 <- indices(gamma_v2)
  # The naive statistic and the estimates are the same in both fits.
  expect_equal(fmr[["chisq"]], fm[["chisq"]])
  expect_equal(nt[["c"]], fm[["df"]] * fm[["chisq.scaling.factor"]],
               tolerance = 1e-8)
  expect_equal(v2[["c"]], fmr[["df"]] * fmr[["chisq.scaling.factor"]],
               tolerance = 1e-8)
  expect_equal(fm[["rmsea.robust"]], nt[["rmsea"]], tolerance = 1e-8)
  expect_equal(fm[["cfi.robust"]], nt[["cfi"]], tolerance = 1e-8)
  expect_equal(fmr[["rmsea.robust"]], v2[["rmsea"]], tolerance = 1e-8)
  expect_equal(fmr[["cfi.robust"]], v2[["cfi"]], tolerance = 1e-8)
  expect_equal(nt[["rmsea"]], 0.0869182751, tolerance = 1e-6)
  expect_equal(v2[["rmsea"]], 0.0859051291, tolerance = 1e-6)
  expect_equal(nt[["cfi"]], 0.9374251254, tolerance = 1e-6)
  expect_equal(v2[["cfi"]], 0.9387864529, tolerance = 1e-6)
  expect_gt(abs(nt[["rmsea"]] - v2[["rmsea"]]), 1e-3)
  # And the labels say which is which.
  expect_identical(.ts_indices_label("robust.two.stage"),
                   "the TS V2 indices of Zhang and Savalei (2023)")
  expect_match(.ts_indices_label("two.stage"),
               "^the normal theory analog of the TS V2 indices")
})

test_that("the residual statistic reads equality constraints kept under ceq.simple = TRUE", {
  .ts_skip()
  skip_on_cran()  # seven fits; the == constraints are projected on CRAN above
  # Under ceq.simple = TRUE lavaan writes no == rows for simple equalities
  # and keeps the map from the estimated to the unconstrained parameters
  # instead; the statistic is the one of ceq.simple = FALSE, the same
  # model, which was checked against semTools and complete data above.
  d <- .ts_o1()
  metric <- lavaan::cfa(.ts_o1_model, data = d, missing = "two.stage",
                        group = "school", group.equal = "loadings",
                        ceq.simple = TRUE)
  expect_true(metric@Model@ceq.simple.only)
  expect_false(any(lavaan::parTable(metric)$op == "=="))
  rm <- .ts_residual_test(metric)
  expect_equal(rm$statistic, 89.8059318463, tolerance = 1e-6)
  expect_equal(rm$df, 54)
  # On complete data it is lavaan's browne.residual.nt under the same
  # setting, and a label shared within one group is read the same way.
  dc <- d[stats::complete.cases(d), ]
  fb <- lavaan::cfa(.ts_o1_model, data = dc, group = "school",
                    group.equal = "loadings", ceq.simple = TRUE,
                    meanstructure = TRUE,
                    test = c("standard", "browne.residual.nt"))
  ft <- lavaan::cfa(.ts_o1_model, data = dc, group = "school",
                    group.equal = "loadings", ceq.simple = TRUE,
                    missing = "two.stage")
  expect_equal(.ts_residual_test(ft)$statistic,
               lavaan::lavInspect(fb, "test")$browne.residual.nt$stat,
               tolerance = 1e-6)
  m1 <- "f =~ a*t6_paragraph_comprehension + a*t7_sentence +
               t9_word_meaning + t1_visual_perception"
  shared <- function(ceq) {
    .ts_residual_test(lavaan::cfa(m1, data = d, std.lv = TRUE,
                                  missing = "two.stage",
                                  ceq.simple = ceq))$statistic
  }
  expect_equal(shared(TRUE), shared(FALSE), tolerance = 1e-6)
  # With a general constraint lavaan ignores ceq.simple and writes every
  # constraint as an == row, which the constraint Jacobian handles.
  m2 <- "f =~ a*t6_paragraph_comprehension + a*t7_sentence +
               b*t9_word_meaning + c*t1_visual_perception
         c == 1.5*b"
  general <- function(ceq) {
    .ts_residual_test(lavaan::cfa(m2, data = d, std.lv = TRUE,
                                  missing = "two.stage",
                                  ceq.simple = ceq))$statistic
  }
  expect_equal(general(TRUE), general(FALSE), tolerance = 1e-8)
})

test_that("rows observed only on an auxiliary variable are what lavaan's two-stage N leaves out", {
  .ts_skip()
  # Six tests with the whole block missing for 77 of 301 pupils, at random
  # given general information (the auxiliary). lavaan's two-stage fit keeps
  # those rows in its data and in its first stage but leaves them out of
  # its N, so each standard error is the two-stage sandwich of every row
  # (Savalei & Bentler, 2009, Eqs. 6 and 13) times sqrt(301/224).
  hs <- holzinger_swineford
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  set.seed(113)
  pr <- stats::plogis(-1.5 + 1.2 *
                        as.numeric(scale(hs$t5_general_information)))
  hs[sample(301, 77, prob = pr), items] <- NA
  m <- "spatial =~ t1_visual_perception + t2_cubes + t4_lozenges
        verbal  =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning"
  fit <- suppressWarnings(lavaan::cfa(m, data = hs, missing = "two.stage",
                                      aux = "t5_general_information",
                                      meanstructure = TRUE))
  expect_identical(nrow(lavaan::lavInspect(fit, "data")), 301L)
  expect_identical(lavaan::lavInspect(fit, "nobs"), 224L)
  Y <- as.matrix(hs[c(items, "t5_general_information")])
  em <- .ts_em_moments(Y)
  A <- .ts_schur_items(.ts_sat_information(Y, em$mean, em$cov),
                       p_all = 7L, p_items = 6L)
  D <- unclass(lavaan::lavInspect(fit, "delta"))
  H <- unclass(lavaan::lavInspect(fit, "wls.v"))
  B <- solve(crossprod(D, H %*% D))
  se_all <- sqrt(diag(B %*% crossprod(D, H) %*% solve(A, H %*% D) %*% B))
  ratio <- sqrt(diag(lavaan::vcov(fit))) / se_all
  expect_equal(unname(ratio), rep(sqrt(301 / 224), length(ratio)),
               tolerance = 1e-6)

  # .ts_rows_read() counts the cases the first stage read, and the fit
  # computed again with that N (.ts_all_rows()) has the sandwich over all
  # of them as its covariance matrix, with the estimates unchanged.
  expect_identical(.ts_rows_read(fit), list(read = 301L, nobs = 224L))
  expect_identical(.ts_aux_only_cases(fit), 77L)
  all_rows <- .ts_all_rows(fit)
  expect_true(all_rows$corrected)
  expect_identical(all_rows$extra, 77L)
  expect_identical(lavaan::lavInspect(all_rows$fit, "nobs"), 301L)
  expect_identical(.ts_aux_only_cases(all_rows$fit), 0L)
  expect_identical(lavaan::coef(all_rows$fit), lavaan::coef(fit))
  expect_false(all_rows$reestimated)
  expect_false(.ts_groups_linked(fit))
  expect_equal(unname(sqrt(diag(lavaan::vcov(all_rows$fit)))), unname(se_all),
               tolerance = 1e-8)
  # The fit keeps the iterations of its estimation and its optimizer, and
  # lavaan's function value goes with the statistic (2 N fmin).
  fm <- lavaan::fitMeasures(all_rows$fit, c("chisq", "fmin", "ntotal"))
  expect_equal(fm[["chisq"]], 2 * 301 * fm[["fmin"]], tolerance = 1e-10)
  expect_identical(lavaan::lavInspect(all_rows$fit, "iterations"),
                   lavaan::lavInspect(fit, "iterations"))
  expect_identical(lavaan::lavInspect(all_rows$fit, "options")$optim.method,
                   lavaan::lavInspect(fit, "options")$optim.method)
  # The data slot counts the cases read as well: lavaan's summary() calls
  # all 301 used, and lavTestLRT() warns that a difference test with a fit
  # lavaan made of the same data on its 224 is invalid (through the first
  # versions of 1.1.0 the summary said 224 were used, and lavTestLRT()
  # reported 35.256 for that pair without a word). The rows of the data,
  # and lavPredict(), are unchanged.
  expect_identical(all_rows$fit@Data@nobs, list(301L))
  out <- utils::capture.output(lavaan::summary(all_rows$fit))
  expect_match(out, "^ +Number of observations +301$", all = FALSE)
  expect_false(any(grepl("Used", out, fixed = TRUE)))
  m1 <- paste("g =~", paste(items, collapse = " + "))
  one <- suppressWarnings(lavaan::cfa(m1, data = hs, missing = "two.stage",
                                      aux = "t5_general_information",
                                      meanstructure = TRUE))
  expect_warning(lavaan::lavTestLRT(one, all_rows$fit),
                 "numbers of observations differ")
  expect_no_warning(lavaan::lavTestLRT(.ts_all_rows(one)$fit, all_rows$fit))
  expect_equal(lavaan::lavPredict(all_rows$fit), lavaan::lavPredict(fit),
               tolerance = 1e-12)
  # A fit without such cases, or without an auxiliary, is returned as it is.
  plain <- suppressWarnings(lavaan::cfa(m, data = holzinger_swineford,
                                        missing = "two.stage",
                                        aux = "t5_general_information"))
  expect_identical(.ts_all_rows(plain)$fit, plain)
  expect_true(is.na(.ts_all_rows(plain)$corrected))
  expect_null(.ts_rows_read(lavaan::cfa(m, data = holzinger_swineford)))
  expect_identical(.ts_dots_aux(list(aux = "z", se = "x")), "z")
  expect_identical(.ts_dots_aux(list(au = c("z", "w"))), c("z", "w"))
  expect_null(.ts_dots_aux(list(se = "x")))
  expect_null(.ts_dots_aux(list()))
})

test_that("in several groups held together by constraints the fit computed again is the two-stage fit with every case read", {
  .ts_skip()
  skip_on_cran()  # the first order information is summed case by case
  # The six tests by school with the whole block missing for 77 pupils at
  # random given general information, and 40 sentence completion scores
  # missing besides. The groups' counts grow in different ratios (130 to
  # 156 and 94 to 145), so in a model with constraints across groups the
  # weights N_g / N of the discrepancy function move, and with them the
  # estimates, the information, and the trace behind the scaling constant;
  # no single factor corrects lavaan's values. Three oracles check the fit
  # computed again. (1) The estimates, the naive statistic, and the log
  # likelihood are those of an ordinary maximum likelihood fit of the
  # stage one moments with the counts of the cases read, the two-stage
  # estimator with every case read; lavaan's own estimates are not
  # (through the first versions of 1.1.0 they were kept, and the values
  # below were computed at them with the new weights). (2) semTools 0.5-9's
  # twostage(), whose first stage counts every case, gives the same naive
  # statistic and log likelihood to the tolerance of its own EM step. (3)
  # The standard errors and the trace are built from the data over every
  # case read, at the estimates of the fit: the asymptotic covariance
  # matrix Gamma_g of the stage one moments, per case, as the inverse of
  # the observed information of the saturated model with the auxiliary (its
  # Schur complement on the items; normal theory) or as the item block of
  # the sandwich of the observed and first order information (robust), and
  # then the two-stage sandwich and tr(U Gamma) with the counts of the
  # cases read (Savalei & Bentler, 2009, Eqs. 6 and 13; Satorra & Bentler,
  # 1994).
  hs <- holzinger_swineford
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
  set.seed(113)
  pr <- stats::plogis(-1.5 + 1.2 *
                        as.numeric(scale(hs$t5_general_information)))
  hs[sample(301, 77, prob = pr), items] <- NA
  set.seed(2)
  hs$t7_sentence[sample(301, 40)] <- NA
  m <- "spatial =~ t1_visual_perception + t2_cubes + t4_lozenges
        verbal  =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning"
  first_order <- function(X, mu, S) {
    p <- ncol(X)
    Dp <- .ts_dup_matrix(p)
    V <- 0
    for (i in seq_len(nrow(X))) {
      o <- which(!is.na(X[i, ]))
      Si <- solve(S[o, o, drop = FALSE])
      r <- X[i, o] - mu[o]
      g_mu <- numeric(p)
      g_mu[o] <- Si %*% r
      G <- matrix(0, p, p)
      G[o, o] <- (Si %*% tcrossprod(r) %*% Si - Si) / 2
      V <- V + tcrossprod(c(g_mu, crossprod(Dp, as.vector(G))))
    }
    V
  }
  oracle <- function(fit, robust) {
    X <- lavaan::lavInspect(fit, "data")
    aux <- .ts_fit_aux(fit)
    H <- lapply(lavaan::lavInspect(fit, "wls.v"), unclass)
    D <- lapply(lavaan::lavInspect(fit, "delta"), unclass)
    J <- unclass(lavaan::lavInspect(fit, "con.jac"))
    q <- qr(t(J))
    K <- qr.Q(q, complete = TRUE)[, -seq_len(q$rank), drop = FALSE]
    n <- numeric(2)
    gam <- vector("list", 2)
    for (g in 1:2) {
      Y <- cbind(as.matrix(X[[g]]), aux$data[[g]])
      Y <- Y[rowSums(!is.na(Y)) > 0, , drop = FALSE]
      n[g] <- nrow(Y)
      em <- .ts_em_moments(Y)
      A <- .ts_sat_information(Y, em$mean, em$cov)
      vi <- which(lower.tri(diag(7), diag = TRUE), arr.ind = TRUE)
      idx <- c(1:6, 7 + which(vi[, 1] <= 6 & vi[, 2] <= 6))
      gam[[g]] <- if (robust) {
        Ai <- solve(A)
        n[g] * (Ai %*% first_order(Y, em$mean, em$cov) %*% Ai)[idx, idx]
      } else {
        n[g] * solve(.ts_schur_items(A, 7L, 6L))
      }
    }
    fg <- n / sum(n)
    E <- Reduce(`+`, lapply(1:2, function(g) {
      fg[g] * crossprod(D[[g]], H[[g]] %*% D[[g]])
    }))
    M <- Reduce(`+`, lapply(1:2, function(g) {
      fg[g] * crossprod(D[[g]], H[[g]] %*% gam[[g]] %*% H[[g]] %*% D[[g]])
    }))
    E_inv <- K %*% solve(crossprod(K, E %*% K), t(K))
    list(n = n,
         se = sqrt(diag(E_inv %*% M %*% E_inv) / sum(n)),
         trace = sum(vapply(1:2, function(g) {
           sum(diag(H[[g]] %*% gam[[g]]))
         }, numeric(1))) - sum(diag(E_inv %*% M)))
  }
  for (missing in c("two.stage", "robust.two.stage")) {
    fit <- suppressWarnings(lavaan::cfa(
      m, data = hs, std.lv = TRUE, missing = missing,
      aux = "t5_general_information", meanstructure = TRUE,
      group = "school", group.equal = c("loadings", "intercepts")))
    expect_identical(.ts_rows_read(fit)$nobs, c(130L, 94L))
    expect_true(.ts_groups_linked(fit))
    a <- .ts_all_rows(fit)
    expect_true(a$corrected)
    expect_true(a$reestimated)
    o <- oracle(a$fit, robust = missing == "robust.two.stage")
    expect_identical(as.numeric(.ts_rows_read(fit)$read), o$n)
    fm <- lavaan::fitMeasures(a$fit)
    # (1) The two-stage estimator with every case read.
    ss <- lavaan::lavInspect(fit, "sampstat")
    ml <- lavaan::cfa(m, sample.cov = lapply(ss, function(x) unclass(x$cov)),
                      sample.mean = lapply(ss, function(x) unclass(x$mean)),
                      sample.nobs = o$n, sample.cov.rescale = FALSE,
                      likelihood = "normal", std.lv = TRUE,
                      meanstructure = TRUE,
                      group.equal = c("loadings", "intercepts"))
    expect_equal(lavaan::coef(a$fit), lavaan::coef(ml), tolerance = 1e-5)
    expect_gt(max(abs(lavaan::coef(fit) - lavaan::coef(a$fit))), 0.5)
    expect_equal(fm[["chisq"]], lavaan::fitMeasures(ml, "chisq")[[1]],
                 tolerance = 1e-9)
    expect_equal(fm[["logl"]], lavaan::fitMeasures(ml, "logl")[[1]],
                 tolerance = 1e-11)
    # lavaan's function value weights the groups by the cases read too
    # (through the first versions of 1.1.0 the fit kept the record of the
    # fit lavaan made, whose function value weights them by lavaan's
    # counts, so 2 N fmin was not the statistic).
    expect_equal(fm[["chisq"]], 2 * sum(o$n) * fm[["fmin"]], tolerance = 1e-10)
    # (2) semTools 0.5-9, twostage(model, data, group = "school", aux =
    # "t5_general_information", fun = "cfa", std.lv = TRUE, group.equal =
    # c("loadings", "intercepts"), meanstructure = TRUE), 2026-11-02: its
    # target fit counts 156 and 145 cases, its naive statistic is
    # 63.2332982986 and its log likelihood -5377.9382727567. At lavaan's
    # estimates with the new weights they would be 63.585 and -5378.114.
    expect_equal(fm[["chisq"]], 63.2332982986, tolerance = 1e-5)
    expect_equal(fm[["logl"]], -5377.9382727567, tolerance = 1e-7)
    # (3) The sandwich and the trace at the fit's estimates.
    expect_equal(unname(sqrt(diag(lavaan::vcov(a$fit)))), o$se,
                 tolerance = 1e-9)
    expect_equal(fm[["chisq.scaling.factor"]] * fm[["df"]], o$trace,
                 tolerance = 1e-9)
    # lavaan's own scaling constant weights the groups by their N.
    expect_gt(abs(lavaan::fitMeasures(fit, "chisq.scaling.factor")[[1]] *
                    fm[["df"]] - o$trace), 0.1)
    # The iterations of lavaan's estimation and of the new one.
    expect_gt(lavaan::lavInspect(a$fit, "iterations"),
              lavaan::lavInspect(fit, "iterations"))
  }
  # Without a constraint across groups the estimates do not depend on the
  # counts and are lavaan's.
  free <- suppressWarnings(lavaan::cfa(
    m, data = hs, std.lv = TRUE, missing = "two.stage",
    aux = "t5_general_information", meanstructure = TRUE, group = "school"))
  af <- .ts_all_rows(free)
  expect_true(af$corrected)
  expect_false(af$reestimated)
  expect_identical(lavaan::coef(af$fit), lavaan::coef(free))
})

test_that("the groups are held together by a shared parameter or a constraint across them", {
  .ts_skip()
  hs <- holzinger_swineford
  m <- "spatial =~ t1_visual_perception + t2_cubes + t4_lozenges"
  fit <- function(model, ...) {
    suppressWarnings(lavaan::cfa(model, data = hs, group = "school",
                                 meanstructure = TRUE, ...))
  }
  expect_false(.ts_groups_linked(fit(m)))
  expect_false(.ts_groups_linked(lavaan::cfa(m, data = hs)))
  expect_true(.ts_groups_linked(fit(m, group.equal = "loadings")))
  expect_true(.ts_groups_linked(fit(m, group.equal = "intercepts",
                                    ceq.simple = TRUE)))
  # A label shared by the groups holds them together (the first loading is
  # the marker, fixed at 1); labels shared within each group do not.
  expect_true(.ts_groups_linked(fit(
    "spatial =~ t1_visual_perception + c(a, a)*t2_cubes + t4_lozenges")))
  expect_false(.ts_groups_linked(fit(paste(
    "spatial =~ t1_visual_perception + c(a1, a2)*t2_cubes +",
    "c(a1, a2)*t4_lozenges"))))
  # A constraint on a defined parameter is taken to hold them together.
  expect_true(.ts_groups_linked(fit(paste(
    "spatial =~ t1_visual_perception + c(a1, a2)*t2_cubes + t4_lozenges",
    "d := a1 - a2", "d == 0", sep = "\n"))))
})

test_that("the residual statistic refuses what it cannot compute", {
  .ts_skip()
  d <- .ts_o1()
  expect_error(.ts_residual_test(lavaan::cfa(.ts_o1_model, data = d,
                                             missing = "ml")),
               "defined for a two-stage fit")
  d$age_c <- d$age_years - mean(d$age_years)
  fx <- lavaan::sem(paste(.ts_o1_model, "\n verbal ~ age_c"), data = d,
                    missing = "two.stage", fixed.x = TRUE)
  expect_error(.ts_residual_test(fx), "fixed.x = FALSE")
  # The refit leaves the estimates as they are only when the exogenous
  # variables are complete: with ages missing, lavaan deletes those cases
  # under fixed.x = TRUE and keeps them under fixed.x = FALSE. Through the
  # first versions of 1.1.0 the message said the estimates stay whatever
  # the data.
  v <- holzinger_swineford[c("t6_paragraph_comprehension", "t7_sentence",
                             "t9_word_meaning", "age_years")]
  set.seed(1)
  v$t7_sentence[sample(301, 90)] <- NA
  v$age_years[c(5, 17, 33, 60, 91, 120, 200, 250)] <- NA
  mv <- "verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
         verbal ~ age_years"
  a <- suppressWarnings(lavaan::sem(mv, data = v, missing = "two.stage"))
  b <- lavaan::sem(mv, data = v, missing = "two.stage", fixed.x = FALSE)
  expect_identical(c(lavaan::lavInspect(a, "nobs"),
                     lavaan::lavInspect(b, "nobs")), c(293L, 301L))
  expect_gt(abs(lavaan::coef(a)[["verbal~age_years"]] -
                  lavaan::coef(b)[["verbal~age_years"]]), 0.005)
  e <- tryCatch(.ts_residual_test(a), error = conditionMessage)
  expect_match(e, paste0("When those variables are complete the refit ",
                         "leaves the estimates as they are; when one has ",
                         "missing values, lavaan has deleted the cases ",
                         "missing it under fixed.x = TRUE, and the refit ",
                         "keeps them, so its estimates differ."),
               fixed = TRUE)
  expect_match(e, "under fixed.x = TRUE it deletes the same cases",
               fixed = TRUE)
  expect_no_match(e, "which leaves the estimates as they are,",
                  fixed = TRUE)
  expect_error(sem_fit_measures(a), "the refit keeps them", fixed = TRUE)
  # The refusal is classed, carries the exogenous variables, and says how
  # to keep the fit's estimates: the refit with fixed.x = FALSE on the rows
  # the fit used, those on which the exogenous variable is observed, has
  # its cases and its estimates.
  caught <- tryCatch(.ts_residual_test(a), dmar_refusal = function(e) e)
  expect_s3_class(caught, c("dmar_two_stage_fixed_x", "dmar_refusal", "error"))
  expect_identical(caught$code, "two_stage_fixed_x")
  expect_identical(caught$variables, "age_years")
  expect_match(conditionMessage(caught), paste0(
    "To keep the fit's estimates, refit with fixed.x = FALSE on the rows the ",
    "fit used, those on which every one of those variables is observed (for ",
    "a lavaan fit, the rows lavaan::lavInspect(fit, \"case.idx\") lists, one ",
    "set per group)."), fixed = TRUE)
  expect_error(sem_fit_measures(a), class = "dmar_two_stage_fixed_x")
  observed <- which(!is.na(v$age_years))
  expect_identical(sort(lavaan::lavInspect(a, "case.idx")), observed)
  kept <- lavaan::sem(mv, data = v[observed, ], missing = "two.stage",
                      fixed.x = FALSE)
  expect_identical(lavaan::lavInspect(kept, "nobs"), 293L)
  shared <- intersect(names(lavaan::coef(a)), names(lavaan::coef(kept)))
  expect_equal(lavaan::coef(kept)[shared], lavaan::coef(a)[shared], tolerance = 1e-4)
  # Under conditional.x = TRUE lavaan does not take fixed.x = FALSE, so
  # the message asks for conditional.x = FALSE as well; the refit on the
  # same rows has the fit's estimates.
  cx <- suppressWarnings(lavaan::sem(mv, data = v, missing = "two.stage",
                                     conditional.x = TRUE))
  caught <- tryCatch(.ts_residual_test(cx), dmar_refusal = function(e) e)
  expect_s3_class(caught, "dmar_two_stage_fixed_x")
  expect_identical(caught$variables, "age_years")
  expect_match(conditionMessage(caught), paste0(
    "Refit with conditional.x = FALSE and fixed.x = FALSE (lavaan does not ",
    "take fixed.x = FALSE beside conditional.x = TRUE)."), fixed = TRUE)
  expect_match(conditionMessage(caught), paste0(
    "lavaan has deleted the cases missing it under conditional.x = TRUE, and ",
    "the refit keeps them"), fixed = TRUE)
  expect_match(conditionMessage(caught), "To keep the fit's estimates, refit that way on",
               fixed = TRUE)
  expect_match(conditionMessage(caught), "under conditional.x = TRUE it deletes the same cases",
               fixed = TRUE)
  expect_error(sem_fit_measures(cx), class = "dmar_two_stage_fixed_x")
  expect_error(lavaan::sem(mv, data = v[observed, ], missing = "two.stage",
                           conditional.x = TRUE, fixed.x = FALSE),
               "conditional")
  joint <- lavaan::sem(mv, data = v[observed, ], missing = "two.stage",
                       conditional.x = FALSE, fixed.x = FALSE)
  shared <- intersect(names(lavaan::coef(cx)), names(lavaan::coef(joint)))
  expect_true("verbal~age_years" %in% shared)
  expect_equal(lavaan::coef(joint)[shared], lavaan::coef(cx)[shared], tolerance = 1e-4)
  expect_equal(.ts_residual_test(b)$statistic, 2.102897, tolerance = 1e-6)
  # lavaan writes a > 10 as a lower bound, on which the loading (about 3)
  # then sits.
  fi <- lavaan::cfa("f =~ a*t6_paragraph_comprehension + t7_sentence +
                          t9_word_meaning + t1_visual_perception
                     a > 10", data = d, std.lv = TRUE, missing = "two.stage")
  expect_error(.ts_residual_test(fi), "on a bound")
  fe <- lavaan::cfa("f =~ a*t6_paragraph_comprehension + b*t7_sentence +
                          t9_word_meaning + t1_visual_perception
                     a > b", data = d, std.lv = TRUE, missing = "two.stage")
  expect_error(.ts_residual_test(fe), "inequality constraints")
})
