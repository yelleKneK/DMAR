## from tests/testthat/test-krippendorff_alpha.R
local({
  # Ratio level with an identical zero pair (HIGH-05). irr expects a
  # raters-by-units matrix (the transpose of DMAR's units-by-raters input).
  ratings <- rbind(c(0, 0), c(1, 1), c(0, 1), c(1, 0))
  oracle  <- irr::kripp.alpha(t(ratings), method = "ratio")$value
  res     <- DMAR::krippendorff_alpha(ratings, level = "ratio", boot = FALSE)
  dmar    <- res$value[res$term == "krippendorff_alpha"]
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-8)))
})

## from tests/testthat/test-audit_regressions.R
local({
  # krippendorff_alpha against irr::kripp.alpha across all four metrics.
  ratings <- matrix(c(
    1, 2, 3, 3, 2, 1, 4, 1, 2, NA, NA, NA,
    1, 2, 3, 3, 2, 2, 4, 1, 2, 5,  NA, 3,
    NA, 3, 3, 3, 2, 3, 4, 2, 2, 5,  1,  NA,
    1, 2, 3, 3, 2, 4, 4, 1, 2, 5,  1,  NA
  ), nrow = 12, ncol = 4)
  for (lv in c("nominal", "ordinal", "interval", "ratio")) {
    dmar   <- DMAR::krippendorff_alpha(ratings, level = lv,
                                       boot = FALSE)$value[1]
    oracle <- irr::kripp.alpha(t(ratings), method = lv)$value
    stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-6)))
  }
})

## from tests/testthat/test-audit_regressions.R
local({
  # gwet_ac coefficient and se against irrCAC::gwet.ac1.raw across weights.
  set.seed(113)
  r1 <- sample(1:5, 60, TRUE, prob = c(.1, .2, .4, .2, .1))
  r2 <- pmin(pmax(r1 + sample(c(-1, 0, 1), 60, TRUE, prob = c(.2, .6, .2)),
                  1), 5)
  for (w in c("unweighted", "linear", "quadratic")) {
    dm <- DMAR::gwet_ac(cbind(r1, r2), weights = w)
    ir <- irrCAC::gwet.ac1.raw(data.frame(r1, r2), weights = w)$est
    dmar   <- c(dm$value[dm$term == "gwet_ac"], dm$value[dm$term == "se"])
    oracle <- c(ir$coeff.val, ir$coeff.se)
    stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-4)))
  }
})

## from tests/testthat/test-cohen_kappa.R
local({
  # Unweighted kappa on the complete 'categories' set (HIGH-08) against
  # irr::kappa2.
  r1 <- c("a", "a", "b", "b", "c", "c", "a", "b")
  r2 <- c("a", "b", "b", "b", "c", "a", "a", "c")
  oracle <- irr::kappa2(cbind(r1, r2))$value
  dmar   <- DMAR::cohen_kappa(r1, r2)$kappa
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-7)))
  dmar_reordered <- DMAR::cohen_kappa(r1, r2,
                                      categories = c("c", "b", "a"))$kappa
  stopifnot(isTRUE(all.equal(dmar_reordered, oracle, tolerance = 1e-7)))
})

## from tests/testthat/test-cohen_kappa.R
local({
  # Unweighted and weighted kappa with the FCE se on Cohen (1968) Table 1
  # against psych::cohen.kappa; psych reports the FCE interval, so its se is
  # backed out of the CI half-width.
  cohen_1968 <- matrix(
    c(88, 10,  2,
      14, 40,  6,
      18, 10, 12),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("Personality disorder", "Neurosis", "Psychosis"),
                    c("Personality disorder", "Neurosis", "Psychosis")))
  v_1968 <- matrix(c(0, 1, 3,
                     1, 0, 6,
                     3, 6, 0), nrow = 3, byrow = TRUE)

  ck_un <- suppressWarnings(psych::cohen.kappa(cohen_1968))
  res_un <- DMAR::cohen_kappa(table = cohen_1968)
  stopifnot(isTRUE(all.equal(res_un$kappa, ck_un$kappa, tolerance = 1e-8)))
  psych_se_un <- (ck_un$confid["unweighted kappa", "upper"] -
                  ck_un$confid["unweighted kappa", "estimate"]) /
                 stats::qnorm(0.975)
  stopifnot(isTRUE(all.equal(res_un$se, psych_se_un, tolerance = 1e-6)))

  w_agree <- 1 - v_1968 / max(v_1968)
  ck_w <- suppressWarnings(psych::cohen.kappa(cohen_1968, w = w_agree))
  res_w <- DMAR::cohen_kappa(table = cohen_1968, weights = v_1968,
                             weight_scaling = "disagreement")
  stopifnot(isTRUE(all.equal(res_w$kappa, ck_w$weighted.kappa,
                             tolerance = 1e-8)))
  psych_se_w <- (ck_w$confid["weighted kappa", "upper"] -
                 ck_w$confid["weighted kappa", "estimate"]) /
                stats::qnorm(0.975)
  stopifnot(isTRUE(all.equal(res_w$se, psych_se_w, tolerance = 1e-6)))
})

## from tests/testthat/test-fleiss_kappa.R
local({
  # Fleiss (1971) Table 1 against irr::kappam.fleiss. irr wants a subjects
  # by raters matrix of categories; expand the counts (rater order within a
  # subject does not matter for Fleiss's kappa).
  fleiss_1971_table1 <- matrix(c(
    0, 0, 0, 6, 0,
    0, 3, 0, 0, 3,
    0, 1, 4, 0, 1,
    0, 0, 0, 0, 6,
    0, 3, 0, 3, 0,
    2, 0, 4, 0, 0,
    0, 0, 4, 0, 2,
    2, 0, 3, 1, 0,
    2, 0, 0, 4, 0,
    0, 0, 0, 0, 6,
    1, 0, 0, 5, 0,
    1, 1, 0, 4, 0,
    0, 3, 3, 0, 0,
    1, 0, 0, 5, 0,
    0, 2, 0, 3, 1,
    0, 0, 5, 0, 1,
    3, 0, 0, 1, 2,
    5, 1, 0, 0, 0,
    0, 2, 0, 4, 0,
    1, 0, 2, 0, 3,
    0, 0, 0, 0, 6,
    0, 1, 0, 5, 0,
    0, 2, 0, 1, 3,
    2, 0, 0, 4, 0,
    1, 0, 0, 4, 1,
    0, 5, 0, 1, 0,
    4, 0, 0, 0, 2,
    0, 2, 0, 4, 0,
    1, 0, 5, 0, 0,
    0, 0, 0, 0, 6
  ), nrow = 30, byrow = TRUE)
  result <- DMAR::fleiss_kappa(fleiss_1971_table1)
  rater_fmt <- t(apply(fleiss_1971_table1, 1, function(cnt) {
    rep(seq_along(cnt), cnt)
  }))
  by_irr <- irr::kappam.fleiss(as.data.frame(rater_fmt))
  dmar   <- c(result$kappa, result$z_value)
  oracle <- c(by_irr$value, by_irr$statistic)
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-12)))
})

## from tests/testthat/test-simple_structure.R
local({
  # Per-item Hoffman complexities against psych::fa on psych::bfi[, 1:10]
  # with a two-factor oblimin solution.
  set.seed(7)
  fa <- psych::fa(psych::bfi[, 1:10], nfactors = 2, rotate = "oblimin",
                  warnings = FALSE)
  L <- unclass(fa$loadings)
  res <- DMAR::simple_structure(L)
  dmar   <- unname(attr(res, "complexity"))
  oracle <- unname(fa$complexity)
  stopifnot(isTRUE(all.equal(dmar, oracle, tolerance = 1e-8)))
})
