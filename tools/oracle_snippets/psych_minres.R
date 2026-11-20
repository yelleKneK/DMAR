# Live oracle comparison against psych, run at release time after DMAR is
# installed (it calls DMAR::fa_exploratory) when psych and GPArotation are installed
# locally, as the other snippets in this folder are.
## Written against psych 2.6.5 and GPArotation (psych's rotation engine).
## from tests/testthat/test-factor_analysis.R
local({
  # fa_exploratory(estimator = "ULS") is the minimum residual solution of Harman and
  # Jones (1966): unweighted least squares on the correlation matrix. psych's
  # fa(fm = "minres") minimizes the same off-diagonal residuals and rotates
  # with GPArotation's oblimin without row weights, as lavaan's efa() does,
  # so the rotated loadings, factor correlations, and communalities agree up
  # to the order and signs of the factors.
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning",
             "t10_addition", "t12_counting_groups_of_dots",
             "t13_straight_and_curved_capitals")
  d <- DMAR::holzinger_swineford[items]
  R <- stats::cor(d)
  oracle <- psych::fa(R, nfactors = 3, fm = "minres", rotate = "oblimin",
                      n.obs = nrow(d))
  dmar <- DMAR::fa_exploratory(d, nfactors = 3, estimator = "ULS", seed = 113)
  L_dmar <- lavaan::lavInspect(attr(dmar, "fit"), "std")$lambda
  L_psych <- unclass(oracle$loadings)
  # Match the factors by the largest absolute column correlation, then align
  # the signs.
  m <- abs(stats::cor(L_psych, L_dmar))
  idx <- apply(m, 1L, which.max)
  stopifnot(!anyDuplicated(idx))
  L_matched <- L_dmar[, idx, drop = FALSE]
  sgn <- sign(colSums(L_psych * L_matched))
  L_matched <- sweep(L_matched, 2L, sgn, `*`)
  stopifnot(isTRUE(all.equal(unname(L_psych), unname(L_matched), tolerance = 1e-4)))
  stopifnot(isTRUE(all.equal(unname(oracle$communality),
                             unname(dmar$estimate[dmar$op == "communality"]),
                             tolerance = 1e-5)))
  Phi_dmar <- lavaan::lavInspect(attr(dmar, "fit"), "std")$psi[idx, idx]
  Phi_dmar <- sweep(sweep(Phi_dmar, 1L, sgn, `*`), 2L, sgn, `*`)
  stopifnot(isTRUE(all.equal(unname(oracle$Phi), unname(Phi_dmar), tolerance = 1e-4)))
  # From a supplied covariance matrix the loadings are the same (the
  # correlation matrix is what ULS fits) and the standard errors are NA.
  from_S <- DMAR::fa_exploratory(S = stats::cov(d), N = nrow(d), nfactors = 3,
                      estimator = "ULS", seed = 113)
  stopifnot(isTRUE(all.equal(dmar$estimate[dmar$op == "communality"],
                             from_S$estimate[from_S$op == "communality"],
                             tolerance = 1e-5)))
  stopifnot(all(is.na(from_S$se)))
})
