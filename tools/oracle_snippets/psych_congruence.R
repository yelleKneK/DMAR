# Live oracle comparison against psych, run at release time after DMAR is
# installed (it calls DMAR::factor_congruence and DMAR::fa_exploratory)
# when psych and GPArotation are installed locally, as the other snippets
# in this folder are.
## Written against psych 2.6.5 and GPArotation (psych's rotation engine).
## from tests/testthat/test-factor_congruence.R
local({
  # psych::fa.congruence() computes Tucker's coefficient between every
  # column of one loading matrix and every column of another, rounded to
  # `digits` (2 by default, so 15 here). It is fed matrices, never lavaan
  # fits: for a fit it reads the unstandardized lambda of @Model@GLIST.
  Lx <- matrix(c(.71, .64, .58, .12, .05, -.08, .20, .11,
                 .09, .14, -.03, .68, .74, .61, .02, .19,
                 .05, -.10, .21, .07, .12, .03, .66, .59), 8, 3,
               dimnames = list(paste0("v", 1:8), c("X1", "X2", "X3")))
  Ly <- matrix(c(.10, .18, .02, .70, .69, .55, .08, .25,
                 -.62, -.70, -.49, -.15, .02, .11, -.30, -.04,
                 .12, -.05, .30, .02, .20, -.06, .58, .71), 8, 3,
               dimnames = list(paste0("v", 1:8), c("Y1", "Y2", "Y3")))
  res <- DMAR::factor_congruence(Lx, Ly)
  stopifnot(isTRUE(all.equal(attr(res, "congruence_matrix"),
                             psych::fa.congruence(Lx, Ly, digits = 15),
                             tolerance = 1e-12)))
  # Items in another order are aligned by name before the comparison.
  shuffled <- Ly[c(5, 2, 8, 1, 4, 7, 3, 6), ]
  stopifnot(isTRUE(all.equal(attr(DMAR::factor_congruence(Lx, shuffled),
                                  "congruence_matrix"),
                             psych::fa.congruence(Lx, Ly, digits = 15),
                             tolerance = 1e-12)))

  # Two solutions of the same items: psych's minimum residual oblimin
  # solution and DMAR's maximum likelihood geomin solution, read from the
  # fa_exploratory() table and from its lavaan fit alike.
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning",
             "t10_addition", "t12_counting_groups_of_dots",
             "t13_straight_and_curved_capitals")
  d <- DMAR::holzinger_swineford[items]
  f_psych <- psych::fa(stats::cor(d), nfactors = 3, fm = "minres",
                       rotate = "oblimin", n.obs = nrow(d))
  f_dmar <- DMAR::fa_exploratory(d, nfactors = 3, rotation = "geomin",
                                 seed = 113)
  L_psych <- unclass(f_psych$loadings)
  L_dmar <- unclass(lavaan::lavInspect(attr(f_dmar, "fit"), "std")$lambda)
  C <- psych::fa.congruence(L_psych, L_dmar, digits = 15)
  for (y in list(f_dmar, attr(f_dmar, "fit"))) {
    stopifnot(isTRUE(all.equal(unname(attr(DMAR::factor_congruence(L_psych, y),
                                           "congruence_matrix")),
                               unname(C), tolerance = 1e-12)))
  }

  # psych's factor.match() pairs each factor of the first solution with
  # the factor of the second that has the largest absolute coefficient
  # and reverses the second's signs to match. Where those best partners
  # are mutual, as here, that is the optimal one-to-one assignment, and
  # the aligned loadings agree; the second solution is psych's own maximum
  # likelihood geomin solution with its factors listed in another order.
  # The factor correlations are not compared: fa.organize(), which
  # factor.match() calls, reverses a factor whose loadings sum below zero
  # without reversing its correlations (psych 2.6.5).
  f_y <- psych::fa(stats::cor(d), nfactors = 3, fm = "ml", rotate = "geominQ",
                   n.obs = nrow(d))
  f_y <- psych::fa.organize(f_y, o = c(3, 1, 2))
  L_y <- unclass(f_y$loadings)
  res <- DMAR::factor_congruence(L_psych, L_y)
  C_y <- psych::fa.congruence(L_psych, L_y, digits = 15)
  stopifnot(isTRUE(all.equal(attr(res, "congruence_matrix"), C_y,
                             tolerance = 1e-12)))
  best <- apply(abs(C_y), 1L, which.max)
  stopifnot(!anyDuplicated(best), !identical(unname(best), 1:3))
  stopifnot(identical(res$factor_y, colnames(C_y)[best]))
  matched <- psych:::factor.match(f_psych, f_y)
  stopifnot(isTRUE(all.equal(unname(attr(res, "loadings_aligned")),
                             unname(unclass(matched$f$loadings)),
                             tolerance = 1e-12)))
})
