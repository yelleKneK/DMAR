# Live oracle comparison against GPArotation, run at release time after DMAR
# is installed (it calls DMAR::fa_exploratory) when GPArotation is installed
# locally, as the other snippets in this folder are.
## Written against GPArotation 2026.8.2 and lavaan 0.7-3.
## from tests/testthat/test-fa_exploratory_target.R
local({
  # fa_exploratory(rotation = "target") is lavaan's rotation to a partially
  # specified target ("pst"), the least squares criterion of Browne (2001,
  # Equation 13) over the specified loadings. GPArotation's targetQ() and
  # targetT() minimize the same criterion, with NA marking an unspecified
  # entry, and its pstQ() and pstT() with the specified entries marked by a
  # weight matrix, by their own gradient projection code, so rotating
  # lavaan's unrotated standardized loadings with them reproduces DMAR's
  # rotated loadings and factor correlations, up to the signs of the
  # factors whose specified values are all zero (a column of zeros leaves
  # its sign free).
  items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
             "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning",
             "t10_addition", "t12_counting_groups_of_dots",
             "t13_straight_and_curved_capitals")
  d <- DMAR::holzinger_swineford[items]
  keys <- list(spatial = items[1:3], verbal = items[4:6], speed = items[7:9])
  B <- matrix(0, 9, 3, dimnames = list(items, names(keys)))
  for (j in seq_along(keys)) B[keys[[j]], j] <- NA
  A <- lavaan::lavInspect(lavaan::efa(data = d, nfactors = 3L, rotation = "none",
                                      missing = "ml", bounds = "none",
                                      output = "lavaan"), "std")$lambda
  align <- function(L, ref) sign(colSums(L * ref))
  compare <- function(dmar, oracle, form, free_sign) {
    L_dmar <- unclass(lavaan::lavInspect(attr(dmar, "fit"), "std")$lambda)
    L_oracle <- unclass(oracle$loadings)
    s <- align(L_oracle, L_dmar)
    # Only a column of zeros may differ in sign; a nonzero value fixes it.
    stopifnot(all(s[!free_sign] == 1))
    L_oracle <- sweep(L_oracle, 2L, s, `*`)
    stopifnot(max(abs(L_oracle - L_dmar)) < 1e-6)
    if (form == "oblique") {
      Phi_oracle <- unclass(oracle$Phi) * tcrossprod(s)
      Phi_dmar <- unclass(lavaan::lavInspect(attr(dmar, "fit"), "std")$psi)
      stopifnot(max(abs(Phi_oracle - Phi_dmar)) < 1e-6)
    }
  }
  # A target of zeros alone, reflected by lavaan.
  for (form in c("oblique", "orthogonal")) {
    dmar <- DMAR::fa_exploratory(d, rotation = if (form == "oblique") "target" else
                                   "target_orthogonal", target = keys, seed = 113)
    set.seed(113)
    oracle <- if (form == "oblique") {
      GPArotation::targetQ(A, Target = B, eps = 1e-10, maxit = 100000,
                           randomStarts = 30)
    } else {
      GPArotation::targetT(A, Target = B, eps = 1e-10, maxit = 100000,
                           randomStarts = 30)
    }
    compare(dmar, oracle, form, free_sign = rep(TRUE, 3L))
  }
  # Nonzero values, which DMAR keeps specified through an explicit mask and
  # which turn lavaan's reflection off: a keyed loading targeted at 0.6 on
  # the first factor and one at 0.8 on the second, the third factor's
  # specified values all zero. pstQ() and pstT() take the same mask as
  # their weight matrix.
  B_nonzero <- B
  B_nonzero["t1_visual_perception", "spatial"] <- 0.6
  B_nonzero["t7_sentence", "verbal"] <- 0.8
  W <- 1 * !is.na(B_nonzero)
  T0 <- B_nonzero
  T0[is.na(T0)] <- 0
  for (form in c("oblique", "orthogonal")) {
    dmar <- DMAR::fa_exploratory(d, rotation = if (form == "oblique") "target" else
                                   "target_orthogonal", target = B_nonzero, seed = 113)
    stopifnot(identical(lavaan::lavInspect(attr(dmar, "fit"),
                                           "options")$rotation.args$reflect, FALSE))
    set.seed(113)
    oracle <- if (form == "oblique") {
      GPArotation::pstQ(A, W = W, Target = T0, eps = 1e-10, maxit = 100000,
                        randomStarts = 30)
    } else {
      GPArotation::pstT(A, W = W, Target = T0, eps = 1e-10, maxit = 100000,
                        randomStarts = 30)
    }
    compare(dmar, oracle, form, free_sign = c(FALSE, FALSE, TRUE))
  }
})
