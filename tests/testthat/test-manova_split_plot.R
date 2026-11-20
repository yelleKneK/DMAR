# manova_split_plot() is a tidy wrapper around car::Anova for a split-plot
# (mixed) design. We verify that (a) the function runs and the output has
# the documented structure, (b) it reproduces the split-plot examples of
# Maxwell, Delaney, and Kelley (2027, Section 14.6), (c) its default Type
# III rows are those of car::Anova with the between-subjects factor coded by
# contr.sum, whatever options("contrasts") says, (d) Type III tests the
# unweighted and Type II the weighted grand means of the within-subjects
# contrasts, with two groups and with three, and the interaction is the
# one-way MANOVA of those contrasts, (e) the values match car's printed
# table, and (f) the result is a DMAR table whose as.data.frame() is the
# plain data frame the function returned before it took that class (with
# the Type III B rows of (c), not those of DMAR 1.0.x), prints its
# p-values to four decimals, and answers tidy() and glance(). The car
# package is in Suggests, so the tests are gated on
# skip_if_not_installed().

mm_fixture <- function() {
  set.seed(113)
  n_per <- 10
  data.frame(
    subject = factor(1:(2 * n_per)),
    group   = factor(rep(c("A", "B"), each = n_per)),
    t1      = c(rnorm(n_per, 0,   1), rnorm(n_per, 0,   1)),
    t2      = c(rnorm(n_per, 0.4, 1), rnorm(n_per, 0.8, 1)),
    t3      = c(rnorm(n_per, 0.8, 1), rnorm(n_per, 1.6, 1))
  )
}

# Reaction times (ms) of 10 young and 10 old participants at angles of 0,
# 4, and 8 degrees: Table 14.9 of Maxwell, Delaney, and Kelley (2027,
# p. 794). Their Table 14.7 (p. 786) is the 0 and 8 degree columns.
mdk_reaction_time <- function() {
  data.frame(
    age     = factor(rep(c("young", "old"), each = 10),
                     levels = c("young", "old")),
    angle_0 = c(450, 390, 570, 450, 510, 360, 510, 510, 510, 510,
                420, 600, 450, 630, 420, 600, 630, 480, 690, 510),
    angle_4 = c(510, 480, 630, 660, 660, 450, 600, 660, 660, 540,
                570, 720, 540, 660, 570, 780, 690, 570, 750, 690),
    angle_8 = c(630, 540, 660, 720, 630, 450, 720, 780, 660, 660,
                690, 810, 690, 780, 780, 870, 870, 720, 900, 810)
  )
}
angles <- c("angle_0", "angle_4", "angle_8")

# The same data with the last four young participants left out, so that the
# groups differ in size (6 and 10) and weighted and unweighted grand means
# differ.
mdk_reaction_time_unequal <- function() mdk_reaction_time()[-(7:10), ]

msp_rows  <- function(res, eff) res[res$effect == eff, ]
msp_wilks <- function(res, eff)
  res$statistic_value[res$effect == eff & res$statistic_name == "Wilks"]

test_that("manova_split_plot() returns the documented columns and effects", {
  skip_if_not_installed("car")
  d <- mm_fixture()
  res <- manova_split_plot(d, within = c("t1", "t2", "t3"), between = "group")

  expect_named(res, c("effect", "statistic_name", "statistic_value",
                      "F_approx", "df_1", "df_2", "p_value"))
  expect_setequal(res$effect,
                  c("A", "B", "interaction", "sum_of_squares_type"))

  tests <- res[res$effect != "sum_of_squares_type", ]
  expect_setequal(tests$statistic_name,
                  c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"))
  # 3 effects x 4 statistics = 12 test rows, plus 1 metadata row.
  expect_equal(nrow(tests), 12L)
  expect_equal(nrow(res), 13L)
})

test_that("manova_split_plot() reports the sum-of-squares type numerically", {
  skip_if_not_installed("car")
  d <- mm_fixture()

  # Default is Type III.
  res <- manova_split_plot(d, within = c("t1", "t2", "t3"), between = "group")
  type_row <- res[res$effect == "sum_of_squares_type", ]
  expect_equal(nrow(type_row), 1L)
  expect_equal(type_row$statistic_value, 3)
  # The value column stays numeric across the whole table.
  expect_type(res$statistic_value, "double")

  # ss_type = 2 changes the car call and is reported as 2 (smoke test
  # that the alternative type runs end to end).
  res2 <- manova_split_plot(d, within = c("t1", "t2", "t3"),
                       between = "group", ss_type = 2)
  expect_equal(
    res2$statistic_value[res2$effect == "sum_of_squares_type"], 2)
  expect_type(res2$statistic_value, "double")

  # The Roman-numeral string form is accepted and normalized.
  res3 <- manova_split_plot(d, within = c("t1", "t2", "t3"),
                       between = "group", ss_type = "III")
  expect_equal(
    res3$statistic_value[res3$effect == "sum_of_squares_type"], 3)

  # Type I is not available for this multivariate path: car computes
  # only Type II and Type III, so requesting it errors clearly.
  expect_error(
    manova_split_plot(d, within = c("t1", "t2", "t3"),
                 between = "group", ss_type = 1),
    "Type I")
  expect_error(
    manova_split_plot(d, within = c("t1", "t2", "t3"),
                 between = "group", ss_type = "I"),
    "Type I")

  # An out-of-range type is rejected before reaching car.
  expect_error(
    manova_split_plot(d, within = c("t1", "t2", "t3"),
                 between = "group", ss_type = 4),
    "ss_type")
})

test_that("manova_split_plot() values match car's printed table", {
  skip_if_not_installed("car")
  d <- mm_fixture()
  res <- manova_split_plot(d, within = c("t1", "t2", "t3"), between = "group")

  # Reference values are those of summary(car::Anova(...)) for
  # lm(cbind(t1, t2, t3) ~ group) with group coded by contr.sum, under
  # R 4.6.1 with car 3.1.5; they only depend on the seed-113 data.
  # When the four statistics agree on s = 1 (the between effect with
  # 1 df), Pillai + Wilks = 1 and Hotelling-Lawley = Roy.
  between <- res[res$effect == "A", ]
  expect_equal(between$F_approx, rep(1.680826, 4L), tolerance = 1e-5)
  expect_equal(between$df_1,     rep(1, 4L))
  expect_equal(between$df_2,     rep(18, 4L))
  expect_equal(round(between$statistic_value[1] +     # Pillai
                     between$statistic_value[2], 6),  # Wilks
               1)
  expect_equal(between$statistic_value[3],            # HL
               between$statistic_value[4],            # Roy
               tolerance = 1e-12)

  # The within-subjects main effect averages over the two groups. DMAR
  # 1.0.x fitted R's default treatment coding and reported F = 5.012657
  # here, the simple effect of time in group A.
  within_eff <- res[res$effect == "B", ]
  expect_equal(within_eff$F_approx, rep(10.46870, 4L), tolerance = 1e-5)
  expect_equal(within_eff$df_1,     rep(2, 4L))
  expect_equal(within_eff$df_2,     rep(17, 4L))
})

test_that("manova_split_plot() reproduces the 2 x 3 split-plot example of Maxwell, Delaney, and Kelley (2027)", {
  skip_if_not_installed("car")
  d <- mdk_reaction_time()
  res <- manova_split_plot(d, within = angles, between = "age")
  # DMAR 1.0.x passed on car's "HF eps > 1" warning here, which belongs to
  # the univariate tests the function does not report.
  expect_no_warning(manova_split_plot(d, within = angles, between = "age"))

  # Age main effect, the one-way ANOVA of each participant's mean
  # (Equation 14.40): F = 7.28 on 1 and 18 df, p = .0147 (p. 794).
  A <- msp_rows(res, "A")
  expect_equal(round(A$F_approx, 2), rep(7.28, 4L))
  expect_equal(A$df_1, rep(1, 4L))
  expect_equal(A$df_2, rep(18, 4L))
  expect_equal(round(A$p_value, 4), rep(0.0147, 4L))

  # Angle main effect (Equation 14.54, p. 798): |E(F)| = 8,323,430,400 and
  # |E(R)| = 116,677,584,000, so Wilks's lambda is their ratio and
  # F = 110.65 on 2 and 17 df, p < .0001. DMAR 1.0.x returned F = 36.22
  # here, the simple effect of angle for the young participants against the
  # pooled error matrix (Equation 14.66, p. 810).
  det_EF   <- 8323430400
  det_ER_B <- 116677584000
  B <- msp_rows(res, "B")
  expect_equal(msp_wilks(res, "B"), det_EF / det_ER_B, tolerance = 1e-10)
  expect_equal(B$F_approx, rep(((det_ER_B - det_EF) / 2) / (det_EF / 17), 4L),
               tolerance = 1e-10)
  expect_equal(round(B$F_approx[1], 2), 110.65)
  expect_equal(B$df_1, rep(2, 4L))
  expect_equal(B$df_2, rep(17, 4L))
  expect_true(all(B$p_value < 1e-4))

  # Age by angle interaction (Equation 14.58, p. 800): |E(R)| =
  # 14,874,192,000 and F = 6.69 on 2 and 17 df, p = .0072 (p. 801), and
  # Pillai's trace V = .4404 (Equation 14.63, p. 804).
  det_ER_AB <- 14874192000
  AB <- msp_rows(res, "interaction")
  expect_equal(msp_wilks(res, "interaction"), det_EF / det_ER_AB,
               tolerance = 1e-10)
  expect_equal(AB$F_approx, rep(((det_ER_AB - det_EF) / 2) / (det_EF / 17), 4L),
               tolerance = 1e-10)
  expect_equal(round(AB$F_approx[1], 2), 6.69)
  expect_equal(round(AB$p_value[1], 4), 0.0072)
  expect_equal(round(AB$statistic_value[AB$statistic_name == "Pillai"], 4),
               0.4404)

  # With equal group sizes the weighted and unweighted grand means agree, so
  # Type II gives the same table.
  res2 <- manova_split_plot(d, within = angles, between = "age", ss_type = 2)
  tests <- res$effect != "sum_of_squares_type"
  expect_equal(res2[tests, ], res[tests, ], tolerance = 1e-10)
})

test_that("manova_split_plot() reproduces the two-angle example of Maxwell, Delaney, and Kelley (2027)", {
  skip_if_not_installed("car")
  d <- mdk_reaction_time()
  res <- manova_split_plot(d, within = c("angle_0", "angle_8"),
                           between = "age")

  # Each effect has one transformed variable, so each F is a univariate F on
  # 1 and 18 df, formed here from the book's printed means and sums of
  # squares (Table 14.8, p. 788). Age (Equation 14.40, p. 787): means of M
  # of 561 and 667.5, grand mean 614.25, within-group sum of squares
  # 110,452.5; F = 9.24, p = .007 (p. 788).
  F_A <- (10 * (561 - 614.25)^2 + 10 * (667.5 - 614.25)^2) / (110452.5 / 18)
  # Interaction (Equation 14.44, p. 789): means of D of 168 and 249, grand
  # mean 208.5, within-group sum of squares 67,050; F = 8.81, p = .008
  # (p. 790).
  F_AB <- (10 * (168 - 208.5)^2 + 10 * (249 - 208.5)^2) / (67050 / 18)
  # Angle (Equation 14.48, p. 792): N times the squared grand mean of D over
  # the same error; F = 233.41. DMAR 1.0.x returned 75.77, the angle effect
  # of the young participants.
  F_B <- 20 * 208.5^2 / (67050 / 18)

  expect_equal(round(c(F_A, F_AB, F_B), 2), c(9.24, 8.81, 233.41))
  expect_equal(msp_rows(res, "A")$F_approx, rep(F_A, 4L), tolerance = 1e-10)
  expect_equal(msp_rows(res, "interaction")$F_approx, rep(F_AB, 4L),
               tolerance = 1e-10)
  expect_equal(msp_rows(res, "B")$F_approx, rep(F_B, 4L), tolerance = 1e-10)
  expect_equal(round(msp_rows(res, "A")$p_value[1], 3), 0.007)
  expect_equal(round(msp_rows(res, "interaction")$p_value[1], 3), 0.008)
})

test_that("manova_split_plot() Type III rows are those of a contr.sum fit, whatever options(\"contrasts\") says", {
  skip_if_not_installed("car")
  d <- mdk_reaction_time_unequal()
  res <- manova_split_plot(d, within = angles, between = "age")

  # The explicit fit: the between-subjects factor coded by contr.sum, Type
  # III tests through car. Wilks's lambda and Pillai's trace are formed from
  # the hypothesis and error matrices car returns for each term.
  Y <- as.matrix(d[, angles])
  fit <- lm(Y ~ age, data = d, contrasts = list(age = contr.sum))
  mv <- summary(car::Anova(fit,
                           idata = data.frame(Time = factor(angles, levels = angles)),
                           idesign = ~ Time, type = 3),
                multivariate = TRUE, univariate = FALSE)$multivariate.tests
  car_terms <- c(A = "age", B = "Time", interaction = "age:Time")
  for (eff in names(car_terms)) {
    H <- mv[[car_terms[[eff]]]]$SSPH
    E <- mv[[car_terms[[eff]]]]$SSPE
    rows <- msp_rows(res, eff)
    expect_equal(rows$statistic_value[rows$statistic_name == "Wilks"],
                 det(E) / det(E + H), tolerance = 1e-10)
    expect_equal(rows$statistic_value[rows$statistic_name == "Pillai"],
                 sum(diag(H %*% solve(E + H))), tolerance = 1e-10)
  }

  # The global contrasts option, and an ordered between-subjects factor,
  # leave the table as it is.
  for (oc in list(c("contr.treatment", "contr.poly"),
                  c("contr.SAS", "contr.poly"),
                  c("contr.helmert", "contr.poly"))) {
    withr::local_options(contrasts = oc)
    expect_equal(manova_split_plot(d, within = angles, between = "age"), res,
                 tolerance = 1e-12)
  }
  d_ordered <- d
  d_ordered$age <- factor(d_ordered$age, levels = c("young", "old"),
                          ordered = TRUE)
  expect_equal(manova_split_plot(d_ordered, within = angles, between = "age"),
               res, tolerance = 1e-12)
})

test_that("manova_split_plot() Type III tests the unweighted and Type II the weighted grand means", {
  skip_if_not_installed("car")
  d <- mdk_reaction_time_unequal()
  res3 <- manova_split_plot(d, within = angles, between = "age")
  res2 <- manova_split_plot(d, within = angles, between = "age", ss_type = 2)

  # By hand, as Maxwell, Delaney, and Kelley (2027, Section 14.6) form the
  # tests: D variables for the linear and quadratic trends of angle, the
  # full model's error matrix E(F) pooled over the groups, and the angle
  # main effect's restricted matrix E(F) + a n~ Dbar_U Dbar_U' for the
  # unweighted grand means Dbar_U, with n~ the harmonic mean of the group
  # sizes (Chapter 14, notes 10 and 14), or E(F) + N Dbar Dbar' for the
  # weighted grand means Dbar (note 11).
  D <- as.matrix(d[, angles]) %*% cbind(c(-1, 0, 1), c(1, -2, 1))
  by_group <- split(as.data.frame(D), d$age)
  n_j <- vapply(by_group, nrow, integer(1L))
  E_F <- Reduce(`+`, lapply(by_group, function(z)
    crossprod(scale(as.matrix(z), scale = FALSE))))
  D_U <- colMeans(t(vapply(by_group, colMeans, numeric(2L))))
  D_W <- colMeans(D)
  a <- length(n_j)
  n_tilde <- a / sum(1 / n_j)
  expect_false(isTRUE(all.equal(D_U, D_W)))

  expect_equal(msp_wilks(res3, "B"),
               det(E_F) / det(E_F + a * n_tilde * tcrossprod(D_U)),
               tolerance = 1e-10)
  expect_equal(msp_wilks(res2, "B"),
               det(E_F) / det(E_F + sum(n_j) * tcrossprod(D_W)),
               tolerance = 1e-10)

  # The between-subjects effect and the interaction do not depend on the
  # type; the within-subjects main effect does.
  for (eff in c("A", "interaction"))
    expect_equal(msp_rows(res2, eff), msp_rows(res3, eff), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(msp_wilks(res2, "B"), msp_wilks(res3, "B"))))
})

test_that("manova_split_plot() forms each test from the D variables with three groups of unequal size", {
  skip_if_not_installed("car")
  # Three groups of 5, 8, and 12 measured on four occasions. With three
  # groups, a times the harmonic mean of the group sizes (22.04) falls
  # further below N = 25 than with two, and the interaction, with a - 1 = 2
  # and b - 1 = 3, has s = 2, where the four statistics give different F
  # values.
  set.seed(113)
  n_j <- c(5L, 8L, 12L)
  g <- factor(rep(c("g1", "g2", "g3"), times = n_j))
  profile <- rbind(c(0,   0.4, 0.8, 1.2),
                   c(0,   0.8, 1.2, 1.4),
                   c(0.5, 0.5, 0.6, 1.5))
  Y <- profile[as.integer(g), ] + matrix(rnorm(sum(n_j) * 4L), ncol = 4L)
  occasions <- paste0("t", 1:4)
  colnames(Y) <- occasions
  d <- data.frame(group = g, Y)
  res3 <- manova_split_plot(d, within = occasions, between = "group")
  res2 <- manova_split_plot(d, within = occasions, between = "group",
                            ss_type = 2)

  # Any b - 1 linearly independent contrasts serve as the D variables. The
  # within-subjects main effect: E(F) + a n~ Dbar_U Dbar_U' for Type III
  # and E(F) + N Dbar Dbar' for Type II, as in the two-group test above.
  D <- Y %*% contr.poly(4L)
  E_F <- crossprod(residuals(lm(D ~ g)))
  D_U <- colMeans(apply(D, 2L, function(z) tapply(z, g, mean)))
  D_W <- colMeans(D)
  a <- length(n_j)
  n_tilde <- a / sum(1 / n_j)
  expect_equal(msp_wilks(res3, "B"),
               det(E_F) / det(E_F + a * n_tilde * tcrossprod(D_U)),
               tolerance = 1e-10)
  expect_equal(msp_wilks(res2, "B"),
               det(E_F) / det(E_F + sum(n_j) * tcrossprod(D_W)),
               tolerance = 1e-10)

  # The interaction is the one-way MANOVA of the D variables over the
  # groups, under either type: each statistic, its F, and its degrees of
  # freedom are those of summary.manova().
  mv <- manova(D ~ g)
  for (st in c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")) {
    ref <- unname(summary(mv, test = st)$stats["g", 2:5])
    for (res in list(res3, res2)) {
      row <- msp_rows(res, "interaction")
      row <- row[row$statistic_name == st, ]
      expect_equal(c(row$statistic_value, row$F_approx, row$df_1, row$df_2),
                   ref, tolerance = 1e-10)
    }
  }

  # The between-subjects effect is the one-way ANOVA of each subject's mean.
  F_M <- anova(lm(rowMeans(Y) ~ g))[["F value"]][1]
  expect_equal(msp_rows(res3, "A")$F_approx, rep(F_M, 4L), tolerance = 1e-10)
  expect_equal(msp_rows(res2, "A"), msp_rows(res3, "A"), tolerance = 1e-12)
})

test_that("manova_split_plot() column names do not reach the model it fits", {
  skip_if_not_installed("car")
  d <- mdk_reaction_time_unequal()
  ref <- manova_split_plot(d, within = angles, between = "age")

  # A between-subjects column named as the model's internal objects, or
  # named in a way that is not syntactic, and within-subjects columns that
  # include a Y. DMAR 1.0.x reported the between-subjects test in the B rows
  # for a column named "Time" and stopped with an error for the others but
  # "group".
  for (nm in c("Time", "Y", "group", "age group")) {
    d_nm <- d
    names(d_nm)[names(d_nm) == "age"] <- nm
    expect_equal(manova_split_plot(d_nm, within = angles, between = nm), ref,
                 tolerance = 1e-12)
  }
  d_xyz <- d
  names(d_xyz)[match(angles, names(d_xyz))] <- c("X", "Y", "Z")
  expect_equal(manova_split_plot(d_xyz, within = c("X", "Y", "Z"),
                                 between = "age"),
               ref, tolerance = 1e-12)
})

test_that("manova_split_plot() returns a DMAR table that holds its plain data frame", {
  skip_if_not_installed("car")
  res <- manova_split_plot(mdk_reaction_time(), within = angles,
                           between = "age")

  # DMAR 1.0.x returned a plain data.frame, whose B rows here held the
  # simple effect of angle among the young participants (F = 36.22; see
  # the test of the 2 x 3 example above). The table is now a dmar_tbl, and
  # as.data.frame() gives the plain data frame back: the same columns,
  # values, and row names, with no attribute beside them.
  expect_identical(class(res), c("dmar_tbl", "data.frame"))
  plain <- as.data.frame(res)
  expect_identical(class(plain), "data.frame")
  expect_identical(unclass(res), unclass(plain))
  expect_setequal(names(attributes(plain)), c("names", "row.names", "class"))
  expect_true(is.data.frame(res))
})

test_that("manova_split_plot() prints its p-values to four decimals", {
  skip_if_not_installed("car")
  res <- manova_split_plot(mdk_reaction_time(), within = angles,
                           between = "age")

  # A width that keeps each row on one line.
  withr::local_options(width = 200)

  # As a plain data frame, the form of DMAR 1.0.x, the table prints the
  # p_value column in scientific notation (1.473401e-02 for age), since the
  # angle effect's p-value is below 1e-9. DMAR 1.0.x printed this table in
  # scientific notation too, with 7.433013e-07 in the angle rows, the
  # p-value of its simple effect among the young participants.
  expect_true(any(grepl("1.473401e-02",
                        capture.output(print(as.data.frame(res))),
                        fixed = TRUE)))

  out <- capture.output(print(res))
  expect_false(any(grepl("[0-9]e[-+][0-9]", out)))
  # Age, F(1, 18) = 7.28, p = .0147 (p. 794); angle, F(2, 17) = 110.65,
  # p < .0001 (p. 798); and their interaction, F(2, 17) = 6.69, p = .0072
  # (p. 801). The statistics show 3 decimal places (7.276, 110.653, and
  # 6.690, which round to the book's values) and the degrees of freedom no
  # decimal part.
  expect_true(any(grepl("^ *A +Pillai +0\\.288 +7\\.276 +1 +18 +0\\.0147 *$",
                        out)))
  expect_true(any(grepl("^ *B +Wilks +0\\.071 +110\\.653 +2 +17 +< 0\\.0001 *$",
                        out)))
  expect_true(any(grepl(
    "^ *interaction +Wilks +0\\.560 +6\\.690 +2 +17 +0\\.0072 *$", out)))
  # The sum-of-squares type row shows the type and nothing else.
  expect_true(any(grepl("^ *sum_of_squares_type +<NA> +3( +<NA>){4} *$",
                        out)))
})

test_that("manova_split_plot() answers tidy() and glance() through the default methods", {
  skip_if_not_installed("car")
  res <- manova_split_plot(mdk_reaction_time(), within = angles,
                           between = "age")

  # One row per row of the table: the term joins the effect and the
  # statistic, the estimate is the statistic, and the other columns follow
  # at full precision.
  td <- generics::tidy(res)
  expect_identical(names(td), c("term", "estimate", "F_approx", "df_1",
                                "df_2", "p_value"))
  expect_identical(td$term, paste(res$effect, res$statistic_name, sep = ":"))
  expect_identical(td$term[1:4], c("A:Pillai", "A:Wilks",
                                   "A:Hotelling-Lawley", "A:Roy"))
  expect_identical(td$estimate, res$statistic_value)
  for (nm in c("F_approx", "df_1", "df_2", "p_value"))
    expect_identical(td[[nm]], res[[nm]])
  # The row that records the sum-of-squares type keeps its place. Its
  # statistic_name is NA, so its term ends in "NA"; its estimate is the
  # type, and its other columns are NA.
  last <- nrow(td)
  expect_identical(td$term[last], "sum_of_squares_type:NA")
  expect_identical(td$estimate[last], 3)
  expect_true(all(is.na(unlist(td[last, c("F_approx", "df_1", "df_2",
                                          "p_value")]))))

  # glance() is one row: the number of rows, with no level and no
  # bootstrap count, since no interval or bootstrap is formed; and
  # results_sentence() has no sentence to write.
  gl <- generics::glance(res)
  expect_identical(names(gl), c("n_terms", "conf_level", "B_used"))
  expect_identical(nrow(gl), 1L)
  expect_identical(gl$n_terms, nrow(res))
  expect_true(is.na(gl$conf_level))
  expect_true(is.na(gl$B_used))
  expect_error(results_sentence(res), "confidence interval")
})

test_that("manova_split_plot() rejects bad input", {
  skip_if_not_installed("car")
  d <- mm_fixture()
  expect_error(manova_split_plot(d, within = "t1", between = "group"),
               "2\\+ column names")
  expect_error(manova_split_plot(d, within = c("t1", "t2"), between = "not_there"),
               "single column name")
  # Between factor with 1 level.
  d_bad <- d
  d_bad$group <- factor("A", levels = "A")
  expect_error(manova_split_plot(d_bad, within = c("t1", "t2"), between = "group"),
               "at least 2 levels")
})
