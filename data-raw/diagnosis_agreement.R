# Build the diagnosis_agreement data set: Cohen's (1968) Table 1, the
# illustrative agreement matrix for weighted kappa.
#
# Source: Cohen, J. (1968). Weighted kappa: Nominal scale agreement with
# provision for scaled disagreement or partial credit. Psychological
# Bulletin, 70(4), 213-220 (Table 1, p. 214).
#
# The table's layout is preserved: Judge B indexes the rows and Judge A
# the columns. Cohen prints, in each cell, the ratio-scaled disagreement
# weight v_ij, the chance-expected proportion (in parentheses), and the
# observed proportion; N = 200. The frequencies below are the observed
# proportions times N. Reconstruction was verified against every
# quantity computed from the table in the paper: q'_o = .90,
# q'_c = 1.38, kappa = .492, kappa_w = .348, and the Formula 10 and 13
# standard errors .0901 and .0916.

categories <- c("Personality disorder", "Neurosis", "Psychosis")

diagnosis_agreement <- data.frame(
  judge_b = factor(rep(categories, each = 3), levels = categories),
  judge_a = factor(rep(categories, times = 3), levels = categories),
  frequency = c(88, 10,  2,
                14, 40,  6,
                18, 10, 12),
  disagreement_weight = c(0, 1, 3,
                          1, 0, 6,
                          3, 6, 0)
)
diagnosis_agreement$observed_proportion <-
  diagnosis_agreement$frequency / sum(diagnosis_agreement$frequency)

# Chance-expected cell proportion: the product of the cell's row (Judge
# B) and column (Judge A) marginal proportions, the parenthetical
# values in Cohen's Table 1.
p_b <- tapply(diagnosis_agreement$observed_proportion,
              diagnosis_agreement$judge_b, sum)
p_a <- tapply(diagnosis_agreement$observed_proportion,
              diagnosis_agreement$judge_a, sum)
diagnosis_agreement$expected_proportion <-
  as.numeric(p_b[diagnosis_agreement$judge_b] *
             p_a[diagnosis_agreement$judge_a])

# Checks against the printed table before saving.
stopifnot(
  sum(diagnosis_agreement$frequency) == 200,
  all.equal(as.numeric(p_b), c(.50, .30, .20)),
  all.equal(as.numeric(p_a), c(.60, .30, .10)),
  all.equal(diagnosis_agreement$expected_proportion,
            c(.30, .15, .05, .18, .09, .03, .12, .06, .02)),
  all.equal(sum(diagnosis_agreement$disagreement_weight *
                diagnosis_agreement$observed_proportion), 0.90),
  all.equal(sum(diagnosis_agreement$disagreement_weight *
                diagnosis_agreement$expected_proportion), 1.38)
)

save(diagnosis_agreement, file = "data/diagnosis_agreement.rda",
     compress = "xz")
