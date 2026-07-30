# Weighted Kappa: Cohen's 1968 Illustration, Worked in Full

``` r
library(DMAR)
```

Cohen’s (1968) weighted kappa generalizes his (1960) kappa to situations
in which some disagreements are graver than others: each cell of the k x
k agreement matrix receives a weight, and the coefficient is the
chance-corrected proportion of weighted agreement. This vignette works
through Cohen’s own illustration from start to finish using
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)
and the `diagnosis_agreement` data set, which reproduces his Table 1
cell for cell. Along the way it covers the two weight scalings, the
standard error question, linear and quadratic weighting for ordinal
categories, Cohen’s connection between weighted kappa and the
product-moment correlation, and, at the end, an orientation question
about the paper’s validity example that the reader can adjudicate with
all computations in view.

## The Data: Cohen’s Table 1

Two judges independently assign *N* = 200 cases to three diagnostic
categories. The `diagnosis_agreement` data set stores one row per cell
of the 3 x 3 matrix, in the paper’s own layout (Judge B indexes the
rows, Judge A the columns), with the three quantities Cohen prints in
each cell: the ratio-scaled disagreement weight, the chance-expected
proportion (his parenthetical values), and the observed proportion.

``` r
data(diagnosis_agreement)
diagnosis_agreement
#>                judge_b              judge_a frequency disagreement_weight
#> 1 Personality disorder Personality disorder        88                   0
#> 2 Personality disorder             Neurosis        10                   1
#> 3 Personality disorder            Psychosis         2                   3
#> 4             Neurosis Personality disorder        14                   1
#> 5             Neurosis             Neurosis        40                   0
#> 6             Neurosis            Psychosis         6                   6
#> 7            Psychosis Personality disorder        18                   3
#> 8            Psychosis             Neurosis        10                   6
#> 9            Psychosis            Psychosis        12                   0
#>   observed_proportion expected_proportion
#> 1                0.44                0.30
#> 2                0.05                0.15
#> 3                0.01                0.05
#> 4                0.07                0.18
#> 5                0.20                0.09
#> 6                0.03                0.03
#> 7                0.09                0.12
#> 8                0.05                0.06
#> 9                0.06                0.02

tab <- xtabs(frequency ~ judge_b + judge_a, data = diagnosis_agreement)
addmargins(tab)
#>                       judge_a
#> judge_b                Personality disorder Neurosis Psychosis Sum
#>   Personality disorder                   88       10         2 100
#>   Neurosis                               14       40         6  60
#>   Psychosis                              18       10        12  40
#>   Sum                                   120       60        20 200
```

The margins are .50/.30/.20 for Judge B (rows) and .60/.30/.10 for Judge
A (columns), and each cell’s expected proportion is the product of its
row and column marginal proportions:

``` r
with(diagnosis_agreement,
     all.equal(expected_proportion,
               as.numeric(tapply(observed_proportion, judge_b, sum)[judge_b] *
                          tapply(observed_proportion, judge_a, sum)[judge_a])))
#> [1] TRUE
```

## Unweighted Kappa

With every disagreement treated alike, kappa is the proportion of
agreement corrected for chance (Cohen, 1960):

``` r
cohen_kappa(table = tab)
```

| weights | kappa | se | lower_limit | upper_limit | z_value | p_value | n | n_categories |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| unweighted | 0.492 | 0.051 | 0.392 | 0.591 | 9.64 | \< 0.0001 | 200 | 3 |

Confidence level: 95%

The estimate is .492: after chance agreement is removed from both the
numerator and the base, about half the judgments agree. The confidence
interval uses the Fleiss, Cohen, and Everitt (1969) asymptotic standard
error.

## Disagreement Scaling and Weighted Kappa

Cohen scales *disagreement*: the diagonal gets 0, a personality
disorder-neurosis confusion gets 1, personality disorder-psychosis gets
3, and neurosis-psychosis, the confusion the illustration treats as
gravest, gets 6. The weights are part of the definition of agreement, so
they must be fixed before the data are collected, and weighted kappa is
invariant to multiplying them by any positive constant.

``` r
v <- unclass(xtabs(disagreement_weight ~ judge_b + judge_a,
                   data = diagnosis_agreement))
v
#>                       judge_a
#> judge_b                Personality disorder Neurosis Psychosis
#>   Personality disorder                    0        1         3
#>   Neurosis                                1        0         6
#>   Psychosis                               3        6         0
#> attr(,"call")
#> xtabs(formula = disagreement_weight ~ judge_b + judge_a, data = diagnosis_agreement)

res <- cohen_kappa(table = tab, weights = v,
                   weight_scaling = "disagreement")
res
```

| weights | kappa | se | lower_limit | upper_limit | z_value | p_value | n | n_categories |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| custom_disagreement | 0.348 | 0.0755 | 0.2 | 0.496 | 4.61 | \< 0.0001 | 200 | 3 |

Confidence level: 95%

Weighted kappa is .348, *smaller* than the unweighted .492. Cohen chose
these values to make a point that is easy to miss: because the same
weights enter both the observed and the chance-expected weighted
disagreement, weighted kappa is fully chance corrected, and “partial
credit” does not push it upward. These judges disagree far less than
chance expectation in the mildly weighted cells and at about the chance
level in the heavily weighted ones; they disagree least where it matters
least. Interchanging the 6 and 1 weights reverses the conclusion:

``` r
v_swap <- v
v_swap[v == 6] <- 1
v_swap[v == 1] <- 6
cohen_kappa(table = tab, weights = v_swap,
            weight_scaling = "disagreement")
```

| weights | kappa | se | lower_limit | upper_limit | z_value | p_value | n | n_categories |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| custom_disagreement | 0.574 | 0.0553 | 0.465 | 0.682 | 10.4 | \< 0.0001 | 200 | 3 |

Confidence level: 95%

The per-cell quantities behind any result travel with it:

``` r
attr(res, "cells")
#>                rater_1              rater_2 observed_proportion
#> 1 Personality disorder Personality disorder                0.44
#> 2 Personality disorder             Neurosis                0.05
#> 3 Personality disorder            Psychosis                0.01
#> 4             Neurosis Personality disorder                0.07
#> 5             Neurosis             Neurosis                0.20
#> 6             Neurosis            Psychosis                0.03
#> 7            Psychosis Personality disorder                0.09
#> 8            Psychosis             Neurosis                0.05
#> 9            Psychosis            Psychosis                0.06
#>   expected_proportion    weight disagreement_weight
#> 1                0.30 1.0000000                   0
#> 2                0.15 0.8333333                   1
#> 3                0.05 0.5000000                   3
#> 4                0.18 0.8333333                   1
#> 5                0.09 1.0000000                   0
#> 6                0.03 0.0000000                   6
#> 7                0.12 0.5000000                   3
#> 8                0.06 0.0000000                   6
#> 9                0.02 1.0000000                   0
```

## Cohen’s Standard Errors and the Modern Interval

Cohen’s Formulas 10 and 13 give large-sample standard errors for
weighted kappa; the arithmetic below reproduces his printed values from
the data set’s columns.

``` r
q_o <- with(diagnosis_agreement,
            sum(disagreement_weight * observed_proportion))
q_c <- with(diagnosis_agreement,
            sum(disagreement_weight * expected_proportion))
c(q_o = q_o, q_c = q_c, kappa_w = 1 - q_o / q_c)
#>       q_o       q_c   kappa_w 
#> 0.9000000 1.3800000 0.3478261

v2_o <- with(diagnosis_agreement,
             sum(disagreement_weight^2 * observed_proportion))
v2_c <- with(diagnosis_agreement,
             sum(disagreement_weight^2 * expected_proportion))
c(se_formula_10 = sqrt((v2_o - q_o^2) / (200 * q_c^2)),
  se_formula_13 = sqrt((v2_c - q_c^2) / (200 * q_c^2)))
#> se_formula_10 se_formula_13 
#>    0.09007104    0.09159718
```

Those are his .0901 and .0916, his 95% interval is .348 plus or minus
1.96 times .0901, and his significance test is z = .348/.0916 = 3.80. A
year later, Fleiss, Cohen, and Everitt (1969) derived the standard error
now in standard use, which
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)
reports; on these data it is noticeably smaller (.0755), so the modern
interval is tighter than the one printed in the 1968 paper.

## Linear and Quadratic Weights for Ordinal Categories

When categories are ordered, two standard weight patterns need no custom
matrix: `weights = "linear"` and `weights = "quadratic"`. Both are
stated on the agreement scale internally, and the quadratic pattern is
the agreement-scale counterpart of disagreement weights proportional to
the squared category distance, $`(i - j)^2`$.

``` r
df <- as.data.frame(as.table(tab))
score_b <- rep(as.integer(df$judge_b), df$Freq)
score_a <- rep(as.integer(df$judge_a), df$Freq)

v_quad <- outer(1:3, 1:3, function(i, j) (i - j)^2)
cohen_kappa(table = tab, weights = v_quad,
            weight_scaling = "disagreement")$kappa
#> [1] 0.4545455
cohen_kappa(score_b, score_a, weights = "quadratic")$kappa
#> [1] 0.4545455
```

Cohen (1968) proves that with quadratic-pattern weights and equal
marginal distributions, weighted kappa *equals* the product-moment
correlation between the category scores, which offers a useful way to
think about what r is: a chance-corrected proportion of agreement with
disagreements weighted by squared distance. With the unequal marginals
of Table 1 the two are close but not equal, with weighted kappa the
smaller, as Cohen notes:

``` r
c(kappa_w_quadratic = cohen_kappa(score_b, score_a,
                                  weights = "quadratic")$kappa,
  pearson_r = cor(score_b, score_a))
#> kappa_w_quadratic         pearson_r 
#>         0.4545455         0.4771653
```

On a table with equal marginals the identity is exact:

``` r
eq <- matrix(c(40, 10, 10,
               10, 60, 10,
               10, 10, 40), nrow = 3, byrow = TRUE)
dfe <- as.data.frame(as.table(as.table(eq)))
t1 <- rep(as.integer(dfe$Var1), dfe$Freq)
t2 <- rep(as.integer(dfe$Var2), dfe$Freq)
c(kappa_w_quadratic = cohen_kappa(table = eq, weights = v_quad,
                                  weight_scaling = "disagreement")$kappa,
  pearson_r = cor(t1, t2))
#> kappa_w_quadratic         pearson_r 
#>               0.5               0.5
```

The quadratic pattern is also the choice under which weighted kappa is
equivalent to the intraclass correlation (Fleiss & Cohen, 1973), which
is why it is the most common weighting for ordinal scales.

## The Validity Example and an Orientation Question

Cohen closes the paper by reinterpreting Table 1 as a validity problem:
Judge A becomes the consensus diagnosis of a panel (the criterion) and
Judge B a computer diagnosis (the predictor). Because a predictor’s two
directions of confusion can carry different costs, the weights may be
asymmetric. His prose says a computer diagnosis of neurosis when the
panel consensus is psychosis (weight 6) is graver than a computer
psychosis when the panel says neurosis (weight 2), and the weight
display printed in the paper, which has Panel labeling its columns and
Computer its rows, matches that prose:

``` r
v_printed <- matrix(c(0, 1, 4,
                      1, 0, 6,
                      2, 2, 0), nrow = 3, byrow = TRUE,
                    dimnames = list(computer = c("D", "N", "P"),
                                    panel    = c("D", "N", "P")))
v_printed
#>         panel
#> computer D N P
#>        D 0 1 4
#>        N 1 0 6
#>        P 2 2 0
```

Here is the question. Table 1 has Judge B (the computer, in this
reinterpretation) on its rows and Judge A (the panel) on its columns,
the same orientation as the printed weight display, so applying the
printed weights cell by cell gives one answer; applying their transpose
gives another. Cohen reports the weighted disagreement sums as .86
(observed) and 1.33 (chance), yielding weighted kappa .353, with Formula
10 and 13 standard errors .0887 and .0915. Both computations follow,
from the same data:

``` r
# As printed, read with its own row and column labels:
res_printed <- cohen_kappa(table = tab, weights = v_printed,
                           weight_scaling = "disagreement")
cells <- attr(res_printed, "cells")
c(sum_v_p_obs = sum(cells$disagreement_weight * cells$observed_proportion),
  sum_v_p_exp = sum(cells$disagreement_weight * cells$expected_proportion),
  kappa_w     = res_printed$kappa)
#> sum_v_p_obs sum_v_p_exp     kappa_w 
#>   0.6200000   1.0700000   0.4205607

# The transpose of the printed display:
res_t <- cohen_kappa(table = tab, weights = t(v_printed),
                     weight_scaling = "disagreement")
cells_t <- attr(res_t, "cells")
c(sum_v_p_obs = sum(cells_t$disagreement_weight * cells_t$observed_proportion),
  sum_v_p_exp = sum(cells_t$disagreement_weight * cells_t$expected_proportion),
  kappa_w     = res_t$kappa)
#> sum_v_p_obs sum_v_p_exp     kappa_w 
#>   0.8600000   1.3300000   0.3533835
```

The transpose reproduces Cohen’s published quantities (.86, 1.33, .353,
and his Formula 10 and 13 values .0887 and .0915 follow from the same
sums); the printed orientation gives .62, 1.07, and .421 instead. So the
published numerical results correspond to the weights applied in the
transpose of the orientation the printed display’s labels state. What to
make of that is a judgment the numbers alone do not settle: the
possibilities include a transposed weight display in typesetting, a
transposed data table relative to the notation used in computation, or a
slip in the computation itself, and nothing in the paper distinguishes
among them. Cohen’s internal arithmetic is consistent throughout, so the
substance of the example, that asymmetric weights make weighted kappa a
validity measure, is untouched either way. The
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)
help page follows the published values (weighted kappa .353), which is
the analysis the paper reports.

## References

Cohen, J. (1960). A coefficient of agreement for nominal scales.
*Educational and Psychological Measurement, 20*(1), 37–46.

Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision for
scaled disagreement or partial credit. *Psychological Bulletin, 70*(4),
213–220.

Fleiss, J. L., & Cohen, J. (1973). The equivalence of weighted kappa and
the intraclass correlation coefficient as measures of reliability.
*Educational and Psychological Measurement, 33*(3), 613–619.

Fleiss, J. L., Cohen, J., & Everitt, B. S. (1969). Large sample standard
errors of kappa and weighted kappa. *Psychological Bulletin, 72*(5),
323–327.
