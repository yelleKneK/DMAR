# Cohen's Kappa Coefficient of Inter-Rater Agreement

Computes Cohen's (1960) kappa coefficient of agreement between two
raters on a categorical variable, with optional weighting (linear,
quadratic, or custom) for ordinal categories, following Cohen's (1968)
weighted kappa. Input can be the two raters' classification vectors or a
published \\k \times k\\ frequency table. A confidence interval and a
Wald test of \\H_0\\: \kappa = 0\\ are returned, using the asymptotic
standard error of Fleiss, Cohen, and Everitt (1969).

## Usage

``` r
cohen_kappa(
  rater_1 = NULL,
  rater_2 = NULL,
  table = NULL,
  weights = "unweighted",
  weight_scaling = c("agreement", "disagreement"),
  categories = NULL,
  conf_level = 0.95,
  ci_method = c("wald", "percentile", "bca"),
  B = 10000L,
  seed = NULL
)
```

## Arguments

- rater_1:

  First rater's classification vector (categorical, factor or coercible
  to factor).

- rater_2:

  Second rater's classification vector (same length as `rater_1`;
  categorical, factor or coercible to factor). Together `rater_1` and
  `rater_2` give the two raters' classifications on the same set of
  subjects; pairs in which either is `NA` are dropped before
  computation.

- table:

  A \\k \times k\\ frequency table of joint assignments (rows are rater
  1's categories, columns rater 2's, in the same order), supplied
  instead of `rater_1` and `rater_2`. This is the form in which
  published agreement studies usually report their data. Row and column
  `dimnames`, when present, must match and are used as the category
  labels.

- weights:

  Either the string `"unweighted"` (the default; appropriate for nominal
  categories), `"linear"`, or `"quadratic"`, or a user-supplied square
  numeric weight matrix. A custom matrix is interpreted on the scale
  named by `weight_scaling`. Asymmetric matrices are allowed; see
  *Details*.

- weight_scaling:

  How a custom `weights` matrix is scaled: `"agreement"` (the default;
  diagonal 1, decreasing off-diagonal entries, typically in \\\[0,
  1\]\\) or `"disagreement"` (Cohen's 1968 ratio-scaled disagreement
  weights: zero diagonal, larger entries for graver disagreements). The
  two scalings give identical \\\kappa_W\\; see *Details*.

- categories:

  Optional character vector listing the full category set in display
  order (used both to set the order of the confusion matrix's
  rows/columns and to map ordinal categories to integers \\1, \ldots,
  k\\ for linear and quadratic weighting). When `NULL` (the default) and
  both raters are factors with the same level set, the factor levels are
  used in their own order; otherwise the sorted union of values observed
  in `rater_1` and `rater_2` is used. When supplied, the set must be
  unique, non-missing, and contain every observed rating; an omitted
  category raises an error, because silently dropping the corresponding
  ratings while still dividing by the original sample size would corrupt
  kappa. With `table` input, supply `categories` only when the table has
  no `dimnames`.

- conf_level:

  Confidence level for the interval (default `0.95`).

- ci_method:

  Interval method: `"wald"` (the default, the asymptotic interval from
  the Fleiss, Cohen, and Everitt standard error), `"percentile"`
  (bootstrap percentile), or `"bca"` (bootstrap bias-corrected and
  accelerated).

- B:

  Number of bootstrap replications when `ci_method` is `"percentile"` or
  `"bca"` (default `10000`; ignored for `"wald"`). The BCa adjustment
  pushes the working quantiles into the tails, so reduce `B` for
  exploration, not for a reported analysis.

- seed:

  Optional integer seed for the bootstrap. The default `NULL` uses the
  current state of the random number generator; a supplied seed is set
  internally and the prior state restored on exit.

## Value

A one-row `data.frame` (class `dmar_tbl`) with columns `weights` (the
form used), `kappa`, `se` (asymptotic standard error), `lower_limit`,
`upper_limit`, `z_value`, `p_value` (Wald test of \\H_0\\: \kappa =
0\\), `n` (number of paired ratings), and `n_categories` (\\k\\).

The per-cell detail behind the coefficient travels with the result as
the `cells` attribute, in the form of Cohen's (1968) Table 1: a
`data.frame` with one row per cell of the confusion matrix giving
`rater_1` and `rater_2` (the cell's categories), `observed_proportion`,
`expected_proportion` (the product of the marginal proportions, the
cell's chance expectation), `weight` (the agreement-scale weight used in
the computation), and, when `weight_scaling = "disagreement"`, the
supplied `disagreement_weight`. Retrieve it with
`attr(result, "cells")`.

## Details

For two raters and a confusion matrix \\P\\ of joint proportions (rater
1 \\\times\\ rater 2), the weighted kappa is \$\$\kappa_W =
\frac{p_o^{(W)} - p_e^{(W)}}{1 - p_e^{(W)}}, \quad p_o^{(W)} =
\sum\_{i,j} W\_{ij}\\P\_{ij}, \quad p_e^{(W)} = \sum\_{i,j}
W\_{ij}\\p\_{i.}\\p\_{.j},\$\$ where \\p\_{i.}\\ and \\p\_{.j}\\ are the
row and column marginals. For `weights = "unweighted"` (the diagonal of
\\W\\ is 1, off- diagonal 0) this collapses to Cohen's original
formulation.

**Disagreement scaling.** Cohen (1968) develops weighted kappa by ratio
scaling *disagreement*: each cell receives a weight \\v\_{ij} \ge 0\\,
zero on the agreement diagonal, with, for example, a weight of 6
representing twice as much disagreement as 3. The weights are part of
the definition of agreement (and of any hypothesis tested about it), so
they must be fixed before the data are collected. \\\kappa_W\\ is
invariant to multiplying the \\v\_{ij}\\ by any positive constant, and a
disagreement matrix is related to an agreement matrix by \\w\_{ij} = 1 -
v\_{ij}/v\_{\max}\\ (Cohen, 1968, Footnote 3), which is the conversion
applied internally when `weight_scaling = "disagreement"`. Either
scaling therefore yields the same \\\kappa_W\\; supply whichever is more
natural.

**Asymmetric weights and validity.** Nothing in \\\kappa_W\\ requires
\\W\_{ij} = W\_{ji}\\. Symmetric weights suit reliability, where the two
sources have equal status; asymmetric weights suit validity, where one
source is a criterion and the other a predictor and the two directions
of a confusion can carry different costs (Cohen, 1968). The examples
reproduce Cohen's computer-diagnosis illustration.

**Standard error.** The Fleiss-Cohen-Everitt (1969) asymptotic variance
for weighted kappa is used: \$\$\mathrm{Var}(\hat\kappa_W) =
\frac{1}{N(1 - p_e^{(W)})^2}\Bigl\[\sum\_{i,j} P\_{ij}\bigl(W\_{ij} -
(\bar W\_{i.} + \bar W\_{.j})(1 - \hat\kappa_W)\bigr)^2 -
\bigl(\hat\kappa_W - p_e^{(W)}(1 - \hat\kappa_W)\bigr)^2\Bigr\],\$\$
with \\\bar W\_{i.} = \sum_j W\_{ij}\\p\_{.j}\\ and \\\bar W\_{.j} =
\sum_i W\_{ij}\\p\_{i.}\\. The Wald confidence interval is \\\hat\kappa
\pm z\_{1-\alpha/2}\\\widehat{\mathrm{SE}}\\. Cohen's (1968) own
Formulas 10 and 13 for the standard error of \\\kappa_W\\ preceded this
result and were superseded by it; the examples reproduce his Table 1
arithmetic for the historical record while the function reports the
Fleiss-Cohen-Everitt interval.

**Choice of weights.** Use `"unweighted"` for nominal categories. For
ordinal categories, `"quadratic"` is the most common choice (and
mathematically equivalent to the intraclass correlation under certain
conditions; Fleiss & Cohen, 1973); `"linear"` is also defensible. Cohen
(1968) further shows that with equal marginals and quadratic-pattern
disagreement weights, \\\kappa_W\\ equals the product-moment correlation
between the category scores.

**Small samples and the bootstrap.** The Wald interval can have poor
coverage for small \\N\\ or extreme values of \\\hat\kappa\\; a
bootstrap interval is more dependable in those regimes (Blackman &
Koval, 2000). With `ci_method = "percentile"` or `"bca"` the subjects
(the rated pairs) are resampled with replacement `B` times, kappa is
recomputed on each resample with the same categories and weights, and
the interval is read off the bootstrap distribution: the percentile
interval takes the empirical quantiles, and the BCa interval adjusts the
quantile positions for median bias (estimated from the bootstrap
distribution) and for acceleration (estimated from the jackknife),
making it second-order accurate where the percentile interval is
first-order accurate (Efron & Tibshirani, 1993). `table` input is
expanded to the equivalent paired ratings and resampled the same way. A
resample on which kappa is undefined (chance agreement 1) is dropped,
and the interval is computed from the replications that return a value;
a single warning reports how many were dropped. The `se`, `z_value`, and
`p_value` columns keep their asymptotic definitions under every
`ci_method`; only the interval changes. Bootstrap results vary from run
to run; supply `seed` for reproducibility (the RNG state is set locally
and the caller's state restored on exit).

## References

Blackman, N. J.-M., & Koval, J. J. (2000). Interval estimation for
Cohen's kappa as a measure of agreement. *Statistics in Medicine,
19*(5), 723–741.

Cohen, J. (1960). A coefficient of agreement for nominal scales.
*Educational and Psychological Measurement, 20*(1), 37–46.

Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision for
scaled disagreement or partial credit. *Psychological Bulletin, 70*(4),
213–220.

Fleiss, J. L., Cohen, J., & Everitt, B. S. (1969). Large sample standard
errors of kappa and weighted kappa. *Psychological Bulletin, 72*(5),
323–327.

Fleiss, J. L., & Cohen, J. (1973). The equivalence of weighted kappa and
the intraclass correlation coefficient as measures of reliability.
*Educational and Psychological Measurement, 33*(3), 613–619.

## See also

[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md)
(Cohen's 1968 Table 1 as a data set),
[`fleiss_kappa`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md)

Other reliability:
[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# ---------------------------------------------------------------------
# Cohen (1968, Table 1): two judges assign N = 200 cases to three
# diagnostic categories. The table ships as the diagnosis_agreement
# data set in the paper's own layout, Judge B in rows and Judge A in
# columns, carrying Cohen's per-cell disagreement weights, observed
# proportions, and chance-expected proportions.
data(diagnosis_agreement)
tab <- xtabs(frequency ~ judge_b + judge_a, data = diagnosis_agreement)
v   <- unclass(xtabs(disagreement_weight ~ judge_b + judge_a,
                     data = diagnosis_agreement))
tab
#>                       judge_a
#> judge_b                Personality disorder Neurosis Psychosis
#>   Personality disorder                   88       10         2
#>   Neurosis                               14       40         6
#>   Psychosis                              18       10        12

# Unweighted kappa, all disagreements equal (Cohen's Formula 4): .492.
cohen_kappa(table = tab)
#>  weights    kappa se    lower_limit upper_limit z_value p_value  n  
#>  unweighted 0.492 0.051 0.392       0.591       9.64    < 0.0001 200
#>  n_categories
#>  3           
#> 
#> Confidence level: 95%

# A bootstrap interval for the same table, which expands it to its
# 200 paired ratings and resamples the subjects. Not run here,
# because 2000 refits of kappa is more than a help page should do;
# the call is:
# cohen_kappa(table = tab, ci_method = "percentile", B = 2000,
#             seed = 113)

# Cohen's ratio-scaled disagreement weights: a neurosis-psychosis
# confusion (weight 6) is six times as grave as a personality
# disorder-neurosis confusion (weight 1). Weighted kappa = .348,
# smaller than the unweighted .492: these judges disagree less than
# chance expectation where it matters little and at about the chance
# level where it matters most.
res <- cohen_kappa(table = tab, weights = v,
                   weight_scaling = "disagreement")
res
#>  weights             kappa se     lower_limit upper_limit z_value p_value  n  
#>  custom_disagreement 0.348 0.0755 0.2         0.496       4.61    < 0.0001 200
#>  n_categories
#>  3           
#> 
#> Confidence level: 95%

# The per-cell quantities of Cohen's Table 1 travel with the result:
# observed and chance-expected proportions and both weight scalings.
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

# Interchanging the 6 and 1 weights reverses the story: kappa_w = .574.
v_swap <- v
v_swap[v == 6] <- 1
v_swap[v == 1] <- 6
cohen_kappa(table = tab, weights = v_swap,
            weight_scaling = "disagreement")
#>  weights             kappa se     lower_limit upper_limit z_value p_value  n  
#>  custom_disagreement 0.574 0.0553 0.465       0.682       10.4    < 0.0001 200
#>  n_categories
#>  3           
#> 
#> Confidence level: 95%

# Cohen's own Table 1 arithmetic (his Formulas 8, 10, and 13),
# computed straight from the data set's per-cell columns and
# reproduced for the historical record. The function reports the
# Fleiss-Cohen-Everitt (1969) standard error, which superseded
# Formulas 10 and 13.
q_o <- with(diagnosis_agreement,
            sum(disagreement_weight * observed_proportion))   # = .90
q_c <- with(diagnosis_agreement,
            sum(disagreement_weight * expected_proportion))   # = 1.38
1 - q_o / q_c                          # kappa_w = .348   (Formula 8)
#> [1] 0.3478261
v2_o <- with(diagnosis_agreement,
             sum(disagreement_weight^2 * observed_proportion))
v2_c <- with(diagnosis_agreement,
             sum(disagreement_weight^2 * expected_proportion))
sqrt((v2_o - q_o^2) / (200 * q_c^2))   # = .0901   (Formula 10)
#> [1] 0.09007104
sqrt((v2_c - q_c^2) / (200 * q_c^2))   # = .0916   (Formula 13)
#> [1] 0.09159718
(1 - q_o / q_c) + c(-1, 1) * qnorm(0.975) * 0.0901  # 95% CI [.171, .524]
#> [1] 0.1712333 0.5244188
(1 - q_o / q_c) / 0.0916               # z = 3.80, p < .001
#> [1] 3.797228

# Cohen's validity reinterpretation: Judge A is a diagnostic panel
# (the criterion), Judge B a computer diagnosis (the predictor), and
# the weights are asymmetric because the two directions of a
# confusion carry different costs. Oriented to the table's layout
# (rows = Judge B = computer, columns = Judge A = panel), the weights
# below reproduce Cohen's published quantities: sum(v * p_o) = .86,
# sum(v * p_c) = 1.33, kappa_w = .353, with his Formulas 10 and 13
# giving .0887 and .0915. The weighted kappa vignette works this
# example in full, including the orientation of the paper's printed
# weight display relative to these values.
v_validity <- matrix(
  c(0, 1, 2,
    1, 0, 2,
    4, 6, 0), nrow = 3, byrow = TRUE)
cohen_kappa(table = tab, weights = v_validity,
            weight_scaling = "disagreement")
#>  weights             kappa se     lower_limit upper_limit z_value p_value  n  
#>  custom_disagreement 0.353 0.0627 0.231       0.476       5.64    < 0.0001 200
#>  n_categories
#>  3           
#> 
#> Confidence level: 95%

# ---------------------------------------------------------------------
# Raw rater vectors and ordinal categories (quadratic agreement
# weights, e.g., Likert-style severity).
set.seed(113)
x <- sample(1:5, 100, replace = TRUE)
y <- pmin(pmax(x + sample(-1:1, 100, replace = TRUE), 1), 5)
cohen_kappa(x, y, weights = "quadratic")
#>  weights   kappa se     lower_limit upper_limit z_value p_value  n  
#>  quadratic 0.85  0.0205 0.81        0.89        41.5    < 0.0001 100
#>  n_categories
#>  5           
#> 
#> Confidence level: 95%
```
