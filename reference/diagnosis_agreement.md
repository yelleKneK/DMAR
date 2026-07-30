# Cohen's (1968) Psychiatric Diagnosis Agreement Table

The illustrative agreement matrix from Cohen's (1968) weighted kappa
paper, Table 1: two judges independently assign *N* = 200 cases to three
diagnostic categories (personality disorder, neurosis, psychosis). The
data set reproduces the printed table cell for cell and in its original
layout, Judge B indexing the rows and Judge A the columns, one row per
cell of the 3 x 3 matrix. Each cell carries the three quantities Cohen
prints: the ratio-scaled disagreement weight, the chance-expected
proportion (his parenthetical values), and the observed proportion; the
raw frequency is the observed proportion times *N*.

## Usage

``` r
diagnosis_agreement
```

## Format

A data frame with 9 observations (one per cell of the 3 x 3 agreement
matrix) on 6 variables.

- `judge_b`:

  Factor: Judge B's diagnostic category (the table's rows), with levels
  `Personality disorder`, `Neurosis`, `Psychosis`.

- `judge_a`:

  Factor: Judge A's diagnostic category (the table's columns), same
  levels.

- `frequency`:

  Number of the 200 cases jointly assigned to the cell.

- `disagreement_weight`:

  Cohen's ratio-scaled disagreement weight \\v\_{ij}\\ for the cell: 0
  on the agreement diagonal, 1 for a personality disorder-neurosis
  confusion, 3 for personality disorder-psychosis, and 6 for
  neurosis-psychosis, the confusion the illustration treats as gravest.

- `observed_proportion`:

  Observed proportion of cases in the cell, `frequency / 200`.

- `expected_proportion`:

  Chance-expected proportion of cases in the cell, the product of the
  cell's row (Judge B) and column (Judge A) marginal proportions; the
  parenthetical values in Cohen's Table 1.

## Source

Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision for
scaled disagreement or partial credit. *Psychological Bulletin, 70*(4),
213–220 (Table 1, p. 214).

## Details

Cohen built this table to make a point that is easy to miss: weighted
kappa is fully chance corrected, and it can be *smaller* than unweighted
kappa on the same data. Here the judges disagree far less than chance
expectation in the mildly weighted personality disorder-neurosis cells
but at about the chance level in the heavily weighted neurosis-psychosis
cells, so \\\kappa = .492\\ while \\\kappa_W = .348\\: they disagree
least where it matters least. Interchanging the 6 and 1 weights reverses
the conclusion (\\\kappa_W = .574\\).

The reconstruction was verified against every quantity Cohen computes
from the table: the marginals (.50/.30/.20 for Judge B, .60/.30/.10 for
Judge A), the chance-expected cell proportions, the weighted
disagreement sums \\q'\_o = .90\\ and \\q'\_c = 1.38\\, \\\kappa =
.492\\, \\\kappa_W = .348\\, and his Formula 10 and 13 standard errors
(.0901 and .0916). The
[`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)
help page replicates the full set of analyses, and the weighted kappa
vignette works the illustration end to end, including the orientation of
the printed weight display in the paper's asymmetric-weight validity
reinterpretation.

## References

Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision for
scaled disagreement or partial credit. *Psychological Bulletin, 70*(4),
213–220. [doi:10.1037/h0026256](https://doi.org/10.1037/h0026256)

## See also

[`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
which analyzes this table in its examples.

Other reliability:
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
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
data(diagnosis_agreement)

# Rebuild Cohen's Table 1 layout (Judge B in rows, Judge A in columns).
xtabs(frequency ~ judge_b + judge_a, data = diagnosis_agreement)
#>                       judge_a
#> judge_b                Personality disorder Neurosis Psychosis
#>   Personality disorder                   88       10         2
#>   Neurosis                               14       40         6
#>   Psychosis                              18       10        12

# Unweighted and weighted kappa, reproducing kappa = .492 and
# kappa_w = .348.
tab <- xtabs(frequency ~ judge_b + judge_a, data = diagnosis_agreement)
v   <- xtabs(disagreement_weight ~ judge_b + judge_a,
             data = diagnosis_agreement)
cohen_kappa(table = tab)
#>  weights    kappa se    lower_limit upper_limit z_value p_value  n  
#>  unweighted 0.492 0.051 0.392       0.591       9.64    < 0.0001 200
#>  n_categories
#>  3           
#> 
#> Confidence level: 95%
cohen_kappa(table = tab, weights = unclass(v),
            weight_scaling = "disagreement")
#>  weights             kappa se     lower_limit upper_limit z_value p_value  n  
#>  custom_disagreement 0.348 0.0755 0.2         0.496       4.61    < 0.0001 200
#>  n_categories
#>  3           
#> 
#> Confidence level: 95%

# Cohen's Formula 8 directly from the per-cell quantities.
with(diagnosis_agreement,
     1 - sum(disagreement_weight * observed_proportion) /
         sum(disagreement_weight * expected_proportion))
#> [1] 0.3478261
```
