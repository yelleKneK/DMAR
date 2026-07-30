# Intraclass Correlation Coefficients With Confidence Intervals

Computes one or more of the six standard intraclass correlation
coefficients (ICC) of Shrout and Fleiss (1979) for an \\n \times k\\
matrix of ratings (rows = subjects, columns = raters/measurements),
along with the *F*-distribution-based confidence interval for each.

## Usage

``` r
icc(x, type = "ICC(2,1)", conf_level = 0.95)
```

## Arguments

- x:

  An \\n \times k\\ numeric matrix or `data.frame`: each row is a
  subject (target) and each column is a rater (or repeated measurement).
  Missing values are not supported; remove or impute first.

- type:

  Which ICC variant(s) to return. Aliases include `"1"`, `"2"`, `"3"`
  (single rater versions of types 1, 2, 3) and `"1k"`, `"2k"`, `"3k"`
  (average-of-\\k\\-raters versions). The full names `"ICC(1,1)"`,
  `"ICC(2,1)"`, `"ICC(3,1)"`, `"ICC(1,k)"`, `"ICC(2,k)"`, `"ICC(3,k)"`
  are also accepted, as is `"all"`. Vectors are accepted; a row is
  returned for each type.

- conf_level:

  Confidence level for the interval (default `0.95`).

## Value

A `data.frame` (class `dmar_tbl`) with one row per requested ICC type
and columns `type`, `value` (point estimate), `lower_limit`,
`upper_limit`, `F_value`, `df_1`, `df_2`, and `p_value` (for the
implicit null \\H_0\\: \mathrm{ICC} = 0\\).

## Details

The six variants follow Shrout and Fleiss's (1979) classification:

- `ICC(1,1)`: one-way random model, single rater. Each subject is rated
  by a different (random) rater; appropriate when raters are not crossed
  with subjects.

- `ICC(1,k)`: one-way random model, average of \\k\\ raters.

- `ICC(2,1)`: two-way random model, single rater, absolute agreement.
  Subjects \\\times\\ raters fully crossed; both effects random.

- `ICC(2,k)`: two-way random model, average of \\k\\ raters, absolute
  agreement.

- `ICC(3,1)`: two-way mixed model, single rater, consistency. Raters are
  fixed (the only ones of interest); rates differences in mean across
  raters are not penalized.

- `ICC(3,k)`: two-way mixed model, average of \\k\\ raters, consistency.

Confidence intervals follow the *F*-distribution-based formulas of
Shrout and Fleiss (1979, pp.\\ 425–426); the `ICC(2,1)` interval uses
their asymmetric approximate-*df* formulation.

## References

Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
assessing rater reliability. *Psychological Bulletin, 86*(2), 420–428.

McGraw, K. O., & Wong, S. P. (1996). Forming inferences about some
intraclass correlation coefficients. *Psychological Methods, 1*(1),
30–46.
[doi:10.1037/1082-989X.1.1.30](https://doi.org/10.1037/1082-989X.1.1.30)

## See also

[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_r.md),
[`descriptives`](https://yelleknek.github.io/DMAR/reference/descriptives.md)

Other reliability:
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
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
# Shrout & Fleiss (1979), Table 1: 4 raters, 6 targets.
shrout_fleiss <- matrix(
  c(9, 2, 5, 8,
    6, 1, 3, 2,
    8, 4, 6, 8,
    7, 1, 2, 6,
    10, 5, 6, 9,
    6, 2, 4, 7),
  nrow = 6, byrow = TRUE,
  dimnames = list(paste0("Target_", 1:6),
                  paste0("Judge_",  1:4))
)

icc(shrout_fleiss, type = "all")
#>  type     value lower_limit upper_limit F_value df_1 df_2 p_value
#>  ICC(1,1) 0.166 -0.133      0.723       1.79    5    18   0.1648 
#>  ICC(2,1) 0.29  0.0188      0.761       11      5    15   0.0001 
#>  ICC(3,1) 0.715 0.342       0.946       11      5    15   0.0001 
#>  ICC(1,k) 0.443 -0.884      0.912       1.79    5    18   0.1648 
#>  ICC(2,k) 0.62  0.0711      0.927       11      5    15   0.0001 
#>  ICC(3,k) 0.909 0.676       0.986       11      5    15   0.0001 
#> 
#> Confidence level: 95%

# Just the two-way mixed-model single-rater consistency ICC:
icc(shrout_fleiss, type = "ICC(3,1)")
#>  type     value lower_limit upper_limit F_value df_1 df_2 p_value
#>  ICC(3,1) 0.715 0.342       0.946       11      5    15   0.0001 
#> 
#> Confidence level: 95%
```
