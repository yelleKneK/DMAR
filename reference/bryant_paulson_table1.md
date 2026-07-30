# Bryant & Paulson (1976) Table 1: critical values of the generalized studentized range

The complete Table 1 of Bryant and Paulson (1976) — the upper 5% and 1%
points \\q\_{\alpha;p,k,\nu}\\ of the generalized studentized range
\\Q_p\\. These are the *original* single-step (Tukey-type) critical
values for simultaneous confidence intervals on, and tests of, contrasts
of covariate-adjusted means in the analysis of covariance when the
covariates are random (Bryant and Paulson, 1976, Equations 7–8). They
are the values returned by
[`cv_bryant_paulson`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md)
(with `procedure = "tukey"`) and
[`qbryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md),
supplied here in tidy long form for convenient lookup.

## Usage

``` r
bryant_paulson_table1
```

## Format

A data frame with 1188 rows (3 covariate counts \\\times\\ 2 error rates
\\\times\\ 18 \\\nu\\ \\\times\\ 11 \\k\\) and 5 variables.

- `p`:

  Number of random covariates (1, 2, or 3).

- `alpha`:

  Upper-tail error rate (`0.05` = Table 1a, `0.01` = Table 1b).

- `nu`:

  Error degrees of freedom \\\nu\\ (one of 2, 3, 4, 5, 6, 7, 8, 10, 12,
  14, 16, 18, 20, 24, 30, 40, 60, 120).

- `k`:

  Number of groups whose adjusted means are compared (one of 2, 3, 4, 5,
  6, 7, 8, 10, 12, 16, 20).

- `critical_value`:

  The critical value \\q\_{\alpha;p,k,\nu}\\, on the studentized-range
  scale, to two decimal places.

## Source

Bryant, J. L., & Paulson, A. S. (1976). An extension of Tukey's method
of multiple comparisons to experimental designs with random concomitant
variables. *Biometrika, 63*(3), 631–638 (Table 1).
[doi:10.1093/biomet/63.3.631](https://doi.org/10.1093/biomet/63.3.631)

## Details

**Provenance of the values.** The numbers are the Bryant–Paulson
critical values *computed* from the distribution (Bryant and Paulson,
1976, Eq. 17; see
[`qbryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md))
and rounded to two decimals on the exact grid of the published Table 1.
They reproduce the printed table to that paper's own stated accuracy of
\\\pm 0.01\\: a cell-by-cell comparison against the scanned 1976 article
agrees for 94.9% of the cleanly readable cells, and the only systematic
departures are at the extreme \\\nu = 2\\, \\\alpha = .01\\ corner,
where the 1976 *numerical tabulation* (not this exact computation) is
the less precise of the two. Computed values are shipped, rather than an
error-prone optical transcription of the half-century-old scan, so that
the table is internally exact and consistent with
[`cv_bryant_paulson`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md).
The package test `test-bryant_paulson_table1.R` checks this data set
both against the function and against a set of anchor cells read by hand
from the printed table.

**Relationship to the other table.** Bryant and Paulson (1976) tabulate
the single-step values here; Bryant and Bruvold (1980) later tabulated
the *stepwise* Duncan multiple-range significant ranges derived from
them (shipped as
[`bryant_bruvold_table2`](https://yelleknek.github.io/DMAR/reference/bryant_bruvold_table2.md))
and showed the same distribution holds when the covariates are not
identically distributed across groups.

## References

Bryant, J. L., & Bruvold, N. T. (1980). Multiple comparison procedures
in the analysis of covariance. *Journal of the American Statistical
Association, 75*(372), 874–880.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

## See also

[`cv_bryant_paulson`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md)
and
[`qbryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)
for computing these values,
[`bryant_bruvold_table2`](https://yelleknek.github.io/DMAR/reference/bryant_bruvold_table2.md)
for the Duncan multiple-range table, and
[`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)
for the confidence intervals they produce.

## Author

Ken Kelley

## Examples

``` r
data(bryant_paulson_table1)
str(bryant_paulson_table1)
#> 'data.frame':    1188 obs. of  5 variables:
#>  $ p             : num  1 1 1 1 1 1 1 1 1 1 ...
#>  $ alpha         : num  0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 ...
#>  $ nu            : num  2 2 2 2 2 2 2 2 2 2 ...
#>  $ k             : num  2 3 4 5 6 7 8 10 12 16 ...
#>  $ critical_value: num  18.9 25.9 30.8 34.5 37.4 ...

# The value Bryant & Bruvold (1980) use in their worked example:
subset(bryant_paulson_table1, p == 1 & alpha == 0.05 & nu == 14 & k == 6)
#>     p alpha nu k critical_value
#> 302 1  0.05 14 6           4.83

# cv_bryant_paulson() returns the same value.
cv_bryant_paulson(0.05, df = 14, groups = 6, covariates = 1,
                  verbose = FALSE)$value
#> [1] 4.829856

# Reshape one sub-table back to the wide form printed in the article.
with(subset(bryant_paulson_table1, p == 1 & alpha == 0.05),
     tapply(critical_value, list(nu = nu, k = k), identity))
#>      k
#> nu       2     3     4     5     6     7     8    10    12    16    20
#>   2   7.95 11.00 12.99 14.46 15.61 16.55 17.35 18.64 19.67 21.25 22.43
#>   3   5.42  7.18  8.32  9.17  9.84 10.39 10.86 11.62 12.22 13.14 13.83
#>   4   4.51  5.84  6.69  7.32  7.82  8.23  8.58  9.15  9.61 10.30 10.82
#>   5   4.06  5.17  5.88  6.40  6.82  7.16  7.45  7.93  8.30  8.88  9.32
#>   6   3.79  4.78  5.40  5.86  6.23  6.53  6.78  7.20  7.53  8.04  8.43
#>   7   3.62  4.52  5.09  5.51  5.84  6.11  6.34  6.72  7.03  7.49  7.84
#>   8   3.49  4.34  4.87  5.26  5.57  5.82  6.03  6.39  6.67  7.10  7.43
#>   10  3.32  4.10  4.58  4.93  5.21  5.43  5.63  5.94  6.19  6.58  6.87
#>   12  3.22  3.95  4.40  4.73  4.98  5.19  5.37  5.67  5.90  6.26  6.53
#>   14  3.15  3.85  4.28  4.59  4.83  5.03  5.20  5.48  5.70  6.03  6.29
#>   16  3.10  3.77  4.19  4.49  4.72  4.91  5.07  5.34  5.55  5.87  6.12
#>   18  3.06  3.72  4.12  4.41  4.63  4.82  4.98  5.23  5.44  5.75  5.98
#>   20  3.03  3.67  4.07  4.35  4.57  4.75  4.90  5.15  5.35  5.65  5.88
#>   24  2.98  3.61  3.99  4.26  4.47  4.65  4.79  5.03  5.22  5.51  5.73
#>   30  2.94  3.55  3.91  4.18  4.38  4.54  4.69  4.91  5.09  5.37  5.58
#>   40  2.89  3.49  3.84  4.09  4.29  4.45  4.58  4.80  4.97  5.23  5.43
#>   60  2.85  3.43  3.77  4.01  4.20  4.35  4.48  4.69  4.85  5.10  5.29
#>   120 2.81  3.37  3.70  3.93  4.11  4.26  4.38  4.58  4.73  4.97  5.15
```
