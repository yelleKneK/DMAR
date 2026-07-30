# Bryant & Bruvold (1980) Table 2: Duncan multiple-range critical values for ANCOVA

The complete Table 2 of Bryant and Bruvold (1980) — the “significant
ranges” for Duncan's multiple-range test extended to the analysis of
covariance with random covariates. These are the stepwise critical
values \\r\_{\alpha;p,k,\nu}\\ of the Bryant–Paulson generalized
studentized range (see
[`cv_bryant_paulson`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md)
with `procedure = "duncan"`). The data set is provided in tidy long form
for convenient lookup and as a reference against which
[`cv_bryant_paulson`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md)
is validated (it reproduces every cell; see the package tests).

## Usage

``` r
bryant_bruvold_table2
```

## Format

A data frame with 948 rows (6 sub-tables \\\times\\ their tabled cells)
and 5 variables.

- `p`:

  Number of random covariates (1, 2, or 3).

- `alpha`:

  Per-step error rate of the multiple-range test (`0.05` or `0.01`).

- `nu`:

  Error degrees of freedom \\\nu\\ (one of 2, 3, 4, 5, 6, 7, 8, 9, 10,
  12, 14, 16, 18, 20, 24, 30, 40, 60, 120).

- `k`:

  Number of groups (steps) spanned by the range (one of 2, 3, 4, 5, 6,
  7, 8, 10, 12, 16, 20).

- `critical_value`:

  The tabled significant range \\r\_{\alpha;p,k,\nu}\\, on the
  studentized-range scale, exactly as printed in the article (two
  decimal places).

## Source

Bryant, J. L., & Bruvold, N. T. (1980). Multiple comparison procedures
in the analysis of covariance. *Journal of the American Statistical
Association, 75*(372), 874–880 (Table 2).
[doi:10.2307/2287175](https://doi.org/10.2307/2287175)

## Details

Bryant and Bruvold's Table 2 is laid out as six sub-tables (covariate
count \\p \in \\1,2,3\\\\ crossed with \\\alpha \in \\.05, .01\\\\),
each a grid of \\\nu\\ (rows) by \\k\\ (columns); not every \\(\nu, k)\\
cell is printed for small \\\nu\\. The values stored here are
transcribed verbatim from the article.
`data-raw/bryant_bruvold_table2.R` holds the six wide matrices exactly
as printed and assembles this long data frame; the construction is
self-contained (no source PDF or digitization step is required).

The Duncan significant range uses variable protection levels \\\alpha_k
= 1 - (1-\alpha)^{k-1}\\ and is built up from the single-step
Bryant–Paulson critical values \\q\_{\alpha;p,k,\nu}\\ by
\\r\_{\alpha;p,2,\nu} = q\_{\alpha;p,2,\nu}\\ and \\r\_{\alpha;p,k,\nu}
= \max\\r\_{\alpha;p,k-1,\nu}, q\_{\alpha_k;p,k,\nu}\\\\;
[`cv_bryant_paulson`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md)
computes exactly these values and reproduces the table to its
two-decimal precision (the only cells off by more than 0.01 are the
extreme \\\nu = 2\\, \\\alpha = .01\\ entries, where the paper's own
accuracy is one unit in the last decimal).

## References

Bryant, J. L., & Paulson, A. S. (1976). An extension of Tukey's method
of multiple comparisons to experimental designs with random concomitant
variables. *Biometrika, 63*, 631–638.

Duncan, D. B. (1955). Multiple range and multiple F tests. *Biometrics,
11*, 1–42.

## See also

[`cv_bryant_paulson`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md)
for computing these critical values,
[`bryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)
for the underlying distribution, and
[`test_market`](https://yelleknek.github.io/DMAR/reference/test_market.md)
for the worked-example data set.

## Author

Ken Kelley

## Examples

``` r
data(bryant_bruvold_table2)
str(bryant_bruvold_table2)
#> 'data.frame':    948 obs. of  5 variables:
#>  $ p             : num  1 1 1 1 1 1 1 1 1 1 ...
#>  $ alpha         : num  0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 ...
#>  $ nu            : num  2 2 3 3 3 4 4 4 4 5 ...
#>  $ k             : num  2 3 2 3 4 2 3 4 5 2 ...
#>  $ critical_value: num  19.1 19.1 10.3 10.3 10.3 ...

# The cell used in the paper's worked example: q_.05;1,6,14-style range.
subset(bryant_bruvold_table2, p == 1 & alpha == 0.05 & nu == 14 & k == 6)
#>     p alpha nu k critical_value
#> 224 1  0.05 14 6            3.5

# cv_bryant_paulson() reproduces a tabled value.
subset(bryant_bruvold_table2, p == 1 & alpha == 0.05 & nu == 20 & k == 2)
#>     p alpha nu k critical_value
#> 251 1  0.05 20 2           3.03
cv_bryant_paulson(0.05, df = 20, groups = 2, covariates = 1,
                  procedure = "duncan", verbose = FALSE)$value
#> [1] 3.027767

# Reshape one sub-table back to the wide form printed in the article.
with(subset(bryant_bruvold_table2, p == 1 & alpha == 0.05),
     tapply(critical_value, list(nu = nu, k = k), identity))
#>      k
#> nu       2    3    4    5    6    7    8   10   12   16   20
#>   2   7.96 7.96   NA   NA   NA   NA   NA   NA   NA   NA   NA
#>   3   5.42 5.42 5.42   NA   NA   NA   NA   NA   NA   NA   NA
#>   4   4.51 4.59 4.60 4.60   NA   NA   NA   NA   NA   NA   NA
#>   5   4.06 4.18 4.22 4.23 4.23   NA   NA   NA   NA   NA   NA
#>   6   3.79 3.92 3.98 4.01 4.02 4.02   NA   NA   NA   NA   NA
#>   7   3.62 3.75 3.82 3.86 3.88 3.89 3.89   NA   NA   NA   NA
#>   8   3.49 3.63 3.71 3.76 3.78 3.80 3.80 3.80   NA   NA   NA
#>   9   3.40 3.54 3.62 3.67 3.71 3.72 3.74 3.74   NA   NA   NA
#>   10  3.32 3.47 3.56 3.61 3.65 3.67 3.68 3.70 3.70   NA   NA
#>   12  3.22 3.37 3.46 3.52 3.56 3.59 3.61 3.63 3.64   NA   NA
#>   14  3.15 3.30 3.39 3.45 3.50 3.53 3.55 3.58 3.60 3.60   NA
#>   16  3.10 3.25 3.34 3.40 3.45 3.48 3.51 3.54 3.56 3.58   NA
#>   18  3.06 3.21 3.30 3.37 3.42 3.45 3.48 3.52 3.54 3.56 3.57
#>   20  3.03 3.18 3.27 3.34 3.39 3.42 3.45 3.49 3.52 3.55 3.56
#>   24  2.98 3.13 3.23 3.30 3.35 3.38 3.42 3.46 3.49 3.53 3.54
#>   30  2.94 3.09 3.18 3.25 3.31 3.35 3.38 3.43 3.46 3.50 3.53
#>   40  2.89 3.04 3.14 3.21 3.26 3.31 3.34 3.39 3.43 3.48 3.51
#>   60  2.85 3.00 3.10 3.17 3.22 3.27 3.30 3.36 3.40 3.46 3.50
#>   120 2.81 2.96 3.06 3.13 3.19 3.23 3.27 3.33 3.37 3.44 3.48
```
