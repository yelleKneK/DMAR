# Convert Between the Standardized Mean Difference and the Correlation

Invertible conversions between a two-group standardized mean difference
(Cohen's *d*) and the (point-biserial) correlation between the outcome
and group membership. `convert_d_r()` maps *d* to *r*; `convert_r_d()`
maps *r* back to *d*. These are the standard conversions used to bring
effect sizes reported in different metrics onto a common scale, for
example when synthesizing a literature in which some studies report mean
differences and others report correlations (Borenstein, Hedges, Higgins,
& Rothstein, 2009, Chapter 7).

## Usage

``` r
convert_d_r(d, n_1 = NULL, n_2 = NULL)

convert_r_d(r, n_1 = NULL, n_2 = NULL)
```

## Arguments

- d:

  The standardized mean difference.

- n_1, n_2:

  Optional per-group sample sizes. When supplied, the conversion uses
  the unequal-group factor \\a = (n_1 + n_2)^2 / (n_1 n_2)\\; when
  omitted, equal group sizes are assumed, for which \\a = 4\\.

- r:

  The point-biserial correlation, in \\(-1, 1)\\.

## Value

A `data.frame` (class `dmar_tbl`) with a single row: term `r` (for
`convert_d_r`) or `smd` (for `convert_r_d`) and its `value`.

## Details

With \\a = (n_1 + n_2)^2/(n_1 n_2)\\ (equal to 4 for equal groups), the
two directions are \$\$r = \frac{d}{\sqrt{d^2 + a}}, \qquad d =
\frac{\sqrt{a}\\ r}{\sqrt{1 - r^2}},\$\$ exact inverses of one another
for a given \\a\\. The same `n_1` and `n_2` must be supplied to both
directions for the round trip to be exact.

## References

Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R.
(2009). *Introduction to meta-analysis*. Wiley.

## See also

[`convert_d_or`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md)
/
[`convert_or_d`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md)
for the odds ratio leg of the same triangle;
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) and
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
for estimating the quantities being converted.

Other parameterization conversions:
[`convert_F_chisq()`](https://yelleknek.github.io/DMAR/reference/convert_F_chisq.md),
[`convert_R2`](https://yelleknek.github.io/DMAR/reference/convert_R2.md),
[`convert_Z_r()`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md),
[`convert_cor_cov()`](https://yelleknek.github.io/DMAR/reference/convert_cor_cov.md),
[`convert_d_or()`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md),
[`convert_r_Z()`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md),
[`convert_t_smd`](https://yelleknek.github.io/DMAR/reference/convert_t_smd.md),
[`convert_z_normal()`](https://yelleknek.github.io/DMAR/reference/convert_z_normal.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Equal groups: d = 0.5 corresponds to r about .243.
convert_d_r(d = 0.5)
#>  term value
#>  r    0.243

# And back, exactly.
convert_r_d(r = convert_d_r(d = 0.5)$value)
#>  term value
#>  smd  0.5  

# Unequal groups change the conversion factor.
convert_d_r(d = 0.5, n_1 = 20, n_2 = 80)
#>  term value
#>  r    0.196
```
