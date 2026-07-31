# Convert a Correlation Coefficient (R) Into the Scale of Fisher's *Z*

This function converts a correlation coefficient into the scale of
Fisher's *Z*, the variance-stabilizing transformation of a correlation.
Many authors call this map the *z*-prime transform and write the
transformed value as *z*'. The capital *Z* is meaningful: Fisher's *Z*
is not a *z*-score (it is not a standardized variate, that is, an
observation centered and divided by a standard deviation). It is the
transform \\Z = \mathrm{atanh}(r)\\ of a correlation coefficient,
applied because the sampling distribution of *Z* is approximately normal
with a variance that does not depend on the population correlation,
which makes *Z* convenient for forming confidence intervals.

## Usage

``` r
convert_r_Z(r)

convert_r_z(r)
```

## Arguments

- r:

  Correlation coefficient (between two variables)

## Value

A 1-row `data.frame` with columns `term` and `value`. The `term` is
`"Z_from_r"` and `value` is Fisher's *Z* corresponding to the supplied
correlation coefficient. The inverse direction is
[`convert_Z_r`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md).

## Details

This function is typically used in the context of forming a confidence
interval for a population correlation coefficient. Note that, in that
situation, the two variables are assumed to follow a bivariate normal
distribution (e.g., Hays, 1994).

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Hays, W. L. (1994). *Statistics* (5th ed.). Fort Worth, TX: Harcourt
Brace College Publishers.

## See also

[`convert_Z_r`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md),
[`ci_cc`](https://yelleknek.github.io/DMAR/reference/ci_cc.md)

Other parameterization conversions:
[`convert_F_chisq()`](https://yelleknek.github.io/DMAR/reference/convert_F_chisq.md),
[`convert_R2`](https://yelleknek.github.io/DMAR/reference/convert_R2.md),
[`convert_Z_r()`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md),
[`convert_cor_cov()`](https://yelleknek.github.io/DMAR/reference/convert_cor_cov.md),
[`convert_d_or()`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md),
[`convert_d_r()`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md),
[`convert_t_smd`](https://yelleknek.github.io/DMAR/reference/convert_t_smd.md),
[`convert_z_normal()`](https://yelleknek.github.io/DMAR/reference/convert_z_normal.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# From Hays (1994, pp. 649--650)
convert_r_Z(.35)
#>  term     value
#>  Z_from_r 0.365

# convert_r_z() is an exported alias of convert_r_Z()
convert_r_z(.35)
#>  term     value
#>  Z_from_r 0.365
```
