# Convert Fisher's *Z* Into the Scale of a Correlation Coefficient (R)

Converts Fisher's *Z* back into the scale of a correlation coefficient
(*r*). Fisher's *Z* is the variance-stabilizing transformation of a
correlation; many authors call it the *z*-prime transform and write the
transformed value as *z*'. The capital *Z* is meaningful: Fisher's *Z*
is not a *z*-score (it is not a standardized variate, that is, an
observation centered and divided by a standard deviation). This function
applies the inverse transform \\r = \mathrm{tanh}(Z)\\ to return to the
scale of a correlation coefficient.

## Usage

``` r
convert_Z_r(Z)

convert_z_r(Z)
```

## Arguments

- Z:

  Fisher's *Z* (the variance-stabilizing transform of a correlation,
  which many authors call *z*')

## Value

A 1-row `data.frame` with columns `term` and `value`. The `term` is
`"r_from_Z"` and `value` is the correlation coefficient corresponding to
the supplied Fisher's *Z*. The inverse direction is
[`convert_r_Z`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md).

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

[`convert_r_Z`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md),
[`ci_cc`](https://yelleknek.github.io/DMAR/reference/ci_cc.md)

Other parameterization conversions:
[`convert_F_chisq()`](https://yelleknek.github.io/DMAR/reference/convert_F_chisq.md),
[`convert_R2`](https://yelleknek.github.io/DMAR/reference/convert_R2.md),
[`convert_cor_cov()`](https://yelleknek.github.io/DMAR/reference/convert_cor_cov.md),
[`convert_d_or()`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md),
[`convert_d_r()`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md),
[`convert_r_Z()`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md),
[`convert_t_smd`](https://yelleknek.github.io/DMAR/reference/convert_t_smd.md),
[`convert_z_normal()`](https://yelleknek.github.io/DMAR/reference/convert_z_normal.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# From Hays (1994, pp. 649--650)
convert_Z_r(0.3654438)
#>  term     value
#>  r_from_Z 0.35 

# convert_z_r() is an exported alias of convert_Z_r()
convert_z_r(0.3654438)
#>  term     value
#>  r_from_Z 0.35 
```
