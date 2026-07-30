# Conversion Functions for Noncentral *t*-distribution

Functions useful for converting a standardized mean difference to a
noncentrality parameter, and vice versa.

## Usage

``` r
convert_delta_lambda(delta, n_1, n_2)

convert_lambda_delta(lambda, n_1, n_2)
```

## Arguments

- delta:

  Population value of the standardized mean difference

- n_1:

  Sample size in group 1

- n_2:

  Sample size in group 2

- lambda:

  noncentral value from a *t*-distribution

## Value

Each function returns a 1-row `data.frame` with columns `term` and
`value`. The `term` entry identifies the conversion (`"delta_lambda"` or
`"lambda_delta"`) and `value` is the converted scalar. The two functions
are exact inverses given the *per-group* sample sizes `n_1` and `n_2`.

## Details

Although `lambda` is the population noncentral value, an estimate of it
is the observed value of a *t*-statistic. Likewise, delta can be
estimated as the observed standardized mean difference. Thus, the
observed standardized mean difference can be converted to the observed
*t*-value. These functions are especially helpful in the context of
forming confidence intervals for the population standardized mean
difference.

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)

Other parameterization conversions:
[`convert_F_chisq()`](https://yelleknek.github.io/DMAR/reference/convert_F_chisq.md),
[`convert_R2`](https://yelleknek.github.io/DMAR/reference/convert_R2.md),
[`convert_Z_r()`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md),
[`convert_cor_cov()`](https://yelleknek.github.io/DMAR/reference/convert_cor_cov.md),
[`convert_d_or()`](https://yelleknek.github.io/DMAR/reference/convert_d_or.md),
[`convert_d_r()`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md),
[`convert_r_Z()`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md),
[`convert_z_normal()`](https://yelleknek.github.io/DMAR/reference/convert_z_normal.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
convert_lambda_delta(lambda = 2, n_1 = 113, n_2 = 113)
#>  term         value
#>  lambda_delta 0.266
convert_delta_lambda(delta = .266076, n_1 = 113, n_2 = 113)
#>  term         value
#>  delta_lambda 2    
```
