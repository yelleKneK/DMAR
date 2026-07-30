# Convert Between *F*, \\R^2\\, and Their Noncentral Parameters

Given values of test statistics (and the appropriate additional
information) the value of the noncentral values can be obtained.
Likewise, given noncentral values (and the appropriate additional
information) the value of the test statistic can be obtained.

## Usage

``` r
convert_R2_f(R2 = NULL, df_1 = NULL, df_2 = NULL, p = NULL, N = NULL)

convert_f_R2(F_value = NULL, df_1 = NULL, df_2 = NULL)

convert_lambda_R2(lambda = NULL, N = NULL)

convert_R2_lambda(R2 = NULL, N = NULL)
```

## Arguments

- R2:

  Squared multiple correlation coefficient (population or observed)

- df_1:

  Degrees of freedom for the numerator of the *F*-distribution

- df_2:

  Degrees of freedom for the denominator of the *F*-distribution

- p:

  Number of predictor variables for `R2`

- N:

  Sample size

- F_value:

  The obtained *F* value from a test of significance for the squared
  multiple correlation coefficient

- lambda:

  The noncentral parameter from an *F*-distribution

## Value

Each of the four functions returns a 1-row `data.frame` with columns
`term` and `value`. The `term` entry identifies the conversion performed
(`"r2_f"`, `"f_r2"`, `"lambda_r2"`, or `"r2_lambda"`) and `value` is the
converted scalar. The conversions are exact inverses of one another
(with the appropriate degrees-of-freedom / sample size inputs supplied),
which is what makes them useful inside the noncentrality-parameter
confidence interval machinery of
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md) and
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md).

## Details

These functions are especially helpful in the search for confidence
intervals for noncentral parameters, as they convert to and from related
quantities.

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

## See also

[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md),
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)

Other parameterization conversions:
[`convert_F_chisq()`](https://yelleknek.github.io/DMAR/reference/convert_F_chisq.md),
[`convert_Z_r()`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md),
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
convert_R2_lambda(R2 = .5, N = 100)
#>  term      value
#>  r2_lambda 100  
```
