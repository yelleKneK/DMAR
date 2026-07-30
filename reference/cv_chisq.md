# Provides the Critical Value(s) for a Chi Square Distribution

Provides the Critical Value(s) for a Chi Square Distribution

## Usage

``` r
cv_chisq(
  alpha_level,
  df,
  alternative = "greater",
  alpha_lower,
  alpha_upper,
  ncp = 0,
  verbose = TRUE
)
```

## Arguments

- alpha_level:

  Type I error rate (i.e., the false positive rate).

- df:

  The number of degrees of freedom (a positive number).

- alternative:

  The type of alternative hypothesis of interest. The default,
  `"greater"`, puts the whole of `alpha_level` in the upper tail, which
  is how the chi square distribution is used to test a model or an
  association (see Details).

- alpha_lower:

  The error rate in the lower tail of the distribution.

- alpha_upper:

  The error rate in the upper tail of the distribution.

- ncp:

  The noncentral parameter (if zero, the default, it is the central chi
  square distribution).

- verbose:

  Provides extra information about areas under the curve.

## Value

Returns the critical value(s), based on the input specifications, in a
output style (a `data.frame` with a row for the lower and the upper
critical value, following the format used by
[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md)).

## Details

Like the *F* distribution and unlike *t* and *z*, the chi square
distribution is not symmetric and takes only non-negative values. Its
common uses are one-sided in the upper tail: a test of association in a
contingency table, a likelihood ratio test, and a test of model fit all
reject for large values, because a poorly fitting model produces a large
discrepancy, never a small one. That is why `alternative` defaults to
`"greater"` here whereas it defaults to `"not_equal"` in
[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md). Maxwell,
Delaney, and Kelley (2027) tabulate these upper-tail values in their
Appendix Table A.9.

Both tails remain available for the situations that need them, most
commonly an interval for a variance, which uses an upper and a lower chi
square quantile. Set `alternative = "not_equal"`, or give `alpha_lower`
and `alpha_upper` directly. When a tail is given zero area its critical
value is the boundary of the support, so `lower_cv` is 0 under the
default.

A noncentral parameter can be supplied, which is what a power analysis
for a test of model fit needs, though it would not be used for a
standard null hypothesis significance test. See
[`conf_limits_nc_chisq`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md)
for confidence limits on the noncentral parameter itself.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (Appendix Table A.9 reports these critical values.)

## See also

[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_f`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`conf_limits_nc_chisq`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md)

Other critical values:
[`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md),
[`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md),
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
[`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The critical value for a test on 3 degrees of freedom at the .05 level.
cv_chisq(alpha_level = .05, df = 3)
#>  term     value area_less area_greater
#>  lower_cv 0     0         1           
#>  upper_cv 7.81  0.95      0.05        

# Simple output.
cv_chisq(alpha_level = .05, df = 3, verbose = FALSE)
#>  term     value
#>  lower_cv 0    
#>  upper_cv 7.81 

# Both tails, as an interval for a variance would need.
cv_chisq(alpha_level = .05, df = 10, alternative = "not_equal")
#>  term     value area_less area_greater
#>  lower_cv 3.25  0.025     0.975       
#>  upper_cv 20.5  0.975     0.025       
```
