# Provides the Critical Value(s) for an *F* Distribution

Provides the Critical Value(s) for an *F* Distribution

## Usage

``` r
cv_f(
  alpha_level,
  df_numerator,
  df_denominator,
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

- df_numerator:

  The numerator degrees of freedom (a positive number). In a model
  comparison this is the difference in the number of parameters between
  the two models.

- df_denominator:

  The denominator (error) degrees of freedom (a positive number).

- alternative:

  The type of alternative hypothesis of interest. The default,
  `"greater"`, puts the whole of `alpha_level` in the upper tail, which
  is how the *F* distribution is used to test a model comparison (see
  Details).

- alpha_lower:

  The error rate in the lower tail of the distribution.

- alpha_upper:

  The error rate in the upper tail of the distribution.

- ncp:

  The noncentral parameter (if zero, the default, it is the central *F*
  distribution).

- verbose:

  Provides extra information about areas under the curve.

## Value

Returns the critical value(s), based on the input specifications, in a
output style (a `data.frame` with a row for the lower and the upper
critical value, following the format used by
[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md)).

## Details

Unlike the *t* and *z* distributions, the *F* distribution is not
symmetric and takes only non-negative values, and the usual test of a
model comparison is one-sided: a restricted model fits worse than a full
model, so evidence against the restriction shows up as a large *F*,
never a small one. That is why `alternative` defaults to `"greater"`
here whereas it defaults to `"not_equal"` in
[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md). Maxwell,
Delaney, and Kelley (2027) tabulate these upper-tail values in their
Appendix Table A.2.

Both tails remain available for the situations that need them, such as
an interval for a ratio of variances, either by setting
`alternative = "not_equal"` or by giving `alpha_lower` and `alpha_upper`
directly. When a tail is given zero area its critical value is the
boundary of the support, so `lower_cv` is 0 under the default.

A noncentral parameter can be supplied, which is what a power analysis
needs, though it would not be used for a standard null hypothesis
significance test. See
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
for confidence limits on the noncentral parameter itself.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3, where the *F* test of a model
comparison is developed; Appendix Table A.2 reports these critical
values.)

## See also

[`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_chisq`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
[`cv_bonferroni_f`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md),
[`cv_scheffe`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)

Other critical values:
[`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md),
[`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md),
[`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
[`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The critical value for a model comparison with 3 numerator and 20
# denominator degrees of freedom, at the .05 level.
cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20)
#>  term     value area_less area_greater
#>  lower_cv 0     0         1           
#>  upper_cv 3.1   0.95      0.05        

# An omnibus test of four groups with 24 participants: a - 1 = 3 and
# N - a = 20 degrees of freedom; simple output.
cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20, verbose = FALSE)
#>  term     value
#>  lower_cv 0    
#>  upper_cv 3.1  

# Both tails, as an interval for a ratio of variances would need.
cv_f(alpha_level = .05, df_numerator = 3, df_denominator = 20,
     alternative = "not_equal")
#>  term     value  area_less area_greater
#>  lower_cv 0.0706 0.025     0.975       
#>  upper_cv 3.86   0.975     0.025       
```
