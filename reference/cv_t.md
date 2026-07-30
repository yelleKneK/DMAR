# Provides the Critical Value(s) for a *t*-distribution

Provides the Critical Value(s) for a *t*-distribution

## Usage

``` r
cv_t(
  alpha_level,
  df,
  alternative = "not_equal",
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

  The number of degrees of freedom (a positive number)

- alternative:

  The type of alternative hypothesis of interest.

- alpha_lower:

  The error rate on the lower (negative) side of the distribution.

- alpha_upper:

  The error rate on the upper (positive) side of the distribution.

- ncp:

  The noncentral parameter (if zero, the default, it is the central
  *t*-distribution).

- verbose:

  Provides extra information about areas under the curve.

## Value

Returns the critical value(s), based on the input specifications, in a
output style.

## Details

Though a noncentral parameter can be included, that would not be done
for a standard null hypothesis significance test.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

Other critical values:
[`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md),
[`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md),
[`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A basic call for finding critical values with equal area in the two tails.
cv_t(alpha_level = .05, df = 13)
#>  term     value area_less area_greater
#>  lower_cv -2.16 0.025     0.975       
#>  upper_cv 2.16  0.975     0.025       

# A basic call for a single-sided confidence interval (for "a greater than" alternative hypothesis)
cv_t(alpha_level = .05, df = 13, alternative = "greater")
#>  term     value area_less area_greater
#>  lower_cv -Inf  0         1           
#>  upper_cv 1.77  0.95      0.05        

# A single-sided confidence interval (for "a greater than" alternative hypothesis); simple output.
cv_t(alpha_lower = 0, alpha_upper = .05, df = 13, verbose = FALSE)
#>  term     value
#>  lower_cv -Inf 
#>  upper_cv 1.77 

# For a nonsymmetric 95% confidence interval.
cv_t(alpha_lower = .01, alpha_upper = .04, df = 13)
#>  term     value area_less area_greater
#>  lower_cv -2.65 0.01      0.99        
#>  upper_cv 1.9   0.96      0.04        
```
