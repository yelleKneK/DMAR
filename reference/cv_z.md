# Provides the Critical Value(s) for the Standard Normal Distribution (the *z*-distribution, With Mean 0 and Variance 1)

Provides the Critical Value(s) for the Standard Normal Distribution (the
*z*-distribution, With Mean 0 and Variance 1)

## Usage

``` r
cv_z(
  alpha_level,
  alternative = "not_equal",
  alpha_lower,
  alpha_upper,
  verbose = TRUE
)
```

## Arguments

- alpha_level:

  Type I error rate (i.e., the false positive rate).

- alternative:

  The type of alternative hypothesis of interest.

- alpha_lower:

  The error rate on the lower (negative) side of the distribution.

- alpha_upper:

  The error rate on the upper (positive) side of the distribution.

- verbose:

  Provides extra information about areas under the curve.

## Value

Returns the critical value(s), based on the input specifications, in a
output style.

## See also

Other critical values:
[`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md),
[`cv_bryant_paulson()`](https://yelleknek.github.io/DMAR/reference/cv_bryant_paulson.md),
[`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
[`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A basic call for finding critical values with equal area in the two tails.
cv_z(alpha_level = .05)
#>  term     value area_less area_greater
#>  lower_cv -1.96 0.025     0.975       
#>  upper_cv 1.96  0.975     0.025       

# A basic call for a single-sided confidence interval (for "a greater than" alternative hypothesis)
cv_z(alpha_level = .05, alternative = "greater")
#>  term     value area_less area_greater
#>  lower_cv -Inf  0         1           
#>  upper_cv 1.64  0.95      0.05        

# A single-sided confidence interval (for "a greater than" alternative hypothesis); simple output.
cv_z(alpha_lower = 0, alpha_upper = .05, verbose = FALSE)
#>  term     value
#>  lower_cv -Inf 
#>  upper_cv 1.64 

# For a nonsymmetric 95% confidence interval.
cv_z(alpha_lower = .01, alpha_upper = .04)
#>  term     value area_less area_greater
#>  lower_cv -2.33 0.01      0.99        
#>  upper_cv 1.75  0.96      0.04        
```
