# Limits of Agreement (Bland-Altman) With Confidence Intervals on the Limits

Computes the limits of agreement (LoA) of Bland and Altman (1986, 1999)
between two methods of measurement applied to the same units, along with
the Carkeet (2015) exact confidence intervals on the LoA themselves. The
exact CIs are based on the noncentral *t* distribution of the sample LoA
and replace the approximate normal CIs originally given by Bland &
Altman (1999).

## Usage

``` r
loa(x, y, coverage = 0.95, conf_level = 0.95)

bland_altman_loa(x, y, coverage = 0.95, conf_level = 0.95)
```

## Arguments

- x, y:

  Paired numeric vectors of equal length (e.g., method A and method B
  applied to the same units).

- coverage:

  Probability content of the limits of agreement. Default `0.95` (the
  conventional 95 is bracketed by \\\pm z\_{(1 + \text{coverage})/2}\\
  standard deviations of the differences.

- conf_level:

  Confidence level for the CIs on the LoA themselves. Default `0.95`.

## Value

A `data.frame` with rows for the mean difference, the SD of differences,
the lower and upper LoA (`loa_lower`, `loa_upper`), and the lower /
upper CI bounds on each LoA.

## Details

**Definition.** For paired observations \\(x_i, y_i)\\, the Bland-Altman
limits of agreement are \$\$\mathrm{LoA}\_\pm \\=\\ \bar d \pm k \cdot
s_d,\$\$ where \\d_i = y_i - x_i\\, \\\bar d\\ is the mean of the
differences, \\s_d\\ is their SD, and \\k = z\_{(1 +
\mathrm{coverage})/2}\\ (for 95% coverage, \\k = 1.96\\). The LoA are
population intervals: they describe the range within which approximately
*coverage*% of *individual* differences are expected to lie if the
differences are normally distributed.

**CIs on the LoA themselves.** The sample LoA are random variables;
their sampling distribution depends on the noncentral- *t*. Carkeet
(2015) derived the exact distribution and the resulting CIs:
\$\$\mathrm{LoA}\_\pm \pm \mathrm{SE}\_{\mathrm{LoA}} \cdot t\_{1 -
\alpha/2, n - 1},\$\$ with the noncentral *t*-based SE replacing the
normal- approximation form (Bland & Altman, 1999) which is biased
downward at small \\n\\. The implementation uses the closed-form
noncentral- *t* CI from Carkeet's equation 6.

**Caveats.** The LoA construction assumes (i) the differences \\d_i\\
are approximately normally distributed, and (ii) the difference does not
systematically depend on the magnitude of the measurement (proportional
bias). Both should be checked, the second by plotting \\d_i\\ against
\\(x_i + y_i)/2\\; a non-flat relationship indicates that a single set
of LoA is inappropriate.

## References

Bland, J. M., & Altman, D. G. (1986). Statistical methods for assessing
agreement between two methods of clinical measurement. *Lancet,
327*(8476), 307–310.

Bland, J. M., & Altman, D. G. (1999). Measuring agreement in method
comparison studies. *Statistical Methods in Medical Research, 8*(2),
135–160.
[doi:10.1191/096228099673819272](https://doi.org/10.1191/096228099673819272)

Carkeet, A. (2015). Exact parametric confidence intervals for
Bland-Altman limits of agreement. *Optometry and Vision Science, 92*(3),
e71–e80.
[doi:10.1097/OPX.0000000000000513](https://doi.org/10.1097/OPX.0000000000000513)

## See also

[`lin_ccc`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md)

Other agreement and measurement:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`content_validity_index()`](https://yelleknek.github.io/DMAR/reference/content_validity_index.md),
[`gwet_ac()`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`krippendorff_alpha()`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md),
[`lin_ccc()`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md),
[`variance_components_mls()`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Two methods that agree well:
set.seed(113)
method_a <- rnorm(40, mean = 100, sd = 15)
method_b <- method_a + rnorm(40, mean = 0, sd = 3)
loa(method_a, method_b)
#>  term                  value
#>  mean_difference       0.259
#>  sd_difference         3.37 
#>  loa_lower             -6.34
#>  loa_lower_lower_limit -8.56
#>  loa_lower_upper_limit -4.82
#>  loa_upper             6.86 
#>  loa_upper_lower_limit 5.34 
#>  loa_upper_upper_limit 9.08 
#> 
#> Confidence level: 95%

# 2. 90% LoA with 95% CIs on the limits:
loa(method_a, method_b, coverage = 0.90, conf_level = 0.95)
#>  term                  value
#>  mean_difference       0.259
#>  sd_difference         3.37 
#>  loa_lower             -5.28
#>  loa_lower_lower_limit -7.25
#>  loa_lower_upper_limit -3.9 
#>  loa_upper             5.8  
#>  loa_upper_lower_limit 4.42 
#>  loa_upper_upper_limit 7.77 
#> 
#> Confidence level: 95%
```
