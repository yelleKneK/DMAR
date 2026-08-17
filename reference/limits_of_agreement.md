# Limits of Agreement (Bland-Altman) With Confidence Intervals on the Limits

Computes the limits of agreement (LoA) of Bland and Altman (1986, 1999)
between two methods of measurement applied to the same units, along with
the Carkeet (2015) exact confidence intervals on the LoA themselves. The
CIs treat the two limits either as a pair (the default), so that the
confidence statement holds for both limits jointly, or individually, one
limit at a time. Both constructions replace the approximate normal CIs
originally given by Bland and Altman (1999), which are too narrow at
small *n*. The short alias `loa()` calls the same function.

## Usage

``` r
limits_of_agreement(
  x,
  y,
  coverage = 0.95,
  conf_level = 0.95,
  method = c("pair", "individual")
)

loa(x, y, coverage = 0.95, conf_level = 0.95, method = c("pair", "individual"))
```

## Arguments

- x, y:

  Paired numeric vectors of equal length (e.g., method A and method B
  applied to the same units).

- coverage:

  Probability content of the limits of agreement. Default `0.95` (the
  conventional 95% LoA). The mean difference is bracketed by \\\pm
  z\_{(1 + \text{coverage})/2}\\ standard deviations of the differences.

- conf_level:

  Confidence level for the CIs on the LoA themselves. Default `0.95`.

- method:

  How the CIs on the LoA are constructed, following Carkeet (2015):
  `"pair"` (the default) treats the two limits as a pair, so that the
  confidence statement holds for both limits jointly; `"individual"`
  treats each limit separately. See Details for when each is
  appropriate.

## Value

A `data.frame` with rows for the mean difference, the SD of differences,
the lower and upper LoA (`loa_lower`, `loa_upper`), and the lower /
upper CI bounds on each LoA. The rows are the same under both methods;
the construction that produced the CI bounds is recorded in the `method`
attribute.

## Details

**Definition.** For paired observations \\(x_i, y_i)\\, the Bland-Altman
limits of agreement are \$\$\mathrm{LoA}\_\pm \\=\\ \bar d \pm k \cdot
s_d,\$\$ where \\d_i = y_i - x_i\\, \\\bar d\\ is the mean of the
differences, \\s_d\\ is their SD, and \\k = z\_{(1 +
\mathrm{coverage})/2}\\ (for 95% coverage, \\k = 1.96\\). The LoA are
population intervals: they describe the range within which approximately
*coverage*% of *individual* differences are expected to lie if the
differences are normally distributed.

**CIs on the LoA themselves.** The sample LoA are random variables, and
Carkeet (2015) derived exact CIs for them in two forms, selected by
`method`. The choice turns on what the agreement claim is about.

**The pair method (the default).** A Bland-Altman analysis is usually
read as a statement about the range of agreement as a whole, the span
from the lower to the upper LoA within which about *coverage*% of
individual differences lie. That claim involves both limits at once, so
the confidence statement should hold for the two limits jointly; this is
the treatment Carkeet (2015) recommends for most situations. Writing
\\k_t(F)\\ for the exact two-sided normal tolerance factor with
confidence \\F\\ and content equal to `coverage` (Odeh, 1978), the CI on
the upper LoA is \$\$\left\[\\ \bar d + k_t(\alpha/2)\\ s_d, \\\\ \bar
d + k_t(1 - \alpha/2)\\ s_d \\\right\],\$\$ with \\\alpha\\ equal to one
minus `conf_level`, and the CI on the lower LoA is its mirror image
about \\\bar d\\. The joint confidence statement runs through the
probability content of the two symmetric intervals: with confidence
`conf_level`, the interval between the inner pair of bounds, \\\bar d
\pm k_t(\alpha/2)\\ s_d\\, captures less than *coverage*% of the
population of differences, while the interval between the outer pair,
\\\bar d \pm k_t(1 - \alpha/2)\\ s_d\\, captures more, so the pair of
population limits is bracketed simultaneously. In the Bland and Altman
(1986) example that Carkeet reanalyzes (\\n = 17\\, \\\bar d = -2.1\\,
\\s_d = 38.8\\), the pair bounds are \\-2.1 \pm 57.81\\ (inner) and
\\-2.1 \pm 119.60\\ (outer).

**The individual method.** When a single limit carries the substantive
question (for example, only the upper limit matters because only
differences in one direction are clinically consequential), each limit
can be treated on its own. Writing \\t\_{p,\\ n - 1}(\delta)\\ for the
\\p\\ quantile of the noncentral *t* distribution with \\n - 1\\ degrees
of freedom and noncentrality parameter \\\delta = k \sqrt{n}\\, the
exact CI on the upper LoA is \$\$\left\[\\ \bar d +
\frac{s_d}{\sqrt{n}}\\ t\_{\alpha/2,\\ n - 1}(\delta), \\\\ \bar d +
\frac{s_d}{\sqrt{n}}\\ t\_{1 - \alpha/2,\\ n - 1}(\delta)
\\\right\],\$\$ and the CI on the lower LoA uses \\-\delta\\ in place of
\\\delta\\. These intervals are asymmetric about the sample LoA, wider
on the side away from the mean difference. In the worked example above,
the individual CI on the upper LoA is \\\[48.9,\\ 120.0\]\\. The
confidence statement is per limit: each limit is covered with
`conf_level` confidence separately, not both at once.

**Numerical accuracy.** The pair tolerance factors are computed by
numerical integration of the Odeh (1978) chi square by normal integral,
which reproduces Carkeet's Table 2 to all four printed decimals. The
individual quantiles come from
[`stats::qt`](https://rdrr.io/r/stats/TDist.html) with a noncentrality
parameter, so at very large *n* their accuracy is bounded by R's
noncentral *t* algorithm: near \\n = 1000\\ the tolerance coefficient
carries an error of about \\3 \times 10^{-4}\\, far past the sample
sizes at which the exact-versus-approximate distinction matters.

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

Odeh, R. E. (1978). Tables of two-sided tolerance factors for a normal
distribution. *Communications in Statistics - Simulation and
Computation, 7*(2), 183–201.

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
# 1. Two methods that agree well; the CIs treat the limits as a
#    pair (the default):
set.seed(113)
method_a <- rnorm(40, mean = 100, sd = 15)
method_b <- method_a + rnorm(40, mean = 0, sd = 3)
limits_of_agreement(method_a, method_b)
#>  term                  value
#>  mean_difference       0.259
#>  sd_difference         3.37 
#>  loa_lower             -6.34
#>  loa_lower_lower_limit -8.33
#>  loa_lower_upper_limit -5.2 
#>  loa_upper             6.86 
#>  loa_upper_lower_limit 5.72 
#>  loa_upper_upper_limit 8.85 
#> 
#> Confidence level: 95%

# 2. Each limit treated individually, for when a single limit
#    carries the substantive question:
limits_of_agreement(method_a, method_b, method = "individual")
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

# 3. 90% LoA with 95% CIs on the limits:
limits_of_agreement(method_a, method_b, coverage = 0.90, conf_level = 0.95)
#>  term                  value
#>  mean_difference       0.259
#>  sd_difference         3.37 
#>  loa_lower             -5.28
#>  loa_lower_lower_limit -6.95
#>  loa_lower_upper_limit -4.32
#>  loa_upper             5.8  
#>  loa_upper_lower_limit 4.84 
#>  loa_upper_upper_limit 7.47 
#> 
#> Confidence level: 95%
```
