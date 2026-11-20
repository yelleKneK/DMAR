# Confidence Interval for a Single Proportion

The Wilson (1927) score interval for a binomial proportion, the
package's default for proportion inference. Unlike the textbook Wald
interval, it stays inside \[0, 1\] and has positive width when the
number of successes is 0 or *N*. For a 95% interval, its coverage
averaged over a uniformly distributed true proportion is .954 at *N* =
10 and .953 at *N* = 20, against .769 and .846 for the Wald interval;
for some true proportions near 0 or 1 the Wilson coverage falls to about
.84. For other intervals for a binomial proportion, see Brown, Cai, and
DasGupta (2001). The Wald interval is available for instruction and
comparison.

## Usage

``` r
ci_proportion(successes, N, conf_level = 0.95, method = c("wilson", "wald"))
```

## Arguments

- successes:

  Number of successes, a single non-negative integer.

- N:

  Number of trials, a single positive integer at least `successes`.

- conf_level:

  Confidence level. Defaults to 0.95.

- method:

  `"wilson"` (default) or `"wald"`.

## Value

A `data.frame` (class `dmar_tbl`) with rows `ci_lower`, `proportion`,
`ci_upper`, `successes`, and `N`, so the point estimate sits between its
confidence limits.

## References

Brown, L. D., Cai, T. T., & DasGupta, A. (2001). Interval estimation for
a binomial proportion. *Statistical Science, 16*(2), 101–133.
[doi:10.1214/ss/1009213286](https://doi.org/10.1214/ss/1009213286)

Wilson, E. B. (1927). Probable inference, the law of succession, and
statistical inference. *Journal of the American Statistical Association,
22*(158), 209–212.

## See also

[`responder_analysis`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
which uses this interval for each group's responder proportion.

Other confidence intervals:
[`ci_crossed()`](https://yelleknek.github.io/DMAR/reference/ci_crossed.md),
[`ci_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/ci_proportion_paired.md),
[`ci_two_phase()`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
ci_proportion(successes = 17, N = 50)
#>  term       value
#>  ci_lower   0.224
#>  proportion 0.340
#>  ci_upper   0.478
#>  successes  17   
#>  N          50   
#> 
#> Confidence level: 95%

# The Wilson interval stays inside [0, 1] even at the boundary.
ci_proportion(successes = 0, N = 20)
#>  term       value
#>  ci_lower   0.000
#>  proportion 0.000
#>  ci_upper   0.161
#>  successes  0    
#>  N          20   
#> 
#> Confidence level: 95%
```
