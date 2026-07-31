# Confidence Interval for a Single Proportion

The Wilson (1927) score interval for a binomial proportion, the
package's default for proportion inference: unlike the textbook Wald
interval it cannot escape \[0, 1\], behaves sensibly at 0 and 1 counts,
and holds close to nominal coverage at small *n* (Brown, Cai, &
DasGupta, 2001, recommend it for general use). The Wald interval is
available for instruction and comparison.

## Usage

``` r
ci_proportion(successes, n, conf_level = 0.95, method = c("wilson", "wald"))
```

## Arguments

- successes:

  Number of successes, a single non-negative integer.

- n:

  Number of trials, a single positive integer at least `successes`.

- conf_level:

  Confidence level. Defaults to 0.95.

- method:

  `"wilson"` (default) or `"wald"`.

## Value

A `data.frame` (class `dmar_tbl`) with rows `lower_limit`, `proportion`,
`upper_limit`, `successes`, and `n`, so the point estimate sits between
its confidence limits.

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

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
ci_proportion(successes = 17, n = 50)
#>  term        value
#>  lower_limit 0.224
#>  proportion  0.34 
#>  upper_limit 0.478
#>  successes   17   
#>  n           50   
#> 
#> Confidence level: 95%

# The Wilson interval stays inside [0, 1] even at the boundary.
ci_proportion(successes = 0, n = 20)
#>  term        value   
#>  lower_limit 1.39e-17
#>  proportion  0       
#>  upper_limit 0.161   
#>  successes   0       
#>  n           20      
#> 
#> Confidence level: 95%
```
