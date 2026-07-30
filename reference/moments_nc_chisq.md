# Moments of the Noncentral Chi Square Distribution

Returns the mean, variance, standard deviation, skewness, and excess
kurtosis of a noncentral chi square distribution with `df` degrees of
freedom and noncentrality parameter `ncp`. The noncentral chi square is
the distribution of a sum of squared independent normals with nonzero
means (\\\sum (Z_i + \mu_i)^2\\, with \\\lambda = \sum \mu_i^2\\); it is
the building block of the noncentral *F* (whose numerator is a
noncentral chi square) and the reference distribution for likelihood
ratio and Wald statistics under the alternative. Unlike the noncentral
*t* and *F*, every moment exists, so none of the returned values is ever
`NA`.

## Usage

``` r
moments_nc_chisq(df, ncp = 0)
```

## Arguments

- df:

  Degrees of freedom, a single positive number (need not be a whole
  number).

- ncp:

  Noncentrality parameter \\\lambda\\, a single non-negative number.
  Defaults to 0, the central chi square.

## Value

A `data.frame` (class `dmar_tbl`) in `term` / `value` layout with the
`mean`, `variance`, `sd`, `skewness`, and `excess_kurtosis`, followed by
the `df` and `ncp` that produced them.

## Details

The cumulants of the noncentral chi square are \\\kappa_n =
2^{n-1}(n-1)!\\(\nu + n\lambda)\\ for \\n \ge 1\\, from which the
moments follow directly: the mean is \\\kappa_1 = \nu + \lambda\\, the
variance is \\\kappa_2 = 2(\nu + 2\lambda)\\, the skewness is \\\kappa_3
/ \kappa_2^{3/2} = \sqrt{8}\\(\nu + 3\lambda)/(\nu + 2\lambda)^{3/2}\\,
and the excess kurtosis is \\\kappa_4 / \kappa_2^{2} = 12(\nu +
4\lambda)/(\nu + 2\lambda)^{2}\\. At \\\lambda = 0\\ these reduce to the
central chi square values: mean \\\nu\\, variance \\2\nu\\, skewness
\\\sqrt{8/\nu}\\, and excess kurtosis \\12/\nu\\.

## References

Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous
univariate distributions* (Vol. 2, 2nd ed., Chapter 29). Wiley.

## See also

[`moments_ncf`](https://yelleknek.github.io/DMAR/reference/moments_ncf.md)
(whose numerator is a noncentral chi square) and
[`moments_nct`](https://yelleknek.github.io/DMAR/reference/moments_nct.md)
for the other noncentral moments;
[`conf_limits_nc_chisq`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md)
for the noncentral chi square confidence limits;
[`dchisq`](https://rdrr.io/r/stats/Chisquare.html) for the density.

Other noncentral distribution moments:
[`moments_ncf()`](https://yelleknek.github.io/DMAR/reference/moments_ncf.md),
[`moments_nct()`](https://yelleknek.github.io/DMAR/reference/moments_nct.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A noncentral chi square with 5 df and noncentrality 3.
moments_nc_chisq(df = 5, ncp = 3)
#>  term            value
#>  mean            8    
#>  variance        22   
#>  sd              4.69 
#>  skewness        1.09 
#>  excess_kurtosis 1.69 
#>  df              5    
#>  ncp             3    

# ncp = 0 is the central chi square: mean df, variance 2 * df.
moments_nc_chisq(df = 5)
#>  term            value
#>  mean            5    
#>  variance        10   
#>  sd              3.16 
#>  skewness        1.26 
#>  excess_kurtosis 2.4  
#>  df              5    
#>  ncp             0    

# Every moment exists for any positive df, so nothing is ever NA.
anyNA(moments_nc_chisq(df = 1, ncp = 10)$value)
#> [1] FALSE
```
