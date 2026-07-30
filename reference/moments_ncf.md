# Moments of the Noncentral *F* Distribution

Returns the mean, variance, standard deviation, skewness, and excess
kurtosis of a noncentral *F* distribution with `df_1` numerator and
`df_2` denominator degrees of freedom and noncentrality parameter `ncp`.
The noncentral *F* is the reference distribution of the *F* statistic
when an effect is present, so its moments describe the sampling behavior
of \\R^2\\, eta squared, and the omnibus *F* test under the alternative.
A central *F* (`ncp = 0`) is the special case.

## Usage

``` r
moments_ncf(df_1, df_2, ncp = 0)
```

## Arguments

- df_1:

  Numerator degrees of freedom, a single positive number.

- df_2:

  Denominator degrees of freedom, a single positive number.

- ncp:

  Noncentrality parameter \\\lambda\\, a single non-negative number.
  Defaults to 0, the central *F*.

## Value

A `data.frame` (class `dmar_tbl`) in `term` / `value` layout with the
`mean`, `variance`, `sd`, `skewness`, and `excess_kurtosis` (any of
which may be `NA` when `df_2` is too small), followed by the `df_1`,
`df_2`, and `ncp` that produced them.

## Details

Writing the noncentral *F* as \\F = (X_1/\nu_1)/(X_2/\nu_2)\\ with \\X_1
\sim \chi^2\_{\nu_1}(\lambda)\\ a noncentral chi square and \\X_2 \sim
\chi^2\_{\nu_2}\\ independent, the raw moments are \$\$\mathrm{E}\[F^k\]
= \Bigl(\tfrac{\nu_2}{\nu_1}\Bigr)^k \mathrm{E}\[X_1^k\]\\
\prod\_{i=1}^{k}\frac{1}{\nu_2 - 2i}, \qquad \nu_2 \> 2k,\$\$ where the
noncentral chi square moments \\\mathrm{E}\[X_1^k\]\\ follow from its
cumulants \\\kappa_n = 2^{n-1}(n-1)!\\(\nu_1 + n\lambda)\\. The mean
exists for \\\nu_2 \> 2\\, the variance for \\\nu_2 \> 4\\, the skewness
for \\\nu_2 \> 6\\, and the excess kurtosis for \\\nu_2 \> 8\\; a moment
whose denominator degrees of freedom condition is not met is returned as
`NA`. The mean reduces to the familiar \\\nu_2(\nu_1 +
\lambda)/\[\nu_1(\nu_2 - 2)\]\\, and the variance to
\\2(\nu_2/\nu_1)^2\[(\nu_1 + \lambda)^2 + (\nu_1 + 2\lambda)(\nu_2 -
2)\] / \[(\nu_2 - 2)^2(\nu_2 - 4)\]\\.

## References

Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous
univariate distributions* (Vol. 2, 2nd ed., Chapter 30). Wiley.

## See also

[`moments_nct`](https://yelleknek.github.io/DMAR/reference/moments_nct.md)
for the noncentral *t*;
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
for the noncentral *F* confidence limits used in effect size intervals;
[`df`](https://rdrr.io/r/stats/Fdist.html) for the density.

Other noncentral distribution moments:
[`moments_nc_chisq()`](https://yelleknek.github.io/DMAR/reference/moments_nc_chisq.md),
[`moments_nct()`](https://yelleknek.github.io/DMAR/reference/moments_nct.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A noncentral F with 3 and 40 df and noncentrality 8.
moments_ncf(df_1 = 3, df_2 = 40, ncp = 8)
#>  term            value
#>  mean            3.86 
#>  variance        5.77 
#>  sd              2.4  
#>  skewness        1.34 
#>  excess_kurtosis 3.12 
#>  df_1            3    
#>  df_2            40   
#>  ncp             8    

# ncp = 0 is the central F: mean df_2 / (df_2 - 2).
moments_ncf(df_1 = 3, df_2 = 40)
#>  term            value
#>  mean            1.05 
#>  variance        0.841
#>  sd              0.917
#>  skewness        1.98 
#>  excess_kurtosis 6.62 
#>  df_1            3    
#>  df_2            40   
#>  ncp             0    

# The variance is undefined for four or fewer denominator df.
moments_ncf(df_1 = 2, df_2 = 4, ncp = 5)
#>  term            value
#>  mean            7    
#>  variance        <NA> 
#>  sd              <NA> 
#>  skewness        <NA> 
#>  excess_kurtosis <NA> 
#>  df_1            2    
#>  df_2            4    
#>  ncp             5    
```
