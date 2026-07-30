# Moments of the Noncentral *t* Distribution

Returns the mean, variance, standard deviation, skewness, and excess
kurtosis of a noncentral *t* distribution with `df` degrees of freedom
and noncentrality parameter `ncp`. These are the closed-form moments
surveyed by Owen (1968); they are the engine behind the bias and
variance of the standardized mean difference (Cohen's *d*), since *d* is
a scaled noncentral *t* variate. A central *t* (`ncp = 0`) is the
special case with mean 0 and the familiar
\\\mathit{df}/(\mathit{df}-2)\\ variance.

## Usage

``` r
moments_nct(df, ncp = 0)
```

## Arguments

- df:

  Degrees of freedom, a single positive number (need not be a whole
  number).

- ncp:

  Noncentrality parameter \\\delta\\, a single number (may be negative,
  which mirrors the distribution about 0). Defaults to 0, the central
  *t*.

## Value

A `data.frame` (class `dmar_tbl`) in `term` / `value` layout with the
`mean`, `variance`, `sd`, `skewness`, and `excess_kurtosis` (any of
which may be `NA` when the degrees of freedom are too small), followed
by the `df` and `ncp` that produced them.

## Details

Writing the noncentral *t* as \\T = (Z + \delta)/\sqrt{W/\nu}\\ with \\Z
\sim N(0, 1)\\ and \\W \sim \chi^2\_\nu\\ independent, the raw moments
are \$\$\mathrm{E}\[T^k\] = \mathrm{E}\[(Z + \delta)^k\]\\
\Bigl(\tfrac{\nu}{2}\Bigr)^{k/2}\\ \frac{\Gamma\\\bigl((\nu -
k)/2\bigr)}{\Gamma(\nu/2)}, \qquad \nu \> k,\$\$ computed here on the
log scale for stability. The mean exists for \\\nu \> 1\\, the variance
for \\\nu \> 2\\, the skewness for \\\nu \> 3\\, and the excess kurtosis
for \\\nu \> 4\\; a moment whose degrees of freedom condition is not met
is returned as `NA`. The mean is
\\\delta\sqrt{\nu/2}\\\Gamma((\nu-1)/2)/\Gamma(\nu/2)\\, the
\\\delta\\-scaled reciprocal of the Hedges (1981) bias-correction factor
that
[`expected_smd`](https://yelleknek.github.io/DMAR/reference/expected_smd.md)
and [`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) use; that
is why the standardized mean difference is upward biased.

## References

Owen, D. B. (1968). A survey of properties and applications of the
noncentral t-distribution. *Technometrics, 10*(3), 445–478.
[doi:10.1080/00401706.1968.10490590](https://doi.org/10.1080/00401706.1968.10490590)

Hedges, L. V. (1981). Distribution theory for Glass's estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128.

## See also

[`moments_ncf`](https://yelleknek.github.io/DMAR/reference/moments_ncf.md)
for the noncentral *F*;
[`expected_smd`](https://yelleknek.github.io/DMAR/reference/expected_smd.md)
and [`var_smd`](https://yelleknek.github.io/DMAR/reference/var_smd.md)
for the same moments specialized to Cohen's *d*;
[`dt`](https://rdrr.io/r/stats/TDist.html) for the density.

Other noncentral distribution moments:
[`moments_nc_chisq()`](https://yelleknek.github.io/DMAR/reference/moments_nc_chisq.md),
[`moments_ncf()`](https://yelleknek.github.io/DMAR/reference/moments_ncf.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A noncentral t with 20 df and noncentrality 2.5.
moments_nct(df = 20, ncp = 2.5)
#>  term            value
#>  mean            2.6  
#>  variance        1.3  
#>  sd              1.14 
#>  skewness        0.393
#>  excess_kurtosis 0.596
#>  df              20   
#>  ncp             2.5  

# ncp = 0 is the central t: mean 0, variance df / (df - 2), no skew.
moments_nct(df = 10)
#>  term            value
#>  mean            0    
#>  variance        1.25 
#>  sd              1.12 
#>  skewness        0    
#>  excess_kurtosis 1    
#>  df              10   
#>  ncp             0    

# The mean is the noncentrality times the Hedges bias factor's reciprocal,
# which is why Cohen's d (a scaled noncentral t) is upward biased.
m <- moments_nct(df = 18, ncp = 1.2)
m$value[m$term == "mean"]
#> [1] 1.253072
```
