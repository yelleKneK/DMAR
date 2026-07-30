# Convert Between the Standardized Mean Difference and the Odds Ratio

Invertible conversions between a two-group standardized mean difference
(Cohen's *d*) and an odds ratio, by the logistic-distribution method of
Hasselblad and Hedges (1995): a continuous outcome split at a threshold
under logistic errors implies \$\$d = \log(\mathrm{OR}) \cdot
\frac{\sqrt{3}}{\pi}, \qquad \mathrm{OR} = \exp\\\bigl(d \cdot \pi /
\sqrt{3}\bigr).\$\$ These conversions let binary-outcome studies enter a
synthesis on the standardized mean difference scale, or mean-difference
studies enter one on the odds ratio scale (Borenstein, Hedges, Higgins,
& Rothstein, 2009, Chapter 7).

## Usage

``` r
convert_d_or(d)

convert_or_d(or)
```

## Arguments

- d:

  The standardized mean difference.

- or:

  The odds ratio, a single positive number.

## Value

A `data.frame` (class `dmar_tbl`) with a single row: term `odds_ratio`
(for `convert_d_or`) or `smd` (for `convert_or_d`) and its `value`.

## References

Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R.
(2009). *Introduction to meta-analysis*. Wiley.

Hasselblad, V., & Hedges, L. V. (1995). Meta-analysis of screening and
diagnostic tests. *Psychological Bulletin, 117*(1), 167–178.
[doi:10.1037/0033-2909.117.1.167](https://doi.org/10.1037/0033-2909.117.1.167)

## See also

[`convert_d_r`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md)
/
[`convert_r_d`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md)
for the correlation leg of the same triangle.

Other parameterization conversions:
[`convert_F_chisq()`](https://yelleknek.github.io/DMAR/reference/convert_F_chisq.md),
[`convert_R2`](https://yelleknek.github.io/DMAR/reference/convert_R2.md),
[`convert_Z_r()`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md),
[`convert_cor_cov()`](https://yelleknek.github.io/DMAR/reference/convert_cor_cov.md),
[`convert_d_r()`](https://yelleknek.github.io/DMAR/reference/convert_d_r.md),
[`convert_r_Z()`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md),
[`convert_t_smd`](https://yelleknek.github.io/DMAR/reference/convert_t_smd.md),
[`convert_z_normal()`](https://yelleknek.github.io/DMAR/reference/convert_z_normal.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# d = 0.5 corresponds to an odds ratio of about 2.48.
convert_d_or(d = 0.5)
#>  term       value
#>  odds_ratio 2.48 

# And back, exactly.
convert_or_d(or = convert_d_or(d = 0.5)$value)
#>  term value
#>  smd  0.5  

# The null maps to the null: d = 0 is an odds ratio of 1.
convert_d_or(d = 0)
#>  term       value
#>  odds_ratio 1    
```
