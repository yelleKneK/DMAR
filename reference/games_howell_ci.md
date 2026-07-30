# Provides Games–Howell Simultaneous Confidence Intervals for All Pairwise Comparisons Without Assuming Homogeneity of Variance

Provides Games–Howell Simultaneous Confidence Intervals for All Pairwise
Comparisons Without Assuming Homogeneity of Variance

## Usage

``` r
games_howell_ci(x, group = NULL, conf_level = 0.95)
```

## Arguments

- x:

  Either (a) a fitted [`lm`](https://rdrr.io/r/stats/lm.html) or
  [`aov`](https://rdrr.io/r/stats/aov.html) object with a single factor
  predictor, or (b) a numeric vector of the outcome, in which case
  `group` must also be supplied.

- group:

  When `x` is a vector, a factor (or coercible to factor) giving the
  group membership of each observation.

- conf_level:

  Family-wise confidence level. Default `0.95`.

## Value

A `data.frame` with one row per pairwise comparison and columns
`contrast`, `mean_difference`, `se`, `df`, `q_statistic`, `lower_limit`,
`upper_limit`, and `p_adjusted`. The `df` column is the
Welch–Satterthwaite degrees of freedom for that pair, which is why it
varies from row to row. The table prints through the
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)
display layer and works with
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) (see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)).

## Details

Tukey's HSD and the Kramer modification for unequal *n*
([`tukey_kramer_ci`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md))
both pool the within-group variances into \\\mathit{MS}\_W\\, so both
assume homogeneity of variance. Neither is robust when that assumption
fails. The Games–Howell procedure drops the assumption: it uses a
separate error term for each pair and a Welch–Satterthwaite degrees of
freedom for each pair, then takes its critical value from the
studentized range.

For groups \\g\\ and \\h\\, the standard error of the difference uses
only those two groups' variances, and the degrees of freedom are
\$\$\mathit{df} = \frac{(s_g^2/n_g +
s_h^2/n_h)^2}{s_g^4/\[n_g^2(n_g-1)\] + s_h^4/\[n_h^2(n_h-1)\]},\$\$ the
same Welch–Satterthwaite expression that base R's
[`t.test`](https://rdrr.io/r/stats/t.test.html) uses by default for two
groups. A pair is declared different when the observed *t* exceeds
\\q/\sqrt{2}\\, with \\q\\ the studentized range critical value
([`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md))
at the pair's degrees of freedom, so the interval for the difference of
means is \\(\bar Y_g - \bar Y_h) \pm
q\_{\alpha;a,\mathit{df}}\sqrt{(s_g^2/n_g + s_h^2/n_h)/2}\\. Maxwell,
Delaney, and Kelley (2027, Chapter 5) develop this as one of the two
modifications of Tukey's HSD for heterogeneous variances (their
Equations 5.13 and 5.14).

**When to use it.** Reach for Games–Howell when the group variances are
not interchangeable and the design is between subjects. It is the
heterogeneity-robust counterpart of
[`tukey_kramer_ci`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md)
and, like it, controls the family-wise error rate across all
\\a(a-1)/2\\ pairs. It handles unequal *n* as a matter of course, so it
does not need a separate unequal-*n* variant.

**When something else is better.** Dunnett (1980) found that
Games–Howell becomes slightly liberal (the family-wise error rate runs
somewhat above the nominal level) when the samples are small. Maxwell,
Delaney, and Kelley (2027, Chapter 5) therefore recommend Games–Howell
for larger samples and Dunnett's T3, which takes its critical value from
the studentized maximum modulus
([`cv_smm`](https://yelleknek.github.io/DMAR/reference/cv_smm.md))
rather than the studentized range, when the groups have fewer than
roughly 50 observations each. When the variances are in fact
homogeneous, use
[`tukey_kramer_ci`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md)
instead: it pools the variances, so it has more error degrees of freedom
and more power. When only treatments are compared to a single control,
use
[`dunnett_ci`](https://yelleknek.github.io/DMAR/reference/dunnett_ci.md).

With \\a = 2\\ groups the procedure is exactly Welch's *t* test: the
interval and the *p*-value equal those from
`t.test(..., var.equal = FALSE)`, because \\q\_{\alpha;2,\mathit{df}} =
\sqrt{2}\\t\_{1-\alpha/2,\mathit{df}}\\.

## References

Games, P. A., & Howell, J. F. (1976). Pairwise multiple comparison
procedures with unequal *n*'s and/or variances: A Monte Carlo study.
*Journal of Educational Statistics, 1*(2), 113–125.
[doi:10.2307/1164979](https://doi.org/10.2307/1164979)

Dunnett, C. W. (1980). Pairwise multiple comparisons in the unequal
variance case. *Journal of the American Statistical Association,
75*(372), 796–800.
[doi:10.2307/2287161](https://doi.org/10.2307/2287161)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 5 on the multiple-comparisons problem,
where the modifications of Tukey's HSD for unequal *n* and unequal
variances are developed.)

## See also

[`tukey_kramer_ci`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md)
for the homogeneity-assuming counterpart,
[`dunnett_ci`](https://yelleknek.github.io/DMAR/reference/dunnett_ci.md)
for many-to-one comparisons,
[`scheffe_ci`](https://yelleknek.github.io/DMAR/reference/scheffe_ci.md)
for arbitrary contrasts,
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md)
for the critical value it uses, and
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)
for the tidy methods.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# PlantGrowth: three groups whose variances are not interchangeable, so
# the pooled error term Tukey's HSD relies on is questionable.
games_howell_ci(PlantGrowth$weight, PlantGrowth$group)
#>  contrast    mean_difference se    df   q_statistic lower_limit upper_limit
#>  trt1 - ctrl -0.371          0.22  16.5 -1.68       -1.17       0.43       
#>  trt2 - ctrl 0.494           0.164 16.8 3.02        -0.101      1.09       
#>  trt2 - trt1 0.865           0.203 14.1 4.26        0.114       1.62       
#>  p_adjusted
#>  0.4746    
#>  0.1129    
#>  0.0237    
#> 
#> Confidence level: 95%

# A fitted one-way model may be passed instead of the two vectors.
games_howell_ci(aov(weight ~ group, data = PlantGrowth))
#>  contrast    mean_difference se    df   q_statistic lower_limit upper_limit
#>  trt1 - ctrl -0.371          0.22  16.5 -1.68       -1.17       0.43       
#>  trt2 - ctrl 0.494           0.164 16.8 3.02        -0.101      1.09       
#>  trt2 - trt1 0.865           0.203 14.1 4.26        0.114       1.62       
#>  p_adjusted
#>  0.4746    
#>  0.1129    
#>  0.0237    
#> 
#> Confidence level: 95%

# With two groups the procedure is Welch's t test, so the limits agree.
pg <- subset(PlantGrowth, group != "trt2")
games_howell_ci(pg$weight, droplevels(pg$group))$lower_limit
#> [1] -1.029516
-t.test(weight ~ group, data = pg)$conf.int[2]
#> [1] -1.029516
```
