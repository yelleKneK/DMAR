# Provides Games–Howell Simultaneous Confidence Intervals for All Pairwise Comparisons Without Assuming Homogeneity of Variance

Provides Games–Howell Simultaneous Confidence Intervals for All Pairwise
Comparisons Without Assuming Homogeneity of Variance

## Usage

``` r
ci_games_howell(x, group = NULL, conf_level = 0.95)
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
([`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md))
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
[`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md)
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
[`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md)
instead: it pools the variances, so it has more error degrees of freedom
and more power. When only treatments are compared to a single control,
use
[`ci_dunnett`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md).

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

[`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md)
for the homogeneity-assuming counterpart,
[`ci_dunnett`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md)
for many-to-one comparisons,
[`ci_scheffe`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md)
for arbitrary contrasts,
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md)
for the critical value it uses, and
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)
for the tidy methods.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The raw drinks_per_week outcome of the drinks_trial data: the
# Standard arm's variance is more than three times that of either CRA
# arm, so the pooled error term Tukey's HSD relies on is questionable.
ci_games_howell(drinks_trial$drinks_per_week, drinks_trial$treatment)
#>  contrast                    mean_difference se   df   q_statistic lower_limit
#>  CRA - Standard              -24.5           14.2 56.4 -1.72       -72.9      
#>  CRA + Disulfiram - Standard -24.7           15.8 53.2 -1.56       -78.6      
#>  CRA + Disulfiram - CRA      -0.216          12.1 35.7 -0.0178     -42.2      
#>  upper_limit p_adjusted
#>  23.9        0.4481    
#>  29.2        0.5158    
#>  41.8        0.9999    
#> 
#> Confidence level: 95%

# A fitted one-way model may be passed instead of the two vectors.
ci_games_howell(aov(drinks_per_week ~ treatment, data = drinks_trial))
#>  contrast                    mean_difference se   df   q_statistic lower_limit
#>  CRA - Standard              -24.5           14.2 56.4 -1.72       -72.9      
#>  CRA + Disulfiram - Standard -24.7           15.8 53.2 -1.56       -78.6      
#>  CRA + Disulfiram - CRA      -0.216          12.1 35.7 -0.0178     -42.2      
#>  upper_limit p_adjusted
#>  23.9        0.4481    
#>  29.2        0.5158    
#>  41.8        0.9999    
#> 
#> Confidence level: 95%

# With two groups the procedure is Welch's t test, so the limits agree;
# the trial's two enrollment cohorts give a two-group comparison.
ci_games_howell(drinks_trial$drinks_per_week, drinks_trial$cohort)$lower_limit
#> [1] -48.36926
-t.test(drinks_per_week ~ cohort, data = drinks_trial)$conf.int[2]
#> [1] -48.36926
```
