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
`contrast`, `mean_difference`, `se`, `df`, `q_statistic`, `ci_lower`,
`ci_upper`, and `p_adjusted`. The `df` column is the Welch–Satterthwaite
degrees of freedom for that pair, which is why it varies from row to
row. The table prints through the
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
Equations 5.13 and 5.14). The adjusted *p*-value is the probability that
the studentized range with \\a\\ means and the pair's degrees of freedom
exceeds the observed \\\|q\|\\.

**How the studentized range is computed.** The per-pair degrees of
freedom are often small and are not whole numbers, and base R's
[`qtukey`](https://rdrr.io/r/stats/Tukey.html) and
[`ptukey`](https://rdrr.io/r/stats/Tukey.html) are not accurate at small
degrees of freedom (with 20 groups and 2 df, `qtukey(.99, 20, 2)` is
44.526 where the exact value is 37.943) and return `NaN` below 2 df,
which a pair of small groups with unequal variances can produce. Where
they are not accurate (fewer than 9 df for up to 6 groups, rising by
steps to fewer than 111 with more than 3000 groups; from 10000 to
\\10^8\\ df; confidence levels above .999; with more than 50 groups,
below the .75 point of the distribution (.80 with more than 1000
groups); and wherever `qtukey` stops short of its root), the quantile
and the *p*-value come from the exact distribution, computed by
integrating the probability integral of the range against the chi square
distribution of the error variance; elsewhere they are `qtukey` and
`ptukey`.
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md)
describes the computation, the regions, and how they were measured.

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
for larger samples and Dunnett's T3
([`ci_dunnett_t3`](https://yelleknek.github.io/DMAR/reference/ci_dunnett_t3.md)),
which takes its critical value from the studentized maximum modulus
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
[`ci_dunnett_t3`](https://yelleknek.github.io/DMAR/reference/ci_dunnett_t3.md)
for Dunnett's T3, the studentized maximum modulus counterpart for
smaller samples,
[`ci_dunnett`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md)
for many-to-one comparisons,
[`ci_scheffe`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md)
for arbitrary contrasts,
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md)
for the critical value it uses, and
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)
for the tidy methods.

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_dunnett_t3()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett_t3.md),
[`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
[`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
[`equivalence_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`factorial_anova()`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`obrien_test()`](https://yelleknek.github.io/DMAR/reference/obrien_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The raw drinks_per_week outcome of the drinks_trial data: the
# Standard arm's variance is more than three times that of either CRA
# arm, so the pooled error term Tukey's HSD relies on is questionable.
ci_games_howell(drinks_trial$drinks_per_week, drinks_trial$treatment)
#>  contrast                    mean_difference se     df     q_statistic ci_lower
#>  CRA - Standard              -24.474         14.219 56.392 -1.721      -72.876 
#>  CRA + Disulfiram - Standard -24.690         15.813 53.210 -1.561      -78.608 
#>  CRA + Disulfiram - CRA      -0.216          12.146 35.653 -0.018      -42.221 
#>  ci_upper p_adjusted
#>  23.929   0.4481    
#>  29.228   0.5158    
#>  41.788   0.9999    
#> 
#> Confidence level: 95%, simultaneous interval (Games-Howell)

# A fitted one-way model may be passed instead of the two vectors.
ci_games_howell(aov(drinks_per_week ~ treatment, data = drinks_trial))
#>  contrast                    mean_difference se     df     q_statistic ci_lower
#>  CRA - Standard              -24.474         14.219 56.392 -1.721      -72.876 
#>  CRA + Disulfiram - Standard -24.690         15.813 53.210 -1.561      -78.608 
#>  CRA + Disulfiram - CRA      -0.216          12.146 35.653 -0.018      -42.221 
#>  ci_upper p_adjusted
#>  23.929   0.4481    
#>  29.228   0.5158    
#>  41.788   0.9999    
#> 
#> Confidence level: 95%, simultaneous interval (Games-Howell)

# With two groups the procedure is Welch's t test, so the limits agree;
# the trial's two enrollment cohorts give a two-group comparison.
ci_games_howell(drinks_trial$drinks_per_week, drinks_trial$cohort)$ci_lower
#> [1] -48.36926
-t.test(drinks_per_week ~ cohort, data = drinks_trial)$conf.int[2]
#> [1] -48.36926
```
