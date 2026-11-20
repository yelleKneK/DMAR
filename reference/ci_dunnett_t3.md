# Provides Dunnett's T3 Simultaneous Confidence Intervals for All Pairwise Comparisons Without Assuming Homogeneity of Variance

Forms simultaneous confidence intervals and multiplicity-adjusted
*p*-values for all \\a(a-1)/2\\ pairwise differences among \\a\\ group
means without assuming that the groups share a variance. Each pair has
its own standard error and Welch–Satterthwaite degrees of freedom, and
its critical value is that of the studentized maximum modulus for the
\\a(a-1)/2\\ comparisons at those degrees of freedom (Dunnett, 1980).
This is the procedure SPSS reports as Dunnett's T3;
[`cv_smm`](https://yelleknek.github.io/DMAR/reference/cv_smm.md) gives
its critical value.

## Usage

``` r
ci_dunnett_t3(x, group = NULL, conf_level = 0.95)
```

## Arguments

- x:

  Either (a) a fitted [`lm`](https://rdrr.io/r/stats/lm.html) or
  [`aov`](https://rdrr.io/r/stats/aov.html) object with a single factor
  predictor, or (b) a numeric vector of the outcome, in which case
  `group` must also be supplied. An infinite outcome value is an error,
  since it leaves its group's mean and variance undefined.

- group:

  When `x` is a vector, a factor (or coercible to factor) giving the
  group membership of each observation. Observations with a missing
  outcome or a missing group are dropped.

- conf_level:

  Family-wise confidence level. Default `0.95`.

## Value

A `data.frame` with one row per pairwise comparison, in the order
[`ci_games_howell`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md)
uses: with the factor levels in their order, each later level is
compared with each earlier one, and the `contrast` label (`"B - A"`)
names the difference taken. The columns are `contrast`,
`mean_difference`, `se` (the separate-variance standard error of the
difference), `df` (the pair's Welch–Satterthwaite degrees of freedom,
which is why it varies from row to row), `t_value` (`mean_difference`
divided by `se`), `ci_lower`, `ci_upper`, and `p_adjusted`. The
`conf_level` attribute records the family-wise level, and the
`ci_adjustment` and `p_adjustment` attributes are `"dunnett_t3"`. The
table prints through the
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)
display layer and works with
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) (see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)).

## Details

**The procedure.** For groups \\g\\ and \\h\\, with means \\\bar Y_g\\
and \\\bar Y_h\\, variances \\s_g^2\\ and \\s_h^2\\, and sizes \\n_g\\
and \\n_h\\, the standard error of the difference uses only those two
groups' variances, \$\$\mathit{SE}\_{gh} = \sqrt{s_g^2/n_g +
s_h^2/n_h},\$\$ and its degrees of freedom are the Welch–Satterthwaite
value \$\$\mathit{df}\_{gh} = \frac{(s_g^2/n_g +
s_h^2/n_h)^2}{s_g^4/\[n_g^2(n_g-1)\] + s_h^4/\[n_h^2(n_h-1)\]},\$\$ as
in Maxwell, Delaney, and Kelley (2027, Equations 5.13 and 5.14) and in
base R's [`t.test`](https://rdrr.io/r/stats/t.test.html) for two groups.
With \\C = a(a-1)/2\\ comparisons in the family and \\\alpha\\ equal to
one minus `conf_level`, the interval for \\\mu_h - \mu_g\\ is \$\$(\bar
Y_h - \bar Y_g) \pm
m\_{\alpha;\\C,\\\mathit{df}\_{gh}}\\\mathit{SE}\_{gh},\$\$ where
\\m\_{\alpha;C,\mathit{df}}\\ is the upper \\\alpha\\ critical value of
the studentized maximum modulus for \\C\\ comparisons on \\\mathit{df}\\
degrees of freedom (Stoline & Ury, 1979), the value
`cv_smm(alpha_level = 1 - conf_level, df, n_comparisons = C)` returns.
The adjusted *p*-value is the probability that a studentized maximum
modulus variate for \\C\\ comparisons on the pair's degrees of freedom
exceeds the observed \\\|t\| = \|\bar Y_h - \bar Y_g\| /
\mathit{SE}\_{gh}\\, so a pair's `p_adjusted` falls below \\\alpha\\
when its interval excludes zero and not otherwise, to numerical
precision: the critical value's tail area matches \\\alpha\\ to a
relative \\10^{-6}\\ or better. Both are computed by numerical
integration rather than simulation, and the degrees of freedom are not
rounded. As with Games–Howell, the family-wise level holds approximately
rather than exactly, since each pair is referred to a distribution
indexed by its own estimated degrees of freedom.

**Dunnett's T3 and Tamhane's T2.** Tamhane (1979) compared procedures
for pairwise comparisons of means with unequal variances. The one known
as Tamhane's T2 uses the same standard errors and degrees of freedom but
takes its critical value from the inequality of Šidák (1967), treating
the \\C\\ *t* statistics as if they were independent: each pair is
tested at the level \\1 - (1 - \alpha)^{1/C}\\, so its critical value is
the \\(1 + (1 - \alpha)^{1/C})/2\\ quantile of *t* on
\\\mathit{df}\_{gh}\\ and its adjusted *p*-value is \\1 - (1 -
p\_{gh})^C\\, with \\p\_{gh}\\ the pair's two-sided Welch *p*-value.
Dunnett (1980) replaced that bound with the studentized maximum modulus
at the same degrees of freedom. The two distribution functions differ
only in where the power is taken. With \\S =
\sqrt{\chi^2\_{\mathit{df}}/\mathit{df}}\\ and \\\Phi\\ the standard
normal distribution function, the probability that all \\C\\ statistics
fall within \\\pm c\\ is \\\\\mathrm{E}\[2\Phi(cS) - 1\]\\^C\\ as T2
computes it and \\\mathrm{E}\\\[2\Phi(cS) - 1\]^C\\\\ under the
studentized maximum modulus, which by Jensen's inequality is never the
smaller of the two. T3's critical value therefore never exceeds T2's at
the same degrees of freedom: its intervals are never wider, and its
adjusted *p*-values never larger. On the reactivated-story scores of
Maxwell, Delaney, and Kelley (2027, Chapter 12, Exercise 25), three
groups of 13, the T2 *p*-values are .045, .001, and .027 and the T3
values .044, .001, and .026, which are the values SPSS reports for
Dunnett's T3. SPSS offers the two as separate procedures, Tamhane's T2
and Dunnett's T3; the T in T3 comes from Tamhane, who developed its
predecessor (Maxwell, Delaney, & Kelley, 2027, Chapter 5). A formula
that raises the single *t* probability \\2F_t(\|t\|) - 1\\, with \\F_t\\
the distribution function of *t* on \\\mathit{df}\_{gh}\\, to the power
\\C\\ computes T2, whatever it is called.

**Dunnett's T3 or Games–Howell.**
[`ci_games_howell`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md)
uses the same standard errors and degrees of freedom but refers
\\\|t\|\sqrt{2}\\ to the studentized range for \\a\\ means. Dunnett
(1980) found that Games–Howell becomes slightly liberal when the samples
are small, and Maxwell, Delaney, and Kelley (2027, Chapter 5) therefore
recommend Dunnett's T3 when the groups have fewer than roughly 50
observations each and Games–Howell for larger samples. When the
variances are in fact homogeneous,
[`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md)
pools them and so has more error degrees of freedom and more power. With
two groups the T3 interval and *p*-value are exactly those of Welch's
*t* test, `t.test(..., var.equal = FALSE)`, because the studentized
maximum modulus for one comparison is the two-sided *t*.

## References

Dunnett, C. W. (1980). Pairwise multiple comparisons in the unequal
variance case. *Journal of the American Statistical Association,
75*(372), 796–800.
[doi:10.2307/2287161](https://doi.org/10.2307/2287161)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 5 on the multiple-comparisons problem,
where the modifications of Tukey's HSD for unequal variances, Dunnett's
T3 and Games–Howell, are developed.)

Šidák, Z. (1967). Rectangular confidence regions for the means of
multivariate normal distributions. *Journal of the American Statistical
Association, 62*(318), 626–633.
[doi:10.1080/01621459.1967.10482935](https://doi.org/10.1080/01621459.1967.10482935)

Stoline, M. R., & Ury, H. K. (1979). Tables of the Studentized maximum
modulus distribution and an application to multiple comparisons among
means. *Technometrics, 21*(1), 87–93.
[doi:10.1080/00401706.1979.10489726](https://doi.org/10.1080/00401706.1979.10489726)

Tamhane, A. C. (1979). A comparison of procedures for multiple
comparisons of means with unequal variances. *Journal of the American
Statistical Association, 74*(366), 471–480.
[doi:10.2307/2286358](https://doi.org/10.2307/2286358)

## See also

[`cv_smm`](https://yelleknek.github.io/DMAR/reference/cv_smm.md) for the
critical value it uses,
[`ci_games_howell`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md)
for the studentized range counterpart,
[`ci_tukey_kramer`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md)
for the homogeneity-assuming counterpart,
[`ci_dunnett`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md)
for many-to-one comparisons, and
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)
for the tidy methods.

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_games_howell()`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md),
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
# The raw drinks_per_week outcome of the drinks_trial data: the arms hold
# 37, 32, and 19 patients, fewer than the roughly 50 per group below
# which Games-Howell runs slightly liberal, and the Standard arm's
# variance is more than three times that of either CRA arm.
ci_dunnett_t3(drinks_trial$drinks_per_week, drinks_trial$treatment)
#>  contrast                    mean_difference se     df     t_value ci_lower
#>  CRA - Standard              -24.474         20.108 56.392 -1.217  -73.898 
#>  CRA + Disulfiram - Standard -24.690         22.363 53.210 -1.104  -79.753 
#>  CRA + Disulfiram - CRA      -0.216          17.178 35.653 -0.013  -43.164 
#>  ci_upper p_adjusted
#>  24.950   0.5366    
#>  30.373   0.6135    
#>  42.732   1.0000    
#> 
#> Confidence level: 95%, simultaneous interval (Dunnett's T3)

# A fitted one-way model may be passed instead of the two vectors.
ci_dunnett_t3(aov(drinks_per_week ~ treatment, data = drinks_trial))
#>  contrast                    mean_difference se     df     t_value ci_lower
#>  CRA - Standard              -24.474         20.108 56.392 -1.217  -73.898 
#>  CRA + Disulfiram - Standard -24.690         22.363 53.210 -1.104  -79.753 
#>  CRA + Disulfiram - CRA      -0.216          17.178 35.653 -0.013  -43.164 
#>  ci_upper p_adjusted
#>  24.950   0.5366    
#>  30.373   0.6135    
#>  42.732   1.0000    
#> 
#> Confidence level: 95%, simultaneous interval (Dunnett's T3)

# Each half-width is the studentized maximum modulus critical value for
# the three comparisons, at that pair's degrees of freedom, times the
# pair's standard error.
t3 <- ci_dunnett_t3(drinks_trial$drinks_per_week, drinks_trial$treatment)
t3$ci_upper[1] - t3$mean_difference[1]
#> [1] 49.42394
cv_smm(alpha_level = .05, df = t3$df[1], n_comparisons = 3)$value * t3$se[1]
#> [1] 49.42394
```
