# Equivalence Test for the Standardized Mean Difference via Two One-Sided Tests (TOST)

Performs a two one-sided tests procedure (Schuirmann, 1987) for
equivalence between two independent groups on the standardized mean
difference (Cohen's *d*) scale. The null hypothesis is that the true
\\\delta\\ lies outside the user-specified equivalence bounds
\\\[-\delta_L, \delta_U\]\\; the alternative is that \\\delta\\ lies
inside them. The test is the joint pair of one-sided *t* tests:
\\H\_{0,L}: \delta \le -\delta_L\\ vs \\H\_{1,L}: \delta \> -\delta_L\\,
and \\H\_{0,U}: \delta \ge \delta_U\\ vs \\H\_{1,U}: \delta \<
\delta_U\\. Equivalence is declared when *both* null hypotheses are
rejected at level \\\alpha\\.

## Usage

``` r
equivalence_smd(
  x = NULL,
  y = NULL,
  smd = NULL,
  n_1 = NULL,
  n_2 = NULL,
  smd_lower = NULL,
  smd_upper = NULL,
  alpha_level = 0.05
)
```

## Arguments

- x, y:

  Numeric vectors of observations from the two groups. Alternatively,
  supply `smd`, `n_1`, `n_2` via the summary-statistic interface.

- smd:

  Observed standardized mean difference (Cohen's *d*). Required if `x`
  and `y` are not supplied.

- n_1, n_2:

  Group sample sizes. Required if `x` and `y` are not supplied.

- smd_lower, smd_upper:

  Equivalence bounds on the *d* scale. Both must be positive; the
  equivalence region is \\\[-\delta_L, +\delta_U\]\\, with `smd_lower`
  as \\\delta_L\\ and `smd_upper` as \\\delta_U\\. If only `smd_upper`
  is supplied, the bounds are symmetric: \\\[-\delta_U, +\delta_U\]\\.

- alpha_level:

  One-sided significance level for each of the two tests. Default
  `0.05`. The TOST is then equivalent to a 100(1 - 2\\\alpha\\)% CI on
  \\\delta\\ lying inside the equivalence bounds.

## Value

A `data.frame` with rows for the observed *d*, the two one-sided test
statistics (`t_lower`, `t_upper`) and their degrees of freedom (`df`),
the two one-sided *p*-values (`p_lower`, `p_upper`), the joint TOST
*p*-value (the larger of the two), the 100(1 - 2\\\alpha\\)% CI on
\\\delta\\, the equivalence bounds (`smd_lower`, stored signed as
\\-\delta_L\\, and `smd_upper`), and a binary decision flag
(`equivalent`: 1 = equivalent, 0 = not). The table names the TOST
*p*-value, `p_tost`, in its `p_sentence` attribute, so that
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
pairs the interval with the test of equivalence rather than with
`p_lower`, the first *p*-value row.

## Details

**Schuirmann's TOST.** The TOST procedure tests \\H_0: \delta \le
-\delta_L \cup \delta \ge \delta_U\\ against \\H_1: -\delta_L \< \delta
\< \delta_U\\. Both component tests are rejected (and equivalence is
declared) when the 100(1 - 2\\\alpha\\)% CI on \\\delta\\ lies entirely
inside \\\[-\delta_L, \delta_U\]\\.

**Critical insight.** A non-significant conventional NHST (\\t\\ test of
\\\delta = 0\\) is *not* evidence of equivalence; it only means we
cannot reject \\\delta = 0\\. TOST inverts the testing logic so that "no
meaningful effect" is the alternative, not the null.

**Choosing \\\delta_L\\ and \\\delta_U\\.** The equivalence bounds
(`smd_lower` and `smd_upper`) must be set *before* data collection,
based on what constitutes the smallest effect size of practical interest
(Lakens, Scheel, & Isager, 2018). Common choices in the literature are
0.2 or 0.3, but the bound should reflect domain-specific meaningfulness.

**Connection to the CI.** The TOST rejection at level \\\alpha\\ is
equivalent to the 100(1 - 2\\\alpha\\)% Cohen's-*d* CI (i.e., 90% for
the default \\\alpha = 0.05\\) lying entirely inside the equivalence
region. This is the "two-one-sided" equivalence and matches the Westlake
(1972) rationale for symmetric bioequivalence CIs.

**Standard error and approximation.** Each one-sided component is a Wald
*t* test, \\t = (\hat d \mp \delta) / \mathrm{SE}(\hat d)\\, referred to
a central *t* distribution on \\n_1 + n_2 - 2\\ degrees of freedom, with
the Hedges and Olkin (1985) large-sample standard error
\\\mathrm{SE}(\hat d) = \sqrt{(n_1 + n_2) / (n_1 n_2) + \hat d^{\\2} /
(2 (n_1 + n_2))}\\. This is the asymptotic form of Schuirmann's (1987)
two one-sided tests applied on the standardized scale; it is accurate in
moderate-to-large samples but is an approximation to the exact
noncentral *t* inversion used by
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md), and
the two can differ in small samples. The 100(1 - 2\\\alpha\\)% interval
reported here, \\\hat d \pm t\_{1 - \alpha,\\ n_1 + n_2 - 2}\\
\mathrm{SE}(\hat d)\\, is usually the wider of the two at the same
\\\hat d\\, group sizes, and level. For 90% intervals from equal groups
and \\\|\hat d\|\\ up to 3, it is wider by 18 to 30 percent with 3 per
group and by less than 1 percent with 100 per group. At the 90, 95, and
99 percent levels it was the wider in every case computed (\\\|\hat
d\|\\ up to 12, 2 to 5,000 per group); at the 80 percent level it can be
the narrower, by less than 1 percent, when \\\|\hat d\|\\ is 5 or more.

**Planning.**
[`ss_aipe_equivalence_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md)
plans the sample size for the interval reported here, and
[`ss_aipe_equivalence_smd_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md)
simulates it, so that the interval a study is planned for, the interval
this function reports, and the interval the simulation scores are one
interval, formed by the same code.

## References

Hedges, L. V., & Olkin, I. (1985). *Statistical methods for
meta-analysis*. Academic Press.

Lakens, D. (2017). Equivalence tests: A practical primer for *t* tests,
correlations, and meta-analyses. *Social Psychological and Personality
Science, 8*(4), 355–362.
[doi:10.1177/1948550617697177](https://doi.org/10.1177/1948550617697177)

Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence testing
for psychological research: A tutorial. *Advances in Methods and
Practices in Psychological Science, 1*(2), 259–269.
[doi:10.1177/2515245918770963](https://doi.org/10.1177/2515245918770963)

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*(6), 657–680.

Westlake, W. J. (1972). Use of confidence intervals in analysis of
comparative bioavailability trials. *Journal of Pharmaceutical Sciences,
61*(8), 1340–1341.

## See also

[`equivalence_r`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`ss_aipe_equivalence_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md)

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_dunnett_t3()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett_t3.md),
[`ci_games_howell()`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md),
[`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
[`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
[`equivalence_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
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

Other equivalence testing:
[`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
[`equivalence_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`power_density_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md),
[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Two groups, raw data, symmetric equivalence bounds of 0.4:
set.seed(113)
x <- rnorm(50, 100, 15); y <- rnorm(50, 101, 15)
equivalence_smd(x = x, y = y, smd_upper = 0.4)
#>  term       value 
#>  smd        0.230 
#>  t_lower    3.142 
#>  t_upper    -0.845
#>  df         98    
#>  p_lower    0.0011
#>  p_upper    0.2001
#>  p_tost     0.2001
#>  ci_lower   -0.103
#>  ci_upper   0.564 
#>  smd_lower  -0.400
#>  smd_upper  0.400 
#>  equivalent 0     
#> 
#> Confidence level: 90%

# 2. Summary-statistic interface: published d = 0.05, n = 40 per group.
equivalence_smd(smd = 0.05, n_1 = 40, n_2 = 40, smd_upper = 0.3)
#>  term       value 
#>  smd        0.050 
#>  t_lower    1.565 
#>  t_upper    -1.118
#>  df         78    
#>  p_lower    0.0608
#>  p_upper    0.1335
#>  p_tost     0.1335
#>  ci_lower   -0.322
#>  ci_upper   0.422 
#>  smd_lower  -0.300
#>  smd_upper  0.300 
#>  equivalent 0     
#> 
#> Confidence level: 90%

# 3. Asymmetric bounds (acceptable -0.2 to +0.5):
equivalence_smd(smd = 0.10, n_1 = 80, n_2 = 80,
         smd_lower = 0.2, smd_upper = 0.5)
#>  term       value 
#>  smd        0.100 
#>  t_lower    1.896 
#>  t_upper    -2.528
#>  df         158   
#>  p_lower    0.0299
#>  p_upper    0.0062
#>  p_tost     0.0299
#>  ci_lower   -0.162
#>  ci_upper   0.362 
#>  smd_lower  -0.200
#>  smd_upper  0.500 
#>  equivalent 1     
#> 
#> Confidence level: 90%
```
