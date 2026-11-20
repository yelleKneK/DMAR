# Equivalence Test for the Population RMSEA (Test of Not Close Fit)

Tests whether the population root mean square error of approximation
(RMSEA) of a structural equation model is below a bound the researcher
chooses in advance, so that a model can be endorsed as fitting within a
stated tolerance rather than merely not rejected. This is the test of
not close fit of MacCallum, Browne, and Sugawara (1996), which Yuan,
Chan, Marcoulides, and Bentler (2016) present as an equivalence test:
the null hypothesis \\H_0: \epsilon \ge \epsilon_0\\ that the population
RMSEA \\\epsilon\\ is at least the bound \\\epsilon_0\\ is tested
against \\H_1: \epsilon \< \epsilon_0\\, and rejecting it supports fit
within the bound. It reverses the test of close fit (Browne & Cudeck,
1993), whose null hypothesis \\\epsilon \le \epsilon_0\\ a well-fitting
model merely fails to reject. The function works from the summary
statistics a paper reports.

## Usage

``` r
equivalence_rmsea(
  chisq,
  df,
  N,
  rmsea_bound,
  n_groups = 1,
  alpha_level = 0.05,
  sample_size_convention = c("N-1", "N")
)
```

## Arguments

- chisq:

  The model's chi square test statistic \\T\\, the likelihood ratio
  statistic of a normal theory maximum likelihood fit.

- df:

  The model's degrees of freedom, a positive whole number.

- N:

  The total sample size, summed over groups: a whole number greater than
  `n_groups`.

- rmsea_bound:

  The bound \\\epsilon_0\\: the largest population RMSEA the researcher
  is prepared to call fit within tolerance, chosen before seeing the
  data. No universal value is assumed.

- n_groups:

  The number of independent groups of a multiple group model, whose
  `chisq`, `df`, and `N` are totals over the groups. Default 1.

- alpha_level:

  The one-sided Type I error rate, in \\(0, 0.5)\\; default 0.05. The
  interval reported has confidence level \\1 - 2\alpha\\, 90 percent at
  the default.

- sample_size_convention:

  The multiplier \\a\\ in \\T = a\hat{F}\\ with which `chisq` was
  formed: `"N-1"` (the default), each group's sample size less one, so
  \\a = N - G\\ for \\G\\ groups, or `"N"`, so \\a = N\\. lavaan's
  default normal likelihood and Mplus form \\T\\ with \\N\\, so their
  statistics take `"N"`.

## Value

A `dmar_tbl` (class `dmar_sem_fit`) with columns `term` and `value` and
the rows `rmsea` (the point estimate \\\sqrt{G \max(T - df, 0)/(a \cdot
df)}\\), `ci_lower` and `ci_upper` (the \\100(1 - 2\alpha)\\ percent
interval for the population RMSEA), `rmsea_bound`, `null_noncentrality`,
`p_value`, `equivalent` (1 when the test rejects \\H_0\\ at
`alpha_level`, 0 otherwise), `chi_square`, `df`, `N`, and `n_groups`.
The `conf_level` attribute is \\1 - 2\alpha\\. The `verdict` attribute
states the inference in a sentence, and the `sem_fit` attribute records
the sample size convention and multiplier, the boundary convention that
applied, and the lines of the printed footer.

## Details

Write \\T\\ for the chi square statistic on \\df\\ degrees of freedom
and \\a\\ for the multiplier with which it was formed. Under a model
whose misfit is small relative to sampling error, \\T\\ is approximately
noncentral chi square with noncentrality \\\lambda = a \cdot df \cdot
\epsilon^2 / G\\ (Steiger & Lind, 1980; Browne & Cudeck, 1993; Steiger,
1998, for \\G\\ groups). At the boundary of the null hypothesis the
noncentrality is \\\lambda_0 = a \cdot df \cdot \epsilon_0^2 / G\\ (row
`null_noncentrality`), and a small \\T\\ is evidence against \\H_0\\, so
the *p*-value is the lower tail, \\P\\\chi^2\_{df}(\lambda_0) \le T\\\\
(MacCallum et al., 1996, Figure 2). The model is declared within
tolerance (row `equivalent` equal to 1) when the *p*-value is below
`alpha_level`. The *p*-value and the interval stay accurate for a
statistic in the millions, where R's
[`pchisq()`](https://rdrr.io/r/stats/Chisquare.html) fails (see
[`ci_nc_chisq`](https://yelleknek.github.io/DMAR/reference/ci_nc_chisq.md)).

The decision is the same as comparing an interval with the bound. The
interval reported is the two-sided \\100(1 - 2\alpha)\\ percent interval
for the population RMSEA, 90 percent at the default, whose upper limit
is the one-sided \\100(1 - \alpha)\\ percent upper bound; the *p*-value
is below \\\alpha\\ exactly when the upper limit is below \\\epsilon_0\\
(MacCallum et al., 1996, Table 1). This is why the RMSEA is
conventionally reported with a 90 percent interval (Browne & Cudeck,
1993): its upper limit settles the test of not close fit at level .05
and its lower limit the test of close fit. The upper limit is the T-size
RMSEA of Yuan et al. (2016, Eq. 4): the data endorse the model at level
\\\alpha\\ for every bound above it, so it is the quantity to report
when no bound was fixed in advance. The point estimate is truncated at
zero; Steiger (2000) shows that the truncation follows from a basic
property of a good estimator and costs little in practice, and the
interval carries the information about precision. Yuan et al. use \\a =
N - 1\\, the default here. The function does not compute their T-size
CFI (their Eqs. 6 to 10, which divide \\\alpha\\ between the model and
baseline statistics) or the adjusted cutoff values for judging a T-size
index, which their Tables 1 and 2 give as regression coefficients in
\\N\\ and \\df\\).

When \\T\\ lies at or below the lower \\\alpha\\ quantile of the central
chi square distribution, no noncentrality is consistent with it at this
level; the upper limit is then 0, as Yuan et al. (2016, Eq. 4) define
the T-size RMSEA and as lavaan reports it, and the model is within
tolerance for every positive bound. That happens with probability
\\\alpha\\ for a model that fits exactly and more often when the normal
theory statistic is deflated, for example by data with light tails, so a
statistic far below its degrees of freedom deserves a look at the data
before it is read as evidence of fit.

Failing to reject \\H_0\\ does not establish poor fit; it says only that
fit within the bound has not been shown at this sample size. The test
rests on the noncentral chi square approximation for a normal theory
statistic. A scaled or robust statistic does not follow that
distribution with this noncentrality, so it must not be supplied here;
for a fitted model, including a robust (MLR, MLM), categorical (WLSMV),
or FIML fit,
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
reports this test as `p_not_close_fit` from lavaan's corrected measures.
For a multiple group model the convention is Steiger's (1998): `chisq`,
`df`, and `N` are totals over the groups, and the RMSEA and its limits
carry the factor \\\sqrt{G}\\.

## References

Browne, M. W., & Cudeck, R. (1993). Alternative ways of assessing model
fit. In K. A. Bollen & J. S. Long (Eds.), *Testing structural equation
models* (pp. 136–162). Sage.

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
analysis and determination of sample size for covariance structure
modeling. *Psychological Methods, 1*(2), 130–149.
[doi:10.1037/1082-989X.1.2.130](https://doi.org/10.1037/1082-989X.1.2.130)

Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA
fit index. *Structural Equation Modeling, 5*(4), 411–419.
[doi:10.1080/10705519809540115](https://doi.org/10.1080/10705519809540115)

Steiger, J. H. (2000). Point estimation, hypothesis testing, and
interval estimation using the RMSEA: Some comments and a reply to Hayduk
and Glaser. *Structural Equation Modeling, 7*(2), 149–162.

Steiger, J. H., & Lind, J. C. (1980). *Statistically-based tests for the
number of common factors*. Paper presented at the annual Spring meeting
of the Psychometric Society, Iowa City, IA.

Yuan, K.-H., Chan, W., Marcoulides, G. A., & Bentler, P. M. (2016).
Assessing structural equation models by equivalence testing with
adjusted fit indexes. *Structural Equation Modeling, 23*(3), 319–330.
[doi:10.1080/10705511.2015.1065414](https://doi.org/10.1080/10705511.2015.1065414)

## See also

[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
for the test of not close fit from a fitted model (robust, categorical,
and FIML fits included) beside the other fit measures,
[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md) for
an interval from a reported RMSEA, and
[`equivalence_r`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md)
and
[`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md)
for equivalence tests of other effect sizes.

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

Other equivalence testing:
[`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
[`equivalence_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`power_density_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md),
[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md)

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
[`ci_bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_adjusted_r()`](https://yelleknek.github.io/DMAR/reference/common_method_adjusted_r.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`compare_sem()`](https://yelleknek.github.io/DMAR/reference/compare_sem.md),
[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
[`fa_bifactor()`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
[`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
[`fa_second_order()`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md),
[`factor_congruence()`](https://yelleknek.github.io/DMAR/reference/factor_congruence.md),
[`factor_scores()`](https://yelleknek.github.io/DMAR/reference/factor_scores.md),
[`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md),
[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`sem_exact_fit_bootstrap()`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md),
[`sem_fit_measures()`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
[`sem_information_criteria()`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md),
[`sem_residual_fit()`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Yuan et al. (2016, Example 1): nine holzinger_swineford
# tests in the Grant-White school (N = 145), the three-factor model with
# one cross-loading, T = 28.098 on 23 df. The 90 percent upper limit,
# their T-size RMSEA, is .083, so fit within .05 is not established.
equivalence_rmsea(chisq = 28.098, df = 23, N = 145, rmsea_bound = 0.05)
#>  term               value 
#>  rmsea              0.039 
#>  ci_lower           0.000 
#>  ci_upper           0.083 
#>  rmsea_bound        0.050 
#>  null_noncentrality 8.280 
#>  p_value            0.3906
#>  equivalent         0     
#>  chi_square         28.098
#>  df                 23    
#>  N                  145   
#>  n_groups           1     
#> 
#> Confidence level: 90%
#> 
#> Not close fit (RMSEA >= 0.05) is not rejected at alpha = 0.05: fit within the bound is not established (upper limit 0.0828), which is not evidence of poor fit.
#> Sample size: N = 145 in 1 group; multiplier a = 144 (sample_size_convention = "N-1")

# A bound of .09 is endorsed at alpha = .05.
equivalence_rmsea(chisq = 28.098, df = 23, N = 145, rmsea_bound = 0.09)
#>  term               value 
#>  rmsea              0.039 
#>  ci_lower           0.000 
#>  ci_upper           0.083 
#>  rmsea_bound        0.090 
#>  null_noncentrality 26.827
#>  p_value            0.0241
#>  equivalent         1     
#>  chi_square         28.098
#>  df                 23    
#>  N                  145   
#>  n_groups           1     
#> 
#> Confidence level: 90%
#> 
#> Not close fit (RMSEA >= 0.09) is rejected at alpha = 0.05: the population RMSEA is below the bound (upper limit 0.0828).
#> Sample size: N = 145 in 1 group; multiplier a = 144 (sample_size_convention = "N-1")

# The statistic lavaan reports for the two-school model of the same
# tests (formed with N, its default), with Steiger's (1998) multiple
# group convention.
equivalence_rmsea(chisq = 115.8513, df = 48, N = 301, rmsea_bound = 0.10,
                  n_groups = 2, sample_size_convention = "N")
#>  term               value  
#>  rmsea              0.097  
#>  ci_lower           0.075  
#>  ci_upper           0.120  
#>  rmsea_bound        0.100  
#>  null_noncentrality 72.240 
#>  p_value            0.4290 
#>  equivalent         0      
#>  chi_square         115.851
#>  df                 48     
#>  N                  301    
#>  n_groups           2      
#> 
#> Confidence level: 90%
#> 
#> Not close fit (RMSEA >= 0.1) is not rejected at alpha = 0.05: fit within the bound is not established (upper limit 0.12), which is not evidence of poor fit.
#> Sample size: N = 301 in 2 groups; multiplier a = 301 (sample_size_convention = "N")
```
