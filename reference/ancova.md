# Analysis of Covariance (ANCOVA)

Fits a one-way analysis of covariance so the covariate-adjusted group
comparison is available from a single call, without assembling the
adjusted means, the omnibus test, and the effect sizes by hand. It
returns the adjusted (covariate-corrected) group means with standard
errors, the covariate-adjusted omnibus *F* for the group effect (Type
III sums of squares), partial \\\eta^2\\ and partial \\\omega^2\\ with
noncentral *F* confidence intervals, and a homogeneity-of-regression
check, returned in one `data.frame`.

## Usage

``` r
ancova(data, outcome, treatment, covariates, conf_level = 0.95)
```

## Arguments

- data:

  A `data.frame` containing the response, the treatment factor, and the
  covariate(s).

- outcome:

  Character name of the response column.

- treatment:

  Character name of the grouping factor column (the groups being
  compared, for example treatment arms); a factor or character column.

- covariates:

  Character vector of one or more covariate column names.

- conf_level:

  Confidence level for the effect size CIs. Default `0.95`.

## Value

A `data.frame` (class `dmar_tbl`) with rows for the omnibus test
(`F_value`, `df_1`, `df_2`, `p_value`), the sum-of-squares type used
(`sum_of_squares_type`; 3 for Type III), the point estimates and
confidence intervals of partial \\\eta^2\\ and partial \\\omega^2\\, the
adjusted group means and their standard errors (one row per level of
`treatment`), and the homogeneity-of-regression *F*-test. The result
carries the `dmar_tbl` class, so it prints to 3 significant figures with
whole numbers (such as the degrees of freedom) shown without a decimal
part and *p*-values to 4 decimal places (a *p*-value below 0.0001 prints
as “\< 0.0001”); the stored values keep full precision. Control the
display with `print(x, digits = )` or globally with
`options(dmar.digits = )` (see
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)).

## Details

**Covariate-adjusted test (Type III sums of squares).** The omnibus *F*
tests the group effect after adjusting for the covariate(s), that is,
the Type III sum of squares for the grouping factor. For a one-way
ANCOVA (one grouping factor, with the covariate slopes held constant)
the Type II and Type III sums of squares for the group effect coincide,
and both equal the sequential sum of squares obtained with the
covariate(s) entered first and the grouping factor last, which is how it
is computed here; the value matches `car::Anova(fit, type = 3)`. The
choice of sum-of-squares type changes the result only in designs with
more than one factor or with interactions among factors (Maxwell,
Delaney, and Kelley, 2027, Chapter 7); for those, the two-way and mixed
analyses report their sum-of-squares type and allow Type I, II, or III.

**Adjusted means.** The adjusted mean for treatment level \\j\\ is the
model-predicted response at \\X = \bar X\\ (the covariate grand mean):
\$\$\hat \mu_j^{\mathrm{adj}} \\=\\ \hat\mu_j - \sum_k \hat\beta_k (\bar
X\_{kj} - \bar X_k),\$\$ where \\\hat\beta_k\\ is the within-cell slope
on covariate \\k\\ and \\\bar X\_{kj}, \bar X_k\\ are the per-cell and
grand means of covariate \\k\\.

**Homogeneity of regression.** The model fit here holds the within-group
covariate slopes \\\beta_k\\ constant across groups. This is a property
of the particular model being fit, not an assumption of analysis of
covariance in general: it is a testable claim. Adding all
group-by-covariate interactions gives an expanded model, and a model
comparison *F*-test of the additive model against the expanded one
([`stats::anova`](https://rdrr.io/r/stats/anova.html)) assesses whether
the slopes differ across groups (Maxwell, Delaney, and Kelley, 2027,
Chapter 9). A large *F* indicates the slopes are not constant, in which
case the single adjusted comparison is not the whole story and the
interaction model should be entertained directly. The check is reported
in the `F_homogeneity_of_regression` rows.

**Effect size CIs.** Partial \\\eta^2\\ and partial \\\omega^2\\ use the
noncentral *F* framework
([`ci_eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)).

## References

Huitema, B. E. (2011). *The analysis of covariance and alternatives*
(2nd ed.). Wiley.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9 on analysis of covariance and Chapter 7
on higher-order designs.)

## See also

[`ci_c_ancova`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ss_aipe_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova.md),
[`omega_squared_partial`](https://yelleknek.github.io/DMAR/reference/omega_squared_partial.md)

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
[`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
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

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Compare the two groups in the Pygmalion data on eighth-grade IQ,
#    adjusting for pre-test IQ.
ancova(outcome = "iq_8", treatment = "treatment",
       covariates = "iq_pre", data = pygmalion)
#>  term                         value   
#>  F_value                      5.38    
#>  df_1                         1       
#>  df_2                         307     
#>  p_value                      0.0210  
#>  sum_of_squares_type          3       
#>  eta_squared_partial          0.0172  
#>  eta_squared_partial_lower    0.000201
#>  eta_squared_partial_upper    0.056   
#>  omega_squared_partial        0.0139  
#>  omega_squared_partial_lower  0.000201
#>  omega_squared_partial_upper  0.056   
#>  adjusted_mean[Control]       107     
#>  adjusted_mean[Bloomer]       111     
#>  se_adjusted_mean[Control]    0.849   
#>  se_adjusted_mean[Bloomer]    1.67    
#>  F_homogeneity_of_regression  3.88    
#>  df_homogeneity_of_regression 1       
#>  p_homogeneity_of_regression  0.0498  
#> 
#> Confidence level: 95%

# 2. The same comparison adjusting for two covariates (pre-test IQ and
#    grade).
ancova(outcome = "iq_8", treatment = "treatment",
       covariates = c("iq_pre", "grade"), data = pygmalion)
#>  term                         value   
#>  F_value                      5.27    
#>  df_1                         1       
#>  df_2                         306     
#>  p_value                      0.0224  
#>  sum_of_squares_type          3       
#>  eta_squared_partial          0.0169  
#>  eta_squared_partial_lower    0.000127
#>  eta_squared_partial_upper    0.0554  
#>  omega_squared_partial        0.0136  
#>  omega_squared_partial_lower  0.000127
#>  omega_squared_partial_upper  0.0554  
#>  adjusted_mean[Control]       107     
#>  adjusted_mean[Bloomer]       111     
#>  se_adjusted_mean[Control]    0.85    
#>  se_adjusted_mean[Bloomer]    1.67    
#>  F_homogeneity_of_regression  4.83    
#>  df_homogeneity_of_regression 2       
#>  p_homogeneity_of_regression  0.0086  
#> 
#> Confidence level: 95%
```
