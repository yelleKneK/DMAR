# Plan Sample Size to Make the Test of a Change in the Squared Multiple Correlation Coefficient Sufficiently Powerful

Determine the necessary sample size, or the statistical power at a
specified sample size, for the model comparison *F*-test that a set of
predictors adds nothing to a regression model that already contains
other predictors. The restricted model has `p_restricted` predictors and
the full model adds one or more predictors to it, for `p_full` in all.
The null hypothesis is that the two models have the same population
squared multiple correlation coefficient, so that the change in \\R^2\\
is zero. This is the test of an increment in \\R^2\\ in hierarchical
regression (Cohen, Cohen, West, & Aiken, 2003), with the effect size
\\f^2\\ of Cohen (1988).
[`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md)
plans the omnibus test of all predictors at once, and
[`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md)
the test of one coefficient with fixed predictors.

## Usage

``` r
ss_power_R2_change(
  true_R2_full = NULL,
  true_R2_restricted = NULL,
  true_R2_change = NULL,
  true_f2 = NULL,
  p_full,
  p_restricted,
  alpha_level = 0.05,
  desired_power = 0.85,
  specified_N = NULL,
  random_predictors = TRUE
)
```

## Arguments

- true_R2_full:

  Population squared multiple correlation coefficient of the outcome
  with all `p_full` predictors of the full model, at least 0 and below
  1.

- true_R2_restricted:

  Population squared multiple correlation coefficient of the outcome
  with the `p_restricted` predictors of the restricted model. The
  restricted model's predictors are a subset of the full model's, so
  `true_R2_restricted` cannot exceed `true_R2_full`.

- true_R2_change:

  Population change in \\R^2\\ when the added predictors join the
  restricted model, `true_R2_full` minus `true_R2_restricted`. It is the
  squared semipartial correlation of the added predictor when one
  predictor is added, and the squared multiple semipartial correlation
  of the added set when several are. Give it with `true_R2_full` or with
  `true_R2_restricted`.

- true_f2:

  Cohen's (1988) effect size for the increment, \\f^2 =
  (\rho^2\_{full} - \rho^2\_{restricted}) / (1 - \rho^2\_{full})\\, a
  nonnegative number. Give it alone, in place of the \\R^2\\ values.

- p_full:

  Number of predictors in the full model, not counting the intercept. A
  whole number of at least 1.

- p_restricted:

  Number of predictors in the restricted model, not counting the
  intercept. A whole number from 0 to `p_full - 1`. Zero is the model
  with an intercept only, for which the test is the omnibus test of
  [`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md).

- alpha_level:

  Type I error rate.

- desired_power:

  Desired statistical power, a single number strictly between 0 and 1 (0
  and 1 themselves are refused, with a condition of class
  `dmar_desired_power_invalid`); used when a sample size is sought. With
  `specified_N` the power is computed and `desired_power` is not used,
  but a value outside (0, 1) is still refused.

- specified_N:

  Total sample size at which to report power, in place of a search for
  the necessary sample size. A whole number of at least `p_full + 2`,
  the smallest size that leaves an error degree of freedom.

- random_predictors:

  Whether the predictors are treated as random (`TRUE`, the default) or
  fixed (`FALSE`). See Details.

## Value

A `data.frame` with columns `term` and `value` and nine rows, in this
order:

- `necessary_N`:

  The smallest total sample size at which power is no less than
  `desired_power`. When `specified_N` was given the row is named
  `specified_N` and holds that value.

- `actual_power`:

  Power at that sample size.

- `noncentral_f_parm`:

  The noncentrality parameter \\N f^2\\ under fixed predictors; `NA`
  under random predictors, where the power does not come from a single
  noncentral *F* with that parameter.

- `true_f2`:

  The planning value of Cohen's \\f^2\\ for the increment: `true_f2` as
  given, or the value the \\R^2\\ planning values imply.

- `true_partial_R2`:

  The population squared multiple partial correlation the planning
  values imply, \\f^2 / (1 + f^2)\\.

- `true_R2_change`:

  The population change in \\R^2\\; `NA` when `true_f2` was given, which
  does not determine it.

- `q`:

  The number of predictors removed from the full model to give the
  restricted model, `p_full` minus `p_restricted`: the predictors being
  tested.

- `df_1`:

  Numerator degrees of freedom, equal to `q`.

- `df_2`:

  Denominator degrees of freedom at the sample size reported, \\N -
  p\_{full} - 1\\.

The result carries the `dmar_ss_power` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) apply, and
the attribute `random_predictors` records the predictor model.

## Details

**The effect.** Write \\\rho^2\_{full}\\ and \\\rho^2\_{restricted}\\
(`true_R2_full` and `true_R2_restricted`) for the two population squared
multiple correlation coefficients, \\q\\ for the number of predictors
removed from the full model to give the restricted model (\\q =
p\_{full} - p\_{restricted}\\, the predictors being tested; row `q` of
the result, and \\u\\ in Cohen, 1988), and \\w\\ for the number in the
restricted model, `p_restricted`. In the partial and semipartial
correlation functions
([`var_partial_r`](https://yelleknek.github.io/DMAR/reference/var_partial_r.md),
[`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md),
and their planners) \\p\\ is `p_full`, the number of variables partialed
out is \\p - q\\, and one predictor is of interest, so \\q = 1\\. The
change \\\rho^2\_{full} - \rho^2\_{restricted}\\ is the squared multiple
correlation of the outcome with the part of the added predictors that
the restricted model's predictors do not explain: the squared multiple
semipartial correlation of the added set (Cohen et al., 2003, Section
5.4.3), and with one added predictor its squared semipartial correlation
(Kelley & Maxwell, 2008, Equation 20, p. 176; Aloe & Becker, 2012,
Equation 1). So the answer to whether the change in \\R^2\\ is the
squared semipartial correlation is yes, for one added predictor, and for
a set it is the multiple version of the same quantity. Dividing the
change by \\1 - \rho^2\_{restricted}\\ gives the squared multiple
partial correlation (the squared partial correlation for one added
predictor). Cohen's effect size divides the change by the full model's
unexplained variance, \$\$f^2 = \frac{\rho^2\_{full} -
\rho^2\_{restricted}}{1 - \rho^2\_{full}},\$\$ which is also the squared
multiple partial correlation divided by one minus itself. Power depends
on the two \\R^2\\ values only through \\f^2\\, so the ways of giving
the effect are equivalent and exactly one is given: any two of the three
\\R^2\\ quantities (`true_R2_full`, `true_R2_restricted`,
`true_R2_change`), or `true_f2` alone. The same change is a larger
effect when the full model explains more: a change of .05 is \\f^2 =
.071\\ when the full model explains .30 and \\f^2 = .25\\ when it
explains .80.

**The test.** The statistic is \$\$F = \frac{(R^2\_{full} -
R^2\_{restricted}) / q}{(1 - R^2\_{full}) / (N - p\_{full} - 1)},\$\$
with \\q\\ and \\N - p\_{full} - 1\\ degrees of freedom (the form of
Cohen et al., 2003, Equation 5.5.1), the general model comparison *F* of
Maxwell, Delaney, and Kelley (2027, Chapter 3) applied to the two
regression models. With one added predictor it is the square of the
*t*-statistic of that predictor's coefficient. In the notation of
Maxwell, Delaney, and Kelley (2027), the restricted model has \\df_R =
N - p\_{restricted} - 1\\ error degrees of freedom and the full model
\\df_F = N - p\_{full} - 1\\; the numerator degrees of freedom are
\\df_R - df_F = q\\ (rows `q` and `df_1`), and the denominator degrees
of freedom are \\df_F\\ (row `df_2`). The function takes the numbers of
predictors, since the error degrees of freedom depend on the sample size
being planned. A categorical predictor with *a* levels counts as *a* - 1
predictors.

**Fixed predictors.** When the predictor values are held constant in
hypothetical replications of the study, *F* follows a noncentral *F*
distribution with \\q\\ and \\N - p\_{full} - 1\\ degrees of freedom and
noncentrality \\\lambda = N f^2\\, where the \\R^2\\ values are those of
the design itself (computed from the design's sums of squares, with
divisor \\N\\). The calculation is exact under normal errors. With one
added predictor this is the square of the noncentrality of Kelley and
Maxwell (2008, Equations 16 and 19, p. 175, and Equation 31, p. 181),
who give \\\lambda_j = f_j \sqrt{N}\\ for the *t*-test of one
coefficient with fixed predictors. That is the noncentrality of
[`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
whose power this function reproduces when one predictor is added: the
noncentral *t* of that function has noncentrality \\\sqrt{N} f\\, whose
square is \\N f^2\\. For any number of added predictors \\N f^2\\ is
also the noncentrality of
[`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md)
with fixed predictors.

**Random predictors.** When the predictors are sampled along with the
outcome from a multivariate normal population, the sample squared
multiple partial correlation of the outcome with the \\q\\ added
predictors, given the \\w\\ others, has the distribution that a sample
squared multiple correlation has with \\q\\ predictors and a sample size
of \\N - w\\ (Gatsonis & Sampson, 1989, p. 518). The critical value is
the same as with fixed predictors, since under the null hypothesis the
statistic has the central *F* distribution in both models; only the
power differs. The power is computed from that distribution with the
noncentral *F* approximation of Lee (1971, Section 5.1) that
[`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md)
and [`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md) use.
Against the exact distribution, the power by Lee's approximation
differed by at most .0043 over the combinations examined (squared
partial correlations up to .95, `alpha_level` from .001 to .10, and one
or more error degrees of freedom) and by at most .0004 with 25 or more
error degrees of freedom; the necessary sample size differed from the
one the exact distribution gives, by one case, in 20 of 1,500 searches,
and in those the power by the exact distribution at the size returned
was within .0011 of `desired_power`. At the levels of power used in
planning, power is lower and the necessary sample size larger than with
fixed predictors, and the difference grows with the number of predictors
in the restricted model relative to the sample size. When power is very
low the order can reverse by a small amount (at most .0005 in 1,920
combinations examined, all with no predictors in the restricted model,
`alpha_level` of .001, and power below .09). Predictors are random in
most studies in the behavioral, educational, and social sciences, so
`random_predictors = TRUE` is the default.

**Relation to the test of one regression coefficient.** With one added
predictor and `random_predictors = FALSE` this function and
[`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md)
plan the same test and return the same sample size and power;
`true_R2_Y_X` there is `true_R2_full` here and `true_R2_Y_X_without_j`
is `true_R2_restricted`, with `p` equal to `p_full` and
`p_restricted = p_full - 1`.
[`ss_power_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md)
has fixed predictors only, so with the default of this function (random
predictors) the sample size here is larger: 56 against 50 for a change
from .30 to .40 with five predictors and power .80 (50 is the example of
Kelley & Maxwell, 2008, pp. 181–182). The two functions report the
effect size differently: `true_f2` here is \\f^2\\, and `true_f` there
is \\f\\.

**Cohen's (1988) noncentrality.** Cohen computes power from \\\lambda =
f^2 (u + v + 1)\\, with \\u\\ the number of predictors tested (\\q\\
here) and \\v\\ the error degrees of freedom (as given by Maxwell, 2000,
Equation 1). For the omnibus test \\q + v + 1 = N\\. For an increment
\\v = N - q - w - 1\\, so Cohen's noncentrality is \\f^2 (N - w)\\,
smaller than \\N f^2\\ by \\w f^2\\. Cohen's value equals the
fixed-predictor calculation for \\q\\ predictors, no restricted model,
and `specified_N` set to \\N - w\\, and the last example reproduces it
that way. That \\N - w\\ is only the number that makes the noncentrality
agree; it is not the study's sample size, which remains \\N\\. The
necessary sample size by Cohen's noncentrality exceeds the
fixed-predictor answer of this function by at most \\w\\ cases. Maxwell
(2000, p. 436) notes that \\N f^2\\ is the noncentrality for fixed
predictors and reads Cohen's form as an adjustment toward random
predictors; the random-predictor calculation here makes that adjustment
from the sampling distribution itself. MBESS's `ss.power.R2()` with a
nonzero `Null.R2` computes this quantity when `p` is given as the number
of predictors added: the sample size it reports is Cohen's total sample
size less \\w\\.

**Error from a larger model.** The error variance of the test can be
taken from a model with further predictors that are not part of the
comparison (the Model 2 error of Cohen et al., 2003, Section 5.5.4).
With fixed predictors that calculation is obtained by giving `true_f2`
with the larger model's unexplained variance in its denominator,
`p_full` as the number of predictors in that larger model, and
`p_restricted` as `p_full` minus the number of predictors tested. Set
`random_predictors = FALSE` for this calculation. The random-predictor
calculation is for the comparison of two nested models only.

**Edge cases.** When the change is zero the null hypothesis is true: at
a `specified_N` the power returned is `alpha_level`, and a search stops
with an error of class `dmar_effect_zero`, because no sample size
attains `desired_power`. A `specified_N` below `p_full + 2` leaves no
error degrees of freedom and is an error, as is one that is not a whole
number. With `p_restricted = 0` the restricted model explains nothing,
so `true_R2_restricted`, if given, must be 0 and `true_R2_change`, if
given, must equal `true_R2_full`; the result is then that of
[`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md).
A search starts at `p_full + 2` and returns that size when it already
attains `desired_power`. The search does not go beyond a sample size of
ten million and stops with an error when the power is not attained by
then. R computes the noncentral *F* distribution accurately for a
noncentrality parameter up to about one million and warns beyond it; no
power is returned from an evaluation at which R gave that warning, and
the error says so.

**Estimation.** To plan for a narrow confidence interval rather than for
power when one predictor is added, see
[`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md);
its `true_r_sp` squared is `true_R2_change` and its `p` is `p_full`
(that function counts the predictor of interest and the variables
partialed from it, and one predictor is added, so
`p_restricted + 1 = p_full`). DMAR has no accuracy in parameter
estimation planner for the change in \\R^2\\ from a set of predictors.

## Note

The sample size returned is the smallest whole number at which power is
no less than `desired_power`, so `actual_power` is slightly above the
target. With random predictors that power is the one Lee's approximation
gives; the paragraph on random predictors in Details says how closely
that approximation follows the exact distribution.

## References

Aloe, A. M., & Becker, B. J. (2012). An effect size for regression
predictors in meta-analysis. *Journal of Educational and Behavioral
Statistics, 37*(2), 278–297.
[doi:10.3102/1076998610396901](https://doi.org/10.3102/1076998610396901)

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum. (See Chapter 9.)

Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). *Applied
multiple regression/correlation analysis for the behavioral sciences*
(3rd ed.). Lawrence Erlbaum. (See Sections 5.4 and 5.5.)

Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact power
and sample size calculations. *Psychological Bulletin, 106*(3), 516–524.

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Lee, Y. S. (1971). Some results on the sampling distribution of the
multiple correlation coefficient. *Journal of the Royal Statistical
Society, Series B, 33*(1), 117–130.

Maxwell, S. E. (2000). Sample size and multiple regression analysis.
*Psychological Methods, 5*(4), 434–458.
[doi:10.1037/1082-989X.5.4.434](https://doi.org/10.1037/1082-989X.5.4.434)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on the general form of the model
comparison *F*-test.)

## See also

[`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
[`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md)

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_power_R2()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
[`ss_power_R2_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md),
[`ss_power_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
[`ss_power_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_c_ancova.md),
[`ss_power_composite_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md),
[`ss_power_composite_ancova_2group()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md),
[`ss_power_composite_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_anova.md),
[`ss_power_composite_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova.md),
[`ss_power_composite_factorial_ancova_het()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova_het.md),
[`ss_power_composite_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_anova.md),
[`ss_power_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md),
[`ss_power_contrast()`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`ss_power_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md),
[`ss_power_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
[`ss_power_r()`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
[`ss_power_rc()`](https://yelleknek.github.io/DMAR/reference/ss_power_rc.md),
[`ss_power_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`ss_power_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md),
[`ss_power_rm_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_sc()`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Five predictors are expected to explain .25 of the variance and four
# added predictors to raise that to .30. Necessary sample size for power
# .80 at alpha .01, with random predictors (the default).
ss_power_R2_change(true_R2_restricted = .25, true_R2_full = .30,
                   p_restricted = 5, p_full = 9, alpha_level = .01,
                   desired_power = .80)
#>  term              value
#>  necessary_N       250  
#>  actual_power      0.802
#>  noncentral_f_parm <NA> 
#>  true_f2           0.071
#>  true_partial_R2   0.067
#>  true_R2_change    0.050
#>  q                 4    
#>  df_1              4    
#>  df_2              240  

# The same plan with predictors fixed by design.
ss_power_R2_change(true_R2_restricted = .25, true_R2_full = .30,
                   p_restricted = 5, p_full = 9, alpha_level = .01,
                   desired_power = .80, random_predictors = FALSE)
#>  term              value 
#>  necessary_N       242   
#>  actual_power      0.802 
#>  noncentral_f_parm 17.286
#>  true_f2           0.071 
#>  true_partial_R2   0.067 
#>  true_R2_change    0.050 
#>  q                 4     
#>  df_1              4     
#>  df_2              232   

# The effect given as the change with the full model's value, and as
# Cohen's f squared. Power at a sample size of 90, fixed predictors.
ss_power_R2_change(true_R2_change = .05, true_R2_full = .30,
                   p_restricted = 5, p_full = 9, alpha_level = .01,
                   specified_N = 90, random_predictors = FALSE)
#>  term              value
#>  specified_N       90   
#>  actual_power      0.241
#>  noncentral_f_parm 6.429
#>  true_f2           0.071
#>  true_partial_R2   0.067
#>  true_R2_change    0.050
#>  q                 4    
#>  df_1              4    
#>  df_2              80   
ss_power_R2_change(true_f2 = .05 / .70, p_restricted = 5, p_full = 9,
                   alpha_level = .01, specified_N = 90,
                   random_predictors = FALSE)
#>  term              value
#>  specified_N       90   
#>  actual_power      0.241
#>  noncentral_f_parm 6.429
#>  true_f2           0.071
#>  true_partial_R2   0.067
#>  true_R2_change    <NA> 
#>  q                 4    
#>  df_1              4    
#>  df_2              80   

# The effect given as the restricted model's value and the change, the
# way a hierarchical regression is often described. Random predictors.
ss_power_R2_change(true_R2_restricted = .25, true_R2_change = .05,
                   p_restricted = 5, p_full = 9, alpha_level = .01,
                   desired_power = .80)
#>  term              value
#>  necessary_N       250  
#>  actual_power      0.802
#>  noncentral_f_parm <NA> 
#>  true_f2           0.071
#>  true_partial_R2   0.067
#>  true_R2_change    0.050
#>  q                 4    
#>  df_1              4    
#>  df_2              240  

# One added predictor with fixed predictors is the test of its regression
# coefficient: the two functions return the same sample size.
ss_power_R2_change(true_R2_restricted = .30, true_R2_full = .40,
                   p_restricted = 4, p_full = 5, desired_power = .80,
                   random_predictors = FALSE)
#>  term              value
#>  necessary_N       50   
#>  actual_power      0.806
#>  noncentral_f_parm 8.333
#>  true_f2           0.167
#>  true_partial_R2   0.143
#>  true_R2_change    0.100
#>  q                 1    
#>  df_1              1    
#>  df_2              44   
ss_power_reg_coef(true_R2_Y_X = .40, true_R2_Y_X_without_j = .30, p = 5,
                  desired_power = .80)
#>  term              value
#>  necessary_N       50   
#>  actual_power      0.806
#>  noncentral_t_parm 2.887
#>  true_f            0.408

# Cohen's (1988) noncentrality for the design of the third example,
# f squared times (N - w) with w = 5: four predictors and no restricted
# model, with specified_N set to 90 - 5 so the noncentrality matches
# Cohen's. The study's sample size is still 90.
ss_power_R2_change(true_f2 = .05 / .70, p_restricted = 0, p_full = 4,
                   alpha_level = .01, specified_N = 85,
                   random_predictors = FALSE)
#>  term              value
#>  specified_N       85   
#>  actual_power      0.223
#>  noncentral_f_parm 6.071
#>  true_f2           0.071
#>  true_partial_R2   0.067
#>  true_R2_change    <NA> 
#>  q                 4    
#>  df_1              4    
#>  df_2              80   
```
