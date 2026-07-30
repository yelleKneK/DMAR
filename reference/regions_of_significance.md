# Regions of Significance for a Covariate by Group Interaction

Finds the values of a covariate at which two groups differ significantly
when the within-group regression slopes are *not* equal, that is, when
there is a covariate-by-group interaction (heterogeneity of regression).
With more than two groups the calculation is carried out for every pair
of groups.

## Usage

``` r
regions_of_significance(
  object,
  data = NULL,
  conf_level = 0.95,
  method = c("simultaneous", "pointwise")
)
```

## Arguments

- object:

  Either a fitted `lm` or `aov` of the form `y ~ x * group`, or a model
  formula of that form. When a formula is supplied, `data` must be
  supplied too and the model is fit with
  [`lm`](https://rdrr.io/r/stats/lm.html).

- data:

  A `data.frame` holding the variables in the formula. Used, and
  required, only when `object` is a formula; ignored with a warning when
  `object` is already a fitted model.

- conf_level:

  Confidence level for the boundaries. Default `0.95`.

- method:

  Character string naming the critical value. `"simultaneous"` (the
  default) uses Potthoff's (1964) simultaneous critical value \\\sqrt{2
  F\_{1 - \alpha; 2, \nu}}\\, which holds the error rate over the whole
  covariate continuum at once. `"pointwise"` uses the classic critical
  value \\t\_{1 - \alpha/2; \nu}\\, which holds it at one covariate
  value chosen in advance.

## Value

A `data.frame` (class `dmar_tbl`) with one row per reported quantity per
pair of groups and columns `pair`, `term`, and a numeric `value`. The
terms, for each pair, are

- `lower_bound`, `upper_bound`:

  The boundaries of the region, sorted, on the scale of the covariate.
  `NA` when a boundary does not exist; when there is a single boundary
  it is reported as `lower_bound` and `upper_bound` is `NA`. Read them
  with `region_code`: the boundaries are the covariate values at which
  the group difference sits exactly on the critical value, and it is
  `region_code` that says on which side of them the groups differ.

- `region_code`:

  How to read the boundaries. `1`: two boundaries, the groups differ
  *outside* them. `2`: two boundaries, the groups differ *between* them.
  `3`: no boundary, the groups differ at every covariate value. `4`: no
  boundary, the groups differ at no covariate value. `5` and `6`: one
  boundary, the groups differ above it (`5`) or below it (`6`).

- `n_boundaries`:

  How many real boundaries exist: 0, 1, or 2.

- `lower_bound_in_range`, `upper_bound_in_range`:

  `1` when the boundary falls inside the covariate values actually
  observed in the two groups, `0` when it falls outside them, `NA` when
  the boundary does not exist. A boundary outside the observed range is
  an extrapolation of the fitted lines and should not be interpreted as
  a covariate value at which anything was or could be seen.

- `difference_intercept`, `difference_slope`:

  The intercept \\d_0\\ and slope \\d_1\\ of the group difference as a
  function of the covariate, so that the estimated difference at
  covariate value \\x\\ is \\d_0 + d_1 x\\.

- `critical_value`:

  The critical value used, \\t\_{crit}\\.

- `df_error`:

  Error degrees of freedom of the fitted model.

- `conf_level`:

  The confidence level.

Everything that is not a number is carried on attributes rather than
forced into `value`: `method`, `outcome`, `covariate`, `group`, the
overall observed `covariate_range`, a `pairs` `data.frame` with the
group labels, the per-pair observed covariate range, and a
plain-language `region` string for each pair, and a `geometry`
`data.frame` holding \\d_0\\, \\d_1\\, and the three elements of their
covariance matrix, which is what
[`plot_regions_of_significance`](https://yelleknek.github.io/DMAR/reference/plot_regions_of_significance.md)
draws.

## Details

**Why the procedure exists.** An analysis of covariance that assumes a
common within-group slope reports one adjusted mean difference, and that
single number is a complete summary of the group comparison only if the
slopes really are common. When the covariate interacts with the group
factor the slopes are not common, the two fitted lines converge or
cross, and there is no such thing as “the” treatment effect: the
difference between the groups depends on where along the covariate you
look. Reporting the adjusted mean difference anyway reports the
difference at one covariate value, the covariate grand mean, and says
nothing about the rest of the range. The question worth answering is
instead *where* on the covariate the groups differ, and that is what
this function answers.

**The calculation.** For two groups, write the estimated difference at
covariate value \\x\\ as the line \$\$\hat D(x) = \hat d_0 + \hat d_1
x,\$\$ where \\\hat d_0\\ is the difference in intercepts and \\\hat
d_1\\ the difference in slopes. Because \\\hat D(x)\\ is a linear
combination of the regression coefficients, its sampling variance
follows from their covariance matrix, \$\$\mathrm{Var}\[\hat D(x)\] =
\mathrm{Var}(\hat d_0) + 2 x \\ \mathrm{Cov}(\hat d_0, \hat d_1) + x^2
\mathrm{Var}(\hat d_1).\$\$ The groups differ significantly at \\x\\
exactly when \\\hat D(x)^2 \> t\_{crit}^2 \\ \mathrm{Var}\[\hat
D(x)\]\\. Setting the two sides equal gives a quadratic in \\x\\,
\$\$(\hat d_1^2 - t\_{crit}^2 \mathrm{Var}(\hat d_1)) x^2 + 2(\hat d_0
\hat d_1 - t\_{crit}^2 \mathrm{Cov}(\hat d_0, \hat d_1)) x + (\hat
d_0^2 - t\_{crit}^2 \mathrm{Var}(\hat d_0)) = 0,\$\$ whose real roots
are the boundaries of the region of significance.

**Every geometry is possible, and the leading coefficient decides
which.** The coefficient on \\x^2\\ is positive exactly when the slope
difference itself clears the critical value. When it is positive the
parabola opens upward and the groups differ *outside* the two
boundaries, the familiar picture of two lines that cross somewhere in
the middle of the covariate and separate at both ends. When it is
negative the parabola opens downward and the groups differ *between* the
boundaries, a middle band of covariate values where the two lines are
far enough apart relative to the precision available there. When there
are no real roots the sign never changes, so the groups differ either
everywhere or nowhere. All of these are reported through `region_code`
rather than being treated as failures.

Of those, “everywhere” is a case the code enumerates but the mathematics
rules out whenever the slopes genuinely differ. The estimated difference
\\\hat D(x)\\ is then a nonconstant line in \\x\\, so it crosses zero at
some covariate value, and at that value the squared difference is zero
while the critical bound is positive. The quadratic therefore always has
two real roots, and the significant set is the pair of tails outside
them or the band between them, never the whole line. Note this is a
statement about the covariate axis extended without limit, not about the
observed data: within the range actually observed the groups may well
differ everywhere, which is why the next paragraph matters.

**Boundaries outside the observed data.** A boundary is a root of an
equation, and the equation is happy to place it far outside the
covariate values that were actually observed. Such a boundary is an
extrapolation of two fitted lines into a region where there are no data
to support them, and it should not be read as a covariate value at which
anything can be claimed. The boundary is still reported, since
suppressing it would hide the shape of the result, but it is flagged in
`lower_bound_in_range` and `upper_bound_in_range` and noted in the
`region` string.

**Simultaneous versus pointwise.** The classic critical value is
\\t\_{1 - \alpha/2}\\, which controls the Type I error rate at a
*single* covariate value fixed in advance. That is not what anyone
actually does: the whole point of the procedure is to scan the covariate
continuum and read off where the groups differ, which is a search over
infinitely many tests. Potthoff (1964) gave the simultaneous critical
value \\\sqrt{2 F\_{1 - \alpha; 2, \nu}}\\, which holds the error rate
over the entire covariate range at once and is therefore the default
here, and the value used in Maxwell, Delaney, and Kelley (2027, Chapter
9). The simultaneous critical value is always the larger of the two, so
its region of significance is always the more conservative one. With
more than two groups the simultaneous guarantee applies to each pair
over the covariate range, not to the family of pairs. For a Bonferroni
protection across the pairs, divide the Type I error rate by the number
of pairs before choosing `conf_level`: with three groups (three pairs)
and a familywise rate of .05, pass `conf_level = 1 - .05 / 3`.

**What the model may contain.** The model must contain exactly one
interaction between a numeric covariate and a grouping factor, and the
grouping factor may not appear in any other term. Additional predictors
that do not interact with the group are allowed and drop out of the
group difference, so `y ~ x * group + block` is fine while
`y ~ x * group * block` is not. A grouping variable stored as a number
(0/1, say) is a numeric predictor to R, not a factor, so convert it with
[`factor`](https://rdrr.io/r/base/factor.html) first.

## References

Johnson, P. O., & Neyman, J. (1936). Tests of certain linear hypotheses
and their application to some educational problems. *Statistical
Research Memoirs, 1*, 57–93.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9 and its extension on heterogeneity of
regression.)

Potthoff, R. F. (1964). On the Johnson-Neyman technique and some
extensions thereof. *Psychometrika, 29*(3), 241–256.
[doi:10.1007/BF02289721](https://doi.org/10.1007/BF02289721)

Rogosa, D. (1980). Comparing nonparallel regression lines.
*Psychological Bulletin, 88*(2), 307–321.
[doi:10.1037/0033-2909.88.2.307](https://doi.org/10.1037/0033-2909.88.2.307)

## See also

[`plot_regions_of_significance`](https://yelleknek.github.io/DMAR/reference/plot_regions_of_significance.md)
to see the group difference and its confidence band across the
covariate;
[`ancova`](https://yelleknek.github.io/DMAR/reference/ancova.md) for the
common-slope analysis and its homogeneity-of-regression test;
[`pygmalion`](https://yelleknek.github.io/DMAR/reference/pygmalion.md)
for the data used below.

Other hypothesis tests:
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
[`dunnett_ci()`](https://yelleknek.github.io/DMAR/reference/dunnett_ci.md),
[`factorial_anova()`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`obrien_test()`](https://yelleknek.github.io/DMAR/reference/obrien_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`scheffe_ci()`](https://yelleknek.github.io/DMAR/reference/scheffe_ci.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`tost_r()`](https://yelleknek.github.io/DMAR/reference/tost_r.md),
[`tost_smd()`](https://yelleknek.github.io/DMAR/reference/tost_smd.md),
[`tukey_kramer_ci()`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# ---- The Pygmalion teacher-expectancy data ----
# Post-test IQ, averaged over the two follow-up assessments, on
# pretest IQ, separately by condition. The slopes differ, so the
# expectancy effect depends on where the child started.
data(pygmalion)
pygmalion$iq_post <- (pygmalion$iq_4 + pygmalion$iq_8) / 2
fit <- lm(iq_post ~ iq_pre * treatment, data = pygmalion)

regions_of_significance(fit)
#>  pair              term                 value
#>  Bloomer - Control lower_bound          104  
#>  Bloomer - Control upper_bound          129  
#>  Bloomer - Control region_code          2    
#>  Bloomer - Control n_boundaries         2    
#>  Bloomer - Control lower_bound_in_range 1    
#>  Bloomer - Control upper_bound_in_range 1    
#>  Bloomer - Control difference_intercept -8.64
#>  Bloomer - Control difference_slope     0.122
#>  Bloomer - Control critical_value       2.46 
#>  Bloomer - Control df_error             306  
#>  Bloomer - Control conf_level           0.95 
#> 
#> Confidence level: 95%

# The plain-language reading of each pair is on an attribute.
attr(regions_of_significance(fit), "pairs")$region
#> [1] "The groups differ for 103.9 < iq_pre < 128.7"

# The pointwise (classic) critical value gives a wider region,
# because it does not pay for scanning the whole covariate.
regions_of_significance(fit, method = "pointwise")
#>  pair              term                 value
#>  Bloomer - Control lower_bound          97.2 
#>  Bloomer - Control upper_bound          175  
#>  Bloomer - Control region_code          2    
#>  Bloomer - Control n_boundaries         2    
#>  Bloomer - Control lower_bound_in_range 1    
#>  Bloomer - Control upper_bound_in_range 0    
#>  Bloomer - Control difference_intercept -8.64
#>  Bloomer - Control difference_slope     0.122
#>  Bloomer - Control critical_value       1.97 
#>  Bloomer - Control df_error             306  
#>  Bloomer - Control conf_level           0.95 
#> 
#> Confidence level: 95%

# ---- Formula interface, and more than two groups ----
set.seed(113)
n <- 150
g <- factor(rep(c("control", "low", "high"), each = n / 3))
x <- rnorm(n, 50, 10)
y <- 2 + 0.5 * x + (g == "high") * (0.4 * x - 15) + rnorm(n, 0, 5)
d <- data.frame(y, x, g)
regions_of_significance(y ~ x * g, data = d)
#>  pair           term                 value 
#>  high - control lower_bound          -114  
#>  high - control upper_bound          45    
#>  high - control region_code          1     
#>  high - control n_boundaries         2     
#>  high - control lower_bound_in_range 0     
#>  high - control upper_bound_in_range 1     
#>  high - control difference_intercept -10.6 
#>  high - control difference_slope     0.295 
#>  high - control critical_value       2.47  
#>  high - control df_error             144   
#>  high - control conf_level           0.95  
#>  low - control  lower_bound          <NA>  
#>  low - control  upper_bound          <NA>  
#>  low - control  region_code          4     
#>  low - control  n_boundaries         0     
#>  low - control  lower_bound_in_range <NA>  
#>  low - control  upper_bound_in_range <NA>  
#>  low - control  difference_intercept 5.58  
#>  low - control  difference_slope     -0.141
#>  low - control  critical_value       2.47  
#>  low - control  df_error             144   
#>  low - control  conf_level           0.95  
#>  low - high     lower_bound          18.7  
#>  low - high     upper_bound          43.2  
#>  low - high     region_code          1     
#>  low - high     n_boundaries         2     
#>  low - high     lower_bound_in_range 0     
#>  low - high     upper_bound_in_range 1     
#>  low - high     difference_intercept 16.1  
#>  low - high     difference_slope     -0.435
#>  low - high     critical_value       2.47  
#>  low - high     df_error             144   
#>  low - high     conf_level           0.95  
#> 
#> Confidence level: 95%
```
