# Adjusted Cell and Marginal Means From a Fitted Linear Model

Given a fitted [`lm`](https://rdrr.io/r/stats/lm.html) or
[`aov`](https://rdrr.io/r/stats/aov.html) object with one or more
factors among its predictors, `adjusted_means()` returns the means the
model actually compares, sometimes called least-squares means or
estimated marginal means. By default the table has one row per cell of
the crossed factor design, each cell's mean being the model's predicted
response at that combination of factor levels with every covariate held
at its sample mean (the adjusted cell means of an ANCOVA; for a model
without covariates, the model-based cell means). Naming one or more
factors in `by` instead returns the marginal means of those factors,
formed by averaging the cell predictions over the remaining factors with
either equal or frequency-proportional weights. Every mean is
accompanied by its standard error and a *t* confidence interval on the
model's residual degrees of freedom.

## Usage

``` r
adjusted_means(
  model,
  by = NULL,
  weights = c("equal", "proportional"),
  conf_level = 0.95
)
```

## Arguments

- model:

  A fitted [`lm`](https://rdrr.io/r/stats/lm.html) or
  [`aov`](https://rdrr.io/r/stats/aov.html) object with one or more
  factors (and optionally covariates) on the right-hand side of the
  formula.

- by:

  `NULL` (default) for the cell means table, or a character vector
  naming one or more of the model's factors for their marginal means.
  The output rows cross the named factors in the order given, first
  factor varying fastest.

- weights:

  Weighting used to average cell predictions into marginal means, so it
  matters only when `by` is supplied and the data are unbalanced.
  `"equal"` (default) weights every combination of the averaged-over
  factors equally; `"proportional"` weights each combination by its
  observed frequency.

- conf_level:

  The confidence level for the intervals (default `0.95`).

## Value

A `data.frame` (class `dmar_tbl`) with one row per cell of the reference
grid or, with `by`, one row per combination of the named factors. The
leading columns give the factor levels; the numeric columns are
`estimate` (the adjusted mean), `se` (its standard error), and
`ci_lower` / `ci_upper` (the *t* confidence limits). The residual
degrees of freedom of the intervals are attached as the `df_residual`
attribute and, when `by` is supplied, the weighting as the `weights`
attribute. The stored values keep full precision; only the display
rounds (see
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)).

## Details

**The reference grid and adjusted cell means.** The reference grid is
the crossing of the model's factor levels, enumerated in the order the
factors appear in the model formula with the first factor varying
fastest (the order
[`expand.grid`](https://rdrr.io/r/base/expand.grid.html) produces). This
is the same cell order
[`contrast_adjusted`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md)
expects, so contrast weights can be read off this table row by row.
Every covariate enters the grid at its sample mean, and a transformed
covariate is evaluated by applying the transformation to the mean of the
raw variable: with `log(x)` in the formula the grid carries `mean(x)`
and the model matrix applies [`log()`](https://rdrr.io/r/base/Log.html),
and a `poly(x, 2)` basis is evaluated at \\\bar{x}\\, matching
[`predict`](https://rdrr.io/r/stats/predict.html) on new data at the
covariate mean. Writing \\L\\ for the matrix whose rows are the
design-matrix rows of the grid cells, the cell means are \\L
\hat{\beta}\\, each standard error is the square root of the
corresponding diagonal element of \\L \\ \mathrm{vcov}(\hat{\beta}) \\
L'\\, and each interval is the *t* interval on the model's residual
degrees of freedom.

**Marginal means and the two weightings.** With `by`, the cell
predictions are averaged over the factors not named there, and the
averaging happens in the coefficient map itself: the marginal mean's
\\L\\ row is the weighted average of its cells' rows, so the estimate
and the standard error both follow from one linear function of the
coefficients. `weights = "equal"` weights every combination of the
averaged-over factors equally; this is the population marginal mean of
Searle, Speed, and Milliken (1980), the mean for a population in which
every cell is equally represented regardless of the sample's cell sizes.
`weights = "proportional"` weights each averaged-over combination by its
observed frequency (in a weighted fit, by its total prior weight), so
the marginal mean targets a population whose margins are shaped like the
sample's. With balanced data the two weightings coincide; with
unbalanced data they generally differ, and the choice between them is a
substantive question about the population of interest, not a technical
one (Maxwell, Delaney, and Kelley, 2027, Chapter 7).

**Nonestimable means.** When the fitted design is rank deficient (for
example an empty factorial cell), the model has no predicted value for
the affected cell, and a marginal mean that averages over such a cell
does not exist either. `adjusted_means()` refuses with an error naming
the affected rows rather than reporting a value contaminated by `lm`'s
arbitrary zero for the aliased coefficient.

**Scope.** The function covers single-stratum `lm` and `aov` fits with a
single response. Multi-stratum `aovlist` fits (within-subjects designs
fit with an `Error()` term) are refused, because a within-subjects
marginal mean takes its standard error from the matching error stratum,
which this function does not compute. Factors must enter the model as
variables in the data, not as conversions inside the formula:
`y ~ factor(g) + x` is refused, so convert `g` in the data first.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 7 on nonorthogonal factorial designs and
Chapter 9 on designs with covariates.)

Searle, S. R., Speed, F. M., & Milliken, G. A. (1980). Population
marginal means in the linear model: An alternative to least squares
means. *The American Statistician, 34*(4), 216–221.

## See also

[`contrast_adjusted`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md)
for a confidence interval on a single contrast of the adjusted cell
means; [`ancova`](https://yelleknek.github.io/DMAR/reference/ancova.md)
for the one-way ANCOVA table;
[`ci_dunnett`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md)
for simultaneous many-to-one comparisons.

Other hypothesis tests:
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
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
# 1. Cell means of a 2 x 3 factorial (no covariate): end-of-study IQ in
#    the pygmalion expectancy experiment, grades 1 through 3. The grade
#    factor is created in the data, not inside the formula.
pyg <- subset(pygmalion, grade <= 3)
pyg$grade <- factor(pyg$grade)
fit <- lm(iq_8 ~ treatment * grade, data = pyg)
adjusted_means(fit)
#>  treatment grade estimate se   ci_lower ci_upper
#>  Control   1     102      2.56 97.3     107     
#>  Bloomer   1     114      6.49 101      127     
#>  Control   2     99.8     2.53 94.8     105     
#>  Bloomer   2     118      4.96 109      128     
#>  Control   3     105      2.72 100      111     
#>  Bloomer   3     104      4.76 94.9     114     
#> 
#> Confidence level: 95%

# 2. Marginal means of grade, averaging the cell means over treatment.
adjusted_means(fit, by = "grade")
#>  grade estimate se   ci_lower ci_upper
#>  1     108      3.49 101      115     
#>  2     109      2.78 104      115     
#>  3     105      2.74 99.5     110     
#> 
#> Confidence level: 95%

# 3. An ANCOVA: each adjusted mean holds the covariate, here the pretest
#    depression score, at its sample mean.
fit_ancova <- lm(bdi_post ~ condition + bdi_pre, data = depression_bdi)
adjusted_means(fit_ancova)
#>  condition estimate se   ci_lower ci_upper
#>  ssri      7.54     1.71 4.03     11      
#>  placebo   12       1.71 8.48     15.5    
#>  wait_list 14       1.71 10.5     17.5    
#> 
#> Confidence level: 95%

# 4. With unbalanced cells (few bloomers in every grade) the two
#    weightings answer different questions.
adjusted_means(fit, by = "treatment")
#>  treatment estimate se   ci_lower ci_upper
#>  Control   103      1.5  99.6     106     
#>  Bloomer   112      3.15 106      118     
#> 
#> Confidence level: 95%
adjusted_means(fit, by = "treatment", weights = "proportional")
#>  treatment estimate se   ci_lower ci_upper
#>  Control   102      1.5  99.5     105     
#>  Bloomer   112      3.13 106      119     
#> 
#> Confidence level: 95%
```
