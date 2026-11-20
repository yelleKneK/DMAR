# Mixed-Design Multivariate ANOVA With All Four Test Statistics

Computes the multivariate analysis of variance for a mixed-design (one
between-subjects factor and one within-subjects factor), returning
Wilks's \\\Lambda\\, Pillai's trace, Hotelling-Lawley trace, and Roy's
largest root, each with the associated *F*-approximation, degrees of
freedom, and *p*-value. The three effects, between-subjects (*A*),
within-subjects (*B*), and the interaction (\\A \times B\\), are tested
separately, following the multivariate approach to the split-plot design
of Maxwell, Delaney, and Kelley (2027, Section 14.6). Wraps
[`Anova`](https://rdrr.io/pkg/car/man/Anova.html) and returns the result
in the tidy DMAR style.

## Usage

``` r
manova_split_plot(data, within, between, ss_type = 3L)
```

## Arguments

- data:

  A `data.frame` in wide format, one row per subject, with the
  within-subjects measurements in separate columns and the
  between-subjects factor as one additional column.

- within:

  Character vector of column names holding the repeated measures values
  (one column per level of the within- subjects factor). Must be in the
  canonical level order.

- between:

  Character name of the between-subjects factor column in `data`.

- ss_type:

  The sum-of-squares type, passed through to
  [`Anova`](https://rdrr.io/pkg/car/man/Anova.html). Accepts the
  integers `1`, `2`, or `3` or the equivalent Roman-numeral strings
  `"I"`, `"II"`, or `"III"`, and defaults to `3L` (Type III). With one
  between-subjects factor the type matters only for the within-subjects
  main effect when the groups differ in size: Type III tests the
  unweighted grand means of the within-subjects contrasts and Type II
  the weighted grand means (see **Details**). The chosen type is
  recorded in the returned table (see **Value**). Type I (`1` or `"I"`)
  is not available for this mixed design, because
  [`Anova`](https://rdrr.io/pkg/car/man/Anova.html) computes only Type
  II and Type III sums of squares for the multivariate repeated measures
  path; requesting it raises an error.

## Value

A `data.frame` (class `dmar_tbl`) with rows for each of the three
effects crossed with each of the four multivariate statistics, plus one
trailing row recording the sum-of-squares type. Columns: `effect`,
`statistic_name`, `statistic_value`, `F_approx`, `df_1`, `df_2`,
`p_value`. The final row has `effect == "sum_of_squares_type"` and
carries the chosen type (`2` or `3`) in its numeric `statistic_value`;
its remaining numeric columns are `NA`, so the `statistic_value` column
stays numeric. The table prints through the
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)
display layer: whole numbers such as the degrees of freedom without a
decimal part, the statistics to 3 decimal places, and *p*-values to 4
decimal places (one below 0.0001 prints as “\< 0.0001”); the stored
values keep full precision. `tidy()` gives one row per row of the table,
its `term` joining `effect` and `statistic_name` and its `estimate` the
`statistic_value`, so the final row is the term
`"sum_of_squares_type:NA"` (its `statistic_name` is `NA`) with the type
as its estimate; `glance()` gives one row, the number of rows
(`n_terms`) with `conf_level` and `B_used` `NA`, since no interval or
bootstrap is formed (see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)).
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the table as a plain `data.frame`.

## Details

**The three tests.** The multivariate approach replaces each subject's
\\b\\ repeated measures with transformed variables: their mean \\M\\,
and \\b - 1\\ contrasts among them, the \\D\\ variables. The
between-subjects effect (*A*) is the one-way ANOVA of \\M\\ across the
\\a\\ groups (Maxwell, Delaney, & Kelley, 2027, Equation 14.40), which
is also the univariate approach's test of *A*. The interaction (\\A
\times B\\) asks whether the groups differ in their means on the \\D\\
variables, a one-way MANOVA of the \\D\\ variables (their Equations
14.62 and 14.63). The within-subjects main effect (*B*) asks whether the
grand means of the \\D\\ variables are all zero, averaging over the
groups (their Equation 14.54). Both within-subjects tests take as their
error matrix the pooled within-group sums of squares and cross-products
of the \\D\\ variables, and no test depends on which \\b - 1\\ contrasts
are formed. A test of *B* within a single group is a simple effect
(their Section 14.7.4), a different question from the main effect.

**Weighted or unweighted grand means.** When the groups differ in size,
the grand mean of a \\D\\ variable can be the unweighted average of its
\\a\\ group means or the average weighted by group size, which is its
mean over all *N* subjects; the two agree when the groups are of equal
size, and the choice affects only the test of *B*. Maxwell, Delaney, and
Kelley (2027, Section 14.6.9, p. 796, and Chapter 14, note 6, p. 854)
generally prefer the unweighted means, which they note SAS and SPSS test
by default, while allowing that the weighted means suit some studies.
The sum-of-squares type selects between them. Type III (the default)
tests the unweighted grand means: the between-subjects factor is coded
with sum-to-zero contrasts, whatever `options("contrasts")` says, so
that the intercept of the regression of each \\D\\ variable on the
groups is the unweighted mean of its group means, and the hypothesis
matrix is \\a \tilde{n} \bar{D}\_U \bar{D}\_U^\top\\, where
\\\bar{D}\_U\\ holds the unweighted grand means and \\\tilde{n}\\ is the
harmonic mean of the group sizes (their notes 10 and 14). Type II tests
the weighted grand means \\\bar{D}\\, with the hypothesis matrix \\N
\bar{D} \bar{D}^\top\\ (their note 11). Under R's default treatment
coding, the intercept would instead be the first group's means, and its
Type III test the simple effect of *B* in the first group (their
Equation 14.66). With one between-subjects factor, the two types give
the same tests of *A* and of the interaction, and with equal group sizes
the same test of *B* as well. Type I, the sequential decomposition, is
not available here: [`Anova`](https://rdrr.io/pkg/car/man/Anova.html)
computes only Type II and Type III for the multivariate repeated
measures path. The type in force is reported in the returned table so
the analysis is self-documenting.

**The four statistics.** For an effect with \\H\\ and \\E\\ hypothesis-
and error-cross-products matrices:

- Wilks's \\\Lambda = \det(E) / \det(E + H)\\

- Pillai's trace \\V = \mathrm{tr}(H (E + H)^{-1})\\

- Hotelling-Lawley trace \\T_0^2 = \mathrm{tr}(H E^{-1})\\

- Roy's largest root \\\theta = \lambda_1(H E^{-1})\\

**When to use which.** For many data sets the four statistics are likely
to reach the same conclusion about statistical significance, but that
agreement is not guaranteed, and without more information about the
population than is usually available none of them can be said to be
best. Wilks's \\\Lambda\\ is historically the most widely used, and
there is some evidence that Pillai's trace is the most robust of the
four (Maxwell, Delaney, & Kelley, 2027, Section 14.6.11). Roy's largest
root uses only the largest eigenvalue of \\H E^{-1}\\, the largest ratio
of hypothesis to error variation over linear combinations of the
transformed variables. The four statistics give the same *F* whenever
\\s = 1\\, where \\s\\ is the smaller of the number of transformed
variables and the effect's degrees of freedom for each of them (Maxwell,
Delaney, & Kelley, 2027, Section 14.6.11): always for *A* and *B*, and
for the interaction when either factor has two levels.

**Dependency.** Requires the car package on CRAN.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 14, Section 14.6, on the multivariate
approach to split-plot designs, and Chapter 7, Sections 7.9 and 7.10, on
weighted and unweighted marginal means.)

Rencher, A. C., & Christensen, W. F. (2012). *Methods of multivariate
analysis* (3rd ed.). Wiley.
[doi:10.1002/9781118391686](https://doi.org/10.1002/9781118391686)

## See also

[`Anova`](https://rdrr.io/pkg/car/man/Anova.html),
[`manova`](https://rdrr.io/r/stats/manova.html),
[`anova_within_two_way`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md)

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
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`factorial_anova()`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md),
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

Other mixed models:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Two groups of ten measured at three times. Both groups start at the
# same place and rise over time, and group B rises twice as fast, so the
# data generating means differ in slope as well as in level.
set.seed(113)
n_per <- 10
d <- data.frame(
  subject = factor(1:(2 * n_per)),
  group   = factor(rep(c("A", "B"), each = n_per)),
  t1      = c(rnorm(n_per, 0,   1), rnorm(n_per, 0,   1)),
  t2      = c(rnorm(n_per, 0.4, 1), rnorm(n_per, 0.8, 1)),
  t3      = c(rnorm(n_per, 0.8, 1), rnorm(n_per, 1.6, 1))
)

# Rows are the between-subjects effect (labeled A), the within-subjects
# main effect (labeled B), and their interaction, each with all four
# multivariate criteria, followed by a row recording the sum-of-squares
# type. The multivariate tests make no sphericity assumption, which is
# what recommends them over the univariate repeated measures F and its
# epsilon corrections.
manova_split_plot(d, within = c("t1", "t2", "t3"), between = "group")
#>  effect              statistic_name   statistic_value F_approx df_1 df_2
#>  A                   Pillai           0.085           1.681    1    18  
#>  A                   Wilks            0.915           1.681    1    18  
#>  A                   Hotelling-Lawley 0.093           1.681    1    18  
#>  A                   Roy              0.093           1.681    1    18  
#>  B                   Pillai           0.552           10.469   2    17  
#>  B                   Wilks            0.448           10.469   2    17  
#>  B                   Hotelling-Lawley 1.232           10.469   2    17  
#>  B                   Roy              1.232           10.469   2    17  
#>  interaction         Pillai           0.082           0.761    2    17  
#>  interaction         Wilks            0.918           0.761    2    17  
#>  interaction         Hotelling-Lawley 0.090           0.761    2    17  
#>  interaction         Roy              0.090           0.761    2    17  
#>  sum_of_squares_type <NA>             3               <NA>     <NA> <NA>
#>  p_value
#>  0.2112 
#>  0.2112 
#>  0.2112 
#>  0.2112 
#>  0.0011 
#>  0.0011 
#>  0.0011 
#>  0.0011 
#>  0.4825 
#>  0.4825 
#>  0.4825 
#>  0.4825 
#>  <NA>   

# With six subjects left in group A, the types differ in the B rows
# alone: Type III (the default) averages the two groups' changes over
# time equally, and Type II weights each group by its size.
d_unequal <- d[-(1:4), ]
manova_split_plot(d_unequal, within = c("t1", "t2", "t3"),
                  between = "group")
#>  effect              statistic_name   statistic_value F_approx df_1 df_2
#>  A                   Pillai           0.278           5.380    1    14  
#>  A                   Wilks            0.722           5.380    1    14  
#>  A                   Hotelling-Lawley 0.384           5.380    1    14  
#>  A                   Roy              0.384           5.380    1    14  
#>  B                   Pillai           0.555           8.119    2    13  
#>  B                   Wilks            0.445           8.119    2    13  
#>  B                   Hotelling-Lawley 1.249           8.119    2    13  
#>  B                   Roy              1.249           8.119    2    13  
#>  interaction         Pillai           0.080           0.567    2    13  
#>  interaction         Wilks            0.920           0.567    2    13  
#>  interaction         Hotelling-Lawley 0.087           0.567    2    13  
#>  interaction         Roy              0.087           0.567    2    13  
#>  sum_of_squares_type <NA>             3               <NA>     <NA> <NA>
#>  p_value
#>  0.0360 
#>  0.0360 
#>  0.0360 
#>  0.0360 
#>  0.0052 
#>  0.0052 
#>  0.0052 
#>  0.0052 
#>  0.5806 
#>  0.5806 
#>  0.5806 
#>  0.5806 
#>  <NA>   
manova_split_plot(d_unequal, within = c("t1", "t2", "t3"),
                  between = "group", ss_type = 2)
#>  effect              statistic_name   statistic_value F_approx df_1 df_2
#>  A                   Pillai           0.278           5.380    1    14  
#>  A                   Wilks            0.722           5.380    1    14  
#>  A                   Hotelling-Lawley 0.384           5.380    1    14  
#>  A                   Roy              0.384           5.380    1    14  
#>  B                   Pillai           0.573           8.717    2    13  
#>  B                   Wilks            0.427           8.717    2    13  
#>  B                   Hotelling-Lawley 1.341           8.717    2    13  
#>  B                   Roy              1.341           8.717    2    13  
#>  interaction         Pillai           0.080           0.567    2    13  
#>  interaction         Wilks            0.920           0.567    2    13  
#>  interaction         Hotelling-Lawley 0.087           0.567    2    13  
#>  interaction         Roy              0.087           0.567    2    13  
#>  sum_of_squares_type <NA>             2               <NA>     <NA> <NA>
#>  p_value
#>  0.0360 
#>  0.0360 
#>  0.0360 
#>  0.0360 
#>  0.0040 
#>  0.0040 
#>  0.0040 
#>  0.0040 
#>  0.5806 
#>  0.5806 
#>  0.5806 
#>  0.5806 
#>  <NA>   
```
