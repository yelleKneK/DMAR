# Correlations Adjusted for Common Method Variance, With Intervals and Tests

When every variable in a study is measured the same way, often by one
respondent on one questionnaire, part of each correlation can come from
the shared method rather than from the constructs. The marker variable
technique of Lindell and Whitney (2001) takes a correlation that should
be zero on substantive grounds, the correlation \\r_M\\ of a marker
variable with the substantive variables, as an estimate of the method
effect common to every correlation, and removes it from each substantive
correlation. `common_method_adjusted_r()` returns, for one or more
values of \\r_M\\, every adjusted correlation with its confidence
interval and Lindell and Whitney's *t* test, so that a reader can see
whether a relationship survives the adjustment and how that conclusion
depends on the marker correlation assumed.

The marker correlation is itself an estimate, and the proxy used when no
marker was measured (the smallest positive correlation among the
substantive variables) is a choice. Reading the adjusted correlations
over a grid of values of \\r_M\\, from no adjustment to a generous one,
is the treatment of that uncertainty: a correlation whose interval
excludes zero across the plausible range does not depend on the value
chosen.
[`common_method_marker`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md)
reports the adjusted matrix for a single marker correlation without
intervals or tests, since it does not take \\N\\; this function takes
\\N\\ and gives both.

## Usage

``` r
common_method_adjusted_r(R, N, marker_r = NULL, conf_level = 0.95)

# S3 method for class 'dmar_common_method_adjusted_r'
tidy(x, ...)

# S3 method for class 'dmar_common_method_adjusted_r'
glance(x, ...)
```

## Arguments

- R:

  A correlation matrix among the substantive variables (not the marker):
  symmetric, with a unit diagonal. Its column (or row) names name the
  variables; an unnamed matrix gets `V1`, `V2`, and so on.

- N:

  The sample size on which `R` rests, at least 4.

- marker_r:

  One or more marker correlations \\r_M\\, each in \\(-1, 1)\\, in the
  order the rows should follow. When `NULL` (the default) the smallest
  positive off-diagonal correlation in `R` is used as the proxy marker,
  as in
  [`common_method_marker`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md).

- conf_level:

  Confidence level of the intervals (default 0.95).

- x:

  A `dmar_common_method_adjusted_r` table returned by
  `common_method_adjusted_r()`.

- ...:

  Unused.

## Value

A `data.frame` (class `dmar_tbl`) with one row per value of `marker_r`
and pair of variables, the pairs in the order (1, 2), (1, 3), ..., (2,
3), ... within each marker value, and the columns `marker_r`,
`variable_1`, `variable_2`, `r` (the unadjusted correlation), `estimate`
(the adjusted correlation), `ci_lower`, `ci_upper`, `t_value`, `df`
(\\N - 3\\), and `p_value`. The attributes are `adjusted` (a list of the
adjusted correlation matrices, one per value of `marker_r`, named by
it), `N`, `conf_level`, and `marker_rule` (`"smallest_positive"` when
the proxy was used, `"supplied"` otherwise). `tidy()` returns one row
per pair and marker value, with `term` naming the pair
(`"variable_1:variable_2"`) and the columns `marker_r`, `r`, `estimate`,
`ci_lower`, `ci_upper`, `statistic` (the *t*), `df`, `p_value`, and
`conf_level`; `glance()` returns `n_pairs`, `n_markers`, `N`, and
`conf_level`.

## Details

**The adjustment.** Each substantive correlation \\r\_{ij}\\ is adjusted
as \$\$r^{A}\_{ij} = \frac{r\_{ij} - r_M}{1 - r_M}\$\$ (Lindell &
Whitney, 2001, Equation 4), the partial correlation of \\i\\ and \\j\\
with the method effect removed when that effect contributes the same
\\r_M\\ to every correlation. The adjusted matrix for each value of
`marker_r` is the `"adjusted"` attribute; for a single value it equals
the `"adjusted"` attribute of
[`common_method_marker`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md).

**The test.** Lindell and Whitney (2001, Equation 5) test each adjusted
correlation as a partial correlation with one variable partialed out,
\$\$t = \frac{r^{A}}{\sqrt{(1 - (r^{A})^2) / (N - 3)}},\$\$ referred to
the *t* distribution with \\N - 3\\ degrees of freedom, with a two-sided
*p*-value.

**The interval.** With \\r_M\\ treated as known, the adjustment is an
increasing linear function of the population correlation \\\rho\_{ij}\\,
so a confidence interval for \\\rho\_{ij}\\ maps into one for the
adjusted correlation with the same coverage: each limit of the Fisher's
*Z* interval for \\\rho\_{ij}\\ (the interval of
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
based on the variance \\1/(N - 3)\\) is adjusted by the same equation as
the correlation. A lower limit that the map carries below \\-1\\ is set
to \\-1\\, the bound of a correlation. Because the map is increasing,
that happens exactly when the lower limit for \\\rho\_{ij}\\ falls below
\\2 r_M - 1\\, whatever the sign of \\r\_{ij}\\ or of the adjusted
correlation: never when \\r_M \le 0\\; for \\0 \< r_M \< .5\\ only when
the lower limit for \\\rho\_{ij}\\ is negative, as it is for a negative
correlation, one near zero, or one from a small sample (with \\r\_{ij} =
.35\\, \\N = 5\\, and \\r_M = .3\\ the adjusted correlation is .071 and
its lower limit is set to \\-1\\); and for \\r_M \ge .5\\ even when that
limit is positive. The interval is conditional on \\r_M\\: it carries
the sampling variability of \\r\_{ij}\\ and none of that of the marker
correlation, which is what the grid over `marker_r` is for.

An adjusted correlation at or beyond \\\pm 1\\ lies outside the range of
a correlation. It arises from a correlation below \\2 r_M - 1\\, which a
method effect of \\r_M\\ common to every correlation cannot produce (or
from a correlation of exactly 1), and it gets no interval and no test:
its `ci_lower`, `ci_upper`, `t_value`, and `p_value` are `NA`.

The interval and the test are different approximations. The interval
excludes zero exactly when the Fisher's *Z* interval for \\\rho\_{ij}\\
excludes \\r_M\\, while the test treats the adjusted correlation as a
sample partial correlation, so near the boundary of significance the two
can disagree. The test is kept for parity with the paper; the interval
is what to report beside the estimate. Neither is adjusted for the
number of pairs or of marker values.

**Reading the results.** The proxy marker is one of the substantive
correlations, so the pair that supplied it is adjusted to exactly zero
and its test says nothing about that pair. The technique presumes the
variables are reflected so that their intercorrelations are positive; a
negative correlation is pushed further from zero rather than toward it,
so reverse code as needed before computing `R`. A correlation that
remains different from zero after the adjustment is consistent with a
relationship that is not only a method artifact, under the assumption
that the method effect is the same for every pair; it is not evidence
that no method variance is present.

## References

Lindell, M. K., & Whitney, D. J. (2001). Accounting for common method
variance in cross-sectional research designs. *Journal of Applied
Psychology, 86*(1), 114–121.
[doi:10.1037/0021-9010.86.1.114](https://doi.org/10.1037/0021-9010.86.1.114)

## See also

[`common_method_marker`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md)
for the adjusted matrix and its summary at one marker correlation;
[`common_method_single_factor`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md)
for the single-factor screen;
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
for the interval of an unadjusted correlation.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
[`ci_bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`compare_sem()`](https://yelleknek.github.io/DMAR/reference/compare_sem.md),
[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
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
# Visual, textual, and speed composites of the nine Holzinger and
# Swineford tests in the primary set, each the mean of three tests.
hs <- holzinger_swineford
scores <- data.frame(
  visual  = rowMeans(hs[c("t1_visual_perception", "t2_cubes",
                          "t4_lozenges")]),
  textual = rowMeans(hs[c("t6_paragraph_comprehension", "t7_sentence",
                          "t9_word_meaning")]),
  speed   = rowMeans(hs[c("t10_addition", "t12_counting_groups_of_dots",
                          "t13_straight_and_curved_capitals")]))
R <- stats::cor(scores)

# The smallest positive correlation as the proxy marker; the pair that
# supplied it (textual with speed) is adjusted to zero. Visual with
# speed shows the interval and the test parting near the boundary: the
# interval includes zero while the p-value is below .05.
common_method_adjusted_r(R, N = nrow(scores))
#>  marker_r variable_1 variable_2 r     estimate ci_lower ci_upper t_value df 
#>  0.237    visual     textual    0.312 0.098    -0.041   0.227    1.700   298
#>  0.237    visual     speed      0.326 0.116    -0.021   0.244    2.022   298
#>  0.237    textual    speed      0.237 0.000    -0.144   0.136    0.000   298
#>  p_value
#>  0.0901 
#>  0.0440 
#>  1.0000 
#> 
#> Confidence level: 95%

# A grid of marker correlations, from no adjustment to a large one.
grid <- common_method_adjusted_r(R, N = nrow(scores),
                                 marker_r = c(0, 0.1, 0.2, 0.3))
grid[grid$variable_1 == "visual" & grid$variable_2 == "textual", ]
#>  marker_r variable_1 variable_2 r     estimate ci_lower ci_upper t_value df 
#>  0.000    visual     textual    0.312 0.312    0.206    0.411    5.670   298
#>  0.100    visual     textual    0.312 0.236    0.118    0.345    4.185   298
#>  0.200    visual     textual    0.312 0.140    0.008    0.263    2.442   298
#>  0.300    visual     textual    0.312 0.017    -0.134   0.158    0.298   298
#>  p_value 
#>  < 0.0001
#>  < 0.0001
#>  0.0152  
#>  0.7662  
#> 
#> Confidence level: 95%
attr(grid, "adjusted")[["0.2"]]
#>            visual   textual     speed
#> visual  1.0000000 0.1400825 0.1575495
#> textual 0.1400825 1.0000000 0.0466276
#> speed   0.1575495 0.0466276 1.0000000
```
