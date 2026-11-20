# Simple Effect F Tests for a Two-Factor Between-Subjects Design

Given a fitted [`aov`](https://rdrr.io/r/stats/aov.html) or
[`lm`](https://rdrr.io/r/stats/lm.html) object for a two-factor
between-subjects design, conventionally written \\Y \sim A \times B\\,
where \\A\\ and \\B\\ are crossed fixed factors, computes the family of
*simple main effects*: the effect of \\A\\ at each level of \\B\\ and/or
the effect of \\B\\ at each level of \\A\\. Each row carries the simple
effect *F* test, its (optionally adjusted) *p*-value, and the partial
\\\eta^2\\ with a noncentrality-based confidence interval under the
pooled error term. The error term can be either the full-model pooled
\\\mathit{MS}\_W\\ (the textbook default) or a Welch–Satterthwaite test
refitted within each conditioning level (robust to within-level
heteroscedasticity).

## Usage

``` r
simple_effects_AB(
  object,
  which = "both",
  error_term = "pooled",
  adjust = "none",
  conf_level = 0.95
)
```

## Arguments

- object:

  A fitted [`aov`](https://rdrr.io/r/stats/aov.html) or
  [`lm`](https://rdrr.io/r/stats/lm.html) object whose right-hand side
  has *exactly two* crossed factors (e.g.\\
  `iq_gain ~ treatment * grade`). The interaction term is strongly
  recommended so that the pooled error is the pure within-cell
  \\\mathit{MS}\_W\\; the function still runs without it but issues a
  warning (the additive-model residual includes interaction variance and
  inflates the error term used for the pooled simple effect *F*).

- which:

  Which family of simple effects to report:

  `"both"` (default)

  :   The \\b\\ tests of \\A\\ at each level of \\B\\ followed by the
      \\a\\ tests of \\B\\ at each level of \\A\\ (\\a + b\\ rows).

  `"A_at_B"`

  :   Only the \\b\\ tests of the first factor at each level of the
      second.

  `"B_at_A"`

  :   Only the \\a\\ tests of the second factor at each level of the
      first.

  The first factor on the right-hand side of the model formula is
  treated as \\A\\; the second as \\B\\.

- error_term:

  Error-term strategy for the simple effect *F*:

  `"pooled"` (default)

  :   Use the full factorial model's \\\mathit{MS}\_W\\ and its residual
      *df*. This is Maxwell, Delaney, and Kelley's preferred default and
      gives the simple effect *F* maximum denominator *df*. Validity
      rests on homogeneity of variance across all \\a \times b\\ cells.

  `"welch"`

  :   Refit a Welch–Satterthwaite one-way test on only the data at the
      conditioning level (using
      [`oneway.test`](https://rdrr.io/r/stats/oneway.test.html) with
      `var.equal = FALSE`). Robust to heteroscedasticity within the
      conditioning level at the cost of fewer (and fractional)
      denominator *df*.

- adjust:

  Multiple-comparison adjustment applied to the *p*-values of the entire
  family of simple effects returned (\\a + b\\ for `which = "both"`).
  One of `"none"` (default), `"bonferroni"`, or any sequential method
  supported by [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html):
  `"holm"`, `"hochberg"`, `"BH"`, `"BY"`.

- conf_level:

  Confidence level for each row's partial \\\eta^2\\ interval. Default
  `0.95`.

## Value

A `data.frame` with one row per simple effect test and columns `effect`,
`focal_factor`, `conditioning_factor`, `conditioning_level`, `F_value`,
`df_effect`, `df_error`, `p_value`, `p_adjusted`, `eta_squared_partial`,
`ci_lower`, `ci_upper`, `n_at_level`. Attributes `error_term`, `adjust`,
`conf_level`, `factor_A`, and `factor_B` record the call options.

## Details

**What a simple effect is.** The simple main effect of \\A\\ at level
\\B = b_j\\ tests whether the \\a\\ cell means at that single level of
\\B\\ differ. It is the one-way analysis of variance of \\Y\\ on \\A\\
restricted to observations with \\B = b_j\\. The counterpart, the simple
effect of \\B\\ at \\A = a_i\\, is defined symmetrically.

**Test statistic.** For the simple effect of \\A\\ at \\B = b_j\\, let
\\\mathit{SS}\_{A\\\|\\b_j} = \sum_i n\_{ij}\\(\bar{Y}\_{ij\cdot} -
\bar{Y}\_{\cdot j\cdot})^2\\ be the between-A sum of squares computed at
that level, and let \\\mathit{MS}\_{A\\\|\\b_j} =
\mathit{SS}\_{A\\\|\\b_j} / (a - 1)\\.

- *Pooled* \\\mathit{MS}\_W\\: \\F = \mathit{MS}\_{A\\\|\\b_j} /
  \mathit{MS}\_W\\ with degrees of freedom \\(a - 1,\\ N - ab)\\, where
  \\\mathit{MS}\_W\\ and its *df* come from the fitted full factorial
  model.

- *Welch*: \\F\\ and its (fractional) denominator *df* come from
  [`oneway.test`](https://rdrr.io/r/stats/oneway.test.html)`(y ~ A, subset = (B == b_j), var.equal = FALSE)`.

Test statistics for \\B\\ at \\a_i\\ are computed by interchanging the
two factors.

**Choosing an error term.** The pooled denominator borrows strength from
all \\N\\ observations and is the textbook default in Maxwell, Delaney,
and Kelley's treatment. Its validity rests on homogeneity of variance
across *all* \\a \times b\\ cells, not merely within the conditioning
level. When that assumption is doubtful , for example, if Levene's or
the Brown–Forsythe test flags heteroscedasticity, or if cell variances
visibly differ, the Welch option provides a level-conditional test that
does not require homogeneity across cells. The trade-off is denominator
*df*: pooled carries the full \\N - ab\\ residual *df*, whereas Welch
carries the Welch–Satterthwaite *df* based on the \\a\\ cell variances
at that level only.

**Partial \\\eta^2\\ and its CI.** The point estimate is
\$\$\hat{\eta}^2_p = \frac{df\_{\text{effect}}\\
F}{df\_{\text{effect}}\\ F + df\_{\text{error}}},\$\$ computed from the
pooled-error *F* and *df*. The confidence interval is built by Steiger's
(2004) transformation principle: a CI for the noncentrality parameter
\\\lambda\\ of the *F* distribution is obtained via
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md) and
then mapped through \\\eta^2_p = \lambda / (\lambda +
N\_{\text{ref}})\\, with \\N\_{\text{ref}}\\ taken to be the total study
*N* for the pooled error term (treating the simple effect as a contrast
within the full factorial design). When the lower limit on \\\lambda\\
is not identified (i.e., the observed *F* is below its one-sided
critical value), the lower limit on \\\eta^2_p\\ is set to 0; when the
upper limit is at infinity, the upper limit on \\\eta^2_p\\ is set to 1.

Under `error_term = "welch"`, partial eta squared and both limits are
`NA`. The transformation of Welch's statistic and fractional degrees of
freedom does not estimate ordinary partial eta squared, and transforming
a noncentral-F interval with the level's sample size targets a different
quantity under unequal variances. The Welch F reference distribution is
itself approximate. The tests remain available; a warning, the
`effect_size_inference` attribute, and the printed explanation identify
the unavailable effect-size inference.

**Multiplicity across the family.** With `which = "both"` the family is
the \\a + b\\ simple effects returned in a single call; the adjustment
is applied to that entire family. If only one direction is wanted, call
the function twice with `which = "A_at_B"` and `which = "B_at_A"` so
each family is adjusted on its own. The `"bonferroni"` adjustment is
\\p\_{\text{adj}} = \min(1, m\\ p)\\ for \\m\\ rows; the sequential
methods (`"holm"`, `"hochberg"`, `"BH"`, `"BY"`) are computed via
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

**When to perform simple effects.** A significant omnibus interaction
need not precede simple effect testing (Maxwell, Delaney, & Kelley,
2027). Simple effects are informative whenever the substantive question
is conditional on a level of the other factor, and the multiplicity
adjustment controls the family-wise error rate independently of any
interaction screen.

**Scope.** Only fixed-effects between-subjects designs with *exactly
two* crossed factors are supported. Within-subjects or mixed designs
(`aovlist` fits) and three- or higher-way designs are out of scope for
this function. Per-cell sample sizes may be unequal.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

Welch, B. L. (1951). On the comparison of several mean values: An
alternative approach. *Biometrika, 38*, 330–336.

## See also

[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)
for within-level pairwise or custom contrasts,
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md)
and
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)
for the omnibus effect size counterparts,
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md) for
the noncentrality machinery,
[`ss_power_factorial_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md)
for power calculations on the omnibus factorial effects.

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
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`obrien_test()`](https://yelleknek.github.io/DMAR/reference/obrien_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 2 x 3 factorial: expectancy treatment (A) x grade (B) on the
# pygmalion data. Grades 4 through 6 are omitted so the family of
# simple effects stays short enough to read at a glance.
pyg <- pygmalion[pygmalion$grade <= 3, ]
pyg$grade <- factor(pyg$grade)
fit <- aov(iq_gain ~ treatment * grade, data = pyg)

# Default: pooled MS_W, both families, no adjustment. The expectancy
# effect on IQ gain is concentrated in grades 1 and 2; at grade 3 the
# F is 0.004, so the lower limit on partial eta squared is clamped to
# 0 and the function notes the clamp in a warning.
simple_effects_AB(fit)
#> Warning: The ci_nc_F() lower-limit clamp fired in 1 of the simple effect rows (observed F below the alpha_lower critical value of the central F-distribution); the corresponding ci_lower on eta_squared_partial is clamped to 0. See ?ci_nc_F for the meaning of the clamp.
#>  effect                      focal_factor conditioning_factor
#>  treatment | grade = 1       treatment    grade              
#>  treatment | grade = 2       treatment    grade              
#>  treatment | grade = 3       treatment    grade              
#>  grade | treatment = Control grade        treatment          
#>  grade | treatment = Bloomer grade        treatment          
#>  conditioning_level F_value df_effect df_error p_value p_adjusted
#>  1                  7.388   1         157      0.0073  0.0073    
#>  2                  5.265   1         157      0.0231  0.0231    
#>  3                  0.004   1         157      0.9473  0.9473    
#>  Control            3.957   2         157      0.0211  0.0211    
#>  Bloomer            6.461   2         157      0.0020  0.0020    
#>  eta_squared_partial ci_lower ci_upper n_at_level
#>  0.045               0.003    0.119    52        
#>  0.032               0.000    0.100    58        
#>  0.000               0.000    0.009    53        
#>  0.048               0.001    0.117    131       
#>  0.076               0.011    0.155    32        
#> 
#> Confidence level: 95%
#> 
#> Pooled-error partial eta squared with a noncentral-F interval.

# Only the simple effects of grade within each treatment level, with
# a Holm adjustment across that family of two tests.
simple_effects_AB(fit, which = "B_at_A", adjust = "holm")
#>  effect                      focal_factor conditioning_factor
#>  grade | treatment = Control grade        treatment          
#>  grade | treatment = Bloomer grade        treatment          
#>  conditioning_level F_value df_effect df_error p_value p_adjusted
#>  Control            3.957   2         157      0.0211  0.0211    
#>  Bloomer            6.461   2         157      0.0020  0.0040    
#>  eta_squared_partial ci_lower ci_upper n_at_level
#>  0.048               0.001    0.117    131       
#>  0.076               0.011    0.155    32        
#> 
#> Confidence level: 95%
#> 
#> Pooled-error partial eta squared with a noncentral-F interval.

# Welch error term: refits a Welch one-way at each conditioning
# level. The Welch denominator df fall well below the pooled 157, so
# more of the lower limits are clamped to 0.
simple_effects_AB(fit, error_term = "welch")
#> Warning: Partial eta squared and its interval are unavailable for the Welch test. Transforming Welch's F and fractional degrees of freedom does not establish the ordinary variance-ratio estimand or a noncentral-F interval for it under unequal variances. The Welch F, degrees of freedom, and p-values remain available.
#>  effect                      focal_factor conditioning_factor
#>  treatment | grade = 1       treatment    grade              
#>  treatment | grade = 2       treatment    grade              
#>  treatment | grade = 3       treatment    grade              
#>  grade | treatment = Control grade        treatment          
#>  grade | treatment = Bloomer grade        treatment          
#>  conditioning_level F_value df_effect df_error p_value p_adjusted
#>  1                  6.868   1         9.027    0.0277  0.0277    
#>  2                  2.968   1         12.467   0.1096  0.1096    
#>  3                  0.007   1         24.514   0.9335  0.9335    
#>  Control            3.211   2         80.697   0.0455  0.0455    
#>  Bloomer            7.483   2         14.520   0.0058  0.0058    
#>  eta_squared_partial ci_lower ci_upper n_at_level
#>  <NA>                <NA>     <NA>     52        
#>  <NA>                <NA>     <NA>     58        
#>  <NA>                <NA>     <NA>     53        
#>  <NA>                <NA>     <NA>     131       
#>  <NA>                <NA>     <NA>     32        
#> 
#> Confidence level: 95%
#> 
#> Partial eta squared and its interval are unavailable for the Welch test. Transforming Welch's F and fractional degrees of freedom does not establish the ordinary variance-ratio estimand or a noncentral-F interval for it under unequal variances. The Welch F, degrees of freedom, and p-values remain available.

# Bonferroni across the full a + b = 5-test family.
simple_effects_AB(fit, adjust = "bonferroni")
#> Warning: The ci_nc_F() lower-limit clamp fired in 1 of the simple effect rows (observed F below the alpha_lower critical value of the central F-distribution); the corresponding ci_lower on eta_squared_partial is clamped to 0. See ?ci_nc_F for the meaning of the clamp.
#>  effect                      focal_factor conditioning_factor
#>  treatment | grade = 1       treatment    grade              
#>  treatment | grade = 2       treatment    grade              
#>  treatment | grade = 3       treatment    grade              
#>  grade | treatment = Control grade        treatment          
#>  grade | treatment = Bloomer grade        treatment          
#>  conditioning_level F_value df_effect df_error p_value p_adjusted
#>  1                  7.388   1         157      0.0073  0.0365    
#>  2                  5.265   1         157      0.0231  0.1154    
#>  3                  0.004   1         157      0.9473  1.0000    
#>  Control            3.957   2         157      0.0211  0.1053    
#>  Bloomer            6.461   2         157      0.0020  0.0101    
#>  eta_squared_partial ci_lower ci_upper n_at_level
#>  0.045               0.003    0.119    52        
#>  0.032               0.000    0.100    58        
#>  0.000               0.000    0.009    53        
#>  0.048               0.001    0.117    131       
#>  0.076               0.011    0.155    32        
#> 
#> Confidence level: 95%
#> 
#> Pooled-error partial eta squared with a noncentral-F interval.

# Simulated 2 x 2 design with a known interaction pattern.
set.seed(113)
d <- simulate_ancova_factorial_data(
  a = 2, b = 2,
  true_mean_Y = c(50, 60, 55, 50),   # crossover at B = 2
  true_mean_X = matrix(10, nrow = 4, ncol = 1),
  true_sd_Y = 8, true_sd_X = 3, true_cor_YX = 0,
  n_per_cell = 30
)
fit_sim <- aov(y ~ A * B, data = d)
simple_effects_AB(fit_sim, conf_level = 0.95)
#>  effect    focal_factor conditioning_factor conditioning_level F_value
#>  A | B = 1 A            B                   1                  31.563 
#>  A | B = 2 A            B                   2                  12.223 
#>  B | A = 1 B            A                   1                  20.424 
#>  B | A = 2 B            A                   2                  21.113 
#>  df_effect df_error p_value  p_adjusted eta_squared_partial ci_lower ci_upper
#>  1         116      < 0.0001 < 0.0001   0.214               0.094    0.330   
#>  1         116      0.0007   0.0007     0.095               0.018    0.201   
#>  1         116      < 0.0001 < 0.0001   0.150               0.048    0.264   
#>  1         116      < 0.0001 < 0.0001   0.154               0.051    0.268   
#>  n_at_level
#>  60        
#>  60        
#>  60        
#>  60        
#> 
#> Confidence level: 95%
#> 
#> Pooled-error partial eta squared with a noncentral-F interval.
```
