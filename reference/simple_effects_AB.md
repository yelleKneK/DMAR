# Simple Effect F Tests for a Two-Factor Between-Subjects Design

Given a fitted [`aov`](https://rdrr.io/r/stats/aov.html) or
[`lm`](https://rdrr.io/r/stats/lm.html) object for a two-factor
between-subjects design, conventionally written \\Y \sim A \times B\\,
where \\A\\ and \\B\\ are crossed fixed factors, computes the family of
*simple main effects*: the effect of \\A\\ at each level of \\B\\ and/or
the effect of \\B\\ at each level of \\A\\. Each row carries the simple
effect *F* test, its (optionally adjusted) *p*-value, and the partial
\\\eta^2\\ with a noncentrality-based confidence interval. The error
term can be either the full-model pooled \\\mathit{MS}\_W\\ (the
textbook default) or a Welch–Satterthwaite test refitted within each
conditioning level (robust to within-level heteroscedasticity).

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
  has *exactly two* crossed factors (e.g.\\ `breaks ~ wool * tension`).
  The interaction term is strongly recommended so that the pooled error
  is the pure within-cell \\\mathit{MS}\_W\\; the function still runs
  without it but issues a warning (the additive-model residual includes
  interaction variance and inflates the error term used for the pooled
  simple effect *F*).

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
`df_effect`, `df_error`, `p_value`, `p_adjusted`, `partial_eta_squared`,
`lower_limit`, `upper_limit`, `n_at_level`. Attributes `error_term`,
`adjust`, `conf_level`, `factor_A`, and `factor_B` record the call
options.

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
*F* and the *df* actually used in the test (so it reflects whichever
error term was chosen). The confidence interval is built by Steiger's
(2004) transformation principle: a CI for the noncentrality parameter
\\\lambda\\ of the *F* distribution is obtained via
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
and then mapped through \\\eta^2_p = \lambda / (\lambda +
N\_{\text{ref}})\\, with \\N\_{\text{ref}}\\ taken to be the total study
*N* for the pooled error term (treating the simple effect as a contrast
within the full factorial design) and the level-conditional sample size
\\n\_{\|b_j}\\ for the Welch error term (since the Welch test uses only
those observations). When the lower limit on \\\lambda\\ is not
identified (i.e., the observed *F* is below its one-sided critical
value), the lower limit on \\\eta^2_p\\ is set to 0; when the upper
limit is at infinity, the upper limit on \\\eta^2_p\\ is set to 1.

**Multiplicity across the family.** With `which = "both"` the family is
the \\a + b\\ simple effects returned in a single call; the adjustment
is applied to that entire family. If only one direction is wanted, call
the function twice with `which = "A_at_B"` and `which = "B_at_A"` so
each family is adjusted on its own. The `"bonferroni"` adjustment is
\\p\_{\text{adj}} = \min(1, m\\ p)\\ for \\m\\ rows; the sequential
methods (`"holm"`, `"hochberg"`, `"BH"`, `"BY"`) are computed via
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

**When to perform simple effects.** It is no longer required that a
significant omnibus interaction precede simple effect testing (Maxwell,
Delaney, & Kelley, 2027). Simple effects are informative whenever the
substantive question is conditional on a level of the other factor, and
the multiplicity adjustment controls the family-wise error rate
independently of any interaction screen.

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
[`eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md)
and
[`ci_eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md)
for the omnibus effect size counterparts,
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
for the noncentrality machinery,
[`ss_power_factorial_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md)
for power calculations on the omnibus factorial effects.

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
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`scheffe_ci()`](https://yelleknek.github.io/DMAR/reference/scheffe_ci.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`tost_r()`](https://yelleknek.github.io/DMAR/reference/tost_r.md),
[`tost_smd()`](https://yelleknek.github.io/DMAR/reference/tost_smd.md),
[`tukey_kramer_ci()`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 2 x 3 factorial: wool (A) x tension (B) on warpbreaks.
fit <- aov(breaks ~ wool * tension, data = warpbreaks)

# Default: pooled MS_W, both families, no adjustment.
simple_effects_AB(fit)
#> Warning: The conf_limits_ncf() lower-limit clamp fired in 3 of the simple effect rows (observed F below the alpha_lower critical value of the central F-distribution); the corresponding lower_limit on partial_eta_squared is clamped to 0. See ?conf_limits_ncf for the meaning of the clamp.
#>  effect             focal_factor conditioning_factor conditioning_level F_value
#>  wool | tension = L wool         tension             L                  10     
#>  wool | tension = M wool         tension             M                  0.858  
#>  wool | tension = H wool         tension             H                  1.26   
#>  tension | wool = A tension      wool                A                  10.3   
#>  tension | wool = B tension      wool                B                  2.37   
#>  df_effect df_error p_value p_adjusted partial_eta_squared lower_limit
#>  1         48       0.0027  0.0027     0.173               0.0216     
#>  1         48       0.3589  0.3589     0.0176              0          
#>  1         48       0.2682  0.2682     0.0255              0          
#>  2         48       0.0002  0.0002     0.301               0.082      
#>  2         48       0.1039  0.1039     0.09                0          
#>  upper_limit n_at_level
#>  0.335       18        
#>  0.134       18        
#>  0.15        18        
#>  0.446       27        
#>  0.229       27        
#> 
#> Confidence level: 95%

# Only the simple effects of tension within each wool level, with a
# Holm adjustment across that family of two tests.
simple_effects_AB(fit, which = "B_at_A", adjust = "holm")
#> Warning: The conf_limits_ncf() lower-limit clamp fired in 1 of the simple effect rows (observed F below the alpha_lower critical value of the central F-distribution); the corresponding lower_limit on partial_eta_squared is clamped to 0. See ?conf_limits_ncf for the meaning of the clamp.
#>  effect             focal_factor conditioning_factor conditioning_level F_value
#>  tension | wool = A tension      wool                A                  10.3   
#>  tension | wool = B tension      wool                B                  2.37   
#>  df_effect df_error p_value p_adjusted partial_eta_squared lower_limit
#>  2         48       0.0002  0.0004     0.301               0.082      
#>  2         48       0.1039  0.1039     0.09                0          
#>  upper_limit n_at_level
#>  0.446       27        
#>  0.229       27        
#> 
#> Confidence level: 95%

# Welch error term: refits a Welch one-way at each conditioning level.
simple_effects_AB(fit, error_term = "welch")
#> Warning: The conf_limits_ncf() lower-limit clamp fired in 3 of the simple effect rows (observed F below the alpha_lower critical value of the central F-distribution); the corresponding lower_limit on partial_eta_squared is clamped to 0. See ?conf_limits_ncf for the meaning of the clamp.
#>  effect             focal_factor conditioning_factor conditioning_level F_value
#>  wool | tension = L wool         tension             L                  5.65   
#>  wool | tension = M wool         tension             M                  1.25   
#>  wool | tension = H wool         tension             H                  2.32   
#>  tension | wool = A tension      wool                A                  4.8    
#>  tension | wool = B tension      wool                B                  5.8    
#>  df_effect df_error p_value p_adjusted partial_eta_squared lower_limit
#>  1         12.4     0.0344  0.0344     0.314               0          
#>  1         15.9     0.2796  0.2796     0.0731              0          
#>  1         11.5     0.1548  0.1548     0.168               0          
#>  2         15.1     0.0243  0.0243     0.389               0.000678   
#>  2         14.3     0.0144  0.0144     0.448               0.0155     
#>  upper_limit n_at_level
#>  0.53        18        
#>  0.348       18        
#>  0.412       18        
#>  0.499       27        
#>  0.535       27        
#> 
#> Confidence level: 95%

# Bonferroni across the full a + b = 5-test family.
simple_effects_AB(fit, adjust = "bonferroni")
#> Warning: The conf_limits_ncf() lower-limit clamp fired in 3 of the simple effect rows (observed F below the alpha_lower critical value of the central F-distribution); the corresponding lower_limit on partial_eta_squared is clamped to 0. See ?conf_limits_ncf for the meaning of the clamp.
#>  effect             focal_factor conditioning_factor conditioning_level F_value
#>  wool | tension = L wool         tension             L                  10     
#>  wool | tension = M wool         tension             M                  0.858  
#>  wool | tension = H wool         tension             H                  1.26   
#>  tension | wool = A tension      wool                A                  10.3   
#>  tension | wool = B tension      wool                B                  2.37   
#>  df_effect df_error p_value p_adjusted partial_eta_squared lower_limit
#>  1         48       0.0027  0.0134     0.173               0.0216     
#>  1         48       0.3589  1.0000     0.0176              0          
#>  1         48       0.2682  1.0000     0.0255              0          
#>  2         48       0.0002  0.0009     0.301               0.082      
#>  2         48       0.1039  0.5193     0.09                0          
#>  upper_limit n_at_level
#>  0.335       18        
#>  0.134       18        
#>  0.15        18        
#>  0.446       27        
#>  0.229       27        
#> 
#> Confidence level: 95%

# Simulated 2 x 2 design with a known interaction pattern.
set.seed(113)
d <- simulate_ancova_factorial_data(
  a = 2, b = 2,
  mu_y    = c(50, 60, 55, 50),   # crossover at B = 2
  mu_x    = matrix(10, nrow = 4, ncol = 1),
  sigma_y = 8, sigma_x = 3, rho_y_x = 0,
  n       = 30
)
fit_sim <- aov(y ~ A * B, data = d)
simple_effects_AB(fit_sim, conf_level = 0.95)
#>  effect    focal_factor conditioning_factor conditioning_level F_value
#>  A | B = 1 A            B                   1                  31.6   
#>  A | B = 2 A            B                   2                  12.2   
#>  B | A = 1 B            A                   1                  20.4   
#>  B | A = 2 B            A                   2                  21.1   
#>  df_effect df_error p_value  p_adjusted partial_eta_squared lower_limit
#>  1         116      < 0.0001 < 0.0001   0.214               0.0935     
#>  1         116      0.0007   0.0007     0.0953              0.0179     
#>  1         116      < 0.0001 < 0.0001   0.15                0.0482     
#>  1         116      < 0.0001 < 0.0001   0.154               0.051      
#>  upper_limit n_at_level
#>  0.33        60        
#>  0.201       60        
#>  0.264       60        
#>  0.268       60        
#> 
#> Confidence level: 95%
```
