# Confidence Interval for Partial Eta Squared (Effect Size for ANOVA)

Computes the point estimate and an exact, noncentrality-based confidence
interval for the population *partial* eta squared (\\\eta^2_p\\).
Accepts either the raw ANOVA summary (*F*, effect df, error df, total
*N*) or a fitted `aov`/`lm`/`aovlist` object, in which case the function
returns one row per effect (with `stratum` identification for
within-subjects fits).

## Usage

``` r
ci_eta_squared_partial(
  object = NULL,
  F_value = NULL,
  df_effect = NULL,
  df_error = NULL,
  N = NULL,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL
)
```

## Arguments

- object:

  Optional. A fitted model object of class
  [`aov`](https://rdrr.io/r/stats/aov.html),
  [`lm`](https://rdrr.io/r/stats/lm.html), or `aovlist`.

- F_value:

  Observed *F*-value (ignored if `object` is supplied).

- df_effect:

  Numerator degrees of freedom for the effect (ignored if `object` is
  supplied).

- df_error:

  Error (residual) degrees of freedom (ignored if `object` is supplied).

- N:

  Total sample size (ignored if `object` is supplied).

- conf_level:

  Desired confidence coverage; default `0.95`.

- alpha_lower, alpha_upper:

  Optional Type I error on the lower and upper side.

## Value

A `data.frame` with one row per effect. Single-stratum fits and the raw
interface return columns `effect`, `eta_squared_partial`, `lower_limit`,
`upper_limit`, `F_value`, `df_effect`, `df_error`, `N`. `aovlist` fits
additionally include a `stratum` column. With the raw-argument interface
`effect` is `"overall"`.

## Details

This is the explicitly-named counterpart of
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md).
The two share point-estimate and CI machinery: in a one-way ANOVA
partial \\\eta^2\\ coincides with \\\eta^2\\; in a factorial or
within-subjects ANOVA both functions return the per-effect *partial*
value computed against that effect's own error stratum. Use
`ci_eta_squared_partial` when you want the function name to make the
partial interpretation explicit.

## References

Cohen, J. (1973). Eta-squared and partial eta-squared in fixed factor
ANOVA designs. *Educational and Psychological Measurement, 33*(1),
107–112.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological
Methods, 17*, 137–152.
[doi:10.1037/a0028086](https://doi.org/10.1037/a0028086)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\\eta^2\\, Chapter 7 on factorial
designs, and Chapter 11 on generalized \\\eta^2\\ for within-subjects
designs.)

Smithson, M. (2001). Correct confidence intervals for various regression
effect sizes and parameters: The importance of noncentral distributions
in computing intervals. *Educational and Psychological Measurement, 61*,
605–632.
[doi:10.1177/00131640121971392](https://doi.org/10.1177/00131640121971392)

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md),
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_cc()`](https://yelleknek.github.io/DMAR/reference/ci_cc.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_r.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_sm()`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Raw-argument interface.
ci_eta_squared_partial(F_value = 11.221, df_effect = 4,
                       df_error = 50, N = 55)
#>  effect  eta_squared_partial lower_limit upper_limit F_value df_effect df_error
#>  overall 0.473               0.226       0.587       11.2    4         50      
#>  N 
#>  55

# Factorial ANOVA: per-effect partial eta squared with CI.
fit <- aov(breaks ~ wool * tension, data = warpbreaks)
ci_eta_squared_partial(fit)
#> Warning: The observed F_value is below the alpha_lower critical value of the central F-distribution; the lower noncentrality limit has been clamped to 0 and the reported 'prob_greater' on the lower_limit row reflects the actual upper-tail probability at lambda = 0.
#>  effect       eta_squared_partial lower_limit upper_limit F_value df_effect
#>  wool         0.0727              0           0.222       3.77    1        
#>  tension      0.261               0.0558      0.411       8.5     2        
#>  wool:tension 0.149               0.00191     0.298       4.19    2        
#>  df_error N 
#>  48       54
#>  48       54
#>  48       54

# Within-subjects ANOVA.
set.seed(113)
n <- 20
rm_data <- data.frame(
  subject = factor(rep(seq_len(n), each = 3)),
  time    = factor(rep(c("Pre", "Mid", "Post"), n),
                   levels = c("Pre", "Mid", "Post")),
  y       = rnorm(n, sd = 1.5)[rep(seq_len(n), each = 3)] +
            0.7 * rep(1:3, n) + rnorm(n * 3, sd = 1.2)
)
fit_rm <- aov(y ~ time + Error(subject/time), data = rm_data)
ci_eta_squared_partial(fit_rm)
#>  effect eta_squared_partial lower_limit upper_limit stratum      F_value
#>  time   0.231               0.0152      0.322       subject:time 5.7    
#>  df_effect df_error N 
#>  2         38       60
```
