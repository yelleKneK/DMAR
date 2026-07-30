# Random Effects Meta-Analysis of Standardized Mean Differences

Pools two-group standardized mean differences across independent
studies. Each study contributes its standardized mean difference and
per-group sample sizes; the function computes the within-study sampling
variances, applies the Hedges small-sample bias correction by default
(the same \\J\\ factor as
[`expected_smd`](https://yelleknek.github.io/DMAR/reference/expected_smd.md)
and [`smd`](https://yelleknek.github.io/DMAR/reference/smd.md)), and
fits the random effects model of
[`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md),
returning the pooled effect with its confidence interval, tau and
tau-squared with intervals, I-squared, Cochran's Q, and a prediction
interval for the effect in a new study.

## Usage

``` r
meta_smd(
  smd,
  n_1,
  n_2,
  unbiased = TRUE,
  method = c("reml", "pm", "dl", "fe"),
  hartung_knapp = TRUE,
  conf_level = 0.95
)
```

## Arguments

- smd:

  Numeric vector of standardized mean differences (Cohen's *d*), one per
  study, positive in the direction of the common hypothesis.

- n_1, n_2:

  Per-group sample sizes for each study.

- unbiased:

  Logical: convert each *d* to Hedges *g* (the small-sample unbiased
  estimator) before pooling? Default `TRUE`. Set `FALSE` to pool the raw
  *d* values, for example when reproducing a historical analysis such as
  Raudenbush (1984) that predates routine use of the correction.

- method, hartung_knapp, conf_level:

  Passed to
  [`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md):
  the tau-squared estimator (`"reml"` default), the Hartung-Knapp
  small-sample adjustment (default `TRUE`), and the confidence level.

## Value

A `data.frame` (class `dmar_tbl`) with the same rows as
[`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md).

## Details

The within-study variance is the standard large-sample form \$\$v_i =
\frac{n\_{1i} + n\_{2i}}{n\_{1i} n\_{2i}} + \frac{g_i^2}{2 (n\_{1i} +
n\_{2i})},\$\$ computed from the bias-corrected \\g_i\\ when
`unbiased = TRUE` (Hedges, 1981; Borenstein, Hedges, Higgins, &
Rothstein, 2009). All reported quantities are in the standardized mean
difference metric.

## References

Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R.
(2009). *Introduction to meta-analysis*. Wiley.

Hedges, L. V. (1981). Distribution theory for Glass's estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128.

Raudenbush, S. W. (1984). Magnitude of teacher expectancy effects on
pupil IQ as a function of the credibility of expectancy induction: A
synthesis of findings from 18 experiments. *Journal of Educational
Psychology, 76*(1), 85–97.

## See also

[`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md) for
the engine and the reported rows;
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) and
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) for the
single-study quantities;
[`plot_forest`](https://yelleknek.github.io/DMAR/reference/plot_forest.md)
for the picture;
[`teacher_expectancy`](https://yelleknek.github.io/DMAR/reference/teacher_expectancy.md)
for the example data.

Other meta-analysis:
[`combine_p()`](https://yelleknek.github.io/DMAR/reference/combine_p.md),
[`meta_contrast()`](https://yelleknek.github.io/DMAR/reference/meta_contrast.md),
[`meta_es()`](https://yelleknek.github.io/DMAR/reference/meta_es.md),
[`meta_r()`](https://yelleknek.github.io/DMAR/reference/meta_r.md),
[`plot_forest()`](https://yelleknek.github.io/DMAR/reference/plot_forest.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Pool the teacher expectancy studies (Raudenbush, 1984). Hedges g and
# the Hartung-Knapp adjustment are on by default; the prediction
# interval shows where a new expectancy study would be expected to land.
data(teacher_expectancy)
meta_smd(smd = teacher_expectancy$d,
         n_1 = teacher_expectancy$n_experimental,
         n_2 = teacher_expectancy$n_control)
#>  term             value  
#>  estimate         0.0544 
#>  se               0.0352 
#>  t                1.55   
#>  p_value          0.1392 
#>  lower_limit      -0.0195
#>  upper_limit      0.128  
#>  prediction_lower -0.0198
#>  prediction_upper 0.129  
#>  tau2             0      
#>  tau2_lower       0      
#>  tau2_upper       0.048  
#>  tau              0      
#>  I2               0      
#>  I2_lower         0      
#>  I2_upper         64.6   
#>  H2               1      
#>  Q                16.7   
#>  Q_df             18     
#>  Q_p              0.5468 
#>  k                19     
#> 
#> Confidence level: 95%
```
