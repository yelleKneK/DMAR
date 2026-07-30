# Random Effects Meta-Analysis of Generic Effect Sizes

Pools independent effect sizes given their sampling variances: the
general engine behind
[`meta_smd`](https://yelleknek.github.io/DMAR/reference/meta_smd.md) and
[`meta_r`](https://yelleknek.github.io/DMAR/reference/meta_r.md),
exposed for any effect metric whose estimates are approximately normal
with known variances. The random effects model is the default and the
fit reports the full uncertainty picture in one table: the pooled
estimate with its confidence interval, the between-study standard
deviation tau and variance tau-squared *with Q-profile confidence
intervals*, I-squared and H-squared *with intervals*, Cochran's Q test,
and, always, a prediction interval for the effect in a new study.
Reporting the prediction interval by default is deliberate: when
heterogeneity is real, the confidence interval for the average effect
understates what the next study will show, and the package treats “where
will the next study land” as part of the answer, not an option.

## Usage

``` r
meta_es(
  yi,
  vi,
  method = c("reml", "pm", "dl", "fe"),
  hartung_knapp = TRUE,
  conf_level = 0.95
)
```

## Arguments

- yi:

  Numeric vector of effect sizes, one per independent study.

- vi:

  Sampling variances of `yi`, one per study.

- method:

  Between-study variance estimator: `"reml"` (restricted maximum
  likelihood, the default), `"pm"` (Paule-Mandel), `"dl"`
  (DerSimonian-Laird), or `"fe"` (a fixed effect / common effect
  analysis, which assumes tau-squared is zero and reports no prediction
  interval).

- hartung_knapp:

  Logical: apply the Hartung-Knapp-Sidik-Jonkman small-sample adjustment
  (the pooled standard error rescaled from the weighted residuals, with
  a *t* reference on \\k - 1\\ degrees of freedom)? Default `TRUE`: with
  the small numbers of studies typical in psychology and education it
  keeps the confidence interval near its nominal coverage, where the
  conventional normal interval is anticonservative. Ignored for
  `method = "fe"`.

- conf_level:

  Confidence level for all intervals. Defaults to 0.95.

## Value

A `data.frame` (class `dmar_tbl`) with rows `estimate`, `se`, the test
statistic (`t` under Hartung-Knapp, `z` otherwise), `p_value`,
`lower_limit` / `upper_limit`, `prediction_lower` / `prediction_upper`,
`tau2` with `tau2_lower` / `tau2_upper`, `tau`, `I2` with limits, `H2`,
`Q` / `Q_df` / `Q_p`, and `k`. The estimator and adjustment are recorded
in the `"method"` and `"hartung_knapp"` attributes.

## Details

The model is \\y_i = \mu + u_i + e_i\\ with \\u_i \sim N(0, \tau^2)\\
and \\e_i \sim N(0, v_i)\\, \\v_i\\ treated as known. The \\\tau^2\\
confidence interval inverts the generalized Q statistic (Viechtbauer,
2007); the I-squared and H-squared intervals map the \\\tau^2\\ interval
through the typical within-study variance of Higgins and Thompson
(2002). The prediction interval follows Higgins, Thompson, and
Spiegelhalter (2009), using *t* with \\k - 2\\ degrees of freedom, and
requires at least three studies.

I-squared is reported because readers expect it, but note its well-known
limitation: it is a *proportion* of variability, not an amount, so the
same tau matched with larger studies yields a larger I-squared. The
quantity with direct scientific meaning is tau (the between-study
standard deviation, in the metric of `yi`) together with the prediction
interval.

## References

DerSimonian, R., & Laird, N. (1986). Meta-analysis in clinical trials.
*Controlled Clinical Trials, 7*(3), 177–188.

Hartung, J., & Knapp, G. (2001). On tests of the overall treatment
effect in meta-analysis with normally distributed responses. *Statistics
in Medicine, 20*(12), 1771–1782.
[doi:10.1002/sim.791](https://doi.org/10.1002/sim.791)

Higgins, J. P. T., & Thompson, S. G. (2002). Quantifying heterogeneity
in a meta-analysis. *Statistics in Medicine, 21*(11), 1539–1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186)

Higgins, J. P. T., Thompson, S. G., & Spiegelhalter, D. J. (2009). A
re-evaluation of random-effects meta-analysis. *Journal of the Royal
Statistical Society: Series A, 172*(1), 137–159.
[doi:10.1111/j.1467-985X.2008.00552.x](https://doi.org/10.1111/j.1467-985X.2008.00552.x)

Viechtbauer, W. (2007). Confidence intervals for the amount of
heterogeneity in meta-analysis. *Statistics in Medicine, 26*(1), 37–52.
[doi:10.1002/sim.2514](https://doi.org/10.1002/sim.2514)

## See also

[`meta_smd`](https://yelleknek.github.io/DMAR/reference/meta_smd.md) and
[`meta_r`](https://yelleknek.github.io/DMAR/reference/meta_r.md) for the
metric- specific front ends;
[`meta_contrast`](https://yelleknek.github.io/DMAR/reference/meta_contrast.md)
for focused moderator contrasts;
[`combine_p`](https://yelleknek.github.io/DMAR/reference/combine_p.md)
for combined significance tests;
[`plot_forest`](https://yelleknek.github.io/DMAR/reference/plot_forest.md)
to see the studies and the pool together.

Other meta-analysis:
[`combine_p()`](https://yelleknek.github.io/DMAR/reference/combine_p.md),
[`meta_contrast()`](https://yelleknek.github.io/DMAR/reference/meta_contrast.md),
[`meta_r()`](https://yelleknek.github.io/DMAR/reference/meta_r.md),
[`meta_smd()`](https://yelleknek.github.io/DMAR/reference/meta_smd.md),
[`plot_forest()`](https://yelleknek.github.io/DMAR/reference/plot_forest.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The teacher expectancy studies (Raudenbush, 1984), pooled in the d
# metric with variances from the standard large-sample formula.
data(teacher_expectancy)
d <- teacher_expectancy$d
n_e <- teacher_expectancy$n_experimental
n_c <- teacher_expectancy$n_control
v <- (n_e + n_c) / (n_e * n_c) + d^2 / (2 * (n_e + n_c))
meta_es(d, v)
#>  term             value  
#>  estimate         0.0549 
#>  se               0.0356 
#>  t                1.54   
#>  p_value          0.1401 
#>  lower_limit      -0.0198
#>  upper_limit      0.13   
#>  prediction_lower -0.0201
#>  prediction_upper 0.13   
#>  tau2             0      
#>  tau2_lower       0      
#>  tau2_upper       0.0513 
#>  tau              0      
#>  I2               0      
#>  I2_lower         0      
#>  I2_upper         66.1   
#>  H2               1      
#>  Q                17.1   
#>  Q_df             18     
#>  Q_p              0.5185 
#>  k                19     
#> 
#> Confidence level: 95%

# A fixed effect (common effect) analysis of the same studies.
meta_es(d, v, method = "fe")
#>  term             value  
#>  estimate         0.0549 
#>  se               0.0365 
#>  z                1.5    
#>  p_value          0.1328 
#>  lower_limit      -0.0167
#>  upper_limit      0.127  
#>  prediction_lower <NA>   
#>  prediction_upper <NA>   
#>  tau2             0      
#>  tau2_lower       <NA>   
#>  tau2_upper       <NA>   
#>  tau              0      
#>  I2               0      
#>  I2_lower         <NA>   
#>  I2_upper         <NA>   
#>  H2               1      
#>  Q                17.1   
#>  Q_df             18     
#>  Q_p              0.5185 
#>  k                19     
#> 
#> Confidence level: 95%
```
