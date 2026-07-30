# Random Effects Meta-Analysis of Correlations

Pools correlations across independent studies on the Fisher *z* scale
and reports the results back in the correlation metric. Optionally, each
study's correlation is first corrected for attenuation due to
measurement error in either or both variables (the Spearman correction
of
[`correction_for_attenuation`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
the basic artifact correction of Hunter and Schmidt's psychometric
meta-analysis), using reliabilities you supply, for example from the
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
family. That combination, synthesis connected to an actual reliability
toolkit, is the measurement-aware path: the pooled quantity is then the
construct-level correlation rather than the attenuated observed one.

## Usage

``` r
meta_r(
  r,
  n,
  reliability_x = NULL,
  reliability_y = NULL,
  method = c("reml", "pm", "dl", "fe"),
  hartung_knapp = TRUE,
  conf_level = 0.95
)
```

## Arguments

- r:

  Numeric vector of observed correlations, one per study, each in (-1,
  1).

- n:

  Per-study sample sizes (integer, at least 4).

- reliability_x, reliability_y:

  Optional per-study reliabilities in (0, 1\] for the two measured
  variables; a single value is recycled across studies. When either is
  supplied, each correlation is disattenuated by \\r_i /
  \sqrt{\rho\_{xx,i}\\ \rho\_{yy,i}}\\ before pooling (a reliability
  left `NULL` is treated as 1). The reliabilities are treated as known.

- method, hartung_knapp, conf_level:

  Passed to
  [`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md).

## Value

A `data.frame` (class `dmar_tbl`) with the same rows as
[`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md): the
`estimate`, `lower_limit` / `upper_limit`, and prediction interval rows
in the correlation metric; the `se`, test statistic, and heterogeneity
rows on the Fisher *z* scale where the model lives.

## Details

Pooling uses \\z_i = \mathrm{atanh}(r_i)\\ with sampling variance \\1 /
(n_i - 3)\\; the pooled estimate, its confidence limits, and the
prediction interval are transformed back through \\\tanh\\. The
heterogeneity quantities (tau, tau-squared, I-squared, H-squared, Q)
remain on the Fisher *z* scale, where the model lives; tau is therefore
the between-study standard deviation of the *z*-scale correlations.

When corrections are applied, the corrected correlation's variance is
computed from its own \\n_i\\ on the *z* scale, the conventional simple
treatment when reliabilities are taken as known constants; the more
elaborate artifact-distribution machinery of Hunter and Schmidt (2004)
is deliberately out of scope here. A corrected correlation that exceeds
1 in magnitude (possible when an observed \\r\\ outruns the supplied
reliabilities) is an error at the pooling stage, unlike the single-study
[`correction_for_attenuation`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
which reports it with a warning: \\\mathrm{atanh}\\ is undefined there.

## References

Hunter, J. E., & Schmidt, F. L. (2004). *Methods of meta-analysis:
Correcting error and bias in research findings* (2nd ed.). Sage.

## See also

[`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md) for
the engine;
[`correction_for_attenuation`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md)
for the single-study correction and its connection to latent variable
modeling;
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
for estimating the reliabilities;
[`convert_r_z`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md)
/
[`convert_z_r`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md)
for the transformation used.

Other meta-analysis:
[`combine_p()`](https://yelleknek.github.io/DMAR/reference/combine_p.md),
[`meta_contrast()`](https://yelleknek.github.io/DMAR/reference/meta_contrast.md),
[`meta_es()`](https://yelleknek.github.io/DMAR/reference/meta_es.md),
[`meta_smd()`](https://yelleknek.github.io/DMAR/reference/meta_smd.md),
[`plot_forest()`](https://yelleknek.github.io/DMAR/reference/plot_forest.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Five validity studies of the same selection instrument.
r <- c(.28, .35, .22, .40, .31)
n <- c(120, 85, 200, 60, 150)
meta_r(r, n)
#>  term             value 
#>  estimate         0.29  
#>  se               0.0319
#>  t                9.34  
#>  p_value          0.0007
#>  lower_limit      0.207 
#>  upper_limit      0.369 
#>  prediction_lower 0.194 
#>  prediction_upper 0.38  
#>  tau2             0     
#>  tau2_lower       0     
#>  tau2_upper       0.0368
#>  tau              0     
#>  I2               0     
#>  I2_lower         0     
#>  I2_upper         80.9  
#>  H2               1     
#>  Q                2.45  
#>  Q_df             4     
#>  Q_p              0.6538
#>  k                5     
#> 
#> Confidence level: 95%

# The same studies corrected for criterion unreliability (reliability
# 0.80 in every study): the construct-level validity.
meta_r(r, n, reliability_y = 0.80)
#>  term             value 
#>  estimate         0.324 
#>  se               0.0367
#>  t                9.18  
#>  p_value          0.0008
#>  lower_limit      0.23  
#>  upper_limit      0.412 
#>  prediction_lower 0.216 
#>  prediction_upper 0.424 
#>  tau2             0     
#>  tau2_lower       0     
#>  tau2_upper       0.052 
#>  tau              0     
#>  I2               0     
#>  I2_lower         0     
#>  I2_upper         85.7  
#>  H2               1     
#>  Q                3.22  
#>  Q_df             4     
#>  Q_p              0.5209
#>  k                5     
#> 
#> Confidence level: 95%
```
