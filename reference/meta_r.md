# Random Effects Meta-Analysis of Correlations

Pools correlations across independent studies on the Fisher's *Z* scale
and reports the results back in the correlation metric. Optionally, each
study's correlation is first corrected for attenuation due to
measurement error in either or both variables (the Spearman correction
of
[`correction_for_attenuation`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
the basic artifact correction of Hunter and Schmidt's psychometric
meta-analysis), using reliabilities you supply, for example from the
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
family. The pooled quantity is then the correlation between the
constructs rather than the attenuated correlation between the observed
scores, and the weight of each study comes from the sampling variance of
its corrected correlation, which the correction enlarges.

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
  left `NULL` is treated as 1), and its sampling variance is that of the
  corrected correlation (see Details). The reliabilities are treated as
  known constants: the uncertainty of a reliability that was estimated
  is not carried into the results.

- method, hartung_knapp, conf_level:

  Passed to
  [`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md).
  Under the Hartung-Knapp adjustment, identical correlations give a zero
  adjusted standard error and a warning of class `dmar_meta_hk_zero_se`;
  see
  [`?meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md).
  In that warning the degenerate limits and the conventional normal
  interval are in the correlation metric, as the reported rows are, and
  the conventional standard error is on the Fisher's *Z* scale, where
  the model is fit, and is labeled so.

## Value

A `data.frame` (class `dmar_meta`, a `dmar_tbl`) with the same rows and
attributes as
[`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md): the
`estimate`, `ci_lower` / `ci_upper`, and prediction interval rows in the
correlation metric; the `se`, `test_statistic`, and heterogeneity rows
on the Fisher's *Z* scale where the model lives. The reference
distribution of `test_statistic` is in the `test_distribution` and
`test_df` attributes (Student's *t* on \\k - 1\\ degrees of freedom
under the Hartung-Knapp adjustment, the standard normal without it).
With reliabilities supplied, every row describes the corrected
correlations, with the sampling variances given in Details. When the
Hartung-Knapp adjustment is applied, the attributes `hartung_knapp_q`
and `se_conventional` hold the weighted residual mean square and the
conventional standard error, both on the Fisher's *Z* scale.

## Details

Without reliabilities, pooling uses \\z_i = \mathrm{atanh}(r_i)\\ with
sampling variance \\1 / (n_i - 3)\\; with them, the effect sizes and
their sampling variances are those given under Corrected correlations
below. Either way, the pooled estimate, its confidence limits, and the
prediction interval are transformed back through \\\tanh\\. The
heterogeneity quantities (tau, tau-squared, I-squared, H-squared, Q)
remain on the Fisher's *Z* scale, where the model lives; tau is
therefore the between-study standard deviation of the Fisher's *Z*
values.

**Corrected correlations.** With \\a_i = \rho\_{xx,i}\\ \rho\_{yy,i}\\,
the corrected correlation is \\r\_{c,i} = r_i / \sqrt{a_i}\\ and the
effect size pooled is \\\mathrm{atanh}(r\_{c,i})\\. The correction
divides the sampling error of \\r_i\\ by \\\sqrt{a_i}\\ along with
\\r_i\\ itself, so a corrected correlation is less precise than the
observed correlation it came from, and the Fisher's *Z* transformation,
which makes the sampling variance of \\\mathrm{atanh}(r_i)\\ free of the
population correlation, does not do the same for
\\\mathrm{atanh}(r\_{c,i})\\. The sampling variance used is the delta
method approximation \$\$v_i = \frac{(1 - r_i^2)^2}{a_i (1 -
r\_{c,i}^2)^2 (n_i - 3)},\$\$ which is \\1 / (n_i - 3)\\, the sampling
variance of \\\mathrm{atanh}(r_i)\\, times the square of the derivative
of \\\mathrm{atanh}(r\_{c,i})\\ with respect to \\\mathrm{atanh}(r_i)\\.
That derivative exceeds 1 whenever \\a_i \< 1\\, and it grows without
bound as \\r\_{c,i}\\ approaches 1 in magnitude. It is evaluated at each
study's observed correlation, so \\v_i\\ depends on that study alone, as
the model of
[`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md) takes
the sampling variances to be known. With no correction (\\a_i = 1\\),
\\v_i\\ is \\1 / (n_i - 3)\\ exactly. The \\v_i\\ enter wherever the
sampling variances do: the weights, the pooled estimate with its
standard error and confidence interval, the between-study variance with
its interval, I-squared, H-squared, Q, and the prediction interval. For
five studies with \\r = 0.40\\, \\n = 100\\, and both reliabilities
0.80, each corrected correlation is 0.50 and each \\v_i\\ is 1.96 times
\\1 / (n_i - 3)\\; the fixed effect interval for the correlation between
the constructs is 0.401 to 0.588.

The reliabilities are treated as known constants. When they are
estimates, from the same samples or from other sources, their sampling
error is a further source of uncertainty in each corrected correlation
that \\v_i\\ does not include, and the intervals are then narrower than
the total uncertainty warrants, by an amount that depends on how
precisely the reliabilities are known. Because \\v_i\\ grows with the
magnitude of the study's own correlation, a study whose corrected
correlation is larger in magnitude receives less weight, which moves the
pooled estimate toward zero. In simulations run for this function (10
studies, a correlation of 0.5 between the constructs, reliabilities from
0.60 to 0.95, 30 to 400 observations per study, and 4,000 sets of
studies per condition, so that a coverage rate near 0.95 has a Monte
Carlo standard error of about 0.0035), the pooled estimate was, on
average, about 0.003 to 0.015 below the construct correlation on the
Fisher's *Z* scale, depending on the condition. When the effects were
homogeneous, the fixed effect interval covered at the nominal 0.95
within Monte Carlo error (0.950 to 0.955; 0.80 to 0.83 with the variance
of an uncorrected correlation). With heterogeneous effects (a
between-study standard deviation of 0.15 on the Fisher's *Z* scale of
the construct correlations), the interval for tau-squared covered at
0.95 to 0.97 (0.79 to 0.85 with the variance of an uncorrected
correlation) and the Hartung-Knapp interval at 0.93 to 0.95. The delta
method approximation is poorer when \\n_i\\ is small or a corrected
correlation is near 1 in magnitude: for a correlation of 0.7 between the
constructs, reliabilities of 0.7, and 60 observations per study, a
corrected correlation reached 1 in 7 percent of the sets of studies, and
in the rest the pooled estimate was, on average, about 0.08 below the
construct correlation on the Fisher's *Z* scale and the intervals
covered at 0.84 to 0.87. The artifact distribution methods of Hunter and
Schmidt (2004) are not implemented. A corrected correlation that reaches
or exceeds 1 in magnitude (possible when an observed \\r\\ outruns the
supplied reliabilities) is an error at the pooling stage, unlike the
single-study
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
[`convert_r_Z`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md)
/
[`convert_Z_r`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md)
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
#>  estimate         0.290 
#>  se               0.032 
#>  test_statistic   9.337 
#>  p_value          0.0007
#>  ci_lower         0.207 
#>  ci_upper         0.369 
#>  prediction_lower 0.194 
#>  prediction_upper 0.380 
#>  tau2             0.000 
#>  tau2_ci_lower    0.000 
#>  tau2_ci_upper    0.037 
#>  tau              0.000 
#>  I2               0.000 
#>  I2_ci_lower      0.000 
#>  I2_ci_upper      80.877
#>  H2               1.000 
#>  Q                2.449 
#>  Q_df             4     
#>  Q_p              0.6538
#>  n_studies        5     
#> 
#> Confidence level: 95%
#> 
#> test_statistic: referred to Student's t on 4 df (Hartung-Knapp adjustment)

# The same studies corrected for criterion unreliability (reliability
# 0.80 in every study): the construct-level validity. Its interval is
# wider than the one above: the correction spreads the corrected
# correlations farther apart and enlarges their sampling variances.
meta_r(r, n, reliability_y = 0.80)
#>  term             value 
#>  estimate         0.323 
#>  se               0.036 
#>  test_statistic   9.185 
#>  p_value          0.0008
#>  ci_lower         0.229 
#>  ci_upper         0.410 
#>  prediction_lower 0.215 
#>  prediction_upper 0.423 
#>  tau2             0.000 
#>  tau2_ci_lower    0.000 
#>  tau2_ci_upper    0.048 
#>  tau              0.000 
#>  I2               0.000 
#>  I2_ci_lower      0.000 
#>  I2_ci_upper      80.837
#>  H2               1.000 
#>  Q                2.430 
#>  Q_df             4     
#>  Q_p              0.6572
#>  n_studies        5     
#> 
#> Confidence level: 95%
#> 
#> test_statistic: referred to Student's t on 4 df (Hartung-Knapp adjustment)
```
