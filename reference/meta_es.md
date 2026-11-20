# Random Effects Meta-Analysis of Generic Effect Sizes

Pools independent effect sizes given their sampling variances: the
general engine behind
[`meta_smd`](https://yelleknek.github.io/DMAR/reference/meta_smd.md) and
[`meta_r`](https://yelleknek.github.io/DMAR/reference/meta_r.md),
exposed for any effect metric whose estimates are approximately normal
with known variances. The random effects model is the default and the
fit reports the full uncertainty picture in one table: the pooled
estimate with its confidence interval, the between-study standard
deviation tau, the between-study variance tau-squared with its Q-profile
confidence interval, I-squared with an interval mapped from the
tau-squared limits, H-squared, Cochran's Q test, and, always, a
prediction interval for the effect in a new study. Reporting the
prediction interval by default is deliberate: when heterogeneity is
real, the confidence interval for the average effect understates what
the next study will show, and the package treats “where will the next
study land” as part of the answer, not an option.

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
  a *t* reference on \\k - 1\\ degrees of freedom)? Default `TRUE`. The
  adjusted standard error is not bounded below by the conventional one,
  so the interval can be narrower than the normal interval; see *The
  Hartung-Knapp adjustment* and *When q is small* in Details. Identical
  effect sizes give a zero adjusted standard error and a warning; see
  *Identical effect sizes*. Ignored for `method = "fe"`.

- conf_level:

  Confidence level for all intervals. Defaults to 0.95.

## Value

A `data.frame` (class `dmar_meta`, a `dmar_tbl`) with rows `estimate`,
`se`, `test_statistic` (the estimate divided by `se`), `p_value`,
`ci_lower` / `ci_upper`, `prediction_lower` / `prediction_upper`, `tau2`
with `tau2_ci_lower` / `tau2_ci_upper`, `tau`, `I2` with limits, `H2`,
`Q` / `Q_df` / `Q_p`, and `n_studies` (the number of studies, \\k\\).
The estimator and adjustment are recorded in the `"method"` and
`"hartung_knapp"` attributes. The reference distribution of
`test_statistic` and `p_value` depends on the adjustment, and the
attributes `test_distribution` and `test_df` record it: `"t"` on \\k -
1\\ degrees of freedom under the Hartung-Knapp adjustment, `"normal"`
(with `test_df` `NA`) without it, as the line under the printed table
says. The attributes `hartung_knapp_q` (the weighted residual mean
square \\q\\) and `se_conventional` are recorded whenever the adjustment
is applied. When the effect sizes are identical and the Hartung-Knapp
adjustment is applied, `se` is zero, the limits equal the estimate, and
the function warns (see *Identical effect sizes* in Details).

## Details

The model is \\y_i = \mu + u_i + e_i\\ with \\u_i \sim N(0, \tau^2)\\
and \\e_i \sim N(0, v_i)\\, \\v_i\\ treated as known. The \\\tau^2\\
confidence interval inverts the generalized Q statistic (Viechtbauer,
2007); the I-squared interval maps the \\\tau^2\\ interval through the
typical within-study variance of Higgins and Thompson (2002). The
prediction interval follows Higgins, Thompson, and Spiegelhalter (2009),
using *t* with \\k - 2\\ degrees of freedom, and requires at least three
studies.

**The Hartung-Knapp adjustment.** With weights \\w_i = 1/(v_i +
\hat{\tau}^2)\\ and pooled estimate \\\hat{\mu}\\, the adjusted standard
error is the conventional one, \\(\sum w_i)^{-1/2}\\, multiplied by
\\\sqrt{q}\\, where \$\$q = \frac{1}{k - 1} \sum\_{i=1}^{k} w_i (y_i -
\hat{\mu})^2\$\$ is the weighted residual mean square on \\k - 1\\
degrees of freedom (Röver, Knapp, & Friede, 2015), and the interval is
the pooled estimate plus and minus the \\1 - \alpha/2\\ quantile of the
*t* distribution on \\k - 1\\ degrees of freedom times that standard
error. The function applies this, the unmodified method, as
`test = "knha"` in metafor does: \\q\\ is not replaced by \\\max(1,
q)\\. With few studies the coverage of the adjusted interval is closer
to the nominal level than that of the normal interval. The prediction
interval uses the standard error of the pooled estimate that the
confidence interval uses, so under the adjustment it narrows and widens
with \\q\\ as well.

**When q is small.** Because \\q\\ is not replaced by \\\max(1, q)\\,
the adjusted standard error is smaller than the conventional one
whenever \\q\\ is below 1, which is common (when the effects are
homogeneous, in most samples). That alone does not make the interval
narrower than the normal interval, since the adjusted interval uses a
*t* quantile on \\k - 1\\ degrees of freedom: it is narrower exactly
when \\q\\ is below the squared ratio of the normal quantile to the *t*
quantile, which for a 95 percent interval is 0.02, 0.21, 0.50, and 0.75
at \\k\\ = 2, 3, 5, and 10. A small \\q\\ is part of the reference
distribution of the method and not a fault in the data: when the weights
are the true ones, the pooled estimate minus the mean effect, divided by
the adjusted standard error, follows the *t* distribution exactly, small
values of \\q\\ included. Röver, Knapp, and Friede (2015) describe the
modification \\\max(1, q)\\, which never gives a standard error below
the conventional one; it is not applied here. In the simulations run for
this function the unmodified interval was at or near the nominal level
when the effects were homogeneous or the studies equally precise, and
below it (0.86 to 0.89 for a nominal 0.95, with the normal interval at
0.70 to 0.87) when one study was far more precise than the rest and the
effects were heterogeneous. When the adjusted interval is much narrower
than the conventional one, report both; `hartung_knapp = FALSE` gives
the second. The attributes `hartung_knapp_q` and `se_conventional` hold
\\q\\ and the conventional standard error for that comparison. No
warning is given for a \\q\\ between 0 and 1, which is the method
working as defined.

**Identical effect sizes.** When every effect size is the same, exactly
or up to rounding error, \\q\\ is 0 and so is the adjusted standard
error. The confidence interval and the prediction interval then have no
width, and the *t* statistic is infinite, or undefined when the common
effect size is zero. That is the arithmetic of the method at its
boundary and not evidence of an effect known exactly: each effect size
still has its sampling variance, and identical values arise from coarse
rounding, from one study entered more than once, or by chance with two
or three studies. The function reports the numbers as computed, as
`test = "knha"` in metafor does, and warns with a condition of class
`dmar_meta_hk_zero_se` that states the conventional standard error and
interval. The condition carries the number of studies as `n_studies`,
the pooled estimate as `estimate`, the conventional standard error as
`se_conventional`, and the limits of the normal interval that standard
error gives as `lower_conventional` and `upper_conventional`. For four
effect sizes of 0.5 with sampling variances of 0.04, the adjusted limits
are 0.5 and 0.5; the conventional standard error is 0.1 and the normal
interval 0.304 to 0.696. The check is made on the effect sizes
themselves: effect sizes that differ in any digit beyond rounding error
(0.5 against \\0.5 + 10^{-9}\\) are not identical, give a small positive
\\q\\, and raise no warning.

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

Raudenbush, S. W. (1984). Magnitude of teacher expectancy effects on
pupil IQ as a function of the credibility of expectancy induction: A
synthesis of findings from 18 experiments. *Journal of Educational
Psychology, 76*(1), 85–97.

Röver, C., Knapp, G., & Friede, T. (2015). Hartung-Knapp-Sidik-Jonkman
approach and its modification for random-effects meta-analysis with few
studies. *BMC Medical Research Methodology, 15*, Article 99.
[doi:10.1186/s12874-015-0091-1](https://doi.org/10.1186/s12874-015-0091-1)

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
#>  estimate         0.055 
#>  se               0.036 
#>  test_statistic   1.544 
#>  p_value          0.1401
#>  ci_lower         -0.020
#>  ci_upper         0.130 
#>  prediction_lower -0.020
#>  prediction_upper 0.130 
#>  tau2             0.000 
#>  tau2_ci_lower    0.000 
#>  tau2_ci_upper    0.051 
#>  tau              0.000 
#>  I2               0.000 
#>  I2_ci_lower      0.000 
#>  I2_ci_upper      66.144
#>  H2               1.000 
#>  Q                17.067
#>  Q_df             18    
#>  Q_p              0.5185
#>  n_studies        19    
#> 
#> Confidence level: 95%
#> 
#> test_statistic: referred to Student's t on 18 df (Hartung-Knapp adjustment)

# A fixed effect (common effect) analysis of the same studies.
meta_es(d, v, method = "fe")
#>  term             value 
#>  estimate         0.055 
#>  se               0.037 
#>  test_statistic   1.503 
#>  p_value          0.1328
#>  ci_lower         -0.017
#>  ci_upper         0.127 
#>  prediction_lower <NA>  
#>  prediction_upper <NA>  
#>  tau2             0.000 
#>  tau2_ci_lower    <NA>  
#>  tau2_ci_upper    <NA>  
#>  tau              0.000 
#>  I2               0.000 
#>  I2_ci_lower      <NA>  
#>  I2_ci_upper      <NA>  
#>  H2               1.000 
#>  Q                17.067
#>  Q_df             18    
#>  Q_p              0.5185
#>  n_studies        19    
#> 
#> Confidence level: 95%
#> 
#> test_statistic: referred to the standard normal

# The adjusted standard error beside the conventional one. Five effect
# sizes that agree closely give a weighted residual mean square q far
# below 1, so the adjusted standard error is a small fraction of the
# conventional one; both intervals are worth reporting.
close <- meta_es(c(.30, .32, .29, .31, .33), c(.04, .05, .03, .06, .04))
attr(close, "hartung_knapp_q")
#> [1] 0.006860532
attr(close, "se_conventional")
#> [1] 0.09128709
meta_es(c(.30, .32, .29, .31, .33), c(.04, .05, .03, .06, .04),
        hartung_knapp = FALSE)
#>  term             value 
#>  estimate         0.308 
#>  se               0.091 
#>  test_statistic   3.376 
#>  p_value          0.0007
#>  ci_lower         0.129 
#>  ci_upper         0.487 
#>  prediction_lower 0.018 
#>  prediction_upper 0.599 
#>  tau2             0.000 
#>  tau2_ci_lower    0.000 
#>  tau2_ci_upper    0.000 
#>  tau              0.000 
#>  I2               0.000 
#>  I2_ci_lower      0.000 
#>  I2_ci_upper      0.000 
#>  H2               1.000 
#>  Q                0.027 
#>  Q_df             4     
#>  Q_p              0.9999
#>  n_studies        5     
#> 
#> Confidence level: 95%
#> 
#> test_statistic: referred to the standard normal
```
