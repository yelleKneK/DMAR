# Consequences of a Design: Power, Sign and Magnitude Errors, and Expected Precision

Evaluates what a design of a given precision will actually deliver,
under both of the package's lenses at once. The *significance lens*: the
`power` of the two-sided test, the `type_s_error` (the probability that
a statistically significant estimate has the wrong sign), and the
`exaggeration_ratio` (Type M: the average factor by which significant
estimates overstate the true effect), following the design analysis of
Gelman and Carlin (2014). The *precision lens*, in the accuracy in
parameter estimation (AIPE) tradition: the expected half-width and full
width of the `conf_level` confidence interval the design will produce,
the spread of that realized width, and, when a target width `w` is
supplied, `pct_ci_less_w`, the probability that the realized interval is
no wider than the target, the same quantities the
`ss_aipe_*_sensitivity()` family estimates by Monte Carlo, here in
closed form.

## Usage

``` r
design_consequences(
  true_effect = NULL,
  se = NULL,
  sd = NULL,
  n_1 = NULL,
  n_2 = NULL,
  alpha_level = 0.05,
  df = NULL,
  conf_level = 0.95,
  w = NULL
)
```

## Arguments

- true_effect:

  The assumed true (population) effect, on the scale of the estimate (a
  mean difference, a regression coefficient, a standardized mean
  difference). May be negative. May be `NULL`, in which case the
  significance-lens rows are returned as `NA` and only the precision
  lens (which does not involve the true effect) is informative.

- se:

  The standard error of the estimate the design will produce, on the
  same scale as `true_effect`. Supply `se` directly, or supply `sd` with
  `n_1` (and `n_2`) and let the function derive it.

- sd, n_1, n_2:

  An alternative to `se` for the two most common cases: with `sd` and
  `n_1` only, the one-sample (or paired-difference) design
  `se = sd / sqrt(n_1)` with `df = n_1 - 1`; with `n_2` as well, the
  two-group design `se = sd * sqrt(1/n_1 + 1/n_2)` with
  `df = n_1 + n_2 - 2` (`sd` is the common within-group standard
  deviation). A `df` supplied explicitly overrides the derived one.

- alpha_level:

  Two-sided Type I error rate of the significance test. Defaults to
  0.05.

- df:

  Degrees of freedom of the reference *t* distribution. Defaults to
  `Inf` (the normal case) unless derived from `n_1` / `n_2`.

- conf_level:

  Confidence level of the interval evaluated by the precision lens.
  Defaults to 0.95.

- w:

  Optional target full width for the confidence interval; when supplied,
  `pct_ci_less_w` reports the probability that the realized interval is
  no wider than `w`.

## Value

A `data.frame` (class `dmar_tbl`) with the significance-lens rows
(`power`, `type_s_error`, `exaggeration_ratio`), the precision-lens rows
(`expected_half_width`, `mean_ci_width`, `median_ci_width`,
`sd_ci_width`, `pct_ci_less_w`, `target_width`; the last two are `NA`
when no `w` is supplied), and the design rows (`true_effect`, `se`,
`df`, `alpha_level`). The confidence level is recorded in the
`"conf_level"` attribute. The schema is constant: rows that do not apply
are `NA`, never dropped.

## Details

Together they answer the two questions a chosen design should be
interrogated with before data collection: *if I run this study and
filter it through a significance test, what will the published record
look like?* and *how precisely will I estimate the effect regardless of
significance?* An underpowered design fails both: its significant
estimates are exaggerated and possibly sign-reversed, and its confidence
intervals are too wide to be informative.

**Significance lens.** Writing \\\lambda = \theta / \mathrm{se}\\ and
\\c\\ for the two-sided critical value, the power and the Type S error
follow from the two tails of the (noncentral *t* or normal) distribution
of the estimate, and the exaggeration ratio is the expected absolute
estimate conditional on significance over the absolute true effect.
Gelman and Carlin's `retrodesign()` computes the exaggeration ratio by
simulation; here it is computed exactly, from truncated normal moments
when `df = Inf` and by numerical integration against the *t* density
otherwise, so no simulation error enters. When `true_effect = 0` the
power equals `alpha_level`, the Type S error is 0.5, and the
exaggeration ratio is undefined (`NA`).

**Precision lens.** The realized interval half-width is
\\t\_{1-\alpha^\*/2,\\\mathit{df}} \cdot \widehat{\mathrm{se}}\\ with
\\\alpha^\* = 1 - \mathtt{conf\\level}\\, and \\\widehat{\mathrm{se}} =
\mathrm{se}\sqrt{W/\mathit{df}}\\ with \\W \sim \chi^2\_{\mathit{df}}\\,
so the width's mean, median, and standard deviation have closed
chi-distribution forms and \$\$P(\mathrm{width} \le w) \\=\\ P\\\left(W
\le \mathit{df}\left\[\frac{w}
{2\\t\\\mathrm{se}}\right\]^{2}\right).\$\$ These are the population
versions of the `mean_ci_width`, `median_ci_width`, `sd_ci_width`, and
`pct_ci_less_w` terms that the `ss_aipe_*_sensitivity()` functions
estimate by Monte Carlo. With `df = Inf` the standard error is treated
as known, the width is deterministic, and `pct_ci_less_w` is a step: 1
when the fixed width is at most `w` and 0 otherwise.

The function complements the planners rather than replacing them:
`ss_power_*()` chooses a sample size for detection, `ss_aipe_*()`
chooses one for precision, and `design_consequences()` interrogates
whatever design came out (or the design a completed study used). The two
lenses are the power and accuracy in parameter estimation approaches to
sample size planning reviewed by Maxwell, Kelley, and Rausch (2008).

## References

Gelman, A., & Carlin, J. (2014). Beyond power calculations: Assessing
Type S (sign) and Type M (magnitude) errors. *Perspectives on
Psychological Science, 9*(6), 641–651.
[doi:10.1177/1745691614551642](https://doi.org/10.1177/1745691614551642)
(Their accompanying `retrodesign()` function obtains the exaggeration
ratio by simulation; the closed-form and *t*-integration computation
here is exact.)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
obtaining precision: Delineating methods of sample size planning.
*Evaluation and the Health Professions, 26*(3), 258–287.
[doi:10.1177/0163278703255242](https://doi.org/10.1177/0163278703255242)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

## See also

[`ss_power_smd`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md)
and
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)
for choosing the sample size by detection or by precision before this
function interrogates the choice;
[`expected_smd`](https://yelleknek.github.io/DMAR/reference/expected_smd.md)
for the unconditional small-sample bias of the standardized mean
difference, a different bias than the significance-filter exaggeration
here.

Other design utilities:
[`deft()`](https://yelleknek.github.io/DMAR/reference/deft.md),
[`effects_coding()`](https://yelleknek.github.io/DMAR/reference/effects_coding.md),
[`helmert_coding()`](https://yelleknek.github.io/DMAR/reference/helmert_coding.md),
[`is_orthogonal_set()`](https://yelleknek.github.io/DMAR/reference/is_orthogonal_set.md),
[`orthogonal_polynomial()`](https://yelleknek.github.io/DMAR/reference/orthogonal_polynomial.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# ---- Both lenses on one underpowered design --------------------------
# True effect 0.1 measured with standard error 0.3 (say, d = .1 with
# about 45 per group): power is 6 percent, a significant result has a
# 17 percent chance of the wrong sign and overstates the truth
# seven-fold, and the 95 percent CI is about 1.2 wide, twelve times the
# effect. Bad for detection, bad for precision.
design_consequences(true_effect = 0.1, se = 0.3)
#>  term                value 
#>  power               0.0628
#>  type_s_error        0.174 
#>  exaggeration_ratio  7.1   
#>  expected_half_width 0.588 
#>  mean_ci_width       1.18  
#>  median_ci_width     1.18  
#>  sd_ci_width         0     
#>  pct_ci_less_w       <NA>  
#>  target_width        <NA>  
#>  true_effect         0.1   
#>  se                  0.3   
#>  df                  Inf   
#>  alpha_level         0.05  
#> 
#> Confidence level: 95%

# ---- The same effect, precisely measured -----------------------------
design_consequences(true_effect = 0.1, se = 0.03)
#>  term                value   
#>  power               0.915   
#>  type_s_error        6.56e-08
#>  exaggeration_ratio  1.05    
#>  expected_half_width 0.0588  
#>  mean_ci_width       0.118   
#>  median_ci_width     0.118   
#>  sd_ci_width         0       
#>  pct_ci_less_w       <NA>    
#>  target_width        <NA>    
#>  true_effect         0.1     
#>  se                  0.03    
#>  df                  Inf     
#>  alpha_level         0.05    
#> 
#> Confidence level: 95%

# ---- From a planned two-group design (sd and per-group n) ------------
# d-type effect 0.4, common sd 1, 60 per group; finite df flows through
# both lenses.
design_consequences(true_effect = 0.4, sd = 1, n_1 = 60, n_2 = 60)
#>  term                value   
#>  power               0.583   
#>  type_s_error        4.99e-05
#>  exaggeration_ratio  1.31    
#>  expected_half_width 0.361   
#>  mean_ci_width       0.722   
#>  median_ci_width     0.721   
#>  sd_ci_width         0.047   
#>  pct_ci_less_w       <NA>    
#>  target_width        <NA>    
#>  true_effect         0.4     
#>  se                  0.183   
#>  df                  118     
#>  alpha_level         0.05    
#> 
#> Confidence level: 95%

# ---- Will the interval beat a target width? --------------------------
# Same design, asking for the probability the realized 95 percent CI is
# no wider than 0.7 (the closed-form pct_ci_less_w that the
# ss_aipe_smd_sensitivity() simulation estimates).
design_consequences(true_effect = 0.4, sd = 1, n_1 = 60, n_2 = 60,
                    w = 0.7)
#>  term                value   
#>  power               0.583   
#>  type_s_error        4.99e-05
#>  exaggeration_ratio  1.31    
#>  expected_half_width 0.361   
#>  mean_ci_width       0.722   
#>  median_ci_width     0.721   
#>  sd_ci_width         0.047   
#>  pct_ci_less_w       0.326   
#>  target_width        0.7     
#>  true_effect         0.4     
#>  se                  0.183   
#>  df                  118     
#>  alpha_level         0.05    
#> 
#> Confidence level: 95%

# ---- Precision lens alone (no effect assumption needed) --------------
design_consequences(true_effect = NULL, sd = 1, n_1 = 60, n_2 = 60,
                    w = 0.7)
#>  term                value
#>  power               <NA> 
#>  type_s_error        <NA> 
#>  exaggeration_ratio  <NA> 
#>  expected_half_width 0.361
#>  mean_ci_width       0.722
#>  median_ci_width     0.721
#>  sd_ci_width         0.047
#>  pct_ci_less_w       0.326
#>  target_width        0.7  
#>  true_effect         <NA> 
#>  se                  0.183
#>  df                  118  
#>  alpha_level         0.05 
#> 
#> Confidence level: 95%

# ---- After an AIPE plan: check the detection side --------------------
# Plan for a full width of 0.5 on the SMD at delta = .4, then ask what
# that design does under the significance filter.
n_plan <- ss_aipe_smd(delta = 0.4, width = 0.5)$value[1]
design_consequences(true_effect = 0.4, sd = 1,
                    n_1 = n_plan, n_2 = n_plan, w = 0.5)
#>  term                value   
#>  power               0.885   
#>  type_s_error        3.06e-07
#>  exaggeration_ratio  1.07    
#>  expected_half_width 0.248   
#>  mean_ci_width       0.496   
#>  median_ci_width     0.496   
#>  sd_ci_width         0.0222  
#>  pct_ci_less_w       0.578   
#>  target_width        0.5     
#>  true_effect         0.4     
#>  se                  0.126   
#>  df                  250     
#>  alpha_level         0.05    
#> 
#> Confidence level: 95%
```
