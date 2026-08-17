# Robust Standardized Mean Difference (Algina-Keselman-Penfield)

Computes the Algina, Keselman, and Penfield (2005) robust standardized
mean difference, which replaces the sample means and pooled SD in
Cohen's *d* with their trimmed-mean and Winsorized-SD counterparts:
\$\$d\_{R} \\=\\ 0.642 \cdot \frac{\bar X\_{t,\\ 1} - \bar X\_{t,\\ 2}}
{s\_{W,\\ p}},\$\$ where \\\bar X\_{t,\\ j}\\ is the trimmed mean of
group \\j\\, \\s\_{W,\\ p}\\ is the pooled Winsorized standard
deviation, and \\0.642\\ is the Algina-Keselman-Penfield (2005) constant
chosen so that \\d_R\\ equals Cohen's \\\delta\\ when the data are
normal. Returns the point estimate, a noncentral *t* confidence
interval, and the trimmed / Winsorized summary statistics.

## Usage

``` r
smd_trimmed(x, y, trim = 0.2, conf_level = 0.95)
```

## Arguments

- x, y:

  Numeric vectors of observations from the two groups.

- trim:

  Proportion to trim from each tail (and Winsorize from each tail). Must
  be in \\\[0, 0.5)\\. Default `0.20` (Wilcox's 2017 recommended
  setting).

- conf_level:

  Confidence level for the CI. Default `0.95`.

## Value

A `data.frame` with rows for the robust *d* estimate, the lower/upper CI
bounds, the per-group trimmed means, the per-group Winsorized SDs, the
pooled Winsorized SD, and the effective sample sizes (after trimming).

## Details

**Why robust.** Under heavy-tailed or skewed marginal distributions, the
conventional Cohen's *d* has very large standard error and biased
coverage. Kelley (2005) documents the coverage distortion of parametric
confidence intervals for the standardized mean difference under
nonnormal distributions. Replacing means by 20%- trimmed means and SD by
20%-Winsorized SD yields an estimator whose efficiency under normality
is roughly 96% (Wilcox, 2017, ch. 5) and whose efficiency under
heavy-tailed contamination is substantially higher than Cohen's *d*.

**The 0.642 constant.** \\0.642 = \mathrm{SD}(X_W) / \mathrm{SD}(X) =
\sqrt{\mathrm{Var}(X_W) / \mathrm{Var}(X)}\\ when \\X \sim N(0, 1)\\ and
\\X_W\\ is the 20%-Winsorized version. Choosing this constant makes
\\d_R = \delta\\ when the data are normal, so the new estimator is on
the same scale as Cohen's *d*.

**CI.** The CI follows the construction of Keselman, Algina, Lix,
Wilcox, and Deering (2008): Yuen's (1974) *t*-statistic on the
trimmed-mean difference (their Equation 8) is referred to a noncentral
*t* distribution with the Yuen-Welch approximate degrees of freedom
(their Equation 9), the noncentrality parameters whose tail
probabilities bracket the observed statistic are located with
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md),
and those limits are rescaled to the \\d_R\\ metric. The degrees of
freedom are reported in the `df_yuen` row of the returned table. At
`trim = 0` the construction reduces to the Welch approximate degrees of
freedom interval; for the exact equal-variance interval on the untrimmed
standardized mean difference use
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md).

## References

Algina, J., Keselman, H. J., & Penfield, R. D. (2005). An alternative to
Cohen's standardized mean difference effect size: A robust parameter and
confidence interval in the two independent groups case. *Psychological
Methods, 10*(3), 317–328.
[doi:10.1037/1082-989X.10.3.317](https://doi.org/10.1037/1082-989X.10.3.317)

Kelley, K. (2005). The effects of nonnormal distributions on confidence
intervals around the standardized mean difference: Bootstrap and
parametric confidence intervals. *Educational and Psychological
Measurement, 65*(1), 51–69.
[doi:10.1177/0013164404264850](https://doi.org/10.1177/0013164404264850)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Keselman, H. J., Algina, J., Lix, L. M., Wilcox, R. R., & Deering, K. N.
(2008). A generally robust approach for testing hypotheses and setting
confidence intervals for effect sizes. *Psychological Methods, 13*(2),
110–129.
[doi:10.1037/1082-989X.13.2.110](https://doi.org/10.1037/1082-989X.13.2.110)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons and Chapter 3
on one-way ANOVA.)

Wilcox, R. R. (2017). *Introduction to robust estimation and hypothesis
testing* (4th ed.). Academic Press.

Yuen, K. K. (1974). The two-sample trimmed *t* for unequal population
variances. *Biometrika, 61*(1), 165–170.

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`var_smd_trimmed`](https://yelleknek.github.io/DMAR/reference/var_smd_trimmed.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
[`eta_squared()`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
[`eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md),
[`expected_partial_r()`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md),
[`expected_r()`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
[`expected_smd()`](https://yelleknek.github.io/DMAR/reference/expected_smd.md),
[`nnt_from_smd()`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`omega_squared_partial()`](https://yelleknek.github.io/DMAR/reference/omega_squared_partial.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Two normal groups: robust d agrees closely with Cohen's d.
set.seed(113)
x <- rnorm(40, 0,   1); y <- rnorm(40, 0.5, 1)
smd_trimmed(x, y)
#>  term                 value 
#>  smd_trimmed          -0.405
#>  lower_limit          -0.883
#>  upper_limit          0.0775
#>  trimmed_mean_x       0.187 
#>  trimmed_mean_y       0.589 
#>  winsorized_sd_x      0.636 
#>  winsorized_sd_y      0.641 
#>  winsorized_sd_pooled 0.638 
#>  h_x                  24    
#>  h_y                  24    
#>  trim                 0.2   
#>  df_yuen              46    
#> 
#> Confidence level: 95%

# 2. Contaminated y: a few outliers; robust d shifts much less
#        than Cohen's d.
set.seed(113)
x <- rnorm(40, 0, 1)
y <- c(rnorm(38, 0.5, 1), 30, -25)
smd_trimmed(x, y)
#>  term                 value 
#>  smd_trimmed          -0.349
#>  lower_limit          -0.826
#>  upper_limit          0.131 
#>  trimmed_mean_x       0.187 
#>  trimmed_mean_y       0.539 
#>  winsorized_sd_x      0.636 
#>  winsorized_sd_y      0.658 
#>  winsorized_sd_pooled 0.647 
#>  h_x                  24    
#>  h_y                  24    
#>  trim                 0.2   
#>  df_yuen              45.9  
#> 
#> Confidence level: 95%
smd(x, y)
#>  term value 
#>  smd  -0.105
```
