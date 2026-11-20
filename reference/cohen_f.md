# Cohen's *f* Effect Size

Computes Cohen's *f*, the standard deviation of the group means divided
by the standard deviation within the groups, from one of three complete
specifications: the group means with the within-group variance (and the
group sizes, when the groups differ in size), eta squared, or the two
standard deviations.

## Usage

``` r
cohen_f(
  means = NULL,
  var_within = NULL,
  n_per_group = NULL,
  eta_squared = NULL,
  sd_means = NULL,
  sd_within = NULL
)
```

## Arguments

- means:

  Numeric vector of group means, one per group and at least two:
  population means for the population *f*, or observed sample means for
  the plug-in value. Give with `var_within`.

- var_within:

  The within-group variance, a single positive number: the population
  variance common to the groups, or, for the plug-in value, the pooled
  within-group mean square (the mean square error of the one-way
  analysis of variance). Give with `means`.

- n_per_group:

  Optional, and used only with `means`: the group sizes, one per group,
  which weight each group mean by its share of the total, \\n_j / N\\. A
  single number states equal groups and gives the same result as `NULL`,
  the default, which weights the means equally (\\1 / k\\ each, for *k*
  groups).

- eta_squared:

  Eta squared, the proportion of the variance that the groups account
  for, a single number in \[0, 1). Give it alone. A population value
  gives the population *f*. A sample eta squared gives the plug-in value
  \\\sqrt{df_1 F / df_2}\\, which is larger than the value from the
  observed means and the pooled mean square (see Details).

- sd_means:

  The standard deviation of the group means, \\\sigma_m\\, a single
  number of at least 0, computed with divisor *k* (or with the weights
  \\n_j / N\\), not with the divisor *k* - 1 that
  [`stats::sd()`](https://rdrr.io/r/stats/sd.html) uses. Give with
  `sd_within`.

- sd_within:

  The within-group standard deviation, \\\sigma\\, a single positive
  number, or, for the plug-in value, the square root of the pooled
  within-group mean square. Give with `sd_means`.

## Value

A one-row `data.frame` (printed through `dmar_tbl`) with columns `term`
and `value`; `term` is `"cohen_f"` and `value` is *f* at full precision.

## Details

Cohen's *f* is defined with population quantities (Cohen, 1988, p. 276):
\\f = \sigma_m / \sigma\\, where \\\sigma_m\\ is the standard deviation
of the population means and \\\sigma\\ is the standard deviation within
the populations. Given population values, the function returns the
population *f*. Given sample summaries, it returns the plug-in value
from the same calculation. The plug-in value depends on which summaries
are given: from the same data, a sample eta squared gives a larger value
than the observed group means with the pooled within-group mean square
(see below). Sample inputs do not make the result an unbiased estimate
of the population *f*: the spread of sample means includes their
sampling error as well as the spread of the population means. Unequal
group sizes determine how the group means are weighted; without
`n_per_group` the means are weighted equally.

A call gives exactly one complete specification: `means` with
`var_within` (and `n_per_group` if wanted), `eta_squared` alone, or
`sd_means` with `sd_within`. Each determines *f* on its own, so an
argument of one specification given with an argument of another is
refused, as is one half of a pair given alone, or `n_per_group` without
`means`; the refusal has class `dmar_argument_invalid`.

From the means, with weights \\w_j = 1 / k\\ by default or \\w_j = n_j /
N\\ when `n_per_group` is given, \$\$\bar\mu = \sum_j w_j \mu_j, \qquad
\sigma_m^2 = \sum_j w_j (\mu_j - \bar\mu)^2, \qquad f = \sqrt{\sigma_m^2
/ \sigma^2}.\$\$ The divisor of \\\sigma_m^2\\ is *k* (or *N*), because
the variance is that of a fixed set of means. From eta squared, \\f =
\sqrt{\eta^2 / (1 - \eta^2)}\\ (Kim, 2016, p. 233), and from the two
standard deviations, \\f = \sigma_m / \sigma\\. Given population values
for the same populations, with the same weighting of the means, the
three give the same *f*.

Given sample summaries, the three agree only when the within-group
variance is defined the same way in each. With the observed means, the
group sizes, and the pooled within-group mean square, \\SS\_{within} /
(N - k)\\, of a one-way design, \\N \sigma_m^2\\ is the between-groups
sum of squares, so the plug-in value equals \\\sqrt{df_1 F / N}\\, with
*F*, \\df_1 = k - 1\\, and \\df_2 = N - k\\ from the analysis of
variance. A sample eta squared, \\SS\_{between} / SS\_{total}\\, pairs
with the within-group variance \\SS\_{within} / N\\ rather than with the
mean square, so the value from it, \\\sqrt{SS\_{between} /
SS\_{within}}\\ \\= \sqrt{df_1 F / df_2}\\, is larger than the value
from the means and the pooled mean square by the factor \\\sqrt{N / (N -
k)}\\. The two are equal when `var_within` is \\SS\_{within} / N\\ (or
`sd_within` is its square root). Without `n_per_group`, the means of
unequal groups are weighted equally and the value differs from both. The
square root of the signal-to-noise ratio of
[`ci_srsnr`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md) is
this *f* with the weights \\n_j / N\\ and the pooled mean square; that
function gives its confidence interval from *F* or from the same
summaries.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Kim, H.-Y. (2016). Statistical notes for clinical researchers: Sample
size calculation 3. Comparison of several means using one-way ANOVA.
*Restorative Dentistry & Endodontics, 41*(3), 231–234.
[doi:10.5395/rde.2016.41.3.231](https://doi.org/10.5395/rde.2016.41.3.231)

## See also

[`ci_srsnr`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`ci_snr`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ss_power_one_way_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Kim (2016, Table 1) plans a four-group study around means of 18.4,
# 22.2, 25.1, and 32.1 with a common standard deviation of 8.6, and
# reports f as .58. At full precision, this f plans a balanced study of
# 40, ten per group, for power .80 at the .05 level, as the article does.
f_kim <- cohen_f(means = c(18.4, 22.2, 25.1, 32.1), var_within = 8.6^2)
f_kim
#>  term    value
#>  cohen_f 0.583
ss_power_one_way_anova(a = 4, true_f = f_kim$value, desired_power = .80)
#>  term              value 
#>  necessary_N       40    
#>  n_per_group       10    
#>  a                 4     
#>  noncentral_f_parm 13.603
#>  actual_power      0.847 

# Cohen (1988, p. 277) places seven population means over half a
# within-population standard deviation, one at each end and five at the
# midpoint, and gives f as .13 (.267 times .5, his Equation 8.2.8).
cohen_f(means = c(-.25, 0, 0, 0, 0, 0, .25), var_within = 1)
#>  term    value
#>  cohen_f 0.134

# The plug-in value from the summaries of an observed one-way design
# with unequal groups: the group means, the pooled within-group mean
# square, and the group sizes of the three arms of drinks_trial.
fit <- stats::lm(log_drinks ~ treatment, data = drinks_trial)
cohen_f(means = tapply(drinks_trial$log_drinks, drinks_trial$treatment, mean),
        var_within = sum(stats::residuals(fit)^2) / stats::df.residual(fit),
        n_per_group = as.vector(table(drinks_trial$treatment)))
#>  term    value
#>  cohen_f 0.311

# The sample eta squared of the same fit gives a larger plug-in value,
# by the factor sqrt(N / (N - k)) = sqrt(88 / 85) (see Details).
ss <- stats::anova(fit)[["Sum Sq"]]
cohen_f(eta_squared = ss[1] / sum(ss))
#>  term    value
#>  cohen_f 0.317

# From eta squared, and from the two standard deviations:
cohen_f(eta_squared = .10)
#>  term    value
#>  cohen_f 0.333
cohen_f(sd_means = 4, sd_within = 8)
#>  term    value
#>  cohen_f 0.500
```
