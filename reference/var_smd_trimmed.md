# Asymptotic Variance of the Robust Trimmed SMD

Computes the asymptotic variance of the Algina-Keselman-Penfield (2005)
robust standardized mean difference under the Yuen (1974) trimmed-mean
framework, suitable for AIPE sample size planning for robust effect
sizes (Keselman, Algina, Lix, Wilcox, & Deering, 2008).

## Usage

``` r
var_smd_trimmed(population_smd_trimmed, n_1, n_2, trim = 0.2)
```

## Arguments

- population_smd_trimmed:

  Anticipated population value of the robust trimmed SMD \\\delta_R\\.

- n_1, n_2:

  Per-group sample sizes.

- trim:

  Proportion to trim and Winsorize. Default `0.20`.

## Value

A 1-row `data.frame` with columns `term` (`"var_smd_trimmed"`) and
`value` (the variance).

## Details

**Variance formula.** Under random sampling with trimming proportion
\\\gamma\\ from each tail, the variance of the trimmed-mean difference
scales by \\1 / h_j\\ (where \\h_j = n_j - 2 \lfloor \gamma n_j
\rfloor\\ is the number of retained observations in group \\j\\) rather
than \\1 / n_j\\. The large-sample variance of the standardized version,
written on the \\d_R\\ scale, is \$\$\mathrm{Var}(\hat d_R) \\\approx\\
\frac{h_1 + h_2}{h_1 h_2} + \frac{\delta_R^2}{2 (h_1 + h_2)}.\$\$ For
\\\gamma = 0\\ this reduces to the standard Hedges-Olkin (1985) variance
of Cohen's *d*.

**When to use.** For AIPE planning of a robust effect size study, use
`var_smd_trimmed()` in place of
[`var_smd()`](https://yelleknek.github.io/DMAR/reference/var_smd.md).
Pair with
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)
for the point estimate and noncentral *t* CI.

## References

Algina, J., Keselman, H. J., & Penfield, R. D. (2005). An alternative to
Cohen's standardized mean difference effect size: A robust parameter and
confidence interval in the two independent groups case. *Psychological
Methods, 10*(3), 317–328.
[doi:10.1037/1082-989X.10.3.317](https://doi.org/10.1037/1082-989X.10.3.317)

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

Yuen, K. K. (1974). The two-sample trimmed *t* for unequal population
variances. *Biometrika, 61*(1), 165–170.

## See also

[`smd_trimmed`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md),
[`var_smd`](https://yelleknek.github.io/DMAR/reference/var_smd.md),
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)

Other variance utilities:
[`var_alpha()`](https://yelleknek.github.io/DMAR/reference/var_alpha.md),
[`var_cv()`](https://yelleknek.github.io/DMAR/reference/var_cv.md),
[`var_ete()`](https://yelleknek.github.io/DMAR/reference/var_ete.md),
[`var_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md),
[`var_omega_squared()`](https://yelleknek.github.io/DMAR/reference/var_omega_squared.md),
[`var_r()`](https://yelleknek.github.io/DMAR/reference/var_r.md),
[`var_smd()`](https://yelleknek.github.io/DMAR/reference/var_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Population delta_R = 0.5, n = 30 per group, 20% trim:
var_smd_trimmed(population_smd_trimmed = 0.5, n_1 = 30, n_2 = 30)
#>  term            value
#>  var_smd_trimmed 0.115

# 2. Variance scales by h / n via the trimming proportion:
var_smd_trimmed(0.5, 30, 30, trim = 0.00)$value
#> [1] 0.06875
var_smd_trimmed(0.5, 30, 30, trim = 0.20)$value
#> [1] 0.1145833
```
