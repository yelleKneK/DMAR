# O'Brien's Test for Homogeneity of Variance

Tests the null hypothesis that two or more groups have equal population
variances using O'Brien's (1981) procedure: each observation is
transformed into a quantity whose expected value equals the group's
variance, and a one-way analysis of variance is then run on those
transformed values. The test is generally regarded as more robust to
non-normality than Bartlett's test while retaining good power.

## Usage

``` r
obrien_test(x, group = NULL, data = NULL, na_action = stats::na.omit)
```

## Arguments

- x:

  Either a numeric vector of observations (in which case `group` must
  also be supplied), or a one-sided formula of the form `y ~ group`, in
  which case `data` is consulted for the variables.

- group:

  A grouping vector or factor of the same length as `x`; used only when
  `x` is a numeric vector.

- data:

  An optional `data.frame` containing the variables named in the
  formula.

- na_action:

  Function specifying how missing values are handled (default
  [`na.omit`](https://rdrr.io/r/stats/na.fail.html)).

## Value

A one-row `data.frame` with columns `statistic` (the *F*-value from the
ANOVA on the transformed scores), `df_1`, `df_2`, `p_value`, `n_groups`,
`n_total`, and `method`.

## Details

Following O'Brien (1981) and the version given in Abdi (2007), each
observation \\Y\_{ij}\\ (the \\j\\th observation in group \\i\\, with
size \\n_i\\ and sample variance \\s_i^2\\) is transformed to
\$\$r\_{ij} = \frac{(n_i - 1.5)\\ n_i\\ (Y\_{ij} - \bar{Y}\_i)^2 - 0.5\\
s_i^2\\ (n_i - 1)}{(n_i - 1)(n_i - 2)}.\$\$ The mean of the \\r\_{ij}\\
within group \\i\\ equals \\s_i^2\\, so a one-way ANOVA on the
\\r\_{ij}\\ tests \\H_0\\: \sigma_1^2 = \cdots = \sigma_k^2\\. Each
group must have at least three observations for the transformation to be
defined.

## References

Abdi, H. (2007). O'Brien's test for homogeneity of variance. In N. J.
Salkind (Ed.), *Encyclopedia of measurement and statistics*. Sage.

O'Brien, R. G. (1981). A simple test for variance effects in
experimental designs. *Psychological Bulletin, 89*(3), 570–574.

## See also

[`bartlett.test`](https://rdrr.io/r/stats/bartlett.test.html),
[`var.test`](https://rdrr.io/r/stats/var.test.html)

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
[`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`factorial_anova()`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Hunter's (1964) "one-is-a-bun" peg-word memory experiment, as discussed
# by Abdi (2007). Sixty-four participants were assigned to a control group
# (no mnemonic instruction) or an experimental group (peg-word mnemonic).
# The score is the number of word pairs (out of 10) recalled. Abdi (2007,
# Table 6) reports F = 1.29 (df = 1, 62) for the O'Brien test of equal
# variances, p = .260 as computed here; the experimental group's apparent
# ceiling effect does not produce statistically detectable variance
# heterogeneity.
hunter_1964 <- data.frame(
  group = factor(
    c(rep("Control", 32), rep("Experimental", 32)),
    levels = c("Control", "Experimental")
  ),
  recall = c(
    # Control group (n = 32):
    5, 5, 5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7,
    8, 8, 8,
    9, 9,
    10, 10,
    # Experimental group (n = 32):
    6,
    7, 7,
    8, 8, 8, 8,
    9, 9, 9, 9, 9, 9, 9, 9, 9,
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
  )
)
obrien_test(recall ~ group, data = hunter_1964)
#>  statistic df_1 df_2 p_value n_groups n_total
#>  1.29      1    62   0.2598  2        64     
#>  method                                    
#>  O'Brien's test for homogeneity of variance

# Comparison against Bartlett's test on the same data.
bartlett.test(recall ~ group, data = hunter_1964)
#> 
#>  Bartlett test of homogeneity of variances
#> 
#> data:  recall by group
#> Bartlett's K-squared = 1.6756, df = 1, p-value = 0.1955
#> 

# Vector / grouping-variable interface, on DMAR's depression_bdi data.
# The wait list variance is about twice the SSRI variance, but with ten
# observations per group the test does not reject equal variances.
obrien_test(depression_bdi$bdi_post, depression_bdi$condition)
#>  statistic df_1 df_2 p_value n_groups n_total
#>  1.29      2    27   0.2918  3        30     
#>  method                                    
#>  O'Brien's test for homogeneity of variance
```
