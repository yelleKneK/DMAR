# Power of Fisher's Exact Test (Noncentral Hypergeometric)

Computes the power of Fisher's exact test (Fisher, 1934) for the 2
\\\times\\ 2 table under Fisher's noncentral hypergeometric
distribution, where the alternative is parameterized by the true odds
ratio \\\psi\\. The power is the probability that the (conditional)
exact test rejects \\H_0: \psi = 1\\ when in fact \\\psi = \psi_1 \ne
1\\.

## Usage

``` r
power_fisher_exact(
  n_1,
  n_2,
  p_1,
  p_2,
  alpha_level = 0.05,
  alternative = c("two_sided", "less", "greater")
)
```

## Arguments

- n_1, n_2:

  Group sample sizes for the two columns of the 2 \\\times\\ 2 table
  (e.g., treatment vs control).

- p_1, p_2:

  Success probabilities in the two groups under the alternative. The
  odds ratio under the alternative is \\\psi = (p_1 / (1 - p_1)) / (p_2
  / (1 - p_2))\\.

- alpha_level:

  Two-sided significance level. Default `0.05`.

- alternative:

  One of `"two_sided"` (default; the base-R spelling `"two.sided"` is
  accepted as an alias), `"less"`, or `"greater"`.

## Value

A `data.frame` with rows for the power, the alternative-odds-ratio
\\\psi_1\\, the alternative success probabilities \\p_1\\ and \\p_2\\,
the expected column-1 total \\E\[S\]\\ across both groups, and the row /
column totals.

## Details

**Setup.** Fisher's exact test conditions on the marginal totals of the
2 \\\times\\ 2 table:

|         |             |                                 |               |
|---------|-------------|---------------------------------|---------------|
|         | Success     | Failure                         | Total         |
| Group 1 | \\X_1\\     | \\n_1 - X_1\\                   | \\n_1\\       |
| Group 2 | \\S - X_1\\ | \\(n_1 + n_2) - n_1 - S + X_1\\ | \\n_2\\       |
| Total   | \\S\\       | \\n_1 + n_2 - S\\               | \\n_1 + n_2\\ |

Under \\H_0: \psi = 1\\, \\X_1\\ given the marginals follows the
*central* hypergeometric. Under the alternative \\\psi_1\\, \\X_1\\
follows *Fisher's noncentral* hypergeometric with odds ratio parameter
\\\psi_1\\ (Fisher, 1935; Fog, 2008), the conditional distribution of
one binomial count given the total of two independent binomials.
(Wallenius' noncentral hypergeometric, which arises from sequential
biased urn sampling, is a different distribution and is not the relevant
one here.)

**Power calculation.** For each possible value of the column-1 total \\S
= 0, 1, \ldots, n_1 + n_2\\:

1.  Determine the rejection region under \\H_0: \psi = 1\\ using the
    central hypergeometric.

2.  Compute \\\Pr(X_1 \in \mathrm{reject} \mid \psi = \psi_1, S)\\ under
    the noncentral hypergeometric.

3.  Weight by \\\Pr(S \mid \psi_1)\\, the marginal probability of total
    column-1 successes under the alternative.

The power is the resulting weighted sum.

## References

Fisher, R. A. (1934). *Statistical methods for research workers* (5th
ed.). Oliver & Boyd.

Fisher, R. A. (1935). The logic of inductive inference. *Journal of the
Royal Statistical Society, 98*(1), 39–82.

Fog, A. (2008). Sampling methods for Wallenius' and Fisher's noncentral
hypergeometric distributions. *Communications in Statistics – Simulation
and Computation, 37*(2), 241–257.
[doi:10.1080/03610910701790236](https://doi.org/10.1080/03610910701790236)

Good, P. I. (2000). *Permutation tests: A practical guide to resampling
methods for testing hypotheses* (2nd ed.). Springer.

O'Brien, R. G. (1998). A tour of UnifyPow: A SAS module/macro for sample
size analysis. *Proceedings of the 23rd SAS Users Group International
Conference*, 1346–1355.

## See also

[`fisher.test`](https://rdrr.io/r/stats/fisher.test.html)

Other sample size for power:
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_power_R2()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
[`ss_power_R2_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md),
[`ss_power_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
[`ss_power_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_c_ancova.md),
[`ss_power_composite_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md),
[`ss_power_composite_ancova_2group()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md),
[`ss_power_composite_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_anova.md),
[`ss_power_composite_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova.md),
[`ss_power_composite_factorial_ancova_het()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova_het.md),
[`ss_power_composite_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_anova.md),
[`ss_power_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md),
[`ss_power_contrast()`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`ss_power_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md),
[`ss_power_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
[`ss_power_r()`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
[`ss_power_rc()`](https://yelleknek.github.io/DMAR/reference/ss_power_rc.md),
[`ss_power_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`ss_power_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md),
[`ss_power_rm_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_sc()`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Power for n_1 = n_2 = 30 when the two population proportions are
#    0.6 and 0.3. The value comes from enumerating the conditional
#    reference set the test itself uses, not from a normal
#    approximation, so it is the power of the test as conducted.
power_fisher_exact(n_1 = 30, n_2 = 30, p_1 = 0.6, p_2 = 0.3)
#>  term           value
#>  power          0.561
#>  odds_ratio_alt 3.5  
#>  p_1            0.6  
#>  p_2            0.3  
#>  expected_s     27   
#>  n_1            30   
#>  n_2            30   
#>  alpha_level    0.05 

# 2. A difference of 0.10 between proportions is much harder to detect:
#    100 per group is not close to enough. Sample size requirements grow
#    quickly as the difference between the two proportions shrinks.
power_fisher_exact(n_1 = 100, n_2 = 100, p_1 = 0.45, p_2 = 0.35)
#>  term           value
#>  power          0.258
#>  odds_ratio_alt 1.52 
#>  p_1            0.45 
#>  p_2            0.35 
#>  expected_s     80   
#>  n_1            100  
#>  n_2            100  
#>  alpha_level    0.05 
```
