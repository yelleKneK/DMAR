# TOST (Two One-Sided Tests) and Noninferiority for a Linear Contrast

Performs the two one-sided tests procedure (Schuirmann, 1987) for
equivalence, and the companion one-sided noninferiority test, for a
linear contrast of group means \\\psi = \sum_j c_j \mu_j\\ in a fixed
effects design with one pooled error term. The equivalence null
hypothesis is that the population contrast lies outside the
user-specified bounds, \\H_0: \psi \le -\delta_L \cup \psi \ge
\delta_U\\, against the alternative \\H_1: -\delta_L \< \psi \<
\delta_U\\. The noninferiority null is \\H_0: \psi \le -\delta_L\\
against \\H_1: \psi \> -\delta_L\\. Following Chattopadhyay,
Bandyopadhyay, Kelley, and Padalunkal (2025), the equivalence decision
is read from the 100(1 - 2\\\alpha\\)% confidence interval: equivalence
is declared when the whole interval lies inside \\(-\delta_L,
\delta_U)\\, and noninferiority when the interval's lower limit exceeds
\\-\delta_L\\.

## Usage

``` r
tost_c(
  means = NULL,
  s_anova = NULL,
  c_weights = NULL,
  n = NULL,
  psi_hat = NULL,
  se = NULL,
  df_error = NULL,
  delta_lower = NULL,
  delta_upper = NULL,
  benchmark = NULL,
  alpha_level = 0.05
)
```

## Arguments

- means:

  A vector of group means. Supply together with `s_anova`, `c_weights`,
  and `n`. Alternatively, use the direct interface via `psi_hat`, `se`,
  and `df_error`.

- s_anova:

  The standard deviation of the errors from the ANOVA model (the square
  root of the mean square error), as in
  [`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md).

- c_weights:

  The contrast weights. For a mean comparison the weights must sum to
  zero, and, so that the bounds are on the raw scale of the response,
  the positive weights must sum to 1 and the negative weights to -1 (use
  fractional values, not integers). When `benchmark` is supplied, the
  weights must instead be nonnegative and sum to 1 (typically a single 1
  selecting one group).

- n:

  Sample sizes per group (if length 1, equal group sizes are assumed).

- psi_hat:

  The estimated contrast, for the direct interface. Supply together with
  `se` and `df_error` when the contrast and its standard error have
  already been computed (for example, from a fitted model with
  covariates).

- se:

  The standard error of `psi_hat`, for the direct interface.

- df_error:

  The error degrees of freedom. On the summary-statistic interface the
  default is \\N - J\\, with \\J\\ the number of groups; it must be
  supplied for designs with additional factors. Required on the direct
  interface.

- delta_lower, delta_upper:

  Equivalence bounds on the raw scale of the response. Both must be
  positive; the equivalence region is \\(-\delta_L, +\delta_U)\\. If
  only `delta_upper` is supplied, the bounds are symmetric.
  Noninferiority uses \\-\delta_L\\ alone.

- benchmark:

  An optional known constant to compare against (for example, a
  normative or regulatory cutoff). When supplied, the contrast is
  \\\sum_j c_j \mu_j - b\\ with nonnegative weights summing to 1, and
  the constant contributes no sampling variability.

- alpha_level:

  One-sided significance level for each of the two tests. Default
  `0.05`, so the interval the decisions are read from is the 90% CI.

## Value

A `data.frame` with rows for the estimated contrast (`psi_hat`), its
standard error (`se`), the error degrees of freedom (`df`), the two
one-sided test statistics (`t_lower`, `t_upper`) and their *p*-values
(`p_lower`, `p_upper`), the joint TOST *p*-value (`p_tost`, the larger
of the two), the noninferiority *p*-value (`p_noninferiority`, equal to
`p_lower` by construction), the 100(1 - 2\\\alpha\\)% confidence limits
(`lower_limit`, `upper_limit`), the bounds (`delta_lower`, stored as the
signed lower bound, and `delta_upper`), and four binary decision flags
(`equivalent`, `noninferior`, `superior`, `inferior`; 1 = declared, 0 =
not). When all four flags are 0, the interval straddles a bound and the
result is inconclusive. The five-way classification is also attached as
the `"verdict"` attribute, one of `"Equivalent"`, `"Superior"`,
`"Inferior"`, `"Non-inferior only"`, or `"Inconclusive"`.

## Details

**One pooled error term.** On the summary-statistic interface the
standard error is \\\mathrm{SE}(\hat\psi) = s\_{\mathrm{anova}}
\sqrt{\sum_j c_j^2 / n_j}\\, the model comparison position of Maxwell,
Delaney, and Kelley (2027): every one-degree-of-freedom contrast is
judged against the same yardstick, the root mean square error of one
model fit to all groups.

**The verdict logic.** Reading the 100(1 - 2\\\alpha\\)% CI against the
bounds: an interval entirely inside \\(-\delta_L, \delta_U)\\ is
*equivalent*; entirely above \\\delta_U\\ is *superior* (which implies
noninferior); entirely below \\-\delta_L\\ is *inferior*; a lower limit
above \\-\delta_L\\ with an upper limit past \\\delta_U\\ is
*noninferior only*; and an interval straddling a bound is
*inconclusive*. An inconclusive result is a statement about precision,
not evidence of a difference: only an interval clearing a bound entirely
licenses a directional claim.

**Why the weights must sum to \\\pm 1\\.** A bound stated in raw units
of the response is only meaningful if the contrast is itself a simple
difference of (weighted) means on that scale, which requires the
positive weights to sum to 1 and the negative weights to -1. A weight
vector such as `c(2, -2)` would silently double the effective bounds, so
it is rejected rather than rescaled.

**Choosing the bounds.** The bounds must be fixed before the data are
examined, on substantive grounds: the smallest difference that would
matter (Serlin & Lapsley, 1985; Lakens, Scheel, & Isager, 2018). They
are never derived from a standard error, which would make the definition
of "close enough" a function of the sample size.

**Agreement with emmeans.** The *p*-values reproduce
`emmeans::test(..., side = "equivalence")` and
`emmeans::test(..., side = "noninferiority")` with `adjust = "none"` to
machine precision. Note two emmeans pitfalls the interface here avoids:
`side = "left"` tests non-superiority, not noninferiority, and
`trt.vs.ctrl` families silently apply a Dunnett-type adjustment unless
`adjust = "none"` is passed.

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, P. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods, 30*(2),
425–439. [doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence testing
for psychological research: A tutorial. *Advances in Methods and
Practices in Psychological Science, 1*(2), 259–269.
[doi:10.1177/2515245918770963](https://doi.org/10.1177/2515245918770963)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons of means.)

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*(6), 657–680.

Serlin, R. C., & Lapsley, D. K. (1985). Rationality in psychological
research: The good-enough principle. *American Psychologist, 40*(1),
73–83.

Wellek, S. (2010). *Testing statistical hypotheses of equivalence and
noninferiority* (2nd ed.). Chapman & Hall/CRC.

## See also

[`tost_smd`](https://yelleknek.github.io/DMAR/reference/tost_smd.md),
[`tost_r`](https://yelleknek.github.io/DMAR/reference/tost_r.md),
[`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`power_equivalence_c`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`ss_power_equivalence_c`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`plot_equivalence`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md)

Other equivalence testing:
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`power_density_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md),
[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`tost_r()`](https://yelleknek.github.io/DMAR/reference/tost_r.md),
[`tost_smd()`](https://yelleknek.github.io/DMAR/reference/tost_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Two of five groups compared against the reference group, with a
#    pooled error term from one model across all five groups.
#    Bounds of 5 raw-scale points; alpha_level = .05, so decisions read
#    from the 90% CI.
tost_c(means = c(70.40, 55.61, 51.91, 65.66, 65.12),
       s_anova = 15.67,
       c_weights = c(-1, 0, 0, 0, 1),
       n = c(113, 74, 76, 80, 61),
       delta_upper = 5)
#>  term             value   
#>  psi_hat          -5.28   
#>  se               2.49    
#>  df               399     
#>  t_lower          -0.112  
#>  t_upper          -4.13   
#>  p_lower          0.5447  
#>  p_upper          < 0.0001
#>  p_tost           0.5447  
#>  p_noninferiority 0.5447  
#>  lower_limit      -9.38   
#>  upper_limit      -1.18   
#>  delta_lower      -5      
#>  delta_upper      5       
#>  equivalent       0       
#>  noninferior      0       
#>  superior         0       
#>  inferior         0       
#> 
#> Confidence level: 90%

# 2. The same contrast through the direct interface, as when the
#    estimate and standard error come from a model with covariates.
res <- tost_c(psi_hat = -5.28, se = 2.49, df_error = 399,
              delta_upper = 5)
res
#>  term             value   
#>  psi_hat          -5.28   
#>  se               2.49    
#>  df               399     
#>  t_lower          -0.112  
#>  t_upper          -4.13   
#>  p_lower          0.5447  
#>  p_upper          < 0.0001
#>  p_tost           0.5447  
#>  p_noninferiority 0.5447  
#>  lower_limit      -9.39   
#>  upper_limit      -1.17   
#>  delta_lower      -5      
#>  delta_upper      5       
#>  equivalent       0       
#>  noninferior      0       
#>  superior         0       
#>  inferior         0       
#> 
#> Confidence level: 90%
attr(res, "verdict")
#> [1] "Inconclusive"

# 3. A group mean against a fixed benchmark of 68: the constant
#    contributes no sampling variability.
tost_c(means = c(70.40, 55.61, 51.91, 65.66, 65.12),
       s_anova = 15.67,
       c_weights = c(1, 0, 0, 0, 0),
       n = c(113, 74, 76, 80, 61),
       benchmark = 68, delta_upper = 5)
#>  term             value   
#>  psi_hat          2.4     
#>  se               1.47    
#>  df               399     
#>  t_lower          5.02    
#>  t_upper          -1.76   
#>  p_lower          < 0.0001
#>  p_upper          0.0393  
#>  p_tost           0.0393  
#>  p_noninferiority < 0.0001
#>  lower_limit      -0.0303 
#>  upper_limit      4.83    
#>  delta_lower      -5      
#>  delta_upper      5       
#>  equivalent       1       
#>  noninferior      1       
#>  superior         0       
#>  inferior         0       
#> 
#> Confidence level: 90%

# 4. Asymmetric bounds: a shortfall of 3 matters, an excess of 8 does.
tost_c(psi_hat = 1.2, se = 1.1, df_error = 120,
       delta_lower = 3, delta_upper = 8)
#>  term             value   
#>  psi_hat          1.2     
#>  se               1.1     
#>  df               120     
#>  t_lower          3.82    
#>  t_upper          -6.18   
#>  p_lower          0.0001  
#>  p_upper          < 0.0001
#>  p_tost           0.0001  
#>  p_noninferiority 0.0001  
#>  lower_limit      -0.623  
#>  upper_limit      3.02    
#>  delta_lower      -3      
#>  delta_upper      8       
#>  equivalent       1       
#>  noninferior      1       
#>  superior         0       
#>  inferior         0       
#> 
#> Confidence level: 90%
```
