# Sample Size Planning for the Standardized Mean Difference (AIPE)

Determines the per-group sample size needed for a two-independent-groups
design so that the (expected) confidence interval for Cohen's *d*, the
population standardized mean difference, denoted \\\delta\\, is no wider
than a user-specified value. This is the Accuracy in Parameter
Estimation (AIPE) framework of Kelley and Rausch (2006), the
standardized-mean-difference companion to power-based planning via
[`ss_power_smd`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md).
Optionally, supplying `assurance` returns the larger sample size needed
so that the realized interval will be at or below the target width with
that probability rather than just on average.

## Usage

``` r
ss_aipe_smd(delta, conf_level = 0.95, width, assurance = NULL)
```

## Arguments

- delta:

  The supposed value of the population standardized mean difference
  \\\delta\\ the sample size is planned against: a value the researcher
  posits, either a minimally important effect or a value believed to be
  true in the population, never a sample estimate. Echoed in the
  returned table as the `supposed_smd` row.

- conf_level:

  Desired confidence level (i.e., \\1-\alpha\\, where \\\alpha\\ is the
  Type I error rate). Default `0.95`.

- width:

  Desired (full) width of the two-sided confidence interval on
  \\\delta\\.

- assurance:

  Optional probability with which the realized confidence interval is to
  be no wider than `width`. When `NULL` (the default), the planning
  targets the *expected* width; when supplied (e.g., 0.80, 0.90, 0.99),
  the procedure returns the larger *N* that guarantees the desired width
  with that assurance. Must be `NULL` or strictly between 0.50 and 1.

## Value

A `data.frame` with columns `term` and `value`. The first row,
`necessary_n_per_group`, is the necessary per-group sample size *N*; the
remaining rows echo the user-supplied planning inputs `supposed_smd` and
`width` (and `assurance` when supplied), so the assumptions the sample
size was planned under travel with the result. The `supposed_smd` row is
the supposed effect the plan is built on: a value the researcher posits,
either a minimally important effect or a value believed to be true in
the population, never a sample estimate. The confidence level is
reported in the printed footer.

## Warning

The returned value is the sample size *per group*.

## References

Anderson, S. F., & Kelley, K. (2024). Sample size planning for
replication studies: The devil is in the design. *Psychological Methods,
29*(5), 844–867.
[doi:10.1037/met0000520](https://doi.org/10.1037/met0000520)

Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
planning for more accurate statistical power: A method adjusting sample
effect sizes for publication bias and uncertainty. *Psychological
Science, 28*(11), 1547–1562.
[doi:10.1177/0956797617723724](https://doi.org/10.1177/0956797617723724)

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and
noncentral distributions. *Educational and Psychological Measurement,
61*(4), 532–574.
[doi:10.1177/0013164401614002](https://doi.org/10.1177/0013164401614002)

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128.

Kelley, K. (2005). The effects of nonnormal distributions on confidence
intervals around the standardized mean difference: Bootstrap and
parametric confidence intervals, *Educational and Psychological
Measurement, 65*, 51–69.
[doi:10.1177/0013164404264850](https://doi.org/10.1177/0013164404264850)

Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
obtaining precision: Delineating methods of sample size planning.
*Evaluation and the Health Professions, 26*(3), 258–287.
[doi:10.1177/0163278703255242](https://doi.org/10.1177/0163278703255242)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons and Chapter 3
on one-way ANOVA.)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`smd_c`](https://yelleknek.github.io/DMAR/reference/smd_c.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md),
[`stats::power.t.test()`](https://rdrr.io/r/stats/power.t.test.html)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
ss_aipe_smd(delta = .5, conf_level = .95, width = .30)
#>  term                  value
#>  necessary_n_per_group 353  
#>  supposed_smd          0.5  
#>  width                 0.3  
#> 
#> Confidence level: 95%
ss_aipe_smd(delta = .5, conf_level = .95, width = .30, assurance = .8)
#>  term                  value
#>  necessary_n_per_group 356  
#>  supposed_smd          0.5  
#>  width                 0.3  
#>  assurance             0.8  
#> 
#> Confidence level: 95%
ss_aipe_smd(delta = .5, conf_level = .95, width = .30, assurance = .95)
#>  term                  value
#>  necessary_n_per_group 359  
#>  supposed_smd          0.5  
#>  width                 0.3  
#>  assurance             0.95 
#> 
#> Confidence level: 95%
```
