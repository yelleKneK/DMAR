# Responder Analysis: Who Cleared the Threshold, by Group

The clinical and behavioral endpoint that mean differences hide: the
proportion of each group whose outcome reaches a meaningful threshold (a
minimal clinically important difference, a remission cut, a mastery
criterion). For each group the function reports the responder count and
proportion with a Wilson confidence interval; with exactly two groups it
adds the risk difference with the Newcombe (1998) score-based hybrid
interval and the number needed to treat; and across any number of groups
it reports the omnibus chi square test of equal responder proportions.
An optional `sweep` repeats the analysis over a grid of thresholds,
making explicit how conclusions depend on where the line is drawn,
disclosing the threshold dependence that any single-threshold claim
leaves implicit.

## Usage

``` r
responder_analysis(
  x,
  group,
  threshold,
  direction = c("ge", "le"),
  conf_level = 0.95,
  sweep = NULL
)
```

## Arguments

- x:

  Numeric vector of outcomes (for example, change scores).

- group:

  Group labels, one per observation (coerced to factor; the first level
  is the reference for the two-group difference).

- threshold:

  The cut defining response.

- direction:

  `"ge"` (default): a responder has `x >= threshold`; `"le"`:
  `x <= threshold` (for outcomes where lower is better).

- conf_level:

  Confidence level for all intervals. Defaults to 0.95.

- sweep:

  Optional numeric vector of additional thresholds; the analysis is
  repeated at each and stacked with a leading `threshold` column.

## Value

A tidy wide `data.frame` (class `dmar_tbl`). One row per group with
`group`, `n`, `responders`, `estimate` (the proportion), `lower_limit`,
`upper_limit`; with two groups, a `difference` row (second level minus
first) and an `nnt` row; and a final `omnibus` row carrying
`chi_square`, `df`, and `p_value` (columns that are `NA` on the other
rows). When `sweep` is supplied, the same table is stacked per threshold
with a leading `threshold` column.

## Details

Per-group intervals are Wilson score intervals
([`ci_proportion`](https://yelleknek.github.io/DMAR/reference/ci_proportion.md)).
The two-group risk difference uses Newcombe's method 10: the difference
interval is assembled from the two Wilson limits, which keeps it inside
\[-1, 1\] and well behaved at boundary counts. The number needed to
treat is \\1/\|\Delta\|\\, with its interval from the inverted
difference limits when the difference interval excludes zero; when it
includes zero the NNT interval is reported as `NA` (the interval is
disjoint and an interval on the NNT scale would mislead; Altman, 1998).
Dichotomizing throws away information, so a responder analysis
complements, never replaces, the analysis of the continuous outcome
(Maxwell, Delaney, & Kelley, 2027).

## References

Altman, D. G. (1998). Confidence intervals for the number needed to
treat. *BMJ, 317*(7168), 1309–1312.
[doi:10.1136/bmj.317.7168.1309](https://doi.org/10.1136/bmj.317.7168.1309)

Newcombe, R. G. (1998). Interval estimation for the difference between
independent proportions: Comparison of eleven methods. *Statistics in
Medicine, 17*(8), 873–890.
[doi:10.1002/(SICI)1097-0258(19980430)17:8\<873::AID-SIM779\>3.0.CO;2-I](https://doi.org/10.1002/%28SICI%291097-0258%2819980430%2917%3A8%3C873%3A%3AAID-SIM779%3E3.0.CO%3B2-I)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ci_proportion`](https://yelleknek.github.io/DMAR/reference/ci_proportion.md)
for the per-group interval;
[`nnt_from_smd`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md)
for the model-based route to the number needed to treat from a
standardized mean difference;
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
and
[`proportion_of_superiority`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md)
for dominance-style effect sizes on the continuous outcome.

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
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A two-arm trial: change scores, response defined as a gain of 10+.
set.seed(113)
change <- c(rnorm(60, 8, 9), rnorm(60, 13, 9))
arm    <- rep(c("control", "treatment"), each = 60)
responder_analysis(change, arm, threshold = 10)
#>  group      n    responders estimate lower_limit upper_limit chi_square df  
#>  control    60   30         0.5      0.377       0.623       <NA>       <NA>
#>  treatment  60   35         0.583    0.457       0.699       <NA>       <NA>
#>  difference <NA> <NA>       0.0833   -0.0925     0.252       <NA>       <NA>
#>  nnt        <NA> <NA>       12       <NA>        <NA>        <NA>       <NA>
#>  omnibus    <NA> <NA>       <NA>     <NA>        <NA>        0.839      1   
#>  p_value
#>  <NA>   
#>  <NA>   
#>  <NA>   
#>  <NA>   
#>  0.3596 
#> 
#> Confidence level: 95%

# How threshold-dependent is that conclusion?
responder_analysis(change, arm, threshold = 10, sweep = c(5, 15))
#>  threshold group      n    responders estimate lower_limit upper_limit
#>  10        control    60   30         0.5      0.377       0.623      
#>  10        treatment  60   35         0.583    0.457       0.699      
#>  10        difference <NA> <NA>       0.0833   -0.0925     0.252      
#>  10        nnt        <NA> <NA>       12       <NA>        <NA>       
#>  10        omnibus    <NA> <NA>       <NA>     <NA>        <NA>       
#>  5         control    60   39         0.65     0.524       0.758      
#>  5         treatment  60   50         0.833    0.72        0.907      
#>  5         difference <NA> <NA>       0.183    0.0263      0.33       
#>  5         nnt        <NA> <NA>       5.45     3.03        38         
#>  5         omnibus    <NA> <NA>       <NA>     <NA>        <NA>       
#>  15        control    60   18         0.3      0.199       0.425      
#>  15        treatment  60   22         0.367    0.256       0.493      
#>  15        difference <NA> <NA>       0.0667   -0.1        0.229      
#>  15        nnt        <NA> <NA>       15       <NA>        <NA>       
#>  15        omnibus    <NA> <NA>       <NA>     <NA>        <NA>       
#>  chi_square df   p_value
#>  <NA>       <NA> <NA>   
#>  <NA>       <NA> <NA>   
#>  <NA>       <NA> <NA>   
#>  <NA>       <NA> <NA>   
#>  0.839      1    0.3596 
#>  <NA>       <NA> <NA>   
#>  <NA>       <NA> <NA>   
#>  <NA>       <NA> <NA>   
#>  <NA>       <NA> <NA>   
#>  5.26       1    0.0218 
#>  <NA>       <NA> <NA>   
#>  <NA>       <NA> <NA>   
#>  <NA>       <NA> <NA>   
#>  <NA>       <NA> <NA>   
#>  0.6        1    0.4386 
#> 
#> Confidence level: 95%
```
