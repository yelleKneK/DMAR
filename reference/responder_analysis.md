# Responder Analysis: Who Cleared the Threshold, by Group

Reports the proportion of each group whose outcome reaches a meaningful
threshold (a minimal clinically important difference, a remission cut, a
mastery criterion), a clinical and behavioral endpoint that a mean
difference does not show. For each group the function reports the
responder count and proportion with a Wilson confidence interval; with
exactly two groups it adds the risk difference with the Newcombe (1998)
score-based hybrid interval and the number needed to treat with its
confidence set; and across any number of groups it reports the omnibus
chi square test of equal responder proportions. An optional `sweep`
repeats the analysis over a grid of thresholds, making explicit how
conclusions depend on where the line is drawn, disclosing the threshold
dependence that any single-threshold claim leaves implicit.

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
  is the reference for the two-group difference). An observation with a
  missing outcome or a missing label is removed, and a level left with
  no observation is dropped.

- threshold:

  The cut defining response: a single number, fixed before the data are
  seen (see *Assumptions*).

- direction:

  `"ge"` (default): a responder has `x >= threshold`; `"le"`:
  `x <= threshold` (for outcomes where lower is better). An outcome
  exactly at the threshold is a responder under either direction.

- conf_level:

  Confidence level for all intervals. Defaults to 0.95.

- sweep:

  Optional numeric vector of additional thresholds; the analysis is
  repeated at each and stacked with a leading `threshold` column.
  `threshold` comes first, then the values of `sweep` in the order
  given, each distinct value once.

## Value

A tidy wide `data.frame` (class `dmar_responder`, a `dmar_tbl`) with
columns `group`, `n`, `responders`, `estimate`, `ci_lower`, `ci_upper`,
`chi_square`, `df`, and `p_value`, and these rows:

- One row per group, named by its level: the group size `n`, the number
  of `responders`, the responder proportion in `estimate`, and its
  Wilson limits. The three test columns are `NA`.

- `difference` (two groups only): the second proportion minus the first
  in `estimate`, with the Newcombe limits. `n`, `responders`, and the
  test columns are `NA`.

- `nnt` (two groups only): \\1/\|d\|\\ in `estimate` (`Inf` when the
  difference is 0), with its limits when the interval for the difference
  excludes zero and `NA` limits otherwise. The other columns are `NA`.

- `omnibus`: `chi_square`, `df`, and `p_value` of the test of equal
  proportions. The other columns are `NA`.

When `sweep` is supplied, the same table is stacked per threshold with a
leading `threshold` column.

With two groups the result carries two attributes. `nnt_set` is
`"interval"` or `"two_intervals"`, one value per threshold.
`nnt_confidence_set` is a `data.frame` with columns `lower` and `upper`
(preceded by `threshold` under `sweep`) and one row per interval of the
set, on the signed scale \\1/\Delta\\: one row, (`1 / ci_upper`,
`1 / ci_lower`) of the `difference` row, for `"interval"`; and for
`"two_intervals"` the rows (`-Inf`, `1 / ci_lower`) and (`1 / ci_upper`,
`Inf`). For two intervals a line printed under the table states the set,
and names each side by the group that responds more often there, from
the `group_labels` attribute (the two levels of `group`, in order). The
attributes belong to the table itself: `tidy()` and `glance()` return
the rows without them, and a table made by
[`rbind()`](https://rdrr.io/r/base/cbind.html) or by taking rows keeps
the attributes of the first table, as every `dmar_tbl` does. The
confidence level is in the `conf_level` attribute.

## Details

**What is estimated.** For each group, the proportion of its population
whose outcome is at or beyond the threshold, \\\pi_g\\. With two groups,
the difference of the two proportions, \\\Delta = \pi_2 - \pi_1\\
(second level of `group` minus the first), and the number needed to
treat, \\1/\|\Delta\|\\: the number of people who would have to be moved
from the first condition to the second for one more of them, on average,
to be a responder. This is the number needed to treat of evidence based
medicine, a function of two response proportions at a stated threshold.
It is not the quantity
[`nnt_from_smd`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md)
returns, which is computed from a standardized mean difference under a
normal model, involves no threshold, and is smaller for a positive
effect (see *Which NNT this is* on that page).

**Assumptions.** The observations are independent, and with two groups
the groups are independent of each other (two arms of a randomized
trial, not two occasions for the same people; for paired proportions see
[`ci_proportion_paired`](https://yelleknek.github.io/DMAR/reference/ci_proportion_paired.md)).
The threshold is fixed in advance. Every interval and the test condition
on the threshold as a known constant, so a threshold chosen after the
data are seen, for instance the one at which the difference is largest,
makes the stated confidence level and the *p*-value wrong, and no
adjustment is made for it here. `sweep` describes how the result depends
on the threshold; it is not a search, and the rows of a sweep are not
adjusted for one another.

**The interval for each proportion.** With \\r\\ responders among \\n\\,
\\p = r/n\\, and \\z\\ the \\1 - \alpha/2\\ quantile of the standard
normal distribution, the Wilson score limits are \$\$\frac{p + z^2/(2n)
\mp z \sqrt{p(1 - p)/n + z^2/(4n^2)}}{1 + z^2/n},\$\$ the same limits
[`ci_proportion`](https://yelleknek.github.io/DMAR/reference/ci_proportion.md)
returns. They lie inside \[0, 1\] for every count, the counts 0 and
\\n\\ included, where the limit at the boundary is exactly 0 or 1 (and
can be computed as a number within about \\10^{-16}\\ of it).

**The interval for the difference.** The two-group risk difference uses
Newcombe's method 10: the difference interval is assembled from the two
Wilson limits, which keeps it inside \[-1, 1\] and well behaved at
boundary counts. With \\(l_1, u_1)\\ and \\(l_2, u_2)\\ the Wilson
limits of the first and second group and \\d = p_2 - p_1\\, the limits
computed are \$\$d - \sqrt{(p_2 - l_2)^2 + (u_1 - p_1)^2} \quad
\mbox{and} \quad d + \sqrt{(u_2 - p_2)^2 + (p_1 - l_1)^2}.\$\$

**The number needed to treat and its confidence set.** The `nnt` row
holds \\1/\|d\|\\, an absolute value; the direction is the sign of the
`difference` row. Its limits come from the limits of the difference, and
what they are depends on where zero falls:

- Zero outside the interval for the difference: the reciprocals of the
  two limits, in absolute value, are the limits of the number needed to
  treat, the smaller in `ci_lower` and the larger in `ci_upper`.

- Zero inside the interval: a difference of zero is a number needed to
  treat of infinity, so the confidence set is every value at least as
  large as the reciprocal of one limit in one direction, together with
  every value at least as large as the reciprocal of the other limit in
  the other direction. These are two intervals that meet at infinity
  (Altman, 1998), which one lower and one upper limit cannot hold, so
  the two limits of the `nnt` row are `NA`.

In both cases the set is recorded in the `nnt_confidence_set` attribute
on the signed scale \\1/\Delta\\, where a negative value means that the
first group responds more often, and when it is two intervals a line
under the table states it, as
[`nnt_from_smd`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md)
does. The set has the coverage of the interval for the difference,
because it is the image of that interval under the reciprocal.

**The omnibus test.** Pearson's chi square test of equal responder
proportions across the groups, without a continuity correction
([`chisq.test`](https://rdrr.io/r/stats/chisq.test.html) with
`correct = FALSE`), on the number of groups minus 1 degrees of freedom.
With two groups it tests \\\Delta = 0\\. It and the Newcombe interval
are different constructions and can disagree near the boundary of
significance.

**Edge cases.**

- Equal observed proportions: the difference is 0, the `nnt` row is
  `Inf`, its limits are `NA`, and the confidence set is two intervals.

- No responder, or only responders, in one group: the Wilson limits stay
  inside \[0, 1\] and every row is reported. With none in one group and
  all in the other the difference is 1 or \\-1\\ and the number needed
  to treat is 1.

- No responder, or only responders, in every group: the proportions and
  their limits are reported, the difference is 0, and the chi square
  statistic is 0/0, so `chi_square` and `p_value` are `NaN`.

- A group with one observation: its Wilson interval is reported (for one
  responder of one at the 95 percent level, 0.207 to 1), and so is
  everything built from it. The chi square approximation is poor with
  expected counts this small; no warning is given.

- Outcomes equal to the threshold count as responders under both
  directions, so the counts under `"ge"` and `"le"` sum to more than the
  group size when there are such ties.

- Missing values in `x` or `group` are removed before anything is
  counted. Fewer than two groups after that is an error, as are fewer
  than four values in `x`.

- More than two groups: one row per group and the `omnibus` row. There
  is no `difference` or `nnt` row and no `nnt_confidence_set` attribute.

**A worked example.** In the first example below, 30 of 60 control
participants and 35 of 60 treated participants gained at least 10
points. The responder proportions are 0.500 (95 percent Wilson interval
0.377 to 0.623) and 0.583 (0.457 to 0.699). The difference is 0.083,
with the Newcombe interval \\-0.093\\ to 0.252. The number needed to
treat is 1/0.083 = 12. Because the interval for the difference includes
zero, the `nnt` row has no limits; the confidence set is the two
intervals whose finite ends are the reciprocals of the limits of the
difference, \\1/(-0.093) = -10.8\\ and \\1/0.252 = 3.97\\: a number
needed to treat of 3.97 or more in favor of treatment, or of 10.8 or
more in favor of control. The chi square statistic is 0.839 on 1 degree
of freedom, *p* = 0.360. The second example repeats the analysis at
thresholds of 5 and 15: at 5 the difference is 0.183 (0.026 to 0.330),
its interval excludes zero, and the number needed to treat is 5.45 with
limits 3.03 and 37.98; at 15 the difference is 0.067 and its interval
includes zero again. One trial thus gives three different answers at
three thresholds, which is the reason to fix the threshold first and to
show the sweep.

**What the threshold costs.** Dichotomizing an outcome discards
information. For two variables that are bivariate normal, MacCallum,
Zhang, Preacher, and Rucker (2002, pp. 24 to 26) show that splitting one
of them at its mean multiplies their correlation by 0.798 and their
squared correlation by 0.637, which for a test of zero correlation is
about the same as discarding 36 percent of the sample, and that the loss
grows as the split moves away from the mean. A responder threshold is
such a split of the outcome. A responder analysis therefore complements,
and never replaces, the analysis of the continuous outcome (Maxwell,
Delaney, & Kelley, 2027).

## References

Altman, D. G. (1998). Confidence intervals for the number needed to
treat. *BMJ, 317*(7168), 1309–1312.
[doi:10.1136/bmj.317.7168.1309](https://doi.org/10.1136/bmj.317.7168.1309)

MacCallum, R. C., Zhang, S., Preacher, K. J., & Rucker, D. D. (2002). On
the practice of dichotomization of quantitative variables.
*Psychological Methods, 7*(1), 19–40.
[doi:10.1037/1082-989X.7.1.19](https://doi.org/10.1037/1082-989X.7.1.19)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Newcombe, R. G. (1998). Interval estimation for the difference between
independent proportions: Comparison of eleven methods. *Statistics in
Medicine, 17*(8), 873–890.
[doi:10.1002/(SICI)1097-0258(19980430)17:8\<873::AID-SIM779\>3.0.CO;2-I](https://doi.org/10.1002/%28SICI%291097-0258%2819980430%2917%3A8%3C873%3A%3AAID-SIM779%3E3.0.CO%3B2-I)

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
[`expected_partial_r()`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md),
[`expected_r()`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
[`expected_smd()`](https://yelleknek.github.io/DMAR/reference/expected_smd.md),
[`nnt_from_smd()`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A two-arm trial: change scores, response defined as a gain of 10+.
# The interval for the difference includes zero, so the NNT row has no
# limits and the line under the table gives the two intervals.
set.seed(113)
change <- c(rnorm(60, 8, 9), rnorm(60, 13, 9))
arm    <- rep(c("control", "treatment"), each = 60)
res <- responder_analysis(change, arm, threshold = 10)
res
#>  group      n    responders estimate ci_lower ci_upper chi_square df   p_value
#>  control    60   30         0.500    0.377    0.623    <NA>       <NA> <NA>   
#>  treatment  60   35         0.583    0.457    0.699    <NA>       <NA> <NA>   
#>  difference <NA> <NA>       0.083    -0.093   0.252    <NA>       <NA> <NA>   
#>  nnt        <NA> <NA>       12.000   <NA>     <NA>     <NA>       <NA> <NA>   
#>  omnibus    <NA> <NA>       <NA>     <NA>     <NA>     0.839      1    0.3596 
#> 
#> Confidence level: 95%
#> 
#> The interval for the difference includes 0, so the confidence set for the NNT, on the signed scale 1 / difference, is two intervals:
#> NNT <= -10.806 (control responds more often) or NNT >= 3.966 (treatment responds more often), joined at infinity (Altman, 1998)
attr(res, "nnt_confidence_set")
#>     lower     upper
#> 1    -Inf -10.80588
#> 2 3.96564       Inf

# How threshold-dependent is that conclusion?
responder_analysis(change, arm, threshold = 10, sweep = c(5, 15))
#>  threshold group      n    responders estimate ci_lower ci_upper chi_square
#>  10        control    60   30         0.500    0.377    0.623    <NA>      
#>  10        treatment  60   35         0.583    0.457    0.699    <NA>      
#>  10        difference <NA> <NA>       0.083    -0.093   0.252    <NA>      
#>  10        nnt        <NA> <NA>       12.000   <NA>     <NA>     <NA>      
#>  10        omnibus    <NA> <NA>       <NA>     <NA>     <NA>     0.839     
#>  5         control    60   39         0.650    0.524    0.758    <NA>      
#>  5         treatment  60   50         0.833    0.720    0.907    <NA>      
#>  5         difference <NA> <NA>       0.183    0.026    0.330    <NA>      
#>  5         nnt        <NA> <NA>       5.455    3.034    37.979   <NA>      
#>  5         omnibus    <NA> <NA>       <NA>     <NA>     <NA>     5.263     
#>  15        control    60   18         0.300    0.199    0.425    <NA>      
#>  15        treatment  60   22         0.367    0.256    0.493    <NA>      
#>  15        difference <NA> <NA>       0.067    -0.100   0.229    <NA>      
#>  15        nnt        <NA> <NA>       15.000   <NA>     <NA>     <NA>      
#>  15        omnibus    <NA> <NA>       <NA>     <NA>     <NA>     0.600     
#>  df   p_value
#>  <NA> <NA>   
#>  <NA> <NA>   
#>  <NA> <NA>   
#>  <NA> <NA>   
#>  1    0.3596 
#>  <NA> <NA>   
#>  <NA> <NA>   
#>  <NA> <NA>   
#>  <NA> <NA>   
#>  1    0.0218 
#>  <NA> <NA>   
#>  <NA> <NA>   
#>  <NA> <NA>   
#>  <NA> <NA>   
#>  1    0.4386 
#> 
#> Confidence level: 95%
#> 
#> The interval for the difference includes 0, so the confidence set for the NNT, on the signed scale 1 / difference, is two intervals at 2 of the 3 thresholds:
#> at threshold 10, NNT <= -10.806 (control responds more often) or NNT >= 3.966 (treatment responds more often);
#> at threshold 15, NNT <= -9.979 (control responds more often) or NNT >= 4.375 (treatment responds more often), joined at infinity (Altman, 1998)
```
