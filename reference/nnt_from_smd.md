# Number Needed to Treat From a Standardized Mean Difference

Converts a standardized mean difference to the success rate difference
(SRD) of Kraemer and Kupfer (2006) and to its reciprocal, which they
define as the number needed to treat (NNT) for a continuous outcome. The
SRD is the probability that a randomly chosen treated person has a
better outcome than a randomly chosen control person, minus the
probability of the reverse. When sample sizes or confidence limits for
the standardized mean difference are given, the function also returns
the confidence interval for the SRD and the corresponding confidence set
for the NNT, which is two intervals when the interval for the
standardized mean difference includes zero.

## Usage

``` r
nnt_from_smd(
  smd,
  n_1 = NULL,
  n_2 = NULL,
  conf_level = 0.95,
  smd_ci_lower = NULL,
  smd_ci_upper = NULL
)
```

## Arguments

- smd:

  The standardized mean difference (Cohen's *d*): the treatment mean
  minus the control mean, with higher scores the better outcome, divided
  by the pooled within-group standard deviation. A single finite number.
  A negative value gives a negative SRD and a negative NNT, read as the
  number needed to harm.

- n_1, n_2:

  Sample sizes of the treatment and control groups. Supply both or
  neither. With both, the noncentral *t* confidence interval for the
  population standardized mean difference is computed
  ([`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md)) and
  carried to the SRD and the NNT.

- conf_level:

  Confidence level of the noncentral *t* interval when `n_1` and `n_2`
  are given. Default `0.95`. With supplied limits, `conf_level` enters
  no computation: it is the caller's statement of the level of those
  limits, and it is recorded and printed only when it is given in the
  call.

- smd_ci_lower, smd_ci_upper:

  Optional confidence limits for the population standardized mean
  difference, computed elsewhere. Supply both or neither, with
  `smd_ci_lower` not above `smd_ci_upper`. When supplied they are used
  as given, the noncentral *t* computation is skipped, and `n_1` and
  `n_2` are not used. An estimate outside its limits is accepted with a
  warning.

## Value

A `data.frame` (class `dmar_nnt`, a `dmar_tbl`) with columns `term` and
`value`. Without an interval the rows are `smd`, `srd` (the success rate
difference), and `nnt`. With an interval the rows are, in this order,
`smd`, `srd`, `nnt`, `smd_ci_lower`, `smd_ci_upper`, `nnt_ci_lower`,
`nnt_ci_upper`, `srd_ci_lower`, and `srd_ci_upper`.

What each row is: `smd` is the standardized mean difference as given.
`srd` is the success rate difference, \\2 \Phi(d / \sqrt{2}) - 1\\.
`nnt` is `1 / srd`. `smd_ci_lower` and `smd_ci_upper` are the confidence
limits for the population standardized mean difference, computed from
`n_1` and `n_2` or as supplied. `nnt_ci_lower` and `nnt_ci_upper` are
the limits of the NNT when its confidence set is one interval.
`srd_ci_lower` and `srd_ci_upper` are the limits of the SRD. The result
is a `dmar_tbl`, so `tidy()`, `glance()`, and
[`rbind()`](https://rdrr.io/r/base/cbind.html) work on it as on any
other. The confidence level is in the `conf_level` attribute: with `n_1`
and `n_2` it is the level of the computed interval, and with supplied
limits it is recorded only when `conf_level` is given in the call. With
an interval the `ci_source` attribute names `smd` as the estimate whose
interval the others are computed from, so that
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
reports the SRD (`term = "nnt"` or `term = "smd"` for another) and
`tidy()` gives the confidence level on the rows of the SRD and the NNT
and not on those of `smd`.

`srd_ci_lower` and `srd_ci_upper` are the SRD at the two limits of the
standardized mean difference and always form an ordinary interval.
`nnt_ci_lower` and `nnt_ci_upper` depend on where zero falls in the
interval for the standardized mean difference:

- Zero outside the interval: `nnt_ci_lower` is `1 / srd_ci_upper` and
  `nnt_ci_upper` is `1 / srd_ci_lower`.

- Lower limit exactly zero: `nnt_ci_lower` is `1 / srd_ci_upper` and
  `nnt_ci_upper` is `Inf`.

- Upper limit exactly zero and lower limit below zero: `nnt_ci_lower` is
  `-Inf` and `nnt_ci_upper` is `1 / srd_ci_lower`.

- Zero strictly inside the interval: the confidence set for the NNT is
  two intervals, which one lower and one upper limit cannot hold, so
  `nnt_ci_lower` and `nnt_ci_upper` are `NA`. The set is given in the
  `nnt_confidence_set` attribute and in a line printed under the table.

In every case, each finite end of the confidence set for the NNT is
`1 / srd_ci_lower` or `1 / srd_ci_upper`.

A result with an interval carries two attributes. `nnt_set` is
`"interval"` or `"two_intervals"`. `nnt_confidence_set` is a
`data.frame` with columns `lower` and `upper` and one row per interval
of the set, in increasing order: one row for `"interval"`, and for
`"two_intervals"` the rows (`-Inf`, `1 / srd_ci_lower`) and
(`1 / srd_ci_upper`, `Inf`). The two attributes belong to the table
itself: `tidy()` and `glance()` return the rows without them, so for a
set of two intervals they show `nnt_ci_lower` and `nnt_ci_upper` as
`NA`; its finite ends are `1 / srd_ci_lower` and `1 / srd_ci_upper`,
which both carry. A table made by
[`rbind()`](https://rdrr.io/r/base/cbind.html) or by taking rows keeps
the attributes of the first table, as every `dmar_tbl` does, so the line
printed beneath it describes that table.

The NNT of no effect is infinite: `nnt` is `Inf` when `smd` is zero, and
the sign of that infinity has no meaning.

## Details

**What is estimated.** Take one person at random from the treated
population and one at random from the control population, with outcomes
\\Y_T\\ and \\Y_C\\, higher being better. The success rate difference is
\$\$\mathrm{SRD} = \Pr(Y_T \> Y_C) - \Pr(Y_C \> Y_T),\$\$ the
probability that the treated person has the better outcome minus the
probability that the control person has, and Kraemer and Kupfer (2006,
p. 994) define the number needed to treat for a scaled outcome as its
reciprocal, \$\$\mathrm{NNT} = 1 / \mathrm{SRD}.\$\$ No response
threshold enters: the comparison is between two people, not between each
person and a cut. An SRD of 0.25 says that in 100 such pairs the treated
person is ahead in about 25 more pairs than the control person is, and
the NNT of 4 says that about 4 people must be treated for one more of
those comparisons to favor treatment. A negative standardized mean
difference gives a negative SRD and a negative NNT, read as the number
needed to harm. A standardized mean difference of zero gives an SRD of
zero and an infinite NNT (no advantage of either group).

**Assumptions, by the part each serves.** The point conversion assumes
two independent normal populations and a standardized mean difference
whose standardizer is the root mean square of the two population
standard deviations, which the pooled standard deviation estimates under
a common variance (or with equal sample sizes). The interval from `n_1`
and `n_2` assumes normality and a common variance, the conditions of the
noncentral *t* interval for the standardized mean difference
([`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md);
Kelley, 2007). Supplied limits carry whatever assumptions produced them.
When the populations are not normal, as with a skewed outcome, a ceiling
or a floor, or an ordinal scale, the same standardized mean difference
corresponds to a different \\\Pr(Y_T \> Y_C)\\, and the SRD and the NNT
computed here are then answers for a normal model that does not describe
the data. How far off they are depends on the two distributions and is
not quantified here.
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
estimates the SRD itself from the two samples with no distributional
model (Cliff's delta is the sample version of the same difference of
probabilities), and
[`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md)
estimates \\\Pr(Y_T \> Y_C)\\, counting half of each tie.

**The conversion.** Under the normal model the difference \\Y_T - Y_C\\
is normal with mean \\d\\ and variance 2 in units of the standardizer,
so \$\$\Pr(Y_T \> Y_C) = \Phi(d / \sqrt{2}),\$\$ \$\$\mathrm{SRD} = 2
\Phi(d / \sqrt{2}) - 1, \qquad \mathrm{NNT} = \frac{1}{2 \Phi(d /
\sqrt{2}) - 1}.\$\$ The function computes the SRD as \\\mathrm{sign}(d)
\Pr(\chi^2_1 \le d^2 / 2)\\, which is the same number written without
the subtraction, so that a standardized mean difference near zero loses
no digits to it. The values below agree with Table 1 of Kraemer and
Kupfer (2006, p. 992) to the decimals printed there, with one exception:
at *d* = 1.0 the table prints an SRD of .521, and the value is 0.5205,
which is 0.520 to three decimals.

|     |       |        |
|-----|-------|--------|
| *d* | SRD   | NNT    |
| 0.1 | 0.056 | 17.739 |
| 0.2 | 0.112 | 8.892  |
| 0.3 | 0.168 | 5.953  |
| 0.4 | 0.223 | 4.490  |
| 0.5 | 0.276 | 3.619  |
| 0.6 | 0.329 | 3.043  |
| 0.7 | 0.379 | 2.636  |
| 0.8 | 0.428 | 2.334  |
| 0.9 | 0.475 | 2.103  |
| 1.0 | 0.520 | 1.921  |

**The confidence set, and why it can be two intervals.** The function
starts from a confidence interval for the population standardized mean
difference, either the noncentral *t* interval computed from `n_1` and
`n_2` or the limits supplied. Because the SRD is an increasing function
of the standardized mean difference, the SRD at the two limits is a
confidence interval for the SRD with the same coverage, and it is always
an ordinary interval. The NNT is the reciprocal of the SRD, and the
reciprocal is not continuous at zero, so what the interval for the SRD
becomes depends on where zero falls:

- Zero outside the interval. The reciprocal reverses the order of the
  limits: the interval for the NNT runs from \\1 /
  \mathrm{SRD}\_{\mathrm{upper}}\\ to \\1 /
  \mathrm{SRD}\_{\mathrm{lower}}\\.

- Lower limit exactly zero. The interval for the NNT runs from \\1 /
  \mathrm{SRD}\_{\mathrm{upper}}\\ to infinity.

- Upper limit exactly zero and lower limit below zero. The interval runs
  from minus infinity to \\1 / \mathrm{SRD}\_{\mathrm{lower}}\\.

- Zero strictly inside the interval. The SRD values from the lower limit
  up to zero are sent to the NNT values from \\1 /
  \mathrm{SRD}\_{\mathrm{lower}}\\ (negative) down to minus infinity,
  and the SRD values from zero up to the upper limit are sent to the NNT
  values from plus infinity down to \\1 /
  \mathrm{SRD}\_{\mathrm{upper}}\\. The confidence set is therefore two
  intervals, \\\mathrm{NNT} \le 1 / \mathrm{SRD}\_{\mathrm{lower}}\\
  (harm) and \\\mathrm{NNT} \ge 1 / \mathrm{SRD}\_{\mathrm{upper}}\\
  (benefit), which meet at infinity, the NNT of no effect (Altman,
  1998).

- No interval requested. Only `smd`, `srd`, and `nnt` are returned.

One lower limit and one upper limit cannot hold a set of two intervals.
Reading the two reciprocals as the ends of a single interval would
describe the values between them, which are exactly the NNT values the
data exclude, and would leave out the NNT of no effect, which the data
include. For that case `nnt_ci_lower` and `nnt_ci_upper` are therefore
`NA`, the set is in the `nnt_confidence_set` attribute, and a line under
the printed table states it. In every case the confidence set for the
NNT is the image of the interval for the standardized mean difference
under a one-to-one transformation, so it holds the population NNT
exactly when that interval holds the population standardized mean
difference, and it has the same coverage. Kraemer and Kupfer (2006) call
confidence intervals for the NNT “not advisable” (pp. 992 to 993) and
prefer the SRD for interval estimation (p. 994), which is why the limits
of the SRD are reported here beside those of the NNT.

**Edge cases.**

- `smd = 0`: the SRD is 0 and `nnt` is `Inf`. The sign of that infinity
  has no meaning.

- A negative `smd`: the SRD and the NNT are negative, and with both
  limits negative the interval for the NNT is an ordered negative
  interval (for \\d = -0.5\\ with 50 per group, \\-17.636\\ to
  \\-2.109\\).

- A limit exactly zero: one end of the interval for the NNT is infinite,
  as above; the set is one interval.

- A very small \\\|d\|\\: the NNT is very large and is computed without
  loss of digits (for \\d = 10^{-8}\\, about \\1.77 \times 10^{8}\\). An
  SRD that underflows to zero gives an infinite NNT with the sign of
  \\d\\.

- Very small samples: the noncentral *t* interval is wide, includes zero
  for most values of \\d\\, and gives two intervals (for \\d = 0.5\\
  with `n_1 = 1` and `n_2 = 2`, limits for \\d\\ of \\-2.07\\ and 2.88).

- An estimate outside the supplied limits is accepted with a warning,
  and the rows are computed from the three numbers as given.

- Supplied limits take precedence over `n_1` and `n_2`, which are then
  neither used nor checked.

- `conf_level` must be a single number in (0, 1) even when it enters no
  computation; `NULL` is an error.

**Which NNT this is.** The NNT returned here depends on the standardized
mean difference alone. The NNT of evidence based medicine is the
reciprocal of a difference between two response proportions at a stated
threshold. Under the same normal model, with control response proportion
\\p_0\\, that NNT is \$\$1 / \bigl\[\Phi\bigl(\Phi^{-1}(p_0) + d\bigr) -
p_0\bigr\],\$\$ which changes with \\p_0\\ and, for a positive *d*, is
larger than the Kraemer and Kupfer NNT at every threshold. Its smallest
value is \\1 / \[2 \Phi(d / 2) - 1\]\\, at \\p_0 = \Phi(-d / 2)\\. At
\\d = 0.5\\ the Kraemer and Kupfer NNT is 3.62; the threshold NNT is at
least 5.07 (at \\p_0 = 0.40\\), 5.22 at \\p_0 = 0.5\\, and 6.01 at \\p_0
= 0.2\\. The two are different quantities and neither is an
approximation to the other, so a report should say which one it gives.
Furukawa and Leucht (2011) compared the two conversions in individual
patient data from antipsychotic trials and found the control event rate
conversion closer to the observed NNTs: it fell inside the observed 95
percent interval in 35 of 36 comparisons, against 20 of 36 for the
conversion used here. For the NNT at a threshold computed from data, see
[`responder_analysis`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md).

**A worked example.** `nnt_from_smd(smd = 0.5, n_1 = 50, n_2 = 50)`
returns an SRD of 0.276 and an NNT of 3.62: in 100 pairs of one treated
and one control person the treated person is ahead in about 28 more
pairs than the control person is, and about 3.6 people must be treated
for one more such comparison to favor treatment. The 95 percent limits
for the standardized mean difference are 0.101 and 0.897, which exclude
zero. The SRD at those limits is 0.057 and 0.474, and the interval for
the NNT is their reciprocals in reverse order, 2.11 to 17.64.
`nnt_from_smd(smd = 0.2, n_1 = 50, n_2 = 50)` returns an SRD of 0.112
and an NNT of 8.89, with limits for the standardized mean difference of
\\-0.193\\ and 0.592, which include zero. The interval for the SRD,
\\-0.109\\ to 0.325, is an ordinary interval and is the one to report.
`nnt_ci_lower` and `nnt_ci_upper` are `NA`, and the line under the table
gives the confidence set for the NNT as “NNT \<= -9.189 (harm) or NNT
\>= 3.079 (benefit)”: the data are compatible with a benefit as large as
one more favorable comparison for every 3.1 people treated, with no
effect (an infinite NNT), and with harm as large as one more unfavorable
comparison for every 9.2 people treated. A report of this result gives
the estimate, 8.89, with both intervals of the set, or gives the SRD
with its interval.

## References

Altman, D. G. (1998). Confidence intervals for the number needed to
treat. *BMJ, 317*(7168), 1309–1312.
[doi:10.1136/bmj.317.7168.1309](https://doi.org/10.1136/bmj.317.7168.1309)

Furukawa, T. A., & Leucht, S. (2011). How to obtain NNT from Cohen's
*d*: Comparison of two methods. *PLoS ONE, 6*(4), e19070.
[doi:10.1371/journal.pone.0019070](https://doi.org/10.1371/journal.pone.0019070)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24. (The noncentral *t* interval for the
standardized mean difference that is carried to the SRD and the NNT
here.)
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kraemer, H. C., & Kupfer, D. J. (2006). Size of treatment effects and
their importance to clinical research and practice. *Biological
Psychiatry, 59*(11), 990–996.
[doi:10.1016/j.biopsych.2005.09.014](https://doi.org/10.1016/j.biopsych.2005.09.014)

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) and
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) for the
standardized mean difference and its interval;
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md) for
\\\Pr(Y_T \> Y_C)\\ under the same normal model;
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
and
[`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md)
for the same comparisons estimated without the normal model;
[`responder_analysis`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md)
for the NNT at a threshold.

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
[`eta_squared()`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
[`expected_partial_r()`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md),
[`expected_r()`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
[`expected_smd()`](https://yelleknek.github.io/DMAR/reference/expected_smd.md),
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Point estimate only.
nnt_from_smd(smd = 0.5)
#>  term value
#>  smd  0.500
#>  srd  0.276
#>  nnt  3.619

# 2. With the noncentral t interval from the sample sizes.
nnt_from_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.95)
#>  term         value 
#>  smd          0.500 
#>  srd          0.276 
#>  nnt          3.619 
#>  smd_ci_lower 0.101 
#>  smd_ci_upper 0.897 
#>  nnt_ci_lower 2.109 
#>  nnt_ci_upper 17.636
#>  srd_ci_lower 0.057 
#>  srd_ci_upper 0.474 
#> 
#> Confidence level: 95%

# 3. With limits for the standardized mean difference computed elsewhere.
nnt_from_smd(smd = 0.5, smd_ci_lower = 0.20, smd_ci_upper = 0.80)
#>  term         value
#>  smd          0.500
#>  srd          0.276
#>  nnt          3.619
#>  smd_ci_lower 0.200
#>  smd_ci_upper 0.800
#>  nnt_ci_lower 2.334
#>  nnt_ci_upper 8.892
#>  srd_ci_lower 0.112
#>  srd_ci_upper 0.428

# 4. An interval for the standardized mean difference that includes zero:
#    two intervals for the NNT, printed under the table.
nnt_from_smd(smd = 0.3, smd_ci_lower = -0.10, smd_ci_upper = 0.70)
#>  term         value 
#>  smd          0.300 
#>  srd          0.168 
#>  nnt          5.953 
#>  smd_ci_lower -0.100
#>  smd_ci_upper 0.700 
#>  nnt_ci_lower <NA>  
#>  nnt_ci_upper <NA>  
#>  srd_ci_lower -0.056
#>  srd_ci_upper 0.379 
#> 
#> The interval for the SRD includes 0, so the confidence set for the NNT is two intervals:
#> NNT <= -17.739 (harm) or NNT >= 2.636 (benefit), joined at infinity (Altman, 1998)

# 5. The same case from sample sizes, with the set as an attribute: one
#    row per interval of the set.
res <- nnt_from_smd(smd = 0.2, n_1 = 50, n_2 = 50)
res
#>  term         value 
#>  smd          0.200 
#>  srd          0.112 
#>  nnt          8.892 
#>  smd_ci_lower -0.193
#>  smd_ci_upper 0.592 
#>  nnt_ci_lower <NA>  
#>  nnt_ci_upper <NA>  
#>  srd_ci_lower -0.109
#>  srd_ci_upper 0.325 
#> 
#> Confidence level: 95%
#> 
#> The interval for the SRD includes 0, so the confidence set for the NNT is two intervals:
#> NNT <= -9.189 (harm) or NNT >= 3.079 (benefit), joined at infinity (Altman, 1998)
attr(res, "nnt_set")
#> [1] "two_intervals"
attr(res, "nnt_confidence_set")
#>      lower     upper
#> 1     -Inf -9.188745
#> 2 3.079331       Inf

# 6. A negative standardized mean difference: the number needed to harm,
#    with an ordered negative interval.
nnt_from_smd(smd = -0.5, n_1 = 50, n_2 = 50)
#>  term         value  
#>  smd          -0.500 
#>  srd          -0.276 
#>  nnt          -3.619 
#>  smd_ci_lower -0.897 
#>  smd_ci_upper -0.101 
#>  nnt_ci_lower -17.636
#>  nnt_ci_upper -2.109 
#>  srd_ci_lower -0.474 
#>  srd_ci_upper -0.057 
#> 
#> Confidence level: 95%
```
