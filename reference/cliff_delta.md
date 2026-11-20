# Cliff's \\\delta\\ Ordinal Effect Size

Takes two independent samples of scores that can be ordered and returns
Cliff's \\\delta\\, the probability that a randomly drawn observation
from the first population exceeds one from the second minus the
probability of the reverse, with a confidence interval and the
proportions of pairs won, lost, and tied. The statistic uses only the
order of the observations, so it needs no assumption of normality or of
equal variances and no transformation of a skewed outcome, and its
interval comes from a variance formula, not from resampling.

## Usage

``` r
cliff_delta(group_1, group_2, conf_level = 0.95)
```

## Arguments

- group_1, group_2:

  Numeric or logical vectors, or ordered factors with the same levels in
  the same order, holding the observations in the two groups. The
  statistic uses only their order. Missing values are removed group by
  group (see the `n_missing` attribute). An unordered factor and a
  character vector are refused; see *What is refused, and why*.

- conf_level:

  Confidence level for the CI. Default `0.95`.

## Value

A `data.frame` (classes `dmar_cliff_delta` and `dmar_tbl`) with columns
`term` and `value` and these rows:

- `cliff_delta`:

  The estimate \\\hat\delta\\, in \\\[-1, 1\]\\.

- `ci_lower`, `ci_upper`:

  The confidence limits at `conf_level`. Both equal the estimate when
  the estimated variance is 0, and are then not a confidence interval.

- `var_cliff_delta`:

  The estimated variance of \\\hat\delta\\ (see *Variance*).

- `p_y1_greater`:

  The proportion of the \\n_1 n_2\\ pairs, one observation from each
  group, in which the observation from group 1 is the larger.

- `p_y1_less`:

  The proportion in which it is the smaller.

- `p_tied`:

  The proportion of tied pairs. The three proportions sum to 1, and
  `cliff_delta` is the first minus the second.

The table carries the attributes `conf_level`; `n_missing`, the numbers
of missing values removed from the two groups, as
`c(group_1 = , group_2 = )`; and `interval_degenerate`, present only
when the estimated variance is 0, with the value `"complete_separation"`
or `"no_variation"`. The `tidy()` and `glance()` verbs (see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md))
add a logical column `interval_degenerate`.

## Details

**Estimand and interval.** For independent observations \\Y_1\\ from the
first population and \\Y_2\\ from the second, \$\$\delta \\=\\ \Pr(Y_1
\> Y_2) - \Pr(Y_1 \< Y_2) \\=\\ 2 \cdot A - 1,\$\$ where \\A\\ is the
Vargha-Delaney (2000) statistic. The sampling units are the \\n_1 +
n_2\\ individuals, each observed once and each in one group. The
interval is a large-sample confidence interval for the population
\\\delta\\, formed on the inverse hyperbolic tangent scale. \\\delta\\
ranges over \\\[-1, 1\]\\; 0 says that neither group tends to score
higher, and 1 or \\-1\\ that the two distributions do not overlap. Feng
and Cliff (2004, p. 322) credit the statistic to Cliff (1993).

**Sample estimate.** \$\$\hat\delta = \frac{\\\\(i,j): y\_{1i} \>
y\_{2j}\\ - \\\\(i,j): y\_{1i} \< y\_{2j}\\}{n_1 n_2},\$\$ the mean of
the \\n_1 n_2\\ signs \\d\_{ij} = \mathrm{sign}(y\_{1i} - y\_{2j})\\. A
tied pair has sign 0 and counts for neither side. The rows
`p_y1_greater`, `p_y1_less`, and `p_tied` give the three proportions of
pairs.

**Assumptions.** The two samples are independent of each other and the
observations within each are independent and identically distributed;
for paired observations use
[`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md).
The response is at least ordinal. The estimate needs nothing more: no
distributional shape, no equality of variances, and ties are allowed.
The interval is a large-sample result; *When the estimated variance is
zero* says how it behaves in small samples with a large effect.

**Variance.** The variance of \\\hat\delta\\ is
\$\$\mathrm{Var}(\hat\delta) \\=\\ \frac{(n_2 - 1) \sigma^2\_{d_1} +
(n_1 - 1) \sigma^2\_{d_2} + \sigma^2_d}{n_1 n_2}\$\$ (Feng & Cliff,
2004, Equation 2, p. 323, who credit Cliff, 1993), where
\\\sigma^2\_{d_1}\\ is the variance of the row means \\d\_{i\cdot}\\,
the mean sign of one observation of the first group against the whole
second population, \\\sigma^2\_{d_2}\\ is the variance of the column
means \\d\_{\cdot j}\\, and \\\sigma^2_d\\ is the variance of the signs
\\d\_{ij}\\. (The equation as printed there attaches \\(n_1 - 1)\\ to
the variance of the first group's row means and \\(n_2 - 1)\\ to the
second's. The coefficients are the other way around, as the variance of
a two-sample U statistic shows: as \\n_2\\ grows, \\\hat\delta\\ becomes
the mean of \\n_1\\ independent row means, whose variance is
\\\sigma^2\_{d_1}/n_1\\. The difference matters only for unequal group
sizes.) The function estimates the three components by their sample
values: the sample variances (divisors \\n_1 - 1\\ and \\n_2 - 1\\) of
the observed row means and column means, and the sample variance
(divisor \\n_1 n_2 - 1\\) of the \\n_1 n_2\\ signs. Feng and Cliff's
Equation 4, an unbiased estimate with a minimum allowed value, is not
used. The result is the row `var_cliff_delta`.

**Confidence interval.** The interval is formed on the scale of
\\\mathrm{arctanh}(\delta)\\, as Fisher's *Z* interval for a correlation
is, and transformed back, so that the limits stay inside \\\[-1, 1\]\\:
\$\$\tanh\\\left\\\mathrm{arctanh}(\hat\delta) \pm z\_{1 - \alpha/2}\\
\frac{\sqrt{\widehat{\mathrm{Var}}(\hat\delta)}}{1 - \hat\delta^2}
\right\\,\$\$ the standard error on the transformed scale coming from
the delta method. Feng and Cliff (2004, Equation 5, p. 324) give an
asymmetric interval of another form, which this function does not
compute.

**Relation to Vargha and Delaney's *A*.** \\A = (\delta + 1)/2\\, and
because \\\mathrm{arctanh}(2A - 1) = \mathrm{logit}(A)/2\\, the interval
here and the logit interval of
[`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md)
are the same construction on the two scales. The two functions estimate
the variance differently, however (the three-term estimate above here;
the two placement variances of DeLong, DeLong, and Clarke-Pearson, 1988,
there), so their limits are not the same after the change of scale. For
`group_1 = 1:5` and `group_2 = 3:7`, \\\hat\delta = -0.64\\ with limits
\\-0.935\\ and 0.177;
[`vargha_delaney_A()`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md)
gives \\\hat A = 0.18\\ with limits 0.034 and 0.576, which are
\\-0.931\\ and 0.151 on the scale of \\\delta\\.

**Scale and missing values.** The observations enter only through their
order, so any increasing rescoring of the response gives the same
result, and an ordinal response may be supplied as an ordered factor.
Missing values (`NA` and `NaN`) are removed group by group; the numbers
removed are in the `n_missing` attribute. Infinite values are ordered as
usual. Dates and date-times (classes `Date` and `POSIXct`) are ordered
by their numeric values, as are two time differences (class `difftime`)
in the same units. A one-column matrix and a one-dimensional array are
read as the vectors they are. Each group needs at least two non-missing
observations; with fewer the function stops, since no variance can be
estimated.

**What is refused, and why.** The statistic uses the order of the
observations and nothing else, so the function must be given objects
whose order is the order meant. It refuses, with a condition of class
`dmar_scores_not_ordered` (see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)),
an unordered factor (its codes follow the order of
[`levels()`](https://rdrr.io/r/base/levels.html), alphabetical unless
set, so the levels low, mid, high are stored as 2, 3, 1, and
\\\hat\delta\\ computed from the codes of the example in *Examples* is
\\-0.19\\ where the ordered levels give 0.31); a character vector
(strings sort alphabetically, `"10"` before `"9"`); an ordered factor
beside anything that is not an ordered factor with the same levels in
the same order (each factor numbers its own levels from 1); time
differences in different units, or one beside plain numbers (each such
object is stored in its own units); a `Date` beside a `POSIXct`; and a
list, data frame, or matrix of several columns. Nothing is computed in
any of these cases. The message says which case applies and what to
supply: an ordered factor with its levels from lowest to highest,
`factor(x, levels = c("low", "mid", "high"), ordered = TRUE)`, or
numeric scores.

**When the estimated variance is zero.** Two kinds of data give a
variance estimate of exactly 0, and they are the only two. *Complete
separation*: every observation in one group exceeds every observation in
the other, so the estimate is at the end of its range (\\\hat\delta\\
\\-1\\ or 1). *No variation*: every observation has the same value, so
the estimate is 0. In both, the two limits are reported equal to the
estimate, `var_cliff_delta` is 0, the table's `interval_degenerate`
attribute names the case, and a warning of class
`dmar_interval_degenerate` is given whose `reason` element is
`"complete_separation"` or `"no_variation"`; the condition also carries
`estimate`, `n` (the two group sizes), and `limit_feng_cliff`. The
reported limits are not a confidence interval with the stated coverage.
The zero describes the sample: the variance estimate measures how much
the comparisons differ from one observation to the next, and without
overlap they do not differ. Complete separation is common when the
effect is large and the samples are small. With five observations per
group from normal populations, it occurred in 19 percent of samples when
the population value of \\\delta\\ was 0.71 and in 59 percent when it
was 0.92, and across all samples the reported 95 percent limits then
covered the population value 79 and about 40 percent of the time,
against 98 and 96 percent among the samples that were not separated
(20,000 simulated samples each). Feng and Cliff (2004, p. 324, Equations
6 to 8) treat this case by assuming the largest variance \\\hat\delta\\
can have, which gives the limit \\(n_b - z^2)/(n_b + z^2)\\ when
\\\hat\delta = 1\\ (and its negative when \\\hat\delta = -1\\), with
\\n_b\\ the smaller group size; they describe the result as
conservative. At the 95 percent level it is \\-0.123\\, 0.020, 0.131,
0.445, and 0.678 for \\n_b\\ of 3, 4, 5, 10, and 20. The warning states
that limit. It is not the limit this function reports.

Particular cases. Infinite values are ordered like any others, so a
group made up of `Inf` beside finite scores is completely separated, and
two groups made up of the same infinite value have no variation. Two
constant groups at different values are completely separated. Groups
that only touch, the largest score of one equal to the smallest of the
other, are not separated: the tied pair gives a positive variance and an
ordinary interval, which is wide (`c(3, 4, 5)` against `c(1, 2, 3)`
gives \\\hat\delta = 0.889\\ with limits \\-0.167\\ and 0.995). The mark
travels with the result:
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
writes the estimate and says that there is no confidence interval, the
printed table says so beneath itself, `tidy()` and `glance()` carry it
in the column `interval_degenerate`,
[`randomization_test`](https://yelleknek.github.io/DMAR/reference/randomization_test.md)
passes this warning on, and
[`ss_aipe_cliff_delta_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md)
counts a separated replication as one without an interval. The
randomization test remains valid in both cases.

**Relation to other measures.** Cliff's \\\delta\\ is a linear function
of the Vargha-Delaney (2000) \\A\\ statistic (\\\delta = 2A - 1\\) and
of the Mann-Whitney \\U\\ statistic (\\U / (n_1 n_2) = A\\). The common
language effect size computed by
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md) is
\\\Pr(Y_1 \> Y_2)\\ for normal distributions with a common variance;
\\(\delta + 1)/2\\ estimates that probability (with ties split) without
those assumptions, which suits skewed outcomes and ordinal scales.

**Reading the output.** In the third example below the 32 participants
given the Community Reinforcement Approach are group 1 and the 37 given
standard care are group 2, so there are 1,184 pairs. The CRA participant
reports more drinks per week in 29.3 percent of the pairs
(`p_y1_greater`), fewer in 59.3 percent (`p_y1_less`), and the same
number in 11.4 percent (`p_tied`), which gives \\\hat\delta = 0.293 -
0.593 = -0.300\\, with 95 percent limits \\-0.535\\ and \\-0.021\\.

## References

Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer
ordinal questions. *Psychological Bulletin, 114*(3), 494–509.
[doi:10.1037/0033-2909.114.3.494](https://doi.org/10.1037/0033-2909.114.3.494)

DeLong, E. R., DeLong, D. M., & Clarke-Pearson, D. L. (1988). Comparing
the areas under two or more correlated receiver operating characteristic
curves: A nonparametric approach. *Biometrics, 44*(3), 837–845.
[doi:10.2307/2531595](https://doi.org/10.2307/2531595)

Feng, D., & Cliff, N. (2004). Monte Carlo evaluation of ordinal d with
improved confidence interval. *Journal of Modern Applied Statistical
Methods, 3*(2), 322–332.
[doi:10.22237/jmasm/1099267560](https://doi.org/10.22237/jmasm/1099267560)

Smith, J. E., Meyers, R. J., & Delaney, H. D. (1998). The community
reinforcement approach with homeless alcohol-dependent individuals.
*Journal of Consulting and Clinical Psychology, 66*(3), 541–548.

Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the
CL common language effect size statistics of McGraw and Wong. *Journal
of Educational and Behavioral Statistics, 25*(2), 101–132.
[doi:10.3102/10769986025002101](https://doi.org/10.3102/10769986025002101)

## See also

[`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md),
[`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`randomization_test`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`nnt_from_smd`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`ss_aipe_cliff_delta`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md)

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
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
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Two groups of different sizes, no ties:
set.seed(113)
a <- rnorm(30, mean = 0, sd = 1)
b <- rnorm(40, mean = 0.5, sd = 1)
cliff_delta(a, b)
#>  term            value 
#>  cliff_delta     -0.223
#>  ci_lower        -0.470
#>  ci_upper        0.056 
#>  var_cliff_delta 0.019 
#>  p_y1_greater    0.388 
#>  p_y1_less       0.612 
#>  p_tied          0.000 
#> 
#> Confidence level: 95%

# 2. With ties (ordinal data):
o1 <- c(1, 2, 2, 3, 3, 3, 4, 4, 5)
o2 <- c(2, 3, 3, 4, 4, 5, 5, 5)
cliff_delta(o1, o2)
#>  term            value 
#>  cliff_delta     -0.403
#>  ci_lower        -0.776
#>  ci_upper        0.179 
#>  var_cliff_delta 0.067 
#>  p_y1_greater    0.194 
#>  p_y1_less       0.597 
#>  p_tied          0.208 
#> 
#> Confidence level: 95%

# 3. Robust to right skew. Cliff's delta on the raw, untransformed
#    drinking outcome from the Smith, Meyers, and Delaney (1998)
#    trial, comparing the Community Reinforcement Approach (CRA)
#    against standard care. Because the statistic uses only ranks it
#    needs no normalizing transformation of the heavily skewed
#    outcome, unlike the standardized mean difference.
data(drinks_trial)
cra <- drinks_trial$drinks_per_week[drinks_trial$treatment == "CRA"]
std <- drinks_trial$drinks_per_week[drinks_trial$treatment == "Standard"]
cliff_delta(cra, std)
#>  term            value 
#>  cliff_delta     -0.300
#>  ci_lower        -0.535
#>  ci_upper        -0.021
#>  var_cliff_delta 0.018 
#>  p_y1_greater    0.293 
#>  p_y1_less       0.593 
#>  p_tied          0.114 
#> 
#> Confidence level: 95%

# 4. An ordinal response is supplied as an ordered factor with its
#    levels stated from lowest to highest. The unordered factors
#    rating_1 and rating_2 themselves are refused, with a message that
#    says why: their codes follow the alphabetical order of the labels
#    (high, low, mid).
rating_1 <- factor(c("low", "mid", "high", "high"))
rating_2 <- factor(c("low", "low", "mid", "high"))
lev <- c("low", "mid", "high")
cliff_delta(factor(rating_1, levels = lev, ordered = TRUE),
            factor(rating_2, levels = lev, ordered = TRUE))
#>  term            value 
#>  cliff_delta     0.312 
#>  ci_lower        -0.517
#>  ci_upper        0.839 
#>  var_cliff_delta 0.170 
#>  p_y1_greater    0.500 
#>  p_y1_less       0.188 
#>  p_tied          0.312 
#> 
#> Confidence level: 95%

# 5. Completely separated groups: the estimated variance is 0, the two
#    limits equal the estimate, and a warning says that they are not a
#    confidence interval and gives the limit of Feng and Cliff (2004).
separated <- cliff_delta(6:10, 1:5)
#> Warning: cliff_delta(): the groups are completely separated (|delta| = 1), so the estimated variance is 0. Both confidence limits are reported equal to the estimate (1, 1). They are not a confidence interval with the stated coverage: the variance estimate is built from how much the comparisons differ from one observation to the next, and in a sample with no overlap they do not differ at all. That describes this sample, not the population. A population with some overlap often produces a sample with none when the samples are small (here 5 and 5 observations), so the limits must not be read as showing that the population value is 1. For this case Feng and Cliff (2004, Equation 8) give the limit (n_b - z^2) / (n_b + z^2) for delta, with n_b the smaller group size and z the normal quantile of the level: 0.131 here at the 95% level, with the other limit 1. Report the estimate with the sample sizes and say that the groups did not overlap; the randomization test (randomization_test()) gives a p-value that remains valid here.
separated
#>  term            value
#>  cliff_delta     1.000
#>  ci_lower        1.000
#>  ci_upper        1.000
#>  var_cliff_delta 0.000
#>  p_y1_greater    1.000
#>  p_y1_less       0.000
#>  p_tied          0.000
#> 
#> Confidence level: 95%
#> 
#> Limits are not a confidence interval: the estimated variance is 0.
results_sentence(separated)
#> [1] "cliff_delta = 1.00 (no confidence interval: the groups are completely separated, so the estimated variance is 0)"
```
