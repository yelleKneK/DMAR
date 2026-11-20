# Vargha and Delaney's *A* (Stochastic-Superiority Effect Size)

Takes two independent samples of scores that can be ordered and returns
Vargha and Delaney's (2000) *A*, the probability that a randomly drawn
observation from the first population exceeds a randomly drawn
observation from the second (a tie counting one half), with its standard
error, a confidence interval, and a test of \\A = 1/2\\. The statistic
uses only the order of the observations, so it needs no assumption of
normality or of equal variances and suits ordinal and skewed outcomes.

## Usage

``` r
vargha_delaney_A(
  x,
  y = NULL,
  data = NULL,
  conf_level = 0.95,
  ci_method = c("logit", "wald"),
  test = c("t", "permutation", "z"),
  exact = NULL,
  n_resamples = 10000L,
  seed = NULL
)
```

## Arguments

- x:

  Either a numeric or logical vector of observations from group 1, an
  ordered factor, or a two-sided formula of the form `outcome ~ group`
  (in which case `data` must be supplied and the grouping variable must
  have exactly two observed levels; its first level is group 1). Two
  ordered factors must have the same levels in the same order. An
  unordered factor and a character vector are refused; see *What is
  refused, and why*.

- y:

  Observations from group 2, of the same kind as `x`. Ignored when `x`
  is a formula.

- data:

  Optional data frame containing the variables named in the formula `x`.

- conf_level:

  Confidence coverage for a symmetric interval (default `0.95`).

- ci_method:

  Either `"logit"` (default; Wald interval on the logit of *A* with
  back-transformation, which keeps the limits inside (0, 1)) or `"wald"`
  (untransformed Wald on the original scale, clipped to \[0, 1\]).

- test:

  Reference distribution of the test of \\A = 1/2\\. `"t"` (default)
  refers \\(\hat A - 1/2)/\mathrm{SE}\\ to Student's *t* with
  Satterthwaite degrees of freedom, the test of Brunner and Munzel
  (2000). `"permutation"` refers the same statistic to its distribution
  over reassignments of the scores to the two groups (a studentized
  permutation test; Neubert and Brunner, 2007). `"z"` refers it to the
  standard normal distribution. See *The test*.

- exact:

  Logical. If `NULL` (default), the permutation test enumerates every
  reassignment when `choose(n_1 + n_2, n_1)` is at most 50,000 and
  samples reassignments otherwise. If `TRUE`, enumeration is forced
  (refused above 1,000,000 reassignments). If `FALSE`, random
  reassignments are forced. Used only when `test = "permutation"`.

- n_resamples:

  A single whole number of at least 1: the number of randomly drawn
  reassignments when enumeration is not used. Default `10000L`. Used
  only when `test = "permutation"`.

- seed:

  Optional integer seed for the randomly drawn reassignments. Default
  `NULL`, under which the reassignments are drawn from, and advance, the
  current random number stream; supply an integer for reproducibility.
  When a seed is supplied the RNG state in place before the call is
  restored on exit. Used only when `test = "permutation"`.

## Value

A one-row `data.frame` (classes `dmar_vargha_delaney_A` and `dmar_tbl`)
with these columns, in this order:

- `A`:

  The estimate \\\hat A\\.

- `se`:

  Its standard error, the square root of the placement variance of
  DeLong, DeLong, and Clarke-Pearson (1988); 0 when the estimated
  variance is 0.

- `ci_lower`, `ci_upper`:

  The confidence limits at `conf_level`, formed as `ci_method` says.
  Both equal `A` when the estimated variance is 0, and are then not a
  confidence interval.

- `test_statistic`:

  The test statistic \\W = (\hat A - 1/2)/\mathrm{SE}\\, the same under
  each choice of `test`; the `test` column names the distribution it is
  referred to. The column has this one name under all three tests. `NA`
  when the estimated variance is 0, except under `test = "permutation"`,
  where it is `Inf` or `-Inf` for completely separated groups and 0 when
  every score is equal.

- `p_value`:

  The two-sided *p*-value of \\A = 1/2\\ (stochastic equality) from the
  reference distribution named in `test`. `NA` when the estimated
  variance is 0, except under `test = "permutation"`.

- `n_1`, `n_2`:

  The numbers of non-missing observations in the two groups.

- `ci_method`:

  The interval method used.

- `df`:

  The Satterthwaite degrees of freedom under `test = "t"`; `NA` under
  the other two tests and when the estimated variance is 0.

- `test`:

  The reference distribution used: `"t"`, `"permutation"`, or `"z"`.

The table carries these attributes: `conf_level`, which the console
prints beneath the table and
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
reports; `n_missing`, the numbers of missing values removed from the two
groups, as `c(group_1 = , group_2 = )`; `interval_degenerate`, present
only when the estimated variance is 0, with the value
`"complete_separation"` or `"no_variation"`; and `resampling`, present
only under `test = "permutation"`, a list with `exact` (whether every
reassignment was enumerated), `n_arrangements`
(`choose(n_1 + n_2, n_1)`), and `n_resamples_used` (the number of
reassignments compared). The `tidy()` and `glance()` verbs (see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md))
give the row in their column convention, with `df` and `test` after
`p_value` and a logical column `interval_degenerate`.

## Details

**Estimand and interval.** For independent observations \\X\\ from the
first population and \\Y\\ from the second, \$\$A = \Pr(X \> Y) +
\tfrac{1}{2}\\\Pr(X = Y).\$\$ The sampling units are the \\n_1 + n_2\\
individuals, each observed once and each in one group. The interval is a
large-sample (Wald) confidence interval for the population *A*, formed
on the logit scale by default. A value of \\A = 0.5\\ is stochastic
equality, and \\A \> 0.5\\ says that the first group tends to score
higher. Stochastic equality (Vargha & Delaney, 2000, p. 102) is weaker
than equality of the two distributions: populations that differ in
spread or shape can have \\A = 0.5\\. The test reported here is of \\A =
0.5\\. It is not a test that the distributions are identical, the null
hypothesis of the Wilcoxon-Mann-Whitney test (pp. 108–109). Qualitative
magnitude labels for *A* are not reported, in keeping with the DMAR
convention of reporting effect sizes as numbers with confidence
intervals.

**Assumptions.** The two samples are independent of each other and the
observations within each are independent and identically distributed;
for paired observations use
[`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md).
The response is at least ordinal. The estimate needs nothing more: no
distributional shape, no equality of variances, and ties are allowed.
The standard error, the interval, and the *t* and *z* tests are
large-sample results; *The test* and *When the estimated variance is
zero* say how they behave in small samples.

**Scale and missing values.** *A* depends only on the order of the
observations (Vargha & Delaney, 2000, pp. 101, 107), so any increasing
rescoring of the response gives the same value, and an ordinal response
may be supplied as an ordered factor. Missing values (`NA` and `NaN`)
are removed group by group; the numbers removed are in the `n_missing`
attribute. Infinite values are ordered as usual. Dates and date-times
(classes `Date` and `POSIXct`) are ordered by their numeric values, as
are two time differences (class `difftime`) in the same units. A
one-column matrix and a one-dimensional array are read as the vectors
they are. Each group needs at least two non-missing observations; with
fewer the function stops, since no variance can be estimated.

**What is refused, and why.** The statistic uses the order of the
observations and nothing else, so the function must be given objects
whose order is the order meant. It refuses, with a condition of class
`dmar_scores_not_ordered` (see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)),
an unordered factor (its codes follow the order of
[`levels()`](https://rdrr.io/r/base/levels.html), alphabetical unless
set, so the levels low, mid, high are stored as 2, 3, 1, and *A*
computed from the codes of the example in *Examples* is 0.406 where the
ordered levels give 0.656); a character vector (strings sort
alphabetically, `"10"` before `"9"`); an ordered factor beside anything
that is not an ordered factor with the same levels in the same order
(each factor numbers its own levels from 1); time differences in
different units, or one beside plain numbers (each such object is stored
in its own units); a `Date` beside a `POSIXct`; and a list, data frame,
or matrix of several columns. Nothing is computed in any of these cases.
The message says which case applies and what to supply: an ordered
factor with its levels from lowest to highest,
`factor(x, levels = c("low", "mid", "high"), ordered = TRUE)`, or
numeric scores.

**Sample estimate.** With \\\bar R_X\\ the mean rank of the first sample
in the pooled ranking, mid-ranks for ties, \$\$\hat A = \frac{\bar R_X -
(n_1 + 1)/2}{n_2}\$\$ (Vargha & Delaney, 2000, Equation 14, p. 107).
Equivalently, \\\hat A = U / (n_1 n_2)\\, with \\U\\ the Mann-Whitney
statistic that counts the pairs with \\X_i \> Y_j\\ and half of the tied
pairs. The function computes the mean rank form, whose rounding can
differ from the exact ratio \\U/(n_1 n_2)\\ in about the thirteenth
decimal place.

**Standard error.** Define the placements \$\$V\_{10}(X_i) =
\frac{1}{n_2}\sum\_{j} \psi(X_i, Y_j), \qquad V\_{01}(Y_j) =
\frac{1}{n_1}\sum\_{i} \psi(X_i, Y_j),\$\$ with \\\psi(x, y) = 1,
\tfrac{1}{2}, 0\\ as \\x \> y, =, \<\\. The variance estimate is
\$\$\widehat{\mathrm{Var}}(\hat A) = \frac{S^2\_{10}}{n_1} +
\frac{S^2\_{01}}{n_2},\$\$ where \\S^2\_{10}\\ and \\S^2\_{01}\\ are the
sample variances (divisors \\n_1 - 1\\ and \\n_2 - 1\\) of the
\\V\_{10}\\ and \\V\_{01}\\ placements. This is the placement
(structural components) estimator of DeLong, DeLong, and Clarke-Pearson
(1988, p. 840) for the area under a single receiver operating
characteristic curve; Schüürhuis, Konietschke, and Brunner (2025) state
that the variance estimator of Brunner and Munzel (2000) is identical to
it. It is consistent for two independent samples from any pair of
distributions, with no parametric or equal-variance assumption. It is
not unbiased in general: its expectation exceeds the variance of \\\hat
A\\ by \\(\xi\_{11} - \xi\_{10} - \xi\_{01})/(n_1 n_2) \ge 0\\, where
\\\xi\_{11}\\ is the variance of \\\psi(X, Y)\\ and \\\xi\_{10}\\ and
\\\xi\_{01}\\ are the covariances of two values of \\\psi\\ that share
their \\X\\ or their \\Y\\. The excess is zero when both variables take
only two values. With three continuous observations per group from one
distribution the variance is 7/108 and the estimator averages 8/108.

**Confidence interval.** `ci_method = "logit"` (the default) forms a
Wald interval for \\\mathrm{logit}(A) = \log\\A/(1-A)\\\\ with the delta
method standard error \\\widehat{\mathrm{SE}}(\hat A)/\\\hat A(1 - \hat
A)\\\\, \$\$\mathrm{logit}(\hat A) \pm z\_{1 - \alpha/2}\\
\frac{\widehat{\mathrm{SE}}(\hat A)}{\hat A(1 - \hat A)},\$\$ and
transforms the two limits back with the inverse logit, so they stay
inside (0, 1). Pohar Perme and Manevski (2019) conclude that the DeLong
variance estimator is reliable and that the interval should be formed on
the logit scale. `ci_method = "wald"` returns \\\hat A \pm z\_{1 -
\alpha/2}\\\widehat{\mathrm{SE}}(\hat A)\\, clipped to \\\[0, 1\]\\.
Both intervals use the normal quantile, whatever `test` is.

**Relation to Cliff's delta.** \\A = (\delta + 1)/2\\, and because
\\\mathrm{arctanh}(2A - 1) = \mathrm{logit}(A)/2\\, the logit interval
here and the interval of
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
are the same construction on the two scales. The two functions estimate
the variance differently, however (placement variances here; there the
sample version of Feng and Cliff's, 2004, Equation 2, which has a third
term), so their limits are not the same after the change of scale. For
`x = 1:5` and `y = 3:7`, \\\hat A = 0.18\\ with logit limits 0.034 and
0.576, which are \\-0.931\\ and 0.151 on the scale of \\\delta\\;
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
gives \\\hat\delta = -0.64\\ with limits \\-0.935\\ and 0.177.

**The test.** The statistic is \$\$W = \frac{\hat A -
1/2}{\widehat{\mathrm{SE}}(\hat A)},\$\$ reported in `test_statistic`.
The three choices of `test` share this statistic and differ in what it
is referred to. `"t"` gives \\p = 2\\\Pr(T_f \ge \|W\|)\\ for Student's
*t* with \$\$f = \frac{\left(S^2\_{10}/n_1 + S^2\_{01}/n_2\right)^2}
{\dfrac{(S^2\_{10}/n_1)^2}{n_1 - 1} + \dfrac{(S^2\_{01}/n_2)^2}{n_2 -
1}}\$\$ degrees of freedom, reported in `df`. The statistic, the degrees
of freedom, and the rejection rule are those of Schüürhuis, Konietschke,
and Brunner (2025, Equations 3 and 7, the display that precedes Equation
8, and Equation 8), written here in placement variances; they present
them as the test of Brunner and Munzel (2000), whose statistic is their
Equation 4.8 (p. 21) and whose degrees of freedom are their Equation 5.9
(p. 22), an approximation they recommend for group sizes of at least 10.
`"z"` gives \\p = 2\\\Phi(-\|W\|)\\. `"permutation"` recomputes \\W\\
for reassignments of the \\n_1 + n_2\\ pooled scores to groups of sizes
\\n_1\\ and \\n_2\\: every one of the `choose(n_1 + n_2, n_1)`
reassignments when there are at most 50,000 (or when `exact = TRUE`),
and `n_resamples` randomly drawn ones otherwise, by the rules of
[`randomization_test`](https://yelleknek.github.io/DMAR/reference/randomization_test.md).
The *p*-value is the share of reassignments with \\\|W_b\| \ge \|W\|\\
under enumeration, and \\(1 + r)/(1 + m)\\ with \\r\\ such reassignments
among \\m\\ random ones. A reassignment whose estimated variance is 0
has \\W_b = 0\\ when its \\\hat A\\ is 1/2 and \\W_b\\ infinite, with
the sign of \\\hat A - 1/2\\, otherwise. Neubert and Brunner (2007, p.
5195) define only the one-sided *p*-value of their studentized
permutation test, the share of reassignments whose statistic is at least
the observed one, with the variance recomputed for each. Schüürhuis et
al. describe its two-sided *p*-value as twice the smaller tail share;
the form used here, the share at least as extreme in absolute value, is
the one
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md)
uses, and the two can differ. The permutation test is exact when the two
distributions are equal, since every reassignment is then equally
likely, and approximate when they differ in spread or shape with \\A =
1/2\\.

The *z* reference rejects a true \\A = 1/2\\ too often in small samples,
which is why it is not the default. Rejection rates at a nominal .05
with normal data and \\A = 1/2\\, from this function (10,000 samples per
row for *z* and *t*, Monte Carlo standard error about .002; 2,000
samples for the permutation test, standard error about .005, with 399
random reassignments where the reassignments were not enumerated); the
ratio is the standard deviation of the first group over that of the
second:

|         |         |          |      |      |             |
|---------|---------|----------|------|------|-------------|
| \\n_1\\ | \\n_2\\ | SD ratio | *z*  | *t*  | permutation |
| 5       | 5       | 1        | .086 | .071 | .047        |
| 10      | 10      | 1        | .073 | .056 | .047        |
| 8       | 20      | 1        | .080 | .062 | .045        |
| 20      | 20      | 1        | .060 | .052 | .046        |
| 5       | 5       | 3        | .084 | .068 | .048        |
| 10      | 10      | 3        | .079 | .051 | .052        |
| 8       | 20      | 3        | .090 | .062 | .063        |
| 20      | 8       | 3        | .065 | .052 | .034        |
| 20      | 20      | 3        | .062 | .050 | .052        |

Over these and the other simulations run for this function, which
include five-category ordinal and lognormal data and group sizes from 5
to 30, the rates were .056 to .099 for *z*, .049 to .071 for *t*, and
.032 to .063 for the permutation test. The *t* reference removes most of
the excess from 10 per group upward and part of it below that. These
ranges describe the conditions that were run. With one group of 5 or
fewer and the other five or more times larger, the *t* reference
rejected a true \\A = 1/2\\ in .07 to .11 of samples when the two
distributions were equal (.073 with 5 and 25, .085 with 4 and 20, .112
with 3 and 30; normal data, 20,000 samples each), where the permutation
test is exact (the *z* reference rejected in .105, .124, and .172 of the
same samples). When the smaller group has fewer than about 6 scores,
`test = "permutation"` is therefore the test to use, unless the two
groups also differ greatly in spread, the case described next. When the
small group also had five or six times the standard deviation of the
large one, so that the distributions were not equal, the permutation
test rejected in .08 of samples with 5 and 25 and in .11 with 4 and 20
(4,000 samples each, 999 reassignments). When the groups are completely
separated the estimated variance is 0 and only the permutation test
gives a *p*-value: \\W\\ is infinite, and the *p*-value is the share of
reassignments that are also completely separated (2/252 for 5 and 5
untied scores).

The test and the default interval are separate procedures. The interval
is formed on the logit scale so that it stays inside 0 and 1, and the
test is on the scale of *A*. They can disagree in small samples: for
`x = c(3.5, 5, 6, 7, 8)` and `y = c(1, 2, 3, 4, 4.2)`, \\\hat A = .92\\,
the logit interval is \[.486, .993\], the *t* test gives *p* = .0033,
and the exact permutation test gives *p* = .032. A test on the logit
scale would agree with the interval but rejects far too rarely (.007 to
.041 at a nominal .05 in simulation) and is not offered.

**When the estimated variance is zero.** Two kinds of data give a
variance estimate of exactly 0, and they are the only two. *Complete
separation*: every observation in one group exceeds every observation in
the other, so the estimate is at the end of its range (\\\hat A\\ 0 or
1). *No variation*: every observation has the same value, so the
estimate is 0.5. In both, the two limits are reported equal to the
estimate under either `ci_method`, `se` is 0, `test_statistic` and
`p_value` are `NA` (under `test = "permutation"` the test is reported,
as described above), the table's `interval_degenerate` attribute names
the case, and a warning of class `dmar_interval_degenerate` is given
whose `reason` element is `"complete_separation"` or `"no_variation"`;
the condition also carries `estimate`, `n` (the two group sizes), and
`limit_feng_cliff`. The reported limits are not a confidence interval
with the stated coverage. The zero describes the sample: the variance
estimate measures how much the comparisons differ from one observation
to the next, and without overlap they do not differ. Complete separation
is common when the effect is large and the samples are small. With five
observations per group from normal populations, it occurred in 19
percent of samples when the population value of *A* was 0.86 and in 59
percent when it was 0.96, and across all samples the reported 95 percent
limits then covered the population value 79 and about 40 percent of the
time, against 98 and 96 percent among the samples that were not
separated (20,000 simulated samples each). Feng and Cliff (2004, p. 324,
Equations 6 to 8) treat this case by assuming the largest variance
\\\hat\delta\\ can have, which gives the limit \\(n_b - z^2)/(n_b +
z^2)\\ when \\\hat\delta = 1\\ (and its negative when \\\hat\delta =
-1\\), with \\n_b\\ the smaller group size; they describe the result as
conservative. On the scale of *A* that limit is \\n_b/(n_b + z^2)\\ when
\\\hat A = 1\\ and \\z^2/(n_b + z^2)\\ when \\\hat A = 0\\: 0.566 and
0.434 at the 95 percent level with five in the smaller group. The
warning states that limit. It is not the limit this function reports.

Particular cases. Infinite values are ordered like any others, so a
group made up of `Inf` beside finite scores is completely separated, and
two groups made up of the same infinite value have no variation. Two
constant groups at different values are completely separated. Groups
that only touch, the largest score of one equal to the smallest of the
other, are not separated: the tied pair gives a positive variance and an
ordinary interval, which is wide (`c(3, 4, 5)` against `c(1, 2, 3)`
gives \\\hat A = 0.944\\ with limits 0.475 and 0.997). The mark travels
with the result:
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
writes the estimate and says that there is no confidence interval, the
printed table says so beneath itself, `tidy()` and `glance()` carry it
in the column `interval_degenerate`,
[`randomization_test`](https://yelleknek.github.io/DMAR/reference/randomization_test.md)
passes on the warning of
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
and
[`ss_aipe_cliff_delta_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md)
counts a separated replication as one without an interval. The
randomization test remains valid in both cases.

**Relation to other measures.** *A* is the area under the receiver
operating characteristic curve (AUC) when the scores are used to tell
the two groups apart, and \\\hat A = U/(n_1 n_2)\\ for the Mann-Whitney
*U*. Cliff's \\\delta\\ is \\2A - 1\\. The common language effect size
of McGraw and Wong (1992), as
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md) computes
it, is \\\Phi(\delta/\sqrt{2})\\ with \\\delta\\ the standardized mean
difference, which is \\\Pr(X \> Y)\\ for normal distributions with a
common variance; *A* estimates that probability (with ties split)
without those assumptions.

**Reading the output.** In the second example below the first level of
`treatment`, Control, is group 1, with 246 children against 64. \\\hat A
= 0.429\\: a randomly drawn control child has the higher `iq_8` than a
randomly drawn child from the bloomer group with estimated probability
.43, ties counting one half. The 95 percent logit interval is \[0.351,
0.511\], and the *t* test of \\A = 1/2\\ has \\W = -1.72\\ on 92.5
degrees of freedom, *p* = .089.

## References

Brunner, E., & Munzel, U. (2000). The nonparametric Behrens-Fisher
problem: Asymptotic theory and a small-sample approximation.
*Biometrical Journal, 42*(1), 17–25.
[doi:10.1002/(SICI)1521-4036(200001)42:1\<17::AID-BIMJ17\>3.0.CO;2-U](https://doi.org/10.1002/%28SICI%291521-4036%28200001%2942%3A1%3C17%3A%3AAID-BIMJ17%3E3.0.CO%3B2-U)

DeLong, E. R., DeLong, D. M., & Clarke-Pearson, D. L. (1988). Comparing
the areas under two or more correlated receiver operating characteristic
curves: A nonparametric approach. *Biometrics, 44*(3), 837–845.
[doi:10.2307/2531595](https://doi.org/10.2307/2531595)

Feng, D., & Cliff, N. (2004). Monte Carlo evaluation of ordinal d with
improved confidence interval. *Journal of Modern Applied Statistical
Methods, 3*(2), 322–332.
[doi:10.22237/jmasm/1099267560](https://doi.org/10.22237/jmasm/1099267560)

McGraw, K. O., & Wong, S. P. (1992). A common language effect size
statistic. *Psychological Bulletin, 111*(2), 361–365.
[doi:10.1037/0033-2909.111.2.361](https://doi.org/10.1037/0033-2909.111.2.361)

Neubert, K., & Brunner, E. (2007). A studentized permutation test for
the non-parametric Behrens-Fisher problem. *Computational Statistics &
Data Analysis, 51*(10), 5192–5204.
[doi:10.1016/j.csda.2006.05.024](https://doi.org/10.1016/j.csda.2006.05.024)

Pohar Perme, M., & Manevski, D. (2019). Confidence intervals for the
Mann-Whitney test. *Statistical Methods in Medical Research, 28*(12),
3755–3768.
[doi:10.1177/0962280218814556](https://doi.org/10.1177/0962280218814556)

Schüürhuis, S., Konietschke, F., & Brunner, E. (2025). A new approach to
the nonparametric Behrens-Fisher problem with compatible confidence
intervals. *arXiv*, 2504.01796. <https://arxiv.org/abs/2504.01796>

Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the
CL common language effect size statistics of McGraw and Wong. *Journal
of Educational and Behavioral Statistics, 25*(2), 101–132.
[doi:10.3102/10769986025002101](https://doi.org/10.3102/10769986025002101)

## See also

[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`randomization_test`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Two numeric vectors.
set.seed(113)
x <- rnorm(40, mean = 0.6)
y <- rnorm(40, mean = 0)
vargha_delaney_A(x, y)
#>  A     se    ci_lower ci_upper test_statistic p_value n_1 n_2 ci_method df    
#>  0.666 0.061 0.538    0.774    2.712          0.0082  40  40  logit     77.962
#>  test
#>  t   
#> 
#> Confidence level: 95%

# Formula interface on the pygmalion field experiment. The first factor
# level (Control) forms group 1, so A below 0.5 says a randomly drawn
# control child tends to score below a child from the bloomer group.
vargha_delaney_A(iq_8 ~ treatment, data = pygmalion)
#>  A     se    ci_lower ci_upper test_statistic p_value n_1 n_2 ci_method df    
#>  0.429 0.041 0.351    0.511    -1.718         0.0892  246 64  logit     92.478
#>  test
#>  t   
#> 
#> Confidence level: 95%

# Wald (untransformed) interval rather than logit.
vargha_delaney_A(x, y, ci_method = "wald")
#>  A     se    ci_lower ci_upper test_statistic p_value n_1 n_2 ci_method df    
#>  0.666 0.061 0.546    0.786    2.712          0.0082  40  40  wald      77.962
#>  test
#>  t   
#> 
#> Confidence level: 95%

# The studentized permutation test and the normal reference. With
# 40 per group the reassignments are sampled, so a seed is given.
# n_resamples = 2000 keeps the example quick; a reported p-value
# deserves the default of 10,000.
vargha_delaney_A(x, y, test = "permutation", n_resamples = 2000,
                 seed = 113)
#>  A     se    ci_lower ci_upper test_statistic p_value n_1 n_2 ci_method df  
#>  0.666 0.061 0.538    0.774    2.712          0.0120  40  40  logit     <NA>
#>  test       
#>  permutation
#> 
#> Confidence level: 95%
#> 
#> Permutation test: 2000 random reassignments
vargha_delaney_A(x, y, test = "z")
#>  A     se    ci_lower ci_upper test_statistic p_value n_1 n_2 ci_method df  
#>  0.666 0.061 0.538    0.774    2.712          0.0067  40  40  logit     <NA>
#>  test
#>  z   
#> 
#> Confidence level: 95%

# An ordinal response is supplied as an ordered factor with its levels
# stated from lowest to highest. The unordered factors rating_1 and
# rating_2 themselves are refused, with a message that says why: their
# codes follow the alphabetical order of the labels (high, low, mid).
rating_1 <- factor(c("low", "mid", "high", "high"))
rating_2 <- factor(c("low", "low", "mid", "high"))
lev <- c("low", "mid", "high")
vargha_delaney_A(factor(rating_1, levels = lev, ordered = TRUE),
                 factor(rating_2, levels = lev, ordered = TRUE))
#>  A     se    ci_lower ci_upper test_statistic p_value n_1 n_2 ci_method df test
#>  0.656 0.209 0.237    0.921    0.748          0.4827  4   4   logit     6  t   
#> 
#> Confidence level: 95%

# Completely separated groups: the estimated variance is 0, the two
# limits equal the estimate, and a warning says that they are not a
# confidence interval. The permutation test still gives a p-value, the
# share of the 252 reassignments that are also completely separated.
separated <- vargha_delaney_A(6:10, 1:5, test = "permutation")
#> Warning: vargha_delaney_A(): the groups are completely separated (A = 1), so the estimated variance is 0. Both confidence limits are reported equal to the estimate (1, 1). They are not a confidence interval with the stated coverage: the variance estimate is built from how much the comparisons differ from one observation to the next, and in a sample with no overlap they do not differ at all. That describes this sample, not the population. A population with some overlap often produces a sample with none when the samples are small (here 5 and 5 observations), so the limits must not be read as showing that the population value is 1. For this case Feng and Cliff (2004, Equation 8) give the limit (n_b - z^2) / (n_b + z^2) for delta = 2A - 1, with n_b the smaller group size and z the normal quantile of the level: 0.131 here at the 95% level, with the other limit 1; on the scale of A, 0.566, with the other limit 1. The permutation test is reported (test_statistic is infinite); the t and z tests are not defined here. Report the estimate with the sample sizes and say that the groups did not overlap; the randomization test (randomization_test()) gives a p-value that remains valid here.
separated
#>  A     se    ci_lower ci_upper test_statistic p_value n_1 n_2 ci_method df  
#>  1.000 0.000 1.000    1.000    Inf            0.0079  5   5   logit     <NA>
#>  test       
#>  permutation
#> 
#> Confidence level: 95%
#> 
#> Limits are not a confidence interval: the estimated variance is 0.
#> 
#> Permutation test: all 252 reassignments
results_sentence(separated)
#> [1] "A = 1.00 (no confidence interval: the groups are completely separated, so the estimated variance is 0), unadjusted p = 0.0079"
```
