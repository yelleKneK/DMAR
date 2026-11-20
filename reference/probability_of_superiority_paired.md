# Probability of Superiority for a Paired-Samples Design

Takes paired observations and returns the probability of superiority for
a paired design (Grissom & Kim, 2005, 2012), the probability that the
second member of a pair scores higher than the first (a tied pair
counting one half), with a confidence interval and the numbers of pairs
won, lost, and tied. Only the order within each pair is used. It is the
paired counterpart of
[`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md)
and
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
which are for two independent groups.

## Usage

``` r
probability_of_superiority_paired(x, y, conf_level = 0.95)
```

## Arguments

- x, y:

  Paired numeric vectors of equal length. \\x\\ and \\y\\ are
  interpreted as repeated measurements on the same units (e.g.,
  pre/post, condition 1 / condition 2, sibling-1 / sibling-2). Logical
  vectors and ordered factors with the same levels in the same order are
  also accepted, since only the order within each pair is used. Pairs
  with a missing member are removed (see the `n_missing` attribute). An
  unordered factor and a character vector are refused; see *What is
  refused, and why*.

- conf_level:

  Confidence level. Default `0.95`.

## Value

A `data.frame` (classes `dmar_probability_of_superiority_paired` and
`dmar_tbl`) with columns `term` and `value` and these rows:

- `probability_of_superiority`:

  The estimate \\\hat P_S\\, the proportion of pairs with \\y_i \> x_i\\
  plus half the proportion of tied pairs.

- `ci_lower`, `ci_upper`:

  The confidence limits at `conf_level`. Both equal the estimate when
  the estimated variance is 0, and are then not a confidence interval.

- `var_ps`:

  The estimated variance of \\\hat P_S\\.

- `wins_y_over_x`:

  The number of complete pairs with \\y_i \> x_i\\.

- `losses_y_under_x`:

  The number with \\y_i \< x_i\\.

- `ties`:

  The number with \\y_i = x_i\\.

The table carries the attributes `conf_level`; `n_missing`, the number
of pairs removed because a member was missing; and
`interval_degenerate`, present only when the estimated variance is 0,
with the value `"complete_separation"` or `"no_variation"`. The `tidy()`
and `glance()` verbs (see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md))
add a logical column `interval_degenerate`.

## Details

**Estimand and interval.** For a pair \\(X, Y)\\ drawn from the
population of pairs, \$\$P_S \\=\\ \Pr(Y \> X) + 0.5 \cdot \Pr(Y =
X).\$\$ The sampling unit is the pair, and the \\n\\ pairs are taken to
be independent and identically distributed; nothing is assumed about the
distribution of the scores or about the dependence within a pair. The
interval is a large-sample confidence interval for \\P_S\\, formed on a
transformed scale that keeps it inside \\\[0, 1\]\\. A value of 0.5 says
that neither member tends to score higher.

**Sample estimate.** With the score of pair \\i\\ \$\$s_i =
\mathrm{I}(y_i \> x_i) + 0.5 \cdot \mathrm{I}(y_i = x_i),\$\$ the
estimate is the mean \\\hat P_S = \bar s\\: the proportion of pairs won
by \\y\\, plus half the proportion tied. Being the mean of \\n\\
independent scores with expectation \\P_S\\, it is unbiased. Each pair
is compared directly, not through the difference \\y_i - x_i\\, so a
pair of equal infinite values is a tie.

**Why a function for pairs.** \\P_S\\ compares the two members of one
pair. Vargha and Delaney's *A* and Cliff's \\\delta\\ compare an
observation of one group with an observation of an unrelated unit in the
other, which is a different quantity when the members of a pair are
correlated, and the variance formulas of
[`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md)
and
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
assume independent groups. Paired data therefore need the estimate and
the variance given here.

**Variance and confidence interval.** The variance of the mean of \\n\\
independent scores is the variance of a score divided by \\n\\, and the
function estimates it by \$\$\widehat{\mathrm{Var}}(\hat P_S) \\=\\
\frac{1}{n^2}\sum\_{i=1}^{n} (s_i - \bar s)^2,\$\$ which uses the
divisor \\n\\ in the variance of the scores (it is \\(n - 1)/n\\ times
the estimate with the divisor \\n - 1\\). This is the row `var_ps`. The
interval is formed for \\d = 2 P_S - 1\\, which lies in \\\[-1, 1\]\\,
on the scale of \\\mathrm{arctanh}(d)\\, as the interval of
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
is: \$\$\tanh\\\left\\\mathrm{arctanh}(\hat d) \pm z\_{1 - \alpha/2}\\
\frac{2\sqrt{\widehat{\mathrm{Var}}(\hat P_S)}}{1 - \hat d^2}
\right\\,\$\$ and the two limits \\d_L\\ and \\d_U\\ are returned as
\\(d_L + 1)/2\\ and \\(d_U + 1)/2\\.

**Scale and missing values.** Only the order within each pair is used,
so any increasing rescoring applied to both members gives the same
result, and an ordinal response may be supplied as two ordered factors
with the same levels. A pair with a missing member (`NA` or `NaN`) is
removed, and the number removed is in the `n_missing` attribute. `x` and
`y` must have the same length. Fewer than three complete pairs is an
error, since no variance can be estimated.

**What is refused, and why.** The statistic uses the order of the
observations and nothing else, so the function must be given objects
whose order is the order meant. It refuses, with a condition of class
`dmar_scores_not_ordered` (see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)),
an unordered factor (its codes follow the order of
[`levels()`](https://rdrr.io/r/base/levels.html), alphabetical unless
set, so the levels low, mid, high are stored as 2, 3, 1); a character
vector (strings sort alphabetically, `"10"` before `"9"`); an ordered
factor beside anything that is not an ordered factor with the same
levels in the same order (each factor numbers its own levels from 1);
time differences in different units, or one beside plain numbers (each
such object is stored in its own units); a `Date` beside a `POSIXct`;
and a list, data frame, or matrix of several columns. Nothing is
computed in any of these cases. The message says which case applies and
what to supply: an ordered factor with its levels from lowest to
highest, `factor(x, levels = c("low", "mid", "high"), ordered = TRUE)`,
or numeric scores.

**When the estimated variance is zero.** Two kinds of data give a
variance estimate of exactly 0, and they are the only two. *Complete
separation*: every pair favors the same member, so the estimate is at
the end of its range (0 or 1). *No variation*: every pair is tied, so
the estimate is 0.5. In both, the two limits are reported equal to the
estimate, `var_ps` is 0, the table's `interval_degenerate` attribute
names the case, and a warning of class `dmar_interval_degenerate` is
given whose `reason` element is `"complete_separation"` or
`"no_variation"`; the condition also carries `estimate` and `n` (the
number of pairs), and its `limit_feng_cliff` is `NA`, since that limit
is for two independent groups. The reported limits are not a confidence
interval with the stated coverage. The zero describes the sample: the
variance estimate measures how much the scores differ from one pair to
the next, and when every pair favors the same member they do not differ.
Complete separation is common when the effect is large and the pairs are
few. Without ties, \\n\\ pairs all favor the same member with
probability \\P_S^n + (1 - P_S)^n\\: 0.47 for five pairs and 0.22 for
ten when \\P_S = 0.86\\, and 0.82 and 0.66 when \\P_S = 0.96\\. In those
samples the reported limits, both 0 or both 1, do not contain a
population value inside (0, 1).

Particular cases. Infinite values are ordered like any others, and a
pair of equal infinite values is a tie. Wins together with ties, and no
losses, give a positive variance and an ordinary interval (three wins
and a tie give 0.875 with limits 0.502 and 0.980). The mark travels with
the result:
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
writes the estimate and says that there is no confidence interval, the
printed table says so beneath itself, and `tidy()` and `glance()` carry
it in the column `interval_degenerate`. The randomization test for
paired observations
([`randomization_test_paired`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md))
remains valid in both cases.

**Reading the output.** In the example below 20 of the 30 pairs have a
higher `post` than `pre` score and none is tied, so \\\hat P_S = 20/30 =
0.667\\, with 95 percent limits 0.484 and 0.810.

## References

Grissom, R. J., & Kim, J. J. (2005). *Effect sizes for research: A broad
practical approach*. Lawrence Erlbaum.

Grissom, R. J., & Kim, J. J. (2012). *Effect sizes for research:
Univariate and multivariate applications* (2nd ed.). Routledge.

Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the
CL common language effect size statistics of McGraw and Wong. *Journal
of Educational and Behavioral Statistics, 25*(2), 101–132.
[doi:10.3102/10769986025002101](https://doi.org/10.3102/10769986025002101)

## See also

[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md),
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`randomization_test_paired`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`proportion_of_superiority`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md)

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
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Paired pre/post data:
set.seed(113)
pre  <- rnorm(30, mean = 100, sd = 15)
post <- pre + rnorm(30, mean =  5, sd = 10)
probability_of_superiority_paired(x = pre, y = post)
#>  term                       value
#>  probability_of_superiority 0.667
#>  ci_lower                   0.484
#>  ci_upper                   0.810
#>  var_ps                     0.007
#>  wins_y_over_x              20   
#>  losses_y_under_x           10   
#>  ties                       0    
#> 
#> Confidence level: 95%

# 2. Every pair favors the same member: the estimated variance is 0,
#    the two limits equal the estimate, and a warning says that they
#    are not a confidence interval.
all_up <- probability_of_superiority_paired(x = c(1, 2, 3, 4),
                                            y = c(4, 5, 6, 7))
#> Warning: probability_of_superiority_paired(): every pair favors the same member (probability of superiority = 1), so the estimated variance is 0. Both confidence limits are reported equal to the estimate (1, 1). They are not a confidence interval with the stated coverage: the variance estimate is built from how much the comparisons differ from one pair to the next, and in a sample in which every pair favors the same member they do not differ at all. That describes this sample, not the population. A population in which some pairs favor the other member often produces a sample with none when there are few pairs (here 4 pairs), so the limits must not be read as showing that the population value is 1. Report the estimate with the sample size and say that every pair favored the same member; the randomization test (randomization_test_paired()) gives a p-value that remains valid here.
all_up
#>  term                       value
#>  probability_of_superiority 1.000
#>  ci_lower                   1.000
#>  ci_upper                   1.000
#>  var_ps                     0.000
#>  wins_y_over_x              4    
#>  losses_y_under_x           0    
#>  ties                       0    
#> 
#> Confidence level: 95%
#> 
#> Limits are not a confidence interval: the estimated variance is 0.
results_sentence(all_up)
#> [1] "probability_of_superiority = 1.00 (no confidence interval: every pair favors the same member, so the estimated variance is 0)"
```
