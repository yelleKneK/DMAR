# Randomization (Permutation) Test for Two Independent Groups

Compares two independent groups by referring the observed statistic to
the distribution of that same statistic over reassignments of the
observed scores to the two groups. That reference distribution, not a
normal or a *t* distribution, supplies the *p*-value, so the test needs
no assumption about the shape of the population. Alongside the test, the
function reports the effect sizes that answer the question the test only
screens: the mean difference and its randomization-based interval, the
standardized mean difference with a noncentral *t* interval, the common
language effect size, and Cliff's delta.

## Usage

``` r
randomization_test(
  x = NULL,
  group = NULL,
  data = NULL,
  group_1 = NULL,
  group_2 = NULL,
  statistic = c("mean", "t"),
  alternative = c("two_sided", "less", "greater"),
  exact = NULL,
  n_resamples = 10000L,
  seed = NULL,
  conf_level = 0.95,
  shift_ci = TRUE
)
```

## Arguments

- x:

  A formula of the form `y ~ group`, or a numeric response vector to be
  paired with `group`. Leave `NULL` when the two samples are supplied
  through `group_1` and `group_2`.

- group:

  A grouping variable the same length as `x`, with exactly two levels
  after unused levels are dropped. Ignored when `x` is a formula.

- data:

  An optional `data.frame` in which to find the variables named in the
  formula.

- group_1, group_2:

  The two samples supplied directly as numeric vectors, an alternative
  to the formula and response-plus-grouping interfaces. Lengths need not
  be equal.

- statistic:

  One of `"mean"` (default) or `"t"`. `"mean"` uses the difference in
  group means. `"t"` uses the studentized (Welch) statistic, which
  divides that difference by its separate-variances standard error and
  is the better choice when the groups may differ in variance (see
  Details).

- alternative:

  One of `"two_sided"` (default; the base-R spelling `"two.sided"` is
  accepted as an alias), `"less"`, or `"greater"`. The direction refers
  to the first group minus the second, the same orientation
  [`t.test`](https://rdrr.io/r/stats/t.test.html) uses.

- exact:

  Logical. If `NULL` (default), the test enumerates every reassignment
  when `choose(N, n_1)` is at most 50,000 and samples reassignments
  otherwise. If `TRUE`, enumeration is forced (refused above 1,000,000
  reassignments). If `FALSE`, Monte Carlo is forced.

- n_resamples:

  Number of randomly drawn reassignments when enumeration is not used.
  Default `10000L`.

- seed:

  Optional integer seed for the Monte Carlo branch. Default `NULL`,
  which leaves the user's current RNG state intact; supply an integer
  for reproducibility. When a seed is supplied the RNG state in place
  before the call is restored on exit.

- conf_level:

  Confidence level for every interval reported, the inverted
  randomization interval included. Default `0.95`.

- shift_ci:

  Logical. Compute the randomization-based interval for the shift by
  inverting the test? Default `TRUE`. Setting it to `FALSE` reports the
  two endpoints as `NA` and skips the inversion, which is the expensive
  part of the call.

## Value

A `data.frame` with a `term` column and a numeric `value` column, in
three blocks.

The test: `mean_difference` (first group minus second), `statistic` (the
statistic actually referred to the reference distribution), `p_value`,
and `p_value_se` (the Monte Carlo standard error of the *p*-value, `NA`
under exact enumeration, which has no Monte Carlo error).

The intervals and effect sizes: `shift_lower_limit` and
`shift_upper_limit` (the randomization interval for the shift, obtained
by inverting the test); `normal_theory_lower_limit` and
`normal_theory_upper_limit` (Welch's *t* interval on the same mean
difference, reported for contrast); `smd` with `smd_lower_limit` and
`smd_upper_limit` from
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md); `cles`
with `cles_lower_limit` and `cles_upper_limit` from
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md); and
`cliff_delta` with `cliff_delta_lower_limit` and
`cliff_delta_upper_limit` from
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md).

The design: `n_1`, `n_2`, `N`, `n_evaluated` (how many reassignments
were actually used), and `exact` (1 if every reassignment was
enumerated, 0 if they were sampled).

Non-numeric information travels on attributes rather than in the `value`
column: `statistic_name`, `method` (“exact enumeration” or “Monte
Carlo”), `alternative`, `group_labels`, `response_name`, `group_name`,
`seed`, `observed_statistic`, and `reference_distribution`, the vector
of statistics over the reassignments that
[`plot_randomization_test`](https://yelleknek.github.io/DMAR/reference/plot_randomization_test.md)
draws.

## Details

**What the randomization distribution is.** Suppose *N* participants
were randomly assigned, \\n_1\\ to one condition and \\n_2\\ to the
other. Under the null hypothesis that the condition a participant
received made no difference to that participant's score, each score
would have been the same number no matter which group the participant
landed in. The assignment actually used was one draw from the
\\\binom{N}{n_1}\\ assignments the randomization could equally well have
produced, so every one of those assignments was equally likely, and each
of them yields a value of the test statistic. Those values are the
randomization distribution. The *p*-value is the proportion of them at
least as extreme as the value the experiment actually produced.

**Why there is no normality assumption.** Nothing in that argument
mentions a population, a normal curve, or a sampling model. The
probability comes from the coin flips the experimenter performed, which
are known exactly because the experimenter performed them. This is the
inferential logic Fisher (1935) used to introduce experimental design,
and it is where Chapter 1 of Maxwell, Delaney, and Kelley (2027) starts,
for the same reason: the validity of the test rests on the randomization
rather than on assumptions a data analyst cannot check.

**What the test does and does not license.** A small *p*-value says the
observed separation between the groups would rarely arise from
reassignment alone, which is evidence that the assignment mattered. It
does not say how much it mattered, and with a large *N* an uninteresting
difference will produce a small *p*-value. It also does not, by itself,
license generalization beyond the participants at hand: randomization
licenses a causal claim about these units, while generalization to a
population is a separate argument that rests on how the units were
recruited. That is why this function reports effect sizes with intervals
rather than a *p*-value alone.

**Why the studentized statistic.** With \\n_1 = n_2\\ and equal
population variances the two statistics give the same *p*-value to
within the discreteness of the reference distribution, because the
denominator of the studentized statistic is then nearly constant across
reassignments. When the variances differ and the groups are unbalanced
they part company. Reassigning scores between groups of unequal size
mixes the two variances in proportions that the observed assignment does
not have, so the reference distribution for the raw mean difference is
built under a null that is false in a second way, and the test's actual
Type I error rate drifts away from the nominal level. The studentized
statistic rescales each reassignment by its own separate-variances
standard error, which removes most of that drift and remains
asymptotically valid under heteroscedasticity (Janssen, 1997; Neuhaus,
1993). Use `statistic = "t"` whenever unequal variances are plausible,
which for unbalanced designs is nearly always.

**Exact or Monte Carlo.** When `choose(N, n_1)` is at most 50,000 every
reassignment is enumerated and the *p*-value is exact: it is a count
divided by a known total, with no approximation anywhere. Above that
threshold `n_resamples` reassignments are drawn at random and the
*p*-value is \\(r + 1) / (m + 1)\\, where *r* counts the sampled
reassignments at least as extreme as the observed one and *m* is
`n_resamples`. Adding one to each part counts the observed assignment,
which is itself a legitimate reassignment; without it a *p*-value of
exactly zero could be reported for a hypothesis the data cannot rule
out, and the test would be anticonservative (Phipson & Smyth, 2010). The
reported `p_value_se` is \\\sqrt{\hat p (1 - \hat p) / m}\\, the
standard error of the resampling itself. It describes how much the
*p*-value would move if the reassignments were drawn again, not how much
it would move in a new experiment. Raising `n_resamples` shrinks it at
the usual \\1/\sqrt{m}\\ rate.

**The randomization interval, and how it differs from the normal theory
one.** Suppose the treatment adds a constant \\\delta\\ to every score
it touches. Subtracting \\\delta\\ from each first-group score should
then leave scores that are exchangeable across groups, so the
randomization test applied to the subtracted data is a test of \\H_0\\:
\mathrm{shift} = \delta\\. The set of \\\delta\\ for which that test
does not reject at level \\1 - \\ `conf_level` is a confidence interval
for the shift, and it is reported as `shift_lower_limit` and
`shift_upper_limit`. Inverting a test this way is the general recipe
(Ernst, 2004); the endpoints are located by bisection on the *p*-value,
using the same reassignments throughout so the interval and the test
agree.

The contrast with `normal_theory_lower_limit` and
`normal_theory_upper_limit`, which are Welch's *t* limits on the same
mean difference, is worth reading whenever both are printed. The
randomization interval is exactly the set of shifts the test being run
does not reject, so the test and the interval can never disagree. The
normal theory interval instead assumes the sampling distribution of the
mean difference has a known shape; it is smooth, symmetric about the
point estimate, and can extend past the range the data can support. The
randomization interval is discrete, need not be symmetric, and in a very
small design is unbounded, a correct statement of how little information
the design carries rather than a defect: with three observations per
group the smallest attainable two-sided *p*-value is 2/20 = 0.10, so no
shift can be rejected at the 5% level and the 95% interval is the whole
real line. The randomization interval also inherits the shift model, so
it answers a narrower question than the test does: the test needs only
exchangeability, while the interval needs the treatment to move every
score by the same amount.

**Effect sizes.** Every effect size reported here comes from the package
function that owns it, so the numbers match a direct call.
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) and
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) supply
the standardized mean difference and its noncentral *t* interval;
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md) supplies
the common language effect size, the probability that a randomly drawn
score from the first group exceeds one from the second, by transforming
those limits through \\\Phi(\cdot/\sqrt 2)\\; and
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
supplies Cliff's delta with its consistent interval. The standardized
mean difference and the common language effect size are normal theory
quantities, so their intervals lean on the assumption the test itself
avoids. Cliff's delta does not: it is a function of the ordering of the
observations alone, which makes it the natural effect size companion to
a randomization test. Reporting all three lets a reader see whether the
distribution-free and normal theory summaries tell the same story.

**Ranks give the Wilcoxon test.** Replacing the scores by their ranks
and running this test with `statistic = "mean"` reproduces the exact
Wilcoxon rank sum test, since the rank sum is a monotone function of the
difference in mean ranks. That equivalence is a useful check and a
reminder of what the rank test is: a randomization test on transformed
data.

## References

Edgington, E. S., & Onghena, P. (2007). *Randomization tests* (4th ed.).
Chapman & Hall/CRC.

Ernst, M. D. (2004). Permutation methods: A basis for exact inference.
*Statistical Science, 19*(4), 676–685.
[doi:10.1214/088342304000000396](https://doi.org/10.1214/088342304000000396)

Fisher, R. A. (1935). *The design of experiments*. Oliver & Boyd.

Janssen, A. (1997). Studentized permutation tests for non-i.i.d.
hypotheses and the generalized Behrens-Fisher problem. *Statistics &
Probability Letters, 36*(1), 9–21.
[doi:10.1016/S0167-7152(97)00043-6](https://doi.org/10.1016/S0167-7152%2897%2900043-6)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 1 on the logic of randomization and the
randomization test.)

Neuhaus, G. (1993). Conditional rank tests for the two-sample problem
under random censorship. *The Annals of Statistics, 21*(4), 1760–1779.
[doi:10.1214/aos/1176349396](https://doi.org/10.1214/aos/1176349396)

Phipson, B., & Smyth, G. K. (2010). Permutation *p*-values should never
be zero: Calculating exact *p*-values when permutations are randomly
drawn. *Statistical Applications in Genetics and Molecular Biology,
9*(1), Article 39.
[doi:10.2202/1544-6115.1585](https://doi.org/10.2202/1544-6115.1585)

Pitman, E. J. G. (1937). Significance tests which may be applied to
samples from any populations. *Supplement to the Journal of the Royal
Statistical Society, 4*(1), 119–130.

## See also

[`plot_randomization_test`](https://yelleknek.github.io/DMAR/reference/plot_randomization_test.md)
for the figure that shows the reference distribution, the observed
statistic, and the rejection region;
[`randomization_test_paired`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md)
for the sign-flip sibling used with paired observations;
[`t.test`](https://rdrr.io/r/stats/t.test.html) and
[`wilcox.test`](https://rdrr.io/r/stats/wilcox.test.html) for the
parametric and rank-based alternatives;
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`cles`](https://yelleknek.github.io/DMAR/reference/cles.md), and
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
for the effect sizes reported here.

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
[`obrien_test()`](https://yelleknek.github.io/DMAR/reference/obrien_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Ten observations, so every one of the choose(10, 5) = 252
#    reassignments is enumerated and the p-value is exact.
treatment <- c(80, 84, 79, 88, 83)
control   <- c(72, 75, 68, 81, 74)
randomization_test(group_1 = treatment, group_2 = control)
#>  term                      value 
#>  mean_difference           8.8   
#>  statistic                 8.8   
#>  p_value                   0.0238
#>  p_value_se                <NA>  
#>  shift_lower_limit         3     
#>  shift_upper_limit         15    
#>  normal_theory_lower_limit 2.6   
#>  normal_theory_upper_limit 15    
#>  smd                       2.1   
#>  smd_lower_limit           0.461 
#>  smd_upper_limit           3.66  
#>  cles                      0.931 
#>  cles_lower_limit          0.628 
#>  cles_upper_limit          0.995 
#>  cliff_delta               0.84  
#>  cliff_delta_lower_limit   -0.117
#>  cliff_delta_upper_limit   0.988 
#>  n_1                       5     
#>  n_2                       5     
#>  N                         10    
#>  n_evaluated               252   
#>  exact                     1     
#> 
#> Confidence level: 95%

# 2. The studentized statistic, preferable when the groups may differ
#    in variance.
randomization_test(group_1 = treatment, group_2 = control,
                   statistic = "t")
#>  term                      value 
#>  mean_difference           8.8   
#>  statistic                 3.32  
#>  p_value                   0.0238
#>  p_value_se                <NA>  
#>  shift_lower_limit         3     
#>  shift_upper_limit         15    
#>  normal_theory_lower_limit 2.6   
#>  normal_theory_upper_limit 15    
#>  smd                       2.1   
#>  smd_lower_limit           0.461 
#>  smd_upper_limit           3.66  
#>  cles                      0.931 
#>  cles_lower_limit          0.628 
#>  cles_upper_limit          0.995 
#>  cliff_delta               0.84  
#>  cliff_delta_lower_limit   -0.117
#>  cliff_delta_upper_limit   0.988 
#>  n_1                       5     
#>  n_2                       5     
#>  N                         10    
#>  n_evaluated               252   
#>  exact                     1     
#> 
#> Confidence level: 95%

# 3. Formula interface: weekly drinking in the two comparable arms of
#    the drinks_trial data, a right-skewed outcome, which is exactly
#    where a distribution-free test earns its keep. With 37 and 32
#    participants there are far too many reassignments to enumerate, so
#    10,000 are drawn and the p-value carries a Monte Carlo standard
#    error.
cra <- droplevels(subset(drinks_trial, treatment != "CRA + Disulfiram"))
set.seed(113)
randomization_test(drinks_per_week ~ treatment, data = cra, seed = 113)
#>  term                      value  
#>  mean_difference           24.5   
#>  statistic                 24.5   
#>  p_value                   0.2784 
#>  p_value_se                0.00448
#>  shift_lower_limit         -13    
#>  shift_upper_limit         63.1   
#>  normal_theory_lower_limit -15.8  
#>  normal_theory_upper_limit 64.7   
#>  smd                       0.282  
#>  smd_lower_limit           -0.195 
#>  smd_upper_limit           0.756  
#>  cles                      0.579  
#>  cles_lower_limit          0.445  
#>  cles_upper_limit          0.704  
#>  cliff_delta               0.3    
#>  cliff_delta_lower_limit   0.0213 
#>  cliff_delta_upper_limit   0.535  
#>  n_1                       37     
#>  n_2                       32     
#>  N                         69     
#>  n_evaluated               10000  
#>  exact                     0      
#> 
#> Confidence level: 95%

# 4. A one-sided test, and the one-sided interval that goes with it.
randomization_test(group_1 = treatment, group_2 = control,
                   alternative = "greater")
#>  term                      value 
#>  mean_difference           8.8   
#>  statistic                 8.8   
#>  p_value                   0.0119
#>  p_value_se                <NA>  
#>  shift_lower_limit         4     
#>  shift_upper_limit         Inf   
#>  normal_theory_lower_limit 3.82  
#>  normal_theory_upper_limit Inf   
#>  smd                       2.1   
#>  smd_lower_limit           0.461 
#>  smd_upper_limit           3.66  
#>  cles                      0.931 
#>  cles_lower_limit          0.628 
#>  cles_upper_limit          0.995 
#>  cliff_delta               0.84  
#>  cliff_delta_lower_limit   -0.117
#>  cliff_delta_upper_limit   0.988 
#>  n_1                       5     
#>  n_2                       5     
#>  N                         10    
#>  n_evaluated               252   
#>  exact                     1     
#> 
#> Confidence level: 95%

# 5. On ranks, the test is the exact Wilcoxon rank sum test.
y <- c(treatment, control)
g <- rep(c("treatment", "control"), each = 5)
res <- randomization_test(rank(y), g)
res$value[res$term == "p_value"]
#> [1] 0.03174603
wilcox.test(treatment, control, exact = TRUE)$p.value
#> [1] 0.03174603
```
