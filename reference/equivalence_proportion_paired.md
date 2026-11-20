# Equivalence and Noninferiority Tests for Two Paired Proportions

Tests whether the difference between two proportions measured on the
same units, \\\psi = p_x - p_y\\, lies inside equivalence bounds or
above a noninferiority margin: whether a cheaper diagnostic procedure
finds a condition about as often as the standard one on the same
patients, or whether a short form of a scale classifies respondents
about as often as the full form. The tests are the score tests of Nam
(1997), which take the place of the *t* tests of
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
and the decisions are read from the score interval of
[`ci_proportion_paired`](https://yelleknek.github.io/DMAR/reference/ci_proportion_paired.md)
(Tango, 1998), with which they agree.

## Usage

``` r
equivalence_proportion_paired(
  x = NULL,
  y = NULL,
  n_10 = NULL,
  n_01 = NULL,
  N = NULL,
  prop_diff_lower = NULL,
  prop_diff_upper = NULL,
  alpha_level = 0.05
)
```

## Arguments

- x, y:

  Numeric or logical vectors of the same length, one element per unit,
  coded 1 (success) and 0 (failure): the outcome of each unit under the
  first condition (`x`) and the second (`y`). A pair with a missing
  value in either vector is dropped. Supply either `x` and `y` or the
  three counts, not both.

- n_10:

  Number of units that succeed under `x` and fail under `y`.

- n_01:

  Number of units that fail under `x` and succeed under `y`.

- N:

  Number of units (pairs), at least `n_10 + n_01`.

- prop_diff_lower, prop_diff_upper:

  Bounds on the scale of the difference in proportions, each a single
  number in (0, 1); the equivalence region is \\(-\delta_L, \delta_U)\\,
  with `prop_diff_lower` as \\\delta_L\\ and `prop_diff_upper` as
  \\\delta_U\\. Supply at least one. With `prop_diff_upper` alone the
  bounds are symmetric. With `prop_diff_lower` alone the call is a
  noninferiority test with margin \\\delta_L\\, and the rows that would
  need an upper bound are `NA` (see Details).

- alpha_level:

  One-sided significance level of each test. Defaults to 0.05, so the
  decisions are read from the 90% score interval.

## Value

A `data.frame` (class `dmar_tbl`) with rows `difference` (\\\hat\psi\\);
the counts `n_10`, `n_01`, and `N`; the score statistics `z_lower`
(\\Z(-\delta_L)\\) and `z_upper` (\\Z(\delta_U)\\) and their one-sided
*p*-values `p_lower` and `p_upper`; the TOST *p*-value `p_tost`; the
noninferiority *p*-value `p_noninferiority` (equal to `p_lower` by
construction); the 100(1 - 2\\\alpha\\)% score limits `ci_lower` and
`ci_upper`; the bounds `prop_diff_lower` (stored signed, as
\\-\delta_L\\) and `prop_diff_upper`; and the flags `equivalent`,
`noninferior`, `superior`, `inferior`, and `nonsuperior` (1 = declared,
0 = not). A noninferiority test alone, with no `prop_diff_upper`, has
the same rows, with `NA` in those that would need an upper bound. The
attributes are `conf_level` (\\1 - 2\alpha\\), the `verdict`, and
`p_sentence`, which names the *p*-value
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
reports beside the interval: `p_tost`, or `p_noninferiority` for a
noninferiority test alone.

## Details

**The tests.** With the score statistic \\Z(\psi_0)\\ described on
[`ci_proportion_paired`](https://yelleknek.github.io/DMAR/reference/ci_proportion_paired.md),
the noninferiority null hypothesis \\H_0: \psi \le -\delta_L\\ is
rejected in favor of \\\psi \> -\delta_L\\ when \\Z(-\delta_L) \>
z\_{1-\alpha}\\, and the nonsuperiority null hypothesis \\H_0: \psi \ge
\delta_U\\ is rejected in favor of \\\psi \< \delta_U\\ when
\\Z(\delta_U) \< -z\_{1-\alpha}\\. Equivalence, \\-\delta_L \< \psi \<
\delta_U\\, is claimed when both are rejected. Because the claim
requires both rejections, the two one-sided tests form an
intersection-union test (Berger & Hsu, 1996): each is carried out at
level \\\alpha\\ and the equivalence claim has level \\\alpha\\, with no
adjustment for there being two tests. The *p*-value of the equivalence
test, `p_tost`, is the larger of the two one-sided *p*-values.

**Why the score statistic and not McNemar's.** McNemar's test estimates
the variance of the per-unit difference as it would be if \\\psi\\ were
0, so it tests no difference and nothing else. A test against a margin
needs the variance at \\\psi = -\delta_L\\ (or \\\delta_U\\), which
depends on a discordant probability that the hypothesis does not fix;
the score statistic replaces that probability with its maximum
likelihood estimate under the hypothesis (Nam, 1997). At \\\psi_0 = 0\\
the score statistic is McNemar's, so the tests here generalize it. Nam
(1997) found that a Wald-type statistic, with the variance estimated
from the observed proportions instead, rejects more often than its
nominal level, while the score test holds close to it.

**Reading the interval.** \\Z\\ decreases in \\\psi_0\\, and the score
interval is the set of differences the score test does not reject. So
the noninferiority test rejects at level \\\alpha\\ exactly when the
lower limit of the 100(1 - 2\\\alpha\\)% score interval exceeds
\\-\delta_L\\, and the nonsuperiority test exactly when its upper limit
is below \\\delta_U\\. The flags are computed from that interval, which
the table reports.

**The flags and the verdict.** The flags read the interval against the
bounds:

- `equivalent`: the interval lies inside \\(-\delta_L, \delta_U)\\.

- `noninferior`: the lower limit exceeds \\-\delta_L\\.

- `superior`: the lower limit exceeds \\\delta_U\\.

- `inferior`: the upper limit is below \\-\delta_L\\.

- `nonsuperior`: the upper limit is below \\\delta_U\\, so a difference
  in favor of `x` as large as \\\delta_U\\ is ruled out. It is the
  mirror image of noninferiority.

With both bounds, the `"verdict"` attribute names the most informative
statement the interval supports: `"Equivalent"`, `"Superior"`,
`"Inferior"`, `"Noninferior only"` (an interval above \\-\delta_L\\ that
contains \\\delta_U\\), `"Not relevantly superior only"` (an interval
below \\\delta_U\\ that contains \\-\delta_L\\), or `"Inconclusive"` (an
interval that contains both bounds). An inconclusive result is a
statement about precision, not evidence of a difference.

**Noninferiority alone.** A call with `prop_diff_lower` and no
`prop_diff_upper` asks only whether \\\psi \> -\delta_L\\, and no upper
bound is assumed for it. The rows that would need one, `z_upper`,
`p_upper`, `p_tost`, `prop_diff_upper`, `equivalent`, `superior`, and
`nonsuperior`, are `NA`, so the table keeps its shape; the verdict is
`"Noninferior"` (the lower limit exceeds \\-\delta_L\\), `"Inferior"`
(the upper limit is below \\-\delta_L\\), or `"Inconclusive"`; and the
*p*-value a results sentence reports is `p_noninferiority`. This is the
noninferiority test of
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
which reads \\-\delta_L\\ alone. A call with `prop_diff_upper` alone
gives the symmetric bounds \\(-\delta_U, \delta_U)\\, as
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md)
does.

**Choosing the bounds.** The bounds are fixed before the data are
examined, on substantive grounds, as the smallest difference in the two
proportions that would matter, and on the scale of the difference
itself: `prop_diff_lower = 0.05` is five percentage points.

## References

Berger, R. L., & Hsu, J. C. (1996). Bioequivalence trials,
intersection-union tests and equivalence confidence sets. *Statistical
Science, 11*(4), 283–319.
[doi:10.1214/ss/1032280304](https://doi.org/10.1214/ss/1032280304)

McNeil, B. J., Keeler, E., & Adelstein, S. J. (1975). Primer on certain
elements of medical decision making. *New England Journal of Medicine,
293*(5), 211–215.
[doi:10.1056/NEJM197507312930501](https://doi.org/10.1056/NEJM197507312930501)

Nam, J.-m. (1997). Establishing equivalence of two treatments and sample
size requirements in matched-pairs design. *Biometrics, 53*(4),
1422–1430. [doi:10.2307/2533508](https://doi.org/10.2307/2533508)

Tango, T. (1998). Equivalence test and confidence interval for the
difference in proportions for the paired-sample design. *Statistics in
Medicine, 17*(8), 891–908.
[doi:10.1002/(SICI)1097-0258(19980430)17:8\<891::AID-SIM780\>3.0.CO;2-B](https://doi.org/10.1002/%28SICI%291097-0258%2819980430%2917%3A8%3C891%3A%3AAID-SIM780%3E3.0.CO%3B2-B)

## See also

[`ci_proportion_paired`](https://yelleknek.github.io/DMAR/reference/ci_proportion_paired.md)
for the interval alone;
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md)
for a linear contrast of means.

Other equivalence testing:
[`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`power_density_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md),
[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md)

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_dunnett_t3()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett_t3.md),
[`ci_games_howell()`](https://yelleknek.github.io/DMAR/reference/ci_games_howell.md),
[`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
[`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`compare_cov_structures()`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`factorial_anova()`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`obrien_test()`](https://yelleknek.github.io/DMAR/reference/obrien_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Liver scan data from McNeil, Keeler, and Adelstein (1975), the
# example of Nam (1997, p. 1427): 344 patients classified by scan (x)
# and by pathology (y), 32 positive by scan only and 27 by pathology
# only. Does the scan understate the rate of disease by less than five
# percentage points? That is a noninferiority question, so only
# prop_diff_lower is given, and the rows that would need an upper bound
# are NA. z_lower is the score statistic Nam reports, 2.795; the table
# rounds it for display and stores the full value.
res <- equivalence_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
                                     prop_diff_lower = 0.05)
res
#>  term             value 
#>  difference       0.015 
#>  n_10             32    
#>  n_01             27    
#>  N                344   
#>  z_lower          2.795 
#>  z_upper          <NA>  
#>  p_lower          0.0026
#>  p_upper          <NA>  
#>  p_tost           <NA>  
#>  p_noninferiority 0.0026
#>  ci_lower         -0.023
#>  ci_upper         0.052 
#>  prop_diff_lower  -0.050
#>  prop_diff_upper  <NA>  
#>  equivalent       <NA>  
#>  noninferior      1     
#>  superior         <NA>  
#>  inferior         0     
#>  nonsuperior      <NA>  
#> 
#> Confidence level: 90%
res$value[res$term == "z_lower"]
#> [1] 2.795403
attr(res, "verdict")
#> [1] "Noninferior"

# Agreement within five points in both directions is a stronger
# claim. The interval rules out a shortfall of five points but not an
# excess of five, so the verdict is noninferiority only.
res <- equivalence_proportion_paired(n_10 = 32, n_01 = 27, N = 344,
                                     prop_diff_lower = 0.05,
                                     prop_diff_upper = 0.05)
attr(res, "verdict")
#> [1] "Noninferior only"

# A short form (x) and the full form (y) of a screening scale given to
# the same 400 respondents: 38 screen positive on the short form only
# and 52 on the full form only. A shortfall of five points is not ruled
# out, but an excess of five is.
res <- equivalence_proportion_paired(n_10 = 38, n_01 = 52, N = 400,
                                     prop_diff_lower = 0.05,
                                     prop_diff_upper = 0.05)
res
#>  term             value 
#>  difference       -0.035
#>  n_10             38    
#>  n_01             52    
#>  N                400   
#>  z_lower          0.632 
#>  z_upper          -3.498
#>  p_lower          0.2636
#>  p_upper          0.0002
#>  p_tost           0.2636
#>  p_noninferiority 0.2636
#>  ci_lower         -0.074
#>  ci_upper         0.004 
#>  prop_diff_lower  -0.050
#>  prop_diff_upper  0.050 
#>  equivalent       0     
#>  noninferior      0     
#>  superior         0     
#>  inferior         0     
#>  nonsuperior      1     
#> 
#> Confidence level: 90%
attr(res, "verdict")
#> [1] "Not relevantly superior only"
```
