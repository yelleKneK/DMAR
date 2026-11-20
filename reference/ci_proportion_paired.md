# Confidence Interval for the Difference Between Two Paired Proportions

Computes a confidence interval for the difference between two
proportions measured on the same units, \\\psi = p_x - p_y\\: two
diagnostic procedures applied to the same patients, the proportion of
students answering an item correctly before and after instruction, or
two raters' rates of assigning a category to the same cases. The two
outcomes of a unit are correlated, so the interval for two independent
proportions does not apply. The default is the score interval of Tango
(1998), which inverts the score test of Nam (1997). It stays inside
\\\[-1, 1\]\\ and is defined when no pair is discordant. When discordant
pairs are rare, as in the case worked through under Details, it covers
the difference more often than its nominal level, erring on the
conservative side, where the familiar Wald interval falls well short of
that level.

## Usage

``` r
ci_proportion_paired(
  x = NULL,
  y = NULL,
  n_10 = NULL,
  n_01 = NULL,
  N = NULL,
  conf_level = 0.95,
  method = c("score", "wald")
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

- conf_level:

  Confidence level. Defaults to 0.95. A `conf_level` of 0.90 gives the
  interval whose limits are the one-sided 95% bounds read in
  noninferiority and equivalence decisions.

- method:

  `"score"` (default; Tango, 1998) or `"wald"`.

## Value

A `data.frame` (class `dmar_tbl`) with rows `ci_lower`, `difference`
(\\\hat\psi\\), `ci_upper`, `n_10`, `n_01`, and `N`, so the point
estimate sits between its confidence limits. The confidence level is the
`conf_level` attribute and the interval's method the `method` attribute.

## Details

**Only the discordant pairs carry the difference.** Classify the \\N\\
units by their pair of outcomes: \\n\_{11}\\ succeed under both
conditions, \\n\_{10}\\ under `x` only, \\n\_{01}\\ under `y` only, and
\\n\_{00}\\ under neither. The two proportions are \\(n\_{11} +
n\_{10})/N\\ and \\(n\_{11} + n\_{01})/N\\, so \\n\_{11}\\ cancels from
their difference, \\\hat\psi = (n\_{10} - n\_{01})/N\\, and the
concordant units enter only through \\N\\. That is why the counts
interface asks for `n_10`, `n_01`, and `N` and nothing more. The
variance of the per-unit difference, \\p\_{10} + p\_{01} - \psi^2\\,
depends on a discordant probability as well as on \\\psi\\, so for
inference on \\\psi\\ the probability \\p\_{01}\\ is a nuisance
parameter.

**The score interval.** For a hypothesized difference \\\psi_0\\ the
score statistic is \$\$Z(\psi_0) = \frac{n\_{10} - n\_{01} - N\psi_0}
{\sqrt{N\\2\tilde p\_{01} + \psi_0(1 - \psi_0)\\}},\$\$ where \\\tilde
p\_{01}\\ is the maximum likelihood estimate of \\p\_{01}\\ under the
restriction \\p\_{10} - p\_{01} = \psi_0\\, the larger root of a
quadratic (Nam, 1997, equations 1 and 3). The denominator is \\N\\ times
the variance of the per-unit difference, evaluated at the restricted
estimates rather than at the observed proportions. With \\\alpha\\ = 1 -
`conf_level`, the lower limit is the value of \\\psi_0\\ at which
\\Z(\psi_0) = z\_{1-\alpha/2}\\ and the upper limit the value at which
\\Z(\psi_0) = -z\_{1-\alpha/2}\\, found by root finding. At \\\psi_0 =
0\\ the statistic is McNemar's, \\(n\_{10} - n\_{01})/\sqrt{n\_{10} +
n\_{01}}\\, so the score interval excludes 0 exactly when McNemar's test
without a continuity correction rejects at level \\\alpha\\: the
interval extends that test to every other value of the difference. With
no discordant pair the limits are \\\mp z^2/(N + z^2)\\, with \\z =
z\_{1-\alpha/2}\\, an interval whose width reflects the \\N\\ units
observed rather than a width of zero.

**The Wald interval.** The interval \\\hat\psi \pm
z\_{1-\alpha/2}\\(n\_{10} + n\_{01})/N - \hat\psi^2\\^{1/2}/\sqrt{N}\\
estimates the variance at the observed proportions instead, and is
available for instruction and comparison. With few discordant pairs it
is too short: it has zero width when no pair is discordant, and it can
extend beyond \\\[-1, 1\]\\. Nam (1997) found by simulation that the
test it inverts rejects more often than its nominal level, while the
score test holds close to the nominal level and tends to the
conservative side. Summing over every possible table with \\p\_{10} =
.03\\, \\p\_{01} = .01\\, and \\N = 50\\, the 95% score interval covers
\\\psi\\ with probability .986 and the 95% Wald interval with
probability .779.

**Noninferiority and equivalence.** The limits of a 100(1 -
2\\\alpha\\)% interval are one-sided 100(1 - \\\alpha\\)% confidence
bounds, so `conf_level = 0.90` gives the interval that
[`equivalence_proportion_paired`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md)
reads against its bounds at `alpha_level = 0.05`.

## References

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

[`equivalence_proportion_paired`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md)
for noninferiority and equivalence tests on the same difference;
[`ci_proportion`](https://yelleknek.github.io/DMAR/reference/ci_proportion.md)
for a single proportion;
[`responder_analysis`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md)
for two independent proportions;
[`mcnemar.test`](https://rdrr.io/r/stats/mcnemar.test.html) for the test
of no difference.

Other confidence intervals:
[`ci_crossed()`](https://yelleknek.github.io/DMAR/reference/ci_crossed.md),
[`ci_proportion()`](https://yelleknek.github.io/DMAR/reference/ci_proportion.md),
[`ci_two_phase()`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Liver scan data from McNeil, Keeler, and Adelstein (1975), the
# example of Nam (1997, p. 1427): 344 patients classified as positive
# or negative by scan (x) and by pathology (y). Only the 32 and 27
# patients on whom the two disagree inform the difference in the
# positive rates.
liver <- matrix(c(231, 27, 32, 54), nrow = 2,
                dimnames = list(scan = c("positive", "negative"),
                                pathology = c("positive", "negative")))
liver
#>           pathology
#> scan       positive negative
#>   positive      231       32
#>   negative       27       54
ci_proportion_paired(n_10 = liver["positive", "negative"],
                     n_01 = liver["negative", "positive"],
                     N = sum(liver))
#>  term       value 
#>  ci_lower   -0.030
#>  difference 0.015 
#>  ci_upper   0.059 
#>  n_10       32    
#>  n_01       27    
#>  N          344   
#> 
#> Confidence level: 95%

# From raw outcomes: two judges assign the same 200 cases to one of
# three diagnoses. diagnosis_agreement holds one row per cell of their
# agreement table, so repeating each row by its frequency gives one row
# per case. How much more often does Judge A (x) than Judge B (y)
# diagnose a personality disorder?
cases <- diagnosis_agreement[rep(seq_len(nrow(diagnosis_agreement)),
                                 diagnosis_agreement$frequency), ]
ci_proportion_paired(x = cases$judge_a == "Personality disorder",
                     y = cases$judge_b == "Personality disorder")
#>  term       value
#>  ci_lower   0.036
#>  difference 0.100
#>  ci_upper   0.166
#>  n_10       32   
#>  n_01       12   
#>  N          200  
#> 
#> Confidence level: 95%

# When the two conditions agree on all 40 units, the Wald interval has
# zero width, as though the difference were known exactly; the score
# interval still reflects how few units were observed.
ci_proportion_paired(n_10 = 0, n_01 = 0, N = 40)
#>  term       value 
#>  ci_lower   -0.088
#>  difference 0.000 
#>  ci_upper   0.088 
#>  n_10       0     
#>  n_01       0     
#>  N          40    
#> 
#> Confidence level: 95%
ci_proportion_paired(n_10 = 0, n_01 = 0, N = 40, method = "wald")
#>  term       value
#>  ci_lower   0.000
#>  difference 0.000
#>  ci_upper   0.000
#>  n_10       0    
#>  n_01       0    
#>  N          40   
#> 
#> Confidence level: 95%
```
