# Sample Size Planning for the Coefficient of Variation Given the Goal of Accuracy in Parameter Estimation Approach to Sample Size Planning

Takes a planning value of the population coefficient of variation and a
desired confidence interval width, and returns the smallest sample size
at which the interval of
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md) is
expected to be that narrow (Kelley, 2007). With `assurance` it returns
the sample size of the published procedure of that article, at which the
realized two-sided interval is that narrow with about the stated
probability. The interval and the plan assume a normal population (see
Details).

The population coefficient of variation may be given directly as
`true_cv` or through `true_mean` and `true_sd`, in which case `true_cv`
is taken as `true_sd / true_mean`. It must be positive.

## Usage

``` r
ss_aipe_cv(
  true_cv = NULL,
  width = NULL,
  conf_level = 0.95,
  assurance = NULL,
  true_mean = NULL,
  true_sd = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL
)
```

## Arguments

- true_cv:

  Population coefficient of variation the sample size is planned for
  (kappa in Kelley, 2007). The planner treats this value as the truth;
  [`ss_aipe_cv_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md)
  shows what the plan delivers when the truth (`true_cv` there) differs
  from the value planned with (`supposed_cv` there).

- width:

  Desired (full) width of the confidence interval

- conf_level:

  Confidence interval coverage; 1-Type I error rate. A level above
  .999999998, whose tail areas are below 1e-9, is refused (see "Small
  tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md))

- assurance:

  Optional. The probability that the realized interval is no wider than
  `width`: a single number strictly between 0 and 1 (for example .80,
  .90, or .99). A value such as .50 may be used; it asks that the
  realized interval be no wider than `width` in about half of the
  studies that could be run, where the plan without `assurance` asks
  that the width at the expected value of the sample coefficient be no
  larger than `width`. Near .50 the realized probability can be several
  points below the request and the sample size can be one smaller than
  without `assurance`; see Details for the procedure and the numbers.
  Any other value is refused with a condition of class
  `dmar_assurance_invalid`. `assurance` cannot be given with
  `alpha_lower = 0`, a one-sided upper confidence bound (see **Edge
  cases**).

- true_mean:

  Population mean, a positive number (specified with `true_sd` when
  `true_cv` is not specified)

- true_sd:

  Population standard deviation, a positive number (specified with
  `true_mean` when `true_cv` is not specified)

- alpha_lower:

  Type I error for the lower confidence limit, given with `alpha_upper`
  when `conf_level` is `NULL`: a single number that is at least 0 and
  below 1

- alpha_upper:

  Type I error for the upper confidence limit. A nonzero tail area below
  1e-9 is refused, as in
  [`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md)

## Value

A `data.frame` with the columns `term` and `value` and one row,
`necessary_N`, the planned sample size, with or without `assurance`.

## Details

**Estimand and interval.** The quantity planned for is the population
coefficient of variation \\\kappa = \sigma / \mu\\ of a normal
population with a positive mean, estimated from \\N\\ independent
observations by \\k = s / M\\. The interval is the one
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md) reports,
from the noncentral *t* distribution of \\\sqrt{N} / k\\, two-sided at
`conf_level` (or with the tail areas `alpha_lower` and `alpha_upper`
when `conf_level` is `NULL`).

**The expected-width plan.** At each candidate sample size \\N\\ the
width is evaluated at the value Kelley (2007, Equation 3) gives for the
expected value of the sample coefficient of variation given `true_cv`
and \\N\\, \$\$E\[k \mid \kappa, N\] \approx \kappa \left\[1 -
\frac{1}{4(N - 1)} + \frac{\kappa^2}{N} + \frac{1}{2(N -
1)^2}\right\],\$\$ and the smallest \\N\\ whose width at that value is
no larger than `width` is returned. The expression is the approximation
of Sokal and Braumann (1980, p. 51, Expression 3) for normal data. It
stops at the third-order term and is above the expected value at small
\\N\\: for a coefficient of .25, by 3.1 percent at \\N = 5\\, 0.6
percent at \\N = 10\\, and 0.1 percent at \\N = 20\\. The function keeps
it because the published tables were computed with it. The search starts
at \\N = 4\\. The interval is computed exactly; for how this relates to
the sample sizes tabled in Kelley (2007), see **The published tables**
below.

The width of the interval is not a linear function of the sample
coefficient of variation, so the width at the value of Equation 3 is not
exactly the expected width. Over 62 conditions (coefficients of
variation of .05 to .50, widths of .025 to .20, and confidence levels of
.90, .95, and .99, with plans from 5 to 1,500) the returned \\N\\ was
the smallest \\N\\ with an expected width no larger than `width` in 47
and one smaller in 15; where it was one smaller the expected width
exceeded `width` by at most 0.7 percent.

**The assurance step.** With `assurance` the function runs the
computation that produced the tables of Kelley (2007). At the
expected-width sample size \\N_0\\ it finds \\\kappa\_\gamma\\, the
value that the sample coefficient exceeds with probability
`1 - assurance`, from the noncentral *t* distribution of Equation 5 of
the article. It then takes the value of Equation 3 at \\\kappa\_\gamma\\
and \\N_0\\ as the planning value and runs the expected-width search for
it, which applies Equation 3 to that value again at each \\N\\. The
quantile and the intervals are computed exactly (see
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)). The
article's text (p. 759) describes the search as run for
\\\kappa\_\gamma\\ itself. Recomputed with R's own noncentral *t*
functions, as in 2007, the computation described here gives the printed
value in 524 of the 528 assurance cells of Tables 1 to 3 (the other four
are in the two rows that show a misprint), where the search for
\\\kappa\_\gamma\\ itself gives the printed value in 365. The two
computations give different plans in 160 cells.

**The published tables.** Kelley (2007) computed Tables 1 to 3 (pp.
760–762) with R's noncentral *t* distribution, which is a normal
approximation above a noncentrality of 37.62. Here the noncentrality is
\\\sqrt{N}/k\\, which narrow intervals pass, and the function computes
there exactly, so it does not return every printed value. Two rows of
the tables show a misprint. Table 1 (90 percent confidence), assurance
.80, width .125, coefficient of variation .40 prints 127, between 69 and
122 in its row; the procedure gives 93. Table 3 (99 percent confidence),
assurance .80, width .010 prints 14,279; 28,321; 50,209; 38,120; and
82,817 for coefficients of variation of .30 to .50, which do not
increase; the procedure gives 14,277; 20,473; 28,320; 38,118; and
50,207, so the printed 28,321 and 50,209 are its values for .40 and .50
(to within the approximation), one and two columns to the left of where
they belong. Apart from the four misprinted cells (127; 28,321; 50,209;
and 82,817), the function returns the printed value in 468 of the 788
cells, one to three fewer cases in 232, and one to three more in 88. It
returns every printed value with a width of .125 or more at 90 and 95
percent confidence, the misprint aside, and 9 of the 120 with a width of
.075 or less at 99 percent. The worked example of the article (pp.
759–760), a coefficient of variation of .25, a width of .10, and 99
percent confidence, gives 101, and 142 with an assurance of .99, where
the article gives 99 and 141: at 99 the interval at the expected value
is .1006 wide, and .0999 by the approximation.

The procedure is a close approximation and not an exact search, so the
probability that the realized interval is no wider than `width` at the
returned sample size is near `assurance` and not exactly it. In the 488
printed cells with an assurance and at most 5,000 cases that probability
was .744 to .902 where .80 was requested (below .80 in 37 of 244 cells)
and .9916 to .9997 where .99 was requested. Away from the published
tables, and more so as the request falls toward .50, it can fall below
the request: for a coefficient of variation of .22, a width of .09, a 95
percent interval, and an assurance of .60 the function returns 57, at
which the probability is .587, and one more case gives .621. In 105
plans computed at 95 percent confidence (coefficients of variation from
.05 to .50, widths of a quarter, a half, and three quarters of the
coefficient) the probability was below the request by at most .038 for a
request of .50, .036 for .55, .058 for .60, .032 for .70, and .001 for
.80, and it was not below the request in any plan for .90 and .99. A
request of .50 does not reproduce the plan without `assurance`: in 12 of
the 15 plans it returned a sample size one smaller than the plan without
`assurance` does, and for a coefficient of variation of .05 and a width
of .0375 it returns 17, at which the probability is .462, where the plan
without `assurance` returns 18. In the other direction, when `width` is
close to `true_cv` and `assurance` is high the returned sample size can
be well above the smallest sufficient one: for a coefficient of
variation of .5, a width of .5, a 95 percent interval, and an assurance
of .99 the function returns 57, and the probability reaches .99 at 37.
[`ss_aipe_cv_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md)
estimates the probability for any plan by simulation.

**Populations that are not normal.** The interval and the plan assume a
normal population. With a skewed population the coverage of the interval
departs from the stated level, and a larger sample does not repair it
(see [`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md)). In
40,000 samples at each of the three plans for a coefficient of variation
of .50, a width of .20, and 95 percent confidence (79 without
`assurance`, 97 with an assurance of .80, and 142 with .99), coverage
was .97 for a gamma population with that coefficient of variation, and
for a lognormal one it was .915, .913, and .905, lowest at the largest
plan. The proportion of intervals no wider than `width` was .55, .88,
and 1.00 for the gamma population and .59, .84, and .99 for the
lognormal one. For a normal population that probability, computed from
the noncentral *t* distribution, is .53, .84, and 1.00.

**Edge cases.** When `width` is already met at \\N = 4\\ the function
returns 4. A `width` that no sample size reaches within 100,000 steps of
the search stops with an error. With `alpha_upper = 0` the interval is a
one-sided lower confidence bound whose upper limit is infinite at every
sample size, so no plan exists, and the call is refused at once with a
condition of class `dmar_argument_invalid`. With `alpha_lower = 0`, an
upper bound, the lower limit is 0 and the width is the bound itself,
which comes down toward the coefficient of variation as \\N\\ grows, so
a `width` larger than `true_cv` can be planned without `assurance`. With
`assurance` the call is refused, with class `dmar_argument_invalid`: the
assurance step was written for a two-sided interval and does not give
the requested probability for a bound. For a coefficient of variation of
.2, a width of .25, `conf_level = NULL`, `alpha_lower = 0`, and
`alpha_upper = .05`, where the plan without `assurance` is 37, the step
gives 33, 91, and 219 for an assurance of .50, .80, and .90. The
probability that the bound is no larger than the width is .470, .878,
and .998 at those sizes, and the smallest sample sizes at which it
reaches the request are 36, 74, and 98. With `alpha_upper` below .5 the
bound at the value of Equation 3 stays above `true_cv` at every sample
size the search can reach: over `true_cv` from .001 to 5 and \\N\\ from
4 to 100,004 the smallest bound was 1.000006 times `true_cv`. (The value
of Equation 3 is itself below `true_cv` when `true_cv` is below
\\\sqrt{N(N - 3)} / (2(N - 1))\\, which is .33 at \\N = 4\\, .46 at \\N
= 10\\, and .49 at \\N = 20\\, and approaches .5.) A `width` at or below
`true_cv` is therefore never reached, and the call is refused at once,
also with class `dmar_argument_invalid`. With a large coefficient of
variation the first sample sizes of the search can be so small that the
upper limit of the interval is infinite (at 99 percent confidence and a
coefficient of .45, \\N = 4\\ is one). An infinite interval is never
narrow enough, so the search goes on, and the interval at the returned
sample size is finite.

**Comparison with MBESS.** MBESS 3.4.0 and later evaluate the width at
`true_cv` multiplied by \\(1 + 1/(4N))\\, the nearly unbiased estimate
of Kelley (2007, Equation 4) applied to the population value, which is
larger than `true_cv` where the expected value is smaller, and use R's
approximation to the noncentral *t* beyond a noncentrality of 37.62.
Their sample sizes are one larger than the tables in Kelley (2007) in
many conditions, and up to three larger with `assurance`; for the worked
example they are 100 and 142, where this function gives 101 and 142.

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research, 51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods, 39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Sokal, R. R., & Braumann, C. A. (1980). Significance tests for
coefficients of variation and variability profiles. *Systematic Zoology,
29*(1), 50–66.
[doi:10.1093/sysbio/29.1.50](https://doi.org/10.1093/sysbio/29.1.50)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on confidence intervals for effect
sizes.)

## See also

[`ss_aipe_cv_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md),
[`cv`](https://yelleknek.github.io/DMAR/reference/cv.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Suppose one wishes to have a confidence interval with an expected width of .10
# for a 99% confidence interval when the population coefficient of variation is .10.
ss_aipe_cv(true_cv = .1, width = .1, conf_level = .99)
#>  term        value
#>  necessary_N 20   
#> 
#> Confidence level: 99%

# The same planning problem parameterized by the population mean and standard
# deviation: true_mean = 10 and true_sd = 1 imply the same coefficient of variation, .10.
ss_aipe_cv(true_mean = 10, true_sd = 1, width = .1, conf_level = .99)
#>  term        value
#>  necessary_N 20   
#> 
#> Confidence level: 99%

# Ensuring that the confidence interval will be sufficiently narrow with a 99%
# certainty for the situation above.
ss_aipe_cv(true_cv = .1, width = .1, conf_level = .99, assurance = .99)
#>  term        value
#>  necessary_N 32   
#> 
#> Confidence level: 99%
```
