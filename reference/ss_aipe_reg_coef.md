# Sample Size Planning for a Single Regression Coefficient (AIPE)

Computes the necessary sample size for the confidence interval on a
targeted regression coefficient \\\beta_j\\, the quantity
`true_reg_coef_j` stands for, in a multiple regression with *p*
predictors to be no wider than a user-specified value. The coefficient
is on the scale `true_sd_Y` and `true_sd_X` describe: with both 1, the
default, its population value is the standardized coefficient. This is
the accuracy in parameter estimation (AIPE) framework of Kelley and
Maxwell (2003), targeted at a specific coefficient rather than at the
omnibus \\R^2\\. The default (`noncentral = FALSE`) plans for the
central *t* interval for the coefficient in the units of the data, the
interval [`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md)
reports; `noncentral = TRUE` plans for the noncentral *t* interval for
the standardized coefficient (Kelley, 2007a), the interval
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md)
reports, with unit standard deviations. The central *t* interval for the
standardized coefficient, the interval of Kelley and Maxwell (2003,
Equation 1), is a third interval, and
[`ss_aipe_src`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src.md)
with `noncentral = FALSE` plans for it; see "The interval planned for"
in Details. Supplying `assurance` returns the smallest *N*, not below
the plan without it (the plan for the width at the expected values of
the sample squared multiple correlations), at which the realized
interval is no wider than `width` with at least that probability,
computed exactly for the central interval; for the noncentral interval
it returns by default the published approximation that MBESS uses, which
does not compute the probability, and with
`assurance_method = "simulation"` the smallest *N* at which a simulated
estimate of the probability reaches the request.

## Usage

``` r
ss_aipe_reg_coef(
  true_R2_Y_X = NULL,
  true_R2_j_X_without_j = NULL,
  p = NULL,
  true_reg_coef_j = NULL,
  width,
  which_width = "Full",
  true_sd_Y = 1,
  true_sd_X = 1,
  true_cor_XX = NULL,
  true_cor_YX = NULL,
  which_predictor = NULL,
  noncentral = FALSE,
  alpha_lower = NULL,
  alpha_upper = NULL,
  conf_level = 0.95,
  assurance = NULL,
  assurance_method = c("analytic", "simulation"),
  G = 20000,
  seed = NULL
)
```

## Arguments

- true_R2_Y_X:

  Population value of \\\rho^2\_{Y \cdot X_1, \ldots, X_p}\\, the
  squared multiple correlation of the outcome *Y* with all *p*
  predictors.

- true_R2_j_X_without_j:

  Population value of \\\rho^2\_{X_j \cdot X\_{-j}}\\, the squared
  multiple correlation when the *j*th predictor is regressed on the
  remaining \\p - 1\\ predictors. Quantifies the multicollinearity faced
  by the targeted coefficient.

- p:

  The number of predictor variables.

- true_reg_coef_j:

  Population value \\\beta_j\\ of the regression coefficient for the
  *j*th predictor, the predictor of interest, on the scale `true_sd_Y`
  and `true_sd_X` describe: with both 1 it is the standardized
  coefficient, and otherwise the unstandardized one, the standardized
  coefficient times `true_sd_Y / true_sd_X` (Kelley & Maxwell, 2003, p.
  316); \\b_j\\ is its sample estimate. On the covariance structure path
  it may be left out; when it is given it must equal the coefficient the
  correlations imply on that scale.

- width:

  Desired (full) width of the two-sided confidence interval on
  \\\beta_j\\.

- which_width:

  Which portion of the confidence interval `width` refers to. Only
  `"Full"` (or an abbreviation of it, such as `"F"`) is currently
  implemented; `"Lower"`, `"Upper"`, and their abbreviations are
  refused.

- true_sd_Y:

  Population standard deviation of *Y*, which with `true_sd_X` sets the
  scale of the coefficient and of `width`.

- true_sd_X:

  Population standard deviation of the *j*th predictor.

- true_cor_XX:

  Population correlation matrix for the *p* predictor variables. If
  supplied with `true_cor_YX`, `true_R2_Y_X` and `true_R2_j_X_without_j`
  are derived from the covariance structure. With one predictor it is
  `matrix(1)` (or the number 1), and `p` and `which_predictor` are then
  1 when not given.

- true_cor_YX:

  Length-*p* vector of population correlations between *Y* and the *p*
  predictors.

- which_predictor:

  Which of the *p* predictors is the targeted coefficient.

- noncentral:

  If `TRUE`, plans for the noncentral *t* interval for the standardized
  coefficient (Kelley, 2007a), the interval of
  [`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md). If
  `FALSE` (the default), plans for the central *t* interval for the
  coefficient in the units of the data, the interval of
  [`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md), even
  when both standard deviations are 1; the central interval for the
  standardized coefficient is planned by
  [`ss_aipe_src`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src.md)
  with `noncentral = FALSE`. The noncentral path requires
  `true_sd_Y = true_sd_X = 1` (*i.e.*, a standardized solution).

- alpha_lower:

  Type I error rate for the lower confidence limit, at least 1e-9 and at
  most .5. A tail area below 1e-9 is refused: the noncentral *t* limits
  lose accuracy there, and the central plan keeps the same floor.

- alpha_upper:

  Type I error rate for the upper confidence limit, at least 1e-9 and at
  most .5.

- conf_level:

  Confidence level (i.e., \\1 - \alpha\\, where \\\alpha\\ is the Type I
  error rate). Default `0.95`. Mutually exclusive with `alpha_lower` and
  `alpha_upper`. A level above .999999998, whose tail areas are below
  1e-9, is refused.

- assurance:

  Optional probability with which the realized confidence interval is to
  be no wider than `width`, a single number above .5 and below 1. When
  `NULL` (the default), the plan is the smallest *N* at which the
  interval is no wider than `width` when the two sample squared multiple
  correlations take their expected values (see "The interval planned
  for" in Details): the width at those expected values, not the mean of
  the width over samples.

- assurance_method:

  How the probability that the realized interval is no wider than
  `width` is found, used only with `assurance` (see "The assurance step"
  in Details). `"analytic"` (the default) computes without simulation:
  for the central interval (`noncentral = FALSE`) the probability is
  exact, and the function returns the smallest sample size, at or above
  the plan without `assurance`, at which it is at least `assurance`,
  reported as `assurance_attained`; for the noncentral interval it
  returns the published approximation that MBESS uses, which does not
  compute the probability (`assurance_attained` is `NA`) and can attain
  less than `assurance`. `"simulation"`, offered for the noncentral
  interval only, returns the smallest sample size, at or above the plan
  without `assurance`, at which the probability estimated from `G`
  simulated samples is at least `assurance`, and reports that estimate.
  With `noncentral = FALSE` `"simulation"` is refused with a condition
  of class `dmar_argument_invalid`, since a simulated estimate of an
  exact probability would only add Monte Carlo error.

- G:

  Number of simulated samples at each sample size the noncentral
  assurance step tries, and in its final estimate of the attained
  probability; a whole number of at least 100 (default 20000). Used and
  checked only with `noncentral = TRUE`, `assurance`, and
  `assurance_method = "simulation"`.

- seed:

  Optional single number. When given, the noncentral assurance step is
  reproducible, and the random number stream of the calling session is
  left as it was. With `NULL` (the default) the draws come from the
  session's current stream. Used only where `G` is.

## Value

A `data.frame` with columns `term` and `value`:

- `necessary_N`:

  the necessary total sample size.

- `assurance`:

  the requested probability; present only when `assurance` was given.

- `assurance_attained`:

  the probability, at `necessary_N`, that the realized interval is no
  wider than `width`; present only when `assurance` was given. Exact for
  the central interval; for the noncentral interval `NA` with
  `assurance_method = "analytic"`, whose approximation does not compute
  it, and with `"simulation"` a Monte Carlo estimate, whose standard
  error is the attribute `assurance_mc_se`.

- `G`:

  the number of simulated samples; present only when the noncentral
  assurance step simulated.

With `assurance` the printed footer names the assurance method and says
how `assurance_attained` was found, or that the approximation does not
compute it; the attributes `assurance_method` (`"analytic"` or
`"simulation"`) and `assurance_probability` (`"exact"`, `"simulation"`,
or `"not computed"`) record the same.

## Details

**Calling conventions.** The function offers several mutually exclusive
ways to supply the population information needed to plan \\N\\; the user
picks the one most aligned with their available planning values. Specify
exactly one of:

- *Covariance structure path.* Supply `true_cor_XX` and `true_cor_YX`
  along with `which_predictor`. The function derives \\\rho^2\_{Y \cdot
  X}\\ and \\\rho^2\_{X_j \cdot X\_{-j}}\\ (the values of `true_R2_Y_X`
  and `true_R2_j_X_without_j`) from the covariance structure, and the
  population coefficient \\\beta_j\\, the *j*th element of
  \\R\_{XX}^{-1} r\_{XY}\\ (\\R\_{XX}\\ is `true_cor_XX` and \\r\_{XY}\\
  is `true_cor_YX`) times `true_sd_Y / true_sd_X`, against which a
  supplied `true_reg_coef_j` is checked.

- *Squared multiple correlations path.* Supply `true_R2_Y_X`,
  `true_R2_j_X_without_j`, `p`, and `true_reg_coef_j` directly when the
  user already has these planning values from prior research and does
  not need to specify the full covariance structure.

- *Standardized solution.* For the `noncentral = TRUE` path, set
  `true_sd_Y = true_sd_X = 1`; the result returns the standardized
  regression coefficient sample size.

**The interval planned for.** With `noncentral = FALSE` the interval is
the central *t* interval for the coefficient in the units of the data,
\\b_j\\ plus the *t* quantiles (on \\N - p - 1\\ degrees of freedom) at
`alpha_lower` and `1 - alpha_upper` times its standard error, the
interval [`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md)
and
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md)
report from unstandardized data, and this is so whatever `true_sd_Y` and
`true_sd_X` are: with both 1 the coefficient is on the standardized
scale, but the interval planned for is still the one computed in the
units of the data, whose standard error comes from the residual sums of
squares of *Y* and of \\X_j\\. The plan without `assurance` is the
smallest *N* at which that interval, computed at the expected values of
the two sample squared multiple correlations, is no wider than `width`
(Kelley & Maxwell, 2003, Equations 1 and 2, with the *t* quantile); the
two tails exchanged give the same plan. The plan is the smallest such
*N* of at least \\p + 2\\, the smallest sample with a residual degree of
freedom, and the expected squared multiple correlation of \\X_j\\ with
the other predictors is that of a regression on those \\p - 1\\
predictors (0 with one predictor; see "Comparison with MBESS"). Kelley
and Maxwell's Equation 1 is the central interval for the standardized
coefficient, whose standard error is formed from the two sample squared
multiple correlations; its plan without `assurance` is this one at unit
standard deviations, but the probability that it is no wider than
`width` is not the probability the assurance step below computes, and
[`ss_aipe_src`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src.md)
with `noncentral = FALSE` plans its assurance. With `noncentral = TRUE`
the interval is the noncentral *t* interval for the standardized
coefficient of
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md): the
limits of the noncentrality parameter at the observed *t* statistic,
times the standard error \\\sqrt{(1 - R^2)/\[(1 - R^2\_{X_j \cdot
X\_{-j}})(N - p - 1)\]}\\, and the plan without `assurance` is the
smallest *N*, of at least \\p + 2\\, at which it is no wider than
`width` at the expected squared multiple correlations and \\\beta_j\\.
That standard error treats the sample standard deviations as constants,
which for random predictors is an approximation (Kelley & Maxwell, 2003,
pp. 315–316). The central path is fast and adequate at moderate to large
*N*; the noncentral path is the one to use when the analysis will report
the noncentral interval for the standardized coefficient. A tail of 0,
or two of .5, is refused, and so is a tail area below 1e-9 or a
`conf_level` above .999999998.

**The assurance step.** A plan without `assurance` makes the interval as
narrow as the target when the sample values equal their expectations; in
roughly half of all samples it is wider. With `assurance` the function
plans so that the realized interval is no wider than `width` with at
least that probability, for random predictors from a multivariate normal
population.

- *The central interval* (`noncentral = FALSE`). The only method is
  `"analytic"`: the function computes \\P_N\\, the probability that the
  realized interval is no wider than `width` at sample size *N*,
  exactly, and returns the smallest *N*, at or above the plan without
  `assurance`, at which \\P_N\\ is at least `assurance`, found by
  doubling a step upward and then bisecting. `assurance_attained` is
  \\P_N\\ at the returned *N*. The squared standard error of \\b_j\\ is
  the residual sum of squares of *Y* on all *p* predictors over \\N -
  p - 1\\, divided by the residual sum of squares of \\X_j\\ on the
  other \\p - 1\\. Under multivariate normality the first, over
  \\\sigma^2_Y (1 - \rho^2\_{Y \cdot X})\\, is chi square on \\N - p -
  1\\ degrees of freedom whatever the predictors are, so it is
  independent of the second, which over \\\sigma^2\_{X_j} (1 -
  \rho^2\_{X_j \cdot X\_{-j}})\\ is chi square on \\N - p\\ degrees of
  freedom. The squared standard error is then \\K F / (N - p)\\, with
  \\K = \sigma^2_Y (1 - \rho^2\_{Y \cdot X}) / \[\sigma^2\_{X_j} (1 -
  \rho^2\_{X_j \cdot X\_{-j}})\]\\ and *F* an *F* variate on \\N - p -
  1\\ and \\N - p\\ degrees of freedom, and with \\q\\ the sum of the
  two *t* quantiles of the interval, \$\$P_N = \Pr(F \le w^2 (N - p) /
  (q^2 K)),\$\$ where \\w\\ is `width`. (With fixed predictors the
  second sum of squares is a constant and the same event is a chi square
  probability; the plan is for random predictors, the model
  [`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md)
  simulates.) The units of *Y* and \\X_j\\ cancel, so changing them,
  with `width` changed to match, leaves the plan as it was, and the sign
  of \\\beta_j\\ does not enter.

- *The noncentral interval* (`noncentral = TRUE`). The width depends on
  the observed *t* statistic and on both sample squared multiple
  correlations, and no closed form for its distribution is used. With
  `assurance_method = "analytic"`, the default, the plan is the
  approximation that MBESS uses (Kelley, 2007a, 2007b): the plan without
  `assurance` made again with a lower confidence limit (at level
  `assurance`) of the squared multiple correlation in place of its
  population value, and not below the first. It does not compute the
  probability it attains, so `assurance_attained` is `NA`, and the
  printed footer says so. Kelley and Maxwell (2008, pp. 186–187) wrote
  that no satisfactory method then existed for the standardized
  coefficient and described finding, by simulation, the smallest sample
  size that attains the requested degree of certainty, which is what
  `assurance_method = "simulation"` does: \\P_N\\ is estimated as the
  proportion of `G` samples, drawn from the exact sampling distribution
  of the sums of squares and products (Wishart, through the entries of
  its triangular factor that the statistics depend on), whose interval
  is no wider than `width`. The same `G` samples serve every sample size
  the search tries (common random numbers: each chi square entry is the
  quantile, at the degrees of freedom of that *N*, of one uniform draw),
  so the estimate changes smoothly with *N*. The search assumes that it
  does not decrease as *N* grows, which held in every set checked
  although a single simulated interval can be wider at \\N + 1\\ than at
  *N* (see below), and under that assumption the step returns the
  smallest *N* at which the estimate is at least `assurance`. The
  distribution of the width depends on the population through *p*,
  \\\beta_j\\, the two squared multiple correlations, and the
  correlation of *Y* with \\X_j\\. On the covariance structure path that
  correlation is given; on the squared multiple correlations path it is
  not, and the step uses the smallest estimated probability over 21
  equally spaced values of it, from the smallest to the largest that the
  three values allow. That is an approximate minimum over the range:
  between two of the values the probability can be a little lower than
  at either, so for a correlation there the plan can fall a little short
  of the request. (With \\p = 2\\ only the two ends of the range are
  possible, and the step uses both.) `assurance_attained` is the
  estimate at the returned *N*, so it is at least `assurance`, and its
  Monte Carlo standard error is the attribute `assurance_mc_se`. The
  exact \\P_N\\ differs from the estimate by that error, so the returned
  *N* can be a case or two from the one the exact probability would
  give, and \\P_N\\ there can be a little below `assurance`; give `seed`
  for a plan that can be reproduced, and a larger `G` for a more precise
  one. Reversing the sign of \\\beta_j\\ (with `alpha_lower` and
  `alpha_upper` exchanged when they differ) gives the same plan.

For both intervals \\P_N\\ did not decrease from one *N* to the next at
or above the plan without `assurance`: for the central interval in 1,530
planning sets (*p* from 1 to 15, squared multiple correlations from 0 to
.9, widths from .05 to 2, and tails of .025 and .025, .005 and .005, and
.01 and .04) over every *N* from that plan to three times it, and for
the noncentral interval in five sets by simulation. With the samples
held fixed across *N*, the estimate did not decrease from one *N* to the
next in seven sets (*p* from 1 to 10), and of 95 million comparisons of
a simulated interval at *N* and at \\N + 1\\, 3208, all at *N* below 60,
were wider at \\N + 1\\. Below the plan without `assurance` \\P_N\\ can
decrease, which is why the search starts there.

The approximations can attain more or less than the request. For two
predictors correlated .3 with correlations of .4 and .3 with *Y*, a
width of .3, and an assurance of .90, the chi square approximation of
Kelley and Maxwell (2003, Equation 3; 2008, Equation 37) for the central
interval, the plan without `assurance` with the variance of \\b_j\\
multiplied by the `assurance` quantile of a chi square on \\N - 1\\
degrees of freedom over \\N - p - 1\\, is 181, at which \\P_N\\ is .852,
and the exact plan is 187 (\\P_N\\ = .902); the noncentral
approximation, the default plan for that interval, is 179, at which the
simulated probability is .948, and `assurance_method = "simulation"`
returns 176 (.908 with `seed = 1`). For the five exchangeable predictors
of Kelley and Maxwell (2003, pp. 312–313), correlated .4 with one
another and .3 with *Y*, at a width of .2 and an assurance of .80, the
chi square approximation is 486 (\\P_N\\ = .765) and the exact plan 491
(.800). Each was checked against the proportion of intervals no wider
than the target among samples of raw multivariate normal data. These
probabilities, and those of the first example, are for the central
interval in the units of the data, the interval this function plans for.
Kelley and Maxwell's example is of the standardized coefficient, with
the central interval for it (their Equation 1), which is computed from
standardized data with the sample standard deviations treated as
constants; its width varies less, and it is no wider than .2 with
probability .877 at 486 and .911 at 491 (their own simulation, p. 317,
reported .816 at their plan of 485, where raw multivariate normal data
give .871 with the *t* quantile of their Equation 1, or .887 with the
normal quantile of their plan). For that interval
[`ss_aipe_src`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src.md)
with `noncentral = FALSE`, `assurance_method = "simulation"`, and
`seed = 1` plans 478 (.807, and .797 at 477), and in the first example
171 (.914, and .895 at 170), where the plans above reach .986 at 181 and
.997 at 187; its plan is simulated, and another seed can give one a case
away (170 for one of seeds 1 to 30 in the first example).

**Comparison with MBESS.** The plans without `assurance` of
`ss.aipe.reg.coef()` in MBESS differ from these in two ways. They take
the expected squared multiple correlation of \\X_j\\ with the other
predictors to be that of a regression on *p* predictors, where that
regression has \\p - 1\\, which overstates it and so the width. And
their search goes up one case at a time from a start, the larger of
\\2p + 1\\ and ten below the normal approximation for the central
interval and six below the central plan for the noncentral one, so it
returns the start when a smaller sample size already meets the width:
for `true_R2_Y_X = .8`, `true_R2_j_X_without_j = 0`, `p = 10`, and
`width = 2`, 21 where 14 does. In 1,530 planning sets the two
differences make the central plan here one case smaller in 1,329 and 2
to 13 cases smaller in 104, and in 810 sets the noncentral plan one
smaller in 596 and 2 to 10 smaller in 177. The approximations start from
these plans (Kelley and Maxwell's example above is 487 in MBESS and 486
here). With unequal tails `ss.aipe.reg.coef()` also takes both *t*
quantiles of its central plans at the average of the two tails, so it
plans for an interval narrower than the one
[`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md) reports;
the central plans here take the quantile at each tail.

## References

Kelley, K. (2007a). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K. (2007b). Methods for the behavioral, educational, and social
sciences: An R package. *Behavior Research Methods, 39*(4), 979–984.
[doi:10.3758/BF03192993](https://doi.org/10.3758/BF03192993)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons of means and
Chapter 6 on trend analysis.)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

## See also

[`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md),
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Covariance structure path: supply the population correlation
#    matrix and the population YX cross-correlations. Five predictors
#    in an exchangeable structure (all pairwise correlations 0.5,
#    all Y-X correlations 0.3).
true_cor_YX <- c(.3, .3, .3, .3, .3)
true_cor_XX <- rbind(c(1, .5, .5, .5, .5), c(.5, 1, .5, .5, .5),
                c(.5, .5, 1, .5, .5), c(.5, .5, .5, 1, .5),
                c(.5, .5, .5, .5, 1))

# Closed-form (central t) planning, targeting the first predictor's
# standardized coefficient.
ss_aipe_reg_coef(width = .10, which_width = "Full",
                 true_sd_Y = 1, true_sd_X = 1,
                 true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
                 which_predictor = 1, noncentral = FALSE,
                 conf_level = .95)
#>  term        value
#>  necessary_N 2184 
#> 
#> Confidence level: 95%

# With an assurance of .85: the smallest N at which the realized
# interval is no wider than 0.10 with probability at least .85. The
# last row is that probability at the returned N, computed exactly.
ss_aipe_reg_coef(width = .10, which_width = "Full",
                 true_sd_Y = 1, true_sd_X = 1,
                 true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
                 which_predictor = 1, noncentral = FALSE,
                 conf_level = .95, assurance = .85)
#>  term               value
#>  necessary_N        2281 
#>  assurance          0.850
#>  assurance_attained 0.850
#> 
#> Confidence level: 95%
#> 
#> Assurance method: analytic; assurance_attained is the exact probability that the interval is no wider than the width

# The noncentral t interval for the standardized coefficient. Required
# N differs at small to moderate samples.
ss_aipe_reg_coef(width = .10, which_width = "Full",
                 true_sd_Y = 1, true_sd_X = 1,
                 true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
                 which_predictor = 1, noncentral = TRUE,
                 conf_level = .95)
#>  term        value
#>  necessary_N 2190 
#> 
#> Confidence level: 95%

# With assurance its default plan is the approximation that MBESS uses,
# which does not compute the probability it attains; by simulation, the
# smallest N at which the estimated probability is at least .85. G = 2000
# keeps the example quick; a reported plan deserves the default G = 20000.
ss_aipe_reg_coef(width = .10, true_cor_XX = true_cor_XX,
                 true_cor_YX = true_cor_YX, which_predictor = 1,
                 noncentral = TRUE, assurance = .85)
#>  term               value
#>  necessary_N        2240 
#>  assurance          0.850
#>  assurance_attained <NA> 
#> 
#> Confidence level: 95%
#> 
#> Assurance method: analytic; the plan is the published approximation, which does not compute the probability that the interval is no wider than the width (assurance_attained is NA) and can attain less than the assurance. assurance_method = "simulation" computes it and returns the smallest plan that attains it
plan <- ss_aipe_reg_coef(width = .10, true_cor_XX = true_cor_XX,
                         true_cor_YX = true_cor_YX, which_predictor = 1,
                         noncentral = TRUE, assurance = .85,
                         assurance_method = "simulation",
                         G = 2000, seed = 113)
plan
#>  term               value
#>  necessary_N        2257 
#>  assurance          0.850
#>  assurance_attained 0.854
#>  G                  2000 
#> 
#> Confidence level: 95%
#> 
#> Assurance method: simulation; assurance_attained is estimated from G = 2000 simulated samples (Monte Carlo SE 0.008)
attr(plan, "assurance_mc_se")
#> [1] 0.007895695

# 2. Squared multiple correlations path: when the user has planning
#    values for true_R2_Y_X and true_R2_j_X_without_j directly (e.g.,
#    from a prior power analysis), without specifying the full
#    covariance structure. true_reg_coef_j is required on this path.
ss_aipe_reg_coef(true_R2_Y_X = 0.30, true_R2_j_X_without_j = 0.20,
                 p = 5, true_reg_coef_j = 0.25,
                 width = .15, which_width = "Full",
                 true_sd_Y = 1, true_sd_X = 1,
                 noncentral = FALSE, conf_level = .95)
#>  term        value
#>  necessary_N 606  
#> 
#> Confidence level: 95%

# 3. One predictor: its correlation matrix is matrix(1), and p and
#    which_predictor follow from it.
ss_aipe_reg_coef(width = .20, true_cor_XX = matrix(1), true_cor_YX = .4)
#>  term        value
#>  necessary_N 327  
#> 
#> Confidence level: 95%
```
