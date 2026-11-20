# Sample Size Necessary for the Accuracy in Parameter Estimation Approach for a Standardized Regression Coefficient of Interest

Plans the sample size from the accuracy in parameter estimation (AIPE)
approach for a standardized regression coefficient of interest: the
smallest *N* at which the confidence interval for the coefficient,
computed at the expected values of the sample squared multiple
correlations, is no wider than `width`, or, with `assurance`, the
published approximation for an interval no wider than `width` with at
least that probability or, by simulation, the smallest *N* at which it
is. Two intervals can be planned for. The default (`noncentral = TRUE`)
is the noncentral *t* interval for the standardized coefficient, the
interval
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md)
reports; `noncentral = FALSE` is the central *t* interval of Kelley and
Maxwell (2003, Equation 1), the interval of
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md)
with `noncentral = FALSE` for standardized data. For the coefficient in
the units of the data use
[`ss_aipe_rc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc.md).

## Usage

``` r
ss_aipe_src(
  true_R2_Y_X = NULL,
  true_R2_j_X_without_j = NULL,
  p = NULL,
  true_std_reg_coef_j = NULL,
  width,
  which_width = "Full",
  true_sd_Y = 1,
  true_sd_X = 1,
  true_cor_XX = NULL,
  true_cor_YX = NULL,
  which_predictor = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL,
  conf_level = 0.95,
  assurance = NULL,
  noncentral = TRUE,
  assurance_method = c("analytic", "simulation"),
  G = 20000,
  seed = NULL
)
```

## Arguments

- true_R2_Y_X:

  Population value of the squared multiple correlation coefficient

- true_R2_j_X_without_j:

  Population value of the squared multiple correlation coefficient
  predicting the *j*th predictor variable from the remaining *p*-1
  predictor variables

- p:

  The number of predictor variables

- true_std_reg_coef_j:

  Population value of the standardized regression coefficient for the
  *j*th predictor variable (the predictor of interest)

- width:

  The desired width of the confidence interval

- which_width:

  Which width (`"Full"`, `"Lower"`, or `"Upper"`) the width refers to
  (at present, only `"Full"` can be specified)

- true_sd_Y:

  The population standard deviation of *Y* (i.e., the dependent
  variable). It must be 1 (the default): the plan is for the
  standardized coefficient. Use
  [`ss_aipe_rc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc.md)
  for the coefficient in the units of the data.

- true_sd_X:

  The population standard deviation of the *j*th *X* variable (i.e., the
  predictor variable of interest). It must be 1 (the default), as
  `true_sd_Y` must.

- true_cor_XX:

  Population correlation matrix for the *p* predictor variables

- true_cor_YX:

  Population *p* length vector of correlation between the dependent
  variable (*Y*) and the *p* independent variables

- which_predictor:

  Identifies which of the *p* predictors is of interest

- alpha_lower:

  Type I error rate for the lower confidence interval limit, at least
  1e-9 and at most .5; a tail area below 1e-9 is refused, as in
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)

- alpha_upper:

  Type I error rate for the upper confidence interval limit, at least
  1e-9 and at most .5

- conf_level:

  Desired level of confidence for the computed interval (i.e., 1 - the
  Type I error rate); a level above .999999998, whose tail areas are
  below 1e-9, is refused

- assurance:

  Optional probability with which the realized confidence interval is to
  be no wider than `width`; when `NULL` (the default) the plan is for
  the width of the interval at the expected values of the sample squared
  multiple correlations (see
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)),
  not for the mean width over samples

- noncentral:

  `TRUE` (the default) plans for the noncentral *t* interval for the
  standardized coefficient, the interval of
  [`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md);
  `FALSE` plans for the central *t* interval of Kelley and Maxwell
  (2003, Equation 1), the estimate plus and minus the *t* quantiles of
  its two tails times the standard error formed from the two sample
  squared multiple correlations.

- assurance_method:

  How the plan with `assurance` is made, used only with `assurance`.
  Neither interval has an exact form for the probability that it is no
  wider than `width`. `"analytic"` (the default), computed without
  simulation, is the published approximation: that of MBESS for the
  noncentral interval, and the chi square approximation of Kelley and
  Maxwell (2003, Equation 3) for the central one; neither computes that
  probability (`assurance_attained` is `NA`), and either can attain more
  or less than `assurance`. `"simulation"` is the smallest sample size,
  at or above the plan without `assurance`, at which the probability
  estimated from `G` simulated samples is at least `assurance`, and
  reports that estimate as `assurance_attained` (see Details).

- G:

  Number of simulated samples at each sample size the assurance step
  tries, and in its final estimate of the attained probability; a whole
  number of at least 100 (default 20000). Used and checked only with
  `assurance` and `assurance_method = "simulation"`

- seed:

  Optional single number. When given, the assurance step is
  reproducible, and the random number stream of the calling session is
  left as it was. With `NULL` (the default) the draws come from the
  session's current stream. Used only where `G` is

## Value

A `data.frame` with columns `term` and `value`: `necessary_N`, the
necessary total sample size; with `assurance`, the rows `assurance` and
`assurance_attained`, which is `NA` for the published approximation of
`assurance_method = "analytic"` and, with `"simulation"`, the simulated
probability at `necessary_N` that the realized interval is no wider than
`width`, with its Monte Carlo standard error in the attribute
`assurance_mc_se` and the row `G`. With `assurance` the printed footer
names the assurance method and says how `assurance_attained` was found,
or that the approximation does not compute it; the attributes
`assurance_method` and `assurance_probability` (`"simulation"` or
`"not computed"`) record the same.

## Details

Not all of the arguments need to be specified, only those that provide
all of the necessary information so that the sample size can be
determined for the conditions specified.

The coefficient planned for is the population standardized regression
coefficient of the *j*th predictor, \\\beta_j\\, the quantity
`true_std_reg_coef_j` stands for; on the covariance structure path it is
the *j*th element of \\R\_{XX}^{-1} r\_{XY}\\, with \\R\_{XX}\\ given by
`true_cor_XX` and \\r\_{XY}\\ by `true_cor_YX`, and a supplied
`true_std_reg_coef_j` must equal it. A condition raised by
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
which makes the plans without `assurance`, names the arguments of this
function (`true_std_reg_coef_j`, not `true_reg_coef_j`).

**The interval planned for.** Both intervals use the standard error
\\\sqrt{(1 - R^2)/\[(1 - R^2\_{X_j \cdot X\_{-j}})(N - p - 1)\]}\\,
formed from the sample squared multiple correlations of *Y* with the *p*
predictors and of \\X_j\\ with the other \\p - 1\\, which treats the
sample standard deviations as constants (see
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md)). With
`noncentral = TRUE` the limits are the limits for the noncentrality
parameter at the observed *t* statistic times that standard error, and
the plan without `assurance` is the noncentral plan of
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
(`noncentral = TRUE`). With `noncentral = FALSE` the interval is the
estimate plus and minus *t* quantiles on \\N - p - 1\\ degrees of
freedom times that standard error (Kelley & Maxwell, 2003, Equation 1),
and the plan without `assurance` is theirs (Equation 2 with the *t*
quantile, at the expected squared multiple correlations), the central
plan of
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
at unit standard deviations, with the *t* quantiles of its two tails.

**The assurance step.** With `assurance` the plan is by default
(`assurance_method = "analytic"`) the published approximation, computed
without simulation: for the noncentral interval that of MBESS, for a
lower confidence limit of the squared multiple correlation, and for the
central interval the chi square approximation of Kelley and Maxwell
(2003, Equation 3), the plan of
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
without `assurance` with the variance of the coefficient multiplied by
the `assurance` quantile of a chi square on \\N - 1\\ degrees of freedom
over \\N - p - 1\\. Neither computes the probability it attains, so
`assurance_attained` is `NA`, and the printed footer says so; both start
from the plan without `assurance`, which can be a case or more smaller
than that of MBESS (see
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)).

With `assurance_method = "simulation"` the plan is the smallest sample
size, not below the plan without `assurance`, at which the interval for
the standardized coefficient is no wider than `width` with at least that
probability, for random predictors from a multivariate normal
population. The width of the noncentral interval depends on the observed
*t* statistic and on both sample squared multiple correlations, and that
of the central interval on the two squared multiple correlations, so for
either interval the probability is estimated from `G` samples drawn from
the exact sampling distribution of the sums of squares and products, the
same samples at every sample size the search tries, so that the estimate
changes smoothly with the sample size. The search assumes that the
estimate does not decrease as the sample size grows, which held in every
set checked (see "The assurance step" in
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md));
`assurance_attained` is the estimate at the returned sample size, at
least `assurance`, with its Monte Carlo standard error in the attribute
`assurance_mc_se`. The distribution of the width also depends on the
correlation of *Y* with the targeted predictor: on the covariance
structure path it is given, and on the squared multiple correlations
path the step uses the smallest probability over 21 equally spaced
values of it, from the smallest to the largest the three planning values
allow (the two ends with \\p = 2\\, the only values possible there), an
approximate minimum over that range. Give `seed` for a plan that can be
reproduced. This is the simulation Kelley and Maxwell (2008, pp.
186–187) described for the standardized coefficient, for which they
wrote that no satisfactory method then existed.

The central interval for the standardized coefficient is not the central
interval for the coefficient in the units of the data that
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
plans for with `noncentral = FALSE`, and the exact probability that
planner computes is that of the interval in the units of the data. The
standard error of that interval is formed from the residual sums of
squares of *Y* and of \\X_j\\; that of the standardized interval divides
each by the total sum of squares of its variable, which gives \\1 -
R^2\\ and \\1 - R^2\_{X_j \cdot X\_{-j}}\\, and these vary less from
sample to sample, so the standardized interval needs fewer cases to be
as narrow with the same probability. For two predictors correlated .3
with correlations of .4 and .3 with *Y*, a width of .3, and an assurance
of .90, this function returns 171 with
`assurance_method = "simulation"`, `seed = 1`, and the default `G`, at
which the central interval for the standardized coefficient is no wider
than .3 with probability .914 (.895 at 170). The plan is simulated, so
another seed can give one a case away: over seeds 1 to 30, 29 gave 171
and one gave 170. The central plan of
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
there is 187 (.997 for this interval), and the chi square approximation
of Kelley and Maxwell (2003, Equation 3, with the *t* quantile in place
of the normal one), the default plan here, is 181 (.986). For their five
exchangeable predictors (pp. 312–313), correlated .4 with one another
and .3 with *Y*, at a width of .2 and an assurance of .80, the simulated
plan is 478 (with each of seeds 1 to 30; .807, and .797 at 477), where
the approximation gives 486 (.877) and
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
491 (.911). Each probability is the proportion of intervals no wider
than the target among samples of raw multivariate normal data (50,000 or
100,000 samples at each *N*, Monte Carlo standard error at most .0014),
standardized and fit by least squares without DMAR. Kelley and Maxwell's
own simulation of that example (p. 317) reported 81.64 percent of the
intervals no wider than .2 at their plan of 485 (made with the normal
quantile) and a mean half-width of .098; raw multivariate normal data
give .871 there with the *t* quantile of their Equation 1 (.887 with the
normal one) and a mean half-width of .097 (50,000 samples, Monte Carlo
standard error at most .0015).

## Note

This function calls upon
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
in DMAR, which plans the coefficient in the units of the data by
default, but has a different naming scheme. See
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
for more details.

## Warning

As discussed in Kelley and Maxwell (2008), the sample size planning
approach from the AIPE perspective used in this function is only an
approximation.

## References

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

## See also

[`ss_aipe_src_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src_sensitivity.md),
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md),
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
[`ss_aipe_rc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Exchangeable correlation structure
true_cor_YX <- c(.3, .3, .3, .3, .3)
true_cor_XX <- rbind(c(1, .5, .5, .5, .5), c(.5, 1, .5, .5, .5), c(.5, .5, 1, .5, .5),
                c(.5, .5, .5, 1, .5), c(.5, .5, .5, .5, 1))

ss_aipe_src(width = .1, which_width = "Full", true_sd_Y = 1, true_sd_X = 1,
            true_cor_XX = true_cor_XX,
            true_cor_YX = true_cor_YX, which_predictor = 1, conf_level = 1 - .05)
#>  term        value
#>  necessary_N 2190 
#> 
#> Confidence level: 95%

# With an assurance of .85: the approximation with a lower confidence
# limit of the squared multiple correlation, and the plan found by
# simulation. G = 2000 keeps the example quick; a reported plan deserves
# the default G = 20000.
ss_aipe_src(width = .1, true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
            which_predictor = 1, assurance = .85)
#>  term               value
#>  necessary_N        2240 
#>  assurance          0.850
#>  assurance_attained <NA> 
#> 
#> Confidence level: 95%
#> 
#> Assurance method: analytic; the plan is the published approximation, which does not compute the probability that the interval is no wider than the width (assurance_attained is NA) and can attain less than the assurance. assurance_method = "simulation" computes it and returns the smallest plan that attains it
plan <- ss_aipe_src(width = .1, which_width = "Full", true_sd_Y = 1, true_sd_X = 1,
                    true_cor_XX = true_cor_XX,
                    true_cor_YX = true_cor_YX, which_predictor = 1, conf_level = 1 - .05,
                    assurance = .85, assurance_method = "simulation",
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

# The central t interval of Kelley and Maxwell (2003, Equation 1): the
# plan without assurance, and with an assurance of .85 by their chi
# square approximation and by simulation.
ss_aipe_src(width = .1, true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
            which_predictor = 1, noncentral = FALSE)
#>  term        value
#>  necessary_N 2184 
#> 
#> Confidence level: 95%
ss_aipe_src(width = .1, true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
            which_predictor = 1, noncentral = FALSE, assurance = .85)
#>  term               value
#>  necessary_N        2258 
#>  assurance          0.850
#>  assurance_attained <NA> 
#> 
#> Confidence level: 95%
#> 
#> Assurance method: analytic; the plan is the published approximation, which does not compute the probability that the interval is no wider than the width (assurance_attained is NA) and can attain less than the assurance. assurance_method = "simulation" computes it and returns the smallest plan that attains it
ss_aipe_src(width = .1, true_cor_XX = true_cor_XX, true_cor_YX = true_cor_YX,
            which_predictor = 1, noncentral = FALSE, assurance = .85,
            assurance_method = "simulation", G = 2000, seed = 113)
#>  term               value
#>  necessary_N        2250 
#>  assurance          0.850
#>  assurance_attained 0.851
#>  G                  2000 
#> 
#> Confidence level: 95%
#> 
#> Assurance method: simulation; assurance_attained is estimated from G = 2000 simulated samples (Monte Carlo SE 0.008)
```
