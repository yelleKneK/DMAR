# Mediation Analysis With Bootstrap Confidence Intervals

Estimates the simple mediation model, a predictor \\X\\ affecting an
outcome \\Y\\ directly and through a mediator \\M\\, and reports the
indirect, direct, and total effects with confidence intervals in one
tidy table. The indirect effect \\a b\\ is the product of the \\X \to
M\\ path and the \\M \to Y\\ path (holding \\X\\), and its sampling
distribution is skewed, which is why the default interval is the
percentile bootstrap rather than a normal approximation; the
bias-corrected and accelerated (BCa) bootstrap, the Monte Carlo
(parametric simulation) interval, and the Sobel normal-theory interval
are available for comparison. This is the analysis counterpart of the
planning functions
[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md)
and
[`ss_power_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md).

## Usage

``` r
mediate(
  data,
  x,
  m,
  y,
  covariates = NULL,
  ci_method = c("boot_percentile", "boot_bca", "monte_carlo", "sobel"),
  B = 2000,
  conf_level = 0.95,
  seed = NULL
)
```

## Arguments

- data:

  A `data.frame` containing the variables.

- x, m, y:

  Names (single character strings) of the predictor, the mediator, and
  the outcome columns in `data`. Each column must be numeric (a Date,
  date-time, or time difference column is read as its number, and a
  one-column matrix, such as
  [`scale()`](https://rdrr.io/r/base/scale.html) returns, as its
  column); a predictor with two categories is supplied coded 0 and 1.
  Any column name is accepted as it is.

- covariates:

  Optional character vector of covariate column names, entered in both
  the mediator and the outcome models. Numeric, logical, factor, and
  character columns are accepted; a factor or character covariate is
  coded once, on the complete cases, with R's contrasts for the session.
  A covariate with a single value in the complete cases adds nothing to
  either model and is left out with a message. A column may fill one
  role only.

- ci_method:

  Confidence interval method for the indirect effect:
  `"boot_percentile"` (default), `"boot_bca"`, `"monte_carlo"` (simulate
  \\a\\ and \\b\\ from their joint normal approximation; MacKinnon,
  Lockwood, & Williams, 2004), or `"sobel"` (the first-order
  normal-theory interval; reported for comparison, not recommended for
  inference). Direct and total effects always carry their ordinary
  *t*-based intervals.

- B:

  Number of bootstrap or Monte Carlo replications, a single whole
  number. Defaults to 2000, because each bootstrap replication refits
  both regressions. The fewest accepted depends on the interval and the
  level: 400 for the percentile and Monte Carlo intervals at the default
  `conf_level` (200 at 0.90, 2000 at 0.99), and 2000 for the BCa
  interval (1000 at 0.90, 10000 at 0.99, which the default does not
  meet). A smaller `B` is refused before anything is computed (class
  `dmar_count_invalid`); see *Bootstrap replications* under *Details*
  and
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).
  Not used, and not checked, for `ci_method = "sobel"`.

- conf_level:

  Confidence level. Defaults to 0.95.

- seed:

  Optional integer seed for the resampling, used locally (the caller's
  random number generator state is restored on exit). Default `NULL`
  draws the resamples from the current random number stream, which
  advances that stream.

## Value

A `data.frame` (class `dmar_tbl`) with rows `indirect_effect` (with its
`ci_method` interval), `direct_effect`, `total_effect`, the paths `a`
and `b`, their standard errors (`se_indirect` per Sobel, `se_a`,
`se_b`), the interval limits (`indirect_effect_ci_lower` /
`indirect_effect_ci_upper`, `direct_effect_ci_lower` /
`direct_effect_ci_upper`, `total_effect_ci_lower` /
`total_effect_ci_upper`), `proportion_mediated` (\\ab/c\\; `NA` when the
total effect is near zero, where the ratio is unstable), `N`, `B`, and
`B_used`, the number of replications the interval of the indirect effect
rests on (equal to `B` for the Monte Carlo interval and `NA`, as `B` is,
for the Sobel interval). The interval method is recorded in the
`"ci_method"` attribute and the confidence level in `"conf_level"`. For
the bootstrap and Monte Carlo intervals the attributes `"B"`
(requested), `"B_used"` (the number the interval rests on),
`"n_dropped"` (their difference), and `"B_required"` (the fewest
replications a request may ask for; that count is held to `B`, not to
`B_used`) hold the counts as integers, which is where
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)
reads them, and `"bootstrap_status"` is one of `"complete"`,
`"dropped"`, `"exploratory"`, `"insufficient_count"`, or
`"insufficient_share"` (the last two mean that the limits of the
indirect effect are `NA`). For the two bootstrap intervals
`"n_covariate_redundant"` is the number of resamples kept in which a
covariate column was redundant (see *Details*).

## Details

The model is the standard pair of regressions \$\$M = i_M + a X +
\mathbf{g}'\mathbf{C} + e_M, \qquad Y = i_Y + c' X + b M +
\mathbf{h}'\mathbf{C} + e_Y,\$\$ with \\\mathbf{C}\\ the optional
covariates. The indirect effect is \\a b\\, the direct effect \\c'\\,
and the total effect \\c = c' + a b\\ (an identity in linear models with
the same cases, which the implementation exploits as an internal
consistency check).

Bootstrap intervals resample cases (rows) with replacement `B` times,
refitting both regressions in each resample (Efron & Tibshirani, 1993).
`"boot_percentile"` takes the interval limits from the empirical
quantiles of the bootstrapped \\a b\\ estimates; it is not forced to be
symmetric about the estimate, which is the point for a skewed sampling
distribution. `"boot_bca"` (the bias-corrected and accelerated interval)
additionally adjusts the two quantile positions for median bias,
estimated from the bootstrap distribution, and for the rate at which the
variance of the estimator changes with the parameter, the acceleration,
estimated by the jackknife (which adds *N* extra pairs of fits); the
adjustments make it second-order accurate where the percentile interval
is first-order accurate (DiCiccio & Efron, 1996). Because the adjusted
quantile positions sit farther into the tails of the bootstrap
distribution, the BCa interval benefits more than the percentile
interval does from a `B` well above the default.

**Bootstrap replications.** Each resample refits the two regressions by
least squares on the rows drawn, with the covariates coded as they are
in the data, and the interval is read from the quantiles of the `B`
values of \\a b\\. A limit read from a quantile is located by the
replications that lie beyond it, so the percentile and Monte Carlo
intervals are run only when `B` puts at least 10 replications beyond
each limit (400 at the 95 percent level), and the BCa interval, whose
limits move farther into the tails, only with 50 beyond each nominal
limit (2000 at the 95 percent level). A resample that leaves a covariate
column redundant (a category with no case in the resample) still gives
unique values of \\a\\ and \\b\\ and is kept; the
`"n_covariate_redundant"` attribute counts these. A resample in which
\\a\\ or \\b\\ is not estimable (the predictor constant, or the
predictor or the mediator a linear combination of the other predictors
there) has no indirect effect. Such a replication is dropped and
counted, never replaced, and the interval is computed from the rest. The
result records `B`, `B_used`, `n_dropped`, and `B_required` as
attributes and `B` and `B_used` as rows. One warning of class
`dmar_bootstrap_replications_dropped` states the counts when any
replication was dropped, with a caution when the dropped replications
are at least the share of `B` beyond one limit (2.5 percent at the 95
percent level), because they could then decide by themselves where that
limit falls. When more than half were dropped, no interval is reported:
the two limits of the indirect effect are `NA`, every other row is as it
would have been, and the warning has class
`dmar_bootstrap_insufficient`. An interval from the replications that
remain describes resamples in which the effect is defined. Dropping
these replications is the rule DMAR applies to a replication that
returns no value; it is a convention of the package, not a result that
gives the interval its nominal coverage. Many dropped replications mean
the data carry little information about a path (a predictor that varies
in few cases), a design in which the case bootstrap is itself not
reliable.

The BCa interval needs the indirect effect with each case left out; when
one of those is not estimable, or the bias correction or the
acceleration is not finite, the BCa interval is refused (class
`dmar_bca_undefined`) and the percentile interval remains available.
After the replications are run, the count of 10 is applied again at the
two positions where the bias correction and the acceleration place the
limits. When fewer than 10 replications lie beyond one of them, no BCa
interval is reported (limits `NA`, warning of class
`dmar_bootstrap_insufficient`, which carries the two positions as
`positions` and a suggested `B` as `B_suggest`).

`options(dmar.bootstrap_exploratory = TRUE)` lets a smaller `B` (2 or
more) run for a first look or a demonstration. Every such call gives a
warning of class `dmar_bootstrap_exploratory`, and the attribute
`"bootstrap_status"` is `"exploratory"`; the interval is not one to
report.

The Sobel standard error is computed by
[`var_indirect_effect`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md).
The proportion mediated is one of the effect size measures for mediation
models surveyed by Preacher and Kelley (2011); its instability when the
total effect is small is why it is reported as `NA` near a zero total
effect. Listwise deletion is applied to the analysis variables; for full
information maximum likelihood under missingness, fit the model in
lavaan (see [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md)
for the package's FIML front end philosophy). The model has one mediator
and no interaction, the cases are treated as independent, and incomplete
cases are removed; for several mediators or missing data see
[`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md).

Mediation language implies causal structure: with observational data the
estimates are conditional associations, and the causal reading requires
the usual no-unmeasured-confounding assumptions for both the \\X \to M\\
and \\M \to Y\\ links (MacKinnon, 2008). The function computes; the
design earns the interpretation.

## References

DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
*Statistical Science, 11*(3), 189–228.

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

MacKinnon, D. P. (2008). *Introduction to statistical mediation
analysis*. Erlbaum.

MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence
limits for the indirect effect: Distribution of the product and
resampling methods. *Multivariate Behavioral Research, 39*(1), 99–128.
[doi:10.1207/s15327906mbr3901_4](https://doi.org/10.1207/s15327906mbr3901_4)

Preacher, K. J., & Hayes, A. F. (2008). Asymptotic and resampling
strategies for assessing and comparing indirect effects in multiple
mediator models. *Behavior Research Methods, 40*(3), 879–891.
[doi:10.3758/BRM.40.3.879](https://doi.org/10.3758/BRM.40.3.879)

Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
models: Quantitative strategies for communicating indirect effects.
*Psychological Methods, 16*(2), 93–115.
[doi:10.1037/a0022658](https://doi.org/10.1037/a0022658)

## See also

[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md)
and
[`ss_power_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md)
for planning the study this function analyzes;
[`var_indirect_effect`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md)
for the Sobel variance.

Other mediation:
[`mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md),
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Simulated mediation: X raises M (a = .5), M raises Y (b = .4), and a
# little direct effect remains (c' = .2).
set.seed(113)
n <- 200
x <- rnorm(n)
m <- 0.5 * x + rnorm(n, 0, sqrt(1 - 0.25))
y <- 0.2 * x + 0.4 * m + rnorm(n, 0, 0.8)
d <- data.frame(x = x, m = m, y = y)

# The percentile bootstrap is the default and is what a reported
# analysis would use. Each resample refits both regressions.
# B = 400 is the fewest replications a 95 percent percentile interval
# is reported from and keeps the example quick; a reported interval
# deserves the default B = 2000 or more.
mediate(d, x = "x", m = "m", y = "y", B = 400, seed = 113)
#>  term                     value
#>  indirect_effect          0.182
#>  indirect_effect_ci_lower 0.116
#>  indirect_effect_ci_upper 0.259
#>  direct_effect            0.194
#>  direct_effect_ci_lower   0.063
#>  direct_effect_ci_upper   0.324
#>  total_effect             0.375
#>  total_effect_ci_lower    0.250
#>  total_effect_ci_upper    0.501
#>  a                        0.449
#>  b                        0.404
#>  se_a                     0.062
#>  se_b                     0.068
#>  se_indirect              0.039
#>  proportion_mediated      0.484
#>  N                        200  
#>  B                        400  
#>  B_used                   400  
#> 
#> Confidence level: 95%
#> 
#> Bootstrap: 400 replications, 400 used (0 dropped; at least 400 must be requested)

# The Monte Carlo interval draws a and b from their joint normal
# approximation instead of refitting, so it stays quick even at a
# large B. Its limits sit close to the bootstrap limits here.
mediate(d, x = "x", m = "m", y = "y", ci_method = "monte_carlo",
        B = 10000, seed = 113)
#>  term                     value
#>  indirect_effect          0.182
#>  indirect_effect_ci_lower 0.109
#>  indirect_effect_ci_upper 0.266
#>  direct_effect            0.194
#>  direct_effect_ci_lower   0.063
#>  direct_effect_ci_upper   0.324
#>  total_effect             0.375
#>  total_effect_ci_lower    0.250
#>  total_effect_ci_upper    0.501
#>  a                        0.449
#>  b                        0.404
#>  se_a                     0.062
#>  se_b                     0.068
#>  se_indirect              0.039
#>  proportion_mediated      0.484
#>  N                        200  
#>  B                        10000
#>  B_used                   10000
#> 
#> Confidence level: 95%
#> 
#> Bootstrap: 10000 replications, 10000 used (0 dropped; at least 400 must be requested)

# The Sobel interval is closed form. It assumes the product ab is
# normally distributed, which is why it is reported for comparison
# rather than used for inference; its B and B_used rows are NA
# because no replications are drawn.
mediate(d, x = "x", m = "m", y = "y", ci_method = "sobel")
#>  term                     value
#>  indirect_effect          0.182
#>  indirect_effect_ci_lower 0.105
#>  indirect_effect_ci_upper 0.259
#>  direct_effect            0.194
#>  direct_effect_ci_lower   0.063
#>  direct_effect_ci_upper   0.324
#>  total_effect             0.375
#>  total_effect_ci_lower    0.250
#>  total_effect_ci_upper    0.501
#>  a                        0.449
#>  b                        0.404
#>  se_a                     0.062
#>  se_b                     0.068
#>  se_indirect              0.039
#>  proportion_mediated      0.484
#>  N                        200  
#>  B                        <NA> 
#>  B_used                   <NA> 
#> 
#> Confidence level: 95%
```
