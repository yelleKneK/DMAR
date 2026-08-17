# Power and Precision for the One-Way ANOVA: A Model Comparison Perspective

This vignette plans, estimates, and reasons about a one-way
between-subjects analysis of variance using the framework of Maxwell,
Delaney, and Kelley (2027): the comparison of a full model against a
restricted model. The same logic that produces the omnibus *F* test
produces its power, its effect size, the confidence interval for that
effect size, the power of a focused contrast, and the sample size needed
for either a decision or a precise estimate. Each quantity is first
worked by hand, so the arithmetic is visible, and then obtained from the
DMAR function that encapsulates it, with a note on what that function is
doing internally.

## Two Models, One Question

Every test in the model comparison perspective asks one question: how
much worse does the data fit once the null hypothesis is imposed? Write
a full model that lets each group have its own mean,

``` math
Y_{ij} = \mu_j + \varepsilon_{ij},
```

and a restricted model that forces a single common mean,

``` math
Y_{ij} = \mu + \varepsilon_{ij}.
```

The omnibus null hypothesis is the restriction
$`\mu_1 = \mu_2 = \cdots =
\mu_a`$, which collapses every group onto the grand mean. The treatment
effects $`\alpha_j = \mu_j - \mu`$ are exactly the per-observation
discrepancies between the two models’ predictions. With $`E_F`$ and
$`E_R`$ the error sums of squares of the full and restricted models, the
test statistic always has the same shape,

``` math
F = \frac{(E_R - E_F)/(df_R - df_F)}{E_F/df_F},
```

and the effect size is the proportional reduction in error,
$`\mathrm{PRE} = (E_R - E_F)/E_R`$, which in the one-way design is the
same number as eta squared. Power, estimation, and sample size planning
all read off these quantities, so design and analysis speak one
language.

## The Four-Group Running Example

A researcher plans a four-condition study. Under the alternative
hypothesis the population means are 90, 92, 88, and 81, homogeneity of
variance is assumed with a common within-group variance of 144 (a
within-group standard deviation of 12), and the plan calls for 20
participants per group.

``` r

mu_j          <- c(90, 92, 88, 81)
sigma_squared <- 144
sigma         <- sqrt(sigma_squared)
n_j           <- 20

a <- length(mu_j)
N <- n_j * a

mu      <- mean(mu_j)        # grand mean: the restricted model's one prediction
alpha_j <- mu_j - mu         # treatment effects: full minus restricted
alpha_j
#> [1]  2.25  4.25  0.25 -6.75
```

The reduction in error from the restricted to the full model is the
between-groups sum of squares, which is the sum over all *N*
observations of the squared discrepancy between the two models’
predictions. In error-variance units this is the noncentrality parameter
that drives power:

``` math
\lambda = \frac{\sum_j n_j\,\alpha_j^2}{\sigma^2}.
```

``` r

reduction <- sum(n_j * alpha_j^2)   # population between-groups sum of squares
lambda    <- reduction / sigma_squared
lambda
#> [1] 9.548611
```

## The Effect Size as a Proportional Reduction in Error

Cohen’s *f* is a population quantity, $`f = \sigma_m / \sigma`$, where
$`\sigma_m`$ is the standard deviation of the population means. The one
subtlety that has caused decades of confusion is the divisor of that
variance. The variance of the means here is a property of a *fixed* set
of parameters, so it divides by *a*, the number of groups, not by
$`a - 1`$:

``` math
\sigma_m^2 = \frac{\sum_j (\mu_j - \mu)^2}{a} = \frac{\sum_j \alpha_j^2}{a}.
```

By hand:

``` r

sigma_squared_m <- sum(alpha_j^2) / a   # divisor a, not a - 1
sigma_m         <- sqrt(sigma_squared_m)
f_by_hand       <- sigma_m / sigma
f_by_hand
#> [1] 0.3454817
```

[`cohen_f()`](https://yelleknek.github.io/DMAR/reference/cohen_f.md)
does this directly from the means and the within-group variance. It uses
the population divisor internally, so the user never has to remember
whether to divide by *a* or by $`a - 1`$:

``` r

cohen_f(mu = mu_j, sigma_squared = sigma_squared)
```

| term    | value |
|:--------|:------|
| cohen_f | 0.345 |

``` r

f <- cohen_f(mu = mu_j, sigma_squared = sigma_squared)$value
```

The same effect, expressed as the population proportional reduction in
error (the population eta squared), is
$`\mathrm{PRE} = f^2 / (1 + f^2)`$. The noncentrality, the proportional
reduction in error, and *f* are three views of one quantity, and the
`convert_*` family moves among them. Here the noncentrality maps to the
population PRE and back:

``` r

PRE <- f^2 / (1 + f^2)
PRE
#> [1] 0.1066305

convert_lambda_R2(lambda = lambda, N = N)   # lambda to population PRE (= R^2)
```

| term      | value |
|:----------|:------|
| lambda_r2 | 0.107 |

``` r

convert_R2_lambda(R2 = PRE,    N = N)       # and back to the noncentrality
```

| term      | value |
|:----------|:------|
| r2_lambda | 9.55  |

Off-the-shelf benchmarks for *f* are a last resort for when nothing is
known about the design; a set of hypothesized means, as used here, is
always the better input, because it ties the effect size to the actual
scientific quantities rather than to a convention.

## Power for the Omnibus Test

The comparison has numerator degrees of freedom
$`df_R - df_F = (N-1) - (N-a) =
a - 1`$ and denominator degrees of freedom $`df_F = N - a`$. Power is
the area of the noncentral *F* distribution, centered at $`\lambda`$,
that lies beyond the critical value taken from the central *F*. By hand:

``` r

df_numerator   <- a - 1
df_denominator <- N - a

critical_F <- qf(1 - .05, df1 = df_numerator, df2 = df_denominator)
pf(critical_F, df1 = df_numerator, df2 = df_denominator,
   ncp = lambda, lower.tail = FALSE)
#> [1] 0.7149725
```

[`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md)
returns the same value. Supplied with a total sample size *N* it reports
the realized power; internally it forms $`\lambda = N f^2`$, sets the
two degrees of freedom, and integrates the noncentral *F* beyond the
central critical value exactly as above:

``` r

ss_power_one_way_anova(a = a, f = f, N = N, alpha_level = .05)
```

| term          | value |
|:--------------|:------|
| specified_N   | 80    |
| a             | 4     |
| noncentrality | 9.55  |
| actual_power  | 0.715 |

## Planning the Sample Size

To plan rather than evaluate, hold the effect fixed and solve for the
sample size that reaches a target power. The function searches over *N*,
reporting the necessary total sample size and the per-group *n*:

``` r

ss_power_one_way_anova(a = a, f = f, desired_power = .80, alpha_level = .05)
```

| term          | value |
|:--------------|:------|
| necessary_N   | 96    |
| n_per_group   | 24    |
| a             | 4     |
| noncentrality | 11.5  |
| actual_power  | 0.803 |

## The Model Comparison, Made Literal

Everything above is a population calculation. To see the model
comparison on actual data, simulate one data set under these parameters
with
[`simulate_anova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_anova_data.md),
fit the restricted and full models by maximum likelihood with
[`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md), and
compare them with the likelihood ratio test. The likelihood ratio test
is the sample analogue of $`E_R - E_F`$: it asks whether freeing the
group means improves the fit by more than sampling noise would explain.

``` r

d <- simulate_anova_data(mu = mu_j, sigma = sigma, a = a, n = n_j)
head(d)
#>   group         y
#> 1     1  91.60025
#> 2     1 106.50266
#> 3     1  98.98458
#> 4     1  74.47379
#> 5     1  83.29475
#> 6     1  69.21061
```

``` r

restricted <- mlmr(y ~ 1,     data = d, ci_method = "wald", effect_sizes = FALSE)
full       <- mlmr(y ~ group, data = d, ci_method = "wald", effect_sizes = FALSE)

# Likelihood ratio test of the nested restricted model against the full model
anova(restricted, full)
#> Likelihood ratio test for nested mlmr fits
#> Model 2: y ~ group
#> Model 1: y ~ 1
#>         Df    AIC    BIC  Chisq Chisq diff Df diff Pr(>Chisq)   
#> Model 2  0 898.46 931.81  0.000                                 
#> Model 1  3 906.97 933.17 14.511     14.511       3   0.002286 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# The sample proportional reduction in error for this same comparison
full$R2
#> [1] 0.1658905
```

The omnibus *F* from the ordinary least squares fit tells the same
story, and its components are the sample $`E_R - E_F`$, $`E_F`$, and the
degrees of freedom of the comparison:

``` r

fit_aov <- aov(y ~ group, data = d)
print_anova(anova(fit_aov))
#> Analysis of Variance Table
#> 
#> Response: y
#> 
#>           Df    Sum Sq  Mean Sq F value Pr(>F)
#> group      3  2543.688 847.8960 5.03838 0.0031
#> Residuals 76 12789.843 168.2874      NA   <NA>
```

## Estimating the Effect Size With a Confidence Interval

Planning supposes an effect; analysis estimates one and states its
uncertainty. Because the effect size is a proportional reduction in
error, the natural deliverable is the estimate together with a
confidence interval, obtained by inverting the noncentral *F*
distribution (Kelley, 2007; Steiger, 2004). Given an observed *F*, DMAR
finds the noncentrality values whose noncentral *F* distributions place
the observed *F* at the chosen lower and upper tail probabilities (this
is what
[`conf_limits_ncf()`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
does), then transforms those noncentrality endpoints into the metric
requested: a proportion of variance, a signal to noise ratio, or its
square root.

The canonical illustration in the literature is Bargman’s (1970)
five-group analysis with 11 per group and an observed *F* of 11.221
(also used by Steiger, 2004). The point estimate of omega squared, a
less biased estimator of the population proportion of variance than eta
squared, travels with its interval:

``` r

os_est <- omega_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                        N = 55)
os_est
#>    effect omega_squared F_value df_effect df_error  N
#> 1 overall     0.4263902  11.221         4       50 55
os_ci <- ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50,
                          N = 55)
os_ci
```

| effect  | omega_squared | lower_limit | upper_limit | F_value | df_effect | df_error | N   |
|:--------|:--------------|:------------|:------------|:--------|:----------|:---------|:----|
| overall | 0.426         | 0.226       | 0.587       | 11.2    | 4         | 50       | 55  |

The population proportion of variance accounted for is estimated at
omega squared = 0.43, with a 95% confidence interval of \[0.23, 0.59\]:
a substantial effect whose width still admits values from roughly a
quarter to three fifths of the variance, a reminder that even a clearly
nonzero effect is estimated with real uncertainty at this sample size.

The same comparison, expressed as a proportion of variance, as the
signal to noise ratio $`f^2`$, and as the square root of the signal to
noise ratio (which is Cohen’s *f*):

``` r

ci_pvaf(F_value  = 11.221, df_1 = 4, df_2 = 50, N = 55)
```

| term            | value | prob_less | prob_greater |
|:----------------|:------|:----------|:-------------|
| lower_limit     | 0.226 | 0.025     | 0.975        |
| pvaf            | 0.473 | NA        | NA           |
| upper_limit     | 0.587 | 0.975     | 0.025        |
| actual_coverage | 0.95  | NA        | NA           |

Confidence level: 95%

``` r

ci_snr(F_value   = 11.221, df_1 = 4, df_2 = 50, N = 55)
```

| term        | value |
|:------------|:------|
| lower_limit | 0.293 |
| upper_limit | 1.42  |

Confidence level: 95%

``` r

ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
```

| term        | value |
|:------------|:------|
| lower_limit | 0.541 |
| upper_limit | 1.19  |

Confidence level: 95%

[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md)
also accepts a design-stage specification. Given the planning means, the
within-group variance, and the per-group sample size, it derives the
implied *F* internally and returns the interval on Cohen’s *f* that the
planned design would produce. For the four-group running example:

``` r

ci_srsnr(means = mu_j, sigma_squared = sigma_squared, n_per_group = n_j)
#> Warning: The observed F_value is below the alpha_lower critical value of the
#> central F-distribution, so the lower confidence limit on the square root of the
#> signal-to-noise ratio is 0.
```

| term        | value |
|:------------|:------|
| lower_limit | 0     |
| upper_limit | 0.542 |

Confidence level: 95%

## Carrying Uncertainty Into the Next Study

A pilot estimate is uncertain, and planning the next study around a
point estimate quietly assumes that uncertainty away. A more defensible
habit is to plan around a conservative bound of the effect. Take a
one-sided 80 percent lower limit for *f* from the Bargman analysis by
putting all of the Type I error in the lower tail, then plan the next
five-group study to that bound:

``` r

lower_f <- ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55,
                    alpha_lower = .20, alpha_upper = 0)
lower_f
```

| term        | value |
|:------------|:------|
| lower_limit | 0.727 |
| upper_limit | Inf   |

Confidence level: 80%

``` r


f_lower <- lower_f$value[lower_f$term == "lower_limit"]
ss_power_one_way_anova(a = 5, f = f_lower, desired_power = .80)
```

| term          | value |
|:--------------|:------|
| necessary_N   | 30    |
| n_per_group   | 6     |
| a             | 5     |
| noncentrality | 15.9  |
| actual_power  | 0.842 |

## A Second Example: Three Groups

The worked example in Maxwell, Delaney, and Kelley (2027) uses three
groups with population means 400, 450, and 500, a within-group standard
deviation of 100, and 10 participants per group. The same spine applies,
now stated compactly with DMAR.

``` r

mu_j_2          <- c(400, 450, 500)
sigma_2         <- 100
sigma_squared_2 <- sigma_2^2
n_j_2           <- 10
a_2             <- length(mu_j_2)
N_2             <- n_j_2 * a_2

# Noncentrality from the model discrepancy, by hand
lambda_2 <- sum(n_j_2 * (mu_j_2 - mean(mu_j_2))^2) / sigma_squared_2
lambda_2
#> [1] 5

# Effect size and realized power at the planned design, from DMAR
f_2 <- cohen_f(mu = mu_j_2, sigma_squared = sigma_squared_2)$value
ss_power_one_way_anova(a = a_2, f = f_2, N = N_2, alpha_level = .05)
```

| term          | value |
|:--------------|:------|
| specified_N   | 30    |
| a             | 3     |
| noncentrality | 5     |
| actual_power  | 0.458 |

``` r


# Sample size needed for power .80
plan2 <- ss_power_one_way_anova(a = a_2, f = f_2, desired_power = .80, alpha_level = .05)
```

The realized power at 10 per group is modest. The sample size needed for
power .80 is about 21 per group, the value the textbook reaches:

``` r

plan2
```

| term          | value |
|:--------------|:------|
| necessary_N   | 63    |
| n_per_group   | 21    |
| a             | 3     |
| noncentrality | 10.5  |
| actual_power  | 0.815 |

## A Planned Contrast as a Model Comparison

The omnibus test compares one mean against *a* means. A focused question
imposes a different restriction: that a single contrast equals zero.
Returning to the four-group example, suppose the substantive question is
whether the average of the first three groups differs from the fourth.
The full model lets the group means be free, so the contrast
$`\psi = \sum_j c_j \mu_j`$ may be nonzero; the restricted model
constrains the means so that $`\psi = 0`$. This is a
one-degree-of-freedom comparison.

The contrast weights must satisfy three conditions, all of which DMAR’s
contrast functions expect for the standardized fractional form. The
weights must sum to zero, which is what makes them a contrast; and to
put $`\psi`$ on the metric of a difference between two composite means,
the positive weights must sum to 1 and the negative weights must sum to
$`-1`$:

``` r

c_weights <- c(1/3, 1/3, 1/3, -1)
sum(c_weights)                  # all weights sum to zero
#> [1] -5.551115e-17
sum(c_weights[c_weights > 0])   # the positive weights sum to 1
#> [1] 1
sum(c_weights[c_weights < 0])   # the negative weights sum to -1
#> [1] -1
```

The population value of the contrast is the discrepancy the restriction
would erase. The single-degree-of-freedom reduction in error is
$`\psi^2 / \sum_j (c_j^2 / n_j)`$, and dividing by the error variance
gives the noncentrality:

``` r

psi <- sum(c_weights * mu_j)
psi
#> [1] 9

reduction_contrast <- psi^2 / sum(c_weights^2 / n_j)
lambda_contrast    <- reduction_contrast / sigma_squared
lambda_contrast
#> [1] 8.4375

# Power by hand, from the noncentral F with one numerator degree of freedom
critical_F_c <- qf(1 - .05, df1 = 1, df2 = N - a)
pf(critical_F_c, df1 = 1, df2 = N - a, ncp = lambda_contrast, lower.tail = FALSE)
#> [1] 0.8180309
```

[`ss_power_contrast()`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md)
returns the same power and will also solve for the sample size. One
point deserves attention, because the sample size argument is not
uniform across the planning family:
[`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md)
takes the *total* sample size, whereas
[`ss_power_contrast()`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md)
takes the *per-group* size in `n_per_group`. Here 20 per group
corresponds to a total of 80:

``` r

# Realized power at 20 per group
ss_power_contrast(c_weights = c_weights, mu = mu_j,
                  sigma_squared = sigma_squared, n_per_group = n_j)
```

| term                  | value |
|:----------------------|:------|
| specified_n_per_group | 20    |
| total_N               | 80    |
| actual_power          | 0.818 |
| noncentral_t_parm     | 2.9   |
| effect_size_f         | 0.325 |

``` r


# Per-group sample size needed for power .90
ss_power_contrast(c_weights = c_weights, mu = mu_j,
                  sigma_squared = sigma_squared, desired_power = .90)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 26    |
| total_N               | 104   |
| actual_power          | 0.907 |
| noncentral_t_parm     | 3.31  |
| effect_size_f         | 0.325 |

The analysis-stage counterpart estimates the contrast with a confidence
interval. [`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md)
forms $`\hat\psi \pm t_{1-\alpha/2,\,N-a}\,\mathrm{SE}`$, where the
standard error is $`s\sqrt{\sum_j c_j^2 / n_j}`$ and *s* is the root
mean square error (the ANOVA standard deviation, not the variance):

``` r

ci_c(means = mu_j, s_anova = sigma, c_weights = c_weights,
     n = n_j, N = N, conf_level = .95)
```

| term        | value |
|:------------|:------|
| lower_limit | 2.83  |
| contrast    | 9     |
| upper_limit | 15.2  |

Confidence level: 95%

From fitted data,
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)
carries the same estimate, its standard error, the *t* statistic, the
degrees of freedom, and a confidence interval, with optional
multiple-comparison adjustments:

``` r

contrast_test(fit_aov,
              contrasts = list("avg of 1, 2, 3 vs 4" = c(1/3, 1/3, 1/3, -1)))
```

| contrast            | estimate | se   | t    | df  | p_value | p_adjusted | ci_lower | ci_upper |
|:--------------------|:---------|:-----|:-----|:----|:--------|:-----------|:---------|:---------|
| avg of 1, 2, 3 vs 4 | 12.2     | 3.35 | 3.64 | 76  | 0.0005  | 0.0005     | 5.52     | 18.9     |

Confidence level: 95%

## Beyond the Reject or Retain Decision

Power is the probability of a correct rejection, but a significant
result is not the end of the story. Following Gelman and Carlin (2014),
[`design_consequences()`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
reports, for the planned contrast, the power together with the Type S
error (the probability that a significant estimate has the wrong sign)
and the Type M error (the factor by which a significant estimate
exaggerates the true effect), and the expected width of the confidence
interval. Feed it the population contrast and its standard error:

``` r

se_psi <- sqrt(sigma_squared * sum(c_weights^2 / n_j))
design_consequences(true_effect = psi, se = se_psi, df = N - a, alpha_level = .05)
```

| term                | value   |
|:--------------------|:--------|
| power               | 0.818   |
| type_s_error        | 8.4e-07 |
| exaggeration_ratio  | 1.11    |
| expected_half_width | 6.15    |
| mean_ci_width       | 12.3    |
| median_ci_width     | 12.3    |
| sd_ci_width         | 0.999   |
| pct_ci_less_w       | NA      |
| target_width        | NA      |
| true_effect         | 9       |
| se                  | 3.1     |
| df                  | 76      |
| alpha_level         | 0.05    |

Confidence level: 95%

A well-powered design has a negligible Type S error and a Type M factor
near 1, meaning a significant estimate neither flips sign nor materially
overstates the effect. Underpowered designs fail on both counts even
when they occasionally reach significance.

## Planning a Contrast for Precision

Power asks whether the test will reject; accuracy in parameter
estimation (AIPE) asks how precisely the contrast will be estimated.
Rather than a target power, the input is a target confidence interval
width, and the planner returns the sample size that achieves it. Here
the goal is an expected 95 percent interval no wider than 6 points, with
90 percent assurance that the realized interval meets that target:

``` r

ss_aipe_c(error_variance = sigma_squared, c_weights = c_weights,
          width = 6, conf_level = .95, assurance = .90)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 91    |

Confidence level: 95%

## Why Not base R’s `power.anova.test()`

The base function
[`power.anova.test()`](https://rdrr.io/r/stats/power.anova.test.html)
takes between- and within-group variances, but its `between.var` is the
*sample* variance of the means (with divisor $`a - 1`$), whereas Cohen’s
*f* and the model comparison use the *population* variance of the means
(with divisor *a*). Passing the population value understates the effect;
the classic fix is to pass `var(mu_j)` instead. DMAR removes the hazard,
because
[`cohen_f()`](https://yelleknek.github.io/DMAR/reference/cohen_f.md)
uses the correct divisor internally:

``` r

sigma_squared_m              # population variance of the means (divisor a)
#> [1] 17.1875
var(mu_j)                    # sample version (divisor a - 1)
#> [1] 22.91667
var(mu_j) / sigma_squared_m  # the ratio is a / (a - 1)
#> [1] 1.333333

cohen_f(mu = mu_j, sigma_squared = sigma_squared)$value
#> [1] 0.3454817
sqrt(sigma_squared_m / sigma_squared)   # identical, computed by hand
#> [1] 0.3454817
```

The chain
[`cohen_f()`](https://yelleknek.github.io/DMAR/reference/cohen_f.md) to
[`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md)
is therefore both correct and unambiguous, which is why this vignette
never reaches for the base function.

## References

Bargman, R. E. (1970). Interpretation and use of a generalized
discriminant function. In R. C. Bose et al. (Eds.), *Essays in
probability and statistics*. University of North Carolina Press.

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Erlbaum.

Gelman, A., & Carlin, J. (2014). Beyond power calculations: Assessing
Type S (sign) and Type M (magnitude) errors. *Perspectives on
Psychological Science, 9*(6), 641–651.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. See Chapter 3 on the one-way design and the analysis of
group differences, and Chapters 4 through 6 on contrasts and focused
comparisons.

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
