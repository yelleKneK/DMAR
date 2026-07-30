# Standardized Mean Difference

Estimates the standardized mean difference (Cohen's *d*), the difference
between two group means divided by the pooled standard deviation, from
either raw data or summary statistics. Expressing the difference in
standard deviation units frees the comparison from the raw measurement
units, so effects can be compared across measures and studies; either
the ordinary or the unbiased (Hedges, 1981) estimate can be returned.

## Usage

``` r
smd(
  group_1 = NULL,
  group_2 = NULL,
  mean_1 = NULL,
  mean_2 = NULL,
  s_1 = NULL,
  s_2 = NULL,
  s = NULL,
  n_1 = NULL,
  n_2 = NULL,
  unbiased = FALSE
)
```

## Arguments

- group_1:

  Raw data for group 1

- group_2:

  Raw data for group 2

- mean_1:

  The mean of group 1

- mean_2:

  The mean of group 2

- s_1:

  The standard deviation of group 1 (i.e., the square root of the
  unbiased estimator of the population variance)

- s_2:

  The standard deviation of group 2 (i.e., the square root of the
  unbiased estimator of the population variance)

- s:

  The pooled group standard deviation (i.e., the square root of the
  unbiased estimator of the population variance)

- n_1:

  The sample size within group 1

- n_2:

  The sample size within group 2

- unbiased:

  Returns the unbiased estimate of the standardized mean difference

## Value

A 1-row `data.frame` with columns `term` (`"smd"`) and `value` (the
estimated standardized mean difference).

## Details

When `unbiased=TRUE`, the unbiased estimate of the standardized mean
difference is returned (Hedges, 1981).

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and
noncentral distributions. *Educational and Psychological Measurement,
61*(4), 532–574.
[doi:10.1177/0013164401614002](https://doi.org/10.1177/0013164401614002)

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128.

Kelley, K. (2005) The effects of nonnormal distributions on confidence
intervals around the standardized mean difference: Bootstrap and
parametric confidence intervals, *Educational and Psychological
Measurement, 65*, 51–69.
[doi:10.1177/0013164404264850](https://doi.org/10.1177/0013164404264850)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons and Chapter 3
on one-way ANOVA.)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`smd_c`](https://yelleknek.github.io/DMAR/reference/smd_c.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md),
[`ss_power_smd`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`plot_smd`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Generate sample data.
set.seed(113)
g.1 <- rnorm(n = 25, mean = .5, sd = 1)
g.2 <- rnorm(n = 25, mean = 0, sd = 1)
smd(group_1 = g.1, group_2 = g.2)
#>  term value
#>  smd  0.399

M.x <- .66745
M.y <- .24878
sd <- 1.048
smd(mean_1 = M.x, mean_2 = M.y, s = sd)
#>  term value
#>  smd  0.399

M.x <- .66745
M.y <- .24878
n1 <- 25
n2 <- 25
sd.1 <- .95817
sd.2 <- 1.1311
smd(mean_1 = M.x, mean_2 = M.y, s_1 = sd.1, s_2 = sd.2, n_1 = n1, n_2 = n2)
#>  term value
#>  smd  0.399

smd(mean_1 = M.x, mean_2 = M.y, s_1 = sd.1, s_2 = sd.2, n_1 = n1, n_2 = n2,
    unbiased = TRUE)
#>  term value
#>  smd  0.393
```
