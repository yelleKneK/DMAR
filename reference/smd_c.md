# Standardized Mean Difference Using the Control Group as the Basis of Standardization

Estimates the standardized mean difference using the control group
standard deviation as the basis of standardization (Glass's *g*), from
either raw data or summary statistics, in ordinary or unbiased form.
Standardizing by the control group alone keeps the scale of the effect
free of any treatment effect on variability.

## Usage

``` r
smd_c(
  group_T = NULL,
  group_C = NULL,
  mean_T = NULL,
  mean_C = NULL,
  s_C = NULL,
  n_C = NULL,
  unbiased = FALSE
)
```

## Arguments

- group_T:

  Raw data for the treatment group

- group_C:

  Raw data for the control group

- mean_T:

  The mean of the treatment group

- mean_C:

  The mean of the control group

- s_C:

  The standard deviation of the control group (i.e., the square root of
  the unbiased estimator of the population variance)

- n_C:

  The sample size of the control group

- unbiased:

  Returns the unbiased estimate of the standardized mean difference
  using the standard deviation of the control group

## Value

A 1-row `data.frame` with columns `term` (`"smd_c"`) and `value` (the
estimated standardized mean difference using the control group standard
deviation as the basis of standardization).

## Details

When `unbiased=TRUE`, the unbiased estimate of the standardized mean
difference (using the control group as the basis of standardization) is
returned (Hedges, 1981). Although the unbiased estimate of the
standardized mean difference is not often reported, at least at the
present time, it is nevertheless made available to those who are
interested in calculating this quantity.

## References

Glass, G. V. (1976). Primary, secondary, and meta-analysis of research.
*Educational Researcher, 5*, 3–8.

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128.

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons and Chapter 3
on one-way ANOVA.)

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Generate sample data.
set.seed(113)
g.T <- rnorm(n = 25, mean = .5, sd = 1)
g.C <- rnorm(n = 25, mean = 0, sd = 1)
smd_c(group_T = g.T, group_C = g.C)
#>  term  value
#>  smd_c 0.37 

M.T <- .66745
M.C <- .24878
sd.c <- 1.1311
n.c <- 25
smd_c(mean_T = M.T, mean_C = M.C, s_C = sd.c)
#>  term  value
#>  smd_c 0.37 
smd_c(mean_T = M.T, mean_C = M.C, s_C = sd.c, n_C = n.c, unbiased = TRUE)
#>  term  value
#>  smd_c 0.358
```
