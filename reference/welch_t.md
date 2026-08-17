# Welch's Separate-Variance *t* Test

Computes Welch's (1947) separate-variance *t* test, the Satterthwaite
(1946) approximation for unequal-variance two-sample inference, and
returns the test statistic, Satterthwaite degrees of freedom, *p*-value,
point estimate of the mean difference, and a confidence interval on the
mean difference, all in a tidy `data.frame`. Unlike Student's
pooled-variance *t* test (`stats::t.test(..., var.equal = TRUE)`),
Welch's test does *not* assume the two populations have equal variances,
and should be the default choice in applied work (Delacre, Lakens, &
Leys, 2017).

## Usage

``` r
welch_t(
  x,
  y,
  mu = 0,
  alternative = c("two_sided", "less", "greater"),
  conf_level = 0.95
)
```

## Arguments

- x, y:

  Numeric vectors of observations from the two groups. The two groups
  are independent and need not be the same length; `NA`s are removed
  from each vector separately.

- mu:

  Null value of the mean difference \\\mu_1 - \mu_2\\. Default `0`.

- alternative:

  One of `"two_sided"` (default; the base-R spelling `"two.sided"` is
  accepted as an alias), `"less"`, or `"greater"`, defining the
  direction of the alternative hypothesis.

- conf_level:

  Confidence level for the CI on the mean difference. Default `0.95`.

## Value

A `data.frame` with rows for the mean difference \\\bar x - \bar y\\,
the Welch *t*-statistic, Satterthwaite degrees of freedom, *p*-value,
the CI lower and upper limits on the mean difference, and the per-group
means, SDs, and *n*.

## Details

**Test statistic.** Welch's *t* is \$\$t \\=\\ \frac{\bar x - \bar y -
\mu_0} {\sqrt{s_1^2 / n_1 + s_2^2 / n_2}},\$\$ which is referred to a
*t* distribution on the Satterthwaite (1946) approximate degrees of
freedom \$\$df \\=\\ \frac{(s_1^2 / n_1 + s_2^2 / n_2)^2} {(s_1^2 /
n_1)^2 / (n_1 - 1) + (s_2^2 / n_2)^2 / (n_2 - 1)}.\$\$

**Why Welch by default.** Student's pooled-variance *t* assumes
\\\sigma_1 = \sigma_2\\; when that assumption fails it has both inflated
and deflated Type I error rates depending on the \\n_1 : n_2\\ ratio
(Ruxton, 2006). Welch's test maintains nominal Type I error across
virtually all combinations of \\\sigma_1 / \sigma_2\\ and \\n_1 / n_2\\,
with no meaningful loss of power when variances are equal. The American
Statistical Association and multiple methodological reviews now
recommend Welch as the default (Delacre et al., 2017; Lakens, 2015).

**Relation to
[`stats::t.test()`](https://rdrr.io/r/stats/t.test.html).** The
numerical results here match `stats::t.test(x, y, var.equal = FALSE)` to
machine precision; this function differs only in returning a tidy
`data.frame` that composes with the rest of DMAR.

## References

Delacre, M., Lakens, D., & Leys, C. (2017). Why psychologists should by
default use Welch's *t*-test instead of Student's *t*-test.
*International Review of Social Psychology, 30*(1), 92–101.
[doi:10.5334/irsp.82](https://doi.org/10.5334/irsp.82)

Lakens, D. (2015, January). Always use Welch's t-test instead of
Student's t-test \[Blog post\]. The 20% Statistician.
<https://daniellakens.blogspot.com/2015/01/always-use-welchs-t-test-instead-of.html>

Ruxton, G. D. (2006). The unequal variance *t*-test is an underused
alternative to Student's *t*-test and the Mann-Whitney *U* test.
*Behavioral Ecology, 17*(4), 688–690.
[doi:10.1093/beheco/ark016](https://doi.org/10.1093/beheco/ark016)

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin, 2*(6), 110–114.

Welch, B. L. (1947). The generalization of "Student's" problem when
several different population variances are involved. *Biometrika,
34*(1/2), 28–35.

## See also

[`t.test`](https://rdrr.io/r/stats/t.test.html),
[`summary_t_test`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md)

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
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Two groups with different variances:
set.seed(113)
x <- rnorm(20, mean = 100, sd = 15)
y <- rnorm(20, mean = 110, sd = 25)
welch_t(x, y)
#>  term            value 
#>  mean_difference -16.9 
#>  t_statistic     -2.41 
#>  df              28.3  
#>  p_value         0.0227
#>  lower_limit     -31.2 
#>  upper_limit     -2.54 
#>  mean_x          100   
#>  mean_y          117   
#>  sd_x            14.3  
#>  sd_y            27.9  
#>  n_x             20    
#>  n_y             20    
#> 
#> Confidence level: 95%

# 2. One-sided test:
welch_t(x, y, alternative = "less")
#>  term            value 
#>  mean_difference -16.9 
#>  t_statistic     -2.41 
#>  df              28.3  
#>  p_value         0.0113
#>  lower_limit     -Inf  
#>  upper_limit     -4.97 
#>  mean_x          100   
#>  mean_y          117   
#>  sd_x            14.3  
#>  sd_y            27.9  
#>  n_x             20    
#>  n_y             20    
#> 
#> Confidence level: 95%

# 3. Side-by-side comparison with base R's stats::t.test().
# The two functions implement the same Welch / Satterthwaite test, so
# the t-statistic, Satterthwaite degrees of freedom, p-value, and the
# CI on the mean difference match exactly. welch_t() differs only in
# what it returns: a data.frame(term, value) rather than a list-
# like htest object. The return composes with dplyr / ggplot2
# pipelines and avoids stringly-typed access like $statistic.
set.seed(113)
a <- rnorm(15, mean = 0,   sd = 1)
b <- rnorm(20, mean = 0.5, sd = 2)

# DMAR (data.frame):
dmar_res <- welch_t(a, b, conf_level = 0.95)
dmar_res
#>  term            value  
#>  mean_difference -1.16  
#>  t_statistic     -2.29  
#>  df              30.1   
#>  p_value         0.0294 
#>  lower_limit     -2.19  
#>  upper_limit     -0.124 
#>  mean_x          -0.0832
#>  mean_y          1.07   
#>  sd_x            1.02   
#>  sd_y            1.93   
#>  n_x             15     
#>  n_y             20     
#> 
#> Confidence level: 95%

# Base R (htest list):
base_res <- stats::t.test(a, b, var.equal = FALSE, conf.level = 0.95)
base_res
#> 
#>  Welch Two Sample t-test
#> 
#> data:  a and b
#> t = -2.2878, df = 30.068, p-value = 0.02935
#> alternative hypothesis: true difference in means is not equal to 0
#> 95 percent confidence interval:
#>  -2.1886355 -0.1241931
#> sample estimates:
#>   mean of x   mean of y 
#> -0.08317142  1.07324287 
#> 

# Verify the four key statistics agree numerically:
pick <- function(term) dmar_res$value[dmar_res$term == term]
stopifnot(
  all.equal(pick("t_statistic"), unname(base_res$statistic)),
  all.equal(pick("df"),          unname(base_res$parameter)),
  all.equal(pick("p_value"),     base_res$p.value),
  all.equal(pick("lower_limit"), base_res$conf.int[1]),
  all.equal(pick("upper_limit"), base_res$conf.int[2])
)
```
