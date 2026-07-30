# Mixed-design multivariate ANOVA with all four test statistics.

Computes the multivariate analysis of variance for a mixed-design (one
between-subjects factor and one within-subjects factor), returning
Wilks's \\\Lambda\\, Pillai's trace, Hotelling-Lawley trace, and Roy's
largest root, each with the associated *F*-approximation, degrees of
freedom, and *p*-value. The three effects, between-subjects (*A*),
within-subjects (*B*), and the interaction (\\A \times B\\), are tested
separately. Wraps [`Anova`](https://rdrr.io/pkg/car/man/Anova.html) and
returns the result in the tidy DMAR style.

## Usage

``` r
manova_mixed(data, within, between)
```

## Arguments

- data:

  A `data.frame` in wide format, one row per subject, with the
  within-subjects measurements in separate columns and the
  between-subjects factor as one additional column.

- within:

  Character vector of column names holding the repeated measures values
  (one column per level of the within- subjects factor). Must be in the
  canonical level order.

- between:

  Character name of the between-subjects factor column in `data`.

## Value

A tidy `data.frame` with rows for each of the three effects crossed with
each of the four multivariate statistics. Columns: `effect`,
`statistic_name`, `statistic_value`, `F_approx`, `df_1`, `df_2`,
`p_value`.

## Details

**The four statistics.** For an effect with \\H\\ and \\E\\ hypothesis-
and error-cross-products matrices:

- Wilks's \\\Lambda = \det(E) / \det(E + H)\\

- Pillai's trace \\V = \mathrm{tr}(H (E + H)^{-1})\\

- Hotelling-Lawley trace \\T_0^2 = \mathrm{tr}(H E^{-1})\\

- Roy's largest root \\\theta = \lambda_1(H E^{-1})\\

**When to use which.** Pillai's trace is the most robust to departures
from the multivariate normal / homogeneous-covariance assumptions.
Wilks's \\\Lambda\\ is the most widely reported. Roy's largest root is
the most powerful when the alternative concentrates on a single
dimension. The four statistics agree exactly when the effect has 1
numerator degree of freedom.

**Dependency.** Requires the car package on CRAN.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 14.)

Rencher, A. C., & Christensen, W. F. (2012). *Methods of multivariate
analysis* (4th ed.). Wiley.

## See also

[`Anova`](https://rdrr.io/pkg/car/man/Anova.html),
[`manova`](https://rdrr.io/r/stats/manova.html),
[`anova_within_two_way`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md)

Other hypothesis tests:
[`ancova`](https://yelleknek.github.io/DMAR/reference/ancova.md)`()`,
[`anova_within`](https://yelleknek.github.io/DMAR/reference/anova_within.md)`()`,
[`compare_cov_structures`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md)`()`,
[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)`()`,
[`correlations_test`](https://yelleknek.github.io/DMAR/reference/correlations_test.md)`()`,
[`dunnett_ci`](https://yelleknek.github.io/DMAR/reference/dunnett_ci.md)`()`,
[`mauchly_test`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md)`()`,
[`mixed_anova`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md)`()`,
[`obrien_test`](https://yelleknek.github.io/DMAR/reference/obrien_test.md)`()`,
[`pairwise_within`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md)`()`,
[`randomization_test_paired`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md)`()`,
[`scheffe_ci`](https://yelleknek.github.io/DMAR/reference/scheffe_ci.md)`()`,
[`simple_effects_AB`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md)`()`,
[`summary_t_test`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md)`()`,
[`tost_r`](https://yelleknek.github.io/DMAR/reference/tost_r.md)`()`,
[`tost_smd`](https://yelleknek.github.io/DMAR/reference/tost_smd.md)`()`,
[`tukey_kramer_ci`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md)`()`,
[`welch_t`](https://yelleknek.github.io/DMAR/reference/welch_t.md)`()`

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
set.seed(113)
n_per <- 10
d <- data.frame(
  subject = factor(1:(2 * n_per)),
  group   = factor(rep(c("A", "B"), each = n_per)),
  t1      = c(rnorm(n_per, 0,   1), rnorm(n_per, 0,   1)),
  t2      = c(rnorm(n_per, 0.4, 1), rnorm(n_per, 0.8, 1)),
  t3      = c(rnorm(n_per, 0.8, 1), rnorm(n_per, 1.6, 1))
)
manova_mixed(d, within = c("t1", "t2", "t3"), between = "group")
#>         effect   statistic_name statistic_value  F_approx df_1 df_2    p_value
#> 1            A           Pillai      0.08540425 1.6808262    1   18 0.21119035
#> 2            A            Wilks      0.91459575 1.6808262    1   18 0.21119035
#> 3            A Hotelling-Lawley      0.09337924 1.6808262    1   18 0.21119035
#> 4            A              Roy      0.09337924 1.6808262    1   18 0.21119035
#> 5            B           Pillai      0.37096011 5.0126566    2   17 0.01944306
#> 6            B            Wilks      0.62903989 5.0126566    2   17 0.01944306
#> 7            B Hotelling-Lawley      0.58972431 5.0126566    2   17 0.01944306
#> 8            B              Roy      0.58972431 5.0126566    2   17 0.01944306
#> 9  interaction           Pillai      0.08217611 0.7610359    2   17 0.48245247
#> 10 interaction            Wilks      0.91782389 0.7610359    2   17 0.48245247
#> 11 interaction Hotelling-Lawley      0.08953364 0.7610359    2   17 0.48245247
#> 12 interaction              Roy      0.08953364 0.7610359    2   17 0.48245247
# }
```
