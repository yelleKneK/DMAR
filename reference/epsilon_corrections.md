# Greenhouse-Geisser, Huynh-Feldt, and Lower-Bound Epsilon Corrections

Computes the three standard sphericity-correction factors for the
univariate within-subjects *F* test. When sphericity holds, all three
equal 1; departures from sphericity reduce them, deflating the effective
degrees of freedom and thus tempering the inflated Type I error rate of
the unadjusted univariate test.

## Usage

``` r
epsilon_corrections(x, id = NULL, time = NULL, outcome = NULL)
```

## Arguments

- x:

  Either an \\n \times k\\ numeric matrix or `data.frame` (rows =
  subjects, columns = repeated measurements); *or* a long-format
  `data.frame` together with `id`, `time`, and `outcome` column names.

- id:

  Column name in `x` identifying the subject when `x` is in long format
  (`NULL` otherwise).

- time:

  Column name in `x` identifying the within-subjects factor level when
  `x` is in long format (`NULL` otherwise).

- outcome:

  Column name in `x` identifying the dependent variable when `x` is in
  long format (`NULL` otherwise).

## Value

A `data.frame` with columns `epsilon_method` (`"Greenhouse-Geisser"`,
`"Huynh-Feldt"`, `"lower_bound"`) and `epsilon` (the correction factor
in \\\[1/(k-1), 1\]\\).

## Details

For an \\(k - 1) \times (k - 1)\\ covariance matrix \\\hat\Sigma_C\\ of
orthonormal contrasts among the \\k\\ repeated measurements (with
eigenvalues \\\lambda_1, \ldots, \lambda\_{k-1}\\):
\$\$\hat\varepsilon\_{\mathrm{GG}} = \frac{(\sum \lambda_i)^2}{(k -
1)\\\sum \lambda_i^2},\$\$ \$\$\hat\varepsilon\_{\mathrm{HF}} =
\min\\\Bigl(1,\\ \frac{n(k - 1)\hat\varepsilon\_{\mathrm{GG}} - 2}{(k -
1)\bigl(n - 1 - (k - 1)\hat\varepsilon\_{\mathrm{GG}}\bigr)}\Bigr),\$\$
\$\$\hat\varepsilon\_{\mathrm{LB}} = \frac{1}{k - 1}.\$\$ The
Greenhouse-Geisser \\\hat\varepsilon\\ tends to be conservative; the
Huynh-Feldt correction adjusts upward to be (approximately) unbiased;
the lower bound is the worst-case adjustment.

## References

Greenhouse, S. W., & Geisser, S. (1959). On methods in the analysis of
profile data. *Psychometrika, 24*(2), 95–112.

Huynh, H., & Feldt, L. S. (1976). Estimation of the Box correction for
degrees of freedom from sample data in randomized block and split-plot
designs. *Journal of Educational Statistics, 1*(1), 69–82.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`mauchly_test`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`anova_within`](https://yelleknek.github.io/DMAR/reference/anova_within.md)

Other within-subjects analysis:
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`anova_within_two_way()`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
Y <- matrix(rnorm(20 * 4), nrow = 20)
epsilon_corrections(Y)
#>  epsilon_method     epsilon
#>  Greenhouse-Geisser 0.787  
#>  Huynh-Feldt        0.906  
#>  lower_bound        0.333  
```
