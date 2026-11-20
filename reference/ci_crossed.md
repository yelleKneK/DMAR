# Confidence Interval for the Mean of a Crossed Design of Units and Raters

The confidence interval for the population mean score of a fully crossed
design, in which every rater scores every unit once, when the raters,
like the units, are a sample from a population. The units are the
objects of measurement (persons, essays, recordings, or the topics of an
information retrieval evaluation) and the raters are whoever or whatever
produces the scores (judges, coders, assessors, or automated scoring
models). Because a different panel of raters would have given a
different mean, the standard error carries the variability of the raters
as well as of the units; treating the raters as fixed, as an interval
computed from the unit means alone does, understates it whenever raters
differ. This is the interval
[`ss_aipe_crossed`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md)
plans for, and its variance component rows are the planning values that
function takes.

## Usage

``` r
ci_crossed(ratings, conf_level = 0.95)
```

## Arguments

- ratings:

  A numeric matrix or `data.frame` with one row per unit and one column
  per rater, holding each rater's score of each unit. At least two units
  and two raters are needed, and every cell must hold a finite score:
  the design is complete, with no missing values. A cell may also hold a
  difference of two scores, as when every rater scores every unit under
  two conditions; the interval is then for the mean difference.

- conf_level:

  Confidence level of the interval (default `0.95`).

## Value

A `data.frame` (class `dmar_tbl`) with rows `ci_lower`, `grand_mean`,
and `ci_upper` (the estimate between its limits), `se` (the estimated
standard error of the grand mean), `df` (its Satterthwaite degrees of
freedom), the estimated variance components `var_unit`, `var_rater`, and
`var_residual` (the first two truncated at zero), and the counts
`n_units` and `n_raters`.

## Details

**Model.** The score of unit \\u\\ by rater \\r\\ is \\y\_{ur} = \mu +
a_u + b_r + e\_{ur}\\, with independent unit effects, rater effects, and
residuals of mean zero and variances \\\sigma^2_u\\, \\\sigma^2_r\\, and
\\\sigma^2_e\\; the residual holds the unit by rater interaction, which
one score per cell cannot separate from error. This is the design Shrout
and Fleiss (1979) call Case 2, a random sample of judges each of whom
rates every target, and it is the crossed design with one random facet
of generalizability theory (Brennan, 2001). The grand mean \\\bar y\\ of
the \\n_u \times n_r\\ array estimates \\\mu\\ with variance \$\$V =
\frac{\sigma^2_u}{n_u} + \frac{\sigma^2_r}{n_r} + \frac{\sigma^2_e}{n_u
n_r}.\$\$

**Estimate and interval.** With the mean squares \\MS_u\\, \\MS_r\\, and
\\MS_e\\ of the two-way analysis of variance for units, raters, and
residual, the ANOVA estimates of the components are \\(MS_u - MS_e) /
n_r\\, \\(MS_r - MS_e) / n_u\\, and \\MS_e\\, each of the first two set
to zero when it is negative. Plugged into \\V\\ they give \$\$\hat V =
\frac{k_u MS_u + k_r MS_r + (1 - k_u - k_r) MS_e}{n_u n_r},\$\$ with
\\k_u\\ equal to 1 when \\MS_u \> MS_e\\ and 0 otherwise, and \\k_r\\
likewise for raters. The interval is \\\bar y \pm t\_{\nu} \sqrt{\hat
V}\\, where the degrees of freedom \\\nu\\ are those Satterthwaite
(1946) gives for a linear combination of independent mean squares,
\$\$\nu = \frac{\left(\sum_k a_k MS_k\right)^2}{\sum_k (a_k MS_k)^2 /
\nu_k},\$\$ with the coefficients \\a_k\\ of the numerator above and
\\\nu_k\\ the degrees of freedom of each mean square: \\n_u - 1\\,
\\n_r - 1\\, and \\(n_u - 1)(n_r - 1)\\. The truncation keeps \\\hat V\\
at or above \\MS_e / (n_u n_r)\\, so it is never negative. When few
raters carry a large rater component, \\\nu\\ is near \\n_r - 1\\ and
the interval is wide: the data then hold little information about how
far the mean of another panel of raters would lie from this one's.

**An information retrieval evaluation.** To compare two retrieval
systems over topics judged by several relevance assessors (people, or
automated judges such as language models), let each row be a topic, each
column an assessor, and each cell the difference between the two
systems' effectiveness scores on that topic computed from that
assessor's judgments. The interval is then for the mean difference
between the systems over the populations of topics and assessors.

**Designs chosen from the data.** The interval assumes the numbers of
units and raters were fixed in advance, as
[`ss_aipe_crossed`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md)
fixes them, and not chosen from the data, as a two-stage sequential
procedure chooses them from a pilot array.

## References

Brennan, R. L. (2001). *Generalizability theory*. Springer.
[doi:10.1007/978-1-4757-3456-0](https://doi.org/10.1007/978-1-4757-3456-0)

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin, 2*(6), 110–114.
[doi:10.2307/3002019](https://doi.org/10.2307/3002019)

Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
assessing rater reliability. *Psychological Bulletin, 86*(2), 420–428.
[doi:10.1037/0033-2909.86.2.420](https://doi.org/10.1037/0033-2909.86.2.420)

## See also

[`ss_aipe_crossed`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md)
to plan the numbers of units and raters,
[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md) for the
intraclass correlations of the same array, and
[`variance_components_mls`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md)
for intervals on variance components.

Other confidence intervals:
[`ci_proportion()`](https://yelleknek.github.io/DMAR/reference/ci_proportion.md),
[`ci_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/ci_proportion_paired.md),
[`ci_two_phase()`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Shrout and Fleiss (1979, Table 2): four judges each rated the same
# six targets.
ratings <- matrix(c( 9, 2, 5, 8,
                     6, 1, 3, 2,
                     8, 4, 6, 8,
                     7, 1, 2, 6,
                    10, 5, 6, 9,
                     6, 2, 4, 7),
                  nrow = 6, byrow = TRUE)
fit <- ci_crossed(ratings)
fit
#>  term         value
#>  ci_lower     1.827
#>  grand_mean   5.292
#>  ci_upper     8.756
#>  se           1.334
#>  df           4.837
#>  var_unit     2.556
#>  var_rater    5.244
#>  var_residual 1.019
#>  n_units      6    
#>  n_raters     4    
#> 
#> Confidence level: 95%

# The judges differ a great deal in leniency, and with only four of
# them the interval for the mean over the population of judges is wide.
# Treating these four judges as the only ones of interest ignores that
# variability; the interval from the six target means is much narrower.
t.test(rowMeans(ratings))$conf.int
#> [1] 3.532362 7.050972
#> attr(,"conf.level")
#> [1] 0.95

# The variance components are planning values for a larger study, here
# one whose 95% interval for the mean is no wider than 2 points.
v <- setNames(fit$value, fit$term)
ss_aipe_crossed(true_var_unit = v[["var_unit"]], true_var_rater = v[["var_rater"]],
                true_var_residual = v[["var_residual"]], width = 2)
#>  term               value 
#>  necessary_n_units  22    
#>  necessary_n_raters 40    
#>  total_judgments    880   
#>  cost               880   
#>  expected_width     1.997 
#>  df                 55.903
#>  rater_floor        20.146
#>  true_var_unit      2.556 
#>  true_var_rater     5.244 
#>  true_var_residual  1.019 
#>  width              2     
#>  conf_level         0.950 
#> 
#> Confidence level: 95%
```
