# Contrast Among Study Effect Sizes (Rosenthal-Rubin)

Tests a focused hypothesis about *differences* among independent study
effect sizes by the method of Rosenthal and Rubin (1982): given effects
\\y_i\\ with sampling variances \\v_i\\ and contrast weights
\\\lambda_i\\ summing to zero, \$\$z \\=\\ \frac{\sum \lambda_i
y_i}{\sqrt{\sum \lambda_i^2 v_i}}\$\$ is referred to the standard
normal. This is how a meta-analyst asks a pointed moderator question
(“do the effects decline with weeks of prior teacher-student contact?”)
rather than the diffuse heterogeneity question (“do the effects differ
at all?”). Raudenbush (1984) used exactly this test for the teacher
expectancy literature, with weights inversely proportional to weeks of
prior contact.

## Usage

``` r
meta_contrast(yi, vi, weights, center = TRUE)
```

## Arguments

- yi:

  Numeric vector of study effect sizes (any metric whose sampling
  distribution is approximately normal; standardized mean differences
  and Fisher-z correlations qualify).

- vi:

  Sampling variances of `yi`, one per study.

- weights:

  Contrast weights, one per study. If they do not already sum to zero
  they are mean-centered (with a message) when `center = TRUE`, the
  convenient route for weights built from a moderator such as
  `1 / (weeks + 2)`.

- center:

  Logical: mean-center `weights` that do not sum to zero? Default
  `TRUE`.

## Value

A `data.frame` (class `dmar_tbl`) with the contrast `estimate` (\\\sum
\lambda_i y_i\\), its `se`, the `z` statistic, the two-sided `p_value`,
and `k`.

## Details

The two-sided *p*-value is reported; halve it for a directional
hypothesis stated in advance (Raudenbush's \\z = 2.75\\ carried the
one-tailed \\p = .003\\). Dividing the squared contrast \\z^2\\ by the
total heterogeneity statistic \\Q\\ from
[`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md) gives
the proportion of between-study heterogeneity the contrast accounts for,
the meta-analytic analog of a contrast's share of the between-group sum
of squares.

## References

Raudenbush, S. W. (1984). Magnitude of teacher expectancy effects on
pupil IQ as a function of the credibility of expectancy induction: A
synthesis of findings from 18 experiments. *Journal of Educational
Psychology, 76*(1), 85–97.

Rosenthal, R., & Rubin, D. B. (1982). Comparing effect sizes of
independent studies. *Psychological Bulletin, 92*(2), 500–504.

## See also

[`meta_es`](https://yelleknek.github.io/DMAR/reference/meta_es.md) for
the pooled effect and the total heterogeneity the contrast partitions;
[`combine_p`](https://yelleknek.github.io/DMAR/reference/combine_p.md)
for combined significance tests;
[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)
for the single-study ANOVA analog.

Other meta-analysis:
[`combine_p()`](https://yelleknek.github.io/DMAR/reference/combine_p.md),
[`meta_es()`](https://yelleknek.github.io/DMAR/reference/meta_es.md),
[`meta_r()`](https://yelleknek.github.io/DMAR/reference/meta_r.md),
[`meta_smd()`](https://yelleknek.github.io/DMAR/reference/meta_smd.md),
[`plot_forest()`](https://yelleknek.github.io/DMAR/reference/plot_forest.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Raudenbush (1984): do expectancy effects decline with weeks of prior
# teacher-student contact? Weights inversely proportional to weeks + 2,
# study-level data (Pellegrini & Hicks merged), d variances from the
# standard large-sample formula.
data(teacher_expectancy)
study <- teacher_expectancy[-c(4, 5), ]
d  <- append(study$d, 0.52, after = 3)
wk <- append(study$weeks, 0, after = 3)
ne <- append(study$n_experimental, 22, after = 3)
nc <- append(study$n_control, 22, after = 3)
v  <- (ne + nc) / (ne * nc) + d^2 / (2 * (ne + nc))
meta_contrast(d, v, weights = 1 / (wk + 2))
#> Contrast weights mean-centered to sum to zero.
#>  term     value 
#>  estimate 0.432 
#>  se       0.156 
#>  z        2.76  
#>  p_value  0.0057
#>  k        18    
# z near 2.75: the better teachers knew their pupils, the smaller the
# expectancy effect (one-tailed p = .003 in the paper).
```
